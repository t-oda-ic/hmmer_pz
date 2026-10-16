/* pz_ssvmulti_utest.c
 *
 * Verify the many-models x many-queries SSV kernel (P7_PZMM) against
 * impl_sse, over a whole pressed profile database.
 *
 * usage: pz_ssvmulti_utest [-n <nmodel>] [-t <threads>] <hmmfile> <seqfile>
 */
#include <p7_config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_stopwatch.h"

#include "hmmer.h"

#ifndef HAVE_PEZY
int main(void) { printf("pz_ssvmulti_utest: built without PEZY support\n"); return 0; }
#else
#include "p7_pz.h"
#include "p7_hmmcache.h"

extern uint8_t get_xE(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om); /* impl_sse/ssvfilter.c */

static ESL_OPTIONS options[] = {
  { "-h", eslARG_NONE,  FALSE, NULL, NULL, NULL, NULL, NULL, "show brief help",                    0 },
  { "-n", eslARG_INT,     "0", NULL, NULL, NULL, NULL, NULL, "check only the first <n> models",    0 },
  { "-t", eslARG_INT,     "4", NULL, NULL, NULL, NULL, NULL, "threads for the table build",        0 },
  { 0,0,0,0,0,0,0,0,0,0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "test the many-model PEZY SSV kernel against impl_sse";

int
main(int argc, char **argv)
{
  ESL_GETOPTS  *go      = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  char         *hmmfile = esl_opt_GetArg(go, 1);
  char         *seqfile = esl_opt_GetArg(go, 2);
  int           ncheck  = esl_opt_GetInteger(go, "-n");
  int           nthread = esl_opt_GetInteger(go, "-t");
  ESL_ALPHABET *abc     = NULL;
  P7_HMMCACHE  *cache   = NULL;
  ESL_SQFILE   *sqfp    = NULL;
  ESL_SQ       *sq      = NULL;
  ESL_SQ      **qs      = NULL;
  P7_PZ        *pz      = NULL;
  P7_PZMM      *mm      = NULL;
  uint8_t      *xE      = NULL;
  ESL_DSQ     **dsq     = NULL;
  int          *L       = NULL;
  ESL_STOPWATCH *w = esl_stopwatch_Create();
  char          errbuf[eslERRBUFSIZE];
  int           nq = 0, i, status;
  long          nraw = 0, nhard = 0, nout = 0, nhost = 0, npair = 0;
  int64_t       cells = 0;
  double        t_dev, t_cpu, t_build;

  impl_Init();
  if (p7_pz_DeviceCount() <= 0) esl_fatal("no PEZY device available");

  abc = esl_alphabet_Create(eslAMINO);
  if (esl_sqfile_OpenDigital(abc, seqfile, eslSQFILE_FASTA, NULL, &sqfp) != eslOK) esl_fatal("cannot open %s", seqfile);
  sq = esl_sq_CreateDigital(abc);
  qs = malloc(sizeof(ESL_SQ *) * 100000);
  while (esl_sqio_Read(sqfp, sq) == eslOK) { qs[nq++] = sq; sq = esl_sq_CreateDigital(abc); }
  esl_sq_Destroy(sq);
  printf("# %d query sequences\n", nq);

  if (p7_hmmcache_Open(hmmfile, &cache, errbuf) != eslOK) esl_fatal("cache: %s", errbuf);
  printf("# %u models\n", cache->n);

  if ((pz = p7_pz_Create(0, errbuf)) == NULL) esl_fatal("p7_pz_Create: %s", errbuf);

  esl_stopwatch_Start(w);
  if ((status = p7_pzmm_Create(pz, cache->list, (int) cache->n, nthread, &mm, errbuf)) != eslOK)
    esl_fatal("p7_pzmm_Create failed (%d): %s", status, errbuf);
  esl_stopwatch_Stop(w); t_build = w->elapsed;
  printf("# tables built+uploaded in %.3f s (%d of %u models on the device, %d host-only)\n",
         t_build, p7_pzmm_NDevice(mm), cache->n, (int) cache->n - p7_pzmm_NDevice(mm));

  dsq = malloc(sizeof(ESL_DSQ *) * nq);
  L   = malloc(sizeof(int) * nq);
  for (i = 0; i < nq; i++) { dsq[i] = qs[i]->dsq; L[i] = (int) qs[i]->n; }
  xE = malloc((size_t) cache->n * nq);

  esl_stopwatch_Start(w);
  if ((status = p7_pzmm_SSV(pz, mm, (const ESL_DSQ **) dsq, L, nq, xE)) != eslOK) esl_fatal("p7_pzmm_SSV failed (%d)", status);
  esl_stopwatch_Stop(w); t_dev = w->elapsed;

  /* ---- reference ---- */
  {
    uint32_t m, mlim = (ncheck > 0 && (uint32_t) ncheck < cache->n) ? (uint32_t) ncheck : cache->n;
    esl_stopwatch_Start(w);
    for (m = 0; m < mlim; m++)
      {
        P7_OPROFILE *om = cache->list[m];
        if (! p7_pzmm_ModelOnDevice(mm, (int) m)) { nhost++; continue; }
        for (i = 0; i < nq; i++)
          {
            uint8_t ref, dev = xE[(size_t) m * nq + i];
            float   sc_ref, sc_dev;
            int     st_ref, st_dev;

            p7_oprofile_ReconfigLength(om, L[i]);
            ref = get_xE(dsq[i], L[i], om);
            npair++; cells += (int64_t) om->M * L[i];

            if (ref != dev) {
              nraw++;
              if (! (ref >= 255 - om->bias_b && dev >= 255 - om->bias_b)) {
                if (nhard < 8) printf("HARD  mismatch: model %u (%s, M=%d) query %d: ref=%u dev=%u\n",
                                      m, om->name, om->M, i, ref, dev);
                nhard++;
              }
            }
            st_ref = p7_SSVFilter(dsq[i], L[i], om, &sc_ref);
            st_dev = p7_pz_SSVOutcome(dev, om, &sc_dev);
            if (st_ref != st_dev || (st_ref == eslOK && memcmp(&sc_ref, &sc_dev, sizeof(float)) != 0)) {
              if (nout < 8) printf("OUTCOME mismatch: model %u (%s) query %d: status %d/%d score %.9g/%.9g\n",
                                   m, om->name, i, st_ref, st_dev, sc_ref, sc_dev);
              nout++;
            }
          }
      }
    esl_stopwatch_Stop(w); t_cpu = w->elapsed;
  }

  printf("# (model,query) pairs   %ld\n", npair);
  printf("# host-only models      %ld\n", nhost);
  printf("# SSV cells             %lld\n", (long long) cells);
  printf("# device SSV            %.3f s   (%.1f Gcells/s)\n", t_dev, cells/t_dev/1e9);
  printf("# cpu reference         %.3f s   (2 passes: get_xE + p7_SSVFilter)\n", t_cpu);
  printf("# raw xE diffs          %ld  (overflow regime only)\n", nraw);
  printf("# HARD xE diffs         %ld\n", nhard);
  printf("# outcome diffs         %ld\n", nout);
  printf("%s\n", (nhard == 0 && nout == 0) ? "ok. many-model PEZY SSV matches impl_sse exactly." : "FAILED");

  free(xE); free(dsq); free(L);
  for (i = 0; i < nq; i++) esl_sq_Destroy(qs[i]);
  free(qs);
  p7_pzmm_Destroy(mm);
  p7_pz_Destroy(pz);
  p7_hmmcache_Close(cache);
  esl_sqfile_Close(sqfp);
  esl_stopwatch_Destroy(w);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return (nhard == 0 && nout == 0) ? 0 : 1;
}
#endif /*HAVE_PEZY*/
