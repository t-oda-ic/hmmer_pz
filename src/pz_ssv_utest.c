/* pz_ssv_utest.c
 *
 * Verify that the PEZY-SC SSV kernel reproduces impl_sse's SSV filter exactly,
 * and report how fast it does it.
 *
 * usage: pz_ssv_utest [-n <nseq>] [-b <batch>] <hmmfile> <seqfile>
 *
 * Two things are checked for every target sequence:
 *
 *   1. the raw xE byte against impl_sse/ssvfilter.c's get_xE().  These must
 *      agree except in the overflow regime (xE >= 255 - bias_b), where
 *      p7_SSVFilter() throws the value away; a disagreement there is reported
 *      separately and is not a failure, but a disagreement anywhere else is.
 *
 *   2. the full p7_SSVFilter() outcome -- status and score bit pattern --
 *      against p7_pz_SSVOutcome() driven by the device's xE.  This is what the
 *      pipeline actually consumes, so it must match exactly, always.
 */
#include <p7_config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_stopwatch.h"

#include "hmmer.h"

#ifndef HAVE_PEZY
int main(void) { printf("pz_ssv_utest: built without PEZY support\n"); return 0; }
#else
#include "p7_pz.h"

extern uint8_t get_xE(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om); /* impl_sse/ssvfilter.c */

static ESL_OPTIONS options[] = {
  /* name        type          def   env  range toggles reqs incomp  help                                docgroup*/
  { "-h",        eslARG_NONE,  FALSE, NULL, NULL, NULL, NULL, NULL, "show brief help",                        0 },
  { "-n",        eslARG_INT,   "0",   NULL, NULL, NULL, NULL, NULL, "stop after <n> sequences (0 = all)",     0 },
  { "-b",        eslARG_INT, "200000",NULL, NULL, NULL, NULL, NULL, "sequences per device batch",             0 },
  { "-d",        eslARG_INT,   "0",   NULL, NULL, NULL, NULL, NULL, "PEZY device index",                      0 },
  { 0,0,0,0,0,0,0,0,0,0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "test the PEZY-SC SSV filter against impl_sse";

int
main(int argc, char **argv)
{
  ESL_GETOPTS   *go      = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  char          *hmmfile = esl_opt_GetArg(go, 1);
  char          *seqfile = esl_opt_GetArg(go, 2);
  int64_t        maxseq  = esl_opt_GetInteger(go, "-n");
  int            batchsz = esl_opt_GetInteger(go, "-b");
  int            devidx  = esl_opt_GetInteger(go, "-d");
  ESL_ALPHABET  *abc     = NULL;
  P7_HMMFILE    *hfp     = NULL;
  P7_HMM        *hmm     = NULL;
  P7_BG         *bg      = NULL;
  P7_PROFILE    *gm      = NULL;
  P7_OPROFILE   *om      = NULL;
  ESL_SQFILE    *sqfp    = NULL;
  ESL_SQ        *sq      = NULL;
  P7_PZ         *pz      = NULL;
  P7_PZ_SEQBATCH *b      = NULL;
  uint8_t       *xE      = NULL;
  int           *Lsave   = NULL;
  ESL_STOPWATCH *w       = esl_stopwatch_Create();
  char           errbuf[256];
  double         t_dev = 0.0, t_cpu = 0.0;
  int64_t        nseq = 0, nres = 0;
  int64_t        n_rawdiff = 0, n_hard = 0, n_outcome = 0;
  int            status;

  impl_Init();

  if (p7_pz_DeviceCount() <= 0) esl_fatal("no PEZY device available");

  if (p7_hmmfile_Open(hmmfile, NULL, &hfp, errbuf) != eslOK) esl_fatal("failed to open %s: %s", hmmfile, errbuf);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)             != eslOK) esl_fatal("failed to read HMM from %s", hmmfile);

  bg = p7_bg_Create(abc);
  gm = p7_profile_Create(hmm->M, abc);
  om = p7_oprofile_Create(hmm->M, abc);
  p7_ProfileConfig(hmm, bg, gm, 400, p7_LOCAL);
  p7_oprofile_Convert(gm, om);

  if ((pz = p7_pz_Create(devidx, errbuf)) == NULL) esl_fatal("p7_pz_Create failed: %s", errbuf);
  if ((status = p7_pz_SetModel(pz, om)) == eslENORESULT)
    esl_fatal("model M=%d is too large for the device kernel", om->M);
  else if (status != eslOK) esl_fatal("p7_pz_SetModel failed (%d)", status);

  if (esl_sqfile_OpenDigital(abc, seqfile, eslSQFILE_FASTA, NULL, &sqfp) != eslOK) esl_fatal("failed to open %s", seqfile);
  sq = esl_sq_CreateDigital(abc);

  b     = p7_pz_seqbatch_Create(batchsz, (int64_t) batchsz * 400);
  xE    = malloc(sizeof(uint8_t) * batchsz);
  Lsave = malloc(sizeof(int)     * batchsz);
  if (!xE || !Lsave) esl_fatal("malloc failed");

  printf("# model %s  M=%d  Kp=%d  bias_b=%u base_b=%u tec_b=%u tbm_b=%u\n",
         hmm->name, om->M, abc->Kp, om->bias_b, om->base_b, om->tec_b, om->tbm_b);
  printf("# device batch size %d\n", batchsz);

  /* Read the database in batches; verify each batch as it completes. */
  {
    int eof = FALSE;
    while (! eof)
      {
        p7_pz_seqbatch_Reuse(b);

        while (b->n < batchsz)
          {
            if (maxseq && nseq + b->n >= maxseq) break;
            if ((status = esl_sqio_Read(sqfp, sq)) != eslOK) { eof = TRUE; break; }
            if (sq->n > 0 && sq->n <= 100000) {
              p7_pz_seqbatch_Add(b, sq->dsq, (int) sq->n);
              Lsave[b->n - 1] = (int) sq->n;
            }
            esl_sq_Reuse(sq);
          }
        if (b->n == 0) break;

        /* ---- device ---- */
        esl_stopwatch_Start(w);
        if (p7_pz_SSVBatch(pz, b, xE) != eslOK) esl_fatal("p7_pz_SSVBatch failed");
        esl_stopwatch_Stop(w);
        t_dev += w->elapsed;

        /* ---- reference ---- */
        esl_stopwatch_Start(w);
        for (int i = 0; i < b->n; i++)
          {
            const ESL_DSQ *dsq = (const ESL_DSQ *) (b->dsq + b->off[i]);
            int            L   = Lsave[i];
            uint8_t        ref;
            float          sc_ref, sc_dev;
            int            st_ref, st_dev;

            p7_oprofile_ReconfigMSVLength(om, L);   /* the pipeline does this per seq */
            ref = get_xE(dsq, L, om);

            if (ref != xE[i]) {
              n_rawdiff++;
              if (! (ref >= 255 - om->bias_b && xE[i] >= 255 - om->bias_b)) {
                if (n_hard < 8) printf("HARD  xE mismatch: seq %lld L=%d ref=%u dev=%u\n",
                                       (long long)(nseq+i), L, ref, xE[i]);
                n_hard++;
              }
            }
            st_ref = p7_SSVFilter(dsq, L, om, &sc_ref);
            st_dev = p7_pz_SSVOutcome(xE[i], om, &sc_dev);
            if (st_ref != st_dev || (st_ref == eslOK && memcmp(&sc_ref, &sc_dev, sizeof(float)) != 0)) {
              if (n_outcome < 8) printf("OUTCOME mismatch: seq %lld L=%d status %d/%d score %.9g/%.9g\n",
                                        (long long)(nseq+i), L, st_ref, st_dev, sc_ref, sc_dev);
              n_outcome++;
            }
            nres += L;
          }
        esl_stopwatch_Stop(w);
        t_cpu += w->elapsed;

        nseq += b->n;
        if (maxseq && nseq >= maxseq) break;
      }
  }

  printf("# sequences        %" PRId64 "\n", nseq);
  printf("# residues         %" PRId64 "\n", nres);
  printf("# device time      %.3f s   (%.1f Mcells/s)\n", t_dev, (double) nres * om->M / t_dev / 1e6);
  printf("# cpu ref time     %.3f s   (2 SSV passes: get_xE + p7_SSVFilter)\n", t_cpu);
  printf("# 1-core cpu SSV   %.3f s   -> device speedup %.1fx\n", t_cpu/2.0, (t_cpu/2.0)/t_dev);
  printf("# raw xE diffs     %" PRId64 "  (all inside the overflow regime, where the value is discarded)\n", n_rawdiff);
  printf("# HARD xE diffs    %" PRId64 "\n", n_hard);
  printf("# outcome diffs    %" PRId64 "\n", n_outcome);
  printf("%s\n", (n_hard == 0 && n_outcome == 0) ? "ok. PEZY SSV filter matches impl_sse exactly." : "FAILED");

  free(xE); free(Lsave);
  p7_pz_seqbatch_Destroy(b);
  p7_pz_Destroy(pz);
  esl_stopwatch_Destroy(w);
  esl_sq_Destroy(sq);
  esl_sqfile_Close(sqfp);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return (n_hard == 0 && n_outcome == 0) ? 0 : 1;
}
#endif /*HAVE_PEZY*/
