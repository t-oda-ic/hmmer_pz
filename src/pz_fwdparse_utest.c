/* pz_fwdparse_utest.c
 *
 * Verify that the PEZY-SC parser-mode Forward kernel reproduces
 * impl_sse/fwdback.c's p7_ForwardParser() bit for bit, and report how fast it
 * does it.
 *
 * usage: pz_fwdparse_utest [-n <nseq>] [-b <batch>] [-d <dev>] <hmmfile> <seqfile>
 *
 * Two things are compared for every target sequence, both bit-exactly:
 *
 *   1. the whole xmx array -- (L+1) rows of {E,N,J,B,C,SCALE} -- because that
 *      is the entire result of a parser-mode pass and everything downstream
 *      (p7_BackwardParser, p7_DomainDecoding) reads it and nothing else.
 *      Compared as raw 32-bit patterns: "close" is not good enough, since the
 *      whole point of the port is that HMMER's output does not move.
 *
 *   2. the Forward score, from p7_pz_FwdParserScore() driven by the device's
 *      xmx, against p7_ForwardParser()'s own opt_sc.  This is what the F3
 *      filter consumes.  It exercises the one operation that stays on the
 *      host: the log() of the scale factors.
 *
 * A mismatch in either is a failure.  Unlike the SSV test there is no benign
 * regime -- there is no overflow value that the caller discards.
 */
#include <p7_config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_stopwatch.h"

#include "hmmer.h"

#ifndef HAVE_PEZY
int main(void) { printf("pz_fwdparse_utest: built without PEZY support\n"); return 0; }
#else
#include "p7_pz.h"
#include "p7_pz_fbparse.h"

static ESL_OPTIONS options[] = {
  /* name      type          def    env  range toggles reqs incomp  help                                  docgroup*/
  { "-h",      eslARG_NONE,  FALSE, NULL, NULL, NULL, NULL, NULL, "show brief help",                          0 },
  { "-n",      eslARG_INT,   "0",   NULL, NULL, NULL, NULL, NULL, "stop after <n> sequences (0 = all)",       0 },
  { "-b",      eslARG_INT,   "8192",NULL, NULL, NULL, NULL, NULL, "sequences per device batch",               0 },
  { "-d",      eslARG_INT,   "0",   NULL, NULL, NULL, NULL, NULL, "PEZY device index",                        0 },
  { "-L",      eslARG_INT,   "400", NULL, NULL, NULL, NULL, NULL, "length the profile is configured for",     0 },
  { "-x",      eslARG_INT,   "0",   NULL, NULL, NULL, NULL, NULL, "skip targets longer than <x> (0 = no cap)", 0 },
  { 0,0,0,0,0,0,0,0,0,0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "test the PEZY-SC parser-mode Forward and Backward against impl_sse";

int
main(int argc, char **argv)
{
  ESL_GETOPTS   *go      = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  char          *hmmfile = esl_opt_GetArg(go, 1);
  char          *seqfile = esl_opt_GetArg(go, 2);
  int64_t        maxseq  = esl_opt_GetInteger(go, "-n");
  int            batchsz = esl_opt_GetInteger(go, "-b");
  int            devidx  = esl_opt_GetInteger(go, "-d");
  int            cfgL    = esl_opt_GetInteger(go, "-L");
  int            capL    = esl_opt_GetInteger(go, "-x");
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
  P7_OMX        *ox      = NULL;
  P7_OMX        *oxb     = NULL;
  float         *it_xf   = NULL;
  int64_t       *xoff    = NULL;
  float         *xmx     = NULL;
  int64_t        xalloc  = 0;
  ESL_STOPWATCH *w       = esl_stopwatch_Create();
  char           errbuf[256];
  double         t_dev = 0.0, t_cpu = 0.0;
  int64_t        nseq = 0, nres = 0, ncell = 0;
  int64_t        n_xmxdiff = 0, n_scdiff = 0, n_stdiff = 0, n_fallback = 0;
  int64_t        n_bxmxdiff = 0, n_bownstate = 0, n_btotdiff = 0, n_bfallback = 0;
  int            maxL = 0;
  int            status;

  impl_Init();

  if (p7_pz_DeviceCount() <= 0) esl_fatal("no PEZY device available");

  if (p7_hmmfile_Open(hmmfile, NULL, &hfp, errbuf) != eslOK) esl_fatal("failed to open %s: %s", hmmfile, errbuf);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)             != eslOK) esl_fatal("failed to read HMM from %s", hmmfile);

  bg = p7_bg_Create(abc);
  gm = p7_profile_Create(hmm->M, abc);
  om = p7_oprofile_Create(hmm->M, abc);
  p7_ProfileConfig(hmm, bg, gm, cfgL, p7_LOCAL);
  p7_oprofile_Convert(gm, om);

  if ((pz = p7_pz_Create(devidx, errbuf)) == NULL) esl_fatal("p7_pz_Create failed: %s", errbuf);
  if (! p7_pz_FBModelIsSupported(om))
    esl_fatal("model M=%d: the DP row does not fit in device local memory", om->M);

  if (esl_sqfile_OpenDigital(abc, seqfile, eslSQFILE_FASTA, NULL, &sqfp) != eslOK) esl_fatal("failed to open %s", seqfile);
  sq = esl_sq_CreateDigital(abc);

  b      = p7_pz_seqbatch_Create(batchsz, (int64_t) batchsz * 400);
  it_xf  = malloc(sizeof(float)   * (size_t) batchsz * P7_PZ_XF_N);
  xoff   = malloc(sizeof(int64_t) * (size_t) batchsz);
  if (!b || !it_xf || !xoff) esl_fatal("malloc failed");

  printf("# model %s  M=%d  Kp=%d  Q=%d  row %zu B/thread\n",
         hmm->name, om->M, abc->Kp, P7_PZ_NQF(om->M),
         sizeof(float) * (size_t) P7_PZ_FBROW_FLOATS(P7_PZ_NQF(om->M)));
  printf("# device batch size %d\n", batchsz);

  {
    int eof = FALSE;
    while (! eof)
      {
        int64_t xn = 0;
        int     i;

        p7_pz_seqbatch_Reuse(b);
        while (b->n < batchsz)
          {
            if (maxseq && nseq + b->n >= maxseq) break;
            if ((status = esl_sqio_Read(sqfp, sq)) != eslOK) { eof = TRUE; break; }
            if (sq->n > 0 && sq->n <= 100000 && (capL == 0 || (int) sq->n <= capL))
              p7_pz_seqbatch_Add(b, sq->dsq, (int) sq->n);
            esl_sq_Reuse(sq);
          }
        if (b->n == 0) break;

        /* xmx block offsets, and the per-item specials.  om->xf depends on the
         * target length, exactly as the pipeline's ReconfigLength makes it. */
        for (i = 0; i < b->n; i++)
          {
            xoff[i] = xn;
            xn     += (int64_t) (b->len[i] + 1) * P7_PZ_NXCELLS;
            if (b->len[i] > maxL) maxL = b->len[i];
            p7_oprofile_ReconfigLength(om, b->len[i]);
            p7_pz_FBPackXf(om, it_xf + (size_t) i * P7_PZ_XF_N);
          }
        if (xn * 2 > xalloc)
          { free(xmx); xmx = malloc(sizeof(float) * (size_t) xn * 2); xalloc = xn * 2;
            if (!xmx) esl_fatal("malloc(xmx, %lld floats) failed", (long long) xn * 2); }

        /* ---- device ---- */
        esl_stopwatch_Start(w);
        if ((status = p7_pz_FBParserBatch(pz, devidx, om, b, it_xf, xoff, xmx, xn, errbuf)) != eslOK)
          esl_fatal("p7_pz_FBParserBatch failed (%d): %s", status, errbuf);
        esl_stopwatch_Stop(w);
        t_dev += w->elapsed;

        /* ---- reference ---- */
        if (ox  == NULL) ox  = p7_omx_Create(om->M, 0, maxL);
        if (oxb == NULL) oxb = p7_omx_Create(om->M, 0, maxL);
        p7_omx_GrowTo(ox,  om->M, 0, maxL);
        p7_omx_GrowTo(oxb, om->M, 0, maxL);

        esl_stopwatch_Start(w);
        for (i = 0; i < b->n; i++)
          {
            const ESL_DSQ *dsq = (const ESL_DSQ *) (b->dsq + b->off[i]);
            int            L   = b->len[i];
            const float   *xf  = it_xf + (size_t) i * P7_PZ_XF_N;
            const float   *dev = xmx + xoff[i];
            const float   *bdev = xmx + P7_PZ_BXOFF(xoff[i], xn);
            float          sc_ref, sc_dev;
            int            st_ref, st_dev;

            p7_oprofile_ReconfigLength(om, L);
            st_ref = p7_ForwardParser(dsq, L, om, ox, &sc_ref);
            st_dev = p7_pz_FwdParserScore(dev, L, xf, &sc_dev, NULL);

            /* The device screens itself out when the host's flush-to-zero
             * could have diverged; those targets are the host's to redo, so
             * there is nothing to compare.  Counted, because how often it
             * happens is the cost of the whole scheme. */
            if (st_dev == eslENORESULT) { n_fallback++; nres += L; ncell += (int64_t) L * om->M; continue; }

            /* the whole xmx array, as raw bits */
            if (memcmp(ox->xmx, dev, sizeof(float) * (size_t) (L + 1) * P7_PZ_NXCELLS) != 0)
              {
                if (n_xmxdiff < 8)
                  {
                    int r, s;
                    printf("XMX mismatch: seq %lld L=%d\n", (long long)(nseq + i), L);
                    for (r = 0; r <= L; r++)
                      for (s = 0; s < P7_PZ_NXCELLS; s++)
                        {
                          float a = ox->xmx[r * P7_PZ_NXCELLS + s], c = dev[r * P7_PZ_NXCELLS + s];
                          if (memcmp(&a, &c, sizeof(float)) != 0)
                            { printf("   first at i=%d s=%d  ref %.9g (%08x)  dev %.9g (%08x)\n",
                                     r, s, a, *(unsigned *)&a, c, *(unsigned *)&c);
                              r = L + 1; break; }
                        }
                  }
                n_xmxdiff++;
              }

            if (st_ref != st_dev)
              { if (n_stdiff < 8) printf("STATUS mismatch: seq %lld L=%d ref %d dev %d\n",
                                         (long long)(nseq + i), L, st_ref, st_dev);
                n_stdiff++; }
            else if (st_ref == eslOK && memcmp(&sc_ref, &sc_dev, sizeof(float)) != 0)
              { if (n_scdiff < 8) printf("SCORE mismatch: seq %lld L=%d ref %.9g (%08x) dev %.9g (%08x)\n",
                                         (long long)(nseq + i), L, sc_ref, *(unsigned *)&sc_ref,
                                         sc_dev, *(unsigned *)&sc_dev);
                n_scdiff++; }

            /* ---- Backward, from the same launch ---- */
            {
              float btot;
              int   own, bst;

              bst = p7_pz_BckParserFinish(dev, bdev, L, &own, &btot);
              if (bst == eslENORESULT) n_bfallback++;
              else
                {
                  p7_omx_GrowTo(oxb, om->M, 0, L);
                  p7_BackwardParser(dsq, L, om, ox, oxb, NULL);
                  if (memcmp(oxb->xmx, bdev, sizeof(float) * (size_t) (L + 1) * P7_PZ_NXCELLS) != 0)
                    {
                      if (n_bxmxdiff < 8)
                        {
                          int r, t;
                          printf("BXMX mismatch: seq %lld L=%d\n", (long long)(nseq + i), L);
                          for (r = L; r >= 0; r--)
                            for (t = 0; t < P7_PZ_NXCELLS; t++)
                              {
                                float a = oxb->xmx[r * P7_PZ_NXCELLS + t], c = bdev[r * P7_PZ_NXCELLS + t];
                                if (memcmp(&a, &c, sizeof(float)) != 0)
                                  { printf("   first at i=%d s=%d  ref %.9g (%08x)  dev %.9g (%08x)\n",
                                           r, t, a, *(unsigned *)&a, c, *(unsigned *)&c);
                                    r = -1; break; }
                              }
                        }
                      n_bxmxdiff++;
                    }
                  if (oxb->has_own_scales != own)
                    { if (n_bownstate < 8) printf("OWNSCALES mismatch: seq %lld L=%d ref %d dev %d\n",
                                                  (long long)(nseq + i), L, oxb->has_own_scales, own);
                      n_bownstate++; }
                  if (memcmp(&oxb->totscale, &btot, sizeof(float)) != 0)
                    { if (n_btotdiff < 8) printf("BTOTSCALE mismatch: seq %lld L=%d ref %.9g dev %.9g\n",
                                                 (long long)(nseq + i), L, oxb->totscale, btot);
                      n_btotdiff++; }
                }
            }

            nres  += L;
            ncell += (int64_t) L * om->M;
          }
        esl_stopwatch_Stop(w);
        t_cpu += w->elapsed;

        nseq += b->n;
        if (maxseq && nseq >= maxseq) break;
      }
  }

  printf("# sequences        %" PRId64 "\n", nseq);
  printf("# residues         %" PRId64 "\n", nres);
  printf("# DP cells (M*L)   %" PRId64 "\n", ncell);
  printf("# device time      %.3f s   (%.2f Gcells/s)\n", t_dev, (double) ncell / t_dev / 1e9);
  printf("# 1-core cpu time  %.3f s   (%.2f Gcells/s)  -> device %.1fx one core\n",
         t_cpu, (double) ncell / t_cpu / 1e9, t_cpu / t_dev);
  printf("# host fallbacks   %" PRId64 "  (%.2f%%: a denormal the host would have flushed)\n",
         n_fallback, 100.0 * (double) n_fallback / (double) (nseq ? nseq : 1));
  printf("# xmx diffs        %" PRId64 "\n", n_xmxdiff);
  printf("# status diffs     %" PRId64 "\n", n_stdiff);
  printf("# score diffs      %" PRId64 "\n", n_scdiff);
  printf("# bck fallbacks    %" PRId64 "\n", n_bfallback);
  printf("# bck xmx diffs    %" PRId64 "\n", n_bxmxdiff);
  printf("# bck ownscales    %" PRId64 "\n", n_bownstate);
  printf("# bck totscale     %" PRId64 "\n", n_btotdiff);
  printf("%s\n", (n_xmxdiff == 0 && n_stdiff == 0 && n_scdiff == 0 &&
                   n_bxmxdiff == 0 && n_bownstate == 0 && n_btotdiff == 0)
                 ? "ok. PEZY parser-mode Forward and Backward match impl_sse bit for bit."
                 : "FAILED");

  free(xmx); free(xoff); free(it_xf);
  if (ox)  p7_omx_Destroy(ox);
  if (oxb) p7_omx_Destroy(oxb);
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
  return (n_xmxdiff == 0 && n_stdiff == 0 && n_scdiff == 0 &&
          n_bxmxdiff == 0 && n_bownstate == 0 && n_btotdiff == 0) ? 0 : 1;
}
#endif /*HAVE_PEZY*/
