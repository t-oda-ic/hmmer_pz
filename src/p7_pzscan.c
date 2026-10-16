/* p7_pzscan.c  --  see p7_pzscan.h */
#include <p7_config.h>

#ifdef HAVE_PEZY

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_gumbel.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"
#include "p7_pz.h"
#include "p7_pzdb.h"
#include "p7_pzscan.h"
#include "p7_pz_fbparse.h"

/* HMMER_PEZY_VERBOSE=1 prints a one-line breakdown of the setup and of each
 * query, on stderr.  Purely diagnostic. */
#include <sys/time.h>
static double pz_now(void) { struct timeval t; gettimeofday(&t,NULL); return t.tv_sec + 1e-6*t.tv_usec; }
/* Threads used for the one-off database read.  The default matches HMMER's own
 * thread budget -- one reader plus <--cpu> workers -- since that is what the
 * unaccelerated program would have had running at this point. */
static int pz_rthreads(int nworker)
{
  const char *v = getenv("HMMER_PEZY_RTHREADS");
  if (v) { int n = atoi(v); if (n > 0) return n; }
  return nworker + 1;
}
static int pz_verbose(void) { const char *v = getenv("HMMER_PEZY_VERBOSE"); return (v && v[0] != '0'); }

static int pz_fb_join(P7_PZSCAN *sc);   /* the finish workers call this; defined below */

/* Sequences handed to a worker per atomic grab.
 *
 * The cheap per-target work (length reconfig, NullOne, the SSV outcome) is
 * ~2 us, so even a chunk of 4 is plenty to amortize one atomic add; the cost
 * that has to be averaged out is the ~2% of targets that pass F1, which are
 * orders of magnitude more expensive and unevenly spread.  So aim for many
 * chunks per worker and cap the chunk small.  Measured on Swiss-Prot x CLPX
 * (M=420, --cpu 128), pipeline phase:
 *
 *   chunk    4     8    16    32    64   128   1000
 *   sec   0.547 0.546 0.566 0.564 0.794 1.251  3.316
 *
 * 1000 is HMMER's own block size and reproduces the old static round-robin
 * time exactly -- that is the factor of 13 this replaces (see p7_pzscan.h).
 */
#define PZ_CHUNK_MAX      32
#define PZ_CHUNK_PERWORKER 128   /* target this many chunks per worker */
static int pz_chunk(int nseq, int nworker)
{
  int n;
  const char *v = getenv("HMMER_PEZY_CHUNK");
  if (v) { n = atoi(v); if (n > 0) return n; }

  n = nseq / (nworker * PZ_CHUNK_PERWORKER);
  if (n < 1)             n = 1;
  if (n > PZ_CHUNK_MAX)  n = PZ_CHUNK_MAX;
  return n;
}

int
p7_pzscan_Create(const char *seqfile, ESL_ALPHABET *abc, int nworker, int blocksize,
                 P7_PZSCAN **ret_sc, char *errbuf)
{
  P7_PZSCAN *sc = NULL;
  int        status;

  *ret_sc = NULL;
  if (nworker < 1) nworker = 1;
  (void) blocksize;   /* HMMER's block size is not what we schedule by any more */

  if (p7_pz_DeviceCount() <= 0)
    { if (errbuf) sprintf(errbuf, "no PEZY device available"); return eslENOTFOUND; }

  ESL_ALLOC(sc, sizeof(P7_PZSCAN));
  memset(sc, 0, sizeof(P7_PZSCAN));
  sc->nworker = nworker;

  { double t0 = pz_now(), t1, t2;
    if ((status = p7_pzdb_Build(seqfile, abc, p7_pzdb_MaxBytes(), pz_rthreads(nworker), &sc->db, errbuf)) != eslOK) goto ERROR;
    t1 = pz_now();
    if ((sc->pz = p7_pz_Create(0, errbuf)) == NULL) { status = eslENOTFOUND; goto ERROR; }
    t2 = pz_now();
    if (pz_verbose())
      fprintf(stderr, "# pezy: db build %.3fs (%d seq, %lld res), device open %.3fs\n",
              t1-t0, sc->db->nseq, (long long) sc->db->nres, t2-t1);
  }

  sc->chunk = pz_chunk(sc->db->nseq, nworker);   /* needs nseq, so after the db build */
  if (pz_verbose())
    fprintf(stderr, "# pezy: %d worker%s, chunk %d\n", nworker, nworker == 1 ? "" : "s", sc->chunk);

  { double _t = pz_now();
  sc->batch = p7_pz_seqbatch_Create(16, 4096);   /* immediately re-pointed at the index */
  if (sc->batch == NULL) { status = eslEMEM; goto ERROR; }
  if ((status = p7_pzdb_ToBatch(sc->db, sc->batch)) != eslOK) goto ERROR;

  ESL_ALLOC(sc->xE, sizeof(uint8_t) * sc->db->nseq);

  /* Every target could in principle reach Forward, so size the candidate array
   * for the whole database once rather than growing it under an atomic. */
  ESL_ALLOC(sc->cand, sizeof(P7_PZ_CAND) * (size_t) sc->db->nseq);
  pthread_mutex_init(&sc->fb_mtx, NULL);
  sc->fb_joined = TRUE;                  /* nothing launched yet */
  if ((sc->fb = p7_pz_seqbatch_Create(4096, 4096L * 400)) == NULL) { status = eslEMEM; goto ERROR; }

  if (pz_verbose()) fprintf(stderr, "# pezy: post-db setup %.3fs (cand %.0f MB)\n",
                            pz_now() - _t, (double) sizeof(P7_PZ_CAND) * sc->db->nseq / 1048576.0);
  }

  *ret_sc = sc;
  return eslOK;

 ERROR:
  p7_pzscan_Destroy(sc);
  return status;
}

void
p7_pzscan_Destroy(P7_PZSCAN *sc)
{
  double t0 = pz_now();

  if (!sc) return;
  if (sc->xE)    free(sc->xE);
  if (sc->cand)  free(sc->cand);
  if (sc->it_xf) free(sc->it_xf);
  if (sc->xoff)  free(sc->xoff);
  if (sc->xmx)   free(sc->xmx);
  if (sc->fb)    p7_pz_seqbatch_Destroy(sc->fb);
  if (sc->batch) p7_pz_seqbatch_Destroy(sc->batch);
  if (sc->pz)    p7_pz_Destroy(sc->pz);
  if (sc->db)    p7_pzdb_Destroy(sc->db);
  free(sc);
  if (pz_verbose()) fprintf(stderr, "# pezy: teardown %.3fs\n", pz_now() - t0);
}

int
p7_pzscan_NSeq(const P7_PZSCAN *sc) { return sc->db->nseq; }


/* Phase 1: every target, from the SSV outcome to just past the Viterbi filter.
 *
 * Everything before Forward needs only the digital sequence and its length,
 * both of which P7_PZDB already holds, so nothing is read off disk here -- not
 * even for the targets that pass.  The ~2% that reach Forward are recorded in
 * the shared candidate array; the ESL_SQ they need for their name and their
 * domain alignments is fetched in phase 3, for the ~1% that get that far.
 */
static void *
pzscan_filter_worker(void *arg)
{
  P7_PZ_WORKER  *w   = (P7_PZ_WORKER *) arg;
  const P7_PZDB *db  = w->db;
  P7_PIPELINE   *pli = w->pli;
  P7_OPROFILE   *om  = w->om;
  P7_BG         *bg  = w->bg;
  ESL_SQ         sqv;          /* a borrowed view: dsq points into the db */
  long           n_msvfallback = 0, n_cand = 0;

  /* Same one-time per-thread setup the ordinary worker does.  This is not
   * optional: impl_Init() turns on SSE flush-to-zero, which changes float
   * results for denormals, so every thread that scores sequences must do it. */
  impl_Init();

  memset(&sqv, 0, sizeof(sqv));
  w->status = eslOK;

  for (;;)
    {
      int blk = __sync_fetch_and_add(w->next, 1);
      int lo  = blk * w->chunk;
      int hi, s;

      if (lo >= db->nseq) break;
      hi = lo + w->chunk;  if (hi > db->nseq) hi = db->nseq;

      for (s = lo; s < hi; s++)
        {
          int    L = db->len[s];
          float  usc, nullsc, filtersc;
          int    st;

          /* the unaccelerated loop counts every record, including empty ones */
          p7_pli_NewSeqLength(pli, L);
          p7_bg_SetLength(bg, L);
          p7_oprofile_ReconfigLength(om, L);
          if (L == 0) { p7_pipeline_Reuse(pli); continue; }   /* p7_Pipeline() skips these */

          sqv.dsq = (ESL_DSQ *) (db->dsq_all + db->dsqoff[s]);
          sqv.n   = L;

          if (p7_pz_SSVOutcome(w->xE[s], om, &usc) == eslENORESULT)
            {
              /* SSV cannot decide (the J state might have been used, or the
               * cost constants leave too little headroom): fall back to the
               * full MSV filter on the host, exactly as p7_MSVFilter() does. */
              n_msvfallback++;
              p7_omx_GrowTo(pli->oxf, om->M, 0, L);
              p7_MSVFilter(sqv.dsq, L, om, pli->oxf, &usc);
            }

          st = p7_Pipeline_ToFwd(pli, om, bg, &sqv, &usc, &nullsc, &filtersc);
          if      (st == eslEOD) { p7_pipeline_Reuse(pli); continue; }   /* filtered out */
          else if (st != eslOK)  { w->status = st; return NULL; }

          { int c = __sync_fetch_and_add(w->ncand, 1);
            w->cand[c].sidx     = s;
            w->cand[c].L        = L;
            w->cand[c].nullsc   = nullsc;
            w->cand[c].filtersc = filtersc;
            w->cand[c].xoff     = -1;
            w->cand[c].bxoff    = -1; }
          n_cand++;

          p7_pipeline_Reuse(pli);
        }
    }
  if (pz_verbose())
    fprintf(stderr, "# pezy: filter worker %d: %ld host MSV fallbacks, %ld past F2\n",
            w->widx, n_msvfallback, n_cand);
  return NULL;
}


/* Phase 3: every candidate, from Forward onwards.
 *
 * <c->xoff> >= 0 means the device computed this target's Forward matrix.  We
 * install it in pli->oxf and hand p7_Pipeline_FromFwd() the score.  Otherwise
 * -- the target was too long for the device batch, or its matrix came back
 * with a value small enough that the host's flush-to-zero could have diverged
 * (p7_pz_FwdParserScore returning something other than eslOK) -- we pass NULL
 * and FromFwd runs p7_ForwardParser() itself, exactly as it always did.
 */
static void *
pzscan_finish_worker(void *arg)
{
  P7_PZ_WORKER  *w   = (P7_PZ_WORKER *) arg;
  const P7_PZDB *db  = w->db;
  P7_PIPELINE   *pli = w->pli;
  P7_OPROFILE   *om  = w->om;
  P7_BG         *bg  = w->bg;

  impl_Init();
  w->status  = eslOK;
  w->n_fbdev = w->n_fbhost = 0;
  w->n_bdev  = w->n_bhost  = 0;

  for (;;)
    {
      int blk = __sync_fetch_and_add(w->next, 1);
      int lo  = blk * w->chunk;
      int hi, i;

      lo += w->cand_lo;
      if (lo >= w->cand_hi) break;
      hi = lo + w->chunk;  if (hi > w->cand_hi) hi = w->cand_hi;

      for (i = lo; i < hi; i++)
        {
          const P7_PZ_CAND *c = &w->cand[i];
          P7_PLI_DEVFB  dev;
          int           have_fwd = FALSE;

          if ((w->status = p7_pzdb_FetchSQ(db, w->sqfp, c->sidx, w->sq)) != eslOK) return NULL;

          p7_bg_SetLength(bg, c->L);
          p7_oprofile_ReconfigLength(om, c->L);

          if (c->xoff >= 0)
            {
              float xf8[P7_PZ_XF_N];

              /* The launches may still be running; whoever gets here first
               * waits for them.  Everyone after sees the flag. */
              if ((w->status = pz_fb_join(w->sc)) != eslOK) return NULL;

              memset(&dev, 0, sizeof(dev));
              dev.fxmx = w->xmx + c->xoff;
              dev.bxmx = w->xmx + c->bxoff;
              p7_pz_FBPackXf(om, xf8);

              /* Only the score is needed here; FromFwd copies the matrices in
               * after the F3 test, which about half of these fail. */
              if (p7_pz_FwdParserScore(dev.fxmx, c->L, xf8, &dev.fwdsc, &dev.ftotscale) == eslOK)
                {
                  have_fwd = TRUE;
                  /* Backward came out of the same launch, but only counts if it
                   * passes its own screen; otherwise FromFwd runs
                   * p7_BackwardParser() over the Forward matrix. */
                  if (p7_pz_BckParserFinish(dev.fxmx, dev.bxmx, c->L, &dev.bown, &dev.btotscale) != eslOK)
                    dev.bxmx = NULL;
                }
            }
          if (have_fwd) { w->n_fbdev++; if (dev.bxmx) w->n_bdev++; else w->n_bhost++; }
          else          { w->n_fbhost++; w->n_bhost++; }

          w->status = p7_Pipeline_FromFwd(pli, om, bg, w->sq, NULL, w->th,
                                          c->nullsc, c->filtersc,
                                          have_fwd ? &dev : NULL);
          if (w->status != eslOK) return NULL;

          esl_sq_Reuse(w->sq);
          p7_pipeline_Reuse(pli);
        }
    }
  if (pz_verbose())
    fprintf(stderr, "# pezy: finish worker %d: fwd %ld dev / %ld host, bck %ld dev / %ld host\n",
            w->widx, w->n_fbdev, w->n_fbhost, w->n_bdev, w->n_bhost);
  return NULL;
}


/* Longest target the device batch will take.
 *
 * A Forward pass cannot be split -- the recursion is sequential in i -- so a
 * launch cannot finish before its longest target does, and that time is dead
 * linear in the length (33.5 us per residue at M=420, i.e. M/12.5e6 s).  So
 * the longest targets stay on the host, and the cap is where we draw the line.
 *
 * The cap that matters is a PERCENTILE of the candidate lengths, not a length.
 * Measured on CLPX (M=420) x Swiss-Prot, --cpu 128, summing the two overlapped
 * phases (host-long + device) and the phase that follows them:
 *
 *   LCAP        0    658   1200   2000   3000
 *   host-long  0.469 0.189 0.076  0.046  0.033
 *   devwait    0.000 0.000 0.000  0.043  0.082
 *   finish     0.000 0.259 0.343  0.364  0.372
 *   total      0.534 0.512 0.484  0.516  0.583   <- best at 1200
 *
 * 1200 is where ~95% of the 17,946 candidates fall below (1000 -> 94.1%,
 * 1500 -> 97.4%).  Below that the host is left with too much of the tail;
 * above it the device's own tail starts showing up as devwait.  A percentile
 * adapts to the length distribution, which is what actually varies between
 * databases, so that is what this computes.
 *
 * The absolute clamp is still there as a backstop: PZ_FB_LCAP_MAX bounds the
 * device tail no matter how long the 95th percentile is, and
 * HMMER_PEZY_FBLCAP overrides the whole thing (0 disables the device path).
 */
#define PZ_FB_LCAP_PCT    90          /* leave the longest 10% to the host     */
#define PZ_FB_LCAP_MIN    256
#define PZ_FB_LCAP_MAX    4096
#define PZ_FB_MINCELLS    2.0e9       /* DP cells below which we stay on the host */

/* Below this much device work, don't launch at all.
 *
 * The offload only pays when there is a lot of it to do.  MEASURED over
 * Swiss-Prot, --cpu 128, wall per query, device on vs off (three runs each):
 *
 *   CLPX, M=420, 17,946 candidates   on 0.89/0.91/0.90   off 0.94/0.94/0.96
 *   12 assorted models, summed       on 6.07/6.12/6.13   off 5.70/5.74/5.78
 *
 * i.e. a 5% win on the one query with a lot of work and a 6.5% loss across a
 * mix of ordinary ones.  What separates them is not the model length -- both
 * sets span it -- but how many targets get past the Viterbi filter, which is
 * what a launch has to amortize its fixed cost (kernel object, eight buffers,
 * the uploads and the readback, on each of four devices) over.
 *
 * HMMER_PEZY_FBMINCELLS overrides it; 0 turns the floor off.
 */
static double
pz_fb_mincells(void)
{
  const char *v = getenv("HMMER_PEZY_FBMINCELLS");
  if (v) { double d = atof(v); if (d >= 0.0) return d; }
  return PZ_FB_MINCELLS;
}

static int
pz_fb_lcap(const P7_PZ_CAND *cand, int n, int M, int ndev)
{
  int  *cnt  = NULL;
  int   maxL = 0, i, v, want, run, lcap;
  const char *e = getenv("HMMER_PEZY_FBLCAP");

  if (e) { int k = atoi(e); if (k >= 0) return k; }
  if (n <= 0 || ndev <= 0) return 0;
  (void) M;

  for (i = 0; i < n; i++) if (cand[i].L > maxL) maxL = cand[i].L;
  if ((cnt = malloc(sizeof(int) * (size_t) (maxL + 2))) == NULL) return PZ_FB_LCAP_MIN;
  memset(cnt, 0, sizeof(int) * (size_t) (maxL + 2));
  for (i = 0; i < n; i++) cnt[cand[i].L]++;

  want = (int) ((double) n * PZ_FB_LCAP_PCT / 100.0);
  run  = 0;  lcap = maxL;
  for (v = 0; v <= maxL; v++) { run += cnt[v]; if (run >= want) { lcap = v; break; } }
  free(cnt);

  if (lcap < PZ_FB_LCAP_MIN) lcap = PZ_FB_LCAP_MIN;
  if (lcap > PZ_FB_LCAP_MAX) lcap = PZ_FB_LCAP_MAX;
  return lcap;
}


/* One device's slice of the Forward batch. */
typedef struct {
  P7_PZSCAN         *sc;
  const P7_OPROFILE *om;
  int                devidx;
  int                i0, i1;      /* items [i0,i1) of sc->fb */
  int                status;
  char               errbuf[256];
} P7_PZ_FBSLICE;

static void *
pz_fb_slice_worker(void *arg)
{
  P7_PZ_FBSLICE  *sl = (P7_PZ_FBSLICE *) arg;
  P7_PZSCAN      *sc = sl->sc;
  P7_PZ_SEQBATCH  sub;
  int64_t        *rel  = NULL;
  int64_t         xbase, xn;
  int             j, n = sl->i1 - sl->i0;

  sl->errbuf[0] = '\0';
  if (n <= 0) { sl->status = eslOK; return NULL; }

  /* A window onto the shared batch.  The dsq buffer is shared and off[] holds
   * offsets into it, so a slice needs no copying there. */
  memset(&sub, 0, sizeof(sub));
  sub.dsq      = sc->fb->dsq;
  sub.off      = sc->fb->off + sl->i0;
  sub.len      = sc->fb->len + sl->i0;
  sub.n        = n;
  sub.nres     = sc->fb->nres;
  sub.borrowed = TRUE;

  /* The xmx window has to be per slice, though.  p7_pz_FwdParserBatch()
   * allocates a device buffer of the size it is given and reads all of it
   * back, so handing every slice the whole shared array would have each
   * device download 124 MB and overwrite the other three slices' results. */
  xbase = sc->xoff[sl->i0];
  xn    = (sl->i1 < sc->fb->n ? sc->xoff[sl->i1] : sc->xmx_n) - xbase;
  if ((rel = malloc(sizeof(int64_t) * (size_t) n)) == NULL) { sl->status = eslEMEM; return NULL; }
  for (j = 0; j < n; j++) rel[j] = sc->xoff[sl->i0 + j] - xbase;

  sl->status = p7_pz_FBParserBatch(sc->pz, sl->devidx, sl->om, &sub,
                                    sc->it_xf + (size_t) sl->i0 * P7_PZ_XF_N,
                                    rel, sc->xmx + 2 * xbase, xn, sl->errbuf);
  free(rel);
  return NULL;
}


/* Wait for the device launches, once.  Called by whichever host worker first
 * reaches a candidate that has a device result; everyone else finds the flag
 * already set.  The mutex is what publishes the launches' writes to xmx to the
 * other workers, as well as guarding the flag. */
static int
pz_fb_join(P7_PZSCAN *sc)
{
  P7_PZ_FBSLICE *sl;
  int            d, status = eslOK;

  if (sc->fb_joined) return eslOK;                 /* fast path, no lock */
  pthread_mutex_lock(&sc->fb_mtx);
  if (! sc->fb_joined)
    {
      for (d = 0; d < sc->fb_nslice; d++) pthread_join(sc->fb_tid[d], NULL);
      sl = (P7_PZ_FBSLICE *) sc->fb_sl;
      for (d = 0; d < sc->fb_nslice; d++)
        if (sl[d].status != eslOK)
          {
            /* A device that would not run is not fatal, and not worth being
             * clever about: drop the whole batch back to the host, which then
             * computes exactly what it always did.  Should not happen, so make
             * it loud. */
            if (pz_verbose())
              fprintf(stderr, "# pezy: fbparse slice %d/%d failed (%d): %s -- falling back to the host\n",
                      d, sc->fb_nslice, sl[d].status, sl[d].errbuf);
            { int i; for (i = 0; i < sc->ncand; i++) sc->cand[i].xoff = -1; }
            break;
          }
      free(sc->fb_tid); sc->fb_tid = NULL;
      free(sc->fb_sl);  sc->fb_sl  = NULL;
      sc->fb_joined = TRUE;
    }
  pthread_mutex_unlock(&sc->fb_mtx);
  return status;
}


/* Phase 2: run Forward for every candidate short enough, on every device. */
static int
pzscan_fwd_batch(P7_PZSCAN *sc, P7_OPROFILE *om, int lcap, int *ret_nhost)
{
  pthread_t *tid = NULL;
  P7_PZ_FBSLICE *sl = NULL;
  int        ndev = p7_pz_NDev(sc->pz);
  int        i, d, status = eslOK;

  *ret_nhost = sc->ncand;      /* until we know better, it is all the host's */

  p7_pz_seqbatch_Reuse(sc->fb);
  sc->xmx_n = 0;

  /* Not enough work to be worth a launch?  Leave it all to the host. */
  { double cells = 0.0;
    for (i = 0; i < sc->ncand; i++)
      if (sc->cand[i].L <= lcap) cells += (double) sc->cand[i].L * (double) om->M;
    if (cells < pz_fb_mincells())
      {
        if (pz_verbose())
          fprintf(stderr, "# pezy: fwd batch skipped: %.2f Gcells under the %.2f Gcell floor\n",
                  cells / 1e9, pz_fb_mincells() / 1e9);
        return eslOK;
      } }

  /* Partition: the LONG ones -- the ones the host must run Forward for itself --
   * go to the FRONT, the device's to the back.  Phase 3 is then a single work
   * queue in that order, so every worker spends the launches' duration on work
   * that needs no device result, and whichever worker first reaches the device
   * range joins the launches (pz_fb_join).  No barrier, so no idle waiting.
   * The order candidates sit in does not affect the output (p7_pzscan.h). */
  { int lo = 0, hi = sc->ncand - 1;
    while (lo <= hi)
      {
        if (sc->cand[lo].L > lcap) lo++;                /* host-only: keep at front */
        else { P7_PZ_CAND t = sc->cand[lo]; sc->cand[lo] = sc->cand[hi]; sc->cand[hi] = t; hi--; }
      }
    *ret_nhost = lo; }

  { int nfb = sc->ncand - *ret_nhost;
    if (nfb == 0) return eslOK;
    if (nfb > sc->fb_alloc)
      {
        float   *nx = realloc(sc->it_xf, sizeof(float)   * (size_t) nfb * P7_PZ_XF_N);
        int64_t *no = realloc(sc->xoff,  sizeof(int64_t) * (size_t) nfb);
        if (nx == NULL || no == NULL) { free(nx); free(no); return eslEMEM; }
        sc->it_xf = nx;  sc->xoff = no;  sc->fb_alloc = nfb;
      } }

  /* Collect the ones that go to the device, in candidate order, and give each
   * an xmx block.  om is reconfigured per target because om->xf depends on the
   * target length -- unlike tfv and rfv, which are the model's. */
  for (i = *ret_nhost; i < sc->ncand; i++)
    {
      const P7_PZ_CAND *c = &sc->cand[i];
      int j;

      if (p7_pz_seqbatch_Add(sc->fb, (const ESL_DSQ *) (sc->db->dsq_all + sc->db->dsqoff[c->sidx]),
                             c->L) != eslOK) break;      /* out of room: the rest go to the host */
      j = sc->fb->n - 1;

      sc->xoff[j] = sc->xmx_n;
      sc->xmx_n  += (int64_t) (c->L + 1) * p7X_NXCELLS;

      p7_oprofile_ReconfigLength(om, c->L);
      p7_pz_FBPackXf(om, sc->it_xf + (size_t) j * P7_PZ_XF_N);
    }

  if (sc->fb->n == 0) return eslOK;

  /* Twice over: the Forward blocks in [0,xmx_n) and the Backward ones in
   * [xmx_n, 2*xmx_n) at the same offsets. */
  if (sc->xmx_n * 2 > sc->xmx_alloc)
    {
      float *nw = realloc(sc->xmx, sizeof(float) * (size_t) sc->xmx_n * 2);
      if (nw == NULL) return eslEMEM;
      sc->xmx = nw;  sc->xmx_alloc = sc->xmx_n * 2;
    }

  ESL_ALLOC(sl,  sizeof(P7_PZ_FBSLICE) * ndev);
  ESL_ALLOC(tid, sizeof(pthread_t)     * ndev);
  /* Each device reads its own output window back, so a slice's blocks have to
   * be contiguous on the host.  Slice d therefore owns
   * [2*xbase_d, 2*xbase_d + 2*xn_d) -- Forward in the first half, Backward in
   * the second -- and since the xbase_d tile [0,xmx_n), the windows tile
   * [0, 2*xmx_n) exactly.  Each candidate carries both of its offsets because
   * they are no longer a fixed distance apart. */
  for (d = 0; d < ndev; d++)
    {
      int64_t xbase, xn_d;
      int     j;

      sl[d].sc     = sc;
      sl[d].om     = om;
      sl[d].devidx = d;
      sl[d].i0     = (int) ((long) sc->fb->n *  d      / ndev);
      sl[d].i1     = (int) ((long) sc->fb->n * (d + 1) / ndev);
      sl[d].status = eslOK;

      if (sl[d].i1 <= sl[d].i0) continue;
      xbase = sc->xoff[sl[d].i0];
      xn_d  = (sl[d].i1 < sc->fb->n ? sc->xoff[sl[d].i1] : sc->xmx_n) - xbase;
      for (j = sl[d].i0; j < sl[d].i1; j++)
        {
          P7_PZ_CAND *c = &sc->cand[*ret_nhost + j];   /* batch item j is candidate nhost+j */
          c->xoff  = 2 * xbase + (sc->xoff[j] - xbase);
          c->bxoff = c->xoff + xn_d;
        }
    }
  /* Every slice goes to a thread, including the first: the caller runs host
   * work while they are busy and joins them lazily (pz_fb_join). */
  for (d = 0; d < ndev; d++)
    if (pthread_create(&tid[d], NULL, pz_fb_slice_worker, &sl[d]) != 0) { status = eslESYS; goto ERROR; }

  if (pz_verbose())
    fprintf(stderr, "# pezy: fwd batch %d/%d candidates on %d dev (lcap %d, %.0f MB xmx)\n",
            sc->fb->n, sc->ncand, ndev, lcap,
            (double) sc->xmx_n * sizeof(float) / 1048576.0);

  sc->fb_tid = tid;  sc->fb_sl = sl;  sc->fb_nslice = ndev;  sc->fb_joined = FALSE;
  return eslOK;

 ERROR:
  if (tid) free(tid);
  if (sl)  free(sl);
  return status;
}


/* Run candidates [lo,hi) through p7_Pipeline_FromFwd() on every worker. */
static int
pzscan_finish_range(P7_PZSCAN *sc, P7_PZ_WORKER *w, pthread_t *tid, int lo, int hi)
{
  int i, status = eslOK;

  if (hi <= lo) return eslOK;
  sc->next = 0;
  for (i = 0; i < sc->nworker; i++)
    {
      w[i].cand_lo = lo;
      w[i].cand_hi = hi;
      w[i].xmx     = sc->xmx;
      w[i].chunk   = 1;          /* candidates are few and wildly uneven */
      w[i].status  = eslOK;
    }
  if (sc->nworker == 1) pzscan_finish_worker(&w[0]);
  else {
    for (i = 0; i < sc->nworker; i++)
      if (pthread_create(&tid[i], NULL, pzscan_finish_worker, &w[i]) != 0) return eslESYS;
    for (i = 0; i < sc->nworker; i++) pthread_join(tid[i], NULL);
  }
  for (i = 0; i < sc->nworker; i++) if (w[i].status != eslOK) status = w[i].status;
  return status;
}


int
p7_pzscan_Run(P7_PZSCAN *sc, P7_PZ_WORKER *w)
{
  pthread_t *tid = NULL;
  int        nhost_cand = 0;
  int        i, lcap;
  int        status;
  double     t0 = pz_now(), t1, t2, t3, t4;

  if ((status = p7_pz_SetModel(sc->pz, w[0].om)) != eslOK) return status;   /* eslENORESULT: too big for the device */
  if ((status = p7_pz_SSVBatch(sc->pz, sc->batch, sc->xE)) != eslOK) return status;
  t1 = pz_now();

  if (sc->nworker > 1) ESL_ALLOC(tid, sizeof(pthread_t) * sc->nworker);

  for (i = 0; i < sc->nworker; i++)
    {
      w[i].db      = sc->db;
      w[i].xE      = sc->xE;
      w[i].widx    = i;
      w[i].nworker = sc->nworker;
      w[i].chunk   = sc->chunk;
      w[i].next    = &sc->next;
      w[i].sc      = sc;
      w[i].cand    = sc->cand;
      w[i].ncand   = &sc->ncand;
      w[i].status  = eslOK;
      if (w[i].sq == NULL && (w[i].sq = esl_sq_CreateDigital(sc->db->abc)) == NULL) { status = eslEMEM; goto ERROR; }
    }

  /* --- phase 1: filters, over every target --- */
  sc->next = 0;  sc->ncand = 0;
  if (sc->nworker == 1) pzscan_filter_worker(&w[0]);
  else {
    for (i = 0; i < sc->nworker; i++)
      if (pthread_create(&tid[i], NULL, pzscan_filter_worker, &w[i]) != 0) { status = eslESYS; goto ERROR; }
    for (i = 0; i < sc->nworker; i++) pthread_join(tid[i], NULL);
  }
  for (i = 0; i < sc->nworker; i++) if (w[i].status != eslOK) { status = w[i].status; goto ERROR; }
  t2 = pz_now();

  /* --- phase 2: Forward+Backward on the devices, for everything short enough ---
   *
   * The launches go to their own threads and nobody waits here.  The candidates
   * the device did NOT take are at the front, so phase 3 -- a single queue over
   * all of them, in that order -- spends the launches' duration on work that
   * needs no device result, and whichever worker first reaches the device range
   * joins them (pz_fb_join).  One queue rather than two ranges with a barrier
   * between: the barrier version measured no better than the synchronous one,
   * because splitting phase 3 in two cost as much in load balance as the
   * overlap gained. */
  lcap = pz_fb_lcap(sc->cand, sc->ncand, w[0].om->M, p7_pz_NDev(sc->pz));
  if (lcap > 0 && p7_pz_FBModelIsSupported(w[0].om))
    {
      if ((status = pzscan_fwd_batch(sc, w[0].om, lcap, &nhost_cand)) != eslOK) goto ERROR;
    }
  t3 = pz_now();

  /* --- phase 3: every candidate, host-only ones first --- */
  status = pzscan_finish_range(sc, w, tid, 0, sc->ncand);
  if (status != eslOK) goto ERROR;
  pz_fb_join(sc);                  /* in case nothing reached the device range */
  t4 = pz_now();

  if (pz_verbose())
    fprintf(stderr, "# pezy: M=%d ssv %.3fs, filters %.3fs, fwd launch %.3fs, finish %.3fs "
                    "(%d host-only + %d device candidates)\n",
            w[0].om->M, t1-t0, t2-t1, t3-t2, t4-t3, nhost_cand, sc->ncand - nhost_cand);

  if (tid) free(tid);
  return eslOK;

 ERROR:
  pz_fb_join(sc);                  /* never leave launch threads unjoined */
  if (tid) free(tid);
  return status;
}

#else /*!HAVE_PEZY*/
int p7_pzscan_dummy_symbol_to_avoid_empty_translation_unit = 0;
#endif /*HAVE_PEZY*/
