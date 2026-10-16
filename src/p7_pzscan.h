/* p7_pzscan.h
 *
 * Run one query model over a resident, pre-digitized target database
 * (P7_PZDB) with the SSV filter on a PEZY-SC device.
 *
 * This replaces the reader-thread / work-queue loop of hmmsearch, phmmer and
 * jackhmmer.  The device scores every target sequence in one pass; the host
 * then only has to materialize the ~2% that survive the F1 threshold and hand
 * them to p7_Pipeline_MSVDone(), which is the ordinary pipeline minus the MSV
 * filter.
 *
 * Sequences are handed to workers in small chunks off a shared atomic counter,
 * i.e. dynamically.  The distribution does not affect the output, which is why
 * we are free to schedule for load:
 *
 *   p7_pli_NewSeqLength() sets pli->Z to the worker's *own* running sequence
 *   count, and p7_pli_TargetReportable() tests exp(lnP) * pli->Z <= pli->E.  A
 *   running count is always <= the final Z, so the in-search test is always
 *   more permissive than the final one, and the hits a worker stores are always
 *   a superset of the hits that survive p7_tophits_Threshold() once the merged
 *   Z is known.  Measured: --cpu 1, 8 and 128 give byte-identical output.
 *
 * Static round-robin over HMMER's 1000-sequence blocks was what this used to
 * do, and it cost a factor of 13: only ~2% of targets pass F1 and the cost of
 * one that does spans orders of magnitude (a multidomain region draws 200
 * stochastic tracebacks), so 4.5 blocks per worker does not average out.  One
 * straggler decided the wall time of the whole phase.
 */
#ifndef p7_PZSCAN_INCLUDED
#define p7_PZSCAN_INCLUDED

#include <p7_config.h>

#ifdef HAVE_PEZY

#include <pthread.h>

#include "easel.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"
#include "p7_pz.h"
#include "p7_pzdb.h"

typedef struct p7_pzscan_s P7_PZSCAN;   /* defined below; the worker needs a pointer */

/* A target that got past the Viterbi filter, i.e. one that Forward has to be
 * run for.  These are what the device batch is built from. */
typedef struct {
  int      sidx;          /* index into the database                          */
  int      L;
  float    nullsc;        /* what p7_Pipeline_ToFwd() returned                */
  float    filtersc;
  int64_t  xoff;          /* offset of this target's Forward xmx block, or -1
                           * if the host is to run Forward for it             */
  int64_t  bxoff;         /* and its Backward one                             */
} P7_PZ_CAND;

typedef struct {
  /* filled in by the caller, one per worker */
  P7_OPROFILE *om;
  P7_BG       *bg;
  P7_PIPELINE *pli;
  P7_TOPHITS  *th;
  ESL_SQFILE  *sqfp;      /* private positionable handle on the database file */

  /* private to the runner */
  ESL_SQ         *sq;
  const P7_PZDB  *db;
  const uint8_t  *xE;
  int             widx, nworker, chunk;
  int            *next;    /* shared chunk counter, bumped with an atomic add */
  int             status;

  /* the two-phase Forward path */
  P7_PZ_CAND     *cand;    /* shared candidate array                          */
  int            *ncand;   /* shared count, bumped with an atomic add         */
  int             cand_lo, cand_hi;    /* the slice of the candidate array this
                                        * pass is to work on                   */
  P7_PZSCAN      *sc;      /* for the lazy join of the device launches       */
  const float    *xmx;     /* shared: the device's matrices; each candidate
                            * carries its own two offsets into it            */
  long            n_fbdev, n_fbhost;   /* diagnostic tallies                  */
  long            n_bdev,  n_bhost;
} P7_PZ_WORKER;

/* One-shot context, created once per process and reused for every query. */
struct p7_pzscan_s {
  P7_PZ          *pz;
  P7_PZDB        *db;
  P7_PZ_SEQBATCH *batch;
  uint8_t        *xE;
  int             nworker;
  int             chunk;    /* sequences handed out per grab */
  int             next;     /* shared chunk counter for the current query */

  /* the two-phase Forward path (see p7_pzscan.c) */
  P7_PZ_CAND     *cand;     /* [nseq] targets that reached Forward           */
  int             ncand;
  P7_PZ_SEQBATCH *fb;       /* the subset that goes to the device            */
  float          *it_xf;    /* [fb->n * P7_PZ_XF_N] per-target specials      */
  int64_t        *xoff;     /* [fb->n] where each one's xmx block lands      */
  float          *xmx;      /* the returned Forward matrices                 */
  int64_t         xmx_n, xmx_alloc;
  int             fb_alloc;

  /* The device launches run in their own threads while the host works on the
   * candidates the device did not take.  Whichever host worker first reaches a
   * candidate that DOES have a device result joins them; the rest just see the
   * flag.  That way there is no barrier between "host-only work" and "device
   * work", so no worker ever sits idle waiting for a launch. */
  pthread_t      *fb_tid;
  void           *fb_sl;      /* P7_PZ_FBSLICE[], opaque here                 */
  int             fb_nslice;
  int             fb_joined;
  pthread_mutex_t fb_mtx;
};

/* Try to set up the accelerated path for <seqfile>.  Returns eslOK with
 * *ret_sc set, or a non-eslOK status (with *ret_sc NULL) if the device or the
 * file is unsuitable -- in which case the caller just uses its normal loop.
 */
extern int  p7_pzscan_Create(const char *seqfile, ESL_ALPHABET *abc, int nworker, int blocksize,
                             P7_PZSCAN **ret_sc, char *errbuf);
extern void p7_pzscan_Destroy(P7_PZSCAN *sc);

/* Number of target sequences, for callers that want to report progress. */
extern int  p7_pzscan_NSeq(const P7_PZSCAN *sc);

/* Run one query.  <w> is an array of <sc->nworker> workers whose om/bg/pli/th
 * the caller has already created.  Returns eslOK, or eslENORESULT if this
 * model cannot go on the device (caller must fall back for this query).
 */
extern int  p7_pzscan_Run(P7_PZSCAN *sc, P7_PZ_WORKER *w);

#endif /*HAVE_PEZY*/
#endif /*p7_PZSCAN_INCLUDED*/
