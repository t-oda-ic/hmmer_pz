/* p7_pz.c
 *
 * PEZY-SC (PZCL) offload of HMMER's SSV filter.
 *
 * Contents:
 *   1. Device discovery and process-wide initialization
 *   2. P7_PZ context (one per worker thread)
 *   3. Model upload (de-striping the SSE SSV score table)
 *   4. Sequence batches
 *   5. Running the filter
 *   6. Host-side tail of p7_SSVFilter()
 */
#include <p7_config.h>

#ifdef HAVE_PEZY

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <sys/time.h>

#include <pzcl/pzcl.h>
#include <pzcl/pzcl_ext.h>

#include "easel.h"
#include "esl_alphabet.h"

#include "hmmer.h"
#include "p7_pz.h"
#include "p7_pz_fpprobe.h"
#include "p7_pz_fbparse.h"

/* The embedded device binary (objcopy'd from p7_pz_kernel.pz). */
extern char _binary_p7_pz_kernel_pz_start;
extern char _binary_p7_pz_kernel_pz_end;

/* Per-thread stack size we request.  What the PE's local memory leaves after
 * eight of these stacks is ours. */
#define P7_PZ_STACKSIZE   1024
#define P7_PZ_LMEM_AVAIL  (24576 - P7_PZ_STACKSIZE * 8)

/* Hard cap: the DP row is per-thread, 8 threads per PE, and is padded up to a
 * multiple of 8 bytes so the SWAR loop can use aligned 64-bit accesses. */
#define P7_PZ_MAXM        (P7_PZ_LMEM_AVAIL / 8 - 8)

/* The parser-Forward kernel has a stack setting of its own, because local
 * memory is what limits how many of a PE's threads can run: each needs a
 * 48*Q-byte DP row out of what the eight stacks leave.  A smaller stack would
 * let more threads fit, but this stays at 1 KB (two threads per PE at M=420).
 * pz_fwdparse_utest is what verifies any change here, since a stack overrun
 * would corrupt the DP row and every result with it. */
#define P7_PZ_FB_STACKSIZE  1024
#define P7_PZ_FB_LMEM_AVAIL (24576 - P7_PZ_FB_STACKSIZE * 8)

/* Keeping a launch shorter than the driver's watchdog
 * ---------------------------------------------------
 * The driver puts a time limit on every kernel launch.  Rather than depend on
 * how a launch that hits it is reported, the host checks for itself that each
 * launch ran to the end, and uses the output only if it did.  Two defences,
 * both required:
 *
 *   1. Split the work into launches of about <P7_PZ_KTIME> seconds each, sized
 *      from a running estimate of the device's cell rate.  counter[0] is
 *      seeded with the first item of the range and <nitem> is passed as the
 *      item one past the last, so no kernel argument had to be added.
 *   2. Verify afterwards that every thread reported completion, in
 *      counter[P7_PZ_CNT_DONE].  Anything else is a failed launch.
 */
#define P7_PZ_CNT_N       64        /* ints in the counter buffer            */
#define P7_PZ_CNT_DONE    16        /* slot the kernel tallies threads in    */
#define P7_PZ_NITEM_ARG    7        /* every kernel takes <nitem> here       */
/* Initial cells/s guess for one card, before any launch has been timed.
 *
 * This sizes the FIRST launch, and for hmmsearch -- one query per process --
 * there is no second query to amortize a bad guess over.  Guessing low is not
 * free: a small launch runs at a much lower rate than a big one, because with
 * few items most of the device's threads get nothing and the time is set by
 * the per-thread work.  MEASURED, one card, Swiss-Prot:
 *
 *   launch size   0.2 Gcells   3.5      16.4     62.8
 *   rate            2.5 G/s    14 G/s   60 G/s   183 G/s
 *
 * and end to end at M=1414 (the largest model in the benchmark), where the old
 * 4.0e10 forced two launches per card:
 *
 *   launches   8            4
 *   ssv        1.082 s      0.747 s      (-31%)
 *   wall       1.49 s       1.19 s       (-20%)
 *
 * So guess near the slowest rate a LARGE launch actually achieves -- 99 G/s
 * per card at M=1414 -- rather than near the rate a small one does.  With
 * P7_PZ_KTIME at 1 s that is a 100 Gcell budget, which at that rate is 1 s,
 * well inside the driver's watchdog limit.  A model slow enough to overrun is
 * caught by the completion tally and retried with half the budget (see the
 * launch loop), so being wrong here costs one launch, not correctness.
 */
#define P7_PZ_CPS0        1.0e11
#define P7_PZ_CPS_MIN     1.0e9
#define P7_PZ_CPS_MAX     2.0e12

static double pz_now(void) { struct timeval t; gettimeofday(&t, NULL); return t.tv_sec + 1e-6 * t.tv_usec; }
static int    pz_verbose(void) { const char *v = getenv("HMMER_PEZY_VERBOSE"); return (v && v[0] != '0'); }

/* Which variant of each kernel to run.
 *
 * A PE's eight threads add no arithmetic throughput of their own; what they
 * buy is latency hiding, and what they cost is cache footprint.  Which way
 * that trade falls depends on whether the eight threads are reading the same
 * score table or eight different ones.
 *
 *   ssv_lmem/gmem   One query model against many sequences: all eight threads
 *                   read the one model's table, so they share it and only hide
 *                   each other's global-load latency.  Eight threads measured
 *                   884 Gcells/s against 345 for one (M=95, Swiss-Prot).
 *                   HMMER_PEZY_TPP=1 runs one thread per PE instead.
 *
 *   ssv_multi       Many models.  Three variants, HMMER_PEZY_MULTI:
 *                     thread  8 threads, each on its own model.  Eight tables
 *                             compete for the same cache and just evict each
 *                             other: 134 Gcells/s.
 *                     pe      1 thread per PE, so one table: 292 Gcells/s.
 *                     share   8 threads all on ONE model, its table copied
 *                             into local memory: both effects at once.  The
 *                             default; needs npq >= 8 to have eight queries of
 *                             one model to hand out, so the <pe> kernel stays
 *                             loaded for smaller batches.
 *                   (Pfam-A x 128 queries.)
 */
#define P7_PZ_MULTI_SHARE  0
#define P7_PZ_MULTI_PE     1
#define P7_PZ_MULTI_THREAD 2

static int
pz_batch_threads(void)
{
  const char *v = getenv("HMMER_PEZY_TPP");
  return (v && atoi(v) == 1) ? 1 : 8;
}

static int
pz_multi_mode(void)
{
  const char *v = getenv("HMMER_PEZY_MULTI");
  if (v && strcmp(v, "pe")     == 0) return P7_PZ_MULTI_PE;
  if (v && strcmp(v, "thread") == 0) return P7_PZ_MULTI_THREAD;
  return P7_PZ_MULTI_SHARE;
}

/* Seconds of device time one launch should aim for.  Well under the driver's
 * watchdog limit, and long enough that per-launch overhead stays invisible. */
static double
pz_ktime(void)
{
  const char *v = getenv("HMMER_PEZY_KTIME");
  double      t = (v ? atof(v) : 1.0);
  if (! (t >= 0.02)) t = 1.0;
  if (t > 3.0)       t = 3.0;
  return t;
}

/* Cells one launch may cover, from the current rate estimate. */
static int64_t
pz_cell_budget(double cps)
{
  double b = cps * pz_ktime();
  if (b < 1.0e7)  b = 1.0e7;
  if (b > 4.0e12) b = 4.0e12;
  return (int64_t) b;
}

/* Fold a launch's measured rate into the estimate: drop to a slower rate at
 * once (that is the safe direction), climb toward a faster one gradually. */
static void
pz_update_cps(double *cps, int64_t cells, double secs)
{
  double r;
  if (secs < 1.0e-3 || cells <= 0) return;     /* too short to measure */
  r = (double) cells / secs;
  *cps = (r < *cps) ? r : (3.0 * *cps + r) / 4.0;
  if (*cps < P7_PZ_CPS_MIN) *cps = P7_PZ_CPS_MIN;
  if (*cps > P7_PZ_CPS_MAX) *cps = P7_PZ_CPS_MAX;
}

/* One device.  Everything that is per-card -- context, queue, program,
 * kernels, buffers, and the rate estimate that sizes its launches -- lives
 * here; P7_PZ owns an array of these plus the host-side staging that all the
 * cards share.
 */
typedef struct {
  int               idx;         /* 0..ndev-1, for diagnostics         */
  pzcl_device_id    dev;
  pzcl_context      ctx;
  pzcl_command_queue q;
  pzcl_program      prog;
  pzcl_kernel       k_lmem;      /* score tables in local memory      */
  pzcl_kernel       k_gmem;      /* score tables in global memory      */
  pzcl_kernel       k_multi;     /* many models x many queries         */
  pzcl_kernel       k_mshare;    /* one model per PE, table in lmem    */
  size_t            gws;

  /* device buffers, grown on demand */
  pzcl_mem  d_sb;    size_t d_sb_sz;
  pzcl_mem  d_dsq;   size_t d_dsq_sz;
  pzcl_mem  d_qoff;  size_t d_qoff_sz;
  pzcl_mem  d_a;     size_t d_a_sz;
  pzcl_mem  d_b;     size_t d_b_sz;
  pzcl_mem  d_xE;    size_t d_xE_sz;
  pzcl_mem  d_cnt;   size_t d_cnt_sz;

  /* This card's slice of the work: items [it0,it1) of the shared list, which
   * cover the <dsq> bytes [dsq_base, dsq_base+dsq_len).  <h_qoff> is the
   * card's own copy of the item offsets, rebased to its slice of <dsq>. */
  int       seq0, seq1;   /* sequences [seq0,seq1) of the batch      */
  int       it0, it1;     /* the items that cover them               */
  int64_t   dsq_base, dsq_len;
  int64_t  *h_qoff;   int qoff_alloc;

  /* what is already sitting on this card, so repeated passes over the same
   * database do not re-upload a few hundred megabytes per query */
  int      up_n;        /* sequences uploaded (0 = nothing)          */
  int64_t  up_nres;
  int64_t  up_base, up_len;   /* the <dsq> slice this card holds       */
  int      up_chunk;    /* chunk length the item list was built for   */
  int      up_it0, up_it1;

  double   cps;         /* running estimate of this card's cells/s */
} P7_PZDEV;

struct p7_pz_s {
  int               ndev;        /* cards this context drives            */
  P7_PZDEV         *dv;          /* [ndev]                               */
  int               tpp_batch;   /* threads per PE: ssv_lmem / ssv_gmem  */
  int               multi_mode;  /* P7_PZ_MULTI_*                        */

  /* host-side work-item staging, grown on demand and reused.  One copy for
   * the whole context: the cards divide the item list, they do not each get
   * their own. */
  int64_t  *h_qoff;
  int32_t  *h_a;
  int32_t  *h_b;
  int32_t  *h_seq;      /* item -> sequence, host only (for the reduction) */
  uint8_t  *h_part;
  int       nitem, item_alloc;
  int       item_chunk;     /* chunk length the item list holds        */
  int       item_n;         /* and the batch it was built from         */
  int64_t   item_nres;
  int       part_ndev;      /* cards the current partition was cut for  */
  int       part_n;         /* the batch it was cut for: sequences...   */
  int64_t   part_nres;      /* ...and residues.  Both: a second batch   */
                            /* of the same length would otherwise reuse */
                            /* the first one's byte offsets.            */
  int       part_nitem;     /* items the item-slices were mapped for    */

  /* current model */
  int   M;
  int   Kp;
  int   use_lmem;   /* TRUE if the sb table fits in local memory */

  /* Work handed out to the cards.  A card takes a range of items, sized from
   * its own rate estimate, runs it, and comes back for more; that balances
   * the cards against each other without assuming they are equally fast or
   * equally free.  <mx> guards <cursor> and <failed>. */
  pthread_mutex_t mx;
  int             cursor;
  int             failed;   /* first non-eslOK status any card returned  */
};

/*****************************************************************
 * 1. Device discovery and process-wide initialization
 *****************************************************************/

static pthread_once_t  pz_once      = PTHREAD_ONCE_INIT;
static int             pz_ndev      = 0;
static pzcl_device_id  pz_devs[16];
static int             pz_enabled    = TRUE;
static pfnPezyExtSetPerThreadStackSize pz_set_stack = NULL;

static void
pz_discover(void)
{
  pzcl_platform_id plat;
  pzcl_uint        np = 0, nd = 0;
  pzcl_device_id   devs[64];
  const char      *s;

  if ((s = getenv("HMMER_PEZY")) != NULL && (s[0] == '0' || s[0] == 'n' || s[0] == 'N'))
    { pz_enabled = FALSE; return; }

  if (pzclGetPlatformIDs(1, &plat, &np) != PZCL_SUCCESS || np == 0) return;
  if (pzclGetDeviceIDs(plat, PZCL_DEVICE_TYPE_DEFAULT, 64, devs, &nd) != PZCL_SUCCESS || nd == 0) return;

  if (nd > 16) nd = 16;
  for (pzcl_uint i = 0; i < nd; i++) pz_devs[i] = devs[i];
  pz_ndev = (int) nd;

  pz_set_stack = (pfnPezyExtSetPerThreadStackSize)
                 pzclGetExtensionFunctionAddress("pezy_set_per_thread_stack_size");
}

int
p7_pz_DeviceCount(void)
{
  pthread_once(&pz_once, pz_discover);
  return pz_enabled ? pz_ndev : 0;
}

int  p7_pz_IsEnabled(void)      { return pz_enabled && p7_pz_DeviceCount() > 0; }
void p7_pz_SetEnabled(int on)   { pz_enabled = on; }

/*****************************************************************
 * 2. P7_PZ context
 *****************************************************************/

static void
pz_free_mem(pzcl_mem *m, size_t *sz)
{
  if (*m) { pzclReleaseMemObject(*m); *m = NULL; }
  *sz = 0;
}

/* Grow a device buffer to at least <need> bytes.  Buffers are only ever
 * enlarged, so steady-state operation does no allocation. */
static int
pz_grow(pzcl_context ctx, pzcl_mem *m, size_t *sz, size_t need)
{
  pzcl_int err;
  if (need == 0) return eslOK;
  if (*m && *sz >= need) return eslOK;
  if (*m) pzclReleaseMemObject(*m);
  need += need / 4;                     /* 25% headroom to damp regrowth */
  *m = pzclCreateBuffer(ctx, PZCL_MEM_READ_WRITE, need, NULL, &err);
  if (err != PZCL_SUCCESS) { *m = NULL; *sz = 0; return eslEMEM; }
  *sz = need;
  return eslOK;
}

/* How many cards to drive.  All of them by default; HMMER_PEZY_NDEV caps it,
 * which is how the single-card path is measured for comparison. */
static int
pz_want_ndev(void)
{
  const char *v = getenv("HMMER_PEZY_NDEV");
  int         n = p7_pz_DeviceCount();
  if (v) { int w = atoi(v); if (w >= 1 && w < n) n = w; }
  return n;
}

/* Bring up one card: context, queue, program, kernels. */
static int
pz_dev_open(P7_PZDEV *dd, pzcl_device_id dev, int idx, int tpp_batch, int multi_mode, char *errbuf)
{
  pzcl_int   err;
  size_t     binsz;
  const unsigned char *bin;

  memset(dd, 0, sizeof(P7_PZDEV));
  dd->idx = idx;
  dd->dev = dev;

  dd->ctx = pzclCreateContext(NULL, 1, &dd->dev, NULL, NULL, &err);
  if (err != PZCL_SUCCESS) { if (errbuf) sprintf(errbuf, "pzclCreateContext(dev %d) failed (%d)", idx, (int) err); return eslFAIL; }

  dd->q = pzclCreateCommandQueue(dd->ctx, dd->dev, 0, &err);
  if (err != PZCL_SUCCESS) { if (errbuf) sprintf(errbuf, "pzclCreateCommandQueue(dev %d) failed (%d)", idx, (int) err); return eslFAIL; }

  binsz = (size_t) (&_binary_p7_pz_kernel_pz_end - &_binary_p7_pz_kernel_pz_start);
  bin   = (const unsigned char *) &_binary_p7_pz_kernel_pz_start;
  /* binary_status is not needed: success is judged from errcode_ret alone. */
  dd->prog = pzclCreateProgramWithBinary(dd->ctx, 1, &dd->dev, &binsz, &bin, NULL, &err);
  if (err != PZCL_SUCCESS) { if (errbuf) sprintf(errbuf, "pzclCreateProgramWithBinary(dev %d) failed (%d)", idx, (int) err); return eslFAIL; }

  /* Kernel names drop the pzc_ prefix.  With one thread per PE we take the
   * _pe variants, in which the other seven threads leave immediately. */
  {
    const char *bsfx  = (tpp_batch == 1) ? "_pe" : "";
    const char *msfx  = (multi_mode == P7_PZ_MULTI_THREAD) ? "" : "_pe";
    char        nm[32];
    int         i;
    struct { pzcl_kernel *slot; const char *base, *sfx; } want[4];
    int         nwant = 3;

    want[0].slot = &dd->k_lmem;   want[0].base = "ssv_lmem";  want[0].sfx = bsfx;
    want[1].slot = &dd->k_gmem;   want[1].base = "ssv_gmem";  want[1].sfx = bsfx;
    want[2].slot = &dd->k_multi;  want[2].base = "ssv_multi"; want[2].sfx = msfx;
    if (multi_mode == P7_PZ_MULTI_SHARE)
      { want[3].slot = &dd->k_mshare; want[3].base = "ssv_multi"; want[3].sfx = "_share"; nwant = 4; }

    for (i = 0; i < nwant; i++)
      {
        snprintf(nm, sizeof(nm), "%s%s", want[i].base, want[i].sfx);
        *want[i].slot = pzclCreateKernel(dd->prog, nm, &err);
        if (err != PZCL_SUCCESS)
          { if (errbuf) sprintf(errbuf, "pzclCreateKernel(%s, dev %d) failed (%d)", nm, idx, (int) err); return eslFAIL; }
        if (pz_set_stack && pz_set_stack(*want[i].slot, P7_PZ_STACKSIZE) != PZCL_SUCCESS)
          { if (errbuf) sprintf(errbuf, "pezy_set_per_thread_stack_size(%s, dev %d) failed", nm, idx); return eslFAIL; }
      }
  }

  dd->gws = 0;
  if (pzclGetDeviceInfo(dd->dev, PZCL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(dd->gws), &dd->gws, NULL) != PZCL_SUCCESS
      || dd->gws == 0)
    { if (errbuf) sprintf(errbuf, "could not query MAX_WORK_GROUP_SIZE (dev %d)", idx); return eslFAIL; }

  dd->up_n = 0; dd->up_nres = 0; dd->up_base = dd->up_len = -1;
  dd->up_chunk = 0; dd->up_it0 = dd->up_it1 = -1;
  dd->cps  = P7_PZ_CPS0;
  return eslOK;
}

/* one card's share of bringing the context up */
typedef struct {
  P7_PZDEV      *dd;
  pzcl_device_id dev;
  int            idx, tpp_batch, multi_mode;
  char           why[eslERRBUFSIZE];
  int            status;
} pz_open_t;

static void *
pz_open_one(void *arg)
{
  pz_open_t *op = (pz_open_t *) arg;
  op->status = pz_dev_open(op->dd, op->dev, op->idx, op->tpp_batch, op->multi_mode, op->why);
  return NULL;
}

static void
pz_dev_close(P7_PZDEV *dd)
{
  pz_free_mem(&dd->d_sb,   &dd->d_sb_sz);
  pz_free_mem(&dd->d_dsq,  &dd->d_dsq_sz);
  pz_free_mem(&dd->d_qoff, &dd->d_qoff_sz);
  pz_free_mem(&dd->d_a,    &dd->d_a_sz);
  pz_free_mem(&dd->d_b,    &dd->d_b_sz);
  pz_free_mem(&dd->d_xE,   &dd->d_xE_sz);
  pz_free_mem(&dd->d_cnt,  &dd->d_cnt_sz);
  if (dd->h_qoff)   free(dd->h_qoff);
  if (dd->k_lmem)   pzclReleaseKernel(dd->k_lmem);
  if (dd->k_gmem)   pzclReleaseKernel(dd->k_gmem);
  if (dd->k_multi)  pzclReleaseKernel(dd->k_multi);
  if (dd->k_mshare) pzclReleaseKernel(dd->k_mshare);
  if (dd->prog)     pzclReleaseProgram(dd->prog);
  if (dd->q)        pzclReleaseCommandQueue(dd->q);
  if (dd->ctx)      pzclReleaseContext(dd->ctx);
  memset(dd, 0, sizeof(P7_PZDEV));
}

/* Create a context over every visible card, starting at <devidx>.
 *
 * <devidx> survives from the single-card interface and now only rotates which
 * physical card is card 0 of this context; every caller in HMMER passes 0.
 */
P7_PZ *
p7_pz_Create(int devidx, char *errbuf)
{
  P7_PZ     *pz = NULL;
  int        ndev = p7_pz_DeviceCount();
  int        i;
  int        status;

  if (ndev <= 0) { if (errbuf) sprintf(errbuf, "no PEZY device available"); return NULL; }

  if (errbuf) errbuf[0] = '\0';

  ESL_ALLOC(pz, sizeof(P7_PZ));
  memset(pz, 0, sizeof(P7_PZ));
  pthread_mutex_init(&pz->mx, NULL);

  pz->tpp_batch  = pz_batch_threads();
  pz->multi_mode = pz_multi_mode();
  pz->ndev       = pz_want_ndev();

  ESL_ALLOC(pz->dv, sizeof(P7_PZDEV) * pz->ndev);
  memset(pz->dv, 0, sizeof(P7_PZDEV) * pz->ndev);

  /* Bring the cards up, one thread each.  Opening a card is about 20 ms --
   * mostly handing the kernel binary to the driver -- and hmmsearch pays it
   * once per process for a search that itself takes a fraction of a second,
   * so four of them in series would be most of what the extra cards win.
   *
   * A card that cannot be opened (one another process is holding, say) is
   * skipped, and we run on what we can get: dropping to three cards is a
   * slowdown, but refusing to start would be a fallback to the CPU for the
   * whole run.
   */
  {
    pz_open_t *op  = NULL;
    pthread_t *tid = NULL;
    int        want = pz->ndev, opened = 0;

    ESL_ALLOC(op, sizeof(pz_open_t) * want);
    for (i = 0; i < want; i++)
      {
        op[i].dd  = &pz->dv[i];
        op[i].dev = pz_devs[((devidx + i) % ndev + ndev) % ndev];
        op[i].idx = i;
        op[i].tpp_batch = pz->tpp_batch; op[i].multi_mode = pz->multi_mode;
        op[i].why[0] = '\0';
        op[i].status = eslOK;
      }

    if (want == 1) pz_open_one(&op[0]);
    else
      {
        if ((tid = malloc(sizeof(pthread_t) * want)) == NULL) { free(op); status = eslEMEM; goto ERROR; }
        for (i = 0; i < want; i++)
          if (pthread_create(&tid[i], NULL, pz_open_one, &op[i]) != 0) { free(op); free(tid); status = eslESYS; goto ERROR; }
        for (i = 0; i < want; i++) pthread_join(tid[i], NULL);
        free(tid);
      }

    /* release any half-open card, then compact the ones that came up */
    for (i = 0; i < want; i++)
      {
        if (op[i].status != eslOK)
          {
            pz_dev_close(&pz->dv[i]);
            if (pz_verbose()) fprintf(stderr, "# pezy: skipping a device (%s)\n", op[i].why);
            if (errbuf && ! errbuf[0]) { strncpy(errbuf, op[i].why, eslERRBUFSIZE-1); errbuf[eslERRBUFSIZE-1] = '\0'; }
            continue;
          }
        if (opened != i) pz->dv[opened] = pz->dv[i];
        pz->dv[opened].idx = opened;
        opened++;
      }
    free(op);
    if (opened == 0) { status = eslENOTFOUND; goto ERROR; }
    pz->ndev = opened;
  }

  if (pz_verbose())
    fprintf(stderr, "# pezy: %d device%s, batch %d thread%s/PE, multi %s\n",
            pz->ndev, pz->ndev == 1 ? "" : "s",
            pz->tpp_batch, pz->tpp_batch == 1 ? "" : "s",
            pz->multi_mode == P7_PZ_MULTI_SHARE  ? "share" :
            pz->multi_mode == P7_PZ_MULTI_PE     ? "pe"    : "thread");

  pz->M = -1;
  return pz;

 ERROR:
  if (pz) p7_pz_Destroy(pz);
  return NULL;
}

void
p7_pz_Destroy(P7_PZ *pz)
{
  int i;
  if (!pz) return;
  if (pz->dv) {
    for (i = 0; i < pz->ndev; i++) pz_dev_close(&pz->dv[i]);
    free(pz->dv);
  }
  if (pz->h_qoff) free(pz->h_qoff);
  if (pz->h_a)    free(pz->h_a);
  if (pz->h_b)    free(pz->h_b);
  if (pz->h_seq)  free(pz->h_seq);
  if (pz->h_part) free(pz->h_part);
  pthread_mutex_destroy(&pz->mx);
  free(pz);
}

/*****************************************************************
 * 3. Model upload
 *****************************************************************/

/* p7_pz_DestripeSSV()
 *
 * impl_sse stores the SSV match scores striped 16 ways:
 *   om->sbv[x][q] byte z  <->  model position k = q+1 + z*nq,  nq = p7O_NQB(M)
 * Undo that so the device can index sb[x*(M+1)+k] directly.
 */
int8_t *
p7_pz_DestripeSSV(const P7_OPROFILE *om)
{
  int      Kp  = om->abc->Kp;
  int      M   = om->M;
  int      nq  = p7O_NQB(M);
  int8_t  *sb  = NULL;
  int      x, k, q, z;
  int      status;

  ESL_ALLOC(sb, (size_t) Kp * (M+1));
  for (x = 0; x < Kp; x++)
    {
      sb[(size_t) x * (M+1)] = 0;                /* k=0 is never read */
      for (k = 1; k <= M; k++)
        {
          q = (k-1) % nq;
          z = (k-1) / nq;
          sb[(size_t) x * (M+1) + k] = ((const int8_t *)(om->sbv[x] + q))[z];
        }
    }
  return sb;

 ERROR:
  if (sb) free(sb);
  return NULL;
}

uint8_t *
p7_pz_BuildSSVTables(const P7_OPROFILE *om, int *ret_Mpad)
{
  int      Kp   = om->abc->Kp;
  int      M    = om->M;
  int      nq   = p7O_NQB(M);
  int      Mpad = (M + 7) & ~7;
  uint8_t *tab  = NULL;
  int      x, k, q, z;
  int      status;

  ESL_ALLOC(tab, (size_t) 2 * Kp * Mpad);

  /* Undo impl_sse's 16-way striping: sbv[x][q] byte z holds model position
   * k = q+1 + z*nq.  Walking (z,q) rather than k matters -- deriving q and z
   * from k costs an integer division per byte, which for a whole pressed
   * database (30,134 models, 4.75M positions x 29 residues) was most of the
   * table-building time. */
  for (x = 0; x < Kp; x++)
    {
      const int8_t *sv = (const int8_t *) om->sbv[x];
      uint8_t      *wp = tab + (size_t) x        * Mpad;
      uint8_t      *wn = tab + (size_t) (Kp + x) * Mpad;

      for (z = 0; z < 16; z++)
        {
          int kbase = z * nq + 1;              /* k at q = 0 */
          int qn    = M - kbase + 1;           /* how many q are still <= M   */
          const int8_t *svz;
          uint8_t      *wpz, *wnz;

          if (qn <= 0) break;
          if (qn > nq) qn = nq;
          svz = sv + z;                        /* stride 16 within one lane   */
          wpz = wp + kbase - 1;                /* contiguous in k             */
          wnz = wn + kbase - 1;
          for (q = 0; q < qn; q++)
            {
              int w = svz[q * 16];
              int m = w >> 31;                 /* -1 if negative, else 0      */
              wpz[q] = (uint8_t) ( w & ~m);
              wnz[q] = (uint8_t) ((-w) &  m);
            }
        }

      /* Padding beyond M: the striped table holds +127 there, which drives the
       * running value to the -infinity floor in one step, so padded lanes never
       * contribute to the maximum.  Reproduce that. */
      for (k = M + 1; k <= Mpad; k++) { wp[k-1] = 127; wn[k-1] = 0; }
    }

  *ret_Mpad = Mpad;
  return tab;

 ERROR:
  if (tab) free(tab);
  return NULL;
}

int
p7_pz_ModelIsSupported(const P7_OPROFILE *om)
{
  return (om->M >= 1 && om->M <= P7_PZ_MAXM) ? TRUE : FALSE;
}

int
p7_pz_SetModel(P7_PZ *pz, const P7_OPROFILE *om)
{
  uint8_t *tab = NULL;
  int      Mpad;
  size_t   n;
  int      i;
  int      status;

  if (! p7_pz_ModelIsSupported(om)) return eslENORESULT;

  pz->M  = om->M;
  pz->Kp = om->abc->Kp;

  if ((tab = p7_pz_BuildSSVTables(om, &Mpad)) == NULL) return eslEMEM;
  n = (size_t) 2 * pz->Kp * Mpad;

  /* the lmem kernel needs both score tables AND 8 DP rows in local memory */
  pz->use_lmem = (n + (size_t) 8 * Mpad <= (size_t) P7_PZ_LMEM_AVAIL) ? TRUE : FALSE;

  /* A few kilobytes to each card; not worth a thread apiece. */
  for (i = 0; i < pz->ndev; i++)
    {
      P7_PZDEV *dd = &pz->dv[i];
      if ((status = pz_grow(dd->ctx, &dd->d_sb, &dd->d_sb_sz, n)) != eslOK) { free(tab); return status; }
      /* Uploads are always explicit writes. */
      if (pzclEnqueueWriteBuffer(dd->q, dd->d_sb, PZCL_TRUE, 0, n, tab, 0, NULL, NULL) != PZCL_SUCCESS)
        { free(tab); return eslFAIL; }
    }
  free(tab);
  return eslOK;
}

/*****************************************************************
 * 4. Sequence batches
 *****************************************************************/

P7_PZ_SEQBATCH *
p7_pz_seqbatch_Create(int nalloc, int64_t resalloc)
{
  P7_PZ_SEQBATCH *b = NULL;
  int status;

  if (nalloc   < 16)   nalloc   = 16;
  if (resalloc < 4096) resalloc = 4096;

  ESL_ALLOC(b, sizeof(P7_PZ_SEQBATCH));
  b->dsq = NULL; b->off = NULL; b->len = NULL;
  b->n = 0; b->nalloc = nalloc; b->nres = 0; b->resalloc = resalloc; b->borrowed = FALSE;
  ESL_ALLOC(b->dsq, (size_t) resalloc);
  ESL_ALLOC(b->off, sizeof(int64_t) * nalloc);
  ESL_ALLOC(b->len, sizeof(int)     * nalloc);
  return b;

 ERROR:
  p7_pz_seqbatch_Destroy(b);
  return NULL;
}

void
p7_pz_seqbatch_Destroy(P7_PZ_SEQBATCH *b)
{
  if (!b) return;
  if (! b->borrowed) {
    if (b->dsq) free(b->dsq);
    if (b->off) free(b->off);
    if (b->len) free(b->len);
  }
  free(b);
}

void
p7_pz_seqbatch_Reuse(P7_PZ_SEQBATCH *b) { b->n = 0; b->nres = 0; }

int
p7_pz_seqbatch_Borrow(P7_PZ_SEQBATCH *b, uint8_t *dsq, int64_t *off, int *len, int n, int64_t nres)
{
  if (! b->borrowed) {
    if (b->dsq) free(b->dsq);
    if (b->off) free(b->off);
    if (b->len) free(b->len);
  }
  b->borrowed = TRUE;
  b->dsq = dsq;  b->off = off;  b->len = len;
  b->n   = n;    b->nres = nres;
  b->nalloc = n; b->resalloc = nres;
  return eslOK;
}

int
p7_pz_seqbatch_Add(P7_PZ_SEQBATCH *b, const ESL_DSQ *dsq, int L)
{
  void *tmp;
  int   status;

  if (b->n == b->nalloc) {
    int na = b->nalloc * 2;
    ESL_RALLOC(b->off, tmp, sizeof(int64_t) * na);
    ESL_RALLOC(b->len, tmp, sizeof(int)     * na);
    b->nalloc = na;
  }
  /* we store dsq[0..L]: the leading sentinel keeps device indexing 1-based */
  if (b->nres + (int64_t)L + 1 > b->resalloc) {
    int64_t ra = b->resalloc * 2;
    while (ra < b->nres + (int64_t)L + 1) ra *= 2;
    ESL_RALLOC(b->dsq, tmp, (size_t) ra);
    b->resalloc = ra;
  }

  b->off[b->n] = b->nres;
  b->len[b->n] = L;
  memcpy(b->dsq + b->nres, dsq, (size_t) L + 1);
  b->nres += (int64_t) L + 1;
  b->n++;
  return eslOK;

 ERROR:
  return status;
}

/*****************************************************************
 * 5. Running the filter
 *****************************************************************/

/* Chunk length, in residues.
 *
 * Long sequences must be split or they set the kernel's critical path: in
 * Swiss-Prot the longest sequence is 35213 residues, several times the ideal
 * per-thread load, so with one sequence per thread most of the device would
 * sit idle waiting for that one.  Each split costs an (M-1)-residue warm-up,
 * so the chunk is scaled with M to keep that overhead near 1/8.
 */
static int
pz_chunk_len(int M)
{
  int c = 8 * M;
  if (c < 1024) c = 1024;
  return c;
}

/* Build the work-item list for <b> at the current M.  Sequences shorter than
 * one chunk produce a single item, so the list is barely longer than the batch
 * for a typical protein database.
 */
static int
pz_build_items(P7_PZ *pz, const P7_PZ_SEQBATCH *b)
{
  int   chunk = pz_chunk_len(pz->M);
  int   need  = 0;
  int   i, a;
  void *tmp;
  int   status;

  for (i = 0; i < b->n; i++) need += (b->len[i] + chunk - 1) / chunk;
  if (need < 1) need = 1;

  if (need > pz->item_alloc) {
    int na = need + need / 4;
    ESL_RALLOC(pz->h_qoff, tmp, sizeof(int64_t) * na);
    ESL_RALLOC(pz->h_a,    tmp, sizeof(int32_t) * na);
    ESL_RALLOC(pz->h_b,    tmp, sizeof(int32_t) * na);
    ESL_RALLOC(pz->h_seq,  tmp, sizeof(int32_t) * na);
    ESL_RALLOC(pz->h_part, tmp, sizeof(uint8_t) * na);
    pz->item_alloc = na;
  }

  pz->nitem = 0;
  for (i = 0; i < b->n; i++)
    for (a = 1; a <= b->len[i]; a += chunk)
      {
        int j = pz->nitem++;
        pz->h_qoff[j] = b->off[i];
        pz->h_a[j]    = a;
        pz->h_b[j]    = (a + chunk - 1 > b->len[i]) ? b->len[i] : a + chunk - 1;
        pz->h_seq[j]  = i;
      }
  return eslOK;

 ERROR:
  return status;
}

/* One launch over the item range [ifirst,ilast).
 *
 * The kernel reads its work cursor from counter[0] and stops at the <nitem>
 * argument, so seeding the cursor with <ifirst> and passing <ilast> as <nitem>
 * restricts a launch to that range without any change to the kernel's
 * signature.  <nitem> is argument P7_PZ_NITEM_ARG in all three kernels; every
 * other argument is expected to be set already.
 *
 * Returns eslOK only if every thread reported completion.  This check is not
 * optional: the output of a launch that did not run to the end must not be
 * used, and this tally is how the host knows -- without it a whole hmmscan run
 * could come out silently wrong.
 */
/* HMMER_PEZY_FAILFIRST=<n> makes the first <n> launches of the process report
 * failure.  The watchdog path is otherwise unreachable on any workload we have
 * -- at the rates this device achieves, no single launch a 1 s budget produces
 * comes close to the driver's limit -- and an untested recovery path is
 * how a silent wrong answer gets shipped.  Diagnostic only; default off. */
static int
pz_failfirst(void)
{
  const char *v = getenv("HMMER_PEZY_FAILFIRST");
  if (v) { int n = atoi(v); if (n > 0) return n; }
  return 0;
}
static int pz_failfirst_left = -1;

static int
pz_run_range(P7_PZDEV *dd, pzcl_kernel k, pzcl_mem d_cnt, int ifirst, int ilast)
{
  int32_t cnt[P7_PZ_CNT_N];
  size_t  gws = dd->gws;

  if (pz_failfirst_left < 0) pz_failfirst_left = pz_failfirst();
  if (pz_failfirst_left > 0)
    {
      pz_failfirst_left--;
      if (pz_verbose())
        fprintf(stderr, "# pezy: dev %d launch [%d,%d) failed on purpose (HMMER_PEZY_FAILFIRST)\n",
                dd->idx, ifirst, ilast);
      return eslFAIL;
    }

  memset(cnt, 0, sizeof(cnt));
  cnt[0] = ifirst;
  if (pzclEnqueueWriteBuffer(dd->q, d_cnt, PZCL_TRUE, 0, sizeof(cnt), cnt, 0, NULL, NULL) != PZCL_SUCCESS) return eslFAIL;
  if (pzclSetKernelArg(k, P7_PZ_NITEM_ARG, sizeof(int), &ilast)                           != PZCL_SUCCESS) return eslFAIL;
  if (pzclEnqueueNDRangeKernel(dd->q, k, 1, NULL, &gws, NULL, 0, NULL, NULL)              != PZCL_SUCCESS) return eslFAIL;
  if (pzclFinish(dd->q)                                                                   != PZCL_SUCCESS) return eslFAIL;
  if (pzclEnqueueReadBuffer(dd->q, d_cnt, PZCL_TRUE, 0, sizeof(cnt), cnt, 0, NULL, NULL)  != PZCL_SUCCESS) return eslFAIL;

  if (cnt[P7_PZ_CNT_DONE] != (int32_t) gws || cnt[0] < ilast)
    {
      if (pz_verbose())
        fprintf(stderr, "# pezy: dev %d launch [%d,%d) did not finish (cursor %d, %d/%d threads done)\n",
                dd->idx, ifirst, ilast, (int) cnt[0], (int) cnt[P7_PZ_CNT_DONE], (int) gws);
      return eslFAIL;
    }
  return eslOK;
}

/* hmmscan's cards share one profile database rather than dividing it (see
 * P7_PZMMDEV), so they pull model ranges off a common cursor.  Handing out
 * 1/ncard of the work per grab is what keeps them from starving each other;
 * more grabs than that only add launch tails. */
#define P7_PZ_SPLIT  1

static int
pz_split(void)
{
  const char *v = getenv("HMMER_PEZY_SPLIT");
  int         n = (v ? atoi(v) : P7_PZ_SPLIT);
  return (n >= 1) ? n : P7_PZ_SPLIT;
}

/* Cut the database into one contiguous slice per card.
 *
 * The alternative -- give every card the whole database and let them pull from
 * a shared cursor -- balances better, but it puts the whole database on every
 * card: four times the bytes over PCIe, paid again in every process, and
 * hmmsearch runs one query per process.  Dividing it instead means the cards
 * upload a quarter each, in parallel, so adding cards costs nothing and a
 * search that ran one launch still runs one launch per card.  What it gives up
 * is the ability of an early-finishing card to help a late one; the cards are
 * identical and the slices are equal, so in practice they land within a few
 * percent of each other.
 *
 * The cut is by cumulative residues, on sequence boundaries.  Every item of
 * one call has the same Mpad and the bands of a sequence tile it exactly, so
 * equal residues *is* equal cells -- and unlike a cut on cells it does not
 * move when M does.  That matters because jackhmmer searches the same database
 * once per iteration with a different M each time: the item list is rebuilt,
 * but the hundreds of megabytes of sequence stay where they are.
 */
static int
pz_partition(P7_PZ *pz, const P7_PZ_SEQBATCH *b)
{
  int64_t nres = 0, acc = 0;
  int     s, c;

  /* Cut on residues, not on <b->nres>.  P7_PZDB lays a record's digital form
   * down at its own file offset, so the array is sparse -- the gaps are the
   * headers and newlines -- and b->nres is the span of that array, half again
   * the residues it holds.  Cutting against the span leaves the last card
   * with nothing. */
  for (s = 0; s < b->n; s++) nres += b->len[s];

  for (c = 0, s = 0; c < pz->ndev; c++)
    {
      P7_PZDEV *dd   = &pz->dv[c];
      int64_t   want = nres * (int64_t) (c + 1) / pz->ndev;

      dd->seq0 = s;
      if (c == pz->ndev - 1) s = b->n;          /* last card takes the rest */
      else while (s < b->n && acc < want) { acc += b->len[s]; s++; }
      dd->seq1 = s;

      if (dd->seq1 <= dd->seq0) { dd->dsq_base = 0; dd->dsq_len = 0; continue; }
      dd->dsq_base = b->off[dd->seq0];
      dd->dsq_len  = b->off[dd->seq1 - 1] + (int64_t) b->len[dd->seq1 - 1] + 2 - dd->dsq_base;
    }

  pz->part_ndev = pz->ndev;
  pz->part_n    = b->n;
  pz->part_nres = b->nres;
  return eslOK;
}

/* Map each card's sequence range onto the current item list, and rebase its
 * item offsets onto its slice of <dsq>.  Redone whenever the item list is
 * rebuilt (that is, whenever M changes the chunk length); the sequence
 * partition itself, and so the uploaded sequence bytes, do not move.
 */
static int
pz_partition_items(P7_PZ *pz)
{
  int   nit = pz->nitem;
  int   i = 0, c, j;
  void *tmp;
  int   status;

  for (c = 0; c < pz->ndev; c++)
    {
      P7_PZDEV *dd = &pz->dv[c];

      dd->it0 = i;
      while (i < nit && pz->h_seq[i] < dd->seq1) i++;
      dd->it1 = i;

      if (dd->it1 - dd->it0 > dd->qoff_alloc)
        {
          int na = (dd->it1 - dd->it0) + (dd->it1 - dd->it0) / 4 + 1;
          ESL_RALLOC(dd->h_qoff, tmp, sizeof(int64_t) * (size_t) na);
          dd->qoff_alloc = na;
        }
      for (j = dd->it0; j < dd->it1; j++) dd->h_qoff[j - dd->it0] = pz->h_qoff[j] - dd->dsq_base;
    }

  pz->part_nitem = nit;
  return eslOK;

 ERROR:
  return status;
}

/* One card's share of a p7_pz_SSVBatch() call. */
typedef struct {
  P7_PZ                *pz;
  P7_PZDEV             *dd;
  const P7_PZ_SEQBATCH *b;
  int64_t               cells;    /* cells this card actually ran       */
  double                tup;      /* seconds spent uploading to it      */
  int                   nlaunch;
  int                   status;
} pz_batch_work_t;

/* Get this card's slice of the database and of the item list onto it, then run
 * the slice in watchdog-sized launches.
 *
 * The upload is inside the worker on purpose: the cards' slices are disjoint,
 * so four threads move four quarters at once instead of one thread moving the
 * whole thing.
 */
static void *
pz_batch_worker(void *arg)
{
  pz_batch_work_t *wk  = (pz_batch_work_t *) arg;
  P7_PZ           *pz  = wk->pz;
  P7_PZDEV        *dd  = wk->dd;
  const P7_PZ_SEQBATCH *b = wk->b;
  pzcl_kernel      k;
  pzcl_int         err;
  int              nc    = dd->it1 - dd->it0;      /* items on this card */
  int              zero  = 0;
  int              chunk = pz_chunk_len(pz->M);
  int              same_dsq, same_it;
  double           tup;
  int              status;

  wk->status = eslOK;
  if (nc <= 0) return NULL;
  tup = pz_now();

  /* The sequence slice and the item slice go stale independently: a new query
   * model changes the chunk length and so the item list, but leaves the
   * hundreds of megabytes of sequence exactly where they are. */
  same_dsq = (dd->up_n == b->n && dd->up_nres == b->nres &&
              dd->up_base == dd->dsq_base && dd->up_len == dd->dsq_len);
  same_it  = (same_dsq && dd->up_chunk == chunk &&
              dd->up_it0 == dd->it0 && dd->up_it1 == dd->it1);

  if ((status = pz_grow(dd->ctx, &dd->d_dsq,  &dd->d_dsq_sz,  (size_t) dd->dsq_len))          != eslOK) goto ERROR;
  if ((status = pz_grow(dd->ctx, &dd->d_qoff, &dd->d_qoff_sz, sizeof(int64_t) * (size_t) nc)) != eslOK) goto ERROR;
  if ((status = pz_grow(dd->ctx, &dd->d_a,    &dd->d_a_sz,    sizeof(int32_t) * (size_t) nc)) != eslOK) goto ERROR;
  if ((status = pz_grow(dd->ctx, &dd->d_b,    &dd->d_b_sz,    sizeof(int32_t) * (size_t) nc)) != eslOK) goto ERROR;
  if ((status = pz_grow(dd->ctx, &dd->d_xE,   &dd->d_xE_sz,   (size_t) nc))                   != eslOK) goto ERROR;
  if ((status = pz_grow(dd->ctx, &dd->d_cnt,  &dd->d_cnt_sz,  sizeof(int32_t) * P7_PZ_CNT_N)) != eslOK) goto ERROR;

  /* Uploads are always explicit writes. */
  if (! same_dsq &&
      pzclEnqueueWriteBuffer(dd->q, dd->d_dsq, PZCL_FALSE, 0, (size_t) dd->dsq_len, b->dsq + dd->dsq_base, 0, NULL, NULL) != PZCL_SUCCESS)
    { status = eslFAIL; goto ERROR; }
  if (! same_it &&
      (pzclEnqueueWriteBuffer(dd->q, dd->d_qoff, PZCL_FALSE, 0, sizeof(int64_t)*(size_t)nc, dd->h_qoff,        0, NULL, NULL) != PZCL_SUCCESS ||
       pzclEnqueueWriteBuffer(dd->q, dd->d_a,    PZCL_FALSE, 0, sizeof(int32_t)*(size_t)nc, pz->h_a + dd->it0, 0, NULL, NULL) != PZCL_SUCCESS ||
       pzclEnqueueWriteBuffer(dd->q, dd->d_b,    PZCL_FALSE, 0, sizeof(int32_t)*(size_t)nc, pz->h_b + dd->it0, 0, NULL, NULL) != PZCL_SUCCESS))
    { status = eslFAIL; goto ERROR; }

  /* The work counter is zeroed from the host rather than in the kernel, which
   * avoids needing a device-wide barrier before the first pz_atomic_add(). */
  if (pzclEnqueueWriteBuffer(dd->q, dd->d_cnt, PZCL_TRUE, 0, sizeof(int32_t), &zero, 0, NULL, NULL) != PZCL_SUCCESS)
    { status = eslFAIL; goto ERROR; }

  dd->up_n = b->n; dd->up_nres = b->nres;
  dd->up_base = dd->dsq_base; dd->up_len = dd->dsq_len;
  dd->up_chunk = chunk; dd->up_it0 = dd->it0; dd->up_it1 = dd->it1;
  wk->tup = pz_now() - tup;

  k = pz->use_lmem ? dd->k_lmem : dd->k_gmem;
  {
    int arg = 0;
    err  = pzclSetKernelArg(k, arg++, sizeof(pzcl_mem), &dd->d_sb);
    err |= pzclSetKernelArg(k, arg++, sizeof(int),      &pz->Kp);
    err |= pzclSetKernelArg(k, arg++, sizeof(int),      &pz->M);
    err |= pzclSetKernelArg(k, arg++, sizeof(pzcl_mem), &dd->d_dsq);
    err |= pzclSetKernelArg(k, arg++, sizeof(pzcl_mem), &dd->d_qoff);
    err |= pzclSetKernelArg(k, arg++, sizeof(pzcl_mem), &dd->d_a);
    err |= pzclSetKernelArg(k, arg++, sizeof(pzcl_mem), &dd->d_b);
    err |= pzclSetKernelArg(k, arg++, sizeof(int),      &nc);   /* per-launch, reset by pz_run_range() */
    err |= pzclSetKernelArg(k, arg++, sizeof(pzcl_mem), &dd->d_xE);
    err |= pzclSetKernelArg(k, arg++, sizeof(pzcl_mem), &dd->d_cnt);
    if (err != PZCL_SUCCESS) { status = eslFAIL; goto ERROR; }
  }

  /* Run the slice in pieces short enough for the driver's watchdog.  One item
   * costs Mpad * its band length, so a piece is as many consecutive items as
   * fit the cell budget -- always at least one, however long that one is.
   * Every extra piece costs a tail (the launch ends when its slowest thread
   * finishes its last grab), so the budget is kept as large as the watchdog
   * allows. */
  {
    int Mpad = (pz->M + 7) & ~7;
    int i0   = 0;

    while (i0 < nc)
      {
        int64_t budget = pz_cell_budget(dd->cps);
        int64_t cells  = 0;
        int     i1     = i0;
        double  t0;

        while (i1 < nc)
          {
            int64_t c = (int64_t) Mpad * (int64_t) (pz->h_b[dd->it0 + i1] - pz->h_a[dd->it0 + i1] + 1);
            if (i1 > i0 && cells + c > budget) break;
            cells += c; i1++;
          }

        t0 = pz_now();
        if ((status = pz_run_range(dd, k, dd->d_cnt, i0, i1)) != eslOK)
          {
            /* Almost certainly the watchdog: the tally came up short.  Halve
             * the rate estimate -- which halves the budget -- and run this
             * range again, smaller.  The retry rewrites the same xE slice, so
             * a killed launch's partial output is harmless.  Give up only when
             * a single item will not run. */
            if (i1 - i0 <= 1) goto ERROR;
            dd->cps *= 0.5;
            if (dd->cps < P7_PZ_CPS_MIN) dd->cps = P7_PZ_CPS_MIN;
            if (pz_verbose())
              fprintf(stderr, "# pezy: dev %d launch [%d,%d) came up short; budget halved to %.3g cells\n",
                      dd->idx, i0, i1, (double) pz_cell_budget(dd->cps));
            continue;
          }
        pz_update_cps(&dd->cps, cells, pz_now() - t0);

        /* the slices are disjoint, so no lock is needed on h_part */
        if (pzclEnqueueReadBuffer(dd->q, dd->d_xE, PZCL_TRUE, (size_t) i0, (size_t) (i1 - i0),
                                  pz->h_part + dd->it0 + i0, 0, NULL, NULL) != PZCL_SUCCESS)
          { status = eslFAIL; goto ERROR; }

        wk->cells += cells; wk->nlaunch++;
        i0 = i1;
      }
  }
  return NULL;

 ERROR:
  wk->status = status;
  return NULL;
}

int
p7_pz_SSVBatch(P7_PZ *pz, const P7_PZ_SEQBATCH *b, uint8_t *xE)
{
  pz_batch_work_t *wk  = NULL;
  pthread_t       *tid = NULL;
  int64_t          total = 0;
  int              n   = b->n;
  int              nit;
  int              Mpad;
  int              i;
  int              status;

  if (n == 0)     return eslOK;
  if (pz->M <= 0) return eslEINVAL;

  /* The item list is host-side and shared by every card.  It only has to be
   * rebuilt when the chunk length changes, which follows M. */
  {
    int chunk = pz_chunk_len(pz->M);
    if (pz->nitem == 0 || pz->item_chunk != chunk || pz->item_n != n || pz->item_nres != b->nres)
      {
        if ((status = pz_build_items(pz, b)) != eslOK) return status;
        pz->item_chunk = chunk; pz->item_n = n; pz->item_nres = b->nres;
        pz->part_nitem = -1;                 /* the item slices moved with it */
      }
  }
  nit  = pz->nitem;
  Mpad = (pz->M + 7) & ~7;
  for (i = 0; i < nit; i++) total += (int64_t) Mpad * (int64_t) (pz->h_b[i] - pz->h_a[i] + 1);

  if (pz->part_ndev != pz->ndev || pz->part_n != n || pz->part_nres != b->nres)
    { if ((status = pz_partition(pz, b))       != eslOK) return status; pz->part_nitem = -1; }
  if (pz->part_nitem != nit)
    { if ((status = pz_partition_items(pz))    != eslOK) return status; }

  ESL_ALLOC(wk, sizeof(pz_batch_work_t) * pz->ndev);
  for (i = 0; i < pz->ndev; i++)
    {
      wk[i].pz = pz; wk[i].dd = &pz->dv[i]; wk[i].b = b;
      wk[i].cells = 0; wk[i].tup = 0.0; wk[i].nlaunch = 0; wk[i].status = eslOK;
    }

  {
    double tall = pz_now();

    if (pz->ndev == 1) pz_batch_worker(&wk[0]);
    else
      {
        ESL_ALLOC(tid, sizeof(pthread_t) * pz->ndev);
        for (i = 0; i < pz->ndev; i++)
          if (pthread_create(&tid[i], NULL, pz_batch_worker, &wk[i]) != 0) { status = eslESYS; goto ERROR; }
        for (i = 0; i < pz->ndev; i++) pthread_join(tid[i], NULL);
        free(tid); tid = NULL;
      }

    for (i = 0; i < pz->ndev; i++) if (wk[i].status != eslOK) { status = wk[i].status; goto ERROR; }

    if (pz_verbose())
      {
        int    nlaunch = 0;
        double tup = 0.0;
        tall = pz_now() - tall;
        for (i = 0; i < pz->ndev; i++) { nlaunch += wk[i].nlaunch; if (wk[i].tup > tup) tup = wk[i].tup; }
        fprintf(stderr, "# pezy: ssv_%s M=%d, %d seqs, %d items, %d dev, %d launches, "
                        "%.1f Gcells in %.3fs (upload %.3fs) = %.0f Gcells/s",
                pz->use_lmem ? "lmem" : "gmem", pz->M, n, nit, pz->ndev, nlaunch,
                total / 1e9, tall, tup, total / tall / 1e9);
        if (pz->ndev > 1)
          for (i = 0; i < pz->ndev; i++) fprintf(stderr, " [d%d %.0f%%]", i, 100.0 * wk[i].cells / (double) total);
        fprintf(stderr, "\n");
      }
  }
  free(wk); wk = NULL;

  /* fold the per-item maxima back into one xE per sequence */
  for (i = 0; i < n;   i++) xE[i] = 0x80;
  for (i = 0; i < nit; i++) {
    int s = pz->h_seq[i];
    if (pz->h_part[i] > xE[s]) xE[s] = pz->h_part[i];
  }
  return eslOK;

 ERROR:
  if (tid) free(tid);
  if (wk)  free(wk);
  return status;
}

/*****************************************************************
 * 6. Host-side tail of p7_SSVFilter()
 *****************************************************************/

/* Verbatim transcription of the tail of p7_SSVFilter() (impl_sse/ssvfilter.c),
 * reading a precomputed xE instead of calling get_xE().  Keeping the float
 * arithmetic here -- and identical -- is what makes the offload exact.
 */
int
p7_pz_SSVOutcome(uint8_t xE8, const P7_OPROFILE *om, float *ret_sc)
{
  uint16_t xE = (uint16_t) xE8;
  uint16_t xJ;

  if (om->tjb_b + om->tbm_b + om->tec_b + om->bias_b >= 127) return eslENORESULT;

  if (xE >= 255 - om->bias_b)
    {
      *ret_sc = eslINFINITY;
      if (om->base_b - om->tjb_b - om->tbm_b < 128) return eslENORESULT;
      return eslERANGE;
    }

  xE += om->base_b - om->tjb_b - om->tbm_b;
  xE -= 128;

  if (xE >= 255 - om->bias_b) { *ret_sc = eslINFINITY; return eslERANGE; }

  xJ = xE - om->tec_b;
  if (xJ > om->base_b) return eslENORESULT;

  *ret_sc  = ((float) (xJ - om->tjb_b) - (float) om->base_b);
  *ret_sc /= om->scale_b;
  *ret_sc -= 3.0;
  return eslOK;
}

#else /*!HAVE_PEZY*/
/* keep the translation unit non-empty for pedantic compilers */
int p7_pz_dummy_symbol_to_avoid_empty_translation_unit = 0;
#endif /*HAVE_PEZY*/

/*****************************************************************
 * 7. Many models x many queries (hmmscan's shape)
 *****************************************************************/

/* One card's copy of the profile database and its per-call scratch.
 *
 * The score tables are replicated rather than partitioned.  Pfam-A comes to
 * about 1.5 GB, a small fraction of a card's memory, so a copy per card costs
 * nothing that matters, and in exchange any card can run any model -- which is what lets
 * the cards pull from one shared model cursor and balance themselves.  A
 * partitioned layout would pin each card to a fixed slice of Pfam.
 */
typedef struct {
  pzcl_mem  d_tab;  size_t d_tab_sz;
  pzcl_mem  d_mdl;  size_t d_mdl_sz;
  pzcl_mem  d_pq;   size_t d_pq_sz;
  pzcl_mem  d_dsq;  size_t d_dsq_sz;
  pzcl_mem  d_xE;   size_t d_xE_sz;
  pzcl_mem  d_cnt;  size_t d_cnt_sz;
  double    cps;    /* running estimate of this card's cells/s */
} P7_PZMMDEV;

struct p7_pzmm_s {
  int       nmodel;
  int       Kp;
  int       rowstride;    /* bytes of local memory per thread's DP row      */
  int       ndev;         /* models actually on the device                  */

  int64_t  *toff;         /* [nmodel] table offset, or -1 if host-only      */
  int      *devof;        /* [nmodel] index in the device list, or -1       */
  int      *devmodel;     /* [ndev]  back to the caller's model index       */

  int64_t  *h_mdl;        /* [2*ndev] (table offset, M) descriptors         */
  int64_t   tabsz;

  int          ncard;     /* cards this P7_PZMM is spread over              */
  P7_PZMMDEV  *cv;        /* [ncard]                                        */

  /* per-call staging, grown and reused */
  int64_t  *h_pq;     int  pq_alloc, npq;
  int      *pq_owner;                        /* pseudo-query -> query index  */
  uint8_t  *h_part;   int  part_alloc;
  P7_PZ_SEQBATCH *qb;                        /* packed query sequences       */
};

void
p7_pzmm_Destroy(P7_PZMM *mm)
{
  int c;
  if (!mm) return;
  if (mm->cv) {
    for (c = 0; c < mm->ncard; c++) {
      P7_PZMMDEV *cc = &mm->cv[c];
      if (cc->d_tab) pzclReleaseMemObject(cc->d_tab);
      if (cc->d_mdl) pzclReleaseMemObject(cc->d_mdl);
      if (cc->d_pq)  pzclReleaseMemObject(cc->d_pq);
      if (cc->d_dsq) pzclReleaseMemObject(cc->d_dsq);
      if (cc->d_xE)  pzclReleaseMemObject(cc->d_xE);
      if (cc->d_cnt) pzclReleaseMemObject(cc->d_cnt);
    }
    free(mm->cv);
  }
  if (mm->toff)     free(mm->toff);
  if (mm->devof)    free(mm->devof);
  if (mm->devmodel) free(mm->devmodel);
  if (mm->h_mdl)    free(mm->h_mdl);
  if (mm->h_pq)     free(mm->h_pq);
  if (mm->pq_owner) free(mm->pq_owner);
  if (mm->h_part)   free(mm->h_part);
  if (mm->qb)       p7_pz_seqbatch_Destroy(mm->qb);
  free(mm);
}

int p7_pzmm_ModelOnDevice(const P7_PZMM *mm, int m) { return (mm->toff[m] >= 0) ? TRUE : FALSE; }
int p7_pzmm_NDevice      (const P7_PZMM *mm)        { return mm->ndev; }

/* one thread's share of the table build */
typedef struct {
  P7_PZMM      *mm;
  P7_OPROFILE **list;
  uint8_t      *tab;
  int           lo, hi;      /* device-list indices [lo,hi) */
  int           bad;
} pzmm_work_t;

static void *
pzmm_build_range(void *arg)
{
  pzmm_work_t *wk = (pzmm_work_t *) arg;
  P7_PZMM     *mm = wk->mm;
  int          d;

  for (d = wk->lo; d < wk->hi; d++)
    {
      int      m    = mm->devmodel[d];
      int      Mpad;
      uint8_t *t    = p7_pz_BuildSSVTables(wk->list[m], &Mpad);
      if (!t) { wk->bad = 1; return NULL; }
      memcpy(wk->tab + mm->toff[m], t, (size_t) 2 * mm->Kp * Mpad);
      free(t);
    }
  return NULL;
}

/* one card's share of the table upload */
typedef struct {
  P7_PZMM     *mm;
  P7_PZDEV    *dd;
  P7_PZMMDEV  *cc;
  const uint8_t *tab;
  int          status;
} pzmm_up_t;

static void *
pzmm_upload_card(void *arg)
{
  pzmm_up_t  *up = (pzmm_up_t *) arg;
  P7_PZMM    *mm = up->mm;
  P7_PZDEV   *dd = up->dd;
  P7_PZMMDEV *cc = up->cc;
  int         status;

  if ((status = pz_grow(dd->ctx, &cc->d_tab, &cc->d_tab_sz, (size_t) mm->tabsz))                     != eslOK) goto ERROR;
  if ((status = pz_grow(dd->ctx, &cc->d_mdl, &cc->d_mdl_sz, sizeof(int64_t) * 2 * (size_t) mm->ndev)) != eslOK) goto ERROR;
  if (pzclEnqueueWriteBuffer(dd->q, cc->d_tab, PZCL_FALSE, 0, (size_t) mm->tabsz,                 up->tab,   0, NULL, NULL) != PZCL_SUCCESS ||
      pzclEnqueueWriteBuffer(dd->q, cc->d_mdl, PZCL_TRUE,  0, sizeof(int64_t)*2*(size_t)mm->ndev, mm->h_mdl, 0, NULL, NULL) != PZCL_SUCCESS)
    { status = eslFAIL; goto ERROR; }
  up->status = eslOK;
  return NULL;

 ERROR:
  up->status = status;
  return NULL;
}

int
p7_pzmm_Create(P7_PZ *pz, P7_OPROFILE **list, int nmodel, int nthreads,
               P7_PZMM **ret_mm, char *errbuf)
{
  P7_PZMM     *mm  = NULL;
  uint8_t     *tab = NULL;
  pzmm_work_t *wk  = NULL;
  pthread_t   *tid = NULL;
  int          cap;                 /* largest padded M a DP row can hold */
  int          m, d, t;
  int          status;

  *ret_mm = NULL;
  if (nmodel <= 0) { if (errbuf) sprintf(errbuf, "no models"); return eslEINVAL; }
  if (nthreads < 1) nthreads = 1;
  if (nthreads > 32) nthreads = 32;

  ESL_ALLOC(mm, sizeof(P7_PZMM));
  memset(mm, 0, sizeof(P7_PZMM));
  mm->nmodel = nmodel;
  mm->Kp     = list[0]->abc->Kp;
  /* 8 threads share the local memory, and ssv_multi_share keeps the last 64
   * bytes of it for the PE's work cursor.  In Pfam-A this rules out one more
   * model (Mpad 2048) than the uncapped 2048 would. */
  cap        = (P7_PZ_LMEM_AVAIL - 64) / 8;

  ESL_ALLOC(mm->toff,     sizeof(int64_t) * nmodel);
  ESL_ALLOC(mm->devof,    sizeof(int)     * nmodel);
  ESL_ALLOC(mm->devmodel, sizeof(int)     * nmodel);

  /* Pass 1: decide who goes on the device, lay out the tables, and set the row
   * stride from the largest model we accept. */
  mm->tabsz = 0; mm->ndev = 0; mm->rowstride = 8;
  for (m = 0; m < nmodel; m++)
    {
      int Mpad = (list[m]->M + 7) & ~7;
      mm->toff[m] = -1; mm->devof[m] = -1;
      if (list[m]->M < 1 || Mpad > cap) continue;      /* host keeps this one */
      mm->toff[m]              = mm->tabsz;
      mm->devof[m]             = mm->ndev;
      mm->devmodel[mm->ndev++] = m;
      mm->tabsz               += (int64_t) 2 * mm->Kp * Mpad;
      if (Mpad > mm->rowstride) mm->rowstride = Mpad;
    }
  if (mm->ndev == 0) { if (errbuf) sprintf(errbuf, "no model fits the device"); status = eslENORESULT; goto ERROR; }

  /* Pass 2: build the tables, in parallel over disjoint ranges. */
  ESL_ALLOC(tab, (size_t) mm->tabsz);
  ESL_ALLOC(wk, sizeof(pzmm_work_t) * nthreads);
  for (t = 0; t < nthreads; t++)
    {
      wk[t].mm = mm; wk[t].list = list; wk[t].tab = tab; wk[t].bad = 0;
      wk[t].lo = (int) ((int64_t) mm->ndev *  t      / nthreads);
      wk[t].hi = (int) ((int64_t) mm->ndev * (t + 1) / nthreads);
    }
  if (nthreads == 1) pzmm_build_range(&wk[0]);
  else
    {
      ESL_ALLOC(tid, sizeof(pthread_t) * nthreads);
      for (t = 0; t < nthreads; t++)
        if (pthread_create(&tid[t], NULL, pzmm_build_range, &wk[t]) != 0) { status = eslESYS; goto ERROR; }
      for (t = 0; t < nthreads; t++) pthread_join(tid[t], NULL);
      free(tid); tid = NULL;
    }
  for (t = 0; t < nthreads; t++)
    if (wk[t].bad) { if (errbuf) sprintf(errbuf, "table build failed"); status = eslEMEM; goto ERROR; }
  free(wk); wk = NULL;

  /* device descriptors: (table offset, M) per device model */
  ESL_ALLOC(mm->h_mdl, sizeof(int64_t) * 2 * mm->ndev);
  for (d = 0; d < mm->ndev; d++)
    {
      m = mm->devmodel[d];
      mm->h_mdl[2*d]     = mm->toff[m];
      mm->h_mdl[2*d + 1] = list[m]->M;
    }

  /* Give every card its own copy.  These are gigabyte-scale transfers, so
   * they go from one thread per card and overlap instead of adding up. */
  mm->ncard = pz->ndev;
  ESL_ALLOC(mm->cv, sizeof(P7_PZMMDEV) * mm->ncard);
  memset(mm->cv, 0, sizeof(P7_PZMMDEV) * mm->ncard);
  for (t = 0; t < mm->ncard; t++) mm->cv[t].cps = P7_PZ_CPS0;

  {
    pzmm_up_t *up = NULL;
    int        bad = 0;

    ESL_ALLOC(up, sizeof(pzmm_up_t) * mm->ncard);
    for (t = 0; t < mm->ncard; t++)
      { up[t].mm = mm; up[t].dd = &pz->dv[t]; up[t].cc = &mm->cv[t]; up[t].tab = tab; up[t].status = eslOK; }

    if (mm->ncard == 1) pzmm_upload_card(&up[0]);
    else
      {
        ESL_ALLOC(tid, sizeof(pthread_t) * mm->ncard);
        for (t = 0; t < mm->ncard; t++)
          if (pthread_create(&tid[t], NULL, pzmm_upload_card, &up[t]) != 0) { free(up); status = eslESYS; goto ERROR; }
        for (t = 0; t < mm->ncard; t++) pthread_join(tid[t], NULL);
        free(tid); tid = NULL;
      }
    for (t = 0; t < mm->ncard; t++) if (up[t].status != eslOK) { bad = up[t].status; break; }
    free(up);
    if (bad) { if (errbuf) sprintf(errbuf, "table upload failed"); status = bad; goto ERROR; }
  }
  free(tab); tab = NULL;

  if ((mm->qb = p7_pz_seqbatch_Create(64, 65536)) == NULL) { status = eslEMEM; goto ERROR; }

  *ret_mm = mm;
  return eslOK;

 ERROR:
  if (tid) free(tid);
  if (wk)  free(wk);
  if (tab) free(tab);
  p7_pzmm_Destroy(mm);
  return status;
}

/* Build the pseudo-query list: one entry per query, except that a query longer
 * than <chunk> is cut into overlapping pieces so no single work item is
 * unbounded.  Each piece is preceded by <rowstride-1> residues of warm-up,
 * which is at least M-1 for every model on the device, so every row is computed
 * exactly by the piece that owns it.  Pieces overlap, but the reduction is a
 * max and a short chain can only underestimate, so duplicates are harmless.
 */
static int
pzmm_build_pq(P7_PZMM *mm, int nq, int chunk)
{
  int   need = 0, i, a;
  void *tmp;
  int   status;

  for (i = 0; i < nq; i++)
    { int L = mm->qb->len[i]; need += (L <= chunk) ? 1 : (L + chunk - 1) / chunk; }
  if (need < 1) need = 1;

  if (need > mm->pq_alloc)
    {
      int na = need + need / 4;
      ESL_RALLOC(mm->h_pq,     tmp, sizeof(int64_t) * 2 * na);
      ESL_RALLOC(mm->pq_owner, tmp, sizeof(int)         * na);
      mm->pq_alloc = na;
    }

  mm->npq = 0;
  for (i = 0; i < nq; i++)
    {
      int L = mm->qb->len[i];
      if (L <= chunk)
        {
          int j = mm->npq++;
          mm->h_pq[2*j]     = mm->qb->off[i];      /* q[1..L] = dsq[1..L] */
          mm->h_pq[2*j + 1] = L;
          mm->pq_owner[j]   = i;
        }
      else
        for (a = 1; a <= L; a += chunk)
          {
            int b     = (a + chunk - 1 > L) ? L : a + chunk - 1;
            int start = a - (mm->rowstride - 1); if (start < 1) start = 1;
            int j     = mm->npq++;
            mm->h_pq[2*j]     = mm->qb->off[i] + start - 1;   /* q[1] = dsq[start] */
            mm->h_pq[2*j + 1] = b - start + 1;
            mm->pq_owner[j]   = i;
          }
    }
  return eslOK;

 ERROR:
  return status;
}

/* One card's share of a p7_pzmm_SSV() call. */
typedef struct {
  P7_PZ      *pz;
  P7_PZDEV   *dd;
  P7_PZMM    *mm;
  P7_PZMMDEV *cc;
  int         nitem;      /* ndev * npq                              */
  int64_t     pqres;      /* residues over all pseudo-queries        */
  int64_t     total;      /* cells in the whole (model x query) grid */
  int64_t     cells;
  int         nlaunch;
  int         share;      /* TRUE if this card ran ssv_multi_share   */
  int         status;
} pzmm_work_run_t;

/* Upload this call's queries to one card, then grab model ranges off the
 * shared cursor until the profile database is exhausted.
 *
 * Item j is (model j/npq, query j%npq), so a contiguous item range is a
 * contiguous range of models -- which is what makes a shared cursor over
 * models the same thing as a shared cursor over items.
 */
static void *
pzmm_run_card(void *arg)
{
  pzmm_work_run_t *wk = (pzmm_work_run_t *) arg;
  P7_PZ           *pz = wk->pz;
  P7_PZDEV        *dd = wk->dd;
  P7_PZMM         *mm = wk->mm;
  P7_PZMMDEV      *cc = wk->cc;
  pzcl_kernel      km;
  pzcl_int         err;
  int              zero = 0;
  int              status;

  wk->status = eslOK;

  if ((status = pz_grow(dd->ctx, &cc->d_dsq, &cc->d_dsq_sz, (size_t) mm->qb->nres))                  != eslOK) goto ERROR;
  if ((status = pz_grow(dd->ctx, &cc->d_pq,  &cc->d_pq_sz,  sizeof(int64_t) * 2 * (size_t) mm->npq)) != eslOK) goto ERROR;
  if ((status = pz_grow(dd->ctx, &cc->d_xE,  &cc->d_xE_sz,  (size_t) wk->nitem))                     != eslOK) goto ERROR;
  if ((status = pz_grow(dd->ctx, &cc->d_cnt, &cc->d_cnt_sz, sizeof(int32_t) * P7_PZ_CNT_N))          != eslOK) goto ERROR;

  if (pzclEnqueueWriteBuffer(dd->q, cc->d_dsq, PZCL_FALSE, 0, (size_t) mm->qb->nres,             mm->qb->dsq, 0, NULL, NULL) != PZCL_SUCCESS ||
      pzclEnqueueWriteBuffer(dd->q, cc->d_pq,  PZCL_FALSE, 0, sizeof(int64_t)*2*(size_t)mm->npq, mm->h_pq,    0, NULL, NULL) != PZCL_SUCCESS ||
      pzclEnqueueWriteBuffer(dd->q, cc->d_cnt, PZCL_TRUE,  0, sizeof(int32_t),                   &zero,       0, NULL, NULL) != PZCL_SUCCESS)
    { status = eslFAIL; goto ERROR; }

  /* ssv_multi_share hands each PE eight consecutive items so its threads share
   * one model's table; with fewer than eight pseudo-queries per model there is
   * nothing to share, so the one-thread-per-PE kernel takes those. */
  km = (mm->npq >= 8 && dd->k_mshare) ? dd->k_mshare : dd->k_multi;
  wk->share = (km == dd->k_mshare);
  {
    int arg = 0;
    err  = pzclSetKernelArg(km, arg++, sizeof(pzcl_mem), &cc->d_mdl);
    err |= pzclSetKernelArg(km, arg++, sizeof(pzcl_mem), &cc->d_pq);
    err |= pzclSetKernelArg(km, arg++, sizeof(int),      &mm->npq);
    err |= pzclSetKernelArg(km, arg++, sizeof(int),      &mm->Kp);
    err |= pzclSetKernelArg(km, arg++, sizeof(int),      &mm->rowstride);
    err |= pzclSetKernelArg(km, arg++, sizeof(pzcl_mem), &cc->d_tab);
    err |= pzclSetKernelArg(km, arg++, sizeof(pzcl_mem), &cc->d_dsq);
    err |= pzclSetKernelArg(km, arg++, sizeof(int),      &wk->nitem); /* per-launch, reset by pz_run_range() */
    err |= pzclSetKernelArg(km, arg++, sizeof(pzcl_mem), &cc->d_xE);
    err |= pzclSetKernelArg(km, arg++, sizeof(pzcl_mem), &cc->d_cnt);
    if (err != PZCL_SUCCESS) { status = eslFAIL; goto ERROR; }
  }

  for (;;)
    {
      int64_t budget, cells = 0;
      int     d0, d1;
      double  t0;

      pthread_mutex_lock(&pz->mx);
      d0 = pz->cursor;
      if (d0 >= mm->ndev || pz->failed != eslOK) { pthread_mutex_unlock(&pz->mx); break; }
      {
        int64_t wb = pz_cell_budget(cc->cps);          /* watchdog ceiling  */
        int64_t sb = (mm->ncard > 1 && wk->total > 0)  /* sharing ceiling   */
                     ? wk->total / ((int64_t) mm->ncard * pz_split()) : 0;
        budget = (sb > 0 && sb < wb) ? sb : wb;
      }
      d1 = d0;
      while (d1 < mm->ndev)
        {
          int     Mp = ((int) mm->h_mdl[2*d1 + 1] + 7) & ~7;
          int64_t c  = (int64_t) Mp * wk->pqres;
          if (d1 > d0 && cells + c > budget) break;
          cells += c; d1++;
        }
      pz->cursor = d1;
      pthread_mutex_unlock(&pz->mx);

      t0 = pz_now();
      if ((status = pz_run_range(dd, km, cc->d_cnt, d0 * mm->npq, d1 * mm->npq)) != eslOK)
        {
          /* Same retry as the single-model path: halve the budget and redo the
           * range.  Unlike there, the cursor has already been advanced past it
           * -- and it is shared with the other cards' threads -- so put it
           * back under the lock before going round again. */
          if (d1 - d0 <= 1) goto ERROR;
          cc->cps *= 0.5;
          if (cc->cps < P7_PZ_CPS_MIN) cc->cps = P7_PZ_CPS_MIN;
          pthread_mutex_lock(&pz->mx);
          if (pz->cursor > d0) pz->cursor = d0;
          pthread_mutex_unlock(&pz->mx);
          if (pz_verbose())
            fprintf(stderr, "# pezy: dev %d model range [%d,%d) came up short; budget halved\n",
                    dd->idx, d0, d1);
          continue;
        }
      pz_update_cps(&cc->cps, cells, pz_now() - t0);

      /* Model ranges are disjoint, so the cards write disjoint slices of
       * h_part and no lock is needed here. */
      if (pzclEnqueueReadBuffer(dd->q, cc->d_xE, PZCL_TRUE,
                                (size_t) d0 * mm->npq, (size_t) (d1 - d0) * mm->npq,
                                mm->h_part + (size_t) d0 * mm->npq, 0, NULL, NULL) != PZCL_SUCCESS)
        { status = eslFAIL; goto ERROR; }

      wk->cells += cells; wk->nlaunch++;
    }
  return NULL;

 ERROR:
  wk->status = status;
  pthread_mutex_lock(&pz->mx);
  if (pz->failed == eslOK) pz->failed = status;
  pthread_mutex_unlock(&pz->mx);
  return NULL;
}

int
p7_pzmm_SSV(P7_PZ *pz, P7_PZMM *mm, const ESL_DSQ **dsq, const int *L, int nq, uint8_t *xE)
{
  pzmm_work_run_t *wk  = NULL;
  pthread_t       *tid = NULL;
  int64_t  pqres = 0, total = 0;
  int      chunk;
  int      nitem;
  int      i, d, j, c;
  int      status;

  if (nq <= 0) return eslOK;

  p7_pz_seqbatch_Reuse(mm->qb);
  for (i = 0; i < nq; i++)
    if ((status = p7_pz_seqbatch_Add(mm->qb, dsq[i], L[i])) != eslOK) return status;

  chunk = 8 * mm->rowstride; if (chunk < 4096) chunk = 4096;
  if ((status = pzmm_build_pq(mm, nq, chunk)) != eslOK) return status;

  nitem = mm->ndev * mm->npq;
  if (nitem > mm->part_alloc)
    { void *tmp; int na = nitem + nitem/8; ESL_RALLOC(mm->h_part, tmp, (size_t) na); mm->part_alloc = na; }

  for (j = 0; j < mm->npq; j++) pqres += mm->h_pq[2*j + 1];
  for (d = 0; d < mm->ndev; d++)
    total += (int64_t) ((((int) mm->h_mdl[2*d + 1]) + 7) & ~7) * pqres;

  pz->cursor = 0;
  pz->failed = eslOK;

  ESL_ALLOC(wk, sizeof(pzmm_work_run_t) * mm->ncard);
  for (c = 0; c < mm->ncard; c++)
    {
      wk[c].pz = pz; wk[c].dd = &pz->dv[c]; wk[c].mm = mm; wk[c].cc = &mm->cv[c];
      wk[c].nitem = nitem; wk[c].pqres = pqres; wk[c].total = total;
      wk[c].cells = 0; wk[c].nlaunch = 0; wk[c].share = FALSE; wk[c].status = eslOK;
    }

  {
    double tall = pz_now();

    if (mm->ncard == 1) pzmm_run_card(&wk[0]);
    else
      {
        ESL_ALLOC(tid, sizeof(pthread_t) * mm->ncard);
        for (c = 0; c < mm->ncard; c++)
          if (pthread_create(&tid[c], NULL, pzmm_run_card, &wk[c]) != 0) { status = eslESYS; goto ERROR; }
        for (c = 0; c < mm->ncard; c++) pthread_join(tid[c], NULL);
        free(tid); tid = NULL;
      }

    for (c = 0; c < mm->ncard; c++) if (wk[c].status != eslOK) { status = wk[c].status; goto ERROR; }

    if (pz_verbose())
      {
        int nlaunch = 0;
        tall = pz_now() - tall;
        for (c = 0; c < mm->ncard; c++) nlaunch += wk[c].nlaunch;
        fprintf(stderr, "# pezy: %s %d queries, %d pseudo-queries, %d dev, %d launches, "
                        "%.1f Gcells in %.3fs = %.0f Gcells/s",
                wk[0].share ? "ssv_multi_share" : "ssv_multi",
                nq, mm->npq, mm->ncard, nlaunch, total / 1e9, tall, total / tall / 1e9);
        if (mm->ncard > 1)
          for (c = 0; c < mm->ncard; c++) fprintf(stderr, " [d%d %.0f%%]", c, 100.0 * wk[c].cells / (double) total);
        fprintf(stderr, "\n");
      }
  }
  free(wk); wk = NULL;

  /* fold the pseudo-queries back into one xE per (model, query) */
  for (d = 0; d < mm->ndev; d++)
    {
      int      m   = mm->devmodel[d];
      uint8_t *row = xE + (size_t) m * nq;
      const uint8_t *part = mm->h_part + (size_t) d * mm->npq;
      for (i = 0; i < nq; i++) row[i] = 0x80;
      for (j = 0; j < mm->npq; j++)
        if (part[j] > row[mm->pq_owner[j]]) row[mm->pq_owner[j]] = part[j];
    }
  return eslOK;

 ERROR:
  if (tid) free(tid);
  if (wk)  free(wk);
  return status;
}


/*****************************************************************
 * 6. Floating-point conformance probe
 *
 * Diagnostic only; nothing in the pipeline calls this.  It exists to decide
 * whether HMMER's Forward/Backward engines could be parallelized per call
 * (one target sequence per thread) without changing the output.  That needs no
 * reassociation -- forward_engine() is a pure function of (dsq, L, om), and
 * --cpu already runs it concurrently with identical results -- so the whole
 * question reduces to whether the SC3 float unit matches x86 SSE bit for bit.
 * See p7_pz_fpprobe.h for the operation list and the reasoning.
 *
 * The kernel object is created here rather than in pz_dev_open() so that the
 * production path neither pays for it nor fails if the probe is removed.
 *****************************************************************/

int p7_pz_NDev(const P7_PZ *pz) { return pz ? pz->ndev : 0; }

int
p7_pz_FPProbe(P7_PZ *pz, int devidx, const float *in, int n, float *out, char *errbuf)
{
  P7_PZDEV   *dd;
  pzcl_kernel k = NULL;
  pzcl_mem    d_in = NULL, d_out = NULL, d_cnt = NULL;
  size_t      insz  = sizeof(float) * (size_t) P7_PZ_FPPROBE_NIN * (size_t) n;
  size_t      outsz = sizeof(float) * (size_t) P7_PZ_FP_NOP      * (size_t) n;
  int32_t     cnt[P7_PZ_CNT_N];
  int         status = eslFAIL;
  pzcl_int    err;
  size_t      gws;

  if (pz == NULL || devidx < 0 || devidx >= pz->ndev || n <= 0) return eslEINVAL;
  dd  = &pz->dv[devidx];
  gws = dd->gws;

  k = pzclCreateKernel(dd->prog, "fpprobe", &err);
  if (err != PZCL_SUCCESS)
    { if (errbuf) sprintf(errbuf, "pzclCreateKernel(fpprobe) failed (%d)", (int) err); goto DONE; }
  if (pz_set_stack && pz_set_stack(k, P7_PZ_STACKSIZE) != PZCL_SUCCESS)
    { if (errbuf) sprintf(errbuf, "pezy_set_per_thread_stack_size(fpprobe) failed"); goto DONE; }

  d_in  = pzclCreateBuffer(dd->ctx, PZCL_MEM_READ_WRITE, insz,        NULL, &err);
  if (err != PZCL_SUCCESS) { if (errbuf) sprintf(errbuf, "pzclCreateBuffer(in) failed (%d)",  (int) err); goto DONE; }
  d_out = pzclCreateBuffer(dd->ctx, PZCL_MEM_READ_WRITE, outsz,       NULL, &err);
  if (err != PZCL_SUCCESS) { if (errbuf) sprintf(errbuf, "pzclCreateBuffer(out) failed (%d)", (int) err); goto DONE; }
  d_cnt = pzclCreateBuffer(dd->ctx, PZCL_MEM_READ_WRITE, sizeof(cnt), NULL, &err);
  if (err != PZCL_SUCCESS) { if (errbuf) sprintf(errbuf, "pzclCreateBuffer(cnt) failed (%d)", (int) err); goto DONE; }

  memset(cnt, 0, sizeof(cnt));
  if (pzclEnqueueWriteBuffer(dd->q, d_in,  PZCL_FALSE, 0, insz,        in,  0, NULL, NULL) != PZCL_SUCCESS ||
      pzclEnqueueWriteBuffer(dd->q, d_cnt, PZCL_TRUE,  0, sizeof(cnt), cnt, 0, NULL, NULL) != PZCL_SUCCESS)
    { if (errbuf) sprintf(errbuf, "upload failed"); goto DONE; }

  if (pzclSetKernelArg(k, 0, sizeof(pzcl_mem), &d_out) != PZCL_SUCCESS ||
      pzclSetKernelArg(k, 1, sizeof(pzcl_mem), &d_in)  != PZCL_SUCCESS ||
      pzclSetKernelArg(k, 2, sizeof(int),      &n)     != PZCL_SUCCESS ||
      pzclSetKernelArg(k, 3, sizeof(pzcl_mem), &d_cnt) != PZCL_SUCCESS)
    { if (errbuf) sprintf(errbuf, "pzclSetKernelArg failed"); goto DONE; }

  if (pzclEnqueueNDRangeKernel(dd->q, k, 1, NULL, &gws, NULL, 0, NULL, NULL) != PZCL_SUCCESS ||
      pzclFinish(dd->q) != PZCL_SUCCESS)
    { if (errbuf) sprintf(errbuf, "launch failed"); goto DONE; }

  if (pzclEnqueueReadBuffer(dd->q, d_cnt, PZCL_TRUE, 0, sizeof(cnt), cnt, 0, NULL, NULL) != PZCL_SUCCESS ||
      pzclEnqueueReadBuffer(dd->q, d_out, PZCL_TRUE, 0, outsz,       out, 0, NULL, NULL) != PZCL_SUCCESS)
    { if (errbuf) sprintf(errbuf, "download failed"); goto DONE; }

  /* Same completion tally as the filters: every thread must have reported. */
  if (cnt[P7_PZ_CNT_DONE] != (int32_t) gws)
    { if (errbuf) sprintf(errbuf, "kernel did not finish (%d/%d threads)", (int) cnt[P7_PZ_CNT_DONE], (int) gws); goto DONE; }

  status = eslOK;

 DONE:
  if (d_in)  pzclReleaseMemObject(d_in);
  if (d_out) pzclReleaseMemObject(d_out);
  if (d_cnt) pzclReleaseMemObject(d_cnt);
  if (k)     pzclReleaseKernel(k);
  return status;
}

/*****************************************************************
 * 7. Parser-mode Forward
 *
 * One target sequence per thread.  See p7_pz_fbparse.h for what crosses PCIe
 * (only the xmx array, (L+1)*6 floats) and why the arithmetic is bit-exact.
 *
 * Kernel choice
 * -------------
 * The DP row is per-thread local memory, P7_PZ_NSCELLS * Q float4 = 48*Q
 * bytes.  With P7_PZ_LMEM_AVAIL = 16 KB that is Q <= 341 (M <= 1364) at one
 * thread per PE, or Q <= 42 (M <= 168) at eight.  Eight threads per PE hides
 * global-load latency, which matters because tfv and rfv stay in global
 * memory, so we take it whenever the row fits; otherwise one thread per PE.
 * Anything that will not fit either way is the caller's problem: the function
 * returns eslENORESULT and the caller runs p7_ForwardParser() on the host.
 *
 * Backward reuses the very same row -- it starts from scratch at row L -- so
 * fusing the two directions into one kernel costs no extra local memory, and
 * it is what makes Backward nearly free to add: backward_engine() needs the
 * Forward matrix ONLY for its scale factors, and those are sitting in global
 * memory where this item's Forward just wrote them.
 *
 * This deliberately allocates its buffers per call, like p7_pz_FPProbe(),
 * rather than keeping them in P7_PZDEV.  Correctness first; the batch is big
 * enough that one allocation per launch is not what will limit it.
 *****************************************************************/

/* Pack om->xf into the order the kernel expects (p7_pz_fbparse.h).
 * om must already be length-reconfigured for the target: xf depends on L. */
void
p7_pz_FBPackXf(const P7_OPROFILE *om, float *xf8)
{
  xf8[P7_PZ_XF_NMOVE] = om->xf[p7O_N][p7O_MOVE];
  xf8[P7_PZ_XF_NLOOP] = om->xf[p7O_N][p7O_LOOP];
  xf8[P7_PZ_XF_CLOOP] = om->xf[p7O_C][p7O_LOOP];
  xf8[P7_PZ_XF_EMOVE] = om->xf[p7O_E][p7O_MOVE];
  xf8[P7_PZ_XF_ELOOP] = om->xf[p7O_E][p7O_LOOP];
  xf8[P7_PZ_XF_JLOOP] = om->xf[p7O_J][p7O_LOOP];
  xf8[P7_PZ_XF_JMOVE] = om->xf[p7O_J][p7O_MOVE];
  xf8[P7_PZ_XF_CMOVE] = om->xf[p7O_C][p7O_MOVE];
}

/* Threads per PE this model can use, or 0 if the row will not fit at all.
 *
 * Each running thread needs its own DP row in local memory, so what fits is
 * whatever 48*Q divides into P7_PZ_FB_LMEM_AVAIL.  But WIDER IS NOT BETTER --
 * extra threads buy only latency hiding, and past a point they cost more in
 * local-memory and cache pressure than they buy.
 *
 * MEASURED, Gcells/s over Swiss-Prot with targets capped at L=400:
 *
 *   M      tpp=1   tpp=2   tpp=4   tpp=8
 *    35     4.03    4.43    3.82    1.63
 *   101     8.65   10.87   10.72    3.82
 *   167    11.48   15.52   15.77    5.70
 *   278    14.21   20.96    --      --
 *   420    22.07   26.14    --      --
 *
 * Two is best or tied-best everywhere, and eight is 3-4x WORSE than two.  This
 * used to take the widest that fit, which meant every model with M <= 168 ran
 * at a third of the speed it could.  HMMER_PEZY_FBTPP overrides the cap.
 */
#define P7_PZ_FB_TPP_MAX  2
static int
pz_fb_tpp(int M)
{
  int         Q   = P7_PZ_NQF(M);
  size_t      row = sizeof(float) * (size_t) P7_PZ_FBROW_FLOATS(Q);
  const char *v   = getenv("HMMER_PEZY_FBTPP");
  int         cap = P7_PZ_FB_TPP_MAX, t;

  if (v) { int w = atoi(v); if (w >= 1 && w <= 8) cap = w; }
  for (t = 8; t >= 1; t >>= 1)
    if (t <= cap && row * (size_t) t <= (size_t) P7_PZ_FB_LMEM_AVAIL) return t;
  return 0;
}

int
p7_pz_FBModelIsSupported(const P7_OPROFILE *om) { return pz_fb_tpp(om->M) > 0; }

/* Order <idx> so that len[idx[0]] >= len[idx[1]] >= ...
 *
 * A counting sort on the length, because target lengths are small integers
 * (the pipeline rejects anything over 100000) and n runs to hundreds of
 * thousands: qsort on that is a measurable slice of a batch that only takes
 * tens of milliseconds on the device. */
static void
pz_sort_desc_by_len(int *idx, const int *len, int n)
{
  int  maxL = 0, i, v;
  int *cnt  = NULL;

  for (i = 0; i < n; i++) if (len[i] > maxL) maxL = len[i];
  if ((cnt = malloc(sizeof(int) * (size_t) (maxL + 2))) == NULL)
    return;                         /* unsorted is slower, not wrong */
  memset(cnt, 0, sizeof(int) * (size_t) (maxL + 2));

  for (i = 0; i < n; i++) cnt[len[i]]++;
  /* prefix sums from the long end, so the longest land first */
  { int run = 0;
    for (v = maxL; v >= 0; v--) { int c = cnt[v]; cnt[v] = run; run += c; } }
  for (i = 0; i < n; i++) idx[cnt[len[i]]++] = i;

  free(cnt);
}

/* The host tail of p7_ForwardParser(), driven by the xmx the device returned.
 *
 * A transcription of fwdback.c:428-484, and it has to stay one: totscale is a
 * float that the reference accumulates with `+= log(xE)` in increasing i, and
 * log() is not guaranteed to round the same way in every math library -- which
 * is exactly why the device returns the raw scale factors instead of their logs.
 * Rows that were not rescaled hold SCALE == 1.0, and log(1.0) is exactly 0.0,
 * so skipping them is bit-identical to summing them.
 *
 * The same pass screens for the FTZ divergence (see p7_pz_fbparse.h): the host
 * scores in flush-to-zero and the kernel does not reproduce that, so a block
 * holding a value small enough for that to matter is handed back for the host
 * to redo.
 *
 * Returns:
 *   eslOK         - *ret_sc is the score; the device block is trustworthy
 *   eslENORESULT  - too small to trust; caller must run p7_ForwardParser()
 *   eslERANGE     - the three conditions the reference itself throws on
 */
int
p7_pz_FwdParserScore(const float *xmx, int L, const float *xf8, float *ret_sc, float *opt_totscale)
{
  float  totscale = 0.0f;
  float  xC       = xmx[L * P7_PZ_NXCELLS + P7_PZ_XC_C];
  int    i, s;

  for (i = 0; i <= L; i++)
    for (s = 0; s < P7_PZ_NXCELLS; s++)
      {
        float v = xmx[i * P7_PZ_NXCELLS + s];
        if (v > 0.0f && v < P7_PZ_FB_TINY) return eslENORESULT;
      }

  for (i = 1; i <= L; i++)
    {
      float sc = xmx[i * P7_PZ_NXCELLS + P7_PZ_XC_SCALE];
      if (sc != 1.0f) totscale += log(sc);
    }

  if      (isnan(xC))          return eslERANGE;
  else if (L > 0 && xC == 0.0) return eslERANGE;
  else if (isinf(xC) == 1)     return eslERANGE;

  if (opt_totscale) *opt_totscale = totscale;
  if (ret_sc)       *ret_sc = totscale + log(xC * xf8[P7_PZ_XF_CMOVE]);
  return eslOK;
}

/* The host tail of p7_BackwardParser(), driven by the two xmx blocks the
 * device returned.
 *
 * Two things the kernel deliberately does not return:
 *
 *   has_own_scales.  backward_engine() flips it when fwd's scale factors turn
 *     out to be insufficient [J3/119], and from then on writes its own.  So it
 *     is TRUE exactly when some Backward scale differs from the Forward one --
 *     and where they agree the flag changes nothing downstream, since
 *     decoding.c only uses it to multiply by fwd/bck, which is then exactly
 *     1.0.
 *
 *   totscale.  A sum of logs, so it is summed here with the host's own log(),
 *     in the reference's own order: log(scale[L]) first, then descending i,
 *     only where scale > 1.
 *
 * Returns eslOK, or eslENORESULT if the block holds a value small enough that
 * the host's flush-to-zero could have diverged (see p7_pz_fbparse.h), in which
 * case the caller must run p7_BackwardParser() itself.
 */
int
p7_pz_BckParserFinish(const float *fxmx, const float *bxmx, int L,
                      int *ret_own_scales, float *ret_totscale)
{
  float totscale;
  int   i, s, own = FALSE;

  for (i = 0; i <= L; i++)
    for (s = 0; s < P7_PZ_NXCELLS; s++)
      {
        float v = bxmx[i * P7_PZ_NXCELLS + s];
        if (v > 0.0f && v < P7_PZ_FB_TINY) return eslENORESULT;
      }

  for (i = 1; i <= L; i++)
    if (bxmx[i * P7_PZ_NXCELLS + P7_PZ_XC_SCALE] != fxmx[i * P7_PZ_NXCELLS + P7_PZ_XC_SCALE])
      { own = TRUE; break; }

  totscale = log(bxmx[L * P7_PZ_NXCELLS + P7_PZ_XC_SCALE]);
  for (i = L-1; i >= 1; i--)
    {
      float sc = bxmx[i * P7_PZ_NXCELLS + P7_PZ_XC_SCALE];
      if (sc > 1.0f) totscale += log(sc);
    }

  *ret_own_scales = own;
  *ret_totscale   = totscale;
  return eslOK;
}


/* Run parser-mode Forward AND Backward for items [0..b->n) of <b> against <om>,
 * on device
 * <devidx>.
 *
 *   it_xf  [b->n * P7_PZ_XF_N]  per-item specials, from p7_pz_FBPackXf() after
 *                              p7_oprofile_ReconfigLength(om, L)
 *   xoff   [b->n]              offset into <xmx> of item j's block
 *   xmx    out, 2*xn floats    the Forward blocks in [0,xn), the Backward ones
 *                              in [xn,2*xn) at the same offsets
 *
 * Returns eslOK, or eslENORESULT if this model's DP row does not fit in local
 * memory (the caller must then use the host engine).
 */
int
p7_pz_FBParserBatch(P7_PZ *pz, int devidx, const P7_OPROFILE *om,
                     const P7_PZ_SEQBATCH *b, const float *it_xf,
                     const int64_t *xoff, float *xmx, int64_t xn,
                     char *errbuf)
{
  P7_PZDEV   *dd;
  pzcl_kernel k = NULL;
  pzcl_mem    d_tfv = NULL, d_rfv = NULL, d_xf = NULL, d_dsq = NULL;
  pzcl_mem    d_off = NULL, d_len = NULL, d_xmx = NULL, d_cnt = NULL;
  int64_t    *off2  = NULL;
  int32_t    *len2  = NULL;
  float      *xf2   = NULL;
  int        *idx   = NULL;
  int32_t     cnt[P7_PZ_CNT_N];
  int         Q     = P7_PZ_NQF(om->M);
  int         Kp    = om->abc->Kp;
  size_t      tfvsz = sizeof(float) * (size_t) 8 * Q * 4;
  size_t      rfvsz = sizeof(float) * (size_t) Kp * Q * 4;
  size_t      xfsz  = sizeof(float) * (size_t) b->n * P7_PZ_XF_N;
  size_t      offsz = sizeof(int64_t) * (size_t) b->n * 3;
  size_t      lensz = sizeof(int32_t) * (size_t) b->n;
  size_t      xmxsz = sizeof(float) * (size_t) xn * 2;   /* Forward then Backward */
  int         tpp, j, status = eslFAIL;
  pzcl_int    err;
  size_t      gws;

  double      t_kern = 0, t_alloc = 0, t_up = 0, t_run = 0, t_down = 0, t_rel = 0, t0;

  if (pz == NULL || devidx < 0 || devidx >= pz->ndev || b->n <= 0) return eslEINVAL;
  if ((tpp = pz_fb_tpp(om->M)) == 0) return eslENORESULT;

  dd  = &pz->dv[devidx];
  /* The full work size either way: the one-thread-per-PE build launches every
   * thread and has 7 of every 8 quit immediately, exactly as ssv_lmem_pe
   * does.  That keeps the completion tally comparable. */
  gws = dd->gws;

  /* Hand the items out longest first.
   *
   * An item is a whole Forward pass and cannot be split -- the recursion is
   * sequential in i -- so the batch cannot finish before its longest sequence
   * does.  Longest-processing-time-first is the classic fix, and it is what
   * this does; xmx offsets travel with their item, so the output still lands
   * in the caller's order whatever order we run in.
   *
   * MEASURED: it buys nothing on its own, because it only helps once there
   * are more items than running threads.  Below that every item
   * starts immediately and the batch time IS the longest item, whatever order
   * they go out in.  The number that matters is 33.5 us per residue of the
   * longest target, dead linear -- see p7_pz_FBMaxL().  Kept because it does help
   * once batches get large, and it costs a counting sort.
   *
   * it_off is one array so the kernel stays at ten arguments: dsq offsets
   * first, then xmx offsets.  See the argument-count note in p7_pz_kernel.pzc. */
  ESL_ALLOC(idx,  sizeof(int) * (size_t) b->n);
  ESL_ALLOC(off2, offsz);
  ESL_ALLOC(len2, lensz);
  ESL_ALLOC(xf2,  xfsz);
  for (j = 0; j < b->n; j++) idx[j] = j;
  pz_sort_desc_by_len(idx, b->len, b->n);
  for (j = 0; j < b->n; j++)
    {
      int o = idx[j];
      off2[j]            = b->off[o];
      off2[b->n + j]     = xoff[o];        /* Forward block  */
      off2[2 * b->n + j] = xoff[o] + xn;   /* Backward block */
      len2[j]            = (int32_t) b->len[o];
      memcpy(xf2 + (size_t) j * P7_PZ_XF_N, it_xf + (size_t) o * P7_PZ_XF_N,
             sizeof(float) * P7_PZ_XF_N);
    }

  t0 = pz_now();
  { char nm[32];
    snprintf(nm, sizeof(nm), "fbparse_t%d", tpp);
    k = pzclCreateKernel(dd->prog, nm, &err);
    if (err != PZCL_SUCCESS)
      { if (errbuf) sprintf(errbuf, "pzclCreateKernel(%s) failed (%d)", nm, (int) err); goto DONE; } }
  if (pz_set_stack && pz_set_stack(k, P7_PZ_FB_STACKSIZE) != PZCL_SUCCESS)
    { if (errbuf) sprintf(errbuf, "pezy_set_per_thread_stack_size(fbparse) failed"); goto DONE; }
  t_kern = pz_now() - t0;

  t0 = pz_now();
#define PZ_FB_BUF(V, SZ, NAME)                                                        \
  do { V = pzclCreateBuffer(dd->ctx, PZCL_MEM_READ_WRITE, (SZ), NULL, &err);           \
       if (err != PZCL_SUCCESS)                                                        \
         { if (errbuf) sprintf(errbuf, "pzclCreateBuffer(%s, %zu B) failed (%d)",      \
                               NAME, (size_t)(SZ), (int) err); goto DONE; } } while (0)
  PZ_FB_BUF(d_tfv, tfvsz,        "tfv");
  PZ_FB_BUF(d_rfv, rfvsz,        "rfv");
  PZ_FB_BUF(d_xf,  xfsz,         "xf");
  PZ_FB_BUF(d_dsq, (size_t) b->nres, "dsq");
  PZ_FB_BUF(d_off, offsz,        "off");
  PZ_FB_BUF(d_len, lensz,        "len");
  PZ_FB_BUF(d_xmx, xmxsz,        "xmx");
  PZ_FB_BUF(d_cnt, sizeof(cnt),  "cnt");
#undef PZ_FB_BUF
  t_alloc = pz_now() - t0;

  t0 = pz_now();
  memset(cnt, 0, sizeof(cnt));
  /* om->rfv[x] = om->rfv[0] + x*Q and om->tfv is 8*Q contiguous vectors, so
   * both upload verbatim -- the kernel indexes exactly the same layout. */
  if (pzclEnqueueWriteBuffer(dd->q, d_tfv, PZCL_FALSE, 0, tfvsz, (const void *) om->tfv,    0, NULL, NULL) != PZCL_SUCCESS ||
      pzclEnqueueWriteBuffer(dd->q, d_rfv, PZCL_FALSE, 0, rfvsz, (const void *) om->rfv[0], 0, NULL, NULL) != PZCL_SUCCESS ||
      pzclEnqueueWriteBuffer(dd->q, d_xf,  PZCL_FALSE, 0, xfsz,  xf2,                       0, NULL, NULL) != PZCL_SUCCESS ||
      pzclEnqueueWriteBuffer(dd->q, d_dsq, PZCL_FALSE, 0, (size_t) b->nres, b->dsq,         0, NULL, NULL) != PZCL_SUCCESS ||
      pzclEnqueueWriteBuffer(dd->q, d_off, PZCL_FALSE, 0, offsz, off2,                      0, NULL, NULL) != PZCL_SUCCESS ||
      pzclEnqueueWriteBuffer(dd->q, d_len, PZCL_FALSE, 0, lensz, len2,                      0, NULL, NULL) != PZCL_SUCCESS ||
      pzclEnqueueWriteBuffer(dd->q, d_cnt, PZCL_TRUE,  0, sizeof(cnt), cnt,                 0, NULL, NULL) != PZCL_SUCCESS)
    { if (errbuf) sprintf(errbuf, "upload failed"); goto DONE; }
  t_up = pz_now() - t0;

  { int M = om->M, nitem = b->n;
    if (pzclSetKernelArg(k, 0, sizeof(pzcl_mem), &d_tfv) != PZCL_SUCCESS ||
        pzclSetKernelArg(k, 1, sizeof(pzcl_mem), &d_rfv) != PZCL_SUCCESS ||
        pzclSetKernelArg(k, 2, sizeof(pzcl_mem), &d_xf)  != PZCL_SUCCESS ||
        pzclSetKernelArg(k, 3, sizeof(pzcl_mem), &d_dsq) != PZCL_SUCCESS ||
        pzclSetKernelArg(k, 4, sizeof(pzcl_mem), &d_off) != PZCL_SUCCESS ||
        pzclSetKernelArg(k, 5, sizeof(pzcl_mem), &d_len) != PZCL_SUCCESS ||
        pzclSetKernelArg(k, 6, sizeof(int),      &M)     != PZCL_SUCCESS ||
        pzclSetKernelArg(k, 7, sizeof(int),      &nitem) != PZCL_SUCCESS ||
        pzclSetKernelArg(k, 8, sizeof(pzcl_mem), &d_xmx) != PZCL_SUCCESS ||
        pzclSetKernelArg(k, 9, sizeof(pzcl_mem), &d_cnt) != PZCL_SUCCESS)
      { if (errbuf) sprintf(errbuf, "pzclSetKernelArg failed"); goto DONE; }
  }

  t0 = pz_now();
  if (pzclEnqueueNDRangeKernel(dd->q, k, 1, NULL, &gws, NULL, 0, NULL, NULL) != PZCL_SUCCESS ||
      pzclFinish(dd->q) != PZCL_SUCCESS)
    { if (errbuf) sprintf(errbuf, "launch failed"); goto DONE; }
  t_run = pz_now() - t0;

  t0 = pz_now();
  if (pzclEnqueueReadBuffer(dd->q, d_cnt, PZCL_TRUE, 0, sizeof(cnt), cnt, 0, NULL, NULL) != PZCL_SUCCESS ||
      pzclEnqueueReadBuffer(dd->q, d_xmx, PZCL_TRUE, 0, xmxsz,       xmx, 0, NULL, NULL) != PZCL_SUCCESS)
    { if (errbuf) sprintf(errbuf, "download failed"); goto DONE; }
  t_down = pz_now() - t0;

  if (pz_verbose())
    /* MEASURED: kern, alloc, up, down and the releases are all 0.000 even for a
     * batch of a few hundred items -- the per-call buffer allocation this
     * function does is not what costs anything.  Kept in the line because that
     * was worth knowing, and cheap to keep knowing. */
    fprintf(stderr, "# pezy: fbparse dev%d M=%d tpp=%d %d items (maxL %d): kern %.3f, alloc %.3f, "
                    "up %.3f (%.1f MB), run %.3f, down %.3f (%.1f MB)\n",
            devidx, om->M, tpp, b->n, (int) len2[0], t_kern, t_alloc,
            t_up,   (double)(tfvsz+rfvsz+xfsz+(size_t)b->nres+offsz+lensz) / 1048576.0,
            t_run,
            t_down, (double) xmxsz / 1048576.0);

  /* Same completion tally as the filters: a launch that did not run to the
   * end must not be used.
   * Threads that quit early (tid != 0 in the one-per-PE build) still report,
   * so the expected total is the full work size either way. */
  if (cnt[P7_PZ_CNT_DONE] != (int32_t) dd->gws)
    { if (errbuf) sprintf(errbuf, "fbparse did not finish (%d/%d threads)",
                          (int) cnt[P7_PZ_CNT_DONE], (int) dd->gws); goto DONE; }
  if (cnt[0] < b->n)
    { if (errbuf) sprintf(errbuf, "fbparse cursor short (%d/%d items)", (int) cnt[0], b->n); goto DONE; }

  status = eslOK;

 DONE:
  t0 = pz_now();
  if (idx)   free(idx);
  if (off2)  free(off2);
  if (len2)  free(len2);
  if (xf2)   free(xf2);
  if (d_tfv) pzclReleaseMemObject(d_tfv);
  if (d_rfv) pzclReleaseMemObject(d_rfv);
  if (d_xf)  pzclReleaseMemObject(d_xf);
  if (d_dsq) pzclReleaseMemObject(d_dsq);
  if (d_off) pzclReleaseMemObject(d_off);
  if (d_len) pzclReleaseMemObject(d_len);
  if (d_xmx) pzclReleaseMemObject(d_xmx);
  if (d_cnt) pzclReleaseMemObject(d_cnt);
  if (k)     pzclReleaseKernel(k);
  t_rel = pz_now() - t0;
  if (pz_verbose() && t_rel > 0.001)
    fprintf(stderr, "# pezy: fbparse dev%d release %.3fs\n", devidx, t_rel);
  return status;

 ERROR:
  if (idx)  free(idx);
  if (off2) free(off2);
  if (len2) free(len2);
  if (xf2)  free(xf2);
  return eslEMEM;
}

/*****************************************************************
 * 8. Device memory bandwidth probe
 *
 * Diagnostic only; see the kernel-side comment in p7_pz_kernel.pzc section 8
 * for why this exists.  Short version: the O(ML) full-matrix chain would move
 * ~70 GB per hmmsearch query, so whether it belongs on the device is decided
 * by the achieved copy bandwidth -- and that is a runtime property of the
 * card, not a constant.
 *****************************************************************/

/* Runs read, write and copy over a buffer of <mb> MB on device <devidx>, and
 * returns the achieved GB/s of each.  Pass mb well above the device's
 * last-level cache or the cache answers instead of device memory. */
int
p7_pz_BWProbe(P7_PZ *pz, int devidx, int mb, double *ret_read, double *ret_write,
              double *ret_copy, char *errbuf)
{
  P7_PZDEV   *dd;
  pzcl_kernel kr = NULL, kw = NULL, kc = NULL;
  pzcl_mem    d_a = NULL, d_b = NULL, d_out = NULL, d_cnt = NULL;
  size_t      nbytes = (size_t) mb * 1048576;
  int64_t     n4     = (int64_t) (nbytes / 16);
  int32_t     cnt[P7_PZ_CNT_N];
  int         status = eslFAIL;
  pzcl_int    err;
  size_t      gws;
  double      t0;

  if (pz == NULL || devidx < 0 || devidx >= pz->ndev || mb <= 0) return eslEINVAL;
  dd  = &pz->dv[devidx];
  gws = dd->gws;

  kr = pzclCreateKernel(dd->prog, "bw_read",  &err);  if (err != PZCL_SUCCESS) goto DONE;
  kw = pzclCreateKernel(dd->prog, "bw_write", &err);  if (err != PZCL_SUCCESS) goto DONE;
  kc = pzclCreateKernel(dd->prog, "bw_copy",  &err);  if (err != PZCL_SUCCESS) goto DONE;
  if (pz_set_stack) { pz_set_stack(kr, 256); pz_set_stack(kw, 256); pz_set_stack(kc, 256); }

  d_a   = pzclCreateBuffer(dd->ctx, PZCL_MEM_READ_WRITE, nbytes, NULL, &err); if (err != PZCL_SUCCESS) goto DONE;
  d_b   = pzclCreateBuffer(dd->ctx, PZCL_MEM_READ_WRITE, nbytes, NULL, &err); if (err != PZCL_SUCCESS) goto DONE;
  d_out = pzclCreateBuffer(dd->ctx, PZCL_MEM_READ_WRITE, sizeof(float) * gws, NULL, &err); if (err != PZCL_SUCCESS) goto DONE;
  d_cnt = pzclCreateBuffer(dd->ctx, PZCL_MEM_READ_WRITE, sizeof(cnt), NULL, &err); if (err != PZCL_SUCCESS) goto DONE;

#define PZ_BW_RUN(K, LABEL, BYTES, SETARGS)                                             \
  do {                                                                                  \
    memset(cnt, 0, sizeof(cnt));                                                        \
    if (pzclEnqueueWriteBuffer(dd->q, d_cnt, PZCL_TRUE, 0, sizeof(cnt), cnt, 0, NULL, NULL) != PZCL_SUCCESS) \
      { if (errbuf) sprintf(errbuf, "%s: counter upload failed", LABEL); goto DONE; }    \
    SETARGS                                                                             \
    t0 = pz_now();                                                                      \
    if (pzclEnqueueNDRangeKernel(dd->q, (K), 1, NULL, &gws, NULL, 0, NULL, NULL) != PZCL_SUCCESS \
        || pzclFinish(dd->q) != PZCL_SUCCESS)                                            \
      { if (errbuf) sprintf(errbuf, "%s: launch failed", LABEL); goto DONE; }            \
    *ret_##K##_gbs = (double)(BYTES) / (pz_now() - t0) / 1e9;                            \
    if (pzclEnqueueReadBuffer(dd->q, d_cnt, PZCL_TRUE, 0, sizeof(cnt), cnt, 0, NULL, NULL) != PZCL_SUCCESS \
        || cnt[P7_PZ_CNT_DONE] != (int32_t) gws)                                         \
      { if (errbuf) sprintf(errbuf, "%s: kernel did not finish (%d/%d)", LABEL,          \
                            (int) cnt[P7_PZ_CNT_DONE], (int) gws); goto DONE; }          \
  } while (0)

  /* read */
  { double *ret_kr_gbs = ret_read;
    PZ_BW_RUN(kr, "bw_read", nbytes,
      if (pzclSetKernelArg(kr, 0, sizeof(pzcl_mem), &d_a)   != PZCL_SUCCESS ||
          pzclSetKernelArg(kr, 1, sizeof(int64_t),  &n4)    != PZCL_SUCCESS ||
          pzclSetKernelArg(kr, 2, sizeof(pzcl_mem), &d_out) != PZCL_SUCCESS ||
          pzclSetKernelArg(kr, 3, sizeof(pzcl_mem), &d_cnt) != PZCL_SUCCESS)
        { if (errbuf) sprintf(errbuf, "bw_read: setarg failed"); goto DONE; } ); }

  /* write */
  { double *ret_kw_gbs = ret_write;
    PZ_BW_RUN(kw, "bw_write", nbytes,
      if (pzclSetKernelArg(kw, 0, sizeof(pzcl_mem), &d_a)   != PZCL_SUCCESS ||
          pzclSetKernelArg(kw, 1, sizeof(int64_t),  &n4)    != PZCL_SUCCESS ||
          pzclSetKernelArg(kw, 2, sizeof(pzcl_mem), &d_cnt) != PZCL_SUCCESS)
        { if (errbuf) sprintf(errbuf, "bw_write: setarg failed"); goto DONE; } ); }

  /* copy: both halves of the traffic count */
  { double *ret_kc_gbs = ret_copy;
    PZ_BW_RUN(kc, "bw_copy", 2 * nbytes,
      if (pzclSetKernelArg(kc, 0, sizeof(pzcl_mem), &d_a)   != PZCL_SUCCESS ||
          pzclSetKernelArg(kc, 1, sizeof(pzcl_mem), &d_b)   != PZCL_SUCCESS ||
          pzclSetKernelArg(kc, 2, sizeof(int64_t),  &n4)    != PZCL_SUCCESS ||
          pzclSetKernelArg(kc, 3, sizeof(pzcl_mem), &d_cnt) != PZCL_SUCCESS)
        { if (errbuf) sprintf(errbuf, "bw_copy: setarg failed"); goto DONE; } ); }
#undef PZ_BW_RUN

  status = eslOK;

 DONE:
  if (d_a)   pzclReleaseMemObject(d_a);
  if (d_b)   pzclReleaseMemObject(d_b);
  if (d_out) pzclReleaseMemObject(d_out);
  if (d_cnt) pzclReleaseMemObject(d_cnt);
  if (kr)    pzclReleaseKernel(kr);
  if (kw)    pzclReleaseKernel(kw);
  if (kc)    pzclReleaseKernel(kc);
  if (status != eslOK && errbuf && ! errbuf[0]) sprintf(errbuf, "bandwidth probe failed");
  return status;
}
