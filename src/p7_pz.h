/* p7_pz.h
 *
 * PEZY-SC offload for HMMER's acceleration pipeline.
 *
 * Only the integer-arithmetic filter stages are offloaded (currently the SSV
 * filter), so device results are bit-identical to the SSE reference
 * implementation.  All floating-point scoring stays on the host.
 */
#ifndef p7_PZ_INCLUDED
#define p7_PZ_INCLUDED

#include <p7_config.h>

#ifdef HAVE_PEZY

#include <stdint.h>
#include "easel.h"
#include "esl_alphabet.h"
#include "hmmer.h"

typedef struct p7_pz_s P7_PZ;

/* A host-side staging area for a batch of digital sequences, packed into one
 * contiguous buffer so it can be uploaded in a single transfer.
 */
typedef struct {
  uint8_t *dsq;        /* packed dsq[0..L] blocks, sentinel included         */
  int64_t *off;        /* [0..n-1] start of sequence i inside <dsq>          */
  int     *len;        /* [0..n-1] L of sequence i                           */
  int      n;          /* number of sequences held                          */
  int      nalloc;
  int64_t  nres;       /* bytes of <dsq> in use                             */
  int64_t  resalloc;
  int      borrowed;   /* TRUE: the three arrays belong to someone else      */
} P7_PZ_SEQBATCH;

/* Point a batch at arrays somebody else owns, instead of copying into it.
 * P7_PZDB already stores its sequences in exactly this layout, and it is
 * hundreds of megabytes, so the copy is worth avoiding. */
extern int p7_pz_seqbatch_Borrow(P7_PZ_SEQBATCH *b, uint8_t *dsq, int64_t *off, int *len,
                                 int n, int64_t nres);

/* runtime availability -------------------------------------------------- */
extern int      p7_pz_DeviceCount(void);
extern int      p7_pz_IsEnabled(void);
extern void     p7_pz_SetEnabled(int on);

/* context --------------------------------------------------------------- */
extern P7_PZ   *p7_pz_Create(int devidx, char *errbuf);
extern void     p7_pz_Destroy(P7_PZ *pz);

/* model ----------------------------------------------------------------- */
/* Upload om's de-striped SSV score table.
 * Returns eslOK, or eslENORESULT if this model cannot be handled on the
 * device (M too large for the DP row to fit in local memory); the caller
 * must then fall back to the CPU filter.
 */
extern int      p7_pz_SetModel(P7_PZ *pz, const P7_OPROFILE *om);
extern int      p7_pz_ModelIsSupported(const P7_OPROFILE *om);

/* sequence batch -------------------------------------------------------- */
extern P7_PZ_SEQBATCH *p7_pz_seqbatch_Create(int nalloc, int64_t resalloc);
extern void            p7_pz_seqbatch_Destroy(P7_PZ_SEQBATCH *b);
extern void            p7_pz_seqbatch_Reuse(P7_PZ_SEQBATCH *b);
extern int             p7_pz_seqbatch_Add(P7_PZ_SEQBATCH *b, const ESL_DSQ *dsq, int L);

/* the filter ------------------------------------------------------------ */
/* Run the SSV filter for every sequence in <b> against the model most
 * recently given to p7_pz_SetModel().  <xE> must have room for b->n bytes and
 * receives the raw-byte xE exactly as impl_sse/ssvfilter.c's get_xE() would
 * have returned it.
 */
extern int      p7_pz_SSVBatch(P7_PZ *pz, const P7_PZ_SEQBATCH *b, uint8_t *xE);

/* Host-side tail of p7_SSVFilter(), driven by a precomputed xE.  This is a
 * verbatim transcription of the code in impl_sse/ssvfilter.c, so that the
 * float arithmetic is unchanged.
 */
extern int      p7_pz_SSVOutcome(uint8_t xE, const P7_OPROFILE *om, float *ret_sc);

/* Build the de-striped SSV score table for <om>: sb[x*(M+1)+k] for k=1..M.
 * Exposed for the unit test.  Caller frees.
 */
extern int8_t  *p7_pz_DestripeSSV(const P7_OPROFILE *om);

/* Build the two device score tables from <om>: the positive and negative
 * halves of the de-striped SSV scores, each Kp rows of <*ret_Mpad> bytes with
 * row x byte (k-1) holding position k.  The kernel works in u = s+128 space
 * where the recursion is an unsigned saturating subtract of the positive half
 * followed by a saturating add of the negative half, so it never has to unpack
 * a sign.  Caller frees.
 */
extern uint8_t *p7_pz_BuildSSVTables(const P7_OPROFILE *om, int *ret_Mpad);


/*****************************************************************
 * Parser-mode Forward (see p7_pz_fbparse.h)
 *
 * One target sequence per thread.  The whole result is the xmx array,
 * (L+1)*p7X_NXCELLS floats -- the DP row never leaves the device.
 *****************************************************************/

/* TRUE if this model's DP row fits in device local memory.  Otherwise the
 * caller must use p7_ForwardParser() on the host. */
extern int  p7_pz_FBModelIsSupported(const P7_OPROFILE *om);

/* Pack om->xf into the kernel's order.  om must already be length-reconfigured
 * for the target: xf depends on L, unlike tfv and rfv. */
extern void p7_pz_FBPackXf(const P7_OPROFILE *om, float *xf8);

/* Run parser-mode Forward AND Backward for every item of <b>, in one launch.
 * <xoff> gives each item's offset into <xmx>, which must have room for
 * 2*<xn> floats: the Forward blocks in [0,xn) and the Backward ones in
 * [xn,2*xn) at the same offsets.  Returns eslOK, or eslENORESULT if the model
 * is too large. */
extern int  p7_pz_FBParserBatch(P7_PZ *pz, int devidx, const P7_OPROFILE *om,
                                const P7_PZ_SEQBATCH *b, const float *it_xf,
                                const int64_t *xoff, float *xmx, int64_t xn,
                                char *errbuf);
/* Where p7_pz_FBParserBatch() puts item <j>'s Backward block, given the
 * Forward offset it was handed and the <xn> of that call. */
#define P7_PZ_BXOFF(xoff, xn)  ((xoff) + (xn))

/* The host tail of p7_ForwardParser(): sums log(scale) in increasing i and
 * finishes the score.  log() stays on the host so that it is the very same
 * log() the reference uses.  Returns eslOK, or eslERANGE on the same three
 * conditions the reference throws on. */
extern int  p7_pz_FwdParserScore(const float *xmx, int L, const float *xf8, float *ret_sc,
                                 float *opt_totscale);

/* The host tail of p7_BackwardParser(): derives has_own_scales and totscale
 * from the two scale columns, and screens the same way. */
extern int  p7_pz_BckParserFinish(const float *fxmx, const float *bxmx, int L,
                                  int *ret_own_scales, float *ret_totscale);


/* Device memory bandwidth probe (diagnostic).  Runs read / write / copy over a
 * buffer of <mb> MB and returns the achieved GB/s of each.  Pass mb well above
 * the device's last-level cache or the cache answers instead of device memory.
 * Exists because whether the O(ML) full-matrix chain belongs on the device is
 * decided by this number, and it is a runtime property of the card, not a
 * constant. */
extern int  p7_pz_BWProbe(P7_PZ *pz, int devidx, int mb, double *ret_read,
                          double *ret_write, double *ret_copy, char *errbuf);


/* Floating-point conformance probe (diagnostic; see p7_pz_fpprobe.h).
 * Runs every operation in the probe list on device <devidx> over <n> input
 * triples.  <in> is 3*n floats (a, then b, then c); <out> must have room for
 * P7_PZ_FP_NOP*n floats and receives op j at [j*n .. j*n+n).
 */
extern int  p7_pz_FPProbe(P7_PZ *pz, int devidx, const float *in, int n, float *out, char *errbuf);
extern int  p7_pz_NDev(const P7_PZ *pz);


/*****************************************************************
 * Many models x many queries (hmmscan's shape)
 *
 * P7_PZMM holds the SSV score tables for a whole profile database on the
 * device, built once, and scores every (model, query) pair in one launch.
 *****************************************************************/

typedef struct p7_pzmm_s P7_PZMM;

/* Build and upload tables for <nmodel> profiles.  <nthreads> is used for the
 * one-off table build only.  Profiles too long for a device DP row are left on
 * the host; ask with p7_pzmm_ModelOnDevice().
 */
extern int  p7_pzmm_Create(P7_PZ *pz, P7_OPROFILE **list, int nmodel, int nthreads,
                           P7_PZMM **ret_mm, char *errbuf);
extern void p7_pzmm_Destroy(P7_PZMM *mm);
extern int  p7_pzmm_ModelOnDevice(const P7_PZMM *mm, int m);
extern int  p7_pzmm_NDevice(const P7_PZMM *mm);

/* Score every (model, query) pair.  <xE> must have room for nmodel*nq bytes and
 * receives the raw-byte xE at [m*nq + q]; entries for models not on the device
 * are left untouched.
 */
extern int  p7_pzmm_SSV(P7_PZ *pz, P7_PZMM *mm, const ESL_DSQ **dsq, const int *L, int nq,
                        uint8_t *xE);

#endif /* HAVE_PEZY */
#endif /* p7_PZ_INCLUDED */
