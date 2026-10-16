/* p7_pz_fbparse.h
 *
 * Layouts shared between the host and the device for the parser-mode
 * Forward/Backward offload.  Included by BOTH p7_pz.c and p7_pz_kernel.pzc, so
 * it must stay free of easel, hmmer and PZCL headers.
 *
 * WHY parser mode is cheap to ship
 * ================================
 * impl_sse/fwdback.c's forward_engine(do_full=FALSE) writes its DP row to
 * ox->dpf[0] every residue -- one row, reused.  Nothing downstream ever reads
 * it back:
 *
 *   - p7_BackwardParser() takes the Forward matrix only for its scale factors
 *     (fwdback.c: `fwd->xmx[L*p7X_NXCELLS+p7X_SCALE]`).
 *   - p7_DomainDecoding() (impl_sse/decoding.c:160) reads oxf->xmx, oxb->xmx
 *     and om->xf, and nothing else.
 *
 * So the whole result of a parser-mode pass is the xmx array: (L+1) rows of
 * p7X_NXCELLS floats.  At L=300 that is 7.2 KB, versus 1.5 MB for a full
 * matrix.  The DP row stays on the device and never crosses PCIe.
 *
 * WHY IT IS BIT-EXACT
 * ===================
 * forward_engine() is a pure function of (dsq, L, om), so one target sequence
 * per thread executes the same operations in the same order as the reference;
 * no reduction is reassociated.  What remains is the arithmetic itself, and
 * pz_fpprobe_utest checks mul/add/div on the device against x86 SSE bit for
 * bit over four input regimes.  Two consequences for this port:
 *
 *   - log() is a library function, not a single rounded operation, so it is
 *     not expected to match and it stays on the host: the device returns the
 *     raw per-row scale factors and the host sums their logs in increasing i,
 *     exactly as `ox->totscale += log(xE)` did.
 *   - the host runs in FTZ+DAZ (impl_sse.h impl_Init()) and the kernel does
 *     not reproduce that.  This one is real and it does bite: xN decays by
 *     xf[N][LOOP] every row and is divided by xE at every rescale, so after a
 *     few hundred residues it underflows.  Measured against p7_ForwardParser()
 *     over Swiss-Prot: ~0.9% of targets came back with a denormal in the N
 *     column where the host had an exact zero.
 *
 *     Emulating FTZ in the kernel means flushing after every multiply AND
 *     every add (the reference's `_mm_add_ps(sv, _mm_mul_ps(..))` rounds
 *     twice, and FTZ applies to both), which roughly doubles the operation
 *     count of the inner loop.  So instead the kernel runs unflushed and the
 *     HOST checks what came back: p7_pz_FwdParserScore() rejects any block
 *     holding a positive value below P7_PZ_FB_TINY, and the caller redoes
 *     those few targets with p7_ForwardParser().  The check costs one pass
 *     over data the host reads anyway.
 *
 *     Why a host-side check on xmx is enough to catch it.  Every value here is
 *     a non-negative probability, so sums only grow: a denormal partial sum
 *     d < 2^-126 can only change the rounded result of d + t when the result
 *     is below about 2^-102, and every such result reaches xmx either directly
 *     (the specials chain, which is where it actually happens) or through xE,
 *     which is a sum of the whole row and so is only perturbed at the 2^-102
 *     level.  P7_PZ_FB_TINY is set at 2^-100 to leave margin, and the claim is
 *     checked empirically by pz_fwdparse_utest over the whole database.
 */
#ifndef p7_PZ_FBPARSE_INCLUDED
#define p7_PZ_FBPARSE_INCLUDED

/* The special-state transitions forward_engine() uses.  om->xf changes with
 * the target length (p7_oprofile_ReconfigLength), so unlike tfv/rfv these are
 * per item, not per model. */
enum p7_pz_xf_e {
  P7_PZ_XF_NMOVE = 0,   /* om->xf[p7O_N][p7O_MOVE] */
  P7_PZ_XF_NLOOP = 1,   /* om->xf[p7O_N][p7O_LOOP] */
  P7_PZ_XF_CLOOP = 2,   /* om->xf[p7O_C][p7O_LOOP] */
  P7_PZ_XF_EMOVE = 3,   /* om->xf[p7O_E][p7O_MOVE] */
  P7_PZ_XF_ELOOP = 4,   /* om->xf[p7O_E][p7O_LOOP] */
  P7_PZ_XF_JLOOP = 5,   /* om->xf[p7O_J][p7O_LOOP] */
  P7_PZ_XF_JMOVE = 6,   /* om->xf[p7O_J][p7O_MOVE] */
  P7_PZ_XF_CMOVE = 7    /* om->xf[p7O_C][p7O_MOVE]; host tail only */
};
#define P7_PZ_XF_N  8

/* Mirrors of the reference's own layout constants, so the kernel does not have
 * to include impl_sse.h.  Kept as separate names on purpose: if HMMER ever
 * changes one, the static asserts in p7_pz.c fail loudly. */
#define P7_PZ_NXCELLS  6   /* p7X_NXCELLS: E N J B C SCALE   */
#define P7_PZ_XC_E     0
#define P7_PZ_XC_N     1
#define P7_PZ_XC_J     2
#define P7_PZ_XC_B     3
#define P7_PZ_XC_C     4
#define P7_PZ_XC_SCALE 5

#define P7_PZ_NSCELLS  3   /* p7X_NSCELLS: M D I, each a float4 */
#define P7_PZ_SC_M     0
#define P7_PZ_SC_D     1
#define P7_PZ_SC_I     2

/* p7O_NQF(M): vectors per DP row */
#define P7_PZ_NQF(M)  ( (((M)-1)/4 + 1) > 2 ? (((M)-1)/4 + 1) : 2 )

/* Floats of local memory one thread's DP row needs. */
#define P7_PZ_FBROW_FLOATS(Q)  ((Q) * P7_PZ_NSCELLS * 4)

/* Below this, a denormal the host would have flushed could have changed the
 * value, so the target is recomputed on the host.  See the FTZ discussion
 * above.  2^-100. */
#define P7_PZ_FB_TINY  7.888609052210118e-31f

#endif /* p7_PZ_FBPARSE_INCLUDED */
