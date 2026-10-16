/* p7_pz_fpprobe.h
 *
 * The operation list for the floating-point conformance probe, shared by the
 * device kernel (p7_pz_kernel.pzc) and the host test (pz_fpprobe_utest.c) so
 * the two cannot drift apart.
 *
 * WHY THIS EXISTS
 * ===============
 * HMMER's Forward/Backward engines (impl_sse/fwdback.c) are the largest
 * remaining share of the hmmsearch profile -- 32% + 24% -- and they are pure
 * float.  Parallelizing them *per call* (one target sequence per thread)
 * needs no reassociation at all: forward_engine() is a pure function of
 * (dsq, L, om), and HMMER already runs it concurrently under --cpu with
 * identical output.  So the only thing that can change the answer is the
 * arithmetic itself: does a PEZY-SC float unit reproduce x86 SSE bit for bit
 * for the operations fwdback.c performs?
 *
 * This probe answers that question directly.  Each op below is computed on the
 * device and on the host from the same inputs, and the raw 32-bit patterns are
 * compared.  The interesting ones:
 *
 *   MADD_OFF / MADD_FAST   fwdback.c:312-316 is _mm_mul_ps then _mm_add_ps:
 *                          two roundings.  If the compiler contracts them
 *                          into a fused multiply-add, the low bits change.
 *                          These two ops pin the compiler's behaviour down
 *                          with #pragma clang fp contract, so we learn both
 *                          what the default does and what the FMA would cost
 *                          us in accuracy.
 *   MADD_ON / MAD_CALL     the other two ways to ask for a multiply-add: the
 *                          llvm.fmuladd intrinsic, which the backend is free
 *                          to fuse, and the PZSDK's mad().  Both are compared
 *                          against two separate roundings.
 *   FMA_CALL               the PZSDK's fma(), compared against the host's
 *                          fmaf().  This is the control: if it matches, a
 *                          single-rounding multiply-add is available, and the
 *                          fact that the MADD_* ops match the two-rounding
 *                          reference means the compiler is genuinely not
 *                          contracting.
 *   RECIP_D / RECIP_F      fwdback.c:424 is _mm_set1_ps(1.0 / xE) -- the 1.0
 *                          is a *double*, so it divides in double and rounds
 *                          once.  fwdback.c:420 (xN / xE) is a float divide.
 *                          A port that writes 1.0f/xE for both is wrong.
 *   FWDSTEP                the five-operation body of the inner loop.
 *   LOGF / LOGD            fwdback.c:432 (totscale += log(xE)) and the final
 *                          score.  Library log() implementations need not
 *                          agree to the last bit, so this may differ -- but it
 *                          is avoidable, by returning the scale factors and
 *                          taking the logs on the host in i order.
 *   CHAIN64                64 dependent mul-adds, to show whether a
 *                          single-ulp difference stays put or compounds.
 */
#ifndef p7_PZ_FPPROBE_INCLUDED
#define p7_PZ_FPPROBE_INCLUDED

/* X(name, human-readable description) */
#define P7_PZ_FPPROBE_OPS(X)                                                  \
  X(ADD,       "a + b")                                                       \
  X(SUB,       "a - b")                                                       \
  X(MUL,       "a * b")                                                       \
  X(DIV,       "a / b")                                                       \
  X(MADD_OFF,  "a*b + c          contract(off): two roundings")               \
  X(MADD_FAST, "a*b + c          contract(fast): fused if the HW has it")     \
  X(MADD_ON,   "a*b + c          contract(on): llvm.fmuladd intrinsic")       \
  X(MAD_CALL,  "mad(a,b,c)       PZSDK mad()")                                \
  X(FMA_CALL,  "fma(a,b,c)       PZSDK fma(); vs host fmaf")                  \
  X(RECIP_D,   "(float)(1.0 / (double)a)   <- fwdback.c:424")                 \
  X(RECIP_F,   "1.0f / a")                                                    \
  X(FWDSTEP,   "forward_engine inner loop, one lane (4 madd + 1 mul)")        \
  X(DMADD,     "(float)((double)a*b + c)   double path")                      \
  X(LOGF,      "log(a)           float")                                      \
  X(LOGD,      "(float)log((double)a)")                                       \
  X(CHAIN64,   "64x r = r*b + c  drift accumulation")

#define P7_PZ_FPPROBE_ENUM(NAME, DESC)  P7_PZ_FP_##NAME,
enum { P7_PZ_FPPROBE_OPS(P7_PZ_FPPROBE_ENUM) P7_PZ_FP_NOP };
#undef P7_PZ_FPPROBE_ENUM

/* input layout: in[0..n) = a, in[n..2n) = b, in[2n..3n) = c
 * output layout: out[op*n + i]
 */
#define P7_PZ_FPPROBE_NIN  3

#endif /*p7_PZ_FPPROBE_INCLUDED*/
