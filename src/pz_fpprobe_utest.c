/* pz_fpprobe_utest.c
 *
 * Does the PEZY-SC float unit reproduce x86 SSE bit for bit, for the
 * operations HMMER's Forward/Backward engines perform?
 *
 * usage: pz_fpprobe_utest [-n <ntrial>] [-d <dev>] [-s <seed>] [-v]
 *
 * WHY THIS TEST EXISTS
 * ====================
 * forward_engine()/backward_engine() (impl_sse/fwdback.c) are 56% of the
 * hmmsearch profile and are pure float, which makes them look un-offloadable.
 * That concern is about reassociating the inner loop -- the xEv accumulation
 * is a *sum*, and changing its order changes the bits.
 *
 * Parallelizing per *call* has no such problem.  forward_engine() is a pure
 * function of (dsq, L, om); one target sequence per thread performs exactly the
 * reference operation sequence, in the reference order.  HMMER already does
 * this on the host: --cpu N runs the pipeline concurrently over sequences and
 * the output is identical.  So the only remaining question is whether the
 * *arithmetic* agrees, and that is a property of the hardware, not of the
 * decomposition.
 *
 * This test answers it by running the same operations on both sides and
 * comparing raw 32-bit patterns.  The host reference uses the same SSE
 * intrinsics fwdback.c uses, so it is the real reference and not a scalar
 * approximation of it.
 *
 * INPUT REGIMES
 * =============
 * Bit-exactness questions are decided by the operands, so the ops are run over
 * four different populations:
 *
 *   forward   log-uniform positive in [1e-12, 1e4].  This is what Forward
 *             actually holds: scaled odds ratios, kept in range by the sparse
 *             rescaling at fwdback.c:415-434.
 *   normal    uniform random mantissa and exponent over the finite normals.
 *             Exercises rounding cases the forward regime rarely reaches.
 *   denormal  exponents pushed low enough that products and sums land in the
 *             denormal range.  x86 SSE keeps denormals unless MXCSR.FTZ is
 *             set; a device that flushes them to zero would diverge here and
 *             nowhere else, which is the worst kind of bug to find late.
 *   tiny      a and b near 1.0 with random low mantissa bits, c small: makes
 *             a*b+c land on or next to a rounding tie, which is where a fused
 *             multiply-add differs from two separate roundings.
 */
#include <p7_config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_random.h"

#include "hmmer.h"

#ifndef HAVE_PEZY
int main(void) { printf("pz_fpprobe_utest: built without PEZY support\n"); return 0; }
#else
#include <xmmintrin.h>
#include <emmintrin.h>

#include "p7_pz.h"
#include "p7_pz_fpprobe.h"

#define P7_PZ_FPPROBE_NAME(NAME, DESC)  #NAME,
static const char *op_name[] = { P7_PZ_FPPROBE_OPS(P7_PZ_FPPROBE_NAME) };
#undef P7_PZ_FPPROBE_NAME

#define P7_PZ_FPPROBE_DESC(NAME, DESC)  DESC,
static const char *op_desc[] = { P7_PZ_FPPROBE_OPS(P7_PZ_FPPROBE_DESC) };
#undef P7_PZ_FPPROBE_DESC

static ESL_OPTIONS options[] = {
  /* name     type          def      env  range toggles reqs incomp  help                                  docgroup*/
  { "-h",     eslARG_NONE,  FALSE,   NULL, NULL, NULL, NULL, NULL, "show brief help",                          0 },
  { "-n",     eslARG_INT,   "1000000",NULL,"n>0",NULL, NULL, NULL, "input triples per regime",                 0 },
  { "-d",     eslARG_INT,   "0",     NULL, NULL, NULL, NULL, NULL, "PEZY device index",                        0 },
  { "-s",     eslARG_INT,   "42",    NULL, NULL, NULL, NULL, NULL, "RNG seed",                                 0 },
  { "-v",     eslARG_NONE,  FALSE,   NULL, NULL, NULL, NULL, NULL, "show example mismatches",                  0 },
  { 0,0,0,0,0,0,0,0,0,0 },
};
static char usage[]  = "[-options]";
static char banner[] = "probe PEZY-SC float conformance against x86 SSE";

/* bit-level view of a float, for reporting */
static uint32_t
bits(float f) { uint32_t u; memcpy(&u, &f, 4); return u; }

/* signed distance in ulps, saturating.  Only meaningful for same-sign finite
 * values, which is all we compare. */
static long
ulps(float a, float b)
{
  int32_t ia, ib;
  if (isnan(a) || isnan(b)) return (isnan(a) && isnan(b)) ? 0 : 1L << 30;
  memcpy(&ia, &a, 4); memcpy(&ib, &b, 4);
  if (ia < 0) ia = (int32_t) 0x80000000 - ia;
  if (ib < 0) ib = (int32_t) 0x80000000 - ib;
  return (long) ia - (long) ib;
}

/* Scalar helpers that go through SSE, so the host reference performs the same
 * instructions fwdback.c does rather than whatever the x87/scalar path would.
 */
static float sse_add(float a, float b) { float r; _mm_store_ss(&r, _mm_add_ss(_mm_set_ss(a), _mm_set_ss(b))); return r; }
static float sse_sub(float a, float b) { float r; _mm_store_ss(&r, _mm_sub_ss(_mm_set_ss(a), _mm_set_ss(b))); return r; }
static float sse_mul(float a, float b) { float r; _mm_store_ss(&r, _mm_mul_ss(_mm_set_ss(a), _mm_set_ss(b))); return r; }
static float sse_div(float a, float b) { float r; _mm_store_ss(&r, _mm_div_ss(_mm_set_ss(a), _mm_set_ss(b))); return r; }

/* The host's answer for every op, computed the way fwdback.c would. */
static void
host_ops(float a, float b, float c, float *o)
{
  float r;
  int   j;

  o[P7_PZ_FP_ADD]       = sse_add(a, b);
  o[P7_PZ_FP_SUB]       = sse_sub(a, b);
  o[P7_PZ_FP_MUL]       = sse_mul(a, b);
  o[P7_PZ_FP_DIV]       = sse_div(a, b);

  /* two roundings, as _mm_mul_ps then _mm_add_ps gives.  Both MADD ops are
   * compared against this: the reference has no FMA. */
  o[P7_PZ_FP_MADD_OFF]  = sse_add(sse_mul(a, b), c);
  o[P7_PZ_FP_MADD_FAST] = o[P7_PZ_FP_MADD_OFF];
  o[P7_PZ_FP_MADD_ON]   = o[P7_PZ_FP_MADD_OFF];
  o[P7_PZ_FP_MAD_CALL]  = o[P7_PZ_FP_MADD_OFF];
  /* the control: a genuinely fused multiply-add, single rounding */
  o[P7_PZ_FP_FMA_CALL]  = fmaf(a, b, c);

  o[P7_PZ_FP_RECIP_D]   = (float) (1.0 / (double) a);   /* fwdback.c:424 */
  o[P7_PZ_FP_RECIP_F]   = sse_div(1.0f, a);

  r = sse_mul(a, b);
  r = sse_add(r, sse_mul(b, c));
  r = sse_add(r, sse_mul(c, a));
  r = sse_add(r, sse_mul(a, a));
  o[P7_PZ_FP_FWDSTEP]   = sse_mul(r, b);

  o[P7_PZ_FP_DMADD]     = (float) ((double) a * (double) b + (double) c);
  o[P7_PZ_FP_LOGF]      = logf(a);
  o[P7_PZ_FP_LOGD]      = (float) log((double) a);

  r = a;
  for (j = 0; j < 64; j++) r = sse_add(sse_mul(r, b), c);
  o[P7_PZ_FP_CHAIN64]   = r;
}

/* --- input regimes ---------------------------------------------------- */

enum { REG_FORWARD = 0, REG_NORMAL, REG_DENORMAL, REG_TINY, REG_N };
static const char *reg_name[REG_N] = { "forward", "normal", "denormal", "tiny" };
static const char *reg_desc[REG_N] = {
  "log-uniform positive 1e-12..1e4 (what Forward holds)",
  "random mantissa and exponent over the finite normals",
  "products and sums land in the denormal range",
  "a,b near 1.0; a*b+c sits on or beside a rounding tie",
};

/* a random float with the given exponent range and a fully random mantissa */
static float
rand_float(ESL_RANDOMNESS *r, int emin, int emax)
{
  uint32_t m = (uint32_t) (esl_random(r) * 8388608.0);          /* 23 bits */
  uint32_t e = (uint32_t) (emin + (int) (esl_random(r) * (emax - emin + 1)));
  uint32_t u = (e << 23) | (m & 0x7FFFFF);
  float    f;
  memcpy(&f, &u, 4);
  return f;
}

static void
fill_regime(ESL_RANDOMNESS *r, int regime, float *in, int n)
{
  float *A = in, *B = in + n, *C = in + 2 * n;
  int    i;

  for (i = 0; i < n; i++)
    {
      switch (regime) {
      case REG_FORWARD:
        /* exp(uniform(-27.6, 9.2)) = 1e-12 .. 1e4 */
        A[i] = (float) exp(-27.6 + 36.8 * esl_random(r));
        B[i] = (float) exp(-27.6 + 36.8 * esl_random(r));
        C[i] = (float) exp(-27.6 + 36.8 * esl_random(r));
        break;

      case REG_NORMAL:
        /* biased exponent 1..254 = the whole finite normal range */
        A[i] = rand_float(r, 1, 254);
        B[i] = rand_float(r, 1, 254);
        C[i] = rand_float(r, 1, 254);
        break;

      case REG_DENORMAL:
        /* a, b in [2^-70, 2^-60): their product lands in [2^-140, 2^-120),
         * i.e. below the smallest normal (2^-126) and above the smallest
         * denormal (2^-149), so it is a denormal and not zero.  c is itself a
         * denormal, so the add sees one on both sides.
         *
         * (Biased exponents, so 2^-70 is 127-70 = 57.  Getting this wrong is
         * how the first version of this test generated NaNs instead.) */
        A[i] = rand_float(r, 127 - 70, 127 - 60);
        B[i] = rand_float(r, 127 - 70, 127 - 60);
        { uint32_t u = (uint32_t) (esl_random(r) * 8388606.0) + 1;   /* exp 0, nonzero mantissa */
          memcpy(&C[i], &u, 4); }
        break;

      case REG_TINY:
      default:
        /* a,b in [1,2) with random low bits: a*b has 46 significant bits, so
         * the product's rounding is nontrivial and c can push it onto a tie. */
        A[i] = rand_float(r, 127, 127);
        B[i] = rand_float(r, 127, 127);
        C[i] = rand_float(r, 127 - 25, 127 - 22);
        break;
      }
    }
}

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go     = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  int             n      = esl_opt_GetInteger(go, "-n");
  int             devidx = esl_opt_GetInteger(go, "-d");
  int             seed   = esl_opt_GetInteger(go, "-s");
  int             be_verbose = esl_opt_GetBoolean(go, "-v");
  ESL_RANDOMNESS *rng    = NULL;
  P7_PZ          *pz     = NULL;
  float          *in = NULL, *dev = NULL, *hst = NULL;
  char            errbuf[256] = "";
  long            bad[REG_N][P7_PZ_FP_NOP];
  long            maxulp[REG_N][P7_PZ_FP_NOP];
  /* controls -- see the report below for why they are needed */
  long            nfma[REG_N];   /* triples where a true FMA differs from two roundings */
  long            nden[REG_N];   /* triples where an operand or the product is denormal */
  int             regime, op, i, status;
  int             nfail = 0;

  if (p7_pz_DeviceCount() == 0) { printf("pz_fpprobe_utest: no PEZY device available\n"); esl_getopts_Destroy(go); return 0; }

  rng = esl_randomness_Create(seed);
  if ((pz = p7_pz_Create(devidx, errbuf)) == NULL) esl_fatal("p7_pz_Create failed: %s", errbuf);

  ESL_ALLOC(in,  sizeof(float) * (size_t) P7_PZ_FPPROBE_NIN * (size_t) n);
  ESL_ALLOC(dev, sizeof(float) * (size_t) P7_PZ_FP_NOP      * (size_t) n);
  ESL_ALLOC(hst, sizeof(float) * (size_t) P7_PZ_FP_NOP);

  printf("# PEZY-SC float conformance vs. x86 SSE\n");
  printf("# %d cards visible, using device %d; %d triples per regime, seed %d\n",
         p7_pz_DeviceCount(), devidx, n, seed);
  printf("#\n");

  memset(bad, 0, sizeof(bad));
  for (regime = 0; regime < REG_N; regime++)
    for (op = 0; op < P7_PZ_FP_NOP; op++) maxulp[regime][op] = 0;
  memset(nfma,  0, sizeof(nfma));
  memset(nden,  0, sizeof(nden));

  for (regime = 0; regime < REG_N; regime++)
    {
      fill_regime(rng, regime, in, n);

      if ((status = p7_pz_FPProbe(pz, 0, in, n, dev, errbuf)) != eslOK)
        esl_fatal("p7_pz_FPProbe failed: %s", errbuf);

      for (i = 0; i < n; i++)
        {
          float a = in[i], b = in[n + i], c = in[2 * n + i];

          /* Control 1: how many of these triples can tell a fused
           * multiply-add apart from two separate roundings at all?  Without
           * this, "MADD_OFF matched" and "FMA_CALL matched" could both be
           * true simply because the inputs never distinguish the two, and the
           * whole test would prove nothing. */
          if (bits(fmaf(a, b, c)) != bits(sse_add(sse_mul(a, b), c))) nfma[regime]++;

          /* Control 2: is the denormal regime actually producing denormals?
           * The first version of this test got the biased-exponent arithmetic
           * wrong and generated NaNs instead, which looked like a pile of
           * mismatches that were really just differing NaN payloads. */
          if (fpclassify(a) == FP_SUBNORMAL || fpclassify(b) == FP_SUBNORMAL ||
              fpclassify(c) == FP_SUBNORMAL || fpclassify(sse_mul(a, b)) == FP_SUBNORMAL)
            nden[regime]++;

          host_ops(a, b, c, hst);
          for (op = 0; op < P7_PZ_FP_NOP; op++)
            {
              float d = dev[(size_t) op * n + i];
              float h = hst[op];
              if (bits(d) == bits(h)) continue;
              /* +-0.0 differ in bits but are the same value; and two NaNs are
               * the same answer as far as HMMER is concerned -- fwdback.c only
               * ever asks isnan().  Neither is a divergence. */
              if (d == 0.0f && h == 0.0f) continue;
              if (isnan(d) && isnan(h))   continue;
              {
                long u = ulps(d, h);
                if (u < 0) u = -u;
                if (u > maxulp[regime][op]) maxulp[regime][op] = u;
              }
              if (bad[regime][op]++ == 0 && be_verbose)
                printf("#   %-8s %-10s a=%.9g b=%.9g c=%.9g  dev=%.9g [%08x]  host=%.9g [%08x]\n",
                       reg_name[regime], op_name[op], in[i], in[n + i], in[2 * n + i],
                       d, bits(d), h, bits(h));
            }
        }
    }

  /* ---- report ---- */
  printf("# mismatching triples out of %d, by operation and input regime\n#\n", n);
  printf("%-11s %10s %10s %10s %10s   %s\n", "op", reg_name[0], reg_name[1], reg_name[2], reg_name[3], "operation");
  printf("%-11s %10s %10s %10s %10s   %s\n", "-----------", "----------", "----------", "----------", "----------", "---------");
  for (op = 0; op < P7_PZ_FP_NOP; op++)
    {
      printf("%-11s", op_name[op]);
      for (regime = 0; regime < REG_N; regime++)
        {
          if (bad[regime][op] == 0) printf(" %10s", "-");
          else                      printf(" %10ld", bad[regime][op]);
        }
      printf("   %s\n", op_desc[op]);
    }

  printf("#\n# worst |ulp| where they differ (0 = identical everywhere)\n#\n");
  printf("%-11s %10s %10s %10s %10s\n", "op", reg_name[0], reg_name[1], reg_name[2], reg_name[3]);
  for (op = 0; op < P7_PZ_FP_NOP; op++)
    {
      printf("%-11s", op_name[op]);
      for (regime = 0; regime < REG_N; regime++) printf(" %10ld", maxulp[regime][op]);
      printf("\n");
    }

  printf("#\n# controls (a clean table above is only meaningful if these are large)\n#\n");
  printf("%-11s %10s %10s %10s %10s   %s\n", "", reg_name[0], reg_name[1], reg_name[2], reg_name[3], "");
  printf("%-11s", "FMA!=mul+add");
  for (regime = 0; regime < REG_N; regime++) printf(" %10ld", nfma[regime]);
  printf("   triples where a true FMA gives a different answer\n");
  printf("%-11s", "denormal");
  for (regime = 0; regime < REG_N; regime++) printf(" %10ld", nden[regime]);
  printf("   triples with a denormal operand or product\n");

  printf("#\n# input regimes\n");
  for (regime = 0; regime < REG_N; regime++) printf("#   %-9s %s\n", reg_name[regime], reg_desc[regime]);

  /* The verdict.  MADD_FAST is *expected* to differ if the hardware has an
   * FMA -- that is the point of asking -- so it does not count as a failure;
   * neither do the two log()s, which are library functions and which the
   * port can keep on the host anyway.  Everything else must be exact for a
   * per-call port of forward_engine() to be bit-identical. */
  printf("#\n");
  for (op = 0; op < P7_PZ_FP_NOP; op++)
    {
      long tot = 0;
      if (op == P7_PZ_FP_MADD_FAST || op == P7_PZ_FP_MADD_ON  ||
          op == P7_PZ_FP_MAD_CALL   || op == P7_PZ_FP_FMA_CALL ||
          op == P7_PZ_FP_LOGF       || op == P7_PZ_FP_LOGD) continue;
      for (regime = 0; regime < REG_N; regime++) tot += bad[regime][op];
      if (tot) { printf("FAIL %s: %ld mismatches\n", op_name[op], tot); nfail++; }
    }
  if (nfail == 0)
    printf("ok - every operation forward_engine() performs is bit-identical to x86 SSE\n"
           "     (excluding log(), which stays on the host, and the deliberate FMA probes)\n");
  else
    printf("%d operation(s) diverge; a per-call port would not be bit-identical\n", nfail);

  free(in); free(dev); free(hst);
  p7_pz_Destroy(pz);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return nfail ? 1 : 0;

 ERROR:
  esl_fatal("allocation failed");
}
#endif /*HAVE_PEZY*/
