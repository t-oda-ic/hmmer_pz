/* p7_fbstat.c
 *
 * Lightweight diagnostic counters (see p7_fbstat.h).
 * Only collects when HMMER_FBSTAT=1, and prints a table to stderr at process exit.
 */
#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "p7_fbstat.h"

int p7_fbstat_on = 0;

#define P7_FBS_NBUCKET 16       /* log2 histogram of L */

static struct {
  uint64_t calls;
  uint64_t cells;               /* sum of M*L */
  uint64_t ns;
  uint64_t sumL;
  uint64_t maxL;
  uint64_t maxM;
  uint64_t maxbytes;            /* largest size in bytes of one matrix requested by a single call */
  uint64_t lhist[P7_FBS_NBUCKET];
} stats[P7_FBS_NSLOTS];

static const char *slotname[P7_FBS_NSLOTS] = {
  "p7_ForwardParser   O(M+L)",
  "p7_Forward         O(ML) ",
  "p7_BackwardParser  O(M+L)",
  "p7_Backward        O(ML) ",
  "p7_Decoding        O(ML) ",
  "p7_DomainDecoding  O(L)  ",
  "p7_OptimalAccuracy O(ML) ",
  "p7_OATrace         O(L)  ",
  "p7_Null2_ByExpect  O(ML) ",
  "p7_Null2_ByTrace   trace ",
  "p7_StochasticTrace rand  ",
  "p7_ViterbiFilter   O(ML) ",
};

/* whether the stage actually allocates an O(ML) matrix (counted in the memory estimate) */
static const int is_full[P7_FBS_NSLOTS] = { 0,1,0,1,1,0,1,0,1,0,0,0 };

uint64_t
p7_fbstat_now(void)
{
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (uint64_t) ts.tv_sec * 1000000000ull + (uint64_t) ts.tv_nsec;
}

void
p7_fbstat_add(int slot, int M, int L, uint64_t ns)
{
  uint64_t bytes;
  int b, i;

  if (slot < 0 || slot >= P7_FBS_NSLOTS) return;
  if (M < 0) M = 0;
  if (L < 0) L = 0;

  /* one row = 3 states (M/I/D) x ceil(M/4) vectors x 16 B; a full matrix has L+1 rows */
  bytes = (uint64_t) 3 * (uint64_t) ((M + 3) / 4 > 0 ? (M + 3) / 4 : 1) * 16;
  if (is_full[slot]) bytes *= (uint64_t) (L + 1);

  __atomic_fetch_add(&stats[slot].calls, 1,                          __ATOMIC_RELAXED);
  __atomic_fetch_add(&stats[slot].cells, (uint64_t) M * (uint64_t) L, __ATOMIC_RELAXED);
  __atomic_fetch_add(&stats[slot].ns,    ns,                         __ATOMIC_RELAXED);
  __atomic_fetch_add(&stats[slot].sumL,  (uint64_t) L,               __ATOMIC_RELAXED);

  for (b = 0, i = L; i > 1 && b < P7_FBS_NBUCKET - 1; i >>= 1) b++;
  __atomic_fetch_add(&stats[slot].lhist[b], 1, __ATOMIC_RELAXED);

  /* the maxima need not be exact (diagnostic only), so no CAS loop */
  if ((uint64_t) L     > stats[slot].maxL)     stats[slot].maxL     = (uint64_t) L;
  if ((uint64_t) M     > stats[slot].maxM)     stats[slot].maxM     = (uint64_t) M;
  if (bytes            > stats[slot].maxbytes) stats[slot].maxbytes = bytes;
}

static void
fbstat_report(void)
{
  uint64_t total_ns = 0, full_ns = 0, parse_ns = 0;
  int s, b;

  if (! p7_fbstat_on) return;

  for (s = 0; s < P7_FBS_NSLOTS; s++) {
    total_ns += stats[s].ns;
    if (is_full[s]) full_ns += stats[s].ns;
  }
  parse_ns = stats[P7_FBS_FWD_PARSE].ns + stats[P7_FBS_BCK_PARSE].ns;
  if (total_ns == 0) return;

  fprintf(stderr, "\n#### HMMER_FBSTAT (pid %ld) ####\n", (long) getpid());
  fprintf(stderr, "%-26s %12s %14s %9s %8s %7s %7s %10s\n",
          "stage", "calls", "cells(M*L)", "sec", "%meas", "meanL", "maxL", "mat/call");
  for (s = 0; s < P7_FBS_NSLOTS; s++) {
    if (stats[s].calls == 0) continue;
    fprintf(stderr, "%-26s %12llu %14llu %9.3f %7.1f%% %7.0f %7llu %8.1fKB\n",
            slotname[s],
            (unsigned long long) stats[s].calls,
            (unsigned long long) stats[s].cells,
            stats[s].ns / 1e9,
            100.0 * stats[s].ns / total_ns,
            stats[s].calls ? (double) stats[s].sumL / stats[s].calls : 0.0,
            (unsigned long long) stats[s].maxL,
            stats[s].maxbytes / 1024.0);
  }
  fprintf(stderr, "%-26s %12s %14s %9.3f %7.1f%%\n", "  [O(ML) full-matrix sum]", "", "",
          full_ns / 1e9, 100.0 * full_ns / total_ns);
  fprintf(stderr, "%-26s %12s %14s %9.3f %7.1f%%\n", "  [O(M+L) parser sum]", "", "",
          parse_ns / 1e9, 100.0 * parse_ns / total_ns);
  fprintf(stderr, "%-26s %9.3f sec\n", "  [all measured stages]", total_ns / 1e9);

  /* L distribution of the full-matrix path (= rows per thread if it ran on the device) */
  fprintf(stderr, "\n# p7_Forward(full) L distribution [log2 buckets]\n");
  for (b = 0; b < P7_FBS_NBUCKET; b++)
    if (stats[P7_FBS_FWD_FULL].lhist[b])
      fprintf(stderr, "#   L in [%6d,%6d) : %llu\n",
              b ? 1 << b : 0, 1 << (b + 1),
              (unsigned long long) stats[P7_FBS_FWD_FULL].lhist[b]);
  fprintf(stderr, "#### end HMMER_FBSTAT ####\n");
}

__attribute__((constructor))
static void
fbstat_init(void)
{
  char *s = getenv("HMMER_FBSTAT");
  if (s != NULL && *s != '\0' && *s != '0') {
    p7_fbstat_on = 1;
    atexit(fbstat_report);
  }
}
