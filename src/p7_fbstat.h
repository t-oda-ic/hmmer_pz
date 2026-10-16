/* p7_fbstat.h
 *
 * Lightweight diagnostic counters.  For each stage around Forward/Backward,
 * tallies the number of calls, DP cells, and elapsed time.
 *
 * Purpose: measure how much of the wall time the O(ML) full-matrix path
 * actually takes, as input to deciding whether that path is worth moving to
 * the device (in chunks, since a whole matrix per thread does not fit).
 *
 * Off by default (costs one branch per call).  HMMER_FBSTAT=1 turns it on,
 * and a table is printed to stderr at process exit.
 */
#ifndef P7_FBSTAT_INCLUDED
#define P7_FBSTAT_INCLUDED

#include <stdint.h>

enum p7_fbstat_slot_e {
  P7_FBS_FWD_PARSE = 0,   /* p7_ForwardParser   : O(M+L) */
  P7_FBS_FWD_FULL,        /* p7_Forward         : O(ML)  */
  P7_FBS_BCK_PARSE,       /* p7_BackwardParser  : O(M+L) */
  P7_FBS_BCK_FULL,        /* p7_Backward        : O(ML)  */
  P7_FBS_DECODING,        /* p7_Decoding        : O(ML)  */
  P7_FBS_DOMDECODE,       /* p7_DomainDecoding  : O(L)   */
  P7_FBS_OPTACC,          /* p7_OptimalAccuracy : O(ML)  */
  P7_FBS_OATRACE,         /* p7_OATrace         : O(L)   traceback */
  P7_FBS_NULL2_EXP,       /* p7_Null2_ByExpectation : O(ML) */
  P7_FBS_NULL2_TRACE,     /* p7_Null2_ByTrace   : depends on the trace */
  P7_FBS_STOTRACE,        /* p7_StochasticTrace : stochastic */
  P7_FBS_VITFILTER,       /* p7_ViterbiFilter   : O(ML) int16 */
  P7_FBS_NSLOTS
};

extern int p7_fbstat_on;

extern uint64_t p7_fbstat_now(void);
extern void     p7_fbstat_add(int slot, int M, int L, uint64_t ns);

#define P7_FBSTAT_T0(v)                 uint64_t v = p7_fbstat_on ? p7_fbstat_now() : 0
#define P7_FBSTAT_T1(v, slot, M, L)                                     \
  do { if (p7_fbstat_on) p7_fbstat_add((slot), (int)(M), (int)(L),      \
                                       p7_fbstat_now() - (v)); } while (0)

#endif /*P7_FBSTAT_INCLUDED*/
