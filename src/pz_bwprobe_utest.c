/* pz_bwprobe_utest.c
 *
 * Report the memory bandwidth this device actually achieves right now.
 *
 * usage: pz_bwprobe_utest [-m <MB>] [-d <dev>] [-r <reps>]
 *
 * Not a test of HMMER; a measurement the port needs in order to decide whether
 * the O(ML) full-matrix chain of p7_domaindef.c belongs on the device.  That
 * chain makes ~8 passes over a 12-byte-per-cell DP matrix -- roughly 70 GB per
 * hmmsearch query -- and the host affords it only because one domain's pair of
 * matrices (~1.3 MB) times the host's threads fits in the CPU's last-level
 * cache.  With a pair per device thread the working set is several GB, far
 * past the device's caches, so it goes to device memory, and the achieved
 * bandwidth there is a runtime property of the card.
 */
#include <p7_config.h>

#include <stdio.h>
#include <stdlib.h>

#include "easel.h"
#include "esl_getopts.h"
#include "hmmer.h"

#ifndef HAVE_PEZY
int main(void) { printf("pz_bwprobe_utest: built without PEZY support\n"); return 0; }
#else
#include "p7_pz.h"

static ESL_OPTIONS options[] = {
  { "-h", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, NULL, "show brief help",                    0 },
  { "-m", eslARG_INT,  "512", NULL, NULL, NULL, NULL, NULL, "buffer size in MB (>> cache size)",  0 },
  { "-d", eslARG_INT,  "-1",  NULL, NULL, NULL, NULL, NULL, "device index (-1 = every device)",   0 },
  { "-r", eslARG_INT,  "3",   NULL, NULL, NULL, NULL, NULL, "repeats; the best of each is shown", 0 },
  { 0,0,0,0,0,0,0,0,0,0 },
};
static char usage[]  = "[-options]";
static char banner[] = "measure the achieved memory bandwidth of the PEZY device";

int
main(int argc, char **argv)
{
  ESL_GETOPTS *go   = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  int          mb   = esl_opt_GetInteger(go, "-m");
  int          only = esl_opt_GetInteger(go, "-d");
  int          reps = esl_opt_GetInteger(go, "-r");
  P7_PZ       *pz   = NULL;
  char         errbuf[256] = "";
  int          d, r, ndev;

  if (p7_pz_DeviceCount() <= 0) esl_fatal("no PEZY device available");
  if ((pz = p7_pz_Create(0, errbuf)) == NULL) esl_fatal("p7_pz_Create failed: %s", errbuf);
  ndev = p7_pz_NDev(pz);

  printf("# %d MB buffer, best of %d\n", mb, reps);
  printf("# %-4s %10s %10s %10s\n", "dev", "read", "write", "copy");
  for (d = 0; d < ndev; d++)
    {
      double br = 0.0, bw = 0.0, bc = 0.0;
      if (only >= 0 && d != only) continue;
      for (r = 0; r < reps; r++)
        {
          double rd, wr, cp;
          errbuf[0] = '\0';
          if (p7_pz_BWProbe(pz, d, mb, &rd, &wr, &cp, errbuf) != eslOK)
            esl_fatal("p7_pz_BWProbe(dev %d) failed: %s", d, errbuf);
          if (rd > br) br = rd;
          if (wr > bw) bw = wr;
          if (cp > bc) bc = cp;
        }
      printf("  %-4d %7.1f GB/s %5.1f GB/s %5.1f GB/s\n", d, br, bw, bc);
    }

  p7_pz_Destroy(pz);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*HAVE_PEZY*/
