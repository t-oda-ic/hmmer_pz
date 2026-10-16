/* pz_pzdb_utest.c
 *
 * Verify that P7_PZDB's fast FASTA reader agrees with esl_sqio, sequence by
 * sequence, over a whole database: same count, same lengths, same digital
 * sequences.  Then verify that p7_pzdb_FetchSQ() reproduces the full ESL_SQ
 * (name, accession, description, dsq) for a sample of records.
 *
 * usage: pz_pzdb_utest [-n <nsample>] <seqfile>
 */
#include <p7_config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_stopwatch.h"

#include "hmmer.h"

#ifndef HAVE_PEZY
int main(void) { printf("pz_pzdb_utest: built without PEZY support\n"); return 0; }
#else
#include "p7_pzdb.h"

static ESL_OPTIONS options[] = {
  { "-h", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, NULL, "show brief help",                       0 },
  { "-n", eslARG_INT, "2000", NULL, NULL, NULL, NULL, NULL, "records to re-fetch and fully compare", 0 },
  { "-t", eslARG_INT,    "3", NULL, NULL, NULL, NULL, NULL, "threads for the reader",                0 },
  { 0,0,0,0,0,0,0,0,0,0 },
};
static char usage[]  = "[-options] <seqfile>";
static char banner[] = "test P7_PZDB against esl_sqio";

static int
streq(const char *a, const char *b)
{
  if (a == NULL && b == NULL) return TRUE;
  if (a == NULL || b == NULL) return FALSE;
  return strcmp(a, b) == 0;
}

int
main(int argc, char **argv)
{
  ESL_GETOPTS   *go      = p7_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char          *seqfile = esl_opt_GetArg(go, 1);
  int            nsample = esl_opt_GetInteger(go, "-n");
  ESL_ALPHABET  *abc     = esl_alphabet_Create(eslAMINO);
  P7_PZDB       *db      = NULL;
  ESL_SQFILE    *sqfp    = NULL;
  ESL_SQ        *sq      = NULL;
  ESL_STOPWATCH *w       = esl_stopwatch_Create();
  char           errbuf[256];
  double         t_fast, t_sqio;
  int            i, status;
  long           nbadlen = 0, nbaddsq = 0, nbadfetch = 0;
  int            nesl    = 0;

  esl_stopwatch_Start(w);
  if ((status = p7_pzdb_Build(seqfile, abc, p7_pzdb_MaxBytes(), esl_opt_GetInteger(go, "-t"), &db, errbuf)) != eslOK)
    esl_fatal("p7_pzdb_Build failed (%d): %s", status, errbuf);
  esl_stopwatch_Stop(w);
  t_fast = w->elapsed;

  printf("# %s\n", seqfile);
  printf("# pzdb : %d sequences, %lld residues, maxL %d, %.3f s  (%.0f MB/s of input)\n",
         db->nseq, (long long) db->nres, db->maxL, t_fast,
         (double) db->fmap_n / t_fast / 1e6);

  /* ---- streaming comparison against esl_sqio ---- */
  if (esl_sqfile_OpenDigital(abc, seqfile, eslSQFILE_FASTA, NULL, &sqfp) != eslOK) esl_fatal("cannot open %s", seqfile);
  sq = esl_sq_CreateDigital(abc);

  esl_stopwatch_Start(w);
  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
    {
      if (nesl < db->nseq)
        {
          if ((int) sq->n != db->len[nesl]) {
            if (nbadlen < 5) printf("LEN  mismatch at %d (%s): sqio %d pzdb %d\n", nesl, sq->name, (int) sq->n, db->len[nesl]);
            nbadlen++;
          } else if (memcmp(sq->dsq, db->dsq_all + db->dsqoff[nesl], (size_t) sq->n + 2) != 0) {
            if (nbaddsq < 5) printf("DSQ  mismatch at %d (%s)\n", nesl, sq->name);
            nbaddsq++;
          }
        }
      nesl++;
      esl_sq_Reuse(sq);
    }
  esl_stopwatch_Stop(w);
  t_sqio = w->elapsed;
  printf("# sqio : %d sequences, %.3f s  ->  pzdb is %.1fx faster\n", nesl, t_sqio, t_sqio / t_fast);
  if (nesl != db->nseq) { printf("COUNT mismatch: sqio %d pzdb %d\n", nesl, db->nseq); nbadlen++; }

  /* ---- full ESL_SQ comparison for a sample of records, via FetchSQ ---- */
  {
    ESL_SQ *sq2 = esl_sq_CreateDigital(abc);
    int     step = (nsample > 0 && db->nseq > nsample) ? db->nseq / nsample : 1;

    esl_sqfile_Position(sqfp, 0);
    for (i = 0; i < db->nseq; i += step)
      {
        esl_sq_Reuse(sq);
        esl_sq_Reuse(sq2);
        /* reference: stream to record i */
        if (esl_sqfile_Position(sqfp, (off_t) db->roff[i]) != eslOK ||
            esl_sqio_Read(sqfp, sq)                        != eslOK) { printf("reference read failed at %d\n", i); nbadfetch++; continue; }
        if (p7_pzdb_FetchSQ(db, sqfp, i, sq2)              != eslOK) { printf("FetchSQ failed at %d\n", i);      nbadfetch++; continue; }

        if (! streq(sq->name, sq2->name) || ! streq(sq->desc, sq2->desc) ||
            ! streq(sq->acc,  sq2->acc)  || sq->n != sq2->n ||
            memcmp(sq->dsq, sq2->dsq, (size_t) sq->n + 2) != 0)
          {
            if (nbadfetch < 5) printf("FETCH mismatch at %d: '%s' vs '%s'\n", i, sq->name, sq2->name);
            nbadfetch++;
          }
      }
    esl_sq_Destroy(sq2);
  }

  printf("# length mismatches %ld / dsq mismatches %ld / fetch mismatches %ld\n", nbadlen, nbaddsq, nbadfetch);
  printf("%s\n", (nbadlen || nbaddsq || nbadfetch) ? "FAILED" : "ok. P7_PZDB agrees with esl_sqio.");

  esl_sq_Destroy(sq);
  esl_sqfile_Close(sqfp);
  p7_pzdb_Destroy(db);
  esl_stopwatch_Destroy(w);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return (nbadlen || nbaddsq || nbadfetch) ? 1 : 0;
}
#endif /*HAVE_PEZY*/
