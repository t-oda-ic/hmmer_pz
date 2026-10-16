/* p7_pzdb.h
 *
 * A resident, pre-digitized target sequence database for the PEZY-SC offload.
 *
 * HMMER's normal path re-reads and re-digitizes the whole FASTA file on every
 * pass over the target database, byte at a time (esl_sqio_ascii.c's seebuf()
 * and addbuf()).  On the benchmark that is ~60% of all CPU time and, because
 * it happens in one reader thread, it is the wall-clock critical path -- so
 * offloading the DP filters alone cannot show up in wall time.
 *
 * P7_PZDB reads the file once with a branch-free digitizer, keeps the digital
 * sequences in memory, and uploads them to the device once.  Every subsequent
 * pass (jackhmmer's iterations, a multi-model query file) is then free on the
 * host side.
 *
 * Only what the pipeline needs before the MSV filter is served from here:
 * the length L of every sequence, and the digital sequence itself.  A
 * sequence that survives the filter is materialized through the ordinary
 * esl_sqio reader, positioned at the record's file offset, so its name,
 * accession, description and dsq are produced by exactly the same code as
 * before.
 */
#ifndef p7_PZDB_INCLUDED
#define p7_PZDB_INCLUDED

#include <p7_config.h>

#ifdef HAVE_PEZY

#include <stdint.h>
#include <sys/types.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"
#include "p7_pz.h"

typedef struct {
  char         *seqfile;
  ESL_ALPHABET *abc;

  int64_t  *roff;         /* [nseq] file offset of the '>' that starts record s */
  int      *len;          /* [nseq] number of residues in record s              */
  int64_t  *dsqoff;       /* [nseq] offset of record s's dsq[0] in <dsq_all>    */
  int       nseq;
  int       nalloc;
  int64_t   nres;

  uint8_t  *dsq_all;      /* dsq[0..L+1] blocks, both sentinels included        */
  int64_t   dsq_n;
  int64_t   dsq_alloc;

  int       maxL;         /* longest record                                     */

  char     *fmap;         /* the mmap'ed file                                   */
  size_t    fmap_n;
  int       fd;
} P7_PZDB;

/* Build the index.  Returns eslOK, or:
 *   eslENOTFOUND / eslESYS  the file cannot be opened or mapped
 *   eslEFORMAT              not something the fast reader will touch
 *                           (not FASTA, or a character it will not guess at)
 *   eslERANGE               larger than <maxbytes>
 * On anything other than eslOK the caller must fall back to the normal reader.
 */
extern int  p7_pzdb_Build(const char *seqfile, ESL_ALPHABET *abc, int64_t maxbytes,
                          int nthreads, P7_PZDB **ret_db, char *errbuf);
extern void p7_pzdb_Destroy(P7_PZDB *db);

/* Fill the P7_PZ sequence batch with the whole database, then upload it.
 * Cheap to call repeatedly: the batch is built once and reused. */
extern int  p7_pzdb_ToBatch(const P7_PZDB *db, P7_PZ_SEQBATCH *b);

/* Materialize record <s> into <sq>, exactly as esl_sqio_Read() would have.
 * <sqfp> must be a positionable handle on the same file, private to the
 * calling thread. */
extern int  p7_pzdb_FetchSQ(const P7_PZDB *db, ESL_SQFILE *sqfp, int s, ESL_SQ *sq);

/* Default resident-size cap, overridable with HMMER_PEZY_MAXDB (bytes). */
extern int64_t p7_pzdb_MaxBytes(void);

#endif /*HAVE_PEZY*/
#endif /*p7_PZDB_INCLUDED*/
