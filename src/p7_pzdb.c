/* p7_pzdb.c
 *
 * A resident, pre-digitized target sequence database (see p7_pzdb.h).
 *
 * Contents:
 *   1. The fast FASTA digitizer
 *   2. Building the index
 *   3. Serving the pipeline
 */
#include <p7_config.h>

#ifdef HAVE_PEZY

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <pthread.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"
#include "p7_pz.h"
#include "p7_pzdb.h"

/* Codes used only inside our own 256-entry map.  Residues are 0..127, so any
 * value above that is a control code; keeping SKIP < BAD lets one comparison
 * per byte detect anything we are not prepared to handle. */
#define PZDB_SKIP 200      /* whitespace: ' ', '\t', '\r', '\n'   */
#define PZDB_BAD  201      /* anything else: bail out to esl_sqio  */

int64_t
p7_pzdb_MaxBytes(void)
{
  const char *s = getenv("HMMER_PEZY_MAXDB");
  if (s) { long long v = atoll(s); if (v > 0) return (int64_t) v; }
  return (int64_t) 8 * 1024 * 1024 * 1024;   /* 8 GB of digital sequence */
}

/*****************************************************************
 * 1. The fast FASTA digitizer
 *****************************************************************/

/* Build our byte map from the alphabet's, using exactly the classification
 * that inmap_fasta() installs for a digital FASTA read:
 *   residue        -> digital code
 *   ' ' '\t' '\r'  -> ignored
 *   '\n'           -> end of line
 *   '>'            -> end of record   (we locate records ourselves, so a '>'
 *                                      inside a sequence line is a bail-out)
 *   '-' and all else -> illegal
 */
static void
pzdb_build_map(const ESL_ALPHABET *abc, uint8_t map[256])
{
  int x;

  for (x = 0; x < 256; x++) map[x] = PZDB_BAD;
  for (x = 0; x < 128; x++)
    if (abc->inmap[x] <= 127) map[x] = (uint8_t) abc->inmap[x];

  map['-']  = PZDB_BAD;      /* FASTA is an ungapped format */
  map['*']  = (abc->inmap['*'] <= 127) ? (uint8_t) abc->inmap['*'] : PZDB_BAD;
  map[' ']  = PZDB_SKIP;
  map['\t'] = PZDB_SKIP;
  map['\r'] = PZDB_SKIP;
  map['\n'] = PZDB_SKIP;
  map['>']  = PZDB_BAD;
}

/* Digitize [src, src+n) into dst[0..], returning the residue count, or -1 if
 * any byte is one we will not guess at.
 *
 * The whole point of this function: the store is unconditional and the
 * increment is a comparison, so there is no data-dependent branch and no
 * mispredict per residue.  esl_sqio's addbuf() branches on every byte, and
 * measures ~90 MB/s on this machine; this runs at GB/s.
 */
static int64_t
pzdb_digitize(const uint8_t map[256], const char *src, size_t n, uint8_t *dst)
{
  size_t  i;
  int64_t w   = 0;
  uint8_t bad = 0;

  for (i = 0; i < n; i++)
    {
      uint8_t x = map[(uint8_t) src[i]];
      dst[w] = x;
      w   += (x < 128);
      bad |= (x > PZDB_SKIP);
    }
  return bad ? -1 : w;
}

/*****************************************************************
 * 2. Building the index
 *****************************************************************/

static int
pzdb_grow_index(P7_PZDB *db, int need)
{
  void *tmp;
  int   status;
  if (need <= db->nalloc) return eslOK;
  {
    int na = (db->nalloc ? db->nalloc * 2 : 65536);
    while (na < need) na *= 2;
    ESL_RALLOC(db->roff,   tmp, sizeof(int64_t) * na);
    ESL_RALLOC(db->len,    tmp, sizeof(int)     * na);
    ESL_RALLOC(db->dsqoff, tmp, sizeof(int64_t) * na);
    db->nalloc = na;
  }
  return eslOK;
 ERROR:
  return status;
}

void
p7_pzdb_Destroy(P7_PZDB *db)
{
  if (!db) return;
  if (db->fmap && db->fmap != MAP_FAILED) munmap(db->fmap, db->fmap_n);
  if (db->fd >= 0) close(db->fd);
  if (db->roff)    free(db->roff);
  if (db->len)     free(db->len);
  if (db->dsqoff)  free(db->dsqoff);
  if (db->dsq_all) free(db->dsq_all);
  if (db->seqfile) free(db->seqfile);
  free(db);
}

/* Find the next record start at or after <p>: a '>' that begins a line.
 * <p> is always strictly inside the file, so p[-1] is valid. */
static const char *
pzdb_next_record(const char *p, const char *end)
{
  while (p < end)
    {
      const char *g = memchr(p, '>', (size_t) (end - p));
      if (!g) return NULL;
      if (g[-1] == '\n') return g;
      p = g + 1;
    }
  return NULL;
}

/* Phase 2 work: digitize a contiguous run of records.
 *
 * Output placement uses the neat property that a record's digital form is
 * never longer than its bytes in the file: '>' + header + '\n' + L residues
 * (+ newlines) is at least L+2 bytes, which is exactly what dsq[0..L+1] needs.
 * So record s can be written at its own file offset, dsqoff[s] = roff[s], and
 * every thread writes into a disjoint range with no coordination and no
 * prefix sum.  The array is sparse -- the gaps are the headers and newlines --
 * and the gaps are zeroed so that nothing uninitialized is ever uploaded.
 */
typedef struct {
  P7_PZDB     *db;
  const char  *fmap;
  size_t       fmap_n;
  const uint8_t *map;
  int          lo, hi;      /* records [lo,hi) */
  int          bad;         /* set if a byte we will not guess at turned up */
} pzdb_work_t;

static void *
pzdb_digitize_range(void *arg)
{
  pzdb_work_t *wk  = (pzdb_work_t *) arg;
  P7_PZDB     *db  = wk->db;
  const char  *end = wk->fmap + wk->fmap_n;
  int          s;

  for (s = wk->lo; s < wk->hi; s++)
    {
      const char *cur = wk->fmap + db->roff[s];
      const char *nxt = (s + 1 < db->nseq) ? wk->fmap + db->roff[s+1] : end;
      const char *hdr_end, *seq_beg;
      int64_t     L, o = db->roff[s];
      int64_t     gap_lo, gap_hi;

      hdr_end = memchr(cur, '\n', (size_t) (nxt - cur));
      seq_beg = hdr_end ? hdr_end + 1 : nxt;

      db->dsq_all[o] = eslDSQ_SENTINEL;
      L = pzdb_digitize(wk->map, seq_beg, (size_t) (nxt - seq_beg), db->dsq_all + o + 1);
      if (L < 0 || L > 100000) { wk->bad = 1; return NULL; }
      db->dsq_all[o + 1 + L] = eslDSQ_SENTINEL;
      db->len[s]    = (int) L;
      db->dsqoff[s] = o;

      /* zero the header/newline bytes we did not write into */
      gap_lo = o + L + 2;
      gap_hi = (nxt - wk->fmap);
      if (gap_hi > gap_lo) memset(db->dsq_all + gap_lo, 0, (size_t) (gap_hi - gap_lo));
    }
  return NULL;
}

int
p7_pzdb_Build(const char *seqfile, ESL_ALPHABET *abc, int64_t maxbytes,
              int nthreads, P7_PZDB **ret_db, char *errbuf)
{
  P7_PZDB     *db = NULL;
  struct stat  st;
  uint8_t      map[256];
  const char  *cur, *end;
  pzdb_work_t *wk  = NULL;
  pthread_t   *tid = NULL;
  int          t, i;
  int          status;

  *ret_db = NULL;
  if (nthreads < 1) nthreads = 1;
  if (nthreads > 32) nthreads = 32;

  ESL_ALLOC(db, sizeof(P7_PZDB));
  memset(db, 0, sizeof(P7_PZDB));
  db->fd  = -1;
  db->abc = abc;

  if ((db->fd = open(seqfile, O_RDONLY)) < 0)
    { if (errbuf) sprintf(errbuf, "cannot open %s", seqfile); status = eslENOTFOUND; goto ERROR; }
  if (fstat(db->fd, &st) < 0 || ! S_ISREG(st.st_mode))
    { if (errbuf) sprintf(errbuf, "%s is not a regular file", seqfile); status = eslEFORMAT; goto ERROR; }
  if ((int64_t) st.st_size > maxbytes)
    { if (errbuf) sprintf(errbuf, "%s is larger than the resident cap", seqfile); status = eslERANGE; goto ERROR; }
  if (st.st_size < 4)
    { if (errbuf) sprintf(errbuf, "%s is too small", seqfile); status = eslEFORMAT; goto ERROR; }

  db->fmap_n = (size_t) st.st_size;
  db->fmap   = mmap(NULL, db->fmap_n, PROT_READ, MAP_PRIVATE, db->fd, 0);
  if (db->fmap == MAP_FAILED)
    { db->fmap = NULL; if (errbuf) sprintf(errbuf, "cannot mmap %s", seqfile); status = eslESYS; goto ERROR; }
  madvise(db->fmap, db->fmap_n, MADV_WILLNEED);

  if (db->fmap[0] != '>')
    { if (errbuf) sprintf(errbuf, "%s does not start with '>'", seqfile); status = eslEFORMAT; goto ERROR; }

  pzdb_build_map(abc, map);

  /* +2 covers the one degenerate case where a record's digital form is longer
   * than its bytes: a bare ">" as the final record needs two sentinels but
   * occupies one byte. */
  db->dsq_alloc = (int64_t) db->fmap_n + 2;
  ESL_ALLOC(db->dsq_all, (size_t) db->dsq_alloc);

  /* ---- phase 1: locate the records (serial; this is just memchr) ---- */
  cur = db->fmap;
  end = db->fmap + db->fmap_n;
  while (cur != NULL)
    {
      if ((status = pzdb_grow_index(db, db->nseq + 1)) != eslOK) goto ERROR;
      db->roff[db->nseq++] = (int64_t) (cur - db->fmap);
      cur = pzdb_next_record(cur + 1, end);
    }

  /* ---- phase 2: digitize, in parallel over disjoint record ranges ---- */
  ESL_ALLOC(wk, sizeof(pzdb_work_t) * nthreads);
  for (t = 0; t < nthreads; t++)
    {
      wk[t].db     = db;
      wk[t].fmap   = db->fmap;
      wk[t].fmap_n = db->fmap_n;
      wk[t].map    = map;
      wk[t].lo     = (int) ((int64_t) db->nseq *  t      / nthreads);
      wk[t].hi     = (int) ((int64_t) db->nseq * (t + 1) / nthreads);
      wk[t].bad    = 0;
    }
  if (nthreads == 1) pzdb_digitize_range(&wk[0]);
  else
    {
      ESL_ALLOC(tid, sizeof(pthread_t) * nthreads);
      for (t = 0; t < nthreads; t++)
        if (pthread_create(&tid[t], NULL, pzdb_digitize_range, &wk[t]) != 0)
          { status = eslESYS; goto ERROR; }
      for (t = 0; t < nthreads; t++) pthread_join(tid[t], NULL);
      free(tid); tid = NULL;
    }
  for (t = 0; t < nthreads; t++)
    if (wk[t].bad)
      { if (errbuf) sprintf(errbuf, "unexpected character or over-long sequence in %s; using the standard reader", seqfile);
        status = eslEFORMAT; goto ERROR; }
  free(wk); wk = NULL;

  db->dsq_n = db->dsq_alloc;
  for (i = 0; i < db->nseq; i++) {
    db->nres += db->len[i];
    if (db->len[i] > db->maxL) db->maxL = db->len[i];
  }

  if ((db->seqfile = strdup(seqfile)) == NULL) { status = eslEMEM; goto ERROR; }

  /* the mapping is only needed for the scan; hits are re-read through sqio */
  munmap(db->fmap, db->fmap_n);
  db->fmap = NULL;

  *ret_db = db;
  return eslOK;

 ERROR:
  if (tid) free(tid);
  if (wk)  free(wk);
  p7_pzdb_Destroy(db);
  return status;
}

/*****************************************************************
 * 3. Serving the pipeline
 *****************************************************************/

/* The index is already laid out the way the batch wants it -- dsq[0..L+1]
 * blocks with an offset table -- so hand the arrays over rather than copying
 * a few hundred megabytes. */
int
p7_pzdb_ToBatch(const P7_PZDB *db, P7_PZ_SEQBATCH *b)
{
  return p7_pz_seqbatch_Borrow(b, db->dsq_all, db->dsqoff, db->len, db->nseq, db->dsq_n);
}

int
p7_pzdb_FetchSQ(const P7_PZDB *db, ESL_SQFILE *sqfp, int s, ESL_SQ *sq)
{
  int status;
  if (s < 0 || s >= db->nseq) return eslEINVAL;
  if ((status = esl_sqfile_Position(sqfp, (off_t) db->roff[s])) != eslOK) return status;
  if ((status = esl_sqio_Read(sqfp, sq))                       != eslOK) return status;
  return eslOK;
}

#else /*!HAVE_PEZY*/
int p7_pzdb_dummy_symbol_to_avoid_empty_translation_unit = 0;
#endif /*HAVE_PEZY*/
