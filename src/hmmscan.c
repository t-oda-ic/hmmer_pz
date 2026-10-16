/* hmmscan: search sequence(s) against a profile HMM database
 *
 * To do:
 *  - in MPI mode, add a check to make sure ncpus >= 2. If 1, then we
 *    only have a master, no workers. See Infernal commit r3972 on the
 *    same point; and same note in hmmsearch.c's to do list.
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

#ifdef HMMER_MPI
#include "mpi.h"
#include "esl_mpi.h"
#endif 

#ifdef HMMER_THREADS
#include <unistd.h>
#include "esl_threads.h"
#include "esl_workqueue.h"
#endif

#include "esl_gumbel.h"
#include "hmmer.h"
#include "p7_hmmcache.h"
#ifdef HAVE_PEZY
#include "p7_pz.h"
#endif
#include <sys/stat.h>
#include <pthread.h>

typedef struct {
#ifdef HMMER_THREADS
  ESL_WORK_QUEUE   *queue;
#endif
  ESL_SQ           *qsq;
  P7_BG            *bg;	         /* null model                              */
  P7_PIPELINE      *pli;         /* work pipeline                           */
  P7_TOPHITS       *th;          /* top hit results                         */
} WORKER_INFO;

typedef struct {
  float  **rows;      /* [n * Kp] row pointers into <data>  */
  float   *data;      /* [n * Kp * fM] odds ratios          */
  int      Kp, fM;
  uint32_t n;
} BIAS_CACHE;

/* What cached_loop() needs to know about the device, if it is in use. */
struct scan_dev_s {
#ifdef HAVE_PEZY
  P7_PZMM       *pzmm;
  const uint8_t *xE;
  int            nq, qi;
#else
  int            unused;
#endif
};


#define REPOPTS     "-E,-T,--cut_ga,--cut_nc,--cut_tc"
#define DOMREPOPTS  "--domE,--domT,--cut_ga,--cut_nc,--cut_tc"
#define INCOPTS     "--incE,--incT,--cut_ga,--cut_nc,--cut_tc"
#define INCDOMOPTS  "--incdomE,--incdomT,--cut_ga,--cut_nc,--cut_tc"
#define THRESHOPTS  "-E,-T,--domE,--domT,--incE,--incT,--incdomE,--incdomT,--cut_ga,--cut_nc,--cut_tc"

#if defined (HMMER_THREADS) && defined (HMMER_MPI)
#define CPUOPTS     "--mpi"
#define MPIOPTS     "--cpu"
#else
#define CPUOPTS     NULL
#define MPIOPTS     NULL
#endif

static ESL_OPTIONS options[] = {
  /* name           type          default  env  range toggles  reqs   incomp                         help                                           docgroup*/
  { "-h",           eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "show brief help on version and usage",                          1 },
  /* Control of output */
  { "-o",           eslARG_OUTFILE, NULL, NULL, NULL,    NULL,  NULL,  NULL,            "direct output to file <f>, not stdout",                         2 },
  { "--tblout",     eslARG_OUTFILE, NULL, NULL, NULL,    NULL,  NULL,  NULL,            "save parseable table of per-sequence hits to file <f>",         2 },
  { "--domtblout",  eslARG_OUTFILE, NULL, NULL, NULL,    NULL,  NULL,  NULL,            "save parseable table of per-domain hits to file <f>",           2 },
  { "--pfamtblout", eslARG_OUTFILE, NULL, NULL, NULL,    NULL,  NULL,  NULL,            "save table of hits and domains to file, in Pfam format <f>",    2 },
  { "--acc",        eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "prefer accessions over names in output",                        2 },
  { "--noali",      eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "don't output alignments, so output is smaller",                 2 },
  { "--notextw",    eslARG_NONE,    NULL, NULL, NULL,    NULL,  NULL, "--textw",        "unlimit ASCII text output line width",                          2 },
  { "--textw",      eslARG_INT,    "120", NULL, "n>=120",NULL,  NULL, "--notextw",      "set max width of ASCII text output lines",                      2 },
  /* Control of reporting thresholds */
  { "-E",           eslARG_REAL,  "10.0", NULL, "x>0",   NULL,  NULL,  REPOPTS,         "report models <= this E-value threshold in output",             4 },
  { "-T",           eslARG_REAL,   FALSE, NULL, NULL,    NULL,  NULL,  REPOPTS,         "report models >= this score threshold in output",               4 },
  { "--domE",       eslARG_REAL,  "10.0", NULL, "x>0",   NULL,  NULL,  DOMREPOPTS,      "report domains <= this E-value threshold in output",            4 },
  { "--domT",       eslARG_REAL,   FALSE, NULL, NULL,    NULL,  NULL,  DOMREPOPTS,      "report domains >= this score cutoff in output",                 4 },
  /* Control of inclusion (significance) thresholds: */
  { "--incE",       eslARG_REAL,  "0.01", NULL, "x>0",   NULL,  NULL,  INCOPTS,         "consider models <= this E-value threshold as significant",      5 },
  { "--incT",       eslARG_REAL,   FALSE, NULL, NULL,    NULL,  NULL,  INCOPTS,         "consider models >= this score threshold as significant",        5 },
  { "--incdomE",    eslARG_REAL,  "0.01", NULL, "x>0",   NULL,  NULL,  INCDOMOPTS,      "consider domains <= this E-value threshold as significant",     5 },
  { "--incdomT",    eslARG_REAL,   FALSE, NULL, NULL,    NULL,  NULL,  INCDOMOPTS,      "consider domains >= this score threshold as significant",       5 },
  /* Model-specific thresholding for both reporting and inclusion */
  { "--cut_ga",     eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  THRESHOPTS,      "use profile's GA gathering cutoffs to set all thresholding",    6 },
  { "--cut_nc",     eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  THRESHOPTS,      "use profile's NC noise cutoffs to set all thresholding",        6 },
  { "--cut_tc",     eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  THRESHOPTS,      "use profile's TC trusted cutoffs to set all thresholding",      6 },
  /* Control of acceleration pipeline */
  { "--max",        eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, "--F1,--F2,--F3", "Turn all heuristic filters off (less speed, more power)",       7 },
  { "--F1",         eslARG_REAL,  "0.02", NULL, NULL,    NULL,  NULL, "--max",          "MSV threshold: promote hits w/ P <= F1",                        7 },
  { "--F2",         eslARG_REAL,  "1e-3", NULL, NULL,    NULL,  NULL, "--max",          "Vit threshold: promote hits w/ P <= F2",                        7 },
  { "--F3",         eslARG_REAL,  "1e-5", NULL, NULL,    NULL,  NULL, "--max",          "Fwd threshold: promote hits w/ P <= F3",                        7 },
  { "--nobias",     eslARG_NONE,    NULL, NULL, NULL,    NULL,  NULL, "--max",          "turn off composition bias filter",                              7 },
  /* Other options */
  { "--nonull2",    eslARG_NONE,    NULL, NULL, NULL,    NULL,  NULL,  NULL,            "turn off biased composition score corrections",                12 },
  { "-Z",           eslARG_REAL,   FALSE, NULL, "x>0",   NULL,  NULL,  NULL,            "set # of comparisons done, for E-value calculation",           12 },
  { "--domZ",       eslARG_REAL,   FALSE, NULL, "x>0",   NULL,  NULL,  NULL,            "set # of significant seqs, for domain E-value calculation",    12 },
  { "--seed",       eslARG_INT,    "42",  NULL, "n>=0",  NULL,  NULL,  NULL,            "set RNG seed to <n> (if 0: one-time arbitrary seed)",          12 },
  { "--qformat",    eslARG_STRING,  NULL, NULL, NULL,    NULL,  NULL,  NULL,            "assert input <seqfile> is in format <s>: no autodetection",    12 },
#ifdef HMMER_THREADS
  { "--cpu",        eslARG_INT,"0","HMMER_NCPU","n>=0",NULL,  NULL, CPUOPTS,            "number of parallel CPU workers to use for multithreads",       12 },  // multithread parallelization off by default. hmmscan is i/o bound on almost all systems.
#endif
#ifdef HMMER_MPI
  { "--stall",      eslARG_NONE,   FALSE, NULL, NULL,    NULL,"--mpi", NULL,            "arrest after start: for debugging MPI under gdb",              12 },  
  { "--mpi",        eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  MPIOPTS,         "run as an MPI parallel program",                               12 },
#endif
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};


/* struct cfg_s : "Global" application configuration shared by all threads/processes
 * 
 * This structure is passed to routines within main.c, as a means of semi-encapsulation
 * of shared data amongst different parallel processes (threads or MPI processes).
 */
struct cfg_s {
  char            *seqfile;           /* query sequence file                             */
  char            *hmmfile;           /* database HMM file                               */

  int              do_mpi;            /* TRUE if we're doing MPI parallelization         */
  int              nproc;             /* how many MPI processes, total                   */
  int              my_rank;           /* who am I, in 0..nproc-1                         */
};

static char usage[]  = "[-options] <hmmdb> <seqfile>";
static char banner[] = "search sequence(s) against a profile database";

static int  serial_master(ESL_GETOPTS *go, struct cfg_s *cfg);
static int  serial_loop  (WORKER_INFO *info, P7_HMMFILE *hfp);
typedef struct scan_dev_s SCAN_DEV;
struct scan_pool_s;                              /* the cached path's persistent workers */
static int  scan_nactive(uint32_t nmodel, int nworker);
static void scan_pool_destroy(struct scan_pool_s *p);
static int  cached_loop  (WORKER_INFO *info, int infocnt, const P7_HMMCACHE *cache, const BIAS_CACHE *bias,
                          SCAN_DEV *dev, int blocksize, struct scan_pool_s **ret_pool);
static BIAS_CACHE *bias_cache_Create (const P7_HMMCACHE *cache, ESL_ALPHABET *abc);
static void        bias_cache_Destroy(BIAS_CACHE *bc);
static int  scan_cache_minq(void);
static int64_t scan_cache_maxbytes(void);
static int  scan_count_queries(ESL_SQFILE *sqfp, ESL_ALPHABET *abc, int upto);
static int64_t scan_cache_estimate(const char *hmmfile);
static int  scan_qbatch(void);

#ifdef HMMER_THREADS
#define BLOCK_SIZE 1000

static int  thread_loop(ESL_THREADS *obj, ESL_WORK_QUEUE *queue, P7_HMMFILE *hfp);
static void pipeline_thread(void *arg);
#endif

#ifdef HMMER_MPI
static int  mpi_master   (ESL_GETOPTS *go, struct cfg_s *cfg);
static int  mpi_worker   (ESL_GETOPTS *go, struct cfg_s *cfg);
#endif

/* process_commandline()
 * 
 * Processes the commandline, filling in fields in <cfg> and creating and returning
 * an <ESL_GETOPTS> options structure. The help page (hmmsearch -h) is formatted
 * here.
 */
static int
process_commandline(int argc, char **argv, ESL_GETOPTS **ret_go, char **ret_hmmfile, char **ret_seqfile)
{
  ESL_GETOPTS *go = esl_getopts_Create(options);
  int          status;

  if (esl_opt_ProcessEnvironment(go)         != eslOK)  { if (printf("Failed to process environment: %s\n", go->errbuf) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK)  { if (printf("Failed to parse command line: %s\n",  go->errbuf) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if (esl_opt_VerifyConfig(go)               != eslOK)  { if (printf("Failed to parse command line: %s\n",  go->errbuf) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
 
  /* help format: */
  if (esl_opt_GetBoolean(go, "-h") == TRUE) 
    {
      p7_banner_help(stdout, argv[0], banner);
      esl_usage(stdout, argv[0], usage);
      if (puts("\nBasic options:")                                           < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1= group; 2 = indentation; 80=textwidth*/

      if (puts("\nOptions controlling output:")                              < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); 

      if (puts("\nOptions controlling reporting thresholds:")                < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80); 

      if (puts("\nOptions controlling inclusion (significance) thresholds:") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 5, 2, 80); 

      if (puts("\nOptions for model-specific thresholding:")                 < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 6, 2, 80); 

      if (puts("\nOptions controlling acceleration heuristics:")             < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 7, 2, 80); 

      if (puts("\nOther expert options:")                                    < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 12, 2, 80); 
      exit(0);
    }

  if (esl_opt_ArgNumber(go)                 != 2)      { if (puts("Incorrect number of command line arguments.")      < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if ((*ret_hmmfile = esl_opt_GetArg(go, 1)) == NULL)  { if (puts("Failed to get <hmmdb> argument on command line")   < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if ((*ret_seqfile = esl_opt_GetArg(go, 2)) == NULL)  { if (puts("Failed to get <seqfile> argument on command line") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }

  /* Validate any attempted use of stdin streams */
  if (strcmp(*ret_hmmfile, "-") == 0) 
    { if (puts("hmmscan cannot read <hmm database> from stdin stream, because it must have hmmpress'ed auxfiles") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");   goto FAILURE;  }

  *ret_go = go;
  return eslOK;
  
 FAILURE:  /* all errors handled here are user errors, so be polite.  */
  esl_usage(stdout, argv[0], usage);
  if (puts("\nwhere most common options are:")                                 < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
  esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1= group; 2 = indentation; 80=textwidth*/
  if (printf("\nTo see more help on available options, do %s -h\n\n", argv[0]) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
  esl_getopts_Destroy(go);
  exit(1);  

 ERROR:
  if (go) esl_getopts_Destroy(go);
  exit(status);
}


static int
output_header(FILE *ofp, ESL_GETOPTS *go, char *hmmfile, char *seqfile)
{
  p7_banner(ofp, go->argv[0], banner);
  
  if (fprintf(ofp, "# query sequence file:             %s\n", seqfile)                                                                              < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (fprintf(ofp, "# target HMM database:             %s\n", hmmfile)                                                                              < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-o")          && fprintf(ofp, "# output directed to file:         %s\n",            esl_opt_GetString(go, "-o"))          < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--tblout")    && fprintf(ofp, "# per-seq hits tabular output:     %s\n",            esl_opt_GetString(go, "--tblout"))    < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--domtblout") && fprintf(ofp, "# per-dom hits tabular output:     %s\n",            esl_opt_GetString(go, "--domtblout")) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--pfamtblout")&& fprintf(ofp, "# pfam-style tabular hit output:   %s\n",            esl_opt_GetString(go, "--pfamtblout")) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--acc")       && fprintf(ofp, "# prefer accessions over names:    yes\n")                                                 < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--noali")     && fprintf(ofp, "# show alignments in output:       no\n")                                                  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--notextw")   && fprintf(ofp, "# max ASCII text line length:      unlimited\n")                                           < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--textw")     && fprintf(ofp, "# max ASCII text line length:      %d\n",            esl_opt_GetInteger(go, "--textw"))    < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");  
  if (esl_opt_IsUsed(go, "-E")          && fprintf(ofp, "# profile reporting threshold:     E-value <= %g\n", esl_opt_GetReal(go, "-E"))            < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-T")          && fprintf(ofp, "# profile reporting threshold:     score >= %g\n",   esl_opt_GetReal(go, "-T"))            < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--domE")      && fprintf(ofp, "# domain reporting threshold:      E-value <= %g\n", esl_opt_GetReal(go, "--domE"))        < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--domT")      && fprintf(ofp, "# domain reporting threshold:      score >= %g\n",   esl_opt_GetReal(go, "--domT"))        < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--incE")      && fprintf(ofp, "# profile inclusion threshold:     E-value <= %g\n", esl_opt_GetReal(go, "--incE"))        < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--incT")      && fprintf(ofp, "# profile inclusion threshold:     score >= %g\n",   esl_opt_GetReal(go, "--incT"))        < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--incdomE")   && fprintf(ofp, "# domain inclusion threshold:      E-value <= %g\n", esl_opt_GetReal(go, "--incdomE"))     < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--incdomT")   && fprintf(ofp, "# domain inclusion threshold:      score >= %g\n",   esl_opt_GetReal(go, "--incdomT"))     < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--cut_ga")    && fprintf(ofp, "# model-specific thresholding:     GA cutoffs\n")                                          < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--cut_nc")    && fprintf(ofp, "# model-specific thresholding:     NC cutoffs\n")                                          < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--cut_tc")    && fprintf(ofp, "# model-specific thresholding:     TC cutoffs\n")                                          < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--max")       && fprintf(ofp, "# Max sensitivity mode:            on [all heuristic filters off]\n")                      < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--F1")        && fprintf(ofp, "# MSV filter P threshold:       <= %g\n",            esl_opt_GetReal(go, "--F1"))          < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--F2")        && fprintf(ofp, "# Vit filter P threshold:       <= %g\n",            esl_opt_GetReal(go, "--F2"))          < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--F3")        && fprintf(ofp, "# Fwd filter P threshold:       <= %g\n",            esl_opt_GetReal(go, "--F3"))          < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--nobias")    && fprintf(ofp, "# biased composition HMM filter:   off\n")                                                 < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--nonull2")   && fprintf(ofp, "# null2 bias corrections:          off\n")                                                 < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-Z")          && fprintf(ofp, "# sequence search space set to:    %.0f\n",          esl_opt_GetReal(go, "-Z"))            < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--domZ")      && fprintf(ofp, "# domain search space set to:      %.0f\n",          esl_opt_GetReal(go, "--domZ"))        < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--seed"))  {
    if (esl_opt_GetInteger(go, "--seed")==0 && fprintf(ofp, "# random number seed:              one-time arbitrary\n")                              < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
    else if (                                  fprintf(ofp, "# random number seed set to:       %d\n",        esl_opt_GetInteger(go, "--seed"))     < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  }
  if (esl_opt_IsUsed(go, "--qformat")   && fprintf(ofp, "# input seqfile format asserted:   %s\n",            esl_opt_GetString(go, "--qformat"))   < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
                                           
#ifdef HMMER_THREADS
  if (esl_opt_IsUsed(go, "--cpu")) {
    if (esl_opt_GetInteger(go, "--cpu") == 0) { if (fprintf(ofp, "# multithread parallelization:     off\n")                                         < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); }
    else                                      { if (fprintf(ofp, "# multithread parallelization:     %d workers\n", esl_opt_GetInteger(go, "--cpu")) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); }
  }
#endif
#ifdef HMMER_MPI
  if (esl_opt_IsUsed(go, "--mpi")       && fprintf(ofp, "# MPI:                             on\n")                                                  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
#endif
  if (fprintf(ofp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n")                                                 < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  return eslOK;
}


int
main(int argc, char **argv)
{
  ESL_GETOPTS     *go  = NULL;	
  struct cfg_s     cfg;         
  int              status   = eslOK;

  impl_Init();			/* processor-specific initialization */
  p7_FLogsumInit();		/* we're going to use table-driven Logsum() approximations at times */

  /* Initialize what we can in the config structure (without knowing the alphabet yet) */
  cfg.hmmfile    = NULL;
  cfg.seqfile    = NULL;
  cfg.do_mpi     = FALSE;	           /* this gets reset below, if we init MPI */
  cfg.nproc      = 0;		           /* this gets reset below, if we init MPI */
  cfg.my_rank    = 0;		           /* this gets reset below, if we init MPI */

  process_commandline(argc, argv, &go, &cfg.hmmfile, &cfg.seqfile);    

  /* Figure out who we are, and send control there: 
   * we might be an MPI master, an MPI worker, or a serial program.
   */
#ifdef HMMER_MPI
  /* pause the execution of the programs execution until the user has a
   * chance to attach with a debugger and send a signal to resume execution
   * i.e. (gdb) signal SIGCONT
   */
  if (esl_opt_GetBoolean(go, "--stall")) pause();

  if (esl_opt_GetBoolean(go, "--mpi")) 
    {
      cfg.do_mpi     = TRUE;
      MPI_Init(&argc, &argv);
      MPI_Comm_rank(MPI_COMM_WORLD, &(cfg.my_rank));
      MPI_Comm_size(MPI_COMM_WORLD, &(cfg.nproc));

      if (cfg.my_rank > 0)  status = mpi_worker(go, &cfg);
      else 		    status = mpi_master(go, &cfg);

      MPI_Finalize();
    }
  else
#endif /*HMMER_MPI*/
    {
      status = serial_master(go, &cfg);
    }

  esl_getopts_Destroy(go);
  return status;
}


/* serial_master()
 * The serial version of hmmsearch.
 * For each query HMM in <hmmdb> search the database for hits.
 * 
 * A master can only return if it's successful. All errors are handled
 * immediately and fatally with p7_Fail().
 */
static int
serial_master(ESL_GETOPTS *go, struct cfg_s *cfg)
{
  FILE            *ofp      = stdout;	         /* output file for results (default stdout)        */
  FILE            *tblfp    = NULL;		 /* output stream for tabular per-seq (--tblout)    */
  FILE            *domtblfp = NULL;	  	 /* output stream for tabular per-seq (--domtblout) */
  FILE            *pfamtblfp= NULL;              /* output stream for pfam tabular output (--pfamtblout)    */
  int              seqfmt   = eslSQFILE_UNKNOWN; /* format of seqfile                               */
  ESL_SQFILE      *sqfp     = NULL;              /* open seqfile                                    */
  P7_HMMFILE      *hfp      = NULL;		 /* open HMM database file                          */
  ESL_ALPHABET    *abc      = NULL;              /* sequence alphabet                               */
  P7_OPROFILE     *om       = NULL;		 /* target profile                                  */
  ESL_STOPWATCH   *w        = NULL;              /* timing                                          */
  ESL_SQ          *qsq      = NULL;		 /* current query (borrowed)                        */
  ESL_SQ          *qsq_own  = NULL;		 /* the streaming path's own buffer                 */
  int              nquery   = 0;
  int              textw;
  int              status   = eslOK;
  int              hstatus  = eslOK;
  int              sstatus  = eslOK;
  int              i;

  int              ncpus    = 0;

  int              infocnt  = 0;
  WORKER_INFO     *info     = NULL;
  int              infoalloc;                    /* info[] entries allocated (>= infocnt) */
  P7_HMMCACHE     *hcache   = NULL;              /* resident profile database, or NULL */
  struct scan_pool_s *scanpool = NULL;           /* persistent workers for the cached path */
  BIAS_CACHE      *bcache   = NULL;              /* per-model bias filter odds ratios  */
#ifdef HAVE_PEZY
  P7_PZ           *pz       = NULL;              /* PEZY-SC context, or NULL           */
  P7_PZMM         *pzmm     = NULL;              /* all models' score tables on device */
  uint8_t         *pz_xE    = NULL;              /* [nmodel * qb_n] SSV results        */
  char             pzerr[eslERRBUFSIZE];
#endif
  ESL_SQ         **qbuf     = NULL;              /* buffered query batch               */
  int              qb_n     = 0, qb_i = 0, qb_alloc = 0;
  char             cerrbuf[eslERRBUFSIZE];
#ifdef HMMER_THREADS
  P7_OM_BLOCK     *block    = NULL;
  ESL_THREADS     *threadObj= NULL;
  ESL_WORK_QUEUE  *queue    = NULL;
#endif
  char             errbuf[eslERRBUFSIZE];

  w = esl_stopwatch_Create();

  if (esl_opt_GetBoolean(go, "--notextw")) textw = 0;
  else                                     textw = esl_opt_GetInteger(go, "--textw");

  /* If caller declared an input format, decode it */
  if (esl_opt_IsOn(go, "--qformat")) {
    seqfmt = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--qformat"));
    if (seqfmt == eslSQFILE_UNKNOWN) p7_Fail("%s is not a recognized input sequence file format\n", esl_opt_GetString(go, "--qformat"));
  }

  /* Open the target profile database to get the sequence alphabet */
  status = p7_hmmfile_Open(cfg->hmmfile, p7_HMMDBENV, &hfp, errbuf);
  if      (status == eslENOTFOUND) p7_Fail("File existence/permissions problem in trying to open HMM file %s.\n%s\n", cfg->hmmfile, errbuf);
  else if (status == eslEFORMAT)   p7_Fail("File format problem, trying to open HMM file %s.\n%s\n",                  cfg->hmmfile, errbuf);
  else if (status != eslOK)        p7_Fail("Unexpected error %d in opening HMM file %s.\n%s\n",               status, cfg->hmmfile, errbuf);  
  if (! hfp->is_pressed)           p7_Fail("Failed to open binary auxfiles for %s: use hmmpress first\n",             hfp->fname);

  hstatus = p7_oprofile_ReadMSV(hfp, &abc, &om);
  if      (hstatus == eslEFORMAT)   p7_Fail("bad format, binary auxfiles, %s:\n%s",     cfg->hmmfile, hfp->errbuf);
  else if (hstatus == eslEINCOMPAT) p7_Fail("HMM file %s contains different alphabets", cfg->hmmfile);
  else if (hstatus != eslOK)        p7_Fail("Unexpected error in reading HMMs from %s", cfg->hmmfile); 

  p7_oprofile_Destroy(om);
  p7_hmmfile_Close(hfp);
  hfp = NULL;             /* this probe leaves <hfp> dangling otherwise, and the
                           * resident-database path passes it to the pipeline */

  /* Open the query sequence database */
  status = esl_sqfile_OpenDigital(abc, cfg->seqfile, seqfmt, NULL, &sqfp);
  if      (status == eslENOTFOUND) p7_Fail("Failed to open sequence file %s for reading\n",      cfg->seqfile);
  else if (status == eslEFORMAT)   p7_Fail("Sequence file %s is empty or misformatted\n",        cfg->seqfile);
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect format of a stdin or .gz seqfile");
  else if (status != eslOK)        p7_Fail("Unexpected error %d opening sequence file %s\n", status, cfg->seqfile);
  qsq_own = esl_sq_CreateDigital(abc);

  /* Open the results output files */
  if (esl_opt_IsOn(go, "-o"))          { if ((ofp      = fopen(esl_opt_GetString(go, "-o"),          "w")) == NULL)  esl_fatal("Failed to open output file %s for writing\n",                 esl_opt_GetString(go, "-o")); }
  if (esl_opt_IsOn(go, "--tblout"))    { if ((tblfp    = fopen(esl_opt_GetString(go, "--tblout"),    "w")) == NULL)  esl_fatal("Failed to open tabular per-seq output file %s for writing\n", esl_opt_GetString(go, "--tblout")); }
  if (esl_opt_IsOn(go, "--domtblout")) { if ((domtblfp = fopen(esl_opt_GetString(go, "--domtblout"), "w")) == NULL)  esl_fatal("Failed to open tabular per-dom output file %s for writing\n", esl_opt_GetString(go, "--domtblout")); }
  if (esl_opt_IsOn(go, "--pfamtblout")){ if ((pfamtblfp = fopen(esl_opt_GetString(go, "--pfamtblout"), "w")) == NULL)  esl_fatal("Failed to open pfam-style tabular output file %s for writing\n", esl_opt_GetString(go, "--pfamtblout")); }

  output_header(ofp, go, cfg->hmmfile, cfg->seqfile);

#ifdef HMMER_THREADS
  /* initialize thread data */
  ncpus = ESL_MIN( esl_opt_GetInteger(go, "--cpu"), esl_threads_GetCPUCount());
  if (ncpus > 0)
    {
      threadObj = esl_threads_Create(&pipeline_thread);
      queue = esl_workqueue_Create(ncpus * 2);
    }
#endif

  infocnt = (ncpus == 0) ? 1 : ncpus;
  infoalloc = infocnt;                 /* infocnt may shrink below; free this many */
  ESL_ALLOC(info, (ptrdiff_t) sizeof(*info) * infocnt);

  for (i = 0; i < infocnt; ++i)
    {
      info[i].bg    = p7_bg_Create(abc);
#ifdef HMMER_THREADS
      info[i].queue = queue;
#endif
    }

#ifdef HMMER_THREADS
  for (i = 0; i < ncpus * 2; ++i)
    {
      block = p7_oprofile_CreateBlock(BLOCK_SIZE);
      if (block == NULL)    esl_fatal("Failed to allocate sequence block");

      status = esl_workqueue_Init(queue, block);
      if (status != eslOK)  esl_fatal("Failed to add block to work queue");
    }
#endif

  /* With more than a handful of queries it is much cheaper to read the profile
   * database once and scan it from memory than to re-read it per query; see the
   * comment above cached_loop().  Anything unusual -- a query file we cannot
   * rewind to count, a database too large to hold -- leaves hcache NULL and the
   * original per-query loop is used. */
  {
    int     minq = scan_cache_minq();
    int     nq   = (minq > 1) ? scan_count_queries(sqfp, abc, minq) : minq;
    int64_t est  = scan_cache_estimate(cfg->hmmfile);

    if (nq >= minq && est > 0 && est <= scan_cache_maxbytes())
      {
        if (p7_hmmcache_Open(cfg->hmmfile, &hcache, cerrbuf) != eslOK) hcache = NULL;
        if (hcache) bcache = bias_cache_Create(hcache, abc);   /* NULL is tolerated */
      }
  }

  /* On the cached path, use only as many workers as there is work for.  Every
   * query creates and merges one pipeline and one hit list PER WORKER, and a
   * pipeline allocates four DP matrices, so an unused worker is not free even
   * if it never runs -- see scan_nactive() for the measurements.  The reported
   * worker count is --cpu itself, so the output does not change. */
  if (hcache != NULL) infocnt = scan_nactive(hcache->n, infocnt);

#ifdef HAVE_PEZY
  /* With the profile database resident we can also put the SSV filter on the
   * device: build every model's score tables once, keep them there, and score a
   * whole batch of queries against all of them in one launch.  Falling back is
   * always safe -- pzmm NULL just means the host runs the filter. */
  if (hcache != NULL && p7_pz_DeviceCount() > 0)
    {
      if ((pz = p7_pz_Create(0, pzerr)) != NULL)
        {
          int nth = (ncpus > 0 ? ncpus : 1) + 1;
          if (p7_pzmm_Create(pz, hcache->list, (int) hcache->n, nth, &pzmm, pzerr) != eslOK) pzmm = NULL;
          if (pzmm) {
            qb_alloc = scan_qbatch();
            ESL_ALLOC(pz_xE, (size_t) hcache->n * qb_alloc);
          }
          if (pzmm == NULL) { p7_pz_Destroy(pz); pz = NULL; }
        }
    }
#endif
  if (hcache != NULL && qb_alloc == 0) qb_alloc = scan_qbatch();
  if (qb_alloc > 0) ESL_ALLOC(qbuf, sizeof(ESL_SQ *) * qb_alloc);

  /* Outside loop: over each query sequence in <seqfile>.
   *
   * When the profile database is resident, queries are taken a batch at a time
   * so the device can score the whole batch against every model in one launch.
   * Results are still emitted one query at a time, in input order.
   */
  for (;;)
    {
      if (qb_alloc > 0)
        {
          if (qb_i == qb_n)
            {
              int k;
              for (k = 0; k < qb_n; k++) esl_sq_Destroy(qbuf[k]);
              qb_n = qb_i = 0;
              while (qb_n < qb_alloc)
                {
                  ESL_SQ *nx = esl_sq_CreateDigital(abc);
                  if ((sstatus = esl_sqio_Read(sqfp, nx)) != eslOK) { esl_sq_Destroy(nx); break; }
                  qbuf[qb_n++] = nx;
                }
              if (qb_n == 0) break;
              sstatus = eslOK;
#ifdef HAVE_PEZY
              if (pzmm)
                {
                  ESL_DSQ *dq[1]; int ll[1];
                  ESL_DSQ **dqs = NULL; int *lls = NULL;
                  int k2;
                  (void) dq; (void) ll;
                  ESL_ALLOC(dqs, sizeof(ESL_DSQ *) * qb_n);
                  ESL_ALLOC(lls, sizeof(int)       * qb_n);
                  for (k2 = 0; k2 < qb_n; k2++) { dqs[k2] = qbuf[k2]->dsq; lls[k2] = (int) qbuf[k2]->n; }
                  /* A launch that does not come back clean leaves scores we
                   * must not use: give the device up and let the ordinary CPU
                   * pipeline finish the run rather than report wrong hits. */
                  if (p7_pzmm_SSV(pz, pzmm, (const ESL_DSQ **) dqs, lls, qb_n, pz_xE) != eslOK)
                    {
                      fprintf(stderr, "\n# WARNING: a PEZY-SC SSV launch did not complete "
                                      "(see pzcl.log); falling back to the CPU filter for the rest of this run.\n");
                      p7_pzmm_Destroy(pzmm); pzmm = NULL;
                      p7_pz_Destroy(pz);     pz   = NULL;
                    }
                  free(dqs); free(lls);
                }
#endif
            }
          qsq = qbuf[qb_i++];
        }
      else
        {
          esl_sq_Reuse(qsq_own);
          if ((sstatus = esl_sqio_Read(sqfp, qsq_own)) != eslOK) break;
          qsq = qsq_own;
        }

      nquery++;
      esl_stopwatch_Start(w);	                          

      /* Open the target profile database (not needed when it is resident) */
      if (hcache == NULL)
        {
          status = p7_hmmfile_Open(cfg->hmmfile, p7_HMMDBENV, &hfp, NULL);
          if (status != eslOK)        p7_Fail("Unexpected error %d in opening hmm file %s.\n",           status, cfg->hmmfile);

#ifdef HMMER_THREADS
          /* if we are threaded, create a lock to prevent multiple readers */
          if (ncpus > 0)
            {
              status = p7_hmmfile_CreateLock(hfp);
              if (status != eslOK) p7_Fail("Unexpected error %d creating lock\n", status);
            }
#endif
        }

      if (fprintf(ofp, "Query:       %s  [L=%ld]\n", qsq->name, (long) qsq->n) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
      if (qsq->acc[0]  != 0 && fprintf(ofp, "Accession:   %s\n", qsq->acc)     < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
      if (qsq->desc[0] != 0 && fprintf(ofp, "Description: %s\n", qsq->desc)    < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");

      for (i = 0; i < infocnt; ++i)
	{
	  /* Create processing pipeline and hit list */
	  info[i].th  = p7_tophits_Create(); 
	  info[i].pli = p7_pipeline_Create(go, 100, 100, FALSE, p7_SCAN_MODELS); /* M_hint = 100, L_hint = 100 are just dummies for now */
	  /* for two-stage input the pipeline needs <hfp>; a resident database is
	   * already complete, so leaving this NULL skips p7_oprofile_ReadRest() */
	  info[i].pli->hfp = hfp;

	  p7_pli_NewSeq(info[i].pli, qsq);
	  info[i].qsq = qsq;

#ifdef HMMER_THREADS
	  /* esl_threads_AddThread() starts its pipeline_thread immediately and
	   * that thread blocks on the work queue forever, so the cached path
	   * must not create them. */
	  if (ncpus > 0 && hcache == NULL) esl_threads_AddThread(threadObj, &info[i]);
#endif
	}

      {
        SCAN_DEV devd;
        SCAN_DEV *devp = NULL;
        memset(&devd, 0, sizeof(devd));
#ifdef HAVE_PEZY
        if (pzmm) { devd.pzmm = pzmm; devd.xE = pz_xE; devd.nq = qb_n; devd.qi = qb_i - 1; devp = &devd; }
#endif
        if (hcache != NULL) hstatus = cached_loop(info, infocnt, hcache, bcache, devp, BLOCK_SIZE, &scanpool);
        else hstatus = eslENORESULT;
      }
      if (hcache == NULL)
        {
#ifdef HMMER_THREADS
          if (ncpus > 0) hstatus = thread_loop(threadObj, queue, hfp);
          else	         hstatus = serial_loop(info, hfp);
#else
          hstatus = serial_loop(info, hfp);
#endif
        }
      switch(hstatus)
	{
	case eslEFORMAT:   p7_Fail("bad file format in HMM file %s",             cfg->hmmfile);	  break;
	case eslEINCOMPAT: p7_Fail("HMM file %s contains different alphabets",   cfg->hmmfile);	  break;
	case eslEOF: 	  /* do nothing */                                                 	  break;
	default: 	   p7_Fail("Unexpected error in reading HMMs from %s",   cfg->hmmfile); 
	}

      /* merge the results of the search results */
      for (i = 1; i < infocnt; ++i)
	{
	  p7_tophits_Merge(info[0].th, info[i].th);
	  p7_pipeline_Merge(info[0].pli, info[i].pli);

	  p7_pipeline_Destroy(info[i].pli);
	  p7_tophits_Destroy(info[i].th);
	}

      /* Print results */
      p7_tophits_SortBySortkey(info->th);
      p7_tophits_Threshold(info->th, info->pli);

      p7_tophits_Targets(ofp, info->th, info->pli, textw); if (fprintf(ofp, "\n\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
      p7_tophits_Domains(ofp, info->th, info->pli, textw); if (fprintf(ofp, "\n\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");

      if (tblfp)     p7_tophits_TabularTargets(tblfp,    qsq->name, qsq->acc, info->th, info->pli, (nquery == 1));
      if (domtblfp)  p7_tophits_TabularDomains(domtblfp, qsq->name, qsq->acc, info->th, info->pli, (nquery == 1));
      if (pfamtblfp) p7_tophits_TabularXfam(pfamtblfp, qsq->name, qsq->acc, info->th, info->pli);

      esl_stopwatch_Stop(w);
      p7_pli_Statistics(ofp, info->pli, w);
      if (fprintf(ofp, "//\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
      fflush(ofp);

      if (hcache == NULL) { p7_hmmfile_Close(hfp); hfp = NULL; }
      p7_pipeline_Destroy(info->pli);
      p7_tophits_Destroy(info->th);
    }
  if      (sstatus == eslEFORMAT) esl_fatal("Parse failed (sequence file %s):\n%s\n",
					    sqfp->filename, esl_sqfile_GetErrorBuf(sqfp));
  else if (sstatus != eslEOF)     esl_fatal("Unexpected error %d reading sequence file %s",
					    sstatus, sqfp->filename);

#ifdef HAVE_PEZY
  if (pz_xE) free(pz_xE);
  if (pzmm)  p7_pzmm_Destroy(pzmm);
  if (pz)    p7_pz_Destroy(pz);
#endif
  if (qbuf) { int k; for (k = 0; k < qb_n; k++) esl_sq_Destroy(qbuf[k]); free(qbuf); }
  if (bcache) bias_cache_Destroy(bcache);
  if (scanpool) scan_pool_destroy(scanpool);
  if (hcache) p7_hmmcache_Close(hcache);

  /* Terminate outputs - any last words?
   */
  if (tblfp)    p7_tophits_TabularTail(tblfp,    "hmmscan", p7_SCAN_MODELS, cfg->seqfile, cfg->hmmfile, go);
  if (domtblfp) p7_tophits_TabularTail(domtblfp, "hmmscan", p7_SCAN_MODELS, cfg->seqfile, cfg->hmmfile, go);
  if (pfamtblfp)p7_tophits_TabularTail(pfamtblfp,"hmmscan", p7_SEARCH_SEQS, cfg->seqfile, cfg->hmmfile, go);
  if (ofp)      { if (fprintf(ofp, "[ok]\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); }

  /* Cleanup - prepare for successful exit
   */
  for (i = 0; i < infoalloc; ++i)
    p7_bg_Destroy(info[i].bg);

#ifdef HMMER_THREADS
  if (ncpus > 0)
    {
      esl_workqueue_Reset(queue);
      while (esl_workqueue_Remove(queue, (void **) &block) == eslOK)
	p7_oprofile_DestroyBlock(block);
      esl_workqueue_Destroy(queue);
      esl_threads_Destroy(threadObj);
    }
#endif

  free(info);

  if (qsq_own) esl_sq_Destroy(qsq_own);
  esl_stopwatch_Destroy(w);
  esl_alphabet_Destroy(abc);
  esl_sqfile_Close(sqfp);

  if (ofp != stdout) fclose(ofp);
  if (tblfp)         fclose(tblfp);
  if (domtblfp)      fclose(domtblfp);
  if (pfamtblfp)     fclose(pfamtblfp);
  return eslOK;

 ERROR:
  return status;
}

#ifdef HMMER_MPI

/* Define common tags used by the MPI master/slave processes */
#define HMMER_ERROR_TAG          1
#define HMMER_HMM_TAG            2
#define HMMER_SEQUENCE_TAG       3
#define HMMER_BLOCK_TAG          4
#define HMMER_PIPELINE_TAG       5
#define HMMER_TOPHITS_TAG        6
#define HMMER_HIT_TAG            7
#define HMMER_TERMINATING_TAG    8
#define HMMER_READY_TAG          9

/* mpi_failure()
 * Generate an error message.  If the clients rank is not 0, a
 * message is created with the error message and sent to the
 * master process for handling.
 */
static void
mpi_failure(char *format, ...)
{
  va_list  argp;
  int      status = eslFAIL;
  int      len;
  int      rank;
  char     str[512];

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  /* format the error mesg */
  va_start(argp, format);
  len = vsnprintf(str, sizeof(str), format, argp);
  va_end(argp);

  /* make sure the error string is terminated */
  str[sizeof(str)-1] = '\0';

  /* if the caller is the master, print the results and abort */
  if (rank == 0)
    {
      if (fprintf(stderr, "\nError: ") < 0) exit(eslEWRITE);
      if (fprintf(stderr, "%s", str)   < 0) exit(eslEWRITE);
      if (fprintf(stderr, "\n")        < 0) exit(eslEWRITE);
      fflush(stderr);

      MPI_Abort(MPI_COMM_WORLD, status);
      exit(1);
    }
  else
    {
      MPI_Send(str, len, MPI_CHAR, 0, HMMER_ERROR_TAG, MPI_COMM_WORLD);
      pause();
    }
}

#define MAX_BLOCK_SIZE (512*1024)

typedef struct {
  uint64_t  offset;
  uint64_t  length;
  uint64_t  count;
} MSV_BLOCK;

typedef struct {
  int        complete;
  int        size;
  int        current;
  int        last;
  MSV_BLOCK *blocks;
} BLOCK_LIST;

/* this routine parses the database keeping track of the blocks
 * offset within the file, number of sequences and the length
 * of the block.  These blocks are passed as work units to the
 * MPI workers.  If multiple hmm's are in the query file, the
 * blocks are reused without parsing the database a second time.
 */
int next_block(P7_HMMFILE *hfp, BLOCK_LIST *list, MSV_BLOCK *block)
{
  P7_OPROFILE   *om       = NULL;
  ESL_ALPHABET  *abc      = NULL;
  int            status   = eslOK;

  /* if the list has been calculated, use it instead of parsing the database */
  if (list->complete)
    {
      if (list->current == list->last)
	{
	  block->offset = 0;
	  block->length = 0;
	  block->count  = 0;

	  status = eslEOF;
	}
      else
	{
	  int inx = list->current++;

	  block->offset = list->blocks[inx].offset;
	  block->length = list->blocks[inx].length;
	  block->count  = list->blocks[inx].count;

	  status = eslOK;
	}

      return status;
    }

  block->offset = 0;
  block->length = 0;
  block->count = 0;

  while (block->length < MAX_BLOCK_SIZE && (status = p7_oprofile_ReadInfoMSV(hfp, &abc, &om)) == eslOK)
    {
      if (block->count == 0) block->offset = om->roff;
      block->length = om->eoff - block->offset + 1;
      block->count++;
      p7_oprofile_Destroy(om);
    }

  if (status == eslEOF && block->count > 0) status = eslOK;
  if (status == eslEOF) list->complete = 1;

  /* add the block to the list of known blocks */
  if (status == eslOK)
    {
      int inx;

      if (list->last >= list->size)
	{
	  void *tmp;
	  list->size += 500;
	  ESL_RALLOC(list->blocks, tmp, sizeof(MSV_BLOCK) * list->size);
	}

      inx = list->last++;
      list->blocks[inx].offset = block->offset;
      list->blocks[inx].length = block->length;
      list->blocks[inx].count  = block->count;
    }

  return status;

 ERROR:
  return eslEMEM;
}

/* mpi_master()
 * The MPI version of hmmscan
 * Follows standard pattern for a master/worker load-balanced MPI program (J1/78-79).
 * 
 * A master can only return if it's successful. 
 * Errors in an MPI master come in two classes: recoverable and nonrecoverable.
 * 
 * Recoverable errors include all worker-side errors, and any
 * master-side error that do not affect MPI communication. Error
 * messages from recoverable messages are delayed until we've cleanly
 * shut down the workers.
 * 
 * Unrecoverable errors are master-side errors that may affect MPI
 * communication, meaning we cannot count on being able to reach the
 * workers and shut them down. Unrecoverable errors result in immediate
 * p7_Fail()'s, which will cause MPI to shut down the worker processes
 * uncleanly.
 */
static int
mpi_master(ESL_GETOPTS *go, struct cfg_s *cfg)
{
  FILE            *ofp      = stdout;	         /* output file for results (default stdout)        */
  FILE            *tblfp    = NULL;		 /* output stream for tabular per-seq (--tblout)    */
  FILE            *domtblfp = NULL;	  	 /* output stream for tabular per-seq (--domtblout) */
  FILE            *pfamtblfp= NULL;              /* output stream for pfam-style tabular output  (--pfamtblout) */
  int              seqfmt   = eslSQFILE_UNKNOWN; /* format of seqfile                               */
  P7_BG           *bg       = NULL;	         /* null model                                      */
  ESL_SQFILE      *sqfp     = NULL;              /* open seqfile                                    */
  P7_HMMFILE      *hfp      = NULL;		 /* open HMM database file                          */
  ESL_ALPHABET    *abc      = NULL;              /* sequence alphabet                               */
  P7_OPROFILE     *om       = NULL;		 /* target profile                                  */
  ESL_STOPWATCH   *w        = NULL;              /* timing                                          */
  ESL_SQ          *qsq      = NULL;		 /* current query (borrowed)                        */
  ESL_SQ          *qsq_own  = NULL;		 /* the streaming path's own buffer                 */
  int              nquery   = 0;
  int              textw;
  int              status   = eslOK;
  int              hstatus  = eslOK;
  int              sstatus  = eslOK;
  int              dest;

  char            *mpi_buf  = NULL;              /* buffer used to pack/unpack structures */
  int              mpi_size = 0;                 /* size of the allocated buffer */
  BLOCK_LIST      *list     = NULL;
  MSV_BLOCK        block;

  int              i;
  int              size;
  MPI_Status       mpistatus;
  char             errbuf[eslERRBUFSIZE];

  w = esl_stopwatch_Create();

  if (esl_opt_GetBoolean(go, "--notextw")) textw = 0;
  else                                     textw = esl_opt_GetInteger(go, "--textw");

  /* If caller declared an input format, decode it */
  if (esl_opt_IsOn(go, "--qformat")) {
    seqfmt = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--qformat"));
    if (seqfmt == eslSQFILE_UNKNOWN) mpi_failure("%s is not a recognized input sequence file format\n", esl_opt_GetString(go, "--qformat"));
  }

  /* Open the target profile database to get the sequence alphabet */
  status = p7_hmmfile_Open(cfg->hmmfile, p7_HMMDBENV, &hfp, errbuf);
  if      (status == eslENOTFOUND) mpi_failure("File existence/permissions problem in trying to open HMM file %s.\n%s\n", cfg->hmmfile, errbuf);
  else if (status == eslEFORMAT)   mpi_failure("File format problem, trying to open HMM file %s.\n%s\n",                  cfg->hmmfile, errbuf);
  else if (status != eslOK)        mpi_failure("Unexpected error %d in opening HMM file %s.\n%s\n",               status, cfg->hmmfile, errbuf);  
  if (! hfp->is_pressed)           mpi_failure("Failed to open binary auxfiles for %s: use hmmpress first\n",             hfp->fname);
  
  hstatus = p7_oprofile_ReadMSV(hfp, &abc, &om);
  if      (hstatus == eslEFORMAT)   mpi_failure("bad format, binary auxfiles, %s",            cfg->hmmfile);
  else if (hstatus == eslEINCOMPAT) mpi_failure("HMM file %s contains different alphabets",   cfg->hmmfile);
  else if (hstatus != eslOK)        mpi_failure("Unexpected error in reading HMMs from %s",   cfg->hmmfile); 

  p7_oprofile_Destroy(om);
  p7_hmmfile_Close(hfp);

  /* Open the query sequence database */
  status = esl_sqfile_OpenDigital(abc, cfg->seqfile, seqfmt, NULL, &sqfp);
  if      (status == eslENOTFOUND) mpi_failure("Failed to open sequence file %s for reading\n",      cfg->seqfile);
  else if (status == eslEFORMAT)   mpi_failure("Sequence file %s is empty or misformatted\n",        cfg->seqfile);
  else if (status == eslEINVAL)    mpi_failure("Can't autodetect format of a stdin or .gz seqfile");
  else if (status != eslOK)        mpi_failure("Unexpected error %d opening sequence file %s\n", status, cfg->seqfile);

  /* Open the results output files */
  if (esl_opt_IsOn(go, "-o")          && (ofp      = fopen(esl_opt_GetString(go, "-o"),          "w")) == NULL)
    mpi_failure("Failed to open output file %s for writing\n",                 esl_opt_GetString(go, "-o"));
  if (esl_opt_IsOn(go, "--tblout")    && (tblfp    = fopen(esl_opt_GetString(go, "--tblout"),    "w")) == NULL)
    mpi_failure("Failed to open tabular per-seq output file %s for writing\n", esl_opt_GetString(go, "--tblfp"));
  if (esl_opt_IsOn(go, "--domtblout") && (domtblfp = fopen(esl_opt_GetString(go, "--domtblout"), "w")) == NULL)  
    mpi_failure("Failed to open tabular per-dom output file %s for writing\n", esl_opt_GetString(go, "--domtblfp"));
  if (esl_opt_IsOn(go, "--pfamtblout") && (pfamtblfp = fopen(esl_opt_GetString(go, "--pfamtblout"), "w")) == NULL)
    mpi_failure("Failed to open pfam-style tabular output file %s for writing\n", esl_opt_GetString(go, "--pfamtblout"));
 
  ESL_ALLOC(list, sizeof(MSV_BLOCK));
  list->complete = 0;
  list->size     = 0;
  list->current  = 0;
  list->last     = 0;
  list->blocks   = NULL;

  output_header(ofp, go, cfg->hmmfile, cfg->seqfile);
  qsq = esl_sq_CreateDigital(abc);
  bg = p7_bg_Create(abc);

  /* Outside loop: over each query sequence in <seqfile>. */
  while ((sstatus = esl_sqio_Read(sqfp, qsq)) == eslOK)
    {
      P7_PIPELINE     *pli     = NULL;		/* processing pipeline                      */
      P7_TOPHITS      *th      = NULL;        	/* top-scoring sequence hits                */

      nquery++;

      esl_stopwatch_Start(w);	                          

      /* seqfile may need to be rewound (multiquery mode) */
      if (nquery > 1) list->current = 0;

      /* Open the target profile database */
      status = p7_hmmfile_Open(cfg->hmmfile, p7_HMMDBENV, &hfp, NULL);
      if (status != eslOK) mpi_failure("Unexpected error %d in opening hmm file %s.\n", status, cfg->hmmfile);  
  
      if (fprintf(ofp, "Query:       %s  [L=%ld]\n", qsq->name, (long) qsq->n) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
      if (qsq->acc[0]  != 0 && fprintf(ofp, "Accession:   %s\n", qsq->acc)     < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
      if (qsq->desc[0] != 0 && fprintf(ofp, "Description: %s\n", qsq->desc)    < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");

      /* Create processing pipeline and hit list */
      th  = p7_tophits_Create(); 
      pli = p7_pipeline_Create(go, 100, 100, FALSE, p7_SCAN_MODELS); /* M_hint = 100, L_hint = 100 are just dummies for now */
      pli->hfp = hfp;  /* for two-stage input, pipeline needs <hfp> */

      p7_pli_NewSeq(pli, qsq);

      /* Main loop: */
      while ((hstatus = next_block(hfp, list, &block)) == eslOK)
	{
	  if (MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &mpistatus) != 0) 
	    mpi_failure("MPI error %d receiving message from %d\n", mpistatus.MPI_SOURCE);

	  MPI_Get_count(&mpistatus, MPI_PACKED, &size);
	  if (mpi_buf == NULL || size > mpi_size) {
	    void *tmp;
	    ESL_RALLOC(mpi_buf, tmp, sizeof(char) * size);
	    mpi_size = size; 
	  }

	  dest = mpistatus.MPI_SOURCE;
	  MPI_Recv(mpi_buf, size, MPI_PACKED, dest, mpistatus.MPI_TAG, MPI_COMM_WORLD, &mpistatus);

	  if (mpistatus.MPI_TAG == HMMER_ERROR_TAG)
	    mpi_failure("MPI client %d raised error:\n%s\n", dest, mpi_buf);
	  if (mpistatus.MPI_TAG != HMMER_READY_TAG)
	    mpi_failure("Unexpected tag %d from %d\n", mpistatus.MPI_TAG, dest);
      
	  MPI_Send(&block, 3, MPI_LONG_LONG_INT, dest, HMMER_BLOCK_TAG, MPI_COMM_WORLD);
	}
      switch(hstatus)
	{
	case eslEFORMAT:
	  mpi_failure("bad file format in HMM file %s",              cfg->hmmfile);
	  break;
	case eslEINCOMPAT:
	  mpi_failure("HMM file %s contains different alphabets",    cfg->hmmfile);
	  break;
	case eslEOF:
	  /* do nothing */
	  break;
	default:
	  mpi_failure("Unexpected error %d in reading HMMs from %s", hstatus, cfg->hmmfile); 
	}

      block.offset = 0;
      block.length = 0;
      block.count  = 0;

      /* wait for all workers to finish up their work blocks */
      for (i = 1; i < cfg->nproc; ++i)
	{
	  if (MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &mpistatus) != 0) 
	    mpi_failure("MPI error %d receiving message from %d\n", mpistatus.MPI_SOURCE);

	  MPI_Get_count(&mpistatus, MPI_PACKED, &size);
	  if (mpi_buf == NULL || size > mpi_size) {
	    void *tmp;
	    ESL_RALLOC(mpi_buf, tmp, sizeof(char) * size);
	    mpi_size = size; 
	  }

	  dest = mpistatus.MPI_SOURCE;
	  MPI_Recv(mpi_buf, size, MPI_PACKED, dest, mpistatus.MPI_TAG, MPI_COMM_WORLD, &mpistatus);

	  if (mpistatus.MPI_TAG == HMMER_ERROR_TAG)
	    mpi_failure("MPI client %d raised error:\n%s\n", dest, mpi_buf);
	  if (mpistatus.MPI_TAG != HMMER_READY_TAG)
	    mpi_failure("Unexpected tag %d from %d\n", mpistatus.MPI_TAG, dest);
	}

      /* merge the results of the search results */
      for (dest = 1; dest < cfg->nproc; ++dest)
	{
	  P7_PIPELINE     *mpi_pli   = NULL;
	  P7_TOPHITS      *mpi_th    = NULL;

	  /* send an empty block to signal the worker they are done */
	  MPI_Send(&block, 3, MPI_LONG_LONG_INT, dest, HMMER_BLOCK_TAG, MPI_COMM_WORLD);

	  /* wait for the results */
	  if ((status = p7_tophits_MPIRecv(dest, HMMER_TOPHITS_TAG, MPI_COMM_WORLD, &mpi_buf, &mpi_size, &mpi_th)) != eslOK)
	    mpi_failure("Unexpected error %d receiving tophits from %d", status, dest);

	  if ((status = p7_pipeline_MPIRecv(dest, HMMER_PIPELINE_TAG, MPI_COMM_WORLD, &mpi_buf, &mpi_size, go, &mpi_pli)) != eslOK)
	    mpi_failure("Unexpected error %d receiving pipeline from %d", status, dest);

	  p7_tophits_Merge(th, mpi_th);
	  p7_pipeline_Merge(pli, mpi_pli);

	  p7_pipeline_Destroy(mpi_pli);
	  p7_tophits_Destroy(mpi_th);
	}

      /* Print the results.  */
      p7_tophits_SortBySortkey(th);
      p7_tophits_Threshold(th, pli);
      p7_tophits_Targets(ofp, th, pli, textw); fprintf(ofp, "\n\n");
      p7_tophits_Domains(ofp, th, pli, textw); fprintf(ofp, "\n\n");

      if (tblfp)     p7_tophits_TabularTargets(tblfp,    qsq->name, qsq->acc, th, pli, (nquery == 1));
      if (domtblfp)  p7_tophits_TabularDomains(domtblfp, qsq->name, qsq->acc, th, pli, (nquery == 1));
      if (pfamtblfp) p7_tophits_TabularXfam(pfamtblfp,   qsq->name, qsq->acc, th, pli);

      esl_stopwatch_Stop(w);
      p7_pli_Statistics(ofp, pli, w);
      if (fprintf(ofp, "//\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");

      p7_hmmfile_Close(hfp);
      p7_pipeline_Destroy(pli);
      p7_tophits_Destroy(th);
      esl_sq_Reuse(qsq);
    }
  if (sstatus == eslEFORMAT) 
    mpi_failure("Parse failed (sequence file %s):\n%s\n", sqfp->filename, esl_sqfile_GetErrorBuf(sqfp));
  else if (sstatus != eslEOF)     
    mpi_failure("Unexpected error %d reading sequence file %s", sstatus, sqfp->filename);

  /* monitor all the workers to make sure they have ended */
  for (i = 1; i < cfg->nproc; ++i)
    {
      if (MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &mpistatus) != 0) 
	mpi_failure("MPI error %d receiving message from %d\n", mpistatus.MPI_SOURCE);

      MPI_Get_count(&mpistatus, MPI_PACKED, &size);
      if (mpi_buf == NULL || size > mpi_size) {
	void *tmp;
	ESL_RALLOC(mpi_buf, tmp, sizeof(char) * size);
	mpi_size = size; 
      }

      dest = mpistatus.MPI_SOURCE;
      MPI_Recv(mpi_buf, size, MPI_PACKED, dest, mpistatus.MPI_TAG, MPI_COMM_WORLD, &mpistatus);

      if (mpistatus.MPI_TAG == HMMER_ERROR_TAG)
	mpi_failure("MPI client %d raised error:\n%s\n", dest, mpi_buf);
      if (mpistatus.MPI_TAG != HMMER_TERMINATING_TAG)
	mpi_failure("Unexpected tag %d from %d\n", mpistatus.MPI_TAG, dest);
    }

 /* Terminate outputs - any last words?
   */
  if (tblfp)    p7_tophits_TabularTail(tblfp,    "hmmscan", p7_SCAN_MODELS, cfg->seqfile, cfg->hmmfile, go);
  if (domtblfp) p7_tophits_TabularTail(domtblfp, "hmmscan", p7_SCAN_MODELS, cfg->seqfile, cfg->hmmfile, go);
  if (pfamtblfp)p7_tophits_TabularTail(pfamtblfp, "hmmscan", p7_SEARCH_SEQS, cfg->seqfile, cfg->hmmfile, go);
  if (ofp)      { if (fprintf(ofp, "[ok]\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); }

  /* Cleanup - prepare for successful exit
   */
  free(list);
  if (mpi_buf != NULL) free(mpi_buf);

  p7_bg_Destroy(bg);
  if (qsq_own) esl_sq_Destroy(qsq_own);
  esl_stopwatch_Destroy(w);
  esl_alphabet_Destroy(abc);
  esl_sqfile_Close(sqfp);

  if (ofp != stdout) fclose(ofp);
  if (tblfp)         fclose(tblfp);
  if (domtblfp)      fclose(domtblfp);
  if (pfamtblfp)     fclose(pfamtblfp);

  return eslOK;

 ERROR:
  return status;
}


static int
mpi_worker(ESL_GETOPTS *go, struct cfg_s *cfg)
{
  int              seqfmt   = eslSQFILE_UNKNOWN; /* format of seqfile                               */
  P7_BG           *bg       = NULL;	         /* null model                                      */
  ESL_SQFILE      *sqfp     = NULL;              /* open seqfile                                    */
  P7_HMMFILE      *hfp      = NULL;		 /* open HMM database file                          */
  ESL_ALPHABET    *abc      = NULL;              /* sequence alphabet                               */
  P7_OPROFILE     *om       = NULL;		 /* target profile                                  */
  ESL_STOPWATCH   *w        = NULL;              /* timing                                          */
  ESL_SQ          *qsq      = NULL;		 /* current query (borrowed)                        */
  ESL_SQ          *qsq_own  = NULL;		 /* the streaming path's own buffer                 */
  int              status   = eslOK;
  int              hstatus  = eslOK;
  int              sstatus  = eslOK;

  char            *mpi_buf  = NULL;              /* buffer used to pack/unpack structures */
  int              mpi_size = 0;                 /* size of the allocated buffer */

  MPI_Status       mpistatus;
  char             errbuf[eslERRBUFSIZE];

  w = esl_stopwatch_Create();

  /* Open the target profile database to get the sequence alphabet */
  status = p7_hmmfile_Open(cfg->hmmfile, p7_HMMDBENV, &hfp, errbuf);
  if      (status == eslENOTFOUND) mpi_failure("File existence/permissions problem in trying to open HMM file %s.\n%s\n", cfg->hmmfile, errbuf);
  else if (status == eslEFORMAT)   mpi_failure("File format problem in trying to open HMM file %s.\n%s\n",                cfg->hmmfile, errbuf);
  else if (status != eslOK)        mpi_failure("Unexpected error %d in opening HMM file %s.\n%s\n",               status, cfg->hmmfile, errbuf);  
  if (! hfp->is_pressed)           mpi_failure("Failed to open binary dbs for HMM file %s: use hmmpress first\n",         hfp->fname);

  hstatus = p7_oprofile_ReadMSV(hfp, &abc, &om);
  if      (hstatus == eslEFORMAT)   mpi_failure("bad file format in HMM file %s",             cfg->hmmfile);
  else if (hstatus == eslEINCOMPAT) mpi_failure("HMM file %s contains different alphabets",   cfg->hmmfile);
  else if (hstatus != eslOK)        mpi_failure("Unexpected error in reading HMMs from %s",   cfg->hmmfile); 

  p7_oprofile_Destroy(om);
  p7_hmmfile_Close(hfp);

  /* Open the query sequence database */
  status = esl_sqfile_OpenDigital(abc, cfg->seqfile, seqfmt, NULL, &sqfp);
  if      (status == eslENOTFOUND) mpi_failure("Failed to open sequence file %s for reading\n",      cfg->seqfile);
  else if (status == eslEFORMAT)   mpi_failure("Sequence file %s is empty or misformatted\n",        cfg->seqfile);
  else if (status == eslEINVAL)    mpi_failure("Can't autodetect format of a stdin or .gz seqfile");
  else if (status != eslOK)        mpi_failure("Unexpected error %d opening sequence file %s\n", status, cfg->seqfile);

  qsq = esl_sq_CreateDigital(abc);
  bg = p7_bg_Create(abc);

  /* Outside loop: over each query sequence in <seqfile>. */
  while ((sstatus = esl_sqio_Read(sqfp, qsq)) == eslOK)
    {
      P7_PIPELINE     *pli     = NULL;		/* processing pipeline                      */
      P7_TOPHITS      *th      = NULL;        	/* top-scoring sequence hits                */

      MSV_BLOCK        block;

      esl_stopwatch_Start(w);

      status = 0;
      MPI_Send(&status, 1, MPI_INT, 0, HMMER_READY_TAG, MPI_COMM_WORLD);

      /* Open the target profile database */
      status = p7_hmmfile_Open(cfg->hmmfile, p7_HMMDBENV, &hfp, NULL);
      if (status != eslOK) mpi_failure("Unexpected error %d in opening hmm file %s.\n", status, cfg->hmmfile);  
  
      /* Create processing pipeline and hit list */
      th  = p7_tophits_Create(); 
      pli = p7_pipeline_Create(go, 100, 100, FALSE, p7_SCAN_MODELS); /* M_hint = 100, L_hint = 100 are just dummies for now */
      pli->hfp = hfp;  /* for two-stage input, pipeline needs <hfp> */

      p7_pli_NewSeq(pli, qsq);

      /* receive a sequence block from the master */
      MPI_Recv(&block, 3, MPI_LONG_LONG_INT, 0, HMMER_BLOCK_TAG, MPI_COMM_WORLD, &mpistatus);
      while (block.count > 0)
	{
	  uint64_t length = 0;
	  uint64_t count  = block.count;

	  hstatus = p7_oprofile_Position(hfp, block.offset);
	  if (hstatus != eslOK) mpi_failure("Cannot position optimized model to %ld\n", block.offset);

	  while (count > 0 && (hstatus = p7_oprofile_ReadMSV(hfp, &abc, &om)) == eslOK)
	    {
	      length = om->eoff - block.offset + 1;

	      p7_pli_NewModel(pli, om, bg);
	      p7_bg_SetLength(bg, qsq->n);
	      p7_oprofile_ReconfigLength(om, qsq->n);
	      
	      p7_Pipeline(pli, om, bg, qsq, NULL, th);
	      
	      p7_oprofile_Destroy(om);
	      p7_pipeline_Reuse(pli);

	      --count;
	    }

	  /* check the status of reading the hmm */

	  /* lets do a little bit of sanity checking here to make sure the blocks are the same */
	  if (count > 0)              
	    {
	      switch(hstatus)
		{
		case eslEFORMAT:
		  mpi_failure("bad file format in HMM file %s",              cfg->hmmfile);
		  break;
		case eslEINCOMPAT:
		  mpi_failure("HMM file %s contains different alphabets",    cfg->hmmfile);
		  break;
		case eslOK:
		case eslEOF:
		  mpi_failure("Block count mismatch - expected %ld found %ld at offset %ld\n", block.count, block.count-count, block.offset);
		  break;
		default:
		  mpi_failure("Unexpected error %d in reading HMMs from %s", hstatus, cfg->hmmfile); 
		}
	    }
	  if (block.length != length) 
	    mpi_failure("Block length mismatch - expected %ld found %ld at offset %ld\n", block.length, length, block.offset);

	  /* inform the master we need another block of sequences */
	  status = 0;
	  MPI_Send(&status, 1, MPI_INT, 0, HMMER_READY_TAG, MPI_COMM_WORLD);

	  /* wait for the next block of sequences */
	  MPI_Recv(&block, 3, MPI_LONG_LONG_INT, 0, HMMER_BLOCK_TAG, MPI_COMM_WORLD, &mpistatus);
	}

      esl_stopwatch_Stop(w);

      /* Send the top hits back to the master. */
      p7_tophits_MPISend(th, 0, HMMER_TOPHITS_TAG, MPI_COMM_WORLD,  &mpi_buf, &mpi_size);
      p7_pipeline_MPISend(pli, 0, HMMER_PIPELINE_TAG, MPI_COMM_WORLD,  &mpi_buf, &mpi_size);

      p7_hmmfile_Close(hfp);
      p7_pipeline_Destroy(pli);
      p7_tophits_Destroy(th);
      esl_sq_Reuse(qsq);
    } /* end outer loop over query HMMs */
  if (sstatus == eslEFORMAT) 
    mpi_failure("Parse failed (sequence file %s):\n%s\n", sqfp->filename, esl_sqfile_GetErrorBuf(sqfp));
  else if (sstatus != eslEOF)
    mpi_failure("Unexpected error %d reading sequence file %s", sstatus, sqfp->filename);

  status = 0;
  MPI_Send(&status, 1, MPI_INT, 0, HMMER_TERMINATING_TAG, MPI_COMM_WORLD);

  if (mpi_buf != NULL) free(mpi_buf);

  p7_bg_Destroy(bg);

  if (qsq_own) esl_sq_Destroy(qsq_own);
  esl_stopwatch_Destroy(w);
  esl_alphabet_Destroy(abc);
  esl_sqfile_Close(sqfp);

  return eslOK;
}
#endif /*HMMER_MPI*/

/*****************************************************************
 * Scanning a resident model database
 *
 * hmmscan re-opens and re-reads the whole pressed profile database inside its
 * per-query loop.  With a multi-sequence query that is by far the dominant
 * cost: on Pfam-A (30,134 models, 531 MB of .h3f) one pass costs ~0.15 s, so a
 * 50-sequence query spends ~7.5 s of its ~9.4 s just re-reading models.
 *
 * Reading the database once into a P7_HMMCACHE and scanning that instead
 * removes all of it.  This is the same thing hmmpgmd's workers do (see
 * hmmdwrkr.c), so the pipeline is already built to score against cached
 * profiles: with pli->hfp left NULL it skips p7_oprofile_ReadRest(), because
 * the cached profiles are already complete.
 *
 * Model blocks are handed to workers round robin by block index.  That matters
 * for more than load: p7_pli_TargetReportable() tests against pli->Z, which
 * mid-scan is the worker's own running model count, so which marginal hits get
 * stored depends on how models were spread over threads.  Using HMMER's own
 * block size and a fixed round robin reproduces what its work queue does when
 * the workers keep up with each other.
 *****************************************************************/

/* The bias filter's two-state HMM depends only on the model (its length and
 * mean composition), yet p7_bg_SetFilter() rebuilds it for every
 * (model, query) pair.  esl_hmm_Configure() -- which derives the emission odds
 * ratios, including averaging over every degenerate residue code -- was 21% of
 * the whole run once the profile database was resident.
 *
 * So derive the odds ratios once per model and let the workers borrow them.
 * Only <eo> is cached: esl_hmm_Configure() writes nothing else, and the rest of
 * the filter HMM is either model-independent or two floats derived from M.
 * Each worker keeps its own t[] and pi[], so nothing shared is ever written.
 */
static void
bias_cache_Destroy(BIAS_CACHE *bc)
{
  if (!bc) return;
  if (bc->rows) free(bc->rows);
  if (bc->data) free(bc->data);
  free(bc);
}

static BIAS_CACHE *
bias_cache_Create(const P7_HMMCACHE *cache, ESL_ALPHABET *abc)
{
  BIAS_CACHE *bc = NULL;
  P7_BG      *bg = NULL;
  uint32_t    m;
  int         x, Kp = abc->Kp, fM;
  int         status;

  if ((bg = p7_bg_Create(abc)) == NULL) return NULL;
  fM = bg->fhmm->M;

  ESL_ALLOC(bc, sizeof(BIAS_CACHE));
  bc->rows = NULL; bc->data = NULL; bc->n = cache->n; bc->Kp = Kp; bc->fM = fM;
  ESL_ALLOC(bc->rows, sizeof(float *) * (size_t) cache->n * Kp);
  ESL_ALLOC(bc->data, sizeof(float)   * (size_t) cache->n * Kp * fM);

  for (m = 0; m < cache->n; m++)
    {
      /* built by the same code that the unaccelerated path uses, so the values
       * are identical by construction */
      p7_bg_SetFilter(bg, cache->list[m]->M, cache->list[m]->compo);
      for (x = 0; x < Kp; x++)
        {
          float *dst = bc->data + ((size_t) m * Kp + x) * fM;
          int    k;
          for (k = 0; k < fM; k++) dst[k] = bg->fhmm->eo[x][k];
          bc->rows[(size_t) m * Kp + x] = dst;
        }
    }
  p7_bg_Destroy(bg);
  return bc;

 ERROR:
  if (bg) p7_bg_Destroy(bg);
  bias_cache_Destroy(bc);
  return NULL;
}

typedef struct {
  WORKER_INFO       *info;
  const P7_HMMCACHE *cache;
  const BIAS_CACHE  *bias;      /* may be NULL */
#ifdef HAVE_PEZY
  P7_PZMM           *pzmm;      /* NULL when the device is not in use */
  const uint8_t     *xE;        /* [nmodel * nq] raw-byte SSV scores   */
  int                nq, qi;    /* batch width, and this query's slot  */
#endif
  int                widx, nworker;
  uint32_t           chunk;     /* models handed out per grab                */
  uint32_t          *next;      /* shared cursor, bumped atomically          */
  int                status;
  struct scan_pool_s *pool;     /* NULL when run directly (infocnt == 1) */
} SCAN_WORKER;

/* How many of the requested workers to actually run.
 *
 * hmmscan runs one query at a time, so a worker's whole share of a query is
 * cache->n / nworker models -- with Pfam-A and --cpu 128 that is 235 models,
 * about 0.5 ms of work, against 0.05-0.1 ms to wake it and put it back to
 * sleep.  MEASURED over 1000 queries against Pfam-A (wall / sys):
 *
 *   --cpu        16        32       128
 *   spawn/query  9.24/2.8  9.08/4.8  12.82/10.9
 *   pool         8.51/2.8  8.42/5.3   9.84/17.1
 *   pool+chunk   7.30/3.2  7.25/8.2  10.07/39.0
 *
 * More threads is worse past about 32, and it is all kernel time.  So cap the
 * number that run by the work available.  The output is invariant to the worker
 * count anyway (measured: --cpu 8 / 32 / 128 differ only in the header line
 * that reports it, and that line prints --cpu itself, not this).
 *
 * The cap is applied to <infocnt> in main() as well, not just to the pool,
 * because every query creates and merges one pipeline and one hit list PER
 * WORKER and a pipeline allocates four DP matrices -- with the cap on the pool
 * alone, --cpu 128 stalled at 8.09 s against 7.25 s for a plain --cpu 32.
 * With both, --cpu 128 lands at 7.23 s.  MEASURED at --cpu 128, 1000 queries:
 *
 *   models/worker   256     512    1024    2048
 *   workers         117      58      29      14
 *   wall           9.87    7.62    7.23    7.57
 *
 * HMMER_HMMCACHE_MINMODELS overrides the divisor.
 */
#define SCAN_MIN_MODELS_PER_WORKER 1024
static int
scan_nactive(uint32_t nmodel, int nworker)
{
  int n;
  const char *v = getenv("HMMER_HMMCACHE_MINMODELS");
  int  minm = SCAN_MIN_MODELS_PER_WORKER;

  if (v) { int k = atoi(v); if (k > 0) minm = k; }
  n = (int) (nmodel / (uint32_t) minm);
  if (n < 1)        n = 1;
  if (n > nworker)  n = nworker;
  return n;
}

/* Models handed to a worker per atomic grab.
 *
 * cached_loop() used to hand out HMMER's 1000-model blocks by a static stride.
 * With Pfam-A (30,134 models) that is only 31 blocks, so with --cpu 128
 * exactly 31 workers got work and 97 woke up, found none, and went back to
 * sleep.  That did not show up while every query was also paying for 128
 * thread spawns; once the pool removed those, it did.
 *
 * Per (model, query) pair the work is ~2 us, so a chunk of a few tens gives
 * plenty of work per atomic while leaving many chunks per worker to average
 * out the ~2% of models that pass F1 and cost orders of magnitude more.
 */
static uint32_t
scan_chunk(uint32_t nmodel, int nworker)
{
  uint32_t n;
  const char *v = getenv("HMMER_HMMCACHE_CHUNK");
  if (v) { int k = atoi(v); if (k > 0) return (uint32_t) k; }

  n = nmodel / ((uint32_t) nworker * 64u);
  if (n <  8u) n =  8u;
  if (n > 64u) n = 64u;
  return n;
}

static void
scan_worker_setup(SCAN_WORKER *sw)
{
  ESL_HMM *fh = sw->info->bg->fhmm;

  /* impl_Init() turns on SSE flush-to-zero, which changes float results for
   * denormals, so every thread that scores must set it. */
  impl_Init();

  if (sw->bias)
    {
      /* the parts of the filter HMM that no model changes.  p7_bg_SetLength()
       * overwrites t[0][0..1] per query, so we need not set them here. */
      fh->t[0][2] = 1.0f;
      fh->t[1][2] = 1.0f;
      fh->pi[0]   = 0.999;
      fh->pi[1]   = 0.001;
    }
}

/* One query's share of the resident database. */
static void
scan_worker_query(SCAN_WORKER *sw)
{
  WORKER_INFO       *info  = sw->info;
  const P7_HMMCACHE *cache = sw->cache;
  const BIAS_CACHE  *bias  = sw->bias;
  ESL_HMM           *fh    = info->bg->fhmm;
  float            **own_eo = fh->eo;

  sw->status = eslOK;

  for (;;)
    {
      uint32_t blk = (uint32_t) __sync_fetch_and_add(sw->next, 1u);
      uint32_t lo  = blk * sw->chunk;
      uint32_t hi, m;

      if (lo >= cache->n) break;
      hi = lo + sw->chunk;  if (hi > cache->n) hi = cache->n;

      for (m = lo; m < hi; m++)
        {
          P7_OPROFILE *om = cache->list[m];

          if (bias)
            {
              float L1 = (float) om->M / 8.0f;     /* mean length of the biased segment */
              fh->t[1][0] = 1.0f / (L1 + 1.0f);
              fh->t[1][1] =   L1 / (L1 + 1.0f);
              fh->eo      = bias->rows + (size_t) m * bias->Kp;   /* borrowed */
              p7_pli_NewModelNoBias(info->pli, om, info->bg);
            }
          else p7_pli_NewModel(info->pli, om, info->bg);

          p7_bg_SetLength(info->bg, info->qsq->n);
          p7_oprofile_ReconfigLength(om, info->qsq->n);

#ifdef HAVE_PEZY
          /* The device already scored every (model, query) pair.  Everything the
           * pipeline does before the MSV filter follows from the query length
           * alone, so we can reproduce its F1 test here and enter the pipeline
           * just past the filter. */
          if (sw->pzmm && p7_pzmm_ModelOnDevice(sw->pzmm, (int) m) && info->qsq->n > 0)
            {
              float  usc, nullsc, seq_score;
              double P;

              p7_omx_GrowTo(info->pli->oxf, om->M, 0, info->qsq->n);
              p7_bg_NullOne(info->bg, info->qsq->dsq, info->qsq->n, &nullsc);

              if (p7_pz_SSVOutcome(sw->xE[(size_t) m * sw->nq + sw->qi], om, &usc) == eslENORESULT)
                p7_MSVFilter(info->qsq->dsq, info->qsq->n, om, info->pli->oxf, &usc);  /* SSV could not decide */

              seq_score = (usc - nullsc) / eslCONST_LOG2;
              P = esl_gumbel_surv(seq_score, om->evparam[p7_MMU], om->evparam[p7_MLAMBDA]);
              if (P <= info->pli->F1 &&
                  (sw->status = p7_Pipeline_MSVDone(info->pli, om, info->bg, info->qsq, NULL, info->th, usc)) != eslOK)
                { fh->eo = own_eo; return; }
            }
          else
#endif
          if ((sw->status = p7_Pipeline(info->pli, om, info->bg, info->qsq, NULL, info->th)) != eslOK)
            { fh->eo = own_eo; return; }

          p7_pipeline_Reuse(info->pli);
          /* the profiles belong to the cache: do not destroy them */
        }
    }
  fh->eo = own_eo;     /* give <bg> its own table back before anyone frees it */
}

/* A persistent pool for the cached (resident-database) path.
 *
 * cached_loop() used to spawn <infocnt> threads and join them for EVERY query,
 * and hmmscan runs one query at a time: with Pfam-A and --cpu 128 that is
 * 128 spawns per query, 563,000 over the E. coli proteome.  MEASURED over 1000
 * queries against Pfam-A: --cpu 128 took wall 12.82 s of which 10.90 s was
 * system time, while --cpu 32 took 9.08 s with 4.78 s -- asking for MORE
 * threads made it slower, and the difference was almost all kernel time.
 *
 * So create the threads once and park them on a condition variable.  Each
 * query bumps a generation counter and broadcasts; the last worker out signals
 * the master.  Everything that varies per query (the pipeline, the hit list,
 * the query sequence, the device's batch slot) is reached through the
 * WORKER_INFO pointer or refreshed in the SCAN_WORKER before the broadcast.
 */
struct scan_pool_s {
  SCAN_WORKER    *sw;
  pthread_t      *tid;
  int             nworker;
  pthread_mutex_t mtx;
  pthread_cond_t  go;         /* master -> workers: a new generation is ready */
  pthread_cond_t  done;       /* last worker -> master                        */
  int             gen;        /* bumped once per query                        */
  int             nrunning;
  int             shutdown;
  uint32_t        next;       /* the models cursor, reset before each query   */
};
typedef struct scan_pool_s SCAN_POOL;

static void *
scan_pool_worker(void *arg)
{
  SCAN_WORKER *sw    = (SCAN_WORKER *) arg;
  SCAN_POOL   *p     = sw->pool;
  int          mygen = 0;

  scan_worker_setup(sw);
  for (;;)
    {
      pthread_mutex_lock(&p->mtx);
      while (p->gen == mygen && ! p->shutdown) pthread_cond_wait(&p->go, &p->mtx);
      if (p->shutdown) { pthread_mutex_unlock(&p->mtx); break; }
      mygen = p->gen;
      pthread_mutex_unlock(&p->mtx);

      scan_worker_query(sw);

      pthread_mutex_lock(&p->mtx);
      if (--p->nrunning == 0) pthread_cond_signal(&p->done);
      pthread_mutex_unlock(&p->mtx);
    }
  return NULL;
}

static void
scan_pool_destroy(SCAN_POOL *p)
{
  int i;
  if (p == NULL) return;
  /* With one worker, cached_loop() runs it inline: there are no threads to
   * join, and the mutex and condition variables were never initialized.  The
   * zeroed mutex happens to behave like PTHREAD_MUTEX_INITIALIZER, so the
   * lock and the broadcast came back fine and it was p->tid[0] -- NULL --
   * that segfaulted.  Skip the whole teardown, not just the join. */
  if (p->tid)
    {
      pthread_mutex_lock(&p->mtx);
      p->shutdown = TRUE;
      pthread_cond_broadcast(&p->go);
      pthread_mutex_unlock(&p->mtx);
      for (i = 0; i < p->nworker; i++) pthread_join(p->tid[i], NULL);
      pthread_cond_destroy(&p->done);
      pthread_cond_destroy(&p->go);
      pthread_mutex_destroy(&p->mtx);
      free(p->tid);
    }
  free(p->sw);
  free(p);
}

static int
cached_loop(WORKER_INFO *info, int infocnt, const P7_HMMCACHE *cache, const BIAS_CACHE *bias,
            SCAN_DEV *dev, int blocksize, SCAN_POOL **ret_pool)
{
  (void) blocksize;            /* HMMER's block size is not what we schedule by */
  SCAN_POOL   *p  = *ret_pool;
  SCAN_WORKER *sw = NULL;
  int          i;
  int          status;

  /* First call: build the pool (or, for one worker, just a SCAN_WORKER). */
  if (p == NULL)
    {
      ESL_ALLOC(p, sizeof(SCAN_POOL));
      memset(p, 0, sizeof(SCAN_POOL));
      p->nworker = scan_nactive(cache->n, infocnt);
      ESL_ALLOC(p->sw, sizeof(SCAN_WORKER) * p->nworker);
      memset(p->sw, 0, sizeof(SCAN_WORKER) * p->nworker);
      for (i = 0; i < p->nworker; i++)
        {
          p->sw[i].info      = &info[i];
          p->sw[i].cache     = cache;
          p->sw[i].bias      = bias;
          p->sw[i].widx      = i;
          p->sw[i].nworker   = p->nworker;
          p->sw[i].chunk     = scan_chunk(cache->n, p->nworker);
          p->sw[i].next      = &p->next;
          p->sw[i].pool      = p;
        }
      if (p->nworker > 1)
        {
          pthread_mutex_init(&p->mtx, NULL);
          pthread_cond_init(&p->go, NULL);
          pthread_cond_init(&p->done, NULL);
          ESL_ALLOC(p->tid, sizeof(pthread_t) * p->nworker);
          for (i = 0; i < p->nworker; i++)
            if (pthread_create(&p->tid[i], NULL, scan_pool_worker, &p->sw[i]) != 0)
              { status = eslESYS; goto ERROR; }
        }
      else scan_worker_setup(&p->sw[0]);
      *ret_pool = p;
    }
  sw = p->sw;

  /* Refresh what changes per query.  info[] itself does not move, and the
   * pipeline / hit list / query sequence are reached through it. */
  for (i = 0; i < p->nworker; i++)
    {
#ifdef HAVE_PEZY
      sw[i].pzmm = dev ? dev->pzmm : NULL;
      sw[i].xE   = dev ? dev->xE   : NULL;
      sw[i].nq   = dev ? dev->nq   : 0;
      sw[i].qi   = dev ? dev->qi   : 0;
#endif
      sw[i].status = eslOK;
    }
  p->next = 0;                 /* every query starts the cursor over */

  if (p->nworker == 1) scan_worker_query(&sw[0]);
  else
    {
      pthread_mutex_lock(&p->mtx);
      p->nrunning = p->nworker;
      p->gen++;
      pthread_cond_broadcast(&p->go);
      while (p->nrunning > 0) pthread_cond_wait(&p->done, &p->mtx);
      pthread_mutex_unlock(&p->mtx);
    }

  status = eslOK;
  for (i = 0; i < p->nworker; i++) if (sw[i].status != eslOK) status = sw[i].status;
  return (status == eslOK) ? eslEOF : status;   /* the callers treat eslEOF as a clean finish */

 ERROR:
  if (p) { free(p->tid); free(p->sw); free(p); }
  *ret_pool = NULL;
  return status;
}

/* Queries per device batch.  The device scores every (model, query) pair in one
 * launch, so a batch amortizes the launch and the per-call uploads; the cap
 * keeps the result array (nmodel x batch bytes) bounded. */
static int
scan_qbatch(void)
{
  const char *v = getenv("HMMER_HMMCACHE_QBATCH");
  if (v) { int n = atoi(v); if (n > 0) return n; }
  return 512;
}

/* How many queries make the one-off cache load worth it, and how big a
 * database we are willing to hold. */
static int
scan_cache_minq(void)
{
  const char *v = getenv("HMMER_HMMCACHE_MINQ");
  if (v) return atoi(v);
  return 4;
}
static int64_t
scan_cache_maxbytes(void)
{
  const char *v = getenv("HMMER_HMMCACHE_MAX");
  if (v) { long long n = atoll(v); if (n > 0) return (int64_t) n; }
  return (int64_t) 64 * 1024 * 1024 * 1024;
}

/* Count query sequences, cheaply, and leave <sqfp> rewound.  Stops counting at
 * <upto> since that is all the decision needs.  Returns -1 if the file cannot
 * be rewound, in which case we must not touch it. */
static int
scan_count_queries(ESL_SQFILE *sqfp, ESL_ALPHABET *abc, int upto)
{
  ESL_SQ *tmp = NULL;
  int     n   = 0;

  if (! esl_sqfile_IsRewindable(sqfp)) return -1;
  if ((tmp = esl_sq_CreateDigital(abc)) == NULL) return -1;
  while (n < upto && esl_sqio_ReadInfo(sqfp, tmp) == eslOK) { n++; esl_sq_Reuse(tmp); }
  esl_sq_Destroy(tmp);
  if (esl_sqfile_Position(sqfp, 0) != eslOK) return -1;
  return n;
}

/* Estimated resident size of the cache: the pressed files are within a few
 * percent of the in-memory profiles they deserialize into. */
static int64_t
scan_cache_estimate(const char *hmmfile)
{
  struct stat st;
  char   *path = NULL;
  int64_t tot  = 0;
  const char *sfx[2] = { ".h3f", ".h3p" };
  int     i;

  for (i = 0; i < 2; i++)
    {
      if (esl_sprintf(&path, "%s%s", hmmfile, sfx[i]) != eslOK) return -1;
      if (stat(path, &st) != 0) { free(path); return -1; }
      tot += (int64_t) st.st_size;
      free(path); path = NULL;
    }
  return tot;
}


static int
serial_loop(WORKER_INFO *info, P7_HMMFILE *hfp)
{
  int            status;

  P7_OPROFILE   *om;
  ESL_ALPHABET  *abc = NULL;
  /* Main loop: */
  while ((status = p7_oprofile_ReadMSV(hfp, &abc, &om)) == eslOK)
    {
      p7_pli_NewModel(info->pli, om, info->bg);
      p7_bg_SetLength(info->bg, info->qsq->n);
      p7_oprofile_ReconfigLength(om, info->qsq->n);

      status = p7_Pipeline(info->pli, om, info->bg, info->qsq, NULL, info->th);
      if (status == eslEINVAL) p7_Fail(info->pli->errbuf);

      p7_oprofile_Destroy(om);
      p7_pipeline_Reuse(info->pli);
    }

  esl_alphabet_Destroy(abc);

  return status;
}

#ifdef HMMER_THREADS
static int
thread_loop(ESL_THREADS *obj, ESL_WORK_QUEUE *queue, P7_HMMFILE *hfp)
{
  int  status   = eslOK;
  int  sstatus  = eslOK;
  int  eofCount = 0;
  P7_OM_BLOCK   *block;
  ESL_ALPHABET  *abc = NULL;
  void          *newBlock;

  esl_workqueue_Reset(queue);
  esl_threads_WaitForStart(obj);

  status = esl_workqueue_ReaderUpdate(queue, NULL, &newBlock);
  if (status != eslOK) esl_fatal("Work queue reader failed");
      
  /* Main loop: */
  while (sstatus == eslOK)
    {
      block = (P7_OM_BLOCK *) newBlock;
      sstatus = p7_oprofile_ReadBlockMSV(hfp, &abc, block);
      if (sstatus == eslEOF)
	{
	  if (eofCount < esl_threads_GetWorkerCount(obj)) sstatus = eslOK;
	  ++eofCount;
	}
	  
      if (sstatus == eslOK)
	{
	  status = esl_workqueue_ReaderUpdate(queue, block, &newBlock);
	  if (status != eslOK) esl_fatal("Work queue reader failed");
	}
    }

  status = esl_workqueue_ReaderUpdate(queue, block, NULL);
  if (status != eslOK) esl_fatal("Work queue reader failed");

  if (sstatus == eslEOF)
    {
      /* wait for all the threads to complete */
      esl_threads_WaitForFinish(obj);
      esl_workqueue_Complete(queue);  
    }
  
  esl_alphabet_Destroy(abc);
  return sstatus;
}

static void 
pipeline_thread(void *arg)
{
  int i;
  int status;
  int workeridx;
  WORKER_INFO   *info;
  ESL_THREADS   *obj;
  P7_OM_BLOCK   *block;
  void          *newBlock;
  
  impl_Init();

  obj = (ESL_THREADS *) arg;
  esl_threads_Started(obj, &workeridx);

  info = (WORKER_INFO *) esl_threads_GetData(obj, workeridx);

  status = esl_workqueue_WorkerUpdate(info->queue, NULL, &newBlock);
  if (status != eslOK) esl_fatal("Work queue worker failed");

  /* loop until all blocks have been processed */
  block = (P7_OM_BLOCK *) newBlock;
  while (block->count > 0)
  {
      /* Main loop: */
    for (i = 0; i < block->count; ++i)
    {
      P7_OPROFILE *om = block->list[i];

      p7_pli_NewModel(info->pli, om, info->bg);
      p7_bg_SetLength(info->bg, info->qsq->n);
      p7_oprofile_ReconfigLength(om, info->qsq->n);

      status = p7_Pipeline(info->pli, om, info->bg, info->qsq, NULL, info->th);
      if (status == eslEINVAL) p7_Fail(info->pli->errbuf);

      p7_oprofile_Destroy(om);
      p7_pipeline_Reuse(info->pli);

      block->list[i] = NULL;
    }

    status = esl_workqueue_WorkerUpdate(info->queue, block, &newBlock);
    if (status != eslOK) esl_fatal("Work queue worker failed");

    block = (P7_OM_BLOCK *) newBlock;
  }

  status = esl_workqueue_WorkerUpdate(info->queue, block, NULL);
  if (status != eslOK) esl_fatal("Work queue worker failed");

  esl_threads_Finished(obj, workeridx);
  return;
}
#endif   /* HMMER_THREADS */



