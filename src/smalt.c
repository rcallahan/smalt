/** SMALT main program */

/****************************************************************************
 ****************************************************************************
 *                                                                          *
 *  Copyright (C) 2010 - 2014 Genome Research Ltd.                          *
 *                                                                          *
 *  Author: Hannes Ponstingl (hp3@sanger.ac.uk)                             *
 *                                                                          *
 *  This file is part of SMALT.                                             *
 *                                                                          *
 *  SMALT is free software: you can redistribute it and/or modify it under  *
 *  the terms of the GNU General Public License as published by the Free    *
 *  Software Foundation, either version 3 of the License, or (at your       *
 *  option) any later version.                                              *
 *                                                                          *
 *  This program is distributed in the hope that it will be useful, but     *
 *  WITHOUT ANY WARRANTY; without even the implied warranty of              *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU        *
 *  General Public License for more details.                                *
 *                                                                          *
 *  You should have received a copy of the GNU General Public License along *
 *  with this program.  If not, see <http://www.gnu.org/licenses/>.         *
 *                                                                          *
 ****************************************************************************
 ****************************************************************************/
//#define smalt_debug
//#define smalt_debug_thread_single
//#define SMALT_DEBUG_XALI
#define SMALT_TIMING

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include <time.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "elib.h"
#include "threads.h"
#include "randef.h"
#include "sequence.h"
#include "infmt.h"
#include "hashidx.h"
#include "rmap.h"
#include "insert.h"
#include "report.h"
#include "results.h"
#include "resultpairs.h"
#include "menu.h"

enum {
  SMALT_REFSEQ_BLKZ = 64*1024*1024,
  SMALT_DEFAULT_SEQSET_FLAGS = SEQSET_COMPRESSED,
  SMALT_TARGET_DEPTH = 512,
  SMALT_MAX_DEPTH = 2048,
#ifdef SMALT_DEBUG_XALI
  SMALT_MAX_REFSEQ_NUM = 1,
#else
  SMALT_MAX_REFSEQ_NUM = 512, /**< If there are more than this number of reference sequences,
				*  lump them all together and find out reference seqence index
				*  later. */
#endif
  SMALT_MAX_THREAD_NUM = 240,  /**< Maximum number of threads to be spawned */
  SMALT_THREAD_ARGFAC_SRT  = 4,/**< Number of buffered thread arguments/results, when output
				* should in order of input, as a multiplicative factor. 
				* I.e. for n>2 threads there are
				* SMALT_THREAD_ARGFAC*n arguments buffered. */
  SMALT_THREAD_ARGFAC_NOSRT = 2,
  SMALT_THREAD_NOARG = -1,    /**< negative thread number signals empty thread stack */
  SMALT_MAXNBITS_PERF = 10,   /**< Maximum number of bits for perfect part of hash key */
  SMALT_MAXNBITS_KEY = 26,    /**< Maximum number of bits for hash key */
  SMALT_NBITS_KEY_MARG = 1,   /**< Minimum number of bits for hash key = perfect bits + SMALT_NBITS_KEY_MARG */
  SMALT_UINT32_NBITS = 32,
  SMALT_INSHISTO_SIZ = 256,   /**< Number of sampling points for histogram of insert size distributions */
  SMALT_INSHISTO_SEEDNUM = SMALT_INSHISTO_SIZ*4, /**< Number of initial seeds (dummy counts) */
  SMALT_LINWIDTH = 80,       /* line width for histogram of insert sizes */
#ifdef smalt_debug
  SMALT_INSHIST_IVAL = 1000, /* report sampled histogram every SMALT_INSHIST_IVAL lines */
#endif
  SMALT_INSSAMPLE_BLKSZ = 2048,     /* Block size for memory allocation of isert size sample */
  SMALT_INSSAMPLE_MINSIZ = 1028,   /* minimum sample size */
  SMALT_NARGS_PER_THREAD = 32,      /* number of per-read buffers in block per thread */
};

enum THREADARG_FLAGS {
  ARGFLG_FILLED = 0x01,   /**< Argument is filled with a read */
  ARGFLG_INUSE = 0x02,    /**< in use by thread */
  ARGFLG_PRINTED = 0x04,  /**< arguments have been read out */
};

#ifdef smalt_debug
enum DEBUG_BUFFTYP {
  BUFFTYP_UNKNOWN,
  BUFFTYP_EMPTY,
  BUFFTYP_LOADED,
  BUFFTYP_MAPPED,
};
#endif

typedef struct _SmaltMapConst { /**< Constant Arguments for mapping function (for multiple threads )*/
  MENUFLG_t menuflg;    /**< Menu flags (combination of MENU_FLAGS) */
  char subprogtyp;     /**< One of MENU_SUBPROG */
  uint8_t inform;      /**< One of MENU_INPUT_FORMATS */
  RMAPFLG_t rmapflg;   /**< Flags tuning mapping strategies (combination of RMAP_FLAGS) */
  RESULTOUTFLG_t rsltouflg;/**< Output flags (combination of RESULTSET_OUTPUT_FLAGS) */
  const char *oufilnam;/**< File name for output */
  REPOUFMT_t outform;     /**< Output format (one of REPORT_OUTPUT_FORMATS) */
  REPMODIFLG_t oumodflg;   /**< Combination of REPORT_MODIFIER_FLAGS */
  SeqCodec *codecp;          /**< Sequence De-/Ent-coder */
  ScorePenalties *scorpltyp; /**< Alignment penalty scores */
  ScoreMatrix *scormtxp;     /**< Scoring matrix (SW alignment scores) */
  HashTable *htp;            /**< Hash Table */
  SeqSet *ssp;               /**< Set of reference sequences */
  ResultFilter *rfp;         /**< Data for filtering output of results */
  InsHist *ihp;              /**< Histogram of insert size distributions */
  int insert_min;            /**< Minimum insert size (paired-end reads) */
  int insert_max;            /**< Maximum insert size (paired-end reads) */
  int nhitmax_tuple;         /**< Maximum number of hits a k-mer word may have across the genome.
			      * k-mer words with hits above that cut-off are disregarded */
  int min_swatscor;          /**< Mimimum Smith-Waterman score for a mapping to be reported */
  int swatscordiff;          /**< Maximum allowed difference to the highest Smith-Waterman score
			      * for a mapping to be reported */
  unsigned char minbasq;     /**< Base quality threshold. K-mer words with bases of qualities below 
			      * this threshold are disregarded. */
  double tupcovmin;          /**< Minimum cover of reads by k-tuples (either as fraction of readlenght
			      * or as number of bases > */
  int readskip;              /**< Step size for skipping greads (when sampling insert sizes). */
  RSLTPAIRLIB_t pairtyp;     /**< Type of read pair library (one of RSLTPAIR_LIB) */
  InFmtReader *ifrp;         /**< Memory allocation point for reader of read-sequences */
  const char *prognam;       /**< Name of main program */
  const char *progversion;   /**< Version of main program */
  int cmdlin_narg;           /**< Number of arguments in original command line */
  char **cmdlin_argv;        /**< Pointer to argument in original command line */
  short threadblksz;         /**< Number of per-read buffers SmaltIOBuffArg in thread block 
			      * SmaltArgBlock */
} SmaltMapConst;

typedef struct _SmaltIOBuffArg { /**< Arguments to be buffered and passed between input,
				  * worker and output threads */
  uint64_t readno;
  unsigned char isPaired;
  SeqFastq *readp;
  SeqFastq *matep;
#ifdef RESULTS_TRACKER
  Track *trkrp;
  Track *trkmp;
#endif
  int isiz;
  RSLTPAIRFLG_t pairflg;
  Report *rep;
} SmaltIOBuffArg;

typedef struct SmaltArgBlock_ {
  short argno;
  short n_iobf;          /**< number of active element 
			  * in block iobfp */
  short n_alloc;         /**< Number of elements allocated in iobfp */        
  SmaltIOBuffArg *iobfp; /**< pointer to first in block */
} SmaltArgBlock;

typedef struct _SmaltMapArgs { /**< Arguments for mapping function
				* (needed a single structure for multiple threads) */
  short threadno;              /**< Thread number */
  const SmaltMapConst *smconstp;/**< stuff that is constant across all reads */
  RMap *rmp;                   /**< Various buffer structures used during mapping */
  int errcode;                 /**< Error code (one of ERRMSG_CODES) for the thread */
#ifdef hashhit_dump_sortarray
  FILE *dumpfp;
  int dumpctr;
#endif
} SmaltMapArgs;

typedef struct _SmaltInput { /**< Arguments for reading function */
  int errcode;
  short threadno;
  MENUFLG_t menuflg;
  uint64_t rctr;   /**< read counter */
  uint64_t pctr;   /**< paired read counter */
  int rival;       /**< sampling interval */
  InFmtReader *ifrp;
} SmaltInput;

typedef struct _SmaltOutput {
  short threadno;
  uint64_t next_readno;   /**< Number of read to be printed next */
  const SmaltMapConst *smcp;
  InsSample *isamp;
  ReportWriter *writerp;
} SmaltOutput;


/*****************************************************************************
 ************************ Globals for multiple threads ***********************
 *****************************************************************************/

#ifdef hashhit_dump_sortarray
static const char smalt_helper_oufilnam_fmt[] = "smalt_%2.2i.bin";
#endif

/*****************************************************************************
 ****************************** misc routines ********************************
 *****************************************************************************/
static REPOUFMT_t convertOutputFormat(REPMODIFLG_t *modiflg, uint8_t form, const MenuOpt *mp)
{
  char rv;
  switch (form) {
  case MENU_OUTFORM_CIGAR:
    rv = REPORTFMT_CIGAR;
    break;
  case MENU_OUTFORM_SAM:
    rv = REPORTFMT_SAM;
    if (!menuTestMapOutputFormatFlags(mp, form, MENU_OUFMTSAM_CLIPPED))
      *modiflg |= REPORTMODIF_SOFTCLIP;
    if (menuTestMapOutputFormatFlags(mp, form, MENU_OUFMTSAM_HEADER))
      *modiflg |= REPORTMODIF_HEADER;
    if (menuTestMapOutputFormatFlags(mp, form, MENU_OUFMTSAM_XCIGAR))
      *modiflg |= REPORTMODIF_XMISMATCH;    
    break;
  case MENU_OUTFORM_BAM:
    rv = REPORTFMT_BAM;
    if (!menuTestMapOutputFormatFlags(mp, form, MENU_OUFMTSAM_CLIPPED))
      *modiflg |= REPORTMODIF_SOFTCLIP;
    if (menuTestMapOutputFormatFlags(mp, form, MENU_OUFMTSAM_HEADER))
      *modiflg |= REPORTMODIF_HEADER;
    if (menuTestMapOutputFormatFlags(mp, form, MENU_OUFMTSAM_XCIGAR))
      *modiflg |= REPORTMODIF_XMISMATCH;    
    break;
  case MENU_OUTFORM_SSAHA:
    rv = REPORTFMT_SSAHA;
    break;
  case MENU_OUTFORM_GFF:
    rv = REPORTFMT_GFF2;
    break;
  default:
    rv = REPORTFMT_CIGAR;
    break;
  }
  return rv;
}

static INFMT_t convertInputFormat(unsigned char form)
{
  INFMT_t rv;
  switch (form) {
  case MENU_INFORM_FASTQ:
    rv = INFMT_FASTQ;
    break;
  case MENU_INFORM_SAM:
    rv = INFMT_SAM;
    break;
  case MENU_INFORM_BAM:
    rv = INFMT_BAM;
    break;
  case MENU_INFORM_UNKNOWN:
  default:
    rv = INFMT_UNKNOWN;
    break;
  }
  return rv;
}

static int selectHashTyp(uint8_t *typ, uint8_t *nbits_key, uint8_t *nbits_perf, 
			 uint8_t wordlen, uint8_t nskip, 
			 const SeqSet *ssp)
{
  SEQNUM_t nseq;
  unsigned short nbk = (unsigned short)(wordlen<<1);
  size_t ntup;
  const SETSIZ_t *offsp;
  uint64_t nkey;

  *nbits_key = 0;
  *nbits_perf = 0;
  *typ = HASHIDXTYP_PERFECT;

  if (nbk > 63)
    return ERRCODE_MAXKTUP;

  if (nskip < 1) 
    nskip = 1;

  nseq = seqSetGetOffsets(ssp, &offsp);
  ntup = offsp[nseq]/nskip;
  nkey = ((uint64_t) 1)<<nbk;

  if (ntup > UINT32_MAX)
    return ERRCODE_SKIPSMALL;

  if (nkey > 2*ntup) {
    uint8_t last_b, i;
    uint32_t t;
    *typ = HASHIDXTYP_HASH32MIX;
   
    /* calculate last_b as the smallest 2^last_b > 2*ntup */
    last_b = (uint8_t) ((ntup & 0x01)? 1:0);
    for (t=ntup, i=0; i<32; i++) {
      t = t>>1;
      if (t&1)
	last_b = i;
    }
    if ((last_b & 1))
      *nbits_key = (uint8_t) (last_b + 1);
    else 
      *nbits_key = last_b;

    if (nbk > SMALT_UINT32_NBITS) {
      *nbits_perf = (uint8_t) (nbk - SMALT_UINT32_NBITS);
      if (*nbits_perf > SMALT_MAXNBITS_PERF)
	return ERRCODE_MAXKTUP;
    } else {
      *nbits_perf = 0;
    }

    if (*nbits_key + *nbits_perf > SMALT_MAXNBITS_KEY) {
      *nbits_key = (uint8_t) (SMALT_MAXNBITS_KEY - *nbits_perf);
    }

    if (*nbits_key < *nbits_perf + SMALT_NBITS_KEY_MARG)
      *nbits_key = (uint8_t) (*nbits_perf + SMALT_NBITS_KEY_MARG);

    if (*nbits_key > SMALT_MAXNBITS_KEY)
      *nbits_key = SMALT_MAXNBITS_KEY;
  }

  return ERRCODE_SUCCESS;
}

/*****************************************************************************
 ************************** calculating hash index ***************************
 *****************************************************************************/

static int buildHashIndex(ErrMsg *errmsgp, const MenuOpt *menup)
{
  int errcode;
  const char is_verbose = 1;
  uint8_t ktup, nskip, nbits_key, nbits_perf, typ;
  const char *hashnam, *filnam;
  HashTable *htp;
  SeqCodec *codecp;
  SeqFastq *sqbufp;
  SeqSet *ssp;

  menuGetFileNames(menup, &filnam, NULL);
  
  if((errcode = menuGetHashParams(menup, &hashnam, &ktup, &nskip)))
    ERRMSGNO(errmsgp, errcode);
  
  if (!(codecp = seqCodecCreate()))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);
  
  if (!(sqbufp = seqFastqCreate(SMALT_REFSEQ_BLKZ, SEQTYP_UNKNOWN)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);
 
  if (!(ssp = seqSetCreate(SMALT_REFSEQ_BLKZ, SMALT_DEFAULT_SEQSET_FLAGS)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  if ((is_verbose)) 
    fprintf(stderr, "# Reading sequences ...\n");

  seqSetAddFromFastqFile(errmsgp, ssp, sqbufp, codecp, 
			 filnam, 0);
  
/*   if ((is_verbose)) */
/*     fprintf(stderr, "# Compressing sequence set ...\n"); */
/*   if ((errcode = seqSetCompress(ssp, codecp))) */
/*     ERRMSGNO(errmsgp, errcode); */

  if ((is_verbose))
    fprintf(stderr, "# Writing sequence set ...\n");
  if ((errcode = seqSetWriteBinFil(ssp, hashnam)))
    ERRMSGNO(errmsgp, errcode);
  
  if ((errcode = selectHashTyp(&typ, &nbits_key, &nbits_perf, ktup, nskip, ssp)))
    ERRMSGNO(errmsgp, errcode);

  if ((is_verbose)) {
    if (typ == HASHIDXTYP_PERFECT)
      fprintf(stderr, "# Setting up perfect hash index ...\n");
    else
      fprintf(stderr, "# Setting up hash index with collisions ...\n");
    
    fprintf(stderr, "# word length = %i bases, skip step = %i bases ...\n", (int) ktup, (int) nskip);
    if (typ != HASHIDXTYP_PERFECT)
      fprintf(stderr, "# number of bits for key = %i with %i perfect bits\n", (int) nbits_key, (int) nbits_perf);
  }
  if (!(htp = hashTableCreate(ktup, nskip, nbits_key, nbits_perf, typ)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  if ((errcode = hashTableSetUp(htp, sqbufp, ssp, NULL, codecp, NULL, 
				is_verbose))) 
    ERRMSGNO(errmsgp, errcode);

  if ((is_verbose)) 
    hashTablePrintStats(stderr, htp);
  
  if ((is_verbose))
    fprintf(stderr, "# Writing table to file ... \n");
  if ((errcode = hashTableWrite(hashnam, htp)))
    ERRMSGNO(errmsgp, errcode);
  
  hashTableDelete(htp);
  seqSetDelete(ssp);
  seqFastqDelete(sqbufp);
  seqCodecDelete(codecp);
  return ERRCODE_SUCCESS;
}

/*****************************************************************************
 ******************** Methods of Private Type SmaltMapConst ******************
 *****************************************************************************/
static void updateInsertBoundariesFromSample(SmaltMapConst *macop, const InsHist *ihp)
{
  int32_t lo, hi;
  insGetHistoData(&lo, &hi, NULL, ihp);
  if (lo < macop->insert_min)
    macop->insert_min = lo;
  if (hi > macop->insert_max)
    macop->insert_max = hi;
  return;
}

static void cleanupMapConst(SmaltMapConst *macop)
{
  infmtDeleteReader(macop->ifrp);
  hashTableDelete(macop->htp);
  seqSetDelete(macop->ssp);
  scoreDeleteMatrix(macop->scormtxp);
  scorePenaltiesDelete(macop->scorpltyp);
  seqCodecDelete(macop->codecp);
  resultSetDeleteFilter(macop->rfp);
  insDeleteHisto(macop->ihp);
  macop->htp = NULL;
  macop->ssp = NULL;
  macop->codecp = NULL;
}

static int initMapConst(SmaltMapConst *smcp,
			const MenuOpt *menup)
{
  int errcode;
  int seed;
  short min_swatscor;
  short swatscordiff;
  double idmin;
  const char *indexnam, *insfilnam, *filnamA, *filnamB, *tmpdir;
  uint8_t menuoutform, informat;
  unsigned char pairtyp = 0;
  int ninfil;
  short nthreads = 0;
  INFMT_t fmt;

  memset(smcp, 0, sizeof(SmaltMapConst));

  smcp->subprogtyp = menuGetSubProgTyp(menup);
  smcp->cmdlin_argv = menuGetCommandLine(menup, &smcp->cmdlin_narg);
  smcp->prognam = menuGetProgramName(&smcp->progversion);
  nthreads = menuGetNumberOfThreads(menup);
  if (nthreads > SMALT_MAX_THREAD_NUM)
    nthreads = SMALT_MAX_THREAD_NUM;
  smcp->threadblksz = (short) ((nthreads > 0)? nthreads*SMALT_NARGS_PER_THREAD: 1);

  if((errcode = menuGetMapParams(menup, &indexnam, 
				 &smcp->nhitmax_tuple, &smcp->tupcovmin,
				 &min_swatscor,
				 &swatscordiff,
				 &smcp->minbasq,
				 &idmin,
				 &seed,
				 &smcp->readskip,
				 &smcp->insert_min, &smcp->insert_max,
				 &menuoutform,
				 &smcp->oufilnam,
				 &insfilnam,
				 &pairtyp)))
    return errcode;
  smcp->menuflg = menuGetFlags(menup);
  smcp->swatscordiff = (int) swatscordiff; 
  smcp->rfp = resultSetCreateFilter();
  if ((errcode = menuGetMapInputFormat(menup, &smcp->inform, NULL)))
    return errcode;

  if (!smcp->rfp)
    return ERRCODE_NOMEM;
  resultSetFilterData(smcp->rfp, (int) min_swatscor, (int) swatscordiff, idmin);
  smcp->outform = convertOutputFormat(&smcp->oumodflg, menuoutform, menup);
  if (smcp->menuflg&MENUFLAG_ALIGNMENT)
    smcp->oumodflg |= REPORTMODIF_ALIOUT;

  if (!(swatscordiff)) {
    smcp->rsltouflg |= RESULTFLG_BEST;
    smcp->rmapflg |= RMAPFLG_BEST;
    if (!(smcp->menuflg & MENUFLAG_RELSCOR)) {
      smcp->rsltouflg |= RESULTFLG_SINGLE;
      if (smcp->menuflg & MENUFLAG_RANDREPEAT) {
	smcp->rsltouflg |= RESULTFLG_RANDSEL;
	RANSEED(seed);
      }
    }
  }

  if (smcp->menuflg&MENUFLAG_SPLITREAD) {
    smcp->rmapflg |= RMAPFLG_SPLIT | RMAPFLG_NOSHRTINFO | RMAPFLG_SENSITIVE;
    //smcp->rmapflg &= ~RMAPFLG_BEST;
    smcp->rsltouflg |= RESULTFLG_SPLIT;
  } 

  if (smcp->menuflg&MENUFLAG_COMPLEXW)
    smcp->rmapflg |= RMAPFLG_CMPLXW;

  if (smcp->menuflg&MENUFLAG_PAIRED && 
      (smcp->inform == MENU_INFORM_FASTQ || smcp->inform == MENU_INFORM_UNKNOWN))
    smcp->rmapflg |= RMAPFLG_PAIRED; /* implies two separate FASTA/FASTQ input files */

  if (pairtyp == MENU_READPAIRTYP_PAIREDEND )
    smcp->pairtyp = RSLTPAIRLIB_PAIREDEND;
  else if (pairtyp == MENU_READPAIRTYP_MATEPAIR)
    smcp->pairtyp = RSLTPAIRLIB_MATEPAIR;
  else if (pairtyp == MENU_READPAIRTYP_SAMESTRAND)
    smcp->pairtyp = RSLTPAIRLIB_SAMESTRAND;
  else if (pairtyp == MENU_READPAIRTYP_UNKNOWN)
    smcp->pairtyp = RSLTPAIRLIB_PAIREDALL;
  else 
    return ERRCODE_MENUE;
  
  if (smcp->menuflg & MENUFLAG_EXHAUSTIVE) {
    smcp->rmapflg |=  RMAPFLG_NOSHRTINFO | RMAPFLG_SENSITIVE | RMAPFLG_ALLPAIR;
  }

  if (!(smcp->scorpltyp = scorePenaltiesCreate()))
    return ERRCODE_NOMEM;

  /* modify sw penalties here */
  if (menuGetMapPenaltyScores(menup, NULL, NULL, NULL, NULL) != 0) {
    int8_t match, subst, gap_open, gap_ext;
    uint8_t bitflg = menuGetMapPenaltyScores(menup, &match, &subst, &gap_open, &gap_ext);
    if (bitflg & 0x1)
      scoreSetPenalty(smcp->scorpltyp, SCORPNLTYP_MATCH, (short) match);
    if (bitflg & 0x2)
      scoreSetPenalty(smcp->scorpltyp, SCORPNLTYP_MISMATCH, (short) subst);
    if (bitflg & 0x4)
      scoreSetPenalty(smcp->scorpltyp, SCORPNLTYP_GAPOPEN, (short) gap_open);
    if (bitflg & 0x8)
      scoreSetPenalty(smcp->scorpltyp, SCORPNLTYP_GAPEXT, (short) gap_ext);
  }

  if (!(smcp->codecp = seqCodecCreate()) ||
      !(smcp->scormtxp = scoreCreateMatrix(smcp->codecp, smcp->scorpltyp)))
    return ERRCODE_NOMEM;

  if ((insfilnam)) {
    if (smcp->menuflg & MENUFLAG_VERBOSE)
      fprintf(stderr, "# Reading distribution of insert sizes from file ...\n");
    smcp->ihp = insReadHisto(&errcode, insfilnam);
    if ((errcode)) {
      smcp->ihp = NULL;
      return errcode;
    }
    if (smcp->menuflg & MENUFLAG_VERBOSE) {
      fprintf(stderr, "# Sampled histogram\n");
      insPrintHisto(stdout, 80, 0, smcp->ihp);
      fprintf(stderr, "# Smoothed histogram\n");
      insPrintHisto(stdout, 80, 1, smcp->ihp);
    }
    updateInsertBoundariesFromSample(smcp, smcp->ihp);
  }

  /* open read files first (should take shorter time than index files) 
     rather that in initInput */
  if ((ninfil = menuGetFileNames(menup, &filnamA, &filnamB)) < 1)
    return ERRCODE_ASSERT;

  if ((errcode = menuGetMapInputFormat(menup, &informat, &tmpdir)))
    return errcode;
  fmt = convertInputFormat(informat);

  if (smcp->menuflg & MENUFLAG_VERBOSE)
    fprintf(stderr, "# Opening read %s ...\n", (ninfil == 2)? "file":"files");
  smcp->ifrp = infmtCreateReader(&errcode, 
				 filnamA, filnamB,
#ifdef HAVE_BAMBAMC
				 tmpdir,
#endif
				 fmt);
  if ((errcode))
    return errcode;

  if (smcp->menuflg & MENUFLAG_VERBOSE)
    fprintf(stderr, "# Reading reference sequences ...\n");
  smcp->ssp = seqSetReadBinFil(&errcode, indexnam);
  if ((errcode)) 
    return errcode;

  if (seqSetGetOffsets(smcp->ssp, NULL) < SMALT_MAX_REFSEQ_NUM)
	smcp->rmapflg |= RMAPFLG_SEQBYSEQ;
  
  if (smcp->menuflg & MENUFLAG_VERBOSE)
    fprintf(stderr, "# Reading hash table ...\n");
  smcp->htp = hashTableRead(&errcode, indexnam);
  if ((errcode))
    return errcode;
  
  if (smcp->menuflg&MENUFLAG_MINSCOR) {
    smcp->min_swatscor = (int) min_swatscor;
  } else {
    /* no minimum score specified on command line -> set default */
    uint8_t nskip =0;
    uint8_t ktup = hashTableGetKtupLen(smcp->htp, &nskip);
    smcp->min_swatscor = (int) ktup + nskip - 1;
  }

  if (smcp->menuflg & MENUFLAG_VERBOSE)
    hashTablePrintStats(stderr, smcp->htp);

  return ERRCODE_SUCCESS;
}

/*****************************************************************************
 ******************** Methods of Private Type SmaltOutput *********************
 *****************************************************************************/
/* THREAD_INITF */
static int initSmaltOutput(void *op, const void *mp, short threadno)
{
  int errcode = ERRCODE_SUCCESS;
  SmaltOutput *dop = (SmaltOutput *) op;
  const SmaltMapConst *mcp= (const SmaltMapConst *) mp;

  dop->threadno = threadno;
  dop->writerp = reportCreateWriter(&errcode,
				    mcp->oufilnam, 
				    mcp->outform, mcp->oumodflg,
				    mcp->ssp,
				    mcp->prognam,
				    mcp->progversion,
				    mcp->cmdlin_argv,
				    mcp->cmdlin_narg
				    );
  if ( (!errcode) && MENU_SAMPLE == mcp->subprogtyp ) {
    dop->isamp = insCreateSample(SMALT_INSSAMPLE_BLKSZ);
    if (dop->isamp == NULL)
      errcode = ERRCODE_NOMEM;
  } else {
    dop->isamp = NULL;
  }
  if (!errcode) {
    dop->next_readno = 0;
    dop->smcp = mcp;
  }

  return errcode;
}

/* THREAD_CLEANF */
static int cleanupSmaltOutput(ErrMsg *errmsgp, void *op)
{
  int errcode = ERRCODE_SUCCESS;
  SmaltOutput *dop = (SmaltOutput *) op;

  insDeleteSample(dop->isamp);
  dop->isamp = NULL;
  reportDeleteWriter(dop->writerp);
  dop->writerp = NULL;
  return errcode;
}
/*****************************************************************************
 ********************* Methods of Private Type SmaltInput ********************
 *****************************************************************************/
/* THREAD_INITF */
static int initInput(void *ip, const void *mp, short threadno)
{
  SmaltInput *inargp = (SmaltInput *) ip;
  const SmaltMapConst *mcp = (const SmaltMapConst *) mp;

  inargp->threadno = threadno;
  inargp->errcode = ERRCODE_SUCCESS;
  inargp->menuflg = mcp->menuflg;
  inargp->rctr = 0LL;
  inargp->pctr = 0LL;
  inargp->rival = 0;
  inargp->ifrp = mcp->ifrp;

  return inargp->errcode;
}

/* THREAD_CLEANF */
static int cleanupInput(ErrMsg *errmsgp, void *inp)
{
  SmaltInput *inargp = (SmaltInput *) inp;
  int errcode = infmtGetReaderStatus(inargp->ifrp);

  if (inargp->menuflg & MENUFLAG_VERBOSE) {
    if (inargp->pctr > 0) {
      fprintf(stderr, "# Processed %llu read pairs", 
	     (unsigned long long) inargp->pctr);
      if (inargp->rctr > inargp->pctr) {
	fprintf (stderr, "\n# and %llu single reads.\n", 
		 (unsigned long long) (inargp->rctr - inargp->pctr));
      } else {
	fprintf (stderr, ".\n");
      }
    } else {
      fprintf (stderr, "# Processed %llu single reads.\n", 
	       (unsigned long long) inargp->rctr);
    }
  } 
      
  if (errcode && errcode != ERRCODE_EOF) {
    inargp->errcode = errcode;
    ERRMSGNO(errmsgp, errcode);
  }

  inargp->ifrp = NULL;
  inargp->rctr = 0LL;
  inargp->pctr = 0LL;
  inargp->rival = 0;
  
  return (inargp->errcode = errcode);
  }

static void setInSampleIntervalInput(SmaltInput *inargp, 
				     const InsSample *isamp) 
{
  if (isamp == NULL) {
    inargp->rival = 0;
  } else {
    insGetSample(NULL, &inargp->rival, isamp);
  }
}

/*****************************************************************************
 ****************** Methods of Private Type SmaltIOBuffArg *******************
 *****************************************************************************/
/* THREAD_INITF callback function */
static int initIOBuffArg(SmaltIOBuffArg *argp, 
			 uint8_t prep_paired)
{
  int errcode = ERRCODE_SUCCESS;

  argp->readp = seqFastqCreate(0, SEQTYP_UNKNOWN);
  if (NULL == argp->readp) 
    errcode = ERRCODE_NOMEM;
  if ((prep_paired) && ERRCODE_SUCCESS == errcode) {
    argp->matep = seqFastqCreate(0, SEQTYP_UNKNOWN);
    if (NULL == argp->matep)
      errcode =  ERRCODE_NOMEM;
  } else {
    argp->matep = NULL;
  }
  if (ERRCODE_SUCCESS == errcode) {
    argp->rep = reportCreate(0);
    if (NULL == argp->rep)
      errcode = ERRCODE_NOMEM;
  }
#ifdef RESULTS_TRACKER
  if (ERRCODE_SUCCESS == errcode) {
    argp->trkrp = trackCreate();
    if (NULL == argp->trkrp)
      errcode = ERRCODE_NOMEM;
    else if ((prep_paired)) {
      argp->trkmp = trackCreate();
      if (NULL == argp->trkmp)
	errcode = ERRCODE_NOMEM;
    }
  }
#endif
  argp->readno = 0LL;
  argp->isiz = 0;
  argp->pairflg = 0;

  return errcode;
}

static int cleanupIOBuffArg(ErrMsg *errmsgp, SmaltIOBuffArg *argp)
{
  if (NULL == argp) {
    ERRMSGNO(errmsgp, ERRCODE_NULLPTR);
  } else {
    seqFastqDelete(argp->readp);
    seqFastqDelete(argp->matep);
    reportDelete(argp->rep);
#ifdef RESULTS_TRACKER
    trackDelete(argp->trkrp);
    trackDelete(argp->trkmp);
#endif
  }

  return ERRCODE_SUCCESS;
}

static int loadIOBuffArg(ErrMsg *errmsgp,
			 SmaltInput *inp, 
			 SmaltIOBuffArg *argp)
{
  int errcode = ERRCODE_SUCCESS;

  assert(inp);
  assert(argp);

  errcode = infmtRead(inp->ifrp, argp->readp, argp->matep, &argp->isPaired);
  
  if ((errcode) && errcode != ERRCODE_EOF)
    ERRMSGNO(errmsgp, errcode);

  while (!(errcode) && 
	 (inp->rival > 0 && (inp->pctr % inp->rival) != 0)) {
    errcode = infmtRead(inp->ifrp, argp->readp, argp->matep, &argp->isPaired);

    if (ERRCODE_SUCCESS == errcode) {
      if ((argp->isPaired)) 
	inp->pctr++;
      inp->rctr++;
    }
  }

  argp->readno = inp->rctr;
  argp->isiz = 0;
  argp->pairflg = 0;
  reportBlank(argp->rep);
  
  if (ERRCODE_SUCCESS == errcode) {
    if ((argp->isPaired)) 
      inp->pctr++;
    inp->rctr++;
  }
#ifdef smalt_debug
  fprintf(stderr, "#SMALT_DBUG:loadBuffArg():loaded read no %llu\n", 
	 argp->readno);
#endif
  return errcode;
}

/* THREAD_PROCF */
static int outputIOBuffArg(ErrMsg *errmsgp, 
			   SmaltOutput *dop, 
			   SmaltIOBuffArg *argp)
{
  int errcode = ERRCODE_SUCCESS;
  const SmaltMapConst *smcp;

  assert(dop != NULL);
  assert(argp != NULL);
  smcp = dop->smcp;
  assert(dop->smcp != NULL);

  if ((argp->pairflg & RSLTPAIRFLG_INSERTSIZ) &&
      MENU_SAMPLE == smcp->subprogtyp &&
      dop->isamp != NULL) {
    insAddSample(dop->isamp, argp->isiz);
  }
  dop->next_readno = argp->readno + 1;
  if (smcp->menuflg & MENUFLAG_RELSCOR && /* output of repetitive mappings */
      (smcp->outform == REPORTFMT_SAM ||
       smcp->outform == REPORTFMT_BAM))
    reportFixMultiplePrimary(argp->rep);
  errcode = reportWrite(dop->writerp, 
#ifdef RESULTS_TRACKER 
			argp->trkrp, argp->trkmp,
#endif
			argp->readp, argp->matep,
			smcp->ssp, smcp->codecp,
			argp->rep);
#ifdef smalt_debug
  fprintf(stderr, "#SMALT_DEBUG:outputBuffArg():wrote read no %llu\n", 
	  (long long unsigned int) argp->readno);
#endif
  if ((errcode))
    ERRMSGNO(errmsgp, errcode);

  return errcode;
}

/*****************************************************************************
 ******************** Methods of Private Type SmaltArgBLock ******************
 *****************************************************************************/

/* THREAD_INITF callback function */
static int initArgBlock(void *ap, const void *ip, short argno)
{
  int errcode = ERRCODE_SUCCESS;
  short i;
  SmaltArgBlock *blockp = (SmaltArgBlock *) ap;
  const SmaltMapConst *mcp = (const SmaltMapConst *) ip;
  uint8_t prep_paired = (uint8_t) ((mcp->rmapflg & RMAPFLG_PAIRED) || 
				   mcp->inform == MENU_INFORM_SAM ||
				   mcp->inform == MENU_INFORM_BAM);
  
  blockp->argno = argno;
  blockp->n_iobf = 0;
  if (NULL == ECALLOCP(mcp->threadblksz, blockp->iobfp)) {
    errcode = ERRCODE_NOMEM;
    blockp->n_alloc = 0;
  } else {
    blockp->n_alloc = mcp->threadblksz;
  }

  for (i=0; i<blockp->n_alloc && !(errcode); i++)
    errcode = initIOBuffArg(blockp->iobfp + i, prep_paired);
  
  return errcode;
}

/* THREAD_CLEANF */
static int cleanupArgBlock(ErrMsg *errmsgp, void *p)
{
  int errcode = ERRCODE_SUCCESS;
  short i;
  SmaltArgBlock *blockp = (SmaltArgBlock *) p;
 
  assert( p != NULL );

  for (i=0; i<blockp->n_alloc && !(errcode); i++)
    errcode = cleanupIOBuffArg(errmsgp, blockp->iobfp + i);

  if (!errcode) {
    free(blockp->iobfp);
    blockp->iobfp = NULL;
    blockp->n_iobf = 0;
    blockp->n_alloc = 0;
  }

  return errcode;
}

/* THREAD_PROCF */
static int loadArgBlock(ErrMsg *errmsgp,
#ifdef THREADS_DEBUG
			   uint64_t *readno,
#endif 
			   void *thargp, void *bufargp)
{
  int errcode = ERRCODE_SUCCESS;
  short i;
  SmaltInput *inp = (SmaltInput *) thargp;
  SmaltArgBlock *blockp = (SmaltArgBlock *) bufargp;

  assert(bufargp != NULL);
  blockp->n_iobf = 0;
  for (i=0; i<blockp->n_alloc; i++) {
    errcode = loadIOBuffArg(errmsgp, 
			    inp,
			    blockp->iobfp+i);
    if (errcode)
      break;
  }
  blockp->n_iobf = i;
#ifdef THREADS_DEBUG
  *readno = (i > 1)? blockp->iobfp->readno: 0;
  threadsPrintDebugMsg(NULL, NULL,
		       "loadArgBlock(): read read no %llu to %llu (%s) "\
		       "(errcode = %i)",
		       (long long unsigned int) *readno,
		       (long long unsigned int) (*readno) + i, 
		       (i == blockp->n_alloc)? "complete": "incomplete",
		       errcode);
#endif

  return (ERRCODE_EOF == errcode && i > 0)? ERRCODE_SUCCESS: errcode;
}

/* THREAD_CHECKF */
static int checkArgBlockReadNo(const void *thargp, const void *buffargp)
{
  const SmaltOutput *dop = (SmaltOutput *) thargp;
  const SmaltArgBlock *blockp = (SmaltArgBlock *) buffargp;
  uint64_t readno;
  assert(thargp != NULL);
  assert(buffargp != NULL);
  assert(blockp->iobfp != NULL);
  readno = blockp->iobfp->readno;
  
#ifdef SMALT_debug
  fprintf(stderr, "#SMALT_DEBUG:checkBuffArgBlockReadNo: %llu (next: %llu)", 
	 (unsigned long long) readno, 
	 (unsigned long long) dop->next_readno);
#endif
  return (blockp->n_iobf > 0 && readno <= dop->next_readno);
}

/* THREAD_CMPF */
static int cmpArgBlockReadNo(const void *buffargAp, 
			     const void *buffargBp)
{
  const SmaltArgBlock *bfAp = (const SmaltArgBlock *) buffargAp;
  const SmaltArgBlock *bfBp = (const SmaltArgBlock *) buffargBp;
  int rv = 0;

  assert(buffargAp != NULL && bfAp->iobfp != NULL);
  assert(buffargBp != NULL && bfBp->iobfp != NULL);

  if (bfAp->iobfp->readno < bfBp->iobfp->readno)
    rv = -1;
  else if (bfAp->iobfp->readno > bfBp->iobfp->readno)
    rv = 1;
  return rv;
}

/* THREAD_PROCF */
static int outputArgBlock(ErrMsg *errmsgp, 
#ifdef THREADS_DEBUG
			  uint64_t *readno,
#endif
			  void *thargp, 
			  void *buffargp)
{
  int errcode = ERRCODE_SUCCESS;
  short i;
  SmaltOutput *dop = (SmaltOutput *) thargp;
  SmaltArgBlock *blockp = (SmaltArgBlock *) buffargp;

  assert(buffargp != NULL);

  for (i=0; i<blockp->n_iobf && !(errcode); i++)
    errcode = outputIOBuffArg(errmsgp, 
			      dop,
			      blockp->iobfp+i);

#ifdef THREADS_DEBUG
  *readno = (i > 1)? blockp->iobfp->readno: 0;
  threadsPrintDebugMsg(NULL, NULL,
		       "outputArgBlock(): wrote read no %llu to %llu (%s)\n",
		       (long long unsigned int) *readno,
		       (long long unsigned int) (*readno) + i, 
		       (i == blockp->n_iobf)? "complete": "incomplete");
#endif 

  return errcode;
}


/*****************************************************************************
 ********************* Methods of Private Type SmaltMapArg *******************
 *****************************************************************************/
/* THREAD_INITF */
static int initMapArgs(void *ap, const void *ip, short threadno)
{
  SmaltMapArgs *map = (SmaltMapArgs *) ap;
  const SmaltMapConst *mcp = (const SmaltMapConst *) ip;
#ifdef hashhit_dump_sortarray
  char filnamstr[FILENAME_MAX];
#endif
  RMAPFLG_t rmapflg = mcp->rmapflg;
  if (mcp->inform == MENU_INFORM_SAM || mcp->inform == MENU_INFORM_BAM)
    rmapflg |= RMAPFLG_PAIRED;
  map->rmp = rmapCreate(mcp->htp, mcp->codecp, mcp->ssp, mcp->scormtxp, rmapflg);
  if (NULL == map->rmp)
    return ERRCODE_NOMEM;
  map->smconstp = mcp;
#ifdef hashhit_dump_sortarray
  sprintf(filnamstr, smalt_helper_oufilnam_fmt, threadno);
  if (NULL == (map->dumpfp = EFOPEN(filnamstr, "wb")))
    return ERRCODE_NOFILE;
  map->dumpctr = 0;
#endif
  map->threadno = threadno;
  map->errcode = ERRCODE_SUCCESS;

  return ERRCODE_SUCCESS;
}

static int cleanupMapArgs(ErrMsg *errmsgp, void *p)
{
  int errcode = ERRCODE_SUCCESS;
  if (NULL == p) {
    ERRMSGNO(errmsgp, ERRCODE_NULLPTR);
  } else {
    SmaltMapArgs *argp = (SmaltMapArgs *) p;
    rmapDelete(argp->rmp);
#ifdef hashhit_dump_sortarray
    errcode = EFCLOSE(argp->dumpfp);
    if (errcode)
       ERRMSGNO(errmsgp, errcode);
#endif
  }
  return errcode;
}

/* THREAD_PROCF */
static int processMapArgs(ErrMsg *errmsgp, 
			  SmaltMapArgs *map, 
			  SmaltIOBuffArg *brgp)
     /**< This routine and functions it calls do not fork threads */
{
  int errcode;
  uint8_t is_frac;
  uint32_t covermin_tuple;
  const ResultSet *rsltp;
  const SmaltMapConst *macop;
  SeqFastq *readp, *matep;

  assert(map != NULL);
  assert(brgp != NULL);

  macop = map->smconstp;
  assert(macop != NULL);

  readp = brgp->readp;
  matep = brgp->matep;

  ERRMSG_READNO(errmsgp, brgp->readno+1);
  ERRMSG_READNAM(errmsgp, seqFastqGetSeqName(readp));

  if (macop->tupcovmin < 0)
    return ERRCODE_ASSERT;

  if ((errcode = seqFastqEncode(readp, macop->codecp))) {
    ERRMSGNO(errmsgp, errcode);
    return errcode;
  }
  
  if (macop->tupcovmin < 1.01) {
    uint32_t readlen;
    seqFastqGetConstSequence(readp, &readlen, NULL);
    if (readlen > INT_MAX)
      ERRMSGNO(errmsgp, ERRCODE_OVERFLOW);
    covermin_tuple = (uint32_t) (macop->tupcovmin*readlen);
    if (covermin_tuple > readlen)
      covermin_tuple = readlen;
    is_frac = 1;
  } else {
    is_frac = 0;
    covermin_tuple = (uint32_t) macop->tupcovmin;
  }

  if (brgp->isPaired) {
    uint32_t covermin_tuple_mate;
    const ResultSet *rslt_matep;
    const ResultPairs *pairp;

    if ((errcode = seqFastqEncode(matep, macop->codecp))) {
      ERRMSGNO(errmsgp, errcode);
      return errcode;
    }

    if (is_frac) {
      uint32_t matelen;
      seqFastqGetConstSequence(matep, &matelen, NULL);
      if (matelen > INT_MAX) 
	ERRMSGNO(errmsgp, ERRCODE_OVERFLOW);
      covermin_tuple_mate = (uint32_t) (macop->tupcovmin*matelen);
      if (covermin_tuple_mate > matelen)
	covermin_tuple_mate = matelen;
    } else {
      covermin_tuple_mate = covermin_tuple;
    }

    rmapPair(errmsgp,
	     map->rmp, readp, matep,
#ifdef RESULTS_TRACKER
	     brgp->trkrp,
	     brgp->trkmp,
#endif
	     &brgp->pairflg,
	     macop->insert_min, macop->insert_max, 
	     macop->pairtyp,
	     macop->nhitmax_tuple,
	     (int) covermin_tuple, (int) covermin_tuple_mate,
	     macop->min_swatscor, macop->minbasq,
	     SMALT_TARGET_DEPTH, SMALT_MAX_DEPTH,
	     (RMAPFLG_t) (macop->rmapflg | RMAPFLG_PAIRED),
	     macop->scormtxp, macop->rfp, 
	     macop->htp, macop->ssp, macop->codecp);
   rmapGetData(&rsltp, 
		&rslt_matep, 
		&pairp,
		NULL, NULL, map->rmp);

    errcode = resultSetAddPairToReport(brgp->rep, 
				       macop->ihp, pairp, 
				       brgp->pairflg, macop->rsltouflg, 
				       rsltp,
				       rslt_matep);
    if ((errcode)) 
      ERRMSGNO(errmsgp, errcode);

    if (MENU_SAMPLE == macop->subprogtyp &&
	ERRCODE_SUCCESS == resultSetInferInsertSize(&brgp->isiz, RSLTSAMSPEC_V1P4, rsltp, rslt_matep))
      brgp->pairflg |= RSLTPAIRFLG_INSERTSIZ;

  } else {
    rmapSingle(errmsgp, 
	       map->rmp, 
#ifdef hashhit_dump_sortarray
	       map->dumpfp,
	       &map->dumpctr,
#endif
	       readp,
#ifdef RESULTS_TRACKER
	       brgp->trkrp,
#endif
	       macop->nhitmax_tuple,
	       covermin_tuple, (int) macop->min_swatscor, 
	       macop->swatscordiff,
	       macop->minbasq, 
	       SMALT_TARGET_DEPTH, SMALT_MAX_DEPTH,
	       (RMAPFLG_t)(macop->rmapflg & ~RMAPFLG_ALLPAIR), 
	       macop->scormtxp, macop->rfp,
	       macop->htp, macop->ssp, macop->codecp);
#ifdef hashhit_dump_sortarray
    if (map->dumpfp != NULL && map->dumpctr > 0) {
      EFCLOSE(map->dumpfp);
      map->dumpfp = NULL;
      errcode = ERRCODE_BREAK;
    }
#endif
    rmapGetData(&rsltp, 
		NULL, NULL,
		NULL, NULL, map->rmp);

    errcode = resultSetAddToReport(brgp->rep, macop->rsltouflg, rsltp);
    if ((errcode)) 
      ERRMSGNO(errmsgp, errcode);
  }
  return errcode;
}

static int processArgBlock(ErrMsg *errmsgp, 
#ifdef THREADS_DEBUG
			   uint64_t *readno,
#endif
			   void *targp, void *bufargp)
     /**< This routine and functions it calls do not fork threads */
{
  int errcode = ERRCODE_SUCCESS;
  short i;
  SmaltMapArgs *map = (SmaltMapArgs *) targp;
  SmaltArgBlock *blockp = (SmaltArgBlock *) bufargp;

  for (i=0; i<blockp->n_iobf && !(errcode); i++)
    errcode = processMapArgs(errmsgp, 
			     map, 
			     blockp->iobfp+i);
#ifdef THREADS_DEBUG
  *readno = (i > 1)? blockp->iobfp->readno: 0;
  threadsPrintDebugMsg(NULL, NULL,
		       "processArgBlock(): wrote read no %llu to %llu (%s)\n",
		       (long long unsigned int) *readno,
		       (long long unsigned int) (*readno) + i, 
		       (i == blockp->n_iobf)? "complete": "incomplete");

#endif
  return errcode;
}

/*****************************************************************************
 ******************* Methods for Histogram of Insert Sizes *******************
 *****************************************************************************/

static int prepSample(ErrMsg *errmsgp, const SmaltMapConst *common_args)
{
  int errcode = ERRCODE_SUCCESS;
  SEQNUM_t nreads = 0;
  SmaltInput *input_args = (SmaltInput *) threadsGetMem(THRTASK_INPUT);
  SmaltArgBlock *blockp = (SmaltArgBlock *) threadsGetMem(THRTASK_ARGBUF);
  SmaltOutput *dop = (SmaltOutput *) threadsGetMem(THRTASK_OUTPUT);
  SmaltIOBuffArg *iobfp;

  assert(common_args != NULL);
  assert(input_args != NULL);
  assert(blockp != NULL);
  assert(dop != NULL);

  iobfp = blockp->iobfp;
  assert(iobfp != NULL);
	 
  errcode = infmtCheckReads(input_args->ifrp, iobfp->readp, iobfp->matep,
			    &nreads, NULL, NULL, errmsgp);

  if ((errcode) && errcode != ERRCODE_RNAMPAIR)
    return errcode;

  if ((common_args->menuflg & MENUFLAG_VERBOSE)) {
    fprintf(stderr,"# Check of read pairs ok ...\n# Mate names %s ...\n", 
	   (errcode == ERRCODE_RNAMPAIR)? "don't match": "match");
  }
  infmtReset(input_args->ifrp);

  insSetSamplingInterval(dop->isamp, nreads, common_args->readskip);
  setInSampleIntervalInput(input_args, dop->isamp);

  return errcode;
}

static int outputHisto(ErrMsg *errmsgp)
{
  int errcode = ERRCODE_SUCCESS;
  SmaltOutput *dop = (SmaltOutput *) threadsGetMem(THRTASK_OUTPUT);
  if (dop->isamp) {
    InsHist *ihistp = insMakeHistoFromSample(dop->isamp);
    FILE *oufp = reportGetWriterStream(dop->writerp);
    if (!oufp) {
      errcode = ERRCODE_ASSERT;
      ERRMSGNO(errmsgp, errcode);
    }
    fprintf(oufp, "# Sampled histogram\n");
    insPrintHisto(oufp, SMALT_LINWIDTH, 0, ihistp);
    fprintf(oufp, "# Smoothed histogram\n");
    insPrintHisto(oufp, SMALT_LINWIDTH, 1, ihistp);
    if ((errcode = insWriteHisto(oufp, 0, ihistp)) &&
	(ERRCODE_FAILURE != errcode || NULL != ihistp))
      ERRMSGNO(errmsgp, errcode);
    insDeleteHisto(ihistp);
  }
  
  return errcode;
}

/*****************************************************************************
 ******************** Main routine for mapping reads *************************
 *****************************************************************************/

static int mapReads(ErrMsg *errmsgp, 
		    const MenuOpt *menup)
/**< with nthreads == 0 no additional thread is forked (single thread execution) 
 */
{
  int errcode = ERRCODE_SUCCESS;
  unsigned char is_verbose;
  short nthreads = menuGetNumberOfThreads(menup);
  SmaltMapConst common_args;
  THREAD_CHECKF *checkf = NULL;
  THREAD_CMPF *cmpf = NULL;
  int arg_fac = SMALT_THREAD_ARGFAC_NOSRT;
#ifdef SMALT_TIMING
  time_t time_start, time_setup, time_stop;
#endif

  if ((errcode = threadsInit())) {
    ERRMSGNO(errmsgp, errcode);
    return errcode;
  }
  if (nthreads < 0) 
    ERRMSGNO(errmsgp, ERRCODE_ASSERT);

  if (nthreads > SMALT_MAX_THREAD_NUM) 
    nthreads = SMALT_MAX_THREAD_NUM;

#ifdef SMALT_TIMING
  time(&time_start);
#endif
  
  if ((errcode = initMapConst(&common_args, menup)))
    ERRMSGNO(errmsgp, errcode);

#ifdef SMALT_TIMING
  time(&time_setup);
#endif

  errcode = threadsSetTask(THRTASK_ARGBUF, 0, 
			   initArgBlock, &common_args, 
			   NULL, cleanupArgBlock, 
			   NULL, NULL,
			   sizeof(SmaltArgBlock));
  if (errcode)
    ERRMSGNO(errmsgp, errcode);

  errcode = threadsSetTask(THRTASK_INPUT, (short) ((nthreads > 0)? 1:0), 
			   initInput, &common_args, 
			   loadArgBlock, cleanupInput, 
			   NULL, NULL,
			   sizeof(SmaltInput));
  if (errcode)
    ERRMSGNO(errmsgp, errcode);
  
  errcode = threadsSetTask(THRTASK_PROC, nthreads, 
			   initMapArgs, &common_args, 
			   processArgBlock, cleanupMapArgs,
			   NULL, NULL,
			   sizeof(SmaltMapArgs));
  if (errcode)
    ERRMSGNO(errmsgp, errcode);

  if (nthreads > 0 && (menuGetFlags(menup) & MENUFLAG_READORDER)) {
    checkf = checkArgBlockReadNo;
    cmpf = cmpArgBlockReadNo;
    arg_fac = SMALT_THREAD_ARGFAC_SRT;
  }
  errcode = threadsSetTask(THRTASK_OUTPUT, 0, 
			   initSmaltOutput, &common_args, 
			   outputArgBlock, cleanupSmaltOutput,
			   checkf, cmpf,
			   sizeof(SmaltOutput));
  
  if (errcode)
    ERRMSGNO(errmsgp, errcode);

  if ((errcode = threadsSetUp(arg_fac)))
    ERRMSGNO(errmsgp, errcode);
  is_verbose = (unsigned char) ((common_args.menuflg & MENUFLAG_VERBOSE) != 0);

  if ( MENU_SAMPLE == menuGetSubProgTyp(menup) ) {
    errcode = prepSample(errmsgp, &common_args);
    common_args.rmapflg |= (RMAPFLG_BEST | RMAPFLG_ALLPAIR);
    if (is_verbose) 
      fprintf(stderr, "# Sampling insert size distribution ...\n");
  } else {
    if (is_verbose) 
      fprintf(stderr, "# Processing query reads ...\n");
  }
  
  errcode = threadsRun(); 
  /* expect this to exit with ERRCODE_EOF */
  if ((errcode) && errcode != ERRCODE_EOF)
    ERRMSGNO(errmsgp, errcode);

  outputHisto(errmsgp);
  
  threadsCleanup();

#ifdef SMALT_TIMING
  time(&time_stop);
#endif
  cleanupMapConst(&common_args);

#ifdef SMALT_TIMING
  menuPrintWallClockTime(stderr, time_start, time_setup,
			 "Time spent setting up hash index");
  menuPrintWallClockTime(stderr, time_setup, time_stop,
			 "Time spent mapping reads");
#endif
  return errcode;
}

/*****************************************************************************
 ********************* Main routine for checking reads ***********************
 *****************************************************************************/

static int checkReads(ErrMsg *errmsgp, const MenuOpt *menup)
{
  int errcode = ERRCODE_SUCCESS;
  const char *filnamA, *filnamB;
  InFmtReader *ifrp;
  SeqFastq *readp = NULL, *matep = NULL;
  int nfil = menuGetFileNames(menup, &filnamA, &filnamB);
  SEQNUM_t n_reads = 0;

  if (nfil < 1 || nfil > 3)
    ERRMSGNO(errmsgp, ERRCODE_ASSERT);

  if (!(readp = seqFastqCreate(0, SEQTYP_UNKNOWN)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  if (filnamB != NULL) {
    if (!(matep = seqFastqCreate(0, SEQTYP_UNKNOWN)))
      ERRMSGNO(errmsgp, ERRCODE_NOMEM);
  }
  ifrp = infmtCreateReader(&errcode, filnamA, filnamB,
#ifdef HAVE_BAMBAMC
			   NULL,
#endif
			   INFMT_FASTQ);
  if (NULL != ifrp) {
    //fprintf(stderr, "# counting %s ...\n", (matep == NULL)? "reads": "read pairs");
    errcode = infmtCheckReads(ifrp, readp, matep, &n_reads, NULL, NULL, errmsgp);
    infmtDeleteReader(ifrp);
  }
  seqFastqDelete(matep);
  seqFastqDelete(readp);
    
  fprintf(stderr, "# checked %llu %s: ", 
	 (long long unsigned int) n_reads,
	 (matep == NULL)? "reads": "read pairs");
  if (errcode == ERRCODE_RNAMPAIR) {
    fprintf(stderr, "ok, but mate names don't match.\n");
  } else if ((errcode)) {
    fprintf(stderr, "failure.\n");
  } else {
    fprintf(stderr, "ok.\n");
  }

  return errcode;
}

/*****************************************************************************
 *********************************** MAIN ************************************
 *****************************************************************************/

int main(int argc, char *argv[])
{
  int errcode = ERRCODE_FAILURE;
  time_t time_start, time_stop;
  MenuOpt *menup = menuCreate();
  ErrMsg *errmsg = 0;
 
  ERRMSG_CREATE(errmsg);

  if (!menup) 
    ERRMSGNO(errmsg, ERRCODE_NOMEM);

  if (!(errcode = menuParseCommandLine(menup, argc, argv))) {
  
    switch (menuGetSubProgTyp(menup)) {
    case MENU_INDEX:
#ifdef smalt_debug
      menuPrint(stderr, menup);
#endif
      time(&time_start);
      buildHashIndex(errmsg, menup);
      time(&time_stop);
      menuPrintWallClockTime(stderr, time_start, time_stop,
			     NULL);
      break;
    case MENU_CHECK:
      errcode = checkReads(errmsg, menup);
      if (errcode == ERRCODE_RNAMPAIR)
	errcode = ERRCODE_SUCCESS;
      break;
    case MENU_SAMPLE: 
    case MENU_MAP:
#ifdef smalt_debug
      menuPrint(stderr, menup);
#endif
      time(&time_start);
      mapReads(errmsg, menup);
      time(&time_stop);
      menuPrintWallClockTime(stderr, time_start, time_stop, 
			     NULL);
      break;
    case MENU_HELP:
      break;
    case MENU_VERSION:
      break;
    default:
      ERRMSG(errmsg, "unknown subprogram", ERRCODE_COMMANDLINE);
      break;
    }
  }
  if ((errcode) && errcode != ERRCODE_MENUE)
    ERRMSGNO(errmsg, errcode);

  menuDelete(menup);
  ERRMSG_END(errmsg);
  return errcode;
}
