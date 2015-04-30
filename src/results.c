/** Processing and output of pairwise alignment results */

/*****************************************************************************
 *****************************************************************************
 *                                                                           *
 *  Copyright (C) 2010 - 2014 Genome Research Ltd.                           *
 *                                                                           *
 *  Author: Hannes Ponstingl (hp3@sanger.ac.uk)                              *
 *                                                                           *
 *  This file is part of SMALT.                                              *
 *                                                                           *
 *  SMALT is free software: you can redistribute it and/or modify it under   *
 *  the terms of the GNU General Public License as published by the Free     *
 *  Software Foundation, either version 3 of the License, or (at your        *
 *  option) any later version.                                               *
 *                                                                           *
 *  This program is distributed in the hope that it will be useful, but      *
 *  WITHOUT ANY WARRANTY; without even the implied warranty of               *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU         *
 *  General Public License for more details.                                 *
 *                                                                           *
 *  You should have received a copy of the GNU General Public License along  *
 *  with this program.  If not, see <http://www.gnu.org/licenses/>.          *
 *                                                                           *
 *****************************************************************************
 *****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <assert.h>

#include "elib.h"
#include "array.h"
#include "randef.h"
#include "sort.h"
#include "diffstr.h"
#include "report.h"
#include "results.h"

//#ifndef results_debug
#define NDEBUG
//#endif

#include <assert.h>

enum RESULT_CONST {
  DEFAULT_BLOCKSIZ_DIFFSTR = 4096,
  DEFAULT_BLOCKSIZ_DIFFSTRBUF = 256,
  DEFAULT_BLOCKSIZ_RESULTSET = 256,
  DEFAULT_BLOCKSIZ_SEGMSET = 32,
  MAPSCOR_MAX = 60,         /**< Maximum mapping score */
  MAPSCOR_DUMMY_COUNT = 3, /**< Dummy count for fraction in log arg */
  MAPSCOR_MAX_RANDOM = 3,  /**< minimum Smith-Waterman score for a unique mapping */
  MAPSCOR_MIN_UNIQ = MAPSCOR_MAX_RANDOM + 1, /**< minimum Smith-Waterman score for a unique mapping */ 
#ifdef results_mapscor_exp
  /* MAPSCOR_MAX_EQUALSW = 30, /\**< Maximum mapping score if Smith-Waterman scores are equal *\/ */
  //MAPSCOR_SCALFAC = 60,    /**< Mapping score scaling factor */
  MAPSCOR_EXPFAC  = 10,       /**< scales the exponent */
  //MAPSCOR_DEFICIT_SCALFAC = 2, /**< reduces the deficit by this factor */
#else
  //MAPSCOR_DUMMY_COUNT = 4,  /**< Dummy count for fraction in log arg */
  //MAXMAPSCOR_SCALFAC = 1,   /**< scaling of maximum mapping score */
  MAPSCOR_SCALFAC = 250,    /**< Mapping score scaling factor */
#endif
  MAPSCOR_THRESH_CONFIDENT = 20, /**< Mapping socore of a confident mapping (pair may be used to
  				  * derrive insert size) */
  MAPSCOR_MAX_TAG = 99,    /**< Maximum mapping score if score is output after the tag,
  			    * separeted by '::' from the tag */
  QUALSCOR_SCAL = 10, /**< scaling factor for quality scores
		       * Q = -QUALSCOR_SCAL*ln10(P) */

  RSLTX_INITVAL = -1,      /**< initialization value for index to primary results
  			    * in _RESULT.rsltx */
  QSEGX_INITVAL = -1,      /**< initialization value for index to complementary segment in query
			    * inf _RESULT.qsegx */
  SAMPLESIZ_MAPQ_RANDOM = 9,/**< Maximum smaple size for random draw up to which a mapping score
  			     * would have to be assigned */
  PAIRMAPSCOR_PROPER_INDEPENDENT = 12, /**< Add this mapping score to pairs that were mapped
  					* independently */
  PAIRMAPSCOR_PROPER_RESTRICTED = 6, /**< Add this mapping score if one mate was mapped in an area
  				      * restricted by the first */
  PAIRMAPSCOR_NBIT = 16,          /** number of bits for pair mapping score */
  PAIRMAPSCOR_MASK = 0x0000ffff,  /** mask out pair mapping score */
  N100PERCENT = 100,
  MIN_QSEGOVERLAP_PERCENT = 80, /**< Percentage of smaller semgent by which a pair segements
				 * of the query read have to overlap in order to be classed
				 * as 'overlapping' as opposed to 'complementary' */
};

enum RSLTSET_STATUS_FLAGS {    /**< Status of set of results */
  RSLTSETFLG_SEQX = 0x01,      /**< Sequence indices have been assigned */
  RSLTSETFLG_SERIALNO = 0x02,  /**< Results have beeb assigned serial numbers */
  RSLTSETFLG_SWSORT = 0x04,    /**< Results have been sorted by Smith-Waterman score,
				*   and have been assigned a rank */
  RSLTSETFLG_SEGIDX = 0x08,   /**< Segment indices have been assigned */
  RSLTSETFLG_MAPQ = 0x10,     /**< Mapping qualities have been assigned */
};

static const double MINLOGARG = 1E-7;
static const float QUALSCOR_LOGBASE = 2.30259;
/**< natural logarithm of 10, Q = -QUALSCOR_SCAL*ln(P)/QUALSCOR_LOGBASE */

typedef unsigned char UCHAR;
typedef int32_t SEGMLEN;
typedef unsigned char BOOL;
#ifdef results_debug
typedef uint32_t SEGIDX;
#endif

#ifdef RESULTS_BASQUAL
typedef struct _BASQ {
  uint64_t *correct;
  uint64_t *error;
  size_t siz;
} BASQ;
#endif

struct _RESULT { /**< Alignmment results */
  short serialno;    /**< Serial number (index) in ResultSet.resr */
#ifdef results_debug
  SEGIDX sgx;        /**< candidate segment index */
  char is_alidir;    /**< added from direct alignment without gaps */
#endif
  RSLTFLG_t status;   /**< combination of RESULT_STATUS_FLAGS */
  int swatscor;      /**< Smith-Waterman alignment score */ 
  int mapscor;       /**< 'Mapping' score */
  double prob;       /**< Mapping score as likelihood of being correct 
		      * (related to mapscor) */
#ifdef results_mscor_calib
  UCHAR mode;        /**< mode by which mapping score was calculated */
#endif
  SEQLEN_t q_start;   /**< Start in query sequence - original, not the potentially reverse complemented
		       * profile - counting from 1 */
  SEQLEN_t q_end;     /**< End in original query sequence, counting from 1 (q_start <= q_end) */
  SETSIZ_t s_start;  /**< if sidx>= 0: Start in subject sequence
		      *   else: Start in concatenated set of subject sequences (counting from 1) */
  SETSIZ_t s_end;    /**< if sidx >= 0: End in subject sequence, counting from 1 ( s_start <= s_end) 
		      *   else: End in concatenated set of subject sequences (counting from 1) */
  SEQNUM_t sidx;        /**< sequence index to which reference segment belongs, if <0 sequence
		      * index is unknown and s_start, s_end represent offset (counting from 0)
		      * in concatenated set of sequences */
  int stroffs;       /**< offset of alignment string. the aligment string is along the reference strand */
  int strlen;        /**< length of the alignment string (excl. termination) */
  short rsltx;       /**< Initialised to RSLTX_INITVAL<0: result is primary alignment. 
		      * >= 0: index of the secondary alignment in the array _ResultSet.resr */
  short qsegx;      /**< Index number of the segment in the set of complementary segments of
		     * the query this result belongs to (used for mapping quality and secondary
		     * alignment */
  short swrank;     /**< 0-based rank by Smith-Waterman score */

};

typedef Result* RESULTARR;     /**< Array of results */
typedef Result** RESULTPTRARR; /**< Array of pointers to results for sorting */

struct _ResultSet { /**< Set of alignment results */
  uint8_t status;     /**< One of RSLTSET_STATUS_FLAGS */
  RESULTARR resr;     /**< Array of alignment results */
  DiffStr *diffstrp;  /**< base pointer for alignment strings */
  int swatscor_max;   /**< maximum Smith-Waterman score in set */
  int swatscor_2ndmax;/**< 2nd largest Smith-Waterman score in set */
  RESULTPTRARR sortr; /**< for sorting overall */
  RESULTPTRARR segsrtr; /**< for sorting per segment */
  short *segnor;         /**< Array for designating the sorted subarray
		       * for segment i of pointers to elements of resr
		       * {segsrtr[segnor[i]], segsrtr[segnor[i]+1], ..., segsrtr[segnor[i+1]-1]}
		       */
  uint32_t *sortidxr;   /**< index for sorting resr */
  uint64_t *sortkeyr;  /**< key for sorting resr */
  DiffStr *diffstrbufp; /**< Buffer for diffstrings - this is used for explicit
			 * alignment output, alignment strings are along the reference strand */
  int n_ali_done;       /**< number of Smith-Waterman alignments executed */
  int n_ali_tot;        /**< total number of Smith-Waterman alignments that would
			 * have to be considered */ 
  short n_ali_max;      /**< Cut-off in the number of Smith-Waterman alignments */
  uint32_t n_hits_used;   /**< Number of seed hits used */
  uint32_t n_hits_tot;    /**< Total number of seed hits for read */
  short qsegno;           /**< number of complementary segments of the read max_i(rsr[i].qsegx) 
			   *  <0: not yet assigned */
#ifdef RESULTS_BASQUAL
  BASQ *basqp;            /**< Statistics on base qualities */
#endif
};

struct _ResultFilter {  /**< Output filter for results */
  int min_swscor;           /**< Minimum Smith-Waterman score */
  int min_swscor_below_max; /**< If >= 0 Filter out mappings with scores smaller than min_swscor_below_max
			     * below the maximum score. If < 0 don't apply this filter. */
  double min_identity;      /**< Minimum identity (number of matches as a fraction of alignment length */
};

struct _ResultBuffer {
  RESULTPTRARR bufr;
  short qsegno;
};

/******************************************************************************
 ************************************ Macros **********************************
 ******************************************************************************/
#define BLANK_RESULT(p) memset((p), 0, sizeof(Result))

#define SETDIFF(gap, typ) (gap) + (((unsigned char) (typ)) << DIFFSTR_TYPSHIFT)
/**< Set the one-character code (typ:gap) for the alignment segment */

#define TEST_RESULT_OVERLAP(r1p, r2p, min_overlap) \
  ((r1p)->q_start + (min_overlap) < (r2p)->q_end && (r2p)->q_start + (min_overlap) < (r1p)->q_end)
/******************************************************************************
 ******************************** Private Methods *****************************
 ******************************************************************************/
static int assignPhredScaledMappingScoreToRandomDraw(int samplesiz)
{
  int mapq;
  if (samplesiz < 1 || samplesiz > SAMPLESIZ_MAPQ_RANDOM) {
    mapq = 0;
  } else if (samplesiz == 1) {
    mapq = MAPSCOR_MAX_RANDOM + 1; /* signals non-random */
  } else {
    mapq = (int) (-QUALSCOR_SCAL*log10(((double)(samplesiz - 1))/samplesiz) + .499);
    if (mapq > MAPSCOR_MAX_RANDOM)
      mapq = MAPSCOR_MAX_RANDOM;
    else if (mapq < 0) {
      mapq = 0;
    }
  }
  return mapq;
}

static int sumQualOverMisMatch(int *qualsum_ali, BOOL with_nonali, 
			       const char *qualstrp, const SEQLEN_t slen, 
			       SEQLEN_t pos_start, SEQLEN_t pos_end, const DIFFSTR_T *dstrp)
     /**< \param spos Start position in sequence (starting fom 1) 
      */
{
  UCHAR q;
  SEQLEN_t spos;
  uint32_t qs;
  DIFFSTR_T gap, typ;
  const DIFFSTR_T *dp;
#define ADD_QVAL_CHECKED() { q = (UCHAR) qualstrp[spos];\
	                     if (q < SEQCOD_QVAL_OFFS) return ERRCODE_QUALVAL;\
	                     qs += (uint32_t) q - SEQCOD_QVAL_OFFS;\
	                     if (qs > INT_MAX) return ERRCODE_OVERFLOW;\
                            }

  if (pos_end < pos_start)
    return ERRCODE_ASSERT;

  //diffStrPrintf(stdout, dstrp,DIFFSTRFORM_RAW,0,0,0);
  //fprintf(stdout, "\n");
  *qualsum_ali = 0;
  qs = 0;

  spos = (pos_start > 0)? pos_start - 1: 0;
  for (dp = dstrp; (*dp); dp++) {
    DIFFSTR_GET((*dp), gap, typ);
    spos += gap;
    if (typ == DIFFCOD_D)
      continue;
    if (typ == DIFFCOD_S) {
      if (!(*(dp+1)))
	continue;
      if (spos < 1 || spos >= slen)
	return ERRCODE_ASSERT;
      ADD_QVAL_CHECKED();
    }
    spos++;     
  }

  if (spos != pos_end) 
    return ERRCODE_ASSERT;
  
  if (with_nonali) {
    for (spos=0; spos < pos_start-1; spos++) 
      ADD_QVAL_CHECKED();
    for (spos=pos_end; spos < slen; spos++) 
      ADD_QVAL_CHECKED();
  }

  *qualsum_ali = (int) qs;

  return ERRCODE_SUCCESS;
}

/******************************************************************************
 ******************************** Public Methods ******************************
 ******************************************************************************/

short resultConvertProbabilityToMappingScore(double p)
{
  short ms;
  double m;
  double isc = 1.0 - p;

  if (isc < MINLOGARG) isc = MINLOGARG;
  m = -QUALSCOR_SCAL * log10(isc);
  if (m > MAPSCOR_MAX) ms = MAPSCOR_MAX;
  else if (m < 0) ms = 0;
  else ms = (short) m;

  return ms;
}


#ifdef RESULTS_BASQUAL
/******************************************************************************
 ************************** Private Methods of Type BASQ **********************
 ******************************************************************************/

static void deleteBASQ(BASQ *p)
{
  if (p != NULL) {
    free(p->correct);
    free(p);
  }
}

static BASQ *createBASQ(void)
{
  size_t ucs = UCHAR_MAX;
  BASQ *p;
  EMALLOCP0(p);
  if (p) {
    ECALLOCP(2*ucs, p->correct);
    if (p->correct == NULL) {
      deleteBASQ(p);
      p = NULL;
    } else {
      p->error = p->correct + ucs;
      p->siz = ucs;
    }
  }
  return p;
}
static int sampleQualOverMisMatch(BASQ *basqp, 
				  const char *qualstrp, SEQLEN_t slen, 
				  SEQLEN_t pos_start, SEQLEN_t pos_end,
				  const DIFFSTR_T *dstrp)
{
  UCHAR g, q;
  SEQLEN_t spos;
  DIFFSTR_T gap, typ;
  const DIFFSTR_T *dp;

  if (pos_end < pos_start)
    return ERRCODE_ASSERT;
  
  spos = (pos_start > 0)? pos_start - 1: 0;
  for (dp = dstrp; (*dp); dp++) {
    DIFFSTR_GET((*dp), gap, typ);
    if (typ == DIFFCOD_M)
      gap++;
    for (g=0; g<gap && spos < pos_end; g++, spos++) {
      if (spos >= slen)
	return ERRCODE_ASSERT;
      q = (UCHAR) qualstrp[spos];
      if (q < SEQCOD_QVAL_OFFS)
	return ERRCODE_ASSERT;
      q -= SEQCOD_QVAL_OFFS;
      basqp->correct[q]++;
    }
    if (typ == DIFFCOD_D)
      continue;

    if (typ == DIFFCOD_S) {
      if (!(*(dp+1)))
	continue;
      if (spos >= slen)
	return ERRCODE_ASSERT;
      q = (UCHAR) qualstrp[spos];
      if (q < SEQCOD_QVAL_OFFS)
	return ERRCODE_ASSERT;
      q -= SEQCOD_QVAL_OFFS;
 
      basqp->error[q]++;
    }
    spos++;
  }

  return (spos == pos_end)? ERRCODE_SUCCESS: ERRCODE_ASSERT;
}

static int fprintBASQ(FILE *fp, const BASQ *basqA, const BASQ *basqB)
{
  UCHAR u;
  int q;
  uint64_t ntot, nerr;

  if (basqA->siz > UCHAR_MAX || 
      ((basqB) && basqA->siz != basqB->siz))
    return ERRCODE_ASSERT;

  for (u=0; u<basqA->siz; u++) {
    ntot = basqA->error[u] + basqA->correct[u];
    nerr = basqA->error[u];
    if (basqB) {
      ntot += basqB->error[u] + basqB->correct[u];
      nerr += basqB->error[u];
    }
    if (ntot < 1) {
      if (nerr != 0)
	return ERRCODE_ASSERT;
      q = 0;
    } else if (nerr < 1) {
      q = (int) UCHAR_MAX;
    } else {
      q = -QUALSCOR_SCAL*log(((double) nerr)/ntot)/QUALSCOR_LOGBASE;
    }
    
    fprintf(fp, "# BASQ [%u] %i\n", (unsigned int) u, q);
  }
  return ERRCODE_SUCCESS;
}
#endif /* def RESULTS_BASQUAL */


/******************************************************************************
 ************************* Private Methods of Type Result *********************
 ******************************************************************************/

#ifdef RESULTS_SUPERFLUOUS
static int cmpOffsRESULTp(const void *p1, const void *p2)
{
  const Result*ap = *((Result**) p1);
  const Result*bp = *((Result**) p2);

  /* compare match starts on subject sequence */
  if (ap->s_start < bp->s_start) return -1;
  if (ap->s_start > bp->s_start) return 1;

  return 0;
}

static int cmpOffsRESULT(const void *p1, const void *p2)
{
  const Result*ap = (Result*) p1;
  const Result*bp = (Result*) p2;

  /* compare match starts on subject sequence */
  if (ap->s_start < bp->s_start) return -1;
  if (ap->s_start > bp->s_start) return 1;

  return 0;
}
#endif

static int cmpRes(const void *p1, const void *p2)
     /**< For sorting results so that large segments are followed by segments
      * they include */
{
  const Result*ap = *((Result**) p1);
  const Result*bp = *((Result**) p2);
  SEQLEN_t da, db;
  /* compare subject sequence index */
  if (ap->sidx < bp->sidx) return -1;
  if (ap->sidx > bp->sidx) return 1;
  
  /* compare sense on subject sequence */
  if (((ap->status)&RSLTFLAG_REVERSE) < ((bp->status)&RSLTFLAG_REVERSE)) return -1;
  if (((ap->status)&RSLTFLAG_REVERSE) > ((bp->status)&RSLTFLAG_REVERSE)) return 1;

  /* compare match starts on subject sequence */
  if (ap->s_start < bp->s_start) return -1;
  if (ap->s_start > bp->s_start) return 1;
    
  /* compare segment lengths on query sequence */
  da = ap->q_end - ap->q_start;
  db = bp->s_end - bp->s_start;
  if (da > db) return -1; /* this is sorts shorter segments later */
  if (da < db) return 1;

  return 0;
}

static int cmpResOutput(const void *p1, const void *p2)
{
  const Result*ap = *((Result**) p1);
  const Result*bp = *((Result**) p2);
  SEQLEN_t da, db;

  /* by decreasing Smith-Waterman score */
  if (ap->swatscor > bp->swatscor) return -1;
  if (ap->swatscor < bp->swatscor) return 1;
  
  /* by strand of the reference sequence: forward 1st */
  if (((ap->status)&RSLTFLAG_REVERSE) < ((bp->status)&RSLTFLAG_REVERSE)) return -1;
  if (((ap->status)&RSLTFLAG_REVERSE) > ((bp->status)&RSLTFLAG_REVERSE)) return 1;

  /* by indreasing index (serial no) of reference sequence */
  if (ap->sidx < bp->sidx) return -1;
  if (ap->sidx > bp->sidx) return 1;

  /* by increasing alignment start points on reference sequence */
  if (ap->s_start < bp->s_start) return -1;
  if (ap->s_start > bp->s_start) return 1;
  
  /* by decreasing segment length */
  da = ap->q_end - ap->q_start;
  db = bp->q_end - bp->q_start;
  if (da > db) return -1;
  if (da < db) return 1;

  return 0;
}

static int cmpResSegSW(const void *p1, const void *p2)
{
  const Result*ap = *((Result**) p1);
  const Result*bp = *((Result**) p2);

  /* by increading serial no of the query segment */
  if (ap->qsegx < bp->qsegx) return -1;
  if (ap->qsegx > bp->qsegx) return 1;

  /* by decreasing Smith-Waterman score */
  if (ap->swatscor > bp->swatscor) return -1;
  if (ap->swatscor < bp->swatscor) return 1;
 
  return 0;
}

static int cmpResSegLen (const void *p1, const void *p2)
{
  const Result*ap = *((Result**) p1);
  const Result*bp = *((Result**) p2);
  SEQLEN_t da, db;

  /* by decreasing Smith-Waterman score */
  if (ap->swatscor > bp->swatscor) return -1;
  if (ap->swatscor < bp->swatscor) return 1;

  /* by decreasing segment lengths on query sequence */
  da = ap->q_end - ap->q_start;
  db = bp->q_end - bp->q_start;
  if (da > db) return -1; /* this is sorts shorter segments later */
  if (da < db) return 1;

  /* by strand of the reference sequence: forward 1st */
  if (((ap->status)&RSLTFLAG_REVERSE) < ((bp->status)&RSLTFLAG_REVERSE)) return -1;
  if (((ap->status)&RSLTFLAG_REVERSE) > ((bp->status)&RSLTFLAG_REVERSE)) return 1;

  /* compare subject sequence index */
  if (ap->sidx < bp->sidx) return -1;
  if (ap->sidx > bp->sidx) return 1;

  /* compare match starts on subject sequence */
  if (ap->s_start < bp->s_start) return -1;
  if (ap->s_start > bp->s_start) return 1;
  
  return 0;
}

static int isIdenticalResult(const Result*ap, const Result*bp)
{
  return (ap->s_start != bp->s_start ||
	  ap->s_end != bp->s_end ||
	  ap->q_start != bp->q_start ||
	  ap->q_end != bp->q_end ||
	  ap->swatscor != bp->swatscor ||
	  ap->sidx != bp->sidx)? 
    0: 1;
}

static int calcRESULTid(const Result*rp, const ResultSet *rsp)
{
  int matchnum = 0;
  diffStrCalcAliLen(&matchnum, rsp->diffstrp->dstrp + rp->stroffs);
  return matchnum;
}

#ifdef results_debug
static void fprintResultDebugInfo(FILE *fp, const Result*rp, const DiffStr *dfsp)
{
  const UCHAR *ucp;
  fprintf(fp, "result.c:RESULT_SWATSCORE %d\n", rp->swatscor);
  fprintf(fp, "result.c:RESULT_SEGMENT_QUERY   %d %d\n", rp->q_start-1, rp->q_end-1);
  fprintf(fp, "result.c:RESULT_SEGMENT_SUBJECT %d %d\n", (int) rp->s_start-1, (int) rp->s_end-1);
  fprintf(fp, "result.c:RESULT_DIFFS_RAW");
  for (ucp = dfsp->dstrp + rp->stroffs; *ucp; ucp++)
    fprintf(fp, " %u", (unsigned int) *ucp);
  fprintf(fp, "\n");
}

static int checkRESULT(SeqFastq *sqbufp, const Result*rp, const DiffStr *dfsp, 
		       const char *qstr, SEQLEN_t qlen,
		       const SeqSet *ssp, const SeqCodec *codecp)
{
  int errcode, s, i;
  char cod;
  const char *seqstr;
  SEQLEN_t seqlen, q, qend;
  const DIFFSTR_T *dp;
  DIFFSTR_T count, typ = 0;

  if (seqCodecType(codecp) != SEQCODTYP_3BITMANGLED)
    return ERRCODE_ASSERT;

  if (rp->q_start < 1 || rp->q_start > rp->q_end)
    return ERRCODE_ASSERT;

  if (rp->s_end < rp->s_start || rp->s_end - rp->s_start >=  INT_MAX)
    return ERRCODE_OVERFLOW;

  if (rp->sidx >= 0) {
    errcode = seqSetFetchSegmentBySequence(sqbufp, rp->sidx, 
					   rp->s_start - 1, 
					   rp->s_end - rp->s_start + 1,
					   ssp, codecp);
  } else {
    SETSIZ_t ss = rp->s_start - 1;
    SETSIZ_t se = rp->s_end - 1;
    errcode = seqSetFetchSegment(sqbufp, &ss, &se,
				 ssp, codecp);
  }
  if ((errcode))
    return errcode;

  seqstr = seqFastqGetConstSequence(sqbufp, &seqlen, &cod);
  if (seqlen != rp->s_end - rp->s_start + 1)
    return ERRCODE_ASSERT;

  if (cod != SEQCOD_MANGLED &&
      (errcode = seqFastqEncode(sqbufp, codecp)))
    return errcode;

  if (seqlen > INT_MAX)
    return ERRCODE_OVERFLOW;

  if (rp->status & RSLTFLAG_REVERSE) {
    q = qlen - rp->q_end;
    qend = qlen - rp->q_start;
  } else {
    q = rp->q_start - 1;
    qend = rp->q_end - 1;
  }
  s = 0;

  for (dp = dfsp->dstrp + rp->stroffs; (*dp); dp++) {
    DIFFSTR_GET((*dp), count, typ);
    if (typ == DIFFCOD_M)
      count++;
    for (i=0; i<count; i++) {
      if ((qstr[q++]&SEQCOD_ALPHA_MASK) != (seqstr[s++]&SEQCOD_ALPHA_MASK))
	return ERRCODE_FAILURE;
    }
    if (typ == DIFFCOD_S) {
      if ((*(dp+1)) && (qstr[q++]&SEQCOD_ALPHA_MASK) == (seqstr[s++]&SEQCOD_ALPHA_MASK))
	return ERRCODE_FAILURE;
    } else if (typ == DIFFCOD_I) {
      q++;
    } else if (typ == DIFFCOD_D) { /* DIFFCOD_D */
      s++;
    }
  }
  if (typ == DIFFCOD_M)
    return ERRCODE_DIFFSTR;
  if (q != qend + 1 || rp->s_end + 1 != rp->s_start + s) {
    diffStrPrintf(stdout, dfsp->dstrp + rp->stroffs, DIFFSTRFORM_RAW,0,0,0);
    errcode = ERRCODE_FAILURE;
  }

  return errcode;
}
#endif

static int sortBySegmentAndSWscor(ResultSet *rsp)
{
  short i, j, nres = (short) ARRLEN(rsp->sortr);
  
  if (nres < 1)
    return ERRCODE_SUCCESS;
    
  /** Allocate 2nd sort array */
  if (((size_t) nres) > ARRNALLOC(rsp->segsrtr)) {
    void *hp = ARREALLOC(rsp->segsrtr, nres);
    if (NULL == hp)
      return ERRCODE_NOMEM;
    rsp->segsrtr = hp;
  }
  if (((size_t) rsp->qsegno) >= ARRNALLOC(rsp->segnor)) {
    void *hp = ARREALLOC(rsp->segnor, rsp->qsegno+1);
    if (NULL == hp)
      return ERRCODE_NOMEM;
    rsp->segnor = hp;
  }
  for (i=0; i<nres; i++)
    rsp->segsrtr[i] = rsp->sortr[i];
  if (nres > 1)
    qsort(rsp->segsrtr, nres, sizeof(Result *), cmpResSegSW);
  ARRLEN(rsp->segsrtr) = nres;
  j=0;
  rsp->segnor[j++] = 0;
  for (i=1; i<nres; i++) {
    if (rsp->segsrtr[i]->qsegx < rsp->segsrtr[i-1]->qsegx)
      return ERRCODE_ASSERT;
    if (rsp->segsrtr[i]->qsegx > rsp->segsrtr[i-1]->qsegx)
      rsp->segnor[j++] = i;
  }
  rsp->segnor[j++] = nres;
  ARRLEN(rsp->segnor) = j;
  return (j == rsp->qsegno + 1)? ERRCODE_SUCCESS: ERRCODE_ASSERT;
}

static int labelComplementarySegments(ResultSet *rsp, 
				      short const min_overlap_percent)
/**< find complementary (on query sequence) segments, link overlapping ones */
{
  int errcode = ERRCODE_SUCCESS;
  RESULTPTRARR rspp = rsp->sortr;
  short i, i_start, n = (short) ARRLEN(rspp);
  double const min_overlap_frac = ((double) min_overlap_percent)/N100PERCENT;
  
  if (n < 1)
    return ERRCODE_SUCCESS;
  else if (n > 1 && !(rsp->status & RSLTSETFLG_SWSORT))
    return ERRCODE_ASSERT;

  for (i=0; i<n; i++) {
    rspp[i]->qsegx = QSEGX_INITVAL;
  }
  i_start = 0;
  rsp->qsegno = 0;

  do {
    Result*r1p = rspp[i_start];
    SEQLEN_t l1 = r1p->q_end - r1p->q_start;
    r1p->qsegx = rsp->qsegno;
    if (i_start >= SHRT_MAX)
      return ERRCODE_OVERFLOW;
    i = (short) (i_start + 1);
    i_start = 0;
    for (; i<n; i++) {
      Result*r2p = rspp[i];
      if (r2p->qsegx < 0) {
	SEQLEN_t l2 = r2p->q_end - r2p->q_start;
	SEQLEN_t min_overlap = (SEQLEN_t) (((l1 < l2)? l1: l2)*min_overlap_frac);
	if (TEST_RESULT_OVERLAP(r1p, r2p, min_overlap)) {
	  r2p->qsegx = rsp->qsegno;
	} else if (i_start == 0) {
	  i_start = i;
	}
      }
    }
    if (rsp->qsegno == SHRT_MAX)
      errcode = ERRCODE_OVERFLOW;
    rsp->qsegno++;
  } while (i_start != 0);

  if (!(errcode) &&
      !(errcode = sortBySegmentAndSWscor(rsp)))
    rsp->status |= RSLTSETFLG_SEGIDX;

  return errcode;
}

static int sortAndPrune(ResultSet *rsp)
     /**< sort results and remove duplicates. 
      * If ssp != Null, assign ref. sequence index and update reference offsets */
{
  short i, nres = (short) ARRLEN(rsp->resr);
  Result**dpp, **prevpp, **endpp;
#ifdef result_debug
  int skipctr=0;
#endif

  ARRLEN(rsp->sortr) = 0;
  for (i=0; i<nres; i++) {
    Result *rp = rsp->resr + i;
    rp->serialno = i;
    rp->swrank = 0;
    if ((rp->status&RSLTFLAG_SELECT)) {
      Result **sortp;
      ARRNEXTP(sortp, rsp->sortr);
      if (!sortp)
	return ERRCODE_NOMEM;
      *sortp = rp;
    }
  }
  rsp->status |= RSLTSETFLG_SERIALNO;

  nres = (short) ARRLEN(rsp->sortr);
  if (nres < 2) {
    rsp->status |= RSLTSETFLG_SWSORT;
    return ERRCODE_SUCCESS;
  }

  qsort(rsp->sortr, nres, sizeof(Result*), cmpRes);

  prevpp = rsp->sortr;
  endpp = rsp->sortr + nres;
  nres = 1;
  for (dpp=rsp->sortr+1; dpp<endpp; dpp++) {
    if ((*dpp)->s_end > (*prevpp)->s_end ||
	(*dpp)->swatscor > (*prevpp)->swatscor ||
	(*dpp)->q_start < (*prevpp)->q_start || 
	(*dpp)->q_end > (*prevpp)->q_end ||
	(*dpp)->sidx != (*prevpp)->sidx ||
	(((*dpp)->status)&RSLTFLAG_REVERSE) != (((*prevpp)->status)&RSLTFLAG_REVERSE)) {
      if (nres == SHRT_MAX)
	return ERRCODE_OVERFLOW;
      nres++;
      if ((++prevpp) < dpp) *prevpp = *dpp;
    } else {
      (*dpp)->status &= ~RSLTFLAG_SELECT;
#ifdef result_debug
      skipctr++;  
#endif
    }
  }
#ifdef result_debug
  printf("result.c:pruneAndSort skipped %i duplicates.\n", skipctr);
#endif

  qsort(rsp->sortr, nres, sizeof(Result*), cmpResOutput);
  ARRLEN(rsp->sortr) = (size_t) nres;

  /* now assign rank (by sw score) */
  if (nres > 0) {
    Result **rspp = rsp->sortr;
    rspp[0]->swrank = 0;
    for (i=1; i<nres; i++) {
      if (rspp[i]->swatscor > rspp[i-1]->swatscor)
	return ERRCODE_ASSERT;
      if (rspp[i]->swatscor < rspp[i-1]->swatscor) {
	rspp[i]->swrank = (short) (rspp[i-1]->swrank + 1);
      } else {
	rspp[i]->swrank = rspp[i-1]->swrank;
      }
    }
  }
  rsp->status |= RSLTSETFLG_SWSORT;

  return ERRCODE_SUCCESS;
}

static BOOL getNumberOfTopSwatRESULTs(short *n_best, RESULTPTRARR rspp)
/**< Get the number of results with best or 2nd best Smith-Waterman scores.
 * \return 0 if there are multiple best scores, 1 other wise (single best score).
 * \param n_best Returns the number of results with the best (0) or best/2nd_best (1)
 *        Smith-Waterman scores.
 * \param rssp Set of results.
 */
{
  BOOL rv;
  short nb, n = (short) ARRLEN(rspp);

  nb = n;
  if (n < 2 || rspp[1]->swatscor != rspp[0]->swatscor) {
    rv = 1;
  } else {
    rv = 0;
  }
 
  if (n > 2) {  
    int scorthresh = rspp[1]->swatscor;
    short i;
    for (i=2; i<n; i++)
      if (rspp[i]->swatscor != scorthresh)
	break;
    nb = i;
  }

  if (n_best) *n_best = nb;

  return rv;
}

/******************************************************************************
 ************************* Public Methods of Type Result **********************
 ******************************************************************************/

int resultGetData(SEQLEN_t *qs, SEQLEN_t *qe,
		  SEQLEN_t *rs, SEQLEN_t *re, SEQNUM_t *rx,
		  int *swscor, RSLTFLG_t *flag,
		  const Result *rp)
{
  int errcode;

  if (NULL == rp) {
    if (qs) *qs = 0;
    if (qe) *qe = 0;
    if (rs) *rs = 0;
    if (re) *re = 0;       
    if (rx) *rx = 0;
    if (swscor) *swscor = 0;
    if (flag) flag = 0;
    errcode = ERRCODE_NULLPTR;
  } else {
    if (qs) *qs = rp->q_start;
    if (qe) *qe = rp->q_end;  
    if (rs) *rs = (SEQLEN_t) rp->s_start;
    if (re) *re = (SEQLEN_t) rp->s_end;
    if (rx) *rx = rp->sidx;
    if (swscor) *swscor = rp->swatscor;
    if (flag) *flag = rp->status;
    errcode = ERRCODE_SUCCESS;
  }

  return errcode;
}

short resultGetFragmentNo(const Result *rp)
{
  return rp->qsegx;
}

short resultGetSWRank(const Result *rp)
{
  return rp->swrank;
}

RSLTFLG_t resultGetStatusFlag(const Result *rp)
{
  RSLTFLG_t flg = 0;
  if (rp != NULL)
    flg = rp->status;
  return flg;
}

int resultGetMapQualScore(double *prob, RSLTFLG_t *flag, const Result * const rp)
{
  int mapq = 0;
  if (rp != NULL) {
    mapq = (rp->mapscor < 0)? 0:rp->mapscor;
    if (prob) *prob = rp->prob;
    if (flag) *flag = rp->status;
  } else {
    if (prob) *prob = 0.0;
    if (flag) *flag = 0;
  }

  return mapq;
}

RSLTPAIRMAPFLG_t resultCalcInsertSize(int *isiz, 
				      unsigned char samspec,
				      const Result *ap,
				      const Result *bp)
{
  RSLTPAIRMAPFLG_t flag = 0; /* 0: (aF,bF), 1: (aR,bF), 2: (aF,bR), 3: (aR,bR) */

  if (ap->status&RSLTFLAG_REVERSE)
    flag |= RSLTPAIRMAPFLG_REVERSE_1st;
  if (bp->status&RSLTFLAG_REVERSE)
    flag |= RSLTPAIRMAPFLG_REVERSE_2nd;
  if (bp->s_start < ap->s_start)
    flag |= RSLTPAIRMAPFLG_LEFTMOST2nd;
  if (ap->sidx < 0 || bp->sidx < 0) {
    flag |= RSLTPAIRMAPFLG_NOCONTIG;
  } else if (ap->sidx == bp->sidx) {
    flag |= RSLTPAIRMAPFLG_SAMECONTIG;
  }

  if (isiz) {
    SETSIZ_t rA, rB;
    if (RSLTSAMSPEC_V1P4 == samspec) { 
      rA = (ap->s_start < bp->s_start)? ap->s_start: bp->s_start;
      rB = (ap->s_end < bp->s_end)? bp->s_end: ap->s_end;
      *isiz = (rA + INT_MAX > rB || rA < rB + INT_MAX)? rB - rA + 1: 0;
      if (flag & RSLTPAIRMAPFLG_LEFTMOST2nd)
	*isiz *= -1;
    } else {
      if (ap->status&RSLTFLAG_REVERSE) {
	rA = ap->s_end + ap->q_start;
      } else {
	rA = ap->s_start - ap->q_start + 1;
      }
      if (bp->status&RSLTFLAG_REVERSE) {    
	rB = bp->s_end + bp->q_start;
      } else {
	rB = bp->s_start - bp->q_start + 1;
      }
      *isiz = (rA + INT_MAX > rB || rA < rB + INT_MAX)? rB - rA: 0;
    }
  }

  return flag;
}

/* /\****************************************************************************** */
/*  ********************** Public Methods of Type ResultBuffer ******************* */
/*  ******************************************************************************\/ */
/* ResultBuffer *resultSetCreateBuffer(void) */
/* { */
/*   ResultBuffer *p; */
/*   EMALLOCP0(p); */
/*   if (p) { */
/*     ARRCREATE(p->bufr, 0); */
/*     if (NULL == p->bufr) { */
/*       resultSetDeleteBuffer(p); */
/*       p = NULL; */
/*     } */
/*   } */
/*   return p; */
/* } */

/* void resultSetDeleteBuffer(ResultBuffer *p) */
/* { */
/*   if (p) { */
/*     ARRDELETE(p->bufr); */
/*   } */
/*   free(p); */
/* } */


/******************************************************************************
 ************************ Private Methods of Type ResultSet *******************
 ******************************************************************************/

#define UPDATE_SWATSCORMAX(rsp,scor) if ((scor)>(rsp)->swatscor_2ndmax) {\
  if ((scor)>(rsp)->swatscor_max) {\
     (rsp)->swatscor_2ndmax = (rsp)->swatscor_max;\
     (rsp)->swatscor_max = (scor);}\
  else if ((scor) < (rsp)->swatscor_max) {\
     (rsp)->swatscor_2ndmax = (scor);}\
}
#ifdef RESULTS_SUPERFLUOUS
static int calcMappingScore(RESULTARR const *rspp)
     /**< Array of pointers into array of results, sorted by Smith-Waterman score
      */
{
  short i, n = ARRLEN(rspp);
  int swatscor_1st = rspp[0]->swatscor;
  int mapscor, swatscor_2nd = 0;
  
  if (swatscor_1st < 1)
    return 0;

  if (n < 2) {
    rspp[0]->mapscor = MAPSCOR_MAX;
    return rspp[0]->mapscor;
  }

  swatscor_2nd = rspp[1]->swatscor;
  for (i=2; i<n; i++) {
    if (rspp[i]->swatscor < swatscor_2nd)
      break;
  }
  /* number of results with 2nd largest score = i-1 */
/*   mapscor = swatscor_1st - swatscor_2nd - 2*(i-2); */
/*   if (mapscor < 0)  */
/*     mapscor = 0; */
  mapscor = (swatscor_1st - swatscor_2nd)*MAPSCOR_MAX/swatscor_1st;
  if (mapscor > MAPSCOR_MAX_TAG) 
    mapscor = MAPSCOR_MAX_TAG;
  
  rspp[0]->mapscor = mapscor;
  for(i=1;i<n;i++)
    rspp[i]->mapscor = 0;
  
  return mapscor;
}

static int calcMappingQuality(const ResultSet *rsetp, const SeqFastq *sqp)
     /**< Array of pointers into array of results, sorted by Smith-Waterman score
      */
{
  int errcode;
  short i;
  SEQLEN_t slen;
  const char *qualstrp;
  RESULTARR const *rspp = rsetp->sortr;
  short n =  ARRLEN(rspp);
  int qn;
  int swatscor_1st;
  int mapscor, swatscor_2nd = 0;
  int qvalsum[2];
  const DIFFSTR_T *dstrp = rsetp->diffstrp->dstrp;
#ifdef results_mscor_calib
  const char *namp;
#endif

  if (n < 1)
    return ERRCODE_SUCCESS;

  swatscor_1st = rspp[0]->swatscor;

  if (swatscor_1st < 1) {
    rspp[0]->mapscor = 0;
    return ERRCODE_SUCCESS;
  }

  if (n < 2) {
    mapscor = MAPSCOR_MAX;
  } else {
    swatscor_2nd = rspp[1]->swatscor;

    qualstrp = seqFastqGetConstQualityFactors(sqp, &slen, NULL);
    if (!qualstrp)
      return ERRCODE_ASSERT;

    if ((errcode = sumQualOverMisMatch(qvalsum, 1, qualstrp, slen, 
				       rspp[0]->q_start, rspp[0]->q_end, dstrp+rspp[0]->stroffs)))
      return errcode;
    if ((errcode = sumQualOverMisMatch(qvalsum+1, 1, qualstrp, slen,
				       rspp[1]->q_start, rspp[1]->q_end, dstrp+rspp[1]->stroffs)))
      return errcode;

    for (i=2; i<n; i++) {
      if (rspp[i]->swatscor < swatscor_2nd)
	break;
    }
    /* number of results with 2nd largest score = i-1 */
    qn = (int) QUALSCOR_SCAL*log((double) (i-1)) / QUALSCOR_LOGBASE;
    mapscor = qvalsum[1] - qvalsum[0] - qn;
    if (mapscor > MAPSCOR_MAX)
      mapscor = MAPSCOR_MAX;
    else if (mapscor < 0) 
      mapscor = 0;

#ifdef results_mscor_calib
    namp = seqFastqGetSeqName(sqp);
    printf("CALDAT_RESULTS %s %i %i %i %hi %i %i %i %i %u %u\n", 
	   namp, mapscor, 
	   qvalsum[1], qvalsum[0],
	   i-1, swatscor_1st, swatscor_2nd,
	   rsetp->n_ali_done, rsetp->n_ali_tot,
	   rsetp->n_hits_used, rsetp->n_hits_tot);
#endif
    
  }
  if (rsetp->n_ali_tot > rsetp->n_ali_done && rsetp->n_ali_tot > 1) {
    
    //const int qfracscor = - QUALSCOR_SCAL * log(1.0 - ((double)rsetp->n_ali_done)/rsetp->n_ali_tot) / QUALSCOR_LOGBASE;
/*     qfracscor = MAPSCOR_MAX * rsetp->n_ali_done / rsetp->n_ali_tot; */
/*     if (mapscor > qfracscor) */
/*       mapscor = qfracscor; */
    const double qfrac = ((double)rsetp->n_ali_done)/rsetp->n_ali_tot;
    mapscor *= qfrac;
  }
  
  rspp[0]->mapscor = (mapscor < MAPSCOR_MAX)? mapscor: MAPSCOR_MAX;
  for(i=1;i<n;i++)
    rspp[i]->mapscor = 0;
  
  return ERRCODE_SUCCESS;
}
#endif

static int calcPhredScaledMappingQuality(short qsegx, 
					 const ResultSet *rsetp, 
					 const SeqFastq *sqp)
     /**< Assigns PHRED scaled mapping score.
      * Degenerate (multiple best) mappings with 
      */
{
  int errcode;
  short i, i_min, n, n_swatscor_2nd = 0;
  SEQLEN_t slen, seglen, seglen_1st;
  const char *qualstrp;
  Result*tmpp;
  int qn, swatscor_1st, swatscor_2nd, mapscor;
  //int mapscor_deficit;
  int maxmapscor;
  double fs, fa;
  int qvalsum_1st=0, qvalsum_2nd=0, qvalsum_ali;
  const DIFFSTR_T *dstrp = rsetp->diffstrp->dstrp;
  RESULTPTRARR rspp;
#ifdef results_mscor_calib
  const char *namp;
  UCHAR mode = 0;
#endif
  //double seedfractlg, alifractlg;

  if (!(rsetp->status & RSLTSETFLG_SEGIDX) || qsegx < 0 || qsegx >= rsetp->qsegno)
    return ERRCODE_ASSERT;
  
  rspp = rsetp->segsrtr + rsetp->segnor[qsegx];
  n = (short) (rsetp->segnor[qsegx + 1] - rsetp->segnor[qsegx]);
  if (n < 1)
    return ERRCODE_SUCCESS;

  swatscor_1st = rspp[0]->swatscor;

  if (swatscor_1st < 1) {
    rspp[0]->mapscor = 0;
    return ERRCODE_SUCCESS;
  }

  
  /*   seedfractlg = (rsetp->n_hits_tot>0 && rsetp->n_hits_used > 0)?  */
  /*     log(((double) rsetp->n_hits_used)/rsetp->n_hits_tot): 0.0; */
  /*   alifractlg = (rsetp->n_ali_done > 0 && rsetp->n_ali_tot > 0)? */
  /*     log(((double) rsetp->n_ali_done)/rsetp->n_ali_tot): 0.0; */
 
  /*   mapscor_deficit = -QUALSCOR_SCAL*(seedfractlg + alifractlg)/QUALSCOR_LOGBASE; */
  /*   if (mapscor_deficit < 0) */
  /*     return ERRCODE_ASSERT; */

  fs = ((double)rsetp->n_hits_used)/(rsetp->n_hits_tot + MAPSCOR_DUMMY_COUNT);
  fa = ((double)rsetp->n_ali_done)/(rsetp->n_ali_tot + MAPSCOR_DUMMY_COUNT);
  if (fs > fa) fs = fa;
  fs = (fs > MINLOGARG)? -QUALSCOR_SCAL * log(fs)/QUALSCOR_LOGBASE: MAPSCOR_MAX;
  maxmapscor = (fs < MAPSCOR_MAX)? MAPSCOR_MAX - (int) fs: 0;

  /*   fs = (rsetp->n_hits_used > 0)?  */
  /*     -QUALSCOR_SCAL * log(((double)rsetp->n_hits_used)/(rsetp->n_hits_tot)/QUALSCOR_LOGBASE: */
  /*     MAPSCOR_MAX; */
  
  /*   fa = (rsetp->n_ali_done > 0)?  */
  /*     -QUALSCOR_SCAL * log(((double)rsetp->n_ali_done)/rsetp->n_ali_tot)/QUALSCOR_LOGBASE: */
  /*     MAPSCOR_MAX; */
  /*   if (fa > fs) fs = fa; */
  /*   if (fs > MAPSCOR_MAX) fs = (double) MAPSCOR_MAX; */
  /*   maxmapscor = MAPSCOR_MAX - fs; */
  /* maxmapscor = MAPSCOR_SCALFAC - MAPSCOR_DEFICIT_SCALFAC*mapscor_deficit; */
  /*   if (maxmapscor < 0)  */
  /*     maxmapscor = 0; */

  /*   if (n < 2) { */
  /*     /\* rspp[0]->mapscor = (maxmapscor< MAPSCOR_MAX)? maxmapscor: MAPSCOR_MAX; *\/ */
  /*     rspp[0]->mapscor = (mapscor_deficit<MAPSCOR_MAX)? MAPSCOR_MAX - mapscor_deficit: 0; */
  /*     return ERRCODE_SUCCESS; */
  /*   } */
  if (n>1) {
      swatscor_2nd = rspp[1]->swatscor;
      for (i=2; i<n && rspp[i]->swatscor == swatscor_2nd; i++); 
      n_swatscor_2nd = (short) (i-1); /* number of results with 2nd largest score = i-1 */
      qn = (int) (QUALSCOR_SCAL*log((double) (n_swatscor_2nd)) / QUALSCOR_LOGBASE);
    } else {
      swatscor_2nd = 0;
      n_swatscor_2nd = 0;
      qn=0;
    }
  if (swatscor_2nd == swatscor_1st && n>1) {
    /* multiple mappings with the same Smith-Waterman score -> select
     * the mapping with the longest segment in the query read. Assign
     * mapping score log(N-1/N) = log(N) - log(N-1)
     * If there are multiple mappings with maximum segment lengths
     * then select the mapping with the lowest sum of base qualities over
     * the mismatches. Assign mapping score qvalsum_2nd - qvalsum_1st - qn
     */
    qsort(rspp, n_swatscor_2nd + 1, sizeof(Result*), cmpResSegLen);
    seglen_1st = rspp[0]->q_end - rspp[0]->q_start;
    seglen     = rspp[1]->q_end - rspp[1]->q_start;
    if (seglen_1st == seglen) {
      qualstrp = seqFastqGetConstQualityFactors(sqp, &slen, NULL);
      if ((qualstrp)) {
	if ((errcode = sumQualOverMisMatch(&qvalsum_1st, 0, qualstrp, slen, 
					   rspp[0]->q_start, rspp[0]->q_end, dstrp+rspp[0]->stroffs)))
	  return errcode;

	if ((errcode = sumQualOverMisMatch(&qvalsum_2nd, 0, qualstrp, slen, 
					   rspp[1]->q_start, rspp[1]->q_end, dstrp+rspp[1]->stroffs)))
	  return errcode;
      
	i_min = 1;
	for (i=2; i<n && rspp[i]->swatscor == swatscor_1st; i++) {
	  seglen = rspp[i]->q_end - rspp[i]->q_start;
	  if (seglen < seglen_1st)
	    break;
	  if ((errcode = sumQualOverMisMatch(&qvalsum_ali, 0, qualstrp, slen,
					     rspp[i]->q_start, rspp[i]->q_end, dstrp+rspp[i]->stroffs)))
	    return errcode;
	  if  (qvalsum_ali < qvalsum_2nd) {
	    qvalsum_2nd = qvalsum_ali;
	    i_min = i;
	  }
	}
	if (qvalsum_1st > qvalsum_2nd) {
	  tmpp = rspp[i_min];
	  rspp[i_min] = rspp[0];
	  rspp[0] = tmpp;
	  /* mapscor = qvalsum_1st - qvalsum_2nd - qn; */
	  mapscor = MAPSCOR_MIN_UNIQ;
#ifdef results_mscor_calib
	  mode = 1;
#endif
	} else {
	  mapscor =  (qvalsum_1st == qvalsum_2nd)? 0: MAPSCOR_MIN_UNIQ;
#ifdef results_mscor_calib
	  mode = 2;
#endif
	}
      } else { /* if ((qualstrp)) */
	mapscor = 0;
#ifdef results_mscor_calib
	mode = 2;
#endif	
      }
    } else { /* if (seglen_1st == seglen) */
      /* mapscor = (int) QUALSCOR_SCAL * (log(n_swatscor_2nd + 1) - log(n_swatscor_2nd))/ QUALSCOR_LOGBASE; */
      mapscor = MAPSCOR_MIN_UNIQ;
#ifdef results_mscor_calib
      mode = 3;
#endif
    }
    if (mapscor < 1)
      qsort(rspp, n_swatscor_2nd + 1, sizeof(Result*), cmpResOutput);
    /* if (mapscor + mapscor_deficit > MAPSCOR_MAX_EQUALSW) */
/*       mapscor = MAPSCOR_MAX_EQUALSW - mapscor_deficit; */
  } else { /* if (swatscor_2nd == swatscor_1st && n>1) */
    /* mapscor = (swatscor_1st - swatscor_2nd)*maxmapscor/swatscor_1st - qn; */
    /* mapscor = (swatscor_1st - swatscor_2nd)*MAPSCOR_SCALFAC/swatscor_1st; */

#ifdef results_mapscor_exp
    SEQLEN_t qlen;
    seqFastqGetConstSequence(sqp, &qlen, NULL);
    mapscor = (int) (MAPSCOR_MAX*(1-exp(((double)(swatscor_2nd - swatscor_1st))*
					MAPSCOR_EXPFAC/qlen)) - qn);
    //mapscor = maxmapscor*(1-exp(((double)(swatscor_2nd - swatscor_1st))*MAPSCOR_EXPFAC/qlen)) - qn;
    //mapscor = MAPSCOR_SCALFAC*(1-exp((swatscor_2nd - swatscor_1st)*MAPSCOR_EXPFAC/swatscor_1st));
    //mapscor -= mapscor_deficit/MAPSCOR_DEFICIT_SCALFAC; 
#else
    SEQLEN_t qlen;
    seqFastqGetConstSequence(sqp, &qlen, NULL);
    mapscor =  MAPSCOR_SCALFAC*swatscor_1st/qlen*(swatscor_1st - swatscor_2nd)/qlen - qn;
#endif
#ifdef results_loscor_capped
    if (mapscor > maxmapscor)
      mapscor = maxmapscor;
    else if (mapscor < 0)
      mapscor = 0;
#else
    //mapscor += MAPSCOR_MIN_UNIQ;
    if (mapscor >= 0)
      mapscor += MAPSCOR_MIN_UNIQ;
    if (mapscor > maxmapscor)
      mapscor = maxmapscor;
/*     else if (mapscor < MAPSCOR_MIN_UNIQ) */
/*       mapscor = MAPSCOR_MIN_UNIQ; */
#endif
    
#ifdef results_mscor_calib
    mode = 4;
#endif
  }
  if (mapscor > MAPSCOR_MAX)
    mapscor = MAPSCOR_MAX;
  else if (mapscor < 0) 
    mapscor = 0;
    
#ifdef results_mscor_calib
  rspp[0]->mode = mode;
  namp = seqFastqGetSeqName(sqp);
  printf("CALDAT_RESULTS %s %i %i %i %i %hi %i %i %i %i %i %u %u\n", 
	 namp, mapscor, (int) mode,
	 qvalsum_1st, qvalsum_2nd,
	 n_swatscor_2nd, qn, swatscor_1st, swatscor_2nd,
	 rsetp->n_ali_done, rsetp->n_ali_tot,
	 rsetp->n_hits_used, rsetp->n_hits_tot);
#endif

  rspp[0]->mapscor = mapscor;
  for(i=1;i<n;i++)
    rspp[i]->mapscor = 0;
  
  return ERRCODE_SUCCESS;
}

static int 
propagateMapQualAsProb
(short qsegx, const ResultSet *rsetp)
/* calcPhredScaledMappingQuality() must have been called beforehand! */
{
  short i, ns, nn, n1=0, n2=0;
  double p1=0.0, p2=0.0;
  RESULTPTRARR rspp;

  if (!(rsetp->status & RSLTSETFLG_SEGIDX) || qsegx < 0 || qsegx >= rsetp->qsegno)
    return ERRCODE_ASSERT;
  
  rspp = rsetp->segsrtr + rsetp->segnor[qsegx];
  nn = (short) (rsetp->segnor[qsegx + 1] - rsetp->segnor[qsegx]);
  if (nn < 1)
    return ERRCODE_SUCCESS;
  
  for (i=1; i<nn && rspp[i]->swatscor == rspp[0]->swatscor; i++);
  n1 = i;
  if (i < nn) {
    assert(rspp[i]->swatscor < rspp[0]->swatscor);
    for (++i;i<nn && rspp[i]->swatscor == rspp[n1]->swatscor; i++);
    n2 = (short) (i - n1);
  }

  /* if (n1 > 1 && rspp[0]->mapscor > MAPSCOR_MAX_RANDOM) { */
  /*   n2 = n1 - 1; */
  /*   n1 = 1; */
  /* } */

  if (n1 == 1) {
    int isc = rspp[0]->mapscor;
    if (isc < 0)
      isc = 0;
    p2 = exp(((double) (-QUALSCOR_LOGBASE * isc))/QUALSCOR_SCAL);
    p1 = 1.0 - p2;
    if (n2 > 1)
      p2 /= n2;
  } else if (n1 > 1) {
    p1 = 1.0/n1;
    p2 = p1;
  }

  /* assign */
  for (i=0; i<n1; i++)
    rspp[i]->prob = p1;
  if (n1 + n2 > SHRT_MAX)
    return ERRCODE_OVERFLOW;
  ns = (short) (n1 + n2);
  for (;i<ns; i++)
    rspp[i]->prob = p2;
  for (;i<nn;i++)
    rspp[i]->prob = 0.0;

  if (1 == n1 && 0 == n2)
    rspp[0]->status |= RSLTFLAG_SINGLE;

  return ERRCODE_SUCCESS;
}

static int 
calcPhredScaledMappingQualityPerQuerySegment
(ResultSet *rsetp, const SeqFastq *sqp)
{
  int errcode = ERRCODE_SUCCESS;
  short qsegx;

  if (!(rsetp->status & RSLTSETFLG_SEGIDX))
    return ERRCODE_ASSERT;

  for (qsegx=0; qsegx<rsetp->qsegno; qsegx++) {
    if ((errcode = calcPhredScaledMappingQuality(qsegx, rsetp, sqp)))
      return errcode;
    if ((errcode = propagateMapQualAsProb(qsegx, rsetp)))
      return errcode;
  }
  if (!errcode)
    rsetp->status |= RSLTSETFLG_MAPQ;

  return errcode;
}

static int findSplitReads(RESULTARR const *rspp)
     /**< Array of pointers into array of results, sorted by 
      * Smith-Waterman core,
      * rspp[i]->rsltidx must have been initialised to 0 
      * This variable is set to the index of the primary alignment
      * in the (sorted) pointer array into results
      */
{
  short i, j, n = (short) ARRLEN(rspp);
  int swatscor_1st = rspp[0]->swatscor;
  short n_split = 0;
  Result*ap, *bp;
 
  for (i=0; i<n; i++) {
    /* only top->scoring hits can be primary alignments */
    ap = rspp[i];
    if (ap->swatscor < swatscor_1st)
      break;
    for (j=(short)(i+1); j<n; j++) {
      bp = rspp[j];
      if (bp->rsltx >= 0)
	continue;
      /* check that there is no overlap between the alignments */
      if (ap->q_end < bp->q_start || 
	  ap->q_start > bp->q_end) {
	bp->rsltx = i;
	ap->status |= RSLTFLAG_HASSECOND;
	n_split++;
	break;
      }
    }
  }

  return n_split;
}

static int splitMultiSpan(RESULTARR *resr, const uint32_t residx,
			  DiffStr *diffstrbufp, DiffStr *diffstrp, 
			  SeqFastq *sqbufp, 
			  SEQNUM_t so, SEQNUM_t eo, 
#ifdef results_debug
			  const char *profiled_seqp,
			  const char *profiled_seqRCp,
#endif
			  const ScoreProfile *scpp,
			  const ScoreProfile *scpRCp,
			  const SeqSet *ssp, 
			  const SeqCodec *codecp)
/**< Split alignment results that span multiple refrence sequences and add them
 * to the array. 
 * \param resr Array of alignment results.
 * \param residx Index of the element in the array that is going to be split.
 * \param so Index (serial no, 0 based) of the first reference sequence in ssp that is part 
 *        of the alignment with the query.
 * \param eo Index (serial no, 0 based) of the first reference sequence after sequence so in ssp
 *        that is not part of the alignment with the query (eo == no sequences in ssp is allowed).
 * \param scpp Profile of the query sequence used for scoring the segment.
 * \param scpRCp Profie of the reverse complement of the query sequence.
 * \param ssp Set of referenc sequences.
 * \param codecp Sequence de/encoder.
 */
{
  int errcode, i, n, idx;
  char code;
  char *seqstr;
  BOOL isReverseComplement;
  int curr_start, curr_end, s_start, s_end, q_start, q_end;
  const SETSIZ_t *ofp;
  SEQNUM_t nseq = seqSetGetOffsets(ssp, &ofp);
  Result*rp;
  const ScoreProfile *scprofp;
  SEQLEN_t seqlen, profiled_seqlen, profiled_seqlenRC;
  scoreGetProfile(NULL, &profiled_seqlen, NULL, NULL, scpp);
  scoreGetProfile(NULL, &profiled_seqlenRC, NULL, NULL, scpRCp);
  if (profiled_seqlen != profiled_seqlenRC)
    return ERRCODE_ASSERT;
#ifdef results_debug
  const char *profdseqp = NULL;
#endif

  if (ARRLEN(*resr) <= residx ||
      so < 0 || eo <= so || eo > nseq ||
      (*resr)[residx].s_start <= ofp[so])
    return ERRCODE_ASSERT;
  
  rp = (*resr) + residx;
  isReverseComplement = (BOOL) ((rp->status & RSLTFLAG_REVERSE) != 0);
  if (isReverseComplement) {
      scprofp = scpRCp; 
#ifdef results_debug
    profdseqp = profiled_seqRCp;
#endif
  } else {
    scprofp = scpp;
#ifdef results_debug
    profdseqp = profiled_seqp;
    if ((errcode = checkRESULT(sqbufp, rp, 
			       diffstrp, 
			       profdseqp, profiled_seqlen,
			       ssp, codecp)))
      return errcode;
#endif
  }

  n = eo - so; /* number of reference squences spanned by the alignment */
  for (i=0; i<n; i++) { 
    Result*hp;
    SEQLEN_t q_start_0based_inprofil;
    /* rp->s_start, rp->s_end start from 1,
     * ofp[i] start from 0. */
    idx = so + i;
    if (rp->s_start > ofp[idx]) {
      curr_start = 0;
    } else {
      curr_start = ofp[idx] - rp->s_start + 1;
    }
    curr_end = ((rp->s_end <= ofp[idx+1])? rp->s_end: ofp[idx+1]) - rp->s_start;
 
    errcode = diffStrSegment(diffstrbufp, 
			     diffstrp->dstrp + rp->stroffs,
			     curr_start, curr_end,
			     &s_start, &s_end,
			     &q_start, &q_end);
    if ( (errcode) ) {
      if (errcode == ERRCODE_NOMATCH)
	continue;
      else
	return errcode;
    }
#ifdef results_debug
    printf("results.c:splitMultiSpan(): curr_start = %i, curr_end = %i\n", 
	   curr_start, curr_end);
    diffStrPrintf(stdout, diffstrp->dstrp + rp->stroffs, DIFFSTRFORM_RAW,
		  0, 0, 0);
    printf("results.c:splitMultiSpan(): s_start = %i, s_end = %i\n", 
	   s_start, s_end);  
    printf("results.c:splitMultiSpan(): q_start = %i, q_end = %i\n", 
	   q_start, q_end);
    printf("results.c:splitMultiSpan(): new diffstr: ");
    diffStrPrintf(stdout, diffstrbufp->dstrp, DIFFSTRFORM_RAW,
		  0, 0, 0);

#endif
    ARRNEXTP(hp, *resr);
    if (!hp) 
      return ERRCODE_NOMEM;
    rp = (*resr) + residx;
    *hp = *rp;

    hp->stroffs = DIFFSTR_LENGTH(diffstrp);
    hp->strlen = DIFFSTR_LENGTH(diffstrbufp);
    if ((errcode = diffStrAppend(diffstrp, diffstrbufp)))
      return errcode;
#ifdef results_debug
    printf("results.c:splitMultiSpan(): appended diffstr: ");
    diffStrPrintf(stdout, diffstrp->dstrp+hp->stroffs, DIFFSTRFORM_RAW,
		  0, 0, 0);
#endif

    if (isReverseComplement) {
      hp->q_start = rp->q_end - q_end;
      hp->q_end = rp->q_end - q_start;
      q_start_0based_inprofil = profiled_seqlen - hp->q_end;
    } else {
      hp->q_start = rp->q_start + q_start;
      hp->q_end = rp->q_start + q_end;
      q_start_0based_inprofil = hp->q_start - 1;
    }
    if (hp->q_start > hp->q_end || hp->q_end > profiled_seqlen)
      return ERRCODE_ASSERT;

    hp->s_start = rp->s_start + s_start - ofp[idx];
    hp->s_end = rp->s_start + s_end - ofp[idx]; 
    if (hp->s_end < hp->s_start || hp->s_end - hp->s_start >= INT_MAX)
      return ERRCODE_OVERFLOW;

    hp->sidx = idx;
    hp->status &= ~RSLTFLAG_NOSEQID;
    hp->status |= RSLTFLAG_SELECT;

#ifdef results_debug
    if ((errcode = checkRESULT(sqbufp, hp, diffstrp, 
			       profdseqp, profiled_seqlen,
			       ssp, codecp)))
      return errcode;
#endif

    /* update Smith-Waterman score */
    if ((errcode = seqSetFetchSegmentBySequence(sqbufp, hp->sidx, 
						(SEQLEN_t) (hp->s_start - 1), 
						(SEQLEN_t) hp->s_end - hp->s_start + 1,
						ssp, codecp)))
      return errcode;
  
    seqstr = seqFastqGetSequence(sqbufp, &seqlen, &code);
    if (code != SEQCOD_MANGLED &&
	(errcode = seqFastqEncode(sqbufp, codecp)))
      return errcode;

    if (seqlen > INT_MAX)
      return ERRCODE_OVERFLOW;

    if ((errcode = aliScoreDiffStr(&hp->swatscor, seqstr, (int) seqlen, 
				   q_start_0based_inprofil,
				   diffstrp->dstrp + hp->stroffs, hp->strlen,
				   scprofp)))
      return errcode;
  }

  return ERRCODE_SUCCESS;
}

#ifdef RESULTS_BASQUAL
static int sampleBasQStats(ResultSet *rsetp, const SeqFastq *sqp)
{
  char cod;
  Result*rp;
  RESULTARR *rspp = rsetp->sortr;
  short n =  ARRLEN(rspp);
  SEQLEN_t slen;
  const char *basqstr = seqFastqGetConstQualityFactors(sqp, &slen, &cod);

  if (n < 1)
    return ERRCODE_SUCCESS;

  if (cod != SEQCOD_ASCII)
    return ERRCODE_SEQCODE;

  rp = rspp[0];
  
  return sampleQualOverMisMatch(rsetp->basqp, basqstr, slen, rp->q_start, rp->q_end,
				rsetp->diffstrp->dstrp + rp->stroffs);
}
#endif


#ifdef results_mapscor_pair
static double getMinProb(const RESULTPTRARR rsr, BOOL is_restricted)
{
  short nr = 1;
  double minprob = 1.0;

  if (!(is_restricted) || rsr[0]->mapscor < 1) {
    if (getNumberOfTopSwatRESULTs(&nr, rsr)) {
      /* single best score */
      minprob_read = exp((-QUALSCOR_LOGBASE * rsr[0]->mapscor)/QUALSCOR_SCAL);
      if (nr > 2) /*  multiple 2nd best scores */
	minprob /= nr-1;
    } else {
      if (nr < 2)
	return ERRCODE_ASSERT;
      minprob = 1.0/nr;
    }
  }
  
  return minprob;
}
#endif

static int assignSequenceIndex(ResultSet *rsp,
			       SeqFastq *sbufp,
#ifdef results_debug
			       const SeqFastq *sqp,
			       const SeqFastq *sqRCp,
#endif
			       const ScoreProfile *scpp,
			       const ScoreProfile *scpRCp,
			       const SeqSet *ssp,
			       const SeqCodec *codecp)
{
  int errcode = ERRCODE_SUCCESS;
  short i, ctr, nres = (short) ARRLEN(rsp->resr);
  const SETSIZ_t *ofp;
  SEQNUM_t s, e, nseq = seqSetGetOffsets(ssp, &ofp);

  ARRLEN(rsp->sortidxr) = 0;
  ARRLEN(rsp->sortkeyr) = 0;
  for (i=0; i<nres; i++) {
    Result *rp = rsp->resr + i;
    if ((rp->status&RSLTFLAG_SELECT) && rp->sidx < 0) {
      uint32_t *sip;
      uint64_t *skp;
      ARRCURR(rsp->sortidxr) = (uint32_t) i;
      ARRNEXTP(sip, rsp->sortidxr);
      ARRCURR(rsp->sortkeyr) = rp->s_start;
      ARRNEXTP(skp, rsp->sortkeyr);
      if (sip == NULL || skp == NULL)
	return ERRCODE_NOMEM;
    }
  }
  
  ctr = (short) ARRLEN(rsp->sortidxr);
  if (ctr > 1) {
    /* qsort(rsp->sortr, ctr, sizeof(Result*), cmpOffsRESULTp); */
    if ((errcode = sortUINT64andUINT32ArraysByQuickSort(ctr, rsp->sortkeyr, rsp->sortidxr)))
      return errcode;
  }
  
  /* assign sequence index */
  for (s=0,i=0; i<ctr && s<nseq; i++) {
    Result *rp = rsp->resr + rsp->sortidxr[i];
    if (rp->status&(RSLTFLAG_NOSEQID | RSLTFLAG_SELECT)) {
      for (;s < nseq && rp->s_start > ofp[s+1]; s++);
      for (e=s+1; e < nseq && rp->s_end > ofp[e]; e++);
      if (rp->s_end > ofp[e])
	return ERRCODE_ASSERT; /* not found */
      if (e > s + 1) {
#ifdef results_debug
	const char *profiled_seqp = seqFastqGetConstSequence(sqp, NULL, NULL);
	const char *profiled_seqRCp = seqFastqGetConstSequence(sqRCp, NULL, NULL);
#endif

	errcode = splitMultiSpan(&rsp->resr, rsp->sortidxr[i],
				 rsp->diffstrbufp, rsp->diffstrp, 
				 sbufp, s, e,
#ifdef results_debug
				 profiled_seqp,
				 profiled_seqRCp,
#endif
				 scpp, scpRCp,
				 ssp, codecp);
	rp = rsp->resr + rsp->sortidxr[i];
	rp->status &= ~RSLTFLAG_SELECT;
	if (errcode)
	  return errcode;
      } else {
	rp->sidx = s;
	rp->s_start -= ofp[s];
	rp->s_end -= ofp[s];
	rp->status &= ~RSLTFLAG_NOSEQID;
      }
    }
  }
  rsp->status &= ~RSLTSETFLG_SWSORT;
  rsp->status |= RSLTSETFLG_SEQX;

  ARRLEN(rsp->sortkeyr) = 0;
  ARRLEN(rsp->sortidxr) = 0;

  return errcode;
}

/******************************************************************************
 ************************ Public Methods of Type ResultSet ********************
 ******************************************************************************/
ResultSet *resultSetCreate(int blocksiz, int blocksiz_diffstr)
{
  ResultSet *rsp;

  EMALLOCP0(rsp);
  if (!rsp) return NULL;

  if (blocksiz < 1) blocksiz = DEFAULT_BLOCKSIZ_RESULTSET;
  if (blocksiz_diffstr < 1) blocksiz = DEFAULT_BLOCKSIZ_DIFFSTR;

  ARRCREATE(rsp->resr, blocksiz);
  ARRCREATE(rsp->sortr, blocksiz);
  ARRCREATE(rsp->segsrtr, blocksiz);
  ARRCREATE(rsp->sortidxr, blocksiz);
  ARRCREATE(rsp->sortkeyr, blocksiz);
  ARRCREATE(rsp->segnor, DEFAULT_BLOCKSIZ_SEGMSET);
  rsp->diffstrp = diffStrCreate(blocksiz_diffstr);
  rsp->diffstrbufp = diffStrCreate(DEFAULT_BLOCKSIZ_DIFFSTRBUF);
#ifdef RESULTS_BASQUAL
  rsp->basqp = createBASQ();
#endif
  if (!((rsp->resr) && (rsp->sortr) && 
	(rsp->segsrtr) && (rsp->segnor) && 
	(rsp->sortidxr) && (rsp->sortidxr) &&
	(rsp->diffstrp) && (rsp->diffstrbufp)
#ifdef RESULTS_BASQUAL
	&& (rsp->basqp)
#endif
	)) {
    resultSetDelete(rsp);
    rsp = 0;
  } else {
    rsp->status = 0;
  }

  return rsp;
}

void resultSetDelete(ResultSet *rsp)
{
  if (rsp) {
    ARRDELETE(rsp->resr);
    ARRDELETE(rsp->sortr);
    ARRDELETE(rsp->segsrtr);
    ARRDELETE(rsp->sortidxr);
    ARRDELETE(rsp->sortkeyr);
    ARRDELETE(rsp->segnor);
    diffStrDelete(rsp->diffstrp);
    diffStrDelete(rsp->diffstrbufp);
#ifdef RESULTS_BASQUAL
    deleteBASQ(rsp->basqp);
#endif
  }
  free(rsp);
}

void resultSetAlignmentStats(ResultSet *rsp, int n_ali_done, int n_ali_tot, short max_depth,
			     uint32_t n_hits_used, uint32_t n_hits_tot)
{
  rsp->n_ali_done = n_ali_done;
  rsp->n_ali_tot = n_ali_tot;
  rsp->n_ali_max = max_depth;
  rsp->n_hits_used = n_hits_used;
  rsp->n_hits_tot = n_hits_tot;
}

unsigned char resultSetAlignmentWasCurtailed(const ResultSet *rsp)
{
  return (unsigned char) (rsp->n_ali_tot > rsp->n_ali_max && rsp->n_ali_done >= rsp->n_ali_max);
}

int resultSetAddFromAli(ResultSet *rsp, const AliRsltSet *arsp, 
			SEQLEN_t soffs, 
			SEQLEN_t qoffs, SEQLEN_t qlen,
			SEQNUM_t seqidx, 
#ifdef results_debug
			SEGIDX sgx,
			SeqFastq *sqbufp,
			const char *profiled_seqp,
			const SeqCodec *codecp,
			const SeqSet *ssp,
#endif
			char is_reverse)
{
  int errcode = ERRCODE_SUCCESS;
  unsigned char is_new;
  short i, nres = aliRsltSetGetSize(arsp);
  int qs, qe, rs, re;
  Result*rp;
  const DiffStr *dfsp;

  if (nres < 1)
    return ERRCODE_SUCCESS;

  ARRNEXTP(rp, rsp->resr);
  if (!rp) return ERRCODE_NOMEM;

  BLANK_RESULT(rp);

  is_new = 0;
  rp->status = 0;
  rsp->status = 0;

  for (i=0; i<nres; i++) {
    if (is_new) {
      ARRNEXTP(rp, rsp->resr);
      if (!rp) return ERRCODE_NOMEM;
      is_new = 0;
      rp->status = 0;
    }
    
    errcode = aliRsltSetFetchData(arsp, i, &rp->swatscor,
				  &qs, &qe, &rs, &re,
				  &dfsp);
    if (errcode)
      break;

    if (is_reverse) {  
      rp->q_start = qoffs + qlen - qe; /* counting from 1 */
      rp->q_end = qoffs + qlen - qs;   /* counting from 1 */
    } else {
      rp->q_start = qs + qoffs + 1; /* counting from 1 */
      rp->q_end   = qe + qoffs + 1; /* counting from 1 */
    } 

    rp->s_start = soffs + rs + 1;
    rp->s_end   = soffs + re + 1;
    rp->sidx = seqidx;
    rp->swrank = 0;

    if (seqidx == RESULTSET_UNKNOWN_SEQIDX)
      rp->status |= RSLTFLAG_NOSEQID;
    is_new = (unsigned char) (ARRLEN(rsp->resr)<2 || !isIdenticalResult(rp, (rp-1)));
    if ((is_new)) {
      /* not a duplicate */
      /* copy alignment string */
      rp->stroffs = DIFFSTR_LENGTH(rsp->diffstrp);
      rp->strlen = DIFFSTR_LENGTH(dfsp);
      if ((errcode = diffStrAppend(rsp->diffstrp, dfsp)))
	return errcode;
     /* copy Smith-Waterman score and update maximum score */
      UPDATE_SWATSCORMAX(rsp, rp->swatscor);
      rp->status |= RSLTFLAG_SELECT;
      if ((is_reverse))
	rp->status |= RSLTFLAG_REVERSE;
      rp->mapscor = 0;
      rp->rsltx = RSLTX_INITVAL; /* primary alignment */
      rp->qsegx = QSEGX_INITVAL; /* complementary segment not yet established */
#ifdef results_debug
      rp->sgx = sgx;
      rp->is_alidir = 0;
      if ((errcode = checkRESULT(sqbufp, rp, rsp->diffstrp, 
				 profiled_seqp, qlen, ssp, codecp)))
	return errcode;
#endif
    } else { //if (is_new)
      if (--ARRLEN(rsp->resr) == 0)
	rsp->status = 0;
    }
  }
  return errcode;
}

int resultSetAddMisMatch(ResultSet *rsp, const int mmoffs[], int mmnum, 
			 SEQLEN_t qoffs, int len, SEQLEN_t soffs, SEQNUM_t sidx, char is_rcpl,
			 int swatscor)
{
  int errcode, diffstrlen;
  size_t newlen;
  int nl;
  Result*rp;

  if (len < 1) return ERRCODE_FAILURE;

  if ((errcode = diffStrGenerateFromMismatches(&diffstrlen, NULL, mmoffs, mmnum, len)))
    return errcode;

  newlen = rsp->diffstrp->len + diffstrlen + 1;
  if (newlen > INT_MAX)
    return ERRCODE_OVERFLOW;
  nl = (int) newlen;

  if (nl > rsp->diffstrp->n_alloc &&
      (errcode = diffStrRealloc(rsp->diffstrp, nl)))
    return errcode;

  ARRNEXTP(rp, rsp->resr);
  if (!(rp))
    return ERRCODE_NOMEM;
  BLANK_RESULT(rp);

  rp->status = RSLTFLAG_SELECT;
  if (is_rcpl) rp->status |= RSLTFLAG_REVERSE;
  rp->q_start = qoffs+1;
  rp->q_end = qoffs + len;
  rp->s_start = soffs+1;
  rp->s_end = soffs + len;
  rp->sidx = sidx;
  rp->swrank = 0;
  if (sidx == RESULTSET_UNKNOWN_SEQIDX) 
    rp->status |= RSLTFLAG_NOSEQID;
  rp->swatscor = swatscor;
  UPDATE_SWATSCORMAX(rsp, rp->swatscor);
  rp->mapscor = 0;
  rp->rsltx = RSLTX_INITVAL; /* primary alignment */
  rp->qsegx = QSEGX_INITVAL; /* complementary segment not yet established */
  rp->stroffs = rsp->diffstrp->len;
  rp->strlen = diffstrlen;

  if ((errcode = diffStrGenerateFromMismatches(&diffstrlen, rsp->diffstrp->dstrp+rp->stroffs, 
					       mmoffs, mmnum, len)))
    return errcode;
#ifdef results_debug
  printf("results_debug::resultSetAddMisMatch(): ");
  diffStrPrintf(stdout, rsp->diffstrp->dstrp + rp->stroffs, DIFFSTRFORM_RAW, 0, 0, 0);
#endif
  rsp->diffstrp->len += diffstrlen;
  rsp->status = 0;

  return ERRCODE_SUCCESS;
} 

void resultSetUpdateFromSegment(ResultSet *rsp, short start_idx, short end_idx,
				SEQLEN_t soffs, SEQLEN_t qoffs, char is_reverse)
{
  short i;
  Result*rp;
  for (i=start_idx; i<end_idx; i++) {
    rp = rsp->resr + i;
    if (rp->status & RSLTFLAG_RAW) {
      rp->s_start += soffs;
      rp->s_end   += soffs;
      rp->q_start += qoffs;
      rp->q_end   += qoffs;
      if ((is_reverse))
	rp->status |= RSLTFLAG_REVERSE;
      rp->status &= ~RSLTFLAG_RAW;
    }
  }
}

int resultSetSortAndAssignSequence(ResultSet *rsp,
				   SeqFastq *sbufp,
				   BOOL search_split,
				   const SeqFastq *sqp,
#ifdef results_debug
				   const SeqFastq *sqRCp,
#endif
				   const ScoreProfile *scpp,
				   const ScoreProfile *scpRCp,
				   const SeqSet *ssp,
				   const SeqCodec *codecp)
{
  int errcode = assignSequenceIndex(rsp, sbufp, 
#ifdef results_debug
				    sqp,
				    sqRCp,
#endif
				    scpp, scpRCp,
				    ssp, codecp);
  if (errcode)
    return errcode;

  if ((errcode = sortAndPrune(rsp)))
    return errcode;

  rsp->qsegno = 0;
  
  if (ARRLEN(rsp->sortr) > 0) {
    if ((errcode = labelComplementarySegments(rsp, MIN_QSEGOVERLAP_PERCENT)))
      return errcode;
   
    if ((errcode = calcPhredScaledMappingQualityPerQuerySegment(rsp, sqp)))
      return errcode;

#ifdef RESULTS_BASQUAL    
    if ((errcode = sampleBasQStats(rsp, sqp)))
      return errcode;
#endif
    if ((search_split))
      findSplitReads(rsp->sortr);
  }
  return ERRCODE_SUCCESS;
}

void resultSetBlank(ResultSet *rsp)
{
  ARRLEN(rsp->resr) = 0;
  ARRLEN(rsp->sortr) = 0;
  rsp->diffstrp->len = 0;
  rsp->swatscor_max = 0;
  rsp->swatscor_2ndmax = 0;
  rsp->n_ali_done = 0;
  rsp->n_ali_tot = 0;
  rsp->n_hits_used = 0;
  rsp->n_hits_tot = 0;
  rsp->qsegno = 0;
  rsp->status = 0;
}

int resultSetGetNumberOfSegments(short *nres, short *nseg, const ResultSet *rsp)
{
  int errcode = ERRCODE_SUCCESS;
  if (rsp == NULL) {
    errcode = ERRCODE_NULLPTR;
  } else if (ARRLEN(rsp->sortr) > SHRT_MAX || 
	     ARRLEN(rsp->segnor) > SHRT_MAX) {
    errcode = ERRCODE_OVERFLOW;
  } else if (ARRLEN(rsp->sortr) < 1) {
    if (nres != NULL) *nres = (short) 0;
    if (nseg != NULL) *nseg = (short) 0;
  } else if (ARRLEN(rsp->segnor) < 2) {
    errcode = ERRCODE_ASSERT;
  } else {
    if (nres != NULL) *nres = (short) ARRLEN(rsp->sortr);
    if (nseg != NULL) {
      *nseg = (short) (ARRLEN(rsp->segnor) - 1);
      if (*nseg < 0) *nseg = 0;
    }
  }
  
  return errcode;
}

int resultSetGetNumberOfResultsInSegment(int segx, const ResultSet *rsp)
{
  int nres = 0;
  
  if (rsp != NULL && segx >= 0 && ((size_t) segx + 1) < ARRLEN(rsp->segnor)) {
    nres = rsp->segnor[segx+1] - rsp->segnor[segx];
    if (nres < 0)
      nres = 0;
  }
  
  return nres;
}

int resultSetGetResultInSegment(const Result **rpp, int segx, int resx, const ResultSet *rsp)
{
  int errcode = ERRCODE_SUCCESS;
  *rpp = NULL;

  if (NULL == rsp) {
    errcode = ERRCODE_NULLPTR;
  } else if (ARRLEN(rsp->sortr) < 1) {
    errcode = ERRCODE_FAILURE;
  } else if (!(rsp->status & RSLTSETFLG_SEGIDX)) {
    errcode = ERRCODE_ASSERT;
  } else if (resx < 0 || segx < 0 || ((size_t) segx+1) >= ARRLEN(rsp->segnor)) {
    errcode = ERRCODE_ARGRANGE;
  } else {
    int xs = rsp->segnor[segx];
    int xe = rsp->segnor[segx+1];
    if (xe < xs) {
      errcode = ERRCODE_ASSERT;
    } else if (resx >= xe - xs) {
      errcode = ERRCODE_ARGRANGE;
    } else {
      *rpp = rsp->segsrtr[xs + resx];
    }
  }

  return errcode;
}

int resultSetGetResultByRank(const Result **rpp, int rank, const ResultSet *rsp)
{
  int errcode = ERRCODE_SUCCESS;
  *rpp = NULL;

  if (NULL == rsp) {
    errcode = ERRCODE_NULLPTR;
  } else if (!(rsp->status & RSLTSETFLG_SWSORT)) {
    errcode = ERRCODE_ASSERT;
  } else if (rank < 0 || ((size_t) rank) >= ARRLEN(rsp->sortr)) {
    errcode = ERRCODE_ARGRANGE;
  } else {
    *rpp = rsp->sortr[rank];
  }
  return errcode;
}

int resultSetGetMaxSwat(const ResultSet *rsp, int *maxswat_2nd)
{
  if (maxswat_2nd) *maxswat_2nd = rsp->swatscor_2ndmax;
  return rsp->swatscor_max;
}

Result *resultSetGetResultBySWrank(const ResultSet *rsp, short rank)
{
  Result *rp = NULL;

  if (rsp != NULL && rank >=0 && ((size_t) rank) < ARRLEN(rsp->sortr))
    rp = rsp->sortr[rank];

  return rp;
}

int resultSetDo(void *argp, ResultSetCallBackf *cbf, 
		const ResultSet *rsp)
{
  int errcode;
  int errc = RSLTSCBRTN_OK;
  short s, nseg, nres;

  if ((errcode = resultSetGetNumberOfSegments(&nres, &nseg, rsp)))
    return errcode;

  if (nres < 1)
    return ERRCODE_SUCCESS;
  
  if (!((rsp->status & RSLTSETFLG_SWSORT) && (rsp->status & RSLTSETFLG_SEGIDX)))
    return ERRCODE_ASSERT; 
  /* must have been sorted by segment number and sw-score */
  
  for (s=0; s<nseg && errc != RSLTSCBRTN_STOP; s++) {
    short r = rsp->segnor[s];
    for (; r<rsp->segnor[s+1]; r++) {
      errc = (*cbf)(&errcode, argp, rsp->segsrtr[r]);
      if (errcode != ERRCODE_SUCCESS)
	return errcode;
      if (errc != RSLTSCBRTN_OK)
	break;
    }
  }
  return errcode;
}

int resultSetAddResultToReport(Report *rep,
			       int pairid,
			       short mapscor,
			       REPMATEFLG_t mateflg,
			       REPPAIRFLG_t pairflg,
			       int isize,
			       const Result*rp,
			       const ResultSet *rsp)
{
  int errcode;

  if (rp == NULL || (rp->status & RSLTFLAG_NOOUTPUT)) {
    errcode = reportAddMap(rep,
			   pairid,
			   0, 0, 
			   0, 0,
			   0, 0, 0, 
			   NULL, 0, 
			   0,
			   mateflg, pairflg);
  } else {
    DIFFSTR_T *dstrp = rsp->diffstrp->dstrp + rp->stroffs;
    mateflg |= REPMATEFLG_MAPPED;

    if (rp->status & RSLTFLAG_REVERSE)
      mateflg |= REPMATEFLG_REVERSE;

    errcode = reportAddMap(rep,
			   pairid,
			   rp->swatscor, 
			   (short) ((pairid<0)? rp->mapscor: mapscor), 
			   rp->q_start, rp->q_end,
			   rp->s_start, rp->s_end, rp->sidx,
			   dstrp, rp->strlen,
			   isize,
			   mateflg, pairflg);
  }
  return errcode;
}

int resultSetAdd2ndaryResultsToReport(Report *rep, 
				      REPMATEFLG_t mateflg,
				      RESULTOUTFLG_t rsltflg,
				      const ResultSet *rsp)
{
  int errcode = ERRCODE_SUCCESS;

  if (NULL != rsp) {
    /* add from remaining segments */
    short qsegx = 0;
    for (; qsegx<rsp->qsegno; qsegx++) {
      short r = rsp->segnor[qsegx];
      int swscor = 0;
      for (; r<rsp->segnor[qsegx+1]; r++) {
	Result *rp = rsp->segsrtr[r];
	if ((rp->status & (RSLTFLAG_NOOUTPUT)))
	  continue;
	if ((rp->status & RSLTFLAG_REPORTED) ||
	    (rp->swatscor < swscor && 
	     ((rsltflg & RESULTFLG_BEST) || 
	      (rp->status & RSLTFLAG_BELOWRELSW))))
	  break;
	if ((errcode = resultSetAddResultToReport(rep, -1, 0, mateflg, 0, 0, rp, rsp)))
	  return errcode;
	rp->status |=  RSLTFLAG_REPORTED;
	swscor = rp->swatscor;
      }
    }
  }
  
  return errcode;
}

int resultSetAddToReport(Report *rep,
			 RESULTOUTFLG_t rsltflg,
			 const ResultSet *rsp)
{
  int errcode = ERRCODE_SUCCESS;
  short i, nsort = (short) ARRLEN(rsp->sortr);
  Result*rp = (nsort < 1)? NULL: rsp->sortr[0];
  REPMATEFLG_t mateflg = 0;
  
  if (rp != NULL) {
    short ns = 0;
    BOOL is_single = getNumberOfTopSwatRESULTs(&ns, rsp->sortr);
    if (rp->mapscor == 0 && !(is_single) && ns > 1 && 
	(rsltflg&RESULTFLG_BEST) && !(rsltflg&RESULTFLG_SPLIT)) {
      mateflg |= REPMATEFLG_MULTI;
      if ((rsltflg&RESULTFLG_RANDSEL)) {
	short r = (short) (RANDRAW_UNIFORM_1()*ns);
	rp = rsp->sortr[r];
	if ((rp))
	  rp->mapscor = assignPhredScaledMappingScoreToRandomDraw(ns);
      } else if ((rsltflg & RESULTFLG_SINGLE)) {
	rp = NULL;
      }     
    }
  }
  if ((errcode = resultSetAddResultToReport(rep, -1, 
					    0, 
					    (REPMATEFLG_t) (mateflg | 
							    REPMATEFLG_PRIMARY), 
					    0, 0, rp, rsp)))
    return errcode;
  if (rp != NULL) 
    rp->status |= RSLTFLAG_REPORTED;

  
  if ((rsltflg & RESULTFLG_SINGLE) && !(rsltflg&RESULTFLG_SPLIT))
    return ERRCODE_SUCCESS;

  /* remaining  best alignments */
  for (i=1; i<nsort; i++) {
    rp = rsp->sortr[i];
    if ((rsltflg&RESULTFLG_BEST) && rp->swatscor < rsp->sortr[i-1]->swatscor)
      break;
    if (!(rp->status & (RSLTFLAG_NOOUTPUT | RSLTFLAG_BELOWRELSW))) {
      if ((errcode = resultSetAddResultToReport(rep, -1, 0, mateflg, 0, 0, rp, rsp)))
	return errcode;
      rp->status |= RSLTFLAG_REPORTED;
#ifdef results_debug
      fprintf(stdout, "results_debug: segment index sgx = %u, is_alidir = %i\n", 
	      rp->sgx, (int) rp->is_alidir);
#endif
    }
  }
  
  /* 2ndary results */
  if ((rsltflg&RESULTFLG_BEST) && (rsltflg&RESULTFLG_SPLIT)) {
    errcode = resultSetAdd2ndaryResultsToReport(rep, 
						(REPMATEFLG_t) (mateflg | 
								REPMATEFLG_PARTIAL), 
						rsltflg,
						rsp);
  }

  return errcode;
}
   
void resultSetPrintDebugInfo(FILE *fp, const ResultSet *rsp)
{
  short i, nres = (short) ARRLEN(rsp->resr);
  const UCHAR *ucp;
  const char *cp;
  const Result*rp;

  for (i=0; i<nres; i++) {
    rp = rsp->resr + i;
    fprintf(fp, "PAIR_START %d\n", i+1);
    fprintf(fp, "PAIR_SEGMENT_A %u %u\n", rp->q_start-1, rp->q_end-1);
    fprintf(fp, "PAIR_SEGMENT_B %llu %llu\n", 
	    (long long unsigned int) rp->s_start-1, 
	    (long long unsigned int) rp->s_end-1);
    fprintf(fp, "PAIR_DIFFS");
    for (ucp = rsp->diffstrp->dstrp + rp->stroffs; *ucp; ucp++)
      fprintf(fp, " %u", (unsigned int) *ucp);
    fprintf(fp, "\n");
    fprintf(fp, "PAIR_DIFFS_B");
    for (cp = (char *) rsp->diffstrp->dstrp + rp->stroffs; *cp; cp++)
      fprintf(fp, " %i", (int) *cp);
    fprintf(fp, "\n");
  }
  fprintf(fp, "PAIR_END\n");
}

int resultSetGetScorStats(const ResultSet *rsp, 
			    int *scor_max, short *num_max, 
			    int *scor_2ndmax, short *num_2ndmax)
{
  short i, j;
  short nsort = (short) ARRLEN(rsp->sortr);
  if ((num_max || num_2ndmax)) {
    for (i=0; i<nsort; i++)
      if (rsp->sortr[i]->swatscor < rsp->swatscor_max)
	break;
    if (num_max) *num_max = i;

    if (num_2ndmax) {
      for (j=i; j<nsort; j++)
	if (rsp->sortr[i]->swatscor < rsp->swatscor_2ndmax)
	  break;
      *num_2ndmax = (short) (j-i);
    }
  }
  if (scor_max) *scor_max = rsp->swatscor_max;
  if (scor_2ndmax) *scor_2ndmax = rsp->swatscor_2ndmax;

  return (int) ARRLEN(rsp->resr);
}

BOOL resultSetGetRankDepth(const ResultSet *rsp, short *depth, short *rank)
{
  short n_max, n_2ndmax;
  
  resultSetGetScorStats(rsp, NULL, &n_max, NULL, &n_2ndmax);
  
  if (n_max < 2) {
    if (depth) *depth = (short) (n_max + n_2ndmax);
    if (rank) *rank = 1;
  } else {
    if (depth) *depth = n_max;
    if (rank) *rank = 0;
  }
  
  return (BOOL) (n_max == 1);
}

int resultSetGetMappingScore(const ResultSet *rsp, int *swscor)
{
  if (ARRLEN(rsp->sortr) < 1)  {
    if (swscor) *swscor = 0;
    return 0;
  }
  if (swscor) *swscor = rsp->sortr[0]->swatscor;

  return rsp->sortr[0]->mapscor;
}

double resultSetGetMapQualAsProb(double *pp2, short *nn1, short *nn2, const ResultSet *rsp)
/**< 
 * \return Probability for mapping with best score.
 * \param p2 Returns probability for mapping with 2nd best score.
 * \param n1 Returns Number of alignments with best score.
 * \param n2 Returns Number of alignments with 2nd best score.
 * \param rsp Set of results.
 */
{
  short n1, n2;
  double p1 = 0.0, p2 = 0.0;

  resultSetGetScorStats(rsp, NULL, &n1, NULL, &n2);
  if (n1 == 1) {
    int isc = rsp->sortr[0]->mapscor;
    if (isc < 0)
      isc = 0;
    p2 = exp(((double) (-QUALSCOR_LOGBASE * isc))/QUALSCOR_SCAL);
    p1 = 1.0 - (p2);
    if (n2 > 1)
      p2 /= n2;
  } else if (n1 > 1) {
    p1 = 1.0/(n1);
    p2 = p1;
  }
 
  if (pp2) *pp2 = p2;
  if (nn1) *nn1 = n1;
  if (nn2) *nn2 = n2;

  return p1;
}

int resultSetInferInsertSize(int *isiz, unsigned char samspec, const ResultSet *rsrp, const ResultSet *rsmp)
{
  int errcode = ERRCODE_FAILURE;
  
  *isiz = 0;
  if (ARRLEN(rsrp->sortr) > 0 && ARRLEN(rsmp->sortr) > 0) {
    const Result*rp = rsrp->sortr[0];
    const Result*mp = rsmp->sortr[0];
    if (rp->mapscor >= MAPSCOR_THRESH_CONFIDENT &&
	mp->mapscor >= MAPSCOR_THRESH_CONFIDENT &&
	rp->sidx == rp->sidx) {
      if (rp->sidx < 0) {
	errcode = ERRCODE_ASSERT;
      } else {
	RSLTPAIRMAPFLG_t flg = resultCalcInsertSize(isiz, samspec, rp, mp);
	if (flg == RSLTPAIRMAPFLG_REVERSE_1st)
	  *isiz *= -1;
	errcode = ERRCODE_SUCCESS;
      }
    }
  }
  
  return errcode;
}

BOOL resultSetTestOverlap(const ResultSet *rs1p, int overlap_percent, const ResultSet *rs2p)
{
  BOOL isOverlap = 0;
  int i;
  int rl1 = (int) ARRLEN(rs1p->sortr);
  int rl2 = (int) ARRLEN(rs2p->sortr);
  SEQLEN_t minlen, seg2len, overlaplen;
  const Result *r1p, *r2p;

  if (rl1 < 1 || rl2 < 1) 
    return 0;

  r1p = rs1p->sortr[0];
  minlen = r1p->q_end - r1p->q_start + 1;

  
  for (i=0; i<rl2 && (!isOverlap); i++) {
    r2p = rs2p->sortr[i];
    if (r2p->swatscor < rs2p->swatscor_2ndmax)
      break;
    seg2len = r2p->q_end - r2p->q_start + 1;
    overlaplen = (seg2len < minlen)? seg2len: minlen;
    overlaplen =  (overlaplen * overlap_percent)/N100PERCENT;
    
    if (TEST_RESULT_OVERLAP(r1p, r2p, overlaplen))
      isOverlap = 1;
  }

  return isOverlap;
}

Result *resultSetGetTopResult(BOOL *is_multi, const BOOL is_randsel, const ResultSet *rsp)
{
  short ntop;
  BOOL is_single = getNumberOfTopSwatRESULTs(&ntop, rsp->sortr);
  Result *toprp = NULL;
  *is_multi = 0;
  
  if (ntop > 0) {
    if ((is_single)) {
     toprp = rsp->sortr[0];
     if (toprp->mapscor < 1) /* best SW but still mapping quality 0 */
       *is_multi = 1;
    } else {
      *is_multi = 1;
    }
    if ((*is_multi) && (is_randsel)) {
      short rsltx = (short) (RANDRAW_UNIFORM_1()*ntop);
      toprp = rsp->sortr[rsltx];
      toprp->mapscor = assignPhredScaledMappingScoreToRandomDraw(ntop);
    }
  }

  return toprp;
}

#ifdef RESULTS_BASQUAL
int resultSetPrintBasqStats(FILE *fp, const ResultSet *rsAp, const ResultSet *rsBp)
{
  return fprintBASQ(fp, rsAp->basqp, (rsBp == NULL)? NULL:rsBp->basqp);
}
#endif

/******************************************************************************
 ********************* Private Methods of Type ResultFilter *******************
 ******************************************************************************/

static SEQLEN_t getFilterIdForRead(const ResultFilter *rfp, const SeqFastq *sqp)
{
  SEQLEN_t rlen, id = 0;

  if (rfp) {
    if (rfp->min_identity <= 1.0 && (sqp)) {
      seqFastqGetConstSequence(sqp, &rlen, NULL);
      id = (SEQLEN_t) (rfp->min_identity * rlen);
    } else {
      id = (SEQLEN_t) rfp->min_identity;
    }
  } 
  return id;
}

/******************************************************************************
 ********************** Public Methods of Type ResultFilter *******************
 ******************************************************************************/

ResultFilter *resultSetCreateFilter(void)
{
  ResultFilter *p;
  EMALLOCP0(p);
  return p;
}

void resultSetDeleteFilter(ResultFilter *p)
{
  free(p);
}

void resultSetFilterData(ResultFilter *p, int sw_abs, int sw_rel, double id_abs)
{
  if (p) {
    p->min_swscor = sw_abs;
    p->min_swscor_below_max = sw_rel;
    p->min_identity = id_abs;
  }
}

int resultSetFilterResults(const ResultSet *rsp,
			   const ResultFilter *rsfp, 
			   const SeqFastq *sqp)
{
  short i, n = (short) ARRLEN(rsp->sortr);
  int maxswscor, minabsswscor, minrelswscor;
  Result*rp;
  SEQLEN_t idthresh = getFilterIdForRead(rsfp, sqp);
  int minid = 0;

  if (n < 1) 
    return ERRCODE_SUCCESS;

  if (idthresh > INT_MAX)
    return ERRCODE_OVERFLOW;
  minid = (int) idthresh;

  maxswscor = rsp->sortr[0]->swatscor;
  minabsswscor = rsfp->min_swscor;
  minrelswscor = 0;
  if (rsfp->min_swscor_below_max >= 0 &&
      minabsswscor + rsfp->min_swscor_below_max < maxswscor)
    minrelswscor = maxswscor - rsfp->min_swscor_below_max;

  for (i=0; i<n; i++) {
    rp = rsp->sortr[i];
    if (rp->swatscor < minabsswscor ||
	calcRESULTid(rp, rsp) < minid) {
      rp->status |= RSLTFLAG_NOOUTPUT;
    } else if (rp->swatscor < minrelswscor)
      rp->status |= RSLTFLAG_BELOWRELSW;
  }
    
  return ERRCODE_SUCCESS;
}
