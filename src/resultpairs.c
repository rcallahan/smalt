/** Alignment results for read mates from the same insert */

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
#include <limits.h>
#include <string.h>

#include "elib.h"
#include "array.h"
#include "randef.h"
#include "resultpairs.h"

enum RESULTPAIRS_CONST {
    DEFAULT_BLOCKSIZ_PAIRS = 256,
    MAXPAIRNUM = 8192,       /**< Maximum number of pairs */
};

enum MAP_FLAGS {  /**< flags for mapping categories */
  MAPFLG_PAIRED = 0x01,  /**< Read(s) was(were) mapped as mate(s) of a pair */
  MAPFLG_CONTIG = 0x02,  /**< Both mates mapped to the same reference sequence */
  MAPFLG_PROPER = 0x04,  /**< Both mates mapped in proper orientation */
  MAPFLG_WITHIN = 0x08,  /**< Both mates mapped with in the distance expected from the insert size */
  MAPFLG_PARTIAL = 0x10, /**< Partial secondary alignment of a read */
  MAPFLG_MULT1ST = 0x20,/**< Multliple possible pair/read placements of 1st mate -> Reported with 
			  * a random assignment */
  MAPFLG_MULT2ND = 0x40, /**< Multliple possible pair/read placements of 2nd mate. */
  MAPFLG_INVALID = 0x80, /**< Read (pair) is flagged out (should not be considered) */
};

enum RSLTPAIRS_STATFLAGS { /**< Status flags for type ResultPairs */
  RSLTPAIRSFLG_INTVAL = 0x01,  /**< Pair intervals are set up */
  RSLTPAIRSFLG_MATES = 0x02,   /**< Array of mate pairs is set up */
  RSLTPAIRSFLG_MAPFLG = 0x04,  /**< The flag MATEPAIR.mapflg is set int the array of mate pairs */
  RSLTPAIRSFLG_CEILING = 0x08, /**< Reached the maximum number of pairs */
  RSLTPAIRSFLG_NPROPER = 0x10, /**< The number of proper pairs is up-to-date */
  RSLTPAIRSFLG_SORTED = 0x20,  /**< The array of mate pairs is sorted */
  RSLTPAIRSFLG_INSLIM = 0x40,  /**< Limits of the insert size range have been set */
};

enum PAIR_CATEGORIES { /**< Categories by which pairs are sorted */
  PAIRCAT_DEFAULT = 0, /**< not on same contig */
  PAIRCAT_SAMECONTIG = 1, /**< on same contig */
  PAIRCAT_WITHINLIMIT = 2, /**< on same contig and within limits */
  PAIRCAT_PROPER = 3,      /**< on same contig, within limits, and in expected orientation */
};

typedef unsigned char UCHAR_t;
typedef unsigned char BOOL_t;
typedef short RSLTIDX_t;
typedef uint8_t MAPFLG_t;     /* combination of MAP_FLAG bit flags */
typedef uint8_t PAIRMAPFLG_t; /* Holds combination of PAIRMAP_FLAGS bit flags */
typedef uint8_t PAIRSTATUSFLG_t; /* Combination of RSLTPAIRS_STATFLAGS bit flags */
typedef Result** RESULTPTRARR;

typedef struct _MATEPAIR { /**< Specifies a mapped mate-pair */
  const Result *ap;  /**< mapping results of a read */
  const Result *bp;  /**< mapping results of its mate */
  int ins;           /**< insert size */
  RSLTPAIRMAPFLG_t flag; /**< combination of RESULT_PAIRMAP_FLAGS */
  MAPFLG_t mapflg;   /**< combination of MAP_FLAGS */
  int32_t mscor;     /**< used for scoring proper mate pairs */ 
  uint32_t combiscor; /**< used for sorting */
  UCHAR_t cls;      /**< One of PAIR_CATEGORIES */
  double pbf;     /**< buffer for normalising probabilities */
  /* The following is for linked lists - linked by read segment (and mate segment) */
  /* struct _MATEPAIR *nextAp; /\**< Next pair that has this read segment in common *\/ */
  /* struct _MATEPAIR *prevAp; /\**< Previous pair with this read segment *\/ */
  /* struct _MATEPAIR *nextBp; /\**< Next pair that has this mate segment in common *\/ */
  /* struct _MATEPAIR *prevBp; /\**< Previous pair with this mate segment *\/ */
} MATEPAIR;

typedef MATEPAIR * MATEPAIRARR; /**< array (array.h) of mate pairs */

typedef struct _OFFSIVAL { /**< Specifies a range of offsets around a mapping of a read
			    * in which the mate may be found. */
  const Result *rp; /**< RESULT that generated this interval */
  RSLTFLG_t status; /**< status flag of the mapping rx that generated this interval */
  SEQNUM_t sidx;     /**< Index of the reference sequence to which the read was mapped */
  SEQLEN_t lower;  /**< lower (start) value of the range of offsets on that sequence */
  SEQLEN_t upper;  /**< upper (end) value of the range of offsets on that sequence */
} OFFSIVAL;

typedef OFFSIVAL * OFFSIVALARR; /**< array (array.h) of distance intervals */

struct _ResultPairs { /**< Results for mate pairs */
  PAIRSTATUSFLG_t status; /**< Combination of RSLTPAIRS_STATFLAGS */
  OFFSIVALARR ivr;  /**< array of offset intervals */
  MATEPAIRARR mpr;  /**< array of mate pairs */
  int n_proper;     /**< Number of proper pairs (proper orientation and within limits) */
  int n_within;     /**< Number of pairs within insert limits but not necessarily in the
		     * proper orientation */
  int dmin;         /**< Lower limit of the insert size range that defines a proper pair,
		     * status &  RSLTPAIRSFLG_INSLIM != 0 if set */
  int dmax;         /**< Upper limit of the insert size range that defines a proper pair,
		     * status &  RSLTPAIRSFLG_INSLIM != 0 if set */
};

static const double MINLOGARG = 1E-7;
static const double CUMULPROB_PROPER_OUTSIDE = 3e-3; /* i.e. 3 stdev of normal -> 99.7% */
static const double CUMULPROB_IMPROPER = 1e-4;

/******************************************************************************
 *********************************** Macros ***********************************
 ******************************************************************************/

#define INIT_MATEPAIR(p) memset((p), 0, sizeof(MATEPAIR))

/******************************************************************************
 ******************************* Private Methods ******************************
 ******************************************************************************/

static MAPFLG_t testProperPair(int isize, RSLTPAIRMAPFLG_t iflag, 
			       int dmin, int dmax, 
			       RSLTPAIRLIB_t libcode)
/** Test wether pair maps within limits in 'proper' orientation and
 * set bits MAPFLG_WITHIN or MAPFLG_PROPER accordingly, iflag contain the code for the paired-read
 * library type (RSLTPAIR_LIB) */
{
  UCHAR_t mapflg = 0;

  if (isize < 0) {
    if (isize <= -dmin && isize >= -dmax) 
      mapflg |= MAPFLG_WITHIN;
    if (libcode == RSLTPAIRLIB_PAIREDALL) {
      mapflg |= MAPFLG_PROPER;
    } else if (libcode == RSLTPAIRLIB_PAIREDEND) {
      if ((iflag & RSLTPAIRMAPFLG_REVERSE_1st) && !(iflag & RSLTPAIRMAPFLG_REVERSE_2nd)
	  && (iflag & RSLTPAIRMAPFLG_LEFTMOST2nd)) 
	mapflg |= MAPFLG_PROPER;
    } else if (libcode == RSLTPAIRLIB_MATEPAIR) {
      if (!(iflag & RSLTPAIRMAPFLG_REVERSE_1st) && (iflag & RSLTPAIRMAPFLG_REVERSE_2nd) &&
	  (iflag & RSLTPAIRMAPFLG_LEFTMOST2nd))
	mapflg |= MAPFLG_PROPER;
    } else if (libcode == RSLTPAIRLIB_SAMESTRAND) {
      if ((iflag & RSLTPAIRMAPFLG_REVERSE_1st) && (iflag & RSLTPAIRMAPFLG_REVERSE_2nd) &&
	  (iflag & RSLTPAIRMAPFLG_LEFTMOST2nd)) 
	mapflg |= MAPFLG_PROPER;
    }
  } else {
    if (isize >= dmin && isize <= dmax) 
      mapflg |= MAPFLG_WITHIN;
    if (libcode == RSLTPAIRLIB_PAIREDALL) {
      mapflg |= MAPFLG_PROPER;
    } else if (libcode == RSLTPAIRLIB_PAIREDEND) {
      if (!(iflag & RSLTPAIRMAPFLG_REVERSE_1st) && (iflag & RSLTPAIRMAPFLG_REVERSE_2nd) &&
	  !(iflag & RSLTPAIRMAPFLG_LEFTMOST2nd)) 
	mapflg |= MAPFLG_PROPER;
    } else if (libcode == RSLTPAIRLIB_MATEPAIR) {
      if ((iflag & RSLTPAIRMAPFLG_REVERSE_1st) && !(iflag & RSLTPAIRMAPFLG_REVERSE_2nd) &&
	  !(iflag & RSLTPAIRMAPFLG_LEFTMOST2nd))
	mapflg |= MAPFLG_PROPER;
    } else if (libcode == RSLTPAIRLIB_SAMESTRAND) {
      if (!(iflag & RSLTPAIRMAPFLG_REVERSE_1st) && !(iflag & RSLTPAIRMAPFLG_REVERSE_2nd) &&
	  !(iflag & RSLTPAIRMAPFLG_LEFTMOST2nd))
	mapflg |= MAPFLG_PROPER;
    }
  }

  return mapflg;
}

/******************************************************************************
 ******************** Functions of type ResultSetCallBackf ********************
 ******************************************************************************/

struct SETUPOFFSIVALARG_ {
  short max_rank;
  SEQLEN_t dmin, dmax;
  OFFSIVALARR oivr;
};

static int setupOFFSIVALcbf(int *errcode, void *argp, const Result *rp)
{
  SEQLEN_t s_start, s_end, q_start, q_end;
  SEQNUM_t sidx;
  RSLTFLG_t status;
  SEQLEN_t r0;
  OFFSIVAL *ivp;
  struct SETUPOFFSIVALARG_ *p = (struct SETUPOFFSIVALARG_ *) argp;

  if (resultGetSWRank(rp) > p->max_rank)
    return RSLTSCBRTN_BREAK;

  *errcode = resultGetData(&q_start, &q_end,
			   &s_start, &s_end, &sidx,
			   NULL, &status, rp);

  if ((*errcode))
    return RSLTSCBRTN_STOP;

  if (status&RSLTFLAG_REVERSE) {
      r0 = s_end + q_start - 2;
  } else {
      r0 = s_start - q_start;
  }

  ARRNEXTP(ivp, p->oivr);
  if (!ivp) {
    *errcode = ERRCODE_NOMEM;
    return RSLTSCBRTN_STOP;
  }
  ivp->rp = rp;
  ivp->sidx = sidx;
  ivp->status = status;
  if (r0 >= p->dmax) {
    ivp->upper = r0 - p->dmin;
    ivp->lower = r0 - p->dmax;
  } else {
    ivp->upper = (r0 > p->dmin)? r0 - p->dmin: 0;
    ivp->lower = 0;
  }

  ARRNEXTP(ivp, p->oivr);
  if (!ivp) {
    *errcode = ERRCODE_NOMEM;
    return RSLTSCBRTN_STOP;
  }

  ivp->rp = rp;
  ivp->sidx = sidx;
  ivp->status = status;
  ivp->upper = r0 + p->dmax;
  ivp->lower = r0 + p->dmin;

  if (ivp->lower <= (ivp-1)->upper) {
    (ivp-1)->upper = ivp->upper;
    ARRLEN(p->oivr)--;
  }

  return RSLTSCBRTN_OK;
}

struct GETPROPERMATEPAIRARG_ {
  short ivalx, max_rank;
  int maxnum;
  int swscor_min;
  RSLTPAIRLIB_t pairlibcode;
  const OFFSIVAL *oivr;
  ResultPairs *pairp;
};

static int getProperMATEPAIRcbf(int *errcode, void *argp, const Result *rp)
{
  SEQLEN_t s_start, s_end, q_start, q_end;
  SEQNUM_t sidx;
  int nival, swscor;
  RSLTFLG_t status;
  struct GETPROPERMATEPAIRARG_ *p = (struct GETPROPERMATEPAIRARG_ *) argp;
  ResultPairs * const pairp = p->pairp;

  if (resultGetSWRank(rp) > p->max_rank)
    return RSLTSCBRTN_BREAK;

  *errcode = resultGetData(&q_start, &q_end,
			   &s_start, &s_end, &sidx,
			   &swscor, &status, rp);

  if ((*errcode))
    return RSLTSCBRTN_STOP;

  if (swscor < p->swscor_min)
    return RSLTSCBRTN_BREAK;
  
  nival = ARRLEN(p->oivr);
  if (nival > SHRT_MAX) {
    *errcode = ERRCODE_ASSERT;
    return RSLTSCBRTN_STOP;
  }

  if (p->ivalx >= nival)
    p->ivalx = 0;
  
  for(;p->ivalx < nival; p->ivalx++) {
    int isiz;
    SEQLEN_t r0;
    MATEPAIR *mp;
    const OFFSIVAL *ivp = p->oivr + p->ivalx;
    if (sidx < ivp->sidx) 
	  break;
    if (sidx > ivp->sidx)
      continue;
    if (status&RSLTFLAG_REVERSE) {
      if (((ivp->status)&RSLTFLAG_REVERSE))
	/* look for 'proper' pair, i.e. for mappings on reverse compl. strand */
	continue;

      r0 = s_end + q_start - 2;
    } else {
      if (!((ivp->status)&RSLTFLAG_REVERSE))
	continue;

      r0 = s_start - q_start;
    } 

    if (r0 > ivp->upper)
      continue;

    if (r0 < ivp->lower)
      break;

    if (ivp->rp == NULL) {
      *errcode = ERRCODE_ASSERT;
      return RSLTSCBRTN_STOP;
    }

    ARRNEXTP(mp, pairp->mpr);
    if (!mp) {
      *errcode = ERRCODE_NOMEM;
      return RSLTSCBRTN_STOP;
    }
    INIT_MATEPAIR(mp);
    mp->ap = ivp->rp;
    mp->bp = rp;  
    mp->flag = resultCalcInsertSize(&mp->ins, RSLTSAMSPEC_V1P4, ivp->rp, rp);
    mp->mapflg = testProperPair(mp->ins, mp->flag, 
				pairp->dmin, pairp->dmax, p->pairlibcode);
    mp->mapflg |= (MAPFLG_PAIRED | MAPFLG_CONTIG);
    mp->mscor = 0;
    isiz = (mp->ins < 0)? -1*mp->ins: mp->ins;
    if (isiz < pairp->dmin || isiz > pairp->dmax)
      ARRLEN(pairp->mpr)--;
    if (ARRLEN(pairp->mpr) >= (size_t) p->maxnum) {
      *errcode = ERRCODE_PAIRNUM;
      return RSLTSCBRTN_STOP;
    }
  }

  return RSLTSCBRTN_OK;
}

struct GETMATEPAIRARG_ {
  short max_rankA;
  short max_rankB;
  RSLTPAIRFLG_t pairflg;
  RSLTPAIRLIB_t pairlibcode;
  const Result *ap;
  const ResultSet *rsltAp;
  const ResultSet *rsltBp;
  ResultPairs *pairp;
};

static int getMATEPAIRcbf_InnerLoop(int *errcode, void *argp, const Result *rp)
{
  struct GETMATEPAIRARG_ *p = (struct GETMATEPAIRARG_ *) argp;
  ResultPairs *pairp = p->pairp;
  MATEPAIR *mp;

  if (resultGetSWRank(rp) > p->max_rankB)
    return RSLTSCBRTN_BREAK;

  ARRNEXTP(mp, pairp->mpr);
  if (!mp) {
    *errcode =  ERRCODE_NOMEM;
    return RSLTSCBRTN_STOP;
  }

  mp->ap = p->ap;
  mp->bp = rp;
  mp->flag = 0;
  mp->ins = 0;
  mp->mapflg = MAPFLG_PAIRED;
  mp->mscor = 0;

  mp->flag = resultCalcInsertSize(&mp->ins, RSLTSAMSPEC_V1P4, mp->ap, mp->bp);
  if ((mp->flag & RSLTPAIRMAPFLG_SAMECONTIG)) { /* on the same contig */
    mp->mapflg = (MAPFLG_t) (mp->mapflg | testProperPair(mp->ins, mp->flag, pairp->dmin, pairp->dmax, 
							 p->pairlibcode));
    if ((mp->mapflg & MAPFLG_WITHIN)) {
      pairp->n_within++;
      if ((mp->mapflg & MAPFLG_PROPER))
	pairp->n_proper++;
    } 
    mp->mapflg |= MAPFLG_CONTIG;
  } 
  
  if (ARRLEN(pairp->mpr) >= MAXPAIRNUM) {
    pairp->status |= RSLTPAIRSFLG_CEILING;
    return RSLTSCBRTN_STOP;
  }

  return RSLTSCBRTN_OK;
}

static int getMATEPAIRcbf_OuterLoop(int *errcode, void *argp, const Result *rp)
{
  struct GETMATEPAIRARG_ *p = (struct GETMATEPAIRARG_ *) argp;
  if (resultGetSWRank(rp) > p->max_rankA)
    return RSLTSCBRTN_BREAK;
  p->ap = rp;
  *errcode = resultSetDo(argp, getMATEPAIRcbf_InnerLoop, 
			 p->rsltBp);
  return (*errcode)? RSLTSCBRTN_STOP: RSLTSCBRTN_OK;
}


/******************************************************************************
 *********************** Private Methods of Type OFFSIVAL *********************
 ******************************************************************************/

static int cmpOFFSIVAL(const void *a, const void *b)
     /**< this has to be coordinated with cmpRes(). The sort by
      * sequence index has to be reversed with respect to cmpRes since
      * the two mates of 'proper' paired-end reads map to opposite
      * strands.
      */
{
  const OFFSIVAL *ap = (OFFSIVAL *) a;
  const OFFSIVAL *bp = (OFFSIVAL *) b;

  if (ap->sidx < bp->sidx) return -1;
  if (ap->sidx > bp->sidx) return 1;

  if (((ap->status)&RSLTFLAG_REVERSE) < ((bp->status)&RSLTFLAG_REVERSE)) return 1;
  if (((ap->status)&RSLTFLAG_REVERSE) > ((bp->status)&RSLTFLAG_REVERSE)) return -1;

  if (ap->lower < bp->lower) return -1;
  if (ap->lower > bp->lower) return 1;

  return 0;
}

static int generateOFFSIVAL(OFFSIVALARR *oivr, 
			    int d_min, int d_max, 
			    const ResultSet *rsp)
     /**< Collect position intervals for the top scoring results of
      * the largest fragment.
      * Array of results must have been sorted by rsr[i].swatscor.
      */
{
  int errcode;
  struct SETUPOFFSIVALARG_ arg;
  
  resultSetGetNumberOfResultsInSegment(0, rsp);
  ARRLEN(*oivr) = (size_t) 0;

  if (resultSetGetNumberOfResultsInSegment(0, rsp) < 1) 
    return ERRCODE_SUCCESS;

  arg.dmin = (SEQLEN_t) (d_min < 0)? 0: d_min;
  arg.dmax = (SEQLEN_t) (d_max < 0)? 0: d_max;
  if (arg.dmin > arg.dmax) 
    return ERRCODE_ASSERT;
  
  arg.oivr = *oivr;
  arg.max_rank = 0;
  if ((errcode = resultSetDo(&arg, setupOFFSIVALcbf, rsp)))
    return errcode;
  *oivr = arg.oivr;

  /* sort intervals by sequence and lower limit */
  qsort(*oivr, ARRLEN(*oivr), sizeof(OFFSIVAL), cmpOFFSIVAL);
  
  return ERRCODE_SUCCESS; 
}

#ifdef results_debug
static void fprintOFFSIVAL(FILE *fp, const OFFSIVALARR oivr)
{
  short i, nival = ARRLEN(oivr);
  const OFFSIVAL *ivp;
  
  for (i=0; i<nival; i++) {
    ivp = oivr + i;
    fprintf(fp, "[%hi] rx = %i, sgx =% i, sidx = %lli, %c, lower = %u, upper = %u\n", 
	    i, (int) ivp->rx, (int) ivp->sgx, (long long signed) ivp->sidx, 
	    (ivp->status&RSLTFLAG_REVERSE)? 'R':'F',
	    ivp->lower, ivp->upper);
  }
}
#endif
/******************************************************************************
 *********************** Private Methods of Type MATEPAIR *********************
 ******************************************************************************/
#ifdef RESULTPAIRS_SUPERFLUOUS
static int cmpMATEPAIRbyScoreDescending(const void *ap, const void *bp)
{
  int32_t a = ((MATEPAIR *) ap)->mscor;
  int32_t b = ((MATEPAIR *) bp)->mscor;
  if (a > b)
    return -1;
  if (a < b) return 1;
  return 0;
}

static int cmpMATEPAIRbyCategoryDescending(const void *ap, const void *bp)
{
  UCHAR_t a = ((MATEPAIR *) ap)->cls;
  UCHAR_t b = ((MATEPAIR *) bp)->cls;
  if (a > b)
    return -1;
  if (a < b) return 1;
  return 0;
}
#endif

static int cmpMATEPAIRbyProbDescending(const void *ap, const void *bp)
{
  double const a = ((MATEPAIR *) ap)->pbf;
  double const b = ((MATEPAIR *) bp)->pbf;
  if (a > b)
    return -1;
  if (a < b) return 1;
  return 0;
}

static int findProperMATEPAIR(ResultPairs * const pairp, 
			      const ResultSet *rsp, 
			      int swscor_min,
			      OFFSIVAL const * const oivr, 
			      RSLTPAIRLIB_t pairlibcode,
			      int maxnum)
     /**< This routine looks for 'proper' mate pairs, i.e. pairs that
      * map to opposite strands within the distance interval. */
{
  int errcode;
  RSLTIDX_t nres;
  struct GETPROPERMATEPAIRARG_ arg;

  ARRLEN(pairp->mpr) = 0;

  if ((errcode = resultSetGetNumberOfSegments(&nres, NULL, rsp)))
    return errcode;
  
  if (nres < 1)
    return ERRCODE_SUCCESS;

  if (maxnum < 1) maxnum = 1;

  if (swscor_min > resultSetGetMaxSwat(rsp, NULL))
    return ERRCODE_SUCCESS;

  arg.ivalx = 0;
  arg.max_rank = 0;
  arg.maxnum = maxnum;
  arg.swscor_min = swscor_min;
  arg.pairlibcode = pairlibcode;
  arg.oivr = oivr;
  arg.pairp = pairp;
  errcode = resultSetDo(&arg, getProperMATEPAIRcbf, rsp);

  return errcode;
}
#ifdef RESULTPAIRS_SUPERFLUOUS
static int findMatePairWithinRange(short *px, int *n_proper, int dmin, int dmax, int idxno,
				   RSLTPAIRLIB_t pairlibcode,
				   const MATEPAIRARR mpr)
     /**< Return number of pairs within range. 
      * \param px Return index >= 0 to a proper pair, if there is a single proper pair. 
      *           Return < 0 otherwise (can be NULL).
      * \param n_proper Returns the number of proper pairs within the range (can be NULL).
      * \param dmin minimum insert size (> 0)
      * \param dmax maximum insert size (> 0)
      * \param idxno If idxno >=0, return index for pair number idxno (used for random sel) 
      *               if there are multiple pairs 
      * \param pairlibcode Code fore read pair library (one of RSLTPAIR_LIB)
      * \param mpr Array of mate pairs.
      */
{
  int i, npair = ARRLEN(mpr);
  int n, npr, idx;
  MAPFLG_t bitflg;

  if (npair < 1)
    return 0;
  n = npr = 0;
  idx = -1;
  for (i=0; i<npair; i++) {
    bitflg = testProperPair(mpr[i].ins, mpr[i].flag, dmin, dmax, pairlibcode);
    if ((bitflg & MAPFLG_WITHIN)) {
      if ((bitflg & MAPFLG_PROPER)) {
	if (idx < 0 || idxno == npr) {
	  idx = i;
	}
	npr++;
      }
      n++;
    }
  }
  
  if (px) *px = (short) ((npr != 1 && idxno < 0)? -1: idx);
  if (n_proper) *n_proper = npr;

  return n;
}
#endif

#ifdef results_pairscor_stats
static int getNumberOfProperPairs(int *n_proper, const MATEPAIRARR mpr)
/**< Returns the number of pairs within range that have the same best 
 * pair score.
 *
 * \param n_proper Returns the number of proper pairs within the range (can be NULL).
 * \param mpr Array of mate pairs.
 *
 * \note Pair mapping score must have been assigned (assignScoresToPairs) and the
 *       pairs sorted (cmpMATEPAIRbyScoreDescending) before this routine gets called.
 */
{
  int i, npair = ARRLEN(mpr);
  int n, npr;
  int32_t maxscor;

  if (npair < 1)
    return 0;

  n = npr = 0;

  for (i=0; i<npair; i++) {
    if (mpr[i].mscor < maxscor)
      break;
    if (mpr[i].mapflg & MAPFLG_WITHIN) {
      if (mpr[i].mapflg & MAPFLG_PROPER)  {
	npr++;
      }
      n++;
    }
  }

  return n;
}
static int fprintMatePairStats(FILE *fp, const MATEPAIRARR mpr, 
			       const RESULTPTRARR rresr, const RESULTPTRARR mresr,
			       int dmin, int dmax, RSLTPAIRLIB_t pairlibcode,
			       const InsHist *ihp)
{
  int i, npair = ARRLEN(mpr);
  int n = 0;
  MAPFLG_t bitflg;
  int32_t count, totnum;
  RESULT *rrp, *mrp;
  for (i=0; i<npair; i++) {
    count = totnum = 0;
    bitflg = testProperPair(mpr[i].ins, mpr[i].flag, dmin, dmax, pairlibcode);
    if ((bitflg & MAPFLG_WITHIN) && (bitflg & MAPFLG_PROPER))
      count = insGetHistoCount(&totnum, mpr[i].ins, 1, ihp);
    rrp = rresr[mpr[i].ax];
    mrp = mresr[mpr[i].bx];
    fprintf(fp, "PAIR[%i] posA=(%i, %llu), posB=(%i, %llu), ins=%i, mapqA=%i, mapqB=%i," \
	    "totnum=%i, counts=%i, mscor=%i\n",
	    n,
	    (int) rrp->sidx,
	    (long long unsigned int) rrp->s_start, 
	    (int) mrp->sidx,
	    (long long unsigned int) mrp->s_start,
	    mpr[i].ins, rrp->mapscor, mrp->mapscor,
	    totnum, count, mpr[i].mscor);    
    n++;
  }

  return n;
}
#endif

/******************************************************************************
 *********************** Pivate Methods of Type ResultPairs *******************
 ******************************************************************************/
#ifdef RESULTPAIRS_SUPERFLUOUS
static int getPairNumbersAndSortByCategory(int *n_proper, int *n_within, const ResultPairs *p)
/**< Return the number of all pairs, n_proper and n_within can be NULL 
 * Assign PAIR_CATEGORIES and sort according to those (proper pairs first)
 */
{
  int i, npr=0, nin=0;
  int n_pairs = ARRLEN(p->mpr);
  MATEPAIR *mp;

  for (i=0; i<n_pairs; i++) {
    mp = p->mpr + i;
    mp->cls = PAIRCAT_DEFAULT;
    if (mp->mapflg & MAPFLG_CONTIG) {
      if (mp->mapflg & MAPFLG_WITHIN) {
	nin++;
	if (mp->mapflg & MAPFLG_PROPER) {
	  mp->cls = PAIRCAT_PROPER;
	  npr++;
	} else {
	  mp->cls = PAIRCAT_WITHINLIMIT;
	}
      } else {
	mp->cls = PAIRCAT_SAMECONTIG;
      }
    }
  }
  if (n_pairs > 1)
    qsort(p->mpr, n_pairs, sizeof(MATEPAIR), cmpMATEPAIRbyCategoryDescending);

  if (n_proper) *n_proper = npr;
  if (n_within) *n_within = nin;

  return n_pairs;
}
#endif

static void resetPairs(ResultPairs *p)
{
  ARRLEN(p->mpr) = 0;
  ARRLEN(p->ivr) = 0;
  p->status = 0;
  p->n_proper = 0;
  p->n_within = 0;
}

static MATEPAIR *drawPairAtRandomByProbability(const MATEPAIRARR mpr)
{
  int i;
  const int n_pairs = ARRLEN(mpr);
  double pthresh;
  double s = 0.0;
  MATEPAIR *mp = NULL;
 
  for (i=0; i<n_pairs; i++)
    s += mpr[i].pbf;

  pthresh = RANDRAW_UNIFORM_1()*s;
  s = 0.0;

  for (i=0; i<n_pairs; i++) {
    s += mpr[i].pbf;
    if (s + MINLOGARG > pthresh) {
      mp = mpr + i;
      break;
    }
  }
  if (mp == NULL && n_pairs > 0) 
    mp = mpr + n_pairs - 1;
  
  return mp;
}

static int assignProbabilityToPairs(MATEPAIRARR mpr,
				   double *psum,
				   double *marga, double *margb,
				   RSLTPAIRFLG_t pairflg,
				   const InsHist *ihistp)
{
  int i, n_pairs = ARRLEN(mpr);
  double const prob_improper = CUMULPROB_IMPROPER;
  double const prob_proper = (1.0 - CUMULPROB_IMPROPER);
  double const prob_out = CUMULPROB_PROPER_OUTSIDE;
  double const prob_in = (1.0 - CUMULPROB_PROPER_OUTSIDE);
  double const prob_allout = prob_improper + prob_proper*prob_out;

  if (n_pairs < 1)
    return ERRCODE_SUCCESS;

  *psum = MINLOGARG;
  *marga = *margb = 0.0;
  for (i=0; i<n_pairs; i++) {
    double pa, pb, iab;
    RSLTFLG_t flga, flgb;
    MATEPAIR *mp = mpr + i;
    resultGetMapQualScore(&pa, &flga, mp->ap);
    resultGetMapQualScore(&pb, &flgb, mp->bp);
    if ((pairflg & RSLTPAIRFLG_RESTRICT_1st)) {
      if (pa > pb) pa = pb;
    } else if ((pairflg & RSLTPAIRFLG_RESTRICT_2nd)) {
      if (pb > pa) pb = pa;
    }
    if ( (mp->mapflg & MAPFLG_PROPER) ) {
      iab = prob_proper;
      if (  (mp->mapflg & MAPFLG_WITHIN) ) {
	if (ihistp == NULL || n_pairs < 2) {
	  iab *= prob_in;
	} else {
	  int32_t count, totnum;
	  double p;
	  count = insGetHistoCountCumulative(&totnum, 
					     (mp->ins < 0)? -1*mp->ins: mp->ins, 
					     1, ihistp); 
	  if (totnum < 1) {
	    totnum = 1;
	    count = 1;
	  }
	  /* the following is a temporary fix so long as insert size range
	     is specified independent of histogram */
	  p = ((double) count)/totnum;
	  if (p >= 0.5)
	    iab = 0.5 - p/2;
	  iab *= p*prob_in + prob_out;
	}
      } else {
	iab *= prob_out;
      }
    } else {
      iab = prob_improper;
    }
    mp->pbf = pa*pb*iab;
    
    *psum += mp->pbf;
    if ((flga &  RSLTFLAG_SINGLE)) {
      double s = (1.0-pa)*(prob_allout)*pb;
      *margb += s;
      *psum += s;
    }
    if ((flgb &  RSLTFLAG_SINGLE)) {
      double s = pa*prob_allout*(1.0-pb);
      *marga += s;
      *psum += s;
    }
  }

  return ERRCODE_SUCCESS;
}

static int scorePairsSimple(const Result **ap,
			    const Result **bp,
			    short *mapqA,
			    short *mapqB,
			    MAPFLG_t *mapflg,
			    int *nmax,
			    MATEPAIRARR mpr,
			    RSLTPAIRFLG_t pairflg,
			    const InsHist *ihistp,
			    RESULTOUTFLG_t rsltouflg,
			    const ResultSet *rsrp,
			    const ResultSet *rsmp)
/**< Select a pair and assign scores to it.
 * \param ap Return pointer to a mapping of the 1st mate.
 * \param bp Return pointer to a mapping of the 2nd mate.
 * \param mapflg returns MAPFLG_MULT1ST |  MAPFLG_MULT2ND if degenerate mapping or 0.
 * \param mpr Set of pairs to be considered.
 * \param pairflg Combination of  RSLTPAIR_FLAG bit flags.
 * \param rsltouflg combination of  RESULTSET_OUTPUT_FLAGS. If (rsltouflg & RESULTFLG_RANDSEL) != 0
 *        make a random selection if there are multiple possible mappings.
 * \param rsrp Mapping results for 1st mate.
 * \param rsmp Mapping results of 2nd mate.
 *
 * This is like scorePairsSimple, except that inserts of 'proper' pairs within the expected
 * range are assigned a certain (prior) probability e.g. 99.7% for 3 standard dev
 * from the normal mean, and all other pairs (incl. pairs in correct orientation
 * on the same chromosome get assigned 0.2% probability.
 * 
 * Top mapping quality scores of each mate are converted to error probablilties (1-P_a), (1-P_b).
 * Mappings of the 2nd best scores get the probability of (1-P_a)/n_a, if n_a is the number of 
 * mappings with the same 2nd best score. 
 * 
 * A pair (A,B) has the probability P(a,b) = P_a * P_b * I_ab
 * (i) for independently and mapped pairs:

 * (I_ab = CUMULPROB_INSIZ_PROPER = 0.998 if proper and within range, I_ab = 0.002 else).
 * (ia) If there are N1_a multipe best mappings of read A (=> mapping score 0 at the moment), assign
 *      P_a = 1/N1_a as probability, for a single mapping P_a = P1_a = 1 - 10^(-M_a/10) with mapping score M_a
 *      
 * (ib) If there is a single best mapping score N1_a == 1, consider the N2_a mappings with the 2nd best
 *      Smith-Waterman score with a probability of P2_a = (1-P1_a)/N2_a 
 *      
 * (ii) One mate (e.g. mate B) was mapped dependent on the other
 *      P_b = P1_b if P1_b < P_a, else P_b = P_a
 *
 * SMALT assigns mapping scores per read mapping (marginal over all placements of the mate).
 *
 * sum_I(Pi(I)) * Pprop + Pother == 1
 * Prop = 1.0 - CUMULPROB_IMPROPER: a priori probability of a 'proper' pair.
 * Pother =  CUMULPROB_IMPROPER;
 * Pabi = pa*pb*(pi(a-b)*Pprop + Pother);
 */
{
  int errcode;
  int i, n_pairs;
  double psum = MINLOGARG, marga = 0.0, margb=0.0;
  double maxprob = 0.0;
  MATEPAIR *mp = NULL;

  *ap = *bp = NULL;
  *mapqA = *mapqB = 0;
  *mapflg = 0;
  *nmax = 0;

  n_pairs = ARRLEN(mpr);
  
  if (n_pairs == 0) {
    /* replace resultSetGetTopIndex for getTopResult */
    *ap = resultSetGetTopResult(mapflg,(unsigned char)((rsltouflg & RESULTFLG_RANDSEL) != 0), rsrp);
    *bp = resultSetGetTopResult(mapflg, (unsigned char)((rsltouflg & RESULTFLG_RANDSEL) != 0), rsmp);
    return ERRCODE_SUCCESS;
  } 

  errcode = assignProbabilityToPairs(mpr, 
				     &psum, 
				     &marga, &margb,
				     pairflg,
				     ihistp);
  if ((errcode))
    return errcode;
    
  if (psum < MINLOGARG)
    psum = MINLOGARG;

  qsort(mpr, n_pairs, sizeof(MATEPAIR), cmpMATEPAIRbyProbDescending);
    
  /* establish, how many pairs there are with the same probability */
  for (i=1; i<n_pairs; i++) {
    if (mpr[i].pbf + MINLOGARG < mpr[0].pbf)
      break;
  }
  *nmax = i;
  mp = mpr;
  maxprob = mp->pbf/psum;
  
  if (maxprob <= 0.6 && n_pairs > 1) {
    /* make a random selection */
    *mapflg = MAPFLG_MULT1ST | MAPFLG_MULT2ND;
    if ((rsltouflg & RESULTFLG_RANDSEL))
      mp = drawPairAtRandomByProbability(mpr); 
    else if (!(rsltouflg & RESULTFLG_SINGLE)) 
      mp = mpr;
    else 
      mp = NULL;
  }
  if (NULL == mp)
    return ERRCODE_SUCCESS;

  *ap = mp->ap;
  *bp = mp->bp;
  *mapflg = (MAPFLG_t) (*mapflg | mp->mapflg);
  /* *mapflg |= mp->mapflg; results in warnings by the intel compiler 10.1 */
  /* calculate marginals for each mate */      
  for (i=0; i<n_pairs; i++) {
    if (mpr[i].ap == mp->ap) {
      marga += mpr[i].pbf;
    }
    if (mpr[i].bp == mp->bp)
      margb += mpr[i].pbf;
  }
  *mapqA = resultConvertProbabilityToMappingScore(marga/psum);
  *mapqB = resultConvertProbabilityToMappingScore(margb/psum);

  return ERRCODE_SUCCESS;
}
#ifdef RESULTPAIRS_SUPERFLUOUS
static int flagOutPairsWithDoubleFragments(const ResultPairs *pairp,
					   const ResultSet *rsrp, 
					   const ResultSet *rsmp)
{
  int errcode;
  MATEPAIRARR mpr = pairp->mpr;
  size_t npairs = ARRLEN(pairp->mpr);
  MAPFLG_t const mask = MAPFLG_PROPER | MAPFLG_WITHIN | MAPFLG_INVALID;
  MAPFLG_t const properflg = MAPFLG_PROPER | MAPFLG_WITHIN;
  int nn, i, j;
  short qsegnoA, qsegnoB;

  if (npairs > INT_MAX)
    return ERRCODE_OVERFLOW;

  nn = (int) npairs;
  if ((errcode = resultSetGetNumberOfSegments(NULL, &qsegnoA, rsrp)))
    return errcode;

  if ((errcode = resultSetGetNumberOfSegments(NULL, &qsegnoB, rsmp)))
    return errcode;

  if (nn < 2 || (qsegnoA < 1 && qsegnoB < 1))
    return ERRCODE_SUCCESS;
  
  for (i=0; i<nn; i++) {
    if ((mpr[i].mapflg&mask) == properflg) {
      /* is first proper pair, flag out all other pairs with complementary fragments */
      int sxA = resultGetFragmentNo(mpr[i].ap);
      int sxB = resultGetFragmentNo(mpr[i].bp);
 
      for (j=0; j<nn; j++) {
	if (!(mpr[j].mapflg&MAPFLG_INVALID) &&
	    j != i &&
	    (mpr[j].mapflg&properflg) != properflg &&
	    (resultGetFragmentNo(mpr[j].ap) == sxA ||
	     resultGetFragmentNo(mpr[j].bp) == sxB)
	    )
	  mpr[j].mapflg |= MAPFLG_INVALID;
      }
    }
  }
  /* skip the flagged out pairs */
  for (j=0; j<nn && !(mpr[j].mapflg&MAPFLG_INVALID); j++);
  for (i=j+1; i<nn; i++) {
    if (!(mpr[i].mapflg&MAPFLG_INVALID)) 
      mpr[j++] = mpr[i];
  }
  ARRLEN(mpr) = (size_t) j;

  return ERRCODE_SUCCESS;
}
#endif

static int addPairResultsToReport(Report *rep, 
				  const MAPFLG_t mapflg, REPMATEFLG_t repmateflg,
				  const Result *rp, short mapqA, const ResultSet *rsrp, 
				  const Result *mp, short mapqB, const ResultSet *rsmp) 
{
  int errcode;
  int isize = 0;
  int pairID = reportNextPairID(rep);
  REPPAIRFLG_t reppairflg = 0;
  REPMATEFLG_t rmAflg, rmBflg;

  if (pairID < 0) 
    return ERRCODE_NOMEM;

  repmateflg |= REPMATEFLG_PAIRED;
  
  if ((mapflg & MAPFLG_PAIRED) && 
      rp != NULL && mp != NULL &&
      !(resultGetStatusFlag(rp)&RSLTFLAG_NOOUTPUT) &&
      !(resultGetStatusFlag(mp)&RSLTFLAG_NOOUTPUT)
      ) {
    reppairflg |= REPPAIR_MAPPED;
    if (mapflg & MAPFLG_CONTIG) {
      reppairflg |= REPPAIR_CONTIG;
      resultCalcInsertSize(&isize, RSLTSAMSPEC_V1P4, rp, mp);
      if (mapflg & MAPFLG_WITHIN)
	reppairflg |= REPPAIR_WITHIN;
      if (mapflg & MAPFLG_PROPER)
	reppairflg |= REPPAIR_PROPER;
    }
  }   

  /* print 1st mate */
  rmAflg = (REPMATEFLG_t) (repmateflg & ~REPMATEFLG_2NDMATE);
  if (mapflg & MAPFLG_MULT1ST)
    rmAflg |= REPMATEFLG_MULTI;
  if ((errcode = resultSetAddResultToReport(rep,
					    pairID,
					    mapqA,
					    rmAflg, 
					    reppairflg,
					    isize,
					    rp, rsrp)))
    return errcode;

  /* print 2nd mate */
  rmBflg = (REPMATEFLG_t) (repmateflg | REPMATEFLG_2NDMATE);
  if (mapflg & MAPFLG_MULT2ND)
    rmBflg |= REPMATEFLG_MULTI;
  if ((errcode = resultSetAddResultToReport(rep,
					    pairID,
					    mapqB,
					    rmBflg, 
					    reppairflg,
					    isize,
					    mp, rsmp)))
    return errcode;

  return ERRCODE_SUCCESS;
}

/******************************************************************************
 *********************** Public Methods of Type ResultPairs *******************
 ******************************************************************************/

ResultPairs *resultSetCreatePairs(short blksz)
{
  ResultPairs *p;

  EMALLOCP0(p);
  if (!p) return 0;
  
  if (blksz < 1) blksz = DEFAULT_BLOCKSIZ_PAIRS;
  ARRCREATE(p->ivr, blksz);
  ARRCREATE(p->mpr, blksz);

  if ((p->ivr) && (p->mpr)) {
    resultSetBlankPairs(p);
  } else {
    resultSetDeletePairs(p);
    p = 0;
  }
  
  return p;
}

void resultSetDeletePairs(ResultPairs *p)
{
  if ((p)) {
    ARRDELETE(p->ivr);
    ARRDELETE(p->mpr);
  }
  free(p);
}

void resultSetBlankPairs(ResultPairs *p)
{
  if ((p)) {
    ARRLEN(p->ivr) = 0;
    ARRLEN(p->mpr) = 0;
    p->status = 0;
    p->n_proper = 0;
    p->n_within = 0;
    p->dmin = 0;
    p->dmax = 0;
  }
}

int resultSetFindPairs(ResultPairs *pairp,
		       RSLTPAIRFLG_t pairflg,
		       RSLTPAIRLIB_t pairlibcode,
		       int dmin, int dmax,
		       const ResultSet *rsltAp, const ResultSet *rsltBp)
{
  int errcode;
  BOOL_t isSingleA, isSingleB;
  struct GETMATEPAIRARG_ arg;
  /* 'A': first read, 'B': 2nd read */
  resetPairs(pairp);
  ARRLEN(pairp->mpr) = 0;
  pairp->n_proper = 0;
  pairp->n_within = 0;
  pairp->status = 0;
  if (dmin > dmax) {
    pairp->dmin = dmax;
    pairp->dmax = dmin;
  } else {
    pairp->dmin = dmin;
    pairp->dmax = dmax;
  }

  isSingleA = resultSetGetRankDepth(rsltAp, NULL, &arg.max_rankA);
  isSingleB = resultSetGetRankDepth(rsltBp, NULL, &arg.max_rankB);
  if ((pairflg & RSLTPAIRFLG_RESTRICT_2nd) && (isSingleA))
    arg.max_rankA = 0;
  else if ((pairflg & RSLTPAIRFLG_RESTRICT_1st) && (isSingleB))
    arg.max_rankB = 0;
   
  arg.pairflg = pairflg;
  arg.pairlibcode = pairlibcode;
  arg.ap = NULL;
  arg.rsltAp = rsltAp;
  arg.rsltBp = rsltBp;
  arg.pairp = pairp;

  errcode = resultSetDo(&arg, getMATEPAIRcbf_OuterLoop, rsltAp);
  if (!(errcode)) {
    pairp->status |= RSLTPAIRSFLG_MATES | RSLTPAIRSFLG_NPROPER | RSLTPAIRSFLG_INSLIM;
    pairp->status &= ~RSLTPAIRSFLG_SORTED;
  }

  return errcode;
}

int resultSetFindProperPairs(ResultPairs *pairp, int dist_lo, int dist_hi,
			     int maxnum,
			     int swscor_min,
			     RSLTPAIRLIB_t pairlibcode,
			     const ResultSet *rsltAp, const ResultSet *rsltBp)
{ 
  int errcode;
  short nresA, nresB;
  /* A: first mate, B: second mate of a pair */
  resultSetGetNumberOfSegments(&nresA, NULL, rsltAp);
  resultSetGetNumberOfSegments(&nresB, NULL, rsltBp);
  if (nresA < 1 || nresB < 1) {
    resetPairs(pairp);
    return ERRCODE_SUCCESS;
  }
  
  if (!(errcode = generateOFFSIVAL(&pairp->ivr, dist_lo, dist_hi,
				   rsltAp))) {
#ifdef results_debug
    fprintOFFSIVAL(stdout, pairp->ivr);
#endif
    pairp->status |= RSLTPAIRSFLG_INTVAL;

    if (swscor_min < 1) {
      int sw2ndmax = 0;
      int swmax = resultSetGetMaxSwat(rsltBp, &sw2ndmax);
      swscor_min = (sw2ndmax > 0)? sw2ndmax: swmax;
    }
    if (dist_lo > dist_hi) {
      pairp->dmin = dist_hi;
      pairp->dmax = dist_lo;
    } else {
      pairp->dmin = dist_lo;
      pairp->dmax = dist_hi;
    }

    errcode = findProperMATEPAIR(pairp, rsltBp,
				 swscor_min, pairp->ivr, pairlibcode,
				 maxnum);

    pairp->n_proper = (int) ARRLEN(pairp->mpr);

    if (errcode == ERRCODE_PAIRNUM) {
      pairp->status |= RSLTPAIRSFLG_CEILING;
      errcode = ERRCODE_SUCCESS;
    }
    
    pairp->status |= (RSLTPAIRSFLG_MATES | RSLTPAIRSFLG_MAPFLG | RSLTPAIRSFLG_INSLIM);
    pairp->status &= ~(RSLTPAIRSFLG_SORTED | RSLTPAIRSFLG_NPROPER);
  }

  return errcode;
}

int resultSetGetNumberOfPairs(int *n_proper, const ResultPairs *pairp)
{
  if (n_proper) *n_proper = pairp->n_proper;
  return (int) ARRLEN(pairp->mpr);
}

int resultSetAddPairToReport(Report *rep,
			     const InsHist *ihistp,
			     const ResultPairs *pairp,
			     RSLTPAIRFLG_t pairflg,
			     RESULTOUTFLG_t rsltouflg,
			     const ResultSet *rsrp,
			     const ResultSet *rsmp
			     )
{
  int errcode = ERRCODE_SUCCESS;
  int n_max = 0;
  short mapqA, mapqB;
  MAPFLG_t mapflg = 0;
  REPMATEFLG_t repmateflg = REPMATEFLG_PAIRED;
  const Result *ap = NULL, *bp = NULL;

#ifdef results_debug
  printf("resultSetPrintPair: dmin=%i, dmax=%i, paiflg=%i, rsltouflg=%i\n",
	 pairp->dmin, pairp->dmax, (int) pairflg, (int) rsltouflg);
#endif
  //flagOutPairsWithDoubleFragments(pairp, rsrp, rsmp);

  if ((errcode = scorePairsSimple(&ap, &bp, &mapqA, &mapqB,
				  &mapflg, &n_max, pairp->mpr,
				  pairflg, ihistp,
				  rsltouflg,
				  rsrp, rsmp)))
    return errcode;
     
  if (n_max > 1 && !(rsltouflg&RESULTFLG_RANDSEL) && (rsltouflg&RESULTFLG_SINGLE)) {
    /* report only unique mate */
    unsigned char isMultiA, isMultiB;
    ap = resultSetGetTopResult(&isMultiA, 0, rsrp);
    bp = resultSetGetTopResult(&isMultiB, 0, rsmp);
    if (!(isMultiA)) {
      bp = NULL;
      mapflg |= MAPFLG_MULT2ND;
    } else if (!(isMultiB)) {
      ap = NULL; 
      mapflg |= MAPFLG_MULT1ST;
    } else {
      mapflg |= MAPFLG_MULT1ST | MAPFLG_MULT2ND;
      ap = NULL;
      bp = NULL;
    }
  }

  if ((errcode = addPairResultsToReport(rep, mapflg, 
					(REPMATEFLG_t) (repmateflg | 
							REPMATEFLG_PRIMARY), 
					ap, mapqA, rsrp, 
					bp, mapqB, rsmp)))
    return errcode;
  if ((mapflg & (MAPFLG_MULT1ST | MAPFLG_MULT2ND)) && 
      !(rsltouflg&RESULTFLG_RANDSEL) && 
      !(rsltouflg&RESULTFLG_SINGLE)) {
    int i;
    for (i=0; i<n_max; i++) {
      MATEPAIR *mp = pairp->mpr + i;
      if (mp->ap != ap || mp->bp != bp) {
	MAPFLG_t mflg = (MAPFLG_t) (mp->mapflg | (mapflg & (MAPFLG_MULT1ST | MAPFLG_MULT2ND)));
	if ((errcode = addPairResultsToReport(rep, mflg, 
					      (REPMATEFLG_t) (repmateflg | 
							      REPMATEFLG_PRIMARY), 
					      mp->ap, mapqA, rsrp, 
					      mp->bp, mapqB, rsmp)))
	  return errcode;
      }
    }
  }

  /* secondary alignments */
  if ((rsltouflg&RESULTFLG_BEST) && (rsltouflg&RESULTFLG_SPLIT)) {
    if ((errcode = resultSetAdd2ndaryResultsToReport(rep, 
						     (REPMATEFLG_t) (repmateflg | 
								     REPMATEFLG_PARTIAL),
						     rsltouflg,
						     rsrp)))
      return errcode;
   
  if ((errcode = resultSetAdd2ndaryResultsToReport(rep, 
						   (REPMATEFLG_t) (repmateflg | 
								   REPMATEFLG_PARTIAL | 
								   REPMATEFLG_2NDMATE), 
						   rsltouflg,
						   rsmp)))
      return errcode;
  }
  
  return errcode;
}
