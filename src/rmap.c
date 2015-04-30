/** K-tuple matching and alignment strategies for read mapping */

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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <limits.h>
#include <stdint.h>

#include "array.h"
#include "interval.h"
#include "hashidx.h"
#include "segment.h"
#include "score.h"
#ifdef SCORE_SIMD
#include "swsimd.h"
#endif
#include "alignment.h"
#include "rmap.h"
#ifdef RESULTS_TRACKER
#include "tracker.h"
#endif

enum {
  HASH_MAXNHITS = 16*1024,      /**< maximum number of k-tuple hits */
  HITINFO_BLKSZ =  1024,        /**< Block size for hit info */
  HITBINS_TARGETNUM = 8192,     /**< Number of hit bins for filtering ktuple hits */
  HITFILTR_BLKSZ = 128,         /**< Block size for hit filter */
  SEGLST_BLOCKSIZ = 64*1024,    /**< memory allocation block size for segment list */ 
  SEGCANDS_BLOCKSIZ = 32*1024,  /**< memory allocation block size for candidate segments */
  RMAPCANDS_BLOCKSIZ = 512,     /**< memory allocation block size for filtered candidate segments */
  EDGELEN_MAX = 500,            /**< maximum edge length */
  UNKNOWN_SEQIDX = -1,          /**< Value of sequence index < 0 if unknown */
  UNKNOWN_SYMBOL = 'N',         /**< Nulceotide symbol for unknown character */
  MAXNUM_CANDSEG = 2000,        /**< maximum number of candidate segments considered */
  MAXN_CONSECUTIVE_MISMATCH = 2,/**< maximum number of consecutive mismatches allowed in
				 * ungapped alignment */
  DEFAULT_SCORE_MATCH = 1,      /**< Default score for a match in ungapped alignment */
  DEFAULT_SCORE_MISMATCH = -2,  /**< Default score for a mismatch in ungapped alignment */
  NUM_STRANDS = 2,              /**< For loops over forward [0] and reverse [1] strands */
  FILTERIVALEXT = 30,           /**< Extenstion (% of read length) for search intervals of
				 *   paired reads.*/
  FILTERIVAL_BLKSZ = 1028,      /**< Block size for momory allocation for array of search 
				 * intervals */
  MAPSCORE_UNIQUE_MAPPED_1ST = 20, /**< Mapping score threshold that defines a 'unique' hit for 
				    * the rare mate which is mapped first. If the mapping score 
				    * satisfies this threshold, the 2nd mate is found via a 
				    * restricted mapping */
  MAXNUM_PAIRS_TOTAL = 1028,     /**< Maximum number of pairs to be considered */
  MAXNUM_PAIRS_PROPER = 2,      /**< Maximum number of pairs to consider when looking for
				 * proper pairs */
  KTUPCOUNT_REDUCTION_FACTOR = 2, /**< Reduce the minimum count of k-tuples by this amount when 
				   * filtering bins */
  KTUPCOUNT_PERBIN_THRESHOLD = 2, /**< Minimum number of ktuple hits for a bin to pass filtering */
  KTUPCOUNT_THRESHOLD_FOR_BINNING = 256, /**< If the number of k-tuple counts is above that threshold,
					  * use binning rather than pure sorting */
  SEED_PERCENTILE = 80,          /**< uses only least frequent seeds as percentage of total */
  MINLEN_QUERY_STRIPED = 32,     /**< minimum read length Smith-Waterman using SIMD */
  BWSCAL_QLEN = 16*3,            /**< The length of the profiled sequence must be less than the bandwidth
				  * times BWSCAL_QLEN so that it pays to run vectorised rather than 
				  * in the band and scalar. */
  MINSCOR_BELOW_MAX_BEST = 0,    /**< look only for best mapping */
  MINSCOR_BELOW_MAX_ALL = -1,    /**< consider all mappings */

#ifdef rmap_finehash_2ndmate
  FINEHASH_WORDLEN = 5,          /**< Word length for fine on-the-fly hashing */
  FINEHASH_SKIPSTEP = 1,         /**< Skip step for fine on-the-fly hashing */
  FINEHASH_MAXKTUPPOS = 128*1024*1024, /**< maximum number of word positions for fine on-the-fly hashing */
#endif
  SPLITREAD_OVERLAP_PERCENT = 60, /**< Max. Percentage of overlap allowed for split reads */
#ifdef RESULTS_TRACKER
  RESULTS_TRACKER_IS1BASED = 1,
#endif
};

enum RMAPCAND_FLAGS {
  RMAPCANDFLG_REVERSE = 0x01,
  RMAPCANDFLG_SCORED = 0x02,
};

typedef unsigned char UCHAR;
typedef signed int SWATSCOR;
typedef unsigned char BOOL;
typedef uint32_t COVERAGE;

typedef struct RMAPCAND_ {  /**< Candidate segment aligned using fast methods */
  RMAPFLG_t flags;  /**< one of RMAPCANDS_FLAGS */
  SEQLEN_t qe;      /**< last base of query segment (0 based)  */
  SEQLEN_t qs;      /**< first base of query segment (0 based) */
  uint64_t rs;      /**< first base of reference segment (0 based) */
  uint64_t re;      /**< last base of query segment (0 based) */
  int band_l;       /**< left edge of alignment band along (reverse complemented query) **/
  int band_r;       /**< right edge of alignment band along (reverse complemented query) **/
  int64_t sqidx;    /* sequence index, if unknown: SEGCAND_UNKNOWN_SEQIDX */
#ifdef results_debug
  int segidx;
#endif
  uint32_t dqo; /**< query offset for direct alignment */
  int dro;    /**< refrence offset for direct alignment */
  SWATSCOR swscor;
} RMAPCAND;

typedef RMAPCAND * RMAPCANDARR;

typedef struct RMAPBUFF_ {  /**< Various buffers for mapping a single read */
  SeqFastq *sqbfp;   /**< Sequence buffer */
  SeqFastq *qbfp;    /**< Sequence buffer for explicit alignment output */
  HashHitList *hhlp; /**< Hit list with k-tuple hits */
  SegLst *sglp;
  SegAliCands *sacp;
  RMAPCANDARR candr;
  AliBuffer *alibufp;
  AliRsltSet *alirsltp;
  SegQMask *qmp;
} RMAPBUFF;

typedef struct RMAPINFO_ {
  HashHitInfo *hhiFp; /* K-tuple stats forward strand */
  HashHitInfo *hhiRp; /* K-tuple stats reverse strand */
} RMAPINFO;

typedef struct RMAPPROF_ {
  SeqFastq *readRCp;  /**< Reverse complement of read */
  ScoreProfile *scorprofp;  /**< Profile for reads */
  ScoreProfile *scorprofRCp; /**< Profile for reverse complement of read */
} RMAPPROF;

struct RMap_ {  
  RMAPBUFF *bfp;  /**< Buffers for mapping a single read */
  RMAPPROF *prp;  /**< Profile for read */
  RMAPPROF *pmp;  /**< Profile for mate */
  RMAPINFO *mrp;  /**< read */
  RMAPINFO *mmp;  /**< mate of read (if paired-end, NULL if single read) */
  RMAPINFO *mr2p; /**< for secondary alignment of read (can be NULL) */
  RMAPINFO *mm2p; /**< for secondary alignment of mate (can be NULL) */
  /* HashHitFilter *hhfp; /\**< for filtering hits of the mate, NULL if single reads *\/ */
  InterVal *ivr; /**< for filtering hits of the mate, NULL if single reads */
  ResultSet *rsrp; /**< Set of alignment results for read */
  ResultSet *rsmp; /**< Set of alignment results for mate */
  ResultPairs *pairp;

#ifdef rmap_finehash_2ndmate
  RMAPINFO *mflyp; /**< mate of read (if paired-end, NULL if single read) */
  HashTable *htflyp;  /**< hash table for on-the-fly hashing */
#endif
};

static const float MINFRACT_MAXSCOR_2ND = 0.8;

/******************************************************************************
 ******************************* Private Methods ******************************
 ******************************************************************************/
static BOOL scorIsAboveFractMax(int scor_read, int scor_mate, float fract, 
				const SeqFastq *readp, const SeqFastq *matep)
{
  SEQLEN_t rlen, mlen;
  seqFastqGetConstSequence(readp, &rlen, NULL);
  seqFastqGetConstSequence(matep, &mlen, NULL);
  return (BOOL) (scor_read >= scor_mate*rlen*fract/mlen);
}

#ifdef RMAP_SUPERFLUOUS_CODE
static int tryUngappedAlignment(int *maxscor,
				const char *qstr, const char *rstr, const RMAPCAND *cp)
     /**< Try ungapped alignment with 1 for match -2 for mismatch
      */
{
  int dlen, maxlen;
  int j, mmconsec_ctr;
  int swscor;

  *maxscor = swscor = 0;

  if ((cp->dqo > cp->qe + 1) ||
      (cp->qe > INT_MAX + cp->dqo))
    return ERRCODE_ASSERT;
  maxlen = (int) (cp->qe - cp->dqo + 1);

  if ((((uint64_t)cp->dro) > cp->re + 1) ||
      (cp->re > (uint64_t) INT_MAX + cp->dro))
    return ERRCODE_ASSERT;
  dlen = cp->re + 1 - cp->dro;

  if (dlen < maxlen)
    maxlen = dlen;

/*   if (maxlen*DEFAULT_SCORE_MATCH < minscor) */
/*     return ERRCODE_SUCCESS; */

  qstr += cp->dqo;
  rstr += cp->dro;
  mmconsec_ctr = 0;
  for (j=0; j<maxlen && (rstr[j]) && (qstr[j]); j++) {
    if (!(rstr[j]&SEQCOD_STDNT_TESTBIT) && !(qstr[j]&SEQCOD_STDNT_TESTBIT)) {
      if ((qstr[j]&SEQCOD_STDNT_MASK) == (rstr[j]&SEQCOD_STDNT_MASK)) {
	if (mmconsec_ctr) 
	  mmconsec_ctr--;
	swscor += DEFAULT_SCORE_MATCH;
	if (swscor > *maxscor)
	  *maxscor = swscor;
      } else {
	if (++mmconsec_ctr > MAXN_CONSECUTIVE_MISMATCH)
	  return ERRCODE_FAILURE;
	swscor += DEFAULT_SCORE_MISMATCH;
	if (swscor < 0)
	  swscor = 0;
      }
    }
  }

  return ERRCODE_SUCCESS;
}
#endif

static uint32_t calcMinKtup(uint32_t *mincover, const HashTable *htp)
{
  UCHAR nskip;
  UCHAR ktup = hashTableGetKtupLen(htp, &nskip); 
  uint32_t minktup = (*mincover >= ktup + nskip)? ((*mincover)-ktup)/nskip: 1;
  *mincover = (minktup-1)*nskip + ktup;
  return minktup;
}

#ifdef RMAP_SUPERFLUOUS_CODE
static uint32_t calcMinCover(int *min_ktup, UCHAR mincov_percent, 
			   const SeqFastq *sqp, const HashTable *htp)
     /** Return the minimum covarage as the number of bases and
      * as the minimum number of k-tuples corresponding to a fraction of 
      * the read length covered by k-tuples.
      */
{
  uint32_t readlen, mincover;
  UCHAR nskip;

  hashTableGetKtupLen(htp, &nskip);
  seqFastqGetConstSequence(sqp, &readlen, NULL); 

  if (mincov_percent > 100)
    mincover = readlen;
  else 
    mincover = readlen*mincov_percent/100;
  *min_ktup = calcMinKtup(&mincover, htp);

  return mincover;
}
#endif

static int collectHits(SegAliCands *sacp, BOOL with_seqidx, 
		       SegLst *slp, HashHitList *hlp, SegQMask *qmp,
#ifdef RESULTS_TRACKER
		       Track *trkp,
#endif
#ifdef hashhit_dump_sortarray
		       FILE *dumpfp,
		       int *dumpctr,
#endif
		       uint32_t n_hit_max, uint32_t n_ktup_min, uint32_t cover_min,
		       const HashHitInfo *hip, const HashTable *htp, const SeqSet *ssp)
{
  int errcode = ERRCODE_SUCCESS;

  if ((with_seqidx)) {
    /* process each sequence separately */
    SEQNUM_t s;
    const SETSIZ_t *soffsp;
    const SEQNUM_t nseq = seqSetGetOffsets(ssp, &soffsp);
    
    for (s=0; s<nseq; s++) {
      hashBlankHitList(hlp);
      if ((errcode = hashCollectHitsForSegment(hlp, 
#ifdef RESULTS_TRACKER
					       trkp,
#endif
#ifdef hashhit_dump_sortarray
					       dumpfp,
					       dumpctr,
#endif
					       soffsp[s], soffsp[s+1], 
#ifndef hashhit_minimise_coverdeficit
					       n_hit_max,
					       1,
#endif
					       hip, htp, NULL)))
	break;
      segLstBlank(slp);
      if ((errcode = segLstFillHits(slp, 
#ifdef RESULTS_TRACKER
				    trkp,
#endif
				    n_ktup_min, hlp)))
	break;
    
      if ((errcode = segAliCandsAddFast(sacp, qmp,
#ifdef RESULTS_TRACKER
					trkp,
#endif 
					slp, cover_min, s)))
	break;
    }
  } else {
    /* treat reference sequence as one (concatenated) sequence
     * find out sequence index later */
    errcode = hashCollectHitsUsingCutoff(hlp, 
#ifdef hashhit_dump_sortarray
					 dumpfp,
					 dumpctr,
#endif
					 n_hit_max, htp, hip);
    if (!errcode) {
      segLstBlank(slp);
      errcode = segLstFillHits(slp, 
#ifdef RESULTS_TRACKER
			       trkp,
#endif
			       n_ktup_min, hlp);
      if (!errcode) 
	errcode = segAliCandsAddFast(sacp, qmp, 
#ifdef RESULTS_TRACKER
				     trkp,
#endif
				     slp, cover_min, SEGCAND_UNKNOWN_SEQIDX);
    }
  }
  
  return errcode;
}


static int setupInterValFromResultSet(InterVal *ivr, int dmin, int dmax, 
				      const SeqFastq *readp, const SeqFastq *matep,
				      const HashTable *htp, const SeqSet *ssp,
				      const ResultSet *rsp)
     /**< Derrive from the best hits in the set of results, an array 
      * of distance intervals that can be used to filter seeds from
      * the hash table.
      * \param ivr Set of intervals
      * \param dmin Minimum insert size.
      * \param dmax Maximum insert size.
      * \param readp Read (for which rsp contains results).
      * \param matep Mate of read.
      * \param htp Hash index of k-mer words.
      * \param ssp Set of reference seqences.
      * \param rsp Results for read of length readlen.
      *
      * \note There is FILTERIVALEXT% of the read length added to allow for insertions.
      *       Assumes that sequence indices have been assigned.
      */
{
  int errcode = ERRCODE_SUCCESS;
  UCHAR nskip;
  UCHAR ktup = hashTableGetKtupLen(htp, &nskip);
  short i, n;
  SEQNUM_t sx;
  uint32_t qs, qe, rs, re, delta, rlen, readlen, matelen;
  int64_t tmp_lo, tmp_hi;
  RSLTFLG_t status;
  const SEQNUM_t nseq = seqSetGetOffsets(ssp, NULL);

#define ADJUST_INTERVAL(tmp, rlen) if ((tmp) >= (rlen)) (tmp) = ((int64_t) (rlen)) - 1; \
                                   if ((tmp) < 1) (tmp) = 0;
  if (dmin > dmax) 
    return ERRCODE_ARGRANGE;
  seqFastqGetConstSequence(readp, &readlen, NULL);
  seqFastqGetConstSequence(matep, &matelen, NULL);
  delta = (((int64_t) matelen)*FILTERIVALEXT)/100;

  resultSetGetScorStats(rsp, NULL, &n, NULL, NULL);
  interValBlank(ivr);

  for (i=0; i<n; i++) {
    const Result *rp;
    if ((errcode = resultSetGetResultByRank(&rp, i, rsp)))
      return errcode;
    if ((errcode = resultGetData(&qs, &qe, &rs, &re, &sx, NULL, &status, rp)))
      return errcode;
    if (!((status&RSLTFLAG_SELECT) && ((re) > (rs)) && (sx >= 0) && (sx < nseq))) {
      errcode = ERRCODE_ASSERT;
      break;
    }
    /* qs, qe, rs, re start from 1 !, sx < 0 means sequence indices have no yet been assigned */
    
    rlen =  seqSetGetSeqDatByIndex(NULL, NULL, sx, ssp);

    /* lower interval: (re-1) + readlen - qe  - dmax;
     *                 (re-1) + readlen - qe  - dmin + matelen - ktup + delta */
    tmp_lo = ((int64_t) re) + readlen - qe - dmax;
    ADJUST_INTERVAL(tmp_lo, rlen);

    tmp_hi = ((int64_t) re) + readlen + matelen + delta - qe - dmin  - ktup;
    ADJUST_INTERVAL(tmp_hi, rlen);

    if (tmp_lo <= tmp_hi &&
	(errcode = interValAppend(ivr, (SEQLEN_t) tmp_lo, (SEQLEN_t) tmp_hi, sx, status))) {
      break;
    }
    /* upper interval: rs - qs + dmin - matelen;
     *                 rs - qs + dmax - ktup + delta;
     */
    tmp_lo = ((int64_t) rs) - qs + dmin - matelen;
    ADJUST_INTERVAL(tmp_lo, rlen);

    tmp_hi = ((int64_t) rs) - qs + dmax - ktup + delta;
    ADJUST_INTERVAL(tmp_hi, rlen);

    if (tmp_lo <= tmp_hi &&
	(errcode = interValAppend(ivr, (SEQLEN_t) tmp_lo, (SEQLEN_t) tmp_hi, sx, status)))
      break;
  }

  return errcode;
}

static int collectHitsFromInterVal(SegAliCands *sacp, SegLst *slp, 
				   HashHitList *hlp, SegQMask *qmp,
#ifdef RESULTS_TRACKER
				   Track *trkp,
#endif
				   uint32_t n_hit_max, uint32_t n_ktup_min, uint32_t cover_min,
				   const HashHitInfo *hip,
				   const HashTable *htp, const SeqSet *ssp, 
				   const InterVal *ivr)
{
  int errcode = ERRCODE_SUCCESS;
  int i, niv = interValNum(ivr);
  SEQLEN_t lo, hi;
  SEQNUM_t sx;
  const SETSIZ_t *soffsp;
  SETSIZ_t offs;

  seqSetGetOffsets(ssp, &soffsp);
  for (i=0; i<niv; i++) {
    if ((errcode = interValGet(&lo, &hi, &sx, NULL, i, ivr)))
      break;

    hashBlankHitList(hlp);
    offs = soffsp[sx];
    if ((errcode = hashCollectHitsForSegment(hlp, 
#ifdef RESULTS_TRACKER
					     trkp,
#endif
#ifdef hashhit_dump_sortarray
					     NULL,
					     NULL,
#endif
					     offs + lo, offs + hi+1, 
#ifndef hashhit_minimise_coverdeficit
					     n_hit_max, 0, 
#endif
					     hip, htp, NULL)))
      break;
    segLstBlank(slp);
    if ((errcode = segLstFillHits(slp, 
#ifdef RESULTS_TRACKER
				  trkp,
#endif
				  n_ktup_min, hlp)))
      break;
    if ((errcode = segAliCandsAddFast(sacp, qmp,
#ifdef RESULTS_TRACKER
				      trkp,
#endif 
				      slp, cover_min, sx)))
      break;
  }

  return errcode;
}

#ifdef rmap_finehash_2ndmate
static int setupFineHashTable(HashTable *htfinep, 
			      SeqFastq *sqbfp,
			      const SeqSet *ssp, 
			      const InterVal *ivr,
			      const HashTable *htp,
			      const SeqCodec *codecp)
{
  uint32_t npos_max = FINEHASH_MAXKTUPPOS;
  int errcode = hashTableSetUp(htfinep, sqbfp, ssp, ivr, codecp, &npos_max, 0);
  if (errcode == ERRCODE_MAXKPOS) {
    int s = npos_max/FINEHASH_MAXKTUPPOS + 1;
    UCHAR nskip;
    UCHAR ktuplen = hashTableGetKtupLen(htp, &nskip);
    if (nskip > s || s > ktuplen) {
      hashTableReset(htfinep, (UCHAR) s);
      errcode = hashTableSetUp(htfinep, sqbfp, ssp, ivr, codecp, NULL, 0);
    }
  } 
  if (ERRCODE_KPOSOFLO == errcode)
    errcode = ERRCODE_MAXKPOS;

  return errcode;
}
#endif
/******************************************************************************
 *********************** Methods of Private Type RMAPCAND *********************
 ******************************************************************************/
#ifdef rmap_debug
static void fprintRMAPCAND(FILE *fp, const RMAPCAND *rcp)
{ 
  if (rcp) {
    fprintf(fp, "RMAPCAND: qe = %u, qs = %u, re = %llu, rs = %llu, sqidx = %i, band_l = %i, band_r = %i\n", 
	    rcp->qe, rcp->qs, (unsigned long long) rcp->rs, (unsigned long long) rcp->re,
	    rcp->sqidx, rcp->band_l, rcp->band_r);
  } else {
    fprintf(fp, "== RMAPCAND (NULL) ==");
  }
}
#endif

static int makeRMAPCANDfromSegment(RMAPCAND *cp, SeqFastq *sqbufp, COVERAGE *coverp,
				   SEQLEN_t qlen, const SeqSet *ssp, const SeqCodec *codecp,
				   uint32_t i, const SegAliCands *sacp)
{
  int errcode;
  uint8_t bitflags;
  SETSIZ_t rs, re;
  errcode = segAliCandsCalcSegmentOffsets(&cp->qs, &cp->qe, 
					  &cp->rs, &cp->re,
					  &cp->band_l, &cp->band_r, 
					  &cp->dqo, &cp->dro,
					  &cp->sqidx, &bitflags,
					  coverp,
#ifdef SCORE_SIMD
					  0,
#else
					  EDGELEN_MAX, 
#endif
					  qlen, 
					  ssp, i, sacp);
  if (errcode) 
    return errcode;

  /* for the moment, offsets still have to be representable by signed/unsigned ints */
  if (cp->qe > INT_MAX || cp->re < cp->rs || cp->re - cp->rs > INT_MAX)
    return ERRCODE_OVERFLOW;

  cp->flags = (RMAPFLG_t) ((bitflags & SEGCANDFLG_REVERSE)? RMAPCANDFLG_REVERSE: 0);

  cp->swscor = 0;
#ifdef results_debug
  cp->segidx = (int) i;
#endif
    

  re =  cp->re;
  rs =  cp->rs;
  if (cp->sqidx == SEGCAND_UNKNOWN_SEQIDX) {
    if ((errcode = seqSetFetchSegment(sqbufp, &rs, &re, ssp, codecp)))
      return errcode;
    cp->rs = rs;
    cp->re = re;
  } else {
    if ((errcode = seqSetFetchSegmentBySequence(sqbufp, cp->sqidx, rs, re-rs+1, ssp, codecp)))
      return errcode;
    /* can throw ERRCODE_SEQOFFS */
  }
    
  errcode = seqFastqEncode(sqbufp, codecp);

  return errcode;
}

static int scoreRMAPCAND(RMAPCANDARR *csr,
			 SWATSCOR *max1scor,
			 SWATSCOR *max2scor,
			 SeqFastq *sqbufp,
			 AliBuffer *alibufp,
			 RMAPFLG_t rmapflag,
			 UCHAR nskip,
#if defined alignment_debug || defined alignment_matrix_debug || defined rmap_debug
			 const char *profiled_seqp,
			 const char *profiled_seqRCp,
#endif
			 const ScoreProfile *profp,
			 const ScoreProfile *profRCp,
			 const SeqSet *ssp, const SeqCodec *codecp,
			 const SegAliCands *sacp)
     /**< Submit the candidate segments to a first pass of alignments using SIMD instructions.
      * csr Array of explicit candidate segments.
      * sqbufp Sequence buffer for obtaining a segment.
      */
{
  int errcode = ERRCODE_SUCCESS;
  RMAPCANDARR hp;
  short mismatchscor, gapinitscor;
  short matchscor = scoreProfileGetAvgPenalties(&mismatchscor, &gapinitscor, NULL, profp);
  short mmscordiff = (short) (matchscor - mismatchscor);
  COVERAGE cover, curr_min_cover, max_cover, min_cover, dcov, cdf;
  COVERAGE cover_deficit[NUM_STRANDS];
  SEQLEN_t i, qlen;
  const uint32_t n_candseg = segAliCandsGetNumberOfSegments(sacp, &curr_min_cover, NULL, 
							    cover_deficit, cover_deficit+1, 
							    NULL);
  RMAPCAND *cp;
#ifdef rmap_stop_candlist_early
  short gapscordiff = matchscor - gapinitscor;
  SWATSCOR max_possible_swscor;
  uint32_t cover_rank, cover_rank_break;
#endif
#ifdef rmap_debug
  const char *dp, *decoderp = seqCodecGetDecoder(codecp, NULL);
#endif
#ifdef SCORE_SIMD
  BOOL isSIMDAliCand;
#endif

  if (n_candseg > INT_MAX || mmscordiff < 1 
#ifdef rmap_stop_candlist_early
      || gapscordiff < 1
#endif
      )
    return ERRCODE_ASSERT;

  if (n_candseg >= ARRNALLOC(*csr)) {
    hp = ARREALLOC(*csr, n_candseg);
    if (hp == NULL) 
      return ERRCODE_NOMEM;
    *csr = hp;
  }

  scoreGetProfile(NULL, &qlen, NULL, NULL, profp);

  if (qlen*matchscor > INT_MAX)
    return ERRCODE_OVERFLOW;

  *max1scor = 0; *max2scor = 0;
  min_cover = 0; max_cover = 0;

#ifdef rmap_stop_candlist_early
  max_possible_swscor = qlen*matchscor;
  cover_rank = 0;
  cover_rank_break = n_candseg;
#endif

  for (cp = *csr, i = 0; i<n_candseg; i++, cp++) { 
    const char *unprofiled_seqp = NULL;
    SEQLEN_t unprofiled_seqlen = 0;
    const ScoreProfile *scprofp = NULL;  /* points to profp or profRCp */

#if defined alignment_debug || defined alignment_matrix_debug
    const char *profdseqp = NULL;        /* points to profiled_seqp or profiled_seqRCp */
#endif
    if ((errcode = makeRMAPCANDfromSegment(cp, sqbufp, &cover, 
					   qlen, ssp, codecp, 
					   i, sacp)))
      return errcode;

#ifdef rmap_stop_candlist_early
    if (cover < curr_min_cover) {
      curr_min_cover = cover;
      cover_rank++;
    }
    if (cover_rank > cover_rank_break)
      break;
#endif

    if ((cp->flags & RMAPCANDFLG_REVERSE)) {
      scprofp = profRCp;
#if defined alignment_debug || defined alignment_matrix_debug
      profdseqp = profiled_seqRCp;
#endif
    } else {
      scprofp = profp;
#if defined alignment_debug || defined alignment_matrix_debug
      profdseqp = profiled_seqp;
#endif
    }
  

    unprofiled_seqp = seqFastqGetConstSequence(sqbufp, &unprofiled_seqlen, NULL);
  
    if (unprofiled_seqlen > INT_MAX || unprofiled_seqlen != cp->re - cp->rs + 1)
      return ERRCODE_ASSERT;

#ifdef rmap_debug
    //printf("rmap_debug:Q::%s\n", qbasp);
    printf("rmap_debug:scoreRMAPCAND:Q::[%i]",i);
    for (dp = profiled_seqp; (*dp); dp++)
      fputc(decoderp[(UCHAR) *dp], stdout);
    fputc('\n', stdout);

    //printf("rmap_debug:S::%s\n", seqp);
    printf("rmap_debug:scoreRMAPCAND:S::[%i]", i);
    for (dp = unprofiled_seqp; *dp; dp++)
      fputc(decoderp[(UCHAR) *dp], stdout);
    fputc('\n', stdout);
#endif

#ifdef SCORE_SIMD
    isSIMDAliCand = (BOOL) 
      (qlen >= MINLEN_QUERY_STRIPED && 
       ((SEQLEN_t) (cp->band_r - cp->band_l)*BWSCAL_QLEN) > qlen &&
       cp->qs == 0 && cp->qe >= qlen-1);
    if ( (isSIMDAliCand) ) {
      errcode = swSIMDAlignStriped(&cp->swscor, 
				   alibufp, 
				   scprofp,
#ifdef  alignment_matrix_debug 
				   codecp,
				   profdseqp,
#endif
				   unprofiled_seqp,
				   unprofiled_seqlen);
    } 
    if ( !(isSIMDAliCand) || (errcode == ERRCODE_SWATEXCEED) ) {
#endif /* #ifdef SCORE_SIMD */
      if (cp->qs > INT_MAX)
	return ERRCODE_ASSERT;
      errcode = aliSmiWatInBandFast(&cp->swscor, 
				    alibufp, 
				    scprofp, 
#if defined alignment_matrix_debug
				    codecp,
				    profdseqp,
#endif
				    unprofiled_seqp, (int) unprofiled_seqlen,
				    cp->band_l, cp->band_r,
				    (int) cp->qs, (int) cp->qe,
				    0, (int) unprofiled_seqlen-1);
#ifdef SCORE_SIMD
    }
#endif
    if (errcode)
      return errcode;
      /*    } */

#ifdef rmap_debug
    printf("rmap_debug:scoreRMAPCAND:score =%i\n", cp->swscor);
#endif    
    cp->flags |= RMAPCANDFLG_SCORED;
    cdf = cover_deficit[(cp->flags & RMAPCANDFLG_REVERSE)? 1 :0];
    if ((rmapflag & RMAPFLG_BEST) && 
	(cover + cdf < min_cover))
	break;
    if (cp->swscor > (*max2scor)) {
      if (cp->swscor > (*max1scor)) {
	*max2scor = *max1scor;
	*max1scor = cp->swscor;
#ifdef rmap_stop_candlist_early
	if ((rmapflag & RMAPFLG_BEST)) {
	  if (*max1scor + gapscordiff > max_possible_swscor) {
	    if (*max1scor + mmscordiff > max_possible_swscor) {
	      cover_rank_break = 1;
	    } else {
	      cover_rank_break = 2;
	    }
	  }
	}
#endif
	if (cover + cdf > max_cover)
	  max_cover = (cover > cdf)? cover - cdf: 0;
      } else {	
	*max2scor = cp->swscor;
      }
      dcov = ((int)((*max1scor - *max2scor)/mmscordiff) + 1)*nskip;
      if (dcov + cdf + min_cover < max_cover)
	min_cover = max_cover - dcov;
    }
  }
  ARRLEN(*csr) = i;

  return ERRCODE_SUCCESS;
}

static int alignRMAPCANDFull(ResultSet *rsp,
			     AliRsltSet *alirsp,
			     AliBuffer *alibufp,
			     SeqFastq *sqbufp,
			     SWATSCOR min_swatscor,
			     int scorlen_min,
			     int bandwidth_min,
			     RMAPFLG_t rmapflag,
			     const ScoreProfile *profp,
			     const ScoreProfile *profRCp,
#if defined alignment_debug || defined alignment_matrix_debug || defined rmap_debug || defined results_debug
			     const char *profiled_seqp,
			     const char *profiled_seqRCp,
#endif
			     const SeqSet *ssp, 
			     const SeqCodec *codecp,
			     const RMAPCANDARR candr)
{
  int i, errcode = ERRCODE_SUCCESS;
  int band_l, band_r, bw;
  const char *unprofiled_seqp;
  const int scrlen = ARRLEN(candr);
  uint32_t unprofiled_seqlen;
  SEQLEN_t profiled_seqlen;
  SETSIZ_t rs, re;
  SWATSCOR swatscor_2ndmax = 0;
#ifdef rmap_debug
  const char *dp, *decoderp = seqCodecGetDecoder(codecp, NULL);
#endif

  for (i=0; i<scrlen; i++) {
    const ScoreProfile *scprofp = NULL;
#if defined alignment_debug || defined alignment_matrix_debug || defined rmap_debug || defined results_debug
    const char *profdseqp = NULL;
#endif
    RMAPCAND *cp = candr + i;
    if ((cp->flags & RMAPCANDFLG_SCORED) &&
	cp->swscor < min_swatscor)
      continue;
    if (cp->re < cp->rs || cp->re - cp->rs  > INT_MAX)
      return ERRCODE_OVERFLOW;
    re = cp->re;
    rs = cp->rs;
    if (cp->sqidx == SEGCAND_UNKNOWN_SEQIDX) {
      if ((errcode = seqSetFetchSegment(sqbufp, &rs, &re, ssp, codecp)))
	return errcode;
      cp->rs = rs;
      cp->re = re;
    } else {
      if ((errcode = seqSetFetchSegmentBySequence(sqbufp, cp->sqidx, rs, re-rs+1, ssp, codecp)))
	return errcode;
    /* can throw ERRCODE_SEQOFFS */
    }
   
    if ((errcode = seqFastqEncode(sqbufp, codecp)))
      return errcode;

    if ((cp->flags & RMAPCANDFLG_REVERSE)) {
      scprofp = profRCp;
#if defined alignment_debug || defined alignment_matrix_debug || defined rmap_debug || defined results_debug
      profdseqp = profiled_seqRCp;
#endif
    } else {
      scprofp = profp;
#if defined alignment_debug || defined alignment_matrix_debug || defined rmap_debug || defined results_debug
      profdseqp = profiled_seqp;
#endif
    }

    unprofiled_seqp = seqFastqGetConstSequence(sqbufp, &unprofiled_seqlen, NULL);
    scoreGetProfile(NULL, &profiled_seqlen, NULL, NULL, scprofp);

    if (unprofiled_seqlen > INT_MAX || unprofiled_seqlen != cp->re - cp->rs + 1 ||
	cp->qs > cp->qe || cp->qe >= profiled_seqlen)
      return ERRCODE_ASSERT;
#ifdef rmap_debug
    //printf("rmap_debug:Q::%s\n", qbasp);
    printf("rmap_debug:alignRMAPCANDFull:Q::[%i]", i);
    for (dp = profdseqp; (*dp); dp++)
      fputc(decoderp[(UCHAR) *dp], stdout);
    fputc('\n', stdout);

    //printf("rmap_debug:S::%s\n", seqp);
    printf("rmap_debug:alignRMAPCANDFull:S::[%i]", i);
    for (dp = unprofiled_seqp; *dp; dp++)
      fputc(decoderp[(UCHAR) *dp], stdout);
    fputc('\n', stdout);

    fprintRMAPCAND(stdout, cp);
#endif

    if (rmapflag & RMAPFLG_BEST) {
      resultSetGetMaxSwat(rsp, &swatscor_2ndmax);
      if (swatscor_2ndmax > min_swatscor)
	min_swatscor = swatscor_2ndmax;
    }
    aliRsltSetReset(alirsp);

    bw = cp->band_r - cp->band_l;
    if (bw < bandwidth_min) {
      bw = (bandwidth_min - bw + 1)/2;
      band_l = cp->band_l - bw;
      band_r = cp->band_r + bw;
    } else {
      band_l = cp->band_l;
      band_r = cp->band_r;
    }

    errcode = aliSmiWatInBand(alirsp, alibufp, scprofp, 
#if defined alignment_debug || defined alignment_matrix_debug
			      codecp,
			      profdseqp,
#endif
			      unprofiled_seqp, (int) unprofiled_seqlen,
			      band_l, band_r,
			      (int) cp->qs, (int) cp->qe,
			      0, (int) unprofiled_seqlen-1, 
			      min_swatscor, scorlen_min);
    
    if (!errcode) {
      errcode = resultSetAddFromAli(rsp, alirsp, cp->rs, 
				    0, profiled_seqlen,
				    (cp->sqidx == SEGCAND_UNKNOWN_SEQIDX)?
				    RESULTSET_UNKNOWN_SEQIDX: cp->sqidx,
#ifdef results_debug 
				    cp->segidx,
				    sqbufp,
				    profdseqp,
				    codecp,
				    ssp,
#endif
				    (char) (cp->flags & RMAPCANDFLG_REVERSE));
    }
    if (errcode)
      break;
  }

  return errcode;
}

/******************************************************************************
 ********************** Methods of Private Type RMAPPROF **********************
 ******************************************************************************/

static void deleteRMAPPROF(RMAPPROF *p)
{
  if (p != NULL) {
    scoreDeleteProfile(p->scorprofRCp);
    scoreDeleteProfile(p->scorprofp);
    seqFastqDelete(p->readRCp);
  }
  free(p);
}

static RMAPPROF *createRMAPPROF(const SeqCodec * const codecp)
/**< \param mode Combination of SCORE_PROFILE_MODES */
{
  RMAPPROF *p;
  const UCHAR mode = SCORPROF_SCALAR
#ifdef SCORE_SIMD
#ifdef SCORE_SIMD_SSE2
    | SCORPROF_STRIPED_8 | SCORPROF_STRIPED_16
#else
    | SCORPROF_STRIPED_32
#endif
#endif
    ;

  EMALLOCP0(p);
  if (NULL == p)
    return NULL;

  p->readRCp = seqFastqCreate(0, SEQTYP_FASTQ);
  p->scorprofp = scoreCreateProfile(0, codecp, mode);
  p->scorprofRCp = scoreCreateProfile(0, codecp, mode);
  
  if (NULL == p->readRCp ||
      NULL == p->scorprofp ||
      NULL == p->scorprofRCp) {
    deleteRMAPPROF(p);
    p = NULL;
  }
  
  return p;
}

static int makeRMAPPROFfromRead(RMAPPROF *prp, 
				const SeqFastq *readp, 
				const ScoreMatrix *scormtxp, 
				const SeqCodec *codecp)
{
  int errcode;

  seqFastqBlank(prp->readRCp);

  errcode = seqFastqAppendSegment(prp->readRCp, readp, 0, 0, 1, codecp);

  if (!(errcode)) {    
    errcode = scoreMakeProfileFromSequence(prp->scorprofp, readp, scormtxp);
  }
   
  if (!(errcode))
    errcode = scoreMakeProfileFromSequence(prp->scorprofRCp, prp->readRCp, scormtxp);

  return errcode;
}

/******************************************************************************
 ********************** Methods of Private Type RMAPINFO **********************
 ******************************************************************************/

static void deleteRMAPINFO(RMAPINFO *rmp)
{
  if (rmp) {
    hashDeleteHitInfo(rmp->hhiFp);
    hashDeleteHitInfo(rmp->hhiRp);
  }
  free(rmp);
}

static RMAPINFO *createRMAPINFO(const HashTable *htp)
{
  RMAPINFO *rmp;
  EMALLOCP0(rmp);
  if (!rmp) return 0;
  
  rmp->hhiFp = hashCreateHitInfo(HITINFO_BLKSZ, htp);
  rmp->hhiRp = hashCreateHitInfo(HITINFO_BLKSZ, htp);

  if (!((rmp->hhiFp) && (rmp->hhiRp))) {
    deleteRMAPINFO(rmp);
    rmp = 0;
  }
 
  return rmp;
}

static int initRMAPINFO(RMAPINFO *rmp,
			UCHAR min_basqval,
			SEQLEN_t seq_start,
			SEQLEN_t seq_end,
			const SeqFastq *readp, const HashTable *htp)
/** Resest Results and collect k-mer stats.
 * If seq_start + htp->ktup < seq_end, do this for a segment (seq_start, seq_end)
 * of the read only.
 */
{
  int errcode;

  errcode = hashCollectHitInfo(rmp->hhiFp, 0, min_basqval, seq_start, seq_end, readp, htp);
  if (!(errcode)) {
    errcode = hashCollectHitInfo(rmp->hhiRp, 1, min_basqval, seq_start, seq_end, readp, htp);
  }
  return errcode;
}

static int initRMAPINFOshort(RMAPINFO *rmp,
			     int maxhit_per_tuple,
			     UCHAR min_basqval,
			     const SeqFastq *readp, const HashTable *htp)
     /**< Initialise per-read data structures using a 'short' list of hits 
      * \param rmp Structure to be initialised.
      * \param maxhit_per_tuple Cut-off of the number of hits per kmer word. Kmer words
      *        with more than maxhit_per_tuple hits are skipped.
      *        If maxhit_per_tuple <= 0, do not apply a cut-off, but count all Kmer words.
      */
{
  int errcode;
  errcode = hashCollectHitInfoShort(rmp->hhiFp, 0, 
				    maxhit_per_tuple, HASH_MAXNHITS,
				    min_basqval,
				    readp, htp);
  if (!(errcode)) {
    errcode = hashCollectHitInfoShort(rmp->hhiRp, 1, 
				      maxhit_per_tuple, HASH_MAXNHITS,
				      min_basqval,
				      readp, htp);
  }
  return errcode;
}

#ifdef RMAP_SUPERFLUOUS_CODE
static void calcCoverDeficit(uint32_t deficit[NUM_STRANDS], const RMAPINFO *rmrp)
{
    deficit[0] = hashCalcHitInfoCoverDeficit(rmrp->hhiFp);
    deficit[1] = hashCalcHitInfoCoverDeficit(rmrp->hhiRp);
}
#endif

static uint32_t calcTotalNumberOfHits(const RMAPINFO *rmrp, int ktuple_maxhit)
{
  uint32_t hitnum = hashCalcHitInfoNumberOfHits(rmrp->hhiFp, ktuple_maxhit);
  return (uint32_t) (hitnum +
		     hashCalcHitInfoNumberOfHits(rmrp->hhiRp, ktuple_maxhit));
}

static uint32_t calcTotalHitNumStats(const RMAPINFO *rmrp, uint32_t *nhit_tot)
{
  uint32_t nhitF_rank, nhitR_rank;
  uint32_t nhitF = hashHitInfoCalcHitNumbers(rmrp->hhiFp, &nhitF_rank);
  uint32_t nhitR = hashHitInfoCalcHitNumbers(rmrp->hhiRp, &nhitR_rank);
  *nhit_tot = nhitF + nhitR;

  return nhitF_rank + nhitR_rank;
}
/******************************************************************************
 ********************** Methods of Private Type RMAPBUFF **********************
 ******************************************************************************/
static void deleteRMAPBUFF(RMAPBUFF *rmp)
{
  if (rmp) {
    seqFastqDelete(rmp->sqbfp);
    seqFastqDelete(rmp->qbfp);
    hashDeleteHitList(rmp->hhlp);
    segLstDelete(rmp->sglp);
    segAliCandsDelete(rmp->sacp);
    ARRDELETE(rmp->candr);
    aliBufferDelete(rmp->alibufp);
    aliRsltSetDelete(rmp->alirsltp);
    segQMaskDelete(rmp->qmp);
  }
  free(rmp);
}

static RMAPBUFF *createRMAPBUFF(const ScoreMatrix *scormtxp)
/** scormtxp is used if rmapflg & RMAPFLG_CMPLXW, is NULL otherwise */
{
  RMAPBUFF *rmp;
  EMALLOCP0(rmp);
  if (!rmp) return 0;

  rmp->sqbfp = seqFastqCreate(0, SEQTYP_FASTQ);
  rmp->qbfp = seqFastqCreate(0, SEQTYP_FASTQ);
  rmp->hhlp = hashCreateHitList(HASH_MAXNHITS);
  rmp->sglp = segLstCreate(SEGLST_BLOCKSIZ);
  rmp->sacp = segAliCandsCreate(SEGCANDS_BLOCKSIZ);
  ARRCREATE(rmp->candr, RMAPCANDS_BLOCKSIZ);
  rmp->alibufp = aliBufferCreate(0);
  rmp->alirsltp = aliRsltSetCreate(scormtxp, 0, 0, 0, 0);
  rmp->qmp = segQMaskCreate(0);

  if (!((rmp->hhlp) && (rmp->sglp) && (rmp->sacp) && (rmp->candr) &&
	(rmp->alibufp) && (rmp->alirsltp) && (rmp->qmp))) {
    deleteRMAPBUFF(rmp);
    rmp = 0;
  } 
  return rmp;
}

static void blankRMAPBUFF(RMAPBUFF *rmp)
{
  if (rmp) {
    seqFastqBlank(rmp->sqbfp);
    if ((rmp->qbfp)) 
      seqFastqBlank(rmp->qbfp);
    hashBlankHitList(rmp->hhlp);
    segLstBlank(rmp->sglp);
    segAliCandsBlank(rmp->sacp);
    ARRLEN(rmp->candr) = 0;
    aliRsltSetReset(rmp->alirsltp);
  }
}

static int fillRMAPBUFF(RMAPBUFF *bufp, const RMAPINFO *rmrp,
#ifdef RESULTS_TRACKER
			Track *trkp,
#endif
#ifdef hashhit_dump_sortarray
			FILE *dumpfp,
			int *dumpctr,
#endif
			BOOL with_seqidx, int ktuple_maxhit, 
			uint32_t min_ktup, uint32_t min_cover,
			const HashTable *htp,
			const SeqSet *ssp,
			const InterVal *ivr)
     /* Gnereate segments for type,
      * \param readp Query read, hast to be in SEQCOD_MANGLED encoding
      */
{
  int errcode;
#ifdef RESULTS_TRACKER
  TRACKFLG_t trk_flg;
  trackGetKtupData(trkp, NULL, NULL, &trk_flg);
#endif

  blankRMAPBUFF(bufp);
  
  if ((ivr)) {
    /* filtered */
    /* forward strand */
    errcode = collectHitsFromInterVal(bufp->sacp, bufp->sglp, bufp->hhlp, bufp->qmp,
#ifdef RESULTS_TRACKER
				      trkp,
#endif
				      ktuple_maxhit, min_ktup, min_cover,
				      rmrp->hhiFp, htp, ssp,
				      ivr);
    if (!(errcode)) {
      /* reverse strand */
      errcode = collectHitsFromInterVal(bufp->sacp, bufp->sglp, bufp->hhlp, bufp->qmp,
#ifdef RESULTS_TRACKER
					trkp,
#endif
					ktuple_maxhit, min_ktup, min_cover,
					rmrp->hhiRp, htp, ssp,
					ivr);
    }
  } else {
    errcode = collectHits(bufp->sacp, with_seqidx, 
			  bufp->sglp, bufp->hhlp, bufp->qmp,
#ifdef RESULTS_TRACKER
			  (trk_flg&TRACKFLG_REVERSE)? NULL:trkp,
#endif
#ifdef hashhit_dump_sortarray
			  dumpfp,
			  dumpctr,
#endif
			  ktuple_maxhit, min_ktup, min_cover,
			  rmrp->hhiFp, htp, ssp);
    if (!(errcode)) {
    /* reverse strand */
      errcode = collectHits(bufp->sacp, with_seqidx,
			    bufp->sglp, bufp->hhlp, bufp->qmp,
#ifdef RESULTS_TRACKER
			    (trk_flg&TRACKFLG_REVERSE)?trkp:NULL,
#endif
#ifdef hashhit_dump_sortarray
			    NULL,
			    NULL,
#endif
			    ktuple_maxhit, min_ktup, min_cover,
			    rmrp->hhiRp, htp, ssp);
    }
  }
  return errcode;
}

static int mapSingleRead(ErrMsg *errmsgp,
			 RMAPBUFF *bufp, 
			 ResultSet *rssp,
#ifdef RESULTS_TRACKER
			 Track *trkp,
#endif
#ifdef hashhit_dump_sortarray
			 FILE *dumpfp,
			 int *dumpctr,
#endif
			 const RMAPINFO *rmrp,
			 const RMAPPROF *rprofp,
			 SeqFastq *readp, 
			 int ktuple_maxhit,
			 uint32_t min_cover, 
			 int min_swatscor, int min_swatscor_below_max,
			 short target_depth, short max_depth, RMAPFLG_t rmapflg, 
			 const HashTable *htp, const SeqSet *ssp, const SeqCodec *codecp,
			 const InterVal *ivr)
{
  int errcode;
  char cod;
  uint32_t rlen, mincov_below_max;
#if defined alignment_matrix_debug || defined alignment_debug || defined rmap_debug || defined results_debug
  const char *qbasp = 0; 
  const char *qbasRCp = 0;
#endif
  UCHAR nskip;
  UCHAR ktup = hashTableGetKtupLen(htp, &nskip);
  int scorlen_min = ktup + nskip;    /* minimum length of aligned query segment */
  short mismatchscor, gapinitscor, gapextscor;
  short matchscor = scoreProfileGetAvgPenalties(&mismatchscor, &gapinitscor, &gapextscor, rprofp->scorprofp);
  short mismatchdiff = (short) (matchscor - mismatchscor);
  int bandwidth_min;
  uint32_t min_ktup = calcMinKtup(&min_cover, htp);
  SWATSCOR max1scor = 0, max2scor = 0, maxscor_perfect;
  uint32_t nseg, nseg_tot, nhit, nhit_tot;

  if (mismatchdiff < 0 || gapextscor >= 0 || mismatchscor >= 0) 
    ERRMSGNO(errmsgp, ERRCODE_ASSERT);

  /* mincov_below max translates the thresholds of the Smith-Waterman score
   * relative to the maximum (min_swatscor_below_max) into a relative threshold
   * of coverage */

  seqFastqGetConstSequence(readp, &rlen, &cod);
  if (rlen < ktup) {         /* consider output of a warning */
    return ERRCODE_SHORTSEQ;
  }

  if (rlen*matchscor > INT_MAX)
    ERRMSGNO(errmsgp, ERRCODE_OVERFLOW);

  maxscor_perfect = rlen*matchscor;

  if (min_swatscor_below_max < 0) {
    mincov_below_max = rlen-1;
  } else {
    mincov_below_max = ((uint32_t) (min_swatscor_below_max/mismatchdiff))*nskip;
  if (mincov_below_max < ktup || (rmapflg&RMAPFLG_BEST))
    mincov_below_max = ktup+2*(nskip-1);
  }

  /* set minimum Smith_Waterman score */
  /* if (min_swatscor < ktup*matchscor) min_swatscor = ktup*matchscor;  */
  
#ifdef rmap_debug
  printf("rmap_debug::rmapSingle(): ");
  printf("min_swatscor = %i, min_swatscor_below_max = %i\n", 
	 min_swatscor, min_swatscor_below_max);
#endif
  
  if ((errcode = fillRMAPBUFF(bufp, rmrp, 
#ifdef RESULTS_TRACKER
			      trkp,
#endif
#ifdef hashhit_dump_sortarray
			      dumpfp,
			      dumpctr,
#endif
			      (BOOL) (rmapflg&RMAPFLG_SEQBYSEQ), 
			      //(rmapflg&RMAPFLG_NOSHRTINFO)? ktuple_maxhit: 0,
			      ktuple_maxhit,
			      min_ktup, min_cover, htp, ssp, ivr)))
    
    ERRMSGNO(errmsgp, errcode);

  if ((errcode = aliBufferInit(bufp->alibufp, rlen)))
    ERRMSGNO(errmsgp, errcode);

  if ((errcode = segAliCandsStats(bufp->sacp,
#ifdef RESULTS_TRACKER
				  trkp,
#endif
				  mincov_below_max, 
				  rmrp->hhiFp,
				  rmrp->hhiRp, 
				  target_depth, max_depth, 
				  (uint8_t)(rmapflg & RMAPFLG_SENSITIVE))))
    ERRMSGNO(errmsgp, errcode);

#ifdef rmap_debug
  segAliCandsPrint(stdout, target_depth, bufp->sacp);
#endif

  nseg = segAliCandsGetNumberOfSegments(bufp->sacp, NULL, NULL, NULL, NULL, &nseg_tot);
  if (nseg > INT_MAX || nseg_tot > INT_MAX)
    ERRMSGNO(errmsgp, ERRCODE_ASSERT);

  nhit = calcTotalHitNumStats(rmrp, &nhit_tot);
  resultSetAlignmentStats(rssp, (int) nseg, (int) nseg_tot, max_depth, nhit, nhit_tot);

#ifdef rmap_mscor_calib
  printf("CALDAT_RMAP %s %u %u\n", seqFastqGetSeqName(readp), nseg, nseg_tot);
#endif
  if (cod == SEQCOD_ASCII &&
      ((errcode = seqFastqEncode(readp, codecp))))
    ERRMSGNO(errmsgp, errcode);

#if defined alignment_matrix_debug || defined alignment_debug || defined rmap_debug || defined results_debug
  qbasp = seqFastqGetConstSequence(readp, NULL, NULL);
  qbasRCp = seqFastqGetConstSequence(rprofp->readRCp, NULL, NULL);
#endif

  errcode = scoreRMAPCAND(&bufp->candr, 
			  &max1scor,
			  &max2scor,
			  bufp->sqbfp,
			  bufp->alibufp,
			  rmapflg,
			  nskip,
#if defined alignment_debug || defined alignment_matrix_debug || defined rmap_debug
			  qbasp,
			  qbasRCp,
#endif
			  rprofp->scorprofp,
			  rprofp->scorprofRCp,
			  ssp, 
			  codecp,
			  bufp->sacp);
    
  if ((errcode) && 
      ((rmapflg&RMAPFLG_SEQBYSEQ) || errcode != ERRCODE_SEQOFFS))
    ERRMSGNO(errmsgp, errcode);

  if (max1scor > maxscor_perfect)
    ERRMSGNO(errmsgp, ERRCODE_ASSERT);

  if (max1scor < 1)
    return ERRCODE_SUCCESS;

  bandwidth_min = (maxscor_perfect - max1scor)/(-1*gapextscor);
 
  if (min_swatscor_below_max >= max1scor)
    min_swatscor_below_max = max1scor;

  if (min_swatscor > max2scor && max2scor > 0)
    min_swatscor = max2scor;
 
  if (min_swatscor_below_max >= 0) {
     SWATSCOR minswc = (max2scor > 0)? max2scor: max1scor;
     if ((rmapflg&RMAPFLG_BEST)) {
       if (minswc > min_swatscor)
	 min_swatscor = minswc;
     } else if (min_swatscor + min_swatscor_below_max < max1scor) {
       min_swatscor =  max1scor - min_swatscor_below_max;
       if (min_swatscor > minswc)
	 min_swatscor = minswc;
     }
  }
  
  if (min_swatscor > scorlen_min*matchscor && matchscor > 0)
    scorlen_min = min_swatscor/matchscor;

  if ((errcode = alignRMAPCANDFull(rssp, 
				   bufp->alirsltp, 
				   bufp->alibufp, 
				   bufp->sqbfp,
				   min_swatscor,
				   scorlen_min,
				   bandwidth_min,
				   rmapflg,
				   rprofp->scorprofp,
				   rprofp->scorprofRCp, 
#if defined alignment_matrix_debug || defined alignment_debug || defined rmap_debug || defined results_debug
				   qbasp,
				   qbasRCp,
#endif
				   ssp, codecp, bufp->candr)))
    ERRMSGNO(errmsgp, errcode);
  
  errcode = resultSetSortAndAssignSequence(rssp,
					   bufp->sqbfp, 
					   0, 
					   readp,
#ifdef results_debug 
					   rprofp->readRCp,
#endif
					   rprofp->scorprofp, 
					   rprofp->scorprofRCp,
					   ssp, codecp);
  if (errcode) 
    ERRMSGNO(errmsgp, errcode);
  
  return errcode;
}

static int mapSecondary(ErrMsg *errmsgp,
			RMAPBUFF *bufp, ResultSet *rssp, 
#ifdef RESULTS_TRACKER
			Track *trkp,
#endif
			RMAPINFO *rmrp,
			const RMAPPROF *rmprp,
			SeqFastq *readp, 
			int ktuple_maxhit,
			uint32_t min_cover,
			int min_swatscor, int min_swatscor_below_max,
			UCHAR min_basqval,
			short target_depth, short max_depth, RMAPFLG_t rmapflg, 
			const HashTable *htp, 
			const SeqSet *ssp,
			const SeqCodec *codecp)
{
  int errcode;
  SEQLEN_t qs, qe, qlen;
  UCHAR ktup, nskip;
  const Result *rp;

  ktup = hashTableGetKtupLen(htp, &nskip);
  seqFastqGetConstSequence(readp, &qlen, NULL);
  errcode = resultSetGetResultInSegment(&rp, 0, 0, rssp);
  if (ERRCODE_SUCCESS != errcode) {
    if (ERRCODE_FAILURE == errcode) {
      return ERRCODE_SUCCESS; /* no alignment found, keep silent */
    } else 
      ERRMSGNO(errmsgp, errcode);
  }

  if ((errcode = resultGetData(&qs, &qe, NULL, NULL, NULL, NULL, NULL, rp)))
    ERRMSGNO(errmsgp, errcode);
  
  if (qe > qlen || qs > qe)
    ERRMSGNO(errmsgp, ERRCODE_ASSERT);

  if (qs + qe > qlen) {
    qe = (qs > 1)? qs - 2: 0;
    qs = 0;
  } else {
    qs =  qe;
    qe = qlen - 1;
  }
  if (qs + ktup + nskip > qe + 1)
    return ERRCODE_SUCCESS;

  errcode = initRMAPINFO(rmrp, min_basqval, qs, qe, readp, htp);
  if (errcode == ERRCODE_SHORTSEQ)
    return ERRCODE_SUCCESS;

  if (errcode)
    ERRMSGNO(errmsgp, errcode);

  mapSingleRead(errmsgp,
		bufp, rssp,
#ifdef RESULTS_TRACKER
		trkp,
#endif
#ifdef hashhit_dump_sortarray
		NULL,
		NULL,
#endif
		rmrp, rmprp, readp, ktuple_maxhit, min_cover, 
		min_swatscor, min_swatscor_below_max,
		target_depth, max_depth, rmapflg, 
		htp, ssp, codecp, NULL);

  return errcode;
}

/******************************************************************************
 ************************* Public Methods of Type RMap ************************
 ******************************************************************************/

RMap *rmapCreate(const HashTable *htp, const SeqCodec *codecp, 
		 const SeqSet *ssp,
		 const ScoreMatrix *scormtxp, RMAPFLG_t rmapflg)
{
  BOOL okflg = 0;
  SETSIZ_t totlen;
  RMap *rmp;
  EMALLOCP0(rmp);
  if (!rmp) return 0;
  
  seqSetGetSeqNumAndTotLen(&totlen, ssp);

  rmp->bfp = createRMAPBUFF((rmapflg & RMAPFLG_CMPLXW)? scormtxp: NULL);
  rmp->prp = createRMAPPROF(codecp);
  rmp->mrp = createRMAPINFO(htp);
  rmp->rsrp = resultSetCreate(0,0);
  okflg = (BOOL) (rmp->bfp != NULL && 
		  rmp->prp != NULL &&
		  rmp->mrp != NULL && 
		  rmp->rsrp != NULL);
  if ((okflg) && (rmapflg & RMAPFLG_PAIRED)) {
#ifdef rmap_finehash_2ndmate
    uint8_t nskip = FINEHASH_SKIPSTEP;
#endif
    rmp->pmp = createRMAPPROF(codecp);
    rmp->mmp = createRMAPINFO(htp);
    rmp->rsmp = resultSetCreate(0,0);
    rmp->ivr = interValCreate(FILTERIVAL_BLKSZ);
    rmp->pairp = resultSetCreatePairs(0);
#ifdef rmap_finehash_2ndmate
    if ((totlen + 1)/nskip >  HASHPOS_MAX) {
      SETSIZ_t ns = (totlen+1)/HASHPOS_MAX + 1;
      if (ns > UINT8_MAX)
	okflg = 0;
      nskip = (uint8_t) ns;
    }
    if (okflg) {
      rmp->htflyp = hashTableCreate(FINEHASH_WORDLEN, nskip, 0, 0, HASHIDXTYP_PERFECT);
      if (rmp->htflyp)
	rmp->mflyp = createRMAPINFO(rmp->htflyp);
    
#endif
      okflg = (BOOL) ((rmp->mmp != NULL) && 
		      (rmp->rsmp != NULL) && 
		      (rmp->ivr != NULL) &&
#ifdef rmap_finehash_2ndmate
		      (rmp->htflyp != NULL) &&
#endif
		      (rmp->pairp != NULL)); 
#ifdef rmap_finehash_2ndmate 
    }
#endif
  } else {
    rmp->pmp = NULL;
    rmp->mmp = NULL;
    rmp->rsmp = NULL;
    rmp->ivr = NULL;
    rmp->pairp = NULL;
#ifdef rmap_finehash_2ndmate
    rmp->htflyp = NULL;
    rmp->mflyp = NULL;
#endif
  }

  if ((okflg) && (rmapflg&RMAPFLG_SPLIT)) {
    rmp->mr2p = createRMAPINFO(htp);
    okflg = (BOOL) (NULL != rmp->mr2p);
    if ((okflg) && (rmapflg & RMAPFLG_PAIRED)) {
      rmp->mm2p = createRMAPINFO(htp);
      okflg = (BOOL) (NULL != rmp->mm2p);
    } else {
      rmp->mm2p = NULL;
    }
  } else {
    rmp->mr2p = NULL;
    rmp->mm2p = NULL;
  }

  if (!(okflg)) {
    rmapDelete(rmp);
    rmp = 0;
  }

  return rmp;
}

void rmapDelete(RMap *rmp)
{
  if (rmp) {
    deleteRMAPBUFF(rmp->bfp);
    deleteRMAPPROF(rmp->prp);
    deleteRMAPPROF(rmp->pmp);
    deleteRMAPINFO(rmp->mrp);
    deleteRMAPINFO(rmp->mr2p);
    deleteRMAPINFO(rmp->mmp);
    deleteRMAPINFO(rmp->mm2p);
#ifdef rmap_finehash_2ndmate
    hashTableDelete(rmp->htflyp);
    deleteRMAPINFO(rmp->mflyp);
#endif
    interValDelete(rmp->ivr);
    resultSetDelete(rmp->rsrp);
    resultSetDelete(rmp->rsmp);
    resultSetDeletePairs(rmp->pairp);
  }

  free(rmp);
}

void rmapBlank(RMap *rmp)
{
  if (rmp) {
    blankRMAPBUFF(rmp->bfp);
    resultSetBlank(rmp->rsrp);
    if (rmp->rsmp != NULL)
      resultSetBlank(rmp->rsmp);
    if (rmp->ivr != NULL) 
      interValBlank(rmp->ivr);
    if (rmp->pairp != NULL)
      resultSetBlankPairs(rmp->pairp);
  }
}

void rmapGetData(const ResultSet **rslt_readp, 
		 const ResultSet **rslt_matep, 
		 const ResultPairs **pairp, 
		 SeqFastq **sbufAp, 
		 SeqFastq **sbufBp,
		 const RMap *rmp)
{
  if (rslt_readp) *rslt_readp = rmp->rsrp;
  if (rslt_matep) *rslt_matep = rmp->rsmp;
  if (pairp) *pairp = rmp->pairp;
  if (sbufAp) *sbufAp = rmp->bfp->sqbfp;
  if (sbufBp) *sbufBp = rmp->bfp->qbfp;
}

int rmapSingle(ErrMsg *errmsgp,
	       RMap *rmp, 
#ifdef hashhit_dump_sortarray
	       FILE *dumpfp,
	       int *dumpctr,
#endif
	       SeqFastq *readp, 
#ifdef RESULTS_TRACKER
	       Track *trkrp,
#endif
	       int ktuple_maxhit,
	       uint32_t min_cover, int min_swatscor, int min_swatscor_below_max,
	       UCHAR min_basqval,
	       short target_depth, short max_depth,
	       RMAPFLG_t rmapflg,
	       const ScoreMatrix *scormtxp,
	       const ResultFilter *rsfp, 
	       const HashTable *htp, 
	       const SeqSet *ssp, 
	       const SeqCodec *codecp)
{
  int errcode;

  rmapBlank(rmp);

#ifdef RESULTS_TRACKER
  if ((errcode = trackMakeFromSequence(trkrp, readp, 
				       RESULTS_TRACKER_IS1BASED, 
				       ssp, htp)))
    ERRMSGNO(errmsgp, errcode);
#endif

  if ((errcode = makeRMAPPROFfromRead(rmp->prp, readp, scormtxp, codecp)))
    ERRMSGNO(errmsgp, errcode);

#ifdef rmap_short_hitinfo
  if (rmapflg & RMAPFLG_NOSHRTINFO) {
#endif
    errcode = initRMAPINFO(rmp->mrp, min_basqval, 0, 0, readp, htp);
#ifdef rmap_short_hitinfo
  } else {
    errcode = initRMAPINFOshort(rmp->mrp, 
				ktuple_maxhit, 
				min_basqval,
				readp, htp);
  }
#endif
  if ((errcode) && errcode != ERRCODE_SHORTSEQ) 
    ERRMSGNO(errmsgp, errcode);

  if (!errcode) {
    mapSingleRead(errmsgp,
		  rmp->bfp, 
		  rmp->rsrp, 
#ifdef RESULTS_TRACKER
		  trkrp,
#endif
#ifdef hashhit_dump_sortarray
		  dumpfp,
		  dumpctr,
#endif
		  rmp->mrp, rmp->prp, readp, 
		  ktuple_maxhit, min_cover, 
		  min_swatscor, min_swatscor_below_max,
		  target_depth, max_depth, rmapflg, 
		  htp, ssp, codecp, NULL);
  }

  if (!(errcode) && (rmapflg & RMAPFLG_SPLIT)) {
    mapSecondary(errmsgp,
		 rmp->bfp, 
		 rmp->rsrp,
#ifdef RESULTS_TRACKER
		 trkrp,
#endif
		 rmp->mr2p, rmp->prp, readp, 
		 ktuple_maxhit, min_cover, 
		 min_swatscor, min_swatscor_below_max, min_basqval,
		 target_depth, max_depth, rmapflg, 
		 htp, ssp, codecp);
  }

  if ((errcode) && errcode != ERRCODE_SHORTSEQ) 
    ERRMSGNO(errmsgp, errcode);

  if (!errcode) 
    errcode = resultSetFilterResults(rmp->rsrp, rsfp, readp);

  if (errcode == ERRCODE_SHORTSEQ)
    errcode = ERRCODE_SUCCESS;
  else if (errcode)
    ERRMSGNO(errmsgp, errcode);

  return errcode;
}

int rmapPair(ErrMsg *errmsgp,
	     RMap *rmp, 
	     SeqFastq *readp, SeqFastq *matep,
#ifdef RESULTS_TRACKER
	     Track *trkrp,
	     Track *trkmp,
#endif
	     RSLTPAIRFLG_t *pairflgp,
	     int d_min, int d_max,
	     RSLTPAIRLIB_t pairlibcode,
	     int ktuple_maxhit, uint32_t mincov_read, uint32_t mincov_mate,
	     int min_swatscor,
	     UCHAR min_basqval,
	     short target_depth, short max_depth,
	     RMAPFLG_t rmapflg,
	     const ScoreMatrix *scormtxp,
	     const ResultFilter *rsfp,
	     const HashTable *htp, const SeqSet *ssp, const SeqCodec *codecp)
{
  int errcode, errcode_read, errcode_mate, n_proper, mapq1, swscor1, swscor2_restricted;
  UCHAR rare_mate;
  uint32_t mincov1, mincov2, nhit_read, nhit_mate;
  RMAPBUFF *bufp = rmp->bfp;
  RMAPINFO *rr1p, *rr2p, *rmrp = rmp->mrp, *rmmp = rmp->mmp;
  RMAPPROF *rp1p, *rp2p;
  ResultSet *rs1p, *rs2p, *rsrp = rmp->rsrp, *rsmp = rmp->rsmp;
  SeqFastq *read1p, *read2p;
#ifdef RESULTS_TRACKER
  Track *trk1p, *trk2p;
#endif


#ifdef rmap_debug
  printf("rmapPair: d_min=%i, dmax=%i, ktuple_maxhit=%i,\n",
	 d_min, d_max, ktuple_maxhit);
  printf("          mincov_read=%u, mincov_mate=%u, min_swatscor=%i,\n",
	 mincov_read, mincov_mate, min_swatscor); 
  printf("          target_depth=%hi, max_depth=%hi, rmapflg=%i\n",
	 target_depth, max_depth, (int) rmapflg);   
#endif
  if (!((rmmp) && (rmp->ivr) && (rmp->pairp)))
    ERRMSGNO(errmsgp, ERRCODE_ASSERT);

  rmapBlank(rmp);
  *pairflgp = RSLTPAIRFLG_PAIRED;  

#ifdef RESULTS_TRACKER
  if ((errcode = trackMakeFromSequence(trkrp, readp, RESULTS_TRACKER_IS1BASED,
				       ssp, htp)))
    ERRMSGNO(errmsgp, errcode);

  if ((errcode = trackMakeFromSequence(trkmp, matep, RESULTS_TRACKER_IS1BASED,
				       ssp, htp)))
     ERRMSGNO(errmsgp, errcode);
#endif

  if ((errcode = makeRMAPPROFfromRead(rmp->prp, readp, scormtxp, codecp)))
    ERRMSGNO(errmsgp, errcode);

  if ((errcode = makeRMAPPROFfromRead(rmp->pmp, matep, scormtxp, codecp)))
    ERRMSGNO(errmsgp, errcode);

#ifdef rmap_short_hitinfo
  if (rmapflg &  RMAPFLG_NOSHRTINFO) {
#endif
    if ((errcode_read = initRMAPINFO(rmrp, min_basqval, 0, 0, readp, htp)) &&
	errcode_read != ERRCODE_SHORTSEQ)
      ERRMSGNO(errmsgp, errcode_read);

    if ((errcode_mate = initRMAPINFO(rmmp, min_basqval, 0, 0, matep, htp)) &&
	errcode_mate != ERRCODE_SHORTSEQ)
      ERRMSGNO(errmsgp, errcode_mate);
#ifdef rmap_short_hitinfo
  } else {
    if ((errcode_read = initRMAPINFOshort(rmrp, 
					  ktuple_maxhit, 
					  min_basqval,
					  readp, htp)) &&
	errcode_read != ERRCODE_SHORTSEQ)
      ERRMSGNO(errmsgp, errcode_read);

    if ((errcode_mate = initRMAPINFOshort(rmmp,
					  ktuple_maxhit,
					  min_basqval,
					  matep, htp)) &&
	errcode_mate != ERRCODE_SHORTSEQ)
      ERRMSGNO(errmsgp, errcode_mate);
  }
#endif
  if ((errcode_read) && (errcode_mate))
    return ERRCODE_SUCCESS;

  if ((errcode_read)) {
    mapSingleRead(errmsgp,
		  bufp,
		  rmp->rsmp,
#ifdef RESULTS_TRACKER
		  trkmp,
#endif
#ifdef hashhit_dump_sortarray
		  NULL,
		  NULL,
#endif
		  rmmp, rmp->pmp, matep, ktuple_maxhit, mincov_mate,
		  min_swatscor, MINSCOR_BELOW_MAX_BEST,
		  target_depth, max_depth, rmapflg, 
		  htp, ssp, codecp, NULL);
  }
  
  if ((errcode_mate)) {
    mapSingleRead(errmsgp,
		  bufp,
		  rmp->rsrp,
#ifdef RESULTS_TRACKER
		  trkrp,
#endif
#ifdef hashhit_dump_sortarray
		  NULL,
		  NULL,
#endif
		  rmrp, rmp->prp, readp, ktuple_maxhit, mincov_read,
		  min_swatscor, MINSCOR_BELOW_MAX_BEST,
		  target_depth, max_depth, rmapflg, 
		  htp, ssp, codecp, NULL);
  }
 
  nhit_read = calcTotalNumberOfHits(rmrp, ktuple_maxhit);
  nhit_mate = calcTotalNumberOfHits(rmmp, ktuple_maxhit);
  if (nhit_read > nhit_mate) {
    *pairflgp |= RSLTPAIRFLG_RAREMATE;
    rare_mate = 1;
    read1p = matep;
    read2p = readp;
    rr1p = rmmp;
    rr2p = rmrp;
    rs1p = rsmp;
    rs2p = rsrp;
    rp1p = rmp->pmp;
    rp2p = rmp->prp;
#ifdef RESULTS_TRACKER
    trk1p = trkmp;
    trk2p = trkrp;
#endif
    mincov1 = mincov_mate;
    mincov2 = mincov_read;
  } else {
    rare_mate = 0;
    read1p = readp;
    read2p = matep;
    rr1p = rmrp;
    rr2p = rmmp;
    rs1p = rsrp;
    rs2p = rsmp;
    rp1p = rmp->prp;
    rp2p = rmp->pmp;
#ifdef RESULTS_TRACKER
    trk1p = trkrp;
    trk2p = trkmp;
#endif
    mincov1 = mincov_read;
    mincov2 = mincov_mate;
  }

  mapSingleRead(errmsgp,
		bufp, 
		rs1p,
#ifdef RESULTS_TRACKER
		trk1p,
#endif
#ifdef hashhit_dump_sortarray
		NULL,
		NULL,
#endif
		rr1p, rp1p, read1p, ktuple_maxhit, mincov1,
		min_swatscor, MINSCOR_BELOW_MAX_BEST,
		target_depth, max_depth, rmapflg, 
		htp, ssp, codecp, NULL);
  
  n_proper = 0;
  mapq1 = resultSetGetMappingScore(rs1p, &swscor1);
  swscor2_restricted = 0;
  /* if (resultSetGetMappingScore(rs1p, &swscor1) >= MAPSCORE_UNIQUE_MAPPED_1ST && */
  /*     !(rmapflg & RMAPFLG_ALLPAIR)) { */
  //if (mapq1 >= MAPSCORE_UNIQUE_MAPPED_1ST && !(rmapflg & RMAPFLG_ALLPAIR)) {
    /* rare mate gives rise to a unique hit. use this to restrict search for
     * 2nd mate */
  /* use 1st mate to restrict 2nd mate */
  /* if (mapq1 >= MAPSCORE_UNIQUE_MAPPED_1ST || */
  /*     !resultSetAlignmentWasCurtailed(rs1p)) { */

  errcode = setupInterValFromResultSet(rmp->ivr, d_min, d_max, read1p, read2p,
				       htp, ssp, rs1p);
  if ((errcode)) 
    ERRMSGNO(errmsgp, errcode);
  interValPrune(rmp->ivr);

  mapSingleRead(errmsgp,
		bufp, 
		rs2p, 
#ifdef RESULTS_TRACKER
		trk2p,
#endif
#ifdef hashhit_dump_sortarray
		NULL,
		NULL,
#endif
		rr2p, rp2p,
		read2p, ktuple_maxhit, mincov2,
		min_swatscor, MINSCOR_BELOW_MAX_BEST,
		target_depth, max_depth, rmapflg,
		htp, ssp, codecp, rmp->ivr);

  errcode = resultSetFindProperPairs(rmp->pairp, d_min, d_max,
				     MAXNUM_PAIRS_TOTAL, 
				     0, pairlibcode,
				     rsrp, rsmp);
  if ((errcode) && errcode != ERRCODE_PAIRNUM)
    ERRMSGNO(errmsgp, errcode);
  resultSetGetMappingScore(rs2p, &swscor2_restricted);
  resultSetGetNumberOfPairs(&n_proper, rmp->pairp);

  if ((rmapflg & RMAPFLG_ALLPAIR) ||
      n_proper < 1 ||
      mapq1 < MAPSCORE_UNIQUE_MAPPED_1ST ||
      //resultSetAlignmentWasCurtailed(rs1p) ||
      !scorIsAboveFractMax(swscor2_restricted, swscor1, MINFRACT_MAXSCOR_2ND, read2p, read1p)) {
    int mapq2;
    int swscor2;
    /* no proper pairs or not confident of 1st mate -> 
     * unrestricted mapping of 2nd mate */
    if (n_proper < 1)
      resultSetBlank(rs2p);
    mapSingleRead(errmsgp,
		  bufp, rs2p, 
#ifdef RESULTS_TRACKER
		  trk2p,
#endif
#ifdef hashhit_dump_sortarray
		  NULL,
		  NULL,
#endif
		  rr2p, rp2p,
		  read2p, ktuple_maxhit, mincov2,
		  min_swatscor, MINSCOR_BELOW_MAX_BEST,
		  target_depth, max_depth, rmapflg,
		  htp, ssp, codecp, NULL);

    /* if best alignment score > than alignment score of first mate,
     * perform restricted mapping with first mate to see whether anything
     * better was missed within boundaries. */
    mapq2 = resultSetGetMappingScore(rs2p, &swscor2);
    if (mapq2 > MAPSCORE_UNIQUE_MAPPED_1ST ||
	 //n_proper < 1 ||
	swscor2 > swscor2_restricted ||
	swscor2 > swscor1) {
      int swscor1_2ndbest = 0;
#ifdef rmap_finehash_2ndmate
      SEQLEN_t rlen = 0;
#endif
      resultSetGetScorStats(rs1p, NULL, NULL, &swscor1_2ndbest, NULL);

      if ((errcode = setupInterValFromResultSet(rmp->ivr, d_min, d_max, read2p, read1p,
						htp, ssp, rs2p)))
	ERRMSGNO(errmsgp, errcode);
      interValPrune(rmp->ivr);

#ifdef rmap_finehash_2ndmate
      seqFastqGetConstSequence(read1p, &rlen, NULL);
      if (hashTableGetKtupLen(htp, NULL) <= rlen) {
	errcode = setupFineHashTable(rmp->htflyp, bufp->sqbfp, 
				     ssp, rmp->ivr, htp, codecp);
	if ((errcode) && errcode != ERRCODE_MAXKPOS)
	  ERRMSGNO(errmsgp, errcode);

	if (!(errcode)) {

	  if ((errcode = initRMAPINFO(rmp->mflyp, min_basqval, 0, 0, 
				    read1p, rmp->htflyp)) &&
	      errcode != ERRCODE_SHORTSEQ)
	    ERRMSGNO(errmsgp, errcode);

	  mapSingleRead(errmsgp,
			bufp, 
			rs1p, 
#ifdef RESULTS_TRACKER
			trk1p,
#endif
#ifdef hashhit_dump_sortarray
			NULL,
			NULL,
#endif
			rmp->mflyp, 
			rp1p, read1p, ktuple_maxhit, mincov1,
			swscor1_2ndbest, MINSCOR_BELOW_MAX_BEST,
			target_depth, max_depth, rmapflg, 
			rmp->htflyp, ssp, codecp, rmp->ivr);
	} else {
	  errcode = ERRCODE_SUCCESS;
#endif //#ifdef rmap_finehash_2ndmate
	  mapSingleRead(errmsgp,
			bufp, 
			rs1p, 
#ifdef RESULTS_TRACKER
			trk1p,
#endif
#ifdef hashhit_dump_sortarray
			NULL,
			NULL,
#endif
			rr1p,
			rp1p, read1p, ktuple_maxhit, mincov1,
			swscor1_2ndbest, MINSCOR_BELOW_MAX_BEST,
			target_depth, max_depth, rmapflg, 
			htp, ssp, codecp, rmp->ivr);
#ifdef rmap_finehash_2ndmate
	}
      }
#endif 
    }
  } else {
    *pairflgp = (RSLTPAIRFLG_t)((*pairflgp) | 
				((rare_mate == 0)? 
				 RSLTPAIRFLG_RESTRICT_2nd: 
				 RSLTPAIRFLG_RESTRICT_1st));
  }

  if (errcode)
    ERRMSGNO(errmsgp, errcode);

  if (rmapflg & RMAPFLG_SPLIT) {
    mapSecondary(errmsgp,
		 bufp, 
		 rmp->rsrp,
#ifdef RESULTS_TRACKER
		 trkrp,
#endif
		 rmp->mr2p,
		 rmp->prp, readp, ktuple_maxhit, mincov_read, 
		 min_swatscor, MINSCOR_BELOW_MAX_BEST, min_basqval,
		 target_depth, max_depth, rmapflg, 
		 htp, ssp, codecp);

    mapSecondary(errmsgp,
		 bufp, 
		 rmp->rsmp,
#ifdef RESULTS_TRACKER
		 trkmp,
#endif 
		 rmp->mm2p,
		 rmp->pmp, matep, ktuple_maxhit, mincov_mate, 
		 min_swatscor, MINSCOR_BELOW_MAX_BEST, min_basqval,
		 target_depth, max_depth, rmapflg, 
		 htp, ssp, codecp);
  }
  
  errcode = resultSetFindPairs(rmp->pairp, *pairflgp, pairlibcode, 
			       d_min, d_max, 
			       rsrp, rsmp);
  if ((errcode))
    ERRMSGNO(errmsgp, errcode);

  if ((errcode = resultSetFilterResults(rsrp, rsfp, readp)))
    ERRMSGNO(errmsgp, errcode);

  if ((errcode = resultSetFilterResults(rsmp, rsfp, matep)))
    ERRMSGNO(errmsgp, errcode);

  return errcode;
}

#ifdef RESULTS_BASQUAL
int rmapPrintObservedBasQ(const RMap *rmp)
{
  return resultSetPrintBasqStats(stdout, rmp->rsrp, (rmp->rsmp == NULL)? NULL: rmp->rsmp);
}
#endif
