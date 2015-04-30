/** Define pairwise segments as candidates for dynamic programming */

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

/* Concept:
 * Overlapping k-tuple hits with identical shift define a SEED.
 * (SEEDs therfore don't overlap).
 *
 * Multiple SEEDs of constant shift define a (constant-shift) SEGMENT.
 * A SEGMENT of more than 1 SEED therefore may contain mismatches
 * between its SEEDS (or indels if insertions and deletions cancel
 * out). SEGMENTs may overlap.  HITREGIONs are defined by maximum
 * shift differences between successive k-tuple hits.
 * 
 * Candidate segments for alignment (SegCand) are defined by segment
 * boundaries in the query and reference sequences and by a shift
 * band, specified as an start shift and a shift range.
 * 
 * HITREGIONs, constant-shift SEGMENTs and SEEDs are not necessary for
 * the simplest algorithm defining candidates for alignment. Depending
 * on sequencing platform and read length they may incur quite a bit
 * of time and memory overhead. Nevertheless, they are helpful as
 * concepts and I am defining them here in order to experiment and
 * gather statistics. Once an optimal strategy crystallizes some or
 * all of them can be removed.
 *
 * HIREGIONs might be helpful during the alignment stage. In cases
 * where there are multiple overlapping alignments one would fetch the
 * segment corresponding to that region only once. They are also used
 * as a filter discarding regions with sparse k-tuple hits (i.e. where
 * the minimum coverage could never be reached by the number of
 * k-tuple hits)
 *
 * SEEDs are useful where there are long exactly matching stretches
 * (Illumina & 454 and next-gen sequencing, but probably not
 * caplillary reads).
 *
 * SEGMENTS are a useful concept where SEEDs are useful and where the
 * strategy for defining candidates is a bit more involved. For
 * example, in order to avoid Smith-Waterman, one might try to
 * directly match a SEGMENT first where the coverage by SEED is so
 * high that max 2 mismatches are possible.
 * 
 * Strategies for defining candidate segments for alignment:
 * 
 * (i) The simplest algorithm would scan along the k-tuple hits (which
 * had been sorted first by shift, 'diagonal search', and then by
 * offset in query read) until the shift difference btween 2
 * successive hits is greater than a cutoff (or grate than the read
 * length), and establish extrema in query and reference offsets as
 * well as in the shift. This would define 1 candidate for alignment
 * per HITREGION. SEEDs, SEGMENTs (and HITREGIONS) would not be
 * needed, but in regions with possible overlapping alignments, only 1
 * alignment would be found. Unless all k-tuples had the same shift,
 * Smith-Waterman would have to be carried out. A possible
 * modification would interrupt the scan where the coverage by k-tuple
 * hits is sufficiently high.
 *
 * (ii) In order to deal with overlapping alignments, we first define
 * SEEDs and SEGMENTs. Then we scan the segments along (increasing)
 * shift and join them until the coverage reaches a threshold and the
 * increase in coverage contributed by the new segment is negligible.
 * The idea is to separate e.g. overlapping segments of sufficiently
 * high coverage into separate Smith-Waterman alignments with narrower
 * bands.
 *
 * (iii) One could start scanning for segments that could potentially
 * be aligned without Smith-Waterman. Where that fails one could exend
 * the segment by joining neigbouring segments (left and right along
 * shift axis) of low coverage until contributed coverage is
 * negligible. Those segments would define a band for
 * Smith-Waterman. Once all those candidates are defined one would
 * attempt to join the remaining SEGMENTs.
 *
 * For the time being I will adopt strategy (ii) and label candidates
 * that could promisingly be aligned without Smith-Waterman as well as
 * remember the shift of that segment.
 *
 * All this may turn out to be overkill and strategy (i) might be
 * quickest. The SSAHA2 software adopts strategy (i) but with
 * ktuple-hit numbers.  This is much faster but doesn't separate
 * overlaps well (important e.g. for P. falciparum). Exact matches can
 * be identified but not nearly exact matches containing 1 or 2
 * mismatches.
 */

#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "elib.h"
#include "array.h"
#include "sort.h"
#include "segment.h"

enum {
  DEFAULT_BLKSZ = 16*4096,
  DEFAULT_QBLKSZ = 1024,   
  /**< Default block size for query mask */
  SEGMENTING_DIFFSHIFT = 3, 
  /**< Maximum shift difference between successive hits
   * in a segment specified in units of ktuple lengths */
  SEEDLEN_FLAGGED = -1,
  /**< value for (struct seed).len that flags out the seed */
  NUMBER_OF_STRANDS = 2,
  /**< For looping over forwared [0] and reverse [1] strands */
  MAXIMUM_DEPTH = 8000, 
  /**< Maximum number of candidate segments allowed for alignment */
  DEFAULT_TARGET_DEPTH = 200, 
  /**< Target depth if none specified */
  EDGE_BAND_FACTOR = 4, 
  /**< divide the length of the largest stretch of query sequence
   * not coverered by seeds by EDGE_BAND_FACTOR to obtain an added 
   * band_width */
  MAX_BANDEDGE_2POW = 4, 
  /**< devide read length by 2^MAX_BANDEDGE_2POW to obtain a maximum
   * added band width */
};

enum SEGLST_FLAGS { /**< Status flags for type SegLst */
  SEGLSTFLG_REVERSE = 0x01,
  SEGLSTFLG_HITREGIONS = 0x02, /**< HITREGIONs for seeds are determined */
  SEGLSTFLG_SEEDS = 0x04,      /**< SEEDs are loaded */
  SEGLSTFLG_SEGMENTS = 0x08,   /**< SEGMENTs are loaded */
};

typedef uint8_t UCHAR;
typedef uint8_t BITFLG;
typedef uint16_t USHORT;
typedef uint32_t SEQOFFS;
typedef int32_t SEQIDX;
typedef uint32_t SEEDIDX;
typedef uint32_t SEGMIDX;
typedef int32_t SEEDLEN;
typedef int32_t SEGMLEN;
typedef uint32_t SEQPOS;
typedef uint8_t BOOL;

typedef struct _SEED { 
  /**< Defines an exactly matching segment 
   * (covered by a run of overlapping tuples)
   * query offset (as the number of bases) 
   *   qo = SEED.sqo&HASHHIT_HALFMASK
   * reference offset (as the number of k-ktuples)
   *   on forward strand: 
   *  <code> rs_F = (SEED.sqo + qo/NSKIP)&HASHHIT_SOFFSMASK </code>;
   *   on reverse strand: <code> rs_R = SEED.sqo - qo/NSKIP </code>;
   *   rs_R is the number of the last matching k-tuple (k-tuple numbers begin
   *   with the first k-tuple on the forward strand)
   *
   * The matching segment on the reverse strand is calculated by
   * reversing the corresponding segment on the forward strand and
   * subsituting the complement bases. The start point of the matching
   * segment on the forward strand (as the number of bases) is
   *
   *   rso_F = rs_F*NSKIP;
   *   rso_R = rs_R*NSKIP - SEED.len + KTUP;
   *
   * NSKIP is the skip step size.
   * KTUP is the k-tuple length.
   */
  /*   int64_t shift; */
  /*   SEQOFFS q_offs; */
  uint64_t sqo;  
  /**< shift (upper 33 bits) and query offset (lower 31 bits) */
  SEEDLEN len; 
  /**< Length  of the exactly matching segment as the number of bases covered.
   * The number of ktuples is (SEED.len-KTUP)/nskip + 1 */
} SEED;
typedef struct _SEED *SEEDARR;

typedef struct _HITREGION { 
  /**< Points to seeds in the same region of the reference
   * sequences. This is used mainly to avoid multiple fetches of
   * overlapping regions on the reference sequence during alignment */
  SEEDIDX idx; 
  /**< start index, points into hit, seed or segment array */
  SEEDLEN num;
  /**< number of hits or seeds */
} HITREGION;
typedef struct _HITREGION *HITREGARR;

typedef struct _SEGMENT { 
  /**< Segment of seeds with identical shift */
  SEEDIDX ix; 
  /**< index of the first seed in the segment */
  SEGMLEN nseed; 
  /**< number of seeds in the segment. Can be < 0, in which case the
   * segment is flagged out and the number of seeds is -nseed */
  SEGCOV_t cover; 
  /**< number of bases covered by seeds */
} SEGMENT;
typedef struct _SEGMENT *SEGMARR;

struct _SegQMask { 
  /**< Essentially just a string of flags for coverage calculations */
  UCHAR *maskp;    /**< String of flags */
  size_t n_alloc;  /**< Allocated memory for coverage calculations */
  int blksiz;      /**< Block size for memory allocation */
};

struct _SegLst {
  HITREGARR hregr; /**< Array of regions of k-tuple hits */
  SEEDARR seedr;   /**< Array of exactly matching segments (seeds) */
  SEGMARR segmr;   
  /**< Array of segments that are covered by seeds of constant shift */
  USHORT dshift_cutoff; /**< segmenting max difference in shift */
  UCHAR nskip;
  UCHAR ktup;
  BITFLG flags;
  SEQPOS qlen;
  SEGCOV_t maxcover; 
  /**< maximum coverage for segments of constant shift in segmr */
};

typedef struct _SEGCAND { /**< candidate for dynamic programming step */
  uint32_t qs;         /**< query segment start */
  uint32_t qe;         /**< query segment end, points to last base of last k-tuple */
  uint32_t rs;         /**< reference segment start as kmer word number */
  uint32_t re;         /**< reference segment end as k-tuple serial number 
			    * points to first base of last k-mer word */
  short shiftoffs;         /**< minimum shift for start of pair s = (rs - qs/nskip) for forward
			    * or reverse strand relative to the situation where both segments
			    * are aligned at the beginning of the reference segment */
  short shift2mm;          /**< shift for Smith-Waterman free alignment with max. 2 mismatches.
			    * Value set if flag&SEGCANDFLG_MMALI is set */
  short srange;            /**< shift range as multiples of the skip step size 
			    * this is == 0 for candiates without indels */
  SEGCOV_t cover;          /**< the number of bases in the query segment
			    * that are covered by seeds, cover */
  SEGBITFLG_t flag;        /**< Combination of SEGCAND_FLAG bits */
  uint32_t segix;      /**< Index of the segment of constant shift (in SegLst) 
			    * this candidate segment corresponds to */
  int32_t nseg;        /**< Number of constant-shift segments (in SegList) that
			    * are part of this segment. 0 indicates that this is a candidate
			    * for mismatch alignment without indels (i.e. Smith-Waterman free) */
  uint32_t hregix;     /**< Index of the hit region this segment comes from */
  int32_t seqidx;      /**< Serial number of the reference sequence.  SEGCAND_UNKNOWN_SEQIDX
			 * if undetermined */
} SEGCAND;

typedef struct _SEGCAND *CANDARR;

struct _SegAliCands { /** Segments as candidates for dynamic programming */
  CANDARR candr;        /**< Array of candidate segments, in order along the concatenated
			 * set of reference sequences (forward and reverse are consecutive) */
  SEGCOV_t *sort_keys;  /**< Keys for sorting */
  SEGNUM_t *sort_idx;   /**< indices into array candr, sorted with the key */
  SEGNUM_t n_alloc;     /**< memory allocated for sort arrays */
  SEGNUM_t alloc_blksz; /**< Block size for memory allocation */
  SEGNUM_t n_sort;      /**< Number of elements in sort arrays */
  SEGCOV_t max_cover;   /**< current maximum coverage */
  SEGCOV_t max2nd_cover;/**< 2nd highest coverage */
  SEGCOV_t cover_deficit[NUMBER_OF_STRANDS]; 
  /**< cover_deficit is the maximum number of bases that might not be
   * covered by k-tuples because those k-tuples were not sampled (for
   * forward [0] and reverse [1] strands). */
  SEGNUM_t n_mincover;  /**< number of segments above the minimum cover threshold */
  UCHAR nskip;          /**< sampling step for k-mers */
  UCHAR ktup;           /**< length k of k-mers */
};

/******************************************************************************
 *********************************** Macros ***********************************
 ******************************************************************************/

#define MIN(a,b) ((a)<(b))? (a):(b)
#define MAX(a,b) ((a)>(b))? (a):(b)

#define INIT_COVERAGE_CALC(segmp, seedr, maskp, qlen)\
{\
  SEGMLEN l;\
  SEEDLEN q;\
  UCHAR *ucp;\
  const SEED *sep = (seedr) + (segmp)->ix;\
  memset((maskp), 0, (qlen));\
  for (l=(segmp)->nseed; l>0; l--, sep++) {\
    ucp = (maskp) + (sep->sqo&HASHHIT_HALFMASK);\
    for(q=0; q<sep->len; q++) ucp[q] = 1;\
  }\
}\
 
#ifdef RESULTS_TRACKER
#define TRACK_SEED(trkpos_hi, trkpos_lo, trk_flg, cover_flag, ktup, nskip, seedp, nseed) \
{\
  if ((trkp) != NULL && !(cover_flag)) {\
    SEQOFFS rs;\
    SEGMLEN l;\
    for(l=0; l<(nseed); l++) {\
      SEQPOS qo = (seedp)[l].sqo&HASHHIT_HALFMASK;\
      SEEDLEN trk_len = ((seedp)[l].len - (ktup))/(nskip);\
      if ((trk_flg)&TRACKFLG_REVERSE) {\
        rs = (((seedp)[l].sqo>>HASHHIT_HALFBIT) - qo/(nskip))&HASHHIT_SOFFSMASK;\
        if (rs <= (trkpos_hi) && rs - trk_len >= (trkpos_lo)) (cover_flag) = 1;\
      } else {\
        rs = (((seedp)[l].sqo>>HASHHIT_HALFBIT) + qo/(nskip))&HASHHIT_SOFFSMASK;\
        if (rs >= (trkpos_lo) && rs + trk_len <= (trkpos_hi)) (cover_flag) = 1; \
      }\
    }\
  }\
}
#endif

#define CALC_COVERAGE(cover_new, segmp, seedr, maskp)\
{\
  SEGMLEN l;\
  SEEDLEN q;\
  UCHAR *ucp;\
  const SEED *sep = (seedr) + (segmp)->ix;\
  (cover_new) = 0;\
  for (l=(segmp)->nseed; l>0; l--, sep++) {\
    ucp = (maskp) + (sep->sqo&HASHHIT_HALFMASK);\
    for(q=0; q<sep->len; q++) if (!(ucp[q])) {(cover_new)++;ucp[q]=1;}\
  }\
}\

/******************************************************************************
 ******************************* Pivate Methods *******************************
 ******************************************************************************/

#ifdef segment_debug
static void fprintfHitRegion(FILE *fp, const HITREGION * const hrgp, 
			     const SEGMARR segmr) 
{
  SEEDLEN i;
  fprintf(fp, "++++ hit region ++++\n");
  fprintf(fp, "number of segments: %u\n", hrgp->num);
  for (i=0; i<hrgp->num; i++) {
    const SEGMENT *sgp = segmr + hrgp->idx + i;
    fprintf(fp, "segment [%u] seeds %u - %u (%i), cover = %u\n", 
	    i, sgp->ix, sgp->ix + sgp->nseed - 1,
	    sgp->nseed, sgp->cover);
  }
}
#endif
 
static SEGCOV_t calcIndelFreeMincover(SEQOFFS slen, UCHAR ktup, UCHAR nskip)
     /**< Return the minimum coverage required for a mapping to be indel free 
      */
{
  USHORT diffcover_1mm = ktup + nskip - 1; /* max coverage knocked out by 1 mismatch */
  USHORT diffcover_2mm = diffcover_1mm<<1;       /* max coverage knocked out by 2 mismatches */
  SEGCOV_t mincover;

  if (slen >= diffcover_2mm+ktup) {
    mincover = slen - diffcover_2mm;
  } else if (slen >= diffcover_1mm+ktup) {
    mincover = slen - diffcover_1mm;
  } else {
    mincover = ktup;
  }
  return mincover;
}

static void printSeed(FILE *fp, const SEED *sp, char is_reverse, UCHAR nskip, UCHAR ktup)
{
  uint32_t qo = sp->sqo&HASHHIT_HALFMASK;
  uint64_t so, shift = sp->sqo>>HASHHIT_HALFBIT;

  if (is_reverse) {
    so = shift - qo/nskip; /* start in reference segment */
    so = so*nskip + ktup - sp->len;
  } else {
    so = (shift + qo/nskip)&HASHHIT_SOFFSMASK;
    so *= nskip;
  }
  /* reference offset is given as k-mer word number */
  fprintf(fp, "SEED shift = %llu, s = %llu, q = %u, l = %i, strand = %i\n", 
	  (long long unsigned int) shift, (long long unsigned int) so, qo, sp->len, (int) is_reverse);

}

static int defineHitRegions(HITREGARR *idxr, int *n_added, USHORT *max_dshift,
			    uint32_t min_ktup, const HashHitList *hhlp)
     /**< Segment the list of hits which were sorted by shift in
      * ascending order. Segmentation is achieved using cut-off in the
      * maximum shift difference between sucessive hits. Segments are
      * discarded that have fewer seeds than min_ktup.
      */
{
  UCHAR ktup, nskip;
  int i, j, nhits;
  uint32_t ds;
  uint64_t dsthresh;
  SEQPOS qlen;
  size_t n0;
  HITREGION *sp;
  
  const uint64_t *shdat = hashGetHitListData(&nhits, NULL,
					   &qlen, &ktup, &nskip, NULL, hhlp);
  if ((n_added)) {
    n0 = ARRLEN(*idxr);
    *n_added = 0;
  } else {
    ARRLEN(*idxr) = 0;
    n0 = 0;
  }

  if (nhits < 1)
    return ERRCODE_SUCCESS;
     
  /* calculate the segmenting shift difference */
  *max_dshift = ktup*SEGMENTING_DIFFSHIFT/nskip;
  ds = (qlen - ktup)/nskip + 1;
  if (ds < *max_dshift) *max_dshift = (USHORT) ds;
  dsthresh = ((uint64_t) *max_dshift)<<HASHHIT_HALFBIT;

  for (i=0; i<nhits;) {
    for (j=i+1; j<nhits; j++) {
      if ((shdat[j]-shdat[j-1]) >= dsthresh)
	break;
    }
    if (((uint32_t) (j-i)) >= min_ktup){
      ARRNEXTP(sp, (*idxr));
      if (!sp)
	return ERRCODE_NOMEM;
      sp->idx = i;
      sp->num = j-i;
    }
    i = j;
  }
  if ((n_added)) {
    n0 = ARRLEN(*idxr) - n0;
    if (n0 > INT_MAX)
      return ERRCODE_OVERFLOW;
    *n_added = (int) n0;
  }

  return ERRCODE_SUCCESS;
}

static int makeSeedsFromHits(SEEDARR *seedr, 
#ifdef RESULTS_TRACKER
			     Track *trkp,
#endif
			     const HITREGARR idxr,
			     int reg_start, int nreg,
			     const HashHitList *hhlp)
     /**< Join overlapping k-mer word hits resulting in seeds,
      * begin with hit region idxr[reg_start] and end with
      * idx[reg_start + nreg - 1] */
{
  UCHAR ktup, nskip;
  int s;
  SEQPOS lastq, qo, qoffs;
/* #ifdef segment_debug */
/*   SEQPOS so; */
/* #endif */
  SEEDIDX i, j, end_idx;
  SEED *sdp;
  uint64_t shift;
  const uint64_t *shdat = hashGetHitListData(NULL, NULL,
					     NULL, &ktup, &nskip, NULL, hhlp);
#ifdef RESULTS_TRACKER
  BOOL trk_is_covering = 0;
  SEQPOS trkpos_lo, trkpos_hi;
  TRACKFLG_t trk_flg;
  trackGetKtupData(trkp, &trkpos_lo, &trkpos_hi, &trk_flg);
#endif

  if (!(nreg)) {
    ARRLEN(*seedr) = 0;
    reg_start = 0;
    nreg = ARRLEN(idxr);
  }
  for (s=reg_start; s<nreg; s++) {
    i=idxr[s].idx;
    end_idx = i + idxr[s].num;
    idxr[s].idx = ARRLEN(*seedr);
    for (; i<end_idx;) {
      ARRNEXTP(sdp, *seedr);
      if (!sdp)
	return ERRCODE_NOMEM;
      sdp->sqo = shdat[i];
      shift = sdp->sqo&~((uint64_t) HASHHIT_HALFMASK);
      qoffs = sdp->sqo&HASHHIT_HALFMASK;
      lastq = qoffs + ktup;
      for (j=i+1; j<end_idx; j++) {
	if ((shdat[j]&~((uint64_t) HASHHIT_HALFMASK)) != shift)
	  break;
	qo = shdat[j]&HASHHIT_HALFMASK;
	if (qo > lastq || 
	    ((qo - qoffs)%nskip) /* out of register */
	    )
	  break;
	lastq = qo + ktup;
      }
      sdp->len = (lastq - qoffs);
/* #ifdef segment_debug */
/*       qo = sdp->sqo&HASHHIT_HALFMASK; */
/*       so = ((sdp->sqo>>HASHHIT_HALFBIT) + qo/nskip)&HASHHIT_SOFFSMASK; */
/*       if (so > 418541720 && so < 418541723) */
/* 	printf("segment_debug: region %i, selected [%i] so = %u, qo = %u, l = %i\n", */
/* 	       s, (int) ARRLEN(*seedr), so, qo, sdp->len); */
/* #endif */
#ifdef RESULTS_TRACKER
      TRACK_SEED(trkpos_hi, trkpos_lo, trk_flg, trk_is_covering, ktup, nskip, sdp, 1);
#endif
      i = j;
    }
    idxr[s].num = ARRLEN(*seedr) - idxr[s].idx;
  }
#ifdef RESULTS_TRACKER
  if (trkp != NULL && (trk_is_covering)) {
    trackSetFlag(trkp, TRACKFLG_SEEDS_COVERED);
  }
#endif

  return ERRCODE_SUCCESS;
}

static int makeSegmentsFromSeeds(SEGMARR *segmr, SEGCOV_t *maxcover,
				 const HITREGARR hregr,
				 SEEDIDX reg_start, SEEDIDX nreg,
				 const SEEDARR seedr, UCHAR nskip)
{
  uint32_t r;
  SEEDIDX i, j, end_i, d;
  SEGMENT *sgp;
  SEQOFFS qoffs;
  uint64_t shift;

  if (!(nreg)) {
    ARRLEN(*segmr) = 0;
    reg_start = 0;
    nreg = ARRLEN(hregr);
  }
  *maxcover = 0;

  for (r=reg_start; r<nreg; r++) {
    i = hregr[r].idx;
    end_i = i + hregr[r].num;
    hregr[r].idx = ARRLEN(*segmr);
    hregr[r].num = 0;
    for (;i<end_i;) {
      ARRNEXTP(sgp, *segmr);
      if (!sgp)
	return ERRCODE_NOMEM;
      hregr[r].num++;
      sgp->ix = i;
      sgp->cover = seedr[i].len;
      shift = seedr[i].sqo&~((uint64_t) HASHHIT_HALFMASK);
      qoffs = seedr[i].sqo&HASHHIT_HALFMASK;
      for (j=i+1; j<end_i; j++) {
	if ((seedr[j].sqo&~((uint64_t) HASHHIT_HALFMASK)) != shift ||
	    ((seedr[j].sqo&HASHHIT_HALFMASK)-qoffs)%nskip)
	  break;
	sgp->cover += seedr[j].len;
      }
      d = j-i;
      if (d > INT_MAX) 
	return ERRCODE_OVERFLOW;
      sgp->nseed = (SEGMLEN) d;
      if (sgp->cover > *maxcover)
	*maxcover = sgp->cover;
      i=j;
    }
  }

  return ERRCODE_SUCCESS;
}


static int calcSegmentOverlap(const SEGMENT *sgAp, const SEGMENT *sgBp, 
			      const SEEDARR seedr)
     /**< Calculate the overlap between two segments of constant shift on the
      * query sequence */
{
  SEQPOS qa_start, qb_start, qa_end, qb_end;
  SEGMLEN a, b;
  SEGCOV_t overlap = 0;
 
  qb_start = seedr[0].sqo&HASHHIT_HALFMASK;
  a=b=0;
  while (a<sgAp->nseed && b<sgBp->nseed) {
    qa_start = seedr[a].sqo&HASHHIT_HALFMASK;
    qa_end = qa_start+seedr[a].len;
    qb_start = seedr[b].sqo&HASHHIT_HALFMASK;
    qb_end = qb_start + seedr[b].len;

    for (a++; qb_start >= qa_end; a++) {
      qa_start = seedr[a].sqo&HASHHIT_HALFMASK;
      qa_end = qa_start+seedr[a].len;
    }
    for (b++; qa_start >= qb_end; b++) {
      qb_start = seedr[b].sqo&HASHHIT_HALFMASK;
      qb_end = qb_start + seedr[b].len;
    }
    overlap += MIN(qa_end, qb_end) - MAX(qa_start, qb_start);
  }
  return overlap;
}

#define CALC_SEGMENT_BOUNDARIES(segcandp,segmentp)\
{\
   const SEED *sp = seedr + (segmentp)->ix;\
   const SEED *ep = sp + (segmentp)->nseed - 1;\
   (segcandp)->qs = sp->sqo&HASHHIT_HALFMASK;\
   (segcandp)->qe = ep->sqo&HASHHIT_HALFMASK + ep->len - 1;\
   if ((is_reverse)) {\
      (segcandp)->rs = (endp->sqo>>HASHHIT_HALFBIT) -\
                       ((endp->sqo&HASHHIT_HALFMASK) - endp->len + ktup)/nskip;\
      (segcandp)->re = ((sp->sqo>>HASHHIT_HALFBIT) +\
                       ((segcandp)->qs/nskip))&HASHHIT_SOFFSMASK;\
   } else { \
      (segcandp)->rs = ((sp->sqo>>HASHHIT_HALFBIT) + (segcandp)->qs/nskip)&HASHHIT_SOFFSMASK; \
      (segcandp)->re = ((endp->sqo>>HASHHIT_HALFBIT) + \
                        ((endp->sqo&HASHHIT_HALFMASK) + endp->len - ktup)/nskip)&HASHHIT_SOFFSMASK; \
   }\
}

static void calcSegmentBoundaries(SEQOFFS *qs, SEQOFFS *qe, SEQOFFS *rs, SEQOFFS *re,
				  const SEGMENT *sgp, const SEEDARR seedr, 
				  UCHAR ktup, UCHAR nskip, 
				  BOOL is_reverse)
     /**< Derrive the outer boundaries of an identical-shift-segment
      * 
      * \param qs Returns the offset of the first base of the segment in the query sequence (counting from 0)
      * \param qe Returns the offset of the last base of the segment in the query sequence (counting from 0)
      * \param rs Returns the offset of the first base of the segment in the set of concatenated 
      *        reference sequences (counting from 0)
      * \param re Returns the offset of the last base of the segment in the set of concatenated 
      *        reference sequences (counting from 0)
      * \param sgp Segment consisting of seeds of identical shift.
      * \param seedr Array of seeds
      * \param ktup K-mer word length.
      * \param nskip Skip step size.
      * \param is_reverse 0: on forward strand, 1: on reverse strand
      */
{
  const SEED *seed_startp = seedr + sgp->ix;
  const SEED *seed_endp = seed_startp + sgp->nseed - 1;
  *qs = seed_startp->sqo&HASHHIT_HALFMASK;
  *qe = (seed_endp->sqo&HASHHIT_HALFMASK) + seed_endp->len - 1;
  if ((is_reverse)) {
    *rs = ((seed_endp->sqo>>HASHHIT_HALFBIT) - (seed_endp->sqo&HASHHIT_HALFMASK)/nskip)&HASHHIT_SOFFSMASK;
    *rs -= (seed_endp->len - ktup)/nskip;
    *re = ((seed_startp->sqo>>HASHHIT_HALFBIT) - (*qs)/nskip)&HASHHIT_SOFFSMASK;
  } else {
    *rs = ((seed_startp->sqo>>HASHHIT_HALFBIT) + (*qs)/nskip)&HASHHIT_SOFFSMASK;
    *re = ((seed_endp->sqo>>HASHHIT_HALFBIT) + (seed_endp->sqo&HASHHIT_HALFMASK)/nskip)&HASHHIT_SOFFSMASK;
    *re += (seed_endp->len - ktup)/nskip;
  }
  return;
}

/******************************************************************************
 ********************** Private Methods of Type SegQMask **********************
 ******************************************************************************/

static int reallocQMask(SegQMask *p, SEQPOS newsiz)
{
  size_t nsz = ((newsiz-1)/p->blksiz+1)*p->blksiz;
  UCHAR *hp;
  hp = EREALLOCP(p->maskp, nsz);
  if (!hp) 
    return ERRCODE_NOMEM;
  p->maskp = hp;
  p->n_alloc = nsz;
  return ERRCODE_SUCCESS;
}
/******************************************************************************
 *********************** Public Methods of Type SegQMask **********************
 ******************************************************************************/

SegQMask *segQMaskCreate(int qblksz)
{
  SegQMask *sqmp;

  EMALLOCP0(sqmp);
  if ((sqmp)) {
    if (qblksz < DEFAULT_QBLKSZ)
      qblksz = DEFAULT_QBLKSZ;
    ECALLOCP(qblksz, sqmp->maskp);
    if (sqmp->maskp) {
      sqmp->n_alloc = (size_t) qblksz;
      sqmp->blksiz = qblksz;
    } else {
      free(sqmp);
      sqmp = 0;
    }
  }

  return sqmp;
}

void segQMaskDelete(SegQMask *sqmp)
{
  if (sqmp)
    free(sqmp->maskp);
  free(sqmp);
}

/******************************************************************************
 ************************ Public Methods of Type SegLst ***********************
 ******************************************************************************/

SegLst *segLstCreate(int blocksiz)
{
  SegLst *sp;

  EMALLOCP0(sp);
  if (!sp) return 0;

  if (blocksiz < DEFAULT_BLKSZ) 
    blocksiz = DEFAULT_BLKSZ;
  ARRCREATE(sp->seedr, blocksiz);
  ARRCREATE(sp->segmr, blocksiz);
  ARRCREATE(sp->hregr, blocksiz/2);
  if (!((sp->hregr) && (sp->seedr) && 
      (sp->segmr))) {
    segLstDelete(sp);
    sp = 0;
  }

  return sp;
}

void segLstDelete(SegLst *sp)
{
  if (sp) {
    ARRDELETE(sp->hregr);
    ARRDELETE(sp->segmr);
    ARRDELETE(sp->seedr);
  }
  free(sp);
}

void segLstBlank(SegLst *sp)
{
  ARRLEN(sp->seedr) = 0;
  ARRLEN(sp->segmr) = 0;
  ARRLEN(sp->hregr) = 0;
  sp->dshift_cutoff = 0;
  sp->nskip = sp->ktup = 0;
  sp->flags = 0;
  sp->maxcover = 0;
}

int segLstFillHits(SegLst *sglp, 
#ifdef RESULTS_TRACKER
		   Track *trkp,
#endif
		   uint32_t min_ktup,
		   const HashHitList *hhlp)
{
  int errcode;
  char is_reverse;
  const char *qmask;

  segLstBlank(sglp);
  hashGetHitListData(NULL, &is_reverse, &sglp->qlen,
		     &sglp->ktup,
		     &sglp->nskip,
		     &qmask, hhlp);
  sglp->flags = (is_reverse)? SEGLSTFLG_REVERSE: 0;

  /* reduce the minimum number of k-tuples by the number of missing k-tuples */
  for (;(*qmask);qmask++) {
    if (*qmask == HITQUAL_NORMHIT) 
      continue;
    if (min_ktup < 2) 
      break;
    min_ktup--;
  }

  if ((errcode = defineHitRegions(&(sglp->hregr), NULL, &sglp->dshift_cutoff, min_ktup, hhlp)))
    return errcode;
  sglp->flags |= SEGLSTFLG_HITREGIONS;

  if ((errcode = makeSeedsFromHits(&(sglp->seedr), 
#ifdef RESULTS_TRACKER
				   trkp,
#endif
				   sglp->hregr, 0, 0, hhlp)))
    return errcode;
  sglp->flags |= SEGLSTFLG_SEEDS;

#ifdef segment_debug
  segLstPrintSeeds(stdout, sglp);
#endif

  errcode = makeSegmentsFromSeeds(&(sglp->segmr), &sglp->maxcover, 
				  sglp->hregr, 0, 0,  sglp->seedr, sglp->nskip);
  if (!errcode) sglp->flags |= SEGLSTFLG_SEGMENTS;
  return errcode;
} 

int segLstAddHits(SegLst *sglp, 
#ifdef RESULTS_TRACKER
		  Track *trkp,
#endif
		  uint32_t min_ktup,
		  const HashHitList *hhlp)
{
  int errcode, nreg;
  char is_reverse;
  const char *qmask;
  size_t reg_start;

  segLstBlank(sglp);
  hashGetHitListData(NULL, &is_reverse, &sglp->qlen,
		     &sglp->ktup,
		     &sglp->nskip,
		     &qmask, hhlp);
  sglp->flags = (is_reverse)? SEGLSTFLG_REVERSE: 0;

  /* reduce the minimum number of k-tuples by the number of missing k-tuples */
  for (;(*qmask);qmask++) {
    if (*qmask == HITQUAL_NORMHIT) 
      continue;
    if (min_ktup < 2) 
      break;
    min_ktup--;
  }
  reg_start = ARRLEN(sglp->hregr);
  if (reg_start > INT_MAX)
    return ERRCODE_OVERFLOW;
  
  if ((errcode = defineHitRegions(&(sglp->hregr), &nreg, &sglp->dshift_cutoff, min_ktup, hhlp)))
    return errcode;
  sglp->flags |= SEGLSTFLG_HITREGIONS;

  if ((errcode = makeSeedsFromHits(&(sglp->seedr),
#ifdef RESULTS_TRACKER
				   trkp,
#endif 
				   sglp->hregr, (int) reg_start, nreg, hhlp)))
    return errcode;
  sglp->flags |= SEGLSTFLG_SEEDS;

#ifdef segment_debug
  segLstPrintSeeds(stdout, sglp);
#endif

  errcode = makeSegmentsFromSeeds(&(sglp->segmr), &sglp->maxcover, 
				  sglp->hregr, (int) reg_start, nreg,  
				  sglp->seedr, sglp->nskip);
  if (!errcode) sglp->flags |= SEGLSTFLG_SEGMENTS;
  return errcode;
} 

BOOL segLstGetStats(const SegLst *sp, SEEDIDX *nhreg, SEEDIDX *nseed, SEEDIDX *nseg)
{
  if (nhreg) *nhreg = ARRLEN(sp->hregr);
  if (nseed) *nseed = ARRLEN(sp->seedr);
  if (nseg) *nseg = ARRLEN(sp->segmr);

  return sp->flags&SEGLSTFLG_REVERSE;
}

SEEDLEN segLstFetchSeed(SEQOFFS *q_offs, SEQOFFS *s_offs, SEEDIDX idx, const SegLst *sp)
{
  SEED *seedp;
  if (idx >= ARRLEN(sp->seedr)) {
    if (q_offs) *q_offs = 0;
    if (s_offs) *s_offs = 0;
    return 0;
  }
  seedp = sp->seedr + idx;
  if (q_offs) *q_offs = seedp->sqo&HASHHIT_HALFMASK;
  if (s_offs) {
    if (sp->flags&SEGLSTFLG_REVERSE ) {
      *s_offs = (seedp->sqo>>HASHHIT_HALFBIT) - (seedp->sqo&HASHHIT_HALFMASK)/sp->nskip;
    } else {
      *s_offs = ((seedp->sqo>>HASHHIT_HALFBIT) + (seedp->sqo&HASHHIT_HALFMASK)/sp->nskip)&HASHHIT_SOFFSMASK;
    }
  }
  return seedp->len;
}


void segLstPrintSeeds(FILE *fp, const SegLst *sglp)
{
  uint32_t i, n_seeds;
  SEED *sp;

  n_seeds = ARRLEN(sglp->seedr);
  fprintf(fp, "=-=-=-=-= Seed List =-=-=-=-=\n");
  fprintf(fp, "%u Seeds\n", n_seeds);
  
  for (i=0; i<n_seeds; i++) {
    sp = sglp->seedr + i;
    fprintf(fp, "[%i] ", i);
    printSeed(fp, sp, sglp->flags&SEGLSTFLG_REVERSE, sglp->nskip, sglp->ktup);
  }
  
  fprintf(fp, "=-=-= End of Seed List =-=-=\n");
}

/******************************************************************************
 *********************** Methods of Private Type SEGCAND **********************
 ******************************************************************************/

static int printSEGCAND(FILE *fp, const SEGCAND *scp)
{
  int nchar = fprintf(fp, "SEGCAND qs = %u, qe = %u, rs = %u, re = %u, shiftoffs = %hi, "\
	  "seqidx = %i, "\
	  "srange = %hi, cover = %u, strand = %c, shift2mm = %hi, mmali = %hi\n",
	  scp->qs, scp->qe, scp->rs, scp->re, scp->shiftoffs, scp->seqidx, scp->srange, 
	  scp->cover, (scp->flag&SEGCANDFLG_REVERSE)? 'R':'F', 
	  scp->shift2mm, (short) (scp->flag & SEGCANDFLG_MMALI));
  return (nchar > 0)? ERRCODE_SUCCESS: ERRCODE_FAILURE;
}

static int derriveSEGCAND(SEGCAND *candp, SEEDLEN segix_start, SEEDLEN nseg,
			  SEGMENT *segmbasp, const SEEDARR seedr, UCHAR ktup, UCHAR nskip, 
			  SEGCOV_t cover, SEGCOV_t mincover_noindel,
			  uint32_t hregix, BOOL is_reverse)
     /**< Derrive a pair of segments as candidates for pairwise alignment.
      *
      * \param candp Returns the candidate segment.
      * \param segix_start Index of the first identical-shift-segment in segmbasp
      * \param nseg Number of identical-shift-segments  for this candidate segment
      * \param segmbasp Array of identical-shift-segments
      * \param seedr Array of seeds.
      * \param ktup Length of k-mer words.
      * \param nskip Skip step with which k-mer words were sampled.
      * \param cover Coverage of the segment.
      * \param mincover_noindel Minimum coverage required for the candindate segment to be flagged
      *        as a candidate for direct alignment.
      * \param hregix index of the hit region this segment belongs to.
      * \param is_reverse Segment is on forward (0) or reverse (1) strand.
      */
{
  SEEDLEN n;
  SEQPOS qs, qe, rs, re;
  uint64_t shift_range;
  UCHAR flag = 0;
  int64_t shift_min, shift_start, diff_shift, shift_2mm;
  SEGCOV_t maxcover;
  SEGMENT *segmp, *segmentp = segmbasp + segix_start;
  uint64_t offbit = ((uint64_t) 1)<<(HASHHIT_HALFBIT+1);
  /* don't declare this as const, too large for 32-bit systems */

  if (segmentp->nseed <0) 
    return ERRCODE_ASSERT;
  /* return boundaries of the first identical-shift-segment */
  calcSegmentBoundaries(&candp->qs, &candp->qe, &candp->rs, &candp->re, 
			segmentp, seedr, ktup, nskip, is_reverse);
  segmentp->nseed *= -1; /* seed flagged out */
  shift_2mm = shift_min = (int64_t) (seedr[segmentp->ix].sqo>>HASHHIT_HALFBIT);
  /* shift_min is smallest shift 
   * shift_2mm eventually will contain the shift of the constant-shift-segment
   * that hast the most bases covered by k-mers fromt the hash-index */
  maxcover = segmentp->cover;
  segmp = segmentp+1;

#ifdef segment_debug
  fprintf(stderr, "#derriveSEGCAND: segix_start = %i, nseg = %i\n",
	  segix_start, nseg);
  fprintf(stderr, "#derriveSEGCAND: [0] q = (%u, %u) r = (%u, %u) RC = %i, "\
	  "shift_min = %lli\n",
	  candp->qs, candp->qe, candp->rs, candp->re, 
	  (int) (is_reverse != 0),
	  (signed long long) shift_min);  
#endif

  /* update boundaries */
  for (n=1; n<nseg; n++, segmp++) {
    if (segmp->nseed <0) 
      return ERRCODE_ASSERT;
    calcSegmentBoundaries(&qs, &qe, &rs, &re, segmp, seedr, ktup, nskip, is_reverse);
    if (segmp->cover > maxcover) {
      /* possible candidate for direct alignment w/o indels -> remember shift */
      shift_2mm = seedr[segmp->ix].sqo>>HASHHIT_HALFBIT;
      maxcover = segmp->cover;
    }

    segmp->nseed *= -1;
    if (qs < candp->qs) 
      candp->qs = qs;
    if (qe > candp->qe) 
      candp->qe = qe;
    if (rs < candp->rs) 
      candp->rs = rs;
    if (re > candp->re)
      candp->re = re;
  }
  segmp--;
  
  /* calculate the reference shift at the beginning (relative to read)
   * of the resulting candidate segment */
  if (is_reverse) {
    flag |= SEGCANDFLG_REVERSE;
    /* calculate the reference shift at the beginning of 
     * the reference segment */
    shift_start = ((int64_t) candp->rs) + (candp->qe-ktup+1)/nskip;
    /* substitute the following if shift_start refers to shift at
     * beginning of read: 
     * shift_start = ((int64_t) candp->re) + candp->qs/nskip; */
  } else {
    shift_start = (((int64_t) candp->rs)|offbit) - candp->qs/nskip;
  }
 
  shift_range = ((int64_t)(seedr[segmp->ix].sqo>>HASHHIT_HALFBIT)) - shift_min;
  diff_shift = shift_min - shift_start;

#ifdef segment_debug
  fprintf(stderr, "#derriveSEGCAND: [0] q = (%u, %u) r = (%u, %u) RC = %i, "\
	  "sr = %lli, so = %lli\n",
	  candp->qs, candp->qe, candp->rs, candp->re, 
	  (int) (is_reverse != 0),
	  shift_range, diff_shift);  
#endif
  
  /* diff_shift is the smallest shift between the reference segment and the query segment
   * relative to segments aligned at the reference segment start (i.e. shift_start == 0) */
  if (shift_range > SHRT_MAX)
    return ERRCODE_OVERFLOW;
  if (diff_shift < SHRT_MIN || diff_shift > SHRT_MAX)
    return ERRCODE_OVERFLOW;

  candp->shiftoffs = (short) diff_shift;
  if (maxcover >= mincover_noindel) {
    /* candidate for alignment w/o indels,
     * shift_2mm is the shift of the constant-shift-segment with the most bases
     * covered by k-mers from the hash-index */
    int64_t ds_2mm = shift_2mm - shift_start;
    flag |= SEGCANDFLG_MMALI;
    if (ds_2mm < SHRT_MIN || ds_2mm > SHRT_MAX)
      return ERRCODE_OVERFLOW;
    candp->shift2mm = (short) ds_2mm;
  } else {
    candp->shift2mm = 0;
  }

  candp->flag = flag;
  candp->srange = (short) shift_range;
  candp->cover = cover;
  candp->nseg = nseg;
  candp->hregix = hregix;
  candp->seqidx = SEGCAND_UNKNOWN_SEQIDX;
  
  return ERRCODE_SUCCESS;
}

static int updateCandBoundaries(SEGCAND *sgcp, const SEGMENT *segmp, 
				const SEEDARR seedr,
				UCHAR ktup, UCHAR nskip)
{
  SEQPOS qs, qe, rs, re;
  int64_t shift0, shift, d_shift;

  calcSegmentBoundaries(&qs, &qe, &rs, &re, segmp, seedr, ktup, nskip, sgcp->flag&SEGCANDFLG_REVERSE);
  shift0 =((int64_t) sgcp->rs) - sgcp->qs/nskip + sgcp->shiftoffs;
  shift = ((int64_t) rs) - qs/nskip;
  d_shift = shift - shift0;
  if (d_shift < 0) {
    if (d_shift < SHRT_MIN) 
      return ERRCODE_OVERFLOW;
    sgcp->shiftoffs += d_shift;
    sgcp->srange -= d_shift;
    shift0 = shift;
  } else if (d_shift > sgcp->srange) {
    if (d_shift > SHRT_MAX)
      return ERRCODE_OVERFLOW;
    sgcp->srange = d_shift;
  }
  /* update boundaries */
  if (qs < sgcp->qs) 
    sgcp->qs = qs;
  if (qe > sgcp->qe) 
    sgcp->qe = qe;
  if (rs < sgcp->rs) 
    sgcp->rs = rs;
  if (re > sgcp->re)
    sgcp->re = re;
  shift = sgcp->rs - sgcp->qs/nskip;
  sgcp->shiftoffs += shift0 - shift;
  return ERRCODE_SUCCESS;
}

static int extendCand(SEGCAND *sgcp, SEGMARR segmr, const HITREGARR hregr, const SEEDARR seedr,
		      UCHAR ktup, UCHAR nskip)
     /**< Incorporate SEGMENTs of constant shift to the left and right if
      * they are complementary (on query) or small */
{
  int errcode, ovl;
  SEGMIDX i, endix;
  SEGMENT *sip, *sap = segmr + sgcp->segix;
  HITREGION *hregp = hregr + sgcp->hregix;
  
  endix = hregp->idx + hregp->num;
  if (sap->ix < hregp->idx || sap->ix >= endix)
    return ERRCODE_ASSERT;

  /* first look to the left ... */
  for (i=sgcp->segix-1; i>=hregp->idx; i--) {
    sip = segmr+i;
    if (sip->nseed < 0) 
      break;
    ovl = calcSegmentOverlap(sip, sap, seedr);
    if (ovl > MIN(sip->cover, sap->cover)/2)
      break;
    if ((errcode = updateCandBoundaries(sgcp, sip, seedr, ktup, nskip)))
      return errcode;
    sip->nseed *= -1; /* flag segment out */
  }
  
  /* ... then to the right */
  for (i=sgcp->segix+1; i<endix; i--) {
    sip = segmr+i;
    if (sip->nseed < 0) 
      break;
    ovl = calcSegmentOverlap(sip, sap, seedr);
    if (ovl > MIN(sip->cover, sap->cover)/2)
      break;
    if ((errcode = updateCandBoundaries(sgcp, sip, seedr, ktup, nskip)))
      return errcode;
    sip->nseed *= -1; /* flag segment out */
  }
  return ERRCODE_SUCCESS;
}


static int addCandsFast(CANDARR *candr, SEGCOV_t *maxcover, SEGCOV_t *max2ndcover,
#ifdef RESULTS_TRACKER
			Track *trkp,
#endif
			const SEGMARR segmr, const SEEDARR seedr,
			const HITREGARR hregr, UCHAR *maskp, SEQPOS qlen,
			UCHAR ktup, UCHAR nskip, BOOL is_reverse, SEGCOV_t mincover,
			SEGCOV_t mincover_noindel, SEQIDX seqidx)
{
  int errcode;
  /*   UCHAR *ucp; */
  uint32_t r, nreg;
  SEEDLEN i,j;
  /*   SEGMLEN l; */
  /*   SEEDLEN q; */
  SEGCOV_t cover, cover_new;
  /* const SEED *seedp; */
  SEGMENT *segmentp, *sgp;
  const HITREGION *hitregp;
  SEGCAND *cdp;

#ifdef RESULTS_TRACKER
  BOOL trk_is_covering = 0;
  SEQPOS trkpos_lo, trkpos_hi;
  TRACKFLG_t trk_flg;
  trackGetKtupData(trkp, &trkpos_lo, &trkpos_hi, &trk_flg);
#endif

  nreg = ARRLEN(hregr);
  for (r=0; r<nreg; r++) {
    hitregp = hregr+r;
    segmentp = segmr + hitregp->idx;
#ifdef segment_debug
    fprintfHitRegion(stderr, hitregp, segmr);
#endif
   for(i=0;i<hitregp->num;) {
      sgp = segmentp+i;
      INIT_COVERAGE_CALC(sgp, seedr, maskp, qlen);      
      cover = sgp->cover;
#ifdef RESULTS_TRACKER
      TRACK_SEED(trkpos_hi, trkpos_lo, trk_flg, trk_is_covering, ktup, nskip, seedr+sgp->ix, sgp->nseed);
#endif
      /* look ahead */
      sgp++;
      for(j=i+1; j<hitregp->num; j++,sgp++) {
	if (sgp->nseed < 0) 
	  break; /* up to already flagged out segment */
	/* calculate coverage overlap and complement */
#ifdef RESULTS_TRACKER
	TRACK_SEED(trkpos_hi, trkpos_lo, trk_flg, trk_is_covering, ktup, nskip, seedr+sgp->ix, sgp->nseed);
#endif
	CALC_COVERAGE(cover_new, sgp, seedr, maskp);
	if ((cover_new<<1) < sgp->cover && cover >= mincover)
	  break;
	cover += cover_new;
      }
      if (cover >= mincover) {
	/* segments i, ..., j-1 define a new candidate for alignment */
	ARRNEXTP(cdp, *candr);
	if (!cdp)
	  return ERRCODE_NOMEM;
	if ((errcode = derriveSEGCAND(cdp, i, j-i, segmentp, seedr, ktup, nskip, 
				      cover, mincover_noindel, r, is_reverse)))
	  return errcode;
	cdp->seqidx = seqidx;
	if (cover > (*max2ndcover)) {
	  if (cover > (*maxcover)) {
	    *max2ndcover = *maxcover;
	    *maxcover = cover;
	    /* TUNING: } else if (cover < *maxcover) { */
	  } else if (cover != (*maxcover)) {
	    *max2ndcover = cover;
	  }
	}
      }			      
      i = j;
    }
  }
#ifdef RESULTS_TRACKER
  if ((trk_is_covering))
    trackSetFlag(trkp, TRACKFLG_SEGMENTS_OVERLAP);
#endif
  return ERRCODE_SUCCESS;
}

static int addCandsFromSortedSegments(CANDARR *candr, SEGMARR segmr, 
				      SEGNUM_t nseg, const SEEDIDX *segmidx,
				      const SEEDARR seedr,
				      const HITREGARR hregr, SEEDIDX hridx,
				      SEGCOV_t mincov_noindel,
				      UCHAR ktup, UCHAR nskip, char is_reverse)
     /**< Add candidate segments for alignment from segments sorted in ascending
      * order by coverage.
      * \param candr Array of candidate segments.
      * \param segmr Array of constant-shift segments.
      * \param nseg  Number of segments to be considered.
      * \param segmidx Segment indices into array segmr so that
      *        segmr[segmidx[0]], ...., segmr[segmidx[nseg-1]] are segments
      *        sorted by coverage in ascending order.
      * \param seedr Array of seeds. Elements of segmr point into that.
      * \param hridx Index of the hit region considered. Gets passed onto candidate
      *              segments.
      * \param mincov_noindel Minimum coverage required for read so that no indel
      *               possible.
      * \param ktup Length K of a k-tuple (k-mer word).
      * \param nskip Equidistant sampling step with which k-tuples are recorded along
      *              the reference sequence.
      * \param is_reverse Flag 0 -> segments in array segmr are on the forward strand
      *                        1 -> segments refer to the reverse complement strand
      */
{
  int errcode = ERRCODE_SUCCESS;
  SEEDIDX i, sgx;
  SEGMENT *sgp;
  SEGCAND *cdp;

  for (i=nseg; i>0; i--) {
    sgx = segmidx[i-1];
    sgp = segmr + sgx;
    if (sgp->nseed <= 0) 
      continue;
    ARRNEXTP(cdp, *candr);
    if (!cdp)
      return ERRCODE_NOMEM;
    calcSegmentBoundaries(&cdp->qs, &cdp->qe, &cdp->rs, &cdp->re,
			  sgp, seedr, ktup, nskip, is_reverse);
    if ((is_reverse)) cdp->flag = SEGCANDFLG_REVERSE;
    cdp->shift2mm = cdp->shiftoffs = 0;
    cdp->srange = 0;
    cdp->cover = sgp->cover;
    if (sgp->cover >= mincov_noindel)
      cdp->flag |= SEGCANDFLG_MMALI;
    cdp->segix = sgx;
    cdp->nseg = 1;
    cdp->hregix = hridx;
    cdp->seqidx = SEGCAND_UNKNOWN_SEQIDX;
    sgp->nseed *= -1; /* flag out by setting sign bit */
    /* look to join neighbouring segments to the left (smaller shifts) */
    if ((errcode = extendCand(cdp, segmr, hregr, seedr, ktup, nskip)))
      return errcode;
  }
  return errcode;
}

static int addNoIndelCandsFromSegmentSeeds(CANDARR *candr, SEGCOV_t *max_cover, 
					   SEGMARR segmr, const SEEDARR seedr,
					   const HITREGARR hregr, SEGCOV_t mincover,
					   UCHAR ktup, UCHAR nskip, char is_reverse)
/**< Make candidate segments which consist of seeds of constant shift and have
      * a coverage greater than the specified minimum coverage. The maximum coverage
      * of such segments is returned.
      */
{
  int s, nseg = ARRLEN(hregr);
  SEGMIDX i, end_i;
  SEGCOV_t cover;
  SEGCAND *cdp;

  for (s=0; s<nseg; s++) {
    i = hregr[s].idx;
    end_i = i + hregr[s].num;
    for (;i<end_i;i++) {
      cover = segmr[i].cover;
      if (cover < mincover || segmr[i].nseed < 0) 
	/* coverage not above threshold or segment flagged out */
	continue;
      if (cover > *max_cover) 
	*max_cover = cover;
      ARRNEXTP(cdp, *candr);
      if (!cdp)
	return ERRCODE_NOMEM;
      calcSegmentBoundaries(&cdp->qs, &cdp->qe, &cdp->rs, &cdp->re,
			    segmr + i, seedr, ktup, nskip, is_reverse);
      if ((is_reverse)) cdp->flag = SEGCANDFLG_REVERSE;
      cdp->shiftoffs = 0;
      cdp->srange = 0;
      cdp->cover = cover;
      cdp->segix = i;
      cdp->nseg = 1;
      cdp->hregix = s;
      cdp->seqidx = SEGCAND_UNKNOWN_SEQIDX;
#ifdef segment_debug
      printSEGCAND(stdout, cdp);
#endif
      segmr[i].nseed *= -1; /* flag out by setting sign bit */
    }
  }
  return ERRCODE_SUCCESS;
}

static void labelSEGCANDWithinInsertRange (SEGCAND *candAp, SEGNUM_t nsegA, 
					   SEGCAND *candBp, SEGNUM_t nsegB, 
					   SEQLEN_t dmax, SEQLEN_t dmin)
/**< Compare two sets of Candidate segments and label those segments
 * SEGCANDFLG_MATEDIST that are placed on the reference within the
 * insert range dmin ... dmax. No attention is paid to reverse
 * complement etc.
 */
{
  SEGNUM_t na, nb;
  SEGCAND *scAp, *scBp;
  BOOL is_swapped = nsegA > nsegB;
     
  if (nsegA == 0 || nsegB == 0)
    return;

  if (is_swapped) {
    SEGNUM_t tmp = nsegA;
    nsegA = nsegB;
    nsegB = tmp;
    scAp = candBp;
    scBp = candAp;
  } else {
    scAp = candAp;
    scBp = candBp;
  }
  
  for (na=nb=0; na < nsegA && nb < nsegB;) {
    if (scAp[na].seqidx < scBp[nb].seqidx) {
      na++;
    } else if (scBp[na].seqidx > scBp[nb].seqidx) {
      nb++;
    } else if (scAp[na].re + dmax < scBp[nb].rs) {
      na++;
    } else if (scAp[na].re + dmin > scBp[nb].rs) {
      nb++;
    } else {
      scAp[na].flag |= SEGCANDFLG_MATEDIST;
      scBp[nb].flag |= SEGCANDFLG_MATEDIST;
      nb++;
    }
  }
  return;
}

static int checkCANDARRsequential(const CANDARR candr, SEGNUM_t *rvc_offs)
/** check that candidate segments are sorted first by forward then reverse, then
 * by position. Return the start of the reverse complement segments. */
{
  int rv = ERRCODE_SUCCESS;
  SEGNUM_t n, ns = ARRLEN(candr);
  SEGBITFLG_t reverseflg;
  if (rvc_offs != NULL) *rvc_offs = 0;

  if (ns < 2) {
    return rv;
  }

  reverseflg = candr[0].flag & SEGCANDFLG_REVERSE;
  for (n=1; n<ns && ERRCODE_SUCCESS == rv; n++) {
    if ((candr[n].flag & SEGCANDFLG_REVERSE) != reverseflg) {
      if ((reverseflg & SEGCANDFLG_REVERSE))
	rv = ERRCODE_FAILURE;
      else if (rvc_offs != NULL)
	*rvc_offs = n;
      reverseflg = (candr[n].flag & SEGCANDFLG_REVERSE);
    }
    if (candr[n].seqidx < candr[n-1].seqidx)
      rv = ERRCODE_FAILURE;
    if (candr[n].rs < candr[n-1].rs)
      rv = ERRCODE_FAILURE;
  }

  return rv;
}

/* static int labelCANDARRWithinDistance(CANDARR candAr, CANDARR candBr,  */
/* 				      SEQLEN_t dmax, SEQLEN_t dmin, */
/* 				      UCHAR pairlibflg) */
/*   /\**< Compare two sets of Candidate segments and label those segments SEGCANDFLG_MATEDIST  */
/*    * that are within the insert range  dist_max placed on the reference. */
/*    * \param candAr Array of segments for the 1st mate. */
/*    * \param candBr Array of segments for the 2nd mate. */
/*    * \param dmax Maximum distance. */
/*    * \param dmin Minimum distance. */
/*    * \param pairlibflg Bit flag with the lower two bits indicating the orientation of mates: */
/*    *        0 (00): mates in same direction (|-1-> |-2->, <-2-| <-1-|) */
/*    *        1 (01): mates pointing towards each other (|-1-> <-2-|, |-2-> <-1-|) */
/*    *        2 (10): mates pointing away from each other (<-1-| |-2->, <-2-| |-1->) */
/*    *\/ */
/* { */
/*   int errcode; */
/*   SEGNUM_t na_rvc, nb_rvc; */

/*   if ((errcode = checkCANDARRsequential(candAr, &na_rvc))) */
/*     return errcode; */
  
/*   if ((errcode = checkCANDARRsequential(candBr, &nb_rvc))) */
/*     return errcode; */

/*   return errcode; */
/* } */

/******************************************************************************
 ********************* Private Methods of Type SegAliCands ********************
 ******************************************************************************/

static int reallocSortArrays(SegAliCands *sacp, SEEDIDX num)
{
  size_t newsiz = ((num-1)/sacp->alloc_blksz + 1)*sacp->alloc_blksz;
  void *hp;

  if (newsiz > UINT_MAX) 
    return ERRCODE_OVERFLOW;
  hp = EREALLOCP(sacp->sort_keys, newsiz);
  if (!hp) 
    return ERRCODE_NOMEM;
  sacp->sort_keys = hp;
  hp = EREALLOCP(sacp->sort_idx, newsiz);
  if (!hp) 
    return ERRCODE_NOMEM;
  sacp->sort_idx = hp;
  sacp->n_alloc = (SEGNUM_t) newsiz;

  return ERRCODE_SUCCESS;
}

static int transferParamFromSegLst(SegAliCands *sacp, const SegLst *sglp)
{
  int errcode = ERRCODE_SUCCESS;

  if (ARRLEN(sacp->candr) == 0) {
    sacp->ktup = sglp->ktup;
    sacp->nskip = sglp->nskip;
  } else if (sglp->ktup != sacp->ktup || sglp->nskip != sacp->nskip)
    errcode = ERRCODE_ASSERT;

  return errcode;
}

/******************************************************************************
 ********************** Public Methods of Type SegAliCands ********************
 ******************************************************************************/

SegAliCands *segAliCandsCreate(int blocksiz)
{
  SegAliCands *sacp;
  
  EMALLOCP0(sacp);
  if (!sacp) return 0;

  if (blocksiz < 1)
    blocksiz = DEFAULT_BLKSZ;

  ARRCREATE(sacp->candr, blocksiz);
  ECALLOCP(blocksiz, sacp->sort_keys);
  ECALLOCP(blocksiz, sacp->sort_idx);
  
  if ((sacp->candr) && (sacp->sort_keys) && (sacp->sort_idx)) {
    sacp->n_alloc = sacp->alloc_blksz = blocksiz;
    sacp->max_cover = 0;
    sacp->max2nd_cover = 0;
    sacp->n_sort = 0;
    sacp->nskip = 0;
    sacp->ktup = 0;
    sacp->cover_deficit[0] = 0;
    sacp->cover_deficit[1] = 0;
  } else {
    segAliCandsDelete(sacp);
    sacp = 0;
  } 
  
  return sacp;
}

void segAliCandsDelete(SegAliCands *sacp)
{
  if ((sacp)) {
    ARRDELETE(sacp->candr);
    free(sacp->sort_keys);
    free(sacp->sort_idx);
  }
  free(sacp);
}

void segAliCandsBlank(SegAliCands *sacp)
{
  if ((sacp)) {
    ARRLEN(sacp->candr) = 0;
    sacp->n_sort = 0;
    sacp->max_cover = 0;
    sacp->max2nd_cover = 0;
    sacp->n_mincover = 0;
    sacp->nskip = 0;
    sacp->ktup = 0;
    sacp->cover_deficit[0] = 0;
    sacp->cover_deficit[1] = 0;
  }
}

int segAliCandsAddFast(SegAliCands *sacp, SegQMask *qmp,
#ifdef RESULTS_TRACKER
		       Track *trkp,
#endif
		       const SegLst *sglp, SEGCOV_t mincover, SEQIDX seqidx)
{
  int errcode = ERRCODE_SUCCESS;
  //SEGCOV_t mincover_noindel = calcIndelFreeMincover(sglp->qlen, sglp->ktup, sglp->nskip);
  
  if (sglp->qlen >= qmp->n_alloc &&
      (errcode = reallocQMask(qmp, sglp->qlen+1)))
    return errcode;

  if ((errcode = transferParamFromSegLst(sacp, sglp)))
    return errcode;

  errcode = addCandsFast(&sacp->candr, &sacp->max_cover, &sacp->max2nd_cover,
#ifdef RESULTS_TRACKER
			 trkp,
#endif
			 sglp->segmr, sglp->seedr, sglp->hregr,
			 qmp->maskp, sglp->qlen, sglp->ktup, sglp->nskip, 
			 sglp->flags&SEGLSTFLG_REVERSE, mincover, 
			 //mincover_noindel,
			 mincover,
			 seqidx);
  return errcode;
}

int segAliCandsAdd(SegAliCands *sacp, const SegLst *sglp)
{
  int errcode = ERRCODE_SUCCESS;
  uint32_t r, nreg;
  SEEDIDX i, i_end;
  SEGNUM_t j;
  HITREGION *hrp;
  SEGCOV_t mincover_noindel = calcIndelFreeMincover(sglp->qlen, sglp->ktup, sglp->nskip);

  if ((errcode = transferParamFromSegLst(sacp, sglp)))
    return errcode;

  nreg = ARRLEN(sglp->hregr);
  for (r=0; r<nreg; r++) {
    hrp = sglp->hregr+r;
    if (((SEGNUM_t) hrp->num) > sacp->n_alloc &&
	(errcode = reallocSortArrays(sacp, hrp->num)))
      return errcode;
    i_end = hrp->idx + hrp->num;
    for (j=0, i= hrp->idx; i<i_end; i++,j++) {
      sacp->sort_keys[j] = sglp->segmr[i].cover;
      sacp->sort_idx[j] = i;
    }
    sacp->n_sort = j;
    if ((errcode = sort2UINTarraysByQuickSort(sacp->n_sort, 
					      sacp->sort_keys, 
					      sacp->sort_idx)))
      return errcode;
    /* sacp->sort_idx are segment indices relative to sglp->segmr + hrp->idx
     * sorted by coverage in ascending order */
    if ((errcode = addCandsFromSortedSegments(&sacp->candr, sglp->segmr, sacp->n_sort, sacp->sort_idx,
					      sglp->seedr, sglp->hregr, r, mincover_noindel, 
					      sglp->ktup, sglp->nskip,
					      sglp->flags&SEGLSTFLG_REVERSE)))
      return errcode;
  }
  return errcode;
}


int segAliCandsAddNoIndel(SegAliCands *sacp, const SegLst *sglp)
{
  int errcode = ERRCODE_SUCCESS;
  SEGCOV_t mincover = calcIndelFreeMincover(sglp->qlen, sglp->ktup, sglp->nskip);

  if (sglp->maxcover >= mincover) {
    if (!(errcode = transferParamFromSegLst(sacp, sglp))) {
      errcode = addNoIndelCandsFromSegmentSeeds(&sacp->candr, &sacp->max_cover, 
						sglp->segmr, sglp->seedr,
						sglp->hregr, mincover,
						sglp->ktup, sglp->nskip,
						sglp->flags&SEGLSTFLG_REVERSE);
    }
  }
  return errcode;
}

int segAliCandsStats(SegAliCands *sacp,
#ifdef RESULTS_TRACKER
		     Track *trkp,
#endif
		     SEGCOV_t min_cover_below_max, 
		     const HashHitInfo *hhiFp,
		     const HashHitInfo *hhiRp,
		     SEGNUM_t target_depth,
		     SEGNUM_t max_depth,
		     BOOL is_sensitive)
{
  int errcode;
  SEGNUM_t i, j, n_cands = ARRLEN(sacp->candr);
  BOOL is_rev;
  UCHAR nskip = sacp->nskip;
  SEGCOV_t min_cover, cdf = 0;
  SEGCOV_t cover_deficit_adjusted[NUMBER_OF_STRANDS];
  SEGCAND *scp;
#ifdef RESULTS_TRACKER
  TRACKFLG_t trk_flg;
  SEQPOS trk_kixs, trk_kixe;
  BOOL isInSortedSet = 0, isInFullSet = 0;
  SEGNUM_t segRank=0, segIdx = 0;
  trackGetKtupData(trkp, &trk_kixs, &trk_kixe, &trk_flg);
#endif

  if (max_depth< 1 || max_depth > MAXIMUM_DEPTH)
    max_depth =  MAXIMUM_DEPTH;
  if (target_depth < 1) 
    target_depth = DEFAULT_TARGET_DEPTH;
  if (target_depth > max_depth)
    target_depth = max_depth;
  
  /* min_cover_below_max depends on the maximimum number of mismatches
   * one would like to capture.  It is the number of bases of the read
   * that might not be covered by k-tuples because mismatches knock
   * them out.
   */

  min_cover = (min_cover_below_max > sacp->max_cover)? 0: sacp->max_cover - min_cover_below_max;
  /* min_cover: abosolute threshold in k-mer cover of the read */

  if (min_cover > sacp->max2nd_cover) {
    cdf = min_cover - sacp->max2nd_cover;
    min_cover = sacp->max2nd_cover;
  }
  /* Make sure that 2nd-best coverage reads are included, even if cover threshold would
   * exclude them. In this case the absolute threshold min_cover is lowered accordingly
   * and the difference saved in cdf.
   * */

  /* cover_deficit: Array of length 2 with the maximum number of bases
   * that might not be covered by k-tuples because those k-tuples were
   * not sampled (for forward [0] and reverse [1] strands. */

  sacp->cover_deficit[0] = hashCalcHitInfoCoverDeficit(hhiFp);
  sacp->cover_deficit[1] = hashCalcHitInfoCoverDeficit(hhiRp);

  for (i=0; i<NUMBER_OF_STRANDS;i++) {
    //cover_deficit_adjusted[i] = sacp->cover_deficit[0] + nskip;
    cover_deficit_adjusted[i] = sacp->cover_deficit[0];
    if (cover_deficit_adjusted[i] > cdf) {
      cover_deficit_adjusted[i] -= cdf;
    } else {
      cover_deficit_adjusted[i] = 0;
    }
    /* if (cover_deficit_adjusted[i] < nskip) */
    /*   cover_deficit_adjusted[i] = nskip; */
  }

  /* Here, min_cover is a threshold in k-mer cover that can be applied to candidate segments
   * while guaranteeing that the best match and all matches with a certain maximum number of 
   * mismatches below the best match are still found.
   * All candidate segments with a coverage higher than min_cover and
   * within nskip of the maximum coverage will be sampled. If this number is less than target_depth
   * candidate segments with a coverage of nskip within the maximum 
   * coverage will be added up to a maximum number of target_depth.
   */
#ifdef segment_debug
  printf("segment.c::seqAliCandsStats(): unsorted candidate segments\n");
  segAliCandsPrintRaw(stdout, n_cands, sacp);
#endif 

  scp = sacp->candr;
  for (i=j=0; i<n_cands; i++) {
    is_rev = (scp[i].flag & SEGCANDFLG_REVERSE)? 1:0;
#ifdef RESULTS_TRACKER
    isInFullSet = 0;
    if (scp[i].re >= trk_kixs && scp[i].rs <= trk_kixe) {
      trackSetFlag(trkp, TRACKFLG_CANDSEG_FULL);
      isInFullSet = 1;
      segIdx = i;
    }
#endif
    if (scp[i].cover + cover_deficit_adjusted[is_rev] < min_cover) {
      continue;
    }
/*     if (scp[i].cover + cover_deficit_adjusted[is_rev] < sacp->max_cover &&  */
/* 	scp[i].cover + nskip < min_cover) */
/*       continue; */
    if (j >= sacp->n_alloc &&
	(errcode = reallocSortArrays(sacp, j+1)))
      return errcode;
    if (scp[i].cover > sacp->max_cover) 
      return ERRCODE_ASSERT;
    sacp->sort_keys[j] = sacp->max_cover - scp[i].cover;
    sacp->sort_idx[j] = i;
#ifdef RESULTS_TRACKER
    if ((isInFullSet)) {
      isInSortedSet = 1;
      segIdx = i;
    }
#endif
    j++;
  }

#ifdef segment_debug
  printf("segment.c::seqAliCandsStats(): maximum coverage = %u\n", 
	 sacp->max_cover);
  printf("segment.c::seqAliCandsStats(): min_cover_below_max = %u, cover_deficit = [%u, %u]\n", 
	 min_cover_below_max, sacp->cover_deficit[0], sacp->cover_deficit[1]);
  printf("segment.c::seqAliCandsStats(): %i candidate segments above threshold  min_cover = %u\n", 
	 j, min_cover);
#endif    
  
  if ((errcode = sort2UINTarraysByQuickSort(j, sacp->sort_keys, sacp->sort_idx)))
    return errcode;
#ifdef RESULTS_TRACKER
  /* find position in sorted array (if target segment is in the array */
  if ((isInSortedSet)) {
    for (i=0; i<j; i++) {
      if (sacp->sort_idx[i] == segIdx) {
	segRank = i;
	trackSetSegCandRank(trkp, i);
	break;
      }
    }
  }
#endif
  sacp->n_mincover = j;
  
  if (j > target_depth) {
    SEGNUM_t maxj = (j < max_depth)? j: max_depth;
    if (is_sensitive) {
      for (j=target_depth; j<maxj; j++) { 
	is_rev = (scp[j].flag & SEGCANDFLG_REVERSE)? 1:0;
	if (sacp->sort_keys[j] >= cover_deficit_adjusted[(scp[j].flag & SEGCANDFLG_REVERSE)? 1:0])
	  break;
      }
      for (; j<sacp->n_mincover && sacp->sort_keys[j] < nskip; j++);
    } else {
      SEGCOV_t cov = sacp->sort_keys[j/2];
      if (cov < nskip)
	cov = nskip;
      for (j=target_depth; j<maxj && sacp->sort_keys[j] < cov; j++);
    }
  }
  sacp->n_sort = j;

#ifdef RESULTS_TRACKER
  if ((isInSortedSet) && (segRank < sacp->n_sort))
    trackSetFlag(trkp, TRACKFLG_CANDSEG_CUTOFF);
#endif
#ifdef segment_debug
  printf("segment.c::seqAliCandsStats(): %i candidate segments within nskip = %hu of maximum = %u\n", 
	 j, (uint16_t) nskip, sacp->max_cover);
#endif
 
  return ERRCODE_SUCCESS;
}

const SEGCAND *segAliCandsGetSegment(const SegAliCands *sacp, SEEDIDX idx)
{
  return (idx < sacp->n_sort)? sacp->candr + sacp->sort_idx[idx]: 0;
}

void segAliCandsPrint(FILE *fp, SEGNUM_t max_depth, const SegAliCands *sacp)
{
  short i, n_cand;

  n_cand = (sacp->n_sort > max_depth)? max_depth: sacp->n_sort;
  printf("======= List of %hi candidate segments =======\n", n_cand);
  for (i=0; i<n_cand; i++) {
    fprintf(fp, "[%hi] %u: ", i, sacp->sort_idx[i]);
    printSEGCAND(fp, sacp->candr + sacp->sort_idx[i]);
  }
  printf("========= End of candidate segments ==========\n");
}

void segAliCandsPrintRaw(FILE *fp, short max_depth, const SegAliCands *sacp)
{
  int n_cand = ARRLEN(sacp->candr);
  short i, i_end = MIN(max_depth, n_cand);
  printf("======= List of %4i raw candidate segments =======\n", n_cand);
  for (i=0; i<i_end; i++) {
    fprintf(fp, "[%hi]", i);
    printSEGCAND(fp, sacp->candr + i);
  }
  printf("========== End of raw candidate segments ==========\n");
}

SEGNUM_t segAliCandsGetNumberOfSegments(const SegAliCands *sacp, 
					SEGCOV_t *max_cover, 
					SEGCOV_t *max2nd_cover,
					SEGCOV_t *cover_deficitF,  
					SEGCOV_t *cover_deficitR,
					SEGNUM_t *n_mincover)
{
  if (max_cover) *max_cover = sacp->max_cover;
  if (max2nd_cover) *max2nd_cover = sacp->max2nd_cover;
  if (cover_deficitF) *cover_deficitF = sacp->cover_deficit[0];
  if (cover_deficitR) *cover_deficitR = sacp->cover_deficit[1];
  if (n_mincover) *n_mincover = sacp->n_mincover;
  return sacp->n_sort;
}

int segAliCandsPrintSegment(FILE *fp, SEEDIDX scidx, const SegAliCands *sacp)
{
  int errcode = ERRCODE_SUCCESS;
  if (scidx < sacp->n_sort) {
    const SEGCAND *scp = sacp->candr + sacp->sort_idx[scidx];
    printSEGCAND(fp, scp);
  } else {
    errcode = ERRCODE_ASSERT;
  }
  return errcode;
}

int segAliCandsGetSegmentData(SEQPOS *qs, SEQPOS *qe, 
			      SEQPOS *rs, SEQPOS *re,
			      SEEDIDX scidx, const SegAliCands *sacp)
{
  const SEGCAND *scp;
  if (scidx >= ARRLEN(sacp->candr) )
    return ERRCODE_ASSERT;
  
  scp = sacp->candr + scidx;
  if (qs) *qs = scp->qs;
  if (qe) *qe = scp->qe;
  if (rs) *rs = scp->rs;
  if (re) *re = scp->re;
 
  return ERRCODE_SUCCESS;
}
 
int segAliCandsCalcSegmentOffsets(SEQLEN_t *qs, SEQLEN_t *qe, 
				  SETSIZ_t *rs, SETSIZ_t *re,
				  int *band_l, int *band_r,
				  SEQLEN_t *qs_directmatch, int *ro_directmatch,
				  SEQNUM_t *seqidx, SEGBITFLG_t *bitflags,
				  SEGCOV_t *cover,
				  short edgelen, uint32_t qlen, 
				  const SeqSet *ssp,
				  SEEDIDX scidx,
				  const SegAliCands *sacp)

{
  int bl, br, nseq, band_offs;
  int ds;
  int q_edge_l, q_edge_r, r_edge_l, r_edge_r, edge_band;
  UCHAR nskip = sacp->nskip;
  UCHAR ktup = sacp->ktup;
  const SETSIZ_t *soffsp;
  SETSIZ_t roffs, rlen;
  const SEGCAND *scandp;

  if (scidx >= sacp->n_sort)
    return ERRCODE_FAILURE;
  scandp = sacp->candr + sacp->sort_idx[scidx];

  *seqidx = scandp->seqidx;
  *bitflags = scandp->flag;
  *cover = scandp->cover;
  nseq = seqSetGetOffsets(ssp, &soffsp);
  if (scandp->seqidx < 0 || scandp->seqidx >= nseq) {
    roffs = 0;
    rlen = soffsp[nseq];
  } else {
    rlen = seqSetGetSeqDatByIndex(&roffs, NULL, scandp->seqidx, ssp);
  }
    
  *rs = ((uint64_t) scandp->rs)*nskip;
  *re = ((uint64_t) scandp->re)*nskip + ktup - 1;

  if (*rs < roffs || *re < *rs)
    return ERRCODE_ASSERT;

  *rs -= roffs;
  *re -= roffs;

  if (*re >= rlen) 
    return ERRCODE_ASSERT; 

  if (scandp->qe < scandp->qs || scandp->qs >= qlen)
    return ERRCODE_ASSERT;

  if (scandp->flag&SEGCANDFLG_REVERSE) {
    *qs = qlen - scandp->qe - 1;
    *qe = qlen - scandp->qs - 1;
  } else {
    *qs = scandp->qs;
    *qe = scandp->qe;  
  }

  /* calculate band */
  edge_band = (int) (qlen - scandp->cover)/EDGE_BAND_FACTOR;
  if (edge_band > nskip) {
    if (edge_band > (int) (qlen>>MAX_BANDEDGE_2POW)) 
      edge_band = (int) (qlen>>MAX_BANDEDGE_2POW);
    edge_band -= nskip-1;
  } else
    edge_band = 0;
 
  br = (-scandp->shiftoffs + 1)*nskip + edge_band + 1;
  bl = br - (scandp->srange + 2)*nskip - 2*edge_band - 2;

  /* edges extend the segment */
  q_edge_l = (*qs >= ((SEQLEN_t) edgelen) && edgelen > 0)? edgelen: (int) *qs;
  q_edge_r = (*qe + edgelen + 1 <= qlen && edgelen > 0)? edgelen: (int) (qlen - *qe - 1);
  *qs -= q_edge_l;
  *qe += q_edge_r;
 
  r_edge_l = q_edge_l + br;
  r_edge_r = q_edge_r - bl;
 
  if (r_edge_l > 0 && *rs < (SETSIZ_t) r_edge_l) {
    r_edge_l = *rs;
    *rs = 0;   
  } else {
    *rs -= r_edge_l; 
  }
  
  if (*re + r_edge_r >= rlen) {
    r_edge_r = rlen - *re - 1;
    *re = rlen -1;
  } else {
    *re += r_edge_r;
  }
  
  if (*re < *rs)
    return ERRCODE_ASSERT;
  
  /* alignment band is reported along profiled query sequence and
   *  relative to profile origin, not the segment start qs */
  band_offs = q_edge_l - r_edge_l;
  ds = scandp->shift2mm*nskip + band_offs;
  *band_l = bl + band_offs + *qs;
  *band_r = br + band_offs + *qs;
   
  if (ds < 0) {
    *qs_directmatch = *qs - ds;
    *ro_directmatch = 0;
  } else {
    *qs_directmatch = *qs;
    *ro_directmatch = ds;
  }

#ifdef segment_debug
  printf("segment_debug::segALiCandsCalcSegmentOffsets():\n");
  printf("   rs = %llu, re = %llu, qs = %u, qe = %u\n",
	 *rs, *re, *qs, *qe);
  printf("   ro_directmatch = %i, qs_directmatch = %u\n",
	 *ro_directmatch, *qs_directmatch);
  printf("   band_l = %i, band_r = %i\n", *band_l, *band_r);
  
  printf("   r_edge_l = %i, r_edge_r = %i, q_edge_l = %i, q_edge_r = %i\n",
	 r_edge_l, r_edge_r, q_edge_l, q_edge_r);
#endif
  return ERRCODE_SUCCESS;
}
