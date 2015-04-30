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

#ifdef __cplusplus
extern "C"
{
#endif

#ifndef SEGMENT_H
#define SEGMENT_H

  //#define segment_debug

#include <stdio.h>
#include <stdint.h>

#include "sequence.h"
#ifdef RESULTS_TRACKER
#include "tracker.h"
#endif
#include "hashhit.h"

  /****************************************************************************
   ****************************** Constants ***********************************
   ****************************************************************************/
 enum SEGCAND_BITFLAGS {  /**< Status flags for type SegCand */
    SEGCANDFLG_REVERSE   = 0x01,  /**< Candidate segment is on reverse complement strand */
    SEGCANDFLG_DISREGARD = 0x02,  /**< Disregard candidate segment (flagged out) */
    SEGCANDFLG_MMALI     = 0x04,  /**< Candidate for mismatch alignment without indels, i.e.
				   * Smith-Waterman free */
    SEGCANDFLG_MATEDIST  = 0x08,  /**< Candidate segment is located within insert distance of 
				   * a segment from a read mate (used to prioritise alignments
				   * of segments from paired reads). */
  };

  enum SEGCAND_CONST {
    SEGCAND_UNKNOWN_SEQIDX = -1,
  };

  /****************************************************************************
   ************************** Transparent types *******************************
   ****************************************************************************/
  
  typedef uint8_t SEGBITFLG_t; /**< Bit flag, combination of SEGCAND_BITFLAGS */
  typedef uint32_t SEGNUM_t; /**< Holds number of segments */
  typedef uint32_t SEGCOV_t; /**< Coverage (num bases covered by k-mer hits */

  /****************************************************************************
   ***************************** Opaque types *********************************
   ****************************************************************************/

  typedef struct _SegQMask SegQMask;
  typedef struct _SegLst SegLst;
  typedef struct _SegAliCands SegAliCands;

  /****************************************************************************
   ************************* Methods of type SegQMask *************************
   ****************************************************************************/ 
  SegQMask *segQMaskCreate(int qblksz);
  /**< Constructor 
   */
  void segQMaskDelete(SegQMask *sqmp);
  /**< Destructor
   */
  
  /****************************************************************************
   ************************** Methods of type SegLst **************************
   ****************************************************************************/ 
  
  SegLst *segLstCreate(int blocksiz);
  /**< Constructor
   */
  void segLstDelete(SegLst *sp);
  /**< Destructor
   */
  void segLstBlank(SegLst *sp);
  /**< Empty the list of segments
   */

  int segLstFillHits(SegLst *sglp, 
#ifdef RESULTS_TRACKER
		     Track *trkp,
#endif
		     uint32_t min_ktup,
		     const HashHitList *hhlp);
  /**< Fill list of segments with runs of exact matches derrived from hash hits.
   * \param sglp List of segments to be filled.
   * \param min_ktup Minimum number of ktuple hits for required for a segment
   * \param hhlp List of k-tuple hits.
   */

  int segLstAddHits(SegLst *sglp, 
#ifdef RESULTS_TRACKER
		    Track *trkp,
#endif
		    uint32_t min_ktup,
		    const HashHitList *hhlp);
  /**< Like sgLstFillHits() but add the segments to existing ones.
   * \param sglp List of segments to be filled.
   * \param min_ktup Minimum number of ktuple hits for required for a segment
   * \param hhlp List of k-tuple hits.
   */ 

  unsigned char segLstGetStats(const SegLst *sp, uint32_t *nhreg,
			       uint32_t *nseed, uint32_t *nseg);
  /**< Return 1 if list is for reverse complement, 0 else
   * \param sp segment list.
   * \param nhreg returns the number of hit regions (can be NULL)
   * \param nseed returns the number of seeds (can be NULL)
   * \param nseg returns the number of segments covered by seeds of constant shift (can be NULL).
   */

  signed int segLstFetchSeed(uint32_t *q_offs, uint32_t *s_offs, 
			     uint32_t idx, const SegLst *sp);
  /**< return seed length and starting positions for seed with index idx.
   */

  void segLstPrintSeeds(FILE *fp, const SegLst *sglp);
  /**< Print the seeds in the segment list 
   */   

  /*****************************************************************************
   ************************* Methods of Type SegAliCands ***********************
   *****************************************************************************/
  
  SegAliCands *segAliCandsCreate(int blocksiz);
  /**< Constructor 
   */
  void segAliCandsDelete(SegAliCands *sacp);
  /**< Destructor
   */
  void segAliCandsBlank(SegAliCands *sacp);
  /**< Empty type without deallocating memory
   */
  int segAliCandsAddFast(SegAliCands *sacp, SegQMask *qmp,
#ifdef RESULTS_TRACKER
			 Track *trkp,
#endif
			 const SegLst *sglp, SEGCOV_t mincover,
			 signed int seqidx);
  /**< Add candidates for alignment from constant-shift segments in type SegLst.
   * Segments are added in sorted order until (i) minimum coverage is reached
   * and (ii) coverage overlap greater than non-overlapping coverage.
   * \param sacp List of candidate segments for alignment.
   * \param qmp Buffer used during coverage calculations.
   * \param sglp List of segments of constant shift from which the candidate segments
   *        are derrived.
   * \param mincover Minimum coverage of k-tuple hits required for a candidate segment.
   * \param seqidx Index of the reference sequence that set of segments was generated
   *        from. Use SEGCAND_UNKNOWN_SEQIDX if undetermined.
   */

  int segAliCandsAdd(SegAliCands *sacp, 
		     const SegLst *sglp);
  /**< Add candidates for alignment from segments in type SegLst.
   * The segments are first sorted by coverage and then joined starting
   * from the segments with the highest coverage.
   */
  int segAliCandsAddNoIndel(SegAliCands *sacp, const SegLst *sglp);
  /**< Add candidate segments from type SegLst that have sufficiently
   * high coverage by seeds so that they probably are free of indels.
   * These might be aligned without Smith-Waterman
   * allowing for 1 or 2 mismatches.
   */
  int segAliCandsStats(SegAliCands *sacp, 
#ifdef RESULTS_TRACKER
		       Track *trkp,
#endif
		       SEGCOV_t min_cover_below_max,
		       const HashHitInfo *hhiFp,
		       const HashHitInfo *hhiRp,  
		       SEGNUM_t target_depth, SEGNUM_t max_depth,
		       uint8_t is_sensitive);
  /**< Collect a subset of candidate sorted by coverage.
   * \param sacp List of candidate segments.
   * \param min_cover_below_max The number of bases of the read that might not be 
   *             covered by k-tuples because of a defined maxmimum number of mismatches
   *             one would like to capture.
   
   * \param hhiFp Hit info for forward strand.
   * \param hhiRp Hit info for reverse strand.
   * \param target_depth minimum number of candidate segments.
   * \param max_depth maximum number of candidate segments.
   * \param is_sensitive If this flag is != 0, all candidate segments are considered 
   *        with a coverage up to and including nskip below the maximum irrespective of
   *        max_depth.
   * \note If min_cover is a coverage threshold applicable to
   * candidate segments while guaranteeing that the best match and all
   * matches with a certain maximum number of mismatches below the
   * best match are still found. If the maximum coverage ovserved among the 
   * segments is max_cover, a maximum value for the threshold is 
   * min_cover = max_cover - min_cover_below_max - cover_dificit;
   * 
   * All candidate segments with a coverage higher than this threshold and within
   * nskip of the max_cover will be sampled.
   * If this number is less than target_depth, candidate segments with a coverage within 
   * MAX(min_cover_below_max, cover_dificit) within max_cover will be added up to a maximum 
   * number of target_depth.
   */

  void segAliCandsPrint(FILE *fp, SEGNUM_t max_depth, const SegAliCands *sacp);
  /**< Print list of candidate segments on stream up to a depth of max_depth
   */
  void segAliCandsPrintRaw(FILE *fp, short max_depth, const SegAliCands *sacp);
   /**< Print list of candidate segments in the order in which they were added 
    * on stream up to a depth of max_depth.
   */

  SEGNUM_t segAliCandsGetNumberOfSegments(const SegAliCands *sacp,
					  SEGCOV_t *max_cover, 
					  SEGCOV_t *max2nd_cover,
					  SEGCOV_t *cover_deficitF, 
					  SEGCOV_t *cover_deficitR,
					  SEGNUM_t *n_mincover);
  /**< Accessor returning the number of candidate segments.
   * \return Number of candidate segments.
   * \param sacp Structure of type SegAliCands.
   * \param max_cover (can be NULL) Returns maximum coverage (number of bases covered by kmer words).
   * \param max2nd_cover (can be NULL) Returns second highest coverage (number of bases covered by kmer words).
   * \param cover_deficitF Maximum mumber of bases that might not be covered by a k-mer hit on forward strand
   *        (can be NULL).
   * \param cover_deficitR Maximum mumber of bases that might not be covered by a k-mer hit on reverse strand
   *        (can be NULL).
   * \param n_mincover (can be NULL) Returns the number of segments that would have to be considered to guarantee
   *        the true hit is found.
   */

  int segAliCandsGetSegmentData(uint32_t *qs, uint32_t *qe, 
				uint32_t *rs, uint32_t *re,
				uint32_t scidx, const SegAliCands *sacp);
  /**< Accessor for debugging purposes
   */

  int segAliCandsPrintSegment(FILE *fp, uint32_t scidx, const SegAliCands *sacp);
  /**< Print data of candidate segment
   * \param fp Output stream.
   * \param scidx Index of the (sorted) candidate segment.
   * \param sacp Structure holding cadidate segments for alignment.
   */

  int segAliCandsCalcSegmentOffsets(SEQLEN_t *qs, SEQLEN_t *qe,
				    SETSIZ_t *rs, SETSIZ_t *re,
				    int *band_l, int *band_r,
				    SEQLEN_t *qo_dirmatch, int *ro_dirmatch,
				    SEQNUM_t *seqidx, SEGBITFLG_t *bitflags,
				    SEGCOV_t *cover,
				    short edgelen, uint32_t qlen, 
				    const SeqSet *ssp,
				    uint32_t scidx,
				    const SegAliCands *sacp);
  /**< Accessor returning the segment boundaries after extension by edgelen to either side.
   * \return error code or ERRCODE_SUCCESS
   * \param qs Returns start position of the segment in (reverse complemented) profile of 
   *           the query sequence (0 based).
   * \param qe Returns last position of the segment in (reverse complemented) profile of
   *           the query sequence (0 based).
   * \param rs Returns start position of the segment in the unprofiled 
   *           (concatenated reference) sequences (0 based).
   * \param re Returns last position of the segment in the unprofiled (concatenated reference)
   *           sequences (0 based).
   * \param band_l Returns the left alignment band limit for banded Smith-Waterman along
   *           (reverse complemented) profile of the query sequence (0 based to start of profile).
   * \param band_r Returns the right alignment band limit for banded Smith-Waterman along
   *           (reverse complemented) profile of the query sequence (0 based to start of profile).
   * \param qo_dirmatch Returns the offset in the query segment for attempting an alignment
   *        without indels (0 based to start of profile).
   * \param ro_dirmatch Returns the offset in the reference segment for attempting an alignment
   *        without indels.
   * \param seqidx Returns the sequence index or SEGCAND_UNKNOWN_SEQIDX if the index is not yet
   *        determined.
   * \param bitflags Returns a combination of SEGCAND_BITFLAGS bit flags.
   * \param cover Number of bases in segment covered by kmer words.
   * \param edgelen Try to extend pairwise segments by this many bases on either end. If 0, extend
   *        to the entire read.
   * \param qlen Length of the query sequence.
   * \param ssp Set of reference sequences.
   * \param scidx Index of the (sorted) candidate segment.
   * \param sacp Structure holding cadidate segments for alignment.
   */


#endif
#ifdef __cplusplus
}
#endif
