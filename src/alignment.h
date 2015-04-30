/**< Pairwise sequence alignment using dynamic programming */

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

#ifndef ALIGNMENT_H
#define ALIGNMENT_H

  //#define alignment_debug
  //#define alignment_debug_limits
  //#define alignment_timing 
  //#define alignment_matrix_debug

#if defined alignment_debug || defined algnment_matrix_debug
#include "sequence.h"
#endif

#include "score.h"
#include "diffstr.h"
#include "alibuffer.h"

  /****************************************************************************
   ******************************* Constants **********************************
   ****************************************************************************/

  /****************************************************************************
   ***************************** Opaque types *********************************
   ****************************************************************************/
  
  typedef struct _AliDiffStr AliDiffStr;
  /**< String capturing the pairwise alignment */

  typedef struct _AliRsltSet AliRsltSet;
  /**< Intermediate alignment results */

  /****************************************************************************
   ****************************** General Methods *****************************
   ****************************************************************************/
 
  int aliScoreDiffStr(int *swscor, const char *unpofiled_seq, int unprofiled_seqlen,
		      unsigned int profiled_offs,
		      const DIFFSTR_T *diffstrp, int diffstrlen, const ScoreProfile *scpp);
  /**< Calculate Smith-Waterman score from diff string.
   * \param swscor Returns the Smith-Waterman score.
   * \param unpofiled_seq Segment of the unprofiled sequence (in SEQCOD_MANGLED code).
   * \param unprofiled_seqlen Length of the segment.
   * \param profiled_offs Starting offset of the segment in the profiled sequence.
   * \param diffstrp Compressed alignment string.
   * \param diffstrlen Length of the compressed alignment string (incl. terminatig 0).
   * \param scpp Sequence profile.
   */

  /******************************************************************************
   ************************* Methods of Type AliRsltSet *************************
   ******************************************************************************
   * stores the (partial) alignment resuls for a pair of segments
   */

  AliRsltSet *aliRsltSetCreate(const ScoreMatrix *smp, short blksz, short diffblksz,
			       int track_blksz, int track_thresh);
  /**< Constructor.
   * \param smp Matrix of alignment scores (can be NULL, in which case the Smith-
   *            Waterman alignment score are *not* complexity weighted.
   * \param blksz Block size (granularity) for memory allocation as the number of 
   *        results. If 0 defaults are used.
   * \param diffblksz Block size (gramularity) for memory allocation as the number 
   *        of bytes. Default value is used if 0.
   * \param track_blksz Granularity for allocation of direction indicators as the
   *        number of bytes. Default value is used if 0.
   * \param track_thresh Threshold (as the number of bytes) above which memory 
   *        allocated for direction indicators is de-allocated after use. If 0
   *        default value is used.
   */

  void aliRsltSetDelete(AliRsltSet *p);
  /**< Destructor.
   */

  void aliRsltSetReset(AliRsltSet *p);
  /**< Reset the set of results.
   */

  short aliRsltSetGetSize(const AliRsltSet *arp);
  /**< Return the number of results 
   */

  int aliRsltSetFetchData(const AliRsltSet *arp, short idx, int *score, 
			 int *ps_start, int *ps_end, int *us_start, int *us_end,
			 const DiffStr **dfspp);
  /**< Accessor for the alignment results.
   * \return ERRCODE_FAILURE if index is out of bounds.
   * \param arp Set of results.
   * \param idx Index of the result for which data should be fetched.
   * \param score Smith-Waterman alignment score.
   * \param ps_start Start of the aligned segment in the profiled sequence.
   * \param ps_end End of the aligned segment in the profiled sequence.
   * \param us_start Start of the aligned segment in the unprofiled sequence.
   * \param us_end End of the aligned segment in the unprofiled sequence.
   * \param dfspp Compressed alignment string.
   * \note Any of the pointers can be NULL.
   */

  /******************************************************************************
   ************************** Public Alignment Methods **************************
   ******************************************************************************/

int aliSmiWatInBand(AliRsltSet *rssp,
		    AliBuffer *bufp,
		    const ScoreProfile *profp,
#if defined alignment_debug || defined algnment_matrix_debug
		    const SeqCodec *codecp,
		    const char *profiled_seqp,
#endif
		    const char *unprofiled_seqp,
		    int unprofiled_seqlen,
		    int l_edge, int r_edge,
		    int profiled_left, int profiled_right,
		    int unprofiled_left, int unprofiled_right,
		    int minscore, int minscorlen);
  /**< Full banded Smith-Waterman.
   * \param rssp Returns set of results for this pairwise alignment.
   * \param bufp Buffer for rows of the dynamic programming matrix
   * \param pltyp Alignment penalty scores.
   * \param profp Sequence profile of sequence psqp.
   * \param unprofiled_seqp Unprofiled sequence in SEQCOD_MANGLED code.
   * \param unprofiled_seqlen Length of the unprofiled sequence.
   * \param l_edge Left edge of the alignment band along profiled sequence (origin 0).
   * \param r_edge Right edge of the alignment band along profiled sequence (origin 0).
   * \param profiled_left Start base of the segment in the profiled sequence.
   * \param profiled_right Last base of the segment in the profiled sequence.
   * \param unprofiled_left Start base of the segment in the unprofiled sequence.
   * \param unprofiled_right Last base of the segment in the unprofiled sequence.
   * \param minscore Minimum alignment score for a partial alignment to be added.
   * \param minscorlen Minimum length for aligned query segment. 
   */

  int aliSmiWatInBandFast(int *maxscor,
			  AliBuffer *bufp,
			  const ScoreProfile *profp,
#if defined algnment_matrix_debug
			  const SeqCodec *codecp,
			  const char *profiled_seqp,
#endif
			  const char *unprofiled_seqp,
			  int unprofiled_seqlen,
			  int l_edge, int r_edge,
			  int profiled_left, int profiled_right,
			  int unprofiled_left, int unprofiled_right);
  /**< Fast banded Smith-Waterman establishing only the maximum score.
   *
   * \param maxscor Returns the maximum score.
   * \param bufp Buffer for rows of the dynamic programming matrix
   * \param profp Sequence profile of sequence psqp.
   * \param unprofiled_seqp Unprofiled sequence in SEQCOD_MANGLED code.
   * \param unprofiled_seqlen Length of the unprofiled sequence.
   * \param l_edge Left edge of the alignment band along profiled sequence (origin 0).
   * \param r_edge Right edge of the alignment band along profiled sequence (origin 0).
   * \param profiled_left Start base of the segment in the profiled sequence.
   * \param profiled_right Last base of the segment in the profiled sequence.
   * \param unprofiled_left Start base of the segment in the unprofiled sequence.
   * \param unprofiled_right Last base of the segment in the unprofiled sequence.
   */

#ifdef alignment_debug
  /******************************************************************************
   ******************* Public Alignment Methods for Debugging *******************
   ******************************************************************************/

int aliDebugFullSmiWat(AliRsltSet *rssp, AliBuffer *bufp,
		       const ScoreProfile *profp,
		       const SeqCodec *codecp,
		       const char *psqp,
		       const char *usqp, int us_len,
		       int l_edge, int r_edge,
		       int ps_start, int ps_end,
		       int us_start, int us_end);
  /**< Performs direct full banded Smith_Waterman alignment for debugging purposes.
   * \param rssp Set of alignment results (used as buffer, no results are actually created).
   * \param profp Sequence profile.
   * \param psqp Profiled sequence (in SEQCOD_MANGLED code).
   * \param usqp Sequence of the unprofiled segment (in SEQCOD_MANGLED code).
   * \param uslen Length of the unprofiled segment.
   * \param l_edge Left edge of the alignment band along the profiled sequence.
   * \param r_edge Right edge of the alignment band along the profiled sequence..
   * \param ps_start Start base of the segment in the profiled sequence.
   * \param ps_end Last base of the segment in the profiled sequence.
   * \param us_start Start base of the segment in the unprofiled sequence.
   * \param us_end Last base of the segment in the unprofiled sequence.
   *
   * \note: start/end base positions are specified with respect to an origin of 0.
   */
#endif

#ifdef alignment_timing
int aliSmiWatInBandDirect(AliRsltSet *rssp, AliBuffer *bufp,
			  const ScoreProfile *profp,
			  const char *psqp,
			  const char *usqp, int us_len,
			  int l_edge, int r_edge,
			  int ps_start, int ps_end,
			  int us_start, int us_end,
			  int *maxscor);

  /**< Performs direct full banded Smith_Waterman alignment for timing.
   * \param rssp Set of alignment results (used as buffer, no results are actually created).
   * \param profp Sequence profile.
   * \param psqp Profiled sequence (in SEQCOD_MANGLED code).
   * \param usqp Sequence of the unprofiled segment (in SEQCOD_MANGLED code).
   * \param uslen Length of the unprofiled segment.
   * \param l_edge Left edge of the alignment band along the profiled sequence.
   * \param r_edge Right edge of the alignment band along the profiled sequence.
   * \param ps_start Start base of the segment in the profiled sequence.
   * \param ps_end Last base of the segment in the profiled sequence.
   * \param us_start Start base of the segment in the unprofiled sequence.
   * \param us_end Last base of the segment in the unprofiled sequence.
   * \param maxscor Returns the maximum score of the alignment.
   *
   * \note: start/end base positions are specified with respect to an origin of 0.
   */

#endif /* #ifdef alignment_timing */

#endif /* ifndef ALIGNMNENT_H */
#ifdef __cplusplus
}
#endif
