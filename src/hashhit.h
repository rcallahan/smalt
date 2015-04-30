/** Gathering hits in hash table */

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

#ifndef HASHHIT_H
#define HASHHIT_H

  //#define hashhit_mergesort /* use merge sort rather than quick sort */
  //#define hashhit_basqual /* quality-ware lookup of k-mer words */
  //#define hashhit_debug
  //#define hashhit_mscor_calib  /* output data for calibration of mapping score */
  //#define hashhit_minimise_coverdeficit /* distribute seeds evenly among frames */
  //#define hashhit_dump_sortarray

#include <stdint.h>

#include "sequence.h"
#ifdef RESULTS_TRACKER
#include "tracker.h"
#endif
#include "hashidx.h"

  /****************************************************************************
   ****************************** Constants ***********************************
   ****************************************************************************/

  /** Hit qualifiers for the k-tuple positions
   * in the query sequence */
  enum HASH_HIT_QUALIFIERS {
    HITQUAL_TERM    = 0, /**< termination, don't change */
    HITQUAL_NORMHIT = 1, /**< k-tuple registered in hit list */
    HITQUAL_MULTIHIT= 2, /**< more hits than cut-off allows */
    HITQUAL_REPEAT  = 3, /**< filtered out by repeat filter */
    HITQUAL_NOHIT   = 4,  /**< no hit in hash table */
    HITQUAL_NONSTDNT = 5, /**< no hit because contains non-standard nucleotide */
    HITQUAL_PARTHIT = 6, /**< Hits only partially sampled */
  };

  enum HASHHIT_CONST {
    HASHHIT_HALFBIT = 31, /**< split 64 bits into 33 bits for shift and 31 bits for query offset */
    HASHHIT_HALFMASK   =  0x7FFFFFFF, /**< mask out lower 31 bits for query offset */
    HASHHIT_SOFFSMASK  =  0xFFFFFFFF, /**< mask out lower 32 bits for subject offset (after >>HASHHIT_HALFBIT)*/
    HASHHIT_SEMIBIT = 32, /**< split 64 bit into 2x32 bits */
  };

  /****************************************************************************
   ******************************** Types *************************************
   ****************************************************************************/

  /****************************************************************************
   ***************************** Opaque types *********************************
   ****************************************************************************/

  typedef struct _HashHitInfo HashHitInfo;
  typedef struct _HashHitList HashHitList;
  typedef struct _HashHitFilter HashHitFilter;
  
  /****************************************************************************
   *********************** Methods of type HashHitInfo ************************
   ****************************************************************************/ 
  
  HashHitInfo *hashCreateHitInfo (int blksz, const HashTable *htp);
  /**< Constructor
   */

  void hashDeleteHitInfo(HashHitInfo *p);
  /**< Destructor
   */

  int hashCollectHitInfo(HashHitInfo *hhip, unsigned char is_reverse, 
			 unsigned char basq_thresh,
			 SEQLEN_t seq_start, SEQLEN_t seq_end,
			 const SeqFastq *seqp, const HashTable *htp);
  /**< Collect stats on k-tuple hits from the read.
   * \param hhip Hit info containing stats.
   * \param is_reverse Flag indicating whether hits are to be collected on the
   *        forward (0) or reverse (1) strand 
   * \param basq_thresh Threshold of the base qualities (based to 0). K-mer words 
   *        with bases below this threshold are disregarded.
   * \param seq_start Specifies the beginning of a segment of the read sqp from
   *        which the k-tuples should be collected. If seq_start >= seq_end, collect
   *        from the entire read.
   * \param seq_end Specifies the end of a segment of the read sqp from
   *        which the k-tuples should be collected. If seq_start >= seq_end, collect
   *        from the entire read.   
   * \param seqp Read on which the k-tuple stats are to be collected. Has to be
   *        encoded in SEQCOD_MANGLED type.
   * \param htp Hash table.
   */

  int hashCollectHitInfoShort(HashHitInfo *hhip, unsigned char is_reverse, 
			      HASHNUM_t maxhit_per_tuple,
			      HASHNUM_t maxhit_total,
			      unsigned char  basq_thresh,
			      const SeqFastq *seqp, const HashTable *htp);
  /**< Collect stats on k-tuples hits from the reads and select the more informative
   * k-tupels for use.
   * \param hhip Hit info containing stats.
   * \param is_reverse Flag indicating whether hits are to be collected on the
   *        forward (0) or reverse (1) strand
   * \param maxhit_per_tuple Cutoff of the maximum number of hits a k-mer word can have.
   *        k-mer words with more hits are not considered. If == 0, no cut-off is applied.
   * \param maxhit_total Cut-off in the total number of k-mer words
   *        hits of the read.  The more abundant k-mer words are removed
   *        so that the total number of hits is less or equal to the cut-off.
   * \param basq_thresh Threshold of the base qualities (based to 0). K-mer words 
   *        with bases below this threshold are disregarded.
   * \param seqp Read on which the k-tuple stats are to be collected. Has to be
   *        encoded in SEQCOD_MANGLED type.
   * \param htp Hash table.
   */
  
       
  int hashSortHitInfo(HashHitInfo *hhip);
  /**< Sort the k-mer words by increasing hit frequency.
   */

  uint32_t hashCalcHitInfoCoverDeficit(const HashHitInfo *hip);
  /**< Calculate the maximum cover deficit due to unsampled k-tuples
   */
  
  uint32_t hashCalcHitInfoNumberOfHits(const HashHitInfo *hip, 
				       HASHNUM_t maxhit_per_tuple);
  /**< Calculate the total number of k-tuple words sampled.
   * \param hip Hit info, i.e. k-mer word statistics.
   *
   * \param maxhit_per_tuple Cut-off of the number of k-mer word
   *        occurrances across the genomic reference seqeunces. K-mer
   *        words that occur more frequently are not counted. If 0, no
   *        cut-off is applied.
   */

  uint32_t hashHitInfoCalcHitNumbers(const HashHitInfo *hhip, uint32_t *nhit_rank);
  /**< Returns the total number of k-mer words seeds and the number of seeds actually used.
   * \return Total number of k-mer word hits.
   * \param hhip Hit info structure.
   * \param nhit_rank Returns the number of k-mer word hits actually used.
   */

  /****************************************************************************
   *********************** Methods of type HashHitList ************************
   ****************************************************************************/ 

  HashHitList *hashCreateHitList(int maxnhits);
  /**< Constructor.
   * \param maxnhits   Maximum total number of hits allowed (allocated) in the list 
   */
  void hashDeleteHitList(HashHitList *hlp);
  /**< Destructor 
   */
  void hashBlankHitList(HashHitList *hlp);
  /**< Resets the hit list without freeing already allocated memory, so that previous
   * contents are subsequently overwritten.
   */
  int hashCollectHits(HashHitList *hlp,
		      const SeqFastq *seqp, const HashTable *htp,
		      char is_reverse, int maxhit_per_tuple, 
		      unsigned char basq_thresh,
		      const HashHitFilter *hhfp);
  /**< Collect k-tuple hits in hit list and sort them. K-tuples are
   * sampled from all positions in the query sequence {0,1, ...,
   * qlen-k-1}.  Where qlen is the length of the query sequence and k
   * is the length of the k-tuple. Hits are sorted by shift and
   * (subsequently) by offset in the query sequence.
   *
   * \param hlp Hit list to be filled. The hit list is overwritten and, if necessary, memory
   *            is re-allocated.
   * \param seqp Sequence from which the k-tuples are sampled.
   * \param htp Hash table containing the k-tuple hits.
   * \param is_reverse Flags a hit list for the reverse complement of the sequences
   *                   in the hash table.
   * \param maxhit_per_tuple Maximum number of hits for a k-tuple to register in the hit list.
   * \param basq_thresh Threshold of the base qualities (based to 0). K-mer words 
   *        with bases below this threshold are disregarded.
   * \param hhfp Filter out hits in particular intervals.
   */
  
  int hashCollectHitsForSegment(HashHitList *hlp,
#ifdef RESULTS_TRACKER
				Track *trkp,
#endif
#ifdef hashhit_dump_sortarray
				FILE *dumpfp,
				int *dumpctr,
#endif
				SETSIZ_t segmoffs_lo, 
				SETSIZ_t segmoffs_hi,
#ifndef hashhit_minimise_coverdeficit
				HASHNUM_t nhit_max,
				unsigned char use_short_hitinfo,
#endif
				const HashHitInfo *hhip,
				const HashTable *htp,
				const HashHitFilter *hhfp);
   /**< Collect k-tuple hits from the HitInfo structure using a cutoff
   * in the number of hits a k-tuple is allowed to have across the
   * genome. Only hits that fall in a particular segment in the
   * concatenated genomic sequences are returned. This can be used,
   * for example, to retrieve hits for each successive chromosomal sequence. The segment
   * positions must be greater than in the previous call and the segment must not overlap
   * with the one of the previous call.
   * \param hlp Hit list to be filled. Contents are overwritten.
   * \param segmoffs_lo Offset (counting from 0) of the first base of the segment 
   *          for which hits should be retrieved.
   * \param segmoffs_hi  Offset (counting from 0) of the last base of the segment plus 1
   *          for which hits should be retrieved.
   * \param nhit_max Cut-off in the number of hits for a k-tuple. K-tuples with a number of hits
   *        greater than the cutoff are regareded as not informative and their positions are
   *        not sampled. If nhit_max == 0, no cut-off is applied.
   * \param use_short_hitinfo if > 0, use only the least frequent k-tuples in hhip. if == 0, use
   *        all k-mers as long as they have fewer than nhit_max hits.
   * \param hhip Hash hit info with the k-tuples of the query sequence.
   * \param htp Hast table.
   * \param hhfp Hit filter allowing the hits to be filtered further by a set of intervals.
   *
   * \note To retrieve hits for each chromosomal sequences. The function has to be called 
   *  successively for each reference sequence in the order in which they where used for the 
   *  construction of the hash table.
   */

  int hashCollectHitsUsingCutoff(HashHitList *hlp, 
#ifdef hashhit_dump_sortarray
				 FILE *dumpfp,
				 int *dumpctr,
#endif
				 HASHNUM_t max_nhit_per_tup, 
				 const HashTable *htp,
				 const HashHitInfo *hip);
  /**< Collect k-tuple hits from the HitInfo structure using a cutoff
   * in the number of hits a k-tuple is allowed to have across the
   * genome.
   */

  void hashPrintHitList(const HashHitList *hlp, FILE *fp);
  int hashCheckHitList(const HashHitList *hlp, const SeqFastq *seqp, 
		       const HashTable *htp, const SeqSet *ssp,
		       const SeqCodec *codep);

  const uint64_t *hashGetHitListData(int *nhits, char *is_reverse,
				     uint32_t *qlen,
				     unsigned char *ktup, 
				     unsigned char *nskip,
				     const char **qmask,
				     const HashHitList *hlp);
  /**< Return the raw k-tuple hits sorted by shift. Each hit is
   * represented by a 64-bit integer. The lower 32-bits represent the
   * offset qo in the query sequence. The upper 32-bits encode the shift as
   * if (is_reverse) ktup_s + qo/nskip 
   * else ktup_s - qo/nskip + shift_offs
   * where nskip is the skip step size and ktup is the ktuple serial number in 
   * the concatenated reference sequences. Hits are sorted in ascending order,
   * i.e. first by shift, then by offset in query sequence.
   *
   * \param nhits Returns the number of seeds.
   * \param is_reverse 0 if forward strand, 1 if reverse strand
   * \param qlen Lenght of the query sequence.
   * \param ktup K-mer word length.
   * \param nskip Skip step for equidistant sampling of k-mer words.
   * \param qmask Returns a sequence of HASH_HIT_QUALIFIERS which allows
   *              tracing of the fate of each k-tuple, i.e. whether it was
   *              excluded by the 'ncut' threshold etc.
   * \param hlp Hit list. 
   */
  uint32_t hashCalcHitListCoverDeficit(const HashHitList *hlp);
  /**< Calculate the maximum cover deficit caused by unsampled k-tuples
   */
  /****************************************************************************
   ********************** Methods of Type HashHitFilter ***********************
   ****************************************************************************/ 
  HashHitFilter *hashCreateHitFilter(short blocksiz);
  /**< Constructor 
   * \param blocksiz Blocks size for memory allocation (as the number of intervals).
   */

  void hashDeleteHitFilter(HashHitFilter *p);
  /**< Destructor 
   */
  void hashResetHitFilter(HashHitFilter *p);
  /**< Empty the hit filter
   */
  
  int hashAddToHitFilter(HashHitFilter *hhfp,
			 uint32_t lo, uint32_t hi);
  /**< Add a pair of distance intervals to the hit filter.
   * 
   * \param hhfp Hit filter to which the distance intevals are to be added.
   * \param lo Lower limit of the interval, absolute k-tuple number in
   *           the sequence of concatenated sequences (starting from 0).
   * \param hi Upper limit of the interval (k-tuple number);
   */

  void hashPruneHitFilter(HashHitFilter *hhfp);
  /**< Sort intervals in the hit filter and remove overlapping intervals
   */
  
#endif
#ifdef __cplusplus
}
#endif
