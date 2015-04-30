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

#ifdef __cplusplus
extern "C"
{
#endif

#ifndef RMAP_H
#define RMAP_H

  //#define rmap_debug
#define rmap_short_hitinfo
  //#define rmap_mscor_calib /* ouptut data for mapping score calibration */
  //#define rmap_stop_candlist_early /* stop aligning candidate segments ealy */ 
#define rmap_finehash_2ndmate

#include <stdio.h>

#include "elib.h"
#include "sequence.h"
#include "hashhit.h"
#include "results.h"
#include "resultpairs.h"

  /****************************************************************************
   ****************************** Constants ***********************************
   ****************************************************************************/
  enum RMAP_FLAGS {
    RMAPFLG_CMPLXW  = 0x01,   /**< Produce complexity weighted Smith-Waterman scores */
    RMAPFLG_BEST    = 0x02,   /**< Look for best alignment only (including second best for mapping score) */
    RMAPFLG_ALLPAIR = 0x04,   /**< Perform exhaustive search for pairs - whether or not
			       * there are multiple hits for the rare mate. */
    RMAPFLG_SPLIT      = 0x08,/**< Look for split reads (overrides RMAPFLG_BEST) */
    RMAPFLG_SEQBYSEQ   = 0x10,/**< Process each reference sequence individually */
    RMAPFLG_NOSHRTINFO = 0x20,/**< Don't use the short but the full hitinfo - includes
			       * a cutoff of the hit frequency */
    RMAPFLG_PAIRED     = 0x40,/**< Align a read pair */
    RMAPFLG_SENSITIVE  = 0x80, /**< Do all alignments where candiate segments are within
				* nskip of maximum. */
 };

  /****************************************************************************
   ******************************** Types *************************************
   ****************************************************************************/

  typedef uint16_t RMAPFLG_t; /**< Holds a combination of RMAP_FLAGS */

  /****************************************************************************
   ***************************** Opaque types *********************************
   ****************************************************************************/

  typedef struct RMap_ RMap;

  /****************************************************************************
   *************************** Methods of type RMap ***************************
   ****************************************************************************/ 

  RMap *rmapCreate(const HashTable *htp, const SeqCodec *codecp,
		   const SeqSet *ssp,
		   const ScoreMatrix *scormtxp, RMAPFLG_t rmapflg);
  /**< Constructor
   * \param htp Hash table.
   * \param codecp Sequence En-/Decoder.
   * \param ssp Set of reference sequences (indexed by htp).
   * \param scormtxp Matrix of alignment scores, must be != NULL if (rmapflg & RMAPFLG_CMPLXW) 
   * \param rmapflg Bit flag, combination of RMAPFLG_CMPLXW and RMAPFLG_PAIRED.
   */

  void rmapDelete(RMap *rmp);
  /**< Destructor
   */

  void rmapBlank(RMap *rmp);
  /**< Resets type to 'empty' state
   */

  void rmapGetData(const ResultSet **rslt_readp,
		   const ResultSet **rslt_matep,
		   const ResultPairs **pairp, 
		   SeqFastq **sbufAp, SeqFastq **sbufBp,
		   const RMap *rmp);
  /**< Accessor.
   *
   * \param rslt_readp Returns set of results for (1st) read (of a pair).
   * \param rslt_matep Returns set of results for (2nd) read.
   * \param pairp Returns Result pairs for paired-end reads.
   * \param sbufAp Sequence buffer.
   * \param sbufBp 2nd sequence buffer.
   * \param rmp Rmap type.
   */

  void rmapSetFilterThresholds(RMap *rmp, int min_swatscor, int min_swatscor_below_max,
			       double min_identity);
  /**< Acessor. Sets thresholds to filter alignments.
   * \param rmp Container.
   * \param min_swatscor Minimum SW score.
   * \param min_swatscor_below_max Minimum SW score relative to maximum.
   * \param min_identity Minimum identity for an alignment to be reported. Fract < 1.0 as
   *        fraction of entire read. >= 1.0 as integer number int(id_abs) of bases.
   */

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
		 int ktuple_maxhit, unsigned int min_cover,
		 int min_swatscor,
		 int min_swatscor_below_max, 
		 unsigned char min_basqval,
		 short target_depth, short max_depth,
		 RMAPFLG_t rmapflg,
		 const ScoreMatrix *scormtxp,
		 const ResultFilter *rsfp,
		 const HashTable *htp, const SeqSet *ssp, const SeqCodec *codecp);
  /**< Map a single read onto a set of reference sequences.
   *
   * \param errmsgp Returns error messages.
   * \param rmp RMap structure containing various buffers.
   * \param readp Query reads, must be encoded as SEQCOD_MANGLED on input.
   *              Is returned decoded as SEQCOD_ASCII.
   * \param ktuple_maxhit Cutoff for the number of hits for a k-mer word in the query sequence.
   *        k-mer words with hit frequencies above that number are discounted (as uninformative).
   * \param min_cover Minimum number of bases in a segment covered by k-mer word hits.
   * \param min_swatscor Minimum Smith-Waterman score for an alignment to be counted.
   * \param min_swatscor_below_max Minimum Smith-Waterman score for an alignmen to be counted.
   *        This is given relative to the maximum Smith-Waterman score of the query read.
   * \param min_basqval Threshold of the base quality (based to 0). K-mer words with bases of 
   *        lower quality are disregarded.
   * \param target_depth Target number of potentially matching pairwise segments that are passed
   *        on to the dynamic programming step. If there a more potentially matching segments than
   *        this number, segments are aligned starting with the segments of the highest coverage until
   *        number target_depth is reached. If the difference to the maximum depth is less than the 
   *        k-word sampling step size, up to max_depth segments may be aligned.
   * \param max_depth Maximum number of potentially matching pairwise segments that are passed
   *        on to the dynamic programming step.
   * \param rmapflg Combination of RMAP_FLAGS.
   * \param scormtxp Matrix of SW-alignment scores.
   * \param rsfp Data like a minimum SW score to filter final list of results.
   * \param htp Hash table.
   * \param ssp Set of reference sequences.
   * \param codecp Sequence En-/Decoder.
   */

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
	       int ktuple_maxhit, 
	       unsigned int mincov_read, unsigned int mincov_mate,
	       int min_swatscor, unsigned char min_basqval,
	       short target_depth, short max_depth, 
	       RMAPFLG_t rmapflg, 
	       const ScoreMatrix *scormtxp,
	       const ResultFilter *rsfp,
	       const HashTable *htp, 
	       const SeqSet *ssp, 
	       const SeqCodec *codecp);
  /**< Map a pair of reads against a set of reference sequences.
   *
   * \param errmsgp Returns error messages.
   * \param rmp RMap structure
   * \param readp 1st mate of the paired-end query read, must be encoded as SEQCOD_MANGLED on input.
   *              Is returned decoded as SEQCOD_ASCII.
   * \param matep 2nd mate of query read, must be encoded as SEQCOD_MANGLED on input.
   *              Is returned decoded as SEQCOD_ASCII.
   * \param pairflgp Returns binary flags, a combination of RSLTPAIR_FLAGS.
   * \param d_min Minimum insert size of pair.
   * \param d_max Maximum insert size of pair.
   * \param pairlibcode Type of read pair library (one of RSLTPAIR_LIB).
   * \param ktuple_maxhit Cutoff for the number of hits for a k-mer word in the query sequence.
   *        k-mer words with hit frequencies above that number are discounted (as uninformative).
   * \param mincov_read Minimum number of bases a segment of the 1st mate that must be covered 
   *              by k-mer word hits.
   * \param mincov_mate Minimum number of bases a segment of the 2nd mate that must be covered 
   *              by k-mer word hits.
   * \param min_swatscor Minimum Smith-Waterman score for an alignment to be considered.
   * \param min_basqval Threshold of the base quality (based to 0). K-mer words with bases of 
   *        lower quality are disregarded.
   * \param target_depth Target number of potentially matching pairwise segments that are passed
   *        on to the dynamic programming step. If there a more potentially matching segments than
   *        this number, segments are aligned starting with the segments of the highest coverage until
   *        number target_depth is reached. If the difference to the maximum depth is less than the 
   *        k-word sampling step size, up to max_depth segments may be aligned.
   * \param max_depth Maximum number of potentially matching pairwise segments that are passed
   *        on to the dynamic programming step.
   * \param rmapflg Combination of RMAP_FLAGS.
   * \param scormtxp Matrix of SW-alignment scores.
   * \param rsfp Data like a minimum SW score to filter final list of results.
   * \param htp Hash table.
   * \param ssp Set of reference sequences.
   * \param codecp Sequence En-/Decoder.
   */
#ifdef RESULTS_BASQUAL
  int rmapPrintObservedBasQ(const RMap *rmp);
  /**< Return phred-scaled base qualties observed from alignments. */
#endif
  
#endif
#ifdef __cplusplus
}
#endif
  
