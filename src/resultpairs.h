/** Processing and output of alignment results from the same insert */

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

#ifndef RESULTPAIRS_H
#define RESULTPAIRS_H

  //#define results_debug
  //#define results_assert    /* self-check of alignment results */
  //#define results_mscor_calib /* output data for calibration of mapping score */

#define results_mapscor_exp    /* exponential scaling of Smith-Waterman score for mapping score */
  //#define results_mapscor_pair   /* 'proper' pair score */
  //#define results_pairscor_stats /* output statistics on pair scores & insert sizes */
  //#define results_loscor_capped

#include <stdint.h>

#include "insert.h"
#include "results.h"
#include "report.h"


  enum RSLTPAIR_FLAGS {
    RSLTPAIRFLG_PAIRED       = 0x01,  /**< Paired reads are mapped (as opposed to single reads) */
    RSLTPAIRFLG_RAREMATE     = 0x02,  /**< Read with fewer k-tuple hits is 1st (0) or 2nd (1) 
				       * read of the pair. This 'rare' mate is also the read 
				       * mapped first */
    RSLTPAIRFLG_RESTRICT_1st = 0x04,  /**< 1st read was mapped in area restricted by 2nd mate */
    RSLTPAIRFLG_RESTRICT_2nd = 0x08,  /**< 2nd read was mapped in area restricted by 1st mate */
    RSLTPAIRFLG_MATEPAIRLIB_1 = 0x10, /**< 1st bit of a 2-bit code for the orientation of mates */
    RSLTPAIRFLG_MATEPAIRLIB_2 = 0x20, /**< 2nd bit of a 2-bit code (RSLTPAIR_LIB) for the orientation of mates */
    RSLTPAIRFLG_INSERTSIZ     = 0x40, /**< insert size was inferred */
    RSLTPAIRFLG_MATEPAIRLIB_BIT = 4,  /**< Start bit of the above 2-bit code. */
    RSLTPAIRFLG_MATEPAIRLIB_MASK = 0x03, /**< For masking out 2-bit library code. */				 

  };

  enum RSLTPAIR_LIB {                 /**< Encodes the orientation of mates in different read pair libraries. */ 
    RSLTPAIRLIB_SINGLE = 0,           /**< Single reads, not paired reads */
    RSLTPAIRLIB_PAIREDEND = 1,        /**< Paired reads with mates on opposite strands, 5' ends
				       * on the outside, 3' ends facing each other (|--><--|) 
				       * as in Illumina paired-end libraries */
    RSLTPAIRLIB_MATEPAIR = 2,         /**< Paired reads with mates on opposite strands, 3' ends
				       * on the outside and the 5' end facing each other (<--| |-->)
				       * as in Illumina mate-pair libraries */
    RSLTPAIRLIB_SAMESTRAND = 3,       /**< Paired reads with mates on the same strand (|--> |-->)
				       * as in 454 libraries. */
    RSLTPAIRLIB_PAIREDALL = 4,        /**< Paired reads in all orientations (used for sampling)
				       */

  };

  /****************************************************************************
   ************************** Transparent Types *******************************
   ****************************************************************************/

  typedef uint8_t RSLTPAIRFLG_t; /**< Holds combination of RSLTPAIR_FLAGS bit flags */
  typedef uint8_t RSLTPAIRLIB_t; /**< Holds code of read pair library */

  /****************************************************************************
   ***************************** Opaque Types *********************************
   ****************************************************************************/

  typedef struct _ResultPairs ResultPairs;
  /**< List of paired read mappings 
   */

  /****************************************************************************
   *********************** Methods of Type ResultPairs ************************
   ****************************************************************************/

  ResultPairs *resultSetCreatePairs(short blksz);
  /**< Constructor.
   * \param blksz Blocks size for memory allocation
   */
  void resultSetDeletePairs(ResultPairs *p);
  /**< Destructor.
   */
  void resultSetBlankPairs(ResultPairs *p);
  /**< Resets type to empty state
   */

  int resultSetFindPairs(ResultPairs *pairp,
			 RSLTPAIRFLG_t pairflg,
			 RSLTPAIRLIB_t pairlibcode,
			 int dist_lo, int dist_hi,
			 const ResultSet *rsltAp, const ResultSet *rsltBp);
  /**< Find paired mappings amongst the top scoring results of both mates
   * irrespective of whether or not the inferred insert size is in the
   * expected range or the pairs are in proper orientation.
   * \param pairp Returns paired mappings.
   * \param pairflg Only bits RSLTPAIRFLG_RESTRICT_1st, RSLTPAIRFLG_RESTRICT_2nd are used.
   *        (indicate restricted mapping of one of the mates).
   * \param dist_lo lower limit of (inferred) insert size range.
   * \param dist_hi upper limit of (inferred) insert size range.
   * \param rsltAp Set of results for 1st read (query) of the pair.
   * \param rsltBp Set of results for 2nd read (mate) of the pair.
   */
  
  int resultSetFindProperPairs(ResultPairs *pairp, int dist_lo, int dist_hi,
			       int maxnum,
			       int swscor_min,
			       RSLTPAIRLIB_t pairlibcode,
			       const ResultSet *rsltAp, const ResultSet *rsltBp);
  /**< Find proper paired mappings among the results of two read mates.
   * \param pairp Paired mappings.
   * \param dist_lo lower limit of (inferred) insert size range.
   * \param dist_hi upper limit of (inferred) insert size range.
   * \param maxnum Maximum number of pairs to be collected.
   * \param swscor_min Minimum Smith-Waterman score for either of the pair.
   * \param pairlibcode Type of read pair library (one of RSLTPAIR_LIB).
   * \param rsltAp Set of results for 1st read (query) of the pair.
   * \param rsltBp Set of results for 2nd read (mate) of the pair.
   */
     
  int resultSetGetNumberOfPairs(int *n_proper, const ResultPairs *pairp);
  /**< Return the total number of pairs and the number of 'proper' pairs.
   * \param n_proper Returns the number of proper pairs (can be NULL).
   * \param pairp Contains all pairs gathered with resultSetFindPairs() 
   *              or resultSetFindProperPairs().
   */

  int resultSetAddPairToReport(Report *rep,
			       const InsHist *ihistp,
			       const ResultPairs *pairp,
			       RSLTPAIRFLG_t pairflg,
			       RESULTOUTFLG_t rsltouflg,
			       const ResultSet *rsrp,
			       const ResultSet *rsmp);

  /**< Output set of results for paired reads in one of a variety of output formats.
   * See also resultSetPrint().
   *
   * \param rep Report to which results should be added.
   * \param ihistp Histogram of insert sizes (can be NULL).
   * \param pairp Set of paired mappings generated, e.g. by resultSetFindBestPairs() or
   *        resultSetFindProperPairs().
   * \param pairflg Flag for paired reads (combination of RSLTPAIR_FLAGS).
   * \param rsltouflg a combination of RESULTSET_OUTPUT_FLAGS
   * \param rsrp Set of results for first read (i.e. 1st mate of paired read as
   *        determined by the order of the input files.
   * \param rsmp Set of results for 2nd read (mate).
   */
  
#endif
#ifdef __cplusplus
}
#endif
