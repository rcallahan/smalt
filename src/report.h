/** Output of alignment results for a template (read) */

/*****************************************************************************
 *****************************************************************************
 *                                                                           *
 *  Copyright (C) 2012 - 2014 Genome Research Ltd.                           *
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

#ifndef REPORT_H
#define REPORT_H

#include <stdio.h>
#include "sequence.h"
#include "diffstr.h"
#ifdef RESULTS_TRACKER
#include "tracker.h"
#endif

  //#define report_debug

  /****************************************************************************
   ******************************** Constants *********************************
   ****************************************************************************/
  enum REPORT_OUTPUT_FORMATS {
    REPORTFMT_CIGAR,  /**< CIGAR line output format */
    REPORTFMT_SSAHA,  /**< SSAHA2 output format */
    REPORTFMT_SAM,    /**< SAM (samtools) output format */
    REPORTFMT_BAM,    /**< BAM ouptut format */
    REPORTFMT_GFF2,   /**< General Feature Format version 2 */
  };

  enum REPORT_MODIFIER_FLAGS {
    REPORTMODIF_ALIOUT   = 0x01,  /**< Produce explicit alignment output */
    REPORTMODIF_SOFTCLIP = 0x02,  /**< Soft clipping (for SAM output) */
    REPORTMODIF_HEADER = 0x04,    /**< Include header (for SAM output) */
    REPORTMODIF_XMISMATCH = 0x08, /**< Produce a cigar string with 'X' for mismatch */
   };

  enum REPORT_MATE_FLAGS {
    REPMATEFLG_MAPPED = 0x01,  /**< Read was mapped */
    REPMATEFLG_REVERSE = 0x02, /**< Read was mapped as reverse complement */
    REPMATEFLG_PAIRED = 0x04,  /**< Read was part of a pair */
    REPMATEFLG_2NDMATE = 0x08, /**< Read was 2nd mate of a pair (implies REPMATEFLG) */
    REPMATEFLG_PRIMARY = 0x10, /**< Primary alignment of a read (implies  REPMATEFLG_PARTIAL not set) */
    REPMATEFLG_PARTIAL = 0x20, /**< Partial secondary alignment of a read (implies REPMATEFLG_PRIMARY not set) */
    REPMATEFLG_MULTI = 0x40,   /**< Multiple possible pair/read placements (may be Reported with a random assignment) */
  };
  
  enum REPORT_PAIR_FLAGS {  /**< flags for mapping categories */
    REPPAIR_MAPPED = 0x01,  /**< Reads are mates of a pair and both mates are mapped */
    REPPAIR_CONTIG = 0x02,  /**< Both mates mapped to the same reference sequence */
    REPPAIR_PROPER = 0x04,  /**< Both mates mapped in proper orientation */
    REPPAIR_WITHIN = 0x08,  /**< Both mates mapped with in the distance expected from the insert size */
  };

  /****************************************************************************
   ****************************** Simple Types ********************************
   ****************************************************************************/

  typedef unsigned char REPOUFMT_t;     /**< Holds one of REPORT_OUTPUT_FORMATS */
  typedef unsigned char REPMODIFLG_t;   /**< Holds a combination of REPORT_MODIFIER_FLAGS */
  typedef unsigned char REPMATEFLG_t;   /**< Holds a combination of REPORT_MATE_FLAGS */
  typedef unsigned char REPPAIRFLG_t;   /**< Holds a combination of REPORT_PAIR_FLAGS */
  

  /****************************************************************************
   ****************************** Opaque Types ********************************
   ****************************************************************************/

  typedef struct _Report Report;
  typedef struct _ReportWriter ReportWriter;

  /****************************************************************************
   ********************** Methods of Type ReportWriter ************************
   ****************************************************************************/
  ReportWriter *reportCreateWriter(int *errcode,
				   const char * const filnam,
				   const REPOUFMT_t outform, 
				   const REPMODIFLG_t modiflg,
				   const SeqSet *ssp,
				   const char *prognam,
				   const char *progversion,
				   char * const *cmdlin_argv,
				   int cmdlin_narg
				   );
  /**< Constructor.
   * \param errcode Returns error code (can be NULL).
   * \param filnam File name (can be NULL which directs to STDOUT).
   * \param outform Output format (one of REPORT_OUTPUT_FORMATS)
   * \param modiflg Format specific flags (combination of REPORT_MODIFIER_FLAGS)
   * \param ssp Set of reference sequences (needed for SAM header & BAM, can
   *        otherwise be NULL).
   * \param prognam Name of the main program.
   * \param progversion Program version.
   * \param cmdlin_argv Pointer to original command line arguments.
   * \param cmdlin_narg Number of command line arguments.
   */

  void reportDeleteWriter(ReportWriter *p);
  /**< Destructor.
   */

  FILE *reportGetWriterStream(const ReportWriter *p);
  /**< Returns the output stream of the writer
   */

  /****************************************************************************
   ************************* Methods of Type Report ***************************
   ****************************************************************************/
 
  Report *reportCreate(int blksz);
  /**< Constructor.
   * \param blksz Block size (as the number of alignments) for memory allocation
   *              If <=0 use default block size.
   */

  void reportDelete(Report *p);
  /**< Destructor.
   */

  void reportBlank(Report *p);
  /**< Reset Report structure
   */

  int reportNextPairID(Report *rep);
  /**< Return an ID >=0 for the next read pair (to be used with reportAddMap).
   * Return < 0 on error.
   */

  int reportAddMap(Report *p, 
		   int pairid,
		   int swatscor, short mapscor,
		   SEQLEN_t q_start, SEQLEN_t q_end,
		   SEQLEN_t s_start, SEQLEN_t r_end, SEQNUM_t s_idx,
		   const DIFFSTR_T *dstrp, int dfslen,
		   int insiz,
		   REPMATEFLG_t flg, REPPAIRFLG_t pairflg);
  /**< Add an alignment result to the Report.
   * \param p Report strucure holding alignment results.
   * \param pairid Add the the alignment result to the pair with the ID pairid.
   *        If < 0, add a single read that is not part of a pair.
   * \param swatscor Smith-Waterman alignment score.
   * \param q_start Start of the aligned segment (base 1) of the query sequence
   * \param q_end End of the the aligned segment (base 1) of the query sequence
   * \param s_start Start of the aligned segment (base 1) of the reference sequence
   * \param s_end End of the the aligned segment (base 1) of the reference sequence
   * \param sidx Sequence inddex of the reference sequence in set of reference sequences
   * \param dstrp Compressed alignment string in the direction of the reference sequence.
   * \param insiz Insert size (as if for 1st mate, i.e. positive).
   * \param flg Comination of REPORT_FLAGS.
   */

  void reportFixMultiplePrimary(Report *rep);
  /**< If there are several pairs with the same read alignment
   * labelled as 'primary', remove these labels.  Use this to ensure
   * that at most one alignment is labelled the 'primary' alignment in
   * SAM/BAM output formats.
   */

  int reportWrite(const ReportWriter *wrp,
#ifdef RESULTS_TRACKER
		  const Track *read_trackp,
		  const Track *mate_trackp,
#endif 
		  const SeqFastq *readp, 
		  const SeqFastq *matep, 
		  const SeqSet *ssp,
		  const SeqCodec *codecp,
		  const Report *rep);
  /**< Write a set of alignment results
   * \param wrp Parameters, buffers and stream for output.
   * \param readp Read (1st mate of a read pair)
   * \param matep 2nd mate of a read pair (can be NULL)
   * \param ssp Set of reference sequences.
   * \param codec Sequence De-/Encoder.
   * \param rep Report to be output.
   */

#endif
#ifdef __cplusplus
}
#endif
