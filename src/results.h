/** Processing and output of (pairwise) alignment results */

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

#ifndef RESULTS_H
#define RESULTS_H

  //#define results_debug
  //#define results_assert    /* self-check of alignment results */
  //#define results_mscor_calib /* output data for calibration of mapping score */

#define results_mapscor_exp    /* exponential scaling of Smith-Waterman score for mapping score */
  //#define results_mapscor_pair   /* 'proper' pair score */
  //#define results_pairscor_stats /* output statistics on pair scores & insert sizes */
  //#define results_loscor_capped

#include <stdint.h>

#include "sequence.h"
#include "hashidx.h"
#include "alignment.h"
#include "insert.h"
#ifdef RESULTS_TRACKER
#include "tracker.h"
#endif
#include "report.h"

  enum RESULTSET_OUTPUT_FLAGS {
    RESULTFLG_BEST         = 0x01,   /**< Look for best mappings only */
    RESULTFLG_SINGLE       = 0x02,   /**< In case of multiple best mappings, don't report all 
				      * of them (if RESULTFLG_RANDSEL is not set). */
    RESULTFLG_SPLIT        = 0x04,   /**< Look also for secondary alignments if they complement
				  * primary alignment (split reads) */
    RESULTFLG_RANDSEL      = 0x08,   /**< In case of multiple best mappings select exactly 
				      * 1 mapping at random */
  };

  enum RESULT_STATUS_FLAGS {
    RSLTFLAG_SELECT  = 0x01,   /**< indicates this result is selected */
    RSLTFLAG_RAW     = 0x02,   /**< raw alignment data: offsets relative to segment */
    RSLTFLAG_REVERSE = 0x04,   /**< indicates reverse complement match */
    RSLTFLAG_NOSEQID = 0x08,   /**< No Sequence index has been assigned and reference 
				* offset is absolute position in the sequence of
				* concatenated sequences. */
    RSLTFLAG_NOOUTPUT = 0x10,  /**< Read has been filtered out for output */
    RSLTFLAG_BELOWRELSW = 0x20,/**< Smith-Waterman score is below relative threshold */
    RSLTFLAG_HASSECOND = 0x40, /**< Has a secondary algnment (split read) */
    RSLTFLAG_PARTIAL   = 0x80, /**< Is the secondary partial alignment of a split read */
    RSLTFLAG_SINGLE    = 0x100,/**< Is the only result obtained */
    RSLTFLAG_REPORTED  = 0x200,/**< Was added to report (used to avoid adding result
				* as 2ndary as well as primary alignment) */
  };

  enum RESULT_PAIRMAP_FLAGS { /**< Flags for paired mapping */
    RSLTPAIRMAPFLG_REVERSE_1st = 0x01, /**< 1st read is mapped on reverse strand */
    RSLTPAIRMAPFLG_REVERSE_2nd = 0x02, /**< 2nd read is mapped on reverse strand */
    RSLTPAIRMAPFLG_SAMECONTIG = 0x04,  /**< Both reads are on the same contig 
					* (RSLTPAIRMAPFLG_NOCONTIG is unset) */
    RSLTPAIRMAPFLG_LEFTMOST2nd = 0x08, /**< Left most base (reference 5' end) of the alignment
					* is convered only by a base of the 2nd mate and not
					* by the first mate. */
    RSLTPAIRMAPFLG_NOCONTIG = 0x10,    /**< If set, RSLTPAIRMAPFLG_SAMECONTIG could not be
					* determined (e.g. because contigs were not assigned) 
					* (implies RSLTPAIRMAPFLG_SAMECONTIG is unset) */
  };

  enum RESULTSET_CONST {
    RESULTSET_UNKNOWN_SEQIDX = -1, /* indicates that sequence index
				    * hasn't yet been determined (must be < 0!)*/
  };

  enum RESULTSET_CALLBACKRTN { /**< Return types for callback function */
    RSLTSCBRTN_OK = 0,
    RSLTSCBRTN_BREAK = 1,  /**< Break the current loop */
    RSLTSCBRTN_STOP = 2,   /**< Break and fall out of all loops w/o error */
  };

  enum RESULTSET_SAMVERSIONS {
    RSLTSAMSPEC_V1P0 = 0,  /**< SAM spec version 1.0 */
    RSLTSAMSPEC_V1P4 = 1,  /**< SAM spec version 1.4 */
  };

  /****************************************************************************
   ************************** Transparent Types *******************************
   ****************************************************************************/

  typedef uint8_t RESULTOUTFLG_t; /**< Holds a combination of RESULT_OUTPUT_FLAGS */
  typedef uint16_t RSLTFLG_t;     /**< Holds a combination of RESULT_STATUS_FLAGS, 
				   * definition of INTERVALFLG_t must be identical */
  typedef uint8_t RSLTPAIRMAPFLG_t; /**< Holds a combination of RESULT_PAIRMAP_FLAGS */
  /****************************************************************************
   ***************************** Opaque Types *********************************
   ****************************************************************************/

  typedef struct _RESULT Result;
  /**< A single alignment result */

  typedef struct _ResultSet ResultSet;
  /**< A list of alignment results 
   */
  typedef struct _ResultFilter ResultFilter;
  /**< Filter for output of results 
   */
  typedef struct _ResultBuffer ResultBuffer;
  /**< For buffering pointers
   */

  /****************************************************************************
   ******************************** Callbacks *********************************
   ****************************************************************************/

  typedef int (ResultSetCallBackf) (int *errcode, void *arg, const Result *rp);
  /**< Hast to return one of RESULTSET_CALLBACKRTN.
   * \param errcode returns one of ERROR_CODE. Check this if 
   *        return value of call back is != RSLTSCBRTN_OK.
   * \param arg Arguments to call back function.
   * \param Pointer to current Result.
   */

  /****************************************************************************
   ********************************* Methods **********************************
   ****************************************************************************/

  short resultConvertProbabilityToMappingScore(double p);
  /**< Return a mapping quality score */

  /****************************************************************************
   ************************* Methods of Type Result ***************************
   ****************************************************************************/
  
  int resultGetData(SEQLEN_t *qs, SEQLEN_t *qe,
		    SEQLEN_t *rs, SEQLEN_t *re, SEQNUM_t *rx,
		    int *swscor, RSLTFLG_t *flag,
		    const Result *rp);
  /**< Accessor.
   * \param qs Returns the (1 based!) offset of the start of the aligment 
   * in the query read (can be NULL).
   * \param qe   Returns the (1 based!) offset of the end of the aligment 
   *                in the query read (can be NULL).
   * \param rs Returns the (1-based!) offset of the start of the aligment 
   * in the reference sequence (can be NULL).
   * \param re   Returns the (1-based) offset of the end of the aligment 
   *             in the reference sequence (can be NULL).
   * \param rx Returns the (0-based) index of the reference sequence.
   * \param swscor Returns Smith-Waterman score.
   * \param flag Returns a combination of RESULT_STATUS_FLAGS.
   * \param rp Result from which to obtain data.
   */
  
  short resultGetFragmentNo(const Result *rp);
  /**< Return the id of the query fragment of this result
   */

  short resultGetSWRank(const Result *rp);
  /**< Return the rank of the Smith-Waterman score
   */

  RSLTFLG_t resultGetStatusFlag(const Result *rp);
  /**< Returns a combination of RESULT_STATUS_FLAGS.
   */

  int resultGetMapQualScore(double *prob, RSLTFLG_t *flag, const Result * const rp);
  /**< Get the mapping quality score.
   * \return Mapping quality score phred scaled.
   * \param prob Mapping quality score likelihood of a correct mapping.
   * \param flag Returns a combination of RESULT_STATUS_FLAGS.
   * \param rp Result from which score should be fetched.
   */

  RSLTPAIRMAPFLG_t resultCalcInsertSize(int *isiz, unsigned char samspec,
					const Result *ap, const Result *bp);
  /**< Calculate insert size for two results according to SAM spec 1.0
   *  and return a flag that determines the orientation of the reads. 
   * \return A combination of RESULT_PAIRMAP_FLAGS (see below).
   * \param isiz Returns insert size (can be NULL).
   * \param samspec Version of the SAM specification (one of RESULTSET_SAMVERSIONS)
   * \param ap Alignment result for first mate.
   * \param bp Alignment result for 2nd mate.
   * \note: The sign of the insert size in SAM Spec 1.4 is determined by the position
   * of the respective insert (fragment) in a series of fragments (not determined
   * by this routine). But we use the convention that the sign is negative if the
   * left-most (reference-5') position of the alignment is covered by the 2nd-mate alone.
   * 
   * The bits determining the orientation are RSLTPAIRMAPFLG_REVERSE_1st and 
   *         RSLTPAIRMAPFLG_REVERSE_2nd
   * RSLTPAIRMAPFLG_SAMECONTIG is set if both alignments are on the same contig.
   * RSLTPAIRMAPFLG_NOCONTIG is set if the routine could not check whether the mapping 
   * results are on the same sequence or not (RSLTPAIRMAPFLG_SAMECONTIG is then unset).
   *
   *  pair orientation : [RSLTPAIRMAPFLG_REVERSE_1st:RSLTPAIRMAPFLG_REVERSE_2nd | 
   *                     sign of isiz (0 pos, 1 neg):RSLTPAIRMAPFLG_LEFTMOST2nd]
   *  |--a-->  <--b--| : 1+ [0:1 | 0:0] : PAIRMAPFLG_REVERSE_2nd set
   *  |--------------> : inferred insert size (+) (SAM spec v1.0)
   *  |--------------| : inferred insert size (SAM spec v1.4)
   *
   *  <--a--|  <--b--| : 3+ [11 | 00]
   *        |--------> : inferred insert size (+) (SAM-1.0)
   *  |--------------| : inferred insert size (SAM-1.4)
   *
   *  |--a-->  |--b--> : 0+ [00 | 00]
   *  |-------->       : inferred insert size (+) (SAM-1.0)
   *  |--------------| : inferred insert size (SAM-1.4)
   *
   *  <--a--|  |--b--> : 2+ [10 | 00]
   *        |-->       : inferred insert size (+) (SAM-1.0)
   *  |--------------| : inferred insert size (SAM-1.4)
   *  
   *  |--b-->  <--a--| : 2- [10 | 11]: PAIRMAPFLG_REVERSE_1st set
   *  <--------------| : inferred insert size (-) (SAM spec v1.0)
   *  |--------------| : inferred insert size (SAM spec v1.4)
   *
   *  <--b--|  <--a--| : 3- [11 | 11]
   *        <--------| : inferred insert size (-) (SAM-1.0)
   *  |--------------| : inferred insert size (SAM-1.4)
   * 
   *  |--b-->  |--a--> : 0- [00 | 11]
   *  <--------|       : inferred insert size (-) (SAM-1.0)
   *  |--------------| : inferred insert size (SAM-1.4)
   *
   *  <--b--|  |--a--> : 1- [01 | 11]
   *        <--|       : inferred insert size (-) (SAM-1.0)
   *  |--------------| : inferred insert size (SAM-1.4)
   *
   * Overlap:
   *     |--a-->
   *  <--b--|           : 1+ [01| 01]
   *     |-->           : inferred insert size (+) (SAM-1.0)
   *  |--------|        : inferred insert size (SAM-1.4)
   *
   *  |--a-->           :
   *  <--b--|           : 1+ [01| 00]
   *  |----->           : inferred insert size (SAM-1.0)
   *  |-----|           : inferred insert size (SAM-1.4)
   * 
   */

  RSLTPAIRMAPFLG_t resultCalcInsertSizeSamSpec1p4(int *isiz, const Result *ap, const Result *bp);
  /**< Calculate size of th insert according to the SAM/BAM spect 1.4
   * The insert size is always >=0.
   *
   *  
   *  pair orientation : RSLTPAIRMAPFLG_REVERSE_1st|RSLTPAIRMAPFLG_REVERSE_2nd|RSLTPAIRMAPFLG_MATEX
   *  |--a-->  <--b--| : 1+ [01 | 0] : PAIRMAPFLG_REVERSE_2nd+
   *  |--------------| | inferred insert size (+)
   *
   *  <--a--|  <--b--| : 3+ [11 | 0]
   *  |--------------| : inferred insert size  (+)
   *
   *  |--a-->  |--b--> : 0+ [00 | 0]
   *  |--------------| : inferred insert size (+)
   *
   *  <--a--|  |--b--> : 2+ [10 | 0]
   *  |--------------| : inferred insert size  (+)
   *  
   *  |--b-->  <--a--| : 2- [10 | 1]: PAIRMAPFLG_REVERSE_1st-
   *  |--------------| : inferred insert size (-)
   *
   *  <--b--|  <--a--| : 3- [11 | 1]
   *  |--------------| : inferred insert size (-)
   * 
   *  |--b-->  |--a--> : 0+ [00 | 1]
   *  |--------------| : inferred insert size (-)
   *
   *  <--b--|  |--a--> : 1- [01 | 1]
   *  |--------------| : inferred insert size (-)
   *
   *
   * Overlap:
   *     |--a-->
   *  <--b--|           : 1- [01| 1]
   *  |--------|        :
   *
   *  |--a-->           : 1+ [01| 0]
   *  <--b--|
   * 
   * RSLTPAIRMAPFLG_MATEX is set if the leftmost aligned point belongs to b
   */


  /* /\**************************************************************************** */
  /*  ********************** Methods of Type ResultBuffer ************************ */
  /*  ****************************************************************************\/ */

  /* ResultBuffer *resultSetCreateBuffer(void); */
  /* /\**< Constructor. */
  /*  *\/ */

  /* void resultSetDeleteBuffer(ResultBuffer *p); */
  /* /\**< Destructor */
  /*  *\/ */

  /****************************************************************************
   ************************ Methods of Type ResultSet *************************
   ****************************************************************************/

  ResultSet *resultSetCreate(int blocksiz, int blocksiz_diffstr);
  /**< Constructor
   * \param blocksiz Memory (re-)allocation in blocks of this size (bytes)
   * \param blocksiz_diffstr Block size for memory (re-) allocation of diff str
   */

  void resultSetDelete(ResultSet *rsp);
  /**< Destructor */
  
  void resultSetAlignmentStats(ResultSet *rsp, 
			       int n_ali_done, int n_ali_tot, short maxdepth,
			       uint32_t n_hits_used, uint32_t n_hits_tot);
  /**< Accessor.
   * Sets the number of alignments performed and the number of alignments that would
   * have to be considered in order to guarantee the best hit is found.
   * \param rsp Set of results.
   * \param n_ali_done Number of aligments actually performed.
   * \param n_ali_tot Total number of alignments that would have to be considered.
   * \param maxdepth Threshold in then number of alignments carried out.
   * \param n_hits_used Number of seed hits considered.
   * \param n_hits_tot Total number of seed hits.
   */

  unsigned char resultSetAlignmentWasCurtailed(const ResultSet *rsp);
  /**< Returns !=O if not all neccessary algnments were performed and
   * therefore a best score might have been missed.
   */

  int resultSetAddFromAli(ResultSet *rsp, const AliRsltSet *arsp,
			  SEQLEN_t soffs,
			  SEQLEN_t qoffs, SEQLEN_t qlen,
			  SEQNUM_t seqidx,
#ifdef results_debug
			  uint32_t sgx,
			  SeqFastq *sqbufp,
			  const char *profiled_seqp,
			  const SeqCodec *codecp,
			  const SeqSet *ssp,
#endif
			  char is_reverse);
  /**< Add alignment results - which are relative to the aligned segments -
   * with absolute offsets for the original sequences.
   * \param rsp Set of results.
   * \param arsp Set of original alignment results
   * \param soffs Offset (counting from 0) of the segment in the subject 
   *        sequence that was submitted to the alignment unprofiled.
   * \param qoffs Offset (counting from 0) of the segment in the query 
   *        sequence that was submitted with its profile to the alignment.
   * \param qlen Length of the query sequence.
   * \param seqidx Reference sequence index, RESULTSET_UNKNOWN_SEQIDX if not
   *        determined.
   * \param is_reverse Flag indicating that the subject segment was aligned as
   *        reverse complement.
   *
   * \note Offsets are on the forward strand (counting from 0).
   */
			    
  int resultSetAddMisMatch(ResultSet *rsp, const int mmoffs[], int mmnum, 
			   SEQLEN_t qoffs, int qlen, SEQLEN_t soffs, SEQNUM_t sidx, 
			   char is_rcpl, int swatscor);
  /**< Add a mapping possibly containing a number of mismatches.
   * \param rsp result set to which the mapping should be added.
   * \param mmoffs Array of mismatch offsets in ascending order. 
   *        Offsets are relative to qoffs starting from 1.
   * \param mmnum Number of mismatches (size of array).
   * \param qoffs Offset of matching segment in query sequence.
   * \param qlen Length of matching segmenet in query sequence.
   * \param soffs Offset of matching segment in reference sequence 
   *        (of concatenated sequences if sidx < 0)
   * \param sidx Sequence index. If sidx == RESULTSET_UNKNOWN_SEQIDX, 
   *        sequence index is unknown and soffs is the offset in the concatenated 
   *        set of sequences.
   * \param is_rcpl Indicates whether (==1) or not (==0) mapping is on reverse complement.
   * \param swatscor is the Smith-Waterman score for the alignment.
   */
  
  void resultSetUpdateFromSegment(ResultSet *rsp, short start_idx, short end_idx,
				  SEQLEN_t soffs, SEQLEN_t qoffs, char is_reverse);
  /**< Update raw sequence offsets (relative to segment) with the 
   * absolute offsets of the segment 
   */
 
  int resultSetSortAndAssignSequence(ResultSet *rsp, 
				     SeqFastq *sbufp,
				     unsigned char search_split,
				     const SeqFastq *sqp,
#ifdef results_debug
				     const SeqFastq *sqRCp,
#endif
				     const ScoreProfile *scpp,
				     const ScoreProfile *scpRCp,
				     const SeqSet *ssp,
				     const SeqCodec *codecp);
  /**< Assign reference sequence index where necessary and prepare for output.
   * \param rsp Result set.
   * \param sbufp Sequence buffer.
   * \param search_split Flag if != 0, look for secondary alignments that
   *        complement the primary alignment (split reads).
   * \param sqp Query (profiled) sequence providing base qualities for mapping score.
   * \param scpp Score profile (query sequence).
   * \param scpRCp Score profile for reverse complement of query sequence.
   * \param ssp Set of reference sequences
   * \param codecp Sequence Ent-/Decoder.
   */

  void resultSetBlank(ResultSet *rsp);
  /**< Reset result set for re-use 
   */

  int resultSetGetNumberOfSegments(short *nres, short *nseg, const ResultSet *rsp);
  /**< Accessor for number of results and number of segments
   * \return ERROR_CODES
   * \param nres Returns the number of selected results.
   * \param nseg Returns the number of segments.
   */

  int resultSetGetNumberOfResultsInSegment(int segx, const ResultSet *rsp);
  /* Accessor for the number of results in segment segx.
   * \return the Number of results.
   * \param segx Index (serial number) of segment.
   * \param rsp Set of results.
   */

  int resultSetGetResultInSegment(const Result **rpp, int segx, int resx, const ResultSet *rsp);
  /**< Accessor for result number resx in segment segx. Returns ERRCODE_FAILURE if
   * there are no alignment results.
   * \param rpp Returns pointer to result. Returns NULL if there are no results found.
   * \param segx Serial number of the query segment.
   * \param resx Serial number of the result in the segment segx.
   * \param rsp Set of results (after calling resultSetSortAndAssignSequence)
   */

  int resultSetGetResultByRank(const Result **rpp, int rank, const ResultSet *rsp);
  /**< Accessor for result of rank resx sorted by Smith-Waterman score.
   * \param rpp Returns pointer to result.
   * \param rank Rank of the result (by descending SW-score).
   * \param rsp Set of results.
   */

  int resultSetGetMaxSwat(const ResultSet *rsp, int *maxswat_2nd);
  /**< Return the maximum Smith-Waterman score in the set.
   * \param rsp Set of results.
   * \param maxswat_2nd Returns second largest Smith-Waterman score (can be NULL)
   */

  Result *resultSetGetResultBySWrank(const ResultSet *rsp, short rank);
  /**< Return Result no rank in the list of results sorted by Smith-Waterman
   * score.
   * \param rsp Set of results.
   * \param rank Position in sorted list.
   */

  int resultSetDo(void *argp, ResultSetCallBackf *cbf, 
		  const ResultSet *rsp);
  /**< Execute a function for all Results in ResultSet.
   * \return one of ERROR_CODES.
   * \param argp Argument to call back function.
   * \param cbf Pointer to call back funcion.
   * \param rsp Set of results.
   */
  
  int resultSetAddResultToReport(Report *rep,
				 int pairid,
				 short mapscor,
				 REPMATEFLG_t mateflg,
				 REPPAIRFLG_t pairflg,
				 int isize,
				 const Result *rp,
				 const ResultSet *rsp);
  /**< Add a single result to a Report.
   * \param rep Report to which result should be added.
   * \param pairid If >= 0, result is part of a pair of ID paird.
   *               If < 0 a single read is added.
   * \param mapscor Mapping score. Overwrite mapping score of Result rp
   *        with this score if adding mate of a pair (pairid >= 0). 
   * \param mateflg A combination of binary flags in REPORT_MATE_FLAGS.
   * \param pairflag A combinatio of binary flags in REPORT_PAIR_FLAGS
   *                 (ignored if pairid < 0).
   * \param isize Insert size (ignored if pairid < 0).
   * \param rp Result to be added.
   * \param rsp Result set to which result rp belongs 
   *             (caution: this is not checked at the moment!) 
   */
  int resultSetAdd2ndaryResultsToReport(Report *rep, 
					REPMATEFLG_t mateflg,
					RESULTOUTFLG_t rsltflg,
					const ResultSet *rsp);
  /**< Add secondary alignments (e.g. alternative mappings, segments of
   *   split alignments) to a report.
   * \param rep Report to which 2ndary results should be added.
   * \param mateflg A combination of binary flags in REPORT_MATE_FLAGS.
   * \param  rsltflg A combination of RESULTSET_OUTPUT_FLAGS.
   * \param rsp Set of results from which secondary alignment are to be added.
   */
  int resultSetAddToReport(Report *rep,
			   RESULTOUTFLG_t rsltflg,
			   const ResultSet *rsp);

  /**< Add set of results to Report.
   * \param rep Report to which set of results should be added.
   * \param rsltflg a combination of RESULTSET_OUTPUT_FLAGS
   * \param rsp Set of results
   */

  void resultSetPrintDebugInfo(FILE *fp, const ResultSet *rsp);
  /**< Print result set for debugging/testing
   */
  
  int resultSetGetScorStats(const ResultSet *rsp, 
			    int *scor_max, short *num_max, 
			    int *scor_2ndmax, short *num_2ndmax);
  /**< Return the number of alignment results and the number of alignments with 
   * best and 2nd best scores.
   * \param rsp Set of results.
   * \param scor_max Returns the maximum score (can be NULL)
   * \param num_max Returs the number of alignments with the maximum score (can be NULL)
   * \param scor_2ndmax Returns the 2nd best score (can be NULL)
   * \param num_2ndmax Returs the number of alignments with the 2nd best score (can be NULL)
   *
   * \note Results must have been sorted.
   */
  
  unsigned char resultSetGetRankDepth(const ResultSet *rsp, short *depth, short *rank);
  /** Return the number of top aligments to be considered.
   * \return 1 if unique best mapping, 0 otherwise.
   * \param rsp Set of results.
   * \param detph Returns number of top alignments to be considered (can be NULL).
   * \param rank Returns the maximum rank to be considered (can be NULL).
   * \note Results must have been sorted.
   */

  int resultSetGetMappingScore(const ResultSet *rsp, int *swscor);
  /**< Return the mapping score for the top mapping. This is 0 if there
   * are no mappings or if there are multiple best mappings with the same
   * Smith-Waterman score.
   * \param rsp Set of alignment results.
   * \param swscor Returns Smith-Waterman score for the alignment (can be NULL).
   */
  
  double resultSetGetMapQualAsProb(double *p2, short *n1, short *n2, const ResultSet * rsp);
  /**< Return top mapping score converted to likelihood.
   * \return likelihood of mapping being correct.
   * \param p2 Returns likelihood of 2nd best mapping being correct.
   * \param n1 Returns number of best mappings with equivalen Smith-Waterman scores.
   * \param n2 Returns number of 2nd best mappings (with equivalen Smith-Waterman scores).
   * \param rsp Set of alignment results.
   */

  int resultSetInferInsertSize(int *isiz, unsigned char samspec,
			       const ResultSet *rsrp, const ResultSet *rsmp);
  /**< Infer insert size if confident mapping. Returns ERRCOD_SUCCESS if 
   * a confident pair is 
   * found.
   *
   * \param isiz Returns insert size.  
   * \param samspec One of RESULTSET_SAMVERSIONS.
   * \param rsrp Result set for first mate.  
   * \param rsmp Result set for 2nd mate.
   */

  unsigned char resultSetTestOverlap(const ResultSet *rs1p, int overlap_percent, const ResultSet *rs2p);
  /**< Test whether top result in first set overlaps significantly
   * with best or 2nd best results in second set. 
   *
   * \param rs1p First set of results.
   * \param overlap as percentage of the length of the top aligned segment
   *        in the first set. 
   * \param rs2p Second set of results.
   */
  
  Result *resultSetGetTopResult(unsigned char *is_multi, 
				const unsigned char is_randsel, 
				const ResultSet *rsp);
  /**< Return pointer to top result. If there are are multiple top results
   * draw at random from those if is_randsel != 0, return < 0 if randsel == 0
   * \return Pointer to (NULL if not available)
   * \param is_multi Returns 0: if a single best value was available, 1 if there were
   *        multiple best values.
   * \param is_randsel == 0: only return index >= 0 if there is single best. !=0 make
   *        a random selection.
   * \param Set of alignment results.
   */

#ifdef RESULTS_BASQUAL
  int resultSetPrintBasqStats(FILE *fp, const ResultSet *rsAp, const ResultSet *rsBp);
  /**< Output phred-scaled base qualities ovserved during the alignment run.
   * \param fp Output stream
   * \param rsAp Set of results for single reads or first mates.
   * \param rsBp Set of results for second mates of paired reads (can be NULL).
   */
#endif
  

  /****************************************************************************
   ********************** Methods of Type ResultFilter ************************
   ****************************************************************************/
  
  ResultFilter *resultSetCreateFilter(void);
  /**< Constructor.
   */

  void resultSetDeleteFilter(ResultFilter *p);
  /**< Destructor
   */

  void resultSetFilterData(ResultFilter *p, int sw_abs, int sw_rel, double id_abs);
  /**< Accessor for seting filter data.
   * \param p Filter structure.
   * \param sw_abs Absolute threshold of the Smith-Waterman score.
   * \param sw_rel Threshold of the Smith-Waterman score relative to maximum.
   * \param id_abs Absolute threshold of the fraction of exact matching nucleotides;
   */
  
  int resultSetFilterResults(const ResultSet *rsp, 
			     const ResultFilter *rsfp, 
			     const SeqFastq *sqp);
  /**< Sort and filter results for output.
   * \param rsp Set of results to be filtered.
   * \param rsfp Filter data.
   * \param sqp Read sequence.
   */
 
#endif
#ifdef __cplusplus
}
#endif
