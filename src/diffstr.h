/** Handling of compressed alignment strings (pairwise alignments) */

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

/*****************************************************************************
 *************************** DIFFSTR FORMAT **********************************
 *****************************************************************************
 *
 * The pairwise alignment is partitioned into segments each of which
 * consists of m matches with 0 <= m <= DIFFSTR_MAXMISMATCH either
 * terminated by and additional 1-nucleotide match (M, i.e. a total of
 * m+1 exact matches) or terminated by a 1-nucleotide mismatch (S for
 * "substitution"), a 1-nucleotide insertion (I) or deletion (D).
 * S also designates the end of a run of matches at the end of
 * the alignment.
 * 
 * Insertions/deletions are with respect to the profiled sequence (query read).
 * I.e. an insertion is one additional nucleotide in the profiled sequence 
 * compared to the unprofiled sequence.
 *
 * Each segment is represented by an 8-bit character that encodes one
 * of those 4 states in the upper (most significant) 2 bits: M(00),
 * I(10), D(01), S(11).  The lower 6 bits encode the number m of exact
 * matches up to the position in the alignment characterise by the
 * state letter. I.e. the number of exactly matching bases n is 
 * n = m for I, D, S and n = m + 1 for M.
 * 
 * The 'diff string' is terminated by 0 (M:0). If the alignment
 * ends in a match, the last state is (S:m) followed by (M:0).
 *
 * Examples:
 * A Match of 100 bases is encoded as (DIFFSTR_MAXMISMATCH = 62)
 * diffs = M:61|S:38|M:0|
 *
 * Profiled sequence (query read)         : AGATCCCCAACAGAGTCCC-AGGATCACACAGAC
 * Unprofiled sequence (reference segment): AGATCCCCAACAGAGTACCAAGGATCACACAGAC
 * Aligment                               : AGATCCCCAACAGAGT CC-AGGATCACACAGAC
 * diffs = S:16|D:2|S:14|M:0|
 *
 * Profiled sequence (query read)         : AGTCCCAAGGATCACACAGA
 * Unprofiled sequence (reference segment): AGTACCAAGGATCACACAGA
 * Aligment                               : AGT CCAAGGATCACACAGA
 * diffs = S:3|S:16|M:0|
 *
 * Profiled sequence (query read)         : AGAAGATCCCCAACAGAGTCCCAGGATCACACAGAC
 * Unprofiled sequence (reference segment): AGAAGATCCCCAACAG---CCC-GGATCACACAGAC
 * Aligment                               : AGAAGATCCCCAACAG---CCC-GGATCACACAGAC
 * diffs = I:16|I:0|I:0|I:3|S:13|M:0|
 */

#ifdef __cplusplus
extern "C"
{
#endif

#ifndef DIFFSTR_H
#define DIFFSTR_H

  //#define diffstr_debug

#include <stdio.h>

  /****************************************************************************
   ******************************* Constants **********************************
   ****************************************************************************/

  enum DIFFSTR_CODES { /**< Diff string codes */
    DIFFCOD_M = 0,  /**< match 'M' */
    DIFFCOD_D = 1,  /**< deletion 'D' (in profiled sequence) */
    DIFFCOD_I = 2,  /**< insertion 'I'(in profiled sequence) */
    DIFFCOD_S = 3,  /**< substitution 'S' */
    DIFFCOD_NUM = 4, /**< number of diffstr codes */
  };

  enum DIFFSTR_CONST {        /**< Diff string constants */
    DIFFSTR_MAXMISMATCH = 61, /**< maximum length of ungapped alignment
			       * for 1 byte/char in diff string */
    DIFFSTR_TYPSHIFT = 6,     /**< Number of bits by which the diff code is 
			       * shifted up in an element of a diff string */
    DIFFSTR_COUNTMASK = 0x3F, /**< Mask out the count in a 
			       * diff string element */
  };

  enum DIFFSTR_OUTPUT_FORMATS { /**< Output formats for compressed alignment string */
    DIFFSTRFORM_RAW = 0,     /**< Raw output format, e.g. "(M:03|I:02|M:10|D:01|M:08|M:00)" */
    DIFFSTRFORM_PLAIN = 1,   /**< Plain raw output format, e.g. "M3I2M10D1M8" */
    DIFFSTRFORM_CIGNORM = 2, /**< GIGAR format */
    DIFFSTRFORM_CIGEXT = 3,  /**< Extended CIGAR format, silent mismatch */
    DIFFSTRFORM_CIGEXT_XMISMATCH = 4,  /**< Extended CIGAR format, with X for mismatch */
  };

  /****************************************************************************
   ********************************* Typedefs *********************************
   ****************************************************************************/
  typedef unsigned char DIFFSTR_T;

  /****************************************************************************
   **************************** Opaque Types **********************************
   ****************************************************************************/

  typedef struct _DiffBuff DiffBuff;
  /**< Contains an explicit alignment representation for each sequence of the pair */

  typedef struct _DiffBlocks DiffBlocks;
  /**< Blocks of ungapped alignment */

  typedef struct _DiffView DiffView;
  /**< Contains formatted alignment string */

  /****************************************************************************
   ************************** Transparent Types *******************************
   ****************************************************************************/

  typedef struct _DiffStr { /**< Container for compressed alignment string */
    DIFFSTR_T *dstrp;/**< Actual alignment string terminated by 0 (M:0) */
    int len;         /**< Length of the string in bytes (including termination) */
    int n_alloc;     /**< Number of bytes allocated */
    int blksz;       /**< Block size for memory allocation */
  } DiffStr;
  
  /****************************************************************************
   ********************************* Macros ***********************************
   ****************************************************************************/
#define DIFFSTR_LENGTH(dfsp) (dfsp)->len
#define DIFFSTR_GET(uc, gap, typ) (typ) = (unsigned char) ((uc) >> DIFFSTR_TYPSHIFT); \
  (gap) = (unsigned char) ((uc) & DIFFSTR_COUNTMASK);
  /* give upper bound for the length of a sring printed in any of the DIFFSTR_OUTPUT_FORMATS
   * this is 5 char per byte plus clipping on either side 2*(ceiling(31*ln2/ln10)+3) */
#define DIFFSTR_MAXPRINTLEN(dfsp) (dfsp->len*5+26)

  /****************************************************************************
   *********************** Methods of Type DiffBlocks *************************
   ****************************************************************************/
  DiffBlocks *diffBlocksCreate(int blksz);
  /**< Constructor.
   * \param blksd Block size for memory allocation
   */

  void diffBlocksDelete(DiffBlocks *p);
  /**< Destructor.
   */

  int diffBlocksGetNumber(const DiffBlocks *p);
  /**< Return the number of blocks
   */

  int diffBlocksGetLen(int *unprof_start, int *prof_start,
		       int blkno, const DiffBlocks *p);
  /**< Return the lenght of an ungapped alignment block number blkno.
   * Return 0 if block not present.
   * \param unprof_start Start coordinate of the block in unprofiled sequence (based to 0).
   * \param prof_start Start coordinate of the block in unprofiled sequence (based to 0).
   * \param Index number of the block.
   * \param p Container of ungapped alignment blocks.
   */


  int diffStrFindBlocks(DiffBlocks *dbp, const DIFFSTR_T *diffstrp);
  /**< Return blocks of ungapped alignment in alignment string.
   * \param dbp Contains block coordinates
   * \param diffstrp Compressed alignement string (terminated by 0)
   */

  /****************************************************************************
   ************************* Methods of Type DiffStr **************************
   ****************************************************************************/
  
  DiffStr *diffStrCreate(int blocksiz);
  /**< Constructor.
   * \param blocksiz Block size for memory allocation. Use default if 0.
   */

  void diffStrDelete(DiffStr *dfsp);
  /**< Destructor.
   */

  int diffStrInit(DiffStr *p, int blocksiz);
  /**< Initialise a compressed alignment string.
   * \param p Compressed alignment string.
   * \param blocksiz Block size for memory allocation. Use default if 0.
   */

  void diffStrCleanUp(DiffStr *p);
  /**< Cleanup contents of container but keep the container itsself
   */

  int diffStrRealloc(DiffStr *dfsp, int n_new);
  /**< Reallocate memory for compressed alignment string.
   * \param dfsp Compressed alignment string.
   * \param n_new Number of bytes to be allocated.
   */

  int diffStrCopy(DiffStr *dfsp, const DIFFSTR_T *diffstrp);
  /**< Copy a compressed alignment string.
   * \param dfsp Container for a compressed alignment string.
   * \param diffstrp Source string.
   */

  int diffStrAdd(DiffStr *top, const DIFFSTR_T *fcp, int len);
  /**< Add to a compressed alignment string.
   * \param dfsp Container for a compressed alignment string.
   * \param fcp Source string.
   * \param len Length of the segement to be copied (in bytes incl termination).
   */

  int diffStrAppend(DiffStr *top, const DiffStr *fromp);
  /**< Concatenate two compressed aligment strings
   * \param top Receiving alignment string.
   * \param fromp String to be appended.
   */

  int diffStrReverse(DiffStr *dfsp, const DIFFSTR_T *diffstrp);
  /**< Return the compressed alignment string for an inversion of the 
   * the aligned sequences.
   * \param dfsp Container receiving rhe reverse compressed alignment string.
   * \param diffstrp Source string.
   */

 int  diffStrCalcSeqLen(int *len_prof, int *len_unprof, const DIFFSTR_T *diffstrp);
  /**< Calculate the lengths of the aligned pair of sequences.
   * \param len_prof Returns the length of the aligned segment of the 
   *        profiled sequence (can be NULL).
   * \param len_unprof Returns the length of the aligned segment of the 
   *        unprofiled sequence (can be NULL).
   * \param diffstrp Compressed alignement string (terminated by 0)
   */

  int diffStrCalcAliLen(int *matchnum, const DIFFSTR_T *diffstrp);
  /**< Calculate the length of the alignment.
   * \return Length of the alignment.
   * \param matchnum Returns the number of exactly matching nucleotides (can be NULL).
   * \param diffstrp Compressed alignement string (terminated by 0)
   */
  
  int diffStrGetDiffStats(int *n_sub, int *n_ins, int *n_del, 
			  const DIFFSTR_T * diffstrp);
  /**< Calculate the number of substitutions, insertions and deletions. 
   * \return error code (!=ERRCODE_SUCCESS on error).
   * \param n_sub Number of substitutions.
   * \param n_ins Number of insertions.
   * \param n_del Number of deletions.
   * \param diffstrp Compressed alignment string.
   */

  int diffStrCalcUngappedBlockLen(int *start_unprof,
				  int *start_prof,
				  const DIFFSTR_T *diffstrp);
  /**< Calculate the length of the next ungapped block in the pairwise alignment.
   * \param start_unprof On input contains minimum position along unprofiled sequence.
   *        On output returns the start of the next block in the unprofiled sequence.
   * \param start_prof On input contains minimum position along profiled sequence.
   *        On output returns the start of the next block in the profiled sequence.
   * \param diffstrp Alignment string.
   */

  int diffStrGenerateFromMismatches(int *dlen, DIFFSTR_T *diffstrp, 
				    const int mmpos[], int mmnum, int qlen);
  /**< Generate a diffstring for a set of ordered mismatches or/and calculate its length.
   * 
   * \return ERRCODE_ASSERT if mismatch positions not in ascending order.
   * \param dlen Returns length of compressed alignment string (incl. termination character). Can be NULL.
   * \param diffstrp Start of the string to be filled (has to have been allocated). Can be NULL.
   * \param mmpos Array of mismatch positions in ascending order
   * \param mmnum Number of mismatches.
   * \param qlen length of query sequence
   */

  int diffStrPrintf(FILE *fp, const DIFFSTR_T *diffstrp, char outform,
		    int clip_start, int clip_end, char is_softclipped);
  /**< Print compressed alignment string in specified output format.
   * \return ERRCODE_ASSERT if unknown output format.
   * \param fp Output stream.
   * \param diffstrp Compressed alignment string.
   * \param outform Output format (one of DIFFSTR_OUTPUT_FORMATS)
   * \param clip_start Start of clipped sequence.
   * \param clip_end End of clipped sequence.
   * \param is_softclipped Flag if != 0 Sequence is 'soft' clipped, i.e. part of the output.
   * \note clip_start, clip_end and is_softclipped are only used with the extended
   * cigar format (outform == DIFFSTRFORM_CIGEXT).
   */

  int diffStrPrintfStr(char *sp, int *nchar, const DIFFSTR_T *diffstrp, 
		       char outform,
		       int clip_start, int clip_end, char is_softclipped);
  /**< Print compressed alignment string in specified output format to a string.
   * \return ERRCODE_ASSERT if unknown output format.
   * \param sp String to which format is written to (must be big enough to hold
   *         the entire resulting string).
   * \param nchar Returns the number of characters written.
   * \param diffstrp Compressed alignment string.
   * \param outform Output format (one of DIFFSTR_OUTPUT_FORMATS)
   * \param clip_start Start of clipped sequence.
   * \param clip_end End of clipped sequence.
   * \param is_softclipped Flag if != 0 Sequence is 'soft' clipped, i.e. part of the output.
   * \note clip_start, clip_end and is_softclipped are only used with the extended
   * cigar format (outform == DIFFSTRFORM_CIGEXT).
   */

  int diffStrParsePlain(DiffStr *dfsp, const char *rawstrp);
  /**< Parse an alignment string written in DIFFSTRFORM_PLAIN format.
   * Stop at end of string or unrecognised character.
   */

  int diffStrParseSimul(DiffStr *dfsp, 
			int *clip_start, int *clip_end, 
			unsigned char isSAMCIGAR, 
			const char *rawstrp);
  /**< Parse an alignment string from a simulated read (see simread.c).
   * Stop at end of string or unrecognised character.
   * \param dfsp Written compressed alignment string.
   * \param clip_start Returns number of clipped bases at beginning of string.
   * \param clip_end Returns the number of clipped bases at the end of the string.
   * \param isSAMCIGAR The CIGAR string is in SAM cigar format with X for mismatches
   * and S or H indicating (soft-/hard-) clipped bases.
   * \param rawstrp CIGAR string as ASCII.
   */

  int diffStrCrop(DIFFSTR_T *dfsp, int *dstrlen, 
		  int start_unprof_target, int end_unprof_target,
		  int *start_unprof, int *end_unprof, 
		  int *start_prof, int *end_prof);
  /**< Crop boundaries of a diff string in place.
   * \param dfsp Compressed alignment string, terminated by 0.
   * \param dstrlen Returns length of resulting alignment string (incl. teminating '\\0')
   * \param start_unprof_target Target start offset on (bases to be cropped) on the 
   *        unprofiled sequence (on forward strand).
   * \param end_unprof_target Target end of the unprofiled sequence (on foward strand).
   * \param start_unprof Number of bases actually cropped at beginning on unprofiled sequence.
   * \param end_unprof Actual end of the cropped segment of unprofiled sequence.
   * \param start_prof Number of bases actually cropped at beginning on profiled sequence.
   * \param end_prof Actual end of the cropped segment of profiled sequence.   
   */
   
  int diffStrFetchSegment(DiffStr *dfsp, const DIFFSTR_T *diffstrp,
			  int start_unprof_target, int end_unprof_target,
			  int *start_unprof, int *end_unprof, 
			  int *start_prof, int *end_prof);
  /**< Fetch segment of a diff string according to a segment of the unprofiled sequence.
   * \param dfsp Container for resulting alignment string.
   * \param diffstrp Original alignment string from which segment is taken.
   * \param start_unprof_target Target start offset on (bases to be cropped) on the 
   *        unprofiled sequence (on forward strand).
   * \param end_unprof_target Target end of the unprofiled sequence (on foward strand).
   * \param start_unprof Number of bases actually cropped at beginning on unprofiled sequence.
   * \param end_unprof Actual end of the cropped segment of unprofiled sequence.
   * \param start_prof Number of bases actually cropped at beginning on profiled sequence.
   * \param end_prof Actual end of the cropped segment of profiled sequence.   
   */

  int diffStrSegment(DiffStr *dfsp, const DIFFSTR_T *diffstrp,
		     int start_unprof_target, int end_unprof_target,
		     int *start_unprof, int *end_unprof, 
		     int *start_prof, int *end_prof);
  /**< Fetch segment of a diff string according to a segment of the unprofiled sequence.
   * \return ERRCODE_NOMATCH if segment does not contain a match. start_unprof, end_unprof,
   *         start_prof, end_prof, dfsp are not defined in that chase.
   * \param dfsp Container for resulting alignment string.
   * \param diffstrp Original alignment string from which segment is taken.
   * \param start_unprof_target Target start offset on the 
   *        unprofiled sequence (on forward strand).
   * \param end_unprof_target Target end of the unprofiled sequence (on foward strand).
   * \param start_unprof Number of bases actually cropped at beginning on unprofiled sequence.
   * \param end_unprof Actual end of the cropped segment of unprofiled sequence.
   * \param start_prof Number of bases actually cropped at beginning on profiled sequence.
   * \param end_prof Actual end of the cropped segment of profiled sequence.   
   */

  int diffStrScore(const DIFFSTR_T *diffstrp, int *swscor, 
		   short match, short mismatch, 
		   short gapopen, short gapextend);
  /**< Calculate the Smith-Waterman score resulting from the diff string 
   */

  int diffStrGetLevenshteinDistance(const DIFFSTR_T *diffstrp);
  /**< Return the Levenshtein or edit distance.
   */

  /******************************************************************************
   ************************** Methods of Type DiffView **************************
   ******************************************************************************/

  DiffView *diffStrCreateView(int blksz);
  /**< Constructor
   * \param Size of the blocks (in bytes) with which memory is allocated.
   */
  
  void diffStrDeleteView(DiffView *p);
  /**< Destructor 
   */

  const char *diffStrGetViewStr(const DiffView *p);
  /**< Returns alignment as character string
   */

  int diffStrAsView(DiffView *dvp, const DIFFSTR_T *dfstrp,
		    char outform,
		    int clip_start, int clip_end, char is_softclipped);
  /**< Generates character string in one of DIFFSTR_OUTPUT_FORMATS (see diffStrPrintfStr).
   * \return ERRCODE_ASSERT if unknown output format.
   * \param dvp Menages memory for string.
   * \param dfstrp Compressed alignment string.
   * \param outform Output format (one of DIFFSTR_OUTPUT_FORMATS)
   * \param clip_start Start of clipped sequence.
   * \param clip_end End of clipped sequence.
   * \param is_softclipped Flag if != 0 Sequence is 'soft' clipped, i.e. part of the output.
   * \note clip_start, clip_end and is_softclipped are only used with the extended
   * cigar format (outform == DIFFSTRFORM_CIGEXT).
   */
  
#endif /* #ifndef DIFFSTR_H */
#ifdef __cplusplus
}
#endif
