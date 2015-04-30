/** Generic interface for reading sequencing reads in different formats */

/****************************************************************************
 ****************************************************************************
 *                                                                          *
 *  Copyright (C) 2010 - 2014 Genome Research Ltd.                          *
 *                                                                          *
 *  Author: Hannes Ponstingl (hp3@sanger.ac.uk)                             *
 *                                                                          *
 *  This file is part of SMALT.                                             *
 *                                                                          *
 *  SMALT is free software: you can redistribute it and/or modify it under  *
 *  the terms of the GNU General Public License as published by the Free    *
 *  Software Foundation, either version 3 of the License, or (at your       *
 *  option) any later version.                                              *
 *                                                                          *
 *  This program is distributed in the hope that it will be useful, but     *
 *  WITHOUT ANY WARRANTY; without even the implied warranty of              *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU        *
 *  General Public License for more details.                                *
 *                                                                          *
 *  You should have received a copy of the GNU General Public License along *
 *  with this program.  If not, see <http://www.gnu.org/licenses/>.         *
 *                                                                          *
 ****************************************************************************
 ****************************************************************************/
#ifdef __cplusplus
extern "C"
{
#endif

#ifndef INFMT_H
#define INFMT_H

#include "elib.h"
#include "sequence.h"

  /****************************************************************************
   ******************************* Constants **********************************
   ****************************************************************************/
  
  enum INFMT_FORMATS {
    INFMT_UNKNOWN,
    INFMT_FASTQ,  /**< FASTA/FASTQ format */
    INFMT_SAM,
    INFMT_BAM,
  };
  
  /****************************************************************************
   ***************************** Simple types *********************************
   ****************************************************************************/

  typedef unsigned char INFMT_t;

  /****************************************************************************
   ***************************** Opaque types *********************************
   ****************************************************************************/

  typedef struct _InFmtReader InFmtReader;

  /****************************************************************************
   ********************* Methods of type InFmtReader **************************
   ****************************************************************************/

  InFmtReader *infmtCreateReader(int *errcode,
				 const char *filnamA, const char *filnamB,
#ifdef HAVE_BAMBAMC
				 const char *tmpdir,
#endif
				 const INFMT_t fmt);
  /**< Constructor.
   * \param errcode Returns the error code, ERRCODE_SUCCESS if successful.
   * \param filnamA Input file name (can be "-" for a pipe).
   * \param filnamB Input file name for 2nd mates ("-" for a pipe). 
   *        Only used for read pairs in FASTA or FASTQ input formats
   *       (can be NULL).
   * \param tmpdir Directory name for writing temporary files to (can be NULL).
   *       Only used for SAM/BAM input formats (if NULL, the calling directory
   *       is used).
   * \param Expected input format (one of INFMT_FORMATS). If set to
   *       INFMT_UNKNOWN, an attempt will be made to determine the input
   *       format.
   */

  void infmtDeleteReader(InFmtReader *p);
  /**< Destructor.
   */

  int infmtGetReaderStatus(const InFmtReader *ifrp);
  /**< Returns one of ERRMSG_CODES, e.g. ERRCODE_SUCCESS, ERRCODE_EOF.
   */

  int infmtRead(InFmtReader *ifrp, SeqFastq *sfqAp, SeqFastq *sfqBp, unsigned char *isPair);
  /**< Read a sequencing read or a pair of sequencing reads from input.
   * \param ifrp Reader.
   * \param sfqAp Sequencing read (nucleotides and base call quality).
   * \param sfqBp 2nd mate of a read pair (can be NULL)
   * \param Returns 1 if a read pair was read, 0 otherwise.
   */

  int infmtCheckReads(InFmtReader *ifrp, SeqFastq *sqbufAp, SeqFastq *sqbufBp,
		      SEQNUM_t *seqnum, SEQLEN_t *maxseqlen,  SEQLEN_t *maxnamlen,
		      ErrMsg *errmsgp);
  /**< Check whether the sequences conform to FASTA/FASTQ format for reading.
   *  Return error code or ERRCODE_SUCCESS. Return ERRCODE_SUCCESS if not FASTA/FASTQ.
   *  If FASTA/FASTQ, the file is scrolled back to the start on return.
   * \note This routine wraps seqIOCheckReads for FASTA/FASTQ files.
   *       If the file format is not FASTA/FASTQ ERRCODE_SUCCESS is returned and 
   *       seqnum, maxseqlen and maxnamlen return 0.
   * \param ifrp Reader structure.
   * \param sqbufAp Sequence buffer, gets overwritten and reallocated.
   * \param sqbufBp Sequence buffer for paired reads (NULL for single reads).
   * \param seqnum Returns the number if sequences if FASTA/FASTQ (can be NULL).
   * \param maxseqlen Returns the maximum sequence length if FASTA/FASTQ format (can be NULL).
   * \param maxnamlen Returns the maximum name length if FASTA/FASTQ format (can be NULL).
   * \param errmsgp Returns error messages.
   */

  int infmtReset(InFmtReader *ifrp);
  /**< Rewind the input streams to beginning if reading FASTA/FASTQ 
   * format.  Return one of ERRMSG_CODES. Return ERRCODE_SUCCESS if
   * not FASTA/FASTQ format.
   */
#endif
#ifdef __cplusplus
}
#endif
