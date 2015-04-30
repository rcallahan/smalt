/** Base quality statistics */

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

#ifndef BASQUAL_H
#define BASQUAL_H

  /* #define basqual_debug */

#include "sequence.h"

  /****************************************************************************
   ***************************** Opaque types *********************************
   ****************************************************************************/

  typedef struct _BasQualFreq_t BasQualFreq_t;
    
  /****************************************************************************
   ********************** Methods of type BasQualFreq *************************
   ****************************************************************************/

  BasQualFreq_t *basQualFreqCreate(uint8_t qmin, uint8_t nq, uint32_t readlen);
  /**< Constructor.
   * \param qmin Minimum base quality (based to SEQCOD_QVAL_OFFS)
   * \param nq number of quality values, i.e. range of base qualities
   *        is {qmin, ..., qmin + nq - 1}
   * \param readlen Read length (maximum)
   */
  void basQualFreqDelete(BasQualFreq_t *p);
  /**< Destructor
   */
  uint32_t *basQualFreqFld(int *errcode, uint32_t r,  uint8_t i, uint8_t j, BasQualFreq_t *bqfp);
  /**< Return a pointer to the field of base quality counts.
   * \param errcode Returns error code (one of ERRMSG_CODES)
   * \param r Position in sequence (based to 0).
   * \param i Quality value at this position r.
   * \param j Quality value of previous position (r-1).
   * \param bqfp Pointer to type instance.
   */
  uint32_t basQualFreqGetData(const BasQualFreq_t *p, const uint32_t **q0p, uint8_t *nq, uint8_t *qmin, 
			      const uint32_t **qtp, const uint64_t **qsp);
  /**< Accessor. 
   * \return (maximum) read length.
   * \param p Pointer to type.
   * \param q0p Returns pointer to array of counts for first quality value.
   *           Array size is *nq.
   * \param nq Returns the maximum quality value
   * \param qmin Returns the minimum quality value (based on SEQCOD_QVAL_OFFS)
   * \param qtp Returns transition counts qtp[r][i][j] (qtp[readlen-1][nq][nq])
   *        for the transition from quality value i to quality value j when going from position r to
   *        r+1 in the sequence. 
   * \param qsp Returns marginals of qtp i.e. qsp[r][i] = Sum_j qtp[r][i][j].
   */ 
  int basQualFreqWrite(const BasQualFreq_t *p, const char *filnam);
  /**< Write base quality stats to file.
   * \param p Structure containing the base quality stats.
   * \param filnam Name of output file (extension will be added).
   */
  BasQualFreq_t *basQualFreqRead(int *errcode, const char *filnam);
  /**< Read base quality stats from file and allocate memory.
   * \param errcode Returns one of ERRMSG_CODES.
   * \param filnam Name of output file (w/o extension).
   */

  int basQualFreqFromFastq(BasQualFreq_t *bqfp, SeqFastq *sqbufp, SeqIO *sifp);
  /**< Obtain frequencies from FASTQ file
   * \param bqfp Base quality frequencies.
   * \param sqbufp Read buffer.
   * \param sifp FASTQ input file.
   */

  void basQualFreqPrint(FILE *fp, const BasQualFreq_t *bqfp);
  /**< Print statistics as text.
   */

  int basQualFreqSum(BasQualFreq_t *bqfp);
  /**< Update the sums of counts.
   */

  int basQualFreqSimulate(char *qualp, uint32_t len, const BasQualFreq_t *bqfp);
  /**< Simulate a string of base qualities.
   * \param qualp String of base qualities (has to have been allocated to length len+1).
   * \param len Number of base quality characters (excl. teminating 0)
   * \param bqfp Base quality count statistics (empirial distribution).
   */

  /****************************************************************************
   **************************** Other Methods *********************************
   ****************************************************************************/

  int basQualFindExtrema(SeqIO *sifp, SeqFastq *sqbufp, uint64_t *nreads, 
			 uint32_t *maxlen, uint32_t *minlen,
			 uint8_t *maxq, uint8_t *minq);
  /**< Find extreme base quality values and reads lengths in a FASTQ file.
   * \param sifp FASTQ file.
   * \param sqbufp Read buffer.
   * \param nreads Returns number of reads (can be NULL).
   * \param maxlen Returns maximum read length (can be NULL).
   * \param minlen Returns minimum read length (can be NULL).
   * \param maxq   Returns maximum base quality (can be NULL).
   * \param minq   Returns minimum base quality (can be NULL).
   */

#endif
#ifdef __cplusplus
}
#endif
