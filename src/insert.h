/** Statistics on insert length of paired reads */

/*****************************************************************************
 *****************************************************************************
 *                                                                           *
 *  Copyright (C) 2011 - 2014 Genome Research Ltd.                           * 
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

#ifndef INSERT_H
#define INSERT_H

  //#define insert_debug

#include <stdio.h>
#include <stdint.h>

  /****************************************************************************
   ***************************** Opaque types *********************************
   ****************************************************************************/

  typedef struct _InsSample InsSample;
  typedef struct _InsHist InsHist;

  /******************************************************************************
   ************************* Methods of Type InsSample **************************
   ******************************************************************************/

  InsSample *insCreateSample(int blksz);
  /**< Constructor.
   * \param blksz block size for memory allocations (as the number of elements).
   */

  void insDeleteSample(InsSample *pInsSample);
  /**< Destructor
   */

  void insSetSamplingInterval(InsSample *pInsSample, uint64_t nreads, int nrskip);
  /**< Specify the interval for the sampling of reads.
   * \param pInsSample Type poninter (can be NULL, in which case True is returned).
   * \param nreads Total number of reads from which insert size will be sampled.
   * \param nrskip Target for the number of reads to be skipped. This will be overwritten
   *        so as to fulfill the minimum sample size.
   */

  int insAddSample(InsSample *pInsSample, int32_t insertsiz);
  /**< Add an insert size to the sample.
   * \param pInsSample Sample of insert sizes.
   * \param insertsiz Inserts size to be added.
   */

  int insIsInSample(const InsSample *pInsSample, uint64_t readno);
  /**< Test whether read is in Sample based on the sequential number.
   * \param pInsSample 
   * \param readno Sequential number of the read.
   */

  size_t insGetSample(int32_t **pSample, int *readival, const InsSample *pInsSample);
  /**< Accessor.
   * \return Size of the sample.
   * \param pSample Returns the array of insert sizes (can be NULL)
   * \param readival Returns the sampling interval (can be NULL)
   */

  /****************************************************************************
   ************************* Methods of type InsHist **************************
   ****************************************************************************/ 

  InsHist *insCreateHisto(int32_t iLen);
  /**< Constructor.
   * \param iLen Number of data points in histogram.
   */
  InsHist *insMakeHistoFromSample(const InsSample *pInsSample);
  /**< Construcor from Sample.
   * \param pInsSample Sample of insert sizes.
   */
  void insDeleteHisto(InsHist *p);
  /**< Destructor.
   */
  void insSeedHistoNormal(InsHist *pHist, int32_t mean, int32_t std, int32_t num);
  /**< Seed the histogram from the normal distribution.
   * \param pHist Histogram of insert sizes.
   * \param mean Mean of the normal distribution.
   * \param std Standard deviation of normal distribution.
   * \param num Number of data points seeded.
   */

  void insUpdateHisto(InsHist *pHist, int32_t insiz);
  /**< Update histogram with a data point.
   * \param pHist Histogram.
   * \param insiz New insert size to add.
   */

  int insSmoothHisto(InsHist *pHist);
  /**< Update smoothed histogram
   */

  double insGetHistoProb(int32_t insiz, const InsHist *pHist);
  /**< Return probability for an insert length.
   * \param insiz Insert size.
   * \param pHist Histogram.
   */

  int32_t insGetHistoCount(int32_t *totnum, int32_t insiz, unsigned char is_smooth, 
			   const InsHist *pHist);
  /**< Return counts for an insert size.
   * \param totnum Return the total number of counts in the histogram (can be NULL).
   * \param insiz Insert size.
   * \param is_smooth If != 0 return counts from smoothed histogram.
   * \param pHist Histogram.
   */

  int32_t insGetHistoCountCumulative(int32_t *totnum, int32_t insiz, unsigned char is_smooth, 
				     const InsHist *pHist);
  /**< Return cumulative counts up to and including an insert size.
   * \param totnum Return the total number of counts in the histogram (can be NULL).
   * \param insiz Insert size.
   * \param is_smooth If != 0 return counts from smoothed histogram.
   * \param pHist Histogram.
   */
  uint64_t insGetHistoData(int32_t *lo, int32_t *hi, int32_t *n_bins, const InsHist *pHist);
  /**< Accessor. Returns the number of insert sizes sampled.
   * \param lo Returns lower bound of insert sizes (can be NULL).
   * \param hi Returns upper bound of insert sizes (can be NULL).
   * \param n_bins Returns the number of bins.
   * \param pHist Histogram.
   */
  int32_t insGetHistoQuartiles(int32_t *qlo, int32_t *qhi, const InsHist *pHist);
  /**< Accessor. Returns median of insert sizes.
   * \param qlo Returns lower quartile of insert lengths  (can be NULL).
   * \param qhi Returns upper quartile of insert lengths (can be NULL).
   * \param pHist Histogram.
   */

  int insPrintHisto(FILE *fp, int linwidth, unsigned char is_smooth, const InsHist *pHist);
  /**< Plot the histogram with horizontal bars.
   * \param fp Output stream.
   * \param linwith Line width for plotting.
   * \param is_smooth If != 0 plot smoothed histogram.
   * \param pHist Histogram.
   */

  int insWriteHisto(FILE *fp, unsigned char is_smooth, const InsHist *pHist);
  /**< Write histogram data to ASCII file.
   * \param filnam Name of output file.
   * \param is_smooth If != 0 plot smoothed histogram.
   * \param pHist Histogram.
   */

  InsHist *insReadHisto(int *errcode, const char *filnam);
  /**< Read Histogram written by insWriteHisto() and allocate memory. Return pointer
   * to allocated structure. Returns NULL if errcode != ERRCODE_SUCCESS.
   * \param errcode Returns ERRCODE_SUCCESS or error code.
   * \param filnam Input file name */
#endif
#ifdef __cplusplus
}
#endif
