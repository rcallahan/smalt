/** Sort routines */

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

#ifndef SORT_H
#define SORT_H

  //#define sort_debug

#include<stdint.h>

  int sort2UINTarraysByQuickSort(int n, uint32_t arr[], uint32_t brr[]);
  /**< sort two uint32_t arrays according to the values of the first
   * in ascending order. Return ERRCODE_SUCCESS on success, ERRCODE_SORTSTACK
   * on stack overflow.
   */
  int sortUINT32arrayByQuickSort(int n, uint32_t arr[]);
  /**< sort a single uint32_t array using quick sort 
   */
  int sortUINT64arrayByQuickSort(int n, uint64_t arr[]);

  int sortUINT64andUINT32ArraysByQuickSort(int n, uint64_t arr[], uint32_t brr[]);

  void sortMultiKey(const char *x[], int num, int maxdepth);
  /**< Sort a set of character strings using the multi key sort algorithm
   * by Bentley & Sedgewick (1997).
   *
   * \param x array of strings to be sorted. Must be terminated by '\\0'.
   * \param num Number of strings in the array.
   * \param maxdepth Stop sorting at a maximum prefix length of \<maxdepth\>.
   */

  void sortSufficesByMultiKeyQuickSort(const char *hstrp, uint32_t *sfxidxp, 
				       uint32_t nsfx, short startdepth, short maxdepth);
  /**< Sort a set of suffices up to a specified depth using the multi key sort algorithm
   * by Bentley & Sedgewick (1997).
   *
   * \param hstrp Character string from which the suffices are taken.
   * \param sfxidxp Array of suffix starting positions in the character string hstrp.
   * \param nsfx Number of suffixes, i.e. size of array sfxidxp.
   * \param startdepth > 0 if suffices are already sorted up to depth \<startdepth\>/
   * \param maxdepth Stop sorting at a maximum prefix length of \<maxdepth\>.
   *
   * \note The the maximum offset in array sfxidxp must be not greater than the number
   * of characters in hstrp minus the maximum depth.
   */

#endif
#ifdef __cplusplus
}
#endif
