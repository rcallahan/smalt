/**< Management of sequence intervals (sequence segments).
 */

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

#ifndef INTERVAL_H
#define INTERVAL_H

#include <stdint.h>
#include "sequence.h"

  /****************************************************************************
   ******************************** Types *************************************
   ****************************************************************************/

  typedef uint16_t INTERVALFLG_t; 
  /**< change this with RLSLTFLG_t */

  /****************************************************************************
   ***************************** Opaque types *********************************
   ****************************************************************************/

  typedef struct V_IVAL_ InterVal;

  /******************************************************************************
   ************************** Methods of Type InterVal **************************
   ******************************************************************************/

  InterVal *interValCreate(int32_t blksz);
  /**< Constructor.
   */
  void interValDelete(InterVal *p);
  /**< Destructor.
   */
  void interValBlank(InterVal *p);
  /**< Blank the set of intevals
   */
  int interValAppend(InterVal *p, SEQLEN_t lo, SEQLEN_t hi, SEQNUM_t sx, INTERVALFLG_t flag);
  /**< Add an interval to the set.
   * \param p Set of intervals to be augmented.
   * \param lo start (lower end) of the interval.
   * \param hi end (higher end) of the interval.
   * \param sx Sequence index the inteval refers to.
   * \param flag Bit flag.
   */
  void interValPrune(InterVal *p);
  /**< Merge overlapping intevals
   */
  int interValNum(const InterVal *p);
  /**< Accessor returning the number of intervals (sequence segments)
   */
  int interValGet(SEQLEN_t *lo, SEQLEN_t *hi, SEQNUM_t *sx, INTERVALFLG_t *flg, int idx, 
		  const InterVal *p);
  /**< Acessor. Returns ERRCODE_ARGRANGE if index out of bounds.
   * \param lo returns start (lower end) of the interval.
   * \param hi returns end (higher end) of the interval.
   * \param sx returns Sequence index the inteval refers to.
   * \param flag returns bit flag.
   * \param idx Index of the interval in sequence.
   * \param p Sequence of intervals.
   */

#endif
#ifdef __cplusplus
}
#endif
