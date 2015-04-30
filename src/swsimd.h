/**< Smith-Waterman alignment using SIMD instructions */

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

#ifndef SWSIMD_H
#define SWSIMD_H

#if defined alignment_matrix_debug
#include "sequence.h"
#endif

#include "score.h"
#include "alibuffer.h"
  
  /****************************************************************************
   ************************* Public Alignment Methods *************************
   ****************************************************************************/

  int swSIMDAlignStriped(int *maxscor, 
			 const AliBuffer *abp,
			 const ScoreProfile *profp, 
#ifdef alignment_matrix_debug 
			 const SeqCodec *codecp,
			 const char *profiled_seqp,
#endif
			 const char *unprofiled_seqp,
			 int unprofiled_seqlen);
  /**< Vectorised Smith-Waterman with striped profile.
   * 
   * \param maxscor Returns the maximum Smith-Waterman score.
   * \param abp Buffer for rows of the dynamic programming matrix.
   * \param profp Sequence profile of sequence psqp. Must have been created
   *              with mode = SCORPROF_STRIPED_8 | SCORPROF_STRIPED_16
   * \param unprofiled_seqp Unprofiled sequence in SEQCOD_MANGLED code.
   * \param unprofiled_seqlen Length of the unprofiled sequence.
   */


#endif /* #ifndef SWSIMD_H */
#ifdef __cplusplus
}
#endif
