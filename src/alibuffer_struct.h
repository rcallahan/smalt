/**< Alignment buffer structure used by various alignment methods */

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

#ifndef ALIBUFFER_STRUCT_H
#define ALIBUFFER_STRUCT_H

#include <config.h>
#include "score.h"

  typedef int ALIDPMSCOR_t;
  
  struct _AliBuffer { /** Row of the dynamic programming matrix */
    int qlen_max;          /**< Max query length for which buffer is currently initialised */
    unsigned char *datap;  /**< point of memory allocation */
    ALIDPMSCOR_t *baseHp;  /**< buffer for maximum cell scores H[i,j] along
			    * DPM row */
    ALIDPMSCOR_t *baseEp;  /**< buffer for maximum gap scores E[i-1, j] -
			    * [gap_ext|gap_init] along DPM row */
    size_t allocsiz;       /**< number of ALIDPMSCOR_t elements currently 
			    * allocated per buffer */
    int blocksiz;          /**< Block size (granularity) for memory allocation
			    * as the number of ALIDPMSCOR_t elements per buffer. */
#ifdef SCORE_SIMD
    SIMDV_t *H1v;
    SIMDV_t *H2v;
    SIMDV_t *Ev;    
#endif
  };

#endif /* ifndef ALIBUFFER_STRUCT_H */
#ifdef __cplusplus
}
#endif
