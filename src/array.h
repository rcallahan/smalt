/** Handling of arrays */

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

#ifndef ARRAY_H
#define ARRAY_H

  //#define array_debug

#include <stdlib.h>
#include "elib.h"

  enum ARRAY_CONST {
    ARRHEADO_NELEM = -1,   /** number of elements (offset from array pointer */
    ARRHEADO_NALLOC =-2,   /** number of elements for which memory is allocated */
    ARRHEADO_ELEMSIZ = -3, /** size of 1 element in bytes */
    ARRHEADO_BLOCKSIZ = -4,/** block size (in bytes) for allocation */
    ARRAY_NHEADER = 4      /**< Number of size_t types in header */
  };

  void *arrayCreate(size_t elemsiz, size_t blocksiz, const char *nam, int lin);

  void *arrayRealloc(void *p, size_t new_nelem, char is_final, const char *nam, int lin);
  /**< Re-allocate memory 
   * if nelem > 0 : for at least nelem elements
   * if nelem == 0: allocate 1 more block
   * if nelem < 1 : free unused memory
   */

  int arrayCast(void **p, size_t elemsiz, const char *nam, int lin);
  /**< Cast array to a new element size */

#define ARRCREATE(p, blocksiz) ((p) = arrayCreate(sizeof(*(p)),(blocksiz),__FILE__, __LINE__))

#define ARRDELETE(p) free(((size_t *) (p)) - ARRAY_NHEADER); (p) = NULL;

#define ARRLEN(p) ((size_t *) ((p)))[ARRHEADO_NELEM]

#define ARRNALLOC(p) ((size_t *) ((p)))[ARRHEADO_NALLOC]

#define ARREALLOC(p, nelem) arrayRealloc((p), (size_t) (nelem), 0, __FILE__, __LINE__)
    
#define ARRNEXTP(elemp, arrp) \
    if (ARRLEN((arrp)) >= ARRNALLOC((arrp))) {\
      (elemp) = ARREALLOC((arrp), 0);\
      if ((elemp) != NULL) {\
        (arrp) = (elemp);\
        (elemp) = (arrp) + ARRLEN((arrp))++;\
      }\
    } else (elemp) = (arrp) + ARRLEN((arrp))++;
  /**< Increment current element of array and allocate memory if necessary, return error code */

#define ARRCURR(p) (p)[ARRLEN((p))]

#define ARRCURP(p) ((p)+ARRLEN((p)))


#endif

#ifdef __cplusplus
}
#endif
