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

#include <stdio.h>
#include "array.h"

enum {
  ARRAY_BLKSZ_DEFAULT = 4096,
};

/* Array is implemented as a pointer to allocated memory
 * which is offset by 4 size_t types containing
 * [ARRHEADO_NELEM   = 0] the number of elements
 * [ARRHEADO_NALLOC  = 1] the number of elements for which memory is allocated
 * [ARRHEADO_ELEMSIZ = 2] the size (in bytes) of one element
 * [ARRHEADO_BLOCKSIZ= 3] the block size (in bytes) of one allocation block
 */

#ifdef array_debug
static void printfHeader(size_t *headp)
{
  short i;
  for (i= -1*ARRAY_NHEADER;i>0; i--)
    printf("ARRAY[%hi]: %llu\n", i, (long long unsigned int) headp[i]);
}
#endif

void *arrayCreate(size_t elemsiz, size_t blocksiz, const char *nam, int lin)
{
  size_t headsiz = ARRAY_NHEADER * sizeof(size_t);
  size_t blksz, *p;

  if (elemsiz < 1) return NULL;
  blksz = (blocksiz > 0)? 
    ((elemsiz*blocksiz-1)/ARRAY_BLKSZ_DEFAULT+1)*ARRAY_BLKSZ_DEFAULT: 
    ARRAY_BLKSZ_DEFAULT;

  p = (size_t *) ecalloc(1, blksz, nam, lin);
  if (!p) return NULL;
  p += ARRAY_NHEADER;
  p[ARRHEADO_NELEM] = 0;
  p[ARRHEADO_NALLOC] = (blksz - headsiz)/elemsiz;
  p[ARRHEADO_ELEMSIZ] = elemsiz;
  p[ARRHEADO_BLOCKSIZ] = blksz;
  
#ifdef array_debug
  printf("array_debug: arrayCreate() in %s:%i\n", nam, lin);
  printfHeader(p);
#endif

  return (void *) p;
}  

void *arrayRealloc(void *p, size_t new_nelem, char is_final, const char *nam, int lin)
{
  size_t *newp =NULL, *headp = (size_t *) p;
  size_t headsiz = ARRAY_NHEADER * sizeof(size_t);
  size_t newsiz, blocksiz = headp[ARRHEADO_BLOCKSIZ], elsiz = headp[ARRHEADO_ELEMSIZ];

#ifdef array_debug
  printf("array_debug: arrayRealloc(%llu) in %s:%i\n", 
	 (long long unsigned int) new_nelem, nam, lin);
  printfHeader(headp);
#endif

  if (new_nelem > 0) {
    /* allocate space for at least new_nelem */
    if (new_nelem < headp[ARRHEADO_NALLOC]) return p;   
    newsiz = (headsiz+new_nelem*elsiz + blocksiz - 1)/blocksiz;
    newsiz *= blocksiz;
  } else if (is_final) {
    /* free unused memory */
    newsiz = headsiz+headp[ARRHEADO_NELEM]*elsiz;
  } else {
    /* allocate one more block */
    newsiz = (headsiz+(headp[ARRHEADO_NALLOC] + 1)*elsiz - 1)/blocksiz + 1;
    newsiz *= blocksiz;
  }
  newp = (size_t *) erealloc((void *) (headp-ARRAY_NHEADER), newsiz, 0, nam, lin);
  if (!newp)
    return newp;
  newp += ARRAY_NHEADER;
  newp[ARRHEADO_NALLOC] = (newsiz-headsiz)/elsiz;
  
  return (void *) newp;
}
