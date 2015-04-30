/**< Buffers for pairwise alignments */

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

#include <stdlib.h>
#include <limits.h>
#include <config.h>

#include "elib.h"
#include "alibuffer.h"
#include "alibuffer_struct.h"

enum {
  SLACK_BYTES = 32,
  ALIBUFFER_DEFAULT_BLKSZ = 256,
};

/****************************************************************************
 ******************************* Macros *************************************
 ****************************************************************************/

/******************************************************************************
 ********************** Public Methods of Type AliBuffer **********************
 ******************************************************************************/

AliBuffer *aliBufferCreate(int blocksiz)
{
  AliBuffer *p;
  EMALLOCP0(p);
  if (!p) return 0;

  if (blocksiz < 1) blocksiz = ALIBUFFER_DEFAULT_BLKSZ;

  ECALLOCP(blocksiz, p->datap);
  if (!p->datap) {
    aliBufferDelete(p);
    return 0;
  }

  p->allocsiz = (size_t ) blocksiz;
  p->blocksiz = blocksiz;
  p->qlen_max = 0; /* not yet initialised (pointers not set) */
  return p;
}

void aliBufferDelete(AliBuffer *p)
{
  if (p)
    free(p->datap);
  free(p);
}

int aliBufferInit(AliBuffer *p, unsigned int qlen)
{
  unsigned int newlen;
  size_t siz, newsiz;
#ifdef SCORE_SIMD
  int nvec;
  const int nvecelem = 
#ifdef SCORE_SIMD_IMIC
    SCORSIMD_NINTS; 
#else
  SCORSIMD_NSHORTS; 
#endif
#endif 

  if (qlen <= (unsigned int) p->qlen_max)
    return ERRCODE_SUCCESS;

  if (qlen > INT_MAX)
    return ERRCODE_SEQLEN;

  newlen = (qlen + p->blocksiz - 1)/p->blocksiz;
  newlen *= p->blocksiz;
  if (newlen > INT_MAX)
    newlen = INT_MAX;

  siz = (
#ifdef SCORE_SIMD
	 SCORSIMD_MEMALIMASK + 1 + /* slack for alignment to 16/64 byte boundary */
#endif
	 (newlen+1)*sizeof(ALIDPMSCOR_t))*2;

#ifdef SCORE_SIMD
  nvec = (newlen + nvecelem - 1)/nvecelem;
  newsiz = (1 + /* slack for alignment to 16/64 byte boundary */
	    3 * nvec) * sizeof(SIMDV_t);
  if (siz < newsiz)
    siz = newsiz;
#endif  

  if (siz > p->allocsiz) {
    void *hp;
    newsiz = (siz + p->blocksiz - 1)/p->blocksiz;
    newsiz *= p->blocksiz;

    hp = EREALLOCP(p->datap, newsiz);
    if (!hp)
      return ERRCODE_NOMEM;
    p->datap = hp;
    p->allocsiz = newsiz;
  }
  p->baseHp = (ALIDPMSCOR_t *) SCORE_ALIGN_MEMORY(p->datap);
  p->baseEp = p->baseHp + newlen + 1;
  p->baseEp = (ALIDPMSCOR_t *) SCORE_ALIGN_MEMORY(p->baseEp);
  

#ifdef SCORE_SIMD
  /* align to 16/64 byte boundary for SIMD register */
  p->H1v = (SIMDV_t *) SCORE_ALIGN_MEMORY(p->datap);
  p->H2v = p->H1v + nvec;
  p->Ev  = p->H2v + nvec;
#endif

  p->qlen_max = (int) newlen;

  return ERRCODE_SUCCESS;
}
