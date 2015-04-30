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
#include <limits.h>

#include "vdef.h"
#include "interval.h"

/******************************************************************************
 ********************************* Constants **********************************
 ******************************************************************************/

enum INTERVAL_CONST {
  IVAL_BLKSZ_DEFAULT = 256, /**< Default block size for memory allocation */
};

/******************************************************************************
 ******************************* Private Types ********************************
 ******************************************************************************/

typedef struct _IVAL {
  SEQLEN_t lo;           /**< lower boundary (incl) */
  SEQLEN_t hi;           /**< upper boundary (incl) */
  INTERVALFLG_t flg;     /**< combination of RESULT_STATUS_FLAGS */
  SEQNUM_t sx;           /**< sequence index */
} IVAL;

/******************************************************************************
 ******************************* Opaque Types *********************************
 ******************************************************************************/

V_STRUCT(IVAL); /**< Sequence of intervals */

/******************************************************************************
 ************************ Methods of private Type IVAL ************************
 ******************************************************************************/

static int cmpIVAL(const void *a, const void *b)
{
  const IVAL *ap = (IVAL *) a;
  const IVAL *bp = (IVAL *) b;
  if (ap->sx < bp->sx) return -1;
  if (ap->sx > bp->sx) return 1;
  if (ap->lo < bp->lo) return -1;
  if (ap->lo > bp->lo)return 1;
  if (ap->hi < bp->hi) return -1;
  if (ap->hi > bp->hi) return 1;
  return 0;
}

/******************************************************************************
 *********************** Public Methods of Type InterVal **********************
 ******************************************************************************/

InterVal *interValCreate(int32_t blksz)
{
  InterVal *p;
  if (blksz < 1)
    blksz = IVAL_BLKSZ_DEFAULT;
  V_CREATE(p, blksz);
  return p;
}

void interValDelete(InterVal *p)
{
  V_DELETE(p);
}

void interValBlank(InterVal *p)
{
  V_LEN(p) = 0;
}

int interValAppend(InterVal *p, SEQLEN_t lo, SEQLEN_t hi, 
		   SEQNUM_t sx, INTERVALFLG_t flag)
{
  IVAL *ivp;
  if (hi < lo)
    return ERRCODE_ASSERT;

  ivp = V_NEXTPTR(p);

  if (NULL == ivp)
    return ERRCODE_NOMEM;

  ivp->hi = hi;
  ivp->lo = lo;
  ivp->sx = sx;
  ivp->flg = flag;

  if (V_LEN(p) > INT_MAX)
    return ERRCODE_OVERFLOW;
  
  return ERRCODE_SUCCESS;
}

void interValPrune(InterVal *p)
{
  int i, j, niv = (int) V_LEN(p);
  IVAL *ivp = V_BASPTR(p);

  if (niv > 0) {
    qsort(ivp, niv, sizeof(IVAL), cmpIVAL);

    for (i=0, j=1; j<niv; j++) {
      if (ivp[j].sx ==  ivp[i].sx && 
	  ivp[j].lo <= ivp[i].hi) {
	/* overlapping intervals */
	if (ivp[j].hi > ivp[i].hi) 
	  ivp[i].hi = ivp[j].hi;
      } else {
	if (++i < j) 
	  ivp[i] = ivp[j];
      }   
    }

    V_LEN(p) = i+1;
  }
}

int interValNum(const InterVal *p)
{
  return (int) V_LEN(p);
}

int interValGet(SEQLEN_t *lo, SEQLEN_t *hi, SEQNUM_t *sx, INTERVALFLG_t *flg, 
		int idx, const InterVal *p)
{
  const IVAL *ivp = V_PTR(p, idx);
  if (!ivp)
    return ERRCODE_ARGRANGE;
  if (ivp->lo > ivp->hi)
    return ERRCODE_ASSERT;
  if (lo) *lo = ivp->lo;
  if (hi) *hi = ivp->hi;
  if (sx) *sx = ivp->sx;
  if (flg) *flg = ivp->flg;

  return ERRCODE_SUCCESS;
}
