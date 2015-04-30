/**< Macros for vectors (arrays)
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

#ifndef VDEF_H
#define VDEF_H

#include <stdlib.h>
#include "elib.h"

  enum VCONST {
    VBLKSZ_DEFAULT = 256,
  };

#define V_TYP(typ) V_##typ
#define V_STRUCT(typ) struct V_##typ##_ {size_t n, n_alloc, n_blksz; typ *vp;}
#define V_INIT(v, blksz) do {			       \
  (v).n_blksz = ((blksz) < 1)? VBLKSZ_DEFAULT:(blksz); \
  (v).n = 0;					       \
  ECALLOCP((v).n_blksz, (v).vp);		       \
  if (NULL == (v).vp) (v).n_alloc = 0;		       \
  else (v).n_alloc = (v).n_blksz;} while(0);
  
#define V_FREE(v) do {free((v).vp);(v).vp = NULL;} while(0);
#define V_DELETE(p) do {if ((p) != NULL) free((p)->vp); free((p)); (p) = NULL;} while(0);

#define V_CREATE(p,blksz) do {						\
    size_t siz = (size_t) ((blksz) < 1)? VBLKSZ_DEFAULT:(blksz);	\
    if (NULL != EMALLOCP0((p)) && NULL != ECALLOCP(siz, (p)->vp)) {	\
      (p)->n_alloc = (p)->n_blksz = siz;				\
      (p)->n = 0;							\
    } else {								\
      V_DELETE((p));							\
    }									\
  } while(0);


#define V_LEN(p) ((p)->n)
#define V_BLKSZ(p) (((size_t) (((p)->n + (p)->n_blksz)/(p)->n_blksz))* \
		    (p)->n_blksz)

#define V_BASPTR(p) (p)->vp
#define V_PTR(p,idx) ((size_t) (idx) < (p)->n)? (p)->vp + (idx): NULL

#define V_NEXTPTR(p) ((p)->n < (p)->n_alloc ||				\
		      (NULL != EREALLOCPP((p)->vp, V_BLKSZ(p)) &&	\
		       ((p)->n_alloc = V_BLKSZ(p)) > (p)->n))?		\
  (p)->vp + (p)->n++: NULL

#define V_APPEND(p,elem) do {			\
    if ((p)->n >= (p)->n_alloc) {		\
      size_t siz = V_BLKSZ(p);			\
      void *hp = EREALLOCP((p)->vp, siz);	\
      if (!hp) return ERRCODE_NOMEM;		\
      (p)->vp = hp; (p)->n_alloc = siz;}	\
    (p)->vp[(p)->n++] = (elem);} while(0)

#define V_DEF_SORTFUNC(typ) int sort_##typ##_quicksort(struct V_##typ##_ *p) {	\
  enum { MINARRSIZE = 7, MAXSTACKSIZE = 60};				\
  int n = (int) p->n;							\
  int i,j, im, il=0, ir=n-1, errcode = ERRCODE_SUCCESS;			\
  int stksz = 0, idxstk[MAXSTACKSIZE];					\
  register typ el;							\
  register typ tmp;							\
  typ *arr = p->vp;							\
  for(;;) {								\
  if (ir-il < MINARRSIZE) {						\
    for (j=il+1;j<=ir;j++) {						\
      el = arr[j];							\
      for (i=j-1; i>=il && arr[i]>el; i--) arr[i+1] = arr[i];		\
      arr[i+1] = el;}							\
    if (!stksz) return errcode;						\
    ir = idxstk[stksz--]; il  = idxstk[stksz--];}			\
  else {								\
    im = (il+ir)>>1;							\
    tmp=arr[im]; arr[im]=arr[il+1]; arr[il+1]=tmp;			\
    if (arr[il] > arr[ir]) {tmp=arr[il]; arr[il]=arr[ir]; arr[ir]=tmp;} \
    if (arr[il+1] > arr[ir]) {tmp=arr[il+1]; arr[il+1]=arr[ir]; arr[ir]=tmp;} \
    if (arr[il] > arr[il+1]) {tmp=arr[il]; arr[il]=arr[il+1]; arr[il+1]=tmp;} \
    i=il+1; j=ir; el = arr[il+1];					\
    for(;;) {								\
      do i++; while (arr[i]<el);					\
      do j--; while (arr[j]>el);					\
      if(j<i) break;							\
      tmp=arr[i]; arr[i]=arr[j]; arr[j]=tmp;}				\
    arr[il+1] = arr[j]; arr[j] = el;					\
    stksz += 2;								\
    if (stksz>MAXSTACKSIZE) return ERRCODE_SORTSTACK;			\
    if (ir-i+1 >= j-il) {						\
      idxstk[stksz] = ir; idxstk[stksz-1] = i; ir = j-1;}		\
    else {								\
      idxstk[stksz] = j-1; idxstk[stksz-1] = il; il = i;}		\
  }}}

#endif
#ifdef __cplusplus
}
#endif
