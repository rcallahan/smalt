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

#define NDEBUG
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "elib.h"
#include "sort.h"

/*****************************************************************************
 **************************** Private Methods ********************************
 *****************************************************************************/
#define SWAPX(p,q) tmp=x[(p)]; x[(p)]=x[(q)]; x[(q)]=tmp
#define INT2CHARX(q) (int) x[(q)][depth]
#define VSWAPX(p,q) for (i=(p), j=(q), k=r; k-->0; i++, j++) {SWAPX(i, j);}
#define MIN(p,q) ((p) < (q))? (p) : (q)

static void mksort(
#ifdef sort_debug
		   int *level,
#endif
		   const char *x[], 
		   int n, int depth, int maxdepth)
     /**< sorts n strings x[0], x[1], ..., x[n-1] by their prefix of length depth+1 if
      * they are known to have been sorted by their prefix of length depth.
      * first call is mksort(x, n, 0, n).
      *
      * Bentley & Sedgewick (1997) Fast Algorithms for sorting and searching strings
      * in Proc. 8th Annual Symposium on Discrete Algorithms, ACM.
      */
{
  int a,b,c,d,i,j,k,r,v;
  const char *tmp;

  if (n <= 1) return;

#ifdef sort_debug
  printf("mksort: level = %i\n", (*level)++);
#endif
  a = rand() % n;
  SWAPX(0,a);
  v = INT2CHARX(0);
  a = b = 1;
  c = d = n-1;
  
  for(;;) {
    while (b <= c && (r = INT2CHARX(b) - v) <= 0 ) {
      if (r == 0) {
	SWAPX(a,b);
	a++;
      }
      b++;
    }
    while ( b <= c && (r = INT2CHARX(c) - v) >= 0) {
      if (r == 0) {
	SWAPX(c, d);
	d--;
      }
      c--;
    }
    if (b > c) break;
    SWAPX(b, c);
    b++;
    c--;
  }

  r = MIN(a, b-a);
  //vswap(0, b-r, r, x);
  VSWAPX(0, b-r);

  r = MIN(d-c, n-d-1);
  //vswap(0, n-r, r, x);
  VSWAPX(b, n-r);  

  r = b-a;
  mksort(
#ifdef sort_debug
	 level,
#endif
	 x, r, depth, maxdepth
	 );

  if (INT2CHARX(r) != 0 && depth < maxdepth) /* termination signal */
    mksort(
#ifdef sort_debug
	 level,
#endif
	 x+r, a+n-d-1, depth+1, maxdepth
	 );
  
  r = d-c;
  mksort(
#ifdef sort_debug
	 level,
#endif
	 x+n-r, r, depth, maxdepth
	 );
}

#define SFX_SWAPX(p,q) tmp=sfxidxp[(p)]; sfxidxp[(p)]=sfxidxp[(q)]; sfxidxp[(q)]=tmp
#define SFX_INT2CHARX(q) (short) hstrp[sfxidxp[(q)]+depth]
#define SFX_VSWAPX(p,q,r) for (ip=sfxidxp+(p), jp=sfxidxp+(q), k=(r); k-->0; ip++, jp++) \
       {tmp = *ip; *ip=*jp; *jp=tmp;}			     
/**< swap vectors x[p], x[p+1], .... x[p+r-1] with x[q], x[q+1], ..., x[q+r-1] */

static void mkeyQSortSuffix(
#ifdef sort_debug
		   int *level,
#endif
		   const char *hstrp,
		   uint32_t *sfxidxp,
		   uint32_t nsfx, short depth, short maxdepth)
     /**< sorts n-maxdepth suffices of a string hstrp up to a suffix length 
      * of maxdepth caracters. The suffices must have been sorted by their prefix of length 
      * depth prior to calling this routine.
      * For a recursive sort of suffices up to a lenght maxdepth first call with depth = 0.
      *
      * \param hstrp Host string must be at least max[sfxidxp[]] + maxdepth characters long.
      * \param sfxidxp Array of nsfx start indices of the suffices in hstrp
      * \param nsfx size of the array sfxidxp
      * \param depth Current sort depth before the routine is called.
      * \param maximum depth of the sort.
      *
      * \Note Bentley & Sedgewick (1997) Fast Algorithms for sorting and searching strings
      * in Proc. 8th Annual Symposium on Discrete Algorithms, ACM.
      */
{
  short r,v;
  uint32_t a,b,c,d,m;
  uint32_t *ip, *jp, k;        /* variables used in macros */
  register uint32_t tmp; /* used in macro */

  if (nsfx <= 1) return;

#ifdef sort_debug
  printf("mkeyQSortSuffix: level = %i\n", (*level)++);
#endif
  a = rand() % nsfx;
  SFX_SWAPX(0,a);
  v = SFX_INT2CHARX(0);
  a = b = 1;
  c = d = nsfx-1;
  
  for(;;) {
    while (b <= c && (r = SFX_INT2CHARX(b) - v) <= 0 ) {
      if (!r) {
	SFX_SWAPX(a,b);
	a++;
      }
      b++;
    }
    while ( b <= c && (r = SFX_INT2CHARX(c) - v) >= 0) {
      if (!r) {
	SFX_SWAPX(c, d);
	d--;
      }
      c--;
    }
    if (b > c) break;
    SFX_SWAPX(b, c);
    b++;
    c--;
  }

  m = MIN(a, b-a);
  //vswap(0, b-r, r, x);
  SFX_VSWAPX(0, b-m, m);

  m = MIN(d-c, nsfx-d-1);
  //vswap(0, n-r, r, x);
  SFX_VSWAPX(b, nsfx-m, m);  

  m = b-a;
  mkeyQSortSuffix(
#ifdef sort_debug
	 level,
#endif
	 hstrp, sfxidxp, m, depth, maxdepth
	 );

  if (SFX_INT2CHARX(m) != 0 && depth < maxdepth) /* termination signal */
    mkeyQSortSuffix (
#ifdef sort_debug
	 level,
#endif
	 hstrp, sfxidxp+m, a+nsfx-d-1, depth+1, maxdepth
	 );
  
  m = d-c;
  mkeyQSortSuffix(
#ifdef sort_debug
	 level,
#endif
	 hstrp, sfxidxp+nsfx-m, m, depth, maxdepth
	 );
}

/*****************************************************************************
 ***************************** Public Methods ********************************
 *****************************************************************************/

#define MAXSTACKSIZE 60
#define MINARRSIZE 7   /* subarrays smaller than this value are sorted by straight insertion */

#define SWAP(a,b) temp=(a);(a)=b;(b)=temp;

int sort2UINTarraysByQuickSort(int n, uint32_t arr[], uint32_t brr[])
     /* sort two uint32_t arrays according to the values of the first
      * in ascending order. Return ERRCODE_SUCCESS on success, ERRCODE_SORTSTACK
      * on stack overflow.
      */

{
  int i,j; /* partitioning indices */
  int i_left=0, i_middle, i_right=n-1; /* left, middle, right element of current subarray */
  register uint32_t partelem_a;  /* partitioning element */
  uint32_t partelem_b;
  int index_stack[MAXSTACKSIZE]; /* stack of index intervals for each subarray */
  int stack_size = 0;
  int errcode = ERRCODE_SUCCESS;
  register uint32_t temp;

  for(;;) {
    if (i_right-i_left < MINARRSIZE) {/* straight insertion */
      for (j=i_left+1;j<=i_right;j++) {
	partelem_a = arr[j];
	partelem_b = brr[j];
	for (i=j-1; i>=i_left && arr[i]>partelem_a; i--) {
	  arr[i+1] = arr[i];
	  brr[i+1] = brr[i];
	}
	arr[i+1] = partelem_a;
	brr[i+1] = partelem_b;
      }
      if (!stack_size) return errcode;
      /* pop index interval for next subarray from stack */
      i_right = index_stack[stack_size--];
      i_left  = index_stack[stack_size--];
    } 
    else {
      /* qicksort, partitioning exchange
       *
       * choose as partitioning element the median of left, middle and
       * right elements of current subarray and arrange
       * arr[i_left]<=arr[i_left+1]<=arr[i_right] */
      i_middle = (i_left+i_right)>>1;
/*       if (i_middle < 0 || i_middle >= n) { */
/* 	printf("SORT.C: i_middle = %i, i_left = %i, i_right = %i, n = %i\n",  */
/* 	       i_middle, i_left, i_right, n); */
/* 	return ERRCODE_FAILURE; */
/*       } */
      SWAP(arr[i_middle], arr[i_left+1]);
      SWAP(brr[i_middle], brr[i_left+1]);
      if (arr[i_left] > arr[i_right]) {
	SWAP(arr[i_left], arr[i_right]);
	SWAP(brr[i_left], brr[i_right]);
      }
      if (arr[i_left+1] > arr[i_right]) {
	SWAP(arr[i_left+1], arr[i_right]);
	SWAP(brr[i_left+1], brr[i_right]);
      }
      if (arr[i_left] > arr[i_left+1]) {
	SWAP(arr[i_left], arr[i_left+1]);
	SWAP(brr[i_left], brr[i_left+1]);
      }

      /* initialise partitioning pointers */
      i=i_left+1;
      j=i_right;
      
      /* pick partitioning element */
      partelem_a = arr[i_left+1];
      partelem_b = brr[i_left+1];

      /* partitioning exchange loop */
      for(;;) {
      	do i++; while (arr[i]<partelem_a);
      	do j--; while (arr[j]>partelem_a);
      	if(j<i) break; /* partitioning pointers crossed */
      	SWAP(arr[i],arr[j]); /* exchange */
      	SWAP(brr[i],brr[j])
      }
      /* insert partitioning element */
      arr[i_left+1] = arr[j];
      brr[i_left+1] = brr[j];
      arr[j] = partelem_a;
      brr[j] = partelem_b;

      /* push the larger of the two partitions on stack */
      stack_size += 2;
      if(stack_size>MAXSTACKSIZE) return ERRCODE_SORTSTACK; /* stack overflow */

      if (i_right-i+1 >= j-i_left) {
	index_stack[stack_size] = i_right;
	index_stack[stack_size-1] = i;
	i_right = j-1;
      } else {
	index_stack[stack_size] = j-1;
	index_stack[stack_size-1] = i_left;
	i_left = i;
      }
    } /* if (i_right-i_left < MINARRSIZE) else */
  } /* loop over all subarrays */
}


int sortUINT32arrayByQuickSort(int n, uint32_t arr[])
     /* sort two uint32_t arrays according to the values of the first
      * in ascending order. Return ERRCODE_SUCCESS on success, ERRCODE_SORTSTACK
      * on stack overflow.
      */

{
  int i,j; /* partitioning indices */
  int i_left=0, i_middle, i_right=n-1; /* left, middle, right element of current subarray */
  register uint32_t partelem;  /* partitioning element */
  int index_stack[MAXSTACKSIZE]; /* stack of index intervals for each subarray */
  int stack_size = 0;
  int errcode = ERRCODE_SUCCESS;
  register uint32_t temp;

  for(;;) {
    if (i_right-i_left < MINARRSIZE) {/* straight insertion */
      for (j=i_left+1;j<=i_right;j++) {
	partelem = arr[j];
	for (i=j-1; i>=i_left && arr[i]>partelem; i--) {
	  arr[i+1] = arr[i];
	}
	arr[i+1] = partelem;
      }
      if (!stack_size) return errcode;
      /* pop index interval for next subarray from stack */
      i_right = index_stack[stack_size--];
      i_left  = index_stack[stack_size--];
    } 
    else {
      /* qicksort, partitioning exchange
       *
       * choose as partitioning element the median of left, middle and
       * right elements of current subarray and arrange
       * arr[i_left]<=arr[i_left+1]<=arr[i_right] */
      i_middle = (i_left+i_right)>>1;
      SWAP(arr[i_middle], arr[i_left+1]);
      if (arr[i_left] > arr[i_right]) {
	SWAP(arr[i_left], arr[i_right]);
      }
      if (arr[i_left+1] > arr[i_right]) {
	SWAP(arr[i_left+1], arr[i_right]);
      }
      if (arr[i_left] > arr[i_left+1]) {
	SWAP(arr[i_left], arr[i_left+1]);
      }

      /* initialise partitioning pointers */
      i=i_left+1;
      j=i_right;
      
      /* pick partitioning element */
      partelem = arr[i_left+1];

      /* partitioning exchange loop */
      for(;;) {
      	do i++; while (arr[i]<partelem);
      	do j--; while (arr[j]>partelem);
      	if(j<i) break; /* partitioning pointers crossed */
      	SWAP(arr[i],arr[j]); /* exchange */
      }
      /* insert partitioning element */
      arr[i_left+1] = arr[j];
      arr[j] = partelem;

      /* push the larger of the two partitions on stack */
      stack_size += 2;
      if(stack_size>MAXSTACKSIZE) return ERRCODE_SORTSTACK; /* stack overflow */

      if (i_right-i+1 >= j-i_left) {
	index_stack[stack_size] = i_right;
	index_stack[stack_size-1] = i;
	i_right = j-1;
      } else {
	index_stack[stack_size] = j-1;
	index_stack[stack_size-1] = i_left;
	i_left = i;
      }
    } /* if (i_right-i_left < MINARRSIZE) else */
  } /* loop over all subarrays */
}

int sortUINT64arrayByQuickSort(int n, uint64_t arr[])
     /* sort two uint32_t arrays according to the values of the first
      * in ascending order. Return ERRCODE_SUCCESS on success, ERRCODE_SORTSTACK
      * on stack overflow.
      */

{
  int i,j; /* partitioning indices */
  int i_left=0, i_middle, i_right=n-1; /* left, middle, right element of current subarray */
  register uint64_t partelem;  /* partitioning element */
  int index_stack[MAXSTACKSIZE]; /* stack of index intervals for each subarray */
  int stack_size = 0;
  int errcode = ERRCODE_SUCCESS;
  register uint64_t temp;

  for(;;) {
    if (i_right-i_left < MINARRSIZE) {/* straight insertion */
      for (j=i_left+1;j<=i_right;j++) {
	partelem = arr[j];
	for (i=j-1; i>=i_left && arr[i]>partelem; i--) {
	  arr[i+1] = arr[i];
	}
	arr[i+1] = partelem;
      }
      if (!stack_size) return errcode;
      /* pop index interval for next subarray from stack */
      i_right = index_stack[stack_size--];
      i_left  = index_stack[stack_size--];
    } 
    else {
      /* qicksort, partitioning exchange
       *
       * choose as partitioning element the median of left, middle and
       * right elements of current subarray and arrange
       * arr[i_left]<=arr[i_left+1]<=arr[i_right] */
      i_middle = (i_left+i_right)>>1;
      SWAP(arr[i_middle], arr[i_left+1]);
      if (arr[i_left] > arr[i_right]) {
	SWAP(arr[i_left], arr[i_right]);
      }
      if (arr[i_left+1] > arr[i_right]) {
	SWAP(arr[i_left+1], arr[i_right]);
      }
      if (arr[i_left] > arr[i_left+1]) {
	SWAP(arr[i_left], arr[i_left+1]);
      }

      /* initialise partitioning pointers */
      i=i_left+1;
      j=i_right;
      
      /* pick partitioning element */
      partelem = arr[i_left+1];

      /* partitioning exchange loop */
      for(;;) {
      	do i++; while (arr[i]<partelem);
      	do j--; while (arr[j]>partelem);
      	if(j<i) break; /* partitioning pointers crossed */
      	SWAP(arr[i],arr[j]); /* exchange */
      }
      /* insert partitioning element */
      arr[i_left+1] = arr[j];
      arr[j] = partelem;

      /* push the larger of the two partitions on stack */
      stack_size += 2;
      if(stack_size>MAXSTACKSIZE) return ERRCODE_SORTSTACK; /* stack overflow */

      if (i_right-i+1 >= j-i_left) {
	index_stack[stack_size] = i_right;
	index_stack[stack_size-1] = i;
	i_right = j-1;
      } else {
	index_stack[stack_size] = j-1;
	index_stack[stack_size-1] = i_left;
	i_left = i;
      }
    } /* if (i_right-i_left < MINARRSIZE) else */
  } /* loop over all subarrays */
}

#define SWAP2(a,b) temp2=(a);(a)=b;(b)=temp2;
int sortUINT64andUINT32ArraysByQuickSort(int n, uint64_t arr[], uint32_t brr[])
     /* sort two uint32_t arrays according to the values of the first
      * in ascending order. Return ERRCODE_SUCCESS on success, ERRCODE_SORTSTACK
      * on stack overflow.
      */

{
  int i,j; /* partitioning indices */
  int i_left=0, i_middle, i_right=n-1; /* left, middle, right element of current subarray */
  register uint64_t partelem_a, temp;  /* partitioning element */
  uint32_t partelem_b;
  register uint32_t temp2;
  int index_stack[MAXSTACKSIZE]; /* stack of index intervals for each subarray */
  int stack_size = 0;
  int errcode = ERRCODE_SUCCESS;

  

  for(;;) {
    if (i_right-i_left < MINARRSIZE) {/* straight insertion */
      for (j=i_left+1;j<=i_right;j++) {
	partelem_a = arr[j];
	partelem_b = brr[j];
	for (i=j-1; i>=i_left && arr[i]>partelem_a; i--) {
	  arr[i+1] = arr[i];
	  brr[i+1] = brr[i];
	}
	arr[i+1] = partelem_a;
	brr[i+1] = partelem_b;
      }
      if (!stack_size) return errcode;
      /* pop index interval for next subarray from stack */
      i_right = index_stack[stack_size--];
      i_left  = index_stack[stack_size--];
    } 
    else {
      /* qicksort, partitioning exchange
       *
       * choose as partitioning element the median of left, middle and
       * right elements of current subarray and arrange
       * arr[i_left]<=arr[i_left+1]<=arr[i_right] */
      i_middle = (i_left+i_right)>>1;
      SWAP(arr[i_middle], arr[i_left+1]);
      SWAP2(brr[i_middle], brr[i_left+1]);
      if (arr[i_left] > arr[i_right]) {
	SWAP(arr[i_left], arr[i_right]);
	SWAP2(brr[i_left], brr[i_right]);
      }
      if (arr[i_left+1] > arr[i_right]) {
	SWAP(arr[i_left+1], arr[i_right]);
	SWAP2(brr[i_left+1], brr[i_right]);
      }
      if (arr[i_left] > arr[i_left+1]) {
	SWAP(arr[i_left], arr[i_left+1]);
	SWAP2(brr[i_left], brr[i_left+1]);
      }

      /* initialise partitioning pointers */
      i=i_left+1;
      j=i_right;
      
      /* pick partitioning element */
      partelem_a = arr[i_left+1];
      partelem_b = brr[i_left+1];

      /* partitioning exchange loop */
      for(;;) {
      	do i++; while (arr[i]<partelem_a);
      	do j--; while (arr[j]>partelem_a);
      	if(j<i) break; /* partitioning pointers crossed */
      	SWAP(arr[i],arr[j]); /* exchange */
      	SWAP2(brr[i],brr[j])
      }
      /* insert partitioning element */
      arr[i_left+1] = arr[j];
      brr[i_left+1] = brr[j];
      arr[j] = partelem_a;
      brr[j] = partelem_b;

      /* push the larger of the two partitions on stack */
      stack_size += 2;
      if(stack_size>MAXSTACKSIZE) return ERRCODE_SORTSTACK; /* stack overflow */

      if (i_right-i+1 >= j-i_left) {
	index_stack[stack_size] = i_right;
	index_stack[stack_size-1] = i;
	i_right = j-1;
      } else {
	index_stack[stack_size] = j-1;
	index_stack[stack_size-1] = i_left;
	i_left = i;
      }
    } /* if (i_right-i_left < MINARRSIZE) else */
  } /* loop over all subarrays */
}


void sortMultiKey(const char *x[], int num, int maxdepth)
{
#ifdef sort_debug
  int level = 0;
#endif
  mksort(
#ifdef sort_debug
	 &level,
#endif
	 x, num, 0, maxdepth);
}

void sortSufficesByMultiKeyQuickSort(const char *hstrp, uint32_t *sfxidxp, 
				     uint32_t nsfx, short startdepth, short maxdepth)
{
#ifdef sort_debug
  int level = 0;
#endif
  mkeyQSortSuffix(
#ifdef sort_debug
	 &level,
#endif
	 hstrp, sfxidxp, nsfx, startdepth, maxdepth);
}
