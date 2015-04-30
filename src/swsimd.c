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

#include <string.h>
#include <limits.h>
#include "swsimd.h"
#include "alibuffer_struct.h"

enum {
#ifdef SCORE_SIMD_IMIC
  NELEM_REGISTER = 16,
  PERMUTMASK = 0xFFFE,
#else
  BYTEMASK = 0x00ff,
  NBITS_PER_BYTE = 8,
#endif
};

typedef unsigned char UCHAR;

#ifdef SCORE_SIMD_IMIC
static __attribute__((align(64))) const int SIMD_LEFTSHIFT_PERMXV[NELEM_REGISTER] = 
  {0,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14};
#else
static short BIAS = ((short) 1)<<15;
#endif

/******************************************************************************
 ********************************** Macros ************************************
 ******************************************************************************/

#ifdef SCORE_SIMD_IMIC

#define SIMD_ZERO_VECTOR(v) \
  (v) = _mm512_setzero_epi32()

#define SIMD_COPY_SCALAR_TO_VECTOR(a, v) \
  (v) = _mm512_set1_epi32((int) (a))

#define SIMD_COPY_VECTOR(a, v)\
  (v) = _mm512_shuffle_epi32((a), 0xe4)

#define SIMD_MAXELEMV(vOrig, vMax, vTmp)	 \
  vTmp = _mm512_permute4f128_epi32(vOrig, 0xB1); \
  vMax = _mm512_max_epi32(vOrig, vTmp);		 \
  vTmp = _mm512_permute4f128_epi32(vMax, 0x1B);	 \
  vMax = _mm512_max_epi32(vMax, vTmp);		 \
  vTmp = _mm512_shuffle_epi32(vMax, 0xB1);	 \
  vMax = _mm512_max_epi32(vMax, vTmp);		 \
  vTmp = _mm512_shuffle_epi32(vMax, 0x1B);	 \
  vMax = _mm512_max_epi32(vMax, vTmp);

#define SIMD_MAXELEMV_DESTRUCT(vMax, vTmp)	\
  SIMD_MAXELEMV(vMax, vMax, vTmp)
/* overwrites vMax */

#else /* #ifdef SCORE_SIMD_IMIC */

#define SIMD_ZERO_VECTOR(v) \
  (v) = _mm_setzero_si128()

/* memset to stop compiler/leak checker comlaining about uninitialised values */
/* #define SIMD_ZERO_VECTOR(v) \ */
/*   memset(&(v), 0, NBYTES_REGISTER); \ */
/*   (v) = _mm_xor_si128 ((v), (v)) */

#define SIMD_COPY_SCALAR_TO_VECTOR(b, v) \
  SIMD_ZERO_VECTOR((v)); \
  (v) = _mm_insert_epi16 ((v), b, 0); \
  (v) = _mm_shufflelo_epi16 ((v), 0); \
  (v) = _mm_shuffle_epi32 ((v), 0);

#define SIMD_NEGBIAS_VECTOR(v) \
  SIMD_ZERO_VECTOR((v)); \
  (v) = _mm_cmpeq_epi16 ((v), (v));\
  (v) = _mm_slli_epi16 ((v), 15); 

#define COPY_BYTE_TO_REGISTER_VALUE(b, v) \
  SIMD_ZERO_VECTOR((v)); \
  intval = ((b) << NBITS_PER_BYTE) | (b & BYTEMASK); \
  (v) = _mm_insert_epi16 ((v), intval, 0); \
  (v) = _mm_shufflelo_epi16 ((v), 0); \
  (v) = _mm_shuffle_epi32 ((v), 0);

#define SIMD_MAXELEMV_DESTRUCT(vMax, vTmp)    \
  vTmp = _mm_srli_si128 (vMax, 8);	      \
  vMax = _mm_max_epi16 (vMax, vTmp);	      \
  vTmp = _mm_srli_si128 (vMax, 4);	      \
  vMax = _mm_max_epi16 (vMax, vTmp);	      \
  vTmp = _mm_srli_si128 (vMax, 2);	      \
  vMax = _mm_max_epi16 (vMax, vTmp);

#define SIMD_MAXELEMV(vOrig, vMax, vTmp)	\
  vMax = _mm_shufflehi_epi16 (vOrig, 0xE4);	\
  SIMD_MAXELEMV_DESTRUCT(vMax, vTmp);

#define SIMD_MAXELEMV_8BIT_DESTRUCT(vMax, vTmp)\
  vTmp = _mm_srli_si128 (vMax, 8);	       \
  vMax = _mm_max_epu8 (vMax, vTmp);	       \
  vTmp = _mm_srli_si128 (vMax, 4);	       \
  vMax = _mm_max_epu8 (vMax, vTmp);	       \
  vTmp = _mm_srli_si128 (vMax, 2);	       \
  vMax = _mm_max_epu8 (vMax, vTmp);	       \
  vTmp = _mm_srli_si128 (vMax, 1);	       \
  vMax = _mm_max_epu8 (vMax, vTmp);

#define SIMD_MAXELEMV_8BIT(vOrig, vMax, vTmp)	\
  vMax = _mm_shufflehi_epi16 (vOrig, 0xE4);	\
  SIMD_MAXELEMV_8BIT_DESTRUCT(vMax, vTmp);

#endif


/******************************************************************************
 ********************** Private SIMD Alignment Methods ************************
 ******************************************************************************/
#ifdef alignment_matrix_debug
#ifdef SCORE_SIMD_IMIC
void printfArr(const int *vp, const char *label, int bias)
{
  int i;
  if (NULL == label)
    printf("V:{%i", vp[0] - bias);
  else
    printf("%s {%i", label, vp[0] - bias);
  for (i=1; i<NELEM_REGISTER; i++)
    printf(",%i", vp[i] - bias);
  printf("}\n");
}

void printfSIMDV(const __m512i *vp, int *oup, const char *label, int bias)
{
  _mm512_store_epi32(oup, *vp);
  printfArr(oup, label, bias);
}
#endif
#ifdef SCORE_SIMD_IMIC
static void printfStripedIntVector(int *vp, unsigned int qlen, int segsiz)
{
  unsigned int j , k, nseg;

  if (segsiz <= 0)
    return;
  nseg = (qlen + segsiz - 1)/segsiz;
  
  for (j=0; j<nseg; j++)
    for (k = (unsigned int) j; k<qlen; k += SCORSIMD_NINTS)
      printf("%3i|", vp[k]);
  //printf("%3u|", k);
}
#else
static void printfStripedShortVector(short *vp, unsigned int qlen, int segsiz)
{
  unsigned int j , k, nseg;

  if (segsiz <= 0)
    return;
  nseg = (qlen + segsiz - 1)/segsiz;
  
  for (j=0; j<nseg; j++)
    for (k = (unsigned int) j; k<qlen; k += SCORSIMD_NSHORTS)
      printf("%3hi|", vp[k]-BIAS);
      //printf("%3u|", k);
}

static void printfStripedByteVector(unsigned char *vp, unsigned int qlen, int segsiz)
{
  unsigned int j , k, nseg;

  if (segsiz <= 0)
    return;
  nseg = (qlen + segsiz - 1)/segsiz;
  
  for (j=0; j<nseg; j++)
    for (k = (unsigned int) j; k<qlen; k += SCORSIMD_NBYTES)
      printf("%3hu|", (unsigned short) vp[k]);
      //printf("%3u|", k);
}
#endif
#endif

#ifdef SCORE_SIMD_IMIC
static int alignSmiWatIntStriped(unsigned int *maxscor,
				 const AliBuffer *abp, 
				 const ScoreProfile *spp,
#ifdef alignment_matrix_debug
				 const char *psqp,
				 const SeqCodec *codecp,
#endif
				 const char *usqp,
				 int uslen)
     /* adapted from M. Farrar (2007) Bioinformatics 23, 156 - 161.
      * http://farrar.michael.googlepages.com/Smith-waterman
      */
{
  int errcode = ERRCODE_SUCCESS;
  int i, j, segsiz, score, cmpval;
  short tmp;
  unsigned short gap_init, gap_ext;
  const __m512i *vProfp;
  __m512i vH, vE, vF, vMax, vMin, vGapI, vGapE, vTmp, *vp;
  __m512i *vEp = abp->Ev;
  __m512i *vHSp = abp->H1v;
  __m512i *vHLp = abp->H2v;
  __mmask16 cmpmsk;
  const __m512i vSLX = _mm512_load_epi32(&SIMD_LEFTSHIFT_PERMXV);
  const __mmask16 mPERM = _mm512_int2mask(PERMUTMASK);
  const __mmask16 mSEL1 = _mm512_int2mask(1);
#ifdef alignment_matrix_debug
  char unprof_symbol;
  const char *decoderp = seqCodecGetDecoder(codecp, NULL);
  unsigned int maximum_score = 0;
  unsigned int q, qlen;
  __m512i vMaxBuf, vTmpBuf;
  /* int lazy_Floop_ctr;  */
  int *oubuf_allocp = (int *) malloc(128);
  int *oubufp = (int *) SCORE_ALIGN_MEMORY(oubuf_allocp);
  if (NULL == oubuf_allocp)
    return ERRCODE_NOMEM;
#endif



  *maxscor = 0;
  if (!((vEp) && (vHSp) && (vHLp)))
    return ERRCODE_ASSERT;

  vProfp = (const __m512i *) 
    scoreGetStripedProfile(NULL, 
#ifdef alignment_matrix_debug
			   &qlen,
#else
			   NULL, 
#endif
			   &gap_init, &gap_ext, NULL, &segsiz,
			   SCORPROF_STRIPED_32, spp);
  if (!vProfp)
    return ERRCODE_SWATSTRIP;

  /* Set SSE constants */
  SIMD_COPY_SCALAR_TO_VECTOR(gap_init, vGapI);
  SIMD_COPY_SCALAR_TO_VECTOR(gap_ext, vGapE);
  

  /* initialize storage vector to 0 */
  SIMD_ZERO_VECTOR(vMax);
  SIMD_ZERO_VECTOR(vMin);

  for (i=0; i<segsiz; i++) {
    _mm512_store_epi32(vEp + i, vMax);
    _mm512_store_epi32(vHSp + i, vMax);
  }

#ifdef alignment_matrix_debug
  printf("allocation address vProfp: %p\n", vProfp);
  printfSIMDV(&vMax, oubufp, "vMax", 0);
  printfSIMDV(&vMin, oubufp, "vMin", 0);
#endif
  for (i=0; i<uslen; i++) {
    const __m512i *vScorep = vProfp + (usqp[i]&SEQCOD_ALPHA_MASK) * segsiz;

    /* zero out F */
    SIMD_ZERO_VECTOR(vF);
    
    /* load the next h value,
     * initialise element on right to conceptual 0 */
    vH = _mm512_load_epi32(vHSp + segsiz - 1);
    vH = _mm512_mask_permutevar_epi32(vMin,
				      mPERM,
				      vSLX,
				      vH);
#ifdef alignment_matrix_debug
    //printfSIMDV(&vH, oubufp, "vH", 0);
#endif
    /* swap the two H vectors */
    vp = vHLp;
    vHLp = vHSp;
    vHSp = vp;

    for (j = 0; j < segsiz; j++) {
      /* load values of vF and vH from previous row (one unit up) */
      vE = _mm512_load_epi32(vEp + j);

      vTmp = _mm512_load_epi32(vScorep + j);

      /* add score to vH */
      vH = _mm512_add_epi32(vH, vTmp);

      /* Update highest score encountered this far */
      vMax = _mm512_max_epi32(vMax, vH);

#ifdef alignment_matrix_debug
      /*       //printfSIMDV(&vMax, oubufp, "vMax", 0); */      
      /* find largest score in the vMax vector */
      SIMD_MAXELEMV(vMax, vMaxBuf, vTmpBuf);
      _mm512_mask_packstorelo_epi32(&score, mSEL1, vMaxBuf);

      /* //printf("sse_score = %i\n", score); */
      if ((unsigned int) score > maximum_score) {
      	maximum_score = (unsigned int) score;
	printf("IMIC_max_scor(%i,%i) %u\n", i, j, maximum_score);
      }
#endif

      /* get max from vH, vE and vF */
      vH = _mm512_max_epi32 (vH, vE);
      vH = _mm512_max_epi32 (vH, vF);
      vH = _mm512_max_epi32 (vH, vMin);

      /* save vH values */
      _mm512_store_epi32 (vHSp + j, vH);

      /* update vE value */
      vH = _mm512_sub_epi32 (vH, vGapI);
      vE = _mm512_sub_epi32 (vE, vGapE);
      vE = _mm512_max_epi32 (vE, vH);

      /* update vF value */
      vF = _mm512_sub_epi32 (vF, vGapE);
      vF = _mm512_max_epi32 (vF, vH);

      /* save vE values */
      _mm512_store_epi32 (vEp + j, vE);

      /* load the next h value */
      vH = _mm512_load_epi32(vHLp + j);
    } /* for (j = 0; j < segsiz; j++) */

    /* reset pointers to the start of the saved data */
    j = 0;
    vH = _mm512_load_epi32(vHSp + j);

    /*  the computed vF value is for the given column.  since */
    /*  we are at the end, we need to shift the vF value over */
    /*  to the next column. */
    vF = _mm512_mask_permutevar_epi32(vMin,
				      mPERM,
				      vSLX,
				      vF);

    vTmp = _mm512_sub_epi32 (vH, vGapI);
    cmpmsk = _mm512_cmp_epi32_mask (vF, vTmp, _MM_CMPINT_GT);
    cmpval = _mm512_mask2int (cmpmsk);

#ifdef alignment_matrix_debug
    /* printf("cmpval[%i,%i]=%x\n", i,j, cmpval); */

    /* lazy_Floop_ctr = 0; */
    while (cmpval != 0x0000
	   /* && lazy_Floop_ctr < 32 */
	   ) { 
#else
    while (cmpval != 0x0000) {
#endif
	/* lazy F-loop */
      vE = _mm512_load_epi32 (vEp + j);
      vH = _mm512_max_epi32 (vH, vF);
      vH = _mm512_max_epi32 (vH, vMin);

      /* save vH values */
      _mm512_store_epi32 (vHSp + j, vH);
      
      /*  update vE in case the new vH value would change it */
      vH = _mm512_sub_epi32 (vH, vGapI);
      vE = _mm512_max_epi32 (vE, vH);
      _mm512_store_epi32 (vEp + j, vE);

      /* update vF value */
      vF = _mm512_sub_epi32 (vF, vGapE);

      j++;
      if (j >= segsiz) {
    	j = 0;
    	/* shift left by 1 element,
    	 * initialise element 0 to conceptual 0 */
    	vF = _mm512_mask_permutevar_epi32(vMin,
    					  mPERM,
    					  vSLX,
    					  vF);
     }

      vH = _mm512_load_epi32 (vHSp + j);

      vTmp = _mm512_sub_epi32 (vH, vGapI);
      cmpmsk = _mm512_cmp_epi32_mask (vF, vTmp, _MM_CMPINT_GT);
      cmpval = _mm512_mask2int (cmpmsk);
#ifdef alignment_matrix_debug
      /* printf("cmpval[%i]=%x\n", lazy_Floop_ctr, cmpval); */
      /* lazy_Floop_ctr++; */
#endif
  }
#ifdef alignment_matrix_debug
    printf("  [%i,%i]: ", i, 0);
    for (q=0; q<qlen; q++)
      printf("%4i", q);
    unprof_symbol = decoderp[(UCHAR) usqp[i]];
    printf("\n%c [%i,%i]: |", unprof_symbol, i, 0);
    for (q=0; q<qlen; q++)
      printf(" %c |", decoderp[(UCHAR) psqp[q]]);
    printf("\nHp[%i,%i]: |", i, 0);
    printfStripedIntVector((int *) vHSp, qlen, segsiz);
    printf("\n\n");
#endif
  }

  /* find largest score in the maxscorv vector */

  SIMD_MAXELEMV_DESTRUCT(vMax, vTmp);
  _mm512_mask_packstorelo_epi32(&score, mSEL1, vMax);
  
  *maxscor = (unsigned int) score;
#ifdef alignment_matrix_debug  
  free(oubuf_allocp);
#endif
  return errcode;
}

#else /* #ifdef SCORE_SIMD_IMIC */
static int alignSmiWatShortStriped(unsigned short *maxscor,
			    const AliBuffer *abp, 
			    const ScoreProfile *spp,
#ifdef alignment_matrix_debug
			    const char *psqp,
			    const SeqCodec *codecp,
#endif
			    const char *usqp,
			    int uslen)
     /* Striped Smith-Waterman using SSE2 instructions.
      * adapted from M. Farrar (2007) Bioinformatics 23, 156 - 161.
      * http://farrar.michael.googlepages.com/Smith-waterman
      */
{
  int errcode, i, j, segsiz, score, cmpval;
  short tmp;
  unsigned short gap_init, gap_ext;
  const SIMDV_t *vScorep, *vProfp;
  SIMDV_t vH, vE, vF, vMax, vMin, vGapI, vGapE, vTmp, *vp;
  SIMDV_t *vEp = abp->Ev;
  SIMDV_t *vHSp = abp->H1v;
  SIMDV_t *vHLp = abp->H2v;
#ifdef alignment_matrix_debug
  char unprof_symbol;
  const char *decoderp = seqCodecGetDecoder(codecp, NULL);
  unsigned short maximum_score = 0;
  unsigned int q, qlen;
  SIMDV_t vMaxBuf, vTmpBuf;
#endif


  *maxscor = 0;
  if (!((vEp) && (vHSp) && (vHLp)))
    return ERRCODE_ASSERT;

  vProfp = (const SIMDV_t *) scoreGetStripedProfile(NULL, 
#ifdef alignment_matrix_debug
				  &qlen,
#else
				  NULL, 
#endif
				  &gap_init, &gap_ext, NULL, &segsiz,
				  SCORPROF_STRIPED_16, spp);
  if (!vProfp)
    return ERRCODE_SWATSTRIP;
#ifdef alignment_matrix_debug
  printf("alignSmiWatShortStriped: segsiz = %i, bias=%hi\n", segsiz, BIAS);
#endif
  /* Set SSE constants */
  /* should be able to do this with _mm_set1_epi16(short b) */
  SIMD_COPY_SCALAR_TO_VECTOR(gap_init, vGapI);
  SIMD_COPY_SCALAR_TO_VECTOR(gap_ext, vGapE);
  
  /*  load vMaxScore with the zeros.  since we are using signed */
  /*  math, we will bias the maxscore to -32768 so we have the */
  /*  full range of the short.
   */
  SIMD_NEGBIAS_VECTOR(vMax)

  /* initialize elements of vMin,
   * conceptual 0 -> byte 0,1 (1st element), 0 for the rest,
   * to be able to shift and set first element of a vector v using or
   *   v = _mm_slli_si128 (v, 2);
   *   v = _mm_or_si128 (v, vMin);
   */
  vMin = _mm_shuffle_epi32 (vMax, 0);
  vMin = _mm_srli_si128 (vMin, 14); /* shift right by 14 bytes, fill with 0s */

  /* initialize storage vector to 0 (biased to -32768) */
  for (i=0; i<segsiz; i++) {
    _mm_store_si128 (vEp + i, vMax);
    _mm_store_si128 (vHSp + i, vMax);
  }

  for (i=0; i<uslen; i++) {
    vScorep = vProfp + (usqp[i]&SEQCOD_ALPHA_MASK) * segsiz;
    
    /* zero out F */
    SIMD_NEGBIAS_VECTOR(vF);
    
    /* load the next h value */
    vH = _mm_load_si128 (vHSp + segsiz - 1);
    vH = _mm_slli_si128 (vH, 2); /* shift left 2 bytes (for short) */
    vH = _mm_or_si128 (vH, vMin); /* initialise new short on right to -32768 (conceptual 0) */
	
    /* swap the two H vectors */
    vp = vHLp;
    vHLp = vHSp;
    vHSp = vp;

    for (j = 0; j < segsiz; j++) {
      /* load values of vF and vH from previous row (one unit up) */
      vE = _mm_load_si128 (vEp + j);

      /* add score to vH */
      vTmp = _mm_load_si128 (vScorep + j);
      vH = _mm_adds_epi16 (vH, vTmp);

      /* Update highest score encountered this far */
      vMax = _mm_max_epi16 (vMax, vH);

#ifdef alignment_matrix_debug
      /* find largest score in the maxscorv vector */
      SIMD_MAXELEMV(vMax, vMaxBuf, vTmpBuf);
      score = _mm_extract_epi16 (vMaxBuf, 0);
      score += BIAS;
      //printf("sse_score = %i\n", score);
      if (score > maximum_score) {
	maximum_score = score;
	printf("sse_max_scor_shrt(%i,%i) = %i\n", i, j, maximum_score);
      }
#endif

      /* get max from vH, vE and vF */
      vH = _mm_max_epi16 (vH, vE);
      vH = _mm_max_epi16 (vH, vF);

      /* save vH values */
      _mm_store_si128 (vHSp + j, vH);

      /* update vE value */
      vH = _mm_subs_epi16 (vH, vGapI);
      vE = _mm_subs_epi16 (vE, vGapE);
      vE = _mm_max_epi16 (vE, vH);

      /* update vF value */
      vF = _mm_subs_epi16 (vF, vGapE);
      vF = _mm_max_epi16 (vF, vH);

      /* save vE values */
      _mm_store_si128 (vEp + j, vE);

      /* load the next h value */
      vH = _mm_load_si128 (vHLp + j);
    }

    /* reset pointers to the start of the saved data */
    j = 0;
    vH = _mm_load_si128 (vHSp + j);

    /*  the computed vF value is for the given column.  since */
    /*  we are at the end, we need to shift the vF value over */
    /*  to the next column. */
    vF = _mm_slli_si128 (vF, 2);
    vF = _mm_or_si128 (vF, vMin);
    vTmp = _mm_subs_epi16 (vH, vGapI);
    vTmp = _mm_cmpgt_epi16(vF, vTmp);
    cmpval  = _mm_movemask_epi8 (vTmp); /* Creates a 16-bit mask from the most significant bits of 
					 * the 16 signed or unsigned 8-bit integers in vTemp and 
					 * zero extends the upper bits. */

    /* lazy F-loop */
    while (cmpval != 0x0000) { /* not all 8-bit integers == 0 */
      vE = _mm_load_si128 (vEp + j);
      vH = _mm_max_epi16 (vH, vF);

      /* save vH values */
      _mm_store_si128 (vHSp + j, vH);
      
      /*  update vE incase the new vH value would change it */
      vH = _mm_subs_epi16 (vH, vGapI);
      vE = _mm_max_epi16 (vE, vH);
      _mm_store_si128 (vEp + j, vE);

      /* update vF value */
      vF = _mm_subs_epi16 (vF, vGapE);

      j++;
      if (j >= segsiz) {
	j = 0;
	vF = _mm_slli_si128 (vF, 2);
	vF = _mm_or_si128 (vF, vMin);
      }

      vH = _mm_load_si128 (vHSp + j);

      vTmp = _mm_subs_epi16 (vH, vGapI);
      vTmp = _mm_cmpgt_epi16 (vF, vTmp);
      cmpval  = _mm_movemask_epi8 (vTmp);
    }
#ifdef alignment_matrix_debug
    printf("  [%i,%i]: ", i, 0);
    for (q=0; q<qlen; q++)
      printf("%4i", q);
    unprof_symbol = decoderp[(UCHAR) usqp[i]];
    printf("\n%c [%i,%i]: |", unprof_symbol, i, 0);
    for (q=0; q<qlen; q++)
      printf(" %c |", decoderp[(UCHAR) psqp[q]]);
    printf("\nHp[%i,%i]: |", i, 0);
    printfStripedShortVector((short *) vHSp, qlen, segsiz);
    printf("\n\n");
#endif
  }

  /* find largest score in the maxscorv vector */

  /* return */
  SIMD_MAXELEMV_DESTRUCT(vMax, vTmp);
  tmp = _mm_extract_epi16 (vMax, 0);
  score = ((int) tmp) - BIAS;

  if (score >= USHRT_MAX)
    errcode = ERRCODE_SWATEXCEED;
  else if (score < 0)
    errcode = ERRCODE_ASSERT;
  else {
    errcode = ERRCODE_SUCCESS;
    *maxscor = (unsigned short) score;
  }
  
  return errcode;
}

static int alignSmiWatByteStriped(UCHAR *maxscor,
			   const AliBuffer *abp,
			   const ScoreProfile *spp,
#ifdef alignment_matrix_debug
			    const char *psqp,
			    const SeqCodec *codecp,
#endif
			   const char *usqp,
			   int uslen)
     /* Striped Smith-Waterman using SSE2 instructions.
      * adapted from M. Farrar (2007) Bioinformatics 23, 156 - 161.
      * http://farrar.michael.googlepages.com/Smith-waterman
      *
      * \param maxscor Returns maximum score
      * \param abp Buffers used for dynamic programming
      * \param ssp Sequence profile (query).
      * \param psqp Profiled (query) sequence in SEQCOD_MANGLED encoding.
      * \param codec Sequence de/encoder.
      * \param usqp Unprofiled (subject) sequence in SEQCOD_MANGLED encoding.
      * \param uslen Length of the unprofiled (subject) sequence.
      */
{
  int errcode, i, j, score = 0, intval, segsiz, cmpval = 0;
  unsigned short gap_init, gap_ext, bias = 0;
  const __m128i *vScorep, *vProfp;
  __m128i vH, vE, vF, vMax, vBias, vZero, vGapI, vGapE, vTmp, *vp;
  __m128i *vEp = abp->Ev;
  __m128i *vHSp = abp->H1v;
  __m128i *vHLp = abp->H2v;
#ifdef alignment_matrix_debug
  char unprof_symbol;
  const char *decoderp = seqCodecGetDecoder(codecp, NULL);
  int maximum_score = 0;
  unsigned int q, qlen;
  __m128i vMaxBuf, vTmpBuf;
#endif

  *maxscor = 0;
  if (!((vEp) && (vHSp) && (vHLp)))
    return ERRCODE_ASSERT;

  vProfp = (const __m128i *) scoreGetStripedProfile(NULL, 
#ifdef alignment_matrix_debug
				  &qlen,
#else
				  NULL, 
#endif
				  &gap_init, &gap_ext, &bias, &segsiz,
				  SCORPROF_STRIPED_8, spp);
  if (!vProfp)
    return ERRCODE_SWATSTRIP;

  /* Load constants to SSE2 register value */
  /* should be able to do this with _mm_set1_epi8(char b) */
  COPY_BYTE_TO_REGISTER_VALUE(bias, vBias);
  COPY_BYTE_TO_REGISTER_VALUE(gap_init, vGapI);
  COPY_BYTE_TO_REGISTER_VALUE(gap_ext, vGapE);
  
  /* variables initialised to 0 */
  SIMD_ZERO_VECTOR(vMax);
  SIMD_ZERO_VECTOR(vZero);

  /* Zero out the storage vector */
  for (i=0; i<segsiz; i++) {
    _mm_store_si128 (vEp + i, vZero);
    _mm_store_si128 (vHSp + i, vZero);
  }

  for (i=0; i<uslen; i++) {
    vScorep = vProfp + (usqp[i]&SEQCOD_ALPHA_MASK) * segsiz;
    
    /* zero out F */
    SIMD_ZERO_VECTOR(vF);
    
    /* load the next h value */
    vH = _mm_load_si128 (vHSp + segsiz - 1);
    vH = _mm_slli_si128 (vH, 1);
	
    /* swap the two H vectors */
    vp = vHLp;
    vHLp = vHSp;
    vHSp = vp;

    for (j = 0; j < segsiz; j++) {
      /* load values of vF and vH from previous row (one unit up) */
      vE = _mm_load_si128 (vEp + j);

      /* add score to vH */
      vTmp = _mm_load_si128 (vScorep + j);
      vH = _mm_adds_epu8 (vH, vTmp);
      vH = _mm_subs_epu8 (vH, vBias);

      /* Update highest score encountered this far */
      vMax = _mm_max_epu8 (vMax, vH);

#ifdef alignment_matrix_debug
      /* find largest score in the maxscorv vector */
      SIMD_MAXELEMV_8BIT(vMax, vMaxBuf, vTmpBuf);
      
      score = _mm_extract_epi16 (vMaxBuf, 0);
      score &= BYTEMASK;
      //printf("sse_score = %i\n", score);
      if (score > maximum_score) {
	maximum_score = score;
	printf("sse_max_scor(%i,%i) = %i\n", i, j, maximum_score);
      }
#endif

      /* get max from vH, vE and vF */
      vH = _mm_max_epu8 (vH, vE);
      vH = _mm_max_epu8 (vH, vF);

      /* save vH values */
      _mm_store_si128 (vHSp + j, vH);

      /* update vE value */
      vH = _mm_subs_epu8 (vH, vGapI);
      vE = _mm_subs_epu8 (vE, vGapE);
      vE = _mm_max_epu8  (vE, vH);

      /* update vF value */
      vF = _mm_subs_epu8 (vF, vGapE);
      vF = _mm_max_epu8 (vF, vH);

      /* save vE values */
      _mm_store_si128 (vEp + j, vE);

      /* load the next h value */
      vH = _mm_load_si128 (vHLp + j);
    } /* for (j = 0; j < segsiz; j++) */

    /* reset pointers to the start of the saved data */
    j = 0;
    vH = _mm_load_si128 (vHSp + j);

    /*  the computed vF value is for the given column.  since */
    /*  we are at the end, we need to shift the vF value over */
    /*  to the next column. */
    vF = _mm_slli_si128 (vF, 1);
    vTmp = _mm_subs_epu8 (vH, vGapI);
    vTmp = _mm_subs_epu8 (vF, vTmp);
    vTmp = _mm_cmpeq_epi8 (vTmp, vZero); /* Where 8-bit integers are equal, all bits are set (0xff)
					  * otherwise unset (0x00). */
    cmpval = _mm_movemask_epi8 (vTmp);  /* Creates a 16-bit mask from the most significant bits of 
					 * the 16 signed or unsigned 8-bit integers in vTemp and 
					 * zero extends the upper bits. */

    /* lazy F-loop */
    while (cmpval != 0xffff) { /* not all 8-bit integers == 0 */
      vE = _mm_load_si128 (vEp + j);
      vH = _mm_max_epu8 (vH, vF);

      /* save vH values */
      _mm_store_si128 (vHSp + j, vH);
      
      /*  update vE incase the new vH value would change it */
      vH = _mm_subs_epu8 (vH, vGapI);
      vE = _mm_max_epu8 (vE, vH);
      _mm_store_si128 (vEp + j, vE);

      /* update vF value */
      vF = _mm_subs_epu8 (vF, vGapE);

      j++;
      if (j >= segsiz) {
	j = 0;
	vF = _mm_slli_si128 (vF, 1);
      }

      vH = _mm_load_si128 (vHSp + j);

      vTmp = _mm_subs_epu8 (vH, vGapI);
      vTmp = _mm_subs_epu8 (vF, vTmp);
      vTmp = _mm_cmpeq_epi8 (vTmp, vZero);
      cmpval  = _mm_movemask_epi8 (vTmp);
    }
#ifdef alignment_matrix_debug
    printf("  [%i,%i]: ", i, 0);
    for (q=0; q<qlen; q++)
      printf("%4i", q);
    unprof_symbol = decoderp[(UCHAR) usqp[i]];
    printf("\n%c [%i,%i]: |", unprof_symbol, i, 0);
    for (q=0; q<qlen; q++)
      printf(" %c |", decoderp[(UCHAR) psqp[q]]);
    printf("\nHp[%i,%i]: |", i, 0);
    printfStripedByteVector((unsigned char *) vHSp, qlen, segsiz);
    printf("\n\n");
#endif
  }

  /* find largest score in the maxscorv vector */
  SIMD_MAXELEMV_8BIT_DESTRUCT(vMax, vTmp);

  /* store in temporary variable */
  score = _mm_extract_epi16 (vMax, 0);
  score &= BYTEMASK;

  /*  check if we might have overflowed */
  if (score + bias >= UCHAR_MAX) {
    errcode = ERRCODE_SWATEXCEED;
  } else {
    errcode = ERRCODE_SUCCESS;
    *maxscor = (UCHAR) score; /* return largest score */
  }
  return errcode;
}
#endif /* #ifdef SWSIMD_IMIC #else */

/******************************************************************************
 ********************** Public SIMD Alignment Methods *************************
 ******************************************************************************/

int swSIMDAlignStriped(ALIDPMSCOR_t *maxscor, 
		       const AliBuffer *abp,
		       const ScoreProfile *profp, 
#ifdef alignment_matrix_debug 
		       const SeqCodec *codecp,
		       const char *profiled_seqp,
#endif
		       const char *unprofiled_seqp,
		       int unprofiled_seqlen)
{
  int errcode;
#ifdef SCORE_SIMD_IMIC
  unsigned int scor;
#else 
  unsigned char bytscor;
  unsigned short shortscor;
#endif

  *maxscor = 0;

#ifdef SCORE_SIMD_IMIC
  errcode = alignSmiWatIntStriped(&scor, 
				  abp, 
				  profp,
#ifdef alignment_matrix_debug
				  profiled_seqp,
				  codecp,
#endif 
				  unprofiled_seqp,
				  unprofiled_seqlen);
  if (!errcode) {
    if (scor > INT_MAX)
      errcode = ERRCODE_SWATEXCEED;
    else
      *maxscor = (int) scor;
  }

#else
  errcode = alignSmiWatByteStriped(&bytscor, 
				   abp, 
				   profp,
#ifdef alignment_matrix_debug
				   profiled_seqp,
				   codecp,
#endif 
				   unprofiled_seqp,
				   unprofiled_seqlen);
  if (errcode == ERRCODE_SWATEXCEED) {
    errcode = alignSmiWatShortStriped(&shortscor, 
				      abp, 
				      profp,
#ifdef alignment_matrix_debug
				      profiled_seqp,
				      codecp,
#endif 
				      unprofiled_seqp,
				      unprofiled_seqlen);
    if (!errcode) 
      *maxscor = (int) shortscor;
  } else if (!errcode) {
    *maxscor = (int) bytscor;
  }
#endif 

  return errcode;
}

