/** Statistics on insert length of paired reads */

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

#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>

#include "elib.h"
#include "vdef.h"
#include "sort.h"
#include "insert.h"

enum {
  INTERVAL_Z = 3,  /**< Insert sizes within INTERVAL_Z standard deviations
		    * of median are considered */
  HISTO_PRINTCHAR = '*', /**< Character for printing histogram */
  PRINTCHAR_NEWLIN = '\n',
  DEFAULT_HISTO_WIDTH = 128,
  KERNEL_CUTOFF_BANDFAC = 3, /**< Cutoff in smoothing kernel in multiples of the
			       * band width */
  KERNEL_MIN_WIDTH = 3,    /**< Minimum kernel band width */
  IOFIL_LINBUFSIZ = 128,   /**< line buffer for text based IO */
  IQR_RANGE_FAC = 3,       /**< Use iqr*IQR_RANGE_FAC as sampling range */
  HISTO_MIN_BINNUM = 16,
  HISTO_MAX_BINNUM = 1028,
  NNEWLINE_CHAR = '\n',
  SAMPLE_TARGETSIZ = 4098,
};

enum INSERT_FLAGS {
  INSFLG_EMPTY = 0,
  INSFLG_SEEDED = 1,
  INSFLG_SMOOTHED = 2,
};

typedef unsigned char BOOL_t;
typedef unsigned char INSFLG_t;

typedef V_STRUCT(int32_t) SAMPLE;

struct _InsSample {
  SAMPLE *pSample;  /**< Sample of reads */
  int readival;
};

struct _InsHist { /**< Contains Histogram of insert size distributions */
  INSFLG_t status;   /**< One of INSERT_FLAGS */
  int32_t *pCounts;  /**< Counts for insert sizes */
  int32_t *pSmoothCounts; /**< Smoothed counts */
  double *pKernelBuf;/**< Buffer for smoothing kernel */
  int32_t iSpan;     /**< Number of cells in histogram (size of array pCounts */
  int32_t iInSizLo;  /**< Insert size of the lowest cell 0. */     
  int32_t iInSizHi;  /**< Insert size of the highest cell iSpan -1. */
  int32_t iScalFac;   /**< Scaling factor for histogram */
  uint64_t iNum;      /**< Total number if insert sizes */
  int32_t median;     /**< Median of the distribution */
  int32_t quart_lo;   /**< lower quartile */
  int32_t quart_hi;   /**< upper quartile */
};

static const char IOFIL_HEADER[] = "# SMALT histogram of insert sizes\n";
static const char IOFIL_KEY_START[] = "HISTO_START";
static const char IOFIL_KEY_FORMAT[] = "HISTO_BINNUM %i\nHISTO_SCALFAC %i\nHISTO_INSIZLO %i\n"\
  "HISTO_INSIZHI %i\nHISTO_TOTNUM %llu\nHISTO_QUARTILES %i %i %i\n";
static const char IOFIL_KEY_END[] = "HISTO_END";
static const char IOFIL_FORMAT[] = "%i %i\n";

/******************************************************************************
 ********************************** Macros ************************************
 ******************************************************************************/

#define CALC_HISTO_IDX(idx, pHist, insiz) if ((insiz) < (pHist)->iInSizLo) { \
    (idx) = 0;								\
  } else if ((insiz) > (pHist)->iInSizHi) {				\
    (idx) = (pHist)->iSpan - 1;						\
  } else {								\
    (idx) = ((insiz) - (pHist)->iInSizLo)/(pHist)->iScalFac;		\
    if ((idx) >= (pHist)->iSpan) {					\
      (idx) = (pHist)->iSpan - 1;					\
    }									\
  }

#ifdef insert_debug
#define CHECK_COUNT_SUM(errcode, pInsHisto, mark) if (!(errcode)) {	\
    uint64_t totnum = 0;						\
    int32_t i;								\
    for (i=0; i<(pInsHisto)->iSpan; i++)				\
      totnum += (pInsHisto)->pCounts[i];				\
    if (totnum != (pInsHisto)->iNum)					\
      (errcode) = ERRCODE_ASSERT;					\
    printf("INSERT_DEBUG: checksum %llu (%llu)\n", (unsigned long long) totnum, (unsigned long long) (pInsHisto)->iNum); \
    if ((errcode)) {							\
      printf("INSERT_DEBUG: checksum failed '%s'\n", (mark));		\
    }									\
  }
#endif
/******************************************************************************
 ****************************** Private Methods *******************************
 ******************************************************************************/

/* static double normal (double x, double mu, double std) */
/* { */
/*   double d = x - mu; */
/*   double v = 2*std*std; */

/*   return exp(-d*d/v)/(std*sqrt(2*M_PI)); */
/* } */

static V_DEF_SORTFUNC(int32_t) /* defines sort_int32_t_quicksort(SAMPLE *) */

static int32_t calcKernelBandWidth(int32_t n, int32_t iqr)
{
  return (n > 0)? (int32_t) (0.9 * pow((double) n, -0.2) * ((double) iqr)/1.34): 0;
}
/******************************************************************************
 ********************** Methods of Private Type KERNEL ************************
 ******************************************************************************/

/******************************************************************************
 *********************** Methods of Private Type SAMPLE ***********************
 ******************************************************************************/

static int getSAMPLEQuartiles(int32_t *med, int32_t *qlo, int32_t *qhi, SAMPLE *pSample)
{ 
  int errcode = sort_int32_t_quicksort(pSample);
  if (errcode) {
    if (med) *med = 0;
    if (qlo) *qlo = 0;
    if (qhi) *qhi = 0;
  } else {
    if (med) *med = V_BASPTR(pSample)[(size_t) (V_LEN(pSample)*.5)];
    if (qlo) *qlo = V_BASPTR(pSample)[(size_t) (V_LEN(pSample)*.25)];
    if (qhi) *qhi = V_BASPTR(pSample)[(size_t) (V_LEN(pSample)*.75)];
  }
  return errcode;
}

/******************************************************************************
 ********************** Public Methods of Type InsSample **********************
 ******************************************************************************/

InsSample *insCreateSample(int blksz)
{
  InsSample *p;
  EMALLOCP0(p);
  if (p) {
    V_CREATE(p->pSample, blksz);
    if (p->pSample == NULL) {
      insDeleteSample(p);
      p = NULL;
    } else {
      p->readival = 0;
    }
  }

  return p;
}

void insDeleteSample(InsSample *pInsSample)
{
  if ((pInsSample)) {
    V_DELETE(pInsSample->pSample);
  }
  free(pInsSample);
}

void insSetSamplingInterval(InsSample *pInsSample, uint64_t nreads, int nrskip)
{
  uint64_t n = nreads/SAMPLE_TARGETSIZ;

  if (n < 1)
    pInsSample->readival = 1;
  else if (n > INT_MAX)
    pInsSample->readival = INT_MAX;
  else
    pInsSample->readival = (int) n;

  if (nrskip < pInsSample->readival && nrskip > 0)
    pInsSample->readival = nrskip;
}

int insAddSample(InsSample *pInsSample, int32_t insertsiz)
{

  V_APPEND(pInsSample->pSample, insertsiz);

  return ERRCODE_SUCCESS;
}

int insIsInSample(const InsSample *pInsSample, uint64_t readno)
{
  return pInsSample == NULL || readno % pInsSample->readival == 0;
}

size_t insGetSample(int32_t **pSample, int *readival, const InsSample *pInsSample)
{
  if (pSample) *pSample = V_BASPTR(pInsSample->pSample);
  if (readival) *readival = pInsSample->readival;
  return V_LEN(pInsSample->pSample);
}
/******************************************************************************
 ********************** Private Methods of Type InsHist ***********************
 ******************************************************************************/

static int findInsHistMax(const InsHist *pHist, int32_t *count_max, 
			  int32_t *range_min, int32_t *range_max)
{
  int32_t i;

  *range_min = *range_max = *count_max = 0;
  for (i=0; i<pHist->iSpan && pHist->pCounts[i] == 0; i++);
  if (i>=pHist->iSpan)
    return ERRCODE_FAILURE;

  *range_min = *range_max = i;
  *count_max = pHist->pCounts[i];

  for (;i<pHist->iSpan; i++)
    if (pHist->pCounts[i] > 0) {
      *range_max = i;
      if (pHist->pCounts[i] > *count_max)
	*count_max = pHist->pCounts[i];
    }
  
  return ERRCODE_SUCCESS;
}

static int smoothGauss(int32_t *targetp, double *Kp, int32_t bw, 
		       const int32_t *sourcep, int32_t n)
{
  int32_t i, j, k, jmax;
  double x;
  int32_t cutoff = KERNEL_CUTOFF_BANDFAC * bw;
  int32_t imax = 2*cutoff + 1;
  const double normfac = sqrt(2*M_PI);
  double tt;

  if (imax > n)
    bw = (n - 1)/(2*KERNEL_CUTOFF_BANDFAC);
  if (bw < KERNEL_MIN_WIDTH)
    bw = KERNEL_MIN_WIDTH;
  cutoff = KERNEL_CUTOFF_BANDFAC * bw;
  imax = 2*cutoff + 1;

  /* calculate Kernel */
  for (i=0; i<imax; i++) {
    x = ((double) (i - cutoff))/bw;
    Kp[i] = exp(-x*x/2)/normfac;
#ifdef insert_debug
    printf("insert_debug:K[%i] %g\n", i, Kp[i]);
#endif
  }
  /* smooth data */
  for (i=0; i<n; i++) {
    targetp[i] = 0;
    if (i > cutoff) {
      j = i - cutoff;
      k = 0;
    } else {
      j = 0;
      k = i;
    }
    if (i + cutoff < n) {
      jmax = i + cutoff;
    } else {
      jmax = n;
    }
    tt = 0.0;
    for (; j<jmax; j++, k++) 
      tt += sourcep[j]*Kp[k];
    targetp[i] = (int32_t) (tt/bw);
  }

  return ERRCODE_SUCCESS;
}
/******************************************************************************
 *********************** Public Methods of Type InsHist ***********************
 ******************************************************************************/

InsHist *insCreateHisto(int32_t iLen)
{
  InsHist *p;

  EMALLOCP0(p);

  if (iLen < 1)
    iLen = DEFAULT_HISTO_WIDTH;
  if (p) {
    ECALLOCP(2*iLen, p->pCounts);
    ECALLOCP(iLen, p->pKernelBuf);
    if ((p->pCounts) && (p->pKernelBuf)) {
      p->pSmoothCounts = p->pCounts + iLen;
      p->iSpan = iLen;
      p->iScalFac = 1;
      p->status = INSFLG_EMPTY;
    } else {
      insDeleteHisto(p);
      p = NULL;
    }
  }

  return p;
}

InsHist *insMakeHistoFromSample(const InsSample *pInsSample)
{
  int32_t med, qlo, qhi, iRange;
  int32_t scf;
  InsHist *pInsHisto = NULL;
  int errcode = getSAMPLEQuartiles(&med, &qlo, &qhi, pInsSample->pSample);
  
  iRange = (qhi - qlo)*IQR_RANGE_FAC*2;
  if (!errcode) {
    int32_t ns = V_LEN(pInsSample->pSample);
    int32_t nbins = (int32_t) (3*sqrt((double) ns));
    if (nbins < HISTO_MIN_BINNUM)
      nbins = HISTO_MIN_BINNUM;
    else if (nbins > HISTO_MAX_BINNUM) 
      nbins = HISTO_MAX_BINNUM;
    scf = iRange/nbins;
    if (scf < 1) {
      nbins = iRange;
      scf = 1;
    } else {
      iRange = scf*nbins;
    }
    pInsHisto = insCreateHisto(nbins);
    if ((pInsHisto)) {
      int32_t i, idx;
      int32_t *p = V_BASPTR(pInsSample->pSample);
      pInsHisto->iScalFac = scf;
      pInsHisto->iInSizLo = med - iRange/2;
      pInsHisto->iInSizHi = pInsHisto->iInSizLo + iRange - 1;
      pInsHisto->median = med;
      pInsHisto->quart_lo = qlo;
      pInsHisto->quart_hi = qhi;

      for (idx=0; idx<pInsHisto->iSpan; idx++)
	pInsHisto->pCounts[idx] = 0;
      for(i=0; i<ns; i++) {
	if (p[i] >= pInsHisto->iInSizLo && p[i] <= pInsHisto->iInSizHi) {
	  CALC_HISTO_IDX(idx, pInsHisto, p[i]);
	  pInsHisto->pCounts[idx]++;
	  pInsHisto->iNum++;
	}
      }
#ifdef insert_debug
      CHECK_COUNT_SUM(errcode, pInsHisto, "insMakeHistoFromSample: before smoothing");
#endif
      errcode = insSmoothHisto(pInsHisto);
#ifdef insert_debug
      CHECK_COUNT_SUM(errcode, pInsHisto, "insMakeHistoFromSample: after smoothing");
#endif
      if (errcode) {
	insDeleteHisto(pInsHisto);
	pInsHisto = NULL;
      }
    }
  }

  return pInsHisto;
}

void insDeleteHisto(InsHist *p)
{
  if (p != NULL) {
    free(p->pCounts);
    free(p->pKernelBuf);
  }
  free(p);
}


/* int insSampleHistoNormal(InsHist *pHist, int32_t mean, int32_t std, int32_t num) */
/* { */
/*   int margin, errcode = ERRCODE_SUCCESS; */
/*   int lo = mean - INTERVAL_Z*std; */
/*   int hi = mean + INTERVAL_Z*std; */
/*   int d = (hi - lo + 1)/pHist->iSpan; */
/*   double *pSample; */
  
/*   if (d >= 1) { */
/*     pHist->iScalFac = d; */
/*   } else { */
/*     pHist->iScalFac = 1; */
/*     margin = pHist->iSpan/d - pHist->iSpan; */
/*     lo -= margin/2; */
/*     hi = lo + pHist->iSpan + margin - 1; */
/*   } */
/*   pHist->iInSizLo = lo; */
/*   pHist->iInSizHi = hi; */
    
/*   ECALLOCP(num, pSample); */
/*   if (pSample == NULL) */
/*     return ERRCODE_NOMEM; */

/*   errcode = rsampleNormal(pSample, num, mean, std*std); */
/*   if (!errcode) { */
/*     int32_t i; */
/*     for (i=0; i<num && !errcode; i++) */
/*       insUpdateHisto(pHist, (int32_t) pSample[i]); */
/*   } */
			  
/*   free(pSample); */
/*   return errcode; */
/* } */

void insSeedHistoNormal(InsHist *pHist, int32_t mean, int32_t std, int32_t num)
{
  int i, margin;
  int lo = mean - INTERVAL_Z*std;
  int hi = mean + INTERVAL_Z*std;
  double x, y, d = ((double) (hi - lo + 1))/pHist->iSpan;
  double var2 = (double) 2*std*std;
  double norm = sqrt(var2*M_PI);
  
  if (d >= 1.0) {
    pHist->iScalFac = (int32_t) d;
  } else {
    pHist->iScalFac = 1;
    margin = (int) (pHist->iSpan/d - pHist->iSpan);
    lo -= margin/2;
    hi = lo + pHist->iSpan + margin - 1;
  }   
  pHist->iInSizLo = lo;
  pHist->iInSizHi = hi;
  
  for (i=0; i<pHist->iSpan; i++) {
    x = pHist->iInSizLo + d*i - mean;
    y = exp(-x*x/var2)/norm;
    pHist->pCounts[i] = (int32_t) (y*num + 0.4999);
  }
  pHist->iNum = num;
  pHist->status = INSFLG_SEEDED;
}

void insUpdateHisto(InsHist *pHist, int32_t insiz)
{
  int idx = (insiz - pHist->iInSizLo)/pHist->iScalFac;

  if (idx >= 0 && idx < pHist->iSpan) {
    pHist->pCounts[idx]++;
    pHist->iNum++;
  }
}

int insSmoothHisto(InsHist *pHist)
{
  int errcode;
  int32_t i, th, kbw;
  int32_t iqr = 0;

  if (pHist->iNum < 2)
    return ERRCODE_FAILURE;

  /* Estimate band width kbw.
   * This would need to be made more robust by combining
   * rank-based stats with stderr (min(std, iqr)) 
   */

  /* determine inter-quartile range (iqr) */
  if (pHist->iSpan > 3) {
    int32_t n=0;
    int q=0;
    uint32_t quart[3];
    th = pHist->iNum/4;
    for (i=0; i<pHist->iSpan && q < 3; i++) {
      n += pHist->pCounts[i];
      if (n > th) {
	quart[q++] = i;
	n -=  pHist->pCounts[i]/2;
	th = pHist->iNum*q/4;
      }
    }
    if (q > 2) 
      iqr = quart[2] - quart[0];
  }

  kbw = calcKernelBandWidth(pHist->iNum, iqr);
  if (kbw < KERNEL_MIN_WIDTH)
    kbw = KERNEL_MIN_WIDTH;

  errcode = smoothGauss(pHist->pSmoothCounts, pHist->pKernelBuf, kbw,
			pHist->pCounts, pHist->iSpan);
  if (!errcode)
    pHist->status = INSFLG_SMOOTHED;

  return errcode;
}

double insGetHistoProb(int32_t insiz, const InsHist *pHist) 
{
  double prob = 0.0;

  if (insiz >= pHist->iInSizLo && insiz <= pHist->iInSizHi &&  pHist->iNum > 0) {
    int32_t idx = (int32_t) (insiz - pHist->iInSizLo)/pHist->iScalFac;
    prob = (pHist->status == INSFLG_SMOOTHED)? 
      (double) pHist->pSmoothCounts[idx]: 
      (double) pHist->pCounts[idx];
    prob /= pHist->iNum;
  }

  return prob;
}

int32_t insGetHistoCount(int32_t *totnum, int32_t insiz, BOOL_t is_smooth, const InsHist *pHist)
{
  int32_t idx, rv = 0;

  if (insiz >= pHist->iInSizLo && insiz <= pHist->iInSizHi) {
    CALC_HISTO_IDX(idx, pHist, insiz);
    rv = ((is_smooth) && pHist->status == INSFLG_SMOOTHED)? pHist->pSmoothCounts[idx]: pHist->pCounts[idx];
  }
  if (totnum) *totnum = pHist->iNum;

  return rv;
}

int32_t insGetHistoCountCumulative(int32_t *totnum, int32_t insiz, BOOL_t is_smooth, const InsHist *pHist)
{
  int32_t i, idx, ccount = 0;
  is_smooth = (BOOL_t) ((is_smooth) && (pHist->status == INSFLG_SMOOTHED));
  if (insiz >= pHist->iInSizLo && insiz <= pHist->iInSizHi) {
    CALC_HISTO_IDX(idx, pHist, insiz);
    for (i=0; i<=idx; i++)
      ccount += (is_smooth)? pHist->pSmoothCounts[i]: pHist->pCounts[i];
  }
  if (totnum) *totnum = pHist->iNum;

  return ccount;
}

uint64_t insGetHistoData(int32_t *lo, int32_t *hi, int32_t *n_bins, const InsHist *pHist)
{
  if (lo) *lo = pHist->iInSizLo;
  if (hi) *hi = pHist->iInSizHi;
  if (n_bins) *n_bins =  pHist->iSpan;

  return pHist->iNum;
}

int32_t insGetHistoQuartiles(int32_t *qlo, int32_t *qhi, const InsHist *pHist)
{
  if (qlo) *qlo = pHist->quart_lo;
  if (qhi) *qhi = pHist->quart_hi;
  return pHist->median;
}

int insPrintHisto(FILE *fp, int linwidth, BOOL_t is_smooth, const InsHist *pHist)
{
  int errcode = ERRCODE_SUCCESS;
  int32_t i, max_count, range_min, range_max;
  double wf;
 
  if ((NULL == pHist) ||
      (errcode = findInsHistMax(pHist, &max_count, &range_min, &range_max))) {
    fprintf(fp, "# Histogram of insert sizes is empty.\n");
    return ERRCODE_FAILURE;
  }
  is_smooth = (BOOL_t) ((is_smooth) && (pHist->status == INSFLG_SMOOTHED));

  wf = ((double) linwidth)/max_count;
  if (wf > 1.0)
    wf = 1.0;
  for (i=range_min; i<=range_max; i++) {
    int32_t j, col_idx = 
      (int32_t) (((is_smooth)? pHist->pSmoothCounts[i]: pHist->pCounts[i])*wf);
    fprintf(fp, "#%5i ", (int) (pHist->iInSizLo + i*pHist->iScalFac));
    for (j=0; j<col_idx; j++)
      if (fputc(HISTO_PRINTCHAR, fp) != HISTO_PRINTCHAR)
	return ERRCODE_WRITEERR;
    fputc(PRINTCHAR_NEWLIN, fp);
  }

  return errcode;
}

int insWriteHisto(FILE *fp, BOOL_t is_smooth, const InsHist *pHist)
{
  int32_t i;
  int32_t *cp;
  uint64_t totnum = 0;
  
  if (pHist == NULL)
    return ERRCODE_FAILURE;

  is_smooth = (BOOL_t) ((is_smooth) && (pHist->status == INSFLG_SMOOTHED));
  cp = (is_smooth)? pHist->pSmoothCounts: pHist->pCounts;
  for (i=0; i<pHist->iSpan; i++) 
    totnum += cp[i];
  if (!is_smooth && totnum != pHist->iNum)
    return ERRCODE_ASSERT;

  fprintf(fp, "%s", IOFIL_HEADER);
  fprintf(fp, "%s\n", IOFIL_KEY_START);
  fprintf(fp, IOFIL_KEY_FORMAT, pHist->iSpan, pHist->iScalFac, 
	  pHist->iInSizLo, pHist->iInSizHi, (unsigned long long) totnum,
	  pHist->quart_lo, pHist->median, pHist->quart_hi);
  for (i=0; i<pHist->iSpan; i++) {
    fprintf(fp, IOFIL_FORMAT, pHist->iInSizLo + i*pHist->iScalFac, cp[i]);
  }
  fprintf(fp, "%s\n", IOFIL_KEY_END);

  return ERRCODE_SUCCESS;
}

InsHist *insReadHisto(int *errcode, const char *filnam)
{
  char linbufp[IOFIL_LINBUFSIZ];
  BOOL_t has_start, has_end;
  int32_t isiz, count, linctr, qlo, qhi, med;
  int binnum, scalfac, insizlo, insizhi;
  unsigned long long totnum;
  size_t slen;
  InsHist *ihp = NULL;
  FILE *fp = EFOPEN(filnam, "r");

  if (fp == NULL) {
    *errcode = ERRCODE_FILEIO;
    return NULL;
  }
  *errcode = ERRCODE_SUCCESS;

  has_start = has_end = 0;
  slen = strlen(IOFIL_KEY_START);
  while (fgets(linbufp, IOFIL_LINBUFSIZ-1, fp) != NULL) {
    if (!strncmp(IOFIL_KEY_START, linbufp, slen)) {
      has_start = 1;
      break;
    }
  }
  if (!has_start || 
      fscanf(fp, IOFIL_KEY_FORMAT, &binnum, &scalfac, 
	     &insizlo, &insizhi, &totnum, &qlo, &med, &qhi) != 8) {
    *errcode = ERRCODE_FILEFORM;
    EFCLOSE(fp);
    return NULL;
  }
  ihp = insCreateHisto(binnum);
  if (ihp == NULL) {
    *errcode = ERRCODE_NOMEM;
    EFCLOSE(fp);
  }
  ihp->iScalFac = scalfac;
  ihp->iInSizLo = insizlo;
  ihp->iInSizHi = insizhi;
  ihp->median = med;
  ihp->quart_lo = qlo;
  ihp->quart_hi = qhi;
  ihp->iNum = 0;

  linctr = 0;
  slen = strlen(IOFIL_KEY_END);
  while (fgets(linbufp, IOFIL_LINBUFSIZ-1, fp) != NULL) {
    if (!strncmp(IOFIL_KEY_END, linbufp, slen)) {
      has_end = 1;
      break;
    }
    if (sscanf(linbufp, IOFIL_FORMAT, &isiz, &count) != 2)
      break;
    if (isiz != insizlo + linctr*scalfac)
      break;
    if (linctr >= binnum)
      break;
    ihp->pCounts[linctr++] = count;
    ihp->iNum += count;
  }

  if (!has_end || ihp->iNum != totnum) {
    *errcode = ERRCODE_FILEFORM;
    EFCLOSE(fp);
  } else {
    *errcode = EFCLOSE(fp);
  }

  if (!(*errcode)) 
      *errcode = insSmoothHisto(ihp);

  return ihp;
}
