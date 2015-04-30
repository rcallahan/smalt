/** Analyse base quality scores of a set of sequencing reads.
 */

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
#include <stdlib.h>
#include <stdint.h>

#include "elib.h"
#include "filio.h"
#include "basqual.h"

enum {
  BASQUALFIL_VERSION = 1,
  BASQUALFIL_HEADSIZ = 3,
};

enum BASQUALFREQ_FLAGS { /**< values for BasQualFreq_t.status flag */
  BASQUALFLG_EMPTY = 0,  /**< Structure is empty */
  BASQUALFLG_COUNTS = 1, /**< Type has counts (but sums have not been updated) */
  BASQUALFLG_SUMS = 3,   /**< Type has counts and up-to-date sums */
};

static const char BASQUALFIL_NAMEXT[] = "smq"; /* file name extension for base quality file */
static const char BASQUALFIL_WRITERRMSG[] = "when writing file with base quality statistics";

struct _BasQualFreq_t {
  uint8_t status;/**< One of BASQUALFREQ_FLAGS */
  uint32_t *q0p; /**< base quality counts for state at beginning */
  uint64_t q0s;  /**< Sum_i q0p[i] */
  uint8_t nq;    /**< number of quality values (maximum quality value) */
  uint8_t qmin;  /**< minimum quality value based on SEQCOD_QVAL_OFFS */
  uint32_t rlen; /**< (maximum) read length */
  
  uint32_t *qtp; /**< qtp[nq][nq][readlen-1] */
  uint64_t *qsp; /**< size qsp[nq][readlen-1], qsp[i][r] = Sum_j qtp[i][j][r] */
};

/******************************************************************************
 *********************************** Macros ***********************************
 ******************************************************************************/

#define DRAW_UNIFORM_1() ((long double) rand()/(((long double) RAND_MAX) + 1))

/******************************************************************************
 ******************** Private Methods of Type BasQualFreq *********************
 ******************************************************************************/

static int allocBasQualFreqSums(BasQualFreq_t *p)
{
  p->q0s = 0;
  if (p->qsp != NULL) {
    free(p->qsp);
    p->qsp = NULL;
  }
  ECALLOCP(p->nq * (p->rlen-1), p->qsp);
  return (p->qsp == NULL)? ERRCODE_NOMEM: ERRCODE_SUCCESS;
}

static int calcSums(BasQualFreq_t *p)
{
  int errcode, i, j;
  int nq = p->nq;
  uint32_t bs, bt, max_bs;
  uint64_t sum;

  if (p->status != BASQUALFLG_COUNTS &&       
      p->status != BASQUALFLG_SUMS)
    return ERRCODE_ASSERT;

  if (p->qsp == NULL && 
      (errcode = allocBasQualFreqSums(p)))
    return errcode;
  
  sum = 0LL;
  for (i=0; i<nq; i++)
    sum += p->q0p[i];
  p->q0s = sum;

  max_bs = nq*(p->rlen-1);
  for (bs=0; bs < max_bs; bs++) {
    sum = 0LL;
    bt = bs*nq;
    for (j=0; j<nq; j++) {
      sum += p->qtp[bt+j];
    }
    p->qsp[bs] = sum;
  }
#ifdef basqual_debug
  for (r=1; r<p->rlen; r++) {
    for (i=0; i<nq; i++) {
      sum = 0LL;
      bs = (r-1)*nq + i;
      bt = bs*nq;
      for (j=0; j<nq; j++) 
	sum += p->qtp[bt+j];
      if (sum != p->qsp[bs]) {
	printf("basqual_debug::calcSums: inconsistent sums (r,i) = (%u,%i)\n", r, i);
	return ERRCODE_ASSERT;
      }
    }
  }
#endif

  p->status = BASQUALFLG_SUMS;

  return ERRCODE_SUCCESS;
}

/******************************************************************************
 ********************* Public Methods of Type BasQualFreq *********************
 ******************************************************************************/

BasQualFreq_t *basQualFreqCreate(uint8_t qmin, uint8_t nq, uint32_t readlen)
{
  BasQualFreq_t *p;

  EMALLOCP0(p);
  if (p != NULL) {
    ECALLOCP(nq + (nq*nq*(readlen-1)), p->q0p);
    if (p->q0p == NULL) {
      basQualFreqDelete(p);
      p = NULL;
    } else {
      p->qtp = p->q0p + nq;
      p->qmin = qmin;
      p->nq = nq;
      p->rlen = readlen;
      p->status = BASQUALFLG_EMPTY;
    }
  }

  return p;
}

void basQualFreqDelete(BasQualFreq_t *p)
{
  if (p != NULL) {
    free(p->q0p);
    free(p->qsp);
    free(p);
  }
}

uint32_t *basQualFreqFld(int *errcode, uint32_t r,  uint8_t i, uint8_t j, BasQualFreq_t *bqfp)
{
  *errcode = ERRCODE_SUCCESS;
  uint32_t *p = NULL;

  if (i > bqfp->nq) {
    *errcode =  ERRCODE_ARGRANGE;
  } else {
    if (r == 0)
      p = bqfp->q0p + i;
    else if (j > bqfp->nq || r > bqfp->rlen)
      *errcode = ERRCODE_ARGRANGE;
    else 
      p = bqfp->q0p + (i*bqfp->nq + j)*(bqfp->rlen-1) + r - 1;
  }
  return p;
}

uint32_t basQualFreqGetData(const BasQualFreq_t *p, const uint32_t **q0p, uint8_t *nq, uint8_t *qmin, 
			    const uint32_t **qtp, const uint64_t **qsp)
{
  if (q0p) *q0p = p->q0p;
  if (nq) *nq = p->nq;
  if (qmin) *qmin = p->qmin;
  if (qtp) *qtp = p->qtp;
  if (qsp) *qsp = p->qsp;
  return p->rlen;
}

int basQualFreqWrite(const BasQualFreq_t *p, const char *filnam)
{
  int errcode = ERRCODE_SUCCESS;
  uint32_t header[BASQUALFIL_HEADSIZ];
  size_t totsiz;
  FILE *fp;
  
  if (p->status == BASQUALFLG_EMPTY)
    return ERRCODE_ASSERT;

  header[0] = p->rlen;
  header[1] = p->nq;
  header[2] = p->qmin;
  totsiz = p->nq + p->nq*p->nq*(p->rlen-1);
  if (totsiz > UINT32_MAX)
    return ERRCODE_OVERFLOW;
  fp = filioOpenForWriting(&errcode, totsiz, FILIOTYP_BASQUAL,
			   BASQUALFIL_VERSION, BASQUALFIL_HEADSIZ,
			   header, filnam, BASQUALFIL_NAMEXT);
  if (errcode)
    return errcode;

  fwrite(p->q0p, sizeof(uint32_t), totsiz, fp);

  if (ferror(fp)) {
    perror(BASQUALFIL_WRITERRMSG);
    errcode = ERRCODE_WRITEERR;
  }
  
  fclose(fp);
  return errcode;
}

BasQualFreq_t *basQualFreqRead(int *errcode, const char *filnam)
{
  uint8_t is_endianid, typ;
  uint32_t totsiz, version;
  size_t siz;
  uint32_t headsiz = BASQUALFIL_HEADSIZ;
  uint32_t header[BASQUALFIL_HEADSIZ];
  BasQualFreq_t *p = NULL;
  FILE *fp = filioOpenForReading(errcode, &is_endianid, &totsiz,
				 &typ, &version,
				 &headsiz, header, filnam, BASQUALFIL_NAMEXT);
  if (*errcode)
    return NULL;

  if (typ != FILIOTYP_BASQUAL)
    *errcode = ERRCODE_FILTYP;
  else if (version != BASQUALFIL_VERSION)
    *errcode = ERRCODE_FILVERSION;
  else if (headsiz > BASQUALFIL_HEADSIZ || 
	   header[0] < 1 ||
	   header[1] < 1 ||
	   header[1] > UINT8_MAX ||
	   header[2] + SEQCOD_QVAL_OFFS > UINT8_MAX)
    *errcode = ERRCODE_FILEFORM;
  if (*errcode) {
    EFCLOSE(fp);
    return NULL;
  }
  
  if (!(p = basQualFreqCreate((uint8_t) header[2],(uint8_t) header[1], header[0]))) {
    *errcode = ERRCODE_NOMEM;
    EFCLOSE(fp);
    return NULL;
  }
  siz = header[1] + (header[1]*header[1]*(header[0] - 1));
  if (siz != totsiz) {
    *errcode =  ERRCODE_ASSERT;
    return p;
  }

  if (fread(p->q0p, sizeof(uint32_t), totsiz, fp) != totsiz)
    *errcode = ERRCODE_FILEFORM;
  else if((ferror(fp)))
    *errcode = ERRCODE_READERR;
  else if (!(is_endianid)) {
    filioSwapEndian(p->q0p, totsiz);
  }

  EFCLOSE(fp);

  p->status = BASQUALFLG_COUNTS;
  return p;
}

int basQualFreqFromFastq(BasQualFreq_t *bqfp, SeqFastq *sqbufp, SeqIO *sifp)
{
  int errcode = ERRCODE_SUCCESS;
  const char *basq;
  int b, b_prev;
  uint32_t r, readlen;
  uint8_t nq;
  uint64_t rctr = 0LL;

  while(!seqIOstatus(sifp)) {
    //printf("Reading sequence %i ...\n", sctr);
    if ((errcode = seqFastqRead(sqbufp, sifp)))
      break;
    basq = seqFastqGetConstQualityFactors(sqbufp, &readlen, NULL);
    
    if (readlen < 1 || readlen > bqfp->rlen)
      return ERRCODE_ASSERT;
    
    b = basq[0] - SEQCOD_QVAL_OFFS - bqfp->qmin;
    if (b < 0)
      continue; /* b < 0 can occur because of applied threshold, but
		 * read won't count */
    if (b >= bqfp->nq) 
      return ERRCODE_ASSERT;
 
    bqfp->q0p[b]++;
    nq = bqfp->nq;
    for (r=1; r<readlen; r++) {
      b_prev = b;
      b = basq[r] - SEQCOD_QVAL_OFFS - bqfp->qmin;
      if (b < 0 || b_prev < 0)
	continue; /* below applied threshold */
      if (b >= bqfp->nq)
	return ERRCODE_ASSERT;
      bqfp->qtp[((r-1)*nq + b_prev)*nq + b]++;
    }
    rctr++;
  }

  if (seqIOstatus(sifp) && 
      seqIOstatus(sifp) != ERRCODE_EOF) 
    errcode = seqIOstatus(sifp);
  
  if (!errcode) 
    bqfp->status = BASQUALFLG_COUNTS;
  return errcode;
}

void basQualFreqPrint(FILE *fp, const BasQualFreq_t *bqfp)
{
  int i, j, nq = (int) bqfp->nq;
  int nq2 = bqfp->nq*bqfp->nq;
  uint32_t r, count;
  uint64_t bas_r, bas_i;

  fprintf(fp, "Base quality | counts\n");
  for (i=0; i<nq; i++)
    fprintf(fp, "%3i %6u\n", (int) bqfp->qmin + i, bqfp->q0p[i]); 

  fprintf(fp, "Transition counts\n");
  fprintf(fp, "Read position | quality | quality at next position| count\n");
  for (r=1; r<bqfp->rlen; r++) {
    bas_r = (r-1)*nq2;
    for (i=0; i<nq; i++) {
      bas_i = bas_r + nq*i;
      for (j=0; j<nq; j++) {
	count = bqfp->qtp[bas_i + j];
	if (count > 0) 
	  fprintf(fp, "%4u %3i %3i %8u\n", r, i + bqfp->qmin, j + bqfp->qmin, count);
      }
    }
  }
}

int basQualFreqSum(BasQualFreq_t *bqfp)
{
  return calcSums(bqfp);
}

int basQualFreqSimulate(char *qualp, uint32_t len, const BasQualFreq_t *bqfp)
{
  int i, j, nq = (int) bqfp->nq;
  unsigned char qbas = (unsigned char) (bqfp->qmin + SEQCOD_QVAL_OFFS);
  uint32_t r, bs, bt;
  uint64_t sum, pivot;

  if (len > bqfp->rlen || bqfp->status != BASQUALFLG_SUMS)
    return ERRCODE_ASSERT;
  
  /* sample start base quality from empirical distribution */
  sum=0;
  pivot = (uint64_t) (bqfp->q0s*DRAW_UNIFORM_1());
  for (i=0; i<bqfp->nq; i++) {
    sum += bqfp->q0p[i];
    if (sum > pivot)
      break;
  }
  qualp[0] = (char) (i + qbas);
  
  /* sample the subsequent base qualities from empirical transition
   * distributions */
  for (r=1; r<len; r++) {
    bs = (r-1)*nq + i;
    bt = bs*nq;
#ifdef basqual_debug
    sum = 0;
    for (j=0; j<nq; j++) 
      sum += bqfp->qtp[bt+j];
    if (sum != bqfp->qsp[bs]) {
      printf("basqual_debug::basQualFreqSimulate: inconsistent sums (r,i) = (%u,%i)\n", r, i);
      return ERRCODE_ASSERT;
    }
#endif
    
    if (bqfp->qsp[bs] > 0) {
      pivot = (uint64_t) (bqfp->qsp[bs]*DRAW_UNIFORM_1());
      sum = 0;
      for (j=0; j<nq; j++) {
	sum += bqfp->qtp[bt+j];
	if (sum > pivot)
	  break;
      }
      if (j >= nq) j = nq - 1;
      qualp[r] = (char) (j + qbas);
      i = j;
    } else { /* if there are no counts for this base quality choose the same base quality */
      qualp[r] = (char) (i + qbas);
    }
  }
  qualp[r] = '\0';
  
  return ERRCODE_SUCCESS;
}
/****************************************************************************
 ************************ Other Public Methods ******************************
 ****************************************************************************/

int basQualFindExtrema(SeqIO *sifp, SeqFastq *sqbufp, uint64_t *nreads, 
		       uint32_t *maxlen, uint32_t *minlen,
		       uint8_t *maxq, uint8_t *minq)
{
  int errcode = ERRCODE_SUCCESS;
  const char *basq;
  uint8_t q_min, q_max;
  uint32_t i, readlen, len_max, len_min;
  uint64_t rctr = 0LL;

  if (nreads) *nreads = 0;
  if (maxlen) *maxlen = 0;
  if (minlen) *minlen = 0;
  if (maxq) *maxq = 0;
  if (minq) *minq = 0;

  len_max = 0;
  len_min = UINT32_MAX;
  q_max = 0;
  q_min = UINT8_MAX;

  while(!seqIOstatus(sifp)) {
    //printf("Reading sequence %i ...\n", sctr);
    if ((errcode = seqFastqRead(sqbufp, sifp)))
      break;
    basq = seqFastqGetConstQualityFactors(sqbufp, &readlen, NULL);
    if (readlen > len_max)
      len_max = readlen;
    else if (readlen < len_min)
      len_min = readlen;
    for (i=0; i<readlen; i++) {
      if (basq[i] > q_max)
	q_max = (uint8_t) basq[i];
      else if (basq[i] < q_min)
	q_min = (uint8_t) basq[i];
    }
    rctr++;
  }
  if (q_max >= q_min) {
    if (maxq) *maxq = (uint8_t) (q_max - SEQCOD_QVAL_OFFS);
    if (minq) *minq = (uint8_t) (q_min - SEQCOD_QVAL_OFFS);
  }
  if (len_max >= len_min) {
    if (maxlen) *maxlen = len_max;
    if (minlen) *minlen = len_min;
  }
  if (nreads) *nreads = rctr;

  if (seqIOstatus(sifp) && 
      seqIOstatus(sifp) != ERRCODE_EOF) 
    errcode = seqIOstatus(sifp);

  return errcode;
} 

