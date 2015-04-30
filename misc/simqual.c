/**< Simulate base qualities for a set of reads. 
 * Original reads are modified according to base qualities.
 *
 * Statistics on base qualities must have been collected using
 * basqcol.c
 */

/****************************************************************************
 ****************************************************************************
 *                                                                          *
 *  Copyright (C) 2010-2014 Genome Research Ltd.                            * 
 *                                                                          *        
 *  Author: Hannes Ponstingl (hp3@sanger.ac.uk)                             *
 *                                                                          *
 *  This file is part of SMALT                                              *
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
#include <limits.h>
#include <stdint.h>
#include <ctype.h>
#include <string.h>
#include <math.h>

#include "elib.h"
#include "randef.h"
#include "sequence.h"
#include "diffstr.h"
#include "basqual.h"

enum {
  NAMBUF_MAXLEN = 4096,
  NAMSIZ_PER_MUTREC = 20, /**< allow this many characters for 1 mutation in string */
  PHRED_SCALE = 10,
  PHRED_LOGBASE = 10,
};

static int mutateNtSeq(SeqFastq *sqp, 
		       int *mutreclen, char *mutchar, int *mutpos, 
		       const SeqCodec *codecp)
     /**< Mutate nuclotide sequence randomly with rates given by the base qualities. 
      * \param sqp Sequence to be mutated (with base qualities)
      * \param mutreclen Returns the number of bases mutated (can be NULL).
      * \param mutchar array returning the original bases 
      *        of each mutation (must have been pre-allocated, can be NULL).
      * \param mutpos array returning the position 
      *        of each mutation (must have been pre_allocated, can be NULL).
      */
{
  char *qualp, *basp, cod;
  short alphlen, encodsiz;
  int modi, stdnt_idx, bq;
  SEQLEN_t i, rlen;
  double errprob, randunit;
  const char *alphabetp = seqCodecGetAlphabet(codecp, &alphlen);
  const unsigned char *encoderp = seqCodecGetEncoder(codecp, &encodsiz);
  const double phredexp = -1*log(PHRED_LOGBASE)/PHRED_SCALE;
  int n_mutrecord = 0;

  if (SEQCOD_STDNT_MASK >= alphlen)
    return ERRCODE_ASSERT;
  if (mutreclen) *mutreclen = 0;
  basp = seqFastqGetSequence(sqp, &rlen, &cod);
  if (cod != SEQCOD_ASCII)
    return ERRCODE_SEQCODE;
  qualp = seqFastqGetQualityFactors(sqp, &rlen, NULL);

  if (rlen > INT_MAX)
    return ERRCODE_OVERFLOW;

  for (i=0; i<rlen; i++) {
    bq = ((int) qualp[i]) - SEQCOD_QVAL_OFFS;
    if (bq < 0) 
      return ERRCODE_ASSERT;
    if (bq == 0) 
      continue;
    errprob = exp(phredexp*bq);
    randunit = RANDRAW_UNIFORM_1();
    if (randunit > errprob)
      continue;
    modi = (int) (randunit*SEQCOD_STDNT_MASK/errprob);
    if (modi + 1 > SEQCOD_STDNT_MASK)
      modi = SEQCOD_STDNT_MASK - 1;
    else if (modi < 0)
      modi = 0;
    
    if (basp[i] >= encodsiz)
      return ERRCODE_ASSERT;

    if (mutchar) mutchar[n_mutrecord] = basp[i];
    if (mutpos) mutpos[n_mutrecord] = (int) i;
    n_mutrecord++;
    stdnt_idx = ((encoderp[(int) basp[i]] & SEQCOD_STDNT_MASK) + modi) % (SEQCOD_STDNT_MASK + 1);
    basp[i] = alphabetp[stdnt_idx];
  }
 
  if (mutreclen) *mutreclen = n_mutrecord;
  return ERRCODE_SUCCESS;
}

int main(int argc, char *argv[])
{
  int errcode = EXIT_SUCCESS;
  char nambuf[NAMBUF_MAXLEN];
  char *qualp, *mutbasp;
  char *filnam_basq, *filnam_fastq_in, *filnam_fastq_out;
  int is_modify;
  uint8_t nqual, qmin;
  int rand_seed;
  SEQLEN_t rlen, rlen_simulated;
  int *mutposp, mutrecnum = 0;
  int64_t tot_mutnum = 0, tot_basnum = 0;
  SeqIO *sifp, *sofp;
  SeqFastq *sqbufp;
  SeqCodec *codecp = NULL;
  BasQualFreq_t *bqfp;
  DiffStr *dfsp = NULL;
  ErrMsg *errmsgp=0;
  
  ERRMSG_CREATE(errmsgp);

  if (argc < 6) {
    fprintf(stderr, "usage: %s <base quality file (in)> <seed> "\
	    "<modify bases [y/n]> <FASTQ file (in)> <FASTQ file (out)>\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  filnam_basq = argv[1];
  rand_seed = atoi(argv[2]);
  is_modify = (toupper(argv[3][0]) == 'Y')? 1: 0;
  filnam_fastq_in = argv[4];
  filnam_fastq_out = argv[5];

  RANSEED(rand_seed);

  if (!(sqbufp = seqFastqCreate(0, SEQTYP_FASTQ)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  if ((is_modify)) {
    if (!(codecp = seqCodecCreate()))
      ERRMSGNO(errmsgp, ERRCODE_NOMEM);

    if (!(dfsp = diffStrCreate(0)))
      ERRMSGNO(errmsgp, ERRCODE_NOMEM);
  }

  printf("\nReading base qualities from file %s...\n", filnam_basq);
  bqfp = basQualFreqRead(&errcode, filnam_basq);
  if ((errcode))
    ERRMSGNO(errmsgp, errcode);

  printf("\nCalculate sums ...\n");
  if ((errcode = basQualFreqSum(bqfp)))
    ERRMSGNO(errmsgp, errcode);
  
  rlen_simulated = basQualFreqGetData(bqfp, NULL, &nqual, &qmin, NULL, NULL);
  printf("Simulated read length = %u\n", (unsigned int) rlen_simulated);
  printf("Minimum quality value = %i\n", (int) qmin);
  printf("Maximum quality value = %i\n", (int) qmin + nqual);

  if (rlen_simulated > INT_MAX)
    ERRMSG(errmsgp, "reads too long for simulation", ERRCODE_OVERFLOW);

  /* allocate helper arrays */
  ECALLOCP(rlen_simulated, mutbasp);
  if (mutbasp == NULL)
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  ECALLOCP(rlen_simulated, mutposp);
  if (mutposp == NULL)
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);
  
  printf("Open FASTQ file for input ...\n");
  sifp = seqIOopen(&errcode, filnam_fastq_in, SEQIO_READ, 0);
  if (errcode) 
    ERRMSGNO(errmsgp, errcode);

  printf("Open FASTQ file for output ...\n");
  sofp = seqIOopen(&errcode, filnam_fastq_out, SEQIO_WRITE_FASTQ, 0);
  if (errcode) 
    ERRMSGNO(errmsgp, errcode);

  while(!seqIOstatus(sifp) && !seqIOstatus(sofp)) {
    //printf("Reading sequence %i ...\n", sctr);
    if ((errcode = seqFastqRead(sqbufp, sifp)))
      break;
    qualp = seqFastqGetQualityFactors(sqbufp, &rlen, NULL);

    if (rlen > rlen_simulated)
      ERRMSG(errmsgp, "read length greater than simulated", ERRCODE_ASSERT);
      
    if ((errcode = basQualFreqSimulate(qualp, rlen, bqfp)))
      ERRMSGNO(errmsgp, errcode);
    
    tot_basnum += rlen;

    if ((is_modify)) {
      int dlen, buflen, nchar=0;
      mutrecnum = 0;
      if ((errcode = mutateNtSeq(sqbufp, &mutrecnum, mutbasp, mutposp, codecp)))
	ERRMSGNO(errmsgp, errcode);
      tot_mutnum += mutrecnum;

      if ((errcode = diffStrGenerateFromMismatches(&dlen, NULL, mutposp, mutrecnum, (int) rlen)))
	ERRMSGNO(errmsgp, errcode);

      if (dlen >= dfsp->n_alloc &&
	  (errcode = diffStrRealloc(dfsp, dlen)))
	ERRMSGNO(errmsgp, errcode);

      diffStrGenerateFromMismatches(NULL, dfsp->dstrp, mutposp, mutrecnum, (int) rlen);

      strncpy(nambuf, seqFastqGetSeqName(sqbufp), NAMBUF_MAXLEN);
      nambuf[NAMBUF_MAXLEN-1] = '\0';
      buflen = (int) strlen(nambuf);
      if (buflen < 1 || buflen*NAMSIZ_PER_MUTREC > NAMBUF_MAXLEN-1)
	ERRMSG(errmsgp, "Name buffer overflow", ERRCODE_ASSERT);

      nambuf[buflen++] = ' ';
      nambuf[buflen] = '\0';
      if ((errcode = diffStrPrintfStr(nambuf + buflen, &nchar, dfsp->dstrp, DIFFSTRFORM_PLAIN, 0, 0, 0)))
	ERRMSGNO(errmsgp, errcode);
      nambuf[NAMBUF_MAXLEN-1] = '\0';
      
      if ((errcode = seqFastqSetAscii(sqbufp, nambuf, NULL, NULL, NULL)))
	ERRMSGNO(errmsgp, errcode);
    }
    
    if ((errcode = seqFastqWrite(sofp, sqbufp, 0)))
      ERRMSGNO(errmsgp, errcode);
  }
  if (seqIOstatus(sofp))
    ERRMSGNO(errmsgp, seqIOstatus(sofp));
  seqIOclose(sofp);

  if (seqIOstatus(sifp) && 
      seqIOstatus(sifp) != ERRCODE_EOF) 
    ERRMSGNO(errmsgp, seqIOstatus(sifp));
  seqIOclose(sifp);

  printf("# simqual: introduced %lli mutations in %lli bases (%Lg%%)\n", 
	 (long long int) tot_mutnum, (long long int) tot_basnum,
	 (((long double) tot_mutnum)/tot_basnum)*100);

  free(mutbasp);
  free(mutposp);

  basQualFreqDelete(bqfp);
  diffStrDelete(dfsp);
  seqCodecDelete(codecp);
  seqFastqDelete(sqbufp);

  ERRMSG_END(errmsgp);
  exit(errcode);
}
