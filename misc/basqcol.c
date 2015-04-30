/** Collect statistics on base qualities
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
#include <math.h>

#include "elib.h"
#include "sequence.h"
#include "basqual.h"

typedef long long unsigned int LLUINT;

int main(int argc, char *argv[])
{
  int argi, iv, argi_start, errcode = EXIT_SUCCESS;
  uint8_t maxq, minq;
  uint8_t maxq_tot = 0, minq_tot = UINT8_MAX, minbasq;
  uint32_t maxlen, minlen;
  uint32_t maxlen_tot = 0, minlen_tot = UINT32_MAX;
  uint64_t nreads, nreads_tot=0;
  char *infilnam, *oufilnam;
  SeqIO *sifp;
  SeqFastq *sqbufp;
  BasQualFreq_t *bqfp;
  ErrMsg *errmsgp=0;
  
  ERRMSG_CREATE(errmsgp);

  if (argc < 4) { 
    fprintf(stderr, "usage: %s <base quality file (output)> <min basqual> <FASTQ file 1> "\
	    "[<FASTQ file 2> <FASTQ file 3> ...]\n", argv[0]);
    
    exit(EXIT_FAILURE);
  }

  oufilnam = argv[1];
  iv = atoi(argv[2]); /* Don't count when base quality below this value */
  if (iv < 0 || iv + SEQCOD_QVAL_OFFS > UINT8_MAX) {
    fprintf(stderr, "base quality threshold must be a number between 0 and %i\n", 
	    UINT8_MAX - SEQCOD_QVAL_OFFS);
    exit(EXIT_FAILURE);
  }
  minbasq = (uint8_t) iv;
  argi_start = 3;

  if (!(sqbufp = seqFastqCreate(0, SEQTYP_FASTQ)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  printf("# Find out Extrema ...\n");
  for (argi=argi_start; argi<argc; argi++) {
    infilnam = argv[argi];
    printf("Processing file %s ...\n", infilnam);
    sifp = seqIOopen(&errcode, infilnam, SEQIO_READ, 0);
    if (errcode) 
      ERRMSGNO(errmsgp, errcode);

    printf("# Find out maximum read length ...\n");
    errcode = basQualFindExtrema(sifp, sqbufp, &nreads, &maxlen, &minlen, &maxq, &minq);
    printf("# Number of reads: %llu\n", (LLUINT) nreads);
    printf("# Maximum read length: %u\n", maxlen);
    printf("# Minimum read length: %u\n", minlen);
    printf("# Maximum quality: %hi\n", maxq);
    printf("# Minimum quality: %hi\n", minq);

    nreads_tot += nreads;
    if (maxlen > maxlen_tot)
      maxlen_tot = maxlen;
    if (minlen < minlen_tot)
      minlen_tot = minlen;
    if (maxq > maxq_tot)
      maxq_tot = maxq;
    if (minq < minq_tot)
      minq_tot = minq;
 
    if (seqIOstatus(sifp) && 
	seqIOstatus(sifp) != ERRCODE_EOF) 
      ERRMSGNO(errmsgp, seqIOstatus(sifp));
    seqIOclose(sifp);
  }
  
  printf("##########################\n");
  printf("# Total number of reads: %llu\n", (LLUINT) nreads_tot);
  printf("# Overall maximum read length: %u\n", maxlen_tot);
  printf("# Overall minimum read length: %u\n", minlen_tot);
  printf("# Overall maximum quality: %hi\n", maxq_tot);
  printf("# Overall minimum quality: %hi\n", minq_tot);
  if (minq_tot < minbasq) {
    printf("# Overall minimum quality threshold applied: %hi\n", minbasq);
    minq_tot = minbasq;
  }
  if (!(bqfp = basQualFreqCreate(minq_tot, (uint8_t) (maxq_tot - minq_tot + 1), maxlen_tot)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);
 
  printf("\nGetting base quality counts ...\n");

  for (argi=argi_start; argi<argc; argi++) {
    infilnam = argv[argi];
    printf("Processing file %s ...\n", infilnam);
    sifp = seqIOopen(&errcode, infilnam, SEQIO_READ, 0);
    if (errcode) 
      ERRMSGNO(errmsgp, errcode);

    if((errcode = basQualFreqFromFastq(bqfp, sqbufp, sifp)))
      ERRMSGNO(errmsgp, errcode);

    if (seqIOstatus(sifp) && 
	seqIOstatus(sifp) != ERRCODE_EOF) 
      ERRMSGNO(errmsgp, seqIOstatus(sifp));
    seqIOclose(sifp);
  }

  printf("\nWriting base qualities to file ...\n");
  if((errcode = basQualFreqWrite(bqfp, oufilnam)))
    ERRMSGNO(errmsgp, errcode);
  
  basQualFreqPrint(stdout, bqfp);

  basQualFreqDelete(bqfp);
  seqFastqDelete(sqbufp);
  
  ERRMSG_END(errmsgp);
  exit(errcode);
}
