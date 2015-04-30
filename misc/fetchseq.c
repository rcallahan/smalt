/** Fetch a segment from type sequence.c:SeqSet.
 *
 * Read a binary file with the set of sequences.
 * Print a segment from that file identified by
 * Sequence name and offset.
 */

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>

#include "../src/elib.h"
#include "../src/sequence.h"

static int cmpStrNonBlank(const char *ap, const char *bp)
/**< compare 2 strings up to white space */
{
  int i;
  for(i=0; i<INT_MAX && ap[i] == bp[i] && ap[i] && !isspace(ap[i]); i++);
  if (ap[i] == bp[i] ||
      (isspace(ap[i]) && !bp[i]) ||
      (isspace(bp[i]) && !ap[i]))
    return 0;
  return (ap[i] < bp[i])? -1: 1;
}

int main (int argc, char *argv[])
{
  int errcode = ERRCODE_SUCCESS;
  SEQNUM_t s, nseq;
  char *infilnam, *seqnam;
  const char *cp, *snam;
  SEQLEN_t qlen, segoffs, seglen;
  SETSIZ_t soffs;
  unsigned long long sso;
  SeqFastq *sqp;
  SeqSet *ssp;
  SeqCodec *codecp;
  ErrMsg *errmsgp = 0;

  ERRMSG_CREATE(errmsgp);

  if (argc != 5) {
    printf("usage: %s <Sequence set binary file (no ext.)> <seq id> <start pos (from 1)> <len>\n", 
	   argv[0]);
    exit(1);
  }

  infilnam = argv[1];
  seqnam = argv[2];
  segoffs = (unsigned int) strtoul(argv[3], NULL, 0);
  if (segoffs < 1) segoffs = 1;
  seglen = (unsigned int) strtoul(argv[4], NULL, 0);

  if (!(codecp = seqCodecCreate()))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  if (!(sqp = seqFastqCreate(0,SEQTYP_UNKNOWN)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  /* re-read the binary file */
  printf("Reading sequence set %s ...\n", infilnam);

  ssp = seqSetReadBinFil(&errcode, infilnam);
  if ((errcode))
    ERRMSGNO(errmsgp, errcode);

  nseq = seqSetGetOffsets(ssp, NULL);
  for (s=0; s<nseq; s++) {
    seqSetGetSeqDatByIndex(&soffs, &snam, s, ssp);
    if (!cmpStrNonBlank(snam, seqnam)) {
      if ((errcode = seqSetFetchSegmentBySequence(sqp, s, segoffs-1, seglen, ssp, codecp)))
	ERRMSGNO(errmsgp, errcode);
      cp = seqFastqGetConstSequence(sqp, &qlen, NULL);
      printf("[%lli] %s %10u %s %-10u\n", (long long) s, snam, segoffs, cp, segoffs + qlen-1);
      if (qlen != seglen)
	printf("Warning sequence length %u not as requested!\n", qlen);

      if ((errcode = seqFastqReverse(sqp, codecp)))
	ERRMSGNO(errmsgp, errcode);
      
      cp = seqFastqGetConstSequence(sqp, &qlen, NULL);
      printf("[%lli] %s %10u %s %-10u\n", (long long) s, snam, segoffs, cp, segoffs + qlen-1);
      sso = soffs + segoffs-1; 
      printf("offset in sequence set: %llu\n", sso);
      break;
    }
  }

  seqSetDelete(ssp);
  seqFastqDelete(sqp);	 
  seqCodecDelete(codecp);
  ERRMSG_END(errmsgp);

  exit(0);
}
