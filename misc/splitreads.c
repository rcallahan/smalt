/** Save a segment of reads of a FASTQ file to another file
 */

/*****************************************************************************
 *****************************************************************************
 *                                                                           *
 *  Copyright (C) 2010-2014 Genome Research Ltd.                             * 
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
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "elib.h"
#include "sequence.h"

enum {
  LINWIDTH = 0,
  SEGMENTSIZ_REPORT = 1000000,
  NCHAR_FILEXT = 6,
};

int main(int argc, char *argv[])
{
  int errcode = ERRCODE_SUCCESS;
  int ctr = 0, filctr = 0;
  int readno_start, readno_end, readnum;
  unsigned char is_partition = 0, is_verbose = 1, isFastaOutput;
  char oufilnam[FILENAME_MAX];
  char *infilnam, *oufilnamroot; 
  SeqIO *sfp_in, *sfp_out;
  SeqFastq *seqp;
  ErrMsg *errmsg = 0;

  ERRMSG_CREATE(errmsg);

  if (argc < 6) {
    printf("usage: %s <fasta/fastq file [in]> <start read no (<0: partition)> "\
	   "<no of reads (0: all)> <convert to fasta [y/n]> <fasta/fastq file [out]>\n", argv[0]);
    exit(0);
  }
  
  infilnam = argv[1];
  readno_start = atoi(argv[2]);
  if (readno_start < 0) {
    is_partition = 1;
    readno_start = 0;
  }
  readnum = atoi(argv[3]);
  readno_end = readno_start + readnum - 1;
  isFastaOutput = (unsigned char) (toupper(argv[4][0]) == 'Y');
  oufilnamroot = argv[5];
  if (strlen(oufilnamroot) + NCHAR_FILEXT > FILENAME_MAX)
    ERRMSGNO(errmsg, ERRCODE_ARGRANGE);

  if (!(seqp = seqFastqCreate(0, (char) ((isFastaOutput)? SEQTYP_FASTA: SEQTYP_UNKNOWN))))
    ERRMSGNO(errmsg, ERRCODE_NOMEM);

  sfp_in = seqIOopen(&errcode, infilnam, SEQIO_READ, 0);
  if (!sfp_in) ERRMSGNO(errmsg, errcode);
  
  while(!(seqIOstatus(sfp_in))) {
    if (is_partition) {
      sprintf(oufilnam, "%s_%3.3i", oufilnamroot, filctr);
    } else {
      strcpy(oufilnam, oufilnamroot);
    }
    sfp_out = seqIOopen(&errcode, oufilnam, SEQIO_WRITE_FASTQ, 0);
    if (!sfp_out) ERRMSGNO(errmsg, errcode);
    if (is_verbose)
      printf("writing file %s ...\n", oufilnam);

    while(!(seqIOstatus(sfp_in) || seqIOstatus(sfp_out))) {
      if ((errcode = seqFastqRead(seqp, sfp_in)))
	ERRMSGNO(errmsg, errcode);
      if ((isFastaOutput) &&
	  (errcode = seqFastqSetType(seqp, SEQTYP_FASTA)))
	ERRMSGNO(errmsg, errcode);
      if ((errcode = seqFastqWrite(sfp_out, seqp, LINWIDTH)))
	ERRMSGNO(errmsg, errcode);
      ctr++;
      if (ctr % SEGMENTSIZ_REPORT == 0)
	printf("%i reads ... \n", ctr);
      if (ctr < readno_start)
	continue;
      if (ctr > readno_end && readnum > 0)
	break;
    }
    if (seqIOstatus(sfp_out)) 
      ERRMSGNO(errmsg, seqIOstatus(sfp_out));
    seqIOclose(sfp_out);
    if (!is_partition)
      break;
    readno_end += readnum;
    filctr++;
  }
  if ((seqIOstatus(sfp_in)) && seqIOstatus(sfp_in) != ERRCODE_EOF) 
    ERRMSGNO(errmsg, seqIOstatus(sfp_in));
  
  seqIOclose(sfp_in);
  seqFastqDelete(seqp);

  ERRMSG_END(errmsg);

  return ERRCODE_SUCCESS;
}
  
