/** Separate mates of paired-end reads into two different files
 * based on '/1' and '/2' extension of the read names.
 */

/*****************************************************************************
 *****************************************************************************
 *                                                                           *
 *  Copyright (C) 2010 -2014Genome Research Ltd.                             * 
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
#include "elib.h"
#include "sequence.h"

enum {
  LINWIDTH = 0,
  MATE_LABEL_SEPARATOR = '/',
  MATE_LABEL_1st = '1',
  MATE_LABEL_2nd = '2',
  SEGMENTSIZ_REPORT = 1000000,
};

int main(int argc, char *argv[])
{
  int errcode = ERRCODE_SUCCESS;
  int ctr = 0;
  char *infilnam, *oufilnamA, *oufilnamB, *filnamprefix; 
  const char *namstrp;
  size_t namlen;
  SeqIO *sfp_in, *sfpA_out, *sfpB_out;
  SeqFastq *seqp;
  ErrMsg *errmsg = 0;

  ERRMSG_CREATE(errmsg);

  if (argc < 3) {
    printf("usage: %s <fasta/fastq file [in]> <prefix (out)>\n", argv[0]);
    exit(0);
  }
  
  infilnam = argv[1];
  filnamprefix = argv[2];
  ESTRCAT(oufilnamA, filnamprefix, "_1.fa");
  ESTRCAT(oufilnamB, filnamprefix, "_2.fa");

  
  if (!(seqp = seqFastqCreate(0,SEQTYP_UNKNOWN)))
    ERRMSGNO(errmsg, ERRCODE_NOMEM);

  sfp_in = seqIOopen(&errcode, infilnam, SEQIO_READ, 0);
  if (!sfp_in) ERRMSGNO(errmsg, errcode);
  
  sfpA_out = seqIOopen(&errcode, oufilnamA, SEQIO_WRITE_FASTQ, 0);
  if (!sfpA_out) ERRMSGNO(errmsg, errcode);
  
  sfpB_out = seqIOopen(&errcode, oufilnamB, SEQIO_WRITE_FASTQ, 0);
  if (!sfpB_out) ERRMSGNO(errmsg, errcode);
  
  while(!(seqIOstatus(sfp_in) || seqIOstatus(sfpA_out) || seqIOstatus(sfpB_out))) {
    if ((errcode = seqFastqRead(seqp, sfp_in)))
      ERRMSGNO(errmsg, errcode);
    ctr++;
    if (ctr % SEGMENTSIZ_REPORT == 0)
      printf("%i reads ... \n", ctr);
    
    namstrp = seqFastqGetSeqName(seqp);
    namlen = strlen(namstrp);
    if (namstrp[namlen-2] == MATE_LABEL_SEPARATOR) {
      if (namstrp[namlen-1] == MATE_LABEL_1st) {
	if ((errcode = seqFastqWrite(sfpA_out, seqp, LINWIDTH)))
	  ERRMSGNO(errmsg, errcode);
      } else if (namstrp[namlen-1] == MATE_LABEL_2nd) {
	if ((errcode = seqFastqWrite(sfpB_out, seqp, LINWIDTH)))
	  ERRMSGNO(errmsg, errcode);
      } else {
	printf("Unrecogised mate label: %s\n", namstrp);
      }
    } else {
      printf("Missing mate label: %s\n", namstrp);
    }
  }
  printf("Processed %i reads.\n", ctr);

  if ((seqIOstatus(sfp_in)) && seqIOstatus(sfp_in) != ERRCODE_EOF) 
    ERRMSGNO(errmsg, seqIOstatus(sfp_in));
  if (seqIOstatus(sfpA_out)) 
    ERRMSGNO(errmsg, seqIOstatus(sfpA_out));
  if (seqIOstatus(sfpB_out)) 
    ERRMSGNO(errmsg, seqIOstatus(sfpB_out));

  seqIOclose(sfpB_out);
  seqIOclose(sfpA_out);
  seqIOclose(sfp_in);
  seqFastqDelete(seqp);
  free(oufilnamA);
  free(oufilnamB);

  ERRMSG_END(errmsg);

  return ERRCODE_SUCCESS;
}
  
