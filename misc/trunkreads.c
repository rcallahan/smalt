/** Trunkate the reads in a FASTQ file
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
#include "elib.h"
#include "sequence.h"

int main(int argc, char *argv[])
{
  int errcode = ERRCODE_SUCCESS;
  int newreadlen, ctr = 0;
  char *infilnam; 
  SEQLEN_t rlen;
  SeqIO *sfp_in, *sfp_out;
  SeqFastq *seqp, *seqbufp;
  ErrMsg *errmsg = 0;

  ERRMSG_CREATE(errmsg);

  if (argc < 3) {
    printf("usage: %s <fasta/fastq file [in]> <new read length>\n", argv[0]);
    exit(0);
  }
  
  infilnam = argv[1];
  newreadlen = atoi(argv[2]);
  
  if (!(seqp = seqFastqCreate(0,SEQTYP_UNKNOWN)))
    ERRMSGNO(errmsg, ERRCODE_NOMEM);
  
  if (!(seqbufp = seqFastqCreate(0,SEQTYP_UNKNOWN)))
    ERRMSGNO(errmsg, ERRCODE_NOMEM);
  
  sfp_in = seqIOopen(&errcode, infilnam, SEQIO_READ, 0);
  if (!sfp_in) ERRMSGNO(errmsg, errcode);
  
  sfp_out = seqIOopen(&errcode, "-", SEQIO_WRITE_FASTQ, 0);
  if (!sfp_out) ERRMSGNO(errmsg, errcode);
    
  while(!(seqIOstatus(sfp_in) || seqIOstatus(sfp_out))) {
    if ((errcode = seqFastqRead(seqp, sfp_in)))
      ERRMSGNO(errmsg, errcode);
    ctr++;
    seqFastqGetSequence(seqp, &rlen, NULL);
    seqFastqBlank(seqbufp);
    if ((newreadlen > 0) && 
	(rlen > ((SEQLEN_t) newreadlen)))
      rlen = newreadlen;
    seqFastqAppendSegment(seqbufp, seqp, 0, rlen, 0, 0);
    seqFastqWrite(sfp_out, seqbufp, 0);
  }

  if ((seqIOstatus(sfp_in)) && seqIOstatus(sfp_in) != ERRCODE_EOF) 
    ERRMSGNO(errmsg, seqIOstatus(sfp_in));
  if (seqIOstatus(sfp_out)) 
    ERRMSGNO(errmsg, seqIOstatus(sfp_out));

  seqIOclose(sfp_out);
  seqIOclose(sfp_in);
  seqFastqDelete(seqbufp);
  seqFastqDelete(seqp);

  ERRMSG_END(errmsg);

  return ERRCODE_SUCCESS;
}
  
