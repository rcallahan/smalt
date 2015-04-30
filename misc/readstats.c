/** Output stats on reads like read length
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
  int ctr = 0;
  char *infilnam;
  const char *seqnam;
  SEQLEN_t seqlen;
  SeqIO *sfp_in;
  SeqFastq *seqp;
  ErrMsg *errmsg = 0;

  ERRMSG_CREATE(errmsg);

  if (argc < 2) {
    printf("usage: %s <fasta/fastq file [in]>\n", argv[0]);
    exit(0);
  }
  
  infilnam = argv[1];
 
  if (!(seqp = seqFastqCreate(0,SEQTYP_UNKNOWN)))
    ERRMSGNO(errmsg, ERRCODE_NOMEM);

  sfp_in = seqIOopen(&errcode, infilnam, SEQIO_READ, 0);
  if (!sfp_in) ERRMSGNO(errmsg, errcode);

  while(!(seqIOstatus(sfp_in))) {
    if ((errcode = seqFastqRead(seqp, sfp_in)))
      ERRMSGNO(errmsg, errcode);
    seqFastqGetSequence(seqp, &seqlen, NULL);
    seqnam = seqFastqGetSeqName(seqp);
    printf("%s %u\n", seqnam, seqlen);
    ctr++;
  }
  if ((seqIOstatus(sfp_in)) && seqIOstatus(sfp_in) != ERRCODE_EOF) 
    ERRMSGNO(errmsg, seqIOstatus(sfp_in));
  
  seqIOclose(sfp_in);
  seqFastqDelete(seqp);

  ERRMSG_END(errmsg);

  return ERRCODE_SUCCESS;
}
  
