/** Mix reads from two sets of fastq files
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

static const char FILEXT_MATE[2][6] = {"_1.fq", "_2.fq"};

int main(int argc, char *argv[])
{
  int errcode = ERRCODE_SUCCESS;
  int i, readActr = 0, readBctr = 0, ctr = 0, bnum;
  char isOK_A = 1, isOK_B = 1, isVerbose = 1;
  char *infilrootnam_A, *infilrootnam_B, *oufilrootnam;
  char *infilnam_A_mate[2], *infilnam_B_mate[2], *oufilnam_mate[2];
  SeqIO *sfp_inA[2], *sfp_inB[2], *sfp_out[2], *sfp;
  SeqFastq *seqp;
  ErrMsg *errmsg = 0;

  ERRMSG_CREATE(errmsg);

  if (argc < 4) {
    printf("usage: %s <root name fastq file A> <root name fastq file B [in] <root name fastq file [out]>\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  
  infilrootnam_A = argv[1];
  infilrootnam_B = argv[2];
  oufilrootnam = argv[3];

  if (!(seqp = seqFastqCreate(0,SEQTYP_UNKNOWN)))
    ERRMSGNO(errmsg, ERRCODE_NOMEM);

  for (i=0; i<2; i++) { 
    if (!ESTRCAT(infilnam_A_mate[i], infilrootnam_A, FILEXT_MATE[i]))
      ERRMSGNO(errmsg, ERRCODE_NOMEM);
    sfp_inA[i] = seqIOopen(&errcode, infilnam_A_mate[i], SEQIO_READ, 0);
    if (!sfp_inA[i]) 
      ERRMSGNO(errmsg, errcode);

    if (!ESTRCAT(infilnam_B_mate[i], infilrootnam_B, FILEXT_MATE[i]))
      ERRMSGNO(errmsg, ERRCODE_NOMEM);
    sfp_inB[i] = seqIOopen(&errcode, infilnam_B_mate[i], SEQIO_READ, 0);
    if (!sfp_inB[i]) 
      ERRMSGNO(errmsg, errcode);

    if (!ESTRCAT(oufilnam_mate[i], oufilrootnam, FILEXT_MATE[i]))
      ERRMSGNO(errmsg, ERRCODE_NOMEM);
    sfp_out[i] = seqIOopen(&errcode, oufilnam_mate[i], SEQIO_WRITE_FASTQ, 0);
    if (!sfp_out[i]) 
      ERRMSGNO(errmsg, errcode);
  }
  
  if ((isVerbose))
    printf("counting reads in input files %s ...\n", infilrootnam_A);
  
  isOK_A=1;
  while((isOK_A)) {
    for (i=0; i<2; i++) {
      if ((errcode = seqFastqRead(seqp, sfp_inA[i])))
	ERRMSGNO(errmsg, errcode);
      if (seqIOstatus(sfp_inA[i]))
	isOK_A = 0;
    }
    readActr++;
  }
  for (i=0; i<2; i++) {
    if ((seqIOstatus(sfp_inA[i])) && seqIOstatus(sfp_inA[i]) != ERRCODE_EOF) 
      ERRMSGNO(errmsg, seqIOstatus(sfp_inA[i]));
  }

  if ((isVerbose))
    printf("%i reads.\n", readActr);

  if ((isVerbose))
    printf("counting reads in input files %s ...\n", infilrootnam_B);
  
  isOK_B=1;
  while((isOK_B)) {
    for (i=0; i<2; i++) {
      if ((errcode = seqFastqRead(seqp, sfp_inB[i])))
	ERRMSGNO(errmsg, errcode);
      if (seqIOstatus(sfp_inB[i]))
	isOK_B = 0;
    }
    readBctr++;
  }
  for (i=0; i<2; i++) {
    if ((seqIOstatus(sfp_inB[i])) && seqIOstatus(sfp_inB[i]) != ERRCODE_EOF) 
      ERRMSGNO(errmsg, seqIOstatus(sfp_inB[i]));
  }
  if ((isVerbose))
    printf("%i reads.\n", readBctr);
 
  bnum = (readActr > readBctr)? readActr/readBctr: readBctr/readActr;
 
  for (i=0; i<2; i++) {
    seqIOReset(sfp_inA[i]);
    seqIOReset(sfp_inB[i]);
  }
 
  if ((isVerbose))
    printf("Inserting reads from file %s every %i reads in file %s ...\n",
	   (readActr > readBctr)? infilrootnam_B: infilrootnam_A,
	   bnum,
	   (readActr > readBctr)? infilrootnam_A: infilrootnam_B);
  
  isOK_A=1;
  isOK_B=1;
  while((isOK_A) || (isOK_B)) {
    for (i=0; i<2; i++) {
      sfp = (readActr > readBctr)? sfp_inA[i]: sfp_inB[i];
      if ((errcode = seqFastqRead(seqp, sfp)))
	ERRMSGNO(errmsg, errcode);
      if (seqIOstatus(sfp))
	isOK_A = 0;
      if ((errcode = seqFastqWrite(sfp_out[i], seqp, 0)))
	ERRMSGNO(errmsg, errcode);
      if (seqIOstatus(sfp_out[i]))
	isOK_A = 0;
    }

    if ((isOK_B) && !((ctr % bnum) && (isOK_A))) {
/*       if ((isVerbose)) */
/* 	printf("insert read %i ...\n", ctr); */
      for (i=0; i<2; i++) {
	sfp = (readActr > readBctr)? sfp_inB[i]: sfp_inA[i];
	if ((errcode = seqFastqRead(seqp, sfp)))
	  ERRMSGNO(errmsg, errcode);
	if (seqIOstatus(sfp))
	  isOK_B = 0;
	if ((errcode = seqFastqWrite(sfp_out[i], seqp, 0)))
	  ERRMSGNO(errmsg, errcode);
	if (seqIOstatus(sfp_out[i]))
	  isOK_B = 0;
      }
    }
    ctr++;
  }

  for (i=0; i<2; i++) {

    if ((seqIOstatus(sfp_inA[i])) && seqIOstatus(sfp_inA[i]) != ERRCODE_EOF) 
      ERRMSGNO(errmsg, seqIOstatus(sfp_inA[i]));
    if ((seqIOstatus(sfp_inB[i])) && seqIOstatus(sfp_inB[i]) != ERRCODE_EOF) 
      ERRMSGNO(errmsg, seqIOstatus(sfp_inB[i])); 
    if (seqIOstatus(sfp_out[i])) 
      ERRMSGNO(errmsg, seqIOstatus(sfp_out[i]));

    seqIOclose(sfp_out[i]);
    seqIOclose(sfp_inA[i]);
    seqIOclose(sfp_inB[i]);

    free(infilnam_A_mate[i]);
    free(infilnam_B_mate[i]);
    free(oufilnam_mate[i]);
  }

  seqFastqDelete(seqp);

  ERRMSG_END(errmsg);

  return ERRCODE_SUCCESS;
}
  
