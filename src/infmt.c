/** Generic interface for reading sequencing reads in different formats */

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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_BAMBAMC
#include <bambamc/BamBam_BamCollatorInterface.h>
#endif

#include "stdlib.h"
#include "infmt.h"

typedef unsigned char BOOL_t;

struct _InFmtReader {
  int errcode; /**< one of ERRMSG_CODES (benign: ERRCODE_SUCCESS, ERRCODE_EOF) */
  INFMT_t fmt;
  SeqIO *sfAp;
  SeqIO *sfBp;
#ifdef HAVE_BAMBAMC
  int fid;
  BamBam_FastQRead rbufA;
  BamBam_FastQRead rbufB;
#endif
};

#ifdef HAVE_BAMBAMC
enum SAMBAM_CONST {
  SAMBAM_LINTERM = '\0', /* don't change that */
};
static char const BAMBAM_FILTYP_BAM[] = "bam";
static char const BAMBAM_FILTYP_SAM[] = "sam";
static char const BAMBAM_TMPDIR_DEFAULT[] = ".";
#endif


/****************************************************************************
 ***************** Private Methods of type InFmtReader **********************
 ****************************************************************************/

static int openINFMTReaderAsFASTQ(InFmtReader *p, 
				  const char *filnamA, const char *filnamB, 
				  BOOL_t doTest, int *errcode_test)
{
  int errcode = ERRCODE_SUCCESS;

  /* fasta/fastq format */
  p->sfAp = seqIOopen(&errcode, filnamA, SEQIO_READ, 0);
    
  
  if (!(errcode) && filnamB != NULL) 
    p->sfBp = seqIOopen(&errcode, filnamB, SEQIO_READ, 0);
 
  if (!(errcode) && (doTest)) {
    SeqFastq *rp = seqFastqCreate(0, SEQTYP_UNKNOWN);
    *errcode_test = seqFastqRead(rp, p->sfAp);
    errcode = seqIOReset(p->sfAp);
    if (!(errcode) && !(*errcode_test) && filnamB != NULL) {
      *errcode_test = seqFastqRead(rp, p->sfBp);
      errcode = seqIOReset(p->sfBp);
    }
    seqFastqDelete(rp);
  }

  return errcode;
}

#ifdef HAVE_BAMBAMC
static int openINFMTReaderAsSAMBAM(InFmtReader *p, const char *filnam,
			    const char *tmpdir, BOOL_t is_sam)
{
  int errcode = ERRCODE_SUCCESS;
  p->fid = BamBam_AllocBamCollator(filnam, 
			    (is_sam)? BAMBAM_FILTYP_SAM: BAMBAM_FILTYP_BAM,
			    (tmpdir)? tmpdir:BAMBAM_TMPDIR_DEFAULT, 
			    1);
  if (p->fid < 0)
    errcode = ERRCODE_NOFILE;

  return errcode;
}
#endif
/****************************************************************************
 ****************** Public Methods of type InFmtReader **********************
 ****************************************************************************/

InFmtReader *infmtCreateReader(int *errcode,
			       const char *filnamA, const char *filnamB,
#ifdef HAVE_BAMBAMC
			       const char *tmpdir,
#endif
			       const INFMT_t fmt)
{
  InFmtReader *p;

  EMALLOCP0(p);

  if (p == NULL) 
    return p;
  p->errcode = ERRCODE_SUCCESS;
  p->fmt = INFMT_UNKNOWN;

#ifdef HAVE_BAMBAMC
  if (NULL == tmpdir) 
    tmpdir = BAMBAM_TMPDIR_DEFAULT;
  p->fid = -1;
#endif

  if (INFMT_UNKNOWN == fmt) {
    /* guess which format this is 
    * */
    int errcode_test = ERRCODE_SUCCESS;
    *errcode = openINFMTReaderAsFASTQ(p, filnamA, filnamB, 1, &errcode_test);
#ifdef HAVE_BAMBAMC
    if (ERRCODE_SUCCESS == *errcode && ERRCODE_SUCCESS == errcode_test) {
	p->fmt = INFMT_FASTQ;
    } else if (ERRCODE_SUCCESS == *errcode || ERRCODE_FASTA == *errcode) {
      *errcode = openINFMTReaderAsSAMBAM(p, filnamA, tmpdir, 0);
      if (ERRCODE_SUCCESS == *errcode)
	p->fmt = INFMT_BAM;
      else
	*errcode = ERRCODE_INFMT;
    }
#endif
  } else if (INFMT_FASTQ == fmt) {
    *errcode = openINFMTReaderAsFASTQ(p, filnamA, filnamB, 0, NULL);
    if (ERRCODE_SUCCESS == *errcode)
      p->fmt = INFMT_FASTQ;
    }
#ifdef HAVE_BAMBAMC
    else if (INFMT_BAM == fmt) {
      *errcode = openINFMTReaderAsSAMBAM(p, filnamA, tmpdir, 0);
    if (ERRCODE_SUCCESS == (*errcode))
      p->fmt = INFMT_BAM;
    }
    else if (INFMT_SAM == fmt) {
      *errcode = openINFMTReaderAsSAMBAM(p, filnamA, tmpdir, 1);
    if (ERRCODE_SUCCESS == (*errcode))
      p->fmt = INFMT_SAM;
    } 
#endif
  else {
    *errcode = ERRCODE_INFMT;
  }
  
  if ((*errcode)) {
    infmtDeleteReader(p);
    p = NULL;
  }

  return p;
}

void infmtDeleteReader(InFmtReader *p)
{
  if (p != NULL) {
    seqIOclose(p->sfAp);
    seqIOclose(p->sfBp);
#ifdef HAVE_BAMBAMC
    if (p->fid >= 0)
      BamBam_FreeBamCollator(p->fid);
#endif
  }
  free(p);
}

int infmtGetReaderStatus(const InFmtReader *ifrp)
{
  return ifrp->errcode;
}

int infmtRead(InFmtReader *ifrp, SeqFastq *sfqAp, SeqFastq *sfqBp, BOOL_t *isPair)
{
  int errcode;
#ifdef HAVE_BAMBAMC
  int typ;
  void *dummyA, *dummyB;
#endif
  *isPair = 0;

  seqFastqBlank(sfqAp);
  seqFastqBlank(sfqBp);

  switch(ifrp->fmt) {
  case INFMT_FASTQ:
    errcode = seqFastqRead(sfqAp, ifrp->sfAp);
    if (!(errcode) &&
	ifrp->sfBp != NULL && sfqBp != NULL) {
      errcode = seqFastqRead(sfqBp, ifrp->sfBp);
      if (ERRCODE_SUCCESS == errcode) {
	int erc = seqIOstatus(ifrp->sfBp);
	*isPair = 1;
	if (ERRCODE_SUCCESS != erc) {
	  if (ERRCODE_EOF == erc) {
	    if (seqIOstatus(ifrp->sfAp) == ERRCODE_EOF) {
	      ifrp->errcode = ERRCODE_EOF;
	    } else { 
	      errcode = ERRCODE_FASTAPAIRNUM;
	    }
	  } else {
	    errcode = erc;
	  }
	}
      }
    }
    if (ERRCODE_SUCCESS == errcode) {
      int erc = seqIOstatus(ifrp->sfAp);
      if (ERRCODE_SUCCESS != erc) {
	if (ERRCODE_SUCCESS == ifrp->errcode)
	  ifrp->errcode = erc;
	if (ERRCODE_EOF != erc)
	  errcode = erc;
      }
    }
   break;
#ifdef HAVE_BAMBAMC
  case INFMT_SAM:
  case INFMT_BAM:
    typ = BamBam_ReadPair(ifrp->fid, &ifrp->rbufA, &ifrp->rbufB, &dummyA, &dummyB,
			  SAMBAM_LINTERM);
    if (BAMBAM_ALIGNMENT_TYPE_NONE == typ) {
      errcode = ERRCODE_EOF;
    } else if (BAMBAM_ALIGNMENT_TYPE_ORPHAN2_PAIR == typ) {
      /* samlt overwrites SAM input flags */
      errcode = seqFastqSetAscii(sfqAp, 
				 ifrp->rbufB.name, ifrp->rbufB.seq, 
				 "", ifrp->rbufB.qual);
    } else {
      errcode = seqFastqSetAscii(sfqAp, 
				 ifrp->rbufA.name, ifrp->rbufA.seq, 
				 "", ifrp->rbufA.qual);
      if (ERRCODE_SUCCESS == errcode && 
	  BAMBAM_ALIGNMENT_TYPE_COMPLETE_PAIR == typ &&
	  NULL != sfqBp &&
	  ifrp->rbufB.seqlength > 0) {
	errcode = seqFastqSetAscii(sfqBp, 
				    ifrp->rbufB.name, ifrp->rbufB.seq, 
				   "", ifrp->rbufB.qual);
	if (ERRCODE_SUCCESS == errcode)
	  *isPair = 1;
      }
    }
    if ((errcode))
      ifrp->errcode = errcode;
    break;
#endif
  default:
    errcode = ERRCODE_INFMT;
  }
  
  if (!errcode) 
    errcode = seqFastqCheck(sfqAp);
  if (!errcode && sfqBp)
    errcode = seqFastqCheck(sfqBp);

  return errcode;
}

int infmtCheckReads(InFmtReader *ifrp, SeqFastq *sqbufAp, SeqFastq *sqbufBp,
		    SEQNUM_t *seqnum, SEQLEN_t *maxseqlen,  SEQLEN_t *maxnamlen, ErrMsg *errmsgp)
{
  int errcode = ERRCODE_SUCCESS;
  if ((seqnum)) *seqnum = 0;
  if ((maxseqlen)) *maxseqlen = 0;
  if ((maxnamlen)) *maxnamlen = 0;
  if (ifrp->fmt == INFMT_FASTQ) {
    errcode = seqIOCheckReads(errmsgp, sqbufAp, ifrp->sfAp, sqbufBp, ifrp->sfBp, 
			      seqnum, maxseqlen, maxnamlen);
  }

  return errcode;
}

int infmtReset(InFmtReader *ifrp)
{
  int errcode = ERRCODE_SUCCESS;

  if (ifrp->fmt == INFMT_FASTQ) {
    errcode = seqIOReset(ifrp->sfAp);
    if (!(errcode) && ifrp->sfBp != NULL)
      errcode = seqIOReset(ifrp->sfBp);
  }

  return errcode;
}
