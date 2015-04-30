/** Store and output of alignment results */

/*****************************************************************************
 *****************************************************************************
 *                                                                           *
 *  Copyright (C) 2012 - 2014 Genome Research Ltd.                           *
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

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <stdint.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_BAMBAMC
#include <bambamc/BamBam_BamWriter.h>
#include <bambamc/BamBam_BamFlagBase.h>
#endif

#include "elib.h"
#include "array.h"
#include "report.h"

enum CONSTANTS {
  REPORT_DEFBLKSZ = 256,     /**< Default block size for memory allocation */
  MAXLINWIDTH_ALI = 120,     /**< Maximum line width for alignment output */
  DEFAULT_LINWIDTH_ALI = 60, /**< Default line width for output */
  DEFAULT_BLOCKSIZ_DIFFSTRBUF = 256,
  DEFAULT_BLKSZ_REPSTR = 512,/**< Block size for string memory allocation */
  NAMEXT_MAXLEN = 16,        /**< Maximum length of a read name extension */
  SEQNAM_SAM_MAXLEN = 512,   /**< Maximum length of a sequence name in the
			      * SAM header */
};

enum OUTPUT_FORMAT_LABELS {
  OUFMT_SENSE_FORWARD = '+',
  OUFMT_SENSE_REVERSE = '-',
  OUFMT_SSAHA_SENSE_FORWARD = 'F',
  OUFMT_SSAHA_SENSE_REVERSE = 'C',
  OUFMT_SENSE_UNKNOWN = '*',
  OUFMT_CLIP_HARD     = 'H',
  OUFMT_CLIP_SOFT     = 'S',
  OUFMT_NAMSTR_MATESEP = '/',
  OUFMT_NAMSTR_MATE1 = '1',
  OUFMT_NAMSTR_MATE2 = '2',
  OUFMT_CIGAR_MAXTAG = 99,    /**< Maximum mapping score if score is output after the tag,
			       * separeted by '::' from the tag */  
  OUFMT_BAM_NOQUAL = -1,      /**< Signals no quality value in BAM format (char) -1 or (unsinged char) 0xff */
};

enum SAM_FLAGS { /**< Flags for SAM output format */
  SAMFLAG_PAIRED     = 0x0001, /**< Signals paired-end read (irrespective of whether or not it is
				* mapped as a pair) */
  SAMFLAG_PROPER     = 0x0002, /**< Read is mapped as a proper pair */
  SAMFLAG_NOMAP      = 0x0004, /**< Query sequence itsself was not mapped */
  SAMFLAG_MATENOMAP  = 0x0008, /**< The mate of the query sequence was not mapped */
  SAMFLAG_STRAND     = 0x0010, /**< Strand of the query sequence (0 forward, 1 reverse) */
  SAMFLAG_MATESTRAND = 0x0020, /**< Strand of the mate (0 forward, 1 reverse) */
  SAMFLAG_1stMATE    = 0x0040, /**< Read is first mate of a pair */
  SAMFLAG_2ndMATE    = 0x0080, /**< Read is 2nd mate of a pair */
  SAMFLAG_NOTPRIMARY = 0x0100, /**< Alignment is not primary (used for split reads) */
};

enum CIGAR_FLAGS { /**< Flags for CIGAR output format */
  CIGFLG_PROPER_INSIDE          = 'A', /**< mates in proper orientations within limits */
  CIGFLG_PROPER_OUTSIDE         = 'B', /**< mates in proper orientations outside limits but
					* on same reference */
  CIGFLG_IMPROPER_SAMEREF       = 'C', /** mates not in proper orientations but on same contig */
  CIGFLG_DIFFREF                = 'D', /** mates on different contigs */
  CIGFLG_NOMAP                  = 'N', /**< unmapped read */
  CIGFLG_SINGLE                 = 'S', /**< either single reads or only mapped read of pair */
  CIGFLG_PARTIAL                = 'P', /**< secondary partial alignment of a split read */
  CIGFLG_MULTI                  = 'R', /**< Read is not placed because of multiple possible mapping locations */
};

enum ALIGNMENT_MATCHTYPES {
  ALIMATCHTYP_MATCH =    ' ',
  ALIMATCHTYP_UNKNOWN =  '?',  /**< unknown nucleotide */
  ALIMATCHTYP_NONSTD =   '!',  /**< non-standard nucleotide code */
  ALIMATCHTYP_SAMETYP =  'i',  /**< purine -> purine or  pyrimidine -> pyrimidine */
  ALIMATCHTYP_SWITCHTYP ='v',  /**< purine -> pyrimidine or pyrimidine -> purine */
  ALIMATCHTYP_GAP =      '-',  /**< alignment gap (insertion or deletion) */
};

typedef uint16_t SAMFLAG_t;
typedef unsigned short USHORT_t;
typedef unsigned char UCHAR_t;
typedef unsigned char BOOL_t;

typedef struct _REPSTR { /**< basic string type */
  char *strp;
  size_t strl;
  size_t n_alloc;
  int blksz;
} REPSTR;

typedef struct _REPNAMBUF {
  REPSTR ref_nam;
  REPSTR mref_nam;
  REPSTR q_nam;
} REPNAMBUF;

typedef struct _REPALI { /**< Alignmment results */
  BOOL_t was_output;  /**< 0 if Alignment was not yet printed */
  REPMATEFLG_t status;/**< Combination of REPORT_MATE_FLAGS */
  int swatscor;      /**< Smith-Waterman alignment score */ 
  short mapscor;       /**< 'Mapping' score */
  SEQLEN_t q_start;  /**< Start in query sequence (counting from 1) */
  SEQLEN_t q_end;    /**< End in query sequence, counting from 1 (q_start <= q_end) */
  SEQNUM_t q_idx;    /**< Index (mate number) of query sequence in template */
  SEQLEN_t s_start;  /**< Start in reference sequence (counting from 1) */
  SEQLEN_t s_end;    /**< if sidx >= 0: End in reference sequence, counting from 1 ( s_start <= s_end) */
  SEQNUM_t s_idx;     /**< Index of reference sequence in set of refrence sequences */
  int dfo;            /**< Offset of compressed alignment string in Report.dfs.dstrp 
		       * alignment in direction of reference [s_start, s_end] */
} REPALI;

typedef struct _REPPAIR {
  REPPAIRFLG_t pairflg; /**< Combination of REPORT_PAIR_FLAGS */
  int isize;            /**< Insert size for pair (relative to 1st mate) */
  int iA;              /**< Index into Report.arAr for 1st mate */
  int iB;               /**< Index into Report.arBr for 2nd mate */
} REPPAIR;

struct _Report {  /**< Holds Results for output, (Report.arAr[0], Report.arBr[0]) 
		   *   is the paired mapping */
  REPPAIR *pairr; /**< Pairs (can be multiple degenerate mappings) */
  REPALI *arAr;   /**< Alignment results 1st mate */
  REPALI *arBr;   /**< Alignment results 2nd mate */
  DiffStr dfs;    /**< Concatenated alignment strings pointed to by REPALI.dfo */
};

struct _ReportWriter {
  REPOUFMT_t oufmt;
  REPMODIFLG_t modflg;
  short linwidth;  /**< used for explicit alignment output */
  SeqFastq *qbufp;
  SeqFastq *sbufp;
  DiffBlocks *dfblkp;  /**< used for some output formats, e.g. GFF2 */
  char namext[NAMEXT_MAXLEN];
  char namext_mate[NAMEXT_MAXLEN];
  REPNAMBUF *nambufp;
  char *filnam;
  FILE *oufp;
#ifdef HAVE_BAMBAMC
  DiffView *dvp; /**< Formatted alignment string (CIGAR) */
  BamBam_BamHeaderInfo *headinfop ;/**< BAM/SAM header */
  BamBam_BamWriter * bamwriterp;
#endif
};

#ifdef HAVE_BAMBAMC
enum SAMBAM_CONST {
  SAMBAM_COMPRESSION_LEVEL = 1,
};

static const char SAMBAM_HEADER_VERSION[] = "1.4";
static const char SAMBAM_SORTORDER_UNKNOWN[] = "unknown";
#endif

static const char OUFMT_FILNAM_STDOUT[] =  "-";

/* format strings for SAM header */
static const char SAMFORM_HEADLINE[] = "@HD\tVN:1.3\tSO:unknown\n";
static const char SAMFORM_REFSEQLINE[] = "@SQ\tSN:%s\tLN:%u\n";
static const char SAMFORM_PROGLINE[] = "@PG\tID:%s\tPN:%s\tVN:%s\tCL:";

/* format strings for SAM line */
static const char OUFMT_SAM_BEFORE[] = "%s\t%hu\t%s\t%i\t%hi\t";
/**< QNAME FLAG RNAME POS MAPQ */
static const char OUFMT_SAM_AFTER[] = "\t%s\t%i\t%i\t%s\t%s\tNM:i:%i\tAS:i:%i\n";
/**< MRNM MPOS ISIZE SEQ QUAL */

/**< TAG VTYPE VALUE */
static const char OUFMT_SAM_NULLSTR[] = "*";
static const char OUFMT_BAM_NULLSTR[] = "";

static const char OUFMT_CIGAR[] = "cigar:%c:%2.2d %s%s %u %u %c %s %u %u + %d ";
static const char OUFMT_ALIGN[] =\
"    QUERY: %10i %s %-10i\n                      %s\n"\
"REFERENCE: %10i %s %-10i\n\n\n";

static const char OUFMT_SSAHA[] = "alignment:%c:%2.2d %-5d %s%s %s %8u %8u %9u %9u   %c %7d %5.2f %u %u\n";
static const char OUFMT_GFF2[] = "gff: %s%s\t" \
  "SMALT\tsimilarity\t%d\t%d\t%d\t%c\t.\t" \
  "Subject \"%s\" %d %d;\t";
static const char OUFMT_GFF2_ALIBLOCK[] = " Align %d %d %d;";
static const char OUFMT_GFF2_NULLSTR[] = "";


/******************************************************************************
 ******************************* Private Methods ******************************
 ******************************************************************************/
static char getMapLabelFromFlag(REPMATEFLG_t mateflg, REPPAIRFLG_t pairflg)
{
  char flagchr;

  if (mateflg & REPMATEFLG_MAPPED) {
    if (mateflg & REPMATEFLG_PARTIAL) {
      flagchr = CIGFLG_PARTIAL;
    } else if (pairflg & REPPAIR_MAPPED) {
      if (pairflg & REPPAIR_CONTIG) {
	if (pairflg & REPPAIR_PROPER) {
	  flagchr = (char) ((pairflg & REPPAIR_WITHIN)? 
			    CIGFLG_PROPER_INSIDE: 
			    CIGFLG_PROPER_OUTSIDE);
	} else {
	  flagchr = CIGFLG_IMPROPER_SAMEREF;
	}
      } else {
	flagchr = CIGFLG_DIFFREF;
      }
    } else {
      flagchr = CIGFLG_SINGLE;
    }
  } else if (mateflg & REPMATEFLG_MULTI) {
    flagchr = CIGFLG_MULTI;
  } else {
    flagchr = CIGFLG_NOMAP;
  }

  return flagchr;
}

static int fprintAlignment(FILE *fp, SeqFastq *sbufp, SeqFastq *qbufp, 
			   const SeqFastq *sqp,
			   USHORT_t linwidth, const DiffStr *dfsp, 
			   const SeqSet *ssp, const SeqCodec *codecp, 
			   const REPALI *rp)
     /**< Print out the alignment explicitly.
      */
{
  int errcode;
  char cod, qbuf[MAXLINWIDTH_ALI], dbuf[MAXLINWIDTH_ALI], sbuf[MAXLINWIDTH_ALI];
  char q_bastyp, s_bastyp;
  const char *qcp, *scp;
  UCHAR_t k, typ, count;
  USHORT_t j;
  SEQLEN_t s, q, len, qlen, slen, q_linstart, s_linstart;
  const DIFFSTR_T *ucp, *diffstrp = dfsp->dstrp + rp->dfo;

  if (!rp)
    return ERRCODE_ASSERT;

  if (seqCodecType(codecp) != SEQCODTYP_3BITMANGLED)
    return ERRCODE_ASSERT;

  seqFastqGetConstSequence(sqp, &qlen, &cod);

  if (rp->q_start < 1 || rp->q_start > rp->q_end || rp->q_end > qlen)
    return ERRCODE_ASSERT;

  if (linwidth >= MAXLINWIDTH_ALI) linwidth = MAXLINWIDTH_ALI-1;
  else if (linwidth < 1)  linwidth = DEFAULT_LINWIDTH_ALI;

  seqFastqBlank(qbufp);
  if (rp->status & REPMATEFLG_REVERSE) {
    if ((errcode = seqFastqAppendSegment(qbufp, sqp, rp->q_start-1, 
					 rp->q_end - rp->q_start + 1,
					 1, codecp)))
      return errcode;
  } else {
    if ((errcode = seqFastqAppendSegment(qbufp, sqp, rp->q_start-1, 
					 rp->q_end - rp->q_start + 1,
					 0, NULL)))
      return errcode;
  }

  if (cod == SEQCOD_MANGLED &&
      (errcode = seqFastqDecode(qbufp, codecp)))
    return errcode;

  qcp = seqFastqGetConstSequence(qbufp, &qlen, &cod);
  
  if (cod != SEQCOD_ASCII)
    return ERRCODE_SEQCODE;

  if (rp->s_start < 1 || rp->s_start > rp->s_end)
    return ERRCODE_ASSERT;
  len = rp->s_end - rp->s_start + 1;
  if (rp->s_idx < 0) 
    return ERRCODE_NOSEQIDX;
  if ((errcode = seqSetFetchSegmentBySequence(sbufp, rp->s_idx, 
					      (SEQLEN_t) rp->s_start-1, 
					      len, ssp, codecp)))
    return errcode;
  scp = seqFastqGetConstSequence(sbufp, &slen, NULL);
  if (slen != len)
    return ERRCODE_ASSERT;

  ucp = diffstrp;
  q_linstart = s_linstart = 0;
  DIFFSTR_GET(*ucp, count, typ);

  for (k=0, q=s=0; (*ucp) && q<=qlen && s<=slen;) {
    for (j=0; j<linwidth && q<=qlen && s<=slen; j++) {  

      if (k++ < count) {
	qbuf[j] = qcp[q++];
	dbuf[j] = ALIMATCHTYP_MATCH;
	sbuf[j] = scp[s++];
	continue;
      } 
      k = 0;
      ucp++;

      if (typ == DIFFCOD_M) {
	qbuf[j] = qcp[q++];
	dbuf[j] = ALIMATCHTYP_MATCH;
	sbuf[j] = scp[s++];
      } else if (typ == DIFFCOD_S) {
	qbuf[j] = qcp[q];
	sbuf[j] = scp[s];
	if (*(ucp)) {
	  /* substitution, figure out which type of substitution */
	  q_bastyp = seqCodecFindBaseClass(qcp[q], codecp);
	  s_bastyp = seqCodecFindBaseClass(scp[s], codecp);
	  if (q_bastyp == SEQBASTYP_NONSTD || s_bastyp == SEQBASTYP_NONSTD)
	    dbuf[j] = ALIMATCHTYP_NONSTD;
	  else if (q_bastyp == SEQBASTYP_UNKNOWN || s_bastyp == SEQBASTYP_UNKNOWN)
	    dbuf[j] = ALIMATCHTYP_UNKNOWN;
	  else if (q_bastyp == s_bastyp)
	    dbuf[j] = ALIMATCHTYP_SAMETYP;
	  else 
	    dbuf[j] = ALIMATCHTYP_SWITCHTYP;
	  q++;
	  s++;
	} else { /* termination signal */
	  dbuf[j] = ALIMATCHTYP_MATCH;
	}
      } else if (typ == DIFFCOD_D) {
	qbuf[j] = ALIMATCHTYP_GAP;
	dbuf[j] = ALIMATCHTYP_GAP;
	sbuf[j] = scp[s++];
      } else if (typ == DIFFCOD_I) {
	qbuf[j] = qcp[q++];
	dbuf[j] = ALIMATCHTYP_GAP;
	sbuf[j] = ALIMATCHTYP_GAP;
      }
   
      if (!(*ucp)) break;
      DIFFSTR_GET(*ucp, count, typ);
    } /* for (j=0; j<linwidth; j++) */

    /* print the line */
    qbuf[j] = dbuf[j] = sbuf[j] = '\0';
    fprintf(fp, OUFMT_ALIGN,				  
	    (rp->status&REPMATEFLG_REVERSE)? 
	    rp->q_end-q_linstart: 
	    rp->q_start+q_linstart, 
	    qbuf,
	    (rp->status&REPMATEFLG_REVERSE)? 
	    rp->q_end - q + 1: 
	    rp->q_start + q - 1, 
	    dbuf,
	    (int) rp->s_start+s_linstart, sbuf,
	    (int) rp->s_start+s-1);

    s_linstart = s;
    q_linstart = q; 
    j = 0;
  } /* (*ucp) */
  
  return ((*ucp))? ERRCODE_DIFFSTR: ERRCODE_SUCCESS;
}

/******************************************************************************
 *********************** Private Methods of Type REPSTR ***********************
 ******************************************************************************/

static void cleanupREPSTR(REPSTR *rsp)
{
  if (rsp) {
    free(rsp->strp);
    rsp->strp = NULL;
    rsp->strl = rsp->n_alloc = 0;
  }
}

static int initREPSTR(REPSTR *rsp, int blksz)
{
  if (NULL == rsp) return ERRCODE_ASSERT;
  if (blksz < 1) blksz = DEFAULT_BLKSZ_REPSTR;
  rsp->blksz = blksz;
  rsp->strl = 0;
  rsp->n_alloc = 0;
  ECALLOCP(blksz, rsp->strp);
  if (NULL == rsp->strp) return ERRCODE_NOMEM;
  rsp->n_alloc = blksz;
  rsp->strp[0] = '\0';
  return ERRCODE_SUCCESS;
}

static int reallocREPSTR(REPSTR *rsp, size_t newlen)
{
  int errcode = ERRCODE_NOMEM;
  size_t newsz = (newlen + rsp->blksz - 1)/rsp->blksz;
  void *hp;
  newsz *= rsp->blksz;

  hp = EREALLOCP(rsp->strp, newsz);
  if (hp != NULL) {
    rsp->strp = hp;
    rsp->n_alloc = newsz;
    errcode = ERRCODE_SUCCESS;
  }

  return errcode;
}

static int copyReadNamStrToREPSTR(REPSTR *rsp, 
				  BOOL_t is_stripped, const char *namp)
{
  int errcode = ERRCODE_SUCCESS;
  int i, c;

  for (i=0; i<INT_MAX-1; i++) {
    c = (int) namp[i];
    if (!c)
      break;
    if (isspace(c))
      break;
    if (((size_t) i) + 1 >= rsp->n_alloc && 
	(errcode = reallocREPSTR(rsp, i+1)))
      break;
    rsp->strp[i] = (char) c;
  }

  /* don't copy the '/1', '/2' extension */
  if ((is_stripped) && i>2 && rsp->strp[i-2] == OUFMT_NAMSTR_MATESEP &&
      (rsp->strp[i-1] == OUFMT_NAMSTR_MATE1 || rsp->strp[i-1] == OUFMT_NAMSTR_MATE2)) 
    i -= 2;

  rsp->strp[i] = '\0';
  rsp->strl = (size_t) i;
  
  return (i < INT_MAX)? errcode: ERRCODE_SEQNAMLEN;
}

static int copyReadNameToREPSTR(REPSTR *rsp, BOOL_t is_stripped, const SeqFastq *sqp)
{
  const char *namp;

  namp = (NULL == sqp)? OUFMT_SAM_NULLSTR: seqFastqGetSeqName(sqp);
 
  return copyReadNamStrToREPSTR(rsp, is_stripped, namp);
}

#ifdef HAVE_BAMBAMC
static int copySAMheaderCommandLineToREPSTR(REPSTR *rsp, 
					    const char *prognam,
					    const char *progversion,
					    int narg, 
					    char * const *argv)
{
  int errcode = ERRCODE_SUCCESS;
  int i, nc;
  size_t maxlen = strlen(SAMFORM_PROGLINE)+1;

  rsp->strp[0] = '\0';
  rsp->strl = 0;
  if (narg < 1 || NULL == argv)
    return ERRCODE_SUCCESS;
  /* calculate max length of line */
  for (i=0; i<narg; i++)
    maxlen += strlen(argv[i]) + 1;

  if (maxlen > rsp->n_alloc &&
      (errcode = reallocREPSTR(rsp, maxlen)))
    return errcode;

  if ((nc = sprintf(rsp->strp, SAMFORM_PROGLINE, prognam, prognam, progversion)) < 1)
    return ERRCODE_WRITEERR; 
  rsp->strl += (size_t) nc;

  if ((nc = sprintf(rsp->strp + rsp->strl, "%s", argv[0])) < 1)
    return ERRCODE_WRITEERR;
  rsp->strl += (size_t) nc;

  for (i=1; i<narg; i++) {
    if ((nc = sprintf(rsp->strp + rsp->strl, " %s", argv[i])) < 1)
	return ERRCODE_WRITEERR;
    rsp->strl += (size_t) nc;
  }
  if ((nc = sprintf(rsp->strp + rsp->strl, "\n")) < 1) 
    return ERRCODE_WRITEERR;
  rsp->strl += (size_t) nc;
	
  return ERRCODE_SUCCESS;
}
#endif
/******************************************************************************
 ********************* Private Methods of Type REPNAMBUF **********************
 ******************************************************************************/
static void deleteREPNAMBUF(REPNAMBUF *p)
{
  if (p) {
    cleanupREPSTR(&p->q_nam);
    cleanupREPSTR(&p->mref_nam);
    cleanupREPSTR(&p->ref_nam);
  }
  free(p);
}

static REPNAMBUF *createREPNAMBUF(void)
{
  REPNAMBUF *p;
  EMALLOCP0(p);
  if ((NULL == p) ||
      (initREPSTR(&p->ref_nam, 0)) ||
      (initREPSTR(&p->mref_nam, 0)) ||
      (initREPSTR(&p->q_nam, 0))) {
    deleteREPNAMBUF(p);
    p = NULL;
  }
  return p;
}

/******************************************************************************
 *********************** Private Methods of Type REPALI ***********************
 ******************************************************************************/
static int findREPALI(const REPALI *rar, 
		      int *idxp,
		      SEQLEN_t q_start, SEQLEN_t q_end,
		      REPMATEFLG_t rmatflg,
		      SEQLEN_t s_start, SEQLEN_t s_end, SEQNUM_t s_idx)
/**< check whether there is an alignment exactly like this already in the 
 * array rar. If there is return ERRCODE_SUCCESS and the indes in idxp,
 * otherwise return ERRCODE_FAILURE and *idxp < 0;
 */
{
  const size_t n = ARRLEN(rar);
  const REPMATEFLG_t mask = REPMATEFLG_REVERSE | REPMATEFLG_2NDMATE;
  int i;
  *idxp = -1;
  if (n < 1)
    return ERRCODE_FAILURE;

  if (n > INT_MAX)
    return ERRCODE_OVERFLOW;

  for (i=(int) n-1; 
       i >=0 &&
	 (s_start != rar[i].s_start ||
	  s_end !=rar[i].s_end ||
	  s_idx != rar[i].s_idx ||
	  q_start != rar[i].q_start ||
	  q_end != rar[i].q_end ||
	  (rmatflg & mask) !=  (rar[i].status & mask));
       i--);
  *idxp = i;

  return (i < 0)? ERRCODE_FAILURE: ERRCODE_SUCCESS;
}
		  
static int fprintREPALIssaha(FILE *fp, const REPALI *rp, short mapscor,
			     REPSTR *q_namp,
			     const SeqFastq *q_sqp, const char *q_namext, 
			     const char *s_nam, SEQLEN_t s_len, const DIFFSTR_T *diffstr,
			     REPPAIRFLG_t pairflg)
{
  int errcode, swatscor;
  SEQLEN_t qs, qe, qlen;
  SETSIZ_t rs, re;
  int matchlen=0, alilen = 0;
  char sensechr, flagchr;
  double idfrac = .0;

  if ((errcode = copyReadNameToREPSTR(q_namp, 0, q_sqp)))
    return errcode;

  seqFastqGetConstSequence(q_sqp, &qlen, NULL);
  
  if ((rp != NULL) && (rp->status & REPMATEFLG_MAPPED)) {
    if (rp->status&REPMATEFLG_REVERSE) {
      qs = rp->q_end;
      qe = rp->q_start;
      sensechr = OUFMT_SSAHA_SENSE_REVERSE;
    } else {
      qs = rp->q_start;
      qe = rp->q_end;
      sensechr = OUFMT_SSAHA_SENSE_FORWARD;
    }
    rs = rp->s_start;
    re = rp->s_end;
    swatscor = rp->swatscor;

    flagchr = getMapLabelFromFlag(rp->status, pairflg);

   if ((diffstr)) {
      alilen = diffStrCalcAliLen(&matchlen, diffstr);
      idfrac = (alilen > 0)? ((double) 100*matchlen)/alilen: .0;
    } 
  } else {
    qs = qe = rs = re = 0;
    sensechr = OUFMT_SENSE_UNKNOWN;
    s_nam = OUFMT_SAM_NULLSTR;
    s_len = 0;
    swatscor = 0;
    mapscor = 0;
    flagchr = (char) (((rp) && (rp->status & REPMATEFLG_MULTI))? CIGFLG_MULTI: CIGFLG_NOMAP);
    idfrac = .0;
  }

#ifdef report_debug
  printf("report_debug: alilen = %i, matchlen = %i: ", alilen, matchlen);
  errcode = diffStrPrintf(stdout, diffstr, DIFFSTRFORM_RAW, 0, 0, 0);
#endif
  fprintf(fp, OUFMT_SSAHA, flagchr,
	  (mapscor > OUFMT_CIGAR_MAXTAG)? OUFMT_CIGAR_MAXTAG:mapscor,
	  swatscor,
	  q_namp->strp, q_namext, 
	  s_nam,
	  (unsigned int) qs, (unsigned int) qe,
	  (unsigned int) rs, (unsigned int) re,
	  sensechr,
	  matchlen,
	  idfrac,
	  qlen,
	  s_len);
  
  return errcode;
}

static int fprintREPALIgff2(FILE *fp, const REPALI *rp,
			    REPSTR *q_namp,
			    const SeqFastq *q_sqp, const char *q_namext, 
			    const char *s_nam, 
			    const DiffBlocks *dfblkp)
{
  int errcode, swatscor;
  int b, blkn = diffBlocksGetNumber(dfblkp);
  SEQLEN_t qs, qe, qlen;
  SETSIZ_t rs, re;
  char sensechr;
  BOOL_t isReverse = (BOOL_t) ((rp) && (rp->status&REPMATEFLG_REVERSE));

  if ((errcode = copyReadNameToREPSTR(q_namp, 0, q_sqp)))
    return errcode;

  seqFastqGetConstSequence(q_sqp, &qlen, NULL);
  
  if ((rp != NULL) && (rp->status & REPMATEFLG_MAPPED)) {
    if (isReverse) {
      qs = rp->q_end;
      qe = rp->q_start;
      sensechr = OUFMT_SENSE_REVERSE;
    } else {
      qs = rp->q_start;
      qe = rp->q_end;
      sensechr = OUFMT_SENSE_FORWARD;
    }
    rs = rp->s_start;
    re = rp->s_end;
    swatscor = rp->swatscor;

  } else {
    qs = qe = rs = re = 0;
    sensechr = OUFMT_SENSE_UNKNOWN;
    s_nam = OUFMT_GFF2_NULLSTR;
    swatscor = 0;
  }

  fprintf(fp, OUFMT_GFF2,
	  q_namp->strp, q_namext, 
	  (unsigned int) qs, (unsigned int) qe,
	  swatscor, sensechr,
	  s_nam,
	  (unsigned int) rs, (unsigned int) re);
  
  for (b=0; b<blkn; b++) {
    int q0 = 0;
    int r0 = 0;
    int len = diffBlocksGetLen(&r0, &q0, b, dfblkp);
    if (len < 1) break;
    if (isReverse) {
      q0 = rp->q_end - rp->q_start - q0;
    } 
    fprintf(fp, OUFMT_GFF2_ALIBLOCK, q0 + 1, r0 + 1, len);
  }
  if (b == 0)
    fprintf(fp, OUFMT_GFF2_ALIBLOCK, 0, 0, 0);
  fprintf(fp, "\n");

  return errcode;
}

static int fprintREPALIcigar(FILE *fp, const REPALI *rp, short mapscor,
			     REPSTR *q_namp,
			     const SeqFastq *q_sqp, const char *q_namext, 
			     const char *s_nam, const DIFFSTR_T *diffstr, 
			     REPPAIRFLG_t pairflg)
{
  int errcode, swatscor;
  SEQLEN_t qs, qe;
  SETSIZ_t rs, re;
  char sensechr, flagchr;

  if ((errcode = copyReadNameToREPSTR(q_namp, 0, q_sqp)))
    return errcode;

  if ((rp != NULL) && (rp->status & REPMATEFLG_MAPPED)) {
    if ((rp->status & REPMATEFLG_REVERSE)) {
      qs = rp->q_end;
      qe = rp->q_start;
      sensechr = OUFMT_SENSE_REVERSE;
    } else {
      qs = rp->q_start;
      qe = rp->q_end;
      sensechr = OUFMT_SENSE_FORWARD;
    }
    rs = rp->s_start;
    re = rp->s_end;
    swatscor = rp->swatscor;

    flagchr = getMapLabelFromFlag(rp->status, pairflg);
  } else {
    qs = qe = rs = re = 0;
    sensechr = OUFMT_SENSE_UNKNOWN;
    s_nam = OUFMT_SAM_NULLSTR;
    swatscor = 0;
    mapscor = 0;
    flagchr = (char) (((rp->status & REPMATEFLG_MULTI))? CIGFLG_MULTI: CIGFLG_NOMAP);
  }
      
  fprintf(fp, OUFMT_CIGAR, flagchr,
	  //mapscor*OUFMT_CIGAR_MAXTAG/MAPSCOR_MAX,
	  (mapscor > OUFMT_CIGAR_MAXTAG)? OUFMT_CIGAR_MAXTAG:mapscor,
	  q_namp->strp, q_namext, qs, qe,
	  sensechr,
	  s_nam, 
	  (unsigned int) rs, (unsigned int) re,
	  swatscor);
  errcode = diffStrPrintf(fp, diffstr, DIFFSTRFORM_CIGNORM, 0, 0, 0);
  fprintf(fp, "\n");
  return errcode;
}

static int fprintREPALIsam(FILE *fp, SeqFastq *sqbufp, 
			   REPSTR *readnamp,
			   short mapscor,
			   const REPALI *rrp, const DIFFSTR_T *diffstr,
			   const SeqFastq *q_sqp,
			   const char *refnam,
			   const REPALI *rmp,
			   const char *mate_refnam,
			   int isize,
			   REPPAIRFLG_t pairflg,
			   REPMODIFLG_t oumodiflg,
			   const SeqCodec *codecp)
{
  int errcode;
  int editdist=0, clip_start=0, clip_end=0;
  char cod;
  SEQLEN_t qlen;
  SAMFLAG_t samflg = 0;
  const char *seqstr, *qualstr, *ms_nam = mate_refnam, *s_nam = refnam;
  int swatscor;
  SEQLEN_t pos=0, mpos=0;

  if (!rrp) 
    return ERRCODE_ASSERT;

  if (!(sqbufp))
    return ERRCODE_NULLPTR;
  seqFastqBlank(sqbufp);

  if ((errcode = copyReadNameToREPSTR(readnamp, 1, q_sqp)))
    return errcode;
    
  seqFastqGetConstSequence(q_sqp, &qlen, &cod);

  if (rrp->status & REPMATEFLG_PAIRED) {
    samflg |= SAMFLAG_PAIRED;
    if ((rrp->status & REPMATEFLG_2NDMATE)) {
      samflg |= SAMFLAG_2ndMATE;
      isize *= -1;
    } else {
      samflg |=  SAMFLAG_1stMATE;
    }
 
    if ((rmp) && (rmp->status & REPMATEFLG_MAPPED)) {
      mpos = (SEQLEN_t) rmp->s_start;
      if (rmp->status&REPMATEFLG_REVERSE)
	samflg |= SAMFLAG_MATESTRAND;
    } else {
      samflg |= SAMFLAG_MATENOMAP;
      isize = 0;
      mpos = 0;
      ms_nam = OUFMT_SAM_NULLSTR;
    }		
  } else {
    ms_nam = OUFMT_SAM_NULLSTR;
  }
  
  if ((rrp->status & REPMATEFLG_MAPPED)) {
    BOOL_t isReverse = (BOOL_t) ((rrp->status & REPMATEFLG_REVERSE)? 1 : 0);
    SEQLEN_t qseg_start, qseg_len;

    if (oumodiflg & REPORTMODIF_SOFTCLIP) {
      qseg_start = 0;
      qseg_len = qlen;
    } else {
      qseg_start = rrp->q_start-1;
      qseg_len = rrp->q_end - rrp->q_start + 1;
    }
    errcode = seqFastqAppendSegment(sqbufp, q_sqp, 
				    qseg_start, qseg_len,
				    (isReverse), (isReverse)? codecp: NULL);
    if (errcode)
      return errcode;
    if (cod == SEQCOD_MANGLED &&
	(errcode = seqFastqDecode(sqbufp, codecp)))
      return errcode;
    
    seqstr = seqFastqGetConstSequence(sqbufp, NULL, &cod);
    if (cod != SEQCOD_ASCII)
      return ERRCODE_SEQCODE;

    qualstr = seqFastqGetConstQualityFactors(sqbufp, NULL, NULL);

    pos = (SEQLEN_t) rrp->s_start;
    if (rrp->q_end > qlen)
      return ERRCODE_ASSERT;

    if (isReverse) {
      samflg |= SAMFLAG_STRAND;
      clip_start = qlen - rrp->q_end;
      clip_end = rrp->q_start - 1;
    } else {
      clip_start = rrp->q_start - 1;
      clip_end = qlen - rrp->q_end;
    }

    if (((pairflg & REPPAIR_PROPER) && (pairflg & REPPAIR_WITHIN)))
      samflg |= SAMFLAG_PROPER;
    if (rrp->status & REPMATEFLG_PARTIAL)
      samflg |= SAMFLAG_NOTPRIMARY;

    swatscor = rrp->swatscor;
  } else { /* if (rrp->status & REPMATEFLG_MAPPED) */
    if (oumodiflg & REPORTMODIF_SOFTCLIP) {
      if ((errcode = seqFastqAppendSegment(sqbufp, q_sqp, 0, 0, 0, NULL)))
	return errcode;
      if (cod == SEQCOD_MANGLED &&
	  (errcode = seqFastqDecode(sqbufp, codecp)))
	return errcode;
      seqstr = seqFastqGetConstSequence(sqbufp, NULL, &cod);
      if (cod != SEQCOD_ASCII)
	return ERRCODE_SEQCODE;
      qualstr = seqFastqGetConstQualityFactors(sqbufp, NULL, NULL);
    } else { 
      seqstr = OUFMT_SAM_NULLSTR;
      qualstr = OUFMT_SAM_NULLSTR;
    }
    samflg |= SAMFLAG_NOMAP;
    s_nam = OUFMT_SAM_NULLSTR;
    swatscor = 0;
    isize = 0;
  }

  if (!qualstr || !qualstr[0])
    qualstr = OUFMT_SAM_NULLSTR;

  fprintf(fp, OUFMT_SAM_BEFORE, 
	  readnamp->strp, samflg, s_nam, pos, mapscor);

  if (rrp->status & REPMATEFLG_MAPPED) {
    errcode = diffStrPrintf(fp, diffstr, 
			    (char) ((oumodiflg & REPORTMODIF_XMISMATCH)? 
				    DIFFSTRFORM_CIGEXT_XMISMATCH: DIFFSTRFORM_CIGEXT), 
				    clip_start, clip_end, 
			    (char) ((oumodiflg & REPORTMODIF_SOFTCLIP) != 0));
    if (!errcode)
      editdist = diffStrGetLevenshteinDistance(diffstr);
  } else {
    fprintf(fp, OUFMT_SAM_NULLSTR);
  }
  fprintf(fp, OUFMT_SAM_AFTER, 
	  ms_nam, mpos, isize, seqstr, qualstr, editdist, swatscor);

  return errcode;
}
#ifdef report_debug
void fprintREPALIraw(FILE *fp, const REPALI *p)
{
  fprintf(fp, "REPALI: q_start=%u, q_end=%u, q_idx=%lli, "\
	  "s_start=%u, s_end=%u, s_idx=%lli\n",
	  p->q_start, p->q_end, (long long signed) p->q_idx,
	  p->s_start, p->s_end, (long long signed) p->s_idx);
}
#endif
#ifdef HAVE_BAMBAMC
static int writeREPALIbam(BamBam_BamWriter *bamwriterp, 
			  SeqFastq *sqbufp, 
			  REPSTR *readnamp,
			  DiffView *dvp,
			  short mapscor,
			  const REPALI *rrp, const DIFFSTR_T *diffstr,
			  const SeqFastq *q_sqp,
			  const REPALI *rmp,
			  int isize,
			  REPPAIRFLG_t pairflg,
			  REPMODIFLG_t oumodiflg,
			  const SeqCodec *codecp)

{
  int errcode, errc;
  int editdist=0, clip_start=0, clip_end=0;
  int32_t s_idx = -1, ms_idx = -1;
  char cod;
  int bamflg = 0;
  int32_t swatscor = 0;
  SEQLEN_t pos = 0, mpos = 0;
  SEQLEN_t qlen;
  const char *seqstr, *qualstr, *cigarstr;

  if (!rrp) 
    return ERRCODE_ASSERT;

  if (!(sqbufp))
    return ERRCODE_NULLPTR;
  seqFastqBlank(sqbufp);

  if ((errcode = copyReadNameToREPSTR(readnamp, 1, q_sqp)))
    return errcode;
    
  if (rrp->s_idx > INT_MAX)
    return ERRCODE_OVERFLOW;
  s_idx = (int) rrp->s_idx;

  seqFastqGetConstSequence(q_sqp, &qlen, &cod);
  
  if (rrp->status & REPMATEFLG_PAIRED) {
    bamflg |= BAMBAMC_FPAIRED;
    if ((rrp->status & REPMATEFLG_2NDMATE)) {
      bamflg |= BAMBAMC_FREAD2;
      isize *= -1;
    } else {
      bamflg |= BAMBAMC_FREAD1;
    }
 
    if ((rmp) && (rmp->status & REPMATEFLG_MAPPED)) {
      mpos = (SEQLEN_t) rmp->s_start;
      if (rmp->s_idx > INT_MAX)
	return ERRCODE_OVERFLOW;
      ms_idx = (int) rmp->s_idx;
      if (rmp->status&REPMATEFLG_REVERSE)
	bamflg |= BAMBAMC_FMREVERSE;
    } else {
      bamflg |= BAMBAMC_FMUNMAP;
      isize = 0;
      mpos = 0;
      ms_idx = -1;
    }		
  } else {
    ms_idx = -1;
  }
  
  if ((rrp->status & REPMATEFLG_MAPPED)) {
    BOOL_t isReverse = (BOOL_t) ((rrp->status & REPMATEFLG_REVERSE)? 1 : 0);
    errcode = seqFastqAppendSegment(sqbufp, q_sqp, 
				    0, 0,
				    (isReverse), (isReverse)? codecp: NULL);
    if (errcode)
      return errcode;
    if (cod == SEQCOD_MANGLED &&
	(errcode = seqFastqDecode(sqbufp, codecp)))
      return errcode;
    
    seqstr = seqFastqGetConstSequence(sqbufp, NULL, &cod);
    if (cod != SEQCOD_ASCII)
      return ERRCODE_SEQCODE;

    qualstr = seqFastqGetConstQualityFactors(sqbufp, NULL, NULL);

    pos = (SEQLEN_t) rrp->s_start;
    if (rrp->q_end > qlen)
      return ERRCODE_ASSERT;

    if (isReverse) {
      bamflg |=  BAMBAMC_FREVERSE;
      clip_start = qlen - rrp->q_end;
      clip_end = rrp->q_start - 1;
    } else {
      clip_start = rrp->q_start - 1;
      clip_end = qlen - rrp->q_end;
    }

    if (((pairflg & REPPAIR_PROPER) && (pairflg & REPPAIR_WITHIN)))
      bamflg |= BAMBAMC_FPROPER_PAIR;
    if (rrp->status & REPMATEFLG_PARTIAL)
      bamflg |= BAMBAMC_FSECONDARY;

    swatscor = rrp->swatscor;

  } else { /* if (rrp->status & REPMATEFLG_MAPPED) */
    if ((errcode = seqFastqAppendSegment(sqbufp, q_sqp, 0, 0, 0, NULL)))
      return errcode;
    if (cod == SEQCOD_MANGLED &&
	(errcode = seqFastqDecode(sqbufp, codecp)))
      return errcode;
    seqstr = seqFastqGetConstSequence(sqbufp, NULL, &cod);
    if (cod != SEQCOD_ASCII)
      return ERRCODE_SEQCODE;
    qualstr = seqFastqGetConstQualityFactors(sqbufp, NULL, NULL);
    bamflg |= BAMBAMC_FUNMAP;
    s_idx = -1;
    swatscor = 0;
    isize = 0;
  }

  if (!qualstr || !qualstr[0]) {
    if ((errcode = seqFastqSetQual(sqbufp, OUFMT_BAM_NOQUAL)))
      return errcode;
    qualstr = seqFastqGetConstQualityFactors(sqbufp, NULL, NULL);
  }

  if (rrp->status & REPMATEFLG_MAPPED) {
    errcode = diffStrAsView(dvp, diffstr, 
			    (char) ((oumodiflg & REPORTMODIF_XMISMATCH)? 
				    DIFFSTRFORM_CIGEXT_XMISMATCH:DIFFSTRFORM_CIGEXT), 
			    clip_start, clip_end, (char) ((oumodiflg & REPORTMODIF_SOFTCLIP) != 0));
    cigarstr = ((errcode))? OUFMT_BAM_NULLSTR: diffStrGetViewStr(dvp);
    /* don't use asterisk here - bambam library can't handle that - use empty string */
    
    if (!errcode)
      editdist = diffStrGetLevenshteinDistance(diffstr);
  } else {
    cigarstr = OUFMT_BAM_NULLSTR;
  }
  
  errc = BamBam_BamWriter_PutAlignment(bamwriterp,
					bamflg,
					s_idx, /* chromosome */
					(uint64_t) ((pos > 0)? pos-1: 0), /* pos, zero based */
					ms_idx, /* mate chromosome */
					(uint64_t) ((mpos > 0)? mpos-1: 0),/* matepos, 0-based!*/
					readnamp->strp, seqstr, qualstr,
					cigarstr, /* cigar string */
					(int32_t) mapscor, isize);

  /* example aux tags */

  if (!errc)
    errc = BamBam_BamWriter_PutAuxNumber(bamwriterp,"NM",'i', &editdist);

  if (!errc)
    errc = BamBam_BamWriter_PutAuxNumber(bamwriterp,"AS",'i', &swatscor);
    //putBamAuxString(bwId,"XX","Test string data");

  /* write alignment */
  if (!errc)
    BamBam_BamWriter_Commit(bamwriterp);

  return (errc)? ERRCODE_WRITEERR: errcode;
}
#endif

static int writeREPALI(
		       FILE *fp,
#ifdef HAVE_BAMBAMC
		       BamBam_BamWriter *bamwriterp,
		       DiffView *dvp,
#endif
		       const REPALI *rp,
		       const DiffStr *rdfsp,
		       const SeqFastq *q_sqp,
#ifdef RESULTS_TRACKER
		       const Track *trackp,
		       BOOL_t *isHit,
#endif
		       const SeqSet *ssp,
		       const char *namext,
		       REPOUFMT_t outform,
		       REPMODIFLG_t oumodiflg,
		       DiffBlocks *dfblkp,
		       REPPAIRFLG_t pairflg,
		       int isize,
		       const REPALI *rsltmp,
		       SeqFastq *sqbufp,
		       REPNAMBUF *nambufp,
		       const SeqCodec *codecp
		       )
     /**< Print one result element 
      * \param fp Stream to which result should be printed.
      * \param rp Mapping result of the read.
      * \param rdfsp Compressed alignment string pointed to by rp->dfo.
      * \param q_sqp Sequence of the mapped read.
      * \param ssp Set of reference sequences.
      * \param namext Name extension of read (used to distinguish paired reads)
      * \param outform Output format, one of REPALI_OUTPUT_FORMATS.
      * \param pairflg Combination of RSLTPAIR_FLAGS.
      * \param mapflg Combination of MAP_FLAG bit flags.
      * \param is_2ndMate 0:rp and q_sqp refer to first read of a pair, 1: 2nd read of a pair
      *        (as determined by the order of input files
      * \param rsltmp Mapping result of the paired mate (can be NULL if single read mapping
      *        of if mate is not mapped.
      * \param matep Sequence of the mate (NULL if single read or mate not mapped).
      * \param sqbufp Buffer (used for SAM output).
      * \param codecp Sequence De-/Entcoder (used for SAM output)
      * 
      */
{ 
  int errcode = ERRCODE_SUCCESS;
  short mapscor =0;
  BOOL_t is_mapped = (BOOL_t) (rp != NULL && (rp->status & REPMATEFLG_MAPPED));
  const char *s_nam, *m_snam;
  DIFFSTR_T *dstrp;
  SEQLEN_t ref_len = 0;
#ifdef results_assert
  int len_prof, len_unprof;
#endif
#ifdef RESULTS_TRACKER
  int target_pos;
  SEQNUM_t target_refidx;
  TRACKFLG_t target_flags;
  const DiffStr *target_dfsp;
#endif

#ifdef report_debug
  printf("report_debug:fprintRESULT()\n");
#endif  

  if (is_mapped) {
    ref_len = seqSetGetSeqDatByIndex(NULL, &s_nam, rp->s_idx, ssp);
    dstrp = rdfsp->dstrp + rp->dfo;

#ifdef results_assert
    if ((errcode = diffStrCalcSeqLen(&len_prof, &len_unprof, dstrp)))
      return errcode;
    if (rp->q_start + len_prof != rp->q_end + 1 ||
	rp->s_start + len_unprof != rp->s_end + 1) {
      fprintf(stderr, "results_assert::fprintRESULT(): diffstr inconsistent with segment length.\n");
      return ERRCODE_DIFFSTR;
    }
#endif
  } else {
    ref_len = 0;
    s_nam = OUFMT_SAM_NULLSTR;
    dstrp = NULL;
 }
  if ((errcode = copyReadNamStrToREPSTR(&nambufp->ref_nam, 0, s_nam)))
    return errcode;

/* set flags */
  if (NULL != rp) {
    mapscor = rp->mapscor;
    if ((is_mapped) && (rsltmp) && rp->s_idx == rsltmp->s_idx) {
      pairflg |= REPPAIR_CONTIG;
    }
  }

  if (!namext) namext = "";
  switch (outform) {
  case REPORTFMT_CIGAR:
    errcode = fprintREPALIcigar(fp, rp, mapscor, &nambufp->q_nam, q_sqp, namext,
				nambufp->ref_nam.strp, dstrp, pairflg);
    break;
  case REPORTFMT_SSAHA:
    errcode = fprintREPALIssaha(fp, rp, mapscor, &nambufp->q_nam, q_sqp, namext,
				nambufp->ref_nam.strp, ref_len, dstrp, pairflg);
    break;
  case  REPORTFMT_GFF2:
    errcode = diffStrFindBlocks(dfblkp, dstrp);
    if (!errcode)
      errcode = fprintREPALIgff2(fp, rp, &nambufp->q_nam, q_sqp, namext,
				 nambufp->ref_nam.strp, dfblkp);
    break;
  case REPORTFMT_SAM:
    if (rsltmp) seqSetGetSeqDatByIndex(NULL, &m_snam, rsltmp->s_idx, ssp);
    else m_snam = OUFMT_SAM_NULLSTR;
    errcode = copyReadNamStrToREPSTR(&nambufp->mref_nam, 0, m_snam);
    if (!errcode)
      errcode = fprintREPALIsam(fp, sqbufp, 
				&nambufp->q_nam,
				mapscor, rp,
				dstrp,
				q_sqp, 
				nambufp->ref_nam.strp, 
				rsltmp, 
				nambufp->mref_nam.strp,
				isize, pairflg,
				oumodiflg, 
				codecp);
    break;
#ifdef HAVE_BAMBAMC
  case REPORTFMT_BAM:
    errcode = writeREPALIbam(bamwriterp, 
			     sqbufp, 
			     &nambufp->q_nam,
			     dvp,
			     mapscor, rp,
			     dstrp,
			     q_sqp, 
			     rsltmp,
			     isize, pairflg,
			     oumodiflg, 
			     codecp);
    break;
#endif
  default:
    errcode = ERRCODE_OUFMT;
    break;
  }
#ifdef RESULTS_TRACKER
  trackGetData(trackp, &target_pos, &target_refidx, &target_flags, NULL, &target_dfsp);
  
  if ((rp) && target_refidx == rp->s_idx) {
    int reflen;
    if ((errcode = diffStrCalcSeqLen(NULL, &reflen, target_dfsp->dstrp)))
      return errcode;
    if ((((SEQLEN_t) target_pos) <= rp->s_start) && 
	(rp->s_end < (SEQLEN_t) (target_pos + reflen))) {
      BOOL_t revflag = (rp->status&REPMATEFLG_REVERSE)? 1:0;
      if (revflag == ((target_flags & TRACKFLG_REVERSE)? 1:0))
	*isHit = 1;
    }
  }
#endif
#ifdef report_debug
  printf("report_debug:writeREPALI(): ");
  if (!errcode) 
    errcode = diffStrPrintf(stdout, dstrp, DIFFSTRFORM_RAW, 0, 0, 0);
  printf("\n");
#endif  
  return errcode;
}

/******************************************************************************
 *********************** Private Methods of Type REPPAIR **********************
 ******************************************************************************/

#define BLANK_REPPAIR(rpp) if ((rpp) != NULL) {		\
    (rpp)->pairflg = 0; (rpp)->isize = 0;		\
    (rpp)->iA = (rpp)->iB = -1;}			\


/******************************************************************************
 ****************** Private Methods Used by Type ReportWriter *****************
 ******************************************************************************/

static int writeSAMHeaderf(FILE *oufp, const SeqSet *ssp, 
		    const char *prognam, const char *progversion,
		    int narg, char * const *argv)
{
  SEQNUM_t s, snum = seqSetGetSeqNumAndTotLen(NULL, ssp);
  int i;
  char nambf[SEQNAM_SAM_MAXLEN];
 
  if (fprintf(oufp, SAMFORM_HEADLINE) < 1)
    return ERRCODE_WRITEERR;
  for (s=0; s<snum; s++) {
    const char *namp;
    SEQLEN_t sl = seqSetGetSeqDatByIndex(NULL, &namp, s, ssp);
    for (i=0; i<SEQNAM_SAM_MAXLEN-1 && namp[i] != '\0' && !isspace((int) namp[i]); i++)
      nambf[i] = namp[i];
    nambf[i] = '\0';

    if (fprintf(oufp, SAMFORM_REFSEQLINE, nambf, (unsigned int) sl) < 1)
      return ERRCODE_WRITEERR;
  }

  if (fprintf(oufp, SAMFORM_PROGLINE, prognam, prognam, progversion) < 0)
    return ERRCODE_WRITEERR;

  if (narg > 0 && argv != NULL) {
    if (fprintf(oufp, "%s", argv[0]) < 0)
       return ERRCODE_WRITEERR;
    for (i=1; i<narg; i++)
      if (fprintf(oufp, " %s", argv[i]) < 0)
	return ERRCODE_WRITEERR;
    if (fprintf(oufp, "\n") < 0)
      return ERRCODE_WRITEERR;
  }
  return ERRCODE_SUCCESS;
}

#ifdef HAVE_BAMBAMC
static int initBAMHeader(BamBam_BamWriter **bamwriterpp, 
			 BamBam_BamHeaderInfo **headinfopp, 
			 REPSTR * sbufp,
			 const char *prognam, 
			 const char *progversion,
			 int narg, 
			 char * const *argv,
			 const char *filnamp,
			 const SeqSet *ssp)
{
  int errcode = ERRCODE_SUCCESS;
  SEQNUM_t s, nseq;

  if (NULL == ssp)
    return ERRCODE_ASSERT;
  
  if ((errcode = copySAMheaderCommandLineToREPSTR(sbufp, 
						  prognam, progversion,
						  narg, argv)))
    return errcode;
 
  
  *headinfopp = BamBam_BamHeaderInfo_New(SAMBAM_HEADER_VERSION, 
					 SAMBAM_SORTORDER_UNKNOWN, 
					 sbufp->strp);

  if (NULL == *headinfopp)
    return ERRCODE_NOMEM;

  nseq = seqSetGetOffsets(ssp, NULL);
  for (s=0; s<nseq; s++) {
    const char *chrnam;
    const SEQLEN_t chrlen = seqSetGetSeqDatByIndex(NULL, &chrnam, s, ssp);

    if ((errcode = copyReadNamStrToREPSTR(sbufp, 0, chrnam)))
      break;

    if (BamBam_BamHeaderInfo_AddChromosome(*headinfopp, sbufp->strp,
					   chrlen)) {
      errcode = ERRCODE_BAMBAM;
      break;
    }
  }
  if (errcode != ERRCODE_SUCCESS)
    return errcode;

  *bamwriterpp = BamBam_BamWriter_New(*headinfopp, 
				      filnamp, 
				      SAMBAM_COMPRESSION_LEVEL);

  return (NULL == *bamwriterpp)? ERRCODE_NOMEM: ERRCODE_SUCCESS;
}
#endif

/******************************************************************************
 ********************* Public Methods of Type ReportWriter ********************
 ******************************************************************************/

ReportWriter *reportCreateWriter(int *errcode,
				 const char * const filnam,
				 const REPOUFMT_t outform, 
				 const REPMODIFLG_t modiflg,
				 const SeqSet *ssp,
				 const char *prognam,
				 const char *progversion,
				 char * const *cmdlin_argv,
				 int cmdlin_narg)
{
  int errc = ERRCODE_SUCCESS;
  ReportWriter *p;
  EMALLOCP0(p);

  if (NULL == p) {
    *errcode = ERRCODE_NOMEM;
    return p;
  }
  if (NULL == filnam) {
    ESTRCPY(p->filnam, OUFMT_FILNAM_STDOUT);
  } else {
    ESTRCPY(p->filnam, filnam);
  }
  if (NULL == p->filnam ||
      (NULL == (p->qbufp = seqFastqCreate(0, SEQTYP_FASTQ))) ||
      (NULL == (p->sbufp = seqFastqCreate(0, SEQTYP_FASTA))) ||
      (NULL == (p->nambufp = createREPNAMBUF())))
    errc = ERRCODE_NOMEM;
  else if (REPORTFMT_GFF2 == outform) {
    p->dfblkp = diffBlocksCreate(0);
    if (NULL == p->dfblkp)
      errc = ERRCODE_NOMEM;
    else 
      p->dfblkp = NULL;
  }

  p->oufmt = outform;
  p->modflg = modiflg;
  p->linwidth = DEFAULT_LINWIDTH_ALI;
  p->namext[0] = '\0';
  p->namext_mate[0] = '\0';
  p->dfblkp = NULL;

  if (ERRCODE_SUCCESS == errc) {   
#ifdef HAVE_BAMBAMC
    if (REPORTFMT_BAM == outform) {
      p->oufp = stdout; /* needed e.g. for explicit alignment output */
      errc = initBAMHeader(&p->bamwriterp, &p->headinfop, 
			   &p->nambufp->ref_nam, 
			   prognam, progversion,
			   cmdlin_narg, cmdlin_argv,
			   p->filnam, 
			   ssp);
     if (!(errc) && 
	  ((p->dvp = diffStrCreateView(0)) == NULL))
	errc = ERRCODE_NOMEM;
    } else {
      p->headinfop = NULL;
      p->bamwriterp = NULL;
      p->dvp = NULL;
#endif
      if (filnam == NULL) {
	p->oufp = stdout;
      } else {
	p->oufp = EFOPEN(p->filnam, "w");
	if (NULL == p->oufp) {
	  errc = ERRCODE_NOFILE;
	} 
      }
      if ((outform == REPORTFMT_SAM) && (modiflg & REPORTMODIF_HEADER)) {
	if (NULL == ssp) {
	  errc = ERRCODE_ASSERT;
	} else {
	  errc = writeSAMHeaderf(p->oufp, ssp,
				 prognam, progversion,
				 cmdlin_narg, cmdlin_argv);
	}
      }
#ifdef HAVE_BAMBAMC
    }
#endif
  }
 
  if (ERRCODE_SUCCESS != errc) {
    reportDeleteWriter(p);
    p = NULL;
  }

  *errcode = errc;

  return p;
}

void reportDeleteWriter(ReportWriter *p)
{
#ifdef HAVE_BAMBAMC
  int deleteStatus = 0;
#endif
  if (p != NULL) {
    free(p->filnam);
    if (p->oufp != NULL && p->oufp != stdout)
      EFCLOSE(p->oufp);
    diffBlocksDelete(p->dfblkp);
    seqFastqDelete(p->sbufp);
    seqFastqDelete(p->qbufp);
    deleteREPNAMBUF(p->nambufp);
#ifdef HAVE_BAMBAMC
    diffStrDeleteView(p->dvp);
    if (p->bamwriterp != NULL)
      BamBam_BamWriter_Delete(p->bamwriterp, &deleteStatus);
    if (p->headinfop != NULL)
      BamBam_BamHeaderInfo_Delete(p->headinfop);
#endif
  }
  free(p);
}

FILE *reportGetWriterStream(const ReportWriter *p)
{
  return p->oufp;
}
/******************************************************************************
 ************************ Private Methods of Type Report **********************
 ******************************************************************************/

static int writeReportForRead(const ReportWriter *wrp,
#ifdef RESULTS_TRACKER
			      const Track *trackp,
#endif
			      const REPALI *ralip,
			      const SeqFastq *readp,
			      const DiffStr *rdfsp,
			      const REPALI *malip,
			      int isize,
			      REPPAIRFLG_t pairflg,
			      const SeqSet *ssp,
			      const SeqCodec *codecp)
{
  int errcode;

#ifdef RESULTS_TRACKER
  BOOL_t isHitRead = 0;
#endif 

  if ((errcode = writeREPALI(wrp->oufp, 
#ifdef HAVE_BAMBAMC
			     wrp->bamwriterp,
			     wrp->dvp,
#endif
			     ralip, rdfsp, readp,
#ifdef RESULTS_TRACKER
			     trackp,
			     &isHitRead,
#endif 
			      ssp,
			      wrp->namext, wrp->oufmt, wrp->modflg, wrp->dfblkp,
			      pairflg, isize, malip, wrp->qbufp, wrp->nambufp,
			      codecp)))
    return errcode;

  if ((wrp->modflg & REPORTMODIF_ALIOUT) != 0 && (ralip) &&
      (ralip->status & REPMATEFLG_MAPPED) != 0 &&
      (errcode = fprintAlignment(wrp->oufp, wrp->sbufp, wrp->qbufp, 
				 readp, wrp->linwidth,
				 rdfsp, ssp, codecp, ralip)))
    return errcode;

#ifdef RESULTS_TRACKER
  if (!isHitRead)
    trackPrintf(wrp->oufp, trackp);
#endif 

  return errcode;
}
/******************************************************************************
 ************************ Private Methods of Type Report **********************
 ******************************************************************************/

/******************************************************************************
 ************************* Public Methods of Type Report **********************
 ******************************************************************************/

Report *reportCreate(int blksz)
{
  Report *p;
  
  if (NULL != EMALLOCP0(p)) {
    if (blksz <= 0) blksz = REPORT_DEFBLKSZ;
    if (NULL == ARRCREATE(p->pairr, blksz) ||
	NULL == ARRCREATE(p->arAr, blksz) ||
	NULL == ARRCREATE(p->arBr, blksz) ||
	diffStrInit(&p->dfs, 0)) {
      reportDelete(p);
      p = NULL;
    } 
  }
  return p;
}

void reportDelete(Report *p)
{
  if (p) {
    ARRDELETE(p->pairr);
    ARRDELETE(p->arAr);
    ARRDELETE(p->arBr);
    diffStrCleanUp(&p->dfs);
  }
  free(p);
}

void reportBlank(Report *p)
{
  if (p) {
    ARRLEN(p->pairr) = 0;
    ARRLEN(p->arAr) = 0;
    ARRLEN(p->arBr) = 0;
    DIFFSTR_LENGTH(&p->dfs) = 0;
  }
}

int reportNextPairID(Report *rep)
{
  int pairid = ARRLEN(rep->pairr);
  REPPAIR *pp;

  ARRNEXTP(pp, rep->pairr);
  if (NULL == pp)
    pairid = -1;
  else {
    BLANK_REPPAIR(rep->pairr + pairid);
  }

  return pairid;
}

int reportAddMap(Report *rep, 
		 int pairid,
		 int swatscor, short mapscor,
		 SEQLEN_t q_start, SEQLEN_t q_end,
		 SEQLEN_t s_start, SEQLEN_t s_end, SEQNUM_t s_idx,
		 const DIFFSTR_T *dstrp, int dfslen,
		 int insiz,
		 REPMATEFLG_t mateflg, REPPAIRFLG_t pairflg)
{
  int errcode = ERRCODE_SUCCESS, errc = ERRCODE_SUCCESS;
  int idx = -1;
  REPALI *rp = NULL;
  REPPAIR *pp = NULL;
  
  if (NULL == dstrp || dfslen < 1) {
    mateflg &= ~REPMATEFLG_MAPPED;
  }

  

  if ((mateflg & REPMATEFLG_PAIRED) && pairid >= 0) {
    if (((size_t) pairid) >= ARRLEN(rep->pairr))
      return ERRCODE_ARGRANGE;
    pp = rep->pairr + pairid;
    if (pp->pairflg == 0) { /* no pair info, yet */
      pp->pairflg = pairflg;
    } else if (pp->pairflg != pairflg) {
      return ERRCODE_ASSERT;
    }
  }

  if (pp != NULL &&
      (mateflg & REPMATEFLG_2NDMATE)) {
    if (pp->iA >= 0) { 
      /* 1st mate already set, so check insert size is consistent */
      if (insiz != pp->isize)
	return ERRCODE_ASSERT;
      errc = findREPALI(rep->arBr, &idx, q_start, q_end, mateflg,
			s_start, s_end, s_idx);
      if (ERRCODE_FAILURE == errc) {
	pp->iB = ARRLEN(rep->arBr);
	ARRNEXTP(rp, rep->arBr)
	if (NULL == rp)
	  errcode = ERRCODE_NOMEM;
      } else if (ERRCODE_SUCCESS == errc) {
	pp->iB = idx;
	rp = rep->arBr + idx;;
      } else {
	errcode = errc;
      }
    } else {
      pp->isize = insiz;
    }
  } else {
    /* here, the 1st mate of a paired read is added or
       a single read that is not part of a pair */
    REPALI **arp = &rep->arAr;
    if (NULL == pp) { 
      /* add single read */
      if (mateflg & REPMATEFLG_2NDMATE)
	arp = &rep->arBr;
    } else {
      /* add 1st mate of a read pair */
      if (pp->iB >= 0) {
	/* 2nd mate already set, so check insert size is consistent */
	if (insiz != pp->isize)
	  return ERRCODE_ASSERT;
      } else {
	pp->isize = insiz;
      }
    }
    errc = findREPALI(*arp, &idx, q_start, q_end, mateflg,
		      s_start, s_end, s_idx);
    if (ERRCODE_FAILURE == errc) {
      if (pp != NULL) 
	pp->iA = ARRLEN(rep->arAr);
      ARRNEXTP(rp, *arp);
      if (NULL == rp) 
	errcode = ERRCODE_NOMEM;
    } else if (ERRCODE_SUCCESS == errc) {
      if (NULL == pp) {
	/* single read/partial mapping is already known -> ignore */
	rp = NULL;
      } else {
	/* 1st mate of a pair, overwrite */
	pp->iA = idx;
	rp = *arp + idx;
      }
    } else {
      errcode = errc;
    }
  }
 
  if (!(errcode) && rp != NULL) {
    rp->status = mateflg;
    rp->dfo = DIFFSTR_LENGTH(&rep->dfs);
    if ((mateflg & REPMATEFLG_MAPPED)) {
      rp->swatscor = swatscor;
      rp->mapscor = mapscor;
      rp->q_start = q_start;
      rp->q_end = q_end;
      rp->s_start = s_start;
      rp->s_end = s_end;
      rp->s_idx = s_idx;
      errcode = diffStrAdd(&rep->dfs, dstrp, dfslen);
    } else {
      rp->swatscor = 0;
      rp->mapscor = 0;
      rp->q_start = 0;
      rp->q_end = 0;
      rp->s_start = 0;
      rp->s_end = 0;
      rp->s_idx = 0;
    }
#ifdef report_debug
  fprintf(stdout, "report_debug::reportAddMap():");
  fprintREPALIraw(stdout, rp);
#endif
  }

  return errcode;
}

void reportFixMultiplePrimary(Report *rep)
{
  int n;
  int np = ARRLEN(rep->pairr);
  int na = ARRLEN(rep->arAr);
  int nb = ARRLEN(rep->arBr);
  int n_primary_A = 0;
  int n_primary_B = 0;

  for (n=0; n<np && (n_primary_A < 2 || n_primary_B < 2); n++) {
     REPPAIR *pp = rep->pairr + n;
     if (rep->arAr[pp->iA].status & REPMATEFLG_PRIMARY)
       n_primary_A++;
     if (rep->arAr[pp->iB].status & REPMATEFLG_PRIMARY)
       n_primary_B++;
  }
  if ( n_primary_A < 2 ) {
    if (n_primary_A > 0) n_primary_A = 0;
    for (n=0; n < na && n_primary_A < 2; n++)
      if (rep->arAr[n].status & REPMATEFLG_PRIMARY)
	n_primary_A++;
  }
  if ( n_primary_B < 2 ) {
    if (n_primary_B > 0) n_primary_B = 0;
    for (n=0; n < nb && n_primary_B < 2; n++)
      if (rep->arBr[n].status & REPMATEFLG_PRIMARY)
	n_primary_B++;
  }
  if ( n_primary_A > 1)
    for (n=0; n < na; n++) 
      rep->arAr[n].status &= ~REPMATEFLG_PRIMARY;
  
   if ( n_primary_B > 1)
    for (n=0; n < nb; n++) 
      rep->arBr[n].status &= ~REPMATEFLG_PRIMARY;
   
   return;   
}

int reportWrite(const ReportWriter *wrp,
#ifdef RESULTS_TRACKER
		const Track *read_trackp,
		const Track *mate_trackp,
#endif 
		const SeqFastq *readp, 
		const SeqFastq *matep, 
		const SeqSet *ssp,
		const SeqCodec *codecp,
		const Report *rep)
{
  int errcode = ERRCODE_SUCCESS;
  int n, na = ARRLEN(rep->arAr);
  int nb = ARRLEN(rep->arBr);
  int np = ARRLEN(rep->pairr);
  REPPAIRFLG_t pairflg;
  REPALI *ap = NULL, *bp = NULL;
  
  /* print all the pairs first */
  for (n=0; n<na; n++)
    rep->arAr[n].was_output = 0;
  for (n=0; n<nb; n++)
    rep->arBr[n].was_output = 0;
  
  for (n=0; n<np; n++) {
    REPPAIR *pp = rep->pairr + n;
    ap = rep->arAr + pp->iA;
    bp = rep->arBr + pp->iB;
    ap->was_output = 1;
    bp->was_output = 1;
    errcode = writeReportForRead(wrp, 
#ifdef RESULTS_TRACKER
				 read_trackp,
#endif				    
				 ap,readp,
				 &rep->dfs,
				 bp,
				 pp->isize,
				 pp->pairflg,
				 ssp, codecp);
    
    if ((errcode)) break;
    
    errcode = writeReportForRead(wrp, 
#ifdef RESULTS_TRACKER
				 mate_trackp,
#endif				    
				 bp,matep,
				 &rep->dfs,
				 ap,
				 pp->isize,
				 pp->pairflg,
				 ssp, codecp);
    if ((errcode))
      break;
  }
  
  if ((errcode))
    return errcode;

  /* ouptput remaining alignments, not yet printed */
  if (n > 0) {
    pairflg = rep->pairr[0].pairflg;
  } else {
   pairflg = 0;
  }
  
  for (n=0; n < na && !(errcode); n++) {
    ap = rep->arAr + n;
    if ((rep->arAr[n].was_output))
      continue;
    errcode = writeReportForRead(wrp, 
#ifdef RESULTS_TRACKER
				 read_trackp,
#endif				    
				 ap,readp,
				 &rep->dfs,
				 NULL,
				 0, 
				 pairflg,
				 ssp, codecp);
  }


  if ((errcode))
    return errcode;


  for (n=0; n < nb && !(errcode); n++) {
    bp = rep->arBr + n;
    if ((rep->arBr[n].was_output))
      continue;
    errcode = writeReportForRead(wrp, 
#ifdef RESULTS_TRACKER

				 mate_trackp,
#endif				    
				 bp, matep,
				 &rep->dfs,
				 NULL,
				 0,
				 pairflg,
				 ssp, codecp);
  }

  return errcode;
}
