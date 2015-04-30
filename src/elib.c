/** Stream/memory management with error messages */

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

#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include"elib.h"

enum ESTRING_CONST {
  ESTRING_DEFAULT_BLKSZ = 128,
};

enum ERRMSG_CONST {
  ERRMSG_EMPTY = 0,     /**< Empty message */
  ERRMSG_UNKNOWN = 1,   /**< Unkown message level */
  ERRMSG_NLEVELS = 5,   /**< number of message levels */
  ERRMSG_MAXNUM  = 10,  /**< maximum number of messages allowed before flush */
  ERRMSG_MAXLEN  = 50,  /**< maximum length of message text */
  ERRMSG_FILELEN = 32,  /**< for file name */
  ERRMSG_MAXNAMLEN = 64,/**< Maximum length of read name */
  ERRMSG_MAXSPACE = 512,/**< Total number of characters allowed */
  ERRMSG_MINSPACE = 32, /**< Less space than that triggers flush of error messages */
};

typedef struct _MESSAGE {
  char *textp;    /**< points to message text (can be NULL) */
  char *filnamp;  /**< points to source file name (can be NULL) */
  char *readnamp; /**< Name of read for which this message was posted */
  size_t readno;  /**< Number of ther read which this message was posted */
  int linno;      /**< line number in source code */
  short num;      /**< Message number */
  int code;       /**< error code (one of ERRMSG_CODES) */
  short level;    /**< Message level, one of ERRMSG_LEVELS (can be ERRMSG_EMPTY) */
} MESSAGE;

struct _ErrMsg {
  MESSAGE messages[ERRMSG_MAXNUM];    /**< Messages */
  char mspace[ERRMSG_MAXSPACE];       /**< Memory for text */
  char currReadNam[ERRMSG_MAXNAMLEN]; /**< Current read name */
  size_t currReadNo;                  /**< Current read number (1-based) */
  short num;                          /**< Number of messages in stack */
  short usedSpace;                    /**< Number of bytes used in mspace */
};

const char ERRMSG_FORMAT[] = "[%i] %s:%d %s: %s %s\n";

enum EF_CONST {
  EF_PIPECHAR = '-'
};

#define ERRMSG_STRCPY(top, fromp, emp)\
  maxfree = ERRMSG_MAXSPACE - (emp)->usedSpace - 1;\
  for (n=0; ((fromp)[n]) && (n < maxfree); n++) (top)[n] = (fromp)[n];\
  (top)[n++] = '\0';\
  emp->usedSpace += n;

const char *errMsgString(int errcode) 
{
  switch(errcode) {
  case ERRCODE_SUCCESS:
    return "success";
  case ERRCODE_FAILURE:
    return "failure";
  case ERRCODE_NOMEM:
    return "memory allocation failed";
  case ERRCODE_PREEOF:
    return "unexpected end of file";
  case ERRCODE_SEQIOTYP:
    return "unknown SeqIO type";
  case ERRCODE_FILTYP:
    return "wrong file type";
  case ERRCODE_FASTA:
    return "wrong FASTQ/FASTA format";
  case ERRCODE_NOPROMPT:
    return "wrong or missing FASTQ prompt";
  case ERRCODE_SEQIOMODE:
    return "wrong I/O mode";
  case ERRCODE_NOFILE:
    return "could not open file";
  case ERRCODE_EOF:
    return "end of file";
  case ERRCODE_WRITEERR:
    return "could not write to file";
  case ERRCODE_NOWRITE:
    return "attempt to write to file opened for reading";
  case ERRCODE_NOREAD:
    return "attempt to read from file opened for writing";
  case ERRCODE_READERR:
    return "could not read from file";
  case ERRCODE_NULLPTR:
    return "unexpected NULL pointer";
  case ERRCODE_MAXKTUP:
    return "hashed word length too long";
  case ERRCODE_MAXKPOS:
    return "maximum number of k-mer positions exceeded";
  case ERRCODE_HASHSEQTYP:
    return "sequence type cannot be hashed";
  case ERRCODE_HASH_NOMORESEQ:
    return "cannot add sequence to hash table because table is already set up";
  case ERRCODE_ARGINVAL:
    return "invalid function argument";
  case ERRCODE_BROKENHASH:
    return "broken hash table";
  case ERRCODE_SEQCODE:
    return "wrong type of sequence encoding";
  case ERRCODE_HASHHITSORT:
    return "list of hits not yet sorted - sequence indices not assigned";
  case ERRCODE_INTRANGE:
    return "integer outside the range of integer type";
  case ERRCODE_HASHREAD:
    return "Read hash files are inconsistent";
  case ERRCODE_SEQNOSTR:
    return "expected sequences represented as character strings";
  case ERRCODE_IOBUFF:
    return "couldn't set I/O buffer size";
  case ERRCODE_COMPRESS:
    return "inconsistency resulting from sequence compression";
  case ERRCODE_ARGRANGE:
    return "arguments outside expected range";
  case ERRCODE_SHORTSEQ:
    return "sequence too short to be hashed";
  case ERRCODE_COMMANDLINE:
    return "command line format error";
  case ERRCODE_ALLOCBOUNDARY:
    return "hit boundary of allocated memory";
  case ERRCODE_FILEFORM:
    return "error in file format";
  case ERRCODE_SORTSTACK:
    return "sort stack overflow";
  case ERRCODE_HITSTATUS:
    return "hit list does not have required status";
  case ERRCODE_SEQLEN:
    return "sequence too long";
  case ERRCODE_NOSEQ:
    return "no sequence to be hashed";
  case ERRCODE_QUALLEN:
    return "quality factor and base sequence lengths differ";
  case ERRCODE_SEQTYP:
    return "wrong sequence type";
  case ERRCODE_FREQSUMS:
    return "rounding error when calculating frequency sums";
  case ERRCODE_SWATEXCEED:
    return "Smith-Waterman score overflow";
  case ERRCODE_ALISCOREXCEED:
    return "Alignment score exceeds target while back tracking";
  case ERRCODE_DIFFCOUNT:
    return "Mutation count inconsistent in diff string";
  case ERRCODE_SWATSCOR:
    return "Inconsistency when calculating Smith-Waterman scores";
  case ERRCODE_HTLSTQLEN:
    return "inconsistent query length between forward and reverse complement hit lists";
  case ERRCODE_MATCHSTAT:
    return "Match list does not have correct status";
  case ERRCODE_ASSERT:
    return "assertion failed";
  case ERRCODE_OVERFLOW:
    return "integer overflow";
  case ERRCODE_FILPOS:
    return "setting file position when error occured";
  case ERRCODE_FILEIO:
    return "file input/output error";
  case ERRCODE_MENUE:
    return "error when parsing command line";   
  case ERRCODE_ENDIAN:
    return "file written under different endianess - conversion failed";
  case ERRCODE_FASTAPAIRNUM:
    return "The two FASTA/FASTQ input file have different numbers of reads";
  case ERRCODE_LONGFILNAM:
    return "file name too long";
  case ERRCODE_FHEADSIZ:
    return "wrong size of file type specific header in binary file";
  case ERRCODE_FILVERSION:
    return "wrong version of binary file format";
  case ERRCODE_LONGSEQ:
    return "sequence too long";
  case ERRCODE_NOSEQIDX:
    return "sequence index missing";
  case ERRCODE_DIFFSTR:
    return "inconsistent diff-string";
  case ERRCODE_SEQOFFS:
    return "sequence offset greater than sequence length";
  case ERRCODE_TERMIT:
    return "termination of iterative calls";
  case ERRCODE_CPLXSCOR:
    return "complexity weighted score exceeds unweighted score";
  case ERRCODE_PTHREAD:
    return "thread returned error";
  case ERRCODE_HITINFO:
    return "HashHitInfo has wrong status";
  case ERRCODE_HITBIN:
    return "K-mer statistics insufficient for binning";
  case ERRCODE_SWATSTRIP:
    return "inconsistency in striping format for Smith-Waterman using SSE2";
  case ERRCODE_QUALVAL:
    return "invalid FASTQ base quality value";
  case ERRCODE_SEGPARTFETCH:
    return "fetched only part of segment requested from sequence set";
  case ERRCODE_NOHASHWORD:
    return "word on not found in hash index";
  case ERRCODE_PTHREADSTACK:
    return "thread stack overflow";
  case ERRCODE_FILSIZ:
    return "file size too big";
  case ERRCODE_SKIPSMALL:
    return "sampling step size to small for total length of reference sequences";
  case ERRCODE_SEQSETSIZ:
    return "total number of nucleotides in sequence set exceeds limit";
  case ERRCODE_PAIRNUM:
    return "number of pairs exceeds limit";
  case ERRCODE_NOMATCH:
    return "Found no match";
  case ERRCODE_RNAMPAIR:
    return "Read names don't match as expected for a read pair";
  case ERRCODE_MATENUMPAIR:
    return "Number of reads differs in the two input files for paired reads";
  case ERRCODE_BREAK:
    return "Early termination triggered by debugging code";
  case ERRCODE_PTHRTERMSIG:
    return "Termination signal for thread";
  case ERRCODE_PTHRPULLCHK:
    return "could not pull thread argument from buffer (check not fulfilled)";
  case ERRCODE_SEQNAMLEN:    
    return "Sequence name too long";
  case ERRCODE_KPOSOFLO:
    return "Number of positions in hash index caused integer overflow";
  case ERRCODE_INFMT:
    return "Unrecognised sequence input format";
  case ERRCODE_OUFMT:
    return "Unrecognised sequence output format";    
  case ERRCODE_SEQNTSMBL:
    return "Nucleotide symbol is not a letter";
  case ERRCODE_SEQNUMSET:
    return "Number of sequences exceeds limit";
  case ERRCODE_BAMBAM:
    return "bambamc library returns error";
  case ERRCODE_SEMOPEN:
    return "named semaphore could not be opened";
  default:
    return "unknown error code";
  }
}

static const char *strMessageLevel(unsigned char level) 
{
  switch(level) {
  case ERRLEVEL_FATAL:
    return "ERROR";
  case ERRLEVEL_WARNING:
    return "WARNING";
  case ERRLEVEL_REPORT:
    return "REPORT";
  default:
    return "UNKNOWN";
  }
}

static void fprintReadNameAndNumber(FILE *fp, const char *rnam, size_t rno)
{
  if (((rnam) && (rnam[0])) ||
      (rno > 1)) {
    fprintf(fp, "  when processing read");
    if (rno > 1) 
      fprintf(fp, " No. %llu", (long long unsigned) rno);
    if ((rnam) && (rnam[0]))
      fprintf(fp, " '%s'", rnam);
    fprintf(fp, "\n");
  }
  return;
}

ErrMsg *errMsgCreate(const char *progfil, int linenum) 
{
  int i;
  ErrMsg *p;

  p = (ErrMsg *) emalloc(sizeof(ErrMsg), progfil, linenum);
  if (!p) 
    exit(EXIT_FAILURE);

  p->num = 0;
  for (i=0; i<ERRMSG_MAXNUM; i++)
    p->messages[i].level = ERRMSG_EMPTY;
  p->currReadNam[0] = '\0';
  p->currReadNo = 0;
  p->usedSpace = 0;
  p->mspace[0] = '\0';
  return p;
}

void errMsgEnd(ErrMsg *p)
{
  if ((p) && p->num > 0) {
    errMsgFlush(stdout, p);
  } 
  free(p);
}

void errMsgSetCurrentReadName(ErrMsg *emp, const char *namp)
{
  strncpy(emp->currReadNam, namp, ERRMSG_MAXNAMLEN);
  emp->currReadNam[ERRMSG_MAXNAMLEN-1] = '\0';
}

void errMsgSetCurrentReadNumber(ErrMsg *emp, size_t rno)
{
  emp->currReadNo = rno;
}

/** Add an error message to the stack */
void errMsgAdd(ErrMsg *emp, const char *message, 
	       const char *progfil, int linenum,
	       int errcode, unsigned char level)
{
  short maxfree, n;
  MESSAGE *mp;

  if (level == ERRLEVEL_FATAL || !(emp)) {
    fprintf(stderr, ERRMSG_FORMAT, 
	    (emp)? emp->num: 0,
	    progfil, linenum, 
	    strMessageLevel(ERRLEVEL_FATAL),
	    errMsgString(errcode),
	    (message)? message: ""); 
    if (emp) {
      fprintReadNameAndNumber(stderr, emp->currReadNam, emp->currReadNo);
      errMsgFlush(stderr, emp);
      exit(EXIT_FAILURE);
    }
  } else {
    mp = emp->messages + emp->num;
    mp->code = errcode;
    mp->linno = linenum;
    mp->num = emp->num;
    mp->level = (level < ERRMSG_NLEVELS)?
      level: ERRMSG_UNKNOWN;

    if (emp->usedSpace < ERRMSG_MAXSPACE) {
      mp->filnamp = emp->mspace + emp->usedSpace;
      ERRMSG_STRCPY(mp->filnamp, progfil, emp);
    } else {
      mp->filnamp = NULL;
    }

    if (emp->usedSpace < ERRMSG_MAXSPACE) {
      mp->textp = emp->mspace + emp->usedSpace;
      ERRMSG_STRCPY(mp->textp, message, emp);
    } else {
      mp->textp = NULL;
    }

    if (emp->usedSpace < ERRMSG_MAXSPACE) {
      mp->readnamp = emp->mspace + emp->usedSpace;
      ERRMSG_STRCPY(mp->readnamp, emp->currReadNam, emp);
    } else {
      mp->readnamp = NULL;
    }
  
    emp->num++;

    if (emp->num >= ERRMSG_MAXNUM || 
	(emp->usedSpace + ERRMSG_MINSPACE > ERRMSG_MAXSPACE)) {
      errMsgFlush(stderr, emp);
      exit(EXIT_FAILURE);
    }
  }
}


/* print Messages on stream */
void errMsgFlush(FILE *fp, ErrMsg *emp)
{
  if ((emp)) {
    MESSAGE *mp;
    while(emp->num--) {
      mp = emp->messages + emp->num;
      fprintf((fp), ERRMSG_FORMAT, 
	      mp->num, mp->filnamp, mp->linno,
	      strMessageLevel(mp->level),
	      errMsgString(mp->code),
	      (mp->textp)? mp->textp: "");
      fprintReadNameAndNumber(fp, mp->readnamp, mp->readno);
    }
  }  
  return;
}
    
void *emalloc(size_t size, const char *progfil, int line)
{
  void *p;
#ifdef elib_debug
  fprintf(stderr, "elib_debug::%s:%i, (malloc) Allocating size = %llu ...\n", 
	 progfil, line, (unsigned long long) size);
  fflush(stderr);
#endif
  if (!(p = malloc(size)))
    fprintf(stderr, "ERROR: malloc(%lu) failed in %s, line %i\n",
	    (unsigned long) size, progfil, line);
  return p;
}

void *ecalloc(size_t nobj, size_t size, const char *progfil, int line)
{
  void *p;
#ifdef elib_debug
  fprintf(stderr, "elib_debug::%s:%i, (ecalloc) Allocating %llu objects of size = %llu ...\n", 
	 progfil, line, (unsigned long long) nobj, (unsigned long long) size);
  fflush(stderr);
#endif
  if (!(p = calloc(nobj, size)))
    fprintf(stderr, "ERROR: calloc(%lu, %lu) failed in %s, line %i\n",
	    (unsigned long) nobj, (unsigned long) size, progfil, line);
  return p;
}

void *erealloc(void *objp, size_t size, size_t size_old, const char *progfil, int line)
     /* initialize newly allocated space if (size_old > 0) */
{
  void *p;
#ifdef elib_debug
  fprintf(stderr, "elib_debug::%s:%i, Reallocating size = %llu ...\n", 
	 progfil, line, (unsigned long long) size);
  fflush(stderr);
#endif
  if (size == size_old) {
    p = objp;
  } else {
    p = realloc(objp, size);
    if (p == NULL) {
      fprintf(stderr, "ERROR: realloc(%p, %lu) failed in %s, line %i\n",
	      objp, (unsigned long) size, progfil, line);
    } else if (size_old > 0 && size > size_old) {
      memset((char *) p + size_old, 0, size - size_old);
    }
  }
  return p;
}

void *ereallocp(void **objp, size_t size, const char *progfil, int line)
{
  void *p; 
#ifdef elib_debug
  fprintf(stderr, "elib_debug::%s:%i, ereallocp size = %llu ...\n", 
	 progfil, line, (unsigned long long) size);
  fflush(stderr);
#endif
  p = realloc(*objp, size);
  if (p == NULL) {
    fprintf(stderr, "ERROR: realloc(%p, %lu) failed in %s, line %i\n",
	    *objp, (unsigned long) size, progfil, line);
  } else {
    *objp = p;
  }

  return p;
}


FILE *efopen(const char *filnam, const char *mode, 
	     const char *progfil, int line)
     /* fopen with error message */
{
  FILE *fp = NULL;

  if (filnam[0] == EF_PIPECHAR && filnam[1] == '\0') {
    if (mode[0] == 'w')
      fp = stdout;
    else if (mode[0] == 'r')
      fp = stdin;
  } else {
    fp = fopen(filnam, mode);
  }
  if (NULL == fp) {
    fprintf(stderr, "ERROR: fopen (%s, \"%s\") failed in %s, line %i\n",
	    filnam, mode, progfil, line);
  }

  return fp;
}

int efclose(FILE *fp, const char *progfil, int line) 
     /* fclose with error checking */
{
  int rv = ERRCODE_SUCCESS;
  if ((fp)) {
    if (ferror(fp)) {
      fprintf(stderr, "ERROR when closing file in %s, line %i\n", progfil, line);
      perror("Error message from <stdlib>:\n");
      rv = ERRCODE_FILEIO;
    }
    if (fp != stdout && fp != stderr && fclose(fp))
      rv = ERRCODE_FILEIO;
  }
  return rv;
}

char *estrcpy(const char *str, const char *progfil, int line)
{
  char *cp = malloc(strlen(str)+1);
  if (cp) strcpy(cp, str);
  else fprintf(stderr, "ERROR: estrcpy() failed in %s, line %i\n",
	       progfil, line);
  return cp;
}

char *estrcat(const char *str1, const char *str2, const char *progfil, int line)
{
  char *cp = malloc(strlen(str1) + strlen(str2) + 1);
  if (cp) {
    strcpy(cp, str1);
    strcat(cp, str2);
  } else {
    fprintf(stderr, "ERROR: estrcat() failed in %s, line %i\n",
	    progfil, line);
  }
  return cp;
}

/*****************************************************************************
 ********************* Private Methods of Type EString ***********************
 *****************************************************************************/

/*****************************************************************************
 ********************** Public Methods of Type EString ***********************
 *****************************************************************************/
char *eStringInit(EString *esp, int blksz, const char *progfil, int line)
{
  if (NULL == esp)
    return NULL;

  esp->len = 0;
  esp->allocsz = 0;
  esp->blksz = (blksz < 1)? ESTRING_DEFAULT_BLKSZ: blksz;
  esp->strp = emalloc(esp->blksz, progfil, line);
  if (esp->strp != NULL) {
    esp->strp[0] = '\0';
    esp->allocsz = (size_t) blksz;
  } 
  return esp->strp;
}

EString *eStringCreate(int blksz, const char *progfil, int line)
{
  EString *esp = emalloc(sizeof(EString), progfil, line);

  if (esp != NULL && 
      NULL == eStringInit(esp, blksz, progfil, line)) {
    eStringDelete(&esp, progfil, line);
    esp = NULL;
  }

  return esp;
}

void eStringDelete(EString **esp, const char *progfil, int line)
{
  if (NULL == esp || NULL == *esp) {
    fprintf(stderr, "ERROR: eStringDelete on NULL pointer in %s, line %i\n",
	    progfil, line);
  } else {
    free((*esp)->strp);
    (*esp)->strp = NULL;
  }
  if (esp != NULL) {
    free(*esp);
    *esp = NULL;
  }
}

int eStringResize(EString *esp, size_t newlen, const char *progfil, int line)
{
  int errcode = ERRCODE_SUCCESS;
  size_t newsz;
  char *hp;

  if (NULL == esp) {
    fprintf(stderr, "ERROR: resize EString on NULL pointer in %s, line %i\n",
	    progfil, line);
    errcode = ERRCODE_NULLPTR;
  } else if (newlen < 1) {
    fprintf(stderr, "ERROR: resize EString to size 0 in %s, line %i\n",
	    progfil, line);
    errcode = ERRCODE_ARGRANGE;
  } else if (esp->blksz < 1) {
    fprintf(stderr, "ERROR: EString has zero blocksize %s, line %i\n",
	    progfil, line);
    errcode = ERRCODE_ASSERT;
  } else {
 
    newsz = (newlen + esp->blksz - 1)/esp->blksz;
    newsz *= esp->blksz;

    hp = erealloc(esp->strp, newsz, 0, progfil, line);
    if (NULL == hp) {
      fprintf(stderr, "ERROR: EString memory re-allocation failed in %s, line %i\n",
	      progfil, line);
      errcode = ERRCODE_NOMEM;
    } else {
      esp->strp = hp;
      esp->allocsz = newsz;
    }
  }
  
  return errcode;
}

int eStringAppend(EString *esp, const char *strp, const char *progfil, int line)
{
  int errcode = ERRCODE_SUCCESS;
  size_t slen = strlen(strp);
  if (esp->len + slen >= esp->allocsz &&
      (errcode = eStringResize(esp, esp->len + slen + 1, progfil, line)))
    return errcode;

  strcpy(esp->strp + esp->len, strp);
  esp->len = strlen(esp->strp);

  return errcode;
}
