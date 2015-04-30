/** Stream/memory magement with error messages */

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

#ifdef __cplusplus
extern "C"
{
#endif

#ifndef ELIB_H
#define ELIB_H

  //#define elib_debug

#include <stdio.h>


  /* **************************************************************************
   * ******************* ERROR MESSAGES/STATUS REPORTS ************************
   * **************************************************************************/

  enum ERRMSG_LEVELS { /**< message levels (0 and 1 reserved) */
    ERRLEVEL_FATAL = 2, /**< This is causes stop of program execution */  
    ERRLEVEL_WARNING = 3, /**< Program continues */
    ERRLEVEL_REPORT = 4,   /**< Program continues */
  };

  enum ERRMSG_CODES {
    ERRCODE_SUCCESS = 0,
    ERRCODE_FAILURE = -1,
    ERRCODE_NOMEM = 2,     /* memory allocation failed */
    ERRCODE_PREEOF = 3,    /* unexpected end of file */
    ERRCODE_SEQIOTYP = 4,  /* unknown SeqIO type */
    ERRCODE_FILTYP = 5,    /* opening wrong file type (e.g. binary when expect text) */
    ERRCODE_FASTA = 6,     /* inconsistent FASTQ/FASTA file format */
    ERRCODE_NOPROMPT = 7,  /* Next non-white char is not a FASTQ prompt */
    ERRCODE_SEQIOMODE = 8, /* wrong I/O mode, expected one of SEQIO_MODES */
    ERRCODE_NOFILE = 9,    /* File not found */
    ERRCODE_EOF = 10,      /* End of file */
    ERRCODE_WRITEERR = 11, /* Could not write to SeqIO or stream */
    ERRCODE_NOWRITE = 12,  /* Attempt to write to read file */
    ERRCODE_NOREAD = 13,   /* Attempt to read from write file */
    ERRCODE_READERR = 14,  /* Could not read from SeqIO */
    ERRCODE_NULLPTR = 15,  /* Unexpected NULL pointr */
    ERRCODE_MAXKTUP = 16,  /* Maximum k-mer length for hashing exceeded */
    ERRCODE_MAXKPOS = 17,  /* Maximum number of k-mer positions exceeded */
    ERRCODE_HASHSEQTYP = 18,/* sequence type cannot be hashed */
    ERRCODE_HASH_NOMORESEQ = 19, /* Cannot add sequences to hash table */
    ERRCODE_ARGINVAL = 20,   /* invalid argument to function */
    ERRCODE_BROKENHASH = 21, /* Hash table not correct */
    ERRCODE_SEQCODE = 22,     /* wrong, unknown or missing sequence encoding */
    ERRCODE_HASHHITSORT = 23, /* hit list not yet sorted, postitions refer to 
			      * concatenated series of sequences */
    ERRCODE_INTRANGE = 24,   /* integer outside the range of integer type */
    ERRCODE_HASHREAD = 25,    /* Read hash files are inconsistent */
    ERRCODE_SEQNOSTR = 26,   /* Expected sequences represented as character strings */
    ERRCODE_IOBUFF = 27,     /* couldn't set I/O buffer size */
    ERRCODE_COMPRESS = 28,   /* inconsistency while compressing sequence */
    ERRCODE_ARGRANGE = 29,    /* arguments outside expected range */
    ERRCODE_SHORTSEQ = 30,     /* sequence too short to be hashed */
    ERRCODE_COMMANDLINE = 31,  /* error in command line format */
    ERRCODE_ALLOCBOUNDARY = 32,/* Hit boundary of allocated memory */
    ERRCODE_FILEFORM = 33,     /* Error in file format */
    ERRCODE_SORTSTACK = 34,    /* sort stack overflow */
    ERRCODE_HITSTATUS = 35,    /* hit list does not have required status */
    ERRCODE_SEQLEN = 36,       /* sequence too long */
    ERRCODE_NOSEQ = 37,        /* No sequence to be hashed */
    ERRCODE_QUALLEN = 38,      /* qality factor sequence length does not
				* agree with base sequence length */
    ERRCODE_SEQTYP = 39,       /* wrong sequence type */
    ERRCODE_FREQSUMS = 40,     /* Rounding error when calculating frequency sums */
    ERRCODE_SWATEXCEED = 41,   /* Overflow of Smith-Waterman score */
    ERRCODE_ALISCOREXCEED = 42,/* Alignment score exceeds target while back tracking */
    ERRCODE_DIFFCOUNT = 43,    /* Mutation count inconsistent in diff string */
    ERRCODE_SWATSCOR = 44,     /* Inconsistency when calculating Smith-Waterman score */
    ERRCODE_HTLSTQLEN = 45,    /* inconsistent query length between hit lists */
    ERRCODE_MATCHSTAT = 46,    /* Match list does not have correct status */
    ERRCODE_ASSERT = 47,       /* assertion failed */
    ERRCODE_OVERFLOW = 48,     /* integer overflow */
    ERRCODE_FILPOS = 49,       /* setting file position */
    ERRCODE_FILEIO = 50,       /* file input/output error */
    ERRCODE_MENUE = 51,        /* error when parsing command line */
    ERRCODE_ENDIAN = 52,       /* file written under different endianess */
    ERRCODE_FASTAPAIRNUM = 53, /* The number of reads in the two FASTA/FASTQ input files differ */
    ERRCODE_LONGFILNAM = 54,   /* file name too long */
    ERRCODE_FHEADSIZ = 55,     /* wrong size of file type specific header in binary file */
    ERRCODE_FILVERSION = 56,   /* wrong format version of binary file */
    ERRCODE_LONGSEQ = 57,      /* sequence too long */
    ERRCODE_NOSEQIDX = 58,     /* sequence index missing */
    ERRCODE_DIFFSTR = 59,      /* Inconsistent diff-string */
    ERRCODE_SEQOFFS = 60,      /* Sequence offset greater than sequence length */
    ERRCODE_TERMIT = 61,       /* Terminate iterative calls to a function */
    ERRCODE_CPLXSCOR = 62,     /* complexity weighted score exceeds unweighted score */
    ERRCODE_PTHREAD = 63,      /* Error occurred in one of the threads */
    ERRCODE_HITINFO = 64,      /* HashHitInfo has wrong status */
    ERRCODE_HITBIN = 65,       /* K-mer statistics insufficient for binning */
    ERRCODE_SWATSTRIP = 66,    /* Inconsistency in striping format for Smith-Waterman using SSE2 */
    ERRCODE_QUALVAL = 67,      /* Invalid FASTQ base quality value */
    ERRCODE_SEGPARTFETCH = 68, /* Fetched only part of requested segment from sequence set */
    ERRCODE_NOHASHWORD = 69,   /* Word not found in hash index */
    ERRCODE_PTHREADSTACK = 70, /* Overflow in thread stack */
    ERRCODE_FILSIZ = 71,       /* File size too big */
    ERRCODE_SKIPSMALL = 72,    /* sampling step size to small for total length of reference sequences */
    ERRCODE_SEQSETSIZ = 73,    /* Total number of bases in set of sequences exceeds limit */
    ERRCODE_PAIRNUM = 74,      /* Number of pairs exceeds limit */
    ERRCODE_NOMATCH = 75,      /* Found no match */
    ERRCODE_RNAMPAIR = 76,     /* Read names don't match as expected for a read pair */
    ERRCODE_MATENUMPAIR = 77,  /* Number of reads differs in the two input files for paired reads */
    ERRCODE_BREAK = 78,        /* Early termination triggered by debugging code */
    ERRCODE_PTHRTERMSIG = 79,  /* Termination signal for thread */
    ERRCODE_PTHRPULLCHK = 80,  /* could not pull thread argument from buffer because
				* check function was not fulfilled */
    ERRCODE_SEQNAMLEN = 81,    /* Sequence name too long */
    ERRCODE_KPOSOFLO = 82,     /* Number of positions in hash index caused
				* integer overflow */
    ERRCODE_INFMT = 83,        /* Unrecognised input format */
    ERRCODE_OUFMT = 84,        /* Unrecognised outpu format */
    ERRCODE_SEQNTSMBL = 85,    /* Nucleotide symbol is not a letter */
    ERRCODE_SEQNUMSET = 86,    /* Number of sequences exceeds limit */
    ERRCODE_BAMBAM = 87,       /* bambamc library returns error */
    ERRCODE_SEMOPEN = 88,      /* named semaphore could not be opened */
  };

  typedef struct _ErrMsg ErrMsg;

  /** Return error description for code.
   * Use this to ouput return status */
  const char *errMsgString(int errcode);

  /**< Constructor */
  ErrMsg *errMsgCreate(const char *progfil, int linenum);
  
  /**< Destructor */
  void errMsgEnd(ErrMsg *p);

  /**< Set the current read name (in order to trace where error occured */
  void errMsgSetCurrentReadName(ErrMsg *emp, const char *namp);

  /** Set the current read number */
  void errMsgSetCurrentReadNumber(ErrMsg *emp, size_t rno);

  /** Add an error message to the stack */
  void errMsgAdd(ErrMsg *emp, const char *message, 
		 const char *progfil, int linenum,
		 int errcode, unsigned char level);

  /* print Message list on stream, delete entries */
  void errMsgFlush(FILE *fp, ErrMsg *emp);

#define ERRMSG_CREATE(emp) ((emp) = errMsgCreate(__FILE__, __LINE__))
#define ERRMSG_END(emp) errMsgEnd(emp); emp = NULL;
#define ERRMSG_FLUSH(emp) errMsgFlush(stdout, emp)
#define ERRMSG_READNAM(emp, readnam) errMsgSetCurrentReadName((emp), (readnam))
#define ERRMSG_READNO(emp, readno) errMsgSetCurrentReadNumber((emp), (readno))
#define ERRMSG(emp, msg, code) errMsgAdd((emp),(msg),__FILE__, __LINE__,(code),ERRLEVEL_FATAL)
#define ERRMSGNO(emp, code) errMsgAdd((emp),"",__FILE__, __LINE__,(code),ERRLEVEL_FATAL)
#define ERRWARNMSG(emp, msg, code) errMsgAdd((emp),(msg),__FILE__, __LINE__,(code),ERRLEVEL_WARNING)
#define ERRNOMSG(emp, code) errMsgAdd((emp),"",__FILE__, __LINE__,(code),ERRLEVEL_WARNING)



  /* **************************************************************************
   * ******* MEMORY ALLOCATION, FILE OPENING WITH ON-SCREEN MESSAGES **********
   * **************************************************************************/

  void *emalloc(size_t, const char *, int);
  void *ecalloc(size_t, size_t, const char *, int);
  void *erealloc(void *, size_t, size_t, const char *, int);
  void *ereallocp(void **, size_t, const char *, int);
  FILE *efopen(const char *, const char *, const char *, int);
  /* fopen with error message */
  int efclose(FILE *, const char *, int line);
  char *estrcpy(const char *, const char *, int);
  char *estrcat(const char *, const char *, const char *, int);
  
#define EMALLOC(size) emalloc((size), __FILE__, __LINE__)
#define EMALLOC0(size) ecalloc(1, (size), __FILE__, __LINE__)
#define EMALLOCP(p) ((p) = emalloc(sizeof (*(p)), __FILE__, __LINE__))
#define EMALLOCP0(p) ((p) = ecalloc(1, sizeof(*(p)), __FILE__, __LINE__))
#define ECALLOC(nobj, size) ecalloc((nobj), (size), __FILE__, __LINE__)
#define ECALLOCP(nobj, p) ((p) = ecalloc((nobj),sizeof(*(p)), __FILE__, __LINE__))
#define EREALLOC(objp, size) erealloc((objp), (size), 0, __FILE__, __LINE__)
#define EREALLOCP(objp, nobj) erealloc((objp), (nobj)*sizeof(*(objp)), 0, __FILE__, __LINE__)
#define EREALLOCPP(objp, nobj) ereallocp((void **) &(objp), (nobj)*sizeof(*(objp)), __FILE__, __LINE__)

  /** reallocate and initialise newly allocated memory */
#define EREALLOCP0(objp, nobj, nobj_old) erealloc((objp), (nobj)*sizeof(*(objp)),\
                                                          (nobj_old)*sizeof(*(objp)), __FILE__, __LINE__)
#define EFOPEN(filnam, mode) \
        efopen((filnam), (mode), __FILE__, __LINE__)
#define EFCLOSE(fp) efclose((fp), __FILE__, __LINE__)
#define ESTRCPY(tostr, fromstr) ((tostr) = estrcpy((fromstr), __FILE__, __LINE__))
#define ESTRCAT(tostr, fromstrA, fromstrB) ((tostr) = estrcat((fromstrA), (fromstrB), __FILE__, __LINE__))


  /****************************************************************************
   *************************** Type EString ***********************************
   ****************************************************************************/
typedef struct EString_ { /** Hold a dynamically allocated string */
  char *strp;
  size_t len;
  size_t allocsz;
  int blksz;
} EString;

#define ESTRING_NEW(p) ((p) = eStringCreate(0, __FILE__, __LINE__))
#define ESTRING_INIT(s) eStringInit(&(s), 0, __FILE__, __LINE__)
#define ESTRING_DELETE(p) eStringDelete(&(p), __FILE__, __LINE__)
#define ESTRING_FREE(s) free((s).strp); memset(&(s),0,sizeof(EString));
#define ESTRING_BLANK(s) (s).strp[0] = '\0';(s).len = 0
#define ESTRING_RESIZE(s,len) eStringResize(&(s), (len), __FILE__, __LINE__) 
#define ESTRING_LENGTH(s) (s).len
#define ESTRING_APPEND(s,cp) eStringAppend(&(s), cp, __FILE__, __LINE__)
#define ESTRING_GETSTR(s) (s).strp

  char *eStringInit(EString *esp, int blksz, const char *progfil, int line);
  EString *eStringCreate(int blksz, const char *progfil, int line);
  void eStringDelete(EString **esp, const char *progfil, int line);
  int eStringResize(EString *esp, size_t newlen, const char *progfil, int line);
  int eStringAppend(EString *esp, const char *strp, const char *progfil, int line);

#endif

#ifdef __cplusplus
}
#endif
