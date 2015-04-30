/** Handling of nucleotide sequences (in FASTQ/FASTA format) */

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

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#ifdef HAVE_ZLIB_H
#include <zlib.h>
#endif
#include "filio.h"
#include "sequence.h"

enum SEQ_CONSTANTS {
  FALSE            = 0,
  TRUE             = 1,
  MODELEN          = 3,
  SIZE_ALPHABET    = 7,          /**< number of letters in 3-bit alphabet (CODEC_ALPHABET) - 1
				    do not count termination signal */
  NBITS_ALPHABET   = 3,          /**< Number of bits for 3-bit alphabet */
  ALPHABET_MASK    = 0x07,       /**< For masking out alphabet index */
  SIZE_STANDARD_ALPHABET = 4,    /**< Number of letters in 2-bit alphabet */
  UNKNOWN_3BIT     = 5,          /**< 3-bit code of the unknown character 'N' */
  NBYTES_PER_INT32 = 4,
  MASK_32BIT = 0xFFFFFFFF,
  NUM_32BIT = 32,
  NBITS_PER_BYTE   = 8,          
  SIZE_CODTAB      = 0x100,      /**< number of 8-bit ASCII characters */
  SIZE_DECODTAB    = 0x100,
  BLOCKSIZE_HEADER = 128,        /* default lenght for FASTA string after header */
  OUTBUFSIZE       = 80,         /* buffer size for output (line width) */
  MAXLEN_HEADER    = 2048,       /* maximum header length */
  MAXN_PER_UNIT    = 10,         /* number of bases per 32-bit unit (3-bit code) */
  MAXN_TERMCHAR    = 8,          /* maximum number of termination characters in a sequence */
  SEQ_BLKSZ_DEFAULT= 256,        /* minimum block size for memory allocation of sequences */
  LINEWIDTH_DEFAULT= 60,
  LINBUFSIZ = 1024,              /**< default size of line buffer for reading */
  SEQ_MAXLEN = UINT32_MAX
};

enum SEQCOD_CONST {
  FASTA_PROMPT = '>',
  FASTQ_PROMPT_SEQ = '@',
  FASTQ_PROMPT_QUAL = '+',
  CODEC_UNKNOWN_CHAR = 'N',
  CODEC_OFFSET_CHAR = 'A', /**< start of offset */
  MASKBIT4 = 0x08,
};

enum SEQSET_CONST {
  SEQSET_BLOCKSIZ_DEFAULT = 4096,    /**< default blocksize for memory allocation in SeqSet */
  SEQSET_NAMBLOCKSIZ_DEFAULT = 1024, /**< default blocksize for memory allocation in SeqSet */
  SEQSET_FORMAT_VERSION = 4,         /**< file format version */
  SEQSET_HEADLEN = 8,                /**< number of UINT32 fields of the file header */
};

enum READNAM_TYPES {
  READNAM_TYP_UNKNOWN = 0,
  READNAM_TYP_ILLUMINA = 1,
  READNAM_TYP_RF = 2,
  READNAM_NUMTYP = 3,
};

enum SEQIO_FLAGS {
  SEQIOFLG_WRITE= 0x01,
  SEQIOFLG_GZIP = 0x02,  /**< gzip compressed */ 
};

static const char READNAM_MATEXT_SEPARATOR[READNAM_NUMTYP] = {'\0', '/', '.'}; /* index corresponds to READNAM_TYPES */
static const char READNAM_MATEXT_ILLUMINA[2][2] = {"1", "2"};
static const char READNAM_MATEXT_FR[2][2] = {"F", "R"};

static const char SEQSET_FILNAMEXT[] = "sma";
static const char CODEC_ALPHABET[] = "ACGTXN";

typedef unsigned char BOOL_t;
typedef unsigned char UCHAR_t;

typedef struct _SEQSEQ {  /**< Sequence container (string, nucleotides, quality factors) */
  char code;              /**< one of SEQ_CODES */
  char *basep;            /**< allocated memory for sequence data (terminated with 0). In
			   * SeqSet.sqp: contains several strings each terminated by 0 */
  int block_size;         /**< blocksize for memory allocation */
  SETSIZ_t size;       /**< sequence length (number of symbols, excl term char). In 
			   * SeqSet.sqp: includes terminating '\\0' except final '\\0' */
  size_t alloc_size;  /**< allocated memory in bytes */
  char nbit_symb;         /**< number of bits for one symbol */
} SEQSEQ;
  
/******************************************************************************
 ******************************* Opaque Types *********************************
 ******************************************************************************/
struct _SeqCodec { /**< Code table */
  char typ;                        /** codec type (one of SEQCODEC_CODETYPES) */
  UCHAR_t alphlen;                   /**< size of the alphabet must be <= SIZE_ALPHABET */
  char alphabet[SIZE_ALPHABET];    /**< Standard nucleotide codes */
  UCHAR_t codtab[SIZE_CODTAB];       /**< Encoding table */
  char decodtab[SIZE_DECODTAB];    /**< Decoding table */
  unsigned char codtab_complement[SIZE_STANDARD_ALPHABET]; /* (encoded) complement table */
};

struct _SeqFastq { /**< Holds nucleotide sequence with base call quality values */
  char type;      /**< one of SEQ_TYPES */
  SEQSEQ *headp;  /**< Header of the nucleotide sequence */
  SEQSEQ *datap;  /**< Nucleotide sequence */
  SEQSEQ *qheadp; /**< Header of quality values (can be NULL if type == SEQTYP_FASTA) */ 
  SEQSEQ *qualp;  /**< Sequence of quality values (can be NULL if type == SEQTYP_FASTA) */ 
};

struct _SeqIO { /**< Wrapper of I/O stream from which data is read */
  UCHAR_t flags;        /**< combination of SEQIO_FLAGS */
  char mode;            /**< One of SEQIO_MODES */
  char *filnam;         /**< file name */
  void *fp;            /**< I/O stream can be FILE or gzFile */
  char fmode[MODELEN+1];/**< fopen() I/O mode */
  size_t bufsize;       /**< I/O buffer size */
  int status;           /**< == 0 on success, otherwise one of ERRMSG_CODES */
  char *linbufp;        /**< line buffer for reading */
};

struct _SeqSet { /**< Holds a set of sequences compressed in memory */
  UCHAR_t statusflag;   /**< combination of SEQSET_FLAGS */
  SEQSEQ *sqp;       /**< Sequence of concatenated sequences */
  SEQSEQ *qqp;       /**< Sequence of concatenated quality values (may be NULL) */
  SEQNUM_t n_seq;    /**< number of sequences */
  SETSIZ_t  *sop; /**< Array of n_seqp+1 offsets: sop[i] is the position of the first base
		      * of sequence i if all seqeuences were concatenated). sop[nseq] is the
		      * total number of bases. */
  SEQNUM_t n_alloc;  /* number of sequences for which space is allocated */
  int blocksiz;      /* block size for memory allocation */
  char *namebasep;   /**< Base address of sequence names */
  SETSIZ_t *namoffs;/**< Array of n_seqp+1 offsets to sequence names allocated at namebasep
			* namoffs[n_seqp+1] is the length of all sequence names in namebasep
			* concatenated. */
  size_t nam_alloc;  /**< size of memory allocated for sequence names */
  int nam_blocksiz;  /**< block size for memory allocation of sequence names */

  SEQNUM_t *sxp;     /**< Array of size n_seq+1 (can be NULL) used to
		      * set sequence indices into arrays of sorted offsets
		      * using sop. This allows assigning sequence
		      * indices to offsets in the sequence of
		      * concatenated sequences and going to offsets per
		      * sequence. */
};

/******************************************************************************
 ********************************** Macros ************************************
 ******************************************************************************/

#define CALC_UNIT_OFFSET(siz, offset, nsu)\
   (offset) = (siz)/MAXN_PER_UNIT;\
   (nsu) = ((offset) + 1)*MAXN_PER_UNIT - (siz);

/* Macro used in methods of type SEQSEQ for dynamic allocation of memory */
#define SEQ_REALLOC(seqptr, baseptr, basectr) {\
 if (((basectr) + 1 >= (seqptr)->alloc_size) &&  \
   reallocSeqBlocks(seqptr, basectr + 2)) return ERRCODE_NOMEM; \
 (baseptr) = (seqptr)->basep + (basectr);\
 }

/******************************************************************************
 ****************************** Private Methods *******************************
 ******************************************************************************/

static int scrollToHeaderLine(
#ifdef HAVE_ZLIB
			      gzFile fp,
#else
			      FILE *fp, 
#endif
			      int *prompt, char bufp[LINBUFSIZ])
{
  char *cp = bufp;
  *prompt = '\0';
  
  do {
    for(;isspace((int) *cp); cp++);
    if (*cp == FASTA_PROMPT ||
	*cp == FASTQ_PROMPT_SEQ) {
      *prompt = *cp;
      return ERRCODE_SUCCESS;
    }
  } while ((cp =
#ifdef HAVE_ZLIB
	    gzgets(fp, bufp, LINBUFSIZ)));
  return (gzeof(fp))?
#else
    fgets(bufp, LINBUFSIZ, fp)));

  return (feof(fp))?
#endif
  ERRCODE_EOF: ERRCODE_READERR;
}

static int cmpPairNamStr(const char *ap, const char *bp, size_t maxlen)
{
  int i, rv = 0;
  BOOL_t isMate1, isMate2;
  const char *extstr[2];
  size_t s;
  for (s=0; s<maxlen && ap[s] == bp[s] && ap[s] != '\0'; s++);

  if (ap[s] != bp[s]) {
    rv = (ap[s] > bp[s])? 1 : -1;
    if (s > 0) {
      for (i=1; i<READNAM_NUMTYP && rv != 0; i++) {
	if (ap[s-1] == READNAM_MATEXT_SEPARATOR[i]) {
	  extstr[0] = NULL;
	  switch (i) {
	  case READNAM_TYP_UNKNOWN:
	    break;
	  case READNAM_TYP_ILLUMINA:
	    extstr[0] = READNAM_MATEXT_ILLUMINA[0];
	    extstr[1] = READNAM_MATEXT_ILLUMINA[1];	
	    break;
	  case READNAM_TYP_RF:
	    extstr[0] = READNAM_MATEXT_FR[0];
	    extstr[1] = READNAM_MATEXT_FR[1];
	  default:
	    break;
	  }
	  if (extstr[0] != NULL) {
	    for (ap += s; *ap == *extstr[0] && *extstr[0] != '\0'; ap++, extstr[0]++);
	    isMate1 = *extstr[0] == '\0' && (*ap == '\0' || isspace(*ap));
	    for (bp += s; *bp == *extstr[1] && *extstr[1] != '\0'; bp++, extstr[1]++);
	    isMate2 = *extstr[1] == '\0' && (*bp == '\0' || isspace(*bp));
	    if ((isMate1) && (isMate2))
	      rv = 0;
	  }
	}
      }
    }
  }
  return rv;
}
/******************************************************************************
 ************************ Pivate Methods of Type SeqCodec *********************
 ******************************************************************************/
static int checkCodec(SeqCodec *codep) 
{
  int errcode = ERRCODE_SUCCESS;
  const char bases[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
  const char *bp;

  for (bp = bases; (*bp); bp++) {
    if (toupper(*bp) != (int) codep->decodtab[codep->codtab[(unsigned char) *bp]] &&
	(toupper(*bp) != 'U' || 
	 (int) codep->decodtab[codep->codtab[(unsigned char) *bp]] != 'T'))  {
      errcode = ERRCODE_FAILURE;
      break;
    }
  }
  if (errcode) {
    printf("ERROR in code table for base: %c, code: %d, decode: %c\n",
	   *bp, codep->codtab[(signed char) *bp], codep->decodtab[codep->codtab[(unsigned char) *bp]]);
  }
  return errcode;
}

static int make3BitMangledCodec(SeqCodec *codp)
{
  int i, offs;
  unsigned char cu, a;
  
  codp->typ = SEQCODTYP_3BITMANGLED;
  codp->alphlen = (UCHAR_t) strlen(CODEC_ALPHABET);
  if (codp->alphlen > SIZE_ALPHABET) 
    return ERRCODE_ASSERT;

  strcpy(codp->alphabet, CODEC_ALPHABET);
  for (i=1; i<SIZE_CODTAB; i++) {
    cu = (unsigned char) toupper(i);
    if (cu == 'U') cu = 'T';
    offs = ((int) cu) - (int) CODEC_OFFSET_CHAR + 1; /* offset from ASCII code for 'A' */
    if (offs > 0 && offs <32) {
      for (a=0; (a<SIZE_STANDARD_ALPHABET) && (cu != codp->alphabet[a]); a++);
      if (a>=SIZE_STANDARD_ALPHABET) a = UNKNOWN_3BIT;
      else codp->codtab_complement[(~a)&SEQCOD_STDNT_MASK] = a + (((unsigned char) offs)<<3);
      codp->codtab[i] = a + (((unsigned char) offs)<<3);
      codp->decodtab[codp->codtab[i]] = cu;
    } else {
      a = UNKNOWN_3BIT; /* unknown character convert to 'N' */
      offs = (int) codp->alphabet[a] - codp->alphabet[0] + 1;
      codp->codtab[i] = a + (((unsigned char) offs)<<3);
    }
  }
  codp->codtab[0] = SEQCOD_TERM;
  codp->decodtab[SEQCOD_TERM] = '\0';
  
  return ERRCODE_SUCCESS;
}

/******************************************************************************
 ************************ Public Methods of Type SeqCodec *********************
 ******************************************************************************/

SeqCodec *seqCodecCreate(void)
     /** Create a table for an 8-bit code where bits 0 and 1 represent
      * the 2-bit code for the standard nucleotides, bits 0-2
      * represent a 3-bit code of the recognised alphabet including
      * 0x07 as sequence termination signal.  Bits 3-7 represent the
      * original letter as the offset of the ASCII code of the
      * corresponding capital letter from the ASCII code 'A' plus 1.
      * => the code is non-zero and char string can be terminated in the usual way.
      * Bit 2 is set if non-standard nucleotide or termination signal
      * The 2-bit code has the property that the complement code is obtained
      * via bitwise complement. All non-standard NT codes are cast to 0x01 ('C').
      */
{
  int errcode;
  SeqCodec *codp;

  EMALLOCP0(codp);
  if (!codp) return NULL;

  errcode = make3BitMangledCodec(codp);
  if (errcode ||
      checkCodec(codp)) {
    seqCodecDelete(codp);
    codp = NULL;
  }

  return codp;
}

SeqCodec *seqCodecCreateFromAlphabet(const char *alphabet, unsigned char code_unknown)
{
  int i, al;
  UCHAR_t cu, a;
  size_t alphsiz = strlen(alphabet);
  SeqCodec *codp;

  if (alphsiz > SIZE_ALPHABET || alphsiz < SIZE_STANDARD_ALPHABET ||
      code_unknown >= alphsiz) return NULL;

  al = (int) alphsiz;

  EMALLOCP0(codp);
  if (!codp) return NULL;
  
  codp->typ = SEQCODTYP_USERDEFINED;
  strcpy(codp->alphabet, alphabet);
  codp->alphlen = (UCHAR_t) al;

  for (i=0; i<SIZE_DECODTAB; i++)
    codp->decodtab[i] = alphabet[code_unknown];

  for (i=0; i<al; i++)
    codp->alphabet[i] = (char) toupper(codp->alphabet[i]);

  for (i=0; i<SIZE_CODTAB; i++) {
    cu = (unsigned char) toupper(i);
    if (cu == 'U') cu = 'T';
    for (a=0; (a<al) && (cu != (UCHAR_t) codp->alphabet[a]); a++);
    if (a < al) {
      if (a < SIZE_STANDARD_ALPHABET) {
	codp->codtab_complement[(~a)&SEQCOD_STDNT_MASK] = a | MASKBIT4;
      }
      a |= MASKBIT4; /* set the 4th bit (so that the (8-bit) code is != '\0' */
      codp->decodtab[a] = cu; 
    } else {
      a = code_unknown | MASKBIT4; /* set the 4th bit (so that the (8-bit) code is != '\0' */
    }
    codp->codtab[i] = a;
  }

  
  return codp;
}

void seqCodecDelete(SeqCodec *codp)    
{
  free(codp);
}  

const char *seqCodecGetAlphabet(const SeqCodec *codp, short *length)
{
  if (length)
    *length = (short) strlen(codp->alphabet);
  return codp->alphabet;
}

const unsigned char *seqCodecGetEncoder(const SeqCodec *codp, short *size)
{
  if (size)
    *size = SIZE_CODTAB;
  return codp->codtab;
}

const char *seqCodecGetDecoder(const SeqCodec *codp, short *size)
{
  if (size)
    *size = SIZE_DECODTAB;
  return codp->decodtab;
}

void seqCodecEncodeString(char *cp, const SeqCodec *codp)
{
  const UCHAR_t *encodp = codp->codtab;
  for (;(*cp); cp++) 
    *cp = (char) encodp[(int) *cp];
}

void seqCodecDecodeString(char *cp, const SeqCodec *codp)
{
  char c;
  const char *decodp = codp->decodtab;
  for (;(*cp); cp++) {
    c = decodp[(UCHAR_t) *cp];
    *cp = c;
  }
}

char seqCodecFindBaseClass(char c, const SeqCodec *codecp)
{
  UCHAR_t cod = codecp->codtab[(int) c];
  char rv;
  if (cod&0x04) {
    if ((cod&ALPHABET_MASK) == UNKNOWN_3BIT) rv = SEQBASTYP_UNKNOWN;
    else rv = SEQBASTYP_NONSTD;
  } 
  else if (cod&0x01) rv = SEQBASTYP_PYRIMIDINE;
  else rv = SEQBASTYP_PURINE;
  return rv;
}

char seqCodecType(const SeqCodec *codecp)
{
  return codecp->typ;
}

/******************************************************************************
 ************************ Public Methods of Type SeqIO ************************
 ******************************************************************************/
SeqIO *seqIOopen(int *errcode, const char *filnam, char mode, unsigned long buffsize)
     /**< Open a file for sequence I/O in FASTA/FASTQ format (7-bit ascii)
      * \param errcode != 0 (one of ERRMSG_CODES) 
      *        if error occurs in which case NULL is returned.
      * \param filnam file name. 
      * \param mode one of SEQIO_MODES 
      *             (SEQIO_READ, SEQIO_WRITE_FASTA, SEQIO_WRITE_FASTQ) 
      * \param buffsize buffer size, if 0 default buffer size is used.
      */
{
  SeqIO *p = 0;
  *errcode = ERRCODE_SUCCESS;
  if (!(filnam) || filnam[0] == '\0' || isspace((int) filnam[0])) {
    *errcode = ERRCODE_ASSERT;
    return NULL;
  }
  EMALLOCP0(p);
  if (!p) {
    *errcode = ERRCODE_NOMEM;
    return NULL;
  }
  p->flags = 0;
  p->bufsize = (buffsize)? buffsize: BUFSIZ;
  p->mode = mode;
  p->linbufp = NULL;
  switch(mode) {
  case SEQIO_READ:
    strcpy(p->fmode, "r");
#ifdef HAVE_ZLIB
  case SEQIO_WRITE_GZIP_FASTA:
  case SEQIO_WRITE_GZIP_FASTQ:
    p->flags = SEQIOFLG_GZIP;
#endif
    break;
  case SEQIO_WRITE_FASTA:
  case SEQIO_WRITE_FASTQ:
    p->flags |= SEQIOFLG_WRITE;
    strcpy(p->fmode, "w");
    break;
  default:
    *errcode = ERRCODE_SEQIOMODE;
    break;
  }

  if (!*errcode) {
    p->filnam = NULL;
    ESTRCPY(p->filnam, filnam);
    if (!p->filnam) *errcode = ERRCODE_NOMEM;
  }
  
  if (!*errcode) {
#ifdef HAVE_ZLIB
    if (mode == SEQIO_READ || 
	mode == SEQIO_WRITE_GZIP_FASTA || mode == SEQIO_WRITE_GZIP_FASTA) {
	p->fp = gzopen(p->filnam, p->fmode);
#else
    if (p->filnam[0] == '-' && 
	(p->filnam[1] == '\0' || isspace(filnam[1]))) {
      p->fp = (mode == SEQIO_READ)? stdin: stdout;
#endif
    } else {
      p->fp = EFOPEN(p->filnam, p->fmode);
    }

    if (!p->fp) *errcode = ERRCODE_NOFILE;
  }

#ifndef HAVE_ZLIB
  if (!(*errcode) &&
      setvbuf(p->fp, NULL, _IOFBF, p->bufsize))
    *errcode = ERRCODE_IOBUFF;
#endif
  
  if (!(*errcode) && mode == SEQIO_READ) {
    ECALLOCP(LINBUFSIZ+1, p->linbufp);
    if (p->linbufp) {
      /* read the first line */
      if (
#ifdef HAVE_ZLIB
	  !gzgets(p->fp, p->linbufp, LINBUFSIZ)
#else
	   !fgets(p->linbufp, LINBUFSIZ, p->fp)
#endif
	  )
	*errcode = ERRCODE_FASTA;
    } else {
      *errcode = ERRCODE_NOMEM;
    }
  }

  if (*errcode) {
    seqIOclose(p);
    p = NULL;
  }
  
  return p;      
}

int seqIOclose(SeqIO *p)
{
  int errcode = ERRCODE_SUCCESS;
  if (p) {
    if (p->fp) {
#ifdef HAVE_ZLIB  
      if ((p->flags & SEQIOFLG_GZIP) || !(p->flags & SEQIOFLG_WRITE)) {
	if (Z_OK != gzclose(p->fp))
	  errcode = ERRCODE_FILEIO;
      } else {
#endif
      if (fflush(p->fp))
	errcode = ERRCODE_FILEIO;
      if (p->filnam[0] != '-' || 
	  (p->filnam[1] != '\0' && !isspace(p->filnam[1])))
	fclose(p->fp);
#ifdef HAVE_ZLIB
      }
#endif
    }

    if (p->linbufp) free(p->linbufp);
    free(p->filnam);
  }
  free(p);

  return errcode;
}

int seqIOReset(SeqIO *p)
{
  if (p->mode == SEQIO_READ && 
      p->status != ERRCODE_SUCCESS && p->status != ERRCODE_EOF)
    return p->status;
  p->status = ERRCODE_SUCCESS;

#ifdef HAVE_ZLIB
  if ((p->flags & (SEQIOFLG_GZIP &~SEQIOFLG_WRITE))) {
    gzrewind(p->fp);
    if (p->mode == SEQIO_READ &&
	!gzgets(p->fp, p->linbufp, LINBUFSIZ))
      p->status = (gzeof(p->fp))? ERRCODE_EOF: ERRCODE_READERR;
  } else {
#endif
    rewind(p->fp);

    /* read the first line */
    if (p->mode == SEQIO_READ &&
	!fgets(p->linbufp, LINBUFSIZ, p->fp))
      p->status = (feof(p->fp))? ERRCODE_EOF: ERRCODE_READERR;
#ifdef HAVE_ZLIB
  }
#endif

  return p->status;
}

int seqIOstatus(SeqIO *p)
{
  if (!p) return ERRCODE_NULLPTR;
  return p->status;
}

 int seqIOCheckReads(ErrMsg *errmsgp, SeqFastq *sqbufp, SeqIO *sfp,
		     SeqFastq *sqbufBp, SeqIO *sfBp,
		     SEQNUM_t *seqnum, SEQLEN_t *maxseqlen, 
		     SEQLEN_t *maxnamlen)
{
  int errcode;
  size_t snum;
  BOOL_t isPaired = sqbufBp != NULL && sfBp != NULL;
  BOOL_t namflg = 1;
  SEQLEN_t mxnaml, mxseql;

  if (sfp->mode != SEQIO_READ) 
    return ERRCODE_FAILURE;
  mxnaml = mxseql = 0;
  for (snum=0;
       !sfp->status && (!isPaired || !sfp->status) && snum < UINT_MAX;
       snum++) {
    errcode = seqFastqRead(sqbufp, sfp);
    ERRMSG_READNO(errmsgp, snum);
    ERRMSG_READNAM(errmsgp, seqFastqGetSeqName(sqbufp));

    if ((errcode)) {
      if (errcode != ERRCODE_EOF)
	ERRMSGNO(errmsgp, errcode);
      return errcode;
    }
 
    if (sqbufp->headp->size > mxnaml) mxnaml = sqbufp->headp->size;
    if (sqbufp->datap->size > mxseql) mxseql = sqbufp->datap->size;
    if ((sqbufp->qheadp) && sqbufp->qheadp->size > mxnaml)
      mxnaml = sqbufp->qheadp->size;

    if (isPaired) {
      if ((errcode = seqFastqRead(sqbufBp, sfBp)))
      return errcode;
      if ( (namflg) ) {
	size_t maxlen = (sqbufp->headp->size > sqbufBp->headp->size)?
	  sqbufBp->headp->size: sqbufp->headp->size;
	if (cmpPairNamStr(sqbufp->headp->basep, sqbufBp->headp->basep, maxlen))
	  namflg = 0;
      }
      if (sqbufBp->headp->size > mxnaml) mxnaml = sqbufBp->headp->size;
      if (sqbufBp->datap->size > mxseql) mxseql = sqbufBp->datap->size;
      if ((sqbufBp->qheadp) && sqbufBp->qheadp->size > mxnaml)
	mxnaml = sqbufBp->qheadp->size;
    }
  }

  if (snum >= UINT_MAX) return ERRCODE_OVERFLOW;
  if (seqnum) *seqnum = (SEQNUM_t) snum;
  if (maxseqlen) *maxseqlen = (unsigned int) mxseql;
  if (maxnamlen) *maxnamlen = (unsigned int) mxnaml;


  if (sfBp == NULL) {
    /* single reads */
    if (sfp->status == ERRCODE_EOF) /* should return EOF */
      errcode = ERRCODE_SUCCESS;
    else if (sfp->status == ERRCODE_SUCCESS) /* did not return EOF (suspicious) */
      errcode = ERRCODE_FASTA;
    else 
      errcode = sfp->status;
    return errcode;
  }

  /* paired reads */
  if (sfp->status == ERRCODE_EOF) {
    if (sfBp->status == ERRCODE_EOF)
      errcode = ERRCODE_SUCCESS;
    else if (sfBp->status == ERRCODE_SUCCESS)
      errcode = ERRCODE_MATENUMPAIR;
    else 
      errcode = sfBp->status;
  } else if (sfp->status == ERRCODE_SUCCESS)
    errcode =  ERRCODE_MATENUMPAIR;
  else
    errcode = sfp->status;
    
  return ( (namflg) )? errcode: ERRCODE_RNAMPAIR;
}

/******************************************************************************
 *********************** Private Methods of Type SEQSEQ ***********************
 ******************************************************************************/

static void deleteSeq(SEQSEQ *sp)
     /**< Destructor */
{
  if (sp != NULL)
    free(sp->basep);
  free(sp);
}

static SEQSEQ *createSeq(uint32_t blocksize)
     /**< Constuctor.
      * \param blocksize Size of the blocks for dynamic memory allocation
      *                  in bytes.
      */
{
  SEQSEQ *sp;

  EMALLOCP0(sp);
  if (sp == NULL) 
    return NULL;

  if (blocksize < 1) blocksize = SEQ_BLKSZ_DEFAULT;
  ECALLOCP(blocksize, sp->basep);
  
  if (sp->basep == NULL) {
    deleteSeq(sp);
    return NULL;
  }
  sp->code = SEQCOD_ASCII;
  *sp->basep = '\0';
  sp->alloc_size = sp->block_size = blocksize;
  sp->size = 0;
  sp->nbit_symb = NBITS_PER_BYTE;

  return sp;
}

static int reallocSeqBlocks(SEQSEQ *sp, size_t minsize)
{
  size_t newsiz = minsize/sp->block_size + 1;
  newsiz *= sp->block_size;
  char *hp = EREALLOCP(sp->basep, newsiz);
  if (hp == NULL) 
    return ERRCODE_NOMEM;
  sp->basep = hp;
  sp->alloc_size = newsiz;

  return ERRCODE_SUCCESS;
}

static void blankSeq(SEQSEQ *sp)
     /**< Reset parameters such as sequence length */
{
  if (sp) {
    if (sp->basep) sp->basep[0] = '\0';
    sp->size = 0;
    sp->code = SEQCOD_ASCII;
    sp->nbit_symb = NBITS_PER_BYTE;
  }
}

static int finaliseSeq(SEQSEQ *sp)
     /**< Free unused memory */
{
  char *hp;
  if (!sp->basep) return ERRCODE_SUCCESS;
  hp = EREALLOCP(sp->basep, sp->size + 1);
  if (!hp) return ERRCODE_NOMEM;
  sp->basep = hp;
  sp->alloc_size = sp->size + 1;
  return ERRCODE_SUCCESS;
}

static int setSeq(SEQSEQ *sp, const char *cp)
     /**< Copy sequence from a string. Memory is dynamically allocated.
      * skip leading and trailing white space */
{
  char *hp;
  SEQLEN_t i;

  if (cp == NULL) return ERRCODE_SUCCESS;
  hp = sp->basep;
  while ((*cp) && isspace((int) *cp)) cp++;
  for(i=0; *cp && i<SEQ_MAXLEN; i++) {
    SEQ_REALLOC(sp, hp, i);
    *hp++ = *cp++;
  }
  while (i>0 && isspace((int) *(hp-1))) {
    hp--;
    i--;
  }
  *hp = '\0';
  if (i>=SEQ_MAXLEN) return ERRCODE_SEQLEN;
  sp->size = i;
  
  return ERRCODE_SUCCESS;
}

static int cropSeq(SEQSEQ *sp, SEQLEN_t s, SEQLEN_t e)
     /**< Crop sequence in-place to new boundaries.
      * \param sp Sequence to be cropped.
      * \param s Start of the new sequence (starting from 0)
      * \param e End of the new sequence (starting from 0)
      */
{
  if (sp->code == SEQCOD_COMPRESSED) 
    return ERRCODE_SEQCODE;
  if (e < s)
    return ERRCODE_ARGRANGE;
  if (sp->size < 2)
    return ERRCODE_SUCCESS;

  if (e >= sp->size) e = sp->size-1;
  sp->size = e - s + 1;

  if (s > 0) {
    char *cp = sp->basep;
    for (; s <= e; s++)
      *cp++ = sp->basep[s];
  }
  sp->basep[e+1] = '\0';

  return ERRCODE_SUCCESS;
}

static int appendSeqSegment(SEQSEQ *top, 
			    SEQLEN_t tcpos[MAXN_TERMCHAR], UCHAR_t *ntc,
			    const SEQSEQ *fromp, 
			    SEQLEN_t start, SEQLEN_t length, 
			    char reverse, char withTerm, 
			    const SeqCodec *codep)
     /** Make a copy of a sequence segment and append it to existing sequence.
      * The sequence must not be compressed.
      *
      * \param top Copy of the sequence
      * \param tcpos returns positions of terminating '\0' ( can be NULL)
      * \param ntc returns the number of terminating '\0' encountered (can be NULL)
      * \param fromp  Source of the sequence segment.
      * \param start  Start position of the segment.
      *               Counting is in bytes and starts from 0.
      * \param length Length of the segment in bytes. If 0
      *               Copy to the end of the sequence.
      * \param reverse If != 0 reverse the segment.
      * \param withTerm If != 0 append after terminating character ('\0')
      * \param codep  En/Decoder. If != NULL and reverse != 0
      *               (fromp->code == SEQCOD_MANGLED) append
      *               reverse complement.
      */
{
  char *hp;
  unsigned char uc;
  const char *cp;
  SEQLEN_t i;
  UCHAR_t nt = 0;
  size_t totsize;
  
  if (fromp->code == SEQCOD_COMPRESSED) return ERRCODE_SEQCODE; 
      
  if (start > fromp->size) return ERRCODE_ARGRANGE;
  if (!(length) || (start + length > fromp->size)) 
    length = fromp->size - start;
  if ((withTerm) && top->size > 0) top->size++;
  if (length + top->size > SETSIZ_MAX)
    return ERRCODE_SEQSETSIZ;
  totsize = length + top->size;
  if (totsize+1 > top->alloc_size) {
    if (reallocSeqBlocks(top, totsize))
      return ERRCODE_NOMEM;
  }
  hp = top->basep + top->size;
  cp = fromp->basep + start;
  if (reverse) {
    cp += length-1;
    if (codep != NULL) { 
      if (fromp->code == SEQCOD_MANGLED) {
	for(i=0;i<length;i++, cp--) {
	  if ((*cp) & SEQCOD_STDNT_TESTBIT) {
	    *hp++ = *cp;
	    if (*cp == SEQCOD_TERM) {
	      if (tcpos && nt < MAXN_TERMCHAR) {
		tcpos[nt++] = i;
	      }
	    }
	  } else {
	    *hp++ = codep->codtab_complement[(*cp)&SEQCOD_STDNT_MASK];
	  }
	}
      } else {
	for(i=0;i<length;i++, cp--) {
	  uc = codep->codtab[(int) *cp];
	  if (uc & SEQCOD_STDNT_TESTBIT) {
	    *hp++ = *cp;
	    if (uc == SEQCOD_TERM) {
	      if (tcpos && nt < MAXN_TERMCHAR) {
		tcpos[nt++] = i;
	      }
	    }
	  } else {
	    *hp++ = codep->decodtab[codep->codtab_complement[uc&SEQCOD_STDNT_MASK]];
	  }
	}
      }
    } else {
      for(i=0;i<length;i++) {
	if (!(*cp) && tcpos && nt < MAXN_TERMCHAR)
	  tcpos[nt++] = i;
	*hp++ = *cp--;
      }
    }
  } else {
    for(i=0;i<length;i++) {
      if (!(*cp) && tcpos && nt < MAXN_TERMCHAR)
	tcpos[nt++] = i;
      *hp++ = *cp++;
    }
  }
  *hp = '\0';
  top->size = top->size + i;
  top->code = fromp->code;
  top->nbit_symb = fromp->nbit_symb;
  if (ntc)  *ntc = nt;

  return ERRCODE_SUCCESS;
}

static int appendMangledToCompressedSeq(SEQSEQ *top, 
					const SEQSEQ *fromp, 
					char withTerm
					)
     /** Append a compressed sequence to another.
      *
      * \param top Sequence to be extended.
      * \param fromp  Sequence to be appended.
      * \param withTerm If != 0 append after terminating code.
      */
{
  const char *fsp;
  int nsu;
  size_t n_alloc;
  SETSIZ_t i, offset, totsiz, n_unit;
  uint32_t to_unit;
  uint32_t *tsp = (uint32_t *) top->basep;
  
  if (fromp->code != SEQCOD_MANGLED) 
    return ERRCODE_SEQCODE; 

  CALC_UNIT_OFFSET(top->size, offset, nsu);

  totsiz = top->size + fromp->size;
  if (top->size < 1) {
    top->code = SEQCOD_COMPRESSED;
    top->nbit_symb = NBITS_ALPHABET;
    to_unit = 0;
  } else {
    if (top->code != SEQCOD_COMPRESSED) 
      return ERRCODE_SEQCODE; 
    to_unit = tsp[offset]>>((nsu-1)*NBITS_ALPHABET);
    if ((to_unit&ALPHABET_MASK) != SEQCOD_TERM)
      return ERRCODE_ASSERT;

    if ((withTerm)) {
      totsiz++;
      nsu--;
      if (nsu < 1) {
	nsu = MAXN_PER_UNIT;
	offset++;
	to_unit = 0;
      }
    } else 
      to_unit >>= NBITS_ALPHABET;
  }

  n_unit = totsiz/MAXN_PER_UNIT + 1;
  n_alloc = n_unit*MAXN_PER_UNIT;
  if (n_alloc >= top->alloc_size &&
      reallocSeqBlocks(top, n_alloc))
    return ERRCODE_NOMEM;
  
  fsp = fromp->basep;
  tsp = (uint32_t *) top->basep;
  for (i=0; i<fromp->size; i++) {
    to_unit = (to_unit<<NBITS_ALPHABET) + (fsp[i]& ALPHABET_MASK);
    if (--nsu < 1) {
      tsp[offset++] = to_unit;
      nsu = MAXN_PER_UNIT;
      to_unit = 0;
    }
  }
  to_unit = (to_unit<<NBITS_ALPHABET) + SEQCOD_TERM;
  if (--nsu < 0)
    return ERRCODE_ASSERT;
  to_unit <<= nsu*NBITS_ALPHABET;
  tsp[offset++] = to_unit;

  if (offset != n_unit)
    return ERRCODE_COMPRESS;

  top->size = totsiz;

  return ERRCODE_SUCCESS;
}

static int reverseComplementSeq(SEQSEQ *sp, const SeqCodec *codecp)
     /** Reverse sequence in-place and, if codep !=NULL, create complement.
      */
{
  unsigned char uc;
  char *cp = sp->basep;
  char *ep = cp + sp->size-1;
  register char tmp;

  if (sp->code == SEQCOD_COMPRESSED) return ERRCODE_SEQCODE;
  if ((codecp)) {
    if (sp->code == SEQCOD_MANGLED) {
      while (cp < ep) {
	tmp = (*cp & SEQCOD_STDNT_TESTBIT)? 
	  *cp: codecp->codtab_complement[(*cp)&SEQCOD_STDNT_MASK];
	*cp++ = (*ep & SEQCOD_STDNT_TESTBIT)? 
	  *ep: codecp->codtab_complement[(*ep)&SEQCOD_STDNT_MASK];
	*ep-- = tmp;
      }
      if (cp == ep && !(*cp & SEQCOD_STDNT_TESTBIT))
	*cp = codecp->codtab_complement[(*cp)&SEQCOD_STDNT_MASK];
    } else { /* SEQCOD_ASCII */
      while (cp < ep) {
	uc = codecp->codtab[(int) *cp];
	tmp = (uc & SEQCOD_STDNT_TESTBIT)? 
	  *cp: codecp->decodtab[codecp->codtab_complement[uc&SEQCOD_STDNT_MASK]];
	uc = codecp->codtab[(int) *ep];
	*cp++ = (uc & SEQCOD_STDNT_TESTBIT)? 
	  *ep: codecp->decodtab[codecp->codtab_complement[uc&SEQCOD_STDNT_MASK]];
	*ep-- = tmp;
      }
      if (cp == ep) {
	uc = codecp->codtab[(int) *cp];
	if (!(uc & SEQCOD_STDNT_TESTBIT))
	  *cp = codecp->decodtab[codecp->codtab_complement[uc&SEQCOD_STDNT_MASK]];
      }
    }
  } else { /* just reverse - don't generate complement */
    while (cp < ep) {
      tmp = *cp;
      *cp++ = *ep;
      *ep-- = tmp;
    }
  }
  return ERRCODE_SUCCESS;
}

static int readHeader(SEQSEQ *sp,
#ifdef HAVE_ZLIB
		      gzFile fp, 
#else
		      FILE *fp,
#endif
		      int *prompt, char bufp[LINBUFSIZ])
     /**< Read the FASTA/FASTQ format header into sequence.
      * Collapse all whitespace into single blanks and strip
      * leading and trailing blanks.
      * \param[out] sp Target sequence.
      * \param[in] fp Data source.
      * \param[out] prompt Prompt encountered.
      *
      * The entire header line is read into the sequence. White space
      * is collapsed into space (' ') chars. The string sp->basep is
      * terminated.
      *
      * \note Reads line by line and assumes that first line is already read into bufp.
      * Expects a prompt as next non-white space character and reads to end of line from there.
      */
{
  int errcode = ERRCODE_SUCCESS;
  char *cp, *top, *nextp;
  BOOL_t was_space = TRUE;
  BOOL_t eol_flag = FALSE;
  SEQLEN_t i;


  blankSeq(sp);
  sp->code = SEQCOD_ASCII;
  sp->nbit_symb = NBITS_PER_BYTE;
  *prompt = (int) '\0';
  
  top = sp->basep;
  nextp = bufp;
  i = 0;

  for(; !eol_flag && nextp;
      nextp = 
#ifdef HAVE_ZLIB
	gzgets(fp, bufp, LINBUFSIZ)
#else
	fgets(bufp, LINBUFSIZ, fp)
#endif
      ) {
    for(cp = nextp; *cp && !eol_flag; cp++) {
      if (was_space) {
	if (isspace((int) *cp)) {
	  eol_flag = *cp == '\n' && (*prompt);
	  continue;
	}
	if (!*prompt) {
	  if (*cp != FASTA_PROMPT &&
	      *cp != FASTQ_PROMPT_SEQ &&
	      *cp != FASTQ_PROMPT_QUAL)
	    return ERRCODE_FASTA;
	  *prompt = *cp;
	  continue;
	}
	was_space = FALSE;
      } else if (isspace((int) *cp)) {
	if ((eol_flag = *cp == '\n' && (*prompt)))
	  continue;
	was_space = TRUE;
      }

      SEQ_REALLOC(sp, top, i);
      *top++ = *cp;
      i++;
      if (i>=SEQ_MAXLEN) return ERRCODE_SEQLEN;
    }
  }
  if (was_space && i > 0) {
    top--;
    i--;
  }
  *top= '\0';
  sp->size = i;
  
  if (!nextp) {
    errcode = 
#ifdef HAVE_ZLIB
      (gzeof(fp))? ERRCODE_EOF:ERRCODE_FILEIO;
#else
    (ferror(fp))? ERRCODE_FILEIO: ERRCODE_EOF;
#endif
  }
  
  return errcode;
}

static SEQLEN_t curtailSeqAtFirstSpace(SEQSEQ *sp)
     /**< Terminate sequence at the first space.
      * This can be used to cut off the rest of a header line after the first
      * white space.*/
{
  SEQLEN_t i;
  for (i=0; i<sp->size; i++)
    if (sp->basep[i] == ' ') {
      sp->basep[i] = '\0';
      sp->size = i;
      break;
    }
  return sp->size;
}

static int readSeq(SEQSEQ *sp, 
#ifdef HAVE_ZLIB
		   gzFile fp, 
#else
		   FILE *fp,
#endif
		   int *prompt)
     /**< Read a chunk of sequence of bases or quality factors 
      * from a FASTA/FASTQ file in plain ASCII. Skip all white
      * space. Stop at a haeder prompt preceeded by whitespace
      * including a new line or stop if EOF. The prompt is returned 
      * to the input stream.
      */
{
  char *cp;
  int c;
  BOOL_t was_newline = FALSE;
  SEQLEN_t i;

  blankSeq(sp);
  sp->code = SEQCOD_ASCII;
  sp->nbit_symb = NBITS_PER_BYTE;
  *prompt = 0;
  cp = sp->basep;
  i=0;

#ifdef HAVE_ZLIB
  for (c=gzgetc(fp);
       c != EOF && i<SEQ_MAXLEN;
       c = gzgetc(fp)) {
#else
  for (c=fgetc(fp);
       c != EOF && i<SEQ_MAXLEN;
       c = fgetc(fp)) {
#endif
    if(isspace(c)) {
      if (c == '\n') was_newline = TRUE;
      continue;
    }
    if (was_newline &&
	(c == FASTA_PROMPT ||
	 c == FASTQ_PROMPT_SEQ ||
	 c == FASTQ_PROMPT_QUAL)) {
      /* prompt preceeded by newline */
      *prompt = c;
#ifdef HAVE_ZLIB
      gzungetc(c, fp);
#else
      ungetc(c, fp);
#endif

      break;
    }
 
    SEQ_REALLOC(sp, cp, i);

    *cp++ = c;
    i++;
  }

  *cp = '\0';
  if (i>SEQ_MAXLEN) return ERRCODE_SEQLEN;
  sp->size = i;
  return (c == EOF)? ERRCODE_EOF: ERRCODE_SUCCESS;
}

static int readSeqFast(SEQSEQ *sp,
#ifdef HAVE_ZLIB
		       gzFile fp, 
#else
		       FILE *fp,
#endif

		       int *prompt, char bufp[LINBUFSIZ], SEQLEN_t minlen)
     /**< Like readSeq but read one line at the time
      * on return, the next header line is in bufp. 
      * \param sp Sequence structure.
      * \param fp Input file.
      * \param prompt Returns the FASTA/FASTQ prompt of the line.
      * \param bufp Line buffer.
      * \param minlen Minimum number of characters to be read.
      */
{
  int errcode = ERRCODE_SUCCESS;
  char *top, *cp, *nextp;
  BOOL_t was_newline = FALSE;
  BOOL_t eos_flag = FALSE;
  SEQLEN_t i;

  blankSeq(sp);
  sp->code = SEQCOD_ASCII;
  sp->nbit_symb = NBITS_PER_BYTE;
  *prompt = 0;
  top = sp->basep;
  i=0;

  /* read at least one line */
  nextp = bufp;
  do {
    for (cp=nextp; *cp && !eos_flag; cp++) {
      if (isspace((int) *cp)) {
	was_newline = *cp == '\n';
      } else {
	if (was_newline) {
	  if (i>=minlen &&
	      (*cp == FASTA_PROMPT ||
	       *cp == FASTQ_PROMPT_SEQ ||
	       *cp == FASTQ_PROMPT_QUAL)) {
	    *prompt = *cp;
	    eos_flag = TRUE;
	    continue;
	  } 
	  was_newline = FALSE;
	}
	SEQ_REALLOC(sp, top, i);
	*top++ = *cp;
	i++;
	if (i>SEQ_MAXLEN) return ERRCODE_SEQLEN;
      }
    }
  } while (!eos_flag && 
	   (nextp = 
#ifdef HAVE_ZLIB
	    gzgets(fp, bufp, LINBUFSIZ)
#else
	    fgets(bufp, LINBUFSIZ, fp)
#endif
	    ));
  *top = '\0';
  if (i>SEQ_MAXLEN) return ERRCODE_SEQLEN;
  sp->size = i; /* sequence length excluding term char */

  if (!nextp)
#ifdef HAVE_ZLIB
    errcode = (gzeof(fp))? ERRCODE_EOF: ERRCODE_FILEIO;
#else
    errcode = (ferror(fp))? ERRCODE_FILEIO: ERRCODE_EOF;
  clearerr(fp);
#endif
 
  return errcode;
}

static int reverseSeqInPlace(SEQSEQ *sp)
     /** reverse the sequence */
{
  char *cp, *tp;
  register char tmp;

  if (sp->code == SEQCOD_COMPRESSED) return ERRCODE_SEQCODE;
  cp = sp->basep;
  tp = cp + sp->size - 1;
  while(tp > cp) {
    tmp = *cp;
    *cp++ = *tp;
    *tp-- = tmp;
  }
  return ERRCODE_SUCCESS;
}

static int complementAsciiSeqInPlace(SEQSEQ *sp, const SeqCodec *codep)
     /**< generate reverse complement */
{
  char *cp;
  unsigned char c;

  if (sp->code != SEQCOD_ASCII) return ERRCODE_SEQCODE;
  for (cp = sp->basep; *cp; cp++) {
    c = (unsigned char) codep->codtab[(int) *cp];
    *cp = (c&SEQCOD_STDNT_TESTBIT)? 
      CODEC_UNKNOWN_CHAR: codep->alphabet[(~c)&SEQCOD_STDNT_MASK];
  }
  return ERRCODE_SUCCESS;
}

static int encodeSeq(SEQSEQ *sp, const SeqCodec *codep)
     /**< Encode sequence in-place. Read through terminating '\0' until
      * sp->size characters are converted */
{
  char *cp;
  SEQLEN_t i;

  if (sp->code != SEQCOD_ASCII) 
    return (sp->code == SEQCOD_MANGLED)? 
      ERRCODE_SUCCESS: ERRCODE_SEQCODE;
  cp = sp->basep;
  for(i=0; i<sp->size; i++) {
    if (cp[i])
      cp[i] = (char) codep->codtab[(int) cp[i]];
  }
/*   for(cp = sp->basep; *cp; cp++) { */
/*     *cp = (char) codep->codtab[(int) *cp]; */
/*   } */
  sp->code = SEQCOD_MANGLED;
  return ERRCODE_SUCCESS;
}

static int compressSeq(SEQSEQ *sp)
     /**< Compress sequence in-place to 3 bits per base (10 bases per 32-bit int -
      * The highest two bits are not used). The sequence is terminated
      * with 3 set bits. The 9 - (sequence_length % 10) lower bits of
      * the last 32-bit integer are zeroed. It is assumed that the
      * input sequence is stored in the 'mangled' way (see type SeqCodec).
      * All sp->size characters are compressed, with a terminating '\0'
      * compressed as SEQCOD_TERM.
      * Unused memory is freed.
      */
{
  char *fromp, *hp;
  int nsu;
  uint32_t *top, to_unit;
  size_t unit_ctr;
  SETSIZ_t i;

  if (sp->code != SEQCOD_MANGLED) return ERRCODE_SEQCODE;
  if (sp->alloc_size < NBYTES_PER_INT32+1) {
    hp=EREALLOCP(sp->basep, NBYTES_PER_INT32+1);
    if (!hp) return ERRCODE_NOMEM;
    sp->basep = hp;
    sp->alloc_size = NBYTES_PER_INT32+1; 
  }
  sp->code = SEQCOD_COMPRESSED;
  sp->nbit_symb = NBITS_ALPHABET;
  top = (uint32_t *) sp->basep;
  unit_ctr=1;
  
  fromp = sp->basep;
  for (i=0, nsu=MAXN_PER_UNIT-1, to_unit=0; i<sp->size; i++) {
    if (fromp[i]) {
      to_unit += ((uint32_t) (fromp[i]) & ALPHABET_MASK)<<(nsu*NBITS_ALPHABET);
    } else {
      to_unit += ((uint32_t) SEQCOD_TERM) << (nsu*NBITS_ALPHABET);
    }
    if(--nsu < 0) {
      *top++ = to_unit;
      to_unit=0;
      nsu = MAXN_PER_UNIT-1;
      unit_ctr++;
    }
  }
  to_unit += ((uint32_t) SEQCOD_TERM) <<(nsu*NBITS_ALPHABET);
  *top++ = to_unit;

  if (unit_ctr != sp->size/MAXN_PER_UNIT + 1) {
#ifdef seq_debug
    printf("seq_debug:compressChunk: length=%llu, number of 32-bit ints = %llu\n", 
	   (unsigned long long) sp->size,
	   (unsigned long long) unit_ctr);
#endif
    return ERRCODE_COMPRESS;
  }

  /* free unused memory */
  hp = EREALLOCP(sp->basep, unit_ctr*NBYTES_PER_INT32+1);
  if (!hp) return ERRCODE_NOMEM;
  sp->basep = hp;
  sp->alloc_size = unit_ctr*NBYTES_PER_INT32;
  hp[sp->alloc_size] = '\0';
  sp->code = SEQCOD_COMPRESSED;

  return ERRCODE_SUCCESS;
}

static int writeCompressedSeq(FILE *fp, SEQSEQ *sp)
{
  size_t n_unit;
  if (sp->code != SEQCOD_COMPRESSED) return ERRCODE_SEQCODE;
  n_unit = sp->size/MAXN_PER_UNIT + 1;
  if (fwrite(sp->basep, sizeof(uint32_t), n_unit, fp) != n_unit) 
    return ERRCODE_FAILURE;
  return ERRCODE_SUCCESS;
}
  
static int readCompressedSeq(SEQSEQ *sp, FILE *fp)
     /**< Read a sequence of bases from a binary file in compressed form up
      * to termination code */
{
  char is_terminated;
  uint32_t *cp;
  short i;
  size_t bytctr, nbytes_uint32_t = sizeof(uint32_t);
  SETSIZ_t l;

  is_terminated=0;
  for (l=0, bytctr=0, cp = (uint32_t *) sp->basep; 
       !is_terminated && l < SEQ_MAXLEN; 
       cp++, bytctr += nbytes_uint32_t) {

    if (bytctr >= sp->alloc_size-nbytes_uint32_t) {
      if (reallocSeqBlocks(sp, bytctr + nbytes_uint32_t))
	return ERRCODE_NOMEM;
      cp = (uint32_t *) (sp->basep + bytctr);
    }
 
    if (fread(cp, nbytes_uint32_t, 1, fp) != 1)
      break;

    for (i=MAXN_PER_UNIT-1; i>=0; i--) {
      if ((((*cp>>(i*3))) & ALPHABET_MASK) == SEQCOD_TERM) {
	is_terminated = 1;
	break;
      }
      l++;
    }
  } 
  sp->size = l;
  sp->code = SEQCOD_COMPRESSED;
  if (!is_terminated) return ERRCODE_FILEFORM;

  return (is_terminated)? ERRCODE_SUCCESS: ERRCODE_FILEFORM;
}

static int readCompressedSeqOfKnownLength(SEQSEQ *sp, FILE *fp, SETSIZ_t length)
     /**< Read a sequence of known length in compressed format 
      * length is the number of nucleotides */
{
  int errcode;
  uint32_t *cp;
  SETSIZ_t n_units = (size_t) length/MAXN_PER_UNIT+1;
  short termoffs;

  if ((errcode = reallocSeqBlocks(sp, n_units*NBYTES_PER_INT32)))
    return errcode;
  if (fread(sp->basep, NBYTES_PER_INT32, n_units, fp) != n_units)
    return ERRCODE_FILEFORM;
  /* check if properly terminated */
  termoffs = n_units * MAXN_PER_UNIT - length - 1;
  cp = (uint32_t *) sp->basep;
  if ((((cp[n_units-1]>>(termoffs*3))) & ALPHABET_MASK) != SEQCOD_TERM) 
    return ERRCODE_FILEFORM;
  sp->size = length;
  sp->code = SEQCOD_COMPRESSED;

  return ERRCODE_SUCCESS;
}

static int uncompressSeq(SEQSEQ *ucp, SEQLEN_t tcpos[MAXN_TERMCHAR], UCHAR_t *ntc,
			 const SEQSEQ *sp, 
			 SETSIZ_t start, SETSIZ_t length, 
			 const SeqCodec *codep)
     /**< uncompress a segment of the sequence back to ASCII 
      * \param tcpos Returns position of maximum MAXN_TERMCHAR positions of 
      *        termination characters (can be null) 
      * \param ntc Returns the number of termination codes found (can be null)
      */
{
  char *bufp;
  int nsu;
  UCHAR_t nt = 0;
  const uint32_t *fromp;
  SETSIZ_t ctr, offset;
  unsigned char code;
  
  if (sp->code != SEQCOD_COMPRESSED) return ERRCODE_SEQCODE;
  if (start > sp->size) return ERRCODE_ARGRANGE;
  if (!(length) ||
      start + length > sp->size) length = sp->size-start;
  if (ucp->alloc_size < length + 1 && 
      (reallocSeqBlocks(ucp, length + 1)))
    return ERRCODE_NOMEM;
  bufp  = ucp->basep;
  CALC_UNIT_OFFSET(start, offset, nsu);
  nsu--;
  fromp = (const uint32_t *) sp->basep + offset;

  for (ctr=0; ctr<length; ctr++) {
    code = (unsigned char) (((*fromp)>>(NBITS_ALPHABET*nsu)) & ALPHABET_MASK);
    if (--nsu < 0) {
      fromp++;
      nsu=MAXN_PER_UNIT-1;
    }
    if (code == SEQCOD_TERM) {
      *bufp++ = '\0';
      if (tcpos && nt < MAXN_TERMCHAR)
	tcpos[nt++] = ctr;
    } else if (code < codep->alphlen) {
      *bufp++ = codep->alphabet[code];
    } else {
      *bufp++ = CODEC_UNKNOWN_CHAR;
    }
  }
  *bufp = '\0';
  ucp->size = ctr;
  ucp->code = SEQCOD_ASCII;
  if (ntc) *ntc = nt;

  return ERRCODE_SUCCESS;
}

static int decodeSeq(SEQSEQ *sp, const SeqCodec *codep)
     /**< convert sequence from mangled encoding back to ASCII characters in-place */
{
  unsigned char *cp;
  if (sp->code != SEQCOD_MANGLED) return ERRCODE_SEQCODE;

  for (cp = (unsigned char *) sp->basep; *cp; cp++)
    *cp = codep->decodtab[*cp];

  sp->code = SEQCOD_ASCII;
  return ERRCODE_SUCCESS;
}

static int decodeSeqAsStandardNt(SEQSEQ *dep, 
				 const SEQSEQ *sp, 
				 SETSIZ_t start, 
				 SETSIZ_t length, 
				 const SeqCodec *codep,
				 char as_rcp)
     /**< Decode sequence as ASCII standard nucleotide codes (ACGT).
      * Sequence encoding must be SEQCOD_MANGLED.
      *
      * \param dep Contains segment, contents are overwritten.
      * \param sp  Seqeunce from which segment is taken.
      * \param start Start of the segment in bases (counting from 0).
      * \param length Length of the segment. If == 0 decode
      *               entire sequence.
      * \param codep En/Decoder
      * \param as_rcp Flag, if !=0 return reverse complement.
      *               Start position is given in forward direction. 
      */
{
  char *bufp;
  unsigned char c;
  const unsigned char *fromp, *endp;

  if (sp->code != SEQCOD_MANGLED) 
    return ERRCODE_SEQCODE;
  if (start > sp->size) return ERRCODE_ARGRANGE;
  if (!(length) || length > sp->size) length = sp->size - start;
  if (dep->alloc_size < length + 1 &&
      (reallocSeqBlocks(dep, length+1)))
    return ERRCODE_NOMEM;

  bufp = dep->basep;
  if (as_rcp) {
    endp = (unsigned char *) sp->basep+start;
    fromp = endp + length;
    while(fromp>endp) {
      c = *--fromp;
      *bufp++ = (c&SEQCOD_STDNT_TESTBIT)? 
	CODEC_UNKNOWN_CHAR: codep->alphabet[(~c)&SEQCOD_STDNT_MASK];
    }
  } else {
    fromp = (unsigned char *) sp->basep+start;
    endp=fromp+length;
    while(fromp<endp)
      *bufp++ = codep->alphabet[(*fromp++)&ALPHABET_MASK];
  }
  *bufp = '\0';
  dep->size = length;
  dep->code = SEQCOD_ASCII;
  return ERRCODE_SUCCESS;
}

static int checkSeqNtSymbolsAreLetters(const SEQSEQ *sp)
{
  int errcode = ERRCODE_SUCCESS;

  if (sp->code == SEQCOD_ASCII) {
    SETSIZ_t i = 0;
    for (i=0; i<sp->size && isalpha(sp->basep[i]); i++);
    if (i < sp->size)
      errcode = ERRCODE_SEQNTSMBL;
  }

  return errcode;
}

static int fprintSeqFastqHeader(void *top,
#ifdef HAVE_ZLIB
				BOOL_t is_gzipped,
#endif
				const SEQSEQ *sp, 
				int prompt)
{
#ifdef HAVE_ZLIB
  gzFile gzfp = (is_gzipped)? (gzFile) top: NULL;
  FILE *fp = (is_gzipped)? NULL: (FILE *) top;
#else
  FILE *fp = (FILE *) top;
#endif 

  if (
#ifdef HAVE_ZLIB
      (is_gzipped)? (gzputc(gzfp, prompt) == EOF):
#endif
      (fputc(prompt, fp) == EOF)
      ) return ERRCODE_WRITEERR;
  if (sp != NULL) {
    if (sp->code != SEQCOD_ASCII) 
      return ERRCODE_SEQCODE;
#ifdef HAVE_ZLIB
    if (is_gzipped) {
      gzputs(gzfp, sp->basep);
    } else {
#endif
    fputs(sp->basep, fp);
#ifdef HAVE_ZLIB
    }
#endif
  }
#ifdef HAVE_ZLIB
  if (is_gzipped) {
    gzputc(gzfp, '\n');
  } else {
#endif
    fputc('\n', fp);
#ifdef HAVE_ZLIB
  }
#endif

  return ERRCODE_SUCCESS;
}

static int fprintSeqFastqSequence(void *top,
#ifdef HAVE_ZLIB
				  BOOL_t is_gzipped,
#endif
				  const SEQSEQ *sp, short linewidth)
     /**< Write sequence to stream. Sequence has to be sp->code == SEQCOD_ASCII
      *
      * \param fp Output stream.
      * \param sp Sequence to be written to stream.
      * \param linewidth Linewidth for output
      */
{
  char buf[OUTBUFSIZE+1];
  SETSIZ_t pos;
#ifdef HAVE_ZLIB
  gzFile gzfp = (is_gzipped)? top: NULL;
  FILE *fp = (is_gzipped)? NULL: (FILE *) top;
#else
  FILE *fp = (FILE *) top;
#endif  


  if (sp->code != SEQCOD_ASCII) return ERRCODE_SEQCODE;
  if (linewidth>OUTBUFSIZE) linewidth = OUTBUFSIZE;
  if (linewidth<1) {
#ifdef HAVE_ZLIB
    if (is_gzipped) {
      gzputs(gzfp, sp->basep);
      if (gzputc(gzfp,'\n') == EOF)
	return ERRCODE_WRITEERR;
    } else {
#endif
      fputs(sp->basep, fp);
      if (fputc('\n', fp) == EOF)
      return ERRCODE_WRITEERR;
#ifdef HAVE_ZLIB
    }
#endif
  } else {
    for (pos=0; pos<sp->size; pos+=linewidth) {
      strncpy(buf, sp->basep+pos, linewidth);
      buf[linewidth] = '\0';
#ifdef HAVE_ZLIB
      if (is_gzipped) {
	gzputs(gzfp, buf);
	if (gzputc(gzfp,'\n') == EOF)
	  return ERRCODE_WRITEERR;
      } else {
#endif
	fputs(buf, fp);
	if (fputc('\n', fp) == EOF) 
	  return ERRCODE_WRITEERR;
#ifdef HAVE_ZLIB
      }
#endif
    }
  }

  return ERRCODE_SUCCESS;
}

   
/******************************************************************************
 ************************ Private Methods of Type SeqFastq ********************
 ******************************************************************************/

static int readQual(SeqFastq *sqp, SeqIO *ifp)
{
  int this_prompt;

  /* read sequence of quality factors if present */
  if (sqp->qheadp == NULL && 
      (sqp->qheadp = createSeq(sqp->headp->block_size)) == NULL)
    return ERRCODE_NOMEM;

  if ((ifp->status = readHeader(sqp->qheadp, ifp->fp, 
				&this_prompt, ifp->linbufp))) 
    return ifp->status;
 
  if (this_prompt != FASTQ_PROMPT_QUAL)  /* expect a FASTQ prompt */
    return (ifp->status = ERRCODE_FASTA);
    
  if (sqp->qualp == NULL &&
      (sqp->qualp = createSeq(sqp->datap->block_size)) == NULL)
    return ifp->status = ERRCODE_NOMEM;

  if (sqp->datap->size > SEQLEN_MAX)
    return ERRCODE_FASTA;

  ifp->status = readSeqFast(sqp->qualp, ifp->fp,
			    &this_prompt, ifp->linbufp, sqp->datap->size);

  sqp->type = SEQTYP_FASTQ;

  return (ifp->status == ERRCODE_EOF)? ERRCODE_SUCCESS: ifp->status;
}
/******************************************************************************
 ************************ Public Methods of Type SeqFastq *********************
 ******************************************************************************/

SeqFastq *seqFastqCreate(int blocksize, char type)
{
  SeqFastq *sqp=NULL;
  EMALLOCP0(sqp);
  if (sqp == NULL) return NULL;
  sqp->type = (type == SEQTYP_FASTA || type == SEQTYP_FASTQ)? type: SEQTYP_UNKNOWN;
  if (blocksize < 1) blocksize = SEQ_BLKSZ_DEFAULT;
  if (!(sqp->headp = createSeq(BLOCKSIZE_HEADER)) ||
      !(sqp->datap = createSeq(blocksize)) ||
      (type == SEQTYP_FASTQ &&
       !(sqp->qualp = createSeq(blocksize)))) {
    seqFastqDelete(sqp);
    sqp = NULL;
  }  
  return sqp;
}

void seqFastqDelete(SeqFastq *sqp)
{
  if (sqp) {
    deleteSeq(sqp->headp);
    deleteSeq(sqp->datap);
    deleteSeq(sqp->qheadp);
    deleteSeq(sqp->qualp);
  }
  free(sqp);
}

void seqFastqBlank(SeqFastq *p)
{
  if (p) {
    blankSeq(p->datap);
    blankSeq(p->qualp);
    blankSeq(p->headp);
    blankSeq(p->qheadp);
  }
}

void seqFastqFreeUnusedMem(SeqFastq *sqp)
{
  if (sqp->headp) finaliseSeq(sqp->headp);
  if (sqp->datap) finaliseSeq(sqp->datap);
  if (sqp->qualp) finaliseSeq(sqp->qualp);
  if (sqp->qheadp) finaliseSeq(sqp->qheadp);
}

int seqFastqCheck(const SeqFastq *sqp)
{
  if (sqp == NULL ||
      sqp->headp == NULL ||
      sqp->datap == NULL)
    return ERRCODE_FAILURE;
  if (sqp->type == SEQTYP_FASTQ &&
      (sqp->qualp == NULL || sqp->datap->size != sqp->qualp->size))
    return ERRCODE_FAILURE;
  
  return checkSeqNtSymbolsAreLetters(sqp->datap);
}

int seqFastqSetType(SeqFastq *sqp, char type)
{
  int errcode = ERRCODE_SUCCESS;

  if (type != SEQTYP_FASTQ && type != SEQTYP_FASTA)
    type = SEQTYP_UNKNOWN;

  if (type != SEQTYP_FASTA && !(sqp->qualp)) {
    SETSIZ_t i = 0LL;
    if (!sqp->datap)
      return ERRCODE_ASSERT;
    if (!(sqp->qualp = createSeq(sqp->datap->block_size)))
      return ERRCODE_NOMEM;
    if (((size_t) sqp->qualp->block_size) < sqp->datap->size &&
	!(errcode = reallocSeqBlocks(sqp->qualp, sqp->datap->size)))
      return errcode;
    for (i=0; i<sqp->datap->size; i++)
      sqp->qualp->basep[i] = SEQCOD_QVAL_OFFS;
  }
  sqp->type = type;

  return errcode;   
}

int seqFastqSetAscii(SeqFastq *sqp,
		     const char *name, const char *seqp, 
		     const char *name_qual, const char *qualp)
{
  int errcode;

  if ((name) && (errcode = setSeq(sqp->headp, name)))
    return errcode;

  if ((seqp) && (errcode = setSeq(sqp->datap, seqp)))
    return errcode;
  
  if ((name_qual)) {
    if (!(sqp->qheadp ||
	  (sqp->qheadp = createSeq(BLOCKSIZE_HEADER))))
      return ERRCODE_NOMEM;
    if ((errcode = setSeq(sqp->qheadp, name_qual)))
      return errcode;
  }

  if ((qualp) && qualp[0] != '\0') {
    sqp->type = SEQTYP_FASTQ;
    if (!(sqp->qualp ||
	  (sqp->qualp = createSeq(sqp->datap->block_size))))
      return ERRCODE_NOMEM;
    if ((errcode = setSeq(sqp->qualp, qualp)))
      return errcode;
 
    if (sqp->datap->size != sqp->qualp->size) return ERRCODE_QUALLEN;
  } else if (sqp->type != SEQTYP_FASTQ || !(sqp->qualp) ||
	     sqp->datap->size != sqp->qualp->size) {
    sqp->type = SEQTYP_FASTA;
  }

  return ERRCODE_SUCCESS;
}

int seqFastqSetQual(SeqFastq *sqp, const char qval)
{
  SETSIZ_t i;

  if (!(sqp->qualp ||
	(sqp->qualp = createSeq(sqp->datap->block_size))))
    return ERRCODE_NOMEM;

  if (sqp->datap->size >= sqp->qualp->alloc_size ||
      reallocSeqBlocks(sqp->qualp, sqp->datap->size))
    return ERRCODE_NOMEM;

  for (i=0; i<sqp->datap->size; i++)
    sqp->qualp->basep[i] = qval;
  sqp->qualp->basep[i] = '\0';
  sqp->qualp->size = i;
  
  sqp->type = SEQTYP_FASTQ;

  return ERRCODE_SUCCESS;
}
  
int seqFastqAppendSegment(SeqFastq *top, const SeqFastq *fromp, 
			  SEQLEN_t start, SEQLEN_t length, 
			  char reverse, const SeqCodec *codep)
{
  int errcode; 
  if ((errcode = appendSeqSegment(top->datap, NULL, NULL, fromp->datap, start, length,
				  reverse, 0, codep)))
    return errcode;
  /* copy the name if blank */
  if (!top->headp->size &&
      (errcode = appendSeqSegment(top->headp, NULL, NULL, fromp->headp, 0, 0, 0, 0, NULL)))
      return errcode;

  if (top->type == SEQTYP_UNKNOWN) top->type = fromp->type; 

  if (fromp->qualp == NULL ||
      fromp->qualp->size < 1) 
    return errcode;
  
  if (top->qualp == NULL &&
      ((top->qualp = createSeq(top->datap->block_size)) == NULL))
    return ERRCODE_NOMEM;

  if ((errcode = appendSeqSegment(top->qualp, NULL, NULL, fromp->qualp, 
				  start, length,
				  reverse, 0, NULL)))
    return errcode;
  top->type = SEQTYP_FASTQ;

  return (top->qualp->size != top->datap->size)?
    ERRCODE_QUALLEN: ERRCODE_SUCCESS;
}

int seqFastqReverse(SeqFastq *sqp, const SeqCodec *codecp)
{
  int errcode = reverseComplementSeq(sqp->datap, codecp);
  if (!(errcode) && sqp->qualp != NULL && sqp->type == SEQTYP_FASTQ)
    errcode = reverseComplementSeq(sqp->qualp, NULL);
  return errcode;
}

int seqFastqRead(SeqFastq *sqp, SeqIO *ifp)
{
  int errcode;
  int this_prompt = FASTQ_PROMPT_QUAL;
  int next_prompt;
  if (ifp->mode != SEQIO_READ) return ERRCODE_NOREAD;
  if (ifp->status) return ifp->status;

  /* skip quality prompts */
  while (this_prompt == FASTQ_PROMPT_QUAL) {

    if ((ifp->status = readHeader(sqp->headp, ifp->fp,
				  &this_prompt, ifp->linbufp))) 
      return ifp->status;
    ifp->status = readSeqFast(sqp->datap, ifp->fp,
			      &next_prompt, ifp->linbufp, 0);
  }

  if (ifp->status ||
      next_prompt != FASTQ_PROMPT_QUAL) {
    sqp->type = SEQTYP_FASTA;
    return (ifp->status == ERRCODE_EOF)? ERRCODE_SUCCESS: ifp->status;
  }
  errcode = readQual(sqp, ifp);

  if (!(errcode) && sqp->type == SEQTYP_FASTQ && 
      sqp->datap->size != sqp->qualp->size)
    errcode = ERRCODE_FASTA;

  return errcode;
}

int seqFastqFind(SeqFastq *sqp, const char *nam, SeqIO *ifp)
{
  int this_prompt = FASTQ_PROMPT_QUAL;
  int next_prompt;
  int cmp = 0;

  if (ifp->mode != SEQIO_READ) return ERRCODE_NOREAD;

  while (!(ifp->status = scrollToHeaderLine(ifp->fp, &this_prompt, ifp->linbufp)) && (cmp)) {
    if ((ifp->status = readHeader(sqp->headp, ifp->fp,
				  &this_prompt, ifp->linbufp)))
      break;
    curtailSeqAtFirstSpace(sqp->headp);
    cmp = strcmp(sqp->headp->basep, nam);
  }
  if (ifp->status) 
    return ifp->status;
  if (cmp) return ERRCODE_FAILURE;

  ifp->status = readSeq(sqp->datap, ifp->fp,
			&next_prompt);
  if (ifp->status ||
      sqp->type == SEQTYP_FASTA || 
      next_prompt != FASTQ_PROMPT_QUAL) {
    sqp->type = SEQTYP_FASTA;
    return (ifp->status == ERRCODE_EOF)? ERRCODE_SUCCESS: ifp->status;
  }

  return readQual(sqp, ifp);
}

int seqFastqWrite(SeqIO *ofp, const SeqFastq *sqp, short linewidth)
{
  int cprompt;
  if (sqp == NULL) return ERRCODE_FAILURE;
  if (ofp->mode ==  SEQIO_READ) return ERRCODE_NOWRITE;

  if (ofp->status) return ofp->status;
    
  if (sqp->type == SEQTYP_FASTA) {
    cprompt = FASTA_PROMPT;
  } else if (sqp->type == SEQTYP_FASTQ) {
    cprompt = FASTQ_PROMPT_SEQ;
  } else {
    return ERRCODE_SEQTYP;
  }

  ofp->status = fprintSeqFastqHeader(ofp->fp, 
#ifdef HAVE_ZLIB
				     (ofp->flags & SEQIOFLG_GZIP) != 0,
#endif
				     sqp->headp, cprompt);
  if (!ofp->status) 
    ofp->status = fprintSeqFastqSequence(ofp->fp, 
#ifdef HAVE_ZLIB
					 (ofp->flags & SEQIOFLG_GZIP) != 0,
#endif
					 sqp->datap, linewidth);
  if (ofp->status || sqp->type != SEQTYP_FASTQ || sqp->qualp == NULL)
    return ofp->status;

  ofp->status = fprintSeqFastqHeader(ofp->fp, 
#ifdef HAVE_ZLIB
				     (ofp->flags & SEQIOFLG_GZIP) != 0,
#endif
				     sqp->qheadp, FASTQ_PROMPT_QUAL);
  if (!ofp->status) 
    ofp->status = fprintSeqFastqSequence(ofp->fp, 
#ifdef HAVE_ZLIB
					 (ofp->flags & SEQIOFLG_GZIP) != 0,
#endif
					 sqp->qualp, linewidth);
  return ofp->status;
}

int seqFastqWriteCompressedToFile(FILE *fp, const SeqFastq *sqp)
{
  if (!sqp) return ERRCODE_NULLPTR;
  return writeCompressedSeq(fp, sqp->datap);
}
  
/* Accessors of type SeqFastq */
const char *seqFastqGetConstSequence(const SeqFastq *sfqp, 
				     SEQLEN_t *length, 
				     char *codtyp) 
{
    if (!sfqp) {
      if (length)
	*length = 0;
	return NULL;
    }
    if ((length) && sfqp->datap->size < SEQLEN_MAX) 
      *length = (SEQLEN_t) sfqp->datap->size;
    if (codtyp) *codtyp = sfqp->datap->code;
    return sfqp->datap->basep;
}

char *seqFastqGetSequence(SeqFastq *sfqp, 
			  SEQLEN_t *length, 
			  char *codtyp) 
{
    if (!sfqp) {
      if (length)
	*length = 0;
	return NULL;
    }
    if (length && sfqp->datap->size < SEQLEN_MAX) 
      *length = (SEQLEN_t) sfqp->datap->size;
    if (codtyp) *codtyp = sfqp->datap->code;
    return sfqp->datap->basep;
}

const char *seqFastqGetConstQualityFactors(const SeqFastq *sfqp, SEQLEN_t *length, char *code)
{
  if (!sfqp) return NULL;
  if (!sfqp->qualp) {
    if (length) *length = 0;
    return NULL;
  }
  if (length && sfqp->qualp->size < SEQLEN_MAX) 
    *length = (SEQLEN_t) sfqp->qualp->size;
  if (code) *code = sfqp->qualp->code;
  return sfqp->qualp->basep;
}

char *seqFastqGetQualityFactors(SeqFastq *sfqp, SEQLEN_t *length, char *code)
{
  if (!sfqp) return NULL;
  if (!sfqp->qualp) {
    if (length) *length = 0;
    return NULL;
  }
  if (length && sfqp->qualp->size < SEQLEN_MAX) 
     *length = (SEQLEN_t) sfqp->qualp->size;
  if (code) *code = sfqp->qualp->code;
  return sfqp->qualp->basep;
}

const char *seqFastqGetSeqName(const SeqFastq *sfqp)
{
  return (sfqp)? sfqp->headp->basep: NULL;
}

void seqFastqCurtailSeqName(SeqFastq *sfqp)
{
  curtailSeqAtFirstSpace(sfqp->headp);
}

const char *seqFastqGetQualName(const SeqFastq *sfqp)
{
  return (sfqp->qheadp)? sfqp->qheadp->basep: NULL;
}

int seqFastqEncode(SeqFastq *sqp, const SeqCodec *codep)
{
  return encodeSeq(sqp->datap, codep);
}

int seqFastqDecode(SeqFastq *sqp, const SeqCodec *codep)
{
  return decodeSeq(sqp->datap, codep);
}

int seqFastqCompress(SeqFastq *sqp)
{
  if (sqp->datap->code != SEQCOD_MANGLED) return ERRCODE_SEQCODE;
  return compressSeq(sqp->datap);
}


int seqFastqUncompress(SeqFastq *ucp, const SeqFastq *sqp,
		       SEQLEN_t start, SEQLEN_t length, 
		       const SeqCodec *codep, char as_rcp)
{
  int errcode;
  if ((errcode = uncompressSeq(ucp->datap, NULL, NULL, sqp->datap, start, length, codep)))
    return errcode;
  ucp->type = SEQTYP_FASTA;
  blankSeq(ucp->headp);
  if ((errcode = appendSeqSegment(ucp->headp, NULL, NULL, sqp->headp, 0, 0, 0, 0, NULL)))
    return errcode;
  if (as_rcp && 
      ((errcode = reverseSeqInPlace(ucp->datap)) ||
       (errcode = complementAsciiSeqInPlace(ucp->datap, codep))))
    return errcode;

  if (sqp->qualp == NULL || sqp->qualp->size <= 0)
    return ERRCODE_SUCCESS;
  if (!(ucp->qualp || (ucp->qualp = createSeq(sqp->qualp->block_size))))
    return ERRCODE_NOMEM;
  ucp->type = SEQTYP_FASTQ;
  blankSeq(ucp->qualp);
  errcode = appendSeqSegment(ucp->qualp, NULL, NULL, sqp->qualp, start, length, 0, 0, NULL);
  if (!errcode && as_rcp)
    errcode = reverseSeqInPlace(ucp->qualp);
  return errcode;
}

int seqFastqReadCompressedBinary(SeqFastq *sqp, FILE *fp, const char *label)
{
  int errcode = readCompressedSeq(sqp->datap, fp);
  if (!errcode && label != NULL)
    errcode = setSeq(sqp->headp, label);
  deleteSeq(sqp->qualp);
  sqp->qualp = NULL;
  sqp->type = SEQTYP_FASTA;
  return errcode;
}

int seqFastqReadCompressedBinaryOfKnownLength(SeqFastq *sqp, FILE *fp, SEQLEN_t len, 
					      const char *label)
{
  int errcode = readCompressedSeqOfKnownLength(sqp->datap, fp, len);
  if (!errcode && label != NULL)
    errcode = setSeq(sqp->headp, label);
  deleteSeq(sqp->qualp);
  sqp->qualp = NULL;
  sqp->type = SEQTYP_FASTA;
  return errcode;
}


int seqFastqDecodeAsStandardNt(SeqFastq *dep, const SeqFastq *sqp, 
			       SEQLEN_t start, SEQLEN_t length, 
			       const SeqCodec *codep,
			       char as_rcp)
     /* threre is still a problem for as_rcp != 0 */
{
  int errcode;
  seqFastqBlank(dep);
  errcode = decodeSeqAsStandardNt(dep->datap, sqp->datap, start, length, codep, as_rcp);
  if (!errcode && dep->type == SEQTYP_FASTQ) {
    if (sqp->type == SEQTYP_FASTQ)
      errcode = appendSeqSegment(dep->qualp, NULL, NULL, sqp->qualp, start, length, as_rcp, 0, NULL);
    else 
      errcode = ERRCODE_SEQTYP;
  }
  return errcode;
}

/******************************************************************************
 ************************ Private Methods of Type SeqSet **********************
 ******************************************************************************/

static int reallocSeqSet(SeqSet *ssp, SEQNUM_t num)
{
  SETSIZ_t nsiz;
  void *hp;

  if (num < 1) return ERRCODE_ARGRANGE;
  nsiz = (((num-1)/ssp->blocksiz)+1)*ssp->blocksiz;
  if (nsiz>INT_MAX-1) return ERRCODE_OVERFLOW;
  hp = EREALLOCP(ssp->sop, nsiz);
  if (!hp) return ERRCODE_NOMEM;
  ssp->sop = (SETSIZ_t *) hp;

  hp = EREALLOCP(ssp->namoffs, nsiz);
  if (!hp) return ERRCODE_NOMEM;
  ssp->namoffs = (SETSIZ_t *) hp;

  ssp->n_alloc = nsiz;

  return ERRCODE_SUCCESS;
}

static int reallocSeqSetName(SeqSet *ssp, size_t len)
{
  size_t nsiz;
  char *hp;

  if (len < 1) 
    return ERRCODE_ARGRANGE;
  nsiz = (((len-1)/ssp->nam_blocksiz)+1)*ssp->nam_blocksiz;
  
  if (nsiz > INT_MAX) 
    return ERRCODE_OVERFLOW;
			
  hp = EREALLOCP(ssp->namebasep, nsiz);
  if (!hp) return ERRCODE_NOMEM;

  ssp->namebasep = hp;
  ssp->nam_alloc = nsiz;

  return ERRCODE_SUCCESS;
}

/******************************************************************************
 ************************ Public Methods of Type SeqSet ***********************
 ******************************************************************************/
SeqSet *seqSetCreate(int blocksiz, UCHAR_t flags)
{
  SeqSet *ssp;

  if (blocksiz < 1) blocksiz = SEQSET_BLOCKSIZ_DEFAULT;
  EMALLOCP0(ssp);
  if (!ssp) return 0;

  ssp->statusflag = flags;
  ssp->sqp = createSeq(blocksiz);
  if ((flags & SEQSET_BASQUAL))
    ssp->qqp = createSeq(blocksiz);
    
  ECALLOCP(blocksiz, ssp->sop);
  ECALLOCP(blocksiz, ssp->namoffs);
  ECALLOCP(SEQSET_NAMBLOCKSIZ_DEFAULT, ssp->namebasep);

  if (!(ssp->sqp) || ((flags & SEQSET_BASQUAL) && !(ssp->qqp)) || 
      !(ssp->sop) || !(ssp->namoffs) || !(ssp->namebasep)) {
    seqSetDelete(ssp);
    return 0;
  }
  ssp->nam_blocksiz = SEQSET_NAMBLOCKSIZ_DEFAULT;
  ssp->n_alloc = ssp->blocksiz = blocksiz;
  ssp->nam_alloc = ssp->nam_blocksiz;
  
  return ssp;
}

void seqSetDelete(SeqSet *ssp)
{
  if (ssp) {
    deleteSeq(ssp->sqp);
    deleteSeq(ssp->qqp);
    free(ssp->sop);
    free(ssp->namoffs);
    free(ssp->namebasep);
    free(ssp->sxp);
  }
  free(ssp);
}

void seqSetBlank(SeqSet *ssp)
{
  blankSeq(ssp->sqp);
  blankSeq(ssp->sqp);
  ssp->n_seq = 0;
}

int seqSetAddSequence(SeqSet *ssp, const SeqFastq *sqp)
{
  int errcode;
  char *cp;
  SETSIZ_t namsiz;
 
  if (sqp->datap->code == SEQCOD_COMPRESSED) 
    return ERRCODE_SEQCODE;

  if ((ssp->statusflag&SEQSET_COMPRESSED)) {
    if (sqp->datap->code != SEQCOD_MANGLED)
      return ERRCODE_SEQCODE;
  } else {
    if ((ssp->n_seq) && sqp->datap->code != ssp->sqp->code)
      return ERRCODE_SEQCODE;
  }

  if (ssp->n_seq+1 >= (SEQNUM_t) SEQNUM_MAX)
    return ERRCODE_SEQNUMSET;

  if ((ssp->statusflag & SEQSET_BASQUAL) && 
      (sqp->qualp == NULL || sqp->type != SEQTYP_FASTQ))
    return ERRCODE_SEQTYP;

  if (ssp->statusflag&SEQSET_COMPRESSED) {
    errcode = appendMangledToCompressedSeq(ssp->sqp, sqp->datap, 
					  (ssp->statusflag & SEQSET_TERMCHAR)>0);
  } else {
    errcode = appendSeqSegment(ssp->sqp, NULL, NULL, sqp->datap, 0, 0, 0, 
			       (ssp->statusflag & SEQSET_TERMCHAR)>0,NULL);
  }
  if ((errcode))
    return errcode;

  if ((ssp->statusflag & SEQSET_BASQUAL) &&
      (errcode = appendSeqSegment(ssp->qqp, NULL, NULL, sqp->qualp, 0, 0, 0, 
				  (ssp->statusflag & SEQSET_TERMCHAR)>0,NULL)))
    return errcode;

  if (ssp->n_seq+1 >= ssp->n_alloc &&
      (errcode = reallocSeqSet(ssp, ssp->n_seq+2)))
    return errcode;
  
  ssp->sop[ssp->n_seq+1] =  ssp->sop[ssp->n_seq] + sqp->datap->size;
  if (ssp->statusflag & SEQSET_TERMCHAR) {
    ssp->sop[ssp->n_seq+1]++;
  }

  namsiz = (ssp->n_seq > 0)? ssp->namoffs[ssp->n_seq]: 0;
  if (sqp->headp->size + namsiz >= ssp->nam_alloc && 
      (errcode = reallocSeqSetName(ssp, sqp->headp->size + namsiz + 1)))
    return errcode;

  cp = ssp->namebasep + namsiz;
  strncpy(cp, sqp->headp->basep, sqp->headp->size);
  cp[sqp->headp->size] = '\0';
  ssp->namoffs[ssp->n_seq + 1] = namsiz + sqp->headp->size + 1;
  ssp->n_seq++;
  return ERRCODE_SUCCESS;
}

int seqSetAddFromFastqFile(ErrMsg *errmsg, SeqSet *ssp, 
			   SeqFastq *sqbufp,
			   const SeqCodec *codecp, 
			   const char *filnam, 
			   char verbose)
{
  int errcode = ERRCODE_SUCCESS;
  SEQNUM_t sctr;
  SeqIO *sfp;

  sfp = seqIOopen(&errcode, filnam, SEQIO_READ, 0);
  if (errcode) {
    seqIOclose(sfp);
    ERRMSGNO(errmsg, errcode);
  }
  if (verbose) 
    printf("File %s opened for reading sequences in FASTA/FASTQ format ...\n",
	   filnam);

  sctr = 0; 
  while(!seqIOstatus(sfp)) {
    if (verbose) 
      printf("Reading sequence %llu ...\n", (unsigned long long) sctr++);
    if ((errcode = seqFastqRead(sqbufp, sfp)))
      ERRMSGNO(errmsg, errcode);
    if ((errcode = seqFastqEncode(sqbufp, codecp)))
      ERRMSGNO(errmsg, errcode);
    if ((errcode = seqSetAddSequence(ssp, sqbufp)))
      ERRMSGNO(errmsg, errcode);
  }
  if (seqIOstatus(sfp) && 
      seqIOstatus(sfp) != ERRCODE_EOF) 
    ERRMSGNO(errmsg, seqIOstatus(sfp));
  seqIOclose(sfp);

  if ((ssp->statusflag & SEQSET_BASQUAL) && 
      (ssp->sqp->size != ssp->qqp->size))
     return ERRCODE_ASSERT;

  if (verbose)
    printf("%llu sequences read from file %s\n", (unsigned long long) sctr, filnam);

  return ERRCODE_SUCCESS;
}

int seqSetCompress(SeqSet *ssp, const SeqCodec *codecp)
{
  int errcode = ERRCODE_SUCCESS;
  if (!(ssp->statusflag & SEQSET_COMPRESSED) && 
      ssp->sqp->code == SEQCOD_ASCII &&
      (errcode = encodeSeq(ssp->sqp, codecp)))
    return errcode;
  errcode = compressSeq(ssp->sqp);
  if (!errcode) ssp->statusflag |= SEQSET_COMPRESSED;
  return errcode;
}

int seqSetWriteBinFil(const SeqSet *ssp, const char *filnam)
{
  int errcode;
  FILE *fp;
  uint32_t headval[SEQSET_HEADLEN];
  uint32_t seqsiz;
  uint32_t *seqlenp;
  SEQNUM_t s;
  uint64_t bigval, totsiz, seqnamsiz;

  if (ssp->n_seq < 1)
    return ERRCODE_FAILURE;

  if (!(ssp->statusflag&SEQSET_COMPRESSED))
    return ERRCODE_SEQCODE;

  /* write header (format version and number of sequences */
  bigval = (uint64_t) ssp->n_seq;
  headval[0] = (uint32_t) (bigval&MASK_32BIT);   /* number of sequences lower 32 bits*/
  headval[1] = (uint32_t) (bigval>>NUM_32BIT);   /* number of sequences remaining bits */
  bigval = (uint64_t) ssp->namoffs[ssp->n_seq];  /* size of name space */
  headval[2] = (uint32_t) (bigval & MASK_32BIT); /* lower 32 bits */
  headval[3] = (uint32_t) (bigval>>NUM_32BIT);   /* upper (most significant) 32 bits */
  bigval = (uint64_t) ssp->sqp->size;            /* total sequence length */
  headval[4] = (uint32_t) (bigval & MASK_32BIT);
  headval[5] = (uint32_t) (bigval>>NUM_32BIT);
  headval[6] = (uint32_t) ssp->statusflag;
  headval[7] = (uint32_t) 0;

  /* calculate space occupied by name offsets as multiple of 
   * 32-bit integer */
  seqnamsiz = (ssp->namoffs[ssp->n_seq]-1)/NBYTES_PER_INT32+1;
  seqsiz = ssp->sqp->size/MAXN_PER_UNIT+1;
  totsiz = SEQSET_HEADLEN + seqsiz + ssp->n_seq + seqnamsiz;
  if (totsiz > UINT32_MAX)
    return ERRCODE_FILSIZ;

  ECALLOCP(ssp->n_seq, seqlenp);
  if (seqlenp == NULL)
    return ERRCODE_NOMEM;

  fp = filioOpenForWriting(&errcode, (uint32_t) totsiz, FILIOTYP_SEQSET,
			   SEQSET_FORMAT_VERSION, SEQSET_HEADLEN,
			   headval, filnam, SEQSET_FILNAMEXT);
  if (errcode) {
    EFCLOSE(fp);
    return errcode;
  }

  fwrite(ssp->namebasep, sizeof(char), ssp->namoffs[ssp->n_seq], fp);
  for (s=0; s<ssp->n_seq; s++) {
    seqlenp[s] = ssp->sop[s+1] - ssp->sop[s];
  }
  fwrite(seqlenp, sizeof(uint32_t), ssp->n_seq, fp);
  free(seqlenp);
  seqlenp = NULL;

  /* write concatenated sequences */
  errcode = writeCompressedSeq(fp, ssp->sqp);
  if (!errcode && ferror(fp))
    errcode = ERRCODE_WRITEERR;

  /* write base qualities */
  if (!errcode && (ssp->statusflag & SEQSET_BASQUAL)) {
    if (ssp->sqp->size != ssp->qqp->size) 
      errcode = ERRCODE_ASSERT;
    else 
      fwrite(ssp->qqp->basep, sizeof(char), ssp->qqp->size+1, fp);
  }

  return EFCLOSE(fp);
}

SeqSet *seqSetReadBinFil(int *errcode, const char *filnam)
{
  UCHAR_t is_endianid, typ;
  uint32_t header[SEQSET_HEADLEN];
  uint32_t totsiz, version;
  uint32_t headsiz = SEQSET_HEADLEN;
  SEQNUM_t i, seqnum_dat;
  
  uint32_t *seqlenp = NULL;
  uint64_t seqsiz, namsiz, seqnum;
  SETSIZ_t j;
  uint32_t statusflg;
  SeqSet *ssp = 0;
  FILE *fp = filioOpenForReading(errcode, &is_endianid, &totsiz,
				 &typ, &version,
				 &headsiz, header, filnam, SEQSET_FILNAMEXT);
  if (*errcode)
    return 0;

  if (typ != FILIOTYP_SEQSET ||
      (version > SEQSET_FORMAT_VERSION || version+3 < SEQSET_FORMAT_VERSION) ||
      headsiz > SEQSET_HEADLEN || headsiz < 4) {
    *errcode = ERRCODE_FILTYP;
    if (version + 2 < SEQSET_FORMAT_VERSION)
      *errcode = ERRCODE_FILVERSION;
    else if (headsiz != SEQSET_HEADLEN)
      *errcode = ERRCODE_FILEFORM;
    fclose(fp);
    return 0;
  }
  
  if (version == SEQSET_FORMAT_VERSION) {
    seqnum = (((uint64_t) header[1])<<NUM_32BIT) + header[0];
    namsiz = (((uint64_t) header[3])<<NUM_32BIT) + header[2];
    seqsiz = (((uint64_t)header[5])<<NUM_32BIT) + header[4];
    statusflg = header[6];
    if (seqnum > SEQNUM_MAX)
      *errcode = ERRCODE_SEQNUMSET;
    seqnum_dat = seqnum;
  } else if (version == 3) {
    seqnum = (uint64_t) header[0];
    namsiz = (((uint64_t) header[2])<<NUM_32BIT) + header[1];
    seqsiz = (((uint64_t)header[4])<<NUM_32BIT) + header[3];
    statusflg = header[5];
    if (seqnum > SEQNUM_MAX)
      *errcode = ERRCODE_SEQNUMSET;
    seqnum_dat = seqnum;
  } else {
    seqnum = (uint64_t) header[0];
    namsiz = header[1];
    seqsiz = header[2];
    statusflg = header[3];
    if (seqnum + 1 > SEQNUM_MAX)
      *errcode = ERRCODE_SEQNUMSET;
    seqnum_dat = seqnum + 1;
  }
  if ((*errcode)) {
    fclose(fp);
    return 0;
  }
  if (seqnum < 1 || seqnum + 1 > SEQNUM_MAX ||
      seqsiz > SETSIZ_MAX ||
      namsiz > SETSIZ_MAX ||
      statusflg > UINT8_MAX) {
    *errcode = ERRCODE_FILEFORM;
    fclose(fp);
    return 0;
  }

  if (seqnum + 1 > SEQSET_BLOCKSIZ_DEFAULT) {
    ssp = seqSetCreate(0, (UCHAR_t) statusflg);
    if (ssp)
      *errcode = reallocSeqSet(ssp, seqnum+1);
  } else {
    ssp = seqSetCreate(seqnum+1, (UCHAR_t) statusflg);
  }
  if (NULL == ssp)
    *errcode = ERRCODE_NOMEM;
  if (!(*errcode))
    *errcode = reallocSeqSetName(ssp, namsiz);
  
  if (!(*errcode) &&
      ECALLOCP(seqnum_dat, seqlenp) == NULL)
    *errcode =  ERRCODE_NOMEM;

  
  if (!(*errcode) &&
      fread(ssp->namebasep, sizeof(char), namsiz, fp) != namsiz)
    *errcode = ERRCODE_FILEFORM;

  if ((*errcode)) {
    seqSetDelete(ssp);
    fclose(fp);
    return 0;
  }
    
  /* set sequence name pointers */
  ssp->namoffs[0] = 0;
  for (i=1, j=0; j<namsiz && i<= (SEQNUM_t) seqnum; j++)
    if (!(ssp->namebasep[j]))
      ssp->namoffs[i++] = j+1;

  if (i-1 != (SEQNUM_t) seqnum || 
      ssp->namoffs[seqnum] != namsiz) {
    *errcode = ERRCODE_ASSERT;
  } else {
    if (fread(seqlenp, sizeof(uint32_t), 
	      seqnum_dat, fp) != (size_t) seqnum_dat) {
      *errcode = ERRCODE_FILEFORM;
    } else {
      *errcode = readCompressedSeqOfKnownLength(ssp->sqp, fp, seqsiz);
      if (!(*errcode) && ferror(fp))
	*errcode = ERRCODE_READERR;
      else if (!(is_endianid)) {
	SETSIZ_t nunit_seqsiz = seqsiz/MAXN_PER_UNIT + 1;
	if (nunit_seqsiz <= UINT32_MAX) {
	  filioSwapEndian(seqlenp, seqnum_dat);
	  filioSwapEndian((uint32_t *) ssp->sqp->basep, (uint32_t) nunit_seqsiz);
	} else {
	  *errcode = ERRCODE_FILEFORM;
	}
      }
    }
    if (!(*errcode)) {
      if (version == SEQSET_FORMAT_VERSION || version == 3) {
	ssp->sop[0] = 0;
	for (i=0; i<seqnum_dat; i++)
	  ssp->sop[i+1] = ssp->sop[i] + seqlenp[i];
      } else {
	for (i=0; i<seqnum_dat; i++)
	  ssp->sop[i] = seqlenp[i];
      }
      if ((ssp->statusflag & SEQSET_BASQUAL)) {
	if (reallocSeqBlocks(ssp->qqp, seqsiz+1))
	  *errcode = ERRCODE_NOMEM;
	else {
	  if (fread(ssp->qqp->basep, sizeof(char), seqsiz+1, fp) != seqsiz+1 ||
	      ferror(fp)) {
	    *errcode = ERRCODE_READERR;
	  } else {
	    ssp->qqp->size = (SETSIZ_t) seqsiz;
	    if (ssp->sqp->size != ssp->qqp->size)
	      *errcode = ERRCODE_FILEFORM;
	  }
	}
	if (!(*errcode) && !(is_endianid))
	  filioSwapEndian((uint32_t *) ssp->qqp->basep, seqsiz + 1);
      }
    }
  }

  EFCLOSE(fp);
  free(seqlenp);

  ssp->n_seq = (SEQNUM_t) seqnum;
  if (statusflg & SEQSET_COMPRESSED)
    ssp->statusflag |= SEQSET_COMPRESSED;

  return ssp;
}

SEQNUM_t seqSetGetOffsets(const SeqSet *ssp, const SETSIZ_t **soffs)
{
  if (soffs) *soffs = ssp->sop;
  return ssp->n_seq;
}
    
int seqSetFetchSegment(SeqFastq *sqp,
		       SETSIZ_t *offs_start, SETSIZ_t *offs_end, 
		       const SeqSet *ssp, const SeqCodec *codecp)
{
  int errcode;
  UCHAR_t ntc = 0;
  BOOL_t with_qual = (ssp->statusflag & SEQSET_BASQUAL) && (sqp->type == SEQTYP_FASTQ);
  SEQLEN_t len, os, oe, tcpos[MAXN_TERMCHAR];
  
  if (*offs_start >= ssp->sop[ssp->n_seq] || *offs_start >= *offs_end ||
      *offs_end >= ssp->sqp->size)
    return ERRCODE_ARGRANGE;
  len = *offs_end - *offs_start + 1;
  blankSeq(sqp->datap);
  if (with_qual) 
    blankSeq(sqp->qualp);

  if ((ssp->statusflag&SEQSET_COMPRESSED)) {
    errcode = uncompressSeq(sqp->datap, tcpos, &ntc, ssp->sqp, *offs_start, len, codecp);
  } else {
    errcode = appendSeqSegment(sqp->datap, tcpos, &ntc, ssp->sqp, *offs_start, len, 0, 0, NULL);
  }
  if ((errcode))
    return errcode;

  if (ntc > 0) {
    if (ntc == 1) {
      if (tcpos[0] > len/2) {
	os = 0;
	oe = tcpos[0]-1;
      } else {
	os = tcpos[0] + 1;
	oe = len - 1;
      }
    } else {
      os = tcpos[0]+1;
      oe = tcpos[1]-1;
    }
    errcode = cropSeq(sqp->datap, os, oe);

    if (!errcode && ntc > 3)
      errcode = ERRCODE_SEGPARTFETCH;
    *offs_end = *offs_start + oe;
    *offs_start += os;
  }
  if (!(errcode) && (with_qual)) {
    len = *offs_end - *offs_start + 1;
    errcode = appendSeqSegment(sqp->qualp, tcpos, &ntc, ssp->qqp, *offs_start, len, 0, 0, NULL);
  }

  return errcode;
}

int seqSetFetchSegmentBySequence(SeqFastq *sqp, SEQNUM_t seqidx,
				 SEQLEN_t offs, SEQLEN_t len, 
				 const SeqSet *ssp, const SeqCodec *codecp)
{
  int errcode;
  BOOL_t with_qual = (ssp->statusflag & SEQSET_BASQUAL) && (sqp->type == SEQTYP_FASTQ);
  SEQLEN_t slen;
  
  if (seqidx >= ssp->n_seq)
    return ERRCODE_ARGRANGE;

  slen = ssp->sop[seqidx+1] - ssp->sop[seqidx];
  if (slen>0 && (ssp->statusflag & SEQSET_TERMCHAR)) slen--;
  if (offs >= slen) 
    return ERRCODE_SEQOFFS;
  if (len < 1) len = slen;
  if (offs + len > slen) {
    len = slen - offs;
  } 
  blankSeq(sqp->datap);
  if ((with_qual))
    blankSeq(sqp->qualp);

  if ((ssp->statusflag & SEQSET_COMPRESSED))
    errcode = uncompressSeq(sqp->datap, NULL, NULL, ssp->sqp, 
			    ssp->sop[seqidx]+offs, len, codecp);
  else
    errcode = appendSeqSegment(sqp->datap, NULL, NULL, ssp->sqp, 
			       ssp->sop[seqidx]+offs, len, 0, 0, NULL);

  if (!(errcode) && (with_qual)) 
    errcode = appendSeqSegment(sqp->qualp, NULL, NULL, ssp->qqp, 
			       ssp->sop[seqidx]+offs, len, 0, 0, NULL);

  return errcode;
}


int seqSetGetIndexAndOffset(SEQNUM_t *seqidx, SEQLEN_t *seqoffs, SETSIZ_t offs, const SeqSet *ssp)
{
  SEQNUM_t s, a, b;

  if (offs >= ssp->sop[ssp->n_seq])
    return ERRCODE_ARGRANGE;

  if (ssp->n_seq < 1) 
    return ERRCODE_FAILURE;

  a = 0;
  b = ssp->n_seq;
  while (a<b-1) {
    s = (b+a)/2;
    if (ssp->sop[s] > offs) {
      b = s;
    } else {
      a = s;
    }
  }
  if (seqidx) *seqidx = a;
  if (seqoffs) *seqoffs = ssp->sop[a];

  return ERRCODE_SUCCESS;
}

SEQLEN_t seqSetGetSeqDatByIndex(SETSIZ_t *offs, const char **name, SEQNUM_t seqidx, const SeqSet *ssp)
{
  SEQLEN_t slen;
  if (seqidx >= ssp->n_seq) {
    if (name) *name = 0;
    return 0;
  }
  if (name) *name = ssp->namebasep + ssp->namoffs[seqidx];
  if (offs) *offs = ssp->sop[seqidx];
  slen = ssp->sop[seqidx+1] - ssp->sop[seqidx];

  return ((ssp->statusflag&SEQSET_TERMCHAR) && slen > 0)? slen-1: slen;
}

SEQNUM_t seqSetGetSeqNumAndTotLen(SETSIZ_t *totseqlen, const SeqSet *ssp)
{
  if (totseqlen) {
    *totseqlen = ssp->sop[ssp->n_seq];
    if (ssp->statusflag & SEQSET_TERMCHAR) {
      if (*totseqlen > (SETSIZ_t) ssp->n_seq)
	*totseqlen -= ssp->n_seq;
      else 
	*totseqlen = 0;
    }
  }

  return ssp->n_seq;
}
