/** Handling of compressed alignment strings (pairwise alignments) */

/*****************************************************************************
 *****************************************************************************
 *                                                                           *
 *  Copyright (C) 2010 - 2014 Genome Research Ltd.                           *
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
#include <ctype.h>
#include <limits.h>

#include "elib.h"
#include "diffstr.h"

/****************************************************************************
 ****************************** Constants ***********************************
 ****************************************************************************/

enum {
  CIGAR_CLIPPED_HARD     = 'H', /**< Hard clipped */
  CIGAR_CLIPPED_SOFT     = 'S', /**< Soft clipped */
  NUMBUF_MAXLEN          = 7,   /**< Maximum buffer length for contsal */
  DEFAULT_BLOCKSIZ = 256,       /**< block size for memory allocation */
  DIFFSTR_MARGIN = 5,
  DIFFBLOCKS_BLKSZ = 64,     /**< Block size for DiffBlocks memory 
				  * allocation */
/**< maximum number of characters to print contents of one byte of a
 * compressed alignment string in any of the output formats
 * (incl. spaces).  This is used to calculate an upper limit for the length
 * of an output string */
  NCHAR_MAX_PRINTLEN = 5,  
  NCHAR_MAX_UNITLEN = 12, /**< Maximum number of characters to print one CIGAR unit 
			   * (i.e. 'I 3', 'M 1028') */
  NCHAR_MAX_CLIPLEN = 13, /**< max string length for clipping (31 bit int)
			   * ceiling(31*ln2/ln10)+3 */
};

enum CIGARSTR_TYPES {
  CIGARSTRTYP_SILENTMISMATCH = 0x01, /* don't mark mismatches */
  CIGARSTRTYP_EXTENDED = 0x02,       /* use extended format */
  CIGARSTRTYP_SOFTCLIPPED = 0x04,    /* use soft clipping */
};

static const char CIGAR_FORM[] = "%c %d ";
static const char CIGAR_EXTF[] = "%d%c"; /**< Extended CIGAR string format */
static const char DIFFSTR_SYMBOLS[] = "MDIS"; /**< Diff string symbols, must be
					       * coordinated with DIFFSTR_CODES.
					       * This is used for raw output */
static const char DIFFSTR_SYMBOLS_X[] = "MDIX"; /**< Diff string symbols, must be
						 * coordinated with DIFFSTR_CODES */
static const char DIFFSTR_NULLSTR[] = "*";

/****************************************************************************
 ******************************* Typedefs ***********************************
 ****************************************************************************/

typedef unsigned char UCHAR;
typedef unsigned char CIGTYP_t; /**< combination of CIGARSTR_TYPES */
typedef signed short SWSCOR;


/****************************************************************************
 **************************** Private types *********************************
 ****************************************************************************/
#ifdef diffstr_debug
typedef struct _DIFFBUFF {
  char *strp;
  int len;
  int n_alloc;
} DIFFBUFF;
#endif

typedef struct _BLOCK { /** A pair of segments that alingn w/o indels */
  int unprof_start;   /**< start of the segement of the unprofiled sequence */
  int prof_start;     /**< start of the segement of the profiled sequence */
  int len;            /**< length of the segments */
} BLOCK;

/****************************************************************************
 ***************************** Opaque types *********************************
 ****************************************************************************/

struct _DiffBlocks {
  BLOCK *blkp;  /**< Array of blocks of aligning w/o indels */
  int nblk;     /**< Number of blocks in array */
  int n_alloc;  /**< Number of blocks for which memory is allocated */
  int ablksz;   /**< Chunk size for memory allocation (a the number of elements) */
};

struct _DiffView { /**< Explicit alignment string representation */
  char *strp;
  size_t strlen;
  size_t n_alloc;
  int blksz;
};

/****************************************************************************
 ***************************** Public types *********************************
 ****************************************************************************/

/****************************************************************************
 ******************************** Macros ************************************
 ****************************************************************************/

#define SETDIFF(gap, typ) (UCHAR) ((gap) + (((UCHAR) (typ)) << DIFFSTR_TYPSHIFT))
/**< Set the one-character code (typ:gap) for the alignment segment */

#define DIFFSTR_GET_TYP(uc, typ) (typ) = (unsigned char) ((uc) >> DIFFSTR_TYPSHIFT);
/******************************************************************************
 ********************* Private Methods of Type DiffView ***********************
 ******************************************************************************/

static int reallocView(DiffView *dvp, size_t newsiz)
{
  int errcode = ERRCODE_SUCCESS;

  newsiz = (newsiz + dvp->blksz - 1)/dvp->blksz;
  newsiz *= dvp->blksz;
  
  if (newsiz != dvp->n_alloc) {
    void *hp = EREALLOCP(dvp->strp, newsiz);
    if (NULL == hp) {
      errcode = ERRCODE_NOMEM;
    } else {
      dvp->strp = hp;
      dvp->n_alloc = newsiz;
    }
  }

  return errcode;
}

/******************************************************************************
 ********************* Private Methods of Type CIGAR_WRITER *******************
 ******************************************************************************/

typedef int (CIGAR_WRITER) (void *, CIGTYP_t, char, int);

inline static int writeCigarToFile(void *top, CIGTYP_t cigtyp, char typc, int ctr)
{
  return (ctr > 0)? (cigtyp & CIGARSTRTYP_EXTENDED)? 
    fprintf((FILE *)top, CIGAR_EXTF, ctr, typc):
    fprintf((FILE *)top, CIGAR_FORM, typc, ctr): 
    fprintf((FILE *)top, "%c", typc);
}

inline static int writeCigarToStr(void *top, CIGTYP_t cigtyp, char typc, int ctr)
{
  return (ctr > 0)? (cigtyp & CIGARSTRTYP_EXTENDED)?  
    sprintf((char *)top, CIGAR_EXTF, ctr, typc):
    sprintf((char *)top, CIGAR_FORM, typc, ctr):
    sprintf((char *)top, "%c", typc);
}

#ifdef DIFFSTR_SUPERFLUOUS
static int writeCigarToDiffView(void *top, CIGTYP_t cigtyp, char typc, int ctr)
{
  int nchar = 0;
  DiffView *dvp = (DiffView *) top;

  if (dvp->strlen + NCHAR_MAX_UNITLEN <= dvp->n_alloc ||
      !reallocView(dvp, dvp->strlen + NCHAR_MAX_UNITLEN))
    nchar = writeCigarToStr(dvp->strp, cigtyp, typc, ctr);
  if (nchar > 0) dvp->strlen += nchar;

  return nchar;
}
#endif
/******************************************************************************
 ********************* Private Methods of Type DIFFBUFF ***********************
 ******************************************************************************/
#ifdef diffstr_debug
int fillDIFFBUFF(DIFFBUFF *p, const DIFFSTR_T *dstrp)
{
  int i, j;
  DIFFSTR_T count, typ;

  diffStrPrintf(stdout, dstrp, DIFFSTRFORM_RAW, 0, 0, 0);
  p->len = 0;
  for (i=0; dstrp[i] != 0 && i<INT_MAX; i++) {
    DIFFSTR_GET(dstrp[i], count, typ);
    for (j=0; j<count; j++) {
      if (p->len >= p->n_alloc) 
	return ERRCODE_OVERFLOW;
      p->strp[p->len++] = DIFFSTR_SYMBOLS[DIFFCOD_M];
    }
    if (p->len >= p->n_alloc) 
      return ERRCODE_OVERFLOW;
    p->strp[p->len++] = DIFFSTR_SYMBOLS[typ];
  }
  if (p->len > 0)
    p->len--; /* the last code should be S */
  if (p->len >= p->n_alloc) 
    return ERRCODE_OVERFLOW;
  p->strp[p->len] = '\0';
  return (typ == DIFFCOD_S)? ERRCODE_SUCCESS: ERRCODE_DIFFSTR;
}

void deleteDIFFBUFF(DIFFBUFF *p)
{
  if (p) {
    free(p->strp);
  }
  free(p);
}

DIFFBUFF *createDIFFBUFF(const DIFFSTR_T *dstrp)
{
  DIFFBUFF *p;
  int alilen = diffStrCalcAliLen(NULL, dstrp);

  if (alilen < 1) return NULL;
  EMALLOCP0(p);
  if (!p) return NULL;

  ECALLOCP(alilen+1, p->strp);
  p->n_alloc = alilen+1;
  if (p->strp == NULL ||
      fillDIFFBUFF(p, dstrp)) {
    deleteDIFFBUFF(p);
    p = NULL;
  }

  return p;
}
#endif
/******************************************************************************
 ******************** Private Methods of Type DIFFSTR_T ***********************
 ******************************************************************************/
static void fprintfDiffStrRaw(FILE *fp, const DIFFSTR_T *dstrp)
{
  const UCHAR *ucp;
  fprintf(fp, "(");
  if ((dstrp)) {
    for (ucp=dstrp; (*ucp); ucp++)
      fprintf(fp, "%c:%2.2hi|", DIFFSTR_SYMBOLS[(*ucp) >> DIFFSTR_TYPSHIFT],
	      (short) ((*ucp) & DIFFSTR_COUNTMASK));
  }
  fprintf(fp, "M:00)\n");
}

static int sprintfDiffStrRaw(char *sp, const DIFFSTR_T *dstrp)
{
  const UCHAR *ucp;
  int nchar = 0;
  nchar += sprintf(sp, "(");
  if ((dstrp)) {
    for (ucp=dstrp; (*ucp); ucp++)
      nchar += sprintf(sp + nchar, "%c:%2.2hi|", DIFFSTR_SYMBOLS[(*ucp) >> DIFFSTR_TYPSHIFT],
		       (short) ((*ucp) & DIFFSTR_COUNTMASK));
  }
  nchar += sprintf(sp + nchar, "M:00)");
  return nchar;
}

static void fprintfDiffStrPlain(FILE *fp, const DIFFSTR_T *dstrp)
{
  const UCHAR *ucp;
  if ((dstrp)) {
    for (ucp=dstrp; (*ucp); ucp++)
      fprintf(fp, "%c%i", DIFFSTR_SYMBOLS[(*ucp) >> DIFFSTR_TYPSHIFT],
	      (short) ((*ucp) & DIFFSTR_COUNTMASK));
  }
}

static int sprintfDiffStrPlain(char *sp, const DIFFSTR_T *dstrp)
{
  const UCHAR *ucp;
  int nchar = 0;
  if ((dstrp)) {
    for (ucp=dstrp; (*ucp); ucp++)
      nchar += sprintf(sp + nchar, "%c%i", DIFFSTR_SYMBOLS[(*ucp) >> DIFFSTR_TYPSHIFT],
		       (short) ((*ucp) & DIFFSTR_COUNTMASK));
  }
  return nchar;
}

static int writeDiffStrCIGAR(void * const top, int *nchar,
			     unsigned char cgt, /* combination of CIGARSTR_TYPES */
			     const DIFFSTR_T *diffstr,
			     int clip_start, 
			     int clip_end,
			     CIGAR_WRITER *writerp)
{
  unsigned int count, prev_count = 0;
  unsigned short typ = DIFFCOD_M, prev_typ = DIFFCOD_M;
  int nc = 0;
  const DIFFSTR_T *ucp;
  const char clipchar = (char) (( 0 != (cgt & CIGARSTRTYP_SOFTCLIPPED) )? 
				CIGAR_CLIPPED_SOFT: CIGAR_CLIPPED_HARD);

#define WRITE_CIGAR(chr, ctr) nc += (*writerp)(writeCigarToStr == *writerp? ((char *)top)+nc: top, \
					       cgt, (chr), (ctr))
  if (!diffstr) {
    WRITE_CIGAR(DIFFSTR_NULLSTR[0], 0);
    *nchar = nc;
    return ERRCODE_SUCCESS;
  }

  if (!(*diffstr)) return ERRCODE_FAILURE;

  if (clip_start > 0 && (cgt & CIGARSTRTYP_EXTENDED))
    WRITE_CIGAR(clipchar, clip_start);
    
  for (ucp = diffstr; *ucp; ucp++) {
    DIFFSTR_GET(*ucp, count, typ);
    
    if (prev_typ == DIFFCOD_M) {
      prev_count += count;
      if (typ == DIFFCOD_M || 
	  (typ == DIFFCOD_S && (cgt & CIGARSTRTYP_SILENTMISMATCH))) {
	prev_count++;
	continue;
      } 
    } else if (typ == prev_typ && count < 1) {
      prev_count++;
      continue;
    }

    if (prev_count > 0)
      WRITE_CIGAR(DIFFSTR_SYMBOLS_X[prev_typ], prev_count);
 
    if (typ == DIFFCOD_M || (typ == DIFFCOD_S && (cgt & CIGARSTRTYP_SILENTMISMATCH))) {
      prev_count = count + 1;
      prev_typ = DIFFCOD_M;
    } else {
      if (count > 0 && prev_typ != DIFFCOD_M)
	 WRITE_CIGAR(DIFFSTR_SYMBOLS_X[DIFFCOD_M], count);
      prev_count = 1;
      prev_typ = typ;
    }
  }

  if (typ != DIFFCOD_S)
    return ERRCODE_DIFFSTR;
  
  if (prev_count > 1) /* may end with mismatch */
    WRITE_CIGAR(DIFFSTR_SYMBOLS_X[(cgt & CIGARSTRTYP_SILENTMISMATCH)?
				  DIFFCOD_M:DIFFCOD_S], 
		prev_count-1);

  if (clip_end > 0 && (cgt & CIGARSTRTYP_EXTENDED))
    WRITE_CIGAR(clipchar, clip_end);

  *nchar += nc;
  return ERRCODE_SUCCESS;
}

static int scrollDiffStr(DIFFSTR_T *dfsp, int pos_uprof_target, UCHAR isEnd,
		  int *pos_uprof, int *pos_prof, int *dfs_offs)
{
  int i, ie;
  UCHAR count, typ;
  int count_uprof, count_prof;
  int count_uprof_lastmatch, count_prof_lastmatch;

  count_uprof = count_prof = 0;
  count_uprof_lastmatch = count_prof_lastmatch = 0;
  ie = 0;
  for (i=0; dfsp[i]; i++) {
    DIFFSTR_GET(dfsp[i], count, typ);
 
    if (typ == DIFFCOD_M || typ == DIFFCOD_S) {
      count_prof += count+1;  
      count_uprof += count+1;
    } else if (typ == DIFFCOD_I) {
      count_prof += count + 1;
      count_uprof += count;
    } else {/* typ == DIFFCOD_D */
      count_prof += count;
      count_uprof += count + 1;
    }
    if (typ == DIFFCOD_S && !dfsp[i+1]) {
      count_uprof--;
      count_prof--;
    }
    if (count > 0) {
      count_uprof_lastmatch = count_uprof;
      count_prof_lastmatch = count_prof;
      ie = i;
    }
    if (isEnd && count_uprof > pos_uprof_target)
      break;
    if (count_uprof_lastmatch > pos_uprof_target)
      break;
  }
  if (dfs_offs) *dfs_offs = ie;
  if (pos_prof) *pos_prof = count_prof_lastmatch;
  if (pos_uprof) *pos_uprof = count_uprof_lastmatch;
  return (count_uprof < pos_uprof_target)? ERRCODE_DIFFSTR: ERRCODE_SUCCESS;   
}

static int scrollDIFFSTRStartEnd 
(			  
 int *start_unprof, 
 int *end_unprof,
 int *start_prof,
 int *end_prof,
 DIFFSTR_T *count_start,
 DIFFSTR_T *count_end,
 DIFFSTR_T *typ_start,
 int *idx_start,
 int *idx_end,
 int start_unprof_target,
 int end_unprof_target,
 const DIFFSTR_T *diffstrp
 )
     /**< Scroll alignment represented by the compressed alignment string to a target start position
      * on the unprofiled sequence. If this position is in a deletion or a mismatch, scroll
      * further up to the next perfect match.
      *
      * Then scroll alignment to a target end position on the
      * unprofiled sequence. If this position is in a deletion or a
      * mismatch, scroll back to the last position that was a perfect
      * match.
      * 
      * \return ERRCODE_NOMATCH if segment selected does not contain a match.
      *
      * \param start_unprof Returns actual start position in the unprofiled seuqence after the scrolling.
      * \param end_unprof Returns actual end position in the unprofiled seuqence after the scrolling.
      * \param start_prof Returns position in profiled sequence corresponding to start_unprof.
      * \param end_prof Returns position in profiled sequence corresponding to end_unprof.
      * \param count_start Returns the number of matches following position (start_unprof, start_prof) up
      *        to the state typ_end. I.e. the starting character of the resulting segment of the 
      *        compressed alignment string is (*typ_end:*count_end) for *typ_end != DIFFCOD_M and
      *        (DIFFCOD_M:*count_end-1) for *typ_end == DIFFCOD_M.
      * \param count_end Returns the number of matches preceding (and including) position 
      *       (pos_unprof, pos_prof) up to the state typ_end. I.e. the last character of the resulting 
      *        segment of the compressed alignment string is (DIFFCOD_S:*count_end).
      * \param typ_start Returns the DIFF_CODE of the compressed diffstr at (start_unprof, start_prof).
      * \param idx_start Returns the index in the alignment string corresponding to (start_unprof, start_prof).
      * \param idx_end Returns the index in the alignment string corresponding to (end_unprof, end_prof).
      * \param start_unprof_target Target start position of the segment in the unprofiled sequence.
      * \param start_unprof_target Target end position of the segment in the unprofiled sequence.
      * \param diffstrp Compressed alignment string terminated by 0.
      */
{
  int i, idx_last, shift = 0, shift_last = 0, pos = 0, pos_last;
  DIFFSTR_T count = 0, count_add = 0, typ = 0;
#ifdef diffstr_debug
  printf("diffstr.c::scrollDIFFSTRStartEnd:");
  fprintfDiffStrRaw(stdout, diffstrp);
  printf("diffstr.c::scrollDIFFSTRStartEnd: start_unprof_target = %i\n", start_unprof_target);
#endif
  
  /* shift = *start_prof - *start_unprof */
  /* scroll to start */
  for (i=0; diffstrp[i] != 0 && i<INT_MAX; i++) {
    DIFFSTR_GET(diffstrp[i], count, typ);
    shift_last = shift;
    if (typ == DIFFCOD_M) {
      count++;
      count_add = 0; /* unprofiled sequence */
    } else if (typ == DIFFCOD_S) {
      count_add = 1;
    } else if (typ == DIFFCOD_I) {
      shift++; 
      count_add = 0;
    } else { /* typ == DIFFCOD_D */
      count_add = 1;
      shift--;
    }
    /* count: number of exact matches
     * count + count_add: number of bases moved along unprofiled sequence
     * shift: pos_profiled - pos_unprofiled
     */
    pos += count;
    if (pos > start_unprof_target && count > 0) 
      break;
    pos += count_add;
  }
  if (i >= INT_MAX)
    return ERRCODE_DIFFSTR;
  
  if (diffstrp[i] == 0)
    return ERRCODE_FAILURE;
  /* pos: position on unprofiled sequence - count_add, 
   * i.e. excluding last mismatch or deletion */
  idx_last = i;
  *count_start = (DIFFSTR_T) (pos - start_unprof_target);
  if (*count_start > count)
    *count_start = count;  
  *start_unprof = pos - *count_start;
  *start_prof = *start_unprof + shift_last;

  pos_last = pos;
  pos += count_add;
  *idx_start = i;
  *typ_start = typ;

#ifdef diffstr_debug
  printf("diffstr.c::scrollDIFFSTRStartEnd: start_unprof = %i, start_prof = %i\n", 
	 *start_unprof, *start_prof);
  printf("diffstr.c::scrollDIFFSTRStartEnd: count_start = %i, pos_last = %i\n", 
	 *count_start, pos_last);
  printf("diffstr.c::scrollDIFFSTRStartEnd: idx_start = %i, typ_start = %c\n", 
	 *idx_start, DIFFSTR_SYMBOLS[*typ_start]);
  printf("diffstr.c::scrollDIFFSTRStartEnd: end_unprof_target = %i\n", end_unprof_target);
  printf("diffstr.c::scrollDIFFSTRStartEnd:");
  fprintfDiffStrRaw(stdout, diffstrp + i + 1);
#endif

  if (*start_unprof > end_unprof_target) {
    /* no match in sequence */
    return ERRCODE_NOMATCH; 
  }
  if (pos <= end_unprof_target) {
    /* scroll to end */
    //shift_last = shift;
    for (i++; diffstrp[i] != 0 && i<INT_MAX; i++) {
      DIFFSTR_GET(diffstrp[i], count, typ);
      if (count > 0)
	shift_last = shift;
      
      if (typ == DIFFCOD_M) {
	count++;
	count_add = 0;
      }	else if (typ == DIFFCOD_S) {
	count_add = 1;
      } else if (typ == DIFFCOD_I) {
	count_add = 0;
	shift++;
      } else { /* typ == DIFFCOD_D */
	count_add = 1;
	shift--;
      }
      
      pos += count;
      /* count is the number of exact matches */

      if (count > 0) {
	pos_last = pos;
	idx_last = i;
      }
      pos += count_add;

      if (pos > end_unprof_target)
	break;

    }
    if (diffstrp[i] == 0)
      i--;
    else if (i >= INT_MAX) 
      return ERRCODE_OVERFLOW;
  }  /* pos <= end_unprof_target */
  if (pos_last > end_unprof_target) {
    *count_end = (DIFFSTR_T) (pos_last - end_unprof_target - 1);
    if (*count_end > count)
      return ERRCODE_ASSERT;
    *count_end = (DIFFSTR_T) (count - *count_end);
    *end_unprof = end_unprof_target;
    *idx_end = i;
  } else {
     DIFFSTR_GET(diffstrp[idx_last], count, typ);
     if (typ == DIFFCOD_M)
       count++;
     *count_end = count;
     *end_unprof = pos_last - 1;
     *idx_end = idx_last;
  }

  *end_prof = *end_unprof + shift_last;
#ifdef diffstr_debug
  printf("diffstr.c::scrollDIFFSTRStartEnd: end_unprof = %i, end_prof = %i\n", 
	 *end_unprof, *end_prof);
  printf("diffstr.c::scrollDIFFSTRStartEnd: count_end = %i, idx_end = %i\n", 
	 *count_end, *idx_end);
#endif

  return ERRCODE_SUCCESS;
}


/******************************************************************************
 ******************* Private Methods of Type DiffBlocks ***********************
 ******************************************************************************/

static int reallocDiffBlocks(DiffBlocks *p, int n_blks)
{
  size_t newsz = (n_blks + p->ablksz - 1)/p->ablksz;
  void *hp;
  newsz *= p->ablksz;
  if (newsz > INT_MAX)
    return ERRCODE_OVERFLOW;
  hp = EREALLOCP(p->blkp, newsz);
  if (hp == NULL)
    return ERRCODE_NOMEM;
  p->blkp = (BLOCK *) hp;
  p->n_alloc = (int) newsz;

  return ERRCODE_SUCCESS;
}

/******************************************************************************
 ******************** Public Methods of Type DiffBlocks ***********************
 ******************************************************************************/

DiffBlocks *diffBlocksCreate(int blksz)
{
  DiffBlocks *p;
  EMALLOCP0(p);
  if (p != NULL) {
    p->nblk = 0;
    if (blksz < 1) blksz = DIFFBLOCKS_BLKSZ;
    ECALLOCP(blksz, p->blkp);
    if (NULL == p->blkp) {
      diffBlocksDelete(p);
      p = NULL;
    } else {
      p->n_alloc = p->ablksz = blksz;
    }
  }
  return p;
}

void diffBlocksDelete(DiffBlocks *p)
{
  if (p)
    free(p->blkp);
  free(p);
}

int diffBlocksGetNumber(const DiffBlocks *p)
{
  return p->nblk;
}

int diffBlocksGetLen(int *unprof_start, int *prof_start,
		     int blkno, const DiffBlocks *p)
{
  int rv = 0;
  if (blkno < p->nblk && blkno >= 0) {
    BLOCK *bp = p->blkp + blkno;
    if (unprof_start) *unprof_start = bp->unprof_start;
    if (prof_start) *prof_start = bp->prof_start;
    rv = bp->len;
  } else {
    if (unprof_start) *unprof_start = 0;
    if (prof_start) *prof_start = 0;
  }
  return rv;
}

int diffStrFindBlocks(DiffBlocks *dbp, const DIFFSTR_T *diffstrp)
{
  int u = 0, p = 0, l = 0;
  DIFFSTR_T count = 0, typ = DIFFCOD_M;
#define ADDBLOCK(dbp, u, p, l) if (l > 0) {			\
    BLOCK *bp;							\
    if ((dbp)->nblk >= (dbp)->n_alloc) {			\
      int errcode = reallocDiffBlocks((dbp), (dbp)->nblk+1);	\
      if (errcode) return errcode;				\
    }								\
    bp = (dbp)->blkp + dbp->nblk++;				\
    bp->unprof_start = (u);					\
    bp->prof_start = (p);					\
    bp->len = (l);						\
    (u) += (l);							\
    (p) += (l);							\
    (l) = 0;							\
}
    
  dbp->nblk = 0;
  if (NULL == diffstrp)
    return ERRCODE_SUCCESS;

  for (; (*diffstrp); diffstrp++) {
    DIFFSTR_GET(diffstrp[0], count, typ);
    l += count;
    if (typ == DIFFCOD_I) {
      ADDBLOCK(dbp, u, p, l);
      p++;
    } else if (typ == DIFFCOD_D) {
      ADDBLOCK(dbp, u, p, l);
      u++;
    } else {
      l++;
    }
  }
  if (typ != DIFFCOD_S)
    return ERRCODE_DIFFSTR;
  l--;
  ADDBLOCK(dbp, u, p, l);

  return ERRCODE_SUCCESS;
}

/******************************************************************************
 ********************* Private Methods of Type DiffStr ************************
 ******************************************************************************/

/******************************************************************************
 ********************* Public Methods of Type DiffStr *************************
 ******************************************************************************/

DiffStr *diffStrCreate(int blocksiz)
{
  DiffStr *dfsp;

  EMALLOCP0(dfsp);
  if (dfsp && diffStrInit(dfsp, blocksiz)) {
    diffStrDelete(dfsp);
    dfsp = 0;
  }

  return dfsp;
}

void diffStrDelete(DiffStr *dfsp)
{
  diffStrCleanUp(dfsp);
  free(dfsp);
}

int diffStrInit(DiffStr *p, int blocksiz)
{
  p->len = p->n_alloc = p->blksz = 0;
  p->dstrp = NULL;
  if (blocksiz < 1) blocksiz = DEFAULT_BLOCKSIZ;
  ECALLOCP(blocksiz, p->dstrp);

  if (!(p->dstrp))
    return ERRCODE_NOMEM;

  p->n_alloc = p->blksz = blocksiz;

  return ERRCODE_SUCCESS;
}

void diffStrCleanUp(DiffStr *p)
{
  if (p) {
    free(p->dstrp);
    p->dstrp = NULL;
    p->len = p->n_alloc = p->blksz = 0;
  }
}

int diffStrRealloc(DiffStr *dfsp, int n_new)
{
  int errcode;
  size_t nsiz;
  DIFFSTR_T *hp;

#ifdef diffstr_debug
  printf("diffStrRealloc(%i, %i)\n", dfsp->n_alloc, n_new);
#endif
  if (!dfsp->dstrp)
    return ERRCODE_NULLPTR;
  if (!dfsp->dstrp &&
      (errcode = diffStrInit(dfsp, 0)))
    return errcode;

  if (!n_new) nsiz = dfsp->n_alloc + dfsp->blksz;
  else nsiz = ((size_t) ((n_new + dfsp->blksz - 1)/dfsp->blksz))*dfsp->blksz;

  if (nsiz > INT_MAX) 
    return ERRCODE_OVERFLOW;

#ifdef diffstr_debug
  printf("diffStrRealloc nsiz = %i\n", (int) nsiz);
#endif

  hp = EREALLOCP(dfsp->dstrp, nsiz);
  if (!hp) return ERRCODE_NOMEM;

  dfsp->dstrp = hp;
  dfsp->n_alloc = (int) nsiz;
  if (dfsp->len > dfsp->n_alloc)
    dfsp->len = dfsp->n_alloc;

  return ERRCODE_SUCCESS;
}

int diffStrCopy(DiffStr *dfsp, const DIFFSTR_T *diffstrp)
{
  int errcode; 
  short l;
  const  DIFFSTR_T *dstrp = diffstrp;
  DIFFSTR_T *ucp;

  for (l=0; (diffstrp[l]); l++)
    if (l >= SHRT_MAX) 
      return ERRCODE_OVERFLOW;
  
  if (l>=dfsp->n_alloc &&
      (errcode = diffStrRealloc(dfsp, l)))
    return errcode;
  
  dfsp->len = l+1; /* length included termination */
  
  for (ucp = dfsp->dstrp; l>=0; l--)
    *ucp++ = *dstrp++;
  
  return ERRCODE_SUCCESS;
}

int diffStrAdd(DiffStr *top, const DIFFSTR_T *fcp, int len)
{
  int l, errcode = ERRCODE_SUCCESS;
  size_t newlen = top->len + ((len > 0)? len: 0);
  DIFFSTR_T *tcp;

  if (len < 1 || NULL == fcp)
    return ERRCODE_SUCCESS;

  if (newlen > INT_MAX)
    return ERRCODE_OVERFLOW;

  if (((int) newlen) > top->n_alloc &&
      (errcode = diffStrRealloc(top, (int) newlen)))
    return errcode;

  tcp = top->dstrp + top->len;
  for (l=0; l<len; l++) 
    *tcp++ = *fcp++;

  top->len = (int) newlen;
  if (newlen > 0 && top->dstrp[newlen-1] != 0)
    errcode = ERRCODE_ASSERT;

  return errcode; 
}

int diffStrAppend(DiffStr *top, const DiffStr *fromp)
{ 
  return diffStrAdd(top, fromp->dstrp, fromp->len);
}

int diffStrReverse(DiffStr *dfsp, const DIFFSTR_T *diffstrp)
{
  int errcode;
  UCHAR typ, count, count_prev;
  short l,u=0;
  DIFFSTR_T *ucp;
#ifdef results_debug
  printf("results_debug::reverseDIFFSTR(): string before ");
  fprintfDiffStrRaw(stdout, diffstrp);
#endif  
  for (l=0; diffstrp[l]; l++)
    if (l >= SHRT_MAX) 
      return ERRCODE_OVERFLOW;
  
  if (l>=dfsp->n_alloc &&
      (errcode = diffStrRealloc(dfsp, l+1)))
    return errcode;
  
  l--;
  DIFFSTR_GET(diffstrp[l], count_prev, typ);
  if (typ != DIFFCOD_S)
    return ERRCODE_DIFFSTR;

  ucp=dfsp->dstrp;
  for (l--;l>=0;l--) {
    DIFFSTR_GET(diffstrp[l], count, typ);
    if (typ == DIFFCOD_M) {
      count_prev = (UCHAR) (count_prev + count + 1);
      if (count_prev > DIFFSTR_MAXMISMATCH) {
	ucp[u++] = SETDIFF(DIFFSTR_MAXMISMATCH, DIFFCOD_M);
	count_prev -= DIFFSTR_MAXMISMATCH + 1;
      }
    } else {
      ucp[u++] = SETDIFF(count_prev, typ);
      count_prev = count; 
    }   
  }
  ucp[u++] = SETDIFF(count_prev, DIFFCOD_S);
  ucp[u++] = SETDIFF(0, DIFFCOD_M);
  dfsp->len = u;
#ifdef results_debug
  printf("results_debug::reverseDIFFSTR(): string after ");
  fprintfDiffStrRaw(stdout, dfsp->dstrp);
#endif  

  return ERRCODE_SUCCESS;
}

int diffStrCalcSeqLen(int *len_prof, int *len_unprof, const DIFFSTR_T *diffstrp)
{
  int errcode = ERRCODE_SUCCESS;
  DIFFSTR_T typ = DIFFCOD_M, count;
  int pl=0, ul=0;
  const int maxlen = INT_MAX - DIFFSTR_MAXMISMATCH;

  for (; (*diffstrp) && pl < maxlen && ul < maxlen; diffstrp++) {
    DIFFSTR_GET(diffstrp[0], count, typ);
    if (typ == DIFFCOD_I) {
      ul += count;
      pl += count + 1;
    } else if (typ == DIFFCOD_D) {
      ul += count + 1;
      pl += count;
    } else {
      ul += count + 1;
      pl += count + 1;
    }
  }
  if (!(*diffstrp)) {
    if (typ == DIFFCOD_S) {
      ul--;
      pl--;
    }
  } else {
    errcode = ERRCODE_DIFFSTR;
  }
  if (len_prof) *len_prof = pl;
  if (len_unprof) *len_unprof = ul;

  return errcode;
}

int diffStrCalcAliLen(int *matchnum, const DIFFSTR_T *diffstrp)
{
  DIFFSTR_T typ=DIFFCOD_M, count;
  int l, m;
  const int maxl = INT_MAX - DIFFSTR_MAXMISMATCH;

  if (matchnum) *matchnum = 0;
  l = m = 0;
  for (; *diffstrp && l < maxl; diffstrp++) {
    DIFFSTR_GET(diffstrp[0], count, typ);
    if (typ == DIFFCOD_M)
      m += count + 1;
    else 
      m += count;
    l += count+1;
  }
  if (matchnum) *matchnum = m;
  if (typ == DIFFCOD_S && !(*diffstrp))
    l--;
  return ((*diffstrp))? 0: l;
}

int diffStrGetDiffStats(int *n_sub, int *n_ins, int *n_del, 
			const DIFFSTR_T * diffstrp)
{
  int ni = 0, nd = 0, ns = 0;
  int errcode = ERRCODE_SUCCESS;
  DIFFSTR_T typ=DIFFCOD_M;

  if (n_sub) *n_sub = 0;
  if (n_ins) *n_ins = 0;
  if (n_del) *n_del = 0;

  for (; !(errcode) && (*diffstrp); diffstrp++) {
    DIFFSTR_GET_TYP(diffstrp[0], typ);
    if (typ == DIFFCOD_I) {
      if (INT_MAX == ni)
	errcode = ERRCODE_OVERFLOW;
      else
	ni++;
    } else if (typ == DIFFCOD_D) {
      if (INT_MAX == nd)
	errcode = ERRCODE_OVERFLOW;
      else
	nd++;
    } else if (typ == DIFFCOD_S) {
      if (diffstrp+1 != NULL) {
	if (INT_MAX == ns)
	  errcode = ERRCODE_OVERFLOW;
	else
	  ns++;
      }
    }
  }
  if (!(errcode) && ((diffstrp) || typ != DIFFCOD_S))
    errcode = ERRCODE_DIFFSTR;
  if (n_sub) *n_sub = ns;
  if (n_ins) *n_ins = ni;
  if (n_del) *n_del = nd;
 
  return errcode;
}

int diffStrGenerateFromMismatches(int *dlen, DIFFSTR_T *diffstrp, 
				  const int mmpos[], int mmnum, int qlen)
     
{
  int i, j, n, ntot;
  int supos;
  UCHAR *dcp = NULL;

  if (dlen) *dlen = 0;
  if (diffstrp) dcp = diffstrp;

  if (mmnum < 1) {
    n = (qlen-1)/DIFFSTR_MAXMISMATCH; /* number of bytes - 1 for exact matching diffstr (excl. termination) */
    supos = qlen;
  } else {
    n = (int) (mmpos[0]>0)? (mmpos[0]-1)/DIFFSTR_MAXMISMATCH: 0;
    supos = mmpos[0];
  }
  if (dcp) {
    for (j=0;j<n;j++)
      *dcp++ = SETDIFF(DIFFSTR_MAXMISMATCH, DIFFCOD_M);
    *dcp++ = SETDIFF(supos-n*DIFFSTR_MAXMISMATCH+1, DIFFCOD_S);
  }
  ntot = n+1;
  if (mmnum>0) {
    for (i=1; i<mmnum; i++) {
      if (mmpos[i] <= mmpos[i-1]) 
	return ERRCODE_ASSERT;
      n = (int) (mmpos[i]-mmpos[i-1]-1)/DIFFSTR_MAXMISMATCH;
      if (dcp) {
	for (j=0; j<n; j++)
	  *dcp++ = SETDIFF(DIFFSTR_MAXMISMATCH, DIFFCOD_M);
	*dcp++ = SETDIFF(mmpos[i]-mmpos[i-1]-n*DIFFSTR_MAXMISMATCH, DIFFCOD_S);
      }
      ntot += n+1;
    }
    if (mmpos[i-1] != qlen-1) {
      n = (int) (qlen-mmpos[i-1]-1)/DIFFSTR_MAXMISMATCH;
      if (dcp) {
	for (j=0; j<n; j++)
	  *dcp++ = SETDIFF(DIFFSTR_MAXMISMATCH, DIFFCOD_M);
	*dcp++ = SETDIFF(qlen-mmpos[i-1]-n*DIFFSTR_MAXMISMATCH, DIFFCOD_S);
      }
      ntot += n+1;
    }
  }
  if ((dcp)) 
    *dcp++ = SETDIFF(0, DIFFCOD_M); /* termination */
  if (dlen) *dlen = ntot+1;

  return ERRCODE_SUCCESS;
}

int diffStrPrintf(FILE *fp, const DIFFSTR_T *diffstrp, char outform, 
		  int clip_start, int clip_end, char is_softclipped)
{
  int errcode = ERRCODE_SUCCESS;
  CIGTYP_t cgt = (CIGTYP_t)((is_softclipped)? CIGARSTRTYP_SOFTCLIPPED: 0);
  int nchar = 0;
  switch (outform) {
  case DIFFSTRFORM_RAW:
    fprintfDiffStrRaw(fp, diffstrp);
    break;
  case DIFFSTRFORM_PLAIN:
    fprintfDiffStrPlain(fp, diffstrp);
    break;
  case DIFFSTRFORM_CIGNORM:
    errcode = writeDiffStrCIGAR(fp, &nchar, 
				CIGARSTRTYP_SILENTMISMATCH, 
				diffstrp, 0, 0, writeCigarToFile);
    break;
  case DIFFSTRFORM_CIGEXT:
    errcode = writeDiffStrCIGAR(fp, &nchar, 
				(CIGTYP_t)(cgt|CIGARSTRTYP_EXTENDED|CIGARSTRTYP_SILENTMISMATCH), 
				diffstrp, clip_start, clip_end, writeCigarToFile);
    break;
  case DIFFSTRFORM_CIGEXT_XMISMATCH:
    errcode = writeDiffStrCIGAR(fp, &nchar, 
				(CIGTYP_t)(cgt|CIGARSTRTYP_EXTENDED), 
				diffstrp, clip_start, clip_end, writeCigarToFile);
   break;
  default:
    errcode = ERRCODE_ASSERT;
    break;
  }

  return errcode;
}

int diffStrPrintfStr(char *sp, int *nchar, const DIFFSTR_T *diffstrp, char outform, 
		     int clip_start, int clip_end, char is_softclipped)
{
  int errcode = ERRCODE_SUCCESS;
  CIGTYP_t cgt = (CIGTYP_t) ((is_softclipped)? CIGARSTRTYP_SOFTCLIPPED: 0);

  *nchar = 0;
  switch (outform) {
  case DIFFSTRFORM_RAW:
    *nchar = sprintfDiffStrRaw(sp, diffstrp);
    break;
  case DIFFSTRFORM_PLAIN:
    *nchar = sprintfDiffStrPlain(sp, diffstrp);
    break;
  case DIFFSTRFORM_CIGNORM:
    errcode = writeDiffStrCIGAR(sp, nchar, 
				0, 
				diffstrp, 0, 0, writeCigarToStr);
    break;  
  case DIFFSTRFORM_CIGEXT:
    errcode = writeDiffStrCIGAR(sp, nchar, 
				(CIGTYP_t) (cgt|CIGARSTRTYP_EXTENDED|CIGARSTRTYP_SILENTMISMATCH), 
				diffstrp, clip_start, clip_end, 
				writeCigarToStr);
    break;
  case DIFFSTRFORM_CIGEXT_XMISMATCH:
    errcode = writeDiffStrCIGAR(sp, nchar, 
				(CIGTYP_t)(cgt|CIGARSTRTYP_EXTENDED), 
				diffstrp, clip_start, clip_end, 
				writeCigarToStr);
    break;
  default:
    errcode = ERRCODE_ASSERT;
    break;
  }

  return errcode;
}

#define DIFFSTR_CHECKMEM()    if (dfsp->len + 1 >= dfsp->n_alloc && \
	(errcode = diffStrRealloc(dfsp, dfsp->len+1))) return errcode

int diffStrParsePlain(DiffStr *dfsp, const char *rawstrp)
{
  int errcode, c;
  short i;
  DIFFSTR_T code;
  int count;
  char numbuf[NUMBUF_MAXLEN];
  const char *cp;

  dfsp->len = 0;

  for (cp = rawstrp; (*cp) && isspace(*cp); cp++);
  while (*cp) {
    for (i=0; (cp[i]) && i<NUMBUF_MAXLEN-1 && isdigit(cp[i]); i++)
      numbuf[i] = cp[i];
    if (i >= NUMBUF_MAXLEN-1)
      return ERRCODE_DIFFSTR;
    
    numbuf[i] = '\0';
    count = atoi(numbuf);
    while (count > DIFFSTR_MAXMISMATCH) {
      DIFFSTR_CHECKMEM();
      dfsp->dstrp[dfsp->len++] = (DIFFCOD_M<<DIFFSTR_TYPSHIFT) + (DIFFSTR_T) (DIFFSTR_MAXMISMATCH & DIFFSTR_COUNTMASK);
      count -= DIFFSTR_MAXMISMATCH;
    }
    
    c = toupper(cp[i]);
    for (code=0; (DIFFSTR_SYMBOLS[code]) && c !=  toupper(DIFFSTR_SYMBOLS[code]); code++);
    if (!DIFFSTR_SYMBOLS[code])
      break;

    DIFFSTR_CHECKMEM();

    dfsp->dstrp[dfsp->len++] = (DIFFSTR_T) ((code<<DIFFSTR_TYPSHIFT) + (DIFFSTR_T) (count & DIFFSTR_COUNTMASK));
    
    cp += i+1;
  }
  DIFFSTR_CHECKMEM();
  dfsp->dstrp[dfsp->len++] = 0;
  
  return ERRCODE_SUCCESS;
}

int diffStrParseSimul(DiffStr *dfsp, 
		      int *clip_start, int *clip_end, 
		      unsigned char isSAMCIGAR, const char *rawstrp)
{
  int errcode = ERRCODE_SUCCESS;
  int i, count, curr_count = 0;
  const char *cp;
  char numbuf[NUMBUF_MAXLEN];
  DIFFSTR_T code;

  *clip_start = 0;
  *clip_end = 0;

  /* skip white space at beginning */
  for (cp = rawstrp; (*cp) && isspace(*cp); cp++);
  
  dfsp->len = 0;

  while (*cp && ERRCODE_SUCCESS == errcode) {
    unsigned char isClip = 0;
    /* load digits into buffer */
    for (i=0; (cp[i]) && i<NUMBUF_MAXLEN-1 && isdigit(cp[i]); i++)
      numbuf[i] = cp[i];
    if (i >= NUMBUF_MAXLEN-1 || !cp[i])
      return ERRCODE_DIFFSTR;
    numbuf[i] = '\0';
    count = atoi(numbuf);
    
    cp += i;
    if (!(*cp))
      return ERRCODE_DIFFSTR;
    /* determine var type */
    if (isSAMCIGAR) {
      for (code=0; (DIFFSTR_SYMBOLS_X[code]) && *cp !=  DIFFSTR_SYMBOLS_X[code]; code++);
      if (!DIFFSTR_SYMBOLS_X[code]) {
	if (*cp == CIGAR_CLIPPED_HARD || *cp == CIGAR_CLIPPED_SOFT) {
	  if (dfsp->len == 0 || *clip_start == 0) {
	    *clip_start = count;
	    isClip = 1;
	  } else if (*clip_end == 0) {
	    *clip_end = count;
	    isClip = 1;
	  } else {
	    errcode = ERRCODE_DIFFSTR;
	  }
	} else {
	  errcode = ERRCODE_DIFFSTR;
	}
      }
    } else {
      for (code=0; (DIFFSTR_SYMBOLS[code]) && *cp !=  DIFFSTR_SYMBOLS[code]; code++);
      if (!DIFFSTR_SYMBOLS[code])
	errcode = ERRCODE_DIFFSTR;
    }

    if (errcode)
      break;
    
    if ((isClip)) {
      if(*clip_end > 0 && (*(cp+1)))
	errcode = ERRCODE_DIFFSTR;
    } else if (code == DIFFCOD_M) {
      curr_count += count;
      while (curr_count > DIFFSTR_MAXMISMATCH+1) {
	DIFFSTR_CHECKMEM();
	dfsp->dstrp[dfsp->len++] = (((DIFFSTR_T) DIFFCOD_M)<<DIFFSTR_TYPSHIFT) +
	  (DIFFSTR_T) DIFFSTR_MAXMISMATCH;
	curr_count -= DIFFSTR_MAXMISMATCH + 1;
      }
    } else {
      while (0 < count--) {
	DIFFSTR_CHECKMEM();
	dfsp->dstrp[dfsp->len++] = (DIFFSTR_T) ((code<<DIFFSTR_TYPSHIFT) +
						(curr_count & DIFFSTR_COUNTMASK));
	curr_count = 0;
      }
    }
    cp++;
  }
  if (curr_count > 0) {
    DIFFSTR_CHECKMEM();
    dfsp->dstrp[dfsp->len++] = (DIFFSTR_T)
      ((((DIFFSTR_T) DIFFCOD_S)<<DIFFSTR_TYPSHIFT) + 
       (curr_count & DIFFSTR_COUNTMASK));
  }
  DIFFSTR_CHECKMEM();
  dfsp->dstrp[dfsp->len++] = 0;
  
  return errcode;
}
 
int diffStrCrop(DIFFSTR_T *diffstrp, int *dstrlen,
		int start_unprof_target, int end_unprof_target,
		int *start_unprof, int *end_unprof, 
		int *start_prof, int *end_prof)
{
  int errcode;
  UCHAR count, typ;
  UCHAR count_over;
  int j, is, ie, dd;
  
  *dstrlen = 0;
  *start_unprof = *end_unprof = 0;
  *start_prof = *end_prof = 0;
#ifdef diffstr_debug
  printf("diffStrCrop: before split");
  fprintfDiffStrRaw(stdout, diffstrp);
#endif

  if ((errcode = scrollDiffStr(diffstrp, start_unprof_target, 
			       0, start_unprof, start_prof, &is)))
    return errcode;
  
  if ((errcode = scrollDiffStr(diffstrp, end_unprof_target, 
			       1, end_unprof, end_prof, &ie)))
    return errcode;
  
  DIFFSTR_GET(diffstrp[is], count, typ);
  if (typ == DIFFCOD_M) count++;
  if (*start_unprof < start_unprof_target)
    return ERRCODE_ASSERT;
  
  if (*start_unprof < start_unprof_target + count) {
    count = (UCHAR) (*start_unprof - start_unprof_target);
  }

  *start_unprof -= count;
  *start_prof -= count;

  j = 0;  
  if (is < ie) {
    if (count < 1 || (typ == DIFFCOD_M && count < 2)) {
      count_over = count;
    } else {
      if (typ == DIFFCOD_M) 
	count--;
      diffstrp[j++] = SETDIFF(count, typ); 
      count_over = 0;
    }
    for (is++; is<ie && diffstrp[is]; is++) {
      DIFFSTR_GET(diffstrp[is], count, typ);
      count = (UCHAR) (count + count_over);
      if (count > DIFFSTR_MAXMISMATCH) {
	diffstrp[j++] = SETDIFF(DIFFSTR_MAXMISMATCH, DIFFCOD_M);
	count -= DIFFSTR_MAXMISMATCH + 1;
	count_over = 0;
      }
      if (count < 1 &&  typ == DIFFCOD_M) {
	count_over = 1;
      } else {
	diffstrp[j++] = SETDIFF(count, typ);
      }
    }
  
    DIFFSTR_GET(diffstrp[ie], count, typ);
    count = (UCHAR) (count + count_over);
    if (typ == DIFFCOD_M) count++;
  }
  if (*end_unprof < end_unprof_target + 1) {
    if (end_unprof_target > *end_unprof + DIFFSTR_MAXMISMATCH)
      return ERRCODE_DIFFSTR;
    count = (UCHAR) (end_unprof_target + 1 - *end_unprof);
    /* e = pos_uprof; boundary changed */
  } else {
    if (*end_unprof > end_unprof_target + count + 1)
      return ERRCODE_DIFFSTR;
    dd = (int) (*end_unprof - end_unprof_target);  
    count = (UCHAR) (count + 1 - dd);
    *end_unprof -= dd;
    *end_prof -= dd;
  }
  diffstrp[j++] = SETDIFF(count, DIFFCOD_S);
  diffstrp[j] = SETDIFF(0, DIFFCOD_M);

#ifdef diffstr_debug
  printf("diffStrCrop: after split");
  fprintfDiffStrRaw(stdout, diffstrp);
#endif  
  
  *dstrlen = j+1; 
  return ERRCODE_SUCCESS;
} 

int diffStrFetchSegment(DiffStr *dfsp, const DIFFSTR_T *diffstrp,
			int start_unprof_target, int end_unprof_target,
			int *start_unprof, int *end_unprof, 
			int *start_prof, int *end_prof)
{
  int errcode;
  if ((errcode = diffStrCopy(dfsp, diffstrp)))
    return errcode;

  return diffStrCrop(dfsp->dstrp, &dfsp->len,
		     start_unprof_target, end_unprof_target,
		     start_unprof, end_unprof,
		     start_prof, end_prof);
}




int diffStrSegment(DiffStr *dfsp, const DIFFSTR_T *diffstrp,
		   int start_unprof_target, int end_unprof_target,
		   int *start_unprof, int *end_unprof, 
		   int *start_prof, int *end_prof)
{
  int errcode;
  int i, idx_start, idx_end;
  int maxsiz;
  int nmatch;
  DIFFSTR_T count, nmatch_start, nmatch_end, typ, typ_start;

#ifdef diffstr_debug
  DIFFBUFF *dstrbfp = createDIFFBUFF(diffstrp);
  if (!dstrbfp)
    return ERRCODE_DIFFSTR;
#endif

  dfsp->len = 0;
  nmatch = nmatch_start = nmatch_end = 0;				    
  errcode = scrollDIFFSTRStartEnd(start_unprof, end_unprof,
				  start_prof, end_prof,
				  &nmatch_start, &nmatch_end,
				  &typ_start, 
				  &idx_start, &idx_end,
				  start_unprof_target, 
				  end_unprof_target,
				  diffstrp);
  if (errcode)
    return errcode; /* hands down ERRCODE_NOMATCH if requested segment does
		     * not contain a match */
  
  maxsiz = idx_end - idx_start + DIFFSTR_MARGIN;
  if (maxsiz > dfsp->n_alloc &&
      (errcode = diffStrRealloc(dfsp, maxsiz)))
    return errcode;
  
  nmatch = 0;
  if (idx_start == idx_end) {
    DIFFSTR_GET(diffstrp[idx_start], count, typ);
    if (typ == DIFFCOD_M) count++;
    nmatch_end = (DIFFSTR_T) (nmatch_end + nmatch_start - count);
  } else {
    if (typ_start == DIFFCOD_M) {
      nmatch = nmatch_start;
    } else if (nmatch_start > 0) {
      dfsp->dstrp[dfsp->len++] = SETDIFF(nmatch_start, typ_start);
      nmatch = 0;
    }
    if (idx_end > idx_start+1) {
      for(i=idx_start+1; i<idx_end && diffstrp[i] != 0;i++) {
	DIFFSTR_GET(diffstrp[i], count, typ);
	if (nmatch + count > INT_MAX)
	  return ERRCODE_DIFFSTR;
	nmatch += count;
	if (typ == DIFFCOD_M) {
	  nmatch++;
	  continue;
	}
	for(; nmatch > DIFFSTR_MAXMISMATCH; nmatch -= DIFFSTR_MAXMISMATCH + 1) {
	  dfsp->dstrp[dfsp->len++] = SETDIFF(DIFFSTR_MAXMISMATCH, DIFFCOD_M);
	  if (dfsp->len >= dfsp->n_alloc)
	    return ERRCODE_ALLOCBOUNDARY;
	}
	dfsp->dstrp[dfsp->len++] = SETDIFF(nmatch, typ);
	if (dfsp->len >= dfsp->n_alloc)
	  return ERRCODE_ALLOCBOUNDARY;
	nmatch = 0;
      }
    }
  }
  nmatch += nmatch_end;
  for(; nmatch > DIFFSTR_MAXMISMATCH + 1; nmatch -= DIFFSTR_MAXMISMATCH + 1) {
    dfsp->dstrp[dfsp->len++] = SETDIFF(DIFFSTR_MAXMISMATCH, DIFFCOD_M);
    if (dfsp->len >= dfsp->n_alloc)
      return ERRCODE_ALLOCBOUNDARY;
  }
  if (dfsp->len + 2 > dfsp->n_alloc)
    return ERRCODE_ALLOCBOUNDARY;
  dfsp->dstrp[dfsp->len++] = SETDIFF(nmatch, DIFFCOD_S);
  dfsp->dstrp[dfsp->len++] = SETDIFF(0, DIFFCOD_M);
 
#ifdef diffstr_debug 
  deleteDIFFBUFF(dstrbfp);
#endif

  return ERRCODE_SUCCESS;
}


int diffStrScore(const DIFFSTR_T *diffstrp, int *swscor, SWSCOR match, SWSCOR mismatch, 
		 SWSCOR gapopen, SWSCOR gapextend)
{
  int errcode = ERRCODE_SUCCESS;
  UCHAR count, typ = DIFFCOD_M;
  UCHAR is_gap_open = 0;
  int i, maxi;

  if (match < 1 || mismatch > 0 || gapopen > gapextend || gapopen > 0) 
    return ERRCODE_ASSERT;

  maxi = INT_MAX - (DIFFSTR_MAXMISMATCH+1)*match;
  *swscor=0;

  for(i=0; (diffstrp[i]) && i < maxi; i++) {
    DIFFSTR_GET(diffstrp[i], count, typ);
    if (typ == DIFFCOD_S) {
      *swscor += match*count + mismatch; 
      is_gap_open = 0;
    } else if (typ == DIFFCOD_M) {
      *swscor += match*(count+1);
      is_gap_open =0 ;
    } else if ((is_gap_open) && count < 1) {
      *swscor += gapextend;
    } else {
      *swscor += match*count + gapopen;
      is_gap_open = 1;
    }
  }
  if ((diffstrp[i]) || typ != DIFFCOD_S)
    errcode = ERRCODE_DIFFSTR;
  else 
    *swscor -= mismatch;

  return errcode;
}

int diffStrGetLevenshteinDistance(const DIFFSTR_T *diffstrp)
{
  UCHAR typ = DIFFCOD_M;
  int i, ed = 0;

  for (i=0; (diffstrp[i]) && i < INT_MAX; i++) {
    DIFFSTR_GET_TYP(diffstrp[i], typ);
    if (typ != DIFFCOD_M)
      ed++;
  }
  if (ed > 0 && i < INT_MAX && typ == DIFFCOD_S)
    ed--; /* don't count terminating S */

  return ed;
}

/******************************************************************************
 ********************** Public Methods of Type DiffView ***********************
 ******************************************************************************/

DiffView *diffStrCreateView(int blksz)
{
  DiffView *p;
  EMALLOCP0(p);

  if (p != NULL) {
    if (blksz < 1) blksz = DEFAULT_BLOCKSIZ;
    ECALLOCP(blksz, p->strp);
    if (NULL == p->strp) {
      diffStrDeleteView(p);
      p = NULL;
    } else {
      p->n_alloc = blksz;
      p->blksz = blksz;
      p->strlen = 0;
      p->strp[0] = '\0';
    }
  }

  return p;
}

void diffStrDeleteView(DiffView *p)
{
  if (p != NULL) 
    free(p->strp);
  free(p);
}

const char *diffStrGetViewStr(const DiffView *p)
{
  return p->strp;
}

int diffStrAsView(DiffView *dvp, const DIFFSTR_T *dstrp,
		  char outform,
		  int clip_start, int clip_end, char is_softclipped)
{
  int errcode = ERRCODE_SUCCESS;
  int i, len = 0;
  size_t maxlen;
  for (i=0; i < INT_MAX && (dstrp[i]); i++);
  if (i >= INT_MAX)
    return ERRCODE_OVERFLOW;
  
  maxlen = ((size_t) i)*NCHAR_MAX_PRINTLEN + 2*NCHAR_MAX_CLIPLEN; 

  if (maxlen >= dvp->n_alloc &&
      (errcode = reallocView(dvp, maxlen+1)))
    return errcode;

  errcode = diffStrPrintfStr(dvp->strp, &len, dstrp, 
			     outform, clip_start, clip_end, is_softclipped);

  if (!(errcode) && ((size_t) len) > maxlen)
    errcode = ERRCODE_ASSERT;

  return errcode; 
}
