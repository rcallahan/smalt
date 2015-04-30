/**< Pairwise sequence alignment using dynamic programming */

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

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "elib.h"
#include "alignment.h"
#include "alibuffer_struct.h"
#if defined SCORE_SIMD && defined alignment_debug
#include "swsimd.h"
#endif

enum {
  FALSE = 0,
  TRUE = 1,
  ALIBKTRK_BLKSZ_DEFAULT = 32768,
  ALIBKTRK_MEMLIM_DEFAULT = 1024*1024,
  ALIRSLTSET_DEFAULT_BLKSZ = 16,
  ALIDIFFSTR_DEFAULT_BLKSZ = 64,
  ALILEN_MIN = 5,                      /**< Minimum alignment length */
};

enum ALITRACK_BITFLAGS{ /**< Encoding of moves in dynamic programmin matrix, where
			 * the profiled sequence is along the fastest running index
			 * i.e. rows of the matrix */
  ALIBKTRKCOD_COL      = 1, /** bit 0: Move along column, from above (in forward direction) */
  ALIBKTRKCOD_ROW      = 2, /** bit 1: Move along row, from left to right (in forward direction) */
  ALIBKTRKCOD_DIA      = 3, /** diagonal move = (ALIBKTRKCOD_COL | ALIBKTRKCOD_ROW) */
};

#ifdef alignment_debug
enum ALIMETA_CONST {
  NUMSTR = 3,             /**< Number of strings in AliMeta */
  ALISTRINGS_DEFAULT_BLKSZ = 256,
  ALIMETA_SPACER = '-',    /**< Spacing character in AliMeta strings */
  ALIMETA_MISMATCH = ' ',  /**< Mismatch character in AliMeta strings */
  ALIMETA_DELIM = '+',     /**< Delimiting character (for Smith-Waterman score of 0) */
};
#endif

static const double LN0P25 = -1.386294; /* natural logarithm of 1/4 */
/****************************************************************************
 ******************************** Types *************************************
 ****************************************************************************/

typedef unsigned char UCHAR;
typedef unsigned char BOOL;
typedef unsigned char ALIDIRECT_t;
typedef signed char ALIMATSCOR_t;

typedef struct _ALICPLX { /**< Stats for complexity scaling of Smith-Waterman score */
  int *countp;   /**< Array of size n_types, used as buffer for complexity scaling of 
		  * alignment score */
  short n_types; /**< Number of distinct nucleotide types (size of array countp). */
  double lambda;
} ALICPLX;

typedef struct _ALIBAND {  /**< delimiters of alignment band */
  int band_width;  /**< Width of alignment band */
  int l_edge_orig; /**< original left (lower) edge of band  (along profiled sequence) */
  int r_edge_orig; /**< original right(higher) edge of band (along profiled sequence) */
  int l_edge;      /**< Absolute left (lower) edge of alignment band along profiled sequence
		    * (adjusted for segment start points s_left and conting from 0) */
  int r_edge;      /**< Absolute right (higer) edge of alignment band along unprofiled sequence 
		    * (adjusted for segment start point s_left and counting from 0) */
  int s_left_orig; /**< original left (lower, starting) position of the subject segment 
		    * (counting starts from 0) */
  int s_left;      /**< starting position of the subject segment adjusted to the 
		    * alignment band (counting from 0) */
  int s_len;       /**< Length of the subject segment adjusted to the alignment band */
  int s_totlen;    /**< Original length of the subject segment */
  int q_left_orig; /**< original left (lower, starting) position of the query  segment 
		    * counting from 0 */
  int q_left;      /**< starting position of the query (profiled) sequence segment 
		    * adjusted to band width (counting from 0) */
  int q_len;       /**< Length of the query (profiled) sequence segment adjusted
		    * to band width */
  int q_totlen;    /**< original length of profiled squence segment */
} ALIBAND;

typedef struct _ALITRACK { /**< Directions for back tracking */
  int max_i;       /**< Position in the unprofiled sequence corresponding to the maximum. This
		    * translates to row index of the DPM[i,j] cell as i = max_i - ALIBAND.s_left */
  int max_j;       /**< Position in the profiled sequence corresponding to the maximum. This
		    * ranslates to the column (fastest running) index of the DPM[i,j] cell as
		    * j = max_j - ALIBAND.q_left_orig */
  ALIDPMSCOR_t max_scor; /**<  Maximum score */
  ALIDIRECT_t *bdp; /**< directions for back tracking bdp[k] is the element of the DPM [i,j] with
		     * k = (i-ALIBAND.s_left)*ALIBAND.bandwidth + (j-ALIBAND.l_edge) - (i-ALIBAND.s_left).
		     * Each element consists of a combination of ALIBACKTRACK_BITFLAGS 
		     * representing direction indicators (see annotation of ALIBACKTRACK_BITFLAGS) */
  size_t blksz;     /**< Block size (granularity) for memory (re-) allocation */
  size_t n_alloc;   /**< Currently allocated memory as the number of elements in bdp */
  size_t n_alloc_thresh; /**< threshold above which memory is (de-/re-)allocated in blocks */
} ALITRACK;

typedef struct _ALIMETA { /**< Contains intermediate alignment result. */
  int prof_start;    /**< Start (counting from 0) of the segment in profiled sequence */
  int prof_end;      /**< End (counting from 0) of the segment in profiled sequence */
  int nonprof_start; /**< Start (counting from 0) of the segment in second, unprofiled
		      * sequence */
  int nonprof_end;   /**< End (counting from 0) of the segment in second, unprofiled
		      * sequence */
  ALIDPMSCOR_t score;

  DiffStr dfs;          /**< short-hand alignment string *in reverse* along profiled sequence */
#ifdef alignment_debug
  char *string[NUMSTR]; /**< Pairwise alignment (in reverse) and consensus.
			 * string[0]: aligned profiled sequence (query). ALIMETA_SPACER for gap
			 * string[1]: aligned unprofiled sequence (reference). ALIMETA_SPACER for gap.
			 * string[2]: consensus ALIMETA_MISMATCH for mismatch, ALIMETA_SPACER for gap.*/
  int string_len;       /**< length of alignment strings */
  unsigned int allocsiz; /**< Memory allocated for alignment strings (as per string). */
  int blocksiz;         /**< Block size for memory allocation. */
  int string_end;       /**< End of alignment strings */
#endif
} ALIMETA;

typedef struct _ALIRESULT { /**< Contains result of pairwise alignment of two segments */
  int score;  /**< Alignment score (e.g. Smith-Waterman score */
  int qs;     /**< Start (counting from 0) of the segment in profiled sequence */
  int qe;     /**< End (counting from 0) of the segment in profiled sequence */
  int rs;     /**< Start (counting from 0) of the segment in second, unprofiled
		 * sequence */
  int re;     /**< End (counting from 0) of the segment in second, unprofiled
		 * sequence */
  DiffStr diffstr; /**< Short-hand alignment string (along profiled sequence) */
} ALIRESULT;

struct _AliRsltSet { /**< Set of aligned segments obtained for one 
		      * pairwise alignment task */
  ALIRESULT *rsp; /**< Array of results */
  short nres;     /**< Number of results (size of array) */
  short n_alloc;  /**< Number of array elements for which memory is allocated */
  short blksz;    /**< Block size for memory allocation */
  short dfblksz;  /**< Block size for compressed alignment strings */
  ALICPLX *cplxp; /**< Buffer used for complexity scaling of Smith-Waterman scores */
  ALIMETA meta;   /**< Buffer for intermediate results */
  ALITRACK track; /**< Buffer for back tracking */
};

/******************************************************************************
 ****************************** Private Methods *******************************
 ******************************************************************************/

/******************************************************************************
 ******************************* Public Methods *******************************
 ******************************************************************************/
int aliScoreDiffStr(int *swscor, const char *unprofiled_seqp, int unprofiled_seqlen, 
		    unsigned int profiled_offs,
		    const DIFFSTR_T *diffstrp, int diffstrlen, const ScoreProfile *scpp)
{
  int i;
  unsigned int profiled_len;
  int rs;
  signed char gap_init, gap_ext;
  UCHAR is_open = 0;
  DIFFSTR_T j, count, typ;
  signed char *const *scorepp = scoreGetProfile(NULL, &profiled_len, &gap_init, &gap_ext, scpp);
  
  *swscor = 0;
  rs = 0;
  for (i=0; i < diffstrlen && (diffstrp[i]); i++) {
    DIFFSTR_GET(diffstrp[i], count, typ);
    if (typ == DIFFCOD_M || (typ == DIFFCOD_S && (diffstrp[i+1])))
      count++;
    if (count > 0) {
      is_open = 0;
      for (j=0; j<count; j++) {
	*swscor += scorepp[(unprofiled_seqp[rs++]&SEQCOD_ALPHA_MASK)][profiled_offs++];
	if (profiled_offs > profiled_len || rs > unprofiled_seqlen)
	  return ERRCODE_ASSERT;
      }
    }
    if (typ == DIFFCOD_I || typ == DIFFCOD_D) {
      if ((is_open)) {
	*swscor -= gap_ext;
      } else {
	*swscor -= gap_init;
	is_open = 1;
      }
      if (typ == DIFFCOD_I) {
	profiled_offs++;
	if (profiled_offs > profiled_len)
	  return ERRCODE_ASSERT;
      } else {
	rs++;
	if (rs > unprofiled_seqlen)
	  return ERRCODE_ASSERT;
      }
    }
  }

  return (diffstrp[i])? ERRCODE_DIFFSTR: ERRCODE_SUCCESS;
}
/******************************************************************************
 ************************ Methods of Private Type ALICPLX *********************
 ******************************************************************************/

static void deleteALICPLX(ALICPLX *p)
{
  if (p) {
    free(p->countp);
  }
  free(p);
}

static ALICPLX *createALICPLX(const ScoreMatrix *smp)
{
  ALICPLX *p;
  short alphabetsiz;

  EMALLOCP0(p);
  if (p) {
    alphabetsiz = scoreMatrixGetAlphabetSize(smp);
    ECALLOCP(alphabetsiz, p->countp);
    if (!(p->countp)) {
      deleteALICPLX(p);
      p = 0;
    } else {
      p->n_types = alphabetsiz;
      p->lambda = scoreMatrixCalcLambda(smp);
    }
  }
  return p;
}

static void blankALICPLX(ALICPLX *p)
{
  short i;
  if (p) {
    for (i=0; i<p->n_types; i++)
      p->countp[i] = 0;
  }
  return;
}

static int scaleALICPLX(int *adj_score, int orig_score, 
			/* int pos_score, */
			const ALICPLX *cplxp)
{
  int i, t_counts, n_letters, count;
  double t_factor, t_sum;

  for (i = t_counts = t_factor = t_sum = n_letters = 0; i < cplxp->n_types; i++) {
    count = cplxp->countp[i];
    if ((count)) {
      t_factor += count * log((double) count);
      t_sum += count * LN0P25;
      t_counts += count;
      n_letters++;
    }
  }
  t_factor -=  t_counts * log((double)t_counts);
  t_sum -= t_factor;
  t_factor /=  t_counts * log((n_letters > 4) ? 1. / n_letters : .25);
  /* old_adj_score = orig_score + t_factor * pos_score - pos_score + .5; */
  *adj_score = orig_score + t_sum / cplxp->lambda + .999;
  
  if ((*adj_score) > orig_score) {
#ifdef alignment_debug
    fprintf(stderr, "\nScore increase: %d %d", orig_score, *adj_score);
    for (i = 0; i < cplxp->n_types; i++) {
      count =  cplxp->countp[i];
      if ((count))
	fprintf(stderr, "\n %c %d", (char)i, count);
    }
#endif
    return ERRCODE_CPLXSCOR;
  }
  if (*adj_score < 0)
    *adj_score = 0;

  return ERRCODE_SUCCESS;
}

/******************************************************************************
 ************************ Methods of Private Type ALIBAND *********************
 ******************************************************************************/
static int initALIBAND(ALIBAND *abp,
		int l_edge, int r_edge,
		int q_left, int q_right, int q_len,
		int s_left, int s_right, int s_len
		)
  /**< Set and adjust band parameters for banded Smith-Waterman alignment
   * Return ERRCODE_FAILURE if limits are inconsistent.
   *
   * Nomenclature:
   * 'query sequence'(q):   the profiled sequence. 
   * 'subject sequence' (s): the sequence to be aligned against the profile.
   * all positions/coordinates are 0-based.
   * 
   * \note Left and right edge of the alignment band are specified along the profiled sequence
   * (0 based). Position counting in sequence starts from 1.
   *
   * \param l_edge  ablsolute left (lower) edge of the alignment band along the profiled sequence 
   *                (0 based).
   * \param r_edge  absolute right (higher) edge of the alignment band along profiled sequence 
   *                (0 based).
   * \param q_left  start of the query (profiled) segment (0 based)
   * \param q_right end of the query (profiled) segment (0 based)
   * \param q_len   total query length (not profiled)
   * \param s_left  start of the subject segment (0 based)
   * \param s_right end of the subject segment (0 based)
   * \param s_len   total subject length (profiled)
   */
{
  /* Segment [s_left, s_right] on unprofiled sequence determines
   * segment [abp->q_left, abp->q_len-1] on the profiled squence
   * (origin (q=0, s=0)). The band [l_edge, r_edge] is specified along
   * profiled sequence with origin (q=0, s=0) for both sequences. 
   *
   * Relevant to the dynamic programming:
   * At the start: 
   * [abp->l_edge, abp->r_edge] are the band limits at the
   * start of the alignment (upper left corner). [abp->s_left,
   * abp->s_len-1] is the segment to be aligned on unprofiled
   * sequence.
   *
   * During the dynamic programming, abp->r_edge gets incremented
   * until full band width is reached, then abp->l_edge and abp->r_edge get
   * incremented (when stepping along unprofiled sequence) until abp->r_edge > abp->q_len-1
   * finally abp->l_edge gets incremented. End is abp->s_len on unprofiled sequence.
   *
   */
  abp->s_len = (s_right<0||s_right>=s_len) ? s_len : s_right + 1;
  abp->q_len = (q_right<0||q_right>=q_len) ? q_len : q_right + 1;

  abp->s_totlen = s_len;
  abp->q_totlen = q_len;

  abp->s_left = abp->s_left_orig = (s_left > 0 && s_left < abp->s_len) ? s_left : 0;
  abp->q_left = abp->q_left_orig = (q_left > 0 && q_left < abp->q_len) ? q_left : 0;

  abp->l_edge_orig = abp->l_edge = l_edge;
  abp->r_edge_orig = abp->r_edge = r_edge;
  abp->band_width = r_edge - l_edge + 1;

  if (abp->band_width <= 0) {
    abp->band_width = 0;
    abp->l_edge = abp->q_left;
    abp->r_edge = abp->q_len - 1;
  } else {
    /* if (abp->q_len - abp->l_edge_orig < abp->slen)  abp->s_len = abp->q_len - abp->l_edge_orig */
    if (abp->l_edge_orig + abp->s_len > abp->q_len) abp->s_len = abp->q_len - abp->l_edge_orig;
    abp->l_edge += abp->s_left;
    if(abp->l_edge >= abp->q_len || abp->r_edge_orig + abp->s_len <= abp->q_left) {
      /* checking that slen - (s_left[projected_on_s] = -abp->r_edge_orig) overlaps with segment (q_left, q_right)
       * at some stage */
      return ERRCODE_FAILURE;
    }
    abp->r_edge += abp->s_left;
    if (abp->r_edge < abp->q_left) {
      /* Band at beginning of unprofiled segment (abp->s_left) does not overlap with
       * profiled segment -> move the band down diagonal until right end overlaps with start 
       * of segment [q_start, q_end] on profiled sequence */
      abp->s_left += abp->q_left - abp->r_edge;
      abp->l_edge += abp->q_left - abp->r_edge;
      abp->r_edge = abp->q_left;
    }
    if (abp->r_edge > abp->q_len - 1) abp->r_edge = abp->q_len - 1;
  }
  abp->band_width = abp->r_edge - abp->l_edge + 1;

  return (abp->band_width >= 0)? ERRCODE_SUCCESS: ERRCODE_FAILURE;
}

#if defined alignment_debug || defined alignment_debug_limits
static void printAliBand(FILE *fp, const ALIBAND *bp)
{
  fprintf(fp, "== ALIBAND ==\n");
  fprintf(fp, "l_edge_orig = %i\n", bp->l_edge_orig);
  fprintf(fp, "l_edge      = %i\n", bp->l_edge);
  fprintf(fp, "r_edge_orig = %i\n", bp->r_edge_orig);
  fprintf(fp, "r_edge      = %i\n", bp->r_edge);
  fprintf(fp, "s_left_orig = %i\n", bp->s_left_orig);
  fprintf(fp, "s_left      = %i\n", bp->s_left);
  fprintf(fp, "q_left_orig = %i\n", bp->q_left_orig);
  fprintf(fp, "q_left      = %i\n", bp->q_left);
  fprintf(fp, "s_len       = %i\n", bp->s_len);
  fprintf(fp, "s_totlen    = %i\n", bp->s_totlen);
  fprintf(fp, "q_len       = %i\n", bp->q_len);
  fprintf(fp, "q_totlen    = %i\n", bp->q_totlen);
  fprintf(fp, "== End of ALIBAND ==\n");
}
#endif
#ifdef alignment_debug_limits
static void printEncSeq(FILE *fp, const char *seqp)
{
  const char const NTCOD[]="ACGT";
  const char *cp;
  for (cp = seqp; (*cp); cp++)
    fputc((int) NTCOD[(int) (*cp)&0x03], fp);
  fputc((int) '\0', fp);
}
#endif

/******************************************************************************
 ********************** Private Methods of Type ALITRACK **********************
 ******************************************************************************/
static void cleanupALITRACK(ALITRACK *p) 
{
  if (p) {
    free(p->bdp);
    p->bdp = 0;
  }
}

static int initALITRACK(ALITRACK *p, size_t blksz, size_t n_alloc_thresh)
{
  if (blksz < 1) blksz = ALIBKTRK_BLKSZ_DEFAULT;
  ECALLOCP(blksz, p->bdp);
  if (!(p->bdp))
    return ERRCODE_NOMEM;
 
  p->n_alloc = p->blksz = blksz;
  p->n_alloc_thresh = (n_alloc_thresh < 1)? ALIBKTRK_MEMLIM_DEFAULT: n_alloc_thresh;
  p->max_i = p->max_j = 0;
  p->max_scor = 0;

  return ERRCODE_SUCCESS;
}

static int setMemALITRACK(ALITRACK *btrkp, ALIBAND *bandp)
{
  void *hp;
  size_t newsiz;

  if (bandp->s_left >= bandp->s_len || bandp->band_width < 0)
    return ERRCODE_ASSERT;

  newsiz = bandp->band_width*(bandp->s_len - bandp->s_left);

  if (newsiz > btrkp->n_alloc || 
      (newsiz + btrkp->blksz < btrkp->n_alloc && btrkp->n_alloc > btrkp->n_alloc_thresh)) {
    newsiz = ((newsiz-1)/btrkp->blksz + 1)*btrkp->blksz;
    hp = EREALLOCP(btrkp->bdp, newsiz);
    if (!hp) return ERRCODE_NOMEM;
    btrkp->bdp = hp;
    btrkp->n_alloc = newsiz;
  }
  btrkp->max_i = btrkp->max_j = 0;
  btrkp->max_scor = 0;

  return ERRCODE_SUCCESS;
}

static void blankALITRACK(ALITRACK *trackp, const ALIBAND *bandp)
{
  int i, li;
  li = bandp->band_width*(bandp->s_len - bandp->s_left);
  for (i=0; i<li; i++)
    trackp->bdp[i] = 0;
  trackp->max_i = trackp->max_j = 0;
  trackp->max_scor = 0;
}

#ifdef alignment_debug
static int fprintALITRACK(FILE *fp, const ALIBAND *bandp, const ALITRACK *trackp)
{
  short ctr, bctr;
  int i, j, il, jl, delta_band;
  ALIDIRECT_t *dp;
  fprintf(fp, "==== ALITRACK ====\n");
  fprintf(fp, "max_scor = %i\n", trackp->max_scor);
  fprintf(fp, "max_i = %i\n", trackp->max_i);
  fprintf(fp, "max_j = %i\n", trackp->max_j);
  if (bandp) {
    if (bandp->q_left > bandp->l_edge) {
      delta_band = bandp->q_left - bandp->l_edge;
      jl = bandp->q_len - bandp->q_left;
    } else {
      delta_band = 0;
      jl = bandp->l_edge - bandp->q_left;
    }
  
    fprintf(fp, "\n     ");
    for (j=0; j<delta_band; j++)
      fprintf(fp, "   ");
    for (bctr=1,ctr=0,j=0; j<jl; j++,ctr++) {
      if (ctr > 9) {
	ctr=0;
	fprintf(fp, " %1.1i ", bctr++);
      } else {
	fprintf(fp, "   ");
      }
    }
    
    fprintf(fp, "\n     ");
    for (j=0; j<delta_band; j++)
      fprintf(fp, "   ");
    for (ctr=0, j=0; j<jl; j++) {
      fprintf(fp, " %1.1i ", ctr++);
      if (ctr > 9) ctr = 0;
    }

    fprintf(fp, "\n     ");
    for (j=0; j<jl+delta_band; j++)
      fprintf(fp, "---");
    fprintf(fp, "\n");
    
    il = bandp->q_len - bandp->q_left;
    dp = trackp->bdp;
    if (il < 0) return ERRCODE_ASSERT;
    for (i=0; i<il; i++) {
      fprintf(fp, "%3i |", i);
      for (j=0; j<i; j++) 
	fprintf(fp, "   ");
      for (j=0; j<bandp->band_width; j++) {
	fprintf(fp, " %1.1i ", (int) (*dp));
	dp++;
      }
      fprintf(fp,"\n");
    }
  }
  return ERRCODE_SUCCESS;
}
#endif

/******************************************************************************
 ************************ Methods of Private Type ALIMETA *********************
 ******************************************************************************/

static void cleanupALIMETA(ALIMETA *p)
{
#ifdef alignment_debug
  short i;
  if (p) {
    for (i=0; i<NUMSTR; i++) {
      free(p->string[i]);
      p->string[i] = 0;
    }
  }
#endif
  diffStrCleanUp(&p->dfs);
}

static int initALIMETA(ALIMETA *p, short dfstr_blksz)
{
  int errcode = diffStrInit(&p->dfs, dfstr_blksz);
#ifdef alignment_debug
  if (!errcode) {
    short i;

    for (i=0; i<NUMSTR; i++) {
      ECALLOCP(ALISTRINGS_DEFAULT_BLKSZ, p->string[i]);
      if (!p->string[i]) {
	errcode = ERRCODE_NOMEM;
	break;
      }
    }
    p->blocksiz = p->allocsiz = ALISTRINGS_DEFAULT_BLKSZ;
  }
#endif
  if (errcode) {
    cleanupALIMETA(p);
  }
  return errcode;
}
#ifdef alignment_debug
static void fprintALIMETA(FILE *fp, const ALIMETA *mp)
{
  short i;
  fprintf(fp, "=== ALIMETA ===\n");
  fprintf(fp, "prof_start = %i\nprof_end = %i\n", mp->prof_start, mp->prof_end);
  fprintf(fp, "nonprof_start = %i\nnonprof_end = %i\n", mp->nonprof_start, mp->nonprof_end);
  fprintf(fp, "score = %i\n", mp->score);
  for (i=0; i<NUMSTR; i++) {
    fprintf(fp, "STRING_%i:%s\n", i, mp->string[i]);
  }
  fprintf(fp, "DiffStr: ");
  diffStrPrintf(fp, mp->dfs.dstrp, DIFFSTRFORM_RAW, 0, 0, 0);
  fprintf(fp, "=== End of ALIMETA ===\n");
}

static int setALIMETA(ALIMETA *p, unsigned int sum_len)
{
  short i;
  if (sum_len >= p->allocsiz) {
    size_t allocsiz;
    void *hp;
    allocsiz = ((int)((sum_len+1)/p->blocksiz + 1))*p->blocksiz;
    for (i=0; i<NUMSTR; i++) {
      hp = EREALLOC(p->string[i], allocsiz*sizeof(char));
      if (!hp) return ERRCODE_NOMEM;
      p->string[i] = hp;
      p->allocsiz = allocsiz;
      p->string[i][0] = '\0';
    }
  }
  p->string_len = 0;
  p->string_end = -1;

  return ERRCODE_SUCCESS;
}
#endif

static int makeMetaFromTrack(ALIMETA *metap,
			     ALICPLX *cplxp,
			     const ALITRACK *tp,
			     const ALIBAND *bp,
			     const ScoreProfile *qp,
#ifdef alignment_debug
			     const SeqCodec *codecp,
			     const char *profiled_seq, 
#endif
			     const char *unprofiled_seq)
     /**< Make a short-hand alignment string from back-tracking directions.
      * if cplxp != NULL complexity-weight alignment score
      */
{
  int errcode = ERRCODE_SUCCESS;
  BOOL is_gap_open = FALSE;
  UCHAR nmatch;
  int i, j;
  ALIDPMSCOR_t s, checksum_score = 0; /* pos_score = 0 */
  const ALIDIRECT_t *dp;
  ALIMATSCOR_t score_gap_open, score_gap_ext;
  ALIMATSCOR_t * const *scorepp = scoreGetProfile(NULL, NULL, &score_gap_open, 
						  &score_gap_ext, qp);
  DiffStr *const dfsp = &metap->dfs;
#ifdef alignment_debug
  char profiled_symbol, unprofiled_symbol;
  char *alistr_noprof, *alistr_prof, *alistr_cons;
  const char *decoderp = seqCodecGetDecoder(codecp, NULL);
  int k=0;

  alistr_prof = metap->string[0];   /* aligned profiled sequence (query), hyphens for gaps */
  alistr_noprof = metap->string[1]; /* aligned unprofiled sequence (subject), hyphens of gaps */
  alistr_cons = metap->string[2];   /* consensus ' ' for mismatch. '+' for match, and '-' for indels (gaps) */
#endif
/* macros */
#define SETDIFF(dfsp, nmatch, typ) \
  if (((dfsp)->len >= (dfsp)->n_alloc) && \
      (errcode = diffStrRealloc((dfsp), (dfsp)->len+1))) \
    return errcode;\
  (dfsp)->dstrp[dfsp->len++] = (nmatch) + (((unsigned char) (typ)) << DIFFSTR_TYPSHIFT);

/* end of macros */ 
  DIFFSTR_LENGTH(dfsp) = 0;
  nmatch = 0;
  if (cplxp) 
    blankALICPLX(cplxp);

  i = tp->max_i - bp->s_left;
  dp = tp->bdp + i * (bp->band_width - 1) + tp->max_j - bp->l_edge;
  i = tp->max_i;
  j = tp->max_j;
  for (i=tp->max_i, j=tp->max_j;  i >= bp->s_left && j >= bp->q_left && (*dp);) { 
    if (*dp == ALIBKTRKCOD_DIA) {
      /* diagonal move */
      s = scorepp[(unprofiled_seq[i]&SEQCOD_ALPHA_MASK)][j];
      if (s > 0) {
	/* match */
	if (nmatch > DIFFSTR_MAXMISMATCH) {
	  SETDIFF(dfsp, DIFFSTR_MAXMISMATCH, DIFFCOD_M);
	  nmatch -= DIFFSTR_MAXMISMATCH;
	} else {
	  nmatch++;
	}
      } else {
	/* substitution (mismatch) */
	SETDIFF(dfsp, nmatch, DIFFCOD_S);
	nmatch = 0;
      }
      checksum_score += s;
      /* pos_score += s; */
#ifdef alignment_debug
      unprofiled_symbol = decoderp[(UCHAR) unprofiled_seq[i]];
      profiled_symbol = decoderp[(UCHAR) profiled_seq[j]];
      alistr_prof[k] = profiled_symbol;
      alistr_noprof[k] = unprofiled_symbol;
      alistr_cons[k] = (profiled_symbol == unprofiled_symbol)? profiled_symbol: ALIMETA_MISMATCH;
      k++;
#endif
      if (cplxp)
	cplxp->countp[unprofiled_seq[i]&SEQCOD_ALPHA_MASK]++;
      is_gap_open = FALSE;
      dp -= bp->band_width;
      i--;
      j--;
      continue;
    }

    if ((is_gap_open)) {
      checksum_score -= score_gap_ext;
    } else {
      checksum_score -= score_gap_open;	
      is_gap_open = TRUE;
    }

    if ((*dp) & ALIBKTRKCOD_COL) {
      /* move along column (gap in profiled sequence, deletion) */
      SETDIFF(dfsp, nmatch, DIFFCOD_D);
      nmatch = 0;
#ifdef alignment_debug
      alistr_prof[k] = ALIMETA_SPACER;
      alistr_noprof[k] = decoderp[(UCHAR) unprofiled_seq[i]];
      alistr_cons[k] = ALIMETA_SPACER;
      k++;
#endif
      dp -= bp->band_width-1;
      i--;
      continue;
    } 
    if (!((*dp) & ALIBKTRKCOD_ROW))
      return ERRCODE_ASSERT;
    
    /* move along row */
    SETDIFF(dfsp, nmatch, DIFFCOD_I);
    nmatch = 0;
#ifdef alignment_debug
    alistr_prof[k] = decoderp[(UCHAR) profiled_seq[j]];
    alistr_noprof[k] = ALIMETA_SPACER;
    alistr_cons[k] = ALIMETA_SPACER;
    k++;
#endif 
    dp--;
    j--;
  }

  SETDIFF(dfsp, nmatch, DIFFCOD_S);
  SETDIFF(dfsp, 0, DIFFCOD_M);

#ifdef alignment_debug
  alistr_prof[k] = '\0';
  alistr_noprof[k] = '\0';
  alistr_cons[k] = '\0';
  metap->string_end = k - 1;
#endif

  metap->nonprof_start = i+1;
  metap->nonprof_end = tp->max_i;
  metap->prof_start = j+1;
  metap->prof_end = tp->max_j;

  if (checksum_score != tp->max_scor) 
    errcode = ERRCODE_SWATSCOR;
  else if (cplxp)
    errcode = scaleALICPLX(&checksum_score, tp->max_scor, 
			   /* pos_score, */
			   cplxp);
   
  metap->score = checksum_score;

#ifdef alignment_debug
  fprintALIMETA(stdout, metap);
#endif

  return errcode;
}

/******************************************************************************
 ***************************** Alignment Methods ******************************
 ******************************************************************************/


static int alignSmiWatBand(ALITRACK *bktp,
			   AliBuffer *bufp,
			   const ALIBAND *bandp,
			   const ScoreProfile *profp,
#ifdef alignment_matrix_debug
			   const char *psqp,
			   const SeqCodec *codecp,
#endif
			   const char *usqp)
{
  int errcode;
  int delta_band_start, delta_band_end;
  int i, j, j_curr_start, j_curr_len;
  int max_i = 0, max_j = 0, max_scor = 0;
  const ALIMATSCOR_t *rowscorp;
  ALIMATSCOR_t *const *scorpp;
  ALIMATSCOR_t gap_ext;
  ALIMATSCOR_t gap_init;
  ALIDIRECT_t *dirp;
  ALIDPMSCOR_t H, F, tmp, currH, *Hp, *Ep;
#ifdef alignment_debug
  unsigned int splen;
#endif
#ifdef alignment_matrix_debug
  char unprof_symbol;
  const char *decoderp = seqCodecGetDecoder(codecp, NULL);
  ALIDIRECT_t *dir_refp;
  ALIDPMSCOR_t *Fp;
#endif

#ifdef alignment_matrix_debug
#define RECORD_MAXIMUM_SCORE(H) if ((H) > max_scor)\
{ max_i = i; \
  max_j = j; \
  max_scor = (H);\
  printf("max_scor(%i,%i) = %i\n", i, j, max_scor);\
}
#else
#define RECORD_MAXIMUM_SCORE(H) if ((H) > max_scor)\
{ max_i = i; \
  max_j = j; \
  max_scor = (H);\
}
#endif
  scorpp = scoreGetProfile(NULL, 
#ifdef alignment_debug
			   &splen,
#else
			   NULL, 
#endif
			   &gap_init, &gap_ext, profp);

#ifdef alignment_debug
  if ((unsigned int) bandp->q_len > splen)
     printf("ERROR: wrong band paramter q_len = %i, should be <= %u\n!",
	    bandp->q_len, splen);
#endif
#ifdef alignment_matrix_debug
  ECALLOCP(bandp->band_width + 1, Fp);
  if (!Fp)
     return ERRCODE_NOMEM;
  for (j=0; j<bandp->band_width+1; j++)
     Fp[j] = 0;
#endif
  if (bandp->q_len >= bufp->qlen_max &&
      (errcode = aliBufferInit(bufp, (unsigned int) bandp->q_len)))
    return errcode;

  if (bandp->q_left > bandp->l_edge) {
    delta_band_start = bandp->q_left - bandp->l_edge;
    j_curr_start = bandp->q_left;
  } else {
    delta_band_start = 0;
    j_curr_start = bandp->l_edge;
  }
  j_curr_len = bandp->r_edge + 1;
  delta_band_end = 0;
  H = currH = 0; /* H[i-1,j-1] */
  dirp = bktp->bdp + delta_band_start;

  Hp = bufp->baseHp;
  Ep = bufp->baseEp;
  /* initialise first row */
  for (j=j_curr_start; j<bandp->q_len; j++)
    Hp[j] = Ep[j] = 0;

  for (i = bandp->s_left; i<bandp->s_len; i++) {
    F = 0; /* at the beginning of the row */
    rowscorp = scorpp[usqp[i]&SEQCOD_ALPHA_MASK];
#ifdef alignment_matrix_debug
    dir_refp = dirp;
#endif
    for (j=j_curr_start; j<j_curr_len; j++, dirp++) {
#ifdef alignment_matrix_debug
      Fp[j-j_curr_start] = F;
#endif
      /* inner loop */
      H = currH + rowscorp[j]; /* H[i-1,j-1] + W[i,j] */
      currH = Hp[j];        
      if (F > 0) { 
	/* F[i,j] = max(F[i,j-1] - gap_ext, H[i,j-1] - gap_init) > 0 */
	if (Ep[j] > 0) { 
	  /****************
	   * E > 0, F > 0 *
	   ****************/	  
	  /* E[i,j] = max(E[i-1,j] - gap_ext, H[i-1,j] - gap_init) > 0 */
	  if (H > Ep[j]) { 
	    /* H[i-1,j-1] + W[i,j] > E */
	    if (H > F) { 
	      /* max(H, E, F, 0) = H[i-1,j-1] + W[i,j] */
	      Hp[j] = H; /* H[i,j] = H[i-1,j-1] + W[i,j] */
	      F -= gap_ext;
	      Ep[j] -= gap_ext;
	      *dirp = ALIBKTRKCOD_DIA;
	      if (H > gap_init) { /* H[i,j] - gap_init > 0 */
		RECORD_MAXIMUM_SCORE(H);
		tmp = H - gap_init;
		if (F < tmp) F = tmp; /* max(F[i,j] - gap_ext, H[i,j] - gap_init) */
		if (Ep[j] < tmp) Ep[j] = tmp; /* E[i,j] = max(E[i,j] - gap_ext, H[i,j] - gap_init) */
	      }
	    } else { 
	      /* F >= H[i-1,j-1] + W[i,j] > E, i.e. F > E */
	      Hp[j] = F;
	      F -= gap_ext;
	      Ep[j] -= gap_ext;
	      *dirp = ALIBKTRKCOD_ROW;
	    }
	  } else {
	    /* H[i-1,j-1] + W[i,j] <= E, max(E,F) undetermined */
	    if (Ep[j] >= F) {
	      Hp[j] = Ep[j];
	      *dirp = ALIBKTRKCOD_COL;
	    } else {
	      /* F > E */
	      Hp[j] = F;
	      *dirp = ALIBKTRKCOD_ROW;
	    }
	    Ep[j] -= gap_ext;
	    F -= gap_ext;
	  }
	} else { 
	  /*****************
	   * E <= 0, F > 0 *
	   *****************/	  
	  if (H > F) {
	    Hp[j] = H;
	    F -= gap_ext;
	    *dirp = ALIBKTRKCOD_DIA;
	    if (H > gap_init) {
	      RECORD_MAXIMUM_SCORE(H);
	      Ep[j] = H - gap_init;
	      if (F < Ep[j]) F = Ep[j]; /* F = max(F[i,j] - gap_ext, H[i,j] - gap_init) */
	    } 	  
	  } else {
	    /* H[i-1,j-1] + W[i,j] <= max(F[i,j-1] - gap_ext, H[i,j-1] - gap_init) */
	    Hp[j] = F;
	    F -= gap_ext;
	    *dirp = ALIBKTRKCOD_ROW;
	  }
	}
      } else if (Ep[j] > 0) {
	/*****************
	 * E > 0, F <= 0 *
	 *****************/	  
	if (H > Ep[j]) {
	  Hp[j] = H;
	  Ep[j] -= gap_ext;
	  *dirp = ALIBKTRKCOD_DIA;
	  if (H > gap_init) {
	    RECORD_MAXIMUM_SCORE(H);
	    F = H - gap_init;
	    if (Ep[j] < F) Ep[j] = F;
	  } 
	} else {
	  /* H[i-1,j-1] + W[i,j] <= max(E[i-1,j] - gap_ext, H[i,j-1] - gap_init) */
	  Hp[j] = Ep[j];
	  Ep[j] -= gap_ext;
	  *dirp = ALIBKTRKCOD_COL;
	}
      } else {
	/******************
	 * E <= 0, F <= 0 *
	 ******************/	  
	if (H > 0) {
	  Hp[j] = H;
	  *dirp = ALIBKTRKCOD_DIA;
	  if (H > gap_init) {
	    RECORD_MAXIMUM_SCORE(H);
	    F = Ep[j] = H - gap_init;
	  }
	} else {
	  Hp[j] =0;
	  *dirp = 0;
	}
      }
    } /* end of j-loop */
#ifdef alignment_matrix_debug
    printf("  [%i,%i]: ", i, j_curr_start);
    for (j=j_curr_start; j<j_curr_len; j++)
      printf("%4i", j);
    unprof_symbol = decoderp[(UCHAR) usqp[i]];
    printf("\n%c [%i,%i]: |", unprof_symbol, i, j_curr_start);
    for (j=j_curr_start; j<j_curr_len; j++)
      printf(" %c |", decoderp[(UCHAR) psqp[j]]);
    printf("\nHp[%i,%i]: |", i, j_curr_start);
    for (j=j_curr_start; j<j_curr_len; j++)
      printf("%3i|", Hp[j]);
    printf("\nEp[%i,%i]: |", i, j_curr_start);
    for (j=j_curr_start; j<j_curr_len; j++)
      printf("%3i|", Ep[j]);  
    printf("\nFp[%i,%i]: |", i, j_curr_start);
    for (j=j_curr_start; j<j_curr_len; j++)
      printf("%3i|", Fp[j-j_curr_start]);      
    printf("\ndi[%i,%i]: |", i, j_curr_start);
    for (; dir_refp < dirp; dir_refp++)
      printf(" %1i |", *dir_refp);
    printf("\n\n");
#endif
    if (delta_band_start > 0) {
      currH = 0;
      dirp += --delta_band_start;
    } else {
      currH = Hp[j_curr_start];
      j_curr_start++;
    }
    if (j_curr_len < bandp->q_len) {
      j_curr_len++;
    } else {
      dirp += delta_band_end++;
    }
      
  }
  bktp->max_i = max_i;
  bktp->max_j = max_j;
  bktp->max_scor = max_scor;
#ifdef alignment_matrix_debug
  free(Fp);
#endif
  return ERRCODE_SUCCESS;
}

static int alignSmiWatBandFast(int *maxswscor,
			       AliBuffer * const bufp,
			       const ALIBAND * const bandp,
			       const ScoreProfile * const profp,
#ifdef alignment_matrix_debug
			       const char * const psqp,
			       const SeqCodec * const codecp,
#endif
			       const char * const usqp)
{
  int errcode;
  int delta_band_start;
  int i, j, j_curr_start, j_curr_len;
  int max_scor = 0;
  const ALIMATSCOR_t *rowscorp;
  ALIMATSCOR_t *const *scorpp;
  ALIMATSCOR_t gap_ext;
  ALIMATSCOR_t gap_init;
  ALIDPMSCOR_t H, F, tmp, currH, *Hp, *Ep;
#ifdef alignment_debug
  unsigned int splen;
#endif
#ifdef alignment_matrix_debug
  char unprof_symbol;
  const char * const decoderp = seqCodecGetDecoder(codecp, NULL);
  ALIDIRECT_t *dir_refp;
  ALIDPMSCOR_t *Fp;
#endif

#define RECORD_MAXIMUM_SCORE_FAST(H) if ((H) > max_scor) max_scor = (H);

  *maxswscor = 0;
  scorpp = scoreGetProfile(NULL, 
#ifdef alignment_debug
			   &splen,
#else
			   NULL, 
#endif
			   &gap_init, &gap_ext, profp);

#ifdef alignment_debug
  if (((unsigned int) bandp->q_len) > splen)
     printf("ERROR: wrong band paramter q_len = %i, should be <= %u\n!",
	    bandp->q_len, splen);
#endif
#ifdef alignment_matrix_debug
  ECALLOCP(bandp->band_width + 1, Fp);
  if (!Fp)
     return ERRCODE_NOMEM;
  for (j=0; j<bandp->band_width+1; j++)
     Fp[j] = 0;
#endif
  if (bandp->q_len >= bufp->qlen_max &&
      (errcode = aliBufferInit(bufp, (unsigned int) bandp->q_len)))
    return errcode;

  if (bandp->q_left > bandp->l_edge) {
    delta_band_start = bandp->q_left - bandp->l_edge;
    j_curr_start = bandp->q_left;
  } else {
    delta_band_start = 0;
    j_curr_start = bandp->l_edge;
  }
  j_curr_len = bandp->r_edge + 1;
  H = currH = 0; /* H[i-1,j-1] */

  Hp = bufp->baseHp;
  Ep = bufp->baseEp;
  /* initialise first row */
  for (j=j_curr_start; j<bandp->q_len; j++)
    Hp[j] = Ep[j] = 0;

  for (i = bandp->s_left; i<bandp->s_len; i++) {
    F = 0; /* at the beginning of the row */
    rowscorp = scorpp[usqp[i]&SEQCOD_ALPHA_MASK];

    for (j=j_curr_start; j<j_curr_len; j++) {
#ifdef alignment_matrix_debug
      Fp[j-j_curr_start] = F;
#endif
      /* inner loop */
      H = currH + rowscorp[j]; /* H[i-1,j-1] + W[i,j] */
      currH = Hp[j];        
      if (F > 0) { 
	/* F[i,j] = max(F[i,j-1] - gap_ext, H[i,j-1] - gap_init) > 0 */
	if (Ep[j] > 0) { 
	  /****************
	   * E > 0, F > 0 *
	   ****************/	  
	  /* E[i,j] = max(E[i-1,j] - gap_ext, H[i-1,j] - gap_init) > 0 */
	  if (H > Ep[j]) { 
	    /* H[i-1,j-1] + W[i,j] > E */
	    if (H > F) { 
	      /* max(H, E, F, 0) = H[i-1,j-1] + W[i,j] */
	      Hp[j] = H; /* H[i,j] = H[i-1,j-1] + W[i,j] */
	      F -= gap_ext;
	      Ep[j] -= gap_ext;
	      if (H > gap_init) { /* H[i,j] - gap_init > 0 */
		RECORD_MAXIMUM_SCORE_FAST(H);
		tmp = H - gap_init;
		if (F < tmp) F = tmp; /* max(F[i,j] - gap_ext, H[i,j] - gap_init) */
		if (Ep[j] < tmp) Ep[j] = tmp; /* E[i,j] = max(E[i,j] - gap_ext, H[i,j] - gap_init) */
	      }
	    } else { 
	      /* F >= H[i-1,j-1] + W[i,j] > E, i.e. F > E */
	      Hp[j] = F;
	      F -= gap_ext;
	      Ep[j] -= gap_ext;
	    }
	  } else {
	    /* H[i-1,j-1] + W[i,j] <= E, max(E,F) undetermined */
	    if (Ep[j] >= F) {
	      Hp[j] = Ep[j];
	    } else {
	      /* F > E */
	      Hp[j] = F;
	    }
	    Ep[j] -= gap_ext;
	    F -= gap_ext;
	  }
	} else { 
	  /*****************
	   * E <= 0, F > 0 *
	   *****************/	  
	  if (H > F) {
	    Hp[j] = H;
	    F -= gap_ext;
	    if (H > gap_init) {
	      RECORD_MAXIMUM_SCORE_FAST(H);
	      Ep[j] = H - gap_init;
	      if (F < Ep[j]) F = Ep[j]; /* F = max(F[i,j] - gap_ext, H[i,j] - gap_init) */
	    } 	  
	  } else {
	    /* H[i-1,j-1] + W[i,j] <= max(F[i,j-1] - gap_ext, H[i,j-1] - gap_init) */
	    Hp[j] = F;
	    F -= gap_ext;
	  }
	}
      } else if (Ep[j] > 0) {
	/*****************
	 * E > 0, F <= 0 *
	 *****************/	  
	if (H > Ep[j]) {
	  Hp[j] = H;
	  Ep[j] -= gap_ext;
	  if (H > gap_init) {
	    RECORD_MAXIMUM_SCORE_FAST(H);
	    F = H - gap_init;
	    if (Ep[j] < F) Ep[j] = F;
	  } 
	} else {
	  /* H[i-1,j-1] + W[i,j] <= max(E[i-1,j] - gap_ext, H[i,j-1] - gap_init) */
	  Hp[j] = Ep[j];
	  Ep[j] -= gap_ext;
	}
      } else {
	/******************
	 * E <= 0, F <= 0 *
	 ******************/	  
	if (H > 0) {
	  Hp[j] = H;
	  if (H > gap_init) {
	    RECORD_MAXIMUM_SCORE_FAST(H);
	    F = Ep[j] = H - gap_init;
	  }
	} else {
	  Hp[j] =0;
	}
      }
    } /* end of j-loop */
#ifdef alignment_matrix_debug
    printf("  [%i,%i]: ", i, j_curr_start);
    for (j=j_curr_start; j<j_curr_len; j++)
      printf("%4i", j);
    unprof_symbol = decoderp[(UCHAR) usqp[i]];
    printf("\n%c [%i,%i]: |", unprof_symbol, i, j_curr_start);
    for (j=j_curr_start; j<j_curr_len; j++)
      printf(" %c |", decoderp[(UCHAR) psqp[j]]);
    printf("\nHp[%i,%i]: |", i, j_curr_start);
    for (j=j_curr_start; j<j_curr_len; j++)
      printf("%3i|", Hp[j]);
    printf("\nEp[%i,%i]: |", i, j_curr_start);
    for (j=j_curr_start; j<j_curr_len; j++)
      printf("%3i|", Ep[j]);  
    printf("\nFp[%i,%i]: |", i, j_curr_start);
    for (j=j_curr_start; j<j_curr_len; j++)
      printf("%3i|", Fp[j-j_curr_start]);      
    printf("\n\n");
#endif
    if (delta_band_start > 0) {
      currH = 0;
    } else {
      currH = Hp[j_curr_start];
      j_curr_start++;
    }
    if (j_curr_len < bandp->q_len) {
      j_curr_len++;
    }       
  }
  *maxswscor = max_scor;
#ifdef alignment_matrix_debug
  free(Fp);
#endif
  return ERRCODE_SUCCESS;
}

/******************************************************************************
 ********************** Private Methods of Type AliRsltSet ********************
 ******************************************************************************/

static int reallocRsltSet(AliRsltSet *arp, short newsiz)
{
  int errcode = ERRCODE_SUCCESS;
  void *hp;
  short i;
  size_t ns = (newsiz + arp->blksz - 1)/arp->blksz;

  if (newsiz < 0) 
    return ERRCODE_ASSERT;

  ns *= arp->blksz;

  if (ns > SHRT_MAX)
    return ERRCODE_OVERFLOW;

  if ((short) ns < newsiz)
    return ERRCODE_ASSERT;

  for (i=ns; i<arp->n_alloc; i++)
    diffStrCleanUp(&(arp->rsp[i].diffstr));

  hp = EREALLOCP(arp->rsp, ns);
  if (!hp)
    return ERRCODE_NOMEM;
  
  arp->rsp = (ALIRESULT *) hp;
  for (i=arp->n_alloc;i<(short) ns; i++) {
    if ((errcode = diffStrInit(&(arp->rsp[i].diffstr), (int) arp->dfblksz)))
      break;
  }
  arp->n_alloc = (short) ns;
  if (arp->nres > arp->n_alloc)
    arp->nres = arp->n_alloc;
  
  return errcode;
}


static int addALIMETAtoRsltSet(AliRsltSet *p)
{
  int errcode;
  ALIRESULT *arp;
  ALIMETA *mp = &p->meta;

  if (p->nres >= p->n_alloc &&
      (errcode = reallocRsltSet(p, p->nres + 1)))
      return errcode;

  arp = p->rsp + p->nres;
  
  arp->score = mp->score;
  arp->qs = mp->prof_start;
  arp->qe = mp->prof_end;
  arp->rs = mp->nonprof_start;
  arp->re = mp->nonprof_end; 
  diffStrReverse(&arp->diffstr, mp->dfs.dstrp);

  p->nres++;
  return ERRCODE_SUCCESS;
}

static int alignSmiWatBandRecursive(
#ifdef alignment_debug
			  short *level,
#endif
			  AliRsltSet *rssp,
			  AliBuffer *bufp,
			  const ScoreProfile *q_profp,
#if defined alignment_debug || defined algnment_matrix_debug
			  const SeqCodec *codecp,
			  const char *q_seqp,
#endif

			  int q_len,
			  const char *s_seqp, int s_len,
			  int l_edge, int r_edge,
			  int q_left, int q_right,
			  int s_left, int s_right,
			  const ALIDPMSCOR_t minscore,
			  const int minscorlen
			  )
{
  int errcode = ERRCODE_SUCCESS;
  int s_start, s_end;
  ALIBAND band;

#ifdef alignment_debug  
  ++*level;
#endif
  if (minscorlen < 2)
    return ERRCODE_ASSERT;

  /* return without raising error if the limits and the alignment
   * band are inconsistent (end of recursion). */
  if (initALIBAND(&band, 
		  l_edge, r_edge,
		  q_left, q_right, q_len,
		  s_left, s_right, s_len
		  ))
    return ERRCODE_SUCCESS;
  
#ifdef alignment_debug_limits
  printf("swat.c::SwatRecursive:\n");
  printf(">Query\n");
  printEncSeq(stdout, q_seqp);
  printf("\n>Subject\n");
  printEncSeq(stdout, s_seqp);
  printf("\ns_len = %i\n", s_len);
  printf("l_edge = %i, r_edge = %i\n", l_edge, r_edge);
  printf("s_left = %i, s_right = %i\n", s_left, s_right);
  printf("q_left = %i, q_right = %i\n", q_left, q_right);
  printAliBand(stdout, &band);
#endif
  if ((errcode = setMemALITRACK(&rssp->track, &band)))
    return errcode;

  if ((errcode = alignSmiWatBand(&rssp->track, bufp, &band, 
				 q_profp, 
#ifdef alignment_matrix_debug
				 q_seqp,
				 codecp,
#endif
				 s_seqp)))
    return errcode;

  if (rssp->track.max_scor < minscore)
    return ERRCODE_SUCCESS;
    
  errcode = makeMetaFromTrack(&rssp->meta, rssp->cplxp, &rssp->track, 
			      &band, q_profp,
#ifdef alignment_debug
			      codecp,
			      q_seqp, 
#endif
			      s_seqp);

  if (errcode)
    return errcode;

  /* exit if aligned segment of the query reads is not at least minscorlen nucleotides long */
  if (rssp->meta.prof_start + minscorlen > rssp->meta.prof_end + 1)
    return ERRCODE_SUCCESS;

  s_start = rssp->meta.nonprof_start;
  s_end = rssp->meta.nonprof_end;
  if (rssp->meta.score >= minscore) {
    if ((errcode = addALIMETAtoRsltSet(rssp)))
      return errcode;
  }

  if (s_left + minscorlen < s_start) {
    errcode = alignSmiWatBandRecursive(
#ifdef alignment_debug
			     level,
#endif
			     rssp, bufp, q_profp,
#if defined alignment_debug || defined algnment_matrix_debug
			     codecp,
			     q_seqp,
#endif
			     q_len,
			     s_seqp, s_len,
			     l_edge, r_edge, 
			     q_left, q_right, 
			     s_left, s_start-1,
			     minscore, 
			     minscorlen
			     );
    if (errcode) 
      return errcode;
  }

  if(s_right > s_end + minscorlen) {
    errcode = alignSmiWatBandRecursive(
#ifdef alignment_debug
			     level,
#endif
			     rssp, bufp, q_profp,
#if defined alignment_debug || defined algnment_matrix_debug
			     codecp,
			     q_seqp,
#endif
			     q_len,
			     s_seqp, s_len,
			     l_edge, r_edge, 
			     q_left, q_right,
			     s_end+1, s_right,
			     minscore,
			     minscorlen
			     );
    if (errcode) 
      return errcode;
  }

  return ERRCODE_SUCCESS;
}

/******************************************************************************
 *********************** Public Methods of Type AliRsltSet ********************
 ******************************************************************************/

AliRsltSet *aliRsltSetCreate(const ScoreMatrix *smp, 
			     short blksz, short diffblksz,
			     int track_blksz, int track_thresh)
{
  int errcode = ERRCODE_SUCCESS;
  short i;
  AliRsltSet *p;
 
  EMALLOCP0(p);
  if (!p) return NULL;

  if (blksz < 1) 
    blksz = ALIRSLTSET_DEFAULT_BLKSZ;

  ECALLOCP(blksz, p->rsp);
  if (!(p->rsp)) {
    aliRsltSetDelete(p);
    return 0;
  }
  for (i=0; i<blksz; i++)
    if((errcode = diffStrInit(&p->rsp[i].diffstr, (int) diffblksz)))
      break;
  if (errcode) {
    aliRsltSetDelete(p);
    return 0;
  }

  p->n_alloc = p->blksz = blksz;

  if ((smp)) {
    p->cplxp = createALICPLX(smp);
    if (!p->cplxp) {
      aliRsltSetDelete(p);
      return 0;
    }
  } else {
    p->cplxp = 0;
  }

  if ((errcode = initALITRACK(&p->track, track_blksz, track_thresh)) ||
      (errcode = initALIMETA(&p->meta, diffblksz))) {
    aliRsltSetDelete(p);
    p = 0;
  }

  return p;
}

void aliRsltSetDelete(AliRsltSet *p)
{
  short i;
  if (p) {
    deleteALICPLX(p->cplxp);
    p->cplxp = 0;
    cleanupALIMETA(&p->meta);
    cleanupALITRACK(&p->track);
    for (i=0; i<p->n_alloc; i++) {
      diffStrCleanUp(&p->rsp[i].diffstr);
    }
    free(p->rsp);
  }
  free(p);
}

void aliRsltSetReset(AliRsltSet *p)
{
  short i;
  for (i=0; i>p->nres;i++)
    p->rsp[i].diffstr.len=0;
  p->nres = 0;
  blankALICPLX(p->cplxp);
}

short aliRsltSetGetSize(const AliRsltSet *arp)
{
  return arp->nres;
}

int aliRsltSetFetchData(const AliRsltSet *arp, short idx, int *score, 
		       int *ps_start, int *ps_end, int *us_start, int *us_end,
		       const DiffStr **dfsp)
{
  const ALIRESULT *rp;
  if (idx >= arp->nres) 
    return ERRCODE_FAILURE;

  rp = arp->rsp + idx;
  
  if (score)
    *score = rp->score;
  if (ps_start)
    *ps_start = rp->qs;
  if (ps_end)
    *ps_end = rp->qe;
  if (us_start)
    *us_start = rp->rs;
  if (us_end)
    *us_end = rp->re;
  if (dfsp)
    *dfsp = &rp->diffstr;

  return ERRCODE_SUCCESS;
}

/******************************************************************************
 ************************* Public Alignment Methods ***************************
 ******************************************************************************/

int aliSmiWatInBand(AliRsltSet *rssp,
		    AliBuffer *bufp,
		    const ScoreProfile *profp,
#if defined alignment_debug || defined algnment_matrix_debug
		    const SeqCodec *codecp,
		    const char *profiled_seqp,
#endif
		    const char *unprofiled_seqp,
		    int unprofiled_seqlen,
		    int l_edge, int r_edge,
		    int profiled_left, int profiled_right,
		    int unprofiled_left, int unprofiled_right,
		    ALIDPMSCOR_t minscore, int minscorlen)
{
#ifdef alignment_debug
  short level = 0;
  int errcode;
#endif
  unsigned int qlen;
  const short matchscor = scoreProfileGetAvgPenalties(NULL, NULL, NULL, profp);

  if (minscore < 1 || matchscor <= 0)
    return ERRCODE_ASSERT;

  if (minscorlen*matchscor < minscore)
    minscorlen =  minscore/matchscor;
  if (minscorlen < ALILEN_MIN)
    return ERRCODE_ASSERT;

  scoreGetProfile(NULL, &qlen, NULL, NULL, profp);
  if (qlen > INT_MAX)
    return ERRCODE_SEQLEN;
  
#ifdef alignment_debug
  if ((errcode = setALIMETA(&rssp->meta, unprofiled_seqlen + qlen + 1)))
    return errcode;
#endif
  return alignSmiWatBandRecursive(
#ifdef alignment_debug
			&level,
#endif
			rssp, bufp, profp,
#if defined alignment_debug || defined algnment_matrix_debug
			codecp,
			profiled_seqp,
#endif
			qlen,
			unprofiled_seqp, unprofiled_seqlen,
			l_edge, r_edge,
			profiled_left, profiled_right,
			unprofiled_left, unprofiled_right,
			minscore, minscorlen
			);
}

int aliSmiWatInBandFast(ALIDPMSCOR_t *maxswscor,
			AliBuffer *bufp,
			const ScoreProfile *profp,
#if defined algnment_matrix_debug
			const SeqCodec *codecp,
			const char *profiled_seqp,
#endif
			const char *unprofiled_seqp,
			int unprofiled_seqlen,
			int l_edge, int r_edge,
			int profiled_left, int profiled_right,
			int unprofiled_left, int unprofiled_right
			)
{
  int errcode;
  unsigned int qlen;
  ALIBAND band;
  
  scoreGetProfile(NULL, &qlen, NULL, NULL, profp);
  if ((errcode = initALIBAND(&band, 
			     l_edge, r_edge,
			     profiled_left, profiled_right, qlen,
			     unprofiled_left, unprofiled_right, unprofiled_seqlen
			     )))
    return errcode;
  
  errcode = alignSmiWatBandFast(maxswscor, bufp, &band, 
				profp, 
#ifdef alignment_matrix_debug
				profiled_seqp,
				codecp,
#endif
				unprofiled_seqp);

  return errcode;
}

#ifdef alignment_debug
int aliDebugFullSmiWat(AliRsltSet *rssp, AliBuffer *bufp,
		       const ScoreProfile *profp,
		       const SeqCodec *codecp,
		       const char *psqp,
		       const char *usqp, int us_len,
		       int l_edge, int r_edge,
		       int ps_start, int ps_end,
		       int us_start, int us_end)
{
  int errcode;
  unsigned int ps_len;
  ALIBAND band;
#ifdef SCORE_SIMD
  int  fastscor;
#endif
  scoreGetProfile(NULL, &ps_len, NULL, NULL, profp);

  if (ps_len > INT_MAX)
    return ERRCODE_SEQLEN;
  
  if ((errcode = setALIMETA(&rssp->meta, us_len + ps_len + 1)))
    return errcode;

  if ((errcode = initALIBAND(&band, 
			     l_edge, r_edge,
			     ps_start, ps_end, (int) ps_len,
			     us_start, us_end, us_len)))
    return errcode;

  
  if (!(errcode = setMemALITRACK(&rssp->track, &band))) {
    errcode = alignSmiWatBand(&rssp->track, bufp, 
			      &band, profp, 
#ifdef alignment_matrix_debug
			      psqp,
			      codecp,
#endif 
			      usqp);
  }
  
  if (errcode)
    return errcode;
  errcode = makeMetaFromTrack(&rssp->meta, NULL, &rssp->track, 
			      &band, profp, 
			      codecp,
			      psqp, usqp);

#ifdef SCORE_SIMD
  if (errcode)
    return errcode;

  if (us_end > us_start){
    errcode = swSIMDAlignStriped(&fastscor, bufp, profp,
#ifdef alignment_matrix_debug
				 codecp,
				 psqp,
#endif
				 usqp+us_start, us_end - us_start + 1);
  } else
    errcode = ERRCODE_ASSERT;
  if (!errcode && fastscor < rssp->track.max_scor) {
    printf("ERROR: scores conventional = %i, sse2 = %i\n", rssp->track.max_scor, fastscor);
    errcode = ERRCODE_ASSERT;
  }
#endif

  return errcode;
}
#endif

#ifdef alignment_timing
int aliSmiWatInBandDirect(AliRsltSet *rssp, AliBuffer *bufp,
			  const ScoreProfile *profp,
			  const char *psqp,
			  const char *usqp, int us_len,
			  int l_edge, int r_edge,
			  int ps_start, int ps_end,
			  int us_start, int us_end,
			  ALIDPMSCOR_t *maxscor)
{
  int errcode;
  unsigned int ps_len;
  ALIBAND band;

  *maxscor = 0;
  scoreGetProfile(NULL, &ps_len, NULL, NULL, profp);

  if (ps_len > INT_MAX)
    return ERRCODE_SEQLEN;

#ifdef alignment_debug  
  if ((errcode = setALIMETA(&rssp->meta, us_len + ps_len + 1)))
    return errcode;
#endif
  if ((errcode = initALIBAND(&band, 
			     l_edge, r_edge,
			     ps_start, ps_end, (int) ps_len,
			     us_start, us_end, us_len)))
    return errcode;

  
  if (!(errcode = setMemALITRACK(&rssp->track, &band))) {
    errcode = alignSmiWatBand(&rssp->track, bufp, 
			      &band, profp, 
#ifdef alignment_matrix_debug
			      psqp,
			      codecp,
#endif 
			      usqp);
  }
  
  if (!errcode)
    *maxscor = rssp->track.max_scor;
  return errcode;
}
#endif //#ifdef alignment_timing
		  
