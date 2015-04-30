/** Score matrix and sequence profile for dynamic programming */

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
#include <math.h>

#include "elib.h"
#include "score.h"

enum {  
  MAXNUM_3BIT = 0x07,             /**< size of a 3-bit alphabet */
  MINALPHABET = 4,                /**< smallest alphabet */
  DEFAULT_BLOCKSIZ_PROFILE = 256, /**<Granularity for memory re-allocation */
};

enum DEFAULT_ALIGNMENT_CODES {
  DEFAULT_MATCH    = 1,
  DEFAULT_MISMATCH = -2,
  DEFAULT_XMATCH   = -3, /**< MISMATCH - 1 */
  DEFAULT_GAPINIT  = -4, /**< MISMATCH - 2 */
  DEFAULT_GAPEXT   = -3, /**< MISMATCH - 1 */
};

static const double LN0P25 = -1.386294; /* natural logarithm of 1/4 */

typedef signed char ALIMATSCOR_t;
typedef int DPMSCOR_t;
typedef unsigned char UCHAR;
typedef unsigned short USHRT;
typedef unsigned int UINT32;

/******************************************************************************
 *************************** Type ScorePenalties ******************************
 ******************************************************************************/

struct ScorePenalties_ {
  ALIMATSCOR_t penalty[SCORPNLTYP_NUM];
};

/******************************************************************************
 ******************** Public Methods of Type ScorePenalties *******************
 ******************************************************************************/

ScorePenalties *scorePenaltiesCreate(void)
{
  ScorePenalties *p;
  EMALLOCP0(p);
  if (p) {
    p->penalty[SCORPNLTYP_MATCH] = DEFAULT_MATCH;
    p->penalty[SCORPNLTYP_MISMATCH] = DEFAULT_MISMATCH;
    p->penalty[SCORPNLTYP_GAPOPEN] = DEFAULT_GAPINIT;
    p->penalty[SCORPNLTYP_GAPEXT] = DEFAULT_GAPEXT;
  } else {
    scorePenaltiesDelete(p);
    p = NULL;
  }
  return p;
}

void scorePenaltiesDelete(ScorePenalties *p)
{
  free(p);
}

int scoreSetPenalty(ScorePenalties *p, short typ, short penalty)
{
  int errcode = ERRCODE_SUCCESS;

  if (typ >= 0 && typ < SCORPNLTYP_NUM) {
    ALIMATSCOR_t pnlty_min, pnlty_max;
    if (SCORPNLTYP_MATCH == typ) {
      pnlty_min = 0;
      pnlty_max = SCORPNLTY_MAXVAL;
    } else {
      pnlty_min = -1*SCORPNLTY_MAXVAL ;
      pnlty_max = 0;
    }
    if (penalty <= pnlty_max && penalty >= pnlty_min)
      p->penalty[typ] = ( ALIMATSCOR_t ) penalty;
    else
      errcode = ERRCODE_SWATEXCEED;
  } else {
    errcode = ERRCODE_SWATEXCEED;
  }

  return errcode;
}

short scoreGetPenalties(const ScorePenalties *p, 
			short *mismatch, short *gapinit, short *gapext)
{
  if (mismatch) *mismatch = (short) p->penalty[SCORPNLTYP_MISMATCH];
  if (gapinit) *gapinit = (short) p->penalty[SCORPNLTYP_GAPOPEN];
  if (gapext) *gapext = (short) p->penalty[SCORPNLTYP_GAPEXT];
  return (short) p->penalty[SCORPNLTYP_MATCH];
}

/******************************************************************************
 ***************************** Type ScoreMatrix *******************************
 ******************************************************************************/

struct _ScoreMatrix {     /**< Alignment scores */
  ALIMATSCOR_t gap_init;  /**< gap opening score */
  ALIMATSCOR_t gap_ext;   /**< gap extension score */
  ALIMATSCOR_t **score;   /**< score matrix for the alphabet defined in SeqCodec */
  short alphabetsiz;      /**< size of the alphabet */
};

/******************************************************************************
 ********************* Private Methods of Type ScoreMatrix ********************
 ******************************************************************************/

static int setScoreMatrix(ScoreMatrix *amp, const SeqCodec *scp,
			  const ScorePenalties *penp)
{
  short i, j;
  short alphabetsiz;
  const char *alphabet;
  short xmatch = penp->penalty[SCORPNLTYP_MISMATCH] - penp->penalty[SCORPNLTYP_MATCH];

  if (xmatch > SCORPNLTY_MAXVAL || xmatch < -1*SCORPNLTY_MAXVAL)
    return ERRCODE_ASSERT;

  if (!scp || !amp) return ERRCODE_NULLPTR;
  alphabet = seqCodecGetAlphabet(scp, &alphabetsiz);
  
  if (alphabetsiz > MAXNUM_3BIT + 1)
    return ERRCODE_ASSERT;

  amp->gap_init = penp->penalty[SCORPNLTYP_GAPOPEN];
  amp->gap_ext  = penp->penalty[SCORPNLTYP_GAPEXT];

  for (i=0; i<=MAXNUM_3BIT; i++) {
    for (j=0; j<=MAXNUM_3BIT; j++) {
      if (i>=alphabetsiz || j>=alphabetsiz ||
	  alphabet[i] == 'N' || alphabet[j] == 'N')
	amp->score[i][j] = 0;
      else if (alphabet[i] == 'X' || alphabet[j] == 'X')
	amp->score[i][j] = (ALIMATSCOR_t) xmatch;
      else if (j == i)
	amp->score[i][j] = penp->penalty[SCORPNLTYP_MATCH];
      else
	amp->score[i][j] = penp->penalty[SCORPNLTYP_MISMATCH];
    }
  }
  amp->alphabetsiz = alphabetsiz;
  return ERRCODE_SUCCESS;
}

static void fprintScoreMatrix(FILE *fp, const ScoreMatrix *amp, const SeqCodec *scp)
{
  short i, j, alphabetsiz;
  const char *alphabet;

  alphabet = seqCodecGetAlphabet(scp, &alphabetsiz);
  fprintf(fp, "== Alignment matrix ==\n   ");
  
  for (j=0; j<alphabetsiz; j++)
    fprintf(fp, "  %c ", alphabet[j]);
  fprintf(fp, "\n");
  for (i=0; i<alphabetsiz; i++) {
    fprintf(fp, " %c ", alphabet[i]);
    for (j=0; j<alphabetsiz; j++)
      fprintf(fp, "%3i ", amp->score[i][j]);
    fprintf(fp, "\n");
  }
  fprintf(fp, "Gap opening score:   %3i\n", amp->gap_init);
  fprintf(fp, "Gap extension score: %3i\n", amp->gap_ext);
  fprintf(fp, "== End of alignment matrix ==\n");
}

/******************************************************************************
 ********************** Public Methods of Type ScoreMatrix ********************
 ******************************************************************************/

ScoreMatrix *scoreCreateMatrix(const SeqCodec *scp, const ScorePenalties *penp)
{
  short i;
  ScoreMatrix *amp;

  if (!scp) return NULL;

  EMALLOCP0(amp);
  if (!amp) return NULL;
    
  ECALLOCP(MAXNUM_3BIT+1, amp->score);
  if (amp->score == NULL) {
    scoreDeleteMatrix(amp);
    return NULL;
  }
  for (i=0; i<=MAXNUM_3BIT; i++) {
    ECALLOCP(MAXNUM_3BIT+1, amp->score[i]);
    if (amp->score[i] == NULL) {
      scoreDeleteMatrix(amp);
      return NULL;
    }
  }
  
  setScoreMatrix(amp, scp, penp);

  return amp;
}

void scoreDeleteMatrix(ScoreMatrix *amp)
{
  short i;
  if ((amp)) {
    if ((amp->score)) {
      for (i=0;i<=MAXNUM_3BIT;i++)
	free(amp->score[i]);
      free(amp->score);
    }
    free(amp);
  }
}

short scoreMatrixGetAlphabetSize(const ScoreMatrix *smp)
{
  return smp->alphabetsiz;
}

#ifdef SMALT_AUTOCONF_ICC_OPT
// limit optimization by 32-bit ICC 10.x to avoid undefined
// references to __svml_exp2
#pragma optimization_level 1
#endif
double scoreMatrixCalcLambda(const ScoreMatrix *smp) {
  int i,j;
  double lambda, sum, lambda_lower, lambda_upper;

#define GETSUM for (i=sum=0; i < MINALPHABET; i++) \
                 for (j = 0; j < MINALPHABET; j++) \
                   sum += exp(lambda * smp->score[i][j]); \
               sum *= 0.0625; /* sum = sum * 1/4 * 1/4 */

  lambda_lower = 0.0;
  lambda = 0.5;
  for (;;) {
    GETSUM;
    if (sum >= 1.0) break;
    lambda_lower = lambda;
    lambda *= 2.0; 
  } 
  lambda_upper = lambda;
  while (lambda_upper - lambda_lower > .00001) {
    lambda = (lambda_lower + lambda_upper) / 2.0;
    GETSUM;
    if (sum >= 1.0) lambda_upper = lambda;
    else lambda_lower = lambda;
  } 
  return lambda;
}

ALIMATSCOR_t scoreMatrixGetMinSubstScore(const ScoreMatrix *smp)
{
  short i, j;
  const ALIMATSCOR_t *sp;
  ALIMATSCOR_t minscor = smp->score[0][0];

  for (i=0; i<=MAXNUM_3BIT; i++) {
    sp = smp->score[i];
    for (j=0; j<=MAXNUM_3BIT; j++) {
      if (minscor > sp[j])
	minscor = sp[j];
    }
  }

  return minscor;
}

ALIMATSCOR_t scoreMatrixGetAvgSubstScores(ALIMATSCOR_t *matchscor, const ScoreMatrix *smp)
{
  short i, j;
  int n_diag = 0, n_offdiag = 0;
  const ALIMATSCOR_t *sp;
  int match = 0, mismatch=0;

  for (i=0; i<MINALPHABET; i++) {
    sp = smp->score[i];
    for (j=0; j<MINALPHABET; j++) {
      if (sp[j] != 0) {
	if (i == j) {
	  match += (int) sp[j];
	  n_diag++;
	} else {
	  mismatch += (int) sp[j];
	  n_offdiag++;
	}
      }
    }
  }

  match /= n_diag;
  if (match > SCORPNLTY_MAXVAL) { 
    match = SCORPNLTY_MAXVAL;
  } else if (match < -1*SCORPNLTY_MAXVAL) {
    match = -1*SCORPNLTY_MAXVAL;
  }

  mismatch /= n_offdiag;
  if (mismatch > SCORPNLTY_MAXVAL) { 
    mismatch = SCORPNLTY_MAXVAL;
  } else if (mismatch < -1*SCORPNLTY_MAXVAL) {
    mismatch = -1*SCORPNLTY_MAXVAL;
  }
  
  if (matchscor) *matchscor = match;
  
  return mismatch;
}

short scoreGetDefaults(short *mismatch, short *gapinit, short *gapext)
{
  if (mismatch) *mismatch = DEFAULT_MISMATCH;
  if (gapinit) *gapinit = DEFAULT_GAPINIT;
  if (gapext) *gapext = DEFAULT_GAPEXT;
  return DEFAULT_MATCH;
}

int scorePrintMatrix(const ScoreMatrix *amp, const SeqCodec *scp)
{
  if (amp == NULL || scp == NULL) return ERRCODE_NULLPTR;
  fprintScoreMatrix(stdout, amp, scp);
  return ERRCODE_SUCCESS;
}

/******************************************************************************
 **************************** Type ScoreProfile *******************************
 ******************************************************************************/

struct _ScoreProfile {
  UCHAR mod;             /**< A combination of SCORE_PROFILE_MODES bit flags specifying
			  * scalar/striped modes */
  short alphabetsiz;     /**< Size of the nucleotide alphabet */
  SEQLEN_t allocsiz;    /**< Memory allocated for profile specified
			  * as equivalent to max.  sequence length
			  * that can be accomodated */
  unsigned int blocksiz; /**< Block size for memory allocation */
  SEQLEN_t length;       /**< length of sequence/profile */
  ALIMATSCOR_t **score;  /**< matrix [alphabetsiz][length+1], such
			  * that score[i][j] is the alignment score for
			  * 2-bit nucleotide code i if aligned to each
			  * position j along the sequence (starting from
			  * 0 for first nucleotide). score is for scalar 
			  * Smith-Waterman. */
#ifdef SCORE_SIMD
  void *striped_datap;    /** < for memory allocation */
  size_t striped_nalloc;  /**< memory allocated for striped profiles */

#ifdef SCORE_SIMD_IMIC
  int *striped_intp;      /**< striped profile 32-bit scores */
#else
  SIMDV_t *striped_bytep; /**< striped profile 8-bit scores */
  SIMDV_t *striped_shortp;/**< striped profile 16-bit scores */
  short bias;             /**< score bias to produce positive scores
			   * because sse2 8-bit operations are unsigned */
#endif
#endif
  ALIMATSCOR_t match_avg; /**< average match score */
  ALIMATSCOR_t mismatch_avg; /**< average mismatch score */
  ALIMATSCOR_t gap_init; /**< gap opening score */
  ALIMATSCOR_t gap_ext;  /**< gap extension score */
};

/****************************************************************************
 ******************************* Macros *************************************
 ****************************************************************************/

/******************************************************************************
 ******************* Private Methods of Type ScoreProfile *********************
 ******************************************************************************/

static int reallocScalarProfile(ScoreProfile *app, int newlen)
{
  short i;
  void *hp;
  unsigned int newsiz = ((int)((newlen+1)/app->blocksiz) + 1)*app->blocksiz;
  for (i=0; i<app->alphabetsiz; i++) {
    hp = EREALLOC(app->score[i], newsiz*sizeof(ALIMATSCOR_t));
    if (!hp) return ERRCODE_NOMEM;
    app->score[i] = (ALIMATSCOR_t *) hp;
  }

  app->allocsiz = newsiz;

  return ERRCODE_SUCCESS;
}

#ifdef score_debug
static int fprintProfile(FILE *fp, const ScoreProfile *app, int s_start, int s_length)
{
  short i;
  int j, s_end;
  
  if (s_start > app->length) return ERRCODE_ARGRANGE;
  if (s_length > app->length - s_start || s_length < 1) {
    s_length = app->length - s_start;
    if (!s_start) s_length += 1;
  }
  
  s_end = s_start + s_length;
  for (i=0;i<app->alphabetsiz;i++) {
    fprintf(fp, "%i: ", i);
    for (j=s_start; j<s_end; j++)
      fprintf(fp, "%2i ", app->score[i][j]);
    fprintf(fp,"\n");
  }
  return ERRCODE_SUCCESS;
}
#endif

#ifdef SCORE_SIMD
static int makeStripedProfileFromSequence(ScoreProfile *app, const char *seq_basp, 
					  SEQLEN_t length, const ScoreMatrix *amp)
{
  size_t len, n_alloc_striped;

#define FILL_STRIPED(cond, typ, nvelem, alimemp, bias)			\
  if ((cond)) {								\
    short i;								\
    SEQLEN_t segsiz, segbyt;						\
    typ *sprofp;							\
    sprofp = (typ *) (alimemp);						\
    segsiz = (length + (nvelem) - 1)/(nvelem);				\
    segbyt = segsiz * nvelem;						\
    for (i=0; i<app->alphabetsiz; i++) {				\
      SEQLEN_t j;							\
      const ALIMATSCOR_t *sp = amp->score[i];				\
      for (j=0; j<segsiz; j++) {					\
	SEQLEN_t k;							\
	for (k=j; k<length; k += segsiz)				\
	  *sprofp++ = (typ) sp[seq_basp[k]&SEQCOD_ALPHA_MASK] - (bias);	\
	for (; k<segbyt; k += segsiz)					\
	  *sprofp++ = 0;						\
      }									\
    }									\
  }

#ifdef SCORE_SIMD_IMIC
  app->striped_intp = 0;
  len = (app->mod & SCORPROF_STRIPED_32)? (length + SCORSIMD_NINTS - 1) / SCORSIMD_NINTS: 0;
#else  
  size_t lenBYT = (app->mod & SCORPROF_STRIPED_8)? (length + SCORSIMD_NBYTES - 1) / SCORSIMD_NBYTES: 0;
  app->striped_bytep = app->striped_shortp = 0;
  len = lenBYT + 
    ((app->mod & SCORPROF_STRIPED_16)? (length + SCORSIMD_NSHORTS - 1) / SCORSIMD_NSHORTS: 0);
#endif

  /* reallocate memory if necessary */
  n_alloc_striped = 1 /* slack for memory alignment */
    + len * app->alphabetsiz; 
  
#ifdef SCORE_SIMD_SSE2
  /* additional slack (for two alignment operations) */
  if ((app->mod & SCORPROF_STRIPED_8) && (app->mod & SCORPROF_STRIPED_16))
    n_alloc_striped += 1;
#endif

  n_alloc_striped *= 
#ifdef SCORE_SIMD_IMIC
    SCORSIMD_NINTS*sizeof(int);
#else
    sizeof(SIMDV_t);
#endif

  if (n_alloc_striped > app->striped_nalloc) {
    void *hp;
    n_alloc_striped = (n_alloc_striped + app->blocksiz - 1)/app->blocksiz;
    n_alloc_striped *= app->blocksiz;
    if (NULL == app->striped_datap) {
      hp = EMALLOC(n_alloc_striped);
    } else {
      hp = EREALLOCP(app->striped_datap, n_alloc_striped);
    }
    if (!hp) 
      return ERRCODE_NOMEM;
    app->striped_datap = hp;
    app->striped_nalloc = n_alloc_striped;
  }

  /* align memory to 16/64 byte boundaries */
#ifdef SCORE_SIMD_IMIC
  if (app->mod & SCORPROF_STRIPED_32) {
    app->striped_intp = (int *) SCORE_ALIGN_MEMORY(app->striped_datap);
  }
#else
  if (app->mod & SCORPROF_STRIPED_8) {
    app->striped_bytep = (SIMDV_t *) SCORE_ALIGN_MEMORY(app->striped_datap);
    app->bias = scoreMatrixGetMinSubstScore(amp);
    if (app->mod & SCORPROF_STRIPED_16)
      app->striped_shortp = app->striped_bytep + lenBYT * app->alphabetsiz;
  } else if (app->mod & SCORPROF_STRIPED_16) {
    app->striped_shortp = (SIMDV_t *) SCORE_ALIGN_MEMORY(app->striped_datap);
  }
#endif

  /* fill striped profile */
#ifdef SCORE_SIMD_IMIC
  FILL_STRIPED(app->mod & SCORPROF_STRIPED_32, 
	       int, SCORSIMD_NINTS,
	       app->striped_intp, 0);
#else
  FILL_STRIPED(app->mod & SCORPROF_STRIPED_8, 
	       UCHAR, SCORSIMD_NBYTES,
	       app->striped_bytep, app->bias);

  FILL_STRIPED(app->mod & SCORPROF_STRIPED_16, 
	       short, SCORSIMD_NSHORTS,
	       app->striped_shortp, 0);
#endif
  return ERRCODE_SUCCESS;
}
#endif //ifdef SCORE_SIMD

/******************************************************************************
 ******************** Public Methods of Type ScoreProfile *********************
 ******************************************************************************/

ScoreProfile *scoreCreateProfile(int blocksize, const SeqCodec *codep, UCHAR mod)
{
  short i;
  short alphabetsiz;
  ScoreProfile *app;

  EMALLOCP0(app);
  if (!app) return NULL;

  seqCodecGetAlphabet(codep, &alphabetsiz);
  if (blocksize < 1) blocksize = DEFAULT_BLOCKSIZ_PROFILE;

#ifdef SCORE_SIMD
  if (!(mod & (SCORPROF_SCALAR |
#ifdef SCORE_SIMD_SSE2
	       SCORPROF_STRIPED_8 | SCORPROF_STRIPED_16
#else
	       SCORPROF_STRIPED_32
#endif
	       )))
    return NULL;
#else
  mod = SCORPROF_SCALAR;
#endif
  if (mod & SCORPROF_SCALAR) {
    ECALLOCP(alphabetsiz, app->score);
    if (!app->score) {
      scoreDeleteProfile(app);
      return NULL;
    }
    for (i=0; i<alphabetsiz; i++) {
      ECALLOCP(blocksize, app->score[i]);
      if (!app->score[i]) {
	scoreDeleteProfile(app);
	return NULL;
      }
    }
  }
#ifdef SCORE_SIMD
    if ((mod & (
#ifdef SCORE_SIMD_SSE2
		 SCORPROF_STRIPED_8 | SCORPROF_STRIPED_16
#else
		 SCORPROF_STRIPED_32
#endif
		))) {
    ECALLOCP(blocksize, app->striped_datap);
    if (!app->striped_datap) {
      scoreDeleteProfile(app);
      return NULL;
    }
    app->striped_nalloc = blocksize;
  }
#endif
  app->allocsiz = blocksize;
  app->blocksiz = blocksize;
  app->alphabetsiz = alphabetsiz;
  app->mod = mod;
  return app;
}

void scoreDeleteProfile(ScoreProfile *app)
{
  short i;
  if (app) {
    if (app->score) {
      for (i=0;i<app->alphabetsiz; i++)
	free(app->score[i]);
      free(app->score);
    }
#ifdef SCORE_SIMD
    free(app->striped_datap);
#endif
    free(app);
  }
}
    
int scoreMakeProfileFromSequence(ScoreProfile *app, const SeqFastq *sqp, 
				 const ScoreMatrix *amp)
{
  int errcode;
  const char *cp, *seq_basp;
  char cod;
  short i;
  SEQLEN_t j, length;
  ALIMATSCOR_t *sc, *hp;

  seq_basp = seqFastqGetConstSequence(sqp, &length, &cod);
  if (cod != SEQCOD_MANGLED) return ERRCODE_SEQCODE;

  if (app->mod & SCORPROF_SCALAR) {
    if (length > app->allocsiz-1 && 
      (errcode = reallocScalarProfile(app, length)))
      return errcode;

    for (i=0; i<app->alphabetsiz;i++) {
      cp = seq_basp;
      hp = app->score[i];
      sc = amp->score[i];
      for (j=length; *cp && j>0; j--)
	*hp++ = sc[(*cp++)&SEQCOD_ALPHA_MASK];
      *hp = 0;
    }
  }
#ifdef SCORE_SIMD
  if ((app->mod & (
#ifdef SCORE_SIMD_SSE2
		   SCORPROF_STRIPED_8 | SCORPROF_STRIPED_16
#else
		   SCORPROF_STRIPED_32
#endif
		   )) && 
      (errcode = makeStripedProfileFromSequence(app, seq_basp, length, amp)))
    return errcode;
#endif

  app->length = length;
  app->mismatch_avg = scoreMatrixGetAvgSubstScores(&app->match_avg, amp);
  app->gap_init = amp->gap_init;
  app->gap_ext = amp->gap_ext;

#ifdef score_debug
  fprintScoreProfile(stdout, app, 0, 0);
  fprintf(stdout, "gap_init: %i\n", app->gap_init);
  fprintf(stdout, "gap_ext: %i\n", app->gap_ext);
#endif

  return ERRCODE_SUCCESS;
}

ALIMATSCOR_t * const *scoreGetProfile(short *alphabetsiz, SEQLEN_t *seqlen, 
				      ALIMATSCOR_t *gap_init, ALIMATSCOR_t *gap_ext,  
				      const ScoreProfile *spp)
{
  if (alphabetsiz) *alphabetsiz = spp->alphabetsiz;
  if (seqlen) *seqlen = spp->length;
  if (gap_init) *gap_init = -1 * spp->gap_init;
  if (gap_ext) *gap_ext = - 1 * spp->gap_ext;

  return spp->score;
}

short scoreProfileGetAvgPenalties(short *mismatch_avg,
				  short *gap_init, short *gap_ext,
				  const ScoreProfile *spp)
{
  if (mismatch_avg) *mismatch_avg = (short) spp->mismatch_avg;
  if (gap_init) *gap_init = (short) spp->gap_init;
  if (gap_ext) *gap_ext = (short) spp->gap_ext;

  return (short) spp->match_avg;
}


#ifdef SCORE_SIMD
const void *scoreGetStripedProfile(short *alphabetsiz, SEQLEN_t *seqlen, 
				   unsigned short *gap_init, unsigned short *gap_ext,
				   unsigned short *bias, int *segsiz, char mod,
				   const ScoreProfile *spp)
{
  const void *p = 0;

  if (alphabetsiz) *alphabetsiz = spp->alphabetsiz;
  if (seqlen) *seqlen = spp->length;
  if (gap_init) *gap_init = (unsigned short) ((spp->gap_init < 0)? -1: 1) * spp->gap_init;
  if (gap_ext) *gap_ext = (unsigned short) ((spp->gap_ext < 0)? -1: 1) * spp->gap_ext; 
  if (bias) *bias = (unsigned short) 
#ifdef SCORE_SIMD_IMIC
	      0;
#else
	      (spp->bias < 0)? -1*spp->bias: 0;
#endif

#ifdef SCORE_SIMD_IMIC
  if (mod  == SCORPROF_STRIPED_32) {
    if (segsiz) 
      *segsiz = (spp->length + SCORSIMD_NINTS - 1) / SCORSIMD_NINTS;
    p = spp->striped_intp; 
  }
#else
  if (mod  == SCORPROF_STRIPED_8) {
    if (segsiz) 
      *segsiz = (spp->length + SCORSIMD_NBYTES - 1) / SCORSIMD_NBYTES;
    p = spp->striped_bytep;
  }
  else if (mod == SCORPROF_STRIPED_16) {
    if (segsiz) 
      *segsiz = (spp->length + SCORSIMD_NSHORTS - 1) / SCORSIMD_NSHORTS;
    p = spp->striped_shortp;
  }
#endif

  return p;
}

#endif
