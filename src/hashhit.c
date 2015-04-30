/** Gathering hits from a hash index */

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
#include <string.h>
#include <limits.h>
#include <math.h>

#include "elib.h"
#include "sort.h"
#include "hashhit.h"

enum {
  FALSE = 0,
  TRUE = 1,                /**< for booleans */
  HASH_QMASKBLOCK = 512,   /**< allocation block for query mask */
  NREPEATS = 4,            /* filter out kmers that repeat within NREPEATS bases */
  MINHIT_PER_TUPLE = 16,   /* value down to which one can adjust mincut */
  DEFAULT_MAXHIT_PER_TUPLE = 10000,/* Default maximum number of hits a seed may have */
  HITLST_MINSIZ = 8192,    /* target length of a hit list */
  HITLST_DEFBLKSZ = 16384, /**< Default block size for hit list memory reallocation */
  //HITLST_QLENSIZ_FACT = 80, /**< Scales the size of the hit list with the read length */
  HITLST_LOGQLENSIZ_FACT = 32, /**< Scales the size of the hit list with nlogn of read length */
  HITLST_MAXSIZ = INT32_MAX, /**< Maximum size of hit list (has to fit a signed int */
  KTUPLIST_DEFBLKSZ = 512, /**< Memory allocation block for type KTUPLST */
  DEFAULT_FILTER_BLOCKSIZ = 1024, /**< Default block size for HashHitFilter */
  MASK32BIT = 0xFFFFFFFF,
  HITINFO_MAXCOVER_PERCENT = 80, /**< Stop collecting seeds if this percentage of the read is covered */
  HITINFO_MINSEEDNUM = 3,        /** Minimum number of seeds */
  HITINFO_MINCOVER_KMER = 2,     /** The rare k-mer words of the read must cover at least
				  * HITINFO_MINCOVER_KMER * kmer_word_length + nskip bases */
#ifdef hashhit_minimise_coverdeficit
  HITINFO_MINCOVER_NKTUP = 3,
  NBITS_PER_BYTE = 8,
#endif
  /* HITINFO_SEED_PERCENTILE = 80, */ /**< use only the least frequent seeds, as percentage of total */  
  NBIT_INT32 = 31,

#ifdef hashhit_debug
  NBITS_HASHWORD = 64,
#endif
#ifdef hashhit_basqual
  NBIT_KMER_NHIT = 24,     /**< Number of bits encoding the number of hits of a k-mer,
			    * don't change that */
  KMER_NHIT_MASK = 0x00ffffff, /** mask set this to (1<<NBIT_KMER_NHIT) - 1; */
  QCLASS_SCALFAC = 10,
  QCLASS_NUM = 4,          /**< number of k-mer classes */ 
#endif
#ifdef hashhit_dump_sortarray
  MIN_ARRLEN = 50000,
#endif
};
#ifdef hashhit_minimise_coverdeficit
enum FRAMECOV_STATUS_FLAGS {
  FRAMECOV_INACTIVE = 0x01, /** used to flag out certain frames */
};
#endif

enum HITINFO_STATUS_FLAGS {
  HITINFO_REVERSE = 0x01,  /**< If set, hit info is for reverse strand */
  HITINFO_SORTED  = 0x02,  /**< If set, kmer words have been sorted by hit frequency
			    * in ascending order */
  HITINFO_RANK    = 0x04,  /**< If set, a rank has been determined in the sorted kmer words
			    * above which k-mer words may not be used because the hit list
			    * gets too long */
  HITINFO_FRAME   = 0x08,  /**< frames FRAMECOV have been set up */
};

enum HASHHIT_STATUS_FLAGS { /**< status of type HashHitList */
  HASHHITSTAT_REVERSE = 0x01, /**< List is for reverse complement */
  HASHHITSTAT_SHIFTS = 0x02,  /**< the array \<HashHitList\>.sqdat contains shifts, i.e. differences
			       * of base offsets between subject and query sequence. The list is
			       * sorted by shifts and query offset */
};

typedef int64_t TDRFILTR;
typedef uint8_t BOOL;
typedef uint8_t UCHAR;
typedef SEQLEN_t COVERAGE;
typedef uint32_t QMASK;
typedef uint32_t SORTKEY_t;
typedef uint64_t HASHPACK_t;

/** The following structure contains k-tupe hits for one query
 * sequence either in forward direction (is_reverse == 0) or for the
 * reverse complement (is_reverse == 1).  
 *
 * When the list is fully set up (status == HASHHIT_SHIFTS) the array
 * sqdat contains k-tuple base offsets in the query sequence (32 bits)
 * and the 'shift' (32 bits) merged into a 64-bit unsigned integer.
 * 
 * The 'shift' shift = s_pos - q_pos is the position of the hit in the
 * subject sequence s_pos relative to the query sequence q_pos. So
 * that a pair of ktuples of the same shift define between them an
 * alignment where the number of insertions equals the number of
 * deletions.
 *
 * For the reverse complement hits (is_reverse == 1) q_pos is the
 * position of the first base of the forward (F-)k-tuple from which
 * the reverse complement RC-k-tuple was generated. s_pos is the first
 * base of the hit of the RC-k-tuple.  I.e. for a k-tuple
 * (b[i],b[i+1], ...., b[i+kmer-1]) with b[i] the base at position i
 * in the query sequence, the reverse complement is (B[i+kmer-1],
 * B[i+kmer-2], ..., B[i]) where B[j] is the complement of base
 * b[j]. For a hit of the reverse complement in the subject sequence
 * (S[j], S[j+1], ..., S[j+kmer-1]) q_pos_rc = i, s_pos_rc = j and
 * shift_rc = q_pos_rc + s_pos_rc;
 * 
 * to make the number stored as shift >= 0, shift' = s_pos - q_pos + q_len - 1
 * for the forward case (is_reverse == 0)
 */
#ifdef hashhit_minimise_coverdeficit
typedef struct _FRAMECOV { /**< contains stuff for calculating the coverage in a frame */
  UCHAR status;            /**< combination of FRAMECOV_STATUS flags */
  SEQLEN_t *sxp;             /**< array of seed indices for this frame */
  SEQLEN_t ns;               /**< number of seeds (size of array sxp) */
  SEQLEN_t ns_used;           /**< maximum index to be used */
  COVERAGE cover;            /**< number of bases covered by k-mer words in this frame */
  QMASK *bufp;               /**< buffer for caclulating coverage in the frame */
} FRAMECOV;
#endif

typedef struct _SEED {
  HASHNUM_t posidx;     /**< Index of the block of word positions */
  HASHNUM_t nhits;      /**< Number of hits */
#ifdef hashit_basqual
  int qclass;
#endif  
#ifdef hashhit_debug
  HASHWORD_t kmer_word;
#endif
  HASHNUM_t  cix;       /**< current offset in ktuples-index. This is used
		      * to feed the seeds to the hit list
		      * sequence-by-sequence remembering the offset of
		      * the last sequence */
  SEQLEN_t qoffs;     /**< position in query sequence */
} SEED;

struct _HashHitInfo {
  UCHAR status;    /**< Combination of binary flags from HITINFO_STATUS_FLAGS */
  SEQLEN_t qlen;   /**< length of query sequence */
  UCHAR *qmaskp;   /**< array of size qlen contains for each k-tuple a
		    * qualifier indicating its fate, in particular
		    * whether or not the k-tuple registered in the hit
		    * list (HASH_HIT_QUALIFIERS)*/
  UCHAR *qbufp;     /**< array of size qlen used as buffer for calculating k-mer word cover */
  UCHAR ktup;      /**< hashed word length */
  UCHAR nskip;     /**< Sampling step for kmer words */
  SEQLEN_t n_seeds;/**< length of the arrays (== length of query
		    * sequence - ktup + 1, unless there are Ns) */
  SEED *seedp;     /**< array of size n_seeds */
  SEQLEN_t *sidxp;   /**< array of size n_seeds, seed indices
		    * for sorting */
  SORTKEY_t *nhitqual_sortkeyp; 
  /**< Array of size n_seeds. The lower NBIT_KMER_NHIT bits of
   * nhitqual_sortkeyp[i] represent the number of hits of the
   * k-mer. If hashhit_basqual is defined, the upper NBIT_UINT32 -
   * NBIT_KMER_NHIT bits encode the quality class of the k-mer. This
   * is used as key for sorting by hit-frequency and base quality
   * simultaneously */

  size_t n_alloc;  /**< number of array elements (qmask, qpos, nhitqual_sortkeyp,
		    * kip) for which memory is allocated */
  int blksz;

#ifdef hashhit_minimise_coverdeficit
  COVERAGE cover_deficit; /**< max difference in cover over the frames */
  FRAMECOV *frcp;     /**< array of size nskip */
  size_t frcbufsz;  /**< Size of buffer frcp[i].bufp */
#else
  uint32_t *coverp;    /** Array of size nskip. coverp[o] is the
		     * coverage achieved by a run of k-mer words with
		     * spacings modulo nskip = 0 where the run has
		     * offset o - can also contain counts dependent on
		     * whether hashCalcHitInfoCover() or
		     * hashCalcHitInfoNumberOfHits(,0) was called
		     * last */
  uint32_t **framep;    /**< Array of size nskip containing for each frame an array of indices of sidxp */
  uint32_t *countp;   /**< Array of size nskip containing the seed counts for each frame 
		     *  this is only used in some applications */
#endif
  uint32_t seed_rank; /**< if > 0 only use the first seed_rank seeds in the 
		     * (sorted) seed array.
		     * seed_rank is determined such that the sum of all seed hits does
		     * not exceed the cut-off nhit_total_cutoff. In addition, 
		     * seed_rank <= n_seeds*seed_percentile/100 as long as the coverage of the seeds
		     * in each frame is at least seed_mincover */
};

struct _HashHitList { /**<contains k-tuple hits for one query sequence */
  UCHAR status;     /**< combination of HASHHIT_STATUS flags */
  int nhits;        /**< total number of hits */
  int nhits_max;    /**< maximum number of hits for query */
  int nhits_alloc;  /**< Memory for this many hits has been allocated */
  int nhits_blksz;  /**< Block size for memory re-allocation */
  HASHPACK_t *sqdat;/**< [nhits] combined shift and query data for k-tuple hit.
		  * for status & HASHHIT_SEQASS == 0: 
		  *      upper 32 bits: ktuple number in subject sequence 
		  *      lower 32 bits: ktuple number == base offset in query sequence
		  * for status & HASHHIT_SHIFTS != 0:
		  *      upper 32 bits: shift
		  *      lower 32 bits: base offset in query sequence
		  */
  UCHAR ktup;    /**< k-tuple length */
  UCHAR nskip;   /**< skip step size */
  SEQLEN_t qlen; /**< query length */
  char *qmask;   /**< contains for each k-tuple a qualifier indicating its fate, in particular
		  * whether or not the k-tuple is registered in the hit list (HASH_HIT_QUALIFIERS)
		  */
  size_t qmask_alloc; /**< currently allocated memory for qmask */
};

typedef struct _FILTERIVAL { /**< Distance interval for filtering hits */
  HASHPOS_t upper;
  HASHPOS_t lower;
} FILTERIVAL;

struct _HashHitFilter { /**< Set of intervals for filtering hit lists */
  FILTERIVAL *ivp;      /**< array of filter intervals */
  short num;            /**< number of filter intervals */ 
  short n_alloc;        /**< allocated memory as number of intervals */
  short blocksiz;
};

/******************************************************************************
 *********************************** Macros ***********************************
 ******************************************************************************/

#define MAKE_NEXT_WORD(ktup_word, data)	\
  if ((is_reverse)) {							\
    (ktup_word) = ((ktup_word)>>2) + (((HASHWORD_t)(((data)^SEQCOD_STDNT_MASK) & SEQCOD_STDNT_MASK))<<(ktup_rc_addpos)); \
  } else {								\
    (ktup_word) = ((ktup_word)<<2) + ((data) & SEQCOD_STDNT_MASK);	\
  }									

#ifdef hashhit_basqual
#define GET_LOWEST_QUAL(minqual_pos, minqual_num, first_pos, last_pos) \
  if ((minqual_pos) < (first_pos)) {					\
    if (minqual_num == 0) return ERRCODE_ASSERT;			\
    SEQLEN_t t;								\
    (minqual_num)--;                                                    \
    if ((minqual_num) > 0) {						\
      for (t = (minqual_pos) + 1; t <= (last_pos) && qualp[t] > qualp[(minqual_pos)]; t++); \
      (minqual_pos) = t;}						\
    else {								\
      (minqual_pos) = (first_pos);					\
      minqual_num = 1;							\
      for (t = (first_pos) + 1; t <= (last_pos); t++) {			\
	if (qualp[t] < qualp[(minqual_pos)]) {(minqual_pos) = t; (minqual_num) = 1;} \
	else if (qualp[t] == qualp[(minqual_pos)]) {(minqual_num)++;}}}} \
  else if (qualp[(last_pos)] < qualp[(minqual_pos)]) {			\
    (minqual_pos) = (last_pos);						\
    (minqual_num) = 1;}							\
  else if (qualp[(last_pos)] == qualp[(minqual_pos)]) {			\
    (minqual_num)++;} 
#endif

#define SET_NEXT_SHIFT(is_reverse, sqdat, pos, q, nskip, offbit)	\
  if ((is_reverse)) {							\
    (sqdat) = ((((uint64_t) (pos)) + (q)/(nskip))<<HASHHIT_HALFBIT) + (q); \
  } else {								\
    (sqdat) = (((((uint64_t) (pos))|offbit) - (q)/(nskip))<<HASHHIT_HALFBIT) + (q); \
  }
    /* if ((is_reverse)) {							\ */
    /* 	  sqdatp[i] = (((uint64_t)(posp[i])+tuplectr/nskip)<<HASHHIT_HALFBIT) + tuplectr; \ */
    /* 	} else {\ */
    /* 	  sqdatp[i] = (((((uint64_t)(posp[i]))|offbit) - tuplectr/nskip)<<HASHHIT_HALFBIT) + tuplectr;\ */
    /* 	} */

#define SET_NEXT_32BIT_SHIFT()\
  	if ((is_reverse)) {\
	  sqdatp[n] = (((uint64_t) ((int) ((posp[n]-segpos_lo)*nskip))+qo)<<HASHHIT_SEMIBIT) + qo;\
	} else {\
	  sqdatp[n] = (((uint64_t) ((int) ((posp[n]-segpos_lo)*nskip))-qo)<<HASHHIT_SEMIBIT) + qo;\
	}

#define CALC_FRAM_BLKSZ(rlen, nskip) (((size_t) ((rlen) + (nskip)))/(nskip))

#define CALC_FRAM_BUFSZ(rlen, nskip)\
  (((size_t) (CALC_FRAM_BLKSZ(rlen,nskip) + sizeof(QMASK)))/(sizeof(QMASK)*NBITS_PER_BYTE))

/******************************************************************************
 ******************************* Private Methods ******************************
 ******************************************************************************/
#ifdef hashhit_debug
static char *writeWord(char *buff, HASHWORD_t key, UCHAR ktup)
{
  const char alphabet[]="ACGT";
  UCHAR i;

  for (i=ktup; i>0; i--) {
    buff[i-1]= alphabet[key & SEQCOD_STDNT_MASK];
    key = key>>2;
  }
  buff[ktup] = '\0';
  return buff;
}
#endif

static BOOL checkForRepeats(TDRFILTR word, TDRFILTR tdrf[])
     /**< Return TRUE the key occurred the last NREPEATS times. 
      *   FALSE else. */
{
  BOOL isRepeat = FALSE;
  UCHAR i;
  
  for (i=0; i<NREPEATS; i++)
    if (word == tdrf[i]) {
      isRepeat = TRUE;
      break;
    }
  memmove(tdrf+1, tdrf, (NREPEATS-1)*sizeof(TDRFILTR));
  tdrf[0] = word;
  return isRepeat;
}

static void initRepeatFilter(TDRFILTR tdrf[])
{
  UCHAR i;
  for (i=0; i<NREPEATS; i++) tdrf[i] = -i-1;
}

/******************************************************************************
 ******************** Private Methods of Type HashHitInfo *********************
 ******************************************************************************/
#ifdef hashhit_minimise_coverdeficit

static void initHitInfoFrames(const HashHitInfo *hip)
/**< Preserve the order of elements in hip->sidxp (may
 * have beed sorted e.g. by increasing number of hits in genome).
 */
{
  UCHAR f;
  uint32_t i;

  for (f=0; f<hip->nskip; f++) {
    FRAMECOV *frp =  hip->frcp + f;
    frp->status = 0;
    frp->ns = 0;
    frp->ns_used = 0;
    frp->cover = 0;
    memset(frp->bufp, 0, hip->frcbufsz*sizeof(QMASK));
  }
  
  for (i=0; i<hip->n_seeds; i++) {
    SEQLEN_t ix = hip->sidxp[i];
    f = (UCHAR) hip->seedp[ix].qoffs%(hip->nskip);
    hip->frcp[f].sxp[hip->frcp[f].ns++] = ix;
  }

  return;
}

static void
resetHitInfoFrames(const HashHitInfo *hip)
/**< Reset FRAMECOV.ns_used to 0 */
{
  UCHAR f;
  for (f=0; f<hip->nskip; f++) {
    hip->frcp[f].status = 0;
    hip->frcp[f].ns_used = 0;
  }
}
#endif //#ifdef hashhit_minimise_coverdeficit

#ifdef hashhit_debug
static void fprintHitInfo(FILE *fp, UCHAR ktup, const HashHitInfo *p)
{
  uint32_t i;
  char keybuf[NBITS_HASHWORD + 1];
  const SEED *sp;

  fprintf(fp, "=================== HashHitInfo ========================\n");
  fprintf(fp, "reverse: %1i\n", (int) (p->status & HITINFO_REVERSE));
  fprintf(fp, "query length: %u\n", p->qlen);
  fprintf(fp, "number of k-tuples: %u\n", p->n_seeds);
  for (i=0; i<p->n_seeds; i++) {
    sp = p->seedp + p->sidxp[i];
    fprintf(fp, "[%u] %s qo = %u nh = %i\n", i, writeWord(keybuf, sp->kmer_word, ktup),
	    sp->qoffs, 
#ifdef hashhit_basqual
	    p->nhitqual_sortkeyp[i]&KMER_NHIT_MASK
#else
	    p->nhitqual_sortkeyp[i]
#endif
	    );
  }
  fprintf(fp, "================End of HashHitInfo =====================\n");
}
#endif


static int reallocHashHitInfo(HashHitInfo *p, uint32_t newlen)
{
  void *hp;
  UCHAR o;
  size_t siz;
#ifndef hashhit_minimise_coverdeficit
  size_t dblk;
#endif

  if (newlen < 1) return ERRCODE_ARGRANGE;

  siz = (newlen + p->blksz - 1)/p->blksz;
  siz *= p->blksz;

  hp = EREALLOCP(p->qmaskp, siz);
  if (!hp) return ERRCODE_NOMEM;
  p->qmaskp = hp;
  
  hp = EREALLOCP(p->qbufp, siz);
  if (!hp) return ERRCODE_NOMEM;
  p->qbufp = hp;
  
  hp = EREALLOCP(p->sidxp, siz);
  if (!hp) return ERRCODE_NOMEM;
  p->sidxp = hp;
  
  hp = EREALLOCP(p->nhitqual_sortkeyp, siz);
  if (!hp) return ERRCODE_NOMEM;
  p->nhitqual_sortkeyp = hp;

  hp = EREALLOCP(p->seedp, siz);
  if (!hp) return ERRCODE_NOMEM;
  p->seedp = hp;
#ifdef hashhit_minimise_coverdeficit
  if (p->frcp) {
    size_t sz = CALC_FRAM_BLKSZ(siz, p->nskip);
    size_t sq = CALC_FRAM_BUFSZ(siz, p->nskip);
    p->frcbufsz = (uint32_t) sq;

    hp = EREALLOCP(p->frcp->sxp, sz*p->nskip);
    if (!hp) return ERRCODE_NOMEM;
    p->frcp->sxp = hp;
    for (o=1; o<p->nskip;o++) p->frcp[o].sxp = p->frcp[0].sxp + o*sz;
   
    hp = EREALLOCP(p->frcp->bufp, sq*p->nskip);
    if (!hp) return ERRCODE_NOMEM;
    p->frcp->bufp = hp;
    for (o=1; o<p->nskip;o++) p->frcp[o].bufp = p->frcp[0].bufp + o*sq;
  }
#else
  hp = EREALLOCP(p->framep[0], siz + p->nskip);
  if (!hp) return ERRCODE_NOMEM;
  p->framep[0] = hp;
  dblk = siz/p->nskip + 1;
  for (o=1; o<p->nskip; o++) {
    p->framep[o] = p->framep[0] + o*dblk;
  }
#endif
  p->n_alloc = siz;
  return ERRCODE_SUCCESS;
}

static int collectHitInfo(HashHitInfo *hip,
			  BOOL is_reverse,
			  HASHNUM_t maxhit_per_tuple,
			  UCHAR basqual_threshold,
			  SEQLEN_t seq_start,
			  SEQLEN_t seq_end,
			  const SeqFastq *seqp,
			  const HashTable *htp)
     /**< Fill HitInfo with ktuple hits in sequence seqp
      * If maxhit_per_ tuple > 0, only count seeds with hits <= maxhit_per_tuple.
      * if seq_end >= seq_start + ktup, only collect hits for that segment in seqp.
      */
{
  int errcode = ERRCODE_SUCCESS;
  char codtyp;
  UCHAR nskip, non_stdnt = 0;
  SEQLEN_t s, seqlen, seedctr, tuplectr;
  const char *datap = seqFastqGetConstSequence(seqp, &seqlen, &codtyp);
  UCHAR *qmaskp;
  const UCHAR *qualp = (const UCHAR *) seqFastqGetConstQualityFactors(seqp, NULL, NULL);
  const UCHAR ktup = hashTableGetKtupLen(htp, &nskip);
  const UCHAR minqval = basqual_threshold + SEQCOD_QVAL_OFFS;
  const UCHAR ktup_rc_addpos = (ktup-1)<<1;
  HASHWORD_t ktup_word = 0LL;
  const HASHWORD_t wordmask = (((HASHWORD_t) 1)<<(ktup<<1))-1;
  SEED *sp;
  HASHNUM_t posidx, nhits;
  TDRFILTR tr_filter[NREPEATS];
#ifdef hashhit_basqual
  int qclass;
  UCHAR  minqual_num;
  SEQLEN_t minqual_pos;
#endif
  hip->status = 0;

  if (ktup > 31 || ktup != hip->ktup || nskip != hip->nskip) 
    return ERRCODE_ASSERT;

  if (((short) basqual_threshold) + SEQCOD_QVAL_OFFS > UCHAR_MAX)
    return ERRCODE_QUALVAL;
 
  if (codtyp != SEQCOD_MANGLED) 
    return ERRCODE_SEQCODE;
  
  if (seqlen < ktup) 
    return ERRCODE_SHORTSEQ;

  /* qualp might be NULL for FASTA input */
  /*   if (qualp == NULL) */
  /*     return ERRCODE_ASSERT; */

  if (seqlen >= hip->n_alloc &&
      (errcode = reallocHashHitInfo(hip, seqlen+1)))
    return errcode;
  if ((is_reverse)) 
    hip->status |= HITINFO_REVERSE;
  hip->qlen = seqlen;

  if (seq_end >= seqlen)
    seq_end = seqlen - 1;
  if (seq_end < (seq_start + ktup - 1)) {
    seq_start = 0;
    seq_end = seqlen - 1;
  }

  qmaskp =  hip->qmaskp;
  for(s = 0; s < seq_start; s++)
    qmaskp[s] = HITQUAL_NOHIT;

  initRepeatFilter(tr_filter);

  seedctr = hip->n_seeds = 0;  

#ifdef hashhit_basqual
  minqual_pos = 0;
  minqual_num = 0;
#endif
  tuplectr = s = seq_start;
  for (; s < seq_start+ktup-1; s++) {
    if ((datap[s]&SEQCOD_STDNT_TESTBIT) || /* do not look up non-standard nt */
	(qualp != NULL && qualp[s] < minqval)) /* do not look up kmer word with low base quality values */
      non_stdnt = ktup;
    else if ((non_stdnt))
      non_stdnt--;
    MAKE_NEXT_WORD(ktup_word, datap[s]);
#ifdef hashhit_basqual
    if (qualp == NULL || qualp[s] < qualp[minqual_pos]) {
	minqual_pos = s;
	minqual_num = 1;
      } else if (qualp[s] == qualp[minqual_pos]) {
	minqual_num++;
      }
    }     
#endif
  }
  
  for(; s <= seq_end; tuplectr++,s++) {
    if ((datap[s]&SEQCOD_STDNT_TESTBIT) || /* do not look up non-standard nt */
	(qualp != NULL && qualp[s] < minqval)) /* do not look up kmer word with low base quality values */
      non_stdnt = ktup;
    else if ((non_stdnt))
      non_stdnt--;
    MAKE_NEXT_WORD(ktup_word, datap[s]);
 
    if ((non_stdnt)) {
      qmaskp[tuplectr] = HITQUAL_NONSTDNT;
      continue;
    }
#ifdef hashhit_basqual
    if (qualp != NULL) {
      GET_LOWEST_QUAL(minqual_pos, minqual_num, tuplectr, s);
    }
#endif
    if (checkForRepeats((TDRFILTR) (ktup_word&wordmask), tr_filter)) {
      qmaskp[tuplectr] = HITQUAL_REPEAT;
      continue;
    }
    
    nhits = hashTableGetKtupleHits(NULL, &posidx, htp, ktup_word);
    
    if (nhits < 1) {
      qmaskp[tuplectr] = HITQUAL_NOHIT;
      continue;
    }

#ifdef hashhit_basqual
    if (nhits > KMER_NHIT_MASK ||
	((maxhit_per_tuple > 0) && nhits > maxhit_per_tuple)) {
      qmaskp[tuplectr] = HITQUAL_MULTIHIT;
      continue;
    }

    if (minqual_num < 1)
      return ERRCODE_ASSERT;

    if (qualp != NULL) {
      qclass = qualp[minqual_pos] - log10((double) minqual_num);
      if (qclass < 0) qclass = 0;
      else {
	qclass = QCLASS_SCALFAC/qclass + 1;
	if (qclass > QCLASS_NUM)
	  qclass = QCLASS_NUM;
      }
      qclass = QCLASS_NUM - qclass;
    } else {
      qclass = 0;
    }
    hip->nhitqual_sortkeyp[seedctr] = (((uint32_t) qclass)<<NBIT_KMER_NHIT) + nhits;
#else
    if (maxhit_per_tuple > 0 && nhits > maxhit_per_tuple) {
      qmaskp[tuplectr] = HITQUAL_MULTIHIT;
      continue;
    }
    hip->nhitqual_sortkeyp[seedctr] = nhits;
#endif 

    qmaskp[tuplectr] = HITQUAL_NORMHIT;
    sp = hip->seedp + seedctr;
    sp->posidx = posidx;
    sp->nhits = nhits;
#ifdef hashhit_basqual
    sp->qclass = qclass;
#endif
#ifdef hashhit_debug
    sp->kmer_word = ktup_word;
#endif
    sp->cix = 0;
    sp->qoffs = tuplectr;
    hip->sidxp[seedctr] = seedctr;
    seedctr++;
  }

  for (;tuplectr<hip->qlen;tuplectr++)
    qmaskp[tuplectr] = HITQUAL_TERM;
  hip->n_seeds = seedctr;

  return ERRCODE_SUCCESS;
}


#if (defined hashhit_debug) && (!defined hashhit_basqual)
/** The following routine works only if hhip->nhitqual_sortkeyp does
 * not encode a base quality class in the upper NBIT_KMER_NHIT bits
 * (hashhit_basqual undefined) */

static int getHitInfoStats(const HashHitInfo *hhip, uint32_t *tot, uint32_t *max, uint32_t *min, 
			   uint32_t *median, uint32_t *quart_lo, uint32_t *quart_hi)
 
{
  const uint32_t ns = hhip->n_seeds;
  uint32_t * const nhp = hhip->nhitqual_sortkeyp;

  if (ns < 1)
    return ERRCODE_ASSERT;

  if (tot) {
    uint32_t i, t;
    for (t=0, i=0; i<ns; i++)
      t += nhp[i];
    *tot = t;
  }

  if (((max) || (min) || (median) || (quart_lo) || (quart_hi))) {
    if (!(hhip->status & HITINFO_SORTED))
	return ERRCODE_HITINFO;

    if (max)
      *max = nhp[ns-1];
    if (min)
      *min = nhp[0];
    if (median)
      *median = nhp[ns>>1];
    if (quart_lo)
      *quart_lo = nhp[ns>>2];
    if (quart_hi)
      *quart_hi = nhp[ns - (ns>>2) - 1];
  }

  return ERRCODE_SUCCESS;
}
#endif

#ifdef hashhit_minimise_coverdeficit
static int setHitInfoFrameRanks(COVERAGE mincover, 
				HASHNUM_t maxnhitsum,
				HASHNUM_t maxnhit_perTuple,
				const HashHitInfo *hip)
/**< Fill SEED indices pointing to hip->seedp in set hip->frcp[f].sxp for 
 * each frame 0 <= f < nskip. Set the number of seeds in hip->frcp[f].ns.
 * Then set hip->frcp[f].ns_used
 */
{
  UCHAR f, nFrameUnused, nFrameWithMinCover;
  const UCHAR nskip = hip->nskip;
  uint64_t nhitsum = 0;

  if (hip->n_seeds > 1 &&
      !(hip->status & HITINFO_SORTED))
    return ERRCODE_HITINFO;

  initHitInfoFrames(hip);
  
  /* make sure,at least mincover bases are covered by k-mer words in each frame */
  for (nFrameUnused=0, nFrameWithMinCover =0;
       (nFrameUnused < nskip) && (nFrameWithMinCover < nskip || nhitsum <= maxnhitsum);
       ) {
    for(f=0; f<nskip; f++) {
      FRAMECOV *frp = hip->frcp + f;
      if (frp->status & FRAMECOV_INACTIVE)
	continue;
      if (frp->ns_used < frp->ns) {      
	const SEED *sp = hip->seedp + frp->sxp[frp->ns_used];
	const uint64_t hsum = nhitsum + sp->nhits;
	if (hsum <= maxnhitsum || 
	    (maxnhit_perTuple = 0 || sp->nhits < maxnhit_perTuple)) {
	  SEQLEN_t q =  sp->qoffs/nskip;
	  SEQLEN_t ru = q/sizeof(QMASK);
	  int rb = q - ru*sizeof(QMASK);
	  QMASK qmask = ((QMASK) 1)<<rb;
	  UCHAR hasMinCover = frp->cover >= mincover;
	  for(; q < sp->qoffs+hip->ktup; q+=nskip) {
	    if (!(frp->bufp[ru] & qmask)) {
	      frp->cover += nskip;
	      frp->bufp[ru] |= qmask;
	    }
	    qmask <<= 1;
	    if (!(qmask)) {
	      qmask = 1;
	      ru++;
	    }
	  }
	  nhitsum = hsum;
	  if (!hasMinCover && frp->cover >= mincover)
	    nFrameWithMinCover++;
	  frp->ns_used++;
	} else {
	  frp->status |= FRAMECOV_INACTIVE;
	  nFrameUnused++;
	}
      } else {
	frp->status |= FRAMECOV_INACTIVE;
	nFrameUnused++;
      }
    }
  }
  return ERRCODE_SUCCESS;
}
#endif

static int getHitInfoMaxRank(uint32_t *n, 
			     COVERAGE mincover, 
			     COVERAGE maxcover,
			     uint32_t maxhit, 
			     const HashHitInfo *hip)
     /**< Return n so that the words {0, 1, ..., n-1} with the lowest
      * hit frequencies cover at least mincover bases of the read in
      * each of the nskip frames. If the total number of hits of these
      * words is greater than maxhit, report n such that the total
      * number of hits is less than or equal to maxhit.
      * n is returned such that n >= min(HITINFO_MINSEEDNUM, hip->n_seeds).
      * K-mer words must have been sorted by increasing hit frequency.
      */
{
  short f;
  const UCHAR nskip = hip->nskip;
  const UCHAR ktup = hip->ktup;
#ifndef hashhit_minimise_coverdeficit
  const size_t qmem = hip->qlen*sizeof(UCHAR);
  UCHAR *qbufp = hip->qbufp;
  SEQLEN_t *ixp;
#endif
  SEQLEN_t i, imax, nmax, ntot;
  COVERAGE cover;

  if (hip->n_seeds < 1 || maxcover < mincover)
    return ERRCODE_ASSERT;

  if (hip->n_seeds > 1 &&
      !(hip->status & HITINFO_SORTED))
    return ERRCODE_HITINFO;
#ifdef hashhit_minimise_coverdeficit
  initHitInfoFrames(hip);
#else
  for (f=0; f<nskip; f++)
    hip->countp[f] = 0;

  for (i=0; i<hip->n_seeds; i++) {
    SEQLEN_t ix = hip->sidxp[i];
    f = hip->seedp[ix].qoffs%nskip; /* frame index f from {0, ..., nksip-1} */
    hip->framep[f][hip->countp[f]++] = i; /* save the rank */
  } 
#endif
  /*
   * hip->countp[f] is the number of seeds in each of the nskip frames
   * hip->framep[f][i] is the (sorted) seed index of the ith seed in frame f
   */
  nmax = 0;            /* maximum rank */
#ifdef hashhit_basqual
  ntot = hip->nhitqual_sortkeyp[0]&KMER_NHIT_MASK;
  for (i=1; i<=hip->n_seeds && ntot <= maxhit; i++)
    ntot += hip->nhitqual_sortkeyp[i]&KMER_NHIT_MASK; 
#else
  ntot = hip->nhitqual_sortkeyp[0];
  for (i=1; i<=hip->n_seeds && ntot <= maxhit; i++)
    ntot += hip->nhitqual_sortkeyp[i]; 
#endif
 
  *n = nmax = i-1;

#ifdef hashhit_minimise_coverdeficit
  for(f=0; f<nskip; f++) {
    FRAMECOV *frp = hip->frcp + f;
    cover = 0;
    imax = hip->frcp[f].ns;  /* number of seeds in frame f */
    if (!imax) 
      continue;
    for (i=0; i < frp->ns && cover <= maxcover && (cover < mincover || frp->sxp[i] <= *n); i++) {
      const SEED *sp = hip->seedp + frp->sxp[i];
      SEQLEN_t q =  sp->qoffs/nskip;
      SEQLEN_t ru = q/sizeof(QMASK);
      int rb = q - ru*sizeof(QMASK);
      QMASK qmask = ((QMASK) 1)<<rb;
      for(; q < sp->qoffs+ktup; q+=nskip) {
	if (frp->bufp[ru] & qmask) {
	  cover++;
	} else {
	  frp->bufp[ru] |= qmask;
	}
	qmask <<= 1;
	if (!(qmask)) {
	  qmask = 1;
	  ru++;
	}
      }
    }
    frp->ns_used = i;
    if (i > 0 && frp->sxp[i-1] > nmax)
      nmax = frp->sxp[i-1];
  }
#else
  for(f=0; f<nskip; f++) {
    imax =  hip->countp[f]; /* number of seeds in frame f */
    if (!imax) 
      continue;

    memset(qbufp, 0, qmem);
    cover = 0;

    ixp = hip->framep[f];   /* array of size imax containing seed indices */
    for (i=0; i < imax && cover <= maxcover && (cover < mincover || ixp[i] <= *n); i++) {
      SEQLEN_t ix = hip->sidxp[ixp[i]];
      const SEED *sp = hip->seedp + ix;
      SEQLEN_t q = sp->qoffs;
      for(; q < sp->qoffs+ktup-1; q++) {
	if (!qbufp[q]) {
	  qbufp[q] = 1;
	  cover++;
	}
      }
    }
    if (i > 0 && ixp[i-1] > nmax)
      nmax = ixp[i-1];
  }
#endif

  if (nmax < HITINFO_MINSEEDNUM) {
    *n = (HITINFO_MINSEEDNUM < hip->n_seeds)? HITINFO_MINSEEDNUM: hip->n_seeds;
  } else {
    *n = nmax;
  }
  return ERRCODE_SUCCESS;
}

/******************************************************************************
 ******************** Public Methods of Type HashHitInfo **********************
 ******************************************************************************/

HashHitInfo *hashCreateHitInfo (int blksz, const HashTable *htp)
{
  UCHAR o;
  HashHitInfo *p;
#ifndef hashhit_minimise_coverdeficit
  uint32_t *sp;
  int dblk;
#endif

  EMALLOCP0(p);
  if (!p) return 0;

  if (blksz < 1) blksz = KTUPLIST_DEFBLKSZ;
  ECALLOCP(blksz, p->qmaskp);
  ECALLOCP(blksz, p->qbufp);
  ECALLOCP(blksz, p->sidxp);
  ECALLOCP(blksz, p->nhitqual_sortkeyp);
  ECALLOCP(blksz, p->seedp);
  if (!((p->qmaskp) && (p->qbufp) && (p->sidxp) && 
	(p->nhitqual_sortkeyp) && (p->seedp))) {
    hashDeleteHitInfo(p);
    return 0;
  }
  p->blksz = blksz;
  p->n_alloc = blksz;
  p->ktup = hashTableGetKtupLen(htp, &p->nskip);

#ifdef hashhit_minimise_coverdeficit  
  ECALLOCP(p->nskip, p->frcp);
  if ((p->frcp)) {
    size_t sz = CALC_FRAM_BLKSZ(blksz, p->nskip);
    size_t sq = CALC_FRAM_BUFSZ(blksz, p->nskip);
 
    ECALLOCP(sz * p->nskip, p->frcp->sxp);
    ECALLOCP(sq * p->nskip, p->frcp->bufp);
    if ((p->frcp->sxp) && (p->frcp->bufp)) {
      for (o=1; o<p->nskip;o++) {
	p->frcp[o].sxp = p->frcp[0].sxp + o*sz;
	p->frcp[o].bufp = p->frcp[0].bufp + o*sq;
      }
    }
    p->frcbufsz = sq;
  }
  if (!(p->frcp) || !(p->frcp->sxp) || !(p->frcp->bufp)) {
    hashDeleteHitInfo(p);
    return 0;
  }
#else
  ECALLOCP(p->nskip, p->coverp);
  ECALLOCP(p->nskip, p->framep);
  ECALLOCP(p->nskip, p->countp);
  ECALLOCP(blksz + p->nskip, sp);
  if (!((p->coverp) && (p->framep) && (p->countp) && (sp))) {
      hashDeleteHitInfo(p);
      return 0;
  }
  p->framep[0] = sp;
  dblk = blksz/p->nskip + 1;
  
  for (o=1; o<p->nskip; o++) 
    p->framep[o] = sp + o*dblk;
#endif
  return p;
}

void hashDeleteHitInfo(HashHitInfo *p)
{
  if (p) {
    free(p->qmaskp);
    free(p->qbufp);
    free(p->sidxp);
    free(p->nhitqual_sortkeyp);
    free(p->seedp);
#ifdef hashhit_minimise_coverdeficit
    if (p->frcp) {
      free(p->frcp->sxp);
      free(p->frcp->bufp);
    }
    free(p->frcp);
#else
    free(p->coverp);
    free(p->countp);
    if (p->framep) 
      free(p->framep[0]);
    free(p->framep);
#endif
  }
  free(p);
}

int hashCollectHitInfo(HashHitInfo *hhip, BOOL is_reverse, UCHAR basq_thresh,
		       SEQLEN_t seq_start, SEQLEN_t seq_end,
		       const SeqFastq *seqp, const HashTable *htp)
{
  int errcode;
#ifdef hashhit_debug
  UCHAR ktup = hashTableGetKtupLen(htp, NULL);
#endif

  errcode = collectHitInfo(hhip, is_reverse, 0, basq_thresh, seq_start, seq_end, seqp, htp);

#ifdef hashhit_debug
  fprintHitInfo(stdout, ktup, hhip); 
#endif

  hhip->seed_rank = 0;

  return errcode;
}

int hashCollectHitInfoShort(HashHitInfo *hhip, BOOL is_reverse,
			    HASHNUM_t maxhit_per_tuple,
			    HASHNUM_t maxhit_total,
			    UCHAR basq_thresh,
			    const SeqFastq *seqp, const HashTable *htp)
{
  int errcode;
  SEQLEN_t slen;
  COVERAGE mincover;
#ifndef hashhit_minimise_coverdeficit
  COVERAGE maxcover;
#endif

  if ((errcode = collectHitInfo(hhip, is_reverse, maxhit_per_tuple, 
				basq_thresh, 0, 0, seqp, htp)))
    return errcode;
#ifdef hashhit_debug
  UCHAR ktup = hashTableGetKtupLen(htp, NULL);
  fprintHitInfo(stdout, ktup, hhip); 
#endif

  if (hhip->n_seeds <= 1) {
    hhip->status |= HITINFO_SORTED;
    hhip->seed_rank = hhip->n_seeds;
    return ERRCODE_SUCCESS;
  }
  

  if ((errcode = sort2UINTarraysByQuickSort(hhip->n_seeds, 
					    hhip->nhitqual_sortkeyp, 
					    hhip->sidxp)))
    return errcode;
  
  hhip->status |= HITINFO_SORTED;
  hhip->seed_rank = 0;
  seqFastqGetConstSequence(seqp, &slen, NULL);

#ifdef hashhit_minimise_coverdeficit
  mincover = HITINFO_MINCOVER_NKTUP*hhip->ktup;
  if (mincover + hhip->ktup > slen) {
    /* this may happen for very short reads
     * => relax the constraints */
    mincover = 0;
  }
  if ((errcode = setHitInfoFrameRanks(mincover, 
				      maxhit_total, maxhit_per_tuple,
				      hhip)))
    return errcode;

  hhip->status |= HITINFO_RANK | HITINFO_FRAME;

#else
  mincover = HITINFO_MINCOVER_KMER*hhip->ktup + hhip->nskip;
  maxcover = slen*HITINFO_MAXCOVER_PERCENT/100;
  if (maxcover < hhip->ktup + hhip->nskip)
    maxcover = hhip->ktup + hhip->nskip;
  else if (maxcover > slen - hhip->nskip)
    maxcover = slen - hhip->nskip;

  if (mincover > maxcover) {
    /* this may happen for very short reads
     * => relax the constraints */
    mincover = 0;
    maxcover = slen;
  }
  if ((errcode = getHitInfoMaxRank(&hhip->seed_rank, 
				   mincover, maxcover, 
				   maxhit_total, hhip)))
    return errcode;
  hhip->status |= HITINFO_RANK;
#endif
  
  return ERRCODE_SUCCESS;
}

int hashSortHitInfo(HashHitInfo *hhip)
{
  int errcode = ERRCODE_SUCCESS;

  if (hhip->n_seeds>1 && !(hhip->status & HITINFO_SORTED)) {
    errcode = sort2UINTarraysByQuickSort(hhip->n_seeds, 
					 hhip->nhitqual_sortkeyp, 
					 hhip->sidxp);
    if (!errcode)
      hhip->status |= HITINFO_SORTED;
  }
  return errcode;
}
  
uint32_t hashCalcHitInfoCoverDeficit(const HashHitInfo *hip)
{
  UCHAR s;
  const UCHAR nskip = hip->nskip;
  uint32_t i, deficit;
  COVERAGE d;
  const UCHAR * const qmaskp = hip->qmaskp;
#ifndef hashhit_minimise_coverdeficit
  const UCHAR ktup = hip->ktup;
#endif

  if (hip->status & HITINFO_RANK) {
    COVERAGE cover, maxcover=0;  
#ifndef hashhit_minimise_coverdeficit
    SEQLEN_t *ixp, ix, imax, q;
    UCHAR *qbufp = hip->qbufp;
    const size_t qmem = hip->qlen*sizeof(UCHAR);
    const SEED *sp;
#endif

    d = hip->qlen;
    for (s=0;s<nskip;s++) {
#ifdef hashhit_minimise_coverdeficit
      cover = hip->frcp[s].cover;
#else
      imax =  hip->countp[s]; /* number of seeds in frame s */
      if (!imax) 
	continue;

      memset(qbufp, 0, qmem);
      cover = 0;
      ixp = hip->framep[s];   /* array of size imax containing seed indices */

      for (i=0; i < imax && ixp[i] < hip->seed_rank; i++) {
	ix = hip->sidxp[ixp[i]];
	sp = hip->seedp + ix;
	for(q=sp->qoffs; q<sp->qoffs + ktup; q++) {
	  if (!qbufp[q]) {
	    qbufp[q] = 1;
	    cover++;
	  }
	}
      }
#endif
      if (cover < d)
	d = cover;
      if (cover > maxcover)
	maxcover = cover;
    }
#ifdef hashhit_mscor_calib
    printf("CALDAT_HASHIT %i\n", 100*maxcover/hip->qlen);
#endif
    deficit = maxcover - d + 1;
  } else {
    UCHAR ctr;
    UCHAR k = hip->ktup/nskip;
    if (k>0) k--;
    deficit = 0;
    for (s=0;s<nskip;s++) {
      d = 0;
      for (ctr=0, i=s; i<hip->qlen; i+=nskip) {
	if (qmaskp[i] == HITQUAL_NORMHIT) 
	  ctr = k;
	else if (ctr)
	  ctr--;
	else
	  d += nskip;
      }
      if (d > deficit)
	deficit = d;
    }
  }
  return deficit;
}

uint32_t hashCalcHitInfoNumberOfHits(const HashHitInfo *hhip, HASHNUM_t maxhit_per_tuple)
{
  uint32_t i, nhit, hnum = 0;
  HASHNUM_t hcutoff;

  hcutoff = (maxhit_per_tuple < 1)? 0: maxhit_per_tuple;
  if ((hcutoff)) {
    for (i=0; i<hhip->n_seeds; i++) {
#ifdef hashhit_basqual
      nhit = hhip->nhitqual_sortkeyp[i]&KMER_NHIT_MASK;
#else
      nhit = hhip->nhitqual_sortkeyp[i];
#endif
      if (nhit <= hcutoff)
	hnum += nhit;
    }
  } else {
    for (i=0; i<hhip->n_seeds; i++) {
#ifdef hashhit_basqual
      hnum += hhip->nhitqual_sortkeyp[i]&KMER_NHIT_MASK;
#else
      hnum += hhip->nhitqual_sortkeyp[i];
#endif
    }
  }

  return hnum;
}

uint32_t hashHitInfoCalcHitNumbers(const HashHitInfo *hhip, uint32_t *nhit_rank)
{
  uint32_t i, ns, nr = 0;
  
  ns = (hhip->seed_rank > 0)? hhip->seed_rank: hhip->n_seeds;
  for (i=0; i<ns; i++)
#ifdef hashhit_basqual
    nr += hhip->nhitqual_sortkeyp[i]&KMER_NHIT_MASK;
#else
  nr += hhip->nhitqual_sortkeyp[i];
#endif
  *nhit_rank = nr;
  for (; i<hhip->n_seeds; i++)
#ifdef hashhit_basqual
    nr += hhip->nhitqual_sortkeyp[i]&KMER_NHIT_MASK;
#else
  nr += hhip->nhitqual_sortkeyp[i];
#endif
  return nr;
}

/******************************************************************************
 ******************** Private Methods of Type HashHitList *********************
 ******************************************************************************/
void blankHitList(HashHitList *hlp)
{
  if (hlp) {
    hlp->nhits=0;
    hlp->status=0;
    memset(hlp->qmask, HITQUAL_NOHIT, hlp->qlen);
  }
}

static int reallocHitList(HashHitList *hlp, size_t newsiz)
{
  size_t nsiz = (newsiz + hlp->nhits_blksz - 1)/hlp->nhits_blksz;
  void *hp;

  nsiz *= hlp->nhits_blksz;
  if (nsiz > INT_MAX)
    return ERRCODE_OVERFLOW;

  hp = EREALLOCP(hlp->sqdat, nsiz);
  if (!hp) return ERRCODE_NOMEM;
  hlp->sqdat = hp;
  hlp->nhits_alloc = (int) nsiz; 

  return ERRCODE_SUCCESS;
}

static int reallocQmask(HashHitList *hlp, uint32_t seqlen)
{
  size_t nsiz = ((uint32_t) (seqlen/HASH_QMASKBLOCK)+1)*HASH_QMASKBLOCK;
  char *qmaskp = EREALLOC(hlp->qmask, nsiz);
  if (!qmaskp) return ERRCODE_NOMEM;
  hlp->qmask = qmaskp;
  hlp->qmask_alloc = nsiz;
  hlp->qlen = seqlen;

  return ERRCODE_SUCCESS;
}

static int initHitList(HashHitList *hlp, const HashHitInfo *hip)
{
  int errcode;
  //size_t target_size = qlen*HITLST_QLENSIZ_FACT;
  size_t target_size = hip->qlen*log((double) hip->qlen)*HITLST_LOGQLENSIZ_FACT;

  hlp->qlen = hip->qlen;
  hlp->ktup = hip->ktup;
  hlp->nskip = hip->nskip;

  if (target_size > HITLST_MAXSIZ)
    target_size = HITLST_MAXSIZ;
  else if (target_size < HITLST_MINSIZ)
    target_size = HITLST_MINSIZ;
  
  if (target_size > INT_MAX)
    target_size = INT_MAX;

  if (((int) target_size) > hlp->nhits_alloc &&
      ((errcode = reallocHitList(hlp, target_size))))
    return errcode;
  
  if (hip->qlen >= hlp->qmask_alloc &&
      ((errcode = reallocQmask(hlp, hip->qlen))))
    return errcode;

  hlp->nhits_max = (int) target_size;

  blankHitList(hlp);
  
  if (hip->status &  HITINFO_REVERSE)
    hlp->status |= HASHHITSTAT_REVERSE;

  return ERRCODE_SUCCESS;
}

#ifdef hashhit_minimise_coverdeficit
static int 
fillHitListFromSegment
(
 HashHitList *hlp,
#ifdef RESULTS_TRACKER
 Track *trkp,
#endif
 HASHPOS_t segpos_lo,
 HASHPOS_t segpos_hi,
 const HashHitInfo *hip,
 const HashTable *htp,
 const HashHitFilter *hhfp
 )
/**< Fill the hit list with seeds in the segment [pos_lo, pos_hi]
 * of the concatenated reference sequences (e.g. a chromosome).
 * This function must be called successively
 * for non-overlapping segments in the order of their position along the
 * concatenated set of reference sequences (Seeds are added starting from where
 * the previous call stopped).
 */
{
  int errcode;
  const UCHAR is_reverse = hip->status & HITINFO_REVERSE;
  UCHAR f;
  const UCHAR nskip = hip->nskip;
  const HASHPACK_t offbit = ((HASHPACK_t) 1)<<(HASHHIT_HALFBIT+1);
#ifdef RESULTS_TRACKER
  HASHPOS_t trkpos_lo, trkpos_hi;
  BOOL is_trk_overlap;
  TRACKFLG_t trk_flg;
  if (trkp != NULL) {
    trackGetKtupData(trkp, &trkpos_lo, &trkpos_hi, &trk_flg);
    is_trk_overlap = trk_flg & TRACKFLG_KMERHITS_OVERLAP;
  }
#endif
 
  if ((errcode = initHitList(hlp, hip)))
    return errcode;

  /* reset the number of already used k-mers */
  for (f=0; f<nskip; f++) {
    FRAMECOV *frcp = hip->frcp + f;
    SEQLEN_t s;
    for (s=0; s<frcp->ns_used; s++) {
      SEED * const seedp = hip->seedp + s;
      HASHPOS_t *posp;
      HASHPACK_t *sqdatp;
      HASHNUM_t i, n, nhits_remaining;
      const HASHNUM_t nhits = hashTableFetchHitPositions(&posp, htp, seedp->posidx);
      
      if (nhits != seedp->nhits)
	return ERRCODE_ASSERT;

      if (seedp->cix >= nhits) {
	if (posp[nhits-1] < segpos_lo)
	  continue;
	seedp->cix = 0;
      }
      
      if (seedp->cix > 0 && posp[seedp->cix-1] > segpos_lo)
	seedp->cix = 0; /* overlapping with last call -> reset helper */
      for (i = seedp->cix; i<nhits && posp[i] < segpos_lo; i++);
      
      posp += i;
      nhits_remaining = nhits - i;
      seedp->cix += i;

      if (hlp->nhits + nhits_remaining > INT_MAX)
	return ERRCODE_OVERFLOW;

      if (((int) (hlp->nhits + nhits_remaining)) > hlp->nhits_alloc &&
	  (errcode = reallocHitList(hlp, hlp->nhits + nhits_remaining)))
	return errcode;
   
      sqdatp = hlp->sqdat + hlp->nhits;
      n = 0;

      if (hhfp == NULL) {
	for (i=0; i<nhits_remaining && posp[i] < segpos_hi; i++) {
#ifdef RESULTS_TRACKER
	  if (trkp != NULL && !is_trk_overlap && posp[i] >= trkpos_lo && posp[i] <= trkpos_hi)
	    is_trk_overlap = 1;
#endif
	  SET_NEXT_SHIFT(is_reverse, sqdatp[i], posp[i], seedp->qoffs, nskip, offbit);
	}
	n = i;
      } else {
	short j;
	for (i=0, j=0; i<nhits_remaining && posp[i] < segpos_hi; i++) {
	  const FILTERIVAL * const ivalp = hhfp->ivp + j;
	  if (posp[i] < ivalp->lower) continue;
	  for(; posp[i] > ivalp->upper && j < hhfp->num; j++);
	  if (j>=hhfp->num) break;
	  if (posp[i] >= ivalp->lower) {
#ifdef RESULTS_TRACKER
	    if (trkp != NULL && !is_trk_overlap && posp[i] >= trkpos_lo && posp[i] <= trkpos_hi)
	      is_trk_overlap = 1;
#endif
	    SET_NEXT_SHIFT(is_reverse, sqdatp[n], posp[i], seedp->qoffs, nskip, offbit);
	    n++;
	  }
	}
      }
      seedp->cix += i;
      hlp->nhits += n;
    }
  }
#ifdef RESULTS_TRACKER
  if (trkp != NULL && (is_trk_overlap))
    trackSetFlag(trkp, TRACKFLG_KMERHITS_OVERLAP);
#endif

  return ERRCODE_SUCCESS;
}

#else //ifdef hashhit_minimise_coverdeficit

static int 
fillHitListFromHitInfoSegment
(
 HashHitList *hlp,
#ifdef RESULTS_TRACKER
 Track *trkp,
#endif
 HASHPOS_t segpos_lo,
 HASHPOS_t segpos_hi,
 HASHNUM_t  maxhit_per_tuple,
 BOOL use_short_hitinfo,
 const HashHitInfo *hip,
 const HashTable *htp,
 const HashHitFilter *hhfp
 )
/**< Fill the hit list with seeds in the segment [pos_lo, pos_hi] of
 * the concatenated reference sequences. This could be e.g. a
 * chromosome. Start with hip->sidxp[0] and leave out ktuples with
 * hits above maxhit_per_tuple. This function must be called
 * successively for non-overlapping semgments in the order of their
 * position along the concatenated set of reference sequences (Seeds
 * are added starting from where the previous call stopped).
 *
 * segsetoffs_lo and segsetoffs_hi are k-tuple numbers.
 */
{
  int errcode;
  const UCHAR is_reverse = hip->status & HITINFO_REVERSE;
  const UCHAR nskip = hip->nskip;
  HASHNUM_t nhits, nh, i, k;
  short j;
  SEQLEN_t n, tuplectr;
  const SEQLEN_t n_seeds = 
    (use_short_hitinfo) && hip->seed_rank > 0 ? 
    hip->seed_rank: hip->n_seeds;
  SEED *seedp;
  HASHPOS_t *posp;
  HASHPACK_t *sqdatp;
  const HASHPACK_t offbit = ((uint64_t) 1)<<(HASHHIT_HALFBIT+1);
  UCHAR * const qmaskp = hip->qmaskp;
  const FILTERIVAL *ivalp = 0;
#ifdef RESULTS_TRACKER
  HASHPOS_t trkpos_lo, trkpos_hi;
  BOOL is_trk_overlap = 0;
  TRACKFLG_t trk_flg;
  if (trkp != NULL) {
    trackGetKtupData(trkp, &trkpos_lo, &trkpos_hi, &trk_flg);
    is_trk_overlap = trk_flg & TRACKFLG_KMERHITS_OVERLAP;
  }
#endif

  if ((errcode = initHitList(hlp, hip)))
    return errcode;
 
  for (n=0; n<n_seeds; n++) {
    seedp = hip->seedp + ((use_short_hitinfo)? hip->sidxp[n]: n);
    if (maxhit_per_tuple > 0 && 
#ifdef hashhit_basqual
	(hip->nhitqual_sortkeyp[n]&KMER_NHIT_MASK) > maxhit_per_tuple
#else
	hip->nhitqual_sortkeyp[n] > maxhit_per_tuple
#endif
	) {
      qmaskp[seedp->qoffs] = HITQUAL_MULTIHIT;
      continue;
    }
    nhits = hashTableFetchHitPositions(&posp, htp, seedp->posidx);
    if (seedp->cix >= nhits) {
      if (posp[nhits-1] < segpos_lo)
	continue;
      seedp->cix = 0;
    }
    if (posp[seedp->cix] > segpos_lo)
      seedp->cix = 0; /* reset helper */
    posp += seedp->cix;
    nh = nhits - seedp->cix;
    for (i=0; i<nh && posp[i] < segpos_lo; i++);
    nh -= i;
    seedp->cix += i;
    posp += i;

    if (hlp->nhits + nh > (HASHNUM_t) hlp->nhits_alloc) {
      if (maxhit_per_tuple > 0)
	return ERRCODE_ALLOCBOUNDARY;
      qmaskp[seedp->qoffs] = HITQUAL_MULTIHIT;
      continue;
    }
     
    tuplectr = seedp->qoffs;
    sqdatp = hlp->sqdat + hlp->nhits;
    k = 0;

    if (hhfp) {
      /* filtered */
      for (i=0, j=0; i<nh && posp[i] < segpos_hi; i++) {
	ivalp = hhfp->ivp+j;
	if (posp[i] < ivalp->lower) continue;
	for(; posp[i] > ivalp->upper && j < hhfp->num; j++)
	  ivalp = hhfp->ivp+j;
	if (j>=hhfp->num) break;
	if (posp[i] >= ivalp->lower) {
#ifdef RESULTS_TRACKER
	  if (trkp != NULL && !is_trk_overlap && posp[i] >= trkpos_lo && posp[i] <= trkpos_hi)
	    is_trk_overlap = 1;
#endif
	  SET_NEXT_SHIFT(is_reverse, sqdatp[k], posp[i], tuplectr, nskip, offbit);
	  k++;
	}
      }
    } else {
      /* unfiltered */
      for (i=0; i<nh && posp[i] < segpos_hi; i++) {
#ifdef RESULTS_TRACKER
	if (trkp != NULL && !is_trk_overlap && posp[i] >= trkpos_lo && posp[i] <= trkpos_hi)
	  is_trk_overlap = 1;
#endif
	SET_NEXT_SHIFT(is_reverse, sqdatp[i], posp[i], tuplectr, nskip, offbit);
      }
      k = i;
    }
    seedp->cix += i;
    hlp->nhits += k;
  }

#ifdef RESULTS_TRACKER
  if (trkp != NULL && (is_trk_overlap))
    trackSetFlag(trkp, TRACKFLG_KMERHITS_OVERLAP);
#endif

  return ERRCODE_SUCCESS;
}
#endif //ifdef hashhit_minimise_coverdeficit

/******************************************************************************
 ******************** Public Methods of Type HashHitList **********************
 ******************************************************************************/

HashHitList *hashCreateHitList(int maxnhits)
{
  HashHitList *hlp = 0;

  EMALLOCP0(hlp);
  if (hlp) {
    if (maxnhits < HITLST_MINSIZ) maxnhits = HITLST_MINSIZ;
    ECALLOCP(maxnhits, hlp->sqdat);
    ECALLOCP(HASH_QMASKBLOCK, hlp->qmask);
    if ((hlp->sqdat) && (hlp->qmask)) {
      hlp->nhits = 0;
      hlp->nhits_max = maxnhits;
      hlp->nhits_alloc = maxnhits;
      hlp->nhits_blksz = HITLST_DEFBLKSZ;
      hlp->qmask_alloc = HASH_QMASKBLOCK;
      hlp->qlen = 0;
      hlp->status = 0;
    } else {
      hashDeleteHitList(hlp);
      hlp = NULL;
    }
  }

  return hlp;
}

void hashDeleteHitList(HashHitList *hlp)
{
  if (hlp) {
    free(hlp->sqdat);
    free(hlp->qmask);
  }
  free(hlp);
}

void hashBlankHitList(HashHitList *hlp)
{
  blankHitList(hlp);
}

int hashCollectHitsUsingCutoff(HashHitList *hlp,
#ifdef hashhit_dump_sortarray
			       FILE *dumpfp,
			       int *dumpctr,
#endif
			       HASHNUM_t max_nhit_per_tup, 
			       const HashTable *htp, const HashHitInfo *hip)
{
  int errcode;
  HASHNUM_t j, nhits;
  BOOL reached_ceiling = 0;
  HASHPOS_t *posp;
  SORTKEY_t nh;
  SEQLEN_t i, q, qo;
  const SEQLEN_t n_seeds = ((hip->seed_rank))? hip->seed_rank: hip->n_seeds;
  const HASHPACK_t offbit = ((uint64_t) 1)<<(HASHHIT_HALFBIT+1);
  HASHPACK_t *dp;
  const SEED *sp;

  if ((errcode = initHitList(hlp, hip)))
    return errcode;

  do { 
    reached_ceiling = 0;
    blankHitList(hlp);
    if (hip->status &  HITINFO_REVERSE)
      hlp->status |= HASHHITSTAT_REVERSE;

    for (i=0 ; i<n_seeds; i++) {
#ifdef hashhit_basqual
      nh = hip->nhitqual_sortkeyp[i]&KMER_NHIT_MASK;
#else
      nh = hip->nhitqual_sortkeyp[i];
#endif
      if (nh < 1) 
	continue;
      sp = hip->seedp + hip->sidxp[i];
      q = sp->qoffs;	

      if (max_nhit_per_tup > 0 && nh > max_nhit_per_tup) {
	hlp->qmask[q] = HITQUAL_MULTIHIT;
	continue;
      }
      if (hlp->nhits + nh > INT_MAX)
	return ERRCODE_OVERFLOW;

      if (((int) (hlp->nhits + nh)) > hlp->nhits_max) {
	reached_ceiling = 1;
	break;
      }

      nhits = hashTableFetchHitPositions(&posp, htp, sp->posidx);
     
      if (nh != nhits)
	return ERRCODE_ASSERT;

      hlp->qmask[q] = HITQUAL_NORMHIT;
      qo = q/hlp->nskip;
      dp = hlp->sqdat + hlp->nhits;
#ifdef hashhit_debug
      printf("hashhit_debug: ==== hits for qo = %u ====\n", q);
      for (j=0; j<nhits; j++)
	printf("hashhit_debug: so = %u\n", posp[j]);
#endif
      if (hlp->status & HASHHITSTAT_REVERSE) {
	for (j=0; j<nhits; j++) 
	  dp[j] = (((HASHPACK_t) posp[j] + qo)<<HASHHIT_HALFBIT) + q;
      } else {
	for (j=0; j<nhits; j++) 
	  dp[j] = (((((HASHPACK_t) posp[j])|offbit) - qo)<<HASHHIT_HALFBIT) + q;
      }
#ifdef hashhit_debug
      for (j=0; j<nhits; j++) {
	if ((dp[j]&HASHHIT_HALFMASK) > hip->qlen) {
	  printf("hashhit_debug: fatal error in query offset %u\n", posp[j]);
	  return ERRCODE_ASSERT;
	}
      }
#endif
      hlp->nhits += nhits;
    }
    max_nhit_per_tup /= 2;
  } while((reached_ceiling) && max_nhit_per_tup > MINHIT_PER_TUPLE);
#ifdef hashhit_dump_sortarray
  if (hlp->nhits >= MIN_ARRLEN) {
    if (dumpfp != NULL)
      if (hlp->nhits != fwrite(hlp->sqdat, sizeof(uint64_t), hlp->nhits, dumpfp))
	return ERRCODE_FILEIO;
    if (dumpctr != NULL)
      (*dumpctr)++;
  }
#endif
  errcode = sortUINT64arrayByQuickSort(hlp->nhits, hlp->sqdat);
  if (!(errcode))
    hlp->status |= HASHHITSTAT_SHIFTS;
  return errcode;
}

int hashCollectHitsForSegment(HashHitList *hlp,
#ifdef RESULTS_TRACKER
			      Track *trkp,
#endif
#ifdef hashhit_dump_sortarray
			      FILE *dumpfp,
			      int *dumpctr,
#endif
			      SETSIZ_t segmoffs_lo,
			      SETSIZ_t segmoffs_hi,
#ifndef hashhit_minimise_coverdeficit
			      HASHNUM_t nhit_max,
			      BOOL use_short_hitinfo,
#endif
			      const HashHitInfo *hhip,
			      const HashTable *htp,
			      const HashHitFilter *hhfp)
{
  int errcode = ERRCODE_SUCCESS;
  UCHAR nskip;

  hashTableGetKtupLen(htp, &nskip);
  segmoffs_lo /= nskip;
  if (segmoffs_lo > HASHPOS_MAX)
    return ERRCODE_ARGRANGE;
  segmoffs_hi /= nskip;
  if (segmoffs_hi > HASHPOS_MAX)
    segmoffs_hi = HASHPOS_MAX;
  /* if (nhit_max < 1) nhit_max = DEFAULT_MAXHIT_PER_TUPLE; */
#ifdef hashhit_minimise_coverdeficit
  errcode = fillHitListFromSegment(hlp, 
#ifdef RESULTS_TRACKER
				   trkp,
#endif
				   (HASHPOS_t) segmoffs_lo, 
				   (HASHPOS_t) segmoffs_hi,
				   hhip, htp, hhfp);  
  if (errcode) return errcode;
#else
  do {
    errcode = fillHitListFromHitInfoSegment(hlp, 
#ifdef RESULTS_TRACKER
					    trkp,
#endif
					    (HASHPOS_t) segmoffs_lo, (HASHPOS_t) segmoffs_hi,
					    nhit_max,
					    use_short_hitinfo,
					    hhip, htp, hhfp);  
    nhit_max /= 2;
  } while (errcode == ERRCODE_ALLOCBOUNDARY && 
	   nhit_max > MINHIT_PER_TUPLE);
  
  if ((errcode) && errcode != ERRCODE_ALLOCBOUNDARY) 
    return errcode;
#endif

#ifdef hashhit_debug
  printf("\n=-=-= Hash Hit list after collection =-=-=-=\n");
  hashPrintHitList(hlp, stdout);
#endif
#ifdef hashhit_dump_sortarray
  if (hlp->nhits >= MIN_ARRLEN) {
    if (dumpfp != NULL)
      if (hlp->nhits != fwrite(hlp->sqdat, sizeof(uint64_t), hlp->nhits, dumpfp))
	return ERRCODE_FILEIO;
    if (dumpctr != NULL)
      (*dumpctr)++;
  }
#endif
#ifdef hashhit_mergesort
  errcode = sortMergeSort(hlp->sqdat, hlp->nhits);
#else
  errcode = sortUINT64arrayByQuickSort(hlp->nhits, hlp->sqdat);
#endif
  if (!(errcode))
    hlp->status |= HASHHITSTAT_SHIFTS;

  return errcode;
}

void hashPrintHitList(const HashHitList *hlp, FILE *fp)
{
  int i;
  uint32_t qo,so;
  fprintf(fp, "=-= List (%s) of %d hits: =-=\n", 
	  (hlp->status & HASHHITSTAT_REVERSE)? "reverse":"forward", hlp->nhits);
  fprintf(fp, " hit_num | qoffs | ktupno offs | shift | \n");
  for(i=0; i<hlp->nhits; i++) { 
    qo = hlp->sqdat[i]&HASHHIT_HALFMASK;
    if (hlp->status & HASHHITSTAT_REVERSE) 
      so = (hlp->sqdat[i]>>HASHHIT_HALFBIT) - qo/hlp->nskip;
    else
      so = ((hlp->sqdat[i]>>HASHHIT_HALFBIT) + qo/hlp->nskip)&MASK32BIT;
    fprintf(fp, "%d %u %u %llu\n", i, qo, (so&HASHHIT_SOFFSMASK),
	    (long long unsigned int) (hlp->sqdat[i]>>HASHHIT_HALFBIT));
  }
  fprintf(fp, "=-= End of list =-=\n\n");
}

int hashCheckHitList(const HashHitList *hlp, const SeqFastq *seqp, 
		     const HashTable *htp, const SeqSet *ssp,
		     const SeqCodec *codep)
     /* Check whether all entries in the hit list are correct.
      * Whether the hit list is also complete is not checked.
      */
{
  char code;
  UCHAR k, ktup, nskip;
  const char *query_datap, *qdp, *sdp;
  uint32_t querylen, qoffs;
  uint64_t soffs, soffs_end, sqdat;
  int j, errcode;
  SeqFastq *ucp = NULL;

  if (!(hlp->status & HASHHITSTAT_SHIFTS))
    return ERRCODE_HITSTATUS;

  if (!(ucp = seqFastqCreate(0, SEQTYP_FASTA)))
    return ERRCODE_NOMEM;
 
  ktup = hashTableGetKtupLen(htp, &nskip);
  if (ktup >= hlp->qlen)
    return ERRCODE_ASSERT;

  query_datap = seqFastqGetConstSequence(seqp, &querylen, &code);
  if (querylen != hlp->qlen) 
    return ERRCODE_ASSERT;

  if (code != SEQCOD_MANGLED) 
    return ERRCODE_SEQCODE;
  
  for (j=0; j<hlp->nhits; j++) {
    sqdat = hlp->sqdat[j];
    qoffs = sqdat&HASHHIT_HALFMASK;
    if (hlp->status & HASHHITSTAT_REVERSE)
      soffs = (sqdat>>HASHHIT_HALFBIT) - qoffs/nskip;
    else
      soffs = ((sqdat>>HASHHIT_HALFBIT) + qoffs/nskip)&MASK32BIT;

    soffs *= nskip;
    if (qoffs >= querylen) {
      printf("Wrong base offsets. qoffs %u , soffs %llu\n",
	     qoffs, (long long unsigned int) soffs);
      return ERRCODE_FAILURE;
    }
    soffs_end = soffs + ktup - 1;
    if ((errcode = seqSetFetchSegment(ucp, &soffs, &soffs_end, ssp, codep)))
      return errcode;
    sdp = seqFastqGetConstSequence(ucp, NULL, &code);
    if (code == SEQCOD_ASCII) {
      if ((errcode = seqFastqEncode(ucp, codep)))
	return errcode;
      sdp = seqFastqGetConstSequence(ucp, NULL, &code);
    }

    if (code != SEQCOD_MANGLED) 
      return ERRCODE_SEQCODE;

    if (hlp->status & HASHHITSTAT_REVERSE) {
      qdp = query_datap + qoffs + ktup - 1;
      for (k=0; k<ktup; k++)
	    if (((*qdp--) & 
		 SEQCOD_STDNT_MASK) != (((*sdp++)^SEQCOD_STDNT_MASK)&SEQCOD_STDNT_MASK))
	      return ERRCODE_FAILURE;
    } else {
      qdp = query_datap + qoffs;
      for (k=0; k<ktup; k++)
	if (((*qdp++) & SEQCOD_STDNT_MASK) != ((*sdp++) & SEQCOD_STDNT_MASK))
	  return ERRCODE_FAILURE;
    }
  }
  seqFastqDelete(ucp);

  return ERRCODE_SUCCESS;
}

const uint64_t *hashGetHitListData(int *nhits, char *is_reverse,
				   uint32_t *qlen,
				   UCHAR *ktup, UCHAR *nskip,
				   const char **qmask,
				   const HashHitList *hlp)
{
  if (nhits) *nhits = hlp->nhits;
  if (is_reverse) *is_reverse = (hlp->status & HASHHITSTAT_REVERSE) != 0;
  if (qlen) *qlen = hlp->qlen;
  if (ktup) *ktup = hlp->ktup;
  if (qmask) *qmask = hlp->qmask;
  if (nskip) *nskip = hlp->nskip;

  return hlp->sqdat;
}

uint32_t hashCalcHitListCoverDeficit(const HashHitList *hlp)
{
  uint32_t i, d, deficit;
  UCHAR s, ctr, nskip = hlp->nskip;
  UCHAR k = hlp->ktup/nskip;
  const char *qmaskp = hlp->qmask;
  if (k>0) k--;
  deficit = 0;
  for (s=0;s<nskip;s++) {
    d = 0;
    for (ctr=0, i=s; i<hlp->qlen; i+=nskip) {
      if (qmaskp[i] == HITQUAL_NORMHIT) 
	ctr = k;
      else if (ctr)
	ctr--;
      else
	d += nskip;
    }
    if (d > deficit)
      deficit = d;
  }

  return deficit;
}

/******************************************************************************
 ********************* Private Methods of Type FILTERIVAL *********************
 ******************************************************************************/
static int cmpFILTERIVAL(const void *this, const void *other)
     /**< Sorting in ascending order by lower bound */
{
  HASHPOS_t a = ((const FILTERIVAL *) this)->lower;
  HASHPOS_t b = ((const FILTERIVAL *) other)->lower;
  if (a<b) return -1;
  if (a>b) return 1;
  return 0;
}

/******************************************************************************
 ******************* Private Methods of Type HashHitFilter ********************
 ******************************************************************************/

static int reallocHitFilter(HashHitFilter *hhfp, short siz)
{
  size_t nsiz;
  void *hp;
  if (siz < 1) return ERRCODE_ARGRANGE;

  nsiz = ((siz-1)/hhfp->blocksiz + 1)*hhfp->blocksiz;
  if (nsiz > SHRT_MAX) return ERRCODE_OVERFLOW;
  if (!(hp = EREALLOCP(hhfp->ivp, nsiz)))
    return ERRCODE_NOMEM;
  hhfp->ivp = hp;
  hhfp->n_alloc = (short) nsiz;
  return ERRCODE_SUCCESS;
}
/******************************************************************************
 ******************** Public Methods of Type HashHitFilter ********************
 ******************************************************************************/
HashHitFilter *hashCreateHitFilter(short blocksiz)
{
  HashHitFilter *p;
  EMALLOCP0(p);
  if (!p) return 0;

  if (blocksiz<1) blocksiz = DEFAULT_FILTER_BLOCKSIZ;
  ECALLOCP(blocksiz, p->ivp);
  if (!p->ivp) {
    hashDeleteHitFilter(p);
    p = 0;
  } else {
    p->n_alloc = blocksiz;
    p->blocksiz = blocksiz;
  }
  return p;
}

void hashDeleteHitFilter(HashHitFilter *p)
{
  if (p) {
    free(p->ivp);
  } 
  free(p);
}

void hashResetHitFilter(HashHitFilter *p)
{
  p->num = 0;
}


int hashAddToHitFilter(HashHitFilter *hhfp,
		       HASHPOS_t lo, HASHPOS_t hi)
{
  int errcode;
  FILTERIVAL *hp;
  if (hhfp->num >= hhfp->n_alloc &&
      (errcode = reallocHitFilter(hhfp,hhfp->num+1)))
    return errcode;
  hp = hhfp->ivp + hhfp->num++;
  hp->upper = hi;
  hp->lower = lo;
  return ERRCODE_SUCCESS;
}

void hashPruneHitFilter(HashHitFilter *hhfp)
{
  short i,j;
  FILTERIVAL *ivp = hhfp->ivp;

  qsort(ivp, hhfp->num, sizeof(FILTERIVAL), cmpFILTERIVAL);
  for (i=0, j=1; j<hhfp->num; j++) {
    if (ivp[j].lower < ivp[i].upper) {
      /* overlapping intervals */
      if (ivp[j].upper > ivp[i].upper)
	ivp[i].upper = ivp[j].upper;
    } else {
      if (++i < j) 
	ivp[i] = ivp[j];
    }   
  }
  hhfp->num = (i<j)? i+1:j;
  return; 
}
