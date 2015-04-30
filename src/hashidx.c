/** Hash table of words of fixed length (perfect hash function) */

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
#include <string.h>

#include "elib.h"
#include "filio.h"
#include "hashidx.h"

enum {
  FALSE = 0,
  TRUE = 1,            /**< for booleans */
  ARRAY_MARGIN = 2,    /**< margin for memory allocation in hash index */
  MAXTESTOUT = 5,
  MAXTUPLEN = 21,         /**< maximum k-tuple length i.e. sizeof(HASHKEY_t)*4 */
  MAX_NBITS_KEY = 32,     /**< maximum number of bits for key */
  HASHTABFIL_OLDVERSION = 2, /**< Previous version number (reading is backwards compatible). */
  HASHTABFIL_VERSION = 3, /**< Version number of the hash table format */
  HASHTABFIL_HEADSIZ = 8, /**< size of the Hash table file specific header as
			   * the number of 32-bit words */
  HASHSEG_BLKSZ_DEFAULT = 256, /**< Default block size for hash segments */
  BLKSIZ_IDXPOS = 8192,   /**< Block size for allocation of position array HashTable.pos */
#ifdef hashidx_debug
  MEGABYT = 1000*1000,
#endif
};

enum HASHSETUP_STATUS {
  HASHSETUP_EMPTY = 0,    /**< Hash table is emtpy (i.e. has just been created) */
  HASHSETUP_IDXALLOC = 1, /**< Index for keys has been allocated */
  HASHSETUP_COMPLETE = 2, /**< Hash table setup is complete */
};

static const char HASHTABFIL_NAMEXT[] = "smi"; /* file name extension for hash table file */
static const char HASHTABFIL_WRITERRMSG[] = "when writing hash table file";
static const char HASHTABFIL_READERRMSG[] = "when reading hash table file";

typedef uint32_t HASHKEY_t;
typedef uint8_t UCHAR_t;
typedef uint32_t HALFWORD_t;
typedef uint32_t WORDNUM_t;

/** Hash table for a series of sequence chunks.  
 * 
 * The sequence chunks are hashed in k-tuples of wordlen successive
 * bases. The k-tuple start bases are nskip bases apart. That means
 * for starting positions k_tuple[i], k_tuple[i+1] of successive
 * k-tuples nskip = k_tuple[i+1] - k_tuple[i].
 *
 * The hash function generates a perfect hash key by encoding each of
 * the wordlen bases b[0], b[1], ..., b[wordlen-1] of the k-tuple in 2
 * bits. The hash key is returned as an unsigned 32-bit integer which
 * dictates 0 < wordlen <= 16. The sequence of bits in the key is such
 * that the least significant bits 0 and 1 of the key encode base b[wordlen-1] and
 * the most significant bits (wordlen-1)*2 and (wordlen-1)*2+1 encode base b[0].
 *
 * K-tuple positions are stored as k-tuple serial numbers taken along
 * the concatenated sequences of type SeqSet. These serial numbers are
 * stored in the array pos[]. Indices for pos[] can be looked up from
 * the hash key in the array idx[] such that for a given hash_key the
 * serial number of the first k-tuple encoded by the key is returned
 * by pos[idx[hash_key]]. The number of such k-tuples is
 * idx[hash_key+1] - idx[hash_key]. The postition of the last k-tuple
 * is pos[idx[hash_key+1]-1].
 * 
 * With the number of distinct k-tuples (hash keys) wordlen_num =
 * 2^(2*wordlen), and n[ktup_num-1] = idx[ktup_num]-idx[ktup_num-1] the
 * number of occurrences of the k-tuple with key (ktup_num-1), the
 * value of the last entry in idx[] is defined by idx[ktup_num] =
 * idx[ktup_num-1] + n[ktup_num-1].
 *
 * \note SSAHA2 calculates the key from right to left, i.e. bits (0,1) encode
 * base[0].
 * \note maximum allowed sequence length is the maximum number that can be stored
 *       in a signed integer.
 */

struct _HashTable {    /**< Hash index*/
  UCHAR_t typ;           /**< One of HASHIDX_TYPES */
  UCHAR_t status;        /**< One of HASHTAB_STATUS */
  UCHAR_t wordlen;       /**< Number of bases in the hashed words */
  UCHAR_t nskip;         /**< number of bases to skip between hashed k-tuples */
  UCHAR_t nbits_key;
  UCHAR_t nbits_lo;    /**< number of bits for the perfect part of the hash key */
  HASHKEY_t nkeys;     /**< number of keys 2^nbits_key (==4^wordlen if 
			* typ == HASHIDXTYP_PERFECT) */
  HASHWORD_t keymask;   /**< nkeys-1 mask out the relevant bits in lookup */
  HASHKEY_t keymask_lo; /**< for masking out the perfect part of the hash key */
  HASHKEY_t keymask_hi; /**< for masking out the imperfect part of the hash key */
  HASHKEY_t keymod;     /**< modulo for calculating the imperfect part
			 * of the hash key, keymod = 2^(nbits_key -
			 * nbits_perf) */
  HASHWORD_t wordmask;    /**< nkeys-1 mask out the relevant bits in lookup */
  HASHWORD_t wordmask_lo; /**< for masking out the part of the hashed
			   * ktuple that forms the perfect part of the
			   * hash key */
  HASHWORD_t wordmask_hi;  /**< for masking out the part of the hashed
			   * ktuple that forms the imperfect part of the
			   * hash key */
  HASHNUM_t *idx;       /**< array of (nkeys+1) indices of array pos[]. 
			*   The number of hits for a given key is idx[key+1]-idx[key].
			*   There is no reason to choose sizeof(HASHNUM_t) > sizeof(SEQPOS)
			*   However SSAHA2 implements HASHNUM_t as 64-bit signed integer and
			*   SEQPOS as 32-bit signed integer */
  HASHPOS_t *pos;       /**< Positions of k-tuples for each key in the sequence.
			*   pos[idx[key]], pos[idx[key]+1], ..., pos[idx[key+1]-1]
			*   the positions are specified as k-tuple serial numbers.
			*   This allows to cover a factor nskip more bases with
			*   a given number of bits for SEQPOS compared to storing
			*   the base offset directly. */
  HASHNUM_t npos;       /**< size of the array pos */
  size_t npos_alloc;   /**< allocated memory for array pos (as number of elements) */
  HASHPOS_t maxpos;     /**< Maximum position in array pos + 1*/
  HALFWORD_t *wordidx;/**< Imperfect part of the hashed words to resolve collisions */
  HASHNUM_t *posidx;    /**< Index points to pos for each element in wordidx */
  WORDNUM_t nwords;     /**< Number of different words */
  WORDNUM_t *wordctr;    /**< array of size nkeys, counting the number of different words for each key 
			* (only used during hash index construction) */
};

typedef int (* KEYWORKERF_t) (const HashTable *, HASHWORD_t, HASHPOS_t);
typedef int (* KEYCHECKERF_t) (const HashTable *, const SeqSet *, HASHWORD_t, SEQNUM_t, SETSIZ_t);

/******************************************************************************
 *********************************** Macros ***********************************
 ******************************************************************************/

#define MAKE_HASHKEY(word, word_hi, key_hi, key)\
               (word_hi) = ((word) & htp->wordmask_hi) >> htp->nbits_lo; \
               (key_hi) = hash32mix(word_hi) % htp->keymod; \
               (key) = ((key_hi) << htp->nbits_lo) + ((word) & htp->wordmask_lo);

/******************************************************************************
 ******************************* Private Methods ******************************
 ******************************************************************************/
static uint32_t hash32mix(uint32_t a)
{
   a = (a+0x7ed55d16) + (a<<12);
   a = (a^0xc761c23c) ^ (a>>19);
   a = (a+0x165667b1) + (a<<5);
   a = (a+0xd3a2646c) ^ (a<<9);
   a = (a+0xfd7046c5) + (a<<3);
   a = (a^0xb55a4f09) ^ (a>>16);
   return a;
}

static char *printWord(char *buff, HASHWORD_t word, uint8_t len)
{
  const char alphabet[]="ACGT";
  uint8_t i;

  for (i=len; i>0; i--) {
    buff[i-1]= alphabet[word & SEQCOD_STDNT_MASK];
    word = word>>2;
  }
  buff[len] = '\0';
  return buff;
}

/******************************************************************************
 ********************* Private function type KEYWORKERF_t *********************
 ******************************************************************************/

static int countPerfectKey(const HashTable *htp, HASHWORD_t word, HASHPOS_t tuplectr)
{
  htp->idx[word & htp->wordmask]++;
  return ERRCODE_SUCCESS;
}

static int setKeyPos(const HashTable *htp, HASHWORD_t word, HASHPOS_t tuplectr)
{
  HASHKEY_t k = word & htp->wordmask;
  if (htp->idx[k+1] - htp->idx[k] > 0) { /* were there any counts? */
    htp->pos[htp->idx[k]++] = tuplectr;
  }
  return ERRCODE_SUCCESS;
}

static int countCollisionKey(const HashTable *htp, HASHWORD_t word, HASHPOS_t tuplectr)
{
  HASHKEY_t key, key_hi;
  uint32_t word_hi;
  
  MAKE_HASHKEY(word, word_hi, key_hi, key);
#ifdef hashidx_debug
  if (key == 0)
    printf("hasidx_debug::countCollisionKey:[%u] word = %llu, tuplectr = %u\n", 
	   htp->idx[key], word&htp->wordmask, tuplectr);
#endif

  htp->idx[key]++;

  return ERRCODE_SUCCESS;
}

static int collectWords(const HashTable *htp, HASHWORD_t word, HASHPOS_t tuplectr)
{
  HASHKEY_t key, key_hi;
  HALFWORD_t word_hi;
  HALFWORD_t *iwp;

  MAKE_HASHKEY(word, word_hi, key_hi, key);
#ifdef hashidx_debug
  if (key == 0)
    printf("hashidx_debug::collectWords:key = %u\n", key);
#endif
  
  iwp = htp->pos + htp->idx[key];

  if ( 0 == htp->wordctr[key] ) {
    htp->wordctr[key] = 1;
    iwp[0] = word_hi;
  } else {
    WORDNUM_t pivot, b, a = 0;
    b = htp->wordctr[key]-1;
    while (a < b) {
      pivot = (a + b) >> 1;
      if (iwp[pivot] < word_hi) {
	/* upper half */
	a = pivot + 1;
      } else {
	/* lower half */
	b = pivot;
      }
    }
    if (a > b || iwp[b] < word_hi) {
      iwp[htp->wordctr[key]++] = word_hi;
    } else if (iwp[b] > word_hi) {
      for (a = htp->wordctr[key]++; a > b; a--)
	iwp[a] = iwp[a-1];
      iwp[b] = word_hi;
    }
  }

  return ERRCODE_SUCCESS;
}

static int countWords(const HashTable *htp, HASHWORD_t word, HASHPOS_t tuplectr)
{
  HASHKEY_t key, key_hi;
  HALFWORD_t word_hi;
  HASHNUM_t pivot, b, a;

  MAKE_HASHKEY(word, word_hi, key_hi, key);
  b = htp->idx[key+1];
  if (b < 1)
    return ERRCODE_NOHASHWORD;
  a = htp->idx[key];
  b--;
  while (a < b) {
    pivot = (a + b) >> 1;
    if (htp->wordidx[pivot] < word_hi) {
      /* upper half */
      a = pivot + 1;
    } else {
      /* lower half */
      b = pivot;
    }
  }
  if (a > b || htp->wordidx[b] != word_hi) 
    return ERRCODE_NOHASHWORD;

  htp->posidx[b]++;
  
  return ERRCODE_SUCCESS;
}

static int setWordPos(const HashTable *htp, HASHWORD_t word, HASHPOS_t tuplectr)
{
  HASHKEY_t key, key_hi;
  HALFWORD_t word_hi;
  HASHNUM_t pivot, b, a;

  MAKE_HASHKEY(word, word_hi, key_hi, key);
  b = htp->idx[key+1];
  if (b < 1)
    return ERRCODE_NOHASHWORD;
  a = htp->idx[key];
  b--;
  while (a < b) {
    pivot = (a + b) >> 1;
    if (htp->wordidx[pivot] < word_hi) {
      /* upper half */
      a = pivot + 1;
    } else {
      /* lower half */
      b = pivot;
    }
  }
  if (a > b || htp->wordidx[b] != word_hi) 
    return ERRCODE_NOHASHWORD;

  htp->pos[htp->posidx[b]++] = tuplectr;
  
  return ERRCODE_SUCCESS;
}

static HASHPOS_t calcKtupOffs(UCHAR_t *ktupoffs, SEQLEN_t offs, SEQNUM_t sx, UCHAR_t nskip, const SeqSet *ssp)
{
  SETSIZ_t soffs;
  HASHPOS_t tc;
  seqSetGetSeqDatByIndex(&soffs, NULL, sx, ssp);
  tc = (soffs + offs)/nskip;
  if (ktupoffs) *ktupoffs = soffs + offs - tc*nskip;
  return tc;
}

/******************************************************************************
 ********************* Private function type KEYCHECKERF_t ********************
 ******************************************************************************/

static int checkPerfectHash(const HashTable *htp, const SeqSet *ssp, HASHWORD_t word, 
			    SEQNUM_t seqidx, SETSIZ_t basidx)
{
  int errcode = ERRCODE_SUCCESS;
  unsigned char flag_found = 0;
  SEQNUM_t si = 0;
  HASHKEY_t key = word & htp->wordmask;
  HASHPOS_t *posp = htp->pos + htp->idx[key];
  HASHPOS_t *endposp = htp->pos + htp->idx[key+1];
  SEQLEN_t offs;
  SETSIZ_t po;
  
  for (; posp<endposp; posp++) {
    po = *posp * htp->nskip;
    if ((errcode = seqSetGetIndexAndOffset(&si, &offs, po, ssp)))
      return errcode;
    if (si == seqidx && po == offs + basidx) {
      flag_found = 1;
      break;
    }
  }
  if (!flag_found) {
    printf("HASH at seqidx %lli, basidx %lu failed with key %u, pos %u, seqidx %lli, offset %u\n",
	   (signed long long) seqidx, (unsigned long) basidx, key, *posp, 
	   (signed long long) si, (unsigned int) offs);
    errcode = ERRCODE_BROKENHASH;
  }

  return errcode;
}

static int checkHashWithCollisions(const HashTable *htp, const SeqSet *ssp, HASHWORD_t word, 
				   SEQNUM_t seqidx, SETSIZ_t basidx)
{
  int errcode = ERRCODE_SUCCESS;
  unsigned char flag_found = 0;
  HASHKEY_t key, key_hi;
  HALFWORD_t word_hi;
  HASHNUM_t pivot, b, a;
  HASHNUM_t i, i_end;
  HASHPOS_t *posp;
  SETSIZ_t po;
  SEQLEN_t offs;
  SEQNUM_t si = 0;

  MAKE_HASHKEY(word, word_hi, key_hi, key);
  b = htp->idx[key+1];
  if (b < 1)
    return ERRCODE_NOHASHWORD;
  a = htp->idx[key];
  b--;
  while (a < b) {
    pivot = (a + b) >> 1;
    if (htp->wordidx[pivot] < word_hi) {
      /* upper half */
      a = pivot + 1;
    } else {
      /* lower half */
      b = pivot;
    }
  }
  if (a > b || htp->wordidx[b] != word_hi) 
    return ERRCODE_NOHASHWORD;

  posp = htp->pos + htp->posidx[b];
  if (htp->posidx[b+1] < htp->posidx[b])
    return ERRCODE_BROKENHASH;
  i_end = htp->posidx[b+1] - htp->posidx[b];
  for (i=0; i < i_end; i++) {
    po = posp[i]*htp->nskip;
    if ((errcode = seqSetGetIndexAndOffset(&si, &offs, po, ssp)))
      return errcode;
    if (si == seqidx && po == offs + basidx) {
      flag_found = 1;
      break;
    }
  }
  
  if (!flag_found) {
    printf("HASH at seqidx %lli, basidx %lu failed with key %u, pos %u, seqidx %llii, offset %u\n",
	   (signed long long) seqidx, (unsigned long) basidx, key, *posp, 
	   (signed long long) si, (unsigned int) offs);
    errcode = ERRCODE_BROKENHASH;
  }

  return errcode;
}
/******************************************************************************
 ********************* Private Methods of Type HashTable **********************
 ******************************************************************************/
static int calcMaxWordNumPerKey(uint32_t *nkey_zero, const HashTable *htp)
{
  HASHKEY_t k;
  int nw, maxnw = 0;
  uint32_t nzero = 0;
  for (k=0; k<htp->nkeys; k++) {
    nw = htp->idx[k+1] - htp->idx[k];
    if (nw > maxnw)
      maxnw = nw;
    else if (nw == 0)
      nzero++;
  }
  if (nkey_zero)
    *nkey_zero = nzero;
  return maxnw;
}

static WORDNUM_t getWordHisto(WORDNUM_t *counts, HASHNUM_t n, const HashTable *htp)
{
  HASHNUM_t j;
  HASHKEY_t i;
  HASHNUM_t nw;
  uint32_t nout = 0;
  for (j=0; j<n; j++)
    counts[j] = 0;
  for (i=0; i<htp->nkeys; i++) {
    nw = htp->idx[i+1] - htp->idx[i];
    if (nw < n) {
      counts[nw]++;
    } else {
      (nout)++;
    }
  }
  return nout;
}

static int doWordsInSeq(HASHPOS_t *tuplectr,
			UCHAR_t *ktup_offs,
			const SeqFastq *sqp,
			const HashTable *htp,
			KEYWORKERF_t keyfunc)
     /**< Count the number of k-tuples for each key in sequence i of the hash table 
      * (do_counting = 1) or write the k-tuple positions into the hash table body
      * and set hash table header accordingly (do_counting == 0).
      * On return, the hash table header contains the updated k-tuple counts.
      * Sequences have to be of code SEQCOD_MANGLED, otherwise ERRCODE_SEQCODE is returned.
      *
      * \param htp      Pointer to hash table. Hash table is updated with the k-tuple counts.
      * \param tuplectr Holds on return the total number of k-tuples in sequence i added
      *                 to the value on input.
      * \param ktup_offs Offset of the first k-tuple in sequence. Returns the offset of 
      *                 first k-tuple in the next sequence calculated as if the two sequences
      *                 were joined together with an extra terminator base in between.
      * \param sqp      Sequence from which k-tuples are read.
      * \param keyfunc  Function to be performed on key.
      *
      * \note: Kmer-words with non-standard nucleotides are discarded.
      */
{
  char codtyp;
  const char *datap;
  UCHAR_t ktup_i, non_stdnt;
  SEQLEN_t seqlen;
  HASHWORD_t word;

  if (*ktup_offs >= htp->nskip)
    return ERRCODE_ARGRANGE;
  datap = seqFastqGetConstSequence(sqp, &seqlen, &codtyp);
  if (codtyp != SEQCOD_MANGLED) 
    return ERRCODE_SEQCODE;
  if (seqlen < htp->wordlen) 
    return ERRCODE_SHORTSEQ; /* sequence too short to be hashed */

  word = 0;
  ktup_i = htp->wordlen + *ktup_offs; /* number of bases to be filled in k */
  non_stdnt= 0;
  while (*datap) {
    if ((*datap)&SEQCOD_STDNT_TESTBIT) { /* do not hash non-standard nucleotides */
      non_stdnt = htp->wordlen;
    } else if ((non_stdnt)) {
      non_stdnt--;
    }
    word = (word<<2) + ((*datap++) & SEQCOD_STDNT_MASK);
    if (--ktup_i > 0) continue;
    if (!(non_stdnt)) {
      (*keyfunc) (htp, word, *tuplectr);
    }
    (*tuplectr)++;
    ktup_i = htp->nskip;
  }

  /* at this stage: total concatenated sequence length 
   * l = (*tuplectr - 1)*htp->nskip + htp->nskip - ktup_i + htp->wordlen
   *   = (*tuplectr)*htp->nskip - ktup_i + htp->wordlen
   */
  *ktup_offs = (htp->wordlen - ktup_i)%htp->nskip;

  if ((*ktup_offs)) *ktup_offs = htp->nskip - *ktup_offs;

  *tuplectr += (htp->wordlen - ktup_i + *ktup_offs)/htp->nskip;

  return ERRCODE_SUCCESS;
}

static int doAllWordsInSeqSet(HASHPOS_t *tuplectrp, SeqFastq *sqbufp, 
			      const HashTable *htp, 
			      const SeqSet *ssp, 
			      const InterVal *ivp,
			      const SeqCodec *codecp,
			      KEYWORKERF_t keyfunc)
{
  int errcode = ERRCODE_SUCCESS;
  char cod;
  UCHAR_t tupleoffs=0;
  SEQNUM_t i, n_seq = seqSetGetOffsets(ssp, NULL);
  SEQLEN_t seqlen;
  if (n_seq<1) return ERRCODE_NOSEQ;

  *tuplectrp = 0;

  if (ivp) {
    int n, ivn = interValNum(ivp);
    SEQLEN_t lo, hi, sl;
    SEQNUM_t sx;
#ifdef hashidx_debug_init
    fprintf(stderr, "hashidx.c::doAllWordsInSeqSet: inv = %i\n", ivn);
#endif
    for (n=0; n<ivn && !errcode; n++) {
      if ((errcode = interValGet(&lo, &hi, &sx, NULL, (int) n, ivp)))
	return errcode;
      *tuplectrp = calcKtupOffs(&tupleoffs, lo, sx, htp->nskip, ssp);

      sl = hi - lo + 1;
      if (sl < htp->wordlen)
	continue;
      if ((errcode = seqSetFetchSegmentBySequence(sqbufp, sx, lo, sl, ssp, codecp)))
	return errcode;
      seqFastqGetConstSequence(sqbufp, &seqlen, &cod);
      if (seqlen != sl)
	return ERRCODE_ASSERT;
      if (seqlen > INT_MAX) 
	return ERRCODE_SEQLEN;
      if (cod == SEQCOD_ASCII &&
	  (errcode = seqFastqEncode(sqbufp, codecp)))
	return errcode;
      errcode = doWordsInSeq(tuplectrp, &tupleoffs, sqbufp, htp, keyfunc);
    }
  } else {
#ifdef hashidx_debug
    UCHAR_t to;
    HASHPOS_t tc;
#endif
#ifdef hashidx_debug_init
    fprintf(stderr, "hashidx.c::doAllWordsInSeqSet: no intervals\n");
#endif
    for (i=0; i<n_seq && !errcode; i++) {
      //if (verbose) printf("... in sequence %i ...\n", i);
      if ((errcode = seqSetFetchSegmentBySequence(sqbufp, i, 0, 0, ssp, codecp)))
	return errcode;
      seqFastqGetConstSequence(sqbufp, &seqlen, &cod);
      if (cod == SEQCOD_ASCII &&
	  (errcode = seqFastqEncode(sqbufp, codecp)))
	return errcode;
      if (seqlen > INT_MAX) 
	errcode = ERRCODE_SEQLEN;
#ifdef hashidx_debug
      tc = calcKtupOffs(&to, 0, i, htp->nskip, ssp);
      if (tc != *tuplectrp || to != tupleoffs)
	errcode = ERRCODE_ASSERT;
#endif
      if (!errcode)
	errcode = doWordsInSeq(tuplectrp, &tupleoffs, sqbufp, htp, keyfunc);
    }
  }
  return errcode;
}

static int checkWordsInSeqSet(SeqFastq *sqbufp, const HashTable *htp, 
			      const SeqSet *ssp, const SeqCodec *codecp,
			      KEYCHECKERF_t checkfunc)
{
  int errcode = ERRCODE_SUCCESS;
  char codtyp;
  const char *datap;
  uint8_t ktup_i, n_nonstdnt;
  SEQLEN_t seqlen;
  SETSIZ_t basidx;
  SEQNUM_t seqidx, n_seq;
  HASHWORD_t word;

  n_seq = seqSetGetOffsets(ssp, NULL);
  ktup_i = 0;
  for (seqidx=0; seqidx<n_seq && !(errcode); seqidx++) {
    printf("Checking sequence %lli ...\n", (signed long long) seqidx);
    if ((errcode = seqSetFetchSegmentBySequence(sqbufp, seqidx, 0, 0, ssp, codecp)))
      return errcode;
    datap = seqFastqGetConstSequence(sqbufp, &seqlen, &codtyp);
    
    if (codtyp == SEQCOD_ASCII) {
      if ((errcode = seqFastqEncode(sqbufp, codecp)))
	return errcode;
    } else if (codtyp != SEQCOD_MANGLED) 
      return ERRCODE_SEQCODE;
    if (seqlen < htp->wordlen) return ERRCODE_SHORTSEQ; /* sequence too short to be hashed */

    word = 0LL;
    basidx = ktup_i;
    ktup_i += htp->wordlen; /* number of bases to be filled in k */
    n_nonstdnt = 0;
    while ((*datap) && !(errcode)) {
      if ((*datap)&SEQCOD_STDNT_TESTBIT) { /* do not hash non-standard nucleotides */
	n_nonstdnt = htp->wordlen;
      } else if ((n_nonstdnt)) {
	n_nonstdnt--;
      }

      word = ((word<<2) + ((*datap++) & SEQCOD_STDNT_MASK));
      if (--ktup_i > 0) continue;
      if (!(n_nonstdnt))
	  errcode = (*checkfunc)(htp, ssp, word, seqidx, basidx);
      ktup_i = htp->nskip;
      basidx += htp->nskip;
    }
    ktup_i = (htp->wordlen - ktup_i)%htp->nskip;
    if ((ktup_i)) ktup_i = htp->nskip - ktup_i;
  }
  return errcode;
}

static int checkQuickPerfectHashIndex(SeqFastq *sqbufp, const HashTable *htp, 
				      const SeqSet *ssp, const SeqCodec *codecp)
{
  int errcode;
  char codtyp;
  UCHAR_t i;
  const char *datap, *cp;
  SEQLEN_t seqlen;
  SETSIZ_t offset, offset_end;
  HASHKEY_t key;
  HASHPOS_t *p, *endp;

  for (key=0; key<htp->nkeys; key++) {
    endp = htp->pos + htp->idx[key+1];
    for (p=htp->pos+htp->idx[key]; p<endp;p++) {
      offset = (*p)*htp->nskip;
      offset_end = offset +  htp->wordlen - 1;
      if ((errcode = seqSetFetchSegment(sqbufp, &offset, &offset_end, ssp, codecp)))
	return errcode;
      datap = seqFastqGetConstSequence(sqbufp, &seqlen, &codtyp);
      if (codtyp == SEQCOD_ASCII) {
	if ((errcode = seqFastqEncode(sqbufp, codecp)))
	  return errcode;
	datap = seqFastqGetConstSequence(sqbufp, &seqlen, &codtyp);
      }
      if (!datap) return ERRCODE_FAILURE;
      if (codtyp != SEQCOD_MANGLED) return ERRCODE_SEQCODE;
      for(i = htp->wordlen, cp = datap; i>0 && *cp;) {
	if (((*cp++) & SEQCOD_STDNT_MASK) != ((key>>((--i)<<1)) & SEQCOD_STDNT_MASK)) {
	  printf("HASH (perfect) failed with key %d, pos %u, offset %lu\n",
		 key, *p, (unsigned long) offset);
	  return ERRCODE_BROKENHASH;
	}
      }
    }
  }

  return ERRCODE_SUCCESS;
}

static int checkQuickHashIndexWithCollisions(SeqFastq *sqbufp, const HashTable *htp, 
					     const SeqSet *ssp, const SeqCodec *codecp)
{
  int errcode;
  char codtyp, word[MAXTUPLEN+1];
  UCHAR_t i;
  const char *pData, *cp;
  SEQLEN_t seqlen;
  HASHNUM_t iha;
  HASHPOS_t ipos;
  SETSIZ_t offset, offset_end;
  HASHKEY_t key;
  HASHWORD_t keyword;

  for (key=0; key<htp->nkeys; key++) {
    for (iha = htp->idx[key]; iha < htp->idx[key+1]; iha++) {
      keyword = (((HASHWORD_t) htp->wordidx[iha])<< htp->nbits_lo) + (key & htp->wordmask_lo);
      for (ipos = htp->posidx[iha]; ipos < htp->posidx[iha+1]; ipos++) { 
	offset = htp->pos[ipos]*htp->nskip;
	offset_end = offset +  htp->wordlen - 1;
	if ((errcode = seqSetFetchSegment(sqbufp, &offset, &offset_end, ssp, codecp)))
	  return errcode;
	pData = seqFastqGetConstSequence(sqbufp, &seqlen, &codtyp);
	if (codtyp == SEQCOD_ASCII) {
	  if ((errcode = seqFastqEncode(sqbufp, codecp)))
	    return errcode;
	  pData = seqFastqGetConstSequence(sqbufp, &seqlen, &codtyp);
	}

	if (!pData) return ERRCODE_FAILURE;
	if (codtyp != SEQCOD_MANGLED) 
	  return ERRCODE_SEQCODE;

	for(i = htp->wordlen, cp = pData; i>0 && *cp;) {
	  if (((*cp++) & SEQCOD_STDNT_MASK) != ((keyword>>((--i)<<1)) & SEQCOD_STDNT_MASK)) {
	    seqFastqDecode(sqbufp, codecp);
	    pData = seqFastqGetConstSequence(sqbufp, &seqlen, &codtyp);
	    printWord(word, keyword, htp->wordlen);
	    printf("HASH (collisions) %d failed with word %s (refseq %s), pos %u, offset %llu\n",
		   key, word, pData, ipos, (long long unsigned int) offset);
	    return ERRCODE_BROKENHASH;
	  }
	}
      }
    }
  }
  return ERRCODE_SUCCESS;
}

/******************************************************************************
 ********************** Public Methods of Type HashTable **********************
 ******************************************************************************/

HashTable *hashTableCreate(UCHAR_t wordlen, UCHAR_t nskip, 
			   UCHAR_t nbits_key, UCHAR_t nbits_perf,
			   UCHAR_t typ)
{
  HashTable *htp;
  
  if (wordlen > MAXTUPLEN)
    return NULL;

  if (typ == HASHIDXTYP_PERFECT) {
    nbits_key = (wordlen<<1);
    nbits_perf = 0;
  }

  if(wordlen > ((MAX_NBITS_KEY + nbits_perf)/2) ||
      nbits_key > MAX_NBITS_KEY ||
      nbits_key <= nbits_perf)
    return NULL;

  EMALLOCP0(htp);
  if (!htp) return NULL;

  if (typ == HASHIDXTYP_PERFECT)
    nbits_key = (wordlen<<1);
  htp->typ = typ;
  htp->status = HASHSETUP_EMPTY;
  htp->wordlen = wordlen;
  htp->nskip = nskip;
  htp->nbits_key = nbits_key;
  htp->nbits_lo = nbits_perf;
  htp->nkeys = ((HASHKEY_t) 1) << nbits_key; /* number of keys */
  htp->keymask = htp->nkeys - 1;
  htp->keymask_lo = (((HASHKEY_t) 1) << nbits_perf)-1;
  htp->wordmask = (((HASHWORD_t) 1) << (wordlen<<1)) - 1;
  htp->wordmask_lo = (((HASHWORD_t) 1) << nbits_perf) - 1;
  htp->wordmask_hi = (~htp->wordmask_lo) & htp->wordmask;
  htp->keymod = ((HASHWORD_t) 1) << (nbits_key - nbits_perf);

  htp->pos = NULL;
  htp->npos = htp->npos_alloc = 0;
  htp->wordidx = NULL;
  htp->wordctr = NULL;

  ECALLOCP(htp->nkeys+ARRAY_MARGIN, htp->idx); /* allocate 4^k+2 */
  if (!htp->idx) {
    hashTableDelete(htp);
    return NULL;
  }
  return htp;
}
  
void hashTableDelete(HashTable *htp)
{
  if (htp) {
    free(htp->idx);
    free(htp->pos);
    free(htp->wordidx);
    free(htp->wordctr);
  }
  free(htp);
}

void hashTableReset(HashTable *htp, UCHAR_t nskip)
{
  if (nskip > 0)
    htp->nskip = nskip;
  if (htp->status != HASHSETUP_EMPTY) {
    memset(htp->idx, 0, (htp->nkeys+ARRAY_MARGIN)*sizeof(uint32_t));
#ifdef hashidx_debug_init
    fprintf(stderr, "hashidx.c: memset(htp->idx,0,%llu) : htp->idx[%llu] = %lu\n",
	    (unsigned long long) (htp->nkeys+ARRAY_MARGIN),
	    (unsigned long long) (htp->nkeys+1),
	    (unsigned long) htp->idx[htp->nkeys + 1]);
#endif
    htp->status = HASHSETUP_EMPTY;
  }
}

int hashTableSetUp(HashTable *htp, SeqFastq *sqbufp, const SeqSet *ssp,
		   const InterVal *ivp,
		   const SeqCodec *codecp, 
		   uint32_t *npos_max,
		   char verbose)
{
  int errcode = ERRCODE_SUCCESS;
  HASHKEY_t w;
  HASHNUM_t i, i_end;
  size_t k;
  uint32_t tuplectr=0;
  size_t npos;

  if (htp->status == HASHSETUP_EMPTY) {
    SETSIZ_t totlen;
    seqSetGetSeqNumAndTotLen(&totlen, ssp);
    if ( htp->nskip < 1 || (totlen+1)/htp->nskip > HASHPOS_MAX )
      return ERRCODE_ASSERT;
#ifdef hashidx_debug_init
    fprintf(stderr, "hashidx.c: htp->status == HASHSETUP_EMPTY\n");
#endif
  } else {
    hashTableReset(htp, 0);
  }
  if (verbose) {
    if (htp->typ == HASHIDXTYP_PERFECT)
      fprintf(stderr, "# Counting k-tuple occurrences ...\n");
    else
      fprintf(stderr, "# Counting keys ...\n");
  }
  htp->idx += 2;

  errcode = doAllWordsInSeqSet(&tuplectr, sqbufp, htp, ssp, ivp, 
			       codecp, 
			       (htp->typ == HASHIDXTYP_PERFECT)? 
			       countPerfectKey: countCollisionKey);
  htp->idx -= 2;
  if (errcode)
    return errcode;

  /* find out the total number of key counts */
  htp->idx[0] = htp->idx[1] = 0;
  npos = 0;
  for (k=2; k<=htp->nkeys+1; k++) {
    npos +=  htp->idx[k];
    htp->idx[k] = (uint32_t) npos;
  }
  if (npos > UINT32_MAX)
    return ERRCODE_KPOSOFLO;

  /* at this point, htp->idx[k+2] - htp->idx[k+1] is the number of occurences of key k */

  if (npos_max != NULL) {
    if (*npos_max > 0 && npos > *npos_max) {
      *npos_max = (uint32_t) npos;
      return ERRCODE_MAXKPOS;
    }
    *npos_max = (uint32_t) npos;
  }

  /* allocate memory for sequence positions accordingly */
  if (htp->pos == NULL || npos > htp->npos_alloc) {
    size_t n_alloc = (npos+1)/BLKSIZ_IDXPOS + 1;
    n_alloc *= BLKSIZ_IDXPOS;

    if (htp->pos == NULL) {
      ECALLOCP(n_alloc, htp->pos);
      if (!htp->pos) return ERRCODE_NOMEM;
    } else {    
      void *hp = EREALLOCP(htp->pos, n_alloc);
      if (!hp) return ERRCODE_NOMEM;
      htp->pos = hp;
    }
    htp->npos_alloc = n_alloc;
  }
  htp->npos = (uint64_t) npos;

  if (htp->typ == HASHIDXTYP_PERFECT) {
    /* fill in the key-positions */
    if(verbose) 
      fprintf(stderr,"# Setting the k-tuple positions in index ...\n");
    htp->idx++;
    if ((errcode = doAllWordsInSeqSet(&tuplectr, sqbufp, htp, ssp, ivp, codecp, setKeyPos)))
      return errcode;
    htp->idx--;
    htp->idx[0] = 0;
    /* at this point, htp->idx[k+1] - htp->idx[k] is the number of keys k */
  } else {
    
    free(htp->wordctr);
    ECALLOCP(htp->nkeys, htp->wordctr);
    if (!htp->wordctr) 
      return ERRCODE_NOMEM;
    if(verbose) 
      fprintf(stderr, "# Counting k-tuples ...\n");
    htp->idx++;
    errcode = doAllWordsInSeqSet(&tuplectr, sqbufp, htp, ssp, ivp, codecp, collectWords);
    htp->idx--;
    if (errcode)
      return errcode;
    /* at this point, htp->wordctr[k] is the number of different words with the same key k,
     * the words are htp->pos[htp->idx[k+1]], ..., htp->pos[htp->idx[k+1]+htp->wordctr[k]-1] 
     * sorted in ascending order.
     */

    if(verbose) 
      fprintf(stderr, "# Allocating k-tuple arrays ...\n");
    for (w=0, k=0; k<htp->nkeys; k++)
      w += htp->wordctr[k];
    htp->nwords = w;

    free(htp->wordidx);
    if ((ECALLOCP(2*(w+ARRAY_MARGIN), htp->wordidx)) == NULL)
      return ERRCODE_NOMEM;
    htp->posidx = htp->wordidx + w + ARRAY_MARGIN;

    w = 0;
    htp->idx[0] = 0;
    for (k=0; k<htp->nkeys; k++) {
      i = htp->idx[k+1];
      i_end = i + htp->wordctr[k];
      htp->idx[k+1] = htp->idx[k] + htp->wordctr[k];
      for (; i<i_end; i++, w++)
	htp->wordidx[w] = htp->pos[i];
    }
    if (w != htp->nwords)
      return ERRCODE_ASSERT;

    free(htp->wordctr);
    htp->wordctr = NULL;
    
    /* at this point, w is the total number of different words, htp->idx[k+1] - htp->idx[k]
     * is the number of different words with the same key htp->wordidx[htp->idx[k]], ..., 
     * htp->wordidx[htp->idx[k+1]-1] */

    if(verbose) 
      fprintf(stderr, "# Counting k-tuple occurrences ...\n");
    htp->posidx += 2;
    errcode = doAllWordsInSeqSet(&tuplectr, sqbufp, htp, ssp, ivp, codecp, countWords);
    htp->posidx -= 2;
    if (errcode)
      return errcode;
    /* at this point htp->wordidx[htp->idx[k]], ..., htp->wordidx[htp->idx[k+1]-1] are the words
     * with key k, and htp->posidx[htp->idx[k]+j+2] is the number of occurrences of word j with key k 
     * where 0 <= j < htp->posidx[htp->idx[k] + 2]. */
    htp->posidx[0] = htp->posidx[1] = 0;
    for (w=2; w<htp->nwords; w++)
      htp->posidx[w+1] += htp->posidx[w];
    /* at this point htp->wordidx[htp->idx[k]], ..., htp->wordidx[htp->idx[k+1]-1] are the words
     * with key k, and htp->wordidx[htp->idx[k]+j+2] - htp->wordidx[htp->idx[k]+j+1] is the number
     * of occurrences of word j with key k where 0 <= j < htp->idx[k+1] - htp->idx[k]. */
    if(verbose) 
      fprintf(stderr, "# Setting the k-tuple positions in index ...\n");
    htp->posidx++;
    errcode = doAllWordsInSeqSet(&tuplectr, sqbufp, htp, ssp, ivp, codecp, setWordPos);
    htp->posidx--;
    if (errcode)
      return errcode;
    /* at this point htp->wordidx[htp->idx[k]], ..., htp->wordidx[htp->idx[k+1]-1] are the words
     * with key k, and 
     * htp->pos[htp->posidx[htp->idx[k]+j]], ..., htp->pos[htp->posidx[htp->idx[k]+j+1]-1] are the
     * positions of word j with key k where 0 <= j < htp->idx[k+1] - htp->idx[k]. */
  }
  htp->maxpos = (tuplectr > 0)? tuplectr-1: 0;
  htp->status = HASHSETUP_COMPLETE;  /* flag that no more sequences can be added */

  if (verbose) 
    fprintf(stderr,"# Hash table is set up.\n");
  return errcode;
  }

int hashTableCheckExtensive(SeqFastq *sqbufp, const HashTable *htp, 
			    const SeqSet *ssp, const SeqCodec *codecp)
{
  return (htp->typ == HASHIDXTYP_PERFECT)?
    checkWordsInSeqSet(sqbufp, htp, ssp, codecp, checkPerfectHash):
    checkWordsInSeqSet(sqbufp, htp, ssp, codecp, checkHashWithCollisions);
}

int hashTableCheckQuick(SeqFastq *sqbufp, const HashTable *htp, 
			const SeqSet *ssp, const SeqCodec *codecp)
{
  if (htp->status != HASHSETUP_COMPLETE)
    return ERRCODE_FAILURE;
  return (htp->typ == HASHIDXTYP_PERFECT)?
    checkQuickPerfectHashIndex(sqbufp, htp, ssp, codecp):
    checkQuickHashIndexWithCollisions(sqbufp, htp, ssp, codecp);
}

UCHAR_t hashTableGetKtupLen(const HashTable *htp, UCHAR_t *nskip)
{
  if (nskip)
    *nskip = htp->nskip;
  return htp->wordlen;
}

HASHPOS_t hashTableGetMaxPos(const HashTable *htp)
{
  return htp->maxpos;
}

void hashTablePrintStats(FILE *fp, const HashTable *htp)
{
  WORDNUM_t maxperkey = 0;
  HASHKEY_t nkey_zero;
#ifdef hashidx_debug
  HASHKEY_t i;
  char buf[MAXTUPLEN+1];
  WORDNUM_t *countsp = NULL;
#endif

  fprintf(fp, "# =-=-=-=-= Hash Index Stats =-=-=-=-=\n");
  if (htp->typ == HASHIDXTYP_PERFECT) {
    fprintf(fp, "# Perfect hash index.\n");
  } else {
    fprintf(fp, "# Hash index with collisions.\n");
  }
  fprintf(fp, "# Word length:              %-d bases\n", htp->wordlen);
  fprintf(fp, "# Skip step:                %-d bases\n", htp->nskip);
  fprintf(fp, "# Number of hash keys:      %u\n", htp->nkeys);
  fprintf(fp, "# Number of word positions: %u\n", htp->npos);
  if (htp->typ == HASHIDXTYP_HASH32MIX) {
    maxperkey = calcMaxWordNumPerKey(&nkey_zero, htp);
    fprintf(fp, "# Number of different words:        %u\n", htp->nwords);
    fprintf(fp, "# Maximum number of words per key:  %i\n", maxperkey);
    fprintf(fp, "# Number of keys without words:     %u\n", nkey_zero);
#ifdef hashidx_debug
    fprintf(fp, "Size of idx:     %-g MB\n", ((double) (htp->nkeys+1))*4/MEGABYT);
    fprintf(fp, "Size of wordidx: %-g MB\n", ((double) (htp->nwords+1))*8/MEGABYT);
    fprintf(fp, "Size of pos:     %-g MB\n", ((double) (htp->npos))*4/MEGABYT);
#endif
  }
  fprintf(fp, "# =-=-= End of Hash Index Stats =-=-=\n");

#ifdef hashidx_debug
  if (htp->typ == HASHIDXTYP_PERFECT) {
    for (i=0; i<htp->nkeys; i++) {
      if (htp->idx[i+1] > htp->idx[i]) {
	printf("[%u] %s %u\n", i, printWord(buf, i, htp->wordlen),
	       htp->idx[i+1] - htp->idx[i]);
      }
    }
  } else if (maxperkey + 1 < INT_MAX) {
    ECALLOCP(maxperkey+1, countsp);
    if (countsp != NULL) {
      WORDNUM_t j, nout = getWordHisto(countsp, maxperkey + 1, htp);
      for (j=0; j<maxperkey+1; j++) {
	if (countsp[j] > 0)
	  printf("[%i] %u\n", j, countsp[j]);
      }
      printf("Number of keys with > %i words: %u\n", maxperkey, nout);
    }
    free(countsp);
  }
#endif
}


int hashTableCmp(const HashTable *ap, const HashTable *bp)
{
  uint32_t i;

  if (ap->nkeys != bp->nkeys) {
    printf("Tables have different key length: (A) %d, (b) %d\n",
	   ap->wordlen, bp->wordlen);
    return ERRCODE_FAILURE;
  }

  if (ap->nskip != bp->nskip) {
    printf("Tables have different skip step size: (A) %d, (b) %d\n",
	   ap->nskip, bp->nskip);
    return ERRCODE_FAILURE;
  }

  if (ap->typ != bp->typ ||
      ap->status != bp->status) {
    printf("Tables differ in type or status.\n");
    return ERRCODE_FAILURE;
  }

  for (i=0; i<ap->nkeys; i++) {
    if (ap->idx[i] != bp->idx[i])
      return ERRCODE_FAILURE;
  }
  for (i=0; i<ap->npos; i++) {
    if (ap->pos[i] != bp->pos[i])
      return ERRCODE_FAILURE;
  }

  if (ap->typ ==  HASHIDXTYP_PERFECT)
    return ERRCODE_SUCCESS;

  if (ap->nbits_key != bp->nbits_key ||
      ap->nbits_lo != bp->nbits_lo ||
      ap->keymask != bp->keymask ||
      ap->keymask_lo != bp->keymask_lo ||
      ap->keymask_hi != bp->keymask_hi ||
      ap->keymod != bp->keymod ||
      ap->wordmask != bp->wordmask ||
      ap->wordmask_lo != bp->wordmask_lo ||
      ap->wordmask_hi != bp->wordmask_hi)
    return ERRCODE_FAILURE;

  if (ap->nwords != bp->nwords) {
    printf("Tables differ in the number of words: (A) %u, (B) %u.\n", ap->nwords, bp->nwords);
    return ERRCODE_FAILURE;
  }

  for (i=0; i<ap->nwords; i++) {
    if (ap->wordidx[i] != bp->wordidx[i] ||
	ap->posidx[i] != bp->posidx[i])
      return ERRCODE_FAILURE;
  }

  return ERRCODE_SUCCESS;
}

HASHNUM_t 
hashTableGetKtupleHits(HASHPOS_t **posp, HASHNUM_t *posidx, 
		       const HashTable *htp, HASHWORD_t word)
{
  HASHKEY_t key;
  HASHNUM_t nhits = 0;

  if (htp->typ == HASHIDXTYP_PERFECT) {
    key = word&htp->wordmask;
    if (posidx) *posidx = key;
    if (key < htp->nkeys) {
      nhits = htp->idx[key+1] - htp->idx[key];
      if (posp) *posp = htp->pos + htp->idx[key];
    }
  } else {
    /* htp->typ != HASHIDXTYP_PERFECT */
    HASHKEY_t key_hi;
    uint32_t word_hi;
    uint32_t pivot, b, a;
    MAKE_HASHKEY(word, word_hi, key_hi, key);
    b = htp->idx[key+1];
    if (b < 1)
      return ERRCODE_SUCCESS;
    a = htp->idx[key];
    b--;
    while (a < b) {
      pivot = (a + b) >> 1;
      if (htp->wordidx[pivot] < word_hi) {
	/* upper half */
	a = pivot + 1;
      } else {
	/* lower half */
	b = pivot;
      }
    }
    if (a == b && htp->wordidx[b] == word_hi) {
      nhits = htp->posidx[b+1] - htp->posidx[b];
      if (posidx) 
	*posidx = b;
      if (posp)
	*posp = htp->pos + htp->posidx[b];
    }
  }

  return nhits;
}

HASHNUM_t 
hashTableFetchHitPositions(HASHPOS_t **posp, const HashTable *htp, HASHNUM_t posidx)
{
  HASHNUM_t nhits = 0;

  *posp = NULL;

  if (htp->typ == HASHIDXTYP_PERFECT) {
    if (posidx < htp->nkeys) {
      nhits = htp->idx[posidx+1] - htp->idx[posidx];
      *posp = htp->pos + htp->idx[posidx];
    }
  } else {
    if (posidx < htp->npos) {
      nhits = htp->posidx[posidx+1] - htp->posidx[posidx];
      *posp = htp->pos + htp->posidx[posidx];
    }
  }
  return nhits;
}

int hashTableWrite(const char *filnam, const HashTable *htp)
{
  int errcode;
  uint32_t header[HASHTABFIL_HEADSIZ];
  size_t totsiz;
  FILE *fp;
  /* calculate size of the table */
  header[0] = htp->wordlen;
  header[1] = htp->nskip;
  header[2] = htp->npos;
  header[3] = htp->maxpos;
  header[4] = htp->typ;
  header[5] = htp->nbits_key;
  header[6] = htp->nbits_lo;
  header[7] = htp->nwords;

  totsiz = htp->npos + htp->nkeys + 1;
  if (htp->typ != HASHIDXTYP_PERFECT)
    totsiz += (htp->nwords + 1)*2;
  if (totsiz > UINT32_MAX)
    return ERRCODE_OVERFLOW;

  fp = filioOpenForWriting(&errcode, (uint32_t) totsiz, FILIOTYP_HASHTAB,
			   HASHTABFIL_VERSION, HASHTABFIL_HEADSIZ,
			   header, filnam, HASHTABFIL_NAMEXT);
  if (errcode)
    return errcode;

  fwrite(htp->idx, sizeof(uint32_t), htp->nkeys+1, fp);
  fwrite(htp->pos, sizeof(uint32_t), htp->npos, fp);
  if (htp->typ != HASHIDXTYP_PERFECT) {
    fwrite(htp->wordidx, sizeof(uint32_t), htp->nwords + 1, fp);
    fwrite(htp->posidx, sizeof(uint32_t), htp->nwords + 1, fp);
  }
  if (ferror(fp)) {
    perror(HASHTABFIL_WRITERRMSG);
    errcode = ERRCODE_WRITEERR;
  }
  
  fclose(fp);
  return errcode;
}

HashTable *hashTableRead(int *errcode, const char *filnam)
{
  UCHAR_t is_endianid, typ, hashtyp, nbits_perf, nbits_key;
  uint32_t header[HASHTABFIL_HEADSIZ];
  uint32_t totsiz, version;
  uint32_t headsiz = HASHTABFIL_HEADSIZ;
  size_t nr;
  HashTable *htp = NULL;
  FILE *fp = filioOpenForReading(errcode, &is_endianid, &totsiz,
				 &typ, &version,
				 &headsiz, header, filnam, HASHTABFIL_NAMEXT);
  if (*errcode)
    return 0;

  if (typ != FILIOTYP_HASHTAB)
    *errcode = ERRCODE_FILTYP;
  else if (version != HASHTABFIL_VERSION && version != HASHTABFIL_OLDVERSION)
    *errcode = ERRCODE_FILVERSION;
  else if (headsiz > HASHTABFIL_HEADSIZ)
    *errcode = ERRCODE_FILEFORM;
  if (*errcode) {
    fclose(fp);
    return 0;
  }

  if (header[0] > MAXTUPLEN ||
      header[1] > UINT8_MAX ||
      (version == HASHTABFIL_VERSION &&
       (header[5] > MAX_NBITS_KEY ||
	header[6] >= header[5]))) {
    *errcode = ERRCODE_FILEFORM;
    fclose(fp);
    return 0;
  }
  
  if (version == HASHTABFIL_VERSION) {
    hashtyp = (UCHAR_t) header[4];
    nbits_key = (UCHAR_t) header[5];
    nbits_perf = (UCHAR_t) header[6];
  } else {
    hashtyp = HASHIDXTYP_PERFECT;
    nbits_key = (UCHAR_t) header[0];
    nbits_perf = 0;
  } 
  htp = hashTableCreate((UCHAR_t) header[0], (UCHAR_t) header[1],
			nbits_key, nbits_perf, hashtyp);
  if (htp == NULL) {
    *errcode = ERRCODE_NOMEM;
  } else {
    htp->npos = header[2];
    htp->maxpos = header[3];
    ECALLOCP(htp->npos, htp->pos);
    if (htp->pos == NULL) {
      *errcode = ERRCODE_NOMEM;
    } else if (htp->typ != HASHIDXTYP_PERFECT) {
      if (version != HASHTABFIL_VERSION) {
	*errcode = ERRCODE_FILEFORM;
      } else {
	htp->nwords = header[7];
	ECALLOCP((htp->nwords + 1)*2, htp->wordidx);
	if (htp->wordidx == NULL)
	  *errcode = ERRCODE_NOMEM;
	else 
	  htp->posidx = htp->wordidx + htp->nwords + 1;
      }
    }
  }
  if ((*errcode)) {
    fclose(fp);
    return 0;
  }

  nr = fread(htp->idx, sizeof(uint32_t), htp->nkeys+1, fp);
  if (nr == htp->nkeys+1) {
    nr = fread(htp->pos, sizeof(uint32_t), htp->npos, fp);
    if (nr == htp->npos) {
      if (htp->typ != HASHIDXTYP_PERFECT) {
	nr = fread(htp->wordidx, sizeof(uint32_t), 2*(htp->nwords) + 1, fp);
	if (nr != 2*(htp->nwords) + 1)
	  *errcode = ERRCODE_FILEFORM;
      }
    } else {
      *errcode =  ERRCODE_FILEFORM;
    }
  } else {
    *errcode =  ERRCODE_FILEFORM;
  }
  if (ferror(fp)) {
    *errcode = ERRCODE_READERR;
  }
  if (!(is_endianid) && !(*errcode)) {
  /* if (fread(htp->idx, sizeof(uint32_t), htp->nkeys+1, fp) != htp->nkeys+1 || */
  /*     fread(htp->pos, sizeof(uint32_t), htp->npos, fp) != htp->npos || */
  /*     (htp->typ != HASHIDXTYP_PERFECT &&  */
  /*      fread(htp->wordidx, sizeof(uint32_t), 2*(htp->nwords) + 1, fp) != 2*(htp->nwords) + 1)) { */
  /*   *errcode = ERRCODE_FILEFORM; */
  /* } else if ((ferror(fp))) { */
  /*   *errcode = ERRCODE_READERR; */
  /* } else if (!(is_endianid)) { */
    filioSwapEndian(htp->idx, htp->nkeys+1);
    filioSwapEndian(htp->pos, htp->npos);
    filioSwapEndian(htp->wordidx, 2*(htp->nwords+1));
  }
  fclose(fp);

  if (!(*errcode))
    htp->status = HASHSETUP_COMPLETE;

  return htp;  
}

  
