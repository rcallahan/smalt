/**< Simulate reads from genomic sequences
 */

/****************************************************************************
 ****************************************************************************
 *                                                                          *
 *  Copyright (C) 2010-2014 Genome Research Ltd.                            * 
 *                                                                          *
 *  Author: Hannes Ponstingl (hp3@sanger.ac.uk)                             *
 *                                                                          *
 *  This file is part of SMALT                                              *
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdint.h>
#include <math.h>

#include "elib.h"
#include "array.h"
#include "randef.h"
#include "rsample.h"
#include "sequence.h"

//#define SIMREAD_DEBUG

/* switch output format */
//#define simread_alistrtyp_explicit

/* Output formats:
 * (([dims])(\d+))+   if simread_alistrtyp_explicit is defined
 * or
 * ((\d+)([dims]))+   else
 * 
 * The number (\d+) represents the number of exact matches
 * d is a 1-base deletion in the read with respect to the reference
 * i is a 1-base insertion in the read with respect to the reference
 * m indicates exact matches up to the end of the alignment
 * s is a 1-base substitution.
 * 
 */

enum {
  LINBUFSIZ = 4096,
  LINBUFSIZ_MARGIN = 16,
  BLOCKSIZ_READ = 256,
  INDEL_BLKSZ = 512*1024,
  DEFAULT_PERCINDEL = 20, /**< Percentage of indels among total variations */
  //  RANDOM_SEED = 0,  /**< Seed for random number generation, 0 -> derived from calendar time. */ 
  MAXLEN_HISTO = 10,/**< Maximum indel length for histogram */
  MAXLEN_INDEL = 16, /**< maximum indel size */
  MAX_QVAL = 99,     /**< maximum base quality value */
  DEFAULT_QUAL = 'I',
  MAX_N_PER_READ = 5,
  ERRCODE_TOOMANY_N = 80,
  ERRCODE_READLEN = 81,
  ERRCODE_SOURCE_LEN = 82,
  ERRCODE_MATEPOS = 83,
  MAXNAMLEN = 256, /**< Maximum length for a read name */
};

enum VARIAT_TYPES {
  VARTYP_UNKNOWN      = 0,
  VARTYP_SUBSTITUTION = 1, /**< Single-base substitution */
  VARTYP_INSERTION    = 2,
  VARTYP_DELETION     = 3,
};

enum VARIAT_CHARS {
  VARCHAR_SUBSTITUTION = 's',
  VARCHAR_INSERTION = 'i',
  VARCHAR_DELETION = 'd',
};

enum ALISTR_TYPES {
  ALISTRTYP_EXPLICIT = 0,
  ALISTRTYP_CIGARLIKE = 1,
};

/*****************************************************************************
 ************************* Distribution parameters ***************************
 *****************************************************************************/

static const double  DEFAULT_INDELSIZ_GEOM_PROB = 0.7;
//static const double INSERTSIZ_STD_1 = 0.1; /* standard deviation of insert sizes */
static const double SMALL_QVAL = 1e-9;
/*****************************************************************************
 ******************************** Type Defs **********************************
 *****************************************************************************/

typedef unsigned char UCHAR;
typedef uint32_t SEQPOS;
typedef long long unsigned int LLUINT;
typedef int INSIZ;

typedef struct _VARIAT {
  UCHAR typ;   /**< one of VARIAT_TYPES */
  UCHAR len;   /**< length of the variation (indel size), is 1 for typ == VARTYP_SUBSTITUTION */
  uint64_t bas;  /**< Number of the generated base at which the variation occurs. Indels are simulated
		* after that base */
} VARIAT;


static char *copyStrField(char *tostr, const char *fromstr, int maxnum)
{
  int i;
  for (i=0; i<maxnum-1; i++) {
    if (fromstr[i] == ' ')
      break;
    tostr[i] = fromstr[i];
  }
  tostr[i] = '\0';
  return tostr;
}

static int cmpVARIAT(const void *ap, const void *bp)
{
  uint64_t a = ((const VARIAT *) ap)->bas;
  uint64_t b = ((const VARIAT *) bp)->bas;
  
  if (a > b)
    return 1;
  if (a < b)
    return -1;
  return 0;
}

static void fprintfVARIAT(FILE *fp, const VARIAT *varp, 
			  const uint32_t varnum_start, const uint32_t varnum_end)
{
  uint32_t v;
  const VARIAT *vp;
  for (v=varnum_start; v<varnum_end; v++) {
    vp = varp + v;
    fprintf(fp, "VARIAT[%2u]: typ = %hu, len = %hu, bas = %llu\n",
	    v, (unsigned short) vp->typ, (unsigned short) vp->len, (LLUINT) vp->bas);
  }
}

static int fprintfVariationStats(FILE *fp, const VARIAT *vp, const uint32_t varnum)
{
  uint32_t v;
  short h;
  UCHAR typ, l, max_ins, max_del;
  uint32_t ins_len_histo[MAXLEN_HISTO+1];
  uint32_t del_len_histo[MAXLEN_HISTO+1];
  uint32_t n_sub, n_ins, n_del;

  n_sub = n_ins = n_del = 0;
  max_ins = max_del = 0;
  for (h=0; h<= MAXLEN_HISTO; h++) {
    ins_len_histo[h] = del_len_histo[h] = 0;
  }

  for (v=0; v<varnum; v++) {
    
    l = vp[v].len;
    typ = vp[v].typ;

    if (typ == VARTYP_SUBSTITUTION) {
      n_sub++;
    } else if (typ == VARTYP_INSERTION) {
      n_ins++;
      if (l <= MAXLEN_HISTO)
	ins_len_histo[l]++;
      if (l > max_ins)
	max_ins = l;
    } else if (typ == VARTYP_DELETION) {
      n_del++;
      if (l <= MAXLEN_HISTO)
	del_len_histo[l]++;
      if (l > max_del)
	max_del = l;
    }
  }
  fprintf(fp, "============ Variation Statistics ===========\n");
  fprintf(fp, "%u variations generated:\n", varnum);
  fprintf(fp, "%u single-base substitutions\n", n_sub);
  fprintf(fp, "%u insertions\n", n_ins);
  fprintf(fp, "%u deletions\n", n_del);
  fprintf(fp, "histogram of insertions | deletions :\n");
  for (h=0; h<= MAXLEN_HISTO; h++)
    fprintf(fp, "[%2i] %4u | %4u\n", h, 
	    ins_len_histo[h], del_len_histo[h]);
  fprintf(fp, "maximum insertion length: %hu\n", (unsigned short) max_ins);
  fprintf(fp, "maximum deletion length: %hu\n", (unsigned short) max_del);
  fprintf(fp, "======== End of Variation Statistics ========\n");

  return ERRCODE_SUCCESS;
}

static int drawRandomVariations(VARIAT *vp, const uint32_t varnum,  
				const uint64_t basnum,
				const int nindel_percent,
				const double indelsz_pgeom)
     /**< Generate random positions for single-base variations and for 
      * indels.
      * \param vp Array of varnum variations
      * \param varnum number of variations (size of array vp)
      * \param array of indel sizes.
      * \param basnum Total number of bases to be simulated (number of reads x read length)
      * \param nindel_percent Number of indels as percentage of total number of variations.
      * \param indelsz_geomp Parameter p for geometric distribution p(1-p)^(k-1) of insert sizes.
      */
{
  int errcode;
  double r;
  int *indelszr;
  uint32_t indelnum, v,n, ctr=0;
  uint32_t n_ins, n_del;

  /* initalize array */
  for (v=0; v<varnum; v++) {
    vp[v].typ = VARTYP_SUBSTITUTION;
    vp[v].len = 1;
    r = RANDRAW_UNIFORM_1();
    vp[v].bas = r*basnum;
  }

  /* determine the number of indels amongst the varnum variations */
  indelnum = (varnum*nindel_percent)/100;
  if (indelnum > varnum) 
    return ERRCODE_ASSERT;

  if (indelnum < 1)
    return ERRCODE_SUCCESS;

  ECALLOCP(indelnum + 1, indelszr);
  if (!indelszr)
    return ERRCODE_NOMEM;

  if ((errcode = rsampleGeometric(indelszr, indelnum + 1, indelsz_pgeom)))
    return errcode;

  n_ins = indelnum/2;       /* number of insertions */
  n_del = indelnum - n_ins; /* number of deletions */

  /* generate insertions */
  for (n=0; n<n_ins; n++) {
    r = RANDRAW_UNIFORM_1();
    v = (uint32_t) (r*varnum);
    vp[v].typ = VARTYP_INSERTION;
    if (indelszr[ctr] > MAXLEN_INDEL || indelszr[ctr] > UCHAR_MAX)
      vp[v].len = 1;
    else
      vp[v].len = (UCHAR) (indelszr[ctr] + 1);
    ctr++;
  }

  /* generate deletions */
  for (n=0; n<n_del; n++) {
    r = RANDRAW_UNIFORM_1();
    v = (uint32_t) (r*varnum);
    vp[v].typ = VARTYP_DELETION;
    if (indelszr[ctr] > MAXLEN_INDEL || indelszr[ctr] > UCHAR_MAX)
      vp[v].len = 1;
    else
      vp[v].len = (UCHAR) (indelszr[ctr] + 1);

    ctr++;
  }
  
  free(indelszr);
  return ERRCODE_SUCCESS;
}

static SETSIZ_t getNumberOfGenomicBases(const SeqSet * const ssp)
{
  const SETSIZ_t *sop;
  SEQNUM_t no = seqSetGetOffsets(ssp, &sop);
  return sop[no];
}

#ifdef SIMREAD_SUPERFLUOUS
static int parseIndelSizes(UCHAR **isizr, const char *filnam)
{
  char linbuf[LINBUFSIZ];
  int il;
  UCHAR *hp;
  FILE *fp = EFOPEN(filnam, "r");

  if (!fp)
    return ERRCODE_NOFILE;
  
  while (fgets(linbuf, LINBUFSIZ, fp)) {
    il = atoi(linbuf);
    if (il < 1 || il > MAXLEN_INDEL)
      il = 1;
    ARRCURR(*isizr) = (UCHAR) il;
    ARRNEXTP(hp, *isizr);
    if (hp == NULL)
      return ERRCODE_NOMEM;
  }

  return EFCLOSE(fp);
}
#endif

static int simulateSingleRead(SeqFastq *readp, SeqFastq *sbufp, 
			      char *targetp, char *target_qualp,
			      const char *readnamprefix,
			      int readnum, UCHAR is_reverse, int mateno,
			      const VARIAT *varp, const uint32_t varnum,
			      const int readlen, const uint64_t basctr, 
			      uint32_t pos, 
			      const SeqSet *ssp, const SeqCodec *codecp)
     /** mateno 0 for single read, 1 for first, 2 for 2nd mate of pair */
{
  int errcode, d, n, t, fetch_len;
  uint32_t s;
  SEQNUM_t sidx;
  char c, *sp;
  char read_name[LINBUFSIZ];
  char alistr[LINBUFSIZ];
  UCHAR typ, cod, l, last_typ;
  short alphlen, encodlen, N_ctr;
  uint64_t exact_len, ctr;
  size_t a;
  uint32_t v, slen, reflen, soffs, so;
  double r;
  const char *alphabetp = seqCodecGetAlphabet(codecp, &alphlen);
  const UCHAR *encoderp = seqCodecGetEncoder(codecp, &encodlen);
  const char *namp;
  char refnamp[MAXNAMLEN];

  read_name[0] = '\0';
  seqFastqBlank(readp);
  if ((errcode = seqSetGetIndexAndOffset(&sidx, &soffs, pos, ssp)))
    return errcode;
  reflen = seqSetGetSeqDatByIndex(NULL, &namp, sidx, ssp);
  copyStrField(refnamp, namp, MAXNAMLEN);

  so = pos - soffs;
  /* find out, how many extra bases have to be fetched */
  fetch_len = readlen;
  exact_len = 0;
  ctr = basctr;
  for (v=0; v<varnum; v++) {
    exact_len += varp[v].bas - ctr + 1;
    ctr = varp[v].bas;
    if (varp[v].typ == VARTYP_DELETION) {
      fetch_len += (int) varp[v].len;
    }
  }
  if (exact_len > INT_MAX)
    return ERRCODE_OVERFLOW;
  if (fetch_len < (int) exact_len) 
    fetch_len = (int) exact_len;

  if ((so + fetch_len) > reflen)
    return ERRCODE_SOURCE_LEN;

  if ((errcode = seqSetFetchSegmentBySequence(sbufp, sidx, 
					      so, (uint32_t) fetch_len,
					      ssp, 
					      codecp)))
    return errcode;
  sp = seqFastqGetSequence(sbufp, &slen, NULL);
  if (((uint32_t) fetch_len) != slen)
    return ERRCODE_READLEN;

  /* make sure there are not too many Ns */
  for (N_ctr=0, s=0; s<slen && N_ctr < MAX_N_PER_READ; s++) {
    if (encoderp[(UCHAR) sp[s]]&SEQCOD_STDNT_TESTBIT)
      N_ctr++;
  }
  if (N_ctr >=  MAX_N_PER_READ)
    return ERRCODE_TOOMANY_N;

  /* generate the simulated sequence */
  s = t = 0;
  n = 0;
  c = sp[0];
  ctr = basctr;
  typ = VARTYP_UNKNOWN;
#ifdef simread_alistrtyp_explicit
  sprintf(alistr, "[");
#else
  alistr[0] = '\0';
#endif
  last_typ = VARTYP_UNKNOWN;
  a = strlen(alistr);
  for (v=0; v<varnum && t<readlen && s<slen; v++) {

    if (a + LINBUFSIZ_MARGIN > LINBUFSIZ)
      return ERRCODE_OVERFLOW;

    if (varp[v].bas < ctr)
      continue;

    exact_len = varp[v].bas - ctr + 1;
    if (exact_len < 1 && 
	((last_typ == VARTYP_DELETION && varp[v].typ == VARTYP_INSERTION) || 
	 (last_typ == VARTYP_INSERTION && varp[v].typ == VARTYP_DELETION)))
      continue;
    n = (int) exact_len; /* number of nucleotides including position of the variation */
    ctr = varp[v].bas+1;
    typ = varp[v].typ;
 
    
    for (d=0; t<readlen && s<slen && exact_len > 0; t++, s++, exact_len--, d++) {
      c = sp[s];
      targetp[t] = c;
    }
    /* d > 0 is the number of exact matches up to variation or end of alignment.
     * If the end of the alignment is reached, d includes the last
     * match. 
     * 
     * t is the current position */

    if (t>=readlen && (exact_len > 0 || typ != VARTYP_SUBSTITUTION)) {
      if (d > 0) {
	sprintf(alistr+a, 
#ifdef simread_alistrtyp_explicit
		"m%hi]", 
#else
		"%him",
#endif
		d);
	a = strlen(alistr);
      }
      break;
    }
    if (typ == VARTYP_SUBSTITUTION) {
      sprintf(alistr+a,
#ifdef simread_alistrtyp_explicit
	      "s%i:", 
#else
	      "%is",
#endif
	      n);
      cod = (UCHAR)(encoderp[(UCHAR) c]&SEQCOD_STDNT_MASK);
      r = RANDRAW_UNIFORM_1();
      cod = (UCHAR)(cod + (SEQCOD_STDNT_MASK)*r + 1);
      cod %= (SEQCOD_STDNT_MASK+1);
      targetp[t-1] = alphabetp[cod];
#ifdef SIMREAD_DEBUG
      fprintf(stderr, "DEBUG::simulateSingleRead: VAR[%u] substitute %llu + %i, pos = %llu\n", 
	     v, (unsigned long long) basctr, t, 
	     (unsigned long long) varp[v].bas);
#endif
    } else if (typ == VARTYP_INSERTION) {
      sprintf(alistr+a, 
#ifdef simread_alistrtyp_explicit
	      "i%i:",
#else
	      "%ii",
#endif 
	      n);
      n += varp[v].len;
      for (l=0; l<varp[v].len && t<readlen; l++, t++) {
	r = RANDRAW_UNIFORM_1();
	cod = (UCHAR) ((SEQCOD_STDNT_MASK+1)*r);
	targetp[t] = alphabetp[cod];
	if (l>0) {
	  a = strlen(alistr);
	  sprintf(alistr+a, 
#ifdef simread_alistrtyp_explicit
		  "i0:"
#else
		  "0i"
#endif
		  );
	}
      }
    } else if (typ == VARTYP_DELETION) {
      sprintf(alistr+a,
#ifdef simread_alistrtyp_explicit 
	      "d%i:", 
#else
	      "%id",
#endif
	      n);
      n += varp[v].len;
      s++;
      for (l=1; l<varp[v].len && s<slen; l++, s++) {
	a = strlen(alistr);
	sprintf(alistr+a, 
#ifdef simread_alistrtyp_explicit 
		  "d0:"
#else
		  "0d"
#endif
		  );
      }
    }
    a = strlen(alistr);
  }

  n = readlen - t;
  if (n>0) {
    sprintf(alistr+a,
#ifdef simread_alistrtyp_explicit 
	    "m%hi]", 
#else
	    "%him",
#endif
	    n);
    for (; t<readlen && s<slen && n > 0; t++, s++, n--) {
      targetp[t] = sp[s];
    }
  } else {
#ifdef simread_alistrtyp_explicit 
    alistr[a++] = ']';
#else
    alistr[a] = '\0';
#endif
  }
  if (t < readlen)
    return ERRCODE_READLEN;
  targetp[t] = '\0';

#ifdef simread_alistrtyp_explicit 
  sprintf(read_name, "%s_%9.9i_%s_%9.9u_%c_%s", readnamprefix, readnum, 
	  refnamp, so+1, (is_reverse)? 'R':'F', alistr);
#else 
  sprintf(read_name, "%s_%9.9i_%s_%9.9u_%lli_%c_%s", readnamprefix, readnum, 
	  refnamp, so+1, (long long signed) sidx, (is_reverse)? 'R':'F', alistr);
#endif
  if (mateno > 0)
    sprintf(read_name + strlen(read_name), "/%1i", mateno);
  errcode = seqFastqSetAscii(readp, read_name, targetp, "", target_qualp);
  if (!errcode && (is_reverse))
    errcode = seqFastqReverse(readp, codecp);
  return errcode;
}

static int simulatePairedRead(SeqFastq *readp, SeqFastq *matep, 
			      SeqFastq *sbufp, 
			      char *target_seqp, char *target_qualp,
			      const char *readnamprefix,
			      int readnum, UCHAR is_reverse,
			      const VARIAT *varp, uint32_t *voffs, const uint32_t varnum,
			      const int readlen, int insertsiz,
			      SETSIZ_t *basctr, uint32_t pos, 
			      const SeqSet *ssp, const SeqCodec *codecp)
/*
 * \param insertsiz length of insert (measured from 5'-ends). 
 * 0 signals single reads.
 */
{
  int errcode;
  int d_insert = insertsiz - readlen;
  UCHAR is_paired = (UCHAR) (insertsiz != 0 && (matep));
  uint32_t vs, ve, mpos;
  uint64_t bctr, next_bctr;

/*   if ((is_paired) && d_insert < 0) */
/*     return ERRCODE_ASSERT; */

  bctr = *basctr;
  for (vs = *voffs; vs < varnum && varp[vs].bas < bctr; vs++);
  next_bctr = bctr+readlen;
  for (ve=vs; ve<varnum && varp[ve].bas < next_bctr; ve++);

  errcode = simulateSingleRead(readp, sbufp,
			       target_seqp, target_qualp,
			       readnamprefix,
			       readnum, is_reverse, (is_paired)? 1:0,
			       varp+vs, ve-vs, readlen, 
			       bctr, pos, ssp, codecp);
  if (is_paired && insertsiz > 0 && !errcode) {
    if ((is_reverse)) {
      if (d_insert >= 0 && pos < ((uint32_t) d_insert)) 
	return ERRCODE_MATEPOS;
      mpos = pos - d_insert;
    } else {
      if (d_insert < 0 && pos < ((uint32_t) (-d_insert)))
	return ERRCODE_MATEPOS;
      mpos = pos + d_insert;
    }
    vs = ve;
    bctr = next_bctr;
    next_bctr += readlen;
    for (;ve<varnum && varp[ve].bas < next_bctr; ve++);
    errcode = simulateSingleRead(matep, sbufp,
				 target_seqp, target_qualp,
				 readnamprefix,
				 readnum, (UCHAR) !is_reverse, 2,
				 varp+vs, ve-vs, readlen, 
				 bctr, mpos, ssp, codecp);
  }
  if (!errcode) {
    *voffs = ve;
    *basctr = next_bctr;
  }

  return errcode;
}

static int generateReads(SeqFastq *readp, SeqFastq *matep, 
			 SeqFastq *sbufp, SeqIO *sfpA, SeqIO *sfpB,
			 char *target_seqp, char *target_qualp,
			 const VARIAT *varp, const uint32_t varnum, 
			 const char *readnamprefix,
			 int readlen, int readnum, 
			 int insertsiz, const double *ins1p,
			 const SeqSet *ssp, const SeqCodec *codecp)
/* 
 * \param readnum Number of reads or read pairs.
 * \param insertsiz mean of the insert size distribution (can be 0 for single reads
 * or -1 for randomly paired reads)
 * \param ins1p random vector of insert sizes (can be NULL) with mean insertsiz
 */
{
  int errcode, i, i_end, rn, pairctr, n_error=0, n_readlen_err=0, n_skipped=0;
  UCHAR is_reverse, is_paired = (UCHAR) ((matep) && insertsiz > 0);
  int isiz;
  uint32_t vs;
  SETSIZ_t pos, reflen = getNumberOfGenomicBases(ssp);
  SETSIZ_t basctr;
  double r, rnfac;

  if (is_paired && insertsiz < readlen)
    return ERRCODE_ASSERT;

  rnfac = (insertsiz < 0)? 3.0: 1.5;
  i_end = (readnum*rnfac > INT_MAX)? INT_MAX:(int) (readnum*rnfac);

  vs = 0;
  basctr = 0;
  pairctr = 0;
  for (i=rn=0; rn<readnum && i<i_end; i++) {
    /* random position */
    r = RANDRAW_UNIFORM_1();
    pos = r*reflen;

    /* reverse or forward? */
    r = RANDRAW_UNIFORM_1();
    is_reverse = (UCHAR)(2*r);

    /* insert size */
    if (is_paired && (ins1p)) {
      isiz = (int) ins1p[rn];
    } else if (insertsiz >= 0) {
      isiz = insertsiz;
    } else {
      isiz = 0;
    }
    errcode = simulatePairedRead(readp, matep, sbufp,
				 target_seqp, target_qualp,
				 readnamprefix,
				 rn, is_reverse,
				 varp, &vs, varnum, 
				 readlen, isiz,
				 &basctr, pos, ssp, codecp);
    if ((errcode)) {
      if (errcode == ERRCODE_TOOMANY_N || errcode == ERRCODE_MATEPOS) {
	n_skipped++;
	if (!(n_skipped % 100))
	  fprintf(stderr, "skipped %i of %i (n_var = %u, error rate = %g) ...\n", 
		  n_skipped, i, vs, ((double) vs)/basctr);
	continue;
      }
      if (errcode == ERRCODE_READLEN || errcode == ERRCODE_SOURCE_LEN) {
	n_readlen_err++;
	fprintf(stderr, "skipped %i of %i, read length\n", n_readlen_err, i);
	continue;
      }
      n_error++;
      fprintf(stderr, "skipped %i of %i, ERROR code: %i ...\n", n_error, i, errcode);
      continue;
    }

    if (insertsiz >= 0) {
      if ((errcode = seqFastqWrite(sfpA, readp, 0)))
	return errcode;
      
      if ((matep) && insertsiz > 0) {
	if ((errcode = seqFastqWrite(sfpB, matep, 0)))
	  return errcode;
      }
    } else {
      if ((errcode = seqFastqWrite(((pairctr))? sfpB: sfpA, readp, 0)))
	return errcode;
    }

    if (insertsiz < 0 && pairctr == 0) {
      pairctr++;
    } else {
      pairctr = 0;
      rn++;
    }
  }
  return ERRCODE_SUCCESS;
}

int main(int argc, char *argv[])
{
  int errcode;
  int i;
  int readnum, readlen, insertsiz, randseed;
  double percerr, *ins1p, stdisiz, insertstd1;
  uint64_t basnum;
  uint32_t varnum;
  SETSIZ_t reflen;
  char target_qual;
  char *binfilnam, *oufilnam, *oufilnamA, *oufilnamB;
  char *target_seqp, *target_qualp, *readnamprefix;
  double indelsz_pgeom = DEFAULT_INDELSIZ_GEOM_PROB;
  int nindel_percent = DEFAULT_PERCINDEL;
  VARIAT *varp;
  SeqCodec *codecp;
  SeqIO *sfpA, *sfpB;
  SeqFastq *readp, *sbufp, *matep = 0;
  SeqSet *ssp;
  ErrMsg *errmsgp=0;
  
  ERRMSG_CREATE(errmsgp);

  if (argc < 11) {
    fprintf(stderr, "usage: %s <compressed sequence file (w/o ext)> ", argv[0]);
    fprintf(stderr, "<read length> <num reads> <error rate[%%]> <with indels [y|n]>");
    fprintf(stderr, "<insert size (0 for single reads, -1 random pairings)> ");
    fprintf(stderr, "<insert size std_1> <seed (0 for calendar derrived)> ");
    fprintf(stderr, "<read name prefix> ");
    fprintf(stderr, "<output file>\n");
    exit(EXIT_FAILURE);
  }

  binfilnam = argv[1];
  readlen = atoi(argv[2]);
  readnum = atoi(argv[3]);
  percerr = atof(argv[4]);
  if ((argv[5][0] != 'y') && (argv[5][0] != 'Y'))
    nindel_percent = 0;

  insertsiz = atoi(argv[6]);
  /* insertsiz == 0: produce single reads, 
   * insertsiz == -1: produce random pairings (no fixed inserts size, different chromosomes)
   */
  insertstd1 = atof(argv[7]);
  randseed = atoi(argv[8]);
  readnamprefix = argv[9];
  oufilnam = argv[10];

  if (readlen < 1) {
    fprintf(stderr, "Invalid read length: %i\n", readlen);
    exit(EXIT_FAILURE);
  }

  if (readnum < 1) {
    fprintf(stderr, "Invalid number of reads: %i\n", readnum);
    exit(EXIT_FAILURE);
  }

  if (percerr < 0.0 || percerr > 100.0) {
    fprintf(stderr, "Invalid error rate: %2.2f %%\n", percerr);
    exit(EXIT_FAILURE);
  }

  ECALLOCP(readlen+1, target_seqp);
  ECALLOCP(readlen+1, target_qualp);

  if (!target_seqp || !target_qualp)
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  if (percerr > SMALL_QVAL*100) {
    double qv = -10*log10(percerr/100);
    if (qv > (double) MAX_QVAL) 
      qv = MAX_QVAL;
    else if (qv < 0)
      qv = 0.0;
    target_qual = (char) (SEQCOD_QVAL_OFFS + qv);
  } else {
    target_qual = DEFAULT_QUAL;
  }
  for (i=0; i<readlen; i++)
    target_qualp[i] = target_qual;
  target_qualp[i] = '\0';

  if (!(codecp = seqCodecCreate()))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  if (!(sbufp = seqFastqCreate(BLOCKSIZ_READ, SEQTYP_FASTQ)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  if (!(readp = seqFastqCreate(BLOCKSIZ_READ, SEQTYP_FASTQ)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);
  
  if (insertsiz > 0) {
    if (!(matep = seqFastqCreate(BLOCKSIZ_READ, SEQTYP_FASTQ)))
      ERRMSGNO(errmsgp, ERRCODE_NOMEM);
  } else {
    matep = NULL;
  }
  
  if (insertsiz != 0) {
    ESTRCAT(oufilnamA, oufilnam, "_1.fq");
    ESTRCAT(oufilnamB, oufilnam, "_2.fq");
    if (!(oufilnamB))
      ERRMSGNO(errmsgp, ERRCODE_NOMEM);
  } else {
    ESTRCAT(oufilnamA, oufilnam, ".fq");
    oufilnamB = 0;
  }
  if (!(oufilnamA)) 
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  sfpA = seqIOopen(&errcode, oufilnamA, SEQIO_WRITE_FASTQ, 0);
  if (errcode) 
    ERRMSGNO(errmsgp, errcode);

  if (insertsiz != 0) {
    sfpB = seqIOopen(&errcode, oufilnamB, SEQIO_WRITE_FASTQ, 0);
    if (errcode) 
      ERRMSGNO(errmsgp, errcode);
  } else {
    sfpB = 0;
  }

  /* seed the random number generator */
  //RANSEED(RANDOM_SEED);
  RANSEED(randseed);

  printf("Simulate %i read%s a %ibp...\n", 
	 readnum, (insertsiz > 0)? "pairs":"s", readlen);

  basnum = ((uint64_t) readlen)*readnum;
  if (insertsiz != 0)
    basnum *= 2;
  printf("Simulate %llu bases ...\n", (LLUINT) basnum);

  varnum = (uint32_t) (basnum*percerr/100);
  printf("Simulate %u errors ...\n", varnum);

  ECALLOCP(varnum, varp);
  if (!varp)
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  if ((errcode = drawRandomVariations(varp, varnum, basnum, nindel_percent, indelsz_pgeom)))
    ERRMSGNO(errmsgp, errcode);

  fprintfVariationStats(stdout, varp, varnum);

  printf("Sorting variations by base number ...\n"); 
  qsort(varp, varnum, sizeof(VARIAT), cmpVARIAT);

  printf("Variations after sort ...\n");
  fprintfVARIAT(stdout, varp, 0, (varnum > 10)? 10: varnum);
  if (varnum>10) {
    printf("...\n");
    fprintfVARIAT(stdout, varp, varnum - 10, varnum);
  }

  if (insertsiz > 0) {
    //stdisiz = INSERTSIZ_STD_1*insertsiz;
    stdisiz = insertstd1*insertsiz;
    ECALLOCP(readnum, ins1p);
    if (!ins1p)
      ERRMSGNO(errmsgp, ERRCODE_NOMEM);
    
    if ((errcode = rsampleNormal(ins1p, readnum, 
				 (double) insertsiz, 
				 stdisiz*stdisiz)))
      ERRMSGNO(errmsgp, errcode);
  } else {
    stdisiz = 0.0;
    ins1p = NULL;
  }


  printf("Reading reference sequences ...\n");
  ssp = seqSetReadBinFil(&errcode, binfilnam);
  if ((errcode))
    ERRMSGNO(errmsgp, errcode);

  reflen = getNumberOfGenomicBases(ssp);
  printf("total length of reference sequences: %llu bp\n", (long long unsigned int) reflen);

  if ((errcode = generateReads(readp, matep, sbufp, sfpA, sfpB,
			       target_seqp, target_qualp,
			       varp, varnum, 
			       readnamprefix,
			       readlen, readnum, 
			       insertsiz, ins1p,
			       ssp, codecp)))
    ERRMSGNO(errmsgp, errcode);

  seqSetDelete(ssp);
  free(ins1p);
  free(varp);
  seqIOclose(sfpB);
  seqIOclose(sfpA);
  seqFastqDelete(sbufp);
  seqFastqDelete(matep);
  seqFastqDelete(readp);
  seqCodecDelete(codecp);
  free(oufilnamB);
  free(oufilnamA);
  free(target_seqp);
  free(target_qualp);
  ERRMSG_END(errmsgp);

  return EXIT_SUCCESS;
}

