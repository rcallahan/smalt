/** Handling of sequences (in FASTQ/FASTA format) */

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

#ifndef SEQUENCE_H
#define SEQUENCE_H

  //#define sequence_debug

#include <stdio.h>
#include <stdint.h>
#include "elib.h"

  /****************************************************************************
   ******************************* Constants **********************************
   ****************************************************************************/
  enum SEQCODEC_CODETYPES {
    SEQCODTYP_3BITMANGLED, /**< default type */
    SEQCODTYP_USERDEFINED, /**< user-defined from alphabet provided */
  };

  enum SEQCODEC_BASETYPES { /**< nucleotide classes */
    SEQBASTYP_UNKNOWN,
    SEQBASTYP_PURINE,
    SEQBASTYP_PYRIMIDINE,
    SEQBASTYP_NONSTD,      /**< non-standard nucleotide */
  };

  enum SEQ_CODES {/**< Type of sequence encoding */
    SEQCOD_ASCII = 0,      /**< Original ASCII characters */
    SEQCOD_MANGLED = 1,    /**< 8-bit code with bits 0,1 representing standard
			    *   nucleotide, bits 0-2 represent a 3-bit code of
			    *   the recognised alphabet including a
			    *   termination signal. Bits 3-7 represent the
			    *   original 7-bit ASCII letter in upper case as
			    *   the offset from the letter 'A'*/
    SEQCOD_COMPRESSED = 2, /**< Each base encoded in 3 bits, 10 bases per
			    * 32-bit integer.  The highest two bits of the
			    * int are unused. The termination signal are 3 set bits.
			    * the remaining bits are zeroed */
  };

  enum SEQIO_MODES { /**< Implemented input/output formats */
    SEQIO_READ        = 0,  /**< Read sequences sequentially, one at a time */
    SEQIO_WRITE_FASTA = 2,  /**< write in FASTA format */
    SEQIO_WRITE_FASTQ = 3,  /**< write in FASTQ format */
#ifdef HAVE_ZLIB
    SEQIO_WRITE_GZIP_FASTA = 4, /**< write gzipped FASTA format */
    SEQIO_WRITE_GZIP_FASTQ = 5, /**< write gzipped FASTQ format */
#endif
  };

  enum SEQ_TYPES {  /**< Implemented sequence types */
    SEQTYP_UNKNOWN = 0, /**< Unknown sequence type */
    SEQTYP_FASTA = 1,   /**< sequence in FASTA format */
    SEQTYP_FASTQ = 2,   /**< sequence in FASTQ format */
  };

  enum SEQSET_FLAGS {
    SEQSET_TERMCHAR = 0x01,   /**< offsets account for an addtional terminating 
			       * (sentinel) character in each sequence. */
    SEQSET_COMPRESSED = 0x02, /**< Sequence Set is compressed */
    SEQSET_BASQUAL = 0x04,    /**< Set contains base qualities */
  };

 
  enum SEQUENCE_CONST {
    SEQCOD_STDNT_TESTBIT = 0x04, /**< Test bit for standard nulcotides. Set bit indicates non-standard NT */
    SEQCOD_STDNT_MASK = 0x03,    /**< For masking out standard nucleotide code */
    SEQCOD_ALPHA_MASK = 0x07,    /**< For masking out alphabet code */
    SEQCOD_TERM = 0x07,          /**< 3-bit code for sequence termination */
    SEQCOD_QVAL_OFFS = 0x21,     /**< Quality values read as ASCII chars can be converted to
				  * PHRED type quality values by subtracting this offset */
  };

  /****************************************************************************
   ***************************** Simple types *********************************
   ****************************************************************************/
  enum SEQUENCE_LIMITS {
    SEQNUM_MAX = INT64_MAX-1,
    SEQLEN_MAX = UINT32_MAX-1,
    SETSIZ_MAX = UINT64_MAX-1,
  };

  typedef int64_t SEQNUM_t;     /**< Holds the number of sequences */
  typedef uint64_t SEQQN_t;     /**< Number of query sequences */
  typedef uint32_t SEQLEN_t;    /**< Length of an individual nucleotide sequence */
  typedef uint64_t SETSIZ_t;    /**< Total length of a set of (concatenated) sequences */

  /****************************************************************************
   ***************************** Opaque types *********************************
   ****************************************************************************/

  /** Sequence encoding/decoding */
  typedef struct _SeqCodec SeqCodec;
  /** Source of the sequence data */
  typedef struct _SeqIO SeqIO;
  /** Holds a sequence and quality string with headers */
  typedef struct _SeqFastq SeqFastq; 
  /** Holds a set of sequences in memory */
  typedef struct _SeqSet SeqSet;  
  
  /****************************************************************************
   *********************** Methods of type SeqCodec ***************************
   ****************************************************************************/
  
  /** \name SeqCodec Methods of type SeqCodec */
  /*@{*/
  SeqCodec *seqCodecCreate(void);
  /**< Constructor of a coder/decoder for 3-bit codes. U is converted to T.
   * \ingroup SeqCodec
   */
  SeqCodec *seqCodecCreateFromAlphabet(const char *alphabet, unsigned char code_unknown);
  /**< Create type SeqCodec from a specified alphabet.
   * \param alphabet Alphabet of at most 7 characters is transformed to upper case
   * \param code_unknown the index for the character that stands for 'unknown' (must be < 7)
   * 
   * \note the first 4 characters must be standard nucleotide codes. Converts to upper case
   * and U (uracyl) to T (threonine).
   * \ingroup SeqCodec
   */
  void seqCodecDelete(SeqCodec *codp);
  /**< Destructor 
   */
  const char *seqCodecGetAlphabet(const SeqCodec *codp, short *length);
  /**< Returns the alphabet used 
   */
  const unsigned char *seqCodecGetEncoder(const SeqCodec *codp, short *size);
  /**< Return pointer to encoding array (mangled code) and its size 
   */
  const char *seqCodecGetDecoder(const SeqCodec *codp, short *size);
  /**< Return pointer to decoding array (from mangled code) and its size 
   */
  void seqCodecEncodeString(char *cp, const SeqCodec *codp);
  /**< Encode a string in-place. String must be terminated by 0 
   */
  void seqCodecDecodeString(char *cp, const SeqCodec *codp);
  /**< Decode an  encoded string in-place. The string must be teminated by 0
   */
  char seqCodecFindBaseClass(char c, const SeqCodec *codecp);
  /**< Return the nucleotide class (one of SEQCODEC_BASETYPES)
   * for a base character.
   * \param c base (i.e. 'A', 'c', etc.)
   * \param codecp En-/Decoder, must be of type SEQCODTYP_3BITMANGLED
   * \note Codec must be of type SEQCODTYP_3BITMANGLED!
   */
  char seqCodecType(const SeqCodec *codecp);
  /**< return the encoding type (one of SEQCODEC_CODETYPES)
   */
  /*@}*/
  /****************************************************************************
   ************************* Methods of type SeqIO ****************************
   ****************************************************************************/

  SeqIO *seqIOopen(int *errcode, const char *filnam, char mode, 
		   unsigned long buffsize);
  /**< Open a file for sequence I/O in FASTA/FASTQ format
   * \param errcode != 0 (one of ERRMSG_CODES) 
   *        if error occurs in which case NULL is returned.
   * \param filnam file name. "-" opens standard I/O.
   * \param mode one of SEQIO_MODES or combination of SEQIOMOD_* flags
   * \param buffsize buffer size, if 0 default buffer size is used.
   */
  
  int seqIOReset(SeqIO *p);
  int seqIOclose(SeqIO *p);
  int seqIOstatus(SeqIO *p);

  int seqIOCheckReads(ErrMsg *errmsgp,
		      SeqFastq *sqbufAp, SeqIO *sfAp,
		      SeqFastq *sqbufBp, SeqIO *sfBp,
		      SEQNUM_t *seqnum, SEQLEN_t *maxseqlen, 
		      SEQLEN_t *maxnamlen);
  /**< Check whether the sequences in SeqIO conform to FASTA/FASTQ format for reading.
   *   Return error code or ERRCODE_SUCCESS. On return the file is scrolled back to the
   *   start.
   * \param errmsgp Returns error messages.
   * \param sqbufAp Sequence buffer, gets overwritten and reallocated, so that
   *              no memory allocation would be neccessary when subsequently using 
   *              this structure for reading a sequence.
   * \param sfAp Structure of type SeqIO. Must have been created for reading, otherwise 
   *          return ERRCODE_FAILURE
   * \param sqbufBp Sequence buffer for paired reads (NULL for single reads).
   * \param sfAp Structure of type SeqIO for paired reads (NULL for single reads).
   * \param[out] seqnum Returns the number if sequences (can be NULL).
   * \param[out] maxseqlen Returns the maximum sequence length (can be NULL).
   * \param[out] maxnamlen Returns the maximum name length (can be NULL).
   */


  /****************************************************************************
   *********************** Methods of type SeqFastq ***************************
   ****************************************************************************/

  SeqFastq *seqFastqCreate(int blocksize, char type);
  /**< Constructor.
   * \param blocksize Memory is dynamically allocated in blocks of blocksize bytes.
   * \param type one of SEQTYPES
   */

  void seqFastqDelete(SeqFastq *sqp);
  /**< Destructor */

  void seqFastqBlank(SeqFastq *p);
  /**< Reset sequence to the empty sequence but don't free allocated memory. */

  int seqFastqSetType(SeqFastq *sqp, char type);
  /**< Convert between sequence types. When converting from SEQTYP_FASTA, quality
   * values are zeroed.
   * \param sqp Sequence structure.
   * \param type one of SEQ_TYPES. 
   */
  
  void seqFastqFreeUnusedMem(SeqFastq *sqp);
  /**< Reallocate memory to the exact sequence length freeing unused memory. */

  int seqFastqCheck(const SeqFastq *sqp);
  /**< Check consistency of the structure */

  int seqFastqSetAscii(SeqFastq *sqp, 
		       const char *name, const char *seqp, 
		       const char *name_qual, const char *qualp);
  /**< Load name, base sequence and quality factors from ASCII strings.
   * \param sqp Structure to be loaded. Contents are overwritten, and allocated
   *            memory extended if necessary to accomodate the strings.
   * \param name Sequence label (can be NULL).
   * \param seqp Sequence of nucleotide codes (can be NULL).
   * \param name_qual Label for quality factors (can be NULL).
   * \param qualp Quality factors as characters (can be NULL).
   */

  int seqFastqSetQual(SeqFastq *sqp, const char qval);
  /**< Set the base qualities to a single value.
   * \param sqp Structur to be modified.
   * \param qval Value to which base qualities should be set
   */

  int seqFastqAppendSegment(SeqFastq *top, const SeqFastq *fromp, 
			    SEQLEN_t start, SEQLEN_t length,
			    char reverse, const SeqCodec *codep);
  /**< Append a segment from a SeqFastq structure to the end of an existing 
   * sequence (which may be empty, e.g. after calling seqFastqBlank().
   * \param top   Target structure to which the segment is appended.
   * \param fromp Source structure from which the segment is taken.
   * \param start Start position of the segment (counting from 0).
   * \param length Length of the segment. If == 0: append entire sequence.
   * \param reverse If != 0 reverse the segment.
   * \param codep  En/Decoder. If != NULL return complement.
   */

  int seqFastqReverse(SeqFastq *sqp, const SeqCodec *codecp);
  /**< Reverse (complement) sequence.
   * \param sqp Fastq sequence.
   * \param codecp Sequence de-/encoder. Can be NULL in which case sequence
   *            is reversed without creating the complement.
   */
  int seqFastqRead(SeqFastq *sqp, SeqIO *ifp);
  /**< Read sequence data in FASTA/FASTQ format.
   * If the sequence data is followed by EOF, ERRCODE_SUCCESS is returned but
   * the status flag of ifp is set to ERRCODE_EOF.
   *
   * \param sqp Structure to be loaded.
   * \param ifp Data source
   */

  int seqFastqFind(SeqFastq *sqp, const char *nam, SeqIO *ifp);
  /**< Find a sequence by (proper) name in header and read in FASTA/FASTQ format
   * If the sequence data is followed by EOF, ERRCODE_SUCCESS is returned but
   * the status flag of ifp is set to ERRCODE_EOF.
   *
   * \param sqp Structure to be loaded.
   * \param nam Sequence name.
   * \param ifp Data source
   */
   
  int seqFastqWrite(SeqIO *ofp, const SeqFastq *sqp, short linewidth);
  /**< write sequence which must be in ASCII encoding in FASTA/FASTQ format.
   * \param ofp Output stream
   * \param sqp Sequence to be written
   * \param linewidth Line width, if == 0 sequence is written in one line.
   */

  int seqFastqWriteCompressedToFile(FILE *fp, const SeqFastq *sqp);
  /**< Write compressed sequence in raw format to file. Nucleotide
   * sequence must be in compressed format. Quality factors are not written.
   */

  const char *seqFastqGetConstSequence(const SeqFastq *sqp, SEQLEN_t *length, char *code);
  /**< Provides access to base sequence raw data.
   * \param      sqp    Sequence structure. 
   * \param[out] length Number of bases in sequence.
   * \param[out] code   The type of encoding used - one of SEQ_CODES.
   */

  char *seqFastqGetSequence(SeqFastq *sqp, SEQLEN_t *length, char *code);
  /**< Provides access to base sequence raw data, so that it can be changed */

  const char *seqFastqGetConstQualityFactors(const SeqFastq *sqp, SEQLEN_t *length, char *code);
  /**< Provides access to sequence of quality factors raw data.
   * \param      sqp    Sequence structure. 
   * \param[out] length Number of quality factors in sequence.
   * \param[out] code   The type of encoding used - one of SEQ_CODES.
   */

  char *seqFastqGetQualityFactors(SeqFastq *sqp, SEQLEN_t *length, char *code);
  /**< Provides access to sequence of quality factors, so that it can be changed */

  const char *seqFastqGetSeqName(const SeqFastq *sqp);
  /**< Accessor for sequence label */

  void seqFastqCurtailSeqName(SeqFastq *sfqp);
  /**< Terminates the header after the first white space (name proper).
   * Subsequent calls of seqFastqGetSeqName will return only the name proper
   */

  const char *seqFastqGetQualName(const SeqFastq *sqp);
  /**< Accessor for quality label */

  int seqFastqEncode(SeqFastq *sqp, const SeqCodec *codep);
  /**< Encode nucleotide sequence in-place */

  int seqFastqDecode(SeqFastq *sqp, const SeqCodec *codep);
  /**< Decode the entire sequence back to ASCII symbols */
  
  int seqFastqCompress(SeqFastq *sqp);
  /**< Compress sequence to 3 bits/base, 10 bases per 32bit int. 
   * Sequence must be encoded (SEQCOD_MANGLED)*/

  int seqFastqUncompress(SeqFastq *ucp, const SeqFastq *sqp,
			 SEQLEN_t start, SEQLEN_t length, 
			 const SeqCodec *codep, char as_rcp);
  /**< Return a segment of a sequence uncompressed 
   * \param ucp Sequence containing the uncompressed segment.
   * \param sqp Compressed sequence from which the segment is taken.
   * \param start Start of the segment in bases (starting from 0).
   * \param length Length of the segment in bases. If 0 return the entire
   *              sequence uncompressed.
   * \param codep Decoder used for (en)coding/(un)compression.
   * \param as_rcp if != 0 return reverse complement
   */

  int seqFastqReadCompressedBinary(SeqFastq *sqp, FILE *fp, const char *label);
  /**< Read a sequence of bases from a binary file in compressed form 
   * (see SEQCOD_COMPRESSED) and label it.
   * \param sqp Sequence container.
   * \param fp input stream.
   * \param label Label the sequence with this name (can be NULL).
   */

  int seqFastqReadCompressedBinaryOfKnownLength(SeqFastq *sqp, FILE *fp, 
						SEQLEN_t len, 
						const char *label);
    /**< Read a sequence of len bases from a binary file in compressed form 
   * (see SEQCOD_COMPRESSED) and label it.
   * \param sqp Sequence container.
   * \param fp input stream.
   * \param len Number of bases in the sequence.
   * \param label Label the sequence with this name (can be NULL).
   */
     
  int seqFastqDecodeAsStandardNt(SeqFastq *dep, const SeqFastq *sqp, 
				 SEQLEN_t start, SEQLEN_t length, 
				 const SeqCodec *codep,
				 char as_rcp);
  /**< Return a segment of an encoded (SEQCOD_MANGLED) sequence cast as standard nucleotides
   * \param dep Sequence structure containing the decoded segment
   * \param sqp Encoded sequence from which the segment is taken.
   * \param start Start of the segment in bases (starting from 0).
   * \param length Lengh of the segment. If 0 decode the entire sequence.
   * \param codep Decoder used for en/decoding.
   * \param as_rcp If != 0 return reverse complement
   */

  /******************************************************************************
   *************************** Methods of Type SeqSet ***************************
   ******************************************************************************/

  /* SeqSet holds a set of sequences in compressed form in
   * memory. Bases can be accessed by specifying either the sequence
   * serial number together with the offset in the sequence or by
   * specifying the offset in the sequence of concatenated sequences.
   */

  SeqSet *seqSetCreate(int blocksiz, unsigned char flags);
  /**< Constructor.  
   * \param blocksiz Block size for memory allocation. 
   * \param flags A combination of SEQSET_FLAGS, e.g
   *        SEQSET_TERMCHAR sets up offsets allowing for an extra termination character,
   *        SEQSET_BASQUAL stores base qualities.
   */

  void seqSetDelete(SeqSet *ssp);
  /**< Destructor.
   */

  void seqSetBlank(SeqSet *ssp);
  /**< Empty sequence set.
   */

  int seqSetAddSequence(SeqSet *ssp, const SeqFastq *sqp);
  /**< Add a sequence to the Set.
   * \param ssp SeqSet to which the sequence should be added.
   * \param sqp Sequence to be added (must be uncompressed).
   */
  int seqSetCompress(SeqSet *ssp, const SeqCodec *codecp);
  /**< Compress the set of sequences. No more sequences can be added.
   * \param ssp Sequence set.
   * \param codecp endoder.
   */
  int seqSetWriteBinFil(const SeqSet *ssp, const char *filnam);
  /**< Write a set of sequences compressed to binary file.
   * \param ssp Sequence set to be written.
   * \param filnam Root file name of binary file. 
   *                Routine adds a constant extension.
   */
  SeqSet *seqSetReadBinFil(int *errcode, const char *filnam);
  /**< Read a set of sequences compressed from binary file.
   * \param errcode Returns one of ERRMSG_CODES
   * \param filnam Root file name of binary file. 
   *                Routine adds a constant extension.
   */
  int seqSetAddFromFastqFile(ErrMsg *errmsg, SeqSet *ssp, 
			     SeqFastq *sqbufp, const SeqCodec *codecp, 
			     const char *filnam, 
			     char verbose);
  /**< Read set of sequences from a FASTQ file.
   */
  SEQNUM_t seqSetGetOffsets(const SeqSet *ssp, const SETSIZ_t **soffs);
  /**< Accessor for the sequence offsets. Returns the number of sequences n_seq.
   * \param ssp Sequence set.
   * \param soffs Returns sequence offsets (can be NULL). *soffs[i] is the position 
   *              of the 1st base of sequence i if all seqeuences were concatenated. 
   *              soffs[n_seq] is the total number of bases.
   */
  int seqSetFetchSegment(SeqFastq *sqp,
			 SETSIZ_t *offs_start, SETSIZ_t *offs_end, 
			 const SeqSet *ssp, const SeqCodec *codecp);
  /**< Return a segment by offset in the sequence of concatenated sequences.
   * \param sqp Structure in which the segment is returned. If *ssp flags SEQSET_BASQUAL
   *            and *sqp is of type SEQTYP_FASTQ, base qualities are returned as well. The
   *            Sequence name remains unchanged.
   * \param offs_start Offset of the segment start in the sequence of concatenated sequences
   *                   (Counting from 0). If the segment stretches sequence boundaries, 
   *                   The offset gets updated to beginning of the first sequence.
   * \param offs_end Offset of the last base of the segment in the sequence of concatenated 
   *                  sequences (Counting from 0). Gets updated if segment spans boundaries 
   *                  of sequences consituting the set.
   * \param ssp Sequence set.
   * \param codecp Encoder/Decoder 
   *
   * \note If the sequences in the set are compressed. The segment is returned as
   * ASCII.
   */
  int seqSetFetchSegmentBySequence(SeqFastq *sqp, SEQNUM_t seqidx,
				   SEQLEN_t offs, SEQLEN_t len, 
				   const SeqSet *ssp, const SeqCodec *codecp);
  /**< Return a segment by sequence index and offset. The segment is returned
   * as ASCII code.
   * \param sqp Structure in which the segment is returned. If *ssp flags SEQSET_BASQUAL
   *            and *sqp is of type SEQTYP_FASTQ, base qualities are returned as well. The
   *            Sequence name remains unchanged.
   * \param seqidx Sequence index of the segment.
   * \param offs Offset of the segment in sequence seqidx.
   * \param len Length of the segment, the returned segment is.
   *            truncated at the sequence end. If 0 return sequence
   *            up to end.
   * \param ssp Sequence set.
   * \param codecp Sequence Encoder/Decoder.
   */
  int seqSetGetIndexAndOffset(SEQNUM_t *seqidx, SEQLEN_t *seqoffs,
			      SETSIZ_t offs, const SeqSet *ssp);
  /**< Return the sequence index and offset in the sequence for the offset in the 
   * concatenated sequences.
   * \param seqidx Returns sequence index (counting from 0).
   * \param seqoffs Returns offset of sequence seqidx in the concatenated
   *        sequence (counting from 0).
   * \param offs Offset in the sequence of concatenated sequences.
   * \param ssp Sequence set.
   */

  SEQLEN_t seqSetGetSeqDatByIndex(SETSIZ_t *offs, const char **name, 
				  SEQNUM_t seqidx, const SeqSet *ssp);
  /**< Accessor for sequence offset, length, name by sequence index. Returns sequence
   * length.
   * \param offs returns sequence offset (can be NULL).
   * \param name return sequence name (can be NULL).
   * \param seqidx Sequence number for which name should be returned.
   * \param ssp Sequence set.
   */

  SEQNUM_t seqSetGetSeqNumAndTotLen(SETSIZ_t *totseqlen, const SeqSet *ssp);
  /**< Return the number of sequences and the total sequnce length.
   * \param totseqlen Returns the total length of all sequences concatenated.
   * \param ssp Sequence set.
   */

#endif
#ifdef __cplusplus
}
#endif
