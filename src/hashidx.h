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

#ifdef __cplusplus
extern "C"
{
#endif

#ifndef HASHIDX_H
#define HASHIDX_H

  //#define hashidx_debug
  //#define hashidx_debug_init

#include <stdint.h>
#include "sequence.h"
#include "interval.h"

  /****************************************************************************
   ****************************** Constants ***********************************
   ****************************************************************************/
  enum HASHIDX_TYPES {         /**< Hash index types */
    HASHIDXTYP_PERFECT = 0,    /**< Perfect hash index (no collisions) */
    HASHIDXTYP_HASH32MIX = 1,  /**< Hash index with collisions and 32-bit mixer */
  };

  enum HASH_STORAGE_MODES { /**< Hash table storage modes */
    HASHDISK_NONE = 0,    /**< Entire hash table is kept in memory */
    HASHDISK_BODY = 1,    /**< The body of the table is left on disk */
    HASHDISK_ALL  = 2     /**< Everthing except the hash table head is left on disk */
  };

  /****************************************************************************
   ***************************** Simple types *********************************
   ****************************************************************************/

  typedef uint32_t HASHNUM_t;  /**< Holds the number of hits */
  typedef uint32_t HASHPOS_t;  /**< Holds the position of a hit (as multiples 
				   * of sampling step (stride) */
  typedef uint64_t HASHWORD_t;    /**< holds an encoded k-mer (word) */

  enum HASH_LIMITS {
    HASHPOS_MAX = UINT32_MAX,
  };

  /****************************************************************************
   ***************************** Opaque types *********************************
   ****************************************************************************/

  typedef struct _HashTable HashTable;

  /****************************************************************************
   *********************** Methods of type HashTable **************************
   ****************************************************************************/

  HashTable *hashTableCreate(uint8_t wordlen, uint8_t nskip,
			     uint8_t nbits_key, uint8_t nbits_perf,
			     uint8_t typ);
  /**< Constructor.
   * \param wordlen Number of bases in word.
   * \param nskip Skip step size, i.e. words are sampled every nskip bases.
   * \param nbits_key Number of bits for the hash key, if == wordlen*2 a perfect
   *        hash function is used.
   * \param nbits_perf Number of bits for the perfect part of the hash key.
   *        nbits_perf < nbits_key. Is ignored if nbits_key == wordlen*2 (perfect hash).
   * \param typ Hash index tupe (one of HASHIDX_TYPES).
   */
  void hashTableDelete(HashTable *htp);
  /**< Destructor 
   */
  void hashTableReset(HashTable *htp, unsigned char nskip);
  /**< Initialise word position table to 0. Set new stride.
   * \param htp Hash table.
   * \param nskip If > 0 set new k-mer sampling stride.
   */

  int hashTableSetUp(HashTable *htp, SeqFastq *sqbufp, const SeqSet *ssp,
		     const InterVal *ivp,
		     const SeqCodec *codecp, 
		     uint32_t *npos_max,
		     char verbose);
  /**< Calculate hash table for a set of sequences. Compress the sequences.
   * \param htp Hash table.
   * \param sqbufp Sequence buffer.
   * \param ssp Set of sequences.
   * \param ivp Set of intevals (segments in the sequences ssp) to which the
   *        hash table should be restricted (can be NULL).
   * \param codecp En-/Decoder for sequences.
   * \param npos_max Returns the number of k-mer positions. Can be NULL.
   *        If pointing to a value > 0 on input the function returns  
   *        ERRCODE_MAXKPOS if the number of elements in the k-mer 
   *        position table is greater than that value. No additional
   *        memory is allocated in that case.
   * \param verbose Flag causing output of messages about progress on 
   *        standard output.
   */
  int hashTableCheckExtensive(SeqFastq *sqbufp, const HashTable *htp,
			      const SeqSet *ssp, const SeqCodec *codecp);
  /**< Check the hash table for consistency by going through k-mer words
   * in the sequences and looking them up in the table.
   * \return ERRCODE_BROKENHASH if discrepancies are found
   * \return ERRCODE_SUCCESS otherwise.
   * \param sqbufp Sequence buffer.
   * \param htp Hash table.
   * \param ssp Set of sequences.
   * \param codecp En-/Decoder for sequences.
   */
  int hashTableCheckQuick(SeqFastq *sqbufp, const HashTable *htp,
			  const SeqSet *ssp, const SeqCodec *codecp);
  /**< Consistency test that goes through all kmer words in the table
   * and checks for each position whether there is the correct k-mer
   * word in the hashed sequence.
   * \return ERRCODE_BROKENHASH if discrepancies are found
   * \return ERRCODE_SUCCESS otherwise.
   * \param sqbufp Sequence buffer.
   * \param htp Hash table.
   * \param ssp Set of sequences.
   * \param codecp En-/Decoder for sequences.
   */

  uint8_t hashTableGetKtupLen(const HashTable *htp, uint8_t *nskip);
  /**< Return the word length of the hash table
   * \param htp Hash table.
   * \param nskip Returns skip step size (can be NULL)
   */

  uint32_t hashTableGetMaxPos(const HashTable *htp);
  /**< Return the maximum k-mer word position
   */

  void hashTablePrintStats(FILE *fp, const HashTable *htp);
  /**< Print hash table statistics.
   */

  int hashTableCmp(const HashTable *htp, const HashTable *ht2p);
  /**< Compare contents of two hash tables.
   */

  HASHNUM_t hashTableGetKtupleHits(HASHPOS_t **posp, HASHNUM_t *posidx, 
				   const HashTable *htp, HASHWORD_t ktup);
  /**< Accessor returning the number of kmer word hits for a k-mer.
   * \return Number of kmer word hits.
   * \param[out] posp Returns pointer to raw tuple positions (can be NULL).
   * \param[out] posidx Returns the index to the block of ktuple positions (can be NULL).
   * \param htp Hash table.
   * \param ktup encoded tuple.
   */
  
  
  HASHNUM_t hashTableFetchHitPositions(HASHPOS_t **posp,
				       const HashTable *htp, HASHNUM_t posidx);
  /**< Accessor returning the raw tuple positions.
   * \return Number of hits.
   * \param[out] posp Return pointer to raw tuple positions.
   * \param htp Hash table.
   * \param[out] posidx Index to the block of ktuple positions 
   *             (e.g. as returned by hashTableGetKtupleHits).
   */

  int hashTableWrite(const char *filnam, const HashTable *htp);
  /**< Write hash table to file filnam
   */

  HashTable *hashTableRead(int *errcode, const char *filnam);
  /**< Read hash table from file filnam. Memory is allocated.
   */

#endif
#ifdef __cplusplus
}
#endif
