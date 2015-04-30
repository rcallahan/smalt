/** Binary file I/O */

/*****************************************************************************
 *****************************************************************************
 *                                                                           *
 *  Copyright (C) 2010 - 2014 Genome Research Ltd.                             * 
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

#ifndef FILIO_H
#define FILIO_H

#include <stdio.h>
#include <stdint.h>

  /****************************************************************************
   ****************************** Constants ***********************************
   ****************************************************************************/

  enum FILIO_TYPES {
    FILIOTYP_UNKNOWN = 0,
    FILIOTYP_SEQSET  = 1,  /**< compressed set of sequences */
    FILIOTYP_HASHTAB = 2,  /**< perfect hash table */
    FILIOTYP_HASHIDX = 3,  /**< hash index with collisions */
    FILIOTYP_BASQUAL = 4,  /**< base quality statistics */
  };

  /****************************************************************************
   ******************************** Methods ***********************************
   ****************************************************************************/

  FILE *filioOpenForWriting(int *errcode, uint32_t siz, uint8_t typ, 
			    uint32_t version, uint32_t headsiz,
			    const uint32_t *header, 
			    const char *filnam, const char *filext);
  /**< Open a binary file for writing and write standard header first.
   * \param errcode Returns one of ERRMSG_CODES.
   * \param siz Size of the file to be written as the number of 32-bit words
   *        excluding the standard header.
   * \param typ File type, one of FILIO_TYPES.
   * \param version Version number of the file format
   * \param headsiz Size of the file type specific header (number of 32-bit words
   *                included in siz).
   * \param header File type specific header.
   * \param filnam Root file name.
   * \param filext File name extension (can be NULL).
   */

  FILE *filioOpenForReading(int *errcode, uint8_t *is_endianid, 
			    uint32_t *siz, uint8_t *typ, 
			    uint32_t *version, uint32_t *headsiz, 
			    uint32_t *header, 
			    const char *filnam, const char *filext);
  /**< Open a binary file for reading.
   * \param errcode Returns one of ERRMSG_CODES.
   * \param is_endianid Flag 0: File was written under different endiannes. 1: File was
   *                            written with the same endiannes as current machine.
   * \param siz Returns the size of the file as number of 32-bit words and excluding
   *            the size of the standard header.
   * \param typ Returns the file type, one of FILIO_TYPES.
   * \param version Returns the file type specific version number.
   * \param headsiz On input: expected size of the file type specific header as the 
   *                number of 32-bit words. On return: Actual size of the file type 
   *                specific header.
   * \param header Returns file type specific header.
   * \param filnam Root file name.
   * \param filext File name extension (can be NULL).
   */

  void filioSwapEndian(uint32_t *p, uint32_t len);
  /**< Swap endiannes of len 32-bit words.
   * \param p Start address.
   * \param len Number of 32-bit words to be swapped.
   */

#endif
#ifdef __cplusplus
}
#endif

