/**< Buffers for pairwise alignments */

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

#ifndef ALIBUFFER_H
#define ALIBUFFER_H

  /****************************************************************************
   ***************************** Opaque Types *********************************
   ****************************************************************************/

  typedef struct _AliBuffer AliBuffer;

  /****************************************************************************
   ************************* Methods of Type AliBuffer ************************
   ****************************************************************************/
  
  AliBuffer *aliBufferCreate(int blocksiz);
  /**< Constructor.
   * \param blocksiz Block size (granularity) as the number of elements per 
   *       buffer for memory allocation.
   */

  void aliBufferDelete(AliBuffer *p);
  /**< Destructor
   */

  int aliBufferInit(AliBuffer *p, unsigned int qlen);
  /**< Initialise buffer and if neccessary reallocate memory for 
   * query length qlen.
   */

#endif /* ifndef ALIBUFFER_H */
#ifdef __cplusplus
}
#endif
