/** Binary file I/O */

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

#include <string.h>

#include "elib.h"
#include "filio.h"

enum {
  MASK8BIT            = 0xFF,
  IOFIL_HEADSIZ       = 12,            /**< header size in units of 32 bits */
  IOFIL_SIGNATURE     = 0x73212173,   /**< endianess invariant signature number */
  IOFIL_ENDIANTESTNUM = 0x6E378A19,   /**< inverse indicates swap in endianess */
};

static const char WRITERRMSG[] = "when writing header of binary file";
static const char READERRMSG[] = "when reading header of binary file";


/******************************************************************************
 *********************************** Macros ***********************************
 ******************************************************************************/

#define SWAP(a,b) tmp = (a); (a) = (b); (b) = tmp;


/******************************************************************************
 ******************************* Private Methods ******************************
 ******************************************************************************/

static int writeHeader(uint32_t filsiz, uint8_t filtyp, uint32_t version, uint32_t headsiz, FILE *fp)
{
  int errcode = ERRCODE_SUCCESS;
  short i;
  uint32_t header[IOFIL_HEADSIZ];

  header[0] = IOFIL_SIGNATURE;
  header[1] = IOFIL_ENDIANTESTNUM;
  header[2] = filsiz;  /* file size in units of 32 bit ints (4 bytes) */
  header[3] = filtyp;
  header[4] = version;
  header[5] = headsiz; /* file type specific header size in units of 32-bit words */
  for (i=6; i<IOFIL_HEADSIZ; i++) header[i] = 0;
  
  fwrite(header, sizeof(uint32_t), IOFIL_HEADSIZ, fp);
  
  if (ferror(fp)) {
    errcode = ERRCODE_WRITEERR;
    perror(WRITERRMSG);
  }
  
  return errcode;
}

static int readHeader(uint32_t *filsiz, uint8_t *filtyp, uint32_t *version,
	       uint8_t *is_endianid, uint32_t *headsiz, FILE *fp)
{
  uint32_t header[IOFIL_HEADSIZ];
  size_t nrobj = fread(header, sizeof(uint32_t), IOFIL_HEADSIZ, fp);

  if ( nrobj != IOFIL_HEADSIZ ) {
    if (ferror(fp))
      perror(READERRMSG);
    return ERRCODE_READERR;
  }
 
  if (header[0] != IOFIL_SIGNATURE) 
    return ERRCODE_FILTYP;

  *is_endianid = header[1] == IOFIL_ENDIANTESTNUM;
  if (!(is_endianid)) {
    filioSwapEndian(header, IOFIL_HEADSIZ);
    if (header[1] != IOFIL_ENDIANTESTNUM) 
      return ERRCODE_ENDIAN;
  }
  *filsiz = header[2];
  *filtyp = header[3]&MASK8BIT;
  *version = header[4];
  *headsiz = header[5];

  return ERRCODE_SUCCESS;
}

/******************************************************************************
 ******************************* Public Methods *******************************
 ******************************************************************************/

FILE *filioOpenForWriting(int *errcode, uint32_t siz, uint8_t typ, 
			  uint32_t version, uint32_t headsiz,
			  const uint32_t *header, 
			  const char *filnam, const char *filext)
{
  char fnam[FILENAME_MAX];
  size_t namlen = strlen(filnam);
  size_t extlen = ((filext))? strlen(filext): 0;
  FILE *fp = 0;
  *errcode = ERRCODE_SUCCESS;
   
  if (FILENAME_MAX < (extlen + namlen + 1))
     *errcode = ERRCODE_LONGFILNAM;

  strcpy(fnam, filnam);
  if ((filext)) {
    fnam[namlen] = '.';
    fnam[namlen+1] = '\0';
    strcat(fnam, filext);
    fnam[namlen+extlen+1] = '\0';
  } else {
    fnam[namlen] = '\0';
  }

  if (!(fp = fopen(fnam, "wb"))) {
    *errcode = ERRCODE_NOFILE;
    return 0;
  }

  if ((*errcode = writeHeader(siz+IOFIL_HEADSIZ, typ, version, headsiz, fp))) {
    fclose(fp);
    return 0;
  }

  fwrite(header, sizeof(uint32_t), headsiz, fp);
  
  if (ferror(fp)) {
    perror(WRITERRMSG);
    fclose(fp);
    *errcode = ERRCODE_WRITEERR;
    fp = 0;
  }
  
  return fp;
}
  
FILE *filioOpenForReading(int *errcode, uint8_t *is_endianid, uint32_t *siz, uint8_t *typ, 
			  uint32_t *version, uint32_t *headsiz,
			  uint32_t *header, const char *filnam, const char *filext)
{
  char fnam[FILENAME_MAX];
  size_t namlen;
  size_t extlen = ((filext))? strlen(filext): 0;
  size_t nrobj;
  uint32_t filsiz=0, hs=0;
  FILE *fp = NULL;
   
  if (NULL == filnam) {
    *errcode = ERRCODE_NULLPTR;
    return fp;
  }
   *errcode = ERRCODE_SUCCESS;
   namlen = strlen(filnam);

  if (FILENAME_MAX < (extlen + namlen + 2))
     *errcode = ERRCODE_LONGFILNAM;

  strcpy(fnam, filnam);
  if (extlen > 0) {
    fnam[namlen] = '.';
    fnam[namlen+1] = '\0';
    strcat(fnam, filext);
    fnam[namlen+extlen+1] = '\0';
  } else {
    fnam[namlen] = '\0';
  }

  if (!(fp = fopen(fnam, "rb"))) {
    *errcode = ERRCODE_NOFILE;
    return 0;
  }

  if ((*errcode = readHeader(&filsiz, typ, version, is_endianid, &hs, fp))) {
    fclose(fp);
    return 0;
  }
  
  if (filsiz < IOFIL_HEADSIZ || hs > *headsiz) {
    *errcode = ERRCODE_FHEADSIZ;
    fclose(fp);
    return 0;
  }
  *siz = filsiz - IOFIL_HEADSIZ;

  nrobj = fread(header, sizeof(uint32_t), hs, fp);
  *headsiz = hs;

  if (nrobj != hs) {
    if (ferror(fp))
      perror(READERRMSG);
    fclose(fp);
    *errcode = ERRCODE_READERR;
    fp = 0;
  }
  
  if (!(is_endianid))
    filioSwapEndian(header, hs);

  return fp;
}

void filioSwapEndian(uint32_t *p, uint32_t len)
{
  char *cp;
  register char tmp;

  while (--len > 0) {
    cp = (char *) p++;
    SWAP(cp[1], cp[4]);
    SWAP(cp[2], cp[3]);
  }
}
