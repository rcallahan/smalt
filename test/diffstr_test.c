/** Testing handling of compressed alignment strings
 */

/*****************************************************************************
 *****************************************************************************
 *                                                                           *
 *  Copyright (C) 2013 - 2014 Genome Research Ltd.                           *
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

#include "elib.h"
#include "diffstr.h"

enum {
  DIFFSTR_BLKSZ = 64,
  CIGAR_MAXLEN = 256,
};

int main(int argc, char *argv[]) 
{
  int errcode = ERRCODE_SUCCESS;
  char *cigstr;
  char cigarout[CIGAR_MAXLEN];
  const unsigned char isSAMCIGAR = 1;
  const char isSoftClipped = 1;
  int clip_start, clip_end, nchar=0;
  DiffStr *dfsp, *dfsbfp;
  ErrMsg *errmsgp = 0;

  if (argc != 2) {
    printf("usage: %s <GIGAR string>\n", argv[0]);
    exit(1);
  }

  cigstr = argv[1];

  ERRMSG_CREATE(errmsgp);

  dfsp = diffStrCreate(DIFFSTR_BLKSZ);
  if (dfsp == NULL)
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  dfsbfp = diffStrCreate(DIFFSTR_BLKSZ);
  if (dfsbfp == NULL)
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);
  
  if ((errcode = diffStrParseSimul(dfsp, &clip_start, &clip_end, isSAMCIGAR, cigstr)))
    ERRMSGNO(errmsgp, errcode);

  printf("CIGAR_RAW ");
  if ((errcode = diffStrPrintf(stdout, dfsp->dstrp, DIFFSTRFORM_RAW, 
			       clip_start, clip_end, isSoftClipped)))
    ERRMSGNO(errmsgp, errcode);

  printf("CIGAR_SSAHA ");
  if ((errcode = diffStrPrintf(stdout, dfsp->dstrp, DIFFSTRFORM_CIGNORM, 
			       clip_start, clip_end, isSoftClipped)))
    ERRMSGNO(errmsgp, errcode);
  
  printf("\nCIGAR_EXT ");

  if ((errcode = diffStrPrintf(stdout, dfsp->dstrp, DIFFSTRFORM_CIGEXT, 
			       clip_start, clip_end, isSoftClipped)))
    ERRMSGNO(errmsgp, errcode);
 
  printf("\nCIGAR_SAM ");

  if ((errcode = diffStrPrintfStr(cigarout, &nchar, dfsp->dstrp, DIFFSTRFORM_CIGEXT, 
				  clip_start, clip_end, isSoftClipped)))
    ERRMSGNO(errmsgp, errcode);

  printf("%s\nCIGAR_SAMX ", cigarout);

  if ((errcode = diffStrPrintfStr(cigarout, &nchar, dfsp->dstrp, DIFFSTRFORM_CIGEXT_XMISMATCH, 
				  clip_start, clip_end, isSoftClipped)))
    ERRMSGNO(errmsgp, errcode);

  printf("%s\n", cigarout);
  
  diffStrDelete(dfsp);

  ERRMSG_END(errmsgp);

  exit(errcode);
}
