/** Command line parser */

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

#ifndef MENU_H
#define MENU_H

#include <time.h>
#include <stdint.h>
  /****************************************************************************
   ****************************** Constants ***********************************
   ****************************************************************************/

  enum MENU_SUBPROG {
    MENU_UNKNOWN,
    MENU_INDEX,    /**< build a binary hash table */
    MENU_MAP,      /**< map reads */
    MENU_CHECK,    /**< map reads, sample insert sizes */
    MENU_SAMPLE,
    MENU_HELP,
    MENU_VERSION,
  };

  enum MENU_FLAGS {
    MENUFLAG_BINHTAB  = 0x01,      /**< if set, a binary hash table is read from disk */
    MENUFLAG_COMPLEXW = 0x02,      /**< if set, Smith-Waterman alignment scores are 
				    * complexity weighted */
    MENUFLAG_VERBOSE  = 0x04,      /**< if set, the output is complemented by comments
				    * about various stages of the program */
    MENUFLAG_SENSITIVE = 0x08,     /**< if submit all potentially matching regions to Smith-Waterman
				    * otherwise, various heuristics are used to curtail the number 
				    * of Smith-Waterman alignments */
    MENUFLAG_PAIRED     = 0x10,     /**< if set expect paired-end reads (for FASTQ input format only) */
    MENUFLAG_HASHDISK   = 0x20,     /**< if set leave body of the hash table on disk */
    MENUFLAG_REFSEQDISK = 0x40,     /**< if set leave the reference sequences on disk */
    MENUFLAG_ALIGNMENT  = 0x80,     /**< if set output explicit alignment */
    MENUFLAG_SPLITREAD  = 0x100,    /**< report partial alignments if they complement each other */
    MENUFLAG_EXHAUSTIVE = 0x200,    /**< perform exhaustive search of mate pairs (MENUFLAG_PAIRED) */
    MENUFLAG_MINSCOR    = 0x400,    /**< a minimum score was specified on the command line */
    MENUFLAG_RELSCOR    = 0x800,    /**< a minimum score relative to the maximum (-d) was specified
                                     * on the command line */
    MENUFLAG_RANDREPEAT = 0x1000,   /**< Report repetitive mappings by picking one at random */
    MENUFLAG_MATEPAIRLIB = 0x2000,  /**< Paired reads are from an illumina mate-pair library (for
				     * long insert sizes) */
    MENUFLAG_READORDER = 0x4000,    /**< Preserve the read order in the output */
  };
 
  enum MENU_OUTPUT_FORMATS {       /**< Output formats */
    MENU_OUTFORM_CIGAR,            /**< Cigar format (default) */
    MENU_OUTFORM_SAM,              /**< SAM (samtools.sourceforge.net) */
    MENU_OUTFORM_BAM,              /**< BAM (binary SAM) */
    MENU_OUTFORM_SSAHA,            /**< SSAHA2 default format */
    MENU_OUTFORM_GFF,              /**< General Feature Format */
  };

  enum MENU_OUFMTSAM_FLAGS {       /**< Flags for SAM output format */
    MENU_OUFMTSAM_HEADER = 0x01,   /**< Output SAM header if set */
    MENU_OUFMTSAM_CLIPPED = 0x02,  /**< Hard clip sequences to CIGAR string of aligned segment */
    MENU_OUFMTSAM_XCIGAR = 0x04,   /**< Cigar string with X for mismatch */
  };

  enum MENU_INPUT_FORMATS {        /**< Expected sequence input formats */
    MENU_INFORM_UNKNOWN,           /**< No particular format expected (to be determined) */
    MENU_INFORM_FASTQ,             /**< FASTA/FASTQ format */
    MENU_INFORM_SAM,               /**< SAM  (samtools.sourceforge.net) */
    MENU_INFORM_BAM,               /**< BAM  (binary SAM) */
  };

  enum MENU_READPAIR_TYPES {  /**< Encodes the orientation of mates in different read pair libraries. */ 
    MENU_READPAIRTYP_UNKNOWN = 0,
    MENU_READPAIRTYP_SINGLE = 1,     /**< Single reads, not paired reads */
    MENU_READPAIRTYP_PAIREDEND = 2,  /**< The 2 mates of a pair are on opposite strands with the 5' ends
				      * on the outside and the 3' ends facing each other, i.e. |--><--|, 
				      * like in the illumina paired-end library */ 
    MENU_READPAIRTYP_MATEPAIR = 3,   /**< The 2 mates of a pair are on opposite strands with the 3' ends
				      * on the outside and the 5' end facing each other, i.e. <--| |-->, 
				      * like in the illumina mate-pair library */ 
    MENU_READPAIRTYP_SAMESTRAND = 4, /**< The 2 mates of a pair are on the same strand, i.e. |--> |-->,
				      * like in (some) 454 libraries. */
  };

  /****************************************************************************
   ********************************* Types ************************************
   ****************************************************************************/
  typedef uint8_t MENUOFMTFLG_t;  /**< contains a combination of output 
				    * format flags (MENU_OUFMTSAM_FLAGS) */
  typedef uint16_t MENUFLG_t;

  /****************************************************************************
   ***************************** Opaque types *********************************
   ****************************************************************************/

  typedef struct _MenuOpt MenuOpt;

  /****************************************************************************
   *********************** Methods of type MenuOpt ****************************
   ****************************************************************************/

  MenuOpt *menuCreate(void);
  /**< Constructor 
   */
  void menuDelete(MenuOpt *mp);
  /**< Destructor 
   */
  int menuParseCommandLine(MenuOpt *mp, int argc, char *argv[]);
  /**< Parse the command line
   * \return ERRCODE_SUCCESS
   * \return ERRCODE_FAILURE on error.
   * \param mp Structure to be filled 
   * \param argc number of arguments (as in main() argument list)
   * \param argv Array of tokens (as in main() argument list)
   */
  void menuPrint(FILE *fp, const MenuOpt *mp);
  /**< Output menue options on stream.
     \param fp Output stream.
     \param mp pointer to structure
  */
  
  char menuGetSubProgTyp(const MenuOpt *mp);
  /**< Returns the type of sub program specified (one of MENU_SUBPROG)
   */

  char **menuGetCommandLine(const MenuOpt *mp, int *argc);
  /**< Return number of and pointer to command line arguments
   */
  const char *menuGetProgramName(const char **version);
  /**< Return Name of main program and version string.
   */

  int menuGetHashParams(const MenuOpt *mp, const char **hashnamp, 
			unsigned char *kmer, unsigned char *nskip);
  /**< Returns options for building the hash index
   * \param mp Menue
   * \param hashnamp Returns root name for hash index (can be NULL)
   * \param kmer Returns hashed word length (can be NULL)
   * \param nskip Returns step size with which words are sampled (can be NULL)
   */
  int menuGetMapParams(const MenuOpt *mp, const char **hashnamp, 
		       int *nhitmax_tuple, double *tupcovmin, 
		       short *minscore, short *scorediff,
		       unsigned char *minbasq,
		       double *minidentity,
		       int *randseed,
		       int *readkip,
		       int *insert_min, int *insert_max, 
		       unsigned char *outform,
		       const char **outfilnam,
		       const char **insfilnam,
		       unsigned char *pairtyp);
  /**< Returns options for read mapping 
   * \param mp Menue
   * \param hashnamp Returns root name for hash index (can be NULL)
   * \param nhitmax_tuple K-mer words with more than nhitmax hits across the genome are
   *        disregarded.
   * \param tupcovmin 'Coverage' of read by k-tuples. Represents absolute number
   *         of bases if >0, else fraction of read length.
   * \param minscore Minimum Smith-Waterman score for an alignment to be counted.
   * \param scorediff Threshold of the Smith-Waterman alignment score relative to the
   *        maximum score. All mappings resulting in Smith-Waterman scores within
   *        of the maximum are reported. Mappings with scores lower than this value
   *        are skipped. < 0 signals no threshold set.
   * \param minbasq Base quality threshold. K-mer words with base qualities below this
   *        threshold are not looked up in the hash index.
   * \param minidentity Minum number (fraction) of identical bases (matches) for an alignment
   *        to be reported.
   * \param randseed Seed for pseudo-random generator.
   * \param readskip Step size for skipping reads (when sampling insert sizes).
   * \param insert_min Minimum insert size for paired read mode.
   * \param insert_max Maximum insert size for paired reads mode.
   * \param outform Returns output format (one of MENU_OUTPUT_FORMATS).
   * \param outfilnam Returns name of output file.
   * \param insfilnam Returns name of file with insert size distribution.
   * \param pairtyp Returns the library type (one of MENU_READPAIR_TYPES).
   */

  int menuTestMapOutputFormatFlags(const MenuOpt *menup, 
				   const unsigned char outform, const MENUOFMTFLG_t flags);
  /**< Test combination of input format specific binary flags 
   * \param menup Menue
   * \param outform Output format (one of MENU_OUTPUT_FORMATS)
   * \param flags Comination of binary flags that are specific to the output format outform
   */

  int menuGetMapInputFormat(const MenuOpt *mp, 
			    uint8_t *infmt, const char **tmpdirnam);
  /**< Return the sequence input format and directory name for temporary files.
   * \param Menu for mapping.
   * \param infmt Returns the sequence input format (one of MENU_INPUT_FORMATS). 
   *        Can be NULL;
   * \param tmpdirnam Returns name of directory to which temporary files are written.
   *        (only used when reading SAM/BAM files). Can be NULL.
   */

  uint8_t menuGetMapPenaltyScores(const MenuOpt *menup, int8_t *match, int8_t *subst, 
				  int8_t *gap_open, int8_t *gap_ext);
  /**< Accessor for alignment penalty scores.
   * \return 4-bit field in the least significant bits with bit set if the
   * corresponding penalty score was set on the command line. The bits are,
   * from most significatn to least significant: match, subst, gap_open, gap_ext.
   * \param mp Menu for mapping.
   * \param match Returns score for a match (can be NULL). Bit 00x1 is set in the returned
   * number if this parameter was set on command line.
   * \param subst Returns score for a mismatch (NULL is allowed). Corresponding bit: 0x02.
   * \param gap_open Returns the penalty for opening a gap (NULL is allowed). Corresponding bit: 0x04.
   * \param gap_ext Returns the penalty for extending a gap (NULL is allowed). Corresponding bit: 0x08.
   */

  MENUFLG_t menuGetFlags(const MenuOpt *mp);
  /**< Returns combination of binary flags MENU_FLAGS (mp can be NULL).
   */

  short menuGetNumberOfThreads(const MenuOpt *mp);
  /**< Returns the number of threads to be forked from the main thread.
   */
  int menuGetFileNames(const MenuOpt *mp, 
		       const char **filnam1p, const char **filnam2p);
  /**< Returns the number and array of file names
   * \param mp Menu.
   * \param filnam1p Returns name of 1st read file (can be NULL).
   * \param filanm2p Returns name of 2nd read file (can be NULL).
   */

  void menuPrintWallClockTime(FILE *fp, 
			      time_t time_start, time_t time_stop,
			      const char *headerp);
  /**< Print the wall clock time passed between 2 time points
   * \param fp Output stream.
   * \param time_start start point.
   * \param time_end end point.
   * \param headerp Text for a header line before time is output.
   */
#endif
#ifdef __cplusplus
}
#endif
