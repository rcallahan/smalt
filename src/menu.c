/** Command line parser */

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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>

#include "elib.h"
#include "menu.h"

enum MENU_OPTION_CONST {
  MENU_OPTION_INDICATOR = '-',
  MENU_PAIR_SEPARATOR = ',',
  MENU_PAIR_RANGENUM = 2,
  TIME24HRSINSECS = 24*60*60,
};

#ifndef PACKAGE_NAME
#define PACKAGE_NAME "smalt"
#endif

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.7.6"
#endif

#ifndef PACKAGE_BUGREPORT
#define PACKAGE_BUGREPORT "hp3@sanger.ac.uk"
#endif
static const char MENU_PROGNAM_SHORT[] = PACKAGE_NAME;
static const char MENU_PROGNAM[] =\
  "              SMALT - Sequence Mapping and Alignment Tool";
static const char MENU_PROGNAM_VERSION_FMT[] = \
  "                             (version: %s)\n";
static const char MENU_RELEASE_VERSION[] = PACKAGE_VERSION;
static const char MENU_RELEASE_DATE[] =  "21-03-2014";
static const char MENU_RELEASE_AUTHORS[] = "Hannes Ponstingl";
static const char MENU_RELEASE_BUGREPORT[] = PACKAGE_BUGREPORT;
static const char MENU_COPYRIGHT_NOTICE[] = "Copyright (C) 2010 - 2014 Genome Research Ltd.";
/* static const char MENU_DEFAULT_HASHNAME[] = "smalt"; */

enum OPTION_TYPES {  /* this is associated with OPTION_TYPSTR */
  OPTYP_FLAG   = 0,
  OPTYP_STRING = 1,
  OPTYP_INT    = 2,
  OPTYP_INT_PAIR = 3,
  OPTYP_FLT = 4
};

/* chage the following in concert with OPTION_TYPES */
static const char *OPTION_TYPSTR[] = {\
  "",
  "STR",
  "INT",
  "INT,INT",
  "FLT"
};
  
typedef struct OPTDOC_ { /**< Documentation for an option */
  char ochr;    /**< option character */
  uint8_t otyp; /**< one of OPTION_TYPES */
  const char *vnam;  /**< Variable name */
  const char *sdesc; /**< short documentation */
  const char *ldesc; /**< long documentation */
} OPTDOC;

typedef struct TASKDOC_ { /**< Documentation for an task */
  const char *synopsis;
  const char *description;
  const OPTDOC *optdoc;
} TASKDOC;

/*****************************************************************************
 ******************************* MENU SUMMARY ********************************
 *****************************************************************************/

static const char MENU_SHORT_DESCRIPTION[] =\
"  Smalt is a pairwise sequence alignment program designed for the mapping of\n" \
"  DNA sequencing reads onto genomic reference sequences.\n"\
"  Running the software involves two steps. First, an index of short words\n"\
"  has to be built for the set of genomic reference sequences (issue \n"\
"  'smalt index -H' for help). Then the sequencing reads are mapped onto the\n" \
"  reference ('smalt map -H' for help).\n\n";

static const char MENU_USAGE_SUMMARY[] =\
"SYNOPSIS:\n"\
"    smalt <task> [TASK_OPTIONS] [<index_name> <file_name_A> [<file_name_B>]]\n\n"\
"Available tasks:\n"\
"    smalt check   - checks FASTA/FASTQ input\n"\
"    smalt help    - prints a brief summary of this software\n"\
"    smalt index   - builds an index of k-mer words for the reference\n"\
"    smalt map     - maps single or paired reads onto the reference\n"\
"    smalt sample  - sample insert sizes for paired reads\n"\
"    smalt version - prints version information\n\n"\
"Help on individual tasks:\n"\
"    smalt <task> -H\n\n";
//"    smalt prep    - checks/splits input files, samples insert sizes\n"


/*****************************************************************************
 ********************************* INDEX TASK ********************************
 *****************************************************************************/

static const OPTDOC MENU_OPTDOC_INDEX[] = {
  {
  ochr:'H',
  otyp:OPTYP_FLAG,
  vnam:"",
  sdesc: "Print more extensive help on options.",
  ldesc:\
  "     Print these instructions.\n"
  },{
  ochr:'k',
  otyp: OPTYP_INT,
  vnam: "wordlen",
  sdesc: "Length of the k-mer words indexed.",
  ldesc:\
  "     Specifies the word length. <wordlen> is an integer within the limits\n"	\
  "     2 < wordlen <= 20. The default word length is 13.\n"
  },{
  ochr:'s',
  otyp: OPTYP_INT,
  vnam: "stepsiz",
  sdesc: "Sample every <stepsiz>-th k-mer word (stride).",
  ldesc: \
  "     Specifies how many bases are skipped between indexed words. With '-s 1'\n"\
  "     every k-mer word along the reference sequences is indexed. With '-s 2'\n"\
  "     every other word is indexed etc. By default the step size is set equal\n"\
  "     to the word length (tiling words).\n"
  },
  {0,0,0,0,0}
};

static const TASKDOC MENU_TASKDOC_INDEX ={
 synopsis:
 "  smalt index [-k <wordlen>] [-s <stepsiz>]  <index_name> <reference_file>\n",
 description:
 "  Generates an index of k-mer words for the genomic reference sequences. The\n" \
 "  words are of fixed length <wordlen> and are sampled at equidistant steps\n" \
 "  <stepsiz> bases apart. The reference sequences are provided in a single\n"
 "  file <reference_file> in FASTA or FASTQ format.\n"			\
 "  Two binary files are output. The file <index_name>.sma contains the \n" \
 "  reference sequences in compressed form. The file <index_name>.smi contains\n" \
 "  the k-mer word index.\n",
 optdoc: MENU_OPTDOC_INDEX,
};

/*****************************************************************************
 ********************************** MAP TASK *********************************
 *****************************************************************************/

#ifdef HAVE_BAMBAMC
static const char MENU_USAGE_MAP_HEADER[] =\
"SYNOPSIS:\n"\
"  smalt map [OPTIONS] <index_name> <query_file> [<mate_file>]"\
"DESCRIPTION:\n"\
"  Map query reads onto the reference sequences. The reads are provided in\n"\
"  FASTA/FASTQ format or in SAM/BAM format in the file <query_file>. If the\n"\
"  name of a second file <mate_file> is specified, both files are in\n"\
"  FASTA/FASTQ format and reads are mapped in pairs. If <query_file> is in\n"\
"  SAM/BAM format, single reads and paired reads can be mixed.\n\n"\
"  The reference sequences and k-mer word index are read from the binary\n"\
"  files <index_name>.sma and <index_name>.smi which must have been created\n"\
"  by the 'index' task (type 'smalt index -H' for help).\n"; 
#else
static const char MENU_USAGE_MAP_HEADER[] =\
"SYNOPSIS:\n"\
"  smalt map [OPTIONS] <index_name> <query_file> [<mate_file>]\n\n"\
"DESCRIPTION:\n"\
"  Map query reads onto the reference sequences. The reads are provided in\n"\
"  FASTA/FASTQ format in the file <query_file>. If the name of a second file\n"\
"  <mate_file> is specified, the reads in both files are mapped in pairs.\n"\
"  The reference sequences and k-mer word index are read from the binary\n"\
"  files <index_name>.sma and <index_name>.smi which must have been created\n"\
"  by the 'index' task (type 'smalt index -H' for help).\n";

#endif
static const OPTDOC MENU_OPTDOC_MAP[] = {
  {
  ochr:'a',
  otyp:OPTYP_FLAG,
  vnam:"",
  sdesc: "Add explicit alignments to output.",
  ldesc:\
  "     Output explicit alignments along with the mapping coordinates.\n"
  },
  {
  ochr:'c',
  otyp:OPTYP_INT,
  vnam:"mincover",
  sdesc: "Threshold of the number of bases covered by k-mer seeds.",
  ldesc:\
  "     Only consider mappings where the k-mer word seeds cover the query read to\n"\
  "     a minimum extent. If <mincover> is an integer or floating point > 1.0, at\n"\
  "     least this many bases of the read must be covered by k-mer word seeds. If\n"\
  "     <mincover> is a floating point <= 1.0, it specifies the fraction of the\n"\
  "     query read length that must be covered by k-mer word seeds. This option\n"\
  "     is only valid in conjunction with the '-x' flag.\n"
  },
  {
  ochr:'d',
  otyp:OPTYP_INT,
  vnam:"scordiff",
  sdesc: "Threshold of the Smith-Waterman score relative to best.",
  ldesc:\
  "     Set a threshold of the Smith-Waterman alignment score relative to the\n"\
  "     maximum score. When mapping single reads, all alignments are reported\n" \
  "     that have Smith-Waterman scores within <scorediff> of the maximum.\n" \
  "     Mappings with lower scores are skipped. If <scorediff> is set to to a\n" \
  "     value < 0, all alignments are printed that have scores above the\n" \
  "     threshold specified with the '-m <minscor>' option.\n"		\
  "     For paired reads, only a value of 0 is supported. With the option '-d 0'\n" \
  "     all aligments (pairings) with the best score are output. By default \n"	\
  "     (without the option '-d 0') single reads/mates with multiple best mappings\n" \
  "     are reported as 'not mapped'.\n"				\
  },
  {
  ochr:'f',
  otyp:OPTYP_STRING,
  vnam:"ouform",
#ifdef HAVE_BAMBAMC
  sdesc: "Output format [sam(default)|bam|cigar|gff|ssaha].\n"\
  "           Ext: [sam|bam]:nohead,x,clip.",
  ldesc:\
  "     Specifies the output format. <ouform> can be either 'sam'(default), 'bam',\n"\
  "     'cigar', 'gff' or 'ssaha'. Optional extension '[sam|bam]:nohead,x,clip'\n" \
  "     (see manual).\n"
#else
  sdesc: "Output format [sam|cigar|gff|saha]. Extension sam:nohead,clip.",
  ldesc:								\
  "     Specifies the output format. <ouform> can be either 'sam'(default),\n"\
  "     'cigar', 'gff' or 'ssaha'. Optional extension 'sam:nohead,x,clip'\n" \
  "     (see manual). Support for BAM format is dependent on additional\n"\
  "     libraries (not installed).\n"
#endif
  },
  {
  ochr:'F',
  otyp:OPTYP_STRING,
  vnam:"inform",
#ifdef HAVE_BAMBAMC
  sdesc: "Input format [fastq (default)|sam|bam].",
  ldesc:\
  "     Specifies the input format. <inform> can be either 'fastq' (default),\n" \
  "     'sam' or 'bam' (see: samtools.sourceforge.net).\n"
#else
  sdesc: "Input format. fastq (default) is the only available format.",
  ldesc:\
  "     Specifies the input format. The only available format is fastq (default).\n"\
  "     Support for BAM and SAM formats (see: samtools.sourceforge.net) depends\n"\
  "     on additional libraries (not installed).\n" 
#endif 
  },
  {
  ochr:'g',
  otyp:OPTYP_STRING,
  vnam:"insfil",
  sdesc: "Reads insert size distribution from file (see 'sample' task).",
  ldesc:\
  "     Use the distribution of insert sizes stored in the file <insfil>. This\n" \
  "     file is in ASCII format and can be generated using the 'sample' task see\n" \
  "     'smalt sample -H' for help).\n"
  },
  {
  ochr:'H',
  otyp:OPTYP_FLAG,
  vnam:"",
  sdesc: "Print more extensive help on options.",
  ldesc:\
  "     Print these instructions.\n"
  },
  {
  ochr:'i',
  otyp:OPTYP_INT,
  vnam:"insert_max",
  sdesc: "Maximum insert size for paired reads (default: 500).",
  ldesc:\
  "     Maximum insert size (only in paired-end mode). The default is 500.\n"
  },
  {
  ochr:'j',
  otyp:OPTYP_INT,
  vnam:"insert_min",
  sdesc: "Minimum insert size for paired reads (default: 0).",
  ldesc:\
  "     Minimum insert size (only in paired-end mode). The default is 0.\n"
  },
  {
  ochr:'l',
  otyp:OPTYP_STRING,
  vnam:"pairtyp",
  sdesc: "Type of paired read library [pe|mp|pp] (default: pe).",
  ldesc:\
  "     Type of read pair library. <pairtyp> can be either 'pe', i.e. for\n" \
  "     the Illumina paired-end library for short inserts (|--> <--|). 'mp'\n"\
  "     for the Illumina mate-pair library for long inserts (<--| |-->) or\n"\
  "     'pp' for mates sequenced on the same strand (|--> |-->). 'pe' is the\n"
  "     default.\n"
  },
  {
  ochr:'m',
  otyp:OPTYP_INT,
  vnam:"minscor",
  sdesc: "Threshold of alignment score.",
  ldesc:\
  "     Sets an absolute threshold of the Smith-Waterman scores. Mappings with\n" \
  "     scores below that threshold will not be reported. The default is\n" \
  "     <minscor> = <wordlen> + <stepsiz> - 1.\n" 
  },
  {
  ochr:'n',
  otyp:OPTYP_INT,
  vnam:"nthreads",
  sdesc: "Number of threads.",
  ldesc:\
  "     Run smalt using mutiple threads. <nthread> is the number of additional\n"\
  "     threads forked. The order of the reads in the input files is not preserved\n" \
  "     for the output unless '-O' is also specified.\n"
  },
  {
  ochr:'o',
  otyp:OPTYP_STRING,
  vnam:"oufilnam",
  sdesc: "Write aligments to specified file (default: stdout).",
  ldesc:\
  "     Write mapping output (e.g. SAM lines) to a separate file. If this option\n" \
  "     is not specified, mappings are written to standard output.\n"
  },
  {
  ochr:'O',
  otyp:OPTYP_FLAG,
  vnam:"",
  sdesc: "Preserve the order of the reads in the output (with '-n').",
  ldesc:\
  "     Output mappings in the order of the reads in the input files when using\n"\
  "     multiple threads (option '-n <nthreads>').\n\n"			\
  },
  {
  ochr:'p',
  otyp:OPTYP_FLAG,
  vnam:"",
  sdesc: "Report split alignments.",
  ldesc:\
  "     Report partial alignments if they are complementary on the read (split\n"\
  "     reads).\n"
  },
  {
  ochr:'q',
  otyp:OPTYP_INT,
  vnam:"minbasq",
  sdesc: "Base quality threshold <= 10 (default 0).",
  ldesc:\
  "     Sets a base quality threshold (0 <= minbasq <= 10, default 0).\n"\
  "     K-mer words of the read with nucleotides that have a base quality below\n" \
  "     this threshold are not looked up in the hash index.\n"
  },
  {
  ochr:'r',
  otyp:OPTYP_INT,
  vnam:"seed",
  sdesc: "Random assignment of degen. mappings (mark 'unmapped' if < 0).",
  ldesc:\
  "     If <seed> >= 0 report an alignment selected at random where there are\n"\
  "     multiple mappings with the same best alignment score. With <seed> = 0\n" \
  "     (default) a seed is derived from the current calendar time. If <seed>\n" \
  "     < 0 reads with multiple best mappings are reported as 'not mapped'.\n"
  },
  {ochr:'S',
   otyp:OPTYP_STRING,
   vnam:"scorspec",
   sdesc: "Set alignment penalties,\n"\
   "           e.g 'match=1,mismatch=-2,gapopen=-4,gapext=-3' (default).",
   ldesc:\
   "     Specify alignment penalty scores for a match or mismatch (substitution),\n"\
   "     or for opening or extending a gap. <scorspec> is a comma speparated\n"\
   "     list of integer assigments to one or more of the following variables:\n"\
   "     match, subst, gapopen, gapext, i.e. 'gapopen=-5,gapext=-4' (no spaces\n"\
   "     allowed in <scorespec>). Default:'match=1,subst=-2,gapopen=-4,gapext=-3'\n"
  },
#ifdef HAVE_BAMBAMC
  {
  ochr:'T',
  otyp:OPTYP_STRING,
  vnam:"tmpdir",
  sdesc: "Write temporary files do specified directory.",
  ldesc:\
  "     Write temporary files to directory <tmpdir> (used with input files in\n" \
  "     SAM/BAM format).\n"
  },
#endif   
  {
  ochr:'w',
  otyp:OPTYP_FLAG,
  vnam:"",
  sdesc: "Use complexity weighted Smith-Waterman scores.",
  ldesc:\
  "     Smith-Waterman scores are complexity weighted.\n"
  },
  {
  ochr:'x',
  otyp:OPTYP_FLAG,
  vnam:"",
  sdesc: "Exhaustive search for alignments (at the cost of speed).",
  ldesc:\
  "     This flag triggers a more exhaustive search for alignments at the cost\n"\
  "     of speed. In paired-end mode each mate is mapped independently.(By\n"\
  "     default the mate with fewer hits in the hash index is mapped first and\n"\
  "     the vicinity is searched for mappings of its mate.)\n"
  },
  {
  ochr:'y',
  otyp:OPTYP_FLT,
  vnam:"minid",
  sdesc: "Identity threshold (default: 0).",
  ldesc:\
  "     Sets an identity threshold for a mapping to be reported (default: 0).\n" \
  "     <minid> specifies the number of exactly matching nucleotides either as\n" \
  "     a positive integer or as a fraction of the read length (<= 1.0).\n"
  },
  {0,0,0,0,0}
};

//"  -c <maxhit>\n"
//"     K-mer words with more than <maxhit> hits across the genome are\n"
//"     disregarded (default <maxhit> = 10000).\n\n"

static const TASKDOC MENU_TASKDOC_MAP ={
 synopsis:
 "  smalt map [OPTIONS] <index_name> <query_file> [<mate_file>]",
 description:
#ifdef HAVE_BAMBAMC
 "  Map query reads onto the reference sequences. The reads are provided in\n" \
 "  FASTA/FASTQ format or in SAM/BAM format in the file <query_file>. If the\n"	\
 "  name of a second file <mate_file> is specified, both files are in\n" \
 "  FASTA/FASTQ format and reads are mapped in pairs. If <query_file> is in\n" \
 "  SAM/BAM format, single reads and paired reads can be mixed.\n\n"	\
 "  The reference sequences and k-mer word index are read from the binary\n" \
 "  files <index_name>.sma and <index_name>.smi which must have been created\n"	\
 "  by the 'index' task (type 'smalt index -H' for help).\n",
#else
 "  Map query reads onto the reference sequences. The reads are provided in\n"\
 "  FASTA/FASTQ format in the file <query_file>. If the name of a second file\n" \
 "  <mate_file> is specified, the reads in both files are mapped in pairs.\n" \
 "  The reference sequences and k-mer word index are read from the binary\n" \
 "  files <index_name>.sma and <index_name>.smi which must have been created\n"	\
 "  by the 'index' task (type 'smalt index -H' for help).\n",
#endif
 optdoc: MENU_OPTDOC_MAP,
};

/*****************************************************************************
 ******************************** CHECK TASK *********************************
 *****************************************************************************/

static const TASKDOC MENU_TASKDOC_CHECK = {	
 synopsis:
 "  smalt check <query_file> [<mate_file>]",
 description:
 "  Check FASTA/FASTQ read files. If <mate_file> is specified, the reads are\n" \
 "  in pairs.\n",
 optdoc: NULL
};

/*****************************************************************************
 ******************************* SAMPLE TASK *********************************
 *****************************************************************************/

static const OPTDOC MENU_OPTDOC_SAMPLE[] = {\
  {
  ochr:'H',
  otyp:OPTYP_FLAG,
  vnam:"",
  sdesc: "Print more extensive help on options.",
  ldesc:\
  "     Print these instructions.\n"
  },
  {
  ochr:'F',
  otyp:OPTYP_STRING,
  vnam:"inform",
#ifdef HAVE_BAMBAMC
  sdesc: "Input format [fastq (default)|sam|bam].",
  ldesc:\
  "     Specifies the input format. <inform> can be either 'fastq' (default),\n" \
  "     'sam' or 'bam' (see: samtools.sourceforge.net).\n"
#else
  sdesc: "Input format. fastq (default) is the only vailable format.",
  ldesc:\
  "     Specifies the input format. The only available format is fastq (default).\n"\
  "     Support for BAM and SAM formats (see: samtools.sourceforge.net) depends\n"\
  "     on additional libraries (not installed).\n" 
#endif 
  },
  {
  ochr:'m',
  otyp:OPTYP_INT,
  vnam:"minscor",
  sdesc: "Threshold of the alingment score.",
  ldesc:\
  "     Sets an absolute threshold of the Smith-Waterman scores. Mappings with\n" \
  "     scores below that threshold will not be reported. The default is\n" \
  "     <minscor> = <wordlen> + <stepsiz> - 1.\n" 
  },
  {
  ochr:'n',
  otyp:OPTYP_INT,
  vnam:"nthreads",
  sdesc: "Run multi-threaded with this number of threads.",
  ldesc:\
  "    Run in multi-threaded mode. <nthread> is the number of threads forked.\n"
  },
  {
  ochr:'o',
  otyp:OPTYP_STRING,
  vnam:"oufilnam",
  sdesc: "Write output to specified file (default: stdout).",
  ldesc:\
  "     Write mapping output (e.g. SAM lines) to a separate file. If this option\n" \
  "     is not specified, mappings are written to standard output.\n"
  },
  {
  ochr:'q',
  otyp:OPTYP_INT,
  vnam:"minbasq",
  sdesc: "Base quality threshold <= 10 (default 0).",
  ldesc:\
  "     Sets a base quality threshold (0 <= minbasq <= 10, default 0).\n"\
  "     K-mer words of the read with nucleotides that have a base quality below\n" \
  "     this threshold are not looked up in the hash index.\n"
  },
  #ifdef HAVE_BAMBAMC
  {
  ochr:'T',
  otyp:OPTYP_STRING,
  vnam:"tmpdir",
  sdesc: "Write temporary files to specified directory.",
  ldesc:\
  "     Write temporary files to directory <tmpdir> (used with input files in\n" \
  "     SAM/BAM format).\n"
  },
#endif 
  {
  ochr:'u',
  otyp:OPTYP_INT,
  vnam:"nreads",
  sdesc: "Map only every <nreads>-th read pair (default 100).",
  ldesc:\
  "     Map only every <nreads>-th read pair (default 100).\n"
  },
  {0, 0, 0, 0, 0}
};

static const TASKDOC MENU_TASKDOC_SAMPLE = {
 synopsis:
 "  smalt sample [OPTIONS] <index_name> <query_file> [<mate_file>]",
 description:
  "  Sample insert size distribution for paired reads. A subset of the read\n" \
  "  pairs is aligned with a reference in order to derrive the distribution of\n" \
  "  insert sizes. The reference sequences and index are read from the files\n"	\
  "  <index_name>.sma and <index_name>.smi created by the 'index' task (type\n"	\
 "  'smalt index -H' for help).\n",
 optdoc: MENU_OPTDOC_SAMPLE
};

enum MENU_OPTION_DEFAULTS {
  MENU_DEFAULTS_DISKLEVEL = 0,  /**< read entire hash table into memory by default */
  MENU_KMERLEN_MAX = 20,        /**< Maximal allowed k-mer length */
  MENU_KMERLEN_MIN = 3,         /**< Minimal allowed k-mer length */
  MENU_DEFAULTS_KMER = 13,      /**< k-tuple length */
  MENU_DEFAULTS_SKIP = 6,       /**< step size */
  MENU_DEFAULTS_MINSCOR = 18,   /**< minimum Smith-Waterman score */
  MENU_DEFAULTS_SCOREDIFF = 0,  /**< best mappings only */
  MENU_DEFAULTS_DEPTH = 500,    /**< target for the maximum number of reported hits and threshold
				 * for potentially matching regions classified as repeats */
  MENU_DEFAULTS_NCUT = 10000,    /**< Threshold in the number of hits a
				  * single k-tuple can have in the hashed
				  * (subject) sequence. K-tuples with hit
				  * frequencies higher than that are
				  * considered non-informative and are
				  * disregarded.*/
  MENU_DEFAULTS_MAXHIT = 100000, /**< Default upper limit for the number of hits 
				  * in a hit list */
  MENU_DEFAULTS_NTHREAD = 0,     /**< Default number of threads forked from main thread */
  MENU_DEFAULTS_MININSERT = 0,   /**< default minimum insert size */
  MENU_DEFAULTS_MAXINSERT = 500, /**< default maximum insert size */
  MENU_SCOREDIFF_ALL = -1,       /**< negative means find all possible alignments */
  MENU_DEFAULTS_RANDSEED = 0,   /**< 0: signals random picking of repeat alignments with 
				 * calendar derrived seed */
  MENU_DEFAULTS_MINBASQUAL = 0,  /**< base quality threshold for k-mer words */
  MENU_DEFAULTS_READSKIP = 100,  /**< map only every MENU_DEFAULTS_READSKIP when sampling insert
				  * sizes. */
  MENU_DEFAULTS_MINMAPQ = 20,    /**< Only use read pairs for the sampling of insert sizes, where 
				  * both mates have a mapping quality score of at least 
				  * MENU_DEFAULTS_MINMAPQ. */
};

static const double MENU_DEFAULTS_MINCOVER = 0.0;
static const double MENU_DEFAULTS_MINIDENTITY = 0.0;

enum OPTION_INTLIMITS {
  OPTLIMIT_MAXDISKLEVEL = 2,
  OPTLIMIT_MINBASQ = 10,          /**< maximum base quality threshold */
};

enum OPTERR {
  OPTERR_SUCCESS  =  0,
  OPTERR_NOARG    = -1, /* option requires an argument */
  OPTERR_NOMEM    = -2, /* could not allocate memory */
  OPTERR_UNKNOWN  = -3, /* unknown option */
  OPTERR_PAIRARG  = -4  /* option requires a pair of integers separated by colon */
};

typedef unsigned char UCHAR;
typedef uint8_t BOOL_t;
typedef uint16_t MENUFLG;
typedef uint8_t OUFMTYP_t;
typedef uint8_t OUFMTSAMFLG_t; /**< holds a combination of OUFMT_SAM_FLAGS */

static const char OUFMT_TYPSEP[] = ":"; /**< Separates output format keywords */
static const char OUFMT_PARSEP[] = ","; /**< Separates output format parameters */
static const char OUFMT_LISTSEP[] = ","; /**< Separates a list of variable-value assignments */

enum {
  OUFMT_ASSIGNMENT_CHAR = '='
};

typedef struct VARLST_ {
  const char *varnam;
  uint8_t typ;
  int8_t val;
  int8_t max;
  int8_t min;
} VARLST;

enum MENU_PENALTY_TYPES {
  MENUPNLTYP_MATCH = 0,
  MENUPNLTYP_MISMATCH = 1,
  MENUPNLTYP_GAPINIT = 2,
  MENUPNLTYP_GAPEXT = 3,
  MENUPNLTYP_NUM = 4,
};

static const VARLST MENU_PENALTY_LST[] = {
  {
  varnam: "match",
  typ: MENUPNLTYP_MATCH,
  val: 1,
  max: 127,
  min: 0,
  },
  {
  varnam: "subst",
  typ: MENUPNLTYP_MISMATCH,
  val: -2,
  max: 0,
  min:-127,
  },
  {
  varnam: "gapopen",
  typ: MENUPNLTYP_GAPINIT,
  val: -4,
  max: 0,
  min:-127
  },
  {
  varnam: "gapext",
  typ: MENUPNLTYP_GAPEXT,
  val: -3,
  max: 0,
  min:-127
  },
  {0,0,0,0,0}
};

typedef struct _OUFMT {
  OUFMTYP_t typ;       /**< one of MENU_OUTPUT_FORMATS */
  MENUOFMTFLG_t flags; /**< combination of output format flags (e.g. MENU_OUFMTSAM_FLAGS) */
} OUFMT;

typedef struct _INDEXMENU {
  unsigned char kmer;
  unsigned char skip;
} INDEXMENU;

typedef struct _MAPMENU {
  OUFMT oufmt;              /**< Ouptut format including parameters */
  UCHAR inform;             /**< Input format (one of MENU_INPUT_FORMATS) */
  short nthread;            /**< Number of threads used minus 1 (main thread) */
  int8_t penalties[MENUPNLTYP_NUM]; /**< alignment penalty scores */
  uint8_t penaltyflags;       /**< Least significant 4 bits indicate which penalty scores
			      * were set on the command line (according to integer values
			      * of MENU_PENALTY_TYPES, highes value->most significant bit. */
  int ncut;                 /**< Cut-off in the number of hits of a k-tuple in the 
			     * hashed (subject) sequence. K-tuples with more than
			     * ncut hits are disregarded. They are not informative. */
  int maxhit;               /**< Maximum total number of hits in hit list */
  double mincover;          /**< Threshold in the minimum coverage of a read by k-tuples
			     * to qualify as a potentially matching segment */
  short target_depth;       /**< Aim to process about \<maxdepth\> potentially matching regions */
  short minscore;           /**< minimum Smith-Waterman score */
  int randseed_repeat;      /**< Seed for random number generator, used for the selection of 
			     * repetitive alignments. */
  short scorediff;          /**< Only report matches with Smith-Waterman scores 
			     * >= (maxscore - scorediff). If set to MENU_SCOREDIFF_ALL 
			     * find as many mappings as you possibly can */
  UCHAR minbasq;            /**< base quality threshold */
  double minidentity;       /**< identity threshold for reported reads */
  char diskuse;             /**< specifies which parts of the hash table are kept on disk:
			     * 0 -> entire hash table in memory (default), 
			     * 1 -> hash table body on disk,
			     * 2 -> hash table body and genomic sequences on disk. */
  int insert_range[MENU_PAIR_RANGENUM]; /**< lower and upper limit of the insert size for paired-end reads */
  char *oufilnam;           /**< Output file */
  char *insfilnam;          /**< File with insert size distribution */
  char *tmpdirnam;          /**< Directory for temporary output files 
			     * (used with SAM/BAM paired read input files) */
  int readskip;            /**< sample only every nskip_read read pair */
  int mapqmin;             /**< Treshold in mapping quality */
  UCHAR pairtyp;           /**< orientation of mates according to library (one of MENU_READPAIR_TYPES) */
  MENUFLG_t flags;         /**< Combination of  MENU_FLAGS */
} MAPMENU;

typedef struct _optflags { /**> for keeping track of specified options */
  unsigned int kmer: 1;
  unsigned int skip: 1;
  unsigned int ncut: 1;
  unsigned int mincover: 1;
  unsigned int minscore: 1;
  unsigned int maxhit: 1;
  unsigned int scorediff: 1;
  unsigned int oufilnam: 1;
  unsigned int insfilnam: 1;
  unsigned int insertmax: 1;
  unsigned int insertmin: 1;
  unsigned int randrepeat: 1;
  unsigned int minbasq: 1;
  unsigned int minidentity: 1;
  unsigned int mapqmin: 1;
  unsigned int readskip: 1;
  unsigned int pairtyp: 1;
  unsigned int penalties: 1;
} OPTFLAGS;

struct _MenuOpt {
  int argc;               /**< Number of command line arguments */
  char **argv;             /**< Pointer to command line arguments */
  char  subprog;           /**< one of MENU_SUBPROG */
  void *paramp;
  char *indir;            /** directory for input files */
  int ninfil;              /** number of input files */
  char **filnams;          /** Beginning of the array of file names */
  EString estrbuf;         /**< Buffer for string operations */
};


typedef int (MENU_PARSER)(MenuOpt *, OPTFLAGS *, int, char **);
typedef int (MENU_CHECKER) (MenuOpt *, const OPTFLAGS *);

/******************************************************************************
 ******************************* Pivate Methods *******************************
 ******************************************************************************/
static void fprintTaskDoc(FILE *oufp, const TASKDOC *tdocp, BOOL_t isLong)
{
  if (NULL == tdocp)
    return;
  if (NULL != tdocp->synopsis)
    fprintf(oufp,"\nSYNOPSIS:\n%s\n", tdocp->synopsis);
  if ((isLong) && NULL != tdocp->description)
    fprintf(oufp,"\nDESCRIPTION:\n%s\n", tdocp->description);
  if (NULL != tdocp->optdoc) {
    const OPTDOC *odp;
    fprintf(oufp, "\nOPTIONS:\n");
    for (odp = tdocp->optdoc; odp->ochr != 0; odp++) {
      fprintf(oufp, "  -%c", odp->ochr);
      if ((isLong)) {
	if (odp->otyp != OPTYP_FLAG)
	  fprintf(oufp, " <%s [%s]>\n", odp->vnam, OPTION_TYPSTR[odp->otyp]);
	fprintf(oufp, "%s\n", odp->ldesc);
      } else {
	if (odp->otyp != OPTYP_FLAG)
	  fprintf(oufp, " [%s] %s\n", OPTION_TYPSTR[odp->otyp], odp->sdesc);
	else
	  fprintf(oufp, "       %s\n", odp->sdesc);
      }
    }
  }
  if (!(isLong))
    fprintf(oufp, "\n");
  return;
}

static void printBlurb(FILE *oufp)
{
  /* size_t l; */
  /*   size_t linlen = strlen(MENU_PROGNAM) + \ */
  /*     strlen(MENU_RELEASE_VERSION) + 13; */

  /*   fputc('\n', oufp); */
  /*   for (l=0; l<linlen; l++) fputc('=', oufp); */
  /*   fprintf(oufp, "\n= %s version %s =\n",MENU_PROGNAM , MENU_RELEASE_VERSION); */
  /*   for (l=0; l<linlen; l++) fputc('=', oufp); */
  /*   fputc('\n', oufp); */
  fprintf(oufp, "\n%s\n", MENU_PROGNAM); 
  fprintf(oufp, "Version: %s\n", MENU_RELEASE_VERSION);
  fprintf(oufp, "Date:    %s\n", MENU_RELEASE_DATE);
  //fprintf(oufp, "\nRelease date: %s\n", MENU_RELEASE_DATE);
  fprintf(oufp, "Author:  %s (%s)\n\n", MENU_RELEASE_AUTHORS, MENU_RELEASE_BUGREPORT);
  //fprintf(oufp, "%s\nAll Rights Reserved.\n\n", MENU_COPYRIGHT_NOTICE);
  fprintf(oufp, "%s\n\n", MENU_COPYRIGHT_NOTICE);
}

static void exitOptionError(const char *option, const char *expl)
     /**< Print error message on screen and exit */
{
  printf("Command line error: ");
  if (option) 
    printf("option '%s' ", option);
  if (expl)
    printf("%s", expl);
  printf("\n");
  exit(EXIT_FAILURE);
}

static void exitPairedReadError(const char *option)
{
  printf("Command line error: ");
  if (option)
    printf("option '%s' ", option);
  printf("implies paired read mapping for which 2 read files are expected as input.\n");
  exit(EXIT_FAILURE);
}

static void exitArgumentRangeError(const char *option, long minval, long maxval)
{
  printf("Command line error: ");
  if (option) {
    printf("option '%s' requires as argument integer values ", option);
    if (maxval > minval)
      printf("between %li and %li", minval, maxval);
    else
      printf(">= %li", minval);
  }

  exit(EXIT_FAILURE);
}

static void parseListOfKeyValueAssignments(int8_t vals[], uint8_t *bitfld, 
					   EString *estrbufp,
					   const VARLST *lst, short nelem, 
					   const char *listp)
{
  char *parkey, *liststr;
  short i;

  *bitfld = 0;
  ESTRING_BLANK(*estrbufp);
  if (ESTRING_APPEND(*estrbufp, listp))
    exit(EXIT_FAILURE);
  liststr = ESTRING_GETSTR(*estrbufp);

  parkey = strtok(liststr,OUFMT_LISTSEP);
  
  while (parkey != NULL) {
    char *value, *cp;
    int ival = 0;
    for (value = parkey; (*value); value++) 
      if (*value == OUFMT_ASSIGNMENT_CHAR) {
	*value++ = '\0';
	break;
      }
    /* check it is all digits */
    cp = value;
    if (*cp == '+' || *cp == '-')
      cp++;
    for (; (*cp) && isdigit((int) *cp); cp++);
    if ( '\0' != *cp || '\0' == *value) {
      printf("Error: wrong format string %s.\n", 
	     parkey);
      exit(EXIT_FAILURE);
    }
    ival = atoi(value);
    if (ival > UINT8_MAX) {
      printf("Error: wrong format string %s.\n", 
	     parkey);
      exit(EXIT_FAILURE);
    }
    for (i=0; i<nelem; i++) {
      if (!strcmp(parkey, lst[i].varnam)) {
	if (ival < lst[i].min || ival > lst[i].max) {
	  printf("Error: value %s=%hi is outside the range [%hi,%hi].\n",
		 parkey, ival, (short) lst[i].min, (short) lst[i].max);
	  exit(EXIT_FAILURE);
	}
	vals[i] = (uint8_t) ival;
	*bitfld = (uint8_t)((*bitfld) | (((uint8_t) 1)<<i));
	*(--value) = OUFMT_ASSIGNMENT_CHAR;
	break;
      }
    }
    if (i == nelem) {
      printf("Error: Keyword '%s' not recognised in list of penalty assignments.\n", 
	     parkey);
      exit(EXIT_FAILURE);
    }
    parkey = strtok(NULL, OUFMT_PARSEP); 
  }
}

static void parseSamParams(MENUOFMTFLG_t *samflg, 
			   const char *formatp, 
			   const char *optstr)
{
  BOOL_t is_err = 0;
  char *parkey, *fmtstr;
  EString *estrbufp;

  if (NULL == ESTRING_NEW(estrbufp))
    exit(EXIT_FAILURE);

  *samflg = MENU_OUFMTSAM_HEADER;

  if (formatp == NULL)
    return;

  ESTRING_BLANK(*estrbufp);
  if (ESTRING_APPEND(*estrbufp, formatp)) {
    ESTRING_DELETE(estrbufp);
    exit(EXIT_FAILURE);
  }
  fmtstr = ESTRING_GETSTR(*estrbufp);

  parkey = strtok(fmtstr, OUFMT_PARSEP);

  while (parkey != NULL && !(is_err)) {
    if (!strcmp(parkey, "nohead")) {
      *samflg &= ~ MENU_OUFMTSAM_HEADER;
    } else if (!strcmp(parkey, "clip")) {
      *samflg |= MENU_OUFMTSAM_CLIPPED;
    } else if (!strcmp(parkey, "x") || !strcmp(parkey, "X")) {
      *samflg |= MENU_OUFMTSAM_XCIGAR;
    } else {
      is_err = 1;
      break;
    }
    parkey = strtok(NULL, OUFMT_PARSEP);
  }

  ESTRING_DELETE(estrbufp);

  if (is_err) {
    printf("Error: Keyword '%s' not recognised in format specification (%s sam).\n", 
	   parkey, optstr);
    exit(EXIT_FAILURE);
  }
}

static void parseOutputFormat(OUFMT *oufmt, EString *estrbufp, const char *optstr, const char *formatp)
{
  BOOL_t is_err = 0;

  ESTRING_BLANK(*estrbufp);
  if (ESTRING_APPEND(*estrbufp, formatp))
    exit(EXIT_FAILURE);

  if (NULL == formatp) {
    is_err = 1;
  } else {
    char *fmtstr = ESTRING_GETSTR(*estrbufp);
    char *fmtkey = strtok(fmtstr, OUFMT_TYPSEP);

    if (!strcmp(fmtkey, "cigar")) {
      oufmt->typ = MENU_OUTFORM_CIGAR;
    } else if (!strcmp(fmtkey, "sam")) {
      char *par = strtok(NULL, OUFMT_TYPSEP);
      oufmt->typ = MENU_OUTFORM_SAM;
      parseSamParams(&oufmt->flags, par, optstr);
    } else if (!strcmp(fmtkey, "samsoft")) {
      /* for backward compatibility */
      oufmt->typ = MENU_OUTFORM_SAM;
      oufmt->flags = MENU_OUFMTSAM_HEADER;
    } else if (!strcmp(fmtkey, "bam")) {
      char *par = strtok(NULL, OUFMT_TYPSEP);
      oufmt->typ = MENU_OUTFORM_BAM;
      parseSamParams(&oufmt->flags, par, optstr);
#ifndef HAVE_BAMBAMC
      is_err = 2;
#endif
    } else if (!strcmp(fmtkey, "ssaha")) {
      oufmt->typ = MENU_OUTFORM_SSAHA;
    } else if (!strcmp(fmtkey, "gff")) {
      oufmt->typ = MENU_OUTFORM_GFF;      
    } else {
      is_err = 1;
    }
  }
  if ((is_err)) {
#ifdef HAVE_BAMBAMC
    exitOptionError(optstr, "requires an output format [sam|bam|cigar|ssaha|gff] as argument.");
#else
    if (2 == is_err) 
      exitOptionError(optstr, "BAM output format not available. Bambam library not installed.");
    else
      exitOptionError(optstr, "requires an output format[sam|cigar|ssaha|gff] as argument.");
#endif
  }
}

static short parseOption(void *optarg, const char *optnam, 
			 char optyp, int narg, char * const *argp)
     /**< Compare the next token in argp to optnam. Return 0 if
      * argp[0] and optnam do not match.  If they match, return the
      * following token with allocated memory in optarg.  Return the
      * number of tokens read from argp.  Return values < 0 indicate
      * an error: -1: no argument found -2: could not allocate memory
      * for string -3: unknown opttyp
      */
{
  short rv;
  char *cp, **strpp;

  if (strcmp(argp[0], optnam)) return OPTERR_UNKNOWN;
  switch (optyp) {
  case OPTYP_STRING: /* assume optarg is pointer to a string */
    if (narg < 2) rv = OPTERR_NOARG;
    else {
      strpp = (char **) optarg;
      free(*strpp);
      ESTRCPY(*strpp, argp[1]);
      rv = (short)((*strpp)? 2: OPTERR_NOMEM);
    }
    break;
  case OPTYP_INT:
    if (narg < 2) rv = OPTERR_NOARG;
    else {
      *(int *) optarg = atoi(argp[1]);
      rv = 2;
    }
    break;
  case OPTYP_INT_PAIR: /* pair of ints separated by colon */
    if (narg < 2) rv = OPTERR_NOARG;
    else {
      /* find colon and replace by 0 */
      for (cp=argp[1]; (*cp); cp++)
	if (*cp == MENU_PAIR_SEPARATOR) {
	  *cp++ = '\0';
	  break;
	}
      if (*cp) {
	*(int *) optarg = atoi(argp[1]);
	*(((int *) optarg)+1) = atoi(cp);
	rv = 2;
      } else {
	rv = OPTERR_PAIRARG;
      }
    }
    break;
  case OPTYP_FLT:
    if (narg < 2) rv = OPTERR_NOARG;
    else {
      *(double *) optarg = atof(argp[1]);
      rv = 2;
    }
    break;
  case OPTYP_FLAG:
    rv = 1;
    break;
  default:
    rv = OPTERR_UNKNOWN;
    break;
  }
    
  return rv;
}

/******************************************************************************
 ********************** Pivate Methods of Type MenuOpt ************************
 ******************************************************************************/
static int initHashOptions(MenuOpt *menup)
{
  INDEXMENU *ip;
  free(menup->paramp);
  EMALLOCP0(ip);
  if (!ip)
    return ERRCODE_NOMEM;
  menup->subprog = MENU_INDEX;
  menup->paramp = ip;
  ip = (INDEXMENU *) menup->paramp;
  ip->kmer = MENU_DEFAULTS_KMER;
  ip->skip = MENU_DEFAULTS_SKIP;
  
  return ERRCODE_SUCCESS;
}

static int initMapOptions(MenuOpt *menup)
{
  MAPMENU *mp;
  int i;
  free(menup->paramp);
  EMALLOCP0(mp);
  if (!mp)
    return ERRCODE_NOMEM;
  if (menup->subprog != MENU_MAP
      && menup->subprog != MENU_SAMPLE)
    return ERRCODE_ASSERT;
  menup->paramp = mp;
  mp = (MAPMENU *) menup->paramp;
  mp->oufmt.typ = MENU_OUTFORM_SAM;
  mp->oufmt.flags = MENU_OUFMTSAM_HEADER;
  mp->inform = MENU_INFORM_FASTQ;
  mp->nthread = MENU_DEFAULTS_NTHREAD;
  for (i=0; i< MENUPNLTYP_NUM; i++)
    mp->penalties[i] = MENU_PENALTY_LST[i].val;
  mp->ncut = MENU_DEFAULTS_NCUT;
  mp->maxhit = MENU_DEFAULTS_MAXHIT;
  mp->mincover = MENU_DEFAULTS_MINCOVER;
  mp->minscore = MENU_DEFAULTS_MINSCOR;
  mp->target_depth = MENU_DEFAULTS_DEPTH;
  mp->scorediff = MENU_DEFAULTS_SCOREDIFF;
  mp->flags = MENUFLAG_VERBOSE | MENUFLAG_RANDREPEAT;
  mp->diskuse = (char) 0;
  mp->oufilnam = NULL;
  mp->insfilnam = NULL;
  mp->tmpdirnam = NULL;
  mp->insert_range[0] = MENU_DEFAULTS_MININSERT;
  mp->insert_range[1] = MENU_DEFAULTS_MAXINSERT;
  mp->randseed_repeat = MENU_DEFAULTS_RANDSEED;
  mp->minbasq = MENU_DEFAULTS_MINBASQUAL;
  mp->minidentity = MENU_DEFAULTS_MINIDENTITY;
  mp->readskip = 0;
  mp->mapqmin = 0;
  mp->pairtyp = MENU_READPAIRTYP_UNKNOWN;
  return ERRCODE_SUCCESS;
}


/* Functions of type MENU_CHECKER */
static int checkHashDefaults(MenuOpt *menup, const OPTFLAGS *optflp)
{
  INDEXMENU *ip;

  if (menup->subprog != MENU_INDEX)
    return ERRCODE_FAILURE;

  ip = (INDEXMENU *) menup->paramp;

  /* if -skip is not specified but -kmer is,
     set menup->skip to menup->kmer (i.e. tiling k-tuples by default) */
  if (!optflp->skip) ip->skip = ip->kmer;
  if(menup->ninfil != 2)
    exitOptionError(NULL, "Expected name of an input file containing the genomic reference sequences"\
		     "in FASTA/FASTQ format.");
  return ERRCODE_SUCCESS;
}

static int checkMapDefaults(MenuOpt *menup, const OPTFLAGS *optflp)
{
  MAPMENU *mp = (MAPMENU *) menup->paramp;
  if (menup->subprog != MENU_MAP)
    return ERRCODE_FAILURE;

  if (menup->ninfil < 2 || menup->ninfil > 3) 
    exitOptionError(NULL, "Expected one or two file names.");

  if (mp->inform == MENU_INFORM_FASTQ) {
    if (menup->ninfil == 3) {
      mp->flags |= MENUFLAG_PAIRED;
    } else {
      mp->pairtyp = MENU_READPAIRTYP_SINGLE;
    if ((optflp->insertmax))
      exitPairedReadError("-i");
    if ((optflp->insertmin)) 
      exitPairedReadError("-j");
    if ((optflp->pairtyp))
      exitPairedReadError("-l");
    }
  }
  if (mp->insert_range[0] > mp->insert_range[1]) {
    int tmp = mp->insert_range[0];
    mp->insert_range[0] = mp->insert_range[1];
    mp->insert_range[1] = tmp;
  }
  if (!optflp->scorediff) mp->scorediff = 0; /* by default only best mappings in paired-end mode */
  if (!optflp->pairtyp) mp->pairtyp = MENU_READPAIRTYP_PAIREDEND;

/*     if (mp->flags&MENUFLAG_EXHAUSTIVE) */
/*       exitPairedReadError("-x"); */
	     
  if (mp->maxhit < 2*mp->ncut) mp->maxhit = 2*mp->ncut;
  if (mp->scorediff != 0) {
    mp->flags &= ~MENUFLAG_RANDREPEAT;
  }

  if ((optflp->mincover) && !(mp->flags & MENUFLAG_EXHAUSTIVE)) {
    exitOptionError("-c", "can only be used in combination with the '-x' flag.");
  }

  return ERRCODE_SUCCESS;
}

static int checkSampleDefaults(MenuOpt *menup, const OPTFLAGS *optflp)
{
  MAPMENU *mp = (MAPMENU *) menup->paramp;
  if (menup->subprog != MENU_SAMPLE)
    return ERRCODE_FAILURE;

  if (mp->inform == MENU_INFORM_FASTQ) {
    if (menup->ninfil != 3) 
      exitOptionError(NULL, "Expected an index name and 2 file names.");
    mp->flags |= MENUFLAG_PAIRED;
  }

  mp->flags |= MENUFLAG_EXHAUSTIVE;
  mp->scorediff = 0; /* by default only best mappings in paired-end mode */
 
  if (mp->maxhit < 2*mp->ncut) mp->maxhit = 2*mp->ncut;
  mp->flags &= ~MENUFLAG_RANDREPEAT;

  if (!optflp->readskip)
    mp->readskip = MENU_DEFAULTS_READSKIP;
  if (!optflp->mapqmin)
    mp->mapqmin = MENU_DEFAULTS_MINMAPQ;
    
  return ERRCODE_SUCCESS;
}

/* Functions of type MENU_PARSER */
static int parseHashOptions(MenuOpt *menup, OPTFLAGS *optflp, int narg, char **argp)
     /**< Parse options concerning hash table build 
      * Return the number of tokens processed. 0 on error */
{
  int n, iv;
  INDEXMENU *ip = (INDEXMENU *) menup->paramp;
  if (menup->subprog != MENU_INDEX)
    return OPTERR_UNKNOWN;

  switch (argp[0][1]) {
  case 'D': /* input directory */
    n = parseOption(&(menup->indir), "-D", OPTYP_STRING, narg, argp);
    if (n == OPTERR_NOARG) exitOptionError(argp[0], "requires a directory name as argument.");
    break;
  case 'k': /* kmer word length */
    n = parseOption(&iv, "-k", OPTYP_INT, narg, argp);
    if (n != OPTERR_UNKNOWN) {
      if (n == OPTERR_NOARG  || iv<MENU_KMERLEN_MIN || iv > MENU_KMERLEN_MAX) 
	exitArgumentRangeError(argp[0], MENU_KMERLEN_MIN, MENU_KMERLEN_MAX);
      ip->kmer = (char) iv;
      optflp->kmer = 1;
    }
    break;
  case 's': /* skip step size */
    n = parseOption(&iv, "-s", OPTYP_INT, narg, argp);
    if (n != OPTERR_UNKNOWN) {
      if (n == OPTERR_NOARG || iv<1 || iv > 127)
	exitArgumentRangeError(argp[0], 1, 127);
      ip->skip = (char) iv;
      optflp->skip = 1;
    }
    break;
    
  default:
    n=OPTERR_UNKNOWN;
    break;
  }

  return n;

}

static int parseMapOptions(MenuOpt *menup, OPTFLAGS *optflp, int narg, char **argp)
     /**< Parse options concerning the read mapping
      * Return the number of tokens processed. 0 on error */
{
  int n, iv, okflg=1, iserr = 0;
  double fv;
  MAPMENU *mp = (MAPMENU *) menup->paramp;
  if (menup->subprog != MENU_MAP)
    return OPTERR_UNKNOWN;

  switch (argp[0][1]) {
  case 'a':
    n = parseOption(NULL, "-a", OPTYP_FLAG, narg, argp);
    if (n != OPTERR_UNKNOWN) {
      mp->flags |= MENUFLAG_ALIGNMENT;
    }
    break;
  case 'D': /* directory for input files */
    n = parseOption(&(menup->indir), "-D", OPTYP_STRING, narg, argp);
    if (n == OPTERR_NOARG) exitOptionError(argp[0], "requires a directory name as argument.");
    break;
/*   case 'c': /\* maximum number of hits before a k-mer word is disregarded *\/ */
/*     n = parseOption(&iv, "-c", OPTYP_INT, narg, argp); */
/*     if (n != OPTERR_UNKNOWN) { */
/*       if (n == OPTERR_NOARG  || iv<1 || iv>SHRT_MAX)  */
/* 	exitArgumentRangeError(argp[0], 1, SHRT_MAX); */
/*       mp->ncut = iv; */
/*       optflp->ncut = 1; */
/*     } */
/*     break; */
  case 'c': /* minimum coverage can be either float or int */
    n = parseOption(&fv, "-c", OPTYP_FLT, narg, argp);
    if (n != OPTERR_UNKNOWN) {
      if (n == OPTERR_NOARG  || fv<-0.0 || fv>SHRT_MAX) 
	exitArgumentRangeError(argp[0], 0, SHRT_MAX);
      mp->mincover = fv;
      optflp->mincover = 1;
    }
    break;
  case 'd':
    n = parseOption(&iv, "-d", OPTYP_INT, narg, argp);
    if (n != OPTERR_UNKNOWN) {
      if (n == OPTERR_NOARG  || iv<SHRT_MIN || iv>SHRT_MAX) 
	exitArgumentRangeError(argp[0], SHRT_MIN, SHRT_MAX);
      mp->scorediff = (short) iv;
      optflp->scorediff = 1;
      mp->flags |= MENUFLAG_RELSCOR;
    }
    break;
  case 'f': /* output format */
    n = 2;
    iserr = 0;
    if (strcmp(argp[0], "-f")) {
      n = OPTERR_UNKNOWN;
    } else if (narg < 2) {
      n = OPTERR_NOARG;
      exitOptionError(argp[0], "requires an output format [cigar|sam|samsoft|ssaha|gff] as argument.");
    } else {
      parseOutputFormat(&mp->oufmt, &menup->estrbuf, argp[0], argp[1]);
    }
    break;
  case 'F': /* expected input format */
    iserr = 0;
    if (strcmp(argp[0], "-F")) {
      n = OPTERR_UNKNOWN;
    } else if (narg < 2) {
      n = OPTERR_NOARG;
      iserr = 1;
      exitOptionError(argp[0], "requires an input format [fastq|sam|bam] as argument.");
    } else {
      n = 2;
      if (!strcmp(argp[1], "fastq")) {
	mp->inform = MENU_INFORM_FASTQ;
      } else if (!strcmp(argp[1], "sam")) {
	mp->inform = MENU_INFORM_SAM;
#ifndef HAVE_BAMBAMC
	iserr = 2;
#endif
      } else if (!strcmp(argp[1], "bam")) {
	mp->inform = MENU_INFORM_BAM;
#ifndef HAVE_BAMBAMC
	iserr = 2;
#endif
      } else {
	iserr = 1;
      }
    }
    if ((iserr)) {
#ifdef HAVE_BAMBAMC
      exitOptionError(argp[0], "requires an input format [fastq|sam|bam] as argument.");
#else
      if (2 == iserr)
	exitOptionError(argp[0], "format not available. Requires bambam library (not installed).");
      else 
	exitOptionError(argp[0], "requires an input format as argument. The only available format is 'fastq'.");
#endif
    }
    break;
  case 'g':
    n = parseOption(&(mp->insfilnam), "-g", OPTYP_STRING, narg, argp);
    if (n == OPTERR_NOARG) exitOptionError(argp[0], "requires a file name as argument.");
    break;
  case 'i':
    n = parseOption(&iv, "-i", OPTYP_INT, narg, argp);
    if (n != OPTERR_UNKNOWN) {
      if (n == OPTERR_NOARG  || iv<0 || iv>INT_MAX) 
	exitArgumentRangeError(argp[0], 0, INT_MAX);
      mp->insert_range[1] = iv;
      optflp->insertmax = 1;
    }
    break;
  case 'j':
    n = parseOption(&iv, "-j", OPTYP_INT, narg, argp);
    if (n != OPTERR_UNKNOWN) {
      if (n == OPTERR_NOARG  || iv<0 || iv>INT_MAX) 
	exitArgumentRangeError(argp[0], 0, INT_MAX);
      mp->insert_range[0] = iv;
      optflp->insertmin = 1;
    }
    break;
  case 'l': /* library format */
    okflg = 1;
    if (strcmp(argp[0], "-l")) {
      n = OPTERR_UNKNOWN;
    } else if (narg < 2) {
      n = OPTERR_NOARG;
      okflg = 0;
    } else {
      n = 2;
      if (!strcmp(argp[1], "pe")) {
	mp->flags &= ~MENUFLAG_MATEPAIRLIB;
	mp->pairtyp = MENU_READPAIRTYP_PAIREDEND;
	optflp->pairtyp = 1;
      } else if (!strcmp(argp[1], "mp")) {
	mp->flags |= MENUFLAG_MATEPAIRLIB;
	mp->pairtyp = MENU_READPAIRTYP_MATEPAIR;
	optflp->pairtyp = 1;
      } else if (!strcmp(argp[1], "pp")) {
	mp->flags &= ~MENUFLAG_MATEPAIRLIB;
	mp->pairtyp = MENU_READPAIRTYP_SAMESTRAND;
	optflp->pairtyp = 1;
      } else 
	okflg = 0;
    }
    if (!okflg) 
      exitOptionError(argp[0], "requires a read pair library type [pe|mp|pp] as argument.");
    break;
  case 'm':
    n = parseOption(&iv, "-m", OPTYP_INT, narg, argp);
    if (n != OPTERR_UNKNOWN) {
      if (n == OPTERR_NOARG  || iv<0 || iv>SHRT_MAX) 
	exitArgumentRangeError(argp[0], 0, SHRT_MAX);
      mp->minscore = (short) iv;
      optflp->minscore = 1;
      mp->flags |= MENUFLAG_MINSCOR;
    }
    break;
  case 'n':
    n = parseOption(&iv, "-n", OPTYP_INT, narg, argp);
    if (n != OPTERR_UNKNOWN) {
      if (n == OPTERR_NOARG  || iv<1 || iv>SHRT_MAX) 
	exitArgumentRangeError(argp[0], 0, SHRT_MAX);
      mp->nthread = (short) iv;
    }
    break;     
  case 'o':
    n = parseOption(&(mp->oufilnam), "-o", OPTYP_STRING, narg, argp);
    if (n == OPTERR_NOARG) exitOptionError(argp[0], "requires a file name as argument.");
    break;
  case 'O':
    n = parseOption(NULL, "-O", OPTYP_FLAG, narg, argp);
    if (n != OPTERR_UNKNOWN) {
      mp->flags |= MENUFLAG_READORDER;
    }
    break;
  case 'p':
    n = parseOption(NULL, "-p", OPTYP_FLAG, narg, argp);
    if (n != OPTERR_UNKNOWN) {
      mp->flags |= MENUFLAG_SPLITREAD;
    }
    break;
  case 'q':
    n = parseOption(&iv, "-q", OPTYP_INT, narg, argp);
    if (n != OPTERR_UNKNOWN) {
      if (n == OPTERR_NOARG  || iv<0 || iv>OPTLIMIT_MINBASQ) 
	exitArgumentRangeError(argp[0], 0, OPTLIMIT_MINBASQ);
      mp->minbasq = (UCHAR) iv;
      optflp->minbasq = 1;
    }
    break; 
  case 'r':
    n = parseOption(&iv, "-r", OPTYP_INT, narg, argp);
    if (n != OPTERR_UNKNOWN) {
      if (n == OPTERR_NOARG  || iv>INT_MAX) 
	exitArgumentRangeError(argp[0], 0, INT_MAX);
      if (iv < 0) {
	mp->flags &= ~MENUFLAG_RANDREPEAT;
      } else {
	mp->randseed_repeat = iv;
	optflp->randrepeat = 1;
	mp->flags |= MENUFLAG_RANDREPEAT;
      }
    }
    break;
  case 'S': /* substitution (penalty scores) */
    n = 2;
    iserr = 0;
    if (strcmp(argp[0], "-S")) {
      n = OPTERR_UNKNOWN;
    } else if (narg < 2) {
      n = OPTERR_NOARG;
      exitOptionError(argp[0], "requires list of assignments as argument "\
		      "(eg. 'mismatch=-2,gapopen=-5,gapext=-4').");
    } else {
      parseListOfKeyValueAssignments(mp->penalties, &mp->penaltyflags,
				     &menup->estrbuf,
				     MENU_PENALTY_LST, MENUPNLTYP_NUM, 
				     argp[1]);
      optflp->penalties = 1;
    }
    break;
#ifdef HAVE_BAMBAMC
  case 'T': /* specify directory for temporary output files */
    n = parseOption(&(mp->tmpdirnam), "-T", OPTYP_STRING, narg, argp);
    if (n == OPTERR_NOARG) exitOptionError(argp[0], "requires a directory name as argument.");
    break;
#endif
  case 'w': /* flags complexity weighted Smith-Waterman scores */
    n = parseOption(NULL, "-w", OPTYP_FLAG, narg, argp);
    if (n != OPTERR_UNKNOWN) {
      mp->flags |= MENUFLAG_COMPLEXW;
    }
    break;
  case 'x': /* exhaustive search for read pairs */
    n = parseOption(NULL, "-x", OPTYP_FLAG, narg, argp);
    if (n != OPTERR_UNKNOWN) {
      mp->flags |= MENUFLAG_EXHAUSTIVE;
    }
    break;
  case 'y': /* minimum identity can be either float or int */
    n = parseOption(&fv, "-y", OPTYP_FLT, narg, argp);
    if (n != OPTERR_UNKNOWN) {
      if (n == OPTERR_NOARG  || fv<-0.0 || fv>SHRT_MAX) 
	exitArgumentRangeError(argp[0], 0, SHRT_MAX);
      mp->minidentity = fv;
      optflp->minidentity = 1;
    }
    break;
  default:
    n=OPTERR_UNKNOWN;
    break;
  }

  return n;
}

static int parseSampleOptions(MenuOpt *menup, OPTFLAGS *optflp, int narg, char **argp)
     /**< Parse options concerning the read mapping
      * Return the number of tokens processed. 0 on error */
{
  int n, iv;
  MAPMENU *mp = (MAPMENU *) menup->paramp;
  if (menup->subprog != MENU_SAMPLE)
    return OPTERR_UNKNOWN;

  switch (argp[0][1]) {
#ifdef HAVE_BAMBAMC
    int iserr;
  case 'F': /* expected input format */
    iserr = 0;
    if (strcmp(argp[0], "-F")) {
      n = OPTERR_UNKNOWN;
    } else if (narg < 2) {
      n = OPTERR_NOARG;
      iserr = 1;
      exitOptionError(argp[0], "requires an input format [fastq|sam|bam] as argument.");
    } else {
      n = 2;
      if (!strcmp(argp[1], "fastq"))
	mp->inform = MENU_INFORM_FASTQ;
      else if (!strcmp(argp[1], "sam")) 
	mp->inform = MENU_INFORM_SAM;
      else if (!strcmp(argp[1], "bam"))
	mp->inform = MENU_INFORM_BAM;
      else 
	iserr = 1;
    }
    if ((iserr))
      exitOptionError(argp[0], "requires an input format [fastq|sam|bam] as argument.");
    break;
#endif
  case 'm':
    n = parseOption(&iv, "-m", OPTYP_INT, narg, argp);
    if (n != OPTERR_UNKNOWN) {
      if (n == OPTERR_NOARG  || iv<0 || iv>SHRT_MAX) 
	exitArgumentRangeError(argp[0], 0, SHRT_MAX);
      mp->minscore = (short) iv;
      optflp->minscore = 1;
      mp->flags |= MENUFLAG_MINSCOR;
    }
    break;
  case 'n':
    n = parseOption(&iv, "-n", OPTYP_INT, narg, argp);
    if (n != OPTERR_UNKNOWN) {
      if (n == OPTERR_NOARG  || iv<1 || iv>SHRT_MAX) 
	exitArgumentRangeError(argp[0], 0, SHRT_MAX);
      mp->nthread = (short) iv;
    }
    break;     
  case 'o':
    n = parseOption(&(mp->oufilnam), "-o", OPTYP_STRING, narg, argp);
    if (n == OPTERR_NOARG) exitOptionError(argp[0], "requires a file name as argument.");
    break;
  case 'q':
    n = parseOption(&iv, "-q", OPTYP_INT, narg, argp);
    if (n != OPTERR_UNKNOWN) {
      if (n == OPTERR_NOARG  || iv<0 || iv>OPTLIMIT_MINBASQ) 
	exitArgumentRangeError(argp[0], 0, OPTLIMIT_MINBASQ);
      mp->minbasq = (UCHAR) iv;
      optflp->minbasq = 1;
    }
    break; 
  case 't': 
    n = parseOption(&iv, "-t", OPTYP_INT, narg, argp);
    if (n != OPTERR_UNKNOWN) {
      if (n == OPTERR_NOARG  || iv<0 || iv>INT_MAX) 
	exitArgumentRangeError(argp[0], 0, INT_MAX);
      mp->mapqmin = iv;
      optflp->mapqmin = 1;
    }
    break;
#ifdef HAVE_BAMBAMC
  case 'T': /* specify directory for temporary output files */
    n = parseOption(&(mp->tmpdirnam), "-T", OPTYP_STRING, narg, argp);
    if (n == OPTERR_NOARG) exitOptionError(argp[0], "requires a directory name as argument.");
    break;
#endif 
  case 'u': 
    n = parseOption(&iv, "-u", OPTYP_INT, narg, argp);
    if (n != OPTERR_UNKNOWN) {
      if (n == OPTERR_NOARG  || iv<0 || iv>INT_MAX) 
	exitArgumentRangeError(argp[0], 0, INT_MAX);
      mp->readskip = iv;
      optflp->readskip = 1;
    }
    break; 
  default:
    n=OPTERR_UNKNOWN;
    break;
  }

  return n;
}

static int fprintHashOptions(FILE *fp, const MenuOpt *menup)
{
  const INDEXMENU *ip = (INDEXMENU *) menup->paramp;
  if (menup->subprog !=  MENU_INDEX)
    return ERRCODE_ASSERT;

  fprintf(fp, "# [-k] kmer: %d\n", ip->kmer);
  fprintf(fp, "# [-s] skip: %d\n", ip->skip);
  return ERRCODE_SUCCESS;
}

static int fprintMapOptions(FILE *fp, const MenuOpt *menup)
{
  const MAPMENU *mp = (MAPMENU *) menup->paramp;
  if (menup->subprog !=  MENU_MAP)
    return ERRCODE_ASSERT;

  fprintf(fp, "# ncut: %d\n", mp->ncut);
  fprintf(fp, "# maxhit: %d\n", mp->maxhit);

  return ERRCODE_SUCCESS;
}

static void fprintCommandLine(FILE *fp, int argc, char *argv[])
{
  int i;
  fprintf(fp, "# Command line: ");
  for (i=0; i<argc; i++) {
    fprintf(fp, " %s", argv[i]);
  }
  fprintf(fp, "\n");
}
/******************************************************************************
 ********************** Public Methods of Type MenuOpt ************************
 ******************************************************************************/

MenuOpt *menuCreate(void)
{
  MenuOpt *mp;
  EMALLOCP0(mp);
  if (mp) {
    if (NULL == ESTRING_INIT(mp->estrbuf)) {
      menuDelete(mp);
      mp = NULL;
    }
    mp->indir = NULL;
  }
  return mp;
}

void menuDelete(MenuOpt *mp)
{
  if (mp) {
    ESTRING_FREE(mp->estrbuf);
    if (mp->subprog == MENU_MAP || mp->subprog == MENU_SAMPLE) {
      MAPMENU *p = (MAPMENU *) mp->paramp;
      free(p->oufilnam);
      free(p->insfilnam);
      free(p->tmpdirnam);
    }
    free(mp->paramp);
    free(mp->indir);
  }
  free(mp);
}

int menuParseCommandLine(MenuOpt *menu, 
			 int argc, char *argv[])
{
  //char *prognam = argv[0];
  char *subprognam = NULL;
  int ntokens = 0;
  struct _optflags optfl;
  static MENU_PARSER *parse_fn=NULL;
  static MENU_CHECKER *check_fn=NULL;
  const TASKDOC *taskdocp = NULL;

  menu->argc = argc;
  menu->argv = argv;
  if (argc < 2) {
    fprintf(stdout, "\n%s\n", MENU_PROGNAM);
    fprintf(stdout, MENU_PROGNAM_VERSION_FMT, MENU_RELEASE_VERSION);
    
    fprintf(stdout, MENU_USAGE_SUMMARY);
    return ERRCODE_MENUE;
  }

  memset(&optfl,0,sizeof(struct _optflags));

  subprognam = argv[1];
  if (!strcmp(subprognam, "index")) {
    taskdocp = &MENU_TASKDOC_INDEX;
    if (argc > 2 && argv[2][0] == MENU_OPTION_INDICATOR &&
	argv[2][1] == 'H') {
      menu->subprog = MENU_HELP;
      fprintTaskDoc(stdout, taskdocp, 1);
      return ERRCODE_SUCCESS;
    }
    menu->subprog = MENU_INDEX;
    fprintCommandLine(stderr, argc, argv);
    initHashOptions(menu);
    parse_fn = parseHashOptions;
    check_fn = checkHashDefaults;
  } else if (!strcmp(subprognam, "map")) {
    taskdocp = &MENU_TASKDOC_MAP;
    if (argc > 2 && argv[2][0] == MENU_OPTION_INDICATOR &&
	argv[2][1] == 'H') {
      menu->subprog = MENU_HELP;
      fprintf(stdout, MENU_USAGE_MAP_HEADER);
      fprintTaskDoc(stdout, taskdocp, 1);
      return ERRCODE_SUCCESS;
    }
    menu->subprog = MENU_MAP;
    fprintCommandLine(stderr, argc, argv);
    initMapOptions(menu);
    parse_fn = parseMapOptions;
    check_fn = checkMapDefaults;
  } else if (!strcmp(subprognam, "check")) {
    menu->subprog = MENU_CHECK;
    taskdocp = &MENU_TASKDOC_CHECK;
    if (argc > 2 && argv[2][0] == MENU_OPTION_INDICATOR &&
	argv[2][1] == 'H') {
      menu->subprog = MENU_HELP;
      fprintTaskDoc(stdout, taskdocp, 1);
      return ERRCODE_SUCCESS;
    }
    fprintCommandLine(stderr, argc, argv);
 } else if (!strcmp(subprognam, "sample")) {
    menu->subprog = MENU_SAMPLE;
    taskdocp = &MENU_TASKDOC_SAMPLE;
    if (argc > 2 && argv[2][0] == MENU_OPTION_INDICATOR &&
	argv[2][1] == 'H') {
      menu->subprog = MENU_HELP;
      fprintTaskDoc(stdout, taskdocp, 1);
      return ERRCODE_SUCCESS;
    }
    fprintCommandLine(stderr, argc, argv);
    initMapOptions(menu);
    parse_fn = parseSampleOptions;
    check_fn = checkSampleDefaults;     
  } else if (!strcmp(subprognam, "help")) {
    menu->subprog = MENU_HELP;
    fprintf(stdout, "\n%s\n\n", MENU_PROGNAM);
    fprintf(stdout, MENU_USAGE_SUMMARY);
    fprintf(stdout, "DESCRIPTION:\n%s", MENU_SHORT_DESCRIPTION);
    return ERRCODE_SUCCESS;
  } else if (!strcmp(subprognam, "version")){
    menu->subprog = MENU_VERSION;
    printBlurb(stdout);
    return ERRCODE_SUCCESS;
  } else {
    menu->subprog = MENU_UNKNOWN;
    printf("ERROR: unknown task switch %s\n\n", subprognam);
    fprintf(stdout, MENU_USAGE_SUMMARY);
    return ERRCODE_MENUE;
  }
  if ( menu->subprog == MENU_UNKNOWN)
    exitOptionError(subprognam, "not yet implemented.");

  argc -= 2;
  argv += 2;
  menu->ninfil=0;

  while (argc > 0)
    if (argv[0][0] == MENU_OPTION_INDICATOR) {
      if (menu->ninfil>0) /* already in list of file names */
	exitOptionError(argv[0], "occurs after file name.");

      if (parse_fn)
	ntokens = (*parse_fn)(menu, &optfl, argc, argv);
      if (ntokens <= 0) {
	if (ntokens == OPTERR_UNKNOWN)
	  exitOptionError(argv[0], "unknown.");
	return ERRCODE_FAILURE;
      }
      argc -= ntokens;
      argv += ntokens;
    } else { /* treat as input file name */
      if (menu->ninfil < 1) menu->filnams = argv;
      menu->ninfil++;
      argv++;
      argc--;
    }
  
  if (menu->ninfil < 2 
      && !(menu->subprog == MENU_CHECK && menu->ninfil == 1)) {
    fprintf(stdout, "no input files specified.\n");
    fprintTaskDoc(stdout, taskdocp, 0);
    exit(EXIT_FAILURE);
  }

  if (menu->ninfil > 3)
    exitOptionError(NULL, "too many input files specified.");
  
  if (check_fn && 
      ((*check_fn)(menu, &optfl)))
      return ERRCODE_FAILURE;
  return ERRCODE_SUCCESS;
}

void menuPrint(FILE *fp, const MenuOpt *menup)
{
  fprintf(fp, "\n========= Menue =========\nSubprogram: ");
  switch (menup->subprog) {
  case MENU_INDEX:
    fprintf(fp, "index\n");
    fprintHashOptions(fp, menup);
    break;
  case MENU_MAP:
    fprintf(fp, "map\n");
    fprintMapOptions(fp, menup);
    break;
  case MENU_UNKNOWN:
  default:
    fprintf(fp, "unknown\n");
    break;
  }
  fprintf(fp, "\nFiles\n");
  fprintf(fp, "dir: %s\n", (menup->indir)? menup->indir: "-");
  fprintf(fp, "infil: %d\n", menup->ninfil);

  fprintf(fp, "====== End of Menue =====\n");
}

char menuGetSubProgTyp(const MenuOpt *mp)
{
  return mp->subprog;
}

char **menuGetCommandLine(const MenuOpt *mp, int *argc)
{
  if (argc) *argc = mp->argc;
  return mp->argv;
}

const char *menuGetProgramName(const char **version)
{
  if (version) *version = MENU_RELEASE_VERSION;
  return MENU_PROGNAM_SHORT;
}

int menuGetHashParams(const MenuOpt *mp, const char **hashnamp, UCHAR *kmer, UCHAR *nskip)
{
  const INDEXMENU *ip = (INDEXMENU *) mp->paramp;
  if (mp->subprog != MENU_INDEX)
    return ERRCODE_FAILURE;
  if (kmer) *kmer = ip->kmer;
  if (nskip) *nskip = ip->skip;
  if (hashnamp) 
    *hashnamp = (mp->ninfil > 0)? mp->filnams[0]: NULL;
    
  return ERRCODE_SUCCESS;
}

int menuGetMapParams(const MenuOpt *menup, const char **hashnamp,
		     int *nhitmax_tuple, double *tupcovmin,
		     short *minscore,
		     short *scorediff,
		     UCHAR *minbasq,
		     double *minidentity,
		     int *randseed,
		     int *readskip,
		     int *insert_min, int *insert_max, 
		     UCHAR *outform,
		     const char **outfilnam,
		     const char **insfilnam,
		     UCHAR *pairtyp)
{
  const MAPMENU *mp = (MAPMENU *) menup->paramp;
  if (menup->subprog != MENU_MAP 
      && menup->subprog != MENU_SAMPLE
      )
    return ERRCODE_FAILURE;
  if (hashnamp)
    *hashnamp = ((	
		  MENU_SAMPLE == menup->subprog ||		 
		  MENU_MAP == menup->subprog
		  ) 
		 && menup->ninfil > 1)?
      menup->filnams[0]: NULL;
  if (nhitmax_tuple) *nhitmax_tuple = mp->ncut;
  if (tupcovmin) *tupcovmin =  mp->mincover;
  if (minscore) *minscore = mp->minscore;
  if (scorediff) *scorediff = mp->scorediff;
  if (minbasq) *minbasq = mp->minbasq;
  if (minidentity) *minidentity = mp->minidentity;
  if (randseed) *randseed = mp->randseed_repeat;
  if (readskip) *readskip = mp->readskip;
  if (insert_min) *insert_min = mp->insert_range[0];
  if (insert_max) *insert_max = mp->insert_range[1];
  if (outform) *outform = mp->oufmt.typ;
  if (outfilnam) *outfilnam = mp->oufilnam;
  if (insfilnam) *insfilnam = mp->insfilnam;
  if (pairtyp) *pairtyp = mp->pairtyp;
  return ERRCODE_SUCCESS;
}

int menuTestMapOutputFormatFlags(const MenuOpt *menup, const UCHAR outform, const MENUOFMTFLG_t flags)
{
  const MAPMENU *mp = (MAPMENU *) menup->paramp;
  if (menup->subprog != MENU_MAP || mp->oufmt.typ != outform)
    return 0;
  return (mp->oufmt.flags & flags) == flags;
}

int menuGetMapInputFormat(const MenuOpt *menup, 
			  unsigned char *infmt, const char **tmpdirnam)
{
  const MAPMENU *mp = (MAPMENU *) menup->paramp;

  if (menup->subprog != MENU_MAP 
      && menup->subprog != MENU_SAMPLE
      )
    return ERRCODE_FAILURE;
  
  if (infmt) *infmt = mp->inform;
  if (tmpdirnam) *tmpdirnam = mp->tmpdirnam;

  return ERRCODE_SUCCESS;
}

uint8_t menuGetMapPenaltyScores(const MenuOpt *menup, int8_t *match, int8_t *subst, 
				int8_t *gap_open, int8_t *gap_ext)
{
  const MAPMENU *mp = (MAPMENU *) menup->paramp;
  if (menup->subprog != MENU_MAP 
      && menup->subprog != MENU_SAMPLE)
    return ERRCODE_FAILURE;

  if (match) *match = mp->penalties[MENUPNLTYP_MATCH];
  if (subst) *subst = mp->penalties[MENUPNLTYP_MISMATCH];
  if (gap_open) *gap_open = mp->penalties[MENUPNLTYP_GAPINIT];
  if (gap_ext) *gap_ext = mp->penalties[MENUPNLTYP_GAPEXT];

  return mp->penaltyflags;
}

MENUFLG_t menuGetFlags(const MenuOpt *menup)
{
  MENUFLG_t flags = 0;
  if ((menup) && (menup->paramp)) {
    if (menup->subprog == MENU_MAP 
	|| menup->subprog == MENU_SAMPLE
	)
      flags = ((const MAPMENU *) menup->paramp)->flags;
  }
  
  return flags;
}

short menuGetNumberOfThreads(const MenuOpt *menup)
{
  const MAPMENU *mp = (MAPMENU *) menup->paramp;
  return mp->nthread;
}

int menuGetFileNames(const MenuOpt *mp, 
		     const char **filnam1p, 
		     const char **filnam2p)
{
  if (filnam1p) *filnam1p = 0;
  if (filnam2p) *filnam2p = 0;

  if (mp->subprog == MENU_MAP 
      || mp->subprog == MENU_SAMPLE
      ) {
    if ((filnam1p) && mp->ninfil > 1) *filnam1p = mp->filnams[1];
    if ((filnam2p) && mp->ninfil > 2) *filnam2p = mp->filnams[2];
  } else if (mp->subprog == MENU_INDEX) {
    if ((filnam1p) && mp->ninfil > 1) *filnam1p = mp->filnams[1];
  }
  else if (mp->subprog == MENU_CHECK) {
    if ((filnam1p) && mp->ninfil > 0) *filnam1p = mp->filnams[0];
    if ((filnam2p) && mp->ninfil > 1) *filnam2p = mp->filnams[1];
  }

  return mp->ninfil;
}

void menuPrintWallClockTime(FILE *fp, time_t time_start, time_t time_stop, const char *headerp)
{
  double secs = difftime(time_stop, time_start);
  short days = (short) (secs/TIME24HRSINSECS);
  short hours = (short) ((secs - days*TIME24HRSINSECS)/3600);
  short mins = (short) ((secs - days*TIME24HRSINSECS - hours*3600)/60);
  double seconds = (secs - days*TIME24HRSINSECS - hours*3600 - mins*60);

  if (NULL == headerp)
    fprintf(fp, "# %s: Total elapsed wall clock time: ", MENU_PROGNAM_SHORT);
  else
    fprintf(fp, "# %s: %s:", MENU_PROGNAM_SHORT, headerp);

  if (days > 0)
    fprintf(fp, "%hi days ", days);
  if (hours > 0)
    fprintf(fp, "%hi hours ", hours);
  if (mins > 0)
    fprintf(fp, "%hi minutes and ", mins);
  fprintf(fp, "%g seconds\n", seconds);
}
