# test splitting of reads aligning across consecutive reference
# sequences. This is relevant when the number of reference sequences >
# SMALT_MAX_REFSEQ_NUM (= 512). In that case the reads are aligned
# against concatenated reference sequences and the index of the
# reference sequence (i.e. name, length) is worked out a-posteriory.
# Splitting the alignment string my result in fragments ending in
# mismatches etc. which have to be cleaned up.

#############################################################################
#############################################################################
#                                                                           #
#  Copyright (C) 2013 - 2014 Genome Research Ltd.                           # 
#                                                                           #
#  Author: Hannes Ponstingl (hp3@sanger.ac.uk)                              #
#                                                                           #
#  This file is part of SMALT.                                              #
#                                                                           #
#  SMALT is free software: you can redistribute it and/or modify it under   #
#  the terms of the GNU General Public License as published by the Free     #
#  Software Foundation, either version 3 of the License, or (at your        #
#  option) any later version.                                               #
#                                                                           #
#  This program is distributed in the hope that it will be useful, but      #
#  WITHOUT ANY WARRANTY; without even the implied warranty of               #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU         #
#  General Public License for more details.                                 #
#                                                                           #
#  You should have received a copy of the GNU General Public License along  #
#  with this program.  If not, see <http://www.gnu.org/licenses/>.          #
#                                                                           #
#############################################################################
#############################################################################

from sys import path
path.append('../misc') # find SAM.py 

PROGNAM = "./smalt_Xali_test"
PROGNAM_REVCOMP = "./sequenceReverseComplement_test"

REFSEQ = (
"ggggaaacagttgccttctcagaaccgtgcggacttataaaaaatgaaagtcagaaactt"\
"gctgctctcggtgatccgtgcggaagtgcaagtgaaaaacgtttttctaaagaacgtgtc"\
"gatgaatatgataagaaaaaaataaaatgtagtaatagtgaaggagcttgcgctccatat",
"agacgattgcacgtatgcgacaaaaatatggtaaaaatggacacaaataatgatagtaaa"\
"gctaaacatgatttattgttggatgtgtgtcttgcagcaaagtatgaaggagagtcatta",
"aaaacatatcgtgcacaatatgatgaacaatatccttcttctgtttctacttttacaatg"\
"tgtactatgttagcacgaagttttgccgatataggtgatattgtcagaggaagagatttg"\
"tatcgtggcaatagaaagaaaaatgaaacaaaaacagaaagagaaaaattagatgataag"
)

READS = (
"aagaaaaaaataaaatgtagtaatagtgaaggagcttgcgctccatat"\
"agacgattgcacgtatgcgacaaaaatatggtaaaaatgg",
"aagaaaaaaataaaatgtagtaatagtgaaggagcttgcgctccatac"\
"agacgattgcacgtatgcgacaaaaatatggtaaaaatgg",
"aagaaaaaaataaaatgtagtaatagtgaaggagcttgcgctccttgc"\
"agacgattgcacgtatgcgacaaaaatatggtaaaaatgg",
"tgtagtaatagtgaaggagcttgcgctccatat"\
"cgacgattgcacgtatgcgacaaaaatatggtaaaaatgg",
"tgtagtaatagtgaaggagcttgcgctccatat"\
"cgcgattgcacgtatgcgacaaaaatatggtaaaaatgg",  
)

# ref no, pos, cigar, edit distance
RESULTS = ((1, 133, "48M40S", 0),
           (1, 133, "47M41S", 0),
           (1, 133, "44M1X1M42S", 1),
           (2, 2, "34S39M", 0),
           (2, 2, "34S1M1D37M", 1))

TMPFIL_PREFIX = "tmp"
KMER = 11
NSKIP = 2

from re import compile
REFNAMNO = compile("^\S+_(\d+)")

def asFASTA(seqnamprefix, sequences):
    oustr = ""
    ctr = 1
    for seq in sequences:
        oustr = oustr + ">%s_%i\n%s\n" % (seqnamprefix, ctr, seq)
        ctr = ctr + 1
    return oustr

def makeFASTAfile(df, prefix, sequences):
    from testdata import openFile

    filnam = df.addTMP("%s%s.fa" % (TMPFIL_PREFIX, prefix))
    oufil = openFile(filnam, 'w')
    oustr = asFASTA(prefix, sequences)
    oufil.write(oustr)
    oufil.close()
    return filnam

def reverseComplement(df, fastq_name_F, fastq_name_R):

    tup = (PROGNAM_REVCOMP,
           fastq_name_F,
           fastq_name_R)
    df.call(tup, "reverse complement failed")

    return

def smalt_index(df, index_name, fasta_name, kmer, nskip):
    tup = (PROGNAM, 'index',
           '-k', '%i' % (int(kmer)),
           '-s', '%i' % (int(nskip)),
           index_name,
           fasta_name)

    df.call(tup, "when indexing")

def smalt_map(df, oufilnam, indexnam, readfil, option=[]):
    tup = [PROGNAM, 'map', '-f', 'sam:x', '-o', oufilnam]
    
    if option:
        tup.extend(option)
    
    tup.extend([indexnam, readfil])
               
    df.call(tup, "when mapping")

def refnoFromNam(refnam):
    m = REFNAMNO.search(refnam)
    return int(m.group(1))

def checkSAM(df, filnam_sam, results):
    from SAM import Sam, openFile

    sam = Sam()
    infil = openFile(filnam_sam, 'r')

    for res in results:
        sam.next(infil)

        if not sam.ok:
            df.exitErr("ERROR: Could not parse file '%s'" % (filnam_sam))

        refno = refnoFromNam(sam.rname)
        if refno != res[0]:
            df.exitErr("ERROR: wrong reference sequence %i '%s' (target: %i)" % \
                       (refno, sam.rname, res[0]))
        if sam.pos != res[1]:
            df.exitErr("ERROR: wrong reference position %i (target:%i)" % (sam.pos, res[1]))

        if sam.cigar != res[2]:
            df.exitErr("ERROR: wrong CIGAR string '%s' (target:'%s')" % (sam.cigar, res[2]))

        (typ, fld) = sam.tags["NM"]
        edist = int(fld)
        if edist != res[3]:
            df.exitErr("ERROR: wrong edit distance %i (target: %i)" % (edist, res[3]))
                          
    infil.close()
    return

if __name__ == '__main__':
    from testdata import DataFiles
    
    df = DataFiles()

    indexnam = df.addIndex(TMPFIL_PREFIX)
    reffilnam = makeFASTAfile(df, "REF", REFSEQ)
    readfilnam = makeFASTAfile(df, "READ", READS)
    rc_readfilnam = df.addTMP(TMPFIL_PREFIX + "RC.fa")
    samfilnam = df.addTMP(TMPFIL_PREFIX + ".sam")
    rc_samfilnam = df.addTMP(TMPFIL_PREFIX + "RC.sam")
    
    smalt_index(df, indexnam, reffilnam, KMER, NSKIP)
    smalt_map(df, samfilnam, indexnam, readfilnam)
    checkSAM(df, samfilnam, RESULTS)
    
    reverseComplement(df, readfilnam, rc_readfilnam)
    smalt_map(df, rc_samfilnam, indexnam, rc_readfilnam)
    checkSAM(df, samfilnam, RESULTS)
    
    df.cleanup()
    exit(0)
       
