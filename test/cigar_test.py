# test cigar strings

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

PROGNAM = "../src/smalt"
FNAM_REF = "cigar_ref.fa.gz"
FNAM_READ1 = "cigar_read1.fq"
FNAM_READ2 = "cigar_read2.fq"

TMPFIL_PREFIX = "TMPcig"

KMER = 13
NSKIP = 2

def smalt_index(df,index_name, fasta_name, kmer, nskip):
    from sys import exit
    from subprocess import call

    tup = (PROGNAM, 'index',
           '-k', '%i' % (int(kmer)),
           '-s', '%i' % (int(nskip)),
           index_name,
           fasta_name)
    df.call(tup, "when indexing")

def smalt_map(df, oufilnam, indexnam, readfil, matefil, typ="fastq", flags=[]):
    from sys import exit
    from subprocess import call
 
    tup = [PROGNAM, 'map']
    if len(flags) > 0:
        tup.extend(flags)
    tup.extend([
           '-f', typ,
           '-o', oufilnam,
           indexnam,
           readfil, matefil])
    df.call(tup, "when mapping")

if __name__ == '__main__':
    from testdata import DataFiles
    
    df = DataFiles()
    
    refnam = df.joinData(FNAM_REF)
    readnamA = df.joinData(FNAM_READ1)
    readnamB = df.joinData(FNAM_READ2)
    indexnam = df.addIndex(TMPFIL_PREFIX)
    oufilnam = df.addTMP(TMPFIL_PREFIX + ".sam")
    
    smalt_index(df,indexnam, refnam, KMER, NSKIP)
    smalt_map(df,oufilnam, indexnam, readnamA, readnamB, "sam", ["-x"])
    
    #print "Test ok."
    
    df.cleanup()
    exit()
