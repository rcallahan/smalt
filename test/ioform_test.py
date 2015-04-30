# tests of input and output formats
PROGNAM = "../src/smalt"
REF_FASTA_NAME = "hs37chrXtrunc.fa.gz"
READ_PREFIX = "hs37l100i300e05q_trunc"

TMPFIL_PREFIX = "TMP"

KMER = 11
NSKIP = 11

SAMTAG_PG = "@PG"

def smalt_index(df, index_name, fasta_name, kmer, nskip):
    from sys import exit
    from subprocess import call

    tup = (PROGNAM, 'index',
           '-k', '%i' % (int(kmer)),
           '-s', '%i' % (int(nskip)),
           index_name,
           fasta_name)
    df.call(tup, "when indexing")
                
def smalt_map(df, oufilnam, indexnam, readfil, typ="fastq", matefil=""):
    from sys import exit
    from subprocess import call
    
    tup = (PROGNAM, 'map',
           '-r', '-1',
           '-F', typ,
           '-o', oufilnam,
           indexnam,
           readfil, matefil)

    df.call(tup, "when mapping")

if __name__ == '__main__':
    from testdata import DataFiles, areFilesIdentical
    
    df = DataFiles()
    refnam = df.joinData(REF_FASTA_NAME)
    readnamA = df.joinData(READ_PREFIX + "_nonam_1.fq.gz")
    readnamB = df.joinData(READ_PREFIX + "_nonam_2.fq.gz")
    bamnam = df.joinData(READ_PREFIX + ".bam")
    samnam = df.unpack(READ_PREFIX + ".sam")
    indexnam = df.addIndex(TMPFIL_PREFIX)
    
    oufilnam_fastq = df.addTMP(TMPFIL_PREFIX + "fastq.sam")
    oufilnam_bam = df.addTMP(TMPFIL_PREFIX + "bam.sam")
    oufilnam_sam = df.addTMP(TMPFIL_PREFIX + "sam.sam")
    oufilnam2_sam = df.addTMP(TMPFIL_PREFIX + "sam2.sam")
    
    smalt_index(df,indexnam, refnam, KMER, NSKIP)
    smalt_map(df,oufilnam_fastq, indexnam, readnamA, "fastq", readnamB)
    smalt_map(df,oufilnam_bam, indexnam, bamnam, "bam")
    
    if not areFilesIdentical(oufilnam_fastq, oufilnam_bam, SAMTAG_PG):
        exit("Output for FASTQ input and BAM input does not agree!")

    #print "SAM and FASTQ comparison ok."
    
    smalt_map(df,oufilnam_sam, indexnam, samnam, "sam")
    if not areFilesIdentical(oufilnam_sam, oufilnam_bam, SAMTAG_PG):
        exit("Output for SAM input and BAM input does not agree!")
        
    #print "SAM and BAM comparison ok."
    # use smalt output as input
    smalt_map(df,oufilnam2_sam, indexnam, oufilnam_sam, "sam")
    if not areFilesIdentical(oufilnam2_sam, oufilnam_bam, SAMTAG_PG):
        exit("Using smalt sam output as input gave inconsistent results!")

    #print "Test ok."
    
    df.cleanup()
    exit()
