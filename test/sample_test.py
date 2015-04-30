# test sampling of insert lengths

PROGNAM = "../src/smalt"
REF_FASTA_NAME = "genome_1.fa.gz"
READ_PREFIX = "gen1l75i300e0"

KMER = 11
NSKIP = 5

TMPFIL_PREFIX = "TMP"
MAXNUM_NONPROPER_PAIRS = 20

def smalt_index(df, index_name, fasta_name, kmer, nskip):
    from sys import exit
    from subprocess import call

    tup = (PROGNAM, 'index',
           '-k', '%i' % (int(kmer)),
           '-s', '%i' % (int(nskip)),
           index_name,
           fasta_name)
    df.call(tup, "when indexing")

def smalt_sample(df,oufilnam, indexnam, readfil, matefil, option=[]):
 
    tup = [PROGNAM, 'sample']
    
    if option:
        tup.extend(option)
    tup.extend(
        ['-o', oufilnam,
         indexnam,
         readfil, matefil]
        )

    df.call(tup, "when sampling")

def smalt_check(df, readfil, matefil=""):
    tup = [PROGNAM, 'check', readfil]
    
    if len(matefil) > 0:
        tup.append(matefil)
    df.call(tup, "when checking")

def smalt_map(df, oufilnam, indexnam, samplenam, readfil, matefil, option=[]):
    from sys import exit
    from subprocess import call
 
    tup = [PROGNAM, 'map']
    
    if option:
        tup.extend(option)
    tup.extend(
        ['-g', samplenam,
        '-f', 'cigar',
        '-o', oufilnam,
         indexnam,
         readfil, matefil]
        )
    df.call(tup, "when mapping")

def assess_mapping(oufilnam):
    from formats import Cigar, getNextCigarPair, openFile

    infil = openFile(oufilnam, 'r')
    
    cigA = Cigar()
    cigB = Cigar()

    pair_ctr = 0
    nonproper_ctr = 0
    (isOK, isEOF) = getNextCigarPair(infil, cigA, cigB)
    while isOK:
        pair_ctr = pair_ctr + 1
        if cigA.mapcls != 'A':
            nonproper_ctr = nonproper_ctr + 1
        if cigB.mapcls != 'A':
            nonproper_ctr = nonproper_ctr + 1
        (isOK, isEOF) = getNextCigarPair(infil, cigA, cigB)
    infil.close()
    
    if pair_ctr != 10000:
        exit("Found %i pairs, but expected 10,000." % pair_ctr)
    if nonproper_ctr > MAXNUM_NONPROPER_PAIRS:
        exit("Found %i non-proper pairs. Expected max. %i" %
             (nonproper_ctr, MAXNUM_NONPROPER_PAIRS))

def compare_mapping(oufilnam1, oufilnam2):
    from formats import Cigar, openFile

    infil1 = openFile(oufilnam1, 'r')
    infil2 = openFile(oufilnam2, 'r')
    
    cig1 = Cigar()
    cig2 = Cigar()

    ctr1 = 0
    ctr2 = 0
    while 1:
        if cig1.next(infil1):
            break
        ctr1 = ctr1 + 1
        
        if cig2.next(infil2):
            break
        ctr2 = ctr2 + 1

        if cmp(cig1.qnam, cig2.qnam):
            exit("readnames don't match: '%s' vs '%s'" % \
                 (cig1.qnam, cig2.qnam))
        if cmp(cig1,cig2) and cig1.mapq > 5 and cig2.mapq > 5:
            exit("mappings don't match for read '%s'" % \
                 cig1.qnam)
    infil2.close()
    infil1.close()

    if ctr1 != ctr1:
        exit("Expected the same number of mates, got %i (A) and %i (B)" % \
             ctr1, ctr2)
    if ctr1 != 20000:
        exit("Expected 20,000 reads, got %i." % ctr1)
    
if __name__ == '__main__':
    from testdata import DataFiles
    
    df = DataFiles()
    
    refnam = df.joinData(REF_FASTA_NAME)
    readnamA = df.joinData(READ_PREFIX + "_1.fq.gz")
    readnamB = df.joinData(READ_PREFIX + "_2.fq.gz")
    indexnam = df.addIndex(TMPFIL_PREFIX)

    samplenam1 = df.addTMP(TMPFIL_PREFIX + ".1.txt")
    samplenam2 = df.addTMP(TMPFIL_PREFIX + ".2.txt")
    
    oufilnam1 = df.addTMP(TMPFIL_PREFIX + ".1.cig")
    oufilnam2 = df.addTMP(TMPFIL_PREFIX + ".2.cig")
    oufilnam3 = df.addTMP(TMPFIL_PREFIX + ".3.cig")

    smalt_check(df,readnamA, readnamB)
    smalt_index(df,indexnam, refnam, KMER, NSKIP)
    smalt_sample(df,samplenam1, indexnam, readnamA, readnamB)
    smalt_map(df,oufilnam1, indexnam, samplenam1, readnamA, readnamB)
    assess_mapping(oufilnam1)

    nthread_tup = ['-n', '4']
    
    smalt_sample(df,samplenam2, indexnam, readnamA, readnamB, nthread_tup)
    smalt_map(df,oufilnam2, indexnam, samplenam2, readnamA, readnamB, nthread_tup)
    assess_mapping(oufilnam2)

    nthread_tup = ['-n', '4', '-O']
    smalt_map(df,oufilnam3, indexnam, samplenam2, readnamA, readnamB, nthread_tup)
    compare_mapping(oufilnam1, oufilnam3)
    
    df.cleanup()
    exit(0)
              
