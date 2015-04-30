# test of multi-threaded execution

PROGNAM = "../src/smalt"

#DATA = ("contigs.fa.gz", "contigs_i3k_1.fq.gz", "contigs_i3k_2.fq.gz")
DATA = ("genome_1.fa.gz", "gen1l75i300e0_1.fq.gz", "gen1l75i300e0_2.fq.gz", int(10000))

KMER = 11
NSKIP = 11
NTHREADS = 2

TMPFIL_PREFIX = "TMPmthread"
MAPQ_THRESH = 6
VERBOSE = False

def smalt_index(df, index_name, fasta_name, kmer, nskip):
    from sys import exit
    from subprocess import call

    tup = (PROGNAM, 'index',
           '-k', '%i' % (int(kmer)),
           '-s', '%i' % (int(nskip)),
           index_name,
           fasta_name)
    df.call(tup, "when indexing")

def smalt_map(df, oufilnam, indexnam, readfil,
              matefil="", added_options=[], ouform="cigar",):
    from sys import exit
    from subprocess import call
    
    tup = [PROGNAM, 'map'];
    tup.extend(added_options);
    tup.extend(['-f', ouform,
                '-o', oufilnam,
                indexnam,
                readfil, matefil])

    df.call(tup, "when mapping")

def cmpCigarFiles(cigfilA, cigfilB, is_verbose=True):
    from formats import Cigar, openFile, getNextCigarPair
    cigA1 = Cigar()
    cigA2 = Cigar()
    cigB1 = Cigar()
    cigB2 = Cigar()
    
    filA = openFile(cigfilA, 'r')
    filB = openFile(cigfilB, 'r')
    ctr = 0
    while 1:
        (isOK, isEOF) = getNextCigarPair(filA, cigA1, cigA2)
        if not isOK:
            break
        (isOK, isEOF) = getNextCigarPair(filB, cigB1, cigB2)
        if not isOK:
            break
        if cigA1 != cigB1:
            if is_verbose:
                print "Not matching:\n%s\n%s" % (cigA1.lin, cigB1.lin)
            if cigA1.mapq > MAPQ_THRESH and cigB1.mapq > MAPQ_THRESH:
                exit("Discrepancy:\n%s\n%s" % (cigA1.lin, cigB1.lin))
        if cigA2 != cigB2:
            if is_verbose:
                print "Not matching:\n%s\n%s" % (cigA2.lin, cigB2.lin) 
            if cigA2.mapq > MAPQ_THRESH and cigB2.mapq > MAPQ_THRESH:
                exit("Discrepancy:\n%s\n%s" % (cigA2.lin, cigB2.lin))
        ctr = ctr + 1
    if not isOK and isEOF:
        isOK = True
    return isOK, ctr
            
        
        
if __name__ == '__main__':
    from testdata import DataFiles, areFilesIdentical

    df = DataFiles()
    refnam = df.joinData(DATA[0])
    readfilnamA = df.joinData(DATA[1])
    readfilnamB = df.joinData(DATA[2])
    n_pairs_expected = DATA[3]
    
    indexnam = df.addIndex(TMPFIL_PREFIX)
    oufilnam_ref = df.addTMP(TMPFIL_PREFIX + ".n0.cig")
    oufilnam_thread = df.addTMP(TMPFIL_PREFIX + ".n%i.cig" % NTHREADS)
    
    smalt_index(df, indexnam, refnam, KMER, NSKIP);

    smalt_map(df, oufilnam_ref, indexnam, readfilnamA, readfilnamB);
    smalt_map(df, oufilnam_thread, indexnam, readfilnamA, readfilnamB,
              ["-n", "%i" % NTHREADS, "-O"]);

    isOK, pairctr = cmpCigarFiles(oufilnam_ref, oufilnam_thread, VERBOSE)
    if VERBOSE:
        print "Test ok=%s, number of pairs = %i" % (isOK, pairctr)
    if not isOK or pairctr != n_pairs_expected:
        exit("Using smalt in multi-threaded mode gave inconsistent results!")

    df.cleanup()
    exit()
