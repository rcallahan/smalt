# test CIGAR output format
# mapping class labels

PROGNAM = "../src/smalt"
RCPROGNAM = "./sequenceReverseComplement_test"

REF_FASTA_NAME = "genome_1.fa.gz"
READ_PREFIX = "gen1l100i500e1"

KMER = 11
NSKIP = 5

TMPFIL_PREFIX = "TMP"

OPTION_LABEL_PAIRS_ORIG = (
    (['-r', '-1'], ("AA", "BB", "DD", "SR", "BB", "RR", "AA")),
    (['-m', '30'], ("AA", "BB", "DD", "??", "BB", "??", "AA")),
    (['-l', 'pp', '-r', '-1'], ("CC", "CC", "DD", "BB", "SR", "RR", "CC")),
    (['-l', 'mp', '-r', '-1'], ("CC", "CC", "DD", "SR", "SR", "RR", "CC")),
)

OPTION_LABEL_PAIRS_RC2ND = (
    (['-r', '-1'], ("CC", "CC", "DD", "BB", "SR", "RR", "CC")),
    (['-m', '30'], ("CC", "CC", "DD", "SN", "??", "??", "CC")),
    (['-l', 'pp', '-r', '-1'], ("AA", "BB", "DD", "SR", "BB", "RR", "AA")),
    (['-l', 'mp', '-r', '-1'], ("CC", "CC", "DD", "SR", "SR", "RR", "CC")),
)

OPTION_LABEL_PAIRS_RC2ND_RC1ST = (
    (['-r', '-1'], ("CC", "CC", "DD", "RS", "RS", "RR", "CC")),
    (['-m', '30'], ("CC", "CC", "DD", "??", "??", "??", "CC")),
    (['-l', 'pp', '-r', '-1'], ("CC", "CC", "DD", "BB", "RS", "RR", "CC")),
    (['-l', 'mp', '-r', '-1'], ("AA", "BB", "DD", "RS", "BB", "RR", "AA")),
)



def checkLabels(cigfilnam, label_pairs, mateno_check=True):
    from formats import Cigar, openFile, getNextCigarPair
    cigA = Cigar()
    cigB = Cigar()
    
    infil = openFile(cigfilnam)
    for lb in label_pairs:
        (isOK, isEOF) = getNextCigarPair(infil, cigA, cigB, mateno_check)
        if (not isOK) or isEOF:
            exit("missing lines in cigar file %s" % cigfilnam)
        if (cigA.mapcls != lb[0] and lb[0] != '?') or \
           (cigB.mapcls != lb[1] and lb[1] != '?'):
            exit("unexpeced cigar mapping labels %s%s (%s) for read pair %s" % \
                 (cigA.mapcls, cigB.mapcls, lb, cigA.qnam))
    return
        
def reverseComplement(df, filnam_in, filnam_out):
    from sys import exit
    from subprocess import call

    df.call([RCPROGNAM, filnam_in, filnam_out],
            "when reverse complementing reads in file '%s'" % (filnam_in))
   
def smalt_index(df, index_name, fasta_name, kmer, nskip):
    from sys import exit
    from subprocess import call

    tup = (PROGNAM, 'index',
           '-k', '%i' % (int(kmer)),
           '-s', '%i' % (int(nskip)),
           index_name,
           fasta_name)
    df.call(tup, "when indexing")

def smalt_map(df,oufilnam, indexnam, readfil, matefil="", option=[]):
    from sys import exit
    from subprocess import call
 
    tup = [PROGNAM, 'map']
    
    if option:
        tup.extend(option)
    tup.extend(
        ['-f', 'cigar',
         '-o', oufilnam,
         indexnam,
         readfil, matefil]
        )
    df.call(tup, "when mapping")

def process(df, indexnam, oufilnam, readnamA, readnamB, option_label_pairs, mateno_check=True):
    for (option, label_pairs) in option_label_pairs:
        smalt_map(df,oufilnam, indexnam, readnamA, readnamB, option)
        checkLabels(oufilnam, label_pairs, mateno_check)

if __name__ == '__main__':
    from testdata import DataFiles
    
    df = DataFiles()
    
    refnam = df.joinData(REF_FASTA_NAME)
    readnamA = df.joinData(READ_PREFIX + "_1.fq")
    readnamB = df.joinData(READ_PREFIX + "_2.fq")
    indexnam = df.addIndex(TMPFIL_PREFIX)
    
    oufilnam = df.addTMP(TMPFIL_PREFIX + ".cig")
    readnamRCA = df.addTMP(TMPFIL_PREFIX + "rc_1.fq")
    readnamRCB = df.addTMP(TMPFIL_PREFIX + "rc_2.fq")
    
    smalt_index(df,indexnam, refnam, KMER, NSKIP)
    process(df,indexnam, oufilnam, readnamA, readnamB, OPTION_LABEL_PAIRS_ORIG)

    #print "reverse complement 2nd read ..."
    reverseComplement(df,readnamB, readnamRCB)
    process(df,indexnam, oufilnam, readnamA, readnamRCB, OPTION_LABEL_PAIRS_RC2ND)

    df.addLog(["reverse complement 1st read and swap with 2nd"])
    # this equates to changing from -l pe to -l mp
    reverseComplement(df,readnamA, readnamRCA)
    process(df,indexnam, oufilnam, readnamRCB, readnamRCA,
            OPTION_LABEL_PAIRS_RC2ND_RC1ST, False)

    df.cleanup()
    exit(0)
