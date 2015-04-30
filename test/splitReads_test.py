# test mapping of split reads (-p flag)

PROGNAM = "../src/smalt"
REF_FASTA_NAME = "genome_1.fa.gz"

KMER = 11
NSKIP = 4

READ_PAIRS = (("ATAAGAACAATTTCCATAAATTATTTGAAGATTATTTAGAAGATTATAAT"\
              "AAGAAGTTTTTTCCACTTAAAAATATATTTATATATATATATATATATTT",
              "GTACGTATTTATATATTTTAGAATAGAGAATTAATGAATAAATGCAAAGA"\
              "ATTCCATAAAAAATATAATGTAGAAATGGTTAAGAAAGTTACCGAATAAT"),)
#cigar:A:36 SPLIT_000000/1 100 61 - MAL14 1060107 1060146 + 40 M 40 
#cigar:A:36 SPLIT_000000/2 28 100 + MAL14 1059974 1060045 + 65 M 60 I 1 M 12 
#cigar:P:47 SPLIT_000000/1 1 58 + MAL14 2573035 2573092 + 58 M 58 
#cigar:P:37 SPLIT_000000/2 25 1 - MAL9 1418385 1418409 + 25 M 25
MAPPED_CIGAR = (('A', 1, (100, 61), 'MAL14', (1060107, 1060146)),
                ('A', 2, (28, 100), 'MAL14', (1059974, 1060045)),
                ('P', 1, (1, 58), 'MAL14', (2573035, 2573092)),
                ('P', 2, (25, 1), 'MAL9', (1418385, 1418409)))
                

FASTA_FMT = ">SPLIT_%6.6i/%1i\n%s\n"
TMPFIL_PREFIX = "tmp"

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

def makeFastqPair(filnamA, filnamB):
    from testdata import openFile

    oufilA = openFile(filnamA, "w")
    oufilB = openFile(filnamB, "w")

    pairctr = 0
    for read_pair in READ_PAIRS:
        oufilA.write(FASTA_FMT % (pairctr, 1, read_pair[0]))
        oufilB.write(FASTA_FMT % (pairctr, 2, read_pair[1]))
        pairctr = pairctr + 1

    oufilB.close()
    oufilA.close()

def checkOutput(cigfilnam, expected_tup):
    from formats import Cigar, openFile
    n_tup = len(expected_tup)
    cig = Cigar()
    okflgs = [False]*n_tup
    
    infil = openFile(cigfilnam)
    linctr = 0
    allok = True

    while not cig.next(infil):
        if linctr > n_tup:
            allok = False
            print "ERROR: %i cigar lines, expected %i" % \
                  (linctr, n_tup)
            break
        mateno = cig.getMateNo()
        observed_tup = (cig.mapcls, mateno, cig.qseg, cig.snam, cig.sseg)
        if observed_tup not in expected_tup:
            print "ERROR: unexpected tuple: ", observed_tup
            allok = False
            break
        
        okflgs[expected_tup.index(observed_tup)] = True  
        linctr = linctr + 1
        
    for i in range(n_tup):
        if not okflgs[i]:
            print "ERROR: could not find tuple: ", expected_tup[i]
            allok = False

    if not allok:
        exit("ERROR when checking. Test not passed.")
    
if __name__ == '__main__':
    from testdata import DataFiles
    
    df = DataFiles()
    
    refnam = df.joinData(REF_FASTA_NAME)
    indexnam = df.addIndex(TMPFIL_PREFIX)
    readfilA = df.addTMP(TMPFIL_PREFIX + "_1.fa")
    readfilB = df.addTMP(TMPFIL_PREFIX + "_2.fa")
    oufilnam = df.addTMP(TMPFIL_PREFIX + ".cig")

    makeFastqPair(readfilA, readfilB)
    smalt_index(df,indexnam, refnam, KMER, NSKIP)
    smalt_map(df,oufilnam, indexnam, readfilA, readfilB, ['-p'])
    checkOutput(oufilnam, MAPPED_CIGAR)
    df.cleanup()
    
    exit(0)
