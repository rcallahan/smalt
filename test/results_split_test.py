# test splitting of alignments that span multiple reference sequences
# note: the number of reference sequences has to be > 256 for this to
# to be triggered.

PROGNAM = "../src/smalt"

REFSEQ = "contigs.fa.gz"

READSEQS = ("@SIM_000000000_contig5121_000000001_5120_R_20m/1\n"\
            "AAATATAGTTAGCCACCTCCAAGACACAGTAGAGAAATAAGCATTCTGAT"\
            "GGGAAATGTTTATGCATTTGTGAAGGTTAGGCCTACAAGTGGACCATGTA\n"\
            "+\n"\
            "Bbf^BcfBbdcfBdcfcddYbaB\dfYcc\ec_efeed\_cTebcacced"\
            "bebf`^fa^cV`fde^\ecc`d`cbccfeTbdbeZcefccf^fccbfBc!\n",
            "@SIM_000000001_contig4209_000030120_4208_R_97m/1\n"\
            "ATACATGCATATCCTTTTCAAATGGAGATGAAAGCACTCAAAATTTGATG"\
            "GCCGCCGTTATTATGAGCAGACGACTTGAAATGATAGTACTGCAAAATGA\n"\
            "+\n"\
            "BeUBdfeYTf]d]BeKaLcddcBeff`e\\bf`f^`f_ebVdbcddfe`Yf"\
            "edff^Yc^`dY`\\fefa^ddcYb^ef`adfdB\`f\`c`defce`adde!\n",
            "@SIM_000000002_contig5043_000004873_5042_F_81m/1\n"\
            "TGTAAGCTGTACTTCTGAGCTGTTGATTTTCATTGAAGGTACCTCACTGT"\
            "GCATCTATTGAGATAATCATGTTGTTTTTTGTAAGTACTCCCACTTGGAC\n"\
            "+\n"\
            "Yfffcc\BYBf`dedBeeaZf^c^adff`BbfeTB_bdf\`YfdeBfBdd"\
            "eB_ae\eVff]dTedc_dcff^cVa_^d\\affdfdfddcadbbBT^fBd!\n",
            "@SIM_000000003_contig1732_000000001_1731_F_55m/1\n"\
            "AACGCCTGGCTCCCACCTAATAAAAATACCTATATTAGATGGAGGGTGAC"\
            "CGACACCAGTAACCCCAGCAGTTGAGGCTCACACCCAACAATCTATGTTT\n"\
            "+\n"\
            "]BBdBYecfbfTX[dQddfeafceDa`^^R^cYdaLfbBfYfadffYdOc"\
            "B``bdZffbeMUcUaTfeceb`V^fecaddde\cfW^de[Bececec`B!\n",
            "@SIM_000000004_contig10737_000030839_10736_F_100m/1\n"\
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGAAACAAAAAA"\
            "AAAAGAAAAAAAAAAAAAAAAAACAAAAAAAAAAAAAAAAAAAAAAAAAA\n"\
            "+\n"\
            "daBZddeYdBfe`fBf`ecfdfcf`c\eacXdbLe^b``ccadK]e]^_T"\
            "fccbcd\Tffb^Yfdfdf^efa\[fbfdad`X\Ue_^eBfcfcfNfdf^!\n"
            )

# cigar:S:16 SIM_000000000_contig5121_000000001_5120_R_20m/1 20 1 - contig5121 1 20 + 20 M 20
RESULTS = ((20, 1, "-", "contig5121", 1, 20, 20),
           (100,1, "-", "contig4209", 30120, 30219, 97),
           (17, 97, "+", "contig5043", 4873, 4953, 81),
           (30, 84, "+", "contig1732", 1, 55, 55),
           (1, 100, "+", "contig10737", 30839, 30938, 97)
           )

KMER = 11
NSKIP = 5

TMPFIL_PREFIX = "TMP"


def smalt_index(df, index_name, fasta_name, kmer, nskip):
    from sys import exit
    from subprocess import call

    tup = (PROGNAM, 'index',
           '-k', '%i' % (int(kmer)),
           '-s', '%i' % (int(nskip)),
           index_name,
           fasta_name)
    df.call(tup, "when indexing")

def smalt_map(df, oufilnam, indexnam, readfilnam):
    from sys import exit
    from subprocess import call
 
    tup = (PROGNAM, 'map',
           '-f', 'cigar',
           '-o', oufilnam,
           indexnam,
           readfilnam)
    df.call(tup, "when mapping")

def compare_result(cigfilnam):
    from formats import Cigar, openFile
    from sys import exit
    
    cig = Cigar()
    infil = openFile(cigfilnam)
    for result in RESULTS:
        cig.next(infil)
        if not cig.ok or \
               cig.qseg != (result[0], result[1]) or \
               cig.sense != result[2] or \
               cig.snam != result[3] or \
               cig.sseg != (result[4], result[5]) or \
               cig.swatscor != result[6]:
            exit("Unexpected result for read '%s'" % cig.qnam)
             
    infil.close()
    
def prep_fasta(filnam):
    from testdata import openFile

    infil = openFile(filnam, 'w')
    for seq in READSEQS:
        for i in range(len(seq)):
            infil.write("%c" % seq[i])
    infil.close()
    
if __name__ == '__main__':
    from testdata import DataFiles

    df = DataFiles()
    #refnam = df.unpack(REFSEQ)
    #refnam = df.addTMP(REFSEQ)
    refnam = df.joinData(REFSEQ)
    indexnam = df.addIndex(TMPFIL_PREFIX)
    readfilnam = df.addTMP(TMPFIL_PREFIX + ".fq")
    oufilnam = df.addTMP(TMPFIL_PREFIX + ".cig")
    prep_fasta(readfilnam)
    smalt_index(df,indexnam, refnam, KMER, NSKIP)
    smalt_map(df,oufilnam, indexnam, readfilnam)
    compare_result(oufilnam)
    df.cleanup()
