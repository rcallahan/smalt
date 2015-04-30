# test cigar string in SAM/BAM output formats

REFSEQ = ("acaaaaaattaaataaaaatattaagaaaaagagaaattacacaaactattacatactat" \
          "aatgcatttttttctttagtgtttatgttattaaaatatatttttttcataatatatatt" \
          "aaatcacgtcatatactattttacaatttatcttatattaattgtattattacttgttct" \
          "aaaattataattctaaatatgttattttaataattatattataataattgctattataat" \
          "aaattacatattatattaaacaagtaacacgatattatttatatactataaaaatacaag" \
          "gaatcgtttatttttgtattttatacgtaattagaatatttaaaaaaaaataataatgga" \
          "atacacatattacattatatgttatatttacatataacatagtaatgcatattttatagt" \
          "aaattagtttgcaaaaccttataaataataataaatataatttaaataatcatcttatac" \
          "ttaataagcaataataaaatccaatcatatataaacttaagcaaggaaatttaaatgagg",
          
          "atagtaaattaacccatttattaaaaaattctcttgaaggcaattgtctagttgtaatga" \
          "tcgcaaatataaacccttctagaacatcctttcaagaatctaataatactcttaaatacg")

READSEQ = [
    ("aaatcacgtcatatactattttacaatttatcttatattaattgtattattacttgttct",
     ("60M", "60M"),
     "NM:i:0"
     ),
    ("acaattataattctaaatatgttatttaataattatattataataattgctattattat",
     ("2S25M1D29M3S", "2S25M1D29M3S"),
            "NM:i:1"
     ),
    ("atatacatattacattatatgttatatttacatatggaacatagcaatgcatattttatagt",
     ("35M2I25M", "3M1X31M2I7M1X17M"),
     "NM:i:4"
     ),
    ("TCCATGATTATTTTTTTTAAATATTCTAATTACGTATAAAAATACAAACATAAACGATTC",
     ("22M1I27M1D10M", "11M1X10M1I27M1D4M1X5M"),
     "NM:i:4"
     )
    ]

READSEQ_PAIR = [
    ("aatgcatttttttctttagtgtttatgttattaaaatatatttttttcataatatatatt",
     ("60M", "60M"),
     "NM:i:0",
     ),
    ("ACTATAAAATGTGCATTACTATGTTATATGTAAATATAACATATAATGTAATATGTGTAT",
     ("60M", "49M1X10M"),
     "NM:i:1",
     )
    ]

PROGNAM = "../src/smalt"
SAMTOOLS_EXEC = "samtools"

# look for SAM.py module in the following directory
from sys import path
path.append('../misc')

KMER = 7
NSKIP = 1

TMPFIL_PREFIX = "TMP"
TMPFIL_PREFIX_PAIRED = TMPFIL_PREFIX + "pair"

REFNAM_PREFIX = "REF_"
READNAM_PREFIX = "READ_"

SAM_TEST_FIELDS = {0: "read name",
                   1: "FLAG field",
                   2: "chromosome (read)",
                   3: "position (read)",
                   4: "mapping quality",
                   5: "cigar string"
                   }
                   
def writeFASTAref(filnam):
    from testdata import openFile
    oufil = openFile(filnam, 'w')
    sn = 0
    for seq in REFSEQ:
        sn = sn + 1
        oufil.write(">%s%i\n%s\n" % (REFNAM_PREFIX, sn, seq))
    oufil.close()

def writeFASTAreads(seqdat, filnam):
    from testdata import openFile
    oufil = openFile(filnam, 'w')
    rctr = 0
    for rdat in seqdat:
        rctr = rctr + 1
        oufil.write(">%s%i\n%s\n" % (READNAM_PREFIX, rctr, rdat[0]))
    oufil.close()

def writeFASTAreadPairs(seqdat, filnamA, filnamB):
    from testdata import openFile
    oufilA = openFile(filnamA, 'w')
    oufilB = openFile(filnamB, 'w')
    rctr = 0
    matectr = 1
    for rdat in seqdat:
        if matectr == 1:
            oufil = oufilA
        else:
            oufil = oufilB
        oufil.write(">READ_%i/%1i\n%s\n" % (rctr, matectr, rdat[0]))
        if matectr > 1:
            matectr = 1
            rctr = rctr + 1
        else:
            matectr = 2

    oufilB.close()
    oufilA.close()

def smalt_index(df, index_name, fasta_name, kmer, nskip):
    tup = (PROGNAM, 'index',
           '-k', '%i' % (int(kmer)),
           '-s', '%i' % (int(nskip)),
           index_name,
           fasta_name)
    df.call(tup, "when indexing")

def smalt_map(df, oufilnam, indexnam, readfil, format="sam", matefil=""): 
    tup = [PROGNAM, 'map']
    
    tup.extend(
        ['-f', format,
         '-o', oufilnam,
         indexnam,
         readfil]
        )
    if matefil and len(matefil) > 0:
        tup.append(matefil)
        
    #print tup
    df.call(tup, "when mapping")

def samtools_bam2sam(infilnam, oufilnam):
    from sys import exit
    from subprocess import call

    tup = [SAMTOOLS_EXEC, 'view', '-h', '-o', oufilnam, infilnam]
    
    rv = call(tup)
    if rv:
        exit("ERROR when converting bam to sam: exited with code %i" % (rv))  

def testSAMfilesAreIdentical(filnamA, filnamB):
    from testdata import openFile

    infilA = openFile(filnamA)
    infilB = openFile(filnamB)

    while 1:
        linA = infilA.readline()
        if not linA or linA[0]!="@":
            break

    while 1:
        linB = infilB.readline()
        if not linB or linB[0]!="@":
            break
        
    okflg = False
    
    while linA and linB:
        fldA = linA.split('\t')
        fldB = linB.split('\t')

        okflg = len(fldA) == len(fldB)

        if not okflg:
            print "lines differ in the number of fields:\n%s\n%s\n" % \
                  (linA.strip(), linB.strip())
            break
        
        for i in SAM_TEST_FIELDS.keys():
            okflg = fldA[i] == fldB[i]
            if not okflg and i == 0:
                # samtools-0.1.18 view -h produces non-printing char
                # directly after header
                okflg = fldA[i][1:] == fldB[i] or fldA[i] == fldB[i][1:]
            if not okflg:
                print "lines differ in %s:\n%s\n%s\n" % \
                      (SAM_TEST_FIELDS[i], linA.strip(), linB.strip())
                break       
        if not okflg:
            break
        linA = infilA.readline()
        linB = infilB.readline()

    if okflg:
        okflg = not linA and not linB
    
    return okflg

def checkSAMfile(df, readdat,samfilnam, xcig=False):
    from SAM import Sam, openFile

    sam = Sam()
    infil = openFile(samfilnam, 'r')

    linctr = 0

    for rdat in readdat:
        if sam.next(infil): break
        if xcig:
            cigar = rdat[1][1]
        else:
            cigar = rdat[1][0]
        if sam.cigar != cigar:
            df.exitErr("Unexpected CIGAR string. Got '%s' - expected '%s'" % (sam.cigar, cigar))
        nmtag = sam.tags["NM"]
        if nmtag:
            tagstr = "NM:%s:%s" % nmtag
            if tagstr != rdat[2]:
                df.exitErr("Unexpected tag '%s' (expected '%s')" % (tagstr, rdat[2]))
        linctr = linctr + 1

    if linctr < len(readdat):
        df.exitErr("CIGAR strings incomplete");
        
    infil.close()
    
if __name__ == '__main__':
    from testdata import DataFiles

    df = DataFiles()
    
    reffilnam = df.addTMP(TMPFIL_PREFIX + "ref.fa")
    writeFASTAref(reffilnam)

    
    indexnam = df.addIndex(TMPFIL_PREFIX)
    smalt_index(df, indexnam, reffilnam, KMER, NSKIP)

    readfilnam = df.addTMP(TMPFIL_PREFIX + "read.fa")
    mate1filnam = df.addTMP(TMPFIL_PREFIX_PAIRED + "mate1.fa")
    mate2filnam = df.addTMP(TMPFIL_PREFIX_PAIRED + "mate2.fa")
    writeFASTAreads(READSEQ, readfilnam)
    writeFASTAreadPairs(READSEQ_PAIR, mate1filnam, mate2filnam)
    
    samoufilnam = df.addTMP(TMPFIL_PREFIX + "out.sam")
    bamoufilnam = df.addTMP(TMPFIL_PREFIX + "out.bam")
    sambamnam = df.addTMP(TMPFIL_PREFIX + "out.bam.sam")

    samoufilnam_x = df.addTMP(TMPFIL_PREFIX + "outx.sam")
    bamoufilnam_x = df.addTMP(TMPFIL_PREFIX + "outx.bam")
    sambamnam_x = df.addTMP(TMPFIL_PREFIX + "outx.bam.sam")
   
    # the same for paired reads
    sam_paired_oufilnam = df.addTMP(TMPFIL_PREFIX_PAIRED + ".out.sam")
    bam_paired_oufilnam = df.addTMP(TMPFIL_PREFIX_PAIRED + ".out.bam")
    sambam_paired_nam = df.addTMP(TMPFIL_PREFIX_PAIRED + "out.bam.sam")

    # cigar strings without X for mismatch, single reads
    smalt_map(df, samoufilnam, indexnam, readfilnam, "sam")
    smalt_map(df, bamoufilnam, indexnam, readfilnam, "bam")
    samtools_bam2sam(bamoufilnam, sambamnam)
    isOK = testSAMfilesAreIdentical(sambamnam, samoufilnam)
    if not isOK:
        exit("ERROR: SAM and BAM files differ!")
    checkSAMfile(df, READSEQ, samoufilnam, xcig=False)

    # cigar strings without X for mismatch, paired reads
    smalt_map(df, sam_paired_oufilnam, indexnam, mate1filnam, "sam", mate2filnam)
    smalt_map(df, bam_paired_oufilnam, indexnam, mate1filnam, "bam", mate2filnam)
    samtools_bam2sam(bam_paired_oufilnam, sambam_paired_nam)
    isOK = testSAMfilesAreIdentical(sambam_paired_nam, sam_paired_oufilnam)
    if not isOK:
        exit("ERROR: SAM and BAM files differ for paired reads!")
    checkSAMfile(df, READSEQ_PAIR,sam_paired_oufilnam, xcig=False)

    # cigar strings with X for mismatch
    smalt_map(df, samoufilnam_x, indexnam, readfilnam, "sam:x")
    smalt_map(df, bamoufilnam_x, indexnam, readfilnam, "bam:x")
    samtools_bam2sam(bamoufilnam_x, sambamnam_x)
    isOK = testSAMfilesAreIdentical(sambamnam_x, samoufilnam_x)
    if not isOK:
        exit("ERROR: SAM and BAM files differ (CIGAR strings with X)!")
    checkSAMfile(df, READSEQ, samoufilnam_x, xcig=True)

    # cigar strings with X for mismatch, paired reads
    smalt_map(df, sam_paired_oufilnam, indexnam, mate1filnam, "sam:x", mate2filnam)
    smalt_map(df, bam_paired_oufilnam, indexnam, mate1filnam, "bam:x", mate2filnam)
    samtools_bam2sam(bam_paired_oufilnam, sambam_paired_nam)
    isOK = testSAMfilesAreIdentical(sambam_paired_nam, sam_paired_oufilnam)
    if not isOK:
        exit("ERROR: SAM and BAM files differ for paired reads!")
    checkSAMfile(df, READSEQ_PAIR,sam_paired_oufilnam, xcig=True)

    #print "Test ok."  
    df.cleanup()
    exit()
