# SAM line parser

class Sam:
    import re
    SAMSTR = re.compile("^(\S+)\t(\d+)\t(\S+)\t(\d+)\t(\d+)\t(\S+)\t(\S+)\t"\
                        "(\d+)\t([\+\-]*\d+)\t(\S+)\t(\S+)"\
                        "(\t.+)?")
    TAGSTR = re.compile("^(\S{2}):([iZ]):(\S+)")
#    SAMSTR = re.compile("^(\S+)\s+(\d+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+"\
#                        "(\d+)\s+([\+\-]*\d+)\s+(\S+)\s+(\S+)(\s+AS:i:\d+){0,1}"\
#                        ".*(\s+NM:i:\d+){0,1}")

    QNAMSTR = re.compile("^(\S+)\/([12])$")
    ALISCOR = re.compile("\tAS:i:(\d+)")
    CLIPSTART = re.compile("^(\d+)([HS])")
    CLIPEND   = re.compile("(\d+)([HS])$")
    QEXT = re.compile("/([12])$")
    
    FLAG_PAIRED     = 0x0001
    FLAG_PROPER     = 0x0002
    FLAG_NOMAP      = 0x0004
    FLAG_MATENOMAP  = 0x0008
    FLAG_STRAND     = 0x0010 # reverse complemented
    FLAG_MATESTRAND = 0x0020
    FLAG_1stMATE    = 0x0040
    FLAG_2ndMATE    = 0x0080
    FLAG_NOTPRIMARY = 0x0100
    FLAG_CHECKFAIL  = 0x0200
    FLAG_DUPLICATE  = 0x0400

    DEFAULT_SWSCOR = 0
    MAPQ_NONRANDOM = 4

    FASTQ_FORMAT = '@%s\n%s\n+\n%s\n'

    COMPLEMENTOR = {'A':'T', 'C':'G', 'G':'C', 'T':'A',
                    'a':'t', 'c':'g', 'g':'c', 't':'a',
                    'N':'N', 'n':'n'}

    def __init__(self):
        self.blank()
        
    def blank(self):
        self.ok = False
        self.qname = ""
        self.flag = 0
        self.rname = ""
        self.pos = 0
        self.mapq = 0
        self.cigar = ""
        self.mrnm = ""
        self.mpos = 0
        self.isize = 0
        self.swscor = 0
        self.nt = "" # nucleotide sequence
        self.qual = "" # base qualities
        self.tags = {}
        
    def parse(self, lin, is_verbose=True):
        self.blank()
        m = Sam.SAMSTR.match(lin)
        if m:
            #print m.groups()
            #print "GROUPS: %i" % len(m.groups())
            self.qname = m.group(1)
            self.flag = int(m.group(2))
            self.rname = m.group(3)
            self.pos = int(m.group(4))
            self.mapq = int(m.group(5))
            self.cigar = m.group(6)
            self.mrnm = m.group(7)
            self.mpos = int(m.group(8))
            self.isize = int(m.group(9))
            self.nt = m.group(10)
            self.qual = m.group(11)
            self.tags = {}

            if len(m.groups()) > 11:
                tags = m.group(12).split('\t')
                #print "TAGS", tags
                for tag in tags:
                    n = Sam.TAGSTR.match(tag)
                    if n:
                        (tagnam, typ, fld) = n.groups()
                        self.tags[tagnam] = (typ, fld)
            self.ok = True
        else:
            if is_verbose: print "NOT PARSED '%s'" % lin
            self.blank()

    def clip(self):
        s = 0
        e = 0
        typ = 'H'
        isCorrect = False
        ms = Sam.CLIPSTART.match(self.cigar)
        if ms:
            s = int(ms.group(1))
            typ = ms.group(2)
        
        me = Sam.CLIPSTART.search(self.cigar)
        if me:
            e = int(me.group(1))
            if ms and ms.group(2) == me.group(2):
                isCorrect = True
        return (typ, isCorrect, s, e)

    def strand(self):
        print "(%i, %s)" % (self.flag, bool(self.flag & Sam.FLAG_STRAND))
        return bool(self.flag & Sam.FLAG_STRAND)
    
    def calcUnclippedStart(self):
        rs = self.pos
        (typ, ok, s, e) = self.clip()
        #print rs
        #print (typ, ok, s, e)
        if typ == 'H' and ok:
            rs = rs - s
        return rs

    def cmpFlags(self, other):
        rv = True
        if self.flag != other.flag:
            f = Sam.FLAG_PROPER | Sam.FLAG_MATENOMAP | Sam.FLAG_MATESTRAND
            if (self.flag & (~f)) != (other.flag & (~f)):
                rv = False
        return rv

    def strFlags(self):
        s = "FLAG:\n"
        if self.flag & Sam.FLAG_PAIRED:
            i = 1
        else:
            i = 0
        s = s + '[%1i] Paired Read\n' % i
        if self.flag & Sam.FLAG_PROPER:
            i = 1
        else:
            i = 0
        s = s + '[%1i] Proper Read\n' % i
        if self.flag & Sam.FLAG_1stMATE:
            i = 1
        else:
            i = 0
        s = s + '[%1i] 1st mate of a pair\n' % i
        if self.flag & Sam.FLAG_2ndMATE:
            i = 1
        else:
            i = 0
        s = s + '[%1i] 2nd mate of a pair\n' % i
        return s
        
    def isMapped(self):
        return  self.name != '*'

    def getMateNam(self):
        #print "getMateNam:%s" % (self.qname)
        m = Sam.QNAMSTR.search(self.qname)
        if m:
            rv = (m.group(1), int(m.group(2)))
        else:
            flg = int(self.flag)
            if flg & Sam.FLAG_1stMATE:
                n = 1
            elif flg & Sam.FLAG_2ndMATE:
                n = 2
            else:
                n = 0
            rv = (self.qname, n)
        return rv
    
    def compare(self, other):
        rv = True
        s = ""
##         if self.cigar != other.cigar:
##             print "CIGAR strings differ! %s::%s" % \
##                   (self.cigar, other.cigar)
##             rv = False
        spos = self.calcUnclippedStart()
        opos = other.calcUnclippedStart()
        if spos != opos:
            if self.cigar == '*':
                s = s + "%s not mapped" % self.qname
                rv = False
            elif other.cigar == '*':
                s = s + "%s not mapped" % other.qname
                rv = False
            else:
                s = s + "Mapping positions differ! (%i:%i, %i:%i)" % \
                    (self.pos, spos, other.pos, opos)
                rv = False
        elif not self.cmpFlags(other):
            s = s + "Flags differ"
            rv = False
        return (rv, s)

    def asFastqStr(self):
        mateno = 0
        seq =""
        qual=""
        if self.flag & Sam.FLAG_1stMATE:
            namstr = '%s/1' % self.qname
        elif self.flag & Sam.FLAG_2ndMATE:
            namstr = '%s/2' % self.qname
        else:
            namstr = self.qname
        if  bool(self.flag & Sam.FLAG_STRAND):
            (seq, qual, okflg) = self.reverseComplement()
            if not okflg:
              print "ERROR: when reverse complementing sequence!"
              exit(1)
        else:
            seq = self.nt
            qual = self.qual

        return namstr, seq, qual
    
    def asFastq(self, clip=0):
        namstr, seq, qual = self.asFastqStr()
        if clip >= len(self.nt):
            print "ERROR: clipped fragment > sequence length\n%s" % (self.qname)
        return Sam.FASTQ_FORMAT % (namstr, seq[clip:], qual[clip:])

    def reverseComplement(self):
        s = list(self.qual)
        s.reverse()
        qual = ""
        for i in range(len(self.qual)):
            qual = qual + s[i]
    
        s = list(self.nt)
        s.reverse()
        seq = ""
        for i in range(len(self.nt)):
            try:
                seq = seq + Sam.COMPLEMENTOR[s[i]]
            except:
                seq = seq + 'N'

        return (seq, qual, len(self.nt) == len(self.qual))

    def next(self, infil, is_verbose=False):
        while 1:
            lin = infil.readline()
            if not lin: break
            self.parse(lin, is_verbose)
            if self.ok: break
            if is_verbose: print "not parsed: %s" % lin
        return not lin # True if EOF
    
        
def fetchNextRead(infil, read):
    isEOF = False
    currnam = read.qname
    while 1:
        lin = infil.readline()
        if not lin:
            isEOF = True
            break
        read.parse(lin)
        if not read.ok:
            print "not parsed: %s" % lin
            continue
        break
        #print read.target.varnum()
        #if read.qname != currnam:
        #    break

    return isEOF

def fetchNextPair(infil, samA, samB):
    # assumes sorted file
    from sys import exit
    switch_flag = 0
    isEOF = fetchNextRead(infil, samA)
    if not isEOF:
        isEOF = fetchNextRead(infil, samB)
        if not isEOF:
            errflg = False
            if (samA.flag & Sam.FLAG_PAIRED) == 0 or (samB.flag & Sam.FLAG_PAIRED) == 0:
                errflg = True
                
            if (samA.flag & Sam.FLAG_1stMATE) == 0:
                switch_flag = 1
                
            if switch_flag == 0:
                if (samB.flag & Sam.FLAG_2ndMATE) == 0:
                    errflg = True
            else:
                if (samB.flag & Sam.FLAG_2ndMATE) == 1:
                    errflg = True
  
            if errflg:
                print "ERROR: unexpected pair flags.\n%s\n%s\n%s\n%s" % \
                      (samA.qname, samA.strFlags(), samB.qname, samB.strFlags())
                exit(1)

    return (isEOF, switch_flag)

def openFile(filnam, mode):
    import gzip

    is_compressed = len(filnam) > 3 and filnam[-3:] == ".gz"
    try:
        if is_compressed:
            oufil = gzip.open(filnam, mode)
        else:
            oufil = open(filnam, mode)
    except:
        print "ERROR when opening file '%s'" % filnam
        exit(1)

    return oufil

if __name__ == '__main__':
    from sys import argv, exit

    if len(argv) < 3:
        print "usage: %s <SAM file (input)> <mapping score threshold>" % argv[0]
        exit(1)

    infilnam = argv[1]
    mapq_min = int(argv[2])

    infil = openFile(infilnam, 'r')
    
    read = Sam()
    readctr = 0
    pairctr = 0
    nomapctr = 0
    lowqctr = 0
    minqctr = 0
    highqctr = 0
    propctr = 0
    chimictr = 0 # counter for chimeric inserts
    old_qnam = ""
    old_rnam = ""
    isizarr = []
    while 1:
        if fetchNextRead(infil, read):
            break;
        flg = int(read.flag)
        readctr = readctr + 1
        if (flg & Sam.FLAG_NOMAP):
            nomapctr = nomapctr + 1
        else:
            if int(read.mapq) >= Sam.MAPQ_NONRANDOM:
                minqctr = minqctr + 1
            if int(read.mapq) >= mapq_min:
                highqctr = highqctr + 1
                if (flg & (Sam.FLAG_PAIRED | Sam.FLAG_NOMAP | Sam.FLAG_MATENOMAP)) == Sam.FLAG_PAIRED:
                    pairctr = pairctr + 1
                    if (flg & Sam.FLAG_PROPER) == 0:
                        # not a proper pair
                        (mnam, mno) = read.getMateNam()
                        if mno == 1: #print mnam
                            old_qnam = mnam
                            old_rnam = read.rname
                        else:
                            if mnam == old_qnam:
                                if read.rname != old_rnam:
                                    chimictr = chimictr + 1
                                    print "[%i, %i]%s: rname1: %s, rname2: %s" % \
                                          (readctr, chimictr, mnam, old_rnam, read.rname)

                    else:
                        # a 'proper' pair
                        print "%s %i %i %x" % (read.qname, abs(read.isize), read.mapq, int(read.flag))
                        propctr =  propctr + 1
                isizarr.append(abs(read.isize))
            else:
                lowqctr = lowqctr + 1
        if readctr % 100000 == 0:
            print "%i reads ..." % readctr
    pairctr = pairctr/2
    
    infil.close()

    print "%i out of a total of %i reads (%5.2f%%) were mapped." % \
          (readctr - nomapctr, readctr, 100*float(readctr - nomapctr)/readctr)
    print "%i out of a total of %i reads (%5.2f%%) mapped with a mapping score > %i" % \
          (minqctr, readctr, 100*float(minqctr)/readctr, Sam.MAPQ_NONRANDOM)
    if Sam.MAPQ_NONRANDOM != mapq_min:
        print "%i out of a total of %i reads (%5.2f%%) mapped with a mapping score > %i" % \
              (highqctr, readctr, 100*float(highqctr)/readctr, mapq_min)
        
    if pairctr > 0:
        print "%i out of a total of %i reads (%5.2f%%) mapped as a proper pair with a mapping score > %i" % \
              (propctr, readctr, 100*float(propctr)/readctr, mapq_min)
        print "%i of %i pairs (%6.3f%%) with a mapping score >= %i mapped to different chromosomes" % \
              (chimictr, pairctr, 200*float(chimictr)/pairctr, mapq_min)
    else:
        print "There were no reads mapped as pairs."
        
    exit(0)
