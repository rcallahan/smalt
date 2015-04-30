# parsing output formats

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

class Cigar:
    from re import compile
    DEBUG = False
    CIGARSTR = compile("^cigar:([ABCDNRSP]):(\d{2})\s+(\S+)\s+(\d+)\s+(\d+)\s+([*+-])\s+" \
                       "(\S+)\s+(\d+)\s+(\d+)\s+\+\s+(\d+)\s+(.+)$")

    def __init__(self):
        self._blank()

    def next(self, infil):
        okflg = False
        isEOF = False
        while not okflg:
            lin = infil.readline()
            if not lin:
                isEOF = True
                break
            self.parse(lin)
            okflg = self.ok
            if Cigar.DEBUG:
                print "DEBUG:Cigar.next %s" % (okflg)
        if Cigar.DEBUG:
            print "DEBUG:Cigar.next returns %s" % (isEOF)
        return isEOF
    
    def parse(self, lin):
        m = Cigar.CIGARSTR.match(lin)
        if m:
            self.lin = lin.strip()
            self.mapcls = m.group(1) # class label (A,B,C,D,N,R,S)
            self.mapq = int(m.group(2)) # mapping quality score
            self.qnam = m.group(3)   # name of query (read)
            self.qseg = (int(m.group(4)), int(m.group(5)))
            self.sense = m.group(6)
            self.snam = m.group(7)   # name of subject (reference)
            self.sseg = (int(m.group(8)), int(m.group(9)))
            self.swatscor = int(m.group(10))
            self.cigar = m.group(11)
            self.cigar.strip()
            self.ok = True
            if Cigar.DEBUG:
                print "DEBUG:Cigar.parse('%s')::'%s'" % (self.lin, self.cigar)
        else:
            lin.strip()
            if lin[0] != '#':
                print "not parsed: %s" % (lin)
            self._blank()

    def getMateNo(self):
        mateno = 0
        if self.qnam[-2:] == "/1":
            mateno = 1
        elif  self.qnam[-2:] == "/2":
            mateno = 2
        return mateno
    
    def _blank(self):
        self.mapcls = ""
        self.mapq = 0
        self.qnam = ""
        self.qseg = (0,0)
        self.sense = ''
        self.snam = ""
        self.sseg = (0,0)
        self.swatscor = 0
        self.lin = ""
        self.ok = False

    def __cmp__(self, other):
        rv = cmp(self.snam, other.snam)
        if rv: return rv
        
        rv = cmp(self.sseg[0], other.sseg[0])
        if rv: return rv

        return cmp(other.sseg[1], self.sseg[1])
                
def getNextCigarPair(infil, cigA, cigB, mateno_check = True):
    isOk = False
    isEOF = cigA.next(infil)

    if not isEOF:
        if cigB.next(infil):
            exit("ERROR EOF when reading 2nd mate ...")
    
        if not cigA.ok:
            exit("ERROR when parsing 1st mate")
            
        if not cigB.ok:
            exit("ERROR: when parsing 2nd mate")
        isOk = True
    if not isEOF and mateno_check:
        if cigA.qnam[-2] != cigB.qnam[-2]:
            print "ERROR: read names don't match %s, %s" % (cigA.qnam, cigB.qnam)
            isOk = False
        elif cigA.getMateNo() != 1 or cigB.getMateNo() != 2:
            print "ERROR in mate number %s [1] and %s [2]" % (cigA.qnam, cigB.qnam)
            isOk = False
    
    return (isOk, isEOF)

def openFile(filnam, mode = 'r'):
    import gzip

    is_compressed = len(filnam) > 3 and filnam[-3:] == ".gz"
    try:
        if is_compressed:
            fil = gzip.open(filnam, mode)
        else:
            fil = open(filnam, mode)
    except:
        exit("ERROR when opening file '%s'" % filnam)

    return fil

if __name__ == "__main__":
    from sys import argv, exit

    if len(argv) < 2:
        exit("usage: %s <cigar file>" % argv[0])

    cig = Cigar()
    
    infil = openFile(argv[1])

    while not cig.next(infil):
        print cig.qnam
    
    infil.close()
