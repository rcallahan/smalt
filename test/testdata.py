# Handling of data files for tests

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

class DataFiles:
    INPUTDIR = "./data"
    WORKDIR = "./tmp"
    LOGFILNAM = "LOGFILE"
    def __init__(self, df=[]):
        from os import access, F_OK, mkdir, path
        if not access(DataFiles.WORKDIR, F_OK):
            mkdir(DataFiles.WORKDIR)
        logfilnam = path.join(DataFiles.WORKDIR, DataFiles.LOGFILNAM)
        self.tmpfiles = [logfilnam + ".out", logfilnam + ".err"]
        self.datafiles = []
        self.logfil_out = openFile(self.tmpfiles[0], 'w')
        self.logfil_err = openFile(self.tmpfiles[1], 'w')
        
        self.addPackedSourceFiles(df)
        
    def __delete_(self):
        self.logfil_out.flush()
        self.logfil_out.close()
        self.logfil_err.flush()
        self.logfil_err.close()
        
    def joinData(self, filnam):
        from os import path
        return path.join(DataFiles.INPUTDIR, filnam)
    
    def addPackedSourceFiles(self, filnams):
        for fn in filnams:
            self.unpack(fn)
        
    def unpack(self, filnam):
        from os import path

        fn_from = path.join(DataFiles.INPUTDIR, filnam)
        if len(filnam) < 4 or filnam[-3:] != '.gz': fn_from = fn_from + '.gz'
        fn_to = path.join(DataFiles.WORKDIR, filnam)
        infil = openFile(fn_from, 'r')
        oufil = openFile(fn_to, 'w')
        self.logfil_out.write("unpacking '%s' -> '%s'\n" % (fn_from, fn_to))
        while 1:
            lin = infil.readline()
            if not lin:
                break
            oufil.write(lin)
        oufil.close()
        infil.close()
        self.datafiles.append(fn_to)

        return fn_to

    def addIndex(self, filnam):
        from os import path
        fn = path.join(DataFiles.WORKDIR, filnam)
        for f in (fn + '.smi', fn + '.sma'):
            if f not in self.tmpfiles:
                self.tmpfiles.append(f)
        return fn
    
    def addTMP(self, filnam):
        from os import path
        fn = path.join(DataFiles.WORKDIR, filnam)
        if fn not in self.tmpfiles:
            self.tmpfiles.append(fn)
        return fn

    def call(self, tup, errmsg="", oufilnam=None):
        from subprocess import call
    
        if oufilnam:
            oufil = openFile(oufilnam, 'w')
        else:
            oufil = self.logfil_out
            
        self.addLog(tup)
        rc = call(tup, stdout=oufil, stderr=self.logfil_err)
        
        if oufilnam:
            oufil.close()

        if rc != 0:
            if not errmsg:
                errmsg = "when calling %s" % (tup[0])
            self.exitErr("ERROR %s: returned with code %i" % (errmsg, rc))
        self.logfil_out.flush()
        self.logfil_err.flush()
                
    def addLog(self, tup):
        if not tup:
            return
        self.logfil_out.write("%s" % tup[0])
        if len(tup) > 1:
            for t in tup[1:]:
                self.logfil_out.write(" %s" % t)
        self.logfil_out.write("\n")
        self.logfil_err.flush()
        return
    
    def addErr(self, msg):
        self.logfil_err.write(msg + "\n")
        self.logfil_out.flush()
        self.logfil_err.flush()

    def exitErr(self, msg):
        from sys import exit
        self.addErr(msg)
        exit(msg)
        
    def cleanup(self, is_verbose=False):
        from os import access, F_OK, remove, removedirs
        self.logfil_out.close()
        self.logfil_err.close()
        for filnam in (self.datafiles + self.tmpfiles):
            if access(filnam, F_OK):
                remove(filnam)
                if is_verbose: print "Removed file '%s'" % filnam
            else:
                print "Could not access file '%s'" % filnam
        self.datafiles = []
        self.tmpfiles = []
        try:
            removedirs(DataFiles.WORKDIR)
        except:
            "Could not remove temporary directory %s" % DataFiles.WORKDIR

def openFile(filnam, mode = 'r'):
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

def areFilesIdentical(filnamA, filnamB, ignore_key=""):
    # compare text files line by line
    # ignore lines that start with ignore_key
    infilA = openFile(filnamA)
    infilB = openFile(filnamB)
    kl = len(ignore_key)
    okflg = True
    while okflg:
        linA = infilA.readline()
        linB = infilB.readline()

        if not linA or not linB:
            okflg = not linA and not linB
            break
        
        if linA != linB and \
               (kl < 1 or \
                linA[:kl] != linB[:kl] or \
                linA[:kl] != ignore_key):
            okflg = False
        
    return okflg

