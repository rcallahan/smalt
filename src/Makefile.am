## Process this file with automake to produce Makefile.in
##
## Mac 64bit binaries: compile with 
## CFLAGS="-Wall -g -O2 -pthread -m64 -D_FILE_OFFSET_BITS=64"
## CXXFLAGS="-m64 -D_FILE_OFFSET_BITS=64"

if ICC_COMPILER
AM_CFLAGS = -Wall
if ICC_MMIC
AM_CFLAGS += -mmic	
endif
else
AM_CFLAGS = -Wall -Wextra
endif

if ICC_NOLIMITH
AM_CFLAGS += -D_GCC_LIMITS_H_
endif

if APPLE_64BIT
AM_CFLAGS += -m64 -D_FILE_OFFSET_BITS=64
APPLELDFLAGS = -pthread 
endif

ALI_SRC = score.c score.h diffstr.c diffstr.h \
	  alibuffer.c alibuffer.h alibuffer_struct.h \
	  alignment.c alignment.h
if WITH_SIMD
ALI_SRC += swsimd.c swsimd.h	
endif

SGM_SRC = sort.c sort.h array.c array.h vdef.h interval.c interval.h \
	  hashidx.c hashidx.h hashhit.c hashhit.h \
          segment.c segment.h 

MAP_SRC = results.c results.h randef.h \
	  resultpairs.c resultpairs.h \
	  report.c report.h \
	  rmap.c rmap.h insert.c insert.h

bin_PROGRAMS = smalt

noinst_LIBRARIES = libseq.a libsgm.a

libseq_a_SOURCES = elib.c elib.h filio.c filio.h sequence.c sequence.h

libsgm_a_SOURCES = ${ALI_SRC} ${SGM_SRC}
libsgm_a_CFLAGS = ${AM_CFLAGS}
if WITH_SIMD
if !ICC_COMPILER
libsgm_a_CFLAGS += -msse2
endif
endif
libsgm_a_LIBADD = libseq.a

smalt_SOURCES = smalt.c menu.h menu.c threads.h threads.c infmt.c infmt.h ${MAP_SRC}
smalt_CFLAGS = $(BAMBAMC_INCL) ${AM_CFLAGS}
smalt_LDFLAGS = $(BAMBAMC_LDFLG) ${APPLELDFLAGS} -static 
smalt_LDADD = libsgm.a libseq.a $(BAMBAMC_DEPS) -lm -lpthread #-lrt
if WITH_ZLIB
smalt_LDADD += -lz	
endif
