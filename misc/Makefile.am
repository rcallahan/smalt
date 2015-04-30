## Process this file with automake to produce Makefile.in
AM_CFLAGS = -I../src/
if ICC_COMPILER
if ICC_MMIC
AM_CFLAGS += -mmic	
endif
AM_CFLAGS += -Wall
if ICC_NOLIMITH
AM_CFLAGS += -D_GCC_LIMITS_H_
endif
else
AM_CFLAGS += -Wall -Wextra
endif

if APPLE_64BIT
AM_CFLAGS += -m64 -D_FILE_OFFSET_BITS=64
AM_CFLAGS += ${APPLECFLAGS}
endif

if WITH_ZLIB
ZLIB_LDFLAG= -lz
else
ZLIB_LDFLAG=
endif

bin_PROGRAMS = simread fetchseq splitreads mixreads \
	splitmates readstats trunkreads basqcol simqual

dist_data_DATA = SAM.py

simread_SOURCES = simread.c rsample.c rsample.h ../src/randef.h
simread_LDADD = ../src/libsgm.a ../src/libseq.a -lm $(ZLIB_LDFLAG)

fetchseq_SOURCES = fetchseq.c 
fetchseq_LDADD = ../src/libseq.a $(ZLIB_LDFLAG)

splitreads_SOURCES = splitreads.c
splitreads_LDADD = ../src/libseq.a $(ZLIB_LDFLAG)

mixreads_SOURCES = mixreads.c
mixreads_LDADD = ../src/libseq.a $(ZLIB_LDFLAG)

readstats_SOURCES = readstats.c
readstats_LDADD = ../src/libseq.a $(ZLIB_LDFLAG)

splitmates_SOURCES = splitmates.c
splitmates_LDADD = ../src/libseq.a $(ZLIB_LDFLAG)

trunkreads_SOURCES = trunkreads.c
trunkreads_LDADD = ../src/libseq.a $(ZLIB_LDFLAG)

basqcol_SOURCES = basqcol.c basqual.c basqual.h
basqcol_LDADD = ../src/libseq.a $(ZLIB_LDFLAG)

simqual_SOURCES = simqual.c basqual.c basqual.h ../src/diffstr.c ../src/diffstr.h \
	 	  ../src/randef.h
simqual_LDADD = ../src/libseq.a -lm $(ZLIB_LDFLAG)

