#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
#
# ############### Structure of this file ###############
#   Autoconf requirements
#   AC_INIT(package, version, bug-report-address)
#   information on the package
#   checks for programs
#   checks for libraries
#   checks for header files
#   checks for types
#   checks for structures
#   checks for compiler characteristics
#   checks for library functions
#   checks for system services
#   AC_CONFIG_FILES([file...])
#   AC_OUTPUT
# ######################################################

AC_PREREQ([2.65])
AC_INIT([smalt],[0.7.6],[hp3@sanger.ac.uk],[smalt],[http://smalt.sourceforge.net])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_CONFIG_SRCDIR([src/elib.c])
AC_CONFIG_HEADERS([config.h:config.hin])

# Checks for programs.
AC_PROG_CC([icc gcc cc])
#AC_PROG_CC

# C++ compiler used for linking with c++ library bambam
#AC_PROG_CXX

# check wether compiler supports '-c -o' together
# this macro is obsolete with autoconf 1.14 
# AM_PROG_CC_C_O

# we use a archiver, so add support for unusual archivers
# AM_PROG_AR

# find pkg-config tests
PKG_PROG_PKG_CONFIG

# set CFLAGS defaults explicitly
#AC_SUBST([AM_CFLAGS],['-O0 -g -Wall'])
#AC_SUBST([CFLAGS],['-O0 -g -Wall'])

# set LDFLAGS defaults explicitly
#AC_SUBST([LDFLAGS],[-static])

AC_PROG_RANLIB
AM_PATH_PYTHON([2.6],[with_testdriver=true],
	[with_testdriver=false;
	AC_MSG_WARN("Could not find python version >= 2.6. Won't run tests.")])

# Checks for programs

# samtools executable is required for some tests
AC_CHECK_PROG([samtoolsexec], [samtools], [yes])

# check if compiling for Intel's Many Integrated Core Architecture (MIC)
AC_ARG_ENABLE([mic], 
              [AS_HELP_STRING([--enable-mic], 
	      [compile for Many Integrated Core Architecture])],
	      AS_IF([test $CC = "icc"], [for_intel_mic=$enableval], [for_intel_mic=no]),
	      [for_intel_mic=no])
#AS_ECHO($for_intel_mic)
# Checks for libraries.

# optional: samtools built as library libbam
AC_ARG_WITH([samtools-dir],
            [AS_HELP_STRING([--with-samtools-dir@<:@=DIR@:>@],
		[search for samtools in DIR])], 
            [],
            [with_samtools=check])
AS_IF([test "x$with_samtools" != xno],
      [AS_IF([test "x$with_samtools" == xcheck],
             [AC_CHECK_LIB([bam],[sam_open],
			   [samtools=yes;],
			   [samtools=no;], -lz)
	     ],
             [ldflg="$LDFLAGS";
	     samtools_ldflg="$withval";
	     samtools_incl="$withval";
	     AC_SUBST([LDFLAGS], ["$LDFLAGS -L$samtools_ldflg"])     
       	     AC_CHECK_LIB([bam],[sam_open],
			  [samtools=yes;
			  AC_SUBST([SAMTOOLS_DEPS], ["-lbam -lz"])
	  		  AC_SUBST([SAMTOOLS_CFLAGS], ["-I$samtools_incl"])
	  		  AC_SUBST([SAMTOOLS_LDFLAGS], ["-L$samtools_ldflg"])
			  ],
      			  [samtools=no;
            		  AC_MSG_FAILURE(
			  [--with-samtools-prefix was specified, but test for samtools library failed])
       	        	  ], -lz)
             AC_SUBST([LDFLAGS], ["$ldflg"])
       	     ]);
        if test "$xsamtool"s == xyes; then
	     AC_DEFINE([HAVE_SAMTOOLS_API], [1],
                       [Define if you have samtools installed as a library])
       	fi])	    
 	    
# optional: German Tischler's bambam library
AC_ARG_WITH([bambamc],
            [AS_HELP_STRING([--with-bambamc],
	      [support BAM format @<:@default=check@:>@])],
            [],
            [with_bambamc=check])
AS_IF([test "x$with_bambamc" != xno],
      [PKG_CHECK_MODULES([BAMBAMC], [bambamc >= 0.0.42], 
         [ bambamc=yes; ],
	 [bambamc=no;
	  if test "x$with_bambamc" = xyes; then
	       if  test -z $BAMBAMC_CFLAGS; then
             	   AC_MSG_FAILURE(
				[bambamc library was not found])
	       else
	          bambamc=yes;
		  AC_MSG_WARN(
			[library not known to pkg-config, will use linker flags])	
	       fi
          else 
             AC_MSG_WARN(
	        [library not found])	
          fi
         ])
	 if test x$bambamc = xyes; then
	    AC_DEFINE([HAVE_BAMBAMC], [1],
                      [Define if you have the package bambamc installed as a library])
 	 fi
	 bbincludes=
	 bblibdirs=
	 bblibdeps=
	 for elem in $BAMBAMC_CFLAGS
	 do
    	   if test ${elem::2} = "-I" || test ${elem::2} = "-f"; then 
	     bbincludes="$bbincludes $elem"
           fi
         done
	 for elem in $BAMBAMC_LIBS
	 do
	   if test ${elem::2} = "-L" || test ${elem::2} = "-f"; then 
	     bblibdirs="$bblibdirs $elem"
           elif test ${elem::2} = "-l"; then
	     bblibdeps="$bblibdeps $elem"
	   fi
	 done
      ])
AC_SUBST([BAMBAMC_INCL], ["$bbincludes"])
AC_SUBST([BAMBAMC_LDFLG], ["$bblibdirs"])
AC_SUBST([BAMBAMC_DEPS],  ["$bblibdeps"])

# check for zlib
AS_IF([test x$for_intel_mic = xno],
	   [ AC_CHECK_LIB([z], [gzread],
             		      [withzlib=true;
              		      AC_DEFINE([HAVE_ZLIB], [1], [Define to 1 if you have zlib])],
	      		      [withzlib=false],) ],
           [withzlib=false;]
)
#AS_ECHO([x$withzlib])
# POSIX threads API with only named but not unnamed semaphores (MacOSX)
#AC_CHECK_LIB
#	(
#	[pthread],
#	[sem_getvalue],
#	[],
#	[AC_DEFINE([HAVE_POSIX_SEAPHOR_NAMEDONLY], [1], [Define to 1 if unnamed semaphores are not supported])
#	)

# Checks for header files.
AC_HEADER_STDC

# AC_CHECK_HEADERS([limits.h stdio.h stdlib.h string.h math.h])
AC_CHECK_HEADERS([math.h pthread.h semaphore.h float.h stddef.h stdint.h])
AS_IF([test x$withzlib = xtrue], AC_CHECK_HEADERS([zlib.h]))

# some versions of icc on 64-bit can't compile limits.h 
AC_CHECK_HEADER([limits.h], [withlimith=true], [withlimith=false], [])

# check for immintrin.h header (intel C compiler icc: prototypes for intrinsics)
AS_IF([test x$for_intel_mic = xyes],
      [AC_CHECK_HEADERS(
         [immintrin.h], 
	 [withintrin=true], 
         [withintrin=false;
	   AC_MSG_WARN([won't compile for Intel MIC])],
         [])],
      [withintrin=false;]
)
#AS_ECHO($withintrin)
# recogising the emmintrin.h header (some gcc versions require '-msse2' set)
AS_IF([test x$withintrin = xtrue],
      [withsse2=false],
      AC_CHECK_HEADERS([emmintrin.h], [withsse2=true], [withsse2=false], [])
)
#AS_ECHO($withsse2)
#establish 32-bit vs 64-bit -> $host_cpu
AC_CANONICAL_HOST
#AS_ECHO($host)
#AS_ECHO($host_cpu)
#AS_ECHO($host_vendor)
#AS_ECHO($host_os)

AC_CANONICAL_BUILD
#AS_ECHO($build)
#AS_ECHO($build_cpu)

# flag MacOSX Darwin which seems to build 32-bit binaries by default
AS_CASE([$build_os], [*darwin*], [build_macosx_darwin=true], [build_macosx_darwin=false])
AS_IF([test x$build_macosx_darwin = xtrue -a $host_cpu = x86_64], 
  [AC_MSG_NOTICE("compiling on darwin for x86_64 ...")])

# Intel's C compiler
# flag icc on 32 bit platform
AS_IF([test $CC = "icc" -a $host_cpu = "i686"], 
	    	AC_DEFINE([SMALT_AUTOCONF_ICC_OPT], [1], 
	    	[ICC on 32bit: avoid undefined references to __svml_exp2]))
AM_CONDITIONAL([ICC_COMPILER], [test $CC = "icc"])
AM_CONDITIONAL([ICC_NOLIMITH], [test $CC = "icc" -a x$withlimith = xfalse])
AM_CONDITIONAL([ICC_MMIC], [test x$for_intel_mic = xyes -a x$withintrin = xtrue])


AM_CONDITIONAL([SSE2], [test x$withsse2 = xtrue])
AM_CONDITIONAL([WITH_SIMD], [test x$withintrin = xtrue -o x$withsse2 = xtrue])
AM_CONDITIONAL([WITH_BAMBAMC], [test x$bambamc = xyes])
AM_CONDITIONAL([WITH_SAMTOOLS], [test x$samtools = xyes])
AM_CONDITIONAL([WITH_SAMTOOLS_EXEC], [test x$samtoolsexec = xyes])
AM_CONDITIONAL([WITH_ZLIB], [test x$withzlib = xtrue])
AM_CONDITIONAL([WITH_TESTALL], [test -e ./test/seqSetIO_test.c])
AM_CONDITIONAL([WITH_TESTDRIVER], [test x$with_testdriver = xtrue])
#AM_CONDITIONAL([WITH_SAMTOOLS], [test -n "$LIBSAMTOOLS"])
#AM_COND_IF([ICC_NOLIMITH], [AC_SUBST([AM_CFLAGS], [-D_GCC_LIMITS_H_])])
#AM_CONDITIONAL([APPLE_DARWIN], [test x$build_macosx_darwin = xtrue])
AM_CONDITIONAL([APPLE_64BIT], [test x$build_macosx_darwin = xtrue -a $host_cpu = x86_64])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
#AC_TYPE_INT8_T


# Checks for functions.
#AC_FUNC_MALLOC
#AC_FUNC_REALLOC
#AC_CHECK_FUNCS([memmove memset strtoul]) 

AC_CONFIG_FILES([Makefile src/Makefile test/Makefile test/data/Makefile misc/Makefile])

AC_OUTPUT
