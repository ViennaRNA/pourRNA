
m4_define([pourRNA_version_major], [1])
m4_define([pourRNA_version_minor], [0])
m4_define([pourRNA_version_patch], [1])

AC_PREREQ([2.59])
# 5 argument version only available with aclocal >= 2.64
AC_INIT( [pourRNA], pourRNA_version_major.pourRNA_version_minor.pourRNA_version_patch, [], [pourRNA], [http://www.bioinf.uni-freiburg.de] )

AC_CANONICAL_HOST
AC_CONFIG_AUX_DIR([.])
AC_CONFIG_SRCDIR([src/pourRNA.cpp])
AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_MACRO_DIR([m4])

# backwards compatibility for autoconf >= 2.64
# PACKAGE_URL should be the fifth argument of AC_INIT
m4_ifndef([AC_PACKAGE_URL],
	[m4_define([AC_PACKAGE_URL], [http://www.bioinf.uni-freiburg.de/])
	AC_DEFINE(PACKAGE_URL, "AC_PACKAGE_URL", [Package URL])
	AC_SUBST(PACKAGE_URL, AC_PACKAGE_URL)])

lt_enable_auto_import=""
case "$host_os" in 
	cygwin* | mingw* | cegcc*)
		AM_LDFLAGS="-Wl,--enable-auto-import $AM_LDFLAGS"
esac

# Checks for programs

# check for C++ compiler
# store current compiler flags to avoid default setup via AC_PROG_CXX and *_CC
OLD_CXXFLAGS=$CXXFLAGS
OLD_CFLAGS=$CFLAGS
# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

AC_CHECK_PROG([HAVE_PKG_CONFIG], [pkg-config], [yes],[no])

# reset compiler flags to initial flags
CXXFLAGS=$OLD_CXXFLAGS
CFLAGS=$OLD_CFLAGS

# automake initialisation (mandatory) and check for minimal automake API version
AM_INIT_AUTOMAKE([1.9])

# use the C++ compiler for the following checks
AC_LANG([C++])

# include openMP for multithreading.
dnl Check for OpenMP and prepare flags
#AX_OPENMP([],[]) # not used

#AC_SUBST(OPENMP_CXXFLAGS)
AC_SUBST(AM_CXXFLAGS,"$OPENMP_CXXFLAGS")

###############################################################################
###############################################################################

############  PARAMETERS  ########################################

###############################################################################
# SETUP FOR COMPILATION WITH 2011 ISO C++ STANDARD 
###############################################################################

AX_CXX_COMPILE_STDCXX_11([noext],[mandatory])

###############################################################################
# DEBUG SUPPORT SETUP
###############################################################################

AC_MSG_CHECKING([whether to build with debug information])
debuger=no
AC_ARG_ENABLE([debug],
	[AS_HELP_STRING([--enable-debug],
	    [enable debug data generation (def=no)])],
  	[debuger="$enableval"])
AC_MSG_RESULT([$debuger])
if test x"$debuger" = x"yes"; then
	AC_DEFINE([DEBUG], [], [Run in DEBUG mode with additional assertions and debug output])
	AM_CXXFLAGS="$AM_CXXFLAGS -pg -g -Wno-uninitialized -O0 -pthread"
	AM_LDFLAGS="$AM_LDFLAGS -g -Wno-uninitialized -O0"
	#AM_CXXFLAGS="$AM_CXXFLAGS -g -O0 -Wall -Werror -Wno-uninitialized -Wno-deprecated"
else
	AC_DEFINE([NDEBUG], [], [Run in normal mode with minimal assertions])
	AM_CXXFLAGS="$AM_CXXFLAGS -O3 -fno-strict-aliasing -Wno-uninitialized -Wno-deprecated -pthread"
fi


###############################################################################
# Vienna RNA package library path support, if not installed in usual directories
###############################################################################
AC_ARG_WITH([RNA],
	[AC_HELP_STRING(
	    [--with-RNA=PREFIX],
	    [alternative prefix path to Vienna RNA library]
	  )],
	  [RNAPATHSET=1],
	  [RNAPATHSET=0]
)
if test  $RNAPATHSET = 1 ; then
	# set compiler and linker flags if needed
	AM_CXXFLAGS="-I$with_RNA/include $AM_CXXFLAGS"
	AM_LDFLAGS="-L$with_RNA/lib $AM_LDFLAGS"
fi

###############################################################################
###############################################################################


############  CHECKS  ############################################

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Checks for header files.
AC_HEADER_STDC


###############################################################################
# BEGIN VIENNA CHECK
###############################################################################
# check for Vienna RNA headers
AC_MSG_CHECKING([for the Vienna RNA package >=2.4.11 headers])
OLD_CXXFLAGS=$CXXFLAGS
OLD_LDFLAGS=$LDFLAGS
CXXFLAGS="$CXXFLAGS $AM_CXXFLAGS"
LDFLAGS="$LDFLAGS $AM_LDFLAGS"
AC_COMPILE_IFELSE(
     [AC_LANG_SOURCE([[#include <ViennaRNA/neighbor.h>]])],
     [ 
		AC_MSG_RESULT([yes])
		AC_CHECK_LIB([RNA], [vrna_neighbors], [RNANOTFOUND=0], [ RNANOTFOUND=1])
     ],
     [
       AC_MSG_RESULT([no])
       RNANOTFOUND=1;
     ]
)
CXXFLAGS=$OLD_CXXFLAGS
LDFLAGS=$OLD_LDFLAGS

# error output if ViennaRNA not found
if test "$RNANOTFOUND" = "1"; then
	AC_MSG_NOTICE()
	AC_MSG_NOTICE([The Vienna RNA C library >= 2.4.11 is required.])
	AC_MSG_NOTICE([ -> It can be obtained from http://www.tbi.univie.ac.at/.])
	AC_MSG_NOTICE()
 	if test "$RNAPATHSET" = "1"; then
		AC_MSG_NOTICE([ -> Can't find the Vienna RNA library in given path '$with_RNA'.])
 	else
		AC_MSG_NOTICE([ -> If installed in a non-standard path, please use '--with-RNA=PREFIX'.])
 	fi
 	DEPENDENCYNOTFOUND=1;
else
	# register Vienna RNA lib for linking
	LIBS="$LIBS -lRNA"
fi


###############################################################################
# END VIENNA CHECK
###############################################################################

###############################################################################
# FINAL DEPENDENCY CHECK AND EXIT IF NEEDED
###############################################################################

# error ABORT if on of the libraries was not found
if test "$DEPENDENCYNOTFOUND" = "1"; then
	AC_MSG_NOTICE()
	AC_MSG_ERROR([Some dependency was not met! See above for errors and relate to './configure --help'.])
fi

##########################################################################

# distribute additional compiler and linker flags
# --> set these variables instead of CXXFLAGS or LDFLAGS
AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_LDFLAGS])
AC_SUBST([LIBS])

# files to generate via autotools (.am or .in source files)
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/Makefile])
AC_CONFIG_FILES([misc/Makefile])

# generate the final Makefile etc.
AC_OUTPUT([packaging/pourRNA.spec packaging/PKGBUILD])

