
m4_define([pourRNA_version_major], [1])
m4_define([pourRNA_version_minor], [3])
m4_define([pourRNA_version_patch], [1])

AC_PREREQ([2.64])
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


# reset compiler flags to initial flags
CXXFLAGS=$OLD_CXXFLAGS
CFLAGS=$OLD_CFLAGS

# automake initialisation (mandatory) and check for minimal automake API version
AM_INIT_AUTOMAKE([1.9])

# use the C++ compiler for the following checks
AC_LANG([C++])

# pkg-config is used to obtain the flags the ViennaRNA library was built with
PKG_PROG_PKG_CONFIG

# OpenMP flags of the C++ compiler; only needed as fallback when linking a
# libRNA that was built with OpenMP but no RNAlib2.pc is available
AC_OPENMP

# gengetopt is only needed when the command line parser is not yet generated
# (e.g. git checkout); release tarballs ship pourRNA_cmdl.c/.h
AC_CHECK_PROG([GENGETOPT], [gengetopt], [gengetopt], [no])
if test "x$GENGETOPT" = "xno" && test ! -f "$srcdir/src/pourRNA_cmdl.c"; then
	AC_MSG_ERROR([gengetopt is required to build from a source checkout (src/pourRNA_cmdl.c is missing). Please install gengetopt.])
fi

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
	[AS_HELP_STRING([--with-RNA=PREFIX],
	    [alternative prefix path to Vienna RNA library])],
	[RNAPATHSET=1],
	[RNAPATHSET=0]
)
if test "$RNAPATHSET" = "1" ; then
	# let pkg-config find RNAlib2.pc of the given installation first
	PKG_CONFIG_PATH="$with_RNA/lib/pkgconfig:$with_RNA/lib64/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"
	export PKG_CONFIG_PATH
fi

###############################################################################
###############################################################################


############  CHECKS  ############################################

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T


###############################################################################
# BEGIN VIENNA CHECK
###############################################################################
# 1) preferred: ask pkg-config for the exact flags libRNA was built with
#    (covers OpenMP, pthreads, GSL, MPFR, ... depending on its configuration)
# 2) fallback: plain -I/-L/-lRNA plus OpenMP and math library
RNANOTFOUND=1
RNA_FROM_PKGCONFIG=no

# libRNA is usually a static library -> request flags for static linking
RNA_SAVED_PKG_CONFIG=$PKG_CONFIG
if test -n "$PKG_CONFIG"; then
	PKG_CONFIG="$PKG_CONFIG --static"
fi
PKG_CHECK_MODULES([RNA], [RNAlib2 >= 2.4.11],
	[RNA_FROM_PKGCONFIG=yes],
	[RNA_FROM_PKGCONFIG=no])
PKG_CONFIG=$RNA_SAVED_PKG_CONFIG

if test "$RNA_FROM_PKGCONFIG" = "no"; then
	AC_MSG_NOTICE([RNAlib2.pc not found, falling back to manual ViennaRNA detection])
	RNA_CFLAGS=""
	RNA_LIBS="-lRNA $OPENMP_CXXFLAGS -lm"
	if test "$RNAPATHSET" = "1" ; then
		RNA_CFLAGS="-I$with_RNA/include"
		RNA_LIBS="-L$with_RNA/lib $RNA_LIBS"
	fi
fi

# verify that a program using ViennaRNA really compiles and links.
# vrna_fold_compound() + vrna_neighbors() pull in most of the library
# (incl. the parts that depend on OpenMP/SVM), unlike a bare symbol check.
AC_MSG_CHECKING([whether ViennaRNA >= 2.4.11 can be used])
OLD_CXXFLAGS=$CXXFLAGS
OLD_LIBS=$LIBS
CXXFLAGS="$CXXFLAGS $RNA_CFLAGS"
LIBS="$RNA_LIBS $LIBS"
AC_LINK_IFELSE(
	[AC_LANG_PROGRAM([[
extern "C" {
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/neighbor.h>
#include <ViennaRNA/mfe.h>
}
]], [[
  vrna_fold_compound_t *fc = vrna_fold_compound("GGGAAACCC", NULL, VRNA_OPTION_DEFAULT);
  float mfe = vrna_mfe(fc, NULL);
  vrna_move_t *n = vrna_neighbors(fc, NULL, VRNA_MOVESET_DEFAULT);
  free(n);
  vrna_fold_compound_free(fc);
  return (int)mfe;
]])],
	[RNANOTFOUND=0; AC_MSG_RESULT([yes (flags from pkg-config: $RNA_FROM_PKGCONFIG)])],
	[RNANOTFOUND=1; AC_MSG_RESULT([no])])
CXXFLAGS=$OLD_CXXFLAGS
LIBS=$OLD_LIBS

# error output if ViennaRNA not found
if test "$RNANOTFOUND" = "1"; then
	AC_MSG_NOTICE()
	AC_MSG_NOTICE([The Vienna RNA C library >= 2.4.11 is required.])
	AC_MSG_NOTICE([ -> It can be obtained from https://www.tbi.univie.ac.at/RNA/])
	AC_MSG_NOTICE()
	if test "$RNAPATHSET" = "1"; then
		AC_MSG_NOTICE([ -> Can't use the Vienna RNA library in given path '$with_RNA'.])
	else
		AC_MSG_NOTICE([ -> If installed in a non-standard path, please use '--with-RNA=PREFIX'])
		AC_MSG_NOTICE([    or set PKG_CONFIG_PATH to the directory containing RNAlib2.pc.])
	fi
	AC_MSG_NOTICE([ -> See config.log for the exact compiler/linker error.])
	DEPENDENCYNOTFOUND=1;
else
	# register Vienna RNA lib for compiling and linking
	AM_CXXFLAGS="$RNA_CFLAGS $AM_CXXFLAGS"
	LIBS="$RNA_LIBS $LIBS"
fi


###############################################################################
# END VIENNA CHECK
###############################################################################

###############################################################################
# help2man check
###############################################################################
AX_PROG_HELP2MAN([src/pourRNA])

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
AC_CONFIG_FILES([doc/Makefile])
AC_CONFIG_FILES([packaging/pourRNA.spec packaging/PKGBUILD])

# generate the final Makefile etc.
AC_OUTPUT

