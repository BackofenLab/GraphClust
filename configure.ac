#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([GraphClust], [0.7.6], [heyne@informatik.uni-freiburg.de], [GraphClust],
        [http://rna.bioinf.uni-freiburg.de/])

# Define here to keep all version information in one place.

AC_CONFIG_MACRO_DIR([m4])
m4_include([m4/openmp1.m4])

LT_INIT                     # for libtool
AC_PROG_LIBTOOL             # deprecated; left for backward compatibility
# AC_SUBST(LIBTOOL_DEPS)


if test -z $CXXFLAGS; then
    CXXFLAGS='-O3'
fi

CXXFLAGS='-O3'
CFLAGS='-O3'

# Checks for programs.
AC_PROG_CXX
AC_PROG_CXXCPP
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL

AC_LANG_CPLUSPLUS

AC_PATH_PROGS([PERL],[perl perl5])

AC_PATH_PROG([HELP2MAN],[help2man],false)
AC_PATH_PROG([POD2MAN],[pod2man],false)

AM_INIT_AUTOMAKE

# Vienna RNA package library path support, if not installed in standard directory/path
VRNA_OK=true;

AC_ARG_WITH([vrna],
    AC_HELP_STRING(
        [--with-vrna=PREFIX],
        [Installation directory of the Vienna RNA library]
      )
)


# register Vienna programs as configuration variables
AC_ARG_VAR([RNAfold],[RNAfold name (def=RNAfold)])
AC_ARG_VAR([RNAplfold],[RNAplfold name (def=RNAplfold)])
AC_ARG_VAR([RNAalifold],[RNAalifold name (def=RNAalifold)])

AC_ARG_VAR([VIENNA_BINPATH],[Vienna RNA Package Version >=2.0 (def="")])
AC_ARG_VAR([RNASHAPES_BINPATH],[RNAshapes (def=RNAshapes)])
AC_ARG_VAR([LOCARNA_BINPATH],[LocARNA Package >=1.7 (def=)])
AC_ARG_VAR([INFERNAL_BINPATH],[Infernal Package (def=)])
AC_ARG_VAR([R_BINPATH],[R Package (def=)])
AC_ARG_VAR([RNAZ_BINPATH],[RNAz >= 2.1 (def=RNAz)])
AC_ARG_VAR([BLASTCLUST_BINPATH],[NCBI BLASTCLUST (def=)])
AC_ARG_VAR([OCTAVE_BINPATH],[GNU Octave (def=)])
AC_ARG_VAR([CMFINDER_BINPATH],[CMfinder (def=)])

AC_PATH_PROG([RNASHAPES_BINPATH],[RNAshapes],[false])
AC_PATH_PROG([LOCARNA_BINPATH],[mlocarna],[false])
AC_PATH_PROG([INFERNAL_BINPATH],[cmsearch],[false])
AC_PATH_PROG([R_BINPATH],[R],[false])
AC_PATH_PROG([RNAZ_BINPATH],[RNAz],[false])
AC_PATH_PROG([BLASTCLUST_BINPATH],[blastclust],[false])
AC_PATH_PROG([OCTAVE_BINPATH],[octave],[false])
AC_PATH_PROG([CMFINDER_BINPATH],[cmfinder],[false])

#ac_cv_path_CMFINDER_BINPATH
#AC_CHECK_PROG

	AC_SUBST(RNAZ_BINPATH,["`echo $RNAZ_BINPATH | sed s/RNAz$//`"])
	AC_SUBST(BLASTCLUST_BINPATH,["`echo $BLASTCLUST_BINPATH | sed s/blastclust$//`"])
	AC_SUBST(LOCARNA_BINPATH,["`echo $LOCARNA_BINPATH | sed s/mlocarna$//`"])
	AC_SUBST(RNASHAPES_BINPATH,["`echo $RNASHAPES_BINPATH | sed s/RNAshapes$//`"])
	AC_SUBST(INFERNAL_BINPATH,["`echo $INFERNAL_BINPATH | sed s/cmsearch$//`"])
	AC_SUBST(R_BINPATH,["`echo $R_BINPATH | sed s/R$//`"])
	AC_SUBST(OCTAVE_BINPATH,["`echo $OCTAVE_BINPATH | sed s/octave$//`"])
	AC_SUBST(CMFINDER_BINPATH,["`echo $CMFINDER_BINPATH | sed s/cmfinder$//`"])


# check for required Vienna programs

if test "$VIENNA_BINPATH" != "" -a "$with_vrna" = "" ; then 
  with_vrna=$VIENNA_BINPATH 
fi

if test "$with_vrna" = "" ; then
   AC_PATH_PROG([RNAfold],[RNAfold],[false])
   AC_PATH_PROG([RNAplfold],[RNAplfold],[false])
   AC_PATH_PROG([RNAalifold],[RNAalifold],[false])
   AC_SUBST(VIENNA_BINPATH,["`echo $RNAfold | sed s/RNAfold$//`"])
   AC_SUBST(with_vrna,["`echo $RNAfold | sed s/RNAfold$//`"])
   AC_SUBST(with_vrna,["`echo $with_vrna | sed s/bin//`"])
else
   AC_SUBST(VIENNA_BINPATH,[$with_vrna/bin/])
   AC_PATH_PROG([RNAfold],[RNAfold],[false],[$with_vrna/bin])
   AC_PATH_PROG([RNAplfold],[RNAplfold],[false],[$with_vrna/bin])
   AC_PATH_PROG([RNAalifold],[RNAalifold],[false],[$with_vrna/bin])
fi



# trigger error if ViennaRNA programs not found
if test "$RNAfold" = "false" -o "$RNAplfold" = "false" -o "$RNAalifold" = "false" ; then
    VRNA_OK=false
fi

if echo "library(evd)" | R --vanilla >/dev/null ; then R_EVD="true" ; else
R_EVD="false" ; fi

AC_SUBST(LDFLAGS,[""])
AC_SUBST(CPPFLAGS,[""])
AC_SUBST(LIBS,[""])
AC_SUBST([with_vrna])


   # if --with-vrna given, add to flags for CPP and LD
#   if test "$with_vrna" != "" ; then
#      VRNA_LDFLAGS="-L$with_vrna/lib"
#      VRNA_CPPFLAGS="-I$with_vrna/include/ViennaRNA"

     # add vienna rna ld and cpp flags
#     LDFLAGS="$VRNA_LDFLAGS $LDFLAGS"
#     CPPFLAGS="$VRNA_CPPFLAGS $CPPFLAGS"
#   fi

#   # Check for Vienna RNA library (primitive check for one header)
#   AC_CHECK_HEADER([fold.h],
#		   [],
# 		   [VRNA_OK=false],[])

   ## link libraries
   LIBS="-lm $LIBS"
   
#   ## link the Vienna library
#   LIBS="-lRNA -lm $LIBS"
#   AC_CHECK_LIB(RNA, [warn_user], [VRNA_OK=true], [VRNA_OK=false])

## openmp is necessary for Vienna RNA 2.0
## and for NSPDK also
AC_OPENMP
AC_SUBST(OPENMP_CXXFLAGS)

CPPFLAGS="$OPENMP_CXXFLAGS $CPPFLAGS"
LDFLAGS="$OPENMP_CXXFLAGS $LDFLAGS"

#################
# Debug option
#
AC_MSG_CHECKING([whether to build with debug information])
debug_it=no
AC_ARG_ENABLE([debug],
    AC_HELP_STRING(
        [--enable-debug],
        [enable debug data generation (def=no)]
    ),
    debug_it="$enableval"
)
AC_MSG_RESULT($debug_it)
if test x"$debug_it" = x"yes"; then
    AC_DEFINE([DEBUG])
    CPPFLAGS="$CPPFLAGS -g -Wall"
else
    AC_DEFINE([NDEBUG])
fi

CPPFLAGS="$CPPFLAGS -Wno-deprecated"


####################
# Static linking
#
AC_MSG_CHECKING([whether to link statically])
link_statically=no
AC_ARG_ENABLE([link-static],
    AC_HELP_STRING(
        [--enable-link-static],
        [link statically (def=no)]
    ),
    link_statically="$enableval"
)

AC_MSG_RESULT($link_statically)
if test x"$link_statically" = x"yes"; then
   LDFLAGS="$LDFLAGS -all-static"
fi

# Checks for typedefs, structures, and compiler characteristics.
#AC_HEADER_STDBOOL
#AC_C_CONST
#AC_TYPE_SIZE_T
#AC_HEADER_STDC
#AC_CHECK_HEADERS(malloc.h strings.h unistd.h)
#AC_INLINE
#AC_CHECK_FUNCS(strdup strstr strchr erand48)
#AC_FUNC_MALLOC

AC_CONFIG_FILES([Makefile])

################################################################################

if test "$R_EVD" = "false"; then

 AC_MSG_NOTICE([=================================================================])
 AC_MSG_NOTICE([Warning: (compilation of GraphClust is still possible)])
 AC_MSG_NOTICE()

 AC_MSG_NOTICE([Cannot find library "evd" for R. Please install it!])
 AC_MSG_NOTICE()
 AC_MSG_NOTICE([Please call R and then type "install.packages("evd")])
 AC_MSG_NOTICE([May be you need to add a path to R_LIBS_USER!])
 AC_MSG_NOTICE()
fi


if test -o "$OCTAVE_BINPATH" = "false"; then

 AC_MSG_NOTICE([=================================================================])
 AC_MSG_NOTICE([Warning: (compilation of GraphClust is still possible)])
 AC_MSG_NOTICE()

 AC_MSG_NOTICE([Cannot find 'octave' tool!])
 AC_MSG_NOTICE()
 AC_MSG_NOTICE([Please install octave if you want to run GraphClust in evaluation mode!])
 AC_MSG_NOTICE([])
 AC_MSG_NOTICE()
fi

if test  "INFERNAL_BINPATH" = "false" -o "$RNASHAPES_BINPATH" = "false" -o "$LOCARNA_BINPATH" = "false" -o "$RNAZ_BINPATH" = "false"; then

 AC_MSG_NOTICE([=================================================================])
 AC_MSG_NOTICE([Error:])
 AC_MSG_NOTICE()

 AC_MSG_NOTICE([Cannot find all necessary external tools:])
 AC_MSG_NOTICE()
 AC_MSG_NOTICE([RNAshapes : $RNASHAPES_BINPATH])
 AC_MSG_NOTICE([Infernal  : $INFERNAL_BINPATH])
 AC_MSG_NOTICE([RNAz      : $RNAZ_BINPATH])
 AC_MSG_NOTICE([LocARNA   : $LOCARNA_BINPATH])
 AC_MSG_NOTICE()
 AC_MSG_NOTICE([Optional tools (only used with specific GraphClust parameters:])
 AC_MSG_NOTICE([BLASTCLUST: $BLASTCLUST_BINPATH])
 AC_MSG_NOTICE([R         : $R_BINPATH])
 AC_MSG_NOTICE([OCTAVE    : $OCTAVE_BINPATH])
 AC_MSG_NOTICE([CMfinder  : $CMFINDER_BINPATH])
 AC_MSG_NOTICE()
 AC_MSG_NOTICE([Please install all tools which were not found!])
 AC_MSG_NOTICE([Maybe you just need to add their path to \$PATH!?])
 AC_MSG_NOTICE()

 	exit -1;

fi


### ----------------------------------------
### Test if Vienna lib available and ok
## don't allow compiling otherwise
#
if test "$VRNA_OK" = "false"; then
   AC_MSG_NOTICE([=================================================================])
   AC_MSG_NOTICE([ERROR:])
   AC_MSG_NOTICE()

   if test "$with_vrna" != "" ; then
      AC_MSG_NOTICE([Cannot find Vienna RNA package in given path '$with_vrna'.])
   else
      AC_MSG_NOTICE([Cannot find Vienna RNA package in standard path. For using a])
      AC_MSG_NOTICE([non-standard path, please specify with --with-vrna=PREFIX.])
   fi

   AC_MSG_NOTICE()
   AC_MSG_NOTICE([The Vienna RNA package (version >= 2.0) is required for compiling.])
   AC_MSG_NOTICE([It is freely available at http://www.tbi.univie.ac.at/~ivo/RNA/.])
   AC_MSG_NOTICE()

   exit -1;
fi


AC_OUTPUT([pgma_src/Makefile NSPDK_src/Makefile Array/Makefile Statistics/Makefile GraphClust_config.pm])
