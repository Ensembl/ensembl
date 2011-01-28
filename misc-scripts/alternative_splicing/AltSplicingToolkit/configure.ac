#                                                       -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

# Copyright (c) 1999-2011 The European Bioinformatics Institute and 
# Genome Research Limited, and others.  All rights reserved. 

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# Ensembl License for more details.

# You should have received a copy of the Ensembl License
# along with this program.  If not, see 
# <http://www.ensembl.org/info/about/code_licence.html>.

dnl --------------------------------
dnl run aclocal first
dnl run autoconf
dnl automake --add-missing
dnl --------------------------------

dnl --------------------------------
dnl Initialization macros.
dnl Every configure script must call
dnl AC_INIT before doing something else
dnl --------------------------------

dnl Require Autoconf 2.60
AC_PREREQ([2.60])
AC_INIT([Ensembl AltSplicingToolkit], [0.5.2], [koscieln@ebi.ac.uk])

dnl -----------------------------------------------
dnl Package name and version number (user defined)
dnl for automake
dnl -----------------------------------------------

AM_INIT_AUTOMAKE
#([1.0 -Wall -Werror foreign])

dnl -----------------------------------------------
dnl Finding configure Input
dnl -----------------------------------------------

AC_CONFIG_SRCDIR([src/Main.cpp])
AC_CONFIG_HEADER([include/config.hpp])

AC_SUBST(VERSION)

ISODATE=`date +%Y-%m-%d`
AC_SUBST(ISODATE)

dnl -----------------------------------------------
dnl Checks for programs.
dnl -----------------------------------------------
# We are using C++
#AC_LANG(C++)
# autoconf2.13: AC_LANG_CPLUSPLUS
#AC_LANG_CPLUSPLUS
#AC_REQUIRE_CPP
AC_PROG_CXX
#CPPFLAGS([-O2 -Wno-deprecated])

AC_PROG_INSTALL
AC_PROG_RANLIB
AC_PROG_CC_C_O

# Turn off shared libraries during beta-testing, since they
# make the build process take too long
#LT_INIT
#AC_PROG_LIBTOOL
#AC_DISABLE_SHARED
#AM_PROG_LIBTOOL
AM_SANITY_CHECK

## ----- ##
## Make. ##
## ----- ##
#AC_PROG_MAKE_CASE_SENSITIVE

dnl These are standard among *NIX systems, but not when cross-compiling
#AC_CHECK_TOOL(RANLIB, ranlib, ':')
#AC_CHECK_TOOL(LD, ld, ld)
#AC_CHECK_TOOL(AR, ar, ar)


dnl -----------------------------------------------
dnl Checks for libraries.
dnl -----------------------------------------------

# Paths to examplelib stuff:
AC_SUBST([LDFLAGS],["${LDFLAGS} -L."])

##
## check for log4cpp
##
AC_ARG_WITH([log4cpp],
            AC_HELP_STRING([--with-log4cpp=LOG4CPP_PATH],
                           [Supply the location of Log for C++ installation]),
            [ql_log4cpp_path=$withval],
            [ql_log4cpp_path=""])
if test [ -n "$ql_log4cpp_path/include/log4cpp" ] ; then
   AC_SUBST([CPPFLAGS],["${CPPFLAGS} -I${ql_log4cpp_path}/include"])
fi
if test [ -n "$ql_log4cpp_path/lib" ] ; then
   AC_SUBST([LDFLAGS],["${LDFLAGS} -L${ql_log4cpp_path}/lib -llog4cpp"])
fi

## check for boost, this segment as well as segment below are 
## from QuantLib's configure.ac
## 
##
## Help the compiler find external software
AC_ARG_WITH([boost-include],
            AC_HELP_STRING([--with-boost-include=BOOST_INCLUDE_PATH],
                           [Supply the location of Boost header files]),
            [ql_boost_include_path=$withval],
            [ql_boost_include_path=""])
if test [ -n "$ql_boost_include_path/boost" ] ; then
   AC_SUBST([CPPFLAGS],["${CPPFLAGS} -I${ql_boost_include_path}"])
   #AC_SUBST([CXXFLAGS],["${CXXFLAGS} -I${ql_boost_include_path}"])
fi
AC_ARG_WITH([boost-lib],
            AC_HELP_STRING([--with-boost-lib=BOOST_LIB_PATH],
                           [Supply the location of Boost libraries]),
            [ql_boost_lib_path=$withval],
            [ql_boost_lib_path=""])
if test [ -n "$ql_boost_lib_path" ] ; then
   AC_SUBST([LDFLAGS],["${LDFLAGS} -L${ql_boost_lib_path} -lboost_program_options -lboost_regex"])
fi



# libboost_regex-gcc41-mt-1_45.a

#PKG_CHECK_MODULES(DEPS, libboost_regex-gcc41-mt-1_45  >= 1.45.0 )
#AC_SUBST(DEPS_CFLAGS)
#AC_SUBST(DEPS_LIBS)

#AC_CHECK_LIB(boost_regex-gcc, main,, [
	#		      AC_CHECK_LIB(boost_regex, main,,
		#	      AC_MSG_ERROR([** cannot find required boost_regex library. Get it from www.boost.org **]))
			#      ])

#AC_CHECK_HEADERS(boost/regex.hpp, , [
#				  echo "######################"
#				  echo " Boost headers needed "
#				  echo "######################"
#				  exit 1
#				  ])

#AC_CHECK_LIB(boost_program_options-gcc, main,, [
#			      AC_CHECK_LIB(boost_program_options, main,,
#			      AC_MSG_ERROR([** cannot find required boost_program_options library. Get it from www.boost.org **]))
#			      ])

#AC_CHECK_HEADERS(boost/program_options.hpp, , [
#				  echo "######################"
#				  echo " Boost headers needed "
#				  echo "######################"
#				  exit 1
#				  ])




#AC_SUBST(CXXFLAGS)
#AC_SUBST(LDFLAGS)
AC_SUBST([ARFLAGS],["${ARFLAGS} -rfs"])
AC_OUTPUT( [Makefile src/Makefile] )
#dnl src/as/Makefile src/gff/Makefile src/util/Makefile] )

