dnl acx_lang_package.m4 --- generic check for packaged software component
dnl                         consisting of header and library in directories
dnl                         that can be used by adding the flags in
dnl                         PACKAGE_LANG_INCLUDE and PACKAGE_LANG_LIB
dnl                         respectively to the language compiler command
dnl
dnl Copyright  (C)  2010  Thomas Jahns <jahns@dkrz.de>
dnl
dnl Version: 1.0
dnl Keywords:
dnl Author: Thomas Jahns <jahns@dkrz.de>
dnl Maintainer: Thomas Jahns <jahns@dkrz.de>
dnl URL: https://www.dkrz.de/redmine/projects/show/scales-ppm
dnl
dnl Redistribution and use in source and binary forms, with or without
dnl modification, are  permitted provided that the following conditions are
dnl met:
dnl
dnl Redistributions of source code must retain the above copyright notice,
dnl this list of conditions and the following disclaimer.
dnl
dnl Redistributions in binary form must reproduce the above copyright
dnl notice, this list of conditions and the following disclaimer in the
dnl documentation and/or other materials provided with the distribution.
dnl
dnl Neither the name of the DKRZ GmbH nor the names of its contributors
dnl may be used to endorse or promote products derived from this software
dnl without specific prior written permission.
dnl
dnl THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
dnl IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
dnl TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
dnl PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
dnl OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
dnl EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
dnl PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
dnl PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
dnl LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
dnl NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
dnl SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
dnl
dnl Commentary:
dnl
dnl
dnl
dnl Code:
dnl
dnl
dnl ACX_GENERIC_PACKAGE(PACKAGE, [INCLUDE], INC_FLAG, [EXTRA-INCLUDES],
dnl [EXTRA-INCLUDEFLAGS], [ACTION-IF_HEADER-NOT-FOUND], [FUNCTION],
dnl LIBFLAG, [LIB-CANDIDATES], [EXTRA-LIBS], [EXTRA-LIBFLAGS],
dnl [ACTION-IF-LIB-NOT-FOUND], [DEFAULT-ROOT])
dnl -------------------------------------------------------------------
dnl Check wether INCLUDE can be compiled and FUNCTION is found in
dnl LIB-CANDIDATES with current language compiler. Sets PACKAGE_LANG_LIB
dnl and PACKAGE_LANG_INCLUDE variables to necessary compiler
dnl switches and arguments such that inclusion/compilation succeed for
dnl program using INCLUDE or FUNCTION respectively.  Also defines
dnl configure --with arguments for PACKAGEROOT, PACKAGE-LIB and
dnl PACKAGE-INCLUDE.
dnl The code sets the shell variable have_PACKAGE__AC_LANG_ABBREV_bindings
dnl to yes if all requested tests succeeded, to no if any test failed.
dnl The library check will not run if the header check failed.
AC_DEFUN([ACX_GENERIC_PACKAGE],
  [AC_REQUIRE([_ASX_TR_ARG_PREPARE])
   AS_VAR_PUSHDEF([acx_pkg_root],[AS_TR_CPP([$1][root])])
   AS_VAR_PUSHDEF([acx_pkg_lib],[AS_TR_CPP([$1_]_AC_LANG_ABBREV[_LIB])])
   AS_VAR_PUSHDEF([acx_pkg_inc],[AS_TR_CPP([$1_]_AC_LANG_ABBREV[_INCLUDE])])
   AS_VAR_PUSHDEF([acx_cv_pkg_inc],
     [AS_TR_SH([acx_cv_]_AC_LANG_ABBREV[_include_]$2)])
   AS_VAR_PUSHDEF([acx_pkg_bindings],
     [AS_TR_SH([have_][$1_]_AC_LANG_ABBREV[_bindings])])
   AS_VAR_SET([acx_pkg_bindings],[yes])
   AC_SUBST(acx_pkg_root)
   AC_ARG_WITH(ASX_TR_ARG([$1-root]),
     [AS_HELP_STRING([--with-]ASX_TR_ARG([$1])[-root],
        [set directory to search for $1 headers and library]m4_ifval([$13],[, @<:@default=$13@:>@]))],
     [AS_VAR_SET([acx_pkg_root],"$AS_TR_SH([with_]ASX_TR_ARG([$1])[_root])")],
     m4_ifval([$13], AS_VAR_SET([acx_pkg_root],[$13])))
   AS_VAR_SET_IF([acx_pkg_root],
     [AS_VAR_SET_IF([acx_pkg_lib],,
        [AS_VAR_COPY([acx_temp],acx_pkg_root)
         AS_VAR_SET([acx_pkg_lib],["$8$acx_temp/lib"])])
      AS_VAR_SET_IF([acx_pkg_inc],,
        [AS_VAR_COPY([acx_temp],acx_pkg_root)
         AS_VAR_SET([acx_pkg_inc],["$3$acx_temp/include"])])])
   m4_ifval([$2],
     [AC_ARG_WITH(ASX_TR_ARG([$1-include]),
        [AS_HELP_STRING([--with-[]ASX_TR_ARG([$1])[]-include],
           [specifically set directory to search for $1 headers, ]dnl
[@<:@default=$]AS_TR_SH(ASX_TR_ARG([with_$1_root]))[/include@:>@])],
        [AS_VAR_SET([acx_pkg_inc],["]$3[$]AS_TR_SH(ASX_TR_ARG([with_$1_include]))["])])
      AC_ARG_VAR(acx_pkg_inc,dnl
[specifically set flags to use when compiling sources
using $1 includes.])
      ACX_LANG_CHECK_INCLUDE_PATHS_IFELSE([$2],[],
        [AS_VAR_COPY([acx_temp],[acx_cv_pkg_inc])
         acx_temp=`echo "$acx_temp" | sed -e 's/^ *//;s/ *$//'`
         AS_VAR_COPY([acx_pkg_inc],[acx_temp])],
        [AS_VAR_SET([acx_pkg_bindings],[no])
         $6],
        [$4],[$]acx_pkg_inc,m4_ifval([$5],[[$5]],[[[]]]))])
   m4_ifval([$7],
     [AC_ARG_WITH(ASX_TR_ARG([$1-lib]),
        [AS_HELP_STRING([--with-]ASX_TR_ARG([$1])[-lib],
        [specifically set directory to search for $1 library, ]dnl
[@<:@default=$]AS_TR_SH(ASX_TR_ARG([with_$1_root]))[/lib@:>@])],
        [AS_VAR_SET([acx_pkg_lib],["$8$AS_TR_SH(ASX_TR_ARG([with_$1_lib]))"])])
      AC_ARG_VAR(acx_pkg_lib,
        [specifically set flags to use when linking $1.])
      AC_SUBST(acx_pkg_lib)
      AS_VAR_IF([acx_pkg_bindings],[yes],
        [AS_VAR_PUSHDEF([ac_Search],[acx_cv_option_search_$7_]_AC_LANG_ABBREV)dnl
         ACX_OPTION_SEARCH_LIBS_MULTI([$7],[$9],,dnl
           [AS_VAR_SET([acx_pkg_bindings],[no])
            $12],[$10],m4_ifval([$11],[$11 ])$[acx_pkg_lib])
         AS_VAR_COPY([acx_temp],[ac_Search])
         acx_temp=`echo "$acx_temp" | sed -e 's/^ *//;s/ *$//'`
         AS_VAR_COPY([acx_pkg_lib],[acx_temp])
         AS_VAR_POPDEF([ac_Search])dnl
        ])
     ])dnl
   AS_VAR_POPDEF([acx_pkg_bindings])
   AS_VAR_POPDEF([acx_cv_pkg_inc])
   AS_VAR_POPDEF([acx_pkg_inc])
   AS_VAR_POPDEF([acx_pkg_lib])
   AS_VAR_POPDEF([acx_pkg_root])
 ])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://www.dkrz.de/redmine/projects/show/scales-ppm"
dnl license-default: "bsd"
dnl End:
