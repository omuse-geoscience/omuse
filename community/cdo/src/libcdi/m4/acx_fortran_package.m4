dnl acx_fortran_package.m4 --- aggregate check for availability of
dnl                            include file(s) and library routines
dnl                            associated with a given package
dnl
dnl Copyright  (C)  2010  Thomas Jahns <jahns@dkrz.de>
dnl
dnl Version: 1.0
dnl Keywords: fortran package availability check
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
dnl ACX_FORTRAN_PACKAGE(PACKAGE,
dnl    [INCLUDE], [EXTRA-INCLUDES], [EXTRA-INCLUDEFLAGS],
dnl       [ACTION-IF_HEADER-NOT-FOUND],
dnl    [FUNCTION], [LIB-CANDIDATES], [EXTRA-LIBS], [EXTRA-LIBFLAGS],
dnl      [ACTION-IF-LIB-NOT-FOUND],
dnl    [DEFAULT-ROOT])
dnl -------------------------------------------------------------------
dnl Check wether INCLUDE can be compiled and FUNCTION is found in
dnl LIB-CANDIDATES. Sets PACKAGE_LIB and PACKAGE_INCLUDE variables to
dnl necessary Fortran compiler switches and arguments such that
dnl inclusion/compilation succeed for program using INCLUDE or
dnl FUNCTION respectively.
dnl Also defines configure --with arguments for PACKAGEROOT,
dnl PACKAGE-LIB and PACKAGE-INCLUDE.
AC_DEFUN([ACX_FORTRAN_PACKAGE],
  [AC_LANG_PUSH([Fortran])
   AC_REQUIRE([AC_PROG_FPP])
   ACX_GENERIC_PACKAGE([$1],[$2],[$FPP_INCOPT],[$3],[$4],[$5],[$6],[-L],[$7],[$8],[$9],[$10],[$11],[$12])
   AC_LANG_POP([Fortran])dnl
  ])
dnl
dnl
dnl ACX_F90_PACKAGE: same as ACX_FORTRAN_PACKAGE but uses module
dnl instead of header
dnl
dnl ACX_F90_PACKAGE(PACKAGE,
dnl    MODULE, [EXTRA-INCLUDES], [EXTRA-INCLUDEFLAGS],
dnl       [ACTION-IF-MODULE-NOT-FOUND],
dnl    SUBROUTINE, [LIB-CANDIDATES], [EXTRA-LIBS], [EXTRA-LIBFLAGS],
dnl      [ACTION-IF-LIB-NOT-FOUND],
dnl    [DEFAULT-ROOT],[EXTRA-DECLARATION],[SUBROUTINE-ARGUMENTS])
dnl -------------------------------------------------------------------

dnl Check wether USE MODULE can be compiled and FUNCTION is found in
dnl LIB-CANDIDATES. Sets PACKAGE_LIB and PACKAGE_MOD variables to
dnl necessary Fortran compiler switches and arguments such that
dnl compilation and linking succeed for programs using MODULE and
dnl calling SUBROUTINE respectively.  Extra declaration will be
dnl inserted between use MODULE and the subroutine call.
dnl
dnl Also defines configure --with arguments for PACKAGEROOT,
dnl PACKAGE_FC_LIB and PACKAGE_FC_MOD.
AC_DEFUN([ACX_F90_PACKAGE],
  [AC_LANG_PUSH([Fortran])
   AC_REQUIRE([AC_PROG_FC])
   AC_REQUIRE([ACX_SL_FC_CHECK_MOD_PATH_FLAG])
   AC_REQUIRE([_ASX_TR_ARG_PREPARE])
   AC_SUBST(AS_TR_CPP([$1][root]))
   AC_ARG_WITH(ASX_TR_ARG([$1-root]),
     [AS_HELP_STRING([--with-]ASX_TR_ARG([$1])[-root],
        [set directory to search for $1 headers and library, @<:@default=$11@:>@])],
        [AS_TR_CPP([$1][root])="$AS_TR_SH([with_]ASX_TR_ARG([$1])[_root])"],
        m4_ifval([$11], [AS_TR_CPP([$1][root])=$11]))
   AS_VAR_SET_IF(AS_TR_CPP([$1][root]),
     [AS_VAR_SET_IF([AS_TR_CPP([$1_]_AC_LANG_ABBREV[_LIB])],,
        [AS_TR_CPP([$1_]_AC_LANG_ABBREV[_LIB])="-L$[]AS_TR_CPP([$1][root])/lib"])
      AS_VAR_SET_IF([AS_TR_CPP([$1_]_AC_LANG_ABBREV[_MOD])],,
        [AS_TR_CPP([$1_]_AC_LANG_ABBREV[_MOD])="$FC_MOD_FLAG$[]AS_TR_CPP([$1][root])/include"])])
   m4_ifval([$2],dnl
     [AC_ARG_WITH(ASX_TR_ARG([$1-]_AC_LANG_ABBREV[-mod]),
        [AS_HELP_STRING([--with-]ASX_TR_ARG([$1-]_AC_LANG_ABBREV)[-mod],
           [specifically set directory to search for $1 module, ]dnl
[@<:@default=$]AS_TR_SH(ASX_TR_ARG([with_$1_root]))[/include@:>@])],
        AS_TR_CPP([$1_]_AC_LANG_ABBREV[_MOD])[="$FC_MOD_FLAG$AS_TR_SH(ASX_TR_ARG([with_$1_]_AC_LANG_ABBREV[_mod]))"],
        [])
      AC_ARG_VAR(AS_TR_CPP([$1_]_AC_LANG_ABBREV[_MOD]),
        [flags to enable 'USE $1' in Fortran program.])
      ACX_FORTRAN_CHECK_MOD_PATHS_IFELSE([$2],["$[]AS_TR_CPP([$1][root])/include" "$[]AS_TR_CPP([$1][root])/lib"],
        AS_TR_CPP([$1_]_AC_LANG_ABBREV[_MOD])[=]AS_VAR_GET([acx_Mod]),
        [$5],
        [$3],m4_ifval([$4],[$4 ])[$[]AS_TR_CPP([$1_]_AC_LANG_ABBREV[_MOD])])])
   m4_ifval([$6],
     [AC_ARG_WITH(ASX_TR_ARG([$1]_AC_LANG_ABBREV[-lib]),
        [AS_HELP_STRING([--with-]ASX_TR_ARG([$1-]_AC_LANG_ABBREV)[-lib],
        [specifically set directory to search for $1 library, ]dnl
[@<:@default=$]AS_TR_SH(ASX_TR_ARG([with_$1_root]))[/lib@:>@])],
        AS_TR_CPP([$1_]_AC_LANG_ABBREV[_LIB])[="-L$AS_TR_SH(ASX_TR_ARG([with_$1_]_AC_LANG_ABBREV[_lib]))"],
        [])
      AS_IF([test x$]AS_TR_SH([have_][$1_]_AC_LANG_ABBREV[_bindings])[ = xyes],
        [AS_VAR_PUSHDEF([ac_Search], [acx_cv_option_search_$6_]_AC_LANG_ABBREV)dnl
         AC_SUBST(AS_TR_CPP([$1_]_AC_LANG_ABBREV[_LIB]))
         acx_save_FCFLAGS="$FCFLAGS"
         FCFLAGS="$[]AS_TR_CPP([$1_]_AC_LANG_ABBREV[_MOD]) $FCFLAGS"
         ACX_OPTION_SEARCH_LIBS_MULTI([$6],[$7],,dnl
           [$10],[$8],m4_ifval([$9],[$9 ])$[]AS_TR_CPP([$1_]_AC_LANG_ABBREV[_LIB]),
           [      use $2]m4_ifval([$12],[[
      $12]]),[      call $6$13])
         FCFLAGS="$acx_save_FCFLAGS"
         AS_TR_CPP([$1_]_AC_LANG_ABBREV[_LIB])="AS_VAR_GET([ac_Search])"
         AS_VAR_POPDEF([ac_Search])])])
   AC_LANG_POP([Fortran])dnl
  ])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://www.dkrz.de/redmine/projects/show/scales-ppm"
dnl license-default: "bsd"
dnl End:
