dnl acx_option_search_libs.m4 --- Search for library needed only by some
dnl                               sub-components and do not add corresponding
dnl                               flags to LIBS but give the user control with
dnl                               respect to how the test result is used
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
dnl _ACX_OPTION_SEARCH_LIBS(FUNCTION, SEARCH-LIBS,
dnl   [OTHER-LIBRARIES], [EXTRA-FLAGS], [PREAMBLE], [CALL-CODE])
dnl same as ACX_OPTION_SEARCH_LIBS but does not cache result
dnl
dnl Uses either AC_LANG_PROGRAM([PREAMBLE],CALL-CODE) or
dnl AC_LANG_CALL([PREAMBLE],FUNCTION), depending on wether CALL-CODE
dnl is given or not.
AC_DEFUN([_ACX_OPTION_SEARCH_LIBS],
  [acx_option_func_search_save_LIBS="$LIBS"
   AC_LANG_CONFTEST([m4_ifval([$6],[AC_LANG_PROGRAM([$5],[$6])],[AC_LANG_CALL([$5], [$1])])])
   for ac_lib in '' $2; do
     AS_IF([test -z "$ac_lib"],
       [ac_res="none required"
        LIBS="m4_ifval([$4],[$4 ])m4_ifnblank($3,[$3 ])$acx_option_func_search_save_LIBS"],
       [ac_res="-l$ac_lib"
        LIBS="m4_ifval([$4],[$4 ])$ac_res m4_ifnblank($3,[$3 ])$acx_option_func_search_save_LIBS"])
     AC_LINK_IFELSE([], [AS_IF([test x"$ac_res" = x"none required"],dnl
        [AS_VAR_SET([ac_Search],["]m4_ifval([$4],[$4 ])[$3"])],dnl
        [AS_VAR_SET([ac_Search],["]m4_ifval([$4],[$4 ])[-l$ac_lib $3"])])])
     AS_VAR_SET_IF([ac_Search], [break])
   done
   rm conftest.$ac_ext
   LIBS="$acx_option_func_search_save_LIBS"])
dnl ACX_OPTION_SEARCH_LIBS(FUNCTION, SEARCH-LIBS,
dnl   [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND], [OTHER-LIBRARIES],
dnl   [EXTRA-FLAGS], [PREAMBLE])
dnl --------------------------------------------------------
dnl Search for a library defining FUNC, if it's not already available.
dnl Do not add said library to default LIBS output variable.
dnl Use $ac_lib or $ac_res in ACTION if needed. Uses OTHER-LIBRARIES
dnl unconditionally, which might provoke linker errors.
AC_DEFUN([ACX_OPTION_SEARCH_LIBS],
  [AS_VAR_PUSHDEF([ac_Search], [acx_cv_option_search_$1_]_AC_LANG_ABBREV)dnl
   AC_CACHE_CHECK([for library containing $1], [ac_Search],
     [_ACX_OPTION_SEARCH_LIBS([$1],[$2],[$5],[$6],[$7])])
   ac_res=AS_VAR_GET([ac_Search])
   AS_VAR_SET_IF([ac_Search],
     [$3],
     [$4])
   AS_VAR_POPDEF([ac_Search])dnl
])
dnl ACX_OPTION_SEARCH_LIBS_MULTI(FUNCTION, SEARCH-LIBS,
dnl   [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND], [OTHER-LIBRARIES],
dnl   [EXTRA-FLAGS], [PREAMBLE], [CALL-CODE])
dnl --------------------------------------------------------
dnl Search for a library defining FUNC, if it's not already available.
dnl Do not add said library to default LIBS output variable.
dnl Use $ac_lib or $ac_res in ACTION if needed.
dnl
dnl If CALL-CODE is present it will be compiled to make the function
dnl call, otherwise AC_LANG_CALL will be used with FUNCTION as
dnl argument.
dnl
dnl Tries first to link without any OTHER-LIBRARY, then successively
dnl tries each library in the list.
AC_DEFUN([ACX_OPTION_SEARCH_LIBS_MULTI],
  [AS_VAR_PUSHDEF([ac_Search], [acx_cv_option_search_$1_]_AC_LANG_ABBREV)dnl
   AC_MSG_CHECKING([for library containing $1])
   AC_CACHE_VAL([ac_Search],dnl
     [while :; do
        m4_if(m4_car($5),[[]],,[_ACX_OPTION_SEARCH_LIBS([$1],[$2],,[$6],[$7],[$8])
        AS_VAR_SET_IF([ac_Search], [break])])
        m4_foreach([ACX_LibSet], [$5],dnl
          [_ACX_OPTION_SEARCH_LIBS([$1],[$2],[ACX_LibSet],[$6],[$7],[$8])
           AS_VAR_SET_IF([ac_Search],[break])
          ])
        break
      done])
   AS_VAR_SET_IF([ac_Search],dnl
     [AS_IF([test x"AS_VAR_GET([ac_Search])" = x],dnl
        [AC_MSG_RESULT([(none required)])],dnl
        [AC_MSG_RESULT([AS_VAR_GET([ac_Search])])])],dnl
     [AC_MSG_RESULT([not found])])
   AS_VAR_SET_IF([ac_Search],[$3],[$4])
   AS_VAR_POPDEF([ac_Search])dnl
])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://www.dkrz.de/redmine/projects/show/scales-ppm"
dnl license-default: "bsd"
dnl End:
