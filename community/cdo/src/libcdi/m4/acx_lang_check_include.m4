dnl acx_lang_check_include.m4 --- generic check for package includes
dnl
dnl Copyright  (C)  2010  Thomas Jahns <jahns@dkrz.de>
dnl
dnl Version: 1.0
dnl Keywords: package include check
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
# ACX_LANG_CHECK_INCLUDE_PATHS_IFELSE(HEADER-FILE, PATH...
#   [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND],
#   [INCLUDES], [INCLUDE-FLAGS])
# ------------------------------------------------------------------
# Check wether the compiler accepts HEADER-FILE.
# The INCLUDES might be defaulted.
AC_DEFUN([ACX_LANG_CHECK_INCLUDE_PATHS_IFELSE],dnl
  [_AC_LANG_DISPATCH([$0],_AC_LANG,$@)])

# ACX_LANG_INCLUDE_PROGRAM(HEADER-FILE,[PRELUDE])
AC_DEFUN([ACX_LANG_INCLUDE_PROGRAM],
  [_AC_LANG_DISPATCH([$0],_AC_LANG,$@)])

# AC_LANG_CHECK_INCLUDE_IFELSE(HEADER-FILE,
#   [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND],
#   [INCLUDES])
# ------------------------------------------------------------------
# Check the compiler accepts HEADER-FILE. The INCLUDES might be defaulted.
AC_DEFUN([_ACX_LANG_CHECK_INCLUDE_IFELSE],
  [AC_COMPILE_IFELSE([ACX_LANG_INCLUDE_PROGRAM([$1],[$4])],
     [AS_VAR_SET([acx_Include], [yes])])
   AS_VAR_SET_IF([acx_Include],[$2],[$3])dnl
])# AC_LANG_CHECK_INCLUDE_IFELSE

# ACX_LANG_CHECK_INCLUDE_IFELSE(HEADER-FILE,
#   [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND],
#   [INCLUDES])
# ------------------------------------------------------------------
# Check the compiler accepts HEADER-FILE. The INCLUDES might be defaulted.
AC_DEFUN([ACX_LANG_CHECK_INCLUDE_IFELSE],
  [AS_VAR_PUSHDEF([acx_Include], [acx_cv_]_AC_LANG_ABBREV[_include_$1])dnl
   AC_CACHE_CHECK([for $1],[acx_Include],dnl
     [_ACX_LANG_CHECK_INCLUDE_IFELSE([$1],[$2],[$3],[$4])])
   AS_VAR_POPDEF([acx_Include])])

# ACX_GENERIC_CHECK_INCLUDE_PATHS_IFELSE(FLAGSVAR, INCOPT, HEADER-FILE, PATH...
#   [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND],
#   [INCLUDES], [INCLUDE-FLAGS], [EXTRA-INCLUDE-FLAG-SETS])
# ------------------------------------------------------------------
# Check wether the the current language compiler accepts HEADER-FILE.
# The INCLUDES might be defaulted.
AC_DEFUN([ACX_GENERIC_CHECK_INCLUDE_PATHS_IFELSE],
  [AS_VAR_PUSHDEF([acx_Include], [acx_cv_]_AC_LANG_ABBREV[_include_$3])dnl
   AS_VAR_PUSHDEF([save_flags],[ac_include_search_$3_SAVE])
   AC_MSG_CHECKING([for $3 extra include path])
   AC_CACHE_VAL([acx_Include],dnl
     [AS_VAR_SET([save_flags],["@S|@$1"])
      while :; do
        m4_foreach([ACX_IncSet],[$9],
          [for ac_incdir in ''m4_ifval([$4],[ $4]); do
             AS_IF([test -z "$ac_incdir"],
               [ac_res="none required"
                $1="m4_ifval(ACX_IncSet,ACX_IncSet )$8 AS_VAR_GET([save_flags])"],
               [ac_res="$2$ac_incdir"
                $1="m4_ifval(ACX_IncSet,ACX_IncSet )$8 $ac_res AS_VAR_GET([save_flags])"])
             _ACX_LANG_CHECK_INCLUDE_IFELSE([$3],dnl
               [AS_IF([test -z "$ac_incdir"],dnl
                  [AS_VAR_SET([acx_Include],["]m4_ifval(ACX_IncSet,ACX_IncSet )$8["])],dnl
                  [AS_VAR_SET([acx_Include],["]m4_ifval(ACX_IncSet,ACX_IncSet )[$8 $2$ac_incdir"])])],,[$7])
             AS_VAR_SET_IF([acx_Include], [break])
           done
           AS_VAR_SET_IF([acx_Include], [break])
          ])
        break
      done
      AS_VAR_COPY([$1],[save_flags])])
   AS_VAR_SET_IF([acx_Include],
     [acx_temp=]AS_VAR_GET([acx_Include])[
      AS_IF([test x"$acx_temp" = x],
        [AC_MSG_RESULT([(none required)])],
        [AC_MSG_RESULT([$acx_temp])])],
     [AC_MSG_RESULT([not found])])
   AS_VAR_SET_IF([acx_Include], [$5], [$6])
   AS_VAR_POPDEF([acx_Include])])dnl
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://www.dkrz.de/redmine/projects/show/scales-ppm"
dnl license-default: "bsd"
dnl End:
