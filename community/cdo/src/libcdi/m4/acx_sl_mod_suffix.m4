dnl acx_sl_mod_suffix.m4 --- determine name of symbol file for Fortran 90 module
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
# ACX_SL_FC_MOD_SUFFIX([ACTION-IF-FOUND],[ACTION-IF-NOT-FOUND])
# -----------------
# Determines the form of the filename of modules produced
# by the Fortran compiler.
# Tests for all forms of file extension I've (TOHW) found in the
# wild. Note that at least one compiler (PGI??) changes the
# case of the basename as well. Whether this happens is
# encoded in the variable ac_fc_mod_uppercase.
#
# This macro depends, of course, on the Fortran compiler producing
# module files. See comment to AC_FC_MOD_PATH_FLAG.
#
# FIXME: This will fail if an F77-only compiler is used.
# Currently we warn and continue. We should maybe error out.
#
AC_DEFUN([ACX_SL_FC_MOD_SUFFIX],
  [AC_MSG_CHECKING([for suffix of module files])
   AC_ARG_VAR([FCMODEXT], [file extension of compiled Fortran module files])
   ac_fc_mod_uppercase=no
   AC_LANG_PUSH([Fortran])
   AC_COMPILE_IFELSE(
[      module conftest
      implicit none
      integer :: i
      end module conftest])
   while :; do
     acx_fc_mod_name=
   m4_foreach([acx_fc_mod_name],dnl
     [[conftest.$FCMODEXT], [conftest.mod], [conftest.MOD], [conftest.M],
      [CONFTEST.MOD]],dnl
     [AS_IF([test -n "acx_fc_mod_name" -a -f "acx_fc_mod_name"],dnl
        [[acx_fc_mod_name]="acx_fc_mod_name" ; break])
     ])
     break
   done
   rm -f conftest*
   AC_LANG_POP([Fortran])
dnl
   AS_CASE(["$acx_fc_mod_name"],dnl
     [conftest.$FCMODEXT], [:],
     [CONFTEST.$FCMODEXT], [ac_fc_mod_uppercase=yes],
     [conftest.mod], [FCMODEXT=mod],
     [conftest.MOD], [FCMODEXT=MOD],
     [conftest.M], [FCMODEXT=M],
     [CONFTEST.MOD], [FCMODEXT=MOD
        ac_fc_mod_uppercase=yes])
   AC_MSG_RESULT([${FCMODEXT-not found}])
   AS_VAR_TEST_SET([FCMODEXT], [$1], [m4_ifval([$2],[$2],dnl
     [AC_MSG_WARN([Could not find Fortran module file extension.])])])
dnl
   AS_IF([test $ac_fc_mod_uppercase = yes],
     [FCMODCASE=uc
      AC_MSG_NOTICE([Fortran module filenames are uppercase.])],
     [FCMODCASE=lc])
   AC_SUBST([FCMODEXT])
   AC_SUBST([FCMODCASE])
])dnl _ACX_SL_FC_MOD_SUFFIX
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://www.dkrz.de/redmine/projects/show/scales-ppm"
dnl license-default: "bsd"
dnl End:
