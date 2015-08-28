dnl acx_sl_fc_mod_path_flag.m4 --- find flag to use module from other directory
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
# ACX_SL_FC_MOD_PATH_FLAG([ACTION-IF-FOUND],[ACTION-IF-NOT-FOUND])
# -------------------------
# Check which flag is necessary to alter the compiler's search path
# for module files.
# This obviously requires that the compiler has some notion of
# module files as separate from object files and some sensible
# method of altering its search path. This will therefore not work
# on early Cray F90 compilers, or on v5 (and 6?) of ifc.
#
# Nearly every compiler I have found uses -Ipath for this purpose;
# Sun F95 v7.1 (at least), uses -Mpath
#
AC_DEFUN([ACX_SL_FC_CHECK_MOD_PATH_FLAG],dnl
  [ACX_ASSERT_LANG_IS_FORTRAN_VARIANT
   AC_REQUIRE([AC_PROG_FC])
   AS_VAR_PUSHDEF([mod_flag],[acx_sl_cv_fc_mod_path_flag_]_AC_LANG_ABBREV)dnl
   ASX_VAR_UNSET([mod_flag])
   AC_CACHE_CHECK([for flag to alter module search path],[mod_flag],
     [mkdir conftestdir
      cd conftestdir
      AC_COMPILE_IFELSE([AC_LANG_SOURCE(
[      module cnftst
      implicit none
      integer :: i
      end module cnftst])],,
        [AC_MSG_ERROR([Cannot compile fortran modules])])
      cd ..
      for i in -I -M -module -p; do
        FCFLAGS_save=$FCFLAGS
        FCFLAGS="$FCFLAGS ${i}conftestdir"
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],
[      use cnftst
      implicit none
      i = 0])],
          [AS_VAR_SET([mod_flag],[$i]) ; FCFLAGS=$FCFLAGS_save ; break],
          [:])
        FCFLAGS=$FCFLAGS_save
      done
      FCFLAGS=$FCFLAGS_save
      rm -rf conftestdir
      AS_VAR_SET_IF([mod_flag],dnl
        [m4_default([$1],[:])],dnl
        [m4_default([$2],[AC_MSG_ERROR([Cannot find flag to alter module search path])])])])
   FC_MOD_FLAG=AS_VAR_GET([mod_flag])
   AC_SUBST([FC_MOD_FLAG])
   AS_VAR_POPDEF([mod_flag])dnl
])# ACX_SL_FC_MOD_PATH_FLAG
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://www.dkrz.de/redmine/projects/show/scales-ppm"
dnl license-default: "bsd"
dnl End:
