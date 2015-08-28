dnl acx_fc_c_link.m4 --- transform library c flags into version
dnl                      that suits the fortran compiler
dnl
dnl Copyright  (C)  2011  Thomas Jahns <jahns@dkrz.de>
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
dnl ACX_FC_XLF_QEXTNAME_ADD_APPENDUS
dnl Test if compiler is xlf and if -qextname is in use.
dnl Add -Dappendus to CPPFLAGS if the above applies.
dnl
AC_DEFUN([ACX_XLF_QEXTNAME_ADD_APPENDUS],
  [AS_CASE([$host],
     [*-ibm-aix*],
     [AC_MSG_CHECKING([if -Dappendus needs to be added to CPPFLAGS for cfortran.h])
      AS_IF([$CC -qversion 2>&1 | grep '^IBM XL C' >/dev/null],
        [acx_temp_qextname_f77flags=`echo "$FFLAGS" | sed -n '/-qextname/{ s/^\(.* \)*-qextname\( .*\)*$/-qextname/;p;}'`
         acx_temp_qextname_fcflags=`echo "$FCFLAGS" | sed -n '/-qextname/{ s/^\(.* \)*-qextname\( .*\)*$/-qextname/;p;}'`
      dnl pretend the same option as for FC was used if F77 isn't used at all
      dnl also in case a non-xl compiler is used it will append an underscore
         AC_PROVIDE_IFELSE([AC_PROG_F77],,
           [AC_PROVIDE_IFELSE([AC_PROG_FC],
              [acx_temp_qextname_f77flags=$acx_temp_qextname_fcflags],
              [m4_fatal([AC_PROG_F77 or AC_PROG_FC must have been invoked prior to ACX_XLF_QEXTNAME_ADD_APPENDUS])])])
      dnl and vice versa for FC
         AC_PROVIDE_IFELSE([AC_PROG_FC],
           [AS_IF([$FC -qversion 2>&1 | grep '^IBM XL Fortran' >/dev/null],,
              [acx_temp_qextname_fcflags=-qextname])],
           [AC_PROVIDE_IFELSE([AC_PROG_F77],
              [acx_temp_qextname_fcflags=$acx_temp_qextname_f77flags])])
         AS_CASE([x"$acx_temp_qextname_fcflags$acx_temp_qextname_f77flags"],
           [x-qextname],
           [AC_MSG_ERROR([Option -qextname must be provided consistently to F77 and FC])],
           [x-qextname-qextname],
           [AC_MSG_RESULT([yes])
            CPPFLAGS="${CPPFLAGS+$CPPFLAGS }-Dappendus"],
           [AC_MSG_RESULT([no])])
        ],[AC_MSG_RESULT([no])])])])
dnl
dnl automate flag elicitation for cfortran.h
AC_DEFUN([_ACX_FIND_CFORTRAN_DEF],
  [AC_REQUIRE([AC_CANONICAL_HOST])
   _AC_FORTRAN_ASSERT
   AC_LANG_CASE([Fortran],[save_FC=$FC ; acx_FC=$FC],
     [Fortran 77],[save_F77=$F77 ; acx_FC=$F77])
   AS_CASE([$host],
     [x86_64-*-linux-*|i*86-*-linux-*|*-apple-darwin*|ia64-*-linux-*|x86_64-*-freebsd*|i*86-*-freebsd*],
     [acx_temp=`$acx_FC -V 2>&1`
      AS_IF([echo "$acx_temp" | grep '^Copyright.*\(The Portland Group\|NVIDIA CORPORATION\)' >/dev/null],
        [AS_VAR_SET([acx_cf_flag],[-DgFortran])],
        [echo "$acx_temp" | grep '^NAG Fortran Compiler Release' >/dev/null],
        [AS_VAR_SET([acx_cf_flag],[-DNAGf90Fortran])],
        [echo "$acx_temp" | grep '^Intel(R) Fortran.*Compiler' >/dev/null],
        [AS_VAR_SET([acx_cf_flag],[-DgFortran])],
        [echo "$acx_temp" | grep '^Cray Fortran' >/dev/null],
        [AS_VAR_SET([acx_cf_flag],[-DgFortran])],
        [acx_temp=`$acx_FC --version 2>&1` \
         && echo $acx_temp | grep '^GNU Fortran' >/dev/null],
        [AS_IF([echo $acx_temp | grep g77 >/dev/null],
           [AS_VAR_SET([acx_cf_flag],[-Dg77Fortran])],
           [dnl assume gfortran
dnl check if compiling with f2c bindings or with default bindings
            AS_IF([echo "]AC_LANG_CASE([Fortran],[$FCFLAGS],
              [Fortran 77],[$FFLAGS])[" | grep '^\(.* \)*-ff2c\( .*\)*$' >/dev/null],
              [AS_VAR_SET([acx_cf_flag],[-Df2cFortran])],
              [AS_VAR_SET([acx_cf_flag],[-DgFortran])])])],
        [acx_temp=`$acx_FC -v 2>&1` \
         && echo $acx_temp | grep '^f2c'],
        [AS_VAR_SET([acx_cf_flag],[-Df2cFortran])])],
     [*-ibm-aix*],
     [dnl xlc set _IBMR2 so nothing needs to be done
      AS_IF([$CC -qversion 2>&1 | grep '^IBM XL C' >/dev/null],,
        [dnl but for other compilers set IBMR2Fortran
         AS_VAR_SET([acx_cf_flag],[-DIBMR2Fortran])])
     ],
     [*-*-hpux*],
     [AS_VAR_SET([acx_cf_flag],[-DhpuxFortran])],
     [sx*-*-*|es*-*-*],
     [dnl fixme: make sure user is actually using sxf90
dnl but currently there is no alternative I know of
      AS_VAR_SET([acx_cf_flag],[-DSXFortran])])])

AC_DEFUN([ACX_FIND_CFORTRAN_DEF],
  [AC_CACHE_CHECK([C preprocessor flags for Fortran calling convention cfortran.h],
     [acx_cv_cf_flag],
     [acx_cv_cf_flag=''
dnl test if user already provided a flag
      AS_FOR([MACRO],[macro],[pgiFortran NAGf90Fortran f2cFortran hpuxFortran apolloFortran sunFortran IBMR2Fortran CRAYFortran PATHSCALE_COMPILER gFortran mipsFortran DECFortran vmsFortran CONVEXFortran PowerStationFortran AbsoftUNIXFortran AbsoftProFortran SXFortran],
        [acx_temp=`echo "$CPPFLAGS $CFLAGS" | sed -n 's/^\(.* \)*-D\('"MACRO"'\)\( .*\)*$/\2/;t print
b
: print
p'`
         AS_IF([test x"$acx_temp" != x],
           [AS_IF([test x"$acx_cv_cf_flag" = x],
              [acx_cv_cf_flag="$acx_temp (user-specified)"],
              [echo ; echo '"'"$acx_cv_cf_flag $acx_temp"'"'
               AC_MSG_ERROR([Multiple specification of cfortran.h flags])])])])
dnl find automatically from machine/compiler
      AS_IF([test x"$acx_cv_cf_flag" = x],
        [AC_PROVIDE_IFELSE([AC_PROG_FC],
           [AS_IF([test -n "$FC" -a X"$FC" != Xno],
              [AC_LANG_PUSH([Fortran])
               AS_VAR_PUSHDEF([acx_cf_flag],[acx_cv_]_AC_LANG_ABBREV[_cf_flag])
               _ACX_FIND_CFORTRAN_DEF
               AS_VAR_POPDEF([acx_cf_flag])
               AC_LANG_POP([Fortran])])])
         AC_PROVIDE_IFELSE([AC_PROG_F77],
           [AS_IF([test -n "$F77" -a X"$F77" != Xno],
              [AC_LANG_PUSH([Fortran 77])
               AS_VAR_PUSHDEF([acx_cf_flag],[acx_cv_]_AC_LANG_ABBREV[_cf_flag])
               _ACX_FIND_CFORTRAN_DEF
               AS_VAR_POPDEF([acx_cf_flag])
               AC_LANG_POP([Fortran 77])])])
dnl check f77 flag matches fc flag
         AC_PROVIDE_IFELSE([AC_PROG_F77],
           [AC_PROVIDE_IFELSE([AC_PROG_FC],
              [dnl both FC and F77 are configured
               AS_IF([test -z "$FC" -o X"$FC" != Xno],
                 [acx_cv_cf_flag="$acx_cv_fc_cf_flag (probed)"],
                 [test -z "$F77" -o X"$F77" != Xno],
                 [acx_cv_cf_flag="$acx_cv_f77_cf_flag (probed)"],
                 [AS_IF([test x"$acx_cv_f77_cf_flag" = x"$acx_cv_fc_cf_flag"],
                    [acx_cv_cf_flag="$acx_cv_f77_cf_flag (probed)"],
                    [AC_MSG_ERROR([cfortran.h flag for $F77 does not match the flag for $FC.
Have you configured compatible compilers?])])])
              ])],[acx_cv_cf_flag="$acx_cv_fc_cf_flag (probed)"])
        ])
     ])
dnl now that flag is established, add (probed) defines to CPPFLAGS
   AS_IF([echo "$acx_cv_cf_flag" | grep ' (probed)$' >/dev/null],
     [CPPFLAGS="${CPPFLAGS+$CPPFLAGS }`echo "$acx_cv_cf_flag" | sed 's/ (probed)$//'`"])
  ])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://www.dkrz.de/redmine/projects/show/scales-ppm"
dnl license-default: "bsd"
dnl End:
