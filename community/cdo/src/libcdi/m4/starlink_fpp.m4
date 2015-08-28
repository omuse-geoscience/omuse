# This file is part of Autoconf.                       -*- Autoconf -*-
# Fortran languages support.
# Copyright (C) 2001, 2003-2005
# Free Software Foundation, Inc.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.
#
# As a special exception, the Free Software Foundation gives unlimited
# permission to copy, distribute and modify the configure scripts that
# are the output of Autoconf.  You need not follow the terms of the GNU
# General Public License when using or distributing such scripts, even
# though portions of the text of Autoconf appear in them.  The GNU
# General Public License (GPL) does govern all other use of the material
# that constitutes the Autoconf program.
#
# Certain portions of the Autoconf source text are designed to be copied
# (in certain cases, depending on the input) into the output of
# Autoconf.  We call these the "data" portions.  The rest of the Autoconf
# source text consists of comments plus executable code that decides which
# of the data portions to output in any given case.  We call these
# comments and executable code the "non-data" portions.  Autoconf never
# copies any of the non-data portions into its output.
#
# This special exception to the GPL applies to versions of Autoconf
# released by the Free Software Foundation.  When you make and
# distribute a modified version of Autoconf, you may extend this special
# exception to the GPL to apply to your modified version as well, *unless*
# your modified version has the potential to copy into its output some
# of the text that was the non-data portion of the version that you started
# with.  (In other words, unless your change moves or copies text from
# the non-data portions to the data portions.)  If your modification has
# such potential, you must delete any notice of this special exception
# to the GPL from your modified version.
#
# Written by David MacKenzie, with help from
# Franc,ois Pinard, Karl Berry, Richard Pixley, Ian Lance Taylor,
# Roland McGrath, Noah Friedman, david d zuhn, and many others.
#
# Fortran preprocessing support written by Martin Wilck, adapted and
# extended by Norman Gray and Toby White.
#
# Reduced to searching the Fortran preprocessor and flags for FC to
# preprocess by Thomas Jahns, 2009


## ------------------------------- ##
## Preprocessor features           ##
## ------------------------------- ##
#
# Supported features are:
#
# include   : correctly process #include directives and -I
# define    : correctly process -D
# substitute: substitute macros in Fortran code
#             (some preprocessors touch only lines starting with #)
# wrap      : wrap lines that become too long through macro substitution
#             fpp is probably the only preprocessor that does this.
# cstyle    : Do not suppress C style comments (-C option in cpp)
# CSTYLE    : *Do* suppress C style comments
#             (e.g. code contains C-style comments, and compiler may not
#             know how to handle them)
# cxxstyle  : Do not suppress C++ style comments (default)
# CXXSTYLE  : *Do* suppress C++ style comments (seems unlikely, but in here
#             for completeness)
#
# Features can be abbreviated: i, in, inc etc. are equivalent to include.
# Features can be deselected (feature not needed) by prepending "no",
#   e.g. nodef (=nodefine), now (=nowrap).
#
# Default for the feature list is
#       [include define substitute nowrap nocstyle noCSTYLE cxxstyle]
# Feature requirements corresponding to the defaults may be omitted
#
# Note that "wrap" implies "substitute", and CSTYLE and cstyle cannot
# be requested at the same time. The macro adjusts this automatically.
#

# -------------------------------------- #
# Feature tests for Preprocessed Fortran #
# -------------------------------------- #

# ------------------#
#   Internal macros #
# ------------------#

# ----------------------------------------------#
#     Some test programs for different features #
# ----------------------------------------------#

# _AC_LANG_PROGRAM_FPP_SIMPLE
# ---------------------------
# The minimum test program - any compiler supporting
# preprocessing should handle this
AC_DEFUN([_ACX_SL_LANG_PROGRAM_FPP_SIMPLE],
         [AC_LANG_PROGRAM(,[@%:@define OK
@%:@ifndef OK
      syntax error
@%:@endif])])#_ACX_SL_LANG_PROGRAM_FPP_SIMPLE


# _ACX_SL_LANG_PROGRAM_FPP_ONLY
# ---------------------------
# Test program for pure preprocessing
# Note that other macros test for literal strings within this, so check
# for those if you have to change anything here.
AC_DEFUN([_ACX_SL_LANG_PROGRAM_FPP_ONLY],
         [AC_LANG_PROGRAM(,[@%:@define OK
@%:@ifdef OK
      REAL A
@%:@else
      syntax error
@%:@endif])])#_ACX_SL_LANG_PROGRAM_FPP_ONLY


# _ACX_SL_LANG_PROGRAM_FPP_D
# ---------------------------
# Like _ACX_SL_LANG_PROGRAM_FPP_SIMPLE, but OK is passed via -D switch
AC_DEFUN([_ACX_SL_LANG_PROGRAM_FPP_D],
[AC_LANG_PROGRAM([],[@%:@ifndef OK
      syntax error
@%:@endif])])#_ACX_SL_LANG_PROGRAM_FPP_D


# _ACX_SL_LANG_PROGRAM_FPP_I
# ---------------------------
# Test for #include statement
# If unsupported, this should give a type error
AC_DEFUN([_ACX_SL_LANG_PROGRAM_FPP_I],
  [AC_LANG_PROGRAM(,[      IMPLICIT CHARACTER (c)
       !     Comments in test programs should be freeform compliant just in case.
       !     conftest.inc contains the Fortran statement "REAL cc"
@%:@include "conftest.inc"
      cc=1.])])#_ACX_SL_LANG_PROGRAM_FPP_I


# _ACX_SL_LANG_PROGRAM_FPP_SUBS
# ---------------------------
# Test whether cpp symbols are expanded in Fortran code lines
# If not, this should give a type error
AC_DEFUN([_ACX_SL_LANG_PROGRAM_FPP_SUBS],
[AC_LANG_PROGRAM(,[@%:@define NM xxxx
      IMPLICIT CHARACTER (n)
      REAL xxxx
      NM=1.])])#_ACX_SL_LANG_PROGRAM_FPP_SUBS


# _ACX_SL_LANG_PROGRAM_FPP_WRAP
# ---------------------------
# Test whether preprocessor breaks lines that become too long due
# to macro substitution.
# If not, this gives an "unterminated character constant" error
AC_DEFUN([_ACX_SL_LANG_PROGRAM_FPP_WRAP],
  [AC_LANG_PROGRAM(,
     [m4_case(_AC_LANG,
        [Fortran],
[@%:@define LONG '901234567890123456789012345678901234567890123456789012345678901234567890'
      CHARACTER(LEN=80) :: A
      A=LONG],
        [Fortran 77],
[@%:@define LONG '901234567890123456789012345678901234567890123456789012345678901234567890'
      CHARACTER*80 A
      A=LONG],
        [m4_fatal([$0: current language is not Fortran: ] _AC_LANG)])])])#_ACX_SL_LANG_PROGRAM_FPP_WRAP


# _ACX_SL_LANG_PROGRAM_FPP_CSTYLE
# ---------------------------
# Test program for C style comments
AC_DEFUN([_ACX_SL_LANG_PROGRAM_FPP_CSTYLE],
[AC_LANG_PROGRAM(,[      A=1. /* C-style comment */])])#_ACX_SL_LANG_PROGRAM_FPP_CSTYLE

# _ACX_SL_LANG_PROGRAM_FPP_CXXSTYLE
# ---------------------------
# Test program for C++ style comments
AC_DEFUN([_ACX_SL_LANG_PROGRAM_FPP_CXXSTYLE],
  [AC_LANG_SOURCE([m4_case(_AC_LANG,
        [Fortran],
[      PROGRAM MAIN
      CHARACTER(LEN=10) :: C
      C = "abcde" // "fghij"; END PROGRAM],
        [Fortran 77],
[      PROGRAM MAIN
      CHARACTER*10 C
      C = "abcde" // "fghij"; END PROGRAM])])])#_ACX_SL_LANG_PROGRAM_FPP_CXXSTYLE

# _ACX_SL_SET_FPP_FEATURE_VARS ([feature list])
# --------------------------------------
# Parse the feature list from configure.in
AC_DEFUN([_ACX_SL_SET_FPP_FEATURE_VARS],
  [# defaults for needed features
   ac_fpp_need_d=yes
   ac_fpp_need_i=yes
   ac_fpp_need_subs=no
   ac_fpp_need_wrap=no
   ac_fpp_need_cstyle=no
   ac_fpp_need_CSTYLE=no
   ac_fpp_need_cxxstyle=yes
   ac_fpp_need_CXXSTYLE=no
   dnl FIXME: this should be feasable within m4 constructs, i.e. without
   dnl using shell loops
   for _t in $1 nil
   do
     AS_CASE([$_t],
       [define], [ac_fpp_need_d=yes],
       [nodefine], [ac_fpp_need_d=no],
       [include], [ac_fpp_need_i=yes],
       [noinclude], [ac_fpp_need_i=no],
       [substitute], [ac_fpp_need_subs=yes],
       [nosubstitute], [ac_fpp_need_subs=no],
       [wrap], [ac_fpp_need_wrap=yes],
       [nowwrap], [ac_fpp_need_wrap=no],
       [cstyle], [ac_fpp_need_cstyle=yes],
       [nocstyle], [ac_fpp_need_cstyle=no],
       [CSTYLE], [ac_fpp_need_CSTYLE=yes],
       [noCSTYLE], [ac_fpp_need_CSTYLE=no],
       [cxxstyle], [ac_fpp_need_cxxstyle=yes],
       [nocxxstyle], [ac_fpp_need_cxxstyle=no],
       [CXXSTYLE], [ac_fpp_need_CXXSTYLE=yes],
       [noCXXSTYLE], [ac_fpp_need_CXXSTYLE=no],
       [nil], [])
   done
   # Wrapping requires substitution
   test $ac_fpp_need_wrap = yes && ac_fpp_need_subs=yes
   # CSTYLE and cstyle are mutually exclusive.
   # CSTYLE takes precedence, since if it is not fulfilled,
   # compile errors may arise
   test $ac_fpp_need_CSTYLE = yes && ac_fpp_need_cstyle=no
   dnl Similarly for cxxstyle
   test $ac_fpp_need_CXXSTYLE = yes && ac_fpp_need_cxxstyle=no
  ])# _ACX_SL_SET_FPP_FEATURE_VARS


# _ACX_SL_TEST_FPP(COMMAND,SUFFIX,[ACTION-IF-SUCCESSFULL],[ACTION-IF-FAILED])
# ------------------------
# A helper macro to test correct fpp behaviour
# Invokes COMMAND <FILE_ONLY_COMPILABLE_WHEN_PREPROC_RUN>
# If the output on stdout is valid input to the Fortran compiler, it sets
#  * acx_sl_prog_fpp
#  * acx_sl_prog_fpp_suffix
# accordingly.
AC_DEFUN([_ACX_SL_TEST_FPP],
  [rm -f conftest*
   ACX_ASSERT_LANG_IS_FORTRAN_VARIANT
   AC_LANG_CONFTEST([_ACX_SL_LANG_PROGRAM_FPP_ONLY])
   ac_tmp=conftest.fppout
   AS_IF([_AC_RUN_LOG([$1 conftest.$2 >$ac_tmp],dnl
       [_AC_DO_ECHO([$1 $FPPFLAGS conftest.$2 >$ac_tmp])])],
     [AS_IF([test -f $ac_tmp \
        && ! cmp conftest.$2 $ac_tmp >/dev/null 2>&1 \
        && grep '^      REAL A' $ac_tmp >/dev/null 2>&1 \
        && ! grep 'syntax error' $ac_tmp >/dev/null 2>&1],
        [# we have Fortran!  See if the file can be compiled:
         mv $ac_tmp conftest.$ac_ext
         AC_COMPILE_IFELSE(, [_AC_MSG_LOG_CONFTEST
            $3],
           [_AC_MSG_LOG_CONFTEST
            $4])],
        [mv $ac_tmp conftest.$ac_ext
         _AC_MSG_LOG_CONFTEST
         $4])])
dnl Preprocessing might fail for one of the following reasons:
dnl * no output was produced (disk full?),
dnl * FPP input and output have the same content,
dnl * the output did not contain the critical part of the source file
dnl * some syntax error crept into the output
dnl indicating that this is a case-insensitive filesystem or
dnl preprocessing failed. So this command doesn't work.
   rm -f conftest*
])# _ACX_SL_TEST_FPP

# _ACX_SL_PROG_FPP([SUFFIX], [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# ------------
# Try to figure out how to preprocess files with the given suffix
# just like the selected Fortran compiler does
#
# Must be run after _ACX_SL_PROG_FC_CPP
AC_DEFUN([_ACX_SL_PROG_FPP],
  [acx_sl_fpp_srcext=m4_default([$1],[${ac_fc_srcext-F}])
   AS_VAR_PUSHDEF([acx_sl_prog_fpp],
     [acx_sl_cv_prog_fpp_]_AC_LANG_ABBREV[_]m4_default([$1],[${ac_fc_srcext-f}]))
   AC_CACHE_CHECK([how to preprocess Fortran files with suffix $acx_sl_fpp_srcext],
     [acx_sl_prog_fpp],
     [AS_VAR_SET([acx_sl_prog_fpp], [])
      # Let the user specify FPP
      AS_IF([test -n "$FPP"],
        [_ACX_SL_TEST_FPP([$FPP],[$acx_sl_fpp_srcext],
           [AS_VAR_SET([acx_sl_prog_fpp], ["$FPP"])],
           [AC_MSG_WARN([user-specified \$FPP ($FPP) does not work])
            FPP=])])
      AS_IF([test -z "$FPP"],
        [for ac_fpp in `cd $srcdir ; pwd`/util/sxpreproc-wrapper \
           `cd $srcdir ; pwd`/util/xlfpreproc-wrapper \
           `cd $srcdir ; pwd`/util/sunf95preproc-wrapper \
           `cd $srcdir ; pwd`/util/crayftnpreproc-wrapper \
           "$FC -F" "$FC -F -fpp" "$FC -E" "$FC -E -cpp" \
           "$FC $FCFLAGS -F" "$FC $FCFLAGS -E" "$FC $FCFLAGS -E" \
           "$FC $FCFLAGS -E -cpp" "$FC $FCFLAGS -x f95-cpp-input -E -P"
         do
           _ACX_SL_TEST_FPP([$ac_fpp],[$acx_sl_fpp_srcext],[FPP="$ac_fpp"
              break])
           _ACX_SL_TEST_FPP([$ac_fpp -P],[$acx_sl_fpp_srcext],
             [FPP="$ac_fpp -P"
              break])
         done
         dnl the above tests might generate an executable
         \rm a.out 2>/dev/null])
      AS_IF([test -z "$FPP"], [$3],
        [AS_VAR_SET([acx_sl_prog_fpp], [$FPP])
         $2])])
   AS_VAR_PUSHDEF([acx_sl_prog_fpp])
])# _ACX_SL_PROG_FPP



# _ACX_SL_PROG_FPP_CSTYLE
# -------------------
# Check whether FPP lets C-style comments through to FC
AC_DEFUN([_ACX_SL_PROG_FPP_CSTYLE],
[AC_CACHE_CHECK([how to pass C-style comments to $FC],
   ac_cv_prog_fpp_cstyle,
[ac_cv_prog_fpp_cstyle=unknown
ACX_ASSERT_LANG_IS_FORTRAN_VARIANT
cat > conftest.$ac_ext << \_ACEOF
_ACX_SL_LANG_PROGRAM_FPP_CSTYLE
_ACEOF

AC_LANG_PUSH(Fortran)
ac_cmd='$FPP $FPPFLAGS conftest.$ac_ext '"$ac_fpp_out"
AS_IF([AC_TRY_EVAL(ac_cmd) && \
   cat conftest.f | grep '[[*]]/.*[[*]]/' >/dev/null 2>&1],
   [ac_cv_prog_fpp_cstyle=],
   [ac_save_FPPFLAGS=$FPPFLAGS
    ac_name=`expr "x$FPP" : 'x\(fpp\)'`
    AS_IF([test "x$ac_name" = xfpp],
      [ac_flag="-c_com=no"],
      [ac_flag="-C"])
   FPPFLAGS="$FPPFLAGS $ac_flag"
   ac_cmd='$FPP $FPPFLAGS conftest.$ac_ext '"$ac_fpp_out"
   AS_IF([AC_TRY_EVAL(ac_cmd) && \
     cat conftest.f | grep '/[[*]].*[[*]]/' >/dev/null 2>&1],
     [ac_cv_prog_fpp_cstyle=$ac_flag])
   FPPFLAGS=$ac_save_FPPFLAGS])
rm -f conftest*
AC_LANG_POP(Fortran)dnl
])
if test "x$ac_cv_prog_fpp_cstyle" = "xunknown"; then
  AC_MSG_WARN([cannot find a way to make $FPP pass C-style comments])
else
  FPPFLAGS="$FPPFLAGS $ac_cv_prog_fpp_cstyle"
fi
])# _ACX_SL_PROG_FPP_CSTYLE

# _ACX_CHECK_FPP_COMPILE([EXTRA-FLAGS], direct|indirect,
#    [ACTION-IF-COMPILABLE],[ACTION-IF-NOT-COMPILABLE])
# helper macro to run the preprocessor or compile directly
# with preprocessor flags
m4_define([_ACX_FPP_COMPILE_IFELSE],
  [m4_if([$2],[direct],
     [cp conftest.$acx_sl_fpp_srcext conftest.${ac_ext}.tmp],
     [_AC_RUN_LOG([$FPP $FPPFLAGS m4_ifval([$1],[$1 ])conftest.$acx_sl_fpp_srcext \
        >conftest.${ac_ext}.tmp],
        [echo Running preprocessor $FPP $FPPFLAGS m4_ifval([$1],[$1 ])conftest.$acx_sl_fpp_srcext])])
   m4_if([$2],[direct],[FCFLAGS="$FPPFLAGS $1 $ac_save_FCFLAGS"])
   mv conftest.${ac_ext}.tmp conftest.${ac_ext}
   AC_COMPILE_IFELSE(,[$3],[$4])])

# ACX_SL_PROG_FC_FPP_FEATURES(FEATURES,METHOD,[SUFFIX],
#  [ACTION-IF-SUCCESS],[ACTION-IF-FAILURE])
# ---------------
# This macro checks whether the chosen preprocessing method
# has all the requested features.
#
# This macro must be called with METHOD set to either
# direct or indirect; it behaves differently accordingly.
#
# In case one of the checks fails, ACTION-IF-FAILED will be called.
#
#FIXME: this is only for fixed form code. Need a separate check for free-form.
#
# NB We are definitely using a suffix of .F in this case. If the filesystem
# is case-insensitive, we may need to force preprocessing.
#
# Sets ac_fpp_ok to "no" if a requested feature is unavailable
#
AC_DEFUN([ACX_SL_PROG_FC_FPP_FEATURES],
  [AS_VAR_PUSHDEF([acx_sl_fpp_ok], [acx_sl_cv_]_AC_LANG_ABBREV[_$2_ok])
   AC_CACHE_CHECK([if ]m4_case([$2],[direct],[compilation with $FC],
[indirect],[preprocessing with $FPP],
[m4_fatal([$0: only direct or indirect method supported: ']$2['])])[ of Fortran source supports required preprocessor features],
[acx_sl_fpp_ok],
     [AS_VAR_SET([acx_sl_fpp_ok], [yes])
      # Set up ac_fpp_need_* flags based on features in $1
      acx_sl_fpp_srcext=m4_default([$3],[${ac_fc_srcext-F}])
      _ACX_SL_SET_FPP_FEATURE_VARS([$1])

      acx_sl_prog_fc_cpp_CSTYLE=no
      acx_sl_prog_fc_cpp_cxxstyle=no

      ACX_ASSERT_LANG_IS_FORTRAN_VARIANT
      ac_save_FCFLAGS=$FCFLAGS

      m4_case([$2],[direct],[],[indirect],[],
         [m4_fatal([$0: only direct or indirect method supported: ] $2)])

      # We need to skip the following tests if we're trying direct compilation
      # and FC won't preprocess.
      AS_IF([test x$ac_fpp_need_d = xyes],
        [acx_sl_prog_fc_cpp_d=no
         _AS_ECHO_LOG([Trying flag to create preprocessor defines.])
         ACX_LANG_OTHER_SUFFIX_CONFTEST([$acx_sl_fpp_srcext],
           [_ACX_SL_LANG_PROGRAM_FPP_D])
         # Normally we expect to be able to define preprocessor macros
         # with -D, but this might be IBM xlf compiler, which needs
         # -WF,-D or Fujitsu VPP700 which needs -Wp,-D
         mv conftest.$acx_sl_fpp_srcext conftest.${acx_sl_fpp_srcext}.bak
         for FPP_DEFOPT in "$FPP_DEFOPT" '-D' '-WF,-D' '-Wp,-D'
         do
           AS_IF([test x"$FPP_DEFOPT" != x ],
             [cp conftest.${acx_sl_fpp_srcext}.bak conftest.$acx_sl_fpp_srcext
              _ACX_FPP_COMPILE_IFELSE([${FPP_DEFOPT}OK],[$2],
                [acx_sl_prog_fc_cpp_d=yes; break])])
         done
         FCFLAGS=$ac_save_FCFLAGS
         AS_IF([test $acx_sl_prog_fc_cpp_d = no],
           [AS_VAR_SET([acx_sl_fpp_ok], [no])])
         AS_IF([test x[]AS_VAR_GET([acx_sl_fpp_ok]) = xyes],
           [_AS_ECHO_LOG([Test successful.])],
           [_AS_ECHO_LOG([Test failed.])])
        ])

      AS_IF([test $ac_fpp_need_i = yes],
        [acx_sl_prog_fc_cpp_i=no
         _AS_ECHO_LOG([Trying flag to add directories to preprocessor search path.])
         AS_MKDIR_P([conftst])
         cd conftst
         ACX_LANG_OTHER_SUFFIX_CONFTEST([inc],
           [AC_LANG_SOURCE([       !     This statement overrides the IMPLICIT statement in the program
      REAL cc])])
         cd ..
         ACX_LANG_OTHER_SUFFIX_CONFTEST([$acx_sl_fpp_srcext],[_ACX_SL_LANG_PROGRAM_FPP_I])
         mv conftest.$acx_sl_fpp_srcext conftest.${acx_sl_fpp_srcext}.bak
         for FPP_INCOPT in "$FPP_INCOPT" '-I' '-WF,-I' '-Wp,-I'
         do
           AS_IF([test x"$FPP_INCOPT" != x ],
             [cp conftest.${acx_sl_fpp_srcext}.bak conftest.$acx_sl_fpp_srcext
              _ACX_FPP_COMPILE_IFELSE([${FPP_INCOPT}conftst],[$2],
                [acx_sl_prog_fc_cpp_i=yes
                 break])])
         done
         FCFLAGS=$ac_save_FCFLAGS
         rm -rf conftst
         AS_IF([test $acx_sl_prog_fc_cpp_i = no],
           [AS_VAR_SET([acx_sl_fpp_ok], [no])])
         AS_IF([test x[]AS_VAR_GET([acx_sl_fpp_ok]) = xyes],
           [_AS_ECHO_LOG([Test successful.])],
           [_AS_ECHO_LOG([Test failed.])])])
dnl
dnl
      AS_IF([test $ac_fpp_need_subs = yes],
        [acx_sl_prog_fc_cpp_subs=no
         _AS_ECHO_LOG([Testing preprocessor expansion in Fortran code.])
         ACX_LANG_OTHER_SUFFIX_CONFTEST([$acx_sl_fpp_srcext],
           [_ACX_SL_LANG_PROGRAM_FPP_SUBS])
         _ACX_FPP_COMPILE_IFELSE(,[$2],
           [acx_sl_prog_fc_cpp_subs=yes],[AS_VAR_SET([acx_sl_fpp_ok], [no])])
         AS_IF([test x[]AS_VAR_GET([acx_sl_fpp_ok]) = xyes],
           [_AS_ECHO_LOG([Test successful.])],
           [_AS_ECHO_LOG([Test failed.])])
        ])
dnl
dnl
      AS_IF([test $ac_fpp_need_wrap = yes],
        [acx_sl_prog_fc_cpp_wrap=no
         _AS_ECHO_LOG([Testing wether preprocessor wraps long lines.])
         ACX_LANG_OTHER_SUFFIX_CONFTEST([$acx_sl_fpp_srcext],
           [_ACX_SL_LANG_PROGRAM_FPP_WRAP])
         _ACX_FPP_COMPILE_IFELSE(,[$2],
           [acx_sl_prog_fc_cpp_wrap=yes], [AS_VAR_SET([acx_sl_fpp_ok], [no])])
         AS_IF([test x[]AS_VAR_GET([acx_sl_fpp_ok]) = xyes],
           [_AS_ECHO_LOG([Test successful.])],
           [_AS_ECHO_LOG([Test failed.])])
        ])
dnl
dnl
      AS_IF([test $ac_fpp_need_CSTYLE = yes],
        [_AS_ECHO_LOG([Testing wether preprocessor removes C-style comments.])
         ACX_LANG_OTHER_SUFFIX_CONFTEST([$acx_sl_fpp_srcext],
           [_ACX_SL_LANG_PROGRAM_FPP_CSTYLE])
         _ACX_FPP_COMPILE_IFELSE(,[$2],
           [acx_sl_prog_fc_cpp_CSTYLE=yes], [AS_VAR_SET([acx_sl_fpp_ok], [no])])
         AS_IF([test x[]AS_VAR_GET([acx_sl_fpp_ok]) = xyes],
           [_AS_ECHO_LOG([Test successful.])],
           [_AS_ECHO_LOG([Test failed.])])
        ])
dnl
      AS_IF([test $ac_fpp_need_cstyle = yes],
        [_AS_ECHO_LOG([Testing wether preprocessor leaves C-style comments in place.])
         ACX_LANG_OTHER_SUFFIX_CONFTEST([$acx_sl_fpp_srcext],
           [_ACX_SL_LANG_PROGRAM_FPP_CSTYLE])
         _ACX_FPP_COMPILE_IFELSE(,[$2],
           [acx_sl_prog_fc_cpp_CSTYLE=yes], [AS_VAR_SET([acx_sl_fpp_ok], [no])])
         AS_IF([test x[]AS_VAR_GET([acx_sl_fpp_ok]) = xyes],
           [_AS_ECHO_LOG([Test failed.])],
           [_AS_ECHO_LOG([Test successful.])])
        ])
dnl
      AS_IF([test $ac_fpp_need_cxxstyle = yes],
        [_AS_ECHO_LOG([Testing if preprocessor leaves C++-style comments in place.])
         ACX_LANG_OTHER_SUFFIX_CONFTEST([$acx_sl_fpp_srcext],
           [_ACX_SL_LANG_PROGRAM_FPP_CXXSTYLE])
         _ACX_FPP_COMPILE_IFELSE(,[$2],
           [acx_sl_prog_fc_cpp_cxxstyle=yes],[AS_VAR_SET([acx_sl_fpp_ok], [no])])
         AS_IF([test x[]AS_VAR_GET([acx_sl_fpp_ok]) = xyes],
           [_AS_ECHO_LOG([Test successful.])],
           [_AS_ECHO_LOG([Test failed.])])
        ])
dnl
      AS_IF([test $ac_fpp_need_CXXSTYLE = yes],
        [_AS_ECHO_LOG([Testing if preprocessor suppresses C++-style comments.])
         ACX_LANG_OTHER_SUFFIX_CONFTEST([$acx_sl_fpp_srcext],
           [_ACX_SL_LANG_PROGRAM_FPP_CXXSTYLE])
         _ACX_FPP_COMPILE_IFELSE(,[$2],
           [acx_sl_prog_fc_cpp_cxxstyle=yes],
           [AS_VAR_SET([acx_sl_fpp_ok], [no])])
         AS_IF([test x[]AS_VAR_GET([acx_sl_fpp_ok]) = xyes],
           [_AS_ECHO_LOG([Test failed.])],
           [_AS_ECHO_LOG([Test successful.])])
        ])
dnl
      FCFLAGS=$ac_save_FCFLAGS
      rm -f conftest.*
      AS_IF([test x[]AS_VAR_GET([acx_sl_fpp_ok]) = xyes],[$4],[$5])
     ])
   #FIXME we should probably do the AC_SUBST somewhere else.
   AC_SUBST([FPP_DEFOPT])
   AC_SUBST([FPP_INCOPT])
   AS_VAR_POPDEF([acx_sl_fpp_ok])
 ])#ACX_SL_PROG_FC_FPP_FEATURES



# _ACX_SL_FC_CHECK_CIFS
# -----------------
# Check whether the filesystem is case-insensitive (eg, HFS+ on
# MacOS X).  Set ac_cv_fc_cifs=yes if so.
AC_DEFUN([_ACX_SL_FC_CHECK_CIFS],
   [AC_CACHE_CHECK([whether the filesystem is case-insensitive],
       [ac_cv_fc_cifs],
       [rm -f conftest.*
        echo wibble >conftest.F
        AS_IF([test -f conftest.f && test "`cat conftest.f`" = wibble],
          [ac_cv_fc_cifs=yes], [ac_cv_fc_cifs=no])
])])# _ACX_SL_FC_CHECK_CIFS



# -----------------------
# User macros
# -----------------------

# AC_PROG_FPP([required features], [SUFFIX])
# --------------------------------------------------
#
# [required features] is a space-separated list of features that the Fortran
# preprocessor must have for the code to compile.
# It is up to the package maintainer to properly set these requirements.
#
# If SUFFIX is set it's assumed that files for the preprocessor have
# this suffix, otherwise .F is assumed. The output of the preprocessor
# is placed in files with the extension .f or whatever was last setup
# with AC_FC_SRCEXT.
#
# This macro will find out how to compile a preprocessable fixed-form
# file, with a .SUFFIX file extension. The type of Fortran code used
# for tests will be determined by the currently selected AC_LANG
# (Fortran or Fortran 77 at the time of this writing)
#
# See the section Preprocessor features for discussion of individual
# features.
#
# This macro sets and substitutes the variables FPP and FPPFLAGS.
#
# The macro depends on both FC and CPP, because we must possibly fall
# back on CPP for preprocessing.
#
AC_DEFUN([AC_PROG_FPP],
  [AC_REQUIRE([AC_PROG_FC])dnl
dnl We are not going to use AC_REQUIRE(AC_PROG_CPP) here for
dnl two reasons:
dnl 1. we don't really need to if FC will preprocess itself
dnl 2. we can't pass in an optional parameter to change the
dnl    default CPP search order, which we need to.
dnl AC_REQUIRE([AC_PROG_CPP([cpp])])dnl

   # Prefer AC_PROG_FC to AC_PROG_F77
   export FC FCFLAGS
   AS_IF([test "X$F77" != X],
     [AC_MSG_WARN([Use A@&t@C_PROG_FC with A@&t@C_PROG_FPP, instead of A@&t@C_PROG_F77])])


   AC_ARG_VAR([FPP], [Command to preprocess Fortran code])
   AC_ARG_VAR([FPPFLAGS], [Flags for the Fortran preprocessor])

   _ACX_SL_PROG_FPP([$2])

   # ACX_SL_PROG_FC_FPP_FEATURES does the actual feature tests,
   # storing results of the checks in non-cv variables like
   # ac_prog_fc_cpp_*, which we copy to cv variables afterwards.  This
   # allows this macro to be reusable for other cv variables (see
   # below)
   ACX_SL_PROG_FC_FPP_FEATURES([$1],[indirect],[$2],,
      [AC_MSG_FAILURE([required Fortran preprocessor not available])])
  ])# AC_PROG_FPP

#
# ACX_FC_INTEGRAL_FPP(FEATURES,[SUFFIX],[ACTION-IF-TRUE],[ACTION-IF-NOT])
# Determines wether $FC can process files containing preprocessor directives
# and processes each of the requested features correctly
# and run the according ACTION-*.
AC_DEFUN([ACX_FC_INTEGRAL_FPP],
  [# On nearly all systems where direct compilation is possible, a
   # .F file containing preprocessable fixed-form will be compiled
   # and preprocessed automatically. However, case-insensitive
   # filesystems (eg HFS+ on MacOSX) may get confused.  Therefore,
   # we must check for cpp flags.
   AC_ARG_VAR([FPPFLAGS], [Flags for the Fortran preprocessor])
   acx_sl_prog_fc_cpp=no
   AC_COMPILE_IFELSE([AC_LANG_PP_MANDATORY_SOURCE],
     [acx_sl_prog_fc_cpp=yes], [AS_VAR_SET([acx_sl_fpp_ok], [no])])
   # It is possible we've failed the previous test because of a
   # Tru64 bug where the compiler fails when called as 'f95' on
   # a .F file. It works when called as f90.
   #FIXME: this does not protect the user's setting of FC, though
   # we set it back if sensible.
   AS_IF([test x$acx_sl_prog_fc_cpp = xno && test $FC = f95],
     [FC=f90
      AC_LINK_IFELSE([AC_LANG_PP_MANDATORY_SOURCE],
        [acx_sl_prog_fc_cpp=yes],[FC=f95
         AS_VAR_SET([acx_sl_fpp_ok], [no])])
     ])
   ACX_SL_PROG_FC_FPP_FEATURES([$1],[indirect],[$2],,
      [AC_MSG_FAILURE([required Fortran preprocessor not available])])
   ac_first_save_FPPFLAGS=$FPPFLAGS
   FPPFLAGS="$FPPFLAGS $FPPFLAGS_F"
   # Don't need to test if $FC removes C++ comments - that
   # way madness lies.
   ACX_SL_PROG_FC_FPP_FEATURES([$1],[direct],[$2],[$3],[$4])
])

AC_DEFUN([_ACX_CHOOSE_FPP_BUILD_RULE],
 [# If so, we don't need to go any further.
if test x$acx_sl_fpp_ok = xyes; then
  ac_cv_fpp_build_rule=direct
  AC_MSG_RESULT([direct])
else
# indirect compilation
  AC_MSG_RESULT([indirect])

# Before we go any further, check that we're not courting disaster,
# here, by using indirect compilation (.F -> .f -> .o) on a
# case-insensitive filesystem.  If we are, there's nothing we can do
# other than fail noisily.
_ACX_SL_FC_CHECK_CIFS
if test $ac_cv_fc_cifs = yes; then
    AC_MSG_ERROR([disaster: this Fortran needs indirect compilation, but we
 have a case-insensitive filesystem, so .F -> .f would fail; further compilation isn't going to work -- consider filing a bug])
fi

# Now we check how to invoke a preprocessor that outputs Fortran code
# that FC can understand
#FIXME: in a joint C/Fortran project, CPP might have already
# been defined. Here we are potentially (probably) redefining it.
# I don't think this matters. Not sure, though.
# In that case, AC_SUBST has already been called on CPP.
# We don't want to fail if we can't find cpp - we might be able
# to fall back on fpp.
#FIXME: actually, we should just prefer cpp to $CPP
# The next macro sets FPP (unless already set by the user)
_ACX_SL_PROG_FPP

# Redefine the compile and link commands for indirect compilation
  ac_fpp_compile='${FPP-fpp} $FPPFLAGS $FPPFLAGS_SRCEXT conftest.$ac_ext '"$ac_fpp_out"' && ${FC-fc} -c $FCFLAGS conftest.f >&AS_MESSAGE_LOG_FD'
  ac_fpp_link='${FPP-fpp} $FPPFLAGS conftest.$ac_ext $FPPFLAGS_SRCEXT '"$ac_fpp_out"' && ${FC-fc} -o conftest${ac_exeext} $FCFLAGS $LDFLAGS conftest.f $LIBS >&AS_MESSAGE_LOG_FD'

  ac_compile=$ac_fpp_compile
  ac_link=$ac_fpp_link
# Redo all the feature checks for indirect compilation.

  ACX_SL_PROG_FC_FPP_FEATURES([$1], [indirect],
     ,[AC_MSG_FAILURE([required fortran preprocessor not available])])

if test $ac_fpp_need_d = yes; then
  AC_CACHE_CHECK([whether $FPP accepts -D],
     ac_cv_prog_fpp_d,
    [ac_cv_prog_fpp_d=$ac_prog_fc_cpp_d])
fi

if test $ac_fpp_need_i = yes; then
  AC_CACHE_CHECK([whether $FPP accepts -I],
     ac_cv_prog_fpp_i,
    [ac_cv_prog_fpp_i=$ac_prog_fc_cpp_i])
fi

if test $ac_fpp_need_subs = yes; then
  AC_CACHE_CHECK([whether $FPP substitutes macros in Fortran code],
     ac_cv_prog_fpp_subs,
    [ac_cv_prog_fpp_subs=$ac_prog_fc_cpp_subs])
fi

if test $ac_fpp_need_wrap = yes; then
  AC_CACHE_CHECK([whether $FPP wraps long lines automatically],
     ac_cv_prog_fpp_wrap,
    [ac_cv_prog_fpp_wrap=$ac_prog_fc_cpp_wrap])
fi

if test $ac_fpp_need_CSTYLE = yes; then
  AC_CACHE_CHECK([whether $FPP suppresses C-style comments],
     ac_cv_prog_fpp_CSTYLE,
    [ac_cv_prog_fpp_CSTYLE=$ac_prog_fc_cpp_CSTYLE])

elif test $ac_fpp_need_cstyle = yes; then
# It only makes sense to test this for indirect compilation,
# i.e., if .f files are generated
    _ACX_SL_PROG_FPP_CSTYLE
fi

if test $ac_fpp_need_cxxstyle = yes; then
  AC_CACHE_CHECK([whether $FPP preserves C++-style comments],
     ac_cv_prog_fpp_cxxstyle,
    [ac_cv_prog_fpp_cxxstyle=$ac_prog_fc_cpp_cxxstyle])
fi

AC_CACHE_CHECK([whether $FPP fulfils requested features],
  ac_cv_prog_fpp_ok,
  [ac_cv_prog_fpp_ok=AS_VAR_GET([acx_sl_fpp_ok])])

  ac_cv_fpp_build_rule=indirect

fi # test acx_sl_fpp_ok != yes
])
