dnl acx_check_cfortran.m4 --- test if a program compiled from
dnl                           a main Fortran program and
dnl                           C functions gives expected results
dnl
dnl Copyright  (C)  2013  Thomas Jahns <jahns@dkrz.de>
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
dnl ACX_CHECK_CFORTRAN([OPTIONAL-CFORTRAN-INC-DIR],[ACTION-IF-SUCCESS],
dnl                    [ACTION-IF-FAILED])
dnl Test if C compiler can produce objects that link fine to Fortran programs
dnl when using cfortran.h.
dnl
AC_DEFUN([ACX_CHECK_CFORTRAN],
  [AC_CACHE_CHECK([if C externals constructed with cfortran.h work],
     [acx_cv_cfortran_works],
     [acx_cv_cfortran_works=no
      save_CPPFLAGS=$CPPFLAGS
      CPPFLAGS="-I]m4_ifval([$1],[$1],[$srcdir/include])[ $CPPFLAGS"
dnl build C function
      AC_LANG_PUSH([C])
      AC_COMPILE_IFELSE([AC_LANG_SOURCE([@%:@include "cfortran.h"
@%:@include <math.h>
@%:@include <stdio.h>
@%:@include <stdlib.h>

PROTOCCALLSFFUN1(FLOAT,CONFTEST_F,conftest_f,FLOAT)
#define conftest_F(v) \
  CCALLSFFUN1(CONFTEST_F,conftest_f,FLOAT,(v));

static float
conftest_C(int i, float v, int *p, float *q)
{
  float f;
  *p = (int)roundf(v * i);
  *q = f = conftest_F(v * i);
  return f;
}

FCALLSCFUN4(FLOAT,conftest_C,CONFTEST_C,conftest_c,INT,FLOAT,PINT,PFLOAT)

/* test string returns */
static const char *
conftest_str_C(void)
{
  static const char msg@<:@100@:>@ = "aaaaaaaaaaaaaaaaaaaa"
    "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
    "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
  return msg;
}

FCALLSCFUN0(STRING,conftest_str_C,CONFTEST_STR_C,conftest_str_c)

/* This function is required simply because some Fortran compilers
 * won't stop with exit code n when encountering STOP n */
static void
errExit(void)
{
  exit(1);
}

FCALLSCSUB0(errExit,ERR_EXIT,err_exit)

])],
        [_AC_RUN_LOG([mv "conftest.$ac_objext" "conftest_c.$ac_objext"],
           [_AS_ECHO_LOG([Renaming C object file.])])
         AC_LANG_PUSH([Fortran])
         AC_COMPILE_IFELSE([AC_LANG_SOURCE(
[      MODULE conftest_data
      IMPLICIT NONE
      PRIVATE
      REAL :: ri
      PUBLIC :: ri
      END MODULE conftest_data

      FUNCTION conftest_f(v) RESULT(r)
      USE conftest_data, ONLY: ri
      REAL, INTENT(in) :: v
      REAL :: r
      r = v * 100.0
      ri = 1.0 / v
      END FUNCTION conftest_f])],
           [_AC_RUN_LOG([mv "conftest.$ac_objext" "conftest_f.$ac_objext"],
              [_AS_ECHO_LOG([Renaming Fortran object file.])])
            save_LIBS=$LIBS
            LIBS="conftest_c.$ac_objext conftest_f.$ac_objext $LIBS"
            AC_RUN_IFELSE([AC_LANG_PROGRAM(,
[      USE conftest_data, ONLY: ri
      IMPLICIT NONE
      INTERFACE
       FUNCTION conftest_c(i, v, p, q) RESULT(f)
         INTEGER, INTENT(in) :: i
         REAL, INTENT(in) :: v
         INTEGER, INTENT(out) :: p
         REAL, INTENT(out) :: q
         REAL :: f
       END FUNCTION conftest_c
       FUNCTION conftest_str_c() result(s)
         CHARACTER(99) :: s
       END FUNCTION conftest_str_c
      END INTERFACE
      REAL, PARAMETER :: eps = 10e-6
      REAL :: foo, boo, too
      INTEGER :: bar, baz, i
      CHARACTER(99) :: aaaaaa
      bar = 5
      foo = 0.3
      too = conftest_c(bar, foo, baz, boo)
      IF (ABS(baz - NINT(bar * foo)) /= 0) THEN
        WRITE (0, '(2(a,i0))') "error checking, when baz, baz=", baz, &
             ", NINT(bar * foo) =", NINT(bar * foo)
        FLUSH(0)
        CALL err_exit
      END IF
      IF (ABS((ri - 1.0 / (bar * foo)) / ABS(ri)) > eps)  THEN
        WRITE (0, '(2(a,g24.15))') "error checking ri, ri=", ri, ", 1.0 / &
             &(bar * foo) = ", 1.0 / (bar * foo)
        FLUSH(0)
        CALL err_exit
      END IF
      IF (ABS((boo - (bar * foo * 100.0))/ABS(boo)) > eps)  THEN
        WRITE (0, '(2(a,g24.15))') "error checking boo, boo=", boo, &
             ", bar * foo * 100.0 = ", bar * foo * 100.0
        FLUSH(0)
        CALL err_exit
      END IF
      IF (too /= boo) THEN
        WRITE (0, '(2(a,g24.15))') "error checking too vs. boo, too=", too, &
             ", boo = ", boo
        FLUSH(0)
        CALL err_exit
      END IF
      aaaaaa = conftest_str_c()
      DO i = 1, 99
        IF (aaaaaa(i:i) /= 'a') THEN
          WRITE (0, '(a,i0,a)') "error checking aaaaaa(", i, ")=", &
              aaaaaa(i:i)
          FLUSH(0)
          CALL err_exit
        END IF
      END DO])],
              [acx_cv_cfortran_works=yes],
              [acx_cv_cfortran_works="error"],
              [AC_MSG_NOTICE([Skipping run test for cfortran.h in cross-compilation mode,])
	       AC_MSG_NOTICE([link test succeeded.])
               acx_cv_cfortran_works=yes])
	    LIBS=$save_LIBS],
           [acx_cv_cfortran_works="error compiling Fortran subroutine"])
         AC_LANG_POP([Fortran])],
        [acx_cv_cfortran_works="compiling with cfortran.h failed"])
      AC_LANG_POP([C])
      CPPFLAGS=$save_CPPFLAGS
     ])
   m4_ifval([$3],
     [AS_IF([test x"$acx_cv_cfortran_works" = xyes],[$2],[$3])],
     [AS_CASE([x"$acx_cv_cfortran_works"],
       [x"error"],
       [AC_MSG_FAILURE([Linking/Running with C EXTERNAL built with cfortran.h does not work!])],
       [x"compiling with cfortran.h failed"],
       [AC_MSG_FAILURE([Compilation with cfortran.h is not working!])],
       [x"error compiling Fortran subroutine"],
       [AC_MSG_FAILURE([compilation of simple Fortran source failed!])],
       [xyes],[$2],
       [AC_MSG_FAILURE([Unexpected error when linking C and Fortran via cfortran.h!])])])
  ])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://www.dkrz.de/redmine/projects/show/scales-ppm"
dnl license-default: "bsd"
dnl End:
