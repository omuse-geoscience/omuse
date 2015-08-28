dnl acx_lt_problems.m4 --- prevent problematic libtool build configurations
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
dnl ACX_LT_PROBLEMS
dnl Test if compiler is able to produce working shared objects
dnl and prevent builds of shared libraries that would not work anyway.
dnl
dnl The following problems are known so far:
dnl
dnl - Intel ifort versions 15.0.1, 15.0.2 or 15.0.3 compile the Fortran
dnl     runtime into any shared object, thus producing various
dnl     interpositioning problems with little likelyness that
dnl     executables linking the shared objects can actually work.
dnl
dnl - IBM mp-compilers (mpcc, mpcc_r, mpCC, mpCC_r, mpxlf, mpxlf_r,
dnl     mpxlf2003_r, mpxlf90, mpxlf90_r, mpxlf95, mpxlf95_r) link
dnl     shared objects with linker option -binitfini:poe_remote_main
dnl     even though this option should only be set for executables, it
dnl     causes a duplicate exit handler that fails and clobbers the
dnl     programs exit status (set to 128).
dnl
AC_DEFUN([_ACX_LT_PROBLEMS],
  [AC_REQUIRE([AC_CANONICAL_HOST])
   AC_LANG_CASE([Fortran],[acx_Comp=$FC],
     [Fortran 77],[acx_Comp=$F77],
     [C],[acx_Comp=$CC])
   _AS_ECHO_LOG([testing if $acx_Comp cannot build working shared objects])
   AS_CASE([$host],
     [*-ibm-aix*],
     [AS_IF([$acx_Comp -G -v 2>&1 | grep ' -binitfini:poe_remote_main ' >/dev/null],
        [acx_cv_disable_shared=yes])],
     [x86_64-*-linux-*|i*86-*-linux-*|*-apple-darwin*|ia64-*-linux-*|x86_64-*-freebsd*|i*86-*-freebsd*],
     [AS_IF([$acx_Comp -V 2>&1 | grep '^Intel(R).*Fortran.*Compiler.*Version 15.0.@<:@123@:>@' >/dev/null],
        [acx_cv_disable_shared=yes])])
   _AS_ECHO_LOG([result: $acx_cv_disable_shared])])
dnl
dnl run test for C and Fortran compilers
AC_DEFUN([ACX_LT_PROBLEMS],
  [AS_IF([test x"$enable_shared" != xno],
     [AC_CACHE_CHECK([any compiler has problems building shared objects],
        [acx_cv_disable_shared],
        [acx_cv_disable_shared=no
         AC_PROVIDE_IFELSE([AC_PROG_FC],
           [AC_LANG_PUSH([Fortran])
            _ACX_LT_PROBLEMS
            AC_LANG_POP([Fortran])])
         AC_PROVIDE_IFELSE([AC_PROG_F77],
           [AC_LANG_PUSH([Fortran 77])
            _ACX_LT_PROBLEMS
            AC_LANG_POP([Fortran 77])])
         AC_PROVIDE_IFELSE([AC_PROG_CC],
           [AC_LANG_PUSH([C])
            _ACX_LT_PROBLEMS
            AC_LANG_POP([C])])])
       AS_IF([test x"$acx_cv_disable_shared" = xyes],
         [enable_shared=no])])])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://www.dkrz.de/redmine/projects/show/scales-ppm"
dnl license-default: "bsd"
dnl End:
