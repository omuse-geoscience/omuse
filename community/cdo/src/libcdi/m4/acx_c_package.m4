dnl
dnl Copyright  (C)  2010  Thomas Jahns <jahns@dkrz.de>
dnl
dnl Version: 1.0
dnl Keywords: configure configure.ac autotools
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
dnl ACX_C_PACKAGE(PACKAGE,
dnl    [INCLUDE], [EXTRA-INCLUDES], [EXTRA-INCLUDEFLAGS],
dnl       [ACTION-IF_HEADER-NOT-FOUND],
dnl    [FUNCTION], [LIB-CANDIDATES], [EXTRA-LIBS], [EXTRA-LIBFLAGS],
dnl      [ACTION-IF-LIB-NOT-FOUND],
dnl    [DEFAULT-ROOT])
dnl -------------------------------------------------------------------
dnl Check wether INCLUDE can be compiled and FUNCTION is found in
dnl LIB-CANDIDATES. Sets PACKAGE_C_LIB and PACKAGE_C_INCLUDE variables to
dnl necessary C compiler switches and arguments such that
dnl inclusion/compilation succeed for program using INCLUDE or
dnl FUNCTION respectively.
dnl Also defines configure --with arguments for PACKAGEROOT,
dnl PACKAGE-LIB and PACKAGE-INCLUDE.
AC_DEFUN([ACX_C_PACKAGE],
  [AC_LANG_PUSH([C])
   ACX_GENERIC_PACKAGE([$1],[$2],[-I],[$3],[$4],[$5],[$6],[-L],[$7],[$8],[$9],[$10],[$11],[$12])
   AC_LANG_POP([C])dnl
  ])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl End:
