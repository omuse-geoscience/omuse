dnl acx_fc_check_strprt_convert.m4 --- unset a shell variable
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
dnl FIXME: this is currently only a placeholder that expects unset to
dnl be present although it's not universally available
dnl
dnl
dnl Code:
dnl
AC_DEFUN([ACX_FC_CHECK_STRPTR_CONVERT],
  [AC_LANG_PUSH([Fortran])
   AC_MSG_CHECKING([if Fortran compiler can handle complex CHARACTER interoperability])
   AC_COMPILE_IFELSE([AC_LANG_SOURCE([module conftest_mod
  use iso_c_binding
  implicit none
  private
  public :: errstr
contains
  function errstr(errno)
    integer, intent(in) :: errno
    interface
      function strerror(errno) bind(c, name='strerror')
        import :: c_int, c_ptr
        integer(c_int), value, intent(in) :: errno
        type(c_ptr) :: strerror
      end function strerror
      function strlen(s) bind(c, name='strlen')
        import :: c_ptr, c_size_t
        type(c_ptr), value, intent(in) :: s
        integer(c_size_t) :: strlen
      end function strlen
    end interface
    type(c_ptr) :: cptr
    character(len=:, kind=c_char), pointer :: errstr

    cptr = strerror(int(errno, c_int))
    errstr => c2f_string(cptr, int(strlen(cptr)))
  end function errstr

  function c2f_string(s, slen)
    type(c_ptr), intent(in) :: s
    integer, intent(in) :: slen
    CHARACTER(len=slen, kind=c_char), POINTER :: c2f_string
    c2f_string => NULL()
    call c_f_pointer(s, c2f_string)
  end function c2f_string
end module conftest_mod

program conftest
  use iso_c_binding
  use conftest_mod, only: errstr
  implicit none
  character(kind=c_char, len=:), pointer :: msg
  msg => errstr(42)
  write (0, '(a)') msg
end program conftest])],
     [AC_MSG_RESULT([yes])
m4_ifval([$1],[$1])],
     [AC_MSG_RESULT([no])
m4_ifval([$2],[$2])])
   AC_LANG_POP([Fortran])])

dnl
dnl Local Variables:
dnl mode: autoconf
dnl End:
dnl
