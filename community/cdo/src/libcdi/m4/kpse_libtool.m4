# Public macros for the TeX Live (TL) tree.
# Copyright (C) 1995 - 2009 Karl Berry <tex-live@tug.org>
# Copyright (C) 2009, 2010 Peter Breitenlohner <tex-live@tug.org>
#
# This file is free software; the copyright holders
# give unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# serial 1

# _KPSE_USE_LIBTOOL()
# Switch link tests over to use libtool so as not to require dependent
# libraries to be listed explicitly.
# Extended for Fortran by Thomas Jahns <jahns@dkrz.de>, 2015
# -------------------
AC_DEFUN([_KPSE_USE_LIBTOOL],
[## $0: Generate a libtool script for use in configure tests
AC_PROVIDE_IFELSE([LT_INIT], ,
                  [m4_fatal([$0: requires libtool])])[]dnl
LT_OUTPUT
m4_append([AC_LANG(C)],
[ac_link="./libtool --mode=link --tag=CC $ac_link"
])[]dnl
AC_PROVIDE_IFELSE([AC_PROG_CXX],
[m4_append([AC_LANG(C++)],
[ac_link="./libtool --mode=link --tag=CXX $ac_link"
])])[]dnl
AC_PROVIDE_IFELSE([AC_PROG_FC],
[m4_append([AC_LANG(Fortran)],
[ac_link="./libtool --mode=link --tag=FC $ac_link"
])])[]dnl
AC_PROVIDE_IFELSE([AC_PROG_F77],
[m4_append([AC_LANG(Fortran 77)],
[ac_link="./libtool --mode=link --tag=F77 $ac_link"
])])[]dnl
AC_LANG(_AC_LANG)[]dnl
]) # _KPSE_USE_LIBTOOL
