#  Process this file with autoconf to produce a configure script.
#
#  automake 1.14
#  autoconf 2.68
#  libtool  2.4.2

AC_INIT([cdi], [1.7.0rc3], [http://mpimet.mpg.de/cdi])

echo "configuring ${PACKAGE_NAME} ${PACKAGE_VERSION}"

CONFIG_ABORT=yes
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([m4])
dnl forbid acx macro names from remaining unexpanded
m4_pattern_forbid([^ACX_])
AC_CANONICAL_HOST
AC_CANONICAL_BUILD


AM_INIT_AUTOMAKE([foreign serial-tests])
AC_CONFIG_HEADERS([src/config.h])
AM_MAINTAINER_MODE([disable])

# Check building environment
AC_CHECK_TOOL([CC],[gcc],[:])
AC_PROG_CC_C99
AM_PROG_CC_C_O
AC_C_RESTRICT
AC_PROG_FC
AS_IF([test -n "$FC" && test "X$FC" != Xno],
  [AC_FC_SRCEXT([f90])
   AC_LANG_PUSH([Fortran])
   AC_PROG_FPP
   AC_LANG_POP([Fortran])])
AC_PROG_F77
AC_PROG_CXX
AC_PROG_INSTALL

# Set up libtool.
AC_MSG_NOTICE([setting up libtool])
dnl before switching on libtool, identify compilers that prevent us from
dnl certain build configurations
ACX_LT_PROBLEMS
LT_INIT([pic-only])
dnl _KPSE_USE_LIBTOOL ensures libtool is also used for configure-time tests,
dnl which deduces dependent libraries automatically
_KPSE_USE_LIBTOOL

# -----------------------------------------------------------------------
# Check endianess of system
AC_C_BIGENDIAN
#  ----------------------------------------------------------------------
# Check large file support on 32 bit systems
AC_SYS_LARGEFILE
#  ----------------------------------------------------------------------
# Checks for library functions.
AC_FUNC_MMAP
#  ----------------------------------------------------------------------
# Checks for library malloc.
AC_CHECK_LIB(malloc, malloc)
#  ----------------------------------------------------------------------
# Checks for structures.
AC_CHECK_MEMBERS([struct stat.st_blksize])
#  ----------------------------------------------------------------------
# Checks for header files
AC_CHECK_HEADERS([malloc.h uuid.h uuid/uuid.h])
AX_EXECINFO
#  ----------------------------------------------------------------------
# Checks for the availability of functions
AC_CHECK_FUNCS([mallinfo])
AC_CHECK_FUNCS([getline])
AC_CHECK_LIB([uuid],[uuid_generate],
  [AC_CHECK_DECLS([uuid_generate],
  [UUID_C_LIB=-luuid],,[AC_INCLUDES_DEFAULT
@%:@include <uuid/uuid.h>])],
  [AC_CHECK_FUNCS([uuid_create],
     [AC_CHECK_DECLS([uuid_create],
        [UUID_C_LIB=""],,[AC_INCLUDES_DEFAULT
@%:@include <uuid.h>])])])
AC_SUBST([UUID_C_LIB])
#  ----------------------------------------------------------------------
# Checks for the availability of ANSI-C99 functions
AC_CHECK_DECLS([isnan],,,[AC_INCLUDES_DEFAULT
@%:@include <math.h>])

#
# Check for non-standard builtin
AC_CHECK_DECLS([__builtin_ctz])

# Check compiler version
case "$CC" in
  pgcc*)  COMP_VERSION=`$CC -V | head -2 | tail -n 1`;;
  *gcc*)  COMP_VERSION=`$CC --version | head -n 1`;;
  g++*)   COMP_VERSION=`$CC --version | head -n 1`;;
  clang*) COMP_VERSION=`$CC --version | head -n 1`;;
  sxc*)   COMP_VERSION=`$CC -V 2>&1   | tail -n 1`;;
  xlc*)   COMP_VERSION=`$CC -qversion 2>&1   | head -n 1`;;
  *)      COMP_VERSION=`$CC -V 2>&1   | head -n 1`;;
esac

if test -z "$COMP_VERSION" ; then COMP_VERSION="unknown"; fi;
AC_DEFINE_UNQUOTED(COMP_VERSION, ["$COMP_VERSION"], [Compiler version])

# Checks for username, hostname and system type
USERNAME=$LOGNAME
if test -z "$USERNAME" ; then USERNAME=$USER; fi;
if test -z "$USERNAME" ; then USERNAME="unknown"; fi;
AC_DEFINE_UNQUOTED([USER_NAME],["$USERNAME"], [User name])
AC_SUBST([USER_NAME],["$USERNAME"])

AS_IF([test -z "$HOST"],
      [HOST=unknown
       AS_IF([test -x /bin/hostname],[HOST=$(hostname)],
             [AS_IF([test -x /bin/uname],
                    [HOST=$(uname -n)])])])
AC_DEFINE_UNQUOTED([HOST_NAME],["$HOST"],[Host name])
AC_SUBST([HOST_NAME],["$HOST"])

AC_DEFINE_UNQUOTED([SYSTEM_TYPE],["$ac_cv_build"], [System type])
AC_SUBST([SYSTEM_TYPE],["$ac_cv_build"])
#  ----------------------------------------------------------------------
#  Check for math library
AC_CHECK_LIB(m, floor)
#  ----------------------------------------------------------------------
# Add basic configure options
ACX_OPTIONS
AM_CONDITIONAL([ENABLE_NETCDF],[test x$ENABLE_NETCDF = xyes])
AC_CHECK_FUNCS([grib_get_length])
#  ----------------------------------------------------------------------
# Compile with MPI support
AC_ARG_ENABLE([mpi],
  AS_HELP_STRING([--enable-mpi],
    [Compile with MPI compiler @<:@default=no@:>@]),
    [enable_mpi=yes
     ACX_C_PACKAGE([MPI],[mpi.h],,,
       [AC_MSG_FAILURE([Required header mpi.h not found or not compilable.])
        enable_MPI=no],[MPI_Waitall],[mpi mpich],,,
       [AC_MSG_FAILURE([Cannot link C MPI programs.])
        enable_MPI=no])
     AS_IF([test -n "$FC" && test "X$FC" != "Xno"],
       [ACX_FORTRAN_PACKAGE([MPI], [mpif.h],,,
          [AC_MSG_FAILURE([Required include mpif.h not found or not compilable.])
           enable_MPI=no],[mpi_waitall],
          [mpi mpi_f90 mpi_f77 mpich],[[-lmpi_f77 -lmpi],[-lmpi]],,
          [AC_MSG_FAILURE([Cannot link Fortran MPI programs.])
           enable_MPI=no],[])])],
    [enable_mpi=no])
AS_IF([test x"${enable_mpi}" = x"yes"],
      [USE_MPI=yes])
HAVE_PARALLEL_NC4=0
enable_ppm=no
AS_IF([test x"$USE_MPI" = xyes],
  [AC_DEFINE([USE_MPI],[1],[parallel I/O requested and available])
   AC_PATH_PROGS([MPI_LAUNCH],[mpirun mpiexec],[true])
   AS_IF([test x"$cross_compiling" = xno],
     [AC_MSG_CHECKING([if $MPI_LAUNCH works])
      saved_CFLAGS=$CFLAGS
      saved_LIBS=$LIBS
      CFLAGS="$CFLAGS $MPI_C_INCLUDE"
      LIBS="$LIBS $MPI_C_LIB"
      AC_LINK_IFELSE([AC_LANG_SOURCE([
@%:@include <stdio.h>
@%:@include <stdlib.h>

@%:@include <mpi.h>

@%:@define xmpi(ret)           \\
  do {                      \\
    if (ret != MPI_SUCCESS) \\
      exit(EXIT_FAILURE);   \\
  } while (0)

int
main(int argc, char **argv)
{
  xmpi(MPI_Init(&argc, &argv));
  char *numarg = argv@<:@1@:>@;
  int cmdnum = atoi(numarg);
  int procnum = 1;
  xmpi(MPI_Allreduce(MPI_IN_PLACE, &procnum, 1, MPI_INT, MPI_SUM,
                     MPI_COMM_WORLD));
  xmpi(MPI_Finalize());
  return (procnum == cmdnum)?EXIT_SUCCESS:EXIT_FAILURE;
}
])],
         [AS_IF([$MPI_LAUNCH -n 4 ./conftest$EXEEXT 4],
            [AC_MSG_RESULT([yes])],
            [AC_MSG_RESULT([no])
             AC_MSG_FAILURE([mpirun doesn't work])])],
         [AC_MSG_FAILURE([Cannot compile simple MPI program])])
      CFLAGS=$saved_CFLAGS
      LIBS=$saved_LIBS
     ])
   AS_IF([test "x$MPI_LAUNCH" = xtrue],
     [AC_MSG_WARN([MPI launch command unavailable])])

   PKG_CHECK_MODULES([YAXT],[yaxt],
     [AC_DEFINE([HAVE_YAXT],,[yaxt library is available])],
     [AC_MSG_FAILURE([Required yaxt library unavailable.])])
   AS_IF([test x"$ENABLE_NC4" = xyes],
     [AC_CHECK_HEADERS([netcdf_par.h],
        [AC_CHECK_DECL([MPI_Bcast],
           [HAVE_PARALLEL_NC4=1],
           [AC_MSG_WARN([The supplied netCDF library does not support MPI parallel invocations])],
           [AC_INCLUDES_DEFAULT
@%:@include <netcdf_par.h>])],
        [AC_CHECK_DECL([MPI_Bcast],
           [HAVE_PARALLEL_NC4=1],
           [AC_MSG_WARN([The supplied netCDF library does not support MPI parallel invocations])],
           [AC_INCLUDES_DEFAULT
@%:@include <netcdf.h>])
])
      dnl parallel netCDF support still requires ScalES-PPM and YAXT to
      dnl re-arrange the data when running with more than one collector
      PKG_CHECK_MODULES([PPM_CORE],[scales-ppm-core],
        [enable_ppm=yes
         AC_DEFINE([HAVE_PPM_CORE],,
           [ScalES PPM C core library is available])
        ],
        [enable_ppm=no])
      dnl if not both scales-ppm and yaxt are available, netcdf can only be
      dnl used in serial mode
      AS_IF([test x$enable_ppm != xyes],
        [HAVE_PARALLEL_NC4=0])
     ])
  ],
  [MPI_LAUNCH="`pwd`/util/serialrun"
   AC_SUBST([MPI_LAUNCH])])
AS_IF([test $HAVE_PARALLEL_NC4 -gt 0],
  [AC_DEFINE([HAVE_PARALLEL_NC4],[1],
             [netCDF library does support MPI parallel invocations])])
AC_SUBST([HAVE_PARALLEL_NC4])
AM_CONDITIONAL([USE_MPI],[test x"$USE_MPI" = xyes])
AM_CONDITIONAL([USE_PPM_CORE],[test $enable_ppm = yes])
AC_SUBST([USE_MPI])
AC_SUBST([ENABLE_MPI],[`test x"$enable_mpi" = xyes && echo true || echo false`])
#  ----------------------------------------------------------------------
#  Create the Fortran Interface via iso_c_binding module (Fortran 2003 Standard)
AC_ARG_ENABLE([iso-c-interface],
              [AS_HELP_STRING([--enable-iso-c-interface],
                 [Create Fortran Interface via iso_c_bindings facility of F2003 [default=no].])],
              [enable_isoc=${enableval}],[enable_isoc=no])
AS_IF([test x"$enable_isoc" = xyes],
  [ACX_FC_CHECK_STRPTR_CONVERT],[enable_isoc=no])
AM_CONDITIONAL([CREATE_ISOC],[test $enable_isoc = 'yes'])
AC_SUBST([ENABLE_F2003_ISOC],[$enable_isoc])
AM_SUBST_NOTMAKE([ENABLE_F2003_ISOC])
# in case the Fortran interface uses the more modern F90-style interface,
# the FCFLAGS for users of the library have to include a switch to use the
# .mod file
CDI_F90_INTERFACE_FCFLAGS=''
AS_IF([test "x${enable_isoc}" = "xyes"],
  [AC_LANG_PUSH([Fortran])
   ACX_SL_FC_CHECK_MOD_PATH_FLAG
   AC_LANG_POP([Fortran])
   CDI_F90_INTERFACE_FCFLAGS="${FC_MOD_FLAG}${includedir}"])
AC_SUBST([CDI_F90_INTERFACE_FCFLAGS])
#  Check the module extension of the fortran compiler
AS_IF([test -n "$FC" && test "X$FC" != "Xno"],
  [ACX_SL_FC_MOD_SUFFIX(,[FCMODEXT=mod])])
# -----------------------------------------------------------------------
# Check for SWIG - Generator for script-language bindings
AC_ARG_ENABLE(swig,
              [AS_HELP_STRING([--enable-swig],[use swig to create extra bindings [default=no] (EXPERIMENTAL)])],
              [AC_CHECK_PROG(SWIG,swig,swig)],
              [enable_swig=no])
AM_CONDITIONAL(ENABLE_SWIG,[test "x$SWIG" != "x"])
#  ----------------------------------------------------------------------
#  Create the Ruby Interface via swig
AC_ARG_ENABLE(ruby,
  [AS_HELP_STRING([--enable-ruby],[ruby language bindings [default=no] (EXPERIMENTAL)])],
  [AC_CHECK_PROG(RUBY,ruby,ruby)
   RUBY_INCLUDES="$($RUBY $srcdir/config/interface.rb)"
   save_CPPFLAGS=$CPPFLAGS
   CPPFLAGS="$CPPFLAGS $RUBY_INCLUDES"
   AC_CHECK_HEADER([ruby.h],,[enable_ruby=no
                              RUBY=
                             ])
   CPPFLAGS=$save_CPPFLAGS
   AS_IF([test "x$RUBY" != "x"],
         [AS_IF([test "x$SWIG" = "x"],
                [AC_MSG_NOTICE([Ruby bindings: Swig is disabled, pregenerated files will be used. Use '--enable-swig' for swig usage])])],
                [AS_IF([test `$SWIG -help 2>&1 | $GREP -c '\-ruby *- Generate'` = 0],
                       [AC_MSG_NOTICE([No (swig does not support -ruby option)])
                        enable_ruby=no
                       ])])],
  [enable_ruby=no])
AM_CONDITIONAL(ENABLE_RUBY,[test "x$enable_ruby" != "xno"])
AS_IF([test "x$enable_ruby" != "xno"],[AC_SUBST([ENABLE_RUBY],[true])],[AC_SUBST([ENABLE_RUBY],[false])])
#  ----------------------------------------------------------------------
#  Create the Python Interface via swig
AC_ARG_ENABLE(python,
  [AS_HELP_STRING([--enable-python],[python language bindings [default=no] (EXPERIMENTAL)])],
  [AC_CHECK_PROG(PYTHON,python,python)
   enable_python=no
   for PYTHON_INCLUDES in \
     "-I`$PYTHON -c 'from distutils import sysconfig;print sysconfig.get_python_inc()'`" \
     "-I`$PYTHON -c 'from distutils import sysconfig;print sysconfig.get_python_lib(1,1)' | $SED 's/lib/include/'`"
   do
     save_CPPFLAGS=$CPPFLAGS
     CPPFLAGS="$CPPFLAGS $PYTHON_INCLUDES"
     AC_CHECK_HEADER([Python.h],[enable_python=yes
       break])
     CPPFLAGS=$save_CPPFLAGS
   done
   AS_IF([test "x$enable_python" != "xno"],
         [AS_IF([test "x$SWIG" = "x"],
                [AC_MSG_NOTICE([Python bindings: Swig is disabled, pregenerated files will be used. Use '--enable-swig' for swig usage])],
                [PYTHON=])],
                [AS_IF([test `$SWIG -help 2>&1 | $GREP -c '\-python *- Generate'` = 0],
                       [AC_MSG_NOTICE([No (swig does not support -python option)])
                        enable_python=no
                       ])])],
  [enable_python=no])
AM_CONDITIONAL(ENABLE_PYTHON, [test "x$enable_python" != "xno"])
AS_IF([test "x$enable_python" != "xno"],[AC_SUBST([ENABLE_PYTHON],[true])],[AC_SUBST([ENABLE_PYTHON],[false])])
#  ----------------------------------------------------------------------
AM_CONDITIONAL(CREATE_INTERFACES, [test  "x$enable_ruby" = "xyes" -o "x$enable_python" = "xyes"])
#  ----------------------------------------------------------------------
#  Create the CDI Fortran77 Interface via cfortran.h
AC_ARG_ENABLE([cf-interface],
              [AS_HELP_STRING([--disable-cf-interface],
                              [Omit building of cfortran Interface])],
              [enable_cfinterface=${enableval}],[enable_cfinterface=yes])
AS_IF([test "x${enable_cfinterface}" = "xyes"], [
dnl ######################################################################
dnl Test whether cfortran.h works correctly
dnl ######################################################################
ACX_FIND_CFORTRAN_DEF
ACX_XLF_QEXTNAME_ADD_APPENDUS
AS_IF([test -n "$FC" -a X"$FC" != Xno -o -n "$F77" -a X"$F77" != Xno],
  [ACX_CHECK_CFORTRAN([$srcdir/src],
    [AC_DEFINE([HAVE_CF_INTERFACE],[1],
       [Defined to 1 if C / Fortran interface cfortran.h works])],
    [AS_CASE([x"$acx_cv_cfortran_works"],
       [x"error"],
       [AC_MSG_NOTICE([Linking/Running with C EXTERNAL built with cfortran.h does not work!])],
       [x"compiling with cfortran.h failed"],
       [AC_MSG_NOTICE([Compilation with cfortran.h is not working!])],
       [x"error compiling Fortran subroutine"],
       [AC_MSG_NOTICE([compilation of simple Fortran source failed!])],
       [AC_MSG_NOTICE([Unexpected error when linking C and Fortran via cfortran.h!])])
     AC_MSG_NOTICE([Disabling cfortran.h bindings generation])
     acx_cv_cfortran_works=no])])
])
AM_CONDITIONAL([USE_FC],[test -n "$FC" -a X"$FC" != Xno -a x"$acx_cv_cfortran_works" = xyes])dnl
dnl
AC_SUBST([CPPFLAGS])dnl
dnl
m4_foreach([build_flag_var],[[BUILD_CFLAGS],[BUILD_FCFLAGS],[BUILD_LDFLAGS],[BUILD_LIBS],[BUILD_MPI_C_LIB],[BUILD_MPI_FC_LIB],[BUILD_C_INCLUDE],[BUILD_FC_INCLUDE]],
  [AC_ARG_VAR(build_flag_var,
     [append to ]m4_bpatsubst(build_flag_var, [BUILD_], [])[ during build but not in configure phase])dnl
AC_CONFIG_COMMANDS_PRE(m4_bpatsubst(build_flag_var, [BUILD_], [])[="$]m4_bpatsubst(build_flag_var, [BUILD_], [])[$]{build_flag_var+ $build_flag_var[}"])dnl
AM_SUBST_NOTMAKE(build_flag_var)])dnl

m4_foreach([build_tool_var],[[BUILD_CC],[BUILD_CXX],[BUILD_FC],[BUILD_F77]],
  [AC_ARG_VAR(build_tool_var,
     [replace ]m4_bpatsubst(build_tool_var, [BUILD_], [])[ with expansion of $]build_tool_var[ during build but not in configure phase])dnl
AC_CONFIG_COMMANDS_PRE(m4_bpatsubst(build_tool_var, [BUILD_], [])[="$]{build_tool_var-$[]m4_bpatsubst(build_tool_var, [BUILD_], [])}["])dnl
AM_SUBST_NOTMAKE(build_tool_var)])dnl
dnl
dnl
dnl
# Checks for compiler
COMPILER="$CC $CFLAGS"
AC_DEFINE_UNQUOTED(COMPILER, ["$COMPILER"], [Compiler])

AC_CONFIG_FILES([tests/test_cksum_grib \
                 tests/test_cksum_nc \
                 tests/test_cksum_nc2 \
                 tests/test_cksum_nc4 \
                 tests/test_cksum_extra \
                 tests/test_cksum_service \
                 tests/test_cksum_ieg \
                 tests/test_chunk_cksum \
                 tests/test_f2003 \
                 tests/pio_write_run \
                 tests/pio_write_deco2d_run \
                 tests/pio_cksum_mpinonb \
                 tests/pio_cksum_fpguard \
                 tests/pio_cksum_asynch \
                 tests/pio_cksum_writer \
                 tests/pio_cksum_cdf \
                 tests/test_resource_copy_mpi_run \
                 tests/test_cdf_transformation \
                 tables/gen_tableheaderfile \
                 util/serialrun],[chmod a+x "$ac_file"])

AC_OUTPUT([Makefile src/Makefile interfaces/Makefile app/Makefile \
	tests/Makefile examples/Makefile cdi.settings \
	examples/pio/Makefile src/pkgconfig/cdi.pc src/pkgconfig/cdipio.pc])

#  ----------------------------------------------------------------------
#  Show configuration
AC_MSG_NOTICE([CDI is configured with the following options:])
./config.status cdi.settings
cat cdi.settings

AS_IF([test "$ac_cv_prog_cc_c99" = "no"],
      [AC_MSG_NOTICE([

  Warning: The C compiler does not accept ANSI C99 source code!
])])
AC_MSG_NOTICE([

  Configuration completed.

    You can now say 'make' to compile the CDI package and 'make install' to install it afterwards.
])
