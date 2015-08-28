AC_DEFUN([ACX_OPTIONS],
[
#  ----------------------------------------------------------------------
#  Checks for multithreaded compiling + linking
AC_ARG_WITH([threads],
            [AC_HELP_STRING([--with-threads=<yes/no/directory>],
                            [Compile + link for multithreading [default=yes]])],
            [],
            [with_threads=yes])
THREADS_INCLUDE=''
THREADS_LIBS=''
AS_CASE([$with_threads],
        [no],[AC_MSG_CHECKING([multithreading])
              AC_MSG_RESULT([suppressed])],
        [yes],[AX_PTHREAD([AC_DEFINE([HAVE_LIBPTHREAD],[1],[Define 1 for multithread support])],[AC_MSG_ERROR([multithreaded settings NOT found])])
               LIBS="$PTHREAD_LIBS $LIBS"
               CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
               CC="$PTHREAD_CC"
               AS_ECHO(["CC:$CC CFLAGS:$CFLAGS LIBS:$LIBS"])],
        [*],[THREADS_ROOT=$with_threads
             LDFLAGS="-L$THREADS_ROOT/lib $LDFLAGS"
             CPPFLAGS="-I$THREADS_ROOT/include $CPPFLAGS "
             AC_CHECK_HEADERS(pthread.h)
             AC_CHECK_LIB([pthread],[pthread_create])
             THREADS_LIBS=" -L$THREADS_ROOT/lib -lpthread"
             THREADS_INCLUDE=" -I$THREADS_ROOT/include"])
AC_SUBST([THREADS_INCLUDE])
AC_SUBST([THREADS_LIBS])
#  ----------------------------------------------------------------------
#  Link application to ZLIB library, needed for netcdf
ZLIB_INCLUDE=''
ZLIB_LIBS=''
AC_ARG_WITH([zlib],
            [AS_HELP_STRING([--with-zlib=<yes|no|directory> (default=yes)],[location of ZLIB compression library (lib and include subdirs), nec. for HDF5/NETCDF4])],
            [AS_CASE(["$with_zlib"],
                     [no],[AC_MSG_CHECKING([for ZLIB library])
                           AC_MSG_RESULT([suppressed])],
                     [yes],[AC_CHECK_HEADERS(zlib.h)
                            AC_SEARCH_LIBS([deflate],[z],[AC_DEFINE([HAVE_LIBZ],[1],[Define 1 for ZLIB support])])
                            ZLIB_LIBS=" -lz"],
                     [*],[ZLIB_ROOT=$with_zlib
                          LDFLAGS="-L$ZLIB_ROOT/lib $LDFLAGS"
                          CPPFLAGS="-I$ZLIB_ROOT/include $CPPFLAGS"
                          AC_CHECK_HEADERS(zlib.h)
                          AC_SEARCH_LIBS([deflate],[z],[AC_DEFINE([HAVE_LIBZ],[1],[Define 1 for ZLIB support])])
                          ZLIB_INCLUDE=" -I$ZLIB_ROOT/include"
                          ZLIB_LIBS=" -L$ZLIB_ROOT/lib -lz"])],
                     [AC_CHECK_HEADERS(zlib.h)
                      AC_SEARCH_LIBS([deflate],[z],[AC_DEFINE([HAVE_LIBZ],[1],[Define 1 for ZLIB support])])
              ZLIB_LIBS=" -lz"])
AC_SUBST([ZLIB_INCLUDE])
AC_SUBST([ZLIB_LIBS])
#  ----------------------------------------------------------------------
#  Compile application with SZLIB library, needed for GRIB1 or for
#  linking against hdf5/netcdf4
SZLIB_INCLUDE=''
SZLIB_LIBS=''
AC_ARG_WITH([szlib],
            [AS_HELP_STRING([--with-szlib=<yes|no|directory> (default=no)],[location of szlib library, optional for GRIB1 and NETCDF4 compression])],
            [AS_CASE(["$with_szlib"],
                     [no],[AC_MSG_CHECKING([for szlib library])
                           AC_MSG_RESULT([suppressed])],
                     [yes],[AC_CHECK_HEADERS(szlib.h)
                            AC_SEARCH_LIBS([SZ_BufftoBuffCompress],
                                           [sz],
                                           [AC_DEFINE([HAVE_LIBSZ],[1],[Define to 1 for SZIP support])],
                                           [AC_MSG_ERROR([Could not link to szlib])])
                            SZLIB_LIBS=" -lsz"],
                     [*],[SZLIB_ROOT=$with_szlib
                          AS_IF([test -d "$SZLIB_ROOT"],
                                [LDFLAGS="-L$SZLIB_ROOT/lib $LDFLAGS"
                                 CPPFLAGS="-I$SZLIB_ROOT/include $CPPFLAGS"
                                 AC_CHECK_HEADERS(szlib.h)
                                 AC_SEARCH_LIBS([SZ_BufftoBuffCompress],
                                                [sz],
                                                [AC_DEFINE([HAVE_LIBSZ],[1],[Define to 1 for SZIP support])],
                                                [AC_MSG_ERROR([Could not link to szlib])])
                                 SZLIB_LIBS=" -L$SZLIB_ROOT/lib -lsz"
                                 SZLIB_INCLUDE=" -I$SZLIB_ROOT/include"],
                                [AC_MSG_NOTICE([$SZLIB_ROOT is not a directory! SZLIB suppressed])])])],
            [AC_MSG_CHECKING([for szlib library])
             AC_MSG_RESULT([suppressed])])
AC_SUBST([SZLIB_INCLUDE])
AC_SUBST([SZLIB_LIBS])
#  ----------------------------------------------------------------------
#  Link application with HDF5 library, required for netcdf4
HDF5_ROOT=''
HDF5_INCLUDE=''
HDF5_LIBS=''
AC_ARG_WITH([hdf5],
            [AS_HELP_STRING([--with-hdf5=<yes|no|directory> (default=no)],[location of hdf5 library, NETCDF4 requires hdf5 high level interface])],
            [AS_CASE(["$with_hdf5"],
                     [no],[AC_MSG_CHECKING([for hdf5 library])
                           AC_MSG_RESULT([suppressed])],
                     [yes],[AC_CHECK_HEADERS([hdf5.h])
                            AC_SEARCH_LIBS([H5Fopen],
                                           [hdf5],
                                           [AC_DEFINE([HAVE_LIBHDF5],[1],[Define to 1 for HDF5 support])],
                                           [AC_MSG_ERROR([Cannot link to hdf5 library! It is required for Netcdf4])])
                            AC_SEARCH_LIBS([H5DSis_scale],
                                           [hdf5_hl],
                                           [have_hdf5_hl=yes],
                                           [AC_MSG_NOTICE([Cannot find hdf5 high level interface! It is required for netCDF4.])
                                            have_hdf5_hl=no])
                            AS_IF([test "x$have_libhdf5_hl" = xyes],
                                  [HDF5_LIBS=" -lhdf5_hl -lhdf5"],
                                  [HDF5_LIBS=" -lhdf5"])
                            ],
                     [*],[AS_IF([test -d "$with_hdf5"],
                                [HDF5_ROOT="$with_hdf5"
                                 LDFLAGS="-L$HDF5_ROOT/lib $LDFLAGS"
                                 CPPFLAGS="-I$HDF5_ROOT/include $CPPFLAGS"
                                 AC_CHECK_HEADERS([hdf5.h])
                                 AC_SEARCH_LIBS([H5Fopen],
                                                [hdf5],
                                                [AC_DEFINE([HAVE_LIBHDF5],[1],[Define to 1 for HDF5 support])],
                                                [AC_MSG_ERROR([Cannot link to hdf5! It is required for netCDF4.])])
                                 AC_SEARCH_LIBS([H5DSis_scale],
                                                [hdf5_hl],
                                                [have_hdf5_hl=yes],
                                                [AC_MSG_NOTICE([Cannot link to hdf5 high level interface! It is required for netCDF4.\
                                                                HDF5 must be built with zlib; the location of zlib must be specified for HDF5 with the \
                                                                --with-zlib option. If HDF5 was also built with szlib, then the location of szlib must also be \
                                                                specified with the --with-szlib option..])
                                                have_hdf5_hl=no])
                                 AS_IF([test "x$have_libhdf5_hl" = 'xyes'],
                                       [HDF5_LIBS=" -L$HDF5_ROOT/lib -lhdf5_hl -lhdf5"],
                                       [HDF5_LIBS=" -L$HDF5_ROOT/lib -lhdf5"])
                                 HDF5_INCLUDE=" -I$HDF5_ROOT/include"],
                                [AC_MSG_NOTICE([$with_hdf5 is not a directory! HDF5 suppressed])])])],
            [AC_MSG_CHECKING([for hdf5 library])
             AC_MSG_RESULT([suppressed])])
AC_SUBST([HDF5_ROOT])
AC_SUBST([HDF5_INCLUDE])
AC_SUBST([HDF5_LIBS])
#  ----------------------------------------------------------------------
#  Compile application with netcdf
NETCDF_ROOT=''
NETCDF_INCLUDE=''
NETCDF_LIBS=''
ENABLE_NETCDF=no
ENABLE_NC2=no
ENABLE_NC4=no
AC_ARG_WITH([netcdf],
            [AS_HELP_STRING([--with-netcdf=<yes|no|directory> (default=no)],[location of netcdf library (lib and include subdirs)])],
            [AS_CASE(["$with_netcdf"],
                     [no],[AC_MSG_CHECKING([for netcdf library])
                           AC_MSG_RESULT([suppressed])],
                     [yes],[AC_CHECK_HEADERS([netcdf.h])
                            AC_SEARCH_LIBS([nc_open],
                                           [netcdf],
                                           [AC_DEFINE([HAVE_LIBNETCDF],[1],[Define to 1 for NETCDF support])
                                            ENABLE_NETCDF=yes],
                                           [AC_MSG_ERROR([Could not link to netcdf library])])
                            NETCDF_LIBS=" -lnetcdf"
                            AC_CHECK_PROG(NC_CONFIG,nc-config,nc-config)
                            AS_IF([test "x$NC_CONFIG" != "x"],
                                  [AC_MSG_CHECKING([netcdf's nc2 support])
                                   AS_IF([test "x$($NC_CONFIG --has-nc2)" = "xyes"],
                                         [AC_DEFINE([HAVE_NETCDF2],[1],[Define to 1 for NETCDF2 support])
                                          ENABLE_NC2=yes
                                          AC_MSG_RESULT([yes])],[AC_MSG_RESULT([no])])
                                   AC_MSG_CHECKING([netcdf's nc4 support])
                                   AS_IF([test "x$($NC_CONFIG --has-nc4)" = "xyes"],
                                         [AC_DEFINE([HAVE_NETCDF4],[1],[Define to 1 for NETCDF4 support])
                                          ENABLE_NC4=yes
                                          AC_MSG_RESULT([yes])],[AC_MSG_RESULT([no])])],
                                  [AS_ECHO([Could not find nc-config! go on with default configuration])])],
                     [*],[AS_IF([test -d "$with_netcdf"],
                                [NETCDF_ROOT=$with_netcdf
                                 LDFLAGS="-L$NETCDF_ROOT/lib $LDFLAGS"
                                 CPPFLAGS="-I$NETCDF_ROOT/include $CPPFLAGS"
                                 AC_CHECK_HEADERS([netcdf.h])
                                 AC_SEARCH_LIBS([nc_open],
                                                [netcdf],
                                                [AC_DEFINE([HAVE_LIBNETCDF],[1],[Define to 1 for NETCDF support])
                                                 ENABLE_NETCDF=yes],
                                                [AC_MSG_ERROR([Could not link to netcdf library])])
                                 NETCDF_LIBS=" -L$NETCDF_ROOT/lib -lnetcdf"
                                 NETCDF_INCLUDE=" -I$NETCDF_ROOT/include"
                                 AC_MSG_CHECKING([nc-config script])
                                 AC_CHECK_PROG(NC_CONFIG,nc-config,[$NETCDF_ROOT/bin/nc-config],,["$NETCDF_ROOT/bin"])
                                 AS_IF([test "x$NC_CONFIG" != "x"],
                                   [AC_MSG_CHECKING([netcdf's OpenDAP support])
                                    AS_IF([test "x$($NC_CONFIG --has-dap)" = "xyes"],
                                          [AC_DEFINE([HAVE_LIBNC_DAP],[1],[Define to 1 for NETCDF OpenDAP])
                                           AC_MSG_RESULT([yes])],[AC_MSG_RESULT([no])])]
                                   [AC_MSG_CHECKING([netcdf's nc2 support])
                                   AS_IF([test "x$($NC_CONFIG --has-nc2)" = "xyes"],
                                         [AC_DEFINE([HAVE_NETCDF2],[1],[Define to 1 for NETCDF2 support])
                                          ENABLE_NC2=yes
                                          AC_MSG_RESULT([yes])],[AC_MSG_RESULT([no])])
                                   AC_MSG_CHECKING([netcdf's nc4 support])
                                   AS_IF([test "x$($NC_CONFIG --has-nc4)" = "xyes"],
                                         [AC_DEFINE([HAVE_NETCDF4],[1],[Define to 1 for NETCDF4 support])
                                          ENABLE_NC4=yes
                                          AC_MSG_RESULT([yes])],[AC_MSG_RESULT([no])])],
                                   [AC_MSG_RESULT([Could not find nc-config! go on with default configuration])])],
                                [AC_MSG_NOTICE([$with_netcdf is not a directory! NETCDF suppressed])])])],
            [AC_MSG_CHECKING([for NETCDF library])
             AC_MSG_RESULT([suppressed])])
AC_SUBST([ENABLE_NETCDF])
AC_SUBST([ENABLE_NC2])
AC_SUBST([ENABLE_NC4])
AC_SUBST([NETCDF_ROOT])
AC_SUBST([NETCDF_INCLUDE])
AC_SUBST([NETCDF_LIBS])
#  ----------------------------------------------------------------------
#  Link application with JASPER library (needed for GRIB2 compression)
JASPER_LIBS=''
AC_ARG_WITH([jasper],
            [AS_HELP_STRING([--with-jasper=<directory>],
                            [Specify location of JASPER library. You must specify its location if GRIB_API was built with JASPER.])],
            [AS_CASE(["$with_jasper"],
                     [no],[AC_MSG_CHECKING([for jasper library])
                           AC_MSG_RESULT([suppressed])],
                     [yes],[AC_CHECK_HEADERS([jasper.h])
                            AC_SEARCH_LIBS([jas_init],[jasper],[AC_DEFINE([HAVE_LIBJASPER],[1],[Define to 1 for JPEG compression for GRIB2])],
                                           [AC_MSG_ERROR([Could not link to jasper library! Required for GRIB_API])])
                            AC_SUBST([JASPER_LIBS],[" -ljasper"])],
                     [*],[JASPER_ROOT=$with_jasper
                          AS_IF([test -d "$JASPER_ROOT"],
                                [LDFLAGS="$LDFLAGS -L$JASPER_ROOT/lib"
                                 CPPFLAGS="$CPPFLAGS -I$JASPER_ROOT/include"
                                 AC_SEARCH_LIBS([jas_stream_memopen],
                                                [jasper],
                                                [AC_DEFINE([HAVE_LIBJASPER],[1],[Define to 1 for JPEG compression for GRIB2])],
                                                [AC_MSG_ERROR([Could not link to jasper library! Required for GRIB_API])])
                                 JASPER_LIBS=" -L$JASPER_ROOT/lib -ljasper"],
                                [AC_MSG_ERROR([$JASPER_ROOT is not a directory! JASPER suppressed])])])],
            [AC_MSG_CHECKING([for the JASPER library])
             AC_MSG_RESULT([suppressed])])
AC_SUBST([JASPER_LIBS])
#  ----------------------------------------------------------------------
#  Link application with openjpeg library (needed for GRIB2 compression)
OPENJPEG_LIBS=''
AC_ARG_WITH([openjpeg],
            [AS_HELP_STRING([--with-openjpeg=<directory>],
                            [Specify location of openjpeg library. You must specify its location if GRIB_API was built with openjpeg.])],
            [AS_CASE(["$with_openjpeg"],
                     [no],[AC_MSG_CHECKING([for openjpeg library])
                           AC_MSG_RESULT([suppressed])],
                     [yes],[AC_CHECK_HEADERS([openjpeg.h])
                            AC_SEARCH_LIBS([opj_image_create],[openjpeg],[AC_DEFINE([HAVE_LIBOPENJPEG],[1],[Define to 1 for JPEG compression for GRIB2])],
                                           [AC_MSG_ERROR([Could not link to openjpeg library! Required for GRIB_API])])
                            AC_SUBST([openjpeg_LIBS],[" -lopenjpeg"])],
                     [*],[OPENJPEG_ROOT=$with_openjpeg
                          AS_IF([test -d "$OPENJPEG_ROOT"],
                                [LDFLAGS="$LDFLAGS -L$OPENJPEG_ROOT/lib"
                                 CPPFLAGS="$CPPFLAGS -I$OPENJPEG_ROOT/include"
                                 AC_SEARCH_LIBS([opj_image_create],
                                                [openjpeg],
                                                [AC_DEFINE([HAVE_LIBOPENJPEG],[1],[Define to 1 for JPEG compression for GRIB2])],
                                                [AC_MSG_ERROR([Could not link to openjpeg library! Required for GRIB_API])])
                                 OPENJPEG_LIBS=" -L$OPENJPEG_ROOT/lib -lopenjpeg"],
                                [AC_MSG_ERROR([$OPENJPEG_ROOT is not a directory! openjpeg suppressed])])])],
            [AC_MSG_CHECKING([for the openjpeg library])
             AC_MSG_RESULT([suppressed])])
AC_SUBST([OPENJPEG_LIBS])
#  ----------------------------------------------------------------------
#  Link application with LIBPNG library (needed for GRIB2 compression)
LIBPNG_LIBS=''
AC_ARG_WITH([libpng],
            [AS_HELP_STRING([--with-libpng=<directory>],
                            [Specify location of LIBPNG library. You must specify its location if GRIB_API was built with LIBPNG.])],
            [AS_CASE(["$with_libpng"],
                     [no],[AC_MSG_CHECKING([for libpng library])
                           AC_MSG_RESULT([suppressed])],
                     [yes],[AC_CHECK_HEADERS([png.h])
                            AC_SEARCH_LIBS([png_warning],[png],[AC_DEFINE([HAVE_LIBLIBPNG],[1],[Define to 1 for PNG compression for GRIB2])],
                                           [AC_MSG_ERROR([Could not link to libpng library! Required for GRIB_API])])
                            AC_SUBST([LIBPNG_LIBS],[" -lpng"])],
                     [*],[LIBPNG_ROOT=$with_libpng
                          AS_IF([test -d "$LIBPNG_ROOT"],
                                [LDFLAGS="$LDFLAGS -L$LIBPNG_ROOT/lib"
                                 CPPFLAGS="$CPPFLAGS -I$LIBPNG_ROOT/include"
                                 AC_SEARCH_LIBS([png_warning],
                                                [png],
                                                [AC_DEFINE([HAVE_LIBLIBPNG],[1],[Define to 1 for PNG compression for GRIB2])],
                                                [AC_MSG_ERROR([Could not link to libpng library! Required for GRIB_API])])
                                 LIBPNG_LIBS=" -L$LIBPNG_ROOT/lib -lpng"],
                                [AC_MSG_ERROR([$LIBPNG_ROOT is not a directory! LIBPNG suppressed])])])],
            [AC_MSG_CHECKING([for the LIBPNG library])
             AC_MSG_RESULT([suppressed])])
AC_SUBST([LIBPNG_LIBS])
#  ----------------------------------------------------------------------
#  Compile application with GRIB_API library (for GRIB2 support)
GRIB_API_INCLUDE=''
GRIB_API_LIBS=''
AC_ARG_WITH([grib_api],
            [AS_HELP_STRING([--with-grib_api=<yes|no|directory>],
                            [library for grib2 compression; if a directory is given, it will be used as a value for --with-jasper-root])],
            [AS_CASE(["$with_grib_api"],
                     [no],[AC_MSG_CHECKING([for GRIB_API library])
                           AC_MSG_RESULT([suppressed])],
                     [yes],[AC_CHECK_HEADERS([grib_api.h])
                            AC_SEARCH_LIBS([grib_get_message],
                                           [grib_api],
                                           [AC_DEFINE([HAVE_LIBGRIB_API],[1],[GRIB_API library is present if defined to 1])],
                                           [AC_MSG_ERROR([Could not link to grib_api library])])],
                     [*],[GRIB_API_ROOT=$with_grib_api
                          AS_IF([test -d "$GRIB_API_ROOT"],
                                [LDFLAGS="-L$GRIB_API_ROOT/lib $LDFLAGS"
                                 CPPFLAGS="-I$GRIB_API_ROOT/include $CPPFLAGS"
                                 AC_CHECK_HEADERS([grib_api.h])
                                 AC_SEARCH_LIBS([grib_get_message],
                                                [grib_api],
                                                [AC_DEFINE([HAVE_LIBGRIB_API],[1],[GRIB_API library is present if defined to 1])],
                                                [AC_MSG_ERROR([Could not link to grib_api library])])
                                 GRIB_API_LIBS=" -L$GRIB_API_ROOT/lib -lgrib_api"
                                 GRIB_API_INCLUDE=" -I$GRIB_API_ROOT/include"],
                                [AC_MSG_ERROR([$GRIB_API_ROOT is not a directory! GRIB_API suppressed])])])],
            [AC_MSG_CHECKING([for the GRIB_API library])
             AC_MSG_RESULT([suppressed])])
AC_SUBST([GRIB_API_INCLUDE])
AC_SUBST([GRIB_API_LIBS])
AM_CONDITIONAL([HAVE_LIBGRIB_API],[test "x$with_grib_api" != 'x' -a "x$with_grib_api" != 'xno' ])
#  ----------------------------------------------------------------------
#  Enable GRIB support
AC_MSG_CHECKING([for GRIB support])
AC_ARG_ENABLE([grib],
              [AS_HELP_STRING([--enable-grib],[GRIB support [default=yes]])],
              [AS_IF([test "x$enable_grib" != 'xno'],
                     [AC_DEFINE(HAVE_LIBGRIB, [1], [Define to 1 for GRIB support])
                      enable_grib=yes])],
              [AC_DEFINE(HAVE_LIBGRIB, [1], [Define to 1 for GRIB support])
               enable_grib=yes])
AC_MSG_RESULT([$enable_grib])
AC_SUBST([ENABLE_GRIB],[$enable_grib])
#  ----------------------------------------------------------------------
#  Compile interface with internal CGRIBEX library
AC_MSG_CHECKING([for CGRIBEX support])
AC_ARG_ENABLE([cgribex],
              [AC_HELP_STRING([--enable-cgribex],[Use the CGRIBEX library [default=yes]])],
              [AS_IF([test "x$enable_cgribex" != 'xno'],
                     [AC_DEFINE(HAVE_LIBCGRIBEX,[1],[Define to 1 for GRIB1 decoding/encoding with cgribex])
                      enable_cgribex=yes])],
              [AC_DEFINE(HAVE_LIBCGRIBEX,[1],[Define to 1 for GRIB1 decoding/encoding with cgribex])
               enable_cgribex=yes])
AC_MSG_RESULT([$enable_cgribex])
AC_SUBST([ENABLE_CGRIBEX],[$enable_cgribex])
#  ----------------------------------------------------------------------
#  Compile interface with internal SERVICE library
AC_MSG_CHECKING([for SERVICE support])
AC_ARG_ENABLE([service],
              [AC_HELP_STRING([--enable-service],[Use the service library [default=yes]])],
              [AS_IF([test "x$enable_service" != 'xno'],
                     [AC_DEFINE(HAVE_LIBSERVICE,[1],[Define to 1 for SERVICE interface])
                      enable_service=yes])],
              [AC_DEFINE(HAVE_LIBSERVICE,[1],[Define to 1 for SERVICE interface])
               enable_service=yes])
AC_MSG_RESULT([$enable_service])
AC_SUBST([ENABLE_SERVICE],[$enable_service])
#  ----------------------------------------------------------------------
#  Compile interface with internal EXTRA library
AC_MSG_CHECKING([for EXTRA support])
AC_ARG_ENABLE([extra],
              [AC_HELP_STRING([--enable-extra],[Use the extra library [default=yes]])],
              [AS_IF([test "x$enable_extra" != 'xno'],
                     [AC_DEFINE(HAVE_LIBEXTRA,[1],[Define to 1 for EXTRA interface])
                      enable_extra=yes])],
              [AC_DEFINE(HAVE_LIBEXTRA,[1],[Define to 1 for EXTRA interface])
               enable_extra=yes])
AC_MSG_RESULT([$enable_extra])
AC_SUBST([ENABLE_EXTRA],[$enable_extra])
#  ----------------------------------------------------------------------
#  Compile interface with internal IEG library
AC_MSG_CHECKING([for IEG support])
AC_ARG_ENABLE([ieg],
              [AC_HELP_STRING([--enable-ieg],[Use the ieg library [default=yes]])],
              [AS_IF([test "x$enable_ieg" != 'xno'],
                     [AC_DEFINE(HAVE_LIBIEG,[1],[Define to 1 for IEG interface])
                      enable_ieg=yes])],
              [AC_DEFINE(HAVE_LIBIEG,[1],[Define to 1 for IEG interface])
               enable_ieg=yes])
AC_MSG_RESULT([$enable_ieg])
AC_SUBST([ENABLE_IEG],[$enable_ieg])
#  ----------------------------------------------------------------------
# At the moment, there are two possible CDI bindings
# (default for CDO) linking directly to CDI convenience library with libtool
# (default for CDI) build and link to a shared CDI library
AS_IF([test "x$CDO_DISABLE_CDILIB" = "x1"],[enable_cdi_lib=no],[enable_cdi_lib=yes])
# save CDI binding mode for later automake use
AM_CONDITIONAL([ENABLE_CDI_LIB],[test x$enable_cdi_lib = 'xyes'])
# create shell variables for the representation of configure results
AS_IF([test x$enable_cdi_lib = 'xno'],[AC_SUBST([ENABLE_CDI_LIB],[false])],[AC_SUBST([ENABLE_CDI_LIB],[true])])
#  ----------------------------------------------------------------------
#  Build a static CDI
AC_MSG_CHECKING([for building an additional static CDI binary])
AC_ARG_ENABLE([all-static],
              [AS_HELP_STRING([--enable-all-static],[build a completely statically linked CDO binary [default=no]])],
              [AS_IF([test "x$enable_all_static" != "xno"],
                     [enable_all_static=yes],
                     [enable_all_static=no])],
              [enable_all_static=no])
AC_MSG_RESULT([$enable_all_static])
AM_CONDITIONAL([ENABLE_ALL_STATIC],[test x$enable_all_static = 'xyes'])
])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl End:
