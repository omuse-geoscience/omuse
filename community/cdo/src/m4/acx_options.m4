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
ENABLE_NC4HDF5=no
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
                                          AC_MSG_RESULT([yes])],[AC_MSG_RESULT([no])])
			           AC_MSG_CHECKING([netcdf's nc4/hdf5 support])
                                   AS_IF([test "x$($NC_CONFIG --has-hdf5)" = "xyes"],
                                         [AC_DEFINE([HAVE_NC4HDF5],[1],[Define to 1 for NETCDF4/HDF5 support])
                                          ENABLE_NC4HDF5=yes
                                          AC_MSG_RESULT([yes])],[AC_MSG_RESULT([no])]) ],
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
                                          AC_MSG_RESULT([yes])],[AC_MSG_RESULT([no])])
			           AC_MSG_CHECKING([netcdf's nc4/hdf5 support])
                                   AS_IF([test "x$($NC_CONFIG --has-hdf5)" = "xyes"],
                                         [AC_DEFINE([HAVE_NC4HDF5],[1],[Define to 1 for NETCDF4/HDF5 support])
                                          ENABLE_NC4HDF5=yes
                                          AC_MSG_RESULT([yes])],[AC_MSG_RESULT([no])]) ],
                                   [AC_MSG_RESULT([Could not find nc-config! go on with default configuration])])],
                                [AC_MSG_NOTICE([$with_netcdf is not a directory! NETCDF suppressed])])])],
            [AC_MSG_CHECKING([for NETCDF library])
             AC_MSG_RESULT([suppressed])])

AS_IF([test "x$ENABLE_NC4HDF5" = "xyes"],
      [AC_SEARCH_LIBS([H5TS_mutex_lock], [netcdf],
               [AC_DEFINE([HAVE_NC4HDF5_THREADSAFE],[1],[Define to 1 for NETCDF4/HDF5 threadsafe support])],,
	       [-lhdf5_hl -lhdf5])])

AC_SUBST([ENABLE_NETCDF])
AC_SUBST([ENABLE_NC2])
AC_SUBST([ENABLE_NC4])
AC_SUBST([ENABLE_NC4HDF5])
AC_SUBST([NETCDF_ROOT])
AC_SUBST([NETCDF_INCLUDE])
AC_SUBST([NETCDF_LIBS])
#  ----------------------------------------------------------------------
#  Link application with CMOR library
CMOR_LIBS=''
AC_ARG_WITH([cmor],
            [AS_HELP_STRING([--with-cmor=<directory>],
                            [Specify location of CMOR library.])],
            [AS_CASE(["$with_cmor"],
                     [no],[AC_MSG_CHECKING([for cmor library])
                           AC_MSG_RESULT([suppressed])],
                     [yes],[AC_CHECK_HEADERS([cmor.h])
                            AC_SEARCH_LIBS([cmor_load_table],[cmor],[AC_DEFINE([HAVE_LIBCMOR],[1],[Define to 1 for CMOR support])],
                                           [AC_MSG_ERROR([Could not link to cmor library!])])
                            AC_SUBST([CMOR_LIBS],[" -lcmor"])],
                     [*],[CMOR_ROOT=$with_cmor
                          AS_IF([test -d "$CMOR_ROOT"],
                                [LDFLAGS="$LDFLAGS -L$CMOR_ROOT/lib"
                                 CPPFLAGS="$CPPFLAGS -I$CMOR_ROOT/include"
                                 AC_SEARCH_LIBS([cmor_load_table],
                                                [cmor],
                                                [AC_DEFINE([HAVE_LIBCMOR],[1],[Define to 1 for CMOR support])],
                                                [AC_MSG_ERROR([Could not link to cmor library!])])
                                 CMOR_LIBS=" -L$CMOR_ROOT/lib -lcmor"],
                                [AC_MSG_ERROR([$CMOR_ROOT is not a directory! CMOR suppressed])])])],
            [AC_MSG_CHECKING([for the CMOR library])
             AC_MSG_RESULT([suppressed])])
AC_SUBST([CMOR_LIBS])
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
                                 AC_SEARCH_LIBS([jas_init],
                                                [jasper],
                                                [AC_DEFINE([HAVE_LIBJASPER],[1],[Define to 1 for JPEG compression for GRIB2])],
                                                [AC_MSG_ERROR([Could not link to jasper library! Required for GRIB_API])])
                                 JASPER_LIBS=" -L$JASPER_ROOT/lib -ljasper"],
                                [AC_MSG_ERROR([$JASPER_ROOT is not a directory! JASPER suppressed])])])],
            [AC_MSG_CHECKING([for the JASPER library])
             AC_MSG_RESULT([suppressed])])
AC_SUBST([JASPER_LIBS])
#  ----------------------------------------------------------------------
#  Compile application with GRIB_API library (for GRIB2 support)
GRIB_API_INCLUDE=''
GRIB_API_LIBS=''
ENABLE_GRIBAPI=no
AC_ARG_WITH([grib_api],
            [AS_HELP_STRING([--with-grib_api=<yes|no|directory>],
                            [library for grib2 compression; if a directory is given, it will be used as a value for --with-jasper-root])],
            [AS_CASE(["$with_grib_api"],
                     [no],[AC_MSG_CHECKING([for GRIB_API library])
                           AC_MSG_RESULT([suppressed])],
                     [yes],[AC_CHECK_HEADERS([grib_api.h])
                            AC_SEARCH_LIBS([grib_get_message],
                                           [grib_api],
                                           [AC_DEFINE([HAVE_LIBGRIB_API],[1],[GRIB_API library is present if defined to 1])
                                            ENABLE_GRIBAPI=yes],
                                           [AC_MSG_ERROR([Could not link to grib_api library])])],
                     [*],[GRIB_API_ROOT=$with_grib_api
                          AS_IF([test -d "$GRIB_API_ROOT"],
                                [LDFLAGS="-L$GRIB_API_ROOT/lib $LDFLAGS"
                                 CPPFLAGS="-I$GRIB_API_ROOT/include $CPPFLAGS"
                                 AC_CHECK_HEADERS([grib_api.h])
                                 AC_SEARCH_LIBS([grib_get_message],
                                                [grib_api],
                                                [AC_DEFINE([HAVE_LIBGRIB_API],[1],[GRIB_API library is present if defined to 1])
                                                 ENABLE_GRIBAPI=yes],
                                                [AC_MSG_ERROR([Could not link to grib_api library])])
                                 GRIB_API_LIBS=" -L$GRIB_API_ROOT/lib -lgrib_api"
                                 GRIB_API_INCLUDE=" -I$GRIB_API_ROOT/include"],
                                [AC_MSG_ERROR([$GRIB_API_ROOT is not a directory! GRIB_API suppressed])])])],
            [AC_MSG_CHECKING([for the GRIB_API library])
             AC_MSG_RESULT([suppressed])])
AC_SUBST([GRIB_API_INCLUDE])
AC_SUBST([GRIB_API_LIBS])
AC_SUBST([ENABLE_GRIBAPI])
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
#  Compile with fftw support
AC_MSG_CHECKING([for FFTW3 support])
AC_ARG_WITH([fftw3],
    [AS_HELP_STRING([--with-fftw3],
      [enable support for fftw3])],
    [],
    [with_fftw3=no])

  AS_IF([test "x$with_fftw3" != xno],
      [AC_CHECK_HEADERS([fftw3.h])
    AC_SEARCH_LIBS([fftw_cleanup],[fftw3],
      [AC_DEFINE([HAVE_LIBFFTW3],[1],[FFTW3 library is present if defined to 1])],
      [AC_MSG_RESULT([Could not link to fftw3 library])])])

#  ----------------------------------------------------------------------
#  Checks for PROJ.4 library
AC_ARG_WITH([proj],
            [AS_HELP_STRING([--with-proj=<directory>],
                            [Specify location of PROJ library for cartographic projections.])],
            [AS_CASE(["$with_proj"],
                     [no],[AC_MSG_CHECKING([for proj library])
                           AC_MSG_RESULT([suppressed])],
                     [yes],[AC_CHECK_HEADERS([proj_api.h])
                            AC_SEARCH_LIBS([pj_init],[proj],[AC_DEFINE([HAVE_LIBPROJ],[1],[Define to 1 for PROJ support])],
                                           [AC_MSG_ERROR([Could not link to PROJ library!])])
                            AC_SUBST([PROJ_LDFLAGS],[" -lproj"])
                            AC_SUBST([PROJ_INCLUDE],[""])],
                     [*],[PROJ_ROOT=$with_proj
                          AS_IF([test -d "$PROJ_ROOT"],
                                [LDFLAGS="-L$PROJ_ROOT/lib $LDFLAGS"
                                 CPPFLAGS="-I$PROJ_ROOT/include $CPPFLAGS"
                                 AC_CHECK_HEADERS([proj_api.h])
                                 AC_SEARCH_LIBS([pj_init],
                                                [proj],
                                                [AC_DEFINE([HAVE_LIBPROJ],[1],[Define to 1 for PROJ support])],
                                                [AC_MSG_ERROR([Could not link to PROJ library!])])
                                 AC_SUBST([PROJ_LDFLAGS],[" -L$PROJ_ROOT/lib -lproj"])
                                 AC_SUBST([PROJ_INCLUDE],[" -I$PROJ_ROOT/include"])],
                                [AC_MSG_ERROR([$PROJ_ROOT is not a directory! PROJ suppressed])])])],
            [AC_MSG_CHECKING([for the PROJ library])
             AC_MSG_RESULT([suppressed])])
#  ----------------------------------------------------------------------
#  Checks for CURL library
AC_ARG_WITH([curl],
            [AS_HELP_STRING([--with-curl=<directory>],
                            [Specify location of CURL library.])],
            [AS_CASE(["$with_curl"],
                     [no],[AC_MSG_CHECKING([for curl library])
                           AC_MSG_RESULT([suppressed])],
                     [yes],[AC_CHECK_HEADERS([curl/curl.h])
                            AC_SEARCH_LIBS([curl_global_init],[curl],[AC_DEFINE([HAVE_LIBCURL],[1],[Define to 1 for CURL support])],
                                           [AC_MSG_ERROR([Could not link to CURL library!])])
                            AC_SUBST([CURL_LDFLAGS],[" -lcurl"])
                            AC_SUBST([CURL_INCLUDE],[""])],
                     [*],[CURL_ROOT=$with_curl
                          AS_IF([test -d "$CURL_ROOT"],
                                [LDFLAGS="-L$CURL_ROOT/lib $LDFLAGS"
                                 CPPFLAGS="-I$CURL_ROOT/include $CPPFLAGS"
                                 AC_CHECK_HEADERS([curl/curl.h])
                                 AC_SEARCH_LIBS([curl_global_init],
                                                [curl],
                                                [AC_DEFINE([HAVE_LIBCURL],[1],[Define to 1 for CURL support])],
                                                [AC_MSG_ERROR([Could not link to CURL library!])])
                                 AC_SUBST([CURL_LDFLAGS],[" -L$CURL_ROOT/lib -lcurl"])
                                 AC_SUBST([CURL_INCLUDE],[" -I$CURL_ROOT/include"])],
                                [AC_MSG_ERROR([$CURL_ROOT is not a directory! CURL suppressed])])])],
            [AC_MSG_CHECKING([for the CURL library])
             AC_MSG_RESULT([suppressed])])
#  ----------------------------------------------------------------------
#  Link application with UDUNITS2 library
AC_ARG_WITH([udunits2],
            [AS_HELP_STRING([--with-udunits2=<directory>],
                            [Specify location of UDUNITS2 library.])],
            [AS_CASE(["$with_udunits2"],
                     [no],[AC_MSG_CHECKING([for udunits2 library])
                           AC_MSG_RESULT([suppressed])],
                     [yes],[AC_CHECK_HEADERS([udunits2.h])
                            AC_SEARCH_LIBS([ut_parse],
                                           [udunits2],
                                           [AC_DEFINE([HAVE_LIBUDUNITS2],[1],[Define to 1 for UDUNITS2 support])],
                                           [AC_MSG_ERROR([Could not link to udunits2 library!])])
                            AC_SUBST([UDUNITS_LDFLAGS],[" -ludunits2"])
                            AC_SUBST([UDUNITS_INCLUDE],[""])],
                     [*],[UDUNITS_ROOT=$with_udunits2
                          AS_IF([test -d "$UDUNITS_ROOT"],
                                [LDFLAGS="$LDFLAGS -L$UDUNITS_ROOT/lib"
                                 CPPFLAGS="$CPPFLAGS -I$UDUNITS_ROOT/include"
                                 AC_CHECK_HEADERS([udunits2.h])
                                 AC_CHECK_HEADERS([udunits2/udunits2.h])
                                 AC_SEARCH_LIBS([ut_parse],
                                                [udunits2],
                                                [AC_DEFINE([HAVE_LIBUDUNITS2],[1],[Define to 1 for UDUNITS2 support])],
                                                [AC_MSG_ERROR([Could not link to udunits2 library!])])
                                 AC_SUBST([UDUNITS_LDFLAGS],[" -L$UDUNITS_ROOT/lib -ludunits2"])
                                 AC_SUBST([UDUNITS_INCLUDE],[" -I$UDUNITS_ROOT/include"])],
                                [AC_MSG_ERROR([$UDUNITS_ROOT is not a directory! UDUNITS2 suppressed])])])],
            [AC_MSG_CHECKING([for the UDUNITS2 library])
             AC_MSG_RESULT([suppressed])])
#  ----------------------------------------------------------------------
#  Compile application with MAGICS (xml required)
MAGICS_ROOT=''
MAGICS_INCLUDE=''
MAGICS_LIBS=''
AC_ARG_WITH([magics],
            [AS_HELP_STRING([--with-magics=<yes|no|directory>],[location of magics library (lib and include subdirs)])],
            [AS_CASE(["$with_magics"],
                     [no],[AC_MSG_CHECKING([for magics library])
                           AC_MSG_RESULT([suppressed])],
                     [yes],[AC_CHECK_HEADERS([magics_api.h])
                            AC_SEARCH_LIBS([mag_open],
                                           [MagPlus],
                                           [AC_DEFINE([HAVE_LIBMAGICS],[1],[Define to 1 for MAGICS support])],
                                           [AC_MSG_ERROR([Could not link to magics library])])
                            AC_SUBST([MAGICS_LIBS],[" -lMagPlus"])],
                     [*],[AS_IF([test -d "$with_magics"],
                                [MAGICS_ROOT=$with_magics
                                 LDFLAGS="-L$MAGICS_ROOT/lib $LDFLAGS"
                                 CPPFLAGS="-I$MAGICS_ROOT/include/magics $CPPFLAGS"
                                 AC_CHECK_HEADERS([magics_api.h])
                                 AC_SEARCH_LIBS([mag_open],
                                                [MagPlus],
                                                [AC_DEFINE([HAVE_LIBMAGICS],[1],[Define to 1 for MAGICS support])],
                                                [AC_MSG_ERROR([Could not link to magics library])])
                                 MAGICS_LIBS=" -L$MAGICS_ROOT/lib -lMagPlus"
                                 MAGICS_INCLUDE=" -I$MAGICS_ROOT/include/magics"],
                                [AC_MSG_NOTICE([$with_magics is not a directory! MAGICS suppressed])])])],
            [AC_MSG_CHECKING([for MAGICS library])
             AC_MSG_RESULT([suppressed])])
AC_SUBST([MAGICS_ROOT])
AC_SUBST([MAGICS_INCLUDE])
AC_SUBST([MAGICS_LIBS])
AC_ARG_WITH([libxml2],
            [AS_HELP_STRING([--with-libxml2=<yes|no|directory>],[location of libxml2 library (lib and include subdirs)])],
            [AS_CASE(["$with_libxml2"],
                     [no],[AC_MSG_CHECKING([for libxml2 library])
                           AC_MSG_RESULT([suppressed])],
                     [yes],[CPPFLAGS="$CPPFLAGS -I/usr/include/libxml2"
                            AC_CHECK_HEADERS([libxml/parser.h])
                            AC_CHECK_HEADERS([libxml/tree.h])
                            AC_SEARCH_LIBS([xmlInitMemory],
                                           [xml2],
                                           [AC_DEFINE([HAVE_LIBXML2],[1],[Define to 1 for XML2 support])],
                                           [AC_MSG_ERROR([Could not link to libxml2 library])])
                            AC_SUBST([XML2_LIBS],[" -lxml2"])],
                     [*],[AS_IF([test -d "$with_libxml2"],
                                [XML2_ROOT=$with_libxml2
                                 LDFLAGS="-L$XML2_ROOT/lib $LDFLAGS"
                                 CPPFLAGS="-I$XML2_ROOT/include/libxml2 $CPPFLAGS"
                                 AC_CHECK_HEADERS([libxml/parser.h])
                                 AC_CHECK_HEADERS([libxml/tree.h])
                                 AC_SEARCH_LIBS([xmlInitMemory],
                                                [xml2],
                                                [AC_DEFINE([HAVE_LIBXML2],[1],[Define to 1 for XML2 support])],
                                                [AC_MSG_ERROR([Could not link to libxml2 library])])
                                 XML2_LIBS=" -L$XML2_ROOT/lib -lxml2"
                                 XML2_INCLUDE=" -I$XML2_ROOT/include/libxml2"],
                                [AC_MSG_NOTICE([$with_libxml2 is not a directory! XML2 suppressed])])])],
            [AC_MSG_CHECKING([for XML2 library])
             AC_MSG_RESULT([suppressed])])
AM_CONDITIONAL([ENABLE_MAGICS],[test ! x$with_magics = 'x' -a ! x$with_libxml2 = 'x'])
#  ----------------------------------------------------------------------
#  How to build CDI into CDO? 
INTERNAL_CDI_DIR=libcdi
# At the moment, there are two possible CDI bindings
# (default)          linking directly to CDI's object files, i.e. a libtool
#                    convenience library
# (--enable-cdi-lib) build and link to a shared CDI library
AC_MSG_CHECKING([for build a separate CDI library and link CDO to it])
AC_ARG_ENABLE([cdi-lib],
              [AS_HELP_STRING([--enable-cdi-lib],[build, install and link to a CDI library [default=no]])],
              [AS_IF([test "x$enable_cdi_lib" != "xno"],
                     [enable_cdi_lib=yes],
                     [enable_cdi_lib=no
                      CDO_DISABLE_CDILIB=1
                      export CDO_DISABLE_CDILIB])],
              [enable_cdi_lib=no
               CDO_DISABLE_CDILIB=1
               export CDO_DISABLE_CDILIB])
AC_MSG_RESULT([$enable_cdi_lib])
# save CDI binding mode for later automake use
AM_CONDITIONAL([ENABLE_CDI_LIB],[test x$enable_cdi_lib = 'xyes'])
# create shell variables for the representation of configure results
AS_IF([test x$enable_cdi_lib = 'xno'],[AC_SUBST([ENABLE_CDI_LIB],[false])],[AC_SUBST([ENABLE_CDI_LIB],[true])])
# scan for CDI as a subproject
AC_CONFIG_SUBDIRS([libcdi])
#  ----------------------------------------------------------------------
#  Build a static CDO
AC_MSG_CHECKING([for building an additional static CDO binary])
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
