/* src/config.h.  Generated from config.h.in by configure.  */
/* src/config.h.in.  Generated from configure.ac by autoheader.  */

/* CDO version */
#define CDO "1.7.0rc5"

/* Compiler */
#define COMPILER "gcc -std=gnu99 -g -O2 -fopenmp "

/* Compiler version */
#define COMP_VERSION "gcc (Ubuntu 4.8.2-19ubuntu1) 4.8.2"

/* Define to 1 for DATA support */
#define ENABLE_DATA 1

/* Define to 1 if you have the <cmor.h> header file. */
/* #undef HAVE_CMOR_H */

/* Define to 1 if you have the <curl/curl.h> header file. */
/* #undef HAVE_CURL_CURL_H */

/* Define to 1 if you have the declaration of `isnan', and to 0 if you don't.
   */
#define HAVE_DECL_ISNAN 1

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define to 1 if you have the <fftw3.h> header file. */
/* #undef HAVE_FFTW3_H */

/* Define to 1 if you have the <fnmatch.h> header file. */
#define HAVE_FNMATCH_H 1

/* Define to 1 if you have the `gethostname' function. */
#define HAVE_GETHOSTNAME 1

/* Define to 1 if you have the `getpagesize' function. */
#define HAVE_GETPAGESIZE 1

/* Define to 1 if you have the `getrlimit' function. */
#define HAVE_GETRLIMIT 1

/* Define to 1 if you have the <glob.h> header file. */
#define HAVE_GLOB_H 1

/* Define to 1 if you have the <grib_api.h> header file. */
/* #undef HAVE_GRIB_API_H */

/* Define to 1 if you have the <hdf5.h> header file. */
/* #undef HAVE_HDF5_H */

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the <jasper.h> header file. */
/* #undef HAVE_JASPER_H */

/* Define to 1 for GRIB1 decoding/encoding with cgribex */
#define HAVE_LIBCGRIBEX 1

/* Define to 1 for CMOR support */
/* #undef HAVE_LIBCMOR */

/* Define to 1 for CURL support */
/* #undef HAVE_LIBCURL */

/* Define to 1 for EXTRA interface */
#define HAVE_LIBEXTRA 1

/* FFTW3 library is present if defined to 1 */
/* #undef HAVE_LIBFFTW3 */

/* Define to 1 for GRIB support */
#define HAVE_LIBGRIB 1

/* GRIB_API library is present if defined to 1 */
/* #undef HAVE_LIBGRIB_API */

/* Define to 1 for HDF5 support */
/* #undef HAVE_LIBHDF5 */

/* Define to 1 for IEG interface */
#define HAVE_LIBIEG 1

/* Define to 1 for JPEG compression for GRIB2 */
/* #undef HAVE_LIBJASPER */

/* Define to 1 if you have the `m' library (-lm). */
#define HAVE_LIBM 1

/* Define to 1 for MAGICS support */
/* #undef HAVE_LIBMAGICS */

/* Define to 1 if you have the `malloc' library (-lmalloc). */
/* #undef HAVE_LIBMALLOC */

/* Define to 1 for NETCDF OpenDAP */
/* #undef HAVE_LIBNC_DAP */

/* Define to 1 for NETCDF support */
#define HAVE_LIBNETCDF 1

/* Define to 1 for PROJ support */
/* #undef HAVE_LIBPROJ */

/* Define to 1 if you have the `pthread' library (-lpthread). */
#define HAVE_LIBPTHREAD 1

/* Define to 1 for SERVICE interface */
#define HAVE_LIBSERVICE 1

/* Define to 1 for SZIP support */
/* #undef HAVE_LIBSZ */

/* Define to 1 for UDUNITS2 support */
/* #undef HAVE_LIBUDUNITS2 */

/* Define to 1 for XML2 support */
/* #undef HAVE_LIBXML2 */

/* Define to 1 if you have the <libxml/parser.h> header file. */
/* #undef HAVE_LIBXML_PARSER_H */

/* Define to 1 if you have the <libxml/tree.h> header file. */
/* #undef HAVE_LIBXML_TREE_H */

/* Define 1 for ZLIB support */
#define HAVE_LIBZ 1

/* Define to 1 if you have the <magics_api.h> header file. */
/* #undef HAVE_MAGICS_API_H */

/* Define to 1 if you have the `mallinfo' function. */
#define HAVE_MALLINFO 1

/* Define to 1 if you have the <malloc.h> header file. */
#define HAVE_MALLOC_H 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if you have a working `mmap' system call. */
#define HAVE_MMAP 1

/* Define to 1 for NETCDF4/HDF5 support */
#define HAVE_NC4HDF5 1

/* Define to 1 for NETCDF4/HDF5 threadsafe support */
#define HAVE_NC4HDF5_THREADSAFE 1

/* Define to 1 for NETCDF2 support */
#define HAVE_NETCDF2 1

/* Define to 1 for NETCDF4 support */
#define HAVE_NETCDF4 1

/* Define to 1 if you have the <netcdf.h> header file. */
#define HAVE_NETCDF_H 1

/* Define to 1 if you have the <proj_api.h> header file. */
/* #undef HAVE_PROJ_API_H */

/* Define to 1 if you have the <pthread.h> header file. */
/* #undef HAVE_PTHREAD_H */

/* Have PTHREAD_PRIO_INHERIT. */
#define HAVE_PTHREAD_PRIO_INHERIT 1

/* Define to 1 if you have the `sqrtl' function. */
#define HAVE_SQRTL 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if `st_blksize' is a member of `struct stat'. */
#define HAVE_STRUCT_STAT_ST_BLKSIZE 1

/* Define to 1 if you have the <sys/param.h> header file. */
#define HAVE_SYS_PARAM_H 1

/* Define to 1 if you have the <sys/resource.h> header file. */
#define HAVE_SYS_RESOURCE_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/times.h> header file. */
#define HAVE_SYS_TIMES_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <szlib.h> header file. */
/* #undef HAVE_SZLIB_H */

/* Define to 1 if you have the <udunits2.h> header file. */
/* #undef HAVE_UDUNITS2_H */

/* Define to 1 if you have the <udunits2/udunits2.h> header file. */
/* #undef HAVE_UDUNITS2_UDUNITS2_H */

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to 1 if you have the <wordexp.h> header file. */
#define HAVE_WORDEXP_H 1

/* Define to 1 if you have the <zlib.h> header file. */
#define HAVE_ZLIB_H 1

/* Host name */
#define HOST_NAME "ben-VirtualBox"

/* Define to the sub-directory in which libtool stores uninstalled libraries.
   */
#define LT_OBJDIR ".libs/"

/* Name of package */
#define PACKAGE "cdo"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "http://mpimet.mpg.de/cdo"

/* Define to the full name of this package. */
#define PACKAGE_NAME "cdo"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "cdo 1.7.0rc5"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "cdo"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "1.7.0rc5"

/* Define to necessary symbol if this constant uses a non-standard name on
   your system. */
/* #undef PTHREAD_CREATE_JOINABLE */

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* System type */
#define SYSTEM_TYPE "i686-pc-linux-gnu"

/* User name */
#define USER_NAME "ben"

/* Version number of package */
#define VERSION "1.7.0rc5"

/* Number of bits in a file offset, on hosts where this is settable. */
#define _FILE_OFFSET_BITS 64

/* Define for large files, on AIX-style hosts. */
/* #undef _LARGE_FILES */

/* Define to the equivalent of the C99 'restrict' keyword, or to
   nothing if this is not supported.  Do not define if restrict is
   supported directly.  */
#define restrict __restrict
/* Work around a bug in Sun C++: it does not support _Restrict or
   __restrict__, even though the corresponding Sun C compiler ends up with
   "#define restrict _Restrict" or "#define restrict __restrict__" in the
   previous line.  Perhaps some future version of Sun C++ will work with
   restrict; if so, hopefully it defines __RESTRICT like Sun C does.  */
#if defined __SUNPRO_CC && !defined __RESTRICT
# define _Restrict
# define __restrict__
#endif
