#if defined(HAVE_CONFIG_H)
#  include "config.h"
#endif

#if defined(HAVE_HDF5_H)
#include <hdf5.h>
#endif

#if defined(HAVE_ZLIB_H)
#include <zlib.h>
#endif

#if defined(HAVE_LIBXML2)
#include <libxml/xmlversion.h>
#endif

#if defined(HAVE_CURL_CURL_H)
#include <curl/curl.h>
#endif

#if defined(HAVE_PROJ_API_H)
#include <proj_api.h>
#endif

#include <stdio.h>
#include <string.h>

#include "cdo_int.h" // HAVE_OPENMP4

void printFeatures(void)
{
  fprintf(stderr, "Features:");
#if defined(HAVE_LIBPTHREAD)
  fprintf(stderr, " PTHREADS");
#endif
#if defined(_OPENMP)
  fprintf(stderr, " OpenMP");
#if defined(HAVE_OPENMP4)
  fprintf(stderr, "4");
#endif
#endif
#if  defined(HAVE_LIBHDF5)
  fprintf(stderr, " HDF5");
#endif
#if  defined(HAVE_NETCDF4)
  fprintf(stderr, " NC4");
#if  defined(HAVE_NC4HDF5)
  fprintf(stderr, "/HDF5");
#if  defined(HAVE_NC4HDF5_THREADSAFE)
  fprintf(stderr, "/threadsafe");
#endif
#endif
#endif
#if  defined(HAVE_LIBNC_DAP)
  fprintf(stderr, " OPeNDAP");
#endif
#if defined(HAVE_LIBSZ)
  fprintf(stderr, " SZ");
#endif
#if defined(HAVE_LIBZ)
  fprintf(stderr, " Z");
#endif
#if defined(HAVE_LIBJASPER)
  fprintf(stderr, " JASPER");
#endif
#if defined(HAVE_LIBUDUNITS2)
  fprintf(stderr, " UDUNITS2");
#endif
#if defined(HAVE_LIBPROJ)
  fprintf(stderr, " PROJ.4");
#endif
#if defined(HAVE_LIBXML2)
  fprintf(stderr, " XML2");
#endif
#if defined(HAVE_LIBMAGICS)
  fprintf(stderr, " MAGICS");
#endif
#if defined(HAVE_LIBDRMAA)
  fprintf(stderr, " DRMAA");
#endif
#if defined(HAVE_LIBCURL)
  fprintf(stderr, " CURL");
#endif
#if defined(HAVE_LIBFFTW3)
  fprintf(stderr, " FFTW3");
#endif
#if defined(__AVX2__)
  fprintf(stderr, " AVX2");
#elif defined(__AVX__)
  fprintf(stderr, " AVX");
#elif defined(__SSE4_2__)
  fprintf(stderr, " SSE4_2");
#elif defined(__SSE4_1__)
  fprintf(stderr, " SSE4_1");
#elif defined(__SSE3__)
  fprintf(stderr, " SSE3");
#elif defined(__SSE2__)
  fprintf(stderr, " SSE2");
#endif 
  fprintf(stderr, "\n");
}


void printLibraries(void)
{
  fprintf(stderr, "Libraries:");
#if  defined(HAVE_LIBHDF5)
  fprintf(stderr, " HDF5");
#if  defined(H5_VERS_MAJOR)
  unsigned h5h_majnum = H5_VERS_MAJOR, h5h_minnum = H5_VERS_MINOR, h5h_relnum = H5_VERS_RELEASE;
  fprintf(stderr, "/%u.%u.%u", h5h_majnum, h5h_minnum, h5h_relnum);

  unsigned h5l_majnum, h5l_minnum, h5l_relnum;
  H5get_libversion(&h5l_majnum, &h5l_minnum, &h5l_relnum);
  if ( (h5h_majnum != h5l_majnum) || (h5h_minnum != h5l_minnum) || (h5h_relnum != h5l_relnum) )
    fprintf(stderr, "(%u.%u.%u)", h5l_majnum, h5l_minnum, h5l_relnum);
#endif
#endif
  /*
#if defined(HAVE_LIBZ)
  {
    fprintf(stderr, " zlib/%s", zlibVersion());
#if defined(ZLIB_VERSION)
    if ( strcmp(ZLIB_VERSION, zlibVersion()) != 0 )
      fprintf(stderr, "(h%s)", ZLIB_VERSION);
#else
    fprintf(stderr, "(header not found)");
#endif
  }
#endif
  */
#if defined(HAVE_LIBPROJ)
  fprintf(stderr, " proj");
#if defined(PJ_VERSION)
  fprintf(stderr, "/%g", PJ_VERSION*0.01);
#endif
#endif

#if defined(HAVE_LIBXML2)
  fprintf(stderr, " xml2");
#if defined(LIBXML_DOTTED_VERSION)
  fprintf(stderr, "/%s", LIBXML_DOTTED_VERSION);
#endif
#endif

#if defined(HAVE_LIBCURL)
  {
    curl_version_info_data *version_data = curl_version_info(CURLVERSION_NOW);
    fprintf(stderr, " curl/%s", version_data->version);
#if defined(LIBCURL_VERSION)
    if ( strcmp(LIBCURL_VERSION, version_data->version) != 0 )
      fprintf(stderr, "(h%s)", LIBCURL_VERSION);
#else
    fprintf(stderr, "(header not found)");
#endif
  }
#endif

  fprintf(stderr, "\n");
}
