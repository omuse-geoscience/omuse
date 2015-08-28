#ifndef PIO_CDF_INT_H
#define PIO_CDF_INT_H

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#ifdef HAVE_LIBNETCDF
#include "cdf_int.h"

void
cdiPioEnableNetCDFParAccess();

#endif
#endif
