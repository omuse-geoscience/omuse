#ifndef PIO_WRITE_H
#define PIO_WRITE_H

#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <stdbool.h>

struct model_config
{
  int nlon, nlat, nts, max_nlev, nvars;
  int filetype, datatype;
  bool compute_checksum;
  const char *suffix;
};

#endif
