#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "cdi.h"

const char *cdiStringError(int cdiErrno)
{
  static const char UnknownError[] = "Unknown Error";
  static const char _EUFTYPE[]     = "Unsupported file type";
  static const char _ELIBNAVAIL[]  = "Unsupported file type (library support not compiled in)";
  static const char _EUFSTRUCT[]   = "Unsupported file structure";
  static const char _EUNC4[]       = "Unsupported netCDF4 structure";
  static const char _ELIMIT[]      = "Internal limits exceeded";

  switch (cdiErrno) {
  case CDI_ESYSTEM:
    {
      const char *cp = strerror(errno);
      if ( cp == NULL ) break;
      return cp;
    }
  case CDI_EUFTYPE:    return _EUFTYPE;
  case CDI_ELIBNAVAIL: return _ELIBNAVAIL;
  case CDI_EUFSTRUCT:  return _EUFSTRUCT;
  case CDI_EUNC4:      return _EUNC4;
  case CDI_ELIMIT:     return _ELIMIT;
  }

  return UnknownError;
}

/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
