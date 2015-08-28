#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <stdio.h>

#include "error.h"
#include "cdi.h"
#include "basetime.h"

#undef  UNDEFID
#define UNDEFID  CDI_UNDEFID

void basetimeInit(basetime_t *basetime)
{
  if ( basetime == NULL )
    Error("Internal problem! Basetime not allocated.");

  basetime->ncvarid       = UNDEFID;
  basetime->ncdimid       = UNDEFID;
  basetime->ncvarboundsid = UNDEFID;
  basetime->leadtimeid    = UNDEFID;
  basetime->lwrf          = 0;
  basetime->timevar_cache = NULL;
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
