#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#if  defined  (HAVE_LIBGRIB_API)
#  include <grib_api.h>
#endif

#include <stdio.h>

#include "cdi.h"
#include "cdi_int.h"
#include "gribapi.h"
#include "dmemory.h"

#define XSTRING(x)	#x
#define STRING(x)	XSTRING(x)

static char gribapi_libvers[64] = "";
#if  defined  (HAVE_LIBGRIB_API)
static int gribapi_libvers_init;
#endif


void gribapiLibraryVersion(int* major_version, int* minor_version, int* revision_version)
{
#if  defined  (HAVE_LIBGRIB_API)
  long version = grib_get_api_version();
  (*major_version)    = (int)(version/10000);
  (*minor_version)    = (int)((version-(*major_version)*10000)/100);
  (*revision_version) = (int)(version-(*major_version)*10000-(*minor_version)*100);
#else
  (*major_version)    = 0;
  (*minor_version)    = 0;
  (*revision_version) = 0;
#endif
}

const char *gribapiLibraryVersionString(void)
{
#if  defined  (HAVE_LIBGRIB_API)
  if (!gribapi_libvers_init)
    {
      int major_version, minor_version, revision_version;

      gribapiLibraryVersion(&major_version, &minor_version, &revision_version);

      sprintf(gribapi_libvers, "%d.%d.%d", major_version, minor_version, revision_version);
      gribapi_libvers_init = 1;
    }
#endif

  return (gribapi_libvers);
}


void gribContainersNew(stream_t * streamptr)
{
  int editionNumber = 2;

  if ( streamptr->filetype == FILETYPE_GRB ) editionNumber = 1;
  (void)editionNumber;
#if  defined  (HAVE_LIBCGRIBEX)
  if ( streamptr->filetype == FILETYPE_GRB )
    {
    }
  else
#endif
    {
      int nvars = streamptr->nvars;

#if defined (GRIBCONTAINER2D)
      gribContainer_t **gribContainers;
      gribContainers = (gribContainer_t **) malloc(nvars*sizeof(gribContainer_t *));

      for ( int varID = 0; varID < nvars; ++varID )
        {
          int nlevs = streamptr->vars[varID].nlevs;
          gribContainers[varID] = (gribContainer_t *) malloc(nlevs*sizeof(gribContainer_t));

          for ( int levelID = 0; levelID < nlevs; ++levelID )
            {
              gribContainers[varID][levelID].gribHandle = gribHandleNew(editionNumber);
              gribContainers[varID][levelID].init = FALSE;
            }
	}

      streamptr->gribContainers = (void **) gribContainers;
#else
      gribContainer_t *gribContainers
        = (gribContainer_t *)xmalloc((size_t)nvars*sizeof(gribContainer_t));

      for ( int varID = 0; varID < nvars; ++varID )
        {
          gribContainers[varID].gribHandle = gribHandleNew(editionNumber);
          gribContainers[varID].init = FALSE;
	}

      streamptr->gribContainers = (void *) gribContainers;
#endif
    }
}


void gribContainersDelete(stream_t * streamptr)
{
  if ( streamptr->gribContainers )
    {
      int nvars = streamptr->nvars;

#if defined (GRIBCONTAINER2D)
      gribContainer_t **gribContainers = (gribContainer_t **) streamptr->gribContainers;

      for ( int varID = 0; varID < nvars; ++varID )
	{
          int nlevs = streamptr->vars[varID].nlevs;
          for ( int levelID = 0; levelID < nlevs; ++levelID )
            {
              gribHandleDelete(gribContainers[varID][levelID].gribHandle);
            }
          free(gribContainers[varID]);
	}
#else
      gribContainer_t *gribContainers = (gribContainer_t *) streamptr->gribContainers;

      for ( int varID = 0; varID < nvars; ++varID )
	{
          gribHandleDelete(gribContainers[varID].gribHandle);
	}
#endif

      free(gribContainers);

      streamptr->gribContainers = NULL;
    }
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
