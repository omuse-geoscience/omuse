#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <string.h>

#include "dmemory.h"
#include "error.h"

#include "cdi.h"
#include "cdi_int.h"


static void streamvar_init_recordtable(stream_t *streamptr, int varID, int isub)
{
  streamptr->vars[varID].recordTable[isub].nlevs    = 0;
  streamptr->vars[varID].recordTable[isub].recordID = NULL;
  streamptr->vars[varID].recordTable[isub].lindex   = NULL;
}


static
void streamvar_init_entry(stream_t *streamptr, int varID)
{
  streamptr->vars[varID].ncvarid      = CDI_UNDEFID;
  streamptr->vars[varID].defmiss      = 0;

  streamptr->vars[varID].subtypeSize  = 0;
  streamptr->vars[varID].recordTable  = NULL;

  streamptr->vars[varID].gridID       = CDI_UNDEFID;
  streamptr->vars[varID].zaxisID      = CDI_UNDEFID;
  streamptr->vars[varID].tsteptype    = CDI_UNDEFID;
  streamptr->vars[varID].subtypeID    = CDI_UNDEFID;
}

static
int streamvar_new_entry(stream_t *streamptr)
{
  int varID = 0;
  int streamvarSize;
  svarinfo_t *streamvar;

  streamvarSize = streamptr->varsAllocated;
  streamvar     = streamptr->vars;
  /*
    Look for a free slot in streamvar.
    (Create the table the first time through).
  */
  if ( ! streamvarSize )
    {
      int i;

      streamvarSize = 2;
      streamvar
        = (svarinfo_t *)xmalloc((size_t)streamvarSize * sizeof(svarinfo_t));
      if ( streamvar == NULL )
	{
          Message("streamvarSize = %d", streamvarSize);
	  SysError("Allocation of svarinfo_t failed");
	}

      for ( i = 0; i < streamvarSize; i++ )
	streamvar[i].isUsed = FALSE;
    }
  else
    {
      while ( varID < streamvarSize )
	{
	  if ( ! streamvar[varID].isUsed ) break;
	  varID++;
	}
    }
  /*
    If the table overflows, double its size.
  */
  if ( varID == streamvarSize )
    {
      int i;

      streamvarSize = 2*streamvarSize;
      streamvar
        = (svarinfo_t *)xrealloc(streamvar,
                                 (size_t)streamvarSize * sizeof (svarinfo_t));
      if ( streamvar == NULL )
	{
          Message("streamvarSize = %d", streamvarSize);
	  SysError("Reallocation of svarinfo_t failed");
	}
      varID = streamvarSize/2;

      for ( i = varID; i < streamvarSize; i++ )
	streamvar[i].isUsed = FALSE;
    }

  streamptr->varsAllocated = streamvarSize;
  streamptr->vars          = streamvar;

  streamvar_init_entry(streamptr, varID);

  streamptr->vars[varID].isUsed = TRUE;
  return (varID);
}


static void
allocate_record_table_entry(stream_t *streamptr, int varID, int subID, int nlevs)
{
  int *level    = (int *)xmalloc((size_t)nlevs * sizeof (int));
  int *lindex   = (int *)xmalloc((size_t)nlevs * sizeof (int));

  for (int levID = 0; levID < nlevs; levID++ )
    {
      level[levID]    = CDI_UNDEFID;
      lindex[levID]   = levID;
    }

  streamptr->vars[varID].recordTable[subID].nlevs    = nlevs;
  streamptr->vars[varID].recordTable[subID].recordID = level;
  streamptr->vars[varID].recordTable[subID].lindex   = lindex;
}


int stream_new_var(stream_t *streamptr, int gridID, int zaxisID, int tilesetID)
{
  if ( CDI_Debug )
    Message("gridID = %d  zaxisID = %d", gridID, zaxisID);

  int varID = streamvar_new_entry(streamptr);
  int nlevs = zaxisInqSize(zaxisID);

  streamptr->nvars++;

  streamptr->vars[varID].gridID  = gridID;
  streamptr->vars[varID].zaxisID = zaxisID;

  int nsub = 1;
  if (tilesetID != CDI_UNDEFID)
    nsub = subtypeInqSize(tilesetID); /* e.g. no of tiles */
  if ( CDI_Debug )
    Message("varID %d: create %d tiles with %d level(s), zaxisID=%d", varID, nsub, nlevs,zaxisID);
  streamptr->vars[varID].recordTable = (sleveltable_t *)xmalloc((size_t)nsub * sizeof (sleveltable_t));
  if( streamptr->vars[varID].recordTable == NULL )
    SysError("Allocation of leveltable failed!");
  streamptr->vars[varID].subtypeSize = nsub;

  for (int isub=0; isub<nsub; isub++) {
    streamvar_init_recordtable(streamptr, varID, isub);
    allocate_record_table_entry(streamptr, varID, isub, nlevs);
    if ( CDI_Debug )
      Message("streamptr->vars[varID].recordTable[isub].recordID[0]=%d",
              streamptr->vars[varID].recordTable[isub].recordID[0]);
  }

  streamptr->vars[varID].subtypeID = tilesetID;

  return (varID);
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
