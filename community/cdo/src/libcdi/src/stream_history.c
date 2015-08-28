#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include "dmemory.h"
#include "cdi.h"
#include "cdi_int.h"
#include "stream_cdf.h"


void streamDefHistory(int streamID, int length, const char *history)
{
#ifdef HAVE_LIBNETCDF
  stream_t *streamptr = stream_to_pointer(streamID);

  if ( streamptr->filetype == FILETYPE_NC  ||
       streamptr->filetype == FILETYPE_NC2 ||
       streamptr->filetype == FILETYPE_NC4 ||
       streamptr->filetype == FILETYPE_NC4C )
    {
      char *histstring;
      size_t len;
      if ( history )
	{
	  len = strlen(history);
	  if ( len )
	    {
              /* FIXME: what's the point of strdupx? Why not use
               * history argument directly? */
	      histstring = strdupx(history);
	      cdfDefHistory(streamptr, length, histstring);
	      free(histstring);
	    }
	}
    }
#else
  (void)streamID; (void)length; (void)history;
#endif
}


int streamInqHistorySize(int streamID)
{
  int size = 0;
#ifdef HAVE_LIBNETCDF
  stream_t *streamptr = stream_to_pointer(streamID);

  if ( streamptr->filetype == FILETYPE_NC  ||
       streamptr->filetype == FILETYPE_NC2 ||
       streamptr->filetype == FILETYPE_NC4 ||
       streamptr->filetype == FILETYPE_NC4C )
    {
      size = cdfInqHistorySize(streamptr);
    }
#else
  (void)streamID;
#endif
  return (size);
}


void streamInqHistoryString(int streamID, char *history)
{
#ifdef HAVE_LIBNETCDF
  stream_t *streamptr = stream_to_pointer(streamID);

  if ( streamptr->filetype == FILETYPE_NC  ||
       streamptr->filetype == FILETYPE_NC2 ||
       streamptr->filetype == FILETYPE_NC4 ||
       streamptr->filetype == FILETYPE_NC4C )
    {
      cdfInqHistoryString(streamptr, history);
    }
#else
  (void)streamID; (void)history;
#endif
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
