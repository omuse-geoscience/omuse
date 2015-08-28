#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>

#include "cdf.h"
#include "cdi.h"
#include "cdi_int.h"
#include "cdf_int.h"


const char *cdfLibraryVersion(void)
{
#if  defined  (HAVE_LIBNETCDF)
  return (nc_inq_libvers());
#else
  return ("library undefined");
#endif
}

#if  defined(HAVE_LIBHDF5)
#if defined(__cplusplus)
extern "C" {
#endif
  int H5get_libversion(unsigned *, unsigned *, unsigned *);
#if defined(__cplusplus)
}
#endif
#endif

const char *hdfLibraryVersion(void)
{
#if  defined(HAVE_LIBHDF5)
  static char hdf_libvers[256];
  unsigned majnum, minnum, relnum;

  H5get_libversion(&majnum, &minnum, &relnum);

  sprintf(hdf_libvers, "%u.%u.%u", majnum, minnum, relnum);

  return (hdf_libvers);
#else
  return ("library undefined");
#endif
}


int CDF_Debug   = 0;    /* If set to 1, debugging           */


void cdfDebug(int debug)
{
  CDF_Debug = debug;

  if ( CDF_Debug )
    Message("debug level %d", debug);
}

#if  defined  (HAVE_LIBNETCDF)
static
void cdfComment(int ncid)
{
  static char comment[256] = "Climate Data Interface version ";
  static int init = 0;

  if ( ! init )
    {
      init = 1;
      const char *libvers = cdiLibraryVersion();
      const char *blank = strchr(libvers, ' ');
      size_t size = blank ? (size_t)(blank - libvers) : 0;

      if ( size == 0 || ! isdigit((int) *libvers) )
	strcat(comment, "??");
      else
	strncat(comment, libvers, size);
      strcat(comment, " (http://mpimet.mpg.de/cdi)");
    }

  cdf_put_att_text(ncid, NC_GLOBAL, "CDI", strlen(comment), comment);
  cdf_put_att_text(ncid, NC_GLOBAL, "Conventions", 6, "CF-1.4");
}
#endif

int cdfOpenFile(const char *filename, const char *mode, int *filetype)
{
  int ncid = -1;
#if  defined  (HAVE_LIBNETCDF)
  int fmode = tolower(*mode);
  int writemode = NC_CLOBBER;
  int readmode = NC_NOWRITE;
  int status;

  if ( filename == NULL )
    ncid = CDI_EINVAL;
  else
    {
      switch (fmode)
	{
	case 'r':
	  status = cdf_open(filename, readmode, &ncid);
	  if ( status > 0 && ncid < 0 ) ncid = CDI_ESYSTEM;
#if  defined  (HAVE_NETCDF4)
	  else
	    {
	      int format;
	      (void) nc_inq_format(ncid, &format);
	      if ( format == NC_FORMAT_NETCDF4_CLASSIC )
		{
		  *filetype = FILETYPE_NC4C;
		}
	    }
#endif
	  break;
	case 'w':
#if  defined  (NC_64BIT_OFFSET)
	  if      ( *filetype == FILETYPE_NC2  ) writemode |= NC_64BIT_OFFSET;
#endif
#if  defined  (HAVE_NETCDF4)
	  if      ( *filetype == FILETYPE_NC4  ) writemode |= NC_NETCDF4;
	  else if ( *filetype == FILETYPE_NC4C ) writemode |= NC_NETCDF4 | NC_CLASSIC_MODEL;
#endif
	  cdf_create(filename, writemode, &ncid);
	  cdfComment(ncid);
	  break;
	case 'a':
	  cdf_open(filename, NC_WRITE, &ncid);
	  break;
	default:
	  ncid = CDI_EINVAL;
	}
    }
#endif

  return (ncid);
}


int cdfOpen(const char *filename, const char *mode)
{
  int fileID = 0;
  int filetype = FILETYPE_NC;

  if ( CDF_Debug )
    Message("Open %s with mode %c", filename, *mode);

  fileID = cdfOpenFile(filename, mode, &filetype);

  if ( CDF_Debug )
    Message("File %s opened with id %d", filename, fileID);

  return (fileID);
}


int cdfOpen64(const char *filename, const char *mode)
{
  int fileID = -1;
  int open_file = TRUE;
  int filetype = FILETYPE_NC2;

  if ( CDF_Debug )
    Message("Open %s with mode %c", filename, *mode);

#if  defined  (HAVE_LIBNETCDF)
#if  ! defined  (NC_64BIT_OFFSET)
  open_file = FALSE;
#endif
#endif

  if ( open_file )
    {
      fileID = cdfOpenFile(filename, mode, &filetype);

      if ( CDF_Debug )
	Message("File %s opened with id %d", filename, fileID);
    }
  else
    {
      fileID = CDI_ELIBNAVAIL;
    }

  return (fileID);
}


int cdf4Open(const char *filename, const char *mode, int *filetype)
{
  int fileID = -1;
  int open_file = FALSE;

  if ( CDF_Debug )
    Message("Open %s with mode %c", filename, *mode);

#if  defined  (HAVE_NETCDF4)
  open_file = TRUE;
#endif

  if ( open_file )
    {
      fileID = cdfOpenFile(filename, mode, filetype);

      if ( CDF_Debug )
	Message("File %s opened with id %d", filename, fileID);
    }
  else
    {
      fileID = CDI_ELIBNAVAIL;
    }

  return (fileID);
}


void cdfCloseFile(int fileID)
{
#if  defined  (HAVE_LIBNETCDF)
  cdf_close(fileID);
#endif
}

void cdfClose(int fileID)
{
  cdfCloseFile(fileID);
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
