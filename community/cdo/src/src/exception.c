#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <errno.h>
#include <cdi.h>
#include "cdo.h"
#include "process.h"
#include "error.h"

#include "config.h"

void cdiOpenError(int cdiErrno, const char *fmt, const char *path)
{	
  printf("\n");
  set_text_color(stderr, RESET, RED);
   fprintf(stderr, "%s: ", processInqPrompt());
  reset_text_color(stderr);
  set_text_color(stderr, RESET, BLACK);
  fprintf(stderr, fmt, path);
  reset_text_color(stderr);
   fprintf(stderr, "\n");

  fprintf(stderr, "%s\n", cdiStringError(cdiErrno));

  if ( cdiErrno == CDI_ELIBNAVAIL )
    {
      int byteorder;
      int filetype = cdiGetFiletype(path, &byteorder);

      switch (filetype)
	{
	case FILETYPE_GRB:
	  {
	    break;
	  }
	case FILETYPE_GRB2:
	  {
	    fprintf(stderr, "To create a CDO application with GRIB2 support use: ./configure --with-netcdf=<GRIB_API root directory> ...\n");
	    break;
	  }
	case FILETYPE_SRV:
	  {
	    break;
	  }
	case FILETYPE_EXT:
	  {
	    break;
	  }
	case FILETYPE_IEG:
	  {
	    break;
	  }
	case FILETYPE_NC:
	case FILETYPE_NC2:
	case FILETYPE_NC4:
	case FILETYPE_NC4C:
	  {
	    const char *ncv = (filetype == FILETYPE_NC4 || filetype == FILETYPE_NC4C) ? "4" : ((filetype == FILETYPE_NC2) ? "2" : "");
#if defined HAVE_LIBNETCDF
	    fprintf(stderr, "CDO was build with a netCDF version which doesn't support netCDF%s data!\n", ncv);
#else
	    fprintf(stderr, "To create a CDO application with netCDF%s support use: ./configure --with-netcdf=<netCDF%s root directory> ...\n", ncv, ncv);
#endif
	    break;
	  }
	default:
	  {
	    break;
	  }
	}
    }
  
  if ( _ExitOnError ) exit(EXIT_FAILURE);
}

void pstreamCloseAll(void);
void cdoAbort(const char *fmt, ...)
{
  va_list args;
	
  va_start(args, fmt);

  printf("\n");
  set_text_color(stderr, RESET, RED);
   fprintf(stderr, "%s (Abort): ", processInqPrompt());
  reset_text_color(stderr);
  set_text_color(stderr, RESET, BLACK);
  vfprintf(stderr, fmt, args);
  reset_text_color(stderr);
   fprintf(stderr, "\n");

  va_end(args);

  if ( _ExitOnError ) pstreamCloseAll();
  if ( _ExitOnError ) exit(EXIT_FAILURE);
}


void cdoWarning(const char *fmt, ...)
{
  if ( _Verbose )
    {
      va_list args;

      va_start(args, fmt);

      set_text_color(stderr, BRIGHT, YELLOW);
      fprintf(stderr, "%s (Warning): ", processInqPrompt());
      reset_text_color(stderr);
      set_text_color(stderr, RESET, BLACK);
      vfprintf(stderr, fmt, args);
      reset_text_color(stderr);
      fprintf(stderr, "\n");

      va_end(args);
    }
}


void cdoPrint(const char *fmt, ...)
{
  va_list args;

  if ( ! cdoSilentMode )
    {
      va_start(args, fmt);

      set_text_color(stderr, RESET, GREEN);
      fprintf(stderr, "%s: ", processInqPrompt());
      reset_text_color(stderr);
      set_text_color(stderr, RESET, BLACK);
      vfprintf(stderr, fmt, args);
      reset_text_color(stderr);
      fprintf(stderr, "\n");

      va_end(args);
    }
}
