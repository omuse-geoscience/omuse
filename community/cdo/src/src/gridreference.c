/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2015 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#if defined(HAVE_CONFIG_H)
#  include "config.h"
#endif

#if defined(HAVE_LIBCURL)
#include <curl/curl.h>
#include <errno.h>
#endif

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"

/*
 * callback function for curl for writing the network retrieved grid file
 */
static
size_t write_data(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
  size_t written;
  written = fwrite(ptr, size, nmemb, stream);
  return written;
}

/* code from grid_tools.2 */
int download_gridfile(const char *restrict uri, const char *restrict basename)
{
  int rval = 1;
#if defined(HAVE_LIBCURL)
  // As curl_easy_init calls non-thread safe curl_global_init the libcurl developer advice
  // to call curl_global_init first and before potential thread spawning.

  CURLcode ret;  
  CURL *hd;
  double length;
  int status;
  int curlflags = CURL_GLOBAL_DEFAULT;

#if defined(CURL_GLOBAL_ACK_EINTR)
  curlflags |= CURL_GLOBAL_ACK_EINTR;
#endif

  ret = curl_global_init(curlflags);
  if(ret != 0)
    {
      fprintf(stderr, "ERROR: %s!\n", curl_easy_strerror(ret));
      return -1;
    }

  hd = curl_easy_init();
  if (hd == NULL)
    {
      fprintf(stderr, "ERROR: could not get curl handler.\n");
      return -1;
    }
  else
    {
      FILE *fp;
      fp = fopen(basename, "w");
      if (fp == NULL)
	{
	  fprintf(stderr, "ERROR: could not open local output file %s. %s.\n", basename, strerror(errno));
	  return -1;
	}

      //curl_easy_setopt(hd, CURLOPT_VERBOSE, 1);
      curl_easy_setopt(hd, CURLOPT_URL, uri);
      curl_easy_setopt(hd, CURLOPT_WRITEFUNCTION, write_data);
      curl_easy_setopt(hd, CURLOPT_WRITEDATA, fp);
      ret = curl_easy_perform(hd);
      fclose(fp);
      if ( ret == 0 ) 
	{
	  /*
	  int ihead;
	  curl_easy_getinfo(hd, CURLINFO_HEADER_SIZE, &ihead);
	  printf("ihead %d\n", ihead);
	  */
	  char *ctype;
	  curl_easy_getinfo(hd, CURLINFO_CONTENT_TYPE, &ctype);

	  if ( strstr(ctype, "html") == NULL ) // no html content
	    {
	      curl_easy_getinfo(hd, CURLINFO_SIZE_DOWNLOAD, &length);
	      if ( cdoVerbose ) cdoPrint("File %s downloaded - size: %.0lf byte", basename, length); 
	      rval = 0;
	    }
	  else
	    {
	      status = remove(basename);
	      if (status == -1) perror(basename);
	      if ( cdoVerbose ) cdoPrint("The requested URL was not found on this server!");
	    }
	}
      else
	{
	  status = remove(basename);
	  if (status == -1) perror(basename);
	  fprintf(stderr, "ERROR: %s. Download %s failed.\n\n", curl_easy_strerror(ret), basename);
	}
      
      curl_easy_cleanup(hd);
    }
#else
  cdoWarning("CURL support not compiled in!");
#endif  

  return rval;
}

#if defined(HAVE_SYS_STAT_H)
#include <sys/stat.h>
#endif

/*
 * Search for filename.
 */
int search_file(const char *restrict directory, const char *restrict filename)
{
#if defined(HAVE_SYS_STAT_H)
  struct stat buf;
  int status;

  status = stat(directory, &buf);

  if ( status == 0 )
    {
      status = stat(filename, &buf);
      if ( status == 0 ) return 0;
    }
  else
    {
      perror(directory);
    }
#endif

  return 1;
}


int referenceToGrid(int gridID1)
{
  int gridID2 = -1;
  int gridsize;
  char griduri[8912];
  char gridpath[8912];

  griduri[0] = 0;
  gridpath[0] = 0;

  if ( gridInqReference(gridID1, NULL) ) gridInqReference(gridID1, griduri);

  if ( griduri[0] == 0 )
    {
      cdoWarning("Reference to horizontal grid not available!");
    }
  else
    {
      int lgriduri = TRUE;
      int status;
      int streamID;
      int number, position;

      char *basename = strrchr(griduri, '/');
      if ( basename == NULL )
	{
	  basename = griduri;
	  lgriduri = FALSE;
	}
      else
	{
	  basename++;
	}

      strcpy(gridpath, "./");
      strcat(gridpath, basename);
      if ( cdoVerbose ) cdoPrint("Search for horizontal grid file \"%s\"", gridpath);
  
      /* scan local directory for file */
      status = search_file("./", gridpath);
      if ( status != 0 )
	{
	  if ( cdoGridSearchDir != NULL)
	    {
	      strcpy(gridpath, cdoGridSearchDir);
	      strcat(gridpath, basename);
	      if ( cdoVerbose ) cdoPrint("Search for horizontal grid file \"%s\"", gridpath);
  
	      /* scan directory given by environment variable */
	      status = search_file(cdoGridSearchDir, gridpath);
	    }

	  if ( status != 0 && lgriduri )
	    {
	      /*
	      strcpy(griduri, "http://icon-downloads.mpimet.mpg.de/grids/public/edzw/icon_grid_0001x_R02B05_R.nc");
	      char *basename = strrchr(griduri, '/') + 1;
	      */
	      if ( cdoVerbose ) cdoPrint("Download horizontal grid file %s to %s", griduri, basename);
	      status = download_gridfile(griduri, basename);
	    }
	}

      if ( status == 0 )
	{
	  if ( cdoVerbose ) cdoPrint("Horizontal grid file used: %s", gridpath);
      
	  gridsize = gridInqSize(gridID1);

	  number = gridInqNumber(gridID1);
	  position = gridInqPosition(gridID1);

	  streamID = streamOpenRead(gridpath);
	  if ( streamID < 0 ) cdiOpenError(streamID, "Open failed on horizontal grid file >%s<", gridpath);

	  int vlistID, gridID = -1;
	  int ngrids;
	  vlistID = streamInqVlist(streamID);
	  ngrids = vlistNgrids(vlistID);
	  if ( position > 0 && position <= ngrids )
	    {
	      gridID = vlistGrid(vlistID, position-1);
	      if ( gridInqSize(gridID) == gridsize )
		gridID2 = gridDuplicate(gridID);
	      else
		cdoWarning("Grid size %d on position %d do not match! Reference=%s", gridsize, position, gridpath);
	    }
	  else if ( position == 0 )
	    {
	      for ( int grididx = 0; grididx < ngrids; ++grididx )
		{
		  gridID = vlistGrid(vlistID, grididx);
		  if ( gridInqSize(gridID) == gridsize )
		    {
		      gridID2 = gridDuplicate(gridID);
		      break;
			}
		}
	    }
	  else
	    cdoWarning("Number of grid in reference %d not available! Reference=%s", position, gridpath);
	  
	  streamClose(streamID);
	}

      if ( gridID2 != -1 )
	{
	  unsigned char uuidOfHGrid1[CDI_UUID_SIZE];
	  unsigned char uuidOfHGrid2[CDI_UUID_SIZE];

	  memset(uuidOfHGrid1, 0, 16);
	  memset(uuidOfHGrid2, 0, 16);

	  gridInqUUID(gridID1, uuidOfHGrid1);
	  gridInqUUID(gridID2, uuidOfHGrid2);
	  
	  if ( uuidOfHGrid1[0] != 0 && uuidOfHGrid1[0] != 0 && memcmp(uuidOfHGrid1, uuidOfHGrid2, 16) != 0 )
	    cdoWarning("UUID of horizontal grids differ!");

	  int number1 = gridInqNumber(gridID1);
	  int number2 = gridInqNumber(gridID2);

	  if ( number1 > 0 && number2 > 0 && number1 != number2 )
	    cdoWarning("Number of grid used of horizontal grids differ!");
	}
    }

  return (gridID2);
}
