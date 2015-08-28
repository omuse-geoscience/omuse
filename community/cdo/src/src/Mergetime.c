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

/*
   This module contains the following operators:

      Merge      mergetime       Merge datasets sorted by date and time
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "util.h"


void *Mergetime(void *argument)
{
  int streamID1, streamID2 = CDI_UNDEFID;
  int tsID2 = 0, recID, varID, levelID;
  int vlistID1, vlistID2;
  int nfiles, fileID;
  int taxisID1, taxisID2 = CDI_UNDEFID;
  int lcopy = FALSE;
  int gridsize;
  int nmiss;
  int vdate, vtime;
  int last_vdate = -1, last_vtime = -1;
  int next_fileID;
  int skip_same_time = FALSE;
  int process_timestep;
  const char *ofilename;
  double *array = NULL;
  typedef struct
  {
    int streamID;
    int vlistID;
    int taxisID;
    int tsID;
    int vdate;
    int vtime;
    int nrecs;
  } sfile_t;
  sfile_t *sf = NULL;

  cdoInitialize(argument);

  {
    char *envstr;
    envstr = getenv("SKIP_SAME_TIME");
    if ( envstr )
      {
	int ival;
	ival = atoi(envstr);
	if ( ival == 1 )
	  {
	    skip_same_time = TRUE;
	    if ( cdoVerbose )
	      cdoPrint("Set SKIP_SAME_TIME to %d", ival);
	  }
      }
  }

  if ( UNCHANGED_RECORD ) lcopy = TRUE;

  nfiles = cdoStreamCnt() - 1;

  sf = (sfile_t*) malloc(nfiles*sizeof(sfile_t));

  for ( fileID = 0; fileID < nfiles; fileID++ )
    {
      if ( cdoVerbose ) cdoPrint("process: %s", cdoStreamName(fileID)->args);

      streamID1 = streamOpenRead(cdoStreamName(fileID));

      vlistID1 = streamInqVlist(streamID1);
      taxisID1 = vlistInqTaxis(vlistID1);

      sf[fileID].streamID = streamID1;
      sf[fileID].vlistID  = vlistID1;
      sf[fileID].taxisID  = taxisID1;
    }

  
  /* check that the contents is always the same */
  for ( fileID = 1; fileID < nfiles; fileID++ )
    vlistCompare(sf[0].vlistID, sf[fileID].vlistID, CMP_ALL);

  /* read the first time step */
  for ( fileID = 0; fileID < nfiles; fileID++ )
    {
      sf[fileID].tsID = 0;
      sf[fileID].nrecs = streamInqTimestep(sf[fileID].streamID, sf[fileID].tsID);
      if ( sf[fileID].nrecs == 0 )
	{
	  streamClose(sf[fileID].streamID);
	  sf[fileID].streamID = -1;
	}
      else
	{
	  sf[fileID].vdate = taxisInqVdate(sf[fileID].taxisID);
	  sf[fileID].vtime = taxisInqVtime(sf[fileID].taxisID);
	}
    }

  ofilename = cdoStreamName(nfiles)->args;

  if ( !cdoSilentMode && !cdoOverwriteMode )
    if ( fileExists(ofilename) )
      if ( !userFileOverwrite(ofilename) )
	cdoAbort("Outputfile %s already exists!", ofilename);

  streamID2 = streamOpenWrite(cdoStreamName(nfiles), cdoFiletype());

  if ( ! lcopy )
    {
      gridsize = vlistGridsizeMax(sf[0].vlistID);
      array = (double*) malloc(gridsize*sizeof(double));
    }

  while ( TRUE )
    {
      process_timestep = TRUE;

      next_fileID = -1;
      vdate = 0;
      vtime = 0;
      for ( fileID = 0; fileID < nfiles; fileID++ )
	{
	  if ( sf[fileID].streamID != -1 )
	    if ( next_fileID == -1 || sf[fileID].vdate < vdate ||
		 (sf[fileID].vdate == vdate && sf[fileID].vtime < vtime) )
	      {
		next_fileID = fileID;
		vdate = sf[fileID].vdate;
		vtime = sf[fileID].vtime;
	      }
	}

      fileID = next_fileID;

      if ( cdoVerbose )
	cdoPrint("nextstep = %d  vdate = %d  vtime = %d", next_fileID, vdate, vtime);

      if ( next_fileID == -1 ) break;

      if ( skip_same_time )
	if ( vdate == last_vdate && vtime == last_vtime )
	  {
	    char vdatestr[32], vtimestr[32];
	    date2str(vdate, vdatestr, sizeof(vdatestr));
	    time2str(vtime, vtimestr, sizeof(vtimestr));
	    cdoPrint("Timestep %4d in stream %d (%s %s) already exists, skipped!",
		     sf[fileID].tsID+1, sf[fileID].streamID, vdatestr, vtimestr);
	    process_timestep = FALSE;
	  }

      if ( process_timestep )
	{
	  if ( tsID2 == 0 )
	    {
	      vlistID1 = sf[0].vlistID;
	      vlistID2 = vlistDuplicate(vlistID1);
	      taxisID1 = vlistInqTaxis(vlistID1);
	      taxisID2 = taxisDuplicate(taxisID1);
	      vlistDefTaxis(vlistID2, taxisID2);
	      
	      streamDefVlist(streamID2, vlistID2);
	    }

	  last_vdate = vdate;
	  last_vtime = vtime;

	  taxisCopyTimestep(taxisID2, sf[fileID].taxisID);

	  streamDefTimestep(streamID2, tsID2);
	       
	  for ( recID = 0; recID < sf[fileID].nrecs; recID++ )
	    {
	      streamInqRecord(sf[fileID].streamID, &varID, &levelID);
	      streamDefRecord(streamID2,  varID,  levelID);
	  
	      if ( lcopy )
		{
		  streamCopyRecord(streamID2, sf[fileID].streamID); 
		}
	      else
		{
		  streamReadRecord(sf[fileID].streamID, array, &nmiss);
		  streamWriteRecord(streamID2, array, nmiss);
		}
	    }

	  tsID2++;
	}

      sf[fileID].nrecs = streamInqTimestep(sf[fileID].streamID, ++sf[fileID].tsID);
      if ( sf[fileID].nrecs == 0 )
	{
	  streamClose(sf[fileID].streamID);
	  sf[fileID].streamID = -1;
	}
      else
	{
	  sf[fileID].vdate = taxisInqVdate(sf[fileID].taxisID);
	  sf[fileID].vtime = taxisInqVtime(sf[fileID].taxisID);
	}
    }

  streamClose(streamID2);

  if ( ! lcopy )
    if ( array ) free(array);

  if ( sf ) free(sf);

  cdoFinish();

  return (0);
}
