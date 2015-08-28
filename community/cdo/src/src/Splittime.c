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

      Splittime  splithour       Split hours
      Splittime  splitday        Split days
      Splittime  splitmon        Split months
      Splittime  splitseas       Split seasons
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "util.h"

#include <time.h>

#define  MAX_STREAMS 32

struct tm datetime_to_tm(int date, int time)
{
  int year, month, day, hour, minute, second;
  cdiDecodeDate(date, &year, &month, &day);
  cdiDecodeTime(time, &hour, &minute, &second);

  struct tm stime;
  memset(&stime, 0, sizeof(struct tm));

  stime.tm_sec  = second;
  stime.tm_min  = minute;
  stime.tm_hour = hour;
  stime.tm_mday = day;
  stime.tm_mon  = month-1;
  stime.tm_year = year-1900;

  return stime;
}

void *Splittime(void *argument)
{
  int streamID2;
  int varID;
  int nrecs;
  int tsID, recID, levelID;
  int  streamIDs[MAX_STREAMS], tsIDs[MAX_STREAMS];
  char filesuffix[32];
  char filename[8192];
  int index = 0;
  int i;
  int vdate, vtime;
  int lcopy = FALSE;
  int gridsize;
  int nmiss;
  int gridID;
  int nlevel;
  double *array = NULL;
  field_t **vars = NULL;
  const char *format = NULL;

  cdoInitialize(argument);

  if ( processSelf() != 0 ) cdoAbort("This operator can't be combined with other operators!");

  if ( UNCHANGED_RECORD ) lcopy = TRUE;

  int SPLITHOUR = cdoOperatorAdd("splithour", func_time, 10000, NULL);
  int SPLITDAY  = cdoOperatorAdd("splitday",  func_date,     1, NULL);
  int SPLITMON  = cdoOperatorAdd("splitmon",  func_date,   100, NULL);
  int SPLITSEAS = cdoOperatorAdd("splitseas", func_date,   100, NULL);

  UNUSED(SPLITDAY);
  UNUSED(SPLITHOUR);

  int operatorID = cdoOperatorID();
  int operfunc   = cdoOperatorF1(operatorID);
  int operintval = cdoOperatorF2(operatorID);

  if ( operatorID == SPLITMON )
    {
      if ( operatorArgc() == 1 ) format = operatorArgv()[0];
    }

  const char *seas_name[4];
  get_season_name(seas_name);

  for ( i = 0; i < MAX_STREAMS; i++ ) streamIDs[i] = -1;
  for ( i = 0; i < MAX_STREAMS; i++ ) tsIDs[i] = 0;

  int streamID1 = streamOpenRead(cdoStreamName(0));

  int vlistID1 = streamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  strcpy(filename, cdoStreamName(1)->args);
  int nchars = strlen(filename);

  const char *refname = cdoStreamName(0)->argv[cdoStreamName(0)->argc-1];
  filesuffix[0] = 0;
  cdoGenFileSuffix(filesuffix, sizeof(filesuffix), streamInqFiletype(streamID1), vlistID1, refname);

  //  if ( ! lcopy )
    {
      gridsize = vlistGridsizeMax(vlistID1);
      if ( vlistNumber(vlistID1) != CDI_REAL ) gridsize *= 2;
      array = (double*) malloc(gridsize*sizeof(double));
    }

  int nvars = vlistNvars(vlistID1);
  int nconst = 0;
  for ( varID = 0; varID < nvars; varID++ )
    if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) nconst++;

  if ( nconst )
    {
      vars = (field_t **) malloc(nvars*sizeof(field_t *));

      for ( varID = 0; varID < nvars; varID++ )
	{
	  if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT )
	    {
	      gridID  = vlistInqVarGrid(vlistID1, varID);
	      nlevel  = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      gridsize = gridInqSize(gridID);
		  
	      vars[varID] = (field_t*) malloc(nlevel*sizeof(field_t));

	      for ( levelID = 0; levelID < nlevel; levelID++ )
		{
		  field_init(&vars[varID][levelID]);
		  vars[varID][levelID].grid = gridID;
		  vars[varID][levelID].ptr  = (double*) malloc(gridsize*sizeof(double));
		}
	    }
	}
    }

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      vdate = taxisInqVdate(taxisID1);
      vtime = taxisInqVtime(taxisID1);

      if ( operfunc == func_date )
	{
	  index = (vdate/operintval)%100;
	  if ( index < 0 ) index = -index;

	  if ( operatorID == SPLITSEAS ) index = month_to_season(index);
	}
      else if ( operfunc == func_time )
	{
	  index = (vtime/operintval)%100;
	}

      if ( index < 0 || index >= MAX_STREAMS )
	cdoAbort("Index out of range!");

      streamID2 = streamIDs[index];
      if ( streamID2 < 0 )
	{
	  if ( operatorID == SPLITSEAS )
	    {
	      sprintf(filename+nchars, "%3s", seas_name[index]);
	      if ( filesuffix[0] )
		sprintf(filename+nchars+3, "%s", filesuffix);
	    }
	  else
	    {
	      size_t slen;
	      char oformat[32];
	      strcpy(oformat, "%02d");

	      if ( operatorID == SPLITMON && format )
		{
		  char sbuf[32];
		  struct tm stime = datetime_to_tm(vdate, vtime);
		  slen = strftime(sbuf, 32, format, &stime);

		  if ( slen ) strcpy(oformat, sbuf);
		}

	      slen = sprintf(filename+nchars, oformat, index);
	      if ( filesuffix[0] )
		sprintf(filename+nchars+slen, "%s", filesuffix);
	    }

	  if ( cdoVerbose ) cdoPrint("create file %s", filename);

	  argument_t *fileargument = file_argument_new(filename);
	  streamID2 = streamOpenWrite(fileargument, cdoFiletype());
	  file_argument_free(fileargument);

	  streamDefVlist(streamID2, vlistID2);

	  streamIDs[index] = streamID2;
	}

      taxisCopyTimestep(taxisID2, taxisID1);
      streamDefTimestep(streamID2, tsIDs[index]);

      if ( tsID > 0 && tsIDs[index] == 0 && nconst )
	{
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT )
		{
		  nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
		  for ( levelID = 0; levelID < nlevel; levelID++ )
		    {
		      streamDefRecord(streamID2, varID, levelID);
		      nmiss = vars[varID][levelID].nmiss;
		      streamWriteRecord(streamID2, vars[varID][levelID].ptr, nmiss);
		    }
		}
	    }
	}

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamDefRecord(streamID2,  varID,  levelID);

	  if ( lcopy && !(tsID == 0 && nconst) )
	    {
	      streamCopyRecord(streamID2, streamID1);
	    }
	  else
	    {
	      streamReadRecord(streamID1, array, &nmiss);
	      streamWriteRecord(streamID2, array, nmiss);

	      if ( tsID == 0 && nconst )
		{
		  if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT )
		    {
		      gridID  = vlistInqVarGrid(vlistID1, varID);
		      gridsize = gridInqSize(gridID);
		      memcpy(vars[varID][levelID].ptr, array, gridsize*sizeof(double));
		      vars[varID][levelID].nmiss = nmiss;
		    }
		}
	    }
	}

      tsIDs[index]++;
      tsID++;
    }

  streamClose(streamID1);

  for ( index = 0; index < MAX_STREAMS; index++ )
    {
      streamID2 = streamIDs[index];
      if ( streamID2 >= 0 ) streamClose(streamID2);
    }
 
  if ( array ) free(array);

  if ( nconst )
    {
      for ( varID = 0; varID < nvars; varID++ )
	{
	  if ( vlistInqVarTsteptype(vlistID2, varID) == TSTEP_CONSTANT )
	    {
	      nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
	      for ( levelID = 0; levelID < nlevel; levelID++ )
		if ( vars[varID][levelID].ptr )
		  free(vars[varID][levelID].ptr);

	      free(vars[varID]);
	    }
	}

      if ( vars ) free(vars);
    }

  vlistDestroy(vlistID2);

  cdoFinish();

  return (0);
}
