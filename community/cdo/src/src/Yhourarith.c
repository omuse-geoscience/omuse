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

      Yhourarith  yhouradd         Add multi-year hourly time series
      Yhourarith  yhoursub         Subtract multi-year hourly time series
      Yhourarith  yhourmul         Multiply multi-year hourly time series
      Yhourarith  yhourdiv         Divide multi-year hourly time series
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


#define  MAX_HOUR  9301  /* 31*12*25 + 1 */

static
int hour_of_year(int vdate, int vtime)
{
  int year, month, day, houroy;
  int hour, minute, second;

  cdiDecodeDate(vdate, &year, &month, &day);
  cdiDecodeTime(vtime, &hour, &minute, &second);
      
  if ( month >= 1 && month <= 12 && day >= 1 && day <=31 && hour >= 0 && hour < 24 )
    houroy = ((month-1)*31 + day - 1)*25 + hour + 1;
  else
    houroy = 0;

  if ( houroy < 0 || houroy >= MAX_HOUR )
    {
      char vdatestr[32], vtimestr[32];
      date2str(vdate, vdatestr, sizeof(vdatestr));
      time2str(vtime, vtimestr, sizeof(vtimestr));
      cdoAbort("Hour of year %d out of range (%s %s)!", houroy, vdatestr, vtimestr);
    }

  return (houroy);
}


void *Yhourarith(void *argument)
{
  int operatorID;
  int operfunc;
  int streamID1, streamID2, streamID3;
  int gridsize;
  int nrecs, nvars, nlev, recID;
  int tsID;
  int varID, levelID;
  int offset;
  int vlistID1, vlistID2, vlistID3;
  int taxisID1, taxisID2, taxisID3;
  int vdate, vtime;
  int houroy;
  field_t field1, field2;
  int **varnmiss2[MAX_HOUR];
  double **vardata2[MAX_HOUR];

  cdoInitialize(argument);

  cdoOperatorAdd("yhouradd", func_add, 0, NULL);
  cdoOperatorAdd("yhoursub", func_sub, 0, NULL);
  cdoOperatorAdd("yhourmul", func_mul, 0, NULL);
  cdoOperatorAdd("yhourdiv", func_div, 0, NULL);

  operatorID = cdoOperatorID();
  operfunc = cdoOperatorF1(operatorID);

  streamID1 = streamOpenRead(cdoStreamName(0));
  streamID2 = streamOpenRead(cdoStreamName(1));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = streamInqVlist(streamID2);
  vlistID3 = vlistDuplicate(vlistID1);

  vlistCompare(vlistID1, vlistID2, CMP_ALL);

  gridsize = vlistGridsizeMax(vlistID1);

  field_init(&field1);
  field_init(&field2);
  field1.ptr = (double*) malloc(gridsize*sizeof(double));
  field2.ptr = (double*) malloc(gridsize*sizeof(double));

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = vlistInqTaxis(vlistID2);
  taxisID3 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID3, taxisID3);

  streamID3 = streamOpenWrite(cdoStreamName(2), cdoFiletype());

  streamDefVlist(streamID3, vlistID3);

  nvars  = vlistNvars(vlistID2);

  for ( houroy = 0; houroy < MAX_HOUR ; ++houroy ) vardata2[houroy] = NULL;

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID2, tsID)) )
    {
      vdate = taxisInqVdate(taxisID2);
      vtime = taxisInqVtime(taxisID2);

      houroy = hour_of_year(vdate, vtime);
      if ( vardata2[houroy] != NULL ) cdoAbort("Hour of year %d already allocatd!", houroy);

      vardata2[houroy]  = (double **) malloc(nvars*sizeof(double *));
      varnmiss2[houroy] = (int **) malloc(nvars*sizeof(int *));

      for ( varID = 0; varID < nvars; varID++ )
	{
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
	  nlev     = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
	  vardata2[houroy][varID]  = (double*) malloc(nlev*gridsize*sizeof(double));
	  varnmiss2[houroy][varID] = (int*) malloc(nlev*sizeof(int));
	}

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID2, &varID, &levelID);

	  gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
	  offset   = gridsize*levelID;

	  streamReadRecord(streamID2, vardata2[houroy][varID]+offset, &field2.nmiss);
	  varnmiss2[houroy][varID][levelID] = field2.nmiss;
	}

      tsID++;
    }


  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      vdate = taxisInqVdate(taxisID1);
      vtime = taxisInqVtime(taxisID1);

      houroy = hour_of_year(vdate, vtime);
      if ( vardata2[houroy] == NULL ) cdoAbort("Hour of year %d not found!", houroy);

      taxisCopyTimestep(taxisID3, taxisID1);

      streamDefTimestep(streamID3, tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamReadRecord(streamID1, field1.ptr, &field1.nmiss);

	  gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
	  offset   = gridsize*levelID;
	  memcpy(field2.ptr, vardata2[houroy][varID]+offset, gridsize*sizeof(double));
	  field2.nmiss   = varnmiss2[houroy][varID][levelID];

	  field1.grid    = vlistInqVarGrid(vlistID1, varID);
	  field1.missval = vlistInqVarMissval(vlistID1, varID);

	  field2.grid    = vlistInqVarGrid(vlistID2, varID);
	  field2.missval = vlistInqVarMissval(vlistID2, varID);

	  farfun(&field1, field2, operfunc);

	  streamDefRecord(streamID3, varID, levelID);
	  streamWriteRecord(streamID3, field1.ptr, field1.nmiss);
	}

      tsID++;
    }

  streamClose(streamID3);
  streamClose(streamID2);
  streamClose(streamID1);

  for ( houroy = 0; houroy < MAX_HOUR; ++houroy )
    if ( vardata2[houroy] )
      {
	for ( varID = 0; varID < nvars; varID++ )
	  {
	    free(vardata2[houroy][varID]);
	    free(varnmiss2[houroy][varID]);
	  }

	free(vardata2[houroy]);
	free(varnmiss2[houroy]);
      }

  if ( field1.ptr ) free(field1.ptr);
  if ( field2.ptr ) free(field2.ptr);

  cdoFinish();

  return (0);
}
