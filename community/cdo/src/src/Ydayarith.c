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

      Ydayarith  ydayadd         Add multi-year daily time series
      Ydayarith  ydaysub         Subtract multi-year daily time series
      Ydayarith  ydaymul         Multiply multi-year daily time series
      Ydayarith  ydaydiv         Divide multi-year daily time series
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


#define  MAX_DAY   1232

void *Ydayarith(void *argument)
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
  int vdate, year, mon, day;
  field_t field1, field2;
  int **varnmiss2[MAX_DAY];
  double **vardata2[MAX_DAY];

  cdoInitialize(argument);

  cdoOperatorAdd("ydayadd", func_add, 0, NULL);
  cdoOperatorAdd("ydaysub", func_sub, 0, NULL);
  cdoOperatorAdd("ydaymul", func_mul, 0, NULL);
  cdoOperatorAdd("ydaydiv", func_div, 0, NULL);

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

  for ( day = 0; day < MAX_DAY ; day++ ) vardata2[day] = NULL;

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID2, tsID)) )
    {
      vdate = taxisInqVdate(taxisID2);

      cdiDecodeDate(vdate, &year, &mon, &day);
      day += mon*100;
      if ( day < 0 || day >= MAX_DAY ) cdoAbort("Day %d out of range!", day);

      if ( vardata2[day] != NULL ) cdoAbort("Day %d already allocatd!", day);

      vardata2[day]  = (double **) malloc(nvars*sizeof(double *));
      varnmiss2[day] = (int **) malloc(nvars*sizeof(int *));

      for ( varID = 0; varID < nvars; varID++ )
	{
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
	  nlev     = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
	  vardata2[day][varID]  = (double*) malloc(nlev*gridsize*sizeof(double));
	  varnmiss2[day][varID] = (int*) malloc(nlev*sizeof(int));
	}

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID2, &varID, &levelID);

	  gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
	  offset   = gridsize*levelID;

	  streamReadRecord(streamID2, vardata2[day][varID]+offset, &field2.nmiss);
	  varnmiss2[day][varID][levelID] = field2.nmiss;
	}

      tsID++;
    }


  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      vdate = taxisInqVdate(taxisID1);

      cdiDecodeDate(vdate, &year, &mon, &day);
      day += mon*100;
      if ( day < 0 || day >= MAX_DAY ) cdoAbort("Day %d out of range!", day);

      taxisCopyTimestep(taxisID3, taxisID1);

      streamDefTimestep(streamID3, tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamReadRecord(streamID1, field1.ptr, &field1.nmiss);

	  gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
	  offset   = gridsize*levelID;
	  if ( vardata2[day] == NULL ) cdoAbort("Day %d not found!", day);
	  memcpy(field2.ptr, vardata2[day][varID]+offset, gridsize*sizeof(double));
	  field2.nmiss = varnmiss2[day][varID][levelID];

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

  for ( day = 0; day < MAX_DAY ; day++ ) 
    if ( vardata2[day] )
      {
	for ( varID = 0; varID < nvars; varID++ )
	  {
	    free(vardata2[day][varID]);
	    free(varnmiss2[day][varID]);
	  }

	free(vardata2[day]);
	free(varnmiss2[day]);
      }

  if ( field1.ptr ) free(field1.ptr);
  if ( field2.ptr ) free(field2.ptr);

  cdoFinish();

  return (0);
}
