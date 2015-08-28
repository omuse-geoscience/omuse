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

      Ymonarith  ymonadd         Add multi-year monthly time series
      Ymonarith  ymonsub         Subtract multi-year monthly time series
      Ymonarith  ymonmul         Multiply multi-year monthly time series
      Ymonarith  ymondiv         Divide multi-year monthly time series
      Ymonarith  yseasadd        Add multi-year seasonal time series
      Ymonarith  yseassub        Subtract multi-year seasonal time series
      Ymonarith  yseasmul        Multiply multi-year seasonal time series
      Ymonarith  yseasdiv        Divide multi-year seasonal time series
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


#define  MAX_MON    12


void *Ymonarith(void *argument)
{
  enum {MONTHLY, SEASONAL};
  int nrecs, nvars, nlev, recID;
  int tsID;
  int varID, levelID;
  int offset;
  int vdate, year, mon, day;
  field_t field1, field2;
  int **varnmiss2[MAX_MON];
  double **vardata2[MAX_MON];
  const char *seas_name[4];

  cdoInitialize(argument);

  cdoOperatorAdd("ymonadd",  func_add, MONTHLY, NULL);
  cdoOperatorAdd("ymonsub",  func_sub, MONTHLY, NULL);
  cdoOperatorAdd("ymonmul",  func_mul, MONTHLY, NULL);
  cdoOperatorAdd("ymondiv",  func_div, MONTHLY, NULL);
  cdoOperatorAdd("yseasadd", func_add, SEASONAL, NULL);
  cdoOperatorAdd("yseassub", func_sub, SEASONAL, NULL);
  cdoOperatorAdd("yseasmul", func_mul, SEASONAL, NULL);
  cdoOperatorAdd("yseasdiv", func_div, SEASONAL, NULL);

  int operatorID = cdoOperatorID();
  int operfunc = cdoOperatorF1(operatorID);
  int opertype = cdoOperatorF2(operatorID);

  int streamID1 = streamOpenRead(cdoStreamName(0));
  int streamID2 = streamOpenRead(cdoStreamName(1));

  int vlistID1 = streamInqVlist(streamID1);
  int vlistID2 = streamInqVlist(streamID2);
  int vlistID3 = vlistDuplicate(vlistID1);

  vlistCompare(vlistID1, vlistID2, CMP_ALL);

  int gridsize = vlistGridsizeMax(vlistID1);

  field_init(&field1);
  field_init(&field2);
  field1.ptr = (double*) malloc(gridsize*sizeof(double));
  field2.ptr = (double*) malloc(gridsize*sizeof(double));

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = vlistInqTaxis(vlistID2);
  int taxisID3 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID3, taxisID3);

  int streamID3 = streamOpenWrite(cdoStreamName(2), cdoFiletype());

  streamDefVlist(streamID3, vlistID3);

  nvars  = vlistNvars(vlistID2);

  if ( opertype == SEASONAL ) get_season_name(seas_name);

  for ( mon = 0; mon < MAX_MON ; mon++ ) vardata2[mon] = NULL;

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID2, tsID)) )
    {
      vdate = taxisInqVdate(taxisID2);

      cdiDecodeDate(vdate, &year, &mon, &day);
      if ( mon < 1 || mon > MAX_MON ) cdoAbort("Month %d out of range!", mon);
      mon--;

      if ( opertype == SEASONAL ) mon = month_to_season(mon+1);

      if ( vardata2[mon] != NULL )
	{
	  if ( opertype == SEASONAL )
	    cdoAbort("Season %s already allocatd!", seas_name[mon]);
	  else
	    cdoAbort("Month %d already allocatd!", mon);
	}

      vardata2[mon]  = (double **) malloc(nvars*sizeof(double *));
      varnmiss2[mon] = (int **) malloc(nvars*sizeof(int *));

      for ( varID = 0; varID < nvars; varID++ )
	{
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
	  nlev     = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
	  vardata2[mon][varID]  = (double*) malloc(nlev*gridsize*sizeof(double));
	  varnmiss2[mon][varID] = (int*) malloc(nlev*sizeof(int));
	}

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID2, &varID, &levelID);

	  gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
	  offset   = gridsize*levelID;

	  streamReadRecord(streamID2, vardata2[mon][varID]+offset, &field2.nmiss);
	  varnmiss2[mon][varID][levelID] = field2.nmiss;
	}

      tsID++;
    }


  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      vdate = taxisInqVdate(taxisID1);

      cdiDecodeDate(vdate, &year, &mon, &day);
      if ( mon < 1 || mon > MAX_MON ) cdoAbort("Month %d out of range!", mon);
      mon--;

      if ( opertype == SEASONAL ) mon = month_to_season(mon+1);

      if ( vardata2[mon] == NULL )
	{
	  if ( opertype == SEASONAL )
	    cdoAbort("Season %s not found!", seas_name[mon]);
	  else
	    cdoAbort("Month %d not found!", mon);
	}

      taxisCopyTimestep(taxisID3, taxisID1);

      streamDefTimestep(streamID3, tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamReadRecord(streamID1, field1.ptr, &field1.nmiss);

	  gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
	  offset   = gridsize*levelID;

	  memcpy(field2.ptr, vardata2[mon][varID]+offset, gridsize*sizeof(double));
	  field2.nmiss = varnmiss2[mon][varID][levelID];

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

  for ( mon = 0; mon < MAX_MON ; mon++ ) 
    if ( vardata2[mon] )
      {
	for ( varID = 0; varID < nvars; varID++ )
	  {
	    free(vardata2[mon][varID]);
	    free(varnmiss2[mon][varID]);
	  }

	free(vardata2[mon]);
	free(varnmiss2[mon]);
      }

  if ( field1.ptr ) free(field1.ptr);
  if ( field2.ptr ) free(field2.ptr);

  cdoFinish();

  return (0);
}
