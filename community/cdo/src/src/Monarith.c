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

      Monarith  monadd         Add monthly time series
      Monarith  monsub         Subtract monthly time series
      Monarith  monmul         Multiply monthly time series
      Monarith  mondiv         Divide monthly time series
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Monarith(void *argument)
{
  int operatorID;
  int operfunc;
  int streamID1, streamID2, streamID3;
  int gridsize;
  int nrecs, nrecs2, nvars, nlev, recID;
  int tsID, tsID2;
  int varID, levelID;
  int offset;
  int vlistID1, vlistID2, vlistID3;
  int taxisID1, taxisID2, taxisID3;
  int vdate;
  int yearmon1, yearmon2 = -1;
  field_t field1, field2;
  int **varnmiss2;
  double **vardata2;

  cdoInitialize(argument);

  cdoOperatorAdd("monadd", func_add, 0, NULL);
  cdoOperatorAdd("monsub", func_sub, 0, NULL);
  cdoOperatorAdd("monmul", func_mul, 0, NULL);
  cdoOperatorAdd("mondiv", func_div, 0, NULL);

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

  vardata2  = (double **) malloc(nvars*sizeof(double *));
  varnmiss2 = (int **) malloc(nvars*sizeof(int *));

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
      nlev     = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
      vardata2[varID]  = (double*) malloc(nlev*gridsize*sizeof(double));
      varnmiss2[varID] = (int*) malloc(nlev*sizeof(int));
    }

  tsID  = 0;
  tsID2 = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      vdate = taxisInqVdate(taxisID1);

      yearmon1 = vdate / 100;

      if ( yearmon1 != yearmon2 )
	{
	  int year1, mon1;

	  year1 = yearmon1/100;
	  mon1  = yearmon1 - (yearmon1/100)*100;

	  if ( cdoVerbose ) cdoPrint("Process: Year = %4d  Month = %2d", year1, mon1);

	  nrecs2 = streamInqTimestep(streamID2, tsID2);
	  if ( nrecs2 == 0 )
	    cdoAbort("Missing year=%4d mon=%2d in %s!", year1, mon1, cdoStreamName(1)->args);

	  vdate = taxisInqVdate(taxisID2);

	  yearmon2 = vdate / 100;

	  if ( yearmon1 != yearmon2 )
	    {
	      int year2, mon2;

	      year2 = yearmon2/100;
	      mon2  = yearmon2 - (yearmon2/100)*100;

	      cdoAbort("Timestep %d in %s has wrong date!\nCurrent year=%4d mon=%2d, expected year=%4d mon=%2d",
		       tsID2+1, cdoStreamName(1)->args, year2, mon2, year1, mon1);
	    }

	  for ( recID = 0; recID < nrecs2; recID++ )
	    {
	      streamInqRecord(streamID2, &varID, &levelID);

	      gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
	      offset   = gridsize*levelID;

	      streamReadRecord(streamID2, vardata2[varID]+offset, &field2.nmiss);
	      varnmiss2[varID][levelID] = field2.nmiss;
	    }

	  tsID2++;
	}

      taxisCopyTimestep(taxisID3, taxisID1);

      streamDefTimestep(streamID3, tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamReadRecord(streamID1, field1.ptr, &field1.nmiss);

	  gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
	  offset   = gridsize*levelID;
	  memcpy(field2.ptr, vardata2[varID]+offset, gridsize*sizeof(double));
	  field2.nmiss = varnmiss2[varID][levelID];

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

  for ( varID = 0; varID < nvars; varID++ )
    {
      free(vardata2[varID]);
      free(varnmiss2[varID]);
    }

  free(vardata2);
  free(varnmiss2);

  if ( field1.ptr ) free(field1.ptr);
  if ( field2.ptr ) free(field2.ptr);

  cdoFinish();

  return (0);
}
