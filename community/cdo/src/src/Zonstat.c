/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.1

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

      Zonstat    zonmin          Zonal minimum
      Zonstat    zonmax          Zonal maximum
      Zonstat    zonrange        Zonal range
      Zonstat    zonsum          Zonal sum
      Zonstat    zonmean         Zonal mean
      Zonstat    zonavg          Zonal average
      Zonstat    zonstd          Zonal standard deviation
      Zonstat    zonstd1         Zonal standard deviation [Divisor is (n-1)]
      Zonstat    zonvar          Zonal variance
      Zonstat    zonvar1         Zonal variance [Divisor is (n-1)]
      Zonstat    zonpctl         Zonal percentiles
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "grid.h"


void *Zonstat(void *argument)
{
  int gridID1 = -1, gridID2 = -1;
  int zongridID = -1;
  int index;
  int recID, nrecs;
  int varID, levelID;
  int pn = 0;

  cdoInitialize(argument);

  cdoOperatorAdd("zonmin",   func_min,   0, NULL);
  cdoOperatorAdd("zonmax",   func_max,   0, NULL);
  cdoOperatorAdd("zonrange", func_range, 0, NULL);
  cdoOperatorAdd("zonsum",   func_sum,   0, NULL);
  cdoOperatorAdd("zonmean",  func_mean,  0, NULL);
  cdoOperatorAdd("zonavg",   func_avg,   0, NULL);
  cdoOperatorAdd("zonvar",   func_var,   0, NULL);
  cdoOperatorAdd("zonvar1",  func_var1,  0, NULL);
  cdoOperatorAdd("zonstd",   func_std,   0, NULL);
  cdoOperatorAdd("zonstd1",  func_std1,  0, NULL);
  cdoOperatorAdd("zonpctl",  func_pctl,  0, NULL);

  int operatorID = cdoOperatorID();
  int operfunc = cdoOperatorF1(operatorID);

  if ( operfunc == func_pctl )
    {
      operatorInputArg("percentile number");
      pn = parameter2int(operatorArgv()[0]);
      
      if ( pn < 1 || pn > 99 )
        cdoAbort("Illegal argument: percentile number %d is not in the range 1..99!", pn);
    }

  int streamID1 = streamOpenRead(cdoStreamName(0));

  int vlistID1 = streamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int ngrids = vlistNgrids(vlistID1);
  for ( index = 0; index < ngrids; index++ )
    {
      if ( gridInqXsize(vlistGrid(vlistID1, index)) > 1 ) 
	{
	  if ( gridID1 == -1 ) gridID1 = vlistGrid(vlistID1, index);
	}
      else
	{
	  if ( zongridID == -1 ) zongridID = vlistGrid(vlistID1, index);
	}
    }

  int ndiffgrids = 0;
  for ( index = 0; index < ngrids; index++ )
    {
      if ( zongridID != -1 && zongridID == vlistGrid(vlistID1, index) ) continue;
      if ( gridID1 != vlistGrid(vlistID1, index) ) ndiffgrids++;
    }
  if ( ndiffgrids > 0 ) cdoAbort("Too many different grids!");

  if ( gridID1 != -1 )
    {
      if ( gridInqType(gridID1) == GRID_LONLAT   ||
	   gridInqType(gridID1) == GRID_GAUSSIAN ||
	   gridInqType(gridID1) == GRID_GENERIC )
	{
	  if ( zongridID != -1 && gridInqYsize(zongridID) == gridInqYsize(gridID1) )
	    gridID2 = zongridID;
	  else
	    gridID2 = gridToZonal(gridID1);
	}
      else
	{
	  cdoAbort("Unsupported gridtype: %s", gridNamePtr(gridInqType(gridID1)));
	}
    }
  else
    {
      gridID2 = zongridID;
      cdoWarning("Input stream contains only zonal data!");
    }

  for ( index = 0; index < ngrids; index++ )
    vlistChangeGridIndex(vlistID2, index, gridID2);

  int streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  gridID1 = vlistInqVarGrid(vlistID1, 0);
  int nlatmax = gridInqYsize(gridID1); /* max nlat ? */

  int lim = vlistGridsizeMax(vlistID1);

  field_t field1, field2;
  field_init(&field2);
  field_init(&field2);
  field1.ptr  = (double*) malloc(lim*sizeof(double));
  field2.ptr  = (double*) malloc(nlatmax*sizeof(double));
  field2.grid = gridID2;

  int tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamReadRecord(streamID1, field1.ptr, &field1.nmiss);

	  field1.grid    = vlistInqVarGrid(vlistID1, varID);
	  field1.missval = vlistInqVarMissval(vlistID1, varID);
	  field2.missval = vlistInqVarMissval(vlistID1, varID);

	  if ( zongridID != -1 && zongridID == field1.grid )
	    {
	      memcpy(field2.ptr, field1.ptr, nlatmax*sizeof(double));
	      field2.nmiss = field1.nmiss;
	    }
	  else
	    {
	      if ( operfunc == func_pctl )
		zonpctl(field1, & field2, pn);
	      else  
		zonfun(field1, &field2, operfunc);
	    }

	  streamDefRecord(streamID2, varID,  levelID);
	  streamWriteRecord(streamID2, field2.ptr, field2.nmiss);
	}
      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  if ( field1.ptr ) free(field1.ptr);
  if ( field2.ptr ) free(field2.ptr);

  cdoFinish();

  return (0);
}
