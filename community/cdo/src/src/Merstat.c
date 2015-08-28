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

      Merstat    mermin          Meridional minimum
      Merstat    mermax          Meridional maximum
      Merstat    mersum          Meridional sum
      Merstat    mermean         Meridional mean
      Merstat    meravg          Meridional average
      Merstat    merstd          Meridional standard deviation
      Merstat    merstd          Meridional standard deviation [Divisor is (n-1)]
      Merstat    mervar          Meridional variance
      Merstat    mervar          Meridional variance [Divisor is (n-1)]
      Merstat    merpctl         Meridional percentiles
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "grid.h"
#include "pstream.h"


void *Merstat(void *argument)
{
  int gridID1, gridID2 = -1, lastgrid = -1;
  int wstatus = FALSE;
  int index;
  int recID, nrecs;
  int varID, levelID;
  int needWeights = FALSE;
  int pn = 0;
  char varname[CDI_MAX_NAME];

  cdoInitialize(argument);

  cdoOperatorAdd("mermin",  func_min,  0, NULL);
  cdoOperatorAdd("mermax",  func_max,  0, NULL);
  cdoOperatorAdd("mersum",  func_sum,  0, NULL);
  cdoOperatorAdd("mermean", func_mean, 0, NULL);
  cdoOperatorAdd("meravg",  func_avg,  0, NULL);
  cdoOperatorAdd("mervar",  func_var,  0, NULL);
  cdoOperatorAdd("mervar1", func_var1, 0, NULL);
  cdoOperatorAdd("merstd",  func_std,  0, NULL);
  cdoOperatorAdd("merstd1", func_std1, 0, NULL);
  cdoOperatorAdd("merpctl", func_pctl, 0, NULL);
 
  int operatorID = cdoOperatorID();
  int operfunc = cdoOperatorF1(operatorID);

  /* RQ */
  if ( operfunc == func_pctl )
    {
      operatorInputArg("percentile number");
      pn = parameter2int(operatorArgv()[0]);
      
      if ( pn < 1 || pn > 99 )
        cdoAbort("Illegal argument: percentile number %d is not in the range 1..99!", pn);
    }
  /* QR */

  if ( operfunc == func_mean || operfunc == func_avg ||
       operfunc == func_var  || operfunc == func_std ||
       operfunc == func_var1 || operfunc == func_std1 )
    needWeights = TRUE;

  int streamID1 = streamOpenRead(cdoStreamName(0));

  int vlistID1 = streamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int ngrids = vlistNgrids(vlistID1);
  int ndiffgrids = 0;
  for ( index = 1; index < ngrids; index++ )
    if ( vlistGrid(vlistID1, 0) != vlistGrid(vlistID1, index))
      ndiffgrids++;

  if ( ndiffgrids > 0 ) cdoAbort("Too many different grids!");

  index = 0;
  gridID1 = vlistGrid(vlistID1, index);

  if ( gridInqType(gridID1) == GRID_LONLAT   ||
       gridInqType(gridID1) == GRID_GAUSSIAN ||
       gridInqType(gridID1) == GRID_GENERIC )
    {
      gridID2 = gridToMeridional(gridID1);
    }
  else
    {
      cdoAbort("Unsupported gridtype: %s", gridNamePtr(gridInqType(gridID1)));
    }

  vlistChangeGridIndex(vlistID2, index, gridID2);

  int streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  gridID1 = vlistInqVarGrid(vlistID1, 0);
  int nlonmax = gridInqXsize(gridID1); /* max nlon ? */
  int lim = vlistGridsizeMax(vlistID1);

  field_t field1, field2;
  field_init(&field1);
  field_init(&field2);

  field1.ptr    = (double*) malloc(lim*sizeof(double));
  field1.weight = NULL;
  if ( needWeights )
    field1.weight = (double*) malloc(lim*sizeof(double));

  field2.ptr  = (double*) malloc(nlonmax*sizeof(double));
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

	  field1.grid = vlistInqVarGrid(vlistID1, varID);
	  if ( needWeights && field1.grid != lastgrid )
	    {
	      lastgrid = field1.grid;
	      wstatus = gridWeights(field1.grid, field1.weight);
	    }
	  if ( wstatus != 0 && tsID == 0 && levelID == 0 )
	    {
	      vlistInqVarName(vlistID1, varID, varname);
	      cdoWarning("Using constant grid cell area weights for variable %s!", varname);
	    }
	  field1.missval = vlistInqVarMissval(vlistID1, varID);
	  field2.missval = vlistInqVarMissval(vlistID1, varID);

	  if ( operfunc == func_pctl )
	    merpctl(field1, & field2, pn);
	  else  
	    merfun(field1, &field2, operfunc);

	  streamDefRecord(streamID2, varID,  levelID);
	  streamWriteRecord(streamID2, field2.ptr, field2.nmiss);
	}
      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  if ( field1.ptr )    free(field1.ptr);
  if ( field1.weight ) free(field1.weight);
  if ( field2.ptr )    free(field2.ptr);

  cdoFinish();

  return (0);
}
