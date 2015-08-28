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

      Arithlat   mulcoslat       Multiply with cos(lat)
      Arithlat   divcoslat       Divide by cos(lat)
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "grid.h"


void *Arithlat(void *argument)
{
  int operatorID;
  int operfunc;
  int streamID1, streamID2;
  int gridtype;
  int gridID, gridID0 = -1;
  int nrecs, recID;
  int tsID;
  int varID, levelID;
  int vlistID1, vlistID2;
  int taxisID1, taxisID2;
  int nmiss;
  long gridsize, i;
  char units[CDI_MAX_NAME];
  double *scale = NULL;
  double *array = NULL;

  cdoInitialize(argument);

  cdoOperatorAdd("mulcoslat", func_mul, 0, NULL);
  cdoOperatorAdd("divcoslat", func_div, 0, NULL);

  operatorID = cdoOperatorID();
  operfunc = cdoOperatorF1(operatorID);

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  gridsize = vlistGridsizeMax(vlistID1);

  array = (double*) malloc(gridsize*sizeof(double));

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamReadRecord(streamID1, array, &nmiss);
	  
	  gridID = vlistInqVarGrid(vlistID1, varID);

	  if ( gridID != gridID0 )
	    {
	      gridID0 = gridID;

	      gridtype = gridInqType(gridID);
	      if ( gridtype == GRID_LONLAT      ||
		   gridtype == GRID_GAUSSIAN    ||
		   gridtype == GRID_LCC )
		{
		  gridID = gridToCurvilinear(gridID, 0);
		}
	      else if ( gridtype == GRID_CURVILINEAR ||
			gridtype == GRID_UNSTRUCTURED )
		{
		  /* No conversion necessary */
		}
	      else if ( gridtype == GRID_GME )
		{
		  gridID = gridToUnstructured(gridID, 0);
		}
	      else
		{
		  if ( gridtype == GRID_GAUSSIAN_REDUCED )
		    cdoAbort("Unsupported grid type: %s, use CDO option -R to convert reduced to regular grid!",
			     gridNamePtr(gridtype));
		  else
		    cdoAbort("Unsupported grid type: %s", gridNamePtr(gridtype));
		}

	      gridsize = gridInqSize(gridID);

	      scale = (double*) realloc(scale, gridsize*sizeof(double));
	      gridInqYvals(gridID, scale);

	      /* Convert lat/lon units if required */
	      
	      gridInqXunits(gridID, units);

	      grid_to_radian(units, gridsize, scale, "grid latitudes");

	      if ( operfunc == func_mul )
		for ( i = 0; i < gridsize; ++i ) scale[i] = cos(scale[i]);
	      else
		for ( i = 0; i < gridsize; ++i ) scale[i] = 1./cos(scale[i]);

	      if ( cdoVerbose ) for ( i = 0; i < 10; ++i ) cdoPrint("coslat  %3d  %g", i+1, scale[i]);
	    }

	  for ( i = 0; i < gridsize; ++i ) array[i] *= scale[i];

	  streamDefRecord(streamID2, varID, levelID);
	  streamWriteRecord(streamID2, array, nmiss);
	}

      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  if ( array ) free(array);
  if ( scale ) free(scale);

  cdoFinish();

  return (0);
}
