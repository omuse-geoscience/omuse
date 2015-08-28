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

      writegrid Write grid
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "grid.h"


void *Writegrid(void *argument)
{
  int streamID;
  int vlistID;
  int gridID;
  long gridsize, i;
  int gridtype;
  int *mask = NULL;

  cdoInitialize(argument);

  streamID = streamOpenRead(cdoStreamName(0));

  vlistID = streamInqVlist(streamID);
  gridID  = vlistGrid(vlistID, 0);

  gridtype = gridInqType(gridID);
  gridsize = gridInqSize(gridID);

  if ( gridtype == GRID_GME ) gridID = gridToUnstructured(gridID, 1);

  if ( gridtype != GRID_CURVILINEAR && gridtype != GRID_UNSTRUCTURED )
    gridID = gridToCurvilinear(gridID, 1);

  if ( gridInqXbounds(gridID, NULL) == 0 || gridInqYbounds(gridID, NULL) == 0 )
    cdoAbort("Grid corner missing!");

  mask = (int*) malloc(gridsize*sizeof(int));

  if ( gridInqMask(gridID, NULL) )
    {
      gridInqMask(gridID, mask);
    }
  else
    {
      for ( i = 0; i < gridsize; i++ ) mask[i] = 1;
    }
      
  writeNCgrid(cdoStreamName(1)->args, gridID, mask);

  streamClose(streamID);

  if ( mask ) free(mask);

  cdoFinish();

  return (0);
}
