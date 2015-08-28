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

      Transpose  transxy         Transpose X/Y
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void transxy(int gridID, double *array1, double *array2)
{
  int i, j;
  int nx, ny;
  double **a2D1, **a2D2;

  nx = gridInqXsize(gridID);
  ny = gridInqYsize(gridID);

  a2D1 = (double **) malloc(ny*sizeof(double *));
  a2D2 = (double **) malloc(nx*sizeof(double *));

  for ( j = 0; j < ny; ++j ) a2D1[j] = array1+j*nx;
  for ( i = 0; i < nx; ++i ) a2D2[i] = array2+i*ny;

  for ( j = 0; j < ny; ++j )
    for ( i = 0; i < nx; ++i )
      a2D2[i][j] = a2D1[j][i];

  free(a2D1);
  free(a2D2);
}


void *Transpose(void *argument)
{
  int streamID1, streamID2;
  int gridsize;
  int ngrids, index;
  int gridID1, gridID2;
  int nx, ny;
  int nrecs, recID;
  int gridID, tsID;
  int varID, levelID;
  int vlistID1, vlistID2;
  int nmiss;
  int taxisID1, taxisID2;
  double *array1, *array2;

  cdoInitialize(argument);

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  ngrids = vlistNgrids(vlistID1);
  for ( index = 0; index < ngrids; index++ )
    {
      gridID1 = vlistGrid(vlistID1, index);
      nx = gridInqXsize(gridID1);
      ny = gridInqYsize(gridID1);

      gridID2 = gridCreate(GRID_GENERIC, nx*ny);
      gridDefXsize(gridID2, ny);
      gridDefYsize(gridID2, nx);

      vlistChangeGridIndex(vlistID2, index, gridID2);
    }

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  gridsize = vlistGridsizeMax(vlistID1);

  array1 = (double*) malloc(gridsize*sizeof(double));
  array2 = (double*) malloc(gridsize*sizeof(double));

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamReadRecord(streamID1, array1, &nmiss);

	  gridID = vlistInqVarGrid(vlistID1, varID);

	  transxy(gridID, array1, array2);

	  streamDefRecord(streamID2, varID, levelID);
	  streamWriteRecord(streamID2, array2, nmiss);
	}
      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  free(array1);
  free(array2);

  cdoFinish();

  return (0);
}
