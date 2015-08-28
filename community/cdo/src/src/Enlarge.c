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

      Enlarge    enlarge         Enlarge fields
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Enlarge(void *argument)
{
  int streamID1, streamID2;
  int vlistID1, vlistID2;
  int gridID1, gridID2;
  int i, index, ngrids;
  int recID, nrecs;
  int tsID, varID, levelID;
  int taxisID1, taxisID2;
  int nmiss;
  int gridsize1, gridsize2;
  int xsize1, ysize1, xsize2, ysize2;
  int ix, iy;
  int linfo = TRUE;
  double missval;
  double *array1, *array2;

  cdoInitialize(argument);

  operatorCheckArgc(1);

  gridID2 = cdoDefineGrid(operatorArgv()[0]);
  xsize2 = gridInqXsize(gridID2);
  ysize2 = gridInqYsize(gridID2);

  if ( cdoVerbose ) fprintf(stderr, "gridID2 %d, xsize2 %d, ysize2 %d\n", gridID2, xsize2, ysize2);

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  gridsize2 = gridInqSize(gridID2);
  if ( gridsize2 < vlistGridsizeMax(vlistID1) )
    cdoAbort("Gridsize of input stream is greater than new gridsize!");

  ngrids = vlistNgrids(vlistID1);
  for ( index = 0; index < ngrids; index++ )
    {
      vlistChangeGridIndex(vlistID2, index, gridID2);
    }

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  array1 = (double*) malloc(gridsize2*sizeof(double));
  array2 = (double*) malloc(gridsize2*sizeof(double));

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamReadRecord(streamID1, array1, &nmiss);

	  missval = vlistInqVarMissval(vlistID1, varID);
	  gridID1 = vlistInqVarGrid(vlistID1, varID);
	  xsize1 = gridInqXsize(gridID1);
	  ysize1 = gridInqYsize(gridID1);
	  gridsize1 = gridInqSize(gridID1);

	  if ( xsize1 == 0 ) xsize1 = 1;
	  if ( ysize1 == 0 ) ysize1 = 1;

	  /* printf("%d %d %d %d\n", xsize1, ysize1, xsize2, ysize2); */
	  if ( xsize1 == 1 && ysize1 == ysize2 && xsize1*ysize1 == gridsize1 )
	    {
	      if ( linfo )
		{
		  cdoPrint("Enlarge zonal");
		  linfo = FALSE;
		}
	      
	      for ( iy = 0; iy < ysize2; iy++ )
		for ( ix = 0; ix < xsize2; ix++ )
		  array2[ix+iy*xsize2] = array1[iy];

	      if ( nmiss ) nmiss *= xsize2;
	    }
	  else if ( ysize1 == 1 && xsize1 == xsize2 && xsize1*ysize1 == gridsize1 )
	    {
	      if ( linfo )
		{
		  cdoPrint("Enlarge meridional");
		  linfo = FALSE;
		}
	      
	      for ( iy = 0; iy < ysize2; iy++ )
		for ( ix = 0; ix < xsize2; ix++ )
		  array2[ix+iy*xsize2] = array1[ix];

	      if ( nmiss ) nmiss *= ysize2;
	    }
	  else
	    {
	      memcpy(array2, array1, gridsize1*sizeof(double));
	      for ( i = gridsize1; i < gridsize2; i++ )
		{
		  array2[i] = array1[gridsize1-1];
		}

	      if ( nmiss && DBL_IS_EQUAL(array1[gridsize1-1], missval) ) nmiss += (gridsize2 - gridsize1);
	    }
	    
	  streamDefRecord(streamID2, varID,  levelID);
	  streamWriteRecord(streamID2, array2, nmiss);
	}

      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  if ( array1 ) free(array1);
  if ( array2 ) free(array2);

  cdoFinish();

  return (0);
}
