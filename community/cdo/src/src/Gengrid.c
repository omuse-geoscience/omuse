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

*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Gengrid(void *argument)
{
  int streamID1, streamID2, streamID3;
  int vlistID1, vlistID2, vlistID3;
  int gridID1, gridID2, gridID3;
  int zaxisID3;
  int datatype;
  int nrecs;
  int tsID, varID, levelID;
  int gridsize, i;
  int xsize, ysize;
  int nmiss1, nmiss2;
  int taxisID3;
  double *array1, *array2, *array3;
  double missval = 0;
  double xminval, xmaxval, yminval, ymaxval;

  cdoInitialize(argument);

  streamID1 = streamOpenRead(cdoStreamName(0));
  streamID2 = streamOpenRead(cdoStreamName(1));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = streamInqVlist(streamID2);

  gridID1 = vlistGrid(vlistID1, 0);
  gridID2 = vlistGrid(vlistID2, 0);

  if ( gridInqSize(gridID1) != gridInqSize(gridID2) )
    cdoAbort("Arrays have different grid size!");

  gridsize = gridInqSize(gridID1);
  xsize = gridInqXsize(gridID1);
  ysize = gridInqYsize(gridID1);

  array1 = (double*) malloc(gridsize*sizeof(double));
  array2 = (double*) malloc(gridsize*sizeof(double));
  array3 = (double*) malloc(gridsize*sizeof(double));

  tsID = 0;
  nrecs = streamInqTimestep(streamID1, tsID);
  nrecs = streamInqTimestep(streamID2, tsID);

  streamInqRecord(streamID1, &varID, &levelID);
  streamReadRecord(streamID1, array1, &nmiss1);
  streamInqRecord(streamID2, &varID, &levelID);
  streamReadRecord(streamID2, array2, &nmiss2);

  datatype = vlistInqVarDatatype(vlistID1, 0);

  streamClose(streamID2);
  streamClose(streamID1);

  if ( nmiss1 || nmiss2 ) cdoAbort("Missing values unsupported!");

  gridID3 = gridCreate(GRID_CURVILINEAR, gridsize);

  if ( cdoVerbose ) cdoPrint("xsize %d  ysize %d", xsize, ysize);
  if ( xsize*ysize != gridsize )
    cdoAbort("xsize*ysize != gridsize");

  gridDefXsize(gridID3, xsize);
  gridDefYsize(gridID3, ysize);
  gridDefXvals(gridID3, array1);
  gridDefYvals(gridID3, array2);

  if ( datatype == DATATYPE_FLT64 )
    gridDefPrec(gridID3, DATATYPE_FLT64);
  else
    gridDefPrec(gridID3, DATATYPE_FLT32);

  xminval = array1[0];
  xmaxval = array1[0];
  yminval = array2[0];
  ymaxval = array2[0];
  for ( i = 1; i < gridsize; ++i )
    {
      if ( array1[i] < xminval ) xminval = array1[i];
      if ( array1[i] > xmaxval ) xmaxval = array1[i];
      if ( array2[i] < yminval ) yminval = array2[i];
      if ( array2[i] > ymaxval ) ymaxval = array2[i];
    }

  if ( cdoVerbose )
    cdoPrint("xminval = %g, xmaxval = %g, yminval = %g, ymaxval = %g",
	     xminval, xmaxval, yminval, ymaxval);

  /* check units */
  if ( xminval > -4 && xmaxval < 8 && yminval > -2 && ymaxval < 2 )
    {
      gridDefXunits(gridID3, "radians");
      gridDefYunits(gridID3, "radians");
    }
  else if ( xminval > -181 && xmaxval < 361 && yminval > -91 && ymaxval < 91 )
    {
      /* default is degrees */
    }
  else
    {
      cdoAbort("Units undefined!");
    }

  zaxisID3 = zaxisCreate(ZAXIS_SURFACE, 1);

  vlistID3 = vlistCreate();
  vlistDefVar(vlistID3, gridID3, zaxisID3, TSTEP_CONSTANT);
  vlistDefVarMissval(vlistID3, 0, missval);
  vlistDefVarName(vlistID3, 0, "dummy");
  vlistDefVarDatatype(vlistID3, 0, DATATYPE_INT8);

  taxisID3 = taxisCreate(TAXIS_ABSOLUTE);

  vlistDefTaxis(vlistID3, taxisID3);

  streamID3 = streamOpenWrite(cdoStreamName(2), cdoFiletype());

  streamDefVlist(streamID3, vlistID3);

  tsID = 0;
  streamDefTimestep(streamID3, tsID);

  for ( i = 0; i < gridsize; ++i ) array3[i] = missval;

  streamDefRecord(streamID3, 0, 0);
  streamWriteRecord(streamID3, array3, gridsize);

  streamClose(streamID3);

  if ( array1 ) free(array1);
  if ( array2 ) free(array2);
  if ( array3 ) free(array3);

  cdoFinish();

  return (0);
}
