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

      Gridcell   gridarea        Grid cell area in m^2
      Gridcell   gridweights     Grid cell weights
      Gridcell   gridmask        Grid mask
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "grid.h"
#include "constants.h"

static
double orthodrome(double px1, double py1, double px2, double py2)
{
  return acos(sin(py1)*sin(py2)+cos(py1)*cos(py2)*cos(px2-px1));
}


void *Gridcell(void *argument)
{
  int GRIDAREA, GRIDWGTS, GRIDMASK, GRIDDX, GRIDDY;
  int operatorID;
  int streamID1, streamID2;
  int vlistID1, vlistID2;
  int gridID, zaxisID;
  int gridtype;
  int status;
  int ngrids;
  int need_radius;
  int tsID, varID, levelID, taxisID;
  long i, gridsize;
  double *array = NULL;

  cdoInitialize(argument);

  GRIDAREA = cdoOperatorAdd("gridarea",     1,  0, NULL);
  GRIDWGTS = cdoOperatorAdd("gridweights",  1,  0, NULL);
  GRIDMASK = cdoOperatorAdd("gridmask",     0,  0, NULL);
  GRIDDX   = cdoOperatorAdd("griddx",       1,  0, NULL);
  GRIDDY   = cdoOperatorAdd("griddy",       1,  0, NULL);

  operatorID = cdoOperatorID();

  need_radius = cdoOperatorF1(operatorID);

  if ( need_radius )
    {
      char *envstr = getenv("PLANET_RADIUS");
      if ( envstr )
	{
	  double fval = atof(envstr);
	  if ( fval > 0 ) PlanetRadius = fval;
	}
    }

  if ( cdoVerbose ) cdoPrint("PlanetRadius: %g", PlanetRadius);

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);

  ngrids = vlistNgrids(vlistID1);

  if ( ngrids > 1 )
    cdoWarning("Found more than 1 grid, using the first one!");

  gridID  = vlistGrid(vlistID1, 0);

  zaxisID = zaxisCreate(ZAXIS_SURFACE, 1);

  vlistID2 = vlistCreate();
  varID    = vlistDefVar(vlistID2, gridID, zaxisID, TSTEP_CONSTANT);
  vlistDefNtsteps(vlistID2, 0);

  if ( operatorID == GRIDAREA )
    {
      vlistDefVarName(vlistID2, varID, "cell_area");
      vlistDefVarStdname(vlistID2, varID, "area");
      vlistDefVarLongname(vlistID2, varID, "area of grid cell");
      vlistDefVarUnits(vlistID2, varID, "m2");
      vlistDefVarDatatype(vlistID2, varID, DATATYPE_FLT64);
    }
  else if ( operatorID == GRIDWGTS )
    {
      vlistDefVarName(vlistID2, varID, "cell_weights");
      vlistDefVarDatatype(vlistID2, varID, DATATYPE_FLT64);
    }
  else if ( operatorID == GRIDMASK )
    {
      vlistDefVarName(vlistID2, varID, "grid_mask");
      vlistDefVarDatatype(vlistID2, varID, DATATYPE_UINT8);
    }
  else if ( operatorID == GRIDDX )
    {
      vlistDefVarName(vlistID2, varID, "dx");
      vlistDefVarLongname(vlistID2, varID, "delta x");
      vlistDefVarUnits(vlistID2, varID, "m");
    }
  else if ( operatorID == GRIDDY )
    {
      vlistDefVarName(vlistID2, varID, "dy");
      vlistDefVarLongname(vlistID2, varID, "delta y");
      vlistDefVarUnits(vlistID2, varID, "m");
    }

  taxisID = taxisCreate(TAXIS_ABSOLUTE);
  vlistDefTaxis(vlistID2, taxisID);


  gridsize = gridInqSize(gridID);
  array = (double*) malloc(gridsize*sizeof(double));


  if ( operatorID == GRIDAREA )
    {
      gridtype = gridInqType(gridID);
      if ( gridtype == GRID_LONLAT      ||
	   gridtype == GRID_GAUSSIAN    ||
	   gridtype == GRID_LCC         ||
	   gridtype == GRID_GME         ||
	   gridtype == GRID_CURVILINEAR ||
	   gridtype == GRID_UNSTRUCTURED )
	{
	  if ( gridHasArea(gridID) )
	    {
	      if ( cdoVerbose ) cdoPrint("Using existing grid cell area!");
	      gridInqArea(gridID, array);
	    }
	  else
	    {
	      status = gridGenArea(gridID, array);
	      if ( status == 1 )
		cdoAbort("Grid corner missing!");
	      else if ( status == 2 )
		cdoAbort("Can't compute grid cell areas for this grid!");

	      for ( i = 0; i < gridsize; ++i )
		array[i] *= PlanetRadius*PlanetRadius;
	    }
	}
      else
	{
	  if ( gridtype == GRID_GAUSSIAN_REDUCED )
	    cdoAbort("Unsupported grid type: %s, use CDO option -R to convert reduced to regular grid!",
		     gridNamePtr(gridtype));
	  else
	    cdoAbort("Unsupported grid type: %s", gridNamePtr(gridtype));
	}
    }
  else if ( operatorID == GRIDWGTS )
    {
      status = gridWeights(gridID, array);
      if ( status != 0 )
	  cdoWarning("Using constant grid cell area weights!");
    }
  else if ( operatorID == GRIDMASK )
    {
      int *mask;
      mask = (int*) malloc(gridsize*sizeof(int));
      if ( gridInqMask(gridID, NULL) )
	{
	  gridInqMask(gridID, mask);
	}
      else
	{
	  for ( i = 0; i < gridsize; ++i ) mask[i] = 1;
	}

      for ( i = 0; i < gridsize; ++i ) array[i] = mask[i];
      free(mask);
    }
  else if ( operatorID == GRIDDX || operatorID == GRIDDY )
    {
      gridtype = gridInqType(gridID);
      if ( gridtype == GRID_LONLAT      ||
	   gridtype == GRID_GAUSSIAN    ||
	   gridtype == GRID_LCC         ||
	   gridtype == GRID_CURVILINEAR )
	{
	  long i, j, xsize, ysize;
	  double *xv, *yv;
	  double len1 = 0, len2 = 0;
	  char units[CDI_MAX_NAME];

	  if ( gridtype != GRID_CURVILINEAR )
	    gridID = gridToCurvilinear(gridID, 1);

	  gridsize = gridInqSize(gridID);
	  xsize = gridInqXsize(gridID);
	  ysize = gridInqYsize(gridID);

	  xv = (double*) malloc(gridsize*sizeof(double));
	  yv = (double*) malloc(gridsize*sizeof(double));

	  gridInqXvals(gridID, xv);
	  gridInqYvals(gridID, yv);

	  /* Convert lat/lon units if required */

	  gridInqXunits(gridID, units);

	  grid_to_radian(units, gridsize, yv, "grid longitudes");
	  grid_to_radian(units, gridsize, yv, "grid latitudes");

	  if ( operatorID == GRIDDX )
	    {
	      for ( j = 0; j < ysize; ++j )
		for ( i = 0; i < xsize; ++i )
		  {
		    if ( i == 0 )
		      {
			len2 = orthodrome(xv[j*xsize+i], yv[j*xsize+i], xv[j*xsize+i+1], yv[j*xsize+i+1]);
			len1 = len2;
		      }
		    else if ( i == (xsize-1) )
		      {
			len1 = orthodrome(xv[j*xsize+i-1], yv[j*xsize+i-1], xv[j*xsize+i], yv[j*xsize+i]);
			len2 = len1;
		      }
		    else
		      {
			len1 = orthodrome(xv[j*xsize+i-1], yv[j*xsize+i-1], xv[j*xsize+i], yv[j*xsize+i]);
			len2 = orthodrome(xv[j*xsize+i], yv[j*xsize+i], xv[j*xsize+i+1], yv[j*xsize+i+1]);
		      }

		    array[j*xsize+i] = 0.5*(len1+len2)*PlanetRadius;
		  }
	    }
	  else
	    {
	      for ( i = 0; i < xsize; ++i )
	        for ( j = 0; j < ysize; ++j )
		  {
		    if ( j == 0 )
		      {
			len2 = orthodrome(xv[j*xsize+i], yv[j*xsize+i], xv[(j+1)*xsize+i], yv[(j+1)*xsize+i]);
			len1 = len2;
		      }
		    else if ( j == (ysize-1) )
		      {
			len1 = orthodrome(xv[(j-1)*xsize+i], yv[(j-1)*xsize+i], xv[j*xsize+i], yv[j*xsize+i]);
			len2 = len1;
		      }
		    else
		      {
			len1 = orthodrome(xv[(j-1)*xsize+i], yv[(j-1)*xsize+i], xv[j*xsize+i], yv[j*xsize+i]);
			len2 = orthodrome(xv[j*xsize+i], yv[j*xsize+i], xv[(j+1)*xsize+i], yv[(j+1)*xsize+i]);
		      }

		    array[j*xsize+i] = 0.5*(len1+len2)*PlanetRadius;
		  }
	    }

	  free(xv);
	  free(yv);
	}
      else
	{
	  if ( gridtype == GRID_GAUSSIAN_REDUCED )
	    cdoAbort("Unsupported grid type: %s, use CDO option -R to convert reduced to regular grid!",
		     gridNamePtr(gridtype));
	  else
	    cdoAbort("Unsupported grid type: %s", gridNamePtr(gridtype));
	}
    }


  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  tsID = 0;
  streamDefTimestep(streamID2, tsID);

  varID   = 0;
  levelID = 0;
  streamDefRecord(streamID2, varID, levelID);
  streamWriteRecord(streamID2, array, 0);

  streamClose(streamID2);
  streamClose(streamID1);

  if ( array ) free(array);

  cdoFinish();

  return (0);
}
