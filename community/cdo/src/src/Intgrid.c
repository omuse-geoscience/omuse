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

      Intgrid    interpolate     PINGO grid interpolation
      Intgrid    intgridbil      Bilinear grid interpolation
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "interpol.h"


int genThinoutGrid(int gridID1, int xinc, int yinc)
{
  int ilon, ilat, olon, olat;
  int gridID2, gridtype;
  int nlon1, nlat1;
  int gridsize2, nlon2, nlat2;
  double *xvals1, *yvals1, *xvals2, *yvals2;

  gridtype = gridInqType(gridID1);
  nlon1 = gridInqXsize(gridID1);
  nlat1 = gridInqYsize(gridID1);

  nlon2 = nlon1/xinc;
  nlat2 = nlat1/yinc;
  if ( nlon1%xinc ) nlon2++;
  if ( nlat1%yinc ) nlat2++;
  gridsize2 = nlon2*nlat2;

  gridID2 = gridCreate(GRID_LONLAT, gridsize2);
  gridDefXsize(gridID2, nlon2);
  gridDefYsize(gridID2, nlat2);

  if ( gridtype == GRID_GAUSSIAN || gridtype == GRID_LONLAT )
    {
      xvals1 = (double*) malloc(nlon1*sizeof(double));
      yvals1 = (double*) malloc(nlat1*sizeof(double));
      xvals2 = (double*) malloc(nlon2*sizeof(double));
      yvals2 = (double*) malloc(nlat2*sizeof(double));
      gridInqXvals(gridID1, xvals1);
      gridInqYvals(gridID1, yvals1);

      olat = 0;
      for ( ilat = 0; ilat < nlat1; ilat+=yinc )
	{
	  yvals2[olat] = yvals1[ilat];
	  olat++;
	}

      olon = 0;
      for ( ilon = 0; ilon < nlon1; ilon+=xinc )
	{
	  xvals2[olon] = xvals1[ilon];
	  olon++;
	}

      gridDefXvals(gridID2, xvals2);
      gridDefYvals(gridID2, yvals2);
    }
  else
    {
      cdoAbort("Unsupported grid: %s", gridNamePtr(gridtype));
    }

  return gridID2;
}


int genBoxavgGrid(int gridID1, int xinc, int yinc)
{
  int i, j, i1;
  int gridID2, gridtype;
  int nlon1, nlat1;
  int gridsize2, nlon2, nlat2;
  double *xvals1, *yvals1, *xvals2, *yvals2;
  double *grid1_corner_lon = NULL, *grid1_corner_lat = NULL;
  double *grid2_corner_lon = NULL, *grid2_corner_lat = NULL;

  gridtype = gridInqType(gridID1);
  nlon1 = gridInqXsize(gridID1);
  nlat1 = gridInqYsize(gridID1);

  nlon2 = nlon1/xinc;
  nlat2 = nlat1/yinc;
  if ( nlon1%xinc ) nlon2++;
  if ( nlat1%yinc ) nlat2++;
  gridsize2 = nlon2*nlat2;

  gridID2 = gridCreate(GRID_LONLAT, gridsize2);
  gridDefXsize(gridID2, nlon2);
  gridDefYsize(gridID2, nlat2);

  if ( gridtype == GRID_GAUSSIAN || gridtype == GRID_LONLAT )
    {
      xvals1 = (double*) malloc(nlon1*sizeof(double));
      yvals1 = (double*) malloc(nlat1*sizeof(double));
      xvals2 = (double*) malloc(nlon2*sizeof(double));
      yvals2 = (double*) malloc(nlat2*sizeof(double));
      gridInqXvals(gridID1, xvals1);
      gridInqYvals(gridID1, yvals1);

      if ( gridInqYbounds(gridID1, NULL) && gridInqXbounds(gridID1, NULL) )
	{
	  grid1_corner_lon = (double*) malloc(2*nlon1*sizeof(double));
	  grid1_corner_lat = (double*) malloc(2*nlat1*sizeof(double));
	  grid2_corner_lon = (double*) malloc(2*nlon2*sizeof(double));
	  grid2_corner_lat = (double*) malloc(2*nlat2*sizeof(double));
	  gridInqXbounds(gridID1, grid1_corner_lon);
	  gridInqYbounds(gridID1, grid1_corner_lat);
	}

      j = 0;
      for ( i = 0; i < nlon1; i += xinc )
	{
	  i1 = i+(xinc-1);
	  if ( i1 >= nlon1-1 ) i1 = nlon1-1; 
	  xvals2[j] = xvals1[i] + (xvals1[i1] - xvals1[i])/2;
	  if ( grid2_corner_lon )
	    {
	      grid2_corner_lon[2*j] = grid1_corner_lon[2*i];
	      grid2_corner_lon[2*j+1] = grid1_corner_lon[2*i1+1];
	    }
	  j++;
	}
      j = 0;
      for ( i = 0; i < nlat1; i += yinc )
	{
	  i1 = i+(yinc-1);
	  if ( i1 >= nlat1-1 ) i1 = nlat1-1; 
	  yvals2[j] = yvals1[i] + (yvals1[i1] - yvals1[i])/2;
	  if ( grid2_corner_lat )
	    {
	      grid2_corner_lat[2*j] = grid1_corner_lat[2*i];
	      grid2_corner_lat[2*j+1] = grid1_corner_lat[2*i1+1];
	    }
	  j++;
	}

      gridDefXvals(gridID2, xvals2);
      gridDefYvals(gridID2, yvals2);

      if ( grid2_corner_lon && grid2_corner_lat )
	{
	  gridDefNvertex(gridID2, 2);
	  gridDefXbounds(gridID2, grid2_corner_lon);
	  gridDefYbounds(gridID2, grid2_corner_lat);

	  free(grid2_corner_lon);
	  free(grid2_corner_lat);
	}
    }
  else
    {
      cdoAbort("Unsupported grid: %s", gridNamePtr(gridtype));
    }

  return gridID2;
}


void boxavg(field_t *field1, field_t *field2, int xinc, int yinc)
{
  int nlon1, nlat1;
  int nlon2, nlat2;
  int ilat, ilon;
  int gridID1, gridID2;
  int nmiss;
  double **xfield1;
  double *array1, *array2;
  double missval;
  int i, j, ii, jj, in;
  double **xfield2;
  /* static int index = 0; */

  gridID1 = field1->grid;
  gridID2 = field2->grid;
  array1  = field1->ptr;
  array2  = field2->ptr;
  missval = field1->missval;

  nlon1 = gridInqXsize(gridID1);
  nlat1 = gridInqYsize(gridID1);

  nlon2 = gridInqXsize(gridID2);
  nlat2 = gridInqYsize(gridID2);

  xfield1 = (double **) malloc(nlat1*sizeof(double *));

  for ( ilat = 0; ilat < nlat1; ilat++ )
    xfield1[ilat] = array1 + ilat*nlon1;


  xfield2 = (double **) malloc(nlat2 * sizeof(double *));

  for ( ilat = 0; ilat < nlat2; ilat++ )
    xfield2[ilat] = array2 + ilat*nlon2;

  for ( ilat = 0; ilat < nlat2; ilat++ )
    for ( ilon = 0; ilon < nlon2; ilon++ )
      {
	xfield2[ilat][ilon] = 0;

	in = 0;
	for ( j = 0; j < yinc; ++j )
	  {
	    jj = ilat*yinc+j;
	    if ( jj >= nlat1 ) break;
	    for ( i = 0; i < xinc; ++i )
	      {
		ii = ilon*xinc+i;
		if ( ii >= nlon1 ) break;
		in++;
		xfield2[ilat][ilon] += xfield1[jj][ii];
	      }
	  }
	xfield2[ilat][ilon] /= in;
      }

  nmiss = 0;
  for ( i = 0; i < nlat2*nlon2; i++ )
    if ( DBL_IS_EQUAL(array2[i], missval) ) nmiss++;

  field2->nmiss = nmiss;

  free(xfield2);
  free(xfield1);
}


void thinout(field_t *field1, field_t *field2, int xinc, int yinc)
{
  int nlon1, nlat1;
  int nlon2, nlat2;
  int ilat, ilon, olat, olon;
  int gridID1, gridID2;
  int nmiss;
  double **xfield1;
  double *array1, *array2;
  double missval;
  int i;
  double **xfield2;

  gridID1 = field1->grid;
  gridID2 = field2->grid;
  array1  = field1->ptr;
  array2  = field2->ptr;
  missval = field1->missval;

  nlon1 = gridInqXsize(gridID1);
  nlat1 = gridInqYsize(gridID1);

  nlon2 = gridInqXsize(gridID2);
  nlat2 = gridInqYsize(gridID2);

  xfield1 = (double **) malloc(nlat1*sizeof(double *));

  for ( ilat = 0; ilat < nlat1; ilat++ )
    xfield1[ilat] = array1 + ilat*nlon1;

  xfield2 = (double **) malloc(nlat2*sizeof(double *));

  for ( ilat = 0; ilat < nlat2; ilat++ )
    xfield2[ilat] = array2 + ilat*nlon2;

  olat = 0;
  for ( ilat = 0; ilat < nlat1; ilat+=yinc )
    {
      olon = 0;
      for ( ilon = 0; ilon < nlon1; ilon+=xinc )
	{
	  xfield2[olat][olon] = xfield1[ilat][ilon];
	  olon++;
	}
      olat++;
    }

  nmiss = 0;
  for ( i = 0; i < nlat2*nlon2; i++ )
    if ( DBL_IS_EQUAL(array2[i], missval) ) nmiss++;
  
  field2->nmiss = nmiss;

  free(xfield2);
  free(xfield1);
}



void *Intgrid(void *argument)
{
  int INTGRIDBIL, INTGRIDCON, INTPOINT, INTERPOLATE, BOXAVG, THINOUT;
  int operatorID;
  int streamID1, streamID2;
  int nrecs, ngrids;
  int index;
  int tsID, recID, varID, levelID;
  int gridsize;
  int vlistID1, vlistID2;
  int gridID1 = -1, gridID2 = -1;
  int nmiss;
  int xinc = 0, yinc = 0;
  double missval;
  double slon, slat;
  double *array1 = NULL, *array2 = NULL;
  field_t field1, field2;
  int taxisID1, taxisID2;

  cdoInitialize(argument);

  INTGRIDBIL  = cdoOperatorAdd("intgridbil",  0, 0, NULL);
  INTGRIDCON  = cdoOperatorAdd("intgridcon",  0, 0, NULL);
  INTPOINT    = cdoOperatorAdd("intpoint",    0, 0, NULL);
  INTERPOLATE = cdoOperatorAdd("interpolate", 0, 0, NULL);
  BOXAVG      = cdoOperatorAdd("boxavg",      0, 0, NULL);
  THINOUT     = cdoOperatorAdd("thinout",     0, 0, NULL);

  operatorID = cdoOperatorID();

  if ( operatorID == INTGRIDBIL || operatorID == INTGRIDCON || operatorID == INTERPOLATE )
    {
      operatorInputArg("grid description file or name");
      gridID2 = cdoDefineGrid(operatorArgv()[0]);
    }
  else if ( operatorID == INTPOINT )
    {
      operatorInputArg("longitude and latitude");
      operatorCheckArgc(2);
      slon = parameter2double(operatorArgv()[0]);
      slat = parameter2double(operatorArgv()[1]);
      gridID2 = gridCreate(GRID_LONLAT, 1);
      gridDefXsize(gridID2, 1);
      gridDefYsize(gridID2, 1);
      gridDefXvals(gridID2, &slon);
      gridDefYvals(gridID2, &slat);
    }
  else if ( operatorID == THINOUT || operatorID == BOXAVG )
    {
      operatorInputArg("xinc, yinc");
      operatorCheckArgc(2);
      xinc = parameter2int(operatorArgv()[0]);
      yinc = parameter2int(operatorArgv()[1]);
    }

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  ngrids = vlistNgrids(vlistID1);
  for ( index = 0; index < ngrids; index++ )
    {
      gridID1 = vlistGrid(vlistID1, index);

      if ( operatorID == BOXAVG || operatorID == THINOUT )
	{
	  if ( index == 0 )
	    {
	      if ( gridInqType(gridID1) != GRID_LONLAT && gridInqType(gridID1) != GRID_GAUSSIAN 
		   /* && gridInqType(gridID1) != GRID_CURVILINEAR */ )
		cdoAbort("Interpolation of %s data unsupported!", gridNamePtr(gridInqType(gridID1)) );

	      if ( operatorID == BOXAVG )
		gridID2 = genBoxavgGrid(gridID1, xinc, yinc);
	      else
		gridID2 = genThinoutGrid(gridID1, xinc, yinc);
	    }
	  else
	    cdoAbort("Too many different grids!");
	}
      else
	{
	  if ( gridInqType(gridID1) != GRID_LONLAT && gridInqType(gridID1) != GRID_GAUSSIAN )
	    cdoAbort("Interpolation of %s data unsupported!", gridNamePtr(gridInqType(gridID1)) );

	  if ( gridIsRotated(gridID1) )
	    cdoAbort("Rotated grids not supported!");
	}

      vlistChangeGridIndex(vlistID2, index, gridID2);
    }

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  gridsize = vlistGridsizeMax(vlistID1);
  array1   = (double*) malloc(gridsize*sizeof(double));

  gridsize = gridInqSize(gridID2);
  array2   = (double*) malloc(gridsize*sizeof(double));

  field_init(&field1);
  field_init(&field2);

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);
	       
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamReadRecord(streamID1, array1, &nmiss);

	  gridID1 = vlistInqVarGrid(vlistID1, varID);
	  missval = vlistInqVarMissval(vlistID1, varID);

	  field1.grid    = gridID1;
	  field1.nmiss   = nmiss;
	  field1.missval = missval;
	  field1.ptr     = array1;
	  field2.grid    = gridID2;
	  field2.ptr     = array2;
	  field2.nmiss   = 0;

	  if ( operatorID == INTGRIDBIL || operatorID == INTPOINT )
	    intgridbil(&field1, &field2);
	  if ( operatorID == INTGRIDCON )
	    intgridcon(&field1, &field2);
	  else if ( operatorID == INTERPOLATE )
	    interpolate(&field1, &field2);
	  else if ( operatorID == BOXAVG )
	    boxavg(&field1, &field2, xinc, yinc);
	  else if ( operatorID == THINOUT )
	    thinout(&field1, &field2, xinc, yinc);

	  nmiss = field2.nmiss;

	  streamDefRecord(streamID2, varID, levelID);
	  streamWriteRecord(streamID2, array2, nmiss);
	}
      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  if ( array2 ) free(array2);
  if ( array1 ) free(array1);

  cdoFinish();

  return (0);
}
