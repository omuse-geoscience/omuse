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

      Invert     invertlat       Invert latitude
      Invert     invertlon       Invert longitude
      Invert     invertlatdes    Invert latitude description
      Invert     invertlondes    Invert longitude description
      Invert     invertlatdata   Invert latitude data
      Invert     invertlondata   Invert longitude data
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "error.h"


static
void invertLonDes(int vlistID)
{
  int index, ngrids;
  int gridID1, gridID2;
  int nlat, nlon, size;
  int ilat, ilon;
  int gridtype, nv, iv;
  double *xv1, *xv2;
  double *xb1, *xb2;

  ngrids = vlistNgrids(vlistID);
  for ( index = 0; index < ngrids; index++ )
    {
      gridID1 = vlistGrid(vlistID, index);
      gridID2 = gridDuplicate(gridID1);

      gridtype = gridInqType(gridID1);

      if ( gridtype != GRID_GENERIC && gridtype != GRID_GAUSSIAN &&
	   gridtype != GRID_LONLAT  && gridtype != GRID_CURVILINEAR )
	cdoAbort("Unsupported gridtype!");

      if ( gridInqXvals(gridID1, NULL) )
	{
	  nlon  = gridInqXsize(gridID1);
	  nlat  = gridInqYsize(gridID1);

	  if ( gridtype == GRID_CURVILINEAR )
	    size = nlon*nlat;
	  else
            size = nlon;

	  xv1 = (double*) malloc(size*sizeof(double));
	  xv2 = (double*) malloc(size*sizeof(double));

	  gridInqXvals(gridID1, xv1);

	  if ( gridtype == GRID_CURVILINEAR )
	    {
	      for ( ilat = 0; ilat < nlat; ilat++ )
		for ( ilon = 0; ilon < nlon; ilon++ )
		  xv2[ilat*nlon + nlon-ilon-1] = xv1[ilat*nlon + ilon];
	    }
	  else
	    {
	      for ( ilon = 0; ilon < nlon; ilon++ )
		xv2[nlon-ilon-1] = xv1[ilon];
	    }

	  gridDefXvals(gridID2, xv2);

	  if ( xv2 ) free(xv2);
	  if ( xv1 ) free(xv1);
	}

      if ( gridInqXbounds(gridID1, NULL) )
	{
	  nlon  = gridInqXsize(gridID1);
	  nlat  = gridInqYsize(gridID1);

	  nv = gridInqNvertex(gridID1);

	  if ( gridtype == GRID_CURVILINEAR )
	    size = nv*nlon*nlat;
	  else
            size = nv*nlon;

	  xb1 = (double*) malloc(size*sizeof(double));
	  xb2 = (double*) malloc(size*sizeof(double));

	  gridInqXbounds(gridID1, xb1);

	  if ( gridtype == GRID_CURVILINEAR )
	    {
	      for ( ilat = 0; ilat < nlat; ilat++ )
		for ( ilon = 0; ilon < nlon; ilon++ )
		  for ( iv = 0; iv < nv; iv++ )
		    xb2[ilat*nlon*nv + (nlon-ilon-1)*nv + iv] = xb1[ilat*nlon*nv + ilon*nv + iv];
	    }
	  else
	    {
		for ( ilon = 0; ilon < nlon; ilon++ )
		  {
		    xb2[nlon*2-ilon*2-1] = xb1[ilon*2];
		    xb2[nlon*2-ilon*2-2] = xb1[ilon*2+1];
		  }
	    }

	  gridDefXbounds(gridID2, xb2);

	  if ( xb2 ) free(xb2);
	  if ( xb1 ) free(xb1);
	}

      vlistChangeGrid(vlistID, gridID1, gridID2);
    }
}

static
void invertLatDes(int vlistID)
{
  int index, ngrids;
  int gridID1, gridID2;
  int nlat, nlon, size;
  int ilat, ilon;
  int gridtype, nv, iv;
  double *yv1, *yv2;
  double *yb1, *yb2;

  ngrids = vlistNgrids(vlistID);
  for ( index = 0; index < ngrids; index++ )
    {
      gridID1 = vlistGrid(vlistID, index);
      gridID2 = gridDuplicate(gridID1);

      gridtype = gridInqType(gridID1);

      if ( gridtype != GRID_GENERIC && gridtype != GRID_GAUSSIAN &&
	   gridtype != GRID_LONLAT  && gridtype != GRID_CURVILINEAR )
	cdoAbort("Unsupported gridtype!");

      if ( gridInqYvals(gridID1, NULL) )
	{
	  nlon  = gridInqXsize(gridID1);
	  nlat  = gridInqYsize(gridID1);

	  if ( gridtype == GRID_CURVILINEAR )
	    size = nlon*nlat;
	  else
            size = nlat;

	  yv1 = (double*) malloc(size*sizeof(double));
	  yv2 = (double*) malloc(size*sizeof(double));


	  if ( gridtype == GRID_CURVILINEAR )
	    {
	      gridInqXvals(gridID1, yv1);

	      for ( ilat = 0; ilat < nlat; ilat++ )
		for ( ilon = 0; ilon < nlon; ilon++ )
		  yv2[(nlat-ilat-1)*nlon + ilon] = yv1[ilat*nlon + ilon];

	      gridDefXvals(gridID2, yv2);

	      gridInqYvals(gridID1, yv1);

	      for ( ilat = 0; ilat < nlat; ilat++ )
		for ( ilon = 0; ilon < nlon; ilon++ )
		  yv2[(nlat-ilat-1)*nlon + ilon] = yv1[ilat*nlon + ilon];

	      gridDefYvals(gridID2, yv2);
	    }
	  else
	    {
	      gridInqYvals(gridID1, yv1);

	      for ( ilat = 0; ilat < nlat; ilat++ )
		yv2[nlat-ilat-1] = yv1[ilat];

	      gridDefYvals(gridID2, yv2);
	    }

	  if ( yv2 ) free(yv2);
	  if ( yv1 ) free(yv1);
	}

      if ( gridInqYbounds(gridID1, NULL) )
	{
	  nlon  = gridInqXsize(gridID1);
	  nlat  = gridInqYsize(gridID1);

	  nv = gridInqNvertex(gridID1);

	  if ( gridtype == GRID_CURVILINEAR )
	    size = nv*nlon*nlat;
	  else
            size = nv*nlat;

	  yb1 = (double*) malloc(size*sizeof(double));
	  yb2 = (double*) malloc(size*sizeof(double));

	  gridInqYbounds(gridID1, yb1);

	  if ( gridtype == GRID_CURVILINEAR )
	    {
	      for ( ilat = 0; ilat < nlat; ilat++ )
		for ( ilon = 0; ilon < nlon; ilon++ )
		  for ( iv = 0; iv < nv; iv++ )
		    yb2[(nlat-ilat-1)*nlon*nv + ilon*nv + iv] = yb1[ilat*nlon*nv + ilon*nv + iv];
	    }
	  else
	    {
		for ( ilat = 0; ilat < nlat; ilat++ )
		  {
		    yb2[nlat*2-ilat*2-1] = yb1[ilat*2];
		    yb2[nlat*2-ilat*2-2] = yb1[ilat*2+1];
		  }
	    }

	  gridDefYbounds(gridID2, yb2);

	  if ( yb2 ) free(yb2);
	  if ( yb1 ) free(yb1);
	}

      vlistChangeGrid(vlistID, gridID1, gridID2);
    }
}

static
void invertLonData(double *array1, double *array2, int gridID1)
{
  int nlat, nlon;
  int ilat, ilon;
  double **field1, **field2;

  nlon = gridInqXsize(gridID1);
  nlat = gridInqYsize(gridID1);

  if ( nlat > 0 )
    {
      field1 = (double **) malloc(nlat*sizeof(double *));
      field2 = (double **) malloc(nlat*sizeof(double *));
  
      for ( ilat = 0; ilat < nlat; ilat++ )
	{
	  field1[ilat] = array1 + ilat*nlon;
	  field2[ilat] = array2 + ilat*nlon;
	}

      for ( ilat = 0; ilat < nlat; ilat++ )
	for ( ilon = 0; ilon < nlon; ilon++ )
	  field2[ilat][nlon-ilon-1] = field1[ilat][ilon];
  
      if ( field1 ) free(field1);
      if ( field2 ) free(field2);
    }
  else
    {
      array2[0] = array1[0];
    }
}

static
void invertLatData(double *array1, double *array2, int gridID1)
{
  int nlat, nlon;
  int ilat;
  double **field1, **field2;

  nlon = gridInqXsize(gridID1);
  nlat = gridInqYsize(gridID1);

  if ( nlat > 0 )
    {
      field1 = (double **) malloc(nlat*sizeof(double *));
      field2 = (double **) malloc(nlat*sizeof(double *));
  
      for ( ilat = 0; ilat < nlat; ilat++ )
	{
	  field1[ilat] = array1 + ilat*nlon;
	  field2[ilat] = array2 + ilat*nlon;
	}

      for ( ilat = 0; ilat < nlat; ilat++ )
	memcpy(field2[nlat-ilat-1], field1[ilat], nlon*sizeof(double));
      
      if ( field1 ) free(field1);
      if ( field2 ) free(field2);
    }
  else
    {
      array2[0] = array1[0];
    }
}


void *Invert(void *argument)
{
  int nrecs;
  int recID, varID, levelID;
  int gridID1;
  int nmiss;

  cdoInitialize(argument);

  cdoOperatorAdd("invertlat",     func_all, func_lat, NULL);
  cdoOperatorAdd("invertlon",     func_all, func_lon, NULL);
  cdoOperatorAdd("invertlatdes",  func_hrd, func_lat, NULL);
  cdoOperatorAdd("invertlondes",  func_hrd, func_lon, NULL);
  cdoOperatorAdd("invertlatdata", func_fld, func_lat, NULL);
  cdoOperatorAdd("invertlondata", func_fld, func_lon, NULL);

  int operatorID = cdoOperatorID();
  int operfunc1 = cdoOperatorF1(operatorID);
  int operfunc2 = cdoOperatorF2(operatorID);

  int streamID1 = streamOpenRead(cdoStreamName(0));

  int vlistID1 = streamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  if ( operfunc1 == func_all || operfunc1 == func_hrd )
    {
      if ( operfunc2 == func_lat )
	invertLatDes(vlistID2);
      else
	invertLonDes(vlistID2);
    }

  int streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  int gridsize = vlistGridsizeMax(vlistID1);

  double *array1 = (double*) malloc(gridsize*sizeof(double));
  double *array2 = (double*) malloc(gridsize*sizeof(double));

  int tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);
	       
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamReadRecord(streamID1, array1, &nmiss);

	  streamDefRecord(streamID2, varID, levelID);

	  if ( operfunc1 == func_all || operfunc1 == func_fld )
	    {
	      gridID1 = vlistInqVarGrid(vlistID1, varID);

	      if ( operfunc2 == func_lat )
		invertLatData(array1, array2, gridID1);
	      else
		invertLonData(array1, array2, gridID1);

	      streamWriteRecord(streamID2, array2, nmiss);     
	    }
	  else
	    {
	      streamWriteRecord(streamID2, array1, nmiss);     
	    }
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
