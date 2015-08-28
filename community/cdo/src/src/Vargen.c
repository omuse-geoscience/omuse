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

      Vargen     const           Create a constant field
      Vargen     random          Field with random values
      Vargen     stdatm          Field values for pressure and temperature for
                                 the standard atmosphere
*/


#if defined(HAVE_CONFIG_H)
#  include "config.h" // ENABLE_DATA
#endif

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "list.h"
#include "grid.h"
#include "constants.h"
#include "stdnametable.h"


#if defined(ENABLE_DATA)
  static double etopo_scale  = 3;
  static double etopo_offset = 11000;
  static const unsigned short etopo[] = {
#include "etopo.h"
  };

  static double temp_scale  =  500;
  static double temp_offset = -220;
  static const unsigned short temp[] = {
#include "temp.h"
  };

  static double mask_scale  =  1;
  static double mask_offset =  0;
  static const unsigned short mask[] = {
#include "mask.h"
  };
#endif

/*  some Constants for creating temperatur and pressure for the standard atmosphere */
#define T_ZERO          (213.0)
#define T_DELTA          (75.0)
#define SCALEHEIGHT   (10000.0)   /* [m] */
#define P_ZERO         (1013.25)  /* surface pressure [hPa] */
#define CC_R             (287.05)  /* specific gas constant for air */
static double TMP4PRESSURE = (C_EARTH_GRAV*SCALEHEIGHT)/(CC_R*T_ZERO);

static double
std_atm_temperatur(double height)
{
  /*
    Compute the temperatur for the given height (in meters) according to the
    solution of the hydrostatic atmosphere
   */
   return (T_ZERO + T_DELTA * exp((-1)*(height/SCALEHEIGHT)));
}

static double
std_atm_pressure(double height)
{
  /*
    Compute the pressure for the given height (in meters) according to the
    solution of the hydrostatic atmosphere
   */
  return (P_ZERO * exp((-1)*TMP4PRESSURE*log((exp(height/SCALEHEIGHT)*T_ZERO + T_DELTA)/(T_ZERO + T_DELTA))));
}

int rect_grid_search(long *ii, long *jj, double x, double y, long nxm, long nym, const double *restrict xm, const double *restrict ym);

void remap_nn_reg2d(int nx, int ny, double *restrict data, int gridID, double *array)
{
  int gridsize = gridInqSize(gridID);
  double *xvals = (double*) malloc(gridsize*sizeof(double));
  double *yvals = (double*) malloc(gridsize*sizeof(double));

  if ( gridInqType(gridID) == GRID_GME ) gridID = gridToUnstructured(gridID, 0);

  if ( gridInqType(gridID) != GRID_UNSTRUCTURED && gridInqType(gridID) != GRID_CURVILINEAR )
    gridID = gridToCurvilinear(gridID, 0);

  gridInqXvals(gridID, xvals);
  gridInqYvals(gridID, yvals);

  /* Convert lat/lon units if required */
  char units[CDI_MAX_NAME];
  gridInqXunits(gridID, units);
  grid_to_degree(units, gridsize, xvals, "grid center lon");
  gridInqYunits(gridID, units);
  grid_to_degree(units, gridsize, yvals, "grid center lat");

  int ii, jj;
  double xval, yval;
  for ( int i = 0; i < gridsize; i++ )
    {
      xval = xvals[i];
      yval = yvals[i];
      if ( xval >=  180 ) xval -= 360;
      if ( xval <  -180 ) xval += 360;
      ii = (xval + 180)*2;
      jj = (yval +  90)*2;
      if ( ii > nx ) ii = nx;
      if ( jj > ny ) jj = ny;
      array[i] = data[jj*nx+ii];
    }

  free(xvals);
  free(yvals);
}


void *Vargen(void *argument)
{
  int ntimesteps, nlevels = 1;
  int varID, varID2 = -1, levelID;
  int i;
  int gridID = -1, gridIDdata = -1, zaxisID;
  int vdate, vtime;
  double rval, rstart = 0, rstop = 0, rinc = 0;
  double rconst = 0;
  double *levels = NULL;
  double lon[720], lat[360];
  int nlon = 720;
  int nlat = 360;

  cdoInitialize(argument);

  int RANDOM  = cdoOperatorAdd("random",  0, 0, "grid description file or name, <seed>");
  int SINCOS  = cdoOperatorAdd("sincos",  0, 0, "grid description file or name");
  int COSHILL = cdoOperatorAdd("coshill", 0, 0, "grid description file or name");
  int CONST   = cdoOperatorAdd("const",   0, 0, "constant value, grid description file or name");
  int FOR     = cdoOperatorAdd("for",     0, 0, "start, end, <increment>");
  int TOPO    = cdoOperatorAdd("topo",    0, 0, NULL);
  int TEMP    = cdoOperatorAdd("temp",    0, 0, NULL);
  int MASK    = cdoOperatorAdd("mask",    0, 0, NULL);
  int STDATM  = cdoOperatorAdd("stdatm",  0, 0, "height levels [m]");

  int operatorID = cdoOperatorID();

  if ( operatorID == RANDOM )
    {
      unsigned int seed = 1;
      operatorInputArg(cdoOperatorEnter(operatorID));
      if ( operatorArgc() < 1 ) cdoAbort("Too few arguments!");
      if ( operatorArgc() > 2 ) cdoAbort("Too many arguments!");
      gridID = cdoDefineGrid(operatorArgv()[0]);
      if ( operatorArgc() == 2 )
        {
          int idum = parameter2int(operatorArgv()[1]);
          if ( idum >= 0 && idum < 0x7FFFFFFF ) seed = idum;
        }
      srand(seed);
    }
  else if ( operatorID == SINCOS || operatorID == COSHILL )
    {
      operatorInputArg(cdoOperatorEnter(operatorID));
      operatorCheckArgc(1);
      gridID = cdoDefineGrid(operatorArgv()[0]);
    }
  else if ( operatorID == CONST )
    {
      operatorInputArg(cdoOperatorEnter(operatorID));
      operatorCheckArgc(2);
      rconst = parameter2double(operatorArgv()[0]);
      gridID = cdoDefineGrid(operatorArgv()[1]);
    }
  else if ( operatorID == TOPO || operatorID == TEMP || operatorID == MASK )
    {
      gridIDdata = gridCreate(GRID_LONLAT, nlon*nlat);
      gridDefXsize(gridIDdata, nlon);
      gridDefYsize(gridIDdata, nlat);

      for ( int i = 0; i < nlon; i++ ) lon[i] = -179.75 + i*0.5;
      for ( int i = 0; i < nlat; i++ ) lat[i] = -89.75 + i*0.5;

      gridDefXvals(gridIDdata, lon);
      gridDefYvals(gridIDdata, lat);

      gridID = gridIDdata;

      if ( operatorArgc() == 1 ) gridID = cdoDefineGrid(operatorArgv()[0]);
    }
  else if ( operatorID == FOR )
    {
      double lon = 0, lat = 0;
      operatorInputArg(cdoOperatorEnter(operatorID));
      if ( operatorArgc() < 2 ) cdoAbort("Too few arguments!");
      if ( operatorArgc() > 3 ) cdoAbort("Too many arguments!");

      rstart = parameter2double(operatorArgv()[0]);
      rstop  = parameter2double(operatorArgv()[1]);
      if ( operatorArgc() == 3 )
        rinc = parameter2double(operatorArgv()[2]);
      else
        rinc = 1;

      if ( DBL_IS_EQUAL(rinc, 0.0) ) cdoAbort("Increment is zero!");

      gridID = gridCreate(GRID_LONLAT, 1);
      gridDefXsize(gridID, 1);
      gridDefYsize(gridID, 1);
      gridDefXvals(gridID, &lon);
      gridDefYvals(gridID, &lat);
    }
  else if ( operatorID == STDATM )
    {
      double lon = 0, lat = 0;
      LIST *flist = listNew(FLT_LIST);

      operatorInputArg(cdoOperatorEnter(operatorID));
      nlevels = args2fltlist(operatorArgc(), operatorArgv(), flist);
      levels  = (double *) listArrayPtr(flist);
      //listDelete(flist);

      if ( cdoVerbose ) for ( i = 0; i < nlevels; ++i ) printf("levels %d: %g\n", i, levels[i]);

      gridID = gridCreate(GRID_LONLAT, 1);
      gridDefXsize(gridID, 1);
      gridDefYsize(gridID, 1);
      gridDefXvals(gridID, &lon);
      gridDefYvals(gridID, &lat);
    }

  if ( operatorID == STDATM )
    {
      zaxisID = zaxisCreate(ZAXIS_HEIGHT, nlevels);
      zaxisDefLevels(zaxisID  , levels);
      zaxisDefName(zaxisID    , "level");
      zaxisDefLongname(zaxisID, "Level");
      zaxisDefUnits(zaxisID   , "m");
    }
  else
    {
      zaxisID = zaxisCreate(ZAXIS_SURFACE, 1);
      nlevels = 1;
    }

  int vlistID = vlistCreate();

  int tsteptype = TSTEP_CONSTANT;
  if ( operatorID == FOR ) tsteptype = TSTEP_INSTANT;

  varID = vlistDefVar(vlistID, gridID, zaxisID, tsteptype);
  /*
     For the standard atmosphere two output variables are generated: pressure and
     temperatur. The first (varID) is pressure, second (varID2) is temperatur.
     Add an additional variable for the standard atmosphere.
   */
  if ( operatorID == STDATM )
    varID2 = vlistDefVar(vlistID, gridID, zaxisID, TSTEP_CONSTANT);

  if ( operatorID == MASK )
    vlistDefVarDatatype(vlistID, varID, DATATYPE_INT8);

  if ( operatorID == STDATM )
    {
      vlistDefVarName(vlistID    , varID , "P");
      vlistDefVarParam(vlistID   , varID , cdiEncodeParam(1, 255, 255));
      vlistDefVarStdname(vlistID , varID , "air_pressure");
      vlistDefVarLongname(vlistID, varID , "pressure");
      vlistDefVarUnits(vlistID   , varID , "hPa");

      vlistDefVarName(vlistID    , varID2, "T");
      vlistDefVarParam(vlistID   , varID2, cdiEncodeParam(130, 128, 255));
      vlistDefVarStdname(vlistID , varID2, var_stdname(air_temperature));
      vlistDefVarLongname(vlistID, varID2, "temperature");
      vlistDefVarUnits(vlistID   , varID2, "K");
    }
  else
    {
      vlistDefVarName(vlistID, varID, cdoOperatorName(operatorID));
      if ( operatorID == TOPO ) vlistDefVarUnits(vlistID, varID , "m");	
      if ( operatorID == TEMP ) vlistDefVarUnits(vlistID, varID , "K");	
    }

  int taxisID = taxisCreate(TAXIS_RELATIVE);
  vlistDefTaxis(vlistID, taxisID);

  if ( operatorID == RANDOM || operatorID == SINCOS || operatorID == COSHILL || operatorID == CONST ||
       operatorID == TOPO || operatorID == TEMP || operatorID == MASK || operatorID == STDATM )
    vlistDefNtsteps(vlistID, 1);

  int streamID = streamOpenWrite(cdoStreamName(0), cdoFiletype());

  streamDefVlist(streamID, vlistID);

  int gridsize = gridInqSize(gridID);
  int datasize = gridsize;
  double *array = (double*) malloc(gridsize*sizeof(double));
  double *data = array;
  if ( gridID != gridIDdata && gridIDdata != -1 )
    {
      datasize = gridInqSize(gridIDdata);
      data = (double*) malloc(datasize*sizeof(double));
    }
  
  if ( operatorID == FOR )
    ntimesteps = 1.001 + ((rstop-rstart)/rinc);
  else
    {
      vlistDefNtsteps(vlistID, 0);
      ntimesteps = 1;
    }

  int julday = date_to_julday(CALENDAR_PROLEPTIC, 10101);

  int nvars = vlistNvars(vlistID);

  for ( int tsID = 0; tsID < ntimesteps; tsID++ )
    {
      rval  = rstart + rinc*tsID;
      vdate = julday_to_date(CALENDAR_PROLEPTIC, julday + tsID);
      vtime = 0;
      taxisDefVdate(taxisID, vdate);
      taxisDefVtime(taxisID, vtime);
      streamDefTimestep(streamID, tsID);

      for ( varID = 0; varID < nvars; varID++ )
        {
          nlevels = zaxisInqSize(vlistInqVarZaxis(vlistID, varID));
          for ( levelID = 0; levelID < nlevels; levelID++ )
            {
              streamDefRecord(streamID, varID, levelID);

              if ( operatorID == RANDOM )
                {
                  for ( i = 0; i < gridsize; i++ )
                    array[i] = rand()/(RAND_MAX+1.0);
                }
              else if ( operatorID == SINCOS || operatorID == COSHILL )
                {
		  double *xvals = (double*) malloc(gridsize*sizeof(double));
		  double *yvals = (double*) malloc(gridsize*sizeof(double));

		  if ( gridInqType(gridID) == GRID_GME ) gridID = gridToUnstructured(gridID, 0);

		  if ( gridInqType(gridID) != GRID_UNSTRUCTURED && gridInqType(gridID) != GRID_CURVILINEAR )
		    gridID = gridToCurvilinear(gridID, 0);

		  gridInqXvals(gridID, xvals);
		  gridInqYvals(gridID, yvals);

		  /* Convert lat/lon units if required */
		  char units[CDI_MAX_NAME];
		  gridInqXunits(gridID, units);
		  grid_to_radian(units, gridsize, xvals, "grid center lon");
		  gridInqYunits(gridID, units);
		  grid_to_radian(units, gridsize, yvals, "grid center lat");

		  if ( operatorID == SINCOS )
		    {
		      for ( i = 0; i < gridsize; i++ )
			array[i] = cos(1.0 * xvals[i]) * sin(2.0 * yvals[i]);
		    }
		  else if ( operatorID == COSHILL )
		    {		     
		      for ( i = 0; i < gridsize; i++ )
			array[i] = 2 - cos(acos(cos(xvals[i]) * cos(yvals[i]))/1.2);
		    }

		  free(xvals);
		  free(yvals);
		}
              else if ( operatorID == CONST )
                {
                  for ( i = 0; i < gridsize; i++ )
                    array[i] = rconst;
                }
              else if ( operatorID == TOPO )
                {
#if defined(ENABLE_DATA)
                  for ( i = 0; i < datasize; i++ )
                    data[i] = etopo[i]/etopo_scale - etopo_offset;
#else
                  cdoAbort("Operator support disabled!");
#endif
                }
              else if ( operatorID == TEMP )
                {
#if defined(ENABLE_DATA)
                  for ( i = 0; i < datasize; i++ )
                    data[i] = temp[i]/temp_scale - temp_offset;
#else
                  cdoAbort("Operator support disabled!");
#endif
                }
              else if ( operatorID == MASK )
                {
#if defined(ENABLE_DATA)
                  for ( i = 0; i < datasize; i++ )
                    data[i] = mask[i]/mask_scale - mask_offset;
#else
                  cdoAbort("Operator support disabled!");
#endif
                }
              else if ( operatorID == FOR )
                {
                  array[0] = rval;
                }
              else if ( operatorID == STDATM )
                {
                  array[0] = (varID == varID2) ? std_atm_temperatur(levels[levelID]) : std_atm_pressure(levels[levelID]);
                }

              if ( gridID != gridIDdata && (operatorID == TOPO || operatorID == TEMP || operatorID == MASK) )
                {
                  remap_nn_reg2d(nlon, nlat, data, gridID, array);
                }

              streamWriteRecord(streamID, array, 0);
            }
        }
    }

  streamClose(streamID);

  vlistDestroy(vlistID);

  if ( gridID != gridIDdata && gridIDdata != -1 ) free(data);
  if ( array ) free(array);
  if ( levels ) free(levels); 

  cdoFinish();

  return (0);
}
