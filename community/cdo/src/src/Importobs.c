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

#include <ctype.h>

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "grid.h"


#define  MAX_VARS   6


static
void init_vars(int vlistID, int gridID, int zaxisID, int nvars)
{
  int  code[]  = {11, 17, 33, 34, 1, 2/*, 3*/};
  char *name[]  = {"temp", "depoint", "u", "v", "height", "pressure" /*, "station"*/};
  char *units[] = {"Celsius", "", "m/s", "m/s", "m", "hPa" /*, ""*/};
  int varID;

  for ( int i = 0; i < nvars; ++i )
    {
      varID = vlistDefVar(vlistID, gridID, zaxisID, TSTEP_INSTANT);
      vlistDefVarParam(vlistID, varID, cdiEncodeParam(code[i], 255, 255));
      vlistDefVarName(vlistID, varID, name[i]);
      vlistDefVarUnits(vlistID, varID, units[i]);
      vlistDefVarDatatype(vlistID, varID, DATATYPE_FLT32);
    }
}

static
void init_data(int vlistID, int nvars, double *data[])
{
  int gridsize;
  double missval;

  for ( int varID = 0; varID < nvars; ++varID )
    {
      gridsize = gridInqSize(vlistInqVarGrid(vlistID, varID));
      missval  = vlistInqVarMissval(vlistID, varID);
      
      for ( int i = 0; i < gridsize; ++i ) data[varID][i] = missval;
    } 
}

static
void write_data(int streamID, int vlistID, int nvars, double *data[])
{
  int nmiss;
  int gridsize;
  double missval;

  for ( int varID = 0; varID < nvars; ++varID )
    {
      gridsize = gridInqSize(vlistInqVarGrid(vlistID, varID));
      missval  = vlistInqVarMissval(vlistID, varID);
      
      streamDefRecord(streamID, varID, 0);

      nmiss = 0;
      for ( int i = 0; i < gridsize; ++i )
	if ( DBL_IS_EQUAL(data[varID][i], missval) ) nmiss++;
      
      streamWriteRecord(streamID, data[varID], nmiss);
    }
}

static
int getDate(const char *name)
{
  char *pname = strchr(name, '_');

  int date = 0;
  if ( pname ) date = atoi(pname+1);

  return(date);
}

#define  MAX_LINE_LEN  4096

void *Importobs(void *argument)
{
  char line[MAX_LINE_LEN];
  int i, j;
  int nvars = MAX_VARS;
  int vtime = 0;
  double *data[MAX_VARS];
  char dummy[32], station[32], datetime[32];
  float lat, lon, height1, pressure, height2, value;
  double latmin = 90, latmax = -90, lonmin = 360, lonmax = -360;
  int code;
  int index;
  double dx, dy;
  char *pstation;

  cdoInitialize(argument);

  cdoOperatorAdd("import_obs",     0, 0, "grid description file or name");

  int operatorID = cdoOperatorID();

  operatorInputArg(cdoOperatorEnter(operatorID));  

  int gridID = cdoDefineGrid(operatorArgv()[0]);

  if ( gridInqType(gridID) != GRID_LONLAT ) 
    cdoAbort("Unsupported grid type: %s", gridNamePtr(gridInqType(gridID)));

  int gridsize = gridInqSize(gridID);
  int xsize = gridInqXsize(gridID);
  int ysize = gridInqYsize(gridID);

  // printf("gridsize=%d, xsize=%d, ysize=%d\n", gridsize, xsize, ysize);

  double *xvals = (double*) malloc(gridsize*sizeof(double));
  double *yvals = (double*) malloc(gridsize*sizeof(double));

  gridInqXvals(gridID, xvals);
  gridInqYvals(gridID, yvals);

  /* Convert lat/lon units if required */
  {
    char units[CDI_MAX_NAME];
    gridInqXunits(gridID, units);
    grid_to_degree(units, gridsize, xvals, "grid center lon");
    gridInqYunits(gridID, units);
    grid_to_degree(units, gridsize, yvals, "grid center lat");
  }

  FILE *fp = fopen(cdoStreamName(0)->args, "r");
  if ( fp == NULL ) { perror(cdoStreamName(0)->args); exit(EXIT_FAILURE); }

  int vdate = getDate(cdoStreamName(0)->args);
  if ( vdate <= 999999 ) vdate = vdate*100 + 1;

  int streamID = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  int zaxisID = zaxisCreate(ZAXIS_SURFACE, 1);

  int taxisID = taxisCreate(TAXIS_ABSOLUTE);

  int vlistID = vlistCreate();
  vlistDefTaxis(vlistID, taxisID);

  for ( i = 0; i < nvars; ++i ) data[i] = (double*) malloc(gridsize*sizeof(double));

  init_vars(vlistID, gridID, zaxisID, nvars);

  streamDefVlist(streamID, vlistID);

  int vdate0 = 0;
  int vtime0 = 0;
  //ntime = 0;
  int tsID = 0;
  while ( readline(fp, line, MAX_LINE_LEN) )
    {
      sscanf(line, "%s %s %s %g %g %g %d %g %g %g", 
             dummy, station, datetime, &lat, &lon, &height1, &code, &pressure, &height2, &value);
      sscanf(datetime, "%d_%d", &vdate, &vtime);

      if ( vdate != vdate0 || vtime != vtime0 )
        {
          if ( tsID > 0 ) write_data(streamID, vlistID, nvars, data);

          vdate0 = vdate;
          vtime0 = vtime;
          /*
            printf("%s %d %d %g %g %g %d %g %g %g\n", 
            station, vdate, vtime, lat, lon, height1, code, pressure, height2, value);	    
          */
          taxisDefVdate(taxisID, vdate);
          taxisDefVtime(taxisID, vtime);
          streamDefTimestep(streamID, tsID);
      
          init_data(vlistID, nvars, data);

          tsID++;
        }

      if ( lon < lonmin ) lonmin = lon;
      if ( lon > lonmax ) lonmax = lon;
      if ( lat < latmin ) latmin = lat;
      if ( lat > latmax ) latmax = lat;

      dy =  yvals[1] - yvals[0];
      for ( j = 0; j < ysize; ++j )
        if ( lat >= (yvals[j]-dy/2) && lat < (yvals[j]+dy/2) )  break;

      dx =  xvals[1] - xvals[0];
      if ( lon < (xvals[0] - dx/2) && lon < 0 ) lon+=360;
      for ( i = 0; i < xsize; ++i )
        if ( lon >= (xvals[i]-dx/2) && lon < (xvals[i]+dx/2) )  break;
	  
      index = -1;
      if ( code == 11 ) index = 0;
      if ( code == 17 ) index = 1;
      if ( code == 33 ) index = 2;
      if ( code == 34 ) index = 3;

      //printf("%d %d %d %g %g %g %g\n", i, j, index, dx, dy, lon, lat);
      if ( i < xsize && j < ysize && index >= 0 )
        {
          pstation = station;
          while (isalpha(*pstation)) pstation++;
          // printf("station %s %d\n", pstation, atoi(pstation));
          data[index][j*xsize+i] = value;
          data[    4][j*xsize+i] = height1;
          data[    5][j*xsize+i] = pressure;
          // data[    6][j*xsize+i] = atoi(pstation);
        }

      /*
        printf("%s %d %d %g %g %g %d %g %g %g\n", 
        station, vdate, vtime, lat, lon, height1, code, pressure, height2, value);
      */
    }

  write_data(streamID, vlistID, nvars, data);

  for ( i = 0; i < nvars; ++i ) free(data[i]);

  if ( cdoVerbose )
    printf("lonmin=%g, lonmax=%g, latmin=%g, latmax=%g\n", lonmin, lonmax, latmin, latmax);

  processDefVarNum(vlistNvars(vlistID), streamID);

  streamClose(streamID);

  fclose(fp);

  vlistDestroy(vlistID);
  gridDestroy(gridID);
  zaxisDestroy(zaxisID);
  taxisDestroy(taxisID);

  free(xvals);
  free(yvals);

  cdoFinish();

  return (0);
}
