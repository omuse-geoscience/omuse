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

      Output     output          ASCII output
      Output     outputf         Formatted output
      Output     outputint       Integer output
      Output     outputsrv       SERVICE output
      Output     outputext       EXTRA output
      Output     outputtab       Table output
*/

#include <ctype.h>
#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "grid.h"
#include "pstream.h"


void *Output(void *argument)
{
  int OUTPUT, OUTPUTINT, OUTPUTSRV, OUTPUTEXT, OUTPUTF, OUTPUTTS, OUTPUTFLD, OUTPUTARR, OUTPUTXYZ, OUTPUTTAB;
  int operatorID;
  int i;
  int indf;
  int varID, recID;
  int gridsize = 0;
  int gridID, zaxisID, code, vdate, vtime;
  int param;
  int gridtype;
  int ngrids;
  int nrecs;
  int levelID;
  int tsID, taxisID;
  int streamID = 0;
  int vlistID;
  int nmiss, nout;
  int nlon, nlat;
  int nelem = 1;
  int len;
  int index;
  int ndiffgrids;
  int lhead = TRUE;
  const char *format = NULL;
  char paramstr[32];
  char vdatestr[32], vtimestr[32];
  double level;
  double *grid_center_lon = NULL, *grid_center_lat = NULL;
  double *array = NULL;
  double xdate;
  double missval;
  double lon, lat;
  char name[CDI_MAX_NAME];
  int year, month, day;
  int *keys = NULL, nkeys = 0, k;
  int nKeys;
  int Keylen[]           = {      0,        8,      11,      4,      8,     6,     6,     6,     6,      4,      4,          6,     10,      8,      5,       2,     2 };
  enum                     {knohead,   kvalue,  kparam,  kcode,  kname,  klon,  klat,  klev,  kbin,  kxind,  kyind,  ktimestep,  kdate,  ktime,  kyear,  kmonth,  kday };
  const char *Keynames[] = {"nohead", "value", "param", "code", "name", "lon", "lat", "lev", "bin", "xind", "yind", "timestep", "date", "time", "year", "month", "day"};


  cdoInitialize(argument);

  OUTPUT    = cdoOperatorAdd("output",    0, 0, NULL);
  OUTPUTINT = cdoOperatorAdd("outputint", 0, 0, NULL);
  OUTPUTSRV = cdoOperatorAdd("outputsrv", 0, 0, NULL);
  OUTPUTEXT = cdoOperatorAdd("outputext", 0, 0, NULL);
  OUTPUTF   = cdoOperatorAdd("outputf",   0, 0, NULL);
  OUTPUTTS  = cdoOperatorAdd("outputts",  0, 0, NULL);
  OUTPUTFLD = cdoOperatorAdd("outputfld", 0, 0, NULL);
  OUTPUTARR = cdoOperatorAdd("outputarr", 0, 0, NULL);
  OUTPUTXYZ = cdoOperatorAdd("outputxyz", 0, 0, NULL);
  OUTPUTTAB = cdoOperatorAdd("outputtab", 0, 0, NULL);

  UNUSED(OUTPUT);

  operatorID = cdoOperatorID();

  if ( operatorID == OUTPUTF )
    {
      operatorInputArg("format and number of elements [optional]");

      if ( operatorArgc() < 1 ) cdoAbort("Too few arguments!");

      format = operatorArgv()[0];
      if ( operatorArgc() == 2 ) nelem = parameter2int(operatorArgv()[1]);
    }
  else if ( operatorID == OUTPUTTAB )
    {
      operatorInputArg("keys to print");
 
      int npar = operatorArgc();
      char **parnames = operatorArgv();

      if ( cdoVerbose )
	for ( i = 0; i < npar; i++ )
	  cdoPrint("key %d = %s", i+1, parnames[i]);

      keys = (int*) malloc(npar*sizeof(int));
      nkeys = 0;
      nKeys = sizeof(Keynames)/sizeof(char *);
      for ( i = 0; i < npar; i++ )
	{
	  for ( k = 0; k < nKeys; ++k )
	    {
	      //	      len = strlen(parnames[i]);
	      len = strlen(Keynames[k]);
	      if ( len < 3 ) len = 3;
	      if ( strncmp(parnames[i], Keynames[k], len) == 0 )
		{
		  if ( k == knohead ) lhead = FALSE;
		  else
		    {
		      keys[nkeys++] = k;
		      if ( parnames[i][len] == ':' && isdigit(parnames[i][len+1]) ) Keylen[k] = atoi(&parnames[i][len+1]);
		    }
		  break;
		}
	    }

	  if ( k == nKeys ) cdoAbort("Key %s unsupported!", parnames[i]);
	}
 
      if ( cdoVerbose )
	for ( k = 0; k < nkeys; ++k )
	  cdoPrint("keynr = %d  keyid = %d  keylen = %d  keyname = %s", k, keys[k], Keylen[keys[k]], Keynames[keys[k]]);

      if ( lhead )
	{
	  fprintf(stdout, "#");
	  for ( k = 0; k < nkeys; ++k )
	    {
	      len = Keylen[keys[k]];
	      //   if ( k == 0 ) len -= 1;
	      fprintf(stdout, "%*s ", len, Keynames[keys[k]]);
	    }
	  fprintf(stdout, "\n");
	}
    }

  for ( indf = 0; indf < cdoStreamCnt(); indf++ )
    {
      streamID = streamOpenRead(cdoStreamName(indf));

      vlistID = streamInqVlist(streamID);

      ngrids = vlistNgrids(vlistID);
      ndiffgrids = 0;
      for ( index = 1; index < ngrids; index++ )
	if ( vlistGrid(vlistID, 0) != vlistGrid(vlistID, index) )
	  ndiffgrids++;

      if ( ndiffgrids > 0 ) cdoAbort("Too many different grids!");

      gridID   = vlistGrid(vlistID, 0);
      gridsize = gridInqSize(gridID);
      gridtype = gridInqType(gridID);

      array = (double*) malloc(gridsize*sizeof(double));

      if ( operatorID == OUTPUTFLD || operatorID == OUTPUTXYZ || operatorID == OUTPUTTAB )
	{
	  if ( gridInqType(gridID) == GRID_GME ) gridID = gridToUnstructured(gridID, 0);

	  if ( gridInqType(gridID) != GRID_UNSTRUCTURED && gridInqType(gridID) != GRID_CURVILINEAR )
	    gridID = gridToCurvilinear(gridID, 0);

	  gridtype = gridInqType(gridID);

	  grid_center_lon = (double*) malloc(gridsize*sizeof(double));
	  grid_center_lat = (double*) malloc(gridsize*sizeof(double));
	  gridInqXvals(gridID, grid_center_lon);
	  gridInqYvals(gridID, grid_center_lat);

	  /* Convert lat/lon units if required */
	  {
	    char units[CDI_MAX_NAME];
	    gridInqXunits(gridID, units);
	    grid_to_degree(units, gridsize, grid_center_lon, "grid center lon");
	    gridInqYunits(gridID, units);
	    grid_to_degree(units, gridsize, grid_center_lat, "grid center lat");
	  }
	}

      tsID = 0;
      taxisID = vlistInqTaxis(vlistID);
      while ( (nrecs = streamInqTimestep(streamID, tsID)) )
	{
	  vdate = taxisInqVdate(taxisID);
	  vtime = taxisInqVtime(taxisID);
	  date2str(vdate, vdatestr, sizeof(vdatestr));
	  time2str(vtime, vtimestr, sizeof(vtimestr));

	  cdiDecodeDate(vdate, &year, &month, &day);

	  for ( recID = 0; recID < nrecs; recID++ )
	    {
	      streamInqRecord(streamID, &varID, &levelID);

	      vlistInqVarName(vlistID, varID, name);
	      param    = vlistInqVarParam(vlistID, varID);
	      code     = vlistInqVarCode(vlistID, varID);
	      gridID   = vlistInqVarGrid(vlistID, varID);
	      zaxisID  = vlistInqVarZaxis(vlistID, varID);
	      missval  = vlistInqVarMissval(vlistID, varID);
	      gridsize = gridInqSize(gridID);
	      nlon     = gridInqXsize(gridID);
	      nlat     = gridInqYsize(gridID);
	      level    = zaxisInqLevel(zaxisID, levelID);
	      
	      cdiParamToString(param, paramstr, sizeof(paramstr));

	      if ( nlon*nlat != gridsize ) { nlon = gridsize; nlat = 1; }

	      streamReadRecord(streamID, array, &nmiss);

	      if ( operatorID == OUTPUTSRV )
		fprintf(stdout, "%4d %8g %8d %4d %8d %8d %d %d\n", code, level, vdate, vtime, nlon, nlat, 0, 0);

	      if ( operatorID == OUTPUTEXT )
		fprintf(stdout, "%8d %4d %8g %8d\n", vdate, code, level, gridsize);
		
	      if ( operatorID == OUTPUTINT )
		{
		  nout = 0;
		  for ( i = 0; i < gridsize; i++ )
		    {
		      if ( nout == 8 )
			{
			  nout = 0;
			  fprintf(stdout, "\n");
			}
		      fprintf(stdout, " %8d", (int) array[i]);
		      nout++;
		    }
		  fprintf(stdout, "\n");
		}
	      else if ( operatorID == OUTPUTF )
		{
		  nout = 0;
		  for ( i = 0; i < gridsize; i++ )
		    {
		      if ( nout == nelem )
			{
			  nout = 0;
			  fprintf(stdout, "\n");
			}
		      fprintf(stdout, format, array[i]);
		      nout++;
		    }
		  fprintf(stdout, "\n");
		}
	      else if ( operatorID == OUTPUTTS )
		{
		  char vdatestr[32], vtimestr[32];
	  
		  if ( gridsize > 1 )
		    cdoAbort("operator works only with one gridpoint!");

		  date2str(vdate, vdatestr, sizeof(vdatestr));
		  time2str(vtime, vtimestr, sizeof(vtimestr));

		  fprintf(stdout, "%s %s %12.12g\n", vdatestr, vtimestr, array[0]);
		}
	      else if ( operatorID == OUTPUTFLD )
		{
		  int hour, minute, second;
		  cdiDecodeTime(vtime, &hour, &minute, &second);
		  xdate  = vdate - (vdate/100)*100 + (hour*3600 + minute*60 + second)/86400.;
		  for ( i = 0; i < gridsize; i++ )
		    if ( !DBL_IS_EQUAL(array[i], missval) )
		      fprintf(stdout, "%g\t%g\t%g\t%g\n", xdate, 
			      grid_center_lat[i], grid_center_lon[i], array[i]);
		}
	      else if ( operatorID == OUTPUTTAB )
		{
		  int xind, yind;
		  int l2d = FALSE;

		  int xsize = gridInqXsize(gridID);
		  // int ysize = gridInqYsize(gridID);
		  if ( gridtype == GRID_CURVILINEAR ) l2d = TRUE;
		      
		  for ( i = 0; i < gridsize; i++ )
		    {
		      yind = i;
		      xind = i;
		      if ( l2d ) { yind /= xsize; xind -= yind*xsize; }
		      lon = grid_center_lon[i];
		      lat = grid_center_lat[i];

		      for ( k = 0; k < nkeys; ++k )
			{
			  len = Keylen[keys[k]];
			  if      ( keys[k] == kvalue    ) fprintf(stdout, "%*g ", len, array[i]);
			  else if ( keys[k] == kparam    ) fprintf(stdout, "%*s ", len, paramstr);
			  else if ( keys[k] == kcode     ) fprintf(stdout, "%*d ", len, code);
			  else if ( keys[k] == kname     ) fprintf(stdout, "%*s ", len, name);
			  else if ( keys[k] == klon      ) fprintf(stdout, "%*g ", len, lon);
			  else if ( keys[k] == klat      ) fprintf(stdout, "%*g ", len, lat);
			  else if ( keys[k] == klev      ) fprintf(stdout, "%*g ", len, level);
			  else if ( keys[k] == kbin      ) fprintf(stdout, "%*g ", len, level);
			  else if ( keys[k] == kxind     ) fprintf(stdout, "%*d ", len, xind+1);
			  else if ( keys[k] == kyind     ) fprintf(stdout, "%*d ", len, yind+1);
			  else if ( keys[k] == ktimestep ) fprintf(stdout, "%*d ", len, tsID+1);
			  else if ( keys[k] == kdate     ) fprintf(stdout, "%*s ", len, vdatestr);
			  else if ( keys[k] == ktime     ) fprintf(stdout, "%*s ", len, vtimestr);
			  else if ( keys[k] == kyear     ) fprintf(stdout, "%*d ", len, year);
			  else if ( keys[k] == kmonth    ) fprintf(stdout, "%*d ", len, month);
			  else if ( keys[k] == kday      ) fprintf(stdout, "%*d ", len, day);
			}
		      fprintf(stdout, "\n");
		    }
		}
	      else if ( operatorID == OUTPUTXYZ )
		{
		  if ( tsID == 0 && recID == 0 )
		    {
		      char *fname = "frontplane.xyz";
		      FILE *fp;
		      double fmin = 0;
		      double dx, x0, y0, z0, x, y, z;
		      for ( i = 0; i < gridsize; i++ )
			if ( !DBL_IS_EQUAL(array[i], missval) )
			  {
			    if ( array[i] < fmin ) fmin = array[i];
			    fprintf(stdout, "%g\t%g\t%g\t%g\n",
				    grid_center_lon[i], grid_center_lat[i], array[i], array[i]);
			  }
		      fp = fopen(fname, "w");
		      if ( fp == NULL ) cdoAbort("Open failed on %s", fname);
		      // first front plane
		      dx = (grid_center_lon[1] - grid_center_lon[0]);
		      x0 = grid_center_lon[0]-dx/2;
		      y0 = grid_center_lat[0]-dx/2;
		      z0 = fmin;
		      fprintf(fp, ">\n");
		      for ( i = 0; i < nlon; ++i )
			{
			  x = x0;  y = y0; z = z0;
			  fprintf(fp, "%g %g %g\n", x, y, z);
			  x = x0;  y = y0; z = array[i];
			  fprintf(fp, "%g %g %g\n", x, y, z);
			  x = x0+dx;  y = y0;
			  fprintf(fp, "%g %g %g\n", x, y, z);
			  x0 = x; /*y0 = y0;*/ z0 = z;
			}
		      x = x0;  y = y0; z = fmin;
		      fprintf(fp, "%g %g %g\n", x, y, z);
		      x = grid_center_lon[0]-dx/2;
		      fprintf(fp, "%g %g %g\n", x, y, z);

		      // second front plane
		      x0 = grid_center_lon[0]-dx/2;
		      y0 = grid_center_lat[0]-dx/2;
		      z0 = fmin;
		      fprintf(fp, ">\n");
		      for ( i = 0; i < nlat; ++i )
			{
			  x = x0;  y = y0; z = z0;
			  fprintf(fp, "%g %g %g\n", x, y, z);
			  x = x0;  y = y0; z = array[i*nlon];
			  fprintf(fp, "%g %g %g\n", x, y, z);
			  x = x0;  y = y0+dx;
			  fprintf(fp, "%g %g %g\n", x, y, z);
			  /*x0 = x0;*/ y0 = y; z0 = z;
			}
		      x = x0;  y = y0; z = fmin;
		      fprintf(fp, "%g %g %g\n", x, y, z);
		      y = grid_center_lat[0]-dx/2;
		      fprintf(fp, "%g %g %g\n", x, y, z);

		      fclose(fp);
		    }
		}
	      else if ( operatorID == OUTPUTARR )
		{
		  for ( i = 0; i < gridsize; i++ )
		    {
		      fprintf(stdout, "  arr[%d] = %12.6g;\n", i, array[i]);
		      nout++;
		    }
		}
	      else
		{
		  double minval, maxval;
		  minval = array[0];
		  maxval = array[0];
		  if ( gridInqType(gridID) == GRID_SPECTRAL && gridsize <= 156 )
		    {
		      for ( i = 1; i < gridsize; i++ )
			{
			  if ( array[i] < minval ) minval = array[i];
			  if ( array[i] > maxval ) maxval = array[i];
			}
		    }

		  if ( gridInqType(gridID) == GRID_SPECTRAL && gridsize <= 156 /* T11 */ &&
		       minval >= -1 && maxval <= 12 )
		    {
		      long m, n, ntr;
		      double *spc = array;
		      ntr = gridInqTrunc(gridID);
		      for ( m = 0; m <= ntr; m++ )
			{
			  for ( n = m; n <= ntr; n++ )
			    {
			      fprintf(stdout, "%3d", (int) *spc++);
			      fprintf(stdout, "%3d", (int) *spc++);
			    }
			  fprintf(stdout, "\n");
			}
		    }
		  else
		    {
		      nout = 0;
		      for ( i = 0; i < gridsize; i++ )
			{
			  if ( nout == 6 )
			    {
			      nout = 0;
			      fprintf(stdout, "\n");
			    }
			  fprintf(stdout, " %12.6g", array[i]);
			  nout++;
			}
		      fprintf(stdout, "\n");
		    }
		}
	    }
	  tsID++;
	}
      streamClose(streamID);

      if ( array ) free(array);
      if ( grid_center_lon ) free(grid_center_lon);
      if ( grid_center_lat ) free(grid_center_lat);
    }

  if ( keys ) free(keys);

  cdoFinish();

  return (0);
}
