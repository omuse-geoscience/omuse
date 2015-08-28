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

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "error.h"


static
void compare_lat_reg2d(int ysize, int gridID1, int gridID2)
{
  if ( ysize > 1 )
    {
      double *yvals1, *yvals2;
      
      yvals1 = (double*) malloc(ysize*sizeof(double));
      yvals2 = (double*) malloc(ysize*sizeof(double));

      gridInqYvals(gridID1, yvals1);
      gridInqYvals(gridID2, yvals2);
		
      if ( IS_EQUAL(yvals1[0], yvals2[ysize-1]) &&
	   IS_EQUAL(yvals1[ysize-1], yvals2[0]) )
	{
	  if ( yvals1[0] > yvals2[0] )
	    cdoWarning("Latitude orientation differ! First grid: N->S; second grid: S->N");
	  else
	    cdoWarning("Latitude orientation differ! First grid: S->N; second grid: N->S");
	}
      else
	{
	  for ( int i = 0; i < ysize; ++i )
	    if ( fabs(yvals1[i] - yvals2[i]) > 1.e-5 )
	      {
		// printf("lat %g %g %g\n", yvals1[i], yvals2[i], yvals1[i] - yvals2[i]);
		cdoWarning("Grid latitudes differ!");
		break;
	      }
	}

      free(yvals1);
      free(yvals2);
    }
}

static
void compare_lon_reg2d(int xsize, int gridID1, int gridID2)
{
  if ( xsize > 1 )
    {
      double *xvals1, *xvals2;

      xvals1 = (double*) malloc(xsize*sizeof(double));
      xvals2 = (double*) malloc(xsize*sizeof(double));

      gridInqXvals(gridID1, xvals1);
      gridInqXvals(gridID2, xvals2);
		  
      for ( int i = 0; i < xsize; ++i )
	if ( fabs(xvals1[i] - xvals2[i]) > 1.e-5 )
	  {
	    cdoWarning("Grid longitudes differ!");
	    break;
	  }
      
      free(xvals1);
      free(xvals2);
    }
}

static
void compare_grid_unstructured(int gridID1, int gridID2)
{
  double *xvals1, *xvals2;
  double *yvals1, *yvals2;
  int gridsize;

  gridsize = gridInqSize(gridID1);

  xvals1 = (double*) malloc(gridsize*sizeof(double));
  xvals2 = (double*) malloc(gridsize*sizeof(double));
  yvals1 = (double*) malloc(gridsize*sizeof(double));
  yvals2 = (double*) malloc(gridsize*sizeof(double));

  gridInqXvals(gridID1, xvals1);
  gridInqXvals(gridID2, xvals2);
  gridInqYvals(gridID1, yvals1);
  gridInqYvals(gridID2, yvals2);
		  
  for ( int i = 0; i < gridsize; ++i )
    if ( fabs(xvals1[i] - xvals2[i]) > 1.e-5 || fabs(yvals1[i] - yvals2[i]) > 1.e-5 )
      {
	cdoWarning("Geographic location of some grid points differ!");
	break;
      }
      
  free(xvals1);
  free(xvals2);
  free(yvals1); 
  free(yvals2); 
}

static
void compareGrids(int gridID1, int gridID2)
{
  /* compare grids of first variable */
  int xsize, ysize;

  if ( gridInqType(gridID1) == gridInqType(gridID2) )
    {
      if ( gridInqType(gridID1) == GRID_GAUSSIAN || gridInqType(gridID1) == GRID_LONLAT )
	{
	  xsize = gridInqXsize(gridID1);
	  ysize = gridInqYsize(gridID1);
		
	  if ( ysize == gridInqYsize(gridID2) )
	    compare_lat_reg2d(ysize, gridID1, gridID2);
	  else
	    cdoWarning("ysize of input grids differ!");
		
	  if ( xsize == gridInqXsize(gridID2) )
	    compare_lon_reg2d(xsize, gridID1, gridID2);
	  else
	    cdoWarning("xsize of input grids differ!");
	}
      else if ( gridInqType(gridID1) == GRID_CURVILINEAR || gridInqType(gridID1) == GRID_UNSTRUCTURED )
	{
	  compare_grid_unstructured(gridID1, gridID2);
	}
    }
  else
    {
      cdoWarning("Grids have different types! First grid: %s; second grid: %s",
		 gridNamePtr(gridInqType(gridID1)), gridNamePtr(gridInqType(gridID2)));
    }
}

static
int cmpnames(const void *s1, const void *s2)
{
  const char *name1 = s1;
  const char *name2 = s2;

  return (strcmp(name1, name2));
}


void vlistCompare(int vlistID1, int vlistID2, int flag)
{
  int varID, nvars;
  int lchecknames = FALSE;

  if ( vlistNvars(vlistID1) != vlistNvars(vlistID2) )
    cdoAbort("Input streams have different number of variables per timestep!");

  if ( vlistNrecs(vlistID1) != vlistNrecs(vlistID2) )
    cdoAbort("Input streams have different number of records per timestep!");

  nvars = vlistNvars(vlistID1);

  for ( varID = 0; varID < nvars; varID++ )
    {
      if ( nvars > 1 )
	{
	  if ( flag & CMP_NAME )
	    {
	      char name1[CDI_MAX_NAME], name2[CDI_MAX_NAME];
	      vlistInqVarName(vlistID1, varID, name1);
	      vlistInqVarName(vlistID2, varID, name2);
	      strtolower(name1);
	      strtolower(name2);
	      if ( strcmp(name1, name2) != 0 )
		{
		  cdoWarning("Input streams have different parameters!");
		  lchecknames = TRUE;
		  flag -= CMP_NAME;
		  //    break;
		}
	    }
	}

      if ( flag & CMP_GRIDSIZE )
	{
	  if ( gridInqSize(vlistInqVarGrid(vlistID1, varID)) !=
	       gridInqSize(vlistInqVarGrid(vlistID2, varID)) )
	    cdoAbort("Grid size of the input parameters do not match!");
	}
      
      if ( flag & CMP_NLEVEL )
	{
	  if ( zaxisInqSize(vlistInqVarZaxis(vlistID1, varID)) !=
	       zaxisInqSize(vlistInqVarZaxis(vlistID2, varID)) )
	    cdoAbort("Number of levels of the input parameters do not match!");
	}
    }

  if ( flag & CMP_GRID )
    {
      int gridID1, gridID2;

      gridID1 = vlistInqVarGrid(vlistID1, 0);
      gridID2 = vlistInqVarGrid(vlistID2, 0);

      compareGrids(gridID1, gridID2);
    }

  if ( lchecknames )
    {
      char names1[nvars][CDI_MAX_NAME], names2[nvars][CDI_MAX_NAME];
      for ( varID = 0; varID < nvars; varID++ )
	vlistInqVarName(vlistID1, varID, names1[varID]);
      for ( varID = 0; varID < nvars; varID++ )
	vlistInqVarName(vlistID2, varID, names2[varID]);

      qsort(names1[0], nvars, CDI_MAX_NAME, cmpnames);

      for ( varID = 0; varID < nvars; varID++ )
	if ( strcmp(names1[varID], names2[varID]) != 0 ) break;

      if ( varID == nvars )
	cdoPrint("Use the CDO option -Q to sort the parameter names, if you have netCDF input files!");
    }
}


int vlistCompareX(int vlistID1, int vlistID2, int flag)
{
  int varID, nvars, nvars2, nlevels2;

  nvars = vlistNvars(vlistID1);
  nvars2 = vlistNvars(vlistID2);
  nlevels2 = zaxisInqSize(vlistInqVarZaxis(vlistID2, 0));

  if ( nvars2 != 1 )
    cdoAbort("Internal problem, vlistCompareX() called with unexpected vlistID2 argument!");

  for ( varID = 0; varID < nvars; varID++ )
    {
      if ( flag & CMP_GRIDSIZE )
	{
	  if ( gridInqSize(vlistInqVarGrid(vlistID1, varID)) !=
	       gridInqSize(vlistInqVarGrid(vlistID2, 0)) )
	    cdoAbort("Grid size of the input parameters do not match!");
	}
      
      if ( flag & CMP_NLEVEL )
	{
	  if ( (zaxisInqSize(vlistInqVarZaxis(vlistID1, varID)) !=
                nlevels2) && nlevels2 > 1 )
	    cdoAbort("Number of levels of the input parameters do not match!");
	}
    }

  if ( flag & CMP_GRID )
    {
      int gridID1, gridID2;

      gridID1 = vlistInqVarGrid(vlistID1, 0);
      gridID2 = vlistInqVarGrid(vlistID2, 0);

      compareGrids(gridID1, gridID2);
    }

  return (nlevels2);
}


int vlistIsSzipped(int vlistID)
{
  int lszip = FALSE;
  int nvars, varID, comptype;

  nvars = vlistNvars(vlistID);

  for ( varID = 0; varID < nvars; varID++ )
    {						
      comptype = vlistInqVarCompType(vlistID, varID);
      if ( comptype == COMPRESS_SZIP )
	{
	  lszip = TRUE;
	  break;
	}
    }      

  return (lszip);
}


int vlistInqNWPV(int vlistID, int varID)
{
  int nwpv; // number of words per value; real:1  complex:2

  if ( vlistInqVarDatatype(vlistID, varID) == DATATYPE_CPX32 || 
       vlistInqVarDatatype(vlistID, varID) == DATATYPE_CPX64 )
    nwpv = 2;
  else
    nwpv = 1;

  return (nwpv);
}


int vlist_check_gridsize(int vlistID)
{
  int lerror = FALSE;
  int ngrids = vlistNgrids(vlistID);
  int gridID = vlistGrid(vlistID, 0);
  int ngp    = gridInqSize(gridID);

  /* check gridsize */
  for ( int index = 0; index < ngrids; ++index )
    {
      gridID = vlistGrid(vlistID, index);
      if ( ngp != gridInqSize(gridID) )
	{
	  lerror = TRUE;
	  break;
	}
    }

  if ( lerror )
    {
      cdoPrint("This operator requires all variables on the same horizontal grid.");
      cdoPrint("Horizontal grids found:");
      for ( int index = 0; index < ngrids; ++index )
	{
	  gridID = vlistGrid(vlistID, index);
	  cdoPrint("  grid=%d  type=%s  points=%d", index+1, gridNamePtr(gridInqType(gridID)), gridInqSize(gridID));
	}
      cdoAbort("The input stream contains variables on different horizontal grids!");
    }

  return ngp;
}
