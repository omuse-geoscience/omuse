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
#include "grid.h"


static
void lonlat_to_xyz(double lon, double lat, double *xyz)
{
  double coslat = cos(lat);
  xyz[0] = coslat * cos(lon);
  xyz[1] = coslat * sin(lon);
  xyz[2] = sin(lat);
}

/*
#include "clipping/grid.h"
#include "clipping/grid_cell.h"
#include "clipping/area.h"

static
double yac_huiliers_area(int num_corners, double *cell_corner_lon, double *cell_corner_lat)
{
  if ( num_corners < 3 ) return 0;

  double coordinates_xyz[num_corners*3];
  enum edge_type edge_types[num_corners];
  struct grid_cell cell =
    {.coordinates_x   = cell_corner_lon,
     .coordinates_y   = cell_corner_lat,
     .coordinates_xyz = coordinates_xyz,
     .edge_type       = edge_types,
     .num_corners     = num_corners};

  for ( int i = 0; i < num_corners; ++i ) edge_types[i] = GREAT_CIRCLE;
  for ( int i = 0; i < num_corners; ++i ) 
    lonlat_to_xyz(cell_corner_lon[i], cell_corner_lat[i], coordinates_xyz+i*3);

  double area = huiliers_area(cell);
  area /= (EarthRadius*EarthRadius);

  return (area);
}
*/
static
void cross_product(const double *restrict a, const double *restrict b, double *restrict c)
{
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}

static
double scalar_product(const double *restrict a, const double *restrict b)
{
  return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
}

static
double norm(const double *restrict a)
{
  return (a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}

static
double mod_cell_area(int num_corners, double *cell_corner_lon, double *cell_corner_lat)
{
  if ( num_corners < 3 ) return 0;

  /* generalised version based on the ICON code, mo_base_geometry.f90
     provided by Luis Kornblueh, MPI-M. */

  int M = num_corners; // number of vertices
  int m; // loop index over number of vertices M
  int i; // loop index over the three dimensions

  double area;
  double s[M];
  double ca[M];
  double a[M];

  double p[M][3];
  double u[M][3];

  /* Convert into cartesian coordinates */

  for ( m = 0; m < M; m++ )
    lonlat_to_xyz(cell_corner_lon[m], cell_corner_lat[m], p[m]);

  /* First, compute cross products Uij = Vi x Vj. */

  for ( m = 0; m < M; m++ )
    cross_product(p[m], p[(m+1)%M], u[m]);

  /*  Normalize Uij to unit vectors. */

  area = 0.0;

  for ( m = 0; m < M; m++ )
    {
      s[m] = norm(u[m]);
      area += s[m];
    }

  /* Test for a degenerated cells associated with collinear vertices. */

  if ( fabs(area) > 0.0 )
    {
      for ( m = 0; m < M; m++ )
	s[m] = sqrt(s[m]);

      for ( m = 0; m < M; m++ )
	for ( i = 0; i < 3; i++ )
	  u[m][i] = u[m][i]/s[m];

      /*  Compute interior angles Ai as the dihedral angles between planes
	  by using the definition of the scalar product

	            ab = |a| |b| cos (phi)

	  As a and b are already normalised this reduces to

                    ab = cos (phi)

          There is no explanation so far for the - in the loop below.
	  But otherwise we don't get the correct results for triangles
	  and cells. Must have something to do with the theorem.
      */

      for ( m = 0; m < M; m++ )
	{
	  ca[m] = - scalar_product(u[m], u[(m+1)%M]);
	  if ( ca[m] < -1.0 ) ca[m] = -1.0;
	  if ( ca[m] >  1.0 ) ca[m] =  1.0;
	  a[m] = acos(ca[m]);
	}

      /*  Compute areas = a1 + a2 + a3 - (M-2) * pi.

	  here for a unit sphere: */

      area = - (double) (M-2) * M_PI;

      for ( m = 0; m < M; m++ )
	area += a[m];

      //   area *= EarthRadius * EarthRadius;

      if ( area < 0.0 ) area = 0.0;
    }

  return (area);
}

/** area of a spherical triangle based on L'Huilier's Theorem
  *
  * source code is taken from code by Robert Oehmke of Earth System Modeling
  * Framework (www.earthsystemmodeling.org)
  *
  * the license statement for this routine is as follows:
  * Earth System Modeling Framework
  * Copyright 2002-2013, University Corporation for Atmospheric Research,
  * Massachusetts Institute of Technology, Geophysical Fluid Dynamics
  * Laboratory, University of Michigan, National Centers for Environmental
  * Prediction, Los Alamos National Laboratory, Argonne National Laboratory,
  * NASA Goddard Space Flight Center.
  * Licensed under the University of Illinois-NCSA License.
  */
static
double mod_tri_area(const double *restrict u, const double *restrict v, const double *restrict w)
{
  double tmp_vec[3];

  cross_product(u, v, tmp_vec);
  double sina = sqrt(norm(tmp_vec));
  double a = asin(sina);

  cross_product(u, w, tmp_vec);
  double sinb = sqrt(norm(tmp_vec));
  double b = asin(sinb);

  cross_product(w, v, tmp_vec);
  double sinc = sqrt(norm(tmp_vec));
  double c = asin(sinc);

  double s = 0.5*(a+b+c);

  double t = tan(s*0.5) * tan((s - a)*0.5) * tan((s - b)*0.5) * tan((s - c)*0.5);

  double area = fabs(4.0 * atan(sqrt(fabs(t))));

  return (area);
}

 /*
  * source code is taken from code by Robert Oehmke of Earth System Modeling
  * Framework (www.earthsystemmodeling.org) and adjusted to CDO data structures
  *
  * the license statement for this routine is as follows:
  * Earth System Modeling Framework
  * Copyright 2002-2013, University Corporation for Atmospheric Research,
  * Massachusetts Institute of Technology, Geophysical Fluid Dynamics
  * Laboratory, University of Michigan, National Centers for Environmental
  * Prediction, Los Alamos National Laboratory, Argonne National Laboratory,
  * NASA Goddard Space Flight Center.
  * Licensed under the University of Illinois-NCSA License.
  */
static
double mod_huiliers_area(int num_corners, double *cell_corner_lon, double *cell_corner_lat)
{
  if ( num_corners < 3 ) return 0;

  // sum areas around cell
  double sum = 0.0;
  double pnt1[3], pnt2[3], pnt3[3];

  lonlat_to_xyz(cell_corner_lon[0], cell_corner_lat[0], pnt1);
  lonlat_to_xyz(cell_corner_lon[1], cell_corner_lat[1], pnt2);

  for ( int i = 2; i < num_corners; i++ )
    {
      // points that make up a side of cell
      lonlat_to_xyz(cell_corner_lon[i], cell_corner_lat[i], pnt3);
 
      // compute angle for pnt2
      sum += mod_tri_area(pnt1, pnt2, pnt3);

      if ( i < (num_corners-1) ) { pnt2[0] = pnt3[0]; pnt2[1] = pnt3[1]; pnt2[2] = pnt3[2]; }
    }

  return (sum);
}

static
double mod_huiliers_area2(int num_corners, double *cell_corner_lon, double *cell_corner_lat, double cell_center_lon, double cell_center_lat)
{
  if ( num_corners < 3 ) return 0;

  // sum areas around cell
  double sum = 0.0;
  double pnt1[3], pnt2[3], pnt3[3];

  lonlat_to_xyz(cell_center_lon, cell_center_lat, pnt1);
  lonlat_to_xyz(cell_corner_lon[0], cell_corner_lat[0], pnt2);

  for ( int i = 1; i < num_corners; i++ )
    {
      if ( IS_EQUAL(cell_corner_lon[i], cell_corner_lon[i-1]) && IS_EQUAL(cell_corner_lat[i], cell_corner_lat[i-1]) ) continue;

      // points that make up a side of cell
      lonlat_to_xyz(cell_corner_lon[i], cell_corner_lat[i], pnt3);
 
      // compute angle for pnt2
      sum += mod_tri_area(pnt1, pnt2, pnt3);

      pnt2[0] = pnt3[0]; pnt2[1] = pnt3[1]; pnt2[2] = pnt3[2];
    }

  if ( !(IS_EQUAL(cell_corner_lon[0], cell_corner_lon[num_corners-1]) && IS_EQUAL(cell_corner_lat[0], cell_corner_lat[num_corners-1])) )
    {
      lonlat_to_xyz(cell_corner_lon[0], cell_corner_lat[0], pnt3);
      sum += mod_tri_area(pnt1, pnt2, pnt3);
    }

  return (sum);
}


int gridGenArea(int gridID, double* area)
{
  int status = 0;
  int gridtype;
  int lgrid_gen_bounds = FALSE;
  int lgriddestroy = FALSE;
  long i;
  long nv, gridsize;
  int* grid_mask = NULL;
  double* grid_center_lon = NULL;
  double* grid_center_lat = NULL;
  double* grid_corner_lon = NULL;
  double* grid_corner_lat = NULL;

  gridsize = gridInqSize(gridID);
  gridtype = gridInqType(gridID);

  if ( gridtype != GRID_LONLAT      &&
       gridtype != GRID_GAUSSIAN    &&
       gridtype != GRID_LCC         &&
       gridtype != GRID_LCC2        &&
       gridtype != GRID_LAEA        &&
       gridtype != GRID_SINUSOIDAL  &&
       gridtype != GRID_GME         &&
       gridtype != GRID_CURVILINEAR &&
       gridtype != GRID_UNSTRUCTURED )
    {
      cdoAbort("Internal error! Unsupported gridtype: %s", gridNamePtr(gridtype)); 
    }

  if ( gridtype != GRID_UNSTRUCTURED && gridtype != GRID_CURVILINEAR )
    {
      if ( gridtype == GRID_GME )
	{
	  lgriddestroy = TRUE;
	  gridID = gridToUnstructured(gridID, 1);
	  grid_mask = (int*) malloc(gridsize*sizeof(int));
	  gridInqMaskGME(gridID, grid_mask);
	}
      else
	{
	  lgriddestroy = TRUE;
	  gridID = gridToCurvilinear(gridID, 1);
	  lgrid_gen_bounds = TRUE;
	}
    }

  if ( gridtype == GRID_UNSTRUCTURED )
    {
      if ( gridInqYvals(gridID, NULL) == 0 || gridInqXvals(gridID, NULL) == 0 )
	{
	  if ( gridInqNumber(gridID) > 0 )
	    {
	      lgriddestroy = TRUE;
	      gridID = referenceToGrid(gridID);
	      if ( gridID == -1 ) return (1);
	    }
	}
    }

  gridtype = gridInqType(gridID);

  if ( gridtype == GRID_UNSTRUCTURED )
    nv = gridInqNvertex(gridID);
  else
    nv = 4;

  if ( gridInqYvals(gridID, NULL) == 0 || gridInqXvals(gridID, NULL) == 0 )
    {
      cdoWarning("Computation of grid cell area weights failed, grid cell center coordinates missing!");
      status = 1;
      return (status);
    }

  if ( nv == 0 )
    {
      cdoWarning("Computation of grid cell area weights failed, grid cell corner coordinates missing!");
      status = 1;
      return (status);
    }

  char xunitstr[CDI_MAX_NAME];
  char yunitstr[CDI_MAX_NAME];
  gridInqXunits(gridID, xunitstr);
  gridInqYunits(gridID, yunitstr);

  grid_center_lon = (double*) malloc(gridsize*sizeof(double));
  grid_center_lat = (double*) malloc(gridsize*sizeof(double));

  gridInqXvals(gridID, grid_center_lon);
  gridInqYvals(gridID, grid_center_lat);

  grid_corner_lon = (double*) malloc(nv*gridsize*sizeof(double));
  grid_corner_lat = (double*) malloc(nv*gridsize*sizeof(double));

  if ( gridInqYbounds(gridID, NULL) && gridInqXbounds(gridID, NULL) )
    {
      gridInqXbounds(gridID, grid_corner_lon);
      gridInqYbounds(gridID, grid_corner_lat);
    }
  else
    {
      if ( lgrid_gen_bounds )
	{
	  int nlon = gridInqXsize(gridID);
	  int nlat = gridInqYsize(gridID);
	  double dlon = 0;
	  if ( nlon == 1 ) dlon = 1;

	  grid_cell_center_to_bounds_X2D(xunitstr, nlon, nlat, grid_center_lon, grid_corner_lon, dlon);
	  grid_cell_center_to_bounds_Y2D(yunitstr, nlon, nlat, grid_center_lat, grid_corner_lat);
	}
      else
	{
	  status = 1;
	  return (status);
	}
    }
  
  grid_to_radian(xunitstr, gridsize,    grid_center_lon, "grid1 center longitudes");
  grid_to_radian(xunitstr, gridsize*nv, grid_corner_lon, "grid1 corner longitudes");

  grid_to_radian(yunitstr, gridsize,    grid_center_lat, "grid1 center latitudes");
  grid_to_radian(yunitstr, gridsize*nv, grid_corner_lat, "grid1 corner latitudes");

  if ( lgriddestroy ) gridDestroy(gridID);

  double findex = 0;

  progressInit();

#if defined(_OPENMP)
#pragma omp parallel for default(none)        \
  shared(findex, gridsize, area, nv, grid_corner_lon, grid_corner_lat, grid_center_lon, grid_center_lat) \
  private(i)
#endif
  for ( i = 0; i < gridsize; ++i )
    {
      int lprogress = 1;
      if ( cdo_omp_get_thread_num() != 0 ) lprogress = 0;

#if defined(_OPENMP)
#include "pragma_omp_atomic_update.h"
#endif
      findex++;
      if ( lprogress ) progressStatus(0, 1, findex/gridsize);

      //area[i] = mod_cell_area(nv, grid_corner_lon+i*nv, grid_corner_lat+i*nv);
      if ( nv <= 4 )
	area[i] = mod_huiliers_area(nv, grid_corner_lon+i*nv, grid_corner_lat+i*nv);
      else
	area[i] = mod_huiliers_area2(nv, grid_corner_lon+i*nv, grid_corner_lat+i*nv, grid_center_lon[i], grid_center_lat[i]);
    }

  if ( cdoVerbose )
    {
      double total_area = 0;
      for ( i = 0; i < gridsize; ++i ) total_area += area[i];
      cdoPrint("Total area = %g steradians", total_area);
    }

  free(grid_center_lon);
  free(grid_center_lat);
  free(grid_corner_lon);
  free(grid_corner_lat);
  if ( grid_mask ) free(grid_mask);

  return (status);
}
