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

#if defined(HAVE_CONFIG_H)
#  include "config.h"
#endif

#if defined(_OPENMP)
#  include <omp.h>
#endif

#include <stdio.h>
#include <stdarg.h> /* va_list */

#if defined(HAVE_LIBPROJ)
#  include "proj_api.h"
#endif

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "error.h"
#include "grid.h"

int referenceToGrid(int gridID1);


static
void scale_vec(double scalefactor, long nvals, double *restrict values)
{
  long n;

#if defined(_OPENMP)
#pragma omp parallel for default(none) shared(nvals, scalefactor, values)
#endif
  for ( n = 0; n < nvals; ++n )
    {
      values[n] *= scalefactor;
    }
}


void grid_to_radian(const char *units, long nvals, double *restrict values, const char *description)
{
  if ( cmpstr(units, "degree") == 0 )
    {
      scale_vec(DEG2RAD, nvals, values);
    }
  else if ( cmpstr(units, "radian") == 0 )
    {
      /* No conversion necessary */
    }
  else
    {
      cdoWarning("Unknown units [%s] supplied for %s; proceeding assuming radians!", units, description);
    }
}


void grid_to_degree(const char *units, long nvals, double *restrict values, const char *description)
{
  if ( cmpstr(units, "radian") == 0 )
    {
      for ( long n = 0; n < nvals; ++n ) values[n] *= RAD2DEG;
    }
  else if ( cmpstr(units, "degree") == 0 )
    {
      /* No conversion necessary */
    }
  else
    {
      cdoWarning("Unknown units [%s] supplied for %s; proceeding assuming degress!", units, description);
    }
}


int gridToZonal(int gridID1)
{
  int gridID2;
  int gridtype, gridsize;
  double  xval = 0;
  double *yvals;

  gridtype = gridInqType(gridID1);
  gridsize = gridInqYsize(gridID1);
  gridID2  = gridCreate(gridtype, gridsize);
	  
  if ( gridtype == GRID_LONLAT   ||
       gridtype == GRID_GAUSSIAN ||
       gridtype == GRID_GENERIC )
    {
      gridDefXsize(gridID2, 1);
      gridDefYsize(gridID2, gridsize);

      gridDefXvals(gridID2, &xval);

      if ( gridInqYvals(gridID1, NULL) )
	{
	  yvals = (double*) malloc(gridsize*sizeof(double));

	  gridInqYvals(gridID1, yvals);
	  gridDefYvals(gridID2, yvals);

	  free(yvals);
	}
    }
  else
    {
      Error("Gridtype %s unsupported!", gridNamePtr(gridtype));
    }

  return (gridID2);
}


int gridToMeridional(int gridID1)
{
  int gridID2;
  int gridtype, gridsize;
  double *xvals;
  double  yval = 0;

  gridtype = gridInqType(gridID1);
  gridsize = gridInqXsize(gridID1);
  gridID2  = gridCreate(gridtype, gridsize);
	  
  if ( gridtype == GRID_LONLAT   ||
       gridtype == GRID_GAUSSIAN ||
       gridtype == GRID_GENERIC )
    {
      gridDefXsize(gridID2, gridsize);
      gridDefYsize(gridID2, 1);

      if ( gridInqXvals(gridID1, NULL) )
	{
	  xvals = (double*) malloc(gridsize*sizeof(double));

	  gridInqXvals(gridID1, xvals);
	  gridDefXvals(gridID2, xvals);

	  free(xvals);
	}

      gridDefYvals(gridID2, &yval);
    }
  else
    {
      Error("Gridtype %s unsupported!", gridNamePtr(gridtype));
    }

  return (gridID2);
}


void grid_gen_corners(long n, const double* restrict vals, double* restrict corners)
{
  long i;

  if ( n == 1 )
    {
      corners[0] = vals[0];
      corners[1] = vals[0];
    }
  else
    {
      for ( i = 0; i < n-1; ++i )
	corners[i+1] = 0.5*(vals[i] + vals[i+1]);

      corners[0] = 2*vals[0] - corners[1];
      corners[n] = 2*vals[n-1] - corners[n-1];
    }
}


void grid_gen_bounds(long n, const double* restrict vals, double* restrict bounds)
{
  long i;

  for ( i = 0; i < n-1; ++i )
    {
      bounds[2*i+1]   = 0.5*(vals[i] + vals[i+1]);
      bounds[2*(i+1)] = 0.5*(vals[i] + vals[i+1]);
    }

  bounds[0]     = 2*vals[0] - bounds[1];
  bounds[2*n-1] = 2*vals[n-1] - bounds[2*(n-1)];
}


void grid_check_lat_borders(int n, double *ybounds)
{
  if ( ybounds[0] > ybounds[n-1] )
    {
      if ( ybounds[0]   >  88 ) ybounds[0]   =  90;
      if ( ybounds[n-1] < -88 ) ybounds[n-1] = -90;
    }
  else
    {
      if ( ybounds[0]   < -88 ) ybounds[0]   = -90;
      if ( ybounds[n-1] >  88 ) ybounds[n-1] =  90;
    }
}


void grid_cell_center_to_bounds_X2D(const char* xunitstr, long xsize, long ysize, const double* restrict grid_center_lon, 
				    double* restrict grid_corner_lon, double dlon)
{
  long i, j, index;
  double minlon, maxlon;

  if ( ! (dlon > 0) ) dlon = 360./xsize;
  /*
  if ( xsize == 1 || (grid_center_lon[xsize-1]-grid_center_lon[0]+dlon) < 359 )
    cdoAbort("Cannot calculate Xbounds for %d vals with dlon = %g", xsize, dlon);
  */
  for ( i = 0; i < xsize; ++i )
    {
      minlon = grid_center_lon[i] - 0.5*dlon;
      maxlon = grid_center_lon[i] + 0.5*dlon;
      for ( j = 0; j < ysize; ++j )
	{
	  index = (j<<2)*xsize + (i<<2);
	  grid_corner_lon[index  ] = minlon;
	  grid_corner_lon[index+1] = maxlon;
	  grid_corner_lon[index+2] = maxlon;
	  grid_corner_lon[index+3] = minlon;
	}
    }
}

static
double genYmin(double y1, double y2)
{
  double ymin, dy;

  dy = y2 - y1;
  ymin = y1 - dy/2;

  if ( y1 < -85 && ymin < -87.5 ) ymin = -90;

  if ( cdoVerbose )
    cdoPrint("genYmin: y1 = %g  y2 = %g  dy = %g  ymin = %g", y1, y2, dy, ymin);

  return (ymin);
}

static
double genYmax(double y1, double y2)
{
  double ymax, dy;

  dy = y1 - y2;
  ymax = y1 + dy/2;

  if ( y1 > 85 && ymax > 87.5 ) ymax = 90;

  if ( cdoVerbose )
    cdoPrint("genYmax: y1 = %g  y2 = %g  dy = %g  ymax = %g", y1, y2, dy, ymax);

  return (ymax);
}



/*****************************************************************************/

void grid_cell_center_to_bounds_Y2D(const char* yunitstr, long xsize, long ysize, const double* restrict grid_center_lat, double* restrict grid_corner_lat)
{
  long i, j, index;
  double minlat, maxlat;
  double firstlat, lastlat;

  firstlat = grid_center_lat[0];
  lastlat  = grid_center_lat[xsize*ysize-1];

  // if ( ysize == 1 ) cdoAbort("Cannot calculate Ybounds for 1 value!");

  for ( j = 0; j < ysize; ++j )
    {
      if ( ysize == 1 )
	{
	  minlat = grid_center_lat[0] - 360./ysize;
	  maxlat = grid_center_lat[0] + 360./ysize;
	}
      else
	{
	  index = j*xsize;
	  if ( firstlat > lastlat )
	    {
	      if ( j == 0 )
		maxlat = genYmax(grid_center_lat[index], grid_center_lat[index+xsize]);
	      else
		maxlat = 0.5*(grid_center_lat[index]+grid_center_lat[index-xsize]);

	      if ( j == (ysize-1) )
		minlat = genYmin(grid_center_lat[index], grid_center_lat[index-xsize]);
	      else
		minlat = 0.5*(grid_center_lat[index]+grid_center_lat[index+xsize]);
	    }
	  else
	    {
	      if ( j == 0 )
		minlat = genYmin(grid_center_lat[index], grid_center_lat[index+xsize]);
	      else
		minlat = 0.5*(grid_center_lat[index]+grid_center_lat[index-xsize]);

	      if ( j == (ysize-1) )
		maxlat = genYmax(grid_center_lat[index], grid_center_lat[index-xsize]);
	      else
		maxlat = 0.5*(grid_center_lat[index]+grid_center_lat[index+xsize]);
	    }
	}

      for ( i = 0; i < xsize; ++i )
	{
	  index = (j<<2)*xsize + (i<<2);
	  grid_corner_lat[index  ] = minlat;
	  grid_corner_lat[index+1] = minlat;
	  grid_corner_lat[index+2] = maxlat;
	  grid_corner_lat[index+3] = maxlat;
	}
    }
}


void gridGenRotBounds(int gridID, int nx, int ny,
		      double *xbounds, double *ybounds, double *xbounds2D, double *ybounds2D)
{
  long i, j, index;
  double minlon, maxlon;
  double minlat, maxlat;
  double xpole, ypole, angle;

  xpole = gridInqXpole(gridID);
  ypole = gridInqYpole(gridID);
  angle = gridInqAngle(gridID);

  for ( j = 0; j < ny; j++ )
    {
      if ( ybounds[0] > ybounds[1] )
	{
	  maxlat = ybounds[2*j];
	  minlat = ybounds[2*j+1];
	}
      else
	{
	  maxlat = ybounds[2*j+1];
	  minlat = ybounds[2*j];
	}

      for ( i = 0; i < nx; i++ )
	{
	  minlon = xbounds[2*i];
	  maxlon = xbounds[2*i+1];

	  index = j*4*nx + 4*i;
	  xbounds2D[index+0] = lamrot_to_lam(minlat, minlon, ypole, xpole, angle);
	  xbounds2D[index+1] = lamrot_to_lam(minlat, maxlon, ypole, xpole, angle);
	  xbounds2D[index+2] = lamrot_to_lam(maxlat, maxlon, ypole, xpole, angle);
	  xbounds2D[index+3] = lamrot_to_lam(maxlat, minlon, ypole, xpole, angle);

	  ybounds2D[index+0] = phirot_to_phi(minlat, minlon, ypole, angle);
	  ybounds2D[index+1] = phirot_to_phi(minlat, maxlon, ypole, angle);
	  ybounds2D[index+2] = phirot_to_phi(maxlat, maxlon, ypole, angle);
	  ybounds2D[index+3] = phirot_to_phi(maxlat, minlon, ypole, angle);
	}
    }
}

static
void gridGenXbounds2D(long nx, long ny, const double* restrict xbounds, double* restrict xbounds2D)
{
  long i, j, index;
  double minlon, maxlon;

#if defined(_OPENMP)
#pragma omp parallel for default(none)        \
  shared(nx, ny, xbounds, xbounds2D)	      \
  private(i, j, minlon, maxlon, index)
#endif
  for ( i = 0; i < nx; ++i )
    {
      minlon = xbounds[2*i  ];
      maxlon = xbounds[2*i+1];

      for ( j = 0; j < ny; ++j )
	{
	  index = j*4*nx + 4*i;
	  xbounds2D[index  ] = minlon;
	  xbounds2D[index+1] = maxlon;
	  xbounds2D[index+2] = maxlon;
	  xbounds2D[index+3] = minlon;
	}
    }
}

static
void gridGenYbounds2D(long nx, long ny, const double* restrict ybounds, double* restrict ybounds2D)
{
  long i, j, index;
  double minlat, maxlat;

#if defined(_OPENMP)
#pragma omp parallel for default(none)        \
  shared(nx, ny, ybounds, ybounds2D)	      \
  private(i, j, minlat, maxlat, index)
#endif
  for ( j = 0; j < ny; ++j )
    {
      if ( ybounds[0] > ybounds[1] )
	{
	  maxlat = ybounds[2*j  ];
	  minlat = ybounds[2*j+1];
	}
      else
	{
	  maxlat = ybounds[2*j+1];
	  minlat = ybounds[2*j  ];
	}

      for ( i = 0; i < nx; ++i )
	{
	  index = j*4*nx + 4*i;
	  ybounds2D[index  ] = minlat;
	  ybounds2D[index+1] = minlat;
	  ybounds2D[index+2] = maxlat;
	  ybounds2D[index+3] = maxlat;
	}
    }
}

static
char *gen_param(const char *fmt, ...)
{
  va_list args;
  char str[256];
  char *rstr;
  int len;

  va_start(args, fmt);

  len = vsprintf(str, fmt, args);

  va_end(args);

  len++;
  rstr = (char*) malloc(len*sizeof(char));
  memcpy(rstr, str, len*sizeof(char));

  return (rstr);
}

static
void lcc_to_geo(int gridID, int gridsize, double *xvals, double *yvals)
{
  double originLon, originLat, lonParY, lat1, lat2, xincm, yincm;
  double zlat, zlon;
  double xi, xj;
  int projflag, scanflag;
  long i;
  proj_info_t proj;

  gridInqLCC(gridID, &originLon, &originLat, &lonParY, &lat1, &lat2, &xincm, &yincm, &projflag, &scanflag);
  /*
    while ( originLon < 0 ) originLon += 360;
    while ( lonParY   < 0 ) lonParY   += 360;
  */
  if ( IS_NOT_EQUAL(xincm, yincm) )
    Warning("X and Y increment must be equal on Lambert Conformal grid (Xinc = %g, Yinc = %g)\n", 
	    xincm, yincm);
  /*
  if ( IS_NOT_EQUAL(lat1, lat2) )
    Warning("Lat1 and Lat2 must be equal on Lambert Conformal grid (Lat1 = %g, Lat2 = %g)\n", 
	    lat1, lat2);
  */
  map_set(PROJ_LC, originLat, originLon, xincm, lonParY, lat1, lat2, &proj);

  for ( i = 0; i < gridsize; i++ )
    {
      xi = xvals[i];
      xj = yvals[i];
      // status = W3FB12(xi, xj, originLat, originLon, xincm, lonParY, lat1, &zlat, &zlon);
      ijll_lc(xi, xj, proj, &zlat, &zlon);
      xvals[i] = zlon;
      yvals[i] = zlat;
    }
}

static
void sinusoidal_to_geo(int gridsize, double *xvals, double *yvals)
{
#if defined(HAVE_LIBPROJ)
  projPJ   *libProj;
  char *params[20];
  int nbpar=0;
  projUV data, res;
  long i;

  nbpar = 0;
  params[nbpar++] = (char*) "proj=sinu";
  params[nbpar++] = (char*) "ellps=WGS84";

  if ( cdoVerbose )
    for ( i = 0; i < nbpar; ++i )
      cdoPrint("Proj.param[%ld] = %s", i+1, params[i]);

  libProj = pj_init(nbpar, params);
  if ( !libProj )
    cdoAbort("proj error: %s", pj_strerrno(pj_errno));

  /* libProj->over = 1; */		/* allow longitude > 180° */

  for ( i = 0; i < gridsize; i++ )
    {
      data.u = xvals[i];
      data.v = yvals[i];
      res = pj_inv(data, libProj);
      xvals[i] = res.u*RAD2DEG;
      yvals[i] = res.v*RAD2DEG;
      if ( xvals[i] < -9000. || xvals[i] > 9000. ) xvals[i] = -9999.;
      if ( yvals[i] < -9000. || yvals[i] > 9000. ) yvals[i] = -9999.;
    }
#else
  cdoAbort("proj4 support not compiled in!");
#endif
}

static
void laea_to_geo(int gridID, int gridsize, double *xvals, double *yvals)
{
#if defined(HAVE_LIBPROJ)
  projPJ   *libProj;
  char *params[20];
  int nbpar=0;
  projUV data, res;
  double a, lon_0, lat_0;
  long i;

  gridInqLaea(gridID, &a , &lon_0, &lat_0);

  nbpar = 0;
  params[nbpar++] = gen_param("proj=laea");
  if ( a > 0 ) params[nbpar++] = gen_param("a=%g", a);
  params[nbpar++] = gen_param("lon_0=%g", lon_0);
  params[nbpar++] = gen_param("lat_0=%g", lat_0);

  if ( cdoVerbose )
    for ( i = 0; i < nbpar; ++i )
      cdoPrint("Proj.param[%d] = %s", i+1, params[i]);

  libProj = pj_init(nbpar, &params[0]);
  if ( !libProj )
    cdoAbort("proj error: %s", pj_strerrno(pj_errno));

  for ( i = 0; i < nbpar; ++i ) free(params[i]);

  /* libProj->over = 1; */		/* allow longitude > 180° */

  for ( i = 0; i < gridsize; i++ )
    {
      data.u = xvals[i];
      data.v = yvals[i];
      res = pj_inv(data, libProj);
      xvals[i] = res.u*RAD2DEG;
      yvals[i] = res.v*RAD2DEG;
      if ( xvals[i] < -9000. || xvals[i] > 9000. ) xvals[i] = -9999.;
      if ( yvals[i] < -9000. || yvals[i] > 9000. ) yvals[i] = -9999.;
    }
#else
  cdoAbort("proj4 support not compiled in!");
#endif
}

static
void lcc2_to_geo(int gridID, int gridsize, double *xvals, double *yvals)
{
#if defined(HAVE_LIBPROJ)
  projPJ   *libProj;
  char *params[20];
  int nbpar=0;
  projUV data, res;
  double a, lon_0, lat_0, lat_1, lat_2;
  long i;

  gridInqLcc2(gridID, &a , &lon_0, &lat_0, &lat_1, &lat_2);

  nbpar = 0;
  params[nbpar++] = gen_param("proj=lcc");
  if ( a > 0 ) params[nbpar++] = gen_param("a=%g", a);
  params[nbpar++] = gen_param("lon_0=%g", lon_0);
  params[nbpar++] = gen_param("lat_0=%g", lat_0);
  params[nbpar++] = gen_param("lat_1=%g", lat_1);
  params[nbpar++] = gen_param("lat_2=%g", lat_2);

  if ( cdoVerbose )
    for ( i = 0; i < nbpar; ++i )
      cdoPrint("Proj.param[%d] = %s", i+1, params[i]);
  
  libProj = pj_init(nbpar, &params[0]);
  if ( !libProj )
    cdoAbort("proj error: %s", pj_strerrno(pj_errno));

  for ( i = 0; i < nbpar; ++i ) free(params[i]);

  /* libProj->over = 1; */		/* allow longitude > 180° */
  
  for ( i = 0; i < gridsize; i++ )
    {
      data.u = xvals[i];
      data.v = yvals[i];
      res = pj_inv(data, libProj);
      xvals[i] = res.u*RAD2DEG;
      yvals[i] = res.v*RAD2DEG;
    }
#else
  cdoAbort("proj4 support not compiled in!");
#endif
}

/*
 * grib_get_reduced_row: code from GRIB_API 1.10.4
 *
 * Description:
 *   computes the number of points within the range lon_first->lon_last and the zero based indexes
 *   ilon_first,ilon_last of the first and last point given the number of points along a parallel (pl)
 *
 */
static
void grib_get_reduced_row(long pl,double lon_first,double lon_last,long* npoints,long* ilon_first, long* ilon_last )
{
  double range=0,dlon_first=0,dlon_last=0;
  long irange;

  range=lon_last-lon_first;
  if (range<0) {range+=360;lon_first-=360;}

  /* computing integer number of points and coordinates without using floating point resolution*/
  *npoints=(range*pl)/360.0+1;
  *ilon_first=(lon_first*pl)/360.0;
  *ilon_last=(lon_last*pl)/360.0;

  irange=*ilon_last-*ilon_first+1;

  if (irange != *npoints) {
    if (irange > *npoints) {
      /* checking if the first point is out of range*/
      dlon_first=((*ilon_first)*360.0)/pl;
      if (dlon_first < lon_first) {(*ilon_first)++;irange--;
      }

      /* checking if the last point is out of range*/
      dlon_last=((*ilon_last)*360.0)/pl;
      if (dlon_last > lon_last) {(*ilon_last)--;irange--;
      }
    } else {
      int ok=0;
      /* checking if the point before the first is in the range*/
      dlon_first=((*ilon_first-1)*360.0)/pl;
      if (dlon_first > lon_first) {(*ilon_first)--;irange++;ok=1;
      }

      /* checking if the point after the last is in the range*/
      dlon_last=((*ilon_last+1)*360.0)/pl;
      if (dlon_last < lon_last) {(*ilon_last)++;irange++;ok=1;
      }

      /* if neither of the two are triggered then npoints is too large */
      if (!ok) {(*npoints)--;
      }
    }

    //   assert(*npoints==irange);
  } else {
	  /* checking if the first point is out of range*/
	  dlon_first=((*ilon_first)*360.0)/pl;
	  if (dlon_first < lon_first) {
		  (*ilon_first)++;(*ilon_last)++;
	  }
  }

  if (*ilon_first<0) *ilon_first+=pl;

  return;
}


int qu2reg3_double(double *pfield, int *kpoint, int klat, int klon,
		   double msval, int *kret, int omisng, int operio, int oveggy);

static
int qu2reg_subarea(int gridsize, int np, double xfirst, double xlast, 
		   double *array, int *rowlon, int ny, double missval, int *iret, int lmiss, int lperio, int lveggy)
{
  int nx = 0;
  /* sub area (longitudes) */
  long ilon_firstx;
  long ilon_first, ilon_last;
  int i, j;
  long row_count;
  int rlon, nwork = 0;
  int np4 = np*4;
  int size = 0;
  int wlen;
  double *work, **pwork;
  int ii;
  double *parray;

  if ( np <= 0 ) cdoAbort("Number of values between pole and equator missing!");

  grib_get_reduced_row(np4, xfirst, xlast, &row_count, &ilon_firstx, &ilon_last);
  nx = row_count;
  // printf("nx %d  %ld %ld lon1 %g lon2 %g\n", nx, ilon_firstx, ilon_last, (ilon_firstx*360.)/np4, (ilon_last*360.)/np4);

  for ( j = 0; j < ny; ++j ) nwork += rowlon[j];

  pwork = (double **) malloc(ny*sizeof(double *));
  work  = (double*) malloc(ny*np4*sizeof(double));
  wlen = 0;
  pwork[0] = work;
  for ( j = 1; j < ny; ++j )
    {
      wlen += rowlon[j-1];
      pwork[j] = work + wlen;
    } 
  // printf(" ny, np4, nwork %d %d %d wlen %d\n", ny, np4, nwork, wlen);

  for ( j = 0; j < ny; ++j )
    {
      rlon = rowlon[j];
      for ( i = 0; i < rlon; ++i ) pwork[j][i] = missval;
    } 

  parray = array;
  for ( j = 0; j < ny; ++j )
    {
      rlon = rowlon[j];
      row_count = 0;
      grib_get_reduced_row(rlon, xfirst, xlast, &row_count, &ilon_first, &ilon_last);
      // printf("j %d xfirst %g xlast %g rowlon %d %ld %ld %ld %g %g\n", j, xfirst, xlast, rlon, row_count, ilon_first, ilon_last, (ilon_first*360.)/rlon, (ilon_last*360.)/rlon);

      for ( i = ilon_first; i < (ilon_first+row_count); ++i )
	{
	  ii = i;
	  if ( ii >= rlon ) ii -= rlon; 
	  pwork[j][ii] = *parray;
	  parray++;
	}
      size += row_count;
    }

  if ( gridsize != size ) cdoAbort("gridsize1 inconsistent!");

  (void) qu2reg3_double(work, rowlon, ny, np4, missval, iret, lmiss, lperio, lveggy);

  wlen = 0;
  pwork[0] = work;
  for ( j = 1; j < ny; ++j )
    {
      wlen += np4;
      pwork[j] = work + wlen;
    } 

  // printf("nx, ilon_firstx %d %ld\n", nx, ilon_firstx);
  parray = array;
  for ( j = 0; j < ny; ++j )
    {
      for ( i = ilon_firstx; i < (ilon_firstx+nx); ++i )
	{
	  ii = i;
	  if ( ii >= np4 ) ii -= np4; 
	  *parray = pwork[j][ii];
	  parray++;
	}
    }

  free(work);
  free(pwork);

  return (nx);
}


void field2regular(int gridID1, int gridID2, double missval, double *array, int nmiss)
{
  int nx = 0, ny, np;
  int gridtype;
  int lmiss, lperio, lveggy;
  int iret;
  int *rowlon;
  double xfirstandlast[2];
  double xfirst, xlast;

  gridtype = gridInqType(gridID1);

  if ( gridtype != GRID_GAUSSIAN_REDUCED ) Error("Not a reduced gaussian grid!");

  lmiss = nmiss > 0;
  lperio = 1;
  lveggy = 0;

  ny = gridInqYsize(gridID1);
  np = gridInqNP(gridID1);

  rowlon = (int*) malloc(ny*sizeof(int));
  gridInqRowlon(gridID1, rowlon);

  xfirstandlast[0] = 0.;
  xfirstandlast[1] = 0.;
  gridInqXvals(gridID1, xfirstandlast);
  xfirst = xfirstandlast[0];
  xlast  = xfirstandlast[1];

  if ( fabs(xfirst) > 0 || (np > 0 && fabs(xlast - (360.0-90.0/np)) > 90.0/np) )
    {
      nx = qu2reg_subarea(gridInqSize(gridID1), np, xfirst, xlast, array, rowlon, ny, missval, &iret, lmiss, lperio, lveggy);
    }
  else
    {
      nx = 2*ny;
      (void) qu2reg3_double(array, rowlon, ny, nx, missval, &iret, lmiss, lperio, lveggy);
    }

  if ( gridInqSize(gridID2) != nx*ny ) Error("Gridsize differ!");

  free(rowlon);
}


int gridToRegular(int gridID1)
{
  int gridID2;
  int gridtype, gridsize;
  int nx = 0, ny, np;
  long i;
  double *xvals = NULL, *yvals = NULL;
  double xfirstandlast[2];
  double xfirst, xlast;

  gridtype = gridInqType(gridID1);

  if ( gridtype != GRID_GAUSSIAN_REDUCED ) Error("Not a reduced gaussian grid!");

  ny = gridInqYsize(gridID1);
  np = gridInqNP(gridID1);

  yvals = (double*) malloc(ny*sizeof(double));
  gridInqYvals(gridID1, yvals);

  xfirstandlast[0] = 0.;
  xfirstandlast[1] = 0.;
  gridInqXvals(gridID1, xfirstandlast);
  xfirst = xfirstandlast[0];
  xlast  = xfirstandlast[1];

  if ( fabs(xfirst) > 0 || (np > 0 && fabs(xlast - (360.0-90.0/np)) > 90.0/np) )
    {
      /* sub area (longitudes) */
      long ilon_first, ilon_last;
      long row_count;
      int np4 = np*4;
      int *rowlon = NULL;

      if ( np <= 0 ) cdoAbort("Number of values between pole and equator missing!");

      grib_get_reduced_row(np4, xfirst, xlast, &row_count, &ilon_first, &ilon_last);

      nx = row_count;
      xvals = (double*) malloc(nx*sizeof(double));
      for ( i = 0; i < nx; ++i )
	{
	  xvals[i] = ((ilon_first+i)*360.)/np4;
	  if ( xfirst > xlast ) xvals[i] -= 360.;
	}

      free(rowlon);
    }
  else
    {
      nx = 2*ny;
      xvals = (double*) malloc(nx*sizeof(double));
      for ( i = 0; i < nx; ++i ) xvals[i] = i * 360./nx;
    }

  gridsize = nx*ny;

  gridID2  = gridCreate(GRID_GAUSSIAN, gridsize);
	  
  gridDefXsize(gridID2, nx);
  gridDefYsize(gridID2, ny);
  
  gridDefXvals(gridID2, xvals);
  gridDefYvals(gridID2, yvals);
  gridDefNP(gridID2, np);

  free(xvals);
  free(yvals);

  return (gridID2);
}

static
void gridCopyMask(int gridID1, int gridID2, long gridsize)
{
  if ( gridInqMask(gridID1, NULL) )
    {
      int *mask;
      mask = (int*) malloc(gridsize*sizeof(int));
      gridInqMask(gridID1, mask);
      gridDefMask(gridID2, mask);
      free(mask);
    }
}

static
int check_range(long n, double *vals, double valid_min, double valid_max)
{
  int status = 0;
  long i;

  for ( i = 0; i < n; ++i )
    {
      if ( vals[i] < valid_min || vals[i] > valid_max )
	{
	  status = 1;
	  break;
	}
    }

  return (status);
}


int gridToCurvilinear(int gridID1, int lbounds)
{
  int gridID2;
  int gridtype;
  size_t gridsize;
  long index;

  gridtype = gridInqType(gridID1);
  gridsize = (size_t) gridInqSize(gridID1);
  gridID2  = gridCreate(GRID_CURVILINEAR, (int) gridsize);
  gridDefPrec(gridID2, DATATYPE_FLT32);
	  
  switch (gridtype)
    {
    case GRID_LONLAT:
    case GRID_GAUSSIAN:
    case GRID_LCC:
    case GRID_LCC2:
    case GRID_LAEA:
    case GRID_SINUSOIDAL:
      {
	long i, j;
	int nx, ny;
	double *xvals = NULL, *yvals = NULL;
	double *xvals2D, *yvals2D;
	double *xbounds = NULL, *ybounds = NULL;
	double *xbounds2D, *ybounds2D;
	char xunits[CDI_MAX_NAME], yunits[CDI_MAX_NAME];
	double xscale = 1, yscale = 1;

	nx = gridInqXsize(gridID1);
	ny = gridInqYsize(gridID1);

	gridInqXunits(gridID1, xunits);
	gridInqYunits(gridID1, yunits);

	if ( gridtype == GRID_LAEA )
	  {
	    int lvalid_xunits = FALSE;
	    int lvalid_yunits = FALSE;
	    int len;
	    len = (int) strlen(xunits);
	    if ( len == 1 && memcmp(xunits, "m",  1) == 0 ) lvalid_xunits = TRUE;
	    if ( len == 2 && memcmp(xunits, "km", 2) == 0 ) lvalid_xunits = TRUE;
	    len = (int) strlen(yunits);
	    if ( len == 1 && memcmp(yunits, "m",  1) == 0 ) lvalid_yunits = TRUE;
	    if ( len == 2 && memcmp(yunits, "km", 2) == 0 ) lvalid_yunits = TRUE;

	    if ( lvalid_xunits == FALSE )
	      cdoWarning("Possibly wrong result! Invalid x-coordinate units: \"%s\" (expected \"m\" or \"km\")", xunits);
	    if ( lvalid_yunits == FALSE )
	      cdoWarning("Possibly wrong result! Invalid y-coordinate units: \"%s\" (expected \"m\" or \"km\")", yunits);
	  }

	if ( memcmp(xunits, "km", 2) == 0 ) xscale = 1000;
	if ( memcmp(yunits, "km", 2) == 0 ) yscale = 1000;

	gridDefXsize(gridID2, nx);
	gridDefYsize(gridID2, ny);

	xvals2D = (double*) malloc(gridsize*sizeof(double));
	yvals2D = (double*) malloc(gridsize*sizeof(double));


	if ( gridtype == GRID_LCC )
	  {
	    for ( j = 0; j < ny; j++ )
	      for ( i = 0; i < nx; i++ )
		{
		  xvals2D[j*nx+i] = i+1;
		  yvals2D[j*nx+i] = j+1;
		}

	    lcc_to_geo(gridID1, gridsize, xvals2D, yvals2D);
	  }
	else
	  {
	    if ( ! (gridInqXvals(gridID1, NULL) && gridInqYvals(gridID1, NULL)) )
	      Error("Grid has no values");

	    xvals = (double*) malloc(nx*sizeof(double));
	    yvals = (double*) malloc(ny*sizeof(double));

	    gridInqXvals(gridID1, xvals);
	    gridInqYvals(gridID1, yvals);

	    if ( gridIsRotated(gridID1) )
	      {
		double xpole, ypole, angle;
		
		xpole = gridInqXpole(gridID1);
		ypole = gridInqYpole(gridID1);
		angle = gridInqAngle(gridID1);
		
		for ( j = 0; j < ny; j++ )
		  for ( i = 0; i < nx; i++ )
		    {
		      xvals2D[j*nx+i] = lamrot_to_lam(yvals[j], xvals[i], ypole, xpole, angle);
		      yvals2D[j*nx+i] = phirot_to_phi(yvals[j], xvals[i], ypole, angle);
		    }	    
	      }
	    else
	      {
		for ( j = 0; j < ny; j++ )
		  for ( i = 0; i < nx; i++ )
		    {
		      xvals2D[j*nx+i] = xscale*xvals[i];
		      yvals2D[j*nx+i] = yscale*yvals[j];
		    }

		if ( gridtype == GRID_SINUSOIDAL )
		  {
		    sinusoidal_to_geo(gridsize, xvals2D, yvals2D);
		    /* correct_sinxvals(nx, ny, xvals2D); */
		  }
		else if ( gridtype == GRID_LAEA )
		  {
		    laea_to_geo(gridID1, gridsize, xvals2D, yvals2D);
		  }
		else if ( gridtype == GRID_LCC2 )
		  {
		    lcc2_to_geo(gridID1, gridsize, xvals2D, yvals2D);
		  }
	      }
	  }

	gridDefXvals(gridID2, xvals2D);
	gridDefYvals(gridID2, yvals2D);

	if ( xvals2D ) free(xvals2D);
	if ( yvals2D ) free(yvals2D);

	if ( !lbounds ) goto NO_BOUNDS;

	if ( gridtype == GRID_LCC )
	  {		
	    xbounds2D = (double*) malloc(4*gridsize*sizeof(double));
	    ybounds2D = (double*) malloc(4*gridsize*sizeof(double));

	    for ( j = 0; j < ny; j++ )
	      for ( i = 0; i < nx; i++ )
		{
		  index = j*4*nx + 4*i;

		  xbounds2D[index+0] = i+1.5;
		  ybounds2D[index+0] = j+1.5;

		  xbounds2D[index+1] = i+0.5;
		  ybounds2D[index+1] = j+1.5;

		  xbounds2D[index+2] = i+0.5;
		  ybounds2D[index+2] = j+0.5;

		  xbounds2D[index+3] = i+1.5;
		  ybounds2D[index+3] = j+0.5;
		}

	    lcc_to_geo(gridID1, 4*gridsize, xbounds2D, ybounds2D);

	    gridDefXbounds(gridID2, xbounds2D);
	    gridDefYbounds(gridID2, ybounds2D);

	    free(xbounds2D);
	    free(ybounds2D);
	  }
	else
	  {
	    if ( gridInqXbounds(gridID1, NULL) )
	      {
		xbounds = (double*) malloc(2*nx*sizeof(double));
		gridInqXbounds(gridID1, xbounds);
		if ( gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN )
		  if ( check_range(2*nx, xbounds, -720, 720) )
		    {
		      cdoWarning("longitude bounds out of range, skipped!");
		      free(xbounds);
		      xbounds = NULL;
		    }
	      }
	    else if ( nx > 1 )
	      {
		xbounds = (double*) malloc(2*nx*sizeof(double));
		grid_gen_bounds(nx, xvals, xbounds);
	      }

	    if ( gridInqYbounds(gridID1, NULL) )
	      {
		ybounds = (double*) malloc(2*ny*sizeof(double));
		gridInqYbounds(gridID1, ybounds);
		if ( gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN )
		  if ( check_range(2*ny, ybounds, -180, 180) )
		    {
		      cdoWarning("latitude bounds out of range, skipped!");
		      free(ybounds);
		      ybounds = NULL;
		    }
	      }
	    else if ( ny > 1 )
	      {
		ybounds = (double*) malloc(2*ny*sizeof(double));
		if ( gridtype == GRID_SINUSOIDAL || 
		     gridtype == GRID_LAEA       || 
		     gridtype == GRID_LCC2 )
		  grid_gen_bounds(ny, yvals, ybounds);
		else
		  {
		    grid_gen_bounds(ny, yvals, ybounds);
		    grid_check_lat_borders(2*ny, ybounds);
		  }
	      }

	    if ( xbounds && ybounds )
	      {
		xbounds2D = (double*) malloc(4*gridsize*sizeof(double));
		ybounds2D = (double*) malloc(4*gridsize*sizeof(double));

		if ( gridIsRotated(gridID1) )
		  {
		    gridGenRotBounds(gridID1, nx, ny, xbounds, ybounds, xbounds2D, ybounds2D);
		  }
		else
		  {
		    if ( gridtype == GRID_SINUSOIDAL ||
			 gridtype == GRID_LAEA       || 
			 gridtype == GRID_LCC2 )
		      {
			for ( j = 0; j < ny; j++ )
			  for ( i = 0; i < nx; i++ )
			    {
			      index = j*4*nx + 4*i;

			      xbounds2D[index+0] = xscale*xbounds[2*i];
			      ybounds2D[index+0] = yscale*ybounds[2*j];

			      xbounds2D[index+1] = xscale*xbounds[2*i];
			      ybounds2D[index+1] = yscale*ybounds[2*j+1];

			      xbounds2D[index+2] = xscale*xbounds[2*i+1];
			      ybounds2D[index+2] = yscale*ybounds[2*j+1];

			      xbounds2D[index+3] = xscale*xbounds[2*i+1];
			      ybounds2D[index+3] = yscale*ybounds[2*j];
			    }
			
			if ( gridtype == GRID_SINUSOIDAL )
			  {
			    sinusoidal_to_geo(4*gridsize, xbounds2D, ybounds2D);
			    /*
			    xvals2D = (double*) malloc(gridsize*sizeof(double));
			    for ( j = 0; j < 4; ++j )
			      {
				for ( i = 0; i < gridsize; ++i ) xvals2D[i] = xbounds2D[i*4+j];
				correct_sinxvals(nx, ny, xvals2D);
				for ( i = 0; i < gridsize; ++i ) xbounds2D[i*4+j] = xvals2D[i];
			      }
			    free(xvals2D);
			    */
			  }
			else if ( gridtype == GRID_LAEA )
			  laea_to_geo(gridID1, 4*gridsize, xbounds2D, ybounds2D);
			else if ( gridtype == GRID_LCC2 )
			  lcc2_to_geo(gridID1, 4*gridsize, xbounds2D, ybounds2D);
		      }
		    else
		      {
			gridGenXbounds2D(nx, ny, xbounds, xbounds2D);
			gridGenYbounds2D(nx, ny, ybounds, ybounds2D);
		      }
		  }
		
		gridDefXbounds(gridID2, xbounds2D);
		gridDefYbounds(gridID2, ybounds2D);
		
		if ( xbounds )  free(xbounds);
		if ( ybounds )  free(ybounds);
		if ( xbounds2D) free(xbounds2D);
		if ( ybounds2D) free(ybounds2D);
	      }
	  }

      NO_BOUNDS:

	if ( xvals ) free(xvals);
	if ( yvals ) free(yvals);

	gridCopyMask(gridID1, gridID2, gridsize);

	break;
      }
    default:
      {
	Error("Grid type >%s< unsupported!", gridNamePtr(gridtype));
	break;
      }
    }

  return (gridID2);
}


int gridToUnstructured(int gridID1, int lbounds)
{
  int gridID2;
  int gridtype, gridsize;

  gridtype = gridInqType(gridID1);
  gridsize = gridInqSize(gridID1);
  gridID2  = gridCreate(GRID_UNSTRUCTURED, gridsize);
  gridDefPrec(gridID2, DATATYPE_FLT32);
	  
  switch (gridtype)
    {
    case GRID_LONLAT:
    case GRID_GAUSSIAN:
      {
	long i, j;
	int nx, ny;
	double *xvals, *yvals;
	double *xvals2D, *yvals2D;

	gridDefXname(gridID2, "lon");
	gridDefYname(gridID2, "lat");
	gridDefXlongname(gridID2, "longitude");
	gridDefYlongname(gridID2, "latitude");
	gridDefXunits(gridID2, "degrees_east");
	gridDefYunits(gridID2, "degrees_north");

	gridDefNvertex(gridID2, 4);

	nx = gridInqXsize(gridID1);
	ny = gridInqYsize(gridID1);
	 
	gridDefXsize(gridID2, gridsize);
	gridDefYsize(gridID2, gridsize);

	xvals = (double*) malloc(nx*sizeof(double));
	yvals = (double*) malloc(ny*sizeof(double));

	xvals2D = (double*) malloc(gridsize*sizeof(double));
	yvals2D = (double*) malloc(gridsize*sizeof(double));

	gridInqXvals(gridID1, xvals);
	gridInqYvals(gridID1, yvals);

	if ( gridIsRotated(gridID1) )
	  {
	    double xpole, ypole, angle;
	    
	    xpole = gridInqXpole(gridID1);
	    ypole = gridInqYpole(gridID1);
	    angle = gridInqAngle(gridID1);
		
	    for ( j = 0; j < ny; j++ )
	      for ( i = 0; i < nx; i++ )
		{
		  xvals2D[j*nx+i] = lamrot_to_lam(yvals[j], xvals[i], ypole, xpole, angle);
		  yvals2D[j*nx+i] = phirot_to_phi(yvals[j], xvals[i], ypole, angle);
		}
	  }
	else
	  {
	    for ( j = 0; j < ny; j++ )
	      for ( i = 0; i < nx; i++ )
		{
		  xvals2D[j*nx+i] = xvals[i];
		  yvals2D[j*nx+i] = yvals[j];
		}
	  }

	gridDefXvals(gridID2, xvals2D);
	gridDefYvals(gridID2, yvals2D);

	free(xvals2D);
	free(yvals2D);

	if ( lbounds )
	  {
	    double *xbounds = NULL, *ybounds = NULL;
	    double *xbounds2D, *ybounds2D;

	    if ( gridInqXbounds(gridID1, NULL) )
	      {
		xbounds = (double*) malloc(2*nx*sizeof(double));
		gridInqXbounds(gridID1, xbounds);
	      }
	    else if ( nx > 1 )
	      {
		xbounds = (double*) malloc(2*nx*sizeof(double));
		grid_gen_bounds(nx, xvals, xbounds);
	      }

	    if ( gridInqYbounds(gridID1, NULL) )
	      {
		ybounds = (double*) malloc(2*ny*sizeof(double));
		gridInqYbounds(gridID1, ybounds);
	      }
	    else if ( ny > 1 )
	      {
		ybounds = (double*) malloc(2*ny*sizeof(double));
		grid_gen_bounds(ny, yvals, ybounds);
		grid_check_lat_borders(2*ny, ybounds);
	      }

	    if ( xbounds && ybounds )
	      {
		xbounds2D = (double*) malloc(4*gridsize*sizeof(double));
		ybounds2D = (double*) malloc(4*gridsize*sizeof(double));

		if ( gridIsRotated(gridID1) )
		  {
		    gridGenRotBounds(gridID1, nx, ny, xbounds, ybounds, xbounds2D, ybounds2D);
		  }
		else
		  {
		    gridGenXbounds2D(nx, ny, xbounds, xbounds2D);
		    gridGenYbounds2D(nx, ny, ybounds, ybounds2D);
		  }

		gridDefXbounds(gridID2, xbounds2D);
		gridDefYbounds(gridID2, ybounds2D);

		free(xbounds);
		free(ybounds);
		free(xbounds2D);
		free(ybounds2D);
	      }
	    }

	free(xvals);
	free(yvals);

	gridCopyMask(gridID1, gridID2, gridsize);

	break;
      }
    case GRID_CURVILINEAR:
      {
	gridID2 = gridDuplicate(gridID1);
	gridChangeType(gridID2, GRID_UNSTRUCTURED);
	gridDefXsize(gridID2, gridsize);
	gridDefYsize(gridID2, gridsize);

	break;
      }
    case GRID_GME:
      {
	int nd, ni, ni2, ni3;
	long i, j;
	int nv = 6;
	int *imask;
	double *xvals, *yvals;
	double *xbounds = NULL, *ybounds = NULL;

	nd  = gridInqGMEnd(gridID1);
	ni  = gridInqGMEni(gridID1);
	ni2 = gridInqGMEni2(gridID1);
	ni3 = gridInqGMEni3(gridID1);

	imask   = (int*) malloc(gridsize*sizeof(int));
	xvals   = (double*) malloc(gridsize*sizeof(double));
	yvals   = (double*) malloc(gridsize*sizeof(double));
	if ( lbounds )
	  {
	    xbounds = (double*) malloc(nv*gridsize*sizeof(double));
	    ybounds = (double*) malloc(nv*gridsize*sizeof(double));
	  }

	gme_grid(lbounds, gridsize, xvals, yvals, xbounds, ybounds, imask, ni, nd, ni2, ni3);
	
	for ( i = 0; i < gridsize; i++ )
	  {
	    xvals[i] *= RAD2DEG;
	    yvals[i] *= RAD2DEG;

	    if ( lbounds )
	      for ( j = 0; j < nv; j++ )
		{
		  xbounds[i*nv + j] *= RAD2DEG;
		  ybounds[i*nv + j] *= RAD2DEG;
		}
	    /* printf("%d %g %g\n", i, xvals[i], yvals[i]); */
	  }
	
	gridDefXsize(gridID2, gridsize);
	gridDefYsize(gridID2, gridsize);

	gridDefXvals(gridID2, xvals);
	gridDefYvals(gridID2, yvals);

	gridDefMaskGME(gridID2, imask);

	gridDefNvertex(gridID2, nv);

	if ( lbounds )
	  {
	    gridDefXbounds(gridID2, xbounds);
	    gridDefYbounds(gridID2, ybounds);
	  }

	gridDefXunits(gridID2, "degrees_east");
	gridDefYunits(gridID2, "degrees_north");

	free (imask);
	free (xvals);
	free (yvals);
	if ( xbounds ) free (xbounds);
	if ( ybounds ) free (ybounds);
	
	gridCopyMask(gridID1, gridID2, gridsize);

	break;
      }
    default:
      {
	Error("Grid type >%s< unsupported!", gridNamePtr(gridtype));
	break;
      }
    }

  return (gridID2);
}


int gridCurvilinearToRegular(int gridID1)
{
  int gridID2 = -1;
  int gridtype, gridsize;
  long i, j;
  int nx, ny;
  int lx = TRUE, ly = TRUE;
  double *xvals = NULL, *yvals = NULL;
  double *xvals2D = NULL, *yvals2D = NULL;
	
  gridtype = gridInqType(gridID1);
  gridsize = gridInqSize(gridID1);

  if ( gridtype != GRID_CURVILINEAR ) return (gridID2);

  nx = gridInqXsize(gridID1);
  ny = gridInqYsize(gridID1);
	
  xvals2D = (double*) malloc(gridsize*sizeof(double));
  yvals2D = (double*) malloc(gridsize*sizeof(double));

  gridInqXvals(gridID1, xvals2D);
  gridInqYvals(gridID1, yvals2D);

  xvals = (double*) malloc(nx*sizeof(double));
  yvals = (double*) malloc(ny*sizeof(double));

  for ( i = 0; i < nx; i++ ) xvals[i] = xvals2D[i];
  for ( j = 0; j < ny; j++ ) yvals[j] = yvals2D[j*nx];

  for ( j = 1; j < ny; j++ )
    for ( i = 0; i < nx; i++ )
      {
	if ( fabs(xvals[i] - xvals2D[j*nx+i]) > 1.e-6 )
	  {
	    lx = FALSE;
	    j = ny;
	    break;
	  }
      }
	
  for ( i = 1; i < nx; i++ )
    for ( j = 0; j < ny; j++ )
      {
	if ( fabs(yvals[j] - yvals2D[j*nx+i]) > 1.e-6 )
	  {
	    ly = FALSE;
	    i = nx;
	    break;
	  }
      }

  free(xvals2D);
  free(yvals2D);

  if ( lx && ly )
    {
      char xunits[CDI_MAX_NAME], yunits[CDI_MAX_NAME];
      
      gridID2  = gridCreate(GRID_LONLAT, gridsize);
      gridDefXsize(gridID2, nx);
      gridDefYsize(gridID2, ny);
      
      //  gridDefPrec(gridID2, DATATYPE_FLT32);

      gridInqXunits(gridID1, xunits);
      gridInqYunits(gridID1, yunits);

      grid_to_degree(xunits, nx, xvals, "grid1 center lon");
      grid_to_degree(yunits, ny, yvals, "grid1 center lat");

      gridDefXvals(gridID2, xvals);
      gridDefYvals(gridID2, yvals);
    }

  free(xvals);
  free(yvals);

  return (gridID2);
}


int gridGenWeights(int gridID, double *grid_area, double *grid_wgts)
{
  int i, nvals, gridsize, gridtype;
  int status = 0;
  int *grid_mask = NULL;
  double total_area;

  gridtype = gridInqType(gridID);
  gridsize = gridInqSize(gridID);
  
  if ( gridtype == GRID_GME )
    {
      gridID = gridToUnstructured(gridID, 1);	  
      grid_mask = (int*) malloc(gridsize*sizeof(int));
      gridInqMaskGME(gridID, grid_mask);
    }

  total_area = 0;
  nvals = 0;
  for ( i = 0; i < gridsize; i++ )
    {
      if ( grid_mask )
	if ( grid_mask[i] == 0 ) continue;
      total_area += grid_area[i];
      nvals++;
    }

  if ( cdoVerbose ) cdoPrint("Total area = %g steradians", total_area);

  for ( i = 0; i < gridsize; i++ )
    {
      if ( grid_mask )
	if ( grid_mask[i] == 0 )
	  {
	    grid_wgts[i] = 0;
	    continue;
	  }
      
      grid_wgts[i] = grid_area[i] / total_area;
    }
  
  if ( grid_mask ) free(grid_mask);

  return (status);
}


int gridWeightsOld(int gridID, double *weights)
{
  int status = FALSE;
  long i, j;
  int len;

  len = gridInqSize(gridID);

  if ( gridHasArea(gridID) )
    {
      gridInqArea(gridID, weights);
    }
  else
    {
      int gridtype = gridInqType(gridID);

      if ( gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN )
	{
	  int     nlat, nlon;
	  int     datapoint;
	  double *lats = NULL, *lons = NULL;
	  double sumw;
	  double phi1, phi2, theta1, theta2, sindphi;

	  nlon = gridInqXsize(gridID);
	  nlat = gridInqYsize(gridID);

	  lons = 1 + (double *) malloc((nlon+2)*sizeof(double));
	  lats = 1 + (double *) malloc((nlat+2)*sizeof(double));

	  gridInqXvals(gridID, lons);
	  gridInqYvals(gridID, lats);

	  /* Interpolate to find latitudes outside boundaries. */
	  lats[-1]   = 2*lats[0] - lats[1];
	  lats[nlat] = 2*lats[nlat-1] - lats[nlat-2];
	  lons[-1]   = 2*lons[0] - lons[1];
	  lons[nlon] = 2*lons[nlon-1] - lons[nlon-2];
  
	  /*  Calculate weights.  */
	  /*  phi 1 and 2 and theta 1 and 2 represent respectively the boundary */
	  /*  latitudes and longitudes of a particular grid square.             */
	  datapoint = 0;
	  sumw = 0;
	  for ( j = 0; j < nlat; j++ )
	    {
	      phi1 = (lats[j-1]+lats[j])/2*DEG2RAD;
	      phi2 = (lats[j+1]+lats[j])/2*DEG2RAD;
	      if ( phi1 < (-1*M_PI/2) ) phi1 = -1*M_PI/2;
	      if ( phi1 > (   M_PI/2) ) phi1 =    M_PI/2;
	      if ( phi2 > (   M_PI/2) ) phi2 =    M_PI/2;
	      if ( phi2 < (-1*M_PI/2) ) phi2 = -1*M_PI/2;
	      sindphi = sin(phi2)-sin(phi1);
	      for( i = 0; i < nlon; i++ )
		{
		  if ( lons[i] >= lons[0]+360 || fabs(lats[j]) > 90 )
		    weights[datapoint] = 0;
		  else
		    {
		      theta1 = (lons[i-1]+lons[i])/2*DEG2RAD;
		      theta2 = (lons[i+1]+lons[i])/2*DEG2RAD;
		      weights[datapoint] = fabs((theta2-theta1)*sindphi);
		      sumw += weights[datapoint];
		    }
		  datapoint++;
		}
	    }

	  /* Normalise weights.  */
	  if( IS_NOT_EQUAL(sumw, 0) )
	    for( i = 0; i < datapoint; i++ ) weights[i] /= sumw;

	  if ( lons-1 ) free(lons-1);
	  if ( lats-1 ) free(lats-1);
	}
      else
	{
	  status = TRUE;

	  for ( i = 0; i < len; i++ ) weights[i] = 1./len;
	}
    }

  return (status);
}


int gridWeights(int gridID, double *grid_wgts)
{
  int i, gridsize, gridtype;
  int a_status, w_status;
  double *grid_area;

  gridtype = gridInqType(gridID);
  gridsize = gridInqSize(gridID);
  
  grid_area = (double*) malloc(gridsize*sizeof(double));

  a_status = 0;

  if ( gridHasArea(gridID) )
    {
      if ( cdoVerbose ) cdoPrint("Using existing grid cell area!");
      gridInqArea(gridID, grid_area);
    }
  else
    {
      if ( gridtype == GRID_LONLAT      ||
	   gridtype == GRID_GAUSSIAN    ||
	   gridtype == GRID_LCC         ||
	   gridtype == GRID_LCC2        ||
	   gridtype == GRID_LAEA        ||
	   gridtype == GRID_SINUSOIDAL  ||
	   gridtype == GRID_GME         ||
	   gridtype == GRID_CURVILINEAR ||
	   gridtype == GRID_UNSTRUCTURED )
	{
	  a_status = gridGenArea(gridID, grid_area);
	}
      else
	{
	  a_status = 1;
	}
    }

  if ( a_status == 0 )
    {
      w_status = gridGenWeights(gridID, grid_area, grid_wgts);
    }
  else
    {
      for ( i = 0; i < gridsize; ++i )
	grid_wgts[i] = 1./gridsize;

      w_status = 1;
    }
  /*
  for ( i = 0; i < gridsize; ++i ) 
    printf("weights: %d %d %d %g %g\n", a_status, w_status, i, grid_area[i], grid_wgts[i]);
  */
  free(grid_area);

  return (w_status);
}
