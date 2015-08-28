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

#define H5_USE_16_API

#if defined(HAVE_LIBHDF5)
#  include "hdf5.h"
#endif

#include <ctype.h>

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


#define  MAX_DSETS  1024

typedef struct {
  char *name;
  char *description;
  char *units;
  char *title;
  char *time;
  int dtype;
  int nx;
  int ny;
  int nz;
  int nt;
  int gridsize;
  int lscale;
  int loffset;
  int lmissval;
  double scale;
  double offset;
  double missval;
  double *array;
}
dset_obj_t;

typedef struct {
  int nsets;
  int mergelevel;
  int lgeoloc;
  int lregion;
  int lprojtype;
  int lmetadata;
  dset_obj_t obj[MAX_DSETS];
}
datasets_t;


#if defined(HAVE_LIBHDF5)
static
void print_filter(hid_t dset_id, char *varname)
{
  hid_t plist;
  H5Z_filter_t filter;
  unsigned int flags;
  int idx;
  unsigned int cd_values;
  int nfilter;
  size_t cd_nelmts = 1;
  size_t pnamelen = 64;
  char pname[64];

  /* get filter */
  plist = H5Dget_create_plist(dset_id);
  nfilter = H5Pget_nfilters(plist);

  for ( idx = 0; idx < nfilter; idx++ )
    {
      filter = H5Pget_filter(plist, idx, &flags, &cd_nelmts, &cd_values,
			     pnamelen, pname);
      cdoPrint("Dataset %s: filter %d =  %s", varname, idx+1, pname);
    }

  H5Pclose(plist);
}

static
void get_grid_info(double c0, double re, int *nrxp, int *nryp, double *r0p, double *s0p, double *cp)
{
  const double pi = M_PI;
  double git, phi, s90;
  double r0, s0, c;
  int nrx, nry;

  git=2.*pi*re*cos(pi/6.)/c0;
  /* number of longitude pixels */
  nrx=2*(int)lround(0.5*git);

  /* central index in longitude */
  r0=nrx/2+0.5;

  /* resolution in km */
  c=2.*pi*re*cos(30.*pi/180.)/nrx;

  phi=pi/2.;
  s90= re/c *sin(phi) / cos(30.*pi/180.);

  nry=(int)floor(s90);
  /* central index in latitude */
  s0=nry+0.5;
  /* number of latitude pixels */
  nry=2*nry;

  *nrxp = nrx;
  *nryp = nry;
  *r0p  = r0;
  *s0p  = s0;
  *cp   = c;
}

static
double det_lon_atovs(double r, double r0, double lts, double c, double re)
{
  const double pi = M_PI;
  double xla;

  xla=(r-r0)*c/re/cos(lts*pi/180.); /* longitude */
  xla=180.*xla/pi;

  return (xla);
}

static
double det_lat_atovs(double s, double s0, double lts, double c, double re)
{
  const double pi = M_PI;
  double siphi;
  double phi;

  siphi=(s-s0)*c*cos(lts*pi/180.)/re;
  phi=180.*asin(siphi)/pi; /* latitude */

  return (phi);
}

static
int defLonLatGrid(int nx, int ny, double c0, double lts, double re)
{
  int gridID;
  int nrx, nry, i;
  double c;
  double r0, s0;
  double r, s;
  double xla, phi;
  double *xvals, *yvals, *xbounds, *ybounds;

  get_grid_info(c0, re, &nrx, &nry, &r0, &s0, &c);

  if ( nx != nrx || ny != nry )
    {
      printf("nrx=%d nry=%d\n", nrx, nry);
      return(-1);
    }

  xvals = (double*) malloc(nx*sizeof(double));
  yvals = (double*) malloc(ny*sizeof(double));
  xbounds = (double*) malloc(nx*2*sizeof(double));
  ybounds = (double*) malloc(nx*2*sizeof(double));

  for ( i = 0; i < nx; ++i )
    {
      r = i+1;
      xla = det_lon_atovs(r, r0, lts, c, re);
      xvals[i] = xla;
      xla = det_lon_atovs(r-0.5, r0, lts, c, re);
      xbounds[2*i] = xla;
      xla = det_lon_atovs(r+0.5, r0, lts, c, re);
      xbounds[2*i+1] = xla;
      /* printf("xla[%d]=%g\n", i, xla); */
    }

  for ( i = 0; i < ny; ++i )
    {
      s = (nry-i-1)+1;
      phi = det_lat_atovs(s, s0, lts, c, re);
      yvals[i] = phi;
      phi = det_lat_atovs(s-0.5, s0, lts, c, re);
      ybounds[2*i] = phi;
      phi = det_lat_atovs(s+0.5, s0, lts, c, re);
      ybounds[2*i+1] = phi;
      /* printf("phi[%d]=%g\n", i, phi); */
    }

  gridID = gridCreate(GRID_LONLAT, nx*ny);
  gridDefXsize(gridID, nx);
  gridDefYsize(gridID, ny);
  gridDefXvals(gridID, xvals);
  gridDefYvals(gridID, yvals);
  /*
  gridDefXbounds(gridID, xbounds);
  gridDefYbounds(gridID, ybounds);
  */
  free(xvals);
  free(yvals);
  free(xbounds);
  free(ybounds);

  return (gridID);
}

static
int defSinusoidalGrid(int nx, int ny, double xmin, double xmax, double ymin, double ymax,
		      double dx, double dy, double p1, double p2, double p3, double p4)
{
  int gridID;
  int i;
  double *xvals, *yvals;

  xvals = (double*) malloc(nx*sizeof(double));
  yvals = (double*) malloc(ny*sizeof(double));

  for ( i = 0; i < nx; ++i )
    {
      xvals[i] = xmin + i*dx + dx/2;
      /* printf("x[%d]=%g\n", i, xvals[i]); */
    }

  for ( i = 0; i < ny; ++i )
    {
      yvals[i] = ymax - i*dx - dx/2;
      /* printf("y[%d]=%g\n", i, yvals[i]); */
    }

  gridID = gridCreate(GRID_SINUSOIDAL, nx*ny);
  gridDefXsize(gridID, nx);
  gridDefYsize(gridID, ny);
  gridDefXvals(gridID, xvals);
  gridDefYvals(gridID, yvals);

  free(xvals);
  free(yvals);

  return (gridID);
}

static
int defLaeaGrid(int nx, int ny, double xmin, double xmax, double ymin, double ymax,
		double dx, double dy, double a, double lon0, double lat0)
{
  int gridID;
  int i;
  double *xvals, *yvals;

  xvals = (double*) malloc(nx*sizeof(double));
  yvals = (double*) malloc(ny*sizeof(double));

  for ( i = 0; i < nx; ++i )
    {
      xvals[i] = xmin + i*dx + dx/2;
      /* printf("x[%d]=%g\n", i, xvals[i]); */
    }

  for ( i = 0; i < ny; ++i )
    {
      yvals[i] = ymax - i*dx - dx/2;
      /* printf("y[%d]=%g\n", i, yvals[i]); */
    }

  gridID = gridCreate(GRID_LAEA, nx*ny);
  gridDefXsize(gridID, nx);
  gridDefYsize(gridID, ny);
  gridDefXvals(gridID, xvals);
  gridDefYvals(gridID, yvals);

  gridDefLaea(gridID, a, lon0, lat0);

  free(xvals);
  free(yvals);

  return (gridID);
}

static
int scan_pcs_def(char *pcs_def, char proj[128], double *a, double *lon0, double *lat0)
{
  char *pcs[64];
  int npcs = 0;
  int i;
  int len;
  int nfound = 0;

  strcpy(proj, "unknown");
  *a = 1;
  *lon0 = 0;
  *lat0 = 0;

  pcs[npcs++] = &pcs_def[0];
  len = (int)strlen(pcs_def);
  for ( i = 0; i < len; ++i )
    if ( pcs_def[i] == ',' && npcs < 64 )
      {
	pcs_def[i] = 0;
	pcs[npcs++] = &pcs_def[i+1];
      }

  for ( i = 0; i < npcs; ++i )
    {
      if ( memcmp(pcs[i], "proj=", 5) == 0 )
	{
	  pcs[i] += 5;
	  strcpy(proj, pcs[i]);
	  nfound++;
	}
      else if ( memcmp(pcs[i], "a=", 2) == 0 )
	{
	  pcs[i] += 2;
	  *a = atof(pcs[i]);
	  nfound++;
	}
      else if ( memcmp(pcs[i], "lon_0=", 6) == 0 )
	{
	  pcs[i] += 6;
	  *lon0 = atof(pcs[i]);
	  nfound++;
	}
      else if ( memcmp(pcs[i], "lat_0=", 6) == 0 )
	{
	  pcs[i] += 6;
	  *lat0 = atof(pcs[i]);
	  nfound++;
	}
    }

  return (nfound);
}

static
int read_geolocation(hid_t loc_id, int nx, int ny, int lprojtype)
{
  int gridID = -1;
  hid_t grp_id;
  hid_t proj_id, region_id;
  hid_t proj_tid, region_tid;
  hid_t str_tid, fltarr_tid;
  hid_t ptype_id;
  herr_t     status;
  hsize_t dims;
  int xsize, ysize;
  typedef struct proj_t {
    char name[64];
    char ellipsoid[64];
    float  parameter[10];
  } proj_t;
  typedef struct region_t {
    float  xmin;
    float  xmax;
    float  ymin;
    float  ymax;
    float  dx;
    float  dy;
  } region_t;

  proj_t proj;
  region_t region;
  char *projection_name = NULL;

  if ( cdoVerbose ) cdoPrint("Read geolocation:");

  if ( lprojtype )
    {
      ptype_id = H5Topen(loc_id, "ProjType");
      if ( ptype_id >= 0 )
	{
	  projection_name = H5Tget_member_name(ptype_id, 0);
	  H5Tclose(ptype_id);
	}
    }

  str_tid = H5Tcopy(H5T_C_S1);
  H5Tset_size(str_tid, 64);
  dims = 10;
  fltarr_tid = H5Tarray_create(H5T_NATIVE_FLOAT, 1, &dims, NULL);

  proj_tid = H5Tcreate(H5T_COMPOUND, sizeof(proj_t));
  if ( projection_name )
    H5Tinsert(proj_tid, projection_name, HOFFSET(proj_t, name), str_tid);
  else
    H5Tinsert(proj_tid, "Projection name", HOFFSET(proj_t, name), str_tid);
  H5Tinsert(proj_tid, "Reference ellipsoid", HOFFSET(proj_t, ellipsoid), str_tid);
  H5Tinsert(proj_tid, "Projection parameter", HOFFSET(proj_t, parameter), fltarr_tid);

  if ( projection_name ) free(projection_name);

  grp_id = H5Gopen(loc_id, "Geolocation");

  proj_id = H5Dopen(grp_id, "Projection");
  if ( proj_id < 0 )
    proj_id = H5Dopen(grp_id, "projection");
  /*
  {
    hid_t tid;
    int nmem;
    int im;

    tid = H5Dget_type(proj_id);
    nmem = H5Tget_nmembers(tid);
    for ( im = 0; im < nmem; ++im )
      {
	printf("%d %s\n", im, H5Tget_member_name(tid, im));
      }
  }
  */
  if ( proj_id < 0 )
    memset(&proj, 0, sizeof(proj_t));
  else
    status = H5Dread(proj_id, proj_tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, &proj);

  H5Dclose(proj_id);
  H5Tclose(proj_tid);
  H5Tclose(str_tid);
  H5Tclose(fltarr_tid);

  if ( cdoVerbose )
    cdoPrint("  Projection: name=%s\n\t\t\tellipsoid=%s\n\t\t\tparameter=%g %g %g %g %g %g",
	     proj.name, proj.ellipsoid,
	     proj.parameter[0], proj.parameter[1], proj.parameter[2],
	     proj.parameter[3], proj.parameter[4], proj.parameter[5]);

  region_tid = H5Tcreate(H5T_COMPOUND, sizeof(region_t));
  H5Tinsert(region_tid, "xmin", HOFFSET(region_t, xmin), H5T_NATIVE_FLOAT);
  H5Tinsert(region_tid, "xmax", HOFFSET(region_t, xmax), H5T_NATIVE_FLOAT);
  H5Tinsert(region_tid, "ymin", HOFFSET(region_t, ymin), H5T_NATIVE_FLOAT);
  H5Tinsert(region_tid, "ymax", HOFFSET(region_t, ymax), H5T_NATIVE_FLOAT);
  H5Tinsert(region_tid, "dx",   HOFFSET(region_t, dx),   H5T_NATIVE_FLOAT);
  H5Tinsert(region_tid, "dy",   HOFFSET(region_t, dy),   H5T_NATIVE_FLOAT);

  region_id = H5Dopen(grp_id, "Region");
  if ( region_id < 0 )
    region_id = H5Dopen(grp_id, "region");

  if ( region_id < 0 )
    memset(&region, 0, sizeof(region_t));
  else
    status = H5Dread(region_id, region_tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, &region);

  H5Dclose(region_id);
  H5Tclose(region_tid);

  if ( region.xmin > region.xmax )
    {
      double xmin = region.xmin;
      region.xmin = region.xmax;
      region.xmax = xmin;
      if ( cdoVerbose ) cdoPrint("  Swap xmin/xmax");
    }

  if ( cdoVerbose )
    cdoPrint("  Region: xmin=%g xmax=%g ymin=%g ymax=%g dx=%g dy=%g",
	     region.xmin, region.xmax, region.ymin, region.ymax, region.dx, region.dy);

  H5Gclose(grp_id);

  /* check region */
  xsize = (int)lround((region.xmax-region.xmin)/region.dx);
  ysize = (int)lround((region.ymax-region.ymin)/region.dy);

  if ( cdoVerbose ) cdoPrint("  Size: xsize=%d  ysize=%d", xsize, ysize);

  /* some CM-SAF files have incorrect entries for some metadata. */
  /* these are corrected in the following sections. */
  /* in case of questions on this, contact frank.kaspar@dwd.de */
  if ( strcmp(proj.ellipsoid, "WSG-84") == 0 ) strcpy(proj.ellipsoid, "WGS-84");

  if ( (int)region.xmin == -8887500 &&  (int)region.xmax == -8887500 &&
       (int)region.ymin == 8887500 && (int)region.ymax == 8887500 &&
       (int)region.dx == 15000 && (int)region.dy  == 15000 )
    {
       region.xmax =  8887500.0;
       region.ymin = -8887500.0;
       if ( cdoVerbose )
         cdoPrint("  Corrected region: xmin=%g xmax=%g ymin=%g ymax=%g dx=%g dy=%g",
                region.xmin, region.xmax, region.ymin, region.ymax, region.dx, region.dy);

       xsize = (int)lround((region.xmax-region.xmin)/region.dx);
       ysize = (int)lround((region.ymax-region.ymin)/region.dy);
       if ( cdoVerbose ) cdoPrint("  Corrected size: xsize=%d  ysize=%d", xsize, ysize);
    }

  if (nx == 298 && ny == 371 &&
      (int)region.xmin == -6709222 && (int)region.xmax == 6709222 &&
      (int)region.ymin == -6664078 && (int)region.ymax == 9984898 &&
      (int)region.dx == 45000 && (int)region.dy  == 45000 )
    {
        region.xmin = -6705000;
        region.xmax =  6705000;
        region.ymin = -6705000;
        region.ymax =  9990000;
        cdoPrint("  Corrected region: xmin=%g xmax=%g ymin=%g ymax=%g dx=%g dy=%g",
               region.xmin, region.xmax, region.ymin, region.ymax, region.dx, region.dy);

        xsize = (int)lround((region.xmax-region.xmin)/region.dx);
        ysize = (int)lround((region.ymax-region.ymin)/region.dy);
        if ( cdoVerbose ) cdoPrint("  Corrected size: xsize=%d  ysize=%d", xsize, ysize);
    }

  if ( strcmp(proj.name, "sinusoidal") != 0  &&
       ( (nx == xsize && ny == ysize &&
	  (int)region.xmin == -8887500 && (int)region.xmax == 8887500 &&
	  (int)region.ymin == -8887500 && (int)region.ymax == 8887500 &&
	  (int)region.dx == 15000 && (int)region.dy  == 15000 )
	 ||
	 (nx == xsize && ny == ysize &&
	  (int)region.xmin == -5827500 && (int)region.xmax == 5827500 &&
	  (int)region.ymin ==  3307500 && (int)region.ymax == 8887500 &&
	  (int)region.dx == 15000 && (int)region.dy  == 15000 )
	 ||
	 (nx == xsize && ny == ysize &&
	  (int)region.xmin == -5827500 && (int)region.xmax == 5827500 &&
	  (int)region.ymin ==  3307500 && (int)region.ymax == 8887500 &&
	  (int)region.dx == 45000 && (int)region.dy  == 45000 )
	 ||
	 (nx == xsize && ny == ysize &&
	  (int)region.xmin == -5827500 && (int)region.xmax == 5827500 &&
	  (int)region.ymin ==  3307500 && (int)region.ymax == 8887500 &&
	  (int)region.dx == 3000 && (int)region.dy  == 3000 )
	 ||
	 (nx == 298 && ny == 371 &&
	  (int)region.xmin == -6709222 && (int)region.xmax == 6709222 &&
	  (int)region.ymin == -6664078 && (int)region.ymax == 9984898 &&
	  (int)region.dx == 45000 && (int)region.dy  == 45000 )
	 ||
	 (nx == xsize && ny == ysize &&
	  (int)region.xmin == -6705000 && (int)region.xmax == 6705000 &&
	  (int)region.ymin == -6705000 && (int)region.ymax == 9990000 &&
	  (int)region.dx == 45000 && (int)region.dy  == 45000 ) ) )
    {
      if ( cdoVerbose ) cdoPrint("Replacing incorrect projection parameters for sinusoidal products:");
      strcpy(proj.ellipsoid, "WGS-84");
      strcpy(proj.name, "sinusoidal");
      proj.parameter[0] = 0.0;
      proj.parameter[1] = 0.0;
      proj.parameter[2] = 0.0;
      proj.parameter[3] = 0.0;
      proj.parameter[4] = -99.99;
      proj.parameter[5] = -99.99;
      if ( cdoVerbose )
	  cdoPrint("proj1 = %g, proj2 = %g, proj3 = %g, proj4 = %g,",
		   proj.parameter[0], proj.parameter[1],
		   proj.parameter[2], proj.parameter[3]);
    }

  if ( nx == xsize && ny == ysize &&
       strcmp(proj.name, "sinusoidal") == 0 &&
       strcmp(proj.ellipsoid, "WGS-84") == 0  )
    {
      gridID = defSinusoidalGrid(nx, ny, region.xmin, region.xmax, region.ymin, region.ymax,
				 region.dx, region.dy,
				 proj.parameter[0], proj.parameter[1],
				 proj.parameter[2], proj.parameter[3]);
    }
  /* modification by Frank Kaspar */
  else if ( nx == xsize && ny == ysize &&
            strcmp(proj.name, "Lambert Azimuthal Equal Area") == 0 &&
            memcmp(proj.ellipsoid, "Sphere", 6) == 0 )
    {
      double a;
      if ( proj.parameter[4] < 0 )
        {
          a=6370997.0;
        }
      else
        {
          a=proj.parameter[4];
        }
      gridID = defLaeaGrid(nx, ny,
                           region.xmin, region.xmax,
                           region.ymin, region.ymax,
                           region.dx, region.dy,
                           a, proj.parameter[2], proj.parameter[3]);
    }
  else if ( memcmp(proj.name, "Cylindrical Equal Area", 22) == 0 &&
            memcmp(proj.ellipsoid, "Sphere", 6) == 0 )
    {
      double c0  = 0.001*sqrt(proj.parameter[5]);       /* nominal spatial resolution */
      double lts = proj.parameter[3];
      double re  = proj.parameter[4]/1000;  /* Earth radius [km]*/
      if ( cdoVerbose ) cdoPrint("  c0 = %g, lts = %g, re = %g", c0, lts, re);
      gridID = defLonLatGrid(nx, ny, c0, lts, re);
    }
  else if ( nx == 386 && ny == 162 )
    {
      double c0  = 90;       /* nominal spatial resolution */
      double lts = 30;
      double re  = 6371.228;  /* Earth radius [km]*/
      if ( cdoVerbose ) cdoPrint("  c0 = %g, lts = %g, re = %g", c0, lts, re);
      gridID = defLonLatGrid(nx, ny, c0, lts, re);
    }

  return (gridID);
}

static
int read_region(hid_t loc_id, int nx, int ny)
{
  int gridID = -1;
  hid_t grp_id;
  hid_t region_id;
  hid_t region_tid;
  hid_t str64_tid, str128_tid, fltarr_tid;
  herr_t     status;
  hsize_t dims;
  typedef struct region_t {
    double area_extent[4];
    int xsize;
    int ysize;
    float xscale;
    float yscale;
    float lat_0;
    float lon_0;
    float lat_ts;
    char id[128];
    char name[128];
    char pcs_id[128];
    char pcs_def[128];
  } region_t;
  region_t region;
  int nfound;
  char proj[128];
  double a, lon0, lat0;
  double xmin, ymin, xmax, ymax;
  double dx, dy;

  if ( cdoVerbose ) cdoPrint("Read region:");

  /*
   * Create a data type for region
   */
  region_tid = H5Tcreate(H5T_COMPOUND, sizeof(region_t));
  dims = 4;
  fltarr_tid = H5Tarray_create(H5T_NATIVE_DOUBLE, 1, &dims, NULL);
  str64_tid = H5Tcopy(H5T_C_S1);
  H5Tset_size(str64_tid, 128);
  str128_tid = H5Tcopy(H5T_C_S1);
  H5Tset_size(str128_tid, 128);

  H5Tinsert(region_tid, "area_extent", HOFFSET(region_t, area_extent), fltarr_tid);
  H5Tinsert(region_tid, "xsize",   HOFFSET(region_t, xsize),   H5T_NATIVE_INT);
  H5Tinsert(region_tid, "ysize",   HOFFSET(region_t, ysize),   H5T_NATIVE_INT);
  H5Tinsert(region_tid, "xscale",  HOFFSET(region_t, xscale),  H5T_NATIVE_FLOAT);
  H5Tinsert(region_tid, "yscale",  HOFFSET(region_t, yscale),  H5T_NATIVE_FLOAT);
  H5Tinsert(region_tid, "lat_0",   HOFFSET(region_t, lat_0),   H5T_NATIVE_FLOAT);
  H5Tinsert(region_tid, "lon_0",   HOFFSET(region_t, lon_0),   H5T_NATIVE_FLOAT);
  H5Tinsert(region_tid, "lat_ts",  HOFFSET(region_t, lat_ts),  H5T_NATIVE_FLOAT);
  H5Tinsert(region_tid, "id",      HOFFSET(region_t, id),      str64_tid);
  H5Tinsert(region_tid, "name",    HOFFSET(region_t, name),    str64_tid);
  H5Tinsert(region_tid, "pcs_id",  HOFFSET(region_t, pcs_id),  str64_tid);
  H5Tinsert(region_tid, "pcs_def", HOFFSET(region_t, pcs_def), str128_tid);

  grp_id = H5Gopen(loc_id, "/");

  region_id = H5Dopen(grp_id, "region");
  /*
  {
    hid_t tid;
    int nmem;
    int im;

    tid = H5Dget_type(proj_id);
    nmem = H5Tget_nmembers(tid);
    for ( im = 0; im < nmem; ++im )
      {
	printf("%d %s\n", im, H5Tget_member_name(tid, im));
      }
  }
  */
  status = H5Dread(region_id, region_tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, &region);

  if ( cdoVerbose )
    {
      printf("area_extent[0] = %g\n", region.area_extent[0]);
      printf("area_extent[1] = %g\n", region.area_extent[1]);
      printf("area_extent[2] = %g\n", region.area_extent[2]);
      printf("area_extent[3] = %g\n", region.area_extent[3]);
      printf("xsize = %d\n", region.xsize);
      printf("ysize = %d\n", region.ysize);
      printf("xscale = %g\n", region.xscale);
      printf("yscale = %g\n", region.yscale);
      printf("lat_0 = %g\n", region.lat_0);
      printf("lon_0 = %g\n", region.lon_0);
      printf("lat_ts = %g\n", region.lat_ts);
      printf("id = %s\n", region.id);
      printf("name = %s\n", region.name);
      printf("pcs_id = %s\n", region.pcs_id);
      printf("pcs_def = %s\n", region.pcs_def);
    }

  H5Dclose(region_id);
  H5Tclose(region_tid);
  H5Tclose(str64_tid);
  H5Tclose(str128_tid);
  H5Tclose(fltarr_tid);

  H5Gclose(grp_id);

  /* check region */

  nfound = scan_pcs_def(region.pcs_def, proj, &a, &lon0, &lat0);

  if ( cdoVerbose )
    {
      printf("proj = %s\n", proj);
      printf("a    = %g\n", a);
      printf("lon0 = %g\n", lon0);
      printf("lat0 = %g\n", lat0);
    }

  xmin =  region.area_extent[0];
  ymin =  region.area_extent[1];
  xmax =  region.area_extent[2];
  ymax =  region.area_extent[3];

  dx = (xmax-xmin) / nx;
  dy = (ymax-ymin) / ny;
  /*
  xsize = (int)lround((region.xmax-region.xmin)/region.dx);
  ysize = (int)lround((region.ymax-region.ymin)/region.dy);

  if ( cdoVerbose ) cdoPrint("  Size: xsize=%d  ysize=%d", xsize, ysize);
  */

  if ( nfound == 4 &&
       nx == region.xsize && ny == region.ysize &&
       strcmp(proj, "laea") == 0 )
    {
      gridID = defLaeaGrid(nx, ny, xmin, xmax, ymin, ymax,
			   dx, dy, a, lon0, lat0);
    }

  return (gridID);
}

static
void read_dataset(hid_t loc_id, const char *name, void *opdata)
{
  hid_t dset_id, type_id;
  hid_t   dataspace;
  hsize_t dims_out[9];  /* dataset dimensions           */
  herr_t  status;	/* Generic return value		*/
  hid_t   attr, atype, atype_mem;
  hid_t   native_type;
  int iattr;
  float fattr;
  double dattr;
  char attname[CDI_MAX_NAME];
  H5T_class_t  type_class;
  H5T_class_t atype_class;
  size_t atype_size;
  int     rank;
  int nx = 0, ny = 0, nz = 0, nt = 0;
  int gridsize, offset;
  double *array;
  double addoffset = 0, scalefactor = 1, missval = cdiInqMissval();
  int laddoffset = 0, lscalefactor = 0, lmissval = 0;
  int nset;
  int i, k;
  int ftype = 0;
  int len;
  int dtype = DATATYPE_FLT32;
  char attstring[4096];     /* Buffer to read string attribute back */
  char varname[CDI_MAX_NAME];
  short *mask = NULL;
  double minval, maxval;
  int nmiss;
  int num_attrs;

  attstring[0] = 0;
  strcpy(varname, name);

  dset_id = H5Dopen(loc_id, varname);

  type_id = H5Dget_type(dset_id);  /* get datatype*/

  type_class = H5Tget_class(type_id);
  if ( type_class < 0 )
    {
      cdoAbort(" Invalid datatype for %s", varname);
    }
  /*
  else {
    if(type_class == H5T_INTEGER)
      puts("   Datatype is 'H5T_NATIVE_INTEGER'.\n");
    if(type_class == H5T_FLOAT)
      puts("   Datatype is 'H5T_NATIVE_FLOAT'.\n");
    if(type_class == H5T_STRING)
      puts("   Datatype is 'H5T_NATIVE_STRING'.\n");
    if(type_class == H5T_BITFIELD)
      puts("   Datatype is 'H5T_NATIVE_BITFIELD'.\n");
    if(type_class == H5T_OPAQUE)
      puts("   Datatype is 'H5T_NATIVE_OPAQUE'.\n");
    if(type_class == H5T_COMPOUND)
      puts("   Datatype is 'H5T_NATIVE_COMPOUND'.\n");
  }
  */
  native_type = H5Tget_native_type(type_id, H5T_DIR_ASCEND);
  if      ( H5Tequal(native_type, H5T_NATIVE_SCHAR)  > 0 ) {ftype=0; dtype = DATATYPE_INT8;}
  else if ( H5Tequal(native_type, H5T_NATIVE_UCHAR)  > 0 ) {ftype=0; dtype = DATATYPE_UINT8;}
  else if ( H5Tequal(native_type, H5T_NATIVE_SHORT)  > 0 ) {ftype=0; dtype = DATATYPE_INT16;}
  else if ( H5Tequal(native_type, H5T_NATIVE_USHORT) > 0 ) {ftype=0; dtype = DATATYPE_UINT16;}
  else if ( H5Tequal(native_type, H5T_NATIVE_INT)    > 0 ) {ftype=0; dtype = DATATYPE_INT32;}
  else if ( H5Tequal(native_type, H5T_NATIVE_UINT)   > 0 ) {ftype=0; dtype = DATATYPE_UINT32;}
  else if ( H5Tequal(native_type, H5T_NATIVE_FLOAT)  > 0 ) {ftype=1; dtype = DATATYPE_FLT32;}
  else if ( H5Tequal(native_type, H5T_NATIVE_DOUBLE) > 0 ) {ftype=1; dtype = DATATYPE_FLT64;}
  else
    {
      cdoWarning("Dataset %s skipped, unsupported native datatype!", varname);
      goto RETURN;
    }
  H5Tclose(native_type);

  dataspace = H5Dget_space(dset_id);    /* dataspace handle */
  rank      = H5Sget_simple_extent_ndims(dataspace);
  status    = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);

  if ( rank == 2 )
    {
      nx = dims_out[1];
      ny = dims_out[0];
      nz = 1;
      nt = 1;
    }
  else if ( rank == 3 )
    {
      nx = dims_out[2];
      ny = dims_out[1];
      nz = 1;
      nt = dims_out[0];
    }
  else
    {
      cdoWarning("Dataset %s skipped, unsupported rank (=%d)!", varname, rank);
      goto RETURN;
    }


  len = (int) strlen(varname);
  if ( ((datasets_t *) opdata)->mergelevel )
    if ( isdigit(varname[len-1]) && memcmp(varname, "Data", 4) != 0 )
      {
	if ( nt > 1 ) cdoAbort("Combination of nlevel > 1 and ntime > 1 not implemented!");

	nz = atoi(&varname[len-1]);
	varname[len-1] = 0;
      }

  gridsize = nx*ny;

  if ( nz == 1 )
    nset = ((datasets_t *) opdata)->nsets;
  else
    {
      for ( nset = 0; nset < ((datasets_t *) opdata)->nsets; ++nset )
	{
	  if ( strcmp(varname, ((datasets_t *) opdata)->obj[nset].name) == 0 ) break;
	}

      if ( nset >= ((datasets_t *) opdata)->nsets )
	cdoAbort("3D var %s not found!", varname);
    }

  if ( nset < MAX_DSETS )
    {
      if ( cdoVerbose ) print_filter(dset_id, varname);

      num_attrs = H5Aget_num_attrs(dset_id);
      for( i = 0; i < num_attrs; i++ )
	{
	  attr = H5Aopen_idx(dset_id, i);
	  atype = H5Aget_type(attr);
	  H5Aget_name(attr, sizeof(attname), attname);

	  if ( strcmp(attname, "CLASS") == 0 ||
	       strcmp(attname, "IMAGE_VERSION") == 0 ||
	       strcmp(attname, "PALETTE") == 0 ) continue;

	  atype_mem = H5Tget_native_type(atype, H5T_DIR_ASCEND);
	  atype_size = H5Tget_size(atype);
	  atype_class = H5Tget_class(atype);

	  len = strlen(attname);
	  for ( k = 0; k < len; ++k ) attname[k] = tolower(attname[k]);

	  if ( strcmp(attname, "intercept") == 0 ||
	       strcmp(attname, "offset")    == 0 )
	    {
	      if ( atype_class == H5T_FLOAT )
		{
		  if ( atype_size == 4 )
		    {
		      status = H5Aread(attr, H5T_NATIVE_FLOAT, &fattr);
		      if ( status >= 0 )
			{
			  addoffset  = fattr;
			  laddoffset = 1;
			}
		    }
		  else
		    {
		      status = H5Aread(attr, H5T_NATIVE_DOUBLE, &dattr);
		      if ( status >= 0 )
			{
			  addoffset  = dattr;
			  laddoffset = 1;
			}
		    }

		  if ( laddoffset == 0 )
		    cdoWarning("Reading of float attribute %s failed!", attname);
		}
	      else if ( atype_class == H5T_INTEGER )
		{
		  status = H5Aread(attr, H5T_NATIVE_INT, &iattr);
		  if ( status >= 0 )
		    {
		      addoffset  = iattr;
		      laddoffset = 1;
		    }
		  else
		    cdoWarning("Reading of integer attribute %s failed!", attname);
		}
	      else
		cdoWarning("Attribute %s has unsupported data type!", attname);
	    }
	  else if ( strcmp(attname, "gain") == 0 ||
		    strcmp(attname, "scaling_factor") == 0 )
	    {
	      if ( atype_class == H5T_FLOAT )
		{
		  if ( atype_size == 4 )
		    {
		      status = H5Aread(attr, H5T_NATIVE_FLOAT, &fattr);
		      if ( status >= 0 )
			{
			  scalefactor  = fattr;
			  lscalefactor = 1;
			}
		    }
		  else
		    {
		      status = H5Aread(attr, H5T_NATIVE_DOUBLE, &dattr);
		      if ( status >= 0 )
			{
			  scalefactor  = dattr;
			  lscalefactor = 1;
			}
		    }

		  if ( lscalefactor == 0 )
		    cdoWarning("Reading of float attribute %s failed!", attname);
		}
	      else if ( atype_class == H5T_INTEGER )
		{
		  status = H5Aread(attr, H5T_NATIVE_INT, &iattr);
		  if ( status >= 0 )
		    {
		      scalefactor  = iattr;
		      lscalefactor = 1;
		    }
		  else
		    cdoWarning("Reading of integer attribute %s failed!", attname);
		}
	      else
		cdoWarning("Attribute %s has unsupported data type!", attname);
	    }
	  else if ( strncmp(attname, "no_data", 7) == 0 ||
		    strncmp(attname, "nodata", 6)  == 0 )
	    {
	      if ( atype_class == H5T_FLOAT )
		{
		  if ( atype_size == 4 )
		    {
		      status = H5Aread(attr, H5T_NATIVE_FLOAT, &fattr);
		      if ( status >= 0 )
			{
			  missval  = fattr;
			  lmissval = 1;
			}
		    }
		  else
		    {
		      status = H5Aread(attr, H5T_NATIVE_DOUBLE, &dattr);
		      if ( status >= 0 )
			{
			  missval  = dattr;
			  lmissval = 1;
			}
		    }

		  if ( lmissval == 0 )
		    cdoWarning("Reading of float attribute %s failed!", attname);
		}
	      else if ( atype_class == H5T_INTEGER )
		{
		  status = H5Aread(attr, H5T_NATIVE_INT, &iattr);
		  if ( status >= 0 )
		    {
		      missval  = iattr;
		      lmissval = 1;
		    }
		  else
		    cdoWarning("Reading of integer attribute %s failed!", attname);
		}
	      else
		cdoWarning("Attribute %s has unsupported data type!", attname);
	    }
	  else if ( strcmp(attname, "description") == 0 )
	    {
	      H5Aread(attr, atype_mem, attstring);
	      if ( ((datasets_t *) opdata)->obj[nset].description )
		free(((datasets_t *) opdata)->obj[nset].description);
	      ((datasets_t *) opdata)->obj[nset].description = strdup(attstring);
	    }
	  else if ( strcmp(attname, "title") == 0 )
	    {
	      H5Aread(attr, atype_mem, attstring);
	      if ( ((datasets_t *) opdata)->obj[nset].title )
		free(((datasets_t *) opdata)->obj[nset].title);
	      ((datasets_t *) opdata)->obj[nset].title = strdup(attstring);
	    }
	  else if ( strcmp(attname, "time") == 0 )
	    {
	      H5Aread(attr, atype_mem, attstring);
	      if ( ((datasets_t *) opdata)->obj[nset].time )
		free(((datasets_t *) opdata)->obj[nset].time);
	      ((datasets_t *) opdata)->obj[nset].time = strdup(attstring);
	    }
	  else if ( strcmp(attname, "unit") == 0 )
	    {
	      H5Aread(attr, atype_mem, attstring);
	      ((datasets_t *) opdata)->obj[nset].units = strdup(attstring);
	    }

	  H5Tclose(atype_mem);
	  H5Aclose(attr);
	  H5Tclose(atype);
	}

      offset = gridsize*(nz-1);
      array  = ((datasets_t *) opdata)->obj[nset].array;
      array  = (double*) realloc(array, gridsize*nz*nt*sizeof(double));
      ((datasets_t *) opdata)->obj[nset].array    = array;
      array  = array+offset;

      if ( ftype )
	{
	  if ( dtype == DATATYPE_FLT32 )
	    {
	      float *farray;
	      int i;
	      farray = (float*) malloc(gridsize*nt*sizeof(float));
	      status = H5Dread(dset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, farray);
	      if ( status < 0 )
		cdoAbort("Reading of NATIVE_FLOAT variable %s failed!", varname);
	      for ( i = 0; i < gridsize*nt; ++i ) array[i] = farray[i];
	      free(farray);
	    }
	  else
	    {
	      status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, array);
	      if ( status < 0 )
		cdoAbort("Reading of NATIVE_DOUBLE variable %s failed!", varname);
	    }
	}
      else
	{
	  int *iarray, i;
	  iarray = (int*) malloc(gridsize*nt*sizeof(int));
	  status = H5Dread(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, iarray);
	  if ( status < 0 )
	    cdoAbort("Reading of NATIVE_INT variable %s failed!", varname);
	  for ( i = 0; i < gridsize*nt; ++i ) array[i] = iarray[i];
	  free(iarray);
	}

      ((datasets_t *) opdata)->obj[nset].name     = strdup(varname);
      ((datasets_t *) opdata)->obj[nset].nx       = nx;
      ((datasets_t *) opdata)->obj[nset].ny       = ny;
      ((datasets_t *) opdata)->obj[nset].nz       = nz;
      ((datasets_t *) opdata)->obj[nset].nt       = nt;
      ((datasets_t *) opdata)->obj[nset].gridsize = gridsize;

      if ( nz > 1 )
	{
	  if ( ((datasets_t *) opdata)->obj[nset].dtype != dtype )
	    cdoWarning("Data type changes over levels!");

	  if ( laddoffset && !DBL_IS_EQUAL(((datasets_t *) opdata)->obj[nset].offset, addoffset) )
	    cdoWarning("Offset changes over levels!");

	  if ( lscalefactor && !DBL_IS_EQUAL(((datasets_t *) opdata)->obj[nset].scale, scalefactor) )
	    cdoWarning("Scalefactor changes over levels!");

	  if ( lmissval && !DBL_IS_EQUAL(((datasets_t *) opdata)->obj[nset].missval, missval) )
	    cdoWarning("Missing value changes over levels!");
	}

      if ( nz == 1 ) ((datasets_t *) opdata)->nsets++;

      mask = (short*) malloc(gridsize*nt*sizeof(short));
      memset(mask, 0, gridsize*nt*sizeof(short));

      nmiss  = 0;

      minval =  1e35;
      maxval = -1e35;
      for ( i = 0; i < gridsize*nt; i++ )
	{
	  if ( array[i] < minval ) minval = array[i];
	  if ( array[i] > maxval ) maxval = array[i];
	}

      if ( cdoVerbose )
	cdoPrint("Dataset %s: missval = %g  addoffset = %g  scalefactor = %g",
		 varname, missval, addoffset, scalefactor);

      if ( cdoVerbose )
	cdoPrint("Dataset %s: dtype = %d  minval = %g  maxval = %g  missval = %g",
		 varname, dtype,  minval, maxval, missval);

      if ( dtype == DATATYPE_UINT8 )
	{
	  if ( minval >= 0 && maxval <= 127 ) dtype = DATATYPE_INT8;
	}
      else if ( dtype == DATATYPE_UINT16 )
	{
	  if ( minval >= 0 && maxval <= 32767 ) dtype = DATATYPE_INT16;
	}

      laddoffset   = IS_NOT_EQUAL(addoffset,   0);
      lscalefactor = IS_NOT_EQUAL(scalefactor, 1);

      if ( laddoffset || lscalefactor )
	{
	  for ( i = 0; i < gridsize*nt; i++ )
	    if ( !DBL_IS_EQUAL(array[i], missval) )
	      {
		mask[i] = 0;

		if ( lscalefactor ) array[i] *= scalefactor;
		if ( laddoffset   ) array[i] += addoffset;
	      }
	    else
	      {
		nmiss++;
		mask[i] = 1;
	      }
	}

      minval =  1e35;
      maxval = -1e35;
      for ( i = 0; i < gridsize*nt; i++ )
	if ( mask[i] == 0 )
	  {
	    if ( array[i] < minval ) minval = array[i];
	    if ( array[i] > maxval ) maxval = array[i];
	  }

      if ( cdoVerbose )
	cdoPrint("Dataset %s: dtype = %d  minval = %g  maxval = %g  missval = %g",
		 varname, dtype,  minval, maxval, missval);

      if ( nmiss > 0 )
	{
	  if ( ! (missval < minval || missval > maxval) )
	    {
	      if ( DBL_IS_EQUAL(missval, 255.) && dtype == DATATYPE_UINT8 )
		{
		  missval = -255;
		  dtype   = DATATYPE_INT16;
		  cdoPrint("Dataset %s: changed missval to %g and datatype to INT16!",
			   varname, missval);

		  for ( i = 0; i < gridsize*nt; i++ )
		    if ( mask[i] ) array[i] = missval;
		}
	      else
		cdoWarning(" Missing value is inside the range of valid values!\n"
                           " Dataset %s,  Missval: %g,  Range: %g - %g",
			   varname, missval, minval, maxval);
	    }
	}

      ((datasets_t *) opdata)->obj[nset].dtype    = dtype;
      ((datasets_t *) opdata)->obj[nset].loffset  = laddoffset;
      ((datasets_t *) opdata)->obj[nset].lscale   = lscalefactor;
      ((datasets_t *) opdata)->obj[nset].lmissval = lmissval;
      ((datasets_t *) opdata)->obj[nset].offset   = addoffset;
      ((datasets_t *) opdata)->obj[nset].scale    = scalefactor;
      ((datasets_t *) opdata)->obj[nset].missval  = missval;

      free(mask);
      mask = NULL;
    }
  else
    {
      cdoWarning("Too many datasets (MAX = %d)!", MAX_DSETS);
      goto RETURN;
    }

  H5Sclose(dataspace);

 RETURN:

  H5Dclose(dset_id);
  H5Tclose(type_id);
}

static herr_t
obj_info(hid_t loc_id, const char *name, void *opdata)
{
  H5G_obj_t obj_type;
  H5G_stat_t statbuf;

  H5Gget_objinfo(loc_id, name, FALSE, &statbuf);

  obj_type = statbuf.type;

  switch (obj_type) {
  case H5G_GROUP:
    if ( cdoVerbose ) cdoPrint(" Object with name %s is a group", name);
    if ( strcmp(name, "Data") == 0 )
      {
	((datasets_t *) opdata)->mergelevel = TRUE;
	H5Giterate(loc_id, name, NULL, obj_info, opdata);
      }
    else if ( strcmp(name, "Geolocation") == 0 )
      {
	((datasets_t *) opdata)->lgeoloc = TRUE;
      }
    else if ( strcmp(name, "Metadata") == 0 )
      {
	((datasets_t *) opdata)->lmetadata = TRUE;
      }
    break;
  case H5G_DATASET:
    if ( cdoVerbose ) cdoPrint(" Object with name %s is a dataset", name);
    if ( strstr(name, "PALETTE") )
      {
	if ( cdoVerbose ) cdoPrint("   Skip dataset: %s", name);
      }
    /*else if ( strstr(name, "egion") ) */
    else if ( strcmp(name, "region") == 0 )
      {
	((datasets_t *) opdata)->lregion = TRUE;
      }
    else
      {
	if ( cdoVerbose ) cdoPrint("   Read dataset: %s", name);
	read_dataset(loc_id, name, opdata);
      }
    break;
  case H5G_TYPE:
    if ( cdoVerbose ) cdoPrint(" Object with name %s is a named datatype", name);
    if ( strcmp(name, "ProjType") == 0 )
      {
	((datasets_t *) opdata)->lprojtype = TRUE;
      }
    break;
  default:
    cdoAbort(" Unable to identify an object %s", name);
    break;
  }

  return 0;
}

static
void get_global_att(hid_t file_id, const char *obj_path, int vlistID)
{
  hid_t   attr, atype, atype_mem, obj_id, grp_id = -1;
  char attname[CDI_MAX_NAME];
  H5T_class_t  type_class;
  int attint;
  double attflt;
  int i, pos;
  int num_attrs;
  char attstring[4096];     /* Buffer to read string attribute back */

  attstring[0] = 0;

  obj_id = H5Gopen(file_id, obj_path);

  num_attrs = H5Aget_num_attrs(obj_id);

  for( i = 0; i < num_attrs; i++ )
    {
      attr = H5Aopen_idx(obj_id, i);
      atype = H5Aget_type(attr);
      H5Aget_name(attr, sizeof(attname), attname);

      /* remove illegal characters */
      for ( pos = 0; pos < (int)strlen(attname); ++pos )
	if ( attname[pos] == '&' ) attname[pos] = '_';

      atype_mem = H5Tget_native_type(atype, H5T_DIR_ASCEND);
      type_class = H5Tget_class(atype);
      if ( type_class == H5T_STRING )
	{
	  H5Aread(attr, atype_mem, attstring);
	  vlistDefAttTxt(vlistID, CDI_GLOBAL, attname, (int)strlen(attstring), attstring);
	}
      else if ( type_class == H5T_INTEGER )
	{
	  H5Aread(attr, H5T_NATIVE_INT, &attint);
	  vlistDefAttInt(vlistID, CDI_GLOBAL, attname, DATATYPE_INT32, 1, &attint);
	}
      else if ( type_class == H5T_FLOAT )
	{
	  H5Aread(attr, H5T_NATIVE_DOUBLE, &attflt);
	  vlistDefAttFlt(vlistID, CDI_GLOBAL, attname, DATATYPE_FLT64, 1, &attflt);
	}
      H5Tclose(atype_mem);
      H5Aclose(attr);
      H5Tclose(atype);
    }

  if ( grp_id >= 0 ) H5Gclose(grp_id);

}

static
int get_vdate(int vlistID)
{
  int vdate = 0;
  int natts;
  int i, len, type;
  char name[CDI_MAX_NAME];
  char attstr[CDI_MAX_NAME];

  vlistInqNatts(vlistID, CDI_GLOBAL, &natts);

  for ( i = 0; i < natts; ++i )
    {
      vlistInqAtt(vlistID, CDI_GLOBAL, i, name, &type, &len);
      if ( type == DATATYPE_TXT )
	{
	  if ( strcmp(name, "DateAndTime") == 0 ||
	       strcmp(name, "Date_Time") == 0 )
	    {
	      vlistInqAttTxt(vlistID, CDI_GLOBAL, name, CDI_MAX_NAME, attstr);
	      if ( len > 8 ) len = 8;
	      attstr[len] = 0;
	      vdate = atoi(attstr);
	      if ( vdate < 999999 ) vdate = vdate*100 + 1;
	    }
	}
    }

  return (vdate);
}

static
void dsets_init(datasets_t *dsets)
{
  int i;

  dsets->nsets      = 0;
  dsets->mergelevel = 0;
  dsets->lgeoloc    = 0;
  dsets->lregion    = 0;
  dsets->lprojtype  = 0;
  dsets->lmetadata  = 0;

  for ( i = 0; i < MAX_DSETS; ++i )
    {
      dsets->obj[i].nx          = 0;
      dsets->obj[i].ny          = 0;
      dsets->obj[i].nz          = 0;
      dsets->obj[i].name        = NULL;
      dsets->obj[i].description = NULL;
      dsets->obj[i].units       = NULL;
      dsets->obj[i].title       = NULL;
      dsets->obj[i].time        = NULL;
      dsets->obj[i].dtype       = cdoDefaultDataType;
      dsets->obj[i].lscale      = 0;
      dsets->obj[i].loffset     = 0;
      dsets->obj[i].lmissval    = 0;
      dsets->obj[i].missval     = cdiInqMissval();
      dsets->obj[i].array       = NULL;
    }
}
#endif

void *Importcmsaf(void *argument)
{
#if defined(HAVE_LIBHDF5)
  int streamID;
  int gridID = -1, zaxisID, taxisID, vlistID;
  int i, offset;
  int nmiss;
  int ivar;
  int varID, levelID, tsID;
  int nx, ny, nz, nt, gridsize;
  double *array;
  double missval, minval, maxval;
  hid_t	  file_id;	/* HDF5 File ID	        	*/
  herr_t  status;	/* Generic return value		*/
  datasets_t dsets;
  int vdate, vtime;
  int *vtimes = NULL;
#endif

  cdoInitialize(argument);

  if ( cdoDefaultFileType == CDI_UNDEFID )
    cdoDefaultFileType = FILETYPE_NC;

#if defined(HAVE_LIBHDF5)
  dsets_init(&dsets);

  /* Open an existing file. */
  file_id = H5Fopen(cdoStreamName(0)->args, H5F_ACC_RDONLY, H5P_DEFAULT);
  if ( file_id < 0 ) cdoAbort("H5Fopen failed on %s", cdoStreamName(0)->args);

  /* cmsaf_type = get_cmsaf_type(file_id); */

  H5Giterate(file_id, "/", NULL, obj_info, (void *) &dsets);

  if ( dsets.nsets == 0 ) cdoAbort("No dataset found!");

  gridsize = dsets.obj[0].gridsize;
  nx = dsets.obj[0].nx;
  ny = dsets.obj[0].ny;
  nz = dsets.obj[0].nz;
  nt = dsets.obj[0].nt;

  for ( ivar = 0; ivar < dsets.nsets; ++ivar )
    if ( dsets.obj[ivar].nt > 1 )
      {
	nt = dsets.obj[ivar].nt;
	break;
      }

  if ( nt > 1 )
    {
      vtimes = (int*) malloc(nt*sizeof(int));

      for ( i = 0; i < nt; ++i ) vtimes[i] = i*10000 + 45*100;

      if ( dsets.obj[ivar].time )
	{
	  long itime;
	  char *pline = dsets.obj[ivar].time;

	  for ( i = 0; i < nt; ++i )
	    {
	      itime = ((int)strtol(pline, &pline, 10))*100;
	      if ( itime < 0 || itime > 240000 )
		{
		  cdoWarning("Wrong time string!");
		  break;
		}
	      vtimes[i] = (int) itime;
	    }
	}
    }

  if ( cdoVerbose )
    for ( ivar = 0; ivar < dsets.nsets; ++ivar )
      cdoPrint(" Var %d %-20s %dx%d nlev = %d nts = %d",
	       ivar, dsets.obj[ivar].name, nx, ny, nz, dsets.obj[ivar].nt);

  for ( ivar = 1; ivar < dsets.nsets; ++ivar )
    {
      if ( nx != dsets.obj[0].nx || ny != dsets.obj[0].ny )
	cdoAbort("Gridsize must not change!");
      if ( nz != dsets.obj[0].nz )
	cdoAbort("Number of levels must not change!");
    }

  if ( dsets.lgeoloc )
    {
      gridID = read_geolocation(file_id, nx, ny, dsets.lprojtype);
    }
  else if ( dsets.lregion )
    {
      gridID = read_region(file_id, nx, ny);
    }

  if ( gridID == -1 )
    {
      gridID = gridCreate(GRID_GENERIC, gridsize);
      gridDefXsize(gridID, nx);
      gridDefYsize(gridID, ny);
    }

  if ( nz == 1 )
    zaxisID = zaxisCreate(ZAXIS_SURFACE, 1);
  else
    {
      double *levels;
      levels = (double*) malloc(nz*sizeof(double));
      for ( i = 0; i < nz; ++i ) levels[i] = i+1;
      zaxisID = zaxisCreate(ZAXIS_GENERIC, nz);
      zaxisDefLevels(zaxisID, levels);
      free(levels);
    }

  vlistID = vlistCreate();

  if ( nt > 1 )
    taxisID = taxisCreate(TAXIS_RELATIVE);
  else
    taxisID = taxisCreate(TAXIS_ABSOLUTE);

  taxisDefCalendar(taxisID, CALENDAR_STANDARD);

  vlistDefTaxis(vlistID, taxisID);

  for ( ivar = 0; ivar < dsets.nsets; ++ivar )
    {
      if ( dsets.obj[ivar].nt > 1 )
	varID = vlistDefVar(vlistID, gridID, zaxisID, TSTEP_INSTANT);
      else
	varID = vlistDefVar(vlistID, gridID, zaxisID, TSTEP_INSTANT);

      vlistDefVarName(vlistID, varID,  dsets.obj[ivar].name);
      if ( dsets.obj[ivar].description )
	vlistDefVarLongname(vlistID, varID,  dsets.obj[ivar].description);
      if ( dsets.obj[ivar].units )
	vlistDefVarUnits(vlistID, varID,  dsets.obj[ivar].units);
      if ( dsets.obj[ivar].title )
	vlistDefAttTxt(vlistID, varID, "title", (int)strlen(dsets.obj[ivar].title),
		       dsets.obj[ivar].title);

      /*
      vlistDefVarUnits(vlistID, varID, units[i]);
      */
      vlistDefVarDatatype(vlistID, varID, dsets.obj[ivar].dtype);
      if ( dsets.obj[ivar].lmissval )
	vlistDefVarMissval(vlistID, varID, dsets.obj[ivar].missval);
      if ( dsets.obj[ivar].lscale )
	vlistDefVarScalefactor(vlistID, varID, dsets.obj[ivar].scale);
      if ( dsets.obj[ivar].loffset )
	vlistDefVarAddoffset(vlistID, varID, dsets.obj[ivar].offset);
    }

  get_global_att(file_id, "/", vlistID);
  if ( dsets.lmetadata ) get_global_att(file_id, "Metadata", vlistID);

  vdate = get_vdate(vlistID);
  if ( vdate == 0 ) vdate = 10101;

  streamID = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID, vlistID);

  for ( tsID = 0; tsID < nt; ++tsID )
    {
      taxisDefVdate(taxisID, vdate);
      vtime = 0;
      if ( vtimes ) vtime = vtimes[tsID];
      taxisDefVtime(taxisID, vtime);
      streamDefTimestep(streamID, tsID);

      for ( ivar = 0; ivar < dsets.nsets; ++ivar )
	{
	  varID   = ivar;

	  if ( tsID > 0 && dsets.obj[ivar].nt == 1 ) continue;

	  gridsize = dsets.obj[ivar].gridsize;
	  missval  = dsets.obj[ivar].missval;

	  for ( levelID = 0; levelID < nz; ++levelID )
	    {
	      offset = gridsize*levelID;
	      if ( nz == 1 ) offset = gridsize*tsID;
	      array  = dsets.obj[ivar].array+offset;

	      nmiss  = 0;
	      minval =  1e35;
	      maxval = -1e35;

	      for ( i = 0; i < gridsize; i++ )
		{
		  if ( !DBL_IS_EQUAL(array[i], missval) )
		    {
		      if ( array[i] < minval ) minval = array[i];
		      if ( array[i] > maxval ) maxval = array[i];
		    }
		}

	      for ( i = 0; i < gridsize; i++ )
		if ( DBL_IS_EQUAL(array[i], missval) ) nmiss++;

	      if ( cdoVerbose )
		cdoPrint(" Write var %d,  level %d, nmiss %d, missval %g, minval %g, maxval %g",
			 varID, levelID, nmiss, missval, minval, maxval);
	      /*
		if ( ! (missval < minval || missval > maxval) )
		cdoWarning(" Missval is inside of valid values! Name: %s  Range: %g - %g  Missval: %g\n",
		dsets.obj[ivar].name, minval, maxval, missval);
	      */
	      streamDefRecord(streamID,  varID,  levelID);
	      streamWriteRecord(streamID, array, nmiss);
	    }
	}
    }

  /* Close file */
  status = H5Fclose(file_id);

  processDefVarNum(vlistNvars(vlistID), streamID);

  streamClose(streamID);

  vlistDestroy(vlistID);
  gridDestroy(gridID);
  zaxisDestroy(zaxisID);
  taxisDestroy(taxisID);

  for ( ivar = 0; ivar < dsets.nsets; ++ivar )
    free(dsets.obj[ivar].array);

  if ( vtimes ) free(vtimes);

  cdoFinish();
#else
  cdoAbort("HDF5 support not compiled in!");
#endif

  return (0);
}
