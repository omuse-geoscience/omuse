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
#include "pstream.h"
#include "interpol.h"
#include "grid.h"
#include "remap.h"

#if defined(HAVE_LIBYAC)
#include "points.h"
#include "grid_reg2d.h"
#include "grid_search.h"
#include "bucket_search.h"
#include "search.h"
#include "clipping.h"
#include "area.h"
#endif

int timer_yar_remap, timer_yar_remap_init, timer_yar_remap_sort, timer_yar_remap_con, timer_yar_remap_bil;


void yar_remap(double *restrict dst_array, double missval, long dst_size, long num_links, double *restrict map_wts, 
	       long num_wts, const int *restrict dst_add, const int *restrict src_add, const double *restrict src_array)
{
  long n;

  for ( n = 0; n < dst_size; ++n ) dst_array[n] = missval;

  for ( n = 0; n < num_links; ++n ) dst_array[dst_add[n]] = 0;

  for ( n = 0; n < num_links; ++n )
    {
      
      //printf("%5d %5d %5d %g # dst_add src_add n\n", dst_add[n], src_add[n], n, map_wts[n]);
      
      //dst_array[dst_add[n]] += src_array[src_add[n]]*map_wts[num_wts*n];
      dst_array[dst_add[n]] += src_array[src_add[n]]*map_wts[n];
    }
  //for ( n = 0; n < 10; ++n ) printf("array1 %d %g\n", n, src_array[n]);
  //for ( n = 0; n < 10; ++n ) printf("array2 %d %g\n", n, dst_array[n]);
}


void yar_store_link_cnsrv(remapvars_t *rv, long add1, long add2, double weight)
{
  /*
    Input variables:
    int  add1         ! address on grid1
    int  add2         ! address on grid2
    double weight     ! remapping weight for this link
  */
  /* Local variables */
  long nlink; /* link index */

  /*  If the weight is ZERO, do not bother storing the link */

  if ( IS_EQUAL(weight, 0) ) return;

  /*
     If the link does not yet exist, increment number of links and 
     check to see if remap arrays need to be increased to accomodate 
     the new link. Then store the link.
  */
  nlink = rv->num_links;

  rv->num_links++;
  if ( rv->num_links >= rv->max_links )
    resize_remap_vars(rv, rv->resize_increment);

  rv->src_cell_add[nlink] = add1;
  rv->tgt_cell_add[nlink] = add2;

  rv->wts[nlink] = weight;
}

static
void gen_xbounds(int nx, double *xvals, double *xbounds)
{
  int i;

  for ( i = 0; i < nx-1; i++ )
    {
      xbounds[i+1]   = 0.5*(xvals[i] + xvals[i+1]);
    }

  xbounds[0]  = 2*xvals[0] - xbounds[1];
  xbounds[nx] = 2*xvals[nx-1] - xbounds[nx-1];
}

static
void gen_ybounds(int ny, double *yvals, double *ybounds)
{
  int i;

  for ( i = 0; i < ny-1; i++ )
    {
      ybounds[i+1]   = 0.5*(yvals[i] + yvals[i+1]);
    }

  ybounds[0]  = 2*yvals[0] - ybounds[1];
  ybounds[ny] = 2*yvals[ny-1] - ybounds[ny-1];

  if ( yvals[0] > yvals[ny-1] )
    {
      if ( ybounds[0]  >  88 ) ybounds[0]  =  90;
      if ( ybounds[ny] < -88 ) ybounds[ny] = -90;
    }
  else
    {
      if ( ybounds[0]  < -88 ) ybounds[0]  = -90;
      if ( ybounds[ny] >  88 ) ybounds[ny] =  90;
    }
}


void set_source_data(double * source_data, double init_value,
                     unsigned size_x, unsigned size_y) {

   for (unsigned i = 0; i < size_x; ++i)
      for (unsigned j = 0; j < size_y; ++j)
         source_data[i + j * size_x] = init_value;
}


/*
  This routine stores the address and weight for four links associated with one destination
  point in the appropriate address and weight arrays and resizes those arrays if necessary.
*/
static
void store_link_bilin(remapvars_t *rv, int dst_add, int src_add[4], double weights[4])
{
  /*
    Input variables:
    int dst_add       ! address on destination grid
    int src_add[4]    ! addresses on source grid
    double weights[4] ! array of remapping weights for these links
  */
  /* link index */
  long nlink = rv->num_links;
  /*
     Increment number of links and check to see if remap arrays need
     to be increased to accomodate the new link. Then store the link.
  */
  rv->num_links += 4;

  if ( rv->num_links >= rv->max_links ) 
    resize_remap_vars(rv, rv->resize_increment);

  for ( long n = 0; n < 4; ++n )
    {
      rv->src_cell_add[nlink+n] = src_add[n];
      rv->tgt_cell_add[nlink+n] = dst_add;
      rv->wts         [nlink+n] = weights[n];
    }

} /* store_link_bilin */

void yar_remap_bil(field_t *field1, field_t *field2)
{
  int nlonIn, nlatIn;
  int nlonOut, nlatOut;
  int ilat, ilon;
  int gridIDin, gridIDout;
  int i, nmiss;
  int gridsize1, gridsize2;
  double *lonIn, *latIn;
  double *lonOut, *latOut;
  double *xlonIn, *xlatIn;
  double *xlonOut, *xlatOut;
  double **fieldIn;
  double **field;
  double *array = NULL;
  double *array1, *array2;
  double missval;
  int testit = 1;
  double dxIn, dxOut;
  remap_t remap;
  /* static int index = 0; */

  gridIDin  = field1->grid;
  gridIDout = field2->grid;
  array1    = field1->ptr;
  array2    = field2->ptr;
  missval   = field1->missval;

  if ( ! (gridInqXvals(gridIDin, NULL) && gridInqYvals(gridIDin, NULL)) )
    cdoAbort("Source grid has no values");

  remap.src_grid.num_srch_bins = 0;
  remap.tgt_grid.num_srch_bins = 0;
  remap.vars.pinit = FALSE;

  if ( cdoTimer ) timer_start(timer_yar_remap_init);
  remap_grids_init(MAP_TYPE_BILINEAR, 0, gridIDin, &remap.src_grid, gridIDout, &remap.tgt_grid);
  remap_vars_init(MAP_TYPE_BILINEAR, remap.src_grid.size, remap.tgt_grid.size, &remap.vars);
  if ( cdoTimer ) timer_stop(timer_yar_remap_init);


  if ( cdoTimer ) timer_start(timer_yar_remap_init);
  nlonIn = gridInqXsize(gridIDin);
  nlatIn = gridInqYsize(gridIDin);
  gridsize1 = gridInqSize(gridIDin);
  lonIn = (double*) malloc(nlonIn*sizeof(double));
  latIn = (double*) malloc(nlatIn*sizeof(double));
  gridInqXvals(gridIDin, lonIn);
  gridInqYvals(gridIDin, latIn);
  for ( int i = 0; i < nlonIn; ++i ) lonIn[i] *= DEG2RAD;
  for ( int i = 0; i < nlatIn; ++i ) latIn[i] *= DEG2RAD;

  xlonIn = (double*) malloc((nlonIn+1)*sizeof(double));
  xlatIn = (double*) malloc((nlatIn+1)*sizeof(double));
  gridInqXvals(gridIDin, xlonIn);
  gridInqYvals(gridIDin, xlatIn);
  dxIn = xlonIn[1] - xlonIn[0];
  for ( int i = 0; i < nlonIn; ++i ) xlonIn[i] -= dxIn/2;
  for ( int i = 0; i < nlatIn; ++i ) xlatIn[i] -= dxIn/2;
  for ( int i = 0; i < nlonIn; ++i ) xlonIn[i] *= DEG2RAD;
  for ( int i = 0; i < nlatIn; ++i ) xlatIn[i] *= DEG2RAD;
  xlonIn[nlonIn] = xlonIn[nlonIn-1]+dxIn*DEG2RAD;
  xlatIn[nlatIn] = -xlatIn[0];
 
  if ( ! (gridInqXvals(gridIDout, NULL) && gridInqYvals(gridIDout, NULL)) )
    cdoAbort("Target grid has no values");

  nlonOut = gridInqXsize(gridIDout);
  nlatOut = gridInqYsize(gridIDout);
  gridsize2 = gridInqSize(gridIDout);
  lonOut = (double*) malloc(nlonOut*sizeof(double));
  latOut = (double*) malloc(nlatOut*sizeof(double));
  gridInqXvals(gridIDout, lonOut);
  gridInqYvals(gridIDout, latOut);
  for ( int i = 0; i < nlonOut; ++i ) lonOut[i] *= DEG2RAD;
  for ( int i = 0; i < nlatOut; ++i ) latOut[i] *= DEG2RAD;

  xlonOut = (double*) malloc((nlonOut+1)*sizeof(double));
  xlatOut = (double*) malloc((nlatOut+1)*sizeof(double));
  gridInqXvals(gridIDout, xlonOut);
  gridInqYvals(gridIDout, xlatOut);
  dxOut = xlonOut[1] - xlonOut[0];
  for ( int i = 0; i < nlonOut; ++i ) xlonOut[i] -= dxOut/2;
  for ( int i = 0; i < nlatOut; ++i ) xlatOut[i] -= dxOut/2;
  for ( int i = 0; i < nlonOut; ++i ) xlonOut[i] *= DEG2RAD;
  for ( int i = 0; i < nlatOut; ++i ) xlatOut[i] *= DEG2RAD;
  xlonOut[nlonOut] = xlonOut[nlonOut-1]+dxOut*DEG2RAD;
  xlatOut[nlatOut] = xlatOut[0];

  printf("source grid: inc = %g  nlon = %d  nlat = %d\n", dxIn, nlonIn, nlatIn);
  printf("target grid: inc = %g  nlon = %d  nlat = %d\n", dxOut, nlonOut, nlatOut);

  printf("lonIn: %g %g %g ... %g %g\n", lonIn[0]/DEG2RAD, lonIn[1]/DEG2RAD, lonIn[2]/DEG2RAD, lonIn[nlonIn-2]/DEG2RAD, lonIn[nlonIn-1]/DEG2RAD);
  printf("latIn: %g %g %g ... %g %g\n", latIn[0]/DEG2RAD, latIn[1]/DEG2RAD, latIn[2]/DEG2RAD, latIn[nlatIn-2]/DEG2RAD, latIn[nlatIn-1]/DEG2RAD);

  printf("lonOut: %g %g %g ... %g %g\n", lonOut[0]/DEG2RAD, lonOut[1]/DEG2RAD, lonOut[2]/DEG2RAD, lonOut[nlonOut-2]/DEG2RAD, lonOut[nlonOut-1]/DEG2RAD);
  printf("latOut: %g %g %g ... %g %g\n", latOut[0]/DEG2RAD, latOut[1]/DEG2RAD, latOut[2]/DEG2RAD, latOut[nlatOut-2]/DEG2RAD, latOut[nlatOut-1]/DEG2RAD);

#if defined(HAVE_LIBYAC)

  //--------------------------------------------
  // define a grid
  //--------------------------------------------
  unsigned num_source_cells[2];
  unsigned num_target_cells[2];
  num_source_cells[0] = nlonIn;
  num_source_cells[1] = nlatIn;
  num_target_cells[0] = nlonOut;
  num_target_cells[1] = nlatOut;

  unsigned cyclic[2] = {0,0};
  struct grid *source_grid;

  source_grid = reg2d_grid_new(xlonIn, xlatIn, num_source_cells, cyclic);

  if ( cdoTimer ) timer_stop(timer_yar_remap_init);

  struct points source_points;

  //--------------------------------------------
  // define points
  //--------------------------------------------
  init_points(&source_points, source_grid, CELL, lonIn, latIn);

  //--------------------------------------------
  // initialise interpolation
  //--------------------------------------------

  struct dep_list deps;
  unsigned search_id;
  //struct interpolation interpolation;

  printf("src num_grid_corners %d\n", get_num_grid_corners(get_point_grid(&source_points)));

  // search_id = search_init(get_point_grid(&source_points));
  struct grid_search *search;
  unsigned const * curr_src_corners;

  if ( cdoTimer ) timer_start(timer_yar_remap_bil);

  search = bucket_search_new(source_grid);
  /*
  search = sphere_part_search_new(source_grid);
  printf("search %p\n", search);
  printf("search->vtable %p\n", search->vtable);
  printf("search->vtable->do_point_search_p3 %d\n", search->vtable->do_point_search_p3);
  */
  /*

  printf("total_num_dependencies: %d\n", get_total_num_dependencies(deps));

  for ( int i = 0; i < 10; ++i )
    {
      printf("num_deps_per_element %d %d\n", i, deps.num_deps_per_element[i]);
      curr_src_corners = get_dependencies_of_element(deps,i);
      for ( int k = 0; k < deps.num_deps_per_element[i]; ++k )
	printf("  curr_src_corners: %d %d\n", k, curr_src_corners[k]);
    }

  for ( int i = 0; i < 10; ++i )
    {
      double lon = lonOut[i];
      double lat = latOut[i];
      printf("num_deps_per_element %d %d\n", i, deps.num_deps_per_element[i]);
      curr_src_corners = get_dependencies_of_element(deps,i);
      for ( int k = 0; k < deps.num_deps_per_element[i]; ++k )
	printf("  curr_src_corners: %d %d\n", k, curr_src_corners[k]);
    }
  */
  for ( int i = 0; i < gridsize2; ++i )
    {
      double wgts[4];
      double plon;
      double plat;
      double iguess = 0.5, jguess = 0.5;         /*  current guess for bilinear coordinate  */
      int ix, iy;
      int dst_add = i;

      iy = dst_add/nlonOut;
      ix = dst_add - iy*nlonOut;

      plon = lonOut[ix];
      plat = latOut[iy];

      // do search
      do_point_search_p3(search, &plon, &plat, 1, &deps, &source_points);

      if ( i < 10 ) printf("total_num_dependencies: %d\n", get_total_num_dependencies(deps));
      //printf("i, ix, iy %d %d %d\n", i, ix, iy);
      curr_src_corners = get_dependencies_of_element(deps,0);
      for ( int k = 0; k < deps.num_deps_per_element[i]; ++k )
	{
	}

      if ( deps.num_deps_per_element[0] == 4 )
	{
	  int sa, src_add[4];                /*  address for the four source points     */
	  double src_lats[4];            /* latitudes  of the four corner points   */
	  double src_lons[4];            /* longitudes of the four corner points   */

	  for ( int k = 0; k < deps.num_deps_per_element[0]; ++k )
	    {
	      sa = curr_src_corners[k];
	      src_add[k] = sa;
	      iy = sa/nlonIn;
	      ix = sa - iy*nlonIn;
	      src_lats[k] = latIn[iy];
	      src_lons[k] = lonIn[ix];
	      if ( i < 3 ) printf("i, k, iy, ix %d %d %d %d\n", i, k, iy, ix);
	    }

	  // try it with do_point_search_p3
	  if ( find_ij_weights(plon, plat, src_lats, src_lons, &iguess, &jguess) )
	    {

	      wgts[0] = (1.-iguess)*(1.-jguess);
	      wgts[1] = iguess*(1.-jguess);
	      wgts[2] = iguess*jguess;
	      wgts[3] = (1.-iguess)*jguess;

#if defined(_OPENMP)
#pragma omp critical
#endif
	      store_link_bilin(&remap.vars, dst_add, src_add, wgts);
	    }
	}
    }

  if ( cdoTimer ) timer_stop(timer_yar_remap_bil);

  if ( cdoTimer ) timer_start(timer_yar_remap);
  yar_remap(array2, missval, gridInqSize(gridIDout), remap.vars.num_links, remap.vars.wts,
	    remap.vars.num_wts, remap.vars.tgt_cell_add, remap.vars.src_cell_add, array1);
  if ( cdoTimer ) timer_stop(timer_yar_remap);

  nmiss = 0;
  for ( int i = 0; i < gridInqSize(gridIDout); ++i )
    if ( DBL_IS_EQUAL(array2[i], missval) ) nmiss++;

  field2->nmiss = nmiss;

  // if (array) free(array);
  //free(lonIn);
  //free(latIn);
  //free(lonOut);
  //free(latOut);
  //free(fieldIn);
#endif
}


void yar_remap_con(field_t *field1, field_t *field2)
{
  int nlonIn, nlatIn;
  int nlonOut, nlatOut;
  int ilat, ilon;
  int gridIDin, gridIDout;
  int i, nmiss;
  int gridsize1, gridsize2;
  double *lonIn, *latIn;
  double *lonOut, *latOut;
  double *xlonIn, *xlatIn;
  double *xlonOut, *xlatOut;
  double **fieldIn;
  double **field;
  double *array1, *array2;
  double missval;
  double dxIn, dxOut;
  remap_t remap;
  /* static int index = 0; */

  gridIDin  = field1->grid;
  gridIDout = field2->grid;
  array1    = field1->ptr;
  array2    = field2->ptr;
  missval   = field1->missval;

  if ( ! (gridInqXvals(gridIDin, NULL) && gridInqYvals(gridIDin, NULL)) )
    cdoAbort("Source grid has no values");

  remap.src_grid.num_srch_bins = 0;
  remap.tgt_grid.num_srch_bins = 0;
  remap.vars.pinit = FALSE;

  if ( cdoTimer ) timer_start(timer_yar_remap_init);
  remap_grids_init(MAP_TYPE_CONSERV, 0, gridIDin, &remap.src_grid, gridIDout, &remap.tgt_grid);
  remap_vars_init(MAP_TYPE_CONSERV, remap.src_grid.size, remap.tgt_grid.size, &remap.vars);
  if ( cdoTimer ) timer_stop(timer_yar_remap_init);


  if ( cdoTimer ) timer_start(timer_yar_remap_init);
  nlonIn = gridInqXsize(gridIDin);
  nlatIn = gridInqYsize(gridIDin);
  gridsize1 = gridInqSize(gridIDin);
  lonIn = (double*) malloc((nlonIn+1)*sizeof(double));
  latIn = (double*) malloc((nlatIn+1)*sizeof(double));
  gridInqXvals(gridIDin, lonIn);
  gridInqYvals(gridIDin, latIn);
  xlonIn = (double*) malloc((nlonIn)*sizeof(double));
  xlatIn = (double*) malloc((nlatIn)*sizeof(double));
  gridInqXvals(gridIDin, xlonIn);
  gridInqYvals(gridIDin, xlatIn);
  dxIn = lonIn[1] - lonIn[0];
  for ( int i = 0; i < nlonIn; ++i ) lonIn[i] -= dxIn/2;
  for ( int i = 0; i < nlatIn; ++i ) latIn[i] -= dxIn/2;
  for ( int i = 0; i < nlonIn; ++i ) lonIn[i] *= DEG2RAD;
  for ( int i = 0; i < nlatIn; ++i ) latIn[i] *= DEG2RAD;
  lonIn[nlonIn] = lonIn[nlonIn-1]+dxIn*DEG2RAD;
  latIn[nlatIn] = -latIn[0];

  if ( ! (gridInqXvals(gridIDout, NULL) && gridInqYvals(gridIDout, NULL)) )
    cdoAbort("Target grid has no values");

  nlonOut = gridInqXsize(gridIDout);
  nlatOut = gridInqYsize(gridIDout);
  gridsize2 = gridInqSize(gridIDout);
  lonOut = (double*) malloc((nlonOut+1)*sizeof(double));
  latOut = (double*) malloc((nlatOut+1)*sizeof(double));
  gridInqXvals(gridIDout, lonOut);
  gridInqYvals(gridIDout, latOut);
  xlonOut = (double*) malloc((nlonOut+1)*sizeof(double));
  xlatOut = (double*) malloc((nlatOut+1)*sizeof(double));
  gridInqXvals(gridIDout, xlonOut);
  gridInqYvals(gridIDout, xlatOut);
  dxOut = lonOut[1] - lonOut[0];
  for ( int i = 0; i < nlonOut; ++i ) lonOut[i] -= dxOut/2;
  for ( int i = 0; i < nlatOut; ++i ) latOut[i] -= dxOut/2;
  for ( int i = 0; i < nlonOut; ++i ) lonOut[i] *= DEG2RAD;
  for ( int i = 0; i < nlatOut; ++i ) latOut[i] *= DEG2RAD;
  lonOut[nlonOut] = lonOut[nlonOut-1]+dxOut*DEG2RAD;
  latOut[nlatOut] = latOut[0];

  printf("source grid: inc = %g  nlon = %d  nlat = %d\n", dxIn, nlonIn, nlatIn);
  printf("target grid: inc = %g  nlon = %d  nlat = %d\n", dxOut, nlonOut, nlatOut);

  printf("lonIn: %g %g %g ... %g %g\n", lonIn[0]/DEG2RAD, lonIn[1]/DEG2RAD, lonIn[2]/DEG2RAD, lonIn[nlonIn-2]/DEG2RAD, lonIn[nlonIn-1]/DEG2RAD);
  printf("latIn: %g %g %g ... %g %g\n", latIn[0]/DEG2RAD, latIn[1]/DEG2RAD, latIn[2]/DEG2RAD, latIn[nlatIn-2]/DEG2RAD, latIn[nlatIn-1]/DEG2RAD);

  printf("lonOut: %g %g %g ... %g %g\n", lonOut[0]/DEG2RAD, lonOut[1]/DEG2RAD, lonOut[2]/DEG2RAD, lonOut[nlonOut-2]/DEG2RAD, lonOut[nlonOut-1]/DEG2RAD);
  printf("latOut: %g %g %g ... %g %g\n", latOut[0]/DEG2RAD, latOut[1]/DEG2RAD, latOut[2]/DEG2RAD, latOut[nlatOut-2]/DEG2RAD, latOut[nlatOut-1]/DEG2RAD);

#if defined(HAVE_LIBYAC)

  //--------------------------------------------
  // define a grid
  //--------------------------------------------
  unsigned num_source_cells[2];
  unsigned num_target_cells[2];
  num_source_cells[0] = nlonIn;
  num_source_cells[1] = nlatIn;
  num_target_cells[0] = nlonOut;
  num_target_cells[1] = nlatOut;

  unsigned cyclic[2] = {1,0};
  struct grid *source_grid, *target_grid;

  source_grid = reg2d_grid_new(lonIn, latIn, num_source_cells, cyclic);
  target_grid = reg2d_grid_new(lonOut, latOut, num_target_cells, cyclic);

  if ( cdoTimer ) timer_stop(timer_yar_remap_init);

  struct points source_points, target_points;

  //--------------------------------------------
  // define points
  //--------------------------------------------
  // init_points(&source_points, &source_grid, CELL, lonIn, latIn);
  // init_points(&target_points, &target_grid, CELL, lonOut, latOut);

  //--------------------------------------------
  // initialise interpolation
  //--------------------------------------------

  struct dep_list deps;
  unsigned search_id;
  //struct interpolation interpolation;

  // printf("src num_grid_corners %d\n", get_num_grid_corners(*get_point_grid(&source_points)));
  //printf("tgt num_grid_corners %d\n", get_num_grid_corners(*get_point_grid(&target_points)));

  struct grid_search *search;

  if ( cdoTimer ) timer_start(timer_yar_remap_con);

  search = bucket_search_new(source_grid);
  // search_id = search_init(&source_grid);
 
  do_cell_search(search, target_grid, &deps);


  printf("total_num_dependencies: %d\n", get_total_num_dependencies(deps));
  int num_elements = deps.num_elements;
  printf("dep num elements: %d\n", deps.num_elements);

  enum edge_type quad_type[] = {GREAT_CIRCLE, GREAT_CIRCLE, GREAT_CIRCLE, GREAT_CIRCLE};

  int n;
  double weight_sum;
  double const epsilon = 1.0e-10; // relative precision 

  double *weight;
  weight = (double*) malloc(gridsize1*sizeof(double));

  double tgt_area;
  double *area;
  area = (double*) malloc(gridsize1*sizeof(double));

  struct grid_cell *SourceCell;
  SourceCell = (struct grid_cell*) malloc(gridsize1  * sizeof(struct grid_cell));

  for ( int n = 0; n <  gridsize1; n++ ) {
    SourceCell[n].num_corners   = 4;
    SourceCell[n].edge_type     = quad_type;
    SourceCell[n].coordinates_x = (double*) malloc( 4 * sizeof(double));
    SourceCell[n].coordinates_y = (double*) malloc( 4 * sizeof(double));
  }

  struct grid_cell  TargetCell;

  TargetCell.num_corners   = 4;
  TargetCell.edge_type     = quad_type;

  TargetCell.coordinates_x = (double*) malloc( 4 * sizeof(double));
  TargetCell.coordinates_y = (double*) malloc( 4 * sizeof(double));

  unsigned const * curr_deps;
  //struct polygons polygons;

  //polygon_create ( &polygons );

  for ( int i = 0; i < num_elements; ++i )
    {
      int index2 = i;
      int ilat2 = index2/nlonOut;
      int ilon2 = index2 - ilat2*nlonOut;

      TargetCell.coordinates_x[0] =  lonOut[ilon2];
      TargetCell.coordinates_y[0] =  latOut[ilat2];
      TargetCell.coordinates_x[1] =  lonOut[ilon2+1];
      TargetCell.coordinates_y[1] =  latOut[ilat2];
      TargetCell.coordinates_x[2] =  lonOut[ilon2+1];
      TargetCell.coordinates_y[2] =  latOut[ilat2+1];
      TargetCell.coordinates_x[3] =  lonOut[ilon2];
      TargetCell.coordinates_y[3] =  latOut[ilat2+1];

      if ( cdoVerbose )
	{
	  printf("target:\n");
	  for ( int n = 0; n < 4; ++n )
	    printf(" %g %g", TargetCell.coordinates_x[n]/DEG2RAD, TargetCell.coordinates_y[n]/DEG2RAD);
	  printf("\n");
	}

      if ( cdoVerbose )
	printf("num_deps_per_element %d %d\n", i, deps.num_deps_per_element[i]);
      int num_deps = deps.num_deps_per_element[i];
      int nSourceCells = num_deps;

      if ( num_deps > 0 ) curr_deps = get_dependencies_of_element(deps, i);
      for ( int k = 0; k < num_deps; ++k )
	{
	  int index1 = curr_deps[k];
	  int ilat1 = index1/nlonIn;
	  int ilon1 = index1 - ilat1*nlonIn;
	  if ( cdoVerbose )
	    printf("  dep: %d %d %d %d %d %d\n", k, nlonOut, nlatOut, index1, ilon1, ilat1);
	
	  SourceCell[k].coordinates_x[0] =  lonIn[ilon1];
	  SourceCell[k].coordinates_y[0] =  latIn[ilat1];
	  SourceCell[k].coordinates_x[1] =  lonIn[ilon1+1];
	  SourceCell[k].coordinates_y[1] =  latIn[ilat1];
	  SourceCell[k].coordinates_x[2] =  lonIn[ilon1+1];
	  SourceCell[k].coordinates_y[2] =  latIn[ilat1+1];
	  SourceCell[k].coordinates_x[3] =  lonIn[ilon1];
	  SourceCell[k].coordinates_y[3] =  latIn[ilat1+1];
	  if ( cdoVerbose )
	    {
	      printf("source: %d\n", k);
	      for ( int n = 0; n < 4; ++n )
		printf(" %g %g", SourceCell[k].coordinates_x[n]/DEG2RAD, SourceCell[k].coordinates_y[n]/DEG2RAD);
	      printf("\n");
	    }
	}
      /*
      polygon_partial_weights(nSourceCells, SourceCell, TargetCell, weight, &polygons);

      correct_weights ( nSourceCells, weight );
      */
      compute_overlap_areas ( nSourceCells, SourceCell, TargetCell, area);

      // tgt_area = huiliers_area(TargetCell);
      tgt_area = cell_area(TargetCell);
      for (n = 0; n < nSourceCells; ++n)
	weight[n] = area[n] / tgt_area;
      /*
      weight_sum = 0.0;
      for ( n = 0; n < nSourceCells; ++n ) weight_sum += weight[n];
      
      if ( fabs(weight_sum-1.0) > epsilon ) {
	printf ("test 2.) weight %lf\n", weight_sum );
	PUT_ERR("test 2.) weight deviates from expected value of 1.0!\n");
      }
      */
      correct_weights ( nSourceCells, weight );
      /*
      weight_sum = 0.0;
      for ( n = 0; n < nSourceCells; ++n ) weight_sum += weight[n];
      
      if ( fabs(weight_sum-1.0) > epsilon ) {
	printf ("test 2.) weight %lf\n", weight_sum );
	PUT_ERR("test 2.) weight deviates from expected value of 1.0!\n");
      }
      */
      for ( int k = 0; k < num_deps; ++k )
	{
	  int index1 = curr_deps[k];
	  int ilat1 = index1/nlonIn;
	  int ilon1 = index1 - ilat1*nlonIn;
	  long add1, add2;

	  add1 = index1;
	  add2 = index2;

	  yar_store_link_cnsrv(&remap.vars, add1, add2, weight[k]);

	  if ( cdoVerbose )
	    printf("  result dep: %d %d %d %d %d %d  %g\n", k, nlonOut, nlatOut, index1, ilon1, ilat1, weight[k]);
	}
      // correct_weights ( nSourceCells, weight );
    }

  //polygon_destroy ( &polygons );
  if ( cdoTimer ) timer_stop(timer_yar_remap_con);

  if ( cdoTimer ) timer_start(timer_yar_remap);
  yar_remap(array2, missval, gridInqSize(gridIDout), remap.vars.num_links, remap.vars.wts,
	    remap.vars.num_wts, remap.vars.tgt_cell_add, remap.vars.src_cell_add, array1);
  if ( cdoTimer ) timer_stop(timer_yar_remap);

  nmiss = 0;
  for ( int i = 0; i < gridInqSize(gridIDout); ++i )
    if ( DBL_IS_EQUAL(array2[i], missval) ) nmiss++;

  field2->nmiss = nmiss;

  //---------------
  // cleanup
  //---------------
  delete_grid_search(search);
  free_dep_list(&deps);
  delete_grid(source_grid);
  delete_grid(target_grid);

  free(weight);
  free(area);

  //free(lonIn);
  //free(latIn);
  //free(lonOut);
  //free(latOut);
  //free(fieldIn);
#endif
}


void *YAR(void *argument)
{
  int YARBIL, YARCON;
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

  if ( cdoTimer )
    {
      timer_yar_remap        = timer_new("yar remap");
      timer_yar_remap_init   = timer_new("yar remap init");
      timer_yar_remap_sort   = timer_new("yar remap sort");
      timer_yar_remap_con    = timer_new("yar remap con");
      timer_yar_remap_bil    = timer_new("yar remap bil");
    }

  cdoInitialize(argument);

  YARBIL = cdoOperatorAdd("yarbil",  0, 0, NULL);
  YARCON = cdoOperatorAdd("yarcon",  0, 0, NULL);

  operatorID = cdoOperatorID();

  operatorInputArg("grid description file or name");
  gridID2 = cdoDefineGrid(operatorArgv()[0]);

  field_init(&field1);
  field_init(&field2);

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

      if ( gridInqType(gridID1) != GRID_LONLAT && gridInqType(gridID1) != GRID_GAUSSIAN )
	cdoAbort("Interpolation of %s data unsupported!", gridNamePtr(gridInqType(gridID1)) );

      if ( gridIsRotated(gridID1) )
	cdoAbort("Rotated grids not supported!");

      vlistChangeGridIndex(vlistID2, index, gridID2);
    }

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  gridsize = vlistGridsizeMax(vlistID1);
  array1   = (double*) malloc(gridsize*sizeof(double));

  gridsize = gridInqSize(gridID2);
  array2   = (double*) malloc(gridsize*sizeof(double));

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

	  if ( operatorID == YARBIL )
	    yar_remap_bil(&field1, &field2);
	  else if ( operatorID == YARCON )
	    yar_remap_con(&field1, &field2);
	  else
	    cdoAbort("Not implemented!");

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
