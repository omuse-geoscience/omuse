/*
  This is a C library of the Fortran SCRIP version 1.4

  ===>>> Please send bug reports to <http://mpimet.mpg.de/cdo> <<<===

  Spherical Coordinate Remapping and Interpolation Package (SCRIP)
  ================================================================

  SCRIP is a software package which computes addresses and weights for
  remapping and interpolating fields between grids in spherical coordinates.
  It was written originally for remapping fields to other grids in a coupled
  climate model, but is sufficiently general that it can be used in other 
  applications as well. The package should work for any grid on the surface
  of a sphere. SCRIP currently supports four remapping options:

  Conservative remapping
  ----------------------
  First- and second-order conservative remapping as described in
  Jones (1999, Monthly Weather Review, 127, 2204-2210).

  Bilinear interpolation
  ----------------------
  Slightly generalized to use a local bilinear approximation
  (only logically-rectangular grids).

  Bicubic interpolation
  ----------------------
  Similarly generalized (only logically-rectangular grids).

  Distance-weighted averaging
  ---------------------------
  Distance-weighted average of a user-specified number of nearest neighbor values.

  Documentation
  =============

  http://climate.lanl.gov/Software/SCRIP/SCRIPusers.pdf

*/
/*
  2013-11-08 Uwe Schulzweida: split remapgrid class to src_grid and tgt_grid
  2012-01-16 Uwe Schulzweida: alloc grid2_bound_box only for conservative remapping
  2011-01-07 Uwe Schulzweida: Changed remap weights from 2D to 1D array
  2009-05-25 Uwe Schulzweida: Changed restrict data type from double to int
  2009-01-11 Uwe Schulzweida: OpenMP parallelization
 */

#if defined(HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <limits.h>
#include <time.h>

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "grid.h"
#include "remap.h"
#include "remap_store_link_cnsrv.h"
#include "util.h"  /* progressStatus */


#define IS_REG2D_GRID(gridID)  (!gridIsRotated(gridID) && (gridInqType(gridID) == GRID_LONLAT || gridInqType(gridID) == GRID_GAUSSIAN))



static int  remap_gen_weights     = TRUE;
static int  remap_write_remap     = FALSE;
static int  remap_num_srch_bins   = 180;
#define  DEFAULT_MAX_ITER  100
long remap_max_iter        = DEFAULT_MAX_ITER;  /* Max iteration count for i, j iteration */

void remap_set_int(int remapvar, int value)
{
  if      ( remapvar == REMAP_STORE_LINK_FAST ) remap_store_link_fast = value;
  else if ( remapvar == REMAP_WRITE_REMAP     ) remap_write_remap     = value;
  else if ( remapvar == REMAP_MAX_ITER        ) remap_max_iter        = value;
  else if ( remapvar == REMAP_NUM_SRCH_BINS   ) remap_num_srch_bins   = value;
  else if ( remapvar == REMAP_GENWEIGHTS      ) remap_gen_weights     = value;
  else      cdoAbort("Unsupported remap variable (%d)!", remapvar);
}

double intlin(double x, double y1, double x1, double y2, double x2);

void remapGridFree(remapgrid_t *grid)
{
  if ( grid->vgpm ) free(grid->vgpm);
  if ( grid->mask ) free(grid->mask);

  if ( grid->reg2d_center_lat ) free(grid->reg2d_center_lat);
  if ( grid->reg2d_center_lon ) free(grid->reg2d_center_lon);
  if ( grid->reg2d_corner_lat ) free(grid->reg2d_corner_lat);
  if ( grid->reg2d_corner_lon ) free(grid->reg2d_corner_lon);

  if ( grid->cell_center_lat ) free(grid->cell_center_lat);
  if ( grid->cell_center_lon ) free(grid->cell_center_lon);
  if ( grid->cell_corner_lat ) free(grid->cell_corner_lat);
  if ( grid->cell_corner_lon ) free(grid->cell_corner_lon);

  if ( grid->cell_area ) free(grid->cell_area);
  if ( grid->cell_frac ) free(grid->cell_frac);

  if ( grid->cell_bound_box ) free(grid->cell_bound_box);

  if ( grid->bin_addr ) free(grid->bin_addr);
  if ( grid->bin_lats ) free(grid->bin_lats);

} /* remapGridFree */

/*****************************************************************************/

void remapVarsFree(remapvars_t *rv)
{
  long i, num_blks;

  if ( rv->pinit == TRUE )
    {
      rv->pinit = FALSE;

      rv->sort_add = FALSE;

      if ( rv->src_cell_add ) free(rv->src_cell_add);
      if ( rv->tgt_cell_add ) free(rv->tgt_cell_add);
      if ( rv->wts ) free(rv->wts);

      if ( rv->links.option == TRUE )
	{
	  rv->links.option = FALSE;

	  if ( rv->links.num_blks )
	    {
	      free(rv->links.num_links);
	      num_blks = rv->links.num_blks;
	      for ( i = 0; i < num_blks; ++i )
		{
		  free(rv->links.src_add[i]);
		  free(rv->links.dst_add[i]);
		  free(rv->links.w_index[i]);
		}
	      free(rv->links.src_add);
	      free(rv->links.dst_add);
	      free(rv->links.w_index);
	    }
	}
    }
  else
    fprintf(stderr, "%s Warning: vars not initialized!\n", __func__);

} /* remapVarsFree */

/*****************************************************************************/

void remapgrid_init(remapgrid_t *grid)
{
  grid->remap_grid_type  = -1;
  grid->num_srch_bins    = remap_num_srch_bins; // only for source grid ?

  grid->num_cell_corners = 0;
  grid->luse_cell_corners  = FALSE;
  grid->lneed_cell_corners = FALSE;

  grid->nvgp             = 0;
  grid->vgpm             = NULL;

  grid->mask             = NULL;

  grid->reg2d_center_lon = NULL;
  grid->reg2d_center_lat = NULL;
  grid->reg2d_corner_lon = NULL;
  grid->reg2d_corner_lat = NULL;

  grid->cell_center_lon  = NULL;
  grid->cell_center_lat  = NULL;
  grid->cell_corner_lon  = NULL;
  grid->cell_corner_lat  = NULL;

  grid->cell_area        = NULL;
  grid->cell_frac        = NULL;

  grid->cell_bound_box   = NULL;

  grid->bin_addr         = NULL;
  grid->bin_lats         = NULL;
}

/*****************************************************************************/

void remapGridRealloc(int map_type, remapgrid_t *grid)
{
  long nalloc;

  if ( grid->nvgp )
    grid->vgpm   = (int*) realloc(grid->vgpm, grid->nvgp*sizeof(int));

  grid->mask     = (int*) realloc(grid->mask, grid->size*sizeof(int));

  if ( remap_write_remap == TRUE || grid->remap_grid_type != REMAP_GRID_TYPE_REG2D )
    {
      grid->cell_center_lon = (double*) realloc(grid->cell_center_lon, grid->size*sizeof(double));
      grid->cell_center_lat = (double*) realloc(grid->cell_center_lat, grid->size*sizeof(double));
    }

  if ( map_type == MAP_TYPE_CONSERV || map_type == MAP_TYPE_CONSERV_YAC )
    {
      grid->cell_area = (double*) realloc(grid->cell_area, grid->size*sizeof(double));
      memset(grid->cell_area, 0, grid->size*sizeof(double));
    }

  grid->cell_frac = (double*) realloc(grid->cell_frac, grid->size*sizeof(double));
  memset(grid->cell_frac, 0, grid->size*sizeof(double));

  if ( grid->lneed_cell_corners )
    {
      if ( grid->num_cell_corners == 0 )
	{
	  cdoAbort("Grid cell corner missing!");
	}
      else
	{
	  nalloc = grid->num_cell_corners*grid->size;

	  grid->cell_corner_lon = (double*) realloc(grid->cell_corner_lon, nalloc*sizeof(double));
	  memset(grid->cell_corner_lon, 0, nalloc*sizeof(double));

	  grid->cell_corner_lat = (double*) realloc(grid->cell_corner_lat, nalloc*sizeof(double));  
	  memset(grid->cell_corner_lat, 0, nalloc*sizeof(double));
	}
    }
}

/*****************************************************************************/
static
void boundbox_from_corners(long size, long nc, const double *restrict corner_lon,
			   const double *restrict corner_lat, restr_t *restrict bound_box)
{
  long i4, inc, i, j;
  restr_t clon, clat;

#if defined(_OPENMP)
#pragma omp parallel for default(none)        \
  shared(bound_box, corner_lat, corner_lon, nc, size)	\
  private(i4, inc, i, j, clon, clat)
#endif
  for ( i = 0; i < size; ++i )
    {
      i4 = i<<2; // *4
      inc = i*nc;
      clat = RESTR_SCALE(corner_lat[inc]);
      clon = RESTR_SCALE(corner_lon[inc]);
      bound_box[i4  ] = clat;
      bound_box[i4+1] = clat;
      bound_box[i4+2] = clon;
      bound_box[i4+3] = clon;
      for ( j = 1; j < nc; ++j )
	{
	  clat = RESTR_SCALE(corner_lat[inc+j]);
	  clon = RESTR_SCALE(corner_lon[inc+j]);
	  if ( clat < bound_box[i4  ] ) bound_box[i4  ] = clat;
	  if ( clat > bound_box[i4+1] ) bound_box[i4+1] = clat;
	  if ( clon < bound_box[i4+2] ) bound_box[i4+2] = clon;
	  if ( clon > bound_box[i4+3] ) bound_box[i4+3] = clon;
	}
    }
}

static
void boundbox_from_center(int lonIsCyclic, long size, long nx, long ny, const double *restrict center_lon,
			  const double *restrict center_lat, restr_t *restrict bound_box)
{
  long n4, i, j, k, n, ip1, jp1;
  long n_add, e_add, ne_add;
  restr_t tmp_lats[4], tmp_lons[4];  /* temps for computing bounding boxes */

#if defined(_OPENMP)
#pragma omp parallel for default(none)        \
  shared(lonIsCyclic, size, nx, ny, center_lon, center_lat, bound_box)	\
  private(n4, i, j, k, n, ip1, jp1, n_add, e_add, ne_add, tmp_lats, tmp_lons)
#endif
  for ( n = 0; n < size; n++ )
    {
      n4 = n<<2;

      /* Find N,S and NE points to this grid point */
      
      j = n/nx;
      i = n - j*nx;

      if ( i < (nx-1) )
	ip1 = i + 1;
      else
	{
	  /* 2009-01-09 Uwe Schulzweida: bug fix */
	  if ( lonIsCyclic )
	    ip1 = 0;
	  else
	    ip1 = i;
	}

      if ( j < (ny-1) )
	jp1 = j + 1;
      else
	{
	  /* 2008-12-17 Uwe Schulzweida: latitute cyclic ??? (bug fix) */
	  jp1 = j;
	}

      n_add  = jp1*nx + i;
      e_add  = j  *nx + ip1;
      ne_add = jp1*nx + ip1;

      /* Find N,S and NE lat/lon coords and check bounding box */

      tmp_lats[0] = RESTR_SCALE(center_lat[n]);
      tmp_lats[1] = RESTR_SCALE(center_lat[e_add]);
      tmp_lats[2] = RESTR_SCALE(center_lat[ne_add]);
      tmp_lats[3] = RESTR_SCALE(center_lat[n_add]);

      tmp_lons[0] = RESTR_SCALE(center_lon[n]);
      tmp_lons[1] = RESTR_SCALE(center_lon[e_add]);
      tmp_lons[2] = RESTR_SCALE(center_lon[ne_add]);
      tmp_lons[3] = RESTR_SCALE(center_lon[n_add]);

      bound_box[n4  ] = tmp_lats[0];
      bound_box[n4+1] = tmp_lats[0];
      bound_box[n4+2] = tmp_lons[0];
      bound_box[n4+3] = tmp_lons[0];

      for ( k = 1; k < 4; k++ )
	{
	  if ( tmp_lats[k] < bound_box[n4  ] ) bound_box[n4  ] = tmp_lats[k];
	  if ( tmp_lats[k] > bound_box[n4+1] ) bound_box[n4+1] = tmp_lats[k];
	  if ( tmp_lons[k] < bound_box[n4+2] ) bound_box[n4+2] = tmp_lons[k];
	  if ( tmp_lons[k] > bound_box[n4+3] ) bound_box[n4+3] = tmp_lons[k];
	}
    }
}

static
void check_lon_range2(long nc, long nlons, double *corners, double *centers)
{
  long n, k;

  assert(corners != NULL);
  /*
#if defined(_OPENMP)
#pragma omp parallel for default(none) shared(nlons, lons)
#endif
  */
  for ( n = 0; n < nlons; ++n )
    {
      bool lneg = false, lpos = false;

      for ( k = 0; k < nc; ++k )
	{
 	  if      ( !lneg && corners[n*nc+k] > -PI && corners[n*nc+k] < 0. ) lneg = true;
	  else if ( !lpos && corners[n*nc+k] <  PI && corners[n*nc+k] > 0. ) lpos = true;
	}

      if ( lneg && lpos )
	{
	  if ( centers[n] > PI )
	    for ( k = 0; k < nc; ++k ) corners[n*nc+k] += PI2;
	}
      else
	{
	  for ( k = 0; k < nc; ++k )
	    {
	      if ( corners[n*nc+k] > PI2  ) corners[n*nc+k] -= PI2;
	      if ( corners[n*nc+k] < ZERO ) corners[n*nc+k] += PI2;
	    }
	}
    }
}

static
void check_lon_range(long nlons, double *lons)
{
  long n;

  assert(lons != NULL);

#if defined(_OPENMP)
#pragma omp parallel for default(none) shared(nlons, lons)
#endif
  for ( n = 0; n < nlons; ++n )
    {
      // remove missing values
      if ( lons[n] < -PI2 ) lons[n] = 0;
      if ( lons[n] > 2*PI2) lons[n] = PI2;

      if ( lons[n] > PI2  ) lons[n] -= PI2;
      if ( lons[n] < ZERO ) lons[n] += PI2;
    }
}

static
void check_lat_range(long nlats, double *lats)
{
  long n;

  assert(lats != NULL);

#if defined(_OPENMP)
#pragma omp parallel for default(none) shared(nlats, lats)
#endif
  for ( n = 0; n < nlats; ++n )
    {
      if ( lats[n] >  PIH ) lats[n] =  PIH;
      if ( lats[n] < -PIH ) lats[n] = -PIH;
    }
}

static
void check_lon_boundbox_range(long nlons, restr_t *bound_box)
{
  long n, n4;

  assert(bound_box != NULL);

#if defined(_OPENMP)
#pragma omp parallel for default(none) shared(nlons, bound_box) private(n4)
#endif
  for ( n = 0; n < nlons; ++n )
    {
      n4 = n<<2;
      if ( RESTR_ABS(bound_box[n4+3] - bound_box[n4+2]) > RESTR_SCALE(PI) )
	{
	  bound_box[n4+2] = RESTR_SCALE(0.);
	  bound_box[n4+3] = RESTR_SCALE(PI2);
	}
    }
}

static
void check_lat_boundbox_range(long nlats, restr_t *restrict bound_box, double *restrict lats)
{
  long n, n4;

  assert(bound_box != NULL);

#if defined(_OPENMP)
#pragma omp parallel for default(none) shared(nlats, bound_box, lats) private(n4)
#endif
  for ( n = 0; n < nlats; ++n )
    {
      n4 = n<<2;
      if ( RESTR_SCALE(lats[n]) < bound_box[n4  ] ) bound_box[n4  ] = RESTR_SCALE(-PIH);
      if ( RESTR_SCALE(lats[n]) > bound_box[n4+1] ) bound_box[n4+1] = RESTR_SCALE( PIH);
    }
}

static
int expand_lonlat_grid(int gridID)
{
  char units[CDI_MAX_NAME];
  int gridIDnew;
  long nx, ny, nxp4, nyp4;
  double *xvals, *yvals;

  nx = gridInqXsize(gridID);
  ny = gridInqYsize(gridID);
  nxp4 = nx+4;
  nyp4 = ny+4;

  xvals = (double*) malloc(nxp4*sizeof(double));
  yvals = (double*) malloc(nyp4*sizeof(double));
  gridInqXvals(gridID, xvals+2);
  gridInqYvals(gridID, yvals+2);

  gridIDnew = gridCreate(GRID_LONLAT, nxp4*nyp4);
  gridDefXsize(gridIDnew, nxp4);
  gridDefYsize(gridIDnew, nyp4);
	      
  gridInqXunits(gridID,    units);
  gridDefXunits(gridIDnew, units);
  gridInqYunits(gridID,    units);
  gridDefYunits(gridIDnew, units);

  xvals[0] = xvals[2] - 2*gridInqXinc(gridID);
  xvals[1] = xvals[2] - gridInqXinc(gridID);
  xvals[nxp4-2] = xvals[nx+1] + gridInqXinc(gridID);
  xvals[nxp4-1] = xvals[nx+1] + 2*gridInqXinc(gridID);

  yvals[0] = yvals[2] - 2*gridInqYinc(gridID);
  yvals[1] = yvals[2] - gridInqYinc(gridID);
  yvals[nyp4-2] = yvals[ny+1] + gridInqYinc(gridID);
  yvals[nyp4-1] = yvals[ny+1] + 2*gridInqYinc(gridID);

  gridDefXvals(gridIDnew, xvals);
  gridDefYvals(gridIDnew, yvals);

  free(xvals);
  free(yvals);

  if ( gridIsRotated(gridID) )
    {
      gridDefXpole(gridIDnew, gridInqXpole(gridID));
      gridDefYpole(gridIDnew, gridInqYpole(gridID));
      gridDefAngle(gridIDnew, gridInqAngle(gridID));
    }

  return(gridIDnew);
}

static
int expand_curvilinear_grid(int gridID)
{
  char units[CDI_MAX_NAME];
  int gridIDnew;
  long gridsize, gridsize_new;
  long nx, ny, nxp4, nyp4;
  long i, j;
  double *xvals, *yvals;

  gridsize = gridInqSize(gridID);
  nx = gridInqXsize(gridID);
  ny = gridInqYsize(gridID);
  nxp4 = nx+4;
  nyp4 = ny+4;
  gridsize_new = gridsize + 4*(nx+2) + 4*(ny+2);

  xvals = (double*) malloc(gridsize_new*sizeof(double));
  yvals = (double*) malloc(gridsize_new*sizeof(double));
  gridInqXvals(gridID, xvals);
  gridInqYvals(gridID, yvals);

  gridIDnew = gridCreate(GRID_CURVILINEAR, nxp4*nyp4);
  gridDefXsize(gridIDnew, nxp4);
  gridDefYsize(gridIDnew, nyp4);

  gridInqXunits(gridID,   units);
  gridDefXunits(gridIDnew, units);
  gridInqYunits(gridID,   units);
  gridDefYunits(gridIDnew, units);

  for ( j = ny-1; j >= 0; j-- )
    for ( i = nx-1; i >= 0; i-- )
      xvals[(j+2)*(nx+4)+i+2] = xvals[j*nx+i];

  for ( j = ny-1; j >= 0; j-- )
    for ( i = nx-1; i >= 0; i-- )
      yvals[(j+2)*(nx+4)+i+2] = yvals[j*nx+i];

  for ( j = 2; j < nyp4-2; j++ )
    {
      xvals[j*nxp4  ] = intlin(3.0, xvals[j*nxp4+3], 0.0, xvals[j*nxp4+2], 1.0);
      xvals[j*nxp4+1] = intlin(2.0, xvals[j*nxp4+3], 0.0, xvals[j*nxp4+2], 1.0); 
      yvals[j*nxp4  ] = intlin(3.0, yvals[j*nxp4+3], 0.0, yvals[j*nxp4+2], 1.0); 
      yvals[j*nxp4+1] = intlin(2.0, yvals[j*nxp4+3], 0.0, yvals[j*nxp4+2], 1.0); 

      xvals[j*nxp4+nxp4-2] = intlin(2.0, xvals[j*nxp4+nxp4-4], 0.0, xvals[j*nxp4+nxp4-3], 1.0); 
      xvals[j*nxp4+nxp4-1] = intlin(3.0, xvals[j*nxp4+nxp4-4], 0.0, xvals[j*nxp4+nxp4-3], 1.0); 
      yvals[j*nxp4+nxp4-2] = intlin(2.0, yvals[j*nxp4+nxp4-4], 0.0, yvals[j*nxp4+nxp4-3], 1.0); 
      yvals[j*nxp4+nxp4-1] = intlin(3.0, yvals[j*nxp4+nxp4-4], 0.0, yvals[j*nxp4+nxp4-3], 1.0); 
    }

  for ( i = 0; i < nxp4; i++ )
    {
      xvals[0*nxp4+i] = intlin(3.0, xvals[3*nxp4+i], 0.0, xvals[2*nxp4+i], 1.0);
      xvals[1*nxp4+i] = intlin(2.0, xvals[3*nxp4+i], 0.0, xvals[2*nxp4+i], 1.0);
      yvals[0*nxp4+i] = intlin(3.0, yvals[3*nxp4+i], 0.0, yvals[2*nxp4+i], 1.0);
      yvals[1*nxp4+i] = intlin(2.0, yvals[3*nxp4+i], 0.0, yvals[2*nxp4+i], 1.0);

      xvals[(nyp4-2)*nxp4+i] = intlin(2.0, xvals[(nyp4-4)*nxp4+i], 0.0, xvals[(nyp4-3)*nxp4+i], 1.0);
      xvals[(nyp4-1)*nxp4+i] = intlin(3.0, xvals[(nyp4-4)*nxp4+i], 0.0, xvals[(nyp4-3)*nxp4+i], 1.0);
      yvals[(nyp4-2)*nxp4+i] = intlin(2.0, yvals[(nyp4-4)*nxp4+i], 0.0, yvals[(nyp4-3)*nxp4+i], 1.0);
      yvals[(nyp4-1)*nxp4+i] = intlin(3.0, yvals[(nyp4-4)*nxp4+i], 0.0, yvals[(nyp4-3)*nxp4+i], 1.0);
    }
  /*
    {
    FILE *fp;
    fp = fopen("xvals.asc", "w");
    for ( i = 0; i < gridsize_new; i++ ) fprintf(fp, "%g\n", xvals[i]);
    fclose(fp);
    fp = fopen("yvals.asc", "w");
    for ( i = 0; i < gridsize_new; i++ ) fprintf(fp, "%g\n", yvals[i]);
    fclose(fp);
    }
  */
  gridDefXvals(gridIDnew, xvals);
  gridDefYvals(gridIDnew, yvals);
  
  free(xvals);
  free(yvals);

  return(gridIDnew);
}

/*****************************************************************************/

static
void grid_check_lat_borders_rad(int n, double *ybounds)
{
  if ( ybounds[0] > ybounds[n-1] )
    {
      if ( RAD2DEG*ybounds[0]   >  PIH ) ybounds[0]   =  PIH;
      if ( RAD2DEG*ybounds[n-1] < -PIH ) ybounds[n-1] = -PIH;
    }
  else
    {
      if ( RAD2DEG*ybounds[0]   < -PIH ) ybounds[0]   = -PIH;
      if ( RAD2DEG*ybounds[n-1] >  PIH ) ybounds[n-1] =  PIH;
    }
}

static
void remap_define_reg2d(int gridID, remapgrid_t *grid)
{
  char unitstr[CDI_MAX_NAME];
  long nx, nxm, ny;
  long nxp1, nyp1;

  nx = grid->dims[0];
  ny = grid->dims[1];

  nxp1 = nx + 1;
  nyp1 = ny + 1;

  nxm = nx;
  if ( grid->is_cyclic ) nxm++;

  if ( grid->size != nx*ny ) cdoAbort("Internal error, wrong dimensions!");

  grid->reg2d_center_lon = (double*) realloc(grid->reg2d_center_lon, nxm*sizeof(double));
  grid->reg2d_center_lat = (double*) realloc(grid->reg2d_center_lat,  ny*sizeof(double));
 
  grid->reg2d_center_lon[0] = 0;
  grid->reg2d_center_lat[0] = 0;
  gridInqXvals(gridID, grid->reg2d_center_lon);
  gridInqYvals(gridID, grid->reg2d_center_lat);

  /* Convert lat/lon units if required */

  gridInqYunits(gridID, unitstr);

  grid_to_radian(unitstr, nx, grid->reg2d_center_lon, "grid reg2d center lon"); 
  grid_to_radian(unitstr, ny, grid->reg2d_center_lat, "grid reg2d center lat"); 

  if ( grid->is_cyclic ) grid->reg2d_center_lon[nx] = grid->reg2d_center_lon[0] + PI2;

  grid->reg2d_corner_lon = (double*) malloc(nxp1*sizeof(double));
  grid->reg2d_corner_lat = (double*) malloc(nyp1*sizeof(double));

  grid_gen_corners(nx, grid->reg2d_center_lon, grid->reg2d_corner_lon);
  grid_gen_corners(ny, grid->reg2d_center_lat, grid->reg2d_corner_lat);
  grid_check_lat_borders_rad(ny+1, grid->reg2d_corner_lat);

  //for ( long i = 0; i < nxp1; ++i ) printf("lon %ld %g\n", i, grid->reg2d_corner_lon[i]);
  //for ( long i = 0; i < nyp1; ++i ) printf("lat %ld %g\n", i, grid->reg2d_corner_lat[i]);

}

static
void remap_define_grid(int map_type, int gridID, remapgrid_t *grid)
{
  char xunitstr[CDI_MAX_NAME];
  char yunitstr[CDI_MAX_NAME];
  long gridsize;
  long i;
  int lgrid_destroy = FALSE;
  int lgrid_gen_bounds = FALSE;
  int gridID_gme = -1;

  if ( gridInqType(grid->gridID) != GRID_UNSTRUCTURED && gridInqType(grid->gridID) != GRID_CURVILINEAR )
    {
      if ( gridInqType(grid->gridID) == GRID_GME )
	{
	  gridID_gme = gridToUnstructured(grid->gridID, 1);
	  grid->nvgp = gridInqSize(gridID_gme);
	  gridID = gridDuplicate(gridID_gme);
	  gridCompress(gridID);
	  grid->luse_cell_corners = TRUE;
	}
      else if ( remap_write_remap == TRUE || grid->remap_grid_type != REMAP_GRID_TYPE_REG2D )
	{
	  lgrid_destroy = TRUE;
	  gridID = gridToCurvilinear(grid->gridID, 1);
	  lgrid_gen_bounds = TRUE;
	}
    }

  gridsize = grid->size = gridInqSize(gridID);

  grid->dims[0] = gridInqXsize(gridID);
  grid->dims[1] = gridInqYsize(gridID);
  if ( gridInqType(grid->gridID) != GRID_UNSTRUCTURED )
    {
      if ( grid->dims[0] == 0 ) cdoAbort("%s grid without longitude coordinates!", gridNamePtr(gridInqType(grid->gridID)));
      if ( grid->dims[1] == 0 ) cdoAbort("%s grid without latitude coordinates!", gridNamePtr(gridInqType(grid->gridID)));
    }

  grid->is_cyclic = gridIsCircular(gridID);

  if ( gridInqType(gridID) == GRID_UNSTRUCTURED )
    grid->rank = 1;
  else
    grid->rank = 2;

  if ( gridInqType(gridID) == GRID_UNSTRUCTURED )
    grid->num_cell_corners = gridInqNvertex(gridID);
  else
    grid->num_cell_corners = 4;

  remapGridRealloc(map_type, grid);

  /* Initialize logical mask */

#if defined(_OPENMP)
#pragma omp parallel for default(none) shared(gridsize, grid)
#endif
  for ( i = 0; i < gridsize; ++i ) grid->mask[i] = TRUE;

  if ( gridInqMask(gridID, NULL) )
    {
      int *mask = (int*) malloc(gridsize*sizeof(int));
      gridInqMask(gridID, mask);
      for ( i = 0; i < gridsize; ++i )
	if ( mask[i] == 0 ) grid->mask[i] = FALSE;
      free(mask);
    }

  if ( remap_write_remap == FALSE && grid->remap_grid_type == REMAP_GRID_TYPE_REG2D ) return;


  gridInqXvals(gridID, grid->cell_center_lon);
  gridInqYvals(gridID, grid->cell_center_lat);

  gridInqXunits(gridID, xunitstr);
  gridInqYunits(gridID, yunitstr);


  if ( grid->lneed_cell_corners )
    {
      if ( gridInqYbounds(gridID, NULL) && gridInqXbounds(gridID, NULL) )
	{
	  gridInqXbounds(gridID, grid->cell_corner_lon);
	  gridInqYbounds(gridID, grid->cell_corner_lat);
	}
      else if ( lgrid_gen_bounds )
	{
	  grid_cell_center_to_bounds_X2D(xunitstr, grid->dims[0], grid->dims[1], grid->cell_center_lon, grid->cell_corner_lon, 0);
	  grid_cell_center_to_bounds_Y2D(yunitstr, grid->dims[0], grid->dims[1], grid->cell_center_lat, grid->cell_corner_lat);
	}
      else
	{
	  cdoAbort("Grid corner missing!");
	}
    }


  if ( gridInqType(grid->gridID) == GRID_GME ) gridInqMaskGME(gridID_gme, grid->vgpm);

  /* Convert lat/lon units if required */

  grid_to_radian(xunitstr, grid->size, grid->cell_center_lon, "grid center lon"); 
  grid_to_radian(yunitstr, grid->size, grid->cell_center_lat, "grid center lat"); 
  /* Note: using units from cell center instead from bounds */
  if ( grid->num_cell_corners && grid->lneed_cell_corners )
    {
      grid_to_radian(xunitstr, grid->num_cell_corners*grid->size, grid->cell_corner_lon, "grid corner lon"); 
      grid_to_radian(yunitstr, grid->num_cell_corners*grid->size, grid->cell_corner_lat, "grid corner lat"); 
    }

  if ( lgrid_destroy ) gridDestroy(gridID);

  /* Convert longitudes to 0,2pi interval */

  check_lon_range(grid->size, grid->cell_center_lon);

  if ( grid->num_cell_corners && grid->lneed_cell_corners )
    {
      //check_lon_range2(grid->num_cell_corners, grid->size, grid->cell_corner_lon, grid->cell_center_lon);
	check_lon_range(grid->num_cell_corners*grid->size, grid->cell_corner_lon);
    }
  /*  Make sure input latitude range is within the machine values for +/- pi/2 */

  check_lat_range(grid->size, grid->cell_center_lat);

  if ( grid->num_cell_corners && grid->lneed_cell_corners )
    check_lat_range(grid->num_cell_corners*grid->size, grid->cell_corner_lat);
}

/*  Compute bounding boxes for restricting future grid searches */
static
void cell_bounding_boxes(remapgrid_t *grid, int remap_grid_basis)
{
  if ( remap_grid_basis == REMAP_GRID_BASIS_SRC || grid->luse_cell_corners )
    grid->cell_bound_box = (restr_t*) realloc(grid->cell_bound_box, 4*grid->size*sizeof(restr_t));

  if ( grid->luse_cell_corners )
    {
      if ( grid->lneed_cell_corners )
	{
	  if ( cdoVerbose ) cdoPrint("Grid: boundbox_from_corners");

	  boundbox_from_corners(grid->size, grid->num_cell_corners, 
				grid->cell_corner_lon, grid->cell_corner_lat, grid->cell_bound_box);
	}
      else /* full grid search */
	{
	  long gridsize;
	  long i, i4;
	  
	  gridsize = grid->size;
  
	  if ( cdoVerbose ) cdoPrint("Grid: bounds missing -> full grid search!");

	  for ( i = 0; i < gridsize; ++i )
	    {
	      i4 = i<<2;
	      grid->cell_bound_box[i4  ] = RESTR_SCALE(-PIH);
	      grid->cell_bound_box[i4+1] = RESTR_SCALE( PIH);
	      grid->cell_bound_box[i4+2] = RESTR_SCALE(0.);
	      grid->cell_bound_box[i4+3] = RESTR_SCALE(PI2);
	    }
	}
    }
  else if ( remap_grid_basis == REMAP_GRID_BASIS_SRC )
    {
      long nx, ny;

      if ( grid->rank != 2 ) cdoAbort("Internal problem, grid rank = %d!", grid->rank);

      nx = grid->dims[0];
      ny = grid->dims[1];

      if ( cdoVerbose ) cdoPrint("Grid: boundbox_from_center");

      boundbox_from_center(grid->is_cyclic, grid->size, nx, ny, 
			   grid->cell_center_lon, grid->cell_center_lat, grid->cell_bound_box);
    }

  if ( remap_grid_basis == REMAP_GRID_BASIS_SRC || grid->lneed_cell_corners )
    check_lon_boundbox_range(grid->size, grid->cell_bound_box);

  /* Try to check for cells that overlap poles */

  if ( remap_grid_basis == REMAP_GRID_BASIS_SRC || grid->lneed_cell_corners )
    check_lat_boundbox_range(grid->size, grid->cell_bound_box, grid->cell_center_lat);
}


void remap_grids_init(int map_type, int lextrapolate, int gridID1, remapgrid_t *src_grid, int gridID2, remapgrid_t *tgt_grid)
{
  int reg2d_src_gridID = gridID1;
  int reg2d_tgt_gridID = gridID2;

  /* Initialize remapgrid structure */
  remapgrid_init(src_grid);
  remapgrid_init(tgt_grid);

  if ( map_type == MAP_TYPE_BILINEAR || map_type == MAP_TYPE_BICUBIC ||
       map_type == MAP_TYPE_DISTWGT  || map_type == MAP_TYPE_CONSERV_YAC )
    {
      if ( IS_REG2D_GRID(gridID1) ) src_grid->remap_grid_type = REMAP_GRID_TYPE_REG2D;
      // src_grid->remap_grid_type = 0;
    }

  if ( src_grid->remap_grid_type == REMAP_GRID_TYPE_REG2D && map_type == MAP_TYPE_CONSERV_YAC )
    {
      if ( IS_REG2D_GRID(gridID2) ) tgt_grid->remap_grid_type = REMAP_GRID_TYPE_REG2D;
      // else src_grid->remap_grid_type = -1;
    }

  if ( remap_gen_weights == FALSE && map_type == MAP_TYPE_BILINEAR )
    {
      if ( IS_REG2D_GRID(gridID2) ) tgt_grid->remap_grid_type = REMAP_GRID_TYPE_REG2D;
      // else src_grid->remap_grid_type = -1;
    }

  if ( lextrapolate > 0 )
    src_grid->lextrapolate = TRUE;
  else
    src_grid->lextrapolate = FALSE;

  if ( map_type == MAP_TYPE_CONSERV || map_type == MAP_TYPE_CONSERV_YAC )
    {
      if ( src_grid->remap_grid_type != REMAP_GRID_TYPE_REG2D )
	{
	  src_grid->luse_cell_corners  = TRUE;
	  src_grid->lneed_cell_corners = TRUE;
	}

      if ( tgt_grid->remap_grid_type != REMAP_GRID_TYPE_REG2D )
	{
	  tgt_grid->luse_cell_corners  = TRUE;
	  tgt_grid->lneed_cell_corners = TRUE;
	}
    }

  src_grid->gridID = gridID1;
  tgt_grid->gridID = gridID2;

  if ( !src_grid->lextrapolate && gridInqSize(src_grid->gridID) > 1 &&
       map_type == MAP_TYPE_DISTWGT &&
       ((gridInqType(gridID1) == GRID_LONLAT && gridIsRotated(gridID1)) ||
	(gridInqType(gridID1) == GRID_LONLAT && src_grid->non_global)) )
    {
      src_grid->gridID = gridID1 = expand_lonlat_grid(gridID1);
      reg2d_src_gridID = gridID1;
    }

  if ( gridInqType(gridID1) == GRID_UNSTRUCTURED )
    {
      if ( gridInqYvals(gridID1, NULL) == 0 || gridInqXvals(gridID1, NULL) == 0 )
	{
	  if ( gridInqNumber(gridID1) > 0 )
	    {
	      src_grid->gridID = gridID1 = referenceToGrid(gridID1);
	      if ( gridID1 == -1 ) cdoAbort("Reference to source grid not found!");
	    }
	}
    }

  if ( gridInqType(gridID2) == GRID_UNSTRUCTURED )
    {
      if ( gridInqYvals(gridID2, NULL) == 0 || gridInqXvals(gridID2, NULL) == 0 )
	{
	  if ( gridInqNumber(gridID2) > 0 )
	    {
	      tgt_grid->gridID = gridID2 = referenceToGrid(gridID2);
	      if ( gridID2 == -1 ) cdoAbort("Reference to target grid not found!");
	    }
	}
    }

  if ( gridInqSize(src_grid->gridID) > 1 && 
       (gridInqType(src_grid->gridID) == GRID_LCC || 
	gridInqType(src_grid->gridID) == GRID_LAEA || 
	gridInqType(src_grid->gridID) == GRID_SINUSOIDAL) )
    {
      src_grid->gridID = gridID1 = gridToCurvilinear(src_grid->gridID, 1);
    }

  if ( !src_grid->lextrapolate && gridInqSize(src_grid->gridID) > 1 &&
       map_type == MAP_TYPE_DISTWGT &&
       (gridInqType(gridID1) == GRID_CURVILINEAR && src_grid->non_global) )
    {
      src_grid->gridID = gridID1 = expand_curvilinear_grid(gridID1);
    }

  if ( map_type == MAP_TYPE_DISTWGT )
    {
      if ( gridInqType(src_grid->gridID) == GRID_UNSTRUCTURED )
	{
	  src_grid->luse_cell_corners  = TRUE;
	  src_grid->lneed_cell_corners = FALSE; /* full grid search */
	}
      /* not used !
      if ( gridInqType(tgt_grid->gridID) == GRID_UNSTRUCTURED )
	{
	  tgt_grid->luse_cell_corners  = TRUE;
	  tgt_grid->lneed_cell_corners = FALSE;
	}
      */
    }

  //if ( src_grid->remap_grid_type != REMAP_GRID_TYPE_REG2D )
    remap_define_grid(map_type, gridID1, src_grid);

  remap_define_grid(map_type, gridID2, tgt_grid);

  if ( src_grid->remap_grid_type == REMAP_GRID_TYPE_REG2D && tgt_grid->remap_grid_type == REMAP_GRID_TYPE_REG2D )
    {
      remap_define_reg2d(reg2d_src_gridID, src_grid);
      remap_define_reg2d(reg2d_tgt_gridID, tgt_grid);
    }
  else if ( src_grid->remap_grid_type == REMAP_GRID_TYPE_REG2D )
    {
      remap_define_reg2d(reg2d_src_gridID, src_grid);
    }
  else
    {
      cell_bounding_boxes(src_grid, REMAP_GRID_BASIS_SRC);
      cell_bounding_boxes(tgt_grid, REMAP_GRID_BASIS_TGT);
      /*
	Set up and assign address ranges to search bins in order to further restrict later searches
      */
      calc_lat_bins(src_grid, tgt_grid, map_type);
    }

}  /* remapGridInit */

/*****************************************************************************/

/*
    This routine initializes some variables and provides an initial
    allocation of arrays (fairly large so frequent resizing unnecessary).
*/
void remap_vars_init(int map_type, long src_grid_size, long tgt_grid_size, remapvars_t *rv)
{
  /* Initialize all pointer */
  if ( rv->pinit == FALSE )
    {
      rv->pinit = TRUE;

      rv->src_cell_add = NULL;
      rv->tgt_cell_add = NULL;
      rv->wts          = NULL;
    }

  /* Determine the number of weights */

#if defined(_OPENMP)
  if ( ompNumThreads > 1 )
    {
      if      ( map_type == MAP_TYPE_CONSERV     ) rv->sort_add = TRUE;
      else if ( map_type == MAP_TYPE_CONSERV_YAC ) rv->sort_add = FALSE;
      else if ( map_type == MAP_TYPE_BILINEAR    ) rv->sort_add = FALSE;
      else if ( map_type == MAP_TYPE_BICUBIC     ) rv->sort_add = FALSE;
      else if ( map_type == MAP_TYPE_DISTWGT     ) rv->sort_add = FALSE;
      else cdoAbort("Unknown mapping method!");
    }
  else
#endif
    {
      if      ( map_type == MAP_TYPE_CONSERV     ) rv->sort_add = TRUE;
      else if ( map_type == MAP_TYPE_CONSERV_YAC ) rv->sort_add = FALSE;
      else if ( map_type == MAP_TYPE_BILINEAR    ) rv->sort_add = FALSE;
      else if ( map_type == MAP_TYPE_BICUBIC     ) rv->sort_add = FALSE;
      else if ( map_type == MAP_TYPE_DISTWGT     ) rv->sort_add = FALSE;
      else cdoAbort("Unknown mapping method!");
    }

  if      ( map_type == MAP_TYPE_CONSERV     ) rv->num_wts = 3;
  else if ( map_type == MAP_TYPE_CONSERV_YAC ) rv->num_wts = 1;
  else if ( map_type == MAP_TYPE_BILINEAR    ) rv->num_wts = 1;
  else if ( map_type == MAP_TYPE_BICUBIC     ) rv->num_wts = 4;
  else if ( map_type == MAP_TYPE_DISTWGT     ) rv->num_wts = 1;
  else cdoAbort("Unknown mapping method!");

   /*
    Initialize num_links and set max_links to four times the largest 
    of the destination grid sizes initially (can be changed later).
    Set a default resize increment to increase the size of link
    arrays if the number of links exceeds the initial size
  */
  rv->num_links = 0;
  rv->max_links = 4 * tgt_grid_size;

  rv->resize_increment = (int) (0.1 * MAX(src_grid_size, tgt_grid_size));

  /*  Allocate address and weight arrays for mapping 1 */
  if ( map_type == MAP_TYPE_CONSERV )
    {
      rv->src_cell_add = (int*) realloc(rv->src_cell_add, rv->max_links*sizeof(int));
      rv->tgt_cell_add = (int*) realloc(rv->tgt_cell_add, rv->max_links*sizeof(int));

      rv->wts = (double*) realloc(rv->wts, rv->num_wts*rv->max_links*sizeof(double));
    }

  rv->links.option    = FALSE;
  rv->links.max_links = 0;
  rv->links.num_blks  = 0;
  rv->links.num_links = NULL;
  rv->links.src_add   = NULL;
  rv->links.dst_add   = NULL;
  rv->links.w_index   = NULL;

} /* remapVarsInit */

/*****************************************************************************/

/*
   This routine resizes remapping arrays by increasing(decreasing) the max_links by increment
*/
void resize_remap_vars(remapvars_t *rv, int increment)
{
  /*
    Input variables:
    int  increment  ! the number of links to add(subtract) to arrays
  */

  /*  Reallocate arrays at new size */

  rv->max_links += increment;

  if ( rv->max_links )
    {
      rv->src_cell_add = (int*) realloc(rv->src_cell_add, rv->max_links*sizeof(int));
      rv->tgt_cell_add = (int*) realloc(rv->tgt_cell_add, rv->max_links*sizeof(int));

      rv->wts = (double*) realloc(rv->wts, rv->num_wts*rv->max_links*sizeof(double));
    }

} /* resize_remap_vars */


/*
  -----------------------------------------------------------------------

  Performs the remapping based on weights computed elsewhere
     
  -----------------------------------------------------------------------
*/
void remap(double *restrict dst_array, double missval, long dst_size, long num_links, double *restrict map_wts, 
	   long num_wts, const int *restrict dst_add, const int *restrict src_add, const double *restrict src_array, 
	   const double *restrict src_grad1, const double *restrict src_grad2, const double *restrict src_grad3,
	   remaplink_t links)
{
  /*
    Input arrays:

    int *dst_add         ! destination address for each link
    int *src_add         ! source      address for each link

    int num_wts          ! num of weights used in remapping

    double *map_wts      ! remapping weights for each link

    double *src_array    ! array with source field to be remapped

    optional:

    double *src_grad1    ! gradient arrays on source grid necessary for
    double *src_grad2    ! higher-order remappings
    double *src_grad3

    output variables:

    double *dst_array    ! array for remapped field on destination grid
  */

  /* Local variables */
  long n;
  int iorder;
  extern int timer_remap;

  /* Check the order of the interpolation */

  if ( src_grad1 )
    iorder = 2;
  else
    iorder = 1;

  for ( n = 0; n < dst_size; ++n ) dst_array[n] = missval;

  if ( cdoTimer ) timer_start(timer_remap);

#ifdef SX
#pragma cdir nodep
#endif
  for ( n = 0; n < num_links; ++n ) dst_array[dst_add[n]] = 0.;

  if ( iorder == 1 )   /* First order remapping */
    {
      if ( links.option == TRUE )
	{
	  long j;
	  for ( j = 0; j < links.num_blks; ++j )
	    {
#ifdef SX
#pragma cdir nodep
#endif
	      for ( n = 0; n < links.num_links[j]; ++n )
		{
		  dst_array[links.dst_add[j][n]] += src_array[links.src_add[j][n]]*map_wts[num_wts*links.w_index[j][n]];
		}
	    }
	}
      else
	{
	  for ( n = 0; n < num_links; ++n )
	    {
	      /*
		printf("%5d %5d %5d %g # dst_add src_add n\n", dst_add[n], src_add[n], n, map_wts[num_wts*n]);
	      */
	      dst_array[dst_add[n]] += src_array[src_add[n]]*map_wts[num_wts*n];
	    }
	}
    }
  else                 /* Second order remapping */
    {
      if ( num_wts == 3 )
	{
	  for ( n = 0; n < num_links; ++n )
	    {
	      dst_array[dst_add[n]] += src_array[src_add[n]]*map_wts[3*n] +
                                       src_grad1[src_add[n]]*map_wts[3*n+1] +
                                       src_grad2[src_add[n]]*map_wts[3*n+2];
	    }
	}
      else if ( num_wts == 4 )
	{
      	  for ( n = 0; n < num_links; ++n )
	    {
              dst_array[dst_add[n]] += src_array[src_add[n]]*map_wts[4*n] +
                                       src_grad1[src_add[n]]*map_wts[4*n+1] +
                                       src_grad2[src_add[n]]*map_wts[4*n+2] +
                                       src_grad3[src_add[n]]*map_wts[4*n+3];
	    }
	}
    }

  if ( cdoTimer ) timer_stop(timer_remap);
}

static
long get_max_add(long num_links, long size, const int *restrict add)
{
  long n, i;
  long max_add;
  int *isum;

  isum = (int*) malloc(size*sizeof(int));
  memset(isum, 0, size*sizeof(int));

  for ( n = 0; n < num_links; ++n ) isum[add[n]]++;

  max_add = 0;
  for ( i = 0; i < size; ++i ) if ( isum[i] > max_add ) max_add = isum[i];
  free(isum);

  return (max_add);
}

static 
long binary_search_int(const int *array, long len, int value)
{       
  long low = 0, high = len - 1, midpoint = 0;
 
  while ( low <= high )
    {
      midpoint = low + (high - low)/2;      
 
      // check to see if value is equal to item in array
      if ( value == array[midpoint] ) return midpoint;

      if ( value < array[midpoint] )
	high = midpoint - 1;
      else
	low  = midpoint + 1;
    }
 
  // item was not found
  return -1L;
}

/*
  -----------------------------------------------------------------------

  Performs the remapping based on weights computed elsewhere
     
  -----------------------------------------------------------------------
*/
void remap_laf(double *restrict dst_array, double missval, long dst_size, long num_links, double *restrict map_wts,
	       long num_wts, const int *restrict dst_add, const int *restrict src_add, const double *restrict src_array)
{
  /*
    Input arrays:

    int *dst_add         ! destination address for each link
    int *src_add         ! source      address for each link

    int num_wts          ! num of weights used in remapping

    double *map_wts      ! remapping weights for each link

    double *src_array    ! array with source field to be remapped

    output variables:

    double *dst_array    ! array for remapped field on destination grid
  */

  /* Local variables */
  long i, n, k, ncls, imax;
  long max_cls;
  double wts;
  double *src_cls;
  double *src_wts;
#if defined(_OPENMP)
  double **src_cls2;
  double **src_wts2;
#endif

  for ( i = 0; i < dst_size; ++i ) dst_array[i] = missval;

  if ( num_links == 0 ) return;

  max_cls = get_max_add(num_links, dst_size, dst_add);

#if defined(_OPENMP)
  src_cls2 = (double **) malloc(ompNumThreads*sizeof(double *));
  src_wts2 = (double **) malloc(ompNumThreads*sizeof(double *));
  for ( i = 0; i < ompNumThreads; ++i )
    {
      src_cls2[i] = (double*) malloc(max_cls*sizeof(double));
      src_wts2[i] = (double*) malloc(max_cls*sizeof(double));
    }
#else
  src_cls = (double*) malloc(max_cls*sizeof(double));
  src_wts = (double*) malloc(max_cls*sizeof(double));
#endif

  for ( n = 0; n < num_links; ++n )
    if ( DBL_IS_EQUAL(dst_array[dst_add[n]], missval) ) dst_array[dst_add[n]] = ZERO;

#if defined(_OPENMP)
#pragma omp parallel for default(none) \
  shared(dst_size, src_cls2, src_wts2, num_links, dst_add, src_add, src_array, map_wts, num_wts, dst_array, max_cls)					\
  private(i, n, k, src_cls, src_wts, ncls, imax, wts) \
  schedule(dynamic,1)
#endif
  for ( i = 0; i < dst_size; ++i )
    {
#if defined(_OPENMP)
      int ompthID = cdo_omp_get_thread_num();
      src_cls = src_cls2[ompthID];
      src_wts = src_wts2[ompthID];
#endif
      memset(src_cls, 0, max_cls*sizeof(double));
      memset(src_wts, 0, max_cls*sizeof(double));
      /*
      ncls = 0;
      for ( n = 0; n < num_links; n++ )
	{
	  if ( i == dst_add[n] )
	    {
	      for ( k = 0; k < ncls; k++ )
		if ( IS_EQUAL(src_array[src_add[n]], src_cls[k]) ) break;
	      
	      if ( k == ncls )
		{
		  src_cls[k] = src_array[src_add[n]];
		  ncls++;
		}
	      
	      src_wts[k] += map_wts[num_wts*n];
	    }
	}
      */
      /* only for sorted dst_add! */
      {
      long min_add = 1, max_add = 0;

      n = binary_search_int(dst_add, num_links, (int)i);

      if ( n >= 0 && n < num_links )
	{
	  min_add = n;
	  
	  for ( n = min_add+1; n < num_links; ++n )
	    if ( i != dst_add[n] ) break;

	  max_add = n;

	  for ( n = min_add; n > 0; --n )
	    if ( i != dst_add[n-1] ) break;

	  min_add = n;
	}

      ncls = 0;
      for ( n = min_add; n < max_add; ++n )
	{
	  for ( k = 0; k < ncls; ++k )
	    if ( IS_EQUAL(src_array[src_add[n]], src_cls[k]) ) break;
	      
	  if ( k == ncls )
	    {
	      src_cls[k] = src_array[src_add[n]];
	      ncls++;
	    }
	      
	  src_wts[k] += map_wts[num_wts*n];
	}
      }
      
      if ( ncls )
	{
	  imax = 0;
	  wts = src_wts[0];
	  for ( k = 1; k < ncls; ++k )
	    {
	      if ( src_wts[k] > wts )
		{
		  wts  = src_wts[k];
		  imax = k;
		}
	    }

	  dst_array[i] = src_cls[imax];
	}
    }

#if defined(_OPENMP)
  for ( i = 0; i < ompNumThreads; ++i )
    {
      free(src_cls2[i]);
      free(src_wts2[i]);
    }

  free(src_cls2);
  free(src_wts2);
#else
  free(src_cls);
  free(src_wts);
#endif
}

/*
  -----------------------------------------------------------------------

  Performs the remapping based on weights computed elsewhere
     
  -----------------------------------------------------------------------
*/
void remap_sum(double *restrict dst_array, double missval, long dst_size, long num_links, double *restrict map_wts,
	       long num_wts, const int *restrict dst_add, const int *restrict src_add, const double *restrict src_array)
{
  /*
    Input arrays:

    int *dst_add         ! destination address for each link
    int *src_add         ! source      address for each link

    int num_wts          ! num of weights used in remapping

    double *map_wts      ! remapping weights for each link

    double *src_array    ! array with source field to be remapped

    output variables:

    double *dst_array    ! array for remapped field on destination grid
  */
  /* Local variables */
  long n;

  for ( n = 0; n < dst_size; ++n ) dst_array[n] = missval;

#ifdef SX
#pragma cdir nodep
#endif
  for ( n = 0; n < num_links; ++n )
    if ( DBL_IS_EQUAL(dst_array[dst_add[n]], missval) ) dst_array[dst_add[n]] = ZERO;

  for ( n = 0; n < num_links; ++n )
    {
      /*
	printf("%5d %5d %5d %g # dst_add src_add n\n", dst_add[n], src_add[n], n, map_wts[num_wts*n]);
      */
      //dst_array[dst_add[n]] += src_array[src_add[n]]*map_wts[num_wts*n];
      dst_array[dst_add[n]] += src_array[src_add[n]]*map_wts[num_wts*n];
      printf("%ld %d %d %g %g %g\n", n, dst_add[n], src_add[n],
	     src_array[src_add[n]], map_wts[num_wts*n], dst_array[dst_add[n]]);
    }
}

/*****************************************************************************/

void remap_stat(int remap_order, remapgrid_t src_grid, remapgrid_t tgt_grid, remapvars_t rv, const double *restrict array1, 
		const double *restrict array2, double missval)
{
  long n, ns, i;
  long idiff, imax, imin, icount;
  int *tgt_count;
  double minval, maxval, sum;
	  
  if ( remap_order == 2 )
    cdoPrint("Second order mapping from grid1 to grid2:");
  else
    cdoPrint("First order mapping from grid1 to grid2:");
  cdoPrint("----------------------------------------");

  ns = 0;
  sum = 0;
  minval =  DBL_MAX;
  maxval = -DBL_MAX;
  for ( n = 0; n < src_grid.size; ++n )
    {
      if ( !DBL_IS_EQUAL(array1[n], missval) )
	{
	  if ( array1[n] < minval ) minval = array1[n];
	  if ( array1[n] > maxval ) maxval = array1[n];
	  sum += array1[n];
	  ns++;
	}
    }
  if ( ns > 0 ) sum /= ns;
  cdoPrint("Grid1 min,mean,max: %g %g %g", minval, sum, maxval);

  ns = 0;
  sum = 0;
  minval =  DBL_MAX;
  maxval = -DBL_MAX;
  for ( n = 0; n < tgt_grid.size; ++n )
    {
      if ( !DBL_IS_EQUAL(array2[n], missval) )
	{
	  if ( array2[n] < minval ) minval = array2[n];
	  if ( array2[n] > maxval ) maxval = array2[n];
	  sum += array2[n];
	  ns++;
	}
    }
  if ( ns > 0 ) sum /= ns;
  cdoPrint("Grid2 min,mean,max: %g %g %g", minval, sum, maxval);

  /* Conservation Test */

  if ( src_grid.cell_area )
    {
      cdoPrint("Conservation:");
      sum = 0;
      for ( n = 0; n < src_grid.size; ++n )
	if ( !DBL_IS_EQUAL(array1[n], missval) )
	  sum += array1[n]*src_grid.cell_area[n]*src_grid.cell_frac[n];
      cdoPrint("Grid1 Integral = %g", sum);

      sum = 0;
      for ( n = 0; n < tgt_grid.size; ++n )
	if ( !DBL_IS_EQUAL(array2[n], missval) )
	  sum += array2[n]*tgt_grid.cell_area[n]*tgt_grid.cell_frac[n];
      cdoPrint("Grid2 Integral = %g", sum);
      /*
      for ( n = 0; n < src_grid.size; n++ )
       fprintf(stderr, "1 %d %g %g %g\n", n, array1[n], src_grid.cell_area[n], src_grid.cell_frac[n]);
      for ( n = 0; n < tgt_grid.size; n++ )
	fprintf(stderr, "2 %d %g %g %g\n", n, array2[n], tgt_grid.cell_area[n], tgt_grid.cell_frac[n]);
      */
    }

  cdoPrint("number of sparse matrix entries %d", rv.num_links);
  cdoPrint("total number of dest cells %d", tgt_grid.size);

  tgt_count = (int*) malloc(tgt_grid.size*sizeof(int));

  for ( n = 0; n < tgt_grid.size; ++n ) tgt_count[n] = 0;

#if defined(SX)
#pragma vdir nodep
#endif
  for ( n = 0; n < rv.num_links; ++n ) tgt_count[rv.tgt_cell_add[n]]++;

  imin = INT_MAX;
  imax = INT_MIN;
  for ( n = 0; n < tgt_grid.size; ++n )
    {
      if ( tgt_count[n] > 0 )
	if ( tgt_count[n] < imin ) imin = tgt_count[n];
      if ( tgt_count[n] > imax ) imax = tgt_count[n];
    }

  idiff =  (imax - imin)/10 + 1;
  icount = 0;
  for ( i = 0; i < tgt_grid.size; ++i )
    if ( tgt_count[i] > 0 ) icount++;

  cdoPrint("number of cells participating in remap %d", icount);

  if ( icount )
    {
      cdoPrint("min no of entries/row = %d", imin);
      cdoPrint("max no of entries/row = %d", imax);

      imax = imin + idiff;
      for ( n = 0; n < 10; ++n )
	{
	  icount = 0;
	  for ( i = 0; i < tgt_grid.size; ++i )
	    if ( tgt_count[i] >= imin && tgt_count[i] < imax ) icount++;

	  if ( icount )
	    cdoPrint("num of rows with entries between %d - %d  %d", imin, imax-1, icount);

	  imin = imin + idiff;
	  imax = imax + idiff;
	}
    }

  free(tgt_count);

  if ( rv.sort_add )
    cdoPrint("Sparse matrix entries are explicitly sorted.");

} /* remap_stat */

/*****************************************************************************/

void remap_gradients(remapgrid_t grid, const double *restrict array, double *restrict grad_lat,
		     double *restrict grad_lon, double *restrict grad_latlon)
{
  long n, nx, ny, grid_size;
  long i, j, ip1, im1, jp1, jm1, in, is, ie, iw, ine, inw, ise, isw;
  double delew, delns;
  double grad_lat_zero, grad_lon_zero;

  if ( grid.rank != 2 )
    cdoAbort("Internal problem (remap_gradients), grid rank = %d!", grid.rank);

  grid_size = grid.size;
  nx = grid.dims[0];
  ny = grid.dims[1];

#if defined(_OPENMP)
#pragma omp parallel for default(none)        \
  shared(grid_size, grad_lat, grad_lon, grad_latlon, grid, nx, ny, array) \
  private(n, i, j, ip1, im1, jp1, jm1, in, is, ie, iw, ine, inw, ise, isw, delew, delns, grad_lat_zero, grad_lon_zero)
#endif
  for ( n = 0; n < grid_size; ++n )
    {
      grad_lat[n] = ZERO;
      grad_lon[n] = ZERO;
      grad_latlon[n] = ZERO;

      if ( grid.mask[n] )
	{
	  delew = HALF;
	  delns = HALF;

	  j = n/nx + 1;
	  i = n - (j-1)*nx + 1;

	  ip1 = i + 1;
	  im1 = i - 1;
	  jp1 = j + 1;
	  jm1 = j - 1;

	  if ( ip1 > nx ) ip1 = ip1 - nx;
	  if ( im1 <  1 ) im1 = nx;
	  if ( jp1 > ny )
	    {
              jp1 = j;
              delns = ONE;
            }
	  if ( jm1 < 1 )
	    {
              jm1 = j;
              delns = ONE;
            }

	  in  = (jp1-1)*nx + i - 1;
	  is  = (jm1-1)*nx + i - 1;
	  ie  = (j  -1)*nx + ip1 - 1;
	  iw  = (j  -1)*nx + im1 - 1;

	  ine = (jp1-1)*nx + ip1 - 1;
	  inw = (jp1-1)*nx + im1 - 1;
	  ise = (jm1-1)*nx + ip1 - 1;
	  isw = (jm1-1)*nx + im1 - 1;

	  /* Compute i-gradient */

	  if ( ! grid.mask[ie] )
	    {
              ie = n;
              delew = ONE;
            }
	  if ( ! grid.mask[iw] )
	    {
              iw = n;
              delew = ONE;
            }
 
	  grad_lat[n] = delew*(array[ie] - array[iw]);

	  /* Compute j-gradient */

	  if ( ! grid.mask[in] )
	    {
              in = n;
              delns = ONE;
            }
	  if ( ! grid.mask[is] )
	    {
              is = n;
              delns = ONE;
            }
 
	  grad_lon[n] = delns*(array[in] - array[is]);

	  /* Compute ij-gradient */

	  delew = HALF;
	  if ( jp1 == j || jm1 == j )
	    delns = ONE;
	  else 
	    delns = HALF;

	  if ( ! grid.mask[ine] )
	    {
              if ( in != n )
		{
		  ine = in;
		  delew = ONE;
		}
              else if ( ie != n )
		{
		  ine = ie;
		  inw = iw;
		  if ( inw == n ) delew = ONE;
		  delns = ONE;
		}
              else
		{
		  ine = n;
		  inw = iw;
		  delew = ONE;
		  delns = ONE;
		}
	    }

	  if ( ! grid.mask[inw] )
	    {
              if ( in != n )
		{
		  inw = in;
		  delew = ONE;
		}
              else if ( iw != n )
		{
		  inw = iw;
		  ine = ie;
		  if ( ie == n ) delew = ONE;
		  delns = ONE;
		}
              else
		{
		  inw = n;
		  ine = ie;
		  delew = ONE;
		  delns = ONE;
		}
	    }

	  grad_lat_zero = delew*(array[ine] - array[inw]);

	  if ( ! grid.mask[ise] )
	    {
              if ( is != n )
		{
		  ise = is;
		  delew = ONE;
		}
              else if ( ie != n )
		{
		  ise = ie;
		  isw = iw;
		  if ( isw == n ) delew = ONE;
		  delns = ONE;
		}
              else
		{
		  ise = n;
		  isw = iw;
		  delew = ONE;
		  delns = ONE;
		}
	    }

	  if ( ! grid.mask[isw] )
	    {
              if ( is != n )
		{
		  isw = is;
		  delew = ONE;
		}
              else if ( iw != n )
		{
		  isw = iw;
		  ise = ie;
		  if ( ie == n ) delew = ONE;
		  delns = ONE;
		}
              else
		{
		  isw = n;
		  ise = ie;
		  delew = ONE;
		  delns = ONE;
		}
	    }

	  grad_lon_zero = delew*(array[ise] - array[isw]);

	  grad_latlon[n] = delns*(grad_lat_zero - grad_lon_zero);
	}
    }
} /* remap_gradients */

/*****************************************************************************/

void reorder_links(remapvars_t *rv)
{
  long j, nval = 0, num_blks = 0;
  long lastval;
  long nlinks;
  long max_links = 0;
  long n;
  long num_links;

  num_links = rv->num_links;

  printf("reorder_links\n");
  printf("  num_links %ld\n", num_links);
  rv->links.option = TRUE;

  lastval = -1;
  for ( n = 0; n < num_links; n++ )
    {
      if ( rv->tgt_cell_add[n] == lastval ) nval++;
      else
	{
	  if ( nval > num_blks ) num_blks = nval;
	  nval = 1;
	  max_links++;
	  lastval = rv->tgt_cell_add[n];
	}
    }

  if ( num_blks )
    {
      rv->links.max_links = max_links;
      rv->links.num_blks  = num_blks;

      printf("num_links %ld  max_links %ld  num_blks %ld\n", rv->num_links, max_links, num_blks);

      rv->links.num_links = (int*) malloc(num_blks*sizeof(int));
      rv->links.dst_add   = (int **) malloc(num_blks*sizeof(int *));
      rv->links.src_add   = (int **) malloc(num_blks*sizeof(int *));
      rv->links.w_index   = (int **) malloc(num_blks*sizeof(int *));
    }

  for ( j = 0; j < num_blks; j++ )
    {
      rv->links.dst_add[j] = (int*) malloc(max_links*sizeof(int));
      rv->links.src_add[j] = (int*) malloc(max_links*sizeof(int));
      rv->links.w_index[j] = (int*) malloc(max_links*sizeof(int));
    }

  for ( j = 0; j < num_blks; j++ )
    {
      nval = 0;
      lastval = -1;
      nlinks = 0;

      for ( n = 0; n < num_links; n++ )
	{
	  if ( rv->tgt_cell_add[n] == lastval ) nval++;
	  else
	    {
	      nval = 1;
	      lastval = rv->tgt_cell_add[n];
	    }
	  
	  if ( nval == j+1 )
	    {
	      rv->links.dst_add[j][nlinks] = rv->tgt_cell_add[n];
	      rv->links.src_add[j][nlinks] = rv->src_cell_add[n];
	      rv->links.w_index[j][nlinks] = n;
	      nlinks++;
	    }
	}

      rv->links.num_links[j] = nlinks;
      printf("loop %ld  nlinks %ld\n", j+1, nlinks);
    }
}
