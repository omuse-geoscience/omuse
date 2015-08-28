#include "cdo.h"
#include "cdo_int.h"
#include "grid.h"
#include "remap.h"
#include "remap_store_link.h"


int rect_grid_search2(long *imin, long *imax, double xmin, double xmax, long nxm, const double *restrict xm);

static
long get_srch_cells_reg2d(const int *restrict src_grid_dims, 
			  const double *restrict src_corner_lat, const double *restrict src_corner_lon,
			  const double *restrict tgt_cell_bound_box, int *srch_add)
{
  long nx = src_grid_dims[0];
  long ny = src_grid_dims[1];
  long num_srch_cells;  /* num cells in restricted search arrays   */
  int lfound;
  long nxp1, nyp1;
  double src_lon_min, src_lon_max;
  int debug = 0;

  nxp1 = nx+1;
  nyp1 = ny+1;

  src_lon_min = src_corner_lon[0];
  src_lon_max = src_corner_lon[nx];

  double bound_lon1, bound_lon2;

  num_srch_cells = 0;

  long imin = nxp1, imax = -1, jmin = nyp1, jmax = -1;
  long im, jm;

  lfound = rect_grid_search2(&jmin, &jmax, tgt_cell_bound_box[0], tgt_cell_bound_box[1], nyp1, src_corner_lat);
  // if ( jmin > 0 ) jmin--;
  // if ( jmax < (ny-2) ) jmax++;
  bound_lon1 = tgt_cell_bound_box[2];
  bound_lon2 = tgt_cell_bound_box[3];
  if ( bound_lon1 <= src_lon_max && bound_lon2 >= src_lon_min )
    {
      if ( debug ) printf("  b1 %g %g\n", bound_lon1*RAD2DEG, bound_lon2*RAD2DEG);
      if ( bound_lon1 < src_lon_min && bound_lon2 > src_lon_min ) bound_lon1 = src_lon_min;
      if ( bound_lon2 > src_lon_max && bound_lon1 < src_lon_max ) bound_lon2 = src_lon_max;
      lfound = rect_grid_search2(&imin, &imax, bound_lon1, bound_lon2, nxp1, src_corner_lon);
      if ( lfound )
	{
	  if ( debug )
	    printf("   %g %g imin %ld  imax %ld  jmin %ld jmax %ld\n", RAD2DEG*src_corner_lon[imin], RAD2DEG*src_corner_lon[imax+1], imin, imax, jmin, jmax);
	  for ( jm = jmin; jm <= jmax; ++jm )
	    for ( im = imin; im <= imax; ++im )
	      srch_add[num_srch_cells++] = jm*nx + im;
	}
    }

  bound_lon1 = tgt_cell_bound_box[2];
  bound_lon2 = tgt_cell_bound_box[3];
  if ( bound_lon1 <= src_lon_min && bound_lon2 >= src_lon_min )
    {
      long imin2 = nxp1, imax2 = -1;
      bound_lon1 += 2*M_PI;
      bound_lon2 += 2*M_PI;
      if ( debug ) printf("  b2 %g %g\n", bound_lon1*RAD2DEG, bound_lon2*RAD2DEG);
      if ( bound_lon1 < src_lon_min && bound_lon2 > src_lon_min ) bound_lon1 = src_lon_min;
      if ( bound_lon2 > src_lon_max && bound_lon1 < src_lon_max ) bound_lon2 = src_lon_max;
      lfound = rect_grid_search2(&imin2, &imax2, bound_lon1, bound_lon2, nxp1, src_corner_lon);
      if ( lfound )
	{
	  if ( imax != -1 && imin2 <= imax ) imin2 = imax+1;
	  if ( imax != -1 && imax2 <= imax ) imax2 = imax+1;
	  if ( imin2 >= 0 && imax2 < nxp1 )
	    {
	      if ( debug )
		printf("   %g %g imin %ld  imax %ld  jmin %ld jmax %ld\n", RAD2DEG*src_corner_lon[imin2], RAD2DEG*src_corner_lon[imax2+1], imin2, imax2, jmin, jmax);
	      for ( jm = jmin; jm <= jmax; ++jm )
		for ( im = imin2; im <= imax2; ++im )
		  srch_add[num_srch_cells++] = jm*nx + im;
	    }
	}
    }

  bound_lon1 = tgt_cell_bound_box[2];
  bound_lon2 = tgt_cell_bound_box[3];
  if ( bound_lon1 <= src_lon_max && bound_lon2 >= src_lon_max )
    {
      long imin3 = nxp1, imax3 = -1;
      bound_lon1 -= 2*M_PI;
      bound_lon2 -= 2*M_PI;
      if ( debug ) printf("  b3 %g %g\n", bound_lon1*RAD2DEG, bound_lon2*RAD2DEG);
      if ( bound_lon1 < src_lon_min && bound_lon2 > src_lon_min ) bound_lon1 = src_lon_min;
      if ( bound_lon2 > src_lon_max && bound_lon1 < src_lon_max ) bound_lon2 = src_lon_max;
      lfound = rect_grid_search2(&imin3, &imax3, bound_lon1, bound_lon2, nxp1, src_corner_lon);
      if ( lfound )
	{
	  if ( imin != nxp1 && imin3 >= imin ) imin3 = imin-1;
	  if ( imax != nxp1 && imax3 >= imin ) imax3 = imin-1;
	  if ( imin3 >= 0 && imin3 < nxp1 )
	    {
	      if ( debug )
		printf("   %g %g imin %ld  imax %ld  jmin %ld jmax %ld\n", RAD2DEG*src_corner_lon[imin3], RAD2DEG*src_corner_lon[imax3+1], imin3, imax3, jmin, jmax);
	      for ( jm = jmin; jm <= jmax; ++jm )
		for ( im = imin3; im <= imax3; ++im )
		  srch_add[num_srch_cells++] = jm*nx + im;
	    }
	}
    }

  if ( debug ) printf(" -> num_srch_cells: %ld\n", num_srch_cells);

  return (num_srch_cells);
}

static
void restrict_boundbox(const double *restrict grid_bound_box, double *restrict bound_box)
{
  if ( bound_box[0] < grid_bound_box[0] && bound_box[1] > grid_bound_box[0] ) bound_box[0] = grid_bound_box[0];
  if ( bound_box[1] > grid_bound_box[1] && bound_box[0] < grid_bound_box[1] ) bound_box[1] = grid_bound_box[1];

  if ( bound_box[2] >= grid_bound_box[3] && (bound_box[3]-2*M_PI) > grid_bound_box[2] ) { bound_box[2] -= 2*M_PI; bound_box[3] -= 2*M_PI; }
  if ( bound_box[3] <= grid_bound_box[2] && (bound_box[2]-2*M_PI) < grid_bound_box[3] ) { bound_box[2] += 2*M_PI; bound_box[3] += 2*M_PI; }
  //  if ( bound_box[2] < grid_bound_box[2] && bound_box[3] > grid_bound_box[2] ) bound_box[2] = grid_bound_box[2];
  //  if ( bound_box[3] > grid_bound_box[3] && bound_box[2] < grid_bound_box[3] ) bound_box[3] = grid_bound_box[3];
}

static
void boundbox_from_corners_reg2d(long grid_add, const int *restrict grid_dims, const double *restrict corner_lon,
				 const double *restrict corner_lat, double *restrict bound_box)
{
  long nx = grid_dims[0];
  long iy = grid_add/nx;
  long ix = grid_add - iy*nx;

  double clat1 = corner_lat[iy  ];
  double clat2 = corner_lat[iy+1];

  if ( clat2 > clat1 )
    {
      bound_box[0] = clat1;
      bound_box[1] = clat2;
    }
  else
    {
      bound_box[0] = clat2;
      bound_box[1] = clat1;
    }

  bound_box[2] = corner_lon[ix  ];
  bound_box[3] = corner_lon[ix+1];
}

static
void boundbox_from_corners1(long ic, long nc, const double *restrict corner_lon,
			    const double *restrict corner_lat, double *restrict bound_box)
{
  long inc, j;
  double clon, clat;

  inc = ic*nc;

  clat = corner_lat[inc];
  clon = corner_lon[inc];

  bound_box[0] = clat;
  bound_box[1] = clat;
  bound_box[2] = clon;
  bound_box[3] = clon;

  for ( j = 1; j < nc; ++j )
    {
      clat = corner_lat[inc+j];
      clon = corner_lon[inc+j];

      if ( clat < bound_box[0] ) bound_box[0] = clat;
      if ( clat > bound_box[1] ) bound_box[1] = clat;
      if ( clon < bound_box[2] ) bound_box[2] = clon;
      if ( clon > bound_box[3] ) bound_box[3] = clon;
    }

  if ( fabs(bound_box[3] - bound_box[2]) > PI )
    {
      bound_box[2] = 0;
      bound_box[3] = PI2;
    }

  /*
  double dlon = fabs(bound_box[3] - bound_box[2]);

  if ( dlon > PI )
    {
      if ( bound_box[3] > bound_box[2] && (bound_box[3]-PI2) < 0. )
	{
	  double tmp = bound_box[2];
	  bound_box[2] = bound_box[3] - PI2;
	  bound_box[3] = tmp;
	}
      else
	{
	  bound_box[2] = 0;
	  bound_box[3] = PI2;
	}
    }
  */
}

static
void boundbox_from_corners1r(long ic, long nc, const double *restrict corner_lon,
			     const double *restrict corner_lat, restr_t *restrict bound_box)
{
  long inc, j;
  restr_t clon, clat;

  inc = ic*nc;

  clat = RESTR_SCALE(corner_lat[inc]);
  clon = RESTR_SCALE(corner_lon[inc]);

  bound_box[0] = clat;
  bound_box[1] = clat;
  bound_box[2] = clon;
  bound_box[3] = clon;

  for ( j = 1; j < nc; ++j )
    {
      clat = RESTR_SCALE(corner_lat[inc+j]);
      clon = RESTR_SCALE(corner_lon[inc+j]);

      if ( clat < bound_box[0] ) bound_box[0] = clat;
      if ( clat > bound_box[1] ) bound_box[1] = clat;
      if ( clon < bound_box[2] ) bound_box[2] = clon;
      if ( clon > bound_box[3] ) bound_box[3] = clon;
    }

  if ( RESTR_ABS(bound_box[3] - bound_box[2]) > RESTR_SCALE(PI) )
    {
      bound_box[2] = 0;
      bound_box[3] = RESTR_SCALE(PI2);
    }
  /*
  if ( RESTR_ABS(bound_box[3] - bound_box[2]) > RESTR_SCALE(PI) )
    {
      if ( bound_box[3] > bound_box[2] && (bound_box[3]-RESTR_SCALE(PI2)) < RESTR_SCALE(0.) )
	{
	  restr_t tmp = bound_box[2];
	  bound_box[2] = bound_box[3] - RESTR_SCALE(PI2);
	  bound_box[3] = tmp;
	}
    }
  */
}

//#if defined(HAVE_LIBYAC)
#include "clipping/clipping.h"
#include "clipping/area.h"
#include "clipping/geometry.h"

static
double gridcell_area(struct grid_cell cell)
{
  return yac_huiliers_area(cell);
}

static
void cdo_compute_overlap_areas(unsigned N,
			       struct grid_cell *overlap_buffer,
			       struct grid_cell *source_cells,
			       struct grid_cell  target_cell,
			       double *partial_areas)
{
  /* Do the clipping and get the cell for the overlapping area */

  yac_cell_clipping(N, source_cells, target_cell, overlap_buffer);

  /* Get the partial areas for the overlapping regions */

  for ( unsigned n = 0; n < N; n++ )
    {
      partial_areas[n] = gridcell_area(overlap_buffer[n]);
    }

#ifdef VERBOSE
  for ( unsigned n = 0; n < N; n++ )
    printf("overlap area : %lf\n", partial_areas[n]);
#endif
}

static double const tol = 1.0e-12;

enum cell_type {
  UNDEF_CELL,
  LON_LAT_CELL,
  LAT_CELL,
  GREAT_CIRCLE_CELL,
  MIXED_CELL
};
/*
static enum cell_type get_cell_type(struct grid_cell target_cell) {

  int count_lat_edges = 0, count_great_circle_edges = 0;

   if ((target_cell.num_corners == 4) &&
       ((target_cell.edge_type[0] == LAT_CIRCLE &&
         target_cell.edge_type[1] == LON_CIRCLE &&
         target_cell.edge_type[2] == LAT_CIRCLE &&
         target_cell.edge_type[3] == LON_CIRCLE) ||
        (target_cell.edge_type[0] == LON_CIRCLE &&
         target_cell.edge_type[1] == LAT_CIRCLE &&
         target_cell.edge_type[2] == LON_CIRCLE &&
         target_cell.edge_type[3] == LAT_CIRCLE)))
      return LON_LAT_CELL;
   else
      for (unsigned i = 0; i < target_cell.num_corners; ++i)
         if (target_cell.edge_type[i] == LON_CIRCLE ||
             target_cell.edge_type[i] == GREAT_CIRCLE)
            count_great_circle_edges++;
         else
            count_lat_edges++;

   if (count_lat_edges && count_great_circle_edges)
      return MIXED_CELL;
   else if (count_lat_edges)
      return LAT_CELL;
   else
      return GREAT_CIRCLE_CELL;
}
*/
static
void cdo_compute_concave_overlap_areas(unsigned N,
				       struct grid_cell *overlap_buffer,
				       struct grid_cell *source_cell,
				       struct grid_cell  target_cell,
				       double target_node_x,
				       double target_node_y,
				       double * partial_areas)
{
  /*
  enum cell_type target_cell_type = UNDEF_CELL;

  if ( target_cell.num_corners > 3 )
    target_cell_type = get_cell_type(target_cell);

  if ( target_cell.num_corners < 4 || target_cell_type == LON_LAT_CELL )
    {
      cdo_compute_overlap_areas(N, overlap_buffer, source_cell, target_cell, partial_areas);
      return;
    }

  if ( target_node_x == NULL || target_node_y == NULL )
    cdoAbort("Internal problem (cdo_compute_concave_overlap_areas): missing target point coordinates!");
  */
  /*
  struct grid_cell target_partial_cell =
    {.coordinates_x   = (double[3]){-1, -1, -1},
     .coordinates_y   = (double[3]){-1, -1, -1},
     .coordinates_xyz = (double[3*3]){-1, -1, -1},
     .edge_type       = (enum yac_edge_type[3]) {GREAT_CIRCLE, GREAT_CIRCLE, GREAT_CIRCLE},
     .num_corners     = 3};
  */
  double coordinates_x[3] = {-1, -1, -1};
  double coordinates_y[3] = {-1, -1, -1};
  double coordinates_xyz[9] = {-1, -1, -1};
  enum yac_edge_type edge_types[3] = {GREAT_CIRCLE, GREAT_CIRCLE, GREAT_CIRCLE};
  struct grid_cell target_partial_cell =
    {.coordinates_x   = coordinates_x,
     .coordinates_y   = coordinates_y,
     .coordinates_xyz = coordinates_xyz,
     .edge_type       = edge_types,
     .num_corners     = 3};

  /* Do the clipping and get the cell for the overlapping area */

  for ( unsigned n = 0; n < N; n++) partial_areas[n] = 0.0;

  // common node point to all partial target cells
  target_partial_cell.coordinates_x[0] = target_node_x;
  target_partial_cell.coordinates_y[0] = target_node_y;

  LLtoXYZ ( target_node_x, target_node_y, target_partial_cell.coordinates_xyz );

  for ( unsigned num_corners = 0; num_corners < target_cell.num_corners; ++num_corners )
    {
      unsigned corner_a = num_corners;
      unsigned corner_b = (num_corners+1)%target_cell.num_corners;

      // skip clipping and area calculation for degenerated triangles
      //
      // If this is not sufficient, instead we can try something like:
      //
      //     struct point_list target_list
      //     init_point_list(&target_list);
      //     generate_point_list(&target_list, target_cell);
      //     struct grid_cell temp_target_cell;
      //     generate_overlap_cell(target_list, temp_target_cell);
      //     free_point_list(&target_list);
      //
      // and use temp_target_cell for triangulation.
      //
      // Compared to the if statement below the alternative seems
      // to be quite costly.

      if ( ( ( fabs(target_cell.coordinates_xyz[0+3*corner_a]-target_cell.coordinates_xyz[0+3*corner_b]) < tol ) &&
	     ( fabs(target_cell.coordinates_xyz[1+3*corner_a]-target_cell.coordinates_xyz[1+3*corner_b]) < tol ) &&
	     ( fabs(target_cell.coordinates_xyz[2+3*corner_a]-target_cell.coordinates_xyz[2+3*corner_b]) < tol ) ) ||
	   ( ( fabs(target_cell.coordinates_xyz[0+3*corner_a]-target_partial_cell.coordinates_xyz[0]) < tol    ) &&
	     ( fabs(target_cell.coordinates_xyz[1+3*corner_a]-target_partial_cell.coordinates_xyz[1]) < tol    ) &&
	     ( fabs(target_cell.coordinates_xyz[2+3*corner_a]-target_partial_cell.coordinates_xyz[2]) < tol    ) ) ||
	   ( ( fabs(target_cell.coordinates_xyz[0+3*corner_b]-target_partial_cell.coordinates_xyz[0]) < tol    ) &&
	     ( fabs(target_cell.coordinates_xyz[1+3*corner_b]-target_partial_cell.coordinates_xyz[1]) < tol    ) &&
	     ( fabs(target_cell.coordinates_xyz[2+3*corner_b]-target_partial_cell.coordinates_xyz[2]) < tol    ) ) )
	continue;

      target_partial_cell.coordinates_x[1] = target_cell.coordinates_x[corner_a];
      target_partial_cell.coordinates_y[1] = target_cell.coordinates_y[corner_a];
      target_partial_cell.coordinates_x[2] = target_cell.coordinates_x[corner_b];
      target_partial_cell.coordinates_y[2] = target_cell.coordinates_y[corner_b];

      target_partial_cell.coordinates_xyz[0+3*1] = target_cell.coordinates_xyz[0+3*corner_a];
      target_partial_cell.coordinates_xyz[1+3*1] = target_cell.coordinates_xyz[1+3*corner_a];
      target_partial_cell.coordinates_xyz[2+3*1] = target_cell.coordinates_xyz[2+3*corner_a];
      target_partial_cell.coordinates_xyz[0+3*2] = target_cell.coordinates_xyz[0+3*corner_b];
      target_partial_cell.coordinates_xyz[1+3*2] = target_cell.coordinates_xyz[1+3*corner_b];
      target_partial_cell.coordinates_xyz[2+3*2] = target_cell.coordinates_xyz[2+3*corner_b];

      yac_cell_clipping(N, source_cell, target_partial_cell, overlap_buffer);

      /* Get the partial areas for the overlapping regions as sum over the partial target cells. */

      for (unsigned n = 0; n < N; n++)
	{
	  partial_areas[n] += gridcell_area(overlap_buffer[n]);
	  // we cannot use pole_area because it is rather inaccurate for great circle
	  // edges that are nearly circles of longitude
	  //partial_areas[n] = pole_area (overlap_buffer[n]);
	}
    }

#ifdef VERBOSE
  for (unsigned n = 0; n < N; n++)
    printf("overlap area %i: %lf \n", n, partial_areas[n]);
#endif
}

//#endif

static
int get_lonlat_circle_index(remapgrid_t *remap_grid)
{
  int lonlat_circle_index = -1;

  if ( remap_grid->num_cell_corners == 4 )
    {
      if ( remap_grid->remap_grid_type == REMAP_GRID_TYPE_REG2D )
	{
	  lonlat_circle_index = 1;
 	}
      else
	{
	  const double* cell_corner_lon = remap_grid->cell_corner_lon;
	  const double* cell_corner_lat = remap_grid->cell_corner_lat;
	  int gridsize = remap_grid->size;
	  int num_i = 0, num_eq0 = 0, num_eq1 = 0;
	  int iadd = gridsize/3-1;

	  if ( iadd == 0 ) iadd++;

	  for ( int i = 0; i < gridsize; i += iadd )
	    {
	      num_i++;

	      if ( IS_EQUAL(cell_corner_lon[i*4+1], cell_corner_lon[i*4+2]) &&
		   IS_EQUAL(cell_corner_lon[i*4+3], cell_corner_lon[i*4+0]) &&
		   IS_EQUAL(cell_corner_lat[i*4+0], cell_corner_lat[i*4+1]) &&
		   IS_EQUAL(cell_corner_lat[i*4+2], cell_corner_lat[i*4+3]) )
		{  
		  num_eq1++;
		}
	      else if ( IS_EQUAL(cell_corner_lon[i*4+0], cell_corner_lon[i*4+1]) &&
			IS_EQUAL(cell_corner_lon[i*4+2], cell_corner_lon[i*4+3]) &&
			IS_EQUAL(cell_corner_lat[i*4+1], cell_corner_lat[i*4+2]) &&
			IS_EQUAL(cell_corner_lat[i*4+3], cell_corner_lat[i*4+0]) )
		{
		  num_eq0++;
		}
	    }

	  if ( num_i == num_eq1 ) lonlat_circle_index = 1;
	  if ( num_i == num_eq0 ) lonlat_circle_index = 0;	      
	}
    }

  //printf("lonlat_circle_index %d\n", lonlat_circle_index);

  return (lonlat_circle_index);
}



static
void normalize_weights(remapgrid_t *tgt_grid, remapvars_t *rv)
{
  /* Include centroids in weights and normalize using destination area if requested */
  long n;
  long num_wts = rv->num_wts;
  long num_links = rv->num_links;
  long tgt_cell_add;       /* current linear address for target grid cell   */
  double norm_factor = 0;  /* factor for normalizing wts */

  if ( rv->norm_opt == NORM_OPT_DESTAREA )
    {
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(_OPENMP)
#pragma omp parallel for default(none) \
  shared(num_wts, num_links, rv, tgt_grid)	\
  private(n, tgt_cell_add, norm_factor)
#endif
      for ( n = 0; n < num_links; ++n )
	{
	  tgt_cell_add = rv->tgt_cell_add[n];

          if ( IS_NOT_EQUAL(tgt_grid->cell_area[tgt_cell_add], 0) )
	    norm_factor = ONE/tgt_grid->cell_area[tgt_cell_add];
          else
            norm_factor = ZERO;

	  rv->wts[n*num_wts] *= norm_factor;
	}
    }
  else if ( rv->norm_opt == NORM_OPT_FRACAREA )
    {
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(_OPENMP)
#pragma omp parallel for default(none) \
  shared(num_wts, num_links, rv, tgt_grid)	\
  private(n, tgt_cell_add, norm_factor)
#endif
      for ( n = 0; n < num_links; ++n )
	{
	  tgt_cell_add = rv->tgt_cell_add[n];

          if ( IS_NOT_EQUAL(tgt_grid->cell_frac[tgt_cell_add], 0) )
	    norm_factor = ONE/tgt_grid->cell_frac[tgt_cell_add];
          else
            norm_factor = ZERO;

	  rv->wts[n*num_wts] *= norm_factor;
	}
    }
  else if ( rv->norm_opt == NORM_OPT_NONE )
    {
    }
}


void remap_weights_conserv(remapgrid_t *src_grid, remapgrid_t *tgt_grid, remapvars_t *rv)
{
  /* local variables */

  int    lcheck = TRUE;

  long   ioffset;
  long   src_num_cell_corners;
  long   tgt_num_cell_corners;
  long   src_cell_add;       /* current linear address for source grid cell   */
  long   tgt_cell_add;       /* current linear address for target grid cell   */
  long   k;                  /* generic counters                        */
  long   nbins;
  long   num_wts;
  long   max_srch_cells;     /* num cells in restricted search arrays  */
  long   num_srch_cells;     /* num cells in restricted search arrays  */
  long   srch_corners;       /* num of corners of srch cells           */
  int*   srch_add;           /* global address of cells in srch arrays */
  int    i;

  /* Variables necessary if segment manages to hit pole */
  long nx = 0, ny = 0;
  int src_remap_grid_type = src_grid->remap_grid_type;
  int tgt_remap_grid_type = tgt_grid->remap_grid_type;
  double src_grid_bound_box[4];
  int lyac = FALSE;
  extern int timer_remap_con;

  if ( cdoVerbose ) cdoPrint("Called %s()", __func__);

  progressInit();

  nbins = src_grid->num_srch_bins;
  num_wts = rv->num_wts;

  if ( cdoTimer ) timer_start(timer_remap_con);

  long src_grid_size = src_grid->size;
  long tgt_grid_size = tgt_grid->size;

  src_num_cell_corners = src_grid->num_cell_corners;
  tgt_num_cell_corners = tgt_grid->num_cell_corners;

  int max_num_cell_corners = src_num_cell_corners;
  if ( tgt_num_cell_corners > max_num_cell_corners ) max_num_cell_corners = tgt_num_cell_corners;

  enum yac_edge_type great_circle_type[32];
  for ( int i = 0; i < max_num_cell_corners; ++i ) great_circle_type[i] = GREAT_CIRCLE;

  enum yac_edge_type lonlat_circle_type[] = {LON_CIRCLE, LAT_CIRCLE, LON_CIRCLE, LAT_CIRCLE, LON_CIRCLE};

  enum yac_edge_type *src_edge_type = great_circle_type;
  enum yac_edge_type *tgt_edge_type = great_circle_type;

  enum cell_type target_cell_type = UNDEF_CELL;

  if ( src_num_cell_corners == 4 )
    {
      int lonlat_circle_index = get_lonlat_circle_index(src_grid);
      if ( lonlat_circle_index >= 0 ) src_edge_type = &lonlat_circle_type[lonlat_circle_index];
    }

  if ( tgt_num_cell_corners == 4 )
    {
      int lonlat_circle_index = get_lonlat_circle_index(tgt_grid);
      if ( lonlat_circle_index >= 0 )
	{
	  target_cell_type = LON_LAT_CELL;
	  tgt_edge_type = &lonlat_circle_type[lonlat_circle_index];
	}
    }

  if ( !(tgt_num_cell_corners < 4 || target_cell_type == LON_LAT_CELL) )
    {
      if ( tgt_grid->cell_center_lon == NULL || tgt_grid->cell_center_lat == NULL )
	cdoAbort("Internal problem (remap_weights_conserv): missing target point coordinates!");
    }

  double tgt_area;

  struct grid_cell* tgt_grid_cell;
  struct grid_cell* tgt_grid_cell2[ompNumThreads];  
  for ( i = 0; i < ompNumThreads; ++i )
    {
      tgt_grid_cell2[i] = (struct grid_cell*) malloc(sizeof(struct grid_cell));
      tgt_grid_cell2[i]->array_size      = tgt_num_cell_corners;
      tgt_grid_cell2[i]->num_corners     = tgt_num_cell_corners;
      tgt_grid_cell2[i]->edge_type       = tgt_edge_type;
      tgt_grid_cell2[i]->coordinates_x   = (double*) malloc(tgt_num_cell_corners*sizeof(double));
      tgt_grid_cell2[i]->coordinates_y   = (double*) malloc(tgt_num_cell_corners*sizeof(double));
      tgt_grid_cell2[i]->coordinates_xyz = (double*) malloc(3*tgt_num_cell_corners*sizeof(double));
    }

  struct grid_cell* src_grid_cells;
  struct grid_cell* overlap_buffer;
  struct grid_cell* src_grid_cells2[ompNumThreads];
  struct grid_cell* overlap_buffer2[ompNumThreads];
  for ( i = 0; i < ompNumThreads; ++i )
    {
      src_grid_cells2[i] = NULL;
      overlap_buffer2[i] = NULL;
    }

  double* partial_areas;
  double* partial_weights;
  double* partial_areas2[ompNumThreads];
  double* partial_weights2[ompNumThreads];
  for ( i = 0; i < ompNumThreads; ++i )
    {
      partial_areas2[i]   = NULL;
      partial_weights2[i] = NULL;
    }

  long max_srch_cells2[ompNumThreads];
  for ( i = 0; i < ompNumThreads; ++i )
    max_srch_cells2[i] = 0;

  int* srch_add2[ompNumThreads];
  for ( i = 0; i < ompNumThreads; ++i )
    srch_add2[i] = (int*) malloc(src_grid_size*sizeof(int));

  srch_corners = src_num_cell_corners;

  if ( src_remap_grid_type == REMAP_GRID_TYPE_REG2D )
    {
      nx = src_grid->dims[0];
      ny = src_grid->dims[1];
     
      src_grid_bound_box[0] = src_grid->reg2d_corner_lat[0];
      src_grid_bound_box[1] = src_grid->reg2d_corner_lat[ny];
      if ( src_grid_bound_box[0] > src_grid_bound_box[1] )
	{
	  src_grid_bound_box[0] = src_grid->reg2d_corner_lat[ny];
	  src_grid_bound_box[1] = src_grid->reg2d_corner_lat[0];
	}
      src_grid_bound_box[2] = src_grid->reg2d_corner_lon[0];
      src_grid_bound_box[3] = src_grid->reg2d_corner_lon[nx];
      //printf("src_grid   lon: %g %g lat: %g %g\n", RAD2DEG*src_grid_bound_box[2],RAD2DEG*src_grid_bound_box[3],RAD2DEG*src_grid_bound_box[0],RAD2DEG*src_grid_bound_box[1] );
    }

  weightlinks_t *weightlinks = (weightlinks_t *) malloc(tgt_grid_size*sizeof(weightlinks_t));
  
  double findex = 0;

  int sum_srch_cells = 0;
  int sum_srch_cells2 = 0;

  /* Loop over destination grid */

#if defined(_OPENMP)
#pragma omp parallel for default(shared) \
  shared(ompNumThreads, lyac, nbins, num_wts, src_remap_grid_type, tgt_remap_grid_type, src_grid_bound_box,	\
	 src_edge_type, tgt_edge_type, partial_areas2, partial_weights2,  \
         rv, cdoVerbose, max_srch_cells2, tgt_num_cell_corners, target_cell_type, \
         weightlinks, \
         srch_corners, src_grid, tgt_grid, tgt_grid_size, src_grid_size,	\
	 overlap_buffer2, src_grid_cells2, srch_add2, tgt_grid_cell2, findex, sum_srch_cells, sum_srch_cells2) \
  private(srch_add, tgt_grid_cell, tgt_area, k, num_srch_cells, max_srch_cells,  \
	  partial_areas, partial_weights, overlap_buffer, src_grid_cells, src_cell_add, tgt_cell_add, ioffset)
#endif
  for ( tgt_cell_add = 0; tgt_cell_add < tgt_grid_size; ++tgt_cell_add )
    {
      double partial_weight;
      long n, num_weights, num_weights_old;
      int ompthID = cdo_omp_get_thread_num();
      int lprogress = 1;
      if ( ompthID != 0 ) lprogress = 0;

#if defined(_OPENMP)
#include "pragma_omp_atomic_update.h"
#endif
      findex++;
      if ( lprogress ) progressStatus(0, 1, findex/tgt_grid_size);

      weightlinks[tgt_cell_add].nlinks = 0;	

      srch_add = srch_add2[ompthID];
      tgt_grid_cell = tgt_grid_cell2[ompthID];

      /* Get search cells */

      if ( src_remap_grid_type == REMAP_GRID_TYPE_REG2D && tgt_remap_grid_type == REMAP_GRID_TYPE_REG2D )
	{
	  double tgt_cell_bound_box[4];
	  boundbox_from_corners_reg2d(tgt_cell_add, tgt_grid->dims, tgt_grid->reg2d_corner_lon, tgt_grid->reg2d_corner_lat, tgt_cell_bound_box);
	  restrict_boundbox(src_grid_bound_box, tgt_cell_bound_box);
	  if ( 0 && cdoVerbose )
	    printf("bound_box %ld  lon: %g %g lat: %g %g\n",
		   tgt_cell_add, RAD2DEG*tgt_cell_bound_box[2],RAD2DEG*tgt_cell_bound_box[3],RAD2DEG*tgt_cell_bound_box[0],RAD2DEG*tgt_cell_bound_box[1] );
	  num_srch_cells = get_srch_cells_reg2d(src_grid->dims, src_grid->reg2d_corner_lat, src_grid->reg2d_corner_lon,
						tgt_cell_bound_box, srch_add);

	  if ( num_srch_cells == 1 && src_grid->dims[0] == 1 && src_grid->dims[1] == 1 &&
	       IS_EQUAL(src_grid->reg2d_corner_lat[0], src_grid->reg2d_corner_lat[1]) && 
	       IS_EQUAL(src_grid->reg2d_corner_lon[0], src_grid->reg2d_corner_lon[1]) ) num_srch_cells = 0;
	}
      else if ( src_remap_grid_type == REMAP_GRID_TYPE_REG2D )
	{
	  double tgt_cell_bound_box[4];
	  boundbox_from_corners1(tgt_cell_add, tgt_num_cell_corners, tgt_grid->cell_corner_lon, tgt_grid->cell_corner_lat, tgt_cell_bound_box);
	  restrict_boundbox(src_grid_bound_box, tgt_cell_bound_box);
	  if ( 0 && cdoVerbose )
	    printf("bound_box %ld  lon: %g %g lat: %g %g\n",
		   tgt_cell_add, RAD2DEG*tgt_cell_bound_box[2],RAD2DEG*tgt_cell_bound_box[3],RAD2DEG*tgt_cell_bound_box[0],RAD2DEG*tgt_cell_bound_box[1] );
	  num_srch_cells = get_srch_cells_reg2d(src_grid->dims, src_grid->reg2d_corner_lat, src_grid->reg2d_corner_lon,
						tgt_cell_bound_box, srch_add);

	  if ( num_srch_cells == 1 && src_grid->dims[0] == 1 && src_grid->dims[1] == 1 &&
	       IS_EQUAL(src_grid->reg2d_corner_lat[0], src_grid->reg2d_corner_lat[1]) && 
	       IS_EQUAL(src_grid->reg2d_corner_lon[0], src_grid->reg2d_corner_lon[1]) ) num_srch_cells = 0;
	}
      else
	{
	  restr_t tgt_cell_bound_box_r[4];
	  boundbox_from_corners1r(tgt_cell_add, tgt_num_cell_corners, tgt_grid->cell_corner_lon, tgt_grid->cell_corner_lat, tgt_cell_bound_box_r);

	  num_srch_cells = get_srch_cells(tgt_cell_add, nbins, tgt_grid->bin_addr, src_grid->bin_addr,
					  tgt_cell_bound_box_r, src_grid->cell_bound_box, src_grid_size, srch_add);
	}

      if ( 0 && cdoVerbose ) sum_srch_cells += num_srch_cells;

      if ( 0 && cdoVerbose )
	printf("tgt_cell_add %ld  num_srch_cells %ld\n", tgt_cell_add, num_srch_cells);

      if ( num_srch_cells == 0 ) continue;

      if ( tgt_remap_grid_type == REMAP_GRID_TYPE_REG2D )
	{
	  long nx = tgt_grid->dims[0];
	  long iy = tgt_cell_add/nx;
	  long ix = tgt_cell_add - iy*nx;

	  tgt_grid_cell->coordinates_x[0] = tgt_grid->reg2d_corner_lon[ix  ];
	  tgt_grid_cell->coordinates_y[0] = tgt_grid->reg2d_corner_lat[iy  ];
	  tgt_grid_cell->coordinates_x[1] = tgt_grid->reg2d_corner_lon[ix+1];
	  tgt_grid_cell->coordinates_y[1] = tgt_grid->reg2d_corner_lat[iy  ];
	  tgt_grid_cell->coordinates_x[2] = tgt_grid->reg2d_corner_lon[ix+1];
	  tgt_grid_cell->coordinates_y[2] = tgt_grid->reg2d_corner_lat[iy+1];
	  tgt_grid_cell->coordinates_x[3] = tgt_grid->reg2d_corner_lon[ix  ];
	  tgt_grid_cell->coordinates_y[3] = tgt_grid->reg2d_corner_lat[iy+1];
	}
      else
	{
	  for ( int ic = 0; ic < tgt_num_cell_corners; ++ic )
	    {
	      tgt_grid_cell->coordinates_x[ic] = tgt_grid->cell_corner_lon[tgt_cell_add*tgt_num_cell_corners+ic];
	      tgt_grid_cell->coordinates_y[ic] = tgt_grid->cell_corner_lat[tgt_cell_add*tgt_num_cell_corners+ic];
	    }
	}
      
      for ( int ic = 0; ic < tgt_num_cell_corners; ++ic )
	LLtoXYZ(tgt_grid_cell->coordinates_x[ic], tgt_grid_cell->coordinates_y[ic], tgt_grid_cell->coordinates_xyz+ic*3);

      //printf("target: %ld\n", tgt_cell_add);
      if ( lyac )
        if ( tgt_cell_add == 174752 )
	  {
	    for ( int n = 0; n < tgt_num_cell_corners; ++n )
	      {
		printf("  TargetCell.coordinates_x[%d] = %g*rad;\n", n, tgt_grid_cell->coordinates_x[n]/DEG2RAD);
		printf("  TargetCell.coordinates_y[%d] = %g*rad;\n", n, tgt_grid_cell->coordinates_y[n]/DEG2RAD);
	      }
	    /*
	    printf("> -Z1\n");
	    for ( int n = 0; n < tgt_num_cell_corners; ++n )
		printf("  %g %g\n", tgt_grid_cell->coordinates_x[n]/DEG2RAD, tgt_grid_cell->coordinates_y[n]/DEG2RAD);
	      printf("  %g %g\n", tgt_grid_cell->coordinates_x[0]/DEG2RAD, tgt_grid_cell->coordinates_y[0]/DEG2RAD);
	    */
	  }
      
      /* Create search arrays */

      max_srch_cells  = max_srch_cells2[ompthID];
      partial_areas   = partial_areas2[ompthID];
      partial_weights = partial_weights2[ompthID];
      overlap_buffer  = overlap_buffer2[ompthID];
      src_grid_cells  = src_grid_cells2[ompthID];

      if ( num_srch_cells > max_srch_cells )
	{
	  partial_areas   = (double*) realloc(partial_areas,   num_srch_cells*sizeof(double));
	  partial_weights = (double*) realloc(partial_weights, num_srch_cells*sizeof(double));

	  overlap_buffer = (struct grid_cell*) realloc(overlap_buffer, num_srch_cells*sizeof(struct grid_cell));
	  src_grid_cells = (struct grid_cell*) realloc(src_grid_cells, num_srch_cells*sizeof(struct grid_cell));

	  for ( n = max_srch_cells; n < num_srch_cells; ++n )
	    {
	      overlap_buffer[n].array_size      = 0;
	      overlap_buffer[n].num_corners     = 0;
	      overlap_buffer[n].edge_type       = NULL;
	      overlap_buffer[n].coordinates_x   = NULL;
	      overlap_buffer[n].coordinates_y   = NULL;
	      overlap_buffer[n].coordinates_xyz = NULL;
	    }

	  for ( n = max_srch_cells; n < num_srch_cells; ++n )
	    {
	      src_grid_cells[n].array_size      = srch_corners;
	      src_grid_cells[n].num_corners     = srch_corners;
	      src_grid_cells[n].edge_type       = src_edge_type;
	      src_grid_cells[n].coordinates_x   = (double*) malloc(srch_corners*sizeof(double));
	      src_grid_cells[n].coordinates_y   = (double*) malloc(srch_corners*sizeof(double));
	      src_grid_cells[n].coordinates_xyz = (double*) malloc(3*srch_corners*sizeof(double));
	    }

	  max_srch_cells = num_srch_cells;

	  max_srch_cells2[ompthID]  = max_srch_cells;
	  partial_areas2[ompthID]   = partial_areas;
	  partial_weights2[ompthID] = partial_weights;
	  overlap_buffer2[ompthID]  = overlap_buffer;
	  src_grid_cells2[ompthID]  = src_grid_cells;
	}

      // printf("  int ii = 0;\n");
      for ( n = 0; n < num_srch_cells; ++n )
	{
	  long srch_corners_new = srch_corners;

	  src_cell_add = srch_add[n];

	  if ( src_remap_grid_type == REMAP_GRID_TYPE_REG2D )
	    {
	      int ix, iy;

	      iy = src_cell_add/nx;
	      ix = src_cell_add - iy*nx;

	      src_grid_cells[n].coordinates_x[0] = src_grid->reg2d_corner_lon[ix  ];
	      src_grid_cells[n].coordinates_y[0] = src_grid->reg2d_corner_lat[iy  ];
	      src_grid_cells[n].coordinates_x[1] = src_grid->reg2d_corner_lon[ix+1];
	      src_grid_cells[n].coordinates_y[1] = src_grid->reg2d_corner_lat[iy  ];
	      src_grid_cells[n].coordinates_x[2] = src_grid->reg2d_corner_lon[ix+1];
	      src_grid_cells[n].coordinates_y[2] = src_grid->reg2d_corner_lat[iy+1];
	      src_grid_cells[n].coordinates_x[3] = src_grid->reg2d_corner_lon[ix  ];
	      src_grid_cells[n].coordinates_y[3] = src_grid->reg2d_corner_lat[iy+1];
	      /*
	      printf("source1: %ld %ld", num_srch_cells, n);
	      for ( k = 0; k < srch_corners; ++k )
		printf(" %g %g", src_grid_cells[n].coordinates_x[k]/DEG2RAD, src_grid_cells[n].coordinates_y[k]/DEG2RAD);
	      printf("\n");
	      */
	    }
	  else
	    {
	      ioffset = src_cell_add*srch_corners;
	      /*
	      for ( k = srch_corners-1; k > 0; --k )
		{
		  if ( IS_NOT_EQUAL(src_grid->cell_corner_lon[ioffset+k], src_grid->cell_corner_lon[ioffset+k-1]) ||
		       IS_NOT_EQUAL(src_grid->cell_corner_lat[ioffset+k], src_grid->cell_corner_lat[ioffset+k-1]) )
		    break;
		}
	      if ( k != srch_corners-1 ) printf("%ld %ld %ld %ld\n", tgt_cell_add, n, srch_corners, k+1);

	      if ( k != srch_corners-1 )
		{
		  srch_corners_new = k+1;
		  src_grid_cells[n].num_corners = srch_corners_new;
		}
	      */
	      for ( k = 0; k < srch_corners_new; ++k )
		{
		  src_grid_cells[n].coordinates_x[k] = src_grid->cell_corner_lon[ioffset+k];
		  src_grid_cells[n].coordinates_y[k] = src_grid->cell_corner_lat[ioffset+k];
		}
	      /*
	      for ( k = 0; k < srch_corners_new; ++k )
		{
		  printf("  SourceCell[ii].coordinates_x[%ld] = %g*rad;\n", k, src_grid_cells[n].coordinates_x[k]/DEG2RAD);
		  printf("  SourceCell[ii].coordinates_y[%ld] = %g*rad;\n", k, src_grid_cells[n].coordinates_y[k]/DEG2RAD);
		}
	      */
	      /*
	      printf("source2: %ld %ld", num_srch_cells, n);
	      for ( k = 0; k < srch_corners_new; ++k )
		printf(" %g %g", src_grid_cells[n].coordinates_x[k]/DEG2RAD, src_grid_cells[n].coordinates_y[k]/DEG2RAD);
	      printf("\n");
	      */
	    }

	  for ( int ic = 0; ic < srch_corners_new; ++ic )
	    LLtoXYZ(src_grid_cells[n].coordinates_x[ic], src_grid_cells[n].coordinates_y[ic], src_grid_cells[n].coordinates_xyz+ic*3);

	  if ( lyac )
	    if ( tgt_cell_add == 174752 )
	    {
	      // printf("n %d\n", (int)n);
	      for ( k = 0; k < srch_corners_new; ++k )
		{
		  printf("  SourceCell[ii].coordinates_x[%ld] = %g*rad;\n", k, src_grid_cells[n].coordinates_x[k]/DEG2RAD);
		  printf("  SourceCell[ii].coordinates_y[%ld] = %g*rad;\n", k, src_grid_cells[n].coordinates_y[k]/DEG2RAD);
		}
	      printf("  ii++;\n");
	      /*
	      printf("> -Z1\n");
	      for ( k = 0; k < srch_corners_new; ++k )
		printf("  %g %g\n", src_grid_cells[n].coordinates_x[k]/DEG2RAD, src_grid_cells[n].coordinates_y[k]/DEG2RAD);
	      printf("  %g %g\n", src_grid_cells[n].coordinates_x[0]/DEG2RAD, src_grid_cells[n].coordinates_y[0]/DEG2RAD);
	      */
	    }
	}

      if ( tgt_num_cell_corners < 4 || target_cell_type == LON_LAT_CELL )
	{
	  cdo_compute_overlap_areas(num_srch_cells, overlap_buffer, src_grid_cells, *tgt_grid_cell, partial_areas);
	}
      else
	{
	  double cell_center_lon = tgt_grid->cell_center_lon[tgt_cell_add];
	  double cell_center_lat = tgt_grid->cell_center_lat[tgt_cell_add];
	  cdo_compute_concave_overlap_areas(num_srch_cells, overlap_buffer, src_grid_cells, *tgt_grid_cell, cell_center_lon, cell_center_lat, partial_areas);
	}

      tgt_area = gridcell_area(*tgt_grid_cell);
      // tgt_area = cell_area(tgt_grid_cell);

      for ( num_weights = 0, n = 0; n < num_srch_cells; ++n )
	{
	  if ( partial_areas[n] > 0 )
	    {
	      //printf(">>>>   %d %d %g %g\n", (int)tgt_cell_add, srch_add[n], tgt_area, partial_areas[n]);
	      partial_areas[num_weights] = partial_areas[n];
	      srch_add[num_weights] = srch_add[n];
	      num_weights++;
	    }
	}

      if ( 0 && cdoVerbose ) sum_srch_cells2 += num_weights;

      for ( n = 0; n < num_weights; ++n )
	partial_weights[n] = partial_areas[n] / tgt_area;

      if ( rv->norm_opt == NORM_OPT_FRACAREA )
	yac_correct_weights((unsigned)num_weights, partial_weights);

      for ( n = 0; n < num_weights; ++n )
	partial_weights[n] *= tgt_area;
      //#endif

      num_weights_old = num_weights;
      for ( num_weights = 0, n = 0; n < num_weights_old; ++n )
	{
	  src_cell_add = srch_add[n];

	  if ( 0 && cdoVerbose )
	    printf("tgt_cell_add %ld, src_cell_add %ld,  partial_weights[n] %g, tgt_area  %g\n", tgt_cell_add, src_cell_add, partial_weights[n], tgt_area);

	  if ( partial_weights[n] <= 0. ) src_cell_add = -1;
	  if ( src_cell_add != -1 )
	    {
	      partial_weights[num_weights] = partial_weights[n];
	      srch_add[num_weights] = src_cell_add;
	      num_weights++;
	    }
	}

      for ( n = 0; n < num_weights; ++n )
	{
	  partial_weight = partial_weights[n];

	  src_cell_add = srch_add[n];

#if defined(_OPENMP)
#pragma omp atomic
#endif
	  src_grid->cell_area[src_cell_add] += partial_weight;
	}


      num_weights_old = num_weights;
      for ( num_weights = 0, n = 0; n < num_weights_old; ++n )
	{
	  src_cell_add = srch_add[n];

	  /*
	    Store the appropriate addresses and weights. 
	    Also add contributions to cell areas.
	    The source grid mask is the master mask
	  */
	  if ( src_grid->mask[src_cell_add] )
	    {
	      partial_weights[num_weights] = partial_weights[n];
	      srch_add[num_weights] = src_cell_add;
	      num_weights++;
	    }
	}

      for ( n = 0; n < num_weights; ++n )
	{
	  partial_weight = partial_weights[n];
	  src_cell_add = srch_add[n];

#if defined(_OPENMP)
#pragma omp atomic
#endif
	  src_grid->cell_frac[src_cell_add] += partial_weight;
		  
	  tgt_grid->cell_frac[tgt_cell_add] += partial_weight;
	}

      store_weightlinks(num_weights, srch_add, partial_weights, tgt_cell_add, weightlinks);

      tgt_grid->cell_area[tgt_cell_add] = tgt_area; 
      // printf("area %d %g %g\n", tgt_cell_add, tgt_grid->cell_area[tgt_cell_add], tgt_area);
    }

  if ( 0 && cdoVerbose )
    {
      printf("sum_srch_cells : %d\n", sum_srch_cells);
      printf("sum_srch_cells2: %d\n", sum_srch_cells2);
    }

  /* Finished with all cells: deallocate search arrays */
  long n;

  for ( i = 0; i < ompNumThreads; ++i )
    {
      for ( n = 0; n < max_srch_cells2[i]; n++ )
	{
	  if ( overlap_buffer2[i][n].array_size > 0 )
	    {
	      free(overlap_buffer2[i][n].coordinates_x);
	      free(overlap_buffer2[i][n].coordinates_y);
	      if ( overlap_buffer2[i][n].coordinates_xyz ) free(overlap_buffer2[i][n].coordinates_xyz);
	      if ( overlap_buffer2[i][n].edge_type ) free(overlap_buffer2[i][n].edge_type);
	    }
	}
      for ( n = 0; n < max_srch_cells2[i]; n++ )
	{
	  free(src_grid_cells2[i][n].coordinates_x);
	  free(src_grid_cells2[i][n].coordinates_y);
	  free(src_grid_cells2[i][n].coordinates_xyz);
	}
      free(src_grid_cells2[i]);
      free(overlap_buffer2[i]);

      free(partial_areas2[i]);
      free(partial_weights2[i]);

      free(tgt_grid_cell2[i]->coordinates_x);
      free(tgt_grid_cell2[i]->coordinates_y);
      free(tgt_grid_cell2[i]->coordinates_xyz);
      free(tgt_grid_cell2[i]);

      free(srch_add2[i]);
    }

  weightlinks2remaplinks(tgt_grid_size, weightlinks, rv);

  if ( weightlinks ) free(weightlinks);

  /* Normalize using destination area if requested */
  normalize_weights(tgt_grid, rv);

  long num_links = rv->num_links;

  if ( cdoVerbose )
    cdoPrint("Total number of links = %ld", rv->num_links);
  
  for ( n = 0; n < src_grid_size; ++n )
    if ( IS_NOT_EQUAL(src_grid->cell_area[n], 0) ) src_grid->cell_frac[n] /= src_grid->cell_area[n];

  for ( n = 0; n < tgt_grid_size; ++n )
    if ( IS_NOT_EQUAL(tgt_grid->cell_area[n], 0) ) tgt_grid->cell_frac[n] /= tgt_grid->cell_area[n];

  /* Perform some error checking on final weights  */

  if ( lcheck )
    {
      for ( n = 0; n < src_grid_size; ++n )
	{
	  if ( src_grid->cell_area[n] < -.01 )
	    cdoPrint("Source grid area error: %d %g", n, src_grid->cell_area[n]);
	}

      for ( n = 0; n < tgt_grid_size; ++n )
	{
	  if ( tgt_grid->cell_area[n] < -.01 )
	    cdoPrint("Target grid area error: %d %g", n, tgt_grid->cell_area[n]);
	}

      for ( n = 0; n < num_links; ++n )
	{
	  src_cell_add = rv->src_cell_add[n];
	  tgt_cell_add = rv->tgt_cell_add[n];

	  if ( rv->wts[n*num_wts] < -0.01 )
	    cdoPrint("Map weight < 0! grid1idx=%d grid2idx=%d nlink=%d wts=%g",
		     src_cell_add, tgt_cell_add, n, rv->wts[n*num_wts]);

	  if ( rv->norm_opt != NORM_OPT_NONE && rv->wts[n*num_wts] > 1.01 )
	    cdoPrint("Map weight > 1! grid1idx=%d grid2idx=%d nlink=%d wts=%g",
		     src_cell_add, tgt_cell_add, n, rv->wts[n*num_wts]);
	}
    } // lcheck

  if ( cdoTimer ) timer_stop(timer_remap_con);

} /* remap_weights_conserv */


void remap_conserv(remapgrid_t *src_grid, remapgrid_t *tgt_grid, const double* restrict src_array, double* restrict tgt_array, double missval)
{
} /* remap_conserv */
