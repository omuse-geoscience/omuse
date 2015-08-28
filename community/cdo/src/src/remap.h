#ifndef _REMAP_H
#define _REMAP_H

#include <math.h>

#ifndef  M_PI
#define  M_PI        3.14159265358979323846264338327950288  /* pi */
#endif

#define  PI       M_PI
#define  PI2      (2.0*PI)
#define  PIH      (0.5*PI)

#define  ZERO     0.0
#define  ONE      1.0
#define  TWO      2.0
#define  THREE    3.0
#define  HALF     0.5
#define  QUART    0.25
#define  BIGNUM   1.e+20
#define  TINY     1.e-14


#define  REMAP_GRID_TYPE_REG2D     1
#define  REMAP_GRID_TYPE_CURVE2D   2
#define  REMAP_GRID_TYPE_UNSTRUCT  3

#define  REMAP_GRID_BASIS_SRC      1
#define  REMAP_GRID_BASIS_TGT      2

#define  RESTR_TYPE  int  /* restrict data types: 0 -> double, float; 1 -> int */

typedef RESTR_TYPE restr_t;
/*
#if RESTR_TYPE == int
#  define RESTR_SCALE(x) ((int) (0.5+100000000*(x)))
#  define RESTR_ABS(x)   abs(x)
#else
#  define RESTR_SCALE(x) (x)
#  define RESTR_ABS(x)   fabs(x)
#endif
*/
/* short
#  define RESTR_SFAC     4000
#  define RESTR_SCALE(x) ((short) (0.5+RESTR_SFAC*(x)))
#  define RESTR_ABS(x)   abs(x)
*/
/* int */
#  define RESTR_SFAC     100000000
#  define RESTR_SCALE(x) ((int) (0.5+RESTR_SFAC*(x)))
#  define RESTR_ABS(x)   abs(x)
/*
#  define RESTR_SFAC     1.
#  define RESTR_SCALE(x) (x)
#  define RESTR_ABS(x)   fabs(x)
*/

#define  TINY_FRAC     1.e-10

#define  NORM_OPT_NONE        1
#define  NORM_OPT_DESTAREA    2
#define  NORM_OPT_FRACAREA    3

#define  MAP_TYPE_CONSERV     1
#define  MAP_TYPE_BILINEAR    2
#define  MAP_TYPE_BICUBIC     3
#define  MAP_TYPE_DISTWGT     4
#define  MAP_TYPE_CONSERV_YAC 5

#define  SUBMAP_TYPE_NONE     0
#define  SUBMAP_TYPE_LAF      1
#define  SUBMAP_TYPE_SUM      2


typedef struct {
  int      lwrite_remap;
  int      gridID;
  int      store_link_fast;
  int      remap_grid_type;
  int      lextrapolate;
  int      non_global;
  int      is_cyclic;
  int      rank;                  /* rank of the grid */
  long     size;                  /* total points on the grid */
  long     num_cell_corners;      /* number of corners for each grid cell */

  int      dims[2];               /* size of grid dimension */

  int      nvgp;                  /* size of vgpm           */
  int*     vgpm;                  /* flag which cells are valid   */

  int*     mask;                  /* flag which cells participate */

  double*  reg2d_center_lon;      /* reg2d lon/lat coordinates for */
  double*  reg2d_center_lat;      /* each grid center in radians   */
  double*  reg2d_corner_lon;      /* reg2d lon/lat coordinates for */
  double*  reg2d_corner_lat;      /* each grid corner in radians   */

  double*  cell_center_lon;       /* lon/lat coordinates for       */
  double*  cell_center_lat;       /* each grid center in radians   */
  double*  cell_corner_lon;       /* lon/lat coordinates for         */
  double*  cell_corner_lat;       /* each grid corner in radians     */

  double*  cell_area;             /* tot area of each grid cell     */
  double*  cell_frac;             /* fractional area of grid cells participating in remapping  */

  int      lneed_cell_corners;
  int      luse_cell_corners;     /* use corners for bounding boxes  */

  restr_t *cell_bound_box;        /* lon/lat bounding box for use    */

  int      num_srch_bins;         /* num of bins for restricted srch */

  int*     bin_addr;              /* min,max adds for grid cells in this lat bin  */

  restr_t* bin_lats;              /* min,max latitude for each search bin   */
}
remapgrid_t;

typedef struct {
  int      option;
  int      max_links;
  int      num_blks;
  int*     num_links;
  int**    src_add;
  int**    dst_add;
  int**    w_index;
}
remaplink_t;

typedef struct {
  int      sort_add;
  int      pinit;            /* TRUE if the pointers are initialized     */
  long     max_links;        /* current size of link arrays              */
  long     num_links;        /* actual number of links for remapping     */
  long     num_wts;          /* num of weights used in remapping         */
  int      map_type;         /* identifier for remapping method          */
  int      norm_opt;         /* option for normalization (conserv only)  */
  int      resize_increment; /* default amount to increase array size    */

  int*     src_cell_add;     /* source grid address for each link        */
  int*     tgt_cell_add;     /* target grid address for each link        */

  double*  wts;              /* map weights for each link [max_links*num_wts] */

  remaplink_t links;
}
remapvars_t;

typedef struct {
  int      gridID;
  int      gridsize;
  int      nmiss;
  remapgrid_t src_grid;
  remapgrid_t tgt_grid;
  remapvars_t vars;
}
remap_t;

#define  REMAP_STORE_LINK_FAST  1
#define  REMAP_WRITE_REMAP      2
#define  REMAP_MAX_ITER         3
#define  REMAP_NUM_SRCH_BINS    4
#define  REMAP_GENWEIGHTS       5

void remap_set_threshhold(double threshhold);
void remap_set_int(int remapvar, int value);


void remap_grids_init(int map_type, int lextrapolate, int gridID1, remapgrid_t *src_grid, int gridID2, remapgrid_t *tgt_grid);
void remap_vars_init(int map_type, long src_grid_size, long tgt_grid_size, remapvars_t *rv);

void remapVarsFree(remapvars_t *rv);
void remapGridFree(remapgrid_t *grid);

void remap(double *restrict dst_array, double missval, long dst_size, long num_links, double *restrict map_wts, 
	   long num_wts, const int *restrict dst_add, const int *restrict src_add, const double *restrict src_array, 
	   const double *restrict src_grad1, const double *restrict src_grad2, const double *restrict src_grad3,
	   remaplink_t links);

void remap_laf(double *restrict dst_array, double missval, long dst_size, long num_links, double *restrict map_wts,
	       long num_wts, const int *restrict dst_add, const int *restrict src_add, const double *restrict src_array);

void remap_sum(double *restrict dst_array, double missval, long dst_size, long num_links, double *restrict map_wts,
	       long num_wts, const int *restrict dst_add, const int *restrict src_add, const double *restrict src_array);

void scrip_remap_weights_bilinear(remapgrid_t *src_grid, remapgrid_t *tgt_grid, remapvars_t *rv);
void scrip_remap_weights_bicubic(remapgrid_t *src_grid, remapgrid_t *tgt_grid, remapvars_t *rv);
void scrip_remap_weights_distwgt(int num_neighbors, remapgrid_t *src_grid, remapgrid_t *tgt_grid, remapvars_t *rv);
void scrip_remap_weights_conserv(remapgrid_t *src_grid, remapgrid_t *tgt_grid, remapvars_t *rv);
void remap_weights_conserv(remapgrid_t *src_grid, remapgrid_t *tgt_grid, remapvars_t *rv);

void scrip_remap_bilinear(remapgrid_t* src_grid, remapgrid_t* tgt_grid, const double* restrict src_array, double* restrict tgt_array, double missval);
void scrip_remap_bicubic(remapgrid_t *src_grid, remapgrid_t *tgt_grid, const double* restrict src_array, double* restrict tgt_array, double missval);
void remap_conserv(remapgrid_t *src_grid, remapgrid_t *tgt_grid, const double* restrict src_array, double* restrict tgt_array, double missval);


void resize_remap_vars(remapvars_t *rv, int increment);

void remap_stat(int remap_order, remapgrid_t src_grid, remapgrid_t tgt_grid, remapvars_t rv, const double *restrict array1, 
		const double *restrict array2, double missval);
void remap_gradients(remapgrid_t grid, const double *restrict array, double *restrict grad_lat,
		     double *restrict grad_lon, double *restrict grad_latlon);

void reorder_links(remapvars_t *rv);

void sort_add(long num_links, long num_wts, int *restrict add1, int *restrict add2, double *restrict weights);
void sort_iter(long num_links, long num_wts, int *restrict add1, int *restrict add2, double *restrict weights, int parent);

void write_remap_scrip(const char *interp_file, int map_type, int submap_type, int num_neighbors,
		       int remap_order, remapgrid_t src_grid, remapgrid_t tgt_grid, remapvars_t rv);
void read_remap_scrip(const char *interp_file, int gridID1, int gridID2, int *map_type, int *submap_type, int *num_neighbors,
		      int *remap_order, remapgrid_t *src_grid, remapgrid_t *tgt_grid, remapvars_t *rv);

void calc_bin_addr(long gridsize, long nbins, const restr_t* restrict bin_lats, const restr_t* restrict cell_bound_box, int* restrict bin_addr);
void calc_lat_bins(remapgrid_t* src_grid, remapgrid_t* tgt_grid, int map_type);
long get_srch_cells(long tgt_cell_add, long nbins, int *bin_addr1, int *bin_addr2,
		    restr_t *tgt_cell_bound_box, restr_t *src_cell_bound_box, long src_grid_size, int *srch_add);

int grid_search_reg2d_nn(long nx, long ny, int *restrict nbr_add, double *restrict nbr_dist, double plat, double plon,
			 const double *restrict src_center_lat, const double *restrict src_center_lon);

int grid_search_reg2d(remapgrid_t *src_grid, int *restrict src_add, double *restrict src_lats, 
		      double *restrict src_lons,  double plat, double plon, const int *restrict src_grid_dims,
		      const double *restrict src_center_lat, const double *restrict src_center_lon);

int grid_search(remapgrid_t *src_grid, int *restrict src_add, double *restrict src_lats, 
		double *restrict src_lons,  double plat, double plon, const int *restrict src_grid_dims,
		const double *restrict src_center_lat, const double *restrict src_center_lon,
		const restr_t *restrict src_grid_bound_box, const int *restrict src_bin_add);

int find_ij_weights(double plon, double plat, double* restrict src_lats, double* restrict src_lons, double *ig, double *jg);
int rect_grid_search(long *ii, long *jj, double x, double y, long nxm, long nym, const double *restrict xm, const double *restrict ym);


#endif  /* _REMAP_H */
