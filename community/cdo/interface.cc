

#include <stdio.h>
#include <string.h>
using namespace std;


#define restrict __restrict__

extern "C" {

#include "src/libcdi/src/cdi.h"
#include "src/src/cdo.h"
#include "src/src/cdo_int.h"
#include "src/src/grid.h"
#include "src/src/remap.h"

}



#define MAX_LENGTH 256



char remap_file_name[MAX_LENGTH];

remapgrid_t src_grid;
remapgrid_t dst_grid;
remapvars_t rv;
void print_remap_grid_info(remapgrid_t *grid);
void print_remap_info(remapvars_t *rv);

//for printing the map type as string
const char *map_type_string[] = { "unknown", "conservative", "bilinear", "bicubic", "distwgt", "conserv_yac" };


int get_remap_file(char** filename) {
    *filename = remap_file_name;
    return 0;

}


int set_remap_file(char* filename) {
//    strncpy(remap_file_name, filename, MAX_LENGTH);
//    remap_file_name = filename;

    //read_remap_scrip(const char *interp_file, int gridID1, int gridID2, int *map_type, int *submap_type, int *num_neighbors,
    //          int *remap_order, remapgrid_t *src_grid, remapgrid_t *tgt_grid, remapvars_t *rv);

    int map_type;
    int submap_type;
    int num_neighbors;
    int remap_order;

    cdoVerbose = 1;

    printf("going to call CDO\n");

    read_remap_scrip(filename, -1, -1, &map_type, &submap_type, &num_neighbors,
              &remap_order, &src_grid, &dst_grid, &rv);

    printf("finished calling CDO\n");

    printf("source grid:\n");
    print_remap_grid_info(&src_grid);

    printf("destination grid:\n");
    print_remap_grid_info(&dst_grid);

    printf("remap info:\n");
    print_remap_info(&rv);


    return 0;
}



int get_src_grid_size(int * size) {
  *size = src_grid.size;
  return 0;
}

int get_dst_grid_size(int * size) {
  *size = dst_grid.size;
  return 0;
}

int get_src_grid_dims(int *x, int *y) {
  *x = src_grid.dims[0];
  *y = src_grid.dims[1];
  return 0;
}
int get_dst_grid_dims(int *x, int *y) {
  *x = dst_grid.dims[0];
  *y = dst_grid.dims[1];
  return 0;
}



/*
 * It would appear that CDO leaves the grid masks empty in the weights
 * file, even if they are supplied when the weights were generated.
 */
int get_src_grid_mask(int *index_i, int *index_j, int *imask, int n) {
  int i;
  int nx = src_grid.dims[0];
  int ny = src_grid.dims[1];

  printf("get_src_grid_mask() called\n");
  printf("nx=%d, ny=%d, n=%d",nx, ny, n);

  for (i=0; i<n; i++) {
      imask[index_i[i]*nx+index_j[i]] = src_grid.vgpm[index_i[i]*nx+index_j[i]];
  }

  return 0;
}





//gradients used for second order and bicubic remappings
//to be computed when the weights are read in or computed
double *src_grad1 = NULL;
double *src_grad2 = NULL;
double *src_grad3 = NULL;


int perform_remap(int *index_i, double *dst, double *src, int n) {

  /*
    This function calls the remapping function in CDO, it also performs some preprocessing and checks that
    would otherwise be performed by the Remap() function that implements the operator

    Arguments: 
        index_i is a list of indexes in the destination grid for which a remapped value from the source grid is requested
        dst is the array that contains the result of this remapping
        src is the array that contains the source grid to be remapped
        n is the number of grid points to be remapped, which should be equal to the number of grid points in the destination grid

    For now partial remappings are not supported n must equal dst grid size.


    void remap(double *restrict dst_array, double missval, long dst_size, long num_links, double *restrict map_wts,
       long num_wts, const int *restrict dst_add, const int *restrict src_add, const double *restrict src_array,
       const double *restrict src_grad1, const double *restrict src_grad2, const double *restrict src_grad3,
       remaplink_t links)
    Input:
    int *dst_add         ! destination address for each link
    int *src_add         ! source      address for each link
    int num_wts          ! num of weights used in remapping
    double *map_wts      ! remapping weights for each link
    double *src_array    ! array with source field to be remapped
    optional inputs:
    double *src_grad1    ! gradient arrays on source grid necessary for
    double *src_grad2    ! higher-order remappings
    double *src_grad3
    output variables:
    Output:
    double *dst_array    ! array for remapped field on destination grid
  */

  double missval = 0.0; //have to find out what this is exactly

  remap(dst, missval, dst_grid.size, rv.num_links, rv.wts,
       rv.num_wts, rv.tgt_cell_add, rv.src_cell_add, src,
       src_grad1, src_grad2, src_grad3,
       rv.links);

  return 0;
}









void print_remap_grid_info(remapgrid_t *grid) {

  printf("lwrite_remap %d\n", grid->lwrite_remap);
  printf("gridID %d\n", grid->gridID);
  printf("store_link_fast %d\n", grid->store_link_fast);
  printf("remap_grid_type %d\n", grid->remap_grid_type);
  printf("lextrapolate %d\n", grid->lextrapolate);
  printf("non_global %d\n", grid->non_global);
  printf("is_cyclic %d\n", grid->is_cyclic);
  printf("rank %d\n", grid->rank);                  /* rank of the grid */
  printf("size %ld\n", grid->size);                  /* total points on the grid */
  printf("num_cell_corners %ld\n", grid->num_cell_corners);      /* number of corners for each grid cell */

  int nx = grid->dims[0];
  int ny = grid->dims[1];

  printf("dims %d, %d\n", nx, ny);               /* size of grid dimension */

  printf("nvgp %d\n", grid->nvgp);                  /* size of vgpm           */

//  printf("vgpm %d\n",*     vgpm;                  /* flag which cells are valid   */

//  printf("mask %d\n",*     mask;                  /* flag which cells participate */

/*
  double*  reg2d_center_lon;      // reg2d lon/lat coordinates for 
  double*  reg2d_center_lat;      // each grid center in radians   

  double*  reg2d_corner_lon;      // reg2d lon/lat coordinates for 
  double*  reg2d_corner_lat;      // each grid corner in radians   

  double*  cell_center_lon;       // lon/lat coordinates for       
  double*  cell_center_lat;       // each grid center in radians   

  double*  cell_corner_lon;       // lon/lat coordinates for         
  double*  cell_corner_lat;       // each grid corner in radians     

  double*  cell_area;             // tot area of each grid cell     
  double*  cell_frac;             // fractional area of grid cells participating in remapping  
*/

  printf("lneed_cell_corners %d\n", grid->lneed_cell_corners);
  printf("luse_cell_corners %d\n", grid->luse_cell_corners);     /* use corners for bounding boxes  */

//  restr_t *cell_bound_box;        /* lon/lat bounding box for use    */

  printf("num_srch_bins %d\n", grid->num_srch_bins);         /* num of bins for restricted srch */

//  printf("bin_addr %d", *bin_addr);              /* min,max adds for grid cells in this lat bin  */

//  restr_t* bin_lats;              /* min,max latitude for each search bin   */

  printf("\n");

}


void print_remap_info(remapvars_t *rv) {

  printf("sort_add %d\n", rv->sort_add);
  printf("pinit %d\n", rv->pinit);            // TRUE if the pointers are initialized     
  printf("max_links %lu\n", rv->max_links);        // current size of link arrays              
  printf("num_links %lu\n", rv->num_links);        // actual number of links for remapping     
  printf("num_wts %lu\n", rv->num_wts);          // num of weights used in remapping         
  printf("map_type %d %s\n", rv->map_type, map_type_string[rv->map_type]);         // identifier for remapping method

  printf("norm_opt %d\n", rv->norm_opt);         // option for normalization (conserv only)  
  printf("resize_increment %d\n", rv->resize_increment); // default amount to increase array size    

//  int*     src_cell_add;     // source grid address for each link        
  printf("rv->src_cell_add %p\n", rv->src_cell_add);

//  int*     tgt_cell_add;     // target grid address for each link        
  printf("rv->tgt_cell_add %p\n", rv->tgt_cell_add);

//  double*  wts;              // map weights for each link [max_links*num_wts] 
  printf("rv->wts %p\n", rv->wts);

  //remaplinks_t links
  printf("links.option %d\n",      rv->links.option);
  printf("links.max_links %d\n",      rv->links.max_links);
  printf("links.num_blks %d\n",      rv->links.num_blks);


  printf("links.num_links %p\n",     rv->links.num_links);


//  int**    rv->links.src_add;
  printf("rv->links.src_add %p\n", rv->links.src_add);
  if (rv->links.src_add) {
    printf("rv->links.src_add %p\n", *(rv->links.src_add));
  }

//  int**    rv->links.dst_add;
  printf("rv->links.dst_add %p\n", rv->links.dst_add);
  if (rv->links.dst_add) {
    printf("rv->links.dst_add %p\n", *(rv->links.dst_add));
  }

//  int**    rv->links.w_index;
  printf("rv->links.w_index %p\n", rv->links.w_index);
  if (rv->links.w_index) {
    printf("rv->links.w_index %p\n", *(rv->links.w_index));
  }


  printf("\n");

}
