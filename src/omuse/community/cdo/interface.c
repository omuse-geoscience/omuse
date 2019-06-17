
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define restrict __restrict__

#include "src/libcdi/src/cdi.h"
#include "src/src/cdo.h"
#include "src/src/cdo_int.h"
#include "src/src/grid.h"
#include "src/src/remap.h"
#include "src/src/griddes.h"

#undef malloc
#undef free

#define MAX_LENGTH 256

#define max(a,b) ((a) > (b) ? (a) : (b))



char weights_file[MAX_LENGTH];
char src_grid_file[MAX_LENGTH];
char dst_grid_file[MAX_LENGTH];

int gridID1 = -1;
int gridID2 = -1;

//remap grid representation
remapgrid_t src_grid;
remapgrid_t dst_grid;
remapvars_t rv;

//CDO grid representation
griddes_t src_griddes;
griddes_t dst_griddes;


void print_remap_grid_info(remapgrid_t *grid);
void print_remap_info(remapvars_t *rv);

extern int remap_order;

//for printing the map type as string
const char *map_type_string[] = { "unknown", "conservative", "bilinear", "bicubic", "distwgt", "conserv_yac" };




#ifdef NULL
  #undef NULL
#endif
#define NULL (double*)0

double *src_grid_values = NULL;
double *dst_grid_values = NULL;

//gradients used for second order and bicubic remappings
//to be computed when the weights are read in or computed
double *src_grad1 = NULL;
double *src_grad2 = NULL;
double *src_grad3 = NULL;

int init_method = 0;

#define INIT_WEIGHTS_FILE 1
#define INIT_GRIDS 2





//for now I'm assuming conservative remapping is what we'll use
//other methods are supported when using weights files
int map_type = MAP_TYPE_CONSERV;
//int map_type = MAP_TYPE_BILINEAR;



int set_weights_file(char* filename) {
    strncpy(weights_file, filename, MAX_LENGTH);
    init_method = INIT_WEIGHTS_FILE;
    return 0;
}
int set_src_grid_file(char* filename) {
    strncpy(src_grid_file, filename, MAX_LENGTH);
    init_method = INIT_GRIDS;
    gridID1 = cdoDefineGrid(src_grid_file);
    return 0;
}
int set_dst_grid_file(char* filename) {
    strncpy(dst_grid_file, filename, MAX_LENGTH);
    init_method = INIT_GRIDS;
    gridID2 = cdoDefineGrid(dst_grid_file);
    return 0;
}
int get_weights_file(char** filename) {
    *filename = weights_file;
    return 0;
}
int get_src_grid_file(char** filename) {
    *filename = src_grid_file;
    return 0;
}
int get_dst_grid_file(char** filename) {
    *filename = dst_grid_file;
    return 0;
}



void read_weights_file() {
    int submap_type;
    int num_neighbors;

    read_remap_scrip(weights_file, -1, -1, &map_type, &submap_type, &num_neighbors,
              &remap_order, &src_grid, &dst_grid, &rv);
}

void init_cdo_griddata(int *gridID, griddes_t griddes, remapgrid_t grid) {

  gridInit(&griddes);

  griddes.size = grid.size;
  griddes.nvertex = grid.num_cell_corners;
  if (grid.rank == 1)  {
    griddes.type = GRID_UNSTRUCTURED;
  } else {
    griddes.type = GRID_CURVILINEAR;
    if (griddes.nvertex != 4) {
      fprintf(stderr, "curvilinear grid with %d corners unsupported",
               griddes.nvertex);
    }
  }

  griddes.xsize = grid.dims[0];
  griddes.ysize = grid.dims[1];

  griddes.prec = DATATYPE_FLT64;

  griddes.xvals   = (double*) malloc(griddes.size*sizeof(double));
  memcpy(griddes.xvals, grid.cell_center_lon, grid.size * sizeof(double));
  griddes.yvals   = (double*) malloc(griddes.size*sizeof(double));
  memcpy(griddes.yvals, grid.cell_center_lat, grid.size * sizeof(double));

  griddes.xbounds = (double*) malloc(griddes.nvertex*griddes.size*sizeof(double));
  memcpy(griddes.xbounds, grid.cell_corner_lon, grid.num_cell_corners*grid.size*sizeof(double));
  griddes.ybounds = (double*) malloc(griddes.nvertex*griddes.size*sizeof(double));
  memcpy(griddes.ybounds, grid.cell_corner_lat, grid.num_cell_corners*grid.size*sizeof(double));

  strncpy(griddes.xunits, "radians", MAX_LENGTH);
  strncpy(griddes.yunits, "radians", MAX_LENGTH);

  griddes.area = NULL;
  griddes.mask = (int *)0;

  *gridID = gridDefine(griddes);

  grid.lneed_cell_corners = TRUE;
  grid.luse_cell_corners = TRUE;

}


void compute_weights() {

  remap_vars_init(map_type, src_grid.size, dst_grid.size, &rv);

  rv.map_type = map_type;
  rv.norm_opt = NORM_OPT_FRACAREA;

  //calling CDO to compute the weights
  if (map_type == MAP_TYPE_CONSERV) {
    scrip_remap_weights_conserv(&src_grid, &dst_grid, &rv);
  } else if (map_type == MAP_TYPE_BILINEAR) {
    scrip_remap_weights_bilinear(&src_grid, &dst_grid, &rv);
  }

  //reduce the allocated space for storing the links
  if ( map_type == MAP_TYPE_CONSERV && rv.num_links != rv.max_links) {
    resize_remap_vars(&rv, rv.num_links-rv.max_links);
  }

}


int initialize_code() {
  cdoVerbose = TRUE;

  gridID1 = -1;
  src_grid.size = 0;
  src_grid.num_cell_corners = 0;
  src_grid.cell_center_lon = NULL;
  src_grid.cell_center_lat = NULL;
  src_grid.cell_corner_lon = NULL;
  src_grid.cell_corner_lat = NULL;

  gridID2 = -1;
  dst_grid.size = 0;
  dst_grid.num_cell_corners = 0;
  dst_grid.cell_center_lon = NULL;
  dst_grid.cell_center_lat = NULL;
  dst_grid.cell_corner_lon = NULL;
  dst_grid.cell_corner_lat = NULL;

  if      ( map_type == MAP_TYPE_CONSERV     ) rv.num_wts = 3;
  else if ( map_type == MAP_TYPE_CONSERV_YAC ) rv.num_wts = 1;
  else if ( map_type == MAP_TYPE_BILINEAR    ) rv.num_wts = 1;
  else if ( map_type == MAP_TYPE_BICUBIC     ) rv.num_wts = 4;
  else if ( map_type == MAP_TYPE_DISTWGT     ) rv.num_wts = 1;
  else cdoAbort("Unknown mapping method!");
  rv.pinit = FALSE;
  rv.max_links = 0;
  rv.num_links = 0;
  rv.src_cell_add = (int *)0;
  rv.tgt_cell_add = (int *)0;
  rv.wts = NULL;

  return 0;
}



int commit_parameters() {

  if (map_type == MAP_TYPE_CONSERV) {
    remap_order = 2;
  }

  //read either a weights file or grid files
  if (init_method == INIT_WEIGHTS_FILE) {
    read_weights_file();
  } else if (init_method == INIT_GRIDS) {

    if (gridID1 == -1) {
      init_cdo_griddata(&gridID1, src_griddes, src_grid);
    }
    if (gridID2 == -1) {
      init_cdo_griddata(&gridID2, dst_griddes, dst_grid);
    }

    remap_set_int(REMAP_NUM_SRCH_BINS, 720); //this is a setting that CDO also sets by default
    int remap_extrapolate = FALSE;
    rv.pinit = FALSE;
    remap_grids_init(map_type, remap_extrapolate, gridID1, &src_grid, gridID2, &dst_grid);

    int i=0;
    //by default CDO ignores the grid mask, we're doing this to reproduce the same results as CDO
    if (src_grid.mask) {
      for (i=0; i<src_grid.size; i++) {
        src_grid.mask[i] = 1;
      }
    }

    src_grid.store_link_fast = FALSE;
    src_grid.lextrapolate = FALSE;
    src_grid.non_global = FALSE;
    src_grid.lwrite_remap = FALSE;

    dst_grid.store_link_fast = FALSE;
    dst_grid.lextrapolate = FALSE;
    dst_grid.non_global = FALSE;
    dst_grid.lwrite_remap = FALSE;

    if (rv.num_links == 0) {
      compute_weights();
    } else {
      //weights have been set through the interface 
      //check for as far as we can if this has been done correctly
      if (rv.src_cell_add == (int *)0) {
        fprintf(stderr, "Warning: Remapper not properly initialized, set_num_links has been called but src_address seems to be missing\n");
      }
      if (rv.tgt_cell_add == (int *)0) {
        fprintf(stderr, "Warning: Remapper not properly initialized, set_num_links has been called but dst_address seems to be missing\n");
      }
      if (rv.wts == NULL) {
        fprintf(stderr, "Warning: Remapper not properly initialized, set_num_links has been called but weights seem to be missing\n");
      }

      //perform the rest of the initialization
      rv.max_links = rv.num_links;

      rv.sort_add = TRUE;
      rv.pinit = TRUE;
      rv.map_type = map_type;
      rv.norm_opt = NORM_OPT_FRACAREA;
      rv.resize_increment = (int) (0.1 * MAX(src_grid.size, dst_grid.size));

      rv.links.option    = FALSE;
      rv.links.max_links = 0;
      rv.links.num_blks  = 0;

    }
  } else {
    fprintf(stderr, "Remapper not properly initialized, set either a weights file or source and destination grid files\n");
    return 1;
  }

  //can not do this earlier since the grid sizes are not yet known
  src_grid_values = (double *)malloc(src_grid.size * sizeof(double));
  dst_grid_values = (double *)malloc(dst_grid.size * sizeof(double));

  return 0;
}




int cleanup_code() {

  if (src_grid_values) { free(src_grid_values); }
  if (dst_grid_values) { free(dst_grid_values); }

  if (src_grad1) { free(src_grad1); }
  if (src_grad2) { free(src_grad2); }
  if (src_grad3) { free(src_grad3); }

  return 0;
}


int set_src_grid_values(int *index_i, double *src_values, int n) {
  int i=0;
  for (i=0; i < n; i++) {
    src_grid_values[index_i[i]] = src_values[i];
  }
  return 0;
}
int set_dst_grid_values(int *index_i, double *dst_values, int n) {
  int i=0;
  for (i=0; i < n; i++) {
    dst_grid_values[index_i[i]] = dst_values[i];
  }
  return 0;
}
int get_src_grid_values(int *index_i, double *src_values, int n) {
  int i=0;
  for (i=0; i < n; i++) {
    src_values[i] = src_grid_values[index_i[i]];
  }
  return 0;
}
int get_dst_grid_values(int *index_i, double *dst_values, int n) {
  int i=0;
  for (i=0; i < n; i++) {
      dst_values[i] = dst_grid_values[index_i[i]];
  }
  return 0;
}



int set_src_grid_center_lat(int *index_i, double *src_center_lat, int n) {
  int i=0;
  if (src_grid.size == 0) {
    fprintf(stderr, "Error: set src grid size first\n");
  } else {
    if (src_grid.cell_center_lat == NULL) {
      src_grid.cell_center_lat = (double *)malloc(src_grid.size * sizeof(double));
    }
  }
  for (i=0; i < n; i++) {
    src_grid.cell_center_lat[index_i[i]] = src_center_lat[i];
  }
  return 0;
}
int set_dst_grid_center_lat(int *index_i, double *dst_center_lat, int n) {
  int i=0;
  if (dst_grid.size == 0) {
    fprintf(stderr, "Error: set dst grid size first\n");
  } else {
    if (dst_grid.cell_center_lat == NULL) {
      dst_grid.cell_center_lat = (double *)malloc(dst_grid.size * sizeof(double));
    }
  }
  for (i=0; i < n; i++) {
    dst_grid.cell_center_lat[index_i[i]] = dst_center_lat[i];
  }
  return 0;
}
int get_src_grid_center_lat(int *index_i, double *src_center_lat, int n) {
  int i=0;
  for (i=0; i < n; i++) {
    src_center_lat[i] = src_grid.cell_center_lat[index_i[i]];
  }
  return 0;
}
int get_dst_grid_center_lat(int *index_i, double *dst_center_lat, int n) {
  int i=0;
  for (i=0; i < n; i++) {
      dst_center_lat[i] = dst_grid.cell_center_lat[index_i[i]];
  }
  return 0;
}

int set_src_grid_center_lon(int *index_i, double *src_center_lon, int n) {
  int i=0;
  if (src_grid.size == 0) {
    fprintf(stderr, "Error: set src grid size first\n");
  } else {
    if (src_grid.cell_center_lon == NULL) {
      src_grid.cell_center_lon = (double *)malloc(src_grid.size * sizeof(double));
    }
  }
  for (i=0; i < n; i++) {
    src_grid.cell_center_lon[index_i[i]] = src_center_lon[i];
  }
  return 0;
}
int set_dst_grid_center_lon(int *index_i, double *dst_center_lon, int n) {
  int i=0;
  if (dst_grid.size == 0) {
    fprintf(stderr, "Error: set src grid size first\n");
  } else {
    if (dst_grid.cell_center_lon == NULL) {
      dst_grid.cell_center_lon = (double *)malloc(dst_grid.size * sizeof(double));
    }
  }
  for (i=0; i < n; i++) {
    dst_grid.cell_center_lon[index_i[i]] = dst_center_lon[i];
  }
  return 0;
}
int get_src_grid_center_lon(int *index_i, double *src_center_lon, int n) {
  int i=0;
  for (i=0; i < n; i++) {
    src_center_lon[i] = src_grid.cell_center_lon[index_i[i]];
  }
  return 0;
}
int get_dst_grid_center_lon(int *index_i, double *dst_center_lon, int n) {
  int i=0;
  for (i=0; i < n; i++) {
      dst_center_lon[i] = dst_grid.cell_center_lon[index_i[i]];
  }
  return 0;
}




int set_src_grid_corner_lat(int *index_i, double *src_corner_lat, int n) {
  int i=0;
  if (src_grid.size == 0 || src_grid.num_cell_corners == 0) {
    fprintf(stderr, "Error: set src grid size and corners first\n");
  } else {
    if (src_grid.cell_corner_lat == NULL) {
      src_grid.cell_corner_lat = (double *)malloc(src_grid.size * src_grid.num_cell_corners * sizeof(double));
    }
  }
  for (i=0; i < n; i++) {
    src_grid.cell_corner_lat[index_i[i]] = src_corner_lat[i];
  }
  return 0;
}
int set_dst_grid_corner_lat(int *index_i, double *dst_corner_lat, int n) {
  int i=0;
  if (dst_grid.size == 0 || dst_grid.num_cell_corners == 0) {
    fprintf(stderr, "Error: set dst grid size and corners first\n");
  } else {
    if (dst_grid.cell_corner_lat == NULL) {
      dst_grid.cell_corner_lat = (double *)malloc(dst_grid.size * dst_grid.num_cell_corners * sizeof(double));
    }
  }
  for (i=0; i < n; i++) {
    dst_grid.cell_corner_lat[index_i[i]] = dst_corner_lat[i];
  }
  return 0;
}
int get_src_grid_corner_lat(int *index_i, double *src_corner_lat, int n) {
  int i=0;
  for (i=0; i < n; i++) {
    src_corner_lat[i] = src_grid.cell_corner_lat[index_i[i]];
  }
  return 0;
}
int get_dst_grid_corner_lat(int *index_i, double *dst_corner_lat, int n) {
  int i=0;
  for (i=0; i < n; i++) {
      dst_corner_lat[i] = dst_grid.cell_corner_lat[index_i[i]];
  }
  return 0;
}

int set_src_grid_corner_lon(int *index_i, double *src_corner_lon, int n) {
  int i=0;
  if (src_grid.size == 0 || src_grid.num_cell_corners == 0) {
    fprintf(stderr, "Error: set src grid size and corners first\n");
  } else {
    if (src_grid.cell_corner_lon == NULL) {
      src_grid.cell_corner_lon = (double *)malloc(src_grid.size * src_grid.num_cell_corners * sizeof(double));
    }
  }
  for (i=0; i < n; i++) {
    src_grid.cell_corner_lon[index_i[i]] = src_corner_lon[i];
  }
  return 0;
}
int set_dst_grid_corner_lon(int *index_i, double *dst_corner_lon, int n) {
  int i=0;
  if (dst_grid.size == 0 || dst_grid.num_cell_corners == 0) {
    fprintf(stderr, "Error: set dst grid size and corners first\n");
  } else {
    if (dst_grid.cell_corner_lon == NULL) {
      dst_grid.cell_corner_lon = (double *)malloc(dst_grid.size * dst_grid.num_cell_corners * sizeof(double));
    }
  }
  for (i=0; i < n; i++) {
    dst_grid.cell_corner_lon[index_i[i]] = dst_corner_lon[i];
  }
  return 0;
}
int get_src_grid_corner_lon(int *index_i, double *src_corner_lon, int n) {
  int i=0;
  for (i=0; i < n; i++) {
    src_corner_lon[i] = src_grid.cell_corner_lon[index_i[i]];
  }
  return 0;
}
int get_dst_grid_corner_lon(int *index_i, double *dst_corner_lon, int n) {
  int i=0;
  for (i=0; i < n; i++) {
      dst_corner_lon[i] = dst_grid.cell_corner_lon[index_i[i]];
  }
  return 0;
}






int get_src_grid_size(int *size) {
  *size = src_grid.size;
  return 0;
}
int get_dst_grid_size(int *size) {
  *size = dst_grid.size;
  return 0;
}

int set_src_grid_size(int size) {
  init_method = INIT_GRIDS;
  src_grid.size = size;
  return 0;
}
int set_dst_grid_size(int size) {
  init_method = INIT_GRIDS;
  dst_grid.size = size;
  return 0;
}

/*
 * In CDO there are always 2 dimensions, in case of unstructured grid y=0
 */
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

int set_src_grid_dims(int x, int y) {
  init_method = INIT_GRIDS;
  src_grid.dims[0] = x;
  src_grid.dims[1] = y;
  if (y <= 1) {
    src_grid.rank = 1;
  } else {
    src_grid.rank = 2;
  }
  int y1 = max(1,y); //y is zero when grid is unstructured
  src_grid.size = x * y1;
  return 0;
}
int set_dst_grid_dims(int x, int y) {
  init_method = INIT_GRIDS;
  dst_grid.dims[0] = x;
  dst_grid.dims[1] = y;
  if (y <= 1) {
    dst_grid.rank = 1;
  } else {
    dst_grid.rank = 2;
  }
  int y1 = max(1,y); //y is zero when grid is unstructured
  dst_grid.size = x * y1;
  return 0;
}




int get_src_grid_corners(int *corners) {
  *corners = src_grid.num_cell_corners;
  return 0;
}
int get_dst_grid_corners(int *corners) {
  *corners = dst_grid.num_cell_corners;
  return 0;
}

int set_src_grid_corners(int corners) {
  src_grid.num_cell_corners = corners;
  return 0;
}
int set_dst_grid_corners(int corners) {
  dst_grid.num_cell_corners = corners;
  return 0;
}






int perform_remap() {

    /*
    This function calls the remapping function in CDO, it also performs some preprocessing and checks that
    would otherwise be performed by the Remap() function that implements the operator

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

    // if remap order is 2 we also need the gradiants of the src field
    if (remap_order == 2) {

        if (!src_grad1) {
            src_grad1 = (double *)malloc(src_grid.size * sizeof(double));
        }
        if (!src_grad2) {
            src_grad2 = (double *)malloc(src_grid.size * sizeof(double));
        }
        if (!src_grad3) {
            src_grad3 = (double *)malloc(src_grid.size * sizeof(double));
        }

        remap_gradients(src_grid, src_grid_values, src_grad1, src_grad2, src_grad3);
    }

  //this value is used if no links exist that contribute to the destination gridpoint in the remapping
  double missval = 0.0;

  remap(dst_grid_values, missval, dst_grid.size, rv.num_links, rv.wts,
       rv.num_wts, rv.tgt_cell_add, rv.src_cell_add, src_grid_values,
       src_grad1, src_grad2, src_grad3,
       rv.links);

  return 0;
}




int get_num_links(int *num_links) {
  *num_links = rv.num_links;
  return 0;
}
int get_remap_links(int *index_i, int *src_address, int *dst_address, double *weight1, double *weight2, double *weight3, int n) {
  int i;

  for (i=0; i<n; i++) {
    src_address[i] = rv.src_cell_add[index_i[i]];
    dst_address[i] = rv.tgt_cell_add[index_i[i]];
    weight1[i] = rv.wts[(index_i[i]*rv.num_wts)+0];
    if (rv.num_wts > 1) {
      weight2[i] = rv.wts[(index_i[i]*rv.num_wts)+1];
      weight3[i] = rv.wts[(index_i[i]*rv.num_wts)+2];
    }
  }

  return 0;
}

int set_num_links(int num_links) {
  rv.num_links = num_links;
  return 0;
}
int set_src_address(int *index_i, int *src_address, int n) {
  int i=0;
  if (rv.num_links == 0) {
    fprintf(stderr, "Error: set num links first\n");
  } else {
    if (rv.src_cell_add == (int *)0) {
      rv.src_cell_add = (int *)malloc(rv.num_links * sizeof(int));
    }
  }
  for (i=0; i < n; i++) {
    rv.src_cell_add[index_i[i]] = src_address[i];
  }
  return 0;
}
int set_dst_address(int *index_i, int *dst_address, int n) {
  int i=0;
  if (rv.num_links == 0) {
    fprintf(stderr, "Error: set num links first\n");
  } else {
    if (rv.tgt_cell_add == (int *)0) {
      rv.tgt_cell_add = (int *)malloc(rv.num_links * sizeof(int));
    }
  }
  for (i=0; i < n; i++) {
    rv.tgt_cell_add[index_i[i]] = dst_address[i];
  }
  return 0;
}
int set_weights(int *index_i, double *weight1, double *weight2, double *weight3, int n) {
  int i=0;
  if (rv.num_links == 0) {
    fprintf(stderr, "Error: set num links first\n");
  } else {
    if (rv.wts == NULL) {
      rv.wts = (double *)malloc(rv.num_links * rv.num_wts * sizeof(double));
    }
  }
  for (i=0; i < n; i++) {
    rv.wts[(index_i[i]*rv.num_wts)+0] = weight1[i];
    if (rv.num_wts > 1) {
        rv.wts[(index_i[i]*rv.num_wts)+1] = weight2[i];
        rv.wts[(index_i[i]*rv.num_wts)+2] = weight3[i];
    }
  }
  return 0;
}




/* ------------- functions below this line are mainly for debugging purposes ----------- */



/* This function does the same as get_remap_links but also tries
 * to verify that a link exists for each destination grid cell
 * when all the remap links are requested
 */
int check_remap_links(int *index_i, int *src_address, int *dst_address, double *weights, int n) {
  int i;

  for (i=0; i<n; i++) {
    src_address[i] = rv.src_cell_add[index_i[i]];
    dst_address[i] = rv.tgt_cell_add[index_i[i]];
    weights[i] = rv.wts[index_i[i] * rv.num_wts];
  }

  if (n == rv.num_links) {
    for (i=0; i<dst_grid.size; i++) {
      dst_grid_values[i] = 2.0;
    }
    for (i=0; i<n; i++) {
      dst_grid_values[dst_address[i]] = 0.0;
    }
    int error = 0;
    for (i=0; i<dst_grid.size; i++) {
      if (dst_grid_values[i] == 2.0) {
        error++;
      }
    }
    if (error > 0) {
      printf("number of missing destination grid cells in remapping = %d\n", error);
    }
  }
  return 0;
}













int print_info() {
  printf("source grid:\n");
  print_remap_grid_info(&src_grid);

  printf("destination grid:\n");
  print_remap_grid_info(&dst_grid);

  printf("remap info:\n");
  print_remap_info(&rv);
  
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

  printf("mask %p\n",    grid->mask);                  /* flag which cells participate */

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
  printf("cell_bound_box = %p\n", grid->cell_bound_box);

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

  int i;
  int min_i=dst_grid.size;
  int max_i=0;
  for (i=0; i < rv->num_links; i++) {
    int dst_i = rv->tgt_cell_add[i];
    if (dst_i < min_i) min_i = dst_i;
    if (dst_i > max_i) max_i = dst_i;
  }
  printf("min destination address = %d, max destination address = %d\n", min_i, max_i);

//  double*  wts;              // map weights for each link [max_links*num_wts] 
  printf("rv->wts %p\n", rv->wts);

/*
  int j;
  printf("first 10 weights (assuming array of structs):\n");
  for (i=0; i < 10 ; i++) {
    for (j=0; j < rv->num_wts; j++) {
      printf("%f ", rv->wts[i*3+j]);
    }
    printf("\n");
  }
*/

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
