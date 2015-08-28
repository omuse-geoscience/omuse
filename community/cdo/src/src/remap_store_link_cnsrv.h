#ifndef _REMAP_STORE_LINK_CNSRV_H
#define _REMAP_STORE_LINK_CNSRV_H

extern int  remap_store_link_fast;

/* used for store_link_fast */

#define BLK_SIZE 4096
#define BLK_NUM(x) (x/grid_store->blk_size)
#define BLK_IDX(x) (x%grid_store->blk_size)
struct grid_layer
{
  int *grid2_link;
  struct grid_layer *next;
};

typedef struct grid_layer grid_layer_t;

typedef struct
{
  int blk_size;
  int max_size;
  int nblocks;
  int *blksize;
  int *nlayers;
  grid_layer_t **layers;
} grid_store_t;

void grid_store_init(grid_store_t *grid_store, long gridsize);
void grid_store_delete(grid_store_t *grid_store);

void store_link_cnsrv_fast(remapvars_t *rv, long add1, long add2, long num_wts, double *weights, grid_store_t *grid_store);
void store_link_cnsrv(remapvars_t *rv, long add1, long add2, double *restrict weights, int *link_add1[2], int *link_add2[2]);

#endif  /* _REMAP_STORE_LINK_CNSRV */
