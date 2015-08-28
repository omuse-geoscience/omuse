#ifndef _REMAP_STORE_LINK_H
#define _REMAP_STORE_LINK_H


typedef struct
{
  int    add;
  double weight;
} addweight_t;

typedef struct
{
  int    add;
  double weight[4];
} addweight4_t;

typedef struct {
  int nlinks;
  int offset;
  addweight_t *addweights;
} weightlinks_t;

typedef struct {
  int nlinks;
  int offset;
  addweight4_t *addweights;
} weightlinks4_t;


void store_weightlinks(long num_weights, int *srch_add, double *weights, long cell_add, weightlinks_t *weightlinks);
void store_weightlinks4(long num_weights, int *srch_add, double weights[4][4], long cell_add, weightlinks4_t *weightlinks);
void weightlinks2remaplinks(long tgt_grid_size, weightlinks_t *weightlinks, remapvars_t *rv);
void weightlinks2remaplinks4(long tgt_grid_size, weightlinks4_t *weightlinks, remapvars_t *rv);
void sort_add_and_wgts(size_t num_weights, int *src_add, double *wgts);
void sort_add_and_wgts4(size_t num_weights, int *src_add, double wgts[4][4]);


#endif  /* _REMAP_STORE_LINK */
