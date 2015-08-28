#ifndef _KDTREE_FAST_H_
#define _KDTREE_FAST_H_
 
#define MAX_DIM 3

struct kd_node_t {
  double xyz[MAX_DIM];
  struct kd_node_t *left, *right;
  unsigned index;
};

struct kd_node_t *make_tree(struct kd_node_t *t, int len, int i, int dim);
void nearest(struct kd_node_t *root, double *xyz, int i, int dim,
             struct kd_node_t **best, double *best_dist);
void *kd_nearest(struct kd_node_t *root, double *xyz);
unsigned kd_item(struct kd_node_t *result);

#endif
