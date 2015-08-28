#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "kdtree_fast.h"

static inline double dist(const double *restrict a_xyz, const double *restrict b_xyz, int dim)
{
  double t, d = 0;
  while (dim--) {
    t = a_xyz[dim] - b_xyz[dim];
    d += t * t;
  }
  return d;
}

static inline void swap(struct kd_node_t *x, struct kd_node_t *y)
{
  double tmp[MAX_DIM];
  memcpy(tmp,    x->xyz, sizeof(tmp));
  memcpy(x->xyz, y->xyz, sizeof(tmp));
  memcpy(y->xyz, tmp,    sizeof(tmp));
  unsigned index = x->index;
  x->index = y->index;
  y->index = index;
}
 
/* see quickselect method */
struct kd_node_t *find_median(struct kd_node_t *start, struct kd_node_t *end, int idx)
{
	if (end <= start) return NULL;
	if (end == start + 1) return start;
 
	struct kd_node_t *p, *store, *md = start + (end - start) / 2;
	double pivot;
	while (1) {
		pivot = md->xyz[idx];
 
		swap(md, end - 1);
		for (store = p = start; p < end; p++) {
			if (p->xyz[idx] < pivot) {
				if (p != store)
					swap(p, store);
				store++;
			}
		}
		swap(store, end - 1);
 
		/* median has duplicate values */
		if (store->xyz[idx] == md->xyz[idx])
			return md;
 
		if (store > md)	end = store;
		else		start = store;
	}
}
 
struct kd_node_t *make_tree(struct kd_node_t *t, int len, int i, int dim)
{
	struct kd_node_t *n;
 
	if (!len) return 0;
 
	if ((n = find_median(t, t + len, i))) {
		i = (i + 1) % dim;
		n->left  = make_tree(t, n - t, i, dim);
		n->right = make_tree(n + 1, t + len - (n + 1), i, dim);
	}
	return n;
}
 
/* global variable, so sue me */
int visited;
 
void nearest(struct kd_node_t *root, double *xyz, int i, int dim,
             struct kd_node_t **best, double *best_dist)
{
  double d, dx, dx2;
 
  if (!root) return;
  d = dist(root->xyz, xyz, dim);
  dx = root->xyz[i] - xyz[i];
  dx2 = dx * dx;
 
  visited ++;
 
  if (!*best || d < *best_dist) {
    *best_dist = d;
    *best = root;
  }
 
  /* if chance of exact match is high */
  if (!*best_dist) return;
 
  if ( ++i >= dim ) i = 0;
 
  nearest(dx > 0 ? root->left : root->right, xyz, i, dim, best, best_dist);
  if ( dx2 >= *best_dist ) return;
  nearest(dx > 0 ? root->right : root->left, xyz, i, dim, best, best_dist);
}


void *kd_nearest(struct kd_node_t *root, double *xyz)
{
  struct kd_node_t *found = NULL;
  double dist;

  nearest(root, xyz, 0, 3, &found, &dist);

  return (void *) found;
}


unsigned kd_item(struct kd_node_t *result)
{
  if ( result )
    return result->index;
  else
    return 0;
}


#ifdef XXX
#define N 1000000
#define rand1()	(rand() / (double)RAND_MAX)
#define rand_pt(v) { v.xyz[0] = rand1(); v.xyz[1] = rand1(); v.xyz[2] = rand1(); }
int main(void)
{
	int i;
	struct kd_node_t wp[] = {
		{{2, 3}}, {{5, 4}}, {{9, 6}}, {{4, 7}}, {{8, 1}}, {{7, 2}}
	};
	struct kd_node_t this = {{9, 2}};
	struct kd_node_t *root, *found, *million;
	double best_dist;
 
	root = make_tree(wp, sizeof(wp) / sizeof(wp[1]), 0, 2);
 
	visited = 0;
	found = 0;
	nearest(root, &this, 0, 2, &found, &best_dist);
 
	printf(">> WP tree\nsearching for (%g, %g)\n"
		"found (%g, %g) dist %g\nseen %d nodes\n\n",
		this.xyz[0], this.xyz[1],
		found->xyz[0], found->xyz[1], sqrt(best_dist), visited);
 
	million = calloc(N, sizeof(struct kd_node_t));
	srand(time(0));
	for (i = 0; i < N; i++) rand_pt(million[i]);
 
	root = make_tree(million, N, 0, 3);
	rand_pt(this);
 
	visited = 0;
	found = 0;
	nearest(root, &this, 0, 3, &found, &best_dist);
 
	printf(">> Million tree\nsearching for (%g, %g, %g)\n"
		"found (%g, %g, %g) dist %g\nseen %d nodes\n",
		this.xyz[0], this.xyz[1], this.xyz[2],
		found->xyz[0], found->xyz[1], found->xyz[2],
		sqrt(best_dist), visited);
 
	/* search many random points in million tree to see average behavior.
	   tree size vs avg nodes visited:
	   	10		~  7
	   	100		~ 16.5
		1000		~ 25.5
		10000		~ 32.8
		100000		~ 38.3
		1000000		~ 42.6
		10000000	~ 46.7				*/
	int sum = 0, test_runs = 100000;
	for (i = 0; i < test_runs; i++) {
		found = 0;
		visited = 0;
		rand_pt(this);
		nearest(root, &this, 0, 3, &found, &best_dist);
		sum += visited;
	}
	printf("\n>> Million tree\n"
		"visited %d nodes for %d random findings (%f per lookup)\n",
		sum, test_runs, sum/(double)test_runs);
 
	// free(million);
 
	return 0;
}
#endif
