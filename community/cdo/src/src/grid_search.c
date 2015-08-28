#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "grid_search.h"

#if defined(KDTREE_FAST)
#include "kdtree_fast.h"
#else
#include "kdtree.h"
#endif

static inline void LLtoXYZ(double lon, double lat, double p_out[])
{
   double cos_lat = cos(lat);
   p_out[0] = cos_lat * cos(lon);
   p_out[1] = cos_lat * sin(lon);
   p_out[2] = sin(lat);
}


struct gridsearch {
  unsigned n;

#if defined(KDTREE_FAST)
  struct kd_node_t *kdarray;
  struct kd_node_t *kdt;
#else
  struct kdtree *kdt;
#endif
};


struct gridsearch *gridsearch_index_create(unsigned n, const double *restrict lons, const double *restrict lats, const unsigned *restrict index)
{
  struct gridsearch *gs = (struct gridsearch *) malloc(sizeof(struct gridsearch));
  
#if defined(KDTREE_FAST)
  struct kd_node_t *kdarray = calloc(n, sizeof(struct kd_node_t));

#pragma omp simd
  for ( unsigned i = 0; i < n; i++ ) 
    {
      LLtoXYZ(lons[index[i]], lats[index[i]], kdarray[i].xyz);
      kdarray[i].index = index[i];
    }

  gs->kdarray = kdarray;
  gs->kdt = make_tree(kdarray, n, 0, 3);

#else
  gs->kdt = kd_create(3);
  double pos[3];

  for ( unsigned i = 0; i < n; i++ ) 
    {
      LLtoXYZ(lons[index[i]], lats[index[i]], pos);
      kd_insert(gs->kdt, pos, index[i]);
    }

#endif

  return gs;
}


struct gridsearch *gridsearch_create(unsigned n, const double *restrict lons, const double *restrict lats)
{
  struct gridsearch *gs = (struct gridsearch *) malloc(sizeof(struct gridsearch));
  
#if defined(KDTREE_FAST)
  struct kd_node_t *kdarray = calloc(n, sizeof(struct kd_node_t));

#pragma omp simd
  for ( unsigned i = 0; i < n; i++ ) 
    {
      LLtoXYZ(lons[i], lats[i], kdarray[i].xyz);
      kdarray[i].index = i;
    }

  gs->kdarray = kdarray;
  gs->kdt = make_tree(kdarray, n, 0, 3);

#else
  gs->kdt = kd_create(3);
  double pos[3];

  for ( unsigned i = 0; i < n; i++ ) 
    {
      LLtoXYZ(lons[i], lats[i], pos);
      kd_insert(gs->kdt, pos, i);
    }

#endif

  return gs;
}


void gridsearch_delete(struct gridsearch *gs)
{
#if defined(KDTREE_FAST)
  free(gs->kdarray);
#else
  kd_free(gs->kdt);
#endif
  
  free(gs);
}


void *gridsearch_nearest(struct gridsearch *gs, double lon, double lat)
{
  double pos[3];
  LLtoXYZ(lon, lat, pos);

#if defined(KDTREE_FAST)
  return (void *) kd_nearest(gs->kdt, pos);
#else
  return (void *) kd_nearest(gs->kdt, pos);
#endif
}


unsigned gridsearch_item(void *gs_result)
{
  unsigned index = 0;

#if defined(KDTREE_FAST)
  index = kd_item((struct kd_node_t *) gs_result);
#else
  index = kd_res_item((struct kdres *) gs_result, NULL);
#endif
  
  return index;
}
