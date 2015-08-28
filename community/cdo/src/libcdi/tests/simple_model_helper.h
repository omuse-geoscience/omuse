#ifndef SIMPLE_MODEL_HELPER_H
#define SIMPLE_MODEL_HELPER_H

#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <inttypes.h>
#include <math.h>
#include <time.h>

void
var_scale(int datatype, double *mscale, double *mrscale);

static inline double
sign_flat(double v)
{
  if (!(v < 0.0 || v > 0.0)) return 0.0;
  //    if (v == 0.0) return 0.0;
  return v;
}

/* data generator function */
static inline double
dg_wobble(double frac_x, double frac_y, double mscale, double mrscale)
{
  double r = sign_flat(round((cos(2.0 * M_PI * frac_x)
                              * sin(2.0 * M_PI * frac_y)) * mscale)) * mrscale;
  return r;
}

time_t
cditime2time_t(int date, int timeofday);
void
time_t2cditime(time_t t, int *date, int *timeofday);

#if defined (USE_MPI) && ! defined(HAVE_PPM_CORE)
struct PPM_extent
{
  int32_t first, size;
};

struct PPM_extent
PPM_uniform_partition(struct PPM_extent set_interval, int nparts,
                      int part_idx);

int
PPM_prime_factorization_32(uint32_t n, uint32_t **factors);

#endif

#endif
