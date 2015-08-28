#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>

#include "cdi.h"
#include "dmemory.h"

#include "simple_model_helper.h"

void
var_scale(int datatype, double *mscale, double *mrscale)
{
  int mant_bits;
  switch (datatype)
    {
    case DATATYPE_PACK8:
      mant_bits = 7;
      break;
    case DATATYPE_PACK16:
      mant_bits = 15;
      break;
    case DATATYPE_PACK24:
      mant_bits = 23;
      break;
    case DATATYPE_FLT32:
      mant_bits = 24;
      break;
    case DATATYPE_FLT64:
      mant_bits = 53;
      break;
    case DATATYPE_INT8:
    case DATATYPE_INT16:
    case DATATYPE_INT32:
    default:
      fprintf(stderr, "Unexpected or unusable content format: %d\n",
              datatype);
      exit(EXIT_FAILURE);
    }
  *mscale = (double)(INT64_C(1) << mant_bits);
  *mrscale = 1.0 / *mscale;
}

/**
 * Compute UNIX epoch-based time_t from CDI's decimal encoding of date.
 */
time_t
cditime2time_t(int date, int timeofday)
{
  struct tm t_s;
  time_t t;
  t_s.tm_year = date / 10000;
  t_s.tm_mon = (date - t_s.tm_year * 10000)/100;
  t_s.tm_mday = date % 100;
  t_s.tm_year -= 1900;
  t_s.tm_hour = timeofday/10000;
  t_s.tm_min = (timeofday%10000)/100;
  t_s.tm_sec = timeofday%100;
  t_s.tm_isdst = 0;
  t = mktime(&t_s);
  return t;
}

/**
 * Build decimal encoding of date from UNIX epoch-based time_t.
 */
void
time_t2cditime(time_t t, int *date, int *timeofday)
{
  struct tm *t_s;
  t_s = localtime(&t);
  *date = (t_s->tm_year + 1900) * 10000 + t_s->tm_mon * 100 + t_s->tm_mday;
  *timeofday = t_s->tm_hour * 10000 + t_s->tm_min * 100 + t_s->tm_sec;
}

#if defined USE_MPI && ! defined HAVE_PPM_CORE

static int32_t
uniform_partition_start(struct PPM_extent set_interval, int nparts,
                        int part_idx);

struct PPM_extent
PPM_uniform_partition(struct PPM_extent set_interval, int nparts,
                      int part_idx)
{
  struct PPM_extent range;
  range.first = uniform_partition_start(set_interval, nparts, part_idx);
  range.size = uniform_partition_start(set_interval, nparts, part_idx + 1)
    - range.first;
  return range;
}

static int32_t
uniform_partition_start(struct PPM_extent set_interval, int nparts,
                        int part_idx)
{
  int32_t part_offset
    = ((int64_t)set_interval.size * (int64_t)part_idx) / (int64_t)nparts;
  int32_t start = set_interval.first + part_offset;
  return start;
}

int
PPM_prime_factorization_32(uint32_t n, uint32_t **factors)
{
  if (n <= 1) return 0;
  uint32_t * restrict pfactors = xrealloc(*factors, 32 * sizeof (pfactors[0]));
  size_t numFactors = 0;
  uint32_t unfactored = n;
  while (!(unfactored & 1))
  {
    pfactors[numFactors] = 2;
    ++numFactors;
    unfactored >>= 1;
  }
  uint32_t divisor = 3;
  while (unfactored > 1)
  {
    while (unfactored % divisor == 0)
    {
      unfactored /= divisor;
      pfactors[numFactors] = divisor;
      ++numFactors;
    }
    divisor += 1;
  }
  *factors = pfactors;
  return numFactors;
}

#endif

