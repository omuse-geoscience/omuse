#ifndef _MERGE_SORT2_H_
#define _MERGE_SORT2_H_

/* MERGE SORT DEFINES */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "cdo.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

void sort_iter_single(long num_links, double *restrict add1, int parent);

#endif
