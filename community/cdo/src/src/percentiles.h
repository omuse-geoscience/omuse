/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2006 Brockmann Consult
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#ifndef PERCENTILES_H_
#define PERCENTILES_H_

#include "field.h"

typedef struct {
  double min;
  double max;
  double step;
  int    nbins;
  int    nsamp;
  void  *ptr;
}
HISTOGRAM;

typedef struct {
  int    nvars;
  int   *nlevels;
  int   *grids;
  HISTOGRAM ***histograms;
}
HISTOGRAM_SET;


HISTOGRAM_SET *hsetCreate(int nvars);
void hsetCreateVarLevels(HISTOGRAM_SET *hset, int varID, int nlevels, int nhists);
void hsetDestroy(HISTOGRAM_SET *hset);

void hsetDefVarLevelBounds(HISTOGRAM_SET *hset, int varID, int levelID, const field_t *min, const field_t *max);
void hsetAddVarLevelValues(HISTOGRAM_SET *histField, int varID, int levelID, const field_t *field);
void hsetGetVarLevelPercentiles(field_t *field, const HISTOGRAM_SET *hset, int varID, int levelID, double pn); 

#endif /*PERCENTILES_H_*/
