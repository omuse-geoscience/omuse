/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2015 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#include "cdo.h"
#include "cdo_int.h"
#include <cdi.h>


void farcfun(field_t *field, double rconst, int function)
{
  if      ( function == func_add ) farcadd(field, rconst);
  else if ( function == func_sub ) farcsub(field, rconst);
  else if ( function == func_mul ) farcmul(field, rconst);
  else if ( function == func_div ) farcdiv(field, rconst);
  else if ( function == func_mod ) farmod(field, rconst);
  else    cdoAbort("%s: function %d not implemented!", __func__, function);
}

void farcmul(field_t *field, double rconst)
{
  int i, len;
  int    nwpv     = field->nwpv;
  int    grid     = field->grid;
  int    nmiss    = field->nmiss;
  double missval1 = field->missval;
  double missval2 = field->missval;
  double *array   = field->ptr;

  if ( nwpv != 2 ) nwpv = 1;

  len    = nwpv*gridInqSize(grid);

  if ( nmiss > 0 )
    {
      for ( i = 0; i < len; i++ ) 
	array[i] = MUL(array[i], rconst);
    }
  else
    {
      /*
#if defined(_OPENMP)
#pragma omp parallel for default(shared) private(i)
#endif
      */
      for ( i = 0; i < len; i++ ) 
	array[i] *= rconst;
    }
}


void farcdiv(field_t *field, double rconst)
{
  int i, len;
  int    grid     = field->grid;
  int    nmiss    = field->nmiss;
  double missval1 = field->missval;
  double missval2 = field->missval;
  double *array   = field->ptr;

  len    = gridInqSize(grid);

  if ( nmiss > 0 || IS_EQUAL(rconst, 0) )
    {
      for ( i = 0; i < len; i++ )
	array[i] = DIV(array[i], rconst);

      if ( IS_EQUAL(rconst, 0) ) field->nmiss = len;
    }
  else
    {
      for ( i = 0; i < len; i++ ) 
	array[i] /= rconst;
    }
}


void farcadd(field_t *field, double rconst)
{
  int i, len;
  int    grid     = field->grid;
  int    nmiss    = field->nmiss;
  double missval1 = field->missval;
  double missval2 = field->missval;
  double *array   = field->ptr;

  len    = gridInqSize(grid);

  if ( nmiss > 0 )
    {
      for ( i = 0; i < len; i++ ) 
	array[i] = ADD(array[i], rconst);
    }
  else
    {
      for ( i = 0; i < len; i++ ) 
	array[i] += rconst;
    }
}


void farcsub(field_t *field, double rconst)
{
  farcadd(field, -rconst);
}


void farinv(field_t *field)
{
  int i, len;
  int    grid     = field->grid;
  double missval1 = field->missval;
  double missval2 = field->missval;
  double *array   = field->ptr;

  len    = gridInqSize(grid);

  for ( i = 0; i < len; i++ ) 
    array[i] = DIV(1.0, array[i]);

  field->nmiss = 0;
  for ( i = 0; i < len; i++ )
    if ( DBL_IS_EQUAL(array[i], missval1) ) field->nmiss++;
}


void farmod(field_t *field, double divisor)
{
  int i, len;
  int    grid     = field->grid;
  double missval1 = field->missval;
  double *array   = field->ptr;

  len    = gridInqSize(grid);

  for ( i = 0; i < len; i++ )
    {
      array[i] = DBL_IS_EQUAL(array[i], missval1) ? missval1 : fmod(array[i], divisor);
    }
}
