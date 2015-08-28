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

#if 0
#ifdef __SSE2__
#include <emmintrin.h>
#endif
#endif

void farfun(field_t *field1, field_t field2, const int function)
{
  if      ( function == func_add   ) faradd(field1, field2);
  else if ( function == func_min   ) farmin(field1, field2);
  else if ( function == func_max   ) farmax(field1, field2);
  else if ( function == func_sum   ) farsum(field1, field2);
  else if ( function == func_mean  ) farsum(field1, field2);
  else if ( function == func_avg   ) faradd(field1, field2);
  else if ( function == func_sub   ) farsub(field1, field2);
  else if ( function == func_mul   ) farmul(field1, field2);
  else if ( function == func_div   ) fardiv(field1, field2);
  else if ( function == func_atan2 ) faratan2(field1, field2);
  else cdoAbort("%s: function %d not implemented!", __func__, function);
}

static
void arradd(const size_t n, double * restrict a, const double * restrict b)
{
  size_t i;
 
  // SSE2 version is 15% faster than the original loop (tested with gcc47)
#if 0
  //#ifdef __SSE2__ /*__SSE2__*/ // bug in this code!!!
  const size_t residual =  n % 8;
  const size_t ofs = n - residual;

  __m128d *av = (__m128d *) a; // assume 16-byte aligned
  __m128d *bv = (__m128d *) b; // assume 16-byte aligned
  for ( i = 0; i < n/2; i+=4 )
    {
      av[i  ] = _mm_add_pd(av[i  ], bv[i  ]);
      av[i+1] = _mm_add_pd(av[i+1], bv[i+1]);
      av[i+2] = _mm_add_pd(av[i+2], bv[i+2]);
      av[i+3] = _mm_add_pd(av[i+3], bv[i+3]);
    }
  printf("residual, ofs, n %ld %ld %ld\n", residual, ofs, n);
  for ( i = 0; i < residual; i++ )  a[ofs+i] += b[ofs+i];

#else

  for ( i = 0; i < n; i++ ) a[i] += b[i];

#endif
}

static
void arraddw(const size_t n, double * restrict a, const double * restrict b, double w)
{
  for ( size_t i = 0; i < n; i++ ) a[i] += w*b[i];
}


void faradd(field_t *field1, field_t field2)
{
  size_t   i, len;
  int          nwpv     = field1->nwpv;
  const int    grid1    = field1->grid;
  const int    nmiss1   = field1->nmiss;
  const double missval1 = field1->missval;
  double *array1  = field1->ptr;
  /*  double *weight1 = field1->weight; */
  const int    grid2    = field2.grid;
  const int    nmiss2   = field2.nmiss;
  const double missval2 = field2.missval;
  double *array2  = field2.ptr;

  if ( nwpv != 2 ) nwpv = 1;

  len = (size_t) (nwpv*gridInqSize(grid1));

  if ( len != (size_t) (nwpv*gridInqSize(grid2)) )
    cdoAbort("Fields have different gridsize (%s)", __func__);

  if ( nmiss1 > 0 || nmiss2 > 0 )
    {
      for ( i = 0; i < len; i++ ) 
	array1[i] = ADD(array1[i], array2[i]);

      field1->nmiss = 0;
      for ( i = 0; i < len; i++ )
	if ( DBL_IS_EQUAL(array1[i], missval1) ) field1->nmiss++;
    }
  else
    {
      arradd(len, array1, array2);
    }
}


void farsum(field_t *field1, field_t field2)
{
  size_t   i, len;
  int          nwpv     = field1->nwpv;
  const int    grid1    = field1->grid;
  const int    nmiss1   = field1->nmiss;
  const double missval1 = field1->missval;
  double *array1  = field1->ptr;
  /*  double *weight1 = field1->weight; */
  const int    grid2    = field2.grid;
  const int    nmiss2   = field2.nmiss;
  const double missval2 = field2.missval;
  double *array2  = field2.ptr;

  if ( nwpv != 2 ) nwpv = 1;

  len = (size_t) (nwpv*gridInqSize(grid1));

  if ( len != (size_t) (nwpv*gridInqSize(grid2)) )
    cdoAbort("Fields have different gridsize (%s)", __func__);

  if ( nmiss1 > 0 || nmiss2 > 0 )
    {
      for ( i = 0; i < len; i++ )
	if ( !DBL_IS_EQUAL(array2[i], missval2) )
	  {
	    if ( !DBL_IS_EQUAL(array1[i], missval1) )
	      array1[i] += array2[i];
	    else
	      array1[i] = array2[i];
	  }

      field1->nmiss = 0;
      for ( i = 0; i < len; i++ )
	if ( DBL_IS_EQUAL(array1[i], missval1) ) field1->nmiss++;
    }
  else
    {
      arradd(len, array1, array2);
    }
}


void farsumw(field_t *field1, field_t field2, double w)
{
  size_t   i, len;
  int          nwpv     = field1->nwpv;
  const int    grid1    = field1->grid;
  const int    nmiss1   = field1->nmiss;
  const double missval1 = field1->missval;
  double *array1  = field1->ptr;
  /*  double *weight1 = field1->weight; */
  const int    grid2    = field2.grid;
  const int    nmiss2   = field2.nmiss;
  const double missval2 = field2.missval;
  double *array2  = field2.ptr;

  if ( nwpv != 2 ) nwpv = 1;

  len = (size_t) (nwpv*gridInqSize(grid1));

  if ( len != (size_t) (nwpv*gridInqSize(grid2)) )
    cdoAbort("Fields have different gridsize (%s)", __func__);

  if ( nmiss1 > 0 || nmiss2 > 0 )
    {
      for ( i = 0; i < len; i++ )
	if ( !DBL_IS_EQUAL(array2[i], missval2) )
	  {
	    if ( !DBL_IS_EQUAL(array1[i], missval1) )
	      array1[i] += w*array2[i];
	    else
	      array1[i] = w*array2[i];
	  }

      field1->nmiss = 0;
      for ( i = 0; i < len; i++ )
	if ( DBL_IS_EQUAL(array1[i], missval1) ) field1->nmiss++;
    }
  else
    {
      arraddw(len, array1, array2, w);
    }
}

/* 
 * Compute the occurrence of values in field, if they do not equal refval.
 * This can be used to compute the lengths of multiple periods in a timeseries.
 * Missing field values are handled like refval, i.e. they stop a running
 * period. If there is missing data in the occurence field, missing fields
 * values do not change anything (they do not start a non-period by setting
 * occurrence to zero).
 */
void farsumtr(field_t *occur, field_t field, const double refval)
{
  size_t   i, len;
  double omissval = occur->missval;
  double  *oarray = occur->ptr;
  double fmissval = field.missval;
  double  *farray = field.ptr;

  len    = (size_t) gridInqSize(occur->grid);

  if ( len != (size_t) gridInqSize(field.grid) )
    cdoAbort("Fields have different gridsize (%s)", __func__);

  if ( occur->nmiss > 0 || field.nmiss > 0 )
    {
#if defined(_OPENMP)
#pragma omp parallel for default(shared) schedule(static)
#endif
      for ( i = 0; i < len; i++ )
	if ( !DBL_IS_EQUAL(farray[i], fmissval) )
	  {
	    if ( !DBL_IS_EQUAL(oarray[i], omissval) )
	      oarray[i] = (DBL_IS_EQUAL(farray[i], refval)) ? 0.0 : oarray[i] + 1.0;
	    else
	      oarray[i] = (DBL_IS_EQUAL(farray[i], refval)) ? 0.0 : 1.0;
	  }
	else
	{
	  if ( !DBL_IS_EQUAL(oarray[i], omissval) )
	    oarray[i] = 0.0;
	}

      occur->nmiss = 0;
      for ( i = 0; i < len; i++ )
	if ( DBL_IS_EQUAL(oarray[i], omissval) ) occur->nmiss++;
    }
  else
    {
#if defined(_OPENMP)
#pragma omp parallel for default(shared)
#endif
      for ( i = 0; i < len; i++ ) 
	oarray[i] = (DBL_IS_EQUAL(farray[i], refval)) ? 0.0 : oarray[i] + 1.0;
    }
}


void farsumq(field_t *field1, field_t field2)
{
  size_t   i, len;
  int          nwpv     = field1->nwpv;
  const int    grid1    = field1->grid;
  const int    nmiss1   = field1->nmiss;
  const double missval1 = field1->missval;
  double *array1  = field1->ptr;
  /*  double *weight1 = field1->weight; */
  const int    grid2    = field2.grid;
  const int    nmiss2   = field2.nmiss;
  const double missval2 = field2.missval;
  double *array2  = field2.ptr;

  if ( nwpv != 2 ) nwpv = 1;

  len = (size_t) (nwpv*gridInqSize(grid1));

  if ( len != (size_t) (nwpv*gridInqSize(grid2)) )
    cdoAbort("Fields have different gridsize (%s)", __func__);

  if ( nmiss1 > 0 || nmiss2 > 0 )
    {
      for ( i = 0; i < len; i++ )
	if ( !DBL_IS_EQUAL(array2[i], missval2) )
	  {
	    if ( !DBL_IS_EQUAL(array1[i], missval1) )
	      array1[i] += array2[i]*array2[i];
	    else
	      array1[i] = array2[i]*array2[i];
	  }

      field1->nmiss = 0;
      for ( i = 0; i < len; i++ )
	if ( DBL_IS_EQUAL(array1[i], missval1) ) field1->nmiss++;
    }
  else
    {
      for ( i = 0; i < len; i++ ) 
	array1[i] += array2[i]*array2[i];
    }
}


void farsumqw(field_t *field1, field_t field2, double w)
{
  size_t   i, len;
  int          nwpv     = field1->nwpv;
  const int    grid1    = field1->grid;
  const int    nmiss1   = field1->nmiss;
  const double missval1 = field1->missval;
  double *array1  = field1->ptr;
  /*  double *weight1 = field1->weight; */
  const int    grid2    = field2.grid;
  const int    nmiss2   = field2.nmiss;
  const double missval2 = field2.missval;
  double *array2  = field2.ptr;

  if ( nwpv != 2 ) nwpv = 1;

  len = (size_t) (nwpv*gridInqSize(grid1));

  if ( len != (size_t) (nwpv*gridInqSize(grid2)) )
    cdoAbort("Fields have different gridsize (%s)", __func__);

  if ( nmiss1 > 0 || nmiss2 > 0 )
    {
      for ( i = 0; i < len; i++ )
	if ( !DBL_IS_EQUAL(array2[i], missval2) )
	  {
	    if ( !DBL_IS_EQUAL(array1[i], missval1) )
	      array1[i] += w*array2[i]*array2[i];
	    else
	      array1[i] = w*array2[i]*array2[i];
	  }

      field1->nmiss = 0;
      for ( i = 0; i < len; i++ )
	if ( DBL_IS_EQUAL(array1[i], missval1) ) field1->nmiss++;
    }
  else
    {
      for ( i = 0; i < len; i++ ) 
	array1[i] += w*array2[i]*array2[i];
    }
}


void farsub(field_t *field1, field_t field2)
{
  size_t   i, len;
  int          nwpv     = field1->nwpv;
  const int    grid1    = field1->grid;
  const int    nmiss1   = field1->nmiss;
  const double missval1 = field1->missval;
  double *array1  = field1->ptr;
  const int    grid2    = field2.grid;
  const int    nmiss2   = field2.nmiss;
  const double missval2 = field2.missval;
  double *array2  = field2.ptr;

  if ( nwpv != 2 ) nwpv = 1;

  len = (size_t) (nwpv*gridInqSize(grid1));

  if ( len != (size_t) (nwpv*gridInqSize(grid2)) )
    cdoAbort("Fields have different gridsize (%s)", __func__);

  if ( nmiss1 > 0 || nmiss2 > 0 )
    {
      for ( i = 0; i < len; i++ ) 
	array1[i] = SUB(array1[i], array2[i]);

      field1->nmiss = 0;
      for ( i = 0; i < len; i++ )
	if ( DBL_IS_EQUAL(array1[i], missval1) ) field1->nmiss++;
    }
  else
    {
      for ( i = 0; i < len; i++ ) 
	array1[i] -= array2[i];
    }
}


void farmul(field_t *field1, field_t field2)
{
  size_t   i, len;
  int          nwpv     = field1->nwpv;
  const int    grid1    = field1->grid;
  const int    nmiss1   = field1->nmiss;
  const double missval1 = field1->missval;
  double *array1  = field1->ptr;
  const int    grid2    = field2.grid;
  const int    nmiss2   = field2.nmiss;
  const double missval2 = field2.missval;
  double *array2  = field2.ptr;

  if ( nwpv != 2 ) nwpv = 1;

  len = (size_t) (nwpv*gridInqSize(grid1));

  if ( len != (size_t) (nwpv*gridInqSize(grid2)) )
    cdoAbort("Fields have different gridsize (%s)", __func__);

  if ( nmiss1 > 0 || nmiss2 > 0 )
    {
      for ( i = 0; i < len; i++ ) 
	array1[i] = MUL(array1[i], array2[i]);

      field1->nmiss = 0;
      for ( i = 0; i < len; i++ )
	if ( DBL_IS_EQUAL(array1[i], missval1) ) field1->nmiss++;
    }
  else
    {
      for ( i = 0; i < len; i++ ) 
	array1[i] *= array2[i];
    }
}


void fardiv(field_t *field1, field_t field2)
{
  size_t   i, len;
  int          nwpv     = field1->nwpv;
  const int    grid1    = field1->grid;
  const double missval1 = field1->missval;
  double *array1  = field1->ptr;
  const int    grid2    = field2.grid;
  const double missval2 = field2.missval;
  double *array2  = field2.ptr;

  if ( nwpv != 2 ) nwpv = 1;

  len = (size_t) (nwpv*gridInqSize(grid1));

  if ( len != (size_t) (nwpv*gridInqSize(grid2)) )
    cdoAbort("Fields have different gridsize (%s)", __func__);

  for ( i = 0; i < len; i++ ) 
    array1[i] = DIV(array1[i], array2[i]);

  field1->nmiss = 0;
  for ( i = 0; i < len; i++ )
    if ( DBL_IS_EQUAL(array1[i], missval1) ) field1->nmiss++;
}


void faratan2(field_t *field1, field_t field2)
{
  size_t   i, len;
  int          nwpv     = field1->nwpv;
  const int    grid1    = field1->grid;
  const double missval1 = field1->missval;
  double *array1  = field1->ptr;
  const int    grid2    = field2.grid;
  const double missval2 = field2.missval;
  double *array2  = field2.ptr;

  if ( nwpv != 2 ) nwpv = 1;

  len = (size_t) (nwpv*gridInqSize(grid1));

  if ( len != (size_t) (nwpv*gridInqSize(grid2)) )
    cdoAbort("Fields have different gridsize (%s)", __func__);

  for ( i = 0; i < len; i++ ) 
    array1[i] = DBL_IS_EQUAL(array1[i],missval1) || DBL_IS_EQUAL(array2[i],missval2) ? missval1 : atan2(array1[i], array2[i]);

  field1->nmiss = 0;
  for ( i = 0; i < len; i++ )
    if ( DBL_IS_EQUAL(array1[i], missval1) ) field1->nmiss++;
}


void farmin(field_t *field1, field_t field2)
{
  size_t   i, len;
  int          nwpv     = field1->nwpv;
  const int    grid1    = field1->grid;
  const int    nmiss1   = field1->nmiss;
  const double missval1 = field1->missval;
  double *array1  = field1->ptr;
  const int    grid2    = field2.grid;
  const int    nmiss2   = field2.nmiss;
  const double missval2 = field2.missval;
  double *array2  = field2.ptr;

  if ( nwpv != 2 ) nwpv = 1;

  len = (size_t) (nwpv*gridInqSize(grid1));

  if ( len != (size_t) (nwpv*gridInqSize(grid2)) )
    cdoAbort("Fields have different gridsize (%s)", __func__);

  if ( nmiss1 > 0 || nmiss2 > 0 )
    {
      for ( i = 0; i < len; i++ )
	{
	  array1[i] = DBL_IS_EQUAL(array2[i], missval2) ? array1[i] :
	              DBL_IS_EQUAL(array1[i], missval1) ? array2[i] :
		      MIN(array1[i], array2[i]);
	}

      field1->nmiss = 0;
      for ( i = 0; i < len; i++ )
	if ( DBL_IS_EQUAL(array1[i], missval1) ) field1->nmiss++;
    }
  else
    {
      for ( i = 0; i < len; i++ )
	array1[i] = MIN(array1[i], array2[i]);
    }
}


void farmax(field_t *field1, field_t field2)
{
  size_t   i, len;
  int          nwpv     = field1->nwpv;
  const int    grid1    = field1->grid;
  const int    nmiss1   = field1->nmiss;
  const double missval1 = field1->missval;
  double *array1  = field1->ptr;
  const int    grid2    = field2.grid;
  const int    nmiss2   = field2.nmiss;
  const double missval2 = field2.missval;
  double *array2  = field2.ptr;

  if ( nwpv != 2 ) nwpv = 1;

  len = (size_t) (nwpv*gridInqSize(grid1));

  if ( len != (size_t) (nwpv*gridInqSize(grid2)) )
    cdoAbort("Fields have different gridsize (%s)", __func__);

  if ( nmiss1 > 0 || nmiss2 > 0 )
    {
      for ( i = 0; i < len; i++ )
	{
	  array1[i] = DBL_IS_EQUAL(array2[i], missval2) ? array1[i] :
	              DBL_IS_EQUAL(array1[i], missval1) ? array2[i] :
		      MAX(array1[i], array2[i]);
	}

      field1->nmiss = 0;
      for ( i = 0; i < len; i++ )
	if ( DBL_IS_EQUAL(array1[i], missval1) ) field1->nmiss++;
    }
  else
    {
      for ( i = 0; i < len; i++ )
	array1[i] = MAX(array1[i], array2[i]);
    }
}

// not used
void farvar0(field_t *field1, field_t field2, field_t field3)
{
  size_t   i, len;
  int          nwpv     = field1->nwpv;
  const int    grid1    = field1->grid;
  const int    nmiss1   = field1->nmiss;
  const double missval1 = field1->missval;
  double *array1  = field1->ptr;
  const int    grid2    = field2.grid;
  const int    nmiss2   = field2.nmiss;
  const double missval2 = field2.missval;
  double *array2  = field2.ptr;
  const int    nmiss3   = field3.nmiss;
  const double missval3 = field3.missval;
  double *array3  = field3.ptr;

  if ( nwpv != 2 ) nwpv = 1;

  len = (size_t) (nwpv*gridInqSize(grid1));

  if ( len != (size_t) (nwpv*gridInqSize(grid2)) )
    cdoAbort("Fields have different gridsize (%s)", __func__);

  if ( nmiss1 > 0 || nmiss2 > 0 || nmiss3 > 0 )
    {
      for ( i = 0; i < len; i++ )
	{
	  if ( !DBL_IS_EQUAL(array1[i], missval1) && !DBL_IS_EQUAL(array2[i], missval2) && !DBL_IS_EQUAL(array3[i], missval3) )
	    {
	      array1[i] = array2[i]*array3[i] - (array1[i]*array3[i])*(array1[i]*array3[i]);
	      if ( array1[i] < 0 && array1[i] > -1.e-5 ) array1[i] = 0;
	    }
	  else
	    array1[i] = missval1;
	}
    }
  else
    {
      for ( i = 0; i < len; i++ )
	{
	  array1[i] = array2[i]*array3[i] - (array1[i]*array3[i])*(array1[i]*array3[i]);
	  if ( array1[i] < 0 && array1[i] > -1.e-5 ) array1[i] = 0;
	}
    }

  field1->nmiss = 0;
  for ( i = 0; i < len; i++ )
    if ( DBL_IS_EQUAL(array1[i], missval1) || array1[i] < 0 )
      {
	array1[i] = missval1;
	field1->nmiss++;
      }
}


void farvar(field_t *field1, field_t field2, field_t field3, const double divisor)
{
  size_t   i, len;
  int          nwpv     = field1->nwpv;
  const int    grid1    = field1->grid;
  const double missval1 = field1->missval;
  double *array1  = field1->ptr;
  const int    grid2    = field2.grid;
  const double missval2 = field2.missval;
  double *array2  = field2.ptr;
  double *array3  = field3.ptr;
  double temp;

  if ( nwpv != 2 ) nwpv = 1;

  len = (size_t) (nwpv*gridInqSize(grid1));

  if ( len != (size_t) (nwpv*gridInqSize(grid2)) )
    cdoAbort("Fields have different gridsize (%s)", __func__);

  for ( i = 0; i < len; i++ )
    {
      temp      = DIV( MUL(array1[i], array1[i]), array3[i]);
      array1[i] = DIV( SUB(array2[i], temp), array3[i]-divisor);
      if ( array1[i] < 0 && array1[i] > -1.e-5 ) array1[i] = 0;
    }

  field1->nmiss = 0;
  for ( i = 0; i < len; i++ )
    if ( DBL_IS_EQUAL(array1[i], missval1) || array1[i] < 0 )
      {
	array1[i] = missval1;
	field1->nmiss++;
      }
}

// not used
void farstd0(field_t *field1, field_t field2, field_t field3)
{
  size_t   i, len;
  int          nwpv     = field1->nwpv;
  const int    grid1    = field1->grid;
  const double missval1 = field1->missval;
  double *array1  = field1->ptr;
  int    grid2    = field2.grid;

  if ( nwpv != 2 ) nwpv = 1;

  len = (size_t) (nwpv*gridInqSize(grid1));

  if ( len != (size_t) (nwpv*gridInqSize(grid2)) )
    cdoAbort("Fields have different gridsize (%s)", __func__);

  farvar0(field1, field2, field3);

  field1->nmiss = 0;
  for ( i = 0; i < len; i++ )
    if ( DBL_IS_EQUAL(array1[i], missval1) || array1[i] < 0 )
      {
	array1[i] = missval1;
	field1->nmiss++;
      }
    else
      {
	array1[i] = IS_NOT_EQUAL(array1[i], 0) ? sqrt(array1[i]) : 0;
      }
}


void farstd(field_t *field1, field_t field2, field_t field3, const double divisor)
{
  size_t   i, len;
  int          nwpv     = field1->nwpv;
  const int    grid1    = field1->grid;
  const double missval1 = field1->missval;
  double *array1  = field1->ptr;
  const int    grid2    = field2.grid;

  if ( nwpv != 2 ) nwpv = 1;

  len = (size_t) (nwpv*gridInqSize(grid1));

  if ( len != (size_t) (nwpv*gridInqSize(grid2)) )
    cdoAbort("Fields have different gridsize (%s)", __func__);

  farvar(field1, field2, field3, divisor);

  field1->nmiss = 0;
  for ( i = 0; i < len; i++ )
    if ( DBL_IS_EQUAL(array1[i], missval1) || array1[i] < 0 )
      {
	array1[i] = missval1;
	field1->nmiss++;
      }
    else
      {
	array1[i] = IS_NOT_EQUAL(array1[i], 0) ? sqrt(array1[i]) : 0;
      }
}

// not used
void farcvar0(field_t *field1, field_t field2, const double rconst1)
{
  size_t   i, len;
  int          nwpv     = field1->nwpv;
  const int    grid1    = field1->grid;
  const int    nmiss1   = field1->nmiss;
  const double missval1 = field1->missval;
  double *array1  = field1->ptr;
  const int    grid2    = field2.grid;
  const int    nmiss2   = field2.nmiss;
  const double missval2 = field2.missval;
  double *array2  = field2.ptr;
  int    nmiss3   = 0;

  if ( nwpv != 2 ) nwpv = 1;

  len = (size_t) (nwpv*gridInqSize(grid1));

  if ( len != (size_t) (nwpv*gridInqSize(grid2)) )
    cdoAbort("Fields have different gridsize (%s)", __func__);

  if ( DBL_IS_EQUAL(rconst1, missval1) ) nmiss3 = 1;

  if ( nmiss1 > 0 || nmiss2 > 0 || nmiss3 > 0 )
    {
      for ( i = 0; i < len; i++ )
	{
	  if ( !DBL_IS_EQUAL(array1[i], missval1) && !DBL_IS_EQUAL(array2[i], missval2) && nmiss3 == 0 )
	    {
	      array1[i] = array2[i]*rconst1 - (array1[i]*rconst1)*(array1[i]*rconst1);
	      if ( array1[i] < 0 && array1[i] > -1.e-5 ) array1[i] = 0;
	    }
	  else
	    array1[i] = missval1;
	}
    }
  else
    {
      for ( i = 0; i < len; i++ )
	{
	  array1[i] = array2[i]*rconst1 - (array1[i]*rconst1)*(array1[i]*rconst1);
	  if ( array1[i] < 0 && array1[i] > -1.e-5 ) array1[i] = 0;
	}
    }

  field1->nmiss = 0;
  for ( i = 0; i < len; i++ )
    if ( DBL_IS_EQUAL(array1[i], missval1) || array1[i] < 0 )
      {
	array1[i] = missval1;
	field1->nmiss++;
      }
}


void farcvar(field_t *field1, field_t field2, const double rconst1, const double divisor)
{
  size_t   i, len;
  int          nwpv     = field1->nwpv;
  const int    grid1    = field1->grid;
  const double missval1 = field1->missval;
  double *array1  = field1->ptr;
  const int    grid2    = field2.grid;
  const double missval2 = field2.missval;
  double *array2  = field2.ptr;
  double temp;

  if ( nwpv != 2 ) nwpv = 1;

  len = (size_t) (nwpv*gridInqSize(grid1));

  if ( len != (size_t) (nwpv*gridInqSize(grid2)) )
    cdoAbort("Fields have different gridsize (%s)", __func__);

  for ( i = 0; i < len; i++ )
    {
      temp      = DIV( MUL(array1[i], array1[i]), rconst1);
      array1[i] = DIV( SUB(array2[i], temp), rconst1-divisor);
      if ( array1[i] < 0 && array1[i] > -1.e-5 ) array1[i] = 0;
    }

  field1->nmiss = 0;
  for ( i = 0; i < len; i++ )
    if ( DBL_IS_EQUAL(array1[i], missval1) || array1[i] < 0 )
      {
	array1[i] = missval1;
	field1->nmiss++;
      }
}

// not used
void farcstd0(field_t *field1, field_t field2, const double rconst1)
{
  size_t   i, len;
  int          nwpv     = field1->nwpv;
  const int    grid1    = field1->grid;
  const double missval1 = field1->missval;
  double *array1  = field1->ptr;
  const int    grid2    = field2.grid;

  if ( nwpv != 2 ) nwpv = 1;

  len = (size_t) (nwpv*gridInqSize(grid1));

  if ( len != (size_t) (nwpv*gridInqSize(grid2)) )
    cdoAbort("Fields have different gridsize (%s)", __func__);

  farcvar0(field1, field2, rconst1);

  field1->nmiss = 0;
  for ( i = 0; i < len; i++ )
    if ( DBL_IS_EQUAL(array1[i], missval1) || array1[i] < 0 )
      {
	array1[i] = missval1;
	field1->nmiss++;
      }
    else
      {
	array1[i] = IS_NOT_EQUAL(array1[i], 0) ? sqrt(array1[i]) : 0;
      }
}


void farcstd(field_t *field1, field_t field2, const double rconst1, const double divisor)
{
  size_t   i, len;
  int          nwpv     = field1->nwpv;
  const int    grid1    = field1->grid;
  const double missval1 = field1->missval;
  double *array1  = field1->ptr;
  const int    grid2    = field2.grid;

  if ( nwpv != 2 ) nwpv = 1;

  len = (size_t) (nwpv*gridInqSize(grid1));

  if ( len != (size_t) (nwpv*gridInqSize(grid2)) )
    cdoAbort("Fields have different gridsize (%s)", __func__);

  farcvar(field1, field2, rconst1, divisor);

  field1->nmiss = 0;
  for ( i = 0; i < len; i++ )
    if ( DBL_IS_EQUAL(array1[i], missval1) || array1[i] < 0 )
      {
	array1[i] = missval1;
	field1->nmiss++;
      }
    else
      {
	array1[i] = IS_NOT_EQUAL(array1[i], 0) ? sqrt(array1[i]) : 0;
      }
}


void farmoq(field_t *field1, field_t field2)
{
  size_t   i, len;
  int          nwpv     = field1->nwpv;
  const int    grid1    = field1->grid;
  const double missval1 = field1->missval;
  double *array1  = field1->ptr;
  const int    grid2    = field2.grid;
  const int    nmiss2   = field2.nmiss;
  double missval2 = field2.missval;
  double *array2  = field2.ptr;

  if ( nwpv != 2 ) nwpv = 1;

  len = (size_t) (nwpv*gridInqSize(grid1));

  if ( len != (size_t) (nwpv*gridInqSize(grid2)) )
    cdoAbort("Fields have different gridsize (%s)", __func__);

  if ( nmiss2 > 0 )
    {
      for ( i = 0; i < len; i++ )
	if ( !DBL_IS_EQUAL(array2[i], missval2) )
	  array1[i] = array2[i]*array2[i];
	else
	  array1[i] = missval1;

      field1->nmiss = 0;
      for ( i = 0; i < len; i++ )
	if ( DBL_IS_EQUAL(array1[i], missval1) ) field1->nmiss++;
    }
  else
    {
      for ( i = 0; i < len; i++ ) 
	array1[i] = array2[i]*array2[i];
    }
}


void farmoqw(field_t *field1, field_t field2, double w)
{
  size_t   i, len;
  int          nwpv     = field1->nwpv;
  const int    grid1    = field1->grid;
  const double missval1 = field1->missval;
  double *array1  = field1->ptr;
  const int    grid2    = field2.grid;
  const int    nmiss2   = field2.nmiss;
  double missval2 = field2.missval;
  double *array2  = field2.ptr;

  if ( nwpv != 2 ) nwpv = 1;

  len = (size_t) (nwpv*gridInqSize(grid1));

  if ( len != (size_t) (nwpv*gridInqSize(grid2)) )
    cdoAbort("Fields have different gridsize (%s)", __func__);

  if ( nmiss2 > 0 )
    {
      for ( i = 0; i < len; i++ )
	if ( !DBL_IS_EQUAL(array2[i], missval2) )
	  array1[i] = w*array2[i]*array2[i];
	else
	  array1[i] = missval1;

      field1->nmiss = 0;
      for ( i = 0; i < len; i++ )
	if ( DBL_IS_EQUAL(array1[i], missval1) ) field1->nmiss++;
    }
  else
    {
      for ( i = 0; i < len; i++ ) 
	array1[i] = w*array2[i]*array2[i];
    }
}


/* RQ */
/**
 * Counts the number of nonmissing values. The result of the operation
 * is computed according to the following rules:
 * 
 * field1  field2  result
 * a       b       a + 1
 * a       miss    a
 * miss    b       1
 * miss    miss    miss
 * 
 * @param field1 the 1st input field, also holds the result
 * @param field2 the 2nd input field
 */  
void farcount(field_t *field1, field_t field2)
{
  size_t   i, len;
  int          nwpv     = field1->nwpv;
  const int    grid1    = field1->grid;
  const int    nmiss1   = field1->nmiss;
  const double missval1 = field1->missval;
  double *array1  = field1->ptr;
  /*  double *weight1 = field1->weight; */
  const int    grid2    = field2.grid;
  const int    nmiss2   = field2.nmiss;
  const double missval2 = field2.missval;
  double *array2  = field2.ptr;

  if ( nwpv != 2 ) nwpv = 1;

  len = (size_t) (nwpv*gridInqSize(grid1));

  if ( len != (size_t) (nwpv*gridInqSize(grid2)) )
    cdoAbort("Fields have different gridsize (%s)", __func__);

  if ( nmiss1 > 0 || nmiss2 > 0 )
    {
      for ( i = 0; i < len; i++ )
	if ( !DBL_IS_EQUAL(array2[i], missval2) )
	  {
	    if ( !DBL_IS_EQUAL(array1[i], missval1) )
	      array1[i] += 1.0;
	    else
	      array1[i] = 1.0;
	  }

      field1->nmiss = 0;
      for ( i = 0; i < len; i++ )
	if ( DBL_IS_EQUAL(array1[i], missval1) ) field1->nmiss++;
    }
  else
    {
      for ( i = 0; i < len; i++ ) 
	array1[i] += 1.0;
    }
}
/* QR */
