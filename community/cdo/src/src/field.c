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
/* RQ */
#include "nth_element.h"
#include "merge_sort2.h"
/* QR */

double crps_det_integrate(double *a, const double d, const size_t n);

double _FADD_(const double x, const double y, const double missval1, const double missval2)
{
  return FADD(x,y);
}

double _FSUB_(const double x, const double y, const double missval1, const double missval2)
{
  return FSUB(x, y);
}

double _FMUL_(const double x, const double y, const double missval1, const double missval2)
{
  return FMUL(x, y);
}

double _FDIV_(const double x, const double y, const double missval1, const double missval2)
{
  return FDIV(x, y);
}

double _FPOW_(const double x, const double y, const double missval1, const double missval2)
{
  return FPOW(x, y);
}

double _FSQRT_(const double x, const double missval1)
{
  return FSQRT(x);
}


double fldfun(field_t field, const int function)
{
  double rval = 0;

  // Why is this not a switch-case statement?
  if      ( function == func_min  )  rval = fldmin(field);
  else if ( function == func_max  )  rval = fldmax(field);
  else if ( function == func_sum  )  rval = fldsum(field);
  else if ( function == func_mean )  rval = fldmean(field);
  else if ( function == func_avg  )  rval = fldavg(field);
  else if ( function == func_std  )  rval = fldstd(field);
  else if ( function == func_std1 )  rval = fldstd1(field);
  else if ( function == func_var  )  rval = fldvar(field);
  else if ( function == func_var1 )  rval = fldvar1(field);
  
  else if ( function == func_crps )  rval = fldcrps(field);
  else if ( function == func_brs )   rval = fldbrs(field);

  else if ( function == func_rank )  rval = fldrank(field);
  else if ( function == func_roc )   rval = fldroc(field);
  else cdoAbort("%s: function %d not implemented!", __func__, function);

  return rval;
}


double fldrank(field_t field) 
{
  double res = 0;
  // Using first value as reference (observation)
  double *array  =  &(field.ptr[1]);
  double val     = array[-1];
  const double missval = field.missval;
  int nmiss      = field.nmiss;
  const size_t len       = field.size-1;
  size_t j;
  
  if ( nmiss )  return(missval);

  sort_iter_single(len,array, 1);

  if ( val > array[len-1] ) 
    res=(double)len;
  else 
    for ( j=0; j<len; j++ )
      if ( array[j] >= val ) {
	res=(double)j; 
	break;
      }

  return res;
}


double fldroc(field_t field) 
{
  return field.missval;
}

double fldcrps(field_t field)
{
  const size_t len     = field.size;
  const int    nmiss   = field.nmiss;
  double *array  = field.ptr;

  if ( nmiss > 0 ) 
    cdoAbort("Missing values not implemented in crps calculation");
  // possible handling of missing values:
  // (1) strip them off, and sort array without missing values
  //     using only (len - 1 - nmiss) values
  // (2) modify merge_sort in a way, that missing values will
  //     always go to the end of the list

  // Use first value as reference
  sort_iter_single(len-1,&array[1],ompNumThreads);

  return crps_det_integrate(&array[1],array[0],len-1);
}


double fldbrs(field_t field) 
{
  const size_t    len   = field.size;
  const int     nmiss   = field.nmiss;
  double *array   = field.ptr;
  const double missval  = field.missval;

  double brs = 0;
  size_t i, count=0;

  // Using first value as reference
  if ( nmiss == 0 ) 
    {
      for ( i=1; i<len; i++ )
	brs += (array[i] - array[0]) * (array[i] - array[0]);
      count = i-1;
    }
  else 
    {
      if ( DBL_IS_EQUAL(array[0], missval) ) return missval;

      for ( i=1; i<len; i++ )
	if ( !DBL_IS_EQUAL(array[i], missval) ) 
	  {
	    brs += (array[i] - array[0]) * (array[i] - array[0]);
	    count ++;
	  }
    }

  return brs/count;
}


double fldmin(field_t field)
{
  size_t   i;
  const size_t len     = field.size;
  const int    nmiss   = field.nmiss;
  const double missval = field.missval;
  const double *restrict array  = field.ptr;
  double rmin = 0;

  if ( nmiss > 0 )
    {
      rmin = DBL_MAX;
      for ( i = 0; i < len; i++ ) 
	if ( !DBL_IS_EQUAL(array[i], missval) )
	  if ( array[i] < rmin ) rmin = array[i];

      if ( IS_EQUAL(rmin, DBL_MAX) )
	rmin = missval;
    }
  else
    {
      rmin = array[0];
      //#pragma simd reduction(min:rmin) 
      for ( i = 1; i < len; i++ ) 
	if ( array[i] < rmin )  rmin = array[i];
    }

  return (rmin);
}


double fldmax(field_t field)
{
  size_t   i;
  const size_t len     = field.size;
  const int    nmiss   = field.nmiss;
  const double missval = field.missval;
  double *array  = field.ptr;
  double rmax = 0;

  if ( nmiss > 0 )
    {
      rmax = -DBL_MAX;
      for ( i = 0; i < len; i++ )
        if ( !DBL_IS_EQUAL(array[i], missval) )
          if ( array[i] > rmax ) rmax = array[i];
      
      if ( IS_EQUAL(rmax, -DBL_MAX) )
        rmax = missval;
    }
  else
    {
      rmax = array[0];
      for ( i = 1; i < len; i++ ) 
        if ( array[i] > rmax )  rmax = array[i];
    }

  return (rmax);
}


double fldsum(field_t field)
{
  size_t   i;
  size_t   nvals   = 0;
  const size_t len     = field.size;
  const int    nmiss   = field.nmiss;
  const double missval = field.missval;
  double *array  = field.ptr;
  double rsum = 0;

  if ( nmiss > 0 )
    {
      for ( i = 0; i < len; i++ ) 
	if ( !DBL_IS_EQUAL(array[i], missval) )
	  {
	    rsum += array[i];
	    nvals++;
	  }

      if ( !nvals ) rsum = missval;
    }
  else
    {
      for ( i = 0; i < len; i++ ) 
	rsum += array[i];
    }

  return (rsum);
}


double fldmean(field_t field)
{
  size_t   i;
  const size_t len      = field.size;
  const int    nmiss    = field.nmiss;
  const double missval1 = field.missval;
  const double missval2 = field.missval;
  double *array   = field.ptr;
  double *w       = field.weight;
  double rsum = 0, rsumw = 0, ravg = 0;

  if ( nmiss > 0 )
    {
      for ( i = 0; i < len; i++ ) 
	if ( !DBL_IS_EQUAL(array[i], missval1) && !DBL_IS_EQUAL(w[i], missval1) )
	  {
	    rsum  += w[i] * array[i];
	    rsumw += w[i];
	  }
    }
  else
    {
      for ( i = 0; i < len; i++ ) 
	{
	  rsum  += w[i] * array[i];
	  rsumw += w[i];
	}
    }

  ravg = DIV(rsum, rsumw);

  return (ravg);
}


double fldavg(field_t field)
{
  size_t   i;
  const size_t len      = field.size;
  const int    nmiss    = field.nmiss;
  const double missval1 = field.missval;
  const double missval2 = field.missval;
  double *array   = field.ptr;
  double *w       = field.weight;
  double rsum = 0, rsumw = 0, ravg = 0;

  if ( nmiss > 0 )
    {
      for ( i = 0; i < len; i++ ) 
	if ( !DBL_IS_EQUAL(w[i], missval1) )
	  {
	    rsum  = ADD(rsum, MUL(w[i], array[i]));
	    rsumw = ADD(rsumw, w[i]);
	  }
    }
  else
    {
      for ( i = 0; i < len; i++ ) 
	{
	  rsum  += w[i] * array[i];
	  rsumw += w[i];
	}
    }

  ravg = DIV(rsum, rsumw);

  return (ravg);
}

static
void prevarsum(const double *restrict array, const double *restrict w, size_t len, int nmiss, 
	       double missval, double *restrict rsum, double *restrict rsumw, double *restrict rsumq, double *restrict rsumwq)
{ 
  size_t i;
  *rsum = *rsumw = 0;
  *rsumq = *rsumwq = 0;

  if ( nmiss > 0 )
    {
      for ( i = 0; i < len; i++ ) 
        if ( !DBL_IS_EQUAL(array[i], missval) && !DBL_IS_EQUAL(w[i], missval) )
          {
            *rsum   += w[i] * array[i];
            *rsumq  += w[i] * array[i] * array[i];
            *rsumw  += w[i];
            *rsumwq += w[i] * w[i];
          }
    }
  else
    {
      for ( i = 0; i < len; i++ ) 
        {
          *rsum   += w[i] * array[i];
          *rsumq  += w[i] * array[i] * array[i];
          *rsumw  += w[i];
          *rsumwq += w[i] * w[i];
        }
    }
}


double fldvar(field_t field)
{
  const size_t len     = field.size;
  const int    nmiss   = field.nmiss;
  const double missval = field.missval;
  double *array  = field.ptr;
  double *w      = field.weight;
  double rsum, rsumw;
  double rsumq, rsumwq;

  prevarsum(array, w, len, nmiss, missval, &rsum, &rsumw, &rsumq, &rsumwq);

  double rvar = IS_NOT_EQUAL(rsumw, 0) ? (rsumq*rsumw - rsum*rsum) / (rsumw*rsumw) : missval;
  if ( rvar < 0 && rvar > -1.e-5 ) rvar = 0;

  return (rvar);
}


double fldvar1(field_t field)
{
  const size_t len     = field.size;
  const int    nmiss   = field.nmiss;
  const double missval = field.missval;
  double *array  = field.ptr;
  double *w      = field.weight;
  double rsum, rsumw;
  double rsumq, rsumwq;

  prevarsum(array, w, len, nmiss, missval, &rsum, &rsumw, &rsumq, &rsumwq);

  double rvar = (rsumw*rsumw > rsumwq) ? (rsumq*rsumw - rsum*rsum) / (rsumw*rsumw - rsumwq) : missval;
  if ( rvar < 0 && rvar > -1.e-5 ) rvar = 0;

  return (rvar);
}


double var_to_std(double rvar, double missval)
{
  double rstd;

  if ( DBL_IS_EQUAL(rvar, missval) || rvar < 0 )
    {
      rstd = missval;
    }
  else
    {
      rstd = IS_NOT_EQUAL(rvar, 0) ? sqrt(rvar) : 0;
    }

  return rstd;
}

double fldstd(field_t field)
{
  return var_to_std(fldvar(field), field.missval);
}


double fldstd1(field_t field)
{
  return var_to_std(fldvar1(field), field.missval);
}


void fldrms(field_t field, field_t field2, field_t *field3)
{
  size_t   i;
  size_t len;
  int    rnmiss = 0;
  int    grid1    = field.grid;
  //  int    nmiss1   = field.nmiss;
  double *array1  = field.ptr;
  int    grid2    = field2.grid;
  //  int    nmiss2   = field2.nmiss;
  double *array2  = field2.ptr;
  const double missval1 = field.missval;
  const double missval2 = field2.missval;
  double *w       = field.weight;
  double rsum = 0, rsumw = 0, ravg = 0;

  len    = gridInqSize(grid1);
  if ( len != (size_t) gridInqSize(grid2) )
    cdoAbort("fields have different size!");

  /*
  if ( nmiss1 > 0 )
  */
    {
      for ( i = 0; i < len; i++ ) 
	if ( !DBL_IS_EQUAL(w[i], missval1) )
	  {
	    rsum  = ADD(rsum, MUL(w[i], MUL(SUB(array2[i], array1[i]),
                                            SUB(array2[i], array1[i]))));
	    rsumw = ADD(rsumw, w[i]);
	  }
    }
    /*
  else
    {
      for ( i = 0; i < len; i++ ) 
	{
	  rsum  += w[i] * array1[i];
	  rsumw += w[i];
	}
    }
    */

  ravg = SQRT(DIV(rsum, rsumw));

  if ( DBL_IS_EQUAL(ravg, missval1) ) rnmiss++;

  field3->ptr[0] = ravg;
  field3->nmiss  = rnmiss;
}


void varrms(field_t field, field_t field2, field_t *field3)
{
  size_t   i, k, nlev, len;
  int    rnmiss = 0;
  int    zaxis    = field.zaxis;
  int    grid1    = field.grid;
  //  int    nmiss1   = field.nmiss;
  double *array1  = field.ptr;
  int    grid2    = field2.grid;
  //  int    nmiss2   = field2.nmiss;
  double *array2  = field2.ptr;
  const double missval1 = field.missval;
  const double missval2 = field2.missval;
  double *w       = field.weight;
  double rsum = 0, rsumw = 0, ravg = 0;

  nlev   = zaxisInqSize(zaxis);
  len    = gridInqSize(grid1);
  if ( len != (size_t) gridInqSize(grid2) )
    cdoAbort("fields have different size!");

  /*
  if ( nmiss1 > 0 )
  */
    {
      for ( k = 0; k < nlev; k++ )
	for ( i = 0; i < len; i++ )
	  /*	  if ( !DBL_IS_EQUAL(w[i], missval1) ) */
	    {
	      rsum  = ADD(rsum, MUL(w[i], MUL(SUB(array2[k*len+i], array1[k*len+i]),
                                              SUB(array2[k*len+i], array1[k*len+i]))));
	      rsumw = ADD(rsumw, w[i]);
	    }
    }
    /*
  else
    {
      for ( i = 0; i < len; i++ )
	{
	  rsum  += w[i] * array1[i];
	  rsumw += w[i];
	}
    }
    */

  ravg = SQRT(DIV(rsum, rsumw));

  if ( DBL_IS_EQUAL(ravg, missval1) ) rnmiss++;

  field3->ptr[0] = ravg;
  field3->nmiss  = rnmiss;
}

/* RQ */
double fldpctl(field_t field, const int p)
{
  const size_t len     = field.size;
  const int    nmiss   = field.nmiss;
  const double missval = field.missval;
  double *array  = field.ptr;
  double *array2;

  size_t i, j;
  double pctl = missval;

  if ( len - nmiss > 0 )
    {
      if ( nmiss > 0 )
        {
          array2 = (double*) malloc((len - nmiss)*sizeof(double));

          for ( i = 0, j = 0; i < len; i++ ) 
            if ( !DBL_IS_EQUAL(array[i], missval) )
              array2[j++] = array[i];

          pctl = nth_element(array2, j, (int)ceil(j*(p/100.0))-1);

          free(array2);
        }
      else
        {
          pctl = nth_element(array, len, (int)ceil(len*(p/100.0))-1);
        }
    }

  return pctl;
}
/* QR */

/*  field_t UTILITIES */
/*  update the number non missing values */
void fldunm(field_t *field)
{
  size_t i;

  field->nmiss = 0;
  for ( i = 0; i < field->size; i++ )
    if ( DBL_IS_EQUAL(field->ptr[i], field->missval) ) field->nmiss++;
}

/*  check for non missval values */
int fldhvs(field_t *fieldPtr, const size_t nlevels)
{
  size_t level;
  field_t field;

  for ( level = 0; level < nlevels; level++)
    {
      field = fieldPtr[level];
      if ( (size_t)field.nmiss != field.size )
        return TRUE;
    }

  return FALSE;
}


double crps_det_integrate(double *a, const double d, const size_t n)
{
  /* *************************************************************************** */
  /* This routine finds the area between the cdf described by the ordered array  */
  /* of doubles (double *a) and the Heavyside function H(d)                      */
  /* INPUT ARGUMENTS:                                                            */
  /*     double *a  - ordered array of doubles describing a cdf                  */
  /*                  as cdf(a[i]) = ( (double)i )/ n                            */
  /*     double d   - describing a reference value                               */
  /*     int n      - the length of array a                                      */
  /* RETURN VALUE:                                                               */
  /*     double     - area under the curve in units of a                         */
  /* *************************************************************************** */

  double area = 0; 
  //  double tmp;
  size_t i;
#if defined(_OPENMP)
#pragma omp parallel for if ( n>10000 ) shared(a) private(i) \
  reduction(+:area) schedule(static,10000) 
#endif                                                     /* **************************** */
  for ( i=1; i<n; i++ ) {                                  /* INTEGRATE CURVE AREA         */
    if ( a[i] < d )                                        /* left of heavyside            */
      area += (a[i]-a[i-1])*(double)i*i/n/n;                   /*                              */
    else if ( a[i-1] > d )                                 /* right of heavyside           */
      area += (a[i]-a[i-1])*(1.-(double)i/n)*(1.-(double)i/n);              /*                              */
    else if ( a[i-1] < d && a[i] > d ) {                   /* hitting jump pf heavyside    */
      area += (d-a[i-1]) * (double)i*i/n/n;                    /* (occurs exactly once!)       */
      area += (a[i]-d) * (1.-(double)i/n)*(1.-(double)i/n);                 /* **************************** */
    }
  }


  return(area);
}

