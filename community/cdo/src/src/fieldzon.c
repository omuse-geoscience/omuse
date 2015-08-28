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
#include "nth_element.h"


void zonfun(field_t field1, field_t *field2, int function)
{
  if      ( function == func_min   )  zonmin(field1, field2);
  else if ( function == func_max   )  zonmax(field1, field2);  
  else if ( function == func_range )  zonrange(field1, field2);  
  else if ( function == func_sum   )  zonsum(field1, field2);  
  else if ( function == func_mean  )  zonmean(field1, field2);  
  else if ( function == func_avg   )  zonavg(field1, field2);  
  else if ( function == func_std   )  zonstd(field1, field2);  
  else if ( function == func_std1  )  zonstd1(field1, field2);  
  else if ( function == func_var   )  zonvar(field1, field2);
  else if ( function == func_var1  )  zonvar1(field1, field2);
  else cdoAbort("function %d not implemented!", function);
}


void zonmin(field_t field1, field_t *field2)
{
  long   i, j, nx, ny;
  int    rnmiss = 0;
  int    grid    = field1.grid;
  int    nmiss   = field1.nmiss;
  double missval = field1.missval;
  double *array  = field1.ptr;
  double rmin = 0;

  nx    = gridInqXsize(grid);
  ny    = gridInqYsize(grid);

  for ( j = 0; j < ny; j++ )
    {
      if ( nmiss > 0 )
	{
	  rmin = DBL_MAX;
	  for ( i = 0; i < nx; i++ )
	    if ( !DBL_IS_EQUAL(array[j*nx+i], missval) )
	      if ( array[j*nx+i] < rmin ) rmin = array[j*nx+i];

	  if ( IS_EQUAL(rmin, DBL_MAX) )
	    {
	      rnmiss++;
	      rmin = missval;
	    }
	}
      else
	{
	  rmin = array[j*nx];
	  for ( i = 1; i < nx; i++ )
	    if ( array[j*nx+i] < rmin )  rmin = array[j*nx+i];
	}

      field2->ptr[j] = rmin;
    }

  field2->nmiss  = rnmiss;
}


void zonmax(field_t field1, field_t *field2)
{
  long   i, j, nx, ny;
  int    rnmiss = 0;
  int    grid    = field1.grid;
  int    nmiss   = field1.nmiss;
  double missval = field1.missval;
  double *array  = field1.ptr;
  double rmax = 0;

  nx    = gridInqXsize(grid);
  ny    = gridInqYsize(grid);

  for ( j = 0; j < ny; j++ )
    {
      if ( nmiss > 0 )
	{
	  rmax = -DBL_MAX;
	  for ( i = 0; i < nx; i++ )
	    if ( !DBL_IS_EQUAL(array[j*nx+i], missval) )
	      if ( array[j*nx+i] > rmax ) rmax = array[j*nx+i];

	  if ( IS_EQUAL(rmax, -DBL_MAX) )
	    {
	      rnmiss++;
	      rmax = missval;
	    }
	}
      else
	{
	  rmax = array[j*nx];
	  for ( i = 1; i < nx; i++ ) 
	    if ( array[j*nx+i] > rmax )  rmax = array[j*nx+i];
	}

      field2->ptr[j] = rmax;
    }

  field2->nmiss  = rnmiss;
}


void zonrange(field_t field1, field_t *field2)
{
  long   i, j, nx, ny;
  int    rnmiss = 0;
  int    grid    = field1.grid;
  int    nmiss   = field1.nmiss;
  double missval = field1.missval;
  double *array  = field1.ptr;
  double rmin = 0;
  double rmax = 0;
  double rrange = 0;

  nx    = gridInqXsize(grid);
  ny    = gridInqYsize(grid);

  for ( j = 0; j < ny; j++ )
    {
      if ( nmiss > 0 )
	{
	  rmin =  DBL_MAX;
	  rmax = -DBL_MAX;
	  for ( i = 0; i < nx; i++ )
	    if ( !DBL_IS_EQUAL(array[j*nx+i], missval) )
	      {
		if      ( array[j*nx+i] < rmin ) rmin = array[j*nx+i];
		else if ( array[j*nx+i] > rmax ) rmax = array[j*nx+i];
	      }

	  if ( IS_EQUAL(rmin, DBL_MAX) || IS_EQUAL(rmax, -DBL_MAX) )
	    {
	      rnmiss++;
	      rrange = missval;
	    }
	  else
	    {
	      rrange = rmax - rmin;
	    }
	}
      else
	{
	  rmin = array[j*nx];
	  rmax = array[j*nx];
	  for ( i = 1; i < nx; i++ )
	    {
	      if      ( array[j*nx+i] < rmin )  rmin = array[j*nx+i];
	      else if ( array[j*nx+i] > rmax )  rmax = array[j*nx+i];
	    }

	  rrange = rmax - rmin;
	}

      field2->ptr[j] = rrange;
    }

  field2->nmiss  = rnmiss;
}


void zonsum(field_t field1, field_t *field2)
{
  long   i, j, nx, ny;
  long   nvals   = 0;
  int    rnmiss  = 0;
  int    grid    = field1.grid;
  int    nmiss   = field1.nmiss;
  double missval = field1.missval;
  double *array  = field1.ptr;
  double rsum = 0;

  nx    = gridInqXsize(grid);
  ny    = gridInqYsize(grid);

  for ( j = 0; j < ny; j++ )
    {
      if ( nmiss > 0 )
	{
	  nvals = 0;
	  rsum = 0;
	  for ( i = 0; i < nx; i++ )
	    if ( !DBL_IS_EQUAL(array[j*nx+i], missval) )
	      {
		rsum += array[j*nx+i];
		nvals++;
	      }

	  if ( !nvals )
	    {
	      rsum = missval;
	      rnmiss++;
	    }
 	}
      else
	{
	  rsum = 0;
	  for ( i = 0; i < nx; i++ )
	    rsum += array[j*nx+i];
	}

      field2->ptr[j] = rsum;
    }

  field2->nmiss  = rnmiss;
}


void zonmean(field_t field1, field_t *field2)
{
  long   i, j, nx, ny;
  int    rnmiss = 0;
  int    grid     = field1.grid;
  int    nmiss    = field1.nmiss;
  double missval1 = field1.missval;
  double missval2 = field1.missval;
  double *array   = field1.ptr;
  double rsum = 0, rsumw = 0, ravg = 0;

  nx    = gridInqXsize(grid);
  ny    = gridInqYsize(grid);

  for ( j = 0; j < ny; j++ )
    {
      rsum  = 0;
      rsumw = 0;
      if ( nmiss > 0 )
	{
	  for ( i = 0; i < nx; i++ )
	    if ( !DBL_IS_EQUAL(array[j*nx+i], missval1) )
	      {
		rsum  += array[j*nx+i];
		rsumw += 1;
	      }
	}
      else
	{
	  for ( i = 0; i < nx; i++ )
	    {
	      rsum  += array[j*nx+i];
	      rsumw += 1;
	    }
	}

      ravg = DIV(rsum, rsumw);

      if ( DBL_IS_EQUAL(ravg, missval1) ) rnmiss++;

      field2->ptr[j] = ravg;
    }

  field2->nmiss  = rnmiss;
}


void zonavg(field_t field1, field_t *field2)
{
  long   i, j, nx, ny;
  int    rnmiss = 0;
  int    grid     = field1.grid;
  int    nmiss    = field1.nmiss;
  double missval1 = field1.missval;
  double missval2 = field1.missval;
  double *array   = field1.ptr;
  double rsum = 0, rsumw = 0, ravg = 0;

  nx    = gridInqXsize(grid);
  ny    = gridInqYsize(grid);

  for ( j = 0; j < ny; j++ )
    {
      rsum  = 0;
      rsumw = 0;
      if ( nmiss > 0 )
	{
	  for ( i = 0; i < nx; i++ )
	    {
	      rsum   = ADD(rsum, array[j*nx+i]);
	      rsumw += 1;
	    }
	}
      else
	{
	  for ( i = 0; i < nx; i++ )
	    {
	      rsum  += array[j*nx+i];
	      rsumw += 1;
	    }
	}

      ravg = DIV(rsum, rsumw);

      if ( DBL_IS_EQUAL(ravg, missval1) ) rnmiss++;

      field2->ptr[j] = ravg;
    }

  field2->nmiss  = rnmiss;
}

static
void prevarsum_zon(const double *restrict array, int nx, int nmiss,  double missval, 
                   double *rsum, double *rsumw, double *rsumq, double *rsumwq)
{
  double w = 1./nx;

  *rsum   = 0;
  *rsumq  = 0;
  *rsumw  = 0;
  *rsumwq = 0;

  if ( nmiss > 0 )
    {
      for ( int i = 0; i < nx; i++ )
        if ( !DBL_IS_EQUAL(array[i], missval) )
          {
            *rsum   += w * array[i];
            *rsumq  += w * array[i] * array[i];
            *rsumw  += w;
            *rsumwq += w * w;
          }
    }
  else
    {
      for ( int i = 0; i < nx; i++ )
        {
          *rsum   += w * array[i];
          *rsumq  += w * array[i] * array[i];
          *rsumw  += w;
          *rsumwq += w * w;
        }
    }
}


void zonvar(field_t field1, field_t *field2)
{
  int    rnmiss = 0;
  int    grid     = field1.grid;
  int    nmiss    = field1.nmiss;
  double missval1 = field1.missval;
  double *array   = field1.ptr;
  double rsum = 0, rsumw = 0, rvar = 0;
  double rsumq = 0, rsumwq = 0;

  int nx = gridInqXsize(grid);
  int ny = gridInqYsize(grid);

  for ( int j = 0; j < ny; j++ )
    {
      prevarsum_zon(array+j*nx, nx, nmiss, missval1, &rsum, &rsumw, &rsumq, &rsumwq);

      rvar = IS_NOT_EQUAL(rsumw, 0) ? (rsumq*rsumw - rsum*rsum) / (rsumw*rsumw) : missval1;
      if ( rvar < 0 && rvar > -1.e-5 ) rvar = 0;

      if ( DBL_IS_EQUAL(rvar, missval1) ) rnmiss++;

      field2->ptr[j] = rvar;
    }

  field2->nmiss  = rnmiss;
}


void zonvar1(field_t field1, field_t *field2)
{
  int    rnmiss = 0;
  int    grid     = field1.grid;
  int    nmiss    = field1.nmiss;
  double missval1 = field1.missval;
  double *array   = field1.ptr;
  double rsum = 0, rsumw = 0, rvar = 0;
  double rsumq = 0, rsumwq = 0;

  int nx = gridInqXsize(grid);
  int ny = gridInqYsize(grid);

  for ( int j = 0; j < ny; j++ )
    {
      prevarsum_zon(array+j*nx, nx, nmiss, missval1, &rsum, &rsumw, &rsumq, &rsumwq);

      rvar = (rsumw*rsumw > rsumwq) ? (rsumq*rsumw - rsum*rsum) / (rsumw*rsumw - rsumwq) : missval1;
      if ( rvar < 0 && rvar > -1.e-5 ) rvar = 0;

      if ( DBL_IS_EQUAL(rvar, missval1) ) rnmiss++;

      field2->ptr[j] = rvar;
    }

  field2->nmiss  = rnmiss;
}


void zonstd(field_t field1, field_t *field2)
{
  int    rnmiss = 0;
  int    grid    = field1.grid;
  double missval = field1.missval;
  double rstd;

  int ny = gridInqYsize(grid);

  zonvar(field1, field2);

  for ( int j = 0; j < ny; j++ )
    {
      rstd = var_to_std(field2->ptr[j], missval);

      if ( DBL_IS_EQUAL(rstd, missval) ) rnmiss++;

      field2->ptr[j] = rstd;
    }

  field2->nmiss = rnmiss;
}


void zonstd1(field_t field1, field_t *field2)
{
  int    rnmiss = 0;
  int    grid    = field1.grid;
  double missval = field1.missval;
  double rstd;

  int ny = gridInqYsize(grid);

  zonvar1(field1, field2);

  for ( int j = 0; j < ny; j++ )
    {
      rstd = var_to_std(field2->ptr[j], missval);

      if ( DBL_IS_EQUAL(rstd, missval) ) rnmiss++;

      field2->ptr[j] = rstd;
    }

  field2->nmiss = rnmiss;
}

/* RQ */
void zonpctl(field_t field1, field_t *field2, int p)
{
  long   i, j, l, nx, ny;
  int    rnmiss = 0;
  int    grid    = field1.grid;
  int    nmiss   = field1.nmiss;
  double missval = field1.missval;
  double *array  = field1.ptr;
  double *array2;

  nx = gridInqXsize(grid);
  ny = gridInqYsize(grid);
  
  if ( nmiss > 0 )
    {
      array2 = (double*) malloc(nx*sizeof(double));
      
      for ( j = 0; j < ny; j++ )
        {
          for ( i = 0, l = 0; i < nx; i++ )
            if ( !DBL_IS_EQUAL(array[j*nx+i], missval) )
              array2[l++] = array[j*nx+i];
	    
          if ( l > 0 )
            {
              field2->ptr[j] = nth_element(array2, l, (int)ceil(l*(p/100.0))-1);
            }
          else
            {
              field2->ptr[j] = missval;
              rnmiss++;
            }
        }
        
      free(array2);
    }
  else
    {
      for ( j = 0; j < ny; j++ )
        {
          if ( nx > 0 )
            {
              field2->ptr[j] = nth_element(&array[j*nx], nx, (int)ceil(nx*(p/100.0))-1);
            }
          else
            {
              field2->ptr[j] = missval;
              rnmiss++;
            }
        }
    }

  field2->nmiss = rnmiss;
}
/* QR */
