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
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include "nth_element.h"

#define SWAP(a, b) {t = (a); (a) = (b); (b) = t;}

double nth_element(double array[], int length, int n)
{
  int i, j, k = 0, l = length-1, m;
  double a,t;

  assert( array != NULL );
  assert( l >= 0 );
  assert( n >= 0 );
  assert( n <= l );

  for ( ; ; )
    {
      if ( l <= k+1 )
        { 
          if ( l == k+1 && array[l] < array[k] ) 
            SWAP(array[k], array[l])
          return array[n];
        } 
      else 
        {
          m = (k+l) >> 1; 
          SWAP(array[m],array[k+1]);
          
          if ( array[k] > array[l] ) 
            SWAP(array[k], array[l])
          if ( array[k+1] > array[l] ) 
            SWAP(array[k+1], array[l])
          if ( array[k] > array[k+1] ) 
            SWAP(array[k], array[k+1])

          i = k+1; 
          j = l;
          a = array[k+1]; 

          for ( ; ; ) 
            { 
              do i++; while ( array[i] < a ); 
              do j--; while ( array[j] > a ); 
              if ( j < i ) break; 
              SWAP(array[i], array[j])
            } 

          array[k+1] = array[j]; 
          array[j] = a;
          if ( j >= n ) l = j-1; 
          if ( j <= n ) k = i; 
        }
    }
}
