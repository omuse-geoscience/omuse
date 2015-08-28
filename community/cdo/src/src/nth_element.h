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

#ifndef NTH_ELEMENT_H_
#define NTH_ELEMENT_H_

/**
 * Finds the nth smallest value of an array of values. The array elements
 * are rearranged as far as necessary to get the nth element to its proper
 * position, with no element less than the nth element placed after it.
 * Based on the Quicksort algorithm.
 * 
 * @param array   the array
 * @param length  the length of the array
 * @param n       the ordinal number of the element to be found. Must be
 *                in the interval [0, length - 1]
 * 
 * @return the nth smallest value in the input array
 * 
 * @author Ralf Quast
 * @version 1.0
 */
double nth_element(double array[], int length, int n);

#endif /*NTH_ELEMENT_H_*/
