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

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <float.h>
#include "cdo_int.h"
#include "dmemory.h"
#include "list.h"


#define  DEFAULT_ALLINC  1024


static void listInit(LIST *list, int type)
{
  list->array  = NULL;
  list->nalloc = 0;
  list->allinc = DEFAULT_ALLINC;
  list->type   = type;
}


LIST *listNew(int type)
{
  LIST *list = NULL;

  if ( type != INT_LIST && type != FLT_LIST )
    {
      fprintf(stderr, "%s: type %d unsupported!\n", __func__, type);
    }
  else
    {
      list = (LIST*) malloc(sizeof(LIST));
      listInit(list, type);
    }

  return (list);
}


void listDelete(LIST *list)
{
  if ( list )
    {
      if ( list->array ) free(list->array);
      free(list);
    }
}


void *listArrayPtr(LIST *list)
{
  return (list->array);
}


static void listCheck(LIST *list, int num)
{
  while ( list->nalloc < (num+1) )
    {
      list->nalloc += list->allinc;
      if ( list->type == INT_LIST )
	list->array = (int*) realloc(list->array, list->nalloc*sizeof(int));
      else
	list->array = (double*) realloc(list->array, list->nalloc*sizeof(double));
    }
}


void listSetInt(LIST *list, int num, int ival)
{
  listCheck(list, num);

  ((int *) list->array)[num] = ival;
}


void listSetFlt(LIST *list, int num, double fval)
{
  listCheck(list, num);

  ((double *) list->array)[num] = fval;
}


int listGetInt(LIST *list, int num)
{
  int ival;

  ival = ((int *) list->array)[num];

  return (ival);
}


double listGetFlt(LIST *list, int num)
{
  double fval;

  fval = ((double *) list->array)[num];

  return (fval);
}

static
int get_ival(const char *intstr, int idefault, int istart, int iend, int *ilast)
{
  int ival = idefault;
  int i;

  for ( i = istart; i < iend; i++ )
    {
      if ( ! (isdigit(intstr[i]) || intstr[i] == '-') )
	{
	  if ( intstr[i] == '/' )
	    ival = parameter2intlist(intstr+i+1);
	  else
	    fprintf(stderr, "Syntax error in >%.*s<! Character %c not allowed.\n", iend, intstr, intstr[i]);
	  break;
	}
    }

  *ilast = i;

  return (ival);
}


void split_intstring(const char *intstr, int *first, int *last, int *inc)
{
  int i, start;
  int istrlen;

  istrlen = strlen(intstr);
  *first = parameter2intlist(intstr);
  *last  = *first;
  *inc   = 1;

  start = 1;
  *last = get_ival(intstr, *first, start, istrlen, &i);

  if ( i < istrlen )
    {
      start = i+1;
      *inc = get_ival(intstr, 1, start, istrlen, &i);
    }
}


int args2intlist(int argc, char **argv, LIST *list)
{
  int nint = 0;
  int ival;
  int first, last, inc;
  int iarg;

  for ( iarg = 0; iarg < argc; iarg++ )
    {
      split_intstring(argv[iarg], &first, &last, &inc);

      if ( inc >= 0 )
	{
	  for ( ival = first; ival <= last; ival += inc )
	    listSetInt(list, nint++, ival);
	}
      else
	{
	  for ( ival = first; ival >= last; ival += inc )
	    listSetInt(list, nint++, ival);
	}
    }

  return (nint);
}


int args2fltlist(int argc, char **argv, LIST *list)
{
  int i, nint = 0;
  int ival;
  int first, last, inc;
  int iarg;
  int len;
  double tmp_val;

  for ( iarg = 0; iarg < argc; iarg++ )
    {
      len = (int) strlen(argv[iarg]);
      for ( i = 0; i < len; i++ )
	if ( argv[iarg][i] != '/' && argv[iarg][i] != '-' && ! isdigit(argv[iarg][i]) ) break;
      
      if ( i != len )
	{
	  /*
	  if      ( strcmp(argv[iarg],  "inf") == 0 )
	    tmp_val =  DBL_MAX;
	  else if ( strcmp(argv[iarg], "-inf") == 0 )
	    tmp_val = -DBL_MAX;
	  else  
	  */                                    
	    tmp_val = parameter2double(argv[iarg]);

	  listSetFlt(list, nint++, tmp_val);
	}
      else
	{
	  split_intstring(argv[iarg], &first, &last, &inc);

	  if ( inc >= 0 )
	    {
	      for ( ival = first; ival <= last; ival += inc )
		listSetFlt(list, nint++, (double) ival);
	    }
	  else
	    {
	      for ( ival = first; ival >= last; ival += inc )
		listSetFlt(list, nint++, (double) ival);
	    }
	}
    }

  return (nint);
}
