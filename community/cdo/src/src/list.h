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

#ifndef _LIST_H
#define _LIST_H

#define  INT_LIST  1
#define  FLT_LIST  2


typedef struct {
  void *array;
  int nalloc;
  int allinc;
  int type;
}
LIST;


LIST *listNew(int type);
void listDelete(LIST *list);
void *listArrayPtr(LIST *list);
void listSetInt(LIST *list, int num, int ival);
void listSetFlt(LIST *list, int num, double fval);
int listGetInt(LIST *list, int num);
double listGetFlt(LIST *list, int num);
int args2intlist(int argc, char **argv, LIST *list);
int args2fltlist(int argc, char **argv, LIST *list);

#endif  /* _LIST_H */
