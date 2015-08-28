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

#ifndef _NAMELIST_H
#define _NAMELIST_H

#include <stdio.h>


#define  NML_INT         1
#define  NML_FLT      2
#define  NML_WORD        3
#define  NML_TEXT        4

#define NML_DEF_INT(name, size, val)  int sel##name[size]; int nsel##name = 0; int name = 0
#define NML_DEF_FLT(name, size, val)  double sel##name[size]; int nsel##name = 0; double name = 0
#define NML_ADD_INT(nml, name) namelistAdd(nml, #name, NML_INT, 0, sel##name, sizeof(sel##name)/sizeof(int))
#define NML_ADD_FLT(nml, name) namelistAdd(nml, #name, NML_FLT, 0, sel##name, sizeof(sel##name)/sizeof(double))
#define NML_NUM(nml, name)   nsel##name = namelistNum(nml, #name)
#define NML_PAR(name)        nsel##name, sel##name, name

#define MAX_NML_ENTRY  256

#define MAX_LINE_LEN  4096

typedef struct
{
  int nptype, namitf, namitl;
  char lineac[MAX_LINE_LEN], lineuc[MAX_LINE_LEN], linelc[MAX_LINE_LEN];
} nml_line_t;

typedef struct
{
  char  *name;
  void  *ptr;
  int    type;
  int    occ;
  int    dis;
  size_t size;
} nml_entry_t;

typedef struct
{
  int          size;
  int          dis;
  char        *name;
  nml_line_t   line;
  nml_entry_t *entry[MAX_NML_ENTRY];
} namelist_t;


namelist_t *namelistNew(const char *name);
void namelistDelete(namelist_t *nml);
void namelistReset(namelist_t *nml);
int  namelistAdd(namelist_t *nml, const char *name, int type, int dis, void *ptr, size_t size);
void namelistPrint(namelist_t *nml);
void namelistRead(FILE *nmlfp, namelist_t *nml);
int  namelistNum(namelist_t *nml, const char *name);

#endif  /* _NAMELIST_H */
