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

/*
   This module contains the following operators:

*/


#include "cdo.h"
#include "cdo_int.h"
#include "namelist.h"


void *Nmltest(void *argument)
{
  namelist_t *nml;
  int i1[5] = {-99, -99, -99, -99, -99};
  int i2    = -99;
  char lop[99] = "";
  double dm = 0;
  char *var[3];

  cdoInitialize(argument);

  nml = namelistNew("SELECT");

  namelistAdd(nml, "i1",  NML_INT,    0, i1,   sizeof(i1)/sizeof(int));
  namelistAdd(nml, "i2",  NML_INT,    1, &i2,  sizeof(i2)/sizeof(int));
  namelistAdd(nml, "lop", NML_TEXT,   2, lop,  sizeof(lop)/sizeof(char));
  namelistAdd(nml, "dm",  NML_FLT, 1, &dm,  sizeof(dm)/sizeof(double));
  namelistAdd(nml, "var", NML_WORD,   0, var,  sizeof(var)/sizeof(char *));

  namelistRead(stdin, nml);

  namelistPrint(nml);

  namelistDelete(nml);

  cdoFinish();

  return (0);
}
