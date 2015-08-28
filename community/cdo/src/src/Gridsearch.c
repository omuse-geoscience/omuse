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

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"



void *Gridsearch(void *argument)
{
  int gridID1, gridID2;

  cdoInitialize(argument);

  cdoOperatorAdd("testpointsearch",  0,   0, NULL);
  cdoOperatorAdd("testcellsearch",   0,   0, NULL);

  operatorInputArg("source and target grid description file or name");
  operatorCheckArgc(2);
  gridID1 = cdoDefineGrid(operatorArgv()[0]);
  gridID2 = cdoDefineGrid(operatorArgv()[1]);

  cdoFinish();

  return (0);
}
