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

      Condc      ifthenc         If then constant
      Condc      ifnotthenc      If not then constant
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Condc(void *argument)
{
  int IFTHENC, IFNOTTHENC;
  int operatorID;
  int streamID1, streamID2;
  int gridsize;
  int nrecs, recID;
  int tsID;
  int varID, levelID;
  int vlistID1, vlistID2;
  int nmiss, nmiss2;
  int i;
  double missval;
  double rc;
  double *array1, *array2;
  int taxisID1, taxisID2;

  cdoInitialize(argument);

  IFTHENC    = cdoOperatorAdd("ifthenc",    0, 0, NULL);
  IFNOTTHENC = cdoOperatorAdd("ifnotthenc", 0, 0, NULL);

  operatorID = cdoOperatorID();

  operatorInputArg("constant value");
  rc = parameter2double(operatorArgv()[0]);

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  nospec(vlistID1);

  gridsize = vlistGridsizeMax(vlistID1);

  array1 = (double*) malloc(gridsize*sizeof(double));
  array2 = (double*) malloc(gridsize*sizeof(double));

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamReadRecord(streamID1, array1, &nmiss);

	  missval  = vlistInqVarMissval(vlistID1, varID);
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));

	  if ( operatorID == IFTHENC )
	    {
	      for ( i = 0; i < gridsize; i++ )
		array2[i] = !DBL_IS_EQUAL(array1[i], missval) && !DBL_IS_EQUAL(array1[i], 0.) ? rc : missval;
	    }
	  else if ( operatorID == IFNOTTHENC )
	    {
	      for ( i = 0; i < gridsize; i++ )
		array2[i] = !DBL_IS_EQUAL(array1[i], missval) && DBL_IS_EQUAL(array1[i], 0.) ? rc : missval;
	    }
	  else
	    {
	      cdoAbort("Operator not implemented!");
	    }

	  nmiss2 = 0;
	  for ( i = 0; i < gridsize; i++ )
	    if ( DBL_IS_EQUAL(array2[i], missval) ) nmiss2++;

	  streamDefRecord(streamID2, varID, levelID);
	  streamWriteRecord(streamID2, array2, nmiss2);
	}
      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  if ( array2 ) free(array2);
  if ( array1 ) free(array1);

  cdoFinish();

  return (0);
}
