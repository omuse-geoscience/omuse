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

     Subtrend   subtrend        Subtract trend
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Subtrend(void *argument)
{
  int gridsize;
  int nrecs;
  int gridID, varID, levelID, recID;
  int tsID;
  int i;
  int streamID1, streamID2, streamID3, streamID4;
  int vlistID1, vlistID2, vlistID3, vlistID4, taxisID1, taxisID4;
  int nmiss;
  double missval, missval1, missval2;
  field_t **vars2, **vars3;
  field_t field1, field4;

  cdoInitialize(argument);

  streamID1 = streamOpenRead(cdoStreamName(0));
  streamID2 = streamOpenRead(cdoStreamName(1));
  streamID3 = streamOpenRead(cdoStreamName(2));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = streamInqVlist(streamID2);
  vlistID3 = streamInqVlist(streamID3);
  vlistID4 = vlistDuplicate(vlistID1);

  vlistCompare(vlistID1, vlistID2, CMP_DIM);
  vlistCompare(vlistID1, vlistID3, CMP_DIM);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID4 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID4, taxisID4);

  streamID4 = streamOpenWrite(cdoStreamName(3), cdoFiletype());

  streamDefVlist(streamID4, vlistID4);

  gridsize = vlistGridsizeMax(vlistID1);

  field_init(&field1);
  field_init(&field4);
  field1.ptr = (double*) malloc(gridsize*sizeof(double));
  field4.ptr = (double*) malloc(gridsize*sizeof(double));

  vars2 = field_malloc(vlistID1, FIELD_PTR);
  vars3 = field_malloc(vlistID1, FIELD_PTR);


  tsID = 0;
  nrecs = streamInqTimestep(streamID2, tsID);

  for ( recID = 0; recID < nrecs; recID++ )
    {
      streamInqRecord(streamID2, &varID, &levelID);
      streamReadRecord(streamID2, vars2[varID][levelID].ptr, &nmiss);
    }

  tsID = 0;
  nrecs = streamInqTimestep(streamID3, tsID);

  for ( recID = 0; recID < nrecs; recID++ )
    {
      streamInqRecord(streamID3, &varID, &levelID);
      streamReadRecord(streamID3, vars3[varID][levelID].ptr, &nmiss);
    }


  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID4, taxisID1);

      streamDefTimestep(streamID4, tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamReadRecord(streamID1, field1.ptr, &nmiss);

	  missval  = vlistInqVarMissval(vlistID1, varID);
	  gridID   = vlistInqVarGrid(vlistID1, varID);
	  gridsize = gridInqSize(gridID);

	  missval1 = missval;
	  missval2 = missval;
	  for ( i = 0; i < gridsize; i++ )
	    field4.ptr[i] = SUB(field1.ptr[i], ADD(vars2[varID][levelID].ptr[i], MUL(vars3[varID][levelID].ptr[i], tsID)));
    
	  nmiss = 0;
	  for ( i = 0; i < gridsize; i++ )
	    if ( DBL_IS_EQUAL(field4.ptr[i], missval) ) nmiss++;

	  streamDefRecord(streamID4, varID, levelID);
	  streamWriteRecord(streamID4, field4.ptr, nmiss);
	}

      tsID++;
    }

  field_free(vars2, vlistID1);
  field_free(vars3, vlistID1);

  if ( field1.ptr ) free(field1.ptr);
  if ( field4.ptr ) free(field4.ptr);

  streamClose(streamID4);
  streamClose(streamID3);
  streamClose(streamID2);
  streamClose(streamID1);

  cdoFinish();

  return (0);
}
