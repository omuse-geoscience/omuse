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


void *Tee(void *argument)
{
  int streamID1, streamID2, streamID3;
  int nrecs;
  int tsID, recID, varID, levelID;
  int lcopy = FALSE;
  int gridsize;
  int vlistID1, vlistID2, vlistID3;
  int nmiss;
  int taxisID1, taxisID2, taxisID3;
  double *array = NULL;

  cdoInitialize(argument);

  if ( UNCHANGED_RECORD ) lcopy = TRUE;

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  taxisID1 = vlistInqTaxis(vlistID1);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  streamID3 = streamOpenWrite(cdoStreamName(2), cdoFiletype());

  vlistID2 = vlistDuplicate(vlistID1);
  vlistID3 = vlistDuplicate(vlistID1);

  taxisID2 = taxisDuplicate(taxisID1);
  taxisID3 = taxisDuplicate(taxisID1);

  vlistDefTaxis(vlistID2, taxisID2);
  vlistDefTaxis(vlistID3, taxisID3);

  streamDefVlist(streamID2, vlistID2);
  streamDefVlist(streamID3, vlistID3);

  gridsize = vlistGridsizeMax(vlistID1);
  array = (double*) malloc(gridsize*sizeof(double));

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      taxisCopyTimestep(taxisID3, taxisID1);

      streamDefTimestep(streamID2, tsID);
      streamDefTimestep(streamID3, tsID);
	       
      for ( recID = 0; recID < nrecs; recID++ )
	{ 
	  if ( lcopy )
	    {
	      streamInqRecord(streamID1, &varID, &levelID);

	      streamDefRecord(streamID2,  varID,  levelID);
	      streamCopyRecord(streamID2, streamID1);

	      streamDefRecord(streamID3,  varID,  levelID);
	      streamCopyRecord(streamID3, streamID1);
	    }
	  else
	    {
	      streamInqRecord(streamID1, &varID, &levelID);
	      streamReadRecord(streamID1, array, &nmiss);

	      streamDefRecord(streamID2,  varID,  levelID);
	      streamWriteRecord(streamID2, array, nmiss);

	      streamDefRecord(streamID3,  varID,  levelID);
	      streamWriteRecord(streamID3, array, nmiss);
	    }
	}

      tsID++;
    }

  streamClose(streamID1);

  streamClose(streamID2);
  streamClose(streamID3);

  if ( array ) free(array);

  cdoFinish();

  return (0);
}
