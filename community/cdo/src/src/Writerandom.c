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

      Writerandom writerandom
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Writerandom(void *argument)
{
  int streamID1, streamID2;
  int vlistID1, vlistID2;
  int gridsize;
  int recID, nrecs;
  int tsID, varID, levelID;
  int index, rindex, ipos;
  double **recdata = NULL;
  int *recvarID, *reclevelID, *recnmiss, *recindex;
  int taxisID1, taxisID2;


  cdoInitialize(argument);

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);

      recdata    = (double**) malloc(nrecs*sizeof(double*));
      recvarID   = (int*) malloc(nrecs*sizeof(int));
      reclevelID = (int*) malloc(nrecs*sizeof(int));
      recnmiss   = (int*) malloc(nrecs*sizeof(int));
      recindex   = (int*) malloc(nrecs*sizeof(int));

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	  recvarID[recID] = varID;
	  reclevelID[recID] = levelID;
	  recdata[recID] = (double*) malloc(gridsize*sizeof(double));
	  streamReadRecord(streamID1, recdata[recID], &recnmiss[recID]);
	}

      for ( recID = 0; recID < nrecs; recID++ ) recindex[recID] = -1;

      for ( rindex = nrecs-1; rindex >= 0; rindex-- )
	{
	  index = (int) (rindex*1.0*rand()/(RAND_MAX+1.0));
	  /*	printf("rindex %d %d\n", rindex, index); */
	  ipos = -1;
	  for ( recID = 0; recID < nrecs; recID++ )
	    {
	      if ( recindex[recID] == -1 ) ipos++;
	      if ( recindex[recID] == -1 && ipos == index )
		{
		  recindex[recID] = rindex;
		  break;
		}
	    }
	}

      /*
      for ( recID = 0; recID < nrecs; recID++ )
	printf("recID %d %d\n", recID, recindex[recID]);
      */
      for ( recID = 0; recID < nrecs; recID++ )
	if ( recindex[recID] == -1 )
	  cdoAbort("Internal problem! Random initialize.");

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  rindex   = recindex[recID];
	  varID    = recvarID[rindex];
	  levelID  = reclevelID[rindex];
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
	  streamDefRecord(streamID2, varID, levelID);
	  streamWriteRecord(streamID2, recdata[rindex], recnmiss[rindex]);
	}

      for ( recID = 0; recID < nrecs; recID++ ) free(recdata[recID]);

      free(recdata);
      free(recvarID);
      free(reclevelID);
      free(recnmiss);
      free(recindex);

      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  cdoFinish();

  return (0);
}
