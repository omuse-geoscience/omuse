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


#define  NALLOC_INC  1024


void *Duplicate(void *argument)
{
  int gridsize;
  int nrecs;
  int gridID, varID, levelID, recID;
  int tsID;
  int nts;
  int nalloc = 0;
  int streamID1, streamID2;
  int vlistID1, vlistID2, taxisID1, taxisID2;
  int nmiss;
  int nvars, nlevel;
  int ntsteps;
  int *vdate = NULL, *vtime = NULL;
  int idup, ndup = 1;
  field_t ***vars = NULL;

  cdoInitialize(argument);

  if ( operatorArgc() > 1 ) cdoAbort("Too many arguments!");
  else if ( operatorArgc() == 1 ) ndup = parameter2int(operatorArgv()[0]);

  if ( cdoVerbose ) cdoPrint("ndup = %d\n", ndup);

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  ntsteps  = vlistNtsteps(vlistID1);
  nvars    = vlistNvars(vlistID1);

  if ( ntsteps == 1 )
    {
      for ( varID = 0; varID < nvars; ++varID )
	if ( vlistInqVarTsteptype(vlistID1, varID) != TSTEP_CONSTANT ) break;

      if ( varID == nvars ) ntsteps = 0;
    }

  if ( ntsteps == 0 )
    {
      for ( varID = 0; varID < nvars; ++varID )
	vlistDefVarTsteptype(vlistID2, varID, TSTEP_INSTANT);
    }
 
  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  nvars = vlistNvars(vlistID1);

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      if ( tsID >= nalloc )
	{
	  nalloc += NALLOC_INC;
	  vdate = (int*) realloc(vdate, nalloc*sizeof(int));
	  vtime = (int*) realloc(vtime, nalloc*sizeof(int));
	  vars  = (field_t ***) realloc(vars, nalloc*sizeof(field_t **));
	}

      vdate[tsID] = taxisInqVdate(taxisID1);
      vtime[tsID] = taxisInqVtime(taxisID1);

      vars[tsID] = field_malloc(vlistID1, FIELD_NONE);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  gridID   = vlistInqVarGrid(vlistID1, varID);
	  gridsize = gridInqSize(gridID);
	  vars[tsID][varID][levelID].ptr = (double*) malloc(gridsize*sizeof(double));
	  streamReadRecord(streamID1, vars[tsID][varID][levelID].ptr, &nmiss);
	  vars[tsID][varID][levelID].nmiss = nmiss;
	}

      tsID++;
    }

  nts = tsID;

  for ( idup = 0; idup < ndup; idup++ )
    {
      for ( tsID = 0; tsID < nts; tsID++ )
	{
	  taxisDefVdate(taxisID2, vdate[tsID]);
	  taxisDefVtime(taxisID2, vtime[tsID]);
	  streamDefTimestep(streamID2, idup*nts+tsID);

	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      for ( levelID = 0; levelID < nlevel; levelID++ )
		{
		  if ( vars[tsID][varID][levelID].ptr )
		    {
		      nmiss = vars[tsID][varID][levelID].nmiss;
		      streamDefRecord(streamID2, varID, levelID);
		      streamWriteRecord(streamID2, vars[tsID][varID][levelID].ptr, nmiss);
		    }
		}
	    }
	}
    }

  for ( tsID = 0; tsID < nts; tsID++ ) field_free(vars[tsID], vlistID1);

  if ( vars  ) free(vars);
  if ( vdate ) free(vdate);
  if ( vtime ) free(vtime);

  streamClose(streamID2);
  streamClose(streamID1);

  cdoFinish();

  return (0);
}
