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

      Cond       ifthen          If then
      Cond       ifnotthen       If not then
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Cond(void *argument)
{
  int IFTHEN, IFNOTTHEN;
  int operatorID;
  enum {FILL_NONE, FILL_TS, FILL_REC};
  int filltype = FILL_NONE;
  int streamID1, streamID2, streamID3;
  int gridsize;
  int ntsteps1, ntsteps2;
  int nrecs, nrecs2, nvars = 0, nlev, recID;
  int tsID;
  int varID, levelID;
  int offset;
  int vlistID1, vlistID2, vlistID3;
  int nmiss1, nmiss2, nmiss3;
  int i;
  double missval1 = -9.E33;
  double missval2 = -9.E33;
  double *array1, *array2, *array3;
  int **varnmiss1 = NULL;
  double **vardata1 = NULL;
  int taxisID2, taxisID3;

  cdoInitialize(argument);

  IFTHEN    = cdoOperatorAdd("ifthen",    0, 0, NULL);
  IFNOTTHEN = cdoOperatorAdd("ifnotthen", 0, 0, NULL);

  operatorID = cdoOperatorID();

  streamID1 = streamOpenRead(cdoStreamName(0));
  streamID2 = streamOpenRead(cdoStreamName(1));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = streamInqVlist(streamID2);
  vlistID3 = vlistDuplicate(vlistID2);

  taxisID2 = vlistInqTaxis(vlistID2);
  taxisID3 = taxisDuplicate(taxisID2);
  vlistDefTaxis(vlistID3, taxisID3);

  ntsteps1 = vlistNtsteps(vlistID1);
  ntsteps2 = vlistNtsteps(vlistID2);
  if ( ntsteps1 == 0 ) ntsteps1 = 1;
  if ( ntsteps2 == 0 ) ntsteps2 = 1;

  if ( vlistNrecs(vlistID1) == 1 && vlistNrecs(vlistID2) != 1 )
    {
      filltype = FILL_REC;
      cdoPrint("Filling up stream1 >%s< by copying the first record.", cdoStreamName(0)->args);
    }

  if ( filltype == FILL_NONE )
    vlistCompare(vlistID1, vlistID2, CMP_DIM);
  
  nospec(vlistID1);
  nospec(vlistID2);

  streamID3 = streamOpenWrite(cdoStreamName(2), cdoFiletype());

  streamDefVlist(streamID3, vlistID3);

  gridsize = vlistGridsizeMax(vlistID2);

  if ( filltype == FILL_REC && gridsize != gridInqSize(vlistGrid(vlistID1, 0)) )
    cdoAbort("Stream1 >%s< has wrong gridsize!", cdoStreamName(0)->args);

  array1 = (double*) malloc(gridsize*sizeof(double));
  array2 = (double*) malloc(gridsize*sizeof(double));
  array3 = (double*) malloc(gridsize*sizeof(double));

  if ( cdoVerbose )
    cdoPrint("Number of timesteps: file1 %d, file2 %d", ntsteps1, ntsteps2);

  if ( filltype == FILL_NONE )
    {
      if ( ntsteps1 == 1 && ntsteps2 != 1 )
	{
	  filltype = FILL_TS;
	  cdoPrint("Filling up stream1 >%s< by copying the first timestep.", cdoStreamName(0)->args);

	  nvars  = vlistNvars(vlistID1);
	  vardata1  = (double **) malloc(nvars*sizeof(double *));
	  varnmiss1 = (int **) malloc(nvars*sizeof(int *));
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	      nlev     = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      vardata1[varID]  = (double*) malloc(nlev*gridsize*sizeof(double));
	      varnmiss1[varID] = (int*) malloc(nlev*sizeof(int));
	    }
	}
    }

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID2, tsID)) )
    {
      if ( tsID == 0 || filltype == FILL_NONE )
	{
	  nrecs2 = streamInqTimestep(streamID1, tsID);
	  if ( nrecs2 == 0 )
	    cdoAbort("Input streams have different number of timesteps!");
	}

      taxisCopyTimestep(taxisID3, taxisID2);

      streamDefTimestep(streamID3, tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID2, &varID, &levelID);
	  streamReadRecord(streamID2, array2, &nmiss2);

	  if ( tsID == 0 || filltype == FILL_NONE )
	    {
	      if ( recID == 0 || filltype != FILL_REC )
		{
		  streamInqRecord(streamID1, &varID, &levelID);
		  streamReadRecord(streamID1, array1, &nmiss1);
		}

	      if ( filltype == FILL_TS )
		{
		  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
		  offset   = gridsize*levelID;
		  memcpy(vardata1[varID]+offset, array1, gridsize*sizeof(double));
		  varnmiss1[varID][levelID] = nmiss1;
		}
	    }
	  else if ( filltype == FILL_TS )
	    {
	      gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	      offset   = gridsize*levelID;
	      memcpy(array1, vardata1[varID]+offset, gridsize*sizeof(double));
	      nmiss1 = varnmiss1[varID][levelID];
	    }

	  gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
	  missval2 = vlistInqVarMissval(vlistID2, varID);
	  if ( recID == 0 || filltype != FILL_REC )
	    {
	      missval1  = vlistInqVarMissval(vlistID1, varID);
	    }

	  if ( operatorID == IFTHEN )
	    {
	      for ( i = 0; i < gridsize; i++ )
		array3[i] = !DBL_IS_EQUAL(array1[i], missval1) && !DBL_IS_EQUAL(array1[i], 0.) ? array2[i] : missval2;
	    }
	  else if ( operatorID == IFNOTTHEN )
	    {
	      for ( i = 0; i < gridsize; i++ )
		array3[i] = !DBL_IS_EQUAL(array1[i], missval1) && DBL_IS_EQUAL(array1[i], 0.) ? array2[i] : missval2;
	    }
	  else
	    {
	      cdoAbort("Operator not implemented!");
	    }

	  nmiss3 = 0;
	  for ( i = 0; i < gridsize; i++ )
	    if ( DBL_IS_EQUAL(array3[i], missval2) ) nmiss3++;

	  streamDefRecord(streamID3, varID, levelID);
	  streamWriteRecord(streamID3, array3, nmiss3);
	}

      tsID++;
    }

  streamClose(streamID3);
  streamClose(streamID2);
  streamClose(streamID1);

  if ( vardata1 )
    {
      for ( varID = 0; varID < nvars; varID++ )
	{
	  free(vardata1[varID]);
	  free(varnmiss1[varID]);
	}

      free(vardata1);
      free(varnmiss1);
    }

  if ( array3 ) free(array3);
  if ( array2 ) free(array2);
  if ( array1 ) free(array1);

  cdoFinish();

  return (0);
}
