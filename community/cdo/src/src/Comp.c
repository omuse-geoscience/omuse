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

      Comp       eq              Equal
      Comp       ne              Not equal
      Comp       le              Less equal
      Comp       lt              Less than
      Comp       ge              Greater equal
      Comp       gt              Greater than
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Comp(void *argument)
{
  int EQ, NE, LE, LT, GE, GT;
  int operatorID;
  enum {FILL_NONE, FILL_TS, FILL_REC};
  int filltype = FILL_NONE;
  int streamIDx1, streamIDx2, streamID1, streamID2, streamID3;
  int gridsize, gridsize1, gridsize2;
  int nrecs, nrecs2, nvars = 0, nlev, recID;
  int tsID;
  int varID, levelID;
  int offset;
  int ntsteps1, ntsteps2;
  int vlistIDx1, vlistIDx2, vlistID1, vlistID2, vlistID3;
  int taxisIDx1, taxisID1, taxisID2, taxisID3;
  int nmiss1, nmiss2, nmiss3;
  int i;
  double missval1, missval2 = 0;
  double *missvalx1, *missvalx2;
  double *arrayx1, *arrayx2, *array1, *array2, *array3;
  double **vardata = NULL;

  cdoInitialize(argument);

  EQ = cdoOperatorAdd("eq", 0, 0, NULL);
  NE = cdoOperatorAdd("ne", 0, 0, NULL);
  LE = cdoOperatorAdd("le", 0, 0, NULL);
  LT = cdoOperatorAdd("lt", 0, 0, NULL);
  GE = cdoOperatorAdd("ge", 0, 0, NULL);
  GT = cdoOperatorAdd("gt", 0, 0, NULL);

  operatorID = cdoOperatorID();

  streamID1 = streamOpenRead(cdoStreamName(0));
  streamID2 = streamOpenRead(cdoStreamName(1));

  streamIDx1 = streamID1;
  streamIDx2 = streamID2;

  missvalx1 = &missval1;
  missvalx2 = &missval2;

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = streamInqVlist(streamID2);
  vlistIDx1 = vlistID1;
  vlistIDx2 = vlistID2;

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = vlistInqTaxis(vlistID2);
  taxisIDx1 = taxisID1;

  ntsteps1 = vlistNtsteps(vlistID1);
  ntsteps2 = vlistNtsteps(vlistID2);
  if ( ntsteps1 == 0 ) ntsteps1 = 1;
  if ( ntsteps2 == 0 ) ntsteps2 = 1;

  if ( vlistNrecs(vlistID1) != 1 && vlistNrecs(vlistID2) == 1 )
    {
      filltype = FILL_REC;
      cdoPrint("Filling up stream2 >%s< by copying the first record.", cdoStreamName(1)->args);
    }
  else if ( vlistNrecs(vlistID1) == 1 && vlistNrecs(vlistID2) != 1 )
    {
      filltype = FILL_REC;
      cdoPrint("Filling up stream1 >%s< by copying the first record.", cdoStreamName(0)->args);
      streamIDx1 = streamID2;
      streamIDx2 = streamID1;
      vlistIDx1 = vlistID2;
      vlistIDx2 = vlistID1;
      taxisIDx1 = taxisID2;
    }

  if ( filltype == FILL_NONE )
    vlistCompare(vlistID1, vlistID2, CMP_ALL);

  nospec(vlistID1);
  nospec(vlistID2);

  gridsize = vlistGridsizeMax(vlistIDx1);

  array1 = (double*) malloc(gridsize*sizeof(double));
  array2 = (double*) malloc(gridsize*sizeof(double));
  array3 = (double*) malloc(gridsize*sizeof(double));

  arrayx1 = array1;
  arrayx2 = array2;

  if ( cdoVerbose )
    cdoPrint("Number of timesteps: file1 %d, file2 %d", ntsteps1, ntsteps2);

  if ( filltype == FILL_NONE )
    {
      if ( ntsteps1 != 1 && ntsteps2 == 1 )
	{
	  filltype = FILL_TS;
	  cdoPrint("Filling up stream2 >%s< by copying the first timestep.", cdoStreamName(1)->args);
	}
      else if ( ntsteps1 == 1 && ntsteps2 != 1 )
	{
	  filltype = FILL_TS;
	  cdoPrint("Filling up stream1 >%s< by copying the first timestep.", cdoStreamName(0)->args);
	  streamIDx1 = streamID2;
          streamIDx2 = streamID1;
	  vlistIDx1 = vlistID2;
	  vlistIDx2 = vlistID1;
	  taxisIDx1 = taxisID2;
	}

      if ( filltype == FILL_TS )
	{
	  nvars  = vlistNvars(vlistIDx2);
	  vardata  = (double **) malloc(nvars*sizeof(double *));
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      gridsize = gridInqSize(vlistInqVarGrid(vlistIDx2, varID));
	      nlev     = zaxisInqSize(vlistInqVarZaxis(vlistIDx2, varID));
	      vardata[varID]  = (double*) malloc(nlev*gridsize*sizeof(double));
	    }
	}
    }

  if ( filltype != FILL_NONE && ntsteps1 == 1 )
    {
      arrayx1 = array2;
      arrayx2 = array1;
      missvalx1 = &missval2;
      missvalx2 = &missval1;
    }

  vlistID3 = vlistDuplicate(vlistIDx1);

  taxisID3 = taxisDuplicate(taxisIDx1);
  vlistDefTaxis(vlistID3, taxisID3);

  streamID3 = streamOpenWrite(cdoStreamName(2), cdoFiletype());

  streamDefVlist(streamID3, vlistID3);

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamIDx1, tsID)) )
    {
      if ( tsID == 0 || filltype == FILL_NONE )
	{
	  nrecs2 = streamInqTimestep(streamIDx2, tsID);
	  if ( nrecs2 == 0 )
	    cdoAbort("Input streams have different number of timesteps!");
	}
	  
      taxisCopyTimestep(taxisID3, taxisIDx1);

      streamDefTimestep(streamID3, tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamIDx1, &varID, &levelID);
	  streamReadRecord(streamIDx1, arrayx1, &nmiss1);

	  if ( tsID == 0 || filltype == FILL_NONE )
	    {
	      if ( recID == 0 || filltype != FILL_REC )
		{
		  streamInqRecord(streamIDx2, &varID, &levelID);
		  streamReadRecord(streamIDx2, arrayx2, &nmiss2);
		}

	      if ( filltype == FILL_TS )
		{
		  gridsize = gridInqSize(vlistInqVarGrid(vlistIDx2, varID));
		  offset   = gridsize*levelID;
		  memcpy(vardata[varID]+offset, arrayx2, gridsize*sizeof(double));
		}
	    }
	  else if ( filltype == FILL_TS )
	    {
	      gridsize = gridInqSize(vlistInqVarGrid(vlistIDx2, varID));
	      offset   = gridsize*levelID;
	      memcpy(arrayx2, vardata[varID]+offset, gridsize*sizeof(double));
	    }

	  gridsize1 = gridInqSize(vlistInqVarGrid(vlistIDx1, varID));
	  *missvalx1 = vlistInqVarMissval(vlistIDx1, varID);

	  if ( filltype == FILL_REC )
	    {
	      gridsize2 = gridInqSize(vlistInqVarGrid(vlistIDx2, 0));
	      *missvalx2 = vlistInqVarMissval(vlistIDx2, 0);
	    }
	  else
	    {
	      gridsize2 = gridInqSize(vlistInqVarGrid(vlistIDx2, varID));
	      *missvalx2 = vlistInqVarMissval(vlistIDx2, varID);
	    }

	  if ( gridsize1 != gridsize2 ) cdoAbort("Streams have different gridsize (gridsize1 = %d; gridsize2 = %d!",
						 gridsize1, gridsize2);

	  gridsize = gridsize1;

	  if ( operatorID == EQ )
	    {
	      for ( i = 0; i < gridsize; i++ )
		array3[i] = (DBL_IS_EQUAL(array1[i], missval1) || DBL_IS_EQUAL(array2[i], missval2) ?
			     missval1 : DBL_IS_EQUAL(array1[i], array2[i]));
	    }
	  else if ( operatorID == NE )
	    {
	      for ( i = 0; i < gridsize; i++ )
		array3[i] = (DBL_IS_EQUAL(array1[i], missval1) || DBL_IS_EQUAL(array2[i], missval2) ?
			     missval1 : !DBL_IS_EQUAL(array1[i], array2[i]));
	    }
	  else if ( operatorID == LE )
	    {
	      for ( i = 0; i < gridsize; i++ )
		array3[i] = (DBL_IS_EQUAL(array1[i], missval1) || DBL_IS_EQUAL(array2[i], missval2) ?
			     missval1 : array1[i] <= array2[i]);
	    }
	  else if ( operatorID == LT )
	    {
	      for ( i = 0; i < gridsize; i++ )
		array3[i] = (DBL_IS_EQUAL(array1[i], missval1) || DBL_IS_EQUAL(array2[i], missval2) ?
			     missval1 : array1[i] < array2[i]);
	    }
	  else if ( operatorID == GE )
	    {
	      for ( i = 0; i < gridsize; i++ )
		array3[i] = (DBL_IS_EQUAL(array1[i], missval1) || DBL_IS_EQUAL(array2[i], missval2) ?
			     missval1 : array1[i] >= array2[i]);
	    }
	  else if ( operatorID == GT )
	    {
	      for ( i = 0; i < gridsize; i++ )
		array3[i] = (DBL_IS_EQUAL(array1[i], missval1) || DBL_IS_EQUAL(array2[i], missval2) ?
			     missval1 : array1[i] > array2[i]);
	    }
	  else
	    {
	      cdoAbort("Operator not implemented!");
	    }

	  nmiss3 = 0;
	  for ( i = 0; i < gridsize; i++ )
	    if ( DBL_IS_EQUAL(array3[i], missval1) ) nmiss3++;

	  streamDefRecord(streamID3, varID, levelID);
	  streamWriteRecord(streamID3, array3, nmiss3);
	}

      tsID++;
    }

  streamClose(streamID3);
  streamClose(streamID2);
  streamClose(streamID1);

  if ( vardata )
    {
      for ( varID = 0; varID < nvars; varID++ ) free(vardata[varID]);
      free(vardata);
    }

  if ( array3 ) free(array3);
  if ( array2 ) free(array2);
  if ( array1 ) free(array1);

  cdoFinish();

  return (0);
}
