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

      Tstepcount  tstepcount  Count number of timesteps
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


#define  NALLOC_INC  1024


static
double tstepcount(long nts, double missval1, double *array1, double refval)
{
  long j;
  long n = 0;

  if ( DBL_IS_EQUAL(refval, missval1) ) return (missval1);

  for ( j = 0; j < nts; j++ )
    {
      n++;
      if ( DBL_IS_EQUAL(array1[j], refval) ) break;  
    }

  if ( j == nts )
    return (missval1);
  else
    return ((double) n);
}


void *Tstepcount(void *argument)
{
  int gridsize;
  int nrecs;
  int gridID, varID, levelID, recID;
  int tsID;
  int i;
  int nts;
  int nalloc = 0;
  int streamID1, streamID2;
  int vlistID1, vlistID2, taxisID1, taxisID2;
  int nmiss;
  int nvars, nlevel;
  int vdate = 0, vtime = 0;
  double missval;
  double refval = 0;
  double count;
  field_t ***vars = NULL;
  typedef struct
  {
    double *array1;
  } memory_t;
  memory_t *mem = NULL;

  cdoInitialize(argument);

  if ( operatorArgc() == 1 ) refval = parameter2double(operatorArgv()[0]);

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  vlistDefNtsteps(vlistID2, 1);

  nvars = vlistNvars(vlistID1);
  for ( varID = 0; varID < nvars; varID++ )
    {
      vlistDefVarUnits(vlistID2, varID, "steps");
    }

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      if ( tsID >= nalloc )
	{
	  nalloc += NALLOC_INC;
	  vars  = (field_t ***) realloc(vars, nalloc*sizeof(field_t **));
	}

      vdate = taxisInqVdate(taxisID1);
      vtime = taxisInqVtime(taxisID1);

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

  mem = (memory_t*) malloc(ompNumThreads*sizeof(memory_t));
  for ( i = 0; i < ompNumThreads; i++ )
    {
      mem[i].array1 = (double*) malloc(nts*sizeof(double));
    }

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridID   = vlistInqVarGrid(vlistID1, varID);
      missval  = vlistInqVarMissval(vlistID1, varID);
      gridsize = gridInqSize(gridID);
      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      for ( levelID = 0; levelID < nlevel; levelID++ )
	{
#if defined(_OPENMP)
#pragma omp parallel for default(shared) private(i, tsID) \
  schedule(dynamic,1)
#endif
	  for ( i = 0; i < gridsize; i++ )
	    {
	      int ompthID = cdo_omp_get_thread_num();

	      for ( tsID = 0; tsID < nts; tsID++ )
		mem[ompthID].array1[tsID] = vars[tsID][varID][levelID].ptr[i];

	      count = tstepcount(nts, missval, mem[ompthID].array1, refval);
	      vars[0][varID][levelID].ptr[i] = count;
	    }
	}
    }

  for ( i = 0; i < ompNumThreads; i++ )
    {
      free(mem[i].array1);
    }
  free(mem);

  taxisDefVdate(taxisID2, vdate);
  taxisDefVtime(taxisID2, vtime);
  streamDefTimestep(streamID2, 0);

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridID   = vlistInqVarGrid(vlistID2, varID);
      missval  = vlistInqVarMissval(vlistID2, varID);
      gridsize = gridInqSize(gridID);
      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));

      for ( levelID = 0; levelID < nlevel; levelID++ )
	{
	  streamDefRecord(streamID2, varID, levelID);

	  nmiss = 0;
	  for ( i = 0; i < gridsize; i++ )
	    if ( DBL_IS_EQUAL(vars[0][varID][levelID].ptr[i], missval) ) nmiss++;

	  streamWriteRecord(streamID2, vars[0][varID][levelID].ptr, nmiss);
	}
    }

  for ( tsID = 0; tsID < nts; tsID++ ) field_free(vars[tsID], vlistID1);

  if ( vars  ) free(vars);

  streamClose(streamID2);
  streamClose(streamID1);

  cdoFinish();

  return (0);
}
