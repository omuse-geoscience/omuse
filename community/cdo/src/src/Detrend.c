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

      Detrend    detrend         Detrend
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


#define  NALLOC_INC  1024


static
void detrend(long nts, double missval1, double *array1, double *array2)
{
  long n;
  long j;
  double zj;
  double sumj, sumjj;
  double sumx, sumjx;
  double work1, work2;
  double missval2 = missval1;

  sumx = sumjx = 0;
  sumj = sumjj = 0;
  n = 0;
  for ( j = 0; j < nts; j++ )
    if ( !DBL_IS_EQUAL(array1[j], missval1) )
      {
        zj = j;
	sumx  += array1[j];
	sumjx += zj * array1[j];
	sumj  += zj;
	sumjj += zj * zj;
	n++;
      }

  work1 = DIV(SUB(sumjx, DIV(MUL(sumx, sumj), n) ),
	      SUB(sumjj, DIV(MUL(sumj, sumj), n)) );
  work2 = SUB(DIV(sumx, n), MUL(work1, DIV(sumj, n)));

  for ( j = 0; j < nts; j++ )
    array2[j] = SUB(array1[j], ADD(work2, MUL(j, work1)));
}


void *Detrend(void *argument)
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
  double missval;
  field_t ***vars = NULL;
  dtlist_type *dtlist = dtlist_new();
  typedef struct
  {
    double *array1;
    double *array2;
  } memory_t;
  memory_t *ompmem = NULL;

  cdoInitialize(argument);

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  nvars = vlistNvars(vlistID1);

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      if ( tsID >= nalloc )
	{
	  nalloc += NALLOC_INC;
	  vars   = (field_t ***) realloc(vars, nalloc*sizeof(field_t **));
	}

      dtlist_taxisInqTimestep(dtlist, taxisID1, tsID);

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

  ompmem = (memory_t*) malloc(ompNumThreads*sizeof(memory_t));
  for ( i = 0; i < ompNumThreads; i++ )
    {
      ompmem[i].array1 = (double*) malloc(nts*sizeof(double));
      ompmem[i].array2 = (double*) malloc(nts*sizeof(double));
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
#pragma omp parallel for default(shared) private(i, tsID)
#endif
	  for ( i = 0; i < gridsize; i++ )
	    {
              int ompthID = cdo_omp_get_thread_num();

	      for ( tsID = 0; tsID < nts; tsID++ )
		ompmem[ompthID].array1[tsID] = vars[tsID][varID][levelID].ptr[i];

	      detrend(nts, missval, ompmem[ompthID].array1, ompmem[ompthID].array2);

	      for ( tsID = 0; tsID < nts; tsID++ )
		vars[tsID][varID][levelID].ptr[i] = ompmem[ompthID].array2[tsID];
	    }
	}
    }

  for ( i = 0; i < ompNumThreads; i++ )
    {
      free(ompmem[i].array1);
      free(ompmem[i].array2);
    }
  free(ompmem);

  for ( tsID = 0; tsID < nts; tsID++ )
    {
      dtlist_taxisDefTimestep(dtlist, taxisID2, tsID);
      streamDefTimestep(streamID2, tsID);

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
		  free(vars[tsID][varID][levelID].ptr);
		  vars[tsID][varID][levelID].ptr = NULL;
		}
	    }
	}

      field_free(vars[tsID], vlistID1);
    }

  if ( vars  ) free(vars);

  dtlist_delete(dtlist);

  streamClose(streamID2);
  streamClose(streamID1);

  cdoFinish();

  return (0);
}
