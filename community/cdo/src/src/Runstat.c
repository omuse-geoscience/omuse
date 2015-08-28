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

      Runstat    runmin          Running minimum
      Runstat    runmax          Running maximum
      Runstat    runsum          Running sum
      Runstat    runmean         Running mean
      Runstat    runavg          Running average
      Runstat    runvar          Running variance
      Runstat    runvar1         Running variance [Divisor is (n-1)]
      Runstat    runstd          Running standard deviation
      Runstat    runstd1         Running standard deviation [Divisor is (n-1)]
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Runstat(void *argument)
{
  int timestat_date = TIMESTAT_MEAN;
  int gridsize;
  int i;
  int varID;
  int recID;
  int nrecs;
  int levelID;
  int tsID;
  int otsID;
  int inp, its;
  int nmiss;
  int nlevel;
  int runstat_nomiss = 0;
  double missval;
  field_t ***vars1 = NULL, ***vars2 = NULL, ***samp1 = NULL;

  cdoInitialize(argument);

  char *envstr = getenv("RUNSTAT_NOMISS");
  if ( envstr )
    {
      char *endptr;
      int envval = (int) strtol(envstr, &endptr, 10);
      if ( envval == 1 ) runstat_nomiss = 1;
    }

  cdoOperatorAdd("runmin",  func_min,  0, NULL);
  cdoOperatorAdd("runmax",  func_max,  0, NULL);
  cdoOperatorAdd("runsum",  func_sum,  0, NULL);
  cdoOperatorAdd("runmean", func_mean, 0, NULL);
  cdoOperatorAdd("runavg",  func_avg,  0, NULL);
  cdoOperatorAdd("runvar",  func_var,  0, NULL);
  cdoOperatorAdd("runvar1", func_var1, 0, NULL);
  cdoOperatorAdd("runstd",  func_std,  0, NULL);
  cdoOperatorAdd("runstd1", func_std1, 0, NULL);

  int operatorID = cdoOperatorID();
  int operfunc = cdoOperatorF1(operatorID);

  operatorInputArg("number of timesteps");
  int ndates = parameter2int(operatorArgv()[0]);

  int lmean   = operfunc == func_mean || operfunc == func_avg;
  int lstd    = operfunc == func_std || operfunc == func_std1;
  int lvarstd = operfunc == func_std || operfunc == func_var || operfunc == func_std1 || operfunc == func_var1;
  double divisor = operfunc == func_std1 || operfunc == func_var1;

  int streamID1 = streamOpenRead(cdoStreamName(0));

  int vlistID1 = streamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  int nvars    = vlistNvars(vlistID1);
  int nrecords = vlistNrecs(vlistID1);

  int *recVarID   = (int*) malloc(nrecords*sizeof(int));
  int *recLevelID = (int*) malloc(nrecords*sizeof(int));

  dtlist_type *dtlist = dtlist_new();
  dtlist_set_stat(dtlist, timestat_date);
  dtlist_set_calendar(dtlist, taxisInqCalendar(taxisID1));

  vars1 = (field_t ***) malloc((ndates+1)*sizeof(field_t **));
  if ( !runstat_nomiss )
    samp1 = (field_t ***) malloc((ndates+1)*sizeof(field_t **));
  if ( lvarstd )
    vars2 = (field_t ***) malloc((ndates+1)*sizeof(field_t **));

  for ( its = 0; its < ndates; its++ )
    {
      vars1[its] = field_malloc(vlistID1, FIELD_PTR);
      if ( !runstat_nomiss )
	samp1[its] = field_malloc(vlistID1, FIELD_PTR);
      if ( lvarstd )
	vars2[its] = field_malloc(vlistID1, FIELD_PTR);
    }

  int gridsizemax = vlistGridsizeMax(vlistID1);
  int *imask = (int*) malloc(gridsizemax*sizeof(int));

  for ( tsID = 0; tsID < ndates; tsID++ )
    {
      nrecs = streamInqTimestep(streamID1, tsID);
      if ( nrecs == 0 ) cdoAbort("File has less then %d timesteps!", ndates);

      dtlist_taxisInqTimestep(dtlist, taxisID1, tsID);
	
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);

	  if ( tsID == 0 )
	    {
	      recVarID[recID]   = varID;
	      recLevelID[recID] = levelID;
	    }
	  
	  streamReadRecord(streamID1, vars1[tsID][varID][levelID].ptr, &nmiss);
	  vars1[tsID][varID][levelID].nmiss = nmiss;

	  if ( runstat_nomiss && nmiss > 0 ) cdoAbort("Missing values supported swichted off!");

	  if ( !runstat_nomiss )
	    {
	      gridsize = gridInqSize(vars1[0][varID][levelID].grid);
	      missval  = vars1[0][varID][levelID].missval;

	      for ( i = 0; i < gridsize; i++ )
		if ( DBL_IS_EQUAL(vars1[tsID][varID][levelID].ptr[i], missval) )
		  imask[i] = 0;
		else
		  imask[i] = 1;

	      for ( i = 0; i < gridsize; i++ )
		samp1[tsID][varID][levelID].ptr[i] = (double) imask[i];

#if defined(_OPENMP)
#pragma omp parallel for default(shared) private(i, inp)
#endif
	      for ( inp = 0; inp < tsID; inp++ )
		{
		  double *ptr = samp1[inp][varID][levelID].ptr;
		  for ( i = 0; i < gridsize; i++ )
		    if ( imask[i] > 0 ) ptr[i]++;
		}
	    }

	  if ( lvarstd )
	    {
	      farmoq(&vars2[tsID][varID][levelID], vars1[tsID][varID][levelID]);
#if defined(_OPENMP)
#pragma omp parallel for default(shared)
#endif
	      for ( inp = 0; inp < tsID; inp++ )
		{
		  farsumq(&vars2[inp][varID][levelID], vars1[tsID][varID][levelID]);
		  farsum(&vars1[inp][varID][levelID], vars1[tsID][varID][levelID]);
		}
	    }
	  else
	    {
#if defined(_OPENMP)
#pragma omp parallel for default(shared)
#endif
	      for ( inp = 0; inp < tsID; inp++ )
		{
		  farfun(&vars1[inp][varID][levelID], vars1[tsID][varID][levelID], operfunc);
		}
	    }
	}
    }

  otsID = 0;
  while ( TRUE )
    {
      if ( lmean )
	for ( varID = 0; varID < nvars; varID++ )
	  {
	    if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;
	    nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	    for ( levelID = 0; levelID < nlevel; levelID++ )
	      {
		if ( runstat_nomiss )
		  farcmul(&vars1[0][varID][levelID], 1.0/ndates);
		else
		  fardiv(&vars1[0][varID][levelID], samp1[0][varID][levelID]);
	      }
	  }
      else if ( lvarstd )
	for ( varID = 0; varID < nvars; varID++ )
	  {
	    if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;
	    nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	    for ( levelID = 0; levelID < nlevel; levelID++ )
	      {
		if ( runstat_nomiss )
		  {
		    if ( lstd )
		      farcstd(&vars1[0][varID][levelID], vars2[0][varID][levelID], ndates, divisor);
		    else
		      farcvar(&vars1[0][varID][levelID], vars2[0][varID][levelID], ndates, divisor);
		  }
		else
		  {
		    if ( lstd )
		      farstd(&vars1[0][varID][levelID], vars2[0][varID][levelID], samp1[0][varID][levelID], divisor);
		    else
		      farvar(&vars1[0][varID][levelID], vars2[0][varID][levelID], samp1[0][varID][levelID], divisor);
		  }
	      }
	  }

      dtlist_stat_taxisDefTimestep(dtlist, taxisID2, ndates);
      streamDefTimestep(streamID2, otsID);

      for ( recID = 0; recID < nrecords; recID++ )
	{
	  varID    = recVarID[recID];
	  levelID  = recLevelID[recID];

	  if ( otsID && vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;

	  streamDefRecord(streamID2, varID, levelID);
	  streamWriteRecord(streamID2, vars1[0][varID][levelID].ptr, vars1[0][varID][levelID].nmiss);
	}

      otsID++;

      dtlist_shift(dtlist);

      vars1[ndates] = vars1[0];
      if ( !runstat_nomiss )
	samp1[ndates] = samp1[0];
      if ( lvarstd )
        vars2[ndates] = vars2[0];

      for ( inp = 0; inp < ndates; inp++ )
	{
	  vars1[inp] = vars1[inp+1];
	  if ( !runstat_nomiss )
	    samp1[inp] = samp1[inp+1];
	  if ( lvarstd )
	    vars2[inp] = vars2[inp+1];
	}

      nrecs = streamInqTimestep(streamID1, tsID);
      if ( nrecs == 0 ) break;

      dtlist_taxisInqTimestep(dtlist, taxisID1, ndates-1);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  
	  streamReadRecord(streamID1, vars1[ndates-1][varID][levelID].ptr, &nmiss);
	  vars1[ndates-1][varID][levelID].nmiss = nmiss;

	  if ( runstat_nomiss && nmiss > 0 ) cdoAbort("Missing values supported swichted off!");

	  if ( !runstat_nomiss )
	    {
	      gridsize = gridInqSize(vars1[0][varID][levelID].grid);
	      missval  = vars1[0][varID][levelID].missval;

	      for ( i = 0; i < gridsize; i++ )
		if ( DBL_IS_EQUAL(vars1[ndates-1][varID][levelID].ptr[i], missval) )
		  imask[i] = 0;
		else
		  imask[i] = 1;

	      for ( i = 0; i < gridsize; i++ )
		samp1[ndates-1][varID][levelID].ptr[i] = (double) imask[i];

#if defined(_OPENMP)
#pragma omp parallel for default(shared) private(i, inp)
#endif
	      for ( inp = 0; inp < ndates-1; inp++ )
		{
		  double *ptr = samp1[inp][varID][levelID].ptr;
		  for ( i = 0; i < gridsize; i++ )
		    if ( imask[i] > 0 ) ptr[i]++;
		}
	    }

	  if ( lvarstd )
	    {
	      farmoq(&vars2[ndates-1][varID][levelID], vars1[ndates-1][varID][levelID]);
#if defined(_OPENMP)
#pragma omp parallel for default(shared)
#endif
	      for ( inp = 0; inp < ndates-1; inp++ )
		{
		  farsumq(&vars2[inp][varID][levelID], vars1[ndates-1][varID][levelID]);
		  farsum(&vars1[inp][varID][levelID], vars1[ndates-1][varID][levelID]);
		}
	    }
	  else
	    {
#if defined(_OPENMP)
#pragma omp parallel for default(shared)
#endif
	      for ( inp = 0; inp < ndates-1; inp++ )
		{
		  farfun(&vars1[inp][varID][levelID], vars1[ndates-1][varID][levelID], operfunc);
		}
	    }
	}

      tsID++;
    }

  for ( its = 0; its < ndates; its++ )
    {
      field_free(vars1[its], vlistID1);
      if ( !runstat_nomiss ) field_free(samp1[its], vlistID1);
      if ( lvarstd ) field_free(vars2[its], vlistID1);
    }

  free(vars1);
  if ( !runstat_nomiss ) free(samp1);
  if ( lvarstd ) free(vars2);

  if ( recVarID   ) free(recVarID);
  if ( recLevelID ) free(recLevelID);
  if ( imask )      free(imask);

  dtlist_delete(dtlist);

  streamClose(streamID2);
  streamClose(streamID1);

  cdoFinish();

  return (0);
}
