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

      Timselstat    timselmin          Time range minimum
      Timselstat    timselmax          Time range maximum
      Timselstat    timselsum          Time range sum
      Timselstat    timselmean         Time range mean
      Timselstat    timselavg          Time range average
      Timselstat    timselvar          Time range variance
      Timselstat    timselvar1         Time range variance [Divisor is (n-1)]
      Timselstat    timselstd          Time range standard deviation
      Timselstat    timselstd1         Time range standard deviation [Divisor is (n-1)]
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Timselstat(void *argument)
{
  int timestat_date = TIMESTAT_MEAN;
  int nrecs = 0;
  int varID, levelID, recID;
  int tsID;
  int nsets;
  int i;
  int nmiss;
  int nlevel;

  cdoInitialize(argument);

  cdoOperatorAdd("timselmin",  func_min,  0, NULL);
  cdoOperatorAdd("timselmax",  func_max,  0, NULL);
  cdoOperatorAdd("timselsum",  func_sum,  0, NULL);
  cdoOperatorAdd("timselmean", func_mean, 0, NULL);
  cdoOperatorAdd("timselavg",  func_avg,  0, NULL);
  cdoOperatorAdd("timselvar",  func_var,  0, NULL);
  cdoOperatorAdd("timselvar1", func_var1, 0, NULL);
  cdoOperatorAdd("timselstd",  func_std,  0, NULL);
  cdoOperatorAdd("timselstd1", func_std1, 0, NULL);

  int operatorID = cdoOperatorID();
  int operfunc = cdoOperatorF1(operatorID);

  operatorInputArg("nsets <noffset <nskip>>");

  int nargc  = operatorArgc();
  int ndates = parameter2int(operatorArgv()[0]);
  int noffset = 0, nskip = 0;
  if ( nargc > 1 ) noffset = parameter2int(operatorArgv()[1]);
  if ( nargc > 2 ) nskip   = parameter2int(operatorArgv()[2]);

  if ( cdoVerbose ) cdoPrint("nsets = %d, noffset = %d, nskip = %d", ndates, noffset, nskip);

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

  int gridsize = vlistGridsizeMax(vlistID1);

  field_t field;
  field_init(&field);
  field.ptr = (double*) malloc(gridsize*sizeof(double));

  field_t **samp1 = field_malloc(vlistID1, FIELD_NONE);
  field_t **vars1 = field_malloc(vlistID1, FIELD_PTR);
  field_t **vars2 = NULL;
  if ( lvarstd ) vars2 = field_malloc(vlistID1, FIELD_PTR);

  for ( tsID = 0; tsID < noffset; tsID++ )
    {
      nrecs = streamInqTimestep(streamID1, tsID);
      if ( nrecs == 0 ) break;

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);

	  if ( tsID == 0 )
	    {
	      recVarID[recID]   = varID;
	      recLevelID[recID] = levelID;
	    }
	}
    }

  if ( tsID < noffset )
    {
      cdoWarning("noffset is larger than number of timesteps!");
      goto LABEL_END;
    }

  int otsID   = 0;
  while ( TRUE )
    {
      for ( nsets = 0; nsets < ndates; nsets++ )
	{
	  nrecs = streamInqTimestep(streamID1, tsID);
	  if ( nrecs == 0 ) break;

	  dtlist_taxisInqTimestep(dtlist, taxisID1, nsets);

	  for ( recID = 0; recID < nrecs; recID++ )
	    {
	      streamInqRecord(streamID1, &varID, &levelID);

	      if ( tsID == 0 )
		{
		  recVarID[recID]   = varID;
		  recLevelID[recID] = levelID;
		}

	      gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));

	      if ( nsets == 0 )
		{
		  streamReadRecord(streamID1, vars1[varID][levelID].ptr, &nmiss);
		  vars1[varID][levelID].nmiss = nmiss;

		  if ( nmiss > 0 || samp1[varID][levelID].ptr )
		    {
		      if ( samp1[varID][levelID].ptr == NULL )
			samp1[varID][levelID].ptr = (double*) malloc(gridsize*sizeof(double));

		      for ( i = 0; i < gridsize; i++ )
			if ( DBL_IS_EQUAL(vars1[varID][levelID].ptr[i],
					  vars1[varID][levelID].missval) )
			  samp1[varID][levelID].ptr[i] = 0;
			else
			  samp1[varID][levelID].ptr[i] = 1;
		    }
		}
	      else
		{
		  streamReadRecord(streamID1, field.ptr, &field.nmiss);
		  field.grid    = vars1[varID][levelID].grid;
		  field.missval = vars1[varID][levelID].missval;

		  if ( field.nmiss > 0 || samp1[varID][levelID].ptr )
		    {
		      if ( samp1[varID][levelID].ptr == NULL )
			{
			  samp1[varID][levelID].ptr = (double*) malloc(gridsize*sizeof(double));
			  for ( i = 0; i < gridsize; i++ )
			    samp1[varID][levelID].ptr[i] = nsets;
			}

		      for ( i = 0; i < gridsize; i++ )
			if ( !DBL_IS_EQUAL(field.ptr[i], vars1[varID][levelID].missval) )
			  samp1[varID][levelID].ptr[i]++;
		    }

		  if ( lvarstd )
		    {
		      farsumq(&vars2[varID][levelID], field);
		      farsum(&vars1[varID][levelID], field);
		    }
		  else
		    {
		      farfun(&vars1[varID][levelID], field, operfunc);
		    }
		}
	    }

	  if ( nsets == 0 && lvarstd )
	    for ( varID = 0; varID < nvars; varID++ )
	      {
		if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;
		nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
		for ( levelID = 0; levelID < nlevel; levelID++ )
		  farmoq(&vars2[varID][levelID], vars1[varID][levelID]);
	      }

	  tsID++;
	}

      if ( nrecs == 0 && nsets == 0 ) break;

      if ( lmean )
	for ( varID = 0; varID < nvars; varID++ )
	  {
	    if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;
	    nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	    for ( levelID = 0; levelID < nlevel; levelID++ )
	      {
		if ( samp1[varID][levelID].ptr == NULL )
		  farcmul(&vars1[varID][levelID], 1.0/nsets);
		else
		  fardiv(&vars1[varID][levelID], samp1[varID][levelID]);
	      }
	  }
      else if ( lvarstd )
	for ( varID = 0; varID < nvars; varID++ )
	  {
	    if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;
	    nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	    for ( levelID = 0; levelID < nlevel; levelID++ )
	      {
		if ( samp1[varID][levelID].ptr == NULL )
		  {
		    if ( lstd )
		      farcstd(&vars1[varID][levelID], vars2[varID][levelID], nsets, divisor);
		    else
		      farcvar(&vars1[varID][levelID], vars2[varID][levelID], nsets, divisor);
		  }
		else
		  {
		    if ( lstd )
		      farstd(&vars1[varID][levelID], vars2[varID][levelID], samp1[varID][levelID], divisor);
		    else
		      farvar(&vars1[varID][levelID], vars2[varID][levelID], samp1[varID][levelID], divisor);
		  }
	      }
	  }

      dtlist_stat_taxisDefTimestep(dtlist, taxisID2, nsets);
      streamDefTimestep(streamID2, otsID);

      for ( recID = 0; recID < nrecords; recID++ )
	{
	  varID   = recVarID[recID];
	  levelID = recLevelID[recID];

	  if ( otsID && vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;

	  streamDefRecord(streamID2, varID, levelID);
	  streamWriteRecord(streamID2, vars1[varID][levelID].ptr,  vars1[varID][levelID].nmiss);
	}

      if ( nrecs == 0 ) break;
      otsID++;

      for ( i = 0; i < nskip; i++ )
	{
	  nrecs = streamInqTimestep(streamID1, tsID);
	  if ( nrecs == 0 ) break;
	  tsID++;
	}

      if ( nrecs == 0 ) break;
    }

 LABEL_END:


  field_free(vars1, vlistID1);
  field_free(samp1, vlistID1);
  if ( lvarstd ) field_free(vars2, vlistID1);

  dtlist_delete(dtlist);

  if ( field.ptr ) free(field.ptr);

  if ( recVarID   ) free(recVarID);
  if ( recLevelID ) free(recLevelID);

  streamClose(streamID2);
  streamClose(streamID1);

  cdoFinish();

  return (0);
}
