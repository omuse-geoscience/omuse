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

      Yearmonstat   yearmonmean        Yearly mean from monthly data
      Yearmonstat   yearmonavg         Yearly average from monthly data
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Yearmonstat(void *argument)
{
  int timestat_date = TIMESTAT_MEAN;
  int gridsize;
  int vdate = 0, vtime = 0;
  int vdate0 = 0, vtime0 = 0;
  int nrecs;
  int varID, levelID, recID;
  int tsID;
  int otsID;
  int i;
  int dpm;
  int year0 = 0, month0 = 0;
  int year, month, day;
  int nmiss;
  int nlevel;
  long nsets;
  double dsets;
  char vdatestr[32], vtimestr[32];
  field_t **vars1 = NULL, **samp1 = NULL;
  field_t field;

  cdoInitialize(argument);

  cdoOperatorAdd("yearmonmean",  func_mean, 0, NULL);
  cdoOperatorAdd("yearmonavg",   func_avg,  0, NULL);

  int operatorID = cdoOperatorID();
  int operfunc = cdoOperatorF1(operatorID);

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

  int calendar = taxisInqCalendar(taxisID1);
  dtlist_type *dtlist = dtlist_new();
  dtlist_set_stat(dtlist, timestat_date);
  dtlist_set_calendar(dtlist, calendar);

  gridsize = vlistGridsizeMax(vlistID1);

  field_init(&field);
  field.ptr = (double*) malloc(gridsize*sizeof(double));

  vars1 = field_malloc(vlistID1, FIELD_PTR);
  samp1 = field_malloc(vlistID1, FIELD_NONE);

  tsID    = 0;
  otsID   = 0;
  while ( TRUE )
    {
      nsets = 0;
      dsets = 0;
      while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
	{
	  dtlist_taxisInqTimestep(dtlist, taxisID1, nsets);
	  vdate = dtlist_get_vdate(dtlist, nsets);
	  vtime = dtlist_get_vtime(dtlist, nsets);
	  cdiDecodeDate(vdate, &year, &month, &day);

	  if ( nsets == 0 ) year0 = year;

	  if ( year != year0 ) break;

	  if ( nsets > 0 && month == month0 )
	    {
	      date2str(vdate0, vdatestr, sizeof(vdatestr));
	      time2str(vtime0, vtimestr, sizeof(vtimestr));
	      cdoWarning("   last timestep: %s %s", vdatestr, vtimestr);
	      date2str(vdate, vdatestr, sizeof(vdatestr));
	      time2str(vtime, vtimestr, sizeof(vtimestr));
	      cdoWarning("current timestep: %s %s", vdatestr, vtimestr);
	      cdoAbort("Month does not change!");
	    }

	  dpm = days_per_month(calendar, year, month);

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

		  farcmul(&vars1[varID][levelID], dpm);

		  if ( nmiss > 0 || samp1[varID][levelID].ptr )
		    {
		      if ( samp1[varID][levelID].ptr == NULL )
			samp1[varID][levelID].ptr = (double*) malloc(gridsize*sizeof(double));

		      for ( i = 0; i < gridsize; i++ )
			if ( DBL_IS_EQUAL(vars1[varID][levelID].ptr[i], vars1[varID][levelID].missval) )
			  samp1[varID][levelID].ptr[i] = 0;
			else
			  samp1[varID][levelID].ptr[i] = dpm;
		    }
		}
	      else
		{
		  streamReadRecord(streamID1, field.ptr, &field.nmiss);
		  field.grid    = vars1[varID][levelID].grid;
		  field.missval = vars1[varID][levelID].missval;

		  farcmul(&field, dpm);

		  if ( field.nmiss > 0 || samp1[varID][levelID].ptr )
		    {
		      if ( samp1[varID][levelID].ptr == NULL )
			{
			  samp1[varID][levelID].ptr = (double*) malloc(gridsize*sizeof(double));
			  for ( i = 0; i < gridsize; i++ )
			    samp1[varID][levelID].ptr[i] = dsets;
			}

		      for ( i = 0; i < gridsize; i++ )
			if ( !DBL_IS_EQUAL(field.ptr[i], vars1[varID][levelID].missval) )
			  samp1[varID][levelID].ptr[i] += dpm;
		    }

		  farfun(&vars1[varID][levelID], field, operfunc);
		}
	    }

	  month0 = month;
	  vdate0 = vdate;
	  vtime0 = vtime;
	  nsets++;
	  dsets += dpm;
	  tsID++;
	}

      if ( nrecs == 0 && nsets == 0 ) break;

      for ( varID = 0; varID < nvars; varID++ )
	{
	  if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;
	  nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	  for ( levelID = 0; levelID < nlevel; levelID++ )
	    {
	      if ( samp1[varID][levelID].ptr == NULL )
		farcmul(&vars1[varID][levelID], 1.0/dsets);
	      else
		fardiv(&vars1[varID][levelID], samp1[varID][levelID]);
	    }
	}

      if ( cdoVerbose )
	{
	  date2str(vdate0, vdatestr, sizeof(vdatestr));
	  time2str(vtime0, vtimestr, sizeof(vtimestr));
	  cdoPrint("%s %s  nsets = %d", vdatestr, vtimestr, nsets);
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
    }


  field_free(vars1, vlistID1);
  field_free(samp1, vlistID1);

  if ( field.ptr ) free(field.ptr);

  if ( recVarID   ) free(recVarID);
  if ( recLevelID ) free(recLevelID);

  dtlist_delete(dtlist);

  streamClose(streamID2);
  streamClose(streamID1);

  cdoFinish();

  return (0);
}
