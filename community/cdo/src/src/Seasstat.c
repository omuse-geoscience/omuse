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

      Seasstat   seasmin         Seasonal minimum
      Seasstat   seasmax         Seasonal maximum
      Seasstat   seassum         Seasonal sum
      Seasstat   seasmean        Seasonal mean
      Seasstat   seasavg         Seasonal average
      Seasstat   seasvar         Seasonal variance
      Seasstat   seasvar1        Seasonal variance [Divisor is (n-1)]
      Seasstat   seasstd         Seasonal standard deviation
      Seasstat   seasstd1        Seasonal standard deviation [Divisor is (n-1)]
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "util.h"


void *Seasstat(void *argument)
{
  int timestat_date = TIMESTAT_MEAN;
  int gridsize;
  int vdate = 0, vtime = 0;
  int vdate0 = 0, vtime0 = 0;
  int vdate1 = 0, vtime1 = 0;
  int nrecs;
  int varID, levelID, recID;
  long nsets;
  int i;
  int year, month, day, seas, seas0 = 0;
  int nmiss;
  int nlevel;
  int newseas, oldmon = 0, newmon;
  int nseason = 0;
  const char *seas_name[4];

  cdoInitialize(argument);

  cdoOperatorAdd("seasmin",  func_min,  0, NULL);
  cdoOperatorAdd("seasmax",  func_max,  0, NULL);
  cdoOperatorAdd("seassum",  func_sum,  0, NULL);
  cdoOperatorAdd("seasmean", func_mean, 0, NULL);
  cdoOperatorAdd("seasavg",  func_avg,  0, NULL);
  cdoOperatorAdd("seasvar",  func_var,  0, NULL);
  cdoOperatorAdd("seasvar1", func_var1, 0, NULL);
  cdoOperatorAdd("seasstd",  func_std,  0, NULL);
  cdoOperatorAdd("seasstd1", func_std1, 0, NULL);

  int operatorID = cdoOperatorID();
  int operfunc = cdoOperatorF1(operatorID);

  int season_start = get_season_start();
  get_season_name(seas_name);

  int lmean   = operfunc == func_mean || operfunc == func_avg;
  int lstd    = operfunc == func_std || operfunc == func_std1;
  int lvarstd = operfunc == func_std || operfunc == func_var || operfunc == func_std1 || operfunc == func_var1;
  double divisor = operfunc == func_std1 || operfunc == func_var1;

  int streamID1 = streamOpenRead(cdoStreamName(0));

  int vlistID1 = streamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  if ( taxisInqType(taxisID2) == TAXIS_FORECAST ) taxisDefType(taxisID2, TAXIS_RELATIVE);
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

  gridsize = vlistGridsizeMax(vlistID1);

  field_t field;
  field_init(&field);
  field.ptr = (double*) malloc(gridsize*sizeof(double));

  field_t **samp1 = field_malloc(vlistID1, FIELD_NONE);
  field_t **vars1 = field_malloc(vlistID1, FIELD_PTR);
  field_t **vars2 = NULL;
  if ( lvarstd ) vars2 = field_malloc(vlistID1, FIELD_PTR);

  int tsID    = 0;
  int otsID   = 0;
  while ( TRUE )
    {
      nsets = 0;
      newseas = FALSE;
      while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
	{
	  dtlist_taxisInqTimestep(dtlist, taxisID1, nsets);
	  vdate = dtlist_get_vdate(dtlist, nsets);
	  vtime = dtlist_get_vtime(dtlist, nsets);

	  cdiDecodeDate(vdate, &year, &month, &day);

	  newmon = month;

	  if ( season_start == START_DEC && newmon == 12 ) newmon = 0;

          seas = month_to_season(month);

	  if ( nsets == 0 )
	    {
	      nseason++;
	      vdate0 = vdate;
	      vtime0 = vtime;
	      seas0  = seas;
	      oldmon = newmon;
	    }

	  if ( newmon < oldmon ) newseas = TRUE;

	  if ( (seas != seas0) || newseas ) break;

	  oldmon = newmon;

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

	  vdate1 = vdate;
	  vtime1 = vtime;
	  nsets++;
	  tsID++;
	}

      if ( nrecs == 0 && nsets == 0 ) break;

      if ( lmean )
	for ( varID = 0; varID < nvars; varID++ )
	  {
	    if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;
	    nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
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
	    nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
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

      if ( cdoVerbose )
	{
	  char vdatestr0[32], vtimestr0[32];
	  char vdatestr1[32], vtimestr1[32];
	  date2str(vdate0, vdatestr0, sizeof(vdatestr0));
	  time2str(vtime0, vtimestr0, sizeof(vtimestr0));
	  date2str(vdate1, vdatestr1, sizeof(vdatestr1));
	  time2str(vtime1, vtimestr1, sizeof(vtimestr1));
	  cdoPrint("season %3d %3s start %s %s end %s %s ntimesteps %d", 
		   nseason, seas_name[seas0], vdatestr0, vtimestr0, vdatestr1, vtimestr1, nsets);
	}

      dtlist_stat_taxisDefTimestep(dtlist, taxisID2, nsets);
      streamDefTimestep(streamID2, otsID);

      if ( nsets < 3 )
	{
	  char vdatestr[32];
	  date2str(vdate0, vdatestr, sizeof(vdatestr));
	  cdoWarning("Season %3d (%s) has only %d input time step%s!", 
		     otsID+1, vdatestr, nsets, nsets == 1 ? "" : "s");
	}

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
  if ( lvarstd ) field_free(vars2, vlistID1);

  if ( field.ptr ) free(field.ptr);

  if ( recVarID   ) free(recVarID);
  if ( recLevelID ) free(recLevelID);

  dtlist_delete(dtlist);

  streamClose(streamID2);
  streamClose(streamID1);

  cdoFinish();

  return (0);
}
