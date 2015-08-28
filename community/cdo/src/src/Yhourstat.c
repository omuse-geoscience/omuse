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

      Yhourstat   yhourmin         Multi-year hourly minimum
      Yhourstat   yhourmax         Multi-year hourly maximum
      Yhourstat   yhoursum         Multi-year hourly sum
      Yhourstat   yhourmean        Multi-year hourly mean
      Yhourstat   yhouravg         Multi-year hourly average
      Yhourstat   yhourvar         Multi-year hourly variance
      Yhourstat   yhourvar1        Multi-year hourly variance [Divisor is (n-1)]
      Yhourstat   yhourstd         Multi-year hourly standard deviation
      Yhourstat   yhourstd1        Multi-year hourly standard deviation [Divisor is (n-1)]
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


#define  MAX_HOUR  9301  /* 31*12*25 + 1 */

static
int hour_of_year(int vdate, int vtime)
{
  int year, month, day, houroy;
  int hour, minute, second;

  cdiDecodeDate(vdate, &year, &month, &day);
  cdiDecodeTime(vtime, &hour, &minute, &second);
      
  if ( month >= 1 && month <= 12 && day >= 1 && day <=31 && hour >= 0 && hour < 24 )
    houroy = ((month-1)*31 + day - 1)*25 + hour + 1;
  else
    houroy = 0;

  if ( houroy < 0 || houroy >= MAX_HOUR )
    {
      char vdatestr[32], vtimestr[32];
      date2str(vdate, vdatestr, sizeof(vdatestr));
      time2str(vtime, vtimestr, sizeof(vtimestr));
      cdoAbort("Hour of year %d out of range (%s %s)!", houroy, vdatestr, vtimestr);
    }

  return (houroy);
}


void *Yhourstat(void *argument)
{
  int operatorID;
  int operfunc;
  int gridsize;
  int i;
  int varID;
  int recID;
  int vdate, vtime;
  int houroy;
  int nrecs, nrecords;
  int levelID;
  int tsID;
  int otsID;
  long nsets[MAX_HOUR];
  int streamID1, streamID2;
  int vlistID1, vlistID2, taxisID1, taxisID2;
  int nmiss;
  int nvars, nlevel;
  int *recVarID, *recLevelID;
  int vdates[MAX_HOUR], vtimes[MAX_HOUR];
  int lmean = FALSE, lvarstd = FALSE, lstd = FALSE;
  double divisor;
  field_t **vars1[MAX_HOUR], **vars2[MAX_HOUR], **samp1[MAX_HOUR];
  field_t field;

  cdoInitialize(argument);

  cdoOperatorAdd("yhourmin",  func_min,  0, NULL);
  cdoOperatorAdd("yhourmax",  func_max,  0, NULL);
  cdoOperatorAdd("yhoursum",  func_sum,  0, NULL);
  cdoOperatorAdd("yhourmean", func_mean, 0, NULL);
  cdoOperatorAdd("yhouravg",  func_avg,  0, NULL);
  cdoOperatorAdd("yhourvar",  func_var,  0, NULL);
  cdoOperatorAdd("yhourvar1", func_var1, 0, NULL);
  cdoOperatorAdd("yhourstd",  func_std,  0, NULL);
  cdoOperatorAdd("yhourstd1", func_std1, 0, NULL);

  operatorID = cdoOperatorID();
  operfunc = cdoOperatorF1(operatorID);

  lmean   = operfunc == func_mean || operfunc == func_avg;
  lstd    = operfunc == func_std || operfunc == func_std1;
  lvarstd = operfunc == func_std || operfunc == func_var || operfunc == func_std1 || operfunc == func_var1;
  divisor = operfunc == func_std1 || operfunc == func_var1;

  for ( houroy = 0; houroy < MAX_HOUR; ++houroy )
    {
      vars1[houroy] = NULL;
      vars2[houroy] = NULL;
      samp1[houroy] = NULL;
      nsets[houroy] = 0;
    }

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  if ( taxisHasBounds(taxisID2) ) taxisDeleteBounds(taxisID2);
  vlistDefTaxis(vlistID2, taxisID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  nvars    = vlistNvars(vlistID1);
  nrecords = vlistNrecs(vlistID1);

  recVarID   = (int*) malloc(nrecords*sizeof(int));
  recLevelID = (int*) malloc(nrecords*sizeof(int));

  gridsize = vlistGridsizeMax(vlistID1);
  field_init(&field);
  field.ptr = (double*) malloc(gridsize*sizeof(double));

  tsID = 0;
  otsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      vdate = taxisInqVdate(taxisID1);
      vtime = taxisInqVtime(taxisID1);

      if ( cdoVerbose ) cdoPrint("process timestep: %d %d %d", tsID+1, vdate, vtime);

      houroy = hour_of_year(vdate, vtime);

      vdates[houroy] = vdate;
      vtimes[houroy] = vtime;

      if ( vars1[houroy] == NULL )
	{
	  vars1[houroy] = field_malloc(vlistID1, FIELD_PTR);
	  samp1[houroy] = field_malloc(vlistID1, FIELD_NONE);
	  if ( lvarstd )
	    vars2[houroy] = field_malloc(vlistID1, FIELD_PTR);
	}

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);

	  if ( tsID == 0 )
	    {
	      recVarID[recID]   = varID;
	      recLevelID[recID] = levelID;
	    }

	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));

	  if ( nsets[houroy] == 0 )
	    {
	      streamReadRecord(streamID1, vars1[houroy][varID][levelID].ptr, &nmiss);
	      vars1[houroy][varID][levelID].nmiss = nmiss;

	      if ( nmiss > 0 || samp1[houroy][varID][levelID].ptr )
		{
		  if ( samp1[houroy][varID][levelID].ptr == NULL )
		    samp1[houroy][varID][levelID].ptr = (double*) malloc(gridsize*sizeof(double));

		  for ( i = 0; i < gridsize; i++ )
		    if ( DBL_IS_EQUAL(vars1[houroy][varID][levelID].ptr[i],
				      vars1[houroy][varID][levelID].missval) )
		      samp1[houroy][varID][levelID].ptr[i] = 0;
		    else
		      samp1[houroy][varID][levelID].ptr[i] = 1;
		}
	    }
	  else
	    {
	      streamReadRecord(streamID1, field.ptr, &field.nmiss);
	      field.grid    = vars1[houroy][varID][levelID].grid;
	      field.missval = vars1[houroy][varID][levelID].missval;

	      if ( field.nmiss > 0 || samp1[houroy][varID][levelID].ptr )
		{
		  if ( samp1[houroy][varID][levelID].ptr == NULL )
		    {
		      samp1[houroy][varID][levelID].ptr = (double*) malloc(gridsize*sizeof(double));
		      for ( i = 0; i < gridsize; i++ )
			samp1[houroy][varID][levelID].ptr[i] = nsets[houroy];
		    }
		  
		  for ( i = 0; i < gridsize; i++ )
		    if ( !DBL_IS_EQUAL(field.ptr[i], vars1[houroy][varID][levelID].missval) )
		      samp1[houroy][varID][levelID].ptr[i]++;
		}

	      if ( lvarstd )
		{
		  farsumq(&vars2[houroy][varID][levelID], field);
		  farsum(&vars1[houroy][varID][levelID], field);
		}
	      else
		{
		  farfun(&vars1[houroy][varID][levelID], field, operfunc);
		}
	    }
	}

      if ( nsets[houroy] == 0 && lvarstd )
	for ( varID = 0; varID < nvars; varID++ )
	  {
	    if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;
	    gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	    nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	    for ( levelID = 0; levelID < nlevel; levelID++ )
	      farmoq(&vars2[houroy][varID][levelID], vars1[houroy][varID][levelID]);
	  }

      nsets[houroy]++;
      tsID++;
    }

  for ( houroy = 0; houroy < MAX_HOUR; ++houroy )
    if ( nsets[houroy] )
      {
	if ( lmean )
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;
	      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      for ( levelID = 0; levelID < nlevel; levelID++ )
		{
		  if ( samp1[houroy][varID][levelID].ptr == NULL )
		    farcmul(&vars1[houroy][varID][levelID], 1.0/nsets[houroy]);
		  else
		    fardiv(&vars1[houroy][varID][levelID], samp1[houroy][varID][levelID]);
		}
	    }
	else if ( lvarstd )
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;
	      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      for ( levelID = 0; levelID < nlevel; levelID++ )
		{
		  if ( samp1[houroy][varID][levelID].ptr == NULL )
		    {
		      if ( lstd )
			farcstd(&vars1[houroy][varID][levelID], vars2[houroy][varID][levelID], nsets[houroy], divisor);
		      else
			farcvar(&vars1[houroy][varID][levelID], vars2[houroy][varID][levelID], nsets[houroy], divisor);
		    }
		  else
		    {
		      if ( lstd )
			farstd(&vars1[houroy][varID][levelID], vars2[houroy][varID][levelID], samp1[houroy][varID][levelID], divisor);
		      else
			farvar(&vars1[houroy][varID][levelID], vars2[houroy][varID][levelID], samp1[houroy][varID][levelID], divisor);
		    }
		}
	    }

	taxisDefVdate(taxisID2, vdates[houroy]);
	taxisDefVtime(taxisID2, vtimes[houroy]);
	streamDefTimestep(streamID2, otsID);

	for ( recID = 0; recID < nrecords; recID++ )
	  {
	    varID    = recVarID[recID];
	    levelID  = recLevelID[recID];

	    if ( otsID && vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;

	    streamDefRecord(streamID2, varID, levelID);
	    streamWriteRecord(streamID2, vars1[houroy][varID][levelID].ptr,
			      vars1[houroy][varID][levelID].nmiss);
	  }

	otsID++;
      }

  for ( houroy = 0; houroy < MAX_HOUR; ++houroy )
    {
      if ( vars1[houroy] != NULL )
	{
	  field_free(samp1[houroy], vlistID1);
	  field_free(vars1[houroy], vlistID1);
	  if ( lvarstd ) field_free(vars2[houroy], vlistID1);
	}
    }

  if ( field.ptr ) free(field.ptr);

  if ( recVarID   ) free(recVarID);
  if ( recLevelID ) free(recLevelID);

  streamClose(streamID2);
  streamClose(streamID1);

  cdoFinish();

  return (0);
}
