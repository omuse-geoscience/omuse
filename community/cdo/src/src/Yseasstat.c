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

      Yseasstat  yseasmin        Multi-year seasonally minimum
      Yseasstat  yseasmax        Multi-year seasonally maximum
      Yseasstat  yseassum        Multi-year seasonally sum
      Yseasstat  yseasmean       Multi-year seasonally mean
      Yseasstat  yseasavg        Multi-year seasonally average
      Yseasstat  yseasvar        Multi-year seasonally variance
      Yseasstat  yseasvar1       Multi-year seasonally variance [Divisor is (n-1)]
      Yseasstat  yseasstd        Multi-year seasonally standard deviation
      Yseasstat  yseasstd1       Multi-year seasonally standard deviation [Divisor is (n-1)]
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "util.h"


#define  NSEAS       4

typedef struct {
  int vdate;
  int vtime;
}
date_time_t;

void set_date(int vdate_new, int vtime_new, date_time_t *datetime)
{
  int year, month, day;

  cdiDecodeDate(vdate_new, &year, &month, &day);
  if ( month == 12 ) vdate_new = cdiEncodeDate(year-1, month, day);

  if ( vdate_new > datetime->vdate )
    {
      datetime->vdate = vdate_new;
      datetime->vtime = vtime_new;
    }
}


void *Yseasstat(void *argument)
{
  int i;
  int varID;
  int recID;
  int vdate, vtime;
  int year, month, day, seas;
  int nrecs;
  int levelID;
  long nsets[NSEAS];
  int nmiss;
  int nlevel;
  date_time_t datetime[NSEAS];
  field_t **vars1[NSEAS], **vars2[NSEAS], **samp1[NSEAS];

  cdoInitialize(argument);

  cdoOperatorAdd("yseasmin",  func_min,  0, NULL);
  cdoOperatorAdd("yseasmax",  func_max,  0, NULL);
  cdoOperatorAdd("yseassum",  func_sum,  0, NULL);
  cdoOperatorAdd("yseasmean", func_mean, 0, NULL);
  cdoOperatorAdd("yseasavg",  func_avg,  0, NULL);
  cdoOperatorAdd("yseasvar",  func_var,  0, NULL);
  cdoOperatorAdd("yseasvar1", func_var1, 0, NULL);
  cdoOperatorAdd("yseasstd",  func_std,  0, NULL);
  cdoOperatorAdd("yseasstd1", func_std1, 0, NULL);

  int operatorID = cdoOperatorID();
  int operfunc = cdoOperatorF1(operatorID);

  for ( seas = 0; seas < NSEAS; seas++ )
    {
      vars1[seas]  = NULL;
      vars2[seas]  = NULL;
      samp1[seas]  = NULL;
      nsets[seas]  = 0;
      datetime[seas].vdate = 0;
      datetime[seas].vtime = 0;
    }

  int lmean   = operfunc == func_mean || operfunc == func_avg;
  int lstd    = operfunc == func_std || operfunc == func_std1;
  int lvarstd = operfunc == func_std || operfunc == func_var || operfunc == func_std1 || operfunc == func_var1;
  double divisor = operfunc == func_std1 || operfunc == func_var1;

  int streamID1 = streamOpenRead(cdoStreamName(0));

  int vlistID1 = streamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  if ( taxisHasBounds(taxisID2) ) taxisDeleteBounds(taxisID2);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  int nvars    = vlistNvars(vlistID1);
  int nrecords = vlistNrecs(vlistID1);

  int *recVarID   = (int*) malloc(nrecords*sizeof(int));
  int *recLevelID = (int*) malloc(nrecords*sizeof(int));

  int gridsize = vlistGridsizeMax(vlistID1);

  field_t field;
  field_init(&field);
  field.ptr = (double*) malloc(gridsize*sizeof(double));

  int tsID = 0;
  int otsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      vdate = taxisInqVdate(taxisID1);
      vtime = taxisInqVtime(taxisID1);
      cdiDecodeDate(vdate, &year, &month, &day);

      seas = month_to_season(month); 

      set_date(vdate, vtime, &datetime[seas]);

      if ( vars1[seas] == NULL )
	{
	  vars1[seas] = field_malloc(vlistID1, FIELD_PTR);
	  samp1[seas] = field_malloc(vlistID1, FIELD_NONE);
	  if ( lvarstd )
	    vars2[seas] = field_malloc(vlistID1, FIELD_PTR);
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

	  if ( nsets[seas] == 0 )
	    {
	      streamReadRecord(streamID1, vars1[seas][varID][levelID].ptr, &nmiss);
	      vars1[seas][varID][levelID].nmiss = nmiss;

	      if ( nmiss > 0 || samp1[seas][varID][levelID].ptr )
		{
		  if ( samp1[seas][varID][levelID].ptr == NULL )
		    samp1[seas][varID][levelID].ptr = (double*) malloc(gridsize*sizeof(double));

		  for ( i = 0; i < gridsize; i++ )
		    if ( DBL_IS_EQUAL(vars1[seas][varID][levelID].ptr[i],
				      vars1[seas][varID][levelID].missval) )
		      samp1[seas][varID][levelID].ptr[i] = 0;
		    else
		      samp1[seas][varID][levelID].ptr[i] = 1;
		}
	    }
	  else
	    {
	      streamReadRecord(streamID1, field.ptr, &field.nmiss);
	      field.grid    = vars1[seas][varID][levelID].grid;
	      field.missval = vars1[seas][varID][levelID].missval;

	      if ( field.nmiss > 0 || samp1[seas][varID][levelID].ptr )
		{
		  if ( samp1[seas][varID][levelID].ptr == NULL )
		    {
		      samp1[seas][varID][levelID].ptr = (double*) malloc(gridsize*sizeof(double));
		      for ( i = 0; i < gridsize; i++ )
			samp1[seas][varID][levelID].ptr[i] = nsets[seas];
		    }
		  
		  for ( i = 0; i < gridsize; i++ )
		    if ( !DBL_IS_EQUAL(field.ptr[i], vars1[seas][varID][levelID].missval) )
		      samp1[seas][varID][levelID].ptr[i]++;
		}

	      if ( lvarstd )
		{
		  farsumq(&vars2[seas][varID][levelID], field);
		  farsum(&vars1[seas][varID][levelID], field);
		}
	      else
		{
		  farfun(&vars1[seas][varID][levelID], field, operfunc);
		}
	    }
	}

      if ( nsets[seas] == 0 && lvarstd )
	for ( varID = 0; varID < nvars; varID++ )
	  {
	    if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;
	    nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	    for ( levelID = 0; levelID < nlevel; levelID++ )
	      farmoq(&vars2[seas][varID][levelID], vars1[seas][varID][levelID]);
	  }

      nsets[seas]++;
      tsID++;
    }

  for ( seas = 0; seas < NSEAS; seas++ )
    if ( nsets[seas] )
      {
	if ( lmean )
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;
	      nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      for ( levelID = 0; levelID < nlevel; levelID++ )
		{
		  if ( samp1[seas][varID][levelID].ptr == NULL )
		    farcmul(&vars1[seas][varID][levelID], 1.0/nsets[seas]);
		  else
		    fardiv(&vars1[seas][varID][levelID], samp1[seas][varID][levelID]);
		}
	    }
	else if ( lvarstd )
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;
	      nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      for ( levelID = 0; levelID < nlevel; levelID++ )
		{
		  if ( samp1[seas][varID][levelID].ptr == NULL )
		    {
		      if ( lstd )
			farcstd(&vars1[seas][varID][levelID], vars2[seas][varID][levelID], nsets[seas], divisor);
		      else
			farcvar(&vars1[seas][varID][levelID], vars2[seas][varID][levelID], nsets[seas], divisor);
		    }
		  else
		    {
		      if ( lstd )
			farstd(&vars1[seas][varID][levelID], vars2[seas][varID][levelID], samp1[seas][varID][levelID], divisor);
		      else
			farvar(&vars1[seas][varID][levelID], vars2[seas][varID][levelID], samp1[seas][varID][levelID], divisor);
		    }
		}
	    }

	taxisDefVdate(taxisID2, datetime[seas].vdate);
	taxisDefVtime(taxisID2, datetime[seas].vtime);
	streamDefTimestep(streamID2, otsID);

	for ( recID = 0; recID < nrecords; recID++ )
	  {
	    varID    = recVarID[recID];
	    levelID  = recLevelID[recID];

	    if ( otsID && vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;

	    streamDefRecord(streamID2, varID, levelID);
	    streamWriteRecord(streamID2, vars1[seas][varID][levelID].ptr,
			      vars1[seas][varID][levelID].nmiss);
	  }

	otsID++;
      }

  for ( seas = 0; seas < NSEAS; seas++ )
    {
      if ( vars1[seas] != NULL )
	{
	  field_free(vars1[seas], vlistID1);
	  field_free(samp1[seas], vlistID1);
	  if ( lvarstd ) free(vars2[seas]);
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
