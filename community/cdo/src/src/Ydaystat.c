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

      Ydaystat   ydaymin         Multi-year daily minimum
      Ydaystat   ydaymax         Multi-year daily maximum
      Ydaystat   ydaysum         Multi-year daily sum
      Ydaystat   ydaymean        Multi-year daily mean
      Ydaystat   ydayavg         Multi-year daily average
      Ydaystat   ydayvar         Multi-year daily variance
      Ydaystat   ydayvar1        Multi-year daily variance [Divisor is (n-1)]
      Ydaystat   ydaystd         Multi-year daily standard deviation
      Ydaystat   ydaystd1        Multi-year daily standard deviation [Divisor is (n-1)]
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


#define  NDAY       373


void *Ydaystat(void *argument)
{
  int i;
  int varID;
  int recID;
  int vdate, vtime;
  int year, month, day, dayoy;
  int nrecs;
  int levelID;
  long nsets[NDAY];
  int nmiss;
  int nlevel;
  int vdates[NDAY], vtimes[NDAY];
  field_t **vars1[NDAY], **vars2[NDAY], **samp1[NDAY];

  cdoInitialize(argument);

  cdoOperatorAdd("ydaymin",  func_min,  0, NULL);
  cdoOperatorAdd("ydaymax",  func_max,  0, NULL);
  cdoOperatorAdd("ydaysum",  func_sum,  0, NULL);
  cdoOperatorAdd("ydaymean", func_mean, 0, NULL);
  cdoOperatorAdd("ydayavg",  func_avg,  0, NULL);
  cdoOperatorAdd("ydayvar",  func_var,  0, NULL);
  cdoOperatorAdd("ydayvar1", func_var1, 0, NULL);
  cdoOperatorAdd("ydaystd",  func_std,  0, NULL);
  cdoOperatorAdd("ydaystd1", func_std1, 0, NULL);

  int operatorID = cdoOperatorID();
  int operfunc = cdoOperatorF1(operatorID);

  int lmean   = operfunc == func_mean || operfunc == func_avg;
  int lstd    = operfunc == func_std || operfunc == func_std1;
  int lvarstd = operfunc == func_std || operfunc == func_var || operfunc == func_std1 || operfunc == func_var1;
  double divisor = operfunc == func_std1 || operfunc == func_var1;

  for ( dayoy = 0; dayoy < NDAY; dayoy++ )
    {
      vars1[dayoy] = NULL;
      vars2[dayoy] = NULL;
      samp1[dayoy] = NULL;
      nsets[dayoy] = 0;
    }

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

      if ( cdoVerbose ) cdoPrint("process timestep: %d %d %d", tsID+1, vdate, vtime);

      cdiDecodeDate(vdate, &year, &month, &day);

      if ( month >= 1 && month <= 12 )
	dayoy = (month-1)*31 + day;
      else
	dayoy = 0;

      if ( dayoy < 0 || dayoy >= NDAY )
	cdoAbort("day of year %d out of range (date=%d)!", dayoy, vdate);

      vdates[dayoy] = vdate;
      vtimes[dayoy] = vtime;

      if ( vars1[dayoy] == NULL )
	{
	  vars1[dayoy] = field_malloc(vlistID1, FIELD_PTR);
	  samp1[dayoy] = field_malloc(vlistID1, FIELD_NONE);
	  if ( lvarstd )
	    vars2[dayoy] = field_malloc(vlistID1, FIELD_PTR);
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

	  if ( nsets[dayoy] == 0 )
	    {
	      streamReadRecord(streamID1, vars1[dayoy][varID][levelID].ptr, &nmiss);
	      vars1[dayoy][varID][levelID].nmiss = nmiss;

	      if ( nmiss > 0 || samp1[dayoy][varID][levelID].ptr )
		{
		  if ( samp1[dayoy][varID][levelID].ptr == NULL )
		    samp1[dayoy][varID][levelID].ptr = (double*) malloc(gridsize*sizeof(double));

		  for ( i = 0; i < gridsize; i++ )
		    if ( DBL_IS_EQUAL(vars1[dayoy][varID][levelID].ptr[i],
				      vars1[dayoy][varID][levelID].missval) )
		      samp1[dayoy][varID][levelID].ptr[i] = 0;
		    else
		      samp1[dayoy][varID][levelID].ptr[i] = 1;
		}
	    }
	  else
	    {
	      streamReadRecord(streamID1, field.ptr, &field.nmiss);
	      field.grid    = vars1[dayoy][varID][levelID].grid;
	      field.missval = vars1[dayoy][varID][levelID].missval;

	      if ( field.nmiss > 0 || samp1[dayoy][varID][levelID].ptr )
		{
		  if ( samp1[dayoy][varID][levelID].ptr == NULL )
		    {
		      samp1[dayoy][varID][levelID].ptr = (double*) malloc(gridsize*sizeof(double));
		      for ( i = 0; i < gridsize; i++ )
			samp1[dayoy][varID][levelID].ptr[i] = nsets[dayoy];
		    }
		  
		  for ( i = 0; i < gridsize; i++ )
		    if ( !DBL_IS_EQUAL(field.ptr[i], vars1[dayoy][varID][levelID].missval) )
		      samp1[dayoy][varID][levelID].ptr[i]++;
		}

	      if ( lvarstd )
		{
		  farsumq(&vars2[dayoy][varID][levelID], field);
		  farsum(&vars1[dayoy][varID][levelID], field);
		}
	      else
		{
		  farfun(&vars1[dayoy][varID][levelID], field, operfunc);
		}
	    }
	}

      if ( nsets[dayoy] == 0 && lvarstd )
	for ( varID = 0; varID < nvars; varID++ )
	  {
	    if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;
	    gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	    nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	    for ( levelID = 0; levelID < nlevel; levelID++ )
	      farmoq(&vars2[dayoy][varID][levelID], vars1[dayoy][varID][levelID]);
	  }

      nsets[dayoy]++;
      tsID++;
    }

  // set the year to the minimum of years found on output timestep
  int outyear = 1e9;
  for ( dayoy = 0; dayoy < NDAY; dayoy++ )
    if ( nsets[dayoy] )
      {
        cdiDecodeDate(vdates[dayoy], &year, &month, &day);
        if ( year < outyear ) outyear = year;
      }
  for ( dayoy = 0; dayoy < NDAY; dayoy++ )
    if ( nsets[dayoy] )
      {
        cdiDecodeDate(vdates[dayoy], &year, &month, &day);
        if ( year > outyear ) vdates[dayoy] = cdiEncodeDate(outyear, month, day);
        //  printf("vdates[%d] = %d  nsets = %d\n", dayoy, vdates[dayoy], nsets[dayoy]);
      }

  for ( dayoy = 0; dayoy < NDAY; dayoy++ )
    if ( nsets[dayoy] )
      {
	if ( lmean )
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;
	      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      for ( levelID = 0; levelID < nlevel; levelID++ )
		{
		  if ( samp1[dayoy][varID][levelID].ptr == NULL )
		    farcmul(&vars1[dayoy][varID][levelID], 1.0/nsets[dayoy]);
		  else
		    fardiv(&vars1[dayoy][varID][levelID], samp1[dayoy][varID][levelID]);
		}
	    }
	else if ( lvarstd )
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;
	      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      for ( levelID = 0; levelID < nlevel; levelID++ )
		{
		  if ( samp1[dayoy][varID][levelID].ptr == NULL )
		    {
		      if ( lstd )
			farcstd(&vars1[dayoy][varID][levelID], vars2[dayoy][varID][levelID], nsets[dayoy], divisor);
		      else
			farcvar(&vars1[dayoy][varID][levelID], vars2[dayoy][varID][levelID], nsets[dayoy], divisor);
		    }
		  else
		    {
		      if ( lstd )
			farstd(&vars1[dayoy][varID][levelID], vars2[dayoy][varID][levelID], samp1[dayoy][varID][levelID], divisor);
		      else
			farvar(&vars1[dayoy][varID][levelID], vars2[dayoy][varID][levelID], samp1[dayoy][varID][levelID], divisor);
		    }
		}
	    }

	taxisDefVdate(taxisID2, vdates[dayoy]);
	taxisDefVtime(taxisID2, vtimes[dayoy]);
	streamDefTimestep(streamID2, otsID);

	for ( recID = 0; recID < nrecords; recID++ )
	  {
	    varID   = recVarID[recID];
	    levelID = recLevelID[recID];

	    if ( otsID && vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;

	    streamDefRecord(streamID2, varID, levelID);
	    streamWriteRecord(streamID2, vars1[dayoy][varID][levelID].ptr,
			      vars1[dayoy][varID][levelID].nmiss);
	  }

	otsID++;
      }

  for ( dayoy = 0; dayoy < NDAY; dayoy++ )
    {
      if ( vars1[dayoy] != NULL )
	{
	  field_free(vars1[dayoy], vlistID1);
	  field_free(samp1[dayoy], vlistID1);
	  if ( lvarstd ) field_free(vars2[dayoy], vlistID1);
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
