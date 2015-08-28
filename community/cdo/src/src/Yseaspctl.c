/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2006 Brockmann Consult
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

      Yseaspctl  yseaspctl       Multi-year seasonally percentiles
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "percentiles.h"
#include "util.h"

#define  NSEAS       4

typedef struct {
  int vdate;
  int vtime;
}
date_time_t;

void set_date(int vdate_new, int vtime_new, date_time_t *datetime);

int getmonthday(int date);

void *Yseaspctl(void *argument)
{
  int varID;
  int recID;
  int gridID;
  int vdate, vtime;
  int year, month, day, seas;
  int nrecs;
  int levelID;
  long nsets[NSEAS];
  int nmiss;
  int nlevels;
  date_time_t datetime1[NSEAS], datetime2[NSEAS];
  field_t **vars1[NSEAS];
  HISTOGRAM_SET *hsets[NSEAS];

  cdoInitialize(argument);
  cdoOperatorAdd("yseaspctl", func_pctl, 0, NULL);

  operatorInputArg("percentile number");
  double pn = parameter2double(operatorArgv()[0]);
      
  if ( !(pn > 0 && pn < 100) )
    cdoAbort("Illegal argument: percentile number %g is not in the range 0..100!", pn);

  for ( seas = 0; seas < NSEAS; seas++ )
    {
      vars1[seas] = NULL;
      hsets[seas] = NULL;
      nsets[seas] = 0;
      datetime1[seas].vdate = 0;
      datetime1[seas].vtime = 0;
      datetime2[seas].vdate = 0;
      datetime2[seas].vtime = 0;
    }

  int streamID1 = streamOpenRead(cdoStreamName(0));
  int streamID2 = streamOpenRead(cdoStreamName(1));
  int streamID3 = streamOpenRead(cdoStreamName(2));

  int vlistID1 = streamInqVlist(streamID1);
  int vlistID2 = streamInqVlist(streamID2);
  int vlistID3 = streamInqVlist(streamID3);
  int vlistID4 = vlistDuplicate(vlistID1);

  vlistCompare(vlistID1, vlistID2, CMP_ALL);
  vlistCompare(vlistID1, vlistID3, CMP_ALL);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = vlistInqTaxis(vlistID2);
  int taxisID3 = vlistInqTaxis(vlistID3);
  /* TODO - check that time axes 2 and 3 are equal */

  int taxisID4 = taxisDuplicate(taxisID1);
  if ( taxisHasBounds(taxisID4) ) taxisDeleteBounds(taxisID4);
  vlistDefTaxis(vlistID4, taxisID4);

  int streamID4 = streamOpenWrite(cdoStreamName(3), cdoFiletype());

  streamDefVlist(streamID4, vlistID4);

  int nvars    = vlistNvars(vlistID1);
  int nrecords = vlistNrecs(vlistID1);

  int *recVarID   = (int*) malloc(nrecords*sizeof(int));
  int *recLevelID = (int*) malloc(nrecords*sizeof(int));

  int gridsize = vlistGridsizeMax(vlistID1);

  field_t field;
  field_init(&field);
  field.ptr = (double*) malloc(gridsize*sizeof(double));

  int tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID2, tsID)) )
    {
      if ( nrecs != streamInqTimestep(streamID3, tsID) )
        cdoAbort("Number of records at time step %d of %s and %s differ!", tsID+1, cdoStreamName(1)->args, cdoStreamName(2)->args);
      
      vdate = taxisInqVdate(taxisID2);
      vtime = taxisInqVtime(taxisID2);
      
      if ( vdate != taxisInqVdate(taxisID3) || vtime != taxisInqVtime(taxisID3) )
        cdoAbort("Verification dates at time step %d of %s and %s differ!", tsID+1, cdoStreamName(1)->args, cdoStreamName(2)->args);
        
      if ( cdoVerbose ) cdoPrint("process timestep: %d %d %d", tsID+1, vdate, vtime);

      cdiDecodeDate(vdate, &year, &month, &day);

      seas = month_to_season(month);

      set_date(vdate, vtime, &datetime2[seas]);

      if ( vars1[seas] == NULL )
	{
	  vars1[seas] = field_malloc(vlistID1, FIELD_PTR);
          hsets[seas] = hsetCreate(nvars);

	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      gridID   = vlistInqVarGrid(vlistID1, varID);
	      nlevels  = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));

              hsetCreateVarLevels(hsets[seas], varID, nlevels, gridID);
	    }
	}
      
      for ( recID = 0; recID < nrecs; recID++ )
        {
          streamInqRecord(streamID2, &varID, &levelID);
	  streamReadRecord(streamID2, vars1[seas][varID][levelID].ptr, &nmiss);
          vars1[seas][varID][levelID].nmiss = nmiss;
        }
      for ( recID = 0; recID < nrecs; recID++ )
        {
          streamInqRecord(streamID3, &varID, &levelID);
	  streamReadRecord(streamID3, field.ptr, &nmiss);
          field.nmiss   = nmiss;
          field.grid    = vars1[seas][varID][levelID].grid;
	  field.missval = vars1[seas][varID][levelID].missval;
	  
	  hsetDefVarLevelBounds(hsets[seas], varID, levelID, &vars1[seas][varID][levelID], &field);
        }
      
      tsID++;
    }

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      vdate = taxisInqVdate(taxisID1);
      vtime = taxisInqVtime(taxisID1);
      
      cdiDecodeDate(vdate, &year, &month, &day);
      if ( month < 0 || month > 16 )
	cdoAbort("Month %d out of range!", month);

      if ( month <= 12 )
	seas = (month % 12) / 3;
      else
	seas = month - 13;
      if ( seas < 0 || seas > 3 )
	cdoAbort("Season %d out of range!", seas+1);

      set_date(vdate, vtime, &datetime1[seas]);

      if ( vars1[seas] == NULL )
        cdoAbort("No data for season %d in %s and %s", seas, cdoStreamName(1)->args, cdoStreamName(2)->args);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);

	  if ( tsID == 0 )
	    {
	      recVarID[recID]   = varID;
	      recLevelID[recID] = levelID;
	    }

	  streamReadRecord(streamID1, vars1[seas][varID][levelID].ptr, &nmiss);
	  vars1[seas][varID][levelID].nmiss = nmiss;
	  
	  hsetAddVarLevelValues(hsets[seas], varID, levelID, &vars1[seas][varID][levelID]);
	}

      nsets[seas]++;
      tsID++;
    }

  int otsID = 0;
  for ( seas = 0; seas < NSEAS; seas++ )
    if ( nsets[seas] )
      {
        if ( getmonthday(datetime1[seas].vdate) != getmonthday(datetime2[seas].vdate) )
          cdoAbort("Verification dates for season %d of %s, %s and %s are different!",
                   seas, cdoStreamName(1)->args, cdoStreamName(2)->args, cdoStreamName(3)->args);

	for ( varID = 0; varID < nvars; varID++ )
	  {
	    if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;
	    nlevels = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	    
	    for ( levelID = 0; levelID < nlevels; levelID++ )
	      hsetGetVarLevelPercentiles(&vars1[seas][varID][levelID], hsets[seas], varID, levelID, pn);
	  }

	taxisDefVdate(taxisID4, datetime1[seas].vdate);
	taxisDefVtime(taxisID4, datetime1[seas].vtime);
	streamDefTimestep(streamID4, otsID);

	for ( recID = 0; recID < nrecords; recID++ )
	  {
	    varID    = recVarID[recID];
	    levelID  = recLevelID[recID];

	    if ( otsID && vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;

	    streamDefRecord(streamID4, varID, levelID);
	    streamWriteRecord(streamID4, vars1[seas][varID][levelID].ptr, vars1[seas][varID][levelID].nmiss);
	  }

	otsID++;
      }

  for ( seas = 0; seas < NSEAS; seas++ )
    {
      if ( vars1[seas] != NULL )
	{
	  field_free(vars1[seas], vlistID1); 
	  hsetDestroy(hsets[seas]);
	}
    }

  if ( field.ptr ) free(field.ptr);

  if ( recVarID   ) free(recVarID);
  if ( recLevelID ) free(recLevelID);

  streamClose(streamID4);
  streamClose(streamID3);
  streamClose(streamID2);
  streamClose(streamID1);

  cdoFinish();

  return (0);
}
