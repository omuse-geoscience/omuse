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

      Ymonpctl   ymonpctl        Multi-year monthly percentiles
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "percentiles.h"

#define  NMONTH     17

static
int getmonth(int date)
{
  int year, month, day;
  cdiDecodeDate(date, &year, &month, &day);
  return month;
}

void *Ymonpctl(void *argument)
{
  int gridsize;
  int varID;
  int recID;
  int gridID;
  int vdate, vtime;
  int year, month, day;
  int nrecs, nrecords;
  int levelID;
  int tsID;
  int otsID;
  long nsets[NMONTH];
  int streamID1, streamID2, streamID3, streamID4;
  int vlistID1, vlistID2, vlistID3, vlistID4, taxisID1, taxisID2, taxisID3, taxisID4;
  int nmiss;
  int nvars, nlevels;
  int *recVarID, *recLevelID;
  int vdates1[NMONTH], vtimes1[NMONTH];
  int vdates2[NMONTH];
  field_t **vars1[NMONTH];
  field_t field;
  double pn;
  HISTOGRAM_SET *hsets[NMONTH];

  cdoInitialize(argument);
  cdoOperatorAdd("ymonpctl", func_pctl, 0, NULL);

  operatorInputArg("percentile number");
  pn = parameter2double(operatorArgv()[0]);
      
  if ( !(pn > 0 && pn < 100) )
    cdoAbort("Illegal argument: percentile number %g is not in the range 0..100!", pn);

  for ( month = 0; month < NMONTH; month++ )
    {
      vars1[month] = NULL;
      hsets[month] = NULL;
      nsets[month] = 0;
    }

  streamID1 = streamOpenRead(cdoStreamName(0));
  streamID2 = streamOpenRead(cdoStreamName(1));
  streamID3 = streamOpenRead(cdoStreamName(2));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = streamInqVlist(streamID2);
  vlistID3 = streamInqVlist(streamID3);
  vlistID4 = vlistDuplicate(vlistID1);

  vlistCompare(vlistID1, vlistID2, CMP_ALL);
  vlistCompare(vlistID1, vlistID3, CMP_ALL);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = vlistInqTaxis(vlistID2);
  taxisID3 = vlistInqTaxis(vlistID3);
  /* TODO - check that time axes 2 and 3 are equal */

  taxisID4 = taxisDuplicate(taxisID1);
  if ( taxisHasBounds(taxisID4) ) taxisDeleteBounds(taxisID4);
  vlistDefTaxis(vlistID4, taxisID4);

  streamID4 = streamOpenWrite(cdoStreamName(3), cdoFiletype());

  streamDefVlist(streamID4, vlistID4);

  nvars    = vlistNvars(vlistID1);
  nrecords = vlistNrecs(vlistID1);

  recVarID   = (int*) malloc(nrecords*sizeof(int));
  recLevelID = (int*) malloc(nrecords*sizeof(int));

  gridsize = vlistGridsizeMax(vlistID1);
  field_init(&field);
  field.ptr = (double*) malloc(gridsize*sizeof(double));

  tsID = 0;
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
      if ( month < 0 || month >= NMONTH )
	cdoAbort("Month %d out of range!", month);

      vdates2[month] = vdate;

      if ( vars1[month] == NULL )
	{
	  vars1[month] = field_malloc(vlistID1, FIELD_PTR);
          hsets[month] = hsetCreate(nvars);

	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      gridID   = vlistInqVarGrid(vlistID1, varID);
	      nlevels  = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));

              hsetCreateVarLevels(hsets[month], varID, nlevels, gridID);
	    }
	}
      
      for ( recID = 0; recID < nrecs; recID++ )
        {
          streamInqRecord(streamID2, &varID, &levelID);
	  streamReadRecord(streamID2, vars1[month][varID][levelID].ptr, &nmiss);
          vars1[month][varID][levelID].nmiss = nmiss;
        }
      for ( recID = 0; recID < nrecs; recID++ )
        {
          streamInqRecord(streamID3, &varID, &levelID);
	  streamReadRecord(streamID3, field.ptr, &nmiss);
          field.nmiss   = nmiss;
          field.grid    = vars1[month][varID][levelID].grid;
	  field.missval = vars1[month][varID][levelID].missval;
	  
	  hsetDefVarLevelBounds(hsets[month], varID, levelID, &vars1[month][varID][levelID], &field);
        }
      
      tsID++;
    }

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      vdate = taxisInqVdate(taxisID1);
      vtime = taxisInqVtime(taxisID1);

      if ( cdoVerbose ) cdoPrint("process timestep: %d %d %d", tsID+1, vdate, vtime);

      cdiDecodeDate(vdate, &year, &month, &day);
      if ( month < 0 || month >= NMONTH )
	cdoAbort("Month %d out of range!", month);

      vdates1[month] = vdate;
      vtimes1[month] = vtime;

      if ( vars1[month] == NULL )
        cdoAbort("No data for month %d in %s and %s", month, cdoStreamName(1)->args, cdoStreamName(2)->args);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);

	  if ( tsID == 0 )
	    {
	      recVarID[recID]   = varID;
	      recLevelID[recID] = levelID;
	    }

	  streamReadRecord(streamID1, vars1[month][varID][levelID].ptr, &nmiss);
	  vars1[month][varID][levelID].nmiss = nmiss;

	  hsetAddVarLevelValues(hsets[month], varID, levelID, &vars1[month][varID][levelID]);
	}

      nsets[month]++;
      tsID++;
    }

  otsID = 0;
  for ( month = 0; month < NMONTH; month++ )
    if ( nsets[month] )
      {
        if ( getmonth(vdates1[month]) != getmonth(vdates2[month]) )
          cdoAbort("Verification dates for month %d of %s, %s and %s are different!",
                   month, cdoStreamName(1)->args, cdoStreamName(2)->args, cdoStreamName(3)->args);

	for ( varID = 0; varID < nvars; varID++ )
	  {
	    if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;
	    nlevels = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      
	    for ( levelID = 0; levelID < nlevels; levelID++ )
	      hsetGetVarLevelPercentiles(&vars1[month][varID][levelID], hsets[month], varID, levelID, pn);
	  }

	taxisDefVdate(taxisID4, vdates1[month]);
	taxisDefVtime(taxisID4, vtimes1[month]);
	streamDefTimestep(streamID4, otsID);

	for ( recID = 0; recID < nrecords; recID++ )
	  {
	    varID    = recVarID[recID];
	    levelID  = recLevelID[recID];

	    if ( otsID && vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;

	    streamDefRecord(streamID4, varID, levelID);
	    streamWriteRecord(streamID4, vars1[month][varID][levelID].ptr, vars1[month][varID][levelID].nmiss);
	  }

	otsID++;
      }

  for ( month = 0; month < NMONTH; month++ )
    {
      if ( vars1[month] != NULL )
	{
	  field_free(vars1[month], vlistID1); 
	  hsetDestroy(hsets[month]);
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
