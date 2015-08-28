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

      Ydrunstat    ydrunmin          Multi-year daily running minimum
      Ydrunstat    ydrunmax          Multi-year daily running maximum
      Ydrunstat    ydrunsum          Multi-year daily running sum
      Ydrunstat    ydrunmean         Multi-year daily running mean
      Ydrunstat    ydrunavg          Multi-year daily running average
      Ydrunstat    ydrunvar          Multi-year daily running variance
      Ydrunstat    ydrunvar1         Multi-year daily running variance [Divisor is (n-1)]
      Ydrunstat    ydrunstd          Multi-year daily running standard deviation
      Ydrunstat    ydrunstd1         Multi-year daily running standard deviation [Divisor is (n-1)]
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


#define NDAY 373


typedef struct {
  int       vdate[NDAY];
  int       vtime[NDAY];  
  field_t **vars1[NDAY]; 
  field_t **vars2[NDAY];
  int       nsets[NDAY];
  int       vlist;
}
YDAY_STATS;


static YDAY_STATS *ydstatCreate(int vlistID);
static void ydstatDestroy(YDAY_STATS *stats);
static void ydstatUpdate(YDAY_STATS *stats, int vdate, int vtime, 
  field_t **vars1, field_t **vars2, int nsets, int operfunc);
static void ydstatFinalize(YDAY_STATS *stats, int operfunc);


void *Ydrunstat(void *argument)
{
  int varID;
  int recID;
  int nrecs;
  int levelID;
  int tsID;
  int otsID;
  int inp, its;
  int nmiss;
  int vdate, vtime;
  int dayoy;
    
  cdoInitialize(argument);

  cdoOperatorAdd("ydrunmin",  func_min,  0, NULL);
  cdoOperatorAdd("ydrunmax",  func_max,  0, NULL);
  cdoOperatorAdd("ydrunsum",  func_sum,  0, NULL);
  cdoOperatorAdd("ydrunmean", func_mean, 0, NULL);
  cdoOperatorAdd("ydrunavg",  func_avg,  0, NULL);
  cdoOperatorAdd("ydrunvar",  func_var,  0, NULL);
  cdoOperatorAdd("ydrunvar1", func_var1, 0, NULL);
  cdoOperatorAdd("ydrunstd",  func_std,  0, NULL);
  cdoOperatorAdd("ydrunstd1", func_std1, 0, NULL);

  int operatorID = cdoOperatorID();
  int operfunc = cdoOperatorF1(operatorID);

  operatorInputArg("number of timesteps");
  int ndates = parameter2int(operatorArgv()[0]);

  int lvarstd = operfunc == func_std || operfunc == func_var || operfunc == func_std1 || operfunc == func_var1;
  
  int streamID1 = streamOpenRead(cdoStreamName(0));

  int vlistID1 = streamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  if ( taxisHasBounds(taxisID2) ) taxisDeleteBounds(taxisID2);
  vlistDefTaxis(vlistID2, taxisID2);

  int calendar = taxisInqCalendar(taxisID1);
  int dpy      = calendar_dpy(calendar);

  int streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  int nrecords = vlistNrecs(vlistID1);

  int *recVarID   = (int*) malloc(nrecords*sizeof(int));
  int *recLevelID = (int*) malloc(nrecords*sizeof(int));

  datetime_t *datetime = (datetime_t*) malloc((ndates+1)*sizeof(datetime_t));
  
  YDAY_STATS *stats = ydstatCreate(vlistID1);
  field_t ***vars1 = (field_t ***) malloc((ndates+1)*sizeof(field_t **));
  field_t ***vars2 = NULL;
  if ( lvarstd )
    vars2 = (field_t ***) malloc((ndates+1)*sizeof(field_t **));
  
  for ( its = 0; its < ndates; its++ )
    {
      vars1[its] = field_malloc(vlistID1, FIELD_PTR);
      if ( lvarstd )
	vars2[its] = field_malloc(vlistID1, FIELD_PTR);
    }
  
  for ( tsID = 0; tsID < ndates; tsID++ )
    {
      nrecs = streamInqTimestep(streamID1, tsID);
      if ( nrecs == 0 )
	cdoAbort("File has less then %d timesteps!", ndates);

      datetime[tsID].date = taxisInqVdate(taxisID1);
      datetime[tsID].time = taxisInqVtime(taxisID1);
	
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

	  if ( lvarstd )
	    {
	      farmoq(&vars2[tsID][varID][levelID], vars1[tsID][varID][levelID]);
	      for ( inp = 0; inp < tsID; inp++ )
		{
		  farsumq(&vars2[inp][varID][levelID], vars1[tsID][varID][levelID]);
		  farsum(&vars1[inp][varID][levelID], vars1[tsID][varID][levelID]);
		}
	    }
	  else
	    {
	      for ( inp = 0; inp < tsID; inp++ )
		{
		  farfun(&vars1[inp][varID][levelID], vars1[tsID][varID][levelID], operfunc);
		}
	    }
	}
    }
  
  while ( TRUE )
    {
      datetime_avg(dpy, ndates, datetime);
      
      vdate = datetime[ndates].date;
      vtime = datetime[ndates].time;
      
      if ( lvarstd )   
        ydstatUpdate(stats, vdate, vtime, vars1[0], vars2[0], ndates, operfunc);
      else
        ydstatUpdate(stats, vdate, vtime, vars1[0], NULL, ndates, operfunc);
        
      datetime[ndates] = datetime[0];
      vars1[ndates] = vars1[0];
      if ( lvarstd )
        vars2[ndates] = vars2[0];

      for ( inp = 0; inp < ndates; inp++ )
	{
	  datetime[inp] = datetime[inp+1];
	  vars1[inp] = vars1[inp+1];
	  if ( lvarstd )
	    vars2[inp] = vars2[inp+1];
	}

      nrecs = streamInqTimestep(streamID1, tsID);
      if ( nrecs == 0 ) break;

      datetime[ndates-1].date = taxisInqVdate(taxisID1);
      datetime[ndates-1].time = taxisInqVtime(taxisID1);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  
	  streamReadRecord(streamID1, vars1[ndates-1][varID][levelID].ptr, &nmiss);
	  vars1[ndates-1][varID][levelID].nmiss = nmiss;

	  if ( lvarstd )
	    {
	      for ( inp = 0; inp < ndates-1; inp++ )
		{
		  farsumq(&vars2[inp][varID][levelID], vars1[ndates-1][varID][levelID]);
		  farsum(&vars1[inp][varID][levelID], vars1[ndates-1][varID][levelID]);
		}
	      farmoq(&vars2[ndates-1][varID][levelID], vars1[ndates-1][varID][levelID]);
	    }
	  else
	    {
	      for ( inp = 0; inp < ndates-1; inp++ )
		{
		  farfun(&vars1[inp][varID][levelID], vars1[ndates-1][varID][levelID], operfunc);
		}
	    }
	}

      tsID++;
    }

  /*
  // set the year to the minimum of years found on output timestep
  int outyear = 1e9;
  int year, month, day;
  for ( dayoy = 0; dayoy < NDAY; dayoy++ )
    if ( stats->nsets[dayoy] )
      {
	cdiDecodeDate(stats->vdate[dayoy], &year, &month, &day);
	if ( year < outyear ) outyear = year;
      }
  for ( dayoy = 0; dayoy < NDAY; dayoy++ )
    if ( stats->nsets[dayoy] )
      {
	cdiDecodeDate(stats->vdate[dayoy], &year, &month, &day);
        // printf("vdates[%d] = %d  nsets = %d\n", dayoy, stats->vdate[dayoy], stats->nsets[dayoy]);
	if ( year > outyear ) stats->vdate[dayoy] = cdiEncodeDate(outyear, month, day);
      }
  */
  ydstatFinalize(stats, operfunc);

  otsID = 0;

  for ( dayoy = 0; dayoy < NDAY; dayoy++ )
    if ( stats->nsets[dayoy] )
      {
	taxisDefVdate(taxisID2, stats->vdate[dayoy]);
	taxisDefVtime(taxisID2, stats->vtime[dayoy]);
	streamDefTimestep(streamID2, otsID);

	for ( recID = 0; recID < nrecords; recID++ )
	  {
	    varID   = recVarID[recID];
	    levelID = recLevelID[recID];

	    if ( otsID && vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;

	    streamDefRecord(streamID2, varID, levelID);
	    streamWriteRecord(streamID2, stats->vars1[dayoy][varID][levelID].ptr,
			      stats->vars1[dayoy][varID][levelID].nmiss);
	  }

	otsID++;
      }
  
  for ( its = 0; its < ndates; its++ )
    {
      field_free(vars1[its], vlistID1);
      if ( lvarstd ) field_free(vars2[its], vlistID1);
    }
  
  ydstatDestroy(stats);
  free(vars1);
  if ( lvarstd ) free(vars2);

  if ( datetime ) free(datetime);

  if ( recVarID   ) free(recVarID);
  if ( recLevelID ) free(recLevelID);

  streamClose(streamID2);
  streamClose(streamID1);

  cdoFinish();

  return (0);
}

static
YDAY_STATS *ydstatCreate(int vlistID)
{
  int dayoy;
  
  YDAY_STATS *stats = (YDAY_STATS*) malloc(sizeof(YDAY_STATS));
  
  for ( dayoy = 0; dayoy < NDAY; dayoy++ )
    {
      stats->vdate[dayoy] = 0;
      stats->vtime[dayoy] = 0;
      stats->vars1[dayoy] = NULL;
      stats->vars2[dayoy] = NULL;
      stats->nsets[dayoy] = 0;
    }
  stats->vlist = vlistID;
  
  return stats;
}

static
void ydstatDestroy(YDAY_STATS *stats)
{
  int dayoy, varID, levelID, nvars, nlevels;
  
  if ( stats != NULL )
    {
      nvars = vlistNvars(stats->vlist);
      
      for ( dayoy = 0; dayoy < NDAY; dayoy++ )
        {
          if ( stats->vars1[dayoy] != NULL )
            {
              for ( varID = 0; varID < nvars; varID++ )
                {
              	  nlevels = zaxisInqSize(vlistInqVarZaxis(stats->vlist, varID));
              	  for ( levelID = 0; levelID < nlevels; levelID++ )
              	    free(stats->vars1[dayoy][varID][levelID].ptr);
              	  free(stats->vars1[dayoy][varID]);
                }
              free(stats->vars1[dayoy]);
            }
          if ( stats->vars2[dayoy] != NULL )
            {
              for ( varID = 0; varID < nvars; varID++ )
                {
              	  nlevels = zaxisInqSize(vlistInqVarZaxis(stats->vlist, varID));
              	  for ( levelID = 0; levelID < nlevels; levelID++ )
              	    free(stats->vars2[dayoy][varID][levelID].ptr);
              	  free(stats->vars2[dayoy][varID]);
                }
              free(stats->vars2[dayoy]);
            }
        }
      free(stats);    
    }
}

static
void ydstatUpdate(YDAY_STATS *stats, int vdate, int vtime, 
		  field_t **vars1, field_t **vars2, int nsets, int operfunc)
{
  int varID, levelID, nvars, nlevels;
  int gridsize;
  int year, month, day, dayoy;
  int lvarstd;

  lvarstd = vars2 != NULL;

  nvars = vlistNvars(stats->vlist);

  cdiDecodeDate(vdate, &year, &month, &day);

  if ( month >= 1 && month <= 12 )
    dayoy = (month - 1) * 31 + day;
  else
    dayoy = 0;

  if ( dayoy < 0 || dayoy >= NDAY )
    cdoAbort("day %d out of range!", dayoy);

  stats->vdate[dayoy] = vdate;
  stats->vtime[dayoy] = vtime;

  if ( stats->vars1[dayoy] == NULL )
    {
      stats->vars1[dayoy] = field_malloc(stats->vlist, FIELD_PTR);
      if ( lvarstd )
	stats->vars2[dayoy] = field_malloc(stats->vlist, FIELD_PTR);
    }

  for ( varID = 0; varID  < nvars; varID++ )
    {
      if ( vlistInqVarTsteptype(stats->vlist, varID) == TSTEP_CONSTANT ) continue;
        
      gridsize = gridInqSize(vlistInqVarGrid(stats->vlist, varID));
      nlevels  = zaxisInqSize(vlistInqVarZaxis(stats->vlist, varID));
          
      for ( levelID = 0; levelID < nlevels; levelID++ )
        {
	  if ( stats->nsets[dayoy] == 0 )
	    {
	      memcpy(stats->vars1[dayoy][varID][levelID].ptr, vars1[varID][levelID].ptr, gridsize * sizeof(double));
	      stats->vars1[dayoy][varID][levelID].nmiss = vars1[varID][levelID].nmiss;
	       
	      if ( lvarstd )
	        {
	          memcpy(stats->vars2[dayoy][varID][levelID].ptr, vars2[varID][levelID].ptr, gridsize * sizeof(double));
	          stats->vars2[dayoy][varID][levelID].nmiss = vars2[varID][levelID].nmiss;
	        }
	    }
	  else
	    {
	      if ( lvarstd )
	        {
		  farsum(&stats->vars1[dayoy][varID][levelID], vars1[varID][levelID]);
		  farsum(&stats->vars2[dayoy][varID][levelID], vars2[varID][levelID]);
		}
	      else
		{
	          farfun(&stats->vars1[dayoy][varID][levelID], vars1[varID][levelID], operfunc);
		}
	    }
        }
    }

  stats->nsets[dayoy] += nsets;
}

static
void ydstatFinalize(YDAY_STATS *stats, int operfunc)
{
  int varID, levelID, nvars, nlevels;
  int dayoy;
  double divisor;

  divisor = operfunc == func_std1 || operfunc == func_var1;

  nvars = vlistNvars(stats->vlist);
  
  for ( dayoy = 0; dayoy < NDAY; dayoy++ )
    if ( stats->nsets[dayoy] )
      {
      	switch ( operfunc )
      	  {
	    case func_avg:
	    case func_mean:
	      for ( varID = 0; varID < nvars; varID++ )
	        {
	          if ( vlistInqVarTsteptype(stats->vlist, varID) == TSTEP_CONSTANT ) continue;
	          nlevels = zaxisInqSize(vlistInqVarZaxis(stats->vlist, varID));
	          for ( levelID = 0; levelID < nlevels; levelID++ )
		    farcmul(&stats->vars1[dayoy][varID][levelID], 1.0 / stats->nsets[dayoy]);
	        }
	      break;
	      
	    case func_std:
	    case func_std1:
	      for ( varID = 0; varID < nvars; varID++ )
	        {
	          if ( vlistInqVarTsteptype(stats->vlist, varID) == TSTEP_CONSTANT ) continue;
	          nlevels = zaxisInqSize(vlistInqVarZaxis(stats->vlist, varID));
	          for ( levelID = 0; levelID < nlevels; levelID++ )
		    farcstd(&stats->vars1[dayoy][varID][levelID], stats->vars2[dayoy][varID][levelID],
                            stats->nsets[dayoy], divisor);
	        }
	      break;
	      
	    case func_var:
	    case func_var1:
	      for ( varID = 0; varID < nvars; varID++ )
	        {
	          if ( vlistInqVarTsteptype(stats->vlist, varID) == TSTEP_CONSTANT ) continue;
	          nlevels = zaxisInqSize(vlistInqVarZaxis(stats->vlist, varID));
	          for ( levelID = 0; levelID < nlevels; levelID++ )
		    farcvar(&stats->vars1[dayoy][varID][levelID], stats->vars2[dayoy][varID][levelID],
			    stats->nsets[dayoy], divisor);
	        }
	      break;
      	  }
      }
}
