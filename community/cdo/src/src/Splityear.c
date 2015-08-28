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

     Splityear  splityear       Split in years
     Splityear  splityearmon    Split in years and month
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


#define MAX_YEARS 99999

void *Splityear(void *argument)
{
  int streamID2 = -1;
  int varID;
  int nrecs;
  int tsID, tsID2, recID, levelID;
  char filesuffix[32];
  char filename[8192];
  int vdate;
  int day;
  int year1, year2;
  int mon1, mon2;
  int lcopy = FALSE;
  int gridsize;
  int ic = 0;
  int cyear[MAX_YEARS];
  int nmiss;
  int gridID;
  int nlevel;
  double *array = NULL;
  field_t **vars = NULL;

  cdoInitialize(argument);

  if ( processSelf() != 0 ) cdoAbort("This operator can't be combined with other operators!");

  if ( UNCHANGED_RECORD ) lcopy = TRUE;

  int SPLITYEAR    = cdoOperatorAdd("splityear",     0, 10000, NULL);
  int SPLITYEARMON = cdoOperatorAdd("splityearmon",  0,   100, NULL);
  
  int operatorID = cdoOperatorID();
  int operintval = cdoOperatorF2(operatorID);

  memset(cyear, 0, MAX_YEARS*sizeof(int));

  int streamID1 = streamOpenRead(cdoStreamName(0));

  int vlistID1 = streamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  strcpy(filename, cdoStreamName(1)->args);
  int nchars = strlen(filename);

  const char *refname = cdoStreamName(0)->argv[cdoStreamName(0)->argc-1];
  filesuffix[0] = 0;
  cdoGenFileSuffix(filesuffix, sizeof(filesuffix), streamInqFiletype(streamID1), vlistID1, refname);

  // if ( ! lcopy )
    {
      gridsize = vlistGridsizeMax(vlistID1);
      if ( vlistNumber(vlistID1) != CDI_REAL ) gridsize *= 2;
      array = (double*) malloc(gridsize*sizeof(double));
    }

  int nvars = vlistNvars(vlistID1);
  int nconst = 0;
  for ( varID = 0; varID < nvars; varID++ )
    if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) nconst++;

  if ( nconst )
    {
      vars = (field_t **) malloc(nvars*sizeof(field_t *));

      for ( varID = 0; varID < nvars; varID++ )
	{
	  if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT )
	    {
	      gridID  = vlistInqVarGrid(vlistID1, varID);
	      nlevel  = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      gridsize = gridInqSize(gridID);
		  
	      vars[varID] = (field_t*) malloc(nlevel*sizeof(field_t));

	      for ( levelID = 0; levelID < nlevel; levelID++ )
		{
		  field_init(&vars[varID][levelID]);
		  vars[varID][levelID].grid = gridID;
		  vars[varID][levelID].ptr  = (double*) malloc(gridsize*sizeof(double));
		}
	    }
	}
    }

  int index1 = - 1<<31;
  int index2;
  year1 = -1;
  mon1  = -1;
  tsID  = 0;
  tsID2 = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      vdate = taxisInqVdate(taxisID1);
      cdiDecodeDate(vdate, &year2, &mon2, &day);

      if ( operatorID == SPLITYEAR )
	{
	  if ( tsID == 0 || year1 != year2 || (year1 == year2 && mon1 > mon2) )
	    {
	      tsID2 = 0;

	      if ( year1 != year2 ) ic = 0;
	      else                  ic++;

	      if ( year2 >= 0 && year2 < MAX_YEARS )
		{
		  ic = cyear[year2];
		  cyear[year2]++;
		}

	      year1 = year2;

	      if ( streamID2 >= 0 ) streamClose(streamID2);

	      sprintf(filename+nchars, "%04d", year1);
	      if ( ic > 0 ) sprintf(filename+strlen(filename), "_%d", ic+1);
	      if ( filesuffix[0] )
		sprintf(filename+strlen(filename), "%s", filesuffix);
	  
	      if ( cdoVerbose ) cdoPrint("create file %s", filename);

	      argument_t *fileargument = file_argument_new(filename);
	      streamID2 = streamOpenWrite(fileargument, cdoFiletype());
	      file_argument_free(fileargument);

	      streamDefVlist(streamID2, vlistID2);
	    }
	  mon1 = mon2;
	}
      else if ( operatorID == SPLITYEARMON )
	{
	  index2 = (vdate/operintval);
	  
	  if ( tsID == 0 || index1 != index2 )
	    {
	      tsID2 = 0;

	      index1 = index2;

	      if ( streamID2 >= 0 ) streamClose(streamID2);

	      sprintf(filename+nchars, "%04d", index1);
	      //if ( ic > 0 ) sprintf(filename+strlen(filename), "_%d", ic+1);
	      if ( filesuffix[0] )
		sprintf(filename+strlen(filename), "%s", filesuffix);
	  
	      if ( cdoVerbose ) cdoPrint("create file %s", filename);

	      argument_t *fileargument = file_argument_new(filename);
	      streamID2 = streamOpenWrite(fileargument, cdoFiletype());
	      file_argument_free(fileargument);

	      streamDefVlist(streamID2, vlistID2);
	    }
	}
      
      taxisCopyTimestep(taxisID2, taxisID1);
      streamDefTimestep(streamID2, tsID2);

      if ( tsID > 0 && tsID2 == 0 && nconst )
	{
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT )
		{
		  nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
		  for ( levelID = 0; levelID < nlevel; levelID++ )
		    {
		      streamDefRecord(streamID2, varID, levelID);
		      nmiss = vars[varID][levelID].nmiss;
		      streamWriteRecord(streamID2, vars[varID][levelID].ptr, nmiss);
		    }
		}
	    }
	}

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamDefRecord(streamID2,  varID,  levelID);

	  if ( lcopy && !(tsID == 0 && nconst) )
	    {
	      streamCopyRecord(streamID2, streamID1);
	    }
	  else
	    {
	      streamReadRecord(streamID1, array, &nmiss);
	      streamWriteRecord(streamID2, array, nmiss);

	      if ( tsID == 0 && nconst )
		{
		  if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT )
		    {
		      gridID  = vlistInqVarGrid(vlistID1, varID);
		      gridsize = gridInqSize(gridID);
		      memcpy(vars[varID][levelID].ptr, array, gridsize*sizeof(double));
		      vars[varID][levelID].nmiss = nmiss;
		    }
		}
	    }
	}

      tsID2++;
      tsID++;
    }

  streamClose(streamID1);
  streamClose(streamID2);
 
  if ( array ) free(array);

  if ( nconst )
    {
      for ( varID = 0; varID < nvars; varID++ )
	{
	  if ( vlistInqVarTsteptype(vlistID2, varID) == TSTEP_CONSTANT )
	    {
	      nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
	      for ( levelID = 0; levelID < nlevel; levelID++ )
		if ( vars[varID][levelID].ptr )
		  free(vars[varID][levelID].ptr);

	      free(vars[varID]);
	    }
	}

      if ( vars ) free(vars);
    }

  cdoFinish();

  return (0);
}
