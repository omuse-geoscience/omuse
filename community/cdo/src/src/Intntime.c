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

      Intntime   intntime        Time interpolation
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "interpol.h"


void *Intntime(void *argument)
{
  int streamID1, streamID2;
  int nrecs, nvars, nlevel;
  int i, nrecords;
  int tsID, tsIDo, recID, varID, levelID;
  int gridsize;
  int vlistID1, vlistID2;
  int taxisID1, taxisID2;
  int vdate, vtime;
  int vdate1, vtime1;
  int vdate2, vtime2;
  int offset;
  int calendar;
  int numts, it;
  int *recVarID, *recLevelID;
  int **nmiss1, **nmiss2, nmiss3;
  double missval1, missval2;
  juldate_t juldate1, juldate2, juldate;
  double fac1, fac2;
  double *array, *single1, *single2;
  double **vardata1, **vardata2, *vardatap;

  cdoInitialize(argument);

  operatorInputArg("number of timesteps between 2 timesteps");
  if ( operatorArgc() < 1 ) cdoAbort("Too few arguments!");

  numts = parameter2int(operatorArgv()[0]);
  if ( numts < 2 ) cdoAbort("parameter must be greater than 1!");

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  nvars    = vlistNvars(vlistID1);
  nrecords = vlistNrecs(vlistID1);

  recVarID   = (int*) malloc(nrecords*sizeof(int));
  recLevelID = (int*) malloc(nrecords*sizeof(int));

  gridsize = vlistGridsizeMax(vlistID1);
  array = (double*) malloc(gridsize*sizeof(double));

  nmiss1   = (int **) malloc(nvars*sizeof(int *));
  nmiss2   = (int **) malloc(nvars*sizeof(int *));
  vardata1 = (double **) malloc(nvars*sizeof(double *));
  vardata2 = (double **) malloc(nvars*sizeof(double *));

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      nmiss1[varID]   = (int*) malloc(nlevel*sizeof(int));
      nmiss2[varID]   = (int*) malloc(nlevel*sizeof(int));
      vardata1[varID] = (double*) malloc(gridsize*nlevel*sizeof(double));
      vardata2[varID] = (double*) malloc(gridsize*nlevel*sizeof(double));
    }

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  if ( taxisHasBounds(taxisID2) ) taxisDeleteBounds(taxisID2);
  vlistDefTaxis(vlistID2, taxisID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  calendar = taxisInqCalendar(taxisID1);

  tsID = 0;
  tsIDo = 0;
  nrecs = streamInqTimestep(streamID1, tsID++);
  vdate1 = taxisInqVdate(taxisID1);
  vtime1 = taxisInqVtime(taxisID1);
  juldate1 = juldate_encode(calendar, vdate1, vtime1);

  taxisCopyTimestep(taxisID2, taxisID1);
  streamDefTimestep(streamID2, tsIDo++);
  for ( recID = 0; recID < nrecs; recID++ )
    {
      streamInqRecord(streamID1, &varID, &levelID);
      gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
      offset   = gridsize*levelID;
      single1  = vardata1[varID] + offset;
      streamReadRecord(streamID1, single1, &nmiss1[varID][levelID]);

      streamDefRecord(streamID2, varID, levelID);
      streamWriteRecord(streamID2, single1, nmiss1[varID][levelID]);
    }

  while ( (nrecs = streamInqTimestep(streamID1, tsID++)) )
    {
      vdate2 = taxisInqVdate(taxisID1);
      vtime2 = taxisInqVtime(taxisID1);
      juldate2 = juldate_encode(calendar, vdate2, vtime2);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);

	  recVarID[recID]   = varID;
	  recLevelID[recID] = levelID;

	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	  offset   = gridsize*levelID;
	  single2  = vardata2[varID] + offset;
	  streamReadRecord(streamID1, single2, &nmiss2[varID][levelID]);
	}

      for ( it = 1; it < numts; it++ )
	{
	  double seconds = it * juldate_to_seconds(juldate_sub(juldate2, juldate1)) / numts;
	  juldate = juldate_add_seconds((int)lround(seconds), juldate1);

	  juldate_decode(calendar, juldate, &vdate, &vtime);

	  if ( cdoVerbose )
	    {
	      char vdatestr[32], vtimestr[32];	  
	      /*
		cdoPrint("juldate1 %f", juldate_to_seconds(juldate1));
		cdoPrint("juldate  %f", juldate_to_seconds(juldate));
		cdoPrint("juldate2 %f", juldate_to_seconds(juldate2));
	      */
	      date2str(vdate, vdatestr, sizeof(vdatestr));
	      time2str(vtime, vtimestr, sizeof(vtimestr));
	      cdoPrint("%s %s", vdatestr, vtimestr);
	    }

	  taxisDefVdate(taxisID2, vdate);
	  taxisDefVtime(taxisID2, vtime);
	  streamDefTimestep(streamID2, tsIDo++);

	  fac1 = juldate_to_seconds(juldate_sub(juldate2, juldate)) / 
	         juldate_to_seconds(juldate_sub(juldate2, juldate1));
	  fac2 = juldate_to_seconds(juldate_sub(juldate,  juldate1)) / 
	         juldate_to_seconds(juldate_sub(juldate2, juldate1));

	  for ( recID = 0; recID < nrecs; recID++ )
	    {
	      varID    = recVarID[recID];
	      levelID  = recLevelID[recID];
	      gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	      offset   = gridsize*levelID;
	      single1  = vardata1[varID] + offset;
	      single2  = vardata2[varID] + offset;

	      nmiss3 = 0;

	      if ( nmiss1[varID][levelID] > 0 || nmiss2[varID][levelID] > 0 )
		{
		  missval1 = vlistInqVarMissval(vlistID1, varID);
		  missval2 = vlistInqVarMissval(vlistID2, varID);

		  for ( i = 0; i < gridsize; i++ )
		    {
		      if ( !DBL_IS_EQUAL(single1[i], missval1) &&
			   !DBL_IS_EQUAL(single2[i], missval2) )
			array[i] = single1[i]*fac1 + single2[i]*fac2;
		      else if ( DBL_IS_EQUAL(single1[i], missval1) &&
				!DBL_IS_EQUAL(single2[i], missval2) && fac2 >= 0.5 )
			array[i] = single2[i];
		      else if ( DBL_IS_EQUAL(single2[i], missval2) &&
				!DBL_IS_EQUAL(single1[i], missval1) && fac1 >= 0.5 )
			array[i] = single1[i];
		      else
			{
			  array[i] = missval1;
			  nmiss3++;
			}
		    }
		}
	      else
		{
		  for ( i = 0; i < gridsize; i++ )
		    array[i] = single1[i]*fac1 + single2[i]*fac2;
		}

	      streamDefRecord(streamID2, varID, levelID);
	      streamWriteRecord(streamID2, array, nmiss3);
	    }
	}

      taxisDefVdate(taxisID2, vdate2);
      taxisDefVtime(taxisID2, vtime2);
      streamDefTimestep(streamID2, tsIDo++);
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  varID   = recVarID[recID];
	  levelID = recLevelID[recID];

	  gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
	  offset   = gridsize*levelID;
	  single2  = vardata2[varID] + offset;

	  streamDefRecord(streamID2, varID, levelID);
	  streamWriteRecord(streamID2, single2, nmiss2[varID][levelID]);
	}

      for ( varID = 0; varID < nvars; varID++ )
	{
	  vardatap        = vardata1[varID];
	  vardata1[varID] = vardata2[varID];
	  vardata2[varID] = vardatap;
	}

      vdate1 = vdate2;
      vtime1 = vtime2;
      juldate1 = juldate2;
    }

  for ( varID = 0; varID < nvars; varID++ )
    {
      free(nmiss1[varID]);
      free(nmiss2[varID]);
      free(vardata1[varID]);
      free(vardata2[varID]);
    }

  free(nmiss1);
  free(nmiss2);
  free(vardata1);
  free(vardata2);

  if ( array )  free(array);

  streamClose(streamID2);
  streamClose(streamID1);

  cdoFinish();

  return (0);
}
