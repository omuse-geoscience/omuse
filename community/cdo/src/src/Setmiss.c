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

      Setmiss    setmissval      Set a new missing value
      Setmiss    setctomiss      Set constant to missing value
      Setmiss    setmisstoc      Set missing value to constant
      Setmiss    setrtomiss      Set range to missing value
      Setmiss    setvrange       Set range of valid value
*/


#if defined(HAVE_CONFIG_H)
#  include "config.h"
#endif

#if defined(HAVE_ISNAN) && ! defined(__cplusplus)
int isnan(const double x);
#endif

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Setmiss(void *argument)
{
  int nrecs, recID;
  int varID, levelID;
  int nmiss;
  int i;
  double missval, missval2 = 0;
  double rconst = 0, rmin = 0, rmax = 0;

  cdoInitialize(argument);

  int SETMISSVAL = cdoOperatorAdd("setmissval", 0, 0, "missing value");
  int SETCTOMISS = cdoOperatorAdd("setctomiss", 0, 0, "constant");
  int SETMISSTOC = cdoOperatorAdd("setmisstoc", 0, 0, "constant");
  int SETRTOMISS = cdoOperatorAdd("setrtomiss", 0, 0, "range (min, max)");
  int SETVRANGE  = cdoOperatorAdd("setvrange",  0, 0, "range (min, max)");

  int operatorID = cdoOperatorID();

  if ( operatorID == SETMISSVAL )
    {
      operatorCheckArgc(1);
      missval2 = parameter2double(operatorArgv()[0]);
    }
  else if ( operatorID == SETCTOMISS || operatorID == SETMISSTOC )
    {
      operatorCheckArgc(1);
      /*
      if ( operatorArgv()[0][0] == 'n' || operatorArgv()[0][0] == 'N' )
	{
#if ! defined(HAVE_ISNAN)
	  cdoWarning("Function >isnan< not available!");
#endif
	  rconst = 0.0/0.0;
	}
      else
      */
      rconst = parameter2double(operatorArgv()[0]);
    }
  else
    {
      operatorCheckArgc(2);
      rmin = parameter2double(operatorArgv()[0]);
      rmax = parameter2double(operatorArgv()[1]);
    }

  int streamID1 = streamOpenRead(cdoStreamName(0));

  int vlistID1 = streamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  if ( operatorID == SETMISSVAL )
    {
      int nvars = vlistNvars(vlistID2);
      for ( varID = 0; varID < nvars; varID++ )
	vlistDefVarMissval(vlistID2, varID, missval2);
    }
  else if ( operatorID == SETMISSTOC )
    {
      int nvars = vlistNvars(vlistID2);
      for ( varID = 0; varID < nvars; varID++ )
	{
	  missval = vlistInqVarMissval(vlistID2, varID);
	  if ( DBL_IS_EQUAL(rconst, missval) )
	    {
	      cdoWarning("Missing value and constant have the same value!");
	      break;
	    }
	}
    }

  /*
  if ( operatorID == SETVRANGE )
    {
      double range[2];
      range[0] = rmin;
      range[1] = rmax;

      nvars = vlistNvars(vlistID2);
      for ( varID = 0; varID < nvars; varID++ )
	vlistDefAttFlt(vlistID2, varID, "valid_range", DATATYPE_FLT64, 2, range);
    }
  */
  int streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  int gridsize = vlistGridsizeMax(vlistID1);

  double *array = (double*) malloc(gridsize*sizeof(double));

  int tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamReadRecord(streamID1, array, &nmiss);

	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	  missval = vlistInqVarMissval(vlistID1, varID);

	  if ( operatorID == SETMISSVAL )
	    {
	      nmiss = 0;
	      for ( i = 0; i < gridsize; i++ )
		if ( DBL_IS_EQUAL(array[i], missval)  || DBL_IS_EQUAL(array[i], (float)missval) ||
		     DBL_IS_EQUAL(array[i], missval2) || DBL_IS_EQUAL(array[i], (float)missval2) )
		  {
		    array[i] = missval2;
		    nmiss++;
		  }
	    }
	  else if ( operatorID == SETCTOMISS )
	    {
#if defined(HAVE_ISNAN)
	      if ( isnan(rconst) )
		{
		  for ( i = 0; i < gridsize; i++ )
		    if ( isnan(array[i]) )
		      {
			array[i] = missval;
			nmiss++;
		      }
		}
	      else
#endif
		{
		  for ( i = 0; i < gridsize; i++ )
		    if ( DBL_IS_EQUAL(array[i], rconst) || DBL_IS_EQUAL(array[i], (float)rconst) )
		      {
			array[i] = missval;
			nmiss++;
		      }
		}
	    }
	  else if ( operatorID == SETMISSTOC )
	    {
	      nmiss = 0;
	      for ( i = 0; i < gridsize; i++ )
		if ( DBL_IS_EQUAL(array[i], missval) || DBL_IS_EQUAL(array[i], (float)missval) )
		  {
		    array[i] = rconst;
		  }
	    }
	  else if ( operatorID == SETRTOMISS )
	    {
	      for ( i = 0; i < gridsize; i++ )
		if ( array[i] >= rmin && array[i] <= rmax )
		  {
		    array[i] = missval;
		    nmiss++;
		  }
	    }
	  else if ( operatorID == SETVRANGE )
	    {
	      for ( i = 0; i < gridsize; i++ )
		if ( array[i] < rmin || array[i] > rmax ) array[i] = missval;

	      nmiss = 0;
	      for ( i = 0; i < gridsize; i++ )
		if ( DBL_IS_EQUAL(array[i], missval) ) nmiss++;
	    }

	  streamDefRecord(streamID2, varID, levelID);
	  streamWriteRecord(streamID2, array, nmiss);
	}

      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  if ( array ) free(array);

  cdoFinish();

  return (0);
}
