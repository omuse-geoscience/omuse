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

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "statistic.h"


void *Tests(void *argument)
{
  int NORMAL, STUDENTT, CHISQUARE, BETA, FISHER;
  int operatorID;
  int streamID1, streamID2 = CDI_UNDEFID;
  int nrecs;
  int i;
  int tsID, recID, varID, levelID;
  int vlistID1, vlistID2;
  int gridsize;
  int nmiss;
  int taxisID1, taxisID2;
  double degree_of_freedom = 0, p = 0, q = 0, n = 0, d = 0;
  double missval;
  double *array1 = NULL, *array2 = NULL;

  cdoInitialize(argument);

  NORMAL    = cdoOperatorAdd("normal",    0, 0, NULL);
  STUDENTT  = cdoOperatorAdd("studentt",  0, 0, "degree of freedom");
  CHISQUARE = cdoOperatorAdd("chisquare", 0, 0, "degree of freedom");
  BETA      = cdoOperatorAdd("beta",      0, 0, "p and q");
  FISHER    = cdoOperatorAdd("fisher",    0, 0, "degree of freedom of nominator and of denominator");

  operatorID = cdoOperatorID();

  if ( operatorID == STUDENTT || operatorID == CHISQUARE )
    {
      operatorInputArg(cdoOperatorEnter(operatorID));

      operatorCheckArgc(1);

      degree_of_freedom = parameter2double(operatorArgv()[0]);

      if ( degree_of_freedom <= 0 )
	cdoAbort("degree of freedom must be positive!");
    }
  else if ( operatorID == BETA )
    {
      operatorInputArg(cdoOperatorEnter(operatorID));

      operatorCheckArgc(2);

      p = parameter2double(operatorArgv()[0]);
      q = parameter2double(operatorArgv()[1]);

      if ( p <= 0 || q <= 0 )
	cdoAbort("p and q must be positive!");
    }
  else if ( operatorID == FISHER )
    {
      operatorInputArg(cdoOperatorEnter(operatorID));

      operatorCheckArgc(2);

      n = parameter2double(operatorArgv()[0]);
      d = parameter2double(operatorArgv()[1]);

      if ( n <= 0 || d <= 0 )
	cdoAbort("both degrees must be positive!");
    }

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  gridsize = vlistGridsizeMax(vlistID1);
  array1 = (double*) malloc(gridsize*sizeof(double));
  array2 = (double*) malloc(gridsize*sizeof(double));

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);
	       
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);	  
	  streamReadRecord(streamID1, array1, &nmiss);

	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	  missval = vlistInqVarMissval(vlistID1, varID);

	  if ( operatorID == NORMAL )
	    {
	      for ( i = 0; i < gridsize; i++ )
		array2[i] = DBL_IS_EQUAL(array1[i], missval) ? missval :
		  normal(array1[i], processInqPrompt());
	    }
	  else if ( operatorID == STUDENTT )
	    {
	      for ( i = 0; i < gridsize; i++ )
		array2[i] = DBL_IS_EQUAL(array1[i], missval) ? missval :
		  student_t(degree_of_freedom, array1[i], processInqPrompt());
	    }
	  else if ( operatorID == CHISQUARE )
	    {
	      for ( i = 0; i < gridsize; i++ )
		array2[i] = DBL_IS_EQUAL(array1[i], missval) ? missval :
		  chi_square(degree_of_freedom, array1[i], processInqPrompt());
	    }
	  else if ( operatorID == BETA )
	    {
	      for ( i = 0; i < gridsize; i++ )
		{
		  if ( array1[i] < 0 || array1[i] > 1 )
		    cdoAbort("Value out of range (0-1)!");

		  array2[i] = DBL_IS_EQUAL(array1[i], missval) ? missval :
		    beta_distr(p, q, array1[i], processInqPrompt());
		}
	    }
	  else if ( operatorID == FISHER )
	    {
	      for ( i = 0; i < gridsize; i++ )
		array2[i] = DBL_IS_EQUAL(array1[i], missval) ? missval :
		  fisher(n, d, array1[i], processInqPrompt());
	    }
	  else
	    {
	      cdoAbort("Internal problem, operator not implemented!");
	    }

	  streamDefRecord(streamID2,  varID,  levelID);
	  streamWriteRecord(streamID2, array2, nmiss);
	}

      tsID++;
    }

  streamClose(streamID1);
  streamClose(streamID2);

  vlistDestroy(vlistID2);

  if ( array1 ) free(array1);
  if ( array2 ) free(array2);

  cdoFinish();

  return (0);
}
