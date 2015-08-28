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

      Intyear    intyear         Year interpolation
*/

#include <ctype.h>

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "interpol.h"
#include "list.h"


void *Intyear(void *argument)
{
  int streamID1, streamID2;
  int nrecs;
  int i, iy;
  int tsID, recID, varID, levelID;
  int gridsize;
  int vlistID1, vlistID2, vlistID3;
  int taxisID1, taxisID2, taxisID3;
  int vtime, vdate1, vdate2, vdate3, year1, year2;
  int nmiss1, nmiss2, nmiss3;
  int *iyears, nyears = 0, *streamIDs = NULL;
  int nchars;
  char filesuffix[32];
  char filename[8192];
  const char *refname;
  double fac1, fac2;
  double missval1, missval2;
  double *array1, *array2, *array3;
  LIST *ilist = listNew(INT_LIST);

  cdoInitialize(argument);

  operatorInputArg("years");

  nyears = args2intlist(operatorArgc(), operatorArgv(), ilist);

  iyears = (int *) listArrayPtr(ilist);

  streamIDs = (int*) malloc(nyears*sizeof(int));

  streamID1 = streamOpenRead(cdoStreamName(0));
  streamID2 = streamOpenRead(cdoStreamName(1));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = streamInqVlist(streamID2);
  vlistID3 = vlistDuplicate(vlistID1);

  vlistCompare(vlistID1, vlistID2, CMP_ALL);

  gridsize = vlistGridsizeMax(vlistID1);
  array1 = (double*) malloc(gridsize*sizeof(double));
  array2 = (double*) malloc(gridsize*sizeof(double));
  array3 = (double*) malloc(gridsize*sizeof(double));

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = vlistInqTaxis(vlistID2);
  taxisID3 = taxisDuplicate(taxisID1);
  if ( taxisHasBounds(taxisID3) ) taxisDeleteBounds(taxisID3);
  vlistDefTaxis(vlistID3, taxisID3);

  strcpy(filename, cdoStreamName(2)->args);
  nchars = strlen(filename);

  refname = cdoStreamName(0)->argv[cdoStreamName(0)->argc-1];
  filesuffix[0] = 0;
  cdoGenFileSuffix(filesuffix, sizeof(filesuffix), streamInqFiletype(streamID1), vlistID1, refname);

  for ( iy = 0; iy < nyears; iy++ )
    {
      sprintf(filename+nchars, "%04d", iyears[iy]);
      if ( filesuffix[0] )
	sprintf(filename+nchars+4, "%s", filesuffix);

      argument_t *fileargument = file_argument_new(filename);
      streamIDs[iy] = streamOpenWrite(fileargument, cdoFiletype());
      file_argument_free(fileargument);

      streamDefVlist(streamIDs[iy], vlistID3);
    }

  tsID = 0;
  while ( TRUE )
    {
      nrecs = streamInqTimestep(streamID1, tsID);
      if ( nrecs == 0 ) break;
      nrecs = streamInqTimestep(streamID2, tsID);
      if ( nrecs == 0 ) cdoAbort("Too few timesteps in second inputfile!");

      vtime  = taxisInqVtime(taxisID1);
      vdate1 = taxisInqVdate(taxisID1);
      year1  = vdate1/10000;
      vdate2 = taxisInqVdate(taxisID2);
      year2  = vdate2/10000;

      for ( iy = 0; iy < nyears; iy++ )
	{
	  if ( iyears[iy] < year1 || iyears[iy] > year2 )
	    cdoAbort("Year %d out of bounds (first year %d; last year %d)!",
		     iyears[iy], year1, year2);
	  vdate3 = vdate1 - year1*10000 + iyears[iy]*10000;
	  taxisDefVdate(taxisID3, vdate3);
	  taxisDefVtime(taxisID3, vtime);
	  streamDefTimestep(streamIDs[iy], tsID);
	}

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamInqRecord(streamID2, &varID, &levelID);

	  streamReadRecord(streamID1, array1, &nmiss1);
	  streamReadRecord(streamID2, array2, &nmiss2);

	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));

	  for ( iy = 0; iy < nyears; iy++ )
	    {
	      fac1 = ((double) year2-iyears[iy]) / (year2-year1);
	      fac2 = ((double) iyears[iy]-year1) / (year2-year1);

	      nmiss3 = 0;

	      if ( nmiss1 > 0 || nmiss2 > 0 )
		{
		  missval1 = vlistInqVarMissval(vlistID1, varID);
		  missval2 = vlistInqVarMissval(vlistID2, varID);

		  for ( i = 0; i < gridsize; i++ )
		    {
		      if ( !DBL_IS_EQUAL(array1[i], missval1) &&
			   !DBL_IS_EQUAL(array2[i], missval2) )
			array3[i] = array1[i]*fac1 + array2[i]*fac2;
		      /* 2010-04-19 Uwe Schulzweida: removed 
		      else if ( DBL_IS_EQUAL(array1[i], missval1) &&
				!DBL_IS_EQUAL(array2[i], missval2) && fac2 >= 0.5 )
			array3[i] = array2[i];
		      else if ( DBL_IS_EQUAL(array2[i], missval2) &&
				!DBL_IS_EQUAL(array1[i], missval1) && fac1 >= 0.5 )
			array3[i] = array1[i];
		      */
		      else
			{
			  array3[i] = missval1;
			  nmiss3++;
			}
		    }
		}
	      else
		{
		  for ( i = 0; i < gridsize; i++ )
		    array3[i] = array1[i]*fac1 + array2[i]*fac2;
		}

	      streamDefRecord(streamIDs[iy], varID, levelID);
	      streamWriteRecord(streamIDs[iy], array3, nmiss3);
	    }
	}

      tsID++;
    }

  for ( iy = 0; iy < nyears; iy++ )
    streamClose(streamIDs[iy]);
  
  streamClose(streamID2);
  streamClose(streamID1);

  if ( array3 )  free(array3);
  if ( array2 )  free(array2);
  if ( array1 )  free(array1);

  free(streamIDs);

  listDelete(ilist);

  cdoFinish();

  return (0);
}
