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

      Replace    replace         Replace variables
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


#define  MAX_VARS  1024

void *Replace(void *argument)
{
  int varID, varID1, varID2;
  int nrecs = 0;
  int recID, levelID, levelID2;
  int nrecs2;
  int nchvars = 0;
  int idx;
  int nlevel1, nlevel2;
  char varname1[CDI_MAX_NAME], varname2[CDI_MAX_NAME];
  int nmiss;
  int gridsize;
  int offset;
  int varlist1[MAX_VARS], varlist2[MAX_VARS];
  int **varlevel = NULL;
  int **varnmiss2 = NULL;
  double **vardata2 = NULL;
  double *parray;

  cdoInitialize(argument);

  int streamID1 = streamOpenRead(cdoStreamName(0));

  int vlistID1 = streamInqVlist(streamID1);
  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID3 = taxisDuplicate(taxisID1);

  int streamID2 = streamOpenRead(cdoStreamName(1));

  int vlistID2 = streamInqVlist(streamID2);

  /* compare all variables in vlistID2 */

  int nvars1 = vlistNvars(vlistID1);
  int nvars2 = vlistNvars(vlistID2);

  for ( varID2 = 0; varID2 < nvars2; varID2++ )
    {
      vlistInqVarName(vlistID2, varID2, varname2);

      for ( varID1 = 0; varID1 < nvars1; varID1++ )
	{
	  vlistInqVarName(vlistID1, varID1, varname1);
	  if ( strcmp(varname1, varname2) == 0 ) break;
	}

      if ( varID1 < nvars1 )
	{
	  int gridsize1, gridsize2, nlevel1, nlevel2;

	  gridsize1 = gridInqSize(vlistInqVarGrid(vlistID1, varID1));
	  nlevel1   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID1));

	  gridsize2 = gridInqSize(vlistInqVarGrid(vlistID2, varID2));
	  nlevel2   = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID2));

	  if ( gridsize1 != gridsize2 )
	    cdoAbort("Variables have different gridsize!");

	  if ( nlevel1 < nlevel2 )
	    cdoAbort("Variables have different number of levels!");

	  if ( cdoVerbose ) cdoPrint("Variable %s replaced.", varname1);

	  varlist1[nchvars] = varID1;
	  varlist2[nchvars] = varID2;
	  nchvars++;
	  if ( nchvars > MAX_VARS ) cdoAbort("Internal problem - too many variables!");
	}
      else
	{
	  cdoPrint("Variable %s not found!", varname2);
	}
    }

  if ( nchvars )
    {
      vardata2  = (double **) malloc(nchvars*sizeof(double *));
      varnmiss2 = (int **) malloc(nchvars*sizeof(int *));
      varlevel  = (int **) malloc(nchvars*sizeof(int *));
      for ( idx = 0; idx < nchvars; idx++ )
	{
	  varID1 = varlist1[idx];
	  varID2 = varlist2[idx];
	  nlevel1  = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID1));
	  nlevel2  = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID2));
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID2));
	  vardata2[idx]  = (double*) malloc(nlevel2*gridsize*sizeof(double));
	  varnmiss2[idx] = (int*) malloc(nlevel2*sizeof(int));
	  varlevel[idx] = (int*) malloc(nlevel1*sizeof(int));
	  /*
	  for ( levelID = 0; levelID < nlevel1; levelID++ )
	    varlevel[idx][levelID] = levelID;
	  */
	  if ( nlevel2 <= nlevel1 )
	    {
	      double *level1 = (double*) malloc(nlevel1*sizeof(double));
	      double *level2 = (double*) malloc(nlevel2*sizeof(double));
	      zaxisInqLevels(vlistInqVarZaxis(vlistID1, varID1), level1);
	      zaxisInqLevels(vlistInqVarZaxis(vlistID2, varID2), level2);

	      for ( levelID = 0; levelID < nlevel1; levelID++ )
		varlevel[idx][levelID] = -1;
	      
	      for ( int l2 = 0; l2 < nlevel2; l2++ )
		{
		  int l1;
		  for ( l1 = 0; l1 < nlevel1; l1++ )
		    if ( IS_EQUAL(level2[l2], level1[l1]) )
		      {
			varlevel[idx][l1] = l2;
			break;
		      }

		  if ( l1 == nlevel1 ) cdoWarning("Level %g not found!", level2[l2]);
		}

	      free(level1);
	      free(level2);
	    }
	}
    }

  int vlistID3 = vlistDuplicate(vlistID1);

  int streamID3 = streamOpenWrite(cdoStreamName(2), cdoFiletype());

  vlistDefTaxis(vlistID3, taxisID3);
  streamDefVlist(streamID3, vlistID3);

  gridsize = vlistGridsizeMax(vlistID1);
  double *array = (double*) malloc(gridsize*sizeof(double));

  int nts2 = vlistNtsteps(vlistID2);

  int tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID3, taxisID1);

      if ( tsID == 0 || (nts2 != 0 && nts2 != 1) )
	{
	  nrecs2 = streamInqTimestep(streamID2, tsID);
	  if ( nrecs2 == 0 )
	    cdoAbort("Input streams have different number of timesteps!");

	  for ( recID = 0; recID < nrecs2; recID++ )
	    {
	      streamInqRecord(streamID2, &varID, &levelID);
	      
	      for ( idx = 0; idx < nchvars; idx++ )
		if ( varlist2[idx] == varID )
		  {
		    gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
		    offset   = gridsize*levelID;
		    parray   = vardata2[idx]+offset;
		    streamReadRecord(streamID2, parray, &nmiss);
		    varnmiss2[idx][levelID] = nmiss;
		    break;
		  }
	    }
	}

      streamDefTimestep(streamID3, tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);

	  parray = array;

	  for ( idx = 0; idx < nchvars; idx++ )
	    if ( varlist1[idx] == varID )
	      {
		levelID2 = varlevel[idx][levelID];
		if ( levelID2 != -1 )
		  {
		    gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
		    offset   = gridsize*levelID2;
		    parray   = vardata2[idx]+offset;
		    nmiss    = varnmiss2[idx][levelID2];
		    break;
		  }
	      }

	  if ( idx == nchvars ) streamReadRecord(streamID1, parray, &nmiss);

	  streamDefRecord(streamID3, varID, levelID);
	  streamWriteRecord(streamID3, parray, nmiss);
	}

      tsID++;
    }

  streamClose(streamID3);
  streamClose(streamID2);
  streamClose(streamID1);
 
  if ( vardata2 )
    {
      for ( idx = 0; idx < nchvars; idx++ )
	{
	  free(vardata2[idx]);
	  free(varnmiss2[idx]);
	  free(varlevel[idx]);
	}

      free(vardata2);
      free(varnmiss2);
      free(varlevel);
    }

  if ( array ) free(array);

  cdoFinish();

  return (0);
}
