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

      Harmonic   harmonic        Harmonic
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Harmonic(void *argument)
{
  int gridsize;
  int nrecs;
  int varID, levelID, recID;
  int tsID;
  int i, j;
  int nts;
  int streamID1, streamID2;
  int *streamIDs;
  int vlistID1, vlistID2, taxisID1, taxisID2;
  int nmiss;
  int nchars;
  int offset;
  int nvars, nlevel;
  int vdate = 0, vtime = 0;
  int n_out, nout, n;
  char filesuffix[32];
  char filename[8192];
  const char *refname;
  double missval;
  double sine, cosine;
  double *array;
  double ***work, ***out;

  cdoInitialize(argument);

  operatorInputArg("wave number and wave length of first harmonic in number of timesteps");

  operatorCheckArgc(2);

  n_out = parameter2int(operatorArgv()[0]);
  n     = parameter2int(operatorArgv()[1]);

  if ( n_out > 9 ) cdoAbort("Maximum number of wave numbers is 9!");

  if ( n < 1 || n < 2 * n_out )
    cdoAbort("The wave length must be positive and smaller than\n"
	     "2 times the number of requested harmonics (=%d)!", n_out);

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisCreate(TAXIS_ABSOLUTE);
  vlistDefTaxis(vlistID2, taxisID2);

  refname = cdoStreamName(0)->argv[cdoStreamName(0)->argc-1];
  filesuffix[0] = 0;
  cdoGenFileSuffix(filesuffix, sizeof(filesuffix), streamInqFiletype(streamID1), vlistID1, refname);

  streamIDs = (int*) malloc(n_out*sizeof(int));

  strcpy(filename, cdoStreamName(1)->args);
  nchars = strlen(filename);

  for ( j = 0; j < n_out; ++j )
    {
      sprintf(filename+nchars, "%1d", j+1);
      if ( filesuffix[0] )
	sprintf(filename+nchars+1, "%s", filesuffix);

      argument_t *fileargument = file_argument_new(filename);
      streamID2 = streamOpenWrite(fileargument, cdoFiletype());
      file_argument_free(fileargument);

      streamIDs[j] = streamID2;

      streamDefVlist(streamID2, vlistID2);
    }

  nvars = vlistNvars(vlistID1);

  out  = (double ***) malloc(n_out*sizeof(double **));
  work = (double ***) malloc(2*n_out*sizeof(double **));

  for ( j = 0; j < n_out; ++j )
    {
      out[j] = (double **) malloc(nvars*sizeof(double *));
      for ( varID = 0; varID < nvars; ++varID )
	{
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	  nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	  out[j][varID] = (double*) malloc(gridsize*nlevel*sizeof(double));
	}
    }

  for ( j = 0; j < n_out*2; ++j )
    {
      work[j] = (double **) malloc(nvars*sizeof(double *));
      for ( varID = 0; varID < nvars; ++varID )
	{
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	  nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	  work[j][varID] = (double*) malloc(gridsize*nlevel*sizeof(double));
	  memset(work[j][varID], 0, gridsize*nlevel*sizeof(double));
	}
    }


  gridsize = vlistGridsizeMax(vlistID1);
  array = (double*) malloc(gridsize*sizeof(double));

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      if ( tsID == 0 )
	{
	  vdate = taxisInqVdate(taxisID1);
	  vtime = taxisInqVtime(taxisID1);
	}

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	  nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	  offset   = gridsize*levelID;

	  streamReadRecord(streamID1, array, &nmiss);

	  if ( nmiss > 0 ) cdoAbort("Missing values are not allowed!");

	  for ( j = 0; j < n_out; ++j )
	    {
	      sine   = sin (2 * M_PI * (((j + 1) * (tsID+1)) % n) / n);
	      cosine = cos (2 * M_PI * (((j + 1) * (tsID+1)) % n) / n);
	      for ( i = 0; i < gridsize; i++ )
		{
		  work[j][varID][i+offset]         += array[i] * sine;
		  work[n_out + j][varID][i+offset] += array[i] * cosine;
		}
	    }
	}

      tsID++;
    }

  nts = tsID;

  if ( array ) free(array);

  streamClose(streamID1);


  if ( nts%n )
    {
      cdoAbort("The length of first harmonic (=%d)"
	       " does not divide the number of timesteps (=%d)!", n, nts);
    }

  for ( j = 0; j < n_out && 2*(j+1) < n; j++ )
    {
      for ( varID = 0; varID < nvars; varID++ )
	{
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
	  nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
	  for ( levelID = 0; levelID < nlevel; levelID++ )
	    {
	      offset = gridsize*levelID;
	      for ( i = 0; i < gridsize; i++ )
		out[j][varID][i+offset] = sqrt(work[j][varID][i+offset] * work[j][varID][i+offset] +
					work[n_out+j][varID][i+offset] * work[n_out+j][varID][i+offset]) * 2 / nts;
	    }
	}
    }

  if ( 2*n_out == n )
    {
      for ( varID = 0; varID < nvars; varID++ )
	{
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
	  nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
	  for ( levelID = 0; levelID < nlevel; levelID++ )
	    {
	      offset = gridsize*levelID;
	      for ( i = 0; i < gridsize; i++ )
		out[n_out - 1][varID][i+offset] = work[2 * n_out - 1][varID][i+offset] / nts;
	    }
	}
    }

  nout = n_out;

  taxisDefVdate(taxisID2, vdate);
  taxisDefVtime(taxisID2, vtime);
  for ( j = 0; j < nout; j++ )
    {
      streamID2 = streamIDs[j];
      streamDefTimestep(streamID2, 0);

      for ( varID = 0; varID < nvars; varID++ )
	{
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
	  nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
	  for ( levelID = 0; levelID < nlevel; levelID++ )
	    {
	      offset = gridsize*levelID;
	      streamDefRecord(streamID2, varID, levelID);
	      streamWriteRecord(streamID2, out[j][varID]+offset, 0);
	    }
	}
    }

  for ( j = 0; j < n_out && 2 * (j + 1) < n; j++ )
    {
      for ( varID = 0; varID < nvars; varID++ )
	{
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
	  nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
	  missval  = vlistInqVarMissval(vlistID2, varID);
	  for ( levelID = 0; levelID < nlevel; levelID++ )
	    {
	      offset = gridsize*levelID;
	      for ( i = 0; i < gridsize; i++ )
		{
		  out[j][varID][i+offset] = work[j][varID][i+offset] || work[n_out+j][varID][i+offset]
		              ? atan2 (work[j][varID][i+offset], work[n_out+j][varID][i+offset]) *
		                n / (j + 1) / 2 / M_PI : missval;
	  
		  if ( out[j][varID][i+offset] < 0 )
		    out[j][varID][i+offset] += n / (j + 1.);
		}
	    }
	}
    }

  nout = n_out;
  if ( 2*n_out == n ) nout -= 1;

  taxisDefVdate(taxisID2, vdate);
  taxisDefVtime(taxisID2, vtime);
  for ( j = 0; j < nout; j++ )
    {
      streamID2 = streamIDs[j];
      streamDefTimestep(streamID2, 1);

      for ( varID = 0; varID < nvars; varID++ )
	{
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
	  nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
	  missval  = vlistInqVarMissval(vlistID2, varID);
	  for ( levelID = 0; levelID < nlevel; levelID++ )
	    {
	      offset = gridsize*levelID;
	      streamDefRecord(streamID2, varID, levelID);
	      nmiss = 0;
	      for ( i = 0; i < gridsize; i++ )
		if ( DBL_IS_EQUAL(out[j][varID][i+offset], missval) ) nmiss++;
	      streamWriteRecord(streamID2, out[j][varID]+offset, nmiss);
	    }
	}
    }

  for ( j = 0; j < n_out; j++ )
    {
      streamID2 = streamIDs[j];
      streamClose(streamID2);
    }

  free(streamIDs);

  for ( j = 0; j < n_out; ++j )
    {
      for ( varID = 0; varID < nvars; ++varID )
	free(out[j][varID]);

      free(out[j]);
    }

  free(out);

  for ( j = 0; j < n_out*2; ++j )
    {
      for ( varID = 0; varID < nvars; ++varID )
	free(work[j][varID]);

      free(work[j]);
    }

  free(work);

  cdoFinish();

  return (0);
}
