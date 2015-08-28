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

      Copy       copy            Copy datasets
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "par_io.h"
#include "pstream.h"


void *Copy(void *argument)
{
  int SELALL, SZIP;
  int operatorID;
  int streamID1, streamID2 = CDI_UNDEFID;
  int nrecs;
  int tsID1, tsID2, recID, varID, levelID;
  int lcopy = FALSE;
  int gridsize;
  int vlistID1, vlistID2 = CDI_UNDEFID;
  int nmiss;
  int streamCnt, nfiles, indf;
  int taxisID1, taxisID2 = CDI_UNDEFID;
  int ntsteps, nvars;
  double *array = NULL;
  par_io_t parIO;

  cdoInitialize(argument);

            cdoOperatorAdd("copy",   0, 0, NULL);
  SELALL  = cdoOperatorAdd("selall", 0, 0, NULL);
  SZIP    = cdoOperatorAdd("szip",   0, 0, NULL);

  if ( UNCHANGED_RECORD ) lcopy = TRUE;

  operatorID = cdoOperatorID();

  if ( operatorID == SZIP )
    {
      cdoCompType  = COMPRESS_SZIP;
      cdoCompLevel = 0;
    }

  streamCnt = cdoStreamCnt();
  nfiles = streamCnt - 1;

  tsID2 = 0;
  for ( indf = 0; indf < nfiles; indf++ )
    {
      if ( cdoVerbose ) cdoPrint("Process file: %s", cdoStreamName(indf)->args);

      streamID1 = streamOpenRead(cdoStreamName(indf));

      vlistID1 = streamInqVlist(streamID1);
      taxisID1 = vlistInqTaxis(vlistID1);

      if ( indf == 0 )
	{
	  streamID2 = streamOpenWrite(cdoStreamName(nfiles), cdoFiletype());

	  vlistID2 = vlistDuplicate(vlistID1);
	  taxisID2 = taxisDuplicate(taxisID1);
	  vlistDefTaxis(vlistID2, taxisID2);

	  ntsteps = vlistNtsteps(vlistID1);
	  nvars   = vlistNvars(vlistID1);

	  if ( ntsteps == 1 )
	    {
	      for ( varID = 0; varID < nvars; ++varID )
		if ( vlistInqVarTsteptype(vlistID1, varID) != TSTEP_CONSTANT ) break;
	      
	      if ( varID == nvars ) ntsteps = 0;
	    }

	  if ( ntsteps == 0 && nfiles > 1 )
	    {	      
	      for ( varID = 0; varID < nvars; ++varID )
		vlistDefVarTsteptype(vlistID2, varID, TSTEP_INSTANT);
	    }

	  streamDefVlist(streamID2, vlistID2);

	  gridsize = vlistGridsizeMax(vlistID1);
	  array = (double*) malloc(gridsize*sizeof(double));
	  if ( cdoParIO )
	    {
	      fprintf(stderr, "Parallel reading enabled!\n");
	      parIO.array = (double*) malloc(gridsize*sizeof(double));
	      parIO.array_size = gridsize;
	    }
	}
      else
	{
	  vlistCompare(vlistID1, vlistID2, CMP_ALL);
	}

      tsID1 = 0;
      while ( (nrecs = streamInqTimestep(streamID1, tsID1)) )
	{
	  taxisCopyTimestep(taxisID2, taxisID1);

	  streamDefTimestep(streamID2, tsID2);
	       
	  for ( recID = 0; recID < nrecs; recID++ )
	    { 
	      if ( lcopy && (operatorID == SELALL || operatorID == SZIP) )
		{
		  streamInqRecord(streamID1, &varID, &levelID);
		  streamDefRecord(streamID2,  varID,  levelID);
		  streamCopyRecord(streamID2, streamID1);
		}
	      else
		{
		  if ( cdoParIO )
		    {
		      parIO.recID = recID; parIO.nrecs = nrecs;
		      /* fprintf(stderr, "in1 streamID %d varID %d levelID %d\n", streamID1, varID, levelID);*/
		      parReadRecord(streamID1, &varID, &levelID, array, &nmiss, &parIO);
		      /* fprintf(stderr, "in2 streamID %d varID %d levelID %d\n", streamID1, varID, levelID);*/
		    }
		  else
		    {
		      streamInqRecord(streamID1, &varID, &levelID);
		      streamReadRecord(streamID1, array, &nmiss);
		    }
		  /*
		  if ( cdoParIO )
		    fprintf(stderr, "out1 %d %d %d\n", streamID2,  varID,  levelID);
		  */
		  streamDefRecord(streamID2,  varID,  levelID);
		  streamWriteRecord(streamID2, array, nmiss);
		  /*
		  if ( cdoParIO )
		    fprintf(stderr, "out2 %d %d %d\n", streamID2,  varID,  levelID);
		  */
		}
	    }
	  tsID1++;
	  tsID2++;
	}

      streamClose(streamID1);
    }

  streamClose(streamID2);

  if ( array ) free(array);
  if ( vlistID2 != CDI_UNDEFID ) vlistDestroy(vlistID2);

  cdoFinish();

  return (0);
}
