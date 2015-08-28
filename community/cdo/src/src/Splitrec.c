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

      Split      splitrec        Split records
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Splitrec(void *argument)
{
  int nchars;
  int streamID1, streamID2;
  int varID;
  int nrecs;
  int tsID, recID, levelID;
  int varID2, levelID2;
  int vlistID1, vlistID2;
  char filesuffix[32];
  char filename[8192];
  const char *refname;
  int index;
  int lcopy = FALSE;
  int gridsize;
  int nmiss;
  double *array = NULL;

  cdoInitialize(argument);

  if ( processSelf() != 0 ) cdoAbort("This operator can't be combined with other operators!");

  if ( UNCHANGED_RECORD ) lcopy = TRUE;

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);

  nrecs  = vlistNrecs(vlistID1);

  strcpy(filename, cdoStreamName(1)->args);
  nchars = strlen(filename);

  refname = cdoStreamName(0)->argv[cdoStreamName(0)->argc-1];
  filesuffix[0] = 0;
  cdoGenFileSuffix(filesuffix, sizeof(filesuffix), streamInqFiletype(streamID1), vlistID1, refname);

  if ( ! lcopy )
    {
      gridsize = vlistGridsizeMax(vlistID1);
      if ( vlistNumber(vlistID1) != CDI_REAL ) gridsize *= 2;
      array = (double*) malloc(gridsize*sizeof(double));
    }

  index = 0;
  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);

	  vlistClearFlag(vlistID1);
	  vlistDefFlag(vlistID1, varID, levelID, TRUE);

	  vlistID2 = vlistCreate();
	  vlistCopyFlag(vlistID2, vlistID1);

	  index++;
	  sprintf(filename+nchars, "%06d", index);
	  if ( filesuffix[0] )
	    sprintf(filename+nchars+6, "%s", filesuffix);

	  if ( cdoVerbose ) cdoPrint("create file %s", filename);

	  argument_t *fileargument = file_argument_new(filename);
	  streamID2 = streamOpenWrite(fileargument, cdoFiletype());
	  file_argument_free(fileargument);

	  streamDefVlist(streamID2, vlistID2);

	  varID2   = vlistFindVar(vlistID2, varID);
	  levelID2 = vlistFindLevel(vlistID2, varID, levelID);

	  streamDefTimestep(streamID2, 0);
	  streamDefRecord(streamID2, varID2, levelID2);
	  if ( lcopy )
	    {
	      streamCopyRecord(streamID2, streamID1);
	    }
	  else
	    {
	      streamReadRecord(streamID1, array, &nmiss);
	      streamWriteRecord(streamID2, array, nmiss);
	    }

	  streamClose(streamID2);
	  vlistDestroy(vlistID2);
	}

      tsID++;
    }

  streamClose(streamID1);

  if ( ! lcopy )
    if ( array ) free(array);

  cdoFinish();

  return (0);
}
