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


void *Seloperator(void *argument)
{
  int streamID1, streamID2 = CDI_UNDEFID;
  int nrecs;
  int tsID, recID, varID, levelID;
  int vlistID1, vlistID2;
  int taxisID1, taxisID2;
  int scode, sltype;
  double slevel = 0, level;
  int nvars, nlevs, code, zaxisID, selfound = FALSE;
  int levID, ltype = 0;
  int varID2, levelID2;
  int sellevel, selcode, selltype;
  int lcopy = FALSE;
  int gridsize, nmiss;
  double *array = NULL;

  cdoInitialize(argument);

  if ( UNCHANGED_RECORD ) lcopy = TRUE;

  operatorInputArg("code, ltype, level");

  scode  = parameter2int(operatorArgv()[0]);
  sltype = parameter2int(operatorArgv()[1]);

  if ( operatorArgc() == 3 )
    slevel = parameter2double(operatorArgv()[2]);

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);

  nvars = vlistNvars(vlistID1);
  for ( varID = 0; varID < nvars; varID++ )
    {
      code    = vlistInqVarCode(vlistID1, varID);
      zaxisID = vlistInqVarZaxis(vlistID1, varID);
      nlevs   = zaxisInqSize(zaxisID);

      ltype = zaxis2ltype(zaxisID);

      for ( levID = 0; levID < nlevs; levID++ )
	{
	  level = zaxisInqLevel(zaxisID, levID);

	  if ( operatorArgc() == 3 )
	    sellevel = IS_EQUAL(level, slevel);
	  else
	    sellevel = TRUE;

	  if ( scode == -1 || scode == code )
            selcode = TRUE;
	  else
            selcode = FALSE;

	  if ( sltype == -1 || sltype == ltype )
            selltype = TRUE;
	  else
            selltype = FALSE;

	  if ( selcode && selltype && sellevel )
	    {
	      vlistDefFlag(vlistID1, varID, levID, TRUE);
	      selfound = TRUE;
	    }
	}
    }

  if ( selfound == FALSE )
    cdoWarning("Code %d, ltype %d, level %g not found!", scode, sltype, slevel);

  vlistID2 = vlistCreate();
  vlistCopyFlag(vlistID2, vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  if ( ! lcopy )
    {
      gridsize = vlistGridsizeMax(vlistID1);
      array = (double*) malloc(gridsize*sizeof(double));
    }

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);
	       
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  if ( vlistInqFlag(vlistID1, varID, levelID) == TRUE )
	    {
	      varID2   = vlistFindVar(vlistID2, varID);
	      levelID2 = vlistFindLevel(vlistID2, varID, levelID);

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
	    }
	}

      tsID++;
    }

  streamClose(streamID1);
  streamClose(streamID2);

  vlistDestroy(vlistID2);

  if ( ! lcopy )
    if ( array ) free(array);

  cdoFinish();

  return (0);
}
