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
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *CDItest(void *argument)
{
  int NCOPY;
  int operatorID;
  int streamID1, streamID2;
  int n;
  int nrecs;
  int tsID1, tsID2, recID, varID, levelID;
  int lcopy = FALSE;
  int gridsize;
  int vlistID1, vlistID2 = -1;
  int nmiss;
  int taxisID1, taxisID2 = CDI_UNDEFID;
  int max_copy = 3;
  double *array = NULL;
  double s_utime, s_stime;
  double e_utime, e_stime;
  double c_cputime = 0, c_usertime = 0, c_systime = 0;

  cdoInitialize(argument);

  NCOPY = cdoOperatorAdd("ncopy",   0, 0, NULL);

  UNUSED(NCOPY);

  //  if ( UNCHANGED_RECORD ) lcopy = TRUE;

  operatorID = cdoOperatorID();

  UNUSED(operatorID);

  //  operatorInputArg("Number of copies");
  if ( operatorArgc() == 1 ) max_copy = parameter2int(operatorArgv()[0]);

  processStartTime(&s_utime, &s_stime);

  n = 0;
  while ( TRUE )
    {
      streamID1 = streamOpenRead(cdoStreamName(0));

      vlistID1 = streamInqVlist(streamID1);
      taxisID1 = vlistInqTaxis(vlistID1);

      streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

      vlistID2 = vlistDuplicate(vlistID1);
      taxisID2 = taxisDuplicate(taxisID1);
      vlistDefTaxis(vlistID2, taxisID2);

      streamDefVlist(streamID2, vlistID2);

      gridsize = vlistGridsizeMax(vlistID1);
      array = (double*) malloc(gridsize*sizeof(double));

      tsID1 = 0;
      tsID2 = 0;
      while ( (nrecs = streamInqTimestep(streamID1, tsID1)) )
	{
	  taxisCopyTimestep(taxisID2, taxisID1);

	  streamDefTimestep(streamID2, tsID2);
	       
	  for ( recID = 0; recID < nrecs; recID++ )
	    {
	      streamInqRecord(streamID1, &varID, &levelID);
	      streamDefRecord(streamID2,  varID,  levelID);
	  
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
	  tsID1++;
	  tsID2++;
	}

      streamClose(streamID1);
      streamClose(streamID2);

      vlistDestroy(vlistID2);
      taxisDestroy(taxisID2);

      if ( array ) free(array);

      n++;

      cdoProcessTime(&e_utime, &e_stime);

      c_usertime = e_utime - s_utime;
      c_systime  = e_stime - s_stime;
      c_cputime  = c_usertime + c_systime;

      s_utime = e_utime;
      s_stime = e_stime;

      cdoPrint("Copy number %d: %.2fs %.2fs %.2fs", n, c_usertime, c_systime, c_cputime);

      if ( n == max_copy ) break;
    }

  cdoFinish();

  return (0);
}
