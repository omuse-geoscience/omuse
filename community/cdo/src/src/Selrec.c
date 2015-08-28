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

      Selrec     selrec          Select records
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"  /* processSelf */
#include "pstream.h"
#include "error.h"
#include "util.h"
#include "list.h"


void *Selrec(void *argument)
{
  int streamID1, streamID2;
  int tsID, nrecs;
  int recID, varID, levelID;
  int *intarr, nsel = 0;
  int vlistID1 = -1, vlistID2 = -1;
  int i;
  int recordID;
  int filetype;
  int taxisID1, taxisID2;
  LIST *ilist = listNew(INT_LIST);

  cdoInitialize(argument);

  if ( processSelf() != 0 ) cdoAbort("This operator can't be combined with other operators!");

  operatorInputArg("records");

  nsel = args2intlist(operatorArgc(), operatorArgv(), ilist);

  intarr = (int *) listArrayPtr(ilist);

  if ( cdoVerbose )
    {
      for ( i = 0; i < nsel; i++ )
	cdoPrint("intarr entry: %d %d", i, intarr[i]);
    }

  streamID1 = streamOpenRead(cdoStreamName(0));

  filetype = streamInqFiletype(streamID1);

  if ( filetype == FILETYPE_NC || filetype == FILETYPE_NC2 || filetype == FILETYPE_NC4 || filetype == FILETYPE_NC4C )
    cdoAbort("This operator does not work on netCDF data!");

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  recordID = 0;
  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);
     
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  recordID++;
	  streamInqRecord(streamID1, &varID, &levelID);

	  for ( i = 0; i < nsel; i++ )
	    {
	      if ( recordID == intarr[i] )
		{
		  streamDefRecord(streamID2, varID, levelID);
		  streamCopyRecord(streamID2, streamID1);

		  break;
		}
	    }
	}
       
      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  listDelete(ilist);

  cdoFinish();

  return (NULL);
}
