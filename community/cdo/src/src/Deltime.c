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

#include <ctype.h>

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Deltime(void *argument)
{
  int DELDAY, DEL29FEB;
  int operatorID;
  int streamID1, streamID2;
  int tsID, tsID2, nrecs;
  int recID, varID, levelID;
  int vlistID1, vlistID2;
  int taxisID1, taxisID2;
  int vdate /*, vtime */;
  int copytimestep;
  int lcopy = FALSE;
  int gridsize;
  int nmiss;
  int nfound;
  int year, month, day;
  int dday, dmon;
  double *array = NULL;
  const char *cmons[]={"", "jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec"};

  cdoInitialize(argument);

  DELDAY   = cdoOperatorAdd("delday",   0, 0, NULL);
  DEL29FEB = cdoOperatorAdd("del29feb", 0, 0, NULL);

  UNUSED(DELDAY);

  operatorID = cdoOperatorID();

  if ( operatorID == DEL29FEB )
    {
      dday = 29;
      dmon = 2;
    }
  else
    {
      int im;
      int nsel;
      char *sarg;
      nsel = operatorArgc();
      if ( nsel < 1 ) cdoAbort("Too few arguments!");
      if ( nsel > 1 ) cdoAbort("Too many arguments!");
      sarg = operatorArgv()[0];
      dday = atoi(sarg);
      dmon = 0;
      while ( isdigit(*sarg) ) sarg++;
      if ( isalpha(*sarg) )
	{
	  char smon[32];
	  strncpy(smon, sarg, sizeof(smon)-1);
	  smon[sizeof(smon)-1] = 0;
	  strtolower(smon);
	  for ( im = 0; im < 12; ++im )
	    if ( memcmp(smon, cmons[im+1], 3) == 0 ) break;

	  if ( im < 12 ) dmon = im + 1;
	}
    }

  if ( cdoVerbose ) cdoPrint("delete day %d%s", dday, cmons[dmon]);

  if ( UNCHANGED_RECORD ) lcopy = TRUE;

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  taxisDefCalendar(taxisID2, CALENDAR_365DAYS);
  vlistDefTaxis(vlistID2, taxisID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  if ( ! lcopy )
    {
      gridsize = vlistGridsizeMax(vlistID1);
      array = (double*) malloc(gridsize*sizeof(double));
    }
      
  nfound = 0;
  tsID  = 0;
  tsID2 = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      vdate = taxisInqVdate(taxisID1);
      // vtime = taxisInqVtime(taxisID1);

      cdiDecodeDate(vdate, &year, &month, &day);

      if ( day == dday && (month == dmon || dmon == 0) )
	{
	  nfound++;
	  copytimestep = FALSE;
	  if ( cdoVerbose )
	    cdoPrint("Delete %4.4d-%2.2d-%2.2d at timestep %d", year, month, day, tsID+1);
	}
      else
	copytimestep = TRUE;

      if ( copytimestep )
	{
	  taxisCopyTimestep(taxisID2, taxisID1);

	  streamDefTimestep(streamID2, tsID2++);

	  for ( recID = 0; recID < nrecs; recID++ )
	    {
	      streamInqRecord(streamID1, &varID, &levelID);
	      streamDefRecord(streamID2, varID, levelID);
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

  streamClose(streamID2);
  streamClose(streamID1);

  if ( nfound == 0 )
    cdoWarning("Day %d%s not found!", dday, cmons[dmon]);

  if ( ! lcopy )
    if ( array ) free(array);

  vlistDestroy(vlistID2);

  cdoFinish();

  return (NULL);
}
