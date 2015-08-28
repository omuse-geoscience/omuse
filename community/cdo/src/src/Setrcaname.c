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

#define  MAX_LINE_LEN  4096

void *Setrcaname(void *argument)
{
  int streamID1, streamID2 = CDI_UNDEFID;
  int nrecs;
  int tsID, recID, varID, levelID;
  int vlistID1, vlistID2;
  int taxisID1, taxisID2;
  char **rcsnames;
  FILE *fp;
  char line[MAX_LINE_LEN];
  char sname[CDI_MAX_NAME], sdescription[CDI_MAX_NAME], sunits[CDI_MAX_NAME];
  int scode, sltype, slevel;
  int nvars;
  int zaxisID, ltype, code, nlev;
  int level;
  int lcopy = FALSE;
  int gridsize, nmiss;
  double *array = NULL;

  cdoInitialize(argument);

  if ( UNCHANGED_RECORD ) lcopy = TRUE;

  operatorInputArg("file name with RCA names");
  rcsnames = operatorArgv();

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  nvars = vlistNvars(vlistID2);

  fp = fopen(rcsnames[0], "r");
  if ( fp != NULL )
    {
      while ( readline(fp, line, MAX_LINE_LEN) )
	{
	  sscanf(line, "%d\t%d\t%d\t%s\t%s\t%s", &scode, &sltype, &slevel, sname, sdescription, sunits);
	  /*
	  printf("%s\n", line);
	  printf("%d:%d:%d:%s:%s:%s\n", scode, sltype, slevel, sname, sdescription, sunits);
	  */
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      code = vlistInqVarCode(vlistID2, varID);
	      zaxisID = vlistInqVarZaxis(vlistID2, varID);
	      nlev = zaxisInqSize(zaxisID);

	      ltype = zaxis2ltype(zaxisID);

	      if ( code == scode )
		{
		  if ( ltype == 105 )
		    {
		      if ( nlev != 1 )
			{
			  cdoWarning("Number of levels should be 1 for level type 105!");
			  cdoWarning("Maybe environment variable SPLIT_LTYPE_105 is not set.");
			  continue;
			}
		      level = (int) zaxisInqLevel(zaxisID, 0);
		      if ( sltype == 105 && slevel == level )
			{
			  vlistDefVarName(vlistID2, varID, sname);
			  vlistDefVarLongname(vlistID2, varID, sdescription);
			  vlistDefVarUnits(vlistID2, varID, sunits);
			  break;
			}
		    }
		  else if ( sltype != 105 )
		    {
		      vlistDefVarName(vlistID2, varID, sname);
		      vlistDefVarLongname(vlistID2, varID, sdescription);
		      vlistDefVarUnits(vlistID2, varID, sunits);
		      break;
		    }
		}
	    }
	}

      fclose(fp);
    }
  else
    {
      perror(rcsnames[0]);
    }

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
