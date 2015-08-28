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


void *Pinfo(void *argument)
{
  int PINFO, PINFOV;
  int operatorID;
  int i;
  int indg;
  int varID, recID;
  int gridsize = 0;
  int gridID, zaxisID, code;
  int nrecs;
  int levelID;
  int tsID;
  int taxisID1, taxisID2;
  int streamID1, streamID2;
  int vlistID1, vlistID2;
  int nmiss;
  int ivals = 0, imiss = 0;
  int vdate, vtime;
  char varname[CDI_MAX_NAME];
  char vdatestr[32], vtimestr[32];	  
  double missval;
  double *array1 = NULL, *array2 = NULL;
  double level;
  double arrmin, arrmax, arrmean, arrvar;

  cdoInitialize(argument);

  PINFO  = cdoOperatorAdd("pinfo",  0, 0, NULL);
  PINFOV = cdoOperatorAdd("pinfov", 0, 0, NULL);

  UNUSED(PINFO);

  operatorID = cdoOperatorID();

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  gridsize = vlistGridsizeMax(vlistID1);

  array1 = (double*) malloc(gridsize*sizeof(double));
  array2 = (double*) malloc(gridsize*sizeof(double));

  indg = 0;
  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      vdate = taxisInqVdate(taxisID1);
      vtime = taxisInqVtime(taxisID1);

      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);

      date2str(vdate, vdatestr, sizeof(vdatestr));
      time2str(vtime, vtimestr, sizeof(vtimestr));

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  if ( tsID == 0 && recID == 0 )
	    {
	      if ( operatorID == PINFOV )
		fprintf(stdout, "   Rec :       Date  Time    Varname     Level    Size    Miss :"
			    "     Minimum        Mean     Maximum\n");
	      else
		fprintf(stdout, "   Rec :       Date  Time    Code  Level    Size    Miss :"
			"     Minimum        Mean     Maximum\n");
	    }

	  streamInqRecord(streamID1, &varID, &levelID);
	  streamReadRecord(streamID1, array1, &nmiss);

	  indg += 1;
	  code     = vlistInqVarCode(vlistID1, varID);
	  gridID   = vlistInqVarGrid(vlistID1, varID);
	  zaxisID  = vlistInqVarZaxis(vlistID1, varID);
	  missval  = vlistInqVarMissval(vlistID1, varID);
	  gridsize = gridInqSize(gridID);

	  if ( operatorID == PINFOV ) vlistInqVarName(vlistID1, varID, varname);

	  if ( operatorID == PINFOV )
	    fprintf(stdout, "%6d :%s %s %-8s ", indg, vdatestr, vtimestr, varname);
	  else
	    fprintf(stdout, "%6d :%s %s %3d", indg, vdatestr, vtimestr, code);

	  level = zaxisInqLevel(zaxisID, levelID);
	  fprintf(stdout, " %7g ", level);

	  fprintf(stdout, "%7d %7d :", gridsize, nmiss);

	  if ( gridInqType(gridID) == GRID_SPECTRAL ||
	       (gridsize == 1 && nmiss == 0) )
	    {
	      fprintf(stdout, "            %#12.5g\n", array1[0]);
	    }
	  else
	    {
	      if ( nmiss > 0 )
		{
		  ivals = 0;
		  arrmean = 0;
		  arrvar  = 0;
		  arrmin  =  1.e300;
		  arrmax  = -1.e300;
		  for ( i = 0; i < gridsize; i++ )
		    {
		      if ( !DBL_IS_EQUAL(array1[i], missval) )
			{
			  if ( array1[i] < arrmin ) arrmin = array1[i];
			  if ( array1[i] > arrmax ) arrmax = array1[i];
			  arrmean += array1[i];
			  arrvar  += array1[i]*array1[i];
			  ivals++;
			}
		    }
		  imiss = gridsize - ivals;
		  gridsize = ivals;
		}
	      else
		{
		  arrmean = array1[0];
		  arrvar  = array1[0];
		  arrmin  = array1[0];
		  arrmax  = array1[0];
		  for ( i = 1; i < gridsize; i++ )
		    {
		      if ( array1[i] < arrmin ) arrmin = array1[i];
		      if ( array1[i] > arrmax ) arrmax = array1[i];
		      arrmean += array1[i];
		      arrvar  += array1[i]*array1[i];
		    }
		}

	      if ( gridsize )
		{
		  arrmean = arrmean/gridsize;
		  arrvar  = arrvar/gridsize - arrmean*arrmean;
		  fprintf(stdout, "%#12.5g%#12.5g%#12.5g\n", arrmin, arrmean, arrmax);
		}
	      else
		{
		  fprintf(stdout, "                     nan\n");
		}

	      if ( imiss != nmiss && nmiss > 0 )
		fprintf(stdout, "Found %d of %d missing values!\n", imiss, nmiss);
	    }

	  for ( i = 0; i < gridsize; i++ ) array2[i] = array1[i];

	  streamDefRecord(streamID2,  varID,  levelID);
	  streamWriteRecord(streamID2, array2, nmiss);
	}
      tsID++;
    }

  streamClose(streamID1);
  streamClose(streamID2);

  if ( array1 ) free(array1);
  if ( array2 ) free(array2);

  cdoFinish();

  return (0);
}
