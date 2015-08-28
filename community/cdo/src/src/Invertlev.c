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

      Invertlev     invertlev       Invert level
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "error.h"



static
void invertLevDes(int vlistID)
{
  int index, nzaxis;
  int zaxisID1, zaxisID2;
  int nlev;
  int ilev;
  int zaxistype;
  double *yv1, *yv2;
  double *yb1, *yb2;

  nzaxis = vlistNzaxis(vlistID);
  for ( index = 0; index < nzaxis; index++ )
    {
      zaxisID1 = vlistZaxis(vlistID, index);
      zaxisID2 = zaxisDuplicate(zaxisID1);

      zaxistype = zaxisInqType(zaxisID1);

      nlev = zaxisInqSize(zaxisID1);

      if ( nlev < 2 || zaxistype == ZAXIS_HYBRID || zaxistype == ZAXIS_HYBRID_HALF ) continue;

      /* if ( zaxisInqLevels(zaxisID1, NULL) ) */
	{

	  yv1 = (double*) malloc(nlev*sizeof(double));
	  yv2 = (double*) malloc(nlev*sizeof(double));

	  zaxisInqLevels(zaxisID1, yv1);
	  for ( ilev = 0; ilev < nlev; ++ilev ) yv2[nlev-ilev-1] = yv1[ilev];
	  zaxisDefLevels(zaxisID2, yv2);

	  if ( yv2 ) free(yv2);
	  if ( yv1 ) free(yv1);
	}

      if ( zaxisInqLbounds(zaxisID1, NULL) && zaxisInqUbounds(zaxisID1, NULL) )
	{
	  yb1 = (double*) malloc(nlev*sizeof(double));
	  yb2 = (double*) malloc(nlev*sizeof(double));

	  zaxisInqLbounds(zaxisID1, yb1);
	  for ( ilev = 0; ilev < nlev; ++ilev ) yb2[nlev-ilev-1] = yb1[ilev];
	  zaxisDefLbounds(zaxisID2, yb2);

	  zaxisInqUbounds(zaxisID1, yb1);
	  for ( ilev = 0; ilev < nlev; ++ilev ) yb2[nlev-ilev-1] = yb1[ilev];
	  zaxisDefUbounds(zaxisID2, yb2);

	  if ( yb2 ) free(yb2);
	  if ( yb1 ) free(yb1);
	}

      vlistChangeZaxis(vlistID, zaxisID1, zaxisID2);
    }
}


void *Invertlev(void *argument)
{
  int nrecs;
  int recID, varID, levelID;
  int nmiss;
  int nlev, nlevel;
  int gridID, zaxisID, zaxistype, offset;
  int lcopy = FALSE;
  int linvert = FALSE;

  cdoInitialize(argument);

  if ( UNCHANGED_RECORD ) lcopy = TRUE;

  cdoOperatorAdd("invertlev",     func_all, 0, NULL);

  int operatorID = cdoOperatorID();
  int operfunc   = cdoOperatorF1(operatorID);

  int streamID1 = streamOpenRead(cdoStreamName(0));

  int vlistID1 = streamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  if ( operfunc == func_all || operfunc == func_hrd ) invertLevDes(vlistID2);

  int streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  int gridsize = vlistGridsizeMax(vlistID1);

  double *array = (double*) malloc(gridsize*sizeof(double));

  int nvars = vlistNvars(vlistID1);

  double **vardata  = (double**) malloc(nvars*sizeof(double*));
  int **varnmiss = (int**) malloc(nvars*sizeof(int*));

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridID    = vlistInqVarGrid(vlistID1, varID);
      zaxisID   = vlistInqVarZaxis(vlistID1, varID);
      gridsize  = gridInqSize(gridID);
      zaxistype = zaxisInqType(zaxisID);
      nlev      = zaxisInqSize(zaxisID);

      if ( nlev < 2 || zaxistype == ZAXIS_HYBRID || zaxistype == ZAXIS_HYBRID_HALF )
	{
	  vardata[varID]  = NULL;
	  varnmiss[varID] = NULL;
	}
      else
	{
	  linvert = TRUE;
	  vardata[varID]  = (double*) malloc(gridsize*nlev*sizeof(double));
	  varnmiss[varID] = (int*) malloc(nlev*sizeof(int));
	}
    }

  if ( linvert == FALSE ) cdoWarning("No variables with invertable levels found!");

  int tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);
	       
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);

	  if ( vardata[varID] )
	    {    
	      gridID   = vlistInqVarGrid(vlistID1, varID);
	      zaxisID  = vlistInqVarZaxis(vlistID1, varID);
	      gridsize = gridInqSize(gridID);
	      nlevel   = zaxisInqSize(zaxisID);
	      offset   = gridsize*levelID;

	      streamReadRecord(streamID1, vardata[varID]+offset, &nmiss);
	      varnmiss[varID][levelID] = nmiss;
	    }
	  else
	    {
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

      for ( varID = 0; varID < nvars; varID++ )
	{
	  if ( vardata[varID] )
	    {
	      gridID   = vlistInqVarGrid(vlistID1, varID);
	      zaxisID  = vlistInqVarZaxis(vlistID1, varID);
	      gridsize = gridInqSize(gridID);
	      nlevel   = zaxisInqSize(zaxisID);
	      for ( levelID = 0; levelID < nlevel; levelID++ )
		{
		  streamDefRecord(streamID2, varID, levelID);

		  offset   = gridsize*(nlevel-levelID-1);

		  nmiss = varnmiss[varID][nlevel-levelID-1];

		  streamWriteRecord(streamID2, vardata[varID]+offset, nmiss);
		}   
	    }
	}

      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  if ( array ) free(array);

  for ( varID = 0; varID < nvars; varID++ )
    {
      if ( vardata[varID] )
	{
	  free(varnmiss[varID]);
	  free(vardata[varID]);
	}
    }

  free(varnmiss);
  free(vardata);

  cdoFinish();

  return (0);
}
