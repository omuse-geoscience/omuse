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


double intlin(double x, double y1, double x1, double y2, double x2);

static
void isosurface(double isoval, long nlev1, double *lev1, field_t *field3D, field_t *field2D)
{
  long i, k, gridsize, nmiss;
  int lmiss1, lmiss2;
  double missval, val1, val2;
  double *data2D, *data3D;

  gridsize = gridInqSize(field3D->grid);
  nmiss    = field3D->nmiss;
  missval  = field3D->missval;
  data3D   = field3D->ptr;
  data2D   = field2D->ptr;

  for ( i = 0; i < gridsize; ++i )
    {
      data2D[i] = missval;

      for ( k = 0; k < (nlev1-1); ++k )
	{
	  val1 = data3D[k*gridsize+i];
	  val2 = data3D[(k+1)*gridsize+i];

	  if ( nmiss > 0 )
	    {
	      lmiss1 = DBL_IS_EQUAL(val1, missval);
	      lmiss2 = DBL_IS_EQUAL(val2, missval);
	      if ( lmiss1 && lmiss2 ) continue;
	      if ( lmiss1 && IS_EQUAL(isoval, val2) ) data2D[i] = lev1[k+1];
	      if ( lmiss2 && IS_EQUAL(isoval, val1) ) data2D[i] = lev1[k]  ;
	      if ( lmiss1 || lmiss2 ) continue;
	    }

	  if ( (isoval >= val1 && isoval <= val2) || (isoval >= val2 && isoval <= val1) )
	    {
	      if ( IS_EQUAL(val1, val2) )
		data2D[i] = lev1[k];
	      else
		data2D[i] = intlin(isoval, lev1[k], val1, lev1[k+1], val2);

	      break;
	    }
	}
    }

  nmiss = 0;
  for ( i = 0; i < gridsize; ++i )
    if ( DBL_IS_EQUAL(data2D[i], missval) ) nmiss++;

  field2D->missval = missval;
  field2D->nmiss   = nmiss;
}


void *Isosurface(void *argument)
{
  int streamID1, streamID2;
  int vlistID1, vlistID2;
  int gridsize, nlevel = 0;
  int recID, nrecs;
  int gridID;
  int nlev1;
  int i, offset;
  int tsID, varID, levelID;
  int nmiss, nvars;
  int zaxisID, zaxisID1 = -1, zaxisIDsfc, nzaxis;
  double missval;
  double isoval = 0;
  int *vars = NULL;
  int *liso = NULL;
  field_t *vars1 = NULL;
  field_t field;
  double *lev1;
  double *single;
  int taxisID1, taxisID2;

  cdoInitialize(argument);

  operatorInputArg("isoval");

  operatorCheckArgc(1);

  isoval = parameter2double(operatorArgv()[0]);

  if ( cdoVerbose ) cdoPrint("Isoval: %g\n", isoval);


  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  nzaxis  = vlistNzaxis(vlistID1);
  for ( i = 0; i < nzaxis; i++ )
    {
      zaxisID = vlistZaxis(vlistID1, i);
      nlevel  = zaxisInqSize(zaxisID);
      if ( zaxisInqType(zaxisID) != ZAXIS_HYBRID && zaxisInqType(zaxisID) != ZAXIS_HYBRID_HALF )
	if ( nlevel > 1 )
	  {
	    zaxisID1 = zaxisID;
	    break;
	  }
    }

  if ( i == nzaxis ) cdoAbort("No processable variable found!");

  nlev1 = nlevel;
  lev1  = (double*) malloc((nlev1)*sizeof(double));
  zaxisInqLevels(zaxisID1, lev1);

  zaxisIDsfc = zaxisCreate(ZAXIS_SURFACE, 1);
  for ( i = 0; i < nzaxis; i++ )
    if ( zaxisID1 == vlistZaxis(vlistID1, i) )
      vlistChangeZaxisIndex(vlistID2, i, zaxisIDsfc);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  gridsize = vlistGridsizeMax(vlistID1);

  field_init(&field);
  field.ptr = (double*) malloc(gridsize*sizeof(double));

  nvars = vlistNvars(vlistID1);

  liso  = (int*)     malloc(nvars*sizeof(int));
  vars  = (int*)     malloc(nvars*sizeof(int));
  vars1 = (field_t*) malloc(nvars*sizeof(field_t));

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridID   = vlistInqVarGrid(vlistID1, varID);
      zaxisID  = vlistInqVarZaxis(vlistID1, varID);
      gridsize = gridInqSize(gridID);
      nlevel   = zaxisInqSize(zaxisID);
      missval  = vlistInqVarMissval(vlistID1, varID);

      if ( zaxisID == zaxisID1 )
	liso[varID] = TRUE;
      else 
	liso[varID] = FALSE;

      field_init(&vars1[varID]);
      vars1[varID].grid    = gridID;
      vars1[varID].zaxis   = zaxisID;
      vars1[varID].nmiss   = 0;
      vars1[varID].missval = missval;
      vars1[varID].ptr     = (double*) malloc(gridsize*nlevel*sizeof(double));
    }

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);

      for ( varID = 0; varID < nvars; varID++ )
	{
	  vars[varID] = FALSE;
	  vars1[varID].nmiss = 0;
	}

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	  nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	  offset   = gridsize*levelID;
	  single   = vars1[varID].ptr + offset;
	  
	  streamReadRecord(streamID1, single, &nmiss);
	  vars1[varID].nmiss += nmiss;
	  vars[varID] = TRUE;
	}

      for ( varID = 0; varID < nvars; varID++ )
	{
	  if ( vars[varID] )
	    {
	      if ( liso )
		{
		  isosurface(isoval, nlev1, lev1, &vars1[varID], &field);

		  streamDefRecord(streamID2, varID, 0);
		  streamWriteRecord(streamID2, field.ptr, field.nmiss);
		}
	      else
		{
		  gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
		  nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
		  missval  = vlistInqVarMissval(vlistID2, varID);

		  for ( levelID = 0; levelID < nlevel; levelID++ )
		    {
		      offset   = gridsize*levelID;
		      single   = vars1[varID].ptr + offset;

		      nmiss = 0;
		      for ( i = 0; i < gridsize; ++i )
			if ( DBL_IS_EQUAL(single[i], missval) ) nmiss++;

		      streamDefRecord(streamID2, varID, levelID);
		      streamWriteRecord(streamID2, single, nmiss);
		    }
		}
	    }
	}

      tsID++;
    }

  for ( varID = 0; varID < nvars; varID++ ) free(vars1[varID].ptr);
  free(vars1);

  free(vars);
  free(liso);
  if (lev1) free(lev1);

  free(field.ptr);

  streamClose(streamID2);
  streamClose(streamID1);

  vlistDestroy(vlistID2);

  cdoFinish();

  return (0);
}
