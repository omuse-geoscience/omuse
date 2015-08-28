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
#include "grid.h"


void trms(field_t field1, field_t field2, double *dp, field_t *field3)
{
  int i, k, nlev, len, rnmiss = 0;
  int    zaxis    = field1.zaxis;
  int    grid1    = field1.grid;
  double *array1  = field1.ptr;
  int    grid2    = field2.grid;
  double *array2  = field2.ptr;
  double missval1 = field1.missval;
  double missval2 = field2.missval;
  double *w       = field1.weight;
  double rsum = 0, rsumw = 0, ravg = 0, wp;

  nlev   = zaxisInqSize(zaxis);
  len    = gridInqSize(grid1);
  if ( len != gridInqSize(grid2) )
    cdoAbort("fields have different size!");

  for ( k = 0; k < nlev; k++ ) 
    for ( i = 0; i < len; i++ ) 
      {
	wp = w[i]*dp[k*len+i];
	rsum  = ADD(rsum, MUL(wp, MUL(SUB(array2[k*len+i], array1[k*len+i]),
				      SUB(array2[k*len+i], array1[k*len+i]))));
	rsumw = ADD(rsumw, wp);
      }

  ravg = SQRT(DIV(rsum, rsumw));

  if ( DBL_IS_EQUAL(ravg, missval1) ) rnmiss++;

  field3->ptr[0] = ravg;
  field3->nmiss  = rnmiss;
}

void *Trms(void *argument)
{
  int streamID1, streamID2, streamID3;
  int vlistID1, vlistID2, vlistID3;
  int gridID1, gridID3, lastgrid = -1;
  int wstatus = FALSE;
  int code = 0, oldcode = 0;
  int zaxisID;
  int index, ngrids, nzaxis;
  int recID, nrecs;
  int nvars, nlevel;
  int gridsize;
  int i, k;
  int nmiss;
  int tsID, varID, levelID;
  int lim;
  int pcode = 152, pvarID = -1;
  int needWeights = FALSE;
  long offset;
  size_t vctsize = 0;
  const double *vct, *va = NULL, *vb = NULL;
  double dp1, dp2;
  double *dp;
  double *single;
  double **vardata1 = NULL, **vardata2 = NULL;
  double slon, slat;
  double sglval;
  field_t field1, field2, field3;
  int taxisID1, taxisID3;

  cdoInitialize(argument);

  needWeights = TRUE;

  streamID1 = streamOpenRead(cdoStreamName(0));
  streamID2 = streamOpenRead(cdoStreamName(1));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = streamInqVlist(streamID2);

  slon = 0;
  slat = 0;
  gridID3 = gridCreate(GRID_LONLAT, 1);
  gridDefXsize(gridID3, 1);
  gridDefYsize(gridID3, 1);
  gridDefXvals(gridID3, &slon);
  gridDefYvals(gridID3, &slat);

  vlistClearFlag(vlistID1);
  nvars    = vlistNvars(vlistID1);
  for ( varID = 0; varID < nvars; varID++ )
    {
      if ( vlistInqVarCode(vlistID1, varID) == pcode )
	pvarID = varID;
      else
	vlistDefFlag(vlistID1, varID, 0, TRUE);
    }

  if ( pvarID == -1 ) cdoAbort("pressure variable missing!");

  vlistID3 = vlistCreate();
  vlistCopyFlag(vlistID3, vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID3 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID3, taxisID3);

  ngrids = vlistNgrids(vlistID1);
  index = 0;
  gridID1 = vlistGrid(vlistID1, index);
  
  if ( needWeights &&
       gridInqType(gridID1) != GRID_LONLAT &&
       gridInqType(gridID1) != GRID_GAUSSIAN )
    cdoAbort("Unsupported gridtype: %s", gridNamePtr(gridInqType(gridID1)));

  vlistChangeGridIndex(vlistID3, index, gridID3);
  if ( ngrids > 1 ) cdoAbort("Too many different grids!");

  nzaxis = vlistNzaxis(vlistID1);
  for ( index = 0; index < nzaxis; index++ )
    {
      zaxisID = vlistZaxis(vlistID1, index);
      if ( zaxisInqType(zaxisID) == ZAXIS_HYBRID )
	{
	  vctsize = zaxisInqVctSize(zaxisID);
	  vct     = zaxisInqVctPtr(zaxisID);
	  va = vct;
	  vb = vct + vctsize/2;
	  /*
	  for ( i = 0; i < vctsize/2; i++ )
	    fprintf(stdout, "%5d %25.17f %25.17f\n", i, vct[i], vct[vctsize/2+i]);
	  for ( i = 0; i < vctsize/2-1; i++ )
	    fprintf(stdout, "%5d %25.17f %25.17f %25.17f\n", i, vct[i], vct[vctsize/2+i],
		    (va[i+1] + vb[i+1]*101300) - (va[i] + vb[i]*101300));
	  */

	  break;
	}
    }

  if ( vctsize == 0 ) cdoAbort("VCT missing!");

  streamID3 = streamOpenWrite(cdoStreamName(2), cdoFiletype());

  streamDefVlist(streamID3, vlistID3);

  vardata1 = (double**) malloc(nvars*sizeof(double*));
  vardata2 = (double**) malloc(nvars*sizeof(double*));

  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, pvarID));
  nlevel   = vctsize/2 - 1;
  dp = (double*) malloc(gridsize*nlevel*sizeof(double));

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      vardata1[varID] = (double*) malloc(gridsize*nlevel*sizeof(double));
      vardata2[varID] = (double*) malloc(gridsize*nlevel*sizeof(double));
    }

  field_init(&field1);
  field_init(&field2);
  field_init(&field3);

  lim = vlistGridsizeMax(vlistID1);
  field1.weight = NULL;
  if ( needWeights )
    field1.weight = (double*) malloc(lim*sizeof(double));

  field2.weight = NULL;

  field3.ptr  = &sglval;
  field3.grid = gridID3;

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      nrecs = streamInqTimestep(streamID2, tsID);

      taxisCopyTimestep(taxisID3, taxisID1);

      streamDefTimestep(streamID3, tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);

	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	  offset   = gridsize*levelID;
	  single   = vardata1[varID] + offset;

	  streamReadRecord(streamID1, single, &nmiss);
	  if ( nmiss ) cdoAbort("Missing values unsupported for this operator!");

	  streamInqRecord(streamID2, &varID, &levelID);

	  single   = vardata2[varID] + offset;

	  streamReadRecord(streamID2, single, &nmiss);
	  if ( nmiss ) cdoAbort("Missing values unsupported for this operator!");
	}

      gridsize = gridInqSize(vlistInqVarGrid(vlistID1, pvarID));
      for ( i = 0; i < gridsize; i++ )
	{
	  vardata1[pvarID][i] = exp(vardata1[pvarID][i]);
	  vardata2[pvarID][i] = exp(vardata2[pvarID][i]);
	}

      nlevel = vctsize/2 - 1;
      for ( k = 0; k < nlevel; k++ )
	{
	  offset = gridsize*k;
	  for ( i = 0; i < gridsize; i++ )
	    {
	      dp1 = (va[k+1] + vb[k+1]*vardata1[pvarID][i]) - (va[k] + vb[k]*vardata1[pvarID][i]);
	      dp2 = (va[k+1] + vb[k+1]*vardata2[pvarID][i]) - (va[k] + vb[k]*vardata2[pvarID][i]);

	      dp[offset+i] = 0.5 * (dp1 + dp2);
	    }
	}

      for ( varID = 0; varID < nvars; varID++ )
	{
	  field1.ptr = vardata1[varID];
	  field2.ptr = vardata2[varID];

	  field1.zaxis   = vlistInqVarZaxis(vlistID1, varID);
	  field1.grid    = vlistInqVarGrid(vlistID1, varID);
	  field2.grid    = vlistInqVarGrid(vlistID2, varID);
	  if ( needWeights && field1.grid != lastgrid )
	    {
	      lastgrid = field1.grid;
	      wstatus = gridWeights(field1.grid, field1.weight);
	    }
	  code = vlistInqVarCode(vlistID1, varID);
	  if ( wstatus != 0 && tsID == 0 && code != oldcode )
	    cdoWarning("Using constant area weights for code %d!", oldcode=code);

	  field1.missval = vlistInqVarMissval(vlistID1, varID);
	  field2.missval = vlistInqVarMissval(vlistID2, varID);
	  field3.missval = vlistInqVarMissval(vlistID3, varID);

	  trms(field1, field2, dp, &field3);

	  streamDefRecord(streamID3, varID, 0);
	  streamWriteRecord(streamID3, &sglval, field3.nmiss);
	}
      tsID++;
    }

  streamClose(streamID3);
  streamClose(streamID2);
  streamClose(streamID1);

  vlistDestroy(vlistID3);

  if ( field1.weight ) free(field1.weight);

  for ( varID = 0; varID < nvars; varID++ )
    {
      free(vardata1[varID]);
      free(vardata2[varID]);
    }

  free(vardata1);
  free(vardata2);
  free(dp);

  cdoFinish();

  return (0);
}
