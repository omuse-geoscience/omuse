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

      Intlevel   intlevel        Linear level interpolation
      Intlevel   intlevel3d      Linear level interpolation on a 3d vertical coordinates variable
*/

#include <ctype.h>

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "list.h"


static
void vert_interp_lev(int gridsize, double missval, double *vardata1, double *vardata2,
		     int nlev2, int *lev_idx1, int *lev_idx2, double *lev_wgt1, double *lev_wgt2)
{
  int i, ilev;
  int idx1, idx2;
  double wgt1, wgt2;
  double w1, w2;
  double *var1L1, *var1L2, *var2;

  for ( ilev = 0; ilev < nlev2; ++ilev )
    {
      idx1 = lev_idx1[ilev];
      idx2 = lev_idx2[ilev];
      wgt1 = lev_wgt1[ilev];
      wgt2 = lev_wgt2[ilev];
      var2 = vardata2+gridsize*ilev;

      if ( cdoVerbose ) cdoPrint("level %d: idx1=%d idx2=%d wgt1=%g wgt2=%g", ilev, idx1, idx2, wgt1, wgt2);

      var1L1 = vardata1+gridsize*idx1;
      var1L2 = vardata1+gridsize*idx2;

#if defined(_OPENMP)
#pragma omp parallel for shared(gridsize, var2, var1L1, var1L2, wgt1, wgt2, missval) private(i, w1, w2)
#endif
      for ( i = 0; i < gridsize; ++i )
	{
	  w1 = wgt1;
	  w2 = wgt2;
	  if ( DBL_IS_EQUAL(var1L1[i], missval) ) w1 = 0;
	  if ( DBL_IS_EQUAL(var1L2[i], missval) ) w2 = 0;

	  if ( IS_EQUAL(w1, 0) && IS_EQUAL(w2, 0) )
	    {
	      var2[i] = missval;
	    }
	  else if ( IS_EQUAL(w1, 0) )
	    {
	      if ( w2 >= 0.5 )
		var2[i] = var1L2[i];
	      else
		var2[i] = missval;	      
	    }
	  else if ( IS_EQUAL(w2, 0) )
	    {
	      if ( w1 >= 0.5 )
		var2[i] = var1L1[i];
	      else
		var2[i] = missval;	      
	    }
	  else
	    {
	      var2[i] = var1L1[i]*w1 + var1L2[i]*w2;
	    }
	}
    }
}

static
void vert_gen_weights(int expol, int nlev1, double *lev1, int nlev2, double *lev2,
		      int *lev_idx1, int *lev_idx2, double *lev_wgt1, double *lev_wgt2)
{
  int i1, i2;
  int    idx1 = 0, idx2 = 0;
  double val1, val2 = 0;

  for ( i2 = 0; i2 < nlev2; ++i2 )
    {
      for ( i1 = 1; i1 < nlev1; ++i1 )
	{
	  if ( lev1[i1-1] < lev1[i1] )
	    {
	      idx1 = i1-1;
	      idx2 = i1;
	    }
	  else
	    {
	      idx1 = i1;
	      idx2 = i1-1;
	    }
	  val1 = lev1[idx1];
	  val2 = lev1[idx2];

	  if ( lev2[i2] > val1 && lev2[i2] <= val2 ) break;
	}

      if ( i1 == nlev1 ) cdoAbort("Level %g not found!", lev2[i2]);
      else
	{
	  if ( i1-1 == 0 )
	    {
	      lev_idx1[i2] = 1;
	      lev_idx2[i2] = 1;
	      lev_wgt1[i2] = 0;
	      if ( expol || IS_EQUAL(lev2[i2], val2) )
		lev_wgt2[i2] = 1;
	      else
		lev_wgt2[i2] = 0;
	    }
	  else if ( i1 == nlev1-1 )
	    {
	      lev_idx1[i2] = nlev1-2;
	      lev_idx2[i2] = nlev1-2;
	      if ( expol || IS_EQUAL(lev2[i2], val2) )
		lev_wgt1[i2] = 1;
	      else
		lev_wgt1[i2] = 0;
	      lev_wgt2[i2] = 0;
	    }
	  else
	    {
	      lev_idx1[i2] = idx1;
	      lev_idx2[i2] = idx2;
	      lev_wgt1[i2] = (lev1[idx2] - lev2[i2]) / (lev1[idx2] - lev1[idx1]);
	      lev_wgt2[i2] = (lev2[i2] - lev1[idx1]) / (lev1[idx2] - lev1[idx1]);
	    }
	  lev_idx1[i2]--;
	  lev_idx2[i2]--;
	  /*
	  printf("%d %g %d %d %g %g %d %d %g %g\n",
	  i2, lev2[i2], idx1, idx2, lev1[idx1], lev1[idx2], 
	  lev_idx1[i2], lev_idx2[i2], lev_wgt1[i2], lev_wgt2[i2]);
	  */
	}
    }
}


void *Intlevel(void *argument)
{
  int gridsize;
  int recID, nrecs;
  int i, offset;
  int tsID, varID, levelID;
  int nvars;
  int nmiss;
  int zaxisID1 = -1, zaxisID2;
  int gridID, zaxisID;
  int nlev1, nlev2, nlevel = 0, maxlev;
  int lup, ldown;
  int **varnmiss = NULL;
  int *varinterp = NULL;
  int *vars = NULL;
  int expol = FALSE;
  double missval;
  double *lev1 = NULL, *lev2 = NULL;
  double *single1, *single2;
  double **vardata1 = NULL, **vardata2 = NULL;
  int *lev_idx1, *lev_idx2;
  double *lev_wgt1, *lev_wgt2;
  LIST *flist = listNew(FLT_LIST);

  cdoInitialize(argument);

  int INTLEVEL   = cdoOperatorAdd("intlevel",  0, 0, NULL);
  int INTLEVELX  = cdoOperatorAdd("intlevelx", 0, 0, NULL);

  UNUSED(INTLEVEL);

  int operatorID = cdoOperatorID();

  if ( operatorID == INTLEVELX ) expol = TRUE;

  operatorInputArg("levels");

  nlev2 = args2fltlist(operatorArgc(), operatorArgv(), flist);
  lev2  = (double *) listArrayPtr(flist);

  if ( cdoVerbose ) for ( i = 0; i < nlev2; ++i ) printf("lev2 %d: %g\n", i, lev2[i]);

  int streamID1 = streamOpenRead(cdoStreamName(0));

  int vlistID1 = streamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int nzaxis  = vlistNzaxis(vlistID1);
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
  lev1  = (double*) malloc((nlev1+2)*sizeof(double));
  zaxisInqLevels(zaxisID1, lev1+1);

  lup = FALSE;
  ldown = FALSE;
  for ( i = 1; i < nlev1; ++i )
    {
      if ( i == 1 )
	{
	  if ( lev1[i+1] > lev1[i] )
	    lup = TRUE;
	  else if ( lev1[i+1] < lev1[i] )
	    ldown = TRUE;	
	}
      else
	{
	  if ( lup )
	    {
	      if ( !(lev1[i+1] > lev1[i]) ) lup = FALSE;
	    }
	  else if ( ldown )
	    {
	      if ( !(lev1[i+1] < lev1[i]) ) ldown = FALSE;
	    }
	}
    }

  if ( lup )
    {
      lev1[0]       = -1.e33; 
      lev1[nlev1+1] =  1.e33;
    }
  else if ( ldown )
    {
      lev1[0]       =  1.e33; 
      lev1[nlev1+1] = -1.e33;
    }
  else
    cdoWarning("Non monotonic zaxis!");

  if ( cdoVerbose ) for ( i = 0; i < nlev1+2; ++i ) printf("lev1 %d: %g\n", i, lev1[i]);

  lev_idx1 = (int*) malloc(nlev2*sizeof(int));
  lev_idx2 = (int*) malloc(nlev2*sizeof(int));
  lev_wgt1 = (double*) malloc(nlev2*sizeof(double));
  lev_wgt2 = (double*) malloc(nlev2*sizeof(double));

  vert_gen_weights(expol, nlev1+2, lev1, nlev2, lev2, lev_idx1, lev_idx2, lev_wgt1, lev_wgt2);

  zaxisID2 = zaxisCreate(zaxisInqType(zaxisID1), nlev2);
  zaxisDefLevels(zaxisID2, lev2);
  {
    char str[256];
    str[0] = 0;
    zaxisInqName(zaxisID1, str);
    zaxisDefName(zaxisID2, str);
    str[0] = 0;
    zaxisInqLongname(zaxisID1, str);
    if ( str[0] ) zaxisDefLongname(zaxisID2, str);
    str[0] = 0;
    zaxisInqUnits(zaxisID1, str);
    if ( str[0] ) zaxisDefUnits(zaxisID2, str);

    zaxisDefPrec(zaxisID2, zaxisInqPrec(zaxisID1));
  }

  for ( i = 0; i < nzaxis; i++ )
    if ( zaxisID1 == vlistZaxis(vlistID1, i) )
      vlistChangeZaxisIndex(vlistID2, i, zaxisID2);

  int streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  nvars = vlistNvars(vlistID1);

  vars      = (int*) malloc(nvars*sizeof(int));
  vardata1  = (double**) malloc(nvars*sizeof(double*));
  vardata2  = (double**) malloc(nvars*sizeof(double*));
  varnmiss  = (int**) malloc(nvars*sizeof(int*));
  varinterp = (int*) malloc(nvars*sizeof(int));

  maxlev   = nlev1 > nlev2 ? nlev1 : nlev2;

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridID   = vlistInqVarGrid(vlistID1, varID);
      zaxisID  = vlistInqVarZaxis(vlistID1, varID);
      gridsize = gridInqSize(gridID);
      nlevel   = zaxisInqSize(zaxisID);

      vardata1[varID] = (double*) malloc(gridsize*nlevel*sizeof(double));

      if ( zaxisID == zaxisID1 )
	{
	  varinterp[varID] = TRUE;
	  vardata2[varID]  = (double*) malloc(gridsize*nlev2*sizeof(double));
	  varnmiss[varID]  = (int*) malloc(maxlev*sizeof(int));
	  memset(varnmiss[varID], 0, maxlev*sizeof(int));
	}
      else
	{
	  varinterp[varID] = FALSE;
	  vardata2[varID]  = vardata1[varID];
	  varnmiss[varID]  = (int*) malloc(nlevel*sizeof(int));
	}
    }


  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      for ( varID = 0; varID < nvars; ++varID ) vars[varID] = FALSE;

      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	  nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	  offset   = gridsize*levelID;
	  single1  = vardata1[varID] + offset;
	  
	  streamReadRecord(streamID1, single1, &varnmiss[varID][levelID]);
	  vars[varID] = TRUE;
	}

      for ( varID = 0; varID < nvars; varID++ )
	{
	  if ( vars[varID] && varinterp[varID] )
	    {
	      gridID   = vlistInqVarGrid(vlistID1, varID);
	      missval  = vlistInqVarMissval(vlistID1, varID);
	      gridsize = gridInqSize(gridID);

	      vert_interp_lev(gridsize, missval, vardata1[varID], vardata2[varID],
			      nlev2, lev_idx1, lev_idx2, lev_wgt1, lev_wgt2);

	      for ( levelID = 0; levelID < nlev2; levelID++ )
		{
		  gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
		  offset   = gridsize*levelID;
		  single2  = vardata2[varID] + offset;
		  nmiss    = 0;
		  for ( i = 0; i < gridsize; ++i )
		    if ( DBL_IS_EQUAL(single2[i], missval) ) nmiss++;
		  varnmiss[varID][levelID] = nmiss;
		}
	    }
	}

      for ( varID = 0; varID < nvars; varID++ )
	{
	  if ( vars[varID] )
	    {
	      nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
	      for ( levelID = 0; levelID < nlevel; levelID++ )
		{
		  gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
		  offset   = gridsize*levelID;
		  single2  = vardata2[varID] + offset;
		  streamDefRecord(streamID2, varID, levelID);
		  streamWriteRecord(streamID2, single2, varnmiss[varID][levelID]);
		}
	    }
	}

      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  for ( varID = 0; varID < nvars; varID++ )
    {
      free(varnmiss[varID]);
      free(vardata1[varID]);
      if ( varinterp[varID] ) free(vardata2[varID]);
    }

  free(varinterp);
  free(varnmiss);
  free(vardata2);
  free(vardata1);
  free(vars);

  free(lev_idx1);
  free(lev_idx2);
  free(lev_wgt1);
  free(lev_wgt2);

  if (lev1) free(lev1);

  listDelete(flist);

  cdoFinish();

  return (0);
}
