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

        Timstat3        varquot2test
        Timstat3        meandiff2test
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "statistic.h"


#define  NIN     2
#define  NOUT    1
#define  NFWORK  4
#define  NIWORK  2

void *Timstat3(void *argument)
{
  int VARQUOT2TEST, MEANDIFF2TEST;
  int streamID[NIN];
  int vlistID[NIN], vlistID2 = -1;
  int gridsize;
  int vdate = 0, vtime = 0;
  int nlevs;
  int i, iw, is;
  int varID, recID, levelID, gridID;
  int nmiss3;
  double rconst, risk;
  double fnvals0, fnvals1;
  double missval, missval1, missval2;
  double fractil_1, fractil_2, statistic;
  double temp0, temp1, temp2, temp3;
  int ***iwork[NIWORK];
  field_t **fwork[NFWORK];
  field_t in[NIN], out[NOUT];
  int reached_eof[NIN];
  int n_in = NIN;


  cdoInitialize(argument);

  VARQUOT2TEST  = cdoOperatorAdd("varquot2test",  0, 0, NULL);
  MEANDIFF2TEST = cdoOperatorAdd("meandiff2test", 0, 0, NULL);

  int operatorID = cdoOperatorID();

  operatorInputArg("constant and risk (e.g. 0.05)");
  operatorCheckArgc(2);
  rconst = parameter2double(operatorArgv()[0]);
  risk   = parameter2double(operatorArgv()[1]);

  if ( operatorID == VARQUOT2TEST )
    {
      if ( rconst <= 0 )
	cdoAbort("Constant must be positive!");
      
      if ( risk <= 0 || risk >= 1 )
	cdoAbort("Risk must be greater than 0 and lower than 1!");
    }

  for ( is = 0; is < NIN; ++is )
    {
      streamID[is] = streamOpenRead(cdoStreamName(is));

      vlistID[is] = streamInqVlist(streamID[is]);
      if ( is > 0 )
	{
	  vlistID2 = streamInqVlist(streamID[is]);
	  vlistCompare(vlistID[0], vlistID2, CMP_ALL);
	}
    }

  int vlistID3 = vlistDuplicate(vlistID[0]);

  gridsize = vlistGridsizeMax(vlistID[0]);
  int nvars = vlistNvars(vlistID[0]);
  int nrecs = vlistNrecs(vlistID[0]);
  int nrecs3 = nrecs;
  int *recVarID   = (int*) malloc(nrecs*sizeof(int));
  int *recLevelID = (int*) malloc(nrecs*sizeof(int));

  int taxisID1 = vlistInqTaxis(vlistID[0]);
  int taxisID3 = taxisDuplicate(taxisID1);
 
  vlistDefTaxis(vlistID3, taxisID3);
  int streamID3 = streamOpenWrite(cdoStreamName(2), cdoFiletype());

  streamDefVlist(streamID3, vlistID3);

  for ( i = 0; i < NIN; ++i ) reached_eof[i] = 0;

  for ( i = 0; i < NIN; ++i )
    {
      field_init(&in[i]);
      in[i].ptr = (double*) malloc(gridsize*sizeof(double));
    }
				 
  for ( i = 0; i < NOUT; ++i )
    {
      field_init(&out[i]);
      out[i].ptr = (double*) malloc(gridsize*sizeof(double));
    }
				 
  for ( iw = 0; iw < NFWORK; ++iw )
    fwork[iw] = (field_t **) malloc(nvars*sizeof(field_t *));
  for ( iw = 0; iw < NIWORK; ++iw )
    iwork[iw] = (int ***) malloc(nvars*sizeof(int **));

  for ( varID = 0; varID < nvars; ++varID )
    {
      gridID   = vlistInqVarGrid(vlistID[0], varID);      
      gridsize = vlistGridsizeMax(vlistID[0]);
      nlevs    = zaxisInqSize(vlistInqVarZaxis(vlistID[0], varID));
      missval  = missval1 = vlistInqVarMissval(vlistID[0], varID);
      missval2 = vlistInqVarMissval(vlistID[1], varID); 

      for ( iw = 0; iw < NFWORK; ++iw )
	fwork[iw][varID] = (field_t*) malloc(nlevs*sizeof(field_t));
      for ( iw = 0; iw < NIWORK; ++iw )
	iwork[iw][varID] = (int **) malloc(nlevs*sizeof(int *));  

      for ( levelID = 0; levelID < nlevs; ++levelID )
	{
	  for ( iw = 0; iw < NFWORK; ++iw )
	    {
	      field_init(&fwork[iw][varID][levelID]);
	      fwork[iw][varID][levelID].grid    = gridID;
	      fwork[iw][varID][levelID].nmiss   = 0;
	      fwork[iw][varID][levelID].missval = missval;
	      fwork[iw][varID][levelID].ptr     = (double*) malloc(gridsize*sizeof(double));
	      memset(fwork[iw][varID][levelID].ptr, 0, gridsize*sizeof(double));
	    }

	  for ( iw = 0; iw < NIWORK; ++iw )
	    {
	      iwork[iw][varID][levelID] = (int*) malloc(gridsize*sizeof(int));
	      memset(iwork[iw][varID][levelID], 0, gridsize*sizeof(int));
	    }
	}
    }
 
  int tsID = 0;
  while ( TRUE )
    {
      for ( is = 0; is < NIN; ++is )
	{
	  if ( reached_eof[is] ) continue;

	  nrecs = streamInqTimestep(streamID[is], tsID);
	  if ( nrecs == 0 )
	    {
	      reached_eof[is] = 1;
	      continue;
	    }

	  vdate = taxisInqVdate(taxisID1);
	  vtime = taxisInqVtime(taxisID1);

	  for ( recID = 0; recID < nrecs; recID++ )
	    {
	      streamInqRecord(streamID[is], &varID, &levelID);

	      gridsize = gridInqSize(vlistInqVarGrid(vlistID[is], varID));

	      in[is].missval = vlistInqVarMissval(vlistID[is], varID);

	      if ( tsID == 0 && is == 0 )
		{
		  recVarID[recID] = varID;
		  recLevelID[recID] = levelID;	     	     
		}	 

	      streamReadRecord(streamID[is], in[is].ptr, &in[is].nmiss);
	      
	      for ( i = 0; i < gridsize; ++i )
		{
		  /*
		  if ( ( ! DBL_IS_EQUAL(array1[i], missval1) ) && 
		  ( ! DBL_IS_EQUAL(array2[i], missval2) ) )
		  */
		    {
		      fwork[NIN*is + 0][varID][levelID].ptr[i] += in[is].ptr[i];
		      fwork[NIN*is + 1][varID][levelID].ptr[i] += in[is].ptr[i] * in[is].ptr[i];
		      iwork[is][varID][levelID][i]++;
		    }
		}	 
	    }
	}

      for ( is = 0; is < NIN; ++is )
	if ( ! reached_eof[is] ) break;

      if ( is == NIN ) break;

      tsID++;
    }


  taxisDefVdate(taxisID3, vdate);
  taxisDefVtime(taxisID3, vtime);
  streamDefTimestep(streamID3, 0);

  for ( recID = 0; recID < nrecs3; recID++ )
    {
      varID    = recVarID[recID];
      levelID  = recLevelID[recID];
  
      missval1 = fwork[0][varID][levelID].missval;
      missval2 = missval1;

      if ( operatorID == VARQUOT2TEST )
	{
	  for ( i = 0; i < gridsize; ++i )
	    {
	      fnvals0 = iwork[0][varID][levelID][i];
	      fnvals1 = iwork[1][varID][levelID][i];

	      temp0 = DIV(MUL(fwork[0][varID][levelID].ptr[i], fwork[0][varID][levelID].ptr[i]), fnvals0);
	      temp1 = DIV(MUL(fwork[2][varID][levelID].ptr[i], fwork[2][varID][levelID].ptr[i]), fnvals1);
	      temp2 = SUB(fwork[1][varID][levelID].ptr[i], temp0);
	      temp3 = SUB(fwork[3][varID][levelID].ptr[i], temp1);
	      statistic = DIV(temp2, ADD(temp2, MUL(rconst, temp3)));
	      
	      if ( fnvals0 <= 1 || fnvals1 <= 1 )
		fractil_1 = fractil_2 = missval1;
	      else
		beta_distr_constants((fnvals0 - 1) / 2,
				     (fnvals1 - 1) / 2, 1 - risk,
				     &fractil_1, &fractil_2, __func__);
	      out[0].ptr[i] = DBL_IS_EQUAL(statistic, missval1) ? missval1 : 
	       	              statistic <= fractil_1 || statistic >= fractil_2 ? 1 : 0;
	    }
	}									       
      else if ( operatorID == MEANDIFF2TEST )
	{
	  int j;
	  double fnvals;
	  double fractil;
	  double mean_factor[NIN], var_factor[NIN];
	  double stddev_estimator, mean_estimator, norm, deg_of_freedom;
	  double tmp;
	  
	  mean_factor[0] = 1;
	  mean_factor[1] = -1;
	  var_factor[0] = var_factor[1] = 1;

	  for ( i = 0; i < gridsize; ++i )
	    {
	      temp0 = 0;
	      deg_of_freedom = -n_in;
	      for ( j = 0; j < n_in; j++ )
		{
		  fnvals = iwork[j][varID][levelID][i];
		  tmp   = DIV(MUL(fwork[2*j][varID][levelID].ptr[i], fwork[2*j][varID][levelID].ptr[i]), fnvals);
		  temp0 = ADD(temp0, DIV(SUB(fwork[2*j+1][varID][levelID].ptr[i], tmp), var_factor[j]));
		  deg_of_freedom = ADD(deg_of_freedom, fnvals);
		}

	      if ( !DBL_IS_EQUAL(temp0, missval1) && temp0 < 0 ) /* This is possible because */
		temp0 = 0;	                                 /* of rounding errors       */

	      stddev_estimator = SQRT(DIV(temp0, deg_of_freedom));
	      mean_estimator = -rconst;
	      for ( j = 0; j < n_in; j++ )
		{
		  fnvals = iwork[j][varID][levelID][i];
		  mean_estimator = ADD(mean_estimator,
				       MUL(mean_factor[j],
					   DIV(fwork[2*j][varID][levelID].ptr[i], fnvals)));
		}

	      temp1 = 0;
	      for ( j = 0; j < n_in; j++ )
		{
		  fnvals = iwork[j][varID][levelID][i];
		  temp1 = ADD(temp1, DIV(MUL(MUL(mean_factor[j], mean_factor[j]),
					     var_factor[j]), fnvals));
		}

	      norm = SQRT(temp1);
	      
	      temp2 = DIV(DIV(mean_estimator, norm), stddev_estimator);
	      fractil = deg_of_freedom < 1 ? missval1 :
		student_t_inv (deg_of_freedom, 1 - risk/2, __func__);

	      out[0].ptr[i] = DBL_IS_EQUAL(temp2, missval1)|| DBL_IS_EQUAL(fractil, missval1) ? 
		              missval1 : fabs(temp2) >= fractil;
	    }
	}

      nmiss3 = 0;
      for ( i = 0; i < gridsize; i++ )
	if ( DBL_IS_EQUAL(out[0].ptr[i], missval1) ) nmiss3++;

      streamDefRecord(streamID3, varID, levelID);
      streamWriteRecord(streamID3, out[0].ptr, nmiss3);
    }

  for ( varID = 0; varID < nvars; varID++ )
    {
      nlevs = zaxisInqSize(vlistInqVarZaxis(vlistID[0], varID));
      for ( levelID = 0; levelID < nlevs; levelID++ )
	{
	  for ( iw = 0; iw < NFWORK; ++iw )
	    free(fwork[iw][varID][levelID].ptr);
	  for ( iw = 0; iw < NIWORK; ++iw )
	    free(iwork[iw][varID][levelID]);
	}
    
      for ( iw = 0; iw < NFWORK; ++iw ) free(fwork[iw][varID]);
      for ( iw = 0; iw < NIWORK; ++iw ) free(iwork[iw][varID]);
    }

  for ( iw = 0; iw < NFWORK; iw++ ) free(fwork[iw]);
  for ( iw = 0; iw < NIWORK; iw++ ) free(iwork[iw]);


  streamClose(streamID3);
  for ( is = 0; is < NIN; ++is )
    streamClose(streamID[is]);

  for ( i = 0; i < NIN; ++i ) free(in[i].ptr);
  for ( i = 0; i < NOUT; ++i ) free(out[i].ptr);

  free(recVarID);
  free(recLevelID);
    
  cdoFinish();   
 
  return (0);
}
