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

     EOF3d        eof3d             3D-EOF in spatial or time space
     EOF3d        eof3dspatial      3D-EOF in spatial space
     EOF3d        eof3dtime         3D-EOF in time space
*/
/*
 * TODO: 
 * Role of the weights for eofs. Should not be mixed up with division with
 * number of contributing values during summation.
 */

#if defined(_OPENMP)
#include <omp.h>
#endif

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "grid.h"
#include "statistic.h"


enum T_EIGEN_MODE get_eigenmode(void);
enum T_WEIGHT_MODE get_weightmode(void);


// NO MISSING VALUE SUPPORT ADDED SO FAR

void *EOF3d(void * argument)
{
  enum {EOF3D_, EOF3D_TIME, EOF3D_SPATIAL};

  int temp_size = 0;
  int i, i2, j, j1, j2, eofID, varID, recID, levelID, tsID;
  int missval_warning=0;
  int nmiss,ngrids,n=0,nlevs=0,npack=0,nts=0;
  int offset;
  int timer_cov = 0, timer_eig = 0;
  int *varID2;

  int calendar = CALENDAR_STANDARD;
  juldate_t juldate;

  double missval=0;
  double sum_w, sum;
  double **cov = NULL;                                /* TODO: covariance matrix / eigenvectors after solving */
  double *eigv;
  double *xvals, *yvals, *zvals;
  double *df1p, *df2p;


  if ( cdoTimer )
    {
      timer_cov  = timer_new("Timeof cov");
      timer_eig  = timer_new("Timeof eig");
    }

  cdoInitialize(argument);
  cdoOperatorAdd("eof3d",        EOF3D_,        0, NULL);
  cdoOperatorAdd("eof3dtime",    EOF3D_TIME,    0, NULL);
  cdoOperatorAdd("eof3dspatial", EOF3D_SPATIAL, 0, NULL);

  int operatorID  = cdoOperatorID();
  int operfunc    = cdoOperatorF1(operatorID);

  operatorInputArg("Number of eigen functions to write out");
  int n_eig       = parameter2int(operatorArgv()[0]);

  enum T_EIGEN_MODE eigen_mode = get_eigenmode();
  enum T_WEIGHT_MODE weight_mode = get_weightmode();

  int streamID1  = streamOpenRead(cdoStreamName(0));
  int vlistID1   = streamInqVlist(streamID1);
  int taxisID1   = vlistInqTaxis(vlistID1);
  int gridID1    = vlistInqVarGrid(vlistID1, 0);
  long gridsize  = vlistGridsizeMax(vlistID1);
  int nvars      = vlistNvars(vlistID1);
  int nrecs      = vlistNrecs(vlistID1);

  double *weight = (double *) malloc(gridsize*sizeof(double));
  for ( i = 0; i < gridsize; ++i ) weight[i] = 1.;

  if ( weight_mode == WEIGHT_ON )
    {
      int wstatus = gridWeights(gridID1, weight);
      if ( wstatus != 0  )
	{
	  weight_mode = WEIGHT_OFF;
	  cdoWarning("Using constant grid cell area weights!");
	}
    }

  /*  eigenvalues */

  if ( operfunc == EOF3D_SPATIAL )
    cdoAbort("Operator not Implemented - use eof3d or eof3dtime instead");

  tsID = 0;

  /* COUNT NUMBER OF TIMESTEPS if EOF3D_ or EOF3D_TIME */
  nts = vlistNtsteps(vlistID1);
  if ( nts == -1 )
    {
      while ( TRUE )
	{
	  nrecs = streamInqTimestep(streamID1, tsID);
	  if ( nrecs == 0 )  break;
	  tsID++;
	}
      
      nts = tsID;
      if ( cdoVerbose ) cdoPrint("Counted %i timeSteps", nts);
    }

  streamClose(streamID1);

  streamID1 = streamOpenRead(cdoStreamName(0));
  vlistID1  = streamInqVlist(streamID1);
  taxisID1  = vlistInqTaxis(vlistID1);

  /* reset the requested number of eigen-function to the maximum if neccessary */
  if ( n_eig > nts )
    {
      cdoWarning("Solving in time-space:");
      cdoWarning("Number of eigen-functions to write out is bigger than number of time-steps.");
      cdoWarning("Setting n_eig to %i.", nts);
      n_eig = nts;
    }

  n = nts;

  if ( cdoVerbose )  cdoPrint("counted %i timesteps",n);

  /* allocation of temporary fields and output structures */
  double *in       = (double *) malloc(gridsize*sizeof(double));
  int **datacounts = (int **) malloc(nvars*sizeof(int*));
  double ***datafields   = (double ***) malloc(nvars*sizeof(double **));
  double ***eigenvectors = (double ***) malloc(nvars*sizeof(double **));
  double ***eigenvalues  = (double ***) malloc(nvars*sizeof(double **));

  for ( varID = 0; varID < nvars; ++varID )
    {
      gridID1             = vlistInqVarGrid(vlistID1, varID);
      gridsize            = vlistGridsizeMax(vlistID1);
      nlevs               = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      temp_size           = gridsize * nlevs;
      missval             = vlistInqVarMissval(vlistID1, varID);

      datacounts[varID]   = (int*) malloc(nlevs*sizeof(int));
      datafields[varID]   = (double **) malloc(nts*sizeof(double *));

      for ( tsID = 0; tsID < nts; tsID++ )
	{
	  datafields[varID][tsID] = (double *) malloc(temp_size*sizeof(double));
	  for ( i = 0; i < temp_size; ++i ) datafields[varID][tsID][i] = 0;
	}
      datacounts[varID] = (int *) malloc(temp_size*sizeof(int));	      
      for( i = 0; i < temp_size; i++) datacounts[varID][i] = 0;
      
      eigenvectors[varID] = (double **) malloc(n_eig*sizeof(double *));
      eigenvalues[varID]  = (double **) malloc(nts*sizeof(double *));

      for ( i = 0; i < n; i++ )
	{
	  if ( i < n_eig )
	    {
	      eigenvectors[varID][i] = (double *) malloc(temp_size*sizeof(double));
	      for ( i2 = 0; i2 < temp_size; ++i2 )
		eigenvectors[varID][i][i2] = missval;
	    }
	  
	  eigenvalues[varID][i]    = (double *) malloc(1*sizeof(double));
	  eigenvalues[varID][i][0] = missval;
	}
    }

  if ( cdoVerbose)
    cdoPrint("allocated eigenvalue/eigenvector with nts=%i, n=%i, gridsize=%i for processing in %s",
	     nts,n,gridsize,"time_space");
  
  tsID = 0;

  /* read the data and create covariance matrices for each var & level */
  while ( TRUE )
    {
      nrecs = streamInqTimestep(streamID1, tsID);
      if ( nrecs == 0 ) break;

      for ( recID = 0; recID < nrecs; recID++ )
        {
          streamInqRecord(streamID1, &varID, &levelID);

          gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));

          missval = vlistInqVarMissval(vlistID1, varID);
          streamReadRecord(streamID1, in, &nmiss);

	  offset = gridsize * levelID;
	  for ( i = 0; i < gridsize; ++i )
	    {
	      if ( ! DBL_IS_EQUAL(in[i], missval ) )
		{
		  datafields[varID][tsID][offset + i] = in[i];
		  datacounts[varID][offset + i]++;
		}
	      else
		{
		  if ( missval_warning == 0 )
		    {
		      cdoWarning("Missing Value Support not Checked for this Operator!");
		      cdoWarning("Does not work with changing locations of missing values in time.");
		      missval_warning = 1;
		    }
		  datafields[varID][tsID][i+offset] = 0;
		}
	    }
        }
      tsID++;
    }

  if ( cdoVerbose ) 
    cdoPrint("Read data for %i variables",nvars);
  
  int *pack = (int*) malloc(temp_size*sizeof(int)); //TODO

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
      nlevs               = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      temp_size = gridsize * nlevs;

      if ( cdoVerbose )  {
	char vname[64];
	vlistInqVarName(vlistID1,varID,&vname[0]);
	cdoPrint("============================================================================");
	cdoPrint("Calculating covariance matrix and SVD for var%i (%s)",varID,vname);
      }

      npack = 0;    // TODO already set to 0

      if ( cdoTimer ) timer_start(timer_cov);
      
      for ( i = 0; i < temp_size ; i++ )
	{
	  if ( datacounts[varID][i] > 1)
	    {
	      pack[npack] = i;
	      npack++;
	    }
	}

      sum_w = 1;
      if ( weight_mode == WEIGHT_ON )
	{
	  sum_w = 0;
	  for ( i = 0; i < npack; i++ )  sum_w += weight[pack[i]];
	}

      if ( npack < 1 ) {
	char vname[64];
	vlistInqVarName(vlistID1,varID,&vname[0]);
	cdoWarning("Refusing to calculate EOF from a single time step for var%i (%s)",varID,&vname[0]);
	continue;
      }

	  
      cov = (double **) malloc(nts*sizeof(double*));
      for ( j1 = 0; j1 < nts; j1++)
	cov[j1] = (double *) malloc(nts*sizeof(double));
      eigv = (double *) malloc(n*sizeof(double));

      if ( cdoVerbose )  {
	cdoPrint("varID %i allocated eigv and cov with nts=%i and n=%i",varID,nts,n);
	cdoPrint("   npack=%i, nts=%i temp_size=%i",npack,nts,temp_size);
      }


#if defined(_OPENMP)
#pragma omp parallel for private(j1,j2,sum,df1p,df2p) default(shared) schedule(static,2000)
#endif 
      for ( j1 = 0; j1 < nts; j1++)
	for ( j2 = j1; j2 < nts; j2++ )
	  {
	    sum = 0;
	    df1p = datafields[varID][j1];
	    df2p = datafields[varID][j2];
	    for ( i = 0; i < npack; i++ )
	      sum += weight[pack[i]%gridsize]*df1p[pack[i]]*df2p[pack[i]];
	    cov[j2][j1] = cov[j1][j2] = sum / sum_w / nts;
	  }
      
      if ( cdoVerbose ) cdoPrint("calculated cov-matrix");

      /* SOLVE THE EIGEN PROBLEM */
      if ( cdoTimer ) timer_stop(timer_cov);


      if ( cdoTimer ) timer_start(timer_eig);

      if ( cdoVerbose ) 
	cdoPrint("Processed correlation matrix for var %2i | npack: %4i",varID,n);

      if ( eigen_mode == JACOBI ) 
	parallel_eigen_solution_of_symmetric_matrix(cov, eigv, n, __func__);
      else 
	eigen_solution_of_symmetric_matrix(cov, eigv, n, __func__);
      /* NOW: cov contains the eigenvectors, eigv the eigenvalues */

      if ( cdoVerbose ) 
	cdoPrint("Processed SVD decomposition for var %i from %i x %i matrix",varID,n,n);

      for( eofID=0; eofID<n; eofID++ )
	eigenvalues[varID][eofID][0] = eigv[eofID];
      
      if ( cdoTimer ) timer_stop(timer_eig);

      for ( eofID = 0; eofID < n_eig; eofID++ )
	{
	  double *eigenvec = eigenvectors[varID][eofID];

#if defined(_OPENMP)
#pragma omp parallel for private(i,j,sum) shared(datafields, eigenvec)
#endif 
	  for ( i = 0; i < npack; i++ )
	    {
	      sum = 0;
	      for ( j = 0; j < nts; j++ )
		sum += datafields[varID][j][pack[i]] * cov[eofID][j];

	      eigenvec[pack[i]] = sum;
	    }
	  // NORMALIZING
	  sum = 0;

#if defined(_OPENMP)
#pragma omp parallel for private(i) default(none) reduction(+:sum) \
  shared(eigenvec,weight,pack,npack,gridsize)
#endif 
	  for ( i = 0; i < npack; i++ )
	    sum +=  weight[pack[i]%gridsize] *
	            eigenvec[pack[i]] * eigenvec[pack[i]];

	  if ( sum > 0 )
	    {
	      sum = sqrt(sum);
#if defined(_OPENMP)
#pragma omp parallel for private(i) default(none) \
  shared(sum,npack,eigenvec,pack)
#endif
	      for( i = 0; i < npack; i++ )
		eigenvec[pack[i]] /= sum;
	    }
	  else
	    {
#if defined(_OPENMP)
#pragma omp parallel for private(i) default(none) \
  shared(eigenvec,pack,missval,npack)
#endif
	      for( i = 0; i < npack; i++ )
		eigenvec[pack[i]] = missval;
	    }
	}     /* for ( eofID = 0; eofID < n_eig; eofID++ )     */

      if ( eigv ) free(eigv);
      for ( i=0; i<n; i++ )
	if ( cov[i] ) 
	  free(cov[i]);
    }         /* for ( varID = 0; varID < nvars; varID++ )    */

  /* write files with eigenvalues (ID3) and eigenvectors (ID2) */

  /*  eigenvalues */
  int streamID2   = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  int taxisID2    = taxisDuplicate(taxisID1);

  int gridID2     = gridCreate(GRID_LONLAT, 1);
  gridDefXsize(gridID2, 1);
  gridDefYsize(gridID2, 1);
  xvals       = (double*) malloc(1*sizeof(double));
  yvals       = (double*) malloc(1*sizeof(double));
  zvals       = (double*) malloc(1*sizeof(double));
  xvals[0]    = 0;
  yvals[0]    = 0;
  zvals[0]    = 0;
  gridDefXvals(gridID2, xvals);
  gridDefYvals(gridID2, yvals);

  int zaxisID2 = zaxisCreate(ZAXIS_GENERIC,1);
  zaxisDefLevels(zaxisID2,zvals);
  zaxisDefName(zaxisID2,"zaxis_Reduced");
  zaxisDefLongname(zaxisID2,"Reduced zaxis from EOF3D - only one eigen value per 3D eigen vector");

  int vlistID2 = vlistCreate();
  taxisDefRdate(taxisID2, 0);
  taxisDefRtime(taxisID2, 0);
  vlistDefTaxis(vlistID2, taxisID2);

  varID2 = (int*) malloc(nvars*sizeof(int));
  for ( varID=0; varID<nvars; varID++ )
    varID2[varID] = vlistDefVar(vlistID2, gridID2, zaxisID2, TSTEP_INSTANT);
  ngrids      = vlistNgrids(vlistID2);
  for ( i = 0; i < ngrids; i++ )
    vlistChangeGridIndex(vlistID2, i, gridID2);

  int streamID3   = streamOpenWrite(cdoStreamName(2), cdoFiletype());

  int vlistID3    = vlistDuplicate(vlistID1);
  int taxisID3    = taxisDuplicate(taxisID1);
  taxisDefRdate(taxisID3, 0);
  taxisDefRtime(taxisID3, 0);
  vlistDefTaxis(vlistID3, taxisID3);

  streamDefVlist(streamID2, vlistID2);
  streamDefVlist(streamID3, vlistID3);

  int vdate = 10101;
  int vtime = 0;
  juldate = juldate_encode(calendar, vdate, vtime);
  for ( tsID = 0; tsID < n; tsID++ )
    {
      juldate = juldate_add_seconds(60, juldate);
      juldate_decode(calendar, juldate, &vdate, &vtime);

      taxisDefVdate(taxisID2, vdate);
      taxisDefVtime(taxisID2, vtime);
      streamDefTimestep(streamID2, tsID);

      if ( tsID < n_eig )
        {
          taxisDefVdate(taxisID3, vdate);
          taxisDefVtime(taxisID3, vtime);
          streamDefTimestep(streamID3, tsID);
        }

      for ( varID = 0; varID < nvars; varID++ )
        {
          nlevs = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
          for ( levelID = 0; levelID < nlevs; levelID++ )
            {
	      offset = levelID * gridsize;
              if ( tsID < n_eig )
                {
                  nmiss = 0;
                  for ( i = 0; i < gridsize; i++ )
                    if ( DBL_IS_EQUAL(eigenvectors[varID][tsID][offset + i], missval) ) nmiss++;

                  streamDefRecord(streamID3, varID, levelID);
                  streamWriteRecord(streamID3, &eigenvectors[varID][tsID][offset], nmiss);
                }
	    }
	  if ( DBL_IS_EQUAL(eigenvalues[varID][tsID][i], missval) ) nmiss = 1;
	  else nmiss = 0;
	  streamDefRecord(streamID2, varID, 0);
	  streamWriteRecord(streamID2, eigenvalues[varID][tsID],nmiss);
        } // for ( varID = 0; ... )
    } // for ( tsID = 0; ... )

  for ( varID = 0; varID < nvars; varID++)
    {
      for( i = 0; i < nts; i++)
	{
	  free(datafields[varID][tsID]);
	  if ( i < n_eig )
	    free(eigenvectors[varID][i]);
	  free(eigenvalues[varID][i]);
	}

      free(datafields[varID]);
      free(datacounts[varID]);
      free(eigenvectors[varID]);
      free(eigenvalues[varID]);
    }

  free(datafields);
  free(datacounts);
  free(eigenvectors);
  free(eigenvalues);
  free(in);

  free(pack);
  free(weight);

  streamClose(streamID3);
  streamClose(streamID2);
  streamClose(streamID1);
  
  cdoFinish();
 
  return (0);
}
