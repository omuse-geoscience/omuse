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

      Intlevel   intlevel3d      Linear level interpolation on a 3d vertical coordinates variable
      Intlevel   intlevelx3d     Linear level interpolation on a 3d vertical coordinates variable with extrapolation
*/

#include <ctype.h>

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "list.h"
#include "after_vertint.h"



void *Intlevel3d(void *argument)
{
  int gridsize,gridSize,gridsizei,gridsizeo;
  int recID, nrecs;
  int i, offset;
  int tsID, varID, levelID;
  int nvars,nvct;
  int nzaxis;
  int nmiss;

  int nlonIn, nlatIn, nlonOut, nlatOut;
  //double *lonIn, *latIn, *lonOut, *latOut;

  int zaxisID1 = -1, zaxisID3;
  int gridID3 = -1, gridID, zaxisID;
  int nlevi, nlevo, nlevel = 0, maxlev;
  int lup, ldown;
  int **varnmiss = NULL;
  int *varinterp = NULL;
  int *vars = NULL;
  int expol = FALSE;
  double missval;
  double *lev1 = NULL, *lev2 = NULL;
  double *single1, *single2;
  double **vardata1 = NULL, **vardata2 = NULL;
  int taxisID1, taxisID3;
  int *lev_idx1, *lev_idx2;
  double *lev_wgt1, *lev_wgt2;
  double *zlevels_in, *zlevels_out;
  int zlevels_in_miss, zlevels_out_miss;
  char varname[10]; 

  cdoInitialize(argument);

  int INTLEVEL3D  = cdoOperatorAdd("intlevel3d",  0, 0, NULL);
  int INTLEVELX3D = cdoOperatorAdd("intlevelx3d",  0, 0, NULL);

  int operatorID = cdoOperatorID();

  if      ( operatorID == INTLEVEL3D )  expol = FALSE;
  else if ( operatorID == INTLEVELX3D ) expol = TRUE;

  operatorInputArg("icoordinate");

  /*  Read filename from Parameter */
  operatorInputArg("filename for vertical source coordinates variable");
  operatorCheckArgc(1);
  argument_t *fileargument = file_argument_new(operatorArgv()[0]);
  int streamID0 = streamOpenRead(fileargument);                     /*  3d vertical input coordinate */
  file_argument_free(fileargument);
  int streamID1 = streamOpenRead(cdoStreamName(0));                 /*  input data */
  int streamID2 = streamOpenRead(cdoStreamName(1));                 /*  3d target vertical coordinate */
  int streamID3 = streamOpenWrite(cdoStreamName(2),cdoFiletype());  /*  output stream */

  int vlistID0 = streamInqVlist(streamID0);
  int vlistID1 = streamInqVlist(streamID1); taxisID1 = vlistInqTaxis(vlistID1);
  int vlistID2 = streamInqVlist(streamID2);
  int vlistID3 = vlistDuplicate(vlistID1);  taxisID3 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID3, taxisID1);

  /*
   * Read 3d input coordinate (streamID0)
   * * two additional levels are added (top + bottom) for later extrapolation checking
   */
  {
    nvars = vlistNvars(vlistID0);
    if (nvars != 1) cdoAbort("Only one single variable is allowed!");

    gridID     = vlistInqVarGrid(vlistID0, 0);
    zaxisID    = vlistInqVarZaxis(vlistID0, 0);
    gridsize   = gridInqSize(gridID);
    nlevel     = zaxisInqSize(zaxisID);

    nlonIn  = gridInqXsize(gridID);
    nlatIn  = gridInqYsize(gridID);
    /*
    lonIn   = (double*) malloc(nlonIn*sizeof(double));
    latIn   = (double*) malloc(nlatIn*sizeof(double));
    gridInqXvals(gridID, lonIn);
    gridInqYvals(gridID, latIn);
    */
    zlevels_in = (double*) malloc(gridsize*(nlevel+2)*sizeof(double));
    nlevi      = nlevel;   /* number of input levels for later use */
    gridsizei  = gridsize; /* horizontal gridsize of input z coordinate */
    nrecs      = streamInqTimestep(streamID0, 0);
    if (cdoVerbose) cdoPrint("%d records input 3d vertical height\n",nrecs);

    for ( recID = 0; recID < nrecs; recID++ )
      {
        streamInqRecord(streamID0, &varID, &levelID);
        gridsize = gridInqSize(vlistInqVarGrid(vlistID0, varID));
        offset   = gridsize + gridsize*levelID;
        single1  = zlevels_in + offset;
        streamReadRecord(streamID0, single1, &zlevels_in_miss);
      }
  }

  /*
   * Read 3d output coordinate (streamID2)
   */
  {
    nvars = vlistNvars(vlistID2);
    if (nvars != 1) cdoAbort("Only one single variable is allowed!");
    gridID      = vlistInqVarGrid(vlistID2, varID);
    gridID3     = gridID;
    zaxisID     = vlistInqVarZaxis(vlistID2, varID);
    gridsize    = gridInqSize(gridID);
    nlevel      = zaxisInqSize(zaxisID);

    nlonOut = gridInqXsize(gridID);
    nlatOut = gridInqYsize(gridID);
    /*
    lonOut  = (double*) malloc(nlonOut*sizeof(double));
    latOut  = (double*) malloc(nlatOut*sizeof(double));
    gridInqXvals(gridID, lonOut);
    gridInqYvals(gridID, latOut);
    */
    zlevels_out = (double*) malloc(gridsize*nlevel*sizeof(double));
    nlevo       = nlevel;  /* number of output levels for later use */
    gridsizeo   = gridsize;/* horizontal gridsize of output z coordinate */
    nrecs       = streamInqTimestep(streamID2, 0);
    if (cdoVerbose) cdoPrint("%d records target 3d vertical height and gridsize %d\n",nrecs,gridsize);

    for ( recID = 0; recID < nrecs; recID++ )
      {
	streamInqRecord(streamID2, &varID, &levelID);
	gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
	offset   = gridsize*levelID;
	single1  = zlevels_out + offset;
	streamReadRecord(streamID2, single1, &zlevels_out_miss);
      }
  }

  /* Missing values are not allowed for coordinate variables */
  if ( 0 != zlevels_in_miss  )
    cdoAbort("Input vertical coordinate variables are not allowd to contain missing values.");
  else
    {
      if ( cdoVerbose ) cdoPrint("Input vertical coordinate has no missing values.");
    }

  if ( 0 != zlevels_out_miss )
    cdoAbort("Output vertical coordinate variables are not allowd to contain missing values.");
  else
    {
      if ( cdoVerbose ) cdoPrint("Output vertical coordinate has no missing values.");
    }


  /*
   * gridsize of input and output vertical coordinate must be equal
   * (later use of gridsizeo ONLY)
   */
  if ( gridsizei != gridsizeo )
    cdoAbort("Input and output vertical coordinate must have the same gridsize");

  gridSize = gridsizeo;

   /* input and output vertical coordinates must have exactly the same horizontal grid */
   if ( nlonIn != nlonOut || 
        nlatIn != nlatOut /*||
        memcmp(lonIn,lonOut,nlonIn*sizeof(double)) ||
        memcmp(latIn,latOut,nlatIn*sizeof(double))*/ )
     {
       /* i =0; printf ( "lonIn:%g latIn:%g lonOut:%g latOut:%g\n",lonIn[i],latIn[i],lonOut[i],latOut[i] ); */
       cdoAbort("Input and output vertical coordinates do NOT exactly have the same horizontal grid.");
     }

  /*
   * Check for the correct vertical axis in the input: Variables with the same
   * number of levels as the input vertical levels from operators parameter
   * (streamID0). Variables with a different z-axis should be copied into output.
   */
  nzaxis  = vlistNzaxis(vlistID1);
  for ( i = 0; i < nzaxis; i++ )
  {
    zaxisID = vlistZaxis(vlistID1, i);
    nlevel  = zaxisInqSize(zaxisID);
    if ( nlevel == nlevi )
    {
      zaxisID1 = zaxisID;
      break;
    }
  }
  if ( i == nzaxis ) cdoAbort("No processable variable found!");

  /*
   * Check monotony of vertical levels
   */
  lup = FALSE;
  ldown = FALSE;
  lev1 = zlevels_in + gridSize;
  for ( i = 0; i < nlevi-1; i++ )
    {
      if ( i == 1 )
	{
	  if ( lev1[(i+1)*gridSize] > lev1[i*gridSize] )
	    lup = TRUE;
	  else if ( lev1[(i+1)*gridSize] < lev1[i*gridSize] )
	    ldown = TRUE;	
	}
      else
	{
	  if ( lup )
	    {
	      if ( !(lev1[(i+1)*gridSize] > lev1[i*gridSize]) ) lup = FALSE;
	    }
	  else if ( ldown )
	    {
	      if ( !(lev1[(i+1)*gridSize] < lev1[i*gridSize]) ) ldown = FALSE;
	    }
	}
    }

  /* Add artificial values for intication of extrapolation areas (lowermost + upmost levels) */
  if ( lup )
    {
      for ( i = 0; i < gridSize ;i++)
      {
        zlevels_in[i]                      = -1.e33;
        zlevels_in[(nlevi+1)*gridSize + i] =  1.e33;
      }
    }
  else if ( ldown )
    {
      for ( i = 0; i < gridSize ;i++)
      {
        zlevels_in[i]                      =  1.e33;
        zlevels_in[(nlevi+1)*gridSize + i] = -1.e33;
      }
    }
  else
    cdoWarning("Non monotonic zaxis!");

  /*
   * Create weights for later interpolation - assumption: input vertical correct is constant in time
   */
  lev_idx1 = (int*) malloc(nlevo*gridSize*sizeof(int));
  lev_idx2 = (int*) malloc(nlevo*gridSize*sizeof(int));
  lev_wgt1 = (double*) malloc(nlevo*gridSize*sizeof(double));
  lev_wgt2 = (double*) malloc(nlevo*gridSize*sizeof(double));

  vert_gen_weights3d(expol, nlevi+2, gridSize, zlevels_in, nlevo, zlevels_out, lev_idx1, lev_idx2, lev_wgt1, lev_wgt2);

  /*
   * Copy z-axis information to output z-axis
   */
  zaxisID3 = zaxisCreate(zaxisInqType(zaxisID1), nlevo);
  lev2 = (double*) malloc(nlevo*sizeof(double));
  /* fill values with its indices */
  for (i=0;i<nlevo;i++)
    lev2[i] = (double) i;
  zaxisDefLevels(zaxisID3, lev2);
  zaxisDefName(zaxisID3, "lev");
  /*  copy VCT from input vlistID1 to output vlistID3 if there is one */
  nvct = zaxisInqVctSize(zaxisID1);
  if ( nvct > 0 ) zaxisDefVct(zaxisID3,zaxisInqVctSize(zaxisID1), zaxisInqVctPtr(zaxisID1));

  for ( i = 0; i < nzaxis; i++ )
    if ( zaxisID1 == vlistZaxis(vlistID1, i) )
      vlistChangeZaxisIndex(vlistID3, i, zaxisID3);
  /* add the vertical output field to the output stream */
  int oz3dvarID = vlistDefVar(vlistID3, gridID3, zaxisID3, TSTEP_INSTANT);
  {
    char str[256];
    str[0] = 0;
    vlistInqVarName(vlistID2,0,str);
    vlistDefVarName(vlistID3,oz3dvarID,str);
    str[0] = 0;
    vlistInqVarLongname(vlistID2,0,str);
    if ( str[0] ) vlistDefVarLongname(vlistID3,oz3dvarID, str);
    str[0] = 0;
    vlistInqVarUnits(vlistID2,0, str);
    if ( str[0] ) vlistDefVarUnits(vlistID3,oz3dvarID, str);
  }

  streamDefVlist(streamID3, vlistID3);

  maxlev   = nlevi > nlevo ? nlevi : nlevo;
  nvars     = vlistNvars(vlistID1);
  vars      = (int*) malloc(nvars*sizeof(int));
  vardata1  = (double**) malloc(nvars*sizeof(double*)); /* input                                         */
  vardata2  = (double**) malloc(nvars*sizeof(double*)); /* output                                        */
  varnmiss  = (int**) malloc(nvars*sizeof(int*));    /* can for missing values of arbitrary variables */
  varinterp = (int*) malloc(nvars*sizeof(int));     /* marker for variables to be interpolated       */

  /* by default no variable should be interpolated */
  for ( i = 0; i < nvars; i++ )
    varinterp[varID] = FALSE;

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridID   = vlistInqVarGrid(vlistID1, varID);
      zaxisID  = vlistInqVarZaxis(vlistID1, varID);
      gridsize = gridInqSize(gridID);
      nlevel   = zaxisInqSize(zaxisID);

      vlistInqVarName(vlistID1, varID, varname);

      vardata1[varID] = (double*) malloc(gridsize*nlevel*sizeof(double));

      /*  variabls for interpolation:
       *  * have the required vertical axis, i.e. the correct number of levels (nlevi)
       *  * have the same number of horizontal grid points (i.e. same gridSize) like the two vertical coordinates
       *  * are NOT the output vertical coordinates itself
       */
      if ( zaxisID == zaxisID1 && varID != oz3dvarID && gridsize == gridSize )
        {
          nlonIn  = gridInqXsize(gridID);
          nlatIn  = gridInqYsize(gridID);
	  /*
          lonIn   = (double*) malloc(nlonIn*sizeof(double));
          latIn   = (double*) malloc(nlatIn*sizeof(double));
          gridInqXvals(gridID, lonIn);
          gridInqYvals(gridID, latIn);
	  */
          if ( nlonIn != nlonOut || 
               nlatIn != nlatOut /*||
               memcmp(lonIn,lonOut,nlonIn*sizeof(double)) ||
               memcmp(latIn,latOut,nlatIn*sizeof(double))*/ )
            {
              varinterp[varID] = FALSE;
              vardata2[varID]  = vardata1[varID];
              varnmiss[varID]  = (int*) malloc(nlevel*sizeof(int));
              if ( cdoVerbose ) cdoPrint("Ignore variable %s with %d levels\n",varname,nlevel);
            }
          else
            {
              varinterp[varID] = TRUE;
              vardata2[varID]  = (double*) malloc(gridsize*nlevo*sizeof(double));
              varnmiss[varID]  = (int*) malloc(maxlev*sizeof(int));
              memset(varnmiss[varID], 0, maxlev*sizeof(int));
            }
        }
      else
        {
          varinterp[varID] = FALSE;
          vardata2[varID]  = vardata1[varID];
          varnmiss[varID]  = (int*) malloc(nlevel*sizeof(int));
          if ( cdoVerbose ) cdoPrint("Ignore variable %s with %d levels\n",varname,nlevel);
        }
      }

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      for ( varID = 0; varID < nvars; ++varID ) vars[varID] = FALSE;

      taxisCopyTimestep(taxisID3, taxisID1);

      streamDefTimestep(streamID3, tsID);

      /*
       * Read the whole 3d data field
       */
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
          vlistInqVarName(vlistID1, varID, varname); 
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	  nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	  offset   = gridsize*levelID;
	  single1  = vardata1[varID] + offset;
          streamReadRecord(streamID1, single1, &varnmiss[varID][levelID]);
	  vars[varID] = TRUE;
	}

      /* Perform the interpolation on all valid data variables */
      for ( varID = 0; varID < nvars; varID++ )
	{
	  if ( vars[varID] && varinterp[varID] )
	    {
	      gridID   = vlistInqVarGrid(vlistID1, varID);
	      missval  = vlistInqVarMissval(vlistID1, varID);
	      gridsize = gridInqSize(gridID);

	      vert_interp_lev3d(gridsize, missval, vardata1[varID], vardata2[varID],
				nlevo, lev_idx1, lev_idx2, lev_wgt1, lev_wgt2);

	      for ( levelID = 0; levelID < nlevo; levelID++ )
		{
		  gridsize = gridInqSize(vlistInqVarGrid(vlistID3, varID));
		  offset   = gridsize*levelID;
		  single2  = vardata2[varID] + offset;
		  nmiss    = 0;
		  for ( i = 0; i < gridsize; ++i )
		    if ( DBL_IS_EQUAL(single2[i], missval) ) nmiss++;
		  varnmiss[varID][levelID] = nmiss;
		}
	    }
          else
          { 
            vlistInqVarName(vlistID1, varID, varname); 
            if ( cdoVerbose ) cdoPrint("Perform no interpolation on variable %s",varname);
          }
	}

      /* write the output */
      for ( varID = 0; varID < nvars; varID++ )
	{
	  if ( vars[varID] )
	    {
	      nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID3, varID));
	      for ( levelID = 0; levelID < nlevel; levelID++ )
		{
		  gridsize = gridInqSize(vlistInqVarGrid(vlistID3, varID));
		  offset   = gridsize*levelID;

		  single2  = vardata2[varID] + offset;
		  streamDefRecord(streamID3, varID, levelID);
		  streamWriteRecord(streamID3, single2, varnmiss[varID][levelID]);

		}
	    }
	}

      /* copy output z coordinate to output stream */
      nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID3, oz3dvarID));
      for ( levelID = 0; levelID < nlevel; levelID++ )
        {
          gridsize = gridInqSize(vlistInqVarGrid(vlistID3, oz3dvarID));
          offset   = gridsize*levelID;
          single2  = zlevels_out + offset;
          streamDefRecord(streamID3, oz3dvarID, levelID);
          streamWriteRecord(streamID3, single2, 0);

        }
      tsID++;
    }

  nvars     = vlistNvars(vlistID1);
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

  streamClose(streamID0);
  streamClose(streamID1);
  streamClose(streamID2);
  streamClose(streamID3);

  vlistDestroy(vlistID3);

  cdoFinish();

  return (0);
}
