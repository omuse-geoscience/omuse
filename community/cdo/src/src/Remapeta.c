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

      Remapeta     remapeta          Model to model level interpolation
*/

#include <ctype.h>

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "after_vertint.h"
#include "list.h"
#include "stdnametable.h"
#include "hetaeta.h"


static 
void setmissval(long nvals, int *imiss, double missval, double *array)
{
  long i;

  if ( imiss )
    {
      for ( i = 0; i < nvals; ++i )
	if ( imiss[i] ) array[i] = missval;
    }
}

static
void corr_hum(long gridsize, double *q, double q_min)
{
  long i;

  for ( i = 0; i < gridsize; ++i )
    {
      if ( q[i] < q_min ) q[i] = q_min;
    }
}
 
static
long ncctop(double cptop, long nlev, long nlevp1, double *vct_a, double *vct_b)
{
  /*
    Description:
    Defines highest level *ncctop* where condensation is allowed.
    
    Author:
  
    E. Roeckner, MPI, October 2001
  */
  /* local variables */
  long nctop = 0;
  long jk;
  double    za, zb, zph[nlevp1], zp[nlev];
  // double    cptop  =  1000.;   /* min. pressure level for cond. */

  /* half level pressure values, assuming 101320. Pa surface pressure */

  for ( jk = 0; jk < nlevp1; ++jk )
    {
      za = vct_a[jk];
      zb = vct_b[jk];
      zph[jk] = za + zb*101320.;
    }

  /* full level pressure */

  for ( jk = 0; jk < nlev; ++jk )
    zp[jk] = (zph[jk] + zph[jk+1])*0.5;

  /* search for pressure level cptop (Pa) */

  for ( jk = 0; jk < nlev; ++jk )
    {
      nctop = jk;
      if ( zp[jk] >= cptop ) break;
    }

  return (nctop);
}


double *vctFromFile(const char *filename, int *nvct)
{
  char line[1024], *pline;
  int num, i = 0;
  int nlevh2, nvct2;
  int maxvct = 8192;
  double *vct2;
  FILE *fp;


  fp = fopen(filename, "r");
  if ( fp == NULL ) { perror(filename); exit(EXIT_FAILURE); }

  vct2 = (double*) malloc(maxvct*sizeof(double));

  while ( readline(fp, line, 1024) )
    {
      if ( line[0] == '#' ) continue;
      if ( line[0] == '\0' ) continue;

      pline = line;
      num = (int) strtod(pline, &pline);
      if ( pline == NULL ) cdoAbort("Format error in VCT file %s!", filename);
      if ( num != i ) cdoWarning("Inconsistent VCT file, entry %d is %d.", i, num);
      
      if ( i+maxvct/2 >= maxvct-1 ) cdoAbort("Too many values in VCT file!");

      vct2[i] = strtod(pline, &pline);
      if ( pline == NULL ) cdoAbort("Format error in VCT file %s!", filename);

      vct2[i+maxvct/2] = strtod(pline, &pline);

      i++;
    }

  fclose(fp);

  nvct2 = 2*i;
  nlevh2 = i - 1;

  for ( i = 0; i < nlevh2+1; ++i )
    vct2[i+nvct2/2] = vct2[i+maxvct/2];
  
  vct2 = (double*) realloc(vct2, nvct2*sizeof(double));

  *nvct = nvct2;

  return (vct2);
}

static
void vert_sum(double *sum, double *var3d, long gridsize, long nlevel)
{
  long i, k;

  for ( i = 0; i < gridsize; ++i ) sum[i] = 0;

  for ( k = 0; k < nlevel; ++k )
    for ( i = 0; i < gridsize; ++i )
      {
	sum[i] += var3d[k*gridsize + i];
      }
}

static
void vert_sumw(double *sum, double *var3d, long gridsize, long nlevel, double *deltap)
{
  long i, k;

  for ( i = 0; i < gridsize; ++i ) sum[i] = 0;

  for ( k = 0; k < nlevel; ++k )
    for ( i = 0; i < gridsize; ++i )
      {
	sum[i] += var3d[k*gridsize + i]*deltap[k*gridsize + i];
      }
}


#define  MAX_VARS3D  1024

void *Remapeta(void *argument)
{
  int REMAPETA, REMAPETAS, REMAPETAZ;
  int operatorID;
  int streamID1, streamID2;
  int vlistID1, vlistID2;
  int nfis2gp = 0;
  int recID, nrecs;
  int i, offset, iv;
  int tsID, varID, levelID;
  int nvars, nvars3D = 0;
  int zaxisID2, zaxisIDh = -1, nzaxis, surfaceID;
  int gridID, zaxisID;
  int nlevel;
  int nvct1, nvct2 = 0;
  int sgeopotID = -1, tempID = -1, sqID = -1, psID = -1, lnpsID = -1, presID = -1;
  int code;
  char varname[CDI_MAX_NAME], stdname[CDI_MAX_NAME];
  double *single2;
  int taxisID1, taxisID2;
  int lhavevct;
  int nhlevf1 = 0, nhlevf2 = 0;
  double *lev2;
  double *vct1 = NULL, *vct2 = NULL;
  double *a1 = NULL, *b1 = NULL, *a2 = NULL, *b2 = NULL;
  double *fis1 = NULL, *ps1 = NULL, *t1 = NULL, *q1 = NULL;
  double *fis2 = NULL, *ps2 = NULL, *t2 = NULL, *q2 = NULL;
  double *tscor = NULL, *pscor = NULL, *secor = NULL;
  int nmiss, nmissout = 0;
  int ltq = FALSE;
  int lfis2 = FALSE;
  int varids[MAX_VARS3D];
  int *imiss = NULL;
  int timer_hetaeta = 0;
  long nctop = 0;
  double *array = NULL;
  double *deltap1 = NULL, *deltap2 = NULL;
  double *half_press1 = NULL, *half_press2 = NULL;
  double *sum1 = NULL, *sum2 = NULL;
  double **vars1 = NULL, **vars2 = NULL;
  double minval, maxval;
  double missval = 0;
  double cconst = 1.E-6;
  const char *fname;
  char *envstr;
  double cptop  = 0; /* min. pressure level for cond. */

  if ( cdoTimer ) timer_hetaeta = timer_new("Remapeta_hetaeta");

  cdoInitialize(argument);

  REMAPETA  = cdoOperatorAdd("remapeta",   0, 0, "VCT file name");
  REMAPETAS = cdoOperatorAdd("remapeta_s", 0, 0, "VCT file name");
  REMAPETAZ = cdoOperatorAdd("remapeta_z", 0, 0, "VCT file name");

  operatorID = cdoOperatorID();

  operatorInputArg(cdoOperatorEnter(operatorID));

  envstr = getenv("REMAPETA_PTOP");
  if ( envstr )
    {
      double fval = atof(envstr);
      if ( fval > 0 )
	{
	  cptop = fval;
	  //	  if ( cdoVerbose )
	  cdoPrint("Set REMAPETA_PTOP to %g", cptop);
	}
    }  

  vct2 = vctFromFile(operatorArgv()[0], &nvct2);
  nhlevf2 = nvct2/2 - 1;

  a2 = vct2;
  b2 = vct2 + nvct2/2;

  if ( cdoVerbose )
    for ( i = 0; i < nhlevf2+1; ++i )
      cdoPrint("vct2: %5d %25.17f %25.17f", i, vct2[i], vct2[nvct2/2+i]);

  streamID1 = streamOpenRead(cdoStreamName(0));

  if ( operatorArgc() == 2 )
    {
      int streamID;

      lfis2 = TRUE;
      fname = operatorArgv()[1];
      
      argument_t *fileargument = file_argument_new(fname);
      streamID = streamOpenRead(fileargument);
      file_argument_free(fileargument);

      vlistID1 = streamInqVlist(streamID);

      streamInqRecord(streamID, &varID, &levelID);
      gridID  = vlistInqVarGrid(vlistID1, varID);
      nfis2gp = gridInqSize(gridID);

      fis2 = (double*) malloc(nfis2gp*sizeof(double));

      streamReadRecord(streamID, fis2, &nmiss);

      if ( nmiss )
	{
	  missval = vlistInqVarMissval(vlistID1, varID);
	  imiss = (int*) malloc (nfis2gp*sizeof(int));
	  for ( i = 0; i < nfis2gp; ++i )
	    {
	      if ( DBL_IS_EQUAL(fis2[i], missval) )
		imiss[i] = 1;
	      else 
		imiss[i] = 0;
	    }

	  nmissout = nmiss;
	}

      /* check range of surface_geopotential */
      minmaxval(nfis2gp, fis2, imiss, &minval, &maxval);
      if ( minval < MIN_FIS || maxval > MAX_FIS )
	cdoWarning("%s out of range (min=%g max=%g)!", var_stdname(surface_geopotential), minval, maxval);

      if ( minval < -1.e10 || maxval > 1.e10 )
	cdoAbort("%s out of range!", var_stdname(surface_geopotential));

      streamClose(streamID); 
    }

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  gridID = vlistGrid(vlistID1, 0);
  if ( gridInqType(gridID) == GRID_SPECTRAL )
    cdoAbort("Spectral data unsupported!");

  int gridsize = vlist_check_gridsize(vlistID1);

  zaxisID2 = zaxisCreate(ZAXIS_HYBRID, nhlevf2);
  lev2 = (double*) malloc(nhlevf2*sizeof(double));
  for ( i = 0; i < nhlevf2; ++i ) lev2[i] = i+1;
  zaxisDefLevels(zaxisID2, lev2);
  free(lev2);

  if ( nvct2 == 0 ) cdoAbort("Internal problem, vct2 undefined!");
  zaxisDefVct(zaxisID2, nvct2, vct2);

  surfaceID = zaxisFromName("surface");

  nzaxis  = vlistNzaxis(vlistID1);
  lhavevct = FALSE;

  if ( cdoVerbose )
    cdoPrint("nzaxis: %d", nzaxis);

  for ( i = 0; i < nzaxis; i++ )
    {
      zaxisID = vlistZaxis(vlistID1, i);
      nlevel  = zaxisInqSize(zaxisID);
      if ( zaxisInqType(zaxisID) == ZAXIS_HYBRID )
	{
	  if ( nlevel > 1 )
	    {
	      nvct1 = zaxisInqVctSize(zaxisID);

              if ( cdoVerbose )
                cdoPrint("i: %d, vct1 size of zaxisID %d = %d", i, zaxisID, nvct1);

	      if ( nlevel == (nvct1/2 - 1) )
		{
		  if ( lhavevct == FALSE )
		    {
		      lhavevct = TRUE;
		      zaxisIDh = zaxisID;
		      nhlevf1  = nlevel;
	      
                      if ( cdoVerbose )
                        cdoPrint("lhavevct=TRUE  zaxisIDh = %d, nhlevf1   = %d", zaxisIDh, nlevel);
 
		      vct1 = (double*) malloc(nvct1*sizeof(double));
		      zaxisInqVct(zaxisID, vct1);
		      
		      vlistChangeZaxisIndex(vlistID2, i, zaxisID2);

		      a1 = vct1;
		      b1 = vct1 + nvct1/2;
		      if ( cdoVerbose )
			for ( i = 0; i < nvct1/2; ++i )
			  cdoPrint("vct1: %5d %25.17f %25.17f", i, vct1[i], vct1[nvct1/2+i]);
		    }
		  else
		    {
		      if ( memcmp(vct1, zaxisInqVctPtr(zaxisID), nvct1*sizeof(double)) == 0 )
			vlistChangeZaxisIndex(vlistID2, i, zaxisID2);
		    }
		}
              else 
                {
		  if ( cdoVerbose )
		    cdoPrint("nlevel = (nvct1/2 - 1): nlevel = %d", nlevel);
		  if ( nlevel < (nvct1/2 - 1) )
		    cdoPrint("z-axis %d has only %d of %d hybrid sigma pressure levels!", i+1, nlevel, (nvct1/2 - 1));
                }
	    }
	  else
	    {
	      vlistChangeZaxisIndex(vlistID2, i, surfaceID);
	    }
	}
      else
        {
	  if ( cdoVerbose )
	    cdoPrint("i: %d, type of zaxisID %d not ZAXIS_HYBRID", i, zaxisID);
        }
    }

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);


  if ( zaxisIDh == -1 )
    cdoWarning("No 3D variable with hybrid sigma pressure coordinate found!");

  nvars = vlistNvars(vlistID1);

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridID  = vlistInqVarGrid(vlistID1, varID);
      zaxisID = vlistInqVarZaxis(vlistID1, varID);
      nlevel  = zaxisInqSize(zaxisID);

      code = vlistInqVarCode(vlistID1, varID);
      /* code = -1; */
      if ( code <= 0 || code == 255 )
	{
	  vlistInqVarName(vlistID1, varID, varname);
	  strtolower(varname);

	  vlistInqVarStdname(vlistID1, varID, stdname);
	  strtolower(stdname);

	  code = echamcode_from_stdname(stdname);

	  if ( code == -1 )
	    {
	      /*                                  ECHAM                            ECMWF       */
	      if      ( sgeopotID == -1 && (strcmp(varname, "geosp") == 0 || strcmp(varname, "z")    == 0) ) code = 129;
	      else if ( tempID    == -1 && (strcmp(varname, "st")    == 0 || strcmp(varname, "t")    == 0) ) code = 130;
	      else if ( psID      == -1 && (strcmp(varname, "aps")   == 0 || strcmp(varname, "ps"  ) == 0) ) code = 134;
	      else if ( lnpsID    == -1 && (strcmp(varname, "lsp")   == 0 || strcmp(varname, "lnsp") == 0) ) code = 152;
	      else if ( sqID      == -1 && (strcmp(varname, "q")     == 0 ) ) code = 133;
	    }
	}

      if      ( code == 129 && nlevel == 1       ) sgeopotID = varID;
      else if ( code == 130 && nlevel == nhlevf1 ) tempID    = varID;
      else if ( code == 133 && nlevel == nhlevf1 ) sqID      = varID;
      else if ( code == 134 && nlevel == 1       ) psID      = varID;
      else if ( code == 152 && nlevel == 1       ) lnpsID    = varID;

      if ( gridInqType(gridID) == GRID_SPECTRAL && zaxisInqType(zaxisID) == ZAXIS_HYBRID )
	cdoAbort("Spectral data on model level unsupported!");

      if ( gridInqType(gridID) == GRID_SPECTRAL )
	cdoAbort("Spectral data unsupported!");

      if ( zaxisInqType(zaxisID) == ZAXIS_HYBRID && zaxisIDh != -1 && nlevel == nhlevf1 )
	{
	  if ( ! (code == 130 || code == 133) )
	    varids[nvars3D++] = varID;
	}
      else
	{
	  if ( code == 130 ) tempID = -1;
	  if ( code == 133 ) sqID   = -1;
	}
    }

  if ( cdoVerbose )
    {
      cdoPrint("Found:");
      if ( tempID    != -1 ) cdoPrint("  %s", var_stdname(air_temperature));
      if ( psID      != -1 ) cdoPrint("  %s", var_stdname(surface_air_pressure));
      if ( lnpsID    != -1 ) cdoPrint("  LOG(%s)", var_stdname(surface_air_pressure));
      if ( sgeopotID != -1 ) cdoPrint("  %s", var_stdname(surface_geopotential));
      if ( sqID      != -1 ) cdoPrint("  %s", var_stdname(specific_humidity));
    }

  if ( tempID != -1 && sqID != -1 )
    {
      ltq = TRUE;
    }
  else
    {
      if ( tempID != -1 ) cdoAbort("Temperature without humidity unsupported!");
      if ( sqID   != -1 ) cdoAbort("Humidity without temperature unsupported!");
    }
  /*
  if ( ltq == FALSE )
    {
      cdoWarning("Temperature and Humidity not found!");
    }
  */
  if ( operatorID == REMAPETA ) {}

  if ( operatorID == REMAPETAS || operatorID == REMAPETAZ)
    {
      sum1 = (double*) malloc(gridsize*sizeof(double));
      sum2 = (double*) malloc(gridsize*sizeof(double));
    }

  if ( operatorID == REMAPETAZ )
    {
      deltap1 = (double*) malloc(gridsize*nhlevf1*sizeof(double));
      deltap2 = (double*) malloc(gridsize*nhlevf2*sizeof(double));
      half_press1 = (double*) malloc(gridsize*(nhlevf1+1)*sizeof(double));
      half_press2 = (double*) malloc(gridsize*(nhlevf2+1)*sizeof(double));
    }

  array = (double*) malloc(gridsize*sizeof(double));

  fis1  = (double*) malloc(gridsize*sizeof(double));
  ps1   = (double*) malloc(gridsize*sizeof(double));

  if ( lfis2 == FALSE ) fis2  = (double*) malloc(gridsize*sizeof(double));
  if ( lfis2 == TRUE && gridsize != nfis2gp ) cdoAbort("Orographies have different grid size!");

  ps2   = (double*) malloc(gridsize*sizeof(double));

  if ( ltq )
    {
      tscor = (double*) malloc(gridsize*sizeof(double));
      pscor = (double*) malloc(gridsize*sizeof(double));
      secor = (double*) malloc(gridsize*sizeof(double));

      t1    = (double*) malloc(gridsize*nhlevf1*sizeof(double));
      q1    = (double*) malloc(gridsize*nhlevf1*sizeof(double));

      t2    = (double*) malloc(gridsize*nhlevf2*sizeof(double));
      q2    = (double*) malloc(gridsize*nhlevf2*sizeof(double));
    }

  if ( nvars3D )
    {
      vars1  = (double**) malloc(nvars*sizeof(double*));
      vars2  = (double**) malloc(nvars*sizeof(double*));

      for ( varID = 0; varID < nvars3D; ++varID )
	{
	  vars1[varID] = (double*) malloc(gridsize*nhlevf1*sizeof(double));
	  vars2[varID] = (double*) malloc(gridsize*nhlevf2*sizeof(double));
	}
    }

  if ( zaxisIDh != -1 && sgeopotID == -1 )
    {
      if ( ltq )
	cdoWarning("%s not found - set to zero!", var_stdname(surface_geopotential));

      memset(fis1, 0, gridsize*sizeof(double));
    }

  presID = lnpsID;
  if ( zaxisIDh != -1 && lnpsID == -1 )
    {
      if ( psID == -1 )
	cdoAbort("%s not found!", var_stdname(surface_air_pressure));
      else
	presID = psID;
    }

  if ( cdoVerbose )
    {
      if ( presID == lnpsID )
	cdoPrint("using LOG(%s)", var_stdname(surface_air_pressure));      
      else
	cdoPrint("using %s", var_stdname(surface_air_pressure));
    }

  if ( cdoVerbose ) cdoPrint("nvars3D = %d   ltq = %d", nvars3D, ltq);

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  zaxisID  = vlistInqVarZaxis(vlistID1, varID);
	  nlevel   = zaxisInqSize(zaxisID);
	  offset   = gridsize*levelID;
	  streamReadRecord(streamID1, array, &nmiss);

	  if ( zaxisIDh != -1 )
	    {
	      if ( varID == sgeopotID )
		memcpy(fis1, array, gridsize*sizeof(double));
	      else if ( varID == presID )
		{
		  if ( lnpsID != -1 )
		    for ( i = 0; i < gridsize; ++i ) ps1[i] = exp(array[i]);
		  else if ( psID != -1 )
		    memcpy(ps1, array, gridsize*sizeof(double));
		}
	      else if ( ltq && varID == tempID )
		memcpy(t1+offset, array, gridsize*sizeof(double));
	      else if ( ltq && varID == sqID )
		memcpy(q1+offset, array, gridsize*sizeof(double));
	      /* else if ( zaxisID == zaxisIDh ) */
	      else if ( zaxisInqType(zaxisID) == ZAXIS_HYBRID && nlevel == nhlevf1 )
		{
		  for ( i = 0; i < nvars3D; ++i )
		    if ( varID == varids[i] ) break;

		  if ( i == nvars3D ) cdoAbort("Internal error, 3D variable not found!");

		  memcpy(vars1[i]+offset, array, gridsize*sizeof(double));
		}
	      else
		{
		  streamDefRecord(streamID2, varID, levelID);
		  streamWriteRecord(streamID2, array, nmiss);
		}
	    }
	  else
	    {
	      streamDefRecord(streamID2, varID, levelID);
	      streamWriteRecord(streamID2, array, nmiss);
	    }
	}

      if ( zaxisIDh != -1 )
	{
	  /* check range of ps_prog */
	  minmaxval(gridsize, ps1, imiss, &minval, &maxval);
	  if ( minval < MIN_PS || maxval > MAX_PS )
	    cdoWarning("Surface pressure out of range (min=%g max=%g)!", minval, maxval);

	  /* check range of geop */
	  minmaxval(gridsize, fis1, imiss, &minval, &maxval);
	  if ( minval < MIN_FIS || maxval > MAX_FIS )
	    cdoWarning("Orography out of range (min=%g max=%g)!", minval, maxval);
	}

      if ( lfis2 == FALSE )
	for ( i = 0; i < gridsize; i++ ) fis2[i] = fis1[i];

      if ( ltq )
	{
	  varID = tempID;
	  nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	  for ( levelID = 0; levelID < nlevel; levelID++ )
	    {
	      offset   = gridsize*levelID;
	      single2  = t1 + offset;

	      minmaxval(gridsize, single2, imiss, &minval, &maxval);
	      if ( minval < MIN_T || maxval > MAX_T )
		cdoWarning("Input temperature at level %d out of range (min=%g max=%g)!",
			   levelID+1, minval, maxval);
	    }

	  varID = sqID;
	  nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	  for ( levelID = 0; levelID < nlevel; levelID++ )
	    {
	      offset   = gridsize*levelID;
	      single2  = q1 + offset;

	      corr_hum(gridsize, single2, MIN_Q);

	      minmaxval(gridsize, single2, imiss, &minval, &maxval);
	      if ( minval < MIN_Q || maxval > MAX_Q )
		cdoWarning("Input humidity at level %d out of range (min=%g max=%g)!",
			   levelID+1, minval, maxval);
	    }
	}

      if ( nvars3D || ltq )
	{
	  if ( cdoTimer ) timer_start(timer_hetaeta);
	  hetaeta(ltq, gridsize, imiss,
		  nhlevf1, a1, b1,
		  fis1, ps1,
		  t1, q1,
		  nhlevf2, a2, b2,
		  fis2, ps2,
		  t2, q2,
		  nvars3D, vars1, vars2,
		  tscor, pscor, secor);
	  if ( cdoTimer ) timer_stop(timer_hetaeta);
	}

      if ( cptop > 0 )
	nctop = ncctop(cptop, (long) nhlevf2, (long) nhlevf2+1, a2, b2);

      if ( zaxisIDh != -1 && sgeopotID != -1 )
	{
	  varID   = sgeopotID;
	  levelID = 0;
	  setmissval(gridsize, imiss, missval, fis2);
	  streamDefRecord(streamID2, varID, levelID);
	  streamWriteRecord(streamID2, fis2, nmissout);
	}

      if ( zaxisIDh != -1 && lnpsID != -1 )
	for ( i = 0; i < gridsize; ++i ) ps2[i] = log(ps2[i]);

      if ( zaxisIDh != -1 && presID != -1 )
	{
	  varID   = presID;
	  levelID = 0;
	  setmissval(gridsize, imiss, missval, ps2);
	  streamDefRecord(streamID2, varID, levelID);
	  streamWriteRecord(streamID2, ps2, nmissout);
	}

      if ( ltq )
	{
	  varID = tempID;
	  nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
	  for ( levelID = 0; levelID < nlevel; levelID++ )
	    {
	      offset   = gridsize*levelID;
	      single2  = t2 + offset;

	      minmaxval(gridsize, single2, imiss, &minval, &maxval);
	      if ( minval < MIN_T || maxval > MAX_T )
		cdoWarning("Output temperature at level %d out of range (min=%g max=%g)!",
			   levelID+1, minval, maxval);

	      setmissval(gridsize, imiss, missval, single2);
	      streamDefRecord(streamID2, varID, levelID);
	      streamWriteRecord(streamID2, single2, nmissout);
	    }

	  varID = sqID;
	  nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
	  for ( levelID = 0; levelID < nlevel; levelID++ )
	    {
	      offset   = gridsize*levelID;
	      single2  = q2 + offset;

	      corr_hum(gridsize, single2, MIN_Q);

	      if ( levelID < nctop )
		for ( i = 0; i < gridsize; ++i ) single2[i] = cconst;

	      minmaxval(gridsize, single2, imiss, &minval, &maxval);
	      if ( minval < MIN_Q || maxval > MAX_Q )
		cdoWarning("Output humidity at level %d out of range (min=%g max=%g)!",
			   levelID+1, minval, maxval);

	      setmissval(gridsize, imiss, missval, single2);
	      streamDefRecord(streamID2, varID, levelID);
	      streamWriteRecord(streamID2, single2, nmissout);
	    }
	}

      for ( iv = 0; iv < nvars3D; ++iv )
	{
	  varID = varids[iv];

	  nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));

	  if ( operatorID == REMAPETAS )
	    {
	      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      vert_sum(sum1, vars1[iv], gridsize, nhlevf1);

	      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
	      vert_sum(sum2, vars2[iv], gridsize, nhlevf2);
	    }
	  else if ( operatorID == REMAPETAZ )
	    {
	      int k;

	      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));

	      presh(NULL, half_press1, vct1, ps1, nhlevf1, gridsize);
	      for ( k = 0; k < nhlevf1; ++k )
		for ( i = 0; i < gridsize; ++i )
		  {
		    deltap1[k*gridsize+i] = half_press1[(k+1)*gridsize+i] - half_press1[k*gridsize+i];
		    deltap1[k*gridsize+i] = log(deltap1[k*gridsize+i]);
		  }
	      vert_sumw(sum1, vars1[iv], gridsize, nhlevf1, deltap1);

	      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));

	      presh(NULL, half_press2, vct2, ps1, nhlevf2, gridsize);
	      for ( k = 0; k < nhlevf2; ++k )
		for ( i = 0; i < gridsize; ++i )
		  {
		    deltap2[k*gridsize+i] = half_press2[(k+1)*gridsize+i] - half_press2[k*gridsize+i];
		    deltap2[k*gridsize+i] = log(deltap2[k*gridsize+i]);
		  }
	      vert_sumw(sum2, vars2[iv], gridsize, nhlevf2, deltap2);
	    }

	  for ( levelID = 0; levelID < nlevel; levelID++ )
	    {
	      offset   = gridsize*levelID;
	      single2  = vars2[iv] + offset;

	      if ( operatorID == REMAPETAS || operatorID == REMAPETAZ )
		{
		  /*
		  for ( i = 0; i < gridsize; ++i )
		    if ( i %100 == 0 )
		      printf("%d %g %g %g %g %g\n",i, single2[i], sum1[i], sum2[i], sum1[i]/sum2[i], single2[i]*sum1[i]/sum2[i]);
		  */
		  for ( i = 0; i < gridsize; ++i )
		    single2[i] = single2[i]*sum1[i]/sum2[i];
		}

	      setmissval(gridsize, imiss, missval, single2);
	      streamDefRecord(streamID2, varID, levelID);
	      streamWriteRecord(streamID2, single2, nmissout);
	    }
	}

      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  if ( nvars3D )
    {
      for ( varID = 0; varID < nvars3D; varID++ )
	{
	  free(vars2[varID]);
	  free(vars1[varID]);
	}
      free(vars2);
      free(vars1);
    }

  if ( ltq )
    {
      free(q2);
      free(t2);
      free(q1);
      free(t1);
      free(secor);
      free(pscor);
      free(tscor);
    }

  if ( imiss ) free(imiss);

  free(ps2);
  free(fis2);
  free(ps1);
  free(fis1);

  if ( sum1 ) free(sum1);
  if ( sum2 ) free(sum2);

  if ( deltap1 ) free(deltap1);
  if ( deltap2 ) free(deltap2);

  if ( half_press1 ) free(half_press1);
  if ( half_press2 ) free(half_press2);

  free(array);
  free(vct2);
  if ( vct1 ) free(vct1);

  cdoFinish();

  return (0);
}
