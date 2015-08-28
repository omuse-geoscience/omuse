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

      Vertint    ml2pl           Model to pressure level interpolation
      Vertint    ml2hl           Model to height level interpolation
*/


#include <ctype.h>

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "after_vertint.h"
#include "list.h"
#include "stdnametable.h"
#include "constants.h"


void *Vertintml(void *argument)
{
  int mode;
  enum {ECHAM_MODE, WMO_MODE};
  enum {func_pl, func_hl};
  enum {type_lin, type_log};
  int gridsize;
  int recID, nrecs;
  int i, k, offset;
  int tsID, varID, levelID;
  int zaxisIDp, zaxisIDh = -1, nzaxis;
  int gridID, zaxisID;
  int nhlev = 0, nhlevf = 0, nhlevh = 0, nlevel;
  int *vert_index = NULL;
  int nvct;
  int sgeopot_needed = FALSE;
  int sgeopotID = -1, geopotID = -1, tempID = -1, psID = -1, lnpsID = -1, presID = -1, gheightID = -1;
  int code, param;
  int pnum, pcat, pdis;
  //int sortlevels = TRUE;
  int *pnmiss = NULL;
  char paramstr[32];
  char varname[CDI_MAX_NAME], stdname[CDI_MAX_NAME];
  double minval, maxval;
  double missval;
  double *vct = NULL;
  double *rvct = NULL; /* reduced VCT for LM */
  double *single1, *single2;
  double *sgeopot = NULL, *ps_prog = NULL, *full_press = NULL, *half_press = NULL;
  double *hyb_press = NULL;
  int Extrapolate = 0;
  int lhavevct;
  int mono_level;
  int instNum, tableNum;
  int useTable;
  gribcode_t gribcodes = {0};
  LIST *flist = listNew(FLT_LIST);

  cdoInitialize(argument);

  int ML2PL     = cdoOperatorAdd("ml2pl",     func_pl, type_lin, "pressure levels in pascal");
  int ML2PLX    = cdoOperatorAdd("ml2plx",    func_pl, type_lin, "pressure levels in pascal");
  int ML2HL     = cdoOperatorAdd("ml2hl",     func_hl, type_lin, "height levels in meter");
  int ML2HLX    = cdoOperatorAdd("ml2hlx",    func_hl, type_lin, "height levels in meter");
  int ML2PL_LP  = cdoOperatorAdd("ml2pl_lp",  func_pl, type_log, "pressure levels in pascal");
  int ML2PLX_LP = cdoOperatorAdd("ml2plx_lp", func_pl, type_log, "pressure levels in pascal");
  int ML2HL_LP  = cdoOperatorAdd("ml2hl_lp",  func_hl, type_log, "height levels in meter");
  int ML2HLX_LP = cdoOperatorAdd("ml2hlx_lp", func_hl, type_log, "height levels in meter");

  int operatorID = cdoOperatorID();
  int operfunc   = cdoOperatorF1(operatorID);
  int opertype   = cdoOperatorF2(operatorID);

  if ( operatorID == ML2PL || operatorID == ML2HL || operatorID == ML2PL_LP || operatorID == ML2HL_LP )
    {
      char *envstr;
      envstr = getenv("EXTRAPOLATE");

      if ( envstr )
	{
	  if ( isdigit((int) envstr[0]) )
	    {
	      Extrapolate = atoi(envstr);
	      if ( Extrapolate == 1 )
		cdoPrint("Extrapolation of missing values enabled!");
	    }
	}
    }
  else if ( operatorID == ML2PLX || operatorID == ML2HLX || operatorID == ML2PLX_LP || operatorID == ML2HLX_LP )
    {
      Extrapolate = 1;
    }

  operatorInputArg(cdoOperatorEnter(operatorID));

  int nplev = 0;
  double *plev = NULL;
  if ( operatorArgc() == 1 && strcmp(operatorArgv()[0], "default") == 0 )
    {
      /*
      double stdlev[] = {100000, 92500, 85000, 77500, 70000, 60000, 50000, 40000, 30000, 25000, 20000,
                          15000, 10000, 7000, 5000, 3000, 2000, 1000, 700, 500, 300, 200, 100, 50, 20, 10};
      */
      double stdlev[] = {100000, 92500, 85000, 70000, 60000, 50000, 40000, 30000, 25000, 20000, 15000,
                          10000,  7000,  5000,  3000,  2000, 1000 };
        nplev = sizeof(stdlev)/sizeof(*stdlev);
      plev  = (double *) malloc(nplev*sizeof(double));
      for ( i = 0; i < nplev; ++i ) plev[i] = stdlev[i];
    }
  else
    {
      nplev = args2fltlist(operatorArgc(), operatorArgv(), flist);
      plev  = (double *) listArrayPtr(flist);
    }
  
  int streamID1 = streamOpenRead(cdoStreamName(0));

  int vlistID1 = streamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  gridsize = vlist_check_gridsize(vlistID1);

  if ( operfunc == func_hl )
    zaxisIDp = zaxisCreate(ZAXIS_HEIGHT, nplev);
  else
    zaxisIDp = zaxisCreate(ZAXIS_PRESSURE, nplev);

  zaxisDefLevels(zaxisIDp, plev);
  nzaxis  = vlistNzaxis(vlistID1);
  lhavevct = FALSE;
  for ( i = 0; i < nzaxis; i++ )
    {
      /* mono_level = FALSE; */
      mono_level = TRUE;
      zaxisID = vlistZaxis(vlistID1, i);
      nlevel  = zaxisInqSize(zaxisID);

      if ( (zaxisInqType(zaxisID) == ZAXIS_HYBRID || zaxisInqType(zaxisID) == ZAXIS_HYBRID_HALF) &&
	   nlevel > 1 )
	{
	  double *level;
	  int l;
	  level = (double*) malloc(nlevel*sizeof(double));
	  zaxisInqLevels(zaxisID, level);
	  for ( l = 0; l < nlevel; l++ )
	    {
	      if ( (l+1) != (int) (level[l]+0.5) ) break;
	    }
	  if ( l == nlevel ) mono_level = TRUE; 
	  free(level);
	}

      if ( (zaxisInqType(zaxisID) == ZAXIS_HYBRID || zaxisInqType(zaxisID) == ZAXIS_HYBRID_HALF) &&
	   nlevel > 1 && mono_level )
	{
	  nvct = zaxisInqVctSize(zaxisID);
	  if ( nlevel == (nvct/2 - 1) )
	    {
	      if ( lhavevct == FALSE )
		{
		  lhavevct = TRUE;
		  zaxisIDh = zaxisID;
		  nhlev    = nlevel;
		  nhlevf   = nhlev;
		  nhlevh   = nhlevf + 1;
	      
		  vct = (double*) malloc(nvct*sizeof(double));
		  zaxisInqVct(zaxisID, vct);

		  vlistChangeZaxisIndex(vlistID2, i, zaxisIDp);
		}
	      else
		{
		  if ( memcmp(vct, zaxisInqVctPtr(zaxisID), nvct*sizeof(double)) == 0 )
		    vlistChangeZaxisIndex(vlistID2, i, zaxisIDp);
		}
	    }
	  else if ( nlevel == (nvct/2) )
	    {
	      if ( lhavevct == FALSE )
		{
		  lhavevct = TRUE;
		  zaxisIDh = zaxisID;
		  nhlev    = nlevel;
		  nhlevf   = nhlev - 1;
		  nhlevh   = nhlev;
	      
		  vct = (double*) malloc(nvct*sizeof(double));
		  zaxisInqVct(zaxisID, vct);

		  vlistChangeZaxisIndex(vlistID2, i, zaxisIDp);
		}
	      else
		{
		  if ( memcmp(vct, zaxisInqVctPtr(zaxisID), nvct*sizeof(double)) == 0 )
		    vlistChangeZaxisIndex(vlistID2, i, zaxisIDp);
		}
	    }
	  else if ( nlevel == (nvct - 4 - 1) )
	    {
	      if ( lhavevct == FALSE )
		{
		  int vctsize;
		  int voff = 4;
		  
		  rvct = (double*) malloc(nvct*sizeof(double));
		  zaxisInqVct(zaxisID, rvct);

		  if ( (int)(rvct[0]+0.5) == 100000 && rvct[voff] < rvct[voff+1] )
		    {
		      lhavevct = TRUE;
		      zaxisIDh = zaxisID;
		      nhlev    = nlevel;
		      nhlevf   = nhlev;
		      nhlevh   = nhlev + 1;

		      vctsize = 2*nhlevh;
		      vct = (double*) malloc(vctsize*sizeof(double));

		      vlistChangeZaxisIndex(vlistID2, i, zaxisIDp);

		      /* calculate VCT for LM */

		      for ( i = 0; i < vctsize/2; i++ )
			{
			  if ( rvct[voff+i] >= rvct[voff] && rvct[voff+i] <= rvct[3] )
			    {
			      vct[i] = rvct[0]*rvct[voff+i];
			      vct[vctsize/2+i] = 0;
			    }
			  else
			    {
			      vct[i] = (rvct[0]*rvct[3]*(1-rvct[voff+i]))/(1-rvct[3]);
			      vct[vctsize/2+i] = (rvct[voff+i]-rvct[3])/(1-rvct[3]);
			    }
			}
		      
		      if ( cdoVerbose )
			{
			  for ( i = 0; i < vctsize/2; i++ )
			    fprintf(stdout, "%5d %25.17f %25.17f\n", i, vct[i], vct[vctsize/2+i]);
			}
		    }
		}
	      else
		{
		  if ( memcmp(rvct, zaxisInqVctPtr(zaxisID), nvct*sizeof(double)) == 0 )
		    vlistChangeZaxisIndex(vlistID2, i, zaxisIDp);
		}
	    }
	}
    }

  int nvars = vlistNvars(vlistID1);

  int vars[nvars];
  double *vardata1[nvars];
  double *vardata2[nvars];
  int *varnmiss[nvars];
  int varinterp[nvars];

  int maxlev   = nhlevh > nplev ? nhlevh : nplev;

  if ( Extrapolate == 0 )
    pnmiss = (int*) malloc(nplev*sizeof(int));

  // check levels
  if ( zaxisIDh != -1 )
    {
      int nlev = zaxisInqSize(zaxisIDh);
      if ( nlev != nhlev ) cdoAbort("Internal error, wrong numner of hybrid level!");
      double levels[nlev];
      zaxisInqLevels(zaxisIDh, levels);

      for ( int ilev = 0; ilev < nlev; ++ilev )
	{
	  if ( (ilev+1) != (int)levels[ilev] )
	    {
	      //sortlevels = FALSE;
	      break;
	    }
	}
    }

  if ( zaxisIDh != -1 && gridsize > 0 )
    {
      vert_index = (int*) malloc(gridsize*nplev*sizeof(int));
      ps_prog    = (double*) malloc(gridsize*sizeof(double));
      full_press = (double*) malloc(gridsize*nhlevf*sizeof(double));
      half_press = (double*) malloc(gridsize*nhlevh*sizeof(double));
    }
  else
    cdoWarning("No 3D variable with hybrid sigma pressure coordinate found!");

  if ( operfunc == func_hl )
    {
      double phlev[nplev];
      height2pressure(phlev, plev, nplev);

      if ( cdoVerbose )
	for ( i = 0; i < nplev; ++i )
	  cdoPrint("level = %d   height = %g   pressure = %g", i+1, plev[i], phlev[i]);

      memcpy(plev, phlev, nplev*sizeof(double));
    }

  if ( opertype == type_log )
    for ( k = 0; k < nplev; k++ ) plev[k] = log(plev[k]);

  useTable = FALSE;
  for ( varID = 0; varID < nvars; varID++ )
    {
      tableNum = tableInqNum(vlistInqVarTable(vlistID1, varID));
      if ( tableNum > 0 && tableNum != 255 )
	{
	  useTable = TRUE;
	  break;
	}
    }

  if ( cdoVerbose && useTable ) cdoPrint("Using code tables!");

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridID   = vlistInqVarGrid(vlistID1, varID);
      zaxisID  = vlistInqVarZaxis(vlistID1, varID);
      // gridsize = gridInqSize(gridID);
      nlevel   = zaxisInqSize(zaxisID);
      instNum  = institutInqCenter(vlistInqVarInstitut(vlistID1, varID));
      tableNum = tableInqNum(vlistInqVarTable(vlistID1, varID));

      code     = vlistInqVarCode(vlistID1, varID);
      param    = vlistInqVarParam(vlistID1, varID);

      cdiParamToString(param, paramstr, sizeof(paramstr));
      cdiDecodeParam(param, &pnum, &pcat, &pdis);
      if ( pdis >= 0 && pdis < 255 ) code = -1;

      if ( useTable )
	{
	  if ( tableNum == 2 )
	    {
	      mode = WMO_MODE;
	      wmo_gribcodes(&gribcodes);
	    }
	  else if ( tableNum == 128 || tableNum == 0 )
	    {
	      mode = ECHAM_MODE;
	      echam_gribcodes(&gribcodes);
	    }
	  else
	    mode = -1;
	}
      else
	{
	  mode = ECHAM_MODE;
	  echam_gribcodes(&gribcodes);
	}

      if ( cdoVerbose )
	cdoPrint("Mode = %d  Center = %d  Code = %d  Param = %s", mode, instNum, code, paramstr);

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
	      if      ( sgeopotID == -1 && (strcmp(varname, "geosp") == 0 || strcmp(varname, "z")    == 0) ) code = gribcodes.geopot;
	      else if ( tempID    == -1 && (strcmp(varname, "st")    == 0 || strcmp(varname, "t")    == 0) ) code = gribcodes.temp;
	      else if ( psID      == -1 && (strcmp(varname, "aps")   == 0 || strcmp(varname, "sp"  ) == 0) ) code = gribcodes.ps;
	      else if ( lnpsID    == -1 && (strcmp(varname, "lsp")   == 0 || strcmp(varname, "lnsp") == 0) ) code = gribcodes.lsp;
	      else if ( geopotID  == -1 && strcmp(stdname, "geopotential_full") == 0 ) code = gribcodes.geopot;
	      /* else if ( strcmp(varname, "geopoth") == 0 ) code = 156; */
	    }
	}

      if ( mode == ECHAM_MODE )
	{
	  if      ( code == gribcodes.geopot  && nlevel == 1      ) sgeopotID = varID;
	  else if ( code == gribcodes.geopot  && nlevel == nhlevf ) geopotID  = varID;
	  else if ( code == gribcodes.temp    && nlevel == nhlevf ) tempID    = varID;
	  else if ( code == gribcodes.ps      && nlevel == 1      ) psID      = varID;
	  else if ( code == gribcodes.lsp     && nlevel == 1      ) lnpsID    = varID;
	  else if ( code == gribcodes.gheight && nlevel == nhlevf ) gheightID = varID;
	}
      else if ( mode == WMO_MODE )
	{
	  if      ( code == gribcodes.geopot  && nlevel == 1      ) sgeopotID = varID;
	  else if ( code == gribcodes.geopot  && nlevel == nhlevf ) geopotID  = varID;
	  else if ( code == gribcodes.temp    && nlevel == nhlevf ) tempID    = varID;
	  else if ( code == gribcodes.ps      && nlevel == 1      ) psID      = varID;
	}

      if ( gridInqType(gridID) == GRID_SPECTRAL && zaxisInqType(zaxisID) == ZAXIS_HYBRID )
	cdoAbort("Spectral data on model level unsupported!");

      if ( gridInqType(gridID) == GRID_SPECTRAL )
	cdoAbort("Spectral data unsupported!");

      if ( varID == gheightID )
	vardata1[varID] = (double*) malloc(gridsize*(nlevel+1)*sizeof(double));
      else
	vardata1[varID] = (double*) malloc(gridsize*nlevel*sizeof(double));

      /* if ( zaxisInqType(zaxisID) == ZAXIS_HYBRID && zaxisIDh != -1 && nlevel == nhlev ) */
      if ( zaxisID == zaxisIDh ||
	   (zaxisInqType(zaxisID) == ZAXIS_HYBRID && zaxisIDh != -1 && (nlevel == nhlevh || nlevel == nhlevf)) )
	{
	  varinterp[varID] = TRUE;
	  vardata2[varID]  = (double*) malloc(gridsize*nplev*sizeof(double));
	  varnmiss[varID]  = (int*) malloc(maxlev*sizeof(int));
	  memset(varnmiss[varID], 0, maxlev*sizeof(int));
	}
      else
	{
	  if ( zaxisInqType(zaxisID) == ZAXIS_HYBRID && zaxisIDh != -1 && nlevel > 1 )
	    cdoWarning("Parameter %d has wrong number of levels, skipped! (param=%s nlevel=%d)",
		       varID+1, paramstr, nlevel);
	  varinterp[varID] = FALSE;
	  vardata2[varID]  = vardata1[varID];
	  varnmiss[varID]  = (int*) malloc(nlevel*sizeof(int));
	}
    }

  if ( cdoVerbose )
    {
      cdoPrint("Found:");
      if ( tempID    != -1 ) cdoPrint("  %s", var_stdname(air_temperature));
      if ( psID      != -1 ) cdoPrint("  %s", var_stdname(surface_air_pressure));
      if ( lnpsID    != -1 ) cdoPrint("  LOG(%s)", var_stdname(surface_air_pressure));
      if ( sgeopotID != -1 ) cdoPrint("  %s", var_stdname(surface_geopotential));
      if ( geopotID  != -1 ) cdoPrint("  %s", var_stdname(geopotential));
      if ( gheightID != -1 ) cdoPrint("  %s", var_stdname(geopotential_height));
    }

  if ( tempID != -1 || gheightID != -1 ) sgeopot_needed = TRUE;

  if ( zaxisIDh != -1 && sgeopot_needed )
    {
      sgeopot = (double*) malloc(gridsize*sizeof(double));
      if ( sgeopotID == -1 )
	{
	  if ( geopotID == -1 )
	    cdoWarning("%s not found - set to zero!", var_stdname(surface_geopotential));
	  else
	    cdoPrint("%s not found - using bottom layer of %s!", var_stdname(surface_geopotential), var_stdname(geopotential));

	  memset(sgeopot, 0, gridsize*sizeof(double));
	}
    }

  if ( zaxisIDh != -1 && gheightID != -1 && tempID == -1 )
    cdoAbort("Temperature not found, needed for vertical interpolation of geopotheight!");

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

  int streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      for ( varID = 0; varID < nvars; ++varID ) vars[varID] = FALSE;

      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  //gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	  zaxisID  = vlistInqVarZaxis(vlistID1, varID);
	  nlevel   = zaxisInqSize(zaxisID);
	  /*
	  if ( sortlevels && zaxisIDh != -1 && zaxisID == zaxisIDh && nlevel == nhlev )
	    {
	      levelID = (int) (zaxisInqLevel(zaxisIDh, levelID)-1);
	      printf("levelID %d\n", levelID);
	    }
	  */
	  offset   = gridsize*levelID;
	  single1  = vardata1[varID] + offset;

	  streamReadRecord(streamID1, single1, &varnmiss[varID][levelID]);
	  vars[varID] = TRUE;
	}

      if ( zaxisIDh != -1 )
	{
	  if ( sgeopot_needed )
	    {
	      if ( sgeopotID != -1 )
		memcpy(sgeopot, vardata1[sgeopotID], gridsize*sizeof(double));
	      else if ( geopotID != -1 )
		memcpy(sgeopot, vardata1[geopotID]+gridsize*(nhlevf-1), gridsize*sizeof(double));

	      /* check range of surface geopot */
	      if ( sgeopotID != -1 || geopotID != -1 )
		{
		  minmaxval(gridsize, sgeopot, NULL, &minval, &maxval);
		  if ( minval < MIN_FIS || maxval > MAX_FIS )
		    cdoWarning("Surface geopotential out of range (min=%g max=%g) [timestep:%d]!", minval, maxval, tsID+1);
		  if ( gridsize > 1 && minval >= 0 && maxval <= 9000 )
		    cdoWarning("Surface geopotential has an unexpected range (min=%g max=%g) [timestep:%d]!", minval, maxval, tsID+1);
		}
	    }

	  if ( lnpsID != -1 )
	    for ( i = 0; i < gridsize; i++ ) ps_prog[i] = exp(vardata1[lnpsID][i]);
	  else if ( psID != -1 )
	    memcpy(ps_prog, vardata1[psID], gridsize*sizeof(double));

	  /* check range of ps_prog */
	  minmaxval(gridsize, ps_prog, NULL, &minval, &maxval);
	  if ( minval < MIN_PS || maxval > MAX_PS )
	    cdoWarning("Surface pressure out of range (min=%g max=%g) [timestep:%d]!", minval, maxval, tsID+1);


	  presh(full_press, half_press, vct, ps_prog, nhlevf, gridsize);

	  if ( opertype == type_log )
	    {
	      for ( i = 0; i < gridsize; i++ ) ps_prog[i] = log(ps_prog[i]);

	      for ( k = 0; k < nhlevh; k++ )
		for ( i = 0; i < gridsize; i++ )
		  half_press[k*gridsize+i] = log(half_press[k*gridsize+i]);

	      for ( k = 0; k < nhlevf; k++ )
		for ( i = 0; i < gridsize; i++ )
		  full_press[k*gridsize+i] = log(full_press[k*gridsize+i]);
	    }

	  genind(vert_index, plev, full_press, gridsize, nplev, nhlevf);

	  if ( Extrapolate == 0 )
	    genindmiss(vert_index, plev, gridsize, nplev, ps_prog, pnmiss);
	}

      for ( varID = 0; varID < nvars; varID++ )
	{
	  if ( vars[varID] )
	    {
	      gridID   = vlistInqVarGrid(vlistID1, varID);
	      zaxisID  = vlistInqVarZaxis(vlistID1, varID);
	      missval  = vlistInqVarMissval(vlistID1, varID);
	      //gridsize = gridInqSize(gridID);
	      nlevel   = zaxisInqSize(zaxisID);
	      if ( varinterp[varID] )
		{
		  /*
		  if ( nlevel == nhlevh )
		    {
		      int i, k;
		      double *vl1, *vl2;

		      for ( k = 1; k < nlevel; k++ )
			{
			  vl1  = vardata1[varID] + gridsize*(k-1);
			  vl2  = vardata1[varID] + gridsize*(k);
			  for ( i = 0; i < gridsize; i++ )
			    vl1[i] = 0.5*(vl1[i] + vl2[i]);
			}
		      
		      nlevel = nhlevf;
		    }
		  */
		  if ( nlevel == nhlevh )
		    {
		      hyb_press = half_press;
		    }
		  else if ( nlevel == nhlevf )
		    {
		      hyb_press = full_press;
		    }
		  else
		    {
		      param = vlistInqVarParam(vlistID1, varID);
		      cdiParamToString(param, paramstr, sizeof(paramstr));
		      cdoAbort("Number of hybrid level differ from full/half level (param=%s)!", paramstr);
		    }

		  for ( levelID = 0; levelID < nlevel; levelID++ )
		    {
		      if ( varnmiss[varID][levelID] )
			cdoAbort("Missing values unsupported for this operator!");
		    }

		  if ( varID == tempID )
		    {
		      if ( nlevel == nhlevh )
			cdoAbort("Temperature on half level unsupported!");

		      if ( opertype == type_log && Extrapolate )
			cdoAbort("Log. extrapolation of temperature unsupported!");

		      interp_T(sgeopot, vardata1[varID], vardata2[varID],
			       full_press, half_press, vert_index,
			       plev, nplev, gridsize, nlevel, missval);
		    }
		  else if ( varID == gheightID )
		    {
		      for ( i = 0; i < gridsize; ++i )
			vardata1[varID][gridsize*nlevel+i] = sgeopot[i]/PlanetGrav;

		      interp_Z(sgeopot, vardata1[varID], vardata2[varID],
			       full_press, half_press, vert_index, vardata1[tempID],
			       plev, nplev, gridsize, nlevel, missval);
		    }
		  else
		    {
		      interp_X(vardata1[varID], vardata2[varID], hyb_press,
			       vert_index, plev, nplev, gridsize, nlevel, missval);
		    }
		  
		  if ( Extrapolate == 0 )
		    memcpy(varnmiss[varID], pnmiss, nplev*sizeof(int));
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
		  //gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
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

  if ( pnmiss     ) free(pnmiss);

  if ( sgeopot    ) free(sgeopot);
  if ( ps_prog    ) free(ps_prog);
  if ( vert_index ) free(vert_index);
  if ( full_press ) free(full_press);
  if ( half_press ) free(half_press);
  if ( vct        ) free(vct);
  if ( rvct       ) free(rvct);

  listDelete(flist);

  cdoFinish();

  return (0);
}
