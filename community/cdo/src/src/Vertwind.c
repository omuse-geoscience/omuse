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

      Vertwind    vertwind      Convert the vertical velocity to [m/s]
*/

#include <ctype.h>

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "after_vertint.h"


#define R  287.07  /* spezielle Gaskonstante fuer Luft */
#define G  9.80665 /* Erdbeschleunigung */

void *Vertwind(void *argument)
{
  int streamID2;
  int vlistID2;
  int taxisID1, taxisID2;
  int gridID, zaxisID, tsID;
  int nlevel, nrecs, recID, code;
  int varID, levelID;
  int nvars, nvct = 0;
  int gridsize, i;
  int offset;
  int nmiss, nmiss_out;
  int temp_code, sq_code, ps_code, omega_code, lsp_code;
  int tempID = -1, sqID = -1, psID = -1, omegaID = -1, lnpsID = -1;
  char varname[CDI_MAX_NAME];
  double *vct = NULL;
  double tv, rho;
  double *level = NULL;
  double *temp = NULL, *sq = NULL, *omega = NULL, *wms = NULL;
  double *fpress = NULL, *hpress = NULL, *ps_prog = NULL;
  double missval_t, missval_sq, missval_wap, missval_out;

  cdoInitialize(argument);

  int streamID1 = streamOpenRead(cdoStreamName(0));

  int vlistID1 = streamInqVlist(streamID1);

  int ngp = vlist_check_gridsize(vlistID1);

  temp_code  = 130;
  sq_code    = 133;
  ps_code    = 134;
  omega_code = 135;
  lsp_code   = 152;

  nvars = vlistNvars(vlistID1);
  for ( varID = 0; varID < nvars; ++varID )
    {
      gridID   = vlistInqVarGrid(vlistID1, varID);
      zaxisID  = vlistInqVarZaxis(vlistID1, varID);
      gridsize = gridInqSize(gridID);
      nlevel   = zaxisInqSize(zaxisID);

      code = vlistInqVarCode(vlistID1, varID);

      if ( code <= 0 )
	{
	  vlistInqVarName(vlistID1, varID, varname);

	  strtolower(varname);

	  if      ( strcmp(varname, "st")    == 0 ) code = 130;
	  else if ( strcmp(varname, "sq")    == 0 ) code = 133;
	  else if ( strcmp(varname, "aps")   == 0 ) code = 134;
	  else if ( strcmp(varname, "omega") == 0 ) code = 135;
	  else if ( strcmp(varname, "lsp")   == 0 ) code = 152;
	}

      if      ( code == temp_code  ) tempID  = varID;
      else if ( code == sq_code    ) sqID    = varID;
      else if ( code == ps_code    ) psID    = varID;
      else if ( code == omega_code ) omegaID = varID;
      else if ( code == lsp_code   ) lnpsID  = varID;
    }

  if ( tempID == -1 || sqID == -1 || omegaID == -1 )
    {
      if ( tempID  == -1 ) cdoWarning("Temperature (code 130) not found!");
      if ( sqID    == -1 ) cdoWarning("Specific humidity (code 133) not found!");
      if ( omegaID == -1 ) cdoWarning("Vertical velocity (code 135) not found!");
      cdoAbort("Parameter not found!");
    }

  /* Get missing values */
  missval_t   = vlistInqVarMissval(vlistID1, tempID);
  missval_sq  = vlistInqVarMissval(vlistID1, sqID);
  missval_wap = vlistInqVarMissval(vlistID1, omegaID);
  missval_out = missval_wap;

  gridID  = vlistInqVarGrid(vlistID1, omegaID);
  zaxisID = vlistInqVarZaxis(vlistID1, omegaID);

  if ( psID == -1 && zaxisInqType(zaxisID) == ZAXIS_HYBRID )
    cdoAbort("Surface pressure (code 134) not found!");

  gridsize = gridInqSize(gridID);
  nlevel = zaxisInqSize(zaxisID);
  level  = (double*) malloc(nlevel*sizeof(double));
  zaxisInqLevels(zaxisID, level);

  temp    = (double*) malloc(gridsize*nlevel*sizeof(double));
  sq      = (double*) malloc(gridsize*nlevel*sizeof(double));
  omega   = (double*) malloc(gridsize*nlevel*sizeof(double));
  wms     = (double*) malloc(gridsize*nlevel*sizeof(double));
  fpress  = (double*) malloc(gridsize*nlevel*sizeof(double));


  if ( zaxisInqType(zaxisID) == ZAXIS_PRESSURE )
    {
      for ( levelID = 0; levelID < nlevel; ++levelID )
	{
	  offset = levelID*gridsize;
	  for ( i = 0; i < gridsize; ++i )
	    fpress[offset+i] = level[levelID];
	}
    }
  else if ( zaxisInqType(zaxisID) == ZAXIS_HYBRID )
    {
      ps_prog = (double*) malloc(gridsize*sizeof(double));
      hpress  = (double*) malloc(gridsize*(nlevel+1)*sizeof(double));
  
      nvct = zaxisInqVctSize(zaxisID);
      if ( nlevel == (nvct/2 - 1) )
	{
	  vct = (double*) malloc(nvct*sizeof(double));
	  zaxisInqVct(zaxisID, vct);
	}
      else
	cdoAbort("Unsupported vertical coordinate table format!");
    }
  else
    cdoAbort("Unsupported Z-Axis type!");


  vlistClearFlag(vlistID1);
  for ( levelID = 0; levelID < nlevel; ++levelID )
    vlistDefFlag(vlistID1, omegaID, levelID, TRUE);

  vlistID2 = vlistCreate();
  vlistCopyFlag(vlistID2, vlistID1);
  vlistDefVarCode(vlistID2, 0, 40);
  vlistDefVarName(vlistID2, 0, "W");
  vlistDefVarLongname(vlistID2, 0, "Vertical velocity");
  vlistDefVarUnits(vlistID2, 0, "m/s");
  vlistDefVarMissval(vlistID2, 0, missval_out);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);
     
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);

	  offset = levelID*gridsize;

	  if      ( varID == tempID )
	    streamReadRecord(streamID1, temp+offset, &nmiss);
	  else if ( varID == sqID )
	    streamReadRecord(streamID1, sq+offset, &nmiss);
	  else if ( varID == omegaID )
	    streamReadRecord(streamID1, omega+offset, &nmiss);
	  else if ( varID == psID && zaxisInqType(zaxisID) == ZAXIS_HYBRID )
	    streamReadRecord(streamID1, ps_prog, &nmiss);
	}

      if ( zaxisInqType(zaxisID) == ZAXIS_HYBRID )
	presh(fpress, hpress, vct, ps_prog, nlevel, gridsize);

      for ( levelID = 0; levelID < nlevel; ++levelID )
	{
	  offset = levelID*gridsize;

	  for ( i = 0; i < gridsize; ++i )
	    {
	      if ( DBL_IS_EQUAL(temp[offset+i],missval_t)    || 
		   DBL_IS_EQUAL(omega[offset+i],missval_wap) ||
		   DBL_IS_EQUAL(sq[offset+i],missval_sq) )
		{
		  wms[offset+i] = missval_out;
		}
	      else
		{
	          /* Virtuelle Temperatur bringt die Feuchteabhaengigkeit hinein */
	          tv = temp[offset+i] * (1. + 0.608*sq[offset+i]);

	          /*
		    Die Dichte erhaelt man nun mit der Gasgleichung rho=p/(R*tv)
		    Level in Pa!
	          */
	          rho = fpress[offset+i] / (R*tv);

	          /*
		    Nun daraus die Vertikalgeschwindigkeit im m/s, indem man die
		    Vertikalgeschwindigkeit in Pa/s durch die Erdbeschleunigung
		    und die Dichte teilt
	          */
	          wms[offset+i] = omega[offset+i]/(G*rho);
	        }
            }
	}

      for ( levelID = 0; levelID < nlevel; ++levelID )
	{
	  nmiss_out = 0;
	  for ( i = 0; i < gridsize; i++ )
            if ( DBL_IS_EQUAL(wms[offset+i],missval_out) )
	      nmiss_out++;

	  offset = levelID*gridsize;

	  streamDefRecord(streamID2, 0, levelID);
	  streamWriteRecord(streamID2, wms+offset, nmiss_out);
	}

      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);
 
  vlistDestroy(vlistID2);

  free(temp);
  free(sq);
  free(omega);
  free(wms);
  free(fpress);

  if ( ps_prog ) free(ps_prog);
  if ( hpress )  free(hpress);
  if ( vct ) free(vct);

  free(level);

  cdoFinish();

  return (0);
}
