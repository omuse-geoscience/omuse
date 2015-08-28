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

      Wind       uv2dv           U and V wind to divergence and vorticity
      Wind       dv2uv           Divergence and vorticity to U and V wind
      Wind       dv2ps           Divergence and vorticity to
                                 velocity potential and stream function
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "specspace.h"
#include "list.h"


void *Wind(void *argument)
{
  int UV2DV, UV2DVL, DV2UV, DV2UVL, DV2PS;
  int operatorID;
  int streamID1, streamID2;
  int nrecs, nvars;
  int tsID, recID, varID, levelID;
  int nlev = 0, gridsize;
  int index, ngrids;
  int vlistID1, vlistID2;
  int gridIDsp = -1, gridIDgp = -1;
  int gridID1 = -1, gridID2 = -1;
  int gridID;
  int nmiss;
  int lcopy = FALSE;
  int taxisID1, taxisID2;
  int nlon, nlat, ntr = -1;
  int code, param;
  int pnum, pcat, pdis;
  int varID1 = -1, varID2 = -1;
  int offset;
  SPTRANS *sptrans = NULL;
  DVTRANS *dvtrans = NULL;
  char varname[CDI_MAX_NAME];
  double *array1 = NULL;
  double *ivar1 = NULL, *ivar2 = NULL, *ovar1 = NULL, *ovar2 = NULL;

  cdoInitialize(argument);

  UV2DV  = cdoOperatorAdd("uv2dv",  0, 0, NULL);
  UV2DVL = cdoOperatorAdd("uv2dvl", 0, 0, NULL);
  DV2UV  = cdoOperatorAdd("dv2uv",  0, 0, NULL);
  DV2UVL = cdoOperatorAdd("dv2uvl", 0, 0, NULL);
  DV2PS  = cdoOperatorAdd("dv2ps",  0, 0, NULL);

  operatorID = cdoOperatorID();

  if ( UNCHANGED_RECORD ) lcopy = TRUE;

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  /* find variables */
  nvars = vlistNvars(vlistID2);
  for ( varID = 0; varID < nvars; varID++ )
    {
      param = vlistInqVarParam(vlistID2, varID);
      cdiDecodeParam(param, &pnum, &pcat, &pdis);
      code = pnum;
      if ( operatorID == UV2DV || operatorID == UV2DVL )
	{
	  /* search for u and v wind */
	  if ( code <= 0 )
	    {
	      vlistInqVarName(vlistID1, varID, varname);
	      strtolower(varname);

	      if      ( strcmp(varname, "u") == 0 ) code = 131;
	      else if ( strcmp(varname, "v") == 0 ) code = 132;
	    }

	  if      ( code == 131 ) varID1 = varID;
	  else if ( code == 132 ) varID2 = varID;
	}
      else if ( operatorID == DV2UV || operatorID == DV2UVL || operatorID == DV2PS )
	{
	  /* search for divergence and vorticity */
	  if ( pdis != 255 ) // GRIB2
	    {
	      vlistInqVarName(vlistID1, varID, varname);
	      strtolower(varname);

	      if      ( strcmp(varname, "d")  == 0 ) code = 155;
	      else if ( strcmp(varname, "vo") == 0 ) code = 138;
	    }
	  else if ( code <= 0 )
	    {
	      vlistInqVarName(vlistID1, varID, varname);
	      strtolower(varname);

	      if      ( strcmp(varname, "sd")  == 0 ) code = 155;
	      else if ( strcmp(varname, "svo") == 0 ) code = 138;
	    }

	  if      ( code == 155 ) varID1 = varID;
	  else if ( code == 138 ) varID2 = varID;
	}
      else
	cdoAbort("Unexpected operatorID %d", operatorID);
    }

  ngrids = vlistNgrids(vlistID1);

  /* find first spectral grid */
  for ( index = 0; index < ngrids; index++ )
    {
      gridID = vlistGrid(vlistID1, index);
      if ( gridInqType(gridID) == GRID_SPECTRAL )
	{
	  gridIDsp = gridID;
	  break;
	}
    }

  /* find first gaussian grid */
  for ( index = 0; index < ngrids; index++ )
    {
      gridID = vlistGrid(vlistID1, index);
      if ( gridInqType(gridID) == GRID_GAUSSIAN )
	{
	  gridIDgp = gridID;
	  break;
	}
    }

  /* define output grid */
  if ( operatorID == UV2DV || operatorID == UV2DVL )
    {
      if ( varID1 == -1 ) cdoWarning("U-wind not found!");
      if ( varID2 == -1 ) cdoWarning("V-wind not found!");

      if ( varID1 != -1 && varID2 != -1 )
	{
	  gridID1 = vlistInqVarGrid(vlistID1, varID1);

	  if ( gridInqType(gridID1) != GRID_GAUSSIAN )
	    cdoAbort("U-wind is not on gaussian grid!");

	  if ( gridID1 != vlistInqVarGrid(vlistID1, varID2) )
	    cdoAbort("U and V wind must have the same grid represention!");

	  if ( gridID1 != -1 )
	    {
	      if ( operatorID == UV2DV )
		ntr = nlat2ntr(gridInqYsize(gridID1));
	      else
		ntr = nlat2ntr_linear(gridInqYsize(gridID1));

	      if ( gridIDsp != -1 )
		if ( ntr != gridInqTrunc(gridIDsp) ) gridIDsp = -1;

	      if ( gridIDsp == -1 )
		{
		  gridIDsp = gridCreate(GRID_SPECTRAL, (ntr+1)*(ntr+2));
		  gridDefTrunc(gridIDsp, ntr);
		  gridDefComplexPacking(gridIDsp, 1);
		}
	    }

	  if ( gridIDsp == -1 && gridInqType(vlistGrid(vlistID1, 0)) == GRID_GAUSSIAN_REDUCED )
	    cdoAbort("Gaussian reduced grid found. Use option -R to convert it to a regular grid!");

	  if ( gridIDsp == -1 ) cdoAbort("No Gaussian grid data found!");

	  gridID2 = gridIDsp;

	  vlistChangeVarGrid(vlistID2, varID1, gridID2);
	  vlistChangeVarGrid(vlistID2, varID2, gridID2);
	  vlistDefVarParam(vlistID2, varID1, cdiEncodeParam(155, 128, 255));
	  vlistDefVarParam(vlistID2, varID2, cdiEncodeParam(138, 128, 255));
	  vlistDefVarName(vlistID2, varID1, "sd");
	  vlistDefVarName(vlistID2, varID2, "svo");
	  vlistDefVarLongname(vlistID2, varID1, "divergence");
	  vlistDefVarLongname(vlistID2, varID2, "vorticity");
	  vlistDefVarUnits(vlistID2, varID1, "1/s");
	  vlistDefVarUnits(vlistID2, varID2, "1/s");

	  nlon = gridInqXsize(gridID1);
	  nlat = gridInqYsize(gridID1);
	  ntr  = gridInqTrunc(gridID2);

	  sptrans = sptrans_new(nlon, nlat, ntr, 1);
	}
    }
  else if ( operatorID == DV2UV || operatorID == DV2UVL )
    {   
      if ( varID1 == -1 ) cdoWarning("Divergence not found!");
      if ( varID2 == -1 ) cdoWarning("Vorticity not found!");

      if ( varID1 != -1 && varID2 != -1 )
	{
	  gridID1 = vlistInqVarGrid(vlistID1, varID2);

	  if ( gridInqType(gridID1) != GRID_SPECTRAL )
	    cdoAbort("Vorticity is not on spectral grid!");

	  if ( gridID1 != vlistInqVarGrid(vlistID1, varID1) )
	    cdoAbort("Divergence and vorticity must have the same grid represention!");

	  if ( gridIDgp != -1 )
	    {
	      if ( operatorID == DV2UV )
		ntr = nlat2ntr(gridInqYsize(gridIDgp));
	      else
		ntr = nlat2ntr_linear(gridInqYsize(gridIDgp));
	      
	      if ( gridInqTrunc(gridIDsp) != ntr ) gridIDgp = -1;
	    }

	  if ( gridIDgp == -1 )
	    {
	      char gridname[20];

	      if ( operatorID == DV2UV )
		sprintf(gridname, "t%dgrid", gridInqTrunc(gridIDsp));
	      else
		sprintf(gridname, "tl%dgrid", gridInqTrunc(gridIDsp));
	  
	      gridIDgp = gridFromName(gridname);
	    }

	  gridID2 = gridIDgp;

	  vlistChangeVarGrid(vlistID2, varID1, gridID2);
	  vlistChangeVarGrid(vlistID2, varID2, gridID2);
	  vlistDefVarParam(vlistID2, varID1, cdiEncodeParam(131, 128, 255));
	  vlistDefVarParam(vlistID2, varID2, cdiEncodeParam(132, 128, 255));
	  vlistDefVarName(vlistID2, varID1, "u");
	  vlistDefVarName(vlistID2, varID2, "v");
	  vlistDefVarLongname(vlistID2, varID1, "u-velocity");
	  vlistDefVarLongname(vlistID2, varID2, "v-velocity");
	  vlistDefVarUnits(vlistID2, varID1, "m/s");
	  vlistDefVarUnits(vlistID2, varID2, "m/s");

	  ntr  = gridInqTrunc(gridID1);
	  nlon = gridInqXsize(gridID2);
	  nlat = gridInqYsize(gridID2);
      
	  sptrans = sptrans_new(nlon, nlat, ntr, 0);
	  dvtrans = dvtrans_new(ntr);
	}
    }
  else if ( operatorID == DV2PS )
    {   
      if ( varID1 == -1 ) cdoWarning("Divergence not found!");
      if ( varID2 == -1 ) cdoWarning("Vorticity not found!");

      if ( varID1 != -1 && varID2 != -1 )
	{
	  gridID1 = vlistInqVarGrid(vlistID1, varID2);

	  if ( gridInqType(gridID1) != GRID_SPECTRAL )
	    cdoAbort("Vorticity is not on spectral grid!");

	  if ( gridID1 != vlistInqVarGrid(vlistID1, varID1) )
	    cdoAbort("Divergence and vorticity must have the same grid represention!");

	  vlistDefVarParam(vlistID2, varID1, cdiEncodeParam(149, 128, 255));
	  vlistDefVarParam(vlistID2, varID2, cdiEncodeParam(148, 128, 255));
	  vlistDefVarName(vlistID2, varID1, "velopot");
	  vlistDefVarName(vlistID2, varID2, "stream");
	  vlistDefVarLongname(vlistID2, varID1, "velocity potential");
	  vlistDefVarLongname(vlistID2, varID2, "streamfunction");
	  vlistDefVarUnits(vlistID2, varID1, "m^2/s");
	  vlistDefVarUnits(vlistID2, varID2, "m^2/s");

	  ntr  = gridInqTrunc(gridID1);
	  gridID2 = gridID1;
	}
    }

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  gridsize = vlistGridsizeMax(vlistID1);
  array1 = (double*) malloc(gridsize*sizeof(double));

  if ( varID1 != -1 && varID2 != -1 )
    {
      nlev     = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID1));

      gridsize = gridInqSize(gridID1);
      ivar1 = (double*) malloc(nlev*gridsize*sizeof(double));
      ivar2 = (double*) malloc(nlev*gridsize*sizeof(double));
  
      gridsize = gridInqSize(gridID2);
      ovar1 = (double*) malloc(nlev*gridsize*sizeof(double));
      ovar2 = (double*) malloc(nlev*gridsize*sizeof(double));
    }

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);
	       
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);

	  if ( (varID1 != -1 && varID2 != -1) && (varID == varID1 || varID == varID2) )
	    {
	      streamReadRecord(streamID1, array1, &nmiss);
	      if ( nmiss ) cdoAbort("Missing values unsupported for spectral data!");

	      gridsize = gridInqSize(gridID1);
	      offset = gridsize*levelID;

	      if      ( varID == varID1 )
		memcpy(ivar1+offset, array1, gridsize*sizeof(double));
	      else if ( varID == varID2 )
	        memcpy(ivar2+offset, array1, gridsize*sizeof(double));
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
		  streamReadRecord(streamID1, array1, &nmiss);
		  streamWriteRecord(streamID2, array1, nmiss);
		}
	    }
	}

      if ( varID1 != -1 && varID2 != -1 )
	{
	  if ( operatorID == UV2DV || operatorID == UV2DVL )
	    trans_uv2dv(sptrans, nlev, gridID1, ivar1, ivar2, gridID2, ovar1, ovar2);
	  else if ( operatorID == DV2UV || operatorID == DV2UVL )
	    trans_dv2uv(sptrans, dvtrans, nlev, gridID1, ivar1, ivar2, gridID2, ovar1, ovar2);
	  else if ( operatorID == DV2PS )
	    {
	      dv2ps(ivar1, ovar1, nlev, ntr);
	      dv2ps(ivar2, ovar2, nlev, ntr);
	    }

	  gridsize = gridInqSize(gridID2);
	  if ( operatorID == UV2DV || operatorID == UV2DVL || operatorID == DV2PS )
	    {
	      for ( levelID = 0; levelID < nlev; levelID++ )
		{
		  offset = gridsize*levelID;
		  streamDefRecord(streamID2, varID2, levelID);
		  streamWriteRecord(streamID2, ovar2+offset, 0);
		}
	      for ( levelID = 0; levelID < nlev; levelID++ )
		{
		  offset = gridsize*levelID;
		  streamDefRecord(streamID2, varID1, levelID);
		  streamWriteRecord(streamID2, ovar1+offset, 0);
		}
	    }
	  else if ( operatorID == DV2UV || operatorID == DV2UVL )
	    {
	      for ( levelID = 0; levelID < nlev; levelID++ )
		{
		  offset = gridsize*levelID;
		  streamDefRecord(streamID2, varID1, levelID);
		  streamWriteRecord(streamID2, ovar1+offset, 0);
		}
	      for ( levelID = 0; levelID < nlev; levelID++ )
		{
		  offset = gridsize*levelID;
		  streamDefRecord(streamID2, varID2, levelID);
		  streamWriteRecord(streamID2, ovar2+offset, 0);
		}
	    }
	}

      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  if ( array1 ) free(array1);

  if ( ivar1 ) free(ivar1);
  if ( ivar2 ) free(ivar2);
  if ( ovar1 ) free(ovar1);
  if ( ovar2 ) free(ovar2);

  sptrans_delete(sptrans);
  if ( dvtrans ) dvtrans_delete(dvtrans);

  cdoFinish();

  return (0);
}
