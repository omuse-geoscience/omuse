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

      Setgrid    setgrid         Set grid
      Setgrid    setgridtype     Set grid type
      Setgrid    setgridarea     Set grid area
      Setgrid    setgridmask     Set grid mask
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "grid.h"


void *Setgrid(void *argument)
{
  int SETGRID, SETGRIDTYPE, SETGRIDAREA, SETGRIDMASK, UNSETGRIDMASK, SETGRIDNUMBER, SETGRIDURI;
  int operatorID;
  int streamID1, streamID2 = CDI_UNDEFID;
  int nrecs;
  int tsID, recID, varID, levelID;
  int vlistID1, vlistID2;
  int taxisID1, taxisID2;
  int gridID1, gridID2 = -1;
  int ngrids, index;
  int gridtype = -1;
  int nmiss;
  int found;
  long i, gridsize;
  long areasize = 0;
  long masksize = 0;
  int lregular = 0;
  int ldereference = 0;
  int ligme = 0;
  int number = 0, position = 0;
  int grid2_nvgp;
  int *grid2_vgpm = NULL;
  char *gridname = NULL;
  char *griduri = NULL;
  double *gridmask = NULL;
  double *areaweight = NULL;
  double *array = NULL;

  cdoInitialize(argument);

  SETGRID       = cdoOperatorAdd("setgrid",       0, 0, "grid description file or name");
  SETGRIDTYPE   = cdoOperatorAdd("setgridtype",   0, 0, "grid type");
  SETGRIDAREA   = cdoOperatorAdd("setgridarea",   0, 0, "filename with area weights");
  SETGRIDMASK   = cdoOperatorAdd("setgridmask",   0, 0, "filename with grid mask");
  UNSETGRIDMASK = cdoOperatorAdd("unsetgridmask", 0, 0, NULL);
  SETGRIDNUMBER = cdoOperatorAdd("setgridnumber", 0, 0, "grid number and optionally grid position");
  SETGRIDURI    = cdoOperatorAdd("setgriduri",    0, 0, "reference URI of the horizontal grid");

  operatorID = cdoOperatorID();

  if ( operatorID != UNSETGRIDMASK )
    operatorInputArg(cdoOperatorEnter(operatorID));  

  if ( operatorID == SETGRID )
    {
      operatorCheckArgc(1);
      gridID2 = cdoDefineGrid(operatorArgv()[0]);
    }
  else if ( operatorID == SETGRIDTYPE )
    {
      operatorCheckArgc(1);
      gridname = operatorArgv()[0];

      if      ( strcmp(gridname, "curvilinear") == 0 )   gridtype = GRID_CURVILINEAR;
      else if ( strcmp(gridname, "cell") == 0 )          gridtype = GRID_UNSTRUCTURED;
      else if ( strcmp(gridname, "unstructured") == 0 )  gridtype = GRID_UNSTRUCTURED;
      else if ( strcmp(gridname, "dereference") == 0 )   ldereference = 1;
      else if ( strcmp(gridname, "lonlat") == 0 )        gridtype = GRID_LONLAT;
      else if ( strcmp(gridname, "gaussian") == 0 )      gridtype = GRID_GAUSSIAN;
      else if ( strcmp(gridname, "regular") == 0 )      {gridtype = GRID_GAUSSIAN; lregular = 1;}
      else cdoAbort("Unsupported grid name: %s", gridname);
    }
  else if ( operatorID == SETGRIDAREA )
    {
      int streamID, vlistID, gridID;
      char *areafile;

      operatorCheckArgc(1);
      areafile = operatorArgv()[0];

      argument_t *fileargument = file_argument_new(areafile);
      streamID = streamOpenRead(fileargument);
      file_argument_free(fileargument);

      vlistID = streamInqVlist(streamID);

      nrecs = streamInqTimestep(streamID, 0);
      streamInqRecord(streamID, &varID, &levelID);

      gridID = vlistInqVarGrid(vlistID, varID);
      areasize = gridInqSize(gridID);
      areaweight = (double*) malloc(areasize*sizeof(double));
  
      streamReadRecord(streamID, areaweight, &nmiss);

      streamClose(streamID);

      if ( cdoVerbose )
	{
	  double arrmean, arrmin, arrmax;

	  arrmean = areaweight[0];
	  arrmin  = areaweight[0];
	  arrmax  = areaweight[0];
	  for ( i = 1; i < areasize; i++ )
	    {
	      if ( areaweight[i] < arrmin ) arrmin = areaweight[i];
	      if ( areaweight[i] > arrmax ) arrmax = areaweight[i];
	      arrmean += areaweight[i];
	    }
	  arrmean = arrmean/areasize;

	  cdoPrint("areaweights: %d %#12.5g%#12.5g%#12.5g", areasize, arrmin, arrmean, arrmax);
	}
    }
  else if ( operatorID == SETGRIDMASK )
    {
      int streamID, vlistID, gridID;
      char *maskfile;
      double missval;

      operatorCheckArgc(1);
      maskfile = operatorArgv()[0];
      argument_t *fileargument = file_argument_new(maskfile);
      streamID = streamOpenRead(fileargument);
      file_argument_free(fileargument);

      vlistID = streamInqVlist(streamID);

      nrecs = streamInqTimestep(streamID, 0);
      streamInqRecord(streamID, &varID, &levelID);

      missval  = vlistInqVarMissval(vlistID, varID);
      gridID   = vlistInqVarGrid(vlistID, varID);
      masksize = gridInqSize(gridID);
      gridmask = (double*) malloc(masksize*sizeof(double));
  
      streamReadRecord(streamID, gridmask, &nmiss);

      streamClose(streamID);

      for ( i = 0; i < masksize; i++ )
	if ( DBL_IS_EQUAL(gridmask[i], missval) ) gridmask[i] = 0;
    }
  else if ( operatorID == SETGRIDNUMBER )
    {
      if ( operatorArgc() >= 1 && operatorArgc() <= 2 )
	{
	  number = parameter2int(operatorArgv()[0]);
	  if ( operatorArgc() == 2 ) position = parameter2int(operatorArgv()[1]);
	}
      else
	{
	  operatorCheckArgc(1);
	}
    }
  else if ( operatorID == SETGRIDURI )
    {
      operatorCheckArgc(1);
      griduri = operatorArgv()[0];
    }

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  if ( operatorID == SETGRID )
    {
      found = 0;
      ngrids = vlistNgrids(vlistID1);
      for ( index = 0; index < ngrids; index++ )
	{
	  gridID1 = vlistGrid(vlistID1, index);

	  if ( gridInqSize(gridID1) == gridInqSize(gridID2) )
	    {
	      vlistChangeGridIndex(vlistID2, index, gridID2);
	      found++;
	    }
	}
      if ( ! found ) cdoWarning("No grid with %d points found!", gridInqSize(gridID2));
    }
  else if ( operatorID == SETGRIDNUMBER || operatorID == SETGRIDURI )
    {
      gridID1 = vlistGrid(vlistID1, 0);

      if ( operatorID == SETGRIDNUMBER )
	{
	  gridID2 = gridCreate(GRID_UNSTRUCTURED, gridInqSize(gridID1));
	  gridDefNumber(gridID2, number);
	  gridDefPosition(gridID2, position);
	}
      else
	{
	  gridID2 = gridDuplicate(gridID1);
	  gridDefReference(gridID2, griduri);
	}

      found = 0;
      ngrids = vlistNgrids(vlistID1);
      for ( index = 0; index < ngrids; index++ )
	{
	  gridID1 = vlistGrid(vlistID1, index);

	  if ( gridInqSize(gridID1) == gridInqSize(gridID2) )
	    {
	      vlistChangeGridIndex(vlistID2, index, gridID2);
	      found++;
	    }
	}
      if ( ! found ) cdoWarning("No horizontal grid with %d cells found!", gridInqSize(gridID2));
    }
  else if ( operatorID == SETGRIDTYPE )
    {
      ngrids = vlistNgrids(vlistID1);
      for ( index = 0; index < ngrids; index++ )
	{
	  gridID1 = vlistGrid(vlistID1, index);
	  gridID2 = -1;

	  if ( gridInqType(gridID1) == GRID_GENERIC && gridInqSize(gridID1) == 1 ) continue;
	  
	  if ( lregular )
	    {
	      if ( gridInqType(gridID1) == GRID_GAUSSIAN_REDUCED )
		{
		  gridID2 = gridToRegular(gridID1);
		}
	    }
	  else if ( ldereference )
	    {
	      gridID2 = referenceToGrid(gridID1);
	      if ( gridID2 == -1 ) cdoAbort("Reference to horizontal grid not found!");
	    }
	  else
	    {
	      if      ( gridtype == GRID_CURVILINEAR  )
		{
		  gridID2 = gridToCurvilinear(gridID1, 1);
		}
	      else if ( gridtype == GRID_UNSTRUCTURED )
		{
		  if ( gridInqType(gridID1) == GRID_GME ) ligme = 1;
		  gridID2 = gridToUnstructured(gridID1, 1);

		  if ( ligme )
		    {
		      grid2_nvgp = gridInqSize(gridID2);
		      grid2_vgpm = (int*) malloc(grid2_nvgp*sizeof(int));
		      gridInqMaskGME(gridID2, grid2_vgpm);
		      gridCompress(gridID2);
		    }
		}
	      else if ( gridtype == GRID_LONLAT && gridInqType(gridID1) == GRID_CURVILINEAR )
		{
		  gridID2 = gridCurvilinearToRegular(gridID1);
		  if ( gridID2 == -1 ) cdoWarning("Conversion of curvilinear grid to regular grid failed!");
 		}
	      else if ( gridtype == GRID_LONLAT && gridInqType(gridID1) == GRID_UNSTRUCTURED )
		{
		  gridID2 = -1;
		  if ( gridID2 == -1 ) cdoWarning("Conversion of unstructured grid to regular grid failed!");
 		}
	      else if ( gridtype == GRID_LONLAT && gridInqType(gridID1) == GRID_GENERIC )
		{
		  gridID2 = -1;
		  if ( gridID2 == -1 ) cdoWarning("Conversion of generic grid to regular grid failed!");
 		}
	      else if ( gridtype == GRID_LONLAT && gridInqType(gridID1) == GRID_LONLAT )
		{
		  gridID2 = gridID1;
		}
	      else cdoAbort("Unsupported grid name: %s", gridname);
	    }

	  if ( gridID2 == -1 ) cdoAbort("Unsupported grid type!");

	  vlistChangeGridIndex(vlistID2, index, gridID2);
	}
    }
  else if ( operatorID == SETGRIDAREA )
    {
      ngrids = vlistNgrids(vlistID1);
      for ( index = 0; index < ngrids; index++ )
	{
	  gridID1  = vlistGrid(vlistID1, index);
	  gridsize = gridInqSize(gridID1);
	  if ( gridsize == areasize )
	    {
	      gridID2 = gridDuplicate(gridID1);
	      gridDefArea(gridID2, areaweight);
	      vlistChangeGridIndex(vlistID2, index, gridID2);
	    }
	}
    }
  else if ( operatorID == SETGRIDMASK )
    {
      ngrids = vlistNgrids(vlistID1);
      for ( index = 0; index < ngrids; index++ )
	{
	  gridID1  = vlistGrid(vlistID1, index);
	  gridsize = gridInqSize(gridID1);
	  if ( gridsize == masksize )
	    {
	      int *mask = (int*) malloc(masksize*sizeof(int));
	      for ( i = 0; i < masksize; i++ )
		{
		  if ( gridmask[i] < 0 || gridmask[i] > 255 )
		    mask[i] = 0;
		  else
		    mask[i] = (int)lround(gridmask[i]);
		}
	      gridID2 = gridDuplicate(gridID1);
	      gridDefMask(gridID2, mask);
	      vlistChangeGridIndex(vlistID2, index, gridID2);
	      free(mask);
	    }
	}
    }
  else if ( operatorID == UNSETGRIDMASK )
    {
      ngrids = vlistNgrids(vlistID1);
      for ( index = 0; index < ngrids; index++ )
	{
	  gridID1  = vlistGrid(vlistID1, index);
	  gridID2 = gridDuplicate(gridID1);
	  gridDefMask(gridID2, NULL);
	  vlistChangeGridIndex(vlistID2, index, gridID2);
	}
    }

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);
  //vlistPrint(vlistID2);

  if ( lregular )
    gridsize = vlistGridsizeMax(vlistID2);
  else
    gridsize = vlistGridsizeMax(vlistID1);

  if ( vlistNumber(vlistID1) != CDI_REAL ) gridsize *= 2;
  array = (double*) malloc(gridsize*sizeof(double));

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);
	       
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamDefRecord(streamID2,  varID,  levelID);
	  
	  streamReadRecord(streamID1, array, &nmiss);

	  gridID1 = vlistInqVarGrid(vlistID1, varID);
	  if ( lregular )
	    {
	      gridID2 = vlistInqVarGrid(vlistID2, varID);
	      if ( gridInqType(gridID1) == GRID_GAUSSIAN_REDUCED )
		{
		  double missval = vlistInqVarMissval(vlistID1, varID);
		  field2regular(gridID1, gridID2, missval, array, nmiss);
		}
	    }
	  else if ( gridInqType(gridID1) == GRID_GME )
	    {
	      int j = 0;
	      gridsize = gridInqSize(gridID1);
	      for ( i = 0; i < gridsize; i++ )
		if ( grid2_vgpm[i] ) array[j++] = array[i];
	    }

	  streamWriteRecord(streamID2, array, nmiss);
	}
      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  if ( gridmask ) free(gridmask);
  if ( areaweight ) free(areaweight);
  if ( array ) free(array);
  if ( grid2_vgpm ) free(grid2_vgpm);

  cdoFinish();

  return (0);
}
