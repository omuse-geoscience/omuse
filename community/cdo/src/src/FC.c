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

      FC         fc2sp           Fourier to spectral
      FC         sp2fc           Spectral to fourier
      FC         fc2gp           Fourier to gridpoint
      FC         gp2fc           Gridpoint to fourier
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "specspace.h"
#include "list.h"


void *FC(void *argument)
{
  int FC2SP, SP2FC, FC2GP, GP2FC;
  int operatorID;
  int streamID1, streamID2;
  int nrecs, nvars;
  int tsID, recID, varID, levelID;
  int gridsize;
  int index, ngrids;
  int vlistID1, vlistID2;
  int gridIDsp = -1, gridIDgp = -1, gridIDfc = -1;
  int gridID1 = -1, gridID2 = -1;
  int gridID;
  int nmiss;
  int *vars;
  int lcopy = FALSE;
  int taxisID1, taxisID2;
  int nlon = 0, nlat = 0, ntr = 0;
  int nsp = 0, nfc = 0;
  double *array1 = NULL, *array2 = NULL;
  SPTRANS *sptrans = NULL;

  cdoInitialize(argument);

  FC2SP  = cdoOperatorAdd("fc2sp",  0, 0, NULL);
  SP2FC  = cdoOperatorAdd("sp2fc", 0, 0, NULL);
  FC2GP  = cdoOperatorAdd("fc2gp",  0, 0, NULL);
  GP2FC  = cdoOperatorAdd("gp2fc", 0, 0, NULL);

  operatorID = cdoOperatorID();

  if ( UNCHANGED_RECORD ) lcopy = TRUE;

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

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
  /* find first fourier grid */
  for ( index = 0; index < ngrids; index++ )
    {
      gridID = vlistGrid(vlistID1, index);
      if ( gridInqType(gridID) == GRID_FOURIER )
	{
	  gridIDfc = gridID;
	  break;
	}
    }

  /* define output grid */
  if ( operatorID == FC2SP )
    {
      if ( gridIDfc == -1 ) cdoWarning("No fourier data found!");

      gridID1 = gridIDfc;

      if ( gridID1 != -1 )
	{
	  nfc  = gridInqSize(gridID1);
	  ntr  = gridInqTrunc(gridID1);
	  nlat = nfc2nlat(nfc, ntr);

	  if ( gridIDsp != -1 )
	    if ( ntr != gridInqTrunc(gridIDsp) ) gridIDsp = -1;

	  if ( gridIDsp == -1 )
	    {
	      nsp = (ntr+1)*(ntr+2);
	      gridIDsp = gridCreate(GRID_SPECTRAL, nsp);
	      gridDefTrunc(gridIDsp, ntr);
	      gridDefComplexPacking(gridIDsp, 1);
	    }

	  gridID2 = gridIDsp;
	  nlon = 2*nlat;
	  ntr  = gridInqTrunc(gridID2);

	  sptrans = sptrans_new(nlon, nlat, ntr, 0);
	}
    }
  else if ( operatorID == SP2FC )
    {   
      if ( gridIDsp == -1 ) cdoWarning("No spectral data found!");

      gridID1 = gridIDsp;

      if ( gridID1 != -1 )
	{
	  ntr  = gridInqTrunc(gridID1);
	  nlat = ntr2nlat(ntr);

	  if ( gridIDfc != -1 )
	    {
	      if ( ntr != gridInqTrunc(gridIDfc) ) gridIDfc = -1;
	    }

	  if ( gridIDfc == -1 )
	    {
	      nfc = 2*nlat*(ntr+1);
              gridIDfc = gridCreate(GRID_FOURIER, nfc);
	      gridDefTrunc(gridIDfc, ntr);
	    }

	  gridID2 = gridIDfc;
	  nlon = 2*nlat;
      
	  sptrans = sptrans_new(nlon, nlat, ntr, 0);
	}
    }
  else if ( operatorID == GP2FC )
    {
      if ( gridIDgp == -1 ) cdoWarning("No Gaussian grid data found!");

      gridID1 = gridIDgp;

      if ( gridID1 != -1 )
	{
	  nlon = gridInqXsize(gridID1);
	  nlat = gridInqYsize(gridID1);
	  ntr  = nlat2ntr(nlat);

	  if ( gridIDfc != -1 )
	    if ( ntr != gridInqTrunc(gridIDfc) ) gridIDfc = -1;

	  if ( gridIDfc == -1 )
	    {
	      nfc = 2*nlat*(ntr+1);
	      gridIDfc = gridCreate(GRID_FOURIER, nfc);
	      gridDefTrunc(gridIDfc, ntr);
	    }

	  gridID2 = gridIDfc;
	  nfc  = gridInqSize(gridID2);

	  sptrans = sptrans_new(nlon, nlat, ntr, 0);
 	}
    }
  else if ( operatorID == FC2GP )
    {   
      if ( gridIDfc == -1 ) cdoWarning("No fourier data found!");

      gridID1 = gridIDfc;

      if ( gridID1 != -1 )
	{
	  nfc  = gridInqSize(gridID1);
	  ntr  = gridInqTrunc(gridID1);
	  nlat = nfc2nlat(nfc, ntr);

	  if ( gridIDgp != -1 )
	    {
	      if ( nlat != gridInqYsize(gridIDgp) ) gridIDgp = -1;
	    }

	  if ( gridIDgp == -1 )
	    {
	      char gridname[20];
	      sprintf(gridname, "t%dgrid", ntr);

	      gridIDgp = gridFromName(gridname);
	    }

	  gridID2 = gridIDgp;
	  nlon = gridInqXsize(gridID2);
	  nlat = gridInqYsize(gridID2);
      
	  sptrans = sptrans_new(nlon, nlat, ntr, 0);
	}
    }

  // printf("nfc %d, ntr %d, nlat %d, nlon %d\n", nfc, ntr, nlat, nlon);

  nvars = vlistNvars(vlistID2);
  vars  = (int*) malloc(nvars*sizeof(int));
  for ( varID = 0; varID < nvars; varID++ )
    {
      if ( gridID1 == vlistInqVarGrid(vlistID1, varID) )
	vars[varID] = TRUE;
      else
	vars[varID] = FALSE;
    }

  if ( gridID1 != -1 ) vlistChangeGrid(vlistID2, gridID1, gridID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  gridsize = vlistGridsizeMax(vlistID1);
  array1 = (double*) malloc(gridsize*sizeof(double));

  if ( gridID2 != -1 )
    {
      gridsize = gridInqSize(gridID2);
      array2 = (double*) malloc(gridsize*sizeof(double));
    }

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);
	       
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);

	  if ( vars[varID] )
	    {
	      streamReadRecord(streamID1, array1, &nmiss);
	      if ( nmiss ) cdoAbort("Missing values unsupported for spectral/fourier data!");

	      gridID1 = vlistInqVarGrid(vlistID1, varID);
	      if ( operatorID == FC2SP )
		four2spec(sptrans, gridID1, array1, gridID2, array2);	      
	      else if ( operatorID == SP2FC )
		spec2four(sptrans, gridID1, array1, gridID2, array2);
	      else if ( operatorID == FC2GP )
		four2grid(sptrans, gridID1, array1, gridID2, array2);	      
	      else if ( operatorID == GP2FC )
		grid2four(sptrans, gridID1, array1, gridID2, array2);

	      streamDefRecord(streamID2, varID, levelID);
	      streamWriteRecord(streamID2, array2, nmiss);  
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
      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  if ( array2 ) free(array2);
  if ( array1 ) free(array1);
  if ( vars )   free(vars);

  sptrans_delete(sptrans);

  cdoFinish();

  return (0);
}
