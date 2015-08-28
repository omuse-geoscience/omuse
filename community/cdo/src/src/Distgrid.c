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

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"

#define  MAX_BLOCKS  16384

static
void genGrids(int gridID1, int *gridIDs, int nxvals, int nyvals, int nxblocks, int nyblocks,
	      int **gridindex, int *ogridsize, int nsplit)
{
  int gridID2;
  int gridtype;
  int nx, ny;
  int gridsize2;
  int index, i, j, ix, iy, offset;
  int *xlsize = NULL, *ylsize = NULL;
  double *xvals = NULL, *yvals = NULL;

  gridtype = gridInqType(gridID1);
  if ( !(gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN || gridtype == GRID_GENERIC) )
    cdoAbort("Unsupported grid type: %s!", gridNamePtr(gridtype));

  nx = gridInqXsize(gridID1);
  ny = gridInqYsize(gridID1);

  xvals = (double*) malloc(nx*sizeof(double));
  yvals = (double*) malloc(ny*sizeof(double));

  xlsize = (int*) malloc(nxblocks*sizeof(int));
  ylsize = (int*) malloc(nyblocks*sizeof(int));

  gridInqXvals(gridID1, xvals);
  gridInqYvals(gridID1, yvals);

  for ( ix = 0; ix < nxblocks; ++ix ) xlsize[ix] = nxvals;
  if ( nx%nxblocks != 0 ) xlsize[nxblocks-1] = nx - (nxblocks-1)*nxvals;
  if ( cdoVerbose ) for ( ix = 0; ix < nxblocks; ++ix ) cdoPrint("xblock %d: %d", ix, xlsize[ix]);

  for ( iy = 0; iy < nyblocks; ++iy ) ylsize[iy] = nyvals;
  if ( ny%nyblocks != 0 ) ylsize[nyblocks-1] = ny - (nyblocks-1)*nyvals;
  if ( cdoVerbose ) for ( iy = 0; iy < nyblocks; ++iy ) cdoPrint("yblock %d: %d", iy, ylsize[iy]);

  index = 0;
  for ( iy = 0; iy < nyblocks; ++iy )
    for ( ix = 0; ix < nxblocks; ++ix )
      {
	offset = iy*nyvals*nx + ix*nxvals;

	gridsize2 = xlsize[ix]*ylsize[iy];
	gridindex[index] = (int*) malloc(gridsize2*sizeof(int));

	gridsize2 = 0;
        // printf("iy %d, ix %d offset %d\n", iy, ix,  offset);
	for ( j = 0; j < ylsize[iy]; ++j )
	  {
	    for ( i = 0; i < xlsize[ix]; ++i )
	      {
	       	// printf(">> %d %d %d\n", j, i, offset + j*nx + i);
		gridindex[index][gridsize2++] = offset + j*nx + i;
	      }
	  }
	// printf("gridsize2 %d\n", gridsize2);

	gridID2 = gridCreate(gridtype, gridsize2);
	gridDefXsize(gridID2, xlsize[ix]);
	gridDefYsize(gridID2, ylsize[iy]);
	gridDefXvals(gridID2, xvals+ix*nxvals);
	gridDefYvals(gridID2, yvals+iy*nyvals);

	gridIDs[index] = gridID2;
	ogridsize[index] = gridsize2;

	index++;
	if ( index > nsplit )
	  cdoAbort("Internal problem, index exceeded bounds!");
      }

  free(xvals);
  free(yvals);
  free(xlsize);
  free(ylsize);
}

static
void window_cell(double *array1, double *array2, long gridsize2, int *cellidx)
{
  long i;

  for ( i = 0; i < gridsize2; ++i )
    array2[i] = array1[cellidx[i]];
}

typedef struct
{
  int gridID;
  int *gridIDs;
  int *gridsize;
  int **gridindex;
} sgrid_t;


void *Distgrid(void *argument)
{
  int nchars;
  int streamID1;
  int *vlistIDs = NULL, *streamIDs = NULL;
  int gridID1 = -1, varID;
  int nrecs, ngrids;
  int tsID, recID, levelID;
  int vlistID1;
  char filesuffix[32];
  char filename[8192];
  const char *refname;
  int index;
  int nsplit;
  int xinc = 1, yinc = 1;
  int gridsize, gridsize2max;
  int gridtype = -1;
  int nmiss;
  int nxblocks = 1, nyblocks = 1;
  int nx, ny, i;
  double missval;
  double *array1 = NULL, *array2 = NULL;
  sgrid_t *grids;

  cdoInitialize(argument);

  operatorInputArg("nxblocks, [nyblocks]");
  if ( operatorArgc() < 1 ) cdoAbort("Too few arguments!");
  if ( operatorArgc() > 2 ) cdoAbort("Too many arguments!");
  nxblocks = parameter2int(operatorArgv()[0]);
  if ( operatorArgc() == 2 ) nyblocks = parameter2int(operatorArgv()[1]);

  if ( nxblocks <= 0 ) cdoAbort("nxblocks has to be greater than 0!");
  if ( nyblocks <= 0 ) cdoAbort("nyblocks has to be greater than 0!");

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);

  ngrids = vlistNgrids(vlistID1);

  for ( index = 0; index < ngrids; index++ )
    {
      gridID1 = vlistGrid(vlistID1, index);
      gridtype = gridInqType(gridID1);
      if ( gridtype == GRID_LONLAT   || gridtype == GRID_GAUSSIAN ||
	  (gridtype == GRID_GENERIC && gridInqXsize(gridID1) > 0 && gridInqYsize(gridID1) > 0) )
	   break;
    }

  if ( index == ngrids )
    cdoAbort("No Lon/Lat, Gaussian or Generic grid found (%s data unsupported)!", gridNamePtr(gridtype));

  gridID1 = vlistGrid(vlistID1, 0);
  gridsize = gridInqSize(gridID1);
  nx = gridInqXsize(gridID1);
  ny = gridInqYsize(gridID1);
  for ( i = 1; i < ngrids; i++ )
    {
      gridID1 = vlistGrid(vlistID1, i);
      if ( gridsize != gridInqSize(gridID1) )
	cdoAbort("Gridsize must not change!");
    }

  if ( nxblocks > nx )
    {
      cdoPrint("nxblocks (%d) greater than nx (%d), set to %d!", nxblocks, nx, nx);
      nxblocks = nx;
    }
  if ( nyblocks > ny )
    {
      cdoPrint("nyblocks (%d) greater than ny (%d), set to %d!", nyblocks, ny, ny);
      nyblocks = ny;
    }

  xinc = nx/nxblocks;
  yinc = ny/nyblocks;

  if ( nx%xinc != 0 ) xinc++;
  if ( ny%yinc != 0 ) yinc++;

  nsplit = nxblocks*nyblocks;
  if ( nsplit > MAX_BLOCKS ) cdoAbort("Too many blocks (max = %d)!", MAX_BLOCKS);

  array1 = (double*) malloc(gridsize*sizeof(double));

  vlistIDs  = (int*) malloc(nsplit*sizeof(int));
  streamIDs = (int*) malloc(nsplit*sizeof(int));

  grids = (sgrid_t*) malloc(ngrids*sizeof(sgrid_t));
  for ( i = 0; i < ngrids; i++ )
    {  
      gridID1 = vlistGrid(vlistID1, i);
      grids[i].gridID    = vlistGrid(vlistID1, i);
      grids[i].gridIDs   = (int*) malloc(nsplit*sizeof(int));
      grids[i].gridsize  = (int*) malloc(nsplit*sizeof(int));
      grids[i].gridindex = (int**) malloc(nsplit*sizeof(int*));

      for ( index = 0; index < nsplit; index++ ) grids[i].gridindex[index] = NULL;
    }

  for ( index = 0; index < nsplit; index++ )
    vlistIDs[index] = vlistDuplicate(vlistID1);

  if ( cdoVerbose ) cdoPrint("ngrids=%d  nsplit=%d", ngrids, nsplit);

  for ( i = 0; i < ngrids; i++ )
    {
      gridID1 = vlistGrid(vlistID1, i);
      genGrids(gridID1, grids[i].gridIDs, xinc, yinc, nxblocks, nyblocks, grids[i].gridindex, grids[i].gridsize, nsplit);
      /*
      if ( cdoVerbose )
	for ( index = 0; index < nsplit; index++ )
	  cdoPrint("Block %d,  gridID %d,  gridsize %d", index+1, grids[i].gridIDs[index], gridInqSize(grids[i].gridIDs[index]));
      */
      for ( index = 0; index < nsplit; index++ )
	vlistChangeGridIndex(vlistIDs[index], i, grids[i].gridIDs[index]);
    }

  gridsize2max = 0;
  for ( index = 0; index < nsplit; index++ )
    if ( grids[0].gridsize[index] > gridsize2max ) gridsize2max = grids[0].gridsize[index];

  array2 = (double*) malloc(gridsize2max*sizeof(double));

  strcpy(filename, cdoStreamName(1)->args);
  nchars = strlen(filename);

  refname = cdoStreamName(0)->argv[cdoStreamName(0)->argc-1];
  filesuffix[0] = 0;
  cdoGenFileSuffix(filesuffix, sizeof(filesuffix), streamInqFiletype(streamID1), vlistID1, refname);

  for ( index = 0; index < nsplit; index++ )
    {
      sprintf(filename+nchars, "%05d", index);
      if ( filesuffix[0] )
	sprintf(filename+nchars+5, "%s", filesuffix);

      argument_t *fileargument = file_argument_new(filename);
      streamIDs[index] = streamOpenWrite(fileargument, cdoFiletype());
      file_argument_free(fileargument);

      streamDefVlist(streamIDs[index], vlistIDs[index]);
    }

  if ( ngrids > 1 ) cdoPrint("Bausstelle: number of different grids > 1!");
  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      for ( index = 0; index < nsplit; index++ )
	streamDefTimestep(streamIDs[index], tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamReadRecord(streamID1, array1, &nmiss);

	  missval = vlistInqVarMissval(vlistID1, varID);

	  for ( index = 0; index < nsplit; index++ )
	    {
	      i = 0;
	      window_cell(array1, array2, grids[i].gridsize[index], grids[i].gridindex[index]);
	      streamDefRecord(streamIDs[index], varID, levelID);
	      if ( nmiss > 0 )
		{
		  nmiss = 0;
		  for ( int k = 0; k < grids[i].gridsize[index]; ++k )
		    if ( DBL_IS_EQUAL(array2[k], missval) ) nmiss++;
		}
	      streamWriteRecord(streamIDs[index], array2, nmiss);
	    }
	}

      tsID++;
    }

  streamClose(streamID1);

  for ( index = 0; index < nsplit; index++ )
    {
      streamClose(streamIDs[index]);
      vlistDestroy(vlistIDs[index]);
    }

  if ( array1 ) free(array1);
  if ( array2 ) free(array2);

  if ( vlistIDs  ) free(vlistIDs);
  if ( streamIDs ) free(streamIDs);

  for ( i = 0; i < ngrids; i++ )
    {
      for ( index = 0; index < nsplit; index++ )
	gridDestroy(grids[i].gridIDs[index]);
      free(grids[i].gridIDs);
      free(grids[i].gridsize);

      for ( index = 0; index < nsplit; index++ )
	free(grids[i].gridindex[index]);
      free(grids[i].gridindex);
    }
  free(grids);

  cdoFinish();

  return (0);
}
