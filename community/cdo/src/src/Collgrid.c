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
#include "util.h"


typedef struct
{
  int streamID;
  int vlistID;
  int gridID;
  double *array;
} ens_file_t;


typedef struct
{
  double x, y;
  int id;
}
xyinfo_t;

static
int cmpx(const void *s1, const void *s2)
{
  int cmp = 0;
  const xyinfo_t *xy1 = (const xyinfo_t *) s1;
  const xyinfo_t *xy2 = (const xyinfo_t *) s2;

  if      ( xy1->x < xy2->x ) cmp = -1;
  else if ( xy1->x > xy2->x ) cmp =  1;

  return (cmp);
}

static
int cmpxy_lt(const void *s1, const void *s2)
{
  int cmp = 0;
  const xyinfo_t *xy1 = (const xyinfo_t *) s1;
  const xyinfo_t *xy2 = (const xyinfo_t *) s2;

  if      ( xy1->y < xy2->y || (!(fabs(xy1->y - xy2->y) > 0) && xy1->x < xy2->x) ) cmp = -1;
  else if ( xy1->y > xy2->y || (!(fabs(xy1->y - xy2->y) > 0) && xy1->x > xy2->x) ) cmp =  1;

  return (cmp);
}

static
int cmpxy_gt(const void *s1, const void *s2)
{
  int cmp = 0;
  const xyinfo_t *xy1 = (const xyinfo_t *) s1;
  const xyinfo_t *xy2 = (const xyinfo_t *) s2;

  if      ( xy1->y > xy2->y || (!(fabs(xy1->y - xy2->y) > 0) && xy1->x < xy2->x) ) cmp = -1;
  else if ( xy1->y < xy2->y || (!(fabs(xy1->y - xy2->y) > 0) && xy1->x > xy2->x) ) cmp =  1;

  return (cmp);
}

static
int genGrid(int nfiles, ens_file_t *ef, int **gridindex, int igrid)
{
  int lsouthnorth = TRUE;
  int fileID;
  int gridID;
  int gridID2 = -1;
  int gridtype = -1;
  int *xsize, *ysize;
  int *xoff, *yoff;
  int xsize2, ysize2;
  int idx;
  int nx, ny, ix, iy, i, j, ij, offset;
  double **xvals, **yvals;
  double *xvals2, *yvals2;
  xyinfo_t *xyinfo;

  gridID   = vlistGrid(ef[0].vlistID, igrid);
  gridtype = gridInqType(gridID);
  if ( gridtype == GRID_GENERIC && gridInqXsize(gridID) == 0 && gridInqYsize(gridID) == 0 )
    return (gridID2);

  xsize = (int*) malloc(nfiles*sizeof(int));
  ysize = (int*) malloc(nfiles*sizeof(int));
  xyinfo = (xyinfo_t*) malloc(nfiles*sizeof(xyinfo_t));
  xvals = (double**) malloc(nfiles*sizeof(double*));
  yvals = (double**) malloc(nfiles*sizeof(double*));

  for ( fileID = 0; fileID < nfiles; fileID++ )
    {
      gridID   = vlistGrid(ef[fileID].vlistID, igrid);
      gridtype = gridInqType(gridID);
      if ( !(gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN ||
	    (gridtype == GRID_GENERIC && gridInqXsize(gridID) > 0 && gridInqYsize(gridID) > 0)) )
	cdoAbort("Unsupported grid type: %s!", gridNamePtr(gridtype));

      xsize[fileID] = gridInqXsize(gridID);
      ysize[fileID] = gridInqYsize(gridID);
      /*
      if ( xsize == 0 ) xsize = gridInqXsize(gridID);
      if ( ysize == 0 ) ysize = gridInqYsize(gridID);
      if ( xsize != gridInqXsize(gridID) ) cdoAbort("xsize differ!");
      if ( ysize != gridInqYsize(gridID) ) cdoAbort("ysize differ!");
      */
      xvals[fileID] = (double*) malloc(xsize[fileID]*sizeof(double));
      yvals[fileID] = (double*) malloc(ysize[fileID]*sizeof(double));
      gridInqXvals(gridID, xvals[fileID]);
      gridInqYvals(gridID, yvals[fileID]);

      // printf("fileID %d, gridID %d\n", fileID, gridID);

      xyinfo[fileID].x  = xvals[fileID][0];
      xyinfo[fileID].y  = yvals[fileID][0];
      xyinfo[fileID].id = fileID;

      if ( ysize[fileID] > 1 )
	{
	  if ( yvals[fileID][0] > yvals[fileID][ysize[fileID]-1] ) lsouthnorth = FALSE;
	}
    }
  /*
  if ( cdoVerbose )
    for ( fileID = 0; fileID < nfiles; fileID++ )
      printf("1 %d %g %g \n",  xyinfo[fileID].id, xyinfo[fileID].x, xyinfo[fileID].y);
  */
  qsort(xyinfo, nfiles, sizeof(xyinfo_t), cmpx);  	      
  /*
  if ( cdoVerbose )
    for ( fileID = 0; fileID < nfiles; fileID++ )
      printf("2 %d %g %g \n",  xyinfo[fileID].id, xyinfo[fileID].x, xyinfo[fileID].y);
  */
  if ( lsouthnorth )
    qsort(xyinfo, nfiles, sizeof(xyinfo_t), cmpxy_lt);  
  else
    qsort(xyinfo, nfiles, sizeof(xyinfo_t), cmpxy_gt);  	      

  if ( cdoVerbose )
    for ( fileID = 0; fileID < nfiles; fileID++ )
      printf("3 %d %g %g \n",  xyinfo[fileID].id, xyinfo[fileID].x, xyinfo[fileID].y);

  nx = 1;
  for ( fileID = 1; fileID < nfiles; fileID++ )
    {
      if ( DBL_IS_EQUAL(xyinfo[0].y, xyinfo[fileID].y) ) nx++;
      else break;
    }
  ny = nfiles/nx;
  if ( cdoVerbose ) cdoPrint("nx %d  ny %d", nx, ny);
  if ( nx*ny != nfiles ) cdoAbort("Number of input files (%d) seems to be incomplete!", nfiles);
 
  xsize2 = 0;
  for ( i = 0; i < nx; ++i ) xsize2 += xsize[xyinfo[i].id];
  ysize2 = 0;
  for ( j = 0; j < ny; ++j ) ysize2 += ysize[xyinfo[j*nx].id];
  if ( cdoVerbose ) cdoPrint("xsize2 %d  ysize2 %d", xsize2, ysize2);

  xvals2 = (double*) malloc(xsize2*sizeof(double));
  yvals2 = (double*) malloc(ysize2*sizeof(double));

  xoff = (int*) malloc((nx+1)*sizeof(int));
  yoff = (int*) malloc((ny+1)*sizeof(int));

  xoff[0] = 0;
  for ( i = 0; i < nx; ++i )
    {
      idx = xyinfo[i].id;
      memcpy(xvals2+xoff[i], xvals[idx], xsize[idx]*sizeof(double));
      xoff[i+1] = xoff[i] + xsize[idx];
    }

  yoff[0] = 0;
  for ( j = 0; j < ny; ++j )
    {
      idx = xyinfo[j*nx].id;
      memcpy(yvals2+yoff[j], yvals[idx], ysize[idx]*sizeof(double));
      yoff[j+1] = yoff[j] + ysize[idx];
    }

  if ( gridindex != NULL )
    {
      for ( fileID = 0; fileID < nfiles; fileID++ )
	{
	  idx = xyinfo[fileID].id;
	  iy = fileID/nx;
	  ix = fileID - iy*nx;

          offset = yoff[iy]*xsize2 + xoff[ix];
	  /*
	  printf("fileID %d %d, iy %d, ix %d, offset %d\n",
		 fileID, xyinfo[fileID].id, iy, ix, offset);
	  */
	  ij = 0;
	  for ( j = 0; j < ysize[idx]; ++j )
	    for ( i = 0; i < xsize[idx]; ++i )
	      {
		gridindex[idx][ij++] = offset+j*xsize2+i;
	      }
	}
    }

  gridID2 = gridCreate(gridtype, xsize2*ysize2);
  gridDefXsize(gridID2, xsize2);
  gridDefYsize(gridID2, ysize2);
  gridDefXvals(gridID2, xvals2);
  gridDefYvals(gridID2, yvals2);

  free(xoff);
  free(yoff);
  free(xsize);
  free(ysize);
  free(xvals2);
  free(yvals2);

  char string[1024];
  string[0] = 0;
  gridID = vlistGrid(ef[0].vlistID, igrid);
  gridInqXname(gridID, string);
  gridDefXname(gridID2, string);
  gridInqYname(gridID, string);
  gridDefYname(gridID2, string);
  gridInqXlongname(gridID, string);
  gridDefXlongname(gridID2, string);
  gridInqYlongname(gridID, string);
  gridDefYlongname(gridID2, string);
  gridInqXunits(gridID, string);
  gridDefXunits(gridID2, string);
  gridInqYunits(gridID, string);
  gridDefYunits(gridID2, string);

  for ( fileID = 0; fileID < nfiles; fileID++ )
    {
      free(xvals[fileID]);
      free(yvals[fileID]);
    }
  free(xvals);
  free(yvals);
  free(xyinfo);

  return (gridID2);
}


void *Collgrid(void *argument)
{
  int varID, recID;
  int nrecs, nrecs0;
  int levelID;
  int nmiss;
  int taxisID1, taxisID2;
  double missval;
  int fileID;

  cdoInitialize(argument);
    
  int nfiles = cdoStreamCnt() - 1;
  const char *ofilename = cdoStreamName(nfiles)->args;

  if ( !cdoSilentMode && !cdoOverwriteMode )
    if ( fileExists(ofilename) )
      if ( !userFileOverwrite(ofilename) )
	cdoAbort("Outputfile %s already exists!", ofilename);

  ens_file_t *ef = (ens_file_t*) malloc(nfiles*sizeof(ens_file_t));

  for ( fileID = 0; fileID < nfiles; fileID++ )
    {
      ef[fileID].streamID = streamOpenRead(cdoStreamName(fileID));
      ef[fileID].vlistID  = streamInqVlist(ef[fileID].streamID);
    }

  int vlistID1 = ef[0].vlistID;
  vlistClearFlag(vlistID1);

  /* check that the contents is always the same */
  for ( fileID = 1; fileID < nfiles; fileID++ )
    vlistCompare(vlistID1, ef[fileID].vlistID, CMP_NAME | CMP_NLEVEL);

  int nvars = vlistNvars(vlistID1);
  int *vars  = (int*) malloc(nvars*sizeof(int));
  for ( varID = 0; varID < nvars; varID++ ) vars[varID] = FALSE;
  int *vars1  = (int*) malloc(nvars*sizeof(int));
  for ( varID = 0; varID < nvars; varID++ ) vars1[varID] = FALSE;

  int nsel = operatorArgc();
  if ( nsel == 0 )
    {
      for ( varID = 0; varID < nvars; varID++ ) vars1[varID] = TRUE;
    }
  else
    {
      char **argnames = operatorArgv();

      if ( cdoVerbose )
	for ( int i = 0; i < nsel; i++ )
	  fprintf(stderr, "name %d = %s\n", i+1, argnames[i]);

      int *selfound = (int*) malloc(nsel*sizeof(int));
      for ( int i = 0; i < nsel; i++ ) selfound[i] = FALSE;

      char varname[CDI_MAX_NAME];
      for ( varID = 0; varID < nvars; varID++ )
	{
	  vlistInqVarName(vlistID1, varID, varname);

	  for ( int isel = 0; isel < nsel; isel++ )
	    {
	      if ( strcmp(argnames[isel], varname) == 0 )
		{
		  selfound[isel] = TRUE;
		  vars1[varID] = TRUE;
		}
	    }
	}

      for ( int isel = 0; isel < nsel; isel++ )
	if ( selfound[isel] == FALSE )
	  cdoAbort("Variable name %s not found!", argnames[isel]);

      free(selfound);
    }

  for ( varID = 0; varID < nvars; varID++ )
    {
      if ( vars1[varID] == TRUE )
	{
	  int zaxisID  = vlistInqVarZaxis(vlistID1, varID);
	  int nlevs    = zaxisInqSize(zaxisID);
	  for ( int levID = 0; levID < nlevs; levID++ )
	    vlistDefFlag(vlistID1, varID, levID, TRUE);
	}
    }

  int gridsize;
  int gridsizemax = 0;
  for ( fileID = 0; fileID < nfiles; fileID++ )
    {
      gridsize = vlistGridsizeMax(ef[fileID].vlistID);
      if ( gridsize > gridsizemax ) gridsizemax = gridsize;
    }
  gridsize = gridsizemax;

  for ( fileID = 0; fileID < nfiles; fileID++ )
    ef[fileID].array = (double*) malloc(gridsizemax*sizeof(double));


  int vlistID2 = vlistCreate();
  vlistCopyFlag(vlistID2, vlistID1);
  /*
  if ( cdoVerbose )
    {
      vlistPrint(vlistID1);
      vlistPrint(vlistID2);
    }
  */
  //int vlistID2 = vlistDuplicate(vlistID1);
  int nvars2 = vlistNvars(vlistID2);
  // int *vars  = (int*) malloc(nvars*sizeof(int));
  //for ( varID = 0; varID < nvars; varID++ ) vars[varID] = FALSE;

  int ngrids1 = vlistNgrids(vlistID1);
  int ngrids2 = vlistNgrids(vlistID2);

  int *gridIDs = (int*) malloc(ngrids2*sizeof(int));
  int **gridindex = (int **) malloc(nfiles*sizeof(int *));
  for ( fileID = 0; fileID < nfiles; fileID++ )
    gridindex[fileID] = (int*) malloc(gridsizemax*sizeof(int));

  int ginit = FALSE;
  for ( int i2 = 0; i2 < ngrids2; ++i2 )
    {
      int i1;
      for ( i1 = 0; i1 < ngrids1; ++i1 )
	if ( vlistGrid(vlistID1, i1) == vlistGrid(vlistID2, i2) ) break;

      //   printf("i1 %d i2 %d\n", i1, i2);

      if ( ginit == FALSE )
	{
	  gridIDs[i2] = genGrid(nfiles, ef, gridindex, i1);
	  if ( gridIDs[i2] != -1 ) ginit = TRUE;
	}
      else
	gridIDs[i2] = genGrid(nfiles, ef, NULL, i1);
    }


  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int gridsize2 = 0;
  for ( int i = 0; i < ngrids2; ++i )
    {
      if ( gridIDs[i] != -1 ) 
	{
	  if ( gridsize2 == 0 ) gridsize2 = gridInqSize(gridIDs[i]);
	  if ( gridsize2 != gridInqSize(gridIDs[i]) ) cdoAbort("gridsize differ!");
	  vlistChangeGridIndex(vlistID2, i, gridIDs[i]);
	}
    }

  for ( varID = 0; varID < nvars2; varID++ )
    {
      int gridID = vlistInqVarGrid(vlistID2, varID);

      for ( int i = 0; i < ngrids2; ++i )
	{
	  if ( gridIDs[i] != -1 ) 
	    {
	      if ( gridID == vlistGrid(vlistID2, i) )
		vars[varID] = TRUE;
	      break;
	    }
	}
    }

  int streamID2 = streamOpenWrite(cdoStreamName(nfiles), cdoFiletype());
      
  streamDefVlist(streamID2, vlistID2);
	  
  double *array2 = (double*) malloc(gridsize2*sizeof(double));

  int tsID = 0;
  do
    {
      nrecs0 = streamInqTimestep(ef[0].streamID, tsID);
      for ( fileID = 1; fileID < nfiles; fileID++ )
	{
	  nrecs = streamInqTimestep(ef[fileID].streamID, tsID);
	  if ( nrecs != nrecs0 )
	    cdoAbort("Number of records at time step %d of %s and %s differ!", tsID+1, cdoStreamName(0)->args, cdoStreamName(fileID)->args);
	}

      taxisCopyTimestep(taxisID2, taxisID1);

      if ( nrecs0 > 0 ) streamDefTimestep(streamID2, tsID);
      
      for ( recID = 0; recID < nrecs0; recID++ )
	{
	  streamInqRecord(ef[0].streamID, &varID, &levelID);
	  if ( cdoVerbose && tsID == 0 ) printf(" tsID, recID, varID, levelID %d %d %d %d\n", tsID, recID, varID, levelID);

	  for ( fileID = 0; fileID < nfiles; fileID++ )
	    {
	      int varIDx, levelIDx;
	      if ( fileID > 0 ) streamInqRecord(ef[fileID].streamID, &varIDx, &levelIDx);
	    }

	  if ( vlistInqFlag(vlistID1, varID, levelID) == TRUE )
	    {
	      int varID2   = vlistFindVar(vlistID2, varID);
	      int levelID2 = vlistFindLevel(vlistID2, varID, levelID);
	      if ( cdoVerbose && tsID == 0 ) printf("varID %d %d levelID %d %d\n", varID, varID2, levelID, levelID2);

	      missval = vlistInqVarMissval(vlistID2, varID2);
	      for ( int i = 0; i < gridsize2; i++ ) array2[i] = missval;

#if defined(_OPENMP)
#pragma omp parallel for default(shared) private(fileID, nmiss)
#endif
	      for ( fileID = 0; fileID < nfiles; fileID++ )
		{
		  streamReadRecord(ef[fileID].streamID, ef[fileID].array, &nmiss);

		  if ( vars[varID2] )
		    {
		      gridsize = gridInqSize(vlistInqVarGrid(ef[fileID].vlistID, varID));
		      for ( int i = 0; i < gridsize; ++i )
			array2[gridindex[fileID][i]] = ef[fileID].array[i];
		    }
		}

	      streamDefRecord(streamID2, varID2, levelID2);

	      if ( vars[varID2] )
		{
		  nmiss = 0;
		  for ( int i = 0; i < gridsize2; i++ )
		    if ( DBL_IS_EQUAL(array2[i], missval) ) nmiss++;

		  streamWriteRecord(streamID2, array2, nmiss);
		}
	      else
		streamWriteRecord(streamID2, ef[0].array, 0);
	    }
	}

      tsID++;
    }
  while ( nrecs0 > 0 );

  for ( fileID = 0; fileID < nfiles; fileID++ )
    streamClose(ef[fileID].streamID);

  streamClose(streamID2);

  for ( fileID = 0; fileID < nfiles; fileID++ )
    if ( ef[fileID].array ) free(ef[fileID].array);

  if ( ef ) free(ef);
  if ( array2 ) free(array2);

  free(gridIDs);
  if ( vars   ) free(vars);
  if ( vars1  ) free(vars1);

  cdoFinish();

  return (0);
}
