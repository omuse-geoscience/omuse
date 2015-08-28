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

*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "list.h"


void *Histogram(void *argument)
{
  int HISTCOUNT, HISTSUM, HISTMEAN, HISTFREQ;
  int operatorID;
  int streamID1, streamID2;
  int nrecs;
  int tsID1, recID, varID, levelID;
  int gridsize;
  int vlistID1, vlistID2;
  int nmiss;
  int taxisID1, taxisID2 = CDI_UNDEFID;
  int nbins;
  int i, nvars;
  int offset;
  int nzaxis, nlevel, zaxisID, zaxisID2, index;
  double *array = NULL;
  double *fltarr = NULL;
  double *bins;
  double missval;
  LIST *flist = listNew(FLT_LIST);
  double **vardata = NULL;
  double **varcount = NULL;
  double **vartcount = NULL;

  cdoInitialize(argument);

  HISTCOUNT = cdoOperatorAdd("histcount", 0, 0, NULL);
  HISTSUM   = cdoOperatorAdd("histsum",   0, 0, NULL);
  HISTMEAN  = cdoOperatorAdd("histmean",  0, 0, NULL);
  HISTFREQ  = cdoOperatorAdd("histfreq",  0, 0, NULL);

  UNUSED(HISTSUM);

  operatorID = cdoOperatorID();

  operatorInputArg("bins");

  nbins = args2fltlist(operatorArgc(), operatorArgv(), flist) - 1;
  if ( nbins < 1 ) cdoAbort("Too few arguments!");
  fltarr = (double *) listArrayPtr(flist);

  if ( cdoVerbose )
    {
      printf("nbins = %d\n", nbins);
      for ( i = 0; i < nbins; i++ )
	printf("flt %d = %g\n", i+1, fltarr[i]);
    }

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  taxisID1 = vlistInqTaxis(vlistID1);

  vlistID2 = vlistDuplicate(vlistID1);

  /* create zaxis for output bins */
  zaxisID2 = zaxisCreate(ZAXIS_GENERIC, nbins);
  bins = (double*) malloc(nbins*sizeof(double));
  /* for ( i = 0; i < nbins; i++ ) bins[i] = (fltarr[i]+fltarr[i+1])/2; */
  for ( i = 0; i < nbins; i++ ) bins[i] = fltarr[i];
  zaxisDefLevels(zaxisID2, bins);
  free(bins);
  zaxisDefLbounds(zaxisID2, fltarr);
  zaxisDefUbounds(zaxisID2, fltarr+1);
  zaxisDefName(zaxisID2, "bin");
  zaxisDefLongname(zaxisID2, "histogram bins");
  zaxisDefUnits(zaxisID2, "level");

  /* check zaxis: only 2D fields allowed */
  nzaxis = vlistNzaxis(vlistID1);
  for ( index = 0; index < nzaxis; index++ )
    {
      zaxisID = vlistZaxis(vlistID1, index);
      nlevel = zaxisInqSize(zaxisID);
      if ( nlevel > 1 )
	cdoAbort("Found 3D field with %d levels. Only 2D fields allowed!", nlevel);
      vlistChangeZaxisIndex(vlistID2, index, zaxisID2);
    }

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  streamDefVlist(streamID2, vlistID2);

  nvars = vlistNvars(vlistID2);
  vardata   = (double **) malloc(nvars*sizeof(double *));
  varcount  = (double **) malloc(nvars*sizeof(double *));
  vartcount = (double **) malloc(nvars*sizeof(double *));
  for ( varID = 0; varID < nvars; varID++ )
    {
      gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
      vardata[varID]  = (double*) malloc(nbins*gridsize*sizeof(double));
      varcount[varID] = (double*) malloc(nbins*gridsize*sizeof(double));
      vartcount[varID] = (double*) malloc(gridsize*sizeof(double));
      memset(vardata[varID], 0, nbins*gridsize*sizeof(double));
      memset(varcount[varID], 0, nbins*gridsize*sizeof(double));
      memset(vartcount[varID], 0, gridsize*sizeof(double));
    }

  gridsize = vlistGridsizeMax(vlistID1);
  array = (double*) malloc(gridsize*sizeof(double));

  tsID1 = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID1)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);
	       
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamReadRecord(streamID1, array, &nmiss);
	  missval = vlistInqVarMissval(vlistID1, varID);

	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));

	  nmiss=0;
	  for ( i = 0; i < gridsize; i++ )
	    {
	      if ( !DBL_IS_EQUAL(array[i], missval) )
		{
		  vartcount[varID][i] += 1;
		  index = 0;
		  while( index < nbins )
		    {
		      offset = gridsize*index;
		      if ( !DBL_IS_EQUAL(vardata[varID][offset+i], missval) &&
			   array[i] >= fltarr[index] && array[i] < fltarr[index+1] )
			{
			  vardata[varID][offset+i]  += array[i];
			  varcount[varID][offset+i] += 1;
			  break;
			}
		      index++;
		    }
		}
	      else { /* missing value */
		nmiss++;
	      }
	    }
	}
      tsID1++;
    }


  streamDefTimestep(streamID2, 0);

  for ( varID = 0; varID < nvars; varID++ )
    {
      missval = vlistInqVarMissval(vlistID2, varID);
      gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));

      /* fix mising values */
      
      for ( index = 0; index < nbins; index++ )
	{
	  nmiss = 0;
	  offset = gridsize*index;

	  for ( i = 0; i < gridsize; i++ )
	    {
	      if ( vartcount[varID][i] > 0 )
		{
		  if ( operatorID == HISTMEAN || operatorID == HISTFREQ )
		    {
		      if ( varcount[varID][offset+i] > 0 ) 
			{
			  if ( operatorID == HISTMEAN )
			    vardata[varID][offset+i] /= varcount[varID][offset+i];	    
			  else 
			    vardata[varID][offset+i] = varcount[varID][offset+i] / vartcount[varID][i];
			} 
		    }
		}
	      else
		{
		  nmiss++;
		  varcount[varID][offset+i] = missval;
		  vardata[varID][offset+i] = missval;
		}
	    }

	  streamDefRecord(streamID2,  varID,  index);

	  if ( operatorID == HISTCOUNT )
	    streamWriteRecord(streamID2, varcount[varID]+offset, nmiss);
	  else
	    streamWriteRecord(streamID2, vardata[varID]+offset, nmiss);
	}
    }
  
  streamClose(streamID1);
  streamClose(streamID2);

  if ( vardata )
    {
      for ( varID = 0; varID < nvars; varID++ )
	{
	  free(vardata[varID]);
	  free(varcount[varID]);
	  free(vartcount[varID]);
	}

      free(vardata);
      free(varcount);
      free(vartcount);
    }

  if ( array ) free(array);

  listDelete(flist);

  cdoFinish();

  return (0);
}
