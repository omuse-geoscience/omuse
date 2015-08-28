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


void *Test(void *argument)
{
  /*
  int streamID1, streamID2;
  */

  cdoInitialize(argument);

  /*
  streamID1 = streamOpenRead(cdoStreamName(0));
  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamClose(streamID2);
  streamClose(streamID1);
  */
  cdoFinish();

  return (0);
}


void *Test2(void *argument)
{
  /*
  int streamID1, streamID2, streamID3;
  */

  cdoInitialize(argument);

  /*
  streamID1 = streamOpenRead(cdoStreamName(0));
  streamID2 = streamOpenRead(cdoStreamName(1));
  streamID3 = streamOpenWrite(cdoStreamName(2), cdoFiletype());

  streamClose(streamID3);
  streamClose(streamID2);
  streamClose(streamID1);
  */
  cdoFinish();

  return (0);
}


void *Testdata(void *argument)
{
  int streamID1, streamID2 = CDI_UNDEFID;
  int nrecs;
  int tsID1, tsID2, recID, varID, levelID;
  int gridsize, i;
  int vlistID1, vlistID2 = -1;
  int nmiss;
  int taxisID1, taxisID2 = CDI_UNDEFID;
  double *array = NULL;
  float *fval;
  int *ival;
  unsigned char *cval;
  unsigned char *cval2;
  FILE *fp;

  cdoInitialize(argument);

  tsID2 = 0;

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  taxisID1 = vlistInqTaxis(vlistID1);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  vlistID2 = vlistDuplicate(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  streamDefVlist(streamID2, vlistID2);

  gridsize = vlistGridsizeMax(vlistID1);
  array = (double*) malloc(gridsize*sizeof(double));
  fval = (float*) malloc(gridsize*sizeof(float));
  ival = (int*) malloc(gridsize*sizeof(int));
  cval = (unsigned char*) malloc(gridsize*sizeof(unsigned char)*4);
  cval2 = (unsigned char*) malloc(gridsize*sizeof(unsigned char)*4);

  fp = fopen("testdata", "w");

  tsID1 = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID1)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID2);
	       
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamDefRecord(streamID2,  varID,  levelID);
	  
	  streamReadRecord(streamID1, array, &nmiss);

	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	  for ( i = 0; i < gridsize; ++i )
	    {
	      fval[i] = (float) array[i];

	      memcpy(&ival[i], &fval[i], 4);
	      memcpy(&cval[i*4], &fval[i], 4);

	      cval2[i+gridsize*0] = cval[i*4+0];
	      cval2[i+gridsize*1] = cval[i*4+1];
	      cval2[i+gridsize*2] = cval[i*4+2];
	      cval2[i+gridsize*3] = cval[i*4+3];

	      if ( tsID1 == 0 && recID == 0 )
	      printf("%4d %3d %3d %3d %3d %d %g\n",
		     i, (unsigned int)cval[4*i+0], (unsigned int)cval[4*i+1], (unsigned int)cval[4*i+2], (unsigned int)cval[4*i+3], ival[i], fval[i]);
	    }

	  streamWriteRecord(streamID2, array, nmiss);

	  fwrite(cval, 4, gridsize, fp);
	}

      tsID1++;
      tsID2++;
    }
  
  fclose(fp);
  streamClose(streamID1);
  streamClose(streamID2);

  if ( array ) free(array);

  cdoFinish();

  return (0);
}
