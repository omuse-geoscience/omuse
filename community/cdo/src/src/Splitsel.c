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

      Splitsel   splitsel        Split time selection
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Splitsel(void *argument)
{
  int gridsize;
  //int vdate = 0, vtime = 0;
  int nrecs = 0;
  int varID, levelID, recID;
  int tsID, tsID2;
  int nsets;
  int streamID1, streamID2;
  int vlistID1, vlistID2, taxisID1, taxisID2;
  int nmiss;
  int gridID;
  int nvars, nlevel;
  int nconst;
/*   int ndates = 0, noffset = 0, nskip = 0, nargc; */
  double ndates, noffset, nskip;
  int i2 = 0;
  int nargc;

  /* from Splittime.c */
  int nchars;
  char filesuffix[32];
  char filename[8192];
  const char *refname;
  int index = 0;
  int lcopy = FALSE;
  double *array = NULL;
  field_t **vars = NULL;

  cdoInitialize(argument);

  if ( processSelf() != 0 ) cdoAbort("This operator can't be combined with other operators!");

  cdoOperatorAdd("splitsel",  0,  0, NULL);

  if ( UNCHANGED_RECORD ) lcopy = TRUE;

  /*  operatorInputArg("nsets <noffset <nskip>>"); */

  nargc = operatorArgc();

  if ( nargc < 1 )
    cdoAbort("Too few arguments! Need %d found %d.", 1, nargc);

/*   ndates = parameter2int(operatorArgv()[0]); */
/*   if ( nargc > 1 ) noffset = parameter2int(operatorArgv()[1]); */
/*   if ( nargc > 2 ) nskip   = parameter2int(operatorArgv()[2]); */
/*   printf("%s %s %s\n", operatorArgv()[0],operatorArgv()[1],operatorArgv()[2]); */
  ndates = noffset = nskip = 0.0;
  ndates = parameter2double(operatorArgv()[0]);
  if ( nargc > 1 ) noffset = parameter2double(operatorArgv()[1]);
  if ( nargc > 2 ) nskip   = parameter2double(operatorArgv()[2]);

  if ( cdoVerbose ) cdoPrint("nsets = %f, noffset = %f, nskip = %f", ndates, noffset, nskip);

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
/*   taxisID2 = taxisCreate(TAXIS_ABSOLUTE); */
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  strcpy(filename, cdoStreamName(1)->args);
  nchars = strlen(filename);

  refname = cdoStreamName(0)->argv[cdoStreamName(0)->argc-1];
  filesuffix[0] = 0;
  cdoGenFileSuffix(filesuffix, sizeof(filesuffix), streamInqFiletype(streamID1), vlistID1, refname);

  //  if ( ! lcopy )
    {
      gridsize = vlistGridsizeMax(vlistID1);
      if ( vlistNumber(vlistID1) != CDI_REAL ) gridsize *= 2;
      array = (double*) malloc(gridsize*sizeof(double));
    }

  nvars = vlistNvars(vlistID1);
  nconst = 0;
  for ( varID = 0; varID < nvars; varID++ )
    if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) nconst++;

  if ( nconst )
    {
      vars = (field_t **) malloc(nvars*sizeof(field_t *));

      for ( varID = 0; varID < nvars; varID++ )
	{
	  if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT )
	    {
	      gridID  = vlistInqVarGrid(vlistID1, varID);
	      nlevel  = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      gridsize = gridInqSize(gridID);
		  
	      vars[varID] = (field_t*) malloc(nlevel*sizeof(field_t));

	      for ( levelID = 0; levelID < nlevel; levelID++ )
		{
		  field_init(&vars[varID][levelID]);
		  vars[varID][levelID].grid    = gridID;
		  vars[varID][levelID].ptr     = (double*) malloc(gridsize*sizeof(double));
		}
	    }
	}
    }


  /* offset */
  for ( tsID = 0; tsID < noffset; tsID++ )
    {
      nrecs = streamInqTimestep(streamID1, tsID);
      if ( nrecs == 0 )
	{
	  cdoWarning("noffset is larger than number of timesteps!");
	  goto LABEL_END;
	}

      if ( tsID == 0 && nconst )
	for ( recID = 0; recID < nrecs; recID++ )
	  {
	    streamInqRecord(streamID1, &varID, &levelID);
	    if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT )
	      streamReadRecord(streamID1, vars[varID][levelID].ptr, &vars[varID][levelID].nmiss);
	  }
    }

  index = 0;
  nsets = 0;
  while ( TRUE )
    {
      sprintf(filename+nchars, "%06d", index);
      sprintf(filename+nchars+6, "%s", filesuffix);
	  
      if ( cdoVerbose ) cdoPrint("create file %s", filename);
      argument_t *fileargument = file_argument_new(filename);
      streamID2 = streamOpenWrite(fileargument, cdoFiletype());
      file_argument_free(fileargument);

      streamDefVlist(streamID2, vlistID2);

      tsID2 = 0;

      for ( ; nsets < (int)(ndates*(index+1)); nsets++ ) 
	{
	  nrecs = streamInqTimestep(streamID1, tsID);
	  if ( nrecs == 0 ) break;

	  /*
	  vdate = taxisInqVdate(taxisID1);
	  vtime = taxisInqVtime(taxisID1);
	  printf("vdate: %d vtime: %d\n", vdate, vtime);
	   */

	  taxisCopyTimestep(taxisID2, taxisID1);
	  streamDefTimestep(streamID2, tsID2);

	  if ( tsID > 0 && tsID2 == 0 && nconst )
	    {
	      for ( varID = 0; varID < nvars; varID++ )
		{
		  if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT )
		    {
		      nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
		      for ( levelID = 0; levelID < nlevel; levelID++ )
			{
			  streamDefRecord(streamID2, varID, levelID);
			  nmiss = vars[varID][levelID].nmiss;
			  streamWriteRecord(streamID2, vars[varID][levelID].ptr, nmiss);
			}
		    }
		}
	    }

	  for ( recID = 0; recID < nrecs; recID++ )
	    {
	      
	      streamInqRecord(streamID1, &varID, &levelID);
	      streamDefRecord(streamID2,  varID,  levelID);
	      if ( lcopy && !(tsID == 0 && nconst) )
		{
		  streamCopyRecord(streamID2, streamID1);
		}
	      else
		{
		  streamReadRecord(streamID1, array, &nmiss);
		  streamWriteRecord(streamID2, array, nmiss);

		  if ( tsID == 0 && nconst )
		    {
		      if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT )
			{
			  gridID  = vlistInqVarGrid(vlistID1, varID);
			  gridsize = gridInqSize(gridID);
			  memcpy(vars[varID][levelID].ptr, array, gridsize*sizeof(double));
			  vars[varID][levelID].nmiss = nmiss;
			}
		    }
		}
	    }
	  
	  tsID++;
	  tsID2++;	  
	}
      
      streamClose(streamID2);
      if ( nrecs == 0 ) break;

      nrecs = streamInqTimestep(streamID1, tsID);
      if ( nrecs == 0 ) break;

      for ( ; i2 < (int)(nskip*(index+1)); i2++ )
	{
	  nrecs = streamInqTimestep(streamID1, tsID);
	  if ( nrecs == 0 ) break;
	  tsID++;
	}

      nrecs = streamInqTimestep(streamID1, tsID);
      if ( nrecs == 0 ) break;

      index++;
    }

 LABEL_END:

  streamClose(streamID1);
 
  if ( array ) free(array);

  if ( nconst )
    {
      for ( varID = 0; varID < nvars; varID++ )
	{
	  if ( vlistInqVarTsteptype(vlistID2, varID) == TSTEP_CONSTANT )
	    {
	      nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
	      for ( levelID = 0; levelID < nlevel; levelID++ )
		if ( vars[varID][levelID].ptr )
		  free(vars[varID][levelID].ptr);

	      free(vars[varID]);
	    }
	}

      if ( vars  ) free(vars);
    }

  vlistDestroy(vlistID2);

  cdoFinish();

  return (0);
}
