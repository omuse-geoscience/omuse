/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2015 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied 1warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

/*
   This module contains the following operators:

      Change_e5slm      change_e5slm          Change ECHAM5 sea land mask
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Change_e5slm(void *argument)
{
  int streamIDslm, streamID1, streamID2;
  const char *fn_slm;
  char name[CDI_MAX_NAME];
  int nrecs, code;
  int tsID, recID, varID, levelID;
  int gridsize;
  int vlistIDslm, vlistID1, vlistID2 = -1;
  int nmiss, nvars;
  int taxisID1, taxisID2;
  short *codes = NULL;
  short *lsea = NULL;
  long i;
  double minval, maxval;
  double *array = NULL;
  double *cland = NULL;

  cdoInitialize(argument);

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  taxisID1 = vlistInqTaxis(vlistID1);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  vlistID2 = vlistDuplicate(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  streamDefVlist(streamID2, vlistID2);


  /* get filename of SLM */
  operatorInputArg("filename of the sea land mask");
  operatorCheckArgc(1);
  fn_slm = operatorArgv()[0];

  /* read SLM */
  argument_t *fileargument = file_argument_new(fn_slm);
  streamIDslm = streamOpenRead(fileargument);
  file_argument_free(fileargument);

  vlistIDslm = streamInqVlist(streamIDslm);

  gridsize = gridInqSize(vlistInqVarGrid(vlistIDslm, 0));

  array = (double*) malloc(gridsize*sizeof(double));
  cland = (double*) malloc(gridsize*sizeof(double));
  lsea  = (short*) malloc(gridsize*sizeof(short));

  streamInqTimestep(streamIDslm, 0);

  streamInqRecord(streamIDslm, &varID, &levelID);
  streamReadRecord(streamIDslm, cland, &nmiss);

  if ( nmiss > 0 ) cdoAbort("SLM with missing values are unsupported!");

  minmaxval(gridsize, cland, NULL, &minval, &maxval);
  if ( minval < 0 || maxval > 1 )
    cdoWarning("Values of SLM out of bounds! (minval=%g, maxval=%g)", minval , maxval);

  streamClose(streamIDslm);

  for ( i = 0; i < gridsize; ++i )
    {
      if ( cland[i] > 0 )
	lsea[i] = FALSE;
      else
	lsea[i] = TRUE;
    }


  nvars = vlistNvars(vlistID1);
  codes = (short*) malloc(nvars*sizeof(short));

  for ( varID = 0; varID < nvars; ++varID )
    {
      if ( gridsize != gridInqSize(vlistInqVarGrid(vlistID1, varID)) )
	cdoAbort("gridsize differ!");

      code = vlistInqVarCode(vlistID1, varID);
      vlistInqVarName(vlistID1, varID, name);

      if ( code < 0 )
	{
	  if      ( strcmp(name, "SLM")       == 0 ) code = 172;
	  else if ( strcmp(name, "ALAKE")     == 0 ) code = 99;
	  else if ( strcmp(name, "WS")        == 0 ) code = 140;
	  else if ( strcmp(name, "AZ0")       == 0 ) code = 173;
	  else if ( strcmp(name, "ALB")       == 0 ) code = 174;
	  else if ( strcmp(name, "VGRAT")     == 0 ) code = 198;
	  else if ( strcmp(name, "FOREST")    == 0 ) code = 212;
	  else if ( strcmp(name, "FAO")       == 0 ) code = 226;
	  else if ( strcmp(name, "WSMX")      == 0 ) code = 229;
	  else if ( strcmp(name, "GLAC")      == 0 ) code = 232;
	  else if ( strcmp(name, "VLTCLIM")   == 0 ) code = 71;
	  else if ( strcmp(name, "VGRATCLIM") == 0 ) code = 70;
	}

      codes[varID] = code;
    }


  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);
      
      for ( recID = 0; recID < nrecs; recID++ )
	{ 
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamReadRecord(streamID1, array, &nmiss);

	  code = codes[varID];
	  if ( code == 172 )
	    {
	      cdoPrint("SLM changed!");
	      for ( i = 0; i < gridsize; ++i )
		array[i] = cland[i];
	    }
	  else if ( code == 99 )
	    {
	      cdoPrint("ALAKE set all values to zero!");
	      for ( i = 0; i < gridsize; ++i )
		array[i] = 0;
	    }
	  else if ( code == 232 )
	    {
	      cdoPrint("GLAC set sea points to %g!", array[0]);
	      for ( i = 0; i < gridsize; ++i )
		if ( cland[i] < 0.5 ) array[i] = array[0];
	    }
	  else if ( code ==  70 || code ==  71 || code == 140 ||
		    code == 173 || code == 174 || code == 198 ||
		    code == 200 || code == 212 || code == 226 ||
		    code == 229 )
	    {
	      cdoPrint("Code %d set sea points to %g!", code, array[0]);
	      for ( i = 0; i < gridsize; ++i )
		if ( lsea[i] ) array[i] = array[0];
	    }

	  streamDefRecord(streamID2,  varID,  levelID);
	  streamWriteRecord(streamID2, array, nmiss);
	}

      tsID++;
    }
  
  streamClose(streamID1);
  streamClose(streamID2);

  free(array);
  free(cland);
  free(lsea);
  free(codes);

  cdoFinish();

  return (0);
}
