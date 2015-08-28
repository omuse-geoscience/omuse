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

      Setzaxis   setzaxis        Set zaxis
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"

void genLayerBounds(int nlev, double *levels, double *lbounds, double *ubounds);

int getkeyval_dp(const char *keyval, const char *key, double *val)
{
  int status = 0;
  size_t keylen = strlen(key);

  if ( strncmp(keyval, key, keylen) == 0 )
    {
      const char *pkv = keyval+keylen;
      if ( pkv[0] == '=' && pkv[1] != 0 )
        {
          *val = parameter2double(&pkv[1]);
          status = 1;
        }
      else
        {
          cdoAbort("Syntax error for parameter %s!", keyval);
        }
    }

  return status;
}


void *Setzaxis(void *argument)
{
  int nrecs;
  int recID, varID, levelID;
  int zaxisID1, zaxisID2 = -1;
  int nzaxis, index;
  int nmiss;
  int found;
  int lztop = FALSE, lzbot = FALSE;
  double ztop = 0, zbot = 0;

  cdoInitialize(argument);

  int SETZAXIS       = cdoOperatorAdd("setzaxis",        0, 0, "zaxis description file");
  int GENLEVELBOUNDS = cdoOperatorAdd("genlevelbounds",  0, 0, NULL);

  int operatorID = cdoOperatorID();

  if ( operatorID == SETZAXIS )
    {
      operatorInputArg(cdoOperatorEnter(operatorID));  
      operatorCheckArgc(1);
      zaxisID2 = cdoDefineZaxis(operatorArgv()[0]);
    }
  else if ( operatorID == GENLEVELBOUNDS )
    { 
      unsigned npar = operatorArgc();
      char **parnames = operatorArgv();

      for ( unsigned i = 0; i < npar; i++ )
        {
          if ( cdoVerbose ) cdoPrint("keyval[%d]: %s", i+1, parnames[i]);
          
          if      ( !lzbot && getkeyval_dp(parnames[i], "zbot", &zbot) ) lzbot = TRUE;
          else if ( !lztop && getkeyval_dp(parnames[i], "ztop", &ztop) ) lztop = TRUE;
          else cdoAbort("Parameter >%s< unsupported! Supported parameter are: zbot, ztop", parnames[i]);
        }
    }

  int streamID1 = streamOpenRead(cdoStreamName(0));

  int vlistID1 = streamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  if ( operatorID == SETZAXIS )
    {
      found = 0;
      nzaxis = vlistNzaxis(vlistID1);
      for ( index = 0; index < nzaxis; index++ )
	{
	  zaxisID1 = vlistZaxis(vlistID1, index);

	  if ( zaxisInqSize(zaxisID1) == zaxisInqSize(zaxisID2) )
	    {
	      vlistChangeZaxisIndex(vlistID2, index, zaxisID2);
	      found++;
	    }
	}
      if ( ! found ) cdoWarning("No zaxis with %d levels found!", zaxisInqSize(zaxisID2));
    }
  else if ( operatorID == GENLEVELBOUNDS )
    {
      nzaxis = vlistNzaxis(vlistID1);
      for ( index = 0; index < nzaxis; index++ )
	{
	  zaxisID1 = vlistZaxis(vlistID1, index);
          int nlev = zaxisInqSize(zaxisID1);
          double *levels  = (double *) malloc(nlev*sizeof(double));
          double *lbounds = (double *) malloc(nlev*sizeof(double));
          double *ubounds = (double *) malloc(nlev*sizeof(double));

          if ( nlev > 1 )
            {
              zaxisInqLevels(zaxisID1, levels);
              zaxisID2 = zaxisDuplicate(zaxisID1);

              genLayerBounds(nlev, levels, lbounds, ubounds);

              if ( lzbot ) lbounds[0] = zbot;
              if ( lztop ) ubounds[nlev-1] = ztop;
              zaxisDefLbounds(zaxisID2, lbounds);
              zaxisDefUbounds(zaxisID2, ubounds);
	      vlistChangeZaxisIndex(vlistID2, index, zaxisID2);
            }

          free(levels);
          free(lbounds);
          free(ubounds);
	}
    }

  streamDefVlist(streamID2, vlistID2);

  int gridsize = vlistGridsizeMax(vlistID1);
  if ( vlistNumber(vlistID1) != CDI_REAL ) gridsize *= 2;
  double *array = (double*) malloc(gridsize*sizeof(double));

  int tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);
	       
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamDefRecord(streamID2,  varID,  levelID);
	  
	  streamReadRecord(streamID1, array, &nmiss);
	  streamWriteRecord(streamID2, array, nmiss);
	}
      tsID++;
    }

  streamClose(streamID1);
  streamClose(streamID2);

  if ( array ) free(array);

  cdoFinish();

  return (0);
}
