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

      Set        setcode         Set code number
      Set        setparam        Set parameter identifier
      Set        setname         Set variable name
      Set        setlevel        Set level
      Set        setltype        Set GRIB level type
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


int stringToParam(const char *paramstr)
{
  int param = 0;
  int pnum = -1, pcat = 255, pdis = 255;

  sscanf(paramstr, "%d.%d.%d", &pnum, &pcat, &pdis);
  
  if ( cdoVerbose ) cdoPrint("pnum, pcat, pdis: %d.%d.%d", pnum, pcat, pdis);

  param = cdiEncodeParam(pnum, pcat, pdis);

  return (param);
}


void *Set(void *argument)
{
  int SETCODE, SETPARAM, SETNAME, SETUNIT, SETLEVEL, SETLTYPE, SETTABNUM;
  int operatorID;
  int streamID1, streamID2 = CDI_UNDEFID;
  int nrecs, nvars, newval = -1, tabnum = 0;
  int tsID1, recID, varID, levelID;
  int vlistID1, vlistID2;
  int taxisID1, taxisID2;
  int nmiss;
  int gridsize;
  int index, zaxisID1, zaxisID2, nzaxis, nlevs;
  int zaxistype;
  int newparam    = 0;
  char *newname   = NULL, *newunit = NULL;
  double newlevel = 0;
  double *levels  = NULL;
  double *array   = NULL;

  cdoInitialize(argument);

  SETCODE    = cdoOperatorAdd("setcode",    0, 0, "code number");
  SETPARAM   = cdoOperatorAdd("setparam",   0, 0, "parameter identifier (format: code[.tabnum] or num[.cat[.dis]])");
  SETNAME    = cdoOperatorAdd("setname",    0, 0, "variable name");
  SETUNIT    = cdoOperatorAdd("setunit",    0, 0, "variable unit");
  SETLEVEL   = cdoOperatorAdd("setlevel",   0, 0, "level");
  SETLTYPE   = cdoOperatorAdd("setltype",   0, 0, "GRIB level type");
  SETTABNUM  = cdoOperatorAdd("settabnum",  0, 0, "GRIB table number");

  operatorID = cdoOperatorID();

  operatorInputArg(cdoOperatorEnter(operatorID));
  if ( operatorID == SETCODE || operatorID == SETLTYPE )
    {
      newval = parameter2int(operatorArgv()[0]);
    }
  else if ( operatorID == SETPARAM )
    {
      newparam = stringToParam(operatorArgv()[0]);
    }
  else if ( operatorID == SETNAME )
    {
      newname = operatorArgv()[0];
    }
  else if ( operatorID == SETUNIT )
    {
      newunit = operatorArgv()[0];
    }
  else if ( operatorID == SETTABNUM )
    {
      tabnum = parameter2int(operatorArgv()[0]);
    }
  else if ( operatorID == SETLEVEL )
    {
      newlevel = parameter2double(operatorArgv()[0]);
    }

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);
  /* vlistPrint(vlistID2);*/

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  if ( operatorID == SETCODE )
    {
      nvars = vlistNvars(vlistID2);
      for ( varID = 0; varID < nvars; varID++ )
	vlistDefVarCode(vlistID2, varID, newval);
    }
  else if ( operatorID == SETPARAM )
    {
      vlistDefVarParam(vlistID2, 0, newparam);
    }
  else if ( operatorID == SETNAME )
    {
      vlistDefVarName(vlistID2, 0, newname);
    }
  else if ( operatorID == SETUNIT )
    {
      vlistDefVarUnits(vlistID2, 0, newunit);
    }
  else if ( operatorID == SETTABNUM )
    {
      int tableID;
      tableID = tableDef(-1, tabnum, NULL);
      nvars = vlistNvars(vlistID2);
      for ( varID = 0; varID < nvars; varID++ )
	vlistDefVarTable(vlistID2, varID, tableID);
    }
  else if ( operatorID == SETLEVEL )
    {
      nzaxis = vlistNzaxis(vlistID2);
      for ( index = 0; index < nzaxis; index++ )
	{
	  zaxisID1 = vlistZaxis(vlistID2, index);
	  zaxisID2 = zaxisDuplicate(zaxisID1);
	  nlevs = zaxisInqSize(zaxisID2);
	  levels = (double*) malloc(nlevs*sizeof(double));
	  zaxisInqLevels(zaxisID2, levels);
	  levels[0] = newlevel;
	  zaxisDefLevels(zaxisID2, levels);
	  vlistChangeZaxis(vlistID2, zaxisID1, zaxisID2);
	  free(levels);
	}
    }
  else if ( operatorID == SETLTYPE )
    {
      nzaxis = vlistNzaxis(vlistID2);
      for ( index = 0; index < nzaxis; index++ )
	{
	  zaxisID1 = vlistZaxis(vlistID2, index);
	  zaxisID2 = zaxisDuplicate(zaxisID1);

	  zaxistype = ZAXIS_GENERIC;
	  zaxisChangeType(zaxisID2, zaxistype);
	  zaxisDefLtype(zaxisID2, newval);
	  vlistChangeZaxis(vlistID2, zaxisID1, zaxisID2);
	}
    }

  /* vlistPrint(vlistID2);*/
  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  gridsize = vlistGridsizeMax(vlistID1);
  if ( vlistNumber(vlistID1) != CDI_REAL ) gridsize *= 2;
  array = (double*) malloc(gridsize*sizeof(double));

  tsID1 = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID1)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID1);
	       
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamDefRecord(streamID2,  varID,  levelID);
	  
	  streamReadRecord(streamID1, array, &nmiss);
	  streamWriteRecord(streamID2, array, nmiss);
	}
      tsID1++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  if ( array ) free(array);

  cdoFinish();

  return (0);
}
