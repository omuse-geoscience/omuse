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

      Change     chcode          Change code number
      Change     chtabnum        Change GRIB1 parameter table number
      Change     chparam         Change parameter identifier
      Change     chname          Change variable name
      Change     chlevel         Change level
      Change     chlevelc        Change level of one code
      Change     chlevelv        Change level of one variable
      Change     chltype         Change GRIB level type 
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"

int stringToParam(const char *paramstr);

#define  MAXARG     16384

void *Change(void *argument)
{
  int CHCODE, CHTABNUM, CHPARAM, CHNAME, CHUNIT, CHLEVEL, CHLEVELC, CHLEVELV, CHLTYPE;  
  int operatorID;
  int streamID1, streamID2 = CDI_UNDEFID;
  int nrecs, nvars;
  int tsID1, recID, varID = 0, levelID;
  int vlistID1, vlistID2;
  int taxisID1, taxisID2;
  int chints[MAXARG], nch = 0;
  char *chnames[MAXARG];
  char varname[CDI_MAX_NAME];
  char *chname = NULL;
  int chcode = 0;
  int param;
  int code, tabnum, i;
  int nmiss;
  int gridsize;
  int nfound;
  int nzaxis, zaxisID1, zaxisID2, k, nlevs, index; 
  double chlevels[MAXARG];
  int  chltypes[MAXARG];              
  double *levels = NULL;
  double *newlevels = NULL;
  double *array = NULL;

  cdoInitialize(argument);

  CHCODE   = cdoOperatorAdd("chcode",   0, 0, "pairs of old and new code numbers");
  CHTABNUM = cdoOperatorAdd("chtabnum", 0, 0, "pairs of old and new GRIB1 table numbers");
  CHPARAM  = cdoOperatorAdd("chparam",  0, 0, "pairs of old and new parameter identifiers");
  CHNAME   = cdoOperatorAdd("chname",   0, 0, "pairs of old and new variable names");
  CHUNIT   = cdoOperatorAdd("chunit",   0, 0, "pairs of old and new variable units");
  CHLEVEL  = cdoOperatorAdd("chlevel",  0, 0, "pairs of old and new levels");
  CHLEVELC = cdoOperatorAdd("chlevelc", 0, 0, "code number, old and new level");
  CHLEVELV = cdoOperatorAdd("chlevelv", 0, 0, "variable name, old and new level");
  CHLTYPE  = cdoOperatorAdd("chltype",  0, 0, "pairs of old and new type");          

  operatorID = cdoOperatorID();

  operatorInputArg(cdoOperatorEnter(operatorID));

  nch = operatorArgc();

  if ( operatorID == CHCODE || operatorID == CHTABNUM )
    {
      if ( nch%2 ) cdoAbort("Odd number of input arguments!");
      for ( i = 0; i < nch; i++ )
	chints[i] = parameter2int(operatorArgv()[i]);
    }
  else if ( operatorID == CHPARAM || operatorID == CHNAME || operatorID == CHUNIT )
    {
      if ( nch%2 ) cdoAbort("Odd number of input arguments!");
      for ( i = 0; i < nch; i++ )
	chnames[i] = operatorArgv()[i];
    }
  else if ( operatorID == CHLEVEL )
    {
      if ( nch%2 ) cdoAbort("Odd number of input arguments!");
      for ( i = 0; i < nch; i++ )
	chlevels[i] = parameter2double(operatorArgv()[i]);
    }
  else if ( operatorID == CHLEVELC )
    {
      operatorCheckArgc(3);
      
      chcode = parameter2int(operatorArgv()[0]);
      chlevels[0] = parameter2double(operatorArgv()[1]);
      chlevels[1] = parameter2double(operatorArgv()[2]);
    }
  else if ( operatorID == CHLEVELV )
    {
      operatorCheckArgc(3);
      
      chname = operatorArgv()[0];
      chlevels[0] = parameter2double(operatorArgv()[1]);
      chlevels[1] = parameter2double(operatorArgv()[2]);
    }
  else if ( operatorID == CHLTYPE )                  
    {
      if ( nch%2 ) cdoAbort("Odd number of input arguments!");
      for ( i = 0; i < nch; i++ )
	chltypes[i] = parameter2int(operatorArgv()[i]);
    }

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  if ( operatorID == CHCODE )
    {
      nvars = vlistNvars(vlistID2);
      for ( varID = 0; varID < nvars; varID++ )
	{
	  code = vlistInqVarCode(vlistID2, varID);
	  for ( i = 0; i < nch; i += 2 )
	    if ( code == chints[i] )
	      vlistDefVarCode(vlistID2, varID, chints[i+1]);
	}
    }
  else if ( operatorID == CHTABNUM )
    {
      int tableID;
      nvars = vlistNvars(vlistID2);
      for ( varID = 0; varID < nvars; varID++ )
	{
	  tabnum = tableInqNum(vlistInqVarTable(vlistID2, varID));
	  for ( i = 0; i < nch; i += 2 )
	    if ( tabnum == chints[i] )
	      {
		tableID = tableDef(-1, chints[i+1], NULL);
		vlistDefVarTable(vlistID2, varID, tableID);
	      }
	}
    }
  else if ( operatorID == CHPARAM )
    {
      nvars = vlistNvars(vlistID2);
      for ( varID = 0; varID < nvars; varID++ )
	{
	  param = vlistInqVarParam(vlistID2, varID);
	  if ( cdoVerbose )
	    {
	      int pnum, pcat, pdis;
	      cdiDecodeParam(param, &pnum, &pcat, &pdis);
	      cdoPrint("pnum, pcat, pdis: %d.%d.%d", pnum, pcat, pdis);
	    }
	  for ( i = 0; i < nch; i += 2 )
	    if ( param == stringToParam(chnames[i]) )
	      vlistDefVarParam(vlistID2, varID, stringToParam(chnames[i+1]));
	}
    }
  else if ( operatorID == CHNAME )
    {
      nvars = vlistNvars(vlistID2);
      for ( varID = 0; varID < nvars; varID++ )
	{
	  vlistInqVarName(vlistID2, varID, varname);
	  for ( i = 0; i < nch; i += 2 )
	    if ( strcmp(varname, chnames[i]) == 0 )
	      vlistDefVarName(vlistID2, varID, chnames[i+1]);
	}
    }
  else if ( operatorID == CHUNIT )
    {
      nvars = vlistNvars(vlistID2);
      for ( varID = 0; varID < nvars; varID++ )
	{

	  vlistInqVarUnits(vlistID2, varID, varname);
	  for ( i = 0; i < nch; i += 2 )
	    if ( strcmp(varname, chnames[i]) == 0 )
	      vlistDefVarUnits(vlistID2, varID, chnames[i+1]);
	}
    }
  else if ( operatorID == CHLEVEL )
    {
      nzaxis = vlistNzaxis(vlistID2);
      for ( index = 0; index < nzaxis; index++ )
	{
	  zaxisID1 = vlistZaxis(vlistID2, index);
	  nlevs = zaxisInqSize(zaxisID1);
	  levels = (double*) malloc(nlevs*sizeof(double));
	  newlevels = (double*) malloc(nlevs*sizeof(double));
	  zaxisInqLevels(zaxisID1, levels);

	  for ( k = 0; k < nlevs; k++ ) newlevels[k] = levels[k];

	  nfound = 0;
	  for ( i = 0; i < nch; i += 2 )
	    for ( k = 0; k < nlevs; k++ )
	      if ( fabs(levels[k] - chlevels[i]) < 0.0001 ) nfound++;

	  if ( nfound )
	    {
	      zaxisID2 = zaxisDuplicate(zaxisID1);
	      for ( i = 0; i < nch; i += 2 )
		for ( k = 0; k < nlevs; k++ )
		  if ( fabs(levels[k] - chlevels[i]) < 0.001 )
		    newlevels[k] = chlevels[i+1];

	      zaxisDefLevels(zaxisID2, newlevels);
	      vlistChangeZaxis(vlistID2, zaxisID1, zaxisID2);
	    }

	  free(levels);
	  free(newlevels);
	}
    }
  else if ( operatorID == CHLEVELC || operatorID == CHLEVELV )
    {
      nvars = vlistNvars(vlistID2);
      if ( operatorID == CHLEVELC )
	{
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      code = vlistInqVarCode(vlistID2, varID);
	      if ( code == chcode ) break;
	    }
	  if ( varID == nvars ) cdoAbort("Code %d not found!", chcode);
	}
      else
	{
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      vlistInqVarName(vlistID2, varID, varname);
	      if ( strcmp(varname, chname) == 0 ) break;
	    }
	  if ( varID == nvars ) cdoAbort("Variable name %s not found!", chname);
	}

      zaxisID1 = vlistInqVarZaxis(vlistID2, varID);
      nlevs = zaxisInqSize(zaxisID1);
      levels = (double*) malloc(nlevs*sizeof(double));
      zaxisInqLevels(zaxisID1, levels);
      nfound = 0;
      for ( k = 0; k < nlevs; k++ )
	if ( fabs(levels[k] - chlevels[0]) < 0.0001 ) nfound++;

      if ( nfound )
	{
	  zaxisID2 = zaxisDuplicate(zaxisID1);
	  for ( k = 0; k < nlevs; k++ )
	    if ( fabs(levels[k] - chlevels[0]) < 0.001 )
	      levels[k] = chlevels[1];

	  zaxisDefLevels(zaxisID2, levels);
	  vlistChangeVarZaxis(vlistID2, varID, zaxisID2);
	}
      else
	cdoAbort("Level %g not found!", chlevels[0]);

      free(levels);
    }
  else if ( operatorID == CHLTYPE )                
    {
      int ltype, ltype1, ltype2;

      nzaxis = vlistNzaxis(vlistID2);
      for ( index = 0; index < nzaxis; index++ )
	{
	  zaxisID1 = vlistZaxis(vlistID2, index);
	  zaxisID2 = zaxisDuplicate(zaxisID1);

	  ltype = zaxisInqLtype(zaxisID1);

	  for ( i = 0; i < nch; i += 2 )
	    {
	      ltype1 = chltypes[i];
	      ltype2 = chltypes[i+1];

	      if ( ltype1 == ltype )
		{
		  zaxisChangeType(zaxisID2, ZAXIS_GENERIC);
		  zaxisDefLtype(zaxisID2, ltype2);
		  vlistChangeZaxis(vlistID2, zaxisID1, zaxisID2);
		}
	    }
	}
    }

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  gridsize = vlistGridsizeMax(vlistID2);
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

  streamClose(streamID1);
  streamClose(streamID2);

  if ( array ) free(array);

  cdoFinish();

  return (0);
}
