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

      Selvar     selparam        Select parameters by identifier (format: code.tabnum  or  pnum.cat.dis)
      Selvar     delparam        Delete parameters by identifier (format: code.tabnum  or  pnum.cat.dis)
      Selvar     selcode         Select parameters by code number
      Selvar     delcode         Delete parameters by code number
      Selvar     selname         Select parameters by name
      Selvar     delname         Delete parameters by name
      Selvar     selstdname      Select parameters by CF standard name
      Selvar     sellevel        Select levels
      Selvar     sellevidx       Select levels by index
      Selvar     selgrid         Select grids
      Selvar     selzaxis        Select zaxis
      Selvar     seltabnum       Select parameter table number
      Selvar     selltype        Select GRIB level type 
*/

#include <ctype.h>  /* isdigit */

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "error.h"
#include "util.h"
#include "list.h"


void *Selvar(void *argument)
{
  int nlevs;
  int code, tabnum, param, gridID, zaxisID, levID;
  int grididx, zaxisidx;
  double level;
  int varID2, levelID2;
  int recID, varID, levelID;
  int *intarr = NULL, nsel = 0;
  int *selfound = NULL;
  double *fltarr = NULL;
  char paramstr[32];
  char varname[CDI_MAX_NAME];
  char stdname[CDI_MAX_NAME];
  char gridname[CDI_MAX_NAME];
  char zaxistypename[CDI_MAX_NAME];
  char zaxisname[CDI_MAX_NAME];
  char **argnames = NULL;
  int isel;
  int i;
  int lcopy = FALSE;
  int gridsize;
  int nmiss;
  int gridnum = 0;
  double *array = NULL;
  LIST *ilist = listNew(INT_LIST);
  LIST *flist = listNew(FLT_LIST);

  cdoInitialize(argument);

# define INVERTS_SELECTION(id) (cdoOperatorF2(id) & 1)
# define TAKES_STRINGS(id) (cdoOperatorF2(id) & 2)
# define TAKES_INTEGERS(id) (cdoOperatorF2(id) & 4)
# define TAKES_FLOATS(id) (cdoOperatorF2(id) & 8)

  int SELPARAM     = cdoOperatorAdd("selparam",     0, 2,   "parameters");
  int SELCODE      = cdoOperatorAdd("selcode",      0, 4,   "code numbers");
  int SELNAME      = cdoOperatorAdd("selname",      0, 2,   "variable names");
  int SELSTDNAME   = cdoOperatorAdd("selstdname",   0, 2,   "standard names");
  int SELLEVEL     = cdoOperatorAdd("sellevel",     0, 8,   "levels");
  int SELLEVIDX    = cdoOperatorAdd("sellevidx",    0, 4,   "index of levels");
  int SELGRID      = cdoOperatorAdd("selgrid",      0, 4|2, "list of grid names or numbers");
  int SELZAXIS     = cdoOperatorAdd("selzaxis",     0, 4|2, "list of zaxis types or numbers");
  int SELZAXISNAME = cdoOperatorAdd("selzaxisname", 0, 2,   "list of zaxis names");
  int SELTABNUM    = cdoOperatorAdd("seltabnum",    0, 4,   "table numbers");
  int DELPARAM     = cdoOperatorAdd("delparam",     0, 2|1, "parameter");
  int DELCODE      = cdoOperatorAdd("delcode",      0, 1,   "code numbers");
  int DELNAME      = cdoOperatorAdd("delname",      0, 2|1, "variable names");
  int SELLTYPE     = cdoOperatorAdd("selltype",     0, 4,   "GRIB level types"); 

  if ( UNCHANGED_RECORD ) lcopy = TRUE;

  int operatorID = cdoOperatorID();

  operatorInputArg(cdoOperatorEnter(operatorID));

  int args_are_numeric = operatorArgc() > 0 && isdigit(*operatorArgv()[0]);

  if ( TAKES_STRINGS(operatorID) && !( TAKES_INTEGERS(operatorID) && args_are_numeric ) )
    {
      nsel     = operatorArgc();
      argnames = operatorArgv();

      if ( cdoVerbose )
	for ( i = 0; i < nsel; i++ )
	  fprintf(stderr, "name %d = %s\n", i+1, argnames[i]);
    }
  else if ( TAKES_FLOATS(operatorID) )
    {
      nsel = args2fltlist(operatorArgc(), operatorArgv(), flist);
      fltarr = (double *) listArrayPtr(flist);

      if ( cdoVerbose )
	for ( i = 0; i < nsel; i++ )
	  printf("flt %d = %g\n", i+1, fltarr[i]);
    }
  else
    {
      nsel = args2intlist(operatorArgc(), operatorArgv(), ilist);
      intarr = (int *) listArrayPtr(ilist);

      if ( cdoVerbose )
	for ( i = 0; i < nsel; i++ )
	  printf("int %d = %d\n", i+1, intarr[i]);
    }

  if ( nsel )
    {
      selfound = (int*) malloc(nsel*sizeof(int));
      for ( i = 0; i < nsel; i++ ) selfound[i] = FALSE;
    }

  /*
  if ( nsel == 0 )
    cdoAbort("missing code argument!");
  */
  int streamID1 = streamOpenRead(cdoStreamName(0));

  int vlistID1 = streamInqVlist(streamID1);
  int nvars = vlistNvars(vlistID1);

  if ( operatorID == SELGRID && !args_are_numeric && nsel == 1 && strncmp(argnames[0], "var=", 4) == 0 )
    {
      const char *gridvarname = argnames[0] + 4;

      if ( *gridvarname == 0 ) cdoAbort("Variable name missing!", gridvarname);
      
      for ( varID = 0; varID < nvars; varID++ )
        {
          vlistInqVarName(vlistID1, varID, varname);
          if ( strcmp(varname, gridvarname) == 0 )
            {
              gridID = vlistInqVarGrid(vlistID1, varID);
              gridnum = 1 + vlistGridIndex(vlistID1, gridID);
              args_are_numeric = TRUE;
              intarr = &gridnum;
              break;
            }
        }

      if ( !gridnum ) cdoAbort("Variable %s not found!", gridvarname);
    }

  vlistClearFlag(vlistID1);
  for ( varID = 0; varID < nvars; varID++ )
    {
      vlistInqVarName(vlistID1, varID, varname);
      vlistInqVarStdname(vlistID1, varID, stdname);
      param    = vlistInqVarParam(vlistID1, varID);
      code     = vlistInqVarCode(vlistID1, varID);
      tabnum   = tableInqNum(vlistInqVarTable(vlistID1, varID));
      gridID   = vlistInqVarGrid(vlistID1, varID);
      grididx  = vlistGridIndex(vlistID1, gridID);
      zaxisID  = vlistInqVarZaxis(vlistID1, varID);
      zaxisidx = vlistZaxisIndex(vlistID1, zaxisID);
      nlevs    = zaxisInqSize(zaxisID);
      gridName(gridInqType(gridID), gridname);
      zaxisInqName(zaxisID, zaxisname);
      zaxisName(zaxisInqType(zaxisID), zaxistypename);

      cdiParamToString(param, paramstr, sizeof(paramstr));

      for ( levID = 0; levID < nlevs; levID++ )
	{
	  level = zaxisInqLevel(zaxisID, levID);

	  if ( operatorID == DELCODE || operatorID == DELNAME || operatorID == DELPARAM )
	    vlistDefFlag(vlistID1, varID, levID, TRUE);

	  for ( isel = 0; isel < nsel; isel++ )
	    {
              int found = 0;
	      if ( operatorID == SELCODE )
		{
		  found = intarr[isel] == code;
		}
	      else if ( operatorID == SELPARAM )
		{
		  found = wildcardmatch(argnames[isel], paramstr) == 0;
		}
	      else if ( operatorID == SELNAME )
		{
		  found = wildcardmatch(argnames[isel], varname) == 0;
		}
	      else if ( operatorID == SELSTDNAME )
		{
		  found = strcmp(argnames[isel], stdname) == 0;
		}
	      else if ( operatorID == SELLEVEL )
		{
		  found = fabs(fltarr[isel] - level) < 0.0001;
		}
	      else if ( operatorID == SELLEVIDX )
		{
		  found = intarr[isel] == (levID+1);
		}
	      else if ( operatorID == SELGRID && args_are_numeric )
		{
		  found = intarr[isel] == (grididx+1);
		}
	      else if ( operatorID == SELGRID && !args_are_numeric )
		{
		  found = memcmp(argnames[isel], gridname, strlen(argnames[isel])) == 0;
		}
	      else if ( operatorID == SELZAXIS && args_are_numeric )
		{
		  found = intarr[isel] == (zaxisidx+1);
		}
	      else if ( operatorID == SELZAXIS && !args_are_numeric )
		{
		  found = memcmp(argnames[isel], zaxistypename, strlen(argnames[isel])) == 0;
		}
	      else if ( operatorID == SELZAXISNAME )
		{
		  found = wildcardmatch(argnames[isel], zaxisname) == 0;
		}
	      else if ( operatorID == SELTABNUM )
		{
		  found = intarr[isel] == tabnum;
		}
	      else if ( operatorID == DELCODE )
		{
		  found = intarr[isel] == code;
		}
	      else if ( operatorID == DELNAME )
		{
		  found = wildcardmatch(argnames[isel], varname) == 0;
		}
	      else if ( operatorID == DELPARAM )
		{
		  found = strcmp(argnames[isel], paramstr) == 0;
		}
	      else if ( operatorID == SELLTYPE )
		{
		  found = intarr[isel] == zaxis2ltype(zaxisID);
		}

	      if ( found )
	        {
		  vlistDefFlag(vlistID1, varID, levID, !INVERTS_SELECTION(operatorID));
		  selfound[isel] = TRUE;
	        }

	    }
	}
    }

  int npar = 0;
  for ( varID = 0; varID < nvars; varID++ )
    {
      zaxisID = vlistInqVarZaxis(vlistID1, varID);
      nlevs   = zaxisInqSize(zaxisID);

      for ( levID = 0; levID < nlevs; levID++ )
	if ( vlistInqFlag(vlistID1, varID, levID) == TRUE ) break;
	      
      if ( levID < nlevs ) npar++;
    }

  for ( isel = 0; isel < nsel; isel++ )
    {
      if ( selfound[isel] == FALSE )
	{
	  if ( operatorID == SELCODE || operatorID == DELCODE )
	    {
	      cdoWarning("Code number %d not found!", intarr[isel]);
	    }
	  else if ( operatorID == SELPARAM || operatorID == DELPARAM )
	    {
	      cdoWarning("Parameter %s not found!", argnames[isel]);
	    }
	  else if ( operatorID == SELNAME || operatorID == DELNAME )
	    {
	      cdoWarning("Variable name %s not found!", argnames[isel]);
	    }
	  else if ( operatorID == SELSTDNAME )
	    {
	      cdoWarning("Variable with standard name %s not found!", argnames[isel]);
	    }
	  else if ( operatorID == SELLEVEL )
	    {
	      cdoWarning("Level %g not found!", fltarr[isel]);
	    }
	  else if ( operatorID == SELLEVIDX )
	    {
	      cdoWarning("Level index %d not found!", intarr[isel]);
	    }
	  else if ( operatorID == SELGRID && args_are_numeric )
	    {
	      cdoWarning("Grid %d not found!", intarr[isel]);
	    }
	  else if ( operatorID == SELGRID && !args_are_numeric )
	    {
	      cdoWarning("Grid name %s not found!", argnames[isel]);
	    }
	  else if ( operatorID == SELZAXIS && args_are_numeric )
	    {
	      cdoWarning("Zaxis %d not found!", intarr[isel]);
	    }
	  else if ( operatorID == SELZAXIS && !args_are_numeric )
	    {
	      cdoWarning("Zaxis type %s not found!", argnames[isel]);
	    }
	  else if ( operatorID == SELZAXISNAME )
	    {
	      cdoWarning("Zaxis name %s not found!", argnames[isel]);
	    }
	  else if ( operatorID == SELTABNUM )
	    {
	      cdoWarning("Table number %d not found!", intarr[isel]);
	    }
	  else if ( operatorID == SELLTYPE )
	    {
	      cdoWarning("GRIB level type %d not found!", intarr[isel]);
	    }
	}
    }

  if ( npar == 0 ) cdoAbort("No variables selected!");

  int vlistID2 = vlistCreate();
  vlistCopyFlag(vlistID2, vlistID1);

  nvars = vlistNvars(vlistID2);
  for ( varID = 0; varID < nvars; ++varID )
    if ( vlistInqVarTsteptype(vlistID2, varID) != TSTEP_CONSTANT ) break;
  if ( varID == nvars ) vlistDefNtsteps(vlistID2, 0);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int nrecs = vlistNrecs(vlistID2);

  int streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  if ( ! lcopy )
    {
      gridsize = vlistGridsizeMax(vlistID1);
      if ( vlistNumber(vlistID1) != CDI_REAL ) gridsize *= 2;
      array = (double*) malloc(gridsize*sizeof(double));
    }

  int tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);
     
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  if ( vlistInqFlag(vlistID1, varID, levelID) == TRUE )
	    {
	      varID2   = vlistFindVar(vlistID2, varID);
	      levelID2 = vlistFindLevel(vlistID2, varID, levelID);

	      streamDefRecord(streamID2, varID2, levelID2);
	      if ( lcopy )
		{
		  streamCopyRecord(streamID2, streamID1);
		}
	      else
		{
		  streamReadRecord(streamID1, array, &nmiss);
		  streamWriteRecord(streamID2, array, nmiss);
		}
     	    }
	}
       
      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);
 
  vlistDestroy(vlistID2);

  if ( ! lcopy )
    if ( array ) free(array);

  if ( selfound ) free(selfound);

  listDelete(ilist);
  listDelete(flist);

  cdoFinish();

  return (NULL);
}
