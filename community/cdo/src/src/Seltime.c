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

      Seltime    seltimestep     Select timesteps
      Seltime    seltime         Select times
      Seltime    selhour         Select hours
      Seltime    selday          Select days
      Seltime    selmon          Select months
      Seltime    selyear         Select years
      Seltime    selseas         Select seasons
      Seltime    seldate         Select dates
      Seltime    selsmon         Select single month
*/

#include <ctype.h>

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "error.h"
#include "util.h"
#include "list.h"


#define  NOPERATORS  32

static
int seaslist(LIST *ilist)
{
  int i;
  char Seas[3];
  int seas[4] = {FALSE, FALSE, FALSE, FALSE};
  int imon[17]; /* 1-16 ! */
  int ival;
  size_t len;

  int season_start = get_season_start();
  int nsel = operatorArgc();
  if ( isdigit(*operatorArgv()[0]))
    for ( i = 0; i < nsel; i++ )
      {
	ival = parameter2int(operatorArgv()[i]);
	if      ( ival == 1 || ival == 13 ) seas[0] = TRUE;
	else if ( ival == 2 || ival == 14 ) seas[1] = TRUE;
	else if ( ival == 3 || ival == 15 ) seas[2] = TRUE;
	else if ( ival == 4 || ival == 16 ) seas[3] = TRUE;
	else cdoAbort("Season %d not available!", ival);
      }
  else
    for ( i = 0; i < nsel; i++ )
      {
	len = strlen(operatorArgv()[i]);
	if ( len > 3 ) len = 3;
	while ( len-- > 0 ) Seas[len] = toupper(operatorArgv()[i][len]);
	if ( season_start == START_DEC )
	  {
	    if      ( memcmp(Seas, "DJF", 3) == 0 ) seas[0] = TRUE;
	    else if ( memcmp(Seas, "MAM", 3) == 0 ) seas[1] = TRUE;
	    else if ( memcmp(Seas, "JJA", 3) == 0 ) seas[2] = TRUE;
	    else if ( memcmp(Seas, "SON", 3) == 0 ) seas[3] = TRUE;
	    else cdoAbort("Season %s not available!", operatorArgv()[i]);
	  }
	else
	  {
	    if      ( memcmp(Seas, "JFM", 3) == 0 ) seas[0] = TRUE;
	    else if ( memcmp(Seas, "AMJ", 3) == 0 ) seas[1] = TRUE;
	    else if ( memcmp(Seas, "JAS", 3) == 0 ) seas[2] = TRUE;
	    else if ( memcmp(Seas, "OND", 3) == 0 ) seas[3] = TRUE;
	    else cdoAbort("Season %s not available!", operatorArgv()[i]);
	  }
      }
  
  for ( i = 0; i < 17; ++i ) imon[i] = 0;
  
  if ( season_start == START_DEC )
    {
      if ( seas[0] ) { imon[12]++; imon[ 1]++; imon[ 2]++; imon[13]++; }
      if ( seas[1] ) { imon[ 3]++; imon[ 4]++; imon[ 5]++; imon[14]++; }
      if ( seas[2] ) { imon[ 6]++; imon[ 7]++; imon[ 8]++; imon[15]++; }
      if ( seas[3] ) { imon[ 9]++; imon[10]++; imon[11]++; imon[16]++; }
    }
  else
    {
      if ( seas[0] ) { imon[ 1]++; imon[ 2]++; imon[ 3]++; imon[13]++; }
      if ( seas[1] ) { imon[ 4]++; imon[ 5]++; imon[ 6]++; imon[14]++; }
      if ( seas[2] ) { imon[ 7]++; imon[ 8]++; imon[ 9]++; imon[15]++; }
      if ( seas[3] ) { imon[10]++; imon[11]++; imon[12]++; imon[16]++; }
    }
  
  nsel = 0;
  for ( i = 1; i < 17; ++i )
    {
      if ( imon[i] )
	listSetInt(ilist, nsel++, i);
    }

  return (nsel);
}


double datestr_to_double(const char *datestr)
{
  int year = 1, month = 1, day = 1, hour = 0, minute = 0, second = 0;
  double fval = 0;

  if ( strchr(datestr, '-') == NULL )
    {
      fval = parameter2double(datestr);
    }
  else if ( strchr(datestr, 'T') )
    {
      sscanf(datestr, "%d-%d-%dT%d:%d:%d", &year, &month, &day, &hour, &minute, &second);
      fval = cdiEncodeTime(hour, minute, second);
      if ( fabs(fval) > 0 ) fval /= 1000000;
      fval += cdiEncodeDate(year, month, day);
    }
  else
    {
      sscanf(datestr, "%d-%d-%d", &year, &month, &day);
      fval = cdiEncodeDate(year, month, day);
    }

  return fval;
}

static
int datelist(LIST *flist)
{
  int set2 = TRUE;
  double fval = 0;

  int nsel = operatorArgc();
  if ( nsel < 1 ) cdoAbort("Too few arguments!");
  if ( nsel > 2 ) cdoAbort("Too many arguments!");

  for ( int i = 0; i < nsel; i++ )
    {
      if ( operatorArgv()[i][0] == '-' && operatorArgv()[i][1] == 0 )
	{
	  if ( i == 0 )
	    fval = -99999999999.;
	  else
	    fval =  99999999999.;

	  listSetFlt(flist, i,  fval);
	}
      else
	{
	  fval = datestr_to_double(operatorArgv()[i]);
	  if ( strchr(operatorArgv()[i], 'T') )
	    set2 = FALSE;
	  else if ( nsel > 1 && i > 0 )
	    fval += 0.999;
	}

      listSetFlt(flist, i, fval);
    }

  if ( nsel == 1 && set2 == TRUE )
    {
      fval += 0.999;
      listSetFlt(flist, nsel, fval);
      nsel = 2;
    }

  return (nsel);
}


void *Seltime(void *argument)
{
  int streamID2 = -1;
  int nrecs;
  int recID, varID, levelID;
  int nsel = 0, selival;
  int vdate, vtime;
  int copytimestep;
  int copy_nts2 = FALSE;
  int i, isel;
  int gridsize;
  int nmiss;
  int ncts = 0, nts, it;
  int *selfound = NULL;
  int hour = 0, minute = 0, second = 0;
  int nts1 = 0, nts2 = 0;
  int its1 = 0, its2 = 0;
  double selfval = 0;
  double *array = NULL;
  LIST *ilist = listNew(INT_LIST);
  LIST *flist = listNew(FLT_LIST);
  int gridID;
  int nlevel;
  int lconstout = FALSE;
  int process_nts1 = FALSE, process_nts2 = FALSE;
  int *vdate_list = NULL, *vtime_list = NULL;
  double *single;
  field_t ***vars = NULL;

  cdoInitialize(argument);

  int SELTIMESTEP = cdoOperatorAdd("seltimestep", func_step,     1, "timesteps");
  int SELDATE     = cdoOperatorAdd("seldate",     func_datetime, 1, "start date and end date (format YYYY-MM-DDThh:mm:ss)");
  int SELTIME     = cdoOperatorAdd("seltime",     func_time,     1, "times (format hh:mm:ss)");
  int SELHOUR     = cdoOperatorAdd("selhour",     func_time, 10000, "hours");
  int SELDAY      = cdoOperatorAdd("selday",      func_date,     1, "days");
  int SELMON      = cdoOperatorAdd("selmon",      func_date,   100, "months");
  int SELYEAR     = cdoOperatorAdd("selyear",     func_date, 10000, "years");
  int SELSEAS     = cdoOperatorAdd("selseas",     func_date,   100, "seasons");
  int SELSMON     = cdoOperatorAdd("selsmon",     func_date,   100, "month[,nts1[,nts2]]");

  int moddat[NOPERATORS];
  moddat[SELTIMESTEP] =          1;
  /*  moddat[SELDATE]     = 1000000000; */
  moddat[SELDATE]     =          0;
  moddat[SELTIME]     =    1000000;
  moddat[SELHOUR]     =      10000;
  moddat[SELDAY]      =        100;
  moddat[SELMON]      =        100;
  moddat[SELYEAR]     = 1000000000;
  moddat[SELSEAS]     =        100;
  moddat[SELSMON]     =        100;

  int operatorID = cdoOperatorID();

  int operfunc = cdoOperatorF1(operatorID);
  int intval   = cdoOperatorF2(operatorID);

  int lcopy = FALSE;
  if ( UNCHANGED_RECORD ) lcopy = TRUE;

  operatorInputArg(cdoOperatorEnter(operatorID));

  if ( operatorID == SELSEAS )
    {
      nsel = seaslist(ilist);
    }
  else if ( operatorID == SELDATE )
    {
      nsel = datelist(flist);
    }
  else if ( operatorID == SELTIME )
    {
      nsel = operatorArgc();
      if ( nsel < 1 ) cdoAbort("Too few arguments!");
      for ( i = 0; i < nsel; i++ )
	{
	  if ( strchr(operatorArgv()[i], ':') )
	    {
	      sscanf(operatorArgv()[i], "%d:%d:%d", &hour, &minute, &second);
	      listSetInt(ilist, i, cdiEncodeTime(hour, minute, second));
	    }
	  else
	    {
	      listSetInt(ilist, i, parameter2int(operatorArgv()[i]));
	    }
	}
    }
  else
    {
      nsel = args2intlist(operatorArgc(), operatorArgv(), ilist);
    }

  if ( nsel < 1 ) cdoAbort("No timestep selected!");

  int *intarr = (int *) listArrayPtr(ilist);
  double *fltarr = (double *) listArrayPtr(flist);

  if ( operatorID == SELSMON )
    {
      if ( nsel > 1 ) nts1 = intarr[1];
      if ( nsel > 2 ) nts2 = intarr[2];
      else            nts2 = nts1;

      if ( nsel > 3 ) cdoAbort("Too many parameters");

      if ( cdoVerbose )
	cdoPrint("mon=%d  nts1=%d  nts2=%d", intarr[0], nts1, nts2);

      nsel = 1;
    }

  if ( nsel )
    {
      selfound = (int*) malloc(nsel*sizeof(int));
      for ( i = 0; i < nsel; i++ ) selfound[i] = FALSE;
    }

  int streamID1 = streamOpenRead(cdoStreamName(0));

  int vlistID1 = streamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  if ( nsel == 1 && operfunc == func_step )  vlistDefNtsteps(vlistID2,  1);
  else                                       vlistDefNtsteps(vlistID2, -1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  if ( ! lcopy )
    {
      gridsize = vlistGridsizeMax(vlistID1);
      if ( vlistNumber(vlistID1) != CDI_REAL ) gridsize *= 2;
      array = (double*) malloc(gridsize*sizeof(double));
    }

  int ntsteps = vlistNtsteps(vlistID1);

  /* add support for negative timestep values */
  if ( operatorID == SELTIMESTEP && ntsteps > 0 )
    {
      for ( i = 0; i < nsel; i++ )
	{
	  if ( intarr[i] < 0 )
	    {
	      if ( cdoVerbose )
		cdoPrint("timestep %d changed to %d", intarr[i], ntsteps + 1 + intarr[i]);
	      intarr[i] = ntsteps + 1 + intarr[i];
	    }
	}
    }

  if ( cdoVerbose )
    {
      for ( i = 0; i < nsel; i++ )
	if ( operatorID == SELDATE )
	  cdoPrint("fltarr entry: %d %14.4f", i+1, fltarr[i]);
	else
	  cdoPrint("intarr entry: %d %d", i+1, intarr[i]);
    }

  int nvars = vlistNvars(vlistID1);
  int nconst = 0;
  for ( varID = 0; varID < nvars; varID++ )
    if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) nconst++;
      
  int lnts1 = operatorID == SELSMON && nts1 > 0;

  if ( lnts1 || nconst )
    {
      if ( lnts1 )
	{
	  vdate_list = (int*) malloc(nts1*sizeof(int));
	  vtime_list = (int*) malloc(nts1*sizeof(int));
	}
      else
	{
	  nts1 = 1;
	}

      vars  = (field_t ***) malloc(nts1*sizeof(field_t **));

      for ( int tsID = 0; tsID < nts1; tsID++ )
	{
	  vars[tsID] = field_malloc(vlistID1, FIELD_NONE);

	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      if ( lnts1 || (vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT) )
		{
		  gridID  = vlistInqVarGrid(vlistID1, varID);
		  nlevel  = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
		  gridsize = gridInqSize(gridID);
		  
		  for ( levelID = 0; levelID < nlevel; levelID++ )
		    {
		      vars[tsID][varID][levelID].ptr = (double*) malloc(gridsize*sizeof(double));
		    }
		}
	    }
	}
    }

  int tsID  = 0;
  int tsID2 = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      vdate = taxisInqVdate(taxisID1);
      vtime = taxisInqVtime(taxisID1);

      copytimestep = FALSE;
      selival = -1;

      if ( operfunc == func_step )
	{
	  selival = tsID + 1;
	  if ( selival > intarr[nsel-1] ) break;
	}
      else if ( operfunc == func_date )
	{
	  selival = (vdate/intval)%moddat[operatorID];
	}
      else if ( operfunc == func_time )
	{
	  selival = (vtime/intval)%moddat[operatorID];
	}
      else if ( operfunc == func_datetime )
	{
	  selfval = vdate + vtime/1000000.;
	}

      if ( operatorID == SELDATE )
	{
	  if ( selfval >= fltarr[0] && selfval <= fltarr[nsel-1] )
	    {
	      copytimestep = TRUE;
	      selfound[0]      = TRUE;
	      selfound[nsel-1] = TRUE;
	    }
	}
      else
	{
	  for ( i = 0; i < nsel; i++ )
	    if ( selival == intarr[i] )
	      {
		copytimestep = TRUE;
		selfound[i] = TRUE;
		break;
	      }
	}

      if ( operatorID == SELSMON && copytimestep == FALSE )
	{
	  copy_nts2 = FALSE;

	  if ( process_nts1 == TRUE )
	    {
	      process_nts2 = TRUE;
	      its2 = 0;
	      process_nts1 = FALSE;
	    }

	  if ( process_nts2 == TRUE )
	    {
	      if ( its2++ < nts2 )
		copy_nts2 = TRUE;
	      else
		process_nts2 = FALSE;
	    }
	}

      if ( copytimestep || copy_nts2 )
	{
	  if ( tsID2 == 0 )
	    {
	      streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

	      streamDefVlist(streamID2, vlistID2);
	    }

	  if ( lnts1 && ncts == 0 )
	    {
	      nts = nts1;
	      if ( its1 < nts1 )
		{
		  nts = its1;
		  cdoWarning("%d timesteps missing before month %d!", nts1-its1, intarr[0]);
		}

	      for ( it = 0; it < nts; it++ )
		{
		  taxisDefVdate(taxisID2, vdate_list[it]);
		  taxisDefVtime(taxisID2, vtime_list[it]);
		  streamDefTimestep(streamID2, tsID2++);
		  
		  for ( varID = 0; varID < nvars; varID++ )
		    {
		      if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT && tsID2 > 1 ) continue;
		      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
		      for ( levelID = 0; levelID < nlevel; levelID++ )
			{
			  streamDefRecord(streamID2, varID, levelID);
			  single = vars[it][varID][levelID].ptr;
			  nmiss  = vars[it][varID][levelID].nmiss;
			  streamWriteRecord(streamID2, single, nmiss);
			}
		    }
		}

	      its1 = 0;
	    }

	  ncts++;
	  if ( process_nts2 == FALSE )
	    {
	      its2 = 0;
	      process_nts1 = TRUE;
	    }

	  taxisCopyTimestep(taxisID2, taxisID1);

	  streamDefTimestep(streamID2, tsID2++);

	  if ( tsID > 0 && lconstout )
	    {
	      lconstout = FALSE;
	      nts = nts1 - 1;
	      for ( varID = 0; varID < nvars; varID++ )
		{
		  if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT )
		    {
		      nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
		      for ( levelID = 0; levelID < nlevel; levelID++ )
			{
			  streamDefRecord(streamID2, varID, levelID);
			  single = vars[nts][varID][levelID].ptr;
			  nmiss  = vars[nts][varID][levelID].nmiss;
			  streamWriteRecord(streamID2, single, nmiss);
			}
		    }
		}
	    }

	  for ( recID = 0; recID < nrecs; recID++ )
	    {
	      streamInqRecord(streamID1, &varID, &levelID);
	      streamDefRecord(streamID2, varID, levelID);
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
      else
	{
	  ncts = 0;

	  if ( lnts1 || tsID == 0 )
	    {
	      if ( tsID == 0 && nconst && (!lnts1) ) lconstout = TRUE;

	      nts = nts1-1;
	      if ( lnts1 )
		{
		  if ( its1 <= nts )
		    nts = its1;
		  else
		    for ( it = 0; it < nts; it++ )
		      {
			vdate_list[it] = vdate_list[it+1];
			vtime_list[it] = vtime_list[it+1];
			for ( varID = 0; varID < nvars; varID++ )
			  {
			    if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;
			    gridID   = vlistInqVarGrid(vlistID1, varID);
			    gridsize = gridInqSize(gridID);
			    nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
			    for ( levelID = 0; levelID < nlevel; levelID++ )
			      {
				memcpy(vars[it][varID][levelID].ptr,
				       vars[it+1][varID][levelID].ptr,
				       gridsize*sizeof(double));
				vars[it][varID][levelID].nmiss = vars[it+1][varID][levelID].nmiss;
			      }
			  }
		      }

		  vdate_list[nts] = taxisInqVdate(taxisID1);
		  vtime_list[nts] = taxisInqVtime(taxisID1);

		  its1++;
		}

	      for ( recID = 0; recID < nrecs; recID++ )
		{
		  streamInqRecord(streamID1, &varID, &levelID);
		  if ( lnts1 || (vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT) )
		    {
		      single = vars[nts][varID][levelID].ptr;
		      streamReadRecord(streamID1, single, &nmiss);
		      vars[nts][varID][levelID].nmiss = nmiss;
		    }
		}
	    }
	}
       
      tsID++;
    }

  if ( streamID2 != -1 ) streamClose(streamID2);
  streamClose(streamID1);
 
  if ( operatorID == SELSMON )
    if ( its2 < nts2 )
      cdoWarning("%d timesteps missing after the last month!", nts2-its2);

  if ( ! lcopy )
    if ( array ) free(array);

  for ( isel = 0; isel < nsel; isel++ )
    {
      if ( selfound[isel] == FALSE )
	{
	 
	  int isel2;
	  int lcont = FALSE;
	  for ( isel2 = isel+1; isel2 < nsel; isel2++ )
	    if ( selfound[isel2] == TRUE ) break;
	  if ( isel2 == nsel && (nsel-isel) > 1 ) lcont = TRUE;

	  if ( operatorID == SELTIMESTEP )
	    {
	      int lcont2 = FALSE;
	      if ( lcont )
		{    
		  for ( isel2 = isel+1; isel2 < nsel; isel2++ )
		    if ( intarr[isel2-1] != intarr[isel2]-1 ) break;
		  if ( isel2 == nsel ) lcont2 = TRUE;
		}

	      if ( lcont2 )
		{
		  cdoWarning("Time steps %d-%d not found!", intarr[isel], intarr[nsel-1]);
		  break;
		}
	      else
		cdoWarning("Time step %d not found!", intarr[isel]);
	    }
	  else if ( operatorID == SELDATE )
	    {
	      if ( isel == 0 )
		cdoWarning("Date between %14.4f and %14.4f not found!", fltarr[0], fltarr[nsel-1]);
	    }
	  else if ( operatorID == SELTIME )
	    {
	      cdoWarning("Time %d not found!", intarr[isel]);
	    }
	  else if ( operatorID == SELHOUR )
	    {
	      cdoWarning("Hour %d not found!", intarr[isel]);
	    }
	  else if ( operatorID == SELDAY )
	    {
	      cdoWarning("Day %d not found!", intarr[isel]);
	    }
	  else if ( operatorID == SELMON )
	    {
	      cdoWarning("Month %d not found!", intarr[isel]);
	    }
	  else if ( operatorID == SELYEAR )
	    {
	      cdoWarning("Year %d not found!", intarr[isel]);
	    }
	  else if ( operatorID == SELSEAS )
	    {
	      if ( isel < 3 )
		cdoWarning("Month %d not found!", intarr[isel]);
	    }
	}
    }

  if ( selfound ) free(selfound);

  listDelete(ilist);

  if ( lnts1 || nconst )
    {
      for ( tsID = 0; tsID < nts1; tsID++ ) field_free(vars[tsID], vlistID2);

      if ( vars  ) free(vars);
      if ( vdate_list ) free(vdate_list);
      if ( vtime_list ) free(vtime_list);
    }

  vlistDestroy(vlistID2);

  if ( tsID2 == 0 ) cdoAbort("No timesteps selected!");

  cdoFinish();

  return (NULL);
}
