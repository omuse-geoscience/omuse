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

      Showinfo   showparam       Show parameters
      Showinfo   showcode        Show code numbers
      Showinfo   showname        Show variable names
      Showinfo   showstdname     Show variable standard names
      Showinfo   showlevel       Show levels
      Showinfo   showyear        Show years
      Showinfo   showmon         Show months
      Showinfo   showdate        Show dates
      Showinfo   showtime        Show timesteps
      Showinfo   showltype       Show level types
      Showinfo   showformat      Show file format
*/


#include <stdio.h>
#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Showinfo(void *argument)
{
  int zaxisID;
  int vdate, vtime;
  int nrecs;
  int nlevs, levelID;
  int ltype;
  int date0 = 0;
  int year, month, day;
  int month0 = 0, year0 = 0;

  cdoInitialize(argument);

  int SHOWYEAR      = cdoOperatorAdd("showyear",      0, 0, NULL);
  int SHOWMON       = cdoOperatorAdd("showmon",       0, 0, NULL);
  int SHOWDATE      = cdoOperatorAdd("showdate",      0, 0, NULL);
  int SHOWTIME      = cdoOperatorAdd("showtime",      0, 0, NULL);
  int SHOWTIMESTAMP = cdoOperatorAdd("showtimestamp", 0, 0, NULL);
  int SHOWCODE      = cdoOperatorAdd("showcode",      0, 0, NULL);
  int SHOWUNIT      = cdoOperatorAdd("showunit",      0, 0, NULL);
  int SHOWPARAM     = cdoOperatorAdd("showparam",     0, 0, NULL);
  int SHOWNAME      = cdoOperatorAdd("showname",      0, 0, NULL);
  int SHOWSTDNAME   = cdoOperatorAdd("showstdname",   0, 0, NULL);
  int SHOWLEVEL     = cdoOperatorAdd("showlevel",     0, 0, NULL);
  int SHOWLTYPE     = cdoOperatorAdd("showltype",     0, 0, NULL);
  int SHOWFORMAT    = cdoOperatorAdd("showformat",    0, 0, NULL);

  int operatorID = cdoOperatorID();

  int streamID = streamOpenRead(cdoStreamName(0));

  int vlistID = streamInqVlist(streamID);

  int nvars   = vlistNvars(vlistID);
  int taxisID = vlistInqTaxis(vlistID);
  int ntsteps = vlistNtsteps(vlistID);

  if ( operatorID == SHOWYEAR )
    {
      // int nyear = 0;
      int tsID = 0;
      if ( ntsteps != 0 )
	while ( (nrecs = streamInqTimestep(streamID, tsID)) )
	  {
	    vdate = taxisInqVdate(taxisID);

	    cdiDecodeDate(vdate, &year, &month, &day);
	 
	    if ( tsID == 0 || year0 != year )
	      {
		// if ( nyear == 10 ) { nyear = 0; fprintf(stdout, "\n"); }
		year0 = year;
		fprintf(stdout, " %4d", year0);
		// nyear++;
	      }

	    tsID++;
	  }
      fprintf(stdout, "\n");
    }
  else if ( operatorID == SHOWMON )
    {
      // int nmonth = 0;
      int tsID = 0;
      if ( ntsteps != 0 )
	while ( (nrecs = streamInqTimestep(streamID, tsID)) )
	  {
	    vdate = taxisInqVdate(taxisID);

	    cdiDecodeDate(vdate, &year, &month, &day);
	 
	    if ( tsID == 0 || month0 != month )
	      {
		// if ( nmonth == 12 ) { nmonth = 0; fprintf(stdout, "\n"); }
		month0 = month;
		fprintf(stdout, " %2d", month0);
		// nmonth++;
	      }

	    tsID++;
	  }
      fprintf(stdout, "\n");
    }
  else if ( operatorID == SHOWDATE )
    {
      char vdatestr[32];
      // int ndate = 0;
      int tsID  = 0;
      if ( ntsteps != 0 )
	while ( (nrecs = streamInqTimestep(streamID, tsID)) )
	  {
	    vdate = taxisInqVdate(taxisID);
	 
	    date2str(vdate, vdatestr, sizeof(vdatestr));

	    if ( tsID == 0 || date0 != vdate )
	      {
		// if ( ndate == 10 ) { ndate = 0; fprintf(stdout, "\n"); }
		date0 = vdate;
		fprintf(stdout, " %s", vdatestr);
		// ndate++;
	      }

	    tsID++;
	  }
      fprintf(stdout, "\n");
    }
  else if ( operatorID == SHOWTIME )
    {
      char vtimestr[32];
      // int nout = 0;
      int tsID = 0;
      if ( ntsteps != 0 )
	while ( (nrecs = streamInqTimestep(streamID, tsID)) )
	  {
	    // if ( nout == 4 ) { nout = 0; fprintf(stdout, "\n"); }
	    vtime = taxisInqVtime(taxisID);

	    time2str(vtime, vtimestr, sizeof(vtimestr));
	    fprintf(stdout, " %s", vtimestr);

	    tsID++;
	    // nout++;
	  }
      fprintf(stdout, "\n");
    }
  else if ( operatorID == SHOWTIMESTAMP )
    {
      char vdatetimestr[64];
      // int nout = 0;
      int tsID = 0;
      if ( ntsteps != 0 )
	while ( (nrecs = streamInqTimestep(streamID, tsID)) )
	  {
	    // if ( nout == 4 ) { nout = 0; fprintf(stdout, "\n"); }
	    vdate = taxisInqVdate(taxisID);
	    vtime = taxisInqVtime(taxisID);

	    datetime2str(vdate, vtime, vdatetimestr, sizeof(vdatetimestr));
	    fprintf(stdout, " %s", vdatetimestr);

	    tsID++;
            // nout++;
	  }
      fprintf(stdout, "\n");
    }
  else if ( operatorID == SHOWCODE )
    {
      // int nout = 0;
      for ( int varID = 0; varID < nvars; varID++ )
	{
	  // if ( nout == 20 ) { nout = 0; fprintf(stdout, "\n"); }
	  fprintf(stdout, " %d", vlistInqVarCode(vlistID, varID));
	  // nout++;
	}
      fprintf(stdout, "\n");
    }
  else if ( operatorID == SHOWUNIT )
    {
      char varunits[CDI_MAX_NAME];
      //int nout = 0;
      for ( int varID = 0; varID < nvars; varID++ )
	{
	  varunits[0] = 0;
	  vlistInqVarUnits(vlistID, varID, varunits);
	  // if ( nout == 20 ) { nout = 0; fprintf(stdout, "\n"); }
          if ( strlen(varunits) ) fprintf(stdout, " %s", varunits);
	  // nout++;
	}
      fprintf(stdout, "\n");
    }
  else if ( operatorID == SHOWPARAM )
    {
      int param;
      char paramstr[32];
      
      for ( int varID = 0; varID < nvars; varID++ )
	{
	  param = vlistInqVarParam(vlistID, varID);
	  cdiParamToString(param, paramstr, sizeof(paramstr));

	  fprintf(stdout, " %s", paramstr);
	}
      fprintf(stdout, "\n");
    }
  else if ( operatorID == SHOWNAME )
    {
      char varname[CDI_MAX_NAME];
      for ( int varID = 0; varID < nvars; varID++ )
	{
	  vlistInqVarName(vlistID, varID, varname);
	  fprintf(stdout, " %s", varname);
	}
      fprintf(stdout, "\n");
    }
  else if ( operatorID == SHOWSTDNAME )
    {
      char stdname[CDI_MAX_NAME];
      for ( int varID = 0; varID < nvars; varID++ )
	{
	  vlistInqVarStdname(vlistID, varID, stdname);
	  if ( stdname[0] != 0 )
	    fprintf(stdout, " %s", stdname);
	  else
	    fprintf(stdout, " unknown");
	}
      fprintf(stdout, "\n");
    }
  else if ( operatorID == SHOWLEVEL )
    {
      for ( int varID = 0; varID < nvars; varID++ )
	{
	  zaxisID = vlistInqVarZaxis(vlistID, varID);
	  nlevs = zaxisInqSize(zaxisID);
	  for ( levelID = 0; levelID < nlevs; levelID++ )
	    fprintf(stdout, " %.9g", zaxisInqLevel(zaxisID, levelID));
	  fprintf(stdout, "\n");
	}
    }
  else if ( operatorID == SHOWLTYPE )
    {
      int nzaxis = vlistNzaxis(vlistID);
      for ( int index = 0; index < nzaxis; index++ )
	{
	  zaxisID = vlistZaxis(vlistID, index);

	  ltype = zaxis2ltype(zaxisID);

	  if ( ltype != -1 ) fprintf(stdout, " %d", ltype);
	}
      fprintf(stdout, "\n"); 
    }
  else if ( operatorID == SHOWFORMAT )
    {
      printFiletype(streamID, vlistID);
    }

  streamClose(streamID);

  cdoFinish();

  return (0);
}
