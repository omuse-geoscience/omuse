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

      Tinfo      tinfo           Time information
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


#define MAX_GAPS   64
#define MAX_NTSM  128
#define LIM_NTSM 1024

enum {TU_SECONDS=0, TU_MINUTES, TU_HOURS, TU_DAYS, TU_MONTHS, TU_YEARS};
char *tunits[] = {"second", "minute", "hour", "day", "month", "year"};
int   iunits[] = {1, 60, 3600, 86400, 1, 12};


void getTimeInc(double jdelta, int vdate0, int vdate1, int *incperiod, int *incunit)
{
  int year0, month0, day0;
  int year1, month1, day1;
  int deltam, deltay;
  int lperiod;
  int sign = 1;

  *incperiod = 0;
  *incunit   = 0;

  if ( jdelta < 0 )
    lperiod = (int)(jdelta-0.5);
  else
    lperiod = (int)(jdelta+0.5);

  if ( lperiod < 0 )
    {
      int tmp;
      tmp = vdate1;
      vdate1 = vdate0;
      vdate0 = tmp;
      lperiod = -lperiod;
      sign = -1;
    }

  // printf("\n%d %d %d\n",lperiod, vdate0, vdate1);

  cdiDecodeDate(vdate0, &year0, &month0, &day0);
  cdiDecodeDate(vdate1, &year1, &month1, &day1);

  deltay = year1-year0;
  deltam = deltay*12 + (month1-month0);

  if ( lperiod/60 > 0 && lperiod/60 < 60 )
    {
      *incperiod = lperiod/60;
      *incunit = TU_MINUTES;
    }
  else if ( lperiod/3600 > 0 && lperiod/3600 < 24 )
    {
      *incperiod = lperiod/3600;
      *incunit = TU_HOURS;
    }
  else if ( lperiod/(3600*24) > 0 && lperiod/(3600*24) < 32 )
    {
      *incperiod = lperiod/(3600*24);
      *incunit = TU_DAYS;
      if ( *incperiod > 27 && deltam == 1 )
	{
	  *incperiod = 1;
	  *incunit = TU_MONTHS;
	}
    }
  else if ( lperiod/(3600*24*30) > 0 && lperiod/(3600*24*30) < 12 )
    {
      *incperiod = deltam;
      *incunit = TU_MONTHS;
    }
  else if ( lperiod/(3600*24*30*12) > 0 )
    {
      *incperiod = deltay;
      *incunit = TU_YEARS;
    }
  else
    {
      *incperiod = lperiod;
      *incunit = TU_SECONDS;
    }

  *incperiod *= sign;
}

static
void printBounds(int taxisID, int calendar)
{
  int vdate0, vdate1;
  int vtime0, vtime1;
  int incperiod = 0, incunit = 0;
  juldate_t juldate1, juldate0;
  double jdelta;
  int i, len;
  char vdatestr[32], vtimestr[32];

  taxisInqVdateBounds(taxisID, &vdate0, &vdate1);
  taxisInqVtimeBounds(taxisID, &vtime0, &vtime1);

  date2str(vdate0, vdatestr, sizeof(vdatestr));
  time2str(vtime0, vtimestr, sizeof(vtimestr));	  
  fprintf(stdout, " %s %s", vdatestr, vtimestr);

  date2str(vdate1, vdatestr, sizeof(vdatestr));
  time2str(vtime1, vtimestr, sizeof(vtimestr));	  
  fprintf(stdout, " %s %s", vdatestr, vtimestr);

  juldate0  = juldate_encode(calendar, vdate0, vtime0);
  juldate1  = juldate_encode(calendar, vdate1, vtime1);
  jdelta    = juldate_to_seconds(juldate_sub(juldate1, juldate0));

  getTimeInc(jdelta, vdate0, vdate1, &incperiod, &incunit);
  
  /* fprintf(stdout, "  %g  %g  %g  %d", jdelta, jdelta/3600, fmod(jdelta,3600), incperiod%3600);*/
  len = fprintf(stdout, " %3d %s%s", incperiod, tunits[incunit], abs(incperiod)>1?"s":"");
  for ( i = 0; i < 11-len; ++i ) fprintf(stdout, " ");
}

static
int fill_gap(int ngaps, int ntsm[MAX_NTSM], int rangetsm[MAX_GAPS][2], 
	     int vdatem[MAX_GAPS][MAX_NTSM], int vtimem[MAX_GAPS][MAX_NTSM],
	     int tsID, int incperiod0, int incunit0, int vdate, int vdate0, int vtime0,
	     int calendar, int day0, juldate_t juldate, juldate_t juldate0)
{
  int its = 0;
  int year, month, day;
  int ndate, ntime;
  int ijulinc = incperiod0 * iunits[incunit0];

  if ( ijulinc > 0 && ngaps < MAX_GAPS )
    {
      rangetsm[ngaps][0] = tsID;
      rangetsm[ngaps][1] = tsID+1;

      if ( incunit0 == TU_MONTHS || incunit0 == TU_YEARS )
	{
	  its = 0;
	  ndate = vdate0;
	  //printf("fill_gap %d\n", ndate);
	  while ( TRUE )
	    {
	      cdiDecodeDate(ndate, &year, &month, &day);
	      
	      month += ijulinc;
				  
	      while ( month > 12 ) { month -= 12; year++; }
	      while ( month <  1 ) { month += 12; year--; }
	      
	      if ( day0 == 31 )
		day = days_per_month(calendar, year, month);
	      
	      ndate = cdiEncodeDate(year, month, day);
	      ntime = vtime0;
	      if ( ndate >= vdate ) break;
	      /* printf("\n1 %d %d\n", ndate, ntime); */
	      if ( its < MAX_NTSM )
		{
		  vdatem[ngaps][its] = ndate;
		  vtimem[ngaps][its] = ntime;
		}
	      else if ( its >= LIM_NTSM )
		break;

	      its++;
	    }
	}
      else
	{
	  its = 0;
	  juldate0 = juldate_add_seconds(ijulinc, juldate0);
	  while ( juldate_to_seconds(juldate0) < juldate_to_seconds(juldate) )
	    {
	      juldate_decode(calendar, juldate0, &ndate, &ntime);
	      juldate0 = juldate_add_seconds(ijulinc, juldate0);
	      if ( its < MAX_NTSM )
		{
		  vdatem[ngaps][its] = ndate;
		  vtimem[ngaps][its] = ntime;
		}
	      else if ( its >= LIM_NTSM )
		break;

	      its++;
	    }
	}			
      ntsm[ngaps] = its;
    }

  return (its);
}


void *Tinfo(void *argument)
{
  int vdate_first = 0, vtime_first = 0;
  int vdate0 = 0, vtime0 = 0;
  int vdate = 0, vtime = 0;
  int fdate = 0, ftime = 0;
  int nrecs, ntsteps;
  int tsID = 0, ntimeout;
  int taxisID;
  int streamID;
  int vlistID;
  int year0, month0, day0;
  int year, month, day;
  int calendar, unit;
  int lforecast = FALSE;
  int incperiod0 = 0, incunit0 = 0;
  int incperiod = 0, incunit = 0;
  int its = 0, igap;
  int ngaps = 0;
  int ntsm[MAX_GAPS];
  int rangetsm[MAX_GAPS][2];
  int vdatem[MAX_GAPS][MAX_NTSM];
  int vtimem[MAX_GAPS][MAX_NTSM];
  juldate_t juldate, juldate0;
  double jdelta = 0, jdelta0 = 0;
  int arrow = 0;
  int i, len;
  char vdatestr[32], vtimestr[32];	  

  cdoInitialize(argument);

  streamID = streamOpenRead(cdoStreamName(0));

  vlistID = streamInqVlist(streamID);

  fprintf(stdout, "\n");

  taxisID = vlistInqTaxis(vlistID);
  ntsteps = vlistNtsteps(vlistID);

  if ( ntsteps != 0 )
    {
      if ( ntsteps == CDI_UNDEFID )
	fprintf(stdout, "   Time axis :  unlimited steps\n");
      else
	fprintf(stdout, "   Time axis :  %d step%s\n", ntsteps, ntsteps == 1 ? "" : "s");

      if ( taxisID != CDI_UNDEFID )
	{
	  if ( taxisInqType(taxisID) != TAXIS_ABSOLUTE )
	    {
	      vdate = taxisInqRdate(taxisID);
	      vtime = taxisInqRtime(taxisID);
	      
	      date2str(vdate, vdatestr, sizeof(vdatestr));
	      time2str(vtime, vtimestr, sizeof(vtimestr));

	      fprintf(stdout, "     RefTime = %s %s", vdatestr, vtimestr);
		      
	      unit = taxisInqTunit(taxisID);
	      if ( unit != CDI_UNDEFID )  fprintf(stdout, "  Units = %s", tunit2str(unit));
	      
	      calendar = taxisInqCalendar(taxisID);
	      if ( calendar != CDI_UNDEFID )  fprintf(stdout, "  Calendar = %s", calendar2str(calendar));

	      if ( taxisHasBounds(taxisID) )
		fprintf(stdout, "  Bounds = true");

	      fprintf(stdout, "\n");

	      if ( taxisInqType(taxisID) == TAXIS_FORECAST )
		{
		  fdate = taxisInqFdate(taxisID);
		  ftime = taxisInqFtime(taxisID);
	      
		  date2str(fdate, vdatestr, sizeof(vdatestr));
		  time2str(ftime, vtimestr, sizeof(vtimestr));

		  fprintf(stdout, "     Forecast RefTime = %s %s", vdatestr, vtimestr);
		      
		  unit = taxisInqForecastTunit(taxisID);
		  if ( unit != CDI_UNDEFID )  fprintf(stdout, "  Units = %s", tunit2str(unit));

		  fprintf(stdout, "\n");

		  lforecast = TRUE;
		}
	    }
	}

      calendar = taxisInqCalendar(taxisID);

      fprintf(stdout, "\n");
      fprintf(stdout, "         Verification Time              ");
      if ( lforecast ) fprintf(stdout, " Forecast Reference Time     ");
      if ( taxisHasBounds(taxisID) )
	fprintf(stdout, " lower bound          upper bound");
      fprintf(stdout, "\n");

      fprintf(stdout, "Timestep YYYY-MM-DD hh:mm:ss   Increment");
      if ( lforecast ) fprintf(stdout, " YYYY-MM-DD hh:mm:ss   Period");
      if ( taxisHasBounds(taxisID) )
	fprintf(stdout, " YYYY-MM-DD hh:mm:ss  YYYY-MM-DD hh:mm:ss  Difference");
      fprintf(stdout, "\n");

      tsID = 0;
      while ( (nrecs = streamInqTimestep(streamID, tsID)) )
	{  
	  vdate = taxisInqVdate(taxisID);
	  vtime = taxisInqVtime(taxisID);
	  
	  cdiDecodeDate(vdate, &year, &month, &day);
	      
	  date2str(vdate, vdatestr, sizeof(vdatestr));
	  time2str(vtime, vtimestr, sizeof(vtimestr));

	  fprintf(stdout, "%6d  %s %s", tsID+1, vdatestr, vtimestr);

	  if ( tsID )
	    {
	      cdiDecodeDate(vdate0, &year0, &month0, &day0);

	      juldate0  = juldate_encode(calendar, vdate0, vtime0);
	      juldate   = juldate_encode(calendar, vdate, vtime);
	      jdelta    = juldate_to_seconds(juldate_sub(juldate, juldate0));

	      getTimeInc(jdelta, vdate0, vdate, &incperiod, &incunit);

	      /* fprintf(stdout, "  %g  %g  %g  %d", jdelta, jdelta/3600, fmod(jdelta,3600), incperiod%3600);*/
	      len = fprintf(stdout, " %3d %s%s", incperiod, tunits[incunit], abs(incperiod)>1?"s":"");
	      for ( i = 0; i < 11-len; ++i ) fprintf(stdout, " ");
	    }
	  else
	    {
	      vdate_first = vdate;
	      vtime_first = vtime;
	      fprintf(stdout, "   --------");
	    }

	  if ( lforecast )
	    {
	      fdate = taxisInqFdate(taxisID);
	      ftime = taxisInqFtime(taxisID);
	      
	      date2str(fdate, vdatestr, sizeof(vdatestr));
	      time2str(ftime, vtimestr, sizeof(vtimestr));

	      fprintf(stdout, " %s %s", vdatestr, vtimestr);

	      double fc_period = taxisInqForecastPeriod(taxisID);
	      fprintf(stdout, " %7g", fc_period);
	    }

	  if ( taxisHasBounds(taxisID) ) printBounds(taxisID, calendar);

	  if (  tsID > 1 && (incperiod != incperiod0 || incunit != incunit0) )
	    {
	      if ( tsID == 2 && (jdelta0 > jdelta) )
		{
		  jdelta0    = jdelta;
		  incperiod0 = incperiod;
		  incunit0   = incunit;

		  its = fill_gap(ngaps, ntsm, rangetsm, vdatem, vtimem,
				 1, incperiod0, incunit0, vdate_first, vdate, vtime,
				 calendar, day, juldate0, 
				 juldate_encode(calendar, vdate_first, vtime_first));

		  arrow = '^';
		}
	      else
		{
		  its = fill_gap(ngaps, ntsm, rangetsm, vdatem, vtimem,
				 tsID, incperiod0, incunit0, vdate, vdate0, vtime0,
				 calendar, day0, juldate, juldate0);

		  arrow = '<';

		  if (  its == 0 && incperiod < 0 )
		    {
		      its = -1;
		      vdate = vdate0;
		      vtime = vtime0;
		    }
		}

	      if ( its > 0 )
		{
		  ngaps++;
		  if ( cdoVerbose )
		    fprintf(stdout, "  %c--- Gap %d, missing %s%d timestep%s",
			    arrow, ngaps, its>=LIM_NTSM?"more than ":"", its, its>1?"s":"");
		}
	      else if ( its < 0 )
		{
		  if ( cdoVerbose )
		    fprintf(stdout, "  %c--- Wrong date/time information, negative increment!", arrow);
		}
	    }

	  if ( tsID == 1 )
	    {
	      jdelta0    = jdelta;
	      incperiod0 = incperiod;
	      incunit0   = incunit;
	    }

	  fprintf(stdout, "\n");

	  vdate0 = vdate;
	  vtime0 = vtime;

	  tsID++;
	}
    }

  streamClose(streamID);

  fprintf(stdout, "\n");

  date2str(vdate_first, vdatestr, sizeof(vdatestr));
  time2str(vtime_first, vtimestr, sizeof(vtimestr));
  fprintf(stdout, " Start date          : %s %s\n", vdatestr, vtimestr);

  date2str(vdate, vdatestr, sizeof(vdatestr));
  time2str(vtime, vtimestr, sizeof(vtimestr));
  fprintf(stdout, " End date            : %s %s\n", vdatestr, vtimestr);

  fprintf(stdout, " Increment           : %3d %s%s\n", 
	  incperiod0, tunits[incunit0], incperiod0>1?"s":"");
  fprintf(stdout, " Number of timesteps : %d\n", tsID);
  fprintf(stdout, " Gaps identified     : %d\n", ngaps);

  if ( cdoVerbose && ngaps )
    {
      fprintf(stdout, "\nFound potentially %d gap%s in the time series", ngaps, ngaps>1?"s":"");
      if ( ngaps >= MAX_GAPS )
	{
	  ngaps = MAX_GAPS;
	  fprintf(stdout, ", here are the first %d", ngaps);
	}
      fprintf(stdout, ":\n");
      for ( igap = 0; igap < ngaps; ++igap )
	{
	  fprintf(stdout, "  Gap %d between timestep %d and %d, missing %d timestep%s",
		  igap+1, rangetsm[igap][0], rangetsm[igap][1], ntsm[igap], ntsm[igap]>1?"s":"");
	  if ( ntsm[igap] >= MAX_NTSM )
	    {
	      ntsm[igap] = MAX_NTSM;
	      fprintf(stdout, ", here are the first %d", ntsm[igap]);
	    }
	  fprintf(stdout, ":\n");
	  
	  ntimeout = 0;
	  for ( its = 0; its < ntsm[igap]; ++its )
	    {
	      if ( ntimeout == 4 )
		{
		  ntimeout = 0;
		  fprintf(stdout, "\n");
		}

	      vdate = vdatem[igap][its];
	      vtime = vtimem[igap][its];

	      date2str(vdate, vdatestr, sizeof(vdatestr));
	      time2str(vtime, vtimestr, sizeof(vtimestr));
	      fprintf(stdout, "  %s %s", vdatestr, vtimestr);

	      ntimeout++;
	      tsID++;
	    }
	  fprintf(stdout, "\n");
	}
    }

  cdoFinish();

  return (0);
}
