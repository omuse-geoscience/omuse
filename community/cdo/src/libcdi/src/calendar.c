#include <limits.h>
#include <stdio.h>

#include "cdi.h"  		/* CALENDAR_ */
#include "error.h"
#include "timebase.h"


static int month_360[12] = {30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30};
static int month_365[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
static int month_366[12] = {31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};


int calendar_dpy(int calendar)
{
  int dpy = 0;

  if      ( calendar == CALENDAR_360DAYS ) dpy = 360;
  else if ( calendar == CALENDAR_365DAYS ) dpy = 365;
  else if ( calendar == CALENDAR_366DAYS ) dpy = 366;

  return (dpy);
}


int days_per_month(int calendar, int year, int month)
{
  int dayspermonth = 0;
  int *dpm = NULL;
  int dpy;

  dpy = calendar_dpy(calendar);

  if      ( dpy == 360 ) dpm = month_360;
  else if ( dpy == 365 ) dpm = month_365;
  else                   dpm = month_366;

  if ( month >= 1 && month <= 12 )
    dayspermonth = dpm[month-1];
  /*
  else
    fprintf(stderr, "days_per_month: month %d out of range\n", month);
  */
  if ( dpy == 0 && month == 2 )
    {
      if ( (year%4 == 0 && year%100 != 0) || year%400 == 0 )
	dayspermonth = 29;
      else
	dayspermonth = 28;
    }

  return (dayspermonth);
}


int days_per_year(int calendar, int year)
{
  int daysperyear;
  int dpy;

  dpy = calendar_dpy(calendar);

  if ( dpy == 0 )
    {
      if ( calendar == CALENDAR_STANDARD )
	{
	  if ( year == 1582 )
	    dpy = 355;
	  else if ( (year%4 == 0 && year%100 != 0) || year%400 == 0 )
	    dpy = 366;
	  else
	    dpy = 365;
	}
      else
	{
	  if ( (year%4 == 0 && year%100 != 0) || year%400 == 0 )
	    dpy = 366;
	  else
	    dpy = 365;
	}
    }

  daysperyear = dpy;
  
  return (daysperyear);
}


static void decode_day(int dpy, int days, int *year, int *month, int *day)
{
  int i = 0;
  int *dpm = NULL;

  *year = (days-1) / dpy;
  days -= (*year*dpy);

  if      ( dpy == 360 ) dpm = month_360;
  else if ( dpy == 365 ) dpm = month_365;
  else if ( dpy == 366 ) dpm = month_366;

  if ( dpm )
    for ( i = 0; i < 12; i++ )
      {
	if ( days > dpm[i] ) days -= dpm[i];
	else break;
      }

  *month = i + 1;
  *day   = days;
}


static int encode_day(int dpy, int year, int month, int day)
{
  int i;
  int *dpm = NULL;
  long rval = (long)dpy * year + day;

  if      ( dpy == 360 ) dpm = month_360;
  else if ( dpy == 365 ) dpm = month_365;
  else if ( dpy == 366 ) dpm = month_366;

  if ( dpm ) for ( i = 0; i < month-1; i++ ) rval += dpm[i];
  if (rval > INT_MAX || rval < INT_MIN)
    Error("Unhandled date: %ld", rval);

  return (int)rval;
}


int date_to_calday(int calendar, int date)
{
  int calday;
  int dpy;
  int year, month, day;

  dpy = calendar_dpy(calendar);

  cdiDecodeDate(date, &year, &month, &day);

  if ( dpy == 360 || dpy == 365 || dpy == 366 )
    calday = encode_day(dpy, year, month, day);
  else
    calday = encode_julday(calendar, year, month, day);

  return (calday);
}


int calday_to_date(int calendar, int calday)
{
  int date;
  int dpy;
  int year, month, day;

  dpy = calendar_dpy(calendar);

  if ( dpy == 360 || dpy == 365 || dpy == 366 )
    decode_day(dpy, calday, &year, &month, &day);
  else
    decode_julday(calendar, calday, &year, &month, &day);

  date = cdiEncodeDate(year, month, day);

  return (date);
}


void encode_caldaysec(int calendar, int year, int month, int day, int hour, int minute, int second,
		      int *julday, int *secofday)
{
  int dpy;

  dpy = calendar_dpy(calendar);

  if ( dpy == 360 || dpy == 365 || dpy == 366 )
    *julday = encode_day(dpy, year, month, day);
  else
    *julday = encode_julday(calendar, year, month, day);

  *secofday = hour*3600 + minute*60 + second;
}


void decode_caldaysec(int calendar, int julday, int secofday, 
		      int *year, int *month, int *day, int *hour, int *minute, int *second)
{
  int dpy;

  dpy = calendar_dpy(calendar);

  if ( dpy == 360 || dpy == 365 || dpy == 366 )
    decode_day(dpy, julday, year, month, day);
  else
    decode_julday(calendar, julday, year, month, day);

  *hour   = secofday/3600;
  *minute = secofday/60 - *hour*60;
  *second = secofday - *hour*3600 - *minute*60;
}


#ifdef TEST
int main(void)
{
  int calendar = CALENDAR_STANDARD;
  int nmin;
  int vdate0, vtime0;
  int vdate, vtime;
  int ijulinc;
  int i, j = 0;
  int year, mon, day, hour, minute, second;
  int calday, secofday;

  /* 1 - Check valid range of years */

  nmin = 11000;
  vdate0 = -80001201;
  vtime0 = 120500;

  printf("start time: %8d %4d\n", vdate0, vtime0);

  for ( i = 0; i < nmin; i++ )
    {
      cdiDecodeDate(vdate0, &year, &mon, &day);
      cdiDecodeTime(vtime0, &hour, &minute, &second);

      calday  = date_to_calday(calendar, vdate0);
      secofday = time_to_sec(vtime0);

      vdate = calday_to_date(calendar, calday);
      vtime = sec_to_time(secofday);

      if ( vdate0 != vdate || vtime0 != vtime )
	printf("%4d %8d %4d %8d %4d %9d %9d\n",
	       ++j, vdate0, vtime0, vdate, vtime, calday, secofday);

      year++;
      vdate0 = cdiEncodeDate(year, mon, day);
      vtime0 = cdiEncodeTime(hour, minute, second);
    }

  printf("stop time: %8d %4d\n", vdate0, vtime0);

  /* 2 - Check time increment of one minute */

  nmin = 120000;
  ijulinc = 60;
  vdate0 = 20001201;
  vtime0 = 0;

  printf("start time: %8d %4d\n", vdate0, vtime0);

  calday = date_to_calday(calendar, vdate0);
  secofday = time_to_sec(vtime0);
  for ( i = 0; i < nmin; i++ )
    {
      cdiDecodeDate(vdate0, &year, &mon, &day);
      cdiDecodeTime(vtime0, &hour, &minute, &second);

      if ( ++minute >= 60 )
	{
	  minute = 0;
	  if ( ++hour >= 24 )
	    {
	      hour = 0;
	      if ( ++day >= 32 )
		{
		  day = 1;
		  if ( ++mon >= 13 )
		    {
		      mon = 1;
		      year++;
		    }
		}
	    }
	}

      vdate0 = cdiEncodeDate(year, mon, day);
      vtime0 = cdiEncodeTime(hour, minute, second);

      julday_add_seconds(ijulinc, &calday, &secofday);

      vdate = calday_to_date(calendar, calday);
      vtime = sec_to_time(secofday);
      if ( vdate0 != vdate || vtime0 != vtime )
	printf("%4d %8d %4d %8d %4d %9d %9d\n",
	       ++j, vdate0, vtime0, vdate, vtime, calday, secofday);
    }

  printf("stop time: %8d %4d\n", vdate0, vtime0);

  return (0);
}
#endif


#ifdef TEST2
int main(void)
{
  int calendar = CALENDAR_STANDARD;
  int i;
  int calday, secofday;
  int year, month, day, hour, minute, second;
  int value = 30;
  int factor = 86400;

  calendar = CALENDAR_360DAYS;

  year=1979; month=1; day=15; hour=12; minute=30; second = 0;

  printf("calendar = %d\n", calendar);
  printf("%d/%02d/%02d %02d:%02d:%02d\n", year, month, day, hour, minute, second);

  encode_caldaysec(calendar, year, month, day, hour, minute, second, &calday, &secofday);

  decode_caldaysec(calendar, calday, secofday, &year, &month, &day, &hour, &minute, &second);
  printf("%d/%02d/%02d %02d:%02d:%02d   %d %d\n", year, month, day, hour, minute, second, calday, secofday);

  for ( i = 0; i < 420; i++ )
    {

      decode_caldaysec(calendar, calday, secofday, &year, &month, &day, &hour, &minute, &second);
      printf("%2d %d/%02d/%02d %02d:%02d:%02d\n", i, year, month, day, hour, minute, second);
      julday_add_seconds(value*factor, &calday, &secofday);
    }

  return (0);
}
#endif
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
