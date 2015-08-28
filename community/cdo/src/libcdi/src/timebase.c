#include <stdio.h>
#include <stdint.h>
#include <math.h>		/* for floor() */

#include "cdi.h"
#include "timebase.h"

/* convert Julian date into year, months, day */
void decode_julday(int calendar,
		   int julday,	/* Julian day number to convert */
		   int *year,	/* Gregorian year (out)         */
		   int *mon,	/* Gregorian month (1-12) (out) */
		   int *day)	/* Gregorian day (1-31) (out)   */
{
  int a = julday;
  double b, c;
  double d, e, f;

  b = floor((a - 1867216.25)/36524.25);
  c = a + b - floor(b/4) + 1525;

  if ( calendar == CALENDAR_STANDARD )
    if ( a < 2299161 )
      {
	c = a + 1524;
      } 

  d = floor((c - 122.1)/365.25);
  e = floor(365.25*d);
  f = floor((c - e)/30.6001);

  *day  = (int)(c - e - floor(30.6001*f));
  *mon  = (int)(f - 1 - 12*floor(f/14));
  *year = (int)(d - 4715 - floor((7 + *mon)/10));
}


/* convert year, month, day into Julian calendar day */
int encode_julday(int calendar, int year, int month, int day)
{
  int iy;
  int im;
  int ib;
  int julday;

  if ( month <= 2 )
    {
      iy = year  - 1;
      im = month + 12;
    }
  else
    {
      iy = year;
      im = month;
    }


  if ( iy < 0 )
    ib = (int)((iy+1)/400) - (int)((iy+1)/100);
  else
    ib = (int)(iy/400) - (int)(iy/100);

  if ( calendar == CALENDAR_STANDARD )
    {
      if ( year > 1582 || (year == 1582 && (month > 10 || (month == 10 && day >= 15))) )
	{
	  /*
	  ** 15th October 1582 AD or later
	  */
	}
      else
	{
	  /*
	  ** 4th October 1582 AD or earlier
	  */
	  ib = -2;
	}
    }

  julday = (int) (floor(365.25*iy) + (int)(30.6001*(im+1)) + ib + 1720996.5 + day + 0.5);

  return (julday);
}


int date_to_julday(int calendar, int date)
{
  int julday;
  int year, month, day;

  cdiDecodeDate(date, &year, &month, &day);

  julday = encode_julday(calendar, year, month, day);

  return (julday);
}


int julday_to_date(int calendar, int julday)
{
  int date;
  int year, month, day;

  decode_julday(calendar, julday, &year, &month, &day);

  date = cdiEncodeDate(year, month, day);

  return (date);
}


int time_to_sec(int time)
{
  int secofday;
  int hour, minute, second;

  cdiDecodeTime(time, &hour, &minute, &second);

  secofday = hour*3600 + minute*60 + second;

  return (secofday);
}


int sec_to_time(int secofday)
{
  int time;
  int hour, minute, second;

  hour   = secofday/3600;
  minute = secofday/60 - hour*60;
  second = secofday - hour*3600 - minute*60;

  time = cdiEncodeTime(hour, minute, second);

  return (time);
}

static
void adjust_seconds(int *julday, int64_t *secofday)
{
  int64_t secperday = 86400;

  while ( *secofday >= secperday ) 
    { 
      *secofday -= secperday; 
      (*julday)++;
    }

  while ( *secofday <  0 ) 
    { 
      *secofday += secperday;
      (*julday)--;
    }
}


void julday_add_seconds(int64_t seconds, int *julday, int *secofday)
{
  int64_t sec_of_day = *secofday;

  sec_of_day += seconds;

  adjust_seconds(julday, &sec_of_day);

  *secofday = (int) sec_of_day;
}

/* add days and secs to julday/secofday */
void julday_add(int days, int secs, int *julday, int *secofday)
{
  int64_t sec_of_day = *secofday;

  sec_of_day += secs;
  *julday    += days;

  adjust_seconds(julday, &sec_of_day);

  *secofday = (int) sec_of_day;
}

/* subtract julday1/secofday1 from julday2/secofday2 and returns the result in seconds */
double julday_sub(int julday1, int secofday1, int julday2, int secofday2, int *days, int *secs)
{
  int64_t sec_of_day;
  int64_t seconds;

  *days = julday2 - julday1;
  *secs = secofday2 - secofday1;

  sec_of_day = *secs;

  adjust_seconds(days, &sec_of_day);

  *secs = (int) sec_of_day;

  seconds = *days * 86400 + sec_of_day;

  return (double)seconds;
}


void encode_juldaysec(int calendar, int year, int month, int day, int hour, int minute, int second, int *julday, int *secofday)
{
  *julday = encode_julday(calendar, year, month, day);

  *secofday = (hour*60 + minute)*60 + second;
}


void decode_juldaysec(int calendar, int julday, int secofday, int *year, int *month, int *day, int *hour, int *minute, int *second)
{
  decode_julday(calendar, julday, year, month, day);

  *hour   = secofday/3600;
  *minute = secofday/60 - *hour*60;
  *second = secofday - *hour*3600 - *minute*60;
}


#ifdef TEST
int main(void)
{
  int nmin;
  int vdate0, vtime0;
  int vdate, vtime;
  int ijulinc;
  int i, j = 0;
  int year, mon, day, hour, minute, second;
  int julday, secofday;
  int calendar = CALENDAR_STANDARD;

  /* 1 - Check valid range of years */

  nmin = 11000;
  vdate0 = -80001201;
  vtime0 = 120500;

  printf("start time: %8d %4d\n", vdate0, vtime0);

  for ( i = 0; i < nmin; i++ )
    {
      cdiDecodeDate(vdate0, &year, &mon, &day);
      cdiDecodeTime(vtime0, &hour, &minute, &second);

      julday  = date_to_julday(calendar, vdate0);
      secofday = time_to_sec(vtime0);

      vdate = julday_to_date(calendar, julday);
      vtime = sec_to_time(secofday);

      if ( vdate0 != vdate || vtime0 != vtime )
	printf("%4d %8d %4d %8d %4d %9d %9d\n",
	       ++j, vdate0, vtime0, vdate, vtime, julday, secofday);

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

  julday = date_to_julday(calendar, vdate0);
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

      julday_add_seconds(ijulinc, &julday, &secofday);

      vdate = julday_to_date(calendar, julday);
      vtime = sec_to_time(secofday);
      if ( vdate0 != vdate || vtime0 != vtime )
	printf("%4d %8d %4d %8d %4d %9d %9d\n",
	       ++j, vdate0, vtime0, vdate, vtime, julday, secofday);
    }

  printf("stop time: %8d %4d\n", vdate0, vtime0);

  return (0);
}
#endif


#ifdef TEST2
int main(void)
{
  int i;
  int julday, secofday;
  int year, month, day, hour, minute, second;
  int value = 30;
  int factor = 86400;
  int calendar = CALENDAR_STANDARD;

  year=1979; month=1; day=15; hour=12; minute=30, second=17;

  printf("%d/%02d/%02d %02d:%02d:%02d\n", year, month, day, hour, minute, second);

  encode_juldaysec(calendar, year, month, day, hour, minute, second, &julday, &secofday);

  decode_juldaysec(calendar, julday, secofday, &year, &month, &day, &hour, &minute, &second);
  printf("%d/%02d/%02d %02d:%02d:%02d   %d %d\n", year, month, day, hour, minute, second, julday, secofday);

  for ( i = 0; i < 420; i++ )
    {

      decode_juldaysec(calendar, julday, secofday, &year, &month, &day, &hour, &minute, &second);
      printf("%2d %d/%02d/%02d %02d:%02d:%02d\n", i, year, month, day, hour, minute, second);
      julday_add_seconds(value*factor, &julday, &secofday);
    }

  return (0);
}
#endif
