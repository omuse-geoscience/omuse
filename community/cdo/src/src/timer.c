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

#if defined(HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#if defined(HAVE_SYS_TIMES_H)
#include <sys/times.h>
#endif


void cdoProcessTime(double *utime, double *stime)
{
#if defined(HAVE_SYS_TIMES_H)
  struct tms tbuf;
  double clock_ticks = 0;

  if ( (int)times(&tbuf) == -1 )
    {
      *utime = 0;
      *stime = 0;
    }
  else
    {
#if defined(_SC_CLK_TCK)
      clock_ticks = (double) sysconf(_SC_CLK_TCK);
#endif
      if ( clock_ticks > 0 )
	{
	  *utime = ((double) tbuf.tms_utime) / clock_ticks; 
	  *stime = ((double) tbuf.tms_stime) / clock_ticks;
	  /*
	    printf("utime  %g\n", ((double) tbuf.tms_utime) / clock_ticks);
	    printf("stime  %g\n", ((double) tbuf.tms_stime) / clock_ticks);
	    printf("cutime %g\n", ((double) tbuf.tms_cutime) / clock_ticks);
	    printf("cstime %g\n", ((double) tbuf.tms_cstime) / clock_ticks);
	  */
	}
      else
	{
	  *utime = 0;
	  *stime = 0;
	}
    }
#else
  *utime = 0;
  *stime = 0;
#endif
}


void util_init_real_time(void);
void util_read_real_time(void *it);
void util_diff_real_time(void *it1, void *it2, double *t);


FILE *rt_unit = NULL;

enum {rt_stat_undef, rt_stat_on, rt_stat_off};

/* minimal internal time needed to do one measurement */

double tm_shift = 0.0;

/* average total overhead for one timer_start & timer_stop pair */

double tm_overhead = 0.0;


#define  MAX_TIMER    128  /* max number of timers allowed */


typedef struct {
  int    reserved;
  int    calls;
  int    stat;
  double tot;
  double min;
  double max;
  double last;
  char   mark1[32]; /* max: 16 on IBM; 8 for double (all other) */
  char   text[128];
}
RT_TYPE;

RT_TYPE rt[MAX_TIMER];
RT_TYPE rt_init = {0, 0, 0, 0, 1.e30, 0, 0, "", "noname"};
int timer_need_init = 1;
int top_timer = 0;

static
void set_time_mark(void *mark)
{
  util_read_real_time(mark);
}

static
double get_time_val(void *mark0)
{
  double dt;
  char mark[32];

  util_read_real_time(mark);
  util_diff_real_time(mark0, mark, &dt);

  dt -= tm_shift;

  return (dt);
}

int ntests = 100; /* tests need about n microsecs on pwr4 */

static
double m1(void)
{
  double dt, dt0;
  int i;
  char mark[32];

  dt0 = 1.0;
  for ( i = 0; i < ntests; i++ )
    {
      set_time_mark(mark);
      dt = get_time_val(mark);
      if ( dt < dt0 ) dt0 = dt;
    }

  return (dt0);
}

static
double m2(void)
{
  char mark1[32], mark2[32];
  double dt1, dt2, dt0;
  int i;

  dt0 = 1.0;
  for ( i = 0; i < ntests; i++ )
    {
      set_time_mark(mark2);
      set_time_mark(mark1);
      dt1 = get_time_val(mark1);
      dt2 = get_time_val(mark2);
      if ( dt2 < dt0 ) dt0 = dt2;
      if ( dt2 < dt1 ) fprintf(rt_unit, "estimate_overhead: internal error\n");
    }

  return (dt0);
}

static
void estimate_overhead(void)
{
  tm_shift    = m1();
  tm_overhead = m2();
}


static
void timer_init(void)
{
  rt_unit = stderr;

  util_init_real_time();

  estimate_overhead();

  timer_need_init = 0;
}


int timer_new(const char *text)
{
  int it;

  if ( timer_need_init ) timer_init();

  if ( top_timer > MAX_TIMER )
    {
      for ( it = 0; it < MAX_TIMER; it++ )
	fprintf(stderr, "timer %3d:  %s\n", it, rt[it].text);

      fprintf(stderr, "timer_new: MAX_TIMER too small!\n");
    }

  it = top_timer;
  top_timer++;

  rt[it] = rt_init;
  rt[it].reserved = 1;

  if ( text ) strcpy(rt[it].text, text);

  return (it);
}

static
void timer_check(int it)
{
  if ( it < 0 || it > (top_timer-1) )
    fprintf(rt_unit, "timer: invalid timer id %d\n", it);
}


double timer_val(int it)
{
  double val, dt;

  timer_check(it);

  val = rt[it].tot;

  if ( rt[it].stat == rt_stat_on )
    {
      dt = get_time_val(rt[it].mark1);
      val += dt;
    }

  return (val);
}

static
void timer_header(void)
{
  fprintf(rt_unit, "\nTimer report:  shift = %g\n", tm_shift);
  fprintf(rt_unit,
	  "timer  calls        t_min    t_average        t_max      t_total  text\n");
}


void timer_report(void)
{
  int it;
  double total, avg;

  timer_header();

  for ( it = 0; it < top_timer; it++ )
    {
      total = timer_val(it);

      avg = rt[it].tot;
      if ( rt[it].calls > 0 ) avg /= rt[it].calls;

      if ( rt[it].stat != rt_stat_undef )
	fprintf(rt_unit, "%4d %7d %12.4g %12.4g %12.4g %12.4g  %s\n",
		it, rt[it].calls, rt[it].min, avg, rt[it].max, total, rt[it].text);
    }

}


void timer_start(int it)
{
  timer_check(it);

  if ( rt[it].stat == rt_stat_on )
    fprintf(rt_unit, "timer_start: timer_stop call missing\n");

  set_time_mark(rt[it].mark1);

  rt[it].stat = rt_stat_on;
}


void timer_stop(int it)
{
  double dt;

  timer_check(it);

  if ( rt[it].stat != rt_stat_on )
    {
      if ( rt[it].stat == rt_stat_off )
        fprintf(rt_unit, "timer_stop: timer_start call missing\n");
      else
        fprintf(rt_unit, "timer_stop: undefined timer >%s<\n", rt[it].text);
    }

  dt = get_time_val(rt[it].mark1);

  rt[it].last  = dt;
  rt[it].tot  += dt;
  rt[it].calls++;
  if ( dt < rt[it].min ) rt[it].min = dt;
  if ( dt > rt[it].max ) rt[it].max = dt;

  rt[it].stat = rt_stat_off;
}



#include "counter.h"

static
void counter_init(counter_t *counter)
{
  counter->cputime = 0;
}


void counter_start(counter_t *counter)
{
  counter_init(counter);

  if ( timer_need_init ) timer_init();

  set_time_mark(counter->mark);
}


void counter_stop(counter_t *counter)
{
  counter->cputime = get_time_val(counter->mark);
}


double counter_cputime(counter_t counter)
{
  return (counter.cputime);
}
