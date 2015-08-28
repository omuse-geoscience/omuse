#include <string.h>
#include <math.h>

#include <cdi.h>
#include "cdo.h"
#include "datetime.h"

static int timestat_date = -1;

static
void get_timestat_date(int *tstat_date)
{
  char *envstr;

  envstr = getenv("CDO_TIMESTAT_DATE");
  if ( envstr == NULL ) envstr = getenv("RUNSTAT_DATE");
  if ( envstr )
    {
      int env_date = -1;
      char envstrl[8];

      memcpy(envstrl, envstr, 8);
      envstrl[7] = 0;
      strtolower(envstrl);

      if      ( memcmp(envstrl, "first", 5)  == 0 )  env_date = TIMESTAT_FIRST;
      else if ( memcmp(envstrl, "last", 4)   == 0 )  env_date = TIMESTAT_LAST;
      else if ( memcmp(envstrl, "middle", 6) == 0 )  env_date = TIMESTAT_MEAN;

      if ( env_date >= 0 )
	{
	  *tstat_date = env_date;

	  if ( cdoVerbose ) cdoPrint("Set CDO_TIMESTAT_DATE to %s", envstr);
	}
    }
}


void dtlist_init(dtlist_type *dtlist)
{
  dtlist->nalloc     = 0;
  dtlist->size       = 0;
  dtlist->calendar   = CALENDAR_STANDARD;
  dtlist->has_bounds = -1;
  dtlist->stat       = TIMESTAT_LAST;
  dtlist->dtinfo     = NULL;

  if ( timestat_date == -1 )
    {
      timestat_date = 0;
      get_timestat_date(&timestat_date);
    }
}


dtlist_type *dtlist_new(void)
{
  dtlist_type *dtlist = (dtlist_type *) malloc(sizeof(dtlist_type));

  dtlist_init(dtlist);

  return dtlist;
}


void dtlist_delete(dtlist_type *dtlist)
{
  if ( dtlist->nalloc > 0 && dtlist->dtinfo )
    free(dtlist->dtinfo);

  free(dtlist);
}


void dtlist_taxisInqTimestep(dtlist_type *dtlist, int taxisID, int tsID)
{
  size_t NALLOC = 128;

  if ( (size_t)tsID >= dtlist->nalloc )
    {
      dtlist->nalloc += NALLOC;
      dtlist->dtinfo = (dtinfo_type *) realloc(dtlist->dtinfo, dtlist->nalloc*sizeof(dtinfo_type));
    }

  if ( (size_t)tsID >= dtlist->size ) dtlist->size = (size_t)tsID + 1;

  dtlist->dtinfo[tsID].v.date = taxisInqVdate(taxisID);
  dtlist->dtinfo[tsID].v.time = taxisInqVtime(taxisID);

  if ( tsID == 0 && dtlist->has_bounds == -1 )
    {
      dtlist->has_bounds = 0;
      if ( taxisHasBounds(taxisID) ) dtlist->has_bounds = 1;
    }

  if ( dtlist->has_bounds )
    {
      taxisInqVdateBounds(taxisID, &(dtlist->dtinfo[tsID].b[0].date), &(dtlist->dtinfo[tsID].b[1].date));
      taxisInqVtimeBounds(taxisID, &(dtlist->dtinfo[tsID].b[0].time), &(dtlist->dtinfo[tsID].b[1].time));
    }
  else
    {
      dtlist->dtinfo[tsID].b[0].date = 0;
      dtlist->dtinfo[tsID].b[1].date = 0;
      dtlist->dtinfo[tsID].b[0].time = 0;
      dtlist->dtinfo[tsID].b[1].time = 0;
    }
}


void dtlist_taxisDefTimestep(dtlist_type *dtlist, int taxisID, int tsID)
{
  if ( tsID < 0 || (size_t)tsID >= dtlist->size )
    cdoAbort("Internal error; tsID out of bounds!");

  taxisDefVdate(taxisID, dtlist->dtinfo[tsID].v.date);
  taxisDefVtime(taxisID, dtlist->dtinfo[tsID].v.time);
  if ( dtlist->has_bounds )
    {
      taxisDefVdateBounds(taxisID, dtlist->dtinfo[tsID].b[0].date, dtlist->dtinfo[tsID].b[1].date);
      taxisDefVtimeBounds(taxisID, dtlist->dtinfo[tsID].b[0].time, dtlist->dtinfo[tsID].b[1].time);
    }
}


void dtlist_mean(dtlist_type *dtlist, int nsteps)
{
  int vdate, vtime;

  if ( nsteps%2 == 0 )
    {
      int calendar = dtlist->calendar;

      vdate = dtlist->dtinfo[nsteps/2-1].v.date;
      vtime = dtlist->dtinfo[nsteps/2-1].v.time;
      juldate_t juldate1 = juldate_encode(calendar, vdate, vtime);

      vdate = dtlist->dtinfo[nsteps/2].v.date;
      vtime = dtlist->dtinfo[nsteps/2].v.time;
      juldate_t juldate2 = juldate_encode(calendar, vdate, vtime);

      double seconds = juldate_to_seconds(juldate_sub(juldate2, juldate1)) / 2;
      juldate_t juldatem = juldate_add_seconds((int)lround(seconds), juldate1);
      juldate_decode(calendar, juldatem, &vdate, &vtime);
    }
  else
    {
      vdate = dtlist->dtinfo[nsteps/2].v.date;
      vtime = dtlist->dtinfo[nsteps/2].v.time;
    }

  dtlist->timestat.v.date = vdate;
  dtlist->timestat.v.time = vtime;
}


void dtlist_stat_taxisDefTimestep(dtlist_type *dtlist, int taxisID, int nsteps)
{
  if ( (size_t)nsteps > dtlist->size )
    cdoAbort("Internal error; unexpected nsteps=%d (limit=%ld)!", nsteps, dtlist->size);

  int stat = dtlist->stat;
  if ( timestat_date > 0 ) stat = timestat_date;

  if      ( stat == TIMESTAT_MEAN  ) dtlist_mean(dtlist, nsteps);
  else if ( stat == TIMESTAT_FIRST ) dtlist->timestat.v = dtlist->dtinfo[0].v;
  else if ( stat == TIMESTAT_LAST  ) dtlist->timestat.v = dtlist->dtinfo[nsteps-1].v;

  if ( dtlist->has_bounds )
    {
      dtlist->timestat.b[0] = dtlist->dtinfo[0].b[0];
      dtlist->timestat.b[1] = dtlist->dtinfo[nsteps-1].b[1];
    }
  else
    {
      dtlist->timestat.b[0] = dtlist->dtinfo[0].v;
      dtlist->timestat.b[1] = dtlist->dtinfo[nsteps-1].v;
    }

  taxisDefVdate(taxisID, dtlist->timestat.v.date);
  taxisDefVtime(taxisID, dtlist->timestat.v.time);
  // if ( dtlist->has_bounds )
    {
      taxisDefVdateBounds(taxisID, dtlist->timestat.b[0].date, dtlist->timestat.b[1].date);
      taxisDefVtimeBounds(taxisID, dtlist->timestat.b[0].time, dtlist->timestat.b[1].time);
    }
}


void dtlist_shift(dtlist_type *dtlist)
{
  for ( size_t inp = 0; inp < dtlist->size-1; inp++ )
    {
      dtlist->dtinfo[inp] = dtlist->dtinfo[inp+1];
    }
}


void dtlist_set_stat(dtlist_type *dtlist, int stat)
{
  dtlist->stat = stat;
}


void dtlist_set_calendar(dtlist_type *dtlist, int calendar)
{
  dtlist->calendar = calendar;
}


int dtlist_get_vdate(dtlist_type *dtlist, int tsID)
{
  if ( tsID < 0 || (size_t)tsID >= dtlist->size )
    cdoAbort("Internal error; tsID out of bounds!");

  return dtlist->dtinfo[tsID].v.date;
}


int dtlist_get_vtime(dtlist_type *dtlist, int tsID)
{
  if ( tsID < 0 || (size_t)tsID >= dtlist->size )
    cdoAbort("Internal error; tsID out of bounds!");

  return dtlist->dtinfo[tsID].v.time;
}


void datetime_avg(int calendar, int ndates, datetime_t *datetime)
{
  int vdate, vtime;
  juldate_t juldate1, juldate2, juldatem;
  double seconds;

  if ( ndates%2 == 0 )
    {
      vdate = datetime[ndates/2-1].date;
      vtime = datetime[ndates/2-1].time;
      juldate1 = juldate_encode(calendar, vdate, vtime);

      vdate = datetime[ndates/2].date;
      vtime = datetime[ndates/2].time;
      juldate2 = juldate_encode(calendar, vdate, vtime);

      seconds = juldate_to_seconds(juldate_sub(juldate2, juldate1)) / 2;
      juldatem = juldate_add_seconds((int)lround(seconds), juldate1);
      juldate_decode(calendar, juldatem, &vdate, &vtime);
    }
  else
    {
      vdate = datetime[ndates/2].date;
      vtime = datetime[ndates/2].time;
    }

  datetime[ndates].date = vdate;
  datetime[ndates].time = vtime;
}
