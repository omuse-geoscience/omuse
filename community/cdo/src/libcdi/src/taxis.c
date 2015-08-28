#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <stddef.h>
#include <string.h>

#include "dmemory.h"

#include "cdi.h"
#include "error.h"
#include "taxis.h"
#include "cdi_cksum.h"
#include "cdi_int.h"
#include "calendar.h"
#include "namespace.h"
#include "serialize.h"
#include "resource_handle.h"
#include "resource_unpack.h"

extern int cdiDefaultCalendar;

static int DefaultTimeType = TAXIS_ABSOLUTE;
static int DefaultTimeUnit = TUNIT_HOUR;


const char *Timeunits[] = {
  "undefined",
  "seconds",
  "minutes",
  "quarters",
  "30minutes",
  "hours",
  "3hours",
  "6hours",
  "12hours",
  "days",
  "months",
  "years",
};


static int    taxisCompareP    ( void * taxisptr1, void * taxisptr2 );
static void   taxisDestroyP    ( void * taxisptr );
static void   taxisPrintKernel(taxis_t *taxisptr, FILE * fp);
static int    taxisGetPackSize ( void * taxisptr, void *context );
static void   taxisPack        ( void * taxisptr, void *buf, int size,
				 int *position, void *context );
static int    taxisTxCode      ( void );

const resOps taxisOps = {
  taxisCompareP,
  taxisDestroyP,
  (void (*)(void *, FILE *))taxisPrintKernel,
  taxisGetPackSize,
  taxisPack,
  taxisTxCode
};

#define container_of(ptr, type, member) \
  ((type *)((unsigned char *)ptr - offsetof(type,member)))

struct refcount_string
{
  int ref_count;
  char string[];
};

static char *
new_refcount_string(size_t len)
{
  struct refcount_string *container
    = xmalloc(sizeof (*container) + len + 1);
  container->ref_count = 1;
  return container->string;
}

static void
delete_refcount_string(void *p)
{
  if (p)
    {
      struct refcount_string *container
        = container_of(p, struct refcount_string, string);
      if (!--(container->ref_count))
        free(container);
    }
}

static char *
dup_refcount_string(char *p)
{
  if (p)
    {
      struct refcount_string *container
        = container_of(p, struct refcount_string, string);
      ++(container->ref_count);
    }
  return p;
}


#undef container_of

static int  TAXIS_Debug = 0;   /* If set to 1, debugging */


const char *tunitNamePtr(int unitID)
{
  const char *name;
  int size = sizeof(Timeunits)/sizeof(*Timeunits);

  if ( unitID > 0 && unitID < size )
    name = Timeunits[unitID];
  else
    name = Timeunits[0];

  return (name);
}

#if 0
static
void taxis_defaults(void)
{
  char *timeunit;

  timeunit = getenv("TIMEUNIT");
  if ( timeunit )
    {
      if ( strcmp(timeunit, "minutes") == 0 )
	DefaultTimeUnit = TUNIT_MINUTE;
      else if ( strcmp(timeunit, "hours") == 0 )
	DefaultTimeUnit = TUNIT_HOUR;
      else if ( strcmp(timeunit, "3hours") == 0 )
	DefaultTimeUnit = TUNIT_3HOURS;
      else if ( strcmp(timeunit, "6hours") == 0 )
	DefaultTimeUnit = TUNIT_6HOURS;
      else if ( strcmp(timeunit, "12hours") == 0 )
	DefaultTimeUnit = TUNIT_12HOURS;
      else if ( strcmp(timeunit, "days") == 0 )
	DefaultTimeUnit = TUNIT_DAY;
      else if ( strcmp(timeunit, "months") == 0 )
	DefaultTimeUnit = TUNIT_MONTH;
      else if ( strcmp(timeunit, "years") == 0 )
	DefaultTimeUnit = TUNIT_YEAR;
      else
	Warning("Unsupported TIMEUNIT %s!", timeunit);
    }
}
#endif

static
void taxisDefaultValue(taxis_t* taxisptr)
{
  taxisptr->self        = CDI_UNDEFID;
  taxisptr->used        = FALSE;
  taxisptr->type        = DefaultTimeType;
  taxisptr->vdate       = 0;
  taxisptr->vtime       = 0;
  taxisptr->rdate       = CDI_UNDEFID;
  taxisptr->rtime       = 0;
  taxisptr->fdate       = CDI_UNDEFID;
  taxisptr->ftime       = 0;
  taxisptr->calendar    = cdiDefaultCalendar;
  taxisptr->unit        = DefaultTimeUnit;
  taxisptr->numavg      = 0;
  taxisptr->climatology = FALSE;
  taxisptr->has_bounds  = FALSE;
  taxisptr->vdate_lb    = 0;
  taxisptr->vtime_lb    = 0;
  taxisptr->vdate_ub    = 0;
  taxisptr->vtime_ub    = 0;
  taxisptr->fc_unit     = DefaultTimeUnit;
  taxisptr->fc_period   = 0;
  taxisptr->name        = NULL;
  taxisptr->longname    = NULL;
}

static taxis_t *
taxisNewEntry(cdiResH resH)
{
  taxis_t *taxisptr = (taxis_t*) xmalloc(sizeof(taxis_t));

  taxisDefaultValue(taxisptr);
  if (resH == CDI_UNDEFID)
    taxisptr->self = reshPut(taxisptr, &taxisOps);
  else
    {
      taxisptr->self = resH;
      reshReplace(resH, taxisptr, &taxisOps);
    }

  return (taxisptr);
}

static
void taxisInit (void)
{
  static int taxisInitialized = 0;
  char *env;

  if ( taxisInitialized ) return;

  taxisInitialized = 1;

  env = getenv("TAXIS_DEBUG");
  if ( env ) TAXIS_Debug = atoi(env);
}

#if 0
static
void taxis_copy(taxis_t *taxisptr2, taxis_t *taxisptr1)
{
  int taxisID2 = taxisptr2->self;
  memcpy(taxisptr2, taxisptr1, sizeof(taxis_t));
  taxisptr2->self = taxisID2;
}
#endif

/*
@Function  taxisCreate
@Title     Create a Time axis

@Prototype int taxisCreate(int taxistype)
@Parameter
    @Item  taxistype  The type of the Time axis, one of the set of predefined CDI time axis types.
                      The valid CDI time axis types are @func{TAXIS_ABSOLUTE} and @func{TAXIS_RELATIVE}.

@Description
The function @func{taxisCreate} creates a Time axis.

@Result
@func{taxisCreate} returns an identifier to the Time axis.

@Example
Here is an example using @func{taxisCreate} to create a relative T-axis
with a standard calendar.

@Source
#include "cdi.h"
   ...
int taxisID;
   ...
taxisID = taxisCreate(TAXIS_RELATIVE);
taxisDefCalendar(taxisID, CALENDAR_STANDARD);
taxisDefRdate(taxisID, 19850101);
taxisDefRtime(taxisID, 120000);
   ...
@EndSource
@EndFunction
*/
int taxisCreate(int taxistype)
{
  if ( CDI_Debug )
    Message("taxistype: %d", taxistype);

  taxisInit ();

  taxis_t *taxisptr = taxisNewEntry(CDI_UNDEFID);
  taxisptr->type = taxistype;

  int taxisID = taxisptr->self;

  if ( CDI_Debug )
    Message("taxisID: %d", taxisID);

  return (taxisID);
}

void taxisDestroyKernel(taxis_t *taxisptr)
{
  int id = taxisptr->self;
  delete_refcount_string(taxisptr->name);
  delete_refcount_string(taxisptr->longname);
  if (id != CDI_UNDEFID)
    reshRemove(id, &taxisOps);
}

/*
@Function  taxisDestroy
@Title     Destroy a Time axis

@Prototype void taxisDestroy(int taxisID)
@Parameter
    @Item  taxisID  Time axis ID, from a previous call to @func{taxisCreate}

@EndFunction
*/
void taxisDestroy(int taxisID)
{
  taxis_t *taxisptr = (taxis_t *)reshGetVal(taxisID, &taxisOps);
  taxisDestroyKernel(taxisptr);
  free(taxisptr);
}


void taxisDestroyP( void * taxisptr )
{
  taxisDestroyKernel((taxis_t *)taxisptr);
  free(taxisptr);
}


int taxisDuplicate(int taxisID1)
{
  taxis_t *taxisptr1 = (taxis_t *)reshGetVal(taxisID1, &taxisOps);
  taxis_t *taxisptr2 = taxisNewEntry(CDI_UNDEFID);

  int taxisID2 = taxisptr2->self;

  if ( CDI_Debug )
    Message("taxisID2: %d", taxisID2);

  ptaxisCopy(taxisptr2, taxisptr1);
  return (taxisID2);
}


void taxisDefType(int taxisID, int type)
{
  taxis_t *taxisptr = ( taxis_t * ) reshGetVal ( taxisID, &taxisOps );

  if (taxisptr->type != type)
    {
      taxisptr->type = type;
      reshSetStatus(taxisID, &taxisOps, RESH_DESYNC_IN_USE);
    }
}

/*
@Function  taxisDefVdate
@Title     Define the verification date

@Prototype void taxisDefVdate(int taxisID, int vdate)
@Parameter
    @Item  taxisID  Time axis ID, from a previous call to @fref{taxisCreate}
    @Item  vdate    Verification date (YYYYMMDD)

@Description
The function @func{taxisDefVdate} defines the verification date of a Time axis.

@EndFunction
*/
void taxisDefVdate(int taxisID, int vdate)
{
  taxis_t *taxisptr = (taxis_t *)reshGetVal(taxisID, &taxisOps);

  if (taxisptr->vdate != vdate)
    {
      taxisptr->vdate = vdate;
      reshSetStatus(taxisID, &taxisOps, RESH_DESYNC_IN_USE);
    }
}

/*
@Function  taxisDefVtime
@Title     Define the verification time

@Prototype void taxisDefVtime(int taxisID, int vtime)
@Parameter
    @Item  taxisID  Time axis ID, from a previous call to @fref{taxisCreate}
    @Item  vtime    Verification time (hhmmss)

@Description
The function @func{taxisDefVtime} defines the verification time of a Time axis.

@EndFunction
*/
void taxisDefVtime(int taxisID, int vtime)
{
  taxis_t *taxisptr = (taxis_t *)reshGetVal(taxisID, &taxisOps);

  if (taxisptr->vtime != vtime)
    {
      taxisptr->vtime = vtime;
      reshSetStatus(taxisID, &taxisOps, RESH_DESYNC_IN_USE);
    }
}

/*
@Function  taxisDefRdate
@Title     Define the reference date

@Prototype void taxisDefRdate(int taxisID, int rdate)
@Parameter
    @Item  taxisID  Time axis ID, from a previous call to @fref{taxisCreate}
    @Item  rdate    Reference date (YYYYMMDD)

@Description
The function @func{taxisDefRdate} defines the reference date of a Time axis.

@EndFunction
*/
void taxisDefRdate(int taxisID, int rdate)
{
  taxis_t *taxisptr = ( taxis_t * ) reshGetVal ( taxisID, &taxisOps );

  if (taxisptr->rdate != rdate)
    {
      taxisptr->rdate = rdate;
      reshSetStatus(taxisID, &taxisOps, RESH_DESYNC_IN_USE);
    }
}

/*
@Function  taxisDefRtime
@Title     Define the reference time

@Prototype void taxisDefRtime(int taxisID, int rtime)
@Parameter
    @Item  taxisID  Time axis ID, from a previous call to @fref{taxisCreate}
    @Item  rtime    Reference time (hhmmss)

@Description
The function @func{taxisDefRtime} defines the reference time of a Time axis.

@EndFunction
*/
void taxisDefRtime(int taxisID, int rtime)
{
  taxis_t *taxisptr = ( taxis_t * ) reshGetVal ( taxisID, &taxisOps );

  if (taxisptr->rtime != rtime)
    {
      taxisptr->rtime = rtime;
      reshSetStatus(taxisID, &taxisOps, RESH_DESYNC_IN_USE);
    }
}

/*
@Function  taxisDefFdate
@Title     Define the forecast reference date

@Prototype void taxisDefFdate(int taxisID, int fdate)
@Parameter
    @Item  taxisID  Time axis ID, from a previous call to @fref{taxisCreate}
    @Item  fdate    Forecast reference date (YYYYMMDD)

@Description
The function @func{taxisDefFdate} defines the forecast reference date of a Time axis.

@EndFunction
*/
void taxisDefFdate(int taxisID, int fdate)
{
  taxis_t *taxisptr = ( taxis_t * ) reshGetVal ( taxisID, &taxisOps );

  if (taxisptr->fdate != fdate)
    {
      taxisptr->fdate = fdate;
      reshSetStatus(taxisID, &taxisOps, RESH_DESYNC_IN_USE);
    }
}

/*
@Function  taxisDefFtime
@Title     Define the forecast reference time

@Prototype void taxisDefFtime(int taxisID, int ftime)
@Parameter
    @Item  taxisID  Time axis ID, from a previous call to @fref{taxisCreate}
    @Item  ftime    Forecast reference time (hhmmss)

@Description
The function @func{taxisDefFtime} defines the forecast reference time of a Time axis.

@EndFunction
*/
void taxisDefFtime(int taxisID, int ftime)
{
  taxis_t *taxisptr = ( taxis_t * ) reshGetVal ( taxisID, &taxisOps );

  if (taxisptr->ftime != ftime)
    {
      taxisptr->ftime = ftime;
      reshSetStatus(taxisID, &taxisOps, RESH_DESYNC_IN_USE);
    }
}

/*
@Function  taxisDefCalendar
@Title     Define the calendar

@Prototype void taxisDefCalendar(int taxisID, int calendar)
@Parameter
    @Item  taxisID  Time axis ID, from a previous call to @fref{taxisCreate}
    @Item  calendar The type of the calendar, one of the set of predefined CDI calendar types.
                    The valid CDI calendar types are @func{CALENDAR_STANDARD}, @func{CALENDAR_PROLEPTIC},
                    @func{CALENDAR_360DAYS}, @func{CALENDAR_365DAYS} and @func{CALENDAR_366DAYS}.

@Description
The function @func{taxisDefCalendar} defines the calendar of a Time axis.

@EndFunction
*/
void taxisDefCalendar(int taxisID, int calendar)
{
  taxis_t *taxisptr = ( taxis_t * ) reshGetVal ( taxisID, &taxisOps );

  if (taxisptr->calendar != calendar)
    {
      taxisptr->calendar = calendar;
      reshSetStatus(taxisID, &taxisOps, RESH_DESYNC_IN_USE);
    }
}


void taxisDefTunit(int taxisID, int unit)
{
  taxis_t *taxisptr = ( taxis_t * ) reshGetVal ( taxisID, &taxisOps );

  if (taxisptr->unit != unit)
    {
      taxisptr->unit = unit;
      reshSetStatus(taxisID, &taxisOps, RESH_DESYNC_IN_USE);
    }
}


void taxisDefForecastTunit(int taxisID, int unit)
{
  taxis_t *taxisptr = (taxis_t *)reshGetVal(taxisID, &taxisOps);

  if (taxisptr->fc_unit != unit)
    {
      taxisptr->fc_unit = unit;
      reshSetStatus(taxisID, &taxisOps, RESH_DESYNC_IN_USE);
    }
}


void taxisDefForecastPeriod(int taxisID, double fc_period)
{
  taxis_t *taxisptr = ( taxis_t * ) reshGetVal ( taxisID, &taxisOps );

  if ( IS_NOT_EQUAL(taxisptr->fc_period, fc_period) )
    {
      taxisptr->fc_period = fc_period;
      reshSetStatus(taxisID, &taxisOps, RESH_DESYNC_IN_USE);
    }
}


void taxisDefNumavg(int taxisID, int numavg)
{
  taxis_t *taxisptr = ( taxis_t * ) reshGetVal ( taxisID, &taxisOps );

  if (taxisptr->numavg != numavg)
    {
      taxisptr->numavg = numavg;
      reshSetStatus(taxisID, &taxisOps, RESH_DESYNC_IN_USE);
    }
}

/*
The type of the time axis, one of the set of predefined CDI time types.
The valid CDI time types are TAXIS_ABSOLUTE and TAXIS_RELATIVE.
*/
int taxisInqType(int taxisID)
{
  taxis_t *taxisptr = (taxis_t *)reshGetVal(taxisID, &taxisOps);

  return (taxisptr->type);
}


int taxisHasBounds(int taxisID)
{
  taxis_t *taxisptr = (taxis_t *)reshGetVal(taxisID, &taxisOps);

  return (taxisptr->has_bounds);
}


void taxisDeleteBounds(int taxisID)
{
  taxis_t *taxisptr = ( taxis_t * ) reshGetVal ( taxisID, &taxisOps );

  if (taxisptr->has_bounds != FALSE)
    {
      taxisptr->has_bounds = FALSE;
      reshSetStatus(taxisID, &taxisOps, RESH_DESYNC_IN_USE);
    }
}


void taxisCopyTimestep(int taxisID2, int taxisID1)
{
  taxis_t *taxisptr1 = (taxis_t *)reshGetVal(taxisID1, &taxisOps),
    *taxisptr2 = (taxis_t *)reshGetVal(taxisID2, &taxisOps);

  reshLock();

  taxisptr2->rdate = taxisptr1->rdate;
  taxisptr2->rtime = taxisptr1->rtime;

  taxisptr2->vdate = taxisptr1->vdate;
  taxisptr2->vtime = taxisptr1->vtime;

  if ( taxisptr2->has_bounds )
    {
      taxisptr2->vdate_lb = taxisptr1->vdate_lb;
      taxisptr2->vtime_lb = taxisptr1->vtime_lb;
      taxisptr2->vdate_ub = taxisptr1->vdate_ub;
      taxisptr2->vtime_ub = taxisptr1->vtime_ub;
    }

  taxisptr2->fdate = taxisptr1->fdate;
  taxisptr2->ftime = taxisptr1->ftime;

  taxisptr2->fc_unit   = taxisptr1->fc_unit;
  taxisptr2->fc_period = taxisptr1->fc_period;

  reshSetStatus(taxisID2, &taxisOps, RESH_DESYNC_IN_USE);
  reshUnlock();
}

/*
@Function  taxisInqVdate
@Title     Get the verification date

@Prototype int taxisInqVdate(int taxisID)
@Parameter
    @Item  taxisID  Time axis ID, from a previous call to @fref{taxisCreate} or @fref{vlistInqTaxis}

@Description
The function @func{taxisInqVdate} returns the verification date of a Time axis.

@Result
@func{taxisInqVdate} returns the verification date.

@EndFunction
*/
int taxisInqVdate(int taxisID)
{
  taxis_t *taxisptr = (taxis_t *)reshGetVal(taxisID, &taxisOps);

  return (taxisptr->vdate);
}


void taxisInqVdateBounds(int taxisID, int *vdate_lb, int *vdate_ub)
{
  taxis_t *taxisptr = (taxis_t *)reshGetVal(taxisID, &taxisOps);

  *vdate_lb = taxisptr->vdate_lb;
  *vdate_ub = taxisptr->vdate_ub;
}


void taxisDefVdateBounds(int taxisID, int vdate_lb, int vdate_ub)
{
  taxis_t *taxisptr = ( taxis_t * ) reshGetVal ( taxisID, &taxisOps );

  if (taxisptr->vdate_lb != vdate_lb
      || taxisptr->vdate_ub != vdate_ub
      || taxisptr->has_bounds != TRUE)
    {
      taxisptr->vdate_lb = vdate_lb;
      taxisptr->vdate_ub = vdate_ub;
      taxisptr->has_bounds = TRUE;
      reshSetStatus(taxisID, &taxisOps, RESH_DESYNC_IN_USE);
    }
}

/*
@Function  taxisInqVtime
@Title     Get the verification time

@Prototype int taxisInqVtime(int taxisID)
@Parameter
    @Item  taxisID  Time axis ID, from a previous call to @fref{taxisCreate} or @fref{vlistInqTaxis}

@Description
The function @func{taxisInqVtime} returns the verification time of a Time axis.

@Result
@func{taxisInqVtime} returns the verification time.

@EndFunction
*/
int taxisInqVtime(int taxisID)
{
  taxis_t *taxisptr = (taxis_t *)reshGetVal(taxisID, &taxisOps);

  return (taxisptr->vtime);
}


void taxisInqVtimeBounds(int taxisID, int *vtime_lb, int *vtime_ub)
{
  taxis_t *taxisptr = (taxis_t *)reshGetVal(taxisID, &taxisOps);

  *vtime_lb = taxisptr->vtime_lb;
  *vtime_ub = taxisptr->vtime_ub;
}


void taxisDefVtimeBounds(int taxisID, int vtime_lb, int vtime_ub)
{
  taxis_t *taxisptr = ( taxis_t * ) reshGetVal ( taxisID, &taxisOps );

  if (taxisptr->vtime_lb != vtime_lb
      || taxisptr->vtime_ub != vtime_ub
      || taxisptr->has_bounds != TRUE)
    {
      taxisptr->vtime_lb = vtime_lb;
      taxisptr->vtime_ub = vtime_ub;
      taxisptr->has_bounds = TRUE;
      reshSetStatus(taxisID, &taxisOps, RESH_DESYNC_IN_USE);
    }
}

/*
@Function  taxisInqRdate
@Title     Get the reference date

@Prototype int taxisInqRdate(int taxisID)
@Parameter
    @Item  taxisID  Time axis ID, from a previous call to @fref{taxisCreate} or @fref{vlistInqTaxis}

@Description
The function @func{taxisInqRdate} returns the reference date of a Time axis.

@Result
@func{taxisInqRdate} returns the reference date.

@EndFunction
*/
int taxisInqRdate(int taxisID)
{
  taxis_t *taxisptr = (taxis_t *)reshGetVal(taxisID, &taxisOps);

  if ( taxisptr->rdate == -1 )
    {
      taxisptr->rdate = taxisptr->vdate;
      taxisptr->rtime = taxisptr->vtime;
      reshSetStatus(taxisID, &taxisOps, RESH_DESYNC_IN_USE);
    }

  return (taxisptr->rdate);
}

/*
@Function  taxisInqRtime
@Title     Get the reference time

@Prototype int taxisInqRtime(int taxisID)
@Parameter
    @Item  taxisID  Time axis ID, from a previous call to @fref{taxisCreate} or @fref{vlistInqTaxis}

@Description
The function @func{taxisInqRtime} returns the reference time of a Time axis.

@Result
@func{taxisInqRtime} returns the reference time.

@EndFunction
*/
int taxisInqRtime(int taxisID)
{
  taxis_t *taxisptr = (taxis_t *)reshGetVal(taxisID, &taxisOps);

  if ( taxisptr->rdate == -1 )
    {
      taxisptr->rdate = taxisptr->vdate;
      taxisptr->rtime = taxisptr->vtime;
      reshSetStatus(taxisID, &taxisOps, RESH_DESYNC_IN_USE);
    }

  return (taxisptr->rtime);
}

/*
@Function  taxisInqFdate
@Title     Get the forecast reference date

@Prototype int taxisInqFdate(int taxisID)
@Parameter
    @Item  taxisID  Time axis ID, from a previous call to @fref{taxisCreate} or @fref{vlistInqTaxis}

@Description
The function @func{taxisInqFdate} returns the forecast reference date of a Time axis.

@Result
@func{taxisInqFdate} returns the forecast reference date.

@EndFunction
*/
int taxisInqFdate(int taxisID)
{
  taxis_t *taxisptr = (taxis_t *)reshGetVal(taxisID, &taxisOps);

  if ( taxisptr->fdate == -1 )
    {
      taxisptr->fdate = taxisptr->vdate;
      taxisptr->ftime = taxisptr->vtime;
    }

  return (taxisptr->fdate);
}

/*
@Function  taxisInqFtime
@Title     Get the forecast reference time

@Prototype int taxisInqFtime(int taxisID)
@Parameter
    @Item  taxisID  Time axis ID, from a previous call to @fref{taxisCreate} or @fref{vlistInqTaxis}

@Description
The function @func{taxisInqFtime} returns the forecast reference time of a Time axis.

@Result
@func{taxisInqFtime} returns the forecast reference time.

@EndFunction
*/
int taxisInqFtime(int taxisID)
{
  taxis_t *taxisptr = (taxis_t *)reshGetVal(taxisID, &taxisOps);

  if ( taxisptr->fdate == -1 )
    {
      taxisptr->fdate = taxisptr->vdate;
      taxisptr->ftime = taxisptr->vtime;
    }

  return (taxisptr->ftime);
}

/*
@Function  taxisInqCalendar
@Title     Get the calendar

@Prototype int taxisInqCalendar(int taxisID)
@Parameter
    @Item  taxisID  Time axis ID, from a previous call to @fref{taxisCreate} or @fref{vlistInqTaxis}

@Description
The function @func{taxisInqCalendar} returns the calendar of a Time axis.

@Result
@func{taxisInqCalendar} returns the type of the calendar,
one of the set of predefined CDI calendar types.
The valid CDI calendar types are @func{CALENDAR_STANDARD}, @func{CALENDAR_PROLEPTIC},
@func{CALENDAR_360DAYS}, @func{CALENDAR_365DAYS} and @func{CALENDAR_366DAYS}.

@EndFunction
*/
int taxisInqCalendar(int taxisID)
{
  taxis_t *taxisptr = ( taxis_t * ) reshGetVal ( taxisID, &taxisOps );

  return (taxisptr->calendar);
}


int taxisInqTunit(int taxisID)
{
  taxis_t *taxisptr = ( taxis_t * ) reshGetVal ( taxisID, &taxisOps );

  return (taxisptr->unit);
}


int taxisInqForecastTunit(int taxisID)
{
  taxis_t *taxisptr = ( taxis_t * ) reshGetVal ( taxisID, &taxisOps );

  return (taxisptr->fc_unit);
}


double taxisInqForecastPeriod(int taxisID)
{
  taxis_t *taxisptr = ( taxis_t * ) reshGetVal ( taxisID, &taxisOps );

  return (taxisptr->fc_period);
}


int taxisInqNumavg(int taxisID)
{
  taxis_t *taxisptr = ( taxis_t * ) reshGetVal ( taxisID, &taxisOps );

  return (taxisptr->numavg);
}


taxis_t *taxisPtr(int taxisID)
{
  taxis_t *taxisptr = (taxis_t *)reshGetVal(taxisID, &taxisOps);

  return (taxisptr);
}

void
ptaxisDefName(taxis_t *taxisptr, const char *name)
{
  if (name)
    {
      size_t len = strlen(name);
      delete_refcount_string(taxisptr->name);
      char *taxisname = taxisptr->name = new_refcount_string(len);
      strcpy(taxisname, name);
    }
}

void
ptaxisDefLongname(taxis_t *taxisptr, const char *longname)
{
  if (longname)
    {
      size_t len = strlen(longname);
      delete_refcount_string(taxisptr->longname);
      char *taxislongname = taxisptr->longname = new_refcount_string(len);
      strcpy(taxislongname, longname);
    }
}


static void
cdiDecodeTimevalue(int timeunit, double timevalue, int *days, int *secs)
{
  static int lwarn = TRUE;

  *days = 0;
  *secs = 0;

  if ( timeunit == TUNIT_MINUTE )
    {
      timevalue *= 60;
      timeunit = TUNIT_SECOND;
    }
  else if ( timeunit == TUNIT_HOUR )
    {
      timevalue /= 24;
      timeunit = TUNIT_DAY;
    }

  if ( timeunit == TUNIT_SECOND )
    {
      *days = (int) (timevalue/86400);
      double seconds = timevalue - *days*86400.;
      *secs = (int)lround(seconds);
      if ( *secs < 0 ) { *days -= 1; *secs += 86400; };
      /*
      {
	double cval = *days*86400. + *secs;
	if ( cval != timevalue )
	  printf("TUNIT_SECOND error: %g %g %d %d\n", timevalue, cval, *days, *secs);
      }
      */
    }
  else if ( timeunit == TUNIT_DAY )
    {
      *days = (int) timevalue;
      double seconds = (timevalue - *days)*86400;
      *secs = (int)lround(seconds);
      if ( *secs < 0 ) { *days -= 1; *secs += 86400; };
      /*
      {
	double cval = *days + *secs/86400.;
	if ( cval != timevalue )
	  printf("TUNIT_DAY error: %g %g %d %d\n", timevalue, cval, *days, *secs);
      }
      */
    }
  else
    {
      if ( lwarn )
	{
	  Warning("timeunit %s unsupported!", tunitNamePtr(timeunit));
	  lwarn = FALSE;
	}
    }
}

static
void cdiEncodeTimevalue(int days, int secs, int timeunit, double *timevalue)
{
  static int lwarn = TRUE;

  if ( timeunit == TUNIT_SECOND )
    {
      *timevalue = days*86400. + secs;
    }
  else if ( timeunit == TUNIT_MINUTE  ||
	    timeunit == TUNIT_QUARTER ||
	    timeunit == TUNIT_30MINUTES )
    {
      *timevalue = days*1440. + secs/60.;
    }
  else if ( timeunit == TUNIT_HOUR   ||
	    timeunit == TUNIT_3HOURS ||
	    timeunit == TUNIT_6HOURS ||
	    timeunit == TUNIT_12HOURS )
    {
      *timevalue = days*24. + secs/3600.;
    }
  else if ( timeunit == TUNIT_DAY )
    {
      *timevalue = days + secs/86400.;
    }
  else
    {
      if ( lwarn )
	{
	  Warning("timeunit %s unsupported!", tunitNamePtr(timeunit));
	  lwarn = FALSE;
	}
    }
}

int days_per_month(int calendar, int year, int month);

void timeval2vtime(double timevalue, taxis_t *taxis, int *vdate, int *vtime)
{
  int year, month, day, hour, minute, second;
  int rdate, rtime;
  int timeunit;
  int calendar;
  int julday, secofday, days, secs;

  *vdate = 0;
  *vtime = 0;

  timeunit = (*taxis).unit;
  calendar = (*taxis).calendar;

  rdate  = (*taxis).rdate;
  rtime  = (*taxis).rtime;

  if ( rdate == 0 && rtime == 0 && DBL_IS_EQUAL(timevalue, 0.) ) return;

  cdiDecodeDate(rdate, &year, &month, &day);
  cdiDecodeTime(rtime, &hour, &minute, &second);

  if ( timeunit == TUNIT_MONTH && calendar == CALENDAR_360DAYS )
    {
      timeunit = TUNIT_DAY;
      timevalue *= 30;
    }

  if ( timeunit == TUNIT_MONTH || timeunit == TUNIT_YEAR )
    {
      int nmon, dpm;
      double fmon;

      if ( timeunit == TUNIT_YEAR ) timevalue *= 12;

      nmon = (int) timevalue;
      fmon = timevalue - nmon;

      month += nmon;

      while ( month > 12 ) { month -= 12; year++; }
      while ( month <  1 ) { month += 12; year--; }

      dpm = days_per_month(calendar, year, month);
      timeunit = TUNIT_DAY;
      timevalue = fmon*dpm;
    }

  encode_caldaysec(calendar, year, month, day, hour, minute, second, &julday, &secofday);

  cdiDecodeTimevalue(timeunit, timevalue, &days, &secs);

  julday_add(days, secs, &julday, &secofday);

  decode_caldaysec(calendar, julday, secofday, &year, &month, &day, &hour, &minute, &second);

  *vdate = cdiEncodeDate(year, month, day);
  *vtime = cdiEncodeTime(hour, minute, second);
}


double vtime2timeval(int vdate, int vtime, taxis_t *taxis)
{
  int ryear, rmonth;
  int year, month, day, hour, minute, second;
  int rdate, rtime;
  double value = 0;
  int timeunit;
  int timeunit0;
  int calendar;
  int julday1, secofday1, julday2, secofday2, days, secs;

  timeunit = (*taxis).unit;
  calendar = (*taxis).calendar;

  rdate = (*taxis).rdate;
  rtime = (*taxis).rtime;
  if ( rdate == -1 )
    {
      rdate  = (*taxis).vdate;
      rtime  = (*taxis).vtime;
    }

  if ( rdate == 0 && rtime == 0 && vdate == 0 && vtime == 0 ) return(value);

  cdiDecodeDate(rdate, &ryear, &rmonth, &day);
  cdiDecodeTime(rtime, &hour, &minute, &second);

  encode_caldaysec(calendar, ryear, rmonth, day, hour, minute, second, &julday1, &secofday1);

  cdiDecodeDate(vdate, &year, &month, &day);
  cdiDecodeTime(vtime, &hour, &minute, &second);

  timeunit0 = timeunit;

  if ( timeunit == TUNIT_MONTH && calendar == CALENDAR_360DAYS )
    {
      timeunit = TUNIT_DAY;
    }

  if ( timeunit == TUNIT_MONTH || timeunit == TUNIT_YEAR )
    {
      int nmonth, dpm;

      value = (year-ryear)*12 - rmonth + month;

      nmonth = (int) value;
      month -= nmonth;

      while ( month > 12 ) { month -= 12; year++; }
      while ( month <  1 ) { month += 12; year--; }

      dpm = days_per_month(calendar, year, month);

      encode_caldaysec(calendar, year, month, day, hour, minute, second, &julday2, &secofday2);

      julday_sub(julday1, secofday1, julday2, secofday2, &days, &secs);

      value += (days+secs/86400.)/dpm;

      if ( timeunit == TUNIT_YEAR ) value = value/12;
    }
  else
    {
      encode_caldaysec(calendar, year, month, day, hour, minute, second, &julday2, &secofday2);

      julday_sub(julday1, secofday1, julday2, secofday2, &days, &secs);

      cdiEncodeTimevalue(days, secs, timeunit, &value);
    }

  if ( timeunit0 == TUNIT_MONTH && calendar == CALENDAR_360DAYS )
    {
      value /= 30;
    }

  return (value);
}


void conv_timeval(double timevalue, int *rvdate, int *rvtime)
{
  int vdate = 0, vtime = 0;
  int hour, minute, second;
  int daysec;

  vdate = (int) timevalue;
  if ( vdate < 0 )
    daysec = (int) (-(timevalue - vdate)*86400 + 0.01);
  else
    daysec = (int) ( (timevalue - vdate)*86400 + 0.01);

  hour   =  daysec / 3600;
  minute = (daysec - hour*3600)/60;
  second =  daysec - hour*3600 - minute*60;
  vtime  = cdiEncodeTime(hour, minute, second);

  *rvdate = vdate;
  *rvtime = vtime;
}


void splitTimevalue(double timevalue, int timeunit, int *date, int *time)
{
  int vdate = 0, vtime = 0;
  int hour, minute, second;
  int year, month, day;
  static int lwarn = TRUE;

  if ( timeunit == TUNIT_SECOND )
    {
      timevalue /= 86400;
      conv_timeval(timevalue, &vdate, &vtime);
    }
  else if ( timeunit == TUNIT_HOUR )
    {
      timevalue /= 24;
      conv_timeval(timevalue, &vdate, &vtime);
    }
  else if ( timeunit == TUNIT_DAY )
    {
      conv_timeval(timevalue, &vdate, &vtime);
    }
  else if ( timeunit == TUNIT_MONTH )
    {
      vdate = (int) timevalue*100;
      vdate += 1;
      vtime = 0;
    }
  else if ( timeunit == TUNIT_YEAR )
    {
      if ( timevalue < -214700 )
	{
	  Warning("Year %g out of range, set to -214700", timevalue);
	  timevalue = -214700;
	}
      else if ( timevalue > 214700 )
	{
	  Warning("Year %g out of range, set to 214700", timevalue);
	  timevalue = 214700;
	}

      vdate = (int) timevalue*10000;
      vdate += 101;
      vtime = 0;
    }
  else
    {
      if ( lwarn )
	{
	  Warning("timeunit %s unsupported!", tunitNamePtr(timeunit));
	  lwarn = FALSE;
	}
    }

  /* verify date and time */

  cdiDecodeDate(vdate, &year, &month, &day);
  cdiDecodeTime(vtime, &hour, &minute, &second);

  if ( month > 17 || day > 31 || hour > 23 || minute > 59 || second > 59 )
    {
      if ( (month  > 17 || day > 31) && (year < -9999 || year > 9999) ) year = 1;
      if ( month  > 17 ) month  = 1;
      if ( day    > 31 ) day    = 1;
      if ( hour   > 23 ) hour   = 0;
      if ( minute > 59 ) minute = 0;
      if ( second > 59 ) second = 0;

      vdate = cdiEncodeDate(year, month, day);
      vtime = cdiEncodeTime(hour, minute, second);

      if ( lwarn )
        {
          lwarn = FALSE;
          Warning("Reset wrong date/time to %4.4d-%2.2d-%2.2d %2.2d:%2.2d:%2.2d!",
                  year, month, day, hour, minute, second);
        }
    }

  *date = vdate;
  *time = vtime;
}


void cdiSetForecastPeriod(double timevalue, taxis_t *taxis)
{
  int year, month, day, hour, minute, second;
  int vdate, vtime;
  int timeunit;
  int calendar;
  int julday, secofday, days, secs;

  (*taxis).fc_period = timevalue;

  timeunit = (*taxis).fc_unit;
  calendar = (*taxis).calendar;

  vdate  = (*taxis).vdate;
  vtime  = (*taxis).vtime;

  if ( vdate == 0 && vtime == 0 && DBL_IS_EQUAL(timevalue, 0.) ) return;

  cdiDecodeDate(vdate, &year, &month, &day);
  cdiDecodeTime(vtime, &hour, &minute, &second);

  if ( timeunit == TUNIT_MONTH && calendar == CALENDAR_360DAYS )
    {
      timeunit = TUNIT_DAY;
      timevalue *= 30;
    }

  if ( timeunit == TUNIT_MONTH || timeunit == TUNIT_YEAR )
    {
      int nmon, dpm;
      double fmon;

      if ( timeunit == TUNIT_YEAR ) timevalue *= 12;

      nmon = (int) timevalue;
      fmon = timevalue - nmon;

      month -= nmon;

      while ( month > 12 ) { month -= 12; year++; }
      while ( month <  1 ) { month += 12; year--; }

      dpm = days_per_month(calendar, year, month);
      timeunit = TUNIT_DAY;
      timevalue = fmon*dpm;
    }

  encode_caldaysec(calendar, year, month, day, hour, minute, second, &julday, &secofday);

  cdiDecodeTimevalue(timeunit, timevalue, &days, &secs);

  julday_add(-days, -secs, &julday, &secofday);

  decode_caldaysec(calendar, julday, secofday, &year, &month, &day, &hour, &minute, &second);

  (*taxis).fdate = cdiEncodeDate(year, month, day);
  (*taxis).ftime = cdiEncodeTime(hour, minute, second);
}


void cdiDecodeTimeval(double timevalue, taxis_t *taxis, int *date, int *time)
{
  if ( taxis->type == TAXIS_ABSOLUTE )
    splitTimevalue(timevalue, taxis->unit, date, time);
  else
    timeval2vtime(timevalue, taxis, date, time);
}


double cdiEncodeTimeval(int date, int time, taxis_t *taxis)
{
  double timevalue;

  if ( taxis->type == TAXIS_ABSOLUTE )
    {
      if ( taxis->unit == TUNIT_YEAR )
	{
	  int year, month, day;
	  cdiDecodeDate(date, &year, &month, &day);

	  timevalue = year;
	}
      else if ( taxis->unit == TUNIT_MONTH )
	{
	  int year, month, day;
	  cdiDecodeDate(date, &year, &month, &day);
	  if ( day == 0 )
	    timevalue = date/100;
	  else
	    timevalue = date/100 + 0.5;
	}
      else
	{
	  int hour, minute, second;
	  cdiDecodeTime(time, &hour, &minute, &second);
	  if ( date < 0 )
	    timevalue = -(-date + (hour*3600 + minute*60 + second)/86400.);
	  else
	    timevalue =    date + (hour*3600 + minute*60 + second)/86400.;
	}
    }
  else
    timevalue = vtime2timeval(date, time, taxis);

  return (timevalue);
}


void ptaxisInit(taxis_t *taxisptr)
{
  taxisDefaultValue ( taxisptr );
}


void ptaxisCopy(taxis_t *dest, taxis_t *source)
{
  reshLock ();

  /* memcpy(dest, source, sizeof(taxis_t)); */
  dest->used        = source->used;
  dest->type        = source->type;
  dest->vdate       = source->vdate;
  dest->vtime       = source->vtime;
  dest->rdate       = source->rdate;
  dest->rtime       = source->rtime;
  dest->fdate       = source->fdate;
  dest->ftime       = source->ftime;
  dest->calendar    = source->calendar;
  dest->unit        = source->unit;
  dest->numavg      = source->numavg;
  dest->climatology = source->climatology;
  dest->has_bounds  = source->has_bounds;
  dest->vdate_lb    = source->vdate_lb;
  dest->vtime_lb    = source->vtime_lb;
  dest->vdate_ub    = source->vdate_ub;
  dest->vtime_ub    = source->vtime_ub;
  dest->fc_unit     = source->fc_unit;
  dest->fc_period   = source->fc_period;

  dest->climatology = source->climatology;
  delete_refcount_string(dest->name);
  delete_refcount_string(dest->longname);
  dest->name = dup_refcount_string(source->name);
  dest->longname = dup_refcount_string(source->longname);
  if (dest->self != CDI_UNDEFID)
    reshSetStatus(dest->self, &taxisOps, RESH_DESYNC_IN_USE);
  reshUnlock ();
}


static void
taxisPrintKernel(taxis_t * taxisptr, FILE * fp)
{
  int vdate_lb, vdate_ub;
  int vtime_lb, vtime_ub;

  taxisInqVdateBounds ( taxisptr->self, &vdate_lb, &vdate_ub);
  taxisInqVtimeBounds ( taxisptr->self, &vtime_lb, &vtime_ub);

  fprintf ( fp, "#\n");
  fprintf ( fp, "# taxisID %d\n", taxisptr->self);
  fprintf ( fp, "#\n");
  fprintf ( fp, "self        = %d\n", taxisptr->self );
  fprintf ( fp, "used        = %d\n", taxisptr->used );
  fprintf ( fp, "type        = %d\n", taxisptr->type );
  fprintf ( fp, "vdate       = %d\n", taxisptr->vdate );
  fprintf ( fp, "vtime       = %d\n", taxisptr->vtime );
  fprintf ( fp, "rdate       = %d\n", taxisptr->rdate );
  fprintf ( fp, "rtime       = %d\n", taxisptr->rtime );
  fprintf ( fp, "fdate       = %d\n", taxisptr->fdate );
  fprintf ( fp, "ftime       = %d\n", taxisptr->ftime );
  fprintf ( fp, "calendar    = %d\n", taxisptr->calendar );
  fprintf ( fp, "unit        = %d\n", taxisptr->unit );
  fprintf ( fp, "numavg      = %d\n", taxisptr->numavg );
  fprintf ( fp, "climatology = %d\n", taxisptr->climatology );
  fprintf ( fp, "has_bounds  = %d\n", taxisptr->has_bounds );
  fprintf ( fp, "vdate_lb    = %d\n", vdate_lb );
  fprintf ( fp, "vtime_lb    = %d\n", vtime_lb );
  fprintf ( fp, "vdate_ub    = %d\n", vdate_ub );
  fprintf ( fp, "vtime_ub    = %d\n", vtime_ub );
  fprintf ( fp, "fc_unit     = %d\n", taxisptr->fc_unit );
  fprintf ( fp, "fc_period   = %g\n", taxisptr->fc_period );
  fprintf ( fp, "\n");
}

void taxisPrint ( int taxisID )
{
  taxis_t * taxisptr;

  taxisptr = ( taxis_t * ) reshGetVal ( taxisID, &taxisOps );
  taxisPrintKernel ( taxisptr, stdout );
}

static int
taxisCompareP(void *taxisptr1, void *taxisptr2)
{
  taxis_t * t1, * t2;

  t1 = ( taxis_t * ) taxisptr1;
  t2 = ( taxis_t * ) taxisptr2;

  xassert ( t1 );
  xassert ( t2 );

  return ! ( t1->used        == t2->used        &&
	     t1->type        == t2->type        &&
	     t1->vdate       == t2->vdate       &&
	     t1->vtime       == t2->vtime       &&
	     t1->rdate       == t2->rdate       &&
	     t1->rtime       == t2->rtime       &&
	     t1->fdate       == t2->fdate       &&
	     t1->ftime       == t2->ftime       &&
	     t1->calendar    == t2->calendar    &&
	     t1->unit        == t2->unit        &&
	     t1->fc_unit     == t2->fc_unit     &&
	     t1->numavg      == t2->numavg      &&
	     t1->climatology == t2->climatology &&
	     t1->has_bounds  == t2->has_bounds  &&
	     t1->vdate_lb    == t2->vdate_lb    &&
	     t1->vtime_lb    == t2->vtime_lb    &&
	     t1->vdate_ub    == t2->vdate_ub    &&
	     t1->vtime_ub    == t2->vtime_ub );
}


static int
taxisTxCode ( void )
{
  return TAXIS;
}

enum { taxisNint = 21 };

static int
taxisGetPackSize(void *p, void *context)
{
  taxis_t *taxisptr = (taxis_t*) p;
  int packBufferSize
    = serializeGetSize(taxisNint, DATATYPE_INT, context)
    + serializeGetSize(1, DATATYPE_UINT32, context)
    + (taxisptr->name ?
       serializeGetSize((int)strlen(taxisptr->name), DATATYPE_TXT, context) : 0)
    + (taxisptr->longname ?
       serializeGetSize((int)strlen(taxisptr->longname), DATATYPE_TXT,
                        context) : 0);
  return packBufferSize;
}

int
taxisUnpack(char * unpackBuffer, int unpackBufferSize, int * unpackBufferPos,
            int originNamespace, void *context, int force_id)
{
  taxis_t * taxisP;
  int intBuffer[taxisNint];
  uint32_t d;
  int idx = 0;

  serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos,
                  intBuffer, taxisNint, DATATYPE_INT, context);
  serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos,
                  &d, 1, DATATYPE_UINT32, context);

  xassert(cdiCheckSum(DATATYPE_INT, taxisNint, intBuffer) == d);

  taxisInit();

  cdiResH targetID = namespaceAdaptKey(intBuffer[idx++], originNamespace);
  taxisP = taxisNewEntry(force_id?targetID:CDI_UNDEFID);

  xassert(!force_id || targetID == taxisP->self);

  taxisP->used        = intBuffer[idx++];
  taxisP->type        = intBuffer[idx++];
  taxisP->vdate       = intBuffer[idx++];
  taxisP->vtime       = intBuffer[idx++];
  taxisP->rdate       = intBuffer[idx++];
  taxisP->rtime       = intBuffer[idx++];
  taxisP->fdate       = intBuffer[idx++];
  taxisP->ftime       = intBuffer[idx++];
  taxisP->calendar    = intBuffer[idx++];
  taxisP->unit        = intBuffer[idx++];
  taxisP->fc_unit     = intBuffer[idx++];
  taxisP->numavg      = intBuffer[idx++];
  taxisP->climatology = intBuffer[idx++];
  taxisP->has_bounds  = intBuffer[idx++];
  taxisP->vdate_lb    = intBuffer[idx++];
  taxisP->vtime_lb    = intBuffer[idx++];
  taxisP->vdate_ub    = intBuffer[idx++];
  taxisP->vtime_ub    = intBuffer[idx++];

  if (intBuffer[idx])
    {
      int len = intBuffer[idx];
      char *name = new_refcount_string((size_t)len);
      serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos,
                      name, len, DATATYPE_TXT, context);
      name[len] = '\0';
      taxisP->name = name;
    }
  idx++;
  if (intBuffer[idx])
    {
      int len = intBuffer[idx];
      char *longname = new_refcount_string((size_t)len);
      serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos,
                      longname, len, DATATYPE_TXT, context);
      longname[len] = '\0';
      taxisP->longname = longname;
    }

  return taxisP->self;
}


static void
taxisPack(void * voidP, void * packBuffer, int packBufferSize, int * packBufferPos,
          void *context)
{
  taxis_t *taxisP = (taxis_t *)voidP;
  int intBuffer[taxisNint];
  uint32_t d;
  int idx = 0;

  intBuffer[idx++] = taxisP->self;
  intBuffer[idx++] = taxisP->used;
  intBuffer[idx++] = taxisP->type;
  intBuffer[idx++] = taxisP->vdate;
  intBuffer[idx++] = taxisP->vtime;
  intBuffer[idx++] = taxisP->rdate;
  intBuffer[idx++] = taxisP->rtime;
  intBuffer[idx++] = taxisP->fdate;
  intBuffer[idx++] = taxisP->ftime;
  intBuffer[idx++] = taxisP->calendar;
  intBuffer[idx++] = taxisP->unit;
  intBuffer[idx++] = taxisP->fc_unit;
  intBuffer[idx++] = taxisP->numavg;
  intBuffer[idx++] = taxisP->climatology;
  intBuffer[idx++] = taxisP->has_bounds;
  intBuffer[idx++] = taxisP->vdate_lb;
  intBuffer[idx++] = taxisP->vtime_lb;
  intBuffer[idx++] = taxisP->vdate_ub;
  intBuffer[idx++] = taxisP->vtime_ub;
  intBuffer[idx++] = taxisP->name ? (int)strlen(taxisP->name) : 0;
  intBuffer[idx++] = taxisP->longname ? (int)strlen(taxisP->longname) : 0;

  serializePack(intBuffer, taxisNint, DATATYPE_INT,
                packBuffer, packBufferSize, packBufferPos, context);
  d = cdiCheckSum(DATATYPE_INT, taxisNint, intBuffer);
  serializePack(&d, 1, DATATYPE_UINT32,
                packBuffer, packBufferSize, packBufferPos, context);
  if (taxisP->name)
    serializePack(taxisP->name, intBuffer[15], DATATYPE_TXT,
                  packBuffer, packBufferSize, packBufferPos, context);
  if (taxisP->longname)
    serializePack(taxisP->longname, intBuffer[16], DATATYPE_TXT,
                  packBuffer, packBufferSize, packBufferPos, context);

}

/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
