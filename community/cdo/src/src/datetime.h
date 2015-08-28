#include <stdio.h>

#define  TIMESTAT_FIRST  1
#define  TIMESTAT_LAST   2
#define  TIMESTAT_MEAN   3


typedef struct {
  int   julday;
  int   secofday;
} juldate_t;


typedef struct {
  int   date;
  int   time;
} datetime_t;


typedef struct {
  int   date;
  int   time;
} datetime_type;

typedef struct
{
  datetime_type v;
  datetime_type b[2];
} dtinfo_type;

typedef struct
{
  size_t       nalloc;
  size_t       size;
  int          has_bounds;
  int          calendar;
  int          stat;
  int          timestat_date;
  dtinfo_type  timestat;
  dtinfo_type *dtinfo;
} dtlist_type;



juldate_t juldate_encode(int calendar, int date, int time);
void      juldate_decode(int calendar, juldate_t juldate, int *date, int *time);
juldate_t juldate_sub(juldate_t juldate2, juldate_t juldate1);
juldate_t juldate_add_seconds(int seconds, juldate_t juldate);
double    juldate_to_seconds(juldate_t juldate);


void    datetime_avg(int dpy, int ndates, datetime_t *datetime);

dtlist_type *dtlist_new(void);
void dtlist_delete(dtlist_type *dtlist);
void dtlist_shift(dtlist_type *dtlist);
void dtlist_set_stat(dtlist_type *dtlist, int stat);
void dtlist_set_calendar(dtlist_type *dtlist, int calendar);
int  dtlist_get_vdate(dtlist_type *dtlist, int tsID);
int  dtlist_get_vtime(dtlist_type *dtlist, int tsID);
void dtlist_taxisInqTimestep(dtlist_type *dtlist, int taxisID, int tsID);
void dtlist_taxisDefTimestep(dtlist_type *dtlist, int taxisID, int tsID);
void dtlist_stat_taxisDefTimestep(dtlist_type *dtlist, int taxisID, int nsteps);
