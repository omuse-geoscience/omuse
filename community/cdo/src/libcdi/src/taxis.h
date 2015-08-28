#ifndef _TAXIS_H
#define _TAXIS_H

#ifndef RESOURCE_HANDLE_H
#include "resource_handle.h"
#endif

typedef struct {
  /* Date format  YYYYMMDD */
  /* Time format    hhmmss */
  int     self;
  int     used;
  int     type;           // time type
  int     vdate;          // verification date
  int     vtime;          // verification time
  int     rdate;          // reference date
  int     rtime;          // reference time
  int     fdate;          // forecast reference date
  int     ftime;          // forecast reference time
  int     calendar;
  int     unit;           // time unit
  int     numavg;
  int     climatology;
  int     has_bounds;
  int     vdate_lb;       // lower bounds of vdate
  int     vtime_lb;       // lower bounds of vtime
  int     vdate_ub;       // upper bounds of vdate
  int     vtime_ub;       // upper bounds of vtime
  int     fc_unit;        // forecast time unit
  double  fc_period;      // forecast time period
  char*   name;
  char*   longname;
}
taxis_t;

void     ptaxisInit(taxis_t* taxis);
void     ptaxisCopy(taxis_t* dest, taxis_t* source);
taxis_t* taxisPtr(int taxisID);
void     cdiSetForecastPeriod(double timevalue, taxis_t *taxis);
void     cdiDecodeTimeval(double timevalue, taxis_t* taxis, int* date, int* time);
double   cdiEncodeTimeval(int date, int time, taxis_t* taxis);
void     timeval2vtime(double timevalue, taxis_t* taxis, int* vdate, int* vtime);
double   vtime2timeval(int vdate, int vtime, taxis_t *taxis);

void    ptaxisDefName(taxis_t *taxisptr, const char *name);
void    ptaxisDefLongname(taxis_t *taxisptr, const char *name);
void    taxisDestroyKernel(taxis_t *taxisptr);
#if !defined (SX)
extern const resOps taxisOps;
#endif

int
taxisUnpack(char *unpackBuffer, int unpackBufferSize, int *unpackBufferPos,
            int originNamespace, void *context, int checkForSameID);

#endif  /* _TAXIS_H */
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
