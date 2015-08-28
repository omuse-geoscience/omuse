#ifndef _BASETIME_H
#define _BASETIME_H

//#define USE_TIMECACHE 1
#define MAX_TIMECACHE_SIZE 1024

typedef struct {
  int size;
  int startid;
  int maxvals;
  double cache[MAX_TIMECACHE_SIZE];
}
timecache_t;

typedef struct {
  int   ncvarid;
  int   ncdimid;
  int   ncvarboundsid;
  int   leadtimeid;
  int   lwrf;     /* TRUE for time axis in WRF format */
  timecache_t *timevar_cache;
}
basetime_t;

void basetimeInit(basetime_t *basetime);

#endif  /* _BASETIME_H */
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
