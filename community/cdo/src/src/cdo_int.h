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

#ifndef _CDO_INT_H
#define _CDO_INT_H

#if defined(HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "compare.h"
#include "timebase.h"
#include "field.h"
#include "functs.h"
#include "dmemory.h"
#include "process.h"
#include "const.h"
#include "util.h"
#include "datetime.h"

#define  OPENMP4  201307
#if defined(_OPENMP) && defined(OPENMP4) && _OPENMP >= OPENMP4
#define  HAVE_OPENMP4  1
#endif


#ifndef strdupx
#ifndef strdup
char *strdup(const char *s);
#endif
#define strdupx  strdup
/*
#define strdupx(s)			          \
({					      	  \
   const char *__old = (s);			  \
   size_t __len = strlen(__old) + 1;		  \
   char *__new = malloc(__len);	  \
   (char *) memcpy(__new, __old, __len);	  \
})
*/
#endif


#define  cmpstr(s1, s2)          (strncmp(s1, s2, strlen(s2)))
#define  cmpstrlen(s1, s2, len)  (strncmp(s1, s2, len = strlen(s2)))


/* sxxxYYYYMMDDhhmm0 */
#define  DATE_LEN  31        /* YYYYMMDDhhmmss allocate DTLEN+1 !!!! */
#define  SET_DATE(dtstr, date, time)      (sprintf(dtstr, "%*d%*d", DATE_LEN-6, date, 6, time))
#define  DATE_IS_NEQ(dtstr1, dtstr2, len) (memcmp(dtstr1, dtstr2, len) != 0)

enum T_WEIGHT_MODE {WEIGHT_OFF, WEIGHT_ON};
enum T_EIGEN_MODE  {JACOBI, DANIELSON_LANCZOS};


#ifndef  M_LN10
#define  M_LN10      2.30258509299404568402  /* log_e 10 */
#endif

#ifndef  M_PI
#define  M_PI        3.14159265358979323846  /* pi */
#endif


#define  IX2D(y,x,nx)  ((y)*(nx)+(x))

#define  MEMTYPE_DOUBLE  1
#define  MEMTYPE_FLOAT   2

#define  CDO_EXP_LOCAL   1
#define  CDO_EXP_REMOTE  2

void print_pthread_info(void);

void cdoProcessTime(double *utime, double *stime);

void    setCommandLine(int argc, char **argv);
char   *commandLine(void);
int     readline(FILE *fp, char *line, int len);

int zaxis2ltype(int zaxisID);


int nfc2nlat(int nfc, int ntr);
int nlat2ntr(int nlat);
int nlat2ntr_linear(int nlat);
int ntr2nlat(int ntr);
int ntr2nlat_linear(int ntr);
int compNlon(int nlat);

void param2str(int param, char *paramstr, int maxlen);
void datetime2str(int date, int time, char *datetimestr, int maxlen);
void date2str(int date, char *datestr, int maxlen);
void time2str(int time, char *timestr, int maxlen);

const char * tunit2str(int tunits);
const char * calendar2str(int calendar);

int     days_per_month(int calendar, int year, int month);
int     days_per_year(int calendar, int year);
int     calendar_dpy(int calendar);

void    defineGrid(const char *gridarg);
void    defineInstitution(char *instarg);
int     defineTable(char *tablearg);

void    cdolog(const char *prompt, double cputime);
void    cdologs(int noper);
void    cdologo(int noper);
void    nospec(int vlistID);
void    gridWrite(FILE *fp, int gridID);

void openLock(void);
void openUnlock(void);

int  cdf_openread(const char *filename);

void printFiletype(int streamID, int vlistID);

void job_submit(const char *expname, const char *jobfilename, const char *jobname, const char *tmppath, const char *ftppath);

void minmaxval(long nvals, double *array, int *imiss, double *minval, double *maxval);

off_t fileSize(const char *restrict filename);

char *expand_filename(const char *string);

double parameter2double(const char *string);
int    parameter2int(const char *string);
int    parameter2intlist(const char *string);


#endif  /* _CDO_INT_H */
