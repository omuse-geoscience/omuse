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

#ifndef _PIPE_H
#define _PIPE_H

#if defined(HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <sys/types.h>

#if defined(HAVE_LIBPTHREAD)

#include <pthread.h>
#include "pthread_debug.h"

#endif

typedef struct {
  short       check_datarange;
  int         gridsize;
  int         datatype;
  double      missval;
  double      addoffset;
  double      scalefactor;
} varlist_t;


typedef struct {
  int            self;
  int            mode;
  int            fileID;
  int            vlistID;
  int            tsID;
  int            filetype;
  int            ispipe;
  int            isopen;
  int            tsID0;
  int            mfiles;
  int            nfiles;
  int            varID;           /* next varID defined with streamDefVar */
  char          *name;
  char         **mfnames;
  varlist_t     *varlist;
#if defined(HAVE_LIBPTHREAD)
  void          *argument;
  struct pipe_s *pipe;
  pthread_t     rthreadID; /* read  thread ID */
  pthread_t     wthreadID; /* write thread ID */
#endif
} pstream_t;


#if defined(HAVE_LIBPTHREAD)

struct pipe_s {
  int     nrecs, EOP;
  int     varID, levelID;
  int     recIDr, recIDw, tsIDr, tsIDw;
  int     hasdata, usedata;
  int     nmiss;
  double *data;
  pstream_t *pstreamptr_in;
  /* unsigned long */ off_t nvals;
  pthread_mutex_t *mutex;
  pthread_cond_t *tsDef, *tsInq, *vlistDef, *isclosed;
  pthread_cond_t *recDef, *recInq;
  pthread_cond_t *writeCond, *readCond;
};

typedef struct pipe_s pipe_t;

pipe_t *pipeNew(void);
void    pipeDelete(pipe_t *pipe);

void  pipeDebug(int debug);

void  pipeDefVlist(pstream_t *pstreamptr, int vlistID);
int   pipeInqVlist(pstream_t *pstreamptr);

void  pipeDefTimestep(pstream_t *pstreamptr, int tsID);
int   pipeInqTimestep(pstream_t *pstreamptr, int tsID);

void  pipeDefRecord(pstream_t *pstreamptr, int  varID, int  levelID);
int   pipeInqRecord(pstream_t *pstreamptr, int *varID, int *levelID);

void  pipeReadRecord(pstream_t *pstreamptr, double *data, int *nmiss);
void  pipeWriteRecord(pstream_t *pstreamptr, double *data, int nmiss);
void  pipeCopyRecord(pstream_t *pstreamptr_dest, pstream_t *pstreamptr_src);

#endif

#endif  /* _PIPE_H */
