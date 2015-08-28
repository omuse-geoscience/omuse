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

#if defined(_OPENMP)
#  include <omp.h>
#endif

#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <errno.h>
#include <sys/stat.h> /* stat */

FILE *popen(const char *command, const char *type);
int pclose(FILE *stream);

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "modules.h"
#include "pstream_int.h"
#include "cdo_int.h"
#include "util.h"
#include "pipe.h"
#include "error.h"
#include "dmemory.h"


static int PSTREAM_Debug = 0;

#define  MAX_PSTREAMS  4096

static int _pstream_max = MAX_PSTREAMS;

static void pstream_initialize(void);

static int _pstream_init = FALSE;

#if defined(HAVE_LIBPTHREAD)
#include <pthread.h>
#include "pthread_debug.h"

// TODO: make threadsafe
static int pthreadScope = 0;

static pthread_mutex_t streamOpenReadMutex  = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t streamOpenWriteMutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t streamMutex          = PTHREAD_MUTEX_INITIALIZER;

static pthread_once_t _pstream_init_thread = PTHREAD_ONCE_INIT;
static pthread_mutex_t _pstream_mutex;

#  define PSTREAM_LOCK()           pthread_mutex_lock(&_pstream_mutex)
#  define PSTREAM_UNLOCK()         pthread_mutex_unlock(&_pstream_mutex)
#  define PSTREAM_INIT()	  \
   if ( _pstream_init == FALSE ) pthread_once(&_pstream_init_thread, pstream_initialize)

#else

#  define PSTREAM_LOCK()
#  define PSTREAM_UNLOCK()
#  define PSTREAM_INIT()	  \
   if ( _pstream_init == FALSE ) pstream_initialize()

#endif


typedef struct _pstreamPtrToIdx {
  int idx;
  pstream_t *ptr;
  struct _pstreamPtrToIdx *next;
} pstreamPtrToIdx;


static pstreamPtrToIdx *_pstreamList  = NULL;
static pstreamPtrToIdx *_pstreamAvail = NULL;


static
void pstream_list_new(void)
{
  assert(_pstreamList == NULL);

  _pstreamList = (pstreamPtrToIdx*) malloc(_pstream_max*sizeof(pstreamPtrToIdx));
}

static
void pstream_list_delete(void)
{
  if ( _pstreamList ) free(_pstreamList);
}

static
void pstream_init_pointer(void)
{  
  for ( int i = 0; i < _pstream_max; ++i )
    {
      _pstreamList[i].next = _pstreamList + i + 1;
      _pstreamList[i].idx  = i;
      _pstreamList[i].ptr  = 0;
    }

  _pstreamList[_pstream_max-1].next = 0;

  _pstreamAvail = _pstreamList;
}

static
pstream_t *pstream_to_pointer(int idx)
{
  pstream_t *pstreamptr = NULL;

  PSTREAM_INIT();

  if ( idx >= 0 && idx < _pstream_max )
    {
      PSTREAM_LOCK();

      pstreamptr = _pstreamList[idx].ptr;

      PSTREAM_UNLOCK();
    }
  else
    Error("pstream index %d undefined!", idx);

  return (pstreamptr);
}

/* Create an index from a pointer */
static
int pstream_from_pointer(pstream_t *ptr)
{
  int idx = -1;

  if ( ptr )
    {
      PSTREAM_LOCK();

      if ( _pstreamAvail )
	{
	  pstreamPtrToIdx *newptr = _pstreamAvail;
	  _pstreamAvail = _pstreamAvail->next;
	  newptr->next  = 0;
	  idx	        = newptr->idx;
	  newptr->ptr   = ptr;
      
	  if ( PSTREAM_Debug )
	    Message("Pointer %p has idx %d from pstream list", ptr, idx);
	}
      else
	Error("Too many open pstreams (limit is %d)!", _pstream_max);

      PSTREAM_UNLOCK();
    }
  else
    Error("Internal problem (pointer %p undefined)", ptr);

  return (idx);
}

static
void pstream_init_entry(pstream_t *pstreamptr)
{
  pstreamptr->self       = pstream_from_pointer(pstreamptr);

  pstreamptr->isopen     = TRUE;
  pstreamptr->ispipe     = FALSE;
  pstreamptr->fileID     = -1;
  pstreamptr->vlistID    = -1;
  pstreamptr->tsID       = -1;
  pstreamptr->filetype   = -1;
  pstreamptr->name       = NULL;
  pstreamptr->tsID0      = 0;
  pstreamptr->mfiles     = 0;
  pstreamptr->nfiles     = 0;
  pstreamptr->varID      = -1;
  pstreamptr->name       = NULL;
  pstreamptr->mfnames    = NULL;
  pstreamptr->varlist    = NULL;
#if defined(HAVE_LIBPTHREAD)
  pstreamptr->argument   = NULL;
  pstreamptr->pipe       = NULL;
  //  pstreamptr->rthreadID  = 0;
  //  pstreamptr->wthreadID  = 0;
#endif
}

static
pstream_t *pstream_new_entry(void)
{
  pstream_t *pstreamptr = (pstream_t*) malloc(sizeof(pstream_t));

  if ( pstreamptr ) pstream_init_entry(pstreamptr);

  return (pstreamptr);
}

static
void pstream_delete_entry(pstream_t *pstreamptr)
{
  int idx = pstreamptr->self;

  PSTREAM_LOCK();

  free(pstreamptr);

  _pstreamList[idx].next = _pstreamAvail;
  _pstreamList[idx].ptr  = 0;
  _pstreamAvail   	 = &_pstreamList[idx];

  PSTREAM_UNLOCK();

  if ( PSTREAM_Debug )
    Message("Removed idx %d from pstream list", idx);
}

static
void pstream_initialize(void)
{
#if defined(HAVE_LIBPTHREAD)
  /* initialize global API mutex lock */
  pthread_mutex_init(&_pstream_mutex, NULL);
#endif

  char *env = getenv("PSTREAM_DEBUG");
  if ( env ) PSTREAM_Debug = atoi(env);

  env = getenv("PSTREAM_MAX");
  if ( env ) _pstream_max = atoi(env);

  if ( PSTREAM_Debug )
    Message("PSTREAM_MAX = %d", _pstream_max);

  pstream_list_new();
  atexit(pstream_list_delete);

  pstream_init_pointer();

  _pstream_init = TRUE;
}

static
int pstreamFindID(const char *name)
{
  pstream_t *pstreamptr;
  int pstreamID;

  for ( pstreamID = 0; pstreamID < _pstream_max; ++pstreamID )
    {
      pstreamptr = pstream_to_pointer(pstreamID);

      if ( pstreamptr )
	if ( pstreamptr->name )
	  if ( strcmp(pstreamptr->name, name) == 0 ) break;
    }

  if ( pstreamID == _pstream_max ) pstreamID = -1;

  return (pstreamID);
}


int pstreamIsPipe(int pstreamID)
{
  pstream_t *pstreamptr = pstream_to_pointer(pstreamID);

  return (pstreamptr->ispipe);
}


int pstreamOpenRead(const argument_t *argument)
{
  PSTREAM_INIT();

  pstream_t *pstreamptr = pstream_new_entry();
  if ( ! pstreamptr ) Error("No memory");

  int pstreamID = pstreamptr->self;

  int ispipe = argument->args[0] == '-';
  /*
  printf("pstreamOpenRead: args >%s<\n", argument->args);
  for ( int i = 0; i < argument->argc; ++i )
    printf("pstreamOpenRead: arg %d >%s<\n", i, argument->argv[i]);
  */
  if ( ispipe )
    {
#if defined(HAVE_LIBPTHREAD)
      char *operatorArg;
      char *operatorName;
      char *newarg;
      char *pipename = (char*) malloc(16);
      int rval;
      pthread_t thrID;
      pthread_attr_t attr;
      // struct sched_param param;
      size_t len;
      size_t stacksize;
      int status;
      argument_t *newargument = (argument_t*) malloc(sizeof(argument_t));

      newargument->argc = argument->argc + 1;
      newargument->argv = (char **) malloc(newargument->argc*sizeof(char *));
      memcpy(newargument->argv, argument->argv, argument->argc*sizeof(char *));

      operatorArg  = argument->argv[0];
      operatorName = getOperatorName(operatorArg);

      len = strlen(argument->args);
      newarg = (char*) malloc(len+16);
      strcpy(newarg, argument->args);
      sprintf(pipename, "(pipe%d.%d)", processSelf() + 1, processInqChildNum() + 1);
      newarg[len] = ' ';
      strcpy(&newarg[len+1], pipename);

      newargument->argv[argument->argc] = pipename;
      newargument->args = newarg;
      /*
      printf("pstreamOpenRead: new args >%s<\n", newargument->args);
      for ( int i = 0; i < newargument->argc; ++i )
	printf("pstreamOpenRead: new arg %d >%s<\n", i, newargument->argv[i]);
      */
      pstreamptr->ispipe    = TRUE;
      pstreamptr->name      = pipename;
      pstreamptr->rthreadID = pthread_self();
      pstreamptr->pipe      = pipeNew();
      pstreamptr->argument  = (void *) newargument;
 
      if ( ! cdoSilentMode )
	cdoPrint("Started child process \"%s\".", newarg+1);

      status = pthread_attr_init(&attr);
      if ( status ) SysError("pthread_attr_init failed for '%s'\n", newarg+1);
      status = pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
      if ( status ) SysError("pthread_attr_setdetachstate failed for '%s'\n", newarg+1);
      /*
      param.sched_priority = 0;
      status = pthread_attr_setschedparam(&attr, &param);
      if ( status ) SysError("pthread_attr_setschedparam failed for '%s'\n", newarg+1);
      */
      /* status = pthread_attr_setinheritsched(&attr, PTHREAD_EXPLICIT_SCHED); */
      /* if ( status ) SysError("pthread_attr_setinheritsched failed for '%s'\n", newarg+1); */

      pthread_attr_getscope(&attr, &pthreadScope);

      /* status = pthread_attr_setscope(&attr, PTHREAD_SCOPE_PROCESS); */
      /* if ( status ) SysError("pthread_attr_setscope failed for '%s'\n", newarg+1); */
      /* If system scheduling scope is specified, then the thread is scheduled against all threads in the system */
      /* pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM); */

      status = pthread_attr_getstacksize(&attr, &stacksize);
      if ( stacksize < 2097152 )
	{
	  stacksize = 2097152;
	  pthread_attr_setstacksize(&attr, stacksize);
	}

      rval = pthread_create(&thrID, &attr, operatorModule(operatorName), newargument);
      if ( rval != 0 )
	{
	  errno = rval;
	  SysError("pthread_create failed for '%s'\n", newarg+1);
	}

      /* free(operatorName); */
      processAddStream(pstreamID);
      /*      pipeInqInfo(pstreamID); */
      if ( PSTREAM_Debug ) Message("pipe %s", pipename);
#else
      cdoAbort("Cannot use pipes, pthread support not compiled in!");
#endif
    }
  else
    {
      size_t len, i;
      int nfiles = 1, j;
      char *filename = NULL;
      const char *pch;

      len = strlen(argument->args);

      for ( i = 0; i < len; i++ )
	if ( argument->args[i] == ':' ) break;

      if ( i < len )
	{
	  pch = &argument->args[i+1];
	  len -= (i+1);
	  if ( len && ( strncmp(argument->args, "filelist:", 9) == 0 || 
			strncmp(argument->args, "flist:", 6) == 0 ) )
	    {
	      for ( i = 0; i < len; i++ ) if ( pch[i] == ',' ) nfiles++;

	      if ( nfiles == 1 )
		{
		  char line[4096];
		  FILE *fp, *fp2;
		  fp = fopen(pch, "r");
		  if ( fp == NULL ) cdoAbort("Open failed on %s", pch);

		  if ( cdoVerbose )
		    cdoPrint("Reading file names from %s", pch);

		  /* find number of files */
		  nfiles = 0;
		  while ( readline(fp, line, 4096) )
		    {
		      if ( line[0] == '#' || line[0] == '\0' || line[0] == ' ' ) continue;

		      fp2 = fopen(line, "r" );
		      if ( fp2 == NULL ) cdoAbort("Open failed on %s", line);
		      fclose(fp2);
		      nfiles++;
		      if ( cdoVerbose )
			cdoPrint("File number %d is %s", nfiles, line);
		    }

		  if ( nfiles == 0 ) cdoAbort("No imput file found in %s", pch);

		  pstreamptr->mfiles = nfiles;
		  pstreamptr->mfnames = (char **) malloc(nfiles*sizeof(char *));
		  
		  rewind(fp);

		  nfiles = 0;
		  while ( readline(fp, line, 4096) )
		    {
		      if ( line[0] == '#' || line[0] == '\0' ||
			   line[0] == ' ' ) continue;

		      pstreamptr->mfnames[nfiles] = strdupx(line);
		      nfiles++;
		    }

		  fclose(fp);
		}
	      else
		{
		  char line[65536];

		  pstreamptr->mfiles = nfiles;
		  pstreamptr->mfnames = (char **) malloc(nfiles*sizeof(char *));
		  
		  strcpy(line, pch);
		  for ( i = 0; i < len; i++ ) if ( line[i] == ',' ) line[i] = 0;

		  i = 0;
		  for ( j = 0; j < nfiles; j++ )
		    {
		      pstreamptr->mfnames[j] = strdupx(&line[i]);
		      i += strlen(&line[i]) + 1;
		    }
		}
	    }
	  else if ( len && strncmp(argument->args, "ls:", 3) == 0 )
	    {
	      char line[4096];
	      char command[4096];
	      char *fnames[16384];
	      FILE *pfp;

	      strcpy(command, "ls ");
	      strcat(command, pch);

	      pfp = popen(command, "r");
	      if ( pfp == 0 )
		SysError("popen %s failed", command);

	      nfiles = 0;
	      while ( readline(pfp, line, 4096) )
		{
		  if ( nfiles >= 16384 ) cdoAbort("Too many input files (limit: 16384)");
		  fnames[nfiles++] = strdupx(line);
		}

	      pclose(pfp);

	      pstreamptr->mfiles = nfiles;
	      pstreamptr->mfnames = (char **) malloc(nfiles*sizeof(char *));

	      for ( j = 0; j < nfiles; j++ )
		pstreamptr->mfnames[j] = fnames[j];
	    }
	}

      if ( pstreamptr->mfiles )
	{
	  len = strlen(pstreamptr->mfnames[0]);
	  filename = (char*) malloc(len+1);
	  strcpy(filename, pstreamptr->mfnames[0]);
	  pstreamptr->nfiles = 1;
	}
      else
	{
	  len = strlen(argument->args);
	  filename = (char*) malloc(len+1);
	  strcpy(filename, argument->args);
	}

      if ( PSTREAM_Debug ) Message("file %s", filename);

#if defined(HAVE_LIBPTHREAD)
      if ( cdoLockIO )
	pthread_mutex_lock(&streamMutex);
      else
	pthread_mutex_lock(&streamOpenReadMutex);
#endif
      int fileID = streamOpenRead(filename);
      if ( fileID < 0 ) cdiOpenError(fileID, "Open failed on >%s<", filename);

      if ( cdoDefaultFileType == CDI_UNDEFID )
	cdoDefaultFileType = streamInqFiletype(fileID);
      /*
      if ( cdoDefaultInstID == CDI_UNDEFID )
	cdoDefaultInstID = streamInqInstID(fileID);
      */
      cdoInqHistory(fileID);
#if defined(HAVE_LIBPTHREAD)
      if ( cdoLockIO )
	pthread_mutex_unlock(&streamMutex);
      else
	pthread_mutex_unlock(&streamOpenReadMutex);
#endif

      pstreamptr->mode   = 'r';
      pstreamptr->name   = filename;
      pstreamptr->fileID = fileID;
    }

  if ( pstreamID < 0 ) cdiOpenError(pstreamID, "Open failed on >%s<", argument->args);
  
  return (pstreamID);
}

static
void query_user_exit(const char *argument)
{
  /* modified code from NCO */
#define USR_RPL_MAX_LNG 10 /* Maximum length for user reply */
#define USR_RPL_MAX_NBR 10 /* Maximum number of chances for user to reply */
  char usr_rpl[USR_RPL_MAX_LNG];
  int usr_rpl_int;
  short nbr_itr=0;
  size_t usr_rpl_lng = 0;

  /* Initialize user reply string */
  usr_rpl[0]='z';
  usr_rpl[1]='\0';

  while ( !(usr_rpl_lng == 1 && 
	    (*usr_rpl == 'o' || *usr_rpl == 'O' || *usr_rpl == 'e' || *usr_rpl == 'E')) )
    {
      if ( nbr_itr++ > USR_RPL_MAX_NBR )
	{
	  (void)fprintf(stdout,"\n%s: ERROR %d failed attempts to obtain valid interactive input.\n",
			processInqPrompt(), nbr_itr-1);
	  exit(EXIT_FAILURE);
	}

      if ( nbr_itr > 1 ) (void)fprintf(stdout,"%s: ERROR Invalid response.\n", processInqPrompt());
      (void)fprintf(stdout,"%s: %s exists ---`e'xit, or `o'verwrite (delete existing file) (e/o)? ",
		    processInqPrompt(), argument);
      (void)fflush(stdout);
      if ( fgets(usr_rpl, USR_RPL_MAX_LNG, stdin) == NULL ) continue;

      /* Ensure last character in input string is \n and replace that with \0 */
      usr_rpl_lng = strlen(usr_rpl);
      if ( usr_rpl_lng >= 1 )
	if ( usr_rpl[usr_rpl_lng-1] == '\n' )
	  {
	    usr_rpl[usr_rpl_lng-1] = '\0';
	    usr_rpl_lng--;
	  }
    }

  /* Ensure one case statement for each exit condition in preceding while loop */
  usr_rpl_int=(int)usr_rpl[0];
  switch(usr_rpl_int)
    {
    case 'E':
    case 'e':
      exit(EXIT_SUCCESS);
      break;
    case 'O':
    case 'o':
      break;
    default:
      exit(EXIT_FAILURE);
      break;
    } /* end switch */
}


int pstreamOpenWrite(const argument_t *argument, int filetype)
{
  int pstreamID = -1;
  pstream_t *pstreamptr;

  PSTREAM_INIT();

  int ispipe = strncmp(argument->args, "(pipe", 5) == 0;

  if ( ispipe )
    {
#if defined(HAVE_LIBPTHREAD)
      if ( PSTREAM_Debug ) Message("pipe %s", argument->args);
      pstreamID = pstreamFindID(argument->args);
      if ( pstreamID == -1 ) Error("%s is not open!", argument->args);

      pstreamptr = pstream_to_pointer(pstreamID);

      pstreamptr->wthreadID = pthread_self();
      pstreamptr->filetype = filetype;
      processAddStream(pstreamID);
#endif
    }
  else
    {
      char *filename = (char*) malloc(strlen(argument->args)+1);

      pstreamptr = pstream_new_entry();
      if ( ! pstreamptr ) Error("No memory");

      pstreamID = pstreamptr->self;
  
      if ( PSTREAM_Debug ) Message("file %s", argument->args);

      if ( filetype == CDI_UNDEFID ) filetype = FILETYPE_GRB;

      if ( cdoInteractive )
	{
	  int rstatus;
	  struct stat stbuf;

	  rstatus = stat(argument->args, &stbuf);
	  /* If permanent file already exists, query user whether to overwrite or exit */
	  if ( rstatus != -1 ) query_user_exit(argument->args);
	}

      if ( processNums() == 1 && ompNumThreads == 1 ) timer_start(timer_write);
#if defined(HAVE_LIBPTHREAD)
      if ( cdoLockIO )
	pthread_mutex_lock(&streamMutex);
      else
	pthread_mutex_lock(&streamOpenWriteMutex);
#endif
      int fileID = streamOpenWrite(argument->args, filetype);
#if defined(HAVE_LIBPTHREAD)
      if ( cdoLockIO )
	pthread_mutex_unlock(&streamMutex);
      else
	pthread_mutex_unlock(&streamOpenWriteMutex);
#endif
      if ( processNums() == 1 && ompNumThreads == 1 ) timer_stop(timer_write);
      if ( fileID < 0 ) cdiOpenError(fileID, "Open failed on >%s<", argument->args);

      cdoDefHistory(fileID, commandLine());

      if ( cdoDefaultByteorder != CDI_UNDEFID )
	streamDefByteorder(fileID, cdoDefaultByteorder);

      if ( cdoCompress )
	{
	  if      ( filetype == FILETYPE_GRB )
	    {
	      cdoCompType  = COMPRESS_SZIP;
	      cdoCompLevel = 0;
	    }
	  else if ( filetype == FILETYPE_NC4 || filetype == FILETYPE_NC4C )
	    {
	      cdoCompType  = COMPRESS_ZIP;
	      cdoCompLevel = 1;
	    }
	}

      if ( cdoCompType != COMPRESS_NONE )
	{
	  streamDefCompType(fileID, cdoCompType);
	  streamDefCompLevel(fileID, cdoCompLevel);

	  if ( cdoCompType == COMPRESS_SZIP &&
	       (filetype != FILETYPE_GRB && filetype != FILETYPE_GRB2 && filetype != FILETYPE_NC4 && filetype != FILETYPE_NC4C) )
	    cdoWarning("SZIP compression not available for non GRIB/netCDF4 data!");

	  if ( cdoCompType == COMPRESS_JPEG && filetype != FILETYPE_GRB2 )
	    cdoWarning("JPEG compression not available for non GRIB2 data!");

	  if ( cdoCompType == COMPRESS_ZIP && (filetype != FILETYPE_NC4 && filetype != FILETYPE_NC4C) )
	    cdoWarning("Deflate compression not available for non netCDF4 data!");
	}
      /*
      if ( cdoDefaultInstID != CDI_UNDEFID )
	streamDefInstID(fileID, cdoDefaultInstID);
      */
      strcpy(filename, argument->args);

      pstreamptr->mode     = 'w';
      pstreamptr->name     = filename;
      pstreamptr->fileID   = fileID;
      pstreamptr->filetype = filetype;
   }

  return (pstreamID);
}


int pstreamOpenAppend(const argument_t *argument)
{
  int pstreamID = -1;

  int ispipe = strncmp(argument->args, "(pipe", 5) == 0;

  if ( ispipe )
    {
      if ( PSTREAM_Debug ) Message("pipe %s", argument->args);
      cdoAbort("this operator doesn't work with pipes!");
    }
  else
    {
      char *filename = (char*) malloc(strlen(argument->args)+1);

      pstream_t *pstreamptr = pstream_new_entry();
      if ( ! pstreamptr ) Error("No memory");

      pstreamID = pstreamptr->self;
  
      if ( PSTREAM_Debug ) Message("file %s", argument->args);

      if ( processNums() == 1 && ompNumThreads == 1 ) timer_start(timer_write);
#if defined(HAVE_LIBPTHREAD)
      if ( cdoLockIO )
	pthread_mutex_lock(&streamMutex);
      else
	pthread_mutex_lock(&streamOpenReadMutex);
#endif
      int fileID = streamOpenAppend(argument->args);
#if defined(HAVE_LIBPTHREAD)
      if ( cdoLockIO )
	pthread_mutex_unlock(&streamMutex);
      else
	pthread_mutex_unlock(&streamOpenReadMutex);
#endif
      if ( processNums() == 1 && ompNumThreads == 1 ) timer_stop(timer_write);
      if ( fileID < 0 ) cdiOpenError(fileID, "Open failed on >%s<", argument->args);
      /*
      cdoInqHistory(fileID);
      cdoDefHistory(fileID, commandLine());
      */
      strcpy(filename, argument->args);

      pstreamptr->mode   = 'a';
      pstreamptr->name   = filename;
      pstreamptr->fileID = fileID;
    }

  return (pstreamID);
}


void pstreamClose(int pstreamID)
{
  pstream_t *pstreamptr = pstream_to_pointer(pstreamID);

  if ( pstreamptr == NULL )
    Error("Internal problem, stream %d not open!", pstreamID);

  if ( pstreamptr->ispipe )
    {
#if defined(HAVE_LIBPTHREAD)
      pipe_t *pipe;
      int lread = FALSE, lwrite = FALSE;
      pthread_t threadID = pthread_self();

      if      ( pthread_equal(threadID, pstreamptr->rthreadID) ) lread  = TRUE;
      else if ( pthread_equal(threadID, pstreamptr->wthreadID) ) lwrite = TRUE;
      else Error("Internal problem! Close pipe %s", pstreamptr->name);

      if ( lread )
	{
	  pipe = pstreamptr->pipe;
	  pthread_mutex_lock(pipe->mutex);
	  pipe->EOP = TRUE;
	  if ( PSTREAM_Debug ) Message("%s read closed", pstreamptr->name);
	  pthread_mutex_unlock(pipe->mutex);     
	  pthread_cond_signal(pipe->tsDef);
	  pthread_cond_signal(pipe->tsInq);
	 
	  pthread_cond_signal(pipe->recInq);
	 
	  pthread_mutex_lock(pipe->mutex);
	  pstreamptr->isopen = FALSE;
	  pthread_mutex_unlock(pipe->mutex);     
	  pthread_cond_signal(pipe->isclosed);

	  pthread_join(pstreamptr->wthreadID, NULL);

	  pthread_mutex_lock(pipe->mutex);
	  if ( pstreamptr->name ) free(pstreamptr->name);
	  if ( pstreamptr->argument )
	    {
	      argument_t *argument = (argument_t *) (pstreamptr->argument);
	      if ( argument->argv ) free(argument->argv);
	      if ( argument->args ) free(argument->args);
	      free(argument);
	    }
	  vlistDestroy(pstreamptr->vlistID);
	  pthread_mutex_unlock(pipe->mutex);

	  processAddNvals(pipe->nvals);
	  pipeDelete(pipe);

	  pstream_delete_entry(pstreamptr);
	}
      else if ( lwrite )
	{
	  pipe = pstreamptr->pipe;
	  pthread_mutex_lock(pipe->mutex);
	  pipe->EOP = TRUE;
	  if ( PSTREAM_Debug ) Message("%s write closed", pstreamptr->name);
	  pthread_mutex_unlock(pipe->mutex);     
	  pthread_cond_signal(pipe->tsDef);
	  pthread_cond_signal(pipe->tsInq);

	  pthread_mutex_lock(pipe->mutex);
	  while ( pstreamptr->isopen )
	    {
	      if ( PSTREAM_Debug ) Message("wait of read close");
	      pthread_cond_wait(pipe->isclosed, pipe->mutex);
	    }
	  pthread_mutex_unlock(pipe->mutex);
	}

      processDelStream(pstreamID);
#else
      cdoAbort("Cannot use pipes, pthread support not compiled in!");
#endif
    }
  else
    {
      if ( PSTREAM_Debug )
	Message("%s fileID %d\n", pstreamptr->name, pstreamptr->fileID);

      if ( pstreamptr->mode == 'r' )
	{
	  processAddNvals(streamNvals(pstreamptr->fileID));
	}

#if defined(HAVE_LIBPTHREAD)
      if ( cdoLockIO ) pthread_mutex_lock(&streamMutex);
#endif
      streamClose(pstreamptr->fileID);
#if defined(HAVE_LIBPTHREAD)
      if ( cdoLockIO ) pthread_mutex_unlock(&streamMutex);
#endif

      if ( cdoExpMode == CDO_EXP_REMOTE )
	{
	  if ( pstreamptr->mode == 'w' )
	    {
	      extern const char *cdojobfiles;
	      FILE *fp = fopen(cdojobfiles, "a");
	      fprintf(fp, "%s\n", pstreamptr->name);
	      fclose(fp);
	    }
	}

      if ( pstreamptr->name )
	{
	  free(pstreamptr->name);
	  pstreamptr->name = NULL;
	}

      if ( pstreamptr->varlist )
	{
	  free(pstreamptr->varlist);
	  pstreamptr->varlist = NULL;
	}

      pstream_delete_entry(pstreamptr);
    }
}


int pstreamInqVlist(int pstreamID)
{
  pstream_t *pstreamptr = pstream_to_pointer(pstreamID);

  int vlistID = -1;

#if defined(HAVE_LIBPTHREAD)
  if ( pstreamptr->ispipe )
    {
      vlistID = pipeInqVlist(pstreamptr);
      if ( vlistID == -1 )
	cdoAbort("Couldn't read data from input stream %s!", pstreamptr->name);
    }
  else
#endif
    {
      if ( processNums() == 1 && ompNumThreads == 1 ) timer_start(timer_read);
#if defined(HAVE_LIBPTHREAD)
      if ( cdoLockIO ) pthread_mutex_lock(&streamMutex);
#endif
      vlistID = streamInqVlist(pstreamptr->fileID);
#if defined(HAVE_LIBPTHREAD)
      if ( cdoLockIO ) pthread_mutex_unlock(&streamMutex);
#endif
      if ( processNums() == 1 && ompNumThreads == 1 ) timer_stop(timer_read);

      int nsubtypes = vlistNsubtypes(vlistID);
      if ( nsubtypes > 0 )
        cdoWarning("Subtypes are unsupported, the processing results are possibly wrong!");

      if ( cdoDefaultTimeType != CDI_UNDEFID )
	taxisDefType(vlistInqTaxis(vlistID), cdoDefaultTimeType);

      pstreamptr->vlistID = vlistID;
    }

  if ( vlistNumber(vlistID) == CDI_COMP && cdoStreamNumber() == CDI_REAL )
    cdoAbort("Complex fields are not supported by this operator!");

  if ( vlistNumber(vlistID) == CDI_REAL && cdoStreamNumber() == CDI_COMP )
    cdoAbort("This operator needs complex fields!");

  processDefVarNum(vlistNvars(vlistID), pstreamID);

  return vlistID;
}

static
const char *cdoComment(void)
{
  static char comment[256];
  static int init = 0;

  if ( ! init )
    {
      init = 1;

      int size = strlen(CDO_Version);

      strncat(comment, CDO_Version, size);
      comment[size] = 0;
    }

  return (comment);
}

static
void pstreamDefVarlist(pstream_t *pstreamptr, int vlistID)
{
  int filetype = pstreamptr->filetype;

  if ( pstreamptr->vlistID != -1 )
    cdoAbort("Internal problem, vlist already defined!");

  if ( pstreamptr->varlist != NULL )
    cdoAbort("Internal problem, varlist already allocated!");

  int nvars = vlistNvars(vlistID);
  varlist_t *varlist = (varlist_t*) malloc(nvars*sizeof(varlist_t));

  for ( int varID = 0; varID < nvars; ++varID )
    {
      varlist[varID].gridsize    = gridInqSize(vlistInqVarGrid(vlistID, varID));
      varlist[varID].datatype    = vlistInqVarDatatype(vlistID, varID);
      varlist[varID].missval     = vlistInqVarMissval(vlistID, varID);
      varlist[varID].addoffset   = vlistInqVarAddoffset(vlistID, varID);
      varlist[varID].scalefactor = vlistInqVarScalefactor(vlistID, varID);

      varlist[varID].check_datarange = FALSE;

      int laddoffset   = IS_NOT_EQUAL(varlist[varID].addoffset, 0);
      int lscalefactor = IS_NOT_EQUAL(varlist[varID].scalefactor, 1);

      int datatype = varlist[varID].datatype;

      if ( filetype == FILETYPE_NC || filetype == FILETYPE_NC2 || filetype == FILETYPE_NC4 || filetype == FILETYPE_NC4C )
	{
	  if ( datatype == DATATYPE_UINT8 && (filetype == FILETYPE_NC || filetype == FILETYPE_NC2) )
	    {
	      datatype = DATATYPE_INT16;
	      varlist[varID].datatype = datatype;
	    }

	  if ( datatype == DATATYPE_UINT16 && (filetype == FILETYPE_NC || filetype == FILETYPE_NC2) )
	    {
	      datatype = DATATYPE_INT32;
	      varlist[varID].datatype = datatype;
	    }

	  if ( laddoffset || lscalefactor )
	    {
	      if ( datatype == DATATYPE_INT8   ||
		   datatype == DATATYPE_UINT8  ||
		   datatype == DATATYPE_INT16  ||
		   datatype == DATATYPE_UINT16 )
		varlist[varID].check_datarange = TRUE;
	    }
	  else if ( cdoCheckDatarange )
	    {
	      varlist[varID].check_datarange = TRUE;
	    }
	}
    }

  pstreamptr->varlist = varlist;
  pstreamptr->vlistID = vlistID; /* used for -r/-a */
}


void pstreamDefVlist(int pstreamID, int vlistID)
{
  pstream_t *pstreamptr = pstream_to_pointer(pstreamID);

#if defined(HAVE_LIBPTHREAD)
  if ( pstreamptr->ispipe )
    {
      int vlistIDcp = vlistDuplicate(vlistID);
      /*    pipeDefVlist(pstreamptr, vlistID);*/
      pipeDefVlist(pstreamptr, vlistIDcp);
    }
  else
#endif
    {
      if ( cdoDefaultDataType != CDI_UNDEFID )
	{
	  int varID, nvars = vlistNvars(vlistID);

	  for ( varID = 0; varID < nvars; ++varID )
	    vlistDefVarDatatype(vlistID, varID, cdoDefaultDataType);

	  if ( cdoDefaultDataType == DATATYPE_FLT64 || cdoDefaultDataType == DATATYPE_FLT32 )
	    {
	      for ( varID = 0; varID < nvars; varID++ )
		{
		  vlistDefVarAddoffset(vlistID, varID, 0.0);
		  vlistDefVarScalefactor(vlistID, varID, 1.0);
		}
	    }
	}

      if ( cdoChunkType != CDI_UNDEFID )
	{
	  int varID, nvars = vlistNvars(vlistID);

	  for ( varID = 0; varID < nvars; ++varID )
	    vlistDefVarChunkType(vlistID, varID, cdoChunkType);
	}

      vlistDefAttTxt(vlistID, CDI_GLOBAL, "CDO", (int)strlen(cdoComment())+1, cdoComment());

#if defined(_OPENMP)
      if ( ompNumThreads > 1 )
	vlistDefAttInt(vlistID, CDI_GLOBAL, "cdo_openmp_thread_number", DATATYPE_INT32, 1, &ompNumThreads);
#endif
      pstreamDefVarlist(pstreamptr, vlistID);

      if ( processNums() == 1 && ompNumThreads == 1 ) timer_start(timer_write);
#if defined(HAVE_LIBPTHREAD)
      if ( cdoLockIO ) pthread_mutex_lock(&streamMutex);
#endif
      streamDefVlist(pstreamptr->fileID, vlistID);
#if defined(HAVE_LIBPTHREAD)
      if ( cdoLockIO ) pthread_mutex_unlock(&streamMutex);
#endif
      if ( processNums() == 1 && ompNumThreads == 1 ) timer_stop(timer_write);
    }
}


int pstreamInqRecord(int pstreamID, int *varID, int *levelID)
{
  pstream_t *pstreamptr;

  pstreamptr = pstream_to_pointer(pstreamID);

#if defined(HAVE_LIBPTHREAD)
  if ( pstreamptr->ispipe )
    pipeInqRecord(pstreamptr, varID, levelID);
  else
#endif
    {
      if ( processNums() == 1 && ompNumThreads == 1 ) timer_start(timer_read);
#if defined(HAVE_LIBPTHREAD)
      if ( cdoLockIO ) pthread_mutex_lock(&streamMutex);
#endif
      streamInqRecord(pstreamptr->fileID, varID, levelID);
#if defined(HAVE_LIBPTHREAD)
      if ( cdoLockIO ) pthread_mutex_unlock(&streamMutex);
#endif
      if ( processNums() == 1 && ompNumThreads == 1 ) timer_stop(timer_read);
    }

  return (0);
}


void pstreamDefRecord(int pstreamID, int varID, int levelID)
{
  pstream_t *pstreamptr;

  pstreamptr = pstream_to_pointer(pstreamID);
  
  pstreamptr->varID = varID;

#if defined(HAVE_LIBPTHREAD)
  if ( pstreamptr->ispipe )
    {
      pipeDefRecord(pstreamptr, varID, levelID);
    }
  else
#endif
    {
      if ( processNums() == 1 && ompNumThreads == 1 ) timer_start(timer_write);
#if defined(HAVE_LIBPTHREAD)
      if ( cdoLockIO ) pthread_mutex_lock(&streamMutex);
#endif
      streamDefRecord(pstreamptr->fileID, varID, levelID);
#if defined(HAVE_LIBPTHREAD)
      if ( cdoLockIO ) pthread_mutex_unlock(&streamMutex);
#endif
      if ( processNums() == 1 && ompNumThreads == 1 ) timer_stop(timer_write);
    }
}


void pstreamReadRecord(int pstreamID, double *data, int *nmiss)
{
  if ( data == NULL ) cdoAbort("Data pointer not allocated (pstreamReadRecord)!");

  pstream_t *pstreamptr = pstream_to_pointer(pstreamID);

#if defined(HAVE_LIBPTHREAD)
  if ( pstreamptr->ispipe )
    pipeReadRecord(pstreamptr, data, nmiss);
  else
#endif
    {
      if ( processNums() == 1 && ompNumThreads == 1 ) timer_start(timer_read);
#if defined(HAVE_LIBPTHREAD)
      if ( cdoLockIO ) pthread_mutex_lock(&streamMutex);
#endif
      streamReadRecord(pstreamptr->fileID, data, nmiss);
#if defined(HAVE_LIBPTHREAD)
      if ( cdoLockIO ) pthread_mutex_unlock(&streamMutex);
#endif
      if ( processNums() == 1 && ompNumThreads == 1 ) timer_stop(timer_read);
    }
}


void pstreamCheckDatarange(pstream_t *pstreamptr, int varID, double *array, int nmiss)
{
  long i;
  long   gridsize    = pstreamptr->varlist[varID].gridsize;
  int    datatype    = pstreamptr->varlist[varID].datatype;
  double missval     = pstreamptr->varlist[varID].missval;
  double addoffset   = pstreamptr->varlist[varID].addoffset;
  double scalefactor = pstreamptr->varlist[varID].scalefactor;

  long ivals   = 0;
  double arrmin  =  1.e300;
  double arrmax  = -1.e300;
  if ( nmiss > 0 )
    {
      for ( i = 0; i < gridsize; ++i )
	{
	  if ( !DBL_IS_EQUAL(array[i], missval) )
	    {
	      if ( array[i] < arrmin ) arrmin = array[i];
	      if ( array[i] > arrmax ) arrmax = array[i];
	      ivals++;
	    }
	}
    }
  else
    {
      for ( i = 0; i < gridsize; ++i )
	{
	  if ( array[i] < arrmin ) arrmin = array[i];
	  if ( array[i] > arrmax ) arrmax = array[i];
	}
      ivals = gridsize;
    }

  if ( ivals > 0 )
    {
      double smin = (arrmin - addoffset)/scalefactor;
      double smax = (arrmax - addoffset)/scalefactor;

      if ( datatype == DATATYPE_INT8  || datatype == DATATYPE_UINT8 ||
	   datatype == DATATYPE_INT16 || datatype == DATATYPE_UINT16 )
	{
	  smin = (int) round(smin);
	  smax = (int) round(smax);
	}

      double vmin = 0, vmax = 0;

      if      ( datatype == DATATYPE_INT8   ) { vmin =        -128.; vmax =        127.; }
      else if ( datatype == DATATYPE_UINT8  ) { vmin =           0.; vmax =        255.; }
      else if ( datatype == DATATYPE_INT16  ) { vmin =      -32768.; vmax =      32767.; }
      else if ( datatype == DATATYPE_UINT16 ) { vmin =           0.; vmax =      65535.; }
      else if ( datatype == DATATYPE_INT32  ) { vmin = -2147483648.; vmax = 2147483647.; }
      else if ( datatype == DATATYPE_UINT32 ) { vmin =           0.; vmax = 4294967295.; }
      else if ( datatype == DATATYPE_FLT32  ) { vmin = -3.40282e+38; vmax = 3.40282e+38; }
      else                                    { vmin =     -1.e+300; vmax =     1.e+300; }

      if ( smin < vmin || smax > vmax )
	cdoWarning("Some data values (min=%g max=%g) are outside the\n"
		   "    valid range (%g - %g) of the used output precision!\n"
		   "    Use the CDO option%s -b 64 to increase the output precision.",
		   smin, smax, vmin, vmax, (datatype == DATATYPE_FLT32) ? "" : " -b 32 or");
    }
}


void pstreamWriteRecord(int pstreamID, double *data, int nmiss)
{
  if ( data == NULL ) cdoAbort("Data pointer not allocated (%s)!", __func__);

  pstream_t *pstreamptr = pstream_to_pointer(pstreamID);

#if defined(HAVE_LIBPTHREAD)
  if ( pstreamptr->ispipe )
    {
      pipeWriteRecord(pstreamptr, data, nmiss);
    }
  else
#endif
    {
      int varID = pstreamptr->varID;
      if ( processNums() == 1 && ompNumThreads == 1 ) timer_start(timer_write);

      if ( pstreamptr->varlist )
	if ( pstreamptr->varlist[varID].check_datarange )
	  pstreamCheckDatarange(pstreamptr, varID, data, nmiss);

#if defined(HAVE_LIBPTHREAD)
      if ( cdoLockIO ) pthread_mutex_lock(&streamMutex);
#endif
      streamWriteRecord(pstreamptr->fileID, data, nmiss);
#if defined(HAVE_LIBPTHREAD)
      if ( cdoLockIO ) pthread_mutex_unlock(&streamMutex);
#endif

      if ( processNums() == 1 && ompNumThreads == 1 ) timer_stop(timer_write);
    }
}


void pstreamWriteRecordF(int pstreamID, float *data, int nmiss)
{
  if ( data == NULL ) cdoAbort("Data pointer not allocated (%s)!", __func__);

  pstream_t *pstreamptr = pstream_to_pointer(pstreamID);

#if defined(HAVE_LIBPTHREAD)
  if ( pstreamptr->ispipe )
    {
      cdoAbort("pipeWriteRecord not implemented for memtype float!");
      //pipeWriteRecord(pstreamptr, data, nmiss);
    }
  else
#endif
    {
      // int varID = pstreamptr->varID;
      if ( processNums() == 1 && ompNumThreads == 1 ) timer_start(timer_write);
      /*
      if ( pstreamptr->varlist )
	if ( pstreamptr->varlist[varID].check_datarange )
	  pstreamCheckDatarange(pstreamptr, varID, data, nmiss);
      */
#if defined(HAVE_LIBPTHREAD)
      if ( cdoLockIO ) pthread_mutex_lock(&streamMutex);
#endif
      streamWriteRecordF(pstreamptr->fileID, data, nmiss);
#if defined(HAVE_LIBPTHREAD)
      if ( cdoLockIO ) pthread_mutex_unlock(&streamMutex);
#endif
      if ( processNums() == 1 && ompNumThreads == 1 ) timer_stop(timer_write);
    }
}


int pstreamInqTimestep(int pstreamID, int tsID)
{
  pstream_t *pstreamptr = pstream_to_pointer(pstreamID);

  int nrecs = 0;

#if defined(HAVE_LIBPTHREAD)
  if ( pstreamptr->ispipe )
    nrecs = pipeInqTimestep(pstreamptr, tsID);
  else
#endif
    {
      if ( pstreamptr->mfiles ) tsID -= pstreamptr->tsID0;

      if ( processNums() == 1 && ompNumThreads == 1 ) timer_start(timer_read);
#if defined(HAVE_LIBPTHREAD)
      if ( cdoLockIO ) pthread_mutex_lock(&streamMutex);
#endif
      nrecs = streamInqTimestep(pstreamptr->fileID, tsID);
#if defined(HAVE_LIBPTHREAD)
      if ( cdoLockIO ) pthread_mutex_unlock(&streamMutex);
#endif
      if ( processNums() == 1 && ompNumThreads == 1 ) timer_stop(timer_read);

      if ( nrecs == 0 && pstreamptr->mfiles &&
	   (pstreamptr->nfiles < pstreamptr->mfiles) )
	{
	  size_t len;
	  int nfile = pstreamptr->nfiles;
	  char *filename = NULL;
	  int fileID;
	  int vlistIDold, vlistIDnew;

	  pstreamptr->tsID0 += tsID;

	  vlistIDold = vlistDuplicate(streamInqVlist(pstreamptr->fileID));
	  streamClose(pstreamptr->fileID);

	  len = strlen(pstreamptr->mfnames[nfile]);
	  filename = (char*) malloc(len+1);
	  strcpy(filename, pstreamptr->mfnames[nfile]);
	  pstreamptr->nfiles++;

#if defined(HAVE_LIBPTHREAD)
	  if ( cdoLockIO )
	    pthread_mutex_lock(&streamMutex);
	  else
	    pthread_mutex_lock(&streamOpenReadMutex);
#endif
	  if ( cdoVerbose ) cdoPrint("Continuation file: %s", filename);

	  if ( processNums() == 1 && ompNumThreads == 1 ) timer_start(timer_read);
	  fileID = streamOpenRead(filename);
	  vlistIDnew = streamInqVlist(fileID);
	  if ( processNums() == 1 && ompNumThreads == 1 ) timer_stop(timer_read);

	  vlistCompare(vlistIDold, vlistIDnew, CMP_HRD);
	  vlistDestroy(vlistIDold);
#if defined(HAVE_LIBPTHREAD)
	  if ( cdoLockIO )
	    pthread_mutex_unlock(&streamMutex);
	  else
	    pthread_mutex_unlock(&streamOpenReadMutex);
#endif
	  if ( fileID < 0 ) cdiOpenError(fileID, "Open failed on >%s<", filename);

	  free(pstreamptr->name);

	  pstreamptr->name   = filename;
	  pstreamptr->fileID = fileID;

	  if ( processNums() == 1 && ompNumThreads == 1 ) timer_start(timer_read);
#if defined(HAVE_LIBPTHREAD)
      if ( cdoLockIO ) pthread_mutex_lock(&streamMutex);
#endif
	  nrecs = streamInqTimestep(pstreamptr->fileID, 0);
#if defined(HAVE_LIBPTHREAD)
      if ( cdoLockIO ) pthread_mutex_unlock(&streamMutex);
#endif
	  if ( processNums() == 1 && ompNumThreads == 1 ) timer_stop(timer_read);
	}

      if ( tsID == 0 && cdoDefaultTimeType != CDI_UNDEFID )
	taxisDefType(vlistInqTaxis(pstreamptr->vlistID), cdoDefaultTimeType);
    }

  if ( nrecs && tsID != pstreamptr->tsID )
    {
      processDefTimesteps(pstreamID);
      pstreamptr->tsID = tsID;
    }

  return (nrecs);
}


void pstreamDefTimestep(int pstreamID, int tsID)
{
  pstream_t *pstreamptr = pstream_to_pointer(pstreamID);

#if defined(HAVE_LIBPTHREAD)
  if ( pstreamptr->ispipe )
    pipeDefTimestep(pstreamptr, tsID);
  else
#endif
    {
      if ( tsID == 0 && cdoDefaultTimeType != CDI_UNDEFID )
	{
	  int taxisID, vlistID;
	  vlistID = pstreamptr->vlistID;
	  taxisID = vlistInqTaxis(vlistID);
      	  taxisDefType(taxisID, cdoDefaultTimeType);
	}

      if ( processNums() == 1 && ompNumThreads == 1 ) timer_start(timer_write);
      /* don't use sync -> very slow on GPFS */
      //  if ( tsID > 0 ) streamSync(pstreamptr->fileID);
#if defined(HAVE_LIBPTHREAD)
      if ( cdoLockIO ) pthread_mutex_lock(&streamMutex);
#endif
      streamDefTimestep(pstreamptr->fileID, tsID);
#if defined(HAVE_LIBPTHREAD)
      if ( cdoLockIO ) pthread_mutex_unlock(&streamMutex);
#endif
      if ( processNums() == 1 && ompNumThreads == 1 ) timer_stop(timer_write);
    }
}


void pstreamCopyRecord(int pstreamIDdest, int pstreamIDsrc)
{
  if ( PSTREAM_Debug )
    Message("pstreamIDdest = %d  pstreamIDsrc = %d", pstreamIDdest, pstreamIDsrc);

  pstream_t *pstreamptr_dest = pstream_to_pointer(pstreamIDdest);
  pstream_t *pstreamptr_src  = pstream_to_pointer(pstreamIDsrc);

  if ( pstreamptr_dest->ispipe || pstreamptr_src->ispipe )
    cdoAbort("This operator can't be combined with other operators!");

  /*
#if defined(HAVE_LIBPTHREAD)
  if ( pstreamptr_dest->ispipe || pstreamptr_src->ispipe )
    {
      pipeCopyRecord(pstreamptr_dest, pstreamptr_src);
    }
  else
#endif
  */
    {
#if defined(HAVE_LIBPTHREAD)
      if ( cdoLockIO ) pthread_mutex_lock(&streamMutex);
#endif
      streamCopyRecord(pstreamptr_dest->fileID, pstreamptr_src->fileID);
#if defined(HAVE_LIBPTHREAD)
      if ( cdoLockIO ) pthread_mutex_unlock(&streamMutex);
#endif
    }
}


void pstreamDebug(int debug)
{
  PSTREAM_Debug = debug;
}


void cdoInitialize(void *argument)
{
#if defined(_OPENMP)
  omp_set_num_threads(ompNumThreads); /* Have to be called for every module (pthread)! */
#endif

  processCreate();

#if defined(HAVE_LIBPTHREAD)
  if ( PSTREAM_Debug )
     Message("process %d  thread %ld", processSelf(), pthread_self());
#endif

  processDefArgument(argument);
}


void pstreamCloseAll(void)
{
  if ( _pstreamList == NULL ) return;

  for ( int i = 0; i < _pstream_max; i++ )
    {
      pstream_t *pstreamptr = _pstreamList[i].ptr;
      if ( pstreamptr && pstreamptr->isopen )
	{
	  if ( !pstreamptr->ispipe )
	    {
	      if ( PSTREAM_Debug )
		Message("Close file %s id %d", pstreamptr->name, pstreamptr->fileID);
	      streamClose(pstreamptr->fileID);
	    }
	}
    }
}

static
void processClosePipes(void)
{
  int nstream = processInqStreamNum();
  for ( int sindex = 0; sindex < nstream; sindex++ )
    {
      int pstreamID = processInqStreamID(sindex);
      pstream_t *pstreamptr = pstream_to_pointer(pstreamID);

      if ( PSTREAM_Debug )
	Message("process %d  stream %d  close streamID %d", processSelf(), sindex, pstreamID);

      if ( pstreamptr ) pstreamClose(pstreamID);
    }
}


void cdoFinish(void)
{
  int processID = processSelf();
  int nvars, ntimesteps;
  char memstring[32] = {""};
  double s_utime, s_stime;
  double e_utime, e_stime;
  double c_cputime = 0, c_usertime = 0, c_systime = 0;
  double p_cputime = 0, p_usertime = 0, p_systime = 0;

#if defined(HAVE_LIBPTHREAD)
  if ( PSTREAM_Debug )
    Message("process %d  thread %ld", processID, pthread_self());
#endif

  int64_t nvals = processInqNvals(processID);
  nvars = processInqVarNum();
  ntimesteps = processInqTimesteps();

  if ( !cdoSilentMode )
    {
      set_text_color(stderr, RESET, GREEN);
      fprintf(stderr, "%s: ", processInqPrompt());
      reset_text_color(stderr);
      if ( nvals > 0 )
	{
	  if ( sizeof(int64_t) > sizeof(long) )
#if defined(_WIN32)
	    fprintf(stderr, "Processed %I64d value%s from %d variable%s",
#else
	    fprintf(stderr, "Processed %jd value%s from %d variable%s",
#endif
		    (intmax_t) nvals, ADD_PLURAL(nvals), nvars, ADD_PLURAL(nvars));
	  else
	    fprintf(stderr, "Processed %ld value%s from %d variable%s",
		    (long) nvals, ADD_PLURAL(nvals), nvars, ADD_PLURAL(nvars));
	}
      else if ( nvars > 0 )
	{
	  fprintf(stderr, "Processed %d variable%s", nvars,  ADD_PLURAL(nvars));
	}

      if ( ntimesteps > 0 )
	fprintf(stderr, " over %d timestep%s", ntimesteps,  ADD_PLURAL(ntimesteps));

      //  fprintf(stderr, ".");
    }
  /*
    fprintf(stderr, "%s: Processed %d variable%s %d timestep%s.",
	    processInqPrompt(), nvars, nvars > 1 ? "s" : "",
	    ntimesteps, ntimesteps > 1 ? "s" : "");
  */
  processStartTime(&s_utime, &s_stime);
  cdoProcessTime(&e_utime, &e_stime);

  c_usertime = e_utime - s_utime;
  c_systime  = e_stime - s_stime;
  c_cputime  = c_usertime + c_systime;

#if defined(HAVE_LIBPTHREAD)
  if ( pthreadScope == PTHREAD_SCOPE_PROCESS )
    {
      c_usertime /= processNums();
      c_systime  /= processNums();
      c_cputime  /= processNums();
    }
#endif

  processDefCputime(processID, c_cputime);  

  processAccuTime(c_usertime, c_systime);

  if ( processID == 0 )
    {
      int mu[] = {'b', 'k', 'm', 'g', 't'};
      int muindex = 0;
      long memmax;

      memmax = memTotal();
      while ( memmax > 9999 )
	{
	  memmax /= 1024;
	  muindex++;
	}

      if ( memmax )
	sprintf(memstring, " %ld%c ", memmax, mu[muindex]);

      processEndTime(&p_usertime, &p_systime);
      p_cputime  = p_usertime + p_systime;

      if ( cdoLogOff == 0 )
	{
	  cdologs(processNums()); 
	  cdologo(processNums()); 
	  cdolog(processInqPrompt(), p_cputime); 
	}
    }

#if defined(HAVE_SYS_TIMES_H)
  if ( cdoBenchmark )
    fprintf(stderr, " ( %.2fs %.2fs %.2fs %s)\n", c_usertime, c_systime, c_cputime, memstring);
  else
    {
      if ( ! cdoSilentMode )
	fprintf(stderr, " ( %.2fs )\n", c_cputime);
    }

  if ( cdoBenchmark && processID == 0 )
    fprintf(stderr, "total: user %.2fs  sys %.2fs  cpu %.2fs  mem%s\n",
	    p_usertime, p_systime, p_cputime, memstring);
#else
  fprintf(stderr, "\n");
#endif

  processClosePipes();

  processDelete();
}


int pstreamInqFiletype(int pstreamID)
{
  pstream_t *pstreamptr = pstream_to_pointer(pstreamID);

  int filetype;
  
#if defined(HAVE_LIBPTHREAD)
  if ( pstreamptr->ispipe )
    filetype = pstreamptr->filetype;
  else
#endif
    filetype = streamInqFiletype(pstreamptr->fileID);

  return (filetype);
}


int pstreamInqByteorder(int pstreamID)
{
  pstream_t *pstreamptr = pstream_to_pointer(pstreamID);

  int byteorder;

#if defined(HAVE_LIBPTHREAD)
  if ( pstreamptr->ispipe )
    byteorder = pstreamptr->filetype;
  else
#endif
    byteorder = streamInqByteorder(pstreamptr->fileID);

  return (byteorder);
}

void pstreamInqGRIBinfo(int pstreamID, int *intnum, float *fltnum, off_t *bignum)
{
  pstream_t *pstreamptr = pstream_to_pointer(pstreamID);

  streamInqGRIBinfo(pstreamptr->fileID, intnum, fltnum, bignum);
}


void cdoVlistCopyFlag(int vlistID2, int vlistID1)
{
#if defined(HAVE_LIBPTHREAD)
  pthread_mutex_lock(&streamMutex);
#endif

  vlistCopyFlag(vlistID2, vlistID1);

#if defined(HAVE_LIBPTHREAD)
  pthread_mutex_unlock(&streamMutex);
#endif

}


void openLock(void)
{
#if defined(HAVE_LIBPTHREAD)
  if ( cdoLockIO )
    pthread_mutex_lock(&streamMutex);
  else
    pthread_mutex_lock(&streamOpenReadMutex);
#endif  
}


void openUnlock(void)
{
#if defined(HAVE_LIBPTHREAD)
  if ( cdoLockIO )
    pthread_mutex_unlock(&streamMutex);
  else
    pthread_mutex_unlock(&streamOpenReadMutex);
#endif
}
