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

#if defined(HAVE_PTHREAD_H)
#  include <pthread.h>
#endif

#include <stdio.h>
#include <string.h>

#if defined(HAVE_GLOB_H)
#include <glob.h>
#endif
#if defined(HAVE_WORDEXP_H)
#include <wordexp.h>
#endif

#include "cdo.h"
#include "cdo_int.h"
#include "error.h"
#include "modules.h"
#include "util.h"
#include "pstream_int.h"
#include "dmemory.h"


#define  MAX_PROCESS    128
#define  MAX_STREAM      64
#define  MAX_OPERATOR   128
#define  MAX_OARGC     4096
#define  MAX_FILES    65536


typedef struct {
  int         f1;
  int         f2;
  const char *name;
  const char *enter;
}
oper_t;

typedef struct {
#if defined(HAVE_LIBPTHREAD)
  pthread_t   threadID;
  int         l_threadID;
#endif
  short       nchild;
  short       nstream;
  short       streams[MAX_STREAM];
  double      s_utime;
  double      s_stime;
  double      a_utime;
  double      a_stime;
  double      cputime;

  off_t       nvals;
  short       nvars;
  int         ntimesteps;
  short       streamCnt;
  argument_t *streamNames;
  char       *xoperator;
  char       *operatorName;
  char       *operatorArg;
  int         oargc;
  char       *oargv[MAX_OARGC];
  char        prompt[64];
  short       noper;
  oper_t      oper[MAX_OPERATOR];
}
process_t;


static process_t Process[MAX_PROCESS];

static int NumProcess = 0;
static int NumProcessActive = 0;

#if defined(HAVE_LIBPTHREAD)
pthread_mutex_t processMutex = PTHREAD_MUTEX_INITIALIZER;
#endif


int processCreate(void)
{
  int processID;

#if defined(HAVE_LIBPTHREAD)
  pthread_mutex_lock(&processMutex);
#endif

  processID = NumProcess++;
  NumProcessActive++;

#if defined(HAVE_LIBPTHREAD)
  pthread_mutex_unlock(&processMutex);  
#endif

  if ( processID >= MAX_PROCESS )
    Error("Limit of %d processes reached!", MAX_PROCESS);

#if defined(HAVE_LIBPTHREAD)
  Process[processID].threadID     = pthread_self();
  Process[processID].l_threadID   = 1;
#endif
  Process[processID].nstream      = 0;
  Process[processID].nchild       = 0;

  cdoProcessTime(&Process[processID].s_utime, &Process[processID].s_stime);
  Process[processID].a_utime      = 0;
  Process[processID].a_stime      = 0;
  Process[processID].cputime      = 0;

  Process[processID].oargc        = 0;
  Process[processID].xoperator    = NULL;
  Process[processID].operatorName = NULL;
  Process[processID].operatorArg  = NULL;

  Process[processID].noper        = 0;

  return (processID);
}


int processSelf(void)
{
  int processID = 0;
#if defined(HAVE_LIBPTHREAD)
  pthread_t thID = pthread_self();

  pthread_mutex_lock(&processMutex);

  for ( processID = 0; processID < NumProcess; processID++ )
    if ( Process[processID].l_threadID )
      if ( pthread_equal(Process[processID].threadID, thID) ) break;

  if ( processID == NumProcess )
    {
      if ( NumProcess > 0 )
        Error("Internal problem, process not found!");
      else
        processID = 0;
    }

  pthread_mutex_unlock(&processMutex);  

#endif

  return (processID);
}


int processNums(void)
{
#if defined(HAVE_LIBPTHREAD)
  pthread_mutex_lock(&processMutex);
#endif

  int pnums = NumProcess;

#if defined(HAVE_LIBPTHREAD)
  pthread_mutex_unlock(&processMutex);  
#endif

  return (pnums);
}


int processNumsActive(void)
{
#if defined(HAVE_LIBPTHREAD)
  pthread_mutex_lock(&processMutex);
#endif

  int pnums = NumProcessActive;

#if defined(HAVE_LIBPTHREAD)
  pthread_mutex_unlock(&processMutex);  
#endif

  return (pnums);
}


void processAddNvals(off_t nvals)
{
  int processID = processSelf();

  Process[processID].nvals += nvals;
}


off_t processInqNvals(int processID)
{
  return (Process[processID].nvals);
}


void processAddStream(int streamID)
{
  int processID = processSelf();
  int sindex;

  if ( pstreamIsPipe(streamID) ) Process[processID].nchild++;
  sindex = Process[processID].nstream++;

  if ( sindex >= MAX_STREAM )
    Error("limit of %d streams per process reached (processID = %d)!", MAX_STREAM, processID);

  Process[processID].streams[sindex] = streamID;
}


void processDelStream(int streamID)
{
}


void processDefCputime(int processID, double cputime)
{
  Process[processID].cputime = cputime;
}


double processInqCputime(int processID)
{
  return (Process[processID].cputime);
}


void processStartTime(double *utime, double *stime)
{
  int processID = processSelf();

  *utime = Process[processID].s_utime;
  *stime = Process[processID].s_stime;
}


void processEndTime(double *utime, double *stime)
{
  *utime = Process[0].a_utime;
  *stime = Process[0].a_stime;
}


void processAccuTime(double utime, double stime)
{
  Process[0].a_utime += utime;
  Process[0].a_stime += stime;
}


int processInqStreamNum(void)
{
  int processID = processSelf();

  return (Process[processID].nstream);
}


int processInqChildNum(void)
{
  int processID = processSelf();

  return (Process[processID].nchild);
}


int processInqStreamID(int streamindex)
{
  int processID = processSelf();

  return (Process[processID].streams[streamindex]);
}


const char *processInqOpername2(int processID)
{
  return (Process[processID].operatorName);
}


const char *processInqOpername(void)
{
  int processID = processSelf();

  return (Process[processID].operatorName);
}


void processDefPrompt(char *opername)
{
  int processID = processSelf();

  if ( processID == 0 )
    sprintf(Process[processID].prompt, "%s %s", Progname, opername);
  else
    sprintf(Process[processID].prompt, "%s(%d) %s", Progname, processID+1, opername);
}


const char *processInqPrompt(void)
{
  int processID = processSelf();

  return (Process[processID].prompt);
}

#if defined(HAVE_GLOB_H)
static
int get_glob_flags(void)
{
  int glob_flags = 0;

#if defined (GLOB_NOCHECK)
  glob_flags |= GLOB_NOCHECK;
#endif
#if defined (GLOB_TILDE)
  glob_flags |= GLOB_TILDE;
#endif

  return (glob_flags);
}
#endif

#if defined(HAVE_WORDEXP_H)
/* Convert a shell pattern into a list of filenames. */
static
argument_t *glob_pattern(const char *restrict string)
{
  size_t cnt, length = 0;
  int flags = WRDE_UNDEF;
  char **p;

  wordexp_t glob_results;
  argument_t *argument = NULL;

  // glob the input argument or do even more shell magic
  wordexp(string, &glob_results, flags);

  // How much space do we need?
  for ( p = glob_results.we_wordv, cnt = glob_results.we_wordc; cnt; p++, cnt-- )
    {
      length += strlen(*p) + 1;
    }

  // Allocate the space and generate the list.
  argument = argument_new(glob_results.we_wordc, length);

  // put all generated filenames into the argument_t data structure
  for ( cnt = 0; cnt < glob_results.we_wordc; cnt++ )
    {
      argument->argv[cnt] = strdupx(glob_results.we_wordv[cnt]);
      strcat(argument->args, glob_results.we_wordv[cnt]);
      if ( cnt < glob_results.we_wordc-1 ) strcat(argument->args, " ");
    }

  wordfree(&glob_results);

  return argument;
}
#endif

int cdoStreamCnt(void)
{
  int processID = processSelf();
  int cnt;

  cnt = Process[processID].streamCnt;

  return (cnt);
}


const argument_t *cdoStreamName(int cnt)
{
  int processID = processSelf();

  if ( cnt > Process[processID].streamCnt || cnt < 0 )
    Error("count %d out of range!", cnt);

  return (&(Process[processID].streamNames[cnt]));
}


const char *processOperator(void)
{
  int processID = processSelf();

  return (Process[processID].xoperator);
}

static
char *getOperatorArg(const char *xoperator)
{
  char *commapos;
  char *operatorArg = NULL;
  size_t len;

  if ( xoperator )
    {
      commapos = strchr(xoperator, ',');

      if ( commapos )
        {
          len = strlen(commapos+1);
          if ( len )
            {
              operatorArg = (char*) malloc(len+1);
              strcpy(operatorArg, commapos+1);
            }
        }
    }

  return (operatorArg);
}

static int skipInputStreams(int argc, char *argv[], int globArgc, int nstreams);

static
int getGlobArgc(int argc, char *argv[], int globArgc)
{
  int streamInCnt;
  int streamOutCnt;
  char *opername;
  char *comma_position;

  opername = &argv[globArgc][1];
  comma_position = strchr(opername, ',');
  if ( comma_position ) *comma_position = 0;

  streamInCnt  = operatorStreamInCnt(opername);
  streamOutCnt = operatorStreamOutCnt(opername);

  if ( streamInCnt == -1 ) streamInCnt = 1;

  if ( streamOutCnt > 1 )
    cdoAbort("More than one output stream not allowed in CDO pipes (Operator %s)!", opername);

  globArgc++;

  if ( streamInCnt > 0 )
    globArgc = skipInputStreams(argc, argv, globArgc, streamInCnt);
  if ( comma_position ) *comma_position = ',';

  return (globArgc);
}

static
int skipInputStreams(int argc, char *argv[], int globArgc, int nstreams)
{
  while ( nstreams > 0 )
    {
      if ( globArgc >= argc )
        {
          cdoAbort("Too few arguments. Check command line!");
          break;
        }
      if ( argv[globArgc][0] == '-' )
        {
          globArgc = getGlobArgc(argc, argv, globArgc);
        }
      else
        globArgc++;

      nstreams--;
    }

  return (globArgc);
}

static
int getStreamCnt(int argc, char *argv[])
{
  int streamCnt = 0;
  int globArgc = 1;

  while ( globArgc < argc )
    {
      if ( argv[globArgc][0] == '-' )
        {
          globArgc = getGlobArgc(argc, argv, globArgc);
        }
      else
        globArgc++;

      streamCnt++;
    }

  return (streamCnt);
}

static
void setStreamNames(int argc, char *argv[])
{
  int processID = processSelf();
  int i, ac;
  int globArgc = 1;
  int globArgcStart;
  char *streamname;
  int len;

  while ( globArgc < argc )
    {
      if ( argv[globArgc][0] == '-' )
        {
          globArgcStart = globArgc;

          globArgc = getGlobArgc(argc, argv, globArgc);
          len = 0;
          for ( i = globArgcStart; i < globArgc; i++ ) len += strlen(argv[i]) + 1;
          streamname = (char*) calloc(1, len);
          for ( i = globArgcStart; i < globArgc; i++ )
            {
              strcat(streamname, argv[i]);
              if ( i < globArgc-1 ) strcat(streamname, " ");
            }
          for ( i = 1; i < len-1; i++ ) if ( streamname[i] == '\0' ) streamname[i] = ' ';
          Process[processID].streamNames[Process[processID].streamCnt].args = streamname;
          ac = globArgc - globArgcStart;
          //printf("setStreamNames:  ac %d  streamname1: %s\n", ac, streamname);
          Process[processID].streamNames[Process[processID].streamCnt].argv = (char **) malloc(ac*sizeof(char *));
          for ( i = 0; i < ac; ++i )
            Process[processID].streamNames[Process[processID].streamCnt].argv[i] = argv[i+globArgcStart];
          Process[processID].streamNames[Process[processID].streamCnt].argc = ac;
          Process[processID].streamCnt++;
          //printf("setStreamNames:  streamname1: %s\n", streamname);
        }
      else
        {
          len = strlen(argv[globArgc]) + 1;
          streamname = (char*) malloc(len);
          strcpy(streamname, argv[globArgc]);
          Process[processID].streamNames[Process[processID].streamCnt].args = streamname;
          ac = 1;
          Process[processID].streamNames[Process[processID].streamCnt].argv = (char **) malloc(ac*sizeof(char *));
          Process[processID].streamNames[Process[processID].streamCnt].argv[0] = argv[globArgc];
          Process[processID].streamNames[Process[processID].streamCnt].argc = ac;
          Process[processID].streamNames[Process[processID].streamCnt].args = streamname;
          Process[processID].streamCnt++;
          //printf("setStreamNames:  streamname2: %s\n", streamname);
          globArgc++;
        }
    }
}

static
int find_wildcard(const char *string, size_t len)
{
  int status = 0;

  if ( len > 0 )
    {
      if ( string[0] == '~' ) status = 1;

      if ( status == 0 )
        {
          for ( size_t i = 0; i < len; ++i )
            if ( string[i] == '?' || string[i] == '*' || string[i] == '[' )
              {
                status = 1;
                break;
              }
        }
    }

  return status;
}


char *expand_filename(const char *string)
{
  char *filename = NULL;

  if ( find_wildcard(string, strlen(string)) )
    {
#if defined(HAVE_GLOB_H)
      int glob_flags = get_glob_flags();
      glob_t glob_results;

      glob(string, glob_flags, 0, &glob_results);

      if ( glob_results.gl_pathc == 1 ) filename = strdupx(glob_results.gl_pathv[0]);

      globfree(&glob_results);
#endif
    }

  return filename;
}

static
int expand_wildcards(int processID, int streamCnt)
{
  const char *streamname0 = Process[processID].streamNames[0].args;

  if ( streamname0[0] == '-' ) return 1;

#if defined(HAVE_WORDEXP_H)
  argument_t *glob_arg = glob_pattern(streamname0);

  // skip if the input argument starts with an operator (starts with -)
  // otherwise adapt streams if there are several files (>1)
  // in case of one filename skip, no adaption needed
  if ( glob_arg->argc > 1 && glob_arg->argv[0][0] != '-' )
    {
      int i;
      streamCnt = streamCnt - 1 + glob_arg->argc;

      free(Process[processID].streamNames[0].argv);
      free(Process[processID].streamNames[0].args);

      Process[processID].streamNames = (argument_t*) realloc(Process[processID].streamNames, streamCnt*sizeof(argument_t));
          
      // move output streams to the end
      for ( i = 1; i < Process[processID].streamCnt; ++i )
        Process[processID].streamNames[i+glob_arg->argc-1] = Process[processID].streamNames[i];

      for ( i = 0; i < glob_arg->argc; ++i )
        {
          Process[processID].streamNames[i].argv    = (char **) malloc(sizeof(char *));
          Process[processID].streamNames[i].argc    = 1;
          Process[processID].streamNames[i].argv[0] = strdupx(glob_arg->argv[i]);
          Process[processID].streamNames[i].args    = strdupx(glob_arg->argv[i]);
        }
      
      Process[processID].streamCnt = streamCnt;
    }

  free(glob_arg);
#endif

  return 1;
}


static
int checkStreamCnt(void)
{
  int processID = processSelf();
  int streamInCnt, streamOutCnt;
  int streamInCnt0;
  int streamCnt = 0;
  int i, j;
  int obase = FALSE;
  int status = 0;

  streamInCnt  = operatorStreamInCnt(Process[processID].operatorName);
  streamOutCnt = operatorStreamOutCnt(Process[processID].operatorName);

  streamInCnt0 = streamInCnt;

  if ( streamOutCnt == -1 )
    {
      streamOutCnt = 1;
      obase = TRUE;
    }

  if ( streamInCnt == -1 && streamOutCnt == -1 )
    cdoAbort("I/O stream counts unlimited no allowed!");
    
  // printf(" streamInCnt, streamOutCnt %d %d\n", streamInCnt, streamOutCnt);
  if ( streamInCnt == -1 )
    {
      streamInCnt = Process[processID].streamCnt - streamOutCnt;
      if ( streamInCnt < 1 ) cdoAbort("Input streams missing!");
    }

  if ( streamOutCnt == -1 )
    {
      streamOutCnt = Process[processID].streamCnt - streamInCnt;
      if ( streamOutCnt < 1 ) cdoAbort("Output streams missing!");
    }
  // printf(" streamInCnt, streamOutCnt %d %d\n", streamInCnt, streamOutCnt);

  streamCnt = streamInCnt + streamOutCnt;
  // printf(" streamCnt %d %d\n", Process[processID].streamCnt, streamCnt);

  if ( Process[processID].streamCnt > streamCnt )
    cdoAbort("Too many streams!"
             " Operator needs %d input and %d output streams.", streamInCnt, streamOutCnt);

  if ( Process[processID].streamCnt < streamCnt )
    cdoAbort("Too few streams specified!"
             " Operator needs %d input and %d output streams.", streamInCnt, streamOutCnt);

  for ( i = streamInCnt; i < streamCnt; i++ )
    {
      if ( Process[processID].streamNames[i].args[0] == '-' )
        {
          cdoAbort("Output file name %s must not begin with \"-\"!\n",
                   Process[processID].streamNames[i].args);
        }
      else if ( !obase )
        {
          for ( j = 0; j < streamInCnt; j++ ) /* does not work with files in pipes */
            if ( strcmp(Process[processID].streamNames[i].args, Process[processID].streamNames[j].args) == 0 )
              cdoAbort("Output file name %s is equal to input file name"
                       " on position %d!\n", Process[processID].streamNames[i].args, j+1);
        }
    }

  if ( streamInCnt == 1 && streamInCnt0 == -1 )
    status = expand_wildcards(processID, streamCnt);

  return (status);
}

static
void setStreams(int argc, char *argv[])
{
  int processID = processSelf();
  int streamCnt;
  int status;
  int i;

  streamCnt = getStreamCnt(argc, argv);

  Process[processID].nvals = 0;
  Process[processID].nvars = 0;
  Process[processID].ntimesteps = 0;

  Process[processID].streamCnt  = 0; /* filled in setStreamNames */
  if ( streamCnt )
    Process[processID].streamNames = (argument_t*) malloc(streamCnt*sizeof(argument_t));
  for ( i = 0; i < streamCnt; i++ )
    {
      Process[processID].streamNames[i].argc = 0;
      Process[processID].streamNames[i].argv = NULL;
      Process[processID].streamNames[i].args = NULL;
    }

  setStreamNames(argc, argv);

  status = checkStreamCnt();

  if ( status == 0 && Process[processID].streamCnt != streamCnt )
    Error("Internal problem with stream count %d %d", Process[processID].streamCnt, streamCnt);
  /*
  for ( i = 0; i < streamCnt; i++ )
    fprintf(stderr, "setStreams: stream %d %s\n", i+1, Process[processID].streamNames[i].args);
  */
}


void processDefArgument(void *vargument)
{
  int processID = processSelf();
  char *operatorArg;
  char *commapos;
  int oargc = 0;
  char **oargv = Process[processID].oargv;
  int argc = ((argument_t *) vargument)->argc;
  char **argv = ((argument_t *) vargument)->argv;

  Process[processID].xoperator    = argv[0];
  Process[processID].operatorName = getOperatorName(Process[processID].xoperator);
  Process[processID].operatorArg  = getOperatorArg(Process[processID].xoperator);
  operatorArg = Process[processID].operatorArg;

  if ( operatorArg )
    {
      oargv[oargc++] = operatorArg;
      //fprintf(stderr, "processDefArgument: %d %s\n", oargc, operatorArg);

      commapos = operatorArg;
      while ( (commapos = strchr(commapos, ',')) != NULL )
        {
          *commapos++ = '\0';
          if ( strlen(commapos) )
            {
              if ( oargc >= MAX_OARGC )
                cdoAbort("Too many parameter (limit=%d)!", MAX_OARGC);

              oargv[oargc++] = commapos;
            }
        }
      Process[processID].oargc = oargc;
    }

  processDefPrompt(Process[processID].operatorName);

  setStreams(argc, argv);
}


void processDefVarNum(int nvars, int streamID)
{
  int processID = processSelf();

  /*  if ( streamID == Process[processID].streams[0] ) */
    Process[processID].nvars += nvars;
}


int processInqVarNum(void)
{
  int processID = processSelf();

  return (Process[processID].nvars);
}


void processDefTimesteps(int streamID)
{
  int processID = processSelf();
  /*
  int i;
  printf("streamID %d %d %d %d\n", streamID, Process[processID].streams[0], Process[processID].streams[1], processID);

  for ( i = 0; i < Process[processID].nstream; i++)
    printf("streamID %d %d %d %d << \n", processID, Process[processID].nstream, i, Process[processID].streams[i]);
  */
  /*  if ( streamID == Process[processID].streams[0] )*/
    Process[processID].ntimesteps++;
}


int processInqTimesteps(void)
{
  int processID = processSelf();

  return (Process[processID].ntimesteps);
}


void processDelete(void)
{
  int processID = processSelf();

  //fprintf(stderr, "delete processID %d\n", processID);
#if defined(HAVE_LIBPTHREAD)
  pthread_mutex_lock(&processMutex);

  Process[processID].l_threadID = 0;
#endif
  NumProcessActive--;

#if defined(HAVE_LIBPTHREAD)
  pthread_mutex_unlock(&processMutex);  
#endif
}


int operatorArgc(void)
{
  int processID = processSelf();

  return (Process[processID].oargc);
}


char **operatorArgv(void)
{
  int processID = processSelf();

  return (Process[processID].oargv);
}


void operatorCheckArgc(int numargs)
{
  int processID = processSelf();
  int argc = Process[processID].oargc;

  if ( argc < numargs )
    cdoAbort("Too few arguments! Need %d found %d.", numargs, argc);
  else if ( argc > numargs )
    cdoAbort("Too many arguments! Need %d found %d.", numargs, argc);
}


void operatorInputArg(const char *enter)
{
  int processID = processSelf();
  int oargc;

  oargc = Process[processID].oargc;

  if ( oargc == 0 )
    {
      char line[1024];
      char *pline = line;
      size_t pos, len, linelen;
      int lreadline = 1;

      if ( enter )
        {
          set_text_color(stderr, BRIGHT, MAGENTA);
          fprintf(stderr, "%-16s : ", processInqPrompt());
          reset_text_color(stderr);
          // set_text_color(stderr, BLINK, BLACK);
          fprintf(stderr, "Enter %s > ", enter);
          // reset_text_color(stderr);
        }

      while ( lreadline )
        {
          readline(stdin, pline, 1024);

          lreadline = 0;
          while ( 1 )
            {
              pos = 0;
              while ( pline[pos] == ' ' || pline[pos] == ',' ) pos++;
              pline += pos;
              linelen = strlen(pline);
              if ( linelen > 0 )
                {
                  if ( pline[0] == '\\' )
                    {
                      lreadline = 1;
                      break;
                    }
                  len = 0;
                  while ( pline[len] != ' '  && pline[len] != ',' &&
                          pline[len] != '\\' && len < linelen ) len++;

                  Process[processID].oargv[oargc] = (char*) malloc(len+1);
                  memcpy(Process[processID].oargv[oargc], pline, len);
                  Process[processID].oargv[oargc][len] = '\0';
                  oargc++;

                  pline += len;
                }
              else
                break;
            }
        }

      Process[processID].oargc = oargc;
    }
}


int cdoOperatorAdd(const char *name, int f1, int f2, const char *enter)
{
  int processID = processSelf();
  int operID = Process[processID].noper;

  if ( operID < MAX_OPERATOR )
    {
      Process[processID].oper[operID].f1     = f1;
      Process[processID].oper[operID].f2     = f2;
      Process[processID].oper[operID].name   = name;
      Process[processID].oper[operID].enter  = enter;

      Process[processID].noper++;
    }
  else
    {
      cdoAbort("Maximum of %d operators reached!", MAX_OPERATOR);
    }

  return (operID);
}


int cdoOperatorID(void)
{
  int processID = processSelf();
  int operID = -1;

  if ( Process[processID].noper > 0 )
    {
      for ( operID = 0; operID < Process[processID].noper; operID++ )
        if ( Process[processID].oper[operID].name )
          if ( strcmp(Process[processID].operatorName, Process[processID].oper[operID].name) == 0 ) break;

      if ( operID == Process[processID].noper )
        cdoAbort("Operator not callable by this name!");
    }
  else
    {
      cdoAbort("Operator not initialized!");
    }

  return (operID);
}


int cdoOperatorF1(int operID)
{
  int processID = processSelf();

  return (Process[processID].oper[operID].f1);
}


int cdoOperatorF2(int operID)
{
  int processID = processSelf();

  return (Process[processID].oper[operID].f2);
}


const char *cdoOperatorName(int operID)
{
  int processID = processSelf();

  return (Process[processID].oper[operID].name);
}


const char *cdoOperatorEnter(int operID)
{
  int processID = processSelf();

  return (Process[processID].oper[operID].enter);
}


int cdoStreamNumber()
{
  int processID = processSelf();

  return (operatorStreamNumber(Process[processID].operatorName));
}
