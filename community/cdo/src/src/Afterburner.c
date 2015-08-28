#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif

#include <cdi.h>

#if defined(CDO)
#include "cdo.h"
#include "cdo_int.h"
#include "pstream_write.h"
#endif

#if defined(AFTERBURNER)
#include "afterdoc.h"
#endif

#include "afterburner.h"
#include "constants.h"
#include "compare.h"
#include "vct_l191.h"

#if  defined  (HAVE_LIBPTHREAD)
#include <pthread.h>
#endif

#if defined (_OPENMP)
#include <omp.h>
#endif

#if ! defined (VERSION)
#define  VERSION  "0.0.1"
#endif

#ifndef CLOCKS_PER_SEC
#define CLOCKS_PER_SEC 1000000
#endif

#if defined(AFTERBURNER)
static double starttime = 0.0;
#endif

void afterInqHistory(int fileID);
void afterDefHistory(int fileID, char *histstring);

int   scan_par_obsolate(char *namelist, char *name, int def);
void  scan_code(char *namelist, struct Variable *vars, int maxCodes, int *numCodes);
int   scan_par(int verbose, char *namelist, char *name, int def);
int   scan_time(int verbose, char *namelist, int *hours, int max_hours);
void  scan_darray(char *namelist, char *name, double *values, int maxValues, int *numValues);

long  get_nfft(void);

char   *zaxisNamePtr(int leveltype);
char   *vlistInqVarNamePtr(int vlistID, int varID);
char   *vlistInqVarLongnamePtr(int vlistID, int varID);
char   *vlistInqVarUnitsPtr(int vlistID, int varID);

typedef struct {
  int lana, nrecs;
  struct Variable *vars;
  struct Control *globs;
}
RARG;

#if defined(AFTERBURNER)
int stdin_is_tty  = 0;  /* true if stdin  is character device */
int stdout_is_tty = 0;  /* true if stdout is character device */
#endif

static int lstdout = 1;
 
static int Source = 0;

static int ofiletype = -1;

static int DataType = -1;

static char *filename;
static char **ifiles;
static char *ifile  = NULL;
static char *ofile  = NULL;
#if defined(AFTERBURNER)
static char *ofile2 = NULL;
#endif

static int ofileidx = 0;

static int specGridID  = -1;
static int gaussGridID = -1;
static int iVertID = -1;
static int oVertID = -1;

static int Lhybrid2pressure = FALSE;

int    TsID;
#if  defined  (HAVE_LIBPTHREAD)
int    ParallelRead = TRUE;
#else
int    ParallelRead = FALSE;
#endif

#define TIMESTEP_INTERVAL  -1
#define MONTHLY_INTERVAL    0
#define DAILY_INTERVAL      1
#define UNLIM_INTERVAL      2

#define MaxHours  24
int nrqh;
int hours[MaxHours+1];

static double *LevelFound;

static
void cdiError(int cdiErrno, const char *fmt, ...)
{
  va_list args;
	
  va_start(args, fmt);

  printf("\n");
  vfprintf(stderr, fmt, args);
   fprintf(stderr, "\n");

  va_end(args);

  fprintf(stderr, "%s\n", cdiStringError(cdiErrno));

  if ( _ExitOnError ) exit(1);
}

static
void lprintf(FILE *fp)
{
  int inum;
  int num = 67;
  int cval = '-';

  fprintf(fp, " ");
  for (inum = 0; inum < num; inum++)
    fprintf(fp, "%c", cval);
  fprintf(fp, "\n");
}

static
void FreeMean(struct Variable *vars)
{
  for ( int code = 0; code < MaxCodes; code++ )
    if ( vars[code].mean )
      {
	free(vars[code].mean);
	vars[code].mean = NULL;
      }
}


static
void after_PostProcess(struct Control *globs)
{
  if ( globs->EndOfInterval )
    {
      if ( lstdout )
	{
	  if      ( globs->OutputInterval == DAILY_INTERVAL )
	    fprintf(stdout, " Processed Day %2d  Month %2d  Year %04d",
		    globs->OldDate.dy, globs->OldDate.mo, globs->OldDate.yr);
	  else if ( globs->OutputInterval == MONTHLY_INTERVAL )
	    fprintf(stdout, " Processed Month %2d  Year %04d", globs->OldDate.mo, globs->OldDate.yr);
	  else if ( globs->OutputInterval == UNLIM_INTERVAL )
	    fprintf(stdout, " Processed range from %6.4d-%2.2d-%2.2d to %6.4d-%2.2d-%2.2d",
		    globs->StartDate.yr, globs->StartDate.mo, globs->StartDate.dy,
		    globs->OldDate.yr, globs->OldDate.mo, globs->OldDate.dy);

	  if ( globs->Mean ) fprintf(stdout, "  (Mean of %3d Terms)\n", globs->MeanCount);
	  else               fprintf(stdout, "   Terms %3d\n", globs->MeanCount);
	}

      globs->EndOfInterval = FALSE;
      globs->MeanCount = 0;
    }
}

/* ================= */
/* switch input file */
/* ================= */
static
void after_SwitchFile(struct Control *globs)
{
  int echam4 = FALSE;
  int i, n;
  char y3, y2, y1, y0;
  char         m1, m0;
  char         d1, d0;

  streamClose(globs->istreamID);

  if ( globs->Multi > 0 )
    {
      i = strlen(ifile);
      if ( i < 10 )
	{
	  fprintf(stderr, " Not a valid filename: %s \n", ifile);
	  exit(1);
	}

      if ( ifile[i-3] == '.' )
	{
	  echam4 = TRUE;
	  y3 = ifile[i-9]; y2 = ifile[i-8];
	  y1 = ifile[i-7]; y0 = ifile[i-6];
	  m1 = ifile[i-5]; m0 = ifile[i-4];
	  d1 = ifile[i-2]; d0 = ifile[i-1];
	}
      else
	{
	  y3 = ifile[i-6]; y2 = ifile[i-5];
	  y1 = ifile[i-4]; y0 = ifile[i-3];
	  m1 = ifile[i-2]; m0 = ifile[i-1];
	  d1 = '0';        d0 = '1'   ;
	}

      for ( n = 0; n < globs->DayIn; n++ )
	{
	  if ( d0 =='9' ) { d0 = '0'; d1++; }
	  else d0++;
	  if ( d1 == '3' && d0 > '0' )
	    {
	      d1 = '0'; d0 = '1';
	      if ( m1 == '0' )
		{
		  if ( m0 == '9' ) { m0 = '0'; m1 = '1'; }
		  else m0++;
		}
	      else
		{
		  if ( m0 < '2' ) m0++;
		  else
		    {
		      m1 = '0';  m0 = '1';  y0++;
		      if ( y0 > '9' ) { y0 = '0'; y1++; }
		      if ( y1 > '9' ) {
			y1 = (char) '0';
			if ( isdigit((int)y2) ) y2++;
			else                    y2 = '1';
			if ( y2 > '9' )
			  {
			    y2 = (char) '0';
			    if ( isdigit((int)y3) ) y3++;
			    else                    y3 = '1';
			  }
		      }
		    }
		}
	    }
	}

      if ( echam4 )
	{
	  ifile[i-9] = y3; ifile[i-8] = y2;
	  ifile[i-7] = y1; ifile[i-6] = y0;
	  ifile[i-5] = m1; ifile[i-4] = m0;
	  ifile[i-2] = d1; ifile[i-1] = d0;
	}
      else
	{
	  ifile[i-6] = y3; ifile[i-5] = y2;
	  ifile[i-4] = y1; ifile[i-3] = y0;
	  ifile[i-2] = m1; ifile[i-1] = m0;
	}

      globs->Multi--;
    }

  if ( globs->Nfiles > 0 ) ifile = ifiles[--globs->Nfiles];

  fprintf(stderr, " Continuation file: %s\n", ifile);

  globs->istreamID = streamOpenRead(ifile);
  if ( globs->istreamID < 0 ) cdiError(globs->istreamID, "Open failed on %s", ifile);

  globs->ivlistID = streamInqVlist(globs->istreamID);
  globs->taxisID  = vlistInqTaxis(globs->ivlistID);
}

static
int after_getDate(struct Date datetime)
{
  return cdiEncodeDate(datetime.yr, datetime.mo, datetime.dy);
}

static
int after_getTime(struct Date datetime)
{
  return cdiEncodeTime(datetime.hr, datetime.mn, 0);
}

static
void after_setDateTime(struct Date *datetime, int date, int time)
{
  int sec;
  cdiDecodeDate(date, &datetime->yr, &datetime->mo, &datetime->dy);
  cdiDecodeTime(time, &datetime->hr, &datetime->mn, &sec);
}

static
void after_printProcessStatus(int tsID)
{
  static int counthead = FALSE;

  if ( tsID == -1 )
    {
      if ( stdout_is_tty )
	{
	  fprintf(stdout, "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
	  fflush(stdout);
	}

      counthead = FALSE;
    }
  else
    {
      if ( counthead == FALSE )
	{
	  if ( stdout_is_tty )
	    fprintf(stdout, " Process timestep :       ");

	  counthead = TRUE;
	}

      if ( stdout_is_tty )
	{
	  fprintf(stdout, "\b\b\b\b\b\b%6d", tsID);
	  fflush(stdout);
	}
    }
}

static
int after_setNextDate(struct Control *globs)
{
  int nrecs;
  int i;
  int vdate, vtime;

  int righttime = FALSE;
  while ( TRUE )
    {
      nrecs = streamInqTimestep(globs->istreamID, TsID);
      if ( nrecs == 0 && ( globs->Multi > 0 || globs->Nfiles > 0 ) )
	{
	  if ( lstdout ) after_printProcessStatus(-1);

	  after_SwitchFile(globs);
	 
	  if ( globs->istreamID >= 0 )
	    {
	      TsID = 0;
	      nrecs = streamInqTimestep(globs->istreamID, TsID);
	    }
	}
      if ( nrecs == 0 ) break;

#if defined(CDO)
      //      processDefTimesteps(globs->istreamID);
#endif
      vdate = taxisInqVdate(globs->taxisID);
      vtime = taxisInqVtime(globs->taxisID);

      after_setDateTime(&globs->NextDate, vdate, vtime);

      for ( i = 0; i < nrqh; i++ )
	if ( hours[i] < 0 || hours[i] == globs->NextDate.hr )
	  {
	    righttime = TRUE;
	    break;
	  }

      if ( righttime )
	break;
      else
	TsID += 1;	  
    }

  return (nrecs);
}


static int num_recs = 0;

static
void *after_readTimestep(void *arg)
{
  int i;
  int recID, varID, gridID, zaxisID, levelID, timeID;
  int code, leveltype;
  int nmiss;
  int analysisData, nrecs;
  struct Variable *vars;
  struct Control *globs;
  RARG *rarg = (RARG *) arg;

  nrecs        = rarg->nrecs;
  analysisData = rarg->lana;
  vars         = rarg->vars;
  globs        = rarg->globs;

  for ( code = 0; code < MaxCodes; code++ ) vars[code].nmiss0 = 0;

  int level = 0;
  int levelOffset = 0;

  for ( recID = 0; recID < nrecs; recID++ )
    {
      streamInqRecord(globs->istreamID, &varID, &levelID);

      code = vlistInqVarCode(globs->ivlistID, varID);
      if ( code <= 0 || code >= MaxCodes ) continue;

      /* Skip records containing unneeded codes */

      if ( ! vars[code].needed0 ) continue;

      vlistInqVar(globs->ivlistID, varID, &gridID, &zaxisID, &timeID);

      leveltype = zaxisInqType(zaxisID);
	  
      /* Skip records with unselected levels */

      levelOffset = -1;
      /*
	if ( vars[code].ozaxisID != vars[code].izaxisID && ! Lhybrid2pressure )
      */
      if ( (vars[code].ozaxisID != vars[code].izaxisID) && (leveltype == ZAXIS_PRESSURE) )
	{
	  level = (int) zaxisInqLevel(zaxisID, levelID);
	  for ( i = 0; i < globs->NumLevelRequest; ++i )
	    {
	      if ( IS_EQUAL(globs->LevelRequest[i], level) )
		{
		  levelOffset = i;
		  break;
		}
	    }

	  if ( levelOffset < 0 ) continue;

	  zaxisID = vars[code].ozaxisID;
	  levelID = levelOffset;
	}

      if ( globs->Debug )
	{
	  fprintf(stderr, "T%d", globs->Truncation);

	  fprintf(stderr, "  Code %3d   Level%6d   %6.4d-%2.2d-%2.2d  %2.2d:%2.2d:00\n",
		  code, (int) zaxisInqLevel(zaxisID, levelID),
		  globs->OldDate.yr, globs->OldDate.mo, globs->OldDate.dy, globs->OldDate.hr, globs->OldDate.mn);
	}

      streamReadRecord(globs->istreamID, globs->Field, &nmiss);

      if ( analysisData )
	after_AnalysisAddRecord(globs, vars, code, gridID, zaxisID, levelID, nmiss);
      else
	after_EchamAddRecord(globs, vars, code, gridID, zaxisID, levelID, nmiss);

      if ( iVertID != -1 && oVertID != -1 && (vars[code].izaxisID == iVertID) )
	vars[code].ozaxisID = oVertID;
    }

  TsID++;
  /*
    printf("%3d  date = %d  time = %04d\n", TsID, vdate, vtime);
  */
  num_recs = after_setNextDate(globs);

  return ((void *) &num_recs);
}

static
void after_defineNextTimestep(struct Control *globs)
{
  static int otsID = 0;
  int vdate = after_getDate(globs->OldDate);
  int vtime = after_getTime(globs->OldDate);
  taxisDefVdate(globs->taxisID2, vdate);
  taxisDefVtime(globs->taxisID2, vtime);

  if ( globs->Mean != 2 )
    {
      if ( otsID == 0 )
	{
	  vlistDefTaxis(globs->ovlistID, globs->taxisID2);
	  streamDefVlist(globs->ostreamID, globs->ovlistID);
	}
      taxisDefNumavg(globs->taxisID2, globs->MeanCount+1);
      streamDefTimestep(globs->ostreamID,  otsID);
    }

#if defined(AFTERBURNER)
  if ( globs->Mean >= 2 )
    {
      if ( otsID == 0 )
	{
	  vlistDefTaxis(globs->ovlistID2, globs->taxisID2);
	  streamDefVlist(globs->ostreamID2, globs->ovlistID2);
	}
      taxisDefNumavg(globs->taxisID2, globs->MeanCount+1);
      streamDefTimestep(globs->ostreamID2, otsID);
    }
#endif

  otsID++;
}

static
void after_setEndOfInterval(struct Control *globs, int nrecs)
{
  if ( nrecs == 0 )
    {
      globs->EndOfInterval = TRUE;
    }
  else
    {
      if      ( globs->OutputInterval == DAILY_INTERVAL )
	globs->EndOfInterval = globs->NewDate.dy != globs->OldDate.dy;
      else if ( globs->OutputInterval == MONTHLY_INTERVAL )
	globs->EndOfInterval = globs->NewDate.mo != globs->OldDate.mo;
      else if ( globs->OutputInterval == UNLIM_INTERVAL )
	globs->EndOfInterval = FALSE;
      else
	Error( "output interval %d not implemented!\n", globs->OutputInterval);
    }
}

static
void after_moveTimestep(struct Variable *vars)
{
  int code;

  for ( code = 0; code < MaxCodes; code++ )
    vars[code].nmiss = vars[code].nmiss0;

  for ( code = 0; code < MaxCodes; code++ )
    if ( vars[code].hybrid0 )
      {
	vars[code].hybrid  = vars[code].hybrid0;
	vars[code].hybrid0 = NULL;
      }

  for ( code = 0; code < MaxCodes; code++ )
    if ( vars[code].spectral0 )
      {
	vars[code].spectral  = vars[code].spectral0;
	vars[code].spectral0 = NULL;
      }

  for ( code = 0; code < MaxCodes; code++ )
    if ( vars[code].grid0 )
      {
	vars[code].grid  = vars[code].grid0;
	vars[code].grid0 = NULL;
      }
}

static
void after_control(struct Control *globs, struct Variable *vars)
{
  int i;
  int tsFirst, nrecs;
  int rdate, rtime;
  int vdate, vtime;
  int code;
  RARG rarg;
  void *statusp = NULL;
#if  defined  (HAVE_LIBPTHREAD)
  pthread_t thrID;
  pthread_attr_t attr;
  int rval;

  if ( ParallelRead )
    {
      size_t stacksize;

      pthread_attr_init(&attr);
      pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
      int status = pthread_attr_getstacksize(&attr, &stacksize);
      if ( stacksize < 2097152 )
	{
	  stacksize = 2097152;
	  pthread_attr_setstacksize(&attr, stacksize);
	}
    }
#endif

  for ( code = 0; code < MaxCodes; code++ )
    vars[code].needed0 = vars[code].needed;

  TsID = 0;

  int righttime = FALSE;
  while ( (nrecs = streamInqTimestep(globs->istreamID, TsID)) > 0 )
    {
      vdate = taxisInqVdate(globs->taxisID);
      vtime = taxisInqVtime(globs->taxisID);
      after_setDateTime(&globs->StartDate, vdate, vtime);
      after_setDateTime(&globs->NewDate, vdate, vtime);

      for ( i = 0; i < nrqh; i++ )
	if ( hours[i] < 0 || hours[i] == globs->NewDate.hr )
	  {
	    righttime = TRUE;
	    break;
	  }

      if ( righttime )
	break;
      else
	TsID++;	  
    }

  if ( taxisInqType(globs->taxisID) == TAXIS_RELATIVE )
    {
      rdate = taxisInqRdate(globs->taxisID);
      rtime = taxisInqRtime(globs->taxisID);
    }
  else
    {
      rdate = after_getDate(globs->StartDate);
      rtime = after_getTime(globs->StartDate);
    }

  if ( ofiletype == FILETYPE_NC || ofiletype == FILETYPE_NC2 || ofiletype == FILETYPE_NC4 )
    {
      taxisDefCalendar(globs->taxisID2, CALENDAR_PROLEPTIC);
      taxisDefType(globs->taxisID2, TAXIS_RELATIVE);
      taxisDefTunit(globs->taxisID2, TUNIT_DAY);
      taxisDefRdate(globs->taxisID2, rdate);
      taxisDefRtime(globs->taxisID2, rtime);
    }

  globs->OldDate = globs->NewDate;

  tsFirst = TRUE;

  while ( nrecs > 0 )
    {
      rarg.nrecs = nrecs;
      rarg.lana  = globs->AnalysisData;
      rarg.vars  = vars;
      rarg.globs = globs;

      if ( tsFirst || ParallelRead == FALSE )
	{
	  if ( ParallelRead == FALSE )
	    {
	      statusp = after_readTimestep(&rarg);
	    }
#if  defined  (HAVE_LIBPTHREAD)
	  else
	    {
	      rval = pthread_create(&thrID, &attr, after_readTimestep, &rarg);
	      if ( rval != 0 ) Error( "pthread_create failed!");
	    }
#endif

	  if ( tsFirst )
	    {
	      if ( globs->Type >  0 ) after_legini_setup(globs, vars);
	    }

#if  defined  (HAVE_LIBPTHREAD)
	  if ( ParallelRead )
	    {
	      pthread_join(thrID, &statusp);
	      if ( *(int *)statusp < 0 )
		Error( "after_readTimestep error! (status = %d)", *(int *)statusp);
	    }
#endif
	  tsFirst = FALSE;
	}
#if  defined  (HAVE_LIBPTHREAD)
      else
	{
	  pthread_join(thrID, &statusp);
	  if ( *(int *)statusp < 0 )
	    Error( "after_readTimestep error! (status = %d)", *(int *)statusp);
	}
#endif
      nrecs = *(int *)statusp;

      globs->MeanCount0 = globs->MeanCount;
      globs->NewDate = globs->NextDate;

      after_moveTimestep(vars);

#if  defined  (HAVE_LIBPTHREAD)
      if ( nrecs && ParallelRead )
	{
	  rval = pthread_create(&thrID, &attr, after_readTimestep, &rarg);
	  if ( rval != 0 ) Error( "pthread_create failed!");
	}
#endif

      after_setEndOfInterval(globs, nrecs);
      
      if ( lstdout ) after_printProcessStatus(TsID);

      if ( lstdout && globs->EndOfInterval ) after_printProcessStatus(-1);
      
      if ( globs->Mean == 0 || globs->EndOfInterval ) after_defineNextTimestep(globs);

      if ( globs->AnalysisData )
	after_processPL(globs, vars);
      else
	after_processML(globs, vars);

      after_PostProcess(globs);
      
      if ( nrecs )
	{
	  if ( globs->AnalysisData )
	    after_AnalysisDependencies(vars, MaxCodes);
	  else
	    after_EchamDependencies(vars, MaxCodes, globs->Type, Source);
	}
	
      globs->OldDate = globs->NewDate;
    }

#if  defined  (HAVE_LIBPTHREAD)
  if ( ParallelRead )
    {
      pthread_attr_destroy(&attr);
    }
#endif
}

static
void after_setLevel(struct Control *globs)
{
  int k, l, found;
  int removeLevel[MaxLevel];
  double level;
  int checkLevel = TRUE;
  int numplevelDefault;  /* default pressure level */
  long plevelDefault[] = { 100000, 92500, 85000, 70000, 60000, 50000, 40000, 30000, 25000, 20000, 15000,
                            10000,  7000,  5000,  3000,  2000, 1000 };
  int numhlevelDefault;  /* default height level */
  long hlevelDefault[] = {  0, 1000, 2000, 5000, 10000, 15000, 20000, 25000, 30000 };

  numplevelDefault = sizeof(plevelDefault)/sizeof(plevelDefault[0]);
  numhlevelDefault = sizeof(hlevelDefault)/sizeof(hlevelDefault[0]);

  if ( iVertID != -1 )
    if ( zaxisInqType(iVertID) == ZAXIS_HYBRID && globs->Type > 20 )
      Lhybrid2pressure = TRUE;

  if ( globs->Verbose ) lprintf(stdout);

  if ( globs->NumLevelRequest == 0 )
    {
      if ( iVertID == -1 )
	{
	  if ( globs->Verbose ) fprintf(stdout," No level detected\n");
	}
      else
	{
	  if ( Lhybrid2pressure )
	    {
	      if ( globs->unitsel == 0 )
		{
		  if ( globs->Verbose ) fprintf(stdout," Default pressure level selected:\n");
		  globs->NumLevelRequest = numplevelDefault;
		  for ( l = 0; l < globs->NumLevelRequest; l++ ) globs->LevelRequest[l] = plevelDefault[l];
		  oVertID = zaxisCreate(ZAXIS_PRESSURE, globs->NumLevelRequest);
		  zaxisDefLevels(oVertID, globs->LevelRequest);
		}
	      else
		{
		  if ( globs->Verbose ) fprintf(stdout," Default height level selected:\n");
		  globs->NumLevelRequest = numhlevelDefault;
		  for ( l = 0; l < globs->NumLevelRequest; l++ ) globs->LevelRequest[l] = hlevelDefault[l];
		  oVertID = zaxisCreate(ZAXIS_HEIGHT, globs->NumLevelRequest);
		  zaxisDefLevels(oVertID, globs->LevelRequest);
		}
	    }
	  else
	    {
	      if ( globs->Verbose )
		{
		  if ( zaxisInqType(iVertID) == ZAXIS_HYBRID )
		    fprintf(stdout," All detected hybrid level selected:\n");
		  else
		    fprintf(stdout," All detected pressure level selected:\n");
		}
	      globs->NumLevelRequest = globs->NumLevelFound;
	      for ( l = 0; l < globs->NumLevelRequest; l++ ) globs->LevelRequest[l] = LevelFound[l];
	      oVertID = iVertID;
	    }
	}
      checkLevel = FALSE;
    }
  else
    {
      if ( iVertID == -1 )
	{
	  if ( globs->Verbose )	fprintf(stdout," No level detected\n");
	  checkLevel = FALSE;
	}
      else if ( globs->NumLevelRequest == 1 && IS_EQUAL(globs->LevelRequest[0], 0) )
	{
	  if ( globs->Verbose ) fprintf(stdout," No level selected\n");
	  globs->NumLevelRequest = 0;
	  checkLevel = FALSE;
	}
      else if ( globs->Verbose )
	{
	  if ( Lhybrid2pressure )
	    {
	      if ( globs->unitsel == 0 )
		fprintf(stdout," Selected pressure level:\n");
	      else
		fprintf(stdout," Selected height level:\n");
	    }
	  else
	    {
	      if ( zaxisInqType(iVertID) == ZAXIS_HYBRID )
		fprintf(stdout," Selected hybrid level:\n");
	      else
		{
		  if ( globs->unitsel == 0 )
		    fprintf(stdout," Selected pressure level:\n");
		  else
		    fprintf(stdout," Selected height level:\n");
		}
	    }
	}
    }

  if ( globs->Verbose && iVertID != -1 )
    for ( l = 0; l < globs->NumLevelRequest; l++ )
      fprintf(stdout, "  Level %2d = %13.4f\n", l+1, globs->LevelRequest[l]);

  if ( checkLevel )
    {
      for ( k = 0; k < globs->NumLevelRequest; k++ ) removeLevel[k] = FALSE;
      for ( k = 0; k < globs->NumLevelRequest; k++ )
	{
	  level = globs->LevelRequest[k];
	  for ( l = k+1; l < globs->NumLevelRequest; l++ )
	    if ( removeLevel[l] == FALSE && IS_EQUAL(level, globs->LevelRequest[l]) )
	      {
		if ( globs->Verbose )
		  fprintf(stdout, "  Level %2d = %13.4f double request\n", l+1, globs->LevelRequest[l]);
		removeLevel[l] = TRUE;
	      }
	}

      l = 0;
      for ( k = 0; k < globs->NumLevelRequest; k++ )
	if ( removeLevel[k] == FALSE )
	  globs->LevelRequest[l++] = globs->LevelRequest[k];

      globs->NumLevelRequest = l;

      if ( globs->AnalysisData || globs->Type < 30 )
	{
	  for ( k = 0; k < globs->NumLevelRequest; k++ ) removeLevel[k] = FALSE;
	  for ( k = 0; k < globs->NumLevelRequest; k++ )
	    {
	      level = globs->LevelRequest[k];
	      found = FALSE;
	      for ( l = 0; l < globs->NumLevelFound; l++ )
		if ( IS_EQUAL(level, LevelFound[l]) ) found = TRUE;

	      if ( ! found )
		{
		  fprintf(stdout, "  Level %2d = %14.4f not in input\n",
			  k+1, globs->LevelRequest[k]);
		  removeLevel[k] = TRUE;
		}
	    }

	  l = 0;
	  for ( k = 0; k < globs->NumLevelRequest; k++ )
	    if ( removeLevel[k] == FALSE )
	      globs->LevelRequest[l++] = globs->LevelRequest[k];

	  if ( l != globs->NumLevelRequest )
	    {
	      extern int labort_after;
	      if ( globs->Verbose ) lprintf(stdout);
	      if ( labort_after )
		Error( "Inconsistent or invalid level list!");
	      else
		Warning( "Inconsistent or invalid level list!");
	    }

	  globs->NumLevelRequest = l;
	}
    }

  if ( globs->Verbose ) lprintf(stdout);
}

static
void after_defineLevel(struct Control *globs, struct Variable *vars)
{
  int code, i;

  /* hybrid, pressure, height */

  switch ( globs->Type ) {
  case  0:
  case 10:
  case 11:
  case 20:
    {
      if ( iVertID == -1 ) break;

      if ( zaxisInqType(iVertID) == ZAXIS_HYBRID )
	{
	  if ( oVertID == -1 )
	    {
	      if ( globs->NumLevelRequest > globs->NumLevelFound )
		Error( "Too much level requested");

	      if ( globs->NumLevelFound == globs->NumLevelRequest )
		{
		  for ( i = 0; i < globs->NumLevelRequest; i++ )
		    if ( IS_NOT_EQUAL(globs->LevelRequest[i], LevelFound[i]) ) break;

		  if ( i == globs->NumLevelRequest )
		    oVertID = iVertID;
		}

	      if ( oVertID == -1 && globs->NumLevelRequest > 0 )
		{
		  oVertID = zaxisCreate(ZAXIS_HYBRID, globs->NumLevelRequest);
		  zaxisDefLevels(oVertID, globs->LevelRequest);
		  zaxisDefVct(oVertID, globs->nvct, globs->vct);
		}
	    }
	      
	  for ( code = 0; code < MaxCodes; code++ )
	    {
	      if ( vars[code].selected )
		{
		  if ( vars[code].izaxisID != -1 )
		    if ( zaxisInqType(vars[code].izaxisID) == ZAXIS_HYBRID &&
			 zaxisInqSize(vars[code].izaxisID) >= globs->NumLevelRequest )
		      vars[code].ozaxisID = oVertID;
		}
	    }
	}
      else
	Error( "%s level data unsupported for TYPE %d",
	      zaxisNamePtr(zaxisInqType(iVertID)), globs->Type);	    

      break;
    }
  case 30:
  case 40:
  case 41:
  case 50:
  case 60:
  case 61:
  case 70:
    {
      if ( iVertID == -1 ) break;

      if ( oVertID == -1 )
	{
	  if ( globs->unitsel == 0 )
	    oVertID = zaxisCreate(ZAXIS_PRESSURE, globs->NumLevelRequest);
	  else
	    oVertID = zaxisCreate(ZAXIS_HEIGHT, globs->NumLevelRequest);

	  zaxisDefLevels(oVertID, globs->LevelRequest);
	}

      for ( code = 0; code < MaxCodes; code++ )
	{
	  if ( vars[code].selected )
	    {
	      if ( vars[code].izaxisID != -1 )
		if ( zaxisInqType(vars[code].izaxisID) == zaxisInqType(iVertID) &&
		     zaxisInqSize(vars[code].izaxisID) == globs->NumLevel &&
		     zaxisInqSize(vars[code].izaxisID) >  1 )
		  vars[code].ozaxisID = oVertID;
	    }
	}

      break;
    }
  default:
    Error( "TYPE %d unsupported", globs->Type);
  }
}

static
void after_defineGrid(struct Control *globs, struct Variable *vars)
{
  int ogridID = -1;
  int code;

  /* spectral, fourier, gauss, zonal mean */

  switch ( globs->Type ) {
  case  0:
  case 50:
    {
      if ( specGridID == -1 )
	{
	  if ( globs->DimSP == 0 )
	    Error( "dim spectral undefined");
	  if ( globs->Truncation == 0 )
	    Error( "truncation undefined");
	    
	  specGridID = gridCreate(GRID_SPECTRAL, globs->DimSP);
	  gridDefTrunc(specGridID, globs->Truncation);
	}

      ogridID = specGridID;
      break;
    }
  case 20:
  case 30:
  case 70:
    {
      if ( gaussGridID == -1 )
	{
	  if ( globs->Longitudes == 0 )
	    Error( "number of longitudes undefined");
	  if ( globs->Latitudes == 0 )
	    Error( "number of latitudes undefined");
	    
	  gaussGridID = gridCreate(GRID_GAUSSIAN, globs->Longitudes*globs->Latitudes);
	  gridDefXsize(gaussGridID, globs->Longitudes);
	  gridDefYsize(gaussGridID, globs->Latitudes);
	}

      ogridID = gaussGridID;
      break;
    }
  case 10:
  case 40:
  case 60:
    {
      if ( globs->Fouriers == 0 )
	Error( "number of fourier coefficients undefined");
      if ( globs->Latitudes == 0 )
	Error( "number of latitudes undefined");
	    
      ogridID = gridCreate(GRID_FOURIER, globs->Fouriers*globs->Latitudes);
      gridDefXsize(ogridID, globs->Latitudes);
      gridDefYsize(ogridID, globs->Fouriers);
      break;
    }
  case 11:
  case 41:
  case 61:
    {
      if ( globs->Latitudes == 0 )
	Error( "Number of latitudes undefined");
	    
      ogridID = gridCreate(GRID_GAUSSIAN, globs->Latitudes);
      gridDefXsize(ogridID, 1);
      gridDefYsize(ogridID, globs->Latitudes);
      break;
    }
  default:
    Error( "TYPE %d unsupported", globs->Type);
  }

  if ( ogridID != -1 )
    for ( code = 0; code < MaxCodes; code++ )
      {
	if ( vars[code].selected )
	  {
	    vars[code].ogridID = ogridID;
	  }
      }

  if ( ogridID == -1 )
    Error( "out grid undefined");
}

static
void after_setCodes(struct Control *globs, struct Variable *vars, int maxCodes, int numCodes)
{
  int code;
  int table, modelID, tableID;
  char *name, *longname;
  int varID;

  if ( globs->Verbose ) lprintf(stdout);

  if ( numCodes == 0 )
    {
      if ( globs->Verbose ) fprintf(stdout, " All detected codes selected:\n");
      
      for ( code = 0; code < maxCodes; code++ )
	if ( vars[code].detected ) vars[code].selected = 1;
    }
  else if ( globs->Verbose )
    fprintf(stdout, " Selected codes:\n");

  if ( globs->Verbose )
    {
      fprintf(stdout, "  Table Code Name              Longname\n");
      fprintf(stdout, "  ----- ---- ----              --------\n");
    }

  for ( code = 0; code < maxCodes; code++ )
    if ( vars[code].selected )
      {
	table    = 0;
	name     = NULL;
	longname = NULL;
	varID    = vars[code].ivarID;

	if ( varID == CDI_UNDEFID )
	  {
	    modelID  = vlistInqVarModel(globs->ivlistID, 0);
	    table    = 128;
	    tableID  = tableInq(modelID, table, NULL);

	    vars[code].tableID = tableID;
	  }
	else
	  {
	    tableID  = vlistInqVarTable(globs->ivlistID, varID);
	    table    = tableInqNum(tableID);
	    name     = vlistInqVarNamePtr(globs->ivlistID, varID);
	    longname = vlistInqVarLongnamePtr(globs->ivlistID, varID);
	  }
	if ( name     == NULL )     name = (char*) tableInqParNamePtr(tableID, code);
	if ( longname == NULL ) longname = (char*) tableInqParLongnamePtr(tableID, code);
	
	if ( globs->Verbose )
	  {
	    fprintf(stdout, " %5d", table);
	    fprintf(stdout, " %4d", code);
	    if ( name == NULL )
	      fprintf(stdout, "  var%d", code);
	    else
	      {
		fprintf(stdout, "  %-16s", name);
	    if ( longname != NULL )
	      fprintf(stdout, "  %s", longname);
	      }
	    fprintf(stdout, "\n");
	  }
      }
}

static
void after_checkNamelist(struct Control *globs)
{
  if ( globs->Mean && globs->Type < 20 )
    {
      Error("Mean is only available for TYPE >= 20!");
    }

  if ( globs->Extrapolate == FALSE && globs->Type >= 30 )
    {
      if ( globs->Type > 30 )
	Error("EXTRAPOLATE = 0 is only available for TYPE = 30!");
      if ( globs->Mean )
	Error("EXTRAPOLATE = 0 is only available with MEAN = 0!");
    }
}

static
void after_usage(void)
{
  fprintf(stderr, "\nafter [options] <InputFiles> <OutputFile> <VarianceFile>\n");
#if defined (_OPENMP)
  fprintf(stderr, "     option -P <nthreads> : Set number of OpenMP threads\n");
#endif
  fprintf(stderr, "     option -a            : Forces analysis data process\n");
  fprintf(stderr, "     option -c            : Print available codes and names\n");
  fprintf(stderr, "     option -d            : Debug mode\n");
  fprintf(stderr, "     option -v <vctfile>  : Read vct from vctfile\n");
  /*  fprintf(stderr, "     option -h : help (this output)\n"); */
  /*  fprintf(stderr, "     option -p : parallel read on\n"); */
  fprintf(stderr, "  <InputFiles> : ECHAM or ECMWF Ana or ReAna files\n");
  fprintf(stderr, "  <OutputFile> : GRIB, netCDF or SERVICE format file\n");
  fprintf(stderr, "<VarianceFile> : GRIB, netCDF or SERVICE format file\n");
  fprintf(stderr, "  namelist is read from <stdin>\n");
  fprintf(stderr, "  output is written to <stdout>\n\n");

  fprintf(stderr, "  default Namelist: \n");
  fprintf(stderr, "  &SELECT\n");
  fprintf(stderr, "    TYPE = 0, CODE = -1, LEVEL = -1, MULTI = 0, DAYIN = 30,\n");
  fprintf(stderr, "    MEAN = 0, TIMESEL = -1, UNITSEL = 0,\n");
  fprintf(stderr, "    FORMAT = 0, PRECISION = 0, SZIP = 0\n");
  fprintf(stderr, "  &END\n");

  exit(1);
}

static
void after_parini(struct Control *globs, struct Variable *vars)
{
  char namelist[65536];

  if ( stdin_is_tty ) 
    {
#if defined(CDO)
      fprintf(stderr, "Default namelist: \n");
      fprintf(stderr, "  TYPE = 0, CODE = -1, LEVEL = -1, MULTI = 0, DAYIN = 30, MEAN = 0, TIMESEL = -1, UNITSEL = 0\n");
#endif
      fprintf(stdout, "Enter namelist parameter:\n");
    }
  else
    {
      fseek(stdin, 0L, SEEK_END);
      long length = ftell(stdin);
      if ( length == 0L )
	{
	  fprintf(stderr,"\n stdin not connected\n");
	  after_usage();
	}
      fseek(stdin, 0L, SEEK_SET);
    }

  int i = 1;
  namelist[0] = ' ';
  int c = getchar();
  while ((c != EOF) && i < (int)(sizeof(namelist) - 1))
    {
           if ((c >= '0' && c <= '9') ||
               (c == '-' || c == '.'))  namelist[i++] = c;
      else if  (c >= 'a' && c <= 'z')   namelist[i++] = c;
      else if  (c >= 'A' && c <= 'Z')   namelist[i++] = tolower(c);
      else c = ' ';

      if (c == ' ' && namelist[i-1] != ' ') namelist[i++] = c;
      c = getchar();
    }
  namelist[i] = 0;

  if ( globs->Debug )
    {
      lprintf(stderr);
      fprintf(stderr,"  Length of namelist:%4d bytes\n", (int) strlen(namelist));

      for (i = 0; i < (int)strlen(namelist); i += 60)
	fprintf(stderr,"  namelist[%02d]=%-60.60s\n", i, namelist+i);
      lprintf(stderr);
    }

  if ( globs->Verbose )
    {
      lprintf(stdout);
      fprintf(stdout, " Namelist:\n");
    }
  
  globs->Type           = scan_par(globs->Verbose, namelist, "type",  0);
  globs->Multi          = scan_par(globs->Verbose, namelist, "multi", 0);
  globs->Mean           = scan_par(globs->Verbose, namelist, "mean",  0);
  globs->OutputInterval = scan_par(globs->Verbose, namelist, "interval", MONTHLY_INTERVAL);

#if defined(CDO)
  if ( globs->Mean >= 2 )
    cdoAbort("Namelist parameter MEAN=%d out of bounds (0:1)", globs->Mean);
#endif

  int fileFormat = scan_par(globs->Verbose, namelist, "format", -1);
  int gribFormat = scan_par_obsolate(namelist, "grib",   0);
  int cdfFormat  = scan_par_obsolate(namelist, "netcdf", 0);

  if ( gribFormat && cdfFormat ) Error( "GRIB or netCDF?");

  switch ( fileFormat )
    {
#if defined(CDO)
    case -1: ofiletype = -1;            break;
#else
    case -1: ofiletype = FILETYPE_SRV;  break;
#endif
    case  0: ofiletype = FILETYPE_SRV;  break;
    case  1: ofiletype = FILETYPE_GRB;  break;
    case  2: ofiletype = FILETYPE_NC;   break;
    case  3: ofiletype = FILETYPE_EXT;  break;
    case  4: ofiletype = FILETYPE_NC2;  break;
    case  6: ofiletype = FILETYPE_NC4;  break;
    default: Error( "unknown file format %d", fileFormat);
    }

  if ( gribFormat )  ofiletype = FILETYPE_GRB;
  if ( cdfFormat  )  ofiletype = FILETYPE_NC;

  int precision = scan_par(globs->Verbose, namelist, "precision", 0);
  if ( precision )
    switch ( precision )
      {
      case  8: DataType = DATATYPE_PACK8;  break;
      case 16: DataType = DATATYPE_PACK16; break;
      case 24: DataType = DATATYPE_PACK24; break;
      case 32: DataType = DATATYPE_FLT32;  break;
      case 64: DataType = DATATYPE_FLT64;  break;
      default: Error( "unsupported data precision %d", precision);
      }


  globs->unitsel        = scan_par(globs->Verbose, namelist, "unitsel",  0);
  globs->DayIn          = scan_par(globs->Verbose, namelist, "dayinc",  30);
  globs->Extrapolate    = scan_par(globs->Verbose, namelist, "extrapolate",  1);
  globs->Szip           = scan_par(globs->Verbose, namelist, "szip",  0);
  int mars              = scan_par_obsolate(namelist, "mars", 0);

  if ( globs->Multi ) --globs->Multi;

  if ( mars )
    {
      extern int Mars;
      Mars = 1;
      PlanetRD     = C_MARS_RD;
      PlanetGrav   = C_MARS_GRAV;
      PlanetRadius = C_MARS_RADIUS;
    }

  nrqh = scan_time(globs->Verbose, namelist, hours, MaxHours);
  scan_code(namelist, vars, MaxCodes, &globs->NumCodesRequest);

  scan_darray(namelist, "level", globs->LevelRequest, MaxLevel, &globs->NumLevelRequest);
  if ( globs->NumLevelRequest == 1 )
    if ( IS_EQUAL(globs->LevelRequest[0], -1) )
      globs->NumLevelRequest = 0;

  if ( globs->Verbose ) lprintf(stdout);

  after_checkNamelist(globs);
}

static
void after_dimcalc(struct Control *globs)
{  
  if ( globs->AnalysisData ) globs->NumLevel = globs->NumLevelRequest;

  if ( globs->Latitudes == 0 )
    {
      globs->Latitudes = 2 * ((globs->Truncation*3 + 3) / 4);
      if ( globs->Truncation == 30 ) globs->Latitudes = 48;
    }

  if ( globs->Longitudes == 0 )
    {
      globs->Longitudes = globs->Latitudes * 2;
      if ( globs->Truncation == 62 ) globs->Longitudes = 192;
    }

  globs->Waves         = globs->Truncation + 1;
  globs->Fouriers      = globs->Waves * 2;
  globs->DimSP         = (globs->Truncation + 1) * (globs->Truncation + 2);
  globs->DimFC         = globs->Latitudes * globs->Fouriers;
  globs->DimGP         = globs->Latitudes * globs->Longitudes;
  globs->Dim3GP        = globs->NumLevel * globs->DimGP;
  globs->Dim3FC        = globs->NumLevel * globs->DimFC;
  globs->Dim3SP        = globs->NumLevel * globs->DimSP;
  globs->HalfLevels    = globs->NumLevel + 1;
  globs->DimSP_half    = globs->DimSP / 2;

  if ( globs->AnalysisData )
    fprintf(stdout, " Found Ana or Re-Ana Data\n");

  if ( globs->Verbose )
    {
      fprintf(stdout, " Dimensions:\n");
      fprintf(stdout, "  Truncation        = %4d\n", globs->Truncation);
      fprintf(stdout, "  Levels            = %4d\n", globs->NumLevel);
      fprintf(stdout, "  Latitudes         = %4d\n", globs->Latitudes);
      fprintf(stdout, "  Longitudes        = %4d\n", globs->Longitudes);
      lprintf(stdout);
    }
}

/* ----------------------------------------------------------- */
/* Extract basic dimension information                         */
/* ----------------------------------------------------------- */
static
void after_precntl(struct Control *globs, struct Variable *vars)
{
  int l;
  int code = 0;
  int gridID, zaxisID, varID, timeID;
  int i, index, leveltype, gridtype;
  int datasize, numlevel;
  int vertfound = 0;
  int nhzaxis = 0;
  int FieldDim = 0;

  int nvars   = vlistNvars(globs->ivlistID);
  int ngrids  = vlistNgrids(globs->ivlistID);
  int nverts  = vlistNzaxis(globs->ivlistID);
  int ntsteps = vlistNtsteps(globs->ivlistID);

  if ( globs->Debug )
    {
      Message( "nvars      = %d", nvars);
      Message( "ngrids     = %d", ngrids);
      Message( "nverts     = %d", nverts);
      Message( "ntsteps    = %d", ntsteps);
    }

  for ( index = 0; index < ngrids; index++ )
    {
      gridID   = vlistGrid(globs->ivlistID, index);
      gridtype = gridInqType(gridID);
      datasize = gridInqSize(gridID);

      if ( datasize > FieldDim ) FieldDim = datasize;

      if ( gridtype == GRID_SPECTRAL && globs->Truncation == 0 )
	{
	  specGridID = gridID;
	  globs->Truncation = gridInqTrunc(gridID);
	}
      else if ( gridtype == GRID_GAUSSIAN && globs->Latitudes == 0 )
	{
	  gaussGridID = gridID;
	  globs->Longitudes  = gridInqXsize(gridID);
	  globs->Latitudes   = gridInqYsize(gridID);
	}
    }

  if ( globs->Truncation == 0 && globs->Latitudes == 0 )
    Error("Unsupported file structure!\n");

  if ( globs->Truncation == 0 )
    {
      if ( globs->Latitudes )
	{
	  switch ( globs->Latitudes ) {
	  case 512: globs->Truncation = 511; break;
	  case 320: globs->Truncation = 213; break;
	  case 192: globs->Truncation = 127; break;
	  case 160: globs->Truncation = 106; break;
	  case 128: globs->Truncation =  85; break;
	  case  96: globs->Truncation =  63; break;
	  case  94: globs->Truncation =  62; break;
	  case  64: globs->Truncation =  42; break;
	  case  48: globs->Truncation =  31; break;
	  case  32: globs->Truncation =  21; break;
	  default :
	    fprintf(stderr,"%d Gaussian latitudes not supported.\n", globs->Latitudes);
	  }
	}
    }

  for ( index = 0; index < nverts; index++ )
    {
      zaxisID   = vlistZaxis(globs->ivlistID, index);
      leveltype = zaxisInqType(zaxisID);
      numlevel  = zaxisInqSize(zaxisID);
      /*
	printf("leveltype : %d %d\n", leveltype, zaxisInqSize(zaxisID));
      */	
      if ( numlevel > 1 )
	{
	  if ( leveltype == ZAXIS_HYBRID || leveltype == ZAXIS_PRESSURE )
	    {
	      if ( leveltype == ZAXIS_HYBRID && globs->nvct == 0 )
		{
		  nhzaxis++;
		  if ( numlevel != (zaxisInqVctSize(zaxisID)/2 - 1) )
		    {
		      if ( ! (numlevel == 191 && zaxisInqVctSize(zaxisID) == 0) )
			{
			  Warning( "Skip %d hybrid level data with %d levels!",
				  (zaxisInqVctSize(zaxisID)/2 - 1), numlevel);
			  continue;
			}
		    }
		}

	      if ( iVertID != - 1 )
		Warning( "More than %d different vertical grid structure found!", vertfound);

	      vertfound++;

	      if ( iVertID != -1 ) continue;

	      iVertID = zaxisID;
	      globs->NumLevelFound = numlevel;
	      LevelFound = (double *) malloc(globs->NumLevelFound*sizeof(double));
	      for ( l = 0; l < globs->NumLevelFound; l++ )
		LevelFound[l] = (int) zaxisInqLevel(zaxisID, l);

	      if ( leveltype == ZAXIS_HYBRID )
		{
		  if ( globs->nvct == 0 )
		    {
		      if ( zaxisInqVctSize(zaxisID) )
			{
			  globs->nvct = zaxisInqVctSize(zaxisID);

			  if ( globs->vct == NULL )
			    {
			      globs->vct = (double *) malloc(globs->nvct*sizeof(double));
			      memcpy(globs->vct, zaxisInqVctPtr(zaxisID), globs->nvct*sizeof(double));
			    }
			}
		      else
			{
			  if ( numlevel == 191 )
			    {
			      fprintf(stderr," Using internal VCT for L191\n");
			      globs->nvct = (191+1)*2;
			      globs->vct = (double *) malloc(globs->nvct*sizeof(double));
			      memcpy(globs->vct, VCT_L191, globs->nvct*sizeof(double));
			      zaxisDefVct(zaxisID, globs->nvct, globs->vct);
			    }
			  else
			    {
			      Error( "VCT not defined in inputfile!");
			    }
			}
		    }

		  if ( numlevel != (globs->nvct/2 - 1) )
		    Error( "Number of hybrid levels %d does not match vct levels %d",
			  numlevel, globs->nvct/2-1);

		  if ( globs->Debug )
		    for ( i = 0; i < globs->nvct/2; i++ )
		      fprintf(stderr," vct: %4d %10.4f %10.4f\n", i, globs->vct[i], globs->vct[i+globs->nvct/2]);
		}

	      if ( leveltype == ZAXIS_PRESSURE ) globs->AnalysisData = TRUE;
	    }
	}
    }

  if ( nhzaxis > 0 && globs->nvct == 0 ) Error( "VCT missing!");

  globs->NumLevel = globs->NumLevelFound;

  if (  specGridID != -1 ) globs->Spectral = TRUE;
  if ( gaussGridID != -1 ) globs->Gaussian = TRUE;

  if ( globs->Debug )
    fprintf(stderr, "   T = %3d   L = %2d\n", globs->Truncation, globs->NumLevelFound);

  if ( globs->Debug )
    fprintf(stderr," CODE CHECK\n");

  if ( globs->Verbose )
    {
      int instID  = vlistInqVarInstitut(globs->ivlistID, 0);
      int modelID = vlistInqVarModel(globs->ivlistID, 0);

      lprintf(stdout);
      fprintf(stdout, " Institute : ");
      if ( instID == CDI_UNDEFID )
	fprintf(stdout, "unknown\n");
      else
	{
	  if ( institutInqLongnamePtr(instID) )
	    fprintf(stdout, "%s\n", institutInqLongnamePtr(instID));
	  else
	    fprintf(stdout, "name unknown\n");
	}

      fprintf(stdout, " Source    : ");
      if ( modelID == CDI_UNDEFID )
	fprintf(stdout, "unknown\n");
      else
	{
	  if ( modelInqNamePtr(modelID) )
	    {
	      if ( strncmp(modelInqNamePtr(modelID), "ECHAM5", 6) == 0 ) Source = S_ECHAM5;
	      fprintf(stdout, "%s\n", modelInqNamePtr(modelID));
	    }
	  else
	    fprintf(stdout, "name unknown\n");
	}
    }
  
  for ( varID = 0; varID < nvars; varID++ )
    {
      vlistInqVar(globs->ivlistID, varID, &gridID, &zaxisID, &timeID);
      code      = vlistInqVarCode(globs->ivlistID, varID);
      if ( code <= 0 || code >= MaxCodes )
	{
	  Warning( "Code number %d out of range, variable ignored!", code);
	  continue;
	}
      gridtype  = gridInqType(gridID);
      numlevel  = zaxisInqSize(zaxisID);
      leveltype = zaxisInqType(zaxisID);

      vars[code].ivarID  = varID;
      vars[code].igridID = gridID;
      vars[code].ogridID = gridID;
      vars[code].izaxisID = zaxisID;
      vars[code].ozaxisID = zaxisID;

      vars[code].detected = TRUE;

      if ( globs->Debug )
	fprintf(stderr,"Code %3d  Levels = %3d  LevelType = %3d  GridType = %3d\n",
		code, numlevel, leveltype, gridtype);
    }

  if ( globs->Debug )
    Message( "FieldDim = %d\n", FieldDim);

  globs->Field = (double *) malloc(FieldDim*sizeof(double));

  if ( globs->Debug )
    for ( code = 0; code < MaxCodes; code++ )
      {
	if ( vars[code].detected )
	  fprintf(stderr," Detected Code %3d with %3d level\n",
		  code, zaxisInqSize(vars[code].izaxisID));
      }
}

/*
 * -----------------------------------------------------------
 * Define output variables
 * -----------------------------------------------------------
 */
static
void after_postcntl(struct Control *globs, struct Variable *vars)
{
  int code = 0;
  int gridID, zaxisID;
  int ovarID, ogridID, ozaxisID;
  int ovarID2;
  int ivarID, instID, modelID, tableID;
  char *name, *longname, *units;
  char histstring[99];
  int datatype;

  sprintf(histstring, "afterburner version %s  type = %d", VERSION, globs->Type);

#if defined(AFTERBURNER)
  afterInqHistory(globs->istreamID);
  if ( globs->Mean != 2 ) afterDefHistory(globs->ostreamID, histstring);
  if ( globs->Mean >= 2 ) afterDefHistory(globs->ostreamID2, histstring);
#endif

  if ( globs->Debug ) lprintf(stdout);
  if ( globs->Debug )
    for ( code = 0; code < MaxCodes; code++ )
      if ( vars[code].detected )
	{
	  gridID = vars[code].igridID;
	  zaxisID = vars[code].izaxisID;
	  fprintf(stderr," Detected Code %3d  grid %-8s size %5d  level %2d %-8s\n",
		  code, gridNamePtr(gridInqType(gridID)), gridInqSize(gridID),
		  zaxisInqSize(zaxisID), zaxisNamePtr(zaxisInqType(zaxisID)));
	}


  if ( globs->Debug ) lprintf(stdout);
  if ( globs->Debug )
    for ( code = 0; code < MaxCodes; code++ )
      if ( vars[code].needed )
	{
	  fprintf(stderr,"   Needed Code %3d\n", code);
	}

  for ( code = 0; code < MaxCodes; code++ )
    if ( vars[code].selected )
      {
	name     = NULL;
	longname = NULL;
	units    = NULL;
	ivarID   = vars[code].ivarID;
	ogridID  = vars[code].ogridID;
	ozaxisID = vars[code].ozaxisID;

	if ( ogridID == -1 )
	  {
	    /*
	    Warning( "undefined grid for code %d", code);
	    */
	    continue;
	  }
	if ( ozaxisID == -1 )
	  {
	    /*
	    Warning( "undefined level for code %d", code);
	    */
	    continue;
	  }

	instID   = vlistInqVarInstitut(globs->ivlistID, ivarID);
	modelID  = vlistInqVarModel(globs->ivlistID, ivarID);
	tableID  = vlistInqVarTable(globs->ivlistID, ivarID);

	vars[code].missval  = vlistInqVarMissval(globs->ivlistID, ivarID);
	vars[code].samp     = NULL;

	if ( DataType != -1 )
	  datatype = DataType;
	else
	  datatype = vlistInqVarDatatype(globs->ivlistID, ivarID);

	if ( vars[code].comp )
	  {
	    tableID = vars[code].tableID;
	  }
	else
	  {
	    name     = vlistInqVarNamePtr(globs->ivlistID, ivarID);
	    longname = vlistInqVarLongnamePtr(globs->ivlistID, ivarID);
	    units    = vlistInqVarUnitsPtr(globs->ivlistID, ivarID);
	  }

	if ( globs->Mean != 2 )
	  {
	    vlistDefTaxis(globs->ovlistID, globs->taxisID2);
	    ovarID = vlistDefVar(globs->ovlistID, ogridID, ozaxisID, TIME_VARIABLE);
	    vlistDefVarCode(globs->ovlistID, ovarID, code);
	    vars[code].ovarID = ovarID;
	    vlistDefVarInstitut(globs->ovlistID, ovarID, instID);
	    vlistDefVarModel(globs->ovlistID, ovarID, modelID);
	    vlistDefVarTable(globs->ovlistID, ovarID, tableID);
	    if ( globs->Mean ) vlistDefVarTimave(globs->ovlistID, ovarID, 1);
	    if ( name )        vlistDefVarName(globs->ovlistID, ovarID, name);
	    if ( longname )    vlistDefVarLongname(globs->ovlistID, ovarID, longname);
	    if ( units )       vlistDefVarUnits(globs->ovlistID, ovarID, units);
	    vlistDefVarDatatype(globs->ovlistID, ovarID, datatype);
	    vlistDefVarMissval(globs->ovlistID, ovarID, vars[code].missval);
	  }

	if ( globs->Mean >= 2 )
	  {
	    vlistDefTaxis(globs->ovlistID2, globs->taxisID2);
	    ovarID2 = vlistDefVar(globs->ovlistID2, ogridID, ozaxisID, TIME_VARIABLE);
	    vlistDefVarCode(globs->ovlistID2, ovarID2, code);
	    vars[code].ovarID2 = ovarID2;
	    vlistDefVarInstitut(globs->ovlistID2, ovarID2, instID);
	    vlistDefVarModel(globs->ovlistID2, ovarID2, modelID);
	    vlistDefVarTable(globs->ovlistID2, ovarID2, tableID);
	    if ( globs->Mean ) vlistDefVarTimave(globs->ovlistID2, ovarID2, 1);
	    if ( name )        vlistDefVarName(globs->ovlistID2, ovarID2, name);
	    if ( longname )    vlistDefVarLongname(globs->ovlistID2, ovarID2, longname);
	    if ( units )       vlistDefVarUnits(globs->ovlistID2, ovarID2, units);
	    vlistDefVarDatatype(globs->ovlistID2, ovarID2, datatype);
	    vlistDefVarMissval(globs->ovlistID2, ovarID2, vars[code].missval);
	  }
      }

  if ( globs->Debug ) lprintf(stdout);
  if ( globs->Debug )
    for ( code = 0; code < MaxCodes; code++ )
      if ( vars[code].selected )
	{
	  gridID  = vars[code].ogridID;
	  zaxisID = vars[code].ozaxisID;
	  fprintf(stderr," Selected Code %3d  grid %-8s size %5d  level %2d %-8s\n",
		  code, gridNamePtr(gridInqType(gridID)), gridInqSize(gridID),
		  zaxisInqSize(zaxisID), zaxisNamePtr(zaxisInqType(zaxisID)));
	}
}

#if defined(AFTERBURNER)
static
void after_readVct(struct Control *globs, const char *vctfile)
{
  char line[1024];
  int i, n;
  double va, vb;

  FILE *fp = fopen(vctfile, "r");
  if ( fp == NULL ) SysError( "Open failed on %s", vctfile);

  while ( fgets(line, 1023, fp) ) globs->nvct++;

  globs->nvct *= 2;
  globs->vct = (double *) malloc(globs->nvct*sizeof(double));

  rewind(fp);
  for ( i = 0; i < globs->nvct/2; i++ )
    {
      fgets(line, 1023, fp);
      sscanf(line, "%d %lg %lg", &n, &va, &vb);
      globs->vct[i]               = va;
      globs->vct[i+globs->nvct/2] = vb;
    }
  fprintf(stdout, "  Reading VCT for %d hybrid levels from file %s\n", globs->nvct/2-1, vctfile);

  fclose(fp);
}
#endif

#if defined(AFTERBURNER)
static
void after_version(void)
{
#if defined (COMPILER)
  fprintf(stderr, "Compiler: %s\n", COMPILER);
#endif
#if defined (COMP_VERSION)
  fprintf(stderr, " version: %s\n", COMP_VERSION);
#endif
#if defined (HAVE_LIBSZ) || defined (_OPENMP)
  fprintf(stderr, "    with:");
#if defined (HAVE_LIBSZ)
  fprintf(stderr, " libsz");
#endif
#if defined (_OPENMP)
  fprintf(stderr, " OpenMP");
#endif
  fprintf(stderr, "\n");
#endif
#if defined (USER_NAME) && defined(HOST_NAME) && defined(SYSTEM_TYPE)
  fprintf(stderr, "Compiled: by %s on %s (%s) %s %s\n",
	  USER_NAME, HOST_NAME, SYSTEM_TYPE, __DATE__, __TIME__);
#endif
  cdiPrintVersion();
  fprintf(stderr, "\n");
}
#endif

static
void after_control_init(struct Control *globs)
{
  memset(globs, 0, sizeof(struct Control));

  globs->AnalysisData = 0; /* 0 = ECHAM Data, 1 = ECMWF Spectral Analyses */
  globs->DayIn       = 0; /* day increment of infiles if Multi = TRUE    */
  globs->Debug       = FALSE;
  globs->Extrapolate = TRUE;
  globs->Szip        = FALSE;

  globs->istreamID   = CDI_UNDEFID;
  globs->ostreamID   = CDI_UNDEFID;
  globs->ostreamID2  = CDI_UNDEFID;
  globs->ivlistID    = CDI_UNDEFID;
  globs->ovlistID    = CDI_UNDEFID;
  globs->ovlistID2   = CDI_UNDEFID;
  globs->taxisID     = -1;
  globs->taxisID2    = -1;
}


static
void after_variable_init(struct Variable *vars)
{
  memset(vars, 0, sizeof(struct Variable));

  vars->ivarID   = -1;
  vars->ovarID   = -1;
  vars->ovarID2  = -1;
  vars->izaxisID = -1;
  vars->ozaxisID = -1;
  vars->igridID  = -1;
  vars->ogridID  = -1;
  vars->tableID  = -1;
}

#if defined(AFTERBURNER)
static
void after_printCodes(void)
{
  int tableID = tableInq(-1, 128, "echam4");
  int ncodes;
  int codes[] = {34,35,36,131,132,135,148,149,151,156,157,259,260,261,262,263,264,268,269,270,271,275};

  ncodes = sizeof(codes)/sizeof(codes[0]);

  lprintf(stdout);

  fprintf(stdout, "  Code Name              Longname\n");
  fprintf(stdout, "  ---- ----              --------\n");

  for ( int i = 0; i < ncodes; i++ )
    {
      int code     = codes[i];
      const char *name     = tableInqParNamePtr(tableID, code);
      const char *longname = tableInqParLongnamePtr(tableID, code);

      fprintf(stdout, " %4d", code);
      if ( name == NULL )
	fprintf(stdout, "  var%d", code);
      else
	{
	  fprintf(stdout, "  %-16s", name);
	  if ( longname != NULL )
	    fprintf(stdout, "  %s", longname);
	}
      fprintf(stdout, "\n");
    }

  lprintf(stdout);
}
#endif

/* =============================================== */
/* procstat   - appends info about memory usage    */
/*              and time consumption               */
/* =============================================== */
#if defined(AFTERBURNER)
static
void after_procstat(char *procpath, int truncation)
{
  FILE *sf;
  double MaxMBytes;
  time_t tp;
  long  yy, mm, dd, hh, mi;
  char mtype[12];
  char *proc;
  char *name;
  char  stat_file[128];
  double CPUTime;

  CPUTime = ((double) clock() - starttime ) / CLOCKS_PER_SEC;

  (void) time(&tp);
  yy    = gmtime(&tp)->tm_year + 1900;
  mm    = gmtime(&tp)->tm_mon + 1;
  dd    = gmtime(&tp)->tm_mday   ;
  hh    = gmtime(&tp)->tm_hour   ;
  mi    = gmtime(&tp)->tm_min    ;
  name  = getpwuid(getuid())->pw_name;

  proc = strrchr(procpath,'/');
  if (proc == 0) proc = procpath;
  else           proc++         ;

  strcpy(stat_file, "/pf/m/m214003/local/log/after.log");

  MaxMBytes = (double) memTotal() / 1048576.;

  sf = fopen(stat_file, "a");
  if ( sf )
    {
      char unknown[] = "";
      char *hostname;

      if ( (hostname = getenv("HOST")) == NULL ) hostname = unknown;

      setvbuf(sf, (char *)NULL, _IONBF, 0);
      fprintf(sf, "%.7s %4.4ld.%2.2ld.%2.2ld %2.2ld:%2.2ld %s "
	      "%-9.9s %7.1f %7.1f T%3.3d %s\n",
	      name,   yy,   mm,   dd,   hh,   mi,   VERSION,
	      proc, MaxMBytes, CPUTime, truncation, hostname);

      fclose(sf);
    }

#if defined (CRAY)
#  if defined (_CRAYMPP)
     strcpy(mtype, " CRAYMPP --");
#  elif (_MAXVL == 64)
     strcpy(mtype, " CRAYVL64 -");
#  elif (_MAXVL == 128)
     strcpy(mtype, " CRAYVL128 ");
#  else
     strcpy(mtype, " CRAY -----");
#  endif
#elif defined (SX)
     strcpy(mtype, " NECSX ----");
#elif defined (__uxp__)
     strcpy(mtype, " FUJI -----");
#elif defined (sun)
     strcpy(mtype, " SUN ------");
#elif defined (i386)
     strcpy(mtype, " i386 -----");
#elif defined (sgi)
     strcpy(mtype, " sgi ------");
#else
     strcpy(mtype, "-----------");
#endif

  fprintf(stdout, "   NORMAL EXIT\n");
  fprintf(stdout, " ------   End    after  -%-11.11s- %7.1f sec", mtype, CPUTime);
  if ( MaxMBytes > 0 )
    fprintf(stdout, " --- %7.1f MB ---\n", MaxMBytes);
  else
    fprintf(stdout, " ----------------\n");
}
#endif

static
void after_processing(struct Control *globs, struct Variable *vars)
{
  int i;

  //#if defined(_PSTREAM_H)
  //  globs->istreamID = streamOpenRead(cdoStreamName(0));
  //#else
  globs->istreamID = streamOpenRead(ifile);
  if ( globs->istreamID < 0 ) cdiError(globs->istreamID, "Open failed on %s", ifile);
  //#endif
  if ( ofiletype == -1 ) ofiletype = streamInqFiletype(globs->istreamID);

  globs->ivlistID = streamInqVlist(globs->istreamID);
  globs->taxisID  = vlistInqTaxis(globs->ivlistID);
  globs->taxisID2 = taxisDuplicate(globs->taxisID);

  if ( globs->Mean != 2 )
    {
#if defined(_PSTREAM_WRITE_H)
      globs->ostreamID = streamOpenWrite(cdoStreamName(ofileidx), ofiletype);
#else
      globs->ostreamID = streamOpenWrite(ofile, ofiletype);
      if ( globs->ostreamID < 0 ) cdiError(globs->ostreamID, "Open failed on %s", ofile);
#endif

      if ( globs->Szip ) streamDefCompType(globs->ostreamID, COMPRESS_SZIP);

      globs->ovlistID = vlistCreate();
    }
#if defined(AFTERBURNER)
  if ( globs->Mean >= 2 )
    {
      globs->ostreamID2 = streamOpenWrite(ofile2, ofiletype);
      if ( globs->ostreamID2 < 0 ) cdiError(globs->ostreamID2, "Open failed on %s", ofile2);

      if ( globs->Szip ) streamDefCompType(globs->ostreamID, COMPRESS_SZIP);

      globs->ovlistID2 = vlistCreate();
    }
#endif

  /* ---------------- */
  /*  pre-processing  */
  /* ---------------- */
  after_precntl(globs, vars);

  /* ----------------- */
  /*  initializations  */
  /* ----------------- */

  after_setCodes(globs, vars, MaxCodes, globs->NumCodesRequest);

  if ( globs->unitsel == 2 )
    for (i = 0; i < globs->NumLevelRequest; i++) globs->LevelRequest[i] = globs->LevelRequest[i] * 1000;

  if ( ! globs->AnalysisData )
    for (i = 0; i < globs->NumLevelRequest; i++)
      {
	if ( (globs->LevelRequest[i] >= 65535) && globs->unitsel && ofiletype == FILETYPE_GRB )
	  {
	    fprintf(stderr,"\n Level %9.2f out of range (max=65535)!\n", globs->LevelRequest[i]);
	    exit(1);
	  }

	if ( !globs->unitsel && globs->Type >= 20 && globs->NumLevelRequest > 1 && IS_EQUAL(globs->LevelRequest[i], 0))
	  {
	    fprintf(stderr,"\n Level %9.2f illegal for Type %d\n", globs->LevelRequest[i], globs->Type);
	    exit(1);
	  }
      }

  after_setLevel(globs);

  after_dimcalc(globs);

  globs->rcoslat          = (double *) malloc(globs->Latitudes*sizeof(double));
  globs->coslat           = (double *) malloc(globs->Latitudes*sizeof(double));
  globs->DerivationFactor = (double *) malloc(globs->Latitudes*sizeof(double));

  if ( globs->Type < 50 && globs->AnalysisData )
    {
      fprintf(stderr," ::::::::::::::::::::::::::::::::::::::::::::::\n");
      fprintf(stderr," -> Type < 50 is not appropriate for Analysis.\n");
      fprintf(stderr," -> Please check wether you can use Type >= 50.\n");
      fprintf(stderr," -> Premature Exit. Sorry.\n");
      exit(1);
    }

  if ( globs->Type == 10 || globs->Type == 40 || globs->Type == 60 )
    {
      if ( ofiletype == FILETYPE_GRB )
	Error("Can't write fourier coefficients to GRIB!");
      else if ( ofiletype == FILETYPE_NC || ofiletype == FILETYPE_NC2 ||
		ofiletype == FILETYPE_NC4 )
	Error("Can't write fourier coefficients to netCDF!");
    }

  filename = strrchr(ifile,'/');
  if (filename == 0) filename = ifile;
  else               filename++ ;

  if ( globs->Type >= 30 && globs->Type < 50 &&
      (vars[DIVERGENCE].selected || vars[VELOPOT].selected ||
       vars[VORTICITY].selected  || vars[STREAM].selected  ||
       globs->AnalysisData) )
    {
      /*
      int newtype = 0;
      */
      if ( globs->Type == 30 ) globs->Type = 70;
      if ( globs->Type == 40 ) globs->Type = 60;
      if ( globs->Type == 41 ) globs->Type = 61;

      if ( globs->AnalysisData )
	fprintf(stderr,"\n TYPE changed to %d (for analysis data)\n", globs->Type);
      else
	fprintf(stderr,"\n TYPE changed to %d (with code %d, %d, %d or %d)\n",
		globs->Type, DIVERGENCE, VELOPOT, VORTICITY, STREAM);
      /*
      if ( globs->Type == 30 ) newtype = 70;
      if ( globs->Type == 40 ) newtype = 60;
      if ( globs->Type == 41 ) newtype = 61;

      if ( globs->AnalysisData )
	fprintf(stderr,"\n Attention: TYPE isn't changed to %d anymore (for analysis data)!!!\n", globs->Type);
      else
	fprintf(stderr,"\n Attention: TYPE isn't changed to %d anymore (with code %d, %d, %d or %d)!!!\n",
		newtype, DIVERGENCE, VELOPOT, VORTICITY, STREAM);
      */
    }

  if ( globs->AnalysisData ) after_AnalysisDependencies(vars, MaxCodes);
  else
    {
      after_EchamDependencies(vars, MaxCodes, globs->Type, Source);
      vars[GEOPOTENTIAL].needed |= globs->Type >= 30 || vars[SLP].comp || vars[GEOPOTHEIGHT].comp;
    }

  /*  if ( vars[U_WIND].needed || vars[V_WIND].needed ) */
  if ( vars[U_WIND].comp || vars[V_WIND].comp )
    {
      globs->dv2uv_f1 = (double *) malloc(globs->DimSP_half*sizeof(double));
      globs->dv2uv_f2 = (double *) malloc(globs->DimSP_half*sizeof(double));
      geninx(globs->Truncation, globs->dv2uv_f1, globs->dv2uv_f2);
    }

  /* --------- */
  /*  Control  */
  /* --------- */

  after_defineLevel(globs, vars);

  after_defineGrid(globs, vars);

  after_postcntl(globs, vars); /* define output variables */

  after_control(globs, vars);

#if defined(_PSTREAM_WRITE_H)
  if ( globs->ostreamID  != CDI_UNDEFID ) pstreamClose(globs->ostreamID);
#else
  if ( globs->ostreamID2 != CDI_UNDEFID ) streamClose(globs->ostreamID2);
  if ( globs->ostreamID  != CDI_UNDEFID ) streamClose(globs->ostreamID);
#endif
#if defined(CDO)
  processDefVarNum(vlistNvars(globs->ivlistID), globs->istreamID);
  processAddNvals(streamNvals(globs->istreamID));
#endif
  streamClose(globs->istreamID);

  if ( globs->rcoslat )          free(globs->rcoslat);
  if ( globs->coslat )           free(globs->coslat);
  if ( globs->DerivationFactor ) free(globs->DerivationFactor);

  if ( globs->Field ) free(globs->Field);

  if ( globs->poli ) free(globs->poli);
  if ( globs->pold ) free(globs->pold);
  if ( globs->pdev ) free(globs->pdev);
  if ( globs->pol2 ) free(globs->pol2);  if ( globs->pol3 ) free(globs->pol3);
}

extern char *optarg;
extern int optind, opterr, optopt;

#if defined(AFTERBURNER)
static
int afterburner(int argc, char *argv[])
{
  int   i, code;
  char *proc = argv[0];
  char  Line[132];
  int c;
  int fargc0, fargcn;
  FILE *fp;
  int numThreads = 0;
  char *Vctfile = NULL;
  extern int dmemory_ExitOnError;

  dmemory_ExitOnError = 1;

  starttime = (double) clock();

#if defined(AFTERBURNER)
  { /* check character device on stdin and stdout */
    struct stat statbuf;
    fstat(0, &statbuf);
    if ( S_ISCHR(statbuf.st_mode) ) stdin_is_tty = 1;  
    fstat(1, &statbuf);
    if ( S_ISCHR(statbuf.st_mode) ) stdout_is_tty = 1;  
  }
#endif

  /* ------------------- */
  /*  print information  */
  /* ------------------- */

  lprintf(stdout);
  fprintf(stdout,"  afterburner version %s\n", VERSION);
  fprintf(stdout,"  ECHAM & analyses postprocessor\n");

  if ( sizeof(double) != 8 || sizeof(int) < 4 )
    {
      fprintf(stderr, "byte size of type double %d\n", (int) sizeof(double));
      fprintf(stderr, "byte size of type int %d\n",    (int) sizeof(int));
      fprintf(stderr, "byte size of type size_t %d\n", (int) sizeof(size_t));
      return(1);
    }

  fp = fopen("/pf/m/m214003/doc/afterburner.doc","r");
  if ( fp )
    {
      do
	{
	  fgets(Line, 130, fp);
	  fprintf(stdout, "%s", &Line[1]);
	}
      while ( ! feof(fp) && Line[0] == '#' );
      fclose(fp);
    }

  struct Control *globs = (struct Control *) malloc(sizeof(struct Control));
  after_control_init(globs);

  globs->Verbose = 1;

  /* --------------------- */
  /*  options & filenames  */
  /* --------------------- */
  extern int labort_after;

  while ( (c = getopt(argc, argv, "P:b:v:acdgpVw")) != EOF )
    switch (c)
      {
      case 'a': globs->AnalysisData = 1; break;
      case 'b': Message( "option -b not longer needed!\n"); break;
      case 'c': after_printCodes(); break;
      case 'd': globs->Debug = 1; break;
      case 'p': ParallelRead = TRUE; break;
      case 'P': numThreads = atoi(optarg); break;
      case 'V': after_version(); break;
      case 'v': Vctfile = optarg; break;
      case 'w': labort_after = FALSE; break;
      default:  Message( "option -%c unsupported!", optopt); after_usage();
      }

#if defined (_OPENMP)
  /* ParallelRead = TRUE; */

  lprintf(stdout);
  if ( numThreads <= 0 ) numThreads = 1;
  omp_set_num_threads(numThreads);
  if ( omp_get_max_threads() > omp_get_num_procs() )
    fprintf(stdout, " Number of threads is greater than number of Cores=%d!\n", omp_get_num_procs());
  fprintf(stdout, " OpenMP:  num_procs = %d  max_threads = %d\n", omp_get_num_procs(), omp_get_max_threads());
#else
  if ( numThreads > 0 )
    {
      fprintf(stderr, "Option -P failed, OpenMP support not compiled in!\n");
      return(-1);
    }
#endif

  if ( ParallelRead )
    {
#if  defined  (HAVE_LIBPTHREAD)
      fprintf(stdout, " Parallel read enabled\n");
#else
      fprintf(stdout, " Parallel read disabled\n");
      ParallelRead = FALSE;
#endif
    }

  fargc0 = optind;
  fargcn = argc;

  if ( optind < argc ) ifile = argv[optind++];
  if ( ! ifile )
    {
      Message( "*** Missing input file ***");
      after_usage();
    }

  struct Variable vars[MaxCodes+5];
  for ( code = 0; code < MaxCodes+5; code++ ) after_variable_init(&vars[code]);

  after_parini(globs, vars); /* read namelist parameter */

  fprintf(stdout, "   Input File: %-25s\n", ifile);
  if ( globs->Mean >= 2 )
    {
      if ( fargcn-fargc0 >= 3 ) ofile2 = argv[--fargcn];
      if ( ! ofile2 )
	{
	  Message( "*** Missing variance file ***");
	  after_usage();
	}
    }

  if ( globs->Mean != 2 )
    {
      if ( optind < argc ) ofile = argv[optind++];
      if ( fargcn-fargc0 >= 2 ) ofile = argv[--fargcn];
      if ( ! ofile )
	{
	  Message( "*** Missing output file ***");
	  after_usage();
	}
      fprintf(stdout, "  Output File: %-25s\n", ofile);
    }

  globs->Nfiles = fargcn-fargc0-1;
  if ( globs->Nfiles > 0 )
    {
      if ( globs->Multi > 0 )
	Error( "Namelist parameter MULTI works only with one inputfile");

      ifiles = (char **) malloc(globs->Nfiles*sizeof(char*));
      for ( i = 0; i < globs->Nfiles; i++ )
	ifiles[i] = argv[--fargcn];
    }

  if ( ofile2 )
    fprintf(stdout, "Variance File: %-25s\n", ofile2);

  if ( globs->Debug )
    {
      extern int afterDebug;
      afterDebug = globs->Debug;
      fprintf(stderr, "* Debug on!                              *\n");
      fprintf(stderr, "  Maximum ffts to run in parallel:  %ld\n", get_nfft());
    }

  /* read option VCT */
  if ( Vctfile ) after_readVct(globs, Vctfile);

  /* --------------------- */
  /*  open in/output file  */
  /* --------------------- */

  cdiDefGlobal("REGULARGRID", 1);

  after_processing(globs, vars);

  after_procstat(proc, globs->Truncation);

  FreeMean(vars);

  free(globs);

  return(0);
}
#endif

#if defined(CDO)
void *Afterburner(void *argument)
{
  cdoInitialize(argument);

  lstdout = !cdoSilentMode;

  struct Control *globs = (struct Control *) malloc(sizeof(struct Control));
  after_control_init(globs);

  globs->Verbose = cdoVerbose;

  struct Variable vars[MaxCodes+5];
  for ( int code = 0; code < MaxCodes+5; code++ ) after_variable_init(&vars[code]);

  after_parini(globs, vars); /* read namelist parameter */

  if ( cdoDefaultFileType != CDI_UNDEFID ) ofiletype = cdoDefaultFileType;

  int streamCnt = cdoStreamCnt();
  int nfiles = streamCnt - 1;

  ofileidx = nfiles;

  ifile = cdoStreamName(0)->args;
  ofile = cdoStreamName(nfiles)->args;

  globs->Nfiles = nfiles-1;
  if ( globs->Nfiles > 0 )
    {
      if ( globs->Multi > 0 )
	Error( "Namelist parameter MULTI works only with one inputfile");

      ifiles = (char **) malloc(globs->Nfiles*sizeof(char*));
      for ( int i = 0; i < globs->Nfiles; ++i )
	ifiles[i] = cdoStreamName(--nfiles)->args;
      for ( int i = 0; i < globs->Nfiles; ++i ) printf("files %d %s\n", i+1, ifiles[i]);
    }

  after_processing(globs, vars);

  FreeMean(vars);

  free(globs);

  cdoFinish();

  return (0);
}
#else
int main(int argc, char *argv[])
{
  return afterburner(argc, argv);
}
#endif
