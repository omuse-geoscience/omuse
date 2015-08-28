#if defined(HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <stdio.h>
#include <string.h>
#include <time.h>
/* #include <pwd.h> */
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>  /* write, close */
#include <fcntl.h>
#include <errno.h>
#include <stdlib.h>  /* qsort */

#include "cdo.h"
#include "process.h"

#if ! defined(VERSION)
#  define  VERSION  "0.0.1"
#endif

#define  MAX_LEN  65536

#define  LOGSSIZE  32
#define  LOGOSIZE  40

#undef HAVE_LOCK
#if defined(F_UNLCK) && defined(F_RDLCK) && defined(F_WRLCK)
#define HAVE_LOCK
#endif

#if defined(HAVE_LOCK)
static
void filestatus(struct flock *lock)
{
  switch(lock->l_type) {
    case F_UNLCK:
      printf("Status: F_UNLCK\n");
      break;
    case F_RDLCK:
      printf("Status: F_RDLCK (PID: %d)\n", lock->l_pid);
      break;
    case F_WRLCK:
      printf("Status: F_WRLCK (PID: %d)\n", lock->l_pid);
      break;
  }
}
#endif


void cdolog(const char *prompt, double cputime)
{
#if defined(HAVE_LOCK)
#if defined(LOGPATH)
#define  XSTRING(x)	#x
#define  STRING(x)	XSTRING(x)
  char logfilename[] = STRING(LOGPATH) "/cdo.log";
  int  logfileno;
  char *username;
  char timestr[30];
  time_t date_and_time_in_sec;
  struct tm *date_and_time;
  int streamID;
  const char *streamName;
  int len, slen, olen, pos;
  int loper;
  char logstring[MAX_LEN];
  int i;
  int status;
  struct flock mylock;

  memset(logstring, 0, MAX_LEN);

  date_and_time_in_sec = time(NULL);
  timestr[0] = 0;

  if ( date_and_time_in_sec != -1 )
    {
      date_and_time = localtime(&date_and_time_in_sec);
      (void) strftime(timestr, sizeof(timestr), "%d/%m/%Y %H:%M", date_and_time);
    }

  username = getenv("LOGNAME");
  if ( username == NULL )
    {
      username = getenv("USER");
      if ( username == NULL ) username = "unknown";
    }

  slen = sprintf(logstring, "%s %-8s %7.2f %s %-3s",
		 timestr,  username, cputime, VERSION, prompt);

  for ( streamID = 0; streamID < cdoStreamCnt(); streamID++ )
    {
      streamName = cdoStreamName(streamID)->args;
      pos = 0;
      while ( pos < (int) strlen(streamName) )
	{
	  len = 0;
	  loper = 0;
	  if ( streamName[pos++] == '-' ) loper = 1;
	  while ( streamName[pos+len] != ' ' &&  streamName[pos+len] != '\0' ) len++;
	  if ( len && loper )
	    {
	      for ( olen = 1; olen < len; olen++ )
		if ( streamName[pos+olen] == ',' ) break;

	      sprintf(logstring+slen, " %.*s", olen, &streamName[pos]);
	      slen += olen + 1;
	    }
	  pos += len + 1;
	}
    }

  sprintf(logstring+slen, "\n");
  slen++;

  errno = 0;
  logfileno = open(logfilename, O_WRONLY | O_APPEND);
  if ( errno )
    {
      errno = 0;
      return;
    }

  mylock.l_type   = F_WRLCK;
  mylock.l_whence = SEEK_SET;
  mylock.l_start  = 0;
  mylock.l_len    = 0;

  /*
  status = fcntl(logfileno, F_SETLKW, &mylock);
  */
  for ( i = 0; i < 100; i++ )
    {
      status = fcntl(logfileno, F_SETLK, &mylock);
      if ( status == 0 ) break;
      usleep(10000);
    }

  if ( status == 0 )
    {
      status = write(logfileno, logstring, slen);

      mylock.l_type   = F_UNLCK;
      status = fcntl(logfileno, F_SETLK, &mylock);
    }

  close(logfileno);

  errno = 0;

  return;

#endif
#endif
}


#include <math.h>
/*
 * convert an IMB float to single precision number v1.0
 *
 *                      Wesley Ebisuzaki
 */
static
float ibm2flt(unsigned char *ibm) {

	int positive, power;
	unsigned int abspower;
	long int mant;
	double value, exp;

	positive = (ibm[0] & 0x80) == 0;
	mant = (ibm[1] << 16) + (ibm[2] << 8) + ibm[3];
	power = (int) (ibm[0] & 0x7f) - 64;
	abspower = power > 0 ? power : -power;


	/* calc exp */
	exp = 16.0;
	value = 1.0;
	while (abspower) {
		if (abspower & 1) {
			value *= exp;
		}
		exp = exp * exp;
		abspower >>= 1;
	}

	if (power < 0) value = 1.0 / value;
	value = value * mant / 16777216.0;
	if (positive == 0) value = -value;
	return (float)value;
}

#if defined(HAVE_LOCK)
#if defined(LOGPATH)
/*
 * convert a float to an IBM single precision number v1.0
 *
 *                      Wesley Ebisuzaki
 *
 * doesn't handle subnormal numbers
 */

static
int flt2ibm(float x, unsigned char *ibm) {

	int sign, exp, i;
	double mant;

	if ( !(fabs((double)x) > 0) ) {
		ibm[0] = ibm[1] = ibm[2] = ibm[3] = 0;
		return 0;
	}

	/* sign bit */
	if (x < 0.0) {
		sign = 128;
		x = -x;
	}
	else sign = 0;

	mant = frexp((double) x, &exp);

	/* round up by adding 2**-24 */
	/* mant = mant + 1.0/16777216.0; */

	if (mant >= 1.0) {
		mant = 0.5;
		exp++;
	}
	while (exp & 3) {
		mant *= 0.5;
		exp++;
	}
	
	exp = exp/4 + 64;

	if (exp < 0) {
		fprintf(stderr,"underflow in flt2ibm\n");
		ibm[0] = ibm[1] = ibm[2] = ibm[3] = 0;
		return 0;
	}
	if (exp > 127) {
		fprintf(stderr,"overflow in flt2ibm\n");
		ibm[0] = sign | 127;
		ibm[1] = ibm[2] = ibm[3] = 255;
		return -1;
	}

	/* normal number */

	ibm[0] = sign | exp;

	mant = mant * 256.0;
	i = (int) floor(mant);
	mant = mant - i;
	ibm[1] = i;

	mant = mant * 256.0;
	i = (int) floor(mant);
	mant = mant - i;
	ibm[2] = i;

	ibm[3] = (int) floor(mant*256.0);

	return 0;
}
#endif
#endif

#define  GET_UINT4(xb)        ((int) (((int)xb[0]<<24) + \
                                      ((int)xb[1]<<16) + \
                                      ((int)xb[2]<<8)  + \
                                      ((int)xb[3])))
#define  GET_UINT8(xb)        ((int64_t) (((int64_t)xb[0]<<56) + \
                                         ((int64_t)xb[1]<<48) + \
                                         ((int64_t)xb[2]<<40) + \
                                         ((int64_t)xb[3]<<32) + \
                                         ((int64_t)xb[4]<<24) + \
                                         ((int64_t)xb[5]<<16) + \
                                         ((int64_t)xb[6]<<8)  + \
					 ((int64_t)xb[7])))

#define  PUT_UINT4(xb, iv)    ((*(xb)   = (iv) >> 24), \
                               (*(xb+1) = (iv) >> 16), \
                               (*(xb+2) = (iv) >>  8), \
                               (*(xb+3) = (iv)))
#define  PUT_UINT8(xb, iv)    ((*(xb)   = (iv) >> 56), \
                               (*(xb+1) = (iv) >> 48), \
                               (*(xb+2) = (iv) >> 40), \
                               (*(xb+3) = (iv) >> 32), \
                               (*(xb+4) = (iv) >> 24), \
                               (*(xb+5) = (iv) >> 16), \
                               (*(xb+6) = (iv) >>  8), \
                               (*(xb+7) = (iv)))


void cdologs(int noper)
{
#if defined(HAVE_LOCK)
#if defined(LOGPATH)
#define  XSTRING(x)	#x
#define  STRING(x)	XSTRING(x)
  char logfilename[] = STRING(LOGPATH) "/cdo.logs";
  int  logfileno;
  char timestr[30];
  time_t date_and_time_in_sec;
  struct tm *date_and_time;
  int i;
  int status;
  int uselock = 1;
  int date = 0, ncdo = 0, nhours = 0;
  int date0 = 0, ncdo0, noper0, nhours0;
  double cputime0;
  int64_t nvals0;
  unsigned char logbuf[LOGSSIZE];
  unsigned char *logdate   =  logbuf;
  unsigned char *logncdo   = &logbuf[4];
  unsigned char *lognoper  = &logbuf[8];
  unsigned char *logctime  = &logbuf[12];
  unsigned char *lognvals  = &logbuf[16];
  unsigned char *lognhours = &logbuf[24];
  const size_t logsize = LOGSSIZE;
  size_t bufsize;
  struct flock mylock;
  struct stat filestat;
  off_t nvals;
  int processID;
  double cputime;

  nvals = 0;
  cputime = 0;
  for ( processID = 0; processID < noper; processID++ )
    {
      cputime += processInqCputime(processID);
      nvals += processInqNvals(processID);
    }

  memset(logbuf, 0, logsize);

  date_and_time_in_sec = time(NULL);

  if ( date_and_time_in_sec != -1 )
    {
      date_and_time = localtime(&date_and_time_in_sec);
      (void) strftime(timestr, sizeof(timestr), "%Y%m%d", date_and_time);
      date = atoi(timestr);
    }

  errno = 0;
  logfileno = open(logfilename, O_RDWR);
  if ( errno )
    {
      errno = 0;
      return;
    }


  if ( (int) lseek(logfileno, (off_t) 0, SEEK_END) == 0 ) uselock = 0;
  status = (int) lseek(logfileno, 0, SEEK_SET);
  
  if ( uselock )
    {
      /*
	memset(&mylock, 0, sizeof(struct flock));
	status = fcntl(logfileno, F_GETLK, &mylock);
	filestatus(&mylock);
      */

      mylock.l_type   = F_WRLCK;
      mylock.l_whence = SEEK_SET;
      mylock.l_start  = 0;
      mylock.l_len    = 0;

      /*
	status = fcntl(logfileno, F_SETLKW, &mylock);
      */
      for ( i = 0; i < 100; i++ )
	{
	  status = fcntl(logfileno, F_SETLK, &mylock);
	  if ( status == 0 ) break;
	  usleep(10000);
	}
      errno = 0;
      if ( status != 0 ) goto endlabel;
    }

  status = fstat(logfileno, &filestat);
  errno = 0;
  if ( status != 0 ) goto endlabel;

  bufsize = (size_t) filestat.st_size;

  if ( bufsize > 0 )
    {
      status = (int) lseek(logfileno, (off_t) (bufsize-logsize), SEEK_SET);
      status = (int) read(logfileno, logbuf, logsize);

      date0    = GET_UINT4(logdate);
      ncdo0    = GET_UINT4(logncdo);
      noper0   = GET_UINT4(lognoper);
      cputime0 = (double) ibm2flt(logctime);
      nvals0   = GET_UINT8(lognvals);
      nhours0  = GET_UINT4(lognhours);

      if ( date == date0 )
	{
	  ncdo = ncdo0;
	  nhours = nhours0;
	  noper += noper0;
	  cputime += cputime0;
	  nvals += (off_t) nvals0;
	  status = (int) lseek(logfileno, (off_t) (bufsize-logsize), SEEK_SET);
	}
    }

  if ( date >= date0 )
    {
      ncdo++;

      while ( cputime >= 3600 ) { cputime -= 3600; nhours++; }

      PUT_UINT4(logdate, date);
      PUT_UINT4(logncdo, ncdo);
      PUT_UINT4(lognoper, noper);
      flt2ibm((float)cputime, logctime);
      PUT_UINT8(lognvals, nvals);
      PUT_UINT4(lognhours, nhours);

      status = (int) write(logfileno, logbuf, logsize);
    }

 endlabel:

  if ( uselock )
    {
      mylock.l_type = F_UNLCK;
      status = fcntl(logfileno, F_SETLK, &mylock);
    }

  close(logfileno);

  errno = 0;

  return;

#endif
#endif
}


void dumplogs(const char *logfilename)
{
#if defined(HAVE_LOCK)
  int  logfileno;
  int status;
  int date0 = 0, ncdo0, noper0, nhours0;
  int nlogs;
  int i;
  double cputime0;
  int64_t nvals0;
  unsigned char logbuf[LOGSSIZE];
  unsigned char *logdate   =  logbuf;
  unsigned char *logncdo   = &logbuf[4];
  unsigned char *lognoper  = &logbuf[8];
  unsigned char *logctime  = &logbuf[12];
  unsigned char *lognvals  = &logbuf[16];
  unsigned char *lognhours = &logbuf[24];
  unsigned char *buffer = NULL;
  const size_t logsize = LOGSSIZE;
  size_t bufsize;
  struct flock mylock;
  struct stat filestat;

  errno = 0;
  logfileno = open(logfilename, O_RDONLY);
  if ( errno )
    {
      cdoAbort("Open failed on %s", logfilename);
      errno = 0;
      return;
    }

  memset(&mylock, 0, sizeof(struct flock));
  status = fcntl(logfileno, F_GETLK, &mylock);
  filestatus(&mylock);

  status = fstat(logfileno, &filestat);
  errno = 0;
  if ( status != 0 ) return;

  bufsize = (size_t) filestat.st_size;

  if ( bufsize > 0 )
    {
      buffer = (unsigned char*) malloc(bufsize);

      status = (int) read(logfileno, buffer, bufsize);

      nlogs = bufsize / logsize;
      for ( i = 0; i < nlogs; i++ )
	{
	  memcpy(logbuf, &buffer[i*logsize], logsize);
	  date0    = GET_UINT4(logdate);
	  ncdo0    = GET_UINT4(logncdo);
	  noper0   = GET_UINT4(lognoper);
	  cputime0 = (double) ibm2flt(logctime);
	  nvals0   = GET_UINT8(lognvals);
	  nhours0  = GET_UINT4(lognhours);

	  if ( sizeof(int64_t) > sizeof(long) )
	    fprintf(stdout, "%8d %10d %10d %19lld %10d %8.2f\n",
		    date0, ncdo0, noper0, (long long)nvals0, nhours0, cputime0);
	  else
	    fprintf(stdout, "%8d %10d %10d %19ld %10d %8.2f\n",
		    date0, ncdo0, noper0, (long)nvals0, nhours0, cputime0);
	}

      free(buffer);
    }

  close(logfileno);

  errno = 0;

  return;
#endif
}


void daylogs(const char *logfilename)
{
  int  logfileno;
  int status;
  int date0 = 0, ncdo0, noper0, nhours0;
  int nlogs;
  int i;
  double cputime0;
  int64_t nvals0;
  unsigned char logbuf[LOGSSIZE];
  unsigned char *logdate   =  logbuf;
  unsigned char *logncdo   = &logbuf[4];
  unsigned char *lognoper  = &logbuf[8];
  unsigned char *logctime  = &logbuf[12];
  unsigned char *lognvals  = &logbuf[16];
  unsigned char *lognhours = &logbuf[24];
  unsigned char *buffer = NULL;
  const size_t logsize = LOGSSIZE;
  size_t bufsize;
  struct stat filestat;

  errno = 0;
  logfileno = open(logfilename, O_RDONLY);
  if ( errno )
    {
      cdoAbort("Open failed on %s", logfilename);
      errno = 0;
      return;
    }

  status = fstat(logfileno, &filestat);
  errno = 0;
  if ( status != 0 ) return;

  bufsize = (size_t) filestat.st_size;

  if ( bufsize > 0 )
    {
      buffer = (unsigned char*) malloc(bufsize);

      status = (int) read(logfileno, buffer, bufsize);

      fprintf(stdout, "# day           noper         size [MB]    time [s]\n");
      nlogs = bufsize / logsize;
      for ( i = 0; i < nlogs; i++ )
	{
	  memcpy(logbuf, &buffer[i*logsize], logsize);
	  date0    = GET_UINT4(logdate);
	  ncdo0    = GET_UINT4(logncdo);
	  noper0   = GET_UINT4(lognoper);
	  cputime0 = (double) ibm2flt(logctime);
	  nvals0   = GET_UINT8(lognvals);
	  nhours0  = GET_UINT4(lognhours);

	  fprintf(stdout, "%8d %12d %12d %12d\n",
		  date0, noper0, (int)(8*nvals0/(1024*1024)), (int) (nhours0*3600.+cputime0));
	}

      free(buffer);
    }

  close(logfileno);

  errno = 0;

  return;
}


void monlogs(const char *logfilename)
{
  int  logfileno;
  int status;
  int date0 = 0, ncdo0, noper0, nhours0;
  int nlogs;
  int i;
  double cputime0;
  int64_t nvals0;
  unsigned char logbuf[LOGSSIZE];
  unsigned char *logdate   =  logbuf;
  unsigned char *logncdo   = &logbuf[4];
  unsigned char *lognoper  = &logbuf[8];
  unsigned char *logctime  = &logbuf[12];
  unsigned char *lognvals  = &logbuf[16];
  unsigned char *lognhours = &logbuf[24];
  unsigned char *buffer = NULL;
  int ymon = 0, ymon0 = 0;
  unsigned int noper = 0;
  double size = 0, cputime = 0;
  const size_t logsize = LOGSSIZE;
  size_t bufsize;
  struct stat filestat;

  errno = 0;
  logfileno = open(logfilename, O_RDONLY);
  if ( errno )
    {
      cdoAbort("Open failed on %s", logfilename);
      errno = 0;
      return;
    }

  status = fstat(logfileno, &filestat);
  errno = 0;
  if ( status != 0 ) return;

  bufsize = (size_t) filestat.st_size;

  if ( bufsize > 0 )
    {
      buffer = (unsigned char*) malloc(bufsize);

      status = (int) read(logfileno, buffer, bufsize);

      fprintf(stdout, "# month         noper         size [GB]    time [h]\n");
      nlogs = bufsize / logsize;
      for ( i = 0; i < nlogs; i++ )
	{
	  memcpy(logbuf, &buffer[i*logsize], logsize);
	  date0    = GET_UINT4(logdate);
	  ncdo0    = GET_UINT4(logncdo);
	  noper0   = GET_UINT4(lognoper);
	  cputime0 = (double) ibm2flt(logctime);
	  nvals0   = GET_UINT8(lognvals);
	  nhours0  = GET_UINT4(lognhours);

	  if ( i == 0 ) ymon0 = date0/100;
	  ymon = date0/100;

	  if ( ymon != ymon0 )
	    {
	      if ( i )
		fprintf(stdout, "%6d   %12u %12d %12d\n",
			ymon0, noper, (int)(8*size/(1024*1024*1024)), (int) (cputime/3600));
		
	      noper = 0;
	      size  = 0;
              cputime = 0;
	    }
	  noper += noper0;
	  size += nvals0;
	  cputime += (nhours0*3600.+cputime0);
	  
	  ymon0 = ymon;
	}
      fprintf(stdout, "%6d   %12u %12d %12d\n",
	      ymon0, noper, (int)(8*size/(1024*1024*1024)), (int) (cputime/3600));
      
      free(buffer);
    }

  close(logfileno);

  errno = 0;

  return;
}


void cdologo(int noper)
{
#if defined(HAVE_LOCK)
#if defined(LOGPATH)
#define  XSTRING(x)	#x
#define  STRING(x)	XSTRING(x)
  char logfilename[] = STRING(LOGPATH) "/cdo.logo";
  int  logfileno;
  int i;
  int status;
  int uselock = 1;
  int nhours = 0;
  int nhours0 = 0;
  double cputime0 = 0;
  int64_t nvals0 = 0;
  unsigned char logbuf[LOGOSIZE];
  unsigned char *logname   =  logbuf;
  unsigned char *lognocc   = &logbuf[16];
  unsigned char *lognvals  = &logbuf[20];
  unsigned char *lognhours = &logbuf[28];
  unsigned char *logctime  = &logbuf[32];
  unsigned char *buffer = NULL;
  const size_t logsize = LOGOSIZE;
  size_t bufsize, newbufsize;
  struct flock mylock;
  struct stat filestat;
  int nocc = 0;
  int nbuf;
  int processID;
  off_t onvals[256];
  double ocputime[256];
  const char *oname[256];

  for ( processID = 0; processID < noper; processID++ )
    {
      ocputime[processID] = processInqCputime(processID);
      onvals[processID]   = processInqNvals(processID);
      oname[processID]    = processInqOpername2(processID);
    }

  memset(logbuf, 0, logsize);

  errno = 0;
  logfileno = open(logfilename, O_RDWR);
  if ( errno )
    {
      errno = 0;
      return;
    }

  if ( (int) lseek(logfileno, (off_t) 0, SEEK_END) == 0 ) uselock = 0;
  status = (int) lseek(logfileno, 0, SEEK_SET);
  
  if ( uselock )
    {
      mylock.l_type   = F_WRLCK;
      mylock.l_whence = SEEK_SET;
      mylock.l_start  = 0;
      mylock.l_len    = 0;

      /*
	status = fcntl(logfileno, F_SETLKW, &mylock);
      */
      for ( i = 0; i < 100; i++ )
	{
	  status = fcntl(logfileno, F_SETLK, &mylock);
	  if ( status == 0 ) break;
	  usleep(10000);
	}
      errno = 0;
      if ( status != 0 ) goto endlabel;
    }

  status = fstat(logfileno, &filestat);
  errno = 0;
  if ( status != 0 ) goto endlabel;

  bufsize = (size_t) filestat.st_size;

  newbufsize = bufsize + noper*logsize;
  buffer = (unsigned char*) malloc(newbufsize);
  if ( bufsize > 0 )
    status = (int) read(logfileno, buffer, bufsize);

  for ( processID = 0; processID < noper; processID++ )
    {
      if ( oname[processID] == NULL ) break; /* baustelle !! */
      nbuf = (int) (bufsize / logsize);

      for ( i = 0; i < nbuf; i++ )
	{
	  memcpy(logbuf, buffer+i*logsize, logsize);
	  nocc     = GET_UINT4(lognocc);
	  nvals0   = GET_UINT8(lognvals);
	  nhours0  = GET_UINT4(lognhours);
	  cputime0 = (double) ibm2flt(logctime);

	  if ( strcmp((const char*)logname, oname[processID]) == 0 ) break;
	}

      if ( i == nbuf )
	{
	  nocc     = 0;
	  nvals0   = 0;
          nhours0  = 0;
          cputime0 = 0;
	  bufsize += logsize;
	  strcpy((char *)logname, oname[processID]);
	}

      nocc++;
	  
      nvals0   += (off_t) onvals[processID];
      nhours0  += nhours;
      cputime0 += ocputime[processID];
      while ( cputime0 >= 3600 ) { cputime0 -= 3600; nhours0++; }

      PUT_UINT4(lognocc, nocc);
      PUT_UINT8(lognvals, nvals0);
      PUT_UINT4(lognhours, nhours0);
      flt2ibm((float)cputime0, logctime);

      memcpy(buffer+i*logsize, logbuf, logsize);
    }

  status = (int) lseek(logfileno, 0, SEEK_SET);
  status = (int) write(logfileno, buffer, bufsize);

  free(buffer);

 endlabel:

  if ( uselock )
    {
      mylock.l_type = F_UNLCK;
      status = fcntl(logfileno, F_SETLK, &mylock);
    }

  close(logfileno);

  errno = 0;

  return;

#endif
#endif
}


typedef struct
{
  int      nocc;
  int64_t  nvals;
  double   time;
  double   perc;
  char     name[128];
}
loginfo_t;

static
int cmplognocc(const void *s1, const void *s2)
{
  int cmp = 0;
  const loginfo_t *x = s1;
  const loginfo_t *y = s2;

  if      ( x->nocc < y->nocc ) cmp =  1;
  else if ( x->nocc > y->nocc ) cmp = -1;

  return (cmp);
}

static
int cmplognvals(const void *s1, const void *s2)
{
  int cmp = 0;
  const loginfo_t *x = s1;
  const loginfo_t *y = s2;

  if      ( x->nvals < y->nvals ) cmp =  1;
  else if ( x->nvals > y->nvals ) cmp = -1;

  return (cmp);
}

static
int cmplogtime(const void *s1, const void *s2)
{
  int cmp = 0;
  const loginfo_t *x = s1;
  const loginfo_t *y = s2;

  if      ( x->time < y->time ) cmp =  1;
  else if ( x->time > y->time ) cmp = -1;

  return (cmp);
}

static
int cmplogperc(const void *s1, const void *s2)
{
  int cmp = 0;
  const loginfo_t *x = s1;
  const loginfo_t *y = s2;

  if      ( x->perc < y->perc ) cmp =  1;
  else if ( x->perc > y->perc ) cmp = -1;

  return (cmp);
}

static
int cmplogname(const void *s1, const void *s2)
{
  const loginfo_t *x = s1;
  const loginfo_t *y = s2;

  return (strcmp(x->name, y->name));
}


void dumplogo(const char *logfilename, int dumptype)
{
#if defined(HAVE_LOCK)
  int logfileno;
  int status;
  int nocc;
  int nhours0;
  int nlogs;
  int i;
  int mem;
  double cputime0;
  int64_t nvals0;
  unsigned char logbuf[LOGOSIZE];
  unsigned char *logname   =  logbuf;
  unsigned char *lognocc   = &logbuf[16];
  unsigned char *lognvals  = &logbuf[20];
  unsigned char *lognhours = &logbuf[28];
  unsigned char *logctime  = &logbuf[32];
  unsigned char *buffer = NULL;
  const size_t logsize = LOGOSIZE;
  size_t bufsize;
  struct flock mylock;
  struct stat filestat;
  loginfo_t **logInfo;

  errno = 0;
  logfileno = open(logfilename, O_RDONLY);
  if ( errno )
    {
      cdoAbort("Open failed on %s", logfilename);
      errno = 0;
      return;
    }

  memset(&mylock, 0, sizeof(struct flock));
  status = fcntl(logfileno, F_GETLK, &mylock);
  filestatus(&mylock);

  status = fstat(logfileno, &filestat);
  errno = 0;
  if ( status != 0 ) return;

  bufsize = (size_t) filestat.st_size;

  if ( bufsize > 0 )
    {
      fprintf(stdout, "# num name                     call        mem [GB]    time [h]     perc [s]\n");
      buffer = (unsigned char*) malloc(bufsize);

      status = (int) read(logfileno, buffer, bufsize);

      nlogs = bufsize / logsize;

      logInfo    = (loginfo_t **) malloc(nlogs*sizeof(loginfo_t *));
      logInfo[0] = (loginfo_t*) malloc(nlogs*sizeof(loginfo_t));
      for ( i = 1; i < nlogs; i++ ) logInfo[i] = logInfo[0] + i;

      for ( i = 0; i < nlogs; i++ )
	{
	  memcpy(logbuf, &buffer[i*logsize], logsize);
	  nocc     = GET_UINT4(lognocc);
	  nvals0   = GET_UINT8(lognvals);
	  nhours0  = GET_UINT4(lognhours);
	  cputime0 = (double) ibm2flt(logctime);

	  strcpy(logInfo[i]->name, (const char *)logname);
	  logInfo[i]->nocc  = nocc;
	  logInfo[i]->nvals = nvals0;
	  logInfo[i]->time  = nhours0 + cputime0/3600;
	  logInfo[i]->perc  = 3600*logInfo[i]->time/nocc;
	}

      if      ( dumptype == 1 )
	qsort(logInfo[0], nlogs, sizeof(loginfo_t), cmplogname);
      else if ( dumptype == 2 )
	qsort(logInfo[0], nlogs, sizeof(loginfo_t), cmplognocc);
      else if ( dumptype == 3 )
	qsort(logInfo[0], nlogs, sizeof(loginfo_t), cmplognvals);
      else if ( dumptype == 4 )
	qsort(logInfo[0], nlogs, sizeof(loginfo_t), cmplogtime);
      else if ( dumptype == 5 )
	qsort(logInfo[0], nlogs, sizeof(loginfo_t), cmplogperc);

      for ( i = 0; i < nlogs; i++ )
	{
	  mem = (int)(8*logInfo[i]->nvals/(1024*1024*1024));
	  if ( sizeof(int64_t) > sizeof(long) )
	    fprintf(stdout, "%4d  %-16s %12d %12d %12.3f %12.3f\n", i+1, logInfo[i]->name, 
		    logInfo[i]->nocc, mem, logInfo[i]->time, logInfo[i]->perc);
	  else
	    fprintf(stdout, "%4d  %-16s %12d %12d %12.3f %12.3f\n", i+1, logInfo[i]->name,
		    logInfo[i]->nocc, mem, logInfo[i]->time, logInfo[i]->perc);
	}

      free(logInfo[0]);
      free(logInfo);
      free(buffer);
    }

  close(logfileno);

  errno = 0;

  return;
#endif
}
