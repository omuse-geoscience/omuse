#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#ifdef HAVE_LIBNETCDF

//#define TEST_GROUPS 1

#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>

#include <netcdf.h>

#include "dmemory.h"
#include "cdi.h"
#include "basetime.h"
#include "gaussgrid.h"
#include "cdi_int.h"
#include "stream_cdf.h"
#include "cdf_int.h"
#include "varscan.h"
#include "vlist.h"

//#define PROJECTION_TEST

#undef  UNDEFID
#define UNDEFID  CDI_UNDEFID


#if defined HAVE_LIBNETCDF
static void cdfDefGlobalAtts(stream_t *streamptr);
static void cdfDefLocalAtts(stream_t *streamptr);
#endif

#define  X_AXIS  1
#define  Y_AXIS  2
#define  Z_AXIS  3
#define  T_AXIS  4

#define  POSITIVE_UP    1
#define  POSITIVE_DOWN  2

typedef struct {
  int     ncvarid;
  int     dimtype;
  size_t  len;
  char    name[CDI_MAX_NAME];
}
ncdim_t;

#define  MAX_COORDVARS  4
#define  MAX_AUXVARS    4

typedef struct {
  int      ncid;
  int      ignore;
  int      isvar;
  int      islon;
  int      islat;
  int      islev;
  int      istime;
  int      warn;
  int      tsteptype;
  int      param;
  int      code;
  int      tabnum;
  int      climatology;
  int      bounds;
  int      gridID;
  int      zaxisID;
  int      gridtype;
  int      zaxistype;
  int      xdim;
  int      ydim;
  int      zdim;
  int      xvarid;
  int      yvarid;
  int      zvarid;
  int      tvarid;
  int      ncoordvars;
  int      coordvarids[MAX_COORDVARS];
  int      nauxvars;
  int      auxvarids[MAX_AUXVARS];
  int      cellarea;
  int      calendar;
  int      tableID;
  int      truncation;
  int      position;
  int      defmissval;
  int      deffillval;
  int      xtype;
  int      ndims;
  int      gmapid;
  int      positive;
  int      dimids[8];
  int      dimtype[8];
  int      chunks[8];
  int      chunked;
  int      chunktype;
  int      natts;
  int     *atts;
  int      deflate;
  int      lunsigned;
  int      lvalidrange;
  size_t   vlen;
  double  *vdata;
  double   missval;
  double   fillval;
  double   addoffset;
  double   scalefactor;
  double   validrange[2];
  char     name[CDI_MAX_NAME];
  char     longname[CDI_MAX_NAME];
  char     stdname[CDI_MAX_NAME];
  char     units[CDI_MAX_NAME];
  char     extra[CDI_MAX_NAME];
  ensinfo_t   *ensdata;    /* Ensemble information */
}
ncvar_t;

static
void strtolower(char *str)
{
  if ( str )
    for (size_t i = 0; str[i]; ++i)
      str[i] = (char)tolower((int)str[i]);
}

static
int get_timeunit(size_t len, const char *ptu)
{
  int timeunit = -1;

  if ( len > 2 )
    {
      if      ( memcmp(ptu, "sec",    3) == 0 )          timeunit = TUNIT_SECOND;
      else if ( memcmp(ptu, "minute", 6) == 0 )          timeunit = TUNIT_MINUTE;
      else if ( memcmp(ptu, "hour",   4) == 0 )          timeunit = TUNIT_HOUR;
      else if ( memcmp(ptu, "day",    3) == 0 )          timeunit = TUNIT_DAY;
      else if ( memcmp(ptu, "month",  5) == 0 )          timeunit = TUNIT_MONTH;
      else if ( memcmp(ptu, "calendar_month", 14) == 0 ) timeunit = TUNIT_MONTH;
      else if ( memcmp(ptu, "year",   4) == 0 )          timeunit = TUNIT_YEAR;
    }
  else if ( len == 1 )
    {
      if ( ptu[0] == 's' ) timeunit = TUNIT_SECOND;
    }

  return (timeunit);
}

static
int isTimeUnits(const char *timeunits)
{
  int status = 0;

  if ( strncmp(timeunits, "sec",    3) == 0 ||
       strncmp(timeunits, "minute", 6) == 0 ||
       strncmp(timeunits, "hour",   4) == 0 ||
       strncmp(timeunits, "day",    3) == 0 ||
       strncmp(timeunits, "month",  5) == 0 ) status = 1;

  return (status);
}

static
int isTimeAxisUnits(const char *timeunits)
{
  char *ptu, *tu;
  int timetype = -1;
  int timeunit;
  int status = FALSE;

  size_t len = strlen(timeunits);
  tu = xmalloc((len+1)*sizeof(char));
  memcpy(tu, timeunits, (len+1) * sizeof(char));
  ptu = tu;

  for (size_t i = 0; i < len; i++ ) ptu[i] = (char)tolower((int)ptu[i]);

  timeunit = get_timeunit(len, ptu);
  if ( timeunit != -1 )
    {

      while ( ! isspace(*ptu) && *ptu != 0 ) ptu++;
      if ( *ptu )
        {
          while ( isspace(*ptu) ) ptu++;

          if ( memcmp(ptu, "as", 2) == 0 )
            timetype = TAXIS_ABSOLUTE;
          else if ( memcmp(ptu, "since", 5) == 0 )
            timetype = TAXIS_RELATIVE;

          if ( timetype != -1 ) status = TRUE;
        }
    }

  free(tu);

  return (status);
}

static
void scanTimeString(const char *ptu, int *rdate, int *rtime)
{
  int year = 1, month = 1, day = 1;
  int hour = 0, minute = 0, second = 0;
  int v1 = 1, v2 = 1, v3 = 1;

  *rdate = 0;
  *rtime = 0;

  if ( *ptu )
    {
      v1 = atoi(ptu);
      if ( v1 < 0 ) ptu++;
      while ( isdigit((int) *ptu) ) ptu++;
      if ( *ptu )
        {
          v2 = atoi(++ptu);
          while ( isdigit((int) *ptu) ) ptu++;
          if ( *ptu )
            {
              v3 = atoi(++ptu);
              while ( isdigit((int) *ptu) ) ptu++;
            }
        }
    }

  if ( v3 > 999 && v1 < 32 )
    { year = v3; month = v2; day = v1; }
  else
    { year = v1; month = v2; day = v3; }

  while ( isspace((int) *ptu) ) ptu++;

  if ( *ptu )
    {
      while ( ! isdigit((int) *ptu) ) ptu++;

      hour = atoi(ptu);
      while ( isdigit((int) *ptu) ) ptu++;
      if ( *ptu == ':' )
        {
          ptu++;
          minute = atoi(ptu);
          while ( isdigit((int) *ptu) ) ptu++;
          if ( *ptu == ':' )
            {
              ptu++;
              second = atoi(ptu);
            }
        }
    }

  *rdate = cdiEncodeDate(year, month, day);
  *rtime = cdiEncodeTime(hour, minute, second);
}

static
int scanTimeUnit(const char *unitstr)
{
  int timeunit = -1;
  size_t len = strlen(unitstr);
  timeunit = get_timeunit(len, unitstr);
  if ( timeunit == -1 )
    Message("Unsupported TIMEUNIT: %s!", unitstr);

  return (timeunit);
}

static
void setForecastTime(const char *timestr, taxis_t *taxis)
{
  int len;

  (*taxis).fdate = 0;
  (*taxis).ftime = 0;

  len = (int) strlen(timestr);
  if ( len == 0 ) return;

  int fdate = 0, ftime = 0;
  scanTimeString(timestr, &fdate, &ftime);

  (*taxis).fdate = fdate;
  (*taxis).ftime = ftime;
}

static
int setBaseTime(const char *timeunits, taxis_t *taxis)
{
  char *ptu, *tu;
  int timetype = TAXIS_ABSOLUTE;
  int rdate = -1, rtime = -1;
  int timeunit;

  size_t len = strlen(timeunits);
  tu = xmalloc((len+1) * sizeof (char));
  memcpy(tu, timeunits, (len+1) * sizeof (char));
  ptu = tu;

  for ( size_t i = 0; i < len; i++ ) ptu[i] = (char)tolower((int) ptu[i]);

  timeunit = get_timeunit(len, ptu);
  if ( timeunit == -1 )
    {
      Message("Unsupported TIMEUNIT: %s!", timeunits);
      return (1);
    }

  while ( ! isspace(*ptu) && *ptu != 0 ) ptu++;
  if ( *ptu )
    {
      while ( isspace(*ptu) ) ptu++;

      if ( memcmp(ptu, "as", 2) == 0 )
        timetype = TAXIS_ABSOLUTE;
      else if ( memcmp(ptu, "since", 5) == 0 )
        timetype = TAXIS_RELATIVE;

      while ( ! isspace(*ptu) && *ptu != 0 ) ptu++;
      if ( *ptu )
        {
          while ( isspace(*ptu) ) ptu++;

          if ( timetype == TAXIS_ABSOLUTE )
            {
              if ( memcmp(ptu, "%y%m%d.%f", 9) != 0 && timeunit == TUNIT_DAY )
                {
                  Message("Unsupported format %s for TIMEUNIT day!", ptu);
                  timeunit = -1;
                }
              else if ( memcmp(ptu, "%y%m.%f", 7) != 0 && timeunit == TUNIT_MONTH )
                {
                  Message("Unsupported format %s for TIMEUNIT month!", ptu);
                  timeunit = -1;
                }
            }
          else if ( timetype == TAXIS_RELATIVE )
            {
              scanTimeString(ptu, &rdate, &rtime);

              (*taxis).rdate = rdate;
              (*taxis).rtime = rtime;

              if ( CDI_Debug )
                Message("rdate = %d  rtime = %d", rdate, rtime);
            }
        }
    }

  (*taxis).type = timetype;
  (*taxis).unit = timeunit;

  free(tu);

  if ( CDI_Debug )
    Message("timetype = %d  unit = %d", timetype, timeunit);

  return (0);
}

static
void cdfGetAttInt(int fileID, int ncvarid, char *attname, int attlen, int *attint)
{
  nc_type atttype;
  size_t nc_attlen;

  *attint = 0;

  cdf_inq_atttype(fileID, ncvarid, attname, &atttype);
  cdf_inq_attlen(fileID, ncvarid, attname, &nc_attlen);

  if ( atttype != NC_CHAR )
    {
      int *pintatt = NULL;

      if ( (int)nc_attlen > attlen )
        pintatt = (int *) malloc(nc_attlen * sizeof (int));
      else
        pintatt = attint;

      cdf_get_att_int(fileID, ncvarid, attname, pintatt);

      if ( (int)nc_attlen > attlen )
        {
          memcpy(attint, pintatt, (size_t)attlen * sizeof (int));
          free(pintatt);
        }
    }
}

static
void cdfGetAttDouble(int fileID, int ncvarid, char *attname, int attlen, double *attdouble)
{
  nc_type atttype;
  size_t nc_attlen;

  *attdouble = 0;

  cdf_inq_atttype(fileID, ncvarid, attname, &atttype);
  cdf_inq_attlen(fileID, ncvarid, attname, &nc_attlen);

  if ( atttype != NC_CHAR )
    {
      double *pdoubleatt = NULL;

      if ( (int)nc_attlen > attlen )
        pdoubleatt = (double *) malloc(nc_attlen * sizeof (double));
      else
        pdoubleatt = attdouble;

      cdf_get_att_double(fileID, ncvarid, attname, pdoubleatt);

      if ( (int)nc_attlen > attlen )
        {
          memcpy(attdouble, pdoubleatt, (size_t)attlen * sizeof (double));
          free(pdoubleatt);
        }
    }
}

static
void cdfGetAttText(int fileID, int ncvarid, char *attname, int attlen, char *atttext)
{
  nc_type atttype;
  size_t nc_attlen;

  cdf_inq_atttype(fileID, ncvarid, attname, &atttype);
  cdf_inq_attlen(fileID, ncvarid, attname, &nc_attlen);

  if ( atttype == NC_CHAR )
    {
      char attbuf[65636];
      if ( nc_attlen < sizeof(attbuf) )
        {
          cdf_get_att_text(fileID, ncvarid, attname, attbuf);

          if ( (int) nc_attlen > (attlen-1) ) nc_attlen = (size_t)(attlen-1);

          attbuf[nc_attlen++] = 0;
          memcpy(atttext, attbuf, nc_attlen);
        }
      else
        {
          atttext[0] = 0;
        }
    }
#if  defined  (HAVE_NETCDF4)
  else if ( atttype == NC_STRING )
    {
      if ( nc_attlen == 1 )
        {
          char *attbuf = NULL;
          cdf_get_att_string(fileID, ncvarid, attname, &attbuf);

          size_t ssize = strlen(attbuf) + 1;

          if ( ssize > (size_t)attlen ) ssize = (size_t)attlen;
          memcpy(atttext, attbuf, ssize);
          atttext[ssize - 1] = 0;
          free(attbuf);
        }
      else
        {
          atttext[0] = 0;
        }
    }
#endif
}

static
int xtypeIsText(int xtype)
{
  int isText = FALSE;

  if ( xtype == NC_CHAR )
    isText = TRUE;
#if  defined  (HAVE_NETCDF4)
  else if ( xtype == NC_STRING )
    isText = TRUE;
#endif

  return (isText);
}

static
int xtypeIsFloat(int xtype)
{
  int isFloat = FALSE;

  if ( xtype == NC_FLOAT || xtype == NC_DOUBLE ) isFloat = TRUE;

  return isFloat;
}

static
int cdfInqDatatype(int xtype, int lunsigned)
{
  int datatype = -1;

#if  defined  (HAVE_NETCDF4)
  if ( xtype == NC_BYTE && lunsigned ) xtype = NC_UBYTE;
#endif

  if      ( xtype == NC_BYTE   )  datatype = DATATYPE_INT8;
  /* else if ( xtype == NC_CHAR   )  datatype = DATATYPE_UINT8; */
  else if ( xtype == NC_SHORT  )  datatype = DATATYPE_INT16;
  else if ( xtype == NC_INT    )  datatype = DATATYPE_INT32;
  else if ( xtype == NC_FLOAT  )  datatype = DATATYPE_FLT32;
  else if ( xtype == NC_DOUBLE )  datatype = DATATYPE_FLT64;
#if  defined  (HAVE_NETCDF4)
  else if ( xtype == NC_UBYTE  )  datatype = DATATYPE_UINT8;
  else if ( xtype == NC_LONG   )  datatype = DATATYPE_INT32;
  else if ( xtype == NC_USHORT )  datatype = DATATYPE_UINT16;
  else if ( xtype == NC_UINT   )  datatype = DATATYPE_UINT32;
  else if ( xtype == NC_INT64  )  datatype = DATATYPE_FLT64;
  else if ( xtype == NC_UINT64 )  datatype = DATATYPE_FLT64;
#endif

  return (datatype);
}

static
int cdfDefDatatype(int datatype, int filetype)
{
  int xtype;

  if ( datatype == DATATYPE_CPX32 || datatype == DATATYPE_CPX64 )
    Error("CDI/netCDF library does not support complex numbers!");

  if ( filetype == FILETYPE_NC4 )
    {
      if      ( datatype == DATATYPE_INT8   ) xtype = NC_BYTE;
      else if ( datatype == DATATYPE_INT16  ) xtype = NC_SHORT;
      else if ( datatype == DATATYPE_INT32  ) xtype = NC_INT;
#if  defined  (HAVE_NETCDF4)
      else if ( datatype == DATATYPE_UINT8  ) xtype = NC_UBYTE;
      else if ( datatype == DATATYPE_UINT16 ) xtype = NC_USHORT;
      else if ( datatype == DATATYPE_UINT32 ) xtype = NC_UINT;
#else
      else if ( datatype == DATATYPE_UINT8  ) xtype = NC_SHORT;
      else if ( datatype == DATATYPE_UINT16 ) xtype = NC_INT;
      else if ( datatype == DATATYPE_UINT32 ) xtype = NC_INT;
#endif
      else if ( datatype == DATATYPE_FLT64  ) xtype = NC_DOUBLE;
      else                                    xtype = NC_FLOAT;
    }
  else
    {
      if      ( datatype == DATATYPE_INT8   ) xtype = NC_BYTE;
      else if ( datatype == DATATYPE_INT16  ) xtype = NC_SHORT;
      else if ( datatype == DATATYPE_INT32  ) xtype = NC_INT;
      else if ( datatype == DATATYPE_UINT8  ) xtype = NC_SHORT;
      else if ( datatype == DATATYPE_UINT16 ) xtype = NC_INT;
      else if ( datatype == DATATYPE_UINT32 ) xtype = NC_INT;
      else if ( datatype == DATATYPE_FLT64  ) xtype = NC_DOUBLE;
      else                                    xtype = NC_FLOAT;
    }

  return (xtype);
}

static inline void *
resizeBuf(void **buf, size_t *bufSize, size_t reqSize)
{
  if (reqSize > *bufSize)
    {
      *buf = xrealloc(*buf, reqSize);
      *bufSize = reqSize;
    }
  return *buf;
}

static
void defineAttributes(int vlistID, int varID, int fileID, int ncvarID)
{
  int natts, iatt;
  int atttype, attlen;
  size_t len;
  char attname[CDI_MAX_NAME+1];
  void *attBuf = NULL;
  size_t attBufSize = 0;

  vlistInqNatts(vlistID, varID, &natts);

  for ( iatt = 0; iatt < natts; iatt++ )
    {
      vlistInqAtt(vlistID, varID, iatt, attname, &atttype, &attlen);

      if ( attlen == 0 ) continue;

      if ( atttype == DATATYPE_TXT )
        {
          size_t attSize = (size_t)attlen*sizeof(char);
          char *atttxt = (char *)resizeBuf(&attBuf, &attBufSize, attSize);
          vlistInqAttTxt(vlistID, varID, attname, attlen, atttxt);
          len = (size_t)attlen;
          cdf_put_att_text(fileID, ncvarID, attname, len, atttxt);
        }
      else if ( atttype == DATATYPE_INT16 || atttype == DATATYPE_INT32 )
        {
          size_t attSize = (size_t)attlen*sizeof(int);
          int *attint = (int *)resizeBuf(&attBuf, &attBufSize, attSize);
          vlistInqAttInt(vlistID, varID, attname, attlen, &attint[0]);
          len = (size_t)attlen;
          cdf_put_att_int(fileID, ncvarID, attname, atttype == DATATYPE_INT16 ? NC_SHORT : NC_INT, len, attint);
        }
      else if ( atttype == DATATYPE_FLT32 || atttype == DATATYPE_FLT64 )
        {
          size_t attSize = (size_t)attlen * sizeof(double);
          double *attflt = (double *)resizeBuf(&attBuf, &attBufSize, attSize);
          vlistInqAttFlt(vlistID, varID, attname, attlen, attflt);
          len = (size_t)attlen;
          cdf_put_att_double(fileID, ncvarID, attname, atttype == DATATYPE_FLT32 ? NC_FLOAT : NC_DOUBLE, len, attflt);
        }
    }
  free(attBuf);
}


void cdfCopyRecord(stream_t *streamptr2, stream_t *streamptr1)
{
  int memtype  = MEMTYPE_DOUBLE;
  int vlistID1 = streamptr1->vlistID;
  int tsID     = streamptr1->curTsID;
  int vrecID   = streamptr1->tsteps[tsID].curRecID;
  int recID    = streamptr1->tsteps[tsID].recIDs[vrecID];
  int ivarID   = streamptr1->tsteps[tsID].records[recID].varID;
  int gridID   = vlistInqVarGrid(vlistID1, ivarID);
  int datasize = gridInqSize(gridID);

  double *data = xmalloc((size_t)datasize * sizeof (double));

  int nmiss;
  cdfReadRecord(streamptr1, data, &nmiss);
  cdf_write_record(streamptr2, memtype, data, nmiss);

  free(data);
}

/* not used
int cdfInqRecord(stream_t *streamptr, int *varID, int *levelID)
{
  int tsID, recID;

  recID = streamptr->tsteps[0].curRecID++;
  printf("cdfInqRecord recID %d %d\n", recID, streamptr->tsteps[0].curRecID);
  printf("cdfInqRecord tsID %d\n", streamptr->curTsID);

  if ( streamptr->tsteps[0].curRecID >= streamptr->tsteps[0].nrecs )
    {
      streamptr->tsteps[0].curRecID = 0;
    }

  *varID   = streamptr->tsteps[0].records[recID].varID;
  *levelID = streamptr->tsteps[0].records[recID].levelID;

  streamptr->record->varID   = *varID;
  streamptr->record->levelID = *levelID;

  if ( CDI_Debug )
    Message("recID = %d  varID = %d  levelID = %d", recID, *varID, *levelID);

  return (recID+1);
}
*/


void cdfDefRecord(stream_t *streamptr)
{
  (void)streamptr;
}


static
void cdfWriteGridTraj(stream_t *streamptr, int gridID)
{
  int vlistID = streamptr->vlistID;
  int fileID  = streamptr->fileID;

  int gridindex = vlistGridIndex(vlistID, gridID);
  int lonID = streamptr->xdimID[gridindex],
    latID = streamptr->ydimID[gridindex];

  double xlon = gridInqXval(gridID, 0),
    xlat = gridInqYval(gridID, 0);
  int tsID = streamptr->curTsID;
  size_t index = (size_t)tsID;

  cdf_put_var1_double(fileID, lonID, &index, &xlon);
  cdf_put_var1_double(fileID, latID, &index, &xlat);
}

static
void cdfReadGridTraj(stream_t *streamptr, int gridID)
{
  int vlistID = streamptr->vlistID;
  int fileID  = streamptr->fileID;

  int gridindex = vlistGridIndex(vlistID, gridID);
  int lonID = streamptr->xdimID[gridindex];
  int latID = streamptr->ydimID[gridindex];

  int tsID = streamptr->curTsID;
  size_t index = (size_t)tsID;

  double xlon, xlat;
  cdf_get_var1_double(fileID, lonID, &index, &xlon);
  cdf_get_var1_double(fileID, latID, &index, &xlat);

  gridDefXvals(gridID, &xlon);
  gridDefYvals(gridID, &xlat);
}


static
void cdfDefVarDeflate(int ncid, int ncvarid, int deflate_level)
{
#if  defined  (HAVE_NETCDF4)
  int retval;
  /* Set chunking, shuffle, and deflate. */
  int shuffle = 1;
  int deflate = 1;

  if ( deflate_level < 1 || deflate_level > 9 ) deflate_level = 1;

  if ((retval = nc_def_var_deflate(ncid, ncvarid, shuffle, deflate, deflate_level)))
    {
      Error("nc_def_var_deflate failed, status = %d", retval);
    }
#else
  static int lwarn = TRUE;

  if ( lwarn )
    {
      lwarn = FALSE;
      Warning("Deflate compression failed, netCDF4 not available!");
    }
#endif
}


#if defined(NC_SZIP_NN_OPTION_MASK)
static
void cdfDefVarSzip(int ncid, int ncvarid)
{
  int retval;
  /* Set options_mask and bits_per_pixel. */
  int options_mask = NC_SZIP_NN_OPTION_MASK;
  int bits_per_pixel = 16;

  if ((retval = nc_def_var_szip(ncid, ncvarid, options_mask, bits_per_pixel)))
    {
      if ( retval == NC_EINVAL )
        {
          static int lwarn = TRUE;

          if ( lwarn )
            {
              lwarn = FALSE;
              Warning("netCDF4/Szip compression not compiled in!");
            }
        }
      else
        Error("nc_def_var_szip failed, status = %d", retval);
    }
}
#endif

static
void cdfDefVarMissval(stream_t *streamptr, int varID, int dtype, int lcheck)
{
  if ( streamptr->vars[varID].defmiss == FALSE )
    {
      int fileID;
      int ncvarid;
      double missval;
      int vlistID;
      int xtype;

      vlistID = streamptr->vlistID;
      fileID  = streamptr->fileID;
      ncvarid = streamptr->vars[varID].ncvarid;
      missval = vlistInqVarMissval(vlistID, varID);

      if ( lcheck && streamptr->ncmode == 2 ) cdf_redef(fileID);

      xtype = cdfDefDatatype(dtype, streamptr->filetype);

      if ( xtype == NC_BYTE && missval > 127 && missval < 256 ) xtype = NC_INT;

      cdf_put_att_double(fileID, ncvarid, "_FillValue", (nc_type) xtype, 1, &missval);
      cdf_put_att_double(fileID, ncvarid, "missing_value", (nc_type) xtype, 1, &missval);

      if ( lcheck && streamptr->ncmode == 2 ) cdf_enddef(fileID);

      streamptr->vars[varID].defmiss = TRUE;
    }
}


void cdf_write_record(stream_t *streamptr, int memtype, const void *data, int nmiss)
{
  int varID;
  int levelID;

  varID   = streamptr->record->varID;
  levelID = streamptr->record->levelID;

  if ( CDI_Debug ) Message("streamID = %d  varID = %d", streamptr->self, varID);

  cdf_write_var_slice(streamptr, varID, levelID, memtype, data, nmiss);
}

void cdfReadRecord(stream_t *streamptr, double *data, int *nmiss)
{
  if ( CDI_Debug ) Message("streamID = %d", streamptr->self);

  int tsID    = streamptr->curTsID;
  int vrecID  = streamptr->tsteps[tsID].curRecID;
  int recID   = streamptr->tsteps[tsID].recIDs[vrecID];
  int varID   = streamptr->tsteps[tsID].records[recID].varID;
  int levelID = streamptr->tsteps[tsID].records[recID].levelID;

  cdfReadVarSliceDP(streamptr, varID, levelID, data, nmiss);
}

static
void cdfDefTimeValue(stream_t *streamptr, int tsID)
{
  int fileID = streamptr->fileID;

  if ( CDI_Debug )
    Message("streamID = %d, fileID = %d", streamptr->self, fileID);

  taxis_t *taxis = &streamptr->tsteps[tsID].taxis;

  if ( streamptr->ncmode == 1 )
    {
      cdf_enddef(fileID);
      streamptr->ncmode = 2;
    }

  size_t index = (size_t)tsID;

  double timevalue = cdiEncodeTimeval(taxis->vdate, taxis->vtime, &streamptr->tsteps[0].taxis);
  if ( CDI_Debug ) Message("tsID = %d  timevalue = %f", tsID, timevalue);

  int ncvarid = streamptr->basetime.ncvarid;
  cdf_put_var1_double(fileID, ncvarid, &index, &timevalue);

  if ( taxis->has_bounds )
    {
      size_t start[2], count[2];

      ncvarid = streamptr->basetime.ncvarboundsid;

      timevalue = cdiEncodeTimeval(taxis->vdate_lb, taxis->vtime_lb, &streamptr->tsteps[0].taxis);
      start[0] = (size_t)tsID; count[0] = 1; start[1] = 0; count[1] = 1;
      cdf_put_vara_double(fileID, ncvarid, start, count, &timevalue);

      timevalue = cdiEncodeTimeval(taxis->vdate_ub, taxis->vtime_ub, &streamptr->tsteps[0].taxis);
      start[0] = (size_t)tsID; count[0] = 1; start[1] = 1; count[1] = 1;
      cdf_put_vara_double(fileID, ncvarid, start, count, &timevalue);
    }

  ncvarid = streamptr->basetime.leadtimeid;
  if ( taxis->type == TAXIS_FORECAST && ncvarid != UNDEFID )
    {
      timevalue = taxis->fc_period;
      cdf_put_var1_double(fileID, ncvarid, &index, &timevalue);
    }

  /*
printf("fileID = %d %d %d %f\n", fileID, time_varid, index, timevalue);
  */
}

static
int cdfDefTimeBounds(int fileID, int nctimevarid, int nctimedimid, char* taxis_name, taxis_t* taxis)
{
  int time_bndsid = -1;
  int dims[2];
  char tmpstr[CDI_MAX_NAME];

  dims[0] = nctimedimid;

  /* fprintf(stderr, "time has bounds\n"); */

  if ( nc_inq_dimid(fileID, "bnds", &dims[1]) != NC_NOERR )
    cdf_def_dim(fileID, "bnds", 2, &dims[1]);

  if ( taxis->climatology )
    {
      strcpy(tmpstr, "climatology");
      strcat(tmpstr, "_bnds");
      cdf_def_var(fileID, tmpstr, NC_DOUBLE, 2, dims, &time_bndsid);

      cdf_put_att_text(fileID, nctimevarid, "climatology", strlen(tmpstr), tmpstr);
    }
  else
    {
      strcpy(tmpstr, taxis_name);
      strcat(tmpstr, "_bnds");
      cdf_def_var(fileID, tmpstr, NC_DOUBLE, 2, dims, &time_bndsid);

      cdf_put_att_text(fileID, nctimevarid, "bounds", strlen(tmpstr), tmpstr);
    }

  return (time_bndsid);
}

static
void cdfDefTimeUnits(char *unitstr, taxis_t* taxis0, taxis_t* taxis)
{
  unitstr[0] = 0;

  if ( taxis0->type == TAXIS_ABSOLUTE )
    {
      if ( taxis0->unit == TUNIT_YEAR )
        sprintf(unitstr, "year as %s", "%Y.%f");
      else if ( taxis0->unit == TUNIT_MONTH )
        sprintf(unitstr, "month as %s", "%Y%m.%f");
      else
        sprintf(unitstr, "day as %s", "%Y%m%d.%f");
    }
  else
    {
      int year, month, day, hour, minute, second;
      int rdate, rtime;
      int timeunit;

      timeunit = taxis->unit;
      if ( timeunit == -1 ) timeunit = TUNIT_HOUR;
      rdate    = taxis->rdate;
      rtime    = taxis->rtime;
      if ( rdate == -1 )
        {
          rdate  = taxis->vdate;
          rtime  = taxis->vtime;
        }

      cdiDecodeDate(rdate, &year, &month, &day);
      cdiDecodeTime(rtime, &hour, &minute, &second);

      if ( timeunit == TUNIT_QUARTER   ) timeunit = TUNIT_MINUTE;
      if ( timeunit == TUNIT_30MINUTES ) timeunit = TUNIT_MINUTE;
      if ( timeunit == TUNIT_3HOURS  ||
	   timeunit == TUNIT_6HOURS  ||
	   timeunit == TUNIT_12HOURS ) timeunit = TUNIT_HOUR;

      sprintf(unitstr, "%s since %d-%02d-%02d %02d:%02d:%02d",
              tunitNamePtr(timeunit), year, month, day, hour, minute, second);
    }
}

static
void cdfDefForecastTimeUnits(char *unitstr, int timeunit)
{
  unitstr[0] = 0;

  if ( timeunit == -1 ) timeunit = TUNIT_HOUR;

  if ( timeunit == TUNIT_QUARTER   ) timeunit = TUNIT_MINUTE;
  if ( timeunit == TUNIT_30MINUTES ) timeunit = TUNIT_MINUTE;
  if ( timeunit == TUNIT_3HOURS  ||
       timeunit == TUNIT_6HOURS  ||
       timeunit == TUNIT_12HOURS ) timeunit = TUNIT_HOUR;

  sprintf(unitstr, "%s", tunitNamePtr(timeunit));
}

static
void cdfDefCalendar(int fileID, int ncvarid, int calendar)
{
  size_t len;
  char calstr[80];

  calstr[0] = 0;

  if      ( calendar == CALENDAR_STANDARD )  strcpy(calstr, "standard");
  else if ( calendar == CALENDAR_PROLEPTIC ) strcpy(calstr, "proleptic_gregorian");
  else if ( calendar == CALENDAR_NONE )      strcpy(calstr, "none");
  else if ( calendar == CALENDAR_360DAYS )   strcpy(calstr, "360_day");
  else if ( calendar == CALENDAR_365DAYS )   strcpy(calstr, "365_day");
  else if ( calendar == CALENDAR_366DAYS )   strcpy(calstr, "366_day");

  len = strlen(calstr);

  if ( len ) cdf_put_att_text(fileID, ncvarid, "calendar", len, calstr);
}

static
void cdfDefTime(stream_t* streamptr)
{
  int fileID;
  int time_varid;
  int time_dimid;
  int time_bndsid = -1;
  char unitstr[CDI_MAX_NAME];
  char tmpstr[CDI_MAX_NAME];
  char default_name[] = "time";
  char* taxis_name = default_name;
  size_t len;
  taxis_t* taxis;

  if ( streamptr->basetime.ncvarid != UNDEFID ) return;

  fileID = streamptr->fileID;

  if ( streamptr->ncmode == 0 ) streamptr->ncmode = 1;

  if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

  taxis = &streamptr->tsteps[0].taxis;

  if ( taxis->name && taxis->name[0] ) taxis_name = taxis->name;

  cdf_def_dim(fileID, taxis_name, NC_UNLIMITED, &time_dimid);
  streamptr->basetime.ncdimid = time_dimid;

  cdf_def_var(fileID, taxis_name, NC_DOUBLE, 1, &time_dimid, &time_varid);

  streamptr->basetime.ncvarid = time_varid;

  strcpy(tmpstr, "time");
  cdf_put_att_text(fileID, time_varid, "standard_name", strlen(tmpstr), tmpstr);

  if ( taxis->longname && taxis->longname[0] )
    cdf_put_att_text(fileID, time_varid, "long_name", strlen(taxis->longname), taxis->longname);

  if ( taxis->has_bounds )
    {
      time_bndsid = cdfDefTimeBounds(fileID, time_varid, time_dimid, taxis_name, taxis);
      streamptr->basetime.ncvarboundsid = time_bndsid;
    }

  cdfDefTimeUnits(unitstr, &streamptr->tsteps[0].taxis, taxis);

  len = strlen(unitstr);
  if ( len )
    {
      cdf_put_att_text(fileID, time_varid, "units", len, unitstr);
      /*
      if ( taxis->has_bounds )
        cdf_put_att_text(fileID, time_bndsid, "units", len, unitstr);
      */
    }

  if ( taxis->calendar != -1 )
    {
      cdfDefCalendar(fileID, time_varid, taxis->calendar);
      /*
      if ( taxis->has_bounds )
        cdfDefCalendar(fileID, time_bndsid, taxis->calendar);
      */
    }

  if ( taxis->type == TAXIS_FORECAST )
    {
      int leadtimeid;

      cdf_def_var(fileID, "leadtime", NC_DOUBLE, 1, &time_dimid, &leadtimeid);

      streamptr->basetime.leadtimeid = leadtimeid;

      strcpy(tmpstr, "forecast_period");
      cdf_put_att_text(fileID, leadtimeid, "standard_name", strlen(tmpstr), tmpstr);

      strcpy(tmpstr, "Time elapsed since the start of the forecast");
      cdf_put_att_text(fileID, leadtimeid, "long_name", strlen(tmpstr), tmpstr);

      cdfDefForecastTimeUnits(unitstr, taxis->fc_unit);

      len = strlen(unitstr);
      if ( len ) cdf_put_att_text(fileID, leadtimeid, "units", len, unitstr);
    }

  if ( streamptr->ncmode == 2 ) cdf_enddef(fileID);
}


void cdfDefTimestep(stream_t *streamptr, int tsID)
{
  int vlistID = streamptr->vlistID;

  if ( vlistHasTime(vlistID) ) cdfDefTime(streamptr);

  cdfDefTimeValue(streamptr, tsID);
}

static
void cdfDefComplex(stream_t *streamptr, int gridID)
{
  char axisname[] = "nc2";
  int dimID = UNDEFID;
  int gridID0, gridtype0;
  int vlistID = streamptr->vlistID;
  int fileID  = streamptr->fileID;

  int ngrids = vlistNgrids(vlistID);

  for ( int index = 0; index < ngrids; index++ )
    {
      if ( streamptr->xdimID[index] != UNDEFID )
        {
          gridID0 = vlistGrid(vlistID, index);
          gridtype0 = gridInqType(gridID0);
          if ( gridtype0 == GRID_SPECTRAL || gridtype0 == GRID_FOURIER )
            {
              dimID = streamptr->xdimID[index];
              break;
            }
        }
    }

  if ( dimID == UNDEFID )
    {
      size_t dimlen = 2;

      if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

      cdf_def_dim(fileID, axisname, dimlen, &dimID);

      cdf_enddef(fileID);
      streamptr->ncmode = 2;
    }

  int gridindex = vlistGridIndex(vlistID, gridID);
  streamptr->xdimID[gridindex] = dimID;
}


static
void cdfDefSP(stream_t *streamptr, int gridID)
{
  /*
  char longname[] = "Spherical harmonic coefficient";
  */
  char axisname[5] = "nspX";
  int index, iz = 0;
  int gridID0, gridtype0, gridindex;
  int dimID = UNDEFID;
  int ngrids;
  int fileID;
  int vlistID;

  vlistID = streamptr->vlistID;
  fileID  = streamptr->fileID;

  ngrids = vlistNgrids(vlistID);

  size_t dimlen = (size_t)gridInqSize(gridID)/2;

  for ( index = 0; index < ngrids; index++ )
    {
      if ( streamptr->ydimID[index] != UNDEFID )
        {
          gridID0 = vlistGrid(vlistID, index);
          gridtype0 = gridInqType(gridID0);
          if ( gridtype0 == GRID_SPECTRAL )
            {
              size_t dimlen0 = (size_t)gridInqSize(gridID0)/2;
              if ( dimlen == dimlen0 )
                {
                  dimID = streamptr->ydimID[index];
                  break;
                }
              else
                iz++;
            }
        }
    }

  if ( dimID == UNDEFID )
    {
      if ( iz == 0 ) axisname[3] = '\0';
      else           sprintf(&axisname[3], "%1d", iz+1);

      if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

      cdf_def_dim(fileID, axisname, dimlen, &dimID);

      cdf_enddef(fileID);
      streamptr->ncmode = 2;
    }

  gridindex = vlistGridIndex(vlistID, gridID);
  streamptr->ydimID[gridindex] = dimID;
}


static
void cdfDefFC(stream_t *streamptr, int gridID)
{
  char axisname[5] = "nfcX";
  int index, iz = 0;
  int gridID0, gridtype0, gridindex;
  int dimID = UNDEFID;
  int ngrids;
  int fileID;
  int vlistID;

  vlistID = streamptr->vlistID;
  fileID  = streamptr->fileID;

  ngrids = vlistNgrids(vlistID);

  size_t dimlen = (size_t)gridInqSize(gridID)/2;

  for ( index = 0; index < ngrids; index++ )
    {
      if ( streamptr->ydimID[index] != UNDEFID )
        {
          gridID0 = vlistGrid(vlistID, index);
          gridtype0 = gridInqType(gridID0);
          if ( gridtype0 == GRID_FOURIER )
            {
              size_t dimlen0 = (size_t)gridInqSize(gridID0)/2;
              if ( dimlen == dimlen0 )
                {
                  dimID = streamptr->ydimID[index];
                  break;
                }
              else
                iz++;
            }
        }
    }

  if ( dimID == UNDEFID )
    {
      if ( iz == 0 ) axisname[3] = '\0';
      else           sprintf(&axisname[3], "%1d", iz+1);

      if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

      cdf_def_dim(fileID, axisname, dimlen, &dimID);

      cdf_enddef(fileID);
      streamptr->ncmode = 2;
    }

  gridindex = vlistGridIndex(vlistID, gridID);
  streamptr->ydimID[gridindex] = dimID;
}


static
void cdfDefTrajLon(stream_t *streamptr, int gridID)
{
  char units[CDI_MAX_NAME];
  char longname[CDI_MAX_NAME];
  char stdname[CDI_MAX_NAME];
  char axisname[CDI_MAX_NAME];
  int gridtype, gridindex;
  int dimID = UNDEFID;
  int fileID;
  int dimlen;
  size_t len;
  int ncvarid;
  int vlistID;
  int xtype = NC_DOUBLE;

  if ( gridInqPrec(gridID) == DATATYPE_FLT32 ) xtype = NC_FLOAT;

  vlistID = streamptr->vlistID;
  fileID  = streamptr->fileID;

  gridtype = gridInqType(gridID);
  dimlen = gridInqXsize(gridID);
  if ( dimlen != 1 ) Error("Xsize isn't 1 for %s grid!", gridNamePtr(gridtype));

  gridindex = vlistGridIndex(vlistID, gridID);
  ncvarid = streamptr->xdimID[gridindex];

  gridInqXname(gridID, axisname);
  gridInqXlongname(gridID, longname);
  gridInqXstdname(gridID, stdname);
  gridInqXunits(gridID, units);

  if ( ncvarid == UNDEFID )
    {
      dimID = streamptr->basetime.ncvarid;

      if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

      cdf_def_var(fileID, axisname, (nc_type) xtype, 1, &dimID, &ncvarid);

      if ( (len = strlen(stdname)) )
        cdf_put_att_text(fileID, ncvarid, "standard_name", len, stdname);
      if ( (len = strlen(longname)) )
        cdf_put_att_text(fileID, ncvarid, "long_name", len, longname);
      if ( (len = strlen(units)) )
        cdf_put_att_text(fileID, ncvarid, "units", len, units);

      cdf_enddef(fileID);
      streamptr->ncmode = 2;
    }

  streamptr->xdimID[gridindex] = ncvarid; /* var ID for trajectory !!! */
}


static
void cdfDefTrajLat(stream_t *streamptr, int gridID)
{
  char units[] = "degrees_north";
  char longname[] = "latitude";
  char stdname[] = "latitude";
  char axisname[] = "tlat";
  int gridtype, gridindex;
  int dimID = UNDEFID;
  int fileID;
  int dimlen;
  size_t len;
  int ncvarid;
  int vlistID;
  int xtype = NC_DOUBLE;

  if ( gridInqPrec(gridID) == DATATYPE_FLT32 ) xtype = NC_FLOAT;

  vlistID = streamptr->vlistID;
  fileID = streamptr->fileID;

  gridtype = gridInqType(gridID);
  dimlen = gridInqYsize(gridID);
  if ( dimlen != 1 ) Error("Ysize isn't 1 for %s grid!", gridNamePtr(gridtype));

  gridindex = vlistGridIndex(vlistID, gridID);
  ncvarid = streamptr->ydimID[gridindex];

  gridInqYname(gridID, axisname);
  gridInqYlongname(gridID, longname);
  gridInqYstdname(gridID, stdname);
  gridInqYunits(gridID, units);

  if ( ncvarid == UNDEFID )
    {
      dimID = streamptr->basetime.ncvarid;

      if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

      cdf_def_var(fileID, axisname, (nc_type) xtype, 1, &dimID, &ncvarid);

      if ( (len = strlen(stdname)) )
        cdf_put_att_text(fileID, ncvarid, "standard_name", len, stdname);
      if ( (len = strlen(longname)) )
        cdf_put_att_text(fileID, ncvarid, "long_name", len, longname);
      if ( (len = strlen(units)) )
        cdf_put_att_text(fileID, ncvarid, "units", len, units);

      cdf_enddef(fileID);
      streamptr->ncmode = 2;
    }

  streamptr->ydimID[gridindex] = ncvarid; /* var ID for trajectory !!! */
}


static
int checkGridName(int type, char *axisname, int fileID, int vlistID, int gridID, int ngrids, int mode)
{
  int iz, index;
  int gridID0;
  int ncdimid;
  char axisname0[CDI_MAX_NAME];
  char axisname2[CDI_MAX_NAME];
  int checkname;
  int status;

  /* check that the name is not already defined */
  checkname = TRUE;
  iz = 0;

  do
    {
      strcpy(axisname2, axisname);
      if ( iz ) sprintf(&axisname2[strlen(axisname2)], "_%d", iz+1);

      //status = nc_inq_varid(fileID, axisname2, &ncvarid);
      if ( type == 'V' ) /* type Var oder Dim */
        status = nc_inq_varid(fileID, axisname2, &ncdimid);
      else
        status = nc_inq_dimid(fileID, axisname2, &ncdimid);

      if ( status != NC_NOERR )
        {
          if ( iz )
            {
              /* check that the name does not exist for other grids */
              for ( index = 0; index < ngrids; index++ )
                {
                  gridID0 = vlistGrid(vlistID, index);
                  if ( gridID != gridID0 )
                    {
                      if ( mode == 'X' ) /* mode X or Y */
                        gridInqXname(gridID0, axisname0);
                      else
                        gridInqYname(gridID0, axisname0);

                      if ( strcmp(axisname0, axisname2) == 0 ) break;
                    }
                }
              if ( index == ngrids ) checkname = FALSE;
            }
          else
            {
              checkname = FALSE;
            }
        }

      if ( checkname ) iz++;
    }
  while (checkname && iz <= 99);


  if ( iz ) sprintf(&axisname[strlen(axisname)], "_%d", iz+1);

  return (iz);
}


static
void cdfDefXaxis(stream_t *streamptr, int gridID, int ndims)
{
  char units[CDI_MAX_NAME];
  char longname[CDI_MAX_NAME];
  char stdname[CDI_MAX_NAME];
  char axisname[CDI_MAX_NAME];
  int index;
  /*  int index2; */
  int gridID0, gridtype0, gridindex;
  int dimID = UNDEFID;
  int dimIDs[2];
  int ngrids = 0;
  int fileID;
  size_t len;
  int ncvarid = UNDEFID, ncbvarid = UNDEFID;
  int nvdimID = UNDEFID;
  int vlistID;
  int xtype = NC_DOUBLE;

  if ( gridInqPrec(gridID) == DATATYPE_FLT32 ) xtype = NC_FLOAT;

  vlistID = streamptr->vlistID;
  fileID  = streamptr->fileID;

  if ( ndims ) ngrids = vlistNgrids(vlistID);

  size_t dimlen = (size_t)gridInqXsize(gridID);
  gridindex = vlistGridIndex(vlistID, gridID);

  gridInqXname(gridID, axisname);
  gridInqXlongname(gridID, longname);
  gridInqXstdname(gridID, stdname);
  gridInqXunits(gridID, units);

  if ( axisname[0] == 0 ) Error("axis name undefined!");

  for ( index = 0; index < ngrids; index++ )
    {
      if ( streamptr->xdimID[index] != UNDEFID )
        {
          gridID0 = vlistGrid(vlistID, index);
          gridtype0 = gridInqType(gridID0);
          if ( gridtype0 == GRID_GAUSSIAN    ||
               gridtype0 == GRID_LONLAT      ||
               gridtype0 == GRID_CURVILINEAR ||
               gridtype0 == GRID_GENERIC )
            {
              size_t dimlen0 = (size_t)gridInqXsize(gridID0);
              if ( dimlen == dimlen0 )
                if ( IS_EQUAL(gridInqXval(gridID0, 0), gridInqXval(gridID, 0)) &&
                     IS_EQUAL(gridInqXval(gridID0, (int)dimlen-1), gridInqXval(gridID, (int)dimlen-1)) )
                  {
                    dimID = streamptr->xdimID[index];
                    break;
                  }
              /*
              for ( index2 = 0; index2 < index; index2++ )
                if ( streamptr->xdimID[index] == streamptr->xdimID[index2] )
                  break;
              if ( index2 == index ) iz++;
              */
            }
        }
    }

  if ( dimID == UNDEFID )
    {
      int status;
      status = checkGridName('V', axisname, fileID, vlistID, gridID, ngrids, 'X');
      if ( status == 0 && ndims )
        status = checkGridName('D', axisname, fileID, vlistID, gridID, ngrids, 'X');

      if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

      if ( ndims )
        {
          cdf_def_dim(fileID, axisname, dimlen, &dimID);

          if ( gridInqXboundsPtr(gridID) || gridInqYboundsPtr(gridID) )
            {
              size_t nvertex = 2;
              if ( nc_inq_dimid(fileID, "bnds", &nvdimID) != NC_NOERR )
                cdf_def_dim(fileID, "bnds", nvertex, &nvdimID);
            }
        }

      if ( gridInqXvalsPtr(gridID) )
        {
          cdf_def_var(fileID, axisname, (nc_type) xtype, ndims, &dimID, &ncvarid);

          if ( (len = strlen(stdname)) )
            cdf_put_att_text(fileID, ncvarid, "standard_name", len, stdname);
          if ( (len = strlen(longname)) )
            cdf_put_att_text(fileID, ncvarid, "long_name", len, longname);
          if ( (len = strlen(units)) )
            cdf_put_att_text(fileID, ncvarid, "units", len, units);

          cdf_put_att_text(fileID, ncvarid, "axis", 1, "X");

          if ( gridInqXboundsPtr(gridID) && nvdimID != UNDEFID )
            {
              strcat(axisname, "_bnds");
              dimIDs[0] = dimID;
              dimIDs[1] = nvdimID;
              cdf_def_var(fileID, axisname, (nc_type) xtype, 2, dimIDs, &ncbvarid);
              cdf_put_att_text(fileID, ncvarid, "bounds", strlen(axisname), axisname);
            }
          /*
          if ( gridIsRotated(gridID) )
            {
              double north_pole = gridInqXpole(gridID);
              cdf_put_att_double(fileID, ncvarid, "north_pole", NC_DOUBLE, 1, &north_pole);
            }
          */
        }

      cdf_enddef(fileID);
      streamptr->ncmode = 2;

      if ( ncvarid  != UNDEFID ) cdf_put_var_double(fileID, ncvarid, gridInqXvalsPtr(gridID));
      if ( ncbvarid != UNDEFID ) cdf_put_var_double(fileID, ncbvarid, gridInqXboundsPtr(gridID));

      if ( ndims == 0 ) streamptr->ncxvarID[gridindex] = ncvarid;
    }

  streamptr->xdimID[gridindex] = dimID;
}


static
void cdfDefYaxis(stream_t *streamptr, int gridID, int ndims)
{
  char units[CDI_MAX_NAME];
  char longname[CDI_MAX_NAME];
  char stdname[CDI_MAX_NAME];
  char axisname[CDI_MAX_NAME];
  int index;
  /*  int index2; */
  int gridID0, gridtype0, gridindex;
  int dimID = UNDEFID;
  int dimIDs[2];
  int ngrids = 0;
  int fileID;
  size_t len;
  int ncvarid = UNDEFID, ncbvarid = UNDEFID;
  int nvdimID = UNDEFID;
  int vlistID;
  int xtype = NC_DOUBLE;

  if ( gridInqPrec(gridID) == DATATYPE_FLT32 ) xtype = NC_FLOAT;

  vlistID = streamptr->vlistID;
  fileID  = streamptr->fileID;

  if ( ndims ) ngrids = vlistNgrids(vlistID);

  size_t dimlen = (size_t)gridInqYsize(gridID);
  gridindex = vlistGridIndex(vlistID, gridID);

  gridInqYname(gridID, axisname);
  gridInqYlongname(gridID, longname);
  gridInqYstdname(gridID, stdname);
  gridInqYunits(gridID, units);

  if ( axisname[0] == 0 ) Error("axis name undefined!");

  for ( index = 0; index < ngrids; index++ )
    {
      if ( streamptr->ydimID[index] != UNDEFID )
        {
          gridID0 = vlistGrid(vlistID, index);
          gridtype0 = gridInqType(gridID0);
          if ( gridtype0 == GRID_GAUSSIAN    ||
               gridtype0 == GRID_LONLAT      ||
               gridtype0 == GRID_CURVILINEAR ||
               gridtype0 == GRID_GENERIC )
            {
              size_t dimlen0 = (size_t)gridInqYsize(gridID0);
              if ( dimlen == dimlen0 )
                if ( IS_EQUAL(gridInqYval(gridID0, 0), gridInqYval(gridID, 0)) &&
                     IS_EQUAL(gridInqYval(gridID0, (int)dimlen-1), gridInqYval(gridID, (int)dimlen-1)) )
                  {
                    dimID = streamptr->ydimID[index];
                    break;
                  }
              /*
              for ( index2 = 0; index2 < index; index2++ )
                if ( streamptr->ydimID[index] == streamptr->ydimID[index2] )
                  break;
              if ( index2 == index ) iz++;
              */
            }
        }
    }

  if ( dimID == UNDEFID )
    {
      int status;
      status = checkGridName('V', axisname, fileID, vlistID, gridID, ngrids, 'Y');
      if ( status == 0 && ndims )
        status = checkGridName('D', axisname, fileID, vlistID, gridID, ngrids, 'Y');

      if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

      if ( ndims )
        {
          cdf_def_dim(fileID, axisname, dimlen, &dimID);

          if ( gridInqXboundsPtr(gridID) || gridInqYboundsPtr(gridID) )
            {
              size_t nvertex = 2;
              if ( nc_inq_dimid(fileID, "bnds", &nvdimID) != NC_NOERR )
                cdf_def_dim(fileID, "bnds", nvertex, &nvdimID);
            }
        }

      if ( gridInqYvalsPtr(gridID) )
        {
          cdf_def_var(fileID, axisname, (nc_type) xtype, ndims, &dimID, &ncvarid);

          if ( (len = strlen(stdname)) )
            cdf_put_att_text(fileID, ncvarid, "standard_name", len, stdname);
          if ( (len = strlen(longname)) )
            cdf_put_att_text(fileID, ncvarid, "long_name", len, longname);
          if ( (len = strlen(units)) )
            cdf_put_att_text(fileID, ncvarid, "units", len, units);

          cdf_put_att_text(fileID, ncvarid, "axis", 1, "Y");

          if ( gridInqYboundsPtr(gridID) && nvdimID != UNDEFID )
            {
              strcat(axisname, "_bnds");
              dimIDs[0] = dimID;
              dimIDs[1] = nvdimID;
              cdf_def_var(fileID, axisname, (nc_type) xtype, 2, dimIDs, &ncbvarid);
              cdf_put_att_text(fileID, ncvarid, "bounds", strlen(axisname), axisname);
            }
          /*
          if ( gridIsRotated(gridID) )
            {
              double north_pole = gridInqYpole(gridID);
              cdf_put_att_double(fileID, ncvarid, "north_pole", NC_DOUBLE, 1, &north_pole);
            }
          */
        }

      cdf_enddef(fileID);
      streamptr->ncmode = 2;

      if ( ncvarid  != UNDEFID ) cdf_put_var_double(fileID, ncvarid, gridInqYvalsPtr(gridID));
      if ( ncbvarid != UNDEFID ) cdf_put_var_double(fileID, ncbvarid, gridInqYboundsPtr(gridID));

      if ( ndims == 0 ) streamptr->ncyvarID[gridindex] = ncvarid;
    }

  streamptr->ydimID[gridindex] = dimID;
}


static
void cdfGridCompress(int fileID, int ncvarid, int gridsize, int filetype, int comptype)
{
#if  defined  (HAVE_NETCDF4)
  if ( gridsize > 1 && comptype == COMPRESS_ZIP && (filetype == FILETYPE_NC4 || filetype == FILETYPE_NC4C) )
    {
      nc_def_var_chunking(fileID, ncvarid, NC_CHUNKED, NULL);
      cdfDefVarDeflate(fileID, ncvarid, 1);
    }
#endif
}


static
void cdfDefCurvilinear(stream_t *streamptr, int gridID)
{
  char xunits[CDI_MAX_NAME];
  char xlongname[CDI_MAX_NAME];
  char xstdname[CDI_MAX_NAME];
  char yunits[CDI_MAX_NAME];
  char ylongname[CDI_MAX_NAME];
  char ystdname[CDI_MAX_NAME];
  char xaxisname[CDI_MAX_NAME];
  char yaxisname[CDI_MAX_NAME];
  char xdimname[4] = "x";
  char ydimname[4] = "y";
  int index;
  int gridID0, gridtype0, gridindex;
  int xdimID = UNDEFID;
  int ydimID = UNDEFID;
  int dimIDs[3];
  int ngrids;
  int fileID;
  size_t len;
  int ncxvarid = UNDEFID, ncyvarid = UNDEFID;
  int ncbxvarid = UNDEFID, ncbyvarid = UNDEFID, ncavarid = UNDEFID;
  int nvdimID = UNDEFID;
  int vlistID;
  int xtype = NC_DOUBLE;

  if ( gridInqPrec(gridID) == DATATYPE_FLT32 ) xtype = NC_FLOAT;

  vlistID = streamptr->vlistID;
  fileID  = streamptr->fileID;

  ngrids = vlistNgrids(vlistID);

  size_t xdimlen = (size_t)gridInqXsize(gridID);
  size_t ydimlen = (size_t)gridInqYsize(gridID);
  gridindex = vlistGridIndex(vlistID, gridID);

  gridInqXname(gridID, xaxisname);
  gridInqXlongname(gridID, xlongname);
  gridInqXstdname(gridID, xstdname);
  gridInqXunits(gridID, xunits);
  gridInqYname(gridID, yaxisname);
  gridInqYlongname(gridID, ylongname);
  gridInqYstdname(gridID, ystdname);
  gridInqYunits(gridID, yunits);

  for ( index = 0; index < ngrids; index++ )
    {
      if ( streamptr->xdimID[index] != UNDEFID )
        {
          gridID0 = vlistGrid(vlistID, index);
          gridtype0 = gridInqType(gridID0);
          if ( gridtype0 == GRID_GAUSSIAN    ||
               gridtype0 == GRID_LONLAT      ||
               gridtype0 == GRID_CURVILINEAR ||
               gridtype0 == GRID_GENERIC )
            {
              size_t dimlen0 = (size_t)gridInqXsize(gridID0);
              if ( xdimlen == dimlen0 )
                if ( IS_EQUAL(gridInqXval(gridID0, 0), gridInqXval(gridID, 0)) &&
                     IS_EQUAL(gridInqXval(gridID0, (int)xdimlen-1), gridInqXval(gridID, (int)xdimlen-1)) )
                  {
                    xdimID = streamptr->xdimID[index];
                    ncxvarid = streamptr->ncxvarID[index];
                    break;
                  }
              dimlen0 = (size_t)gridInqYsize(gridID0);
              if ( ydimlen == dimlen0 )
                if ( IS_EQUAL(gridInqYval(gridID0, 0), gridInqYval(gridID, 0)) &&
                     IS_EQUAL(gridInqYval(gridID0, (int)xdimlen-1), gridInqYval(gridID, (int)xdimlen-1)) )
                  {
                    ydimID = streamptr->ydimID[index];
                    ncyvarid = streamptr->ncyvarID[index];
                    break;
                  }
            }
        }
    }

  if ( xdimID == UNDEFID || ydimID == UNDEFID )
    {
      checkGridName('V', xaxisname, fileID, vlistID, gridID, ngrids, 'X');
      checkGridName('V', yaxisname, fileID, vlistID, gridID, ngrids, 'Y');
      checkGridName('D', xdimname, fileID, vlistID, gridID, ngrids, 'X');
      checkGridName('D', ydimname, fileID, vlistID, gridID, ngrids, 'Y');

      if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

      cdf_def_dim(fileID, xdimname, xdimlen, &xdimID);
      cdf_def_dim(fileID, ydimname, ydimlen, &ydimID);

      if ( gridInqXboundsPtr(gridID) || gridInqYboundsPtr(gridID) )
        {
          size_t nvertex = 4;
          if ( nc_inq_dimid(fileID, "vertices", &nvdimID) != NC_NOERR )
            cdf_def_dim(fileID, "vertices", nvertex, &nvdimID);
        }

      dimIDs[0] = ydimID;
      dimIDs[1] = xdimID;

      if ( gridInqXvalsPtr(gridID) )
        {
          cdf_def_var(fileID, xaxisname, (nc_type) xtype, 2, dimIDs, &ncxvarid);
          cdfGridCompress(fileID, ncxvarid, (int)(xdimlen*ydimlen), streamptr->filetype, streamptr->comptype);

          if ( (len = strlen(xstdname)) )
            cdf_put_att_text(fileID, ncxvarid, "standard_name", len, xstdname);
          if ( (len = strlen(xlongname)) )
            cdf_put_att_text(fileID, ncxvarid, "long_name", len, xlongname);
          if ( (len = strlen(xunits)) )
            cdf_put_att_text(fileID, ncxvarid, "units", len, xunits);

          /* attribute for Panoply */
          cdf_put_att_text(fileID, ncxvarid, "_CoordinateAxisType", 3, "Lon");

          if ( gridInqXboundsPtr(gridID) && nvdimID != UNDEFID )
            {
              strcat(xaxisname, "_vertices");
              dimIDs[0] = ydimID;
              dimIDs[1] = xdimID;
              dimIDs[2] = nvdimID;
              cdf_def_var(fileID, xaxisname, (nc_type) xtype, 3, dimIDs, &ncbxvarid);
              cdfGridCompress(fileID, ncbxvarid, (int)(xdimlen*ydimlen), streamptr->filetype, streamptr->comptype);

              cdf_put_att_text(fileID, ncxvarid, "bounds", strlen(xaxisname), xaxisname);
            }
        }

      if ( gridInqYvalsPtr(gridID) )
        {
          cdf_def_var(fileID, yaxisname, (nc_type) xtype, 2, dimIDs, &ncyvarid);
          cdfGridCompress(fileID, ncyvarid, (int)(xdimlen*ydimlen), streamptr->filetype, streamptr->comptype);

          if ( (len = strlen(ystdname)) )
            cdf_put_att_text(fileID, ncyvarid, "standard_name", len, ystdname);
          if ( (len = strlen(ylongname)) )
            cdf_put_att_text(fileID, ncyvarid, "long_name", len, ylongname);
          if ( (len = strlen(yunits)) )
            cdf_put_att_text(fileID, ncyvarid, "units", len, yunits);

          /* attribute for Panoply */
          cdf_put_att_text(fileID, ncyvarid, "_CoordinateAxisType", 3, "Lat");

          if ( gridInqYboundsPtr(gridID) && nvdimID != UNDEFID )
            {
              strcat(yaxisname, "_vertices");
              dimIDs[0] = ydimID;
              dimIDs[1] = xdimID;
              dimIDs[2] = nvdimID;
              cdf_def_var(fileID, yaxisname, (nc_type) xtype, 3, dimIDs, &ncbyvarid);
              cdfGridCompress(fileID, ncbyvarid, (int)(xdimlen*ydimlen), streamptr->filetype, streamptr->comptype);

              cdf_put_att_text(fileID, ncyvarid, "bounds", strlen(yaxisname), yaxisname);
            }
        }

      if ( gridInqAreaPtr(gridID) )
        {
          char yaxisname[] = "cell_area";
          char units[] = "m2";
          char longname[] = "area of grid cell";
          char stdname[] = "cell_area";

          cdf_def_var(fileID, yaxisname, (nc_type) xtype, 2, dimIDs, &ncavarid);

          cdf_put_att_text(fileID, ncavarid, "standard_name", strlen(stdname), stdname);
          cdf_put_att_text(fileID, ncavarid, "long_name", strlen(longname), longname);
          cdf_put_att_text(fileID, ncavarid, "units", strlen(units), units);
        }

      cdf_enddef(fileID);
      streamptr->ncmode = 2;

      if ( ncxvarid  != UNDEFID ) cdf_put_var_double(fileID, ncxvarid,  gridInqXvalsPtr(gridID));
      if ( ncbxvarid != UNDEFID ) cdf_put_var_double(fileID, ncbxvarid, gridInqXboundsPtr(gridID));
      if ( ncyvarid  != UNDEFID ) cdf_put_var_double(fileID, ncyvarid,  gridInqYvalsPtr(gridID));
      if ( ncbyvarid != UNDEFID ) cdf_put_var_double(fileID, ncbyvarid, gridInqYboundsPtr(gridID));
      if ( ncavarid  != UNDEFID ) cdf_put_var_double(fileID, ncavarid,  gridInqAreaPtr(gridID));
    }

  streamptr->xdimID[gridindex] = xdimID;
  streamptr->ydimID[gridindex] = ydimID;
  streamptr->ncxvarID[gridindex] = ncxvarid;
  streamptr->ncyvarID[gridindex] = ncyvarid;
  streamptr->ncavarID[gridindex] = ncavarid;
}


static
void cdfDefRgrid(stream_t *streamptr, int gridID)
{
  char axisname[7] = "rgridX";
  int index, iz = 0;
  int gridID0, gridtype0, gridindex;
  int dimID = UNDEFID;
  int ngrids;
  int fileID;
  int vlistID;
  int lwarn = TRUE;

  vlistID = streamptr->vlistID;
  fileID  = streamptr->fileID;

  ngrids = vlistNgrids(vlistID);

  size_t dimlen = (size_t)gridInqSize(gridID);

  for ( index = 0; index < ngrids; index++ )
    {
      if ( streamptr->xdimID[index] != UNDEFID )
        {
          gridID0 = vlistGrid(vlistID, index);
          gridtype0 = gridInqType(gridID0);
          if ( gridtype0 == GRID_GAUSSIAN_REDUCED )
            {
              size_t dimlen0 = (size_t)gridInqSize(gridID0);

              if ( dimlen == dimlen0 )
                {
                  dimID = streamptr->xdimID[index];
                  break;
                }
              else
                iz++;
            }
        }
    }

  if ( dimID == UNDEFID )
    {
      if ( lwarn )
        {
          Warning("Creating a netCDF file with data on a gaussian reduced grid.");
          Warning("The further processing of the resulting file is unsupported!");
          lwarn = FALSE;
        }

      if ( iz == 0 ) axisname[5] = '\0';
      else           sprintf(&axisname[5], "%1d", iz+1);

      if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

      cdf_def_dim(fileID, axisname, dimlen, &dimID);

      cdf_enddef(fileID);
      streamptr->ncmode = 2;
    }

  gridindex = vlistGridIndex(vlistID, gridID);
  streamptr->xdimID[gridindex] = dimID;
}


static
void cdfDefGdim(stream_t *streamptr, int gridID)
{
  int index, iz = 0;
  int gridID0, gridtype0, gridindex;
  int dimID = UNDEFID;
  int ngrids;
  int fileID;
  int vlistID;

  vlistID = streamptr->vlistID;
  fileID  = streamptr->fileID;

  ngrids = vlistNgrids(vlistID);

  size_t dimlen = (size_t)gridInqSize(gridID);

  if ( gridInqYsize(gridID) == 0 )
    for ( index = 0; index < ngrids; index++ )
      {
        if ( streamptr->xdimID[index] != UNDEFID )
          {
            gridID0 = vlistGrid(vlistID, index);
            gridtype0 = gridInqType(gridID0);
            if ( gridtype0 == GRID_GENERIC )
              {
                size_t dimlen0 = (size_t)gridInqSize(gridID0);
                if ( dimlen == dimlen0 )
                  {
                    dimID = streamptr->xdimID[index];
                    break;
                  }
                else
                  iz++;
              }
          }
      }

  if ( gridInqXsize(gridID) == 0 )
    for ( index = 0; index < ngrids; index++ )
      {
        if ( streamptr->ydimID[index] != UNDEFID )
          {
            gridID0 = vlistGrid(vlistID, index);
            gridtype0 = gridInqType(gridID0);
            if ( gridtype0 == GRID_GENERIC )
              {
                size_t dimlen0 = (size_t)gridInqSize(gridID0);
                if ( dimlen == dimlen0 )
                  {
                    dimID = streamptr->ydimID[index];
                    break;
                  }
                else
                  iz++;
              }
          }
      }

  if ( dimID == UNDEFID )
    {
      char axisname[CDI_MAX_NAME];
      strcpy(axisname, "gsize");

      checkGridName('D', axisname, fileID, vlistID, gridID, ngrids, 'X');
      /*
      if ( iz == 0 ) axisname[5] = '\0';
      else           sprintf(&axisname[5], "%1d", iz+1);
      */
      if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

      //printf("axisname, dimlen %s %d\n", axisname, dimlen);
      cdf_def_dim(fileID, axisname, dimlen, &dimID);

      cdf_enddef(fileID);
      streamptr->ncmode = 2;
    }

  gridindex = vlistGridIndex(vlistID, gridID);
  streamptr->xdimID[gridindex] = dimID;
}


static
void cdfDefGridReference(stream_t *streamptr, int gridID)
{
  int fileID  = streamptr->fileID;
  int number = gridInqNumber(gridID);

  if ( number > 0 )
    {
      cdf_put_att_int(fileID, NC_GLOBAL, "number_of_grid_used", NC_INT, 1, &number);
    }

  if ( gridInqReference(gridID, NULL) )
    {
      char gridfile[8912];
      gridInqReference(gridID, gridfile);

      if ( gridfile[0] != 0 )
        cdf_put_att_text(fileID, NC_GLOBAL, "grid_file_uri", strlen(gridfile), gridfile);
    }
}

static
void cdfDefGridUUID(stream_t *streamptr, int gridID)
{
  unsigned char uuidOfHGrid[CDI_UUID_SIZE];

  gridInqUUID(gridID, uuidOfHGrid);
  if ( !cdiUUIDIsNull(uuidOfHGrid) )
    {
      char uuidOfHGridStr[37];
      uuid2str(uuidOfHGrid, uuidOfHGridStr);
      if ( uuidOfHGridStr[0] != 0 && strlen(uuidOfHGridStr) == 36 )
        {
          int fileID  = streamptr->fileID;
          //if ( streamptr->ncmode == 2 ) cdf_redef(fileID);
          cdf_put_att_text(fileID, NC_GLOBAL, "uuidOfHGrid", 36, uuidOfHGridStr);
          //if ( streamptr->ncmode == 2 ) cdf_enddef(fileID);
        }
    }
}

static
void cdfDefZaxisUUID(stream_t *streamptr, int zaxisID)
{
  unsigned char uuidOfVGrid[CDI_UUID_SIZE];
  zaxisInqUUID(zaxisID, uuidOfVGrid);

  if ( uuidOfVGrid[0] != 0 )
    {
      char uuidOfVGridStr[37];
      uuid2str(uuidOfVGrid, uuidOfVGridStr);
      if ( uuidOfVGridStr[0] != 0 && strlen(uuidOfVGridStr) == 36 )
        {
          int fileID  = streamptr->fileID;
          if ( streamptr->ncmode == 2 ) cdf_redef(fileID);
          cdf_put_att_text(fileID, NC_GLOBAL, "uuidOfVGrid", 36, uuidOfVGridStr);
          if ( streamptr->ncmode == 2 ) cdf_enddef(fileID);
        }
    }
}

static
void cdfDefUnstructured(stream_t *streamptr, int gridID)
{
  char xunits[CDI_MAX_NAME];
  char xlongname[CDI_MAX_NAME];
  char xstdname[CDI_MAX_NAME];
  char yunits[CDI_MAX_NAME];
  char ylongname[CDI_MAX_NAME];
  char ystdname[CDI_MAX_NAME];
  char xaxisname[CDI_MAX_NAME];
  char yaxisname[CDI_MAX_NAME];
  int index;
  int gridID0, gridtype0, gridindex;
  int dimID = UNDEFID;
  int ngrids;
  int fileID;
  size_t len;
  int ncxvarid = UNDEFID, ncyvarid = UNDEFID;
  int ncbxvarid = UNDEFID, ncbyvarid = UNDEFID, ncavarid = UNDEFID;
  int nvdimID = UNDEFID;
  int vlistID;
  int xtype = NC_DOUBLE;

  if ( gridInqPrec(gridID) == DATATYPE_FLT32 ) xtype = NC_FLOAT;

  vlistID = streamptr->vlistID;
  fileID  = streamptr->fileID;

  ngrids = vlistNgrids(vlistID);

  size_t dimlen = (size_t)gridInqSize(gridID);
  gridindex = vlistGridIndex(vlistID, gridID);

  gridInqXname(gridID, xaxisname);
  gridInqXlongname(gridID, xlongname);
  gridInqXstdname(gridID, xstdname);
  gridInqXunits(gridID, xunits);
  gridInqYname(gridID, yaxisname);
  gridInqYlongname(gridID, ylongname);
  gridInqYstdname(gridID, ystdname);
  gridInqYunits(gridID, yunits);

  for ( index = 0; index < ngrids; index++ )
    {
      if ( streamptr->xdimID[index] != UNDEFID )
        {
          gridID0 = vlistGrid(vlistID, index);
          gridtype0 = gridInqType(gridID0);
          if ( gridtype0 == GRID_UNSTRUCTURED )
            {
              size_t dimlen0 = (size_t)gridInqSize(gridID0);
              if ( dimlen == dimlen0 )
		if ( gridInqNvertex(gridID0) == gridInqNvertex(gridID) &&
		     IS_EQUAL(gridInqXval(gridID0, 0), gridInqXval(gridID, 0)) &&
                     IS_EQUAL(gridInqXval(gridID0, (int)dimlen-1), gridInqXval(gridID, (int)dimlen-1)) )
		  {
		    dimID = streamptr->xdimID[index];
                    ncxvarid = streamptr->ncxvarID[index];
                    ncyvarid = streamptr->ncyvarID[index];
                    ncavarid = streamptr->ncavarID[index];
		    break;
		  }
            }
        }
    }

  if ( dimID == UNDEFID )
    {
      char axisname[CDI_MAX_NAME];
      char vertname[CDI_MAX_NAME];
      strcpy(axisname, "ncells");
      strcpy(vertname, "vertices");

      checkGridName('V', xaxisname, fileID, vlistID, gridID, ngrids, 'X');
      checkGridName('V', yaxisname, fileID, vlistID, gridID, ngrids, 'Y');
      checkGridName('D', axisname, fileID, vlistID, gridID, ngrids, 'X');
      checkGridName('D', vertname, fileID, vlistID, gridID, ngrids, 'X');

      if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

      cdf_def_dim(fileID, axisname, dimlen, &dimID);

      size_t nvertex = (size_t)gridInqNvertex(gridID);
      if ( nvertex > 0 ) cdf_def_dim(fileID, vertname, nvertex, &nvdimID);

      cdfDefGridReference(streamptr, gridID);

      cdfDefGridUUID(streamptr, gridID);

      if ( gridInqXvalsPtr(gridID) )
        {
          cdf_def_var(fileID, xaxisname, (nc_type) xtype, 1, &dimID, &ncxvarid);
          cdfGridCompress(fileID, ncxvarid, (int)dimlen, streamptr->filetype, streamptr->comptype);

          if ( (len = strlen(xstdname)) )
            cdf_put_att_text(fileID, ncxvarid, "standard_name", len, xstdname);
          if ( (len = strlen(xlongname)) )
            cdf_put_att_text(fileID, ncxvarid, "long_name", len, xlongname);
          if ( (len = strlen(xunits)) )
            cdf_put_att_text(fileID, ncxvarid, "units", len, xunits);

          if ( gridInqXboundsPtr(gridID) && nvdimID != UNDEFID )
            {
              int dimIDs[2];
              dimIDs[0] = dimID;
              dimIDs[1] = nvdimID;
              strcat(xaxisname, "_vertices");
              cdf_def_var(fileID, xaxisname, (nc_type) xtype, 2, dimIDs, &ncbxvarid);
              cdfGridCompress(fileID, ncbxvarid, (int)dimlen, streamptr->filetype, streamptr->comptype);

              cdf_put_att_text(fileID, ncxvarid, "bounds", strlen(xaxisname), xaxisname);
            }
        }

      if ( gridInqYvalsPtr(gridID) )
        {
          cdf_def_var(fileID, yaxisname, (nc_type) xtype, 1, &dimID, &ncyvarid);
          cdfGridCompress(fileID, ncyvarid, (int)dimlen, streamptr->filetype, streamptr->comptype);

          if ( (len = strlen(ystdname)) )
            cdf_put_att_text(fileID, ncyvarid, "standard_name", len, ystdname);
          if ( (len = strlen(ylongname)) )
            cdf_put_att_text(fileID, ncyvarid, "long_name", len, ylongname);
          if ( (len = strlen(yunits)) )
            cdf_put_att_text(fileID, ncyvarid, "units", len, yunits);

          if ( gridInqYboundsPtr(gridID) && nvdimID != UNDEFID )
            {
              int dimIDs[2];
              dimIDs[0] = dimID;
              dimIDs[1] = nvdimID;
              strcat(yaxisname, "_vertices");
              cdf_def_var(fileID, yaxisname, (nc_type) xtype, 2, dimIDs, &ncbyvarid);
              cdfGridCompress(fileID, ncbyvarid, (int)dimlen, streamptr->filetype, streamptr->comptype);

              cdf_put_att_text(fileID, ncyvarid, "bounds", strlen(yaxisname), yaxisname);
            }
        }

      if ( gridInqAreaPtr(gridID) )
        {
          char yaxisname[] = "cell_area";
          char units[] = "m2";
          char longname[] = "area of grid cell";
          char stdname[] = "cell_area";

          cdf_def_var(fileID, yaxisname, (nc_type) xtype, 1, &dimID, &ncavarid);

          cdf_put_att_text(fileID, ncavarid, "standard_name", strlen(stdname), stdname);
          cdf_put_att_text(fileID, ncavarid, "long_name", strlen(longname), longname);
          cdf_put_att_text(fileID, ncavarid, "units", strlen(units), units);
        }

      cdf_enddef(fileID);
      streamptr->ncmode = 2;

      if ( ncxvarid  != UNDEFID ) cdf_put_var_double(fileID, ncxvarid,  gridInqXvalsPtr(gridID));
      if ( ncbxvarid != UNDEFID ) cdf_put_var_double(fileID, ncbxvarid, gridInqXboundsPtr(gridID));
      if ( ncyvarid  != UNDEFID ) cdf_put_var_double(fileID, ncyvarid,  gridInqYvalsPtr(gridID));
      if ( ncbyvarid != UNDEFID ) cdf_put_var_double(fileID, ncbyvarid, gridInqYboundsPtr(gridID));
      if ( ncavarid  != UNDEFID ) cdf_put_var_double(fileID, ncavarid,  gridInqAreaPtr(gridID));
    }

  streamptr->xdimID[gridindex] = dimID;
  streamptr->ncxvarID[gridindex] = ncxvarid;
  streamptr->ncyvarID[gridindex] = ncyvarid;
  streamptr->ncavarID[gridindex] = ncavarid;
}


static
void cdfDefVCT(stream_t *streamptr, int zaxisID)
{
  int type;

  type = zaxisInqType(zaxisID);
  if ( type == ZAXIS_HYBRID || type == ZAXIS_HYBRID_HALF )
    {
      int i;
      int fileID;
      int ilev = zaxisInqVctSize(zaxisID)/2;
      int mlev = ilev - 1;
      size_t start;
      size_t count = 1;
      int ncdimid, ncdimid2;
      int hyaiid, hybiid, hyamid, hybmid;
      double mval;
      char tmpname[CDI_MAX_NAME];

      if ( streamptr->vct.ilev > 0 )
        {
          if ( streamptr->vct.ilev != ilev )
            Error("more than one VCT for each file unsupported!");
          return;
        }

      if ( ilev == 0 )
        {
          Warning("VCT missing");
          return;
        }

      fileID = streamptr->fileID;

      if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

      cdf_def_dim(fileID, "nhym", (size_t)mlev, &ncdimid);
      cdf_def_dim(fileID, "nhyi", (size_t)ilev, &ncdimid2);

      streamptr->vct.mlev   = mlev;
      streamptr->vct.ilev   = ilev;
      streamptr->vct.mlevID = ncdimid;
      streamptr->vct.ilevID = ncdimid2;

      cdf_def_var(fileID, "hyai", NC_DOUBLE, 1, &ncdimid2, &hyaiid);
      cdf_def_var(fileID, "hybi", NC_DOUBLE, 1, &ncdimid2, &hybiid);
      cdf_def_var(fileID, "hyam", NC_DOUBLE, 1, &ncdimid,  &hyamid);
      cdf_def_var(fileID, "hybm", NC_DOUBLE, 1, &ncdimid,  &hybmid);

      strcpy(tmpname, "hybrid A coefficient at layer interfaces");
      cdf_put_att_text(fileID, hyaiid, "long_name", strlen(tmpname), tmpname);
      strcpy(tmpname, "Pa");
      cdf_put_att_text(fileID, hyaiid, "units", strlen(tmpname), tmpname);
      strcpy(tmpname, "hybrid B coefficient at layer interfaces");
      cdf_put_att_text(fileID, hybiid, "long_name", strlen(tmpname), tmpname);
      strcpy(tmpname, "1");
      cdf_put_att_text(fileID, hybiid, "units", strlen(tmpname), tmpname);
      strcpy(tmpname, "hybrid A coefficient at layer midpoints");
      cdf_put_att_text(fileID, hyamid, "long_name", strlen(tmpname), tmpname);
      strcpy(tmpname, "Pa");
      cdf_put_att_text(fileID, hyamid, "units", strlen(tmpname), tmpname);
      strcpy(tmpname, "hybrid B coefficient at layer midpoints");
      cdf_put_att_text(fileID, hybmid, "long_name", strlen(tmpname), tmpname);
      strcpy(tmpname, "1");
      cdf_put_att_text(fileID, hybmid, "units", strlen(tmpname), tmpname);

      cdf_enddef(fileID);
      streamptr->ncmode = 2;

      const double *vctptr = zaxisInqVctPtr(zaxisID);

      cdf_put_var_double(fileID, hyaiid, vctptr);
      cdf_put_var_double(fileID, hybiid, vctptr+ilev);

      for ( i = 0; i < mlev; i++ )
        {
          start = (size_t)i;
          mval = (vctptr[i] + vctptr[i+1]) * 0.5;
          cdf_put_vara_double(fileID, hyamid, &start, &count, &mval);
          mval = (vctptr[ilev+i] + vctptr[ilev+i+1]) * 0.5;
          cdf_put_vara_double(fileID, hybmid, &start, &count, &mval);
        }
    }
}


static
void cdfDefZaxis(stream_t *streamptr, int zaxisID)
{
  /*  char zaxisname0[CDI_MAX_NAME]; */
  char axisname[CDI_MAX_NAME];
  char stdname[CDI_MAX_NAME];
  char longname[CDI_MAX_NAME];
  char units[CDI_MAX_NAME];
  char tmpname[CDI_MAX_NAME];
  int index;
  int zaxisID0;
  int dimID = UNDEFID;
  int dimIDs[2];
  int fileID;
  size_t len;
  int ncvarid = UNDEFID, ncbvarid = UNDEFID;
  int nvdimID = UNDEFID;
  int type;
  int nzaxis;
  int ilevel = 0;
  int vlistID;
  int zaxisindex;
  int xtype = NC_DOUBLE;
  int positive;

  if ( zaxisInqPrec(zaxisID) == DATATYPE_FLT32 ) xtype = NC_FLOAT;

  vlistID = streamptr->vlistID;
  fileID  = streamptr->fileID;

  zaxisindex = vlistZaxisIndex(vlistID, zaxisID);

  nzaxis = vlistNzaxis(vlistID);

  size_t dimlen = (size_t)zaxisInqSize(zaxisID);
  type   = zaxisInqType(zaxisID);

  if (dimlen == 1)
    switch (type)
      {
      case ZAXIS_SURFACE:
      case ZAXIS_CLOUD_BASE:
      case ZAXIS_CLOUD_TOP:
      case ZAXIS_ISOTHERM_ZERO:
      case ZAXIS_TOA:
      case ZAXIS_SEA_BOTTOM:
      case ZAXIS_ATMOSPHERE:
      case ZAXIS_MEANSEA:
      case ZAXIS_LAKE_BOTTOM:
      case ZAXIS_SEDIMENT_BOTTOM:
      case ZAXIS_SEDIMENT_BOTTOM_TA:
      case ZAXIS_SEDIMENT_BOTTOM_TW:
      case ZAXIS_MIX_LAYER:
        return;
      }

  zaxisInqName(zaxisID, axisname);
  /*
  for ( index = 0; index < nzaxis; index++ )
    {
      if ( streamptr->zaxisID[index] != UNDEFID )
        {
          zaxisID0 = vlistZaxis(vlistID, index);
          zaxisInqName(zaxisID0, zaxisname0);
          if ( strcmp(zaxisname0, axisname) == 0 ) ilevel++;
        }
    }
  */
  if ( dimID == UNDEFID )
    {
      char axisname0[CDI_MAX_NAME];
      char axisname2[CDI_MAX_NAME];
      int checkname = FALSE;
      int status;

      /* check that the name is not already defined */
      checkname = TRUE;
      ilevel = 0;

      while ( checkname )
        {
          strcpy(axisname2, axisname);
          if ( ilevel ) sprintf(&axisname2[strlen(axisname2)], "_%d", ilevel+1);

          status = nc_inq_varid(fileID, axisname2, &ncvarid);
          if ( status != NC_NOERR )
            {
              if ( ilevel )
                {
                  /* check that the name does not exist for other grids */
                  for ( index = 0; index < nzaxis; index++ )
                    {
                      zaxisID0 = vlistZaxis(vlistID, index);
                      if ( zaxisID != zaxisID0 )
                        {
                          zaxisInqName(zaxisID0, axisname0);
                          if ( strcmp(axisname0, axisname2) == 0 ) break;
                        }
                    }
                  if ( index == nzaxis ) checkname = FALSE;
                }
              else
                {
                  checkname = FALSE;
                }
            }

          if ( checkname ) ilevel++;

          if ( ilevel > 99 ) break;
        }

      if ( ilevel ) sprintf(&axisname[strlen(axisname)], "_%1d", ilevel+1);

      if ( type == ZAXIS_REFERENCE )
	cdfDefZaxisUUID(streamptr, zaxisID);

      if ( type == ZAXIS_HYBRID || type == ZAXIS_HYBRID_HALF )
        {
          if ( type == ZAXIS_HYBRID )
            {
	      if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

	      cdf_def_dim(fileID, axisname, dimlen, &dimID);
	      cdf_def_var(fileID, axisname, (nc_type) xtype, 1, &dimID,  &ncvarid);

	      strcpy(tmpname, "hybrid_sigma_pressure");
	      cdf_put_att_text(fileID, ncvarid, "standard_name", strlen(tmpname), tmpname);
	      strcpy(tmpname, "hybrid level at layer midpoints");
	      cdf_put_att_text(fileID, ncvarid, "long_name", strlen(tmpname), tmpname);
	      strcpy(tmpname, "level");
	      cdf_put_att_text(fileID, ncvarid, "units", strlen(tmpname), tmpname);
	      strcpy(tmpname, "down");
	      cdf_put_att_text(fileID, ncvarid, "positive", strlen(tmpname), tmpname);
	      strcpy(tmpname, "hyam hybm (mlev=hyam+hybm*aps)");
	      cdf_put_att_text(fileID, ncvarid, "formula", strlen(tmpname), tmpname);
	      strcpy(tmpname, "ap: hyam b: hybm ps: aps");
	      cdf_put_att_text(fileID, ncvarid, "formula_terms", strlen(tmpname), tmpname);
	      /*
	      strcpy(tmpname, "ilev");
	      cdf_put_att_text(fileID, ncvarid, "borders", strlen(tmpname), tmpname);
	      */
	      cdf_enddef(fileID);
	      streamptr->ncmode = 2;

	      cdf_put_var_double(fileID, ncvarid, zaxisInqLevelsPtr(zaxisID));
            }

          if ( type == ZAXIS_HYBRID_HALF )
            {
	      if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

	      cdf_def_dim(fileID, axisname, dimlen, &dimID);
	      cdf_def_var(fileID, axisname, (nc_type) xtype, 1, &dimID,  &ncvarid);

	      strcpy(tmpname, "hybrid_sigma_pressure");
	      cdf_put_att_text(fileID, ncvarid, "standard_name", strlen(tmpname), tmpname);
	      strcpy(tmpname, "hybrid level at layer interfaces");
	      cdf_put_att_text(fileID, ncvarid, "long_name", strlen(tmpname), tmpname);
	      strcpy(tmpname, "level");
	      cdf_put_att_text(fileID, ncvarid, "units", strlen(tmpname), tmpname);
	      strcpy(tmpname, "down");
	      cdf_put_att_text(fileID, ncvarid, "positive", strlen(tmpname), tmpname);
	      strcpy(tmpname, "hyai hybi (ilev=hyai+hybi*aps)");
	      cdf_put_att_text(fileID, ncvarid, "formula", strlen(tmpname), tmpname);
	      strcpy(tmpname, "ap: hyai b: hybi ps: aps");
	      cdf_put_att_text(fileID, ncvarid, "formula_terms", strlen(tmpname), tmpname);

	      cdf_enddef(fileID);
	      streamptr->ncmode = 2;

	      cdf_put_var_double(fileID, ncvarid, zaxisInqLevelsPtr(zaxisID));
            }

          cdfDefVCT(streamptr, zaxisID);

          if ( dimID == UNDEFID )
            {
              if ( type == ZAXIS_HYBRID )
                streamptr->zaxisID[zaxisindex] = streamptr->vct.mlevID;
              else
                streamptr->zaxisID[zaxisindex] = streamptr->vct.ilevID;
            }
        }
      else
        {
          if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

          cdf_def_dim(fileID, axisname, dimlen, &dimID);

          zaxisInqLongname(zaxisID, longname);
          zaxisInqUnits(zaxisID, units);
          zaxisInqStdname(zaxisID, stdname);

          cdf_def_var(fileID, axisname, (nc_type) xtype, 1, &dimID, &ncvarid);

          if ( (len = strlen(stdname)) )
            cdf_put_att_text(fileID, ncvarid, "standard_name", len, stdname);
          if ( (len = strlen(longname)) )
            cdf_put_att_text(fileID, ncvarid, "long_name", len, longname);
          if ( (len = strlen(units)) )
            cdf_put_att_text(fileID, ncvarid, "units", len, units);

	  positive = zaxisInqPositive(zaxisID);
	  if ( positive == POSITIVE_UP )
	    {
	      strcpy(tmpname, "up");
	      cdf_put_att_text(fileID, ncvarid, "positive", strlen(tmpname), tmpname);
	    }
	  else if ( positive == POSITIVE_DOWN )
	    {
	      strcpy(tmpname, "down");
	      cdf_put_att_text(fileID, ncvarid, "positive", strlen(tmpname), tmpname);
	    }

          cdf_put_att_text(fileID, ncvarid, "axis", 1, "Z");

	  if ( zaxisInqLbounds(zaxisID, NULL) && zaxisInqUbounds(zaxisID, NULL) )
            {
              size_t nvertex = 2;
	      if ( nc_inq_dimid(fileID, "bnds", &nvdimID) != NC_NOERR )
		cdf_def_dim(fileID, "bnds", nvertex, &nvdimID);

	      if ( nvdimID != UNDEFID )
		{
		  strcat(axisname, "_bnds");
		  dimIDs[0] = dimID;
		  dimIDs[1] = nvdimID;
		  cdf_def_var(fileID, axisname, (nc_type) xtype, 2, dimIDs, &ncbvarid);
		  cdf_put_att_text(fileID, ncvarid, "bounds", strlen(axisname), axisname);
		}
	    }

          cdf_enddef(fileID);
          streamptr->ncmode = 2;

          cdf_put_var_double(fileID, ncvarid, zaxisInqLevelsPtr(zaxisID));

          if ( ncbvarid != UNDEFID )
	    {
	      double *zbounds, *lbounds, *ubounds;

	      lbounds = (double *) malloc(dimlen*sizeof(double));
	      ubounds = (double *) malloc(dimlen*sizeof(double));
	      zbounds = (double *) malloc(2*dimlen*sizeof(double));

	      zaxisInqLbounds(zaxisID, lbounds);
	      zaxisInqUbounds(zaxisID, ubounds);

	      for (size_t i = 0; i < dimlen; ++i )
		{
		  zbounds[2*i  ] = lbounds[i];
		  zbounds[2*i+1] = ubounds[i];
		}

	      cdf_put_var_double(fileID, ncbvarid, zbounds);

	      free(zbounds);
	      free(ubounds);
	      free(lbounds);
	    }
        }
    }

  if ( dimID != UNDEFID )
    streamptr->zaxisID[zaxisindex] = dimID;
}


static
void cdfDefPole(stream_t *streamptr, int gridID)
{
  int fileID;
  int ncvarid = UNDEFID;
  int ncerr;
  double xpole, ypole, angle;
  char varname[] = "rotated_pole";
  char mapname[] = "rotated_latitude_longitude";

  fileID  = streamptr->fileID;

  ypole = gridInqYpole(gridID);
  xpole = gridInqXpole(gridID);
  angle = gridInqAngle(gridID);

  cdf_redef(fileID);

  ncerr = nc_def_var(fileID, varname, (nc_type) NC_CHAR, 0, NULL, &ncvarid);
  if ( ncerr == NC_NOERR )
    {
      cdf_put_att_text(fileID, ncvarid, "grid_mapping_name", strlen(mapname), mapname);
      cdf_put_att_double(fileID, ncvarid, "grid_north_pole_latitude", NC_DOUBLE, 1, &ypole);
      cdf_put_att_double(fileID, ncvarid, "grid_north_pole_longitude", NC_DOUBLE, 1, &xpole);
      if ( IS_NOT_EQUAL(angle, 0) )
        cdf_put_att_double(fileID, ncvarid, "north_pole_grid_longitude", NC_DOUBLE, 1, &angle);
    }

  cdf_enddef(fileID);
}


static
void cdfDefMapping(stream_t *streamptr, int gridID)
{
  int fileID;
  int ncvarid = UNDEFID;
  int ncerr;

  if ( gridInqType(gridID) == GRID_SINUSOIDAL )
    {
      char varname[] = "sinusoidal";
      char mapname[] = "sinusoidal";

      fileID  = streamptr->fileID;

      cdf_redef(fileID);

      ncerr = nc_def_var(fileID, varname, (nc_type) NC_CHAR, 0, NULL, &ncvarid);
      if ( ncerr == NC_NOERR )
        {
          cdf_put_att_text(fileID, ncvarid, "grid_mapping_name", strlen(mapname), mapname);
          /*
          cdf_put_att_double(fileID, ncvarid, "grid_north_pole_latitude", NC_DOUBLE, 1, &ypole);
          cdf_put_att_double(fileID, ncvarid, "grid_north_pole_longitude", NC_DOUBLE, 1, &xpole);
          */
        }

      cdf_enddef(fileID);
    }
  else if ( gridInqType(gridID) == GRID_LAEA )
    {
      char varname[] = "laea";
      char mapname[] = "lambert_azimuthal_equal_area";

      fileID  = streamptr->fileID;

      cdf_redef(fileID);

      ncerr = nc_def_var(fileID, varname, (nc_type) NC_CHAR, 0, NULL, &ncvarid);
      if ( ncerr == NC_NOERR )
        {
          double a, lon_0, lat_0;

          gridInqLaea(gridID, &a, &lon_0, &lat_0);

          cdf_put_att_text(fileID, ncvarid, "grid_mapping_name", strlen(mapname), mapname);
          cdf_put_att_double(fileID, ncvarid, "earth_radius", NC_DOUBLE, 1, &a);
          cdf_put_att_double(fileID, ncvarid, "longitude_of_projection_origin", NC_DOUBLE, 1, &lon_0);
          cdf_put_att_double(fileID, ncvarid, "latitude_of_projection_origin", NC_DOUBLE, 1, &lat_0);
        }

      cdf_enddef(fileID);
    }
  else if ( gridInqType(gridID) == GRID_LCC2 )
    {
      char varname[] = "Lambert_Conformal";
      char mapname[] = "lambert_conformal_conic";

      fileID  = streamptr->fileID;

      cdf_redef(fileID);

      ncerr = nc_def_var(fileID, varname, (nc_type) NC_CHAR, 0, NULL, &ncvarid);
      if ( ncerr == NC_NOERR )
        {
          double radius, lon_0, lat_0, lat_1, lat_2;

          gridInqLcc2(gridID, &radius, &lon_0, &lat_0, &lat_1, &lat_2);

          cdf_put_att_text(fileID, ncvarid, "grid_mapping_name", strlen(mapname), mapname);
          if ( radius > 0 )
            cdf_put_att_double(fileID, ncvarid, "earth_radius", NC_DOUBLE, 1, &radius);
          cdf_put_att_double(fileID, ncvarid, "longitude_of_central_meridian", NC_DOUBLE, 1, &lon_0);
          cdf_put_att_double(fileID, ncvarid, "latitude_of_projection_origin", NC_DOUBLE, 1, &lat_0);
          if ( IS_EQUAL(lat_1, lat_2) )
            cdf_put_att_double(fileID, ncvarid, "standard_parallel", NC_DOUBLE, 1, &lat_1);
          else
            {
              double lat_1_2[2];
              lat_1_2[0] = lat_1;
              lat_1_2[1] = lat_2;
              cdf_put_att_double(fileID, ncvarid, "standard_parallel", NC_DOUBLE, 2, lat_1_2);
            }
        }

      cdf_enddef(fileID);
    }
}


static
void cdfDefGrid(stream_t *streamptr, int gridID)
{
  int gridtype, size;
  int gridindex;
  int vlistID;

  vlistID = streamptr->vlistID;
  gridindex = vlistGridIndex(vlistID, gridID);
  if ( streamptr->xdimID[gridindex] != UNDEFID ) return;

  gridtype = gridInqType(gridID);
  size     = gridInqSize(gridID);

  if ( CDI_Debug )
    Message("gridtype = %d  size = %d", gridtype, size);

  if ( gridtype == GRID_GAUSSIAN ||
       gridtype == GRID_LONLAT   ||
       gridtype == GRID_GENERIC )
    {
      if ( gridtype == GRID_GENERIC )
        {
          if ( size == 1 && gridInqXsize(gridID) == 0 && gridInqYsize(gridID) == 0 )
            {
              /* no grid information */
            }
          else
            {
              int lx = 0, ly = 0;
              if ( gridInqXsize(gridID) > 0 /*&& gridInqXvals(gridID, NULL) > 0*/ )
                {
                  cdfDefXaxis(streamptr, gridID, 1);
                  lx = 1;
                }

              if ( gridInqYsize(gridID) > 0 /*&& gridInqYvals(gridID, NULL) > 0*/ )
                {
                  cdfDefYaxis(streamptr, gridID, 1);
                  ly = 1;
                }

              if ( lx == 0 && ly == 0 ) cdfDefGdim(streamptr, gridID);
            }
        }
      else
        {
          int ndims = 1;
          if ( gridtype == GRID_LONLAT && size == 1 && gridInqHasDims(gridID) == FALSE )
            ndims = 0;

          if ( gridInqXsize(gridID) > 0 ) cdfDefXaxis(streamptr, gridID, ndims);
          if ( gridInqYsize(gridID) > 0 ) cdfDefYaxis(streamptr, gridID, ndims);
        }

      if ( gridIsRotated(gridID) ) cdfDefPole(streamptr, gridID);
    }
  else if ( gridtype == GRID_CURVILINEAR )
    {
      cdfDefCurvilinear(streamptr, gridID);
    }
  else if ( gridtype == GRID_UNSTRUCTURED )
    {
      cdfDefUnstructured(streamptr, gridID);
    }
  else if ( gridtype == GRID_GAUSSIAN_REDUCED )
    {
      cdfDefRgrid(streamptr, gridID);
    }
  else if ( gridtype == GRID_SPECTRAL )
    {
      cdfDefComplex(streamptr, gridID);
      cdfDefSP(streamptr, gridID);
    }
  else if ( gridtype == GRID_FOURIER )
    {
      cdfDefComplex(streamptr, gridID);
      cdfDefFC(streamptr, gridID);
    }
  else if ( gridtype == GRID_TRAJECTORY )
    {
      cdfDefTrajLon(streamptr, gridID);
      cdfDefTrajLat(streamptr, gridID);
    }
  else if ( gridtype == GRID_SINUSOIDAL || gridtype == GRID_LAEA || gridtype == GRID_LCC2 )
    {
      cdfDefXaxis(streamptr, gridID, 1);
      cdfDefYaxis(streamptr, gridID, 1);

      cdfDefMapping(streamptr, gridID);
    }
  /*
  else if ( gridtype == GRID_LCC )
    {
      cdfDefLcc(streamptr, gridID);
    }
  */
  else
    {
      Error("Unsupported grid type: %s", gridNamePtr(gridtype));
    }
}


static
int cdfDefVar(stream_t *streamptr, int varID)
{
  int ncvarid = -1;
  int fileID;
  int xid = UNDEFID, yid = UNDEFID, zid = UNDEFID, tid = UNDEFID;
  size_t xsize = 0, ysize = 0;
  int code, param, gridID, zaxisID;
  int pnum, pcat, pdis;
  char varname[CDI_MAX_NAME];
  const char *name = NULL;
  const char *longname = NULL;
  const char *stdname = NULL;
  const char *units = NULL;
  int dims[4];
  int lchunk = FALSE;
  int chunktype;
  size_t chunks[4] = {0,0,0,0};
  int tableID;
  int ndims = 0;
  int tsteptype;
  int xtype, dtype;
  int gridtype, gridsize;
  int gridindex, zaxisindex;
  int tablenum;
  int vlistID;
  int dimorder[3];
  size_t iax = 0;
  char axis[5];
  int ensID, ensCount, forecast_type;
  int retval;

  fileID  = streamptr->fileID;

  if ( CDI_Debug )
    Message("streamID = %d, fileID = %d, varID = %d", streamptr->self, fileID, varID);

  if ( streamptr->vars[varID].ncvarid != UNDEFID )
    return (streamptr->vars[varID].ncvarid);

  vlistID   = streamptr->vlistID;
  gridID    = vlistInqVarGrid(vlistID, varID);
  zaxisID   = vlistInqVarZaxis(vlistID, varID);
  tsteptype = vlistInqVarTsteptype(vlistID, varID);
  code      = vlistInqVarCode(vlistID, varID);
  param     = vlistInqVarParam(vlistID, varID);
  cdiDecodeParam(param, &pnum, &pcat, &pdis);

  chunktype = vlistInqVarChunkType(vlistID, varID);

  vlistInqVarDimorder(vlistID, varID, &dimorder);

  gridsize  = gridInqSize(gridID);
  if ( gridsize > 1 ) lchunk = TRUE;
  gridtype  = gridInqType(gridID);
  gridindex = vlistGridIndex(vlistID, gridID);
  if ( gridtype != GRID_TRAJECTORY )
    {
      xid = streamptr->xdimID[gridindex];
      yid = streamptr->ydimID[gridindex];
      if ( xid != UNDEFID ) cdf_inq_dimlen(fileID, xid, &xsize);
      if ( yid != UNDEFID ) cdf_inq_dimlen(fileID, yid, &ysize);
    }

  zaxisindex = vlistZaxisIndex(vlistID, zaxisID);
  zid = streamptr->zaxisID[zaxisindex];

  if ( dimorder[0] != 3 ) lchunk = FALSE; /* ZYX and ZXY */

  if ( ((dimorder[0]>0)+(dimorder[1]>0)+(dimorder[2]>0)) < ((xid!=UNDEFID)+(yid!=UNDEFID)+(zid!=UNDEFID)) )
    {
      printf("zid=%d  yid=%d  xid=%d\n", zid, yid, xid);
      Error("Internal problem, dimension order missing!");
    }

  tid = streamptr->basetime.ncdimid;

  if ( tsteptype != TSTEP_CONSTANT )
    {
      if ( tid == UNDEFID ) Error("Internal problem, time undefined!");
      chunks[ndims] = 1;
      dims[ndims++] = tid;
      axis[iax++] = 'T';
    }
  /*
  if ( zid != UNDEFID ) axis[iax++] = 'Z';
  if ( zid != UNDEFID ) chunks[ndims] = 1;
  if ( zid != UNDEFID ) dims[ndims++] = zid;

  if ( yid != UNDEFID ) chunks[ndims] = ysize;
  if ( yid != UNDEFID ) dims[ndims++] = yid;

  if ( xid != UNDEFID ) chunks[ndims] = xsize;
  if ( xid != UNDEFID ) dims[ndims++] = xid;
  */
  for ( int id = 0; id < 3; ++id )
    {
      if ( dimorder[id] == 3 && zid != UNDEFID )
        {
          axis[iax++] = 'Z';
          chunks[ndims] = 1;
          dims[ndims] = zid;
          ndims++;
        }
      else if ( dimorder[id] == 2 && yid != UNDEFID )
        {
          if ( chunktype == CHUNK_LINES )
            chunks[ndims] = 1;
          else
            chunks[ndims] = ysize;
          dims[ndims] = yid;
          ndims++;
        }
      else if ( dimorder[id] == 1 && xid != UNDEFID )
        {
          chunks[ndims] = xsize;
          dims[ndims] = xid;
          ndims++;
        }
    }

  if ( CDI_Debug )
    fprintf(stderr, "chunktype %d  chunks %d %d %d %d\n", chunktype, (int)chunks[0], (int)chunks[1], (int)chunks[2], (int)chunks[3]);

  tableID  = vlistInqVarTable(vlistID, varID);

  name     = vlistInqVarNamePtr(vlistID, varID);
  longname = vlistInqVarLongnamePtr(vlistID, varID);
  stdname  = vlistInqVarStdnamePtr(vlistID, varID);
  units    = vlistInqVarUnitsPtr(vlistID, varID);

  if ( name     == NULL )     name = tableInqParNamePtr(tableID, code);
  if ( longname == NULL ) longname = tableInqParLongnamePtr(tableID, code);
  if ( units    == NULL )    units = tableInqParUnitsPtr(tableID, code);
  if ( name )
    {
      int checkname;
      int iz;
      int status;

      sprintf(varname, "%s", name);

      checkname = TRUE;
      iz = 0;

      while ( checkname )
        {
          if ( iz ) sprintf(varname, "%s_%d", name, iz+1);

          status = nc_inq_varid(fileID, varname, &ncvarid);
          if ( status != NC_NOERR )
            {
              checkname = FALSE;
            }

          if ( checkname ) iz++;

          if ( iz >= CDI_MAX_NAME ) Error("Double entry of variable name '%s'!", name);
        }

      if ( strcmp(name, varname) != 0 )
        {
          if ( iz == 1 )
            Warning("Changed double entry of variable name '%s' to '%s'!", name, varname);
          else
            Warning("Changed multiple entry of variable name '%s' to '%s'!", name, varname);
        }

      name = varname;
    }
  else
    {
      int checkname;
      int iz;
      int status;
      char *varname2;

      if ( code < 0 ) code = -code;
      if ( pnum < 0 ) pnum = -pnum;

      if ( pdis == 255 )
	sprintf(varname, "var%d", code);
      else
	sprintf(varname, "param%d.%d.%d", pnum, pcat, pdis);

      varname2 = varname+strlen(varname);

      checkname = TRUE;
      iz = 0;

      while ( checkname )
        {
          if ( iz ) sprintf(varname2, "_%d", iz+1);

          status = nc_inq_varid(fileID, varname, &ncvarid);
          if ( status != NC_NOERR ) checkname = FALSE;

          if ( checkname ) iz++;

          if ( iz >= CDI_MAX_NAME ) break;
        }

      name = varname;
      code = 0;
      pdis = 255;
    }

  /* if ( streamptr->ncmode == 2 ) cdf_redef(fileID); */

  dtype = vlistInqVarDatatype(vlistID, varID);
  xtype = cdfDefDatatype(dtype, streamptr->filetype);

  cdf_def_var(fileID, name, (nc_type) xtype, ndims, dims, &ncvarid);

#if  defined  (HAVE_NETCDF4)
  if ( lchunk && (streamptr->filetype == FILETYPE_NC4 || streamptr->filetype == FILETYPE_NC4C) )
    {
      if ( chunktype == CHUNK_AUTO )
        retval = nc_def_var_chunking(fileID, ncvarid, NC_CHUNKED, NULL);
      else
        retval = nc_def_var_chunking(fileID, ncvarid, NC_CHUNKED, chunks);

      if ( retval ) Error("nc_def_var_chunking failed, status = %d", retval);
    }
#endif

  if ( streamptr->comptype == COMPRESS_ZIP )
    {
      if ( lchunk && (streamptr->filetype == FILETYPE_NC4 || streamptr->filetype == FILETYPE_NC4C) )
        {
          cdfDefVarDeflate(fileID, ncvarid, streamptr->complevel);
        }
      else
        {
          if ( lchunk )
            {
              static int lwarn = TRUE;

              if ( lwarn )
                {
                  lwarn = FALSE;
                  Warning("Deflate compression is only available for netCDF4!");
                }
            }
        }
    }

  if ( streamptr->comptype == COMPRESS_SZIP )
    {
      if ( lchunk && (streamptr->filetype == FILETYPE_NC4 || streamptr->filetype == FILETYPE_NC4C) )
        {
#if defined (NC_SZIP_NN_OPTION_MASK)
          cdfDefVarSzip(fileID, ncvarid);
#else
          static int lwarn = TRUE;

          if ( lwarn )
            {
              lwarn = FALSE;
              Warning("netCDF4/SZIP compression not available!");
            }
#endif
        }
      else
        {
          static int lwarn = TRUE;

          if ( lwarn )
            {
              lwarn = FALSE;
              Warning("SZIP compression is only available for netCDF4!");
            }
        }
    }

  if ( stdname && *stdname )
    cdf_put_att_text(fileID, ncvarid, "standard_name", strlen(stdname), stdname);

  if ( longname && *longname )
    cdf_put_att_text(fileID, ncvarid, "long_name", strlen(longname), longname);

  if ( units && *units )
    cdf_put_att_text(fileID, ncvarid, "units", strlen(units), units);

  if ( code > 0 && pdis == 255 )
    cdf_put_att_int(fileID, ncvarid, "code", NC_INT, 1, &code);

  if ( pdis != 255 )
    {
      char paramstr[32];
      cdiParamToString(param, paramstr, sizeof(paramstr));
      cdf_put_att_text(fileID, ncvarid, "param", strlen(paramstr), paramstr);
    }

  if ( tableID != UNDEFID )
    {
      tablenum = tableInqNum(tableID);
      if ( tablenum > 0 )
        cdf_put_att_int(fileID, ncvarid, "table", NC_INT, 1, &tablenum);
    }

  if ( gridtype != GRID_GENERIC && gridtype != GRID_LONLAT  && gridtype != GRID_CURVILINEAR )
    {
      size_t len = strlen(gridNamePtr(gridtype));
      if ( len > 0 )
        cdf_put_att_text(fileID, ncvarid, "grid_type", len, gridNamePtr(gridtype));
    }

  if ( gridIsRotated(gridID) )
    {
      char mapping[] = "rotated_pole";
      cdf_put_att_text(fileID, ncvarid, "grid_mapping", strlen(mapping), mapping);
    }

  if ( gridtype == GRID_SINUSOIDAL )
    {
      char mapping[] = "sinusoidal";
      cdf_put_att_text(fileID, ncvarid, "grid_mapping", strlen(mapping), mapping);
    }
  else if ( gridtype == GRID_LAEA )
    {
      char mapping[] = "laea";
      cdf_put_att_text(fileID, ncvarid, "grid_mapping", strlen(mapping), mapping);
    }
  else if ( gridtype == GRID_LCC2 )
    {
      char mapping[] = "Lambert_Conformal";
      cdf_put_att_text(fileID, ncvarid, "grid_mapping", strlen(mapping), mapping);
    }
  else if ( gridtype == GRID_TRAJECTORY )
    {
      cdf_put_att_text(fileID, ncvarid, "coordinates", 9, "tlon tlat" );
    }
  else if ( gridtype == GRID_LONLAT && xid == UNDEFID && yid == UNDEFID && gridsize == 1 )
    {
      char coordinates[CDI_MAX_NAME] = "";
      int ncxvarID, ncyvarID;
      int gridindex;
      size_t len;

      gridindex = vlistGridIndex(vlistID, gridID);
      ncxvarID = streamptr->ncxvarID[gridindex];
      ncyvarID = streamptr->ncyvarID[gridindex];
      if ( ncxvarID != CDI_UNDEFID )
        cdf_inq_varname(fileID, ncxvarID, coordinates);
      len = strlen(coordinates);
      if ( ncyvarID != CDI_UNDEFID )
        {
          if ( len ) coordinates[len++] = ' ';
          cdf_inq_varname(fileID, ncyvarID, coordinates+len);
        }
      len = strlen(coordinates);
      if ( len )
        cdf_put_att_text(fileID, ncvarid, "coordinates", len, coordinates);
    }
  else if ( gridtype == GRID_UNSTRUCTURED || gridtype == GRID_CURVILINEAR )
    {
      char coordinates[CDI_MAX_NAME] = "";
      char cellarea[CDI_MAX_NAME] = "area: ";
      int ncxvarID, ncyvarID, ncavarID;
      int gridindex;
      size_t len;

      gridindex = vlistGridIndex(vlistID, gridID);
      ncxvarID = streamptr->ncxvarID[gridindex];
      ncyvarID = streamptr->ncyvarID[gridindex];
      ncavarID = streamptr->ncavarID[gridindex];
      if ( ncxvarID != CDI_UNDEFID )
        cdf_inq_varname(fileID, ncxvarID, coordinates);
      len = strlen(coordinates);
      if ( ncyvarID != CDI_UNDEFID )
        {
          if ( len ) coordinates[len++] = ' ';
          cdf_inq_varname(fileID, ncyvarID, coordinates+len);
        }
      len = strlen(coordinates);
      if ( len )
        cdf_put_att_text(fileID, ncvarid, "coordinates", len, coordinates);

      if ( ncavarID != CDI_UNDEFID )
        {
          len = strlen(cellarea);
          cdf_inq_varname(fileID, ncavarID, cellarea+len);
          len = strlen(cellarea);
          cdf_put_att_text(fileID, ncvarid, "cell_measures", len, cellarea);
        }

      if ( gridtype == GRID_UNSTRUCTURED )
        {
          int position = gridInqPosition(gridID);
          if ( position > 0 )
            cdf_put_att_int(fileID, ncvarid, "number_of_grid_in_reference", NC_INT, 1, &position);
        }
    }
  else if ( gridtype == GRID_SPECTRAL || gridtype == GRID_FOURIER )
    {
      int gridTruncation = gridInqTrunc(gridID);

      axis[iax++] = '-';
      axis[iax++] = '-';
      cdf_put_att_text(fileID, ncvarid, "axis", iax, axis);
      cdf_put_att_int(fileID, ncvarid, "truncation", NC_INT, 1, &gridTruncation);
    }

  /*  if ( xtype == NC_BYTE || xtype == NC_SHORT || xtype == NC_INT ) */
    {
      int laddoffset, lscalefactor;
      double addoffset, scalefactor;
      int astype = NC_DOUBLE;

      addoffset    = vlistInqVarAddoffset(vlistID, varID);
      scalefactor  = vlistInqVarScalefactor(vlistID, varID);
      laddoffset   = IS_NOT_EQUAL(addoffset, 0);
      lscalefactor = IS_NOT_EQUAL(scalefactor, 1);

      if ( laddoffset || lscalefactor )
        {
          if ( IS_EQUAL(addoffset,   (double) ((float) addoffset)) &&
               IS_EQUAL(scalefactor, (double) ((float) scalefactor)) )
            {
              astype = NC_FLOAT;
            }

          if ( xtype == (int) NC_FLOAT ) astype = NC_FLOAT;

          cdf_put_att_double(fileID, ncvarid, "add_offset",   (nc_type) astype, 1, &addoffset);
          cdf_put_att_double(fileID, ncvarid, "scale_factor", (nc_type) astype, 1, &scalefactor);
        }
    }

  if ( dtype == DATATYPE_UINT8 && xtype == NC_BYTE )
    {
      int validrange[2] = {0, 255};
      cdf_put_att_int(fileID, ncvarid, "valid_range", NC_SHORT, 2, validrange);
      cdf_put_att_text(fileID, ncvarid, "_Unsigned", 4, "true");
    }

  streamptr->vars[varID].ncvarid = ncvarid;

  if ( vlistInqVarMissvalUsed(vlistID, varID) )
    cdfDefVarMissval(streamptr, varID, vlistInqVarDatatype(vlistID, varID), 0);

  if ( zid == -1 )
    {
      if ( zaxisInqType(zaxisID) == ZAXIS_CLOUD_BASE          ||
           zaxisInqType(zaxisID) == ZAXIS_CLOUD_TOP           ||
           zaxisInqType(zaxisID) == ZAXIS_ISOTHERM_ZERO       ||
           zaxisInqType(zaxisID) == ZAXIS_TOA                 ||
           zaxisInqType(zaxisID) == ZAXIS_SEA_BOTTOM          ||
           zaxisInqType(zaxisID) == ZAXIS_LAKE_BOTTOM         ||
           zaxisInqType(zaxisID) == ZAXIS_SEDIMENT_BOTTOM     ||
           zaxisInqType(zaxisID) == ZAXIS_SEDIMENT_BOTTOM_TA  ||
           zaxisInqType(zaxisID) == ZAXIS_SEDIMENT_BOTTOM_TW  ||
           zaxisInqType(zaxisID) == ZAXIS_MIX_LAYER           ||
           zaxisInqType(zaxisID) == ZAXIS_ATMOSPHERE )
        {
          zaxisInqName(zaxisID, varname);
          cdf_put_att_text(fileID, ncvarid, "level_type", strlen(varname), varname);
        }
    }

  if ( vlistInqVarEnsemble( vlistID,  varID, &ensID, &ensCount, &forecast_type ) )
    {
      /* void cdf_put_att_int(  int ncid, int varid, const char *name, nc_type xtype,
	                        size_t len, const int *ip )
       */
	cdf_put_att_int(fileID, ncvarid, "realization", NC_INT, 1, &ensID);
	cdf_put_att_int(fileID, ncvarid, "ensemble_members", NC_INT, 1, &ensCount);
	cdf_put_att_int(fileID, ncvarid, "forecast_init_type", NC_INT, 1, &forecast_type);

#ifdef DBG
	if( DBG )
	  {
	    fprintf( stderr, "cdfDefVar :\n EnsID  %d\n Enscount %d\n Forecast init type %d\n",  ensID,
		     ensCount,  forecast_type );
	  }
#endif
    }

  /* Attributes */
  defineAttributes(vlistID, varID, fileID, ncvarid);

  /* if ( streamptr->ncmode == 2 ) cdf_enddef(fileID); */

  return (ncvarid);
}

static
void scale_add(size_t size, double *data, double addoffset, double scalefactor)
{
  int laddoffset;
  int lscalefactor;

  laddoffset   = IS_NOT_EQUAL(addoffset, 0);
  lscalefactor = IS_NOT_EQUAL(scalefactor, 1);

  if ( laddoffset || lscalefactor )
    {
      for (size_t i = 0; i < size; ++i )
        {
          if ( lscalefactor ) data[i] *= scalefactor;
          if ( laddoffset )   data[i] += addoffset;
        }
    }
}

static
void cdfGetSlapDescription(stream_t *streamptr, int varID, size_t (*start)[4], size_t (*count)[4])
{
  int vlistID = streamptr->vlistID;
  int tsID = streamptr->curTsID;
  int gridID    = vlistInqVarGrid(vlistID, varID);
  int zaxisID   = vlistInqVarZaxis(vlistID, varID);
  int tsteptype = vlistInqVarTsteptype(vlistID, varID);
  int gridindex = vlistGridIndex(vlistID, gridID);

  if ( CDI_Debug ) Message("tsID = %d", tsID);

  int xid = UNDEFID, yid = UNDEFID;
  if ( gridInqType(gridID) == GRID_TRAJECTORY )
    {
      cdfReadGridTraj(streamptr, gridID);
    }
  else
    {
      xid = streamptr->xdimID[gridindex];
      yid = streamptr->ydimID[gridindex];
    }
  int zaxisindex = vlistZaxisIndex(vlistID, zaxisID);
  int zid = streamptr->zaxisID[zaxisindex];

  int ndims = 0;
#define addDimension(startCoord, length) do \
    { \
      (*start)[ndims] = startCoord; \
      (*count)[ndims] = length; \
      ndims++; \
    } while(0)
  if ( tsteptype != TSTEP_CONSTANT ) addDimension((size_t)tsID, 1);
  if ( zid != UNDEFID ) addDimension(0, (size_t)zaxisInqSize(zaxisID));
  if ( yid != UNDEFID ) addDimension(0, (size_t)gridInqYsize(gridID));
  if ( xid != UNDEFID ) addDimension(0, (size_t)gridInqXsize(gridID));
#undef addDimension

  assert(ndims <= (int)(sizeof(*start)/sizeof(**start)));
  assert(ndims <= (int)(sizeof(*count)/sizeof(**count)));

  if ( CDI_Debug )
    for (int idim = 0; idim < ndims; idim++)
      Message("dim = %d  start = %d  count = %d", idim, start[idim], count[idim]);
}

//Scans the data array for missVals, optionally applying first a scale factor and then an offset.
//Returns the number of missing + out-of-range values encountered.
static
size_t cdfDoInputDataTransformationDP(size_t valueCount, double *data, bool haveMissVal, double missVal, double scaleFactor, double offset, double validMin, double validMax)
 {
  const bool haveOffset   = IS_NOT_EQUAL(offset, 0);
  const bool haveScaleFactor = IS_NOT_EQUAL(scaleFactor, 1);
  size_t missValCount = 0;

  if ( IS_EQUAL(validMin, VALIDMISS) ) validMin = DBL_MIN;
  if ( IS_EQUAL(validMax, VALIDMISS) ) validMax = DBL_MAX;

  bool haveRangeCheck = (IS_NOT_EQUAL(validMax, DBL_MAX)) | (IS_NOT_EQUAL(validMin,DBL_MIN));
  assert(!haveRangeCheck || haveMissVal);

  switch ((int)haveMissVal | ((int)haveScaleFactor << 1)
          | ((int)haveOffset << 2) | ((int)haveRangeCheck << 3))
    {
    case 15: /* haveRangeCheck & haveMissVal & haveScaleFactor & haveOffset */
      for ( size_t i = 0; i < valueCount; i++ )
        {
          int outOfRange = data[i] < validMin || data[i] > validMax;
          int isMissVal = DBL_IS_EQUAL(data[i], missVal);
          missValCount += (size_t)(outOfRange | isMissVal);
          data[i] = outOfRange ? missVal
            : isMissVal ? data[i] : data[i] * scaleFactor + offset;
        }
      break;
    case 13: /* haveRangeCheck & haveMissVal & haveOffset */
      for ( size_t i = 0; i < valueCount; i++ )
        {
          int outOfRange = data[i] < validMin || data[i] > validMax;
          int isMissVal = DBL_IS_EQUAL(data[i], missVal);
          missValCount += (size_t)(outOfRange | isMissVal);
          data[i] = outOfRange ? missVal
            : isMissVal ? data[i] : data[i] + offset;
        }
      break;
    case 11: /* haveRangeCheck & haveMissVal & haveScaleFactor */
      for ( size_t i = 0; i < valueCount; i++ )
        {
          int outOfRange = data[i] < validMin || data[i] > validMax;
          int isMissVal = DBL_IS_EQUAL(data[i], missVal);
          missValCount += (size_t)(outOfRange | isMissVal);
          data[i] = outOfRange ? missVal
            : isMissVal ? data[i] : data[i] * scaleFactor;
        }
      break;
    case 9: /* haveRangeCheck & haveMissVal */
      for ( size_t i = 0; i < valueCount; i++ )
        {
          int outOfRange = data[i] < validMin || data[i] > validMax;
          int isMissVal = DBL_IS_EQUAL(data[i], missVal);
          missValCount += (size_t)(outOfRange | isMissVal);
          data[i] = outOfRange ? missVal : data[i];
        }
      break;
    case 7: /* haveMissVal & haveScaleFactor & haveOffset */
      for ( size_t i = 0; i < valueCount; i++ )
        if ( DBL_IS_EQUAL(data[i], missVal) )
          missValCount++;
        else
          data[i] = data[i] * scaleFactor + offset;
      break;
    case 6: /* haveOffset & haveScaleFactor */
      for ( size_t i = 0; i < valueCount; i++ )
        data[i] = data[i] * scaleFactor + offset;
      break;
    case 5: /* haveMissVal & haveOffset */
      for ( size_t i = 0; i < valueCount; i++ )
        if ( DBL_IS_EQUAL(data[i], missVal) )
          missValCount++;
        else
          data[i] += offset;
      break;
    case 4: /* haveOffset */
      for ( size_t i = 0; i < valueCount; i++ )
        data[i] += offset;
      break;
    case 3: /* haveMissVal & haveScaleFactor */
      for ( size_t i = 0; i < valueCount; i++ )
        if ( DBL_IS_EQUAL(data[i], missVal) )
          missValCount++;
        else
          data[i] *= scaleFactor;
      break;
    case 2: /* haveScaleFactor */
      for ( size_t i = 0; i < valueCount; i++ )
        data[i] *= scaleFactor;
      break;
    case 1: /* haveMissVal */
      for ( size_t i = 0; i < valueCount; i++ )
        missValCount += (unsigned)DBL_IS_EQUAL(data[i], missVal);
      break;
    }

  return missValCount;
}

static
size_t cdfDoInputDataTransformationSP(size_t valueCount, float *data, bool haveMissVal, double missVal, double scaleFactor, double offset, double validMin, double validMax)
 {
  const bool haveOffset   = IS_NOT_EQUAL(offset, 0);
  const bool haveScaleFactor = IS_NOT_EQUAL(scaleFactor, 1);
  size_t missValCount = 0;

  if ( IS_EQUAL(validMin, VALIDMISS) ) validMin = DBL_MIN;
  if ( IS_EQUAL(validMax, VALIDMISS) ) validMax = DBL_MAX;

  bool haveRangeCheck = (IS_NOT_EQUAL(validMax, DBL_MAX)) | (IS_NOT_EQUAL(validMin,DBL_MIN));
  assert(!haveRangeCheck || haveMissVal);

  switch ((int)haveMissVal | ((int)haveScaleFactor << 1)
          | ((int)haveOffset << 2) | ((int)haveRangeCheck << 3))
    {
    case 15: /* haveRangeCheck & haveMissVal & haveScaleFactor & haveOffset */
      for ( size_t i = 0; i < valueCount; i++ )
        {
          int outOfRange = data[i] < validMin || data[i] > validMax;
          int isMissVal = DBL_IS_EQUAL(data[i], missVal);
          missValCount += (size_t)(outOfRange | isMissVal);
          data[i] = outOfRange ? (float)missVal
            : isMissVal ? data[i] : (float)(data[i] * scaleFactor + offset);
        }
      break;
    case 13: /* haveRangeCheck & haveMissVal & haveOffset */
      for ( size_t i = 0; i < valueCount; i++ )
        {
          int outOfRange = data[i] < validMin || data[i] > validMax;
          int isMissVal = DBL_IS_EQUAL(data[i], missVal);
          missValCount += (size_t)(outOfRange | isMissVal);
          data[i] = outOfRange ? (float)missVal
            : isMissVal ? data[i] : (float)(data[i] + offset);
        }
      break;
    case 11: /* haveRangeCheck & haveMissVal & haveScaleFactor */
      for ( size_t i = 0; i < valueCount; i++ )
        {
          int outOfRange = data[i] < validMin || data[i] > validMax;
          int isMissVal = DBL_IS_EQUAL(data[i], missVal);
          missValCount += (size_t)(outOfRange | isMissVal);
          data[i] = outOfRange ? (float)missVal
            : isMissVal ? data[i] : (float)(data[i] * scaleFactor);
        }
      break;
    case 9: /* haveRangeCheck & haveMissVal */
      for ( size_t i = 0; i < valueCount; i++ )
        {
          int outOfRange = data[i] < validMin || data[i] > validMax;
          int isMissVal = DBL_IS_EQUAL(data[i], missVal);
          missValCount += (size_t)(outOfRange | isMissVal);
          data[i] = outOfRange ? (float)missVal : data[i];
        }
      break;
    case 7: /* haveMissVal & haveScaleFactor & haveOffset */
      for ( size_t i = 0; i < valueCount; i++ )
        if ( DBL_IS_EQUAL(data[i], missVal) )
          missValCount++;
        else
          data[i] = (float)(data[i] * scaleFactor + offset);
      break;
    case 6: /* haveOffset & haveScaleFactor */
      for ( size_t i = 0; i < valueCount; i++ )
        data[i] = (float)(data[i] * scaleFactor + offset);
      break;
    case 5: /* haveMissVal & haveOffset */
      for ( size_t i = 0; i < valueCount; i++ )
        if ( DBL_IS_EQUAL(data[i], missVal) )
          missValCount++;
        else
          data[i] = (float)(data[i] + offset);
      break;
    case 4: /* haveOffset */
      for ( size_t i = 0; i < valueCount; i++ )
        data[i] = (float)(data[i] + offset);
      break;
    case 3: /* haveMissVal & haveScaleFactor */
      for ( size_t i = 0; i < valueCount; i++ )
        if ( DBL_IS_EQUAL(data[i], missVal) )
          missValCount++;
        else
          data[i] = (float)(data[i] * scaleFactor);
      break;
    case 2: /* haveScaleFactor */
      for ( size_t i = 0; i < valueCount; i++ )
        data[i] = (float)(data[i] * scaleFactor);
      break;
    case 1: /* haveMissVal */
      for ( size_t i = 0; i < valueCount; i++ )
        missValCount += (unsigned)DBL_IS_EQUAL(data[i], missVal);
      break;
    }

  return missValCount;
}

static void
cdf_write_var_data(int fileID, int vlistID, int varID, int ncvarid, int dtype, long nvals, size_t xsize, size_t ysize, int swapxy, size_t *start, size_t *count, int memtype, const void *data, int nmiss)
{
  const double *pdata_dp = (const double *) data;
  double *mdata_dp = NULL;
  double *sdata_dp = NULL;
  const float *pdata_sp = (const float *) data;
  float *mdata_sp = NULL;
  float *sdata_sp = NULL;
  extern int CDF_Debug;

  /*  if ( dtype == DATATYPE_INT8 || dtype == DATATYPE_INT16 || dtype == DATATYPE_INT32 ) */
    {
      bool laddoffset, lscalefactor;
      double addoffset, scalefactor;
      double missval;

      addoffset    = vlistInqVarAddoffset(vlistID, varID);
      scalefactor  = vlistInqVarScalefactor(vlistID, varID);
      laddoffset   = IS_NOT_EQUAL(addoffset, 0);
      lscalefactor = IS_NOT_EQUAL(scalefactor, 1);

      missval      = vlistInqVarMissval(vlistID, varID);

      if ( laddoffset || lscalefactor )
        {
          if ( memtype == MEMTYPE_FLOAT )
            {
              mdata_sp = (float *) malloc((size_t)nvals*sizeof(float));
              memcpy(mdata_sp, pdata_sp, (size_t)nvals * sizeof (float));
              pdata_sp = mdata_sp;

              if ( nmiss > 0 )
                {
                  for (long i = 0; i < nvals; i++ )
                    {
                      double temp = mdata_sp[i];
                      if ( !DBL_IS_EQUAL(temp, missval) )
                        {
                          if ( laddoffset )   temp -= addoffset;
                          if ( lscalefactor ) temp /= scalefactor;
                          mdata_sp[i] = (float)temp;
                        }
                    }
                }
              else
                {
                  for (long i = 0; i < nvals; i++ )
                    {
                      double temp = mdata_sp[i];
                      if ( laddoffset )   temp -= addoffset;
                      if ( lscalefactor ) temp /= scalefactor;
                      mdata_sp[i] = (float)temp;
                    }
                }
            }
          else
            {
              mdata_dp = (double *) malloc((size_t)nvals * sizeof(double));
              memcpy(mdata_dp, pdata_dp, (size_t)nvals * sizeof(double));
              pdata_dp = mdata_dp;

              if ( nmiss > 0 )
                {
                  for (long i = 0; i < nvals; i++ )
                    {
                      if ( !DBL_IS_EQUAL(mdata_dp[i], missval) )
                        {
                          if ( laddoffset )   mdata_dp[i] -= addoffset;
                          if ( lscalefactor ) mdata_dp[i] /= scalefactor;
                        }
                    }
                }
              else
                {
                  for (long i = 0; i < nvals; i++ )
                    {
                      if ( laddoffset )   mdata_dp[i] -= addoffset;
                      if ( lscalefactor ) mdata_dp[i] /= scalefactor;
                    }
                }
            }
        }

      if ( dtype == DATATYPE_UINT8 || dtype == DATATYPE_INT8 ||
           dtype == DATATYPE_INT16 || dtype == DATATYPE_INT32 )
        {
          if ( memtype == MEMTYPE_FLOAT )
            {
              if ( mdata_sp == NULL )
                {
                  mdata_sp = (float *) malloc((size_t)nvals * sizeof(float));
                  memcpy(mdata_sp, pdata_sp, (size_t)nvals * sizeof(float));
                  pdata_sp = mdata_sp;
                }

              for (long i = 0; i < nvals; i++ ) mdata_sp[i] = roundf(mdata_sp[i]);

              if ( dtype == DATATYPE_UINT8 )
                {
                  nc_type xtype;
                  cdf_inq_vartype(fileID, ncvarid, &xtype);
                  if ( xtype == NC_BYTE )
                    {
                      for ( long i = 0; i < nvals; ++i )
                        if ( mdata_sp[i] > 127 ) mdata_sp[i] -= 256;
                    }
                }
            }
          else
            {
              if ( mdata_dp == NULL )
                {
                  mdata_dp = (double *) malloc((size_t)nvals * sizeof (double));
                  memcpy(mdata_dp, pdata_dp, (size_t)nvals * sizeof (double));
                  pdata_dp = mdata_dp;
                }

              for (long i = 0; i < nvals; i++ ) mdata_dp[i] = round(mdata_dp[i]);

              if ( dtype == DATATYPE_UINT8 )
                {
                  nc_type xtype;
                  cdf_inq_vartype(fileID, ncvarid, &xtype);
                  if ( xtype == NC_BYTE )
                    {
                      for (long i = 0; i < nvals; ++i )
                        if ( mdata_dp[i] > 127 ) mdata_dp[i] -= 256;
                    }
                }
            }
        }

      if ( CDF_Debug && memtype != MEMTYPE_FLOAT )
        {
          double fmin, fmax;
          fmin =  1.0e200;
          fmax = -1.0e200;
          for ( long i = 0; i < nvals; ++i )
            {
              if ( !DBL_IS_EQUAL(pdata_dp[i], missval) )
                {
                  if ( pdata_dp[i] < fmin ) fmin = pdata_dp[i];
                  if ( pdata_dp[i] > fmax ) fmax = pdata_dp[i];
                }
            }
          Message("nvals = %d, nmiss = %d, missval = %g, minval = %g, maxval = %g",
                  nvals, nmiss, missval, fmin, fmax);
        }
    }

  if ( swapxy )
    {
      if ( memtype == MEMTYPE_FLOAT )
        {
          /* malloc and the loop imply nvals >= ysize * xsize,
           * but that is not checked and the types don't match */
          sdata_sp = (float *)xmalloc((size_t)nvals * sizeof (float));
          for ( size_t j = 0; j < ysize; ++j )
            for ( size_t i = 0; i < xsize; ++i )
              sdata_sp[i*ysize+j] = pdata_sp[j*xsize+i];
          pdata_sp = sdata_sp;
        }
      else
        {
          sdata_dp = (double *)xmalloc((size_t)nvals * sizeof (double));
          for ( size_t j = 0; j < ysize; ++j )
            for ( size_t i = 0; i < xsize; ++i )
              sdata_dp[i*ysize+j] = pdata_dp[j*xsize+i];
          pdata_dp = sdata_dp;
        }
    }

  if ( memtype == MEMTYPE_FLOAT )
    cdf_put_vara_float(fileID, ncvarid, start, count, pdata_sp);
  else
    cdf_put_vara_double(fileID, ncvarid, start, count, pdata_dp);

  if ( mdata_dp ) free(mdata_dp);
  if ( sdata_dp ) free(sdata_dp);
  if ( mdata_sp ) free(mdata_sp);
  if ( sdata_sp ) free(sdata_sp);
}


void cdf_write_var(stream_t *streamptr, int varID, int memtype, const void *data, int nmiss)
{
  int fileID;
  int gridID;
  int zaxisID;
  int xid = UNDEFID, yid = UNDEFID, zid = UNDEFID;
  int ncvarid;
  size_t xsize = 0, ysize = 0;
  size_t size;
  size_t start[5];
  size_t count[5];
  long nvals;
  int swapxy = FALSE;
  int ndims = 0;
  int idim;
  int tsteptype;
  int gridindex, zaxisindex;
  int dtype;
  int vlistID;

  if ( CDI_Debug ) Message("streamID = %d  varID = %d", streamptr->self, varID);

  vlistID = streamptr->vlistID;
  fileID  = streamptr->fileID;

  long ntsteps = streamptr->ntsteps;
  if ( CDI_Debug ) Message("ntsteps = %ld", ntsteps);

  if ( vlistHasTime(vlistID) ) cdfDefTime(streamptr);

  ncvarid = cdfDefVar(streamptr, varID);

  gridID    = vlistInqVarGrid(vlistID, varID);
  zaxisID   = vlistInqVarZaxis(vlistID, varID);
  tsteptype = vlistInqVarTsteptype(vlistID, varID);

  gridindex = vlistGridIndex(vlistID, gridID);
  if ( gridInqType(gridID) == GRID_TRAJECTORY )
    {
      cdfWriteGridTraj(streamptr, gridID);
    }
  else
    {
      xid = streamptr->xdimID[gridindex];
      yid = streamptr->ydimID[gridindex];
    }

  zaxisindex = vlistZaxisIndex(vlistID, zaxisID);
  zid = streamptr->zaxisID[zaxisindex];

  if ( tsteptype != TSTEP_CONSTANT )
    {
      start[ndims] = (size_t)ntsteps - 1;
      count[ndims] = 1;
      ndims++;
    }
  if ( zid != UNDEFID )
    {
      start[ndims] = 0;
      count[ndims] = (size_t)zaxisInqSize(zaxisID);
      ndims++;
    }
  if ( yid != UNDEFID )
    {
      start[ndims] = 0;
      cdf_inq_dimlen(fileID, yid, &size);
      /*      count[ndims] = gridInqYsize(gridID); */
      count[ndims] = size;
      ndims++;
    }
  if ( xid != UNDEFID )
    {
      start[ndims] = 0;
      cdf_inq_dimlen(fileID, xid, &size);
      /*      count[ndims] = gridInqXsize(gridID); */
      count[ndims] = size;
      ndims++;
    }

  if ( CDI_Debug )
    for (idim = 0; idim < ndims; idim++)
      Message("dim = %d  start = %d  count = %d", idim, start[idim], count[idim]);

  if ( streamptr->ncmode == 1 )
    {
      cdf_enddef(fileID);
      streamptr->ncmode = 2;
    }

  dtype = vlistInqVarDatatype(vlistID, varID);

  if ( nmiss > 0 ) cdfDefVarMissval(streamptr, varID, dtype, 1);

  nvals = gridInqSize(gridID)*zaxisInqSize(zaxisID);

  cdf_write_var_data(fileID, vlistID, varID, ncvarid, dtype, nvals, xsize, ysize, swapxy, start, count, memtype, data, nmiss);

}


void cdf_write_var_chunk(stream_t *streamptr, int varID, int memtype,
                         const int rect[][2], const void *data, int nmiss)
{
  int fileID;
  int gridID;
  int zaxisID;
  int xid = UNDEFID, yid = UNDEFID, zid = UNDEFID;
  int ncvarid;
  size_t xsize = 0, ysize = 0;
  size_t start[5];
  size_t count[5];
  long nvals;
  int swapxy = FALSE;
  int ndims = 0;
  int idim;
  int tsteptype;
  int gridindex, zaxisindex;
  int dtype;
  int vlistID;
  int streamID = streamptr->self;

  if ( CDI_Debug )
    Message("streamID = %d  varID = %d", streamID, varID);

  vlistID = streamInqVlist(streamID);
  fileID  = streamInqFileID(streamID);

  long ntsteps = streamptr->ntsteps;
  if ( CDI_Debug )
    Message("ntsteps = %ld", ntsteps);

  if ( vlistHasTime(vlistID) ) cdfDefTime(streamptr);

  ncvarid = cdfDefVar(streamptr, varID);

  gridID    = vlistInqVarGrid(vlistID, varID);
  zaxisID   = vlistInqVarZaxis(vlistID, varID);
  tsteptype = vlistInqVarTsteptype(vlistID, varID);

  gridindex = vlistGridIndex(vlistID, gridID);
  if ( gridInqType(gridID) == GRID_TRAJECTORY )
    {
      cdfWriteGridTraj(streamptr, gridID);
    }
  else
    {
      xid = streamptr->xdimID[gridindex];
      yid = streamptr->ydimID[gridindex];
    }

  zaxisindex = vlistZaxisIndex(vlistID, zaxisID);
  zid = streamptr->zaxisID[zaxisindex];

  if ( tsteptype != TSTEP_CONSTANT )
    {
      start[ndims] = (size_t)ntsteps - 1;
      count[ndims] = 1;
      ndims++;
    }
  if ( zid != UNDEFID )
    {
      int size = zaxisInqSize(zaxisID);
      xassert(rect[2][0] >= 0 && rect[2][0] <= rect[2][1]
              && rect[2][1] <= size);
      start[ndims] = (size_t)rect[2][0];
      count[ndims] = (size_t)rect[2][1] - (size_t)rect[2][0] + 1;
      ndims++;
    }
  if ( yid != UNDEFID )
    {
      size_t size;
      cdf_inq_dimlen(fileID, yid, &size);
      xassert(rect[1][0] >= 0 && rect[1][0] <= rect[1][1]
              && (size_t)rect[1][1] <= size);
      start[ndims] = (size_t)rect[1][0];
      count[ndims] = (size_t)rect[1][1] - (size_t)rect[1][0] + 1;
      ndims++;
    }
  if ( xid != UNDEFID )
    {
      size_t size;
      cdf_inq_dimlen(fileID, xid, &size);
      xassert(rect[0][0] >= 0 && rect[0][0] <= rect[0][1]
              && (size_t)rect[0][1] <= size);
      start[ndims] = (size_t)rect[0][0];
      count[ndims] = (size_t)rect[0][1] - (size_t)rect[0][0] + 1;
      ndims++;
    }

  if ( CDI_Debug )
    for (idim = 0; idim < ndims; idim++)
      Message("dim = %d  start = %d  count = %d", idim, start[idim], count[idim]);

  if ( streamptr->ncmode == 1 )
    {
      cdf_enddef(fileID);
      streamptr->ncmode = 2;
    }

  dtype = vlistInqVarDatatype(vlistID, varID);

  if ( nmiss > 0 ) cdfDefVarMissval(streamptr, varID, dtype, 1);

  nvals = gridInqSize(gridID)*zaxisInqSize(zaxisID);

  cdf_write_var_data(fileID, vlistID, varID, ncvarid, dtype, nvals,
                     xsize, ysize, swapxy, start, count, memtype, data, nmiss);
}

static
size_t min_size(size_t a, size_t b)
{
  return a < b ? a : b;
}

static
void transpose2dArrayDP(size_t inWidth, size_t inHeight, double* data)
{
  const size_t cacheBlockSize = 256;    // Purely an optimization parameter. Current value of 32 means we are handling 8kB blocks,
                                       // which should be a decent compromise on many architectures.
  double (*temp)[inWidth] = malloc(inHeight*sizeof(*temp));
  double (*out)[inHeight] = (double (*)[inHeight])data;
  memcpy(temp, data, inHeight*sizeof(*temp));
  /*
  for ( size_t y = 0; y < inHeight; ++y )
    for ( size_t x = 0; x < inWidth; ++x )
      out[x][y] = temp[y][x];
  */
  for ( size_t yBlock = 0; yBlock < inHeight; yBlock += cacheBlockSize )
    {
      for ( size_t xBlock = 0; xBlock < inWidth; xBlock += cacheBlockSize )
        {
          for ( size_t y = yBlock, yEnd = min_size(yBlock + cacheBlockSize, inHeight); y < yEnd; y++ )
            {
              for ( size_t x = xBlock, xEnd = min_size(xBlock + cacheBlockSize, inWidth); x < xEnd; x++ )
                {
                  out[x][y] = temp[y][x];
                }
            }
        }
    }

  free(temp);
}

static
void transpose2dArraySP(size_t inWidth, size_t inHeight, float* data)
{
  const size_t cacheBlockSize = 256;    // Purely an optimization parameter. Current value of 32 means we are handling 8kB blocks,
                                       // which should be a decent compromise on many architectures.
  float (*temp)[inWidth] = malloc(inHeight*sizeof(*temp));
  float (*out)[inHeight] = (float (*)[inHeight])data;
  memcpy(temp, data, inHeight*sizeof(*temp));
  /*
  for ( size_t y = 0; y < inHeight; ++y )
    for ( size_t x = 0; x < inWidth; ++x )
      out[x][y] = temp[y][x];
  */
  for ( size_t yBlock = 0; yBlock < inHeight; yBlock += cacheBlockSize )
    {
      for ( size_t xBlock = 0; xBlock < inWidth; xBlock += cacheBlockSize )
        {
          for ( size_t y = yBlock, yEnd = min_size(yBlock + cacheBlockSize, inHeight); y < yEnd; y++ )
            {
              for ( size_t x = xBlock, xEnd = min_size(xBlock + cacheBlockSize, inWidth); x < xEnd; x++ )
                {
                  out[x][y] = temp[y][x];
                }
            }
        }
    }

  free(temp);
}

static
void cdfInqDimIds(stream_t *streamptr, int varId, int (*outDimIds)[3])
{
  int gridId = vlistInqVarGrid(streamptr->vlistID, varId);
  int gridindex = vlistGridIndex(streamptr->vlistID, gridId);

  (*outDimIds)[0] = (*outDimIds)[1] = (*outDimIds)[2] = UNDEFID;
  switch ( gridInqType(gridId) )
    {
      case GRID_TRAJECTORY:
        cdfReadGridTraj(streamptr, gridId);
        break;

      case GRID_UNSTRUCTURED:
        (*outDimIds)[0] = streamptr->xdimID[gridindex];
        break;

      default:
        (*outDimIds)[0] = streamptr->xdimID[gridindex];
        (*outDimIds)[1] = streamptr->ydimID[gridindex];
        break;
    }

  int zaxisID = vlistInqVarZaxis(streamptr->vlistID, varId);
  int zaxisindex = vlistZaxisIndex(streamptr->vlistID, zaxisID);
  (*outDimIds)[2] = streamptr->zaxisID[zaxisindex];
}

static
int cdfGetSkipDim(int fileId, int ncvarid, int (*dimIds)[3])
{
  if((*dimIds)[0] != UNDEFID) return 0;
  if((*dimIds)[1] != UNDEFID) return 0;
  int nvdims;
  cdf_inq_varndims(fileId, ncvarid, &nvdims);
  if(nvdims != 3) return 0;

  int varDimIds[3];
  cdf_inq_vardimid(fileId, ncvarid, varDimIds);
  size_t size = 0;
  if ( (*dimIds)[2] == varDimIds[2] )
    {
      cdf_inq_dimlen(fileId, varDimIds[1], &size);
      if ( size == 1 ) return 1;
    }
  else if ( (*dimIds)[2] == varDimIds[1] )
    {
      cdf_inq_dimlen(fileId, varDimIds[2], &size);
      if ( size == 1 ) return 2;
    }
  return 0;
}


static
void cdfGetSliceSlapDescription(stream_t *streamptr, int varId, int levelId, bool *outSwapXY, size_t (*start)[4], size_t (*count)[4])
{
  int tsID = streamptr->curTsID;
  if ( CDI_Debug ) Message("tsID = %d", tsID);

  int fileId = streamptr->fileID;
  int vlistId = streamptr->vlistID;
  int ncvarid = streamptr->vars[varId].ncvarid;

  int gridId = vlistInqVarGrid(vlistId, varId);
  int tsteptype = vlistInqVarTsteptype(vlistId, varId);
  int gridsize = gridInqSize(gridId);

  streamptr->numvals += gridsize;

  int dimIds[3];    //this array joins the old variables xid, yid, and zid
  cdfInqDimIds(streamptr, varId, &dimIds);

  int skipdim = cdfGetSkipDim(fileId, ncvarid, &dimIds);

  int dimorder[3];
  vlistInqVarDimorder(vlistId, varId, &dimorder);

  *outSwapXY = (dimorder[2] == 2 || dimorder[0] == 1) && dimIds[0] != UNDEFID && dimIds[1] != UNDEFID ;

  int ndims = 0;

#define addDimension(startIndex, extent) do {   \
      (*start)[ndims] = startIndex; \
      (*count)[ndims] = extent; \
      ndims++; \
  } while(0)

  if ( tsteptype != TSTEP_CONSTANT ) addDimension((size_t)tsID, 1);
  if ( skipdim == 1 ) addDimension(0, 1);

  for ( int id = 0; id < 3; ++id )
    {
      size_t size;
      int curDimId = dimIds[dimorder[id]-1];
      if ( curDimId == UNDEFID ) continue;
      switch ( dimorder[id] )
        {
          Error("Internal errror: Malformed dimension order encountered. Please report this bug.\n");
          case 1:
          case 2:
            cdf_inq_dimlen(fileId, curDimId, &size);
            addDimension(0, size);
            break;

          case 3:
            addDimension((size_t)levelId, 1);
            break;

          default:
            Error("Internal errror: Malformed dimension order encountered. Please report this bug.\n");
        }
    }

  if ( skipdim == 2 ) addDimension(0, 1);

  assert(ndims <= (int)(sizeof(*start)/sizeof(**start)));
  assert(ndims <= (int)(sizeof(*count)/sizeof(**count)));

#undef addDimension

  if ( CDI_Debug )
    for (int idim = 0; idim < ndims; idim++)
      Message("dim = %d  start = %d  count = %d", idim, start[idim], count[idim]);

  int nvdims;
  cdf_inq_varndims(fileId, ncvarid, &nvdims);

  if ( nvdims != ndims )
    Error("Internal error, variable %s has an unsupported array structure!", vlistInqVarNamePtr(vlistId, varId));
}

void cdfReadVarDP(stream_t *streamptr, int varID, double *data, int *nmiss)
{
  if ( CDI_Debug ) Message("streamID = %d  varID = %d", streamptr->self, varID);

  int vlistID = streamptr->vlistID;
  int fileID  = streamptr->fileID;

  int ncvarid = streamptr->vars[varID].ncvarid;

  int gridID  = vlistInqVarGrid(vlistID, varID);
  int zaxisID = vlistInqVarZaxis(vlistID, varID);

  size_t start[4];
  size_t count[4];
  cdfGetSlapDescription(streamptr, varID, &start, &count);

  cdf_get_vara_double(fileID, ncvarid, start, count, data);

  size_t size = (size_t)gridInqSize(gridID)*(size_t)zaxisInqSize(zaxisID);
  double missval = vlistInqVarMissval(vlistID, varID);
  const bool haveMissVal = vlistInqVarMissvalUsed(vlistID, varID);
  double validRange[2];
  if (!(haveMissVal && vlistInqVarValidrange(vlistID, varID, validRange)))
    validRange[0] = DBL_MIN, validRange[1] = DBL_MAX;
  double addoffset   = vlistInqVarAddoffset(vlistID, varID);
  double scalefactor = vlistInqVarScalefactor(vlistID, varID);
  size_t nmiss_ = cdfDoInputDataTransformationDP(size, data, haveMissVal, missval, scalefactor, addoffset, validRange[0], validRange[1]);
  assert(nmiss_ <= INT_MAX);
  *nmiss = (int)nmiss_;
}

void cdfReadVarSP(stream_t *streamptr, int varID, float *data, int *nmiss)
{
  if ( CDI_Debug ) Message("streamID = %d  varID = %d", streamptr->self, varID);

  int vlistID = streamptr->vlistID;
  int fileID  = streamptr->fileID;

  int ncvarid = streamptr->vars[varID].ncvarid;

  int gridID  = vlistInqVarGrid(vlistID, varID);
  int zaxisID = vlistInqVarZaxis(vlistID, varID);

  size_t start[4];
  size_t count[4];
  cdfGetSlapDescription(streamptr, varID, &start, &count);

  cdf_get_vara_float(fileID, ncvarid, start, count, data);

  size_t size = (size_t)gridInqSize(gridID) * (size_t)zaxisInqSize(zaxisID);
  double missval = vlistInqVarMissval(vlistID, varID);
  const bool haveMissVal = vlistInqVarMissvalUsed(vlistID, varID);
  double validRange[2];
  if (!(haveMissVal && vlistInqVarValidrange(vlistID, varID, validRange)))
    validRange[0] = DBL_MIN, validRange[1] = DBL_MAX;
  double addoffset   = vlistInqVarAddoffset(vlistID, varID);
  double scalefactor = vlistInqVarScalefactor(vlistID, varID);
  size_t nmiss_ = cdfDoInputDataTransformationSP(size, data, haveMissVal, missval, scalefactor, addoffset, validRange[0], validRange[1]);
  assert(nmiss_ <= INT_MAX);
  *nmiss = (int)nmiss_;
}

void cdfReadVarSliceDP(stream_t *streamptr, int varID, int levelID, double *data, int *nmiss)
{
  size_t start[4];
  size_t count[4];

  if ( CDI_Debug )
    Message("streamID = %d  varID = %d  levelID = %d", streamptr->self, varID, levelID);

  int vlistID = streamptr->vlistID;
  int fileID = streamptr->fileID;

  bool swapxy;
  cdfGetSliceSlapDescription(streamptr, varID, levelID, &swapxy, &start, &count);

  int ncvarid = streamptr->vars[varID].ncvarid;
  int gridId = vlistInqVarGrid(vlistID, varID);
  size_t gridsize = (size_t)gridInqSize(gridId);
  size_t xsize = (size_t)gridInqXsize(gridId);
  size_t ysize = (size_t)gridInqYsize(gridId);

  if ( vlistInqVarDatatype(vlistID, varID) == DATATYPE_FLT32 )
    {
      float *data_fp = (float *)xmalloc(gridsize*sizeof(*data_fp));
      cdf_get_vara_float(fileID, ncvarid, start, count, data_fp);
      for ( size_t i = 0; i < gridsize; i++ )
        data[i] = (double) data_fp[i];
      free(data_fp);
    }
  else if ( vlistInqVarDatatype(vlistID, varID) == DATATYPE_UINT8 )
    {
      nc_type xtype;
      cdf_inq_vartype(fileID, ncvarid, &xtype);
      if ( xtype == NC_BYTE )
        {
          for ( size_t i = 0; i < gridsize; i++ )
            if ( data[i] < 0 ) data[i] += 256;
        }
    }
  else
    {
      cdf_get_vara_double(fileID, ncvarid, start, count, data);
    }

  if ( swapxy ) transpose2dArrayDP(ysize, xsize, data);

  double missval = vlistInqVarMissval(vlistID, varID);
  const bool haveMissVal = vlistInqVarMissvalUsed(vlistID, varID);
  double validRange[2];
  if (!(haveMissVal && vlistInqVarValidrange(vlistID, varID, validRange)))
    validRange[0] = DBL_MIN, validRange[1] = DBL_MAX;
  double addoffset   = vlistInqVarAddoffset(vlistID, varID);
  double scalefactor = vlistInqVarScalefactor(vlistID, varID);
  size_t nmiss_ = cdfDoInputDataTransformationDP(gridsize, data, haveMissVal, missval, scalefactor, addoffset, validRange[0], validRange[1]);
  assert(nmiss_ <= INT_MAX);
  *nmiss = (int)nmiss_;
}


void cdfReadVarSliceSP(stream_t *streamptr, int varID, int levelID, float *data, int *nmiss)
{
  size_t start[4];
  size_t count[4];

  if ( CDI_Debug )
    Message("streamID = %d  varID = %d  levelID = %d", streamptr->self, varID, levelID);

  int vlistID = streamptr->vlistID;
  int fileID = streamptr->fileID;

  bool swapxy;
  cdfGetSliceSlapDescription(streamptr, varID, levelID, &swapxy, &start, &count);

  int ncvarid = streamptr->vars[varID].ncvarid;
  int gridId = vlistInqVarGrid(vlistID, varID);
  size_t gridsize = (size_t)gridInqSize(gridId);
  size_t xsize = (size_t)gridInqXsize(gridId);
  size_t ysize = (size_t)gridInqYsize(gridId);

  if ( vlistInqVarDatatype(vlistID, varID) == DATATYPE_FLT64 )
    {
      double *data_dp = (double *)xmalloc(gridsize*sizeof(*data_dp));
      cdf_get_vara_double(fileID, ncvarid, start, count, data_dp);
      for ( size_t i = 0; i < gridsize; i++ )
        data[i] = (float) data_dp[i];
      free(data_dp);
    }
  else if ( vlistInqVarDatatype(vlistID, varID) == DATATYPE_UINT8 )
    {
      nc_type xtype;
      cdf_inq_vartype(fileID, ncvarid, &xtype);
      if ( xtype == NC_BYTE )
        {
          for ( size_t i = 0; i < gridsize; i++ )
            if ( data[i] < 0 ) data[i] += 256;
        }
    }
  else
    {
      cdf_get_vara_float(fileID, ncvarid, start, count, data);
    }

  if ( swapxy ) transpose2dArraySP(ysize, xsize, data);

  double missval = vlistInqVarMissval(vlistID, varID);
  bool haveMissVal = vlistInqVarMissvalUsed(vlistID, varID);
  double validRange[2];
  if (!(haveMissVal && vlistInqVarValidrange(vlistID, varID, validRange)))
    validRange[0] = DBL_MIN, validRange[1] = DBL_MAX;
  double addoffset   = vlistInqVarAddoffset(vlistID, varID);
  double scalefactor = vlistInqVarScalefactor(vlistID, varID);
  size_t nmiss_ = cdfDoInputDataTransformationSP(gridsize, data, haveMissVal, missval, scalefactor, addoffset, validRange[0], validRange[1]);
  assert(nmiss_ <= INT_MAX);
  *nmiss = (int)nmiss_;
}


void cdf_write_var_slice(stream_t *streamptr, int varID, int levelID, int memtype, const void *data, int nmiss)
{
  size_t xsize = 0, ysize = 0;
  size_t start[5];
  size_t count[5];
  int dimorder[3];
  int xid = UNDEFID, yid = UNDEFID, zid = UNDEFID;

  if ( CDI_Debug ) Message("streamID = %d  varID = %d", streamptr->self, varID);

  int vlistID = streamptr->vlistID;
  int fileID  = streamptr->fileID;

  long ntsteps = streamptr->ntsteps;
  if ( CDI_Debug ) Message("ntsteps = %ld", ntsteps);

  if ( vlistHasTime(vlistID) ) cdfDefTime(streamptr);

  int ncvarid = cdfDefVar(streamptr, varID);

  int gridID    = vlistInqVarGrid(vlistID, varID);
  int zaxisID   = vlistInqVarZaxis(vlistID, varID);
  int tsteptype = vlistInqVarTsteptype(vlistID, varID);
  vlistInqVarDimorder(vlistID, varID, &dimorder);


  if ( gridInqType(gridID) == GRID_TRAJECTORY )
    {
      cdfWriteGridTraj(streamptr, gridID);
    }
  else
    {
      int gridindex = vlistGridIndex(vlistID, gridID);
      xid = streamptr->xdimID[gridindex];
      yid = streamptr->ydimID[gridindex];
    }
  {
    int zaxisindex = vlistZaxisIndex(vlistID, zaxisID);
    zid = streamptr->zaxisID[zaxisindex];
  }

  int swapxy = (dimorder[2] == 2 || dimorder[0] == 1) && xid != UNDEFID && yid != UNDEFID;
  /*
  printf("swapxy %d\n", swapxy);
  printf("dimorder: %d %d %d\n", dimorder[0], dimorder[1], dimorder[2]);
  */

  size_t ndims = 0;
  if ( tsteptype != TSTEP_CONSTANT )
    {
      start[ndims] = (size_t)ntsteps - 1;
      count[ndims] = 1;
      ndims++;
    }

  for ( int id = 0; id < 3; ++id )
    {
      if ( dimorder[id] == 3 && zid != UNDEFID )
        {
          start[ndims] = (size_t)levelID;
          count[ndims] = 1;
          ndims++;
        }
      else if ( dimorder[id] == 2 && yid != UNDEFID )
        {
          start[ndims] = 0;
          cdf_inq_dimlen(fileID, yid, &ysize);
          count[ndims] = ysize;
          ndims++;
        }
      else if ( dimorder[id] == 1 && xid != UNDEFID )
        {
          start[ndims] = 0;
          cdf_inq_dimlen(fileID, xid, &xsize);
          count[ndims] = xsize;
          ndims++;
        }
    }

  if ( CDI_Debug )
    for (size_t idim = 0; idim < ndims; idim++)
      Message("dim = %d  start = %d  count = %d", idim, start[idim], count[idim]);

  int dtype = vlistInqVarDatatype(vlistID, varID);

  if ( nmiss > 0 ) cdfDefVarMissval(streamptr, varID, dtype, 1);

  long nvals = gridInqSize(gridID);

  cdf_write_var_data(fileID, vlistID, varID, ncvarid, dtype, nvals, xsize, ysize, swapxy, start, count, memtype, data, nmiss);

}


static
void cdfCreateRecords(stream_t *streamptr, int tsID)
{
  int varID, levelID, recID, vrecID, zaxisID;
  int nlev, nvrecs;

  int vlistID  = streamptr->vlistID;

  if ( tsID < 0 || (tsID >= streamptr->ntsteps && tsID > 0) ) return;

  if ( streamptr->tsteps[tsID].nallrecs > 0 ) return;

  tsteps_t* sourceTstep = streamptr->tsteps;
  tsteps_t* destTstep = sourceTstep + tsID;

  int nvars = vlistNvars(vlistID);
  int nrecs = vlistNrecs(vlistID);

  if ( nrecs <= 0 ) return;

  if ( tsID == 0 )
    {
      nvrecs = nrecs; /* use all records at first timestep */

      streamptr->nrecs += nrecs;

      destTstep->records    = (record_t *)xmalloc((size_t)nrecs*sizeof(record_t));
      destTstep->nrecs      = nrecs;
      destTstep->nallrecs   = nrecs;
      destTstep->recordSize = nrecs;
      destTstep->curRecID   = UNDEFID;
      destTstep->recIDs     = (int *)xmalloc((size_t)nvrecs*sizeof (int));;
      for ( recID = 0; recID < nvrecs; recID++ ) destTstep->recIDs[recID] = recID;

      record_t *records = destTstep->records;

      recID = 0;
      for ( varID = 0; varID < nvars; varID++ )
        {
          zaxisID = vlistInqVarZaxis(vlistID, varID);
          nlev    = zaxisInqSize(zaxisID);
          for ( levelID = 0; levelID < nlev; levelID++ )
            {
              recordInitEntry(&records[recID]);
              records[recID].varID   = (short)varID;
              records[recID].levelID = (short)levelID;
              recID++;
            }
        }
    }
  else if ( tsID == 1 )
    {
      nvrecs = 0;
      for ( varID = 0; varID < nvars; varID++ )
        {
          if ( vlistInqVarTsteptype(vlistID, varID) != TSTEP_CONSTANT )
            {
              zaxisID = vlistInqVarZaxis(vlistID, varID);
              nvrecs += zaxisInqSize(zaxisID);
            }
        }

      streamptr->nrecs += nvrecs;

      destTstep->records    = (record_t *) xmalloc((size_t)nrecs*sizeof(record_t));
      destTstep->nrecs      = nvrecs;
      destTstep->nallrecs   = nrecs;
      destTstep->recordSize = nrecs;
      destTstep->curRecID   = UNDEFID;

      memcpy(destTstep->records, sourceTstep->records, (size_t)nrecs*sizeof(record_t));

      if ( nvrecs )
        {
          destTstep->recIDs = (int *) xmalloc((size_t)nvrecs * sizeof (int));
          vrecID = 0;
          for ( recID = 0; recID < nrecs; recID++ )
            {
              varID = destTstep->records[recID].varID;
              if ( vlistInqVarTsteptype(vlistID, varID) != TSTEP_CONSTANT )
                {
                  destTstep->recIDs[vrecID++] = recID;
                }
            }
        }
    }
  else
    {
      if ( streamptr->tsteps[1].records == 0 ) cdfCreateRecords(streamptr, 1);

      nvrecs = streamptr->tsteps[1].nrecs;

      streamptr->nrecs += nvrecs;

      destTstep->records    = (record_t *) xmalloc((size_t)nrecs*sizeof(record_t));
      destTstep->nrecs      = nvrecs;
      destTstep->nallrecs   = nrecs;
      destTstep->recordSize = nrecs;
      destTstep->curRecID   = UNDEFID;

      memcpy(destTstep->records, sourceTstep->records, (size_t)nrecs*sizeof(record_t));

      destTstep->recIDs     = (int *) xmalloc((size_t)nvrecs * sizeof(int));

      memcpy(destTstep->recIDs, streamptr->tsteps[1].recIDs, (size_t)nvrecs*sizeof(int));
    }
}


static
int cdfTimeDimID(int fileID, int ndims, int nvars)
{
  int dimid = UNDEFID;
  int timedimid = UNDEFID;
  char dimname[80];
  char timeunits[CDI_MAX_NAME];
  char attname[CDI_MAX_NAME];
  char name[CDI_MAX_NAME];
  nc_type xtype;
  int nvdims, nvatts;
  int dimids[9];
  int varid, iatt;

  for ( dimid = 0; dimid < ndims; dimid++ )
    {
      cdf_inq_dimname(fileID, dimid, dimname);
      if ( memcmp(dimname, "time", 4) == 0 )
        {
          timedimid = dimid;
          break;
        }
    }

  if ( timedimid == UNDEFID )
    {
      for ( varid = 0; varid < nvars; varid++ )
        {
          cdf_inq_var(fileID, varid, name, &xtype, &nvdims, dimids, &nvatts);
          if ( nvdims == 1 )
            {
              for ( iatt = 0; iatt < nvatts; iatt++ )
                {
                  cdf_inq_attname(fileID, varid, iatt, attname);
                  if ( strncmp(attname, "units", 5) == 0 )
                    {
                      cdfGetAttText(fileID, varid, "units", sizeof(timeunits), timeunits);
                      strtolower(timeunits);

                      if ( isTimeUnits(timeunits) )
                        {
                          timedimid = dimids[0];
                          break;
                        }
                    }
                }
            }
        }
    }

  return (timedimid);
}

static
void init_ncdims(long ndims, ncdim_t *ncdims)
{
  long ncdimid;

  for ( ncdimid = 0; ncdimid < ndims; ncdimid++ )
    {
      ncdims[ncdimid].ncvarid      = UNDEFID;
      ncdims[ncdimid].dimtype      = UNDEFID;
      ncdims[ncdimid].len          = 0;
      ncdims[ncdimid].name[0]      = 0;
    }
}

static
void init_ncvars(long nvars, ncvar_t *ncvars)
{
  long ncvarid;

  for ( ncvarid = 0; ncvarid < nvars; ++ncvarid )
    {
      ncvars[ncvarid].ncid            = UNDEFID;
      ncvars[ncvarid].ignore          = FALSE;
      ncvars[ncvarid].isvar           = UNDEFID;
      ncvars[ncvarid].islon           = FALSE;
      ncvars[ncvarid].islat           = FALSE;
      ncvars[ncvarid].islev           = FALSE;
      ncvars[ncvarid].istime          = FALSE;
      ncvars[ncvarid].warn            = FALSE;
      ncvars[ncvarid].tsteptype       = TSTEP_CONSTANT;
      ncvars[ncvarid].param           = UNDEFID;
      ncvars[ncvarid].code            = UNDEFID;
      ncvars[ncvarid].tabnum          = 0;
      ncvars[ncvarid].calendar        = FALSE;
      ncvars[ncvarid].climatology     = FALSE;
      ncvars[ncvarid].bounds          = UNDEFID;
      ncvars[ncvarid].gridID          = UNDEFID;
      ncvars[ncvarid].zaxisID         = UNDEFID;
      ncvars[ncvarid].gridtype        = UNDEFID;
      ncvars[ncvarid].zaxistype       = UNDEFID;
      ncvars[ncvarid].xdim            = UNDEFID;
      ncvars[ncvarid].ydim            = UNDEFID;
      ncvars[ncvarid].zdim            = UNDEFID;
      ncvars[ncvarid].xvarid          = UNDEFID;
      ncvars[ncvarid].yvarid          = UNDEFID;
      ncvars[ncvarid].zvarid          = UNDEFID;
      ncvars[ncvarid].tvarid          = UNDEFID;
      ncvars[ncvarid].ncoordvars      = 0;
      for ( int i = 0; i < MAX_COORDVARS; ++i )
        ncvars[ncvarid].coordvarids[i]  = UNDEFID;
      ncvars[ncvarid].nauxvars      = 0;
      for ( int i = 0; i < MAX_AUXVARS; ++i )
        ncvars[ncvarid].auxvarids[i]  = UNDEFID;
      ncvars[ncvarid].cellarea        = UNDEFID;
      ncvars[ncvarid].tableID         = UNDEFID;
      ncvars[ncvarid].xtype           = 0;
      ncvars[ncvarid].ndims           = 0;
      ncvars[ncvarid].gmapid          = UNDEFID;
      ncvars[ncvarid].vlen            = 0;
      ncvars[ncvarid].vdata           = NULL;
      ncvars[ncvarid].truncation      = 0;
      ncvars[ncvarid].position        = 0;
      ncvars[ncvarid].positive        = 0;
      ncvars[ncvarid].chunked         = 0;
      ncvars[ncvarid].chunktype       = UNDEFID;
      ncvars[ncvarid].defmissval      = 0;
      ncvars[ncvarid].deffillval      = 0;
      ncvars[ncvarid].missval         = 0;
      ncvars[ncvarid].fillval         = 0;
      ncvars[ncvarid].addoffset       = 0;
      ncvars[ncvarid].scalefactor     = 1;
      ncvars[ncvarid].name[0]         = 0;
      ncvars[ncvarid].longname[0]     = 0;
      ncvars[ncvarid].stdname[0]      = 0;
      ncvars[ncvarid].units[0]        = 0;
      ncvars[ncvarid].extra[0]        = 0;
      ncvars[ncvarid].natts           = 0;
      ncvars[ncvarid].atts            = NULL;
      ncvars[ncvarid].deflate         = 0;
      ncvars[ncvarid].lunsigned       = 0;
      ncvars[ncvarid].lvalidrange     = 0;
      ncvars[ncvarid].validrange[0]   = VALIDMISS;
      ncvars[ncvarid].validrange[1]   = VALIDMISS;
      ncvars[ncvarid].ensdata         = NULL;
    }
}

static
int isLonAxis(const char *units, const char *stdname)
{
  int status = FALSE;
  char lc_units[16];

  memcpy(lc_units, units, 15);
  lc_units[15] = 0;
  strtolower(lc_units);

  if ( ((memcmp(lc_units, "degree", 6) == 0 || memcmp(lc_units, "radian", 6) == 0) &&
        (memcmp(stdname, "grid_longitude", 14) == 0 || memcmp(stdname, "longitude", 9) == 0)) )
    {
      status = TRUE;
    }

  if ( status == FALSE &&
       memcmp(stdname, "grid_latitude", 13) && memcmp(stdname, "latitude", 8) &&
       memcmp(lc_units, "degree", 6) == 0 )
    {
      int ioff = 6;
      if ( lc_units[ioff] == 's' ) ioff++;
      if ( lc_units[ioff] == '_' ) ioff++;
      if ( lc_units[ioff] == 'e' ) status = TRUE;
    }

  return (status);
}

static
int isLatAxis(const char *units, const char *stdname)
{
  int status = FALSE;
  char lc_units[16];

  memcpy(lc_units, units, 15);
  lc_units[15] = 0;
  strtolower(lc_units);

  if ( ((memcmp(lc_units, "degree", 6) == 0 || memcmp(lc_units, "radian", 6) == 0) &&
        (memcmp(stdname, "grid_latitude", 13) == 0 || memcmp(stdname, "latitude", 8) == 0)) )
    {
      status = TRUE;
    }

  if ( status == FALSE &&
       memcmp(stdname, "grid_longitude", 14) && memcmp(stdname, "longitude", 9) &&
       memcmp(lc_units, "degree", 6) == 0 )
    {
      int ioff = 6;
      if ( lc_units[ioff] == 's' ) ioff++;
      if ( lc_units[ioff] == '_' ) ioff++;
      if ( lc_units[ioff] == 'n' ) status = TRUE;
    }

  return (status);
}

static
int isDBLAxis(/*const char *units,*/ const char *longname)
{
  int status = FALSE;

  if ( strcmp(longname, "depth below land")         == 0 ||
       strcmp(longname, "depth_below_land")         == 0 ||
       strcmp(longname, "levels below the surface") == 0 )
    {
      /*
      if ( strcmp(ncvars[ncvarid].units, "cm") == 0 ||
           strcmp(ncvars[ncvarid].units, "dm") == 0 ||
           strcmp(ncvars[ncvarid].units, "m")  == 0 )
      */
        status = TRUE;
    }

  return (status);
}

static
int unitsIsMeter(const char *units)
{
  return (units[0] == 'm' && (!units[1] || strncmp(units, "meter", 5) == 0));
}

static
int isDepthAxis(const char *stdname, const char *longname)
{
  int status = FALSE;

  if ( strcmp(stdname, "depth") == 0 ) status = TRUE;

  if ( status == FALSE )
    if ( strcmp(longname, "depth_below_sea") == 0 ||
         strcmp(longname, "depth below sea") == 0 )
      {
        status = TRUE;
      }

  return (status);
}

static
int isHeightAxis(const char *stdname, const char *longname)
{
  int status = FALSE;

  if ( strcmp(stdname, "height") == 0 ) status = TRUE;

  if ( status == FALSE )
    if ( strcmp(longname, "height") == 0 ||
         strcmp(longname, "height above the surface") == 0 )
      {
        status = TRUE;
      }

  return (status);
}

static
int unitsIsPressure(const char *units)
{
  int status = FALSE;

  if ( memcmp(units, "millibar", 8) == 0 ||
       memcmp(units, "mb", 2)       == 0 ||
       memcmp(units, "hectopas", 8) == 0 ||
       memcmp(units, "hPa", 3)      == 0 ||
       memcmp(units, "Pa", 2)       == 0 )
    {
      status = TRUE;
    }

  return (status);
}

static
int isGaussGrid(size_t ysize, double yinc, double *yvals)
{
  int lgauss = FALSE;
  double *yv, *yw;

  if ( IS_EQUAL(yinc, 0) && ysize > 2 ) /* check if gaussian */
    {
      size_t i;
      yv = (double *) malloc(ysize*sizeof(double));
      yw = (double *) malloc(ysize*sizeof(double));
      gaussaw(yv, yw, ysize);
      free(yw);
      for ( i = 0; i < ysize; i++ )
        yv[i] = asin(yv[i])/M_PI*180.0;

      for ( i = 0; i < ysize; i++ )
        if ( fabs(yv[i] - yvals[i]) >
             ((yv[0] - yv[1])/500) ) break;

      if ( i == ysize ) lgauss = TRUE;

      /* check S->N */
      if ( lgauss == FALSE )
        {
          for ( i = 0; i < ysize; i++ )
            if ( fabs(yv[i] - yvals[ysize-i-1]) >
                 ((yv[0] - yv[1])/500) ) break;

          if ( i == ysize ) lgauss = TRUE;
        }

      free(yv);
    }

  return (lgauss);
}

static
void cdfSetVar(ncvar_t *ncvars, int ncvarid, int isvar)
{
  if ( isvar != TRUE && isvar != FALSE )
    Error("Internal problem! var %s undefined", ncvars[ncvarid].name);

  if ( ncvars[ncvarid].isvar != UNDEFID &&
       ncvars[ncvarid].isvar != isvar   &&
       ncvars[ncvarid].warn  == FALSE )
    {
      if ( ! ncvars[ncvarid].ignore )
        Warning("Inconsistent variable definition for %s!", ncvars[ncvarid].name);

      ncvars[ncvarid].warn = TRUE;
      isvar = FALSE;
    }

  ncvars[ncvarid].isvar = isvar;
}

static
void cdfSetDim(ncvar_t *ncvars, int ncvarid, int dimid, int dimtype)
{
  if ( ncvars[ncvarid].dimtype[dimid] != UNDEFID &&
       ncvars[ncvarid].dimtype[dimid] != dimtype )
    {
      Warning("Inconsistent dimension definition for %s! dimid = %d;  type = %d;  newtype = %d",
              ncvars[ncvarid].name, dimid, ncvars[ncvarid].dimtype[dimid], dimtype);
    }

  ncvars[ncvarid].dimtype[dimid] = dimtype;
}

static
void printNCvars(ncvar_t *ncvars, int nvars, const char *oname)
{
  char axis[7];
  int ncvarid, i;
  int ndim;
  static const char iaxis[] = {'t', 'z', 'y', 'x'};

  fprintf(stderr, "%s:\n", oname);

  for ( ncvarid = 0; ncvarid < nvars; ncvarid++ )
    {
      ndim = 0;
      if ( ncvars[ncvarid].isvar )
        {
          axis[ndim++] = 'v';
          axis[ndim++] = ':';
          for ( i = 0; i < ncvars[ncvarid].ndims; i++ )
            {/*
              if      ( ncvars[ncvarid].tvarid != -1 ) axis[ndim++] = iaxis[0];
              else if ( ncvars[ncvarid].zvarid != -1 ) axis[ndim++] = iaxis[1];
              else if ( ncvars[ncvarid].yvarid != -1 ) axis[ndim++] = iaxis[2];
              else if ( ncvars[ncvarid].xvarid != -1 ) axis[ndim++] = iaxis[3];
              else
             */
              if      ( ncvars[ncvarid].dimtype[i] == T_AXIS ) axis[ndim++] = iaxis[0];
              else if ( ncvars[ncvarid].dimtype[i] == Z_AXIS ) axis[ndim++] = iaxis[1];
              else if ( ncvars[ncvarid].dimtype[i] == Y_AXIS ) axis[ndim++] = iaxis[2];
              else if ( ncvars[ncvarid].dimtype[i] == X_AXIS ) axis[ndim++] = iaxis[3];
              else                                             axis[ndim++] = '?';
            }
        }
      else
        {
          axis[ndim++] = 'c';
          axis[ndim++] = ':';
          if      ( ncvars[ncvarid].istime ) axis[ndim++] = iaxis[0];
          else if ( ncvars[ncvarid].islev  ) axis[ndim++] = iaxis[1];
          else if ( ncvars[ncvarid].islat  ) axis[ndim++] = iaxis[2];
          else if ( ncvars[ncvarid].islon  ) axis[ndim++] = iaxis[3];
          else                               axis[ndim++] = '?';
        }

      axis[ndim++] = 0;

      fprintf(stderr, "%3d %3d  %-6s %s\n", ncvarid, ndim-3, axis, ncvars[ncvarid].name);
    }
}

static
void cdfScanVarAttributes(int nvars, ncvar_t *ncvars, ncdim_t *ncdims,
                          int timedimid, int modelID, int format)
{
  int ncid;
  int ncdimid;
  int nvdims, nvatts;
  int *dimidsp;
  int iatt;
  int i;
  int tablenum;
  nc_type xtype, atttype;
  size_t attlen;
  char name[CDI_MAX_NAME];
  char attname[CDI_MAX_NAME];
  const int attstringlen = 8192; char attstring[8192];

  int nchecked_vars = 0;
  char *checked_vars[9];
  for ( i = 0; i < 9; ++i ) checked_vars[i] = NULL;

  for ( int ncvarid = 0; ncvarid < nvars; ncvarid++ )
    {
      ncid    = ncvars[ncvarid].ncid;
      dimidsp = ncvars[ncvarid].dimids;

      cdf_inq_var(ncid, ncvarid, name, &xtype, &nvdims, dimidsp, &nvatts);
      strcpy(ncvars[ncvarid].name, name);

      for ( ncdimid = 0; ncdimid < nvdims; ncdimid++ )
        ncvars[ncvarid].dimtype[ncdimid] = -1;

      ncvars[ncvarid].xtype = xtype;
      ncvars[ncvarid].ndims = nvdims;

#if  defined  (HAVE_NETCDF4)
      if ( format == NC_FORMAT_NETCDF4_CLASSIC || format == NC_FORMAT_NETCDF4 )
        {
          char buf[CDI_MAX_NAME];
          int shuffle, deflate, deflate_level;
          size_t chunks[nvdims];
          int storage_in;
          nc_inq_var_deflate(ncid, ncvarid, &shuffle, &deflate, &deflate_level);
          if ( deflate > 0 ) ncvars[ncvarid].deflate = 1;

          if ( nc_inq_var_chunking(ncid, ncvarid, &storage_in, chunks) == NC_NOERR )
            {
              if ( storage_in == NC_CHUNKED )
                {
                  ncvars[ncvarid].chunked = 1;
                  for ( int i = 0; i < nvdims; ++i ) ncvars[ncvarid].chunks[i] = (int)chunks[i];
                  if ( CDI_Debug )
                    {
                      fprintf(stderr, "%s: chunking %d %d %d  chunks ", name, storage_in, NC_CONTIGUOUS, NC_CHUNKED);
                      for ( int i = 0; i < nvdims; ++i ) fprintf(stderr, "%ld ", chunks[i]);
                      fprintf(stderr, "\n");
                    }
                  strcat(ncvars[ncvarid].extra, "chunks=");
                  for ( int i = nvdims-1; i >= 0; --i )
                    {
                      sprintf(buf, "%ld", (long) chunks[i]);
                      strcat(ncvars[ncvarid].extra, buf);
                      if ( i > 0 ) strcat(ncvars[ncvarid].extra, "x");
                    }
                  strcat(ncvars[ncvarid].extra, " ");
                }
            }
        }
#endif

      if ( nvdims > 0 )
        {
          if ( timedimid == dimidsp[0] )
            {
              ncvars[ncvarid].tsteptype = TSTEP_INSTANT;
              cdfSetDim(ncvars, ncvarid, 0, T_AXIS);
            }
          else
            {
              for ( ncdimid = 1; ncdimid < nvdims; ncdimid++ )
                {
                  if ( timedimid == dimidsp[ncdimid] )
                    {
                      Warning("Time must be the first dimension! Unsupported array structure, skipped variable %s!", ncvars[ncvarid].name);
                      ncvars[ncvarid].isvar = FALSE;
                    }
                }
            }
        }

      for ( iatt = 0; iatt < nvatts; iatt++ )
        {
          cdf_inq_attname(ncid, ncvarid, iatt, attname);
          cdf_inq_atttype(ncid, ncvarid, attname, &atttype);
          cdf_inq_attlen(ncid, ncvarid, attname, &attlen);

          if ( strcmp(attname, "long_name") == 0 && xtypeIsText(atttype) )
            {
              cdfGetAttText(ncid, ncvarid, attname, CDI_MAX_NAME, ncvars[ncvarid].longname);
            }
          else if ( strcmp(attname, "standard_name") == 0 && xtypeIsText(atttype) )
            {
              cdfGetAttText(ncid, ncvarid, attname, CDI_MAX_NAME, ncvars[ncvarid].stdname);
            }
          else if ( strcmp(attname, "units") == 0 && xtypeIsText(atttype) )
            {
              cdfGetAttText(ncid, ncvarid, attname, CDI_MAX_NAME, ncvars[ncvarid].units);
            }
          else if ( strcmp(attname, "calendar") == 0 )
            {
              ncvars[ncvarid].calendar = TRUE;
            }
          else if ( strcmp(attname, "param") == 0 && xtypeIsText(atttype) )
            {
	      char paramstr[32];
	      int pnum = 0, pcat = 255, pdis = 255;
              cdfGetAttText(ncid, ncvarid, attname, sizeof(paramstr), paramstr);
	      sscanf(paramstr, "%d.%d.%d", &pnum, &pcat, &pdis);
	      ncvars[ncvarid].param = cdiEncodeParam(pnum, pcat, pdis);
              cdfSetVar(ncvars, ncvarid, TRUE);
            }
          else if ( strcmp(attname, "code") == 0 && !xtypeIsText(atttype) )
            {
              cdfGetAttInt(ncid, ncvarid, attname, 1, &ncvars[ncvarid].code);
              cdfSetVar(ncvars, ncvarid, TRUE);
            }
          else if ( strcmp(attname, "table") == 0 && !xtypeIsText(atttype) )
            {
              cdfGetAttInt(ncid, ncvarid, attname, 1, &tablenum);
              if ( tablenum > 0 )
                {
                  ncvars[ncvarid].tabnum = tablenum;
                  ncvars[ncvarid].tableID = tableInq(modelID, tablenum, NULL);
                  if ( ncvars[ncvarid].tableID == CDI_UNDEFID )
                    ncvars[ncvarid].tableID = tableDef(modelID, tablenum, NULL);
                }
              cdfSetVar(ncvars, ncvarid, TRUE);
            }
          else if ( strcmp(attname, "trunc_type") == 0 && xtypeIsText(atttype) )
            {
              cdfGetAttText(ncid, ncvarid, attname, attstringlen, attstring);
              if ( memcmp(attstring, "Triangular", attlen) == 0 )
                ncvars[ncvarid].gridtype = GRID_SPECTRAL;
            }
          else if ( strcmp(attname, "grid_type") == 0 && xtypeIsText(atttype) )
            {
              cdfGetAttText(ncid, ncvarid, attname, attstringlen, attstring);
              strtolower(attstring);

              if      ( strcmp(attstring, "gaussian reduced") == 0 )
                ncvars[ncvarid].gridtype = GRID_GAUSSIAN_REDUCED;
              else if ( strcmp(attstring, "gaussian") == 0 )
                ncvars[ncvarid].gridtype = GRID_GAUSSIAN;
              else if ( strncmp(attstring, "spectral", 8) == 0 )
                ncvars[ncvarid].gridtype = GRID_SPECTRAL;
              else if ( strncmp(attstring, "fourier", 7) == 0 )
                ncvars[ncvarid].gridtype = GRID_FOURIER;
              else if ( strcmp(attstring, "trajectory") == 0 )
                ncvars[ncvarid].gridtype = GRID_TRAJECTORY;
              else if ( strcmp(attstring, "generic") == 0 )
                ncvars[ncvarid].gridtype = GRID_GENERIC;
              else if ( strcmp(attstring, "cell") == 0 )
                ncvars[ncvarid].gridtype = GRID_UNSTRUCTURED;
              else if ( strcmp(attstring, "unstructured") == 0 )
                ncvars[ncvarid].gridtype = GRID_UNSTRUCTURED;
              else if ( strcmp(attstring, "curvilinear") == 0 )
                ncvars[ncvarid].gridtype = GRID_CURVILINEAR;
              else if ( strcmp(attstring, "sinusoidal") == 0 )
                ;
              else if ( strcmp(attstring, "laea") == 0 )
                ;
              else if ( strcmp(attstring, "lcc2") == 0 )
                ;
              else if ( strcmp(attstring, "linear") == 0 ) // ignore grid type linear
                ;
              else
                {
                  static int warn = TRUE;
                  if ( warn )
                    {
                      warn = FALSE;
                      Warning("netCDF attribute grid_type='%s' unsupported!", attstring);
                    }
                }

              cdfSetVar(ncvars, ncvarid, TRUE);
            }
          else if ( strcmp(attname, "level_type") == 0 && xtypeIsText(atttype) )
            {
              cdfGetAttText(ncid, ncvarid, attname, attstringlen, attstring);
              strtolower(attstring);

              if      ( strcmp(attstring, "toa") == 0 )
                ncvars[ncvarid].zaxistype = ZAXIS_TOA;
              else if ( strcmp(attstring, "cloudbase") == 0 )
                ncvars[ncvarid].zaxistype = ZAXIS_CLOUD_BASE;
              else if ( strcmp(attstring, "cloudtop") == 0 )
                ncvars[ncvarid].zaxistype = ZAXIS_CLOUD_TOP;
              else if ( strcmp(attstring, "isotherm0") == 0 )
                ncvars[ncvarid].zaxistype = ZAXIS_ISOTHERM_ZERO;
              else if ( strcmp(attstring, "seabottom") == 0 )
                ncvars[ncvarid].zaxistype = ZAXIS_SEA_BOTTOM;
              else if ( strcmp(attstring, "lakebottom") == 0 )
                ncvars[ncvarid].zaxistype = ZAXIS_LAKE_BOTTOM;
              else if ( strcmp(attstring, "sedimentbottom") == 0 )
                ncvars[ncvarid].zaxistype = ZAXIS_SEDIMENT_BOTTOM;
              else if ( strcmp(attstring, "sedimentbottomta") == 0 )
                ncvars[ncvarid].zaxistype = ZAXIS_SEDIMENT_BOTTOM_TA;
              else if ( strcmp(attstring, "sedimentbottomtw") == 0 )
                ncvars[ncvarid].zaxistype = ZAXIS_SEDIMENT_BOTTOM_TW;
              else if ( strcmp(attstring, "mixlayer") == 0 )
                ncvars[ncvarid].zaxistype = ZAXIS_MIX_LAYER;
              else if ( strcmp(attstring, "atmosphere") == 0 )
                ncvars[ncvarid].zaxistype = ZAXIS_ATMOSPHERE;
              else
                {
                  static int warn = TRUE;
                  if ( warn )
                    {
                      warn = FALSE;
                      Warning("netCDF attribute level_type='%s' unsupported!", attstring);
                    }
                }

              cdfSetVar(ncvars, ncvarid, TRUE);
            }
          else if ( strcmp(attname, "trunc_count") == 0 && !xtypeIsText(atttype) )
            {
              cdfGetAttInt(ncid, ncvarid, attname, 1, &ncvars[ncvarid].truncation);
            }
          else if ( strcmp(attname, "truncation") == 0 && !xtypeIsText(atttype) )
            {
              cdfGetAttInt(ncid, ncvarid, attname, 1, &ncvars[ncvarid].truncation);
            }
          else if ( strcmp(attname, "number_of_grid_in_reference") == 0 && !xtypeIsText(atttype) )
            {
              cdfGetAttInt(ncid, ncvarid, attname, 1, &ncvars[ncvarid].position);
            }
          else if ( strcmp(attname, "add_offset") == 0 && !xtypeIsText(atttype) )
            {
	      cdfGetAttDouble(ncid, ncvarid, attname, 1, &ncvars[ncvarid].addoffset);
	      /*
		if ( atttype != NC_BYTE && atttype != NC_SHORT && atttype != NC_INT )
		if ( ncvars[ncvarid].addoffset != 0 )
		Warning("attribute add_offset not supported for atttype %d", atttype);
	      */
	      /* (also used for lon/lat) cdfSetVar(ncvars, ncvarid, TRUE); */
            }
          else if ( strcmp(attname, "scale_factor") == 0 && !xtypeIsText(atttype) )
            {
	      cdfGetAttDouble(ncid, ncvarid, attname, 1, &ncvars[ncvarid].scalefactor);
	      /*
		if ( atttype != NC_BYTE && atttype != NC_SHORT && atttype != NC_INT )
		if ( ncvars[ncvarid].scalefactor != 1 )
		Warning("attribute scale_factor not supported for atttype %d", atttype);
	      */
	      /* (also used for lon/lat) cdfSetVar(ncvars, ncvarid, TRUE); */
            }
          else if ( strcmp(attname, "climatology") == 0 && xtypeIsText(atttype) )
            {
              int status, ncboundsid;

              cdfGetAttText(ncid, ncvarid, attname, attstringlen, attstring);

              status = nc_inq_varid(ncid, attstring, &ncboundsid);

              if ( status == NC_NOERR )
                {
                  ncvars[ncvarid].climatology = TRUE;
                  ncvars[ncvarid].bounds = ncboundsid;
                  cdfSetVar(ncvars, ncvars[ncvarid].bounds, FALSE);
                  cdfSetVar(ncvars, ncvarid, FALSE);
                }
              else
                Warning("%s - %s", nc_strerror(status), attstring);
            }
          else if ( strcmp(attname, "bounds") == 0 && xtypeIsText(atttype) )
            {
              int status, ncboundsid;

              cdfGetAttText(ncid, ncvarid, attname, attstringlen, attstring);

              status = nc_inq_varid(ncid, attstring, &ncboundsid);

              if ( status == NC_NOERR )
                {
                  ncvars[ncvarid].bounds = ncboundsid;
                  cdfSetVar(ncvars, ncvars[ncvarid].bounds, FALSE);
                  cdfSetVar(ncvars, ncvarid, FALSE);
                }
              else
                Warning("%s - %s", nc_strerror(status), attstring);
            }
          else if ( strcmp(attname, "cell_measures") == 0 && xtypeIsText(atttype) )
            {
              char *pstring, *cell_measures = NULL, *cell_var = NULL;

              cdfGetAttText(ncid, ncvarid, attname, attstringlen, attstring);
              pstring = attstring;

              while ( isspace((int) *pstring) ) pstring++;
              cell_measures = pstring;
              while ( isalnum((int) *pstring) ) pstring++;
              *pstring++ = 0;
              while ( isspace((int) *pstring) ) pstring++;
              cell_var = pstring;
              while ( ! isspace((int) *pstring) && *pstring != 0 ) pstring++;
              *pstring++ = 0;
              /*
              printf("cell_measures >%s<\n", cell_measures);
              printf("cell_var >%s<\n", cell_var);
              */
              if ( memcmp(cell_measures, "area", 4) == 0 )
                {
                  int status;
                  int nc_cell_id;

                  status = nc_inq_varid(ncid, cell_var, &nc_cell_id);
                  if ( status == NC_NOERR )
                    {
                      ncvars[ncvarid].cellarea = nc_cell_id;
                      /* ncvars[nc_cell_id].isvar = UNDEFID; */
                      cdfSetVar(ncvars, nc_cell_id, FALSE);
                    }
                  else
                    Warning("%s - %s", nc_strerror(status), cell_var);
                }
              else
                {
                  Warning("%s has an unexpected contents: %s", attname, cell_measures);
                }
              cdfSetVar(ncvars, ncvarid, TRUE);
            }
          /*
          else if ( strcmp(attname, "coordinates") == 0 )
            {
              char *pstring, *xvarname = NULL, *yvarname = NULL;

              cdfGetAttText(ncid, ncvarid, attname, attstringlen, attstring);
              pstring = attstring;

              while ( isspace((int) *pstring) ) pstring++;
              xvarname = pstring;
              while ( isgraph((int) *pstring) ) pstring++;
              *pstring++ = 0;
              while ( isspace((int) *pstring) ) pstring++;
              yvarname = pstring;
              while ( isgraph((int) *pstring) ) pstring++;
              *pstring++ = 0;

              cdf_inq_varid(ncid, xvarname, &ncvars[ncvarid].xvarid);
              cdf_inq_varid(ncid, yvarname, &ncvars[ncvarid].yvarid);

              cdfSetVar(ncvars, ncvars[ncvarid].xvarid, FALSE);
              cdfSetVar(ncvars, ncvars[ncvarid].yvarid, FALSE);
              cdfSetVar(ncvars, ncvarid, TRUE);
            }
          */
          else if ( (strcmp(attname, "associate")  == 0 || strcmp(attname, "coordinates") == 0) && xtypeIsText(atttype) )
            {
              int status;
              char *pstring, *varname = NULL;
              int lstop = FALSE;
              int dimvarid;
              extern int cdiIgnoreAttCoordinates;

              cdfGetAttText(ncid, ncvarid, attname, attstringlen, attstring);
              pstring = attstring;

              for ( i = 0; i < MAX_COORDVARS; i++ )
                {
                  while ( isspace((int) *pstring) ) pstring++;
                  if ( *pstring == 0 ) break;
                  varname = pstring;
                  while ( !isspace((int) *pstring) && *pstring != 0 ) pstring++;
                  if ( *pstring == 0 ) lstop = TRUE;
                  *pstring++ = 0;

                  status = nc_inq_varid(ncid, varname, &dimvarid);
                  if ( status == NC_NOERR )
                    {
                      cdfSetVar(ncvars, dimvarid, FALSE);
                      if ( cdiIgnoreAttCoordinates == FALSE )
                        {
                          ncvars[ncvarid].coordvarids[i] = dimvarid;
                          ncvars[ncvarid].ncoordvars++;
                        }
                    }
                  else
                    {
                      int k;
                      for ( k = 0; k < nchecked_vars; ++k )
                        if ( strcmp(checked_vars[k], varname) == 0 ) break;

                      if ( k == nchecked_vars )
                        {
                          if ( nchecked_vars < 9 ) checked_vars[nchecked_vars++] = strdup(varname);
                          Warning("%s - %s", nc_strerror(status), varname);
                        }
                    }

                  if ( lstop ) break;
                }

              cdfSetVar(ncvars, ncvarid, TRUE);
            }
          else if ( (strcmp(attname, "auxiliary_variable") == 0) && xtypeIsText(atttype) )
            {
              int status;
              char *pstring, *varname = NULL;
              int lstop = FALSE;
              int dimvarid;
              extern int cdiIgnoreAttCoordinates;

              cdfGetAttText(ncid, ncvarid, attname, attstringlen, attstring);
              pstring = attstring;

              for ( i = 0; i < MAX_AUXVARS; i++ )
                {
                  while ( isspace((int) *pstring) ) pstring++;
                  if ( *pstring == 0 ) break;
                  varname = pstring;
                  while ( !isspace((int) *pstring) && *pstring != 0 ) pstring++;
                  if ( *pstring == 0 ) lstop = TRUE;
                  *pstring++ = 0;

                  status = nc_inq_varid(ncid, varname, &dimvarid);
                  if ( status == NC_NOERR )
                    {
                      cdfSetVar(ncvars, dimvarid, FALSE);
                      //  if ( cdiIgnoreAttCoordinates == FALSE )
                        {
                          ncvars[ncvarid].auxvarids[i] = dimvarid;
                          ncvars[ncvarid].nauxvars++;
                        }
                    }
                  else
                    Warning("%s - %s", nc_strerror(status), varname);

                  if ( lstop ) break;
                }

              cdfSetVar(ncvars, ncvarid, TRUE);
            }
          else if ( strcmp(attname, "grid_mapping") == 0 && xtypeIsText(atttype) )
            {
              int status;
              int nc_gmap_id;

              cdfGetAttText(ncid, ncvarid, attname, attstringlen, attstring);

              status = nc_inq_varid(ncid, attstring, &nc_gmap_id);
              if ( status == NC_NOERR )
                {
                  ncvars[ncvarid].gmapid = nc_gmap_id;
                  cdfSetVar(ncvars, ncvars[ncvarid].gmapid, FALSE);
                }
              else
                Warning("%s - %s", nc_strerror(status), attstring);

              cdfSetVar(ncvars, ncvarid, TRUE);
            }
          else if ( strcmp(attname, "positive") == 0 && xtypeIsText(atttype) )
            {
              cdfGetAttText(ncid, ncvarid, attname, attstringlen, attstring);
              strtolower(attstring);

              if    ( memcmp(attstring, "down", 4) == 0 ) ncvars[ncvarid].positive = POSITIVE_DOWN;
              else if ( memcmp(attstring, "up", 2) == 0 ) ncvars[ncvarid].positive = POSITIVE_UP;

              if ( ncvars[ncvarid].ndims == 1 )
                {
                  cdfSetVar(ncvars, ncvarid, FALSE);
                  cdfSetDim(ncvars, ncvarid, 0, Z_AXIS);
                  ncdims[ncvars[ncvarid].dimids[0]].dimtype = Z_AXIS;
                }
            }
          else if ( strcmp(attname, "_FillValue") == 0 && !xtypeIsText(atttype) )
            {
	      cdfGetAttDouble(ncid, ncvarid, attname, 1, &ncvars[ncvarid].fillval);
	      ncvars[ncvarid].deffillval = TRUE;
	      /* cdfSetVar(ncvars, ncvarid, TRUE); */
            }
          else if ( strcmp(attname, "missing_value") == 0 && !xtypeIsText(atttype) )
            {
	      cdfGetAttDouble(ncid, ncvarid, attname, 1, &ncvars[ncvarid].missval);
	      ncvars[ncvarid].defmissval = TRUE;
	      /* cdfSetVar(ncvars, ncvarid, TRUE); */
            }
          else if ( strcmp(attname, "valid_range") == 0 && attlen == 2 )
            {
              if ( ncvars[ncvarid].lvalidrange == FALSE )
                {
                  extern int cdiIgnoreValidRange;
                  int lignore = FALSE;
                  if ( xtypeIsFloat(atttype) != xtypeIsFloat(xtype) ) lignore = TRUE;
                  if ( cdiIgnoreValidRange == FALSE && lignore == FALSE )
                    {
                      cdfGetAttDouble(ncid, ncvarid, attname, 2, ncvars[ncvarid].validrange);
                      ncvars[ncvarid].lvalidrange = TRUE;
                      if ( ((int)ncvars[ncvarid].validrange[0]) == 0 && ((int)ncvars[ncvarid].validrange[1]) == 255 )
                        ncvars[ncvarid].lunsigned = TRUE;
                      /* cdfSetVar(ncvars, ncvarid, TRUE); */
                    }
                  else if ( lignore )
                    {
                      Warning("Inconsistent data type for attribute %s:valid_range, ignored!", name);
                    }
                }
            }
          else if ( strcmp(attname, "valid_min") == 0 && attlen == 1 )
            {
              if ( ncvars[ncvarid].lvalidrange == FALSE )
                {
                  extern int cdiIgnoreValidRange;
                  int lignore = FALSE;
                  if ( xtypeIsFloat(atttype) != xtypeIsFloat(xtype) ) lignore = TRUE;
                  if ( cdiIgnoreValidRange == FALSE && lignore == FALSE )
                    {
                      cdfGetAttDouble(ncid, ncvarid, attname, 1, &(ncvars[ncvarid].validrange)[0]);
                      ncvars[ncvarid].lvalidrange = TRUE;
                    }
                  else if ( lignore )
                    {
                      Warning("Inconsistent data type for attribute %s:valid_min, ignored!", name);
                    }
                }
            }
          else if ( strcmp(attname, "valid_max") == 0 && attlen == 1 )
            {
              if ( ncvars[ncvarid].lvalidrange == FALSE )
                {
                  extern int cdiIgnoreValidRange;
                  int lignore = FALSE;
                  if ( xtypeIsFloat(atttype) != xtypeIsFloat(xtype) ) lignore = TRUE;
                  if ( cdiIgnoreValidRange == FALSE && lignore == FALSE )
                    {
                      cdfGetAttDouble(ncid, ncvarid, attname, 1, &(ncvars[ncvarid].validrange)[1]);
                      ncvars[ncvarid].lvalidrange = TRUE;
                    }
                  else if ( lignore )
                    {
                      Warning("Inconsistent data type for attribute %s:valid_max, ignored!", name);
                    }
                }
            }
          else if ( strcmp(attname, "_Unsigned") == 0 && xtypeIsText(atttype) )
            {
              cdfGetAttText(ncid, ncvarid, attname, attstringlen, attstring);
              strtolower(attstring);

              if ( memcmp(attstring, "true", 4) == 0 )
                {
                  ncvars[ncvarid].lunsigned = TRUE;
                  /*
                  ncvars[ncvarid].lvalidrange = TRUE;
                  ncvars[ncvarid].validrange[0] = 0;
                  ncvars[ncvarid].validrange[1] = 255;
                  */
                }
	      /* cdfSetVar(ncvars, ncvarid, TRUE); */
            }
          else if ( strcmp(attname, "cdi") == 0 && xtypeIsText(atttype) )
            {
	      cdfGetAttText(ncid, ncvarid, attname, attstringlen, attstring);
	      strtolower(attstring);

	      if ( memcmp(attstring, "ignore", 6) == 0 )
		{
		  ncvars[ncvarid].ignore = TRUE;
		  cdfSetVar(ncvars, ncvarid, FALSE);
		}
            }
          else if ( strcmp(attname, "axis") == 0 && xtypeIsText(atttype) )
            {
              cdfGetAttText(ncid, ncvarid, attname, attstringlen, attstring);
	      attlen = strlen(attstring);

	      if ( (int) attlen > nvdims )
		{
		  if ( nvdims > 0 )
		    Warning("Unexpected axis attribute length for %s, ignored!", name);
		}
	      else
		{
		  strtolower(attstring);
		  for ( i = 0; i < (int)attlen; ++i )
		    {
		      if ( attstring[i] != '-' && attstring[i] != 't' && attstring[i] != 'z' &&
			   attstring[i] != 'y' && attstring[i] != 'x' )
			{
			  Warning("Unexpected character in axis attribute for %s, ignored!", name);
			  break;
			}
		    }

		  if ( i == (int) attlen && (int) attlen == nvdims)
		    {
		      while ( attlen-- )
			{
			  if ( (int) attstring[attlen] == 't' )
			    {
			      if ( attlen != 0 ) Warning("axis attribute 't' not on first position");
			      cdfSetDim(ncvars, ncvarid, (int)attlen, T_AXIS);
			    }
			  else if ( (int) attstring[attlen] == 'z' )
			    {
			      ncvars[ncvarid].zdim = dimidsp[attlen];
			      cdfSetDim(ncvars, ncvarid, (int)attlen, Z_AXIS);

			      if ( ncvars[ncvarid].ndims == 1 )
				{
				  cdfSetVar(ncvars, ncvarid, FALSE);
				  ncdims[ncvars[ncvarid].dimids[0]].dimtype = Z_AXIS;
				}
			    }
			  else if ( (int) attstring[attlen] == 'y' )
			    {
			      ncvars[ncvarid].ydim = dimidsp[attlen];
			      cdfSetDim(ncvars, ncvarid, (int)attlen, Y_AXIS);

			      if ( ncvars[ncvarid].ndims == 1 )
				{
				  cdfSetVar(ncvars, ncvarid, FALSE);
				  ncdims[ncvars[ncvarid].dimids[0]].dimtype = Y_AXIS;
				}
			    }
			  else if ( (int) attstring[attlen] == 'x' )
			    {
			      ncvars[ncvarid].xdim = dimidsp[attlen];
			      cdfSetDim(ncvars, ncvarid, (int)attlen, X_AXIS);

			      if ( ncvars[ncvarid].ndims == 1 )
				{
				  cdfSetVar(ncvars, ncvarid, FALSE);
				  ncdims[ncvars[ncvarid].dimids[0]].dimtype = X_AXIS;
				}
			    }
			}
		    }
		}
	    }
	  else if ( ( strcmp(attname, "realization") == 0 )         ||
	            ( strcmp(attname, "ensemble_members") == 0 )    ||
	            ( strcmp(attname, "forecast_init_type") == 0 )    )
	    {
	      int temp;

	      if( ncvars[ncvarid].ensdata == NULL )
		ncvars[ncvarid].ensdata = (ensinfo_t *) malloc( sizeof( ensinfo_t ) );

	      cdfGetAttInt(ncid, ncvarid, attname, 1, &temp);

	      if( strcmp(attname, "realization") == 0 )
		ncvars[ncvarid].ensdata->ens_index = temp;
	      else if( strcmp(attname, "ensemble_members") == 0 )
		ncvars[ncvarid].ensdata->ens_count = temp;
	      else if( strcmp(attname, "forecast_init_type") == 0 )
		ncvars[ncvarid].ensdata->forecast_init_type = temp;

	      cdfSetVar(ncvars, ncvarid, TRUE);
	    }
	  else
	    {
	      if ( ncvars[ncvarid].natts == 0 )
		ncvars[ncvarid].atts
                  = (int *)xmalloc((size_t)nvatts * sizeof (int));

	      ncvars[ncvarid].atts[ncvars[ncvarid].natts++] = iatt;
	      /*
	      int attrint;
	      double attrflt;
	      nc_type attrtype;
	      cdf_inq_attlen(ncid, ncvarid, attname, &attlen);
	      cdf_inq_atttype(ncid, ncvarid, attname, &attrtype);
	      if ( attlen == 1 && (attrtype == NC_INT || attrtype == NC_SHORT) )
		{
		  cdfGetAttInt(ncid, ncvarid, attname, 1, &attrint);
		  printf("int: %s.%s = %d\n", ncvars[ncvarid].name, attname, attrint);
		}
	      else if ( attlen == 1 && (attrtype == NC_FLOAT || attrtype == NC_DOUBLE) )
		{
		  cdfGetAttDouble(ncid, ncvarid, attname, 1, &attrflt);
		  printf("flt: %s.%s = %g\n", ncvars[ncvarid].name, attname, attrflt);
		}
	      else if ( attrtype == NC_CHAR )
		{
		  cdfGetAttText(ncid, ncvarid, attname, attstringlen, attstring);
		  attstring[attlen] = 0;
		  printf("txt: %s.%s = %s\n", ncvars[ncvarid].name, attname, attstring);
		}
	      else
		printf("att: %s.%s = unknown\n", ncvars[ncvarid].name, attname);
	      */
	    }
	}
    }

  for ( i = 0; i < 9; ++i ) if ( checked_vars[i] ) free(checked_vars[i]);
}

static
void setDimType(int nvars, ncvar_t *ncvars, ncdim_t *ncdims)
{
  int ndims;
  int ncvarid, ncdimid;
  int i;

  for ( ncvarid = 0; ncvarid < nvars; ncvarid++ )
    {
      if ( ncvars[ncvarid].isvar == TRUE )
	{
	  int lxdim = 0, lydim = 0, lzdim = 0/* , ltdim = 0 */;
	  ndims = ncvars[ncvarid].ndims;
	  for ( i = 0; i < ndims; i++ )
	    {
	      ncdimid = ncvars[ncvarid].dimids[i];
	      if      ( ncdims[ncdimid].dimtype == X_AXIS ) cdfSetDim(ncvars, ncvarid, i, X_AXIS);
	      else if ( ncdims[ncdimid].dimtype == Y_AXIS ) cdfSetDim(ncvars, ncvarid, i, Y_AXIS);
	      else if ( ncdims[ncdimid].dimtype == Z_AXIS ) cdfSetDim(ncvars, ncvarid, i, Z_AXIS);
	      else if ( ncdims[ncdimid].dimtype == T_AXIS ) cdfSetDim(ncvars, ncvarid, i, T_AXIS);
	    }

	  if ( CDI_Debug )
	    {
	      Message("var %d %s", ncvarid, ncvars[ncvarid].name);
	      for ( i = 0; i < ndims; i++ )
		printf("  dim%d type=%d  ", i, ncvars[ncvarid].dimtype[i]);
	      printf("\n");
	    }

	  for ( i = 0; i < ndims; i++ )
	    {
	      if      ( ncvars[ncvarid].dimtype[i] == X_AXIS ) lxdim = TRUE;
	      else if ( ncvars[ncvarid].dimtype[i] == Y_AXIS ) lydim = TRUE;
	      else if ( ncvars[ncvarid].dimtype[i] == Z_AXIS ) lzdim = TRUE;
	      /* else if ( ncvars[ncvarid].dimtype[i] == T_AXIS ) ltdim = TRUE; */
	    }

          if ( lxdim == FALSE && ncvars[ncvarid].xvarid != UNDEFID )
            {
              if (  ncvars[ncvars[ncvarid].xvarid].ndims == 0 ) lxdim = TRUE;
            }

          if ( lydim == FALSE && ncvars[ncvarid].yvarid != UNDEFID )
            {
              if (  ncvars[ncvars[ncvarid].yvarid].ndims == 0 ) lydim = TRUE;
            }

          //   if ( ndims > 1 )
            for ( i = ndims-1; i >= 0; i-- )
              {
                if ( ncvars[ncvarid].dimtype[i] == -1 )
                  {
                    if ( lxdim == FALSE )
                      {
                        cdfSetDim(ncvars, ncvarid, i, X_AXIS);
                        lxdim = TRUE;
                      }
                    else if ( lydim == FALSE && ncvars[ncvarid].gridtype != GRID_UNSTRUCTURED )
                      {
                        cdfSetDim(ncvars, ncvarid, i, Y_AXIS);
                        lydim = TRUE;
                      }
                    else if ( lzdim == FALSE )
                      {
                        cdfSetDim(ncvars, ncvarid, i, Z_AXIS);
                        lzdim = TRUE;
                      }
                  }
              }
	}
    }
}

/* verify coordinate vars - first scan (dimname == varname) */
static
void verify_coordinate_vars_1(int ndims, ncdim_t *ncdims, ncvar_t *ncvars, int timedimid)
{
  int ncdimid, ncvarid;

  for ( ncdimid = 0; ncdimid < ndims; ncdimid++ )
    {
      ncvarid = ncdims[ncdimid].ncvarid;
      if ( ncvarid != -1 )
	{
	  if ( ncvars[ncvarid].dimids[0] == timedimid )
	    {
              ncvars[ncvarid].istime = TRUE;
	      ncdims[ncdimid].dimtype = T_AXIS;
	      continue;
	    }

	  if ( ncvars[ncvarid].units[0] != 0 )
	    {
	      if ( isLonAxis(ncvars[ncvarid].units, ncvars[ncvarid].stdname) )
		{
		  ncvars[ncvarid].islon = TRUE;
		  cdfSetVar(ncvars, ncvarid, FALSE);
		  cdfSetDim(ncvars, ncvarid, 0, X_AXIS);
		  ncdims[ncdimid].dimtype = X_AXIS;
		}
	      else if ( isLatAxis(ncvars[ncvarid].units, ncvars[ncvarid].stdname) )
		{
		  ncvars[ncvarid].islat = TRUE;
		  cdfSetVar(ncvars, ncvarid, FALSE);
		  cdfSetDim(ncvars, ncvarid, 0, Y_AXIS);
		  ncdims[ncdimid].dimtype = Y_AXIS;
		}
	      else if ( unitsIsPressure(ncvars[ncvarid].units) )
		{
		  ncvars[ncvarid].zaxistype = ZAXIS_PRESSURE;
		}
	      else if ( strcmp(ncvars[ncvarid].units, "level") == 0 || strcmp(ncvars[ncvarid].units, "1") == 0 )
		{
		  if      ( strcmp(ncvars[ncvarid].longname, "hybrid level at layer midpoints") == 0 )
		    ncvars[ncvarid].zaxistype = ZAXIS_HYBRID;
		  else if ( strncmp(ncvars[ncvarid].longname, "hybrid level at midpoints", 25) == 0 )
		    ncvars[ncvarid].zaxistype = ZAXIS_HYBRID;
		  else if ( strcmp(ncvars[ncvarid].longname, "hybrid level at layer interfaces") == 0 )
		    ncvars[ncvarid].zaxistype = ZAXIS_HYBRID_HALF;
		  else if ( strncmp(ncvars[ncvarid].longname, "hybrid level at interfaces", 26) == 0 )
		    ncvars[ncvarid].zaxistype = ZAXIS_HYBRID_HALF;
		  else if ( strcmp(ncvars[ncvarid].units, "level") == 0 )
		    ncvars[ncvarid].zaxistype = ZAXIS_GENERIC;
		}
	      else if ( isDBLAxis(ncvars[ncvarid].longname) )
                {
                  ncvars[ncvarid].zaxistype = ZAXIS_DEPTH_BELOW_LAND;
		}
	      else if ( unitsIsMeter(ncvars[ncvarid].units) )
		{
		  if ( isDepthAxis(ncvars[ncvarid].stdname, ncvars[ncvarid].longname) )
		    ncvars[ncvarid].zaxistype = ZAXIS_DEPTH_BELOW_SEA;
		  else if ( isHeightAxis(ncvars[ncvarid].stdname, ncvars[ncvarid].longname) )
		    ncvars[ncvarid].zaxistype = ZAXIS_HEIGHT;
		}
	    }

	  if ( ncvars[ncvarid].islon == FALSE && ncvars[ncvarid].longname[0] != 0 &&
               ncvars[ncvarid].islat == FALSE && ncvars[ncvarid].longname[1] != 0 )
	    {
	      if ( memcmp(ncvars[ncvarid].longname+1, "ongitude", 8) == 0 )
		{
		  ncvars[ncvarid].islon = TRUE;
		  cdfSetVar(ncvars, ncvarid, FALSE);
		  cdfSetDim(ncvars, ncvarid, 0, X_AXIS);
		  ncdims[ncdimid].dimtype = X_AXIS;
		  continue;
		}
	      else if ( memcmp(ncvars[ncvarid].longname+1, "atitude", 7) == 0 )
		{
		  ncvars[ncvarid].islat = TRUE;
		  cdfSetVar(ncvars, ncvarid, FALSE);
		  cdfSetDim(ncvars, ncvarid, 0, Y_AXIS);
		  ncdims[ncdimid].dimtype = Y_AXIS;
		  continue;
		}
	    }

	  if ( ncvars[ncvarid].zaxistype != UNDEFID )
	    {
              ncvars[ncvarid].islev = TRUE;
	      cdfSetVar(ncvars, ncvarid, FALSE);
	      cdfSetDim(ncvars, ncvarid, 0, Z_AXIS);
	      ncdims[ncdimid].dimtype = Z_AXIS;
	    }
	}
    }
}

/* verify coordinate vars - second scan (all other variables) */
static
void verify_coordinate_vars_2(int nvars, ncvar_t *ncvars)
{
  int ncvarid;

  for ( ncvarid = 0; ncvarid < nvars; ncvarid++ )
    {
      if ( ncvars[ncvarid].isvar == 0 )
	{
	  if ( ncvars[ncvarid].units[0] != 0 )
	    {
	      if ( isLonAxis(ncvars[ncvarid].units, ncvars[ncvarid].stdname) )
		{
		  ncvars[ncvarid].islon = TRUE;
		  continue;
		}
	      else if ( isLatAxis(ncvars[ncvarid].units, ncvars[ncvarid].stdname) )
		{
		  ncvars[ncvarid].islat = TRUE;
		  continue;
		}
	      else if ( unitsIsPressure(ncvars[ncvarid].units) )
		{
		  ncvars[ncvarid].zaxistype = ZAXIS_PRESSURE;
		  continue;
		}
	      else if ( strcmp(ncvars[ncvarid].units, "level") == 0 || strcmp(ncvars[ncvarid].units, "1") == 0 )
		{
		  if      ( strcmp(ncvars[ncvarid].longname, "hybrid level at layer midpoints") == 0 )
		    ncvars[ncvarid].zaxistype = ZAXIS_HYBRID;
		  else if ( strncmp(ncvars[ncvarid].longname, "hybrid level at midpoints", 25) == 0 )
		    ncvars[ncvarid].zaxistype = ZAXIS_HYBRID;
		  else if ( strcmp(ncvars[ncvarid].longname, "hybrid level at layer interfaces") == 0 )
		    ncvars[ncvarid].zaxistype = ZAXIS_HYBRID_HALF;
		  else if ( strncmp(ncvars[ncvarid].longname, "hybrid level at interfaces", 26) == 0 )
		    ncvars[ncvarid].zaxistype = ZAXIS_HYBRID_HALF;
		  else if ( strcmp(ncvars[ncvarid].units, "level") == 0 )
		    ncvars[ncvarid].zaxistype = ZAXIS_GENERIC;
		  continue;
		}
	      else if ( isDBLAxis(ncvars[ncvarid].longname) )
		{
                  ncvars[ncvarid].zaxistype = ZAXIS_DEPTH_BELOW_LAND;
		  continue;
		}
	      else if ( unitsIsMeter(ncvars[ncvarid].units) )
		{
		  if ( isDepthAxis(ncvars[ncvarid].stdname, ncvars[ncvarid].longname) )
		    ncvars[ncvarid].zaxistype = ZAXIS_DEPTH_BELOW_SEA;
		  else if ( isHeightAxis(ncvars[ncvarid].stdname, ncvars[ncvarid].longname) )
		    ncvars[ncvarid].zaxistype = ZAXIS_HEIGHT;
		  continue;
		}
            }

	  /* not needed anymore for rotated grids */
	  if ( ncvars[ncvarid].islon == FALSE && ncvars[ncvarid].longname[0] != 0 &&
               ncvars[ncvarid].islat == FALSE && ncvars[ncvarid].longname[1] != 0 )
	    {
	      if ( memcmp(ncvars[ncvarid].longname+1, "ongitude", 8) == 0 )
		{
		  ncvars[ncvarid].islon = TRUE;
		  continue;
		}
	      else if ( memcmp(ncvars[ncvarid].longname+1, "atitude", 7) == 0 )
		{
		  ncvars[ncvarid].islat = TRUE;
		  continue;
		}
	    }
	}
    }
}

#if defined (PROJECTION_TEST)
static
void copy_numeric_projatts(int gridID, int ncvarID, int ncfileID)
{
  int iatt, nvatts;
  size_t attlen;
  char attname[CDI_MAX_NAME];
  nc_type xtype;

  cdf_inq_varnatts(ncfileID, ncvarID, &nvatts);

  for ( iatt = 0; iatt < nvatts; iatt++ )
    {
      cdf_inq_attname(ncfileID, ncvarID, iatt, attname);
      cdf_inq_atttype(ncfileID, ncvarID, attname, &xtype);
      cdf_inq_attlen(ncfileID, ncvarID, attname, &attlen);

      //  printf("%s %d\n", attname, (int)attlen);
    }

}
#endif

static
void grid_set_chunktype(grid_t *grid, ncvar_t *ncvar)
{
  if ( ncvar->chunked )
    {
      int ndims = ncvar->ndims;

      if ( grid->type == GRID_UNSTRUCTURED )
        {
          if ( ncvar->chunks[ndims-1] == grid->size )
            ncvar->chunktype = CHUNK_GRID;
          else
            ncvar->chunktype = CHUNK_AUTO;
        }
      else
        {
          if ( grid->xsize > 1 && grid->ysize > 1 && ndims > 1 &&
               grid->xsize == ncvar->chunks[ndims-1] &&
               grid->ysize == ncvar->chunks[ndims-2] )
            ncvar->chunktype = CHUNK_GRID;
          else if ( grid->xsize > 1 && grid->xsize == ncvar->chunks[ndims-1] )
            ncvar->chunktype = CHUNK_LINES;
          else
            ncvar->chunktype = CHUNK_AUTO;
        }
    }
}

/* define all input grids */
static
void define_all_grids(stream_t *streamptr, int vlistID, ncdim_t *ncdims, int nvars, ncvar_t *ncvars, int timedimid, unsigned char *uuidOfHGrid, char *gridfile, int number_of_grid_used)
{
  int ncvarid, ncvarid2;
  int ndims;
  int nbdims;
  int i;
  int nvatts;
  int skipvar;
  size_t nvertex;
  grid_t grid;
  grid_t proj;
  int gridindex;
  size_t size = 0, xsize, ysize, np;
  char name[CDI_MAX_NAME];
  int iatt;
  int ltwarn = TRUE;
  size_t attlen;
  char attname[CDI_MAX_NAME];
  double datt;

  for ( ncvarid = 0; ncvarid < nvars; ++ncvarid )
    {
      if ( ncvars[ncvarid].isvar && ncvars[ncvarid].gridID == UNDEFID )
	{
	  int xdimids[2] = {-1,-1}, ydimids[2] = {-1,-1};
	  int xdimid = -1, ydimid = -1;
	  int xvarid = -1, yvarid = -1;
	  int islon = 0, islat = 0;
	  int nxdims = 0, nydims = 0;
	  double xinc = 0, yinc = 0;

	  xsize = 0;
	  ysize = 0;
          np    = 0;

	  ndims = ncvars[ncvarid].ndims;
	  for ( i = 0; i < ndims; i++ )
	    {
	      if ( ncvars[ncvarid].dimtype[i] == X_AXIS && nxdims < 2 )
		{
		  xdimids[nxdims] = ncvars[ncvarid].dimids[i];
		  nxdims++;
		}
	      else if ( ncvars[ncvarid].dimtype[i] == Y_AXIS && nydims < 2 )
		{
		  ydimids[nydims] = ncvars[ncvarid].dimids[i];
		  nydims++;
		}
	    }

	  if ( nxdims == 2 )
	    {
	      xdimid = xdimids[1];
	      ydimid = xdimids[0];
	    }
	  else if ( nydims == 2 )
	    {
	      xdimid = ydimids[1];
	      ydimid = ydimids[0];
	    }
	  else
	    {
	      xdimid = xdimids[0];
	      ydimid = ydimids[0];
	    }

	  if ( ncvars[ncvarid].xvarid != UNDEFID )
	    xvarid = ncvars[ncvarid].xvarid;
	  else if ( xdimid != UNDEFID )
	    xvarid = ncdims[xdimid].ncvarid;

	  if ( ncvars[ncvarid].yvarid != UNDEFID )
	    yvarid = ncvars[ncvarid].yvarid;
	  else if ( ydimid != UNDEFID )
	    yvarid = ncdims[ydimid].ncvarid;

	  /*
	  if ( xdimid != UNDEFID )
	    xvarid = ncdims[xdimid].ncvarid;
	  if ( xvarid == UNDEFID && ncvars[ncvarid].xvarid != UNDEFID )
	    xvarid = ncvars[ncvarid].xvarid;

	  if ( ydimid != UNDEFID )
	    yvarid = ncdims[ydimid].ncvarid;
	  if ( yvarid == UNDEFID && ncvars[ncvarid].yvarid != UNDEFID )
	    yvarid = ncvars[ncvarid].yvarid;
	  */

	  if ( xdimid != UNDEFID ) xsize = ncdims[xdimid].len;
	  if ( ydimid != UNDEFID ) ysize = ncdims[ydimid].len;

	  if ( ydimid == UNDEFID && yvarid != UNDEFID )
	    {
	      if ( ncvars[yvarid].ndims == 1 )
		{
		  ydimid = ncvars[yvarid].dimids[0];
		  ysize  = ncdims[ydimid].len;
		}
	    }

	  if ( ncvars[ncvarid].gridtype == UNDEFID || ncvars[ncvarid].gridtype == GRID_GENERIC )
	    if ( xdimid != UNDEFID && xdimid == ydimid ) ncvars[ncvarid].gridtype = GRID_UNSTRUCTURED;

	  grid_init(&grid);
	  grid_init(&proj);

	  grid.prec  = DATATYPE_FLT64;
	  grid.trunc = ncvars[ncvarid].truncation;

	  if ( ncvars[ncvarid].gridtype == GRID_TRAJECTORY )
	    {
	      if ( ncvars[ncvarid].xvarid == UNDEFID )
		Error("Longitude coordinate undefined for %s!", name);
	      if ( ncvars[ncvarid].yvarid == UNDEFID )
		Error("Latitude coordinate undefined for %s!", name);
	    }
	  else
	    {
	      size_t start[3], count[3];
	      int ltgrid = FALSE;

	      if ( xvarid != UNDEFID && yvarid != UNDEFID )
		{
		  if ( ncvars[xvarid].ndims != ncvars[yvarid].ndims )
		    {
		      Warning("Inconsistent grid structure for variable %s!", ncvars[ncvarid].name);
		      ncvars[ncvarid].xvarid = UNDEFID;
		      ncvars[ncvarid].yvarid = UNDEFID;
		      xvarid = UNDEFID;
		      yvarid = UNDEFID;
		    }

		  if ( ncvars[xvarid].ndims > 2 || ncvars[yvarid].ndims > 2 )
		    {
		      if ( ncvars[xvarid].ndims == 3 && ncvars[xvarid].dimids[0] == timedimid &&
			   ncvars[yvarid].ndims == 3 && ncvars[yvarid].dimids[0] == timedimid )
			{
			  if ( ltwarn )
			    Warning("Time varying grids unsupported, using grid at time step 1!");
			  ltgrid = TRUE;
			  ltwarn = FALSE;
			  start[0] = start[1] = start[2] = 0;
			  count[0] = 1; count[1] = ysize; count[2] = xsize;
			}
		      else
			{
			  Warning("Unsupported grid structure for variable %s (grid dims > 2)!", ncvars[ncvarid].name);
			  ncvars[ncvarid].xvarid = UNDEFID;
			  ncvars[ncvarid].yvarid = UNDEFID;
			  xvarid = UNDEFID;
			  yvarid = UNDEFID;
			}
		    }
		}

              if ( xvarid != UNDEFID )
                {
                  if ( ncvars[xvarid].ndims > 3 || (ncvars[xvarid].ndims == 3 && ltgrid == FALSE) )
                    {
                      Warning("Coordinate variable %s has to many dimensions (%d), skipped!", ncvars[xvarid].name, ncvars[xvarid].ndims);
                      //ncvars[ncvarid].xvarid = UNDEFID;
                      xvarid = UNDEFID;
                    }
                }

              if ( yvarid != UNDEFID )
                {
                  if ( ncvars[yvarid].ndims > 3 || (ncvars[yvarid].ndims == 3 && ltgrid == FALSE) )
                    {
                      Warning("Coordinate variable %s has to many dimensions (%d), skipped!", ncvars[yvarid].name, ncvars[yvarid].ndims);
                      //ncvars[ncvarid].yvarid = UNDEFID;
                      yvarid = UNDEFID;
                    }
                }

              if ( xvarid != UNDEFID )
		{
                  skipvar = TRUE;
		  islon = ncvars[xvarid].islon;
		  ndims = ncvars[xvarid].ndims;
		  if ( ndims == 2 || ndims == 3 )
		    {
		      ncvars[ncvarid].gridtype = GRID_CURVILINEAR;
		      size = xsize*ysize;
		      /* Check size of 2 dimensional coordinate variables */
		      {
			int dimid;
			size_t dimsize1, dimsize2;
			dimid = ncvars[xvarid].dimids[ndims-2];
			dimsize1 = ncdims[dimid].len;
			dimid = ncvars[xvarid].dimids[ndims-1];
			dimsize2 = ncdims[dimid].len;
			if ( dimsize1*dimsize2 == size ) skipvar = FALSE;
		      }
		    }
		  else if ( ndims == 1 )
		    {
		      size = xsize;
		      /* Check size of 1 dimensional coordinate variables */
		      {
			int dimid;
			size_t dimsize;
			dimid = ncvars[xvarid].dimids[0];
			dimsize = ncdims[dimid].len;
			if ( dimsize == size ) skipvar = FALSE;
		      }
		    }
		  else if ( ndims == 0 && xsize == 0 )
		    {
                      xsize = 1;
		      size = xsize;
                      skipvar = FALSE;
		    }

                  if ( skipvar )
                    {
                      Warning("Unsupported array structure, skipped variable %s!", ncvars[ncvarid].name);
                      ncvars[ncvarid].isvar = -1;
                      continue;
                    }

		  if ( ncvars[xvarid].xtype == NC_FLOAT ) grid.prec = DATATYPE_FLT32;
		  grid.xvals = (double *) malloc(size*sizeof(double));

		  if ( ltgrid )
		    cdf_get_vara_double(ncvars[xvarid].ncid, xvarid, start, count, grid.xvals);
		  else
		    cdf_get_var_double(ncvars[xvarid].ncid, xvarid, grid.xvals);

                  scale_add(size, grid.xvals, ncvars[xvarid].addoffset, ncvars[xvarid].scalefactor);

		  strcpy(grid.xname, ncvars[xvarid].name);
		  strcpy(grid.xlongname, ncvars[xvarid].longname);
		  strcpy(grid.xunits, ncvars[xvarid].units);
		  /* don't change the name !!! */
		  /*
		  if ( (len = strlen(grid.xname)) > 2 )
		    if ( grid.xname[len-2] == '_' && isdigit((int) grid.xname[len-1]) )
		      grid.xname[len-2] = 0;
		  */
		  if ( islon && xsize > 1 )
		    {
		      xinc = fabs(grid.xvals[0] - grid.xvals[1]);
		      for ( i = 2; i < (int) xsize; i++ )
			if ( (fabs(grid.xvals[i-1] - grid.xvals[i]) - xinc) > (xinc/1000) ) break;

		      if ( i < (int) xsize ) xinc = 0;
		    }
		}

	      if ( yvarid != UNDEFID )
		{
                  skipvar = TRUE;
		  islat = ncvars[yvarid].islat;
		  ndims = ncvars[yvarid].ndims;
		  if ( ndims == 2 || ndims == 3 )
		    {
		      ncvars[ncvarid].gridtype = GRID_CURVILINEAR;
		      size = xsize*ysize;
		      /* Check size of 2 dimensional coordinate variables */
		      {
			int dimid;
			size_t dimsize1, dimsize2;
			dimid = ncvars[yvarid].dimids[ndims-2];
			dimsize1 = ncdims[dimid].len;
			dimid = ncvars[yvarid].dimids[ndims-1];
			dimsize2 = ncdims[dimid].len;
			if ( dimsize1*dimsize2 == size ) skipvar = FALSE;
		      }
		    }
		  else if ( ndims == 1 )
		    {
		      if ( (int) ysize == 0 ) size = xsize;
		      else                    size = ysize;

		      /* Check size of 1 dimensional coordinate variables */
		      {
			int dimid;
			size_t dimsize;
			dimid = ncvars[yvarid].dimids[0];
			dimsize = ncdims[dimid].len;
			if ( dimsize == size ) skipvar = FALSE;
		      }
		    }
		  else if ( ndims == 0 && ysize == 0 )
		    {
                      ysize = 1;
		      size = ysize;
                      skipvar = FALSE;
		    }

                  if ( skipvar )
                    {
                      Warning("Unsupported array structure, skipped variable %s!", ncvars[ncvarid].name);
                      ncvars[ncvarid].isvar = -1;
                      continue;
                    }

		  if ( ncvars[yvarid].xtype == NC_FLOAT ) grid.prec = DATATYPE_FLT32;
		  grid.yvals = (double *) malloc(size*sizeof(double));

		  if ( ltgrid )
		    cdf_get_vara_double(ncvars[yvarid].ncid, yvarid, start, count, grid.yvals);
		  else
		    cdf_get_var_double(ncvars[yvarid].ncid, yvarid, grid.yvals);

                  scale_add(size, grid.yvals, ncvars[yvarid].addoffset, ncvars[yvarid].scalefactor);

		  strcpy(grid.yname, ncvars[yvarid].name);
		  strcpy(grid.ylongname, ncvars[yvarid].longname);
		  strcpy(grid.yunits, ncvars[yvarid].units);
		  /* don't change the name !!! */
		  /*
		  if ( (len = strlen(grid.yname)) > 2 )
		    if ( grid.yname[len-2] == '_' && isdigit((int) grid.yname[len-1]) )
		      grid.yname[len-2] = 0;
		  */
		  if ( islon && (int) ysize > 1 )
		    {
		      yinc = fabs(grid.yvals[0] - grid.yvals[1]);
		      for ( i = 2; i < (int) ysize; i++ )
			if ( (fabs(grid.yvals[i-1] - grid.yvals[i]) - yinc) > (yinc/1000) ) break;

		      if ( i < (int) ysize ) yinc = 0;
		    }
		}

	      if      ( (int) ysize == 0 ) size = xsize;
	      else if ( (int) xsize == 0 ) size = ysize;
	      else if ( ncvars[ncvarid].gridtype == GRID_UNSTRUCTURED ) size = xsize;
	      else                         size = xsize*ysize;
	    }

	  if ( ncvars[ncvarid].gridtype == UNDEFID ||
	       ncvars[ncvarid].gridtype == GRID_GENERIC )
	    {
	      if ( islat && islon )
		{
		  if ( isGaussGrid(ysize, yinc, grid.yvals) )
                    {
                      ncvars[ncvarid].gridtype = GRID_GAUSSIAN;
                      np = ysize/2;
                    }
                  else
		    ncvars[ncvarid].gridtype = GRID_LONLAT;
		}
	      else if ( islat && !islon && xsize == 0 )
		{
		  if ( isGaussGrid(ysize, yinc, grid.yvals) )
                    {
                      ncvars[ncvarid].gridtype = GRID_GAUSSIAN;
                      np = ysize/2;
                    }
                  else
		    ncvars[ncvarid].gridtype = GRID_LONLAT;
		}
	      else if ( islon && !islat && ysize == 0 )
		{
		  ncvars[ncvarid].gridtype = GRID_LONLAT;
		}
	      else
		ncvars[ncvarid].gridtype = GRID_GENERIC;
	    }

	  switch (ncvars[ncvarid].gridtype)
	    {
	    case GRID_GENERIC:
	    case GRID_LONLAT:
	    case GRID_GAUSSIAN:
	    case GRID_UNSTRUCTURED:
	    case GRID_CURVILINEAR:
	      {
		grid.size  = (int)size;
		grid.xsize = (int)xsize;
		grid.ysize = (int)ysize;
                grid.np    = (int)np;
		if ( xvarid != UNDEFID )
		  {
		    grid.xdef  = 1;
		    if ( ncvars[xvarid].bounds != UNDEFID )
		      {
			nbdims = ncvars[ncvars[xvarid].bounds].ndims;
			if ( nbdims == 2 || nbdims == 3 )
			  {
			    nvertex = ncdims[ncvars[ncvars[xvarid].bounds].dimids[nbdims-1]].len;
			    grid.nvertex = (int) nvertex;
			    grid.xbounds = (double *) malloc(nvertex*size*sizeof(double));
			    cdf_get_var_double(ncvars[xvarid].ncid, ncvars[xvarid].bounds, grid.xbounds);
			  }
		      }
		  }
		if ( yvarid != UNDEFID )
		  {
		    grid.ydef  = 1;
		    if ( ncvars[yvarid].bounds != UNDEFID )
		      {
			nbdims = ncvars[ncvars[yvarid].bounds].ndims;
			if ( nbdims == 2 || nbdims == 3 )
			  {
			    nvertex = ncdims[ncvars[ncvars[yvarid].bounds].dimids[nbdims-1]].len;
			    /*
			    if ( nvertex != grid.nvertex )
			      Warning("nvertex problem! nvertex x %d, nvertex y %d",
				      grid.nvertex, (int) nvertex);
			    */
			    grid.ybounds = (double *) malloc(nvertex*size*sizeof(double));
			    cdf_get_var_double(ncvars[yvarid].ncid, ncvars[yvarid].bounds, grid.ybounds);
			  }
		      }
		  }

		if ( ncvars[ncvarid].cellarea != UNDEFID )
		  {
		    grid.area = (double *) malloc(size*sizeof(double));
		    cdf_get_var_double(ncvars[ncvarid].ncid, ncvars[ncvarid].cellarea, grid.area);
		  }

		break;
	      }
	    case GRID_SPECTRAL:
	      {
		grid.size = (int)size;
		grid.lcomplex = 1;
		break;
	      }
	    case GRID_FOURIER:
	      {
		grid.size = (int)size;
		break;
	      }
	    case GRID_TRAJECTORY:
	      {
		grid.size = 1;
		break;
	      }
	    }

	  grid.type = ncvars[ncvarid].gridtype;

	  if ( grid.size == 0 )
	    {
	      if ( (ncvars[ncvarid].ndims == 1 && ncvars[ncvarid].dimtype[0] == T_AXIS) ||
		   (ncvars[ncvarid].ndims == 1 && ncvars[ncvarid].dimtype[0] == Z_AXIS) ||
		   (ncvars[ncvarid].ndims == 2 && ncvars[ncvarid].dimtype[0] == T_AXIS && ncvars[ncvarid].dimtype[1] == Z_AXIS) )
		{
		  grid.type  = GRID_GENERIC;
		  grid.size  = 1;
		  grid.xsize = 0;
		  grid.ysize = 0;
		}
	      else
		{
		  Warning("Variable %s has an unsupported grid, skipped!", ncvars[ncvarid].name);
		  ncvars[ncvarid].isvar = -1;
		  continue;
		}
	    }

	  if ( number_of_grid_used != UNDEFID && (grid.type == UNDEFID || grid.type == GRID_GENERIC) )
            grid.type   = GRID_UNSTRUCTURED;

	  if ( number_of_grid_used != UNDEFID && grid.type == GRID_UNSTRUCTURED )
            grid.number = number_of_grid_used;

	  if ( ncvars[ncvarid].gmapid >= 0 && ncvars[ncvarid].gridtype != GRID_CURVILINEAR )
	    {
	      cdf_inq_varnatts(ncvars[ncvarid].ncid, ncvars[ncvarid].gmapid, &nvatts);

	      for ( iatt = 0; iatt < nvatts; iatt++ )
		{
		  cdf_inq_attname(ncvars[ncvarid].ncid, ncvars[ncvarid].gmapid, iatt, attname);
		  cdf_inq_attlen(ncvars[ncvarid].ncid, ncvars[ncvarid].gmapid, attname, &attlen);

		  if ( strcmp(attname, "grid_mapping_name") == 0 )
		    {
                      enum {
                        attstringlen = 8192,
                      };
                      char attstring[attstringlen];

		      cdfGetAttText(ncvars[ncvarid].ncid, ncvars[ncvarid].gmapid, attname, attstringlen, attstring);
		      strtolower(attstring);

		      if ( strcmp(attstring, "rotated_latitude_longitude") == 0 )
			grid.isRotated = TRUE;
		      else if ( strcmp(attstring, "sinusoidal") == 0 )
			grid.type = GRID_SINUSOIDAL;
		      else if ( strcmp(attstring, "lambert_azimuthal_equal_area") == 0 )
			grid.type = GRID_LAEA;
		      else if ( strcmp(attstring, "lambert_conformal_conic") == 0 )
			grid.type = GRID_LCC2;
		      else if ( strcmp(attstring, "lambert_cylindrical_equal_area") == 0 )
			{
			  proj.type = GRID_PROJECTION;
			  proj.name = strdup(attstring);
			}
		    }
		  else if ( strcmp(attname, "earth_radius") == 0 )
		    {
		      cdfGetAttDouble(ncvars[ncvarid].ncid, ncvars[ncvarid].gmapid, attname, 1, &datt);
		      grid.laea_a = datt;
		      grid.lcc2_a = datt;
		    }
		  else if ( strcmp(attname, "longitude_of_projection_origin") == 0 )
		    {
		      cdfGetAttDouble(ncvars[ncvarid].ncid, ncvars[ncvarid].gmapid, attname, 1, &grid.laea_lon_0);
		    }
		  else if ( strcmp(attname, "longitude_of_central_meridian") == 0 )
		    {
		      cdfGetAttDouble(ncvars[ncvarid].ncid, ncvars[ncvarid].gmapid, attname, 1, &grid.lcc2_lon_0);
		    }
		  else if ( strcmp(attname, "latitude_of_projection_origin") == 0 )
		    {
		      cdfGetAttDouble(ncvars[ncvarid].ncid, ncvars[ncvarid].gmapid, attname, 1, &datt);
		      grid.laea_lat_0 = datt;
		      grid.lcc2_lat_0 = datt;
		    }
		  else if ( strcmp(attname, "standard_parallel") == 0 )
		    {
		      if ( attlen == 1 )
			{
			  cdfGetAttDouble(ncvars[ncvarid].ncid, ncvars[ncvarid].gmapid, attname, 1, &datt);
			  grid.lcc2_lat_1 = datt;
			  grid.lcc2_lat_2 = datt;
			}
		      else
			{
			  double datt2[2];
			  cdfGetAttDouble(ncvars[ncvarid].ncid, ncvars[ncvarid].gmapid, attname, 2, datt2);
			  grid.lcc2_lat_1 = datt2[0];
			  grid.lcc2_lat_2 = datt2[1];
			}
		    }
		  else if ( strcmp(attname, "grid_north_pole_latitude") == 0 )
		    {
		      cdfGetAttDouble(ncvars[ncvarid].ncid, ncvars[ncvarid].gmapid, attname, 1, &grid.ypole);
		    }
		  else if ( strcmp(attname, "grid_north_pole_longitude") == 0 )
		    {
		      cdfGetAttDouble(ncvars[ncvarid].ncid, ncvars[ncvarid].gmapid, attname, 1, &grid.xpole);
		    }
		  else if ( strcmp(attname, "north_pole_grid_longitude") == 0 )
		    {
		      cdfGetAttDouble(ncvars[ncvarid].ncid, ncvars[ncvarid].gmapid, attname, 1, &grid.angle);
		    }
		}
	    }

          if ( grid.type == GRID_UNSTRUCTURED )
            {
              int zdimid = UNDEFID;
              int xdimidx = -1, ydimidx = -1;

              for ( i = 0; i < ndims; i++ )
                {
                  if      ( ncvars[ncvarid].dimtype[i] == X_AXIS ) xdimidx = i;
                  else if ( ncvars[ncvarid].dimtype[i] == Y_AXIS ) ydimidx = i;
                  else if ( ncvars[ncvarid].dimtype[i] == Z_AXIS ) zdimid = ncvars[ncvarid].dimids[i];
                }

              if ( xdimid != UNDEFID && ydimid != UNDEFID && zdimid == UNDEFID )
                {
                  if ( grid.xsize > grid.ysize && grid.ysize < 1000 )
                    {
                      ncvars[ncvarid].dimtype[ydimidx] = Z_AXIS;
                      ydimid = UNDEFID;
                      grid.size  = grid.xsize;
                      grid.ysize = 0;
                    }
                  else if ( grid.ysize > grid.xsize && grid.xsize < 1000 )
                    {
                      ncvars[ncvarid].dimtype[xdimidx] = Z_AXIS;
                      xdimid = ydimid;
                      ydimid = UNDEFID;
                      grid.size  = grid.ysize;
                      grid.xsize = grid.ysize;
                      grid.ysize = 0;
                    }
                }

              if ( grid.size != grid.xsize )
                {
                  Warning("Unsupported array structure, skipped variable %s!", ncvars[ncvarid].name);
                  ncvars[ncvarid].isvar = -1;
                  continue;
                }

              if ( ncvars[ncvarid].position > 0 ) grid.position = ncvars[ncvarid].position;
              if ( uuidOfHGrid[0] != 0 ) memcpy(grid.uuid, uuidOfHGrid, 16);
            }

#if defined (PROJECTION_TEST)
	  if ( proj.type == GRID_PROJECTION )
	    {
	      if ( grid.type == GRID_GENERIC )
		{
		  grid.type = GRID_CURVILINEAR;
		}

	      if ( grid.type == GRID_CURVILINEAR )
		{
		  proj.size  = grid.size;
		  proj.xsize = grid.xsize;
                  proj.ysize = grid.ysize;
		}

	      //  grid.proj = gridGenerate(proj);
	    }
#endif

	  if ( CDI_Debug )
	    {
	      Message("grid: type = %d, size = %d, nx = %d, ny %d",
		      grid.type, grid.size, grid.xsize, grid.ysize);
	      Message("proj: type = %d, size = %d, nx = %d, ny %d",
		      proj.type, proj.size, proj.xsize, proj.ysize);
	    }

#if defined (PROJECTION_TEST)
	  if ( proj.type == GRID_PROJECTION )
	    {
	      ncvars[ncvarid].gridID = varDefGrid(vlistID, &proj, 1);
	      copy_numeric_projatts(ncvars[ncvarid].gridID, ncvars[ncvarid].gmapid, ncvars[ncvarid].ncid);
	    }
	  else
#endif
	    ncvars[ncvarid].gridID = varDefGrid(vlistID, &grid, 1);

          if ( grid.type == GRID_UNSTRUCTURED )
            {
              if ( gridfile[0] != 0 ) gridDefReference(ncvars[ncvarid].gridID, gridfile);
            }

          if ( ncvars[ncvarid].chunked ) grid_set_chunktype(&grid, &ncvars[ncvarid]);

	  gridindex = vlistGridIndex(vlistID, ncvars[ncvarid].gridID);
	  streamptr->xdimID[gridindex] = xdimid;
	  streamptr->ydimID[gridindex] = ydimid;
          if ( xdimid == -1 && ydimid == -1 && grid.size == 1 )
            gridDefHasDims(ncvars[ncvarid].gridID, FALSE);

	  if ( CDI_Debug )
	    Message("gridID %d %d %s", ncvars[ncvarid].gridID, ncvarid, ncvars[ncvarid].name);

	  for ( ncvarid2 = ncvarid+1; ncvarid2 < nvars; ncvarid2++ )
	    if ( ncvars[ncvarid2].isvar == TRUE && ncvars[ncvarid2].gridID == UNDEFID )
	      {
		int xdimid2 = UNDEFID, ydimid2 = UNDEFID, zdimid2 = UNDEFID;
                int xdimidx = -1, ydimidx = -1;
		int ndims2 = ncvars[ncvarid2].ndims;

		for ( i = 0; i < ndims2; i++ )
		  {
		    if ( ncvars[ncvarid2].dimtype[i] == X_AXIS )
		      { xdimid2 = ncvars[ncvarid2].dimids[i]; xdimidx = i; }
		    else if ( ncvars[ncvarid2].dimtype[i] == Y_AXIS )
		      { ydimid2 = ncvars[ncvarid2].dimids[i]; ydimidx = i; }
		    else if ( ncvars[ncvarid2].dimtype[i] == Z_AXIS )
		      { zdimid2 = ncvars[ncvarid2].dimids[i]; }
		  }

                if ( ncvars[ncvarid2].gridtype == UNDEFID && grid.type == GRID_UNSTRUCTURED )
                  {
                    if ( xdimid == xdimid2 && ydimid2 != UNDEFID && zdimid2 == UNDEFID )
                      {
                        ncvars[ncvarid2].dimtype[ydimidx] = Z_AXIS;
                        ydimid2 = UNDEFID;
                      }

                    if ( xdimid == ydimid2 && xdimid2 != UNDEFID && zdimid2 == UNDEFID )
                      {
                        ncvars[ncvarid2].dimtype[xdimidx] = Z_AXIS;
                        xdimid2 = ydimid2;
                        ydimid2 = UNDEFID;
                      }
                  }

                if ( xdimid == xdimid2 &&
		    (ydimid == ydimid2 || (xdimid == ydimid && ydimid2 == UNDEFID)) )
		  {
		    int same_grid = TRUE;
                    /*
		    if ( xvarid != -1 && ncvars[ncvarid2].xvarid != UNDEFID &&
			 xvarid != ncvars[ncvarid2].xvarid ) same_grid = FALSE;

		    if ( yvarid != -1 && ncvars[ncvarid2].yvarid != UNDEFID &&
			 yvarid != ncvars[ncvarid2].yvarid ) same_grid = FALSE;
                    */
		    if ( ncvars[ncvarid].xvarid != ncvars[ncvarid2].xvarid ) same_grid = FALSE;
		    if ( ncvars[ncvarid].yvarid != ncvars[ncvarid2].yvarid ) same_grid = FALSE;

		    if ( ncvars[ncvarid].position != ncvars[ncvarid2].position ) same_grid = FALSE;

		    if ( same_grid )
		      {
			if ( CDI_Debug )
			  Message("Same gridID %d %d %s", ncvars[ncvarid].gridID, ncvarid2, ncvars[ncvarid2].name);
			ncvars[ncvarid2].gridID = ncvars[ncvarid].gridID;
			ncvars[ncvarid2].chunktype = ncvars[ncvarid].chunktype;
		      }
		  }
	      }

	  grid_free(&grid);
	  grid_free(&proj);
	}
    }
}

/* define all input zaxes */
static
void define_all_zaxes(stream_t *streamptr, int vlistID, ncdim_t *ncdims, int nvars, ncvar_t *ncvars,
		      size_t vctsize, double *vct, unsigned char *uuidOfVGrid)
{
  int ncvarid, ncvarid2;
  int i, ilev, ndims;
  int zaxisindex;
  int zprec;
  int nbdims, nvertex, nlevel;
  int positive = 0;
  char *pname, *plongname, *punits;

  for ( ncvarid = 0; ncvarid < nvars; ncvarid++ )
    {
      if ( ncvars[ncvarid].isvar == TRUE && ncvars[ncvarid].zaxisID == UNDEFID )
	{
	  int with_bounds = FALSE;
	  int zdimid = UNDEFID;
	  int zvarid = UNDEFID;
	  int zsize = 1;
	  double *lbounds = NULL;
	  double *ubounds = NULL;
	  int zaxisType;

          positive = 0;

	  ndims = ncvars[ncvarid].ndims;
	  for ( i = 0; i < ndims; i++ )
	    {
	      if ( ncvars[ncvarid].dimtype[i] == Z_AXIS )
		zdimid = ncvars[ncvarid].dimids[i];
	    }

	  if ( zdimid != UNDEFID )
	    {
	      zvarid = ncdims[zdimid].ncvarid;
	      zsize  = (int)ncdims[zdimid].len;
	    }

	  if ( CDI_Debug ) Message("nlevs = %d", zsize);

	  double *zvar = (double *)xmalloc((size_t)zsize * sizeof (double));

	  zaxisType = UNDEFID;

	  if ( zvarid != UNDEFID ) zaxisType = ncvars[zvarid].zaxistype;

	  if ( zaxisType == UNDEFID )  zaxisType = ZAXIS_GENERIC;

	  zprec = DATATYPE_FLT64;

	  if ( zvarid != UNDEFID )
	    {
	      positive  = ncvars[zvarid].positive;
	      pname     = ncvars[zvarid].name;
	      plongname = ncvars[zvarid].longname;
	      punits    = ncvars[zvarid].units;
	      if ( ncvars[zvarid].xtype == NC_FLOAT ) zprec = DATATYPE_FLT32;
	      /* don't change the name !!! */
	      /*
	      if ( (len = strlen(pname)) > 2 )
		if ( pname[len-2] == '_' && isdigit((int) pname[len-1]) )
		  pname[len-2] = 0;
	      */
	      cdf_get_var_double(ncvars[zvarid].ncid, zvarid, zvar);

	      if ( ncvars[zvarid].bounds != UNDEFID )
		{
		  nbdims = ncvars[ncvars[zvarid].bounds].ndims;
		  if ( nbdims == 2 )
		    {
		      nlevel  = (int)ncdims[ncvars[ncvars[zvarid].bounds].dimids[0]].len;
		      nvertex = (int)ncdims[ncvars[ncvars[zvarid].bounds].dimids[1]].len;
		      if ( nlevel == zsize && nvertex == 2 )
			{
			  double *zbounds;
			  with_bounds = TRUE;
			  zbounds = (double *) malloc(2*(size_t)nlevel*sizeof(double));
			  lbounds = (double *) malloc((size_t)nlevel*sizeof(double));
			  ubounds = (double *) malloc((size_t)nlevel*sizeof(double));
			  cdf_get_var_double(ncvars[zvarid].ncid, ncvars[zvarid].bounds, zbounds);
			  for ( i = 0; i < nlevel; ++i )
			    {
			      lbounds[i] = zbounds[i*2];
			      ubounds[i] = zbounds[i*2+1];
			    }
			  free(zbounds);
			}
		    }
		}
	    }
	  else
	    {
	      pname     = NULL;
	      plongname = NULL;
	      punits    = NULL;

	      if ( zsize == 1 )
		{
                  if ( ncvars[ncvarid].zaxistype != UNDEFID )
                    zaxisType = ncvars[ncvarid].zaxistype;
                  else
                    zaxisType = ZAXIS_SURFACE;

		  zvar[0] = 0;
		  /*
		  if ( zdimid == UNDEFID )
		    zvar[0] = 9999;
		  else
		    zvar[0] = 0;
		  */
		}
	      else
		{
		  for ( ilev = 0; ilev < zsize; ilev++ ) zvar[ilev] = ilev + 1;
		}
	    }

      	  ncvars[ncvarid].zaxisID = varDefZaxis(vlistID, zaxisType, (int) zsize, zvar, with_bounds, lbounds, ubounds,
						(int)vctsize, vct, pname, plongname, punits, zprec, 1, 0);

	  if ( uuidOfVGrid[0] != 0 )
            {
              // printf("uuidOfVGrid: defined\n");
              zaxisDefUUID(ncvars[ncvarid].zaxisID, uuidOfVGrid);
            }

          if ( positive > 0 ) zaxisDefPositive(ncvars[ncvarid].zaxisID, positive);

	  free(zvar);
	  free(lbounds);
	  free(ubounds);

	  zaxisindex = vlistZaxisIndex(vlistID, ncvars[ncvarid].zaxisID);
	  streamptr->zaxisID[zaxisindex]  = zdimid;

	  if ( CDI_Debug )
	    Message("zaxisID %d %d %s", ncvars[ncvarid].zaxisID, ncvarid, ncvars[ncvarid].name);

	  for ( ncvarid2 = ncvarid+1; ncvarid2 < nvars; ncvarid2++ )
	    if ( ncvars[ncvarid2].isvar == TRUE && ncvars[ncvarid2].zaxisID == UNDEFID && ncvars[ncvarid2].zaxistype == UNDEFID )
	      {
		int zdimid2 = -1;
		ndims = ncvars[ncvarid2].ndims;
		for ( i = 0; i < ndims; i++ )
		  {
		    if ( ncvars[ncvarid2].dimtype[i] == Z_AXIS )
		      zdimid2 = ncvars[ncvarid2].dimids[i];
		  }
		if ( zdimid == zdimid2 )
		  {
		    if ( CDI_Debug )
		      Message("zaxisID %d %d %s",
			      ncvars[ncvarid].zaxisID, ncvarid2, ncvars[ncvarid2].name);
		    ncvars[ncvarid2].zaxisID = ncvars[ncvarid].zaxisID;
		  }
	      }
	}
    }
}

struct varinfo
{
  int      ncvarid;
  const char *name;
};

static
int cmpvarname(const void *s1, const void *s2)
{
  const struct varinfo *x = (const struct varinfo *)s1,
    *y = (const struct varinfo *)s2;
  return (strcmp(x->name, y->name));
}

/* define all input data variables */
static
void define_all_vars(stream_t *streamptr, int vlistID, int instID, int modelID, int *varids, int nvars, int num_ncvars, ncvar_t *ncvars)
{
  if ( streamptr->sortname )
    {
      struct varinfo *varInfo
        = (struct varinfo *)xmalloc((size_t)nvars * sizeof (struct varinfo));

      for ( int varID = 0; varID < nvars; varID++ )
	{
	  int ncvarid = varids[varID];
	  varInfo[varID].ncvarid = ncvarid;
	  varInfo[varID].name = ncvars[ncvarid].name;
	}
      qsort(varInfo, (size_t)nvars, sizeof(varInfo[0]), cmpvarname);
      for ( int varID = 0; varID < nvars; varID++ )
	{
	  varids[varID] = varInfo[varID].ncvarid;
	}
      free(varInfo);
    }

  for ( int varID1 = 0; varID1 < nvars; varID1++ )
    {
      int gridID, zaxisID;

      int ncvarid = varids[varID1];
      gridID  = ncvars[ncvarid].gridID;
      zaxisID = ncvars[ncvarid].zaxisID;

      stream_new_var(streamptr, gridID, zaxisID, CDI_UNDEFID);
      int varID = vlistDefVar(vlistID, gridID, zaxisID, ncvars[ncvarid].tsteptype);

#if  defined  (HAVE_NETCDF4)
      if ( ncvars[ncvarid].deflate )
	vlistDefVarCompType(vlistID, varID, COMPRESS_ZIP);

      if ( ncvars[ncvarid].chunked && ncvars[ncvarid].chunktype != UNDEFID )
        vlistDefVarChunkType(vlistID, varID, ncvars[ncvarid].chunktype);
#endif

      streamptr->vars[varID1].defmiss = 0;
      streamptr->vars[varID1].ncvarid = ncvarid;

      vlistDefVarName(vlistID, varID, ncvars[ncvarid].name);
      if ( ncvars[ncvarid].param != UNDEFID ) vlistDefVarParam(vlistID, varID, ncvars[ncvarid].param);
      if ( ncvars[ncvarid].code != UNDEFID )  vlistDefVarCode(vlistID, varID, ncvars[ncvarid].code);
      if ( ncvars[ncvarid].code != UNDEFID )
	{
	  int param;
	  param = cdiEncodeParam(ncvars[ncvarid].code, ncvars[ncvarid].tabnum, 255);
	  vlistDefVarParam(vlistID, varID, param);
	}
      if ( ncvars[ncvarid].longname[0] )  vlistDefVarLongname(vlistID, varID, ncvars[ncvarid].longname);
      if ( ncvars[ncvarid].stdname[0] )   vlistDefVarStdname(vlistID, varID, ncvars[ncvarid].stdname);
      if ( ncvars[ncvarid].units[0] )     vlistDefVarUnits(vlistID, varID, ncvars[ncvarid].units);

      if ( ncvars[ncvarid].lvalidrange )
        vlistDefVarValidrange(vlistID, varID, ncvars[ncvarid].validrange);

      if ( IS_NOT_EQUAL(ncvars[ncvarid].addoffset, 0) )
	vlistDefVarAddoffset(vlistID, varID, ncvars[ncvarid].addoffset);
      if ( IS_NOT_EQUAL(ncvars[ncvarid].scalefactor, 1) )
	vlistDefVarScalefactor(vlistID, varID, ncvars[ncvarid].scalefactor);

      vlistDefVarDatatype(vlistID, varID, cdfInqDatatype(ncvars[ncvarid].xtype, ncvars[ncvarid].lunsigned));

      vlistDefVarInstitut(vlistID, varID, instID);
      vlistDefVarModel(vlistID, varID, modelID);
      if ( ncvars[ncvarid].tableID != UNDEFID )
	vlistDefVarTable(vlistID, varID, ncvars[ncvarid].tableID);

      if ( ncvars[ncvarid].deffillval == FALSE && ncvars[ncvarid].defmissval == TRUE )
        {
          ncvars[ncvarid].deffillval = TRUE;
          ncvars[ncvarid].fillval    = ncvars[ncvarid].missval;
        }

      if ( ncvars[ncvarid].deffillval == TRUE )
        vlistDefVarMissval(vlistID, varID, ncvars[ncvarid].fillval);

      if ( CDI_Debug )
	Message("varID = %d  gridID = %d  zaxisID = %d", varID,
		vlistInqVarGrid(vlistID, varID), vlistInqVarZaxis(vlistID, varID));

      int gridindex = vlistGridIndex(vlistID, gridID);
      int xdimid = streamptr->xdimID[gridindex];
      int ydimid = streamptr->ydimID[gridindex];

      int zaxisindex = vlistZaxisIndex(vlistID, zaxisID);
      int zdimid = streamptr->zaxisID[zaxisindex];

      int ndims = ncvars[ncvarid].ndims;
      int iodim = 0;
      int ixyz = 0;
      int ipow10[4] = {1, 10, 100, 1000};

      if ( ncvars[ncvarid].tsteptype != TSTEP_CONSTANT ) iodim++;

      if ( gridInqType(gridID) == GRID_UNSTRUCTURED && ndims-iodim <= 2 && ydimid == xdimid )
        {
          if ( xdimid == ncvars[ncvarid].dimids[ndims-1] )
            {
              ixyz = 321;
            }
          else
            {
              ixyz = 213;
            }
        }
      else
        {
          for ( int idim = iodim; idim < ndims; idim++ )
            {
              if      ( xdimid == ncvars[ncvarid].dimids[idim] )
                ixyz += 1*ipow10[ndims-idim-1];
              else if ( ydimid == ncvars[ncvarid].dimids[idim] )
                ixyz += 2*ipow10[ndims-idim-1];
              else if ( zdimid == ncvars[ncvarid].dimids[idim] )
                ixyz += 3*ipow10[ndims-idim-1];
            }
        }

      vlistDefVarXYZ(vlistID, varID, ixyz);
      /*
      printf("ixyz %d\n", ixyz);
      printf("ndims %d\n", ncvars[ncvarid].ndims);
      for ( int i = 0; i < ncvars[ncvarid].ndims; ++i )
        printf("dimids: %d %d\n", i, ncvars[ncvarid].dimids[i]);
      printf("xdimid, ydimid %d %d\n", xdimid, ydimid);
      */
      if ( ncvars[ncvarid].ensdata != NULL )
        {
          vlistDefVarEnsemble( vlistID, varID, ncvars[ncvarid].ensdata->ens_index,
                               ncvars[ncvarid].ensdata->ens_count,
                               ncvars[ncvarid].ensdata->forecast_init_type );
          free(ncvars[ncvarid].ensdata);
          ncvars[ncvarid].ensdata = NULL;
        }

      if ( ncvars[ncvarid].extra[0] != 0 )
        {
          vlistDefVarExtra(vlistID, varID, ncvars[ncvarid].extra);
        }
    }

  for ( int varID = 0; varID < nvars; varID++ )
    {
      int ncvarid = varids[varID];
      int ncid = ncvars[ncvarid].ncid;

      if ( ncvars[ncvarid].natts )
	{
	  int attnum;
	  int iatt;
	  nc_type attrtype;
	  size_t attlen;
	  char attname[CDI_MAX_NAME];
	  const int attstringlen = 8192; char attstring[8192];
	  int nvatts = ncvars[ncvarid].natts;

	  for ( iatt = 0; iatt < nvatts; iatt++ )
	    {
	      attnum = ncvars[ncvarid].atts[iatt];
	      cdf_inq_attname(ncid, ncvarid, attnum, attname);
	      cdf_inq_attlen(ncid, ncvarid, attname, &attlen);
	      cdf_inq_atttype(ncid, ncvarid, attname, &attrtype);

	      if ( attrtype == NC_SHORT || attrtype == NC_INT )
		{
		  int *attint = (int *) malloc(attlen*sizeof(int));
		  cdfGetAttInt(ncid, ncvarid, attname, (int)attlen, attint);
		  if ( attrtype == NC_SHORT )
		    vlistDefAttInt(vlistID, varID, attname, DATATYPE_INT16, (int)attlen, attint);
		  else
		    vlistDefAttInt(vlistID, varID, attname, DATATYPE_INT32, (int)attlen, attint);
		  free(attint);
		}
	      else if ( attrtype == NC_FLOAT || attrtype == NC_DOUBLE )
		{
		  double *attflt = (double *) malloc(attlen*sizeof(double));
		  cdfGetAttDouble(ncid, ncvarid, attname, (int)attlen, attflt);
		  if ( attrtype == NC_FLOAT )
		    vlistDefAttFlt(vlistID, varID, attname, DATATYPE_FLT32, (int)attlen, attflt);
		  else
		    vlistDefAttFlt(vlistID, varID, attname, DATATYPE_FLT64, (int)attlen, attflt);
		  free(attflt);
		}
	      else if ( xtypeIsText(attrtype) )
		{
		  cdfGetAttText(ncid, ncvarid, attname, attstringlen, attstring);
		  vlistDefAttTxt(vlistID, varID, attname, (int)attlen, attstring);
		}
	      else
		{
		  if ( CDI_Debug ) printf("att: %s.%s = unknown\n", ncvars[ncvarid].name, attname);
		}
	    }

	  free(ncvars[ncvarid].atts);
          ncvars[ncvarid].atts = NULL;
	}
    }

  /* release mem of not freed attributes */
  for ( int ncvarid = 0; ncvarid < num_ncvars; ncvarid++ )
    if ( ncvars[ncvarid].atts ) free(ncvars[ncvarid].atts);

  if ( varids ) free(varids);

  for ( int varID = 0; varID < nvars; varID++ )
    {
      if ( vlistInqVarCode(vlistID, varID) == -varID-1 )
	{
	  const char *pname = vlistInqVarNamePtr(vlistID, varID);
	  size_t len = strlen(pname);
	  if ( len > 3 && isdigit((int) pname[3]) )
	    {
	      if ( memcmp("var", pname, 3) == 0 )
		{
		  vlistDefVarCode(vlistID, varID, atoi(pname+3));
		  vlistDestroyVarName(vlistID, varID);
		}
	    }
	  else if ( len > 4 && isdigit((int) pname[4]) )
	    {
	      if ( memcmp("code", pname, 4) == 0 )
		{
		  vlistDefVarCode(vlistID, varID, atoi(pname+4));
		  vlistDestroyVarName(vlistID, varID);
		}
	    }
	  else if ( len > 5 && isdigit((int) pname[5]) )
	    {
	      if ( memcmp("param", pname, 5) == 0 )
		{
		  int pnum = -1, pcat = 255, pdis = 255;
		  sscanf(pname+5, "%d.%d.%d", &pnum, &pcat, &pdis);
		  vlistDefVarParam(vlistID, varID, cdiEncodeParam(pnum, pcat, pdis));
		  vlistDestroyVarName(vlistID, varID);
		}
	    }
	}
    }

  for ( int varID = 0; varID < nvars; varID++ )
    {
      int instID  = vlistInqVarInstitut(vlistID, varID);
      int modelID = vlistInqVarModel(vlistID, varID);
      int tableID = vlistInqVarTable(vlistID, varID);
      int code = vlistInqVarCode(vlistID, varID);
      if ( cdiDefaultTableID != UNDEFID )
	{
	  if ( tableInqParNamePtr(cdiDefaultTableID, code) )
	    {
	      vlistDestroyVarName(vlistID, varID);
	      vlistDestroyVarLongname(vlistID, varID);
	      vlistDestroyVarUnits(vlistID, varID);

	      if ( tableID != UNDEFID )
		{
		  vlistDefVarName(vlistID, varID, tableInqParNamePtr(cdiDefaultTableID, code));
		  if ( tableInqParLongnamePtr(cdiDefaultTableID, code) )
		    vlistDefVarLongname(vlistID, varID, tableInqParLongnamePtr(cdiDefaultTableID, code));
		  if ( tableInqParUnitsPtr(cdiDefaultTableID, code) )
		    vlistDefVarUnits(vlistID, varID, tableInqParUnitsPtr(cdiDefaultTableID, code));
		}
	      else
		{
		  tableID = cdiDefaultTableID;
		}
	    }

	  if ( cdiDefaultModelID != UNDEFID ) modelID = cdiDefaultModelID;
	  if ( cdiDefaultInstID  != UNDEFID ) instID  = cdiDefaultInstID;
	}
      if ( instID  != UNDEFID ) vlistDefVarInstitut(vlistID, varID, instID);
      if ( modelID != UNDEFID ) vlistDefVarModel(vlistID, varID, modelID);
      if ( tableID != UNDEFID ) vlistDefVarTable(vlistID, varID, tableID);
    }
}

static
void scan_global_attributes(int fileID, int vlistID, stream_t *streamptr, int ngatts, int *instID, int *modelID, int *ucla_les, unsigned char *uuidOfHGrid, unsigned char *uuidOfVGrid, char *gridfile, int *number_of_grid_used)
{
  nc_type xtype;
  size_t attlen;
  char attname[CDI_MAX_NAME];
  enum { attstringlen = 65636 };
  char attstring[attstringlen];
  int iatt;

  for ( iatt = 0; iatt < ngatts; iatt++ )
    {
      cdf_inq_attname(fileID, NC_GLOBAL, iatt, attname);
      cdf_inq_atttype(fileID, NC_GLOBAL, attname, &xtype);
      cdf_inq_attlen(fileID, NC_GLOBAL, attname, &attlen);

      if ( xtypeIsText(xtype) )
	{
	  cdfGetAttText(fileID, NC_GLOBAL, attname, attstringlen, attstring);

          size_t attstrlen = strlen(attstring);

	  if ( attlen > 0 && attstring[0] != 0 )
	    {
	      if ( strcmp(attname, "history") == 0 )
		{
		  streamptr->historyID = iatt;
		}
	      else if ( strcmp(attname, "institution") == 0 )
		{
		  *instID = institutInq(0, 0, NULL, attstring);
		  if ( *instID == UNDEFID )
		    *instID = institutDef(0, 0, NULL, attstring);
		}
	      else if ( strcmp(attname, "source") == 0 )
		{
		  *modelID = modelInq(-1, 0, attstring);
		  if ( *modelID == UNDEFID )
		    *modelID = modelDef(-1, 0, attstring);
		}
	      else if ( strcmp(attname, "Source") == 0 )
		{
		  if ( strncmp(attstring, "UCLA-LES", 8) == 0 )
		    *ucla_les = TRUE;
		}
	      /*
	      else if ( strcmp(attname, "Conventions") == 0 )
		{
		}
	      */
	      else if ( strcmp(attname, "CDI") == 0 )
		{
		}
	      else if ( strcmp(attname, "CDO") == 0 )
		{
		}
              /*
	      else if ( strcmp(attname, "forecast_reference_time") == 0 )
		{
                  memcpy(fcreftime, attstring, attstrlen+1);
		}
              */
	      else if ( strcmp(attname, "grid_file_uri") == 0 )
		{
                  memcpy(gridfile, attstring, attstrlen+1);
		}
	      else if ( strcmp(attname, "uuidOfHGrid") == 0 && attstrlen == 36 )
		{
                  attstring[36] = 0;
                  str2uuid(attstring, uuidOfHGrid);
                  //   printf("uuid: %d %s\n", attlen, attstring);
		}
	      else if ( strcmp(attname, "uuidOfVGrid") == 0 && attstrlen == 36 )
		{
                  attstring[36] = 0;
                  str2uuid(attstring, uuidOfVGrid);
		}
	      else
		{
                  if ( strcmp(attname, "ICON_grid_file_uri") == 0 && gridfile[0] == 0 )
                    {
                      memcpy(gridfile, attstring, attstrlen+1);
                    }

		  vlistDefAttTxt(vlistID, CDI_GLOBAL, attname, (int)attstrlen, attstring);
		}
	    }
	}
      else if ( xtype == NC_SHORT || xtype == NC_INT )
	{
	  if ( strcmp(attname, "number_of_grid_used") == 0 )
	    {
	      (*number_of_grid_used) = UNDEFID;
	      cdfGetAttInt(fileID, NC_GLOBAL, attname, 1, number_of_grid_used);
	    }
 	  else
            {
              int *attint;
              attint = (int *) malloc(attlen*sizeof(int));
              cdfGetAttInt(fileID, NC_GLOBAL, attname, (int)attlen, attint);
              if ( xtype == NC_SHORT )
                vlistDefAttInt(vlistID, CDI_GLOBAL, attname, DATATYPE_INT16, (int)attlen, attint);
              else
                vlistDefAttInt(vlistID, CDI_GLOBAL, attname, DATATYPE_INT32, (int)attlen, attint);
              free(attint);
            }
        }
      else if ( xtype == NC_FLOAT || xtype == NC_DOUBLE )
	{
	  double *attflt;
	  attflt = (double *) malloc(attlen*sizeof(double));
	  cdfGetAttDouble(fileID, NC_GLOBAL, attname, (int)attlen, attflt);
	  if ( xtype == NC_FLOAT )
	    vlistDefAttFlt(vlistID, CDI_GLOBAL, attname, DATATYPE_FLT32, (int)attlen, attflt);
	  else
	    vlistDefAttFlt(vlistID, CDI_GLOBAL, attname, DATATYPE_FLT64, (int)attlen, attflt);
	  free(attflt);
	}
    }
}

static
int find_leadtime(int nvars, ncvar_t *ncvars)
{
  int leadtime_id = UNDEFID;

  for ( int ncvarid = 0; ncvarid < nvars; ncvarid++ )
    {
      if ( ncvars[ncvarid].stdname[0] )
        {
          if ( strcmp(ncvars[ncvarid].stdname, "forecast_period") == 0 )
            {
              leadtime_id = ncvarid;
              break;
            }
        }
    }

  return (leadtime_id);
}

static
void find_time_vars(int nvars, ncvar_t *ncvars, ncdim_t *ncdims, int timedimid, stream_t *streamptr,
                    int *time_has_units, int *time_has_bounds, int *time_climatology)
{
  int ncvarid;

  if ( timedimid == UNDEFID )
    {
      char timeunits[CDI_MAX_NAME];

      for ( ncvarid = 0; ncvarid < nvars; ncvarid++ )
        {
          if ( ncvars[ncvarid].ndims == 0 && strcmp(ncvars[ncvarid].name, "time") == 0 )
            {
              if ( ncvars[ncvarid].units[0] )
                {
                  strcpy(timeunits, ncvars[ncvarid].units);
                  strtolower(timeunits);

                  if ( isTimeUnits(timeunits) )
                    {
                      streamptr->basetime.ncvarid = ncvarid;
                      break;
                    }
                }
            }
        }
    }
  else
    {
      int ltimevar = FALSE;

      if ( ncdims[timedimid].ncvarid != UNDEFID )
        {
          streamptr->basetime.ncvarid = ncdims[timedimid].ncvarid;
          ltimevar = TRUE;
        }

      for ( ncvarid = 0; ncvarid < nvars; ncvarid++ )
        if ( ncvarid != streamptr->basetime.ncvarid &&
             ncvars[ncvarid].ndims == 1 &&
             timedimid == ncvars[ncvarid].dimids[0] &&
             !xtypeIsText(ncvars[ncvarid].xtype) &&
             isTimeAxisUnits(ncvars[ncvarid].units) )
          {
            ncvars[ncvarid].isvar = FALSE;

            if ( !ltimevar )
              {
                streamptr->basetime.ncvarid = ncvarid;
                ltimevar = TRUE;
                if ( CDI_Debug )
                  fprintf(stderr, "timevar %s\n", ncvars[ncvarid].name);
              }
            else
              {
                Warning("Found more than one time variable, skipped variable %s!", ncvars[ncvarid].name);
              }
          }

      if ( ltimevar == FALSE ) /* search for WRF time description */
        {
          for ( ncvarid = 0; ncvarid < nvars; ncvarid++ )
            if ( ncvarid != streamptr->basetime.ncvarid &&
                 ncvars[ncvarid].ndims == 2 &&
                 timedimid == ncvars[ncvarid].dimids[0] &&
                 !xtypeIsText(ncvars[ncvarid].xtype) &&
                 ncdims[ncvars[ncvarid].dimids[1]].len == 19 )
              {
                streamptr->basetime.ncvarid = ncvarid;
                streamptr->basetime.lwrf    = TRUE;
                break;
              }
        }

      /* time varID */
      ncvarid = streamptr->basetime.ncvarid;

      if ( ncvarid == UNDEFID )
        {
          Warning("Time variable >%s< not found!", ncdims[timedimid].name);
        }
    }

  /* time varID */
  ncvarid = streamptr->basetime.ncvarid;

  if ( ncvarid != UNDEFID && streamptr->basetime.lwrf == FALSE )
    {
      if ( ncvars[ncvarid].units[0] != 0 ) *time_has_units = TRUE;

      if ( ncvars[ncvarid].bounds != UNDEFID )
        {
          int nbdims = ncvars[ncvars[ncvarid].bounds].ndims;
          if ( nbdims == 2 )
            {
              int len = (int) ncdims[ncvars[ncvars[ncvarid].bounds].dimids[nbdims-1]].len;
              if ( len == 2 && timedimid == ncvars[ncvars[ncvarid].bounds].dimids[0] )
                {
                  *time_has_bounds = TRUE;
                  streamptr->basetime.ncvarboundsid = ncvars[ncvarid].bounds;
                  if ( ncvars[ncvarid].climatology ) *time_climatology = TRUE;
                }
            }
        }
    }
}

int cdfInqContents(stream_t *streamptr)
{
  int ndims, nvars, ngatts, unlimdimid;
  int ncvarid;
  int ncdimid;
  size_t ntsteps;
  int timedimid = -1;
  int *varids;
  int nvarids;
  int time_has_units = FALSE;
  int time_has_bounds = FALSE;
  int time_climatology = FALSE;
  int leadtime_id = UNDEFID;
  int nvars_data;
  int nvcth_id = UNDEFID, vcta_id = UNDEFID, vctb_id = UNDEFID;
  size_t vctsize = 0;
  double *vct = NULL;
  int instID  = UNDEFID;
  int modelID = UNDEFID;
  int taxisID;
  int i;
  int calendar = UNDEFID;
  ncdim_t *ncdims;
  ncvar_t *ncvars = NULL;
  int format = 0;
  int ucla_les = FALSE;
  unsigned char uuidOfHGrid[CDI_UUID_SIZE];
  unsigned char uuidOfVGrid[CDI_UUID_SIZE];
  char gridfile[8912];
  char fcreftime[CDI_MAX_NAME];
  int number_of_grid_used = UNDEFID;

  memset(uuidOfHGrid, 0, CDI_UUID_SIZE);
  memset(uuidOfVGrid, 0, CDI_UUID_SIZE);
  gridfile[0] = 0;
  fcreftime[0] = 0;

  int vlistID = streamptr->vlistID;
  int fileID  = streamptr->fileID;

  if ( CDI_Debug ) Message("streamID = %d, fileID = %d", streamptr->self, fileID);

#if  defined  (HAVE_NETCDF4)
  nc_inq_format(fileID, &format);
#endif

  cdf_inq(fileID, &ndims , &nvars, &ngatts, &unlimdimid);

  if ( CDI_Debug )
    Message("root: ndims %d, nvars %d, ngatts %d", ndims, nvars, ngatts);

  if ( ndims == 0 )
    {
      Warning("ndims = %d", ndims);
      return (CDI_EUFSTRUCT);
    }

  /* alloc ncdims */
  ncdims = (ncdim_t *)xmalloc((size_t)ndims * sizeof (ncdim_t));
  init_ncdims(ndims, ncdims);

  if ( nvars > 0 )
    {
      /* alloc ncvars */
      ncvars = (ncvar_t *)xmalloc((size_t)nvars * sizeof (ncvar_t));
      init_ncvars(nvars, ncvars);

      for ( ncvarid = 0; ncvarid < nvars; ++ncvarid )
        ncvars[ncvarid].ncid = fileID;
    }

#if  defined  (TEST_GROUPS)
#if  defined  (HAVE_NETCDF4)
  if ( format == NC_FORMAT_NETCDF4 )
    {
      int ncid;
      int numgrps;
      int ncids[NC_MAX_VARS];
      char name1[CDI_MAX_NAME];
      int gndims, gnvars, gngatts, gunlimdimid;
      nc_inq_grps(fileID, &numgrps, ncids);
      for ( int i = 0; i < numgrps; ++i )
        {
          ncid = ncids[i];
          nc_inq_grpname (ncid, name1);
          cdf_inq(ncid, &gndims , &gnvars, &gngatts, &gunlimdimid);

          if ( CDI_Debug )
            Message("%s: ndims %d, nvars %d, ngatts %d", name1, gndims, gnvars, gngatts);

          if ( gndims == 0 )
            {
            }
        }
    }
#endif
#endif

  if ( nvars == 0 )
    {
      Warning("nvars = %d", nvars);
      return (CDI_EUFSTRUCT);
    }

  /* scan global attributes */
  scan_global_attributes(fileID, vlistID, streamptr, ngatts, &instID, &modelID, &ucla_les, 
                         uuidOfHGrid, uuidOfVGrid, gridfile, &number_of_grid_used);

  /* find time dim */
  if ( unlimdimid >= 0 )
    timedimid = unlimdimid;
  else
    timedimid = cdfTimeDimID(fileID, ndims, nvars);

  streamptr->basetime.ncdimid = timedimid;

  if ( timedimid != UNDEFID )
    cdf_inq_dimlen(fileID, timedimid, &ntsteps);
  else
    ntsteps = 0;

  if ( CDI_Debug ) Message("Number of timesteps = %d", ntsteps);
  if ( CDI_Debug ) Message("Time dimid = %d", streamptr->basetime.ncdimid);

  /* read ncdims */
  for ( ncdimid = 0; ncdimid < ndims; ncdimid++ )
    {
      cdf_inq_dimlen(fileID, ncdimid, &ncdims[ncdimid].len);
      cdf_inq_dimname(fileID, ncdimid, ncdims[ncdimid].name);
      if ( timedimid == ncdimid )
	ncdims[ncdimid].dimtype = T_AXIS;
    }

  if ( CDI_Debug ) printNCvars(ncvars, nvars, "cdfScanVarAttributes");

  /* scan attributes of all variables */
  cdfScanVarAttributes(nvars, ncvars, ncdims, timedimid, modelID, format);


  if ( CDI_Debug ) printNCvars(ncvars, nvars, "find coordinate vars");

  /* find coordinate vars */
  for ( ncdimid = 0; ncdimid < ndims; ncdimid++ )
    {
      for ( ncvarid = 0; ncvarid < nvars; ncvarid++ )
	{
	  if ( ncvars[ncvarid].ndims == 1 )
	    {
	      if ( timedimid != UNDEFID && timedimid == ncvars[ncvarid].dimids[0] )
		{
		  if ( ncvars[ncvarid].isvar != FALSE ) cdfSetVar(ncvars, ncvarid, TRUE);
		}
	      else
		{
                  //  if ( ncvars[ncvarid].isvar != TRUE ) cdfSetVar(ncvars, ncvarid, FALSE);
		}
	      // if ( ncvars[ncvarid].isvar != TRUE ) cdfSetVar(ncvars, ncvarid, FALSE);

	      if ( ncdimid == ncvars[ncvarid].dimids[0] && ncdims[ncdimid].ncvarid == UNDEFID )
		if ( strcmp(ncvars[ncvarid].name, ncdims[ncdimid].name) == 0 )
		  {
		    ncdims[ncdimid].ncvarid = ncvarid;
		    ncvars[ncvarid].isvar = FALSE;
		  }
	    }
	}
    }

  /* find time vars */
  find_time_vars(nvars, ncvars, ncdims, timedimid, streamptr, &time_has_units, &time_has_bounds, &time_climatology);

  leadtime_id = find_leadtime(nvars, ncvars);
  if ( leadtime_id != UNDEFID ) ncvars[leadtime_id].isvar = FALSE;

  /* check ncvars */
  for ( ncvarid = 0; ncvarid < nvars; ncvarid++ )
    {
      if ( timedimid != UNDEFID )
	if ( ncvars[ncvarid].isvar == -1 &&
	     ncvars[ncvarid].ndims > 1   &&
	     timedimid == ncvars[ncvarid].dimids[0] )
	  cdfSetVar(ncvars, ncvarid, TRUE);

      if ( ncvars[ncvarid].isvar == -1 && ncvars[ncvarid].ndims == 0 )
	cdfSetVar(ncvars, ncvarid, FALSE);

      //if ( ncvars[ncvarid].isvar == -1 && ncvars[ncvarid].ndims > 1 )
      if ( ncvars[ncvarid].isvar == -1 && ncvars[ncvarid].ndims >= 1 )
	cdfSetVar(ncvars, ncvarid, TRUE);

      if ( ncvars[ncvarid].isvar == -1 )
	{
	  ncvars[ncvarid].isvar = 0;
	  Warning("Variable %s has an unknown type, skipped!", ncvars[ncvarid].name);
	  continue;
	}

      if ( ncvars[ncvarid].ndims > 4 )
	{
	  ncvars[ncvarid].isvar = 0;
	  Warning("%d dimensional variables are not supported, skipped variable %s!",
		ncvars[ncvarid].ndims, ncvars[ncvarid].name);
	  continue;
	}

      if ( ncvars[ncvarid].ndims == 4 && timedimid == UNDEFID )
	{
	  ncvars[ncvarid].isvar = 0;
	  Warning("%d dimensional variables without time dimension are not supported, skipped variable %s!",
		ncvars[ncvarid].ndims, ncvars[ncvarid].name);
	  continue;
	}

      if ( xtypeIsText(ncvars[ncvarid].xtype) )
	{
	  ncvars[ncvarid].isvar = 0;
	  continue;
	}

      if ( cdfInqDatatype(ncvars[ncvarid].xtype, ncvars[ncvarid].lunsigned) == -1 )
	{
	  ncvars[ncvarid].isvar = 0;
	  Warning("Variable %s has an unsupported data type, skipped!", ncvars[ncvarid].name);
	  continue;
	}

      if ( timedimid != UNDEFID && ntsteps == 0 && ncvars[ncvarid].ndims > 0 )
	{
	  if ( timedimid == ncvars[ncvarid].dimids[0] )
	    {
	      ncvars[ncvarid].isvar = 0;
	      Warning("Number of time steps undefined, skipped variable %s!", ncvars[ncvarid].name);
	      continue;
	    }
	}
    }

  /* verify coordinate vars - first scan (dimname == varname) */
  verify_coordinate_vars_1(ndims, ncdims, ncvars, timedimid);

  /* verify coordinate vars - second scan (all other variables) */
  verify_coordinate_vars_2(nvars, ncvars);

  if ( CDI_Debug ) printNCvars(ncvars, nvars, "verify_coordinate_vars");

  if ( ucla_les == TRUE )
    {
      for ( ncdimid = 0; ncdimid < ndims; ncdimid++ )
	{
	  ncvarid = ncdims[ncdimid].ncvarid;
	  if ( ncvarid != -1 )
	    {
	      if ( ncdims[ncdimid].dimtype == UNDEFID && ncvars[ncvarid].units[0] == 'm' )
		{
		  if      ( ncvars[ncvarid].name[0] == 'x' ) ncdims[ncdimid].dimtype = X_AXIS;
		  else if ( ncvars[ncvarid].name[0] == 'y' ) ncdims[ncdimid].dimtype = Y_AXIS;
		  else if ( ncvars[ncvarid].name[0] == 'z' ) ncdims[ncdimid].dimtype = Z_AXIS;
		}
	    }
	}
    }
  /*
  for ( ncdimid = 0; ncdimid < ndims; ncdimid++ )
    {
      ncvarid = ncdims[ncdimid].ncvarid;
      if ( ncvarid != -1 )
	{
	  printf("coord var %d %s %s\n", ncvarid, ncvars[ncvarid].name, ncvars[ncvarid].units);
	  if ( ncdims[ncdimid].dimtype == X_AXIS )
	    printf("coord var %d %s is x dim\n", ncvarid, ncvars[ncvarid].name);
	  if ( ncdims[ncdimid].dimtype == Y_AXIS )
	    printf("coord var %d %s is y dim\n", ncvarid, ncvars[ncvarid].name);
	  if ( ncdims[ncdimid].dimtype == Z_AXIS )
	    printf("coord var %d %s is z dim\n", ncvarid, ncvars[ncvarid].name);
	  if ( ncdims[ncdimid].dimtype == T_AXIS )
	    printf("coord var %d %s is t dim\n", ncvarid, ncvars[ncvarid].name);

	  if ( ncvars[ncvarid].islon )
	    printf("coord var %d %s is lon\n", ncvarid, ncvars[ncvarid].name);
	  if ( ncvars[ncvarid].islat )
	    printf("coord var %d %s is lat\n", ncvarid, ncvars[ncvarid].name);
	  if ( ncvars[ncvarid].islev )
	    printf("coord var %d %s is lev\n", ncvarid, ncvars[ncvarid].name);
	}
    }
  */

  /* Set coordinate varids (att: associate)  */
  for ( ncvarid = 0; ncvarid < nvars; ncvarid++ )
    {
      if ( ncvars[ncvarid].isvar == TRUE && ncvars[ncvarid].ncoordvars )
	{
	  /* ndims = ncvars[ncvarid].ndims; */
	  ndims = ncvars[ncvarid].ncoordvars;
	  for ( i = 0; i < ndims; i++ )
	    {
	      if ( ncvars[ncvars[ncvarid].coordvarids[i]].islon )
		ncvars[ncvarid].xvarid = ncvars[ncvarid].coordvarids[i];
	      else if ( ncvars[ncvars[ncvarid].coordvarids[i]].islat )
		ncvars[ncvarid].yvarid = ncvars[ncvarid].coordvarids[i];
	      else if ( ncvars[ncvars[ncvarid].coordvarids[i]].islev )
		ncvars[ncvarid].zvarid = ncvars[ncvarid].coordvarids[i];
	    }
	}
    }

  /* set dim type */
  setDimType(nvars, ncvars, ncdims);

  /* find VCT */
  for ( ncvarid = 0; ncvarid < nvars; ncvarid++ )
    {
      if ( ncvars[ncvarid].ndims == 1 )
	{
	  if ( memcmp(ncvars[ncvarid].name, "hyai", 4) == 0 )
	    {
	      vcta_id = ncvarid;
	      nvcth_id = ncvars[ncvarid].dimids[0];
              ncvars[ncvarid].isvar = FALSE;
	      continue;
	    }
	  if ( memcmp(ncvars[ncvarid].name, "hybi", 4) == 0 )
	    {
	      vctb_id = ncvarid;
	      nvcth_id = ncvars[ncvarid].dimids[0];
              ncvars[ncvarid].isvar = FALSE;
	      continue;
	    }

	  if      ( memcmp(ncvars[ncvarid].name, "hyam", 4) == 0 ) ncvars[ncvarid].isvar = FALSE;
	  else if ( memcmp(ncvars[ncvarid].name, "hybm", 4) == 0 ) ncvars[ncvarid].isvar = FALSE;
	}
    }

  if ( CDI_Debug ) printNCvars(ncvars, nvars, "define_all_grids");


  /* define all grids */
  define_all_grids(streamptr, vlistID, ncdims, nvars, ncvars, timedimid, uuidOfHGrid, gridfile, number_of_grid_used);


  /* read VCT */
  if ( nvcth_id != UNDEFID && vcta_id != UNDEFID && vctb_id != UNDEFID )
    {
      vctsize = ncdims[nvcth_id].len;
      vctsize *= 2;
      vct = (double *) malloc(vctsize*sizeof(double));
      cdf_get_var_double(fileID, vcta_id, vct);
      cdf_get_var_double(fileID, vctb_id, vct+vctsize/2);
    }


  /* define all zaxes */
  define_all_zaxes(streamptr, vlistID, ncdims, nvars, ncvars, vctsize, vct, uuidOfVGrid);


  if ( vct ) free(vct);

  /* select vars */
  varids = (int *)xmalloc((size_t)nvars * sizeof (int));
  nvarids = 0;
  for ( ncvarid = 0; ncvarid < nvars; ncvarid++ )
    if ( ncvars[ncvarid].isvar == TRUE ) varids[nvarids++] = ncvarid;

  nvars_data = nvarids;

  if ( CDI_Debug ) Message("time varid = %d", streamptr->basetime.ncvarid);
  if ( CDI_Debug ) Message("ntsteps = %d", ntsteps);
  if ( CDI_Debug ) Message("nvars_data = %d", nvars_data);


  if ( nvars_data == 0 )
    {
      streamptr->ntsteps = 0;
      return (CDI_EUFSTRUCT);
    }

  if ( ntsteps == 0 && streamptr->basetime.ncdimid == UNDEFID && streamptr->basetime.ncvarid != UNDEFID )
    ntsteps = 1;

  streamptr->ntsteps = (long)ntsteps;

  /* define all data variables */
  define_all_vars(streamptr, vlistID, instID, modelID, varids, nvars_data, nvars, ncvars);


  cdiCreateTimesteps(streamptr);

  /* time varID */
  ncvarid = streamptr->basetime.ncvarid;

  if ( time_has_units )
    {
      taxis_t *taxis = &streamptr->tsteps[0].taxis;

      if ( setBaseTime(ncvars[ncvarid].units, taxis) == 1 )
        {
          ncvarid = UNDEFID;
          streamptr->basetime.ncvarid = UNDEFID;
        }

      if ( leadtime_id != UNDEFID && taxis->type == TAXIS_RELATIVE )
        {
          streamptr->basetime.leadtimeid = leadtime_id;
          taxis->type = TAXIS_FORECAST;

          int timeunit = -1;
          if ( ncvars[leadtime_id].units[0] != 0 ) timeunit = scanTimeUnit(ncvars[leadtime_id].units);
          if ( timeunit == -1 ) timeunit = taxis->unit;
          taxis->fc_unit = timeunit;

          setForecastTime(fcreftime, taxis);
        }
    }

  if ( time_has_bounds )
    {
      streamptr->tsteps[0].taxis.has_bounds = TRUE;
      if ( time_climatology ) streamptr->tsteps[0].taxis.climatology = TRUE;
    }

  if ( ncvarid != UNDEFID )
    {
      taxis_t *taxis = &streamptr->tsteps[0].taxis;
      ptaxisDefName(taxis, ncvars[ncvarid].name);
      if ( ncvars[ncvarid].longname[0] )
        ptaxisDefLongname(taxis, ncvars[ncvarid].longname);
    }

  if ( ncvarid != UNDEFID )
    if ( ncvars[ncvarid].calendar == TRUE )
      {
        enum {
          attstringlen = 8192,
        };
        char attstring[attstringlen];

	cdfGetAttText(fileID, ncvarid, "calendar", attstringlen, attstring);
	strtolower(attstring);

	if ( memcmp(attstring, "standard", 8)  == 0 ||
	     memcmp(attstring, "gregorian", 9) == 0 )
	  calendar = CALENDAR_STANDARD;
	else if ( memcmp(attstring, "none", 4) == 0 )
	  calendar = CALENDAR_NONE;
	else if ( memcmp(attstring, "proleptic", 9) == 0 )
	  calendar = CALENDAR_PROLEPTIC;
	else if ( memcmp(attstring, "360", 3) == 0 )
	  calendar = CALENDAR_360DAYS;
	else if ( memcmp(attstring, "365", 3) == 0 ||
		  memcmp(attstring, "noleap", 6)  == 0 )
	  calendar = CALENDAR_365DAYS;
	else if ( memcmp(attstring, "366", 3)  == 0 ||
		  memcmp(attstring, "all_leap", 8) == 0 )
	  calendar = CALENDAR_366DAYS;
	else
	  Warning("calendar >%s< unsupported!", attstring);
      }

  if ( streamptr->tsteps[0].taxis.type == TAXIS_FORECAST )
    {
      taxisID = taxisCreate(TAXIS_FORECAST);
    }
  else if ( streamptr->tsteps[0].taxis.type == TAXIS_RELATIVE )
    {
      taxisID = taxisCreate(TAXIS_RELATIVE);
    }
  else
    {
      taxisID = taxisCreate(TAXIS_ABSOLUTE);
      if ( !time_has_units )
	{
	  taxisDefTunit(taxisID, TUNIT_DAY);
	  streamptr->tsteps[0].taxis.unit = TUNIT_DAY;
	}
    }


  if ( calendar == UNDEFID && streamptr->tsteps[0].taxis.type != TAXIS_ABSOLUTE )
    {
      calendar = CALENDAR_STANDARD;
    }

#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ > 5)
#pragma GCC diagnostic push
#pragma GCC diagnostic warning "-Wstrict-overflow"
#endif
  if ( calendar != UNDEFID )
    {
      taxis_t *taxis = &streamptr->tsteps[0].taxis;
      taxis->calendar = calendar;
      taxisDefCalendar(taxisID, calendar);
    }
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ > 5)
#pragma GCC diagnostic pop
#endif

  vlistDefTaxis(vlistID, taxisID);

  streamptr->curTsID = 0;
  streamptr->rtsteps = 1;

  (void) cdfInqTimestep(streamptr, 0);

  cdfCreateRecords(streamptr, 0);

  /* free ncdims */
  free(ncdims);

  /* free ncvars */
  free(ncvars);

  return (0);
}

static
void wrf_read_timestep(int fileID, int nctimevarid, int tsID, taxis_t *taxis)
{
  size_t start[2], count[2];
  char stvalue[32];
  start[0] = (size_t) tsID; start[1] = 0;
  count[0] = 1; count[1] = 19;
  stvalue[0] = 0;
  cdf_get_vara_text(fileID, nctimevarid, start, count, stvalue);
  stvalue[19] = 0;
  {
    int year = 1, month = 1, day = 1 , hour = 0, minute = 0, second = 0;
    if ( strlen(stvalue) == 19 )
      sscanf(stvalue, "%d-%d-%d_%d:%d:%d", &year, &month, &day, &hour, &minute, &second);
    taxis->vdate = cdiEncodeDate(year, month, day);
    taxis->vtime = cdiEncodeTime(hour, minute, second);
    taxis->type = TAXIS_ABSOLUTE;
  }
}

static
double get_timevalue(int fileID, int nctimevarid, int tsID, timecache_t *tcache)
{
  double timevalue = 0;
  size_t index = (size_t) tsID;

  if ( tcache )
    {
      if ( tcache->size == 0 || (tsID < tcache->startid || tsID > (tcache->startid+tcache->size-1)) )
        {
          int maxvals = MAX_TIMECACHE_SIZE;
          tcache->startid = (tsID/MAX_TIMECACHE_SIZE)*MAX_TIMECACHE_SIZE;
          if ( (tcache->startid + maxvals) > tcache->maxvals ) maxvals = (tcache->maxvals)%MAX_TIMECACHE_SIZE;
          tcache->size = maxvals;
          index = (size_t) tcache->startid;
          // fprintf(stderr, "fill time cache: %d %d %d %d %d\n", tcache->maxvals, tsID, tcache->startid, tcache->startid+maxvals-1, maxvals);
          for ( int ival = 0; ival < maxvals; ++ival )
            {
              cdf_get_var1_double(fileID, nctimevarid, &index, &timevalue);
              if ( timevalue >= NC_FILL_DOUBLE || timevalue < -NC_FILL_DOUBLE ) timevalue = 0;
              tcache->cache[ival] = timevalue;
              index++;
            }
        }

      timevalue = tcache->cache[tsID%MAX_TIMECACHE_SIZE];
    }
  else
    {
      cdf_get_var1_double(fileID, nctimevarid, &index, &timevalue);
      if ( timevalue >= NC_FILL_DOUBLE || timevalue < -NC_FILL_DOUBLE ) timevalue = 0;
    }

  return timevalue;
}

int cdfInqTimestep(stream_t * streamptr, int tsID)
{
  long nrecs = 0;
  double timevalue;
  int fileID;
  taxis_t *taxis;

  if ( CDI_Debug ) Message("streamID = %d  tsID = %d", streamptr->self, tsID);

  if ( tsID < 0 ) Error("unexpected tsID = %d", tsID);

  if ( tsID < streamptr->ntsteps && streamptr->ntsteps > 0 )
    {
      cdfCreateRecords(streamptr, tsID);

      taxis = &streamptr->tsteps[tsID].taxis;
      if ( tsID > 0 )
	ptaxisCopy(taxis, &streamptr->tsteps[0].taxis);

      timevalue = tsID;

      int nctimevarid = streamptr->basetime.ncvarid;
      if ( nctimevarid != UNDEFID )
	{
	  fileID = streamptr->fileID;
	  size_t index  = (size_t)tsID;

	  if ( streamptr->basetime.lwrf )
	    {
              wrf_read_timestep(fileID, nctimevarid, tsID, taxis);
	    }
	  else
	    {
#if defined (USE_TIMECACHE)
              if ( streamptr->basetime.timevar_cache == NULL )
                {
                  streamptr->basetime.timevar_cache = (timecache_t *) malloc(MAX_TIMECACHE_SIZE*sizeof(timecache_t));
                  streamptr->basetime.timevar_cache->size = 0;
                  streamptr->basetime.timevar_cache->maxvals = streamptr->ntsteps;
                }
#endif
              timevalue = get_timevalue(fileID, nctimevarid, tsID, streamptr->basetime.timevar_cache);
	      cdiDecodeTimeval(timevalue, taxis, &taxis->vdate, &taxis->vtime);
	    }

	  int nctimeboundsid = streamptr->basetime.ncvarboundsid;
	  if ( nctimeboundsid != UNDEFID )
	    {
	      size_t start[2], count[2];
              start[0] = index; count[0] = 1; start[1] = 0; count[1] = 1;
	      cdf_get_vara_double(fileID, nctimeboundsid, start, count, &timevalue);
              if ( timevalue >= NC_FILL_DOUBLE || timevalue < -NC_FILL_DOUBLE ) timevalue = 0;

	      cdiDecodeTimeval(timevalue, taxis, &taxis->vdate_lb, &taxis->vtime_lb);

              start[0] = index; count[0] = 1; start[1] = 1; count[1] = 1;
	      cdf_get_vara_double(fileID, nctimeboundsid, start, count, &timevalue);
              if ( timevalue >= NC_FILL_DOUBLE || timevalue < -NC_FILL_DOUBLE ) timevalue = 0;

	      cdiDecodeTimeval(timevalue, taxis, &taxis->vdate_ub, &taxis->vtime_ub);
	    }

          int leadtimeid = streamptr->basetime.leadtimeid;
          if ( leadtimeid != UNDEFID )
            {
              timevalue = get_timevalue(fileID, leadtimeid, tsID, NULL);
              cdiSetForecastPeriod(timevalue, taxis);
            }
	}
    }

  streamptr->curTsID = tsID;
  nrecs = streamptr->tsteps[tsID].nrecs;

  return ((int) nrecs);
}


void cdfEndDef(stream_t *streamptr)
{
  int varID;
  int nvars;
  int fileID;

  fileID  = streamptr->fileID;

  cdfDefGlobalAtts(streamptr);
  cdfDefLocalAtts(streamptr);

  if ( streamptr->accessmode == 0 )
    {
      nvars =  streamptr->nvars;

      if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

      for ( varID = 0; varID < nvars; varID++ )
	cdfDefVar(streamptr, varID);

      if ( streamptr->ncmode == 2 )
        {
          extern size_t CDI_netcdf_hdr_pad;

          if ( CDI_netcdf_hdr_pad == 0UL )
            cdf_enddef(fileID);
          else
            cdf__enddef(fileID, CDI_netcdf_hdr_pad);
        }

      streamptr->accessmode = 1;
    }
}


static void cdfDefInstitut(stream_t *streamptr)
{
  int fileID, instID;
  size_t len;
  int vlistID;

  vlistID = streamptr->vlistID;
  fileID  = streamptr->fileID;
  instID  = vlistInqInstitut(vlistID);

  if ( instID != UNDEFID )
    {
      const char *longname = institutInqLongnamePtr(instID);
      if ( longname )
	{
	  len = strlen(longname);
	  if ( len > 0 )
	    {
	      if ( streamptr->ncmode == 2 ) cdf_redef(fileID);
	      cdf_put_att_text(fileID, NC_GLOBAL, "institution", len, longname);
	      if ( streamptr->ncmode == 2 ) cdf_enddef(fileID);
	    }
	}
    }
}


static void cdfDefSource(stream_t *streamptr)
{
  int vlistID = streamptr->vlistID;
  int fileID  = streamptr->fileID;
  int modelID = vlistInqModel(vlistID);

  if ( modelID != UNDEFID )
    {
      const char *longname = modelInqNamePtr(modelID);
      if ( longname )
	{
          size_t len = strlen(longname);
	  if ( len > 0 )
	    {
	      if ( streamptr->ncmode == 2 ) cdf_redef(fileID);
	      cdf_put_att_text(fileID, NC_GLOBAL, "source", len, longname);
	      if ( streamptr->ncmode == 2 ) cdf_enddef(fileID);
	    }
	}
    }
}


static void cdfDefGlobalAtts(stream_t *streamptr)
{
  int natts;

  if ( streamptr->globalatts ) return;

  int vlistID = streamptr->vlistID;
  int fileID  = streamptr->fileID;

  cdfDefSource(streamptr);
  cdfDefInstitut(streamptr);

  vlistInqNatts(vlistID, CDI_GLOBAL, &natts);

  if ( natts > 0 && streamptr->ncmode == 2 ) cdf_redef(fileID);

  defineAttributes(vlistID, CDI_GLOBAL, fileID, NC_GLOBAL);

  if ( natts > 0 && streamptr->ncmode == 2 ) cdf_enddef(fileID);

  streamptr->globalatts = 1;
}


static void cdfDefLocalAtts(stream_t *streamptr)
{
  int vlistID = streamptr->vlistID;
  int fileID  = streamptr->fileID;

  if ( streamptr->localatts ) return;
  if ( vlistInqInstitut(vlistID) != UNDEFID ) return;

  streamptr->localatts = 1;

  if ( streamptr->ncmode == 2 ) cdf_redef(fileID);

  for ( int varID = 0; varID < streamptr->nvars; varID++ )
    {
      int instID = vlistInqVarInstitut(vlistID, varID);
      if ( instID != UNDEFID )
	{
          int ncvarid = streamptr->vars[varID].ncvarid;
  	  const char *name = institutInqNamePtr(instID);
	  if ( name )
	    {
              size_t len = strlen(name);
	      cdf_put_att_text(fileID, ncvarid, "institution", len, name);
	    }
	}
      }

  if ( streamptr->ncmode == 2 ) cdf_enddef(fileID);
}


void cdfDefHistory(stream_t *streamptr, int size, const char *history)
{
  int ncid = streamptr->fileID;
  cdf_put_att_text(ncid, NC_GLOBAL, "history", (size_t) size, history);
}


int cdfInqHistorySize(stream_t *streamptr)
{
  size_t size = 0;
  int ncid = streamptr->fileID;
  if ( streamptr->historyID != UNDEFID )
    cdf_inq_attlen(ncid, NC_GLOBAL, "history", &size);

  return ((int) size);
}


void cdfInqHistoryString(stream_t *streamptr, char *history)
{
  int ncid = streamptr->fileID;
  if ( streamptr->historyID != UNDEFID )
    cdf_get_att_text(ncid, NC_GLOBAL, "history", history);
}


void cdfDefVars(stream_t *streamptr)
{
  /* int  nvars, ncvarid; */

  int vlistID = streamptr->vlistID;
  if ( vlistID == UNDEFID )
    Error("Internal problem! vlist undefined for streamptr %p", streamptr);

  /* nvars  = vlistNvars(vlistID); */
  int ngrids = vlistNgrids(vlistID);
  int nzaxis = vlistNzaxis(vlistID);
  /*
  if ( vlistHasTime(vlistID) ) cdfDefTime(streamptr);
  */
  if (ngrids > 0)
    for (size_t index = 0; index < (size_t)ngrids; index++ )
      {
        int gridID = vlistGrid(vlistID, (int)index);
        cdfDefGrid(streamptr, gridID);
      }

  if (nzaxis > 0)
    for (size_t index = 0; index < (size_t)nzaxis; index++ )
      {
        int zaxisID = vlistZaxis(vlistID, (int)index);
        if ( streamptr->zaxisID[index] == UNDEFID ) cdfDefZaxis(streamptr, zaxisID);
      }
  /*
    define time first!!!
  for (varID = 0; varID < nvars; varID++ )
    {
      ncvarid = cdfDefVar(streamptr, varID);
    }
  */
}
#endif
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
