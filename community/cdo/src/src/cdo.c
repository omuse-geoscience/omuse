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

#ifndef _XOPEN_SOURCE
//#define _XOPEN_SOURCE 600 /* gethostname */
#endif

#include <ctype.h>
/*#include <malloc.h>*/ /* mallopt and malloc_stats */
#include <sys/stat.h>
#if defined(HAVE_GETRLIMIT)
#if defined(HAVE_SYS_RESOURCE_H)
#include <sys/time.h>       /* getrlimit */
#include <sys/resource.h>   /* getrlimit */
#endif
#endif
#include <unistd.h>         /* sysconf, gethostname */

#if defined(SX)
#define RLIM_T  long long
#else
#define RLIM_T  rlim_t
#endif

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"

#include "cdo_getopt.h"

#if defined(HAVE_LIBPTHREAD)
#include "pstream_int.h"
#include "pthread_debug.h"
#endif

#include "modules.h"
#include "util.h"
#include "error.h"

#if defined(_OPENMP)
#  include <omp.h>
#endif

#if ! defined(VERSION)
#  define  VERSION  "0.0.1"
#endif

#define MAX_NUM_VARNAMES 256

static int Debug = 0;
static int Version = 0;
static int Help = 0;
static int DebugLevel = 0;
static int numThreads = 0;
static int timer_total;
static int CDO_netcdf_hdr_pad = 0;
static int CDO_Rusage = 0;


#define PRINT_RLIMIT(resource) \
      { \
        int status; \
        struct rlimit rlim; \
        status = getrlimit(resource, &rlim); \
        if ( status == 0 ) \
          { \
            if ( sizeof(RLIM_T) > sizeof(long) ) \
              { \
                fprintf(stderr, "CUR %-15s = %llu\n", #resource, (long long) rlim.rlim_cur); \
                fprintf(stderr, "MAX %-15s = %llu\n", #resource, (long long) rlim.rlim_max); \
              } \
            else \
              { \
                fprintf(stderr, "CUR %-15s = %lu\n", #resource, (long) rlim.rlim_cur); \
                fprintf(stderr, "MAX %-15s = %lu\n", #resource, (long) rlim.rlim_max); \
              } \
          } \
      }


static
void cdo_version(void)
{
  const int   filetypes[] = {FILETYPE_SRV, FILETYPE_EXT, FILETYPE_IEG, FILETYPE_GRB, FILETYPE_GRB2, FILETYPE_NC, FILETYPE_NC2, FILETYPE_NC4, FILETYPE_NC4C};
  const char* typenames[] = {        "srv",        "ext",        "ieg",        "grb",        "grb2",        "nc",        "nc2",        "nc4",        "nc4c"};

  fprintf(stderr, "%s\n", CDO_Version);
#if defined(USER_NAME) && defined(HOST_NAME) && defined(SYSTEM_TYPE)
  fprintf(stderr, "Compiled: by %s on %s (%s) %s %s\n", USER_NAME, HOST_NAME, SYSTEM_TYPE, __DATE__, __TIME__);
#endif
#if defined(COMPILER)
  fprintf(stderr, "Compiler: %s\n", COMPILER);
#endif
#if defined(COMP_VERSION)
  fprintf(stderr, " version: %s\n", COMP_VERSION);
#endif

  printFeatures();
  printLibraries();

  fprintf(stderr, "Filetypes: ");
  set_text_color(stderr, BRIGHT, GREEN);
  for ( size_t i = 0; i < sizeof(filetypes)/sizeof(int); ++i )
    if ( cdiHaveFiletype(filetypes[i]) ) fprintf(stderr, "%s ", typenames[i]);
  reset_text_color(stderr);
  fprintf(stderr, "\n");

  cdiPrintVersion();
  fprintf(stderr, "\n");
}

static
void cdo_usage(void)
{
  const char *name;

  /*  fprintf(stderr, "%s\n", CDO_Version);*/
  /*  fprintf(stderr, "\n");*/
  fprintf(stderr, "usage : cdo  [Options]  Operator1  [-Operator2  [-OperatorN]]\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  Options:\n");
  set_text_color(stderr, RESET, BLUE);
  fprintf(stderr, "    -a             Generate an absolute time axis\n");
  fprintf(stderr, "    -b <nbits>     Set the number of bits for the output precision\n");
  fprintf(stderr, "                   (I8/I16/I32/F32/F64 for nc/nc2/nc4/nc4c; F32/F64 for grb2/srv/ext/ieg; P1 - P24 for grb/grb2)\n");
  fprintf(stderr, "                   Add L or B to set the byteorder to Little or Big endian\n");
  fprintf(stderr, "    -f, --format <format>\n");
  fprintf(stderr, "                   Format of the output file. (grb/grb2/nc/nc2/nc4/nc4c/srv/ext/ieg)\n");
  fprintf(stderr, "    -g <grid>      Set default grid name or file. Available grids: \n");
  fprintf(stderr, "                   n<N>, t<RES>, tl<RES>, global_<DXY>, r<NX>x<NY>, g<NX>x<NY>, gme<NI>, lon=<LON>/lat=<LAT>\n");
  fprintf(stderr, "    -h, --help     Help information for the operators\n");
  fprintf(stderr, "    --history      Do not append to netCDF \"history\" global attribute\n");
  fprintf(stderr, "    --netcdf_hdr_pad, --hdr_pad, --header_pad <nbr>\n");
  fprintf(stderr, "                   Pad netCDF output header with nbr bytes\n");
  /*
  fprintf(stderr, "    -i <inst>      Institution name/file\n");
  fprintf(stderr, "                   Predefined instituts: ");
  for ( int id = 0; id < institutInqNumber; id++ )
    if ( (name = institutInqNamePtr(id)) )
      fprintf(stderr, " %s", name);
  fprintf(stderr, "\n");
  */
  /* fprintf(stderr, "    -l <level>     Level file\n"); */
  fprintf(stderr, "    -k <chunktype> NetCDF4 chunk type: auto, grid or lines\n");
  fprintf(stderr, "    -L             Lock IO (sequential access)\n");
  fprintf(stderr, "    -M             Switch to indicate that the I/O streams have missing values\n");
  fprintf(stderr, "    -m <missval>   Set the default missing value (default: %g)\n", cdiInqMissval());
  fprintf(stderr, "    --no_warnings  Inhibit warning messages\n");
  fprintf(stderr, "    -O             Overwrite existing output file, if checked\n");
#if defined(_OPENMP)
  fprintf(stderr, "    -P <nthreads>  Set number of OpenMP threads\n");
#endif
  fprintf(stderr, "    -Q             Alphanumeric sorting of netCDF parameter names\n");
  fprintf(stderr, "    --reduce_dim   Reduce netCDF dimensions (module: TIMSTAT, FLDSTAT)\n");
  fprintf(stderr, "    -R, --regular  Convert GRIB1 data from reduced to regular grid (only with cgribex)\n");
  fprintf(stderr, "    -r             Generate a relative time axis\n");
  fprintf(stderr, "    -S             Create an extra output stream for the module TIMSTAT. This stream\n");
  fprintf(stderr, "                   contains the number of non missing values for each output period.\n");
  fprintf(stderr, "    -s, --silent   Silent mode\n");
  fprintf(stderr, "    -t <partab>    Set default parameter table name or file\n");
  fprintf(stderr, "                   Predefined tables: ");
  for ( int id = 0; id < tableInqNumber(); id++ )
    if ( (name = tableInqNamePtr(id)) )
      fprintf(stderr, " %s", name);
  fprintf(stderr, "\n");

  fprintf(stderr, "    -V, --version  Print the version number\n");
  fprintf(stderr, "    -v, --verbose  Print extra details for some operators\n");
  fprintf(stderr, "    -W             Print extra warning messages\n");
  fprintf(stderr, "    -z szip        SZIP compression of GRIB1 records\n");
  fprintf(stderr, "       jpeg        JPEG compression of GRIB2 records\n");
  fprintf(stderr, "        zip[_1-9]  Deflate compression of netCDF4 variables\n");
  reset_text_color(stderr);
  fprintf(stderr, "\n");

  fprintf(stderr, "  Operators:\n");
  set_text_color(stderr, RESET, GREEN);
  operatorPrintAll();
  reset_text_color(stderr);

  fprintf(stderr, "\n");
  fprintf(stderr, "  CDO version %s, Copyright (C) 2003-2015 Uwe Schulzweida\n", VERSION);
  //  fprintf(stderr, "  Available from <http://mpimet.mpg.de/cdo>\n");
  fprintf(stderr, "  This is free software and comes with ABSOLUTELY NO WARRANTY\n");
  fprintf(stderr, "  Report bugs to <http://mpimet.mpg.de/cdo>\n");
}

static
void cdo_init_is_tty(void)
{
  struct stat statbuf;
  fstat(0, &statbuf);
  if ( S_ISCHR(statbuf.st_mode) ) stdin_is_tty = 1;  
  fstat(1, &statbuf);
  if ( S_ISCHR(statbuf.st_mode) ) stdout_is_tty = 1;  
  fstat(2, &statbuf);
  if ( S_ISCHR(statbuf.st_mode) ) stderr_is_tty = 1;  
}

static
void cdoPrintHelp(char *phelp[]/*, char *xoperator*/)
{
  if ( phelp == NULL )
    fprintf(stderr, "No help available for this operator!\n");
  else
    {
      int lprint;
      while ( *phelp )
        {
          lprint = TRUE;
          if ( *phelp[0] == '\0' )
            if ( *(phelp+1) )
              if ( *(phelp+1)[0] == ' ' ) lprint = FALSE;
          
          if ( lprint )
            {
              if ( COLOR_STDOUT )
                {
                  if ( (strcmp(*phelp, "NAME")        == 0) ||
                       (strcmp(*phelp, "SYNOPSIS")    == 0) ||
                       (strcmp(*phelp, "DESCRIPTION") == 0) ||
                       (strcmp(*phelp, "OPERATORS")   == 0) ||
                       (strcmp(*phelp, "NAMELIST")    == 0) ||
                       (strcmp(*phelp, "PARAMETER")   == 0) ||
                       (strcmp(*phelp, "ENVIRONMENT") == 0) ||
                       (strcmp(*phelp, "EXAMPLES")    == 0) )
                    {
                      set_text_color(stdout, BRIGHT, BLACK);
                      fprintf(stdout, "%s", *phelp);
                      reset_text_color(stdout);
                      fprintf(stdout, "\n");
                    }
                  else
                    fprintf(stdout, "%s\n", *phelp);
                }
              else
                {
                  fprintf(stdout, "%s\n", *phelp);
                }
            }

          phelp++;
        }
    }
}


static
void cdoSetDebug(int level)
{
  /*
    level   0: off
    level   1: on
    level   2: cdi
    level   4: memory
    level   8: file
    level  16: format
    level  32: cdo
    level  64: stream
    level 128: pipe
    level 256: pthread
   */
  cdiDebug(level);

  if ( level == 1 || (level &  32) ) cdoDebug = 1;
  if ( level == 1 || (level &  64) ) pstreamDebug(1);
#if defined(HAVE_LIBPTHREAD)
  if ( level == 1 || (level & 128) ) pipeDebug(1);
  if ( level == 1 || (level & 256) ) Pthread_debug(1);
#endif
}

#undef  IsBigendian
#define IsBigendian()  ( u_byteorder.c[sizeof(long) - 1] )

static
void setDefaultDataType(char *datatypestr)
{
  static union {unsigned long l; unsigned char c[sizeof(long)];} u_byteorder = {1};
  int nbits = -1;
  enum {D_UINT, D_INT, D_FLT, D_CPX};
  int dtype = -1;

  if      ( *datatypestr == 'i' || *datatypestr == 'I' )
    {
      dtype = D_INT;
      datatypestr++;
    }
  else if ( *datatypestr == 'u' || *datatypestr == 'U' )
    {
      dtype = D_UINT;
      datatypestr++;
    }
  else if ( *datatypestr == 'f' || *datatypestr == 'F' )
    {
      dtype = D_FLT;
      datatypestr++;
    }
  else if ( *datatypestr == 'c' || *datatypestr == 'C' )
    {
      dtype = D_CPX;
      datatypestr++;
    }
  else if ( *datatypestr == 'p' || *datatypestr == 'P' )
    {
      datatypestr++;
    }

  if ( isdigit((int) *datatypestr) )
    {
      nbits = atoi(datatypestr);
      if ( nbits < 10 )
        datatypestr += 1;
      else
        datatypestr += 2;

      if ( dtype == -1 )
        {
          if      ( nbits > 0 && nbits < 32 ) cdoDefaultDataType = nbits;
          else if ( nbits == 32 )
            {
              if ( cdoDefaultFileType == FILETYPE_GRB )
                cdoDefaultDataType = DATATYPE_PACK32;
              else
                cdoDefaultDataType = DATATYPE_FLT32;
            }
          else if ( nbits == 64 ) cdoDefaultDataType = DATATYPE_FLT64;
          else
            {
              fprintf(stderr, "Unsupported number of bits %d!\n", nbits);
              fprintf(stderr, "Use I8/I16/I32/F32/F64 for nc/nc2/nc4/nc4c; F32/F64 for grb2/srv/ext/ieg; P1 - P24 for grb/grb2.\n");
              exit(EXIT_FAILURE);
            }
        }
      else
        {
          if ( dtype == D_INT )
            {
              if      ( nbits ==  8 ) cdoDefaultDataType = DATATYPE_INT8;
              else if ( nbits == 16 ) cdoDefaultDataType = DATATYPE_INT16;
              else if ( nbits == 32 ) cdoDefaultDataType = DATATYPE_INT32;
              else
                {
                  fprintf(stderr, "Unsupported number of bits = %d for datatype INT!\n", nbits);
                  exit(EXIT_FAILURE);
                }
            }
          else if ( dtype == D_UINT )
            {
              if      ( nbits ==  8 ) cdoDefaultDataType = DATATYPE_UINT8;
              else if ( nbits == 16 ) cdoDefaultDataType = DATATYPE_UINT16;
              else if ( nbits == 32 ) cdoDefaultDataType = DATATYPE_UINT32;
              else
                {
                  fprintf(stderr, "Unsupported number of bits = %d for datatype UINT!\n", nbits);
                  exit(EXIT_FAILURE);
                }
            }
          else if ( dtype == D_FLT )
            {
              if      ( nbits == 32 ) cdoDefaultDataType = DATATYPE_FLT32;
              else if ( nbits == 64 ) cdoDefaultDataType = DATATYPE_FLT64;
              else
                {
                  fprintf(stderr, "Unsupported number of bits = %d for datatype FLT!\n", nbits);
                  exit(EXIT_FAILURE);
                }
            }
          else if ( dtype == D_CPX )
            {
              if      ( nbits == 32 ) cdoDefaultDataType = DATATYPE_CPX32;
              else if ( nbits == 64 ) cdoDefaultDataType = DATATYPE_CPX64;
              else
                {
                  fprintf(stderr, "Unsupported number of bits = %d for datatype CPX!\n", nbits);
                  exit(EXIT_FAILURE);
                }
            }
        }
    }

  if ( *datatypestr != 0 )
    {
      if ( *datatypestr == 'l' || *datatypestr == 'L' )
        {
          if ( IsBigendian() ) cdoDefaultByteorder = CDI_LITTLEENDIAN;
          datatypestr++;
        }
      else if ( *datatypestr == 'b' || *datatypestr == 'B' )
        {
          if ( ! IsBigendian() ) cdoDefaultByteorder = CDI_BIGENDIAN;
          datatypestr++;
        }
      else
        {
          fprintf(stderr, "Unsupported character in number of bytes: >%s< !\n", datatypestr);
          exit(EXIT_FAILURE);
        }
    }
}

static
void setDefaultDataTypeByte(char *datatypestr)
{
  static union {unsigned long l; unsigned char c[sizeof(long)];} u_byteorder = {1};
  int datatype = -1;

  if ( isdigit((int) *datatypestr) )
    {
      datatype = atoi(datatypestr);
      datatypestr++;

      if      ( datatype == 1 ) cdoDefaultDataType = DATATYPE_PACK8;
      else if ( datatype == 2 ) cdoDefaultDataType = DATATYPE_PACK16;
      else if ( datatype == 3 ) cdoDefaultDataType = DATATYPE_PACK24;
      else if ( datatype == 4 ) cdoDefaultDataType = DATATYPE_FLT32;
      else if ( datatype == 8 ) cdoDefaultDataType = DATATYPE_FLT64;
      else
        {
          fprintf(stderr, "Unsupported datatype %d!\n", datatype);
          fprintf(stderr, "Use 4/8 for filetype nc/srv/ext/ieg and 1/2/3 for grb/grb2.\n");
          exit(EXIT_FAILURE);
        }
    }

  if ( *datatypestr != 0 )
    {
      if ( *datatypestr == 'l' || *datatypestr == 'L' )
        {
          if ( IsBigendian() ) cdoDefaultByteorder = CDI_LITTLEENDIAN;
          datatypestr++;
        }
      else if ( *datatypestr == 'b' || *datatypestr == 'B' )
        {
          if ( ! IsBigendian() ) cdoDefaultByteorder = CDI_BIGENDIAN;
          datatypestr++;
        }
      else
        {
          fprintf(stderr, "Unsupported character in number of bytes: %s!\n", datatypestr);
          exit(EXIT_FAILURE);
        }
    }
}

static
void setDefaultFileType(char *filetypestr, int labort)
{
  if ( filetypestr )
    {
      char *ftstr = filetypestr;
      size_t len;

      if      ( cmpstrlen(filetypestr, "grb2", len)  == 0 ) { ftstr += len; cdoDefaultFileType = FILETYPE_GRB2;}
      else if ( cmpstrlen(filetypestr, "grb1", len)  == 0 ) { ftstr += len; cdoDefaultFileType = FILETYPE_GRB; }
      else if ( cmpstrlen(filetypestr, "grb",  len)  == 0 ) { ftstr += len; cdoDefaultFileType = FILETYPE_GRB; }
      else if ( cmpstrlen(filetypestr, "nc2",  len)  == 0 ) { ftstr += len; cdoDefaultFileType = FILETYPE_NC2; }
      else if ( cmpstrlen(filetypestr, "nc4c", len)  == 0 ) { ftstr += len; cdoDefaultFileType = FILETYPE_NC4C;}
      else if ( cmpstrlen(filetypestr, "nc4",  len)  == 0 ) { ftstr += len; cdoDefaultFileType = FILETYPE_NC4; }
      else if ( cmpstrlen(filetypestr, "nc",   len)  == 0 ) { ftstr += len; cdoDefaultFileType = FILETYPE_NC;  }
      else if ( cmpstrlen(filetypestr, "srv",  len)  == 0 ) { ftstr += len; cdoDefaultFileType = FILETYPE_SRV; }
      else if ( cmpstrlen(filetypestr, "ext",  len)  == 0 ) { ftstr += len; cdoDefaultFileType = FILETYPE_EXT; }
      else if ( cmpstrlen(filetypestr, "ieg",  len)  == 0 ) { ftstr += len; cdoDefaultFileType = FILETYPE_IEG; }
      else
        {
          if ( labort )
            {
              fprintf(stderr, "Unsupported filetype %s!\n", filetypestr);
              fprintf(stderr, "Available filetypes: grb/grb2/nc/nc2/nc4/nc4c/srv/ext/ieg\n");
              exit(EXIT_FAILURE);
            }
          else
            {
              return;
            }
        }

      if ( cdoDefaultFileType != CDI_UNDEFID && *ftstr != 0 )
        {
          if ( *ftstr == '_' )
            {
              ftstr++;

              setDefaultDataType(ftstr);
            }
          else
            {
              fprintf(stderr, "Unexpected character >%c< in file type >%s<!\n", *ftstr, filetypestr);
              fprintf(stderr, "Use format[_nbits] with:\n");
              fprintf(stderr, "    format = grb, grb2, nc, nc2, nc4, nc4c, srv, ext or ieg\n");
              fprintf(stderr, "    nbits  = 32/64 for grb2/nc/nc2/nc4/nc4c/srv/ext/ieg; 1 - 24 for grb/grb2\n");
              exit(EXIT_FAILURE);
            }
        }
    }
}

#if defined(malloc)
#undef malloc
#undef free
#endif

#define NTESTS 11
#include <inttypes.h>
static
int getMemAlignment(void)
{
  int ma = -1;
  int i, k;
  double *ptr[NTESTS];
  int64_t iptr;
  size_t tsize[NTESTS] = {1, 3, 5, 9, 17, 33, 69, 121, 251, 510, 1025};
  size_t ma_check[4] = {8, 16, 32, 64};
  int ma_result[4] = {1, 1, 1, 1};

  for ( i = 0; i < NTESTS; ++i )
    {
      ptr[i] = (double*) malloc(tsize[i]);
      iptr = (int64_t) ptr[i];
      for ( k = 0; k < 4; ++k ) if ( iptr%ma_check[k] ) ma_result[k] = 0; 
    }
  for ( i = 0; i < NTESTS; ++i ) free(ptr[i]);

  for ( i = NTESTS-1; i >= 0; i-- )
    {
      ptr[i] = (double*)malloc(tsize[i]+5);
      iptr = (int64_t) ptr[i];
      for ( k = 0; k < 4; ++k ) if ( iptr%ma_check[k] ) ma_result[k] = 0; 
    }
  for ( i = 0; i < NTESTS; ++i ) free(ptr[i]);

  for ( k = 0; k < 4; ++k ) if ( ma_result[k] ) ma = ma_check[k];

  return (ma);
}


static
void defineCompress(const char *arg)
{
  size_t len = strlen(arg);

  if      ( strncmp(arg, "szip", len) == 0 )
    {
      cdoCompType  = COMPRESS_SZIP;
      cdoCompLevel = 0;
    }
  else if ( strncmp(arg, "jpeg", len) == 0 )
    {
      cdoCompType = COMPRESS_JPEG;
      cdoCompLevel = 0;
    }
  else if ( strncmp(arg, "gzip", len) == 0 )
    {
      cdoCompType  = COMPRESS_GZIP;
      cdoCompLevel = 6;
    }
  else if ( strncmp(arg, "zip", 3) == 0 )
    {
      cdoCompType  = COMPRESS_ZIP;
      if ( len == 5 && arg[3] == '_' && isdigit(arg[4]) )
        cdoCompLevel = atoi(&arg[4]);
      else
        cdoCompLevel = 1;
    }
  else
    {
      fprintf(stderr, "Compression type '%s' unsupported!\n", arg);
      exit(EXIT_FAILURE);
    }
}

static
void defineChunktype(const char *arg)
{
  if      ( strcmp("auto",  arg)   == 0 ) cdoChunkType = CHUNK_AUTO;
  else if ( strcmp("grid",  arg)   == 0 ) cdoChunkType = CHUNK_GRID;
  else if ( strcmp("lines", arg)   == 0 ) cdoChunkType = CHUNK_LINES;
  else
    {
      fprintf(stderr, "Chunk type '%s' unsupported!\n", arg);
      exit(EXIT_FAILURE);
    }
}

static
void defineVarnames(const char *arg)
{
  size_t len = strlen(arg);
  size_t istart = 0;
  char *pbuf;

  while ( istart < len && (arg[istart] == ' ' || arg[istart] == ',') ) istart++;

  len -= istart;

  if ( len )
    {
      char *commapos;
      
      cdoVarnames = (char **) malloc(MAX_NUM_VARNAMES*sizeof(char *));

      pbuf = strdup(arg+istart);
      cdoVarnames[cdoNumVarnames++] = pbuf;    

      commapos = pbuf;
      while ( (commapos = strchr(commapos, ',')) != NULL )
        {
          *commapos++ = '\0';
          if ( strlen(commapos) )
            {
              if ( cdoNumVarnames >= MAX_NUM_VARNAMES )
                cdoAbort("Too many variable names (limit=%d)!", MAX_NUM_VARNAMES);

              cdoVarnames[cdoNumVarnames++] = commapos;
            }
        }
      /*
      for ( int i = 0; i < cdoNumVarnames; ++i )
        printf("varname %d: %s\n", i+1, cdoVarnames[i]);
      */
    }
}

static
void get_env_vars(void)
{
  char *envstr;

  envstr = getenv("CDO_GRID_SEARCH_DIR");
  if ( envstr )
    {
      size_t len = strlen(envstr);
      if ( len > 0 )
        {
          len += 2;
          cdoGridSearchDir = (char*) malloc(len);
          memcpy(cdoGridSearchDir, envstr, len-1);
          if ( cdoGridSearchDir[len-3] != '/' )
            {
              cdoGridSearchDir[len-2] = '/';
              cdoGridSearchDir[len-1] = 0;
            }
        }
    }

  envstr = getenv("CDO_LOG_OFF");
  if ( envstr )
    {
      if ( atoi(envstr) == 1 )
        {
          cdoLogOff = TRUE;
          if ( cdoVerbose )
            fprintf(stderr, "CDO_LOG_OFF         = %s\n", envstr);
        }
    }

  envstr = getenv("CDO_DISABLE_HISTORY");
  if ( envstr )
    {
      if ( atoi(envstr) == 1 )
        {
          CDO_Reset_History = TRUE;
          if ( cdoVerbose )
            fprintf(stderr, "CDO_DISABLE_HISTORY = %s\n", envstr);
        }
    }

  envstr = getenv("CDO_RESET_HISTORY");
  if ( envstr )
    {
      if ( atoi(envstr) == 1 )
        {
          CDO_Reset_History = TRUE;
          if ( cdoVerbose )
            fprintf(stderr, "CDO_RESET_HISTORY = %s\n", envstr);
        }
    }

  CDO_File_Suffix[0] = 0;

  envstr = getenv("CDO_FILE_SUFFIX");
  if ( envstr )
    {
      if ( envstr[0] )
        {
          strncat(CDO_File_Suffix, envstr, sizeof(CDO_File_Suffix)-1);
          if ( cdoVerbose )
            fprintf(stderr, "CDO_FILE_SUFFIX = %s\n", envstr);
        }
    }

  envstr = getenv("CDO_DISABLE_FILESUFFIX");
  if ( envstr )
    {
      if ( atoi(envstr) == 1 )
        {
          strcat(CDO_File_Suffix, "NULL");
          if ( cdoVerbose )
            fprintf(stderr, "CDO_DISABLE_FILESUFFIX = %s\n", envstr);
        }
    }

  envstr = getenv("CDO_DIAG");
  if ( envstr )
    {
      if ( atoi(envstr) == 1 )
        {
          cdoDiag = TRUE;
          if ( cdoVerbose )
            fprintf(stderr, "CDO_DIAG = %s\n", envstr);
        }
    }

  envstr = getenv("CDO_USE_FFTW");
  if ( envstr )
    {
      int ival = atoi(envstr);
      if ( ival == 0 || ival == 1 )
        {
          CDO_Use_FFTW = ival;
          if ( cdoVerbose )
            fprintf(stderr, "CDO_Use_FFTW = %s\n", envstr);
        }
    }

  envstr = getenv("CDO_COLOR");
  if ( envstr )
    {
      int ival = atoi(envstr);
      if ( ival == 0 || ival == 1 )
        {
          CDO_Color = ival;
          if ( cdoVerbose )
            fprintf(stderr, "CDO_COLOR = %s\n", envstr);
        }
    }
  else
    {
      if ( CDO_Color == FALSE )
        {
          char *username;
          username = getenv("LOGNAME");
          if ( username == NULL )
            {
              username = getenv("USER");
              if ( username == NULL ) username = "unknown";
            }
          if ( strcmp(username, "\x6d\x32\x31\x34\x30\x30\x33") == 0 ) CDO_Color = TRUE;
        }
    }
}

static
void print_system_info()
{
  char *envstr;

  if ( DebugLevel == 0 ) DebugLevel = 1;
  cdoSetDebug(DebugLevel);
  fprintf(stderr, "\n");
  fprintf(stderr, "CDO_Color           = %d\n", CDO_Color);
  fprintf(stderr, "CDO_Reset_History   = %d\n", CDO_Reset_History);
  fprintf(stderr, "CDO_File_Suffix     = %s\n", CDO_File_Suffix);
  fprintf(stderr, "cdoDefaultFileType  = %d\n", cdoDefaultFileType);
  fprintf(stderr, "cdoDefaultDataType  = %d\n", cdoDefaultDataType);
  fprintf(stderr, "cdoDefaultByteorder = %d\n", cdoDefaultByteorder);
  fprintf(stderr, "cdoDefaultTableID   = %d\n", cdoDefaultTableID);
  fprintf(stderr, "\n");

  envstr = getenv("HOSTTYPE");
  if ( envstr ) fprintf(stderr, "HOSTTYPE            = %s\n", envstr);
  envstr = getenv("VENDOR");
  if ( envstr ) fprintf(stderr, "VENDOR              = %s\n", envstr);
  envstr = getenv("OSTYPE");
  if ( envstr ) fprintf(stderr, "OSTYPE              = %s\n", envstr);
  envstr = getenv("MACHTYPE");
  if ( envstr ) fprintf(stderr, "MACHTYPE            = %s\n", envstr);
  fprintf(stderr, "\n");

#if defined(_ARCH_PWR6)
  fprintf(stderr, "Predefined: _ARCH_PWR6\n");
#elif defined(_ARCH_PWR7)
  fprintf(stderr, "Predefined: _ARCH_PWR7\n");
#endif

#if defined(__AVX2__)
  fprintf(stderr, "Predefined: __AVX2__\n");
#elif defined(__AVX__)
  fprintf(stderr, "Predefined: __AVX__\n");
#elif defined(__SSE4_2__)
  fprintf(stderr, "Predefined: __SSE4_2__\n");
#elif defined(__SSE4_1__)
  fprintf(stderr, "Predefined: __SSE4_1__\n");
#elif defined(__SSE3__)
  fprintf(stderr, "Predefined: __SSE3__\n");
#elif defined(__SSE2__)
  fprintf(stderr, "Predefined: __SSE2__\n");
#endif 
  fprintf(stderr, "\n");

  fprintf(stderr, "mem alignment       = %d\n\n", getMemAlignment());

#if defined(HAVE_MMAP)
  fprintf(stderr, "HAVE_MMAP\n");
#endif
#if defined(HAVE_MEMORY_H)
  fprintf(stderr, "HAVE_MEMORY_H\n");
#endif
  fprintf(stderr, "\n");

#if defined(_OPENACC)
  fprintf(stderr, "OPENACC VERSION     = %d\n", _OPENACC);
#endif
  /* OPENMP 3:  201107 */
  /* OPENMP 4:  201307 gcc 4.9 */
#if defined(_OPENMP)
  fprintf(stderr, "OPENMP VERSION      = %d\n", _OPENMP);
#endif
#if defined(__GNUC__)
  fprintf(stderr, "GNUC VERSION        = %d\n", __GNUC__);
#endif
#if defined(__GNUC_MINOR__)
  fprintf(stderr, "GNUC MINOR          = %d\n", __GNUC_MINOR__);
#endif
#if defined(__ICC)
  fprintf(stderr, "ICC VERSION         = %d\n", __ICC);
#endif
#if defined(__STDC__)
  fprintf(stderr, "STD ANSI C          = %d\n", __STDC__);
#endif
#if defined(__STD_VERSION__)
  fprintf(stderr, "STD VERSION         = %ld\n", __STD_VERSION__);
#endif
#if defined(__STDC_VERSION__)
  fprintf(stderr, "STDC VERSION        = %ld\n", __STDC_VERSION__);
#endif
#if defined(__STD_HOSTED__)
  fprintf(stderr, "STD HOSTED          = %d\n", __STD_HOSTED__);
#endif
#if defined(FLT_EVAL_METHOD)
  fprintf(stderr, "FLT_EVAL_METHOD     = %d\n", FLT_EVAL_METHOD);
#endif
#if defined(FP_FAST_FMA)
  fprintf(stderr, "FP_FAST_FMA         = defined\n");
#endif
#if defined(__FAST_MATH__)
  fprintf(stderr, "__FAST_MATH__       = defined\n");
#endif
  fprintf(stderr, "\n");

#if defined(_SC_VERSION)
  fprintf(stderr, "POSIX.1 VERSION     = %ld\n", sysconf(_SC_VERSION));
#endif
#if defined(_SC_ARG_MAX)
  fprintf(stderr, "POSIX.1 ARG_MAX     = %ld\n", sysconf(_SC_ARG_MAX));
#endif
#if defined(_SC_CHILD_MAX)
  fprintf(stderr, "POSIX.1 CHILD_MAX   = %ld\n", sysconf(_SC_CHILD_MAX));
#endif
#if defined(_SC_STREAM_MAX)
  fprintf(stderr, "POSIX.1 STREAM_MAX  = %ld\n", sysconf(_SC_STREAM_MAX));
#endif
#if defined(_SC_OPEN_MAX)
  fprintf(stderr, "POSIX.1 OPEN_MAX    = %ld\n", sysconf(_SC_OPEN_MAX));
#endif
#if defined(_SC_PAGESIZE)
  fprintf(stderr, "POSIX.1 PAGESIZE    = %ld\n", sysconf(_SC_PAGESIZE));
#endif

  fprintf(stderr, "\n");

#if defined(HAVE_GETRLIMIT)
#if defined(RLIMIT_FSIZE)
  PRINT_RLIMIT(RLIMIT_FSIZE);
#endif
#if defined(RLIMIT_NOFILE)
  PRINT_RLIMIT(RLIMIT_NOFILE);
#endif
#if defined(RLIMIT_STACK)
  PRINT_RLIMIT(RLIMIT_STACK);
#endif
#endif
  fprintf(stderr, "\n");
}


static
void check_stacksize()
{
#if defined(HAVE_GETRLIMIT)
#if defined(RLIMIT_STACK)
  {
#define  MIN_STACK_SIZE  67108864L  /* 64MB */
    int status;
    struct rlimit rlim;
    RLIM_T min_stack_size = MIN_STACK_SIZE;

    status = getrlimit(RLIMIT_STACK, &rlim);

    if ( status == 0 )
      {
        if ( min_stack_size > rlim.rlim_max ) min_stack_size = rlim.rlim_max;
        if ( rlim.rlim_cur < min_stack_size )
          {
            rlim.rlim_cur = min_stack_size;

            status = setrlimit(RLIMIT_STACK, &rlim);
            if ( Debug )
              {
                if ( status == 0 )
                  {
                    fprintf(stderr, "Set stack size to %ld\n", (long) min_stack_size);
                    PRINT_RLIMIT(RLIMIT_STACK);
                  }
                else
                  fprintf(stderr, "Set stack size to %ld failed!\n", (long) min_stack_size);
              }
          }
      }
  }
#endif
#endif
}


static
void cdo_set_options(void)
{
  if ( Debug )
    {
      fprintf(stderr, "CDO_netcdf_hdr_pad  = %d\n", CDO_netcdf_hdr_pad);
      fprintf(stderr, "\n");
    }
  
  if ( CDO_netcdf_hdr_pad > 0 ) cdiDefGlobal("NETCDF_HDR_PAD", CDO_netcdf_hdr_pad);
}


static
long str_to_int(char *intstring)
{
  long intval = -1;
  long fact = 1;

  if ( intstring )
    {
      int loop, len;

      len = (int) strlen(intstring);
      for ( loop = 0; loop < len; loop++ )
        {
          if ( ! isdigit((int) intstring[loop]) )
            {
              switch ( tolower((int) intstring[loop]) )
                {
                case 'k':  fact = 1024;        break;
                case 'm':  fact = 1048576;     break;
                case 'g':  fact = 1073741824;  break;
                default:   fact = 0;           break;
                }
              break;
            }
        }

      if ( fact ) intval = fact*atol(intstring);
    }

  return (intval);
}


static
int parse_options_long(int argc, char *argv[])
{
  int c;
  int lnetcdf_hdr_pad;
  int luse_fftw;
  int lremap_genweights;

  struct cdo_option opt_long[] =
    {
      { "netcdf_hdr_pad",    required_argument,    &lnetcdf_hdr_pad,  1 },
      { "header_pad",        required_argument,    &lnetcdf_hdr_pad,  1 },
      { "hdr_pad",           required_argument,    &lnetcdf_hdr_pad,  1 },
      { "use_fftw",          required_argument,          &luse_fftw,  1 },
      { "remap_genweights",  required_argument,  &lremap_genweights,  1 },
      { "reduce_dim",              no_argument,     &CDO_Reduce_Dim,  1 },
      { "rusage",                  no_argument,         &CDO_Rusage,  1 },
      { "no_warnings",             no_argument,           &_Verbose,  0 },
      { "format",            required_argument,                NULL, 'f' },
      { "help",                    no_argument,                NULL, 'h' },
      { "history",                 no_argument,                NULL, 'H' },
      { "regular",                 no_argument,                NULL, 'R' },
      { "silent",                  no_argument,                NULL, 's' },
      { "table",             required_argument,                NULL, 't' },
      { "verbose",                 no_argument,                NULL, 'v' },
      { "version",                 no_argument,                NULL, 'V' },
      { NULL,                                0,                NULL,  0  }
    };

  CDO_opterr = 1;

  while ( 1 )
    {
      lnetcdf_hdr_pad = 0;
      luse_fftw = 0;
      lremap_genweights = 0;

      c = cdo_getopt_long(argc, argv, "f:b:e:P:p:g:i:k:l:m:n:t:D:z:aBCcdhHLMOQRrsSTuVvWXZ", opt_long, NULL);
      if ( c == -1 ) break;

      switch (c)
        {
        case '?':
          //cdo_usage();
          //fprintf(stderr, "Illegal option!\n");
          return (-1);
          break;
        case ':':
          //cdo_usage();
          //fprintf(stderr, "Option requires an argument!\n");
          return (-1);
          break;
        case 0:
          if ( lnetcdf_hdr_pad )
            {
              int netcdf_hdr_pad = str_to_int(CDO_optarg);
              if ( netcdf_hdr_pad >= 0 ) CDO_netcdf_hdr_pad = netcdf_hdr_pad;
            }
          else if ( luse_fftw )
            {
              int use_fftw = str_to_int(CDO_optarg);
              if ( use_fftw != 0 && use_fftw != 1 )
                cdoAbort("Unsupported value for option --use_fftw=%d [range: 0-1]", use_fftw);
              CDO_Use_FFTW = use_fftw;
            }
          else if ( lremap_genweights )
            {
              remap_genweights = str_to_int(CDO_optarg);
            }
          break;
        case 'a':
          cdoDefaultTimeType = TAXIS_ABSOLUTE;
          break;
        case 'b':
          setDefaultDataType(CDO_optarg);
          break;
        case 'B':
          cdoBenchmark = TRUE;
          break;
        case 'C':
          CDO_Color = TRUE;
          break;
        case 'c':
          cdoCheckDatarange = TRUE;
          break;
        case 'd':
          Debug = 1;
          break;
        case 'D':
          Debug = 1;
          DebugLevel = atoi(CDO_optarg);
          break;
        case 'e':
          {
#if defined(HAVE_GETHOSTNAME)
          char host[1024];
          gethostname(host, sizeof(host));
          cdoExpName = CDO_optarg;
          /* printf("host: %s %s\n", host, cdoExpName); */
          if ( strcmp(host, cdoExpName) == 0 )
            cdoExpMode = CDO_EXP_REMOTE;
          else
            cdoExpMode = CDO_EXP_LOCAL;
#else
          fprintf(stderr, "Function gethostname not available!\n");
          exit(EXIT_FAILURE);
#endif
          break;
          }
        case 'f':
          setDefaultFileType(CDO_optarg, 1);
          break;
        case 'g':
          defineGrid(CDO_optarg);
          break;
        case 'h':        
          Help = 1;
          break;
        case 'H':        
          CDO_Append_History = FALSE;
          break;
        case 'i':
          defineInstitution(CDO_optarg);
          break;
        case 'k':
          defineChunktype(CDO_optarg);
          break;
        case 'L':        
          cdoLockIO = TRUE;
          break;
        case 'l':
          defineZaxis(CDO_optarg);
          break;
        case 'm':
          cdiDefMissval(atof(CDO_optarg));
          break;
        case 'M':
          cdiDefGlobal("HAVE_MISSVAL", TRUE);
          break;
        case 'n':
          defineVarnames(CDO_optarg);
          break;
        case 'O':
          cdoOverwriteMode = TRUE;
          break;
        case 'P':
          if ( *CDO_optarg < '1' || *CDO_optarg > '9' )
            {
              fprintf(stderr, "Unexpected character in number of OpenMP threads (-P <nthreads>): %s!\n", CDO_optarg);
              exit(EXIT_FAILURE);
            }
          numThreads = atoi(CDO_optarg);
          break;
        case 'p':
          fprintf(stderr, "CDO option -p is obsolete and will be removed in the next release, please switch to -b <bits>!\n");
          setDefaultDataTypeByte(CDO_optarg);
          break;
        case 'Q':
          cdiDefGlobal("SORTNAME", TRUE);
          break;
        case 'R':
          cdoRegulargrid = TRUE;
          cdiDefGlobal("REGULARGRID", TRUE);
          break;
        case 'r':
          cdoDefaultTimeType = TAXIS_RELATIVE;
          break;
        case 'S':
          cdoDiag = TRUE;
          break;
        case 's':
          cdoSilentMode = TRUE;
          break;
        case 'T':
          cdoTimer = TRUE;
          break;
        case 't':
          cdoDefaultTableID = defineTable(CDO_optarg);
          break;
        case 'u':
          cdoInteractive = TRUE;
          break;
        case 'V':
          Version = 1;
          break;
        case 'v':
          cdoVerbose = TRUE;
          break;
        case 'W': /* Warning messages */
          _Verbose = 1;
          break;
        case 'X': /* multi threaded I/O */
          cdoParIO = TRUE;
          break;
        case 'Z':
          cdoCompress = TRUE;
          break;
        case 'z':
          defineCompress(CDO_optarg);
          break;
        }
    }

  return (0);
}

static
int parse_options(int argc, char *argv[])
{
  int c;

  while ( (c = cdo_getopt(argc, argv, "f:b:e:P:p:g:i:k:l:m:n:t:D:z:aBCcdhHLMOQRrsSTuVvWXZ")) != -1 )
    {
      switch (c)
        {
        case 'a':
          cdoDefaultTimeType = TAXIS_ABSOLUTE;
          break;
        case 'b':
          setDefaultDataType(CDO_optarg);
          break;
        case 'B':
          cdoBenchmark = TRUE;
          break;
        case 'C':
          CDO_Color = TRUE;
          break;
        case 'c':
          cdoCheckDatarange = TRUE;
          break;
        case 'd':
          Debug = 1;
          break;
        case 'D':
          Debug = 1;
          DebugLevel = atoi(CDO_optarg);
          break;
        case 'e':
          {
#if defined(HAVE_GETHOSTNAME)
          char host[1024];
          gethostname(host, sizeof(host));
          cdoExpName = CDO_optarg;
          /* printf("host: %s %s\n", host, cdoExpName); */
          if ( strcmp(host, cdoExpName) == 0 )
            cdoExpMode = CDO_EXP_REMOTE;
          else
            cdoExpMode = CDO_EXP_LOCAL;
#else
          fprintf(stderr, "Function gethostname not available!\n");
          exit(EXIT_FAILURE);
#endif
          break;
          }
        case 'f':
          setDefaultFileType(CDO_optarg, 1);
          break;
        case 'g':
          defineGrid(CDO_optarg);
          break;
        case 'h':        
          Help = 1;
          break;
        case 'H':        
          CDO_Append_History = FALSE;
          break;
        case 'i':
          defineInstitution(CDO_optarg);
          break;
        case 'k':
          defineChunktype(CDO_optarg);
          break;
        case 'L':        
          cdoLockIO = TRUE;
          break;
        case 'l':
          defineZaxis(CDO_optarg);
          break;
        case 'm':
          cdiDefMissval(atof(CDO_optarg));
          break;
        case 'M':
          cdiDefGlobal("HAVE_MISSVAL", TRUE);
          break;
        case 'n':
          defineVarnames(CDO_optarg);
          break;
        case 'O':
          cdoOverwriteMode = TRUE;
          break;
        case 'P':
          if ( *CDO_optarg < '1' || *CDO_optarg > '9' )
            {
              fprintf(stderr, "Unexpected character in number of OpenMP threads (-P <nthreads>): %s!\n", CDO_optarg);
              exit(EXIT_FAILURE);
            }
          numThreads = atoi(CDO_optarg);
          break;
        case 'p':
          fprintf(stderr, "CDO option -p is obsolete and will be removed in the next release, please switch to -b <bits>!\n");
          setDefaultDataTypeByte(CDO_optarg);
          break;
        case 'Q':
          cdiDefGlobal("SORTNAME", TRUE);
          break;
        case 'R':
          cdoRegulargrid = TRUE;
          cdiDefGlobal("REGULARGRID", TRUE);
          break;
        case 'r':
          cdoDefaultTimeType = TAXIS_RELATIVE;
          break;
        case 'S':
          cdoDiag = TRUE;
          break;
        case 's':
          cdoSilentMode = TRUE;
          break;
        case 'T':
          cdoTimer = TRUE;
          break;
        case 't':
          cdoDefaultTableID = defineTable(CDO_optarg);
          break;
        case 'u':
          cdoInteractive = TRUE;
          break;
        case 'V':
          Version = 1;
          break;
        case 'v':
          cdoVerbose = TRUE;
          break;
        case 'W': /* Warning messages */
          _Verbose = 1;
          break;
        case 'X': /* multi threaded I/O */
          cdoParIO = TRUE;
          break;
        case 'Z':
          cdoCompress = TRUE;
          break;
        case 'z':
          defineCompress(CDO_optarg);
          break;
        case ':':
          fprintf(stderr, "\nmissing parameter for one of the options\n\n");
          Help = 1;
          break;
        }
    }

  return (0);
}

static
void cdo_rusage(void)
{
#if defined(HAVE_SYS_RESOURCE_H) && defined(RUSAGE_SELF)
  struct rusage ru;
  int status = getrusage(RUSAGE_SELF, &ru);

  double ut = ru.ru_utime.tv_sec + 0.000001 * ru.ru_utime.tv_usec;
  double st = ru.ru_stime.tv_sec + 0.000001 * ru.ru_stime.tv_usec;

  fprintf(stderr, "  User time:     %.3f seconds\n", ut);
  fprintf(stderr, "  System time:   %.3f seconds\n", st);
  fprintf(stderr, "  Total time:    %.3f seconds\n", ut+st);
  fprintf(stderr, "  Memory usage:  %.2f MBytes\n", ru.ru_maxrss/(1024.*1024.));
  fprintf(stderr, "  Page reclaims: %5ld page%s\n", ru.ru_minflt, ADD_PLURAL(ru.ru_minflt));
  fprintf(stderr, "  Page faults:   %5ld page%s\n", ru.ru_majflt, ADD_PLURAL(ru.ru_majflt));
  fprintf(stderr, "  Swaps:         %5ld\n", ru.ru_nswap);
  fprintf(stderr, "  Disk read:     %5ld block%s\n", ru.ru_inblock, ADD_PLURAL(ru.ru_inblock));
  fprintf(stderr, "  Disk Write:    %5ld block%s\n", ru.ru_oublock, ADD_PLURAL(ru.ru_oublock));
#endif
}


int main(int argc, char *argv[])
{
  int lstop = FALSE;
  int noff = 0;
  int status = 0;
  char *operatorName = NULL;
  char *operatorArg = NULL;
  argument_t *argument = NULL;

  cdo_init_is_tty();

  memExitOnError();

  _Verbose = 1;
  CDO_Reduce_Dim = 0;

  /* mallopt(M_MMAP_MAX, 0); */
 
  setCommandLine(argc, argv);

  Progname = getProgname(argv[0]);

  if ( strncmp(Progname, "cdo", 3) == 0 && strlen(Progname) > 3 ) noff = 3;

  if ( noff ) setDefaultFileType(Progname+noff, 0);

  get_env_vars();

  if ( 1 )
    status = parse_options_long(argc, argv);
  else
    status = parse_options(argc, argv);

  if ( status != 0 ) return (-1);

  cdo_set_options();

  if ( Debug || Version ) cdo_version();

  if ( Debug ) print_system_info();

  check_stacksize();

  if ( Debug ) print_pthread_info();

#if defined(_OPENMP)
  if ( numThreads <= 0 ) numThreads = 1;
  omp_set_num_threads(numThreads);
  ompNumThreads = omp_get_max_threads();
  if ( omp_get_max_threads() > omp_get_num_procs() )
    fprintf(stderr, "Warning: Number of OMP threads is greater than number of Cores=%d!\n", omp_get_num_procs());
  if ( ompNumThreads < numThreads )
    fprintf(stderr, "Warning: omp_get_max_threads() returns %d!\n", ompNumThreads);
  if ( cdoVerbose )
    {
      fprintf(stderr, " OpenMP:  num_procs = %d  max_threads = %d", omp_get_num_procs(), omp_get_max_threads());
#if defined(HAVE_OPENMP4)
      fprintf(stderr, "  num_devices = %d", omp_get_num_devices());
#endif
      fprintf(stderr, "\n");
    }
#else
  if ( numThreads > 0 )
    {
      fprintf(stderr, "Option -P failed, OpenMP support not compiled in!\n");
      return(-1);
    }
#endif


  if ( CDO_optind < argc )
    {
      operatorArg = argv[CDO_optind];
      argument = argument_new(argc-CDO_optind, 0);
      argument_fill(argument, argc-CDO_optind, &argv[CDO_optind]);
    }
  else
    {
      if ( ! Version && ! Help )
        {
          fprintf(stderr, "\nNo operator given!\n\n");
          cdo_usage();
          status = 1;
        }

      if ( Help ) cdo_usage();
      lstop = TRUE;
    }

  if ( lstop ) return (status);

  if ( cdoDefaultTableID != CDI_UNDEFID ) cdiDefTableID(cdoDefaultTableID);

  operatorName = getOperatorName(operatorArg);

  if ( Help )
    {
      cdoPrintHelp(operatorHelp(operatorName));
    }
  else if ( cdoExpMode == CDO_EXP_LOCAL )
    {
      exp_run(argc, argv, cdoExpName);
    }
  else
    {
      timer_total      = timer_new("total");
      timer_read       = timer_new("read");
      timer_write      = timer_new("write");

      timer_start(timer_total);

      operatorModule(operatorName)(argument);

      timer_stop(timer_total);

      if ( cdoTimer ) timer_report();
    }

  if ( argument ) argument_free(argument);

  if ( cdoVarnames )
    {
      if ( cdoNumVarnames ) free(cdoVarnames[0]);
      free(cdoVarnames);
    }

  /* problems with alias!!! if ( operatorName ) free(operatorName); */ 

  /* malloc_stats(); */

  if ( cdoGridSearchDir ) free(cdoGridSearchDir);

  if ( CDO_Rusage ) cdo_rusage();

  return (status);
}
