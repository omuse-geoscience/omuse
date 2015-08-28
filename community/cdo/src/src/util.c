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

#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 600 /* ftello */
#endif

#if defined(HAVE_CONFIG_H)
#  include "config.h"
#endif

#if defined(_OPENMP)
#  include <omp.h>
#endif

#if defined(HAVE_FNMATCH_H)
#include <fnmatch.h>
#endif


#include <stdio.h>
#include <string.h>
#include <ctype.h>   /* tolower */

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "modules.h"
#include "util.h"


#if ! defined(VERSION)
#  define  VERSION  "0.0.1"
#endif
 

/* refactor: moved here from *.c */

int CDO_opterr = 0;      // refactor: moved here from cdo_getopt.c
char *CDO_optarg = NULL; // refactor: moved here from cdo_getopt.c
int CDO_optind = 1;      // refactor: moved here from cdo_getopt.c


/* refactor: moved here from cdo.c */

char *Progname;
char CDO_Version[] = "Climate Data Operators version "VERSION" (http://mpimet.mpg.de/cdo)"; // refactor: moved here from cdo.c

int ompNumThreads = 1;

int stdin_is_tty  = 0;
int stdout_is_tty = 0;
int stderr_is_tty = 0;

char* cdoGridSearchDir   = NULL;

int cdoDefaultFileType   = CDI_UNDEFID;
int cdoDefaultDataType   = CDI_UNDEFID;
int cdoDefaultByteorder  = CDI_UNDEFID;
int cdoDefaultTableID    = CDI_UNDEFID;
int cdoDefaultInstID     = CDI_UNDEFID;     // moved here from institution.c, was UNDEFID
int cdoDefaultTimeType   = CDI_UNDEFID;

int cdoLockIO            = FALSE;
int cdoCheckDatarange    = FALSE;

int CDO_Color            = FALSE;
int CDO_Use_FFTW         = TRUE;
int cdoDiag              = FALSE;

int CDO_Reduce_Dim       = FALSE;
int CDO_Append_History   = TRUE;
int CDO_Reset_History    = FALSE;

int cdoCompType          = COMPRESS_NONE;  // compression type
int cdoCompLevel         = 0;              // compression level
int cdoDebug             = 0;
int cdoChunkType         = CDI_UNDEFID;
int cdoLogOff            = FALSE;
int cdoSilentMode        = FALSE;
int cdoOverwriteMode     = FALSE;
int cdoBenchmark         = FALSE;
int cdoTimer             = FALSE;
int cdoVerbose           = FALSE;
int cdoCompress          = FALSE;
int cdoInteractive       = FALSE;
int cdoParIO             = FALSE;
int cdoRegulargrid       = FALSE;

int cdoNumVarnames       = 0;
char **cdoVarnames       = NULL;

char CDO_File_Suffix[32];

int cdoExpMode           = -1;
char *cdoExpName         = NULL;

int timer_read, timer_write;


#if defined(HAVE_FNMATCH_H)
int wildcardmatch(const char *pattern, const char *string)
{
  return fnmatch(pattern, string, 0);
}
#else
// The wildcardmatch function checks if two given strings match. 
// The first string may contain wildcard characters
// * --> Matches with 0 or more instances of any character or set of characters.
// ? --> Matches with any one character.
// source code from http://www.geeksforgeeks.org/wildcard-character-matching/
int wildcardmatch(const char *w, const char *s)
{
    // If we reach at the end of both strings, we are done
    if ( *w == '\0' && *s == '\0' ) return 0;
 
    // Make sure that the characters after '*' are present in second string.
    // This function assumes that the first string will not contain two consecutive '*'
    if ( *w == '*' && *(w+1) != '\0' && *s == '\0' ) return 1;
 
    // If the first string contains '?', or current characters of both strings match
    if ( (*w == '?' && *s != '\0') || *w == *s ) return wildcardmatch(w+1, s+1);
 
    // If there is *, then there are two possibilities
    // a) We consider current character of second string
    // b) We ignore current character of second string.
    if ( *w == '*' ) return wildcardmatch(w+1, s) || wildcardmatch(w, s+1);

    return 1;
}
#endif

int cdo_omp_get_thread_num(void)
{
  int threadnum = 0;

#if defined(_OPENMP)
  threadnum = omp_get_thread_num();
#endif

  return (threadnum);
}


char *getProgname(char *string)
{
  char *progname;

#if defined(_WIN32)
  /*  progname = strrchr(string, '\\'); */
  progname = " cdo";
#else
  progname = strrchr(string, '/');
#endif

  if ( progname == NULL ) progname = string;
  else                    progname++;

  return (progname);
}

char *getOperator(const char *argument)
{
  char *operatorArg = NULL;
  size_t len;

  if ( argument )
    {
      len = 1 + strlen(argument);

      operatorArg = (char*) malloc(len);

      memcpy(operatorArg, argument, len);
    }

  return (operatorArg);
}

char *operatorAlias(char *operatorName);

char *getOperatorName(const char *operatorArg)
{
  char *commapos;
  char *operatorName = NULL;
  size_t len;

  if ( operatorArg )
    {
      if ( operatorArg[0] == '-' ) operatorArg++;

      commapos = strchr(operatorArg, ',');

      if ( commapos )
        len = commapos - operatorArg;
      else
        len = strlen(operatorArg);

      operatorName = (char*) malloc(len+1);

      memcpy(operatorName, operatorArg, len);
      operatorName[len] = '\0';
    }

  /*  return (operatorName); */
  return (operatorAlias(operatorName));
}


argument_t *file_argument_new(const char *filename)
{
  argument_t *argument;

  argument = (argument_t*) calloc(1, sizeof(argument_t));

  argument->argc = 1;
  argument->argv = (char **) calloc(1, sizeof(char *));
  argument->argv[0] = (char *) filename;
  argument->args = (char *) filename;

  return (argument);
}


void file_argument_free(argument_t *argument)
{
  if ( argument )
    {
      if ( argument->argc )
        {
          assert(argument->argc == 1);
          free(argument->argv);
        }
      free(argument);
    }
}


argument_t *argument_new(size_t argc, size_t len)
{
  argument_t *argument;

  argument = (argument_t*) calloc(1, sizeof(argument_t));

  if ( argc > 0 )
    {
      argument->argc = argc;
      argument->argv = (char **) calloc(argc, sizeof(char *));
    }

  if ( len > 0 )
    argument->args = (char*) calloc(len, sizeof(char));

  return (argument);
}


void argument_free(argument_t *argument)
{
  if ( argument )
    {
      if ( argument->argc )
        {
          int argc =  argument->argc;
          for ( int i = 0; i < argc; ++i )
            {
              if ( argument->argv[i] )
                {
                  free(argument->argv[i]);
                  argument->argv[i] = NULL;
                }
            }

          free(argument->argv);
          argument->argv = NULL;
          argument->argc = 0;
        }

      if ( argument->args )
        {
          free(argument->args);
          argument->args = NULL;
        }

      free(argument);
    }
}


void argument_fill(argument_t *argument, int argc, char *argv[])
{
  int iarg;

  assert(argument->argc == argc);

  for ( iarg = 0; iarg < argc; ++iarg )
    argument->argv[iarg] = strdup(argv[iarg]);
}


char *getFileArg(char *argument)
{
  char *fileArg = NULL;
  char *parg;
  char *blankpos;
  size_t len;

  if ( argument )
    {
      blankpos = strchr(argument, ' ');

      if ( blankpos )
        {
          parg = blankpos + 1;
          len = strlen(parg);
          fileArg = (char*) malloc(len+1);
          strcpy(fileArg, parg);
        }
    }

  return (fileArg);
}


void input_int(char *arg, int intarr[], int maxint, int *nintfound)
{
  int nint = 0;

  intarr[nint++] = atoi(arg);

  while ( (arg = strchr(arg, ',')) && (nint < maxint) )
    intarr[nint++] = atoi(++arg);
    
  *nintfound = nint;
}


void strtolower(char *str)
{
  int i, len;

  if ( str )
    {
      len = (int) strlen(str);
      for ( i = 0; i < len; i++ )
        str[i] = tolower((int) str[i]);
    }
}


double parameter2double(const char *string)
{
  char *endptr = NULL;

  double fval = strtod(string, &endptr);

  if ( *endptr != 0 )
    cdoAbort("Float parameter >%s< contains invalid character at position %d!",
	     string, (int)(endptr-string+1));

  return (fval);
}


int parameter2int(const char *string)
{
  char *endptr = NULL;

  int ival = (int) strtol(string, &endptr, 10);

  if ( *endptr != 0 )
    cdoAbort("Integer parameter >%s< contains invalid character at position %d!",
	     string, (int)(endptr-string+1));

  return (ival);
}


int parameter2intlist(const char *string)
{
  char *endptr = NULL;

  int ival = (int) strtol(string, &endptr, 10);

  if ( *endptr != 0 && *endptr != '/' && (endptr - string) == 0 )
    cdoAbort("Integer parameter >%s< contains invalid character at position %d!",
	     string, (int)(endptr-string+1));

  return (ival);
}


const char *seas_name_dec[4] = {"DJF", "MAM", "JJA", "SON"};
const char *seas_name_jan[4] = {"JFM", "AMJ", "JAS", "OND"};

static int season_start = START_DEC;

int get_season_start(void)
{
  static int lgetenv = TRUE;

  if ( lgetenv )
    {
      lgetenv = FALSE;
  
      char *envstr = getenv("CDO_SEASON_START");
      if ( envstr )
        {
          if      ( strcmp(envstr, "DEC") == 0 ) season_start = START_DEC;
          else if ( strcmp(envstr, "JAN") == 0 ) season_start = START_JAN;
      
          if ( cdoVerbose )
            {
              if      ( season_start == START_DEC )
                cdoPrint("Set SEASON_START to December");
              else if ( season_start == START_JAN )
                cdoPrint("Set SEASON_START to January");
            }
        }
    }

  return season_start;
}


void get_season_name(const char *seas_name[])
{
  long i;

  if ( get_season_start() == START_DEC )
    for ( i = 0; i < 4; ++i ) seas_name[i] = seas_name_dec[i];
  else
    for ( i = 0; i < 4; ++i ) seas_name[i] = seas_name_jan[i];
}


int month_to_season(int month)
{
  int season_start = get_season_start();
  int seas = -1;

  if ( month < 0 || month > 16 ) cdoAbort("Month %d out of range!", month);

  if ( season_start == START_DEC )
    {
      if ( month <= 12 )
        seas = (month % 12) / 3;
      else
        seas = month - 13;
    }
  else
    {
      if ( month <= 12 )
        seas = (month - 1) / 3;
      else
        seas = month - 13;
    }

  if ( seas < 0 || seas > 3 ) cdoAbort("Season %d out of range!", seas+1);

  return seas;
}

//#include <sys/types.h>
#include <sys/stat.h>
//#include <unistd.h>

int fileExists(const char *restrict filename)
{
  int status = 0;
  struct stat buf;

  if ( stat(filename, &buf) == 0 )
    {
      if ( S_ISREG(buf.st_mode) && buf.st_size > 0 ) status = 1;
    }

  return status;
}


int userFileOverwrite(const char *restrict filename)
{
  int status = 0, len;
  char line[1024], *pline;

  fprintf(stderr, "File %s already exists, overwrite? (yes/no): ", filename);
  readline(stdin, line, 1024);
  pline = line;
  while ( isspace((int) *pline) ) pline++;
  len = strlen(pline);
  if ( len == 3 )
    {
      if ( pline[0] == 'y' && pline[1] == 'e' && pline[2] == 's' )
        status = 1;
      else if ( pline[0] == 'Y' && pline[1] == 'E' && pline[2] == 'S' )
        status = 1;
    }
  else if ( len == 1 )
    {
      if ( pline[0] == 'y' ) status = 1;
    }

  return (status);
}


int ps_lhead = FALSE;
int ps_nch   = 0;
int ps_cval  = -1;

void progressInit(void)
{
  ps_lhead = FALSE;
  ps_nch   = 0;
  ps_cval  = -1;
}


void progressStatus(double offset, double refval, double curval)
{
  int ival;

  if ( cdoSilentMode ) return;
  if ( !stdout_is_tty ) return;

  offset = offset < 0 ? 0: offset;
  offset = offset > 1 ? 1: offset;
  refval = refval < 0 ? 0: refval;
  refval = refval > 1 ? 1: refval;
  curval = curval < 0 ? 0: curval;
  curval = curval > 1 ? 1: curval;

  ival = (offset + refval*curval)*100;

  if ( ps_cval == -1 )
    {
      ps_nch = fprintf(stdout, "%s: %3d%%", processInqPrompt(), 0);
      fflush(stdout);
      ps_lhead = TRUE;
    }

  if ( ival != ps_cval )
    {
      ps_cval = ival;
      fprintf(stdout, "\b\b\b\b%3d%%", ps_cval);
      fflush(stdout);
    }

  if ( ps_cval == 100 && ps_lhead )
    {
      ps_lhead = FALSE;
      while ( ps_nch-- ) fprintf(stdout, "\b \b");
      fflush(stdout);
    }
}


int datatype2str(int datatype, char *datatypestr)
{
  int status = 0;

  if      ( datatype == DATATYPE_PACK   ) strcpy(datatypestr, "P0");
  else if ( datatype > 0 && datatype <= 32  ) sprintf(datatypestr, "P%d", datatype);
  else if ( datatype == DATATYPE_CPX32  ) strcpy(datatypestr, "C32");
  else if ( datatype == DATATYPE_CPX64  ) strcpy(datatypestr, "C64");
  else if ( datatype == DATATYPE_FLT32  ) strcpy(datatypestr, "F32");
  else if ( datatype == DATATYPE_FLT64  ) strcpy(datatypestr, "F64");
  else if ( datatype == DATATYPE_INT8   ) strcpy(datatypestr, "I8");
  else if ( datatype == DATATYPE_INT16  ) strcpy(datatypestr, "I16");
  else if ( datatype == DATATYPE_INT32  ) strcpy(datatypestr, "I32");
  else if ( datatype == DATATYPE_UINT8  ) strcpy(datatypestr, "U8");
  else if ( datatype == DATATYPE_UINT16 ) strcpy(datatypestr, "U16");
  else if ( datatype == DATATYPE_UINT32 ) strcpy(datatypestr, "U32");
  else                                  { strcpy(datatypestr, "-1"); status = -1;}

  return (status);
}


int str2datatype(const char *datatypestr)
{
  int datatype = -1;
  size_t len;

  len = strlen(datatypestr);

  if ( len > 1 )
    {
      int ilen = atoi(datatypestr+1);
      if      ( strncmp(datatypestr, "P0",  len) == 0 ) datatype = DATATYPE_PACK;
      else if ( strncmp(datatypestr, "P",     1) == 0 &&
                ilen > 0 && ilen <= 32 )               datatype = atoi(datatypestr+1);
      else if ( strncmp(datatypestr, "C32", len) == 0 ) datatype = DATATYPE_CPX32;
      else if ( strncmp(datatypestr, "C64", len) == 0 ) datatype = DATATYPE_CPX64;
      else if ( strncmp(datatypestr, "F32", len) == 0 ) datatype = DATATYPE_FLT32;
      else if ( strncmp(datatypestr, "F64", len) == 0 ) datatype = DATATYPE_FLT64;
      else if ( strncmp(datatypestr, "I8",  len) == 0 ) datatype = DATATYPE_INT8;
      else if ( strncmp(datatypestr, "I16", len) == 0 ) datatype = DATATYPE_INT16;
      else if ( strncmp(datatypestr, "I32", len) == 0 ) datatype = DATATYPE_INT32;
      else if ( strncmp(datatypestr, "U8",  len) == 0 ) datatype = DATATYPE_UINT8;
      else if ( strncmp(datatypestr, "U16", len) == 0 ) datatype = DATATYPE_UINT16;
      else if ( strncmp(datatypestr, "U32", len) == 0 ) datatype = DATATYPE_UINT32;
      else if ( strncmp(datatypestr, "real",   len) == 0 ) datatype = DATATYPE_FLT32;
      else if ( strncmp(datatypestr, "double", len) == 0 ) datatype = DATATYPE_FLT64;
    }

  return (datatype);
}


off_t fileSize(const char *restrict filename)
{
  off_t filesize = 0;

  if ( filename[0] == '(' && filename[1] == 'p' )
    {
    }
  else
    {
      struct stat buf;
      if ( stat(filename, &buf) == 0 ) filesize = buf.st_size;
    }
  
  return filesize;
}


/* 
 * Return the filetype extension (const char)
 * for a given filetype (int)
 * TODO: handle lists of extensions i.e. grb and grb2 for GRIB2-format
 */
const char *filetypeext(int filetype)
{
  switch ( filetype )
    {
    case FILETYPE_GRB:
    case FILETYPE_GRB2: return (".grb");   break;
    case FILETYPE_NC:
    case FILETYPE_NC2:
    case FILETYPE_NC4:
    case FILETYPE_NC4C: return (".nc");    break;
    case FILETYPE_SRV:  return (".srv");   break;
    case FILETYPE_EXT:  return (".ext");   break;
    case FILETYPE_IEG:  return (".ieg");   break;
    default:            return ("");
    }
}


/*
 * Remove file extension:
 * -------------------------------------------------
 * Remove file extension if it is the expected one
 * Do nothing otherwise
 */
void rm_filetypeext(char *file, const char *ext)
{
  // length of filename
  int namelen = (int) strlen(file);
  // length of the original file extension
  int extlen =  (int) strlen(ext);

  // delete original extension if it is the expected one
  if ( strcmp(&file[namelen-extlen], ext) == 0 )
      file[namelen-extlen] = 0;
}


/*
 * Replace or just add file extension:
 * -------------------------------------------------
 * Replace file extension with new one
 * or just add the new file extension 
 * if the original extension is not the expected one
 */
void repl_filetypeext(char file[], const char *oldext, const char *newext)
{
  // delete original extension if it is the expected one
  rm_filetypeext(file, oldext);

  // add new file extension
  strcat(file, newext);
}


void cdoGenFileSuffix(char *filesuffix, size_t maxlen, int filetype, int vlistID, const char *refname)
{
  if ( strncmp(CDO_File_Suffix, "NULL", 4) != 0 )
    {
      if ( CDO_File_Suffix[0] != 0 )
        {
          strncat(filesuffix, CDO_File_Suffix, maxlen-1);
        }
      else
        {
          int lready = FALSE;
          int lcompsz = FALSE;
          
          if ( filetype == cdoDefaultFileType && cdoDefaultDataType == -1 && cdoDefaultByteorder == -1 )
            {
              size_t len = 0;
              if ( refname != NULL && *refname != 0 && *refname != '-' && *refname != '.' ) len = strlen(refname);

              if ( len > 2 )
                {
                  char *result = strrchr(refname, '.');
                  if ( result != NULL && result[1] != 0 )
                    {
                      int firstchar = tolower(result[1]);
                      switch (firstchar)
                        {
                        case 'g':
                          if ( cdoDefaultFileType == FILETYPE_GRB || cdoDefaultFileType == FILETYPE_GRB2 ) lready = TRUE;
                          break;
                        case 'n':
                          if ( cdoDefaultFileType == FILETYPE_NC || cdoDefaultFileType == FILETYPE_NC2 ||
                               cdoDefaultFileType == FILETYPE_NC4 || cdoDefaultFileType == FILETYPE_NC4C ) lready = TRUE;
                          break;
                        case 's':
                          if ( cdoDefaultFileType == FILETYPE_SRV ) lready = TRUE;
                          break;
                        case 'e':
                          if ( cdoDefaultFileType == FILETYPE_EXT ) lready = TRUE;
                          break;
                        case 'i':
                          if ( cdoDefaultFileType == FILETYPE_IEG ) lready = TRUE;
                          break;
                        }
                    }

                  //if ( lready )  strncat(filesuffix, result, maxlen-1);
		  if ( lready && ((len=strlen(result)) < (maxlen-1)) )
		    {
		      while ( len-- )
			{
			  if ( *result == '.' || isalnum(*result) ) 
			    strncat(filesuffix, result, 1);
			  result++;
			}
		    }
                }
            }

          if ( !lready )
            {
              strncat(filesuffix, streamFilesuffix(cdoDefaultFileType), maxlen-1);
              if ( cdoDefaultFileType == FILETYPE_GRB && vlistIsSzipped(vlistID) ) lcompsz = TRUE;
            }

          if ( cdoDefaultFileType == FILETYPE_GRB && cdoCompType == COMPRESS_SZIP ) lcompsz = TRUE;
          if ( lcompsz ) strncat(filesuffix, ".sz", maxlen-1);
        }
    }
}


int cdoFiletype(void)
{
  if ( cdoDefaultFileType == CDI_UNDEFID )
    {
      cdoDefaultFileType = FILETYPE_GRB;
      if ( ! cdoSilentMode )
        cdoPrint("Set default filetype to GRIB");
    }

  return (cdoDefaultFileType);
}

