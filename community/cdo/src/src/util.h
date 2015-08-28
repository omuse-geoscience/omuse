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

#ifndef _UTIL_H
#define _UTIL_H


/* dummy use of unused parameters to silence compiler warnings */
#define  UNUSED(x) (void)x

#undef   TRUE
#define  TRUE   1
#undef   FALSE
#define  FALSE  0

#undef   MIN
#define  MIN(a,b)  ((a) < (b) ? (a) : (b))
#undef   MAX
#define  MAX(a,b)  ((a) > (b) ? (a) : (b))

#define  ADD_PLURAL(n)  ((n)!=1 ? "s" : "")

#define  UNCHANGED_RECORD  (processSelf() == 0 && cdoStreamName(0)->argv[0][0] != '-' && cdoRegulargrid == FALSE && cdoDefaultFileType == -1 && cdoDefaultDataType == -1 && cdoDefaultByteorder == -1 )

extern char *Progname;
extern char *cdoGridSearchDir;
extern int CDO_Reduce_Dim;
extern int CDO_Append_History;
extern int CDO_Reset_History;
extern int timer_read, timer_write; // refactor: both pstream.c and CDIread.c CDIwrite.c defined in cdo.c

extern int   CDO_optind;
extern char *CDO_optarg;
extern int CDO_opterr;
extern int remap_genweights;

extern char *cdoExpName;
extern int ompNumThreads;

extern int stdin_is_tty;
extern int stdout_is_tty;
extern int stderr_is_tty;

extern int cdoDefaultFileType;
extern int cdoDefaultDataType;
extern int cdoDefaultByteorder;
extern int cdoDefaultTableID;
extern int cdoDefaultInstID;
extern int cdoDefaultTimeType;
extern int cdoLogOff;

extern int cdoLockIO;
extern int cdoCheckDatarange;

extern int cdoSilentMode;
extern int cdoOverwriteMode;
extern int cdoRegulargrid;
extern int cdoBenchmark;
extern int cdoTimer;
extern int cdoVerbose;
extern int cdoDebug;
extern int cdoCompress;
extern int cdoInteractive;
extern int cdoParIO;

extern int cdoCompType;
extern int cdoCompLevel;

extern int cdoChunkType;

extern int cdoExpMode;

extern int CDO_Color;
extern int CDO_Use_FFTW;
extern int cdoDiag;

extern int cdoNumVarnames;
extern char **cdoVarnames;
extern char CDO_File_Suffix[32]; // refactor: added keyword extern


/* moved here from *.c */
extern char CDO_Version[]; // refactor: moved here from pstream.c


typedef struct {
  int    argc;
  int    argl;
  char **argv;
  char  *args;
} argument_t;

argument_t *file_argument_new(const char *filename);
void        file_argument_free(argument_t *argument);
argument_t *argument_new(size_t argc, size_t len);
void        argument_free(argument_t *argument);
void        argument_fill(argument_t *argument, int argc, char *argv[]);

char *getProgname(char *string);
char *getOperator(const char *argument);
char *getOperatorName(const char *xoperator);

argument_t makeArgument(int argc, char *argv[]);
char *getFileArg(char *argument);

enum {START_DEC, START_JAN};
int get_season_start(void);
void get_season_name(const char *seas_name[]);
int month_to_season(int month);

void init_is_tty(void);

void progressInit(void);
void progressStatus(double offset, double refval, double curval);

int fileExists(const char *filename);
int userFileOverwrite(const char *filename);

/* convert a CDI datatype to string */
int datatype2str(int datatype, char *datatypestr);
int str2datatype(const char *datatypestr);

/* filename manipulation */
const char *filetypeext(int filetype);
void rm_filetypeext(char *file, const char *ext);
void repl_filetypeext(char file[], const char *oldext, const char *newext);


/* moved here from cdo.h */
void    cdiOpenError(int cdiErrno, const char *fmt, const char *path);
void    cdoAbort(const char *fmt, ...);
void    cdoWarning(const char *fmt, ...);
void    cdoPrint(const char *fmt, ...);

int  timer_new(const char *text);
void timer_report(void);
void timer_start(int it);
void timer_stop(int it);
double timer_val(int it);


void    operatorInputArg(const char *enter);
int     operatorArgc(void);
char  **operatorArgv(void);
void    operatorCheckArgc(int numargs);

const argument_t *cdoStreamName(int cnt);

void    cdoInitialize(void *argument);
void    cdoFinish(void);

int     cdoStreamNumber(void);
int     cdoStreamCnt(void);
int     cdoOperatorAdd(const char *name, int func, int intval, const char *enter);
int     cdoOperatorID(void);
int     cdoOperatorF1(int operID);
int     cdoOperatorF2(int operID);
const char *cdoOperatorName(int operID);
const char *cdoOperatorEnter(int operID);

int     cdoFiletype(void);

void    cdoInqHistory(int fileID);
void    cdoDefHistory(int fileID, char *histstring);

int     cdoDefineGrid(const char *gridfile);
int     cdoDefineZaxis(const char *zaxisfile);

int     vlistInqNWPV(int vlistID, int varID);
int     vlistIsSzipped(int vlistID);
int     vlist_check_gridsize(int vlistID);


void cdoGenFileSuffix(char *filesuffix, size_t maxlen, int filetype, int vlistID, const char *refname);

void writeNCgrid(const char *gridfile, int gridID, int *imask);
void defineZaxis(const char *zaxisarg);
void cdiDefTableID(int tableID);

int gridFromName(const char *gridname);
int zaxisFromName(const char *zaxisname);

/* refactor: moved here from cdo.h */
int cdo_omp_get_thread_num(void);
void strtolower(char *str);

/* refactor: moved here from cdo.c */
void exp_run(int argc, char *argv[], char *cdoExpName); // job.c
void printFeatures(void); // features.c
void printLibraries(void);  // features.c  

int wildcardmatch(const char *w, const char *s);

#endif  /* _UTIL_H */
