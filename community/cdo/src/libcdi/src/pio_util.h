#ifndef PIO_UTIL_
#define PIO_UTIL_

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <mpi.h>

#ifndef _ERROR_H
#include "error.h"
#endif

#define MAXDEBUG           3

#define ddebug             0

#define debugString "#####"

void
cdiAbortC_MPI(const char * caller, const char * filename,
              const char *functionname, int line,
              const char * errorString, va_list ap)
  __attribute__((noreturn));

void cdiPioWarning(const char *caller, const char *fmt, va_list ap);

static inline int
callsToMPIAreAllowed()
{
  int init_flag = 0, finished_flag = 0;
  return MPI_Initialized(&init_flag) == MPI_SUCCESS && init_flag
    && MPI_Finalized(&finished_flag) == MPI_SUCCESS && !finished_flag;
}

static inline int
getMPICommWorldRank()
{
  int rank = -1;
  if (callsToMPIAreAllowed())
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank;
}

#define xdebug(fmt, ...)                                                \
  if ( ddebug ){                                                        \
    int rank = getMPICommWorldRank();                                   \
    fprintf ( stderr, "%s pe%d in %s, %s, line %d: " fmt "\n",          \
              debugString, rank,  __func__, __FILE__,  __LINE__,        \
              __VA_ARGS__ );                                            \
  }

#define xdebug3(fmt, ...)                                               \
  if ( ddebug == MAXDEBUG ){                                            \
    int rank = getMPICommWorldRank();                                   \
    fprintf ( stderr, "pe%d in %s, %s, line %d: " fmt "\n",             \
              rank,  __func__, __FILE__,  __LINE__,                     \
              __VA_ARGS__ );                                            \
  }

void pcdiXMPI(int iret, const char *, int);
#define xmpi(ret) do {                                  \
    int tmpIRet = (ret);                                   \
    if (tmpIRet != MPI_SUCCESS)                            \
      pcdiXMPI(tmpIRet, __FILE__, __LINE__ );              \
  } while(0)

void pcdiXMPIStat ( int, const char *, int, MPI_Status * );
#define xmpiStat(ret,stat) pcdiXMPIStat ( ret, __FILE__, __LINE__, stat )

void pcdiDebugComm ( const char *filename, const char *functionname, int line, \
                     MPI_Comm *comm );
#define xdebugComm(comm)\
  if ( ddebug ) pcdiDebugComm (  __FILE__, __func__, __LINE__, comm )

void pcdiDebugMsg ( const char * cdiDebugString, const char *filename, const char *functionname, int line, \
                    int tag, int source, int nfinished );
#define xdebugMsg(tag,source,nfinished) \
  if ( ddebug ) \
      pcdiDebugMsg ( debugString, __FILE__, __func__, __LINE__, tag, source, nfinished )

void pcdiDebugMsg2 ( const char *filename, const char *functionname, int line, \
                    int tag, int source, char * text );
#define xdebugMsg2(tag,source,text) \
  if ( ddebug ) pcdiDebugMsg ( __FILE__, __func__,  __LINE__, tag, source, text )

static inline int
sum_int(size_t n, int *a)
{
  int sum = 0;
  for (size_t i = 0; i < n; ++i)
    sum += a[i];
  return sum;
}


void printArray ( const char *, char *, const void *, int, int, const char *, const char *, int );
#define xprintArray(ps,array,n,datatype)                                \
  if ( ddebug )                                                         \
      printArray ( debugString, ps, array, n, datatype,  __func__, __FILE__, __LINE__ )
 
#define xprintArray3(ps,array,n,datatype)                                \
  if ( ddebug == MAXDEBUG )                                                         \
      printArray ( debugString, ps, array, n, datatype,  __func__, __FILE__, __LINE__ )

/**
 * @return number of dimensions
 */
int
cdiPioQueryVarDims(int varShape[3], int vlistID, int varID);

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
