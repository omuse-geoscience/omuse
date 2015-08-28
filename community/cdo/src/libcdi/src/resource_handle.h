#ifndef RESOURCE_HANDLE_H
#define RESOURCE_HANDLE_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>

/*
 * CDI internal handling of resource handles given to user code
 */

/*
 * for reasons of compatibility with cfortran.h, the handle type is: int
 */
typedef int cdiResH;

/* return 0 on equality, not 0 otherwise */
typedef int    ( * valCompareFunc     )( void *, void * );
typedef void   ( * valDestroyFunc     )( void * );
typedef void   ( * valPrintFunc       )( void *, FILE * );
typedef int    ( * valGetPackSizeFunc )( void *, void *context );
typedef void   ( * valPackFunc        )( void *, void *buf, int size, int *pos, void *context );
typedef int    ( * valTxCodeFunc      )( void );

typedef struct {
  valCompareFunc     valCompare;
  valDestroyFunc     valDestroy;
  valPrintFunc       valPrint;
  valGetPackSizeFunc valGetPackSize;
  valPackFunc        valPack;
  valTxCodeFunc      valTxCode;
}resOps;

enum {
  RESH_IN_USE_BIT = 1 << 0,
  RESH_SYNC_BIT = 1 << 1,
  /* resource holds no value */
  RESH_UNUSED = 0,
  /* resource was deleted and needs to be synced */
  RESH_DESYNC_DELETED
    = RESH_SYNC_BIT,
  /* resource is synchronized */
  RESH_IN_USE
    = RESH_IN_USE_BIT,
  /* resource is in use, desynchronized and needs to be synced */
  RESH_DESYNC_IN_USE
    = RESH_IN_USE_BIT | RESH_SYNC_BIT,
};

void   reshListCreate(int namespaceID);
void   reshListDestruct(int namespaceID);
int    reshPut ( void *, const resOps * );
void reshReplace(cdiResH resH, void *p, const resOps *ops);
void   reshRemove ( cdiResH, const resOps * );
/*> doesn't check resource type */
void reshDestroy(cdiResH);

unsigned reshCountType(const resOps *resTypeOps);

void * reshGetValue(const char* caller, const char* expressionString, cdiResH id, const resOps* ops);
#define reshGetVal(resH, ops)  reshGetValue(__func__, #resH, resH, ops)

int reshEntryExists(cdiResH id, const resOps* ops);
#define reshExists(resH, ops)  reshEntryExists(resH, ops)

void reshGetResHListOfType(unsigned numIDs, int IDs[], const resOps *ops);

enum cdiApplyRet {
  CDI_APPLY_ERROR = -1,
  CDI_APPLY_STOP,
  CDI_APPLY_GO_ON,
};
enum cdiApplyRet
cdiResHApply(enum cdiApplyRet (*func)(int id, void *res, const resOps *p,
                                      void *data), void *data);
enum cdiApplyRet
cdiResHFilterApply(const resOps *p,
                   enum cdiApplyRet (*func)(int id, void *res,
                                            void *data),
                   void *data);

void   reshPackBufferCreate ( char **, int *, void *context );
void   reshPackBufferDestroy ( char ** );
int    reshResourceGetPackSize_intern(int resh, const resOps *ops, void *context, const char* caller, const char* expressionString);
#define reshResourceGetPackSize(resh, ops, context) reshResourceGetPackSize_intern(resh, ops, context, __func__, #resh)
void   reshPackResource_intern(int resh, const resOps *ops, void *buf, int buf_size, int *position, void *context, const char* caller, const char* expressionString);
#define reshPackResource(resh, ops, buf, buf_size, position, context) reshPackResource_intern(resh, ops, buf, buf_size, position, context, __func__, #resh)

void   reshSetStatus ( cdiResH, const resOps *, int );
int    reshGetStatus ( cdiResH, const resOps * );

void   reshLock   ( void );
void   reshUnlock ( void );

enum reshListMismatch {
  cdiResHListOccupationMismatch,
  cdiResHListResourceTypeMismatch,
  cdiResHListResourceContentMismatch,
};

int reshListCompare(int nsp0, int nsp1);
void reshListPrint(FILE *fp);

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
