#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 600 /* PTHREAD_MUTEX_RECURSIVE */
#endif

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#if defined (HAVE_EXECINFO_H)
#include <execinfo.h>
#endif

static
void show_stackframe()
{
#if defined HAVE_EXECINFO_H && defined backtrace_size_t && defined HAVE_BACKTRACE
  void *trace[16];
  backtrace_size_t trace_size = backtrace(trace, 16);
  char **messages = backtrace_symbols(trace, trace_size);

  fprintf(stderr, "[bt] Execution path:\n");
  if ( messages ) {
    for ( backtrace_size_t i = 0; i < trace_size; ++i )
      fprintf(stderr, "[bt] %s\n", messages[i]);
    free(messages);
  }
#endif
}

#include "dmemory.h"
#include "resource_handle.h"
#include "namespace.h"
#include "serialize.h"
#include "cdi.h"
#include "error.h"
#include "file.h"
#include "resource_unpack.h"
#include "institution.h"
#include "model.h"

enum { MIN_LIST_SIZE = 128 };

static void listInitialize(void);

typedef struct listElem {
  union
  {
    /* free-list management data */
    struct
    {
      int next, prev;
    } free;
    /* holding an actual value */
    struct
    {
      const resOps *ops;
      void         *val;//ptr
    } v;
  } res;
  int           status;
} listElem_t;

struct resHList_t
{
  int size, freeHead, hasDefaultRes;
  listElem_t *resources;
};

static struct resHList_t *resHList;

static int resHListSize = 0;

#if  defined  (HAVE_LIBPTHREAD)
#  include <pthread.h>

static pthread_once_t  listInitThread = PTHREAD_ONCE_INIT;
static pthread_mutex_t listMutex;

#  define LIST_LOCK()         pthread_mutex_lock(&listMutex)
#  define LIST_UNLOCK()       pthread_mutex_unlock(&listMutex)
#  define LIST_INIT(init0)         do {                         \
    pthread_once(&listInitThread, listInitialize);              \
    pthread_mutex_lock(&listMutex);                             \
    if ((init0) && (!resHList || !resHList[0].resources))       \
      reshListCreate(0);                                        \
    pthread_mutex_unlock(&listMutex);                           \
  } while (0)



#else

static int listInit = 0;

#  define LIST_LOCK()
#  define LIST_UNLOCK()
#  define LIST_INIT(init0)        do {                          \
  if ( !listInit )                                              \
    {                                                           \
      listInitialize();                                         \
      if ((init0) && (!resHList || !resHList[0].resources))     \
        reshListCreate(0);                                      \
      listInit = 1;                                             \
    }                                                           \
  } while(0)

#endif

/**************************************************************/

static void
listInitResources(int nsp)
{
  xassert(nsp < resHListSize && nsp >= 0);
  int size = resHList[nsp].size = MIN_LIST_SIZE;
  xassert(resHList[nsp].resources == NULL);
  resHList[nsp].resources = (listElem_t*) xcalloc(MIN_LIST_SIZE, sizeof(listElem_t));
  listElem_t *p = resHList[nsp].resources;

  for (int i = 0; i < size; i++ )
    {
      p[i].res.free.next = i + 1;
      p[i].res.free.prev = i - 1;
      p[i].status = RESH_UNUSED;
    }

  p[size-1].res.free.next = -1;
  resHList[nsp].freeHead = 0;
  int oldNsp = namespaceGetActive();
  namespaceSetActive(nsp);
  instituteDefaultEntries();
  modelDefaultEntries();
  namespaceSetActive(oldNsp);
}

static inline void
reshListClearEntry(int i)
{
  resHList[i].size = 0;
  resHList[i].resources = NULL;
  resHList[i].freeHead = -1;
}

void
reshListCreate(int namespaceID)
{
  LIST_INIT(namespaceID != 0);
  LIST_LOCK();
  if (resHListSize <= namespaceID)
    {
      resHList = (struct resHList_t *)xrealloc(resHList, (size_t)(namespaceID + 1) * sizeof (resHList[0]));
      for (int i = resHListSize; i <= namespaceID; ++i)
        reshListClearEntry(i);
      resHListSize = namespaceID + 1;
    }
  listInitResources(namespaceID);
  LIST_UNLOCK();
}


/**************************************************************/

void
reshListDestruct(int namespaceID)
{
  LIST_LOCK();
  xassert(resHList && namespaceID >= 0 && namespaceID < resHListSize);
  int callerNamespaceID = namespaceGetActive();
  namespaceSetActive(namespaceID);
  if (resHList[namespaceID].resources)
    {
      for ( int j = 0; j < resHList[namespaceID].size; j++ )
        {
          listElem_t *listElem = resHList[namespaceID].resources + j;
          if (listElem->status & RESH_IN_USE_BIT)
            listElem->res.v.ops->valDestroy(listElem->res.v.val);
        }
      free(resHList[namespaceID].resources);
      resHList[namespaceID].resources = NULL;
      reshListClearEntry(namespaceID);
    }
  if (resHList[callerNamespaceID].resources)
    namespaceSetActive(callerNamespaceID);
  LIST_UNLOCK();
}


static void listDestroy ( void )
{
  LIST_LOCK();
  for (int i = resHListSize; i > 0; --i)
    if (resHList[i-1].resources)
      namespaceDelete(i-1);
  resHListSize = 0;
  free(resHList);
  resHList = NULL;
  cdiReset();
  LIST_UNLOCK();
}

/**************************************************************/

static
void listInitialize ( void )
{
#if  defined  (HAVE_LIBPTHREAD)
  pthread_mutexattr_t ma;
  pthread_mutexattr_init(&ma);
  pthread_mutexattr_settype(&ma, PTHREAD_MUTEX_RECURSIVE);
  /* initialize global API mutex lock */
  pthread_mutex_init ( &listMutex, &ma);
  pthread_mutexattr_destroy(&ma);
#endif
  /* file is special and has its own table, which needs to be
   * created, before we register the listDestroy exit handler */
  {
    int null_id;
    null_id = fileOpen_serial("/dev/null", "r");
    if (null_id != -1)
      fileClose_serial(null_id);
  }
  atexit ( listDestroy );
}

/**************************************************************/

static
void listSizeExtend()
{
  int nsp = namespaceGetActive ();
  int oldSize = resHList[nsp].size;
  size_t newListSize = (size_t)oldSize + MIN_LIST_SIZE;

  resHList[nsp].resources = (listElem_t*) xrealloc(resHList[nsp].resources,
                                                   newListSize * sizeof(listElem_t));

  listElem_t *r = resHList[nsp].resources;
  for (size_t i = (size_t)oldSize; i < newListSize; ++i)
    {
      r[i].res.free.next = (int)i + 1;
      r[i].res.free.prev = (int)i - 1;
      r[i].status = RESH_UNUSED;
    }

  if (resHList[nsp].freeHead != -1)
    r[resHList[nsp].freeHead].res.free.prev = (int)newListSize - 1;
  r[newListSize-1].res.free.next = resHList[nsp].freeHead;
  r[oldSize].res.free.prev = -1;
  resHList[nsp].freeHead = oldSize;
  resHList[nsp].size = (int)newListSize;
}

/**************************************************************/

static void
reshPut_(int nsp, int entry, void *p, const resOps *ops)
{
  listElem_t *newListElem = resHList[nsp].resources + entry;
  int next = newListElem->res.free.next,
    prev = newListElem->res.free.prev;
  if (next != -1)
    resHList[nsp].resources[next].res.free.prev = prev;
  if (prev != -1)
    resHList[nsp].resources[prev].res.free.next = next;
  else
    resHList[nsp].freeHead = next;
  newListElem->res.v.val = p;
  newListElem->res.v.ops = ops;
  newListElem->status = RESH_DESYNC_IN_USE;
}

int reshPut ( void *p, const resOps *ops )
{
  xassert ( p && ops );

  LIST_INIT(1);

  LIST_LOCK();

  int nsp = namespaceGetActive ();

  if ( resHList[nsp].freeHead == -1) listSizeExtend();
  int entry = resHList[nsp].freeHead;
  cdiResH resH = namespaceIdxEncode2(nsp, entry);
  reshPut_(nsp, entry, p, ops);

  LIST_UNLOCK();

  return resH;
}

/**************************************************************/

static void
reshRemove_(int nsp, int idx)
{
  int curFree = resHList[nsp].freeHead;
  listElem_t *r = resHList[nsp].resources;
  r[idx].res.free.next = curFree;
  r[idx].res.free.prev = -1;
  if (curFree != -1)
    r[curFree].res.free.prev = idx;
  r[idx].status = RESH_DESYNC_DELETED;
  resHList[nsp].freeHead = idx;
}

void reshDestroy(cdiResH resH)
{
  int nsp;
  namespaceTuple_t nspT;

  LIST_LOCK();

  nsp = namespaceGetActive ();

  nspT = namespaceResHDecode ( resH );

  xassert ( nspT.nsp == nsp
            && nspT.idx >= 0
            && nspT.idx < resHList[nsp].size
            && resHList[nsp].resources[nspT.idx].res.v.ops);

  if (resHList[nsp].resources[nspT.idx].status & RESH_IN_USE_BIT)
    reshRemove_(nsp, nspT.idx);

  LIST_UNLOCK();
}

void reshRemove ( cdiResH resH, const resOps * ops )
{
  int nsp;
  namespaceTuple_t nspT;

  LIST_LOCK();

  nsp = namespaceGetActive ();

  nspT = namespaceResHDecode ( resH );

  xassert ( nspT.nsp == nsp
            && nspT.idx >= 0
            && nspT.idx < resHList[nsp].size
            && (resHList[nsp].resources[nspT.idx].status & RESH_IN_USE_BIT)
            && resHList[nsp].resources[nspT.idx].res.v.ops
            && resHList[nsp].resources[nspT.idx].res.v.ops == ops );

  reshRemove_(nsp, nspT.idx);

  LIST_UNLOCK();
}

/**************************************************************/

void reshReplace(cdiResH resH, void *p, const resOps *ops)
{
  xassert(p && ops);
  LIST_INIT(1);
  LIST_LOCK();
  int nsp = namespaceGetActive();
  namespaceTuple_t nspT = namespaceResHDecode(resH);
  while (resHList[nsp].size <= nspT.idx)
    listSizeExtend();
  listElem_t *q = resHList[nsp].resources + nspT.idx;
  if (q->status & RESH_IN_USE_BIT)
    {
      q->res.v.ops->valDestroy(q->res.v.val);
      reshRemove_(nsp, nspT.idx);
    }
  reshPut_(nsp, nspT.idx, p, ops);
  LIST_UNLOCK();
}


static listElem_t *
reshGetElem(const char *caller, const char* expressionString, cdiResH resH, const resOps *ops)
{
  listElem_t *listElem;
  int nsp;
  namespaceTuple_t nspT;
  xassert ( ops );

  LIST_INIT(1);

  LIST_LOCK();

  nsp = namespaceGetActive ();

  nspT = namespaceResHDecode ( resH );
  assert(nspT.idx >= 0);

  if (nspT.nsp == nsp &&
      nspT.idx < resHList[nsp].size)
    {
      listElem = resHList[nsp].resources + nspT.idx;
      LIST_UNLOCK();
    }
  else
    {
      LIST_UNLOCK();
      show_stackframe();

      if ( resH == CDI_UNDEFID )
        {
          xabortC(caller, "Error while trying to resolve the ID \"%s\" in `%s()`: the value is CDI_UNDEFID (= %d).\n\tThis is most likely the result of a failed earlier call. Please check the IDs returned by CDI.", expressionString, caller, resH);
        }
      else
        {
          xabortC(caller, "Error while trying to resolve the ID \"%s\" in `%s()`: the value is garbage (= %d, which resolves to namespace = %d, index = %d).\n\tThis is either the result of using an uninitialized variable,\n\tof using a value as an ID that is not an ID,\n\tor of using an ID after it has been invalidated.", expressionString, caller, resH, nspT.nsp, nspT.idx);
        }
    }

  if ( !(listElem && listElem->res.v.ops == ops) )
    {
      show_stackframe();

      xabortC(caller, "Error while trying to resolve the ID \"%s\" in `%s()`: list element not found. The failed ID is %d", expressionString, caller, (int)resH);
    }

  return listElem;
}

void *reshGetValue(const char * caller, const char* expressionString, cdiResH resH, const resOps * ops)
{
  return reshGetElem(caller, expressionString, resH, ops)->res.v.val;
}

int
reshEntryExists(cdiResH resH, const resOps * ops)
{
  int nsp;
  namespaceTuple_t nspT;
  listElem_t *listElem = NULL;
  xassert ( ops );
  LIST_INIT(1);
  LIST_LOCK();
  nsp = namespaceGetActive ();
  nspT = namespaceResHDecode ( resH );
  LIST_UNLOCK();

  int found = 0;
  if (nspT.nsp == nsp &&
      nspT.idx < resHList[nsp].size)
      listElem = resHList[nsp].resources + nspT.idx;

  LIST_UNLOCK();
  if ( listElem && (listElem->res.v.ops == ops) )  found = 1;

  return found;
}


/**************************************************************/

void reshGetResHListOfType(unsigned numIDs, int resHs[], const resOps *ops)
{
  xassert ( resHs && ops );

  LIST_INIT(1);

  LIST_LOCK();

  int nsp = namespaceGetActive();
  unsigned j = 0;
  for (int i = 0; i < resHList[nsp].size && j < numIDs; i++ )
    if ((resHList[nsp].resources[i].status & RESH_IN_USE_BIT)
        && resHList[nsp].resources[i].res.v.ops == ops)
      resHs[j++] = namespaceIdxEncode2(nsp, i);

  LIST_UNLOCK();
}

enum cdiApplyRet
cdiResHApply(enum cdiApplyRet (*func)(int id, void *res, const resOps *p,
                                      void *data), void *data)
{
  xassert(func);

  LIST_INIT(1);

  LIST_LOCK();

  int nsp = namespaceGetActive ();
  enum cdiApplyRet ret = CDI_APPLY_GO_ON;
  for (int i = 0; i < resHList[nsp].size && ret > 0; ++i)
    if (resHList[nsp].resources[i].status & RESH_IN_USE_BIT)
      ret = func(namespaceIdxEncode2(nsp, i),
                 resHList[nsp].resources[i].res.v.val,
                 resHList[nsp].resources[i].res.v.ops, data);
  LIST_UNLOCK();
  return ret;
}


enum cdiApplyRet
cdiResHFilterApply(const resOps *p,
                   enum cdiApplyRet (*func)(int id, void *res, void *data),
                   void *data)
{
  xassert(p && func);

  LIST_INIT(1);

  LIST_LOCK();

  int nsp = namespaceGetActive ();
  enum cdiApplyRet ret = CDI_APPLY_GO_ON;
  listElem_t *r = resHList[nsp].resources;
  for (int i = 0; i < resHList[nsp].size && ret > 0; ++i)
    if ((r[i].status & RESH_IN_USE_BIT) && r[i].res.v.ops == p)
      ret = func(namespaceIdxEncode2(nsp, i), r[i].res.v.val,
                 data);
  LIST_UNLOCK();
  return ret;
}




/**************************************************************/

unsigned reshCountType(const resOps *ops)
{
  unsigned countType = 0;

  xassert(ops);

  LIST_INIT(1);

  LIST_LOCK();

  int nsp = namespaceGetActive ();

  listElem_t *r = resHList[nsp].resources;
  size_t len = (size_t)resHList[nsp].size;
  for (size_t i = 0; i < len; i++ )
    countType += ((r[i].status & RESH_IN_USE_BIT) && r[i].res.v.ops == ops);

  LIST_UNLOCK();

  return countType;
}

/**************************************************************/

int
reshResourceGetPackSize_intern(int resH, const resOps *ops, void *context, const char* caller, const char* expressionString)
{
  listElem_t *curr = reshGetElem(caller, expressionString, resH, ops);
  return curr->res.v.ops->valGetPackSize(curr->res.v.val, context);
}

void
reshPackResource_intern(int resH, const resOps *ops, void *buf, int buf_size, int *position, void *context,
                        const char* caller, const char* expressionString)
{
  listElem_t *curr = reshGetElem(caller, expressionString, resH, ops);
  curr->res.v.ops->valPack(curr->res.v.val, buf, buf_size, position, context);
}

enum {
  resHPackHeaderNInt = 2,
  resHDeleteNInt = 2,
};

static int getPackBufferSize(void *context)
{
  int intpacksize, packBufferSize = 0;

  int nsp = namespaceGetActive ();

  /* pack start marker, namespace and sererator marker */
  packBufferSize += resHPackHeaderNInt * (intpacksize = serializeGetSize(1, DATATYPE_INT, context));

  /* pack resources, type marker and seperator marker */
  listElem_t *r = resHList[nsp].resources;
  for ( int i = 0; i < resHList[nsp].size; i++)
    if (r[i].status & RESH_SYNC_BIT)
      {
        if (r[i].status == RESH_DESYNC_DELETED)
          {
            packBufferSize += resHDeleteNInt * intpacksize;
          }
        else if (r[i].status == RESH_DESYNC_IN_USE)
          {
            xassert ( r[i].res.v.ops );
            /* packed resource plus 1 int for type */
            packBufferSize +=
              r[i].res.v.ops->valGetPackSize(r[i].res.v.val, context)
              + intpacksize;
          }
      }
  /* end marker */
  packBufferSize += intpacksize;

  return packBufferSize;
}

/**************************************************************/

void reshPackBufferDestroy ( char ** buffer )
{
  if ( buffer ) free ( *buffer );
}

/**************************************************************/

void reshPackBufferCreate(char **packBuffer, int *packBufferSize, void *context)
{
  int i, packBufferPos = 0;
  int end = END;

  xassert ( packBuffer );

  LIST_LOCK();

  int nsp = namespaceGetActive ();

  int pBSize = *packBufferSize = getPackBufferSize(context);
  char *pB = *packBuffer = (char *)xcalloc(1, (size_t)pBSize);

  {
    int header[resHPackHeaderNInt] = { START, nsp };
    serializePack(header, resHPackHeaderNInt,  DATATYPE_INT, pB, pBSize, &packBufferPos, context);
  }

  listElem_t *r = resHList[nsp].resources;
  for ( i = 0; i < resHList[nsp].size; i++ )
    if (r[i].status & RESH_SYNC_BIT)
      {
        if (r[i].status == RESH_DESYNC_DELETED)
          {
            int temp[resHDeleteNInt]
              = { RESH_DELETE, namespaceIdxEncode2(nsp, i) };
            serializePack(temp, resHDeleteNInt, DATATYPE_INT,
                          pB, pBSize, &packBufferPos, context);
          }
        else
          {
            listElem_t * curr = r + i;
            xassert ( curr->res.v.ops );
            int type = curr->res.v.ops->valTxCode();
            if ( ! type ) continue;
            serializePack(&type, 1, DATATYPE_INT, pB,
                          pBSize, &packBufferPos, context);
            curr->res.v.ops->valPack(curr->res.v.val,
                                     pB, pBSize, &packBufferPos, context);
          }
        r[i].status &= ~RESH_SYNC_BIT;
      }

  LIST_UNLOCK();

  serializePack(&end, 1,  DATATYPE_INT, pB, pBSize, &packBufferPos, context);
}

/**************************************************************/

/* for thread safety this feature would have to be integrated in reshPut */

void reshSetStatus ( cdiResH resH, const resOps * ops, int status )
{
  int nsp;
  namespaceTuple_t nspT;
  listElem_t * listElem;

  xassert(ops && (status & RESH_IN_USE_BIT));

  LIST_INIT(1);

  LIST_LOCK();

  nsp = namespaceGetActive ();

  nspT = namespaceResHDecode ( resH );

  xassert ( nspT.nsp == nsp &&
            nspT.idx >= 0 &&
            nspT.idx < resHList[nsp].size );

  xassert ( resHList[nsp].resources );
  listElem = resHList[nsp].resources + nspT.idx;

  xassert ( listElem->res.v.ops == ops );

  listElem->status = status;

  LIST_UNLOCK();
}

/**************************************************************/

int reshGetStatus ( cdiResH resH, const resOps * ops )
{
  int nsp;
  namespaceTuple_t nspT;

  xassert ( ops );

  LIST_INIT(1);

  LIST_LOCK();

  nsp = namespaceGetActive ();

  nspT = namespaceResHDecode ( resH );

  xassert ( nspT.nsp == nsp &&
            nspT.idx >= 0 &&
            nspT.idx < resHList[nsp].size );

  listElem_t *listElem = resHList[nsp].resources + nspT.idx;

  const resOps *elemOps = listElem->res.v.ops;

  LIST_UNLOCK();

  xassert(listElem && elemOps == ops);

  return listElem->status;
}

/**************************************************************/

void reshLock ()
{
  LIST_LOCK();
}

/**************************************************************/

void reshUnlock ()
{
  LIST_UNLOCK();
}

/**************************************************************/

int reshListCompare ( int nsp0, int nsp1 )
{
  LIST_INIT(1);
  LIST_LOCK();

  xassert(resHListSize > nsp0 && resHListSize > nsp1 &&
          nsp0 >= 0 && nsp1 >= 0);

  int valCompare = 0;
  int i, listSizeMin = (resHList[nsp0].size <= resHList[nsp1].size)
    ? resHList[nsp0].size : resHList[nsp1].size;
  listElem_t *resources0 = resHList[nsp0].resources,
    *resources1 = resHList[nsp1].resources;
  for (i = 0; i < listSizeMin; i++)
    {
      int occupied0 = (resources0[i].status & RESH_IN_USE_BIT) != 0,
        occupied1 = (resources1[i].status & RESH_IN_USE_BIT) != 0;
      /* occupation mismatch ? */
      int diff = occupied0 ^ occupied1;
      valCompare |= (diff << cdiResHListOccupationMismatch);
      if (!diff && occupied0)
        {
          /* both occupied, do resource types match? */
          diff = (resources0[i].res.v.ops != resources1[i].res.v.ops
                  || resources0[i].res.v.ops == NULL);
          valCompare |= (diff << cdiResHListResourceTypeMismatch);
          if (!diff)
            {
              /* types match, does content match also? */
              diff
                = resources0[i].res.v.ops->valCompare(resources0[i].res.v.val,
                                                      resources1[i].res.v.val);
              valCompare |= (diff << cdiResHListResourceContentMismatch);
            }
        }
    }
  /* find resources in nsp 0 beyond end of nsp 1 */
  for (int j = listSizeMin; j < resHList[nsp0].size; ++j)
    valCompare |= (((resources0[j].status & RESH_IN_USE_BIT) != 0)
                   << cdiResHListOccupationMismatch);
  /* find resources in nsp 1 beyond end of nsp 0 */
  for (; i < resHList[nsp1].size; ++i)
    valCompare |= (((resources1[i].status & RESH_IN_USE_BIT) != 0)
                   << cdiResHListOccupationMismatch);

  LIST_UNLOCK();

  return valCompare;
}

/**************************************************************/

void reshListPrint(FILE *fp)
{
  int i, j, temp;
  listElem_t * curr;

  LIST_INIT(1);


  temp = namespaceGetActive ();

  fprintf ( fp, "\n\n##########################################\n#\n#  print " \
            "global resource list \n#\n" );

  for ( i = 0; i < namespaceGetNumber (); i++ )
    {
      namespaceSetActive ( i );

      fprintf ( fp, "\n" );
      fprintf ( fp, "##################################\n" );
      fprintf ( fp, "#\n" );
      fprintf ( fp, "# namespace=%d\n", i );
      fprintf ( fp, "#\n" );
      fprintf ( fp, "##################################\n\n" );

      fprintf ( fp, "resHList[%d].size=%d\n", i, resHList[i].size );

      for ( j = 0; j < resHList[i].size; j++ )
        {
          curr = resHList[i].resources + j;
          if (!(curr->status & RESH_IN_USE_BIT))
            {
              curr->res.v.ops->valPrint(curr->res.v.val, fp);
              fprintf ( fp, "\n" );
            }
        }
    }
  fprintf ( fp, "#\n#  end global resource list" \
            "\n#\n##########################################\n\n" );

  namespaceSetActive ( temp );
}


/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
