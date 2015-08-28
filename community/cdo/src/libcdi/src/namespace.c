#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 600 /* PTHREAD_MUTEX_RECURSIVE */
#endif

#include <limits.h>
#include <stdlib.h>
#include <stdio.h>

#include "cdi.h"
#include "dmemory.h"
#include "namespace.h"
#include "resource_handle.h"
#include "serialize.h"
#include "error.h"
#include "cdf_int.h"
#include "file.h"
#include "cdi_int.h"
#include "stream_cdf.h"

static unsigned nNamespaces = 1;
static int activeNamespace = 0;

#ifdef HAVE_LIBNETCDF
#define CDI_NETCDF_SWITCHES                     \
  { .func = (void (*)()) nc__create },          \
  { .func = (void (*)()) cdf_def_var_serial },  \
  { .func = (void (*)()) cdfDefTimestep },      \
  { .func = (void (*)()) cdfDefVars }

#else
#define CDI_NETCDF_SWITCHES
#endif

#define defaultSwitches {                                   \
    { .func = (void (*)()) cdiAbortC_serial },              \
    { .func = (void (*)()) cdiWarning },                    \
    { .func = (void (*)()) serializeGetSizeInCore },        \
    { .func = (void (*)()) serializePackInCore },           \
    { .func = (void (*)()) serializeUnpackInCore },         \
    { .func = (void (*)()) fileOpen_serial },               \
    { .func = (void (*)()) fileWrite },                     \
    { .func = (void (*)()) fileClose_serial },              \
    { .func = (void (*)()) cdiStreamOpenDefaultDelegate },  \
    { .func = (void (*)()) cdiStreamDefVlist_ },            \
    { .func = (void (*)()) cdiStreamWriteVar_ },            \
    { .func = (void (*)()) cdiStreamwriteVarChunk_ },       \
    { .func = (void (*)()) 0 },                             \
    { .func = (void (*)()) 0 },                             \
    { .func = (void (*)()) cdiStreamCloseDefaultDelegate }, \
    { .func = (void (*)()) cdiStreamDefTimestep_ }, \
    { .func = (void (*)()) cdiStreamSync_ },                \
    CDI_NETCDF_SWITCHES                        \
    }

#if defined (SX)
static const union namespaceSwitchValue
  defaultSwitches_[NUM_NAMESPACE_SWITCH] = defaultSwitches;
#endif

static struct Namespace
{
  statusCode resStage;
  union namespaceSwitchValue switches[NUM_NAMESPACE_SWITCH];
} initialNamespace = {
  .resStage = STAGE_DEFINITION,
  .switches = defaultSwitches
};

static struct Namespace *namespaces = &initialNamespace;

static unsigned namespacesSize = 1;

#if  defined  (HAVE_LIBPTHREAD)
#  include <pthread.h>

static pthread_once_t  namespaceOnce = PTHREAD_ONCE_INIT;
static pthread_mutex_t namespaceMutex;

static void
namespaceInitialize(void)
{
  pthread_mutexattr_t ma;
  pthread_mutexattr_init(&ma);
  pthread_mutexattr_settype(&ma, PTHREAD_MUTEX_RECURSIVE);
  pthread_mutex_init(&namespaceMutex, &ma);
  pthread_mutexattr_destroy(&ma);
}

#  define NAMESPACE_LOCK()         pthread_mutex_lock(&namespaceMutex)
#  define NAMESPACE_UNLOCK()       pthread_mutex_unlock(&namespaceMutex)
#  define NAMESPACE_INIT()         pthread_once(&namespaceOnce, \
                                                namespaceInitialize)


#else

#  define NAMESPACE_INIT() do { } while (0)
#  define NAMESPACE_LOCK()
#  define NAMESPACE_UNLOCK()

#endif


enum {
  intbits = sizeof(int) * CHAR_BIT,
  nspbits = 4,
  idxbits = intbits - nspbits,
  nspmask = (( 1 << nspbits ) - 1) << idxbits,
  idxmask = ( 1 << idxbits ) - 1,
};

enum {
  NUM_NAMESPACES = 1 << nspbits,
  NUM_IDX = 1 << idxbits,
};


int namespaceIdxEncode ( namespaceTuple_t tin )
{
  xassert ( tin.nsp < NUM_NAMESPACES && tin.idx < NUM_IDX);
  return ( tin.nsp << idxbits ) + tin.idx;
}

int namespaceIdxEncode2 ( int nsp, int idx )
{
  xassert(nsp < NUM_NAMESPACES && idx < NUM_IDX);
  return ( nsp << idxbits ) + idx;
}


namespaceTuple_t namespaceResHDecode ( int resH )
{
  namespaceTuple_t tin;

  tin.idx = resH & idxmask;
  tin.nsp = (int)(((unsigned)( resH & nspmask )) >> idxbits);

  return tin;
}

int
namespaceNew()
{
  int newNamespaceID = -1;
  NAMESPACE_INIT();
  NAMESPACE_LOCK();
  if (namespacesSize > nNamespaces)
    {
      /* namespace is already available and only needs reinitialization */
      for (unsigned i = 0; i < namespacesSize; ++i)
        if (namespaces[i].resStage == STAGE_UNUSED)
          {
            newNamespaceID = (int)i;
            break;
          }
    }
  else if (namespacesSize == 1)
    {
      /* make room for additional namespace */
      struct Namespace *newNameSpaces
        = (struct Namespace *)xmalloc(((size_t)namespacesSize + 1) * sizeof (namespaces[0]));
      memcpy(newNameSpaces, namespaces, sizeof (namespaces[0]));
      namespaces = newNameSpaces;
      ++namespacesSize;
      newNamespaceID = 1;
    }
  else if (namespacesSize < NUM_NAMESPACES)
    {
      /* make room for additional namespace */
      newNamespaceID = (int)namespacesSize;
      namespaces
        = (struct Namespace *)xrealloc(namespaces, ((size_t)namespacesSize + 1) * sizeof (namespaces[0]));
      ++namespacesSize;
    }
  else /* implicit: namespacesSize >= NUM_NAMESPACES */
    {
      NAMESPACE_UNLOCK();
      return -1;
    }
  xassert(newNamespaceID >= 0 && newNamespaceID < NUM_NAMESPACES);
  ++nNamespaces;
  namespaces[newNamespaceID].resStage = STAGE_DEFINITION;
#if defined (SX)
  memcpy(namespaces[newNamespaceID].switches,
         defaultSwitches_,
         sizeof (namespaces[newNamespaceID].switches));
#else
  memcpy(namespaces[newNamespaceID].switches,
         (union namespaceSwitchValue[NUM_NAMESPACE_SWITCH])defaultSwitches,
         sizeof (namespaces[newNamespaceID].switches));
#endif
  reshListCreate(newNamespaceID);
  NAMESPACE_UNLOCK();
  return newNamespaceID;
}

void
namespaceDelete(int namespaceID)
{
  NAMESPACE_INIT();
  NAMESPACE_LOCK();
  xassert(namespaceID >= 0 && (unsigned)namespaceID < namespacesSize
          && nNamespaces);
  reshListDestruct(namespaceID);
  namespaces[namespaceID].resStage = STAGE_UNUSED;
  --nNamespaces;
  NAMESPACE_UNLOCK();
}

int namespaceGetNumber ()
{
  return (int)nNamespaces;
}


void namespaceSetActive ( int nId )
{
  xassert((unsigned)nId < namespacesSize
          && namespaces[nId].resStage != STAGE_UNUSED);
  activeNamespace = nId;
}


int namespaceGetActive ()
{
  return activeNamespace;
}

int namespaceAdaptKey ( int originResH, int originNamespace )
{
  namespaceTuple_t tin;
  int nsp;

  if ( originResH == CDI_UNDEFID ) return CDI_UNDEFID;

  tin.idx = originResH & idxmask;
  tin.nsp = (int)(((unsigned)( originResH & nspmask )) >> idxbits);

  xassert ( tin.nsp == originNamespace );

  nsp = namespaceGetActive ();

  return namespaceIdxEncode2 ( nsp, tin.idx );
}


int namespaceAdaptKey2 ( int originResH )
{
  namespaceTuple_t tin;
  int nsp;

  if ( originResH == CDI_UNDEFID ) return CDI_UNDEFID;

  tin.idx = originResH & idxmask;
  tin.nsp = (int)(((unsigned)( originResH & nspmask )) >> idxbits);

  nsp = namespaceGetActive ();

  return namespaceIdxEncode2 ( nsp, tin.idx );
}


void namespaceDefResStatus ( statusCode argResStatus )
{
  int nsp = namespaceGetActive ();
  namespaces[nsp].resStage = argResStatus;
}


statusCode namespaceInqResStatus ( void )
{
  int nsp = namespaceGetActive ();
  return namespaces[nsp].resStage;
}

void namespaceSwitchSet(enum namespaceSwitch sw, union namespaceSwitchValue value)
{
  xassert(sw > NSSWITCH_NO_SUCH_SWITCH && sw < NUM_NAMESPACE_SWITCH);
  int nsp = namespaceGetActive();
  namespaces[nsp].switches[sw] = value;
}

union namespaceSwitchValue namespaceSwitchGet(enum namespaceSwitch sw)
{
  xassert(sw > NSSWITCH_NO_SUCH_SWITCH && sw < NUM_NAMESPACE_SWITCH);
  int nsp = namespaceGetActive();
  return namespaces[nsp].switches[sw];
}

void cdiReset(void)
{
  NAMESPACE_INIT();
  NAMESPACE_LOCK();
  for (unsigned namespaceID = 0; namespaceID < namespacesSize; ++namespaceID)
    if (namespaces[namespaceID].resStage != STAGE_UNUSED)
      namespaceDelete((int)namespaceID);
  if (namespaces != &initialNamespace)
    {
      free(namespaces);
      namespaces = &initialNamespace;
    }
  namespacesSize = 1;
  nNamespaces = 1;
  activeNamespace = 0;
  NAMESPACE_UNLOCK();
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
