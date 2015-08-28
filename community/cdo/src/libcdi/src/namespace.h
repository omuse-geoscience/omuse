#ifndef NAMESPACE_H
#define NAMESPACE_H

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif


typedef enum {
  STAGE_DEFINITION = 0,
  STAGE_TIMELOOP   = 1,
  STAGE_CLEANUP    = 2,
  STAGE_UNUSED     = 3,
} statusCode;

typedef struct {
  int idx;
  int nsp;
} namespaceTuple_t;

enum namespaceSwitch
{
  NSSWITCH_NO_SUCH_SWITCH = -1,
  NSSWITCH_ABORT,
  NSSWITCH_WARNING,
  NSSWITCH_SERIALIZE_GET_SIZE,
  NSSWITCH_SERIALIZE_PACK,
  NSSWITCH_SERIALIZE_UNPACK,
  NSSWITCH_FILE_OPEN,
  NSSWITCH_FILE_WRITE,
  NSSWITCH_FILE_CLOSE,
  NSSWITCH_STREAM_OPEN_BACKEND,
  NSSWITCH_STREAM_DEF_VLIST_,
  NSSWITCH_STREAM_WRITE_VAR_,
  NSSWITCH_STREAM_WRITE_VAR_CHUNK_,
  NSSWITCH_STREAM_WRITE_VAR_PART_,
  NSSWITCH_STREAM_WRITE_SCATTERED_VAR_PART_,
  NSSWITCH_STREAM_CLOSE_BACKEND,
  NSSWITCH_STREAM_DEF_TIMESTEP_,
  NSSWITCH_STREAM_SYNC,
#ifdef HAVE_LIBNETCDF
  NSSWITCH_NC__CREATE,
  NSSWITCH_CDF_DEF_VAR,
  NSSWITCH_CDF_DEF_TIMESTEP,
  NSSWITCH_CDF_STREAM_SETUP,
#endif
  NUM_NAMESPACE_SWITCH,
};

union namespaceSwitchValue
{
  void *data;
  void (*func)();
};

#define NSSW_FUNC(p) ((union namespaceSwitchValue){ .func = (void (*)())(p) })
#define NSSW_DATA(p) ((union namespaceSwitchValue){ .data = (void *)(p) })

int              namespaceNew();
void             namespaceDelete(int namespaceID);
void             namespaceCleanup      ( void );
int              namespaceGetNumber    ( void );
void namespaceSetActive(int namespaceID);
int              namespaceGetActive    ( void );
int              namespaceIdxEncode    ( namespaceTuple_t );
int              namespaceIdxEncode2   ( int, int );
namespaceTuple_t namespaceResHDecode   ( int );
int              namespaceAdaptKey     ( int originResH, int originNamespace);
int              namespaceAdaptKey2    ( int );
void             namespaceDefResStatus ( statusCode );
statusCode       namespaceInqResStatus ( void );
void namespaceSwitchSet(enum namespaceSwitch sw,
                        union namespaceSwitchValue value);
union namespaceSwitchValue namespaceSwitchGet(enum namespaceSwitch sw);


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
