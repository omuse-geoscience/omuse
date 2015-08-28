#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <limits.h>

#include "dmemory.h"
#include "cdi.h"
#include "cdi_int.h"
#include "resource_handle.h"
#include "resource_unpack.h"
#include "namespace.h"
#include "serialize.h"

#undef  UNDEFID
#define UNDEFID -1

int ECHAM4 = UNDEFID;
int ECHAM5 = UNDEFID;
int COSMO  = UNDEFID;

typedef struct
{
  int      self;
  int      used;
  int      instID;
  int      modelgribID;
  char    *name;
}
model_t;


static int  MODEL_Debug = 0;   /* If set to 1, debugging */

static void modelInit(void);


static int modelCompareP(void *modelptr1, void *modelptr2);
static void   modelDestroyP ( void * modelptr );
static void   modelPrintP   ( void * modelptr, FILE * fp );
static int    modelGetSizeP ( void * modelptr, void *context);
static void   modelPackP    ( void * modelptr, void * buff, int size,
                              int *position, void *context);
static int    modelTxCode   ( void );

static const resOps modelOps = {
  modelCompareP,
  modelDestroyP,
  modelPrintP,
  modelGetSizeP,
  modelPackP,
  modelTxCode
};

static
void modelDefaultValue ( model_t *modelptr )
{
  modelptr->self        = UNDEFID;
  modelptr->used        = 0;
  modelptr->instID      = UNDEFID;
  modelptr->modelgribID = UNDEFID;
  modelptr->name        = NULL;
}

static model_t *
modelNewEntry(cdiResH resH, int instID, int modelgribID, const char *name)
{
  model_t *modelptr;

  modelptr = (model_t *) xmalloc(sizeof(model_t));
  modelDefaultValue ( modelptr );
  if (resH == CDI_UNDEFID)
    modelptr->self = reshPut(modelptr, &modelOps);
  else
    {
      modelptr->self = resH;
      reshReplace(resH, modelptr, &modelOps);
    }
  modelptr->used = 1;
  modelptr->instID = instID;
  modelptr->modelgribID = modelgribID;
  if ( name && *name ) modelptr->name = strdupx(name);

  return (modelptr);
}

void modelDefaultEntries ( void )
{
  int instID, i;
  enum { nDefModels = 10 };
  cdiResH resH[nDefModels];

  instID  = institutInq(  0,   0, "ECMWF", NULL);
  /* (void)    modelDef(instID, 131, "ERA15"); */
  /* (void)    modelDef(instID, 199, "ERA40"); */
  instID  = institutInq(  0,   0, "MPIMET", NULL);

  resH[0] = ECHAM5  = modelDef(instID,  64, "ECHAM5.4");
  resH[1] = modelDef(instID,  63, "ECHAM5.3");
  resH[2] = modelDef(instID,  62, "ECHAM5.2");
  resH[3] = modelDef(instID,  61, "ECHAM5.1");

  instID  = institutInq( 98, 255, "MPIMET", NULL);
  resH[4] = modelDef(instID,  60, "ECHAM5.0");
  resH[5] = ECHAM4  = modelDef(instID,  50, "ECHAM4");
  resH[6] = modelDef(instID, 110, "MPIOM1");

  instID  = institutInq(  0,   0, "DWD", NULL);
  resH[7] = modelDef(instID, 149, "GME");

  instID  = institutInq(  0,   0, "MCH", NULL);
  //(void)  = modelDef(instID, 137, "COSMO");
  resH[8] = COSMO   = modelDef(instID, 255, "COSMO");

  instID  = institutInq(  0,   1, "NCEP", NULL);
  resH[9] = modelDef(instID,  80, "T62L28MRF");

  /* pre-defined models are not synchronized */
  for ( i = 0; i < nDefModels ; i++ )
    reshSetStatus(resH[i], &modelOps, RESH_IN_USE);
}

static
void modelInit(void)
{
  static int modelInitialized = 0;

  if (modelInitialized) return;

  modelInitialized = 1;
  char *env = getenv("MODEL_DEBUG");
  if ( env ) MODEL_Debug = atoi(env);
}

struct modelLoc
{
  char *name;
  int instID, modelgribID, resID;
};

static enum cdiApplyRet
findModelByID(int resID, void *res, void *data)
{
  model_t *modelptr = (model_t*) res;
  struct modelLoc *ret = (struct modelLoc*) data;
  int instID = ret->instID, modelgribID = ret->modelgribID;
  if (modelptr->used
      && modelptr->instID == instID
      && modelptr->modelgribID == modelgribID)
    {
      ret->resID = resID;
      return CDI_APPLY_STOP;
    }
  else
    return CDI_APPLY_GO_ON;
}

static enum cdiApplyRet
findModelByName(int resID, void *res, void *data)
{
  model_t *modelptr = (model_t*) res;
  struct modelLoc *ret = (struct modelLoc*) data;
  int instID = ret->instID, modelgribID = ret->modelgribID;
  const char *name = ret->name;
  if (modelptr->used
      && (instID == -1 || modelptr->instID == instID)
      && (modelgribID == 0 || modelptr->modelgribID == modelgribID)
      && modelptr->name)
    {
      const char *p = name, *q = modelptr->name;
      while (*p != '\0' && *p == *q)
        ++p, ++q;
      if (*p == '\0' || *q == '\0')
        {
          ret->resID = resID;
          return CDI_APPLY_STOP;
        }
    }
  return CDI_APPLY_GO_ON;
}

int modelInq(int instID, int modelgribID, char *name)
{
  modelInit ();

  struct modelLoc searchState = { .name = name, .instID = instID,
                                  .modelgribID = modelgribID,
                                  .resID = UNDEFID };
  if (name && *name)
    cdiResHFilterApply(&modelOps, findModelByName, &searchState);
  else
    cdiResHFilterApply(&modelOps, findModelByID, &searchState);
  return searchState.resID;
}


int modelDef(int instID, int modelgribID, const char *name)
{
  model_t *modelptr;

  modelInit ();

  modelptr = modelNewEntry(CDI_UNDEFID, instID, modelgribID, name);

  return modelptr->self;
}


int modelInqInstitut(int modelID)
{
  model_t *modelptr = NULL;

  modelInit ();

  if ( modelID != UNDEFID )
    modelptr = ( model_t * ) reshGetVal ( modelID, &modelOps );

  return modelptr ? modelptr->instID : UNDEFID;
}


int modelInqGribID(int modelID)
{
  model_t *modelptr = NULL;

  modelInit ();

  if ( modelID != UNDEFID )
    modelptr = ( model_t * ) reshGetVal ( modelID, &modelOps );

  return modelptr ? modelptr->modelgribID : UNDEFID;
}


const char *modelInqNamePtr(int modelID)
{
  model_t *modelptr = NULL;

  modelInit ();

  if ( modelID != UNDEFID )
    modelptr = ( model_t * ) reshGetVal ( modelID, &modelOps );

  return modelptr ? modelptr->name : NULL;
}


static int
modelCompareP(void *modelptr1, void *modelptr2)
{
  model_t *model1 = modelptr1, *model2 = modelptr2;
  int diff = (namespaceResHDecode(model1->instID).idx
              != namespaceResHDecode(model2->instID).idx)
    | (model1->modelgribID != model2->modelgribID)
    | (strcmp(model1->name, model2->name) != 0);
  return diff;
}


void modelDestroyP ( void * modelptr )
{
  model_t *mp = (model_t*) modelptr;
  if (mp->name)
    free(mp->name);
  free(mp);
}


void modelPrintP   ( void * modelptr, FILE * fp )
{
  model_t *mp = (model_t*) modelptr;

  if ( !mp ) return;

  fprintf ( fp, "#\n");
  fprintf ( fp, "# modelID %d\n", mp->self);
  fprintf ( fp, "#\n");
  fprintf ( fp, "self          = %d\n", mp->self );
  fprintf ( fp, "used          = %d\n", mp->used );
  fprintf ( fp, "instID        = %d\n", mp->instID );
  fprintf ( fp, "modelgribID   = %d\n", mp->modelgribID );
  fprintf ( fp, "name          = %s\n", mp->name ? mp->name : "NN" );
}


static int
modelTxCode ( void )
{
  return MODEL;
}

enum {
  model_nints = 4,
};


static int modelGetSizeP(void * modelptr, void *context)
{
  model_t *p = (model_t*)modelptr;
  size_t txsize = (size_t)serializeGetSize(model_nints, DATATYPE_INT, context)
    + (size_t)serializeGetSize(p->name?(int)strlen(p->name) + 1:0, DATATYPE_TXT, context);
  xassert(txsize <= INT_MAX);
  return (int)txsize;
}


static void modelPackP(void * modelptr, void * buf, int size, int *position, void *context)
{
  model_t *p = (model_t*) modelptr;
  int tempbuf[model_nints];
  tempbuf[0] = p->self;
  tempbuf[1] = p->instID;
  tempbuf[2] = p->modelgribID;
  tempbuf[3] = p->name ? (int)strlen(p->name) + 1 : 0;
  serializePack(tempbuf, model_nints, DATATYPE_INT, buf, size, position, context);
  if (p->name)
    serializePack(p->name, tempbuf[3], DATATYPE_TXT, buf, size, position, context);
}

int
modelUnpack(void *buf, int size, int *position, int originNamespace, void *context,
            int force_id)
{
  int tempbuf[model_nints];
  char *name;
  serializeUnpack(buf, size, position, tempbuf, model_nints, DATATYPE_INT, context);
  if (tempbuf[3] != 0)
    {
      name = (char *)xmalloc((size_t)tempbuf[3]);
      serializeUnpack(buf, size, position,
                      name, tempbuf[3], DATATYPE_TXT, context);
    }
  else
    {
      name = "";
    }
  int targetID = namespaceAdaptKey(tempbuf[0], originNamespace);
  model_t *mp = modelNewEntry(force_id?targetID:CDI_UNDEFID,
                              namespaceAdaptKey(tempbuf[1], originNamespace),
                              tempbuf[2], name);
  if (tempbuf[3] != 0)
    free(name);
  xassert(!force_id
          || (mp->self == namespaceAdaptKey(tempbuf[0], originNamespace)));
  return mp->self;
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
