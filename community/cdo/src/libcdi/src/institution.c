#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <assert.h>
#include <limits.h>

#include "dmemory.h"
#include "cdi.h"
#include "cdi_int.h"
#include "resource_handle.h"
#include "resource_handle.h"
#include "resource_unpack.h"
#include "namespace.h"
#include "serialize.h"
#include "institution.h"

#undef  UNDEFID
#define UNDEFID  -1

int ECMWF  = UNDEFID;
int MPIMET = UNDEFID;
int DWD    = UNDEFID;
int MCH    = UNDEFID;

typedef struct
{
  int    self;
  int    used;
  int    center;
  int    subcenter;
  char  *name;
  char  *longname;
}
institute_t;


static int instituteCompareKernel(institute_t *ip1, institute_t *ip2);
static void instituteDestroyP(institute_t *instituteptr);
static void   institutePrintP(institute_t *instituteptr, FILE * fp);
static int instituteGetPackSize(institute_t *instituteptr, void *context);
static void   institutePackP    ( void * instituteptr, void *buf, int size, int *position, void *context );
static int    instituteTxCode   ( void );

static const resOps instituteOps = {
  (int (*)(void *, void *))instituteCompareKernel,
  (void (*)(void *))instituteDestroyP,
  (void (*)(void *, FILE *))institutePrintP,
  (int (*)(void *, void *))instituteGetPackSize,
  institutePackP,
  instituteTxCode
};

static
void instituteDefaultValue ( institute_t * instituteptr )
{
  instituteptr->self       = UNDEFID;
  instituteptr->used       = 0;
  instituteptr->center     = UNDEFID;
  instituteptr->subcenter  = UNDEFID;
  instituteptr->name       = NULL;
  instituteptr->longname   = NULL;
}

void instituteDefaultEntries ( void )
{
  cdiResH resH[]
    = { ECMWF   = institutDef( 98,   0, "ECMWF",     "European Centre for Medium-Range Weather Forecasts"),
        MPIMET  = institutDef( 98, 232, "MPIMET",    "Max-Planck-Institute for Meteorology"),
        institutDef( 98, 255, "MPIMET",    "Max-Planck-Institute for Meteorology"),
        institutDef( 98, 232, "MPIMET",    "Max-Planck Institute for Meteorology"),
        institutDef( 78,   0, "DWD",       "Deutscher Wetterdienst"),
        institutDef( 78, 255, "DWD",       "Deutscher Wetterdienst"),
        MCH     = institutDef(215, 255, "MCH",       "MeteoSwiss"),
        institutDef(  7,   0, "NCEP",      "National Centers for Environmental Prediction"),
        institutDef(  7,   1, "NCEP",      "National Centers for Environmental Prediction"),
        institutDef( 60,   0, "NCAR",      "National Center for Atmospheric Research"),
        institutDef( 74,   0, "METOFFICE", "U.K. Met Office"),
        institutDef( 97,   0, "ESA",       "European Space Agency"),
        institutDef( 99,   0, "KNMI",      "Royal Netherlands Meteorological Institute"),
  };
  /*     (void) institutDef(  0,   0, "IPSL", "IPSL (Institut Pierre Simon Laplace, Paris, France)"); */

  size_t n = sizeof(resH)/sizeof(*resH);

  for (size_t i = 0; i < n ; i++ )
    reshSetStatus(resH[i], &instituteOps, RESH_IN_USE);
}


static int
instituteCompareKernel(institute_t *  ip1, institute_t * ip2)
{
  int differ = 0;
  size_t len1, len2;

  if ( ip1->name )
    {
      if ( ip1->center    > 0 && ip2->center    != ip1->center )    differ = 1;
      if ( ip1->subcenter > 0 && ip2->subcenter != ip1->subcenter ) differ = 1;

      if ( !differ )
        {
          if ( ip2->name )
            {
              len1 = strlen(ip1->name);
              len2 = strlen(ip2->name);
              if ( (len1 != len2) || memcmp(ip2->name, ip1->name, len2) ) differ = 1;
            }
        }
    }
  else if ( ip1->longname )
    {
      if ( ip2->longname )
        {
          len1 = strlen(ip1->longname);
          len2 = strlen(ip2->longname);
          if ( (len1 < len2) || memcmp(ip2->longname, ip1->longname, len2) ) differ = 1;
        }
    }
  else
    {
      if ( !( ip2->center    == ip1->center &&
              ip2->subcenter == ip1->subcenter )) differ = 1;
    }

  return differ;
}


struct instLoc
{
  institute_t *ip;
  int id;
};

static enum cdiApplyRet
findInstitute(int id, void *res, void *data)
{
  institute_t * ip1 = ((struct instLoc *)data)->ip;
  institute_t * ip2 = (institute_t*) res;
  if (ip2->used && !instituteCompareKernel(ip1, ip2))
    {
      ((struct instLoc *)data)->id = id;
      return CDI_APPLY_STOP;
    }
  else
    return CDI_APPLY_GO_ON;
}


int institutInq(int center, int subcenter, const char *name, const char *longname)
{
  institute_t * ip_ref = (institute_t *)xmalloc(sizeof (*ip_ref));
  ip_ref->self       = UNDEFID;
  ip_ref->used       = 0;
  ip_ref->center     = center;
  ip_ref->subcenter  = subcenter;
  ip_ref->name       = name && name[0] ? (char *)name : NULL;
  ip_ref->longname   = longname && longname[0] ? (char *)longname : NULL;

  struct instLoc state = { .ip = ip_ref, .id = UNDEFID };
  cdiResHFilterApply(&instituteOps, findInstitute, &state);

  free(ip_ref);

  return state.id;
}

static
institute_t *instituteNewEntry(cdiResH resH, int center, int subcenter,
                               const char *name, const char *longname)
{
  institute_t *instituteptr = (institute_t*) xmalloc(sizeof(institute_t));
  instituteDefaultValue(instituteptr);
  if (resH == CDI_UNDEFID)
    instituteptr->self = reshPut(instituteptr, &instituteOps);
  else
    {
      instituteptr->self = resH;
      reshReplace(resH, instituteptr, &instituteOps);
    }
  instituteptr->used = 1;
  instituteptr->center = center;
  instituteptr->subcenter = subcenter;
  if ( name && *name )
    instituteptr->name = strdupx(name);
  if (longname && *longname)
    instituteptr->longname = strdupx(longname);
  return  instituteptr;
}


int institutDef(int center, int subcenter, const char *name, const char *longname)
{
  institute_t * instituteptr
    = instituteNewEntry(CDI_UNDEFID, center, subcenter, name, longname);
  return instituteptr->self;
}


int institutInqCenter(int instID)
{
  institute_t * instituteptr = NULL;

  if ( instID != UNDEFID )
    instituteptr = ( institute_t * ) reshGetVal ( instID, &instituteOps );

  return  instituteptr ? instituteptr->center : UNDEFID;
}


int institutInqSubcenter(int instID)
{
  institute_t * instituteptr = NULL;

  if ( instID != UNDEFID )
    instituteptr = ( institute_t * ) reshGetVal ( instID, &instituteOps );

  return instituteptr ? instituteptr->subcenter: UNDEFID;
}


const char *institutInqNamePtr(int instID)
{
  institute_t * instituteptr = NULL;

  if ( instID != UNDEFID )
    instituteptr = ( institute_t * ) reshGetVal ( instID, &instituteOps );

  return instituteptr ? instituteptr->name : NULL;
}


const char *institutInqLongnamePtr(int instID)
{
  institute_t * instituteptr = NULL;

  if ( instID != UNDEFID )
    instituteptr = ( institute_t * ) reshGetVal ( instID, &instituteOps );

  return instituteptr ? instituteptr->longname : NULL;
}

static enum cdiApplyRet
activeInstitutes(int id, void *res, void *data)
{
  (void)id;
  if (res && ((institute_t *)res)->used)
    ++(*(int *)data);
  return CDI_APPLY_GO_ON;
}

int institutInqNumber(void)
{
  int instNum = 0;

  cdiResHFilterApply(&instituteOps, activeInstitutes, &instNum);
  return instNum;
}


static void
instituteDestroyP(institute_t *instituteptr)
{
  xassert(instituteptr);

  int instituteID = instituteptr->self;
  free(instituteptr->name);
  free(instituteptr->longname);
  reshRemove(instituteID, &instituteOps);
  free(instituteptr);
}


static void institutePrintP(institute_t *ip, FILE * fp )
{
  if (ip)
    fprintf(fp, "#\n"
            "# instituteID %d\n"
            "#\n"
            "self          = %d\n"
            "used          = %d\n"
            "center        = %d\n"
            "subcenter     = %d\n"
            "name          = %s\n"
            "longname      = %s\n",
            ip->self, ip->self, ip->used, ip->center, ip->subcenter,
            ip->name ? ip->name : "NN",
            ip->longname ? ip->longname : "NN");
}


static int
instituteTxCode ( void )
{
  return INSTITUTE;
}

enum {
  institute_nints = 5,
};

static int instituteGetPackSize(institute_t *ip, void *context)
{
  size_t namelen = strlen(ip->name), longnamelen = strlen(ip->longname);
  xassert(namelen < INT_MAX && longnamelen < INT_MAX);
  size_t txsize = (size_t)serializeGetSize(institute_nints, DATATYPE_INT, context)
    + (size_t)serializeGetSize((int)namelen + 1, DATATYPE_TXT, context)
    + (size_t)serializeGetSize((int)longnamelen + 1, DATATYPE_TXT, context);
  xassert(txsize <= INT_MAX);
  return (int)txsize;
}

static void institutePackP(void * instituteptr, void *buf, int size, int *position, void *context)
{
  institute_t *p = (institute_t*) instituteptr;
  int tempbuf[institute_nints];
  tempbuf[0] = p->self;
  tempbuf[1] = p->center;
  tempbuf[2] = p->subcenter;
  tempbuf[3] = (int)strlen(p->name) + 1;
  tempbuf[4] = (int)strlen(p->longname) + 1;
  serializePack(tempbuf, institute_nints, DATATYPE_INT, buf, size, position, context);
  serializePack(p->name, tempbuf[3], DATATYPE_TXT, buf, size, position, context);
  serializePack(p->longname, tempbuf[4], DATATYPE_TXT, buf, size, position, context);
}

int instituteUnpack(void *buf, int size, int *position, int originNamespace,
                    void *context, int force_id)
{
  int tempbuf[institute_nints];
  int instituteID;
  char *name, *longname;
  serializeUnpack(buf, size, position, tempbuf, institute_nints, DATATYPE_INT, context);
  name = (char *)xmalloc((size_t)tempbuf[3] + (size_t)tempbuf[4]);
  longname = name + tempbuf[3];
  serializeUnpack(buf, size, position, name, tempbuf[3], DATATYPE_TXT, context);
  serializeUnpack(buf, size, position, longname, tempbuf[4], DATATYPE_TXT, context);
  int targetID = namespaceAdaptKey(tempbuf[0], originNamespace);
  institute_t *ip = instituteNewEntry(force_id?targetID:CDI_UNDEFID,
                                      tempbuf[1], tempbuf[2], name, longname);
  instituteID = ip->self;
  xassert(!force_id || instituteID == targetID);
  free(name);
  return instituteID;
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
