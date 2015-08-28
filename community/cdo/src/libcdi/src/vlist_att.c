#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>

#include "dmemory.h"

#include "cdi.h"
#include "cdi_int.h"
#include "vlist.h"
#include "error.h"
#include "serialize.h"

static
cdi_atts_t *get_attsp(vlist_t *vlistptr, int varID)
{
  cdi_atts_t *attsp = NULL;

  if ( varID == CDI_GLOBAL )
    {
      attsp = &vlistptr->atts;
    }
  else
    {
      if ( varID >= 0 && varID < vlistptr->nvars )
	attsp = &(vlistptr->vars[varID].atts);
    }

  return (attsp);
}

static
cdi_att_t *find_att(cdi_atts_t *attsp, const char *name)
{
  xassert(attsp != NULL);

  if ( attsp->nelems == 0 ) return NULL;

  size_t slen = strlen(name);
  if ( slen > CDI_MAX_NAME ) slen = CDI_MAX_NAME;

  cdi_att_t *atts = attsp->value;
  for ( size_t attid = 0; attid < attsp->nelems; attid++ )
    {
      cdi_att_t *attp = atts + attid;
      if ( attp->namesz == slen && memcmp(attp->name, name, slen) == 0 )
        return (attp); /* Normal return */
    }

  return (NULL);
}

static
cdi_att_t *new_att(cdi_atts_t *attsp, const char *name)
{
  cdi_att_t *attp;
  size_t slen;

  xassert(attsp != NULL);
  xassert(name  != NULL);

  if ( attsp->nelems == attsp->nalloc ) return (NULL);

  attp = &(attsp->value[attsp->nelems]);
  attsp->nelems++;

  slen = strlen(name);
  if ( slen > CDI_MAX_NAME ) slen = CDI_MAX_NAME;

  attp->name = (char *) malloc(slen+1);
  memcpy(attp->name, name, slen+1);
  attp->namesz = slen;
  attp->xvalue = NULL;

  return (attp);
}

static
void fill_att(cdi_att_t *attp, int indtype, int exdtype, size_t nelems, size_t xsz, const void *xvalue)
{
  xassert(attp != NULL);

  attp->xsz = xsz;
  attp->indtype = indtype;
  attp->exdtype = exdtype;
  attp->nelems  = nelems;

  if ( xsz > 0 )
    {
      attp->xvalue = xrealloc(attp->xvalue, xsz);
      memcpy(attp->xvalue, xvalue, xsz);
    }
}

/*
@Function  vlistInqNatts
@Title     Get number of variable attributes

@Prototype int vlistInqNatts(int vlistID, int varID, int *nattsp)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate} or @fref{streamInqVlist}.
    @Item  varID    Variable identifier, or @func{CDI_GLOBAL} for a global attribute.
    @Item  nattsp   Pointer to location for returned number of variable attributes.

@Description
The function @func{vlistInqNatts} gets the number of variable attributes assigned to this variable.

@EndFunction
*/
int vlistInqNatts(int vlistID, int varID, int *nattsp)
{
  int status = CDI_NOERR;
  vlist_t *vlistptr;
  cdi_atts_t *attsp;

  vlistptr = vlist_to_pointer(vlistID);

  attsp = get_attsp(vlistptr, varID);
  xassert(attsp != NULL);

  *nattsp = (int)attsp->nelems;

  return (status);
}

/*
@Function  vlistInqAtt
@Title     Get information about an attribute

@Prototype int vlistInqAtt(int vlistID, int varID, int attnum, char *name, int *typep, int *lenp)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate} or @fref{streamInqVlist}.
    @Item  varID    Variable identifier, or @func{CDI_GLOBAL} for a global attribute.
    @Item  attnum   Attribute number (from 0 to natts-1).
    @Item  name     Pointer to the location for the returned attribute name. The caller must allocate space for the
                    returned string. The maximum possible length, in characters, of
                    the string is given by the predefined constant @func{CDI_MAX_NAME}.
    @Item  typep    Pointer to location for returned attribute type.
    @Item  lenp     Pointer to location for returned attribute number.

@Description
The function @func{vlistInqAtt} gets information about an attribute.

@EndFunction
*/
int vlistInqAtt(int vlistID, int varID, int attnum, char *name, int *typep, int *lenp)
{
  int status = CDI_NOERR;
  cdi_att_t *attp = NULL;

  xassert(name != NULL);

  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  cdi_atts_t *attsp = get_attsp(vlistptr, varID);
  xassert(attsp != NULL);

  if ( attnum >= 0 && attnum < (int)attsp->nelems )
    attp = &(attsp->value[attnum]);

  if ( attp != NULL ) /* name in use */
    {
      memcpy(name, attp->name, attp->namesz+1);
      *typep  = attp->exdtype;
      *lenp   = (int)attp->nelems;
    }
  else
    {
      name[0] =  0;
      *typep  = -1;
      *lenp   =  0;
      status  = -1;
    }

  return (status);
}


int vlistDelAtts(int vlistID, int varID)
{
  int status = CDI_NOERR;
  vlist_t *vlistptr;
  cdi_att_t *attp = NULL;
  cdi_atts_t *attsp;
  int attid;

  vlistptr = vlist_to_pointer(vlistID);

  attsp = get_attsp(vlistptr, varID);
  xassert(attsp != NULL);

  for ( attid = 0; attid < (int)attsp->nelems; attid++ )
    {
      attp = &(attsp->value[attid]);
      if ( attp->name   ) free(attp->name);
      if ( attp->xvalue ) free(attp->xvalue);
    }

  attsp->nelems = 0;

  return (status);
}


int vlistDelAtt(int vlistID, int varID, const char *name)
{
  int status = CDI_NOERR;

  UNUSED(vlistID);
  UNUSED(varID);
  UNUSED(name);

  fprintf(stderr, "vlistDelAtt not implemented!\n");

  return (status);
}

static
int vlist_def_att(int indtype, int exdtype, int vlistID, int varID, const char *name, size_t len, size_t xsz, const void *xp)
{
  int status = CDI_NOERR;
  vlist_t *vlistptr;
  cdi_att_t *attp;
  cdi_atts_t *attsp;

  if ( len != 0 && xp == NULL ) /* Null arg */
    {
      return (CDI_EINVAL);
    }

  vlistptr = vlist_to_pointer(vlistID);

  attsp = get_attsp(vlistptr, varID);
  xassert(attsp != NULL);

  attp = find_att(attsp, name);
  if ( attp == NULL )
    attp = new_att(attsp, name);

  if ( attp != NULL )
    fill_att(attp, indtype, exdtype, len, xsz, xp);

  return (status);
}

static
int vlist_inq_att(int indtype, int vlistID, int varID, const char *name, size_t mxsz, void *xp)
{
  int status = CDI_NOERR;
  vlist_t *vlistptr;
  cdi_att_t *attp;
  cdi_atts_t *attsp;
  size_t xsz;

  if ( mxsz != 0 && xp == NULL ) /* Null arg */
    {
      return (CDI_EINVAL);
    }

  vlistptr = vlist_to_pointer(vlistID);

  attsp = get_attsp(vlistptr, varID);
  xassert(attsp != NULL);

  attp = find_att(attsp, name);
  if ( attp != NULL ) /* name in use */
    {
      if ( attp->indtype == indtype )
	{
	  xsz = attp->xsz;
	  if ( mxsz < xsz ) xsz = mxsz;
	  if ( xsz > 0 )
	    memcpy(xp, attp->xvalue, xsz);
	}
      else
	{
	  Warning("Attribute %s has wrong data type!", name);
          status = -2;
	}
    }
  else
    {
      //Warning("Internal problem, attribute %s not found!", name);
      status = -1;
    }

  return (status);
}


int vlistCopyVarAtts(int vlistID1, int varID_1, int vlistID2, int varID_2)
{
  int status = CDI_NOERR;
  vlist_t *vlistptr1;
  cdi_att_t *attp = NULL;
  cdi_atts_t *attsp1;
  int attid;

  vlistptr1 = vlist_to_pointer(vlistID1);

  attsp1 = get_attsp(vlistptr1, varID_1);
  xassert(attsp1 != NULL);

  for ( attid = 0; attid < (int)attsp1->nelems; attid++ )
    {
      attp = &(attsp1->value[attid]);
      vlist_def_att(attp->indtype, attp->exdtype, vlistID2, varID_2, attp->name, attp->nelems, attp->xsz, attp->xvalue);
    }

  return (status);
}

/*
@Function  vlistDefAttInt
@Title     Define an integer attribute

@Prototype int vlistDefAttInt(int vlistID, int varID, const char *name, int type, int len, const int *ip)

@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  varID    Variable identifier, or @func{CDI_GLOBAL} for a global attribute.
    @Item  name     Attribute name.
    @Item  type     External data type (@func{DATATYPE_INT16} or @func{DATATYPE_INT32}).
    @Item  len      Number of values provided for the attribute.
    @Item  ip       Pointer to one or more integer values.

@Description
The function @func{vlistDefAttInt} defines an integer attribute.

@EndFunction
*/
int vlistDefAttInt(int vlistID, int varID, const char *name, int type, int len, const int *ip)
{
  return vlist_def_att(DATATYPE_INT, type, vlistID, varID, name, (size_t)len, (size_t)len * sizeof (int), ip);
}

/*
@Function  vlistDefAttFlt
@Title     Define a floating point attribute

@Prototype int vlistDefAttFlt(int vlistID, int varID, const char *name, int type, int len, const double *dp)

@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  varID    Variable identifier, or @func{CDI_GLOBAL} for a global attribute.
    @Item  name     Attribute name.
    @Item  type     External data type (@func{DATATYPE_FLT32} or @func{DATATYPE_FLT64}).
    @Item  len      Number of values provided for the attribute.
    @Item  dp       Pointer to one or more floating point values.

@Description
The function @func{vlistDefAttFlt} defines a floating point attribute.

@EndFunction
*/
int vlistDefAttFlt(int vlistID, int varID, const char *name, int type, int len, const double *dp)
{
  return vlist_def_att(DATATYPE_FLT, type, vlistID, varID, name, (size_t)len, (size_t)len * sizeof (double), dp);
}

/*
@Function  vlistDefAttTxt
@Title     Define a text attribute

@Prototype int vlistDefAttTxt(int vlistID, int varID, const char *name, int len, const char *tp)

@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  varID    Variable identifier, or @func{CDI_GLOBAL} for a global attribute.
    @Item  name     Attribute name.
    @Item  len      Number of values provided for the attribute.
    @Item  tp       Pointer to one or more character values.

@Description
The function @func{vlistDefAttTxt} defines a text attribute.

@EndFunction
*/
int vlistDefAttTxt(int vlistID, int varID, const char *name, int len, const char *tp)
{
  return vlist_def_att(DATATYPE_TXT, DATATYPE_TXT, vlistID, varID, name, (size_t)len, (size_t)len, tp);
}

/*
@Function  vlistInqAttInt
@Title     Get the value(s) of an integer attribute

@Prototype int vlistInqAttInt(int vlistID, int varID, const char *name, int mlen, int *ip)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate} or @fref{streamInqVlist}.
    @Item  varID    Variable identifier, or @func{CDI_GLOBAL} for a global attribute.
    @Item  name     Attribute name.
    @Item  mlen     Number of allocated values provided for the attribute.
    @Item  ip       Pointer location for returned integer attribute value(s).

@Description
The function @func{vlistInqAttInt} gets the values(s) of an integer attribute.

@EndFunction
*/
int vlistInqAttInt(int vlistID, int varID, const char *name, int mlen, int *ip)
{
  return vlist_inq_att(DATATYPE_INT, vlistID, varID, name, (size_t)mlen * sizeof (int), ip);
}

/*
@Function  vlistInqAttFlt
@Title     Get the value(s) of a floating point attribute

@Prototype int vlistInqAttFlt(int vlistID, int varID, const char *name, int mlen, double *dp)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate} or @fref{streamInqVlist}.
    @Item  varID    Variable identifier, or @func{CDI_GLOBAL} for a global attribute.
    @Item  name     Attribute name.
    @Item  mlen     Number of allocated values provided for the attribute.
    @Item  dp       Pointer location for returned floating point attribute value(s).

@Description
The function @func{vlistInqAttFlt} gets the values(s) of a floating point attribute.

@EndFunction
*/
int vlistInqAttFlt(int vlistID, int varID, const char *name, int mlen, double *dp)
{
  return vlist_inq_att(DATATYPE_FLT, vlistID, varID, name, (size_t)mlen * sizeof (double), dp);
}

/*
@Function  vlistInqAttTxt
@Title     Get the value(s) of a text attribute

@Prototype int vlistInqAttTxt(int vlistID, int varID, const char *name, int mlen, char *tp)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate} or @fref{streamInqVlist}.
    @Item  varID    Variable identifier, or @func{CDI_GLOBAL} for a global attribute.
    @Item  name     Attribute name.
    @Item  mlen     Number of allocated values provided for the attribute.
    @Item  tp       Pointer location for returned text attribute value(s).

@Description
The function @func{vlistInqAttTxt} gets the values(s) of a text attribute.

@EndFunction
*/
int vlistInqAttTxt(int vlistID, int varID, const char *name, int mlen, char *tp)
{
  return vlist_inq_att(DATATYPE_TXT, vlistID, varID, name, (size_t)mlen * sizeof (char), tp);
}

enum {
  vlist_att_nints = 4,          /* namesz, exdtype, indtype, nelems */
};

static inline int
vlistAttTypeLookup(cdi_att_t *attp)
{
  int type;
  switch (attp->indtype)
  {
  case DATATYPE_FLT:
    type = DATATYPE_FLT64;
    break;
  case DATATYPE_INT:
  case DATATYPE_TXT:
    type = attp->indtype;
    break;
  default:
    xabort("Unknown datatype encountered in attribute %s: %d\n",
            attp->name, attp->indtype);
  }
  return type;
}


int vlist_att_compare(vlist_t *a, int varIDA, vlist_t *b, int varIDB,
                      int attnum)
{
  cdi_atts_t *attspa = get_attsp(a, varIDA),
    *attspb = get_attsp(b, varIDB);
  if (attspa == NULL && attspb == NULL)
    return 0;
  xassert(attnum >= 0 && attnum < (int)attspa->nelems
          && attnum < (int)attspb->nelems);
  cdi_att_t *attpa = attspa->value + attnum,
    *attpb = attspb->value + attnum;
  size_t len;
  if ((len = attpa->namesz) != attpb->namesz)
    return 1;
  int diff;
  if ((diff = memcmp(attpa->name, attpb->name, len)))
    return 1;
  if (attpa->indtype != attpb->indtype
      || attpa->exdtype != attpb->exdtype
      || attpa->nelems != attpb->nelems)
    return 1;
  return memcmp(attpa->xvalue, attpb->xvalue, attpa->xsz);
}


static int
vlistAttGetSize(vlist_t *vlistptr, int varID, int attnum, void *context)
{
  cdi_atts_t *attsp;
  cdi_att_t *attp;

  xassert(attsp = get_attsp(vlistptr, varID));
  xassert(attnum >= 0 && attnum < (int)attsp->nelems);
  attp = &(attsp->value[attnum]);
  int txsize = serializeGetSize(vlist_att_nints, DATATYPE_INT, context)
    + serializeGetSize((int)attp->namesz, DATATYPE_TXT, context);
  txsize += serializeGetSize((int)attp->nelems, vlistAttTypeLookup(attp), context);
  return txsize;
}

int
vlistAttsGetSize(vlist_t *p, int varID, void *context)
{
  cdi_atts_t *attsp = get_attsp(p, varID);
  int txsize = serializeGetSize(1, DATATYPE_INT, context);
  size_t numAtts = attsp->nelems;
  for (size_t i = 0; i < numAtts; ++i)
    txsize += vlistAttGetSize(p, varID, (int)i, context);
  return txsize;
}

static void
vlistAttPack(vlist_t *vlistptr, int varID, int attnum,
             void * buf, int size, int *position, void *context)
{
  cdi_atts_t *attsp;
  cdi_att_t *attp;
  int tempbuf[vlist_att_nints];

  xassert(attsp = get_attsp(vlistptr, varID));
  xassert(attnum >= 0 && attnum < (int)attsp->nelems);
  attp = &(attsp->value[attnum]);
  tempbuf[0] = (int)attp->namesz;
  tempbuf[1] = attp->exdtype;
  tempbuf[2] = attp->indtype;
  tempbuf[3] = (int)attp->nelems;
  serializePack(tempbuf, vlist_att_nints, DATATYPE_INT, buf, size, position, context);
  serializePack(attp->name, (int)attp->namesz, DATATYPE_TXT, buf, size, position, context);
  serializePack(attp->xvalue, (int)attp->nelems, vlistAttTypeLookup(attp),
                buf, size, position, context);
}

void
vlistAttsPack(vlist_t *p, int varID,
              void * buf, int size, int *position, void *context)
{
  cdi_atts_t *attsp = get_attsp(p, varID);
  size_t numAtts = attsp->nelems;
  int numAttsI = (int)numAtts;
  xassert(numAtts <= INT_MAX);
  serializePack(&numAttsI, 1, DATATYPE_INT, buf, size, position, context);
  for (size_t i = 0; i < numAtts; ++i)
    vlistAttPack(p, varID, (int)i, buf, size, position, context);
}

static void
vlistAttUnpack(int vlistID, int varID,
               void * buf, int size, int *position, void *context)
{
  int tempbuf[vlist_att_nints];

  serializeUnpack(buf, size, position,
                  tempbuf, vlist_att_nints, DATATYPE_INT, context);
  char *attName = (char *)xmalloc((size_t)tempbuf[0] + 1);
  serializeUnpack(buf, size, position, attName, tempbuf[0], DATATYPE_TXT, context);
  attName[tempbuf[0]] = '\0';
  int attVDt;
  size_t elemSize;
  switch (tempbuf[2])
  {
  case DATATYPE_FLT:
    attVDt = DATATYPE_FLT64;
    elemSize = sizeof(double);
    break;
  case DATATYPE_INT:
    attVDt = DATATYPE_INT;
    elemSize = sizeof(int);
    break;
  case DATATYPE_TXT:
    attVDt = DATATYPE_TXT;
    elemSize = 1;
    break;
  default:
    xabort("Unknown datatype encountered in attribute %s: %d\n",
           attName, tempbuf[2]);
  }
  void *attData = (void *)xmalloc(elemSize * (size_t)tempbuf[3]);
  serializeUnpack(buf, size, position, attData, tempbuf[3], attVDt, context);
  vlist_def_att(tempbuf[2], tempbuf[1], vlistID, varID, attName,
                (size_t)tempbuf[3], (size_t)tempbuf[3] * elemSize, attData);
  free(attName);
  free(attData);
}

void
vlistAttsUnpack(int vlistID, int varID,
                void * buf, int size, int *position, void *context)
{
  int numAtts, i;
  serializeUnpack(buf, size, position, &numAtts, 1, DATATYPE_INT, context);
  for (i = 0; i < numAtts; ++i)
  {
    vlistAttUnpack(vlistID, varID, buf, size, position, context);
  }
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
