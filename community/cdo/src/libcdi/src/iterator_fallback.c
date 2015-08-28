#include "iterator_fallback.h"

#include "cdi.h"
#include "cdi_int.h"
#include "dmemory.h"
#include "vlist.h"      //Required for vlist_t, which we require because there is no safe function available to access a variable name.

#include <assert.h>
#include <limits.h>
#include <stdlib.h>

//For more information on the condestruct() pattern, see comment in src/iterator_grib.c
static CdiFallbackIterator* cdiFallbackIterator_condestruct(CdiFallbackIterator* me, const char* path, int filetype)
{
  if(me) goto destruct;

  me = xmalloc(sizeof(*me));
  baseIterConstruct(&me->super, filetype);

  me->streamId = streamOpenRead(path);
  if(me->streamId == CDI_UNDEFID) goto destructSuper;
  me->vlistId = streamInqVlist(me->streamId);
  if(me->vlistId == CDI_UNDEFID) goto closeStream;
  me->variableCount = vlistNvars(me->vlistId);
  if(me->variableCount <= 0) goto closeStream;
  me->curLevelCount = -1;        //Will be set in cdiFallbackIterator_nextField()

  //These values are chosen so that the natural increment at the start of cdiFallbackIterator_nextField() will correctly position us at the first slice.
  me->curTimestep = 0;
  if(streamInqTimestep(me->streamId, me->curTimestep) <= 0) goto closeStream;
  me->curVariable = 0;
  me->curLevel = -1;
  me->path = strdup(path);
  if(!me->path) goto closeStream;

  return me;

// ^        constructor code        ^
// |                                |
// v destructor/error-cleanup code  v

destruct:
  free(me->path);
closeStream:
  streamClose(me->streamId);
destructSuper:
  baseIterDestruct(&me->super);
  free(me);
  return NULL;
}

CdiIterator* cdiFallbackIterator_new(const char* path, int filetype)
{
  return &cdiFallbackIterator_condestruct(NULL, path, filetype)->super;
}

//Fetches the info that is published by the variables in the base class from the current field.
static void fetchSuperInfo(CdiFallbackIterator* me)
{
  me->super.datatype = vlistInqVarDatatype(me->vlistId, me->curVariable);
  me->super.timesteptype = vlistInqVarTsteptype(me->vlistId, me->curVariable);
  me->super.gridId = vlistInqVarGrid(me->vlistId, me->curVariable);
  int param = vlistInqVarParam(me->vlistId, me->curVariable);
  cdiDecodeParam(param, &me->super.param.number, &me->super.param.category, &me->super.param.discipline);
}

CdiFallbackIterator* cdiFallbackIterator_clone(CdiIterator* super)
{
  CdiFallbackIterator* me = (CdiFallbackIterator*)super;

  //Make another stream for this file. This yields an unadvanced iterator.
  CdiFallbackIterator* clone = cdiFallbackIterator_condestruct(NULL, me->path, me->super.filetype);
  if(!clone) return NULL;

  //Point the clone to the same position in the file.
  clone->variableCount = me->variableCount;
  clone->curVariable = me->curVariable;
  clone->curLevelCount = me->curLevelCount;
  clone->curLevel = me->curLevel;
  clone->curTimestep = me->curTimestep;

  clone->super.isAdvanced = super->isAdvanced;
  if(super->isAdvanced) fetchSuperInfo(clone);

  return clone;
}

char* cdiFallbackIterator_serialize(CdiIterator* super)
{
  CdiFallbackIterator* me = (CdiFallbackIterator*)super;

  char* escapedPath = cdiEscapeSpaces(me->path);
  char* result = xmalloc(strlen(escapedPath)
                         + 5 * (3 * sizeof (int) * CHAR_BIT / 8 + 1) + 1);
  sprintf(result, "%s %d %d %d %d %d", escapedPath, me->variableCount, me->curVariable, me->curLevelCount, me->curLevel, me->curTimestep);
  free(escapedPath);
  return result;
}

CdiFallbackIterator* cdiFallbackIterator_deserialize(const char* description)
{
  CdiFallbackIterator* me = xmalloc(sizeof(*me));
  if(!me) goto fail;

  description = baseIter_constructFromString(&me->super, description);

  while(*description == ' ') description++;
  me->path = cdiUnescapeSpaces(description, &description);
  if(!me->path) goto destructSuper;

  me->streamId = streamOpenRead(me->path);
  if(me->streamId == CDI_UNDEFID) goto freePath;
  me->vlistId = streamInqVlist(me->streamId);
  if(me->vlistId == CDI_UNDEFID) goto closeStream;

  //This reads one variable from the description string, does error checking, and advances the given string pointer.
  #define decodeValue(variable, description) do \
    { \
      const char* savedStart = description; \
      long long decodedValue = strtoll(description, (char**)&description, 0);   /*The cast is a workaround for the wrong signature of strtoll().*/ \
      variable = (int)decodedValue; \
      if(savedStart == description) goto closeStream; \
      if((long long)decodedValue != (long long)variable) goto closeStream; \
    } while(0)
  decodeValue(me->variableCount, description);
  decodeValue(me->curVariable, description);
  decodeValue(me->curLevelCount, description);
  decodeValue(me->curLevel, description);
  decodeValue(me->curTimestep, description);
  #undef decodeValue

  if(streamInqTimestep(me->streamId, me->curTimestep) <= 0) goto closeStream;
  if(me->super.isAdvanced) fetchSuperInfo(me);

  return me;

closeStream:
  streamClose(me->streamId);
freePath:
  free(me->path);
destructSuper:
  baseIterDestruct(&me->super);
  free(me);
fail:
  return NULL;
}

static int advance(CdiFallbackIterator* me)
{
  me->curLevel++;
  if(me->curLevel == me->curLevelCount)
    {
      me->curLevel = 0;
      me->curVariable++;
      if(me->curVariable == me->variableCount)
        {
          me->curVariable = 0;
          me->curTimestep++;
          if(streamInqTimestep(me->streamId, me->curTimestep) <= 0) return CDI_EEOF;
        }
    }
  return CDI_NOERR;
}

int cdiFallbackIterator_nextField(CdiIterator* super)
{
  CdiFallbackIterator* me = (CdiFallbackIterator*)super;
  int result = advance(me);
  if(result) return result;

  if(!me->curLevel)
    { //Fetch the information that may have changed (we are processing a new variable/timestep if this point is reached).
      fetchSuperInfo(me);
      me->curLevelCount = zaxisInqSize(vlistInqVarZaxis(me->vlistId, me->curVariable));
    }
  return CDI_NOERR;
}

char* cdiFallbackIterator_inqTime(CdiIterator* super, bool getEndTime)
{
  CdiFallbackIterator* me = (CdiFallbackIterator*)super;
  if(getEndTime) return NULL;   //The stream interface does not export the start/end times of statistical fields, so we treat all data as point of time data, returning the validity time as the start time.
  int taxisId = vlistInqTaxis(me->vlistId);
  int date = taxisInqVdate(taxisId);
  int time = taxisInqVtime(taxisId);
  int year, month, day, hour, minute, second;
  cdiDecodeDate(date, &year, &month, &day);
  cdiDecodeTime(time, &hour, &minute, &second);
  char *result = xmalloc(   4+1 +2+1 +2+1 +2+1 +2+1 +2+4+1);
  sprintf(result, "%04d-%02d-%02dT%02d:%02d:%02d.000", year, month, day, hour, minute, second);
  return result;
}

int cdiFallbackIterator_levelType(CdiIterator* super, int levelSelector, char** outName, char** outLongName, char** outStdName, char** outUnit)
{
  CdiFallbackIterator* me = (CdiFallbackIterator*)super;
  int zaxisId = vlistInqVarZaxis(me->vlistId, me->curVariable);
  (void)levelSelector;
  #define copyString(outPointer, function) do \
    { \
      if(outPointer) \
        { \
          char tempBuffer[CDI_MAX_NAME]; \
          function(zaxisId, tempBuffer); \
          *outPointer = strdup(tempBuffer); \
        } \
    } \
  while(0)
  copyString(outName, zaxisInqName);    //FIXME: zaxisInqName is unsafe.
  copyString(outLongName, zaxisInqLongname);    //FIXME: zaxisInqLongname is unsafe.
  copyString(outStdName, zaxisInqStdname);    //FIXME: zaxisInqStdname is unsafe.
  copyString(outUnit, zaxisInqUnits);    //FIXME: zaxisInqUnits is unsafe.
  #undef copyString
  return zaxisInqLtype(zaxisId);
}

int cdiFallbackIterator_level(CdiIterator* super, int levelSelector, double* outValue1, double* outValue2)
{
  CdiFallbackIterator* me = (CdiFallbackIterator*)super;
  int zaxisId = vlistInqVarZaxis(me->vlistId, me->curVariable);

  //handle NULL pointers once and for all
  double trash;
  if(!outValue1) outValue1 = &trash;
  if(!outValue2) outValue2 = &trash;

  //get the level value
  if(levelSelector)
    {
      *outValue1 = (zaxisInqLbounds(zaxisId, NULL))
                 ? zaxisInqLbound(zaxisId, me->curLevel)
                 : zaxisInqLevel(zaxisId, me->curLevel);
    }
  else
    {
      *outValue1 = (zaxisInqUbounds(zaxisId, NULL))
                 ? zaxisInqUbound(zaxisId, me->curLevel)
                 : zaxisInqLevel(zaxisId, me->curLevel);
    }
  *outValue2 = 0.0;

  //if this is a hybrid zaxis, lookup the coordinates in the vertical coordinate table
  ssize_t intLevel = (ssize_t)(2**outValue1);
  if(0 <= intLevel && intLevel < zaxisInqVctSize(zaxisId) - 1)
    {
      const double* coordinateTable = zaxisInqVctPtr(zaxisId);
      *outValue1 = coordinateTable[intLevel];
      *outValue2 = coordinateTable[intLevel + 1];
    }
  return CDI_NOERR;
}

int cdiFallbackIterator_zaxisUuid(CdiIterator* super, int* outVgridNumber, int* outLevelCount, unsigned char (*outUuid)[16])
{
  CdiFallbackIterator* me = (CdiFallbackIterator*)super;
  int zaxisId = vlistInqVarZaxis(me->vlistId, me->curVariable);
  if(zaxisInqLtype(zaxisId) != ZAXIS_HYBRID) return CDI_EINVAL;
  if(outVgridNumber) *outVgridNumber = zaxisInqNumber(zaxisId);
  if(outLevelCount) *outLevelCount = zaxisInqNlevRef(zaxisId);
  if(outUuid) zaxisInqUUID(zaxisId, *outUuid);
  return CDI_NOERR;
}

char* cdiFallbackIterator_copyVariableName(CdiIterator* super)
{
  CdiFallbackIterator* me = (CdiFallbackIterator*)super;
  return vlistCopyVarName(me->vlistId, me->curVariable);
}

void cdiFallbackIterator_readField(CdiIterator* super, double* buffer, size_t* nmiss)
{
  CdiFallbackIterator* me = (CdiFallbackIterator*)super;
  int missingValues = 0;
  streamReadVarSlice(me->streamId, me->curVariable, me->curLevel, buffer, &missingValues);
  if(nmiss) *nmiss = (size_t)missingValues;
}

void cdiFallbackIterator_readFieldF(CdiIterator* super, float* buffer, size_t* nmiss)
{
  CdiFallbackIterator* me = (CdiFallbackIterator*)super;
  int missingValues = 0;
  streamReadVarSliceF(me->streamId, me->curVariable, me->curLevel, buffer, &missingValues);
  if(nmiss) *nmiss = (size_t)missingValues;
}

void cdiFallbackIterator_delete(CdiIterator* super)
{
  CdiFallbackIterator* me = (CdiFallbackIterator*)super;
  cdiFallbackIterator_condestruct(me, NULL, 0);
}
