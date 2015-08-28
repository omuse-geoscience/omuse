#include "iterator_grib.h"

#include "cdi_int.h"
#include "cgribex.h"
#include "dmemory.h"
#include "error.h"
#include "gribapi.h"
#include "gribapi_utilities.h"
#include "stream_grb.h"
#include "zaxis.h"

#include <assert.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>


#ifdef HAVE_LIBGRIB_API

//Since the error handling in constructors is usually very closely related to the workings of a destructor,
//this function combines both functions in one, using a centralized exit.
//The mode of operation depends on whether me is a NULL pointer on entry:
//If it is NULL, a new object is allocated and constructed, which is returned if construction is successful.
//If a non-NULL pointer is passed in, the object is destructed and NULL is returned. In this case, the other arguments are ignored.
static CdiGribIterator* cdiGribIterator_condestruct(CdiGribIterator* me, const char* path, int filetype)
{
  #define super() (&me->super)
  if(me) goto destruct;
  me = xmalloc(sizeof(*me));
  baseIterConstruct(super(), filetype);

  me->file = cdiInputFile_make(path);
  if(!me->file) goto destructSuper;
  me->fileOffset = 0;
  me->gribHandle = NULL;
  me->gribBuffer = NULL;
  me->bufferSize = me->curRecordSize = 0;
  me->super.gridId = CDI_UNDEFID;

  goto success;

// ^        constructor code        ^
// |                                |
// v destructor/error-cleanup code  v

destruct:
  if(me->super.gridId != CDI_UNDEFID) gridDestroy(me->super.gridId);
  if(me->gribHandle) grib_handle_delete(me->gribHandle);
  free(me->gribBuffer);
  cdiRefObject_release(&me->file->super);
destructSuper:
  baseIterDestruct(super());
  free(me);
  me = NULL;

success:
  return me;
  #undef super
}

CdiIterator* cdiGribIterator_new(const char* path, int filetype)
{
  return &cdiGribIterator_condestruct(NULL, path, filetype)->super;
}

CdiGribIterator* cdiGribIterator_makeClone(CdiIterator* super)
{
  CdiGribIterator* me = (CdiGribIterator*)super;

  //Allocate memory and copy data. (operations that may fail)
  CdiGribIterator* result = xmalloc(sizeof(*result));
  if(!result) goto fail;
  *result = (CdiGribIterator)
    {
      .file = me->file,
      .fileOffset = me->fileOffset,
      .gribBuffer = NULL,
      .bufferSize = me->bufferSize,
      .curRecordSize = me->curRecordSize,
      .gribHandle = NULL
    };
  if(me->gribBuffer)
    {
      result->gribBuffer = xmalloc(me->bufferSize);
      if(!result->gribBuffer) goto freeResult;
      memcpy(result->gribBuffer, me->gribBuffer, me->curRecordSize);
    }
  if(me->gribHandle)
    {
      result->gribHandle = grib_handle_new_from_message(NULL, result->gribBuffer, result->curRecordSize);
      if(!result->gribHandle) goto freeBuffer;
    }
  if(super->gridId != CDI_UNDEFID)
    {
      result->super.gridId = gridDuplicate(super->gridId);
      if(result->super.gridId == CDI_UNDEFID) goto deleteGribHandle;
    }

  //Finish construction. (operations that cannot fail)
  baseIterConstruct(&result->super, super->filetype);
  result->super.datatype = super->datatype;
  result->super.timesteptype = super->timesteptype;
  result->super.param = super->param;
  cdiRefObject_retain(&result->file->super);

  return result;

  //Error handling.
deleteGribHandle:
  if(result->gribHandle) grib_handle_delete(result->gribHandle);
freeBuffer:
  free(result->gribBuffer);
freeResult:
  free(result);
fail:
  return NULL;
}

char* cdiGribIterator_serialize(CdiIterator* super)
{
  CdiGribIterator* me = (CdiGribIterator*)super;

  const char* path = cdiInputFile_getPath(me->file);
  char* escapedPath = cdiEscapeSpaces(path);
  char* result = xmalloc(strlen(escapedPath) + 3 * sizeof(int) * CHAR_BIT/8);
  sprintf(result, "%s %zu", escapedPath, me->fileOffset);
  free(escapedPath);
  return result;
}


CdiGribIterator* cdiGribIterator_deserialize(const char* description)
{
  char* path;
  CdiGribIterator* me = xmalloc(sizeof(*me));
  if(!me) goto fail;

  description = baseIter_constructFromString(&me->super, description);

  while(*description == ' ') description++;
  path = cdiUnescapeSpaces(description, &description);
  if(!path) goto destructSuper;

  me->file = cdiInputFile_make(path);
  free(path);
  if(!me->file) goto destructSuper;

  {
    const char* savedStart = description;
    long long decodedOffset = strtoll(description, (char**)&description, 0);    //The cast is a workaround for the wrong signature of strtoll() (it should have been `long long strtoll(const char*, const char**, int)`, not `long long strtoll(const char*, char**, int)`.
    me->fileOffset = (off_t)decodedOffset;
    if(savedStart == description) goto closeFile;
    if((unsigned long long)decodedOffset > (unsigned long long)me->fileOffset) goto closeFile;
  }

  me->gribBuffer = NULL;
  me->bufferSize = me->curRecordSize = 0;
  me->gribHandle = NULL;
  me->super.gridId = CDI_UNDEFID;
  if(me->super.isAdvanced && cdiGribIterator_nextField(&me->super)) goto closeFile;

  return me;


closeFile:
  cdiRefObject_release(&me->file->super);
destructSuper:
  baseIterDestruct(&me->super);
  free(me);
fail:
  return NULL;
}

static void cdiGribIterator_ensureBuffer(CdiGribIterator* me, size_t requiredSize)
{
  if(me->bufferSize < requiredSize)
    {
      me->bufferSize *= 2;
      if(me->bufferSize < requiredSize) me->bufferSize = requiredSize;
      me->gribBuffer = xrealloc(me->gribBuffer, me->bufferSize);
    }
}

static bool isGrib1DualLevel(int levelType)
{
  switch(levelType)
    {
      case 101: case 104: case 106: case 108: case 110: case 112:
      case 114: case 116: case 120: case 121: case 128: case 141:   //This is the complete list after grib_api-1.12.3/definitions/grib1/sections.1.def:106-117:, the code in cdi/src/stream_gribapi.c:grib1GetLevel() seems to be incomplete.
        return true;
      default:
        return false;
    }
}

static const unsigned char* positionOfGribMarker(const unsigned char* data, size_t size)
{
  for(const unsigned char* currentPosition = data, *end = data + size; currentPosition < end; currentPosition++)
    {
      currentPosition = memchr(currentPosition, 'G', size - (size_t)(currentPosition - data) - 3);      //-3 to ensure that we don't overrun the buffer during the strncmp() call.
      if(!currentPosition) return NULL;
      if(!strncmp((const char*)currentPosition, "GRIB", 4)) return currentPosition;
    }
  return NULL;
}

//This clobbers the contents of the gribBuffer!
//Returns the file offset of the next 'GRIB' marker.
static ssize_t scanToGribMarker(CdiGribIterator* me)
{
  cdiGribIterator_ensureBuffer(me, 8*1024);
  const size_t kMaxScanSize = 16*1024*1024;
  for(size_t scannedBytes = 0, scanSize; scannedBytes < kMaxScanSize; scannedBytes += scanSize)
    {
      //Load a chunk of data into our buffer.
      scanSize = me->bufferSize;
      if(scannedBytes + scanSize > kMaxScanSize) scanSize = kMaxScanSize - scannedBytes;
      assert(scanSize <= me->bufferSize);
      int status = cdiInputFile_read(me->file, me->fileOffset + (off_t)scannedBytes, scanSize, &scanSize, me->gribBuffer);
      if(status != CDI_NOERR && status != CDI_EEOF) return status;

      const unsigned char* startPosition = positionOfGribMarker(me->gribBuffer, scanSize);
      if(startPosition)
        {
          return (ssize_t)(me->fileOffset + (off_t)scannedBytes + (off_t)(startPosition - me->gribBuffer));
        }

      //Get the offset for the next iteration if there is a next iteration.
      scanSize -= 3;        //so that we won't miss a 'GRIB' sequence that happens to be cut off
      scannedBytes += scanSize;
      scannedBytes &= ~(size_t)0xf; //make 16 bytes aligned
    }
  return -1;
}

static unsigned decode24(void* beData)
{
  unsigned char* bytes = beData;
  return ((unsigned)bytes[0] << 16) + ((unsigned)bytes[1] << 8) + (unsigned)bytes[2];
}

static uint64_t decode64(void* beData)
{
  unsigned char* bytes = beData;
  uint64_t result = 0;
  for(size_t i = 0; i < 8; i++) result = (result << 8) + bytes[i];
  return result;
}

//Determine the size of the GRIB record that begins at the given file offset.
static int getRecordSize(CdiGribIterator* me, off_t gribFileOffset, size_t* outRecordSize)
{
  char buffer[16];
  size_t readSize;
  int status = cdiInputFile_read(me->file, gribFileOffset, sizeof(buffer), &readSize, buffer);
  if(status != CDI_NOERR && status != CDI_EEOF) return status;
  if(readSize < sizeof(buffer)) return CDI_EEOF;
  *outRecordSize = 0;
  switch(buffer[7])
    {
      case 1:
        *outRecordSize = decode24(&buffer[4]);
        if(*outRecordSize & (1 << 23))
          {
            *outRecordSize = 120*(*outRecordSize & ((1 << 23) - 1));    //Rescaling for long records.
            //The corresponding code in cgribexlib.c:4532-4570: is much more complicated
            //due to the fact that it subtracts the padding bytes that are inserted after section 4.
            //However, we are only interested in the total size of data we need to read here,
            //so we can ignore the presence of some padding bytes.
          }
        return CDI_NOERR;

      case 2:
        *outRecordSize =  decode64(&buffer[8]);
        return CDI_NOERR;

      default:
        return CDI_EUFTYPE;
    }
}

#if 0
static void hexdump(void* data, size_t size)
{
  unsigned char* charData = data;
  for(size_t offset = 0; offset < size; )
    {
      printf("%016zx:", offset);
      for(size_t i = 0; i < 64 && offset < size; i++, offset++)
        {
          if((i & 63) && !(i & 15)) printf(" |");
          if((i & 15) && !(i & 3)) printf("  ");
          printf(" %02x", charData[offset]);
        }
      printf("\n");
    }
}
#endif

//Read a record into memory and wrap it in a grib_handle.
//XXX: I have omitted checking for szip compression as it is done in grbReadVarDP() & friends since that appears to be a non-standard extension of the GRIB1 standard: bit 1 in octet 14 of the binary data section which is used to signal szip compressio is defined to be reserved in the standard. As such, it seems prudent not to support this and to encourage people with such szip compressed files to switch to the GRIB2/JPEG2000 format. However, in the case that this reasoning is wrong, this function is probably the place to add the check for zsip compression.
static int readMessage(CdiGribIterator* me)
{
  //Destroy the old grib_handle.
  if(me->gribHandle) grib_handle_delete(me->gribHandle), me->gribHandle = NULL;
  me->fileOffset += (off_t)me->curRecordSize;

  //Find the next record and determine its size.
  ssize_t gribFileOffset = scanToGribMarker(me);
  int result = CDI_EEOF;
  if(gribFileOffset < 0) goto fail;
  result = getRecordSize(me, gribFileOffset, &me->curRecordSize);
  if(result) goto fail;

  //Load the whole record into our buffer and create a grib_handle for it.
  cdiGribIterator_ensureBuffer(me, me->curRecordSize);
  result = cdiInputFile_read(me->file, gribFileOffset, me->curRecordSize, NULL, me->gribBuffer);
  if(result) goto fail;
  me->gribHandle = grib_handle_new_from_message(NULL, me->gribBuffer, me->curRecordSize);
  result = CDI_EUFSTRUCT;
  if(!me->gribHandle) goto fail;

  return CDI_NOERR;

fail:
  me->curRecordSize = 0;        //This ensures that we won't jump to an uncontrolled file position if cdiGribIterator_nextField() is called another time after it has returned an error.
  return result;
}

int cdiGribIterator_nextField(CdiIterator* super)
{
  CdiGribIterator* me = (CdiGribIterator*)super;

  if(super->gridId != CDI_UNDEFID) gridDestroy(super->gridId), super->gridId = CDI_UNDEFID;

  //Get the next GRIB message into our buffer.
  int result = readMessage(me);
  if(result) return result;

  //Get the metadata that's published as variables in the base class.
  super->datatype = gribGetDatatype(me->gribHandle);
  super->timesteptype = gribapiGetTsteptype(me->gribHandle);
  cdiDecodeParam(gribapiGetParam(me->gribHandle), &super->param.number, &super->param.category, &super->param.discipline);
  grid_t grid;
  gribapiGetGrid(me->gribHandle, &grid);
  super->gridId = gridGenerate(&grid);

  return CDI_NOERR;
}

char* cdiGribIterator_inqTime(CdiIterator* super, bool getEndTime)
{
  CdiGribIterator* me = (CdiGribIterator*)super;
  return gribMakeTimeString(me->gribHandle, getEndTime);
}

int cdiGribIterator_levelType(CdiIterator* super, int levelSelector, char** outName, char** outLongName, char** outStdName, char** outUnit)
{
  CdiGribIterator* me = (CdiGribIterator*)super;

  //First determine the zaxis type corresponding to the given level.
  int zaxisType = ZAXIS_GENERIC;
  if(gribEditionNumber(me->gribHandle) <= 1)
    {
      int levelType = (int)gribGetLongDefault(me->gribHandle, "indicatorOfTypeOfLevel", 255);
      if(levelSelector && !isGrib1DualLevel(levelType)) levelType = 255;
      zaxisType = grib1ltypeToZaxisType(levelType);
    }
  else
    {
      int levelType = (int)gribGetLongDefault(me->gribHandle, levelSelector ? "typeOfSecondFixedSurface" : "typeOfFirstFixedSurface", 255);
      zaxisType = grib2ltypeToZaxisType(levelType);
    }

  //Then lookup the requested names.
  const char* name, *longName, *stdName, *unit;
  zaxisGetTypeDescription(zaxisType, NULL, &name, &longName, &stdName, &unit);
  if(outName) *outName = strdup(name);
  if(outLongName) *outLongName = strdup(longName);
  if(outStdName) *outStdName = strdup(stdName);
  if(outUnit) *outUnit = strdup(unit);

  return zaxisType;
}

static double logicalLevelValue2(long gribType, long storedValue, long power)
{
  double factor = 1;
  while(power--) factor *= 10;      //this is precise up to factor == 22.
  switch(gribType)
    {
      case GRIB2_LTYPE_LANDDEPTH:
      case GRIB2_LTYPE_ISOBARIC:
      case GRIB2_LTYPE_SIGMA:
        return (double)storedValue * (1000.0/factor);      //The evaluation order allows the factors of ten to cancel out before rounding.

      case 255:
        return 0;

      default:
        return (double)storedValue/factor;
    }
}

//The output values must be preinitialized, this function does not always write them.
static int readLevel2(grib_handle* gribHandle, const char* levelTypeKey, const char* powerKey, const char* valueKey, double* outValue1, double* outValue2)
{
  assert(levelTypeKey && powerKey && valueKey && outValue1 && outValue2);

  long levelType = gribGetLongDefault(gribHandle, levelTypeKey, 255);   //1 byte
  switch(levelType)
    {
      case 255: break;

      case 105: case 113:
        {
          unsigned long value = (unsigned long)gribGetLongDefault(gribHandle, valueKey, 0);
          unsigned long coordinateCount = (unsigned long)gribGetLongDefault(gribHandle, "numberOfCoordinatesValues", 0);
          if(value >= coordinateCount/2)
            {
              Error("Invalid level coordinate: Level has the hybrid coordinate index %lu, but only %lu coordinate pairs are present.", value, coordinateCount/2);
              return CDI_EUFSTRUCT;
            }
          int status;
          //XXX: I'm not 100% sure about how the coordinate pairs are stored in the file.
          //     I'm assuming an array of pairs due to the example code in grib_api-1.12.3/examples/F90/set_pv.f90, but that may be wrong.
          if((status = (int)grib_get_double_element(gribHandle, "pv", (int)value*2    , outValue1))) return status;
          if((status = (int)grib_get_double_element(gribHandle, "pv", (int)value*2 + 1, outValue2))) return status;
          break;
        }

      default:
        {
          long power = gribGetLongDefault(gribHandle, powerKey, 0);  //1 byte
          if(power == 255) power = 0;
          long value = gribGetLongDefault(gribHandle, valueKey, 0);   //4 bytes
          *outValue1 = logicalLevelValue2(levelType, value, power);
        }
    }
  return CDI_NOERR;
}

int cdiGribIterator_level(CdiIterator* super, int levelSelector, double* outValue1, double* outValue2)
{
  CdiGribIterator* me = (CdiGribIterator*)super;
  double trash;
  if(!outValue1) outValue1 = &trash;
  if(!outValue2) outValue2 = &trash;
  *outValue1 = *outValue2 = 0;

  if(gribEditionNumber(me->gribHandle) > 1)
    {
      if(levelSelector)
        {
          return readLevel2(me->gribHandle, "typeOfFirstFixedSurface", "scaleFactorOfFirstFixedSurface", "scaledValueOfFirstFixedSurface", outValue1, outValue2);
        }
      else
        {
          return readLevel2(me->gribHandle, "typeOfSecondFixedSurface", "scaleFactorOfSecondFixedSurface", "scaledValueOfSecondFixedSurface", outValue1, outValue2);
        }
    }
  else
    {
      long levelType = (uint8_t)gribGetLongDefault(me->gribHandle, "indicatorOfTypeOfLevel", -1);    //1 byte
      if(levelType == 255)
        {}
      else if(isGrib1DualLevel((int)levelType))
        {
          *outValue1 = (double)gribGetLongDefault(me->gribHandle, (levelSelector ? "bottomLevel" : "topLevel"), 0);
        }
      else if(levelType == 100)
        {
          *outValue1 = 100 * (double)gribGetLongDefault(me->gribHandle, "level", 0);        //2 bytes
        }
      else
        {
          *outValue1 = (double)gribGetLongDefault(me->gribHandle, "level", 0);        //2 bytes
        }
    }
  return CDI_NOERR;
}

int cdiGribIterator_zaxisUuid(CdiIterator* super, int* outVgridNumber, int* outLevelCount, unsigned char (*outUuid)[16])
{
  CdiGribIterator* me = (CdiGribIterator*)super;

  if(outVgridNumber)
    {
      long temp;
      if(grib_get_long(me->gribHandle, "numberOfVGridUsed", &temp)) return CDI_EINVAL;
      *outVgridNumber = (int)temp;
    }
  if(outLevelCount)
    {
      long temp;
      if(grib_get_long(me->gribHandle, "nlev", &temp)) return CDI_EINVAL;
      *outLevelCount = (int)temp;
    }
  if(outUuid)
    {
      size_t size = sizeof(*outUuid);
      if(grib_get_bytes(me->gribHandle, "uuidOfVGrid", *outUuid, &size)) return CDI_EINVAL;
      if(size != sizeof(*outUuid)) return CDI_EUFSTRUCT;
    }

  return CDI_NOERR;
}

char* cdiGribIterator_copyVariableName(CdiIterator* super)
{
  CdiGribIterator* me = (CdiGribIterator*)super;
  return gribCopyString(me->gribHandle, "shortName");
}

void cdiGribIterator_readField(CdiIterator* super, double* buffer, size_t* nmiss)
{
  CdiGribIterator* me = (CdiGribIterator*)super;

  gribGetDoubleArray(me->gribHandle, "values", buffer);
  long gridType = gribGetLong(me->gribHandle, "gridDefinitionTemplateNumber");
  if(nmiss)
    {
      *nmiss = (gridType >= 50 && gridType <= 53) ? (size_t)0 : (size_t)gribGetLong(me->gribHandle, "numberOfMissing");        //The condition excludes harmonic data.
    }
}

void cdiGribIterator_readFieldF(CdiIterator* super, float* buffer, size_t* nmiss)
{
  CdiGribIterator* me = (CdiGribIterator*)super;

  size_t valueCount = gribGetArraySize(me->gribHandle, "values");
  double* temp = malloc(valueCount*sizeof(*temp));
  cdiGribIterator_readField(super, temp, nmiss);
  for(size_t i = valueCount; i--; ) buffer[i] = (float)temp[i];
  free(temp);
}
#endif

/**
@Function cdiGribIterator_delete
@Title Dispose off a CdiGribIterator instance.

@Prototype void cdiGribIterator_delete(CdiGribIterator* me)
@Parameter
    @item me The iterator to delete.

@Description
    Combined destructor and deallocator. Make sure to match every call to cdiGribIterator_clone() with a call to this function.
*/
void cdiGribIterator_delete(CdiGribIterator* me)
{
#ifdef HAVE_LIBGRIB_API
  if(me) cdiGribIterator_condestruct(me, NULL, 0);
#else
  (void)me;
  xabort("CDI was compiled without GribAPI support, so you can't possibly have a valid CdiGribIterator* to call this function with");
#endif
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// callthroughs to provide direct access to the grib keys //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

/**
@Function cdiGribIterator_inqEdition
@Title Get the version of the GRIB standard that is used

@Prototype int cdiGribIterator_inqEdition(CdiGribIterator* me)
@Parameter
    @item me The iterator to operate on.

@Result The GRIB version.

@Description
    Returns the version of the file format.
*/
int cdiGribIterator_inqEdition(CdiGribIterator* me)
{
#ifdef HAVE_LIBGRIB_API
  return (int)gribEditionNumber(me->gribHandle);
#else
  (void)me;
  xabort("CDI was compiled without GribAPI support, so you can't possibly have a valid CdiGribIterator* to call this function with");
#endif
}

/**
@Function cdiGribIterator_getLong
@Title Access to grib_get_long()

@Prototype int cdiGribIterator_getLong(CdiGribIterator* me, const char* key, long* result)
@Parameter
    @item me The iterator to operate on.
    @item ... The arguments to the underlying GRIB-API function.

@Result An error code.

@Description
    Callthrough to grib_get_long().
*/
int cdiGribIterator_getLong(CdiGribIterator* me, const char* key, long* result)
{
#ifdef HAVE_LIBGRIB_API
  return grib_get_long(me->gribHandle, key, result);
#else
  (void)me;
  (void)key;
  (void)result;
  xabort("CDI was compiled without GribAPI support, so you can't possibly have a valid CdiGribIterator* to call this function with");
#endif
}

/**
@Function cdiGribIterator_getLength
@Title Access to grib_get_length()

@Prototype int cdiGribIterator_getLength(CdiGribIterator* me, const char* key, size_t* result)
@Parameter
    @item me The iterator to operate on.
    @item ... The arguments to the underlying GRIB-API function.

@Result An error code.

@Description
    Callthrough to grib_get_length().
*/
int cdiGribIterator_getLength(CdiGribIterator* me, const char* key, size_t* result)
{
#ifdef HAVE_GRIB_GET_LENGTH
  return grib_get_length(me->gribHandle, key, result);
#elif defined(HAVE_LIBGRIB_API)
  (void)me;
  (void)key;
  (void)result;
  Error("grib_get_length() is not available, so cdiGribIterator_getLength() can't be used");
  return -1;
#else
  (void)me;
  (void)key;
  (void)result;
  xabort("CDI was compiled without GribAPI support, so you can't possibly have a valid CdiGribIterator* to call this function with");
#endif
}

/**
@Function cdiGribIterator_getString
@Title Access to grib_get_string()

@Prototype int cdiGribIterator_getString(CdiGribIterator* me, const char* key, char* result, size_t* length)
@Parameter
    @item me The iterator to operate on.
    @item ... The arguments to the underlying GRIB-API function.

@Result An error code.

@Description
    Callthrough to grib_get_string().
*/
int cdiGribIterator_getString(CdiGribIterator* me, const char* key, char* result, size_t* length)
{
#ifdef HAVE_LIBGRIB_API
  return grib_get_string(me->gribHandle, key, result, length);
#else
  (void)me;
  (void)key;
  (void)result;
  (void)length;
  xabort("CDI was compiled without GribAPI support, so you can't possibly have a valid CdiGribIterator* to call this function with");
#endif
}

/**
@Function cdiGribIterator_inqLongValue
@Title Get the value of a GRIB-API key as a long

@Prototype long cdiGribIterator_inqLongValue(CdiGribIterator* me, const char* key)
@Parameter
    @item me The iterator to operate on.
    @item key The GRIB-API key to retrieve.

@Result The value of the key.

@Description
    Use this to fetch a grib value if you are certain that the given key must be present.
    This will abort the process if the key cannot be retrieved.
*/
long cdiGribIterator_inqLongValue(CdiGribIterator* me, const char* key)
{
#ifdef HAVE_LIBGRIB_API
  return gribGetLong(me->gribHandle, key);
#else
  (void)me;
  (void)key;
  xabort("CDI was compiled without GribAPI support, so you can't possibly have a valid CdiGribIterator* to call this function with");
#endif
}

/**
@Function cdiGribIterator_inqLongDefaultValue
@Title Get the value of a GRIB-API key as a long

@Prototype long cdiGribIterator_inqLongDefaultValue(CdiGribIterator* me, const char* key, long defaultValue)
@Parameter
    @item me The iterator to operate on.
    @item key The GRIB-API key to retrieve.
    @item defaultValue The value to return if the key is not present.

@Result The value of the key or the given default value.

@Description
    Use this if you can handle failure to fetch the key by supplying a default value.
    This function cannot fail.
*/
long cdiGribIterator_inqLongDefaultValue(CdiGribIterator* me, const char* key, long defaultValue)
{
#ifdef HAVE_LIBGRIB_API
  return gribGetLongDefault(me->gribHandle, key, defaultValue);
#else
  (void)me;
  (void)key;
  (void)defaultValue;
  xabort("CDI was compiled without GribAPI support, so you can't possibly have a valid CdiGribIterator* to call this function with");
#endif
}

/**
@Function cdiGribIterator_inqStringValue
@Title Safely retrieve a GRIB-API key with a string value

@Prototype char* cdiGribIterator_inqStringValue(CdiGribIterator* me, const char* key)
@Parameter
    @item me The iterator to operate on.
    @item key The GRIB-API key to retrieve.

@Result A malloc'ed string or NULL.

@Description
    This will first call grib_get_length() to inquire the actual size of the string,
    allocate memory accordingly, call grib_get_string(), and return the pointer to the new string.
    Returns NULL on failure.
*/
char* cdiGribIterator_inqStringValue(CdiGribIterator* me, const char* key)
{
#ifdef HAVE_LIBGRIB_API
  return gribCopyString(me->gribHandle, key);
#else
  (void)me;
  (void)key;
  xabort("CDI was compiled without GribAPI support, so you can't possibly have a valid CdiGribIterator* to call this function with");
#endif
}

/**
@Function cdiGribIterator_getDouble
@Title Access to grib_get_double()

@Prototype int cdiGribIterator_getDouble(CdiGribIterator* me, const char* key, double* result)
@Parameter
    @item me The iterator to operate on.
    @item ... The arguments to the underlying GRIB-API function.

@Result An error code.

@Description
    Callthrough to grib_get_double().
*/
int cdiGribIterator_getDouble(CdiGribIterator* me, const char* key, double* result)
{
#ifdef HAVE_LIBGRIB_API
  return grib_get_double(me->gribHandle, key, result);
#else
  (void)me;
  (void)key;
  (void)result;
  xabort("CDI was compiled without GribAPI support, so you can't possibly have a valid CdiGribIterator* to call this function with");
#endif
}

/**
@Function cdiGribIterator_getSize
@Title Access to grib_get_size()

@Prototype int cdiGribIterator_getSize(CdiGribIterator* me, const char* key, size_t* result)
@Parameter
    @item me The iterator to operate on.
    @item ... The arguments to the underlying GRIB-API function.

@Result An error code.

@Description
    Callthrough to grib_get_size().
*/
int cdiGribIterator_getSize(CdiGribIterator* me, const char* key, size_t* result)
{
#ifdef HAVE_LIBGRIB_API
  return grib_get_size(me->gribHandle, key, result);
#else
  (void)me;
  (void)key;
  (void)result;
  xabort("CDI was compiled without GribAPI support, so you can't possibly have a valid CdiGribIterator* to call this function with");
#endif
}

/**
@Function cdiGribIterator_getLongArray
@Title Access to grib_get_long_array()

@Prototype int cdiGribIterator_getLongArray(CdiGribIterator* me, const char* key, long* result, size_t* size)
@Parameter
    @item me The iterator to operate on.
    @item ... The arguments to the underlying GRIB-API function.

@Result An error code.

@Description
    Callthrough to grib_get_long_array().
*/
int cdiGribIterator_getLongArray(CdiGribIterator* me, const char* key, long* result, size_t* size)
{
#ifdef HAVE_LIBGRIB_API
  return grib_get_long_array(me->gribHandle, key, result, size);
#else
  (void)me;
  (void)key;
  (void)result;
  (void)size;
  xabort("CDI was compiled without GribAPI support, so you can't possibly have a valid CdiGribIterator* to call this function with");
#endif
}

/**
@Function cdiGribIterator_getDoubleArray
@Title Access to grib_get_double_array()

@Prototype int cdiGribIterator_getDoubleArray(CdiGribIterator* me, const char* key, double* result, size_t* size)
@Parameter
    @item me The iterator to operate on.
    @item ... The arguments to the underlying GRIB-API function.

@Result An error code.

@Description
    Callthrough to grib_get_double_array().
*/
int cdiGribIterator_getDoubleArray(CdiGribIterator* me, const char* key, double* result, size_t* size)
{
#ifdef HAVE_LIBGRIB_API
  return grib_get_double_array(me->gribHandle, key, result, size);
#else
  (void)me;
  (void)key;
  (void)result;
  (void)size;
  xabort("CDI was compiled without GribAPI support, so you can't possibly have a valid CdiGribIterator* to call this function with");
#endif
}

/**
@Function cdiGribIterator_inqDoubleValue
@Title Get the value of a GRIB-API key as a double

@Prototype double cdiGribIterator_inqDoubleValue(CdiGribIterator* me, const char* key)
@Parameter
    @item me The iterator to operate on.
    @item key The GRIB-API key to retrieve.

@Result The value of the key.

@Description
    Use this to fetch a grib value if you are certain that the given key must be present.
    This will abort the process if the key cannot be retrieved.
*/
double cdiGribIterator_inqDoubleValue(CdiGribIterator* me, const char* key)
{
#ifdef HAVE_LIBGRIB_API
  return gribGetDouble(me->gribHandle, key);
#else
  (void)me;
  (void)key;
  xabort("CDI was compiled without GribAPI support, so you can't possibly have a valid CdiGribIterator* to call this function with");
#endif
}

/**
@Function cdiGribIterator_inqDoubleDefaultValue
@Title Get the value of a GRIB-API key as a double

@Prototype double cdiGribIterator_inqDoubleDefaultValue(CdiGribIterator* me, const char* key, double defaultValue)
@Parameter
    @item me The iterator to operate on.
    @item key The GRIB-API key to retrieve.
    @item defaultValue The value to return if the key is not present.

@Result The value of the key or the given default value.

@Description
    Use this if you can handle failure to fetch the key by supplying a default value.
    This function cannot fail.
*/
double cdiGribIterator_inqDoubleDefaultValue(CdiGribIterator* me, const char* key, double defaultValue)
{
#ifdef HAVE_LIBGRIB_API
  return gribGetDoubleDefault(me->gribHandle, key, defaultValue);
#else
  (void)me;
  (void)key;
  (void)defaultValue;
  xabort("CDI was compiled without GribAPI support, so you can't possibly have a valid CdiGribIterator* to call this function with");
#endif
}
