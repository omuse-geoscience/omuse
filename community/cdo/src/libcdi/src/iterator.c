#include "cdi.h"
#include "dmemory.h"
#include "iterator.h"
#include "iterator_fallback.h"
#include "iterator_grib.h"
#include "cdi_int.h"

#include <assert.h>
#include <ctype.h>

const char* const kUnexpectedFileTypeMessage = "Internal error: Unexpected file type encountered in iterator.\nThis is either due to an illegal memory access by the application or an internal logical error in CDI (unlikely, but possible).";
const char* const kAdvancedString = "advanced";
const char* const kUnadvancedString = "unadvanced";

//Returns a static string.
static const char* fileType2String(int fileType)
{
  switch(fileType)
    {
      #ifdef HAVE_LIBGRIB_API
        case FILETYPE_GRB: return "CDI::Iterator::GRIB1";
        case FILETYPE_GRB2: return "CDI::Iterator::GRIB2";
      #endif
      #ifdef HAVE_LIBNETCDF
        case FILETYPE_NC: return "CDI::Iterator::NetCDF";
        case FILETYPE_NC2: return "CDI::Iterator::NetCDF2";
        case FILETYPE_NC4: return "CDI::Iterator::NetCDF4";
        case FILETYPE_NC4C: return "CDI::Iterator::NetCDF4C";
      #endif
      #ifdef HAVE_LIBSERVICE
        case FILETYPE_SRV: return "CDI::Iterator::SRV";
      #endif
      #ifdef HAVE_LIBEXTRA
        case FILETYPE_EXT: return "CDI::Iterator::EXT";
      #endif
      #ifdef HAVE_LIBIEG
        case FILETYPE_IEG: return "CDI::Iterator::IEG";
      #endif

      default: return NULL;
    }
}

static int string2FileType(const char* fileType, const char** outRestString)
{
  //This first part unconditionally checks all known type strings, and only if the given string matches one of these strings, we use fileType2string() to check whether support for this type has been compiled in. This is to avoid throwing "invalid type string" errors when we just have a library version mismatch.
  #define check(givenString, typeString, typeConstant) do \
    { \
      if(givenString == strstr(givenString, typeString)) \
        { \
          if(outRestString) *outRestString = givenString + strlen(typeString); \
          if(fileType2String(typeConstant)) return typeConstant; \
          Error("Support for " typeString " not compiled in. Please check that the result of `cdiIterator_serialize()` is only passed to a `cdiIterator_deserialize()` implementation of the same CDI library version."); \
          return FILETYPE_UNDEF; \
        } \
    } while(0)
  check(fileType, "CDI::Iterator::GRIB1", FILETYPE_GRB);
  check(fileType, "CDI::Iterator::GRIB2", FILETYPE_GRB2);
  check(fileType, "CDI::Iterator::NetCDF", FILETYPE_NC);
  check(fileType, "CDI::Iterator::NetCDF2", FILETYPE_NC2);
  check(fileType, "CDI::Iterator::NetCDF4", FILETYPE_NC4);
  check(fileType, "CDI::Iterator::NetCDF4C", FILETYPE_NC4C);
  check(fileType, "CDI::Iterator::SRV", FILETYPE_SRV);
  check(fileType, "CDI::Iterator::EXT", FILETYPE_EXT);
  check(fileType, "CDI::Iterator::IEG", FILETYPE_IEG);
  #undef check

  //If this point is reached, the given string does not seem to be produced by a cdiIterator_serialize() call.
  Error("The string \"%s\" does not start with a valid iterator type. Please check the source of this string.", fileType);
  *outRestString = fileType;
  return FILETYPE_UNDEF;
}

/**
@Function cdiIterator_new
@Title Create an iterator for an input file

@Prototype CdiIterator* cdiIterator_new(const char* path)
@Parameter
    @item path Path to the file that is to be read.

@Result An iterator for the given file.

@Description
    Combined allocator and constructor for CdiIterator.

    The returned iterator does not point to the first field yet,
    it must first be advanced once before the first field can be introspected.
    This design decision has two benefits: 1. Empty files require no special
    cases, 2. Users can start a while(!cdiIterator_nextField(iterator)) loop
    right after the call to cdiIterator_new().
*/
CdiIterator* cdiIterator_new(const char* path)
{
  int trash;
  int filetype = cdiGetFiletype(path, &trash);
  switch(filetype)
    {
      case FILETYPE_UNDEF:
        Warning("Can't open file \"%s\": unknown format\n", path);
        return NULL;

      #ifdef HAVE_LIBGRIB_API
        case FILETYPE_GRB:
        case FILETYPE_GRB2:
          return cdiGribIterator_new(path, filetype);
      #endif

      #ifdef HAVE_LIBNETCDF
        case FILETYPE_NC:
        case FILETYPE_NC2:
        case FILETYPE_NC4:
        case FILETYPE_NC4C:
      #endif
      #ifdef HAVE_LIBSERVICE
        case FILETYPE_SRV:
      #endif
      #ifdef HAVE_LIBEXTRA
        case FILETYPE_EXT:
      #endif
      #ifdef HAVE_LIBIEG
        case FILETYPE_IEG:
      #endif
          return cdiFallbackIterator_new(path, filetype);

      default:
        Warning("the file \"%s\" is of type %s, but support for this format is not compiled in!", path, strfiletype(filetype));
        return NULL;
    }
}

void baseIterConstruct(CdiIterator* me, int filetype)
{
  me->filetype = filetype;
  me->isAdvanced = false;
}

const char* baseIter_constructFromString(CdiIterator* me, const char* description)
{
  const char* result = description;
  me->filetype = string2FileType(result, &result);
  assert(me->filetype != FILETYPE_UNDEF && "Please report this error.");        //This condition should have been checked for in a calling function.
  for(; *result && isspace(*result); result++);
  if(result == strstr(result, kAdvancedString))
    {
      me->isAdvanced = true;
      result += strlen(kAdvancedString);
    }
  else if(result == strstr(result, kUnadvancedString))
    {
      me->isAdvanced = false;
      result += strlen(kUnadvancedString);
    }
  else
    {
      Error("Invalid iterator description string \"%s\". Please check the origin of this string.", description);
      return NULL;
    }
  return result;
}

#define sanityCheck(me) do { \
    if(!me) xabort("NULL was passed to %s as an iterator. Please check the return value of cdiIterator_new().", __func__); \
    if(!me->isAdvanced) xabort("Calling %s is not allowed without calling cdiIterator_nextField() first.", __func__); \
} while(0)

/**
@Function cdiIterator_clone
@Title Make a copy of an iterator

@Prototype CdiIterator* cdiIterator_clone(CdiIterator* me)
@Parameter
    @item iterator The iterator to copy.

@Result The clone.

@Description
    Clones the given iterator. Make sure to call cdiIterator_delete() on both
    the copy and the original.

    This is not a cheap operation: Depending on the type of the file, it will
    either make a copy of the current field in memory (GRIB files), or reopen
    the file (all other file types). Use it sparingly. And if you do, try to
    avoid keeping too many clones around: their memory footprint is
    significant.
*/
CdiIterator* cdiIterator_clone(CdiIterator* me)
{
  sanityCheck(me);
  switch(me->filetype)
    {
      #ifdef HAVE_LIBGRIB_API
        case FILETYPE_GRB:
        case FILETYPE_GRB2:
          return &cdiGribIterator_clone(me)->super;
      #endif

      #ifdef HAVE_LIBNETCDF
        case FILETYPE_NC:
        case FILETYPE_NC2:
        case FILETYPE_NC4:
        case FILETYPE_NC4C:
      #endif
      #ifdef HAVE_LIBSERVICE
        case FILETYPE_SRV:
      #endif
      #ifdef HAVE_LIBEXTRA
        case FILETYPE_EXT:
      #endif
      #ifdef HAVE_LIBIEG
        case FILETYPE_IEG:
      #endif
          return &cdiFallbackIterator_clone(me)->super;

      default:
        Error(kUnexpectedFileTypeMessage);
        return NULL;
    }
}

/**
@Function cdiGribIterator_clone
@Title Gain access to GRIB specific functionality

@Prototype CdiGribIterator* cdiGribIterator_clone(CdiIterator* me)
@Parameter
    @item iterator The iterator to operate on.

@Result A clone that allows access to GRIB specific functionality, or NULL if the underlying file is not a GRIB file.

@Description
    Clones the given iterator iff the underlying file is a GRIB file, the returned iterator allows access to GRIB specific functionality.
    Make sure to check that the return value is not NULL, and to call cdiGribIterator_delete() on the copy.

    This is not a cheap operation: It will make a copy of the current field in memory. Use it sparingly. And if you do, try to avoid keeping too many clones around, their memory footprint is significant.
*/
CdiGribIterator* cdiGribIterator_clone(CdiIterator* me)
{
  sanityCheck(me);
  switch(me->filetype)
    {
      #ifdef HAVE_LIBGRIB_API
        case FILETYPE_GRB:
        case FILETYPE_GRB2:
          return cdiGribIterator_makeClone(me);
      #endif

      default:
        return NULL;
    }
}

/**
@Function cdiIterator_serialize
@Title Serialize an iterator for sending it to another process

@Prototype char* cdiIterator_serialize(CdiIterator* me)
@Parameter
    @item iterator The iterator to operate on.

@Result A malloc'ed string that contains the full description of the iterator.

@Description
    Make sure to call free() on the resulting string.
*/
char* cdiIterator_serialize(CdiIterator* me)
{
  if(!me) xabort("NULL was passed to %s as an iterator. Please check the return value of cdiIterator_new().", __func__); \
  char* subclassDescription = NULL;
  switch(me->filetype)
    {
      #ifdef HAVE_LIBGRIB_API
        case FILETYPE_GRB:
        case FILETYPE_GRB2:
          subclassDescription = cdiGribIterator_serialize(me);
          break;
      #endif

      #ifdef HAVE_LIBNETCDF
        case FILETYPE_NC:
        case FILETYPE_NC2:
        case FILETYPE_NC4:
        case FILETYPE_NC4C:
      #endif
      #ifdef HAVE_LIBSERVICE
        case FILETYPE_SRV:
      #endif
      #ifdef HAVE_LIBEXTRA
        case FILETYPE_EXT:
      #endif
      #ifdef HAVE_LIBIEG
        case FILETYPE_IEG:
      #endif
          subclassDescription = cdiFallbackIterator_serialize(me);
          break;

      default:
        Error(kUnexpectedFileTypeMessage);
        return NULL;
    }

  const char *ftypeStr = fileType2String(me->filetype),
    *advStr = me->isAdvanced ? kAdvancedString : kUnadvancedString;
  char* result = xmalloc(strlen(ftypeStr) + 1 + strlen(advStr) + 1
                         + strlen(subclassDescription) + 1);
  sprintf(result, "%s %s %s", ftypeStr, advStr, subclassDescription);
  free(subclassDescription);
  return result;
}

/**
@Function cdiIterator_deserialize
@Title Recreate an iterator from its textual description

@Prototype CdiIterator* cdiIterator_deserialize(const char* description)
@Parameter
    @item description The result of a call to cdiIterator_serialize().

@Result A clone of the original iterator.

@Description
    A pair of cdiIterator_serialize() and cdiIterator_deserialize() is functionally equivalent to a call to cdiIterator_clone()

    This function will reread the current field from disk, so don't expect immediate return.
*/
//This only checks the type of the iterator and calls the corresponding subclass function,
//the real deserialization is done in baseIter_constructFromString().
CdiIterator* cdiIterator_deserialize(const char* description)
{
  switch(string2FileType(description, NULL))
    {
      #ifdef HAVE_LIBGRIB_API
        case FILETYPE_GRB:
        case FILETYPE_GRB2:
          return &cdiGribIterator_deserialize(description)->super;
      #endif

      #ifdef HAVE_LIBNETCDF
        case FILETYPE_NC:
        case FILETYPE_NC2:
        case FILETYPE_NC4:
        case FILETYPE_NC4C:
      #endif
      #ifdef HAVE_LIBSERVICE
        case FILETYPE_SRV:
      #endif
      #ifdef HAVE_LIBEXTRA
        case FILETYPE_EXT:
      #endif
      #ifdef HAVE_LIBIEG
        case FILETYPE_IEG:
      #endif
          return &cdiFallbackIterator_deserialize(description)->super;

      default:
        Error(kUnexpectedFileTypeMessage);
        return NULL;
    }
}


/**
@Function cdiIterator_print
@Title Print a textual description of the iterator to a stream

@Prototype void cdiIterator_print(CdiIterator* iterator, FILE* stream);
@Parameter
    @item iterator The iterator to print.
    @item stream The stream to print to.

@Description
    Use for debugging output.
*/
void cdiIterator_print(CdiIterator* me, FILE* stream)
{
  char* description = cdiIterator_serialize(me);
  fprintf(stream, "%s\n", description);
  free(description);
}


/**
@Function cdiIterator_nextField
@Title Advance an iterator to the next field in the file

@Prototype int cdiIterator_nextField(CdiIterator* iterator)
@Parameter
    @item iterator The iterator to operate on.

@Result An error code. May be one of:
  * CDI_NOERR: The iterator has successfully been advanced to the next field.
  * CDI_EEOF: No more fields to read in this file.

@Description
    One call to cdiIterator_nextField() is required before the metadata of the first field can be examined.
    Usually, it will be used directly as the condition for a while() loop.
*/
int cdiIterator_nextField(CdiIterator* me)
{
  if(!me) xabort("NULL was passed in as an iterator. Please check the return value of cdiIterator_new().");
  me->isAdvanced = true;
  switch(me->filetype)
    {
      #ifdef HAVE_LIBGRIB_API
        case FILETYPE_GRB:
        case FILETYPE_GRB2:
          return cdiGribIterator_nextField(me);
      #endif

      #ifdef HAVE_LIBNETCDF
        case FILETYPE_NC:
        case FILETYPE_NC2:
        case FILETYPE_NC4:
        case FILETYPE_NC4C:
      #endif
      #ifdef HAVE_LIBSERVICE
        case FILETYPE_SRV:
      #endif
      #ifdef HAVE_LIBEXTRA
        case FILETYPE_EXT:
      #endif
      #ifdef HAVE_LIBIEG
        case FILETYPE_IEG:
      #endif
          return cdiFallbackIterator_nextField(me);

      default:
        Error(kUnexpectedFileTypeMessage);
        return CDI_EINVAL;
    }
}

static char* cdiIterator_inqTime(CdiIterator* me, bool getEndTime)
{
  sanityCheck(me);
  switch(me->filetype)
    {
      #ifdef HAVE_LIBGRIB_API
        case FILETYPE_GRB:
        case FILETYPE_GRB2:
          return cdiGribIterator_inqTime(me, getEndTime);
      #endif

      #ifdef HAVE_LIBNETCDF
        case FILETYPE_NC:
        case FILETYPE_NC2:
        case FILETYPE_NC4:
        case FILETYPE_NC4C:
      #endif
      #ifdef HAVE_LIBSERVICE
        case FILETYPE_SRV:
      #endif
      #ifdef HAVE_LIBEXTRA
        case FILETYPE_EXT:
      #endif
      #ifdef HAVE_LIBIEG
        case FILETYPE_IEG:
      #endif
          return cdiFallbackIterator_inqTime(me, getEndTime);

      default:
        Error(kUnexpectedFileTypeMessage);
        return NULL;
    }
}

/**
@Function cdiIterator_inqStartTime
@Title Get the start time of a measurement

@Prototype char* cdiIterator_inqStartTime(CdiIterator* me)
@Parameter
    @item iterator The iterator to operate on.

@Result A malloc'ed string containing the (start) time of the current field in the format "YYYY-MM-DDTHH:MM:SS.mmm".

@Description
The returned time is either the time of the data (fields defined at a time point),
or the start time of an integration time range (statistical fields).

Converts the time to the ISO-8601 format and returns it in a newly allocated buffer.
The caller is responsible to free() the resulting string.

If the file is a GRIB file, the calendar that is used to resolve the relative times is the proleptic calendar
as it is implemented by the standard C mktime() function.
This is due to the fact that GRIB-API version 1.12.3 still does not implement the calendar identification fields.
*/
char* cdiIterator_inqStartTime(CdiIterator* me)
{
  return cdiIterator_inqTime(me, false);
}

/**
@Function cdiIterator_inqEndTime
@Title Get the end time of a measurement

@Prototype char* cdiIterator_inqEndTime(CdiIterator* me)
@Parameter
    @item iterator The iterator to operate on.

@Result A malloc'ed string containing the end time of the current field in the format "YYYY-MM-DDTHH:MM:SS.mmm", or NULL if no such time is defined.

@Description
The returned time is the end time of an integration period if such a time exists (statistical fields).
Otherwise, NULL is returned.

Converts the time to the ISO-8601 format and returns it in a newly allocated buffer.
The caller is responsible to free() the resulting string.

If the file is a GRIB file, the calendar that is used to resolve the relative times is the proleptic calendar
as it is implemented by the standard C mktime() function.
This is due to the fact that GRIB-API version 1.12.3 still does not implement the calendar identification fields.
*/
char* cdiIterator_inqEndTime(CdiIterator* me)
{
  return cdiIterator_inqTime(me, true);
}

/**
@Function cdiIterator_inqVTime
@Title Get the validity time of the current field

@Prototype char* cdiIterator_inqVTime(CdiIterator* me)
@Parameter
    @item iterator The iterator to operate on.

@Result A malloc'ed string containing the validity time of the current field in the format "YYYY-MM-DDTHH:MM:SS.mmm".

@Description
The returned time is the validity time as it is returned by taxisInqVtime(), only more precise.
That is, if the field is a time point, its time is returned,
if it is a statistical field with an integration period, the end time of the integration period is returned.

Converts the time to the ISO-8601 format and returns it in a newly allocated buffer.
The caller is responsible to free() the resulting string.

If the file is a GRIB file, the calendar that is used to resolve the relative times is the proleptic calendar
as it is implemented by the standard C mktime() function.
This is due to the fact that GRIB-API version 1.12.3 still does not implement the calendar identification fields.
*/
char* cdiIterator_inqVTime(CdiIterator* me)
{
  char* result = cdiIterator_inqEndTime(me);
  return (result) ? result : cdiIterator_inqStartTime(me);
}

/**
@Function cdiIterator_inqLevelType
@Title Get the type of a level

@Prototype int cdiIterator_inqLevelType(CdiIterator* me, int levelSelector, char** outName = NULL, char** outLongName = NULL, char** outStdName = NULL, char** outUnit = NULL)
@Parameter
    @item iterator The iterator to operate on.
    @item levelSelector Zero for the top level, one for the bottom level
    @item outName Will be set to a malloc()'ed string with the name of the level if not NULL.
    @item outLongName Will be set to a malloc()'ed string with the long name of the level if not NULL.
    @item outStdName Will be set to a malloc()'ed string with the standard name of the level if not NULL.
    @item outUnit Will be set to a malloc()'ed string with the unit of the level if not NULL.

@Result An integer indicating the type of the level.

@Description
Find out some basic information about the given level, the levelSelector selects the function of the requested level.
If the requested level does not exist, this returns CDI_UNDEFID.
*/
int cdiIterator_inqLevelType(CdiIterator* me, int levelSelector, char** outName, char** outLongName, char** outStdName, char** outUnit)
{
  sanityCheck(me);
  switch(me->filetype)
    {
      #ifdef HAVE_LIBGRIB_API
        case FILETYPE_GRB:
        case FILETYPE_GRB2:
          return cdiGribIterator_levelType(me, levelSelector, outName, outLongName, outStdName, outUnit);
      #endif

      #ifdef HAVE_LIBNETCDF
        case FILETYPE_NC:
        case FILETYPE_NC2:
        case FILETYPE_NC4:
        case FILETYPE_NC4C:
      #endif
      #ifdef HAVE_LIBSERVICE
        case FILETYPE_SRV:
      #endif
      #ifdef HAVE_LIBEXTRA
        case FILETYPE_EXT:
      #endif
      #ifdef HAVE_LIBIEG
        case FILETYPE_IEG:
      #endif
          return cdiFallbackIterator_levelType(me, levelSelector, outName, outLongName, outStdName, outUnit);

      default:
        Error(kUnexpectedFileTypeMessage);
        return CDI_UNDEFID;
    }
}

/**
@Function cdiIterator_inqLevel
@Title Get the value of the z-coordinate

@Prototype void cdiIterator_inqLevel(CdiIterator* me, int levelSelector, double* outValue1, double* outValue2 = NULL)
@Parameter
    @item iterator The iterator to operate on.
    @item levelSelector Zero for the top level, one for the bottom level
    @item outValue1 For "normal" levels this returns the value, for hybrid levels the first coordinate, for generalized levels the level number.
    @item outValue2 Zero for "normal" levels, for hybrid levels, this returns the second coordinate, for generalized levels the level count.

@Result An error code.

@Description
Returns the value of the z-coordinate, whatever that may be.
*/
int cdiIterator_inqLevel(CdiIterator* me, int levelSelector, double* outValue1, double* outValue2)
{
  sanityCheck(me);
  switch(me->filetype)
    {
      #ifdef HAVE_LIBGRIB_API
        case FILETYPE_GRB:
        case FILETYPE_GRB2:
          return cdiGribIterator_level(me, levelSelector, outValue1, outValue2);
      #endif

      #ifdef HAVE_LIBNETCDF
        case FILETYPE_NC:
        case FILETYPE_NC2:
        case FILETYPE_NC4:
        case FILETYPE_NC4C:
      #endif
      #ifdef HAVE_LIBSERVICE
        case FILETYPE_SRV:
      #endif
      #ifdef HAVE_LIBEXTRA
        case FILETYPE_EXT:
      #endif
      #ifdef HAVE_LIBIEG
        case FILETYPE_IEG:
      #endif
          return cdiFallbackIterator_level(me, levelSelector, outValue1, outValue2);

      default:
        Error(kUnexpectedFileTypeMessage);
        return CDI_EINVAL;
    }
}

/**
@Function cdiIterator_inqLevelUuid
@Title Get the UUID of the z-axis used by this field

@Prototype int cdiIterator_inqLevelUuid(CdiIterator* me, int levelSelector, unsigned char (*outUuid)[16])
@Parameter
    @item iterator The iterator to operate on.
    @item outVgridNumber The number of the associated vertical grid description.
    @item outLevelCount The number of levels in the associated vertical grid description.
    @item outUuid A pointer to a user supplied buffer of 16 bytes to store the UUID in.

@Result An error code.

@Description
Returns identifying information for the external z-axis description. May only be called for generalized levels.
*/
int cdiIterator_inqLevelUuid(CdiIterator* me, int* outVgridNumber, int* outLevelCount, unsigned char (*outUuid)[16])
{
  sanityCheck(me);
  switch(me->filetype)
    {
      #ifdef HAVE_LIBGRIB_API
        case FILETYPE_GRB:
        case FILETYPE_GRB2:
          return cdiGribIterator_zaxisUuid(me, outVgridNumber, outLevelCount, outUuid);
      #endif

      #ifdef HAVE_LIBNETCDF
        case FILETYPE_NC:
        case FILETYPE_NC2:
        case FILETYPE_NC4:
        case FILETYPE_NC4C:
      #endif
      #ifdef HAVE_LIBSERVICE
        case FILETYPE_SRV:
      #endif
      #ifdef HAVE_LIBEXTRA
        case FILETYPE_EXT:
      #endif
      #ifdef HAVE_LIBIEG
        case FILETYPE_IEG:
      #endif
          return cdiFallbackIterator_zaxisUuid(me, outVgridNumber, outLevelCount, outUuid);

      default:
        Error(kUnexpectedFileTypeMessage);
        return CDI_ELIBNAVAIL;
    }
}

/**
@Function cdiIterator_inqParam
@Title Get discipline, category, and number

@Prototype CdiParam cdiIterator_inqParam(CdiIterator* iterator)
@Parameter
    @item iterator The iterator to operate on.

@Result A struct containing the requested information.

@Description
    Simple metadata inspection function.
*/
CdiParam cdiIterator_inqParam(CdiIterator* me)
{
  sanityCheck(me);
  return me->param;
}

/**
@Function cdiIterator_inqDatatype
@Title Get the datatype of the current field

@Prototype int cdiIterator_inqDatatype(CdiIterator* iterator)
@Parameter
    @item iterator The iterator to operate on.

@Result The datatype that is used to store this field on disk.

@Description
    Simple metadata inspection function.
*/
int cdiIterator_inqDatatype(CdiIterator* me)
{
  sanityCheck(me);
  return me->datatype;
}

/**
@Function cdiIterator_inqTsteptype
@Title Get the timestep type

@Prototype int cdiIterator_inqTsteptype(CdiIterator* iterator)
@Parameter
    @item iterator The iterator to operate on.

@Result The timestep type.

@Description
    Simple metadata inspection function.
*/
int cdiIterator_inqTsteptype(CdiIterator* me)
{
  sanityCheck(me);
  return me->timesteptype;
}

/**
@Function cdiIterator_inqVariableName
@Title Get the variable name of the current field

@Prototype char* cdiIterator_inqVariableName(CdiIterator* iterator)
@Parameter
    @item iterator The iterator to operate on.

@Result A pointer to a C-string containing the name. The storage for this string is allocated with malloc(), and it is the responsibility of the caller to free() it.

@Description
    Allocates a buffer to hold the string, copies the current variable name into this buffer, and returns the buffer.
    The caller is responsible to make the corresponding free() call.
*/
char* cdiIterator_inqVariableName(CdiIterator* me)
{
  sanityCheck(me);
  switch(me->filetype)
    {
      #ifdef HAVE_LIBGRIB_API
        case FILETYPE_GRB:
        case FILETYPE_GRB2:
          return cdiGribIterator_copyVariableName(me);
      #endif

      #ifdef HAVE_LIBNETCDF
        case FILETYPE_NC:
        case FILETYPE_NC2:
        case FILETYPE_NC4:
        case FILETYPE_NC4C:
      #endif
      #ifdef HAVE_LIBSERVICE
        case FILETYPE_SRV:
      #endif
      #ifdef HAVE_LIBEXTRA
        case FILETYPE_EXT:
      #endif
      #ifdef HAVE_LIBIEG
        case FILETYPE_IEG:
      #endif
          return cdiFallbackIterator_copyVariableName(me);

      default:
        Error(kUnexpectedFileTypeMessage);
        return NULL;
    }
}

/**
@Function cdiIterator_inqGridId
@Title Get the ID of the current grid

@Prototype int cdiIterator_inqGridId(CdiIterator* iterator)
@Parameter
    @item iterator The iterator to operate on.

@Result A gridId that can be used for further introspection.

@Description
    This provides access to the grid related metadata.
    The resulting ID is only valid until the next time cdiIterator_nextField() is called.
*/
int cdiIterator_inqGridId(CdiIterator* me)
{
  sanityCheck(me);
  return me->gridId;
}

/**
@Function cdiIterator_readField
@Title Read the whole field into a double buffer

@Prototype void cdiIterator_readField(CdiIterator* me, double* buffer, size_t* nmiss)
@Parameter
    @item iterator The iterator to operate on.
    @item buffer A pointer to the double array that the data should be written to.
    @item nmiss A pointer to a variable where the count of missing values will be stored. May be NULL.

@Description
    It is assumed that the caller first analyses the return value of cdiIterator_inqGridId to determine the required size of the buffer.
    Failing to do so results in undefined behavior. You have been warned.
*/
void cdiIterator_readField(CdiIterator* me, double* buffer, size_t* nmiss)
{
  sanityCheck(me);
  if(!buffer) xabort("NULL was passed in a buffer. Please provide a suitable buffer.");
  switch(me->filetype)
    {
      #ifdef HAVE_LIBGRIB_API
        case FILETYPE_GRB:
        case FILETYPE_GRB2:
          cdiGribIterator_readField(me, buffer, nmiss);
	  return;
      #endif

      #ifdef HAVE_LIBNETCDF
        case FILETYPE_NC:
        case FILETYPE_NC2:
        case FILETYPE_NC4:
        case FILETYPE_NC4C:
      #endif
      #ifdef HAVE_LIBSERVICE
        case FILETYPE_SRV:
      #endif
      #ifdef HAVE_LIBEXTRA
        case FILETYPE_EXT:
      #endif
      #ifdef HAVE_LIBIEG
        case FILETYPE_IEG:
      #endif
          cdiFallbackIterator_readField(me, buffer, nmiss);
          return;
      default:
        Error(kUnexpectedFileTypeMessage);
    }
}

/**
@Function cdiIterator_readFieldF
@Title Read the whole field into a double buffer

@Prototype void cdiIterator_readFieldF(CdiIterator* me, float* buffer, size_t* nmiss)
@Parameter
    @item iterator The iterator to operate on.
    @item buffer A pointer to the double array that the data should be written to.
    @item nmiss A pointer to a variable where the count of missing values will be stored. May be NULL.

@Description
    It is assumed that the caller first analyses the return value of cdiIterator_inqGridId to determine the required size of the buffer.
    Failing to do so results in undefined behavior. You have been warned.
*/
void cdiIterator_readFieldF(CdiIterator* me, float* buffer, size_t* nmiss)
{
  sanityCheck(me);
  if(!buffer) xabort("NULL was passed in a buffer. Please provide a suitable buffer.");
  switch(me->filetype)
    {
      #ifdef HAVE_LIBGRIB_API
        case FILETYPE_GRB:
        case FILETYPE_GRB2:
          cdiGribIterator_readFieldF(me, buffer, nmiss);
	  return;
      #endif

      #ifdef HAVE_LIBNETCDF
        case FILETYPE_NC:
        case FILETYPE_NC2:
        case FILETYPE_NC4:
        case FILETYPE_NC4C:
      #endif
      #ifdef HAVE_LIBSERVICE
        case FILETYPE_SRV:
      #endif
      #ifdef HAVE_LIBEXTRA
        case FILETYPE_EXT:
      #endif
      #ifdef HAVE_LIBIEG
        case FILETYPE_IEG:
      #endif
          cdiFallbackIterator_readFieldF(me, buffer, nmiss);
          return; 
      default:
        Error(kUnexpectedFileTypeMessage);
    }
}

/**
@Function cdiIterator_delete
@Title Destroy an iterator

@Prototype void cdiIterator_delete(CdiIterator* iterator)
@Parameter
    @item iterator The iterator to operate on.

@Description
    Combined destructor & deallocator.
*/
void cdiIterator_delete(CdiIterator* me)
{
  if(!me) xabort("NULL was passed in as an iterator. Please check the return value of cdiIterator_new().");
  switch(me->filetype)
    {
      #ifdef HAVE_LIBGRIB_API
        case FILETYPE_GRB:
        case FILETYPE_GRB2:
          cdiGribIterator_delete((CdiGribIterator*)me);
          break;
      #endif

      #ifdef HAVE_LIBNETCDF
        case FILETYPE_NC:
        case FILETYPE_NC2:
        case FILETYPE_NC4:
        case FILETYPE_NC4C:
      #endif
      #ifdef HAVE_LIBSERVICE
        case FILETYPE_SRV:
      #endif
      #ifdef HAVE_LIBEXTRA
        case FILETYPE_EXT:
      #endif
      #ifdef HAVE_LIBIEG
        case FILETYPE_IEG:
      #endif
          cdiFallbackIterator_delete(me);
          break;

      default:
        Error(kUnexpectedFileTypeMessage);
    }
}

void baseIterDestruct(CdiIterator* me)
{
  /*currently empty, but that's no reason not to call it*/
  (void)me;
}
