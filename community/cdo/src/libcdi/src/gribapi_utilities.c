#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#ifdef HAVE_LIBGRIB_API

#include "gribapi_utilities.h"

#include "cdi.h"
#include "dmemory.h"
#include "error.h"
#include "gribapi.h"

#include <assert.h>
#include <time.h>

#define FAIL_ON_GRIB_ERROR(function, gribHandle, key, ...) do\
{\
  int errorCode = (int)function(gribHandle, key, __VA_ARGS__);  \
  if(errorCode)\
    {\
      fprintf(stderr, "%s:%d: Error in function `%s`: `%s` returned error code %d for key \"%s\"", __FILE__, __LINE__, __func__, #function, errorCode, key);\
      exit(errorCode);\
    }\
} while(0)

//A simple wrapper for grib_get_string() that returns a newly allocated string.
char* gribCopyString(grib_handle* gribHandle, const char* key)
{
  char* result = NULL;
  size_t length;
#ifdef HAVE_GRIB_GET_LENGTH
  if(!grib_get_length(gribHandle, key, &length))
  {
    char* result = xmalloc(length);
    if(!grib_get_string(gribHandle, key, result, &length))
      result = xrealloc(result, length);
    else
    {
      free(result);
      result = NULL;
    }
  }
#else
  length = 1024;         /* there's an implementation limit
                          * that makes strings longer than
                          * this unlikely in grib_api versions
                          * not providing grib_get_length */
  int rc;
  result = xmalloc(length);
  while ((rc = grib_get_string(gribHandle, key, result, &length))
         == GRIB_BUFFER_TOO_SMALL || rc == GRIB_ARRAY_TOO_SMALL)
  {
    if (length <= 1024UL * 1024UL)
    {
      length *= 2;
      result = xrealloc(result, length);
    }
    else
      break;
  }
  if (!rc)
  {
    result = xrealloc(result, length);
    return result;
  }
  free(result);
#endif
  return NULL;
}

//A simple wrapper for grib_get_string() for the usecase that the result is only compared to a given constant string.
//Returns true if the key exists and the value is equal to the given string.
bool gribCheckString(grib_handle* gribHandle, const char* key, const char* expectedValue)
{
  size_t expectedLength = strlen(expectedValue) + 1;
#ifdef HAVE_GRIB_GET_LENGTH
  size_t length;
  if(grib_get_length(gribHandle, key, &length)) return false;
  if(length != expectedLength) return false;
  char *value = xmalloc(length);
  if(grib_get_string(gribHandle, key, value, &length)) return false;
  int rc = !strcmp(value, expectedValue);
  free(value);
#else
  char *value = gribCopyString(gribHandle, key);
  int rc;
  if (value)
  {
    rc = strlen(value) + 1 == expectedLength ?
      !strcmp(value, expectedValue)
      : false;
  }
  else
    rc = false;
  free(value);
#endif
  return rc;
}

//A simple wrapper for grib_get_long() for the usecase that the result is only compared to a given constant value.
//Returns true if the key exists and the value is equal to the given one.
bool gribCheckLong(grib_handle* gribHandle, const char* key, long expectedValue)
{
  long value;
  if(grib_get_long(gribHandle, key, &value)) return false;
  return value == expectedValue;
}

//A simple wrapper for grib_get_long() for the usecase that failure to fetch the value is fatal.
long gribGetLong(grib_handle* gh, const char* key)
{
  long result;
  FAIL_ON_GRIB_ERROR(grib_get_long, gh, key, &result);
  return result;
}

//A simple wrapper for grib_get_long() for the usecase that a default value is used in the case that the operation fails.
long gribGetLongDefault(grib_handle* gribHandle, const char* key, long defaultValue)
{
  long result;
  if(grib_get_long(gribHandle, key, &result)) return defaultValue;
  if(result == GRIB_MISSING_LONG) return defaultValue;
  return result;
}

//A simple wrapper for grib_get_double() for the usecase that failure to fetch the value is fatal.
double gribGetDouble(grib_handle* gh, const char* key)
{
  double result;
  FAIL_ON_GRIB_ERROR(grib_get_double, gh, key, &result);
  return result;
}

//A sample wrapper for grib_get_double() for the usecase that a default value is used in the case that the operation fails.
double gribGetDoubleDefault(grib_handle* gribHandle, const char* key, double defaultValue)
{
  double result;
  if(grib_get_double(gribHandle, key, &result)) return defaultValue;
  if(IS_EQUAL(result, GRIB_MISSING_DOUBLE)) return defaultValue;
  return result;
}

//A simple wrapper for grib_get_size() for the usecase that failure to fetch the value is fatal.
size_t gribGetArraySize(grib_handle* gribHandle, const char* key)
{
  size_t result;
  FAIL_ON_GRIB_ERROR(grib_get_size, gribHandle, key, &result);
  return result;
}

//A simple wrapper for grib_get_double_array() for the usecase that failure to fetch the data is fatal.
void gribGetDoubleArray(grib_handle* gribHandle, const char* key, double* array)
{
  size_t valueCount = gribGetArraySize(gribHandle, key);
  FAIL_ON_GRIB_ERROR(grib_get_double_array, gribHandle, key, array, &valueCount);
}

//A simple wrapper for grib_get_long_array() for the usecase that failure to fetch the data is fatal.
void gribGetLongArray(grib_handle* gribHandle, const char* key, long* array)
{
  size_t valueCount = gribGetArraySize(gribHandle, key);
  FAIL_ON_GRIB_ERROR(grib_get_long_array, gribHandle, key, array, &valueCount);
}


//We need the edition number so frequently, that it's convenient to give it its own function.
long gribEditionNumber(grib_handle* gh)
{
  return gribGetLong(gh, "editionNumber");
}

//This return value of this should be passed to a call to resetTz(), it is a malloc'ed string with the content of the TZ environment variable before the call (or NULL if that was not set).
static char* setUtc()
{
  char* temp = getenv("TZ"), *result = NULL;
  if(temp) result = strdup(temp);
  setenv("TZ", "UTC", 1);
  return result;
}

//Undoes the effect of setUtc(), pass to it the return value of the corresponding setUtc() call, it will free the string.
static void resetTz(char* savedTz)
{
  if(savedTz)
    {
      setenv("TZ", savedTz, 1);
      free(savedTz);
    }
  else
    {
      unsetenv("TZ");
    }
}

//This function uses the system functions to normalize the date representation according to the gregorian calendar.
//Returns zero on success.
static int normalizeDays(struct tm* me)
{
  char* savedTz = setUtc();     //Ensure that mktime() does not interprete the date according to our local time zone.

  int result = mktime(me) == (time_t)-1;        //This does all the heavy lifting.

  resetTz(savedTz);
  return result;
}

//Returns zero on success.
static int addSecondsToDate(struct tm* me, long long amount)
{
  //It is irrelevant here whether days are zero or one based, the correction would have be undone again so that it is effectless.
  long long seconds = ((me->tm_mday*24ll + me->tm_hour)*60 + me->tm_min)*60 + me->tm_sec;    //The portion of the date that uses fixed increments.
  seconds += amount;
  me->tm_mday = (int)(seconds / 24 / 60 / 60);
  seconds -= (long long)me->tm_mday * 24 * 60 * 60;
  me->tm_hour = (int)(seconds / 60 / 60);
  seconds -= (long long)me->tm_hour * 60 * 60;
  me->tm_min = (int)(seconds / 60);
  seconds -= (long long)(me->tm_min * 60);
  me->tm_sec = (int)seconds;
  return normalizeDays(me);
}

static void addMonthsToDate(struct tm* me, long long amount)
{
  long long months = me->tm_year*12ll + me->tm_mon;
  months += amount;
  me->tm_year = (int)(months/12);
  months -= (long long)me->tm_year*12;
  me->tm_mon = (int)months;
}

//unit is a value according to code table 4.4 of the GRIB2 specification, returns non-zero on error
static int addToDate(struct tm* me, long long amount, long unit)
{
  switch(unit)
    {
      case 0: return addSecondsToDate(me,       60*amount);   // minute
      case 1: return addSecondsToDate(me,    60*60*amount);   // hour
      case 2: return addSecondsToDate(me, 24*60*60*amount);   // day

      case 3: addMonthsToDate(me,        amount); return 0;   // month
      case 4: addMonthsToDate(me,     12*amount); return 0;   // year
      case 5: addMonthsToDate(me,  10*12*amount); return 0;   // decade
      case 6: addMonthsToDate(me,  30*12*amount); return 0;   // normal
      case 7: addMonthsToDate(me, 100*12*amount); return 0;   // century

      case 10: return addSecondsToDate(me,  3*60*60*amount);  // eighth of a day
      case 11: return addSecondsToDate(me,  6*60*60*amount);  // quarter day
      case 12: return addSecondsToDate(me, 12*60*60*amount);  // half day
      case 13: return addSecondsToDate(me,          amount);  // second

      default: return 1;        //reserved, unknown, or missing
    }
}

static char* makeDateString(struct tm* me)
{
  char *result
    = xmalloc(       4+1+ 2+1+ 2+1+ 2+1+ 2+1+ 2+ 4+ 1);
  sprintf(result, "%04d-%02d-%02dT%02d:%02d:%02d.000", me->tm_year + 1900, me->tm_mon + 1, me->tm_mday, me->tm_hour, me->tm_min, me->tm_sec);
  return result;
}

//FIXME: This ignores any calendar definition that might be present.
//XXX: Identification templates are not implemented in grib_api-1.12.3, so even if I implemented the other calendars now, it wouldn't be possible to use them.
static int getAvailabilityOfRelativeTimes(grib_handle* gh, bool* outHaveForecastTime, bool* outHaveTimeRange)
{
  switch(gribGetLong(gh, "productDefinitionTemplateNumber"))
    {
      case 20: case 30: case 31: case 254: case 311: case 2000:
        *outHaveForecastTime = false, *outHaveTimeRange = false;
        return 0;

      case 0: case 1: case 2: case 3: case 4: case 5: case 6: case 7: case 15: case 32: case 33: case 40: case 41: case 44: case 45: case 48: case 51: case 53: case 54: case 60: case 1000: case 1002: case 1100: case 40033:
        *outHaveForecastTime = true, *outHaveTimeRange = false;
        return 0;

      case 8: case 9: case 10: case 11: case 12: case 13: case 14: case 34: case 42: case 43: case 46: case 47: case 61: case 91: case 1001: case 1101: case 40034:
        *outHaveForecastTime = true, *outHaveTimeRange = true;
        return 0;

      default:
        return 1;
    }
}

char* gribMakeTimeString(grib_handle* gh, bool getEndTime)
{
  //Get the parts of the reference date.
  struct tm date = {
    .tm_mon = (int)gribGetLong(gh, "month") - 1,   //months are zero based in struct tm and one based in GRIB
    .tm_mday = (int)gribGetLong(gh, "day"),
    .tm_hour = (int)gribGetLong(gh, "hour"),
    .tm_min = (int)gribGetLong(gh, "minute")
  };
  if(gribEditionNumber(gh) == 1)
    {
      date.tm_year = (int)gribGetLong(gh, "yearOfCentury");  //years are -1900 based both in struct tm and GRIB1
    }
  else
    {
      date.tm_year = (int)gribGetLong(gh, "year") - 1900;   //years are -1900 based in struct tm and zero based in GRIB2
      date.tm_sec = (int)gribGetLong(gh, "second");

      //Determine whether we have a forecast time and a time range.
      bool haveForecastTime, haveTimeRange;
      if(getAvailabilityOfRelativeTimes(gh, &haveForecastTime, &haveTimeRange)) return NULL;
      if(getEndTime && !haveTimeRange) return NULL;     //tell the caller that the requested time does not exist

      //If we have relative times, apply them to the date
      if(haveForecastTime)
        {
          long offset = gribGetLongDefault(gh, "forecastTime", 0);  //if(stepUnits == indicatorOfUnitOfTimeRange) assert(startStep == forecastTime)
          long offsetUnit = gribGetLongDefault(gh, "indicatorOfUnitOfTimeRange", 255);
          if(addToDate(&date, offset, offsetUnit)) return NULL;
          if(getEndTime)
            {
              assert(haveTimeRange);
              long range = gribGetLongDefault(gh, "lengthOfTimeRange", 0);       //if(stepUnits == indicatorOfUnitForTimeRange) assert(endStep == startStep + lengthOfTimeRange)
              long rangeUnit = gribGetLongDefault(gh, "indicatorOfUnitForTimeRange", 255);
              if(addToDate(&date, range, rangeUnit)) return NULL;
            }
        }
    }

  //Bake the date into a string.
  return makeDateString(&date);
}

int gribapiTimeIsFC(grib_handle *gh)
{
  if(gribEditionNumber(gh) <= 1) return true;

  long sigofrtime;
  FAIL_ON_GRIB_ERROR(grib_get_long, gh, "significanceOfReferenceTime", &sigofrtime);
  return sigofrtime != 3;
}

//Fetches the value of the "stepType" key and converts it into a constant in the TSTEP_* range.
int gribapiGetTsteptype(grib_handle *gh)
{
  int tsteptype = TSTEP_INSTANT;
  static bool lprint = true;

  if ( gribapiTimeIsFC(gh) )
    {
      int status;
      size_t len = 256;
      char stepType[256];

      status = grib_get_string(gh, "stepType", stepType, &len);
      if ( status == 0 && len > 1 && len < 256 )
        {
          if      ( strncmp("instant", stepType, len) == 0 ) tsteptype = TSTEP_INSTANT;
          else if ( strncmp("avg",     stepType, len) == 0 ) tsteptype = TSTEP_AVG;
          else if ( strncmp("accum",   stepType, len) == 0 ) tsteptype = TSTEP_ACCUM;
          else if ( strncmp("max",     stepType, len) == 0 ) tsteptype = TSTEP_MAX;
          else if ( strncmp("min",     stepType, len) == 0 ) tsteptype = TSTEP_MIN;
          else if ( strncmp("diff",    stepType, len) == 0 ) tsteptype = TSTEP_DIFF;
          else if ( strncmp("rms",     stepType, len) == 0 ) tsteptype = TSTEP_RMS;
          else if ( strncmp("sd",      stepType, len) == 0 ) tsteptype = TSTEP_SD;
          else if ( strncmp("cov",     stepType, len) == 0 ) tsteptype = TSTEP_COV;
          else if ( strncmp("ratio",   stepType, len) == 0 ) tsteptype = TSTEP_RATIO;
          else if ( lprint )
            {
              Message("Time stepType %s unsupported, set to instant!", stepType);
              lprint = false;
            }

          // printf("stepType: %s %ld %d\n", stepType, len, tsteptype);
        }
    }

  return (tsteptype);
}

int gribGetDatatype(grib_handle* gribHandle)
{
  int datatype;
  if(gribEditionNumber(gribHandle) > 1 && gribCheckString(gribHandle, "packingType", "grid_ieee"))
    {
      datatype = gribCheckLong(gribHandle, "precision", 1) ? DATATYPE_FLT32 : DATATYPE_FLT64;
    }
  else
    {
      long bitsPerValue;
      datatype = (!grib_get_long(gribHandle, "bitsPerValue", &bitsPerValue) && bitsPerValue > 0 && bitsPerValue <= 32) ? (int)bitsPerValue : DATATYPE_PACK;
    }
  return datatype;
}

int gribapiGetParam(grib_handle *gh)
{
  long pdis, pcat, pnum;
  if ( gribEditionNumber(gh) <= 1 )
    {
      pdis = 255;
      FAIL_ON_GRIB_ERROR(grib_get_long, gh, "table2Version", &pcat);
      FAIL_ON_GRIB_ERROR(grib_get_long, gh, "indicatorOfParameter", &pnum);
    }
  else
    {
      FAIL_ON_GRIB_ERROR(grib_get_long, gh, "discipline", &pdis);
      if(grib_get_long(gh, "parameterCategory", &pcat)) pcat = 0;
      if(grib_get_long(gh, "parameterNumber", &pnum)) pnum = 0;
    }
  return cdiEncodeParam((int)pnum, (int)pcat, (int)pdis);
}

int gribapiGetGridType(grib_handle *gh)
{
  int gridtype = GRID_GENERIC;
  switch (gribGetLongDefault(gh, "gridDefinitionTemplateNumber", -1))
    {
      case  GRIB2_GTYPE_LATLON:
        gridtype = ( gribGetLong(gh, "Ni") == (long) GRIB_MISSING_LONG ) ? GRID_GENERIC : GRID_LONLAT;
        break;

      case  GRIB2_GTYPE_GAUSSIAN:
        gridtype = ( gribGetLong(gh, "Ni") == (long) GRIB_MISSING_LONG ) ? GRID_GAUSSIAN_REDUCED : GRID_GAUSSIAN;
        break;

      case GRIB2_GTYPE_LATLON_ROT:   gridtype = GRID_LONLAT; break;
      case GRIB2_GTYPE_LCC:          gridtype = GRID_LCC; break;
      case GRIB2_GTYPE_SPECTRAL:     gridtype = GRID_SPECTRAL; break;
      case GRIB2_GTYPE_GME:          gridtype = GRID_GME; break;
      case GRIB2_GTYPE_UNSTRUCTURED: gridtype = GRID_UNSTRUCTURED; break;
    }

  return gridtype;
}

static
int gribapiGetIsRotated(grib_handle *gh)
{
  return gribGetLongDefault(gh, "gridDefinitionTemplateNumber", -1) == GRIB2_GTYPE_LATLON_ROT;
}

//TODO: Simplify by use of the convenience functions (gribGetLong(), gribGetLongDefault(), etc.).
void gribapiGetGrid(grib_handle *gh, grid_t *grid)
{
  long editionNumber = gribEditionNumber(gh);
  int gridtype = gribapiGetGridType(gh);
  /*
  if ( streamptr->unreduced && gridtype == GRID_GAUSSIAN_REDUCED )
    {
      gridtype = GRID_GAUSSIAN;
      ISEC2_NumLon = 2*ISEC2_NumLat;
      ISEC4_NumValues = ISEC2_NumLon*ISEC2_NumLat;
    }
  */
  memset(grid, 0, sizeof(grid_t));

  size_t datasize;
  FAIL_ON_GRIB_ERROR(grib_get_size, gh, "values", &datasize);
  long numberOfPoints;
  FAIL_ON_GRIB_ERROR(grib_get_long, gh, "numberOfPoints", &numberOfPoints);

  switch (gridtype)
    {
    case GRID_LONLAT:
    case GRID_GAUSSIAN:
      {
        long lpar;
        FAIL_ON_GRIB_ERROR(grib_get_long, gh, "Ni", &lpar);
        /* FIXME: assert(lpar <= INT_MAX && lpar >= INT_MIN) */
        int nlon = (int)lpar;
        FAIL_ON_GRIB_ERROR(grib_get_long, gh, "Nj", &lpar);
        /* FIXME: assert(lpar <= INT_MAX && lpar >= INT_MIN) */
        int nlat = (int)lpar;

        if ( gridtype == GRID_GAUSSIAN )
          {
            FAIL_ON_GRIB_ERROR(grib_get_long, gh, "numberOfParallelsBetweenAPoleAndTheEquator", &lpar);
            grid->np = (int)lpar;
          }

        if ( numberOfPoints != nlon*nlat )
          Error("numberOfPoints (%ld) and gridSize (%d) differ!", numberOfPoints, nlon*nlat);

        /* FIXME: assert(numberOfPoints <= INT_MAX && numberOfPoints >= INT_MIN) */
        grid->size  = (int)numberOfPoints;
        grid->xsize = nlon;
        grid->ysize = nlat;
        grid->xinc  = 0;
        grid->yinc  = 0;
        grid->xdef  = 0;
        FAIL_ON_GRIB_ERROR(grib_get_double, gh, "longitudeOfFirstGridPointInDegrees", &grid->xfirst);
        FAIL_ON_GRIB_ERROR(grib_get_double, gh, "longitudeOfLastGridPointInDegrees",  &grid->xlast);
        FAIL_ON_GRIB_ERROR(grib_get_double, gh, "latitudeOfFirstGridPointInDegrees",  &grid->yfirst);
        FAIL_ON_GRIB_ERROR(grib_get_double, gh, "latitudeOfLastGridPointInDegrees",   &grid->ylast);
        FAIL_ON_GRIB_ERROR(grib_get_double, gh, "iDirectionIncrementInDegrees", &grid->xinc);
        if ( gridtype == GRID_LONLAT )
          FAIL_ON_GRIB_ERROR(grib_get_double, gh, "jDirectionIncrementInDegrees", &grid->yinc);

        if ( grid->xinc < -999 || grid->xinc > 999 ) grid->xinc = 0;
        if ( grid->yinc < -999 || grid->yinc > 999 ) grid->yinc = 0;

        /* if ( IS_NOT_EQUAL(grid->xfirst, 0) || IS_NOT_EQUAL(grid->xlast, 0) ) */
          {
            if ( grid->xsize > 1 )
              {
                if ( (grid->xfirst >= grid->xlast) && (grid->xfirst >= 180) ) grid->xfirst -= 360;

                if ( editionNumber <= 1 )
                  {
                    /* correct xinc if necessary */
                    if ( IS_EQUAL(grid->xfirst, 0) && grid->xlast > 354 )
                      {
                        double xinc = 360. / grid->xsize;

                        if ( fabs(grid->xinc-xinc) > 0.0 )
                          {
                            grid->xinc = xinc;
                            if ( CDI_Debug ) Message("set xinc to %g", grid->xinc);
                          }
                      }
                  }
              }
            grid->xdef = 2;
          }
        grid->ydef = 0;
        /* if ( IS_NOT_EQUAL(grid->yfirst, 0) || IS_NOT_EQUAL(grid->ylast, 0) ) */
          {
            if ( grid->ysize > 1 )
              {
                if ( editionNumber <= 1 )
                  {
                  }
              }
            grid->ydef = 2;
          }
        break;
      }
    case GRID_GAUSSIAN_REDUCED:
      {
        size_t dummy;
        long *pl;

        long lpar;
        FAIL_ON_GRIB_ERROR(grib_get_long, gh, "numberOfParallelsBetweenAPoleAndTheEquator", &lpar);
        /* FIXME: assert(lpar <= INT_MAX && lpar >= INT_MIN) */
        grid->np = (int)lpar;

        FAIL_ON_GRIB_ERROR(grib_get_long, gh, "Nj", &lpar);
        /* FIXME: assert(lpar <= INT_MAX && lpar >= INT_MIN) */
        int nlat = (int)lpar;

        /* FIXME: assert(numberOfPoints <= INT_MAX && numberOfPoints >= INT_MIN) */
        grid->size   = (int)numberOfPoints;

        grid->rowlon = (int *) malloc((size_t)nlat * sizeof (int));
        pl          = (long *) malloc((size_t)nlat * sizeof (long));
        dummy       = (size_t)nlat;
        FAIL_ON_GRIB_ERROR(grib_get_long_array, gh, "pl", pl, &dummy);
        /* FIXME: assert(pl[i] >= INT_MIN && pl[i] <= INT_MIN) */
        for (int i = 0; i < nlat; ++i ) grid->rowlon[i] = (int)pl[i];
        free(pl);

        grid->ysize  = nlat;
        grid->xinc   = 0;
        grid->yinc   = 0;
        grid->xdef   = 0;
        FAIL_ON_GRIB_ERROR(grib_get_double, gh, "longitudeOfFirstGridPointInDegrees", &grid->xfirst);
        FAIL_ON_GRIB_ERROR(grib_get_double, gh, "longitudeOfLastGridPointInDegrees",  &grid->xlast);
        FAIL_ON_GRIB_ERROR(grib_get_double, gh, "latitudeOfFirstGridPointInDegrees",  &grid->yfirst);
        FAIL_ON_GRIB_ERROR(grib_get_double, gh, "latitudeOfLastGridPointInDegrees",   &grid->ylast);
        FAIL_ON_GRIB_ERROR(grib_get_double, gh, "iDirectionIncrementInDegrees", &grid->xinc);

        if ( IS_EQUAL(grid->xinc, GRIB_MISSING_DOUBLE) ) grid->xinc = 0;

        /* if ( IS_NOT_EQUAL(grid->xfirst, 0) || IS_NOT_EQUAL(grid->xlast, 0) ) */
          {
            if ( grid->xsize > 1 )
              {
                if ( (grid->xfirst > grid->xlast) && (grid->xfirst >= 180) ) grid->xfirst -= 360;

                if ( editionNumber <= 1 )
                  {
                    /* correct xinc if necessary */
                    if ( IS_EQUAL(grid->xfirst, 0) && grid->xlast > 354 )
                      {
                        double xinc = 360. / grid->xsize;

                        if ( fabs(grid->xinc-xinc) > 0.0 )
                          {
                            grid->xinc = xinc;
                            if ( CDI_Debug ) Message("set xinc to %g", grid->xinc);
                          }
                      }
                  }
              }
            grid->xdef = 2;
          }
        grid->ydef  = 0;
        /* if ( IS_NOT_EQUAL(grid->yfirst, 0) || IS_NOT_EQUAL(grid->ylast, 0) ) */
          {
            if ( grid->ysize > 1 )
              {
                if ( editionNumber <= 1 )
                  {
                  }
              }
            grid->ydef = 2;
          }
        break;
      }
    case GRID_LCC:
      {
        int nlon, nlat;
        long lpar;

        FAIL_ON_GRIB_ERROR(grib_get_long, gh, "Nx", &lpar);
        nlon = (int)lpar;
        FAIL_ON_GRIB_ERROR(grib_get_long, gh, "Ny", &lpar);
        nlat = (int)lpar;

        if ( numberOfPoints != nlon*nlat )
          Error("numberOfPoints (%d) and gridSize (%d) differ!", (int)numberOfPoints, nlon*nlat);

        grid->size  = (int)numberOfPoints;
        grid->xsize = nlon;
        grid->ysize = nlat;

        FAIL_ON_GRIB_ERROR(grib_get_double, gh, "DxInMetres", &grid->lcc_xinc);
        FAIL_ON_GRIB_ERROR(grib_get_double, gh, "DyInMetres", &grid->lcc_yinc);
        FAIL_ON_GRIB_ERROR(grib_get_double, gh, "longitudeOfFirstGridPointInDegrees", &grid->lcc_originLon);
        FAIL_ON_GRIB_ERROR(grib_get_double, gh, "latitudeOfFirstGridPointInDegrees", &grid->lcc_originLat);
        FAIL_ON_GRIB_ERROR(grib_get_double, gh, "LoVInDegrees", &grid->lcc_lonParY);
        FAIL_ON_GRIB_ERROR(grib_get_double, gh, "Latin1InDegrees", &grid->lcc_lat1);
        FAIL_ON_GRIB_ERROR(grib_get_double, gh, "Latin2InDegrees", &grid->lcc_lat2);

        if ( editionNumber <= 1 )
          {
            FAIL_ON_GRIB_ERROR(grib_get_long, gh, "projectionCenterFlag", &lpar);
            grid->lcc_projflag  = (int) lpar;
            FAIL_ON_GRIB_ERROR(grib_get_long, gh, "scanningMode", &lpar);
            grid->lcc_scanflag  = (int) lpar;
          }

        grid->xdef   = 0;
        grid->ydef   = 0;

        break;
      }
    case GRID_SPECTRAL:
      {
        size_t len = 256;
        char typeOfPacking[256];
        FAIL_ON_GRIB_ERROR(grib_get_string, gh, "packingType", typeOfPacking, &len);
        grid->lcomplex = 0;
        if ( strncmp(typeOfPacking, "spectral_complex", len) == 0 ) grid->lcomplex = 1;

        /* FIXME: assert(datasize >= INT_MIN && datasize <= INT_MAX) */
        grid->size  = (int)datasize;
        long lpar;
        FAIL_ON_GRIB_ERROR(grib_get_long, gh, "J", &lpar);
        /* FIXME: assert(lpar >= INT_MIN && lpar <= INT_MAX) */
        grid->trunc = (int)lpar;

        break;
      }
    case GRID_GME:
      {
        /* FIXME: assert(numberOfPoints <= INT_MAX && numberOfPoints >= INT_MIN) */
        grid->size  = (int)numberOfPoints;
        long lpar;
        /* FIXME: assert(lpar >= INT_MIN && lpar <= INT_MAX) */
        if ( grib_get_long(gh, "nd", &lpar) == 0 ) grid->nd  = (int)lpar;
        /* FIXME: assert(lpar >= INT_MIN && lpar <= INT_MAX) */
        if ( grib_get_long(gh, "Ni", &lpar) == 0 ) grid->ni  = (int)lpar;
        /* FIXME: assert(lpar >= INT_MIN && lpar <= INT_MAX) */
        if ( grib_get_long(gh, "n2", &lpar) == 0 ) grid->ni2 = (int)lpar;
        /* FIXME: assert(lpar >= INT_MIN && lpar <= INT_MAX) */
        if ( grib_get_long(gh, "n3", &lpar) == 0 ) grid->ni3 = (int)lpar;

        break;
      }
    case GRID_UNSTRUCTURED:
      {
        unsigned char uuid[CDI_UUID_SIZE];
        /*
        char reference_link[8192];
        size_t len = sizeof(reference_link);
        reference_link[0] = 0;
        */

        /* FIXME: assert(numberOfPoints <= INT_MAX && numberOfPoints >= INT_MIN) */
            grid->size  = (int)numberOfPoints;
        long lpar;
        if ( grib_get_long(gh, "numberOfGridUsed", &lpar) == 0 )
          {
            /* FIXME: assert(lpar <= INT_MAX && lpar >= INT_MIN) */
            grid->number   = (int)lpar;
            /* FIXME: assert(lpar <= INT_MAX && lpar >= INT_MIN) */
            if ( grib_get_long(gh, "numberOfGridInReference", &lpar) == 0 )
              grid->position = (int)lpar;
            /*
            if ( grib_get_string(gh, "gridDescriptionFile", reference_link, &len) == 0 )
              {
                if ( strncmp(reference_link, "file://", 7) == 0 )
                  grid->reference = strdupx(reference_link);
              }
            */
            size_t len = (size_t)CDI_UUID_SIZE;
            if ( grib_get_bytes(gh, "uuidOfHGrid", uuid, &len) == 0)
              {
                memcpy(grid->uuid, uuid, CDI_UUID_SIZE);
              }
          }
        break;
      }
    case GRID_GENERIC:
      {
        int nlon = 0, nlat = 0;
        long lpar;
        /* FIXME: assert(lpar <= INT_MAX && lpar >= INT_MIN) */
        if ( grib_get_long(gh, "Ni", &lpar) == 0 ) nlon = (int)lpar;
        /* FIXME: assert(lpar <= INT_MAX && lpar >= INT_MIN) */
        if ( grib_get_long(gh, "Nj", &lpar) == 0 ) nlat = (int)lpar;

        /* FIXME: assert(numberOfPoints <= INT_MAX && numberOfPoints >= INT_MIN) */
        grid->size  = (int)numberOfPoints;

        if ( nlon > 0 && nlat > 0 && nlon*nlat == grid->size )
          {
            grid->xsize = nlon;
            grid->ysize = nlat;
          }
        else
          {
            grid->xsize = 0;
            grid->ysize = 0;
          }

        break;
      }
    default:
      {
        Error("Unsupported grid type: %s", gridNamePtr(gridtype));
        break;
      }
    }

  grid->isRotated = FALSE;
  if ( gribapiGetIsRotated(gh) )
    {
      grid->isRotated = TRUE;
      FAIL_ON_GRIB_ERROR(grib_get_double, gh, "latitudeOfSouthernPoleInDegrees",  &grid->ypole);
      FAIL_ON_GRIB_ERROR(grib_get_double, gh, "longitudeOfSouthernPoleInDegrees", &grid->xpole);
      FAIL_ON_GRIB_ERROR(grib_get_double, gh, "angleOfRotation", &grid->angle);
      /* change from south to north pole */
      grid->ypole = -grid->ypole;
      grid->xpole =  grid->xpole - 180;
    }

  grid->xvals = NULL;
  grid->yvals = NULL;
  grid->type  = gridtype;
}
#endif
