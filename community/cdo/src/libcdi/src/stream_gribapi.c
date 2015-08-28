#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <stdio.h>

#include "dmemory.h"
#include "cdi.h"
#include "cdi_int.h"
#include "file.h"
#include "gribapi_utilities.h"
#include "stream_grb.h"
#include "varscan.h"
#include "datetime.h"
#include "vlist.h"
#include "stream_grb.h"
#include "calendar.h"
#include "subtype.h"


#if  defined  (HAVE_LIBGRIB_API)
#  include "cgribex.h"      /* gribGetSize, gribRead, gribGetZip, GRIB1_LTYPE_99 */
#  include "gribapi.h"

#  include <grib_api.h>
#endif

extern int cdiInventoryMode;

#if  defined  (HAVE_LIBGRIB_API)
static const var_tile_t dummy_tiles = { -1, -1, -1, -1, -1, -1 };
#endif

typedef struct {
  int param;
  int level1;
  int level2;
  int ltype;
  int tsteptype;
  char name[32];

  var_tile_t tiles;

} compvar2_t;


#if  defined  (HAVE_LIBGRIB_API)
static
int my_grib_set_double(grib_handle* h, const char* key, double val)
{
  if ( cdiGribApiDebug )
    fprintf(stderr, "grib_set_double(\tgrib_handle* h, \"%s\", %f)\n", key, val);

  return grib_set_double(h, key, val);
}

static
int my_grib_set_long(grib_handle* h, const char* key, long val)
{
  if ( cdiGribApiDebug )
    fprintf(stderr, "grib_set_long(  \tgrib_handle* h, \"%s\", %ld)\n", key, val);

  return grib_set_long(h, key, val);
}

static
int my_grib_set_string(grib_handle* h, const char* key, const char* val, size_t* length)
{
  if ( cdiGribApiDebug )
    fprintf(stderr, "grib_set_string(\tgrib_handle* h, \"%s\", \"%s\")\n", key, val);

  return grib_set_string(h, key, val, length);
}

static
int gribapiGetZaxisType(long editionNumber, int grib_ltype)
{
  int zaxistype = ZAXIS_GENERIC;

  if ( editionNumber <= 1 )
    {
      zaxistype = grib1ltypeToZaxisType(grib_ltype);
    }
  else
    {
      zaxistype = grib2ltypeToZaxisType(grib_ltype);
    }

  return (zaxistype);
}

static
int getTimeunits(long unitsOfTime)
{
  int timeunits = -1;

  switch (unitsOfTime)
    {
    case 13:  timeunits = TUNIT_SECOND;  break;
    case  0:  timeunits = TUNIT_MINUTE;  break;
    case  1:  timeunits = TUNIT_HOUR;    break;
    case 10:  timeunits = TUNIT_3HOURS;  break;
    case 11:  timeunits = TUNIT_6HOURS;  break;
    case 12:  timeunits = TUNIT_12HOURS; break;
    case  2:  timeunits = TUNIT_DAY;     break;
    default:  timeunits = TUNIT_HOUR;    break;
    }

  return (timeunits);
}

static
double timeunit_factor(int tu1, int tu2)
{
  double factor = 1;

  if ( tu2 == TUNIT_HOUR )
    {
      switch (tu1)
        {
        case TUNIT_SECOND:  factor = 3600;   break;
        case TUNIT_MINUTE:  factor = 60;     break;
        case TUNIT_HOUR:    factor = 1;      break;
        case TUNIT_3HOURS:  factor = 1./3;   break;
        case TUNIT_6HOURS:  factor = 1./6;   break;
        case TUNIT_12HOURS: factor = 1./12;  break;
        case TUNIT_DAY:     factor = 1./24;  break;
        }
    }

  return (factor);
}

static
int gribapiGetTimeUnits(grib_handle *gh)
{
  int timeunits = -1;
  long unitsOfTime = -1;

  grib_get_long(gh, "indicatorOfUnitOfTimeRange", &unitsOfTime);

  GRIB_CHECK(my_grib_set_long(gh, "stepUnits", unitsOfTime), 0);

  timeunits = getTimeunits(unitsOfTime);

  return (timeunits);
}

static
int gribapiGetEndStep(grib_handle *gh, int startStep, int timeunits)
{
  int endStep = startStep;
  int timeunits2 = timeunits;

  long unitsOfTime;
  int status = grib_get_long(gh, "stepUnits", &unitsOfTime);
  if ( status == 0 ) timeunits2 = getTimeunits(unitsOfTime);
  //timeunits2 = gribapiGetTimeUnits(gh);

  long lpar;
  status = grib_get_long(gh, "endStep", &lpar);

  if ( status == 0 )
    endStep = (int) (((double)lpar * timeunit_factor(timeunits, timeunits2)) + 0.5);
  // printf("%d %d %d %d %d %g\n", startStep, endStep, lpar, timeunits, timeunits2, timeunit_factor(timeunits, timeunits2));

  return (endStep);
}

static
void gribapiGetDataDateTime(grib_handle *gh, int *datadate, int *datatime)
{
  long lpar;

  GRIB_CHECK(grib_get_long(gh, "dataDate", &lpar), 0);
  *datadate = (int) lpar;
  GRIB_CHECK(grib_get_long(gh, "dataTime", &lpar), 0);  //FIXME: This looses the seconds in GRIB2 files.
  *datatime = (int) lpar*100;
}

static
void gribapiSetDataDateTime(grib_handle *gh, int datadate, int datatime)
{
  GRIB_CHECK(my_grib_set_long(gh, "dataDate", datadate), 0);
  GRIB_CHECK(my_grib_set_long(gh, "dataTime", datatime/100), 0);
}

static
int gribapiGetValidityDateTime(grib_handle *gh, int *vdate, int *vtime)
{
  int rdate, rtime;
  int timeUnits, startStep = 0, endStep;
  int tstepRange = 0;
  int range;
  int status;
  long lpar;
  long sigofrtime = 3;

  if ( gribEditionNumber(gh) > 1 )
    {
      GRIB_CHECK(grib_get_long(gh, "significanceOfReferenceTime", &sigofrtime), 0);
    }
  else
    {
      GRIB_CHECK(grib_get_long(gh, "timeRangeIndicator", &sigofrtime), 0);
    }

  if ( sigofrtime == 3 )        //XXX: This looks like a bug to me, because timeRangeIndicator == 3 does not seem to have the same meaning as significanceOfReferenceTime == 3. I would recommend replacing this condition with `if(!gribapiTimeIsFC())`.
    {
      gribapiGetDataDateTime(gh, vdate, vtime);
    }
  else
    {
      gribapiGetDataDateTime(gh, &rdate, &rtime);

      status = grib_get_long(gh, "forecastTime", &lpar);
      if ( status == 0 ) startStep = (int) lpar;
      timeUnits = gribapiGetTimeUnits(gh);
      endStep = gribapiGetEndStep(gh, startStep, timeUnits);

      range = endStep - startStep;

      if ( range > 0 )
	{
	  if ( startStep == 0 ) tstepRange = -1;
	  else                  tstepRange =  1;
	}

      {
	static int lprint = TRUE;
	extern int grib_calendar;
	int ryear, rmonth, rday, rhour, rminute, rsecond;
	int julday, secofday;
	int64_t time_period = endStep;
        int64_t addsec;

	cdiDecodeDate(rdate, &ryear, &rmonth, &rday);
	cdiDecodeTime(rtime, &rhour, &rminute, &rsecond);

        if ( rday > 0 )
          {
            encode_caldaysec(grib_calendar, ryear, rmonth, rday, rhour, rminute, rsecond, &julday, &secofday);

            addsec = 0;
            switch ( timeUnits )
              {
              case TUNIT_SECOND:  addsec =         time_period; break;
              case TUNIT_MINUTE:  addsec =    60 * time_period; break;
              case TUNIT_HOUR:    addsec =  3600 * time_period; break;
              case TUNIT_3HOURS:  addsec = 10800 * time_period; break;
              case TUNIT_6HOURS:  addsec = 21600 * time_period; break;
              case TUNIT_12HOURS: addsec = 43200 * time_period; break;
              case TUNIT_DAY:     addsec = 86400 * time_period; break;
              default:
                if ( lprint )
                  {
                    Warning("Time unit %d unsupported", timeUnits);
                    lprint = FALSE;
                  }
                break;
              }

            julday_add_seconds(addsec, &julday, &secofday);

            decode_caldaysec(grib_calendar, julday, secofday, &ryear, &rmonth, &rday, &rhour, &rminute, &rsecond);
          }

	*vdate = cdiEncodeDate(ryear, rmonth, rday);
	*vtime = cdiEncodeTime(rhour, rminute, rsecond);
      }
    }

  return (tstepRange);
}

static
void grib1GetLevel(grib_handle *gh, int *leveltype, int *lbounds, int *level1, int *level2)
{
  *leveltype = 0;
  *lbounds   = 0;
  *level1    = 0;
  *level2    = 0;

  long lpar;
  if(!grib_get_long(gh, "indicatorOfTypeOfLevel", &lpar))       //1 byte
    {
      *leveltype = (int) lpar;

      switch (*leveltype)
	{
	case GRIB1_LTYPE_SIGMA_LAYER:
	case GRIB1_LTYPE_HYBRID_LAYER:
	case GRIB1_LTYPE_LANDDEPTH_LAYER:
	  { *lbounds = 1; break; }
	}

      if ( *lbounds )
	{
	  GRIB_CHECK(grib_get_long(gh, "topLevel", &lpar), 0);  //1 byte
	  *level1 = (int)lpar;
	  GRIB_CHECK(grib_get_long(gh, "bottomLevel", &lpar), 0);       //1 byte
	  *level2 = (int)lpar;
	}
      else
	{
          double dlevel;
	  GRIB_CHECK(grib_get_double(gh, "level", &dlevel), 0); //2 byte
	  if ( *leveltype == 100 ) dlevel *= 100;
	  if ( dlevel < -2.e9 || dlevel > 2.e9 ) dlevel = 0;
	  if ( *leveltype == GRIB1_LTYPE_99 ) *leveltype = 100;

	  *level1 = (int) dlevel;
	  *level2 = 0;
	}
    }
}

static
double grib2ScaleFactor(long factor)
{
  switch(factor)
    {
      case GRIB_MISSING_LONG: return 1;
      case 0: return 1;
      case 1: return 0.1;
      case 2: return 0.01;
      case 3: return 0.001;
      case 4: return 0.0001;
      case 5: return 0.00001;
      case 6: return 0.000001;
      case 7: return 0.0000001;
      case 8: return 0.00000001;
      case 9: return 0.000000001;
      default: return 0;
    }
}

static
int calcLevel(int level_sf, long factor, long level)
{
  double result = 0;
  if(level != GRIB_MISSING_LONG) result = (double)level*grib2ScaleFactor(factor);
  if(level_sf) result *= level_sf;
  return (int)result;
}

static
void grib2GetLevel(grib_handle *gh, int *leveltype1, int *leveltype2, int *lbounds, int *level1, 
                   int *level2, int *level_sf, int *level_unit)
{
  int status;
  long lpar;
  long factor;

  *leveltype1 = 0;
  *leveltype2 = -1;
  *lbounds    = 0;
  *level1     = 0;
  *level2     = 0;
  *level_sf   = 0;
  *level_unit = 0;

  status = grib_get_long(gh, "typeOfFirstFixedSurface", &lpar); //1 byte
  if ( status == 0 )
    {
      long llevel;

      *leveltype1 = (int) lpar;

      status = grib_get_long(gh, "typeOfSecondFixedSurface", &lpar); //1 byte
      /* FIXME: assert(lpar >= INT_MIN && lpar <= INT_MAX) */
      if ( status == 0 ) *leveltype2 = (int)lpar;

      if ( *leveltype1 != 255 && *leveltype2 != 255 && *leveltype2 > 0 ) *lbounds = 1;
      switch(*leveltype1)
        {
          case GRIB2_LTYPE_REFERENCE:
            if(*leveltype2 == 1) *lbounds = 0;
            break;

          case GRIB2_LTYPE_LANDDEPTH:
            *level_sf = 1000;
            *level_unit = CDI_UNIT_M;
            break;

          case GRIB2_LTYPE_ISOBARIC:
            *level_sf = 1000;
            *level_unit = CDI_UNIT_PA;
            break;

          case GRIB2_LTYPE_SIGMA:
            *level_sf = 1000;
            *level_unit = 0;
            break;
        }

      GRIB_CHECK(grib_get_long(gh, "scaleFactorOfFirstFixedSurface", &factor), 0);      //1 byte
      GRIB_CHECK(grib_get_long(gh, "scaledValueOfFirstFixedSurface", &llevel), 0);      //4 byte
      *level1 = calcLevel(*level_sf, factor, llevel);

      if ( *lbounds )
        {
          GRIB_CHECK(grib_get_long(gh, "scaleFactorOfSecondFixedSurface", &factor), 0); //1 byte
          GRIB_CHECK(grib_get_long(gh, "scaledValueOfSecondFixedSurface", &llevel), 0); //4 byte
          *level2 = calcLevel(*level_sf, factor, llevel);
        }
    }
}

static
void gribGetLevel(grib_handle *gh, int* leveltype1, int* leveltype2, int* lbounds, int* level1, int* level2, int* level_sf, int* level_unit, var_tile_t* tiles)
{
  if ( gribEditionNumber(gh) <= 1 )
    {
      grib1GetLevel(gh, leveltype1, lbounds, level1, level2);
      *leveltype2 = -1;
      *level_sf = 0;
      *level_unit = 0;
    }
  else
    {
      grib2GetLevel(gh, leveltype1, leveltype2, lbounds, level1, level2, level_sf, level_unit);

      /* read in tiles attributes (if there are any) */
      tiles->tileindex = (int)gribGetLongDefault(gh, cdiSubtypeAttributeName[SUBTYPE_ATT_TILEINDEX], -1);
      tiles->totalno_of_tileattr_pairs = (int)gribGetLongDefault(gh, cdiSubtypeAttributeName[SUBTYPE_ATT_TOTALNO_OF_TILEATTR_PAIRS], -1);
      tiles->tileClassification = (int)gribGetLongDefault(gh, cdiSubtypeAttributeName[SUBTYPE_ATT_TILE_CLASSIFICATION], -1);
      tiles->numberOfTiles = (int)gribGetLongDefault(gh, cdiSubtypeAttributeName[SUBTYPE_ATT_NUMBER_OF_TILES], -1);
      tiles->numberOfAttributes = (int)gribGetLongDefault(gh, cdiSubtypeAttributeName[SUBTYPE_ATT_NUMBER_OF_ATTR], -1);
      tiles->attribute = (int)gribGetLongDefault(gh, cdiSubtypeAttributeName[SUBTYPE_ATT_TILEATTRIBUTE], -1);
    }
}

static
void gribapiGetString(grib_handle *gh, const char *key, char *string, size_t length)
{
  string[0] = 0;

  int ret = grib_get_string(gh, key, string, &length);
  if (ret != 0)
    {
      fprintf(stderr, "grib_get_string(gh, \"%s\", ...) failed!\n", key);
      GRIB_CHECK(ret, 0);
    }
  if      ( length == 8 && memcmp(string, "unknown", length) == 0 ) string[0] = 0;
  else if ( length == 2 && memcmp(string, "~", length)       == 0 ) string[0] = 0;
}

#if  defined  (HAVE_LIBGRIB_API)
static
void gribapiAddRecord(stream_t * streamptr, int param, grib_handle *gh,
                      size_t recsize, off_t position, int datatype, int comptype, const char *varname,
                      int leveltype1, int leveltype2, int lbounds, int level1, int level2, int level_sf, int level_unit,
                      const var_tile_t *tiles, int lread_additional_keys)
{
  int zaxistype;
  int gridID = CDI_UNDEFID, varID;
  int levelID = 0;
  int tsID, recID;
  int numavg;
  int tsteptype;
  record_t *record;
  grid_t grid;
  int vlistID;
  long lpar;
  int status;
  char stdname[CDI_MAX_NAME], longname[CDI_MAX_NAME], units[CDI_MAX_NAME];
  size_t vlen;
  long ens_index = 0, ens_count = 0, ens_forecast_type = 0;

  vlistID = streamptr->vlistID;
  tsID    = streamptr->curTsID;
  recID   = recordNewEntry(streamptr, tsID);
  record  = &streamptr->tsteps[tsID].records[recID];

  tsteptype = gribapiGetTsteptype(gh);
  // numavg  = ISEC1_AvgNum;
  numavg  = 0;

  // fprintf(stderr, "param %d %d %d %d\n", param, level1, level2, leveltype1);

  (*record).size      = recsize;
  (*record).position  = position;
  (*record).param     = param;
  (*record).ilevel    = level1;
  (*record).ilevel2   = level2;
  (*record).ltype     = leveltype1;
  (*record).tsteptype = tsteptype;
  if ( tiles ) (*record).tiles = *tiles;
  else         (*record).tiles = dummy_tiles;

  //FIXME: This may leave the variable name unterminated (which is the behavior that I found in the code).
  //       I don't know precisely how this field is used, so I did not change this behavior to avoid regressions,
  //       but I think that it would be better to at least add a line
  //
  //           record->varname[sizeof(record->varname) - 1] = 0;`
  //
  //       after the `strncpy()` call.
  //
  //       I would consider using strdup() (that requires POSIX-2008 compliance, though), or a similar homebrew approach.
  //       I. e. kick the fixed size array and allocate enough space, whatever that may be.
  strncpy(record->varname, varname, sizeof(record->varname));

  gribapiGetGrid(gh, &grid);

  gridID = varDefGrid(vlistID, &grid, 0);

  zaxistype = gribapiGetZaxisType(gribEditionNumber(gh), leveltype1);

  switch (zaxistype)
    {
    case ZAXIS_HYBRID:
    case ZAXIS_HYBRID_HALF:
      {
        size_t vctsize;
        size_t dummy;
        double *vctptr;

        GRIB_CHECK(grib_get_long(gh, "NV", &lpar), 0);
        /* FIXME: assert(lpar >= 0) */
        vctsize = (size_t)lpar;
        if ( vctsize > 0 )
          {
            vctptr = (double *) malloc(vctsize*sizeof(double));
            dummy = vctsize;
            GRIB_CHECK(grib_get_double_array(gh, "pv", vctptr, &dummy), 0);
            varDefVCT(vctsize, vctptr);
            free(vctptr);
          }
        break;
      }
    case ZAXIS_REFERENCE:
      {
        size_t len;
        unsigned char uuid[CDI_UUID_SIZE];
        long ltmp;
        long nhlev, nvgrid;

        GRIB_CHECK(grib_get_long(gh, "NV", &lpar), 0);
        if ( lpar != 6 )
          {
            fprintf(stderr, "Warning ...\n");
          }
        GRIB_CHECK(grib_get_long(gh, "nlev", &ltmp), 0);
        nhlev = ltmp;
        GRIB_CHECK(grib_get_long(gh, "numberOfVGridUsed", &ltmp), 0);
        nvgrid = ltmp;
        len = (size_t)CDI_UUID_SIZE;
        memset(uuid, 0, CDI_UUID_SIZE);
        GRIB_CHECK(grib_get_bytes(gh, "uuidOfVGrid", uuid, &len), 0);
        varDefZAxisReference((int) nhlev, (int) nvgrid, uuid);
        break;
      }
    }

  // if ( datatype > 32 ) datatype = DATATYPE_PACK32;
  if ( datatype <  0 ) datatype = DATATYPE_PACK;

  stdname[0] = 0;
  longname[0] = 0;
  units[0] = 0;

  if ( varname[0] != 0 )
    {
      vlen = CDI_MAX_NAME;
      gribapiGetString(gh, "name", longname, vlen);
      vlen = CDI_MAX_NAME;
      gribapiGetString(gh, "units", units, vlen);

      {
        vlen = CDI_MAX_NAME;
        status = grib_get_string(gh, "cfName", stdname, &vlen);
        if ( status != 0 || vlen <= 1 ) stdname[0] = 0;
        else if ( strncmp(stdname, "unknown", 7) == 0 ) stdname[0] = 0;
      }
    }
  // fprintf(stderr, "param %d name %s %s %s\n", param, name, longname, units);

  /* add the previously read record data to the (intermediate) list of records */
  int tile_index = -1;
  varAddRecord(recID, param, gridID, zaxistype, lbounds, level1, level2, level_sf, level_unit,
	       datatype, &varID, &levelID, tsteptype, numavg, leveltype1, leveltype2,
	       varname, stdname, longname, units, tiles, &tile_index);

  (*record).varID   = (short)varID;
  (*record).levelID = (short)levelID;

  varDefCompType(varID, comptype);

  /*
    Get the ensemble Info from the grib-2 Tables and update the intermediate datastructure.
    Further update to the "vlist" is handled in the same way as for GRIB-1 by "cdi_generate_vars"
  */
  status = grib_get_long(gh, "typeOfEnsembleForecast", &ens_forecast_type );
  if ( status == 0 )
    {
      GRIB_CHECK(grib_get_long(gh, "numberOfForecastsInEnsemble", &ens_count ), 0);
      GRIB_CHECK(grib_get_long(gh, "perturbationNumber", &ens_index ), 0);
    }

  if ( ens_index > 0 )
    varDefEnsembleInfo(varID, (int)ens_index, (int)ens_count, (int)ens_forecast_type);

  long typeOfGeneratingProcess = 0;
  status = grib_get_long(gh, "typeOfGeneratingProcess", &typeOfGeneratingProcess);
  if ( status == 0 )
    varDefTypeOfGeneratingProcess(varID, (int) typeOfGeneratingProcess);

  long productDefinitionTemplate = 0;
  status = grib_get_long(gh, "productDefinitionTemplateNumber", &productDefinitionTemplate);
  if ( status == 0 )
    varDefProductDefinitionTemplate(varID, (int) productDefinitionTemplate);

  int    i;
  long   lval;
  double dval;

  if (lread_additional_keys)
    for ( i = 0; i < cdiNAdditionalGRIBKeys; i++ )
      {
        /* note: if the key is not defined, we do not throw an error! */
        if ( grib_get_long(gh, cdiAdditionalGRIBKeys[i], &lval) == 0 )
          varDefOptGribInt(varID, tile_index, lval, cdiAdditionalGRIBKeys[i]);
        if ( grib_get_double(gh, cdiAdditionalGRIBKeys[i], &dval) == 0 )
          varDefOptGribDbl(varID, tile_index, dval, cdiAdditionalGRIBKeys[i]);
      }

  if ( varInqInst(varID) == CDI_UNDEFID )
    {
      long center, subcenter;
      int instID;
      GRIB_CHECK(grib_get_long(gh, "centre", &center), 0);
      GRIB_CHECK(grib_get_long(gh, "subCentre", &subcenter), 0);
      instID    = institutInq((int)center, (int)subcenter, NULL, NULL);
      if ( instID == CDI_UNDEFID )
	instID = institutDef((int)center, (int)subcenter, NULL, NULL);
      varDefInst(varID, instID);
    }

  if ( varInqModel(varID) == CDI_UNDEFID )
    {
      int modelID;
      long processID;
      status = grib_get_long(gh, "generatingProcessIdentifier", &processID);
      if ( status == 0 )
	{
          /* FIXME: assert(processID >= INT_MIN && processID <= INT_MAX) */
	  modelID = modelInq(varInqInst(varID), (int)processID, NULL);
	  if ( modelID == CDI_UNDEFID )
	    modelID = modelDef(varInqInst(varID), (int)processID, NULL);
	  varDefModel(varID, modelID);
	}
    }

  if ( varInqTable(varID) == CDI_UNDEFID )
    {
      int pdis, pcat, pnum;

      cdiDecodeParam(param, &pnum, &pcat, &pdis);

      if ( pdis == 255 )
	{
	  int tableID;
	  int tabnum = pcat;

	  tableID = tableInq(varInqModel(varID), tabnum, NULL);

	  if ( tableID == CDI_UNDEFID )
	    tableID = tableDef(varInqModel(varID), tabnum, NULL);
	  varDefTable(varID, tableID);
	}
    }

  streamptr->tsteps[tsID].nallrecs++;
  streamptr->nrecs++;

  if ( CDI_Debug )
    Message("varID = %d  param = %d  zaxistype = %d  gridID = %d  levelID = %d",
	    varID, param, zaxistype, gridID, levelID);
}
#endif

static compvar2_t gribapiVarSet(int param, int level1, int level2, int leveltype, 
                                int tsteptype, char *name, var_tile_t tiles_data)
{
  compvar2_t compVar;
  size_t maxlen = sizeof(compVar.name);
  size_t len = strlen(name);
  if ( len > maxlen ) len = maxlen;

  compVar.param     = param;
  compVar.level1    = level1;
  compVar.level2    = level2;
  compVar.ltype     = leveltype;
  compVar.tsteptype = tsteptype;
  memset(compVar.name, 0, maxlen);
  memcpy(compVar.name, name, len);

  compVar.tiles = tiles_data;
  return (compVar);
}
#endif

#ifdef HAVE_LIBGRIB_API
static
int gribapiVarCompare(compvar2_t compVar, record_t record, int flag)
{
  compvar2_t compVar0;
  size_t maxlen = sizeof(compVar.name);

  compVar0.param     = record.param;
  compVar0.level1    = record.ilevel;
  compVar0.level2    = record.ilevel2;
  compVar0.ltype     = record.ltype;
  compVar0.tsteptype = record.tsteptype;
  memcpy(compVar0.name, record.varname, maxlen);

  if ( flag == 0 )
    {
      if ( compVar0.tsteptype == TSTEP_INSTANT  && compVar.tsteptype == TSTEP_INSTANT3 ) compVar0.tsteptype = TSTEP_INSTANT3;
      if ( compVar0.tsteptype == TSTEP_INSTANT3 && compVar.tsteptype == TSTEP_INSTANT  ) compVar0.tsteptype = TSTEP_INSTANT;
    }

  compVar0.tiles = record.tiles;

  int rstatus = memcmp(&compVar0, &compVar, sizeof(compvar2_t));

  return (rstatus);
}

static void ensureBufferSize(size_t requiredSize, size_t* curSize, unsigned char** buffer) {
  if ( *curSize < requiredSize )
    {
      *curSize = requiredSize;
      *buffer = realloc(*buffer, *curSize);
    }
}

static
grib_handle* gribapiGetDiskRepresentation(size_t recsize, size_t* buffersize, unsigned char** gribbuffer, int* outDatatype, int* outCompressionType, long* outUnzipsize)
{
  int lieee = FALSE;

  grib_handle* gh = grib_handle_new_from_message(NULL, (void *) *gribbuffer, recsize);
  if(gribEditionNumber(gh) > 1)
    {
      size_t len = 256;
      char typeOfPacking[256];

      if ( grib_get_string(gh, "packingType", typeOfPacking, &len) == 0 )
        {
          // fprintf(stderr, "packingType %d %s\n", len, typeOfPacking);
          if      ( strncmp(typeOfPacking, "grid_jpeg", len) == 0 ) *outCompressionType = COMPRESS_JPEG;
          else if ( strncmp(typeOfPacking, "grid_ccsds", len) == 0 ) *outCompressionType = COMPRESS_SZIP;
          else if ( strncmp(typeOfPacking, "grid_ieee", len) == 0 ) lieee = TRUE;
        }
    }
  else
    {
      if( gribGetZip((long)recsize, *gribbuffer, outUnzipsize) > 0 )
        {
          *outCompressionType = COMPRESS_SZIP;
          ensureBufferSize((size_t)*outUnzipsize + 100, buffersize, gribbuffer);
        }
      else
        {
          *outCompressionType = COMPRESS_NONE;
        }
    }

  if ( lieee )
    {
      *outDatatype = DATATYPE_FLT64;
      long precision;
      int status = grib_get_long(gh, "precision", &precision);
      if ( status == 0 && precision == 1 ) *outDatatype = DATATYPE_FLT32;
    }
  else
    {
      *outDatatype = DATATYPE_PACK;
      long bitsPerValue;
      if ( grib_get_long(gh, "bitsPerValue", &bitsPerValue) == 0 )
        {
          if ( bitsPerValue > 0 && bitsPerValue <= 32 ) *outDatatype = (int)bitsPerValue;
        }
    }
  return gh;
}
#endif

#if  defined  (HAVE_LIBGRIB_API)
typedef enum { CHECKTIME_OK, CHECKTIME_SKIP, CHECKTIME_STOP, CHECKTIME_INCONSISTENT } checkTimeResult;
static checkTimeResult checkTime(stream_t* streamptr, compvar2_t compVar, const DateTime* verificationTime, const DateTime* expectedVTime) {
  //First determine whether the current record exists already.
  int recID = 0;
  for ( ; recID < streamptr->nrecs; recID++ )
    {
      if ( gribapiVarCompare(compVar, streamptr->tsteps[0].records[recID], 1) == 0 ) break;
    }
  int recordExists = recID < streamptr->nrecs;

  //Then we need to know whether the verification time is consistent.
  int consistentTime = !memcmp(verificationTime, expectedVTime, sizeof(*verificationTime));

  //Finally, we make a decision.
  if ( cdiInventoryMode == 1 )
    {
      if ( recordExists ) return CHECKTIME_STOP;
      if ( !consistentTime ) return CHECKTIME_INCONSISTENT;
    }
  else
    {
      if ( !consistentTime ) return CHECKTIME_STOP;
      if ( recordExists ) return CHECKTIME_SKIP;
    }
  return CHECKTIME_OK;
}
#endif

#define gribWarning(text, nrecs, timestep, varname, param, level1, level2) do \
  { \
    char paramstr[32]; \
    cdiParamToString(param, paramstr, sizeof(paramstr)); \
    Warning("Record %2d (name=%s id=%s lev1=%d lev2=%d) timestep %d: %s", nrecs, varname, paramstr, level1, level2, timestep, text); \
  } \
while(0)

#if  defined  (HAVE_LIBGRIB_API)
int gribapiScanTimestep1(stream_t * streamptr)
{
  off_t recpos = 0;
  unsigned char *gribbuffer = NULL;
  size_t buffersize = 0;
  DateTime datetime0 = { .date = 10101, .time = 0 };
  int nrecs_scanned = 0;        //Only used for debug output.
  int warn_time = TRUE;
  // int warn_numavg = TRUE;
  int rdate = 0, rtime = 0, tunit = 0, fcast = 0;
  grib_handle *gh = NULL;

  streamptr->curTsID = 0;

  int tsID  = tstepsNewEntry(streamptr);
  taxis_t *taxis = &streamptr->tsteps[tsID].taxis;

  if ( tsID != 0 )
    Error("Internal problem! tstepsNewEntry returns %d", tsID);

  int fileID = streamptr->fileID;

  unsigned nrecs = 0;
  while ( TRUE )
    {
      int level1 = 0, level2 = 0;
      size_t recsize = (size_t)gribGetSize(fileID);
      recpos  = fileGetPos(fileID);

      if ( recsize == 0 )
        {
          streamptr->ntsteps = 1;
          break;
        }
      ensureBufferSize(recsize, &buffersize, &gribbuffer);

      size_t readsize = recsize;
      int rstatus = gribRead(fileID, gribbuffer, &readsize);        //Search for next 'GRIB', read the following record, and position file offset after it.
      if ( rstatus ) break;

      int datatype, comptype = 0;
      long unzipsize;
      gh = gribapiGetDiskRepresentation(recsize, &buffersize, &gribbuffer, &datatype, &comptype, &unzipsize);

      nrecs_scanned++;
      GRIB_CHECK(my_grib_set_double(gh, "missingValue", cdiDefaultMissval), 0);

      int param = gribapiGetParam(gh);
      int leveltype1 = -1, leveltype2 = -1, lbounds, level_sf, level_unit;
      var_tile_t tiles = dummy_tiles;
      gribGetLevel(gh, &leveltype1, &leveltype2, &lbounds, &level1, &level2, &level_sf, &level_unit, &tiles);

      char varname[256];
      varname[0] = 0;
      gribapiGetString(gh, "shortName", varname, sizeof(varname));

      int tsteptype = gribapiGetTsteptype(gh);
      int vdate = 0, vtime = 0;
      gribapiGetValidityDateTime(gh, &vdate, &vtime);
      DateTime datetime = { .date = vdate, .time = vtime };
      /*
      printf("%d %d %d\n", vdate, vtime, leveltype1);
      */

      if( datetime0.date == 10101 && datetime0.time == 0 )
        {
          if( memcmp(&datetime, &datetime0, sizeof(datetime)) || !nrecs )       //Do we really need this condition? I have included it in order not to change the number of times gribapiGetDataDateTime() etc. get called. But if those are sideeffect-free, this condition should be removed.
            {
              datetime0 = datetime;

              gribapiGetDataDateTime(gh, &rdate, &rtime);

              fcast = gribapiTimeIsFC(gh);
              if ( fcast ) tunit = gribapiGetTimeUnits(gh);
            }
        }

      if(nrecs)
        {
          checkTimeResult result = checkTime(streamptr, gribapiVarSet(param, level1, level2, leveltype1, tsteptype, varname, tiles), &datetime, &datetime0);
          if(result == CHECKTIME_STOP)
            {
              break;
            }
          else if(result == CHECKTIME_SKIP)
            {
              gribWarning("Parameter already exist, skipped!", nrecs_scanned, tsID+1, varname, param, level1, level2);
              continue;
            }
          else if(result == CHECKTIME_INCONSISTENT && warn_time)
            {
              gribWarning("Inconsistent verification time!", nrecs_scanned, tsID+1, varname, param, level1, level2);
              warn_time = FALSE;
            }
          assert(result == CHECKTIME_OK || result == CHECKTIME_INCONSISTENT);
        }
      /*
      if ( ISEC1_AvgNum )
        {
          if (  taxis->numavg && warn_numavg && (taxis->numavg != ISEC1_AvgNum) )
            {
              Message("Change numavg from %d to %d not allowed!",
                      taxis->numavg, ISEC1_AvgNum);
              warn_numavg = FALSE;
            }
          else
            {
              taxis->numavg = ISEC1_AvgNum;
            }
        }
      */
      nrecs++;

      if ( CDI_Debug )
        {
          char paramstr[32];
          cdiParamToString(param, paramstr, sizeof(paramstr));
          Message("%4u %8d name=%s id=%s ltype=%d lev1=%d lev2=%d vdate=%d vtime=%d",
                nrecs, (int)recpos, varname, paramstr, leveltype1, level1, level2, vdate, vtime);
        }

      var_tile_t *ptiles = NULL;
      if ( memcmp(&tiles, &dummy_tiles, sizeof(var_tile_t)) != 0 ) ptiles = &tiles;
      gribapiAddRecord(streamptr, param, gh, recsize, recpos, datatype, comptype, varname,
                       leveltype1, leveltype2, lbounds, level1, level2, level_sf, level_unit, ptiles, 1);

      grib_handle_delete(gh);
      gh = NULL;
    }

  if ( gh ) grib_handle_delete(gh);

  streamptr->rtsteps = 1;

  if ( nrecs == 0 ) return (CDI_EUFSTRUCT);

  cdi_generate_vars(streamptr);

  int taxisID = -1;
  if ( fcast )
    {
      taxisID = taxisCreate(TAXIS_RELATIVE);
      taxis->type  = TAXIS_RELATIVE;
      taxis->rdate = rdate;
      taxis->rtime = rtime;
      taxis->unit  = tunit;
    }
  else
    {
      taxisID = taxisCreate(TAXIS_ABSOLUTE);
      taxis->type  = TAXIS_ABSOLUTE;
    }

  taxis->vdate = (int)datetime0.date;
  taxis->vtime = (int)datetime0.time;

  int vlistID = streamptr->vlistID;
  vlistDefTaxis(vlistID, taxisID);

  int nrecords = streamptr->tsteps[0].nallrecs;
  if ( nrecords < streamptr->tsteps[0].recordSize )
    {
      streamptr->tsteps[0].recordSize = nrecords;
      streamptr->tsteps[0].records =
        (record_t *)xrealloc(streamptr->tsteps[0].records, (size_t)nrecords*sizeof(record_t));
    }

  streamptr->tsteps[0].recIDs = (int *)xmalloc((size_t)nrecords*sizeof(int));
  streamptr->tsteps[0].nrecs = nrecords;
  for ( int recID = 0; recID < nrecords; recID++ )
    streamptr->tsteps[0].recIDs[recID] = recID;

  streamptr->record->buffer     = gribbuffer;
  streamptr->record->buffersize = buffersize;

  if ( streamptr->ntsteps == -1 )
    {
      tsID = tstepsNewEntry(streamptr);
      if ( tsID != streamptr->rtsteps )
        Error("Internal error. tsID = %d", tsID);

      streamptr->tsteps[tsID-1].next   = TRUE;
      streamptr->tsteps[tsID].position = recpos;
    }

  if ( streamptr->ntsteps == 1 )
    {
      if ( taxis->vdate == 0 && taxis->vtime == 0 )
        {
          streamptr->ntsteps = 0;
          for ( int varID = 0; varID < streamptr->nvars; varID++ )
            {
              vlistDefVarTsteptype(vlistID, varID, TSTEP_CONSTANT);
            }
        }
    }

  return (0);
}
#endif


#ifdef HAVE_LIBGRIB_API
int gribapiScanTimestep2(stream_t * streamptr)
{
  int rstatus = 0;
  off_t recpos = 0;
  unsigned char *gribbuffer = NULL;
  size_t buffersize = 0;
  int fileID;
  DateTime datetime0;
  // int gridID;
  int recID;
  //  int warn_numavg = TRUE;
  grib_handle *gh = NULL;

  streamptr->curTsID = 1;

  fileID  = streamptr->fileID;
  int vlistID = streamptr->vlistID;
  int taxisID = vlistInqTaxis(vlistID);

  gribbuffer = (unsigned char *) streamptr->record->buffer;
  buffersize = streamptr->record->buffersize;

  int tsID = streamptr->rtsteps;
  if ( tsID != 1 )
    Error("Internal problem! unexpected timestep %d", tsID+1);

  taxis_t *taxis = &streamptr->tsteps[tsID].taxis;

  fileSetPos(fileID, streamptr->tsteps[tsID].position, SEEK_SET);

  cdi_create_records(streamptr, tsID);

  int nrecords = streamptr->tsteps[tsID].nallrecs;
  streamptr->tsteps[1].recIDs = (int *) malloc((size_t)nrecords*sizeof(int));
  streamptr->tsteps[1].nrecs = 0;
  for ( recID = 0; recID < nrecords; recID++ )
    streamptr->tsteps[1].recIDs[recID] = -1;

  for ( recID = 0; recID < nrecords; recID++ )
    {
      streamptr->tsteps[tsID].records[recID].position = streamptr->tsteps[0].records[recID].position;
      streamptr->tsteps[tsID].records[recID].size     = streamptr->tsteps[0].records[recID].size;
    }

  int nrecs_scanned = nrecords; //Only used for debug output
  int rindex = 0;
  while ( TRUE )
    {
      if ( rindex > nrecords ) break;

      size_t recsize = (size_t)gribGetSize(fileID);
      recpos  = fileGetPos(fileID);
      if ( recsize == 0 )
	{
	  streamptr->ntsteps = 2;
	  break;
	}
      ensureBufferSize(recsize, &buffersize, &gribbuffer);

      size_t readsize = recsize;
      rstatus = gribRead(fileID, gribbuffer, &readsize);
      if ( rstatus ) break;

      long unzipsize;
      if ( gribGetZip((long)recsize, gribbuffer, &unzipsize) > 0 )
        ensureBufferSize((size_t)unzipsize + 100, &buffersize, &gribbuffer);

      nrecs_scanned++;
      gh = grib_handle_new_from_message(NULL, (void *) gribbuffer, recsize);
      GRIB_CHECK(my_grib_set_double(gh, "missingValue", cdiDefaultMissval), 0);

      int param = gribapiGetParam(gh);
      int level1 = 0, level2 = 0, leveltype1, leveltype2, lbounds, level_sf, level_unit;
      var_tile_t tiles = dummy_tiles;
      gribGetLevel(gh, &leveltype1, &leveltype2, &lbounds, &level1, &level2, &level_sf, &level_unit, &tiles);

      char varname[256];
      varname[0] = 0;
      gribapiGetString(gh, "shortName", varname, sizeof(varname));

      int vdate = 0, vtime = 0;
      gribapiGetValidityDateTime(gh, &vdate, &vtime);

      if ( rindex == 0 )
	{
	  if ( taxisInqType(taxisID) == TAXIS_RELATIVE )
	    {
	      taxis->type  = TAXIS_RELATIVE;

              gribapiGetDataDateTime(gh, &(taxis->rdate), &(taxis->rtime));

	      taxis->unit  = gribapiGetTimeUnits(gh);
	    }
	  else
	    {
	      taxis->type  = TAXIS_ABSOLUTE;
	    }
	  taxis->vdate = vdate;
	  taxis->vtime = vtime;

	  datetime0.date = vdate;
	  datetime0.time = vtime;
	}

      int tsteptype = gribapiGetTsteptype(gh);
      /*
      if ( ISEC1_AvgNum )
	{
	  if (  taxis->numavg && warn_numavg &&
		(taxis->numavg != ISEC1_AvgNum) )
	    {
	      warn_numavg = FALSE;
	    }
	  else
	    {
	      taxis->numavg = ISEC1_AvgNum;
	    }
	}
      */
      DateTime datetime = {
        .date = vdate,
        .time = vtime
      };

      compvar2_t compVar = gribapiVarSet(param, level1, level2, leveltype1, tsteptype, varname, tiles);

      for ( recID = 0; recID < nrecords; recID++ )
        if ( gribapiVarCompare(compVar, streamptr->tsteps[tsID].records[recID], 0) == 0 ) break;

      if ( recID == nrecords )
	{
	  gribWarning("Parameter not defined at timestep 1!", nrecs_scanned, tsID+1, varname, param, level1, level2);
	  return (CDI_EUFSTRUCT);
	}

      if ( streamptr->tsteps[tsID].records[recID].used )
        {
          if ( cdiInventoryMode == 1 ) break;
          else
	    {
	      if ( memcmp(&datetime, &datetime0, sizeof(DateTime)) != 0 ) break;

              gribWarning("Parameter already exist, skipped!", nrecs_scanned, tsID+1, varname, param, level1, level2);
	      continue;
	    }
	}

      streamptr->tsteps[tsID].records[recID].used = TRUE;
      streamptr->tsteps[tsID].recIDs[rindex] = recID;

      if ( CDI_Debug )
        {
          char paramstr[32];
          cdiParamToString(param, paramstr, sizeof(paramstr));
          Message("%4d %8d name=%s id=%s ltype=%d lev1=%d lev2=%d vdate=%d vtime=%d",
                  nrecs_scanned, (int)recpos, varname, paramstr, leveltype1, level1, level2, vdate, vtime);
        }

      streamptr->tsteps[tsID].records[recID].size = recsize;

      if ( gribapiVarCompare(compVar, streamptr->tsteps[tsID].records[recID], 0) != 0 )
	{
	  Message("tsID = %d recID = %d param = %3d new %3d  level = %3d new %3d",
		  tsID, recID,
		  streamptr->tsteps[tsID].records[recID].param, param,
		  streamptr->tsteps[tsID].records[recID].ilevel, level1);
	  return (CDI_EUFSTRUCT);
	}

      streamptr->tsteps[1].records[recID].position = recpos;
      int varID = streamptr->tsteps[tsID].records[recID].varID;
      /*
      gridID = vlistInqVarGrid(vlistID, varID);
      if ( gridInqSize(gridID) == 1 && gridInqType(gridID) == GRID_LONLAT )
	{
	  if ( IS_NOT_EQUAL(gridInqXval(gridID, 0),ISEC2_FirstLon*0.001) ||
	       IS_NOT_EQUAL(gridInqYval(gridID, 0),ISEC2_FirstLat*0.001) )
	    gridChangeType(gridID, GRID_TRAJECTORY);
	}
      */
      if ( tsteptype != vlistInqVarTsteptype(vlistID, varID) )
	vlistDefVarTsteptype(vlistID, varID, tsteptype);

      grib_handle_delete(gh);
      gh = NULL;

      rindex++;
    }

  if ( gh ) grib_handle_delete(gh);

  int nrecs = 0;
  for ( recID = 0; recID < nrecords; recID++ )
    {
      if ( ! streamptr->tsteps[tsID].records[recID].used )
	{
	  int varID = streamptr->tsteps[tsID].records[recID].varID;
	  vlistDefVarTsteptype(vlistID, varID, TSTEP_CONSTANT);
	}
      else
	{
	  nrecs++;
	}
    }
  streamptr->tsteps[tsID].nrecs = nrecs;

  streamptr->rtsteps = 2;

  if ( streamptr->ntsteps == -1 )
    {
      tsID = tstepsNewEntry(streamptr);
      if ( tsID != streamptr->rtsteps )
	Error("Internal error. tsID = %d", tsID);

      streamptr->tsteps[tsID-1].next   = TRUE;
      streamptr->tsteps[tsID].position = recpos;
    }

  streamptr->record->buffer     = gribbuffer;
  streamptr->record->buffersize = buffersize;

  return (rstatus);
}
#endif


#if  defined  (HAVE_LIBGRIB_API)
int gribapiScanTimestep(stream_t * streamptr)
{
  int vrecID, recID;
  //int warn_numavg = TRUE;
  int nrecs = 0;
  int vlistID = streamptr->vlistID;

  if ( CDI_Debug )
    {
      Message("streamID = %d", streamptr->self);
      Message("cts = %d", streamptr->curTsID);
      Message("rts = %d", streamptr->rtsteps);
      Message("nts = %d", streamptr->ntsteps);
    }

  int tsID  = streamptr->rtsteps;
  taxis_t *taxis = &streamptr->tsteps[tsID].taxis;

  if ( streamptr->tsteps[tsID].recordSize == 0 )
    {
      unsigned char* gribbuffer = (unsigned char *) streamptr->record->buffer;
      size_t buffersize = streamptr->record->buffersize;

      cdi_create_records(streamptr, tsID);

      nrecs = streamptr->tsteps[1].nrecs;

      streamptr->tsteps[tsID].nrecs = nrecs;
      streamptr->tsteps[tsID].recIDs = (int *) malloc((size_t)nrecs*sizeof(int));
      for ( recID = 0; recID < nrecs; recID++ )
	streamptr->tsteps[tsID].recIDs[recID] = streamptr->tsteps[1].recIDs[recID];

      int fileID = streamptr->fileID;

      fileSetPos(fileID, streamptr->tsteps[tsID].position, SEEK_SET);

      int nrecs_scanned = streamptr->tsteps[0].nallrecs + streamptr->tsteps[1].nrecs*(tsID-1);    //Only used for debug output.
      int rindex = 0;
      off_t recpos = 0;
      DateTime datetime0;
      grib_handle *gh = NULL;
      char varname[256];
      while ( TRUE )
	{
	  if ( rindex > nrecs ) break;

	  size_t recsize = (size_t)gribGetSize(fileID);
	  recpos  = fileGetPos(fileID);
	  if ( recsize == 0 )
	    {
	      streamptr->ntsteps = streamptr->rtsteps + 1;
	      break;
	    }

	  if ( rindex >= nrecs ) break;

          ensureBufferSize(recsize, &buffersize, &gribbuffer);

	  size_t readsize = recsize;
	  if (gribRead(fileID, gribbuffer, &readsize))
	    {
	      Warning("Inconsistent timestep %d (GRIB record %d/%d)!", tsID+1, rindex+1,
		      streamptr->tsteps[tsID].recordSize);
	      break;
	    }

          long unzipsize;
	  if ( gribGetZip((long)recsize, gribbuffer, &unzipsize) > 0 )
            ensureBufferSize((size_t)unzipsize + 100, &buffersize, &gribbuffer);

          nrecs_scanned++;
	  gh = grib_handle_new_from_message(NULL, (void *) gribbuffer, recsize);
	  GRIB_CHECK(my_grib_set_double(gh, "missingValue", cdiDefaultMissval), 0);

          int param = gribapiGetParam(gh);
          int level1 = 0, level2 = 0, leveltype1, leveltype2 = -1, lbounds, level_sf, level_unit;
          var_tile_t tiles = dummy_tiles;
          gribGetLevel(gh, &leveltype1, &leveltype2, &lbounds, &level1, &level2, &level_sf, &level_unit, &tiles);

          varname[0] = 0;
	  gribapiGetString(gh, "shortName", varname, sizeof(varname));

          int vdate = 0, vtime = 0;
	  gribapiGetValidityDateTime(gh, &vdate, &vtime);

	  if ( rindex == nrecs ) break;

	  if ( rindex == 0 )
	    {
              int taxisID = vlistInqTaxis(vlistID);
	      if ( taxisInqType(taxisID) == TAXIS_RELATIVE )
		{
		  taxis->type  = TAXIS_RELATIVE;

                  gribapiGetDataDateTime(gh, &(taxis->rdate), &(taxis->rtime));

		  taxis->unit  = gribapiGetTimeUnits(gh);
		}
	      else
		{
		  taxis->type  = TAXIS_ABSOLUTE;
		}
	      taxis->vdate = vdate;
	      taxis->vtime = vtime;

	      datetime0.date = vdate;
	      datetime0.time = vtime;
	    }
	  /*
	  if ( ISEC1_AvgNum )
	    {
	      if (  taxis->numavg && warn_numavg &&
		   (taxis->numavg != ISEC1_AvgNum) )
		{
		  warn_numavg = FALSE;
		}
	      else
		{
		  taxis->numavg = ISEC1_AvgNum;
		}
	    }
	  */
          DateTime datetime = {
            .date  = vdate,
            .time  = vtime
          };

          int tsteptype = gribapiGetTsteptype(gh);

          compvar2_t compVar = gribapiVarSet(param, level1, level2, leveltype1, tsteptype, varname, tiles);

	  for ( vrecID = 0; vrecID < nrecs; vrecID++ )
	    {
	      recID   = streamptr->tsteps[1].recIDs[vrecID];
	      if ( gribapiVarCompare(compVar, streamptr->tsteps[tsID].records[recID], 0) == 0 ) break;
	    }

	  if ( vrecID == nrecs )
	    {
	      gribWarning("Parameter not defined at timestep 1!", nrecs_scanned, tsID+1, varname, param, level1, level2);

	      if ( cdiInventoryMode == 1 )
		return (CDI_EUFSTRUCT);
	      else
		continue;
	    }

	  if ( cdiInventoryMode != 1 )
	    {
	      if ( streamptr->tsteps[tsID].records[recID].used )
		{
		  if ( memcmp(&datetime, &datetime0, sizeof(DateTime)) != 0 ) break;

		  if ( CDI_Debug )
                    gribWarning("Parameter already exist, skipped!", nrecs_scanned, tsID+1, varname, param, level1, level2);

		  continue;
		}
	    }

          streamptr->tsteps[tsID].records[recID].used = TRUE;
          streamptr->tsteps[tsID].recIDs[rindex] = recID;

	  if ( CDI_Debug )
	    Message("%4d %8d %4d %8d %8d %6d", rindex+1, (int)recpos, param, level1, vdate, vtime);

	  if ( gribapiVarCompare(compVar, streamptr->tsteps[tsID].records[recID], 0) != 0 )
	    {
	      Message("tsID = %d recID = %d param = %3d new %3d  level = %3d new %3d",
		      tsID, recID,
		      streamptr->tsteps[tsID].records[recID].param, param,
		      streamptr->tsteps[tsID].records[recID].ilevel, level1);
	      Error("Invalid, unsupported or inconsistent record structure");
	    }

	  streamptr->tsteps[tsID].records[recID].position = recpos;
	  streamptr->tsteps[tsID].records[recID].size = recsize;

	  if ( CDI_Debug )
	    Message("%4d %8d %4d %8d %8d %6d", rindex, (int)recpos, param, level1, vdate, vtime);

	  grib_handle_delete(gh);
	  gh = NULL;

	  rindex++;
	}

      if ( gh ) grib_handle_delete(gh);

      for ( vrecID = 0; vrecID < nrecs; vrecID++ )
	{
	  recID   = streamptr->tsteps[tsID].recIDs[vrecID];
	  if ( ! streamptr->tsteps[tsID].records[recID].used ) break;
	}

      if ( vrecID < nrecs )
	{
	  gribWarning("Paramameter not found!", nrecs_scanned, tsID+1, varname, streamptr->tsteps[tsID].records[recID].param,
                      streamptr->tsteps[tsID].records[recID].ilevel, streamptr->tsteps[tsID].records[recID].ilevel2);
	  return (CDI_EUFSTRUCT);
	}

      streamptr->rtsteps++;

      if ( streamptr->ntsteps != streamptr->rtsteps )
	{
	  tsID = tstepsNewEntry(streamptr);
	  if ( tsID != streamptr->rtsteps )
	    Error("Internal error. tsID = %d", tsID);

	  streamptr->tsteps[tsID-1].next   = 1;
	  streamptr->tsteps[tsID].position = recpos;
	}

      fileSetPos(fileID, streamptr->tsteps[tsID].position, SEEK_SET);
      streamptr->tsteps[tsID].position = recpos;

      streamptr->record->buffer     = gribbuffer;
      streamptr->record->buffersize = buffersize;
    }

  if ( nrecs > 0 && nrecs < streamptr->tsteps[tsID].nrecs )
    {
      Warning("Incomplete timestep. Stop scanning at timestep %d.", tsID);
      streamptr->ntsteps = tsID;
    }

  return (int)streamptr->ntsteps;
}
#endif

#ifdef gribWarning
#undef gribWarning
#endif

#ifdef HAVE_LIBGRIB_API
int gribapiDecode(unsigned char *gribbuffer, int gribsize, double *data, int gridsize,
		  int unreduced, int *nmiss, double missval, int vlistID, int varID)
{
  int status = 0;
  long lpar;
  long numberOfPoints;
  size_t datasize, dummy, recsize;
  grib_handle *gh = NULL;

  UNUSED(vlistID);
  UNUSED(varID);

  if ( unreduced )
    {
      static int lwarn = 1;

      if ( lwarn )
	{
	  lwarn = 0;
	  Warning("Conversion of gaussian reduced grids unsupported!");
	}
    }

  recsize = (size_t)gribsize;
  gh = grib_handle_new_from_message(NULL, (void *) gribbuffer, recsize);
  GRIB_CHECK(my_grib_set_double(gh, "missingValue", missval), 0);

  /* get the size of the values array*/
  GRIB_CHECK(grib_get_size(gh, "values", &datasize), 0);
  GRIB_CHECK(grib_get_long(gh, "numberOfPoints", &numberOfPoints), 0);

  // printf("values_size = %d  numberOfPoints = %ld\n", datasize, numberOfPoints);

  if ( gridsize != (int) datasize )
    Error("Internal problem: gridsize(%d) != datasize(%d)!", gridsize, datasize);
  dummy = datasize;
  GRIB_CHECK(grib_get_double_array(gh, "values", data, &dummy), 0);

  int gridtype;
  GRIB_CHECK(grib_get_long(gh, "gridDefinitionTemplateNumber", &lpar), 0);
  gridtype = (int) lpar;

  *nmiss = 0;
  if ( gridtype < 50 || gridtype > 53 )
    {
      GRIB_CHECK(grib_get_long(gh, "numberOfMissing", &lpar), 0);
      *nmiss = (int) lpar;
      // printf("gridtype %d, nmiss %d\n", gridtype, nmiss);
    }

  grib_handle_delete(gh);

  return (status);
}
#endif


#if  defined  (HAVE_LIBGRIB_API)
static
void gribapiDefInstitut(grib_handle *gh, int vlistID, int varID)
{
  int instID;

  if ( vlistInqInstitut(vlistID) != CDI_UNDEFID )
    instID = vlistInqInstitut(vlistID);
  else
    instID = vlistInqVarInstitut(vlistID, varID);

  if ( instID != CDI_UNDEFID )
    {
      long center, subcenter;
      long center0, subcenter0;

      center    = institutInqCenter(instID);
      subcenter = institutInqSubcenter(instID);

      GRIB_CHECK(grib_get_long(gh, "centre", &center0), 0);
      GRIB_CHECK(grib_get_long(gh, "subCentre", &subcenter0), 0);

      if ( center != center0 )
	GRIB_CHECK(my_grib_set_long(gh, "centre", center), 0);
      if ( subcenter != subcenter0 )
	GRIB_CHECK(my_grib_set_long(gh, "subCentre", subcenter), 0);
    }
}

static
void gribapiDefModel(grib_handle *gh, int vlistID, int varID)
{
  int modelID;

  if ( vlistInqModel(vlistID) != CDI_UNDEFID )
    modelID = vlistInqModel(vlistID);
  else
    modelID = vlistInqVarModel(vlistID, varID);

  if ( modelID != CDI_UNDEFID )
    GRIB_CHECK(my_grib_set_long(gh, "generatingProcessIdentifier", modelInqGribID(modelID)), 0);
}

static
void gribapiDefParam(int editionNumber, grib_handle *gh, int param, const char *name, const char *stdname)
{
  bool ldefined = false;

  int pdis, pcat, pnum;
  cdiDecodeParam(param, &pnum, &pcat, &pdis);

  if ( pnum < 0 )
    {
      size_t len;
      len = strlen(stdname);
      if ( len )
        {
          int status = my_grib_set_string(gh, "cfName", stdname, &len);
          if ( status == 0 ) ldefined = true;
          else Warning("grib_api: No match for cfName=%s", stdname);
        }

      if ( ldefined == false )
        {
          len = strlen(name);
          int status = my_grib_set_string(gh, "shortName", name, &len);
          if ( status == 0 ) ldefined = true;
          else Warning("grib_api: No match for shortName=%s", name);
        }
    }

  if ( ldefined == false )
    {
      if ( pnum < 0 ) pnum = -pnum;

      static bool lwarn_pnum = true;
      if ( pnum > 255 && lwarn_pnum )
        {
          Warning("Parameter number %d out of range (1-255), set to %d!", pnum, pnum%256);
          lwarn_pnum = false;
          pnum = pnum%256;
        }

      if ( editionNumber <= 1 )
	{
          static bool lwarn_pdis = true;
	  if ( pdis != 255 && lwarn_pdis )
	    {
	      char paramstr[32];
	      cdiParamToString(param, paramstr, sizeof(paramstr));
	      Warning("Can't convert GRIB2 parameter ID (%s) to GRIB1, set to %d.%d!", paramstr, pnum, pcat);
              lwarn_pdis = false;
	    }

	  GRIB_CHECK(my_grib_set_long(gh, "table2Version",        pcat), 0);
	  GRIB_CHECK(my_grib_set_long(gh, "indicatorOfParameter", pnum), 0);
	}
      else
	{
	  GRIB_CHECK(my_grib_set_long(gh, "discipline",        pdis), 0);
	  GRIB_CHECK(my_grib_set_long(gh, "parameterCategory", pcat), 0);
	  GRIB_CHECK(my_grib_set_long(gh, "parameterNumber",   pnum), 0);
	}
    }

  // printf("param: %d.%d.%d %s\n", pnum, pcat, pdis, name);
}

static
int getTimeunitFactor(int timeunit)
{
  int factor = 1;

  switch (timeunit)
    {
    case TUNIT_SECOND:  factor =     1;  break;
    case TUNIT_MINUTE:  factor =    60;  break;
    case TUNIT_HOUR:    factor =  3600;  break;
    case TUNIT_3HOURS:  factor = 10800;  break;
    case TUNIT_6HOURS:  factor = 21600;  break;
    case TUNIT_12HOURS: factor = 43200;  break;
    case TUNIT_DAY:     factor = 86400;  break;
    default:            factor =  3600;  break;
    }

  return (factor);
}

static
void gribapiDefStepUnits(grib_handle *gh, int timeunit, int proDefTempNum, int gcinit)
{
  long unitsOfTime;

  switch (timeunit)
    {
    case TUNIT_SECOND:  unitsOfTime = 13;  break;
    case TUNIT_MINUTE:  unitsOfTime =  0;  break;
    case TUNIT_HOUR:    unitsOfTime =  1;  break;
    case TUNIT_3HOURS:  unitsOfTime = 10;  break;
    case TUNIT_6HOURS:  unitsOfTime = 11;  break;
    case TUNIT_12HOURS: unitsOfTime = 12;  break;
    case TUNIT_DAY:     unitsOfTime =  2;  break;
    default:            unitsOfTime =  1;  break;
    }

  if ( !gcinit )
    {
      GRIB_CHECK(my_grib_set_long(gh, "stepUnits", unitsOfTime), 0);
      if ( proDefTempNum == 8 || proDefTempNum == 11 )
        GRIB_CHECK(my_grib_set_long(gh, "indicatorOfUnitForTimeRange", unitsOfTime), 0);
      GRIB_CHECK(my_grib_set_long(gh, "indicatorOfUnitOfTimeRange", unitsOfTime), 0);
    }
}

static
int gribapiDefSteptype(int editionNumber, grib_handle *gh, int productDefinitionTemplate, int typeOfGeneratingProcess, int tsteptype, int gcinit)
{
  long proDefTempNum = 0;
  size_t len = 64;
  char stepType[len];

  switch ( tsteptype )
    {
    case TSTEP_AVG:      strcpy(stepType, "avg");     proDefTempNum = 8; break;
    case TSTEP_ACCUM:    strcpy(stepType, "accum");   proDefTempNum = 8; break;
    case TSTEP_MAX:      strcpy(stepType, "max");     proDefTempNum = 8; break;
    case TSTEP_MIN:      strcpy(stepType, "min");     proDefTempNum = 8; break;
    case TSTEP_DIFF:     strcpy(stepType, "diff");    proDefTempNum = 8; break;
    case TSTEP_RMS:      strcpy(stepType, "rms");     proDefTempNum = 8; break;
    case TSTEP_SD:       strcpy(stepType, "sd");      proDefTempNum = 8; break;
    case TSTEP_COV:      strcpy(stepType, "cov");     proDefTempNum = 8; break;
    case TSTEP_RATIO:    strcpy(stepType, "ratio");   proDefTempNum = 8; break;
    case TSTEP_INSTANT:  strcpy(stepType, "instant"); proDefTempNum = 0; break;
    default:             strcpy(stepType, "instant"); proDefTempNum = 0; break;
    }

  if ( typeOfGeneratingProcess == 4 )
    {
      if ( proDefTempNum == 8 ) proDefTempNum = 11;
      else                      proDefTempNum = 1;
    }

  if ( productDefinitionTemplate != -1 ) proDefTempNum = productDefinitionTemplate;

  if ( !gcinit )
    {
      if ( editionNumber > 1 ) GRIB_CHECK(my_grib_set_long(gh, "productDefinitionTemplateNumber", proDefTempNum), 0);
      len = strlen(stepType);
      GRIB_CHECK(my_grib_set_string(gh, "stepType", stepType, &len), 0);
    }

  return ((int)proDefTempNum);
}

static
void gribapiDefDateTimeAbs(int editionNumber, grib_handle *gh, int date, int time, int productDefinitionTemplate, int typeOfGeneratingProcess, int tsteptype, int gcinit)
{
  (void ) gribapiDefSteptype(editionNumber, gh, productDefinitionTemplate, typeOfGeneratingProcess, tsteptype, gcinit);

  if ( editionNumber > 1 ) GRIB_CHECK(my_grib_set_long(gh, "significanceOfReferenceTime", 0), 0);
  if ( editionNumber > 1 ) GRIB_CHECK(my_grib_set_long(gh, "stepRange", 0), 0);

  if ( date == 0 ) date = 10101;
  gribapiSetDataDateTime(gh, date, time);
}

static
int gribapiDefDateTimeRel(int editionNumber, grib_handle *gh, int rdate, int rtime, int vdate, int vtime,
                          int productDefinitionTemplate, int typeOfGeneratingProcess, int tsteptype, int timeunit, int calendar, int gcinit)
{
  int status = -1;
  int year, month, day, hour, minute, second;
  int julday1, secofday1, julday2, secofday2, days, secs;
  int factor;
  long startStep = 0, endStep;

  cdiDecodeDate(rdate, &year, &month, &day);
  cdiDecodeTime(rtime, &hour, &minute, &second);
  encode_juldaysec(calendar, year, month, day, hour, minute, second, &julday1, &secofday1);

  if ( vdate == 0 && vtime == 0 ) { vdate = rdate; vtime = rtime; }

  cdiDecodeDate(vdate, &year, &month, &day);
  cdiDecodeTime(vtime, &hour, &minute, &second);
  encode_juldaysec(calendar, year, month, day, hour, minute, second, &julday2, &secofday2);

  (void) julday_sub(julday1, secofday1, julday2, secofday2, &days, &secs);

  factor = getTimeunitFactor(timeunit);

  if ( !(int) fmod(days*86400.0 + secs, factor) )
    {
      int proDefTempNum = gribapiDefSteptype(editionNumber, gh, productDefinitionTemplate, typeOfGeneratingProcess, tsteptype, gcinit);

      gribapiDefStepUnits(gh, timeunit, proDefTempNum, gcinit);

      endStep = (int) ((days*86400.0 + secs)/factor);

      if ( editionNumber > 1 ) GRIB_CHECK(my_grib_set_long(gh, "significanceOfReferenceTime", 1), 0);
      if ( editionNumber > 1 ) GRIB_CHECK(my_grib_set_long(gh, "stepRange", 0), 0);

      if ( rdate == 0 ) rdate = 10101;
      gribapiSetDataDateTime(gh, rdate, rtime);

      // printf(">>>>> tsteptype %d  startStep %ld  endStep %ld\n", tsteptype, startStep, endStep);

      // Product Definition Template Number: defined in GRIB_API file 4.0.table
      // point in time products:
      if ( (proDefTempNum >= 0 && proDefTempNum <=  7) || 
           proDefTempNum == 55 || proDefTempNum == 40055 ) // Tile
        startStep = endStep;

      if ( editionNumber > 1 ) GRIB_CHECK(my_grib_set_long(gh, "forecastTime", startStep), 0);
      GRIB_CHECK(my_grib_set_long(gh, "endStep", endStep), 0);

      status = 0;
    }

  return (status);
}

static
void gribapiDefTime(int editionNumber, int productDefinitionTemplate, int typeOfGeneratingProcess, grib_handle *gh,
                    int vdate, int vtime, int tsteptype, int numavg, int taxisID, int gcinit)
{
  int taxistype = -1;

  UNUSED(numavg);

  if ( taxisID != -1 ) taxistype = taxisInqType(taxisID);

  if ( typeOfGeneratingProcess == 196 )
    {
      vdate = 10101;
      vtime = 0;
      taxistype = TAXIS_ABSOLUTE;
    }
  /*
  else if ( typeOfGeneratingProcess == 9 )
    {
    }
  */

  if ( taxistype == TAXIS_RELATIVE )
    {
      int status;
      int calendar = taxisInqCalendar(taxisID);
      int rdate    = taxisInqRdate(taxisID);
      int rtime    = taxisInqRtime(taxisID);
      int timeunit = taxisInqTunit(taxisID);

      status = gribapiDefDateTimeRel(editionNumber, gh, rdate, rtime, vdate, vtime,
                                     productDefinitionTemplate, typeOfGeneratingProcess, tsteptype, timeunit, calendar, gcinit);

      if ( status != 0 ) taxistype = TAXIS_ABSOLUTE;
    }

  if ( taxistype == TAXIS_ABSOLUTE )
    {
      gribapiDefDateTimeAbs(editionNumber, gh, vdate, vtime, productDefinitionTemplate, typeOfGeneratingProcess, tsteptype, gcinit);
    }
}

static
void gribapiDefGrid(int editionNumber, grib_handle *gh, int gridID, int comptype, int lieee, int datatype, int nmiss, int gcinit)
{
  int gridtype;
  int status;
  static short lwarn = TRUE;
  size_t len;
  char *mesg;

  UNUSED(nmiss);

  gridtype = gridInqType(gridID);

  if ( editionNumber <= 1 )
    if ( gridtype == GRID_GME || gridtype == GRID_UNSTRUCTURED )
      gridtype = -1;

  if ( gridtype == GRID_GENERIC )
    {
      int xsize, ysize, gridsize;

      gridsize = gridInqSize(gridID);
      xsize = gridInqXsize(gridID);
      ysize = gridInqYsize(gridID);

      if ( (ysize ==  32 || ysize ==  48 || ysize ==  64 ||
	    ysize ==  96 || ysize == 160 || ysize == 192 ||
	    ysize == 240 || ysize == 320 || ysize == 384 ||
	    ysize == 480 || ysize == 768 ) &&
	   (xsize == 2*ysize || xsize == 1) )
	{
	  gridtype = GRID_GAUSSIAN;
	  gridChangeType(gridID, gridtype);
	}
      else if ( gridsize == 1 )
	{
	  gridtype = GRID_LONLAT;
	  gridChangeType(gridID, gridtype);
	}
      else if ( gridInqXvals(gridID, NULL) && gridInqYvals(gridID, NULL) )
	{
	  gridtype = GRID_LONLAT;
	  gridChangeType(gridID, gridtype);
	}
    }
  else if ( gridtype == GRID_CURVILINEAR )
    {
      if ( lwarn && gridInqSize(gridID) > 1 )
	{
	  lwarn = FALSE;
	  Warning("Curvilinear grids are unsupported in GRIB format! Created wrong GDS!");
	}
      gridtype = GRID_LONLAT;
    }

  if ( gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN )
    {
      if ( editionNumber != 2 || lieee ) { comptype = 0; }

      if ( comptype )
        {
          if ( comptype == COMPRESS_JPEG )
            {
              mesg = "grid_jpeg"; len = strlen(mesg);
              GRIB_CHECK(my_grib_set_string(gh, "packingType", mesg, &len), 0);
            }
          else if ( comptype == COMPRESS_SZIP )
            {
              mesg = "grid_ccsds"; len = strlen(mesg);
              GRIB_CHECK(my_grib_set_string(gh, "packingType", mesg, &len), 0);
            }
          else
            {
              mesg = "grid_simple"; len = strlen(mesg);
              GRIB_CHECK(my_grib_set_string(gh, "packingType", mesg, &len), 0);
            }
        }
    }

  if ( gcinit ) return;

  switch (gridtype)
    {
    case GRID_LONLAT:
    case GRID_GAUSSIAN:
    case GRID_GAUSSIAN_REDUCED:
    case GRID_TRAJECTORY:
      {
	int nlon = 0, nlat;
	double xfirst = 0, xlast = 0, xinc = 0;
	double yfirst = 0, ylast = 0, yinc = 0;
	double latIncr;

	if ( gridtype == GRID_GAUSSIAN )
	  {
	    mesg = "regular_gg"; len = strlen(mesg);
	    GRIB_CHECK(my_grib_set_string(gh, "gridType", mesg, &len), 0);
	  }
	else if ( gridtype == GRID_GAUSSIAN_REDUCED )
	  {
	    mesg = "reduced_gg"; len = strlen(mesg);
	    GRIB_CHECK(my_grib_set_string(gh, "gridType", mesg, &len), 0);
	  }
	else if ( gridtype == GRID_LONLAT && gridIsRotated(gridID) )
	  {
	    mesg = "rotated_ll"; len = strlen(mesg);
	    GRIB_CHECK(my_grib_set_string(gh, "gridType", mesg, &len), 0);
	  }
	else
	  {
	    mesg = "regular_ll"; len = strlen(mesg);
	    GRIB_CHECK(my_grib_set_string(gh, "gridType", mesg, &len), 0);
	  }

	nlon = gridInqXsize(gridID);
	nlat = gridInqYsize(gridID);

	if ( gridtype == GRID_GAUSSIAN_REDUCED )
	  {
	    int *rowlon, i;
	    long *pl = NULL;

	    nlon = 0;

	    rowlon = (int *)xmalloc((size_t)nlat*sizeof(int));
	    pl     = (long *)xmalloc((size_t)nlat*sizeof(long));
	    gridInqRowlon(gridID, rowlon);
	    for ( i = 0; i < nlat; ++i ) pl[i] = rowlon[i];

	    // GRIB_CHECK(my_grib_set_long_array(gh, "pl", pl, nlat), 0);

	    free(pl);
	    free(rowlon);
	  }
	else
	  {
	    if ( nlon == 0 )
	      {
		nlon = 1;
	      }
	    else
	      {
		xfirst = gridInqXval(gridID,      0);
		xlast  = gridInqXval(gridID, nlon-1);
		xinc   = gridInqXinc(gridID);
	      }
	  }

	if ( nlat == 0 )
	  {
	    nlat = 1;
	  }
	else
	  {
	    yfirst = gridInqYval(gridID,      0);
	    ylast  = gridInqYval(gridID, nlat-1);
	    yinc   = gridInqYinc(gridID);
	  }

	GRIB_CHECK(my_grib_set_long(gh, "Ni", nlon), 0);
	GRIB_CHECK(my_grib_set_long(gh, "Nj", nlat), 0);
	GRIB_CHECK(my_grib_set_double(gh, "longitudeOfFirstGridPointInDegrees", xfirst), 0);
	GRIB_CHECK(my_grib_set_double(gh, "longitudeOfLastGridPointInDegrees",  xlast), 0);
	GRIB_CHECK(my_grib_set_double(gh, "latitudeOfFirstGridPointInDegrees",  yfirst), 0);
	GRIB_CHECK(my_grib_set_double(gh, "latitudeOfLastGridPointInDegrees",   ylast), 0);
	GRIB_CHECK(my_grib_set_double(gh, "iDirectionIncrementInDegrees", xinc), 0);

        {
          long jscan = 0;
          if ( yfirst < ylast ) jscan = 1;
          GRIB_CHECK(my_grib_set_long(gh, "jScansPositively", jscan), 0);
        }
	/*
	if ( fabs(xinc*1000 - ISEC2_LonIncr) > FLT_EPSILON )
	  ISEC2_LonIncr = 0;
	*/
	if ( gridtype == GRID_GAUSSIAN || gridtype == GRID_GAUSSIAN_REDUCED )
          {
            int np = gridInqNP(gridID);
            if ( np == 0 ) np = nlat/2;
            GRIB_CHECK(my_grib_set_long(gh, "numberOfParallelsBetweenAPoleAndTheEquator", np), 0);
          }
	else
	  {
	    latIncr = yinc;
	    if ( latIncr < 0 ) latIncr = -latIncr;
	    GRIB_CHECK(my_grib_set_double(gh, "jDirectionIncrementInDegrees", latIncr), 0);
	    /*
	    if ( fabs(yinc*1000 - ISEC2_LatIncr) > FLT_EPSILON )
	      ISEC2_LatIncr = 0;
	    */
	  }
	/*
	if ( ISEC2_NumLon > 1 && ISEC2_NumLat == 1 )
	  if ( ISEC2_LonIncr != 0 && ISEC2_LatIncr == 0 ) ISEC2_LatIncr = ISEC2_LonIncr;

	if ( ISEC2_NumLon == 1 && ISEC2_NumLat > 1 )
	  if ( ISEC2_LonIncr == 0 && ISEC2_LatIncr != 0 ) ISEC2_LonIncr = ISEC2_LatIncr;

	if ( ISEC2_LatIncr == 0 || ISEC2_LonIncr == 0 )
	  ISEC2_ResFlag = 0;
	else
	  ISEC2_ResFlag = 128;
	*/
	if ( gridIsRotated(gridID) )
	  {
	    double xpole, ypole, angle;
	    xpole = gridInqXpole(gridID);
	    ypole = gridInqYpole(gridID);
	    angle = gridInqAngle(gridID);
	    /* change from north to south pole */
	    ypole = -ypole;
	    xpole =  xpole + 180;
	    GRIB_CHECK(my_grib_set_double(gh, "latitudeOfSouthernPoleInDegrees",  ypole), 0);
	    GRIB_CHECK(my_grib_set_double(gh, "longitudeOfSouthernPoleInDegrees", xpole), 0);
	    GRIB_CHECK(my_grib_set_double(gh, "angleOfRotation", angle), 0);
	  }

	/* East -> West */
	//if ( ISEC2_LastLon < ISEC2_FirstLon ) ISEC2_ScanFlag += 128;

	/* South -> North */
	//if ( ISEC2_LastLat > ISEC2_FirstLat ) ISEC2_ScanFlag += 64;

        if ( editionNumber != 2 ) { lieee = 0; comptype = 0; }

        if ( lieee )
          {
            mesg = "grid_ieee"; len = strlen(mesg);
            GRIB_CHECK(my_grib_set_string(gh, "packingType", mesg, &len), 0);

	    if ( datatype == DATATYPE_FLT64 )
	      GRIB_CHECK(my_grib_set_long(gh, "precision", 2), 0);
	    else
	      GRIB_CHECK(my_grib_set_long(gh, "precision", 1), 0);
          }
        else
	  {
            if ( comptype == COMPRESS_JPEG )
              {
                mesg = "grid_jpeg"; len = strlen(mesg);
                GRIB_CHECK(my_grib_set_string(gh, "packingType", mesg, &len), 0);
              }
            else if ( comptype == COMPRESS_SZIP )
              {
                mesg = "grid_ccsds"; len = strlen(mesg);
                GRIB_CHECK(my_grib_set_string(gh, "packingType", mesg, &len), 0);
              }
            else
              {
                mesg = "grid_simple"; len = strlen(mesg);
                GRIB_CHECK(my_grib_set_string(gh, "packingType", mesg, &len), 0);
              }
	  }

	break;
      }
    case GRID_LCC:
      {
	double originLon, originLat, lonParY, lat1, lat2, xincm, yincm;
	int xsize, ysize;
	int projflag, scanflag;

	xsize = gridInqXsize(gridID);
	ysize = gridInqYsize(gridID);

	gridInqLCC(gridID, &originLon, &originLat, &lonParY, &lat1, &lat2, &xincm, &yincm,
		   &projflag, &scanflag);

        mesg = "lambert"; len = strlen(mesg);
        GRIB_CHECK(my_grib_set_string(gh, "gridType", mesg, &len), 0);

	GRIB_CHECK(my_grib_set_long(gh, "Nx", xsize), 0);
	GRIB_CHECK(my_grib_set_long(gh, "Ny", ysize), 0);

        /* FIXME: lround should probably be round here */
	GRIB_CHECK(my_grib_set_double(gh, "DxInMetres", (double)lround(xincm)), 0);
        /* FIXME: lround should probably be round here */
	GRIB_CHECK(my_grib_set_double(gh, "DyInMetres", (double)lround(yincm)), 0);
	GRIB_CHECK(my_grib_set_double(gh, "longitudeOfFirstGridPointInDegrees", originLon), 0);
	GRIB_CHECK(my_grib_set_double(gh, "latitudeOfFirstGridPointInDegrees", originLat), 0);
	GRIB_CHECK(my_grib_set_double(gh, "LoVInDegrees", lonParY), 0);
	GRIB_CHECK(my_grib_set_double(gh, "Latin1InDegrees", lat1), 0);
	GRIB_CHECK(my_grib_set_double(gh, "Latin2InDegrees", lat2), 0);

        if ( editionNumber <= 1 )
          {
            GRIB_CHECK(my_grib_set_long(gh, "projectionCenterFlag", projflag), 0);
            GRIB_CHECK(my_grib_set_long(gh, "scanningMode", scanflag), 0);
          }
        /*
	ISEC2_Lambert_LatSP  = 0;
	ISEC2_Lambert_LatSP  = 0;
        */
	break;
      }
    case GRID_SPECTRAL:
      {
	int trunc = gridInqTrunc(gridID);

	mesg = "sh"; len = strlen(mesg);
	GRIB_CHECK(my_grib_set_string(gh, "gridType", mesg, &len), 0);

	GRIB_CHECK(my_grib_set_long(gh, "J", trunc), 0);
	GRIB_CHECK(my_grib_set_long(gh, "K", trunc), 0);
	GRIB_CHECK(my_grib_set_long(gh, "M", trunc), 0);

	// GRIB_CHECK(my_grib_set_long(gh, "numberOfDataPoints", gridInqSize(gridID)), 0);
        /*
        if ( lieee )
          {
            printf("spectral_ieee\n");
            if ( editionNumber == 2 ) GRIB_CHECK(my_grib_set_long(gh, "numberOfValues", gridInqSize(gridID)), 0);
            mesg = "spectral_ieee"; len = strlen(mesg);
            GRIB_CHECK(my_grib_set_string(gh, "packingType", mesg, &len), 0);
          }
        else */ if ( gridInqComplexPacking(gridID) )
	  {
	    if ( editionNumber == 2 ) GRIB_CHECK(my_grib_set_long(gh, "numberOfValues", gridInqSize(gridID)), 0);
	    mesg = "spectral_complex"; len = strlen(mesg);
	    GRIB_CHECK(my_grib_set_string(gh, "packingType", mesg, &len), 0);

	    GRIB_CHECK(my_grib_set_long(gh, "JS", 20), 0);
	    GRIB_CHECK(my_grib_set_long(gh, "KS", 20), 0);
	    GRIB_CHECK(my_grib_set_long(gh, "MS", 20), 0);
	  }
	else
	  {
	    mesg = "spectral_simple"; len = strlen(mesg);
	    GRIB_CHECK(my_grib_set_string(gh, "packingType", mesg, &len), 0);
	  }

	break;
      }
    case GRID_GME:
      {
	GRIB_CHECK(my_grib_set_long(gh, "gridDefinitionTemplateNumber", GRIB2_GTYPE_GME), 0);

	GRIB_CHECK(my_grib_set_long(gh, "nd", gridInqGMEnd(gridID)), 0);
	GRIB_CHECK(my_grib_set_long(gh, "Ni", gridInqGMEni(gridID)), 0);
	GRIB_CHECK(my_grib_set_long(gh, "n2", gridInqGMEni2(gridID)), 0);
	GRIB_CHECK(my_grib_set_long(gh, "n3", gridInqGMEni3(gridID)), 0);
	GRIB_CHECK(my_grib_set_long(gh, "latitudeOfThePolePoint", 90000000), 0);
	GRIB_CHECK(my_grib_set_long(gh, "longitudeOfThePolePoint", 0), 0);

	GRIB_CHECK(my_grib_set_long(gh, "numberOfDataPoints", gridInqSize(gridID)), 0);
	GRIB_CHECK(my_grib_set_long(gh, "totalNumberOfGridPoints", gridInqSize(gridID)), 0);

        if ( comptype == COMPRESS_SZIP )
          {
            mesg = "grid_ccsds"; len = strlen(mesg);
            GRIB_CHECK(my_grib_set_string(gh, "packingType", mesg, &len), 0);
          }

	break;
      }
    case GRID_UNSTRUCTURED:
      {
	static int warning = 1;

	status = my_grib_set_long(gh, "gridDefinitionTemplateNumber", GRIB2_GTYPE_UNSTRUCTURED);
	if ( status != 0 && warning )
	  {
	    warning = 0;
	    Warning("Can't write reference grid!");
	    Warning("gridDefinitionTemplateNumber %d not found (grib2/template.3.%d.def)!",
		    GRIB2_GTYPE_UNSTRUCTURED, GRIB2_GTYPE_UNSTRUCTURED);
	  }
	else
	  {
            unsigned char uuid[CDI_UUID_SIZE];
            int position = gridInqPosition(gridID);
            int number = gridInqNumber(gridID);
            if ( position < 0 ) position = 0;
            if ( number < 0 ) number = 0;
	    GRIB_CHECK(my_grib_set_long(gh, "numberOfGridUsed", number), 0);
	    GRIB_CHECK(my_grib_set_long(gh, "numberOfGridInReference", position), 0);
            len = CDI_UUID_SIZE;
            gridInqUUID(gridID, uuid);
	    if (grib_set_bytes(gh, "uuidOfHGrid", uuid, &len) != 0)
	      Warning("Can't write UUID!");
	  }

        if ( comptype == COMPRESS_SZIP )
          {
            mesg = "grid_ccsds"; len = strlen(mesg);
            GRIB_CHECK(my_grib_set_string(gh, "packingType", mesg, &len), 0);
          }

	break;
      }
    default:
      {
	Error("Unsupported grid type: %s", gridNamePtr(gridtype));
	break;
      }
    }
}

static
void getLevelFactor(double level, long *factor, long *out_scaled_value)
{
  double scaled_value  = level;
  /* FIXME: lround might be better here */
  long   iscaled_value = (long) round(scaled_value);
  long   i;

  const double eps = 1.e-8;
  for ( i=0; (fabs(scaled_value - (double) iscaled_value) >= eps) && i < 7; i++ )
    {
      scaled_value *= 10.;
      /* FIXME: lround might be better here */
      iscaled_value = (long)round(scaled_value);
    }

  (*factor)           = i;
  (*out_scaled_value) = iscaled_value;
}

static
void gribapiDefLevelType(grib_handle *gh, int gcinit, const char *keyname, long leveltype)
{
  if ( !gcinit ) GRIB_CHECK(my_grib_set_long(gh, keyname, leveltype), 0);
}

static
void grib2DefLevel(grib_handle *gh, int gcinit, long leveltype1, long leveltype2, int lbounds, double level, double dlevel1, double dlevel2)
{
  long scaled_level;
  long factor;

  gribapiDefLevelType(gh, gcinit, "typeOfFirstFixedSurface", leveltype1);
  if ( lbounds ) gribapiDefLevelType(gh, gcinit, "typeOfSecondFixedSurface", leveltype2);

  if ( !lbounds ) dlevel1 = level;

  getLevelFactor(dlevel1, &factor, &scaled_level);
  GRIB_CHECK(my_grib_set_long(gh, "scaleFactorOfFirstFixedSurface", factor), 0);
  GRIB_CHECK(my_grib_set_long(gh, "scaledValueOfFirstFixedSurface", scaled_level), 0);

  if ( lbounds )
    {
      getLevelFactor(dlevel2, &factor, &scaled_level);
      GRIB_CHECK(my_grib_set_long(gh, "scaleFactorOfSecondFixedSurface", factor), 0);
      GRIB_CHECK(my_grib_set_long(gh, "scaledValueOfSecondFixedSurface", scaled_level), 0);
    }
}

static
void gribapiDefLevel(int editionNumber, grib_handle *gh, int param, int zaxisID, int levelID, int gcinit)
{
  int lbounds = 0;
  static int warning = 1;
  double dlevel1 = 0, dlevel2 = 0;

  int zaxistype = zaxisInqType(zaxisID);
  int ltype = zaxisInqLtype(zaxisID);
  int ltype2 = zaxisInqLtype2(zaxisID);
  double level = zaxisInqLevel(zaxisID, levelID);

  if ( zaxisInqLbounds(zaxisID, NULL) && zaxisInqUbounds(zaxisID, NULL) )
    {
      lbounds = 1;
      dlevel1 = zaxisInqLbound(zaxisID, levelID);
      dlevel2 = zaxisInqUbound(zaxisID, levelID);
    }
  else
    {
      dlevel1 = level;
      dlevel2 = 0;
    }

  if ( zaxistype == ZAXIS_GENERIC && ltype == 0 )
    {
      Message("Changed zaxis type from %s to %s", zaxisNamePtr(zaxistype), zaxisNamePtr(ZAXIS_PRESSURE));
      zaxistype = ZAXIS_PRESSURE;
      zaxisChangeType(zaxisID, zaxistype);
      zaxisDefUnits(zaxisID, "Pa");
    }

  int grib2ltype = zaxisTypeToGrib2ltype(zaxistype);

  switch (zaxistype)
    {
    case ZAXIS_SURFACE:
    case ZAXIS_MEANSEA:
    case ZAXIS_HEIGHT:
    case ZAXIS_ALTITUDE:
    case ZAXIS_SIGMA:
    case ZAXIS_DEPTH_BELOW_SEA:
    case ZAXIS_ISENTROPIC:
      {
	if ( editionNumber <= 1 )
          {
            gribapiDefLevelType(gh, gcinit, "indicatorOfTypeOfLevel", zaxisTypeToGrib1ltype(zaxistype));
            GRIB_CHECK(my_grib_set_long(gh, "level", (long)level), 0);
          }
        else
          {
            grib2DefLevel(gh, gcinit, grib2ltype, grib2ltype, lbounds, level, dlevel1, dlevel2);
          }

	break;
      }
    case ZAXIS_CLOUD_BASE:
    case ZAXIS_CLOUD_TOP:
    case ZAXIS_ISOTHERM_ZERO:
    case ZAXIS_TOA:
    case ZAXIS_SEA_BOTTOM:
    case ZAXIS_LAKE_BOTTOM:
    case ZAXIS_SEDIMENT_BOTTOM:
    case ZAXIS_SEDIMENT_BOTTOM_TA:
    case ZAXIS_SEDIMENT_BOTTOM_TW:
    case ZAXIS_MIX_LAYER:
    case ZAXIS_ATMOSPHERE:
      {
        if ( editionNumber <= 1 )
          {
            gribapiDefLevelType(gh, gcinit, "indicatorOfTypeOfLevel", zaxisTypeToGrib1ltype(zaxistype));
            if ( lbounds )
              {
                GRIB_CHECK(my_grib_set_long(gh, "topLevel", (long) dlevel1), 0);
                GRIB_CHECK(my_grib_set_long(gh, "bottomLevel", (long) dlevel2), 0);
              }
            else
              {
                GRIB_CHECK(my_grib_set_long(gh, "level", (long) level), 0);
              }
          }
        else
          {
            grib2DefLevel(gh, gcinit, grib2ltype, grib2ltype, lbounds, level, dlevel1, dlevel2);
          }

        break;
      }
    case ZAXIS_HYBRID:
    case ZAXIS_HYBRID_HALF:
      {
        if ( editionNumber <= 1 )
          {
            if ( lbounds )
              {
                gribapiDefLevelType(gh, gcinit, "indicatorOfTypeOfLevel", GRIB1_LTYPE_HYBRID_LAYER);
                GRIB_CHECK(my_grib_set_long(gh, "topLevel", (long) dlevel1), 0);
                GRIB_CHECK(my_grib_set_long(gh, "bottomLevel", (long) dlevel2), 0);
              }
            else
              {
                gribapiDefLevelType(gh, gcinit, "indicatorOfTypeOfLevel", GRIB1_LTYPE_HYBRID);
                GRIB_CHECK(my_grib_set_long(gh, "level", (long) level), 0);
              }
          }
        else
          {
            grib2DefLevel(gh, gcinit, GRIB2_LTYPE_HYBRID, GRIB2_LTYPE_HYBRID, lbounds, level, dlevel1, dlevel2);
          }

        if ( !gcinit )
          {
            int vctsize = zaxisInqVctSize(zaxisID);
            if ( vctsize == 0 && warning )
              {
                char paramstr[32];
                cdiParamToString(param, paramstr, sizeof(paramstr));
                Warning("VCT missing ( param = %s, zaxisID = %d )", paramstr, zaxisID);
                warning = 0;
              }
            GRIB_CHECK(my_grib_set_long(gh, "PVPresent", 1), 0);
            GRIB_CHECK(grib_set_double_array(gh, "pv", zaxisInqVctPtr(zaxisID), (size_t)vctsize), 0);
          }

	break;
      }
    case ZAXIS_PRESSURE:
      {
	double dum;
	char units[128];

	if ( level < 0 ) Warning("Pressure level of %f Pa is below zero!", level);

	zaxisInqUnits(zaxisID, units);
	if ( memcmp(units, "Pa", 2) != 0 )
          {
            level   *= 100;
            dlevel1 *= 100;
            dlevel2 *= 100;
          }

        if ( editionNumber <= 1 )
          {
            long leveltype = GRIB1_LTYPE_ISOBARIC;

            if ( level < 32768 && (level < 100 || modf(level/100, &dum) > 0) )
              leveltype = GRIB1_LTYPE_99;
            else
              level /= 100;

            gribapiDefLevelType(gh, gcinit, "indicatorOfTypeOfLevel", leveltype);
            GRIB_CHECK(my_grib_set_double(gh, "level", level), 0);
	  }
	else
	  {
            if ( ltype2 == -1 ) ltype2 = GRIB2_LTYPE_ISOBARIC;
            grib2DefLevel(gh, gcinit, GRIB2_LTYPE_ISOBARIC, ltype2, lbounds, level, dlevel1, dlevel2);
	  }

	break;
      }
    case ZAXIS_SNOW:
      {
        if ( editionNumber <= 1 )
          ; // not available
	else
          {
            grib2DefLevel(gh, gcinit, GRIB2_LTYPE_SNOW, GRIB2_LTYPE_SNOW, lbounds, level, dlevel1, dlevel2);
          }

	break;
      }
    case ZAXIS_DEPTH_BELOW_LAND:
      {
	char units[128];

	zaxisInqUnits(zaxisID, units);

	if ( editionNumber <= 1 )
	  {
            double scalefactor;
	    if      ( memcmp(units, "mm", 2) == 0 ) scalefactor =   0.1;
	    else if ( memcmp(units, "cm", 2) == 0 ) scalefactor =   1; // cm
	    else if ( memcmp(units, "dm", 2) == 0 ) scalefactor =  10;
	    else                                    scalefactor = 100;

	    gribapiDefLevelType(gh, gcinit, "indicatorOfTypeOfLevel", GRIB1_LTYPE_LANDDEPTH);
	    GRIB_CHECK(my_grib_set_double(gh, "level", level*scalefactor), 0);
	  }
	else
	  {
            double scalefactor;
	    if      ( memcmp(units, "mm", 2) == 0 ) scalefactor = 0.001;
	    else if ( memcmp(units, "cm", 2) == 0 ) scalefactor = 0.01;
	    else if ( memcmp(units, "dm", 2) == 0 ) scalefactor = 0.1;
	    else                                    scalefactor = 1; // meter

            level   *= scalefactor;
            dlevel1 *= scalefactor;
            dlevel2 *= scalefactor;

            grib2DefLevel(gh, gcinit, GRIB2_LTYPE_LANDDEPTH, GRIB2_LTYPE_LANDDEPTH, lbounds, level, dlevel1, dlevel2);
	  }

	break;
      }
    case ZAXIS_REFERENCE:
      {
        unsigned char uuid[CDI_UUID_SIZE];

        if ( !gcinit )
          {
            GRIB_CHECK(my_grib_set_long(gh, "genVertHeightCoords", 1), 0);
          }

        if ( lbounds )
          {
            if ( editionNumber <= 1 )
              ; // not available
            else
              {
                int number = zaxisInqNumber(zaxisID);
                gribapiDefLevelType(gh, gcinit, "typeOfFirstFixedSurface", GRIB2_LTYPE_REFERENCE);
                gribapiDefLevelType(gh, gcinit, "typeOfSecondFixedSurface", GRIB2_LTYPE_REFERENCE);
                GRIB_CHECK(my_grib_set_long(gh, "NV", 6), 0);
                GRIB_CHECK(my_grib_set_long(gh, "nlev", zaxisInqNlevRef(zaxisID)), 0);
                GRIB_CHECK(my_grib_set_long(gh, "numberOfVGridUsed", number), 0);
                size_t len = CDI_UUID_SIZE;
                zaxisInqUUID(zaxisID, uuid);
                if (grib_set_bytes(gh, "uuidOfVGrid", uuid, &len) != 0)
                  {
                    Warning("Can't write UUID!");
                  }
                GRIB_CHECK(my_grib_set_long(gh, "topLevel", (long) dlevel1), 0);
                GRIB_CHECK(my_grib_set_long(gh, "bottomLevel", (long) dlevel2), 0);
              }
          }
        else
          {
            if ( editionNumber <= 1 )
              ; // not available
            else
              {
                int number = zaxisInqNumber(zaxisID);
                gribapiDefLevelType(gh, gcinit, "typeOfFirstFixedSurface", GRIB2_LTYPE_REFERENCE);
                GRIB_CHECK(my_grib_set_long(gh, "NV", 6), 0);
                GRIB_CHECK(my_grib_set_long(gh, "nlev", zaxisInqNlevRef(zaxisID)), 0);
                GRIB_CHECK(my_grib_set_long(gh, "numberOfVGridUsed", number), 0);
                size_t len = CDI_UUID_SIZE;
                zaxisInqUUID(zaxisID, uuid);
                if (grib_set_bytes(gh, "uuidOfVGrid", uuid, &len) != 0)
                  {
                    Warning("Can't write UUID!");
                  }
                GRIB_CHECK(my_grib_set_double(gh, "level", level), 0);
              }
          }

        break;
      }
    case ZAXIS_GENERIC:
      {
	if ( editionNumber <= 1 )
          gribapiDefLevelType(gh, gcinit, "indicatorOfTypeOfLevel", ltype);
        else
          gribapiDefLevelType(gh, gcinit, "typeOfFirstFixedSurface", ltype);

	GRIB_CHECK(my_grib_set_double(gh, "level", level), 0);

	break;
      }
    default:
      {
	Error("Unsupported zaxis type: %s", zaxisNamePtr(zaxistype));
	break;
      }
    }
}
#endif

/* #define GRIBAPIENCODETEST 1 */

#ifdef HAVE_LIBGRIB_API
size_t gribapiEncode(int varID, int levelID, int vlistID, int gridID, int zaxisID,
		     int vdate, int vtime, int tsteptype, int numavg,
		     long datasize, const double *data, int nmiss, unsigned char **gribbuffer, size_t *gribbuffersize,
		     int comptype, void *gribContainer)
{
  size_t recsize = 0;
  void *dummy = NULL;
  int lieee = FALSE;
  /*  int ensID, ensCount, forecast_type; *//* Ensemble Data */
  int typeOfGeneratingProcess;
  int productDefinitionTemplate;
  long bitsPerValue;
  long editionNumber = 2;
  char name[256];
  char stdname[256];
  gribContainer_t *gc = (gribContainer_t *) gribContainer;
  // extern unsigned char _grib_template_GRIB2[];

  int param    = vlistInqVarParam(vlistID, varID);
  int datatype = vlistInqVarDatatype(vlistID, varID);
  typeOfGeneratingProcess = vlistInqVarTypeOfGeneratingProcess(vlistID, varID);
  productDefinitionTemplate = vlistInqVarProductDefinitionTemplate(vlistID, varID);

  vlistInqVarName(vlistID, varID, name);
  vlistInqVarStdname(vlistID, varID, stdname);

#if defined(GRIBAPIENCODETEST)
  grib_handle *gh = (grib_handle *) gribHandleNew(editionNumber);
#else
  grib_handle *gh = gc->gribHandle;
#endif
  GRIB_CHECK(grib_get_long(gh, "editionNumber", &editionNumber), 0);

  if ( editionNumber == 2 )
    {
      if ( typeOfGeneratingProcess == -1 ) typeOfGeneratingProcess = 0;
      if ( ! gc->init ) GRIB_CHECK(my_grib_set_long(gh, "typeOfGeneratingProcess", typeOfGeneratingProcess), 0);
    }

  /*
  if( vlistInqVarEnsemble( vlistID,  varID, &ensID, &ensCount, &forecast_type ) )
    {
      GRIB_CHECK(my_grib_set_long(gh, "typeOfEnsembleForecast", forecast_type ), 0);
      GRIB_CHECK(my_grib_set_long(gh, "numberOfForecastsInEnsemble", ensCount ), 0);
      GRIB_CHECK(my_grib_set_long(gh, "perturbationNumber", ensID ), 0);
    }
  */

  gribapiDefTime((int)editionNumber, productDefinitionTemplate, typeOfGeneratingProcess, gh, vdate, vtime, tsteptype, numavg, vlistInqTaxis(vlistID), gc->init);

  if ( ! gc->init ) gribapiDefInstitut(gh, vlistID, varID);
  if ( ! gc->init ) gribapiDefModel(gh, vlistID, varID);

  if ( ! gc->init ) gribapiDefParam((int)editionNumber, gh, param, name, stdname);

  if ( editionNumber == 2 && (datatype == DATATYPE_FLT32 || datatype == DATATYPE_FLT64) ) lieee = TRUE;

  /* bitsPerValue have to be defined before call to DefGrid (complex packing) */
  //  if ( lieee == FALSE )
    {
      bitsPerValue = grbBitsPerValue(datatype);
      GRIB_CHECK(my_grib_set_long(gh, "bitsPerValue", bitsPerValue), 0);
    }

  gribapiDefGrid((int)editionNumber, gh, gridID, comptype, lieee, datatype, nmiss, gc->init);

  gribapiDefLevel((int)editionNumber, gh, param, zaxisID, levelID, gc->init);

  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  //if (!gc->init)
  {
    int ret = 0;

    /* NOTE: Optional key/value pairs: Note that we do not distinguish
     *       between tiles here! */

    for ( int i=0; i<vlistptr->vars[varID].opt_grib_nentries; i++ )
      {
        if ( vlistptr->vars[varID].opt_grib_kvpair[i].update )
          {
            //DR: Fix for multi-level fields (otherwise only the 1st level is correct)
            if ( zaxisInqSize(zaxisID)==(levelID+1) )
              vlistptr->vars[varID].opt_grib_kvpair[i].update = FALSE;

            if (vlistptr->vars[varID].opt_grib_kvpair[i].data_type == t_double)
              {
                if ( CDI_Debug )
                  Message("key \"%s\"  :   double value = %g\n",
                          vlistptr->vars[varID].opt_grib_kvpair[i].keyword,
                          vlistptr->vars[varID].opt_grib_kvpair[i].dbl_val);
                my_grib_set_double(gh, vlistptr->vars[varID].opt_grib_kvpair[i].keyword,
                                   vlistptr->vars[varID].opt_grib_kvpair[i].dbl_val);
                GRIB_CHECK(ret, 0);
                }
            if (vlistptr->vars[varID].opt_grib_kvpair[i].data_type == t_int)
              {
                if ( CDI_Debug )
                  Message("key \"%s\"  :   integer value = %d\n",
                          vlistptr->vars[varID].opt_grib_kvpair[i].keyword,
                          vlistptr->vars[varID].opt_grib_kvpair[i].int_val);
                my_grib_set_long(gh, vlistptr->vars[varID].opt_grib_kvpair[i].keyword,
                                 (long) vlistptr->vars[varID].opt_grib_kvpair[i].int_val);
                GRIB_CHECK(ret, 0);
              }
          }
      }
  }

  if ( nmiss > 0 )
    {
      GRIB_CHECK(my_grib_set_long(gh, "bitmapPresent", 1), 0);
      GRIB_CHECK(my_grib_set_double(gh, "missingValue", vlistInqVarMissval(vlistID, varID)), 0);
    }

  GRIB_CHECK(grib_set_double_array(gh, "values", data, (size_t)datasize), 0);

  /* get the size of coded message  */
  GRIB_CHECK(grib_get_message(gh, (const void **)&dummy, &recsize), 0);
  recsize += 512; /* add some space for possible filling */
  *gribbuffersize = recsize;
  *gribbuffer = (unsigned char *) malloc(*gribbuffersize);

  /* get a copy of the coded message */
  GRIB_CHECK(grib_get_message_copy(gh, *gribbuffer, &recsize), 0);

#if defined(GRIBAPIENCODETEST)
  gribHandleDelete(gh);
#endif

  gc->init = TRUE;

  return (recsize);
}
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
