#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <stdarg.h>
#include <ctype.h>

#include "cdf.h"
#include "cdi.h"
#include "cdi_int.h"
#include "dmemory.h"
#include "file.h"
#include "gribapi.h"
#ifdef HAVE_LIBNETCDF
#include "stream_cdf.h"
#endif
#include "namespace.h"
#include "resource_handle.h"

#if  defined  (HAVE_LIBCGRIBEX)
#include "cgribex.h"
#endif

extern int cdiPioSerialOpenFileMap(int streamID);

int cdiDefaultCalendar = CALENDAR_PROLEPTIC;

int cdiDefaultInstID   = CDI_UNDEFID;
int cdiDefaultModelID  = CDI_UNDEFID;
int cdiDefaultTableID  = CDI_UNDEFID;
//int cdiNcMissingValue  = CDI_UNDEFID;
int cdiNcChunksizehint = CDI_UNDEFID;
int cdiChunkType       = CHUNK_GRID;
int cdiSplitLtype105   = CDI_UNDEFID;

int cdiIgnoreAttCoordinates = FALSE;
int cdiIgnoreValidRange     = FALSE;
int cdiSkipRecords          = 0;
int cdiInventoryMode        = 1;
size_t CDI_netcdf_hdr_pad   = 0UL;

char *cdiPartabPath   = NULL;
int   cdiPartabIntern = 1;

double cdiDefaultMissval = -9.E33;

const char Filetypes[][9] = {
  "UNKNOWN",
  "GRIB",
  "GRIB2",
  "netCDF",
  "netCDF2",
  "netCDF4",
  "netCDF4c",
  "SERVICE",
  "EXTRA",
  "IEG",
  "HDF5",
};

#undef  UNDEFID
#define UNDEFID  CDI_UNDEFID


int CDI_Debug   = 0;    /* If set to 1, debugging           */
int CDI_Recopt = 0;

int cdiGribApiDebug     = 0;
int cdiDefaultLeveltype = -1;
int cdiDataUnreduced = 0;
int cdiSortName = 0;
int cdiHaveMissval = 0;


static long cdiGetenvInt(char *envName)
{
  char *envString;
  long envValue = -1;
  long fact = 1;

  envString = getenv(envName);

  if ( envString )
    {
      int loop, len;

      len = (int) strlen(envString);
      for ( loop = 0; loop < len; loop++ )
	{
	  if ( ! isdigit((int) envString[loop]) )
	    {
	      switch ( tolower((int) envString[loop]) )
		{
		case 'k':  fact = 1024;        break;
		case 'm':  fact = 1048576;     break;
		case 'g':  fact = 1073741824;  break;
		default:
		  fact = 0;
		  Message("Invalid number string in %s: %s", envName, envString);
		  Warning("%s must comprise only digits [0-9].",envName);
		  break;
		}
	      break;
	    }
	}

      if ( fact ) envValue = fact*atol(envString);

      if ( CDI_Debug ) Message("set %s to %ld", envName, envValue);
    }

  return (envValue);
}

static void
cdiPrintDefaults(void)
{
  fprintf(stderr, "default instID     :  %d\n"
          "default modelID    :  %d\n"
          "default tableID    :  %d\n"
          "default missval    :  %g\n", cdiDefaultInstID,
          cdiDefaultModelID, cdiDefaultTableID, cdiDefaultMissval);
}

void cdiPrintVersion(void)
{
  fprintf(stderr, "     CDI library version : %s\n", cdiLibraryVersion());
#if  defined  (HAVE_LIBCGRIBEX)
  fprintf(stderr, " CGRIBEX library version : %s\n", cgribexLibraryVersion());
#endif
#if  defined  (HAVE_LIBGRIB_API)
  fprintf(stderr, "GRIB_API library version : %s\n", gribapiLibraryVersionString());
#endif
#if  defined  (HAVE_LIBNETCDF)
  fprintf(stderr, "  netCDF library version : %s\n", cdfLibraryVersion());
#endif
#if  defined  (HAVE_LIBHDF5)
  fprintf(stderr, "    HDF5 library version : %s\n", hdfLibraryVersion());
#endif
#if  defined  (HAVE_LIBSERVICE)
  fprintf(stderr, " SERVICE library version : %s\n", srvLibraryVersion());
#endif
#if  defined  (HAVE_LIBEXTRA)
  fprintf(stderr, "   EXTRA library version : %s\n", extLibraryVersion());
#endif
#if  defined  (HAVE_LIBIEG)
  fprintf(stderr, "     IEG library version : %s\n", iegLibraryVersion());
#endif
  fprintf(stderr, "    FILE library version : %s\n", fileLibraryVersion());
}

void cdiDebug(int level)
{
  if ( level == 1 || (level &  2) ) CDI_Debug = 1;

  if ( CDI_Debug ) Message("debug level %d", level);

  if ( level == 1 || (level &  4) ) memDebug(1);

  if ( level == 1 || (level &  8) ) fileDebug(1);

  if ( level == 1 || (level & 16) )
    {
#if  defined  (HAVE_LIBCGRIBEX)
      gribSetDebug(1);
#endif
#if  defined  (HAVE_LIBNETCDF)
      cdfDebug(1);
#endif
#if  defined  (HAVE_LIBSERVICE)
      srvDebug(1);
#endif
#if  defined  (HAVE_LIBEXTRA)
      extDebug(1);
#endif
#if  defined  (HAVE_LIBIEG)
      iegDebug(1);
#endif
    }

  if ( CDI_Debug )
    {
      cdiPrintDefaults();
      cdiPrintDatatypes();
    }
}


int cdiHaveFiletype(int filetype)
{
  int status = 0;

  switch (filetype)
    {
#if  defined  (HAVE_LIBSERVICE)
    case FILETYPE_SRV:  { status = 1; break; }
#endif
#if  defined  (HAVE_LIBEXTRA)
    case FILETYPE_EXT:  { status = 1; break; }
#endif
#if  defined  (HAVE_LIBIEG)
    case FILETYPE_IEG:  { status = 1; break; }
#endif
#if  defined  (HAVE_LIBGRIB)
#if  defined  (HAVE_LIBGRIB_API) || defined  (HAVE_LIBCGRIBEX)
    case FILETYPE_GRB:  { status = 1; break; }
#endif
#if  defined  (HAVE_LIBGRIB_API)
    case FILETYPE_GRB2: { status = 1; break; }
#endif
#endif
#if  defined  (HAVE_LIBNETCDF)
    case FILETYPE_NC:   { status = 1; break; }
#if  defined  (HAVE_NETCDF2)
    case FILETYPE_NC2:  { status = 1; break; }
#endif
#if  defined  (HAVE_NETCDF4)
    case FILETYPE_NC4:  { status = 1; break; }
    case FILETYPE_NC4C: { status = 1; break; }
#endif
#endif
    default: { status = 0; break; }
    }

  return (status);
}

void cdiDefTableID(int tableID)
{
  cdiDefaultTableID = tableID;
  int modelID = cdiDefaultModelID = tableInqModel(tableID);
  cdiDefaultInstID = modelInqInstitut(modelID);
}

static
void cdiSetChunk(const char *chunkAlgo)
{
  //char *pch;
  //size_t len = strlen(chunkAlgo);
  int algo = -1;

  if      ( strcmp("auto",  chunkAlgo)   == 0 ) algo = CHUNK_AUTO;
  else if ( strcmp("grid",  chunkAlgo)   == 0 ) algo = CHUNK_GRID;
  else if ( strcmp("lines", chunkAlgo)   == 0 ) algo = CHUNK_LINES;
  /*
  else if ( (pch = strstr(chunkAlgo,"x")) != 0 )
    {
      int ix, iy;
      ix = atoi(chunkAlgo);
      iy = atoi(pch+1);
      if ( ix > 0 && iy > 0 )
        {
          cdiChunkX = ix;
          cdiChunkY = iy;
          algo = CHUNK_USER;
        }
      else
        Warning("Invalid environment variable CDI_CHUNK_ALGO: %s", chunkAlgo);
    }
  */
  else
    Warning("Invalid environment variable CDI_CHUNK_ALGO: %s", chunkAlgo);

  if ( algo != -1 )
    {
      cdiChunkType = algo;
      if ( CDI_Debug ) Message("set ChunkAlgo to %s", chunkAlgo);
    }
}


void cdiInitialize(void)
{
  static int Init_CDI = FALSE;
  char *envString;
  long value;

  if ( ! Init_CDI )
    {
      Init_CDI = TRUE;

#if  defined  (HAVE_LIBCGRIBEX)
      gribFixZSE(1);   // 1: Fix ZeroShiftError of simple packed spherical harmonics
      gribSetConst(1); // 1: Don't pack constant fields on regular grids
#endif

      value = cdiGetenvInt("CDI_DEBUG");
      if ( value >= 0 ) CDI_Debug = (int) value;

      value = cdiGetenvInt("CDI_GRIBAPI_DEBUG");
      if ( value >= 0 ) cdiGribApiDebug = (int) value;

      value = cdiGetenvInt("CDI_RECOPT");
      if ( value >= 0 ) CDI_Recopt = (int) value;

      value = cdiGetenvInt("CDI_REGULARGRID");
      if ( value >= 0 ) cdiDataUnreduced = (int) value;

      value = cdiGetenvInt("CDI_SORTNAME");
      if ( value >= 0 ) cdiSortName = (int) value;

      value = cdiGetenvInt("CDI_HAVE_MISSVAL");
      if ( value >= 0 ) cdiHaveMissval = (int) value;

      value = cdiGetenvInt("CDI_LEVELTYPE");
      if ( value >= 0 ) cdiDefaultLeveltype = (int) value;

      value = cdiGetenvInt("CDI_NETCDF_HDR_PAD");
      if ( value >= 0 ) CDI_netcdf_hdr_pad = (size_t) value;

      envString = getenv("CDI_MISSVAL");
      if ( envString ) cdiDefaultMissval = atof(envString);
      /*
      envString = getenv("NC_MISSING_VALUE");
      if ( envString ) cdiNcMissingValue = atoi(envString);
      */
      envString = getenv("NC_CHUNKSIZEHINT");
      if ( envString ) cdiNcChunksizehint = atoi(envString);

      envString = getenv("CDI_CHUNK_ALGO");
      if ( envString ) cdiSetChunk(envString);

      envString = getenv("SPLIT_LTYPE_105");
      if ( envString ) cdiSplitLtype105 = atoi(envString);

      envString = getenv("IGNORE_ATT_COORDINATES");
      if ( envString ) cdiIgnoreAttCoordinates = atoi(envString);

      envString = getenv("IGNORE_VALID_RANGE");
      if ( envString ) cdiIgnoreValidRange = atoi(envString);

      envString = getenv("CDI_SKIP_RECORDS");
      if ( envString )
	{
	  cdiSkipRecords = atoi(envString);
	  cdiSkipRecords = cdiSkipRecords > 0 ? cdiSkipRecords : 0;
	}

      envString = getenv("CDI_INVENTORY_MODE");
      if ( envString )
	{
	  if ( strncmp(envString, "time", 4) == 0 )
	    {
	      cdiInventoryMode = 2;
	      if ( CDI_Debug )
		Message("Inventory mode was set to timestep!");
	    }
	}

      envString = getenv("CDI_CALENDAR");
      if ( envString )
	{
	  if      ( strncmp(envString, "standard", 8) == 0 )
	    cdiDefaultCalendar = CALENDAR_STANDARD;
	  else if ( strncmp(envString, "proleptic", 9) == 0 )
	    cdiDefaultCalendar = CALENDAR_PROLEPTIC;
	  else if ( strncmp(envString, "360days", 7) == 0 )
	    cdiDefaultCalendar = CALENDAR_360DAYS;
	  else if ( strncmp(envString, "365days", 7) == 0 )
	    cdiDefaultCalendar = CALENDAR_365DAYS;
	  else if ( strncmp(envString, "366days", 7) == 0 )
	    cdiDefaultCalendar = CALENDAR_366DAYS;
	  else if ( strncmp(envString, "none", 4) == 0 )
	    cdiDefaultCalendar = CALENDAR_NONE;

	  if ( CDI_Debug )
	    Message("Default calendar set to %s!", envString);
	}
#if  defined  (HAVE_LIBCGRIBEX)
      gribSetCalendar(cdiDefaultCalendar);
#endif

      envString = getenv("PARTAB_INTERN");
      if ( envString ) cdiPartabIntern = atoi(envString);

      envString = getenv("PARTAB_PATH");
      if ( envString ) cdiPartabPath = strdup(envString);
    }
}


const char *strfiletype(int filetype)
{
  const char *name;
  int size = (int) (sizeof(Filetypes)/sizeof(char *));

  if ( filetype > 0 && filetype < size )
    name = Filetypes[filetype];
  else
    name = Filetypes[0];

  return (name);
}


void cdiDefGlobal(const char *string, int val)
{
  if      ( strcmp(string, "REGULARGRID")      == 0 ) cdiDataUnreduced = val;
  else if ( strcmp(string, "GRIBAPI_DEBUG")    == 0 ) cdiGribApiDebug = val;
  else if ( strcmp(string, "SORTNAME")         == 0 ) cdiSortName = val;
  else if ( strcmp(string, "HAVE_MISSVAL")     == 0 ) cdiHaveMissval = val;
  else if ( strcmp(string, "NC_CHUNKSIZEHINT") == 0 ) cdiNcChunksizehint = val;
  else if ( strcmp(string, "NETCDF_HDR_PAD")   == 0 ) CDI_netcdf_hdr_pad = (size_t) val;
  else Warning("Unsupported global key: %s", string);
}


void cdiDefMissval(double missval)
{
  cdiInitialize();

  cdiDefaultMissval = missval;
}


double cdiInqMissval(void)
{
  cdiInitialize();

  return (cdiDefaultMissval);
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

