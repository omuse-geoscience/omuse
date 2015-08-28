#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 600
#endif

#include <ctype.h>
#include <stdio.h>
#include <string.h>

#include "cdi.h"
#include "cdi_int.h"
#include "cdi_cksum.h"
#include "cdf.h"
#include "dmemory.h"
#include "error.h"
#include "stream_grb.h"
#include "stream_cdf.h"
#include "stream_srv.h"
#include "stream_ext.h"
#include "stream_ieg.h"
#include "file.h"
#include "cgribex.h"
#include "gribapi.h"
#include "cdf.h"
#include "service.h"
#include "extra.h"
#include "ieg.h"
#include "vlist.h"
#include "serialize.h"
#include "resource_handle.h"
#include "resource_unpack.h"

#include "namespace.h"


static stream_t *stream_new_entry(int resH);
static void stream_delete_entry(stream_t *streamptr);
static int streamCompareP(void * streamptr1, void * streamptr2);
static void streamDestroyP(void * streamptr);
static void streamPrintP(void * streamptr, FILE * fp);
static int streamGetPackSize(void * streamptr, void *context);
static void streamPack(void * streamptr, void * buff, int size, int * position, void *context);
static int streamTxCode(void);

const resOps streamOps = {
  streamCompareP,
  streamDestroyP,
  streamPrintP,
  streamGetPackSize,
  streamPack,
  streamTxCode
};




#undef  IsBigendian
#define IsBigendian()  ( u_byteorder.c[sizeof(long) - 1] )


static
int getByteorder(int byteswap)
{
  static union {unsigned long l; unsigned char c[sizeof(long)];} u_byteorder = {1};
  int byteorder = -1;

  if ( IsBigendian() )
    {
      if ( byteswap ) byteorder = CDI_LITTLEENDIAN;
      else            byteorder = CDI_BIGENDIAN;
    }
  else
    {
      if ( byteswap ) byteorder = CDI_BIGENDIAN;
      else            byteorder = CDI_LITTLEENDIAN;
    }

  return (byteorder);
}

// used also in CDO
int cdiGetFiletype(const char *filename, int *byteorder)
{
  int filetype = CDI_EUFTYPE;
  int swap = 0;
  int version;
  long recpos;
  char buffer[8];

  int fileID = fileOpen(filename, "r");

  if ( fileID == CDI_UNDEFID )
    {
      if ( strncmp(filename, "http:", 5) == 0 || strncmp(filename, "https:", 6) == 0 )
	return (FILETYPE_NC);
      else
	return (CDI_ESYSTEM);
    }

  if ( fileRead(fileID, buffer, 8) != 8 ) return (CDI_EUFTYPE);

  fileRewind(fileID);

  if ( memcmp(buffer, "GRIB", 4) == 0 )
    {
      version = buffer[7];
      if ( version <= 1 )
	{
	  filetype = FILETYPE_GRB;
	  if ( CDI_Debug ) Message("found GRIB file = %s, version %d", filename, version);
	}
      else if ( version == 2 )
	{
	  filetype = FILETYPE_GRB2;
	  if ( CDI_Debug ) Message("found GRIB2 file = %s", filename);
	}
    }
  else if ( memcmp(buffer, "CDF\001", 4) == 0 )
    {
      filetype = FILETYPE_NC;
      if ( CDI_Debug ) Message("found CDF1 file = %s", filename);
    }
  else if ( memcmp(buffer, "CDF\002", 4) == 0 )
    {
      filetype = FILETYPE_NC2;
      if ( CDI_Debug ) Message("found CDF2 file = %s", filename);
    }
  else if ( memcmp(buffer+1, "HDF", 3) == 0 )
    {
      filetype = FILETYPE_NC4;
      if ( CDI_Debug ) Message("found HDF file = %s", filename);
    }
#if  defined  (HAVE_LIBSERVICE)
  else if ( srvCheckFiletype(fileID, &swap) )
    {
      filetype = FILETYPE_SRV;
      if ( CDI_Debug ) Message("found SRV file = %s", filename);
    }
#endif
#if  defined  (HAVE_LIBEXTRA)
  else if ( extCheckFiletype(fileID, &swap) )
    {
      filetype = FILETYPE_EXT;
      if ( CDI_Debug ) Message("found EXT file = %s", filename);
    }
#endif
#if  defined  (HAVE_LIBIEG)
  else if ( iegCheckFiletype(fileID, &swap) )
    {
      filetype = FILETYPE_IEG;
      if ( CDI_Debug ) Message("found IEG file = %s", filename);
    }
#endif
  else if ( gribCheckSeek(fileID, &recpos, &version) == 0 )
    {
      if ( version <= 1 )
	{
	  filetype = FILETYPE_GRB;
	  if ( CDI_Debug ) Message("found seeked GRIB file = %s", filename);
	}
      else if ( version == 2 )
	{
	  filetype = FILETYPE_GRB2;
	  if ( CDI_Debug ) Message("found seeked GRIB2 file = %s", filename);
	}
    }

  fileClose(fileID);

  *byteorder = getByteorder(swap);

  return (filetype);
}

/*
@Function  streamInqFiletype
@Title     Get the filetype

@Prototype int streamInqFiletype(int streamID)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenRead} or @fref{streamOpenWrite}.

@Description
The function @func{streamInqFiletype} returns the filetype of a stream.

@Result
@func{streamInqFiletype} returns the type of the file format,
one of the set of predefined CDI file format types.
The valid CDI file format types are @func{FILETYPE_GRB}, @func{FILETYPE_GRB2}, @func{FILETYPE_NC}, @func{FILETYPE_NC2},
@func{FILETYPE_NC4}, @func{FILETYPE_NC4C}, @func{FILETYPE_SRV}, @func{FILETYPE_EXT} and @func{FILETYPE_IEG}.

@EndFunction
*/
int streamInqFiletype(int streamID)
{
  stream_t *streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  return (streamptr->filetype);
}


int getByteswap(int byteorder)
{
  static union {unsigned long l; unsigned char c[sizeof(long)];} u_byteorder = {1};
  int byteswap = 0;

  if ( IsBigendian() )
    {
      if ( byteorder == CDI_LITTLEENDIAN ) byteswap = TRUE;
    }
  else
    {
      if ( byteorder == CDI_BIGENDIAN ) byteswap = TRUE;
    }

  return (byteswap);
}

/*
@Function  streamDefByteorder
@Title     Define the byte order

@Prototype void streamDefByteorder(int streamID, int byteorder)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenWrite}.
    @Item  byteorder The byte order of a dataset, one of the CDI constants @func{CDI_BIGENDIAN} and
                     @func{CDI_LITTLEENDIAN}.

@Description
The function @func{streamDefByteorder} defines the byte order of a binary dataset
with the file format type @func{FILETYPE_SRV}, @func{FILETYPE_EXT} or @func{FILETYPE_IEG}.

@EndFunction
*/
void streamDefByteorder(int streamID, int byteorder)
{
  int filetype;
  stream_t *streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  streamptr->byteorder = byteorder;
  filetype = streamptr->filetype;

  switch (filetype)
    {
#if  defined  (HAVE_LIBSERVICE)
    case FILETYPE_SRV:
      {
	srvrec_t *srvp = (srvrec_t*) streamptr->record->exsep;
	srvp->byteswap = getByteswap(byteorder);

	break;
      }
#endif
#if  defined  (HAVE_LIBEXTRA)
    case FILETYPE_EXT:
      {
	extrec_t *extp = (extrec_t*) streamptr->record->exsep;
	extp->byteswap = getByteswap(byteorder);

	break;
      }
#endif
#if  defined  (HAVE_LIBIEG)
    case FILETYPE_IEG:
      {
	iegrec_t *iegp = (iegrec_t*) streamptr->record->exsep;
	iegp->byteswap = getByteswap(byteorder);

	break;
      }
#endif
    }
  reshSetStatus(streamID, &streamOps, RESH_DESYNC_IN_USE);
}

/*
@Function  streamInqByteorder
@Title     Get the byte order

@Prototype int streamInqByteorder(int streamID)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenRead} or @fref{streamOpenWrite}.

@Description
The function @func{streamInqByteorder} returns the byte order of a binary dataset
with the file format type @func{FILETYPE_SRV}, @func{FILETYPE_EXT} or @func{FILETYPE_IEG}.

@Result
@func{streamInqByteorder} returns the type of the byte order.
The valid CDI byte order types are @func{CDI_BIGENDIAN} and @func{CDI_LITTLEENDIAN}

@EndFunction
*/
int streamInqByteorder(int streamID)
{
  stream_t *streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  return (streamptr->byteorder);
}


const char *streamFilesuffix(int filetype)
{
  // static char *fileSuffix[] = {"", ".grb", ".g2", ".nc", ".nc", ".nc4", ".nc4", ".srv", ".ext", ".ieg"};
  /* note: the 2nd dimenstion of the fileSuffix array must be equal to or
   * larger than the length of the longest suffix (dot and \0 terminator
   * included) */
  static const char fileSuffix[][5] = {"", ".grb", ".grb", ".nc", ".nc", ".nc", ".nc", ".srv", ".ext", ".ieg"};
  int size = (int)(sizeof(fileSuffix)/sizeof(fileSuffix[0]));

  if ( filetype > 0 && filetype < size )
    return (fileSuffix[filetype]);
  else
    return (fileSuffix[0]);
}


const char *streamFilename(int streamID)
{
  stream_t *streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  return (streamptr->filename);
}

static
long cdiInqTimeSize(int streamID)
{
  int tsID = 0, nrecs;
  stream_t *streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  long ntsteps = streamptr->ntsteps;

  if ( ntsteps == (long)CDI_UNDEFID )
    while ( (nrecs = streamInqTimestep(streamID, tsID++)) )
      ntsteps = streamptr->ntsteps;

  return (ntsteps);
}

static
int cdiInqContents(stream_t * streamptr)
{
  int status = 0;

  int filetype = streamptr->filetype;

  switch (filetype)
    {
#if  defined  (HAVE_LIBGRIB)
    case FILETYPE_GRB:
    case FILETYPE_GRB2:
      {
        status = grbInqContents(streamptr);
	break;
      }
#endif
#if  defined  (HAVE_LIBSERVICE)
    case FILETYPE_SRV:
      {
        status = srvInqContents(streamptr);
	break;
      }
#endif
#if  defined  (HAVE_LIBEXTRA)
    case FILETYPE_EXT:
      {
        status = extInqContents(streamptr);
	break;
      }
#endif
#if  defined  (HAVE_LIBIEG)
    case FILETYPE_IEG:
      {
        status = iegInqContents(streamptr);
	break;
      }
#endif
#if  defined  (HAVE_LIBNETCDF)
    case FILETYPE_NC:
    case FILETYPE_NC2:
    case FILETYPE_NC4:
    case FILETYPE_NC4C:
      {
        status = cdfInqContents(streamptr);
	break;
      }
#endif
    default:
      {
	if ( CDI_Debug )
	  Message("%s support not compiled in!", strfiletype(filetype));

	status = CDI_ELIBNAVAIL;
        break;
      }
    }

  if ( status == 0 )
    {
      int vlistID = streamptr->vlistID;
      int taxisID = vlistInqTaxis(vlistID);
      if ( taxisID != CDI_UNDEFID )
        {
          taxis_t *taxisptr1 = &streamptr->tsteps[0].taxis;
          taxis_t *taxisptr2 = taxisPtr(taxisID);
          ptaxisCopy(taxisptr2, taxisptr1);
        }
    }

  return (status);
}

int cdiStreamOpenDefaultDelegate(const char *filename, const char *filemode,
                                 int filetype, stream_t *streamptr,
                                 int recordBufIsToBeCreated)
{
  int fileID;
  switch (filetype)
    {
#if  defined  (HAVE_LIBGRIB)
    case FILETYPE_GRB:
    case FILETYPE_GRB2:
      {
        fileID = gribOpen(filename, filemode);
        if ( fileID < 0 ) fileID = CDI_ESYSTEM;
        if (recordBufIsToBeCreated)
          {
            streamptr->record = (Record *) malloc(sizeof(Record));
            streamptr->record->buffer = NULL;
          }
        break;
      }
#endif
#if  defined  (HAVE_LIBSERVICE)
    case FILETYPE_SRV:
      {
        fileID = fileOpen(filename, filemode);
        if ( fileID < 0 ) fileID = CDI_ESYSTEM;
        if (recordBufIsToBeCreated)
          {
            streamptr->record = (Record *) malloc(sizeof(Record));
            streamptr->record->buffer = NULL;
            streamptr->record->exsep  = srvNew();
          }
        break;
      }
#endif
#if  defined  (HAVE_LIBEXTRA)
    case FILETYPE_EXT:
      {
        fileID = fileOpen(filename, filemode);
        if ( fileID < 0 ) fileID = CDI_ESYSTEM;
        if (recordBufIsToBeCreated)
          {
            streamptr->record = (Record *) malloc(sizeof(Record));
            streamptr->record->buffer = NULL;
            streamptr->record->exsep  = extNew();
          }
        break;
      }
#endif
#if  defined  (HAVE_LIBIEG)
    case FILETYPE_IEG:
      {
        fileID = fileOpen(filename, filemode);
        if ( fileID < 0 ) fileID = CDI_ESYSTEM;
        if (recordBufIsToBeCreated)
          {
            streamptr->record = (Record *) malloc(sizeof(Record));
            streamptr->record->buffer = NULL;
            streamptr->record->exsep   = iegNew();
          }
        break;
      }
#endif
#if  defined  (HAVE_LIBNETCDF)
    case FILETYPE_NC:
      {
        fileID = cdfOpen(filename, filemode);
        break;
      }
    case FILETYPE_NC2:
      {
        fileID = cdfOpen64(filename, filemode);
        break;
      }
    case FILETYPE_NC4:
    case FILETYPE_NC4C:
      {
        fileID = cdf4Open(filename, filemode, &filetype);
        break;
      }
#endif
    default:
      {
        if ( CDI_Debug ) Message("%s support not compiled in!", strfiletype(filetype));
        return (CDI_ELIBNAVAIL);
      }
    }

  streamptr->filetype = filetype;

  return fileID;
}


static int
streamOpenID(const char *filename, const char *filemode, int filetype,
             int resH)
{
  int fileID = CDI_UNDEFID;
  int status;

  if ( CDI_Debug )
    Message("Open %s mode %c file %s", strfiletype(filetype), (int) *filemode,
            filename?filename:"(NUL)");

  if ( ! filename || ! filemode || filetype < 0 ) return (CDI_EINVAL);

  stream_t *streamptr = stream_new_entry(resH);
  int streamID = CDI_ESYSTEM;

  {
    int (*streamOpenDelegate)(const char *filename, const char *filemode,
                              int filetype, stream_t *streamptr, int recordBufIsToBeCreated)
      = (int (*)(const char *, const char *, int, stream_t *, int))
      namespaceSwitchGet(NSSWITCH_STREAM_OPEN_BACKEND).func;

    fileID = streamOpenDelegate(filename, filemode, filetype, streamptr, 1);
  }

  if (fileID < 0)
    {
      free(streamptr->record);
      stream_delete_entry(streamptr);
      streamID = fileID;
    }
  else
    {
      streamID  = streamptr->self;

      if ( streamID < 0 ) return (CDI_ELIMIT);

      streamptr->filemode = tolower(*filemode);
      streamptr->filename = strdupx(filename);
      streamptr->fileID   = fileID;

      if ( streamptr->filemode == 'r' )
	{
	  int vlistID = vlistCreate();
	  if ( vlistID < 0 ) return(CDI_ELIMIT);

	  streamptr->vlistID = vlistID;
	  /* cdiReadByteorder(streamID); */
	  status = cdiInqContents(streamptr);
	  if ( status < 0 ) return (status);
	  vlist_t *vlistptr = vlist_to_pointer(streamptr->vlistID);
	  vlistptr->ntsteps = streamptr->ntsteps;
	}
    }

  return (streamID);
}

int streamOpen(const char *filename, const char *filemode, int filetype)
{
  return streamOpenID(filename, filemode, filetype, CDI_UNDEFID);
}

static int streamOpenA(const char *filename, const char *filemode, int filetype)
{
  int fileID = CDI_UNDEFID;
  int streamID = CDI_ESYSTEM;
  int status;
  stream_t *streamptr = stream_new_entry(CDI_UNDEFID);
  vlist_t *vlistptr;

  if ( CDI_Debug )
    Message("Open %s file (mode=%c); filename: %s", strfiletype(filetype), (int) *filemode, filename);
  if ( CDI_Debug ) printf("streamOpenA: %s\n", filename); // seg fault without this line on thunder/squall with "cdo cat x y"

  if ( ! filename || ! filemode || filetype < 0 ) return (CDI_EINVAL);

  {
    int (*streamOpenDelegate)(const char *filename, const char *filemode,
                              int filetype, stream_t *streamptr, int recordBufIsToBeCreated)
      = (int (*)(const char *, const char *, int, stream_t *, int))
      namespaceSwitchGet(NSSWITCH_STREAM_OPEN_BACKEND).func;

    fileID = streamOpenDelegate(filename, "r", filetype, streamptr, 1);
  }

  if ( fileID == CDI_UNDEFID || fileID == CDI_ELIBNAVAIL || fileID == CDI_ESYSTEM ) return (fileID);

  streamID = streamptr->self;

  streamptr->filemode = tolower(*filemode);
  streamptr->filename = strdupx(filename);
  streamptr->fileID   = fileID;

  streamptr->vlistID = vlistCreate();
  /* cdiReadByteorder(streamID); */
  status = cdiInqContents(streamptr);
  if ( status < 0 ) return (status);
  vlistptr = vlist_to_pointer(streamptr->vlistID);
  vlistptr->ntsteps = (int)cdiInqTimeSize(streamID);

  {
    void (*streamCloseDelegate)(stream_t *streamptr, int recordBufIsToBeDeleted)
      = (void (*)(stream_t *, int))
      namespaceSwitchGet(NSSWITCH_STREAM_CLOSE_BACKEND).func;

    streamCloseDelegate(streamptr, 0);
  }

  switch (filetype)
    {
#if  defined  (HAVE_LIBGRIB)
    case FILETYPE_GRB:
    case FILETYPE_GRB2:
      {
        fileID = gribOpen(filename, filemode);
	break;
      }
#endif
#if  defined  (HAVE_LIBSERVICE)
    case FILETYPE_SRV:
      {
        fileID = fileOpen(filename, filemode);
	break;
      }
#endif
#if  defined  (HAVE_LIBEXTRA)
    case FILETYPE_EXT:
      {
        fileID = fileOpen(filename, filemode);
	break;
      }
#endif
#if  defined  (HAVE_LIBIEG)
    case FILETYPE_IEG:
      {
        fileID = fileOpen(filename, filemode);
	break;
      }
#endif
#if  defined  (HAVE_LIBNETCDF)
    case FILETYPE_NC:
      {
	fileID = cdfOpen(filename, filemode);
	streamptr->ncmode = 2;
	break;
      }
    case FILETYPE_NC2:
      {
	fileID = cdfOpen64(filename, filemode);
	streamptr->ncmode = 2;
	break;
      }
    case FILETYPE_NC4:
    case FILETYPE_NC4C:
      {
	fileID = cdf4Open(filename, filemode, &filetype);
	streamptr->ncmode = 2;
	break;
      }
#endif
    default:
      {
	if ( CDI_Debug ) Message("%s support not compiled in!", strfiletype(filetype));
	return (CDI_ELIBNAVAIL);
      }
    }

  if ( fileID == CDI_UNDEFID )
    streamID = CDI_UNDEFID;
  else
    streamptr->fileID = fileID;

  return (streamID);
}

/*
@Function  streamOpenRead
@Title     Open a dataset for reading

@Prototype int streamOpenRead(const char *path)
@Parameter
    @Item  path  The name of the dataset to be read.

@Description
The function @func{streamOpenRead} opens an existing dataset for reading.

@Result
Upon successful completion @func{streamOpenRead} returns an identifier to the
open stream. Otherwise, a negative number with the error status is returned.

@Errors
@List
   @Item  CDI_ESYSTEM     Operating system error.
   @Item  CDI_EINVAL      Invalid argument.
   @Item  CDI_EUFILETYPE  Unsupported file type.
   @Item  CDI_ELIBNAVAIL  Library support not compiled in.
@EndList

@Example
Here is an example using @func{streamOpenRead} to open an existing netCDF
file named @func{foo.nc} for reading:

@Source
#include "cdi.h"
   ...
int streamID;
   ...
streamID = streamOpenRead("foo.nc");
if ( streamID < 0 ) handle_error(streamID);
   ...
@EndSource
@EndFunction
*/
int streamOpenRead(const char *filename)
{
  cdiInitialize();

  int byteorder = 0;
  int filetype = cdiGetFiletype(filename, &byteorder);

  if ( filetype < 0 ) return (filetype);

  int streamID = streamOpen(filename, "r", filetype);

  if ( streamID >= 0 )
    {
      stream_t *streamptr = stream_to_pointer(streamID);
      streamptr->byteorder = byteorder;
    }

  return (streamID);
}


int streamOpenAppend(const char *filename)
{
  cdiInitialize();

  int byteorder = 0;
  int filetype = cdiGetFiletype(filename, &byteorder);

  if ( filetype < 0 ) return (filetype);

  int streamID = streamOpenA(filename, "a", filetype);

  if ( streamID >= 0 )
    {
      stream_t *streamptr = stream_to_pointer(streamID);
      streamptr->byteorder = byteorder;
    }

  return (streamID);
}

/*
@Function  streamOpenWrite
@Title     Create a new dataset

@Prototype int streamOpenWrite(const char *path, int filetype)
@Parameter
    @Item  path      The name of the new dataset.
    @Item  filetype  The type of the file format, one of the set of predefined CDI file format types.
                     The valid CDI file format types are @func{FILETYPE_GRB}, @func{FILETYPE_GRB2}, @func{FILETYPE_NC},
                     @func{FILETYPE_NC2}, @func{FILETYPE_NC4}, @func{FILETYPE_NC4C}, @func{FILETYPE_SRV},
                     @func{FILETYPE_EXT} and @func{FILETYPE_IEG}.

@Description
The function @func{streamOpenWrite} creates a new datset.
@Result
Upon successful completion @func{streamOpenWrite} returns an identifier to the
open stream. Otherwise, a negative number with the error status is returned.

@Errors
@List
   @Item  CDI_ESYSTEM     Operating system error.
   @Item  CDI_EINVAL      Invalid argument.
   @Item  CDI_EUFILETYPE  Unsupported file type.
   @Item  CDI_ELIBNAVAIL  Library support not compiled in.
@EndList

@Example
Here is an example using @func{streamOpenWrite} to create a new netCDF file named @func{foo.nc} for writing:

@Source
#include "cdi.h"
   ...
int streamID;
   ...
streamID = streamOpenWrite("foo.nc", FILETYPE_NC);
if ( streamID < 0 ) handle_error(streamID);
   ...
@EndSource
@EndFunction
*/
int streamOpenWrite(const char *filename, int filetype)
{
  cdiInitialize();

  return (streamOpen(filename, "w", filetype));
}

static
void streamDefaultValue ( stream_t * streamptr )
{
  int i;

  streamptr->self              = CDI_UNDEFID;
  streamptr->accesstype        = CDI_UNDEFID;
  streamptr->accessmode        = 0;
  streamptr->filetype          = FILETYPE_UNDEF;
  streamptr->byteorder         = CDI_UNDEFID;
  streamptr->fileID            = 0;
  streamptr->filemode          = 0;
  streamptr->numvals           = 0;
  streamptr->filename          = NULL;
  streamptr->record            = NULL;
  streamptr->varsAllocated     = 0;
  streamptr->nrecs             = 0;
  streamptr->nvars             = 0;
  streamptr->vars              = NULL;
  streamptr->ncmode            = 0;
  streamptr->curTsID           = CDI_UNDEFID;
  streamptr->rtsteps           = 0;
  streamptr->ntsteps           = CDI_UNDEFID;
  streamptr->tsteps            = NULL;
  streamptr->tstepsTableSize   = 0;
  streamptr->tstepsNextID      = 0;
  streamptr->historyID         = CDI_UNDEFID;
  streamptr->vlistID           = CDI_UNDEFID;
  streamptr->globalatts        = 0;
  streamptr->localatts         = 0;
  streamptr->vct.ilev          = 0;
  streamptr->vct.mlev          = 0;
  streamptr->vct.ilevID        = CDI_UNDEFID;
  streamptr->vct.mlevID        = CDI_UNDEFID;
  streamptr->unreduced         = cdiDataUnreduced;
  streamptr->sortname          = cdiSortName;
  streamptr->have_missval      = cdiHaveMissval;
  streamptr->comptype          = COMPRESS_NONE;
  streamptr->complevel         = 0;

  basetimeInit(&streamptr->basetime);

  for ( i = 0; i < MAX_GRIDS_PS; i++ ) streamptr->xdimID[i]   = CDI_UNDEFID;
  for ( i = 0; i < MAX_GRIDS_PS; i++ ) streamptr->ydimID[i]   = CDI_UNDEFID;
  for ( i = 0; i < MAX_ZAXES_PS; i++ ) streamptr->zaxisID[i]  = CDI_UNDEFID;
  for ( i = 0; i < MAX_GRIDS_PS; i++ ) streamptr->ncxvarID[i] = CDI_UNDEFID;
  for ( i = 0; i < MAX_GRIDS_PS; i++ ) streamptr->ncyvarID[i] = CDI_UNDEFID;
  for ( i = 0; i < MAX_GRIDS_PS; i++ ) streamptr->ncavarID[i] = CDI_UNDEFID;

  streamptr->gribContainers    = NULL;
  streamptr->vlistIDorig       = CDI_UNDEFID;
}


static stream_t *stream_new_entry(int resH)
{
  stream_t *streamptr;

  cdiInitialize(); /* ***************** make MT version !!! */

  streamptr = (stream_t *) xmalloc(sizeof(stream_t));
  streamDefaultValue ( streamptr );
  if (resH == CDI_UNDEFID)
    streamptr->self = reshPut(streamptr, &streamOps);
  else
    {
      streamptr->self = resH;
      reshReplace(resH, streamptr, &streamOps);
    }

  return streamptr;
}


void
cdiStreamCloseDefaultDelegate(stream_t *streamptr, int recordBufIsToBeDeleted)
{
  int fileID   = streamptr->fileID;
  int filetype = streamptr->filetype;
  if ( fileID == CDI_UNDEFID )
    Warning("File %s not open!", streamptr->filename);
  else
    switch (filetype)
      {
#if  defined  (HAVE_LIBGRIB)
      case FILETYPE_GRB:
      case FILETYPE_GRB2:
        {
          gribClose(fileID);
          if (recordBufIsToBeDeleted)
            gribContainersDelete(streamptr);
          break;
        }
#endif
#if  defined  (HAVE_LIBSERVICE)
      case FILETYPE_SRV:
        {
          fileClose(fileID);
          if (recordBufIsToBeDeleted)
            srvDelete(streamptr->record->exsep);
          break;
        }
#endif
#if  defined  (HAVE_LIBEXTRA)
      case FILETYPE_EXT:
        {
          fileClose(fileID);
          if (recordBufIsToBeDeleted)
            extDelete(streamptr->record->exsep);
          break;
        }
#endif
#if  defined  (HAVE_LIBIEG)
      case FILETYPE_IEG:
        {
          fileClose(fileID);
          if (recordBufIsToBeDeleted)
            iegDelete(streamptr->record->exsep);
          break;
        }
#endif
#if  defined  (HAVE_LIBNETCDF)
      case FILETYPE_NC:
      case FILETYPE_NC2:
      case FILETYPE_NC4:
      case FILETYPE_NC4C:
        {
          cdfClose(fileID);
          break;
        }
#endif
      default:
        {
          Error("%s support not compiled in (fileID = %d)!", strfiletype(filetype), fileID);
          break;
        }
      }
}


static void deallocate_sleveltable_t(sleveltable_t *entry)
{
  if (entry->recordID) free(entry->recordID);
  if (entry->lindex)   free(entry->lindex);
  entry->recordID = NULL;
  entry->lindex   = NULL;
}


/*
@Function  streamClose
@Title     Close an open dataset

@Prototype  void streamClose(int streamID)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenRead} or @fref{streamOpenWrite}.

@Description
The function @func{streamClose} closes an open dataset.

@EndFunction
*/
void streamClose(int streamID)
{
  int index;
  int vlistID;
  stream_t *streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  if ( CDI_Debug )
    Message("streamID = %d filename = %s", streamID, streamptr->filename);

  vlistID  = streamptr->vlistID;

  void (*streamCloseDelegate)(stream_t *streamptr, int recordBufIsToBeDeleted)
    = (void (*)(stream_t *, int))
    namespaceSwitchGet(NSSWITCH_STREAM_CLOSE_BACKEND).func;

  if ( streamptr->filetype != -1 ) streamCloseDelegate(streamptr, 1);

  if ( streamptr->record )
    {
      if ( streamptr->record->buffer )
        free(streamptr->record->buffer);

      free(streamptr->record);
    }

  streamptr->filetype = 0;
  if ( streamptr->filename ) free(streamptr->filename);

  for ( index = 0; index < streamptr->nvars; index++ )
    {
      sleveltable_t *pslev = streamptr->vars[index].recordTable;
      unsigned nsub = streamptr->vars[index].subtypeSize >= 0
        ? (unsigned)streamptr->vars[index].subtypeSize : 0U;
      for (size_t isub=0; isub < nsub; isub++)
        {
          deallocate_sleveltable_t(pslev + isub);
        }
      if (pslev) free(pslev);
    }
  free(streamptr->vars);
  streamptr->vars = NULL;

  for ( index = 0; index < streamptr->ntsteps; ++index )
    {
      if ( streamptr->tsteps[index].records )
	free(streamptr->tsteps[index].records);
      if ( streamptr->tsteps[index].recIDs )
	free(streamptr->tsteps[index].recIDs);
      taxisDestroyKernel(&streamptr->tsteps[index].taxis);
    }

  if ( streamptr->tsteps ) free(streamptr->tsteps);

  if ( streamptr->basetime.timevar_cache ) free(streamptr->basetime.timevar_cache);

  if ( vlistID != -1 )
    {
      if ( streamptr->filemode != 'w' )
	if ( vlistInqTaxis(vlistID) != -1 )
	  {
	    taxisDestroy(vlistInqTaxis(vlistID));
	  }

      vlistDestroy(vlistID);
      /* decrease lock counter of the original vlist by 1 */
      if ( streamptr->vlistIDorig != CDI_UNDEFID ) {
        /* Note: Here we have to check if the original vlist still
         * exists. If, for example, the garbage collection routine
         * reshListDestruct takes care of the destruction of objects,
         * then the original vlist might have been destroyed
         * beforehand. */
        if (reshExists(streamptr->vlistIDorig, &vlistOps) != 0)
          vlist_unlock(streamptr->vlistIDorig);
      }
    }

  stream_delete_entry(streamptr);
}

static void stream_delete_entry(stream_t *streamptr)
{
  int idx;

  xassert ( streamptr );

  idx = streamptr->self;
  free ( streamptr );
  reshRemove ( idx, &streamOps );

  if ( CDI_Debug )
    Message("Removed idx %d from stream list", idx);
}


void cdiStreamSync_(stream_t *streamptr)
{
  int fileID   = streamptr->fileID;
  int filetype = streamptr->filetype;
  int vlistID  = streamptr->vlistID;
  int nvars    = vlistNvars(vlistID);

  if ( fileID == CDI_UNDEFID )
    Warning("File %s not open!", streamptr->filename);
  else if ( vlistID == CDI_UNDEFID )
    Warning("Vlist undefined for file %s!", streamptr->filename);
  else if ( nvars == 0 )
    Warning("No variables defined!");
  else
    {
      if ( streamptr->filemode == 'w' || streamptr->filemode == 'a' )
	{
	  switch (filetype)
	    {
#if  defined  (HAVE_LIBNETCDF)
	    case FILETYPE_NC:
	    case FILETYPE_NC2:
	    case FILETYPE_NC4:
	    case FILETYPE_NC4C:
	      {
		void cdf_sync(int ncid);
		if ( streamptr->ncmode == 2 ) cdf_sync(fileID);
		break;
	      }
#endif
	    default:
	      {
		fileFlush(fileID);
		break;
	      }
	    }
	}
    }
}

/*
@Function  streamSync
@Title     Synchronize an Open Dataset to Disk

@Prototype  void streamSync(int streamID)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenWrite}.

@Description
The function @func{streamSync} offers a way to synchronize the disk copy of a dataset with in-memory buffers.

@EndFunction
*/
void streamSync(int streamID)
{
  stream_t *streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  void (*myStreamSync_)(stream_t *streamptr)
    = (void (*)(stream_t *))namespaceSwitchGet(NSSWITCH_STREAM_SYNC).func;
  myStreamSync_(streamptr);
}


int cdiStreamDefTimestep_(stream_t *streamptr, int tsID)
{
  int taxisID;

  if ( CDI_Debug )
    Message("streamID = %d  tsID = %d", streamptr->self, tsID);

  stream_check_ptr(__func__, streamptr);

  int vlistID = streamptr->vlistID;

  int time_is_varying = vlistHasTime(vlistID);

  if ( time_is_varying )
    {
      taxisID = vlistInqTaxis(vlistID);
      if ( taxisID == CDI_UNDEFID )
        {
          Warning("taxisID undefined for fileID = %d! Using absolute time axis.", streamptr->self);
          taxisID = taxisCreate(TAXIS_ABSOLUTE);
          vlistDefTaxis(vlistID, taxisID);
        }
    }

  int newtsID = tstepsNewEntry(streamptr);

  if ( tsID != newtsID )
    Error("Internal problem: tsID = %d newtsID = %d", tsID, newtsID);

  streamptr->curTsID = tsID;

  if ( time_is_varying )
    {
      taxis_t *taxisptr1 = taxisPtr(taxisID);
      taxis_t *taxisptr2 = &streamptr->tsteps[tsID].taxis;
      ptaxisCopy(taxisptr2, taxisptr1);
    }

  streamptr->ntsteps = tsID + 1;

#ifdef HAVE_LIBNETCDF
  if ((streamptr->filetype == FILETYPE_NC  ||
       streamptr->filetype == FILETYPE_NC2 ||
       streamptr->filetype == FILETYPE_NC4 ||
       streamptr->filetype == FILETYPE_NC4C)
      && vlistHasTime(vlistID))
    {
      void (*myCdfDefTimestep)(stream_t *streamptr, int tsID)
        = (void (*)(stream_t *, int))
        namespaceSwitchGet(NSSWITCH_CDF_DEF_TIMESTEP).func;
      myCdfDefTimestep(streamptr, tsID);
    }
#endif

  cdi_create_records(streamptr, tsID);

  return (int)streamptr->ntsteps;
}

/*
@Function  streamDefTimestep
@Title     Define time step

@Prototype int streamDefTimestep(int streamID, int tsID)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenWrite}.
    @Item  tsID      Timestep identifier.

@Description
The function @func{streamDefTimestep} defines the time step of a stream.

@Result
@func{streamDefTimestep} returns the number of records of the time step.

@EndFunction
*/
int streamDefTimestep(int streamID, int tsID)
{
  stream_t *streamptr = stream_to_pointer(streamID);
  int (*myStreamDefTimestep_)(stream_t *streamptr, int tsID)
    = (int (*)(stream_t *, int))
    namespaceSwitchGet(NSSWITCH_STREAM_DEF_TIMESTEP_).func;
  return myStreamDefTimestep_(streamptr, tsID);
}

int streamInqCurTimestepID(int streamID)
{
  stream_t *streamptr = stream_to_pointer(streamID);
  return streamptr->curTsID;
}


/*
@Function  streamInqTimestep
@Title     Get time step

@Prototype int streamInqTimestep(int streamID, int tsID)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenRead} or @fref{streamOpenWrite}.
    @Item  tsID      Timestep identifier.

@Description
The function @func{streamInqTimestep} returns the time step of a stream.

@Result
@func{streamInqTimestep} returns the number of records of the time step.

@EndFunction
*/
int streamInqTimestep(int streamID, int tsID)
{
  int nrecs = 0;
  int taxisID;
  stream_t *streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  int vlistID = streamptr->vlistID;

  if ( tsID < streamptr->rtsteps )
    {
      streamptr->curTsID = tsID;
      nrecs = streamptr->tsteps[tsID].nrecs;
      streamptr->tsteps[tsID].curRecID = CDI_UNDEFID;
      taxisID = vlistInqTaxis(vlistID);
      if ( taxisID == -1 )
	Error("Timestep undefined for fileID = %d", streamID);
      ptaxisCopy(taxisPtr(taxisID), &streamptr->tsteps[tsID].taxis);

      return (nrecs);
    }

  if ( tsID >= streamptr->ntsteps && streamptr->ntsteps != CDI_UNDEFID )
    {
      return (0);
    }

  int filetype = streamptr->filetype;

  if ( CDI_Debug )
    Message("streamID = %d  tsID = %d  filetype = %d", streamID, tsID, filetype);

  switch (filetype)
    {
#if  defined  (HAVE_LIBGRIB)
    case FILETYPE_GRB:
    case FILETYPE_GRB2:
      {
        nrecs = grbInqTimestep(streamptr, tsID);
	break;
      }
#endif
#if  defined  (HAVE_LIBSERVICE)
    case FILETYPE_SRV:
      {
        nrecs = srvInqTimestep(streamptr, tsID);
	break;
      }
#endif
#if  defined  (HAVE_LIBEXTRA)
    case FILETYPE_EXT:
      {
        nrecs = extInqTimestep(streamptr, tsID);
	break;
      }
#endif
#if  defined  (HAVE_LIBIEG)
    case FILETYPE_IEG:
      {
        nrecs = iegInqTimestep(streamptr, tsID);
	break;
      }
#endif
#if  defined  (HAVE_LIBNETCDF)
    case FILETYPE_NC:
    case FILETYPE_NC2:
    case FILETYPE_NC4:
    case FILETYPE_NC4C:
      {
        nrecs = cdfInqTimestep(streamptr, tsID);
	break;
      }
#endif
    default:
      {
	Error("%s support not compiled in!", strfiletype(filetype));
	break;
      }
    }

  taxisID = vlistInqTaxis(vlistID);
  if ( taxisID == -1 )
    Error("Timestep undefined for fileID = %d", streamID);

  ptaxisCopy(taxisPtr(taxisID), &streamptr->tsteps[tsID].taxis);

  return (nrecs);
}

/* the single image implementation */
static
void cdiStreamReadVar(int streamID, int varID, int memtype, void *data, int *nmiss)
{
  if ( CDI_Debug ) Message("streamID = %d  varID = %d", streamID, varID);

  check_parg(data);
  check_parg(nmiss);

  stream_t *streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  int filetype = streamptr->filetype;

  *nmiss = 0;

  switch (filetype)
    {
#if  defined  (HAVE_LIBGRIB)
    case FILETYPE_GRB:
    case FILETYPE_GRB2:
      {
        if ( memtype == MEMTYPE_FLOAT ) Error("grbReadVar not implemented for memtype float!");
        grbReadVarDP(streamptr, varID, data, nmiss);
	break;
      }
#endif
#if  defined  (HAVE_LIBSERVICE)
    case FILETYPE_SRV:
      {
        if ( memtype == MEMTYPE_FLOAT ) Error("srvReadVar not implemented for memtype float!");
        srvReadVarDP(streamptr, varID, data, nmiss);
	break;
      }
#endif
#if  defined  (HAVE_LIBEXTRA)
    case FILETYPE_EXT:
      {
        if ( memtype == MEMTYPE_FLOAT ) Error("extReadVar not implemented for memtype float!");
        extReadVarDP(streamptr, varID, data, nmiss);
	break;
      }
#endif
#if  defined  (HAVE_LIBIEG)
    case FILETYPE_IEG:
      {
        if ( memtype == MEMTYPE_FLOAT ) Error("iegReadVar not implemented for memtype float!");
        iegReadVarDP(streamptr, varID, data, nmiss);
	break;
      }
#endif
#if  defined  (HAVE_LIBNETCDF)
    case FILETYPE_NC:
    case FILETYPE_NC2:
    case FILETYPE_NC4:
    case FILETYPE_NC4C:
      {
        if ( memtype == MEMTYPE_FLOAT )
          cdfReadVarSP(streamptr, varID, data, nmiss);
        else
          cdfReadVarDP(streamptr, varID, data, nmiss);

	break;
      }
#endif
    default:
      {
	Error("%s support not compiled in!", strfiletype(filetype));
	break;
      }
    }
}

/*
@Function  streamReadVar
@Title     Read a variable

@Prototype void streamReadVar(int streamID, int varID, double *data, int *nmiss)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenRead}.
    @Item  varID     Variable identifier.
    @Item  data      Pointer to the location into which the data values are read.
                     The caller must allocate space for the returned values.
    @Item  nmiss     Number of missing values.

@Description
The function streamReadVar reads all the values of one time step of a variable
from an open dataset.
@EndFunction
*/
void streamReadVar(int streamID, int varID, double *data, int *nmiss)
{
  cdiStreamReadVar(streamID, varID, MEMTYPE_DOUBLE, data, nmiss);
}

/*
@Function  streamReadVarF
@Title     Read a variable

@Prototype void streamReadVar(int streamID, int varID, float *data, int *nmiss)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenRead}.
    @Item  varID     Variable identifier.
    @Item  data      Pointer to the location into which the data values are read.
                     The caller must allocate space for the returned values.
    @Item  nmiss     Number of missing values.

@Description
The function streamReadVar reads all the values of one time step of a variable
from an open dataset.
@EndFunction
*/
void streamReadVarF(int streamID, int varID, float *data, int *nmiss)
{
  cdiStreamReadVar(streamID, varID, MEMTYPE_FLOAT, data, nmiss);
}

/* the single image implementation */
void cdiStreamWriteVar_(int streamID, int varID, int memtype, const void *data, int nmiss)
{
  if ( CDI_Debug ) Message("streamID = %d varID = %d", streamID, varID);

  check_parg(data);

  stream_t *streamptr = stream_to_pointer(streamID);
  if (subtypeInqActiveIndex(streamptr->vars[varID].subtypeID) != 0)
    Error("Writing of non-trivial subtypes not yet implemented!");

  stream_check_ptr(__func__, streamptr);

  // check taxis
  if ( streamptr->curTsID == CDI_UNDEFID ) streamDefTimestep(streamID, 0);

  int filetype = streamptr->filetype;

  switch (filetype)
    {
#if  defined  (HAVE_LIBGRIB)
    case FILETYPE_GRB:
    case FILETYPE_GRB2:
      {
        grb_write_var(streamptr, varID, memtype, data, nmiss);
	break;
      }
#endif
#if  defined  (HAVE_LIBSERVICE)
    case FILETYPE_SRV:
      {
        if ( memtype == MEMTYPE_FLOAT ) Error("srvWriteVar not implemented for memtype float!");
        srvWriteVarDP(streamptr, varID, data);
	break;
      }
#endif
#if  defined  (HAVE_LIBEXTRA)
    case FILETYPE_EXT:
      {
        if ( memtype == MEMTYPE_FLOAT ) Error("extWriteVar not implemented for memtype float!");
        extWriteVarDP(streamptr, varID, data);
	break;
      }
#endif
#if  defined  (HAVE_LIBIEG)
    case FILETYPE_IEG:
      {
        if ( memtype == MEMTYPE_FLOAT ) Error("iegWriteVar not implemented for memtype float!");
        iegWriteVarDP(streamptr, varID, data);
	break;
      }
#endif
#if  defined  (HAVE_LIBNETCDF)
    case FILETYPE_NC:
    case FILETYPE_NC2:
    case FILETYPE_NC4:
    case FILETYPE_NC4C:
      {
	if ( streamptr->accessmode == 0 ) cdfEndDef(streamptr);
        cdf_write_var(streamptr, varID, memtype, data, nmiss);
	break;
      }
#endif
    default:
      {
	Error("%s support not compiled in!", strfiletype(filetype));
	break;
      }
    }
}

/*
@Function  streamWriteVar
@Title     Write a variable

@Prototype void streamWriteVar(int streamID, int varID, const double *data, int nmiss)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenWrite}.
    @Item  varID     Variable identifier.
    @Item  data      Pointer to a block of double precision floating point data values to be written.
    @Item  nmiss     Number of missing values.

@Description
The function streamWriteVar writes the values of one time step of a variable to an open dataset.
The values are converted to the external data type of the variable, if necessary.
@EndFunction
*/
void streamWriteVar(int streamID, int varID, const double *data, int nmiss)
{
  void (*myCdiStreamWriteVar_)(int streamID, int varID, int memtype,
                               const void *data, int nmiss)
    = (void (*)(int, int, int, const void *, int))
    namespaceSwitchGet(NSSWITCH_STREAM_WRITE_VAR_).func;
  myCdiStreamWriteVar_(streamID, varID, MEMTYPE_DOUBLE, data, nmiss);
}

/*
@Function  streamWriteVarF
@Title     Write a variable

@Prototype void streamWriteVarF(int streamID, int varID, const float *data, int nmiss)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenWrite}.
    @Item  varID     Variable identifier.
    @Item  data      Pointer to a block of single precision floating point data values to be written.
    @Item  nmiss     Number of missing values.

@Description
The function streamWriteVarF writes the values of one time step of a variable to an open dataset.
The values are converted to the external data type of the variable, if necessary.
Only support for netCDF was implemented in this function.
@EndFunction
*/
void streamWriteVarF(int streamID, int varID, const float *data, int nmiss)
{
  void (*myCdiStreamWriteVar_)(int streamID, int varID, int memtype,
                               const void *data, int nmiss)
    = (void (*)(int, int, int, const void *, int))
    namespaceSwitchGet(NSSWITCH_STREAM_WRITE_VAR_).func;
  myCdiStreamWriteVar_(streamID, varID, MEMTYPE_FLOAT, data, nmiss);
}

static
int cdiStreamReadVarSlice(int streamID, int varID, int levelID, int memtype, void *data, int *nmiss)
{
  // May fail if memtype == MEMTYPE_FLOAT and the file format does not support single precision reading.
  // A value > 0 is returned in this case, otherwise it returns zero.
  int status = 0;

  if ( CDI_Debug ) Message("streamID = %d  varID = %d", streamID, varID);

  check_parg(data);
  check_parg(nmiss);

  stream_t *streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  int filetype = streamptr->filetype;

  *nmiss = 0;

  switch (filetype)
    {
#if  defined  (HAVE_LIBGRIB)
    case FILETYPE_GRB:
    case FILETYPE_GRB2:
      {
        if ( memtype == MEMTYPE_FLOAT ) return 1;
        grbReadVarSliceDP(streamptr, varID, levelID, data, nmiss);
	break;
      }
#endif
#if  defined  (HAVE_LIBSERVICE)
    case FILETYPE_SRV:
      {
        if ( memtype == MEMTYPE_FLOAT ) return 1;
        srvReadVarSliceDP(streamptr, varID, levelID, data, nmiss);
	break;
      }
#endif
#if  defined  (HAVE_LIBEXTRA)
    case FILETYPE_EXT:
      {
        if ( memtype == MEMTYPE_FLOAT ) return 1;
        extReadVarSliceDP(streamptr, varID, levelID, data, nmiss);
	break;
      }
#endif
#if  defined  (HAVE_LIBIEG)
    case FILETYPE_IEG:
      {
        if ( memtype == MEMTYPE_FLOAT ) return 1;
        iegReadVarSliceDP(streamptr, varID, levelID, data, nmiss);
	break;
      }
#endif
#if  defined  (HAVE_LIBNETCDF)
    case FILETYPE_NC:
    case FILETYPE_NC2:
    case FILETYPE_NC4:
    case FILETYPE_NC4C:
      {
        if ( memtype == MEMTYPE_FLOAT )
          cdfReadVarSliceSP(streamptr, varID, levelID, data, nmiss);
        else
          cdfReadVarSliceDP(streamptr, varID, levelID, data, nmiss);
        break;
      }
#endif
    default:
      {
	Error("%s support not compiled in!", strfiletype(filetype));
        status = 2;
	break;
      }
    }

  return status;
}

/*
@Function  streamReadVarSlice
@Title     Read a horizontal slice of a variable

@Prototype void streamReadVarSlice(int streamID, int varID, int levelID, double *data, int *nmiss)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenRead}.
    @Item  varID     Variable identifier.
    @Item  levelID   Level identifier.
    @Item  data      Pointer to the location into which the data values are read.
                     The caller must allocate space for the returned values.
    @Item  nmiss     Number of missing values.

@Description
The function streamReadVarSlice reads all the values of a horizontal slice of a variable
from an open dataset.
@EndFunction
*/
void streamReadVarSlice(int streamID, int varID, int levelID, double *data, int *nmiss)
{
  if ( cdiStreamReadVarSlice(streamID, varID, levelID, MEMTYPE_DOUBLE, data, nmiss) )
    {
      Warning("Unexpected error returned from cdiStreamReadVarSlice()!");
      size_t elementCount = (size_t)gridInqSize(vlistInqVarGrid(streamInqVlist(streamID), varID));
      memset(data, 0, elementCount * sizeof(*data));
    }
}

/*
@Function  streamReadVarSliceF
@Title     Read a horizontal slice of a variable

@Prototype void streamReadVarSliceF(int streamID, int varID, int levelID, float *data, int *nmiss)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenRead}.
    @Item  varID     Variable identifier.
    @Item  levelID   Level identifier.
    @Item  data      Pointer to the location into which the data values are read.
                     The caller must allocate space for the returned values.
    @Item  nmiss     Number of missing values.

@Description
The function streamReadVarSliceF reads all the values of a horizontal slice of a variable
from an open dataset.
@EndFunction
*/
void streamReadVarSliceF(int streamID, int varID, int levelID, float *data, int *nmiss)
{
  if ( cdiStreamReadVarSlice(streamID, varID,levelID, MEMTYPE_FLOAT, data, nmiss) )
    {
      // In case the file format does not support single precision reading,
      // we fall back to double precision reading, converting the data on the fly.
      size_t elementCount = (size_t)gridInqSize(vlistInqVarGrid(streamInqVlist(streamID), varID));
      double* conversionBuffer = malloc(elementCount * sizeof(*conversionBuffer));
      streamReadVarSlice(streamID, varID, levelID, conversionBuffer, nmiss);
      for (size_t i = elementCount; i--; ) data[i] = (float)conversionBuffer[i];
      free(conversionBuffer);
    }
}

static
void cdiStreamWriteVarSlice(int streamID, int varID, int levelID, int memtype, const void *data, int nmiss)
{
  if ( CDI_Debug ) Message("streamID = %d varID = %d", streamID, varID);

  check_parg(data);

  stream_t *streamptr = stream_to_pointer(streamID);
  if (subtypeInqActiveIndex(streamptr->vars[varID].subtypeID) != 0)
    Error("Writing of non-trivial subtypes not yet implemented!");

  stream_check_ptr(__func__, streamptr);

  // check taxis
  if ( streamptr->curTsID == CDI_UNDEFID ) streamDefTimestep(streamID, 0);

  int filetype = streamptr->filetype;

  switch (filetype)
    {
#if  defined  (HAVE_LIBGRIB)
    case FILETYPE_GRB:
    case FILETYPE_GRB2:
      {
        grb_write_var_slice(streamptr, varID, levelID, memtype, data, nmiss);
	break;
      }
#endif
#if  defined  (HAVE_LIBSERVICE)
    case FILETYPE_SRV:
      {
        if ( memtype == MEMTYPE_FLOAT ) Error("srvWriteVarSlice not implemented for memtype float!");
        srvWriteVarSliceDP(streamptr, varID, levelID, data);
	break;
      }
#endif
#if  defined  (HAVE_LIBEXTRA)
    case FILETYPE_EXT:
      {
        if ( memtype == MEMTYPE_FLOAT ) Error("extWriteVarSlice not implemented for memtype float!");
        extWriteVarSliceDP(streamptr, varID, levelID, data);
	break;
      }
#endif
#if  defined  (HAVE_LIBIEG)
    case FILETYPE_IEG:
      {
        if ( memtype == MEMTYPE_FLOAT ) Error("iegWriteVarSlice not implemented for memtype float!");
        iegWriteVarSliceDP(streamptr, varID, levelID, data);
	break;
      }
#endif
#if  defined  (HAVE_LIBNETCDF)
    case FILETYPE_NC:
    case FILETYPE_NC2:
    case FILETYPE_NC4:
    case FILETYPE_NC4C:
      if ( streamptr->accessmode == 0 ) cdfEndDef(streamptr);
      cdf_write_var_slice(streamptr, varID, levelID, memtype, data, nmiss);
      break;
#endif
    default:
      {
	Error("%s support not compiled in!", strfiletype(filetype));
	break;
      }
    }
}

/*
@Function  streamWriteVarSlice
@Title     Write a horizontal slice of a variable

@Prototype void streamWriteVarSlice(int streamID, int varID, int levelID, const double *data, int nmiss)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenWrite}.
    @Item  varID     Variable identifier.
    @Item  levelID   Level identifier.
    @Item  data      Pointer to a block of double precision floating point data values to be written.
    @Item  nmiss     Number of missing values.

@Description
The function streamWriteVarSlice writes the values of a horizontal slice of a variable to an open dataset.
The values are converted to the external data type of the variable, if necessary.
@EndFunction
*/
void streamWriteVarSlice(int streamID, int varID, int levelID, const double *data, int nmiss)
{
  cdiStreamWriteVarSlice(streamID, varID, levelID, MEMTYPE_DOUBLE, data, nmiss);
}

/*
@Function  streamWriteVarSliceF
@Title     Write a horizontal slice of a variable

@Prototype void streamWriteVarSliceF(int streamID, int varID, int levelID, const float *data, int nmiss)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenWrite}.
    @Item  varID     Variable identifier.
    @Item  levelID   Level identifier.
    @Item  data      Pointer to a block of single precision floating point data values to be written.
    @Item  nmiss     Number of missing values.

@Description
The function streamWriteVarSliceF writes the values of a horizontal slice of a variable to an open dataset.
The values are converted to the external data type of the variable, if necessary.
Only support for netCDF was implemented in this function.
@EndFunction
*/
void streamWriteVarSliceF(int streamID, int varID, int levelID, const float *data, int nmiss)
{
  cdiStreamWriteVarSlice(streamID, varID, levelID, MEMTYPE_FLOAT, data, nmiss);
}


void
streamWriteVarChunk(int streamID, int varID,
                    const int rect[][2], const double *data, int nmiss)
{
  void (*myCdiStreamWriteVarChunk_)(int streamID, int varID, int memtype,
                                    const int rect[][2], const void *data,
                                    int nmiss)
    = (void (*)(int, int, int, const int [][2], const void *, int))
    namespaceSwitchGet(NSSWITCH_STREAM_WRITE_VAR_CHUNK_).func;
  myCdiStreamWriteVarChunk_(streamID, varID, MEMTYPE_DOUBLE, rect, data, nmiss);
}

/* single image implementation */
void
cdiStreamwriteVarChunk_(int streamID, int varID, int memtype,
                        const int rect[][2], const void *data, int nmiss)
{
  if ( CDI_Debug ) Message("streamID = %d varID = %d", streamID, varID);

  stream_t *streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  // streamDefineTaxis(streamID);

  int filetype = streamptr->filetype;

  switch (filetype)
    {
#if defined (HAVE_LIBGRIB)
    case FILETYPE_GRB:
    case FILETYPE_GRB2:
#endif
#if defined (HAVE_LIBSERVICE)
    case FILETYPE_SRV:
#endif
#if defined (HAVE_LIBEXTRA)
    case FILETYPE_EXT:
#endif
#if defined (HAVE_LIBIEG)
    case FILETYPE_IEG:
#endif
#if  defined (HAVE_LIBGRIB) || defined (HAVE_LIBSERVICE)      \
  || defined (HAVE_LIBEXTRA) || defined (HAVE_LIBIEG)
      xabort("streamWriteVarChunk not implemented for filetype %s!",
             strfiletype(filetype));
      break;
#endif
#if  defined  (HAVE_LIBNETCDF)
    case FILETYPE_NC:
    case FILETYPE_NC2:
    case FILETYPE_NC4:
    case FILETYPE_NC4C:
      if ( streamptr->accessmode == 0 ) cdfEndDef(streamptr);
      cdf_write_var_chunk(streamptr, varID, memtype, rect, data, nmiss);
      break;
#endif
    default:
      Error("%s support not compiled in!", strfiletype(filetype));
      break;
    }
}

void streamWriteContents(int streamID, char *cname)
{
  stream_t *streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  int vlistID = streamptr->vlistID;

  FILE *cnp = fopen(cname, "w");

  if ( cnp == NULL ) SysError(cname);

  fprintf(cnp, "#CDI library version %s\n", cdiLibraryVersion());
  fprintf(cnp, "#\n");

  fprintf(cnp, "filename: %s\n", streamptr->filename);
  int filetype = streamptr->filetype;
  fprintf(cnp, "filetype: %s\n", strfiletype(filetype));

  fprintf(cnp, "#\n");
  fprintf(cnp, "#grids:\n");

  int ngrids = vlistNgrids(vlistID);
  for ( int i = 0; i < ngrids; i++ )
    {
      int gridID   = vlistGrid(vlistID, i);
      int gridtype = gridInqType(gridID);
      int xsize    = gridInqXsize(gridID);
      int ysize    = gridInqYsize(gridID);
      fprintf(cnp, "%4d:%4d:%4d:%4d\n", i+1, gridtype, xsize, ysize);
    }

  fprintf(cnp, "#\n");

  fprintf(cnp, "varID:code:gridID:zaxisID:tsteptype:datatype\n");

  int nvars = vlistNvars(vlistID);
  for ( int varID = 0; varID < nvars; varID++ )
    {
      int code      = vlistInqVarCode(vlistID, varID);
      int gridID    = vlistInqVarGrid(vlistID, varID);
      int zaxisID   = vlistInqVarZaxis(vlistID, varID);
      int tsteptype = vlistInqVarTsteptype(vlistID, varID);
      int datatype  = vlistInqVarDatatype(vlistID, varID);
      fprintf(cnp, "%4d:%4d:%4d:%4d:%4d:%4d:\n",
	      varID+1, code, gridID, zaxisID, tsteptype, datatype);
    }

  fprintf(cnp, "#\n");

  fprintf(cnp, "tsID:nrecs:date:time\n");

  int tsID = 0;
  while (1)
    {
      int nrecs      = streamptr->tsteps[tsID].nallrecs;
      int date       = streamptr->tsteps[tsID].taxis.vdate;
      int time       = streamptr->tsteps[tsID].taxis.vtime;
      off_t position = streamptr->tsteps[tsID].position;

      fprintf(cnp, "%4d:%4d:%4d:%4d:%ld\n",
	      tsID, nrecs, date, time, (long) position);

      if ( streamptr->tsteps[tsID].next )
	tsID++;
      else
	break;
    }

  fprintf(cnp, "#\n");

  fprintf(cnp, "tsID:recID:varID:levID:size:pos\n");

  tsID = 0;
  while (1)
    {
      int nrecs = streamptr->tsteps[tsID].nallrecs;
      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  int varID   = streamptr->tsteps[tsID].records[recID].varID;
	  int levelID = streamptr->tsteps[tsID].records[recID].levelID;
	  off_t recpos = streamptr->tsteps[tsID].records[recID].position;
	  long recsize = (long)streamptr->tsteps[tsID].records[recID].size;
	  fprintf(cnp, "%4d:%4d:%4d:%4d:%4ld:%ld\n",
		  tsID, recID, varID, levelID, recsize, (long) recpos);
	}

      if ( streamptr->tsteps[tsID].next )
	tsID++;
      else
	break;
    }

  fclose(cnp);
}

// This function is used in CDO!
off_t streamNvals(int streamID)
{
  stream_t *streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  return (streamptr->numvals);
}

/*
@Function  streamDefVlist
@Title     Define the variable list

@Prototype void streamDefVlist(int streamID, int vlistID)
@Parameter
    @Item  streamID Stream ID, from a previous call to @fref{streamOpenWrite}.
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.

@Description
The function @func{streamDefVlist} defines the variable list of a stream.

@EndFunction
*/
void streamDefVlist(int streamID, int vlistID)
{
  void (*myStreamDefVlist)(int streamID, int vlistID)
    = (void (*)(int, int))namespaceSwitchGet(NSSWITCH_STREAM_DEF_VLIST_).func;
  myStreamDefVlist(streamID, vlistID);
}

/* the single image implementation of streamDefVlist */
void
cdiStreamDefVlist_(int streamID, int vlistID)
{
  stream_t *streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  if ( streamptr->vlistID == CDI_UNDEFID )
    cdiStreamSetupVlist(streamptr, vlistDuplicate(vlistID), vlistID);
  else
    Warning("vlist already defined for %s!", streamptr->filename);
}

/*
@Function  streamInqVlist
@Title     Get the variable list

@Prototype int streamInqVlist(int streamID)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenRead} or @fref{streamOpenWrite}.

@Description
The function @func{streamInqVlist} returns the variable list of a stream.

@Result
@func{streamInqVlist} returns an identifier to the variable list.

@EndFunction
*/
int streamInqVlist(int streamID)
{
  stream_t *streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  return (streamptr->vlistID);
}


int streamInqVlistIDorig(int streamID)
{
  stream_t *streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  return (streamptr->vlistIDorig);
}


void streamDefCompType(int streamID, int comptype)
{
  stream_t *streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  if (streamptr->comptype != comptype)
    {
      streamptr->comptype = comptype;
      reshSetStatus(streamID, &streamOps, RESH_DESYNC_IN_USE);
    }
}


void streamDefCompLevel(int streamID, int complevel)
{
  stream_t *streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  if (streamptr->complevel != complevel)
    {
      streamptr->complevel = complevel;
      reshSetStatus(streamID, &streamOps, RESH_DESYNC_IN_USE);
    }
}


int streamInqCompType(int streamID)
{
  stream_t *streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  return (streamptr->comptype);
}


int streamInqCompLevel(int streamID)
{
  stream_t *streamptr = stream_to_pointer(streamID);

  stream_check_ptr(__func__, streamptr);

  return (streamptr->complevel);
}

int streamInqFileID(int streamID)
{
  stream_t *streamptr;

  streamptr = ( stream_t *) reshGetVal ( streamID, &streamOps );

  return (streamptr->fileID);
}


void cdiDefAccesstype(int streamID, int type)
{
  stream_t *streamptr = reshGetVal(streamID, &streamOps);

  if ( streamptr->accesstype == CDI_UNDEFID )
    {
      streamptr->accesstype = type;
    }
  else if ( streamptr->accesstype != type )
    Error("Changing access type from %s not allowed!",
          streamptr->accesstype == TYPE_REC ? "REC to VAR" : "VAR to REC");
}


int cdiInqAccesstype(int streamID)
{
  stream_t *streamptr = (stream_t *) reshGetVal ( streamID, &streamOps );

  return (streamptr->accesstype);
}

static
int streamTxCode(void)
{
  return STREAM;
}


void cdiStreamSetupVlist(stream_t *streamptr, int vlistID, int vlistIDorig)
{
  int nvars = vlistNvars(vlistID);
  streamptr->vlistID = vlistID;
  streamptr->vlistIDorig = vlistIDorig;

  /* increase lock counter of the original vlist by 1 */
  if ( CDI_Debug ) Message("attach vlist=%d to stream, stream vlist=%d", vlistIDorig, vlistID);
  vlist_lock(vlistIDorig);

  for (int varID = 0; varID < nvars; varID++ )
    {
      int gridID    = vlistInqVarGrid(vlistID, varID);
      int zaxisID   = vlistInqVarZaxis(vlistID, varID);
      int tilesetID = vlistInqVarSubtype(vlistID, varID);
      stream_new_var(streamptr, gridID, zaxisID, tilesetID);
      if ( streamptr->have_missval )
        vlistDefVarMissval(vlistID, varID,
                           vlistInqVarMissval(vlistID, varID));
    }

  if (streamptr->filemode == 'w' )
    switch (streamptr->filetype)
      {
#ifdef HAVE_LIBNETCDF
      case FILETYPE_NC:
      case FILETYPE_NC2:
      case FILETYPE_NC4:
      case FILETYPE_NC4C:
        {
          void (*myCdfDefVars)(stream_t *streamptr)
            = (void (*)(stream_t *))
            namespaceSwitchGet(NSSWITCH_CDF_STREAM_SETUP).func;
          myCdfDefVars(streamptr);
        }
        break;
#endif
      case FILETYPE_GRB:
      case FILETYPE_GRB2:
        gribContainersNew(streamptr);
        break;
      }
}


void cdiStreamGetIndexList(unsigned numIDs, int IDs[numIDs])
{
  reshGetResHListOfType(numIDs, IDs, &streamOps);
}

int streamInqNvars ( int streamID )
{
  stream_t *streamptr = reshGetVal(streamID, &streamOps);
  return streamptr->nvars;
}


static int streamCompareP(void * streamptr1, void * streamptr2)
{
  stream_t * s1 = ( stream_t * ) streamptr1;
  stream_t * s2 = ( stream_t * ) streamptr2;
  enum {
    differ = -1,
    equal  = 0,
  };

  xassert ( s1 );
  xassert ( s2 );

  if ( s1->filetype  != s2->filetype  ) return differ;
  if (  namespaceAdaptKey2 ( s1->vlistIDorig ) !=
	namespaceAdaptKey2 ( s2->vlistIDorig )) return differ;
  if ( s1->byteorder != s2->byteorder ) return differ;
  if ( s1->comptype  != s2->comptype  ) return differ;
  if ( s1->complevel != s2->complevel ) return differ;

  if ( s1->filename )
    {
      if (strcmp(s1->filename, s2->filename))
	return differ;
    }
  else if ( s2->filename )
    return differ;

  return equal;
}


void streamDestroyP ( void * streamptr )
{
  stream_t * sp = ( stream_t * ) streamptr;

  xassert ( sp );

  int id = sp->self;
  streamClose ( id );
}


void streamPrintP   ( void * streamptr, FILE * fp )
{
  stream_t * sp = ( stream_t * ) streamptr;

  if ( !sp ) return;

  fprintf ( fp, "#\n");
  fprintf ( fp, "# streamID %d\n", sp->self);
  fprintf ( fp, "#\n");
  fprintf ( fp, "self          = %d\n", sp->self );
  fprintf ( fp, "accesstype    = %d\n", sp->accesstype );
  fprintf ( fp, "accessmode    = %d\n", sp->accessmode );
  fprintf ( fp, "filetype      = %d\n", sp->filetype );
  fprintf ( fp, "byteorder     = %d\n", sp->byteorder );
  fprintf ( fp, "fileID        = %d\n", sp->fileID );
  fprintf ( fp, "filemode      = %d\n", sp->filemode );
  fprintf ( fp, "//off_t numvals;\n" );
  fprintf ( fp, "filename      = %s\n", sp->filename );
  fprintf ( fp, "//Record   *record;\n" );
  fprintf ( fp, "nrecs         = %d\n", sp->nrecs );
  fprintf ( fp, "nvars         = %d\n", sp->nvars );
  fprintf ( fp, "//svarinfo_t *vars;\n" );
  fprintf ( fp, "varsAllocated = %d\n", sp->varsAllocated );
  fprintf ( fp, "curTsID       = %d\n", sp->curTsID );
  fprintf ( fp, "rtsteps       = %d\n", sp->rtsteps );
  fprintf ( fp, "//long ntsteps;\n" );
  fprintf ( fp, "//  tsteps_t   *tsteps;\n" );
  fprintf ( fp, "tstepsTableSize= %d\n", sp->tstepsTableSize );
  fprintf ( fp, "tstepsNextID  = %d\n", sp->tstepsNextID );
  fprintf ( fp, "//basetime_t  basetime;\n" );
  fprintf ( fp, "ncmode        = %d\n", sp->ncmode );
  fprintf ( fp, "vlistID       = %d\n", sp->vlistID );
  fprintf ( fp, "//  int       xdimID[MAX_GRIDS_PS];\n" );
  fprintf ( fp, "//  int       ydimID[MAX_GRIDS_PS];\n" );
  fprintf ( fp, "//  int       zaxisID[MAX_ZAXES_PS];\n" );
  fprintf ( fp, "//  int       ncxvarID[MAX_GRIDS_PS];\n" );
  fprintf ( fp, "//  int       ncyvarID[MAX_GRIDS_PS];\n" );
  fprintf ( fp, "//  int       ncavarID[MAX_GRIDS_PS];\n" );
  fprintf ( fp, "historyID     = %d\n", sp->historyID );
  fprintf ( fp, "globalatts    = %d\n", sp->globalatts );
  fprintf ( fp, "localatts     = %d\n", sp->localatts );
  fprintf ( fp, "//  VCT       vct;\n" );
  fprintf ( fp, "unreduced     = %d\n", sp->unreduced );
  fprintf ( fp, "sortname      = %d\n", sp->sortname );
  fprintf ( fp, "have_missval  = %d\n", sp->have_missval );
  fprintf ( fp, "ztype         = %d\n", sp->comptype );
  fprintf ( fp, "zlevel        = %d\n", sp->complevel );
  fprintf ( fp, "//  void    **gribContainers;\n" );
  fprintf ( fp, "vlistIDorig   = %d\n", sp->vlistIDorig );
}

enum {
  streamNint = 11,
};

static int
streamGetPackSize(void * voidP, void *context)
{
  stream_t * streamP = ( stream_t * ) voidP;
  int packBufferSize
    = serializeGetSize(streamNint, DATATYPE_INT, context)
    + serializeGetSize(2, DATATYPE_UINT32, context)
    + serializeGetSize((int)strlen(streamP->filename) + 1,
                       DATATYPE_TXT, context)
    + serializeGetSize(1, DATATYPE_FLT64, context);
  return packBufferSize;
}


static void
streamPack(void * streamptr, void * packBuffer, int packBufferSize,
           int * packBufferPos, void *context)
{
  stream_t * streamP = ( stream_t * ) streamptr;
  int intBuffer[streamNint];

  intBuffer[0]  = streamP->self;
  intBuffer[1]  = streamP->filetype;
  intBuffer[2]  = (int)strlen(streamP->filename) + 1;
  intBuffer[3]  = streamP->vlistID;
  intBuffer[4]  = streamP->vlistIDorig;
  intBuffer[5]  = streamP->byteorder;
  intBuffer[6]  = streamP->comptype;
  intBuffer[7]  = streamP->complevel;
  intBuffer[8]  = streamP->unreduced;
  intBuffer[9]  = streamP->sortname;
  intBuffer[10] = streamP->have_missval;

  serializePack(intBuffer, streamNint, DATATYPE_INT, packBuffer, packBufferSize, packBufferPos, context);
  uint32_t d = cdiCheckSum(DATATYPE_INT, streamNint, intBuffer);
  serializePack(&d, 1, DATATYPE_UINT32, packBuffer, packBufferSize, packBufferPos, context);

  serializePack(&cdiDefaultMissval, 1, DATATYPE_FLT64, packBuffer, packBufferSize, packBufferPos, context);
  serializePack(streamP->filename, intBuffer[2], DATATYPE_TXT, packBuffer, packBufferSize, packBufferPos, context);
  d = cdiCheckSum(DATATYPE_TXT, intBuffer[2], streamP->filename);
  serializePack(&d, 1, DATATYPE_UINT32, packBuffer, packBufferSize, packBufferPos, context);
}

struct streamAssoc
streamUnpack(char * unpackBuffer, int unpackBufferSize,
             int * unpackBufferPos, int originNamespace, void *context)
{
  int intBuffer[streamNint], streamID;
  uint32_t d;
  char filename[CDI_MAX_NAME];

  serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos,
                  intBuffer, streamNint, DATATYPE_INT, context);
  serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos,
                  &d, 1, DATATYPE_UINT32, context);
  xassert(cdiCheckSum(DATATYPE_INT, streamNint, intBuffer) == d);

  serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos,
                  &cdiDefaultMissval, 1, DATATYPE_FLT64, context);
  serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos,
                  &filename, intBuffer[2], DATATYPE_TXT, context);
  serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos,
                  &d, 1, DATATYPE_UINT32, context);
  xassert(d == cdiCheckSum(DATATYPE_TXT, intBuffer[2], filename));
  int targetStreamID = namespaceAdaptKey(intBuffer[0], originNamespace);
  streamID = streamOpenID(filename, "w", intBuffer[1], targetStreamID);
  xassert(streamID >= 0 && targetStreamID == streamID);
  streamDefByteorder(streamID, intBuffer[5]);
  streamDefCompType(streamID, intBuffer[6]);
  streamDefCompLevel(streamID, intBuffer[7]);
  stream_t *streamptr = stream_to_pointer(streamID);
  streamptr->unreduced = intBuffer[8];
  streamptr->sortname = intBuffer[9];
  streamptr->have_missval = intBuffer[10];
  struct streamAssoc retval = { streamID, intBuffer[3], intBuffer[4] };
  return retval;
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
