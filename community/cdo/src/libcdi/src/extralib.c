#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>

#include "dmemory.h"

#include "extra.h"
#include "error.h"
#include "file.h"
#include "binary.h"
#include "stream_fcommon.h"
#include "swap.h"


enum {
  EXT_HEADER_LEN = 4,
};


static int initExtLib       = 0;
static int extDefaultPrec   = 0;
static int extDefaultNumber = EXT_REAL;


/*
 * A version string.
 */

#undef  LIBVERSION
#define LIBVERSION      1.3.2
#define XSTRING(x)	#x
#define STRING(x)	XSTRING(x)
static const char ext_libvers[] = STRING(LIBVERSION) " of "__DATE__" "__TIME__;

const char *extLibraryVersion(void)
{
  return (ext_libvers);
}


static int EXT_Debug = 0;    /* If set to 1, debugging */


void extDebug(int debug)
{
  EXT_Debug = debug;

  if ( EXT_Debug )
    Message("debug level %d", debug);
}


static void extLibInit()
{
  char *envString;
  char *envName = "EXT_PRECISION";


  envString = getenv(envName);
  if ( envString )
    {
      int pos = 0;

      if ( strlen(envString) == 2  )
	{
	  switch ( tolower((int) envString[pos]) )
	    {
	    case 'r':
	      {
		extDefaultNumber = EXT_REAL;
		switch ( (int) envString[pos+1] )
		  {
		  case '4': extDefaultPrec = SINGLE_PRECISION; break;
		  case '8': extDefaultPrec = DOUBLE_PRECISION; break;
		  default:
		    Message("Invalid digit in %s: %s", envName, envString);
		  }
		break;
	      }
	    case 'c':
	      {
		extDefaultNumber = EXT_COMP;
		switch ( (int) envString[pos+1] )
		  {
		  case '4': extDefaultPrec = SINGLE_PRECISION; break;
		  case '8': extDefaultPrec = DOUBLE_PRECISION; break;
		  default:
		    Message("Invalid digit in %s: %s", envName, envString);
		  }
		break;
	      }
	    default:
              {
                Message("Invalid character in %s: %s", envName, envString);
                break;
              }
            }
	}
    }

  initExtLib = 1;
}

static
void extInit(extrec_t *extp)
{
  extp->checked    = 0;
  extp->byteswap   = 0;
  extp->prec       = 0;
  extp->number     = extDefaultNumber;
  extp->datasize   = 0;
  extp->buffersize = 0;
  extp->buffer     = NULL;
}


void *extNew(void)
{
  extrec_t *extp;

  if ( ! initExtLib ) extLibInit();

  extp = (extrec_t *) malloc(sizeof(extrec_t));

  extInit(extp);

  return ((void*)extp);
}


void extDelete(void *ext)
{
  extrec_t *extp = (extrec_t *) ext;

  if ( extp )
    {
      if ( extp->buffer ) free(extp->buffer);
      free(extp);
    }
}


int extCheckFiletype(int fileID, int *swap)
{
  size_t blocklen = 0, fact = 0;
  size_t sblocklen = 0;
  size_t data =  0;
  size_t dimxy = 0;
  int found = 0;
  unsigned char buffer[40], *pbuf;

  if ( fileRead(fileID, buffer, 4) != 4 ) return (found);

  blocklen  = (size_t) get_UINT32(buffer);
  sblocklen = (size_t) get_SUINT32(buffer);

  if ( EXT_Debug )
    Message("blocklen = %d sblocklen = %d", blocklen, sblocklen);

  if ( blocklen == 16 )
    {
     *swap = 0;
      fact = blocklen/4;
      if ( fileRead(fileID, buffer, blocklen+8) != blocklen+8 ) return (found);
      pbuf = buffer+3*fact;      dimxy = (size_t) get_UINT32(pbuf);
      pbuf = buffer+blocklen+4;  data  = (size_t) get_UINT32(pbuf);
    }
  else if ( blocklen == 32 )
    {
     *swap = 0;
      fact = blocklen/4;
      if ( fileRead(fileID, buffer, blocklen+8) != blocklen+8 ) return (found);
      pbuf = buffer+3*fact;      dimxy = (size_t) get_UINT64(pbuf);
      pbuf = buffer+blocklen+4;  data  = (size_t) get_UINT32(pbuf);
    }
  else if ( sblocklen == 16 )
    {
     *swap = 1;
      fact = sblocklen/4;
      if ( fileRead(fileID, buffer, sblocklen+8) != sblocklen+8 ) return (found);
      pbuf = buffer+3*fact;       dimxy = (size_t) get_SUINT32(pbuf);
      pbuf = buffer+sblocklen+4;  data  = (size_t) get_SUINT32(pbuf);
    }
  else if ( sblocklen == 32 )
    {
     *swap = 1;
      fact = sblocklen/4;
      if ( fileRead(fileID, buffer, sblocklen+8) != sblocklen+8 ) return (found);
      pbuf = buffer+3*fact;       dimxy = (size_t) get_SUINT64(pbuf);
      pbuf = buffer+sblocklen+4;  data  = (size_t) get_SUINT32(pbuf);
    }

  fileRewind(fileID);

  if      ( data && dimxy*fact   == data ) found = 1;
  else if ( data && dimxy*fact*2 == data ) found = 1;

  if ( EXT_Debug )
    {
      Message("swap = %d fact = %d", *swap, fact);
      Message("dimxy = %lu data = %lu", dimxy, data);
    }

  return (found);
}


int extInqHeader(void *ext, int *header)
{
  extrec_t *extp = (extrec_t *) ext;
  size_t i;

  for ( i = 0; i < EXT_HEADER_LEN; i++ )
    header[i] = extp->header[i];

  if ( EXT_Debug ) Message("datasize = %lu", extp->datasize);

  return (0);
}


int extDefHeader(void *ext, const int *header)
{
  extrec_t *extp = (extrec_t *) ext;
  size_t i;

  for ( i = 0; i < EXT_HEADER_LEN; i++ )
    extp->header[i] = header[i];

  extp->datasize = (size_t)header[3];
  if ( extp->number == EXT_COMP ) extp->datasize *= 2;

  if ( EXT_Debug ) Message("datasize = %lu", extp->datasize);

  return (0);
}


static int extInqData(extrec_t *extp, int prec, void *data)
{
  size_t datasize;
  size_t i;
  int ierr = 0;
  int rprec;
  void *buffer;
  int byteswap = extp->byteswap;

  datasize = extp->datasize;
  buffer   = extp->buffer;
  rprec    = extp->prec;

  switch ( rprec )
    {
    case SINGLE_PRECISION:
      {
	if ( sizeof(FLT32) == 4 )
	  {
	    if ( byteswap ) swap4byte(buffer, datasize);

	    if ( rprec == prec )
	      memcpy(data, buffer, datasize*sizeof(FLT32));
	    else
	      for ( i = 0; i < datasize; ++i )
		((double *) data)[i] = (double) ((float *) buffer)[i];
	  }
	else
	  {
	    Error("not implemented for %d byte float", sizeof(FLT32));
	  }
	break;
      }
    case DOUBLE_PRECISION:
	if ( sizeof(FLT64) == 8 )
	  {
	    if ( byteswap ) swap8byte(buffer, datasize);

	    if ( rprec == prec )
	      memcpy(data, buffer, datasize*sizeof(FLT64));
	    else
	      for ( i = 0; i < datasize; ++i )
		((float *) data)[i] = (float) ((double *) buffer)[i];
	  }
	else
	  {
	    Error("not implemented for %d byte float", sizeof(FLT64));
	  }
	break;
    default:
      {
	Error("unexpected data precision %d", rprec);
	break;
      }
    }

  return (ierr);
}


int extInqDataSP(void *ext, float *data)
{
  return (extInqData(ext, SINGLE_PRECISION, (void *) data));
}


int extInqDataDP(void *ext, double *data)
{
  return (extInqData(ext, DOUBLE_PRECISION, (void *) data));
}


int extDefData(void *ext, int prec, const void *data)
{
  extrec_t *extp = (extrec_t *) ext;
  size_t datasize;
  size_t blocklen;
  size_t buffersize;
  size_t i;
  int rprec;
  int *header;
  void *buffer;

  if ( extDefaultPrec ) rprec = extDefaultPrec;
  else                  rprec = extp->prec;

  if ( ! rprec ) rprec = prec;

  extp->prec = rprec;

  header = extp->header;

  datasize = (size_t)header[3];
  if ( extp->number == EXT_COMP ) datasize *= 2;
  blocklen = datasize * (size_t)rprec;

  extp->datasize = datasize;

  buffersize = extp->buffersize;

  if ( buffersize != blocklen )
    {
      buffersize = blocklen;
      buffer = extp->buffer;
      buffer = realloc(buffer, buffersize);
      extp->buffer = buffer;
      extp->buffersize = buffersize;
    }
  else
    buffer = extp->buffer;

  switch ( rprec )
    {
    case SINGLE_PRECISION:
      {
	if ( rprec == prec )
	  memcpy(buffer, data, datasize*sizeof(FLT32));
	else
	  for (i = 0; i < datasize; i++)
	    ((float *) buffer)[i] = (float) ((double *) data)[i];

	break;
      }
    case DOUBLE_PRECISION:
      {
	if ( rprec == prec )
	  memcpy(buffer, data, datasize*sizeof(FLT64));
	else
	  for (i = 0; i < datasize; i++)
	    ((double *) buffer)[i] = (double) ((float *) data)[i];

	break;
      }
    default:
      {
	Error("unexpected data precision %d", rprec);
        break;
      }
    }

  return (0);
}


int extDefDataSP(void *ext, const float *data)
{
  return (extDefData(ext, SINGLE_PRECISION, (void *) data));
}


int extDefDataDP(void *ext, const double *data)
{
  return (extDefData(ext, DOUBLE_PRECISION, (void *) data));
}


int extRead(int fileID, void *ext)
{
  extrec_t *extp = (extrec_t *) ext;
  size_t blocklen, blocklen2;
  size_t i;
  void *buffer;
  int byteswap;
  int status;

  if ( ! extp->checked )
    {
      status = extCheckFiletype(fileID, &extp->byteswap);
      if ( status == 0 ) Error("Not a EXTRA file!");
      extp->checked = 1;
    }

  byteswap = extp->byteswap;

  /* read header record */
  blocklen = binReadF77Block(fileID, byteswap);

  if ( fileEOF(fileID) ) return (-1);

  if ( EXT_Debug )
    Message("blocklen = %lu", blocklen);

  size_t hprec = blocklen / EXT_HEADER_LEN;

  extp->prec = (int)hprec;

  switch ( hprec )
    {
    case SINGLE_PRECISION:
      {
        INT32 tempheader[4];
	binReadInt32(fileID, byteswap, EXT_HEADER_LEN, tempheader);

	for ( i = 0; i < EXT_HEADER_LEN; i++ )
          extp->header[i] = (int)tempheader[i];

	break;
      }
    case DOUBLE_PRECISION:
      {
        INT64 tempheader[4];
	binReadInt64(fileID, byteswap, EXT_HEADER_LEN, tempheader);

	for ( i = 0; i < EXT_HEADER_LEN; i++ )
          extp->header[i] = (int)tempheader[i];

	break;
      }
    default:
      {
	Error("Unexpected header precision %d", hprec);
        break;
      }
    }

  blocklen2 = binReadF77Block(fileID, byteswap);

  if ( blocklen2 != blocklen )
    {
      Warning("Header blocklen differ (blocklen1=%d; blocklen2=%d)!", blocklen, blocklen2);
      if ( blocklen2 != 0 ) return (-1);
    }

  extp->datasize = (size_t)extp->header[3];

  if ( EXT_Debug ) Message("datasize = %lu", extp->datasize);

  blocklen = binReadF77Block(fileID, byteswap);

  size_t buffersize = (size_t)extp->buffersize;

  if ( buffersize < blocklen )
    {
      buffersize = blocklen;
      buffer = extp->buffer;
      buffer = realloc(buffer, buffersize);
      extp->buffer = buffer;
      extp->buffersize = buffersize;
    }
  else
    buffer = extp->buffer;

  size_t dprec = blocklen / extp->datasize;

  if ( dprec == hprec )
    {
      extp->number = EXT_REAL;
    }
  else if ( dprec == 2*hprec )
    {
      dprec /= 2;
      extp->datasize *= 2;
      extp->number = EXT_COMP;
    }

  if ( dprec != SINGLE_PRECISION && dprec != DOUBLE_PRECISION )
    {
      Warning("Unexpected data precision %d", dprec);
      return (-1);
    }

  fileRead(fileID, buffer, blocklen);

  blocklen2 = binReadF77Block(fileID, byteswap);

  if ( blocklen2 != blocklen )
    {
      Warning("Data blocklen differ (blocklen1=%d; blocklen2=%d)!", blocklen, blocklen2);
      if ( blocklen2 != 0 ) return (-1);
    }

  return (0);
}


int extWrite(int fileID, void *ext)
{
  extrec_t *extp = (extrec_t *) ext;
  size_t datasize;
  size_t blocklen;
  size_t i;
  int rprec, number;
  char tempheader[32];
  int *header;
  void *buffer;
  int byteswap = extp->byteswap;


  rprec  = extp->prec;
  number = extp->number;
  header = extp->header;

  /* write header record */
  blocklen = EXT_HEADER_LEN * (size_t)rprec;

  binWriteF77Block(fileID, byteswap, blocklen);

  switch ( rprec )
    {
    case SINGLE_PRECISION:
      {
	for (i = 0; i < EXT_HEADER_LEN; i++)
          ((INT32 *) tempheader)[i] = (INT32) header[i];

	binWriteInt32(fileID, byteswap, EXT_HEADER_LEN, (INT32 *) tempheader);

	break;
      }
    case DOUBLE_PRECISION:
      {
	for (i = 0; i < EXT_HEADER_LEN; i++)
          ((INT64 *) tempheader)[i] = (INT64) header[i];

	binWriteInt64(fileID, byteswap, EXT_HEADER_LEN, (INT64 *) tempheader);

	break;
      }
    default:
      {
	Error("unexpected header precision %d", rprec);
        break;
      }
    }

  binWriteF77Block(fileID, byteswap, blocklen);

  datasize = (size_t)header[3];
  if ( number == EXT_COMP ) datasize *= 2;
  blocklen = datasize * (size_t)rprec;

  binWriteF77Block(fileID, byteswap, blocklen);

  extp->datasize = datasize;

  buffer = extp->buffer;

  switch ( rprec )
    {
    case SINGLE_PRECISION:
      {
	binWriteFlt32(fileID, byteswap, datasize, (FLT32 *) buffer);
	break;
      }
    case DOUBLE_PRECISION:
      {
	binWriteFlt64(fileID, byteswap, datasize, (FLT64 *) buffer);
	break;
      }
    default:
      {
	Error("unexpected data precision %d", rprec);
        break;
      }
    }

  binWriteF77Block(fileID, byteswap, blocklen);

  return (0);
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
