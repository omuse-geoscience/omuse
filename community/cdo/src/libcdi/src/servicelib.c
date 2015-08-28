#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#ifdef HAVE_LIBSERVICE

#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>

#include "dmemory.h"

#include "service.h"
#include "error.h"
#include "file.h"
#include "binary.h"
#include "stream_fcommon.h"
#include "swap.h"


enum {
  SRV_HEADER_LEN = 8,
};


static int initSrvLib      = 0;
static int srvDefaultHprec = 0;
static int srvDefaultDprec = 0;


/*
 * A version string.
 */

#undef  LIBVERSION
#define LIBVERSION      1.3.2
#define XSTRING(x)	#x
#define STRING(x)	XSTRING(x)
static const char srv_libvers[] = STRING(LIBVERSION) " of "__DATE__" "__TIME__;

const char *srvLibraryVersion(void)
{
  return (srv_libvers);
}


int SRV_Debug = 0;    /* If set to 1, debugging */


void srvDebug(int debug)
{
  SRV_Debug = debug;

  if ( SRV_Debug )
    Message("debug level %d", debug);
}


void srvLibInit()
{
  char *envString;
  char *envName = "SRV_PRECISION";


  envString = getenv(envName);
  if ( envString )
    {
      int pos;
      int nrun;
      if ( strlen(envString) == 2 ) nrun = 1;
      else                          nrun = 2;

      pos = 0;
      while ( nrun-- )
	{
	  switch ( tolower((int) envString[pos]) )
	    {
	    case 'i':
	      {
		switch ( (int) envString[pos+1] )
		  {
		  case '4': srvDefaultHprec = SINGLE_PRECISION; break;
		  case '8': srvDefaultHprec = DOUBLE_PRECISION; break;
		  default:
		    Message("Invalid digit in %s: %s", envName, envString);
		  }
		break;
	      }
	    case 'r':
	      {
		switch ( (int) envString[pos+1] )
		  {
		  case '4': srvDefaultDprec = SINGLE_PRECISION; break;
		  case '8': srvDefaultDprec = DOUBLE_PRECISION; break;
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
	  pos += 2;
	}
    }

  initSrvLib = 1;
}


void srvInit(srvrec_t *srvp)
{
  srvp->checked    = 0;
  srvp->byteswap   = 0;
  srvp->hprec      = 0;
  srvp->dprec      = 0;
  srvp->datasize   = 0;
  srvp->buffersize = 0;
  srvp->buffer     = NULL;
}


srvrec_t *srvNew(void)
{
  srvrec_t *srvp;

  if ( ! initSrvLib ) srvLibInit();

  srvp = (srvrec_t *) malloc(sizeof(srvrec_t));

  srvInit(srvp);

  return (srvp);
}


void srvDelete(srvrec_t *srvp)
{
  if ( srvp )
    {
      if ( srvp->buffer ) free(srvp->buffer);
      free(srvp);
    }
}


int srvCheckFiletype(int fileID, int *swap)
{
  size_t blocklen = 0;
  size_t sblocklen = 0;
  size_t data = 0;
  size_t dimx = 0, dimy = 0;
  size_t fact = 0;
  int found = 0;
  unsigned char buffer[72], *pbuf;

  if ( fileRead(fileID, buffer, 4) != 4 ) return (found);

  blocklen  = (size_t) get_UINT32(buffer);
  sblocklen = (size_t) get_SUINT32(buffer);

  if ( SRV_Debug )
    Message("blocklen = %d sblocklen = %d", blocklen, sblocklen);

  if ( blocklen == 32 )
    {
     *swap = 0;
      fact = blocklen>>3;
      if ( fileRead(fileID, buffer, blocklen+8) != blocklen+8 ) return (found);
      pbuf = buffer+4*fact;      dimx = (size_t) get_UINT32(pbuf);
      pbuf = buffer+5*fact;      dimy = (size_t) get_UINT32(pbuf);
      pbuf = buffer+blocklen+4;  data = (size_t) get_UINT32(pbuf);
    }
  else if ( blocklen == 64 )
    {
     *swap = 0;
      fact = blocklen>>3;
      if ( fileRead(fileID, buffer, blocklen+8) != blocklen+8 ) return (found);
      pbuf = buffer+4*fact;      dimx = (size_t) get_UINT64(pbuf);
      pbuf = buffer+5*fact;      dimy = (size_t) get_UINT64(pbuf);
      pbuf = buffer+blocklen+4;  data = (size_t) get_UINT32(pbuf);
    }
  else if ( sblocklen == 32 )
    {
     *swap = 1;
      fact = sblocklen>>3;
      if ( fileRead(fileID, buffer, sblocklen+8) != sblocklen+8 ) return (found);
      pbuf = buffer+4*fact;       dimx = (size_t) get_SUINT32(pbuf);
      pbuf = buffer+5*fact;       dimy = (size_t) get_SUINT32(pbuf);
      pbuf = buffer+sblocklen+4;  data = (size_t) get_SUINT32(pbuf);
    }
  else if ( sblocklen == 64 )
    {
     *swap = 1;
      fact = sblocklen>>3;
      if ( fileRead(fileID, buffer, sblocklen+8) != sblocklen+8 ) return (found);
      pbuf = buffer+4*fact;       dimx = (size_t) get_SUINT64(pbuf);
      pbuf = buffer+5*fact;       dimy = (size_t) get_SUINT64(pbuf);
      pbuf = buffer+sblocklen+4;  data = (size_t) get_SUINT32(pbuf);
    }

  fileRewind(fileID);

  if      ( data && dimx*dimy*fact == data ) found = 1;
  else if ( data && dimx*dimy*8    == data ) found = 1;

  if ( SRV_Debug )
    {
      Message("swap = %d fact = %d", *swap, fact);
      Message("dimx = %lu dimy = %lu data = %lu", dimx, dimy, data);
    }

  return (found);
}


int srvInqHeader(srvrec_t *srvp, int *header)
{
  size_t i;

  for ( i = 0; i < SRV_HEADER_LEN; i++ )
    header[i] = srvp->header[i];
  
  if ( SRV_Debug )
    Message("datasize = %lu", srvp->datasize);

  return (0);
}


int srvDefHeader(srvrec_t *srvp, const int *header)
{
  size_t i;

  for ( i = 0; i < SRV_HEADER_LEN; i++ )
    srvp->header[i] = header[i];

  srvp->datasize = (size_t)(header[4] * header[5]);

  if ( SRV_Debug )
    Message("datasize = %lu", srvp->datasize);

  return (0);
}


int srvInqData(srvrec_t *srvp, int prec, void *data)
{
  size_t datasize;
  size_t i;
  int ierr = 0;
  int dprec;
  void *buffer;
  int byteswap = srvp->byteswap;

  datasize = srvp->datasize;

  buffer = srvp->buffer;

  dprec = srvp->dprec;

  switch ( dprec )
    {
    case SINGLE_PRECISION:
      {
	if ( sizeof(FLT32) == 4 )
	  {
	    if ( byteswap ) swap4byte(buffer, datasize);

	    if ( dprec == prec )
	      memcpy(data, buffer, datasize*sizeof(FLT32));
	    else
	      for (i = 0; i < datasize; i++)
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

	    if ( dprec == prec )
	      memcpy(data, buffer, datasize*sizeof(FLT64));
	    else
	      for (i = 0; i < datasize; i++)
		((float *) data)[i] = (float) ((double *) buffer)[i];
	  }
	else
	  {
	    Error("not implemented for %d byte float", sizeof(FLT64));
	  }
	break;
    default:
      {
	Error("unexpected data precision %d", dprec);
        break;
      }
    }

  return (ierr);
}


int srvInqDataSP(srvrec_t *srvp, float *data)
{
  return (srvInqData(srvp, SINGLE_PRECISION, (void *) data));
}


int srvInqDataDP(srvrec_t *srvp, double *data)
{
  return (srvInqData(srvp, DOUBLE_PRECISION, (void *) data));
}


int srvDefData(srvrec_t *srvp, int prec, const void *data)
{
  size_t datasize;
  size_t blocklen;
  size_t buffersize;
  size_t i;
  int dprec, hprec;
  int *header;
  void *buffer;

  if ( srvDefaultDprec ) dprec = srvDefaultDprec;
  else                   dprec = srvp->dprec;

  if ( ! dprec ) dprec = prec;

  srvp->dprec = dprec;

  if ( srvDefaultHprec ) hprec = srvDefaultHprec;
  else                   hprec = srvp->hprec;

  if ( ! hprec ) hprec = dprec;

  srvp->hprec = hprec;

  header = srvp->header;

  datasize = (size_t)(header[4] * header[5]);
  blocklen = datasize * (size_t)dprec;

  srvp->datasize = datasize;

  buffersize = srvp->buffersize;

  if ( buffersize != blocklen )
    {
      buffersize = blocklen;
      buffer = srvp->buffer;
      buffer = realloc(buffer, buffersize);
      srvp->buffer = buffer;
      srvp->buffersize = buffersize;
    }
  else
    buffer = srvp->buffer;

  switch ( dprec )
    {
    case SINGLE_PRECISION:
      {
	if ( dprec == prec )
	  memcpy(buffer, data, datasize*sizeof(FLT32));
	else
	  for (i = 0; i < datasize; i++)
	    ((float *) buffer)[i] = (float) ((double *) data)[i];

	break;
      }
    case DOUBLE_PRECISION:
      {
	if ( dprec == prec )
	  memcpy(buffer, data, datasize*sizeof(FLT64));
	else
	  for (i = 0; i < datasize; i++)
	    ((double *) buffer)[i] = (double) ((float *) data)[i];

	break;
      }
    default:
      {
	Error("unexpected data precision %d", dprec);
        break;
      }
    }

  return (0);
}


int srvDefDataSP(srvrec_t *srvp, const float *data)
{
  return (srvDefData(srvp, SINGLE_PRECISION, (void *) data));
}


int srvDefDataDP(srvrec_t *srvp, const double *data)
{
  return (srvDefData(srvp, DOUBLE_PRECISION, (void *) data));
}


int srvRead(int fileID, srvrec_t *srvp)
{
  size_t datasize;
  size_t blocklen, blocklen2;
  size_t i;
  char tempheader[64];
  void *buffer;
  int byteswap;
  int status;

  if ( ! srvp->checked )
    {
      status = srvCheckFiletype(fileID, &srvp->byteswap);
      if ( status == 0 ) Error("Not a SERVICE file!");
      srvp->checked = 1;
    }

  byteswap = srvp->byteswap;

  /* read header record */
  blocklen = binReadF77Block(fileID, byteswap);

  if ( fileEOF(fileID) ) return (-1);

  if ( SRV_Debug )
    Message("blocklen = %lu", blocklen);

  size_t hprec = blocklen / SRV_HEADER_LEN;

  srvp->hprec = (int)hprec;

  switch ( hprec )
    {
    case SINGLE_PRECISION:
      {
	binReadInt32(fileID, byteswap, SRV_HEADER_LEN, (INT32 *) tempheader);

	for ( i = 0; i < SRV_HEADER_LEN; i++ )
          srvp->header[i] = (int) ((INT32 *) tempheader)[i];

	break;
      }
    case DOUBLE_PRECISION:
      {
	binReadInt64(fileID, byteswap, SRV_HEADER_LEN, (INT64 *) tempheader);

	for ( i = 0; i < SRV_HEADER_LEN; i++ )
          srvp->header[i] = (int) ((INT64 *) tempheader)[i];

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

  srvp->datasize = (size_t)(srvp->header[4] * srvp->header[5]);

  if ( SRV_Debug )
    Message("datasize = %lu", srvp->datasize);

  blocklen = binReadF77Block(fileID, byteswap);

  size_t buffersize = srvp->buffersize;

  if ( buffersize < blocklen )
    {
      buffersize = blocklen;
      buffer = srvp->buffer;
      buffer = realloc(buffer, buffersize);
      srvp->buffer = buffer;
      srvp->buffersize = buffersize;
    }
  else
    buffer = srvp->buffer;

  datasize = srvp->datasize;

  size_t dprec = blocklen / datasize;

  srvp->dprec = (int)dprec;

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


void srvWrite(int fileID, srvrec_t *srvp)
{
  size_t datasize;
  size_t blocklen;
  size_t i;
  int dprec, hprec;
  char tempheader[64];
  int *header;
  void *buffer;
  int byteswap = srvp->byteswap;

  dprec  = srvp->dprec;
  hprec  = srvp->hprec;
  header = srvp->header;

  /* write header record */
  blocklen = SRV_HEADER_LEN * (size_t)hprec;

  binWriteF77Block(fileID, byteswap, blocklen);

  switch ( hprec )
    {
    case SINGLE_PRECISION:
      {
	for (i = 0; i < SRV_HEADER_LEN; i++)
          ((INT32 *) tempheader)[i] = (INT32) header[i];

	binWriteInt32(fileID, byteswap, SRV_HEADER_LEN, (INT32 *) tempheader);

	break;
      }
    case DOUBLE_PRECISION:
      {
	for (i = 0; i < SRV_HEADER_LEN; i++)
          ((INT64 *) tempheader)[i] = (INT64) header[i];

	binWriteInt64(fileID, byteswap, SRV_HEADER_LEN, (INT64 *) tempheader);

	break;
      }
    default:
      {
	Error("unexpected header precision %d", hprec);
        break;
      }
    }

  binWriteF77Block(fileID, byteswap, blocklen);

  datasize = (size_t)(header[4] * header[5]);
  blocklen = datasize * (size_t)dprec;

  binWriteF77Block(fileID, byteswap, blocklen);

  srvp->datasize = datasize;

  buffer = srvp->buffer;

  switch ( dprec )
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
	Error("unexpected data precision %d", dprec);
        break;
      }
    }

  binWriteF77Block(fileID, byteswap, blocklen);
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
#endif  /* HAVE_LIBSERVICE */
