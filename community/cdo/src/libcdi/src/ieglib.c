#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>

#include "dmemory.h"

#include "ieg.h"
#include "error.h"
#include "file.h"
#include "binary.h"
#include "stream_fcommon.h"
#include "swap.h"


static int initIegLib      = 0;
static int iegDefaultDprec = 0;


/*
 * A version string.
 */

#undef  LIBVERSION
#define LIBVERSION      1.3.3
#define XSTRING(x)	#x
#define STRING(x)	XSTRING(x)
static const char ieg_libvers[] = STRING(LIBVERSION) " of "__DATE__" "__TIME__;

const char *iegLibraryVersion(void)
{
  return (ieg_libvers);
}


int IEG_Debug = 0;    /* If set to 1, debugging */



void iegLibInit()
{
  char *envString;
  char *envName = "IEG_PRECISION";


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
	    case 'r':
	      {
		switch ( (int) envString[pos+1] )
		  {
		  case '4': iegDefaultDprec = SINGLE_PRECISION; break;
		  case '8': iegDefaultDprec = DOUBLE_PRECISION; break;
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

  initIegLib = 1;
}


void iegDebug(int debug)
{
  IEG_Debug = debug;

  if ( IEG_Debug )
    Message("debug level %d", debug);
}


void iegInit(iegrec_t *iegp)
{
  iegp->checked    = 0;
  iegp->byteswap   = 0;
  iegp->dprec      = 0;
  iegp->refval     = 0;
  iegp->datasize   = 0;
  iegp->buffersize = 0;
  iegp->buffer     = NULL;
}


void iegInitMem(iegrec_t *iegp)
{
  memset(iegp->ipdb, 0, sizeof(iegp->ipdb));
  memset(iegp->igdb, 0, sizeof(iegp->igdb));
  memset(iegp->vct,  0, sizeof(iegp->vct));
}


iegrec_t *iegNew(void)
{
  iegrec_t *iegp;

  if ( ! initIegLib ) iegLibInit();

  iegp = (iegrec_t *) malloc(sizeof(iegrec_t));

  iegInit(iegp);
  iegInitMem(iegp);

  return (iegp);
}


void iegDelete(iegrec_t *iegp)
{
  if ( iegp )
    {
      if ( iegp->buffer ) free(iegp->buffer);
      free(iegp);
    }
}


int iegCheckFiletype(int fileID, int *swap)
{
  size_t blocklen = 0;
  size_t sblocklen = 0;
  size_t data = 0;
  size_t dimx = 0, dimy = 0;
  size_t fact = 0;
  unsigned char buffer[1048], *pbuf;

  if ( fileRead(fileID, buffer, 4) != 4 ) return (0);

  blocklen  = get_UINT32(buffer);
  sblocklen = get_SUINT32(buffer);

  if ( IEG_Debug )
    Message("blocklen = %d sblocklen = %d", blocklen, sblocklen);

  if ( blocklen == 636 || blocklen == 640 )
    {
     *swap = 0;
      fact = 4;
      if ( fileRead(fileID, buffer, blocklen+8) != blocklen+8 ) return (0);
      pbuf = buffer+(37+4)*4;    dimx = (size_t) get_UINT32(pbuf);
      pbuf = buffer+(37+5)*4;    dimy = (size_t) get_UINT32(pbuf);
      pbuf = buffer+blocklen+4;  data = (size_t) get_UINT32(pbuf);
    }
  else if ( blocklen == 1040 || blocklen == 1036 )
    {
     *swap = 0;
      fact = 8;
      if ( fileRead(fileID, buffer, blocklen+8) != blocklen+8 ) return (0);
      pbuf = buffer+(37+4)*4;    dimx = (size_t) get_UINT32(pbuf);
      pbuf = buffer+(37+5)*4;    dimy = (size_t) get_UINT32(pbuf);
      pbuf = buffer+blocklen+4;  data = (size_t) get_UINT32(pbuf);
    }
  else if ( sblocklen == 636 || sblocklen == 640 )
    {
     *swap = 1;
      fact = 4;
      if ( fileRead(fileID, buffer, sblocklen+8) != sblocklen+8 ) return (0);
      pbuf = buffer+(37+4)*4;     dimx = (size_t) get_SUINT32(pbuf);
      pbuf = buffer+(37+5)*4;     dimy = (size_t) get_SUINT32(pbuf);
      pbuf = buffer+sblocklen+4;  data = (size_t) get_SUINT32(pbuf);
    }
  else if ( sblocklen == 1040 || sblocklen == 1036 )
    {
     *swap = 1;
      fact = 8;
      if ( fileRead(fileID, buffer, sblocklen+8) != sblocklen+8 ) return (0);
      pbuf = buffer+(37+4)*4;     dimx = (size_t) get_SUINT32(pbuf);
      pbuf = buffer+(37+5)*4;     dimy = (size_t) get_SUINT32(pbuf);
      pbuf = buffer+sblocklen+4;  data = (size_t) get_SUINT32(pbuf);
    }

  fileRewind(fileID);

  int found = data && (dimx*dimy*fact == data || dimx*dimy*8 == data);

  if ( IEG_Debug )
    {
      Message("swap = %d fact = %d", *swap, fact);
      Message("dimx = %lu dimy = %lu data = %lu", dimx, dimy, data);
    }

  return (found);
}


void iegCopyMeta(iegrec_t *diegp, iegrec_t *siegp)
{
  /*  diegp->byteswap = siegp->byteswap; */
  diegp->dprec    = siegp->dprec;
  diegp->refval   = siegp->refval;

  memcpy(diegp->ipdb, siegp->ipdb, sizeof(siegp->ipdb));
  memcpy(diegp->igdb, siegp->igdb, sizeof(siegp->igdb));
  memcpy(diegp->vct,  siegp->vct,  sizeof(siegp->vct));
}


int iegInqData(iegrec_t *iegp, int prec, void *data)
{
  size_t datasize;
  size_t i;
  int ierr = 0;
  int dprec;
  void *buffer;
  int byteswap = iegp->byteswap;


  datasize = iegp->datasize;

  buffer = iegp->buffer;

  dprec = iegp->dprec;

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


int iegInqDataSP(iegrec_t *iegp, float *data)
{
  return (iegInqData(iegp, SINGLE_PRECISION, (void *) data));
}


int iegInqDataDP(iegrec_t *iegp, double *data)
{
  return (iegInqData(iegp, DOUBLE_PRECISION, (void *) data));
}


int iegDefData(iegrec_t *iegp, int prec, const void *data)
{
  size_t datasize;
  size_t blocklen;
  size_t buffersize;
  size_t i;
  int dprec;
  void *buffer;


  if ( iegDefaultDprec ) dprec = iegDefaultDprec;
  else                   dprec = iegp->dprec;

  if ( ! dprec ) dprec = prec;

  iegp->dprec = dprec;

  datasize = (size_t)IEG_G_NumLon(iegp->igdb) * (size_t)IEG_G_NumLat(iegp->igdb);
  blocklen = datasize * (size_t)dprec;

  iegp->datasize = datasize;

  buffersize = iegp->buffersize;

  if ( buffersize != blocklen )
    {
      buffersize = blocklen;
      buffer = iegp->buffer;
      buffer = realloc(buffer, buffersize);
      iegp->buffer = buffer;
      iegp->buffersize = buffersize;
    }
  else
    buffer = iegp->buffer;

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


int iegDefDataSP(iegrec_t *iegp, const float *data)
{
  return (iegDefData(iegp, SINGLE_PRECISION, (void *) data));
}


int iegDefDataDP(iegrec_t *iegp, const double *data)
{
  return (iegDefData(iegp, DOUBLE_PRECISION, (void *) data));
}


int iegRead(int fileID, iegrec_t *iegp)
{
  size_t datasize;
  size_t blocklen, blocklen2;
  size_t i;
  char tmpbuffer[800], *tmpbuf = tmpbuffer;
  int dprec = 0;
  void *buffer;
  int byteswap;
  int status;

  if ( ! iegp->checked )
    {
      status = iegCheckFiletype(fileID, &iegp->byteswap);
      if ( status == 0 ) Error("Not a IEG file!");
      iegp->checked = 1;
    }

  byteswap = iegp->byteswap;

  /* read header record */
  blocklen = binReadF77Block(fileID, byteswap);

  if ( fileEOF(fileID) ) return (-1);

  if ( IEG_Debug )
    Message("blocklen = %lu", blocklen);

  if ( blocklen == 636 || blocklen == 640 )
    dprec = 4;
  else if ( blocklen == 1040 || blocklen == 1036 )
    dprec = 8;
  else
    {
      Warning("unexpecteted header size %d!", (int) blocklen);
      return (-1);
    }

  iegp->dprec = dprec;

  binReadInt32(fileID, byteswap, 37, (INT32 *) tmpbuf);
  for ( i = 0; i < 37; i++ ) iegp->ipdb[i] = (int) ((INT32 *) tmpbuf)[i];

  binReadInt32(fileID, byteswap, 18, (INT32 *) tmpbuf);
  for ( i = 0; i < 18; i++ ) iegp->igdb[i] = (int) ((INT32 *) tmpbuf)[i];

  if ( blocklen == 636 || blocklen == 1036 )
    {
      fileRead(fileID, tmpbuf, 4);
      if ( byteswap ) swap4byte(tmpbuf, 1);
      iegp->refval = (double) ((float *) tmpbuf)[0];
    }
  else
    {
      fileRead(fileID, tmpbuf, 8);
      if ( byteswap ) swap8byte(tmpbuf, 1);
      iegp->refval = (double) ((double *) tmpbuf)[0];
    }

  binReadInt32(fileID, byteswap, 3, (INT32 *) tmpbuf);
  for ( i = 0; i < 3; i++ ) iegp->igdb[18+i] = (int) ((INT32 *) tmpbuf)[i];

  if ( dprec == SINGLE_PRECISION )
    {
      fileRead(fileID, tmpbuf, 400);
      if ( byteswap ) swap4byte(tmpbuf, 100);
      for ( i = 0; i < 100; i++ )
	iegp->vct[i] = (double) ((float *) tmpbuf)[i];
    }
  else
    {
      fileRead(fileID, tmpbuf, 800);
      if ( byteswap ) swap8byte(tmpbuf, 100);
      for ( i = 0; i < 100; i++ )
	iegp->vct[i] = (double) ((double *) tmpbuf)[i];
    }

  /*
  fprintf(stderr, "refval %g\n", iegp->refval);

  for ( i = 0; i < 100; i++ )
    fprintf(stderr, "%3d %g\n", i, iegp->vct[i]);

  {
    int i;
    for ( i = 0; i < 37; i++ )
      fprintf(stderr, "pdb: %d %d\n", i, iegp->ipdb[i]);
    for ( i = 0; i < 22; i++ )
      fprintf(stderr, "gdb: %d %d\n", i, iegp->igdb[i]);
  }
  */
  blocklen2 = binReadF77Block(fileID, byteswap);

  if ( blocklen2 != blocklen )
    {
      Warning("header blocklen differ!");
      return (-1);
    }

  iegp->datasize = (size_t)IEG_G_NumLon(iegp->igdb)
    * (size_t)IEG_G_NumLat(iegp->igdb);

  if ( IEG_Debug )
    Message("datasize = %lu", iegp->datasize);

  blocklen = binReadF77Block(fileID, byteswap);

  size_t buffersize = iegp->buffersize;

  if ( buffersize < blocklen )
    {
      buffersize = blocklen;
      buffer = iegp->buffer;
      buffer = realloc(buffer, buffersize);
      iegp->buffer = buffer;
      iegp->buffersize = buffersize;
    }
  else
    buffer = iegp->buffer;

  datasize = iegp->datasize;

  if ( dprec != (int) (blocklen/datasize) )
    {
      Warning("data precision differ! (h = %d; d = %d)",
	      (int) dprec, (int) (blocklen/datasize));
      return (-1);
    }

  fileRead(fileID, buffer, blocklen);

  blocklen2 = binReadF77Block(fileID, byteswap);

  if ( blocklen2 != blocklen )
    {
      Warning("data blocklen differ!");
      return (-1);
    }

  return (0);
}


int iegWrite(int fileID, iegrec_t *iegp)
{
  size_t datasize;
  size_t blocklen;
  size_t i;
  int dprec;
  float refvalf;
  double refval;
  char tmpbuf[800];
  float fvct[100];
  void *buffer;
  int byteswap = iegp->byteswap;


  dprec  = iegp->dprec;

  /* write header record */
  if ( dprec == SINGLE_PRECISION )
    blocklen = 636;
  else
    blocklen = 1040;

  binWriteF77Block(fileID, byteswap, blocklen);

  for ( i = 0; i < 37; i++ ) ((INT32 *) tmpbuf)[i] = (INT32) iegp->ipdb[i];
  binWriteInt32(fileID, byteswap, 37, (INT32 *) tmpbuf);

  for ( i = 0; i < 18; i++ ) ((INT32 *) tmpbuf)[i] = (INT32) iegp->igdb[i];
  binWriteInt32(fileID, byteswap, 18, (INT32 *) tmpbuf);

  refval = iegp->refval;
  refvalf = (float) refval;
  if ( dprec == SINGLE_PRECISION )
    binWriteFlt32(fileID, byteswap, 1, (FLT32 *) &refvalf);
  else
    binWriteFlt64(fileID, byteswap, 1, (FLT64 *) &refval);

  for ( i = 0; i < 3; i++ ) ((INT32 *) tmpbuf)[i] = (INT32) iegp->igdb[18+i];
  binWriteInt32(fileID, byteswap, 3, (INT32 *) tmpbuf);

  if ( dprec == SINGLE_PRECISION )
    {
      for ( i = 0; i < 100; i++ ) fvct[i] = (float) iegp->vct[i];
      binWriteFlt32(fileID, byteswap, 100, fvct);
    }
  else
    {
      binWriteFlt64(fileID, byteswap, 100, iegp->vct);
    }

  binWriteF77Block(fileID, byteswap, blocklen);

  datasize = (size_t)iegp->igdb[4] * (size_t)iegp->igdb[5];
  blocklen = datasize * (size_t)dprec;

  binWriteF77Block(fileID, byteswap, blocklen);

  iegp->datasize = datasize;

  buffer = iegp->buffer;

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
