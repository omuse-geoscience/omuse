#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#include <ctype.h>
#include <yaxt.h>

#include "file.h"
#include "cdi_int.h"
#include "namespace.h"

#include "pio.h"
#include "cdi.h"
#include "cdipio.h"
#include "pio_comm.h"
#include "pio_impl.h"
#include "pio_interface.h"
#include "pio_rpc.h"
#include "pio_util.h"

const char *const cdiPioCmdStrTab[] = {
  "IO_Open_file",
  "IO_Close_file",
  "IO_Get_fp",
  "IO_Set_fp",
  "IO_Send_buffer",
  "IO_Finalize"
};

char *token = "%";

/***************************************************************/

size_t
cdiPioFileWrite(int fileID, const void *restrict buffer, size_t len, int tsID)
{
  size_t iret = 0;

  switch ( commInqIOMode ())
    {
    case PIO_MPI:
      iret = fwMPINONB ( fileID, tsID, buffer, len );
      break;
#ifndef _SX
    case PIO_ASYNCH:
#endif
    case PIO_WRITER:
      iret = pioSendWrite(fileID, tsID, buffer, len);
      break;
    case PIO_FPGUARD:
      iret = fwPOSIXFPGUARDSENDRECV ( fileID, tsID, buffer, len );
      break;
    }

  return iret;
}

/***************************************************************/

int pioFileClose ( int id )
{
  int iret = CDI_UNDEFID;
  switch ( commInqIOMode ())
    {
    case PIO_MPI:
      iret = fcMPINONB ( id );
      break;
#ifndef _SX
    case PIO_ASYNCH:
#endif
    case PIO_WRITER:
      iret = pioSendClose(id);
      break;
    case PIO_FPGUARD:
      iret = fcPOSIXFPGUARDSENDRECV ( id );
      break;
    }

  return iret;
}

/***************************************************************/

int pioFileOpen(const char *filename, const char *mode)
{
  int iret = CDI_UNDEFID;

  if ((mode[0] != 'w' && mode[0] != 'W') || mode[0] == 0 || mode[1] != 0)
    xabort("Unsupported mode \"%s\" in parallel file open.", mode);

  switch ( commInqIOMode ())
    {
    case PIO_MPI:
      iret = fowMPINONB ( filename );
      break;
#ifndef _SX
    case PIO_ASYNCH:
#endif
    case PIO_WRITER:
      iret = pioSendOpen(filename);
      break;
    case PIO_FPGUARD:
      iret = fowPOSIXFPGUARDSENDRECV ( filename );
      break;
    }

  return iret;
}

/***************************************************************/

void cdiPioFileWritingInit(void (*postCommSetupActions)(void))
{
  int IOMode = commInqIOMode ();

  commDefCommNode ();

  xassert ( IOMode != PIO_NONE  || commInqSizeNode () == 1 );

  switch ( IOMode )
    {
    case PIO_NONE:
      commDefCommColl ( 1 );
      commSendNodeInfo ();
      commRecvNodeMap ();
      commDefCommsIO ();
      break;
    case PIO_MPI:
      initMPINONB(postCommSetupActions);
      break;
#ifndef _SX
    case PIO_ASYNCH:
#endif
    case PIO_WRITER:
      pioSendInitialize(postCommSetupActions);
      break;
    case PIO_FPGUARD:
      initPOSIXFPGUARDSENDRECV(postCommSetupActions);
      break;
    }
}

/***************************************************************/

void cdiPioFileWritingFinalize(void)
{
  int IOMode = commInqIOMode ();
  switch ( IOMode )
    {
    case PIO_NONE:
      break;
    case PIO_MPI:
      finalizeMPINONB ();
      break;
#ifndef _SX
    case PIO_ASYNCH:
#endif
    case PIO_WRITER:
      pioSendFinalize();
      break;
    case PIO_FPGUARD:
      finalizePOSIXFPGUARDSENDRECV ();
      break;
    default:
      xdebug("%s", " BACKENDCLEANUP FUNCTION NOT IMPLEMENTED YET.");
    }
}

/***************************************************************/

/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
