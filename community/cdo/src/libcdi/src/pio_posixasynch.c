/* 
   todo  
   README: specialRank Pe closes down, when all output files are closed    
*/
#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif


#ifndef _SX

#include <aio.h>
#include <errno.h>
#include <fcntl.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <mpi.h>

#include "pio.h"
#include "pio_comm.h"
#include "pio_impl.h"
#include "pio_util.h"
#include "dmemory.h"

extern char *token;

typedef struct
{
  struct dBuffer *fb;
  struct aiocb *ctrlBlks;
  off_t offset;
  int currOpIndex;
  int nextOpIndex;
  int prefIndex;
  int activeCollectors;
  int handle, fileID;
  char name[];
} bFiledataPA;

static int
fileIDTest(void *a, void *fileID)
{
  return ((bFiledataPA *)a)->fileID == (int)(intptr_t)fileID;
}

static int nPrefStreams = 4;

/***************************************************************/

static bFiledataPA *
initBFiledataPA(char *filename, size_t bs, int nc)
{
  bFiledataPA *bfd;

  xdebug ( "filename=%s, buffersize=%zu, ncollectors=%d, nPrefetchStreams=%d",
           filename, bs, nc, nPrefStreams );

  bfd = (bFiledataPA*) xmalloc( sizeof (*bfd) + strlen(filename) + 1);
  strcpy(bfd->name, filename);

  if (( bfd->handle = open ( bfd->name, O_CREAT | O_WRONLY, 0666 )) == -1 )
    xabort("Failed to open %s", bfd->name);

  dbuffer_init(&(bfd->fb), (size_t)nPrefStreams * bs);

  bfd->ctrlBlks = xcalloc((size_t)nPrefStreams, sizeof (bfd->ctrlBlks[0]));

  for (int i = 0; i < nPrefStreams; i++ )
    {
      bfd->ctrlBlks[i].aio_fildes     = bfd->handle;
      bfd->ctrlBlks[i].aio_buf = bfd->fb->buffer + (size_t)i * bs;
      bfd->ctrlBlks[i].aio_reqprio    = 0;
      bfd->ctrlBlks[i].aio_sigevent.sigev_notify = SIGEV_NONE;   
    }
  
  bfd->nextOpIndex = 0;
  bfd->prefIndex = 0; 
  bfd->offset = 0;
  bfd->activeCollectors = nc;

  xdebug ( "filename=%s, opened file, return", bfd->name );

  return bfd;
}

/***************************************************************/

static int
destroyBFiledataPA ( void *v )
{
  bFiledataPA *bfd = (bFiledataPA * ) v;
  const struct aiocb *ccBP[1];
  int iret = 0;
  ssize_t ssiret;
  int nextFinishOp = (bfd->nextOpIndex - bfd->prefIndex + nPrefStreams)
    % nPrefStreams;

  xdebug ( "filename=%s, cleanup and close file", bfd->name );

  /* close file */

  for (; bfd->prefIndex > 0 ; --(bfd->prefIndex))
    {
      xdebug("file: %s, prefIndex=%d", bfd->name, (int)bfd->prefIndex);
      ccBP[0] = ( bfd->ctrlBlks + nextFinishOp );
      do
	{
	  iret = aio_suspend ( ccBP, 1, NULL );
	  if ( iret < 0 && errno != EINTR ) xabort ( "aio_suspend () failed" );
	}
      while ( iret != 0 );

      iret = aio_error(bfd->ctrlBlks + nextFinishOp);
      if (( ssiret = aio_return ( bfd->ctrlBlks + nextFinishOp )) == -1 )
	xabort("aio_return () failed: %s", strerror(iret));

      nextFinishOp = ( nextFinishOp + 1 ) % nPrefStreams;
    }

  if ((iret = ftruncate(bfd->handle, bfd->offset)) == -1)
    xabort("failed to truncate file %s: %s", bfd->name, strerror(errno));
  if (( iret = close ( bfd->handle )) == -1 )
    xabort("failed to close %s", bfd->name);

  /* file closed, cleanup */

  dbuffer_cleanup ( &( bfd->fb ));

  free(bfd->ctrlBlks);
  free(bfd);

  xdebug("%s", "closed file and cleaned up, return");

  return iret;
}

/***************************************************************/

static bool
compareNamesBPA(void *v1, void *v2)
{
  bFiledataPA *bfd1 = v1, *bfd2 = v2;

  return !strcmp(bfd1->name, bfd2->name);
}

/***************************************************************/

static void
writePA(bFiledataPA *bfd, size_t amount)
{
  const struct aiocb *ccBP[1];
  ssize_t iret;

  xdebug ( "file %s, in", bfd->name );
  
  bfd->ctrlBlks[bfd->currOpIndex].aio_nbytes = amount;
  bfd->ctrlBlks[bfd->currOpIndex].aio_offset = bfd->offset;

  xdebug ( " before aio_write(), file %s, aio_nbytes=%zu, aio_offset=%zu",
           bfd->name, bfd->ctrlBlks[bfd->currOpIndex].aio_nbytes,
           bfd->ctrlBlks[bfd->currOpIndex].aio_offset );

  iret = aio_write ( bfd->ctrlBlks + bfd->currOpIndex );

  xdebug ( "after aio_write(), file %s, aio_nbytes=%zu, aio_offset=%zu,"
           "iret=aio_write()=%d",
           bfd->name, bfd->ctrlBlks[bfd->currOpIndex].aio_nbytes,
           bfd->ctrlBlks[bfd->currOpIndex].aio_offset, (int)iret );
   
  if ( iret == -1 ) 
    {
      xabort ( "did not succeed writing buffer" );
    }
  else
    xdebug ( "buffer written to %s",  bfd->name );
     
  bfd->offset += (off_t)amount;
  bfd->prefIndex ++;

  if ( bfd->prefIndex >= nPrefStreams ) 
    {
      ccBP[0] = ( bfd->ctrlBlks + bfd->nextOpIndex );
      do
	{
	  iret = aio_suspend ( ccBP, 1, NULL );
	  if ( iret < 0 && errno != EINTR )
	    xabort ( "aio_suspend () failed" );
	} while ( iret != 0 );

      if (( iret = aio_return ( bfd->ctrlBlks + bfd->nextOpIndex )) == -1 )
	xabort ( "aio_return () failed" );

      bfd->prefIndex --;
    }

  xdebug ( "filename=%s, prefIndex=%d, return", bfd->name, bfd->prefIndex );
}

/***************************************************************/
static void
elemCheck(void *q, void *nm)
{
  bFiledataPA *bfd = q;
  const char *name = nm;

  if (!strcmp(name, bfd->name))
    xabort("Filename %s has already been inserted\n", name);
}

/***************************************************************/

void pioWriterAIO(void)
{
  bFiledataPA *bfd; 
  listSet * bibBFiledataPA;
  MPI_Comm commNode = commInqCommNode ();
  int nProcsCollNode = commInqSizeNode () - commInqSizeColl ();
  bool * sentFinalize, doFinalize;

  if ( nPrefStreams < 1 ) xabort("USAGE: # PREFETCH STREAMS >= 1");
  xdebug ( "nProcsCollNode=%d on this node", nProcsCollNode );
 
  bibBFiledataPA = listSetNew(destroyBFiledataPA, compareNamesBPA);
  sentFinalize = xcalloc((size_t)nProcsCollNode, sizeof (sentFinalize[0]));
  
  for ( ;; )
    {
      MPI_Status status;
      xmpiStat(MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, commNode,
                         &status ), &status);

      int source = status.MPI_SOURCE;
      int tag = status.MPI_TAG;
      struct fileOpTag rtag = decodeFileOpTag(tag);

      int messagesize;
      xmpi(MPI_Get_count(&status, MPI_UNSIGNED_CHAR, &messagesize));

      xdebug ( "receive message from source=%d, id=%d, command=%d ( %s ), "
               "messagesize=%d", source, rtag.id, rtag.command,
               cdiPioCmdStrTab[rtag.command], messagesize);

      switch (rtag.command)
	{
      	case IO_Open_file:
          {
            char *messageBuffer = xmalloc((size_t)messagesize
                                          * sizeof (messageBuffer[0]));
            char *pMB = messageBuffer;

            xmpi(MPI_Recv(messageBuffer, messagesize, MPI_UNSIGNED_CHAR, source,
                          tag, commNode, &status ));

            char *filename = strtok(pMB, token);
            pMB += (strlen(filename) + 1);
            char *temp = strtok(pMB, token);
            long buffersize =  strtol(temp, NULL, 16);
            pMB += (strlen(temp) + 1);
            size_t amount = (size_t)(messageBuffer + messagesize - pMB);

            xdebug("command  %s, filename=%s, buffersize=%ld, amount=%zd",
                   cdiPioCmdStrTab[rtag.command], filename, buffersize, amount);

            if (!(bfd = listSetGet(bibBFiledataPA, fileIDTest,
                                   (void *)(intptr_t)rtag.id)))
              {
                listSetForeach(bibBFiledataPA, elemCheck, filename);
                bfd = initBFiledataPA(filename, (size_t)buffersize, nProcsCollNode);
                int id;
                if ((id = listSetAdd(bibBFiledataPA, bfd)) < 0)
                  xabort("fileID=%d not unique", rtag.id);
                bfd->fileID = id;
              }
            else
              if (strcmp(filename, bfd->name) != 0)
                xabort("filename is not consistent, fileID=%d", rtag.id);

            bfd->currOpIndex = bfd->nextOpIndex;
            bfd->nextOpIndex = ( bfd->nextOpIndex + 1 ) % nPrefStreams;

            memcpy((void *)bfd->ctrlBlks[bfd->currOpIndex].aio_buf,
                   pMB, (size_t)amount);

            writePA(bfd, amount);

            free(messageBuffer);
          }
	  break;

	case IO_Send_buffer:

          if (!(bfd = listSetGet(bibBFiledataPA, fileIDTest,
                               (void *)(intptr_t)rtag.id)))
            xabort("fileID=%d is not in set", rtag.id);

	  xdebug("command: %s, id=%d, name=%s",
                 cdiPioCmdStrTab[rtag.command], rtag.id, bfd->name );

	  bfd->currOpIndex = bfd->nextOpIndex;
	  bfd->nextOpIndex = ( bfd->nextOpIndex + 1 ) % nPrefStreams;

	  xmpi(MPI_Recv((void *)bfd->ctrlBlks[bfd->currOpIndex].aio_buf,
                        messagesize, MPI_UNSIGNED_CHAR, source, tag, commNode,
                        &status));
	  writePA(bfd, (size_t)messagesize);
	  break;

	case IO_Close_file:

          if (!(bfd = listSetGet(bibBFiledataPA, fileIDTest,
                               (void *)(intptr_t)rtag.id)))
            xabort("fileID=%d is not in set", rtag.id);

	  xdebug(" command %s, id=%d, name=%s",
                 cdiPioCmdStrTab[rtag.command], rtag.id, bfd->name);

	  bfd->currOpIndex = bfd->nextOpIndex;

	  bfd->nextOpIndex = ( bfd->nextOpIndex + 1 ) % nPrefStreams;

	  xmpi(MPI_Recv((void *)bfd->ctrlBlks[bfd->currOpIndex].aio_buf,
                        messagesize, MPI_UNSIGNED_CHAR,
                        source, tag, commNode, &status));

	  writePA(bfd, (size_t)messagesize);

	  if ( ! --(bfd->activeCollectors))
	    {
              xdebug ( "all are finished with file %d, delete node", rtag.id);
              listSetRemove(bibBFiledataPA, fileIDTest,
                            (void *)(intptr_t)rtag.id);
	    }
          break;
        case IO_Finalize:
          {
            int buffer, collID;

            xmpi ( MPI_Recv ( &buffer, 1, MPI_INT, source, tag, commNode, &status ));
            sentFinalize[source] = true;
            doFinalize = true;
            for ( collID = 0; collID < nProcsCollNode; collID++ )
              doFinalize &= sentFinalize[collID];
            if ( doFinalize )
              {
                if (!listSetIsEmpty(bibBFiledataPA))
                  xabort("Set bibBfiledataP is not empty.");
                else
                  {
                    xdebug("%s", "all files are finished, destroy set,"
                           " return");
                    listSetDelete(bibBFiledataPA);
                  }
                free(sentFinalize);
                return;
              }
          }

          break;
        default:
          xabort ( "COMMAND NOT IMPLEMENTED" );
	}
    }
}



/***************************************************************/


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
