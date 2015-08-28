/* 
   todo 
   build in control, for consistance of pairs filename / filenumber 
   ( pioOpenFile member name, recv in tmpbuffer, if(!uniqueName(q,v,n))abort )
*/

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <inttypes.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>

#include <mpi.h>

#include "cdi.h"
#include "dmemory.h"

#include "pio.h"
#include "pio_comm.h"
#include "pio_impl.h"
#include "pio_util.h"

typedef struct
{
  struct dBuffer *db1;
  struct dBuffer *db2;
  struct dBuffer *db;
  FILE *fp;
  enum IO_Server_command command;
  int tsID, fileID;
  char name[];
} aFiledataPF;

static int
fileIDTestA(void *a, void *fileID)
{
  return ((aFiledataPF *)a)->fileID == (int)(intptr_t)fileID;
}

typedef struct
{
  long offset;
  bool finished;
  int fileID;
  bool nfinished[];
} bFiledataPF;

static int
fileIDTestB(void *a, void *fileID)
{
  return ((bFiledataPF *)a)->fileID == (int)(intptr_t)fileID;
}

static bool
fileIDCmpB(void *a, void *b)
{
  return ((bFiledataPF *)a)->fileID == ((bFiledataPF *)b)->fileID;
}

static listSet *bibAFiledataPF;

/***************************************************************/
  
static aFiledataPF *initAFiledataPF ( const char *key, size_t bs)
{
  aFiledataPF *afd;
  size_t len;
  int iret;

  len = strlen(key);
  afd = xcalloc(1, sizeof (*afd) + len + 1);
  strcpy(afd->name, key);
  afd->tsID = 0;

  /* init output buffer */

  xdebug ( " name=%s, init output buffer",  afd->name );
   
  iret = dbuffer_init(&(afd->db1), bs);
  iret += dbuffer_init(&(afd->db2), bs);

  if ( iret > 0 )
    xabort("dbuffer_init did not succeed");

  afd->db = afd->db1;

  /* open file */ 
  xdebug ( "name=%s, open file",  afd->name );

  if ( ( afd->fp = fopen ( afd->name, "w" )) == NULL ) 
    xabort("Failed to open %s", afd->name);

  afd->command = IO_Open_file;
  return afd;
}

/***************************************************************/
static bFiledataPF *
initBFiledataPF(int fileID, int nc)
{
  bFiledataPF *bfd;
  size_t bfdSize = sizeof (bFiledataPF) + (size_t)nc * sizeof (bool);
  bfd = xcalloc(1, bfdSize);
  bfd->offset = 0;
  bfd->finished = false;
  bfd->fileID = fileID;

  return bfd;
}

/***************************************************************/

static int
destroyAFiledataPF(void *v)
{
  int iret = 0;
  aFiledataPF *afd = ( aFiledataPF * ) v;

  /* close file */
  xdebug("name=%s, close file", afd->name);
  if ((iret = fclose(afd->fp)) == EOF)
    xabort("Failed to close %s", afd->name);

  /* file closed, cleanup */
  xdebug("name=%s, file closed, cleanup ...",  afd->name);
  dbuffer_cleanup(&(afd->db1));
  dbuffer_cleanup(&(afd->db2));

  free(afd);

  return iret;
}

/***************************************************************/

static int
destroyBFiledataPF(void *v)
{
  int iret = 0;
  bFiledataPF *bfd = (bFiledataPF * ) v;
  
  free ( bfd );

  return iret;
}

/***************************************************************/

static bool
compareNamesAPF(void *v1, void *v2)
{
  aFiledataPF *afd1 = v1, *afd2 = v2;

  return !strcmp(afd1->name, afd2->name);
}

/***************************************************************/

static void
fpgPOSIXFPGUARDSENDRECV(void)
{
  int i, source, iret;
  struct fileOpTag rtag;
  MPI_Status status;
  bFiledataPF *bfd; 
  listSet *bibBFiledataPF;
  long amount;
  MPI_Comm commNode = commInqCommNode ();
  int nProcsCollNode =  commInqSizeNode () - commInqSizeColl ();
  bool * sentFinalize, doFinalize = false;

  xdebug ( "ncollectors=%d on this node", nProcsCollNode );
  
  bibBFiledataPF = listSetNew( destroyBFiledataPF, fileIDCmpB);
  sentFinalize = xcalloc((size_t)nProcsCollNode, sizeof (sentFinalize[0]));

  for ( ;; )
    {
      xmpi ( MPI_Probe ( MPI_ANY_SOURCE, MPI_ANY_TAG, commNode, &status ));
      source = status.MPI_SOURCE;
      rtag = decodeFileOpTag(status.MPI_TAG);
      
      xdebug("receive message from source=%d, id=%d, command=%d ( %s )",
             source, rtag.id, rtag.command, cdiPioCmdStrTab[rtag.command]);

      switch (rtag.command)
      	{
      	case IO_Open_file:

          if (!(bfd = listSetGet(bibBFiledataPF, fileIDTestB,
                                 (void *)(intptr_t)rtag.id)))
	    {
	      bfd = initBFiledataPF(rtag.id, nProcsCollNode);

	      if ((iret = listSetAdd(bibBFiledataPF, bfd)) < 0)
		xabort("fileID=%d not unique", rtag.id);
              bfd->fileID = iret;
	    }

          xdebug("id=%d, command=%d ( %s ), send offset=%ld", rtag.id,
                 rtag.command, cdiPioCmdStrTab[rtag.command], bfd->offset);

	  xmpi ( MPI_Sendrecv ( &( bfd->offset ), 1, MPI_LONG, source,  status.MPI_TAG,
                                &amount, 1, MPI_LONG, source,  status.MPI_TAG,
                                commNode, &status ));

	  bfd->offset += amount; 
 
          xdebug("id=%d, command=%d ( %s ), recv amount=%ld, set offset=%ld",
                 rtag.id, rtag.command, cdiPioCmdStrTab[rtag.command], amount,
                 bfd->offset);

	  break;

	case IO_Set_fp:

          if (!(bfd = listSetGet(bibBFiledataPF, fileIDTestB,
                                 (void *)(intptr_t)rtag.id)))
            xabort("fileId=%d not in set", rtag.id);

          xdebug("id=%d, command=%d ( %s ), send offset=%ld", rtag.id,
                 rtag.command, cdiPioCmdStrTab[rtag.command], bfd->offset);

	  xmpi ( MPI_Sendrecv ( &( bfd->offset ), 1, MPI_LONG, source,  status.MPI_TAG,
                                &amount, 1, MPI_LONG, source,  status.MPI_TAG,
                                commNode, &status ));

	  bfd->offset += amount;

          xdebug("id=%d, command=%d ( %s ), recv amount=%ld, set offset=%ld",
                 rtag.id, rtag.command, cdiPioCmdStrTab[rtag.command], amount,
                 bfd->offset);

	  break;

	case IO_Close_file:

          if (!(bfd = listSetGet(bibBFiledataPF, fileIDTestB,
                               (void *)(intptr_t)rtag.id)))
            xabort("fileId=%d not in set", rtag.id);

          xdebug("id=%d, command=%d ( %s )), send offset=%ld", rtag.id,
                 rtag.command, cdiPioCmdStrTab[rtag.command], bfd->offset);

	  xmpi ( MPI_Sendrecv ( &( bfd->offset ), 1, MPI_LONG, source,  status.MPI_TAG,
                                &amount, 1, MPI_LONG, source,  status.MPI_TAG,
                                commNode, &status ));

	  bfd->offset += amount;

          xdebug("id=%d, command=%d ( %s ), recv amount=%ld, set offset=%ld",
                 rtag.id, rtag.command, cdiPioCmdStrTab[rtag.command], amount,
                 bfd->offset);


	  bfd->nfinished[source] = true;  
	  bfd->finished          = true;
	  
	  for ( i = 0; i < nProcsCollNode; i++ )
	    if ( !( bfd->nfinished[i] ))
	      {
		bfd->finished = false;
		break;
	      }

	  if ( bfd->finished )
            listSetRemove(bibBFiledataPF, fileIDTestB,
                          (void *)(intptr_t)rtag.id);
          break;
        case IO_Finalize:
          {  
            int buffer = CDI_UNDEFID, collID; 

            xmpi ( MPI_Recv ( &buffer, 1, MPI_INT, source, status.MPI_TAG,
                              commNode, &status ));
            sentFinalize[source] = true;
            doFinalize = true;
            for ( collID = 0; collID < nProcsCollNode; collID++ )
              doFinalize &= sentFinalize[collID];
            if ( doFinalize )
              {
                if (!listSetIsEmpty(bibBFiledataPF))
                  xabort("set bibBFiledataM not empty");
                else
                  {
                    xdebug("%s", "destroy set");
                    listSetDelete(bibBFiledataPF);
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

//*******************************************************

static void
writePF(aFiledataPF *afd)
{
  size_t amount, written;
  long offset;
  long amountL;
  int error, tag;
  MPI_Status status;
  int specialRank = commInqSpecialRankNode ();
  MPI_Comm commNode = commInqCommNode ();
  
  /* send buffersize, recv offset */

  amount = dbuffer_data_size ( afd->db );
  amountL = ( long ) amount;
  int id = afd->fileID;
  tag = encodeFileOpTag(id, afd->command);
  
  xmpi ( MPI_Sendrecv ( &amountL, 1, MPI_LONG, specialRank, tag,
                        &offset, 1, MPI_LONG, specialRank, tag,
                        commNode, &status ));
  xdebug ( "id=%d, command=%d, amount=%zu, send amountL=%ld, recv offset=%ld", 
           id, afd->command, amount, amountL, offset );
  
  /* write buffer */
  
  if (( error = fseek ( afd->fp, offset, SEEK_SET )) != 0 )
    xabort ( "did not succeed seeking fp" );

  if (( written = 
	fwrite ( afd->db->buffer, sizeof ( char ), amount, afd->fp )) !=
      amount )
    xabort("fileId=%d, expect to write %zu byte, written %zu byte",
           id, amount, written);
 
  xdebug("written %zu bytes in file %d with offset %ld",
         written, id, offset);
  
  /* change outputBuffer */
  
  dbuffer_reset ( afd->db );
  
  if ( afd->db == afd->db1 )
    {
      xdebug ( "id=%d, change to buffer 2 ...", id );
      afd->db =  afd->db2;
    }
  else 
    {
      xdebug ( "id=%d, change to buffer 1 ...", id );
      afd->db =  afd->db1;
    }
  
  afd->command = IO_Set_fp;
}


/***************************************************************/

static void
defTimestepPF(aFiledataPF *afd, int tsID)
{
  if ( afd == NULL || tsID < 0 || tsID != afd->tsID + 1 ) 
    xabort ( " defTimestepPF() didn't succeed." );
  afd->tsID = tsID;
}


/***************************************************************/

static void
flushOp(aFiledataPF *a, int tsID)
{
  writePF(a);
  defTimestepPF(a, tsID);
}


size_t fwPOSIXFPGUARDSENDRECV( int fileID, int tsID, const void *buffer, size_t len )
{
  int error = 0;
  int filled = 0;
  aFiledataPF *afd
    = listSetGet(bibAFiledataPF, fileIDTestA, (void *)(intptr_t)fileID);

  bool flush = tsID != afd->tsID;

  if (flush)
    {
      xdebug("fileID %d, tsID = %d, flush buffer", fileID, tsID);
      flushOp(afd, tsID);
      xmpi ( MPI_Barrier ( commInqCommColl ())); 
    }

  filled = dbuffer_push(afd->db, ( unsigned char * ) buffer, len);

  xdebug ( "fileID = %d, tsID = %d, pushed data on buffer, filled = %d", 
           fileID, tsID, filled ); 

  if ( filled == 1 ) 
    {
      if ( flush )
	error = filled;
      else
	{
	  writePF(afd);
     
	  error = dbuffer_push ( afd->db, ( unsigned char * ) buffer, len );
	}
    }
  
  if ( error == 1 )
    xabort("did not succeed filling output buffer, fileID=%d", fileID);
  
  return len;
}

/***************************************************************/

int fcPOSIXFPGUARDSENDRECV ( int id )
{
  aFiledataPF *afd;
  int iret;

  xdebug("write buffer, close file %d and cleanup", id);

  afd = listSetGet(bibAFiledataPF, fileIDTestA, (void *)(intptr_t)id);

  afd->command = IO_Close_file;

  writePF(afd);

  /* remove file element */
  iret = listSetRemove(bibAFiledataPF, fileIDTestA, (void *)(intptr_t)id);
  /* make sure the file is closed on all collectors before proceeding */
  xmpi(MPI_Barrier(commInqCommColl()));
  return iret;
}

/***************************************************************/
static void
elemCheck(void *q, void *nm)
{
  aFiledataPF *afd = q;
  const char *name = nm;

  if (!strcmp(name, afd->name))
    xabort("Filename %s has already been added to set\n", name);
}

int fowPOSIXFPGUARDSENDRECV ( const char *filename )
{
  int id;
  enum {
    bcastRoot = 0
  };
  aFiledataPF *afd;
  static unsigned long buffersize = 0;

  /* broadcast buffersize to collectors */
  if (!buffersize)
    {
      if (commInqRankColl() == bcastRoot)
        buffersize = findWriteAccumBufsize();
      xmpi(MPI_Bcast(&buffersize, 1, MPI_UNSIGNED_LONG, bcastRoot,
                     commInqCommColl()));
    }

  /* init and add file element */
  listSetForeach(bibAFiledataPF, elemCheck, (void *)filename);

  afd = initAFiledataPF(filename, (size_t)buffersize);

  if ((id = listSetAdd(bibAFiledataPF, afd)) < 0)
    xabort("filename %s not unique", afd->name);
  afd->fileID = id;
  xdebug("name=%s, init and add aFiledataPF, return id = %d",
         filename, id);
  {
    long offset, amount = 0L;
    int tag = encodeFileOpTag(afd->fileID, afd->command);
    int specialRank = commInqSpecialRankNode ();
    MPI_Status status;
    MPI_Comm commNode = commInqCommNode ();
    xmpi(MPI_Sendrecv(&amount, 1, MPI_LONG, specialRank, tag,
                      &offset, 1, MPI_LONG, specialRank, tag,
                      commNode, &status));
  }
  afd->command = IO_Set_fp;
  return id;
}

/***************************************************************/

void
finalizePOSIXFPGUARDSENDRECV(void)
{
  int buffer = 0, tag = encodeFileOpTag(0, IO_Finalize);

  xmpi(MPI_Send(&buffer, 1, MPI_INT, commInqSpecialRankNode(),
                tag, commInqCommNode()));

  if (!listSetIsEmpty(bibAFiledataPF))
    xabort("set bibAFiledataM not empty");
  else
    {
      xdebug("%s", "destroy set");
      listSetDelete(bibAFiledataPF);
    }
}

/***************************************************************/

void
initPOSIXFPGUARDSENDRECV(void (*postCommSetupActions)(void))
{
  if ( commInqSizeNode () < 2 ) 
    xabort ( "USAGE: # IO PROCESSES ON A PHYSICAL NODE >= 2" );

  int isCollector = commInqRankNode () != commInqSpecialRankNode ();
  commDefCommColl(isCollector);
  commSendNodeInfo ();
  commRecvNodeMap ();
  commDefCommsIO ();
  postCommSetupActions();
  if (!isCollector)
    fpgPOSIXFPGUARDSENDRECV ();
  else
    bibAFiledataPF = listSetNew( destroyAFiledataPF, compareNamesAPF );
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
