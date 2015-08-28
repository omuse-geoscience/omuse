#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <inttypes.h>
#include <stdbool.h>
#include <stdio.h>
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
  size_t size;
  struct dBuffer *db1;
  struct dBuffer *db2;
  struct dBuffer *db;
  MPI_File fh;
  MPI_Request request;
  int fileID;
  int tsID;
  bool finished;
  char name[];
} aFiledataM;

static listSet *bibAFiledataM;

static int
fileIDTest(void *a, void *fileID)
{
  return ((aFiledataM *)a)->fileID == (int)(intptr_t)fileID;
}


/***************************************************************/

static aFiledataM *initAFiledataMPINONB ( const char *filename, size_t bs )
{
  aFiledataM *of = NULL;
  int iret;
  MPI_Comm commNode = commInqCommNode ();

  of = (aFiledataM*) xmalloc(sizeof (*of) + strlen(filename) + 1);

  strcpy(of->name, filename);
  of->size = bs;
  of->db1 = NULL;
  of->db2 = NULL;

  /* init output buffer */

  iret = dbuffer_init ( &( of->db1 ), of->size );
  iret += dbuffer_init ( &( of->db2 ), of->size );

  if ( iret > 0 ) xabort ( "dbuffer_init did not succeed" );

  of->db = of->db1;

  of->tsID = CDI_UNDEFID;

  /* open file */
  xmpi(MPI_File_open(commNode, of->name, MPI_MODE_CREATE|MPI_MODE_WRONLY,
                     MPI_INFO_NULL, &( of->fh )));
  of->request = MPI_REQUEST_NULL;
  of->finished = false;

  return of;
}

/***************************************************************/

static int
destroyAFiledataMPINONB(void *v)
{
  int iret = 0;
  aFiledataM *of;
  MPI_Status status;
  int rankNode = commInqRankNode ();
  MPI_Offset endpos;

  of = (aFiledataM * ) v;

  xdebug ( "IOPE%d: close file %d, name=\"%s\"",
           rankNode, of->fileID, of->name );

  /* close file */
  xmpi(MPI_Wait(&of->request, &status));
  xmpi(MPI_Barrier(commInqCommNode()));
  xmpi(MPI_File_get_position_shared(of->fh, &endpos));
  xmpi(MPI_File_set_size(of->fh, endpos));
  iret = MPI_File_close ( & ( of->fh ));

  /* file closed, cleanup */

  dbuffer_cleanup ( & ( of->db1 ));
  dbuffer_cleanup ( & ( of->db2 ));

  free ( of );

  xdebug ( "IOPE%d: closed file, cleaned up, return",
           rankNode );

  return iret == MPI_SUCCESS ? 0 : -1;
}

/***************************************************************/

static bool
compareNamesMPINONB(void *v1, void *v2)
{
  aFiledataM *afm1 = v1, *afm2 = v2;
  return !strcmp(afm1->name, afm2->name);
}

/***************************************************************/

static void
writeMPINONB(aFiledataM *of)
{
  int amount;
  MPI_Status status;
  int rankNode = commInqRankNode ();
  int fileID = of->fileID;

  /* write buffer */

  amount = ( int ) dbuffer_data_size ( of->db );

  if ( amount == 0 ) return;

  xdebug3 ( "IOPI%d: Write buffer, size %d bytes, in",
           rankNode, amount );

  xmpi ( MPI_Wait ( & ( of->request ), &status ));
  xmpi(MPI_File_iwrite_shared(of->fh, of->db->buffer, amount, MPI_UNSIGNED_CHAR,
                              &of->request));
  xdebug("%d bytes written for fileID=%d", amount, fileID);

  /* change outputBuffer */

  dbuffer_reset ( of->db );

  if ( of->db == of->db1 )
    {
        xdebug3 ( "IOPE%d: fileID=%d, change to buffer 2 ...",
                 rankNode, fileID );
      of->db =  of->db2;
    }
  else
    {
        xdebug3 ( "IOPE%d: fileID=%d, change to buffer 1 ...",
                  rankNode, fileID );
      of->db =  of->db1;
    }

  return;
}

/***************************************************************/

size_t fwMPINONB ( int fileID, int tsID, const void *buffer, size_t len )
{
  int error = 0;
  int filled = 0;
  aFiledataM *of;
  int rankNode = commInqRankNode ();

  of = listSetGet(bibAFiledataM, fileIDTest, (void *)(intptr_t)fileID);
  xassert(of);

  bool flush = tsID != of->tsID;

  if (flush)
    {
      xdebug3("IOPE%d: tsID = %d, flush buffer", rankNode, tsID);
      writeMPINONB(of);
      of->tsID = tsID;
      MPI_Status status;
      xmpi(MPI_Wait(&(of->request), &status));
      xmpi(MPI_Barrier(commInqCommNode()));
    }

  filled = dbuffer_push ( of->db, ( unsigned char * ) buffer, len );

  xdebug3 ( "IOPE%d: fileID = %d, tsID = %d,"
           " pushed data on buffer, filled = %d",
           rankNode, fileID, tsID, filled );

  if ( filled == 1 )
    {
      if ( flush )
        error = filled;
      else
        {
          writeMPINONB(of);

          error = dbuffer_push ( of->db, ( unsigned char * ) buffer, len );
        }
    }

  if ( error == 1 )
    xabort("did not succeed filling output buffer, fileID=%d", fileID);

  return len;
}

/***************************************************************/

int fcMPINONB ( int fileID )
{
  aFiledataM *of;
  int rankNode = commInqRankNode ();

  xdebug("IOPE%d: write buffer, close file and cleanup, in %d",
         rankNode, fileID );

  if (!(of = listSetGet(bibAFiledataM, fileIDTest, (void *)(intptr_t)fileID)))
    xabort("listSet, fileID=%d not found", fileID);

  writeMPINONB(of);

  /* remove file element */
  int iret = listSetRemove(bibAFiledataM, fileIDTest, (void *)(intptr_t)fileID);
  return iret;
}

/***************************************************************/
static void
elemCheck(void *q, void *nm)
{
  aFiledataM *afm = q;
  const char *name = nm;

  if (!strcmp(name, afm->name))
    xabort("Filename %s has already been added to set\n", name);
}


int fowMPINONB ( const char *filename )
{
  static aFiledataM *of;
  static unsigned long buffersize = 0;
  int id;
  enum {
    bcastRoot = 0
  };
  MPI_Comm commNode = commInqCommNode ();
  int rankNode = commInqRankNode ();

  /* broadcast buffersize to collectors ( just once, for all files )*/

  if (!buffersize)
    {
      if (rankNode == bcastRoot)
        buffersize = findWriteAccumBufsize();
      xmpi(MPI_Bcast(&buffersize, 1, MPI_UNSIGNED_LONG, bcastRoot, commNode));
    }

  xdebug("buffersize=%ld", buffersize);

  listSetForeach(bibAFiledataM, elemCheck, (void *)filename);
  of = initAFiledataMPINONB(filename, (size_t)buffersize);

  if ((id = listSetAdd(bibAFiledataM, of)) < 0 )
    xabort("filename %s not unique", of->name);

  xdebug("IOPE%d: name=%s, init and added aFiledataM, return id = %d",
         rankNode, filename, id);
  of->fileID = id;
  return id;
}

/***************************************************************/

void finalizeMPINONB(void)
{
  if (!listSetIsEmpty(bibAFiledataM))
    xabort("set bibAFiledataM not empty");
  else
    {
      xdebug("%s", "destroy set");
      listSetDelete(bibAFiledataM);
    }
}

/***************************************************************/

void
initMPINONB(void (*postCommSetupActions)(void))
{
  commDefCommColl ( 1 );
  commSendNodeInfo ();
  commRecvNodeMap ();
  commDefCommsIO ();
  postCommSetupActions();
  bibAFiledataM = listSetNew( destroyAFiledataMPINONB, compareNamesMPINONB );

  if ( bibAFiledataM == NULL )
    xabort ( "listSetNew did not succeed" );
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
