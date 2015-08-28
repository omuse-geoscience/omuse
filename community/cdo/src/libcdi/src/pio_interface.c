#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

#include <mpi.h>
#include <yaxt.h>

#include "cdi.h"
#include "cdipio.h"
#include "cdi_int.h"
#include "dmemory.h"
#include "namespace.h"
#include "pio.h"
#include "pio_client.h"
#include "pio_serialize.h"
#include "pio_interface.h"
#include "pio_comm.h"
#include "pio_rpc.h"
#include "pio_server.h"
#include "pio_util.h"
#include "resource_handle.h"
#include "vlist.h"
#include "vlist_var.h"


static struct rdmaWin
{
  size_t size;
  unsigned char *buffer, *head;
  MPI_Win win;
  int postSet, refuseFuncCall;
  MPI_Group ioGroup;
  int dictSize, dictDataUsed, dictRPCUsed;
} *txWin = NULL;


const char * const funcMap[numRPCFuncs] = {
  "streamOpen",
  "streamDefVlist",
  "streamClose",
  "streamDefTimestep",
};

float cdiPIOpartInflate_;

/****************************************************/

static int cmp ( const void * va, const void * vb )
{
    const int ** a, ** b;

  a = ( const int ** ) va;
  b = ( const int ** ) vb;

  return (( **a < **b ) - ( **a > **b ));
}

void memcpyPackFunc(void *dataDesc, void *buf, int size, int *pos,
                    void *context)
{
  struct memCpyDataDesc *p = dataDesc;
  (void)context;
  xassert(size >= *pos && (size_t)(size - *pos) >= p->obj_size);
  memcpy((unsigned char *)buf + *pos, p->obj, p->obj_size);
  *pos += (int)p->obj_size;
}

/****************************************************/

static void
mapProblems(int problemSizes[], int * problemMapping, int nProblems,
            int nWriter, double * w)
{
  int *ip[nProblems];
  int dummy[nProblems];
  int buckets[nWriter];
  int currCapacity, nextCapacity;
  double meanBucket[nWriter];
  int sum = 0;
  int writerIdx = 0;
  int i, j;


  for ( i = 0; i < nProblems; i++ )
    {
      ip[i] = &problemSizes[i];
      sum += problemSizes[i];
    }

  qsort(ip, (size_t)nProblems, sizeof ( int * ), cmp );

  for ( i = 0; i < nProblems; i++ )
    dummy[i] = (int)(ip[i] - problemSizes);

  for ( j = 0; j < nWriter; j++ )
    meanBucket[j] = ( double ) sum * ( * ( w + j ));

  memset(buckets, 0, sizeof (buckets));

  for ( i = 0; i < nProblems; i++ )
    {
      currCapacity = INT_MIN;

      for ( j = 0; j < nWriter; j++ )
	{
	  nextCapacity = (int)meanBucket[j] - ( buckets[j] + ( *ip[i] ));

	  if ( nextCapacity > currCapacity )
	    {
	      currCapacity = nextCapacity;
	      writerIdx = j;
	    }
	}
      problemMapping[ dummy[i] ] = writerIdx;
      buckets[writerIdx] +=  *ip[i];
    }

  if ( ddebug )
    {
      xprintArray3 (  "problemSizes = ", problemSizes, nProblems, DATATYPE_INT );
      xprintArray3 ( "vector of indices, qsort of problemSizes", dummy,
                    nProblems, DATATYPE_INT );
      xprintArray3 ( "problemMapping", problemMapping, nProblems, DATATYPE_INT );
      xprintArray3 ( "meanBucket", meanBucket, nWriter, DATATYPE_FLT );
      xprintArray3 ( "actual buckets", buckets, nWriter, DATATYPE_INT );
    }
}

/****************************************************/

/**
   @brief is encapsulated in CDI library.

   @param vSizes array with number of levels for all var_t 's in order of streams

   @param sSizes array with number of var_t for all stream_t 's

   @param varMapping return value, array with ranks of I/O PEs assigned to var_t 's
in order of vSizes

   @param nStreams number of stream_t 's

   @param nodeSizes array of number of I/O PEs on each physical nodes, increasing
   order of ranks assumed

   @param nNodes number of physical nodes hosting I/O PEs

   @return
*/

static void
varMapGen(int *vSizes, int *sSizes, int *varMapping,
          unsigned nStreams, int *nodeSizes, int nNodes)
{

  int weightsStreams[nStreams];
  int streamMapping[nStreams];
  int nPEs = 0, nVars = 0;
  int i, j, k, offset = 0, offsetN = 0;

  int * weightsVarsNode;
  int * varMappingNode;
  int nVarsNode, summandRank = 0;
  int nProcsColl = commInqNProcsColl ();

  int buckets[nProcsColl];

  for (unsigned i = 0; i < nStreams; i++ )
    {
      nVars += sSizes[i];
      weightsStreams[i] = 0;
      for ( j = 0; j < * ( sSizes + i ); j++ )
	weightsStreams[i] += * ( vSizes + offset++ );
    }

  double *w = (double *)xmalloc((size_t)nNodes * sizeof ( double ));
  for ( j = 0; j < nNodes; j++ )
    nPEs += * ( nodeSizes + j );

  for ( j = 0; j < nNodes; j++ )
    w[j] = ( double ) * ( nodeSizes + j ) / ( double ) nPEs;

  mapProblems(weightsStreams, streamMapping, (int)nStreams, nNodes, w);
  free ( w );

  for ( i = 0; i < nNodes; i++ )
    {
      nVarsNode = 0;
      for (unsigned j = 0; j < nStreams; j++)
	if ( * ( streamMapping + j ) == i )
	  nVarsNode += * ( sSizes + j );

      weightsVarsNode = xmalloc((size_t)nVarsNode * sizeof (int));
      varMappingNode = xmalloc((size_t)nVarsNode * sizeof ( int ));
      w = xmalloc((size_t)nodeSizes[i] * sizeof (double));
      offset = 0;
      offsetN = 0;

      for (unsigned j = 0; j < nStreams; j++ )
	if (streamMapping[j] == i)
	  for ( k = 0; k < * ( sSizes + j ); k ++ )
	    weightsVarsNode[offsetN++] = vSizes[offset++];
	else
	  offset += sSizes[j];

      for ( j = 0; j < * ( nodeSizes + i ); j ++ )
	w[j] = 1.0 / ( double ) * ( nodeSizes + i );

      mapProblems ( weightsVarsNode, varMappingNode, nVarsNode,
		    * ( nodeSizes + i ),  w );

      offset = 0;
      offsetN = 0;

      for (unsigned j = 0; j < nStreams; j++)
	if ( * ( streamMapping + j ) == i )
	  for ( k = 0; k < * ( sSizes + j ); k ++ )
	    * ( varMapping + offset ++ ) =
              commCollID2RankGlob ( * ( varMappingNode + offsetN ++ ) +
                                    summandRank );
	else
	  offset += * ( sSizes + j );

      summandRank += * ( nodeSizes + i );

      free ( w );
      free ( varMappingNode );
      free ( weightsVarsNode );
    }

  if ( ddebug )
    {
      xprintArray ( "varMapping", varMapping, nVars, DATATYPE_INT  );
      for ( i = 0; i < nProcsColl; i++ )
	buckets[i] = 0;
      for ( i = 0; i < nVars; i ++ )
	buckets[commRankGlob2CollID ( *(varMapping + i ))] += * ( vSizes + i );
      xprintArray ( "buckets", buckets, nProcsColl, DATATYPE_INT );
    }
}

/************************************************************************/

static void
varsMapNDeco(int nNodes, int *nodeSizes)
{
  int nVars, * resHs, * streamSizes, * varSizes, * varMapping,
    * collectsData;
  int k = 0;
  int nProcsColl = commInqNProcsColl ();

  xdebug ( "START, nProcsColl=%d", nProcsColl );

  unsigned nStreams = reshCountType(&streamOps);

  resHs       = xmalloc(nStreams * sizeof (resHs[0]));
  streamSizes = xmalloc(nStreams * sizeof (streamSizes[0]));
  collectsData = xmalloc((size_t)nProcsColl * sizeof (collectsData[0]));
  cdiStreamGetIndexList(nStreams, resHs);

  for (unsigned i = 0; i < nStreams; i++ )
    streamSizes[i] = streamInqNvars(resHs[i]);

  nVars = sum_int((size_t)nStreams, streamSizes);
  varSizes   = xcalloc((size_t)nVars, sizeof (varSizes[0]));
  varMapping = xmalloc((size_t)nVars * sizeof (varMapping[0]));

  for (unsigned i = 0; i < nStreams; i++ )
    {
      int vlistID = streamInqVlist(resHs[i]);
      for (int j = 0; j < streamSizes[i]; j++ )
        varSizes[k++] += vlistInqVarSize(vlistID, j);
    }

  xassert ( k == nVars );

  varMapGen ( varSizes, streamSizes, varMapping,
	      nStreams, nodeSizes, nNodes );

  k = 0;
  for (unsigned i = 0; i < nStreams; i++ )
    for (int j = 0; j < * ( streamSizes + i ); j++ )
      {
        vlistDefVarIOrank ( streamInqVlist ( * ( resHs + i )), j,
                            * ( varMapping + k ));
        vlistDefVarIOrank ( streamInqVlistIDorig ( * ( resHs + i )), j,
                            * ( varMapping + k ));
        collectsData[commRankGlob2CollID ( varMapping[k++] )] = 1;
      }

  for (int j = 0; j < nProcsColl; j++ )
    if ( collectsData[j] == 0 )
      xabort("AT LEAST ONE COLLECTOR PROCESS IDLES, "
             "CURRENTLY NOT COVERED: "
             "PE%d collects no data",
             commCollID2RankGlob(j));

  if ( varMapping )   free ( varMapping );
  if ( varSizes )     free ( varSizes );
  if ( collectsData ) free ( collectsData );
  if ( streamSizes )  free ( streamSizes );
  if ( resHs )        free ( resHs );

  xdebug("%s", "RETURN");
}

/************************************************************************/

static
void modelWinCleanup ( void )
{
  int collID;

  xdebug("%s", "START");
  if (txWin != NULL)
    {
      for ( collID = 0; collID < commInqNProcsColl (); collID ++ )
        {
          if (txWin[collID].postSet)
            xmpi(MPI_Win_wait(txWin[collID].win));
          xmpi(MPI_Win_free(&txWin[collID].win));
          xmpi ( MPI_Free_mem ( txWin[collID].buffer ));
          xmpi(MPI_Group_free(&txWin[collID].ioGroup));
        }
      free(txWin);
    }

  xdebug("%s", "RETURN. CLEANED UP MPI_WIN'S");
}

/************************************************************************/

struct collDesc
{
  int numDataRecords, numRPCRecords;
};

static void
modelWinDefBufferSizes(void)
{
  int collID, * streamIndexList, nvars, varID;
  size_t sumWinBufferSize = 0;
  int nProcsColl  = commInqNProcsColl ();
  int rankGlob    = commInqRankGlob ();
  int root = commInqRootGlob ();
  struct collDesc *collIndex;

  xdebug("%s", "START");
  xassert(txWin != NULL);

  unsigned nstreams = reshCountType ( &streamOps );
  streamIndexList = xmalloc((size_t)nstreams * sizeof (streamIndexList[0]));
  collIndex = xcalloc((size_t)nProcsColl, sizeof (collIndex[0]));
  reshGetResHListOfType ( nstreams, streamIndexList, &streamOps );
  for (unsigned streamNo = 0; streamNo < nstreams; streamNo++ )
    {
      // memory required for data
      int streamID = streamIndexList[streamNo];
      int vlistID = streamInqVlist(streamID);
      nvars = vlistNvars ( vlistID );
      for ( varID = 0; varID < nvars; varID++ )
        {
          int collID = commRankGlob2CollID(vlistInqVarIOrank(vlistID, varID));
          size_t collIDchunk;
          {
            int varSize = vlistInqVarSize(vlistID, varID);
            int nProcsModel = commInqNProcsModel();
            collIDchunk = (size_t)ceilf(cdiPIOpartInflate_
                                        * (float)(varSize + nProcsModel - 1)
                                        / (float)nProcsModel);
          }
          xassert ( collID != CDI_UNDEFID && collIDchunk > 0 );
          collIndex[collID].numDataRecords += 2;
          txWin[collID].size += (size_t)collIDchunk * sizeof (double)
            /* re-align chunks to multiple of double size */
            + sizeof (double) - 1
            /* one header for data record, one for corresponding part
             * descriptor*/
            + 2 * sizeof (struct winHeaderEntry)
            /* FIXME: heuristic for size of packed Xt_idxlist */
            + sizeof (Xt_int) * collIDchunk * 3;
        }

      // memory required for the function calls encoded
      // for remote execution
      // once per stream and timestep for all collprocs only on the modelproc root
      if ( rankGlob == root )
        for ( collID = 0; collID < nProcsColl; collID++ )
          {
            collIndex[collID].numRPCRecords += numRPCFuncs;
            txWin[collID].size +=
              numRPCFuncs * sizeof (struct winHeaderEntry)
              + MAXDATAFILENAME
              /* data part of streamDefTimestep */
              + (2 * CDI_MAX_NAME + sizeof (taxis_t));
          }
    }
  for (collID = 0; collID < nProcsColl; ++collID)
    {
      int numRecords = 1 + collIndex[collID].numDataRecords
        + collIndex[collID].numRPCRecords;
      txWin[collID].dictSize = numRecords;
      txWin[collID].dictDataUsed = 1;
      txWin[collID].dictRPCUsed = 0;
      /* account for size header */
      txWin[collID].size += sizeof (struct winHeaderEntry);
      txWin[collID].size = roundUpToMultiple(txWin[collID].size,
                                             PIO_WIN_ALIGN);
      sumWinBufferSize += (size_t)txWin[collID].size;
    }
  free(collIndex);
  free ( streamIndexList );

  xdebug("sumWinBufferSize=%zu, MAXWINBUFFERSIZE=%zu", sumWinBufferSize,
         (size_t)MAXWINBUFFERSIZE);
  xassert ( sumWinBufferSize <= (size_t)MAXWINBUFFERSIZE );
  xdebug("%s", "RETURN");
}


/************************************************************************/


static
  void modelWinFlushBuffer ( int collID )
{
  int nProcsColl = commInqNProcsColl ();

  xassert ( collID                >= 0         &&
            collID                < nProcsColl &&
            txWin != NULL      &&
            txWin[collID].buffer     != NULL      &&
            txWin[collID].size <= MAXWINBUFFERSIZE);
  memset(txWin[collID].buffer, 0, txWin[collID].size);
  txWin[collID].head = txWin[collID].buffer
    + (size_t)txWin[collID].dictSize * sizeof (struct winHeaderEntry);
  txWin[collID].refuseFuncCall = 0;
  txWin[collID].dictDataUsed = 1;
  txWin[collID].dictRPCUsed = 0;
}


/************************************************************************/


static
void modelWinCreate ( void )
{
  int collID, ranks[1];
  int nProcsColl = commInqNProcsColl ();

  xdebug("%s", "START");
  txWin = xcalloc((size_t)nProcsColl, sizeof (txWin[0]));

  modelWinDefBufferSizes ();
  ranks[0] = commInqNProcsModel ();

  MPI_Info no_locks_info;
  xmpi(MPI_Info_create(&no_locks_info));
  xmpi(MPI_Info_set(no_locks_info, "no_locks", "true"));
  for ( collID = 0; collID < nProcsColl; collID ++ )
    {
      xassert(txWin[collID].size > 0);
      txWin[collID].buffer = NULL;
      xmpi(MPI_Alloc_mem((MPI_Aint)txWin[collID].size, MPI_INFO_NULL,
                         &txWin[collID].buffer));
      xassert ( txWin[collID].buffer != NULL );
      txWin[collID].head = txWin[collID].buffer
        + (size_t)txWin[collID].dictSize * sizeof (struct winHeaderEntry);
      xmpi(MPI_Win_create(txWin[collID].buffer, (MPI_Aint)txWin[collID].size, 1,
                          no_locks_info, commInqCommsIO(collID),
                          &txWin[collID].win));
      MPI_Group commGroup;
      xmpi(MPI_Comm_group(commInqCommsIO(collID), &commGroup));
      xmpi(MPI_Group_incl(commGroup, 1, ranks, &txWin[collID].ioGroup));
      xmpi(MPI_Group_free(&commGroup));
    }

  xmpi(MPI_Info_free(&no_locks_info));

  xdebug("%s", "RETURN, CREATED MPI_WIN'S");
}

/************************************************************************/

static void
modelWinEnqueue(int collID,
                struct winHeaderEntry header, const void *data,
                valPackFunc packFunc)
{
  struct winHeaderEntry *winDict
    = (struct winHeaderEntry *)txWin[collID].buffer;
  int targetEntry;
  if (header.id > 0 || header.id == PARTDESCMARKER)
    targetEntry = (txWin[collID].dictDataUsed)++;
  else
    targetEntry = txWin[collID].dictSize - ++(txWin[collID].dictRPCUsed);
  if (header.id > 0)
    {
      int offset = header.offset
        = (int)roundUpToMultiple((size_t)(txWin[collID].head
                                          - txWin[collID].buffer),
                                 sizeof (double));
      MPI_Comm comm = commInqCommsIO(collID);
      packFunc((void *)data, txWin[collID].buffer, (int)txWin[collID].size,
               &offset, &comm);
      txWin[collID].head = txWin[collID].buffer + offset;
    }
  else if (header.id == PARTDESCMARKER)
    {
      Xt_uid uid = header.specific.partDesc.uid;
      int offset = -1;
      /* search if same uid entry has already been enqueued */
      for (int entry = 2; entry < targetEntry; entry += 2)
        {
          xassert(winDict[entry].id == PARTDESCMARKER);
          if (winDict[entry].specific.partDesc.uid == uid)
            {
              offset = winDict[entry].offset;
              break;
            }
        }
      if (offset == -1)
        {
          /* not yet used partition descriptor, serialize at
           * current position */
          int position = header.offset
            = (int)(txWin[collID].head - txWin[collID].buffer);
          MPI_Comm comm = commInqCommsIO(collID);
          packFunc((void *)data, txWin[collID].buffer, (int)txWin[collID].size,
                   &position, &comm);
          txWin[collID].head = txWin[collID].buffer + position;
        }
      else
        /* duplicate entries are copied only once per timestep */
        header.offset = offset;
    }
  else
    {
      int position = header.offset
        = (int)(txWin[collID].head - txWin[collID].buffer);
      MPI_Comm comm = commInqCommsIO(collID);
      packFunc((void *)data, txWin[collID].buffer, (int)txWin[collID].size,
               &position, &comm);
      txWin[collID].head = txWin[collID].buffer + position;
    }
  winDict[targetEntry] = header;
}

static void
cdiPio_xt_idxlist_pack_wrap(void *data, void *buf, int size, int *pos,
                            void *context)
{
  MPI_Comm comm = *(MPI_Comm *)context;
  size_t pack_size = xt_idxlist_get_pack_size((Xt_idxlist)data, comm);
  xassert(size >= *pos && pack_size <= (size_t)(size - *pos));
  xt_idxlist_pack((Xt_idxlist)data, (unsigned char *)buf,
                  size, pos, comm);
}

static inline void
collWait(int collID)
{
  if (txWin[collID].postSet)
    {
      xmpi(MPI_Win_wait(txWin[collID].win));
      txWin[collID].postSet = 0;
      modelWinFlushBuffer(collID);
    }
}

static inline void
collProbe(int collID)
{
  if (txWin[collID].postSet)
    {
      int flag;
      xmpi(MPI_Win_test(txWin[collID].win, &flag));
      if (flag)
        {
          txWin[collID].postSet = 0;
          modelWinFlushBuffer(collID);
        }
    }
}

void
cdiPioRDMAProgress()
{
  int nProcsColl = commInqNProcsColl();
  for (int collID = 0; collID < nProcsColl; collID++)
    collProbe(collID);
}


static void
pioBufferPartData_(int streamID, int varID,
                   const void *packData, valPackFunc packDataFunc,
                   int nmiss, Xt_idxlist partDesc)
{
  int vlistID, collID = CDI_UNDEFID;

  vlistID  = streamInqVlist ( streamID );
  collID   = commRankGlob2CollID ( vlistInqVarIOrank ( vlistID, varID ));
  xassert ( collID         >= 0                    &&
            collID         <  commInqNProcsColl () &&
            txWin != NULL);

  collWait(collID);


  struct winHeaderEntry dataHeader
    = { .id = streamID, .specific.dataRecord = { varID, nmiss }, .offset = -1 };
  modelWinEnqueue(collID, dataHeader, packData, packDataFunc);
  {
    struct winHeaderEntry partHeader
      = { .id = PARTDESCMARKER,
          .specific.partDesc = { .uid = xt_idxlist_get_uid(partDesc) },
          .offset = 0 };
    modelWinEnqueue(collID, partHeader, partDesc, cdiPio_xt_idxlist_pack_wrap);
  }

  txWin[collID].refuseFuncCall = 1;
}

void
pioBufferPartData(int streamID, int varID, const double *data,
                  int nmiss, Xt_idxlist partDesc)
{
  int chunk = xt_idxlist_get_num_indices(partDesc);
  xassert(chunk <= INT_MAX);
  pioBufferPartData_(streamID, varID,
                     &(struct memCpyDataDesc){data, (size_t)chunk * sizeof (data[0])},
                     memcpyPackFunc,
                     nmiss, partDesc);
}

struct scatterGatherDesc
{
  void *data;
  const int *blocklengths, *displacements;
  size_t elemSize;
  unsigned numBlocks;
  unsigned numElems;
};

static void
scatterGatherPackFunc(void *dataDesc, void *buf, int size, int *pos,
                      void *context)
{
  (void)context;
  const struct scatterGatherDesc *p = dataDesc;
  unsigned numBlocks = p->numBlocks;
  const int *bls = p->blocklengths, *disps = p->displacements;
  int pos_ = *pos;
  unsigned char *dstBuf = buf + pos_, *bufEnd = (unsigned char *)buf + size;
  size_t elemSize = p->elemSize;
  xassert(elemSize <= SSIZE_MAX);
  const unsigned char *data = p->data;
  unsigned copyCount = 0, numElems = p->numElems;
  for (unsigned j = 0; j < numBlocks && copyCount < numElems; ++j)
    {
      int bl = bls[j];
      if (bl > 0)
        {
          if ((unsigned)bl + copyCount > numElems)
            {
              bl = (int)(numElems - copyCount);
              Warning("%s: %s", "streamWriteScatteredVarPart",
                      "blocks longer than number of elements in index list!");
            }
          size_t bsize = (size_t)bl * elemSize;
          xassert(dstBuf + bsize <= bufEnd);
          memcpy(dstBuf, data + (ssize_t)elemSize * (ssize_t)disps[j], bsize);
          dstBuf += bsize;
        }
    }
  *pos = (int)(dstBuf - (unsigned char *)buf);
}


void
cdiPioBufferPartDataGather(int streamID, int varID, const double *data,
                           int numBlocks, const int blocklengths[],
                           const int displacements[],
                           int nmiss, Xt_idxlist partDesc)
{
  xassert(numBlocks >= 0);
  pioBufferPartData_(streamID, varID,
                     &(struct scatterGatherDesc)
                     { .data = (void *)data, .blocklengths = blocklengths,
                       .displacements = displacements,
                       .elemSize = sizeof (data[0]),
                       .numBlocks = (unsigned)numBlocks,
                       .numElems
                         = (unsigned)xt_idxlist_get_num_indices(partDesc) },
                     scatterGatherPackFunc,
                     nmiss, partDesc);
}


/************************************************************************/

void pioBufferFuncCall(struct winHeaderEntry header,
                       const void *data, valPackFunc dataPackFunc)
{
  int rankGlob = commInqRankGlob ();
  int root = commInqRootGlob ();
  int collID, nProcsColl = commInqNProcsColl ();
  int funcID = header.id;

  xassert(funcID >= MINFUNCID && funcID <= MAXFUNCID);
  xdebug("%s, func: %s", "START", funcMap[(-1 - funcID)]);

  if ( rankGlob != root ) return;

  xassert(txWin != NULL);

  for (collID = 0; collID < nProcsColl; ++collID)
    {
      collWait(collID);
      xassert(txWin[collID].dictRPCUsed + txWin[collID].dictDataUsed
              < txWin[collID].dictSize);
      xassert(txWin[collID].refuseFuncCall == 0);
      modelWinEnqueue(collID, header, data, dataPackFunc);
    }

  xdebug("%s", "RETURN");
}


void
cdiPioNoPostCommSetup(void)
{
}

/*****************************************************************************/

/* pioInit definition must currently compile even in non-MPI configurations */
/**
   @brief initializes the MPI_Communicators needed for the
  communication between the calculator PEs and the I/O PEs and within the
  group of I/O PEs.

  commGlob: all PEs

  commPIO: I/O PEs, PEs with highest ranks in commGlob

  commModel: calculating PEs, no I/O PEs

  commsIO[i]:

  Collective call

  @param comm MPI_Communicator of all calling PEs
  @param nIOP number of I/O PEs
  @param partInflate allow for array partitions on comute
  PE that are at most sized \f$ partInflate * \lceil arraySize /
  numComputePEs \rceil \f$
  @param postSetupActions function which is called by all I/O servers
  after communicator split
  @return int indicating wether the calling PE is a calcutator (1) or not (0)
*/

static int pioNamespace_ = -1;
static int xtInitByCDI = 0;

MPI_Comm
pioInit(MPI_Comm commGlob, int nProcsIO, int IOMode,
        int *pioNamespace, float partInflate,
        void (*postCommSetupActions)(void))
{
  int sizeGlob;

  namespaceSwitchSet(NSSWITCH_WARNING, NSSW_FUNC(cdiPioWarning));

  if ( IOMode < PIO_MINIOMODE || IOMode > PIO_MAXIOMODE )
    xabort ( "IOMODE IS NOT VALID." );

#ifdef _SX
  if ( IOMode ==  PIO_ASYNCH )
    xabort ( "PIO_ASYNCH DOES NOT WORK ON SX." );
#endif

  if ((xtInitByCDI = (!xt_initialized() || xt_finalized())))
    xt_initialize(commGlob);
  commInit ();
  commDefCommGlob ( commGlob );
  sizeGlob = commInqSizeGlob ();

  if (((IOMode != PIO_NONE && (nProcsIO <= 0 || nProcsIO > sizeGlob - 1)))
      || (IOMode == PIO_NONE && nProcsIO != 1))
    xabort("DISTRIBUTION OF TASKS ON PROCS IS NOT VALID.\n"
           "nProcsIO=%d, sizeGlob=%d\n", nProcsIO, sizeGlob);

  commDefNProcsIO ( nProcsIO );
  commDefIOMode   ( IOMode );
  commDefCommPio  ();

  xassert(partInflate >= 1.0);
  cdiPIOpartInflate_ = partInflate;

  // JUST FOR TEST CASES WITH ONLY ONE MPI TASK
  if ( commInqSizeGlob () == 1 )
    {
      pioNamespace_ = *pioNamespace = namespaceNew();
      return commInqCommGlob ();
    }

  if ( commInqIsProcIO ())
    {
      cdiPioSerializeSetMPI();
      namespaceSwitchSet(NSSWITCH_ABORT, NSSW_FUNC(cdiAbortC_MPI));
      namespaceSwitchSet(NSSWITCH_FILE_OPEN, NSSW_FUNC(pioFileOpen));
      namespaceSwitchSet(NSSWITCH_FILE_CLOSE, NSSW_FUNC(pioFileClose));
      cdiPioServer(postCommSetupActions);
      namespaceNew();
      commDestroy ();
      if (xtInitByCDI)
        xt_finalize();
      return MPI_COMM_NULL;
    }
  else
    cdiPioClientSetup(&pioNamespace_, pioNamespace);

  xdebug ( "nProcsGlob=%d, RETURN", sizeGlob );
  return commInqCommModel ();
}

/*****************************************************************************/

void  pioEndDef ( void )
{
  char   * buffer;
  int bufferSize;
  int rankGlob = commInqRankGlob ();

  xdebug("%s", "START");

  varsMapNDeco ( commInqNNodes (), commInqNodeSizes ());

  if ( rankGlob < commInqNProcsColl ())
    {
      MPI_Comm comm = commInqCommsIO ( rankGlob );
      reshPackBufferCreate(&buffer, &bufferSize, &comm);

      xmpi ( MPI_Send ( buffer, bufferSize, MPI_PACKED, commInqNProcsModel (),
                        RESOURCES, commInqCommsIO ( rankGlob )));

      xdebug("%s", "SENT MESSAGE WITH TAG \"RESOURCES\"");

      reshPackBufferDestroy ( &buffer );
    }

  modelWinCreate ();
  namespaceDefResStatus ( STAGE_TIMELOOP );
  xdebug("%s", "RETURN");
}

/************************************************************************/

void  pioEndTimestepping ( void )
{
  xdebug("%s", "START");
  namespaceDefResStatus ( STAGE_CLEANUP );
  xdebug("%s", "RETURN");
}


/****************************************************/


/**
  @brief is invoked by the calculator PEs, to inform
  the I/O PEs that no more data will be written.

  @param

  @return
*/

void pioFinalize ( void )
{
  xdebug("%s", "START");

  /* pioNamespace_ is unchanged on I/O servers */
  if (pioNamespace_ == -1)
    return;
  namespaceDelete(pioNamespace_);
  for (int collID = 0; collID < commInqNProcsColl (); collID++ )
    {
      xmpi(MPI_Send(NULL, 0, MPI_INT, commInqNProcsModel(),
                    FINALIZE, commInqCommsIO ( collID )));
      xdebug("%s", "SENT MESSAGE WITH TAG \"FINALIZE\"");
    }
  modelWinCleanup ();
  commDestroy ();
  if (xtInitByCDI)
    xt_finalize();
  xdebug("%s", "RETURN");
}

 /************************************************************************/

void pioWriteTimestep()
{
  int collID, iAssert = MPI_MODE_NOPUT;
  /* int tokenEnd = END; */
  int rankGlob = commInqRankGlob ();
  int nProcsColl = commInqNProcsColl ();
  int nProcsModel = commInqNProcsModel ();

  xdebug("%s", "START");

  xassert(txWin != NULL);

  if ( rankGlob < nProcsColl )
    {
      xmpi(MPI_Send(NULL, 0, MPI_INT, nProcsModel,
                    WRITETS, commInqCommsIO(rankGlob)));
      xdebug("%s", "SENT MESSAGE WITH TAG \"WRITETS\"");
    }

  for ( collID = 0; collID < nProcsColl; collID++ )
    {
      collWait(collID);
      struct winHeaderEntry header
        = { .id = HEADERSIZEMARKER,
            .specific.headerSize
            = { .numDataEntries = txWin[collID].dictDataUsed,
                .numRPCEntries = txWin[collID].dictRPCUsed } };
      struct winHeaderEntry *winDict
        = (struct winHeaderEntry *)txWin[collID].buffer;
      winDict[0] = header;

      xmpi(MPI_Win_post(txWin[collID].ioGroup, iAssert, txWin[collID].win));
      txWin[collID].postSet = 1;
    }

  xdebug("%s", "RETURN. messages sent, windows posted");
}

void
streamWriteVarPart(int streamID, int varID, const void *data,
                   int nmiss, Xt_idxlist partDesc)
{
  if ( CDI_Debug ) Message("streamID = %d varID = %d", streamID, varID);

  int chunk = xt_idxlist_get_num_indices(partDesc);
  xassert(chunk == 0 || data);

  void (*myStreamWriteVarPart)(int streamID, int varID, const void *data,
                               int nmiss, Xt_idxlist partDesc)
    = (void (*)(int, int, const void *, int, Xt_idxlist))
    namespaceSwitchGet(NSSWITCH_STREAM_WRITE_VAR_PART_).func;

  if (!myStreamWriteVarPart)
    xabort("local part writing is unsupported!");

  myStreamWriteVarPart(streamID, varID, data, nmiss, partDesc);
}

void
streamWriteScatteredVarPart(int streamID, int varID, const void *data,
                            int numBlocks, const int blocklengths[],
                            const int displacements[],
                            int nmiss, Xt_idxlist partDesc)
{
  if ( CDI_Debug ) Message("streamID = %d varID = %d", streamID, varID);

  int chunk = xt_idxlist_get_num_indices(partDesc);
  xassert(chunk == 0 || data);

  void (*myStreamWriteScatteredVarPart)(int streamID, int varID,
                                        const void *data,
                                        int numBlocks, const int blocklengths[],
                                        const int displacements[],
                                        int nmiss, Xt_idxlist partDesc)
    = (void (*)(int, int, const void *, int, const int[], const int[], int,
                Xt_idxlist))
    namespaceSwitchGet(NSSWITCH_STREAM_WRITE_SCATTERED_VAR_PART_).func;

  if (!myStreamWriteScatteredVarPart)
    xabort("local part writing is unsupported!");

  myStreamWriteScatteredVarPart(streamID, varID, data,
                                numBlocks, blocklengths, displacements,
                                nmiss, partDesc);
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
