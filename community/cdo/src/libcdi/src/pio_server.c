/** @file ioServer.c
*/
#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include "pio_server.h"


#include <limits.h>
#include <stdlib.h>
#include <stdio.h>

#ifdef HAVE_PARALLEL_NC4
#include <core/ppm_combinatorics.h>
#include <core/ppm_rectilinear.h>
#include <ppm/ppm_uniform_partition.h>
#endif
#include <yaxt.h>

#include "cdi.h"
#include "cdipio.h"
#include "dmemory.h"
#include "namespace.h"
#include "taxis.h"
#include "pio.h"
#include "pio_comm.h"
#include "pio_interface.h"
#include "pio_rpc.h"
#include "pio_util.h"
#include "cdi_int.h"
#ifndef HAVE_NETCDF_PAR_H
#define MPI_INCLUDED
#endif
#include "pio_cdf_int.h"
#include "resource_handle.h"
#include "resource_unpack.h"
#include "stream_cdf.h"
#include "vlist_var.h"


static struct
{
  size_t size;
  unsigned char *buffer;
  int dictSize;
} *rxWin = NULL;

static MPI_Win getWin = MPI_WIN_NULL;
static MPI_Group groupModel = MPI_GROUP_NULL;

#ifdef HAVE_PARALLEL_NC4
/* prime factorization of number of pio collectors */
static uint32_t *pioPrimes;
static int numPioPrimes;
#endif

/************************************************************************/

static
void serverWinCleanup ()
{
  if (getWin != MPI_WIN_NULL)
    xmpi(MPI_Win_free(&getWin));
  if (rxWin)
    {
      free(rxWin[0].buffer);
      free(rxWin);
    }

  xdebug("%s", "cleaned up mpi_win");
}

 /************************************************************************/

static size_t
collDefBufferSizes()
{
  int *streamIndexList, vlistID, nvars, varID, iorank;
  int modelID;
  size_t sumGetBufferSizes = 0;
  int rankGlob = commInqRankGlob ();
  int nProcsModel = commInqNProcsModel ();
  int root = commInqRootGlob ();

  xassert(rxWin != NULL);

  unsigned nstreams = reshCountType ( &streamOps );
  streamIndexList = xmalloc((size_t)nstreams * sizeof (streamIndexList[0]));
  reshGetResHListOfType ( nstreams, streamIndexList, &streamOps );
  for (unsigned streamNo = 0; streamNo < nstreams; streamNo++)
    {
      // space required for data
      vlistID = streamInqVlist ( streamIndexList[streamNo] );
      nvars = vlistNvars ( vlistID );
      for ( varID = 0; varID < nvars; varID++ )
        {
          iorank = vlistInqVarIOrank ( vlistID, varID );
          xassert ( iorank != CDI_UNDEFID );
          if ( iorank == rankGlob )
            {
              for ( modelID = 0; modelID < nProcsModel; modelID++ )
                {
                  int decoChunk;
                  {
                    int varSize = vlistInqVarSize(vlistID, varID);
                    int nProcsModel = commInqNProcsModel();
                    decoChunk =
                      (int)ceilf(cdiPIOpartInflate_
                                 * (float)(varSize + nProcsModel - 1)
                                 / (float)nProcsModel);
                  }
                  xassert ( decoChunk > 0 );
                  rxWin[modelID].size += (size_t)decoChunk * sizeof (double)
                    /* re-align chunks to multiple of double size */
                    + sizeof (double) - 1
                    /* one header for data record, one for
                     * corresponding part descriptor*/
                    + 2 * sizeof (struct winHeaderEntry)
                    /* FIXME: heuristic for size of packed Xt_idxlist */
                    + sizeof (Xt_int) * (size_t)decoChunk * 3;
                  rxWin[modelID].dictSize += 2;
                }
            }
        }
      // space required for the 3 function calls streamOpen, streamDefVlist, streamClose 
      // once per stream and timestep for all collprocs only on the modelproc root
      rxWin[root].size += numRPCFuncs * sizeof (struct winHeaderEntry)
        /* serialized filename */
        + MAXDATAFILENAME
        /* data part of streamDefTimestep */
        + (2 * CDI_MAX_NAME + sizeof (taxis_t));
      rxWin[root].dictSize += numRPCFuncs;
    }
  free ( streamIndexList );

  for ( modelID = 0; modelID < nProcsModel; modelID++ )
    {
      /* account for size header */
      rxWin[modelID].dictSize += 1;
      rxWin[modelID].size += sizeof (struct winHeaderEntry);
      rxWin[modelID].size = roundUpToMultiple(rxWin[modelID].size,
                                              PIO_WIN_ALIGN);
      sumGetBufferSizes += (size_t)rxWin[modelID].size;
    }
  xassert ( sumGetBufferSizes <= MAXWINBUFFERSIZE );
  return sumGetBufferSizes;
}

 /************************************************************************/

static void
serverWinCreate(void)
{
  int ranks[1], modelID;
  MPI_Comm commCalc = commInqCommCalc ();
  MPI_Group groupCalc;
  int nProcsModel = commInqNProcsModel ();
  MPI_Info no_locks_info;
  xmpi(MPI_Info_create(&no_locks_info));
  xmpi(MPI_Info_set(no_locks_info, "no_locks", "true"));

  xmpi(MPI_Win_create(MPI_BOTTOM, 0, 1, no_locks_info, commCalc, &getWin));

  /* target group */
  ranks[0] = nProcsModel;
  xmpi ( MPI_Comm_group ( commCalc, &groupCalc ));
  xmpi ( MPI_Group_excl ( groupCalc, 1, ranks, &groupModel ));

  rxWin = xcalloc((size_t)nProcsModel, sizeof (rxWin[0]));
  size_t totalBufferSize = collDefBufferSizes();
  rxWin[0].buffer = (unsigned char*) xmalloc(totalBufferSize);
  size_t ofs = 0;
  for ( modelID = 1; modelID < nProcsModel; modelID++ )
    {
      ofs += rxWin[modelID - 1].size;
      rxWin[modelID].buffer = rxWin[0].buffer + ofs;
    }

  xmpi(MPI_Info_free(&no_locks_info));

  xdebug("%s", "created mpi_win, allocated getBuffer");
}

/************************************************************************/

static void
readFuncCall(struct winHeaderEntry *header)
{
  int root = commInqRootGlob ();
  int funcID = header->id;
  union funcArgs *funcArgs = &(header->specific.funcArgs);

  xassert(funcID >= MINFUNCID && funcID <= MAXFUNCID);
  switch ( funcID )
    {
    case STREAMCLOSE:
      {
        int streamID
          = namespaceAdaptKey2(funcArgs->streamChange.streamID);
        streamClose(streamID);
        xdebug("READ FUNCTION CALL FROM WIN:  %s, streamID=%d,"
               " closed stream",
               funcMap[(-1 - funcID)], streamID);
      }
      break;
    case STREAMOPEN:
      {
        size_t filenamesz = (size_t)funcArgs->newFile.fnamelen;
        xassert ( filenamesz > 0 && filenamesz < MAXDATAFILENAME );
        const char *filename
          = (const char *)(rxWin[root].buffer + header->offset);
        xassert(filename[filenamesz] == '\0');
        int filetype = funcArgs->newFile.filetype;
        int streamID = streamOpenWrite(filename, filetype);
        xassert(streamID != CDI_ELIBNAVAIL);
        xdebug("READ FUNCTION CALL FROM WIN:  %s, filenamesz=%zu,"
               " filename=%s, filetype=%d, OPENED STREAM %d",
               funcMap[(-1 - funcID)], filenamesz, filename,
               filetype, streamID);
      }
      break;
    case STREAMDEFVLIST:
      {
        int streamID
          = namespaceAdaptKey2(funcArgs->streamChange.streamID);
        int vlistID = namespaceAdaptKey2(funcArgs->streamChange.vlistID);
        streamDefVlist(streamID, vlistID);
        xdebug("READ FUNCTION CALL FROM WIN:  %s, streamID=%d,"
               " vlistID=%d, called streamDefVlist ().",
               funcMap[(-1 - funcID)], streamID, vlistID);
      }
      break;
    case STREAMDEFTIMESTEP:
      {
        MPI_Comm commCalc = commInqCommCalc ();
        int streamID = funcArgs->streamNewTimestep.streamID;
        int originNamespace = namespaceResHDecode(streamID).nsp;
        streamID = namespaceAdaptKey2(streamID);
        int oldTaxisID
          = vlistInqTaxis(streamInqVlist(streamID));
        int position = header->offset;
        int changedTaxisID
          = taxisUnpack((char *)rxWin[root].buffer, (int)rxWin[root].size,
                        &position, originNamespace, &commCalc, 0);
        taxis_t *oldTaxisPtr = taxisPtr(oldTaxisID);
        taxis_t *changedTaxisPtr = taxisPtr(changedTaxisID);
        ptaxisCopy(oldTaxisPtr, changedTaxisPtr);
        taxisDestroy(changedTaxisID);
        streamDefTimestep(streamID, funcArgs->streamNewTimestep.tsID);
      }
      break;
    default:
      xabort ( "REMOTE FUNCTIONCALL NOT IMPLEMENTED!" );
    }
}

/************************************************************************/

static void
resizeVarGatherBuf(int vlistID, int varID, double **buf, int *bufSize)
{
  int size = vlistInqVarSize(vlistID, varID);
  if (size <= *bufSize) ; else
    *buf = xrealloc(*buf, (size_t)(*bufSize = size) * sizeof (buf[0][0]));
}

static void
gatherArray(int root, int nProcsModel, int headerIdx,
            int vlistID,
            double *gatherBuf, int *nmiss)
{
  struct winHeaderEntry *winDict
    = (struct winHeaderEntry *)rxWin[root].buffer;
  int streamID = winDict[headerIdx].id;
  int varID = winDict[headerIdx].specific.dataRecord.varID;
  int varShape[3] = { 0, 0, 0 };
  cdiPioQueryVarDims(varShape, vlistID, varID);
  Xt_int varShapeXt[3];
  static const Xt_int origin[3] = { 0, 0, 0 };
  for (unsigned i = 0; i < 3; ++i)
    varShapeXt[i] = varShape[i];
  int varSize = varShape[0] * varShape[1] * varShape[2];
  struct Xt_offset_ext *partExts
    = xmalloc((size_t)nProcsModel * sizeof (partExts[0]));
  Xt_idxlist *part = xmalloc((size_t)nProcsModel * sizeof (part[0]));
  MPI_Comm commCalc = commInqCommCalc();
  {
    int nmiss_ = 0;
    for (int modelID = 0; modelID < nProcsModel; modelID++)
      {
        struct dataRecord *dataHeader
          = &((struct winHeaderEntry *)
              rxWin[modelID].buffer)[headerIdx].specific.dataRecord;
        int position =
          ((struct winHeaderEntry *)rxWin[modelID].buffer)[headerIdx + 1].offset;
        xassert(namespaceAdaptKey2(((struct winHeaderEntry *)
                                    rxWin[modelID].buffer)[headerIdx].id)
                == streamID
                && dataHeader->varID == varID
                && ((struct winHeaderEntry *)
                    rxWin[modelID].buffer)[headerIdx + 1].id == PARTDESCMARKER
                && position > 0
                && ((size_t)position
                    >= sizeof (struct winHeaderEntry) * (size_t)rxWin[modelID].dictSize)
                && ((size_t)position < rxWin[modelID].size));
        part[modelID] = xt_idxlist_unpack(rxWin[modelID].buffer,
                                          (int)rxWin[modelID].size,
                                          &position, commCalc);
        unsigned partSize = (unsigned)xt_idxlist_get_num_indices(part[modelID]);
        size_t charOfs = (size_t)((rxWin[modelID].buffer
                                   + ((struct winHeaderEntry *)
                                      rxWin[modelID].buffer)[headerIdx].offset)
                                  - rxWin[0].buffer);
        xassert(charOfs % sizeof (double) == 0
                && charOfs / sizeof (double) + partSize <= INT_MAX);
        int elemOfs = (int)(charOfs / sizeof (double));
        partExts[modelID].start = elemOfs;
        partExts[modelID].size = (int)partSize;
        partExts[modelID].stride = 1;
        nmiss_ += dataHeader->nmiss;
      }
    *nmiss = nmiss_;
  }
  Xt_idxlist srcList = xt_idxlist_collection_new(part, nProcsModel);
  for (int modelID = 0; modelID < nProcsModel; modelID++)
    xt_idxlist_delete(part[modelID]);
  free(part);
  Xt_xmap gatherXmap;
  {
    Xt_idxlist dstList
      = xt_idxsection_new(0, 3, varShapeXt, varShapeXt, origin);
    struct Xt_com_list full = { .list = dstList, .rank = 0 };
    gatherXmap = xt_xmap_intersection_new(1, &full, 1, &full, srcList, dstList,
                                          MPI_COMM_SELF);
    xt_idxlist_delete(dstList);
  }
  xt_idxlist_delete(srcList);

  struct Xt_offset_ext gatherExt = { .start = 0, .size = varSize, .stride = 1 };
  Xt_redist gatherRedist
    = xt_redist_p2p_ext_new(gatherXmap, nProcsModel, partExts, 1, &gatherExt,
                            MPI_DOUBLE);
  xt_xmap_delete(gatherXmap);
  xt_redist_s_exchange1(gatherRedist, rxWin[0].buffer, gatherBuf);
  free(partExts);
  xt_redist_delete(gatherRedist);
}

struct xyzDims
{
  int sizes[3];
};

static inline int
xyzGridSize(struct xyzDims dims)
{
  return dims.sizes[0] * dims.sizes[1] * dims.sizes[2];
}

#ifdef HAVE_PARALLEL_NC4
static void
queryVarBounds(struct PPM_extent varShape[3], int vlistID, int varID)
{
  varShape[0].first = 0;
  varShape[1].first = 0;
  varShape[2].first = 0;
  int sizes[3];
  cdiPioQueryVarDims(sizes, vlistID, varID);
  for (unsigned i = 0; i < 3; ++i)
    varShape[i].size = sizes[i];
}

/* compute distribution of collectors such that number of collectors
 * <= number of variable grid cells in each dimension */
static struct xyzDims
varDimsCollGridMatch(const struct PPM_extent varDims[3])
{
  xassert(PPM_extents_size(3, varDims) >= commInqSizeColl());
  struct xyzDims collGrid = { { 1, 1, 1 } };
  /* because of storage order, dividing dimension 3 first is preferred */
  for (int i = 0; i < numPioPrimes; ++i)
    {
      for (int dim = 2; dim >=0; --dim)
        if (collGrid.sizes[dim] * pioPrimes[i] <= varDims[dim].size)
          {
            collGrid.sizes[dim] *= pioPrimes[i];
            goto nextPrime;
          }
      /* no position found, retrack */
      xabort("Not yet implemented back-tracking needed.");
      nextPrime:
      ;
    }
  return collGrid;
}

static void
myVarPart(struct PPM_extent varShape[3], struct xyzDims collGrid,
          struct PPM_extent myPart[3])
{
  int32_t myCollGridCoord[3];
  {
    struct PPM_extent collGridShape[3];
    for (int i = 0; i < 3; ++i)
      {
        collGridShape[i].first = 0;
        collGridShape[i].size = collGrid.sizes[i];
      }
    PPM_lidx2rlcoord_e(3, collGridShape, commInqRankColl(), myCollGridCoord);
    xdebug("my coord: (%d, %d, %d)", myCollGridCoord[0], myCollGridCoord[1],
           myCollGridCoord[2]);
  }
  PPM_uniform_partition_nd(3, varShape, collGrid.sizes,
                           myCollGridCoord, myPart);
}
#elif defined (HAVE_LIBNETCDF)
/* needed for writing when some files are only written to by a single process */
/* cdiOpenFileMap(fileID) gives the writer process */
int cdiPioSerialOpenFileMap(int streamID)
{
  return stream_to_pointer(streamID)->ownerRank;
}
/* for load-balancing purposes, count number of files per process */
/* cdiOpenFileCounts[rank] gives number of open files rank has to himself */
static int *cdiSerialOpenFileCount = NULL;

static int
cdiPioNextOpenRank()
{
  xassert(cdiSerialOpenFileCount != NULL);
  int commCollSize = commInqSizeColl();
  int minRank = 0, minOpenCount = cdiSerialOpenFileCount[0];
  for (int i = 1; i < commCollSize; ++i)
    if (cdiSerialOpenFileCount[i] < minOpenCount)
      {
        minOpenCount = cdiSerialOpenFileCount[i];
        minRank = i;
      }
  return minRank;
}

static void
cdiPioOpenFileOnRank(int rank)
{
  xassert(cdiSerialOpenFileCount != NULL
          && (unsigned)rank < (unsigned)commInqSizeColl());
  ++(cdiSerialOpenFileCount[rank]);
}

static void
cdiPioCloseFileOnRank(int rank)
{
  xassert(cdiSerialOpenFileCount != NULL
          && rank >= 0 && rank < commInqSizeColl());
  xassert(cdiSerialOpenFileCount[rank] > 0);
  --(cdiSerialOpenFileCount[rank]);
}

static void
cdiPioServerCdfDefVars(stream_t *streamptr)
{
  int rank, rankOpen;
  if (commInqIOMode() == PIO_NONE
      || ((rank = commInqRankColl())
          == (rankOpen = cdiPioSerialOpenFileMap(streamptr->self))))
    cdfDefVars(streamptr);
}

#endif

struct streamMapping {
  int streamID, filetype;
  int firstHeaderIdx, lastHeaderIdx;
  int numVars, *varMap;
};

struct streamMap
{
  struct streamMapping *entries;
  int numEntries;
};

static int
smCmpStreamID(const void *a_, const void *b_)
{
  const struct streamMapping *a = a_, *b = b_;
  int streamIDa = a->streamID, streamIDb = b->streamID;
  return (streamIDa > streamIDb) - (streamIDa < streamIDb);
}

static inline int
inventorizeStream(struct streamMapping *streamMap, int numStreamIDs,
                  int *sizeStreamMap_, int streamID, int headerIdx)
{
  int sizeStreamMap = *sizeStreamMap_;
  if (numStreamIDs < sizeStreamMap) ; else
    {
      streamMap = xrealloc(streamMap,
                           (size_t)(sizeStreamMap *= 2)
                           * sizeof (streamMap[0]));
      *sizeStreamMap_ = sizeStreamMap;
    }
  streamMap[numStreamIDs].streamID = streamID;
  streamMap[numStreamIDs].firstHeaderIdx = headerIdx;
  streamMap[numStreamIDs].lastHeaderIdx = headerIdx;
  streamMap[numStreamIDs].numVars = -1;
  int filetype = streamInqFiletype(streamID);
  streamMap[numStreamIDs].filetype = filetype;
  if (filetype == FILETYPE_NC || filetype == FILETYPE_NC2
      || filetype == FILETYPE_NC4)
    {
      int vlistID = streamInqVlist(streamID);
      int nvars = vlistNvars(vlistID);
      streamMap[numStreamIDs].numVars = nvars;
      streamMap[numStreamIDs].varMap
        = xmalloc(sizeof (streamMap[numStreamIDs].varMap[0]) * (size_t)nvars);
      for (int i = 0; i < nvars; ++i)
        streamMap[numStreamIDs].varMap[i] = -1;
    }
  return numStreamIDs + 1;
}

static inline int
streamIsInList(struct streamMapping *streamMap, int numStreamIDs,
               int streamIDQuery)
{
  int p = 0;
  for (int i = 0; i < numStreamIDs; ++i)
    p |= streamMap[i].streamID == streamIDQuery;
  return p;
}

static struct streamMap
buildStreamMap(struct winHeaderEntry *winDict)
{
  int streamIDOld = CDI_UNDEFID;
  int oldStreamIdx = CDI_UNDEFID;
  int filetype = FILETYPE_UNDEF;
  int sizeStreamMap = 16;
  struct streamMapping *streamMap
    = xmalloc((size_t)sizeStreamMap * sizeof (streamMap[0]));
  int numDataEntries = winDict[0].specific.headerSize.numDataEntries;
  int numStreamIDs = 0;
  /* find streams written on this process */
  for (int headerIdx = 1; headerIdx < numDataEntries; headerIdx += 2)
    {
      int streamID
        = winDict[headerIdx].id
        = namespaceAdaptKey2(winDict[headerIdx].id);
      xassert(streamID > 0);
      if (streamID != streamIDOld)
        {
          for (int i = numStreamIDs - 1; i >= 0; --i)
            if ((streamIDOld = streamMap[i].streamID) == streamID)
              {
                oldStreamIdx = i;
                goto streamIDInventorized;
              }
          oldStreamIdx = numStreamIDs;
          streamIDOld = streamID;
          numStreamIDs = inventorizeStream(streamMap, numStreamIDs,
                                           &sizeStreamMap, streamID, headerIdx);
        }
      streamIDInventorized:
      filetype = streamMap[oldStreamIdx].filetype;
      streamMap[oldStreamIdx].lastHeaderIdx = headerIdx;
      if (filetype == FILETYPE_NC || filetype == FILETYPE_NC2
          || filetype == FILETYPE_NC4)
        {
          int varID = winDict[headerIdx].specific.dataRecord.varID;
          streamMap[oldStreamIdx].varMap[varID] = headerIdx;
        }
    }
  /* join with list of streams written to in total */
  {
    int *streamIDs, *streamIsWritten;
    unsigned numTotalStreamIDs = reshCountType(&streamOps);
    streamIDs = xmalloc(2 * sizeof (streamIDs[0]) * (size_t)numTotalStreamIDs);
    cdiStreamGetIndexList(numTotalStreamIDs, streamIDs);
    streamIsWritten = streamIDs + numTotalStreamIDs;
    for (unsigned i = 0; i < numTotalStreamIDs; ++i)
      streamIsWritten[i] = streamIsInList(streamMap, numStreamIDs,
                                          streamIDs[i]);
    /* Find what streams are written to at all on any process */
    xmpi(MPI_Allreduce(MPI_IN_PLACE, streamIsWritten, (int)numTotalStreamIDs,
                       MPI_INT, MPI_BOR, commInqCommColl()));
    /* append streams written to on other tasks to mapping */
    for (unsigned i = 0; i < numTotalStreamIDs; ++i)
      if (streamIsWritten[i] && !streamIsInList(streamMap, numStreamIDs,
                                                streamIDs[i]))
        numStreamIDs = inventorizeStream(streamMap, numStreamIDs,
                                         &sizeStreamMap, streamIDs[i], -1);

    free(streamIDs);
  }
  /* sort written streams by streamID */
  streamMap = xrealloc(streamMap, sizeof (streamMap[0]) * (size_t)numStreamIDs);
  qsort(streamMap, (size_t)numStreamIDs, sizeof (streamMap[0]), smCmpStreamID);
  return (struct streamMap){ .entries = streamMap, .numEntries = numStreamIDs };
}

static void
writeGribStream(struct winHeaderEntry *winDict, struct streamMapping *mapping,
                double **data_, int *currentDataBufSize, int root,
                int nProcsModel)
{
  int streamID = mapping->streamID;
  int headerIdx, lastHeaderIdx = mapping->lastHeaderIdx;
  int vlistID = streamInqVlist(streamID);
  if (lastHeaderIdx < 0)
    {
      /* write zero bytes to trigger synchronization code in fileWrite */
      cdiPioFileWrite(streamInqFileID(streamID), NULL, 0,
                      streamInqCurTimestepID(streamID));
    }
  else
    for (headerIdx = mapping->firstHeaderIdx;
         headerIdx <= lastHeaderIdx;
         headerIdx += 2)
      if (streamID == winDict[headerIdx].id)
        {
          int varID = winDict[headerIdx].specific.dataRecord.varID;
          int size = vlistInqVarSize(vlistID, varID);
          int nmiss;
          resizeVarGatherBuf(vlistID, varID, data_, currentDataBufSize);
          double *data = *data_;
          gatherArray(root, nProcsModel, headerIdx,
                      vlistID, data, &nmiss);
          streamWriteVar(streamID, varID, data, nmiss);
          if ( ddebug > 2 )
            {
              char text[1024];
              sprintf(text, "streamID=%d, var[%d], size=%d",
                      streamID, varID, size);
              xprintArray(text, data, size, DATATYPE_FLT);
            }
        }
}

#ifdef HAVE_NETCDF4
static void
buildWrittenVars(struct streamMapping *mapping, int **varIsWritten_,
                 int myCollRank, MPI_Comm collComm)
{
  int nvars = mapping->numVars;
  int *varMap = mapping->varMap;
  int *varIsWritten = *varIsWritten_
    = xrealloc(*varIsWritten_, sizeof (*varIsWritten) * (size_t)nvars);
  for (int varID = 0; varID < nvars; ++varID)
    varIsWritten[varID] = ((varMap[varID] != -1)
                           ?myCollRank+1 : 0);
  xmpi(MPI_Allreduce(MPI_IN_PLACE, varIsWritten, nvars,
                     MPI_INT, MPI_BOR, collComm));
}
#endif

static void readGetBuffers()
{
  int nProcsModel = commInqNProcsModel ();
  int root        = commInqRootGlob ();
#ifdef HAVE_NETCDF4
  int myCollRank = commInqRankColl();
  MPI_Comm collComm = commInqCommColl();
#endif
  xdebug("%s", "START");

  struct winHeaderEntry *winDict
    = (struct winHeaderEntry *)rxWin[root].buffer;
  xassert(winDict[0].id == HEADERSIZEMARKER);
  {
    int dictSize = rxWin[root].dictSize,
      firstNonRPCEntry = dictSize - winDict[0].specific.headerSize.numRPCEntries - 1,
      headerIdx,
      numFuncCalls = 0;
    for (headerIdx = dictSize - 1;
         headerIdx > firstNonRPCEntry;
         --headerIdx)
      {
        xassert(winDict[headerIdx].id >= MINFUNCID
                && winDict[headerIdx].id <= MAXFUNCID);
        ++numFuncCalls;
        readFuncCall(winDict + headerIdx);
      }
    xassert(numFuncCalls == winDict[0].specific.headerSize.numRPCEntries);
  }
  /* build list of streams, data was transferred for */
  {
    struct streamMap map = buildStreamMap(winDict);
    double *data = NULL;
#ifdef HAVE_NETCDF4
    int *varIsWritten = NULL;
#endif
#if defined (HAVE_PARALLEL_NC4)
    double *writeBuf = NULL;
#endif
    int currentDataBufSize = 0;
    for (int streamIdx = 0; streamIdx < map.numEntries; ++streamIdx)
      {
        int streamID = map.entries[streamIdx].streamID;
        int vlistID = streamInqVlist(streamID);
        int filetype = map.entries[streamIdx].filetype;

        switch (filetype)
          {
          case FILETYPE_GRB:
          case FILETYPE_GRB2:
            writeGribStream(winDict, map.entries + streamIdx,
                            &data, &currentDataBufSize,
                            root, nProcsModel);
            break;
#ifdef HAVE_NETCDF4
          case FILETYPE_NC:
          case FILETYPE_NC2:
          case FILETYPE_NC4:
#ifdef HAVE_PARALLEL_NC4
            /* HAVE_PARALLE_NC4 implies having ScalES-PPM and yaxt */
            {
              int nvars = map.entries[streamIdx].numVars;
              int *varMap = map.entries[streamIdx].varMap;
              buildWrittenVars(map.entries + streamIdx, &varIsWritten,
                               myCollRank, collComm);
              for (int varID = 0; varID < nvars; ++varID)
                if (varIsWritten[varID])
                  {
                    struct PPM_extent varShape[3];
                    queryVarBounds(varShape, vlistID, varID);
                    struct xyzDims collGrid = varDimsCollGridMatch(varShape);
                    xdebug("writing varID %d with dimensions: "
                           "x=%d, y=%d, z=%d,\n"
                           "found distribution with dimensions:"
                           " x=%d, y=%d, z=%d.", varID,
                           varShape[0].size, varShape[1].size, varShape[2].size,
                           collGrid.sizes[0], collGrid.sizes[1],
                           collGrid.sizes[2]);
                    struct PPM_extent varChunk[3];
                    myVarPart(varShape, collGrid, varChunk);
                    int myChunk[3][2];
                    for (int i = 0; i < 3; ++i)
                      {
                        myChunk[i][0] = PPM_extent_start(varChunk[i]);
                        myChunk[i][1] = PPM_extent_end(varChunk[i]);
                      }
                    xdebug("Writing chunk { { %d, %d }, { %d, %d },"
                           " { %d, %d } }", myChunk[0][0], myChunk[0][1],
                           myChunk[1][0], myChunk[1][1], myChunk[2][0],
                           myChunk[2][1]);
                    Xt_int varSize[3];
                    for (int i = 0; i < 3; ++i)
                      varSize[2 - i] = varShape[i].size;
                    Xt_idxlist preRedistChunk, preWriteChunk;
                    /* prepare yaxt descriptor for current data
                       distribution after collect */
                    int nmiss;
                    if (varMap[varID] == -1)
                      {
                        preRedistChunk = xt_idxempty_new();
                        xdebug("%s", "I got none\n");
                      }
                    else
                      {
                        Xt_int preRedistStart[3] = { 0, 0, 0 };
                        preRedistChunk
                          = xt_idxsection_new(0, 3, varSize, varSize,
                                              preRedistStart);
                        resizeVarGatherBuf(vlistID, varID, &data,
                                           &currentDataBufSize);
                        int headerIdx = varMap[varID];
                        gatherArray(root, nProcsModel, headerIdx,
                                    vlistID, data, &nmiss);
                        xdebug("%s", "I got all\n");
                      }
                    MPI_Bcast(&nmiss, 1, MPI_INT, varIsWritten[varID] - 1,
                              collComm);
                    /* prepare yaxt descriptor for write chunk */
                    {
                      Xt_int preWriteChunkStart[3], preWriteChunkSize[3];
                      for (int i = 0; i < 3; ++i)
                        {
                          preWriteChunkStart[2 - i] = varChunk[i].first;
                          preWriteChunkSize[2 - i] = varChunk[i].size;
                        }
                      preWriteChunk = xt_idxsection_new(0, 3, varSize,
                                                        preWriteChunkSize,
                                                        preWriteChunkStart);
                    }
                    /* prepare redistribution */
                    {
                      Xt_xmap xmap = xt_xmap_all2all_new(preRedistChunk,
                                                         preWriteChunk,
                                                         collComm);
                      Xt_redist redist = xt_redist_p2p_new(xmap, MPI_DOUBLE);
                      xt_idxlist_delete(preRedistChunk);
                      xt_idxlist_delete(preWriteChunk);
                      xt_xmap_delete(xmap);
                      writeBuf = (double*) xrealloc(writeBuf,
                                                    sizeof (double)
                                                    * PPM_extents_size(3, varChunk));
                      xt_redist_s_exchange1(redist, data, writeBuf);
                      xt_redist_delete(redist);
                    }
                    /* write chunk */
                    streamWriteVarChunk(streamID, varID,
                                        (const int (*)[2])myChunk, writeBuf,
                                        nmiss);
                  }
            }
#else
            /* determine process which has stream open (writer) and
             * which has data for which variable (var owner)
             * three cases need to be distinguished */
            {
              int nvars = map.entries[streamIdx].numVars;
              int *varMap = map.entries[streamIdx].varMap;
              buildWrittenVars(map.entries + streamIdx, &varIsWritten,
                               myCollRank, collComm);
              int writerRank;
              if ((writerRank = cdiPioSerialOpenFileMap(streamID))
                  == myCollRank)
                {
                  for (int varID = 0; varID < nvars; ++varID)
                    if (varIsWritten[varID])
                      {
                        int nmiss;
                        int size = vlistInqVarSize(vlistID, varID);
                        resizeVarGatherBuf(vlistID, varID, &data,
                                           &currentDataBufSize);
                        int headerIdx = varMap[varID];
                        if (varIsWritten[varID] == myCollRank + 1)
                          {
                            /* this process has the full array and will
                             * write it */
                            xdebug("gathering varID=%d for direct writing",
                                   varID);
                            gatherArray(root, nProcsModel, headerIdx,
                                        vlistID, data, &nmiss);
                          }
                        else
                          {
                            /* another process has the array and will
                             * send it over */
                            MPI_Status stat;
                            xdebug("receiving varID=%d for writing from"
                                   " process %d",
                                   varID, varIsWritten[varID] - 1);
                            xmpiStat(MPI_Recv(&nmiss, 1, MPI_INT,
                                              varIsWritten[varID] - 1,
                                              COLLBUFNMISS,
                                              collComm, &stat), &stat);
                            xmpiStat(MPI_Recv(data, size, MPI_DOUBLE,
                                              varIsWritten[varID] - 1,
                                              COLLBUFTX,
                                              collComm, &stat), &stat);
                          }
                        streamWriteVar(streamID, varID, data, nmiss);
                      }
                }
              else
                for (int varID = 0; varID < nvars; ++varID)
                  if (varIsWritten[varID] == myCollRank + 1)
                    {
                      /* this process has the full array and another
                       * will write it */
                      int nmiss;
                      int size = vlistInqVarSize(vlistID, varID);
                      resizeVarGatherBuf(vlistID, varID, &data,
                                         &currentDataBufSize);
                      int headerIdx = varMap[varID];
                      gatherArray(root, nProcsModel, headerIdx,
                                  vlistID, data, &nmiss);
                      MPI_Request req;
                      MPI_Status stat;
                      xdebug("sending varID=%d for writing to"
                             " process %d",
                             varID, writerRank);
                      xmpi(MPI_Isend(&nmiss, 1, MPI_INT,
                                     writerRank, COLLBUFNMISS,
                                     collComm, &req));
                      xmpi(MPI_Send(data, size, MPI_DOUBLE,
                                    writerRank, COLLBUFTX,
                                    collComm));
                      xmpiStat(MPI_Wait(&req, &stat), &stat);
                    }
            }
#endif
            break;
#endif
          default:
            xabort("unhandled filetype in parallel I/O.");
          }
      }
#ifdef HAVE_NETCDF4
    free(varIsWritten);
#ifdef HAVE_PARALLEL_NC4
    free(writeBuf);
#endif
#endif
    free(map.entries);
    free(data);
  }
  xdebug("%s", "RETURN");
} 

/************************************************************************/


static
void clearModelWinBuffer(int modelID)
{
  int nProcsModel = commInqNProcsModel ();

  xassert((unsigned)modelID < (unsigned)nProcsModel &&
          rxWin != NULL && rxWin[modelID].buffer != NULL &&
          rxWin[modelID].size > 0 &&
          rxWin[modelID].size <= MAXWINBUFFERSIZE );
  memset(rxWin[modelID].buffer, 0, rxWin[modelID].size);
}


/************************************************************************/


static
void getTimeStepData()
{
  int modelID;
  char text[1024];
  int nProcsModel = commInqNProcsModel ();
  void *getWinBaseAddr;
  int attrFound;

  xdebug("%s", "START");

  for ( modelID = 0; modelID < nProcsModel; modelID++ )
    clearModelWinBuffer(modelID);
  // todo put in correct lbs and ubs
  xmpi(MPI_Win_start(groupModel, 0, getWin));
  xmpi(MPI_Win_get_attr(getWin, MPI_WIN_BASE, &getWinBaseAddr, &attrFound));
  xassert(attrFound);
  for ( modelID = 0; modelID < nProcsModel; modelID++ )
    {
      xdebug("modelID=%d, nProcsModel=%d, rxWin[%d].size=%zu,"
             " getWin=%p, sizeof(int)=%u",
             modelID, nProcsModel, modelID, rxWin[modelID].size,
             getWinBaseAddr, (unsigned)sizeof(int));
      /* FIXME: this needs to use MPI_PACK for portability */
      xmpi(MPI_Get(rxWin[modelID].buffer, (int)rxWin[modelID].size,
                   MPI_UNSIGNED_CHAR, modelID, 0,
                   (int)rxWin[modelID].size, MPI_UNSIGNED_CHAR, getWin));
    }
  xmpi ( MPI_Win_complete ( getWin ));

  if ( ddebug > 2 )
    for ( modelID = 0; modelID < nProcsModel; modelID++ )
      {
        sprintf(text, "rxWin[%d].size=%zu from PE%d rxWin[%d].buffer",
                modelID, rxWin[modelID].size, modelID, modelID);
        xprintArray(text, rxWin[modelID].buffer,
                    (int)(rxWin[modelID].size / sizeof (double)),
                    DATATYPE_FLT);
      }
  readGetBuffers();

  xdebug("%s", "RETURN");
}

/************************************************************************/

#if defined (HAVE_LIBNETCDF) && ! defined (HAVE_PARALLEL_NC4)
static int
cdiPioStreamCDFOpenWrap(const char *filename, const char *filemode,
                        int filetype, stream_t *streamptr,
                        int recordBufIsToBeCreated)
{
  switch (filetype)
    {
    case FILETYPE_NC4:
    case FILETYPE_NC4C:
      {
        /* Only needs initialization to shut up gcc */
        int rank = -1, fileID;
        int ioMode = commInqIOMode();
        if (ioMode == PIO_NONE
            || commInqRankColl() == (rank = cdiPioNextOpenRank()))
          fileID = cdiStreamOpenDefaultDelegate(filename, filemode, filetype,
                                                streamptr,
                                                recordBufIsToBeCreated);
        else
          streamptr->filetype = filetype;
        if (ioMode != PIO_NONE)
          xmpi(MPI_Bcast(&fileID, 1, MPI_INT, rank, commInqCommColl()));
        streamptr->ownerRank = rank;
        cdiPioOpenFileOnRank(rank);
        return fileID;
      }
    default:
      return cdiStreamOpenDefaultDelegate(filename, filemode, filetype,
                                          streamptr, recordBufIsToBeCreated);
    }
}

static void
cdiPioStreamCDFCloseWrap(stream_t *streamptr, int recordBufIsToBeDeleted)
{
  int fileID   = streamptr->fileID;
  int filetype = streamptr->filetype;
  if ( fileID == CDI_UNDEFID )
    Warning("File %s not open!", streamptr->filename);
  else
    switch (filetype)
      {
      case FILETYPE_NC:
      case FILETYPE_NC2:
      case FILETYPE_NC4:
      case FILETYPE_NC4C:
        {
          int rank, rankOpen = cdiPioSerialOpenFileMap(streamptr->self);
          if (commInqIOMode() == PIO_NONE
              || ((rank = commInqRankColl()) == rankOpen))
            cdiStreamCloseDefaultDelegate(streamptr, recordBufIsToBeDeleted);
          cdiPioCloseFileOnRank(rankOpen);
          break;
        }
      default:
        cdiStreamCloseDefaultDelegate(streamptr, recordBufIsToBeDeleted);
      }
}

static void
cdiPioCdfDefTimestep(stream_t *streamptr, int tsID)
{
  int rank, rankOpen, streamID = streamptr->self;
  if (commInqIOMode() == PIO_NONE
      || ((rank = commInqRankColl())
          == (rankOpen = cdiPioSerialOpenFileMap(streamID))))
    cdfDefTimestep(streamptr, tsID);
}

#endif

/**
  @brief is encapsulated in CDI library and run on I/O PEs.

  @param

  @return
*/

void cdiPioServer(void (*postCommSetupActions)(void))
{
  int nProcsModel=commInqNProcsModel();
  static int nfinished = 0;
  MPI_Comm commCalc;
  MPI_Status status;

  xdebug("%s", "START");

  cdiPioFileWritingInit(postCommSetupActions);
  if (commInqRankNode() == commInqSpecialRankNode())
    return;
  commCalc = commInqCommCalc ();
#ifdef HAVE_PARALLEL_NC4
  cdiPioEnableNetCDFParAccess();
  numPioPrimes = PPM_prime_factorization_32((uint32_t)commInqSizeColl(),
                                            &pioPrimes);
#elif defined (HAVE_LIBNETCDF)
  cdiSerialOpenFileCount = xcalloc(sizeof (cdiSerialOpenFileCount[0]),
                                   (size_t)commInqSizeColl());
  namespaceSwitchSet(NSSWITCH_STREAM_OPEN_BACKEND,
                     NSSW_FUNC(cdiPioStreamCDFOpenWrap));
  namespaceSwitchSet(NSSWITCH_STREAM_CLOSE_BACKEND,
                     NSSW_FUNC(cdiPioStreamCDFCloseWrap));
  namespaceSwitchSet(NSSWITCH_CDF_DEF_TIMESTEP,
                     NSSW_FUNC(cdiPioCdfDefTimestep));
  namespaceSwitchSet(NSSWITCH_CDF_STREAM_SETUP,
                     NSSW_FUNC(cdiPioServerCdfDefVars));
#endif
  namespaceSwitchSet(NSSWITCH_FILE_WRITE,
                     NSSW_FUNC(cdiPioFileWrite));

  for ( ;; )
    {
      xmpi ( MPI_Probe ( MPI_ANY_SOURCE, MPI_ANY_TAG, commCalc, &status ));
      
      int source = status.MPI_SOURCE;
      int tag = status.MPI_TAG;
      
      switch ( tag )
        {
        case FINALIZE:
          xdebugMsg(tag, source, nfinished);
          xmpi(MPI_Recv(NULL, 0, MPI_INT, source, tag, commCalc, &status));
          xdebug("%s", "RECEIVED MESSAGE WITH TAG \"FINALIZE\"");
          nfinished++;
          xdebug("nfinished=%d, nProcsModel=%d", nfinished, nProcsModel);
          if ( nfinished == nProcsModel )
            {
              {
                unsigned nStreams = reshCountType(&streamOps);

                if ( nStreams > 0 )
                  {
                    int *resHs = xmalloc(nStreams * sizeof (resHs[0]));
                    cdiStreamGetIndexList(nStreams, resHs);
                    for (unsigned streamNo = 0; streamNo < nStreams; ++streamNo)
                      streamClose(resHs[streamNo]);
                    free(resHs);
                  }
              }
              cdiPioFileWritingFinalize();
              serverWinCleanup();
#ifdef HAVE_PARALLEL_NC4
              free(pioPrimes);
#endif
              /* listDestroy(); */
              xdebug("%s", "RETURN");
              return;
            }
	  
          break;
          
	case RESOURCES:
          {
            int size;
            xdebugMsg(tag, source, nfinished);
            xmpi(MPI_Get_count(&status, MPI_CHAR, &size));
            char *buffer = xmalloc((size_t)size);
            xmpi(MPI_Recv(buffer, size, MPI_PACKED, source,
                          tag, commCalc, &status));
            xdebug("%s", "RECEIVED MESSAGE WITH TAG \"RESOURCES\"");
            reshUnpackResources(buffer, size, &commCalc);
            xdebug("%s", "");
            free(buffer);
            int rankGlob = commInqRankGlob();
            if ( ddebug > 0 && rankGlob == nProcsModel)
              {
                static const char baseName[] = "reshListIOServer.",
                  suffix[] = ".txt";
                /* 9 digits for rank at most */
                char buf[sizeof (baseName) + 9 + sizeof (suffix) + 1];
                snprintf(buf, sizeof (buf), "%s%d%s", baseName, rankGlob,
                         suffix);
                FILE *fp = fopen(buf, "w");
                xassert(fp);
                reshListPrint(fp);
                fclose(fp);
              }
          }
          serverWinCreate ();
	  break;
	case WRITETS:
          {
            xdebugMsg(tag, source, nfinished);
            xmpi(MPI_Recv(NULL, 0, MPI_INT, source,
                          tag, commCalc, &status));
            xdebug("RECEIVED MESSAGE WITH TAG \"WRITETS\": source=%d",
                   source);
            getTimeStepData();
          }
	  break;

	default:
	  xabort ( "TAG NOT DEFINED!" );
	}
    }
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
