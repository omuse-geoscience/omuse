#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>

#include "pio_comm.h"
#include "cdi.h"
#include "dmemory.h"
#include "pio_util.h"

#include "cdipio.h"

typedef struct {
  int IOMode;
  int nProcsIO;
  int nProcsModel;
  int isProcIO;

  MPI_Comm commGlob;
  int sizeGlob;
  int rankGlob;
  int root;

  MPI_Comm commPio;
  int sizePio;
  int rankPio;

  MPI_Comm commNode;
  int sizeNode;
  int rankNode;
  char hostname [ MPI_MAX_PROCESSOR_NAME + 1 ];
  nodeInfo_t nodeInfo;
  int specialRankNode;

  MPI_Comm commColl;
  int sizeColl;
  int rankColl;

  MPI_Comm commCalc;
  int rankCalc, sizeCalc;

  MPI_Comm * commsIO;
  int nProcsColl;
  int * procsCollMap;
  int * nodeSizes;
  int * nodeMap;
} pioInfo_t;


static pioInfo_t * info = NULL;


static
void pioInfoInit ( pioInfo_t * p )
{
  p->IOMode               = CDI_UNDEFID;
  p->nProcsIO             = CDI_UNDEFID;
  p->nProcsModel          = CDI_UNDEFID;
  p->isProcIO             = CDI_UNDEFID;

  p->commGlob             = MPI_COMM_NULL;
  p->sizeGlob             = CDI_UNDEFID;
  p->rankGlob             = CDI_UNDEFID;
  p->root                 = CDI_UNDEFID;

  p->commPio              = MPI_COMM_NULL;
  p->sizePio              = CDI_UNDEFID;
  p->rankPio              = -1;

  p->commNode             = MPI_COMM_NULL;
  p->sizeNode             = CDI_UNDEFID;
  p->rankNode             = CDI_UNDEFID;
  p->hostname[0]          = 0;
  p->specialRankNode      = CDI_UNDEFID;

  p->commColl             = MPI_COMM_NULL;
  p->sizeColl             = CDI_UNDEFID;
  p->rankColl             = CDI_UNDEFID;

  p->commCalc             = MPI_COMM_NULL;
  p->rankCalc             = -1;
  p->sizeCalc             = -1;

  p->commsIO              = NULL;
  p->nodeInfo.hostID      = CDI_UNDEFID;
  p->nodeInfo.isProcColl  = CDI_UNDEFID;
  p->nodeInfo.nNodes      = CDI_UNDEFID;
  p->nProcsColl           = CDI_UNDEFID;
  p->procsCollMap         = NULL;
  p->nodeSizes            = NULL;
  p->nodeMap              = NULL;
}


void commInit ( void )
{
  xassert ( info == 0 );
  info = xcalloc(1, sizeof (pioInfo_t));
  pioInfoInit ( info );
}


void commDestroy ( void )
{
  int collID;

  xassert ( info != NULL );

  if ( info->nodeMap   != NULL ) free ( info->nodeMap );
  if ( info->nodeSizes != NULL ) free ( info->nodeSizes );

  if ( info->commsIO != NULL )
    {
      for ( collID = 0; collID < info->nProcsColl; collID++ )
	if ( info->commsIO[collID] != MPI_COMM_NULL )
	  xmpi(MPI_Comm_free(info->commsIO + collID));
      free ( info->commsIO );
      info->commsIO = NULL;
    }

  if ( info->commColl != MPI_COMM_NULL )
    {
      xmpi ( MPI_Comm_free ( &info->commColl ));
      info->commColl = MPI_COMM_NULL;
    }

  free(info->procsCollMap);

  if ( info->commNode != MPI_COMM_NULL )
    {
      xmpi ( MPI_Comm_free ( &info->commNode ));
      info->commNode = MPI_COMM_NULL;
    }

  if ( info->commPio != MPI_COMM_NULL )
    {
      xmpi ( MPI_Comm_free ( &info->commPio ));
      info->commPio = MPI_COMM_NULL;
    }

  free ( info );
  info = NULL;
}


void commDefCommGlob ( MPI_Comm c )
{
  xassert ( info != NULL &&
	   c != MPI_COMM_NULL );

  info->commGlob = c;
  xmpi ( MPI_Comm_size ( c, &info->sizeGlob ));
  xmpi ( MPI_Comm_rank ( c, &info->rankGlob ));
  info->root = 0;
}
 

MPI_Comm commInqCommGlob ( void )
{
  xassert ( info != NULL &&
           info->commGlob != MPI_COMM_NULL );
  return info->commGlob;
}


int commInqSizeGlob ( void )
{
  xassert ( info != NULL  &&
           info->sizeGlob != CDI_UNDEFID );
  return info->sizeGlob;
}


int commInqRankGlob ( void )
{
  xassert ( info != NULL &&
           info->rankGlob != CDI_UNDEFID );
  return info->rankGlob;
}


int commInqRootGlob ( void )
{
  xassert ( info != NULL &&
           info->root != CDI_UNDEFID );
  return info->root;
}


void commDefNProcsIO ( int n )
{
  xassert ( info != NULL &&
	   n >= 0 &&
	   n < MAXNPROCSIO &&
           info->commGlob != MPI_COMM_NULL );
  info->nProcsIO = n;
  info->nProcsModel = info->sizeGlob - info->nProcsIO;
}


int commInqNProcsIO ( void )
{
  xassert ( info != NULL &&
           info->nProcsIO != CDI_UNDEFID );
  return info->nProcsIO;
}


int      commInqIsProcIO     ( void )
{
  xassert ( info != NULL &&
           info->isProcIO != CDI_UNDEFID );
  return info->isProcIO;
}


int commInqNProcsModel ( void )
{
  xassert ( info != NULL &&
           info->nProcsModel != CDI_UNDEFID );
  return info->nProcsModel;
}


void     commDefIOMode  ( int IOMode )
{
  xassert(info != NULL && IOMode >= PIO_MINIOMODE && IOMode <= PIO_MAXIOMODE );
  info->IOMode = IOMode;
}


int      commInqIOMode  ( void )
{
  if (info != NULL)
    {
      xassert ( info->IOMode != CDI_UNDEFID );
      return info->IOMode;
    }
  else
    return PIO_NONE;
}


void     commDefCommPio  ( void )
{
  int nProcsCalc;
  MPI_Group grpAll, grpCalc, grpPio;
  int calcRange[3], pioRange[3];

  xassert ( info != NULL &&
           info->commGlob != MPI_COMM_NULL &&
           info->IOMode  != CDI_UNDEFID &&
           info->nProcsIO != CDI_UNDEFID );

  nProcsCalc = info->sizeGlob - info->nProcsIO;

  info->isProcIO = info->rankGlob >= nProcsCalc ? 1 : 0;

  xmpi(MPI_Comm_group(info->commGlob, &grpAll));

  calcRange[0] = 0;
  calcRange[1] = info->sizeGlob - info->nProcsIO - 1;
  calcRange[2] = 1;
  xmpi(MPI_Group_range_incl(grpAll, info->nProcsIO==info->sizeGlob?0:1,
                            &calcRange, &grpCalc));
  xmpi(MPI_Comm_create(info->commGlob, grpCalc, &info->commCalc));

  pioRange[0] = info->sizeGlob - info->nProcsIO;
  pioRange[1] = info->sizeGlob - 1;
  pioRange[2] = 1;
  xmpi(MPI_Group_range_incl(grpAll, 1, &pioRange, &grpPio));
  xmpi(MPI_Comm_create(info->commGlob, grpPio, &info->commPio));

  if (info->commPio != MPI_COMM_NULL)
    {
      xmpi(MPI_Comm_size(info->commPio, &info->sizePio));
      xmpi(MPI_Comm_rank(info->commPio, &info->rankPio));
      info->sizeCalc = nProcsCalc;
      info->rankCalc = -1;
    }
  else
    {
      info->sizePio = nProcsCalc;
      info->rankPio = -1;
      xmpi(MPI_Comm_size(info->commCalc, &info->sizeCalc));
      xmpi(MPI_Comm_rank(info->commCalc, &info->rankCalc));
    }
  xmpi(MPI_Group_free(&grpCalc));
  xmpi(MPI_Group_free(&grpPio));
  xmpi(MPI_Group_free(&grpAll));
}


MPI_Comm commInqCommPio  ( void )
{
  xassert ( info != NULL &&
           info->commPio != MPI_COMM_NULL &&
           info->isProcIO == 1 );
  return info->commPio;
}


MPI_Comm commInqCommModel ( void )
{
  xassert ( info != NULL &&
           info->commCalc != MPI_COMM_NULL &&
           info->isProcIO == 0 );
  return info->commCalc;
}


int      commInqRankPio       ( void )
{ 
  xassert ( info != NULL &&
           info->rankPio >= 0 &&
           info->isProcIO == 1 );
  return info->rankPio;
}


int      commInqRankModel       ( void )
{ 
  xassert ( info != NULL &&
           info->rankCalc >= 0 &&
           info->isProcIO == 0 );
  return info->rankCalc;
}


void     commDefCommColl  ( int isProcColl )
{
  xassert ( info != NULL &&
           info->commNode != MPI_COMM_NULL &&
           info->commColl == MPI_COMM_NULL );

  info->nodeInfo.isProcColl = isProcColl;
  xmpi ( MPI_Comm_split ( info->commNode, info->nodeInfo.isProcColl, 0, 
                          &info->commColl ));
  xmpi ( MPI_Comm_size ( info->commColl, &info->sizeColl ));
  xmpi ( MPI_Comm_rank ( info->commColl, &info->rankColl )); 
}


MPI_Comm commInqCommColl    ( void )
{
  xassert ( info != NULL &&
           info->commColl != MPI_COMM_NULL );
  return info->commColl;
}


int      commInqSizeColl    ( void )
{
  xassert ( info != NULL &&
           info->sizeColl != CDI_UNDEFID );
  return info->sizeColl; 
}


int      commInqRankColl    ( void )
{
  xassert ( info != NULL &&
           info->rankColl != CDI_UNDEFID );
  return info->rankColl; 
}

static int
cmpstringp(const void *p1, const void *p2)
{
  return strcmp(* (char * const *) p1, * (char * const *) p2);
}

void commDefCommNode ( void )
{
  int size;
  char * myHost, (*allHosts)[MPI_MAX_PROCESSOR_NAME], **sortedHosts;

  xassert ( info != NULL &&
           info->commPio != MPI_COMM_NULL );

  size = info->sizePio;

  myHost = (char*) xmalloc(MPI_MAX_PROCESSOR_NAME);
  {
    int len;
    xmpi ( MPI_Get_processor_name ( myHost, &len ));
    xassert ( myHost[0] != '\0' );
    strncpy(info->hostname, myHost, (size_t)len);
    info->hostname[len] = '\0';
  }

  allHosts = xmalloc((size_t)size * MPI_MAX_PROCESSOR_NAME);
  sortedHosts = xmalloc((size_t)size * sizeof (sortedHosts[0]));

  for (int i = 0; i < size; ++i)
    sortedHosts[i] = allHosts[i];

  xmpi(MPI_Allgather(myHost, MPI_MAX_PROCESSOR_NAME, MPI_CHAR,
                     allHosts[0], MPI_MAX_PROCESSOR_NAME,
                     MPI_CHAR, info->commPio ));

  qsort(sortedHosts, (size_t)size, sizeof (sortedHosts[0]),
        (int (*)(const void *, const void *))cmpstringp);

  {
    int i = 1, nHosts = 1;
    if (!strcmp(myHost, sortedHosts[0]))
      {
        info->nodeInfo.hostID = 1;
      }
    else
      {
        for (; i < size && strcmp(myHost, sortedHosts[i]); ++i)
          if (strcmp(sortedHosts[i - 1], sortedHosts[i]))
            {
              /* new host seen, might be ours */
              nHosts += 1;
              if (!strcmp(sortedHosts[i], myHost))
                info->nodeInfo.hostID = nHosts - 1;
            }
      }
    for (; i < size && strcmp(myHost, sortedHosts[i]); ++i)
      if (strcmp(sortedHosts[i - 1], sortedHosts[i]))
        nHosts += 1;
    info->nodeInfo.nNodes = nHosts;
  }

  xassert ( info->nodeInfo.hostID != CDI_UNDEFID );

  xmpi ( MPI_Comm_split ( info->commPio, info->nodeInfo.hostID, 0, 
                          &info->commNode ));
  xmpi ( MPI_Comm_size ( info->commNode, &info->sizeNode ));
  xmpi ( MPI_Comm_rank ( info->commNode, &info->rankNode ));
  if ( info->IOMode >= PIO_MINIOMODEWITHSPECIALPROCS )
    info->specialRankNode = info->sizeNode - 1;

  free(sortedHosts);
  free(allHosts);
  free(myHost);

  return;
}


MPI_Comm commInqCommNode    ( void )
{
  xassert ( info != NULL &&
           info->commNode != MPI_COMM_NULL );
  return info->commNode;
} 


int commInqSizeNode ( void )
{
  xassert ( info != NULL &&
           info->sizeNode != CDI_UNDEFID );
  return info->sizeNode;
}


int commInqRankNode ( void )
{
  xassert ( info != NULL &&
           info->rankNode != CDI_UNDEFID );
  return info->rankNode;
}


int commInqSpecialRankNode ( void )
{
  xassert ( info != NULL );
  return info->specialRankNode;
}


void commSendNodeInfo ( void )
{
  xassert ( info != NULL &&
           info->root                != CDI_UNDEFID &&
           info->nodeInfo.hostID     != CDI_UNDEFID &&
           info->nodeInfo.isProcColl != CDI_UNDEFID &&
           info->nodeInfo.nNodes     != CDI_UNDEFID &&
           info->commGlob            != MPI_COMM_NULL );
  
  xmpi ( MPI_Send ( &info->nodeInfo, sizeNodeInfo, MPI_INTEGER, 
                    info->root, NODEINFO, info->commGlob ));
}


void     commRecvNodeMap    ( void )
{
  MPI_Status status;
  int source; 

  xassert ( info != NULL &&
           info->commGlob    != MPI_COMM_NULL &&
           info->rankGlob    != CDI_UNDEFID &&
           info->nProcsModel != CDI_UNDEFID &&
           info->nodeMap     == NULL );

  source = info->root;

  xmpi ( MPI_Probe ( source, NODEMAP, info->commGlob, &status ));
  xmpi ( MPI_Get_count ( &status, MPI_INTEGER, &info->nProcsColl ));

  xdebug ( "info->nProcsColl=%d", info->nProcsColl );

  info->nodeMap = (int *)xmalloc((size_t)info->nProcsColl
                                 * sizeof (info->nodeMap[0]));

  xmpi ( MPI_Recv ( info->nodeMap, info->nProcsColl, MPI_INTEGER, 
                    source, NODEMAP, info->commGlob, &status ));
}


void     commDefNNodes      ( int nNodes )
{
  xassert ( info != NULL );
  info->nodeInfo.nNodes = nNodes;
}


int      commInqNNodes      ( void )
{
  xassert ( info != NULL &&
           info->nodeInfo.nNodes != CDI_UNDEFID );
  return info->nodeInfo.nNodes;
}


void     commEvalPhysNodes  ( void )
{
  nodeInfo_t * nodeInfo;
  MPI_Status status;
  int i, IOID, collID, size, nNodes = CDI_UNDEFID;
  int ** p1, ** p2,  idx;

  xassert ( info != NULL &&
           info->root        != CDI_UNDEFID &&
           info->nProcsIO    != CDI_UNDEFID &&
           info->nProcsModel != CDI_UNDEFID &&
           info->sizeGlob    != CDI_UNDEFID &&
           info->rankGlob    != CDI_UNDEFID &&
           info->commGlob    != MPI_COMM_NULL );

  size = info->nProcsIO * sizeNodeInfo;

  nodeInfo = (nodeInfo_t *)xmalloc((size_t)size * sizeof (int));

  if ( info->rankGlob == info->root )
    {
      xassert ( info->rankCalc == info->root );
 
      for ( i = info->nProcsModel; i < info->sizeGlob; i++ )
        {
          IOID = i - info->nProcsModel;
          xmpi ( MPI_Recv ( nodeInfo + IOID, sizeNodeInfo, MPI_INTEGER, i, 
                            NODEINFO, info->commGlob, &status ));
        }
    }

  xmpi ( MPI_Bcast ( nodeInfo, size, MPI_INTEGER, info->root, 
                     info->commCalc ));

  // consistency check, count collectors
  for ( IOID = 0; IOID < info->nProcsIO; IOID++ )
    {
      if ( IOID == 0 )
        {
          xassert ( nodeInfo[IOID].nNodes > 0 &&
                   nodeInfo[IOID].nNodes <= info->nProcsIO );
          nNodes = nodeInfo[IOID].nNodes;
        }
      else xassert ( nodeInfo[IOID].nNodes == nNodes );         
      xassert ( nodeInfo[IOID].hostID > 0 &&
               nodeInfo[IOID].hostID <= nNodes );
      if ( nodeInfo[IOID].isProcColl ) 
        {
          if ( info->nProcsColl == CDI_UNDEFID ) info->nProcsColl = 0;
          info->nProcsColl++;
        }
    }

  xdebug ( "info->nProcsColl=%d", info->nProcsColl );

  xassert ( info->nProcsColl <= info->nProcsModel );

  info->procsCollMap = (int *)xmalloc((size_t)info->nProcsColl
                                      * sizeof (info->procsCollMap[0]));

  // define nodeSizes
  info->nodeInfo.nNodes = nNodes; 
  info->nodeSizes = xcalloc((size_t)info->nodeInfo.nNodes,
                            sizeof (info->nodeSizes[0]));
  collID = 0;
  for ( IOID = 0; IOID < info->nProcsIO; IOID++ )
    if ( nodeInfo[IOID].isProcColl )
      {
        info->nodeSizes[nodeInfo[IOID].hostID - 1]++;
        info->procsCollMap[collID++] = IOID + info->nProcsModel;
      }

  // define nodeMap
  info->nodeMap = (int *)xmalloc((size_t)info->nProcsColl
                                 * sizeof (info->nodeMap[0]));
  // helpers
  p1 = (int **)xmalloc((size_t)info->nodeInfo.nNodes * sizeof (p1[0]));
  p2 = (int **)xmalloc((size_t)info->nodeInfo.nNodes * sizeof (p2[0]));
  idx = 0;
  for ( i = 0; i < info->nodeInfo.nNodes; i++ )
    {
      xassert ( idx >= 0 && idx < info->nProcsColl );
      p1[i] = &info->nodeMap[idx];
      p2[i] = p1[i] + info->nodeSizes[i]; 
      idx += info->nodeSizes[i];
    }
  xassert ( idx == info->nProcsColl );

  // rankGlob in nodeMap
  for ( IOID = 0; IOID < info->nProcsIO; IOID++ )
    {
      if ( nodeInfo[IOID].isProcColl )
        {
          xassert ( p1[nodeInfo[IOID].hostID - 1] <  
                   p2[nodeInfo[IOID].hostID - 1] );
          * p1[nodeInfo[IOID].hostID - 1]++ = IOID + info->nProcsModel;
        }
    }

  if ( info->rankGlob == info->root )
    for ( IOID = info->nProcsModel; IOID < info->sizeGlob; IOID++ )
      xmpi ( MPI_Send ( info->nodeMap, info->nProcsColl, MPI_INTEGER, 
                        IOID, NODEMAP, info->commGlob ));
  free ( p2 );
  free ( p1 );
  free ( nodeInfo );
}


int      commCollID2RankGlob      ( int collID ) 
{
  xassert ( info != NULL &&
           info->nProcsColl != CDI_UNDEFID &&
           collID >= 0 &&
           collID < info->nProcsColl &&
           info->nodeMap != NULL );
  return info->nodeMap[collID];
}


int      commRankGlob2CollID      ( int rankGlob)
{
  int out = CDI_UNDEFID, collID;

  xassert ( info != NULL &&
           info->nProcsColl != CDI_UNDEFID &&
           info->nProcsModel != CDI_UNDEFID &&
           info->sizeGlob != CDI_UNDEFID &&
           rankGlob >= info->nProcsModel &&
           rankGlob <= info->sizeGlob );

  for ( collID = 0; collID < info->nProcsColl; collID++ )
    if ( info->nodeMap[collID] == rankGlob ) out = collID;

  return out;
} 


int *    commInqNodeSizes   ( void )
{
  xassert ( info != NULL &&
           info->nodeSizes != NULL );
  return info->nodeSizes;
}


// collective call
void     commDefCommsIO     ( void )
{
  MPI_Group groupGlob;
  int collID, * ranks, i, currIORank;
  char name[MAXCOMMIONAME];

  xassert ( info != NULL &&
           info->nProcsColl != CDI_UNDEFID &&
           info->commsIO == NULL &&
           info->nProcsModel != CDI_UNDEFID &&
           info->commGlob != MPI_COMM_NULL );

  info->commsIO = (MPI_Comm *)xmalloc((size_t)info->nProcsColl
                                      * sizeof (info->commsIO[0]));
  for ( collID = 0; collID < info->nProcsColl; collID++ )
    info->commsIO[collID] = MPI_COMM_NULL;

  strncpy ( name, "COMMSIO_", 8 );
  name[MAXCOMMIONAME - 1] = '\0';

  ranks = (int *)xmalloc(((size_t)info->nProcsModel + 1) * sizeof (ranks[0]));
  for ( i = 0; i < info->nProcsModel; i++ )
    ranks[i] = i;

  xmpi ( MPI_Comm_group ( info->commGlob, &groupGlob )); 

  for ( collID = 0; collID < info->nProcsColl; collID++ )
    {
      currIORank = info->nodeMap[collID];
      ranks[info->nProcsModel] = currIORank;
      MPI_Group currGroupIO;
      xmpi(MPI_Group_incl(groupGlob, info->nProcsModel + 1,
                          ranks, &currGroupIO));
      xmpi(MPI_Comm_create(info->commGlob, currGroupIO,
                           info->commsIO + collID));
      xmpi(MPI_Group_free(&currGroupIO));
      if ( info->rankGlob == currIORank )
	info->commCalc = info->commsIO[collID];

      // set names for debugging
      if ( info->rankGlob < info->nProcsModel || 
           info->rankGlob == currIORank )
	{
	  sprintf ( &name[8], "%d", collID ); 
	  xmpi ( MPI_Comm_set_name ( info->commsIO[collID], name ));
	}
    }

  if ( ddebug >= 2 )
    {
      if ( info->rankGlob < info->nProcsModel )
	for ( collID = 0; collID < info->nProcsColl; collID++ ) 
          xdebugComm ( info->commsIO + collID );
      else if ( info->nodeInfo.isProcColl )
	xdebugComm ( &info->commCalc );
    }

  xmpi ( MPI_Group_free ( &groupGlob ));
  free ( ranks );

  commPrint ( stdout );
}


MPI_Comm commInqCommsIO     ( int collID )
{
  xassert ( collID >= 0 &&
           info != NULL &&
           info->nProcsColl != CDI_UNDEFID &&
           collID < info->nProcsColl &&
           info->commsIO != NULL &&
           info->commsIO[collID] != MPI_COMM_NULL );

  return info->commsIO[collID];
}


MPI_Comm commInqCommCalc    ( void )
{
  xassert ( info != NULL &&
           info->commCalc != MPI_COMM_NULL );
  return info->commCalc;
}


int      commInqNProcsColl  ( void )
{
  xassert ( info != NULL &&
           info->nProcsColl != CDI_UNDEFID );
  return info->nProcsColl;
}


void commPrint ( FILE * fp )
{
  int i;

  if ( ddebug == 0 ) return;

  xassert ( info != NULL );

  fprintf ( fp, "\n" );
  fprintf ( fp, "######## pioinfo PE%d ###########\n", info->rankGlob );
  fprintf ( fp, "#\n" );
  fprintf ( fp, "# IOMode      = %d\n", info->IOMode );
  fprintf ( fp, "# nProcsIO    = %d\n", info->nProcsIO );
  fprintf ( fp, "# nProcsModel = %d\n", info->nProcsModel );
  fprintf ( fp, "# isProcIO    = %d\n", info->isProcIO );
  fprintf ( fp, "#\n" );
  fprintf ( fp, "# commGlob    = %d\n", (int)MPI_Comm_c2f(info->commGlob));
  fprintf ( fp, "# sizeGlob    = %d\n", info->sizeGlob );
  fprintf ( fp, "# rankGlob    = %d\n", info->rankGlob );
  fprintf ( fp, "# root        = %d\n", info->root );
  fprintf ( fp, "#\n" );
  fprintf ( fp, "# commPio     = %d\n", (int)MPI_Comm_c2f(info->commPio));
  fprintf ( fp, "# sizePio     = %d\n", info->sizePio );
  fprintf ( fp, "# rankPio     = %d\n", info->rankPio );
  fprintf ( fp, "#\n" );
  fprintf ( fp, "# commNode    = %d\n", (int)MPI_Comm_c2f(info->commNode));
  fprintf ( fp, "# sizeNode    = %d\n", info->sizeNode );
  fprintf ( fp, "# rankNode    = %d\n", info->rankNode );
  fprintf ( fp, "# hostname    = %s\n", info->hostname[0] == 0 ? "nn" : 
	    info->hostname );
  fprintf ( fp, "# specialRankNode = %d\n", info->specialRankNode );
  fprintf ( fp, "#\n" );
  fprintf ( fp, "# commColl    = %d\n", (int)MPI_Comm_c2f(info->commColl));
  fprintf ( fp, "# sizeColl    = %d\n", info->sizeColl );
  fprintf ( fp, "# rankColl    = %d\n", info->rankColl );
  fprintf ( fp, "#\n" );
  fprintf ( fp, "# commCalc    = %d\n", (int)MPI_Comm_c2f(info->commCalc));
  if ( info->commsIO != NULL )
    {
      fprintf ( fp, "#\n" );
      for ( i = 0; i < info->nProcsColl; i++ )
	fprintf(fp, "# commsIO[%d]  = %d\n", i,
                (int)MPI_Comm_c2f(info->commsIO[i]));
    }
  else
    fprintf ( fp, "# commsIO     = NULL\n" );
  fprintf ( fp, "# hostID      = %d\n", info->nodeInfo.hostID );
  fprintf ( fp, "# isProcColl  = %d\n", info->nodeInfo.isProcColl );
  fprintf ( fp, "# nNodes      = %d\n", info->nodeInfo.nNodes );
  fprintf ( fp, "# nProcsColl  = %d\n", info->nProcsColl );
  if ( info->procsCollMap != NULL )
    {
      fprintf ( fp, "# procsCollMap   = " );
      xassert ( info->nProcsColl != CDI_UNDEFID );
      for ( i = 0; i < info->nProcsColl; i++ )
        fprintf ( fp, "%d ", info->procsCollMap[i] );
      fprintf ( fp, "\n" );
    }
  else     
    fprintf ( fp, "# procsCollMap   = NULL\n" );
  if ( info->nodeSizes != NULL )
    {
      fprintf ( fp, "# nodeSizes   = " );
      xassert ( info->nodeInfo.nNodes != CDI_UNDEFID );
      for ( i = 0; i < info->nodeInfo.nNodes; i++ )
        fprintf ( fp, "%d ", info->nodeSizes[i] );
      fprintf ( fp, "\n" );
    }
  else     
    fprintf ( fp, "# nodeSizes   = NULL\n" );
  if ( info->nodeMap != NULL )
    {
      fprintf ( fp, "# nodeMap     = " );
      xassert ( info->nProcsColl != CDI_UNDEFID );
      for ( i = 0; i < info->nProcsColl; i++ )
        fprintf ( fp, "%d ", info->nodeMap[i] );
      fprintf ( fp, "\n" );
    }
  else     
    fprintf ( fp, "# nodeMap     = NULL\n" );
  fprintf ( fp, "#\n" );
  fprintf ( fp, "############################\n" );
  fprintf ( fp, "\n" );
}

/************************************************************************/
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
