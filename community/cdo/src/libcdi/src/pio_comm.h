#ifndef PIO_INFO_
#define PIO_INFO_

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <mpi.h>
#include <stdio.h>


typedef struct {
  int hostID;
  int isProcColl;
  int nNodes;
} nodeInfo_t;

enum { MAXNPROCSIO   = 256,
       MAXCOMMIONAME = 12 };

enum { sizeNodeInfo  = 3,
       NODEINFO      = 1111, 
       NODEMAP       = 2222};

void     commInit               ( void );
void     commDestroy            ( void );

void     commDefCommGlob        ( MPI_Comm ); 
MPI_Comm commInqCommGlob        ( void );
int      commInqSizeGlob        ( void );
int      commInqRankGlob        ( void );
int      commInqRootGlob        ( void );

void     commDefNProcsIO        ( int );
int      commInqNProcsIO        ( void );
int      commInqNProcsModel     ( void );
int      commInqIsProcIO        ( void );
void     commDefIOMode          ( int );
int      commInqIOMode          ( void );

void     commDefCommPio         ( void );
MPI_Comm commInqCommPio         ( void );
MPI_Comm commInqCommModel       ( void );
int      commInqRankPio         ( void );
int      commInqRankModel       ( void );

void     commDefCommNode        ( void );
MPI_Comm commInqCommNode        ( void );
int      commInqSizeNode        ( void );
int      commInqRankNode        ( void );
int      commInqSpecialRankNode ( void );

void     commDefCommColl        ( int ); 
MPI_Comm commInqCommColl        ( void );
int      commInqSizeColl        ( void ); 
int      commInqRankColl        ( void ); 
                                
void     commSendNodeInfo       ( void );
void     commRecvNodeMap        ( void );// todo switch to gatherNodeInfo inside commpio
int *    commInqNodeSizes       ( void );
int      commInqNNodes          ( void );
void     commEvalPhysNodes      ( void );
int      commCollID2RankGlob    ( int );
int      commRankGlob2CollID    ( int );
void     commDefCommsIO         ( void ); 
MPI_Comm commInqCommsIO         ( int );
MPI_Comm commInqCommCalc        ( void );
int      commInqNProcsColl      ( void );

void     commPrint              ( FILE * );

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
