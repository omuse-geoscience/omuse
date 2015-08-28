#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <stdio.h>
#include <string.h>

#ifdef USE_MPI
#include <mpi.h>
#else
typedef int MPI_Comm;
#define MPI_COMM_NULL 0
#endif

#include "cdi.h"
#include "error.h"

#ifdef USE_MPI
#include "cdipio.h"
#include "pio_util.h"
#endif

static void hoursPassingHack ( int * vdate, int * vtime, int hoursPassed )
{
  int sum, days, hours, oldDays;

  xassert ( hoursPassed % 10000 == 0 );
  sum = * vtime + hoursPassed;
  days = sum / 240000;
  hours = sum % 240000;
  oldDays = * vdate % 100;
  if ( oldDays + days > 28 ) xabort ( "UNEXPECTED USE" );
  * vtime = hours;
  * vdate = * vdate + days;
}

#ifdef USE_MPI
static int
uniform_partition_start(int set_interval[2], int nparts, int part_idx);
#endif

static void modelRun(MPI_Comm commModel)
{
  enum {
    filetype    = FILETYPE_GRB,
    nStreams    = 5,
    MAXNSTREAMS = 25,
    ntfiles     = 2,
    ntsteps     = 3,
    nVars       = 5,
    nlon        = 12,
    nlat        = 6,
    MAXNLEV     = 5 };

  static int nlev[nStreams][nVars] =
    {{1,1,5,5,2},{3,5,2,2,1},{3,5,2,2,1},{5,2,2,2,1}, {3,3,3,3,3}};

  static char nameExp[] = "example";
  static int asciiA     = 65;
  char filename[1024];

  int gridID, zaxisID[nStreams][nVars], taxisID;
  int streamID[nStreams], vlistID[nStreams], varID[nStreams][nVars], tsID, tfID = 0;
  int i, j, nmiss = 0;
  double lons[nlon] = {0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330};
  double lats[nlat] = {-75, -45, -15, 15, 45, 75};
  double levs[MAXNLEV] = {101300, 92500, 85000, 50000, 20000};
  double *var;
  int vdate = 19850101, vtime = 120000, hourStep = 20000;
  size_t varSize[nStreams][nVars];
#if USE_MPI
  int rank, comm_size;
  struct var1DDeco {
    int chunkSize, start;
    Xt_idxlist partDesc;
  } varDeco[nStreams][nVars];
#endif
  xassert ( nStreams < MAXNSTREAMS );

  gridID = gridCreate ( GRID_LONLAT, nlon*nlat );
  gridDefXsize ( gridID, nlon );
  gridDefYsize ( gridID, nlat );
  gridDefXvals ( gridID, lons );
  gridDefYvals ( gridID, lats );

  for ( i = 0; i < nStreams; i++ )
    {
      for ( j = 0; j < nVars; j++ )
        {
          xassert ( nlev[i][j] > 0 && nlev[i][j] <= MAXNLEV );
          zaxisID[i][j] = zaxisCreate ( ZAXIS_PRESSURE, nlev[i][j] );
          zaxisDefLevels ( zaxisID[i][j], levs );
        }
      vlistID[i] = vlistCreate ();
    }

#if USE_MPI
  xmpi ( MPI_Comm_rank ( commModel, &rank ));
  xmpi ( MPI_Comm_size ( commModel, &comm_size ));
#endif
  {
    int maxChunkSize = 0;
    for ( i = 0; i < nStreams; i++ )
      for ( j = 0; j < nVars; j++ )
      {
        varID[i][j] = vlistDefVar(vlistID[i], gridID, zaxisID[i][j],
                                  TIME_VARIABLE);
        varSize[i][j] = nlon * nlat * (size_t)nlev[i][j];
#ifdef USE_MPI
        {
          int start = uniform_partition_start((int [2]){ 0, (int)varSize[i][j] - 1 },
                                              comm_size, rank),
            chunkSize = uniform_partition_start((int [2]){ 0, (int)varSize[i][j] - 1 },
                                                comm_size, rank + 1) - start;
          if (maxChunkSize < chunkSize)
            maxChunkSize = chunkSize;
          fprintf(stderr, "%d: start=%d, chunkSize = %d\n", rank,
                  start, chunkSize);
          Xt_idxlist idxlist
            = xt_idxstripes_new(&(struct Xt_stripe){ .start = start,
                  .nstrides = chunkSize, .stride = 1 }, 1);
          varDeco[i][j] = (struct var1DDeco){
            .start = start,
            .chunkSize = chunkSize,
            .partDesc = idxlist
          };
        }
#else
        if (maxChunkSize < varSize[i][j])
          maxChunkSize = (int)varSize[i][j];
#endif
      }
    var = (double *)malloc((size_t)maxChunkSize * sizeof (var[0]));
  }
  taxisID = taxisCreate ( TAXIS_ABSOLUTE );
  for ( i = 0; i < nStreams; i++ )
    vlistDefTaxis ( vlistID[i], taxisID );

  for ( i = 0; i < nStreams; i++ )
    {
      memset ( filename, 0, 1024 );
      sprintf ( &filename[0], "%s%c_%d.grb", nameExp, asciiA + i, tfID );
      streamID[i] = streamOpenWrite ( filename, filetype );
      xassert ( streamID[i] >= 0 );
      streamDefVlist ( streamID[i], vlistID[i]);
    }

#ifdef USE_MPI
  pioEndDef ();
#endif

  for ( tfID = 0; tfID < ntfiles; tfID++ )
    {
      if ( tfID > 0 )
	{
          for ( i = 0; i < nStreams; i++ )
            {
              streamClose ( streamID[i] );
              sprintf ( &filename[0], "%s%c_%d.grb", nameExp, asciiA + i, tfID );
              streamID[i] = streamOpenWrite ( filename, filetype );
              xassert ( streamID[i] >= 0 );
              streamDefVlist ( streamID[i], vlistID[i] );
            }
	}

      for ( tsID = 0; tsID < ntsteps; tsID++ )
	{
          hoursPassingHack ( &vdate, &vtime, hourStep );
	  taxisDefVdate ( taxisID, vdate );
	  taxisDefVtime ( taxisID, vtime );
          /* temporary fix for problem calling streamDefTimestep after
           * streamWriteVarPart
           * FIXME: this can be merged with the loop below again once
           * per-stream RDMA windows are realized */
          for ( i = 0; i < nStreams; i++ )
            streamDefTimestep(streamID[i], tsID);
          for ( i = 0; i < nStreams; i++ )
            {
              for ( j = 0; j < nVars; j++ )
                {
#ifdef USE_MPI
                  int start = varDeco[i][j].start;
                  int chunk = varDeco[i][j].chunkSize;
#else
                  int start = 0, chunk = (int)varSize[i][j];
#endif
                  for(int k = 0; k < chunk; k++)
                    var[k] = 3.3 * (double)(k + start);
#ifdef USE_MPI
                  streamWriteVarPart(streamID[i], varID[i][j], var, nmiss,
                                     varDeco[i][j].partDesc);
#else
                  streamWriteVar(streamID[i], varID[i][j], var, nmiss);
#endif
                }
            }
#ifdef USE_MPI
	  pioWriteTimestep();
#endif
	}
#ifdef USE_MPI
      MPI_Barrier(commModel);
#endif
    }

#ifdef USE_MPI
  pioEndTimestepping ();
  for (int streamID = 0; streamID < nStreams; ++streamID)
    for (int varID = 0; varID < nVars; ++varID)
      xt_idxlist_delete(varDeco[streamID][varID].partDesc);
#endif

  for ( i = 0; i < nStreams; i++ )
    {
      streamClose ( streamID[i] );
      vlistDestroy ( vlistID[i] );
    }
  taxisDestroy ( taxisID );
  for ( i = 0; i < nStreams; i++ )
    for ( j = 0; j < nVars; j++ )
      zaxisDestroy ( zaxisID[i][j] );
  gridDestroy ( gridID );
  free(var);
}


int main (int argc, char *argv[])
{
  enum {
    nProcsIODef    = 3,
    //IOModeDef       = PIO_NONE,
    //IOModeDef       = PIO_MPI,
#ifdef USE_MPI
    IOModeDef       = PIO_FPGUARD,
#endif
    //IOModeDef       = PIO_ASYNCH,
    //IOModeDef       = PIO_WRITER,
  };

  MPI_Comm commModel = MPI_COMM_NULL;
#ifdef USE_MPI
  MPI_Comm commGlob;
  int sizeGlob;
  int rankGlob;
  int IOMode, nProcsIO, pioNamespace;

  xmpi ( MPI_Init ( &argc, &argv));
  commGlob = MPI_COMM_WORLD;
  xt_initialize(commGlob);
  xmpi ( MPI_Comm_set_errhandler ( commGlob, MPI_ERRORS_RETURN ));
  xmpi ( MPI_Comm_size ( commGlob, &sizeGlob ));
  xmpi ( MPI_Comm_rank ( commGlob, &rankGlob ));

  if (argc > 1)
    {
      xassert ( argc >= 3 );
      IOMode = atoi(argv[1]);
      nProcsIO = atoi(argv[2]);
      xassert ( IOMode >= PIO_MINIOMODE && IOMode <= PIO_MAXIOMODE &&
                nProcsIO >= 0 && nProcsIO <= sizeGlob );
      printf ( "IOMode=%d, nProcsIO=%d", IOMode, nProcsIO );
    }
  else
    {
      IOMode = IOModeDef;
      nProcsIO = nProcsIODef;
    }

  commModel = pioInit(commGlob, nProcsIO, IOMode, &pioNamespace, 1.0f,
                      cdiPioNoPostCommSetup);
  if (commModel != MPI_COMM_NULL)
    {
      namespaceSetActive(pioNamespace);
#endif

      modelRun(commModel);

#ifdef USE_MPI
    }
  pioFinalize ();
  xt_finalize();
  MPI_Finalize ();
#endif
  return 0;
}

#ifdef USE_MPI
static int
uniform_partition_start(int set_interval[2], int nparts, int part_idx)
{
  int part_offset
    = (int)((((long long)set_interval[1] - (long long)set_interval[0] + 1LL)
             * (long long)part_idx) / (long long)nparts);
  int start = set_interval[0] + part_offset;
  return start;
}
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
