#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef USE_MPI
#include <unistd.h>
#include <mpi.h>
#include <yaxt.h>

#include "pio_util.h"
#else
typedef int MPI_Comm;
#define MPI_COMM_NULL 0
#endif

#include "cdi.h"
#include "error.h"

#ifdef USE_MPI
#include "cdipio.h"

static int
uniform_partition_start(int set_interval[2], int nparts, int part_idx);
#endif

static void modelRun(MPI_Comm commModel)
{
  enum {
    filetype    = FILETYPE_GRB,
    ntfiles     = 2,
    ntsteps     = 3,
    nVars       = 5,
    nlon        = 12,
    nlat        = 6,
    maxlev      = 5 };

  static int nlev[nVars]    = {1,1,5,5,2};
  static char * name        = "example";

  int gridID, zaxisID[nVars], taxisID;
  int vlistID, varID[nVars], streamID, tsID, tfID = 0;
  int i, j, nmiss = 0;
  double lons[nlon] = {0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330};
  double lats[nlat] = {-75, -45, -15, 15, 45, 75};
  double levs[maxlev] = {101300, 92500, 85000, 50000, 20000};
  double var[nlon*nlat*maxlev];
  int vdate = 19850101, vtime = 120000;
  char filename[1024];
  size_t varSize[nVars];
#if USE_MPI
  int rank, comm_size;
  struct var1DDeco {
    int chunkSize, start;
    Xt_idxlist partDesc;
  } varDeco[nVars];
#endif

  gridID = gridCreate ( GRID_LONLAT, nlon*nlat );
  gridDefXsize ( gridID, nlon );
  gridDefYsize ( gridID, nlat );
  gridDefXvals ( gridID, lons );
  gridDefYvals ( gridID, lats );

  for ( i = 0; i < nVars; i++ )
    {
      zaxisID[i] = zaxisCreate ( ZAXIS_PRESSURE, nlev[i] );
      zaxisDefLevels ( zaxisID[i], levs );
      varSize[i] = nlon * nlat * (size_t)nlev[i];
    }
  vlistID = vlistCreate ();

#if USE_MPI
  xmpi ( MPI_Comm_rank ( commModel, &rank ));
  xmpi ( MPI_Comm_size ( commModel, &comm_size ));
#endif
  for ( i = 0; i < nVars; i++ )
    {
      varID[i] = vlistDefVar ( vlistID, gridID, zaxisID[i], TIME_VARIABLE);
#ifdef USE_MPI
      {
        int start = uniform_partition_start((int [2]){ 0, (int)varSize[i] - 1 },
                                             comm_size, rank),
          chunkSize = uniform_partition_start((int [2]){ 0, (int)varSize[i] - 1 },
                                              comm_size, rank + 1) - start;
        fprintf(stderr, "%d: start=%d, chunkSize = %d\n", rank,
                start, chunkSize);
        Xt_idxlist idxlist
          = xt_idxstripes_new(&(struct Xt_stripe){ .start = start,
                .nstrides = chunkSize, .stride = 1 }, 1);
        varDeco[i] = (struct var1DDeco){
          .start = start,
          .chunkSize = chunkSize,
          .partDesc = idxlist
        };
      }
#endif
    }
  taxisID = taxisCreate ( TAXIS_ABSOLUTE );
  vlistDefTaxis ( vlistID, taxisID );

  sprintf ( &filename[0], "%s_%d.grb", name, tfID );
  streamID = streamOpenWrite ( filename, filetype );
  xassert ( streamID >= 0 );
  streamDefVlist ( streamID, vlistID);

#ifdef USE_MPI
  pioEndDef ();
#endif
  for ( tfID = 0; tfID < ntfiles; tfID++ )
    {
      /* if ( tfID > 0 ) */
	{
	  streamClose ( streamID );
	  sprintf ( &filename[0], "%s_%d.grb", name, tfID );
	  streamID = streamOpenWrite ( filename, filetype );
	  xassert ( streamID >= 0 );
	  streamDefVlist ( streamID, vlistID );
	}
      for ( tsID = 0; tsID < ntsteps; tsID++ )
	{
	  taxisDefVdate ( taxisID, vdate );
	  taxisDefVtime ( taxisID, vtime );
	  streamDefTimestep ( streamID, tsID );
	  for ( i = 0; i < nVars; i++ )
	    {
#ifdef USE_MPI
              int chunk = varDeco[i].chunkSize;
#else
              int chunk = (int)varSize[i];
#endif
	      for (j = 0; j < chunk; ++j) var[j] = 2.2;
#ifdef USE_MPI
              streamWriteVarPart(streamID, varID[i], var, nmiss,
                                 varDeco[i].partDesc);
#else
	      streamWriteVar ( streamID, varID[i], var, nmiss );
#endif
	    }
#ifdef USE_MPI
	  pioWriteTimestep();
#endif
	}
    }
#ifdef USE_MPI
  pioEndTimestepping ();
#endif
  streamClose ( streamID );
  vlistDestroy ( vlistID );
  taxisDestroy ( taxisID );
  for ( i = 0; i < nVars; i++ )
    zaxisDestroy ( zaxisID[i] );
  gridDestroy ( gridID );
#ifdef USE_MPI
  for (int varID = 0; varID < nVars; ++varID)
    xt_idxlist_delete(varDeco[varID].partDesc);
  MPI_Barrier(commModel);
#endif
}

#ifdef USE_MPI
static struct {
  char *text;
  int mode;
} mode_map[] = {
  { "PIO_MPI", PIO_MPI },
  { "PIO_FPGUARD", PIO_FPGUARD },
  { "PIO_ASYNCH", PIO_ASYNCH },
  { "PIO_WRITER", PIO_WRITER },
  { "PIO_FPGUARD", PIO_FPGUARD},
};

static inline int
search_iomode_str(const char *modestr)
{
  int retval = -1;
  for (size_t i = 0;
       i < sizeof (mode_map) / sizeof (mode_map[0]);
       ++i)
    if (!strcmp(modestr, mode_map[i].text))
      {
        retval = (int)i;
        break;
      }
  return retval;
}
#endif


int main (int argc, char *argv[])
{
  MPI_Comm commModel = MPI_COMM_NULL;
#ifdef USE_MPI
  MPI_Comm commGlob;
  int sizeGlob, rankGlob, pioNamespace;
  int IOMode = PIO_WRITER;
  int nProcsIO = 2;

  xmpi ( MPI_Init ( &argc, &argv));
  commGlob = MPI_COMM_WORLD;
  xt_initialize(commGlob);
  xmpi ( MPI_Comm_set_errhandler ( commGlob, MPI_ERRORS_RETURN ));
  xmpi ( MPI_Comm_size ( commGlob, &sizeGlob ));
  xmpi ( MPI_Comm_rank ( commGlob, &rankGlob ));

  {
    int opt;
    while ((opt = getopt(argc, argv, "p:w:")) != -1)
      switch (opt) {
      case 'p':
        {
          int entry = search_iomode_str(optarg);
          if (entry < 0)
            {
              fprintf(stderr, "Unsupported PIO mode requested: %s\n", optarg);
              exit(EXIT_FAILURE);
            }
          IOMode = mode_map[entry].mode;
        }
        break;
      case 'w':
        {
          long temp = strtol(optarg, NULL, 0);
          if (temp < 0 || temp > INT_MAX/2)
            {
              fprintf(stderr, "Unsupported number of I/O servers: %ld\n", temp);
              exit(EXIT_FAILURE);
            }
          nProcsIO = (int)temp;
          break;
        }
      }
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
