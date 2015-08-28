#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifdef USE_MPI
#include <mpi.h>
#include <yaxt.h>
#else
typedef int MPI_Comm;
#endif

#include "cdi.h"
#ifdef USE_MPI
#include "cdipio.h"
#include "pio_util.h"
#ifdef HAVE_PPM_CORE
#include <ppm/ppm_uniform_partition.h>
#endif
#endif

#include "cksum.h"
#include "dmemory.h"
#include "error.h"
#include "pio_write.h"

#include "simple_model_helper.h"
#include "create_uuid.h"

enum {
  ntfiles     = 2,
};


static void
modelRegionCompute(double region[], size_t offset, size_t len,
                   int nlev, int nlat, int nlon,
                   int tsID, const double lons[], const double lats[],
                   double mscale, double mrscale)
{
  size_t local_pos;
  (void)nlev;
  for (local_pos = 0; local_pos < len; ++local_pos)
    {
      size_t global_pos = offset + local_pos;
      int k = (int)(global_pos / (size_t)(nlon * nlat));
      int j = (int)(global_pos % (size_t)(nlon * nlat) / (size_t)nlon);
      int i = (int)(global_pos % (size_t)nlon);
      region[local_pos]
        = dg_wobble((lons[(i + tsID)%nlon] - lons[0])
                    / (lons[nlon-1] - lons[0]),
                    (lats[(j + k)%nlat] - lats[0])
                    / (lats[nlat-1] - lats[0]), mscale, mrscale);
    }
}

void
modelRun(struct model_config setup, MPI_Comm comm)
{
  static const char * const fname_prefix        = "example";

  struct
  {
    size_t size;
    int nlev, zaxisID, id, code;
    uint32_t checksum_state;
#if USE_MPI
    int chunkSize, start;
    Xt_idxlist partDesc;
#endif
  } *varDesc;
  int gridID, taxisID, vlistID, streamID, tsID, tfID = 0;
  int i, nmiss = 0;
  double *lons, *lats, *levs;
  double *var = NULL, *varslice = NULL;
  double mscale, mrscale;
  time_t current_time;
  int vdate = 19850101, vtime = 120000;
  int rank = 0;
  char filename[1024];
  int nlon = setup.nlon, nlat = setup.nlat;
  int nVars = setup.nvars;
  size_t varslice_size = 0;
#if USE_MPI
  int *chunks = NULL, *displs = NULL, comm_size = 1;
#endif

#if USE_MPI
  xmpi ( MPI_Comm_rank ( comm, &rank ));
  xmpi ( MPI_Comm_size ( comm, &comm_size ));
  if (rank == 0 && setup.compute_checksum)
    {
      chunks = xmalloc((size_t)comm_size * sizeof (chunks[0]));
      displs = xmalloc((size_t)comm_size * sizeof (displs[0]));
      var = xmalloc((size_t)nlon * (size_t)nlat
                    * (size_t)setup.max_nlev * sizeof(var[0]));
    }
#endif

  var_scale(setup.datatype, &mscale, &mrscale);

  gridID = gridCreate ( GRID_LONLAT, nlon*nlat );
  gridDefXsize ( gridID, nlon );
  gridDefYsize ( gridID, nlat );
  lons = xmalloc((size_t)nlon * sizeof (lons[0]));
  for (i = 0; i < nlon; ++i)
    lons[i] = ((double)(i * 360))/nlon;
  lats = xmalloc((size_t)nlat * sizeof (lats[0]));
  for (i = 0; i < nlat; ++i)
    lats[i] = ((double)(i * 180))/nlat - 90.0;
  gridDefXvals ( gridID, lons );
  gridDefYvals ( gridID, lats );
  {
    unsigned char uuid[CDI_UUID_SIZE];
    if (rank == 0)
      create_uuid(uuid);
#if USE_MPI
    MPI_Bcast(uuid, CDI_UUID_SIZE, MPI_UNSIGNED_CHAR, 0, comm);
#endif
    gridDefUUID(gridID, uuid);
  }
  levs = xmalloc((size_t)setup.max_nlev * sizeof (levs[0]));
  {
    double lscale = 1.0/(double)(setup.max_nlev - 1);
    for (i = 0; i < setup.max_nlev; ++i)
      levs[i] = 101300.0 - 13000.0 * expm1(2.173 * (double)i * lscale);
  }
  vlistID = vlistCreate ();

  varDesc = xmalloc((size_t)nVars * sizeof (varDesc[0]));
  for (int varIdx = 0; varIdx < nVars; varIdx++ )
    {
      int varLevs = (int)random()%4;
      switch (varLevs)
        {
        case 1:
          varLevs = setup.max_nlev / 3;
          break;
        case 2:
          varLevs = setup.max_nlev >= 11 ? 11 : setup.max_nlev / 2;
          break;
        case 3:
          varLevs = setup.max_nlev - 1;
          break;
        }
      ++varLevs;
      varDesc[varIdx].nlev = varLevs;
      for (size_t i = 0; i < (size_t)varIdx; ++i)
        if (varDesc[i].nlev == varLevs)
          {
            varDesc[varIdx].zaxisID = varDesc[i].zaxisID;
            goto zaxisIDset;
          }
      if (varLevs == 1)
        varDesc[varIdx].zaxisID = zaxisCreate(ZAXIS_SURFACE, 1);
      else
        {
          varDesc[varIdx].zaxisID
            = zaxisCreate(ZAXIS_PRESSURE, varDesc[varIdx].nlev);
          zaxisDefLevels(varDesc[varIdx].zaxisID, levs);
        }
      {
        unsigned char uuid[16];
        if (rank == 0)
          create_uuid(uuid);
#if USE_MPI
        MPI_Bcast(uuid, CDI_UUID_SIZE, MPI_UNSIGNED_CHAR, 0, comm);
#endif
        zaxisDefUUID(varDesc[varIdx].zaxisID, uuid);
      }
      zaxisIDset:
      varDesc[varIdx].id
        = vlistDefVar(vlistID, gridID, varDesc[varIdx].zaxisID, TIME_VARIABLE);
      varDesc[varIdx].size
        = (size_t)nlon * (size_t)nlat * (size_t)varDesc[varIdx].nlev;
#ifdef USE_MPI
      {
        struct PPM_extent range
          = PPM_uniform_partition((struct PPM_extent){ 0,
                (int32_t)varDesc[varIdx].size }, comm_size, rank);
        int start = range.first;
        int chunkSize = range.size;
        Xt_idxlist idxlist
          = xt_idxstripes_new(&(struct Xt_stripe){ .start = start,
                .nstrides = chunkSize, .stride = 1 }, 1);
        varDesc[varIdx].start = start;
        varDesc[varIdx].chunkSize = chunkSize;
        varDesc[varIdx].partDesc = idxlist;
      }
#endif
      varDesc[varIdx].code = 129 + varIdx;
      vlistDefVarCode(vlistID, varDesc[varIdx].id, varDesc[varIdx].code);
      vlistDefVarDatatype(vlistID, varDesc[varIdx].id, setup.datatype);
    }

  taxisID = taxisCreate ( TAXIS_ABSOLUTE );
  vlistDefTaxis ( vlistID, taxisID );

  sprintf ( &filename[0], "%s_%d.%s", fname_prefix, tfID, setup.suffix );
  streamID = streamOpenWrite ( filename, setup.filetype );
  xassert ( streamID >= 0 );
  streamDefVlist ( streamID, vlistID);

#ifdef USE_MPI
  pioEndDef ();
#endif

  for ( tfID = 0; tfID < ntfiles; tfID++ )
    {
      for (int varIdx = 0; varIdx < nVars; ++varIdx)
        varDesc[varIdx].checksum_state = 0;
      if ( tfID > 0 )
	{
	  streamClose ( streamID );
	  sprintf ( &filename[0], "%s_%d.%s", fname_prefix, tfID, setup.suffix );
	  streamID = streamOpenWrite ( filename, setup.filetype );
	  xassert ( streamID >= 0 );
	  streamDefVlist ( streamID, vlistID );
	}
      vdate = 19850101;
      vtime = 120000;
      current_time = cditime2time_t(vdate, vtime);
      for ( tsID = 0; tsID < setup.nts; tsID++ )
	{
          int vdatetime[2];
          time_t2cditime(current_time, &vdatetime[1], &vdatetime[0]);
	  taxisDefVdate(taxisID, vdatetime[1]);
	  taxisDefVtime(taxisID, vdatetime[0]);
	  streamDefTimestep ( streamID, tsID );
          if (setup.filetype == FILETYPE_EXT)
            {
              /* EXTRA doesn't store time, only date
               * set the value to 0 before checksumming, because a
               * time field of 0 is what reading an EXTRA file will
               * return */
              vdatetime[0] = 0;
            }
	  for (int varID = 0; varID < nVars; ++varID)
	    {
#ifdef USE_MPI
              int start = varDesc[varID].start;
              int chunk = varDesc[varID].chunkSize;
#else
              int chunk = (int)varDesc[varID].size;
              int start = 0;
#endif
              if (varslice_size < (size_t)chunk)
                {
                  varslice = xrealloc(varslice, (size_t)chunk * sizeof (var[0]));
                  varslice_size = (size_t)chunk;
                }
              modelRegionCompute(varslice, (size_t)start, (size_t)chunk,
                                 varDesc[varID].nlev, nlat, nlon,
                                 tsID, lons, lats,
                                 mscale, mrscale);
              if (setup.compute_checksum)
                {
#if USE_MPI
                  xmpi(MPI_Gather(&chunk, 1, MPI_INT,
                                  chunks, 1, MPI_INT, 0, comm));
                  if (rank == 0)
                    {
                      displs[0] = 0;
                      for (i = 1; i < comm_size; ++i)
                        displs[i] = displs[i - 1] + chunks[i - 1];
                    }
                  xmpi(MPI_Gatherv(varslice, chunk, MPI_DOUBLE,
                                   var, chunks, displs, MPI_DOUBLE, 0, comm));
#else
                  var = varslice;
#endif
                }
              if (rank == 0 && setup.compute_checksum)
                {
                  memcrc_r(&varDesc[varID].checksum_state,
                           (const unsigned char *)vdatetime,
                           sizeof (vdatetime));
                  memcrc_r(&varDesc[varID].checksum_state,
                           (const unsigned char *)var,
                           varDesc[varID].size * sizeof (var[0]));
                }

#ifdef USE_MPI
	      streamWriteVarPart(streamID, varDesc[varID].id, varslice, nmiss,
                                 varDesc[varID].partDesc);
#else
	      streamWriteVar(streamID, varDesc[varID].id, varslice, nmiss);
#endif
	    }
          current_time += 86400;
#ifdef USE_MPI
	  pioWriteTimestep ( tsID, vdate, vtime );
#endif
	}
      if (rank == 0 && setup.compute_checksum)
        {
          FILE *tablefp;
          {
            sprintf(filename, "%s_%d.cksum", fname_prefix, tfID);
            if (!(tablefp = fopen(filename, "w")))
              {
                perror("failed to open table file");
                exit(EXIT_FAILURE);
              }
            for (i = 0; i < nVars; ++i)
              {
                uint32_t cksum;
                int code;
                cksum = memcrc_finish(&varDesc[i].checksum_state,
                                      (off_t)((varDesc[i].size
                                               * sizeof (var[0])
                                               + sizeof (int) * 2)
                                              * (size_t)setup.nts));
                code = vlistInqVarCode(vlistID, varDesc[i].id);
                if (fprintf(tablefp, "%08lx %d\n", (unsigned long)cksum,
                            code) < 0)
                  {
                    perror("failed to write table file");
                    exit(EXIT_FAILURE);
                  }
              }
            fclose(tablefp);
          }
        }
    }
  free(varslice);
#ifdef USE_MPI
  pioEndTimestepping ();
#endif
  streamClose ( streamID );
  vlistDestroy ( vlistID );
  taxisDestroy ( taxisID );
  for (int varID = 0; varID < nVars; varID++ )
    {
      int zID = varDesc[varID].zaxisID;
      if (zID != CDI_UNDEFID)
        {
          zaxisDestroy(zID);
          for (int j = varID + 1; j < nVars; ++j)
            if (zID == varDesc[j].zaxisID)
              varDesc[j].zaxisID = CDI_UNDEFID;
        }
    }
  gridDestroy ( gridID );
#if USE_MPI
  for (int varID = 0; varID < nVars; ++varID)
    xt_idxlist_delete(varDesc[varID].partDesc);
  free(displs);
  free(chunks);
  free(var);
#endif
  free(varDesc);
  free(levs);
  free(lats);
  free(lons);
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
