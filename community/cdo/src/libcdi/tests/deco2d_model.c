#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <inttypes.h>
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
#include <core/ppm_combinatorics.h>
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
  nproma = 16,
};

static void
modelRegionCompute(double region[], int nlev, int nlat, int nlon,
                   const int chunkStart[3], const int chunkSize[3],
                   int tsID, const double lons[], const double lats[],
                   double mscale, double mrscale)
{
  (void)nlev;
  unsigned is = (unsigned)chunkStart[0],
    js = (unsigned)chunkStart[1],
    ks = (unsigned)chunkStart[2],
    m = (unsigned)chunkSize[0],
    n = (unsigned)chunkSize[1],
    o = (unsigned)chunkSize[2],
    jstride = (unsigned)chunkSize[0],
    kstride = ((jstride * (unsigned)chunkSize[1] + nproma - 1)/nproma)*nproma;

  for (unsigned k = 0; k < o; ++k)
    for (unsigned j = 0; j < n; ++j)
      for (unsigned i = 0; i < m; ++i)
        region[k * kstride + j * jstride + i]
          = dg_wobble((lons[(i + is + (unsigned)tsID)%(unsigned)nlon] - lons[0])
                      / (lons[nlon-1] - lons[0]),
                      (lats[(j + js + k + ks)%(unsigned)nlat] - lats[0])
                      / (lats[nlat-1] - lats[0]),
                      mscale, mrscale);
}

#ifdef USE_MPI
static void
findPartition2D(int npart[2], int num_parts);
#endif

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
    int chunkSize[2], start[2];
    Xt_idxlist partDesc;
    Xt_redist redist4gather;
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
  int comm_size = 1;
  int npart[2], rank_coord[2];
  int *blk_displ, *blk_lens;
#endif

#if USE_MPI
  xmpi ( MPI_Comm_rank ( comm, &rank ));
  xmpi ( MPI_Comm_size ( comm, &comm_size ));
#endif

  if (rank == 0 && setup.compute_checksum)
    {
      var = xmalloc((size_t)nlon * (size_t)nlat
                    * (size_t)setup.max_nlev * sizeof(var[0]));
    }

#if USE_MPI
  if (comm_size == 1)
    {
      npart[0] = 1;
      npart[1] = 1;
      rank_coord[0] = 0;
      rank_coord[1] = 0;
    }
  else
    {
      findPartition2D(npart, comm_size);
      rank_coord[0] = rank % npart[0],
        rank_coord[1] = rank / npart[0];
    }
  blk_displ = xmalloc((size_t)setup.max_nlev * sizeof (blk_displ[0]) * 2);
  blk_lens = blk_displ + setup.max_nlev;
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
        int start[2], chunkSize[3], varSize[2] = { nlon, nlat };
        for (size_t i = 0; i < 2; ++i)
          {
            struct PPM_extent range
              = PPM_uniform_partition((struct PPM_extent){ 0, varSize[i] },
                                      npart[i], rank_coord[i]);
            start[i] = range.first;
            chunkSize[i] = range.size;
            fprintf(stderr, "%d: start[%zu]=%d, chunkSize[%zu] = %d\n", rank,
                    i, start[i], i, chunkSize[i]);
            varDesc[varIdx].start[i] = range.first;
            varDesc[varIdx].chunkSize[i] = range.size;
          }
        Xt_int varSizeXt[3] = { (Xt_int)nlon, (Xt_int)nlat, (Xt_int)varLevs };
        chunkSize[2] = varLevs;
        Xt_int varStartXt[3] = { start[0], start[1], 0 };
        for (int i = 0; i < varIdx; ++i)
          if (varDesc[i].nlev == varLevs)
            {
              varDesc[varIdx].redist4gather = varDesc[i].redist4gather;
              varDesc[varIdx].partDesc = varDesc[i].partDesc;
              goto gatherRedistSet;
            }
        Xt_idxlist part_idxlist
          = xt_idxsection_new(0, (varLevs > 1 ? 3 : 2), varSizeXt,
                              chunkSize, varStartXt),
          gather_idxlist;
        varDesc[varIdx].partDesc = part_idxlist;
        if (setup.compute_checksum)
          {
            if (rank == 0)
              {
                gather_idxlist
                  = xt_idxstripes_new(&(struct Xt_stripe){.start = 0,
                        .stride = 1, .nstrides = (int)varDesc[varIdx].size }, 1);
              }
            else
              gather_idxlist = xt_idxempty_new();
            Xt_xmap xmap4gather
              = xt_xmap_all2all_new(part_idxlist, gather_idxlist, comm);
            xt_idxlist_delete(gather_idxlist);
            struct Xt_offset_ext *src_blocks = xmalloc((size_t)varLevs
                                                       * sizeof (*src_blocks));
            struct Xt_offset_ext dst_block = { .start = 0,
                                               .size = nlon * nlat * varLevs,
                                               .stride = 1 };
            size_t levStride
              = (((size_t)chunkSize[0] * (size_t)chunkSize[1] + nproma - 1)
                 / nproma) * nproma;
            for (size_t i = 0; i < (size_t)varLevs; ++i)
              src_blocks[i] = (struct Xt_offset_ext)
                { .start = (int)(i * levStride),
                  .size = chunkSize[0] * chunkSize[1],
                  .stride = 1 };
            varDesc[varIdx].redist4gather
              = xt_redist_p2p_ext_new(xmap4gather,
                                      varLevs, src_blocks, 1, &dst_block,
                                      MPI_DOUBLE);
            free(src_blocks);
            xt_xmap_delete(xmap4gather);
          }
        gatherRedistSet: ;
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
	  for (int varID = 0; varID < nVars; ++varID)
	    {
              size_t varLevs = (size_t)varDesc[varID].nlev;
#ifdef USE_MPI
              int start[3] = { varDesc[varID].start[0],
                               varDesc[varID].start[1],
                               0 };
              int chunk[3] = { varDesc[varID].chunkSize[0],
                               varDesc[varID].chunkSize[1],
                               (int)varLevs };
#else
              int chunk[3] = { nlon, nlat, (int)varLevs };
              int start[3] = { 0, 0, 0 };
#endif
              size_t chunkSize
                = (((size_t)chunk[0] * (size_t)chunk[1] + (size_t)(nproma - 1))
                   / (size_t)nproma) * (size_t)nproma * varLevs;
              if (varslice_size < chunkSize)
                {
                  varslice = xrealloc(varslice, chunkSize * sizeof (var[0]));
                  varslice_size = chunkSize;
                }
              modelRegionCompute(varslice, (int)varLevs, nlat, nlon,
                                 start, chunk, tsID, lons, lats,
                                 mscale, mrscale);
              if (setup.compute_checksum)
                {
#if USE_MPI
                  xt_redist_s_exchange1(varDesc[varID].redist4gather,
                                        varslice, var);
                  size_t layerSize = (size_t)(chunk[0] * chunk[1]);
                  size_t nblk = (layerSize + nproma - 1)/nproma - 1;
                  for (size_t k = 0; k < varLevs; ++k)
                    {
                      blk_displ[k] = (int)(k * (nblk + 1) * nproma);
                      blk_lens[k] = (int)layerSize;
                    }
#else
                  size_t layerSize = (size_t)(chunk[0] * chunk[1]);
                  size_t nblk = (layerSize + nproma - 1)/nproma - 1;
                  size_t npromz = layerSize - nblk * nproma;
                  for (size_t k = 0; k < varLevs; ++k)
                    {
                      for (size_t j = 0; j < nblk; ++j)
                        for (size_t i = 0; i < nproma; ++i)
                          var[k * layerSize + j * nproma + i] =
                            varslice[k * (nblk + 1) * nproma + j * nproma + i];
                      for (size_t i = 0; i < npromz; ++i)
                        var[k * layerSize + nblk * nproma + i] =
                          varslice[k * (nblk + 1) * nproma + nblk * nproma + i];
                    }
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
	      streamWriteScatteredVarPart(streamID, varDesc[varID].id,
                                          varslice,
                                          (int)varLevs, blk_lens, blk_displ,
                                          nmiss, varDesc[varID].partDesc);
#else
	      streamWriteVar(streamID, varDesc[varID].id, var, nmiss);
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
#if USE_MPI
          xt_idxlist_delete(varDesc[varID].partDesc);
          xt_redist_delete(varDesc[varID].redist4gather);
#endif
          for (int j = varID + 1; j < nVars; ++j)
            if (zID == varDesc[j].zaxisID)
              varDesc[j].zaxisID = CDI_UNDEFID;
        }
    }
  gridDestroy ( gridID );
  free(var);
#if USE_MPI
  free(blk_displ);
#endif
  free(varDesc);
  free(levs);
  free(lats);
  free(lons);
}

#ifdef USE_MPI
static void
findPartition2D(int npart[2], int num_parts)
{
  const uint64_t rscale = 256;
  uint32_t *factors = NULL;
  xassert(num_parts > 0);
  int numFactors
    = PPM_prime_factorization_32((uint32_t)num_parts, &factors);
  /* try to distribute prime factors on dimensions to get
   * approx. 2 times as many parts in x dim than y dim */
  const uint64_t optimumRatio = rscale * 2;
  npart[0] = num_parts, npart[1] = 1;
  uint_fast32_t npart_attempt[2];
  uint64_t bestRatio = (uint64_t)num_parts * rscale,
    bestDiff = (uint64_t)llabs((long long)(bestRatio - optimumRatio));
  /* test all assignments of factors to dimensions, starting with
   * only one assigned to x dim (omitting 0 because that would
   * always give npart[1] > npart[0] */
  for (int assign2X = 1; assign2X <= numFactors; ++assign2X)
    {
      uint_fast32_t pattern = (UINT32_C(1) << assign2X) - 1,
        lastPattern = pattern << (numFactors - assign2X);
      do {
        npart_attempt[0] = 1;
        npart_attempt[1] = 1;
        /* loop over all factors */
        for (uint_fast32_t i = 0; i < (uint_fast32_t)numFactors; ++i)
          {
            uint_fast32_t dim_idx = (pattern >> i) & 1;
            npart_attempt[dim_idx] *= factors[i];
          }
        uint64_t ratio = ((uint64_t)npart_attempt[0] * rscale)
          / (uint64_t)npart_attempt[1];
        uint64_t diff = (uint64_t)llabs((long long)(ratio - optimumRatio));
        if (diff < bestDiff)
          {
            npart[0] = (int)npart_attempt[0];
            npart[1] = (int)npart_attempt[1];
            bestDiff = diff;
            bestRatio = ratio;
          }
        {
          uint_fast32_t t;
#if HAVE_DECL___BUILTIN_CTZ
          t = pattern | (pattern - 1);
          pattern = (t + 1)
            | (((~t & -~t) - 1) >> (__builtin_ctz((unsigned)pattern) + 1));
#else
          t = (pattern | (pattern - 1)) + 1;
          pattern = t | ((((t & -t) / (pattern & -pattern)) >> 1) - 1);
#endif
        }
      } while (pattern <= lastPattern);
    }
  free(factors);
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
