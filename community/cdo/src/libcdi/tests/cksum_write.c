#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <errno.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "cdi.h"

#include "cksum.h"
#include "simple_model_helper.h"

static int
parse_intarg(const char msg[])
{
  char *end;
  long temp = strtol(optarg, &end, 0);
  if ((errno == ERANGE && (temp == LONG_MAX || temp == LONG_MIN))
      || (errno != 0 && temp == 0)) {
    perror(msg);
    exit(EXIT_FAILURE);
  }
  if (temp > INT_MAX || temp < INT_MIN)
  {
    fprintf(stderr, "range error: %ld\n", temp);
    exit(EXIT_FAILURE);
  }
  return (int)temp;
}

static void
allocError(const char *msg)
{
  perror(msg);
  exit(EXIT_FAILURE);
}

static char *
fname_create(const char *prefix, const char *suffix)
{
  size_t prefix_len, suffix_len;
  char *fname;
  if (!(fname = malloc((prefix_len = strlen(prefix)) + 1
                       + (suffix_len = strlen(suffix)) + 1)))
    allocError("cannot create string");
  strcpy(fname, prefix);
  fname[prefix_len] = '.';
  strcpy(fname + prefix_len + 1, suffix);
  return fname;
}

static inline void
check_positive(int v, const char *msg)
{
  if (v < 1)
    {
      fprintf(stderr, "error: number of %s must be positive!\n", msg);
      exit(EXIT_FAILURE);
    }
}

#ifdef TEST_CHUNK_WRITE
static void
get_chunk(double *chunkBuf, double *var, int varShape[3],
          int chunk[3][2])
{
  size_t ofs = 0;
  unsigned start_k = (unsigned)chunk[2][0], start_j = (unsigned)chunk[1][0],
    start_i = (unsigned)chunk[0][0];
  unsigned size_k = (unsigned)chunk[2][1] - (unsigned)chunk[2][0] + 1,
    size_j = (unsigned)chunk[1][1] - (unsigned)chunk[1][0] + 1,
    size_i = (unsigned)chunk[0][1] - (unsigned)chunk[0][0] + 1;
  size_t stride_k = (size_t)varShape[0] * (size_t)varShape[1],
    stride_j = (size_t)varShape[0];
  for (unsigned k = 0; k < size_k ; ++k)
    for (unsigned j = 0; j < size_j; ++j)
      for (unsigned i = 0; i < size_i; ++i)
        chunkBuf[ofs++] = var[(k + start_k) * stride_k
                              + (j + start_j) * stride_j + (i + start_i)];
}
#endif

#ifndef TEST_CHUNK_WRITE
static const struct {
  char suffix[4];
  int type, defaultDT, defaultGrid;
} suffix2type[] = {
  { "nc", FILETYPE_NC, DATATYPE_FLT64, GRID_LONLAT },
  { "grb",  FILETYPE_GRB, DATATYPE_PACK24, GRID_LONLAT },
  { "nc2", FILETYPE_NC2, DATATYPE_FLT64, GRID_LONLAT },
  { "nc4", FILETYPE_NC4, DATATYPE_FLT64, GRID_LONLAT },
  { "ext", FILETYPE_EXT, DATATYPE_FLT64, GRID_GENERIC, },
  { "svc", FILETYPE_SRV, DATATYPE_FLT64, GRID_GENERIC, },
  { "ieg", FILETYPE_IEG, DATATYPE_FLT64, GRID_LONLAT },
};
#endif
enum {
  nvars = 2,
};

static const int varCodes[nvars] = { 42, 55 };

int
main(int argc, char *argv[])
{
  int gridID, zaxisID[nvars], taxisID;
  int vlistID, varID[nvars], streamID, tsID;
  int nlon = 12, //!< Number of longitudes
    nlat = 6, //!< Number of latitudes
    nlev = 5, //!< Number of levels
    nts = 3; //!< Number of time steps
  int i, j, k, nmiss = 0;
  double *lons, *lats, *var[nvars], *levs, mscale, mrscale;
  size_t varSize[nvars];
  char *varName[nvars] = { "varname1", "varname2" };
#ifndef TEST_CHUNK_WRITE
  const char *suffix = "grb", *prefix = "example";
  int grid = GRID_LONLAT;
  int filetype = FILETYPE_GRB, datatype = DATATYPE_PACK24;
#else
  const char *suffix = "nc", *prefix = "example";
  int grid = GRID_LONLAT;
  int filetype = FILETYPE_NC, datatype = DATATYPE_FLT64;
#endif
  {
    int opt;
    while ((opt = getopt(argc, argv,
#ifndef TEST_CHUNK_WRITE
                         "f:"
#endif
                         "m:n:o:t:")) != -1)
      switch (opt) {
#ifndef TEST_CHUNK_WRITE
      case 'f':
        {
          int found = 0;
          for (size_t i = 0;
               i < sizeof (suffix2type) / sizeof (suffix2type[0]);
               ++i)
            if (!strcmp(optarg, suffix2type[i].suffix))
              {
                found = 1;
                filetype = suffix2type[i].type;
                suffix = suffix2type[i].suffix;
                datatype = suffix2type[i].defaultDT;
                break;
              }
          if (!found)
            {
              fprintf(stderr, "Unsupported format requested: %s\n", optarg);
              exit(EXIT_FAILURE);
            }
        }
        break;
#endif
      case 'm':
        nlon = parse_intarg("error parsing number of longitudes");
        check_positive(nlon, "longitudes");
#ifdef TEST_CHUNK_WRITE
        if (nlon < 2)
          {
            fputs("number of longitudes must be larger 1 for "
                  "chunk write test\n", stderr);
            exit(EXIT_FAILURE);
          }
#endif
        break;
      case 'n':
        check_positive(nlat = parse_intarg("error parsing number of latitudes"),
                       "latitudes");
        break;
      case 'o':
        check_positive(nlev = parse_intarg("error parsing number of levels"),
                       "levels");
        break;
      case 't':
        check_positive(nts = parse_intarg("error parsing number of timesteps"),
                       "timesteps");
        break;
      default: /* '?' */
        fprintf(stderr, "Usage: %s [-m nlon] [-n nlat] [-o nlev] [-t nts]\n", argv[0]);
        exit(EXIT_FAILURE);
      }
  }

  lons = malloc((size_t)nlon * sizeof (lons[0]));
  for (i = 0; i < nlon; ++i)
    lons[i] = ((double)(i * 360))/nlon;
  lats = malloc((size_t)nlat * sizeof (lats[0]));
  for (i = 0; i < nlat; ++i)
    lats[i] = ((double)(i * 180))/nlat - 90.0;
  levs = malloc((size_t)nlev * sizeof (levs[0]));
  for (i = 0; i < nlev; ++i)
    levs[i] = 101300 - floor(3940.3 * (exp(2.3579 * (double)(i)/(nlev - 1)) - 1.0));

  varSize[0] = (size_t)nlon * (size_t)nlat;
  varSize[1] = (size_t)nlon * (size_t)nlat * (size_t)nlev;

  // Create a regular lon/lat grid
  gridID = gridCreate(grid, nlon*nlat);
  gridDefXsize(gridID, nlon);
  gridDefYsize(gridID, nlat);
  gridDefXvals(gridID, lons);
  gridDefYvals(gridID, lats);

  // Create a surface level Z-axis
  zaxisID[0] = zaxisCreate(ZAXIS_SURFACE, 1);

  // Create a pressure level Z-axis
  zaxisID[1] = zaxisCreate(ZAXIS_PRESSURE, nlev);
  zaxisDefLevels(zaxisID[1], levs);

  /* add uuids to zaxis and grid */
  {
    unsigned char uuid[16];
    int str2uuid(const char *uuidstr, unsigned char *uuid);

    static char gridUUIDTxt[] = "107d7a5b-348c-4d1a-90a9-d745914f2fb6";

    str2uuid(gridUUIDTxt, uuid);
    gridDefUUID(gridID, uuid);

    static char zaxisUUIDTxt[2][37] = {
      { "d157f399-5496-4097-a3d8-437a6dda6311" },
      { "6f784a65-bce8-48c9-afa4-4c40130709c7" }
    };

    for (int i = 0; i < 2; ++i)
      {
        str2uuid(zaxisUUIDTxt[i], uuid);
        zaxisDefUUID(zaxisID[i], uuid);
      }
  }

  // Create a Time axis
  taxisID = taxisCreate(TAXIS_ABSOLUTE);

  // Create a variable list
  vlistID = vlistCreate();

  for (i = 0; i < nvars; ++i)
    {
      // Define the variables
      varID[i] = vlistDefVar(vlistID, gridID, zaxisID[i], TIME_VARIABLE);
      // Define the variable names,
      vlistDefVarName(vlistID, varID[i], varName[i]);
      // the codes
      vlistDefVarCode(vlistID, varID[i], varCodes[i]);
      // and set the data type
      vlistDefVarDatatype(vlistID, varID[i], datatype);
      // create memory for variables
      var[i] = malloc(varSize[i] * sizeof (var[i][0]));
    }

  var_scale(datatype, &mscale, &mrscale);

  // Assign the Time axis to the variable list
  vlistDefTaxis(vlistID, taxisID);

  // Create a dataset
  {
    char *fname = fname_create(prefix, suffix);
    if ((streamID = streamOpenWrite(fname, filetype)) < 0)
      {
        fprintf(stderr, "error opening output file %s: %s\n",
                fname, cdiStringError(streamID));
        exit(EXIT_FAILURE);
      }
    free(fname);
  }

  // Assign the variable list to the dataset
  streamDefVlist(streamID, vlistID);

  {
    uint32_t checksum_state[nvars] = { 0, 0 };
    // Loop over the number of time steps
    for ( tsID = 0; tsID < nts; tsID++ )
      {
        int vdatetime[2] = { 120000, 19850101+tsID };
        // Set the verification date to 1985-01-01 + tsID
        taxisDefVdate(taxisID, vdatetime[1]);
        // Set the verification time to 12:00:00
        taxisDefVtime(taxisID, vdatetime[0]);
        // Define the time step
        streamDefTimestep(streamID, tsID);

        // Init var1 and var2
        for (j = 0; j < nlat; j++)
          for (i = 0; i < nlon; i++)
            var[0][i+j*nlon]
              = sign_flat(round(
                   (sin(2.0 * M_PI * (lons[(i + tsID)%nlon] - lons[0])
                        / (lons[nlon-1] - lons[0]))
                    * cos(2.0 * M_PI * (lats[j] - lats[0])
                          / (lons[nlat-1] - lats[0]))
                    ) * mscale)) * mrscale;
        for (k = 0; k < nlev; ++k)
          for (j = 0; j < nlat; j++)
            for (i = 0; i < nlon; i++)
              var[1][i+j*nlon+k*nlon*nlat]
                = sign_flat(round(
                     (cos(2.0 * M_PI * (lons[(i + tsID)%nlon] - lons[0])
                          / (lons[nlon-1] - lons[0]))
                      * sin(2.0 * M_PI * (lats[j] - lats[0])
                            / (lons[nlat-1] - lats[0]))
                      ) * mscale)) * mrscale;

        if (filetype == FILETYPE_EXT)
          {
            /* EXTRA doesn't store time, only date
             * set the value to 0 before checksumming, because a
             * time field of 0 is what reading an EXTRA file will
             * return */
            vdatetime[0] = 0;
          }
        memcrc_r(&checksum_state[0], (const unsigned char *)vdatetime, sizeof (vdatetime));
        memcrc_r(&checksum_state[0], (const unsigned char *)var[0], varSize[0] * sizeof (var[0][0]));
        memcrc_r(&checksum_state[1], (const unsigned char *)vdatetime, sizeof (vdatetime));
        memcrc_r(&checksum_state[1], (const unsigned char *)var[1], varSize[1] * sizeof (var[1][0]));

        // Write var1 and var2
#ifdef TEST_CHUNK_WRITE
        {
          size_t maxChunkSize
            = ((size_t)nlon + 1)/2 * (size_t)nlat * (size_t)nlev;
          double *chunkBuf = malloc(maxChunkSize * sizeof (double));
          int varShape[2][3] = { { nlon, nlat, 1 }, { nlon, nlat, nlev } },
            chunk[3][2] = { { 0, nlon/2 - 1 }, { 0, nlat - 1 }, { 0, 0 } };
          if (!chunkBuf)
              allocError("cannot create chunk buffer");
          chunk[0][0] = 0; chunk[0][1] = nlon/2 - 1;
          chunk[2][1] = 0;
          get_chunk(chunkBuf, var[0], varShape[0], chunk);
          streamWriteVarChunk(streamID, varID[0], (const int (*)[2])chunk,
                              chunkBuf, nmiss);
          chunk[2][1] = nlev - 1;
          get_chunk(chunkBuf, var[1], varShape[1], chunk);
          streamWriteVarChunk(streamID, varID[1], (const int (*)[2])chunk,
                              chunkBuf, nmiss);
          chunk[0][0] = chunk[0][1] + 1; chunk[0][1] = nlon - 1;
          chunk[2][1] = 0;
          get_chunk(chunkBuf, var[0], varShape[0], chunk);
          streamWriteVarChunk(streamID, varID[0], (const int (*)[2])chunk,
                              chunkBuf, nmiss);
          chunk[2][1] = nlev - 1;
          get_chunk(chunkBuf, var[1], varShape[1], chunk);
          streamWriteVarChunk(streamID, varID[1], (const int (*)[2])chunk,
                              chunkBuf, nmiss);
          free(chunkBuf);
        }
#else
        streamWriteVar(streamID, varID[0], var[0], nmiss);
        streamWriteVar(streamID, varID[1], var[1], nmiss);
#endif
      }
    // write checksums to table file
    {
      FILE *tablefp;
      {
        char *fname = fname_create(prefix, "cksum");
        if (!(tablefp = fopen(fname, "w")))
          {
            perror("failed to open table file");
            exit(EXIT_FAILURE);
          }
        free(fname);
      }
      for (i = 0; i < nvars; ++i)
        {
          uint32_t cksum;
          int code;
          cksum
            = memcrc_finish(&checksum_state[i],
                            (off_t)((varSize[i] * sizeof (var[i][0])
                                     + sizeof (int) * 2) * (size_t)nts));
          code = vlistInqVarCode(vlistID, varID[i]);
          if (fprintf(tablefp, "%08lx %d\n", (unsigned long)cksum, code) < 0)
            {
              perror("failed to write table file");
              exit(EXIT_FAILURE);
            }
        }
      fclose(tablefp);
    }
  }

  // Close the output stream
  streamClose(streamID);

  // Destroy the objects
  for (i = 0; i < nvars; ++i)
    free(var[i]);
  vlistDestroy(vlistID);
  taxisDestroy(taxisID);
  zaxisDestroy(zaxisID[0]);
  zaxisDestroy(zaxisID[1]);
  gridDestroy(gridID);
  free(levs);
  free(lats);
  free(lons);
  return 0;
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
