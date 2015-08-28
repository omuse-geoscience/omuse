#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <errno.h>
#include <inttypes.h>
#include <limits.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#ifdef USE_MPI
#include <mpi.h>
#include <yaxt.h>
#else
typedef int MPI_Comm;
#endif

#include "cdi.h"
#include "dmemory.h"
#include "pio_write.h"
#ifdef USE_MPI
#include "cdipio.h"
#include "pio_util.h"
#endif

void
modelRun(struct model_config setup, MPI_Comm comm);

static const struct model_config default_setup
  = { .nlon = 12, .nts = 3, .nlat = 6, .nvars = 5,
      .filetype = FILETYPE_GRB, .datatype = DATATYPE_PACK24,
      .compute_checksum = 1,
      .suffix = "grb",
      .max_nlev = 5,
};

#ifdef USE_MPI
static const struct {
  char *text;
  int mode;
} mode_map[] = {
  { "PIO_MPI", PIO_MPI },
  { "PIO_FPGUARD", PIO_FPGUARD },
  { "PIO_ASYNCH", PIO_ASYNCH },
  { "PIO_WRITER", PIO_WRITER }
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

int main(int argc, char *argv[])
{
  struct model_config setup = default_setup;

  MPI_Comm commModel;
#ifdef USE_MPI
  MPI_Comm commGlob;
  int sizeGlob;
  int rankGlob;
  int IOMode = PIO_MPI;
  int nProcsIO = 2;

  xmpi ( MPI_Init ( &argc, &argv));
  commGlob = MPI_COMM_WORLD;
  xt_initialize(commGlob);
  xmpi ( MPI_Comm_set_errhandler ( commGlob, MPI_ERRORS_RETURN ));
  xmpi ( MPI_Comm_size ( commGlob, &sizeGlob ));
  xmpi ( MPI_Comm_rank ( commGlob, &rankGlob ));
#endif

  {
    int opt;
    while ((opt = getopt(argc, argv, "f:m:n:z:t:y:c"
#ifdef USE_MPI
                         "p:w:"
#endif
                         )) != -1)
      switch (opt) {
#ifdef USE_MPI
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
        }
        break;
#endif
      case 'f':
        {
          int found = 0;
          for (size_t i = 0;
               i < sizeof (suffix2type) / sizeof (suffix2type[0]);
               ++i)
            if (!strcmp(optarg, suffix2type[i].suffix))
              {
                found = 1;
                setup.filetype = suffix2type[i].type;
                setup.suffix = suffix2type[i].suffix;
                setup.datatype = suffix2type[i].defaultDT;
                break;
              }
          if (!found)
            {
              fprintf(stderr, "Unsupported format requested: %s\n", optarg);
              exit(EXIT_FAILURE);
            }
        }
        break;
      case 'm':
        setup.nlon = parse_intarg("error parsing number of longitudes");
        break;
      case 'n':
        setup.nlat = parse_intarg("error parsing number of latitudes");
        break;
      case 'y':
        setup.nvars = parse_intarg("error parsing number of variables");
        if (setup.nvars < 1)
          {
            fputs("number of levels must be greater than zero!\n",
                  stderr);
            exit(EXIT_FAILURE);
          }
        if (setup.nvars > 127)
          {
            fputs("number of variables must not exceed 127!\n",
                  stderr);
            exit(EXIT_FAILURE);
          }
        break;
      case 'z':
        setup.max_nlev = parse_intarg("error parsing number of levels");
        if (setup.max_nlev < 1)
          {
            fputs("number of levels must be greater than zero!\n",
                  stderr);
            exit(EXIT_FAILURE);
          }
        break;
      case 't':
        setup.nts = parse_intarg("error parsing number of timesteps");
        break;
      case 'c':
        setup.compute_checksum = 0;
        break;
      default: /* '?' */
        fprintf(stderr, "Usage: %s "
                "[-m nlon] [-n nlat] [-z nlev] [-t nts] [-y num_vars]"
#ifdef USE_MPI
                " [-p PIO_MODE] [-w NIOSERVERS] [-c]"
#endif
                "\n", argv[0]);
        exit(EXIT_FAILURE);
      }

  }

#ifdef USE_MPI
  int pioNamespace;
  commModel = pioInit(commGlob, nProcsIO, IOMode, &pioNamespace, 1.0,
                      cdiPioNoPostCommSetup);
  if (commModel != MPI_COMM_NULL)
    {
      namespaceSetActive(pioNamespace);
#else
      commModel = -1;
#endif

      modelRun (setup, commModel);

#ifdef USE_MPI
    }
  pioFinalize ();
  xt_finalize();
  MPI_Finalize ();
#endif
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
