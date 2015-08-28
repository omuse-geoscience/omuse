#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <stdio.h>

#include "cdi.h"
#include "cksum.h"
#include "stream_cksum.h"
#include "dmemory.h"

struct cksum_table *
cksum_stream(const char *fname, size_t *table_len)
{
  int taxisID, vlistID, streamID;
  int nvars, ngrids, nzaxis;
  int i;
  uint32_t *checksum_state = NULL;
  struct
  {
    int x, y, z;
    size_t chars;
  } *varSize = NULL;
  size_t var_size_max_chars = 0;
  double *buf = NULL;
  struct cksum_table *file_vars = NULL;

  do {
    // Open the dataset
    streamID = streamOpenRead(fname);
    if ( streamID < 0 )
      {
        fprintf(stderr, "Cannot open data input file %s: %s\n",
                fname, cdiStringError(streamID));
        nvars = -1;
        break;
      }

    // Get the variable list of the dataset
    vlistID = streamInqVlist(streamID);

    nvars = vlistNvars(vlistID);
    ngrids = vlistNgrids(vlistID);
    nzaxis = vlistNzaxis(vlistID);
    if (nzaxis < 0 || ngrids < 0)
      {
        fprintf(stderr, "Error in grid/zaxis count query %d:%d\n",
                ngrids, nzaxis);
        nvars = -1;
        break;
      }
    checksum_state = xcalloc((size_t)nvars, sizeof (checksum_state[0]));
    varSize = xmalloc((size_t)nvars * sizeof (varSize[0]));

    for (i = 0; i < nvars; ++i)
      {
        int grid = vlistInqVarGrid(vlistID, i), gridType;
        int zaxis = vlistInqVarZaxis(vlistID, i);
        if (grid == CDI_UNDEFID || zaxis == CDI_UNDEFID)
          {
            fputs("error in axis/grid inquiry\n", stderr);
            nvars = -1;
            break;
          }
      if ((varSize[i].z = zaxisInqSize(zaxis)) <= 0)
        {
          fputs("invalid Z-axis found\n", stderr);
          nvars = -1;
          break;
        }
      if ((gridType = gridInqType(grid)) != GRID_LONLAT
          && gridType != GRID_GENERIC)
        {
          fprintf(stderr, "unexpected non-lonlat grid found: %d\n",
                  gridType);
          nvars = -1;
          break;
        }
      if ((varSize[i].x = gridInqXsize(grid)) < 0)
        {
          fprintf(stderr, "invalid X-size found: %d\n", varSize[i].x);
          nvars = -1;
          break;
        }
      if (varSize[i].x == 0) varSize[i].x = 1;
      if ((varSize[i].y = gridInqYsize(grid)) < 0)
        {
          fprintf(stderr, "invalid Y-size found: %d\n", varSize[i].y);
          nvars = -1;
          break;
        }
      if (varSize[i].y == 0) varSize[i].y = 1;
      varSize[i].chars = (size_t)varSize[i].x * (size_t)varSize[i].y
        * (size_t)varSize[i].z * sizeof (buf[0]);
      if (var_size_max_chars < varSize[i].chars)
        var_size_max_chars = varSize[i].chars;
    }
    buf = xmalloc(var_size_max_chars);

    if (nvars == -1)
      break;

    // Get the Time axis from the variable list
    taxisID = vlistInqTaxis(vlistID);

    int tsID = 0;
    // Inquire the time step
    while (streamInqTimestep(streamID, tsID))
      {
        // Get the verification date and time
        int vdatetime[2] = { taxisInqVtime(taxisID), taxisInqVdate(taxisID) };
        // Read var1 and var2
        for (i = 0; i < nvars; ++i)
          {
            int nmiss;
            streamReadVar(streamID, i, buf, &nmiss);
            memcrc_r(checksum_state + i, (const unsigned char *)vdatetime, sizeof (vdatetime));
            memcrc_r(checksum_state + i, (const unsigned char *)buf,
                     varSize[i].chars);
          }
        ++tsID;
      }

    file_vars = xmalloc((size_t)nvars * sizeof (file_vars[0]));
    for (i = 0; i < nvars; ++i)
      {
        file_vars[i].code = vlistInqVarCode(vlistID, i);
        file_vars[i].cksum
          = memcrc_finish(checksum_state + i,
                          (off_t)((varSize[i].chars + sizeof (int) * 2) * (size_t)tsID));
      }
    // Close the input stream
    streamClose(streamID);

  } while (0);

  // free resources
  free(checksum_state);
  free(varSize);
  free(buf);
  *table_len = (size_t)nvars;

  return file_vars;
}

