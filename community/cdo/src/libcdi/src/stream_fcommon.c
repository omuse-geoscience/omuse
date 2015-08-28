#include <stdlib.h>

#include "cdi_int.h"
#include "dmemory.h"
#include "file.h"
#include "stream_fcommon.h"

void streamFCopyRecord(stream_t *streamptr2, stream_t *streamptr1,
                       const char *container_name)
{

  int fileID1 = streamptr1->fileID;
  int fileID2 = streamptr2->fileID;

  int tsID    = streamptr1->curTsID;
  int vrecID  = streamptr1->tsteps[tsID].curRecID;
  int recID   = streamptr1->tsteps[tsID].recIDs[vrecID];
  off_t recpos  = streamptr1->tsteps[tsID].records[recID].position;
  size_t recsize = streamptr1->tsteps[tsID].records[recID].size;

  if (fileSetPos(fileID1, recpos, SEEK_SET) != 0)
    Error("Cannot seek input file for %s record copy!", container_name);

  char *buffer = xmalloc(recsize);

  if (fileRead(fileID1, buffer, recsize) != recsize)
    Error("Failed to read record from %s file for copying!", container_name);

  if (fileWrite(fileID2, buffer, recsize) != recsize)
    Error("Failed to write record to %s file when copying!", container_name);

  free(buffer);
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
