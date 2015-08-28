#ifndef _STREAM_IEG_H
#define _STREAM_IEG_H

#ifndef _IEG_H
#  include "ieg.h"
#endif

int    iegInqContents(stream_t *streamptr);
int    iegInqTimestep(stream_t *streamptr, int tsID);

int    iegInqRecord(stream_t *streamptr, int *varID, int *levelID);
void   iegDefRecord(stream_t *streamptr);
void   iegCopyRecord(stream_t *streamptr2, stream_t *streamptr1);
void   iegReadRecord(stream_t *streamptr, double *data, int *nmiss);
void   iegWriteRecord(stream_t *streamptr, const double *data);

void   iegReadVarDP (stream_t *streamptr, int varID,       double *data, int *nmiss);
void   iegWriteVarDP(stream_t *streamptr, int varID, const double *data);

void   iegReadVarSliceDP (stream_t *streamptr, int varID, int levelID,       double *data, int *nmiss);
void   iegWriteVarSliceDP(stream_t *streamptr, int varID, int levelID, const double *data);

#endif  /* _STREAM_IEG_H */
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
