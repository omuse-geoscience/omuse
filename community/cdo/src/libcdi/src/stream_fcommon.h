#ifndef STREAM_FCOMMON_H
#define STREAM_FCOMMON_H

#ifndef  _CDI_INT_H
#include "cdi_int.h"
#endif

enum {
  SINGLE_PRECISION = 4,
  DOUBLE_PRECISION = 8,
};

void streamFCopyRecord(stream_t *streamptr2, stream_t *streamptr1,
                       const char *container_name);

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
