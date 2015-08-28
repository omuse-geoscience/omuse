#ifndef PIO_INTERFACE_
#define PIO_INTERFACE_

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <mpi.h>
#include <yaxt.h>

#include "resource_handle.h"
#include "pio_rpc.h"

void
pioBufferPartData(int streamID, int varID, const double *data,
                  int nmiss, Xt_idxlist partDesc);
void
cdiPioBufferPartDataGather(int streamID, int varID, const double *data,
                           int numBlocks, const int blocklengths[],
                           const int displacements[],
                           int nmiss, Xt_idxlist partDesc);

void pioBufferFuncCall(struct winHeaderEntry header,
                       const void *data, valPackFunc dataPackFunc);


struct memCpyDataDesc
{
  const void *obj;
  size_t obj_size;
};

void memcpyPackFunc(void *dataDesc, void *buf, int size, int *pos, void *context);

extern float cdiPIOpartInflate_;

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
