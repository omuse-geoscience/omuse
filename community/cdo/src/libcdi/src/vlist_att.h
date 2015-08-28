#ifndef VLIST_ATT_H
#define VLIST_ATT_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

int
vlistAttsGetSize(vlist_t *p, int varID, void *context);

void
vlistAttsPack(vlist_t *p, int varID,
              void * buf, int size, int *position, void *context);

void
vlistAttsUnpack(int vlistID, int varID,
                void * buf, int size, int *position, void *context);


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
