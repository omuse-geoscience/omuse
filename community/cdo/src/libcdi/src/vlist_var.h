#ifndef VLIST_VAR_H
#define VLIST_VAR_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef _VLIST_H
#include "vlist.h"
#endif

int  vlistVarGetPackSize(vlist_t *p, int varID, void *context);
void vlistVarPack(vlist_t *p, int varID,
                  char * buffer, int bufferSize, int * pos, void *context);
void vlistVarUnpack(int vlistID,
                    char * buf, int size, int *position, int, void *context);
int vlistVarCompare(vlist_t *a, int varIDA, vlist_t *b, int varIDB);
void vlistDefVarIOrank    ( int, int, int );
int  vlistInqVarIOrank    ( int, int );

void cdiVlistCreateVarLevInfo(vlist_t *vlistptr, int varID);

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
