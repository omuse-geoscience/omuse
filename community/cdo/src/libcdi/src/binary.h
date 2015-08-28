#ifndef _BINARY_H
#define _BINARY_H

#include <stdio.h>

#ifndef _DTYPES_H
#include "dtypes.h"
#endif


UINT32 get_UINT32(unsigned char *x);
UINT32 get_SUINT32(unsigned char *x);
UINT64 get_UINT64(unsigned char *x);
UINT64 get_SUINT64(unsigned char *x);


size_t binReadF77Block(int fileID, int byteswap);
void   binWriteF77Block(int fileID, int byteswap, size_t blocksize);

int binReadInt32(int fileID, int byteswap, size_t size, INT32 *ptr);
int binReadInt64(int fileID, int byteswap, size_t size, INT64 *ptr);

int binWriteInt32(int fileID, int byteswap, size_t size, INT32 *ptr);
int binWriteInt64(int fileID, int byteswap, size_t size, INT64 *ptr);

int binReadFlt32(int fileID, int byteswap, size_t size, FLT32 *ptr);
int binReadFlt64(int fileID, int byteswap, size_t size, FLT64 *ptr);

int binWriteFlt32(int fileID, int byteswap, size_t size, FLT32 *ptr);
int binWriteFlt64(int fileID, int byteswap, size_t size, FLT64 *ptr);

#endif  /* _BINARY_H */
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
