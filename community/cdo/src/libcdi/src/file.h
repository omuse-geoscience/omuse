#ifndef _FILE_H
#define _FILE_H

#include <stdio.h>
#include <sys/types.h>


#define  FILE_UNDEFID      -1

#define  FILE_TYPE_OPEN     1
#define  FILE_TYPE_FOPEN    2

/* buffer types for FILE_TYPE_OPEN */
#define  FILE_BUFTYPE_STD   1
#define  FILE_BUFTYPE_MMAP  2

const
char  *fileLibraryVersion(void);

void   fileDebug(int debug);

void  *filePtr(int fileID);

int    fileSetBufferType(int fileID, int type);
void   fileSetBufferSize(int fileID, long buffersize);

int    fileOpen(const char *filename, const char *mode);
int    fileOpen_serial(const char *filename, const char *mode);
int    fileClose(int fileID);
int    fileClose_serial(int fileID);

char  *fileInqName(int fileID);
int    fileInqMode(int fileID);

int    fileFlush(int fileID);
void   fileClearerr(int fileID);
int    fileEOF(int fileID);
int    filePtrEOF(void *fileptr);
void   fileRewind(int fileID);

off_t  fileGetPos(int fileID);
int    fileSetPos(int fileID, off_t offset, int whence);

int    fileGetc(int fileID);
int    filePtrGetc(void *fileptr);

size_t filePtrRead(void *fileptr, void *restrict ptr, size_t size);
size_t fileRead(int fileID, void *restrict ptr, size_t size);
size_t fileWrite(int fileID, const void *restrict ptr, size_t size);

#endif  /* _FILE_H */
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
