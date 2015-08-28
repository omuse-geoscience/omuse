#ifndef INCLUDE_GUARD_CDI_GRIB_FILE_H
#define INCLUDE_GUARD_CDI_GRIB_FILE_H

#include "referenceCounting.h"

/*
CdiInputFile is a file abstraction that allows accessing an input file through any number of channels:
It is reference counted, so that it is closed at the right place,
and it is stateless, so that accesses from different callers cannot interfere with each other.
Once the reference counting code is threadsafe, CdiInputFile will also be threadsafe.
*/
typedef struct CdiInputFile {
  //public:
    CdiReferencedObject super;

  //private:
    char* path;
    int fileDescriptor;
} CdiInputFile;

//Final class, the constructor is private and not defined here.
CdiInputFile* cdiInputFile_make(const char* path);   //The caller is responsible to call cdiRefObject_release() on the returned object.
int cdiInputFile_read(const CdiInputFile* me, off_t readPosition, size_t readSize, size_t* outActualReadSize, void* buffer);       //Returns one of CDI_EINVAL, CDI_ESYSTEM, CDI_EEOF, OR CDI_NOERR.
/* Returns path string, don't use after destruction of CdiInputFile
 * object */
const char* cdiInputFile_getPath(const CdiInputFile* me);
//Destructor is private as well.

#endif
