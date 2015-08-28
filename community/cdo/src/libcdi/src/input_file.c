#define _XOPEN_SOURCE 600
#include "input_file.h"

#include "cdi.h"
#include "dmemory.h"

#include <errno.h>
#include <fcntl.h>
#include <pthread.h>
#include <string.h>
#include <unistd.h>

static void cdiInputFile_destruct(CdiInputFile* me);

//For an explanation of the condestruct() pattern, see the comment in iterator_grib.c
//path != NULL -> construction
//path = NULL -> destruction
static CdiInputFile* cdiInputFile_condestruct(CdiInputFile* me, const char* path)
{
  #define super() (&me->super)
  if(!path) goto destruct;
  cdiRefObject_construct(super());
  me->path = strdup(path);
  if(!me->path) goto destructSuper;
  do
    {
      me->fileDescriptor = open(me->path, O_RDONLY);
    }
  while(me->fileDescriptor == -1 && (errno == EINTR || errno == EAGAIN));
  if(me->fileDescriptor == -1) goto freePath;
  //construction successfull, now we can set our own destructor
  super()->destructor = (void(*)(CdiReferencedObject*))cdiInputFile_destruct;
  goto success;

// ^        constructor code       ^
// |                               |
// v destructor/error-cleanup code v

destruct:
  close(me->fileDescriptor);
freePath:
  free(me->path);
destructSuper:
  cdiRefObject_destruct(super());
  me = NULL;

success:
  return me;
  #undef super
}

static CdiInputFile** openFileList = NULL;
static size_t openFileCount = 0, openFileListSize = 0;
static pthread_mutex_t openFileListLock = PTHREAD_MUTEX_INITIALIZER;

//This either returns a new object, or retains and returns a preexisting open file.
CdiInputFile* cdiInputFile_make(const char* path)
{
  CdiInputFile* result = NULL;
  xassert(path);
  int error = pthread_mutex_lock(&openFileListLock);
  xassert(!error);
    {
      //Check the list of open files for the given path.
      for(size_t i = openFileCount; i-- && !result; )
        {
          if(!strcmp(path, openFileList[i]->path)) result = openFileList[i];
        }
      //If no open file was found, we open one, otherwise we just retain the existing one one more time.
      if(result)
        {
          cdiRefObject_retain(&result->super);
        }
      else
        {
          result = xmalloc(sizeof(*result));
          if(!cdiInputFile_condestruct(result, path))
            {
              //An error occured during construction, avoid a memory leak.
              free(result);
              result = NULL;
            }
          else
            {
              //Add the new file to the list of open files.
              if(openFileCount == openFileListSize)
                {
                  openFileListSize *= 2;
                  if(openFileListSize < 16) openFileListSize = 16;
                  openFileList = xrealloc(openFileList, openFileListSize);
                }
              xassert(openFileCount < openFileListSize);
              openFileList[openFileCount++] = result;
            }
        }
    }
  error = pthread_mutex_unlock(&openFileListLock);
  xassert(!error);
  return result;
}

int cdiInputFile_read(const CdiInputFile* me, off_t readPosition, size_t readSize, size_t* outActualReadSize, void* buffer)
{
  char* byteBuffer = buffer;
  size_t trash;
  if(!outActualReadSize) outActualReadSize = &trash;
  *outActualReadSize = 0;
  while(readSize)
    {
      ssize_t bytesRead = pread(me->fileDescriptor, byteBuffer, readSize, readPosition);
      if(bytesRead == -1) return (errno == EINVAL) ?  CDI_EINVAL : CDI_ESYSTEM;
      if(bytesRead == 0) return CDI_EEOF;
      byteBuffer += bytesRead;
      readPosition += bytesRead;
      readSize -= (size_t)bytesRead;
      *outActualReadSize += (size_t)bytesRead;
    }
  return CDI_NOERR;
}

const char* cdiInputFile_getPath(const CdiInputFile* me)
{
  return me->path;
}

void cdiInputFile_destruct(CdiInputFile* me)
{
  int error = pthread_mutex_lock(&openFileListLock);
  xassert(!error);
    {
      //Find the position of me in the list of open files.
      ssize_t position = (ssize_t)openFileCount;
      while (position > 0 && openFileList[--position] != me);
      //Remove me from the list
      openFileList[position] = openFileList[--openFileCount];
    }
  error = pthread_mutex_unlock(&openFileListLock);
  xassert(!error);
  cdiInputFile_condestruct(me, NULL);
}
