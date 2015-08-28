#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <fcntl.h>
#include <limits.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>  // gettimeofday()

#include "dmemory.h"
#include "error.h"
#include "file.h"

#include "namespace.h"

#if ! defined (O_BINARY)
#define O_BINARY 0
#endif

#ifndef strdupx
#ifndef strdup
char *strdup(const char *s);
#endif
#define strdupx  strdup
/*
#define strdupx(s)                                \
({                                                \
   const char *__old = (s);                       \
   size_t __len = strlen(__old) + 1;              \
   char *__new = (char *) malloc(__len);          \
   (char *) memcpy(__new, __old, __len);          \
})
*/
#endif


#if defined (HAVE_MMAP)
#  include <sys/mman.h> /* mmap() is defined in this header */
#endif


#if ! defined   (FALSE)
#  define  FALSE  0
#endif

#if ! defined   (TRUE)
#  define  TRUE   1
#endif

/* #define  MAX_FILES  FOPEN_MAX */
#define  MAX_FILES  4096

static int _file_max = MAX_FILES;

static void file_initialize(void);

static int _file_init = FALSE;

#if  defined  (HAVE_LIBPTHREAD)
#include <pthread.h>

static pthread_once_t  _file_init_thread = PTHREAD_ONCE_INIT;
static pthread_mutex_t _file_mutex;

#  define FILE_LOCK()         pthread_mutex_lock(&_file_mutex)
#  define FILE_UNLOCK()       pthread_mutex_unlock(&_file_mutex)
#  define FILE_INIT()        \
   if ( _file_init == FALSE ) pthread_once(&_file_init_thread, file_initialize)

#else

#  define FILE_LOCK()
#  define FILE_UNLOCK()
#  define FILE_INIT()        \
   if ( _file_init == FALSE ) file_initialize()

#endif


typedef struct
{
  int        self;
  int        flag;           /* access and error flag         */
  int        eof;            /* end of file flag              */
  int        fd;             /* file descriptor used for read */
  FILE      *fp;             /* FILE pointer used for write   */
  int        mode;           /* file access mode              */
  char      *name;           /* file name                     */
  off_t      size;           /* file size                     */
  off_t      position;       /* file position                 */
  long       access;         /* file access                   */
  off_t      byteTrans;      /*                               */
  size_t     blockSize;      /* file block size               */
  int        type;           /* file type ( 1:open 2:fopen )  */
  int        bufferType;     /* buffer type ( 1:std 2:mmap )  */
  size_t     bufferSize;     /* file buffer size              */
  size_t     mappedSize;     /* mmap buffer size              */
  char      *buffer;         /* file buffer                   */
  long       bufferNumFill;  /* number of buffer fill         */
  char      *bufferPtr;      /* file buffer pointer           */
  off_t      bufferPos;
  off_t      bufferStart;
  off_t      bufferEnd;
  size_t     bufferCnt;
  double     time_in_sec;
}
bfile_t;


enum F_I_L_E_Flags
  {
    FILE_READ  =  01,
    FILE_WRITE =  02,
    FILE_UNBUF =  04,
    FILE_EOF   = 010,
    FILE_ERROR = 020
  };


static int FileInfo  = FALSE;


#if ! defined (MIN_BUF_SIZE)
#  define  MIN_BUF_SIZE  131072L
#endif


static size_t FileBufferSizeMin = MIN_BUF_SIZE;
static long   FileBufferSizeEnv = -1;
static int    FileBufferTypeEnv =  0;

static int    FileTypeRead  = FILE_TYPE_OPEN;
static int    FileTypeWrite = FILE_TYPE_FOPEN;
static int    FileFlagWrite = 0;

static int    FILE_Debug = 0;   /* If set to 1, debugging */


static void file_table_print(void);

/*
 * A version string.
 */
#undef   LIBVERSION
#define  LIBVERSION      1.8.2
#define  XSTRING(x)	 #x
#define  STRING(x) 	 XSTRING(x)
const char file_libvers[] = STRING(LIBVERSION) " of "__DATE__" "__TIME__;

/*
  21/05/2004  1.3.2 set min I/O Buffersize to 128k
  31/05/2005  1.4.0 replace fileTable by _fileList
  26/08/2005  1.4.1 fileClose with return value
                    checks for all fileptr
  01/09/2005  1.5.0 thread safe version
  06/11/2005  1.5.1 add filePtrEOF, filePtr, filePtrGetc
  03/02/2006  1.5.2 ansi C: define getpagesize and strdupx
  27/12/2007  1.6.0 add FILE_TYPE_FOPEN
  24/03/2008  1.6.1 add O_BINARY if available
                    remove default HAVE_MMAP
                    use HAVE_STRUCT_STAT_ST_BLKSIZE
  22/08/2010  1.7.0 refactor
  11/11/2010  1.7.1 update for changed interface of error.h
  02/02/2012  1.8.0 cleanup
  16/11/2012  1.8.1 added support for unbuffered write
  27/06/2013  1.8.2 added env. var. FILE_TYPE_WRITE (1:open; 2:fopen)
 */


typedef struct _filePtrToIdx {
  int idx;
  bfile_t *ptr;
  struct _filePtrToIdx *next;
} filePtrToIdx;


static filePtrToIdx *_fileList  = NULL;
static filePtrToIdx *_fileAvail = NULL;

static
void file_list_new(void)
{
  assert(_fileList == NULL);

  _fileList = (filePtrToIdx *)xmalloc((size_t)_file_max * sizeof (filePtrToIdx));
}

static
void file_list_delete(void)
{
  if ( _fileList )
    {
      free(_fileList);
      _fileList = NULL;
    }
}

static
void file_init_pointer(void)
{
  int  i;

  for ( i = 0; i < _file_max; i++ )
    {
      _fileList[i].next = _fileList + i + 1;
      _fileList[i].idx  = i;
      _fileList[i].ptr  = 0;
    }

  _fileList[_file_max-1].next = 0;

  _fileAvail = _fileList;
}

static
bfile_t *file_to_pointer(int idx)
{
  bfile_t *fileptr = NULL;

  FILE_INIT();

  if ( idx >= 0 && idx < _file_max )
    {
      FILE_LOCK();

      fileptr = _fileList[idx].ptr;

      FILE_UNLOCK();
    }
  else
    Error("file index %d undefined!", idx);

  return (fileptr);
}

/* Create an index from a pointer */
static
int file_from_pointer(bfile_t *ptr)
{
  int      idx = -1;
  filePtrToIdx *newptr;

  if ( ptr )
    {
      FILE_LOCK();

      if ( _fileAvail )
	{
	  newptr       = _fileAvail;
	  _fileAvail   = _fileAvail->next;
	  newptr->next = 0;
	  idx	       = newptr->idx;
	  newptr->ptr  = ptr;

	  if ( FILE_Debug )
	    Message("Pointer %p has idx %d from file list", ptr, idx);
	}
      else
	Warning("Too many open files (limit is %d)!", _file_max);

      FILE_UNLOCK();
    }
  else
    Error("Internal problem (pointer %p undefined)", ptr);

  return (idx);
}

static
void file_init_entry(bfile_t *fileptr)
{
  fileptr->self          = file_from_pointer(fileptr);

  fileptr->flag          = 0;
  fileptr->fd            = -1;
  fileptr->fp            = NULL;
  fileptr->mode          = 0;
  fileptr->size          = 0;
  fileptr->name          = NULL;
  fileptr->access        = 0;
  fileptr->position      = 0;
  fileptr->byteTrans     = 0;
  fileptr->type          = 0;
  fileptr->bufferType    = 0;
  fileptr->bufferSize    = 0;
  fileptr->mappedSize    = 0;
  fileptr->buffer        = NULL;
  fileptr->bufferNumFill = 0;
  fileptr->bufferStart   = 0;
  fileptr->bufferEnd     = -1;
  fileptr->bufferPos     = 0;
  fileptr->bufferCnt     = 0;
  fileptr->bufferPtr     = NULL;
  fileptr->time_in_sec   = 0.0;
}

static
bfile_t *file_new_entry(void)
{
  bfile_t *fileptr;

  fileptr = (bfile_t *) malloc(sizeof(bfile_t));

  if ( fileptr ) file_init_entry(fileptr);

  return (fileptr);
}

static
void file_delete_entry(bfile_t *fileptr)
{
  int idx;

  idx = fileptr->self;

  FILE_LOCK();

  free(fileptr);

  _fileList[idx].next = _fileAvail;
  _fileList[idx].ptr  = 0;
  _fileAvail   	      = &_fileList[idx];

  FILE_UNLOCK();

  if ( FILE_Debug )
    Message("Removed idx %d from file list", idx);
}


const char *fileLibraryVersion(void)
{
  return (file_libvers);
}


#ifndef POSIXIO_DEFAULT_PAGESIZE
#define POSIXIO_DEFAULT_PAGESIZE 4096
#endif

static
int pagesize(void)
{
#if defined(_SC_PAGESIZE)
  return ((int) sysconf(_SC_PAGESIZE));
#else
  return ((int) POSIXIO_DEFAULT_PAGESIZE);
#endif
}

static
double file_time()
{
  double tseconds = 0.0;
  struct timeval mytime;
  gettimeofday(&mytime, NULL);
  tseconds = (double) mytime.tv_sec + (double) mytime.tv_usec*1.0e-6;
  return (tseconds);
}

void fileDebug(int debug)
{
  FILE_Debug = debug;

  if ( FILE_Debug )
    Message("Debug level %d", debug);
}


void *filePtr(int fileID)
{
  bfile_t *fileptr;

  fileptr = file_to_pointer(fileID);

  return (fileptr);
}

static
void file_pointer_info(const char *caller, int fileID)
{
  if ( FILE_Debug )
    {
      fprintf(stdout, "%-18s : ", caller);
      fprintf(stdout, "The fileID %d underlying pointer is not valid!", fileID);
      fprintf(stdout, "\n");
    }
}


int fileSetBufferType(int fileID, int type)
{
  int ret = 0;
  bfile_t *fileptr;

  fileptr = file_to_pointer(fileID);

  if ( fileptr )
    {
      switch (type)
	{
	case FILE_BUFTYPE_STD:
	case FILE_BUFTYPE_MMAP:
	  fileptr->bufferType = type;
	  break;
	default:
	  Error("File type %d not implemented!", type);
	}
    }

#if ! defined (HAVE_MMAP)
  if ( type == FILE_BUFTYPE_MMAP ) ret = 1;
#endif

  return (ret);
}


int fileGetBufferType(int fileID)
{
  bfile_t *fileptr;
  int bufferType = 0;

  fileptr = file_to_pointer(fileID);

  if ( fileptr ) bufferType = fileptr->bufferType;

  return (bufferType);
}


int fileFlush(int fileID)
{
  bfile_t *fileptr;
  int retval = 0;

  fileptr = file_to_pointer(fileID);

  if ( fileptr ) retval = fflush(fileptr->fp);

  return (retval);
}


void fileClearerr(int fileID)
{
  bfile_t *fileptr;

  fileptr = file_to_pointer(fileID);

  if ( fileptr )
    {
      if ( fileptr->mode != 'r' )
	clearerr(fileptr->fp);
    }
}


int filePtrEOF(void *vfileptr)
{
  bfile_t *fileptr = (bfile_t *) vfileptr;
  int retval = 0;

  if ( fileptr ) retval = (fileptr->flag & FILE_EOF) != 0;

  return (retval);
}


int fileEOF(int fileID)
{
  bfile_t *fileptr;
  int retval = 0;

  fileptr = file_to_pointer(fileID);

  if ( fileptr ) retval = (fileptr->flag & FILE_EOF) != 0;

  return (retval);
}


int fileError(int fileID)
{
  bfile_t *fileptr;
  int retval = 0;

  fileptr = file_to_pointer(fileID);

  if ( fileptr ) retval = (fileptr->flag & FILE_ERROR) != 0;

  return (retval);
}


void fileRewind(int fileID)
{
  fileSetPos(fileID, (off_t) 0, SEEK_SET);
  fileClearerr(fileID);
}


off_t fileGetPos(int fileID)
{
  off_t filepos = 0;
  bfile_t *fileptr;

  fileptr = file_to_pointer(fileID);

  if ( fileptr )
    {
      if ( fileptr->mode == 'r' && fileptr->type == FILE_TYPE_OPEN )
	filepos = fileptr->position;
      else
	filepos = ftell(fileptr->fp);
    }

  if ( FILE_Debug ) Message("Position %ld", filepos);

  return (filepos);
}


int fileSetPos(int fileID, off_t offset, int whence)
{
  int status = 0;
  bfile_t *fileptr;

  fileptr = file_to_pointer(fileID);

  if ( FILE_Debug ) Message("Offset %8ld  Whence %3d", (long) offset, whence);

  if ( fileptr == 0 )
    {
      file_pointer_info(__func__, fileID);
      return (1);
    }

  switch (whence)
    {
    case SEEK_SET:
      if ( fileptr->mode == 'r' && fileptr->type == FILE_TYPE_OPEN )
	{
	  off_t position = offset;
	  fileptr->position = position;
	  if ( position < fileptr->bufferStart || position > fileptr->bufferEnd )
	    {
	      if ( fileptr->bufferType == FILE_BUFTYPE_STD )
		fileptr->bufferPos = position;
	      else
		fileptr->bufferPos = position - position % pagesize();

	      fileptr->bufferCnt = 0;
	      fileptr->bufferPtr = NULL;
	    }
	  else
	    {
	      if ( fileptr->bufferPos != fileptr->bufferEnd + 1 )
		{
		  if ( FILE_Debug )
		    Message("Reset buffer pos from %ld to %ld",
			    fileptr->bufferPos, fileptr->bufferEnd + 1);

		  fileptr->bufferPos = fileptr->bufferEnd + 1;
		}
	      fileptr->bufferCnt = (size_t)(fileptr->bufferEnd - position) + 1;
	      fileptr->bufferPtr = fileptr->buffer + position - fileptr->bufferStart;
	    }
	}
      else
	{
	  status = fseek(fileptr->fp, offset, whence);
	}
      break;
    case SEEK_CUR:
      if ( fileptr->mode == 'r' && fileptr->type == FILE_TYPE_OPEN )
	{
	  fileptr->position += offset;
	  off_t position = fileptr->position;
	  if ( position < fileptr->bufferStart || position > fileptr->bufferEnd )
	    {
	      if ( fileptr->bufferType == FILE_BUFTYPE_STD )
		fileptr->bufferPos = position;
	      else
		fileptr->bufferPos = position - position % pagesize();

	      fileptr->bufferCnt = 0;
	      fileptr->bufferPtr = NULL;
	    }
	  else
	    {
	      if ( fileptr->bufferPos != fileptr->bufferEnd + 1 )
		{
		  if ( FILE_Debug )
		    Message("Reset buffer pos from %ld to %ld",
			    fileptr->bufferPos, fileptr->bufferEnd + 1);

		  fileptr->bufferPos = fileptr->bufferEnd + 1;
		}
	      fileptr->bufferCnt -= (size_t)offset;
	      fileptr->bufferPtr += offset;
	    }
	}
      else
	{
	  status = fseek(fileptr->fp, offset, whence);
	}
      break;
    default:
      Error("Whence = %d not implemented", whence);
    }

  if ( fileptr->position < fileptr->size )
    if ( (fileptr->flag & FILE_EOF) != 0 )
      fileptr->flag -= FILE_EOF;

  return (status);
}

static
void file_table_print(void)
{
  int fileID;
  int lprintHeader = 1;
  bfile_t *fileptr;

  for ( fileID = 0; fileID < _file_max; fileID++ )
    {
      fileptr = file_to_pointer(fileID);

      if ( fileptr )
	{
	  if ( lprintHeader )
	    {
	      fprintf(stderr, "\nFile table:\n");
	      fprintf(stderr, "+-----+---------+");
	      fprintf(stderr, "----------------------------------------------------+\n");
	      fprintf(stderr, "|  ID |  Mode   |");
	      fprintf(stderr, "  Name                                              |\n");
	      fprintf(stderr, "+-----+---------+");
	      fprintf(stderr, "----------------------------------------------------+\n");
	      lprintHeader = 0;
	    }

	  fprintf(stderr, "| %3d | ", fileID);

	  switch ( fileptr->mode )
	    {
	    case 'r':
	      fprintf(stderr, "read   ");
	      break;
	    case 'w':
	      fprintf(stderr, "write  ");
	      break;
	    case 'a':
	      fprintf(stderr, "append ");
	      break;
	    default:
	      fprintf(stderr, "unknown");
	    }

          fprintf(stderr, " | %-51s|\n", fileptr->name);
	}
    }

  if ( lprintHeader == 0 )
    {
      fprintf(stderr, "+-----+---------+");
      fprintf(stderr, "----------------------------------------------------+\n");
    }
}


char *fileInqName(int fileID)
{
  bfile_t *fileptr;
  char *name = NULL;

  fileptr = file_to_pointer(fileID);

  if ( fileptr ) name = fileptr->name;

  return (name);
}


int fileInqMode(int fileID)
{
  bfile_t *fileptr;
  int mode = 0;

  fileptr = file_to_pointer(fileID);

  if ( fileptr ) mode = fileptr->mode;

  return (mode);
}

static
long file_getenv(const char *envName)
{
  char *envString;
  long envValue = -1;
  long fact = 1;

  envString = getenv(envName);

  if ( envString )
    {
      int loop;

      for ( loop = 0; loop < (int) strlen(envString); loop++ )
	{
	  if ( ! isdigit((int) envString[loop]) )
	    {
	      switch ( tolower((int) envString[loop]) )
		{
		case 'k':  fact =       1024;  break;
		case 'm':  fact =    1048576;  break;
		case 'g':  fact = 1073741824;  break;
		default:
		  fact = 0;
		  Message("Invalid number string in %s: %s", envName, envString);
		  Warning("%s must comprise only digits [0-9].",envName);
		}
	      break;
	    }
	}

      if ( fact ) envValue = fact*atol(envString);

      if ( FILE_Debug ) Message("Set %s to %ld", envName, envValue);
    }

  return (envValue);
}

static
void file_initialize(void)
{
  long value;
  char *envString;

#if  defined  (HAVE_LIBPTHREAD)
  /* initialize global API mutex lock */
  pthread_mutex_init(&_file_mutex, NULL);
#endif

  value = file_getenv("FILE_DEBUG");
  if ( value >= 0 ) FILE_Debug = (int) value;

  value = file_getenv("FILE_MAX");
  if ( value >= 0 ) _file_max = (int) value;

  if ( FILE_Debug )
    Message("FILE_MAX = %d", _file_max);

  FileInfo  = (int) file_getenv("FILE_INFO");

  value  = file_getenv("FILE_BUFSIZE");
  if ( value >= 0 ) FileBufferSizeEnv = value;
  else
    {
      value  = file_getenv("GRIB_API_IO_BUFFER_SIZE");
      if ( value >= 0 ) FileBufferSizeEnv = value;
    }

  value = file_getenv("FILE_TYPE_READ");
  if ( value > 0 )
    {
      switch (value)
	{
	case FILE_TYPE_OPEN:
	case FILE_TYPE_FOPEN:
	  FileTypeRead = (int)value;
	  break;
	default:
	  Warning("File type %d not implemented!", value);
	}
    }

  value = file_getenv("FILE_TYPE_WRITE");
  if ( value > 0 )
    {
      switch (value)
	{
	case FILE_TYPE_OPEN:
	case FILE_TYPE_FOPEN:
	  FileTypeWrite = (int)value;
	  break;
	default:
	  Warning("File type %d not implemented!", value);
	}
    }

#if defined (O_NONBLOCK)
  FileFlagWrite = O_NONBLOCK;
#endif
  envString = getenv("FILE_FLAG_WRITE");
  if ( envString )
    {
#if defined (O_NONBLOCK)
      if ( strcmp(envString, "NONBLOCK") == 0 ) FileFlagWrite = O_NONBLOCK;
#endif
    }

  value = file_getenv("FILE_BUFTYPE");
#if ! defined (HAVE_MMAP)
  if ( value == FILE_BUFTYPE_MMAP )
    {
      Warning("MMAP not available!");
      value = 0;
    }
#endif
  if ( value > 0 )
    {
      switch (value)
	{
	case FILE_BUFTYPE_STD:
	case FILE_BUFTYPE_MMAP:
	  FileBufferTypeEnv = (int)value;
	  break;
	default:
	  Warning("File buffer type %d not implemented!", value);
	}
    }

  file_list_new();
  atexit(file_list_delete);

  FILE_LOCK();

  file_init_pointer();

  FILE_UNLOCK();

  if ( FILE_Debug ) atexit(file_table_print);

  _file_init = TRUE;
}

static
void file_set_buffer(bfile_t *fileptr)
{
  size_t buffersize = 0;

  if ( fileptr->mode == 'r' )
    {
      if ( FileBufferTypeEnv )
	fileptr->bufferType = FileBufferTypeEnv;
      else if ( fileptr->bufferType == 0 )
	fileptr->bufferType = FILE_BUFTYPE_STD;

      if ( FileBufferSizeEnv >= 0 )
	buffersize = (size_t) FileBufferSizeEnv;
      else if ( fileptr->bufferSize > 0 )
	buffersize = fileptr->bufferSize;
      else
	{
	  buffersize = fileptr->blockSize * 4;
	  if ( buffersize < FileBufferSizeMin ) buffersize = FileBufferSizeMin;
	}

      if ( (size_t) fileptr->size < buffersize )
	buffersize = (size_t) fileptr->size;

      if ( fileptr->bufferType == FILE_BUFTYPE_MMAP )
	{
	  size_t blocksize = (size_t) pagesize();
	  size_t minblocksize = 4 * blocksize;
	  buffersize = buffersize - buffersize % minblocksize;

	  if ( buffersize < (size_t) fileptr->size && buffersize < minblocksize )
	    buffersize = minblocksize;
	}

      if ( buffersize == 0 ) buffersize = 1;
    }
  else
    {
      fileptr->bufferType = FILE_BUFTYPE_STD;

      if ( FileBufferSizeEnv >= 0 )
	buffersize = (size_t) FileBufferSizeEnv;
      else if ( fileptr->bufferSize > 0 )
	buffersize = fileptr->bufferSize;
      else
	{
	  buffersize = fileptr->blockSize * 4;
	  if ( buffersize < FileBufferSizeMin ) buffersize = FileBufferSizeMin;
	}
    }

  if ( fileptr->bufferType == FILE_BUFTYPE_STD || fileptr->type == FILE_TYPE_FOPEN )
    {
      if ( buffersize > 0 )
        {
          fileptr->buffer = (char *) malloc(buffersize);
          if ( fileptr->buffer == NULL )
            SysError("Allocation of file buffer failed!");
        }
    }

  if ( fileptr->type == FILE_TYPE_FOPEN )
    if ( setvbuf(fileptr->fp, fileptr->buffer, fileptr->buffer ? _IOFBF : _IONBF, buffersize) )
      SysError("setvbuf failed!");

  fileptr->bufferSize = buffersize;
}

static
int file_fill_buffer(bfile_t *fileptr)
{
  ssize_t nread;
  int fd;
  long offset = 0;
  off_t retseek;

  if ( FILE_Debug )
    Message("file ptr = %p  Cnt = %ld", fileptr, fileptr->bufferCnt);

  if ( (fileptr->flag & FILE_EOF) != 0 ) return (EOF);

  if ( fileptr->buffer == NULL ) file_set_buffer(fileptr);

  if ( fileptr->bufferSize == 0 ) return (EOF);

  fd = fileptr->fd;

#if defined (HAVE_MMAP)
  if ( fileptr->bufferType == FILE_BUFTYPE_MMAP )
    {
      if ( fileptr->bufferPos >= fileptr->size )
	{
	  nread = 0;
	}
      else
	{
          xassert(fileptr->bufferSize <= SSIZE_MAX);
	  nread = (ssize_t)fileptr->bufferSize;
	  if ( (nread + fileptr->bufferPos) > fileptr->size )
	    nread = fileptr->size - fileptr->bufferPos;

	  if ( fileptr->buffer )
	    {
              int ret;
	      ret = munmap(fileptr->buffer, fileptr->mappedSize);
	      if ( ret == -1 ) SysError("munmap error for read %s", fileptr->name);
	      fileptr->buffer = NULL;
	    }

	  fileptr->mappedSize = (size_t)nread;

	  fileptr->buffer = (char*) mmap(NULL, (size_t) nread, PROT_READ, MAP_PRIVATE, fd, fileptr->bufferPos);

	  if ( fileptr->buffer == MAP_FAILED ) SysError("mmap error for read %s", fileptr->name);

	  offset = fileptr->position - fileptr->bufferPos;
	}
    }
  else
#endif
    {
      retseek = lseek(fileptr->fd, fileptr->bufferPos, SEEK_SET);
      if ( retseek == (off_t)-1 )
	SysError("lseek error at pos %ld file %s", (long) fileptr->bufferPos, fileptr->name);

      nread = read(fd, fileptr->buffer, fileptr->bufferSize);
    }

  if ( nread <= 0 )
    {
      if ( nread == 0 )
	fileptr->flag |= FILE_EOF;
      else
	fileptr->flag |= FILE_ERROR;

      fileptr->bufferCnt = 0;
      return (EOF);
    }

  fileptr->bufferPtr = fileptr->buffer;
  fileptr->bufferCnt = (size_t)nread;

  fileptr->bufferStart = fileptr->bufferPos;
  fileptr->bufferPos  += nread;
  fileptr->bufferEnd   = fileptr->bufferPos - 1;

  if ( FILE_Debug )
    {
      Message("fileID = %d  Val     = %d",  fileptr->self, (int) fileptr->buffer[0]);
      Message("fileID = %d  Start   = %ld", fileptr->self, fileptr->bufferStart);
      Message("fileID = %d  End     = %ld", fileptr->self, fileptr->bufferEnd);
      Message("fileID = %d  nread   = %ld", fileptr->self, nread);
      Message("fileID = %d  offset  = %ld", fileptr->self, offset);
      Message("fileID = %d  Pos     = %ld", fileptr->self, fileptr->bufferPos);
      Message("fileID = %d  postion = %ld", fileptr->self, fileptr->position);
    }

  if ( offset > 0 )
    {
      if ( offset > nread )
	Error("Internal problem with buffer handling. nread = %d offset = %d", nread, offset);

      fileptr->bufferPtr += offset;
      fileptr->bufferCnt -= (size_t)offset;
    }

  fileptr->bufferNumFill++;

  return ((unsigned char) *fileptr->bufferPtr);
}

static
void file_copy_from_buffer(bfile_t *fileptr, void *ptr, size_t size)
{
  if ( FILE_Debug )
    Message("size = %ld  Cnt = %ld", size, fileptr->bufferCnt);

  if ( fileptr->bufferCnt < size )
    Error("Buffer too small. bufferCnt = %d", fileptr->bufferCnt);

  if ( size == 1 )
    {
      ((char *)ptr)[0] = fileptr->bufferPtr[0];

      fileptr->bufferPtr++;
      fileptr->bufferCnt--;
    }
  else
    {
      memcpy(ptr, fileptr->bufferPtr, size);

      fileptr->bufferPtr += size;
      fileptr->bufferCnt -= size;
    }
}

static
size_t file_read_from_buffer(bfile_t *fileptr, void *ptr, size_t size)
{
  size_t nread, rsize;
  size_t offset = 0;

  if ( FILE_Debug )
    Message("size = %ld  Cnt = %ld", size, (long) fileptr->bufferCnt);

  if ( ((long)fileptr->bufferCnt) < 0L )
    Error("Internal problem. bufferCnt = %ld", (long) fileptr->bufferCnt);

  rsize = size;

  while ( fileptr->bufferCnt < rsize )
    {
      nread = fileptr->bufferCnt;
      /*
      fprintf(stderr, "rsize = %d nread = %d\n", (int) rsize, (int) nread);
      */
      if ( nread > (size_t) 0 )
	file_copy_from_buffer(fileptr, (char *)ptr+offset, nread);
      offset += nread;
      if ( nread < rsize )
	rsize -= nread;
      else
	rsize = 0;

      if ( file_fill_buffer(fileptr) == EOF ) break;
    }

  nread = size - offset;

  if ( fileptr->bufferCnt < nread ) nread = fileptr->bufferCnt;

  if ( nread > (unsigned) 0 )
    file_copy_from_buffer(fileptr, (char *)ptr+offset, nread);

  return (nread+offset);
}


void fileSetBufferSize(int fileID, long buffersize)
{
  bfile_t *fileptr = file_to_pointer(fileID);
  xassert(buffersize >= 0);
  if ( fileptr ) fileptr->bufferSize = (size_t)buffersize;
}

/*
 *   Open a file. Returns file ID, or -1 on error
 */
int fileOpen(const char *filename, const char *mode)
{
  int (*myFileOpen)(const char *filename, const char *mode)
    = (int (*)(const char *, const char *))
    namespaceSwitchGet(NSSWITCH_FILE_OPEN).func;
  return myFileOpen(filename, mode);
}

int fileOpen_serial(const char *filename, const char *mode)
{
  FILE *fp = NULL;    /* file pointer    (used for write) */
  int fd = -1;        /* file descriptor (used for read)  */
  int fileID = FILE_UNDEFID;
  int fmode = 0;
  struct stat filestat;
  bfile_t *fileptr = NULL;

  FILE_INIT();

  fmode = tolower((int) mode[0]);

  switch ( fmode )
    {
    case 'r':
      if ( FileTypeRead == FILE_TYPE_FOPEN )
	fp = fopen(filename, "rb");
      else
	fd =  open(filename, O_RDONLY | O_BINARY);
      break;
    case 'x':  fp = fopen(filename, "rb");      break;
    case 'w':
      if ( FileTypeWrite == FILE_TYPE_FOPEN )
        fp = fopen(filename, "wb");
      else
	fd =  open(filename, O_CREAT | O_TRUNC | O_WRONLY | O_BINARY | FileFlagWrite, 0666);
      break;
    case 'a':  fp = fopen(filename, "ab");      break;
    default:   Error("Mode %c unexpected!", fmode);
    }

  if ( FILE_Debug )
    if ( fp == NULL && fd == -1 )
      Message("Open failed on %s mode %c errno %d", filename, fmode, errno);

  if ( fp )
    {
      if ( stat(filename, &filestat) != 0 ) return (fileID);

      fileptr = file_new_entry();
      if ( fileptr )
	{
	  fileID = fileptr->self;
	  fileptr->fp = fp;
	}
    }
  else if ( fd >= 0 )
    {
      if ( fstat(fd, &filestat) != 0 ) return (fileID);

      fileptr = file_new_entry();
      if ( fileptr )
	{
	  fileID = fileptr->self;
	  fileptr->fd = fd;
	}
    }

  if ( fileID >= 0 )
    {
      fileptr->mode = fmode;
      fileptr->name = strdupx(filename);

#if defined (HAVE_STRUCT_STAT_ST_BLKSIZE)
      fileptr->blockSize = (size_t) filestat.st_blksize;
#else
      fileptr->blockSize = (size_t) 4096;
#endif

      if ( fmode == 'r' )
        fileptr->type = FileTypeRead;
      else if ( fmode == 'w' )
        fileptr->type = FileTypeWrite;
      else
	fileptr->type = FILE_TYPE_FOPEN;

      if ( fmode == 'r' ) fileptr->size = filestat.st_size;

      if ( fileptr->type == FILE_TYPE_FOPEN ) file_set_buffer(fileptr);

      if ( FILE_Debug )
	Message("File %s opened with ID %d", filename, fileID);
    }

  return (fileID);
}

/*
 *   Close a file.
 */
int fileClose(int fileID)
{
  int (*myFileClose)(int fileID)
    = (int (*)(int))namespaceSwitchGet(NSSWITCH_FILE_CLOSE).func;
  return myFileClose(fileID);
}

int fileClose_serial(int fileID)
{
  char *name;
  int ret;
  char *fbtname[] = {"unknown", "standard", "mmap"};
  char *ftname[] = {"unknown", "open", "fopen"};
  bfile_t *fileptr = file_to_pointer(fileID);
  double rout = 0;

  if ( fileptr == NULL )
    {
      file_pointer_info(__func__, fileID);
      return (1);
    }

  name = fileptr->name;

  if ( FILE_Debug )
    Message("fileID = %d  filename = %s", fileID, name);

  if ( FileInfo > 0 )
    {
      fprintf(stderr, "____________________________________________\n");
      fprintf(stderr, " file ID          : %d\n",  fileID);
      fprintf(stderr, " file name        : %s\n",  fileptr->name);
      fprintf(stderr, " file type        : %d (%s)\n", fileptr->type, ftname[fileptr->type]);

      if ( fileptr->type == FILE_TYPE_FOPEN )
	fprintf(stderr, " file pointer     : %p\n",  (void *) fileptr->fp);
      else
        {
          fprintf(stderr, " file descriptor  : %d\n",  fileptr->fd);
          fprintf(stderr, " file flag        : %d\n", FileFlagWrite);
        }
      fprintf(stderr, " file mode        : %c\n",  fileptr->mode);

      if ( sizeof(off_t) > sizeof(long) )
	{
#if defined (_WIN32)
	  fprintf(stderr, " file size        : %I64d\n", (long long) fileptr->size);
	  if ( fileptr->type == FILE_TYPE_OPEN )
	    fprintf(stderr, " file position    : %I64d\n", (long long) fileptr->position);
	  fprintf(stderr, " bytes transfered : %I64d\n", (long long) fileptr->byteTrans);
#else
	  fprintf(stderr, " file size        : %lld\n", (long long) fileptr->size);
	  if ( fileptr->type == FILE_TYPE_OPEN )
	    fprintf(stderr, " file position    : %lld\n", (long long) fileptr->position);
	  fprintf(stderr, " bytes transfered : %lld\n", (long long) fileptr->byteTrans);
#endif
	}
      else
	{
	  fprintf(stderr, " file size        : %ld\n", (long) fileptr->size);
	  if ( fileptr->type == FILE_TYPE_OPEN )
	    fprintf(stderr, " file position    : %ld\n", (long) fileptr->position);
	  fprintf(stderr, " bytes transfered : %ld\n", (long) fileptr->byteTrans);
	}

      if ( fileptr->time_in_sec > 0 )
        {
          rout = (double)fileptr->byteTrans / (1024.*1024.*fileptr->time_in_sec);
        }

      fprintf(stderr, " wall time [s]    : %.2f\n", fileptr->time_in_sec);
      fprintf(stderr, " data rate [MB/s] : %.1f\n", rout);

      fprintf(stderr, " file access      : %ld\n", fileptr->access);
      if ( fileptr->mode == 'r' && fileptr->type == FILE_TYPE_OPEN )
	{
	  fprintf(stderr, " buffer type      : %d (%s)\n", fileptr->bufferType, fbtname[fileptr->bufferType]);
	  fprintf(stderr, " num buffer fill  : %ld\n", fileptr->bufferNumFill);
	}
      fprintf(stderr, " buffer size      : %lu\n", (unsigned long) fileptr->bufferSize);
      fprintf(stderr, " block size       : %lu\n", (unsigned long) fileptr->blockSize);
      fprintf(stderr, " page size        : %d\n",  pagesize());
      fprintf(stderr, "--------------------------------------------\n");
    }

  if ( fileptr->type == FILE_TYPE_FOPEN )
    {
      ret = fclose(fileptr->fp);
      if ( ret == EOF )
	SysError("EOF returned for close of %s!", name);
    }
  else
    {
#if defined (HAVE_MMAP)
      if ( fileptr->buffer && fileptr->mappedSize )
	{
	  ret = munmap(fileptr->buffer, fileptr->mappedSize);
	  if ( ret == -1 ) SysError("munmap error for close %s", fileptr->name);
	  fileptr->buffer = NULL;
	}
#endif
      ret = close(fileptr->fd);
      if ( ret == -1 )
	SysError("EOF returned for close of %s!", name);
    }

  if ( fileptr->name )    free((void*) fileptr->name);
  if ( fileptr->buffer )  free((void*) fileptr->buffer);

  file_delete_entry(fileptr);

  return (0);
}


int filePtrGetc(void *vfileptr)
{
  int ivalue = EOF;
  int fillret = 0;
  bfile_t *fileptr = (bfile_t *) vfileptr;

  if ( fileptr )
    {
      if ( fileptr->mode == 'r' && fileptr->type == FILE_TYPE_OPEN )
	{
	  if ( fileptr->bufferCnt == 0 ) fillret = file_fill_buffer(fileptr);

	  if ( fillret >= 0 )
	    {
	      ivalue = (unsigned char) *fileptr->bufferPtr++;
	      fileptr->bufferCnt--;
	      fileptr->position++;

	      fileptr->byteTrans++;
	      fileptr->access++;
	    }
	}
      else
	{
	  ivalue = fgetc(fileptr->fp);
	  if ( ivalue >= 0 )
	    {
	      fileptr->byteTrans++;
	      fileptr->access++;
	    }
	  else
	    fileptr->flag |= FILE_EOF;
	}
    }

  return (ivalue);
}


int fileGetc(int fileID)
{
  int ivalue;
  bfile_t *fileptr;

  fileptr = file_to_pointer(fileID);

  ivalue = filePtrGetc((void *)fileptr);

  return (ivalue);
}


size_t filePtrRead(void *vfileptr, void *restrict ptr, size_t size)
{
  size_t nread = 0;
  bfile_t *fileptr = (bfile_t *) vfileptr;

  if ( fileptr )
    {
      if ( fileptr->mode == 'r' && fileptr->type == FILE_TYPE_OPEN )
	nread = file_read_from_buffer(fileptr, ptr, size);
      else
	{
	  nread = fread(ptr, 1, size, fileptr->fp);
	  if ( nread != size )
	    {
	      if ( nread == 0 )
		fileptr->flag |= FILE_EOF;
	      else
		fileptr->flag |= FILE_ERROR;
	    }
	}

      fileptr->position  += (off_t)nread;
      fileptr->byteTrans += (off_t)nread;
      fileptr->access++;
    }

  if ( FILE_Debug ) Message("size %ld  nread %ld", size, nread);

  return (nread);
}


size_t fileRead(int fileID, void *restrict ptr, size_t size)
{
  size_t nread = 0;
  bfile_t *fileptr;

  fileptr = file_to_pointer(fileID);

  if ( fileptr )
    {
      double t_begin = 0.0;

      if ( FileInfo ) t_begin = file_time();

      if ( fileptr->type == FILE_TYPE_OPEN )
	nread = file_read_from_buffer(fileptr, ptr, size);
      else
	{
	  nread = fread(ptr, 1, size, fileptr->fp);
	  if ( nread != size )
	    {
	      if ( nread == 0 )
		fileptr->flag |= FILE_EOF;
	      else
		fileptr->flag |= FILE_ERROR;
	    }
	}

      if ( FileInfo ) fileptr->time_in_sec += file_time() - t_begin;

      fileptr->position  += (off_t)nread;
      fileptr->byteTrans += (off_t)nread;
      fileptr->access++;
    }

  if ( FILE_Debug ) Message("size %ld  nread %ld", size, nread);

  return (nread);
}


size_t fileWrite(int fileID, const void *restrict ptr, size_t size)
{
  size_t nwrite = 0;
  bfile_t *fileptr;

  fileptr = file_to_pointer(fileID);

  if ( fileptr )
    {
      double t_begin = 0.0;

      /* if ( fileptr->buffer == NULL ) file_set_buffer(fileptr); */

      if ( FileInfo ) t_begin = file_time();

      if ( fileptr->type == FILE_TYPE_FOPEN )
        nwrite = fwrite(ptr, 1, size, fileptr->fp);
      else
        {
          ssize_t temp = write(fileptr->fd, ptr, size);
          if (temp == -1)
            {
              perror("error writing to file");
              nwrite = 0;
            }
          else
            nwrite = (size_t)temp;
        }

      if ( FileInfo ) fileptr->time_in_sec += file_time() - t_begin;

      fileptr->position  += (off_t)nwrite;
      fileptr->byteTrans += (off_t)nwrite;
      fileptr->access++;
    }

  return (nwrite);
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
