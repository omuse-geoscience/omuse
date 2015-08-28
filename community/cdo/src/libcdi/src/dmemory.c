#if  defined(HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdarg.h>
#include <errno.h>

#if !defined(HAVE_CONFIG_H) && !defined(HAVE_MALLOC_H) && defined(SX)
#  define  HAVE_MALLOC_H
#endif

#if  defined(HAVE_MALLOC_H)
#  include <malloc.h>
#endif

//There is no point in avoiding to include our own header, it is likely to be included before this file in the backed cdilib.c.
//As such, we may as well always include it and let the compiler check our function signatures.
#include "dmemory.h"

//However, we need to avoid clobbering our own `malloc()` calls, so we ensure that our own malloc calls cannot be interpreted as macro calls.
#define protected_malloc (malloc)
#define protected_calloc (calloc)
#define protected_realloc (realloc)
#define protected_free (free)

enum             {MALLOC_FUNC=0, CALLOC_FUNC, REALLOC_FUNC, FREE_FUNC};
static const char *memfunc[] = {"Malloc", "Calloc", "Realloc", "Free"};

#undef   MEM_UNDEFID
#define  MEM_UNDEFID  -1

#define  MEM_MAXNAME  32   /* Min = 8, for  "unknown" ! */

int dmemory_ExitOnError = 0;

typedef struct
{
  void     *ptr;
  int       item;
  size_t    size;
  size_t    nobj;
  int       mtype;
  int       line;
  char      file[MEM_MAXNAME];
  char      caller[MEM_MAXNAME];
}
MemTable_t;

static MemTable_t *memTable;
static size_t  memTableSize  = 0;
static long    memAccess     = 0;

static size_t  MemObjs       = 0;
static size_t  MaxMemObjs    = 0;
static size_t  MemUsed       = 0;
static size_t  MaxMemUsed    = 0;

static int     MEM_Debug     = 0;   /* If set to 1, debugging */
static int     MEM_Info      = 0;   /* If set to 1, print mem table at exit */

void memDebug(int debug)
{
  MEM_Debug = debug;
}

static
void memInternalProblem(const char *caller, const char *fmt, ...)
{
  va_list args;

  va_start(args, fmt);

  printf("\n");
   fprintf(stderr, "Internal problem (%s) : ", caller);
  vfprintf(stderr, fmt, args);
   fprintf(stderr, "\n");

  va_end(args);

  exit(EXIT_FAILURE);
}

static
void memError(const char *caller, const char *file, int line, size_t size)
{
  printf("\n");
  fprintf(stderr, "Error (%s) : Allocation of %zu bytes failed. [ line %d file %s ]\n",
	  caller, size, line, file);

  if ( errno ) perror("System error message ");

  exit(EXIT_FAILURE);
}

static
void memListPrintEntry(int mtype, int item, size_t size, void *ptr,
		       const char *caller, const char *file, int line)
{
  fprintf(stderr, "[%-7s ", memfunc[mtype]);

  fprintf(stderr, "memory item %3d ", item);
  fprintf(stderr, "(%6zu byte) ", size);
  fprintf(stderr, "at %p", ptr);
  if ( file != NULL )
    {
      fprintf(stderr, " line %4d", line);
      fprintf(stderr, " file %s", file);
    }
  if ( caller != NULL )
    fprintf(stderr, " (%s)", caller);
  fprintf(stderr, "]\n");
}

static
void memListPrintTable(void)
{
  if ( MemObjs ) fprintf(stderr, "\nMemory table:\n");

  for ( size_t memID = 0; memID < memTableSize; memID++ )
    {
      if ( memTable[memID].item != MEM_UNDEFID )
        memListPrintEntry(memTable[memID].mtype, memTable[memID].item,
                          memTable[memID].size*memTable[memID].nobj,
                          memTable[memID].ptr, memTable[memID].caller,
                          memTable[memID].file, memTable[memID].line);
    }

  if ( MemObjs )
    {
      fprintf(stderr, "  Memory access             : %6u\n", (unsigned) memAccess);
      fprintf(stderr, "  Maximum objects           : %6zu\n", memTableSize);
      fprintf(stderr, "  Objects used              : %6u\n", (unsigned) MaxMemObjs);
      fprintf(stderr, "  Objects in use            : %6u\n", (unsigned) MemObjs);
      fprintf(stderr, "  Memory allocated          : ");
      if (MemUsed > 1024*1024*1024)
	fprintf(stderr, " %5d GB\n",   (int) (MemUsed/(1024*1024*1024)));
      else if (MemUsed > 1024*1024)
	fprintf(stderr, " %5d MB\n",   (int) (MemUsed/(1024*1024)));
      else if (MemUsed > 1024)
	fprintf(stderr, " %5d KB\n",   (int) (MemUsed/(1024)));
      else
	fprintf(stderr, " %5d Byte\n", (int)  MemUsed);
    }

  if ( MaxMemUsed )
    {
      fprintf(stderr, "  Maximum memory allocated  : ");
      if (MaxMemUsed > 1024*1024*1024)
	fprintf(stderr, " %5d GB\n",   (int) (MaxMemUsed/(1024*1024*1024)));
      else if (MaxMemUsed > 1024*1024)
	fprintf(stderr, " %5d MB\n",   (int) (MaxMemUsed/(1024*1024)));
      else if (MaxMemUsed > 1024)
	fprintf(stderr, " %5d KB\n",   (int) (MaxMemUsed/(1024)));
      else
	fprintf(stderr, " %5d Byte\n", (int)  MaxMemUsed);
    }
}

static
void memGetDebugLevel(void)
{
  const char *envstr;

  envstr = getenv("MEMORY_INFO");
  if ( envstr && isdigit((int) envstr[0]) ) MEM_Info = atoi(envstr);

  envstr = getenv("MEMORY_DEBUG");

  if ( envstr && isdigit((int) envstr[0]) ) MEM_Debug = atoi(envstr);

  if ( MEM_Debug && !MEM_Info ) MEM_Info = 1;

  if ( MEM_Info ) atexit(memListPrintTable);
}

static
void memInit(void)
{
  static int initDebugLevel = 0;

  if ( ! initDebugLevel )
    {
      memGetDebugLevel();
      initDebugLevel = 1;
    }
}

static
int memListDeleteEntry(void *ptr, size_t *size)
{
  int item = MEM_UNDEFID;
  size_t memID;

  for ( memID = 0; memID < memTableSize; memID++ )
    {
      if ( memTable[memID].item == MEM_UNDEFID ) continue;
      if ( memTable[memID].ptr == ptr ) break;
    }

  if ( memID != memTableSize )
    {
      MemObjs--;
      MemUsed -= memTable[memID].size * memTable[memID].nobj;
      *size = memTable[memID].size * memTable[memID].nobj;
      item = memTable[memID].item;
      memTable[memID].item = MEM_UNDEFID;
    }

  return item;
}

static
void memTableInitEntry(size_t memID)
{
  if ( memID >= memTableSize )
    memInternalProblem(__func__, "memID %d undefined!", memID);

  memTable[memID].ptr    = NULL;
  memTable[memID].item   = MEM_UNDEFID;
  memTable[memID].size   = 0;
  memTable[memID].nobj   = 0;
  memTable[memID].mtype  = MEM_UNDEFID;
  memTable[memID].line   = MEM_UNDEFID;
}

static
int memListNewEntry(int mtype, void *ptr, size_t size, size_t nobj,
		    const char *caller, const char *file, int line)
{
  static int item = 0;
  size_t memSize = 0;
  size_t memID = 0;

  /*
    Look for a free slot in memTable.
    (Create the table the first time through).
  */
  if ( memTableSize == 0 )
    {
      memTableSize = 8;
      memSize  = memTableSize * sizeof(MemTable_t);
      memTable = (MemTable_t *) protected_malloc(memSize);
      if( memTable == NULL ) memError(__func__, __FILE__, __LINE__, memSize);

      for ( size_t i = 0; i < memTableSize; i++ )
	memTableInitEntry(i);
    }
  else
    {
      while ( memID < memTableSize )
	{
	  if ( memTable[memID].item == MEM_UNDEFID ) break;
	  memID++;
	}
    }
  /*
    If the table overflows, double its size.
  */
  if ( memID == memTableSize )
    {
      memTableSize = 2*memTableSize;
      memSize  = memTableSize*sizeof(MemTable_t);
      memTable = (MemTable_t*) protected_realloc(memTable, memSize);
      if ( memTable == NULL ) memError(__func__, __FILE__, __LINE__, memSize);

      for ( size_t i = memID; i < memTableSize; i++ )
	memTableInitEntry(i);
    }

  memTable[memID].item  = item;
  memTable[memID].ptr   = ptr;
  memTable[memID].size  = size;
  memTable[memID].nobj  = nobj;
  memTable[memID].mtype = mtype;
  memTable[memID].line  = line;

  if ( file )
    {
      size_t len = strlen(file);
      if ( len > MEM_MAXNAME-1 ) len = MEM_MAXNAME-1;

      (void) memcpy(memTable[memID].file, file, len);
      memTable[memID].file[len] = '\0';
    }
  else
    {
      (void) strcpy(memTable[memID].file, "unknown");
    }

  if ( caller )
    {
      size_t len = strlen(caller);
      if ( len > MEM_MAXNAME-1 ) len = MEM_MAXNAME-1;

      (void) memcpy(memTable[memID].caller, caller, len);
      memTable[memID].caller[len] = '\0';
    }
  else
    {
      (void) strcpy(memTable[memID].caller, "unknown");
    }

  MaxMemObjs++;
  MemObjs++;
  MemUsed += size*nobj;
  if ( MemUsed > MaxMemUsed ) MaxMemUsed = MemUsed;

  return item++;
}

static
int memListChangeEntry(void *ptrold, void *ptr, size_t size,
		       const char *caller, const char *file, int line)
{
  int item = MEM_UNDEFID;
  size_t memID = 0;
  size_t sizeold;

  while( memID < memTableSize )
    {
      if ( memTable[memID].item != MEM_UNDEFID )
	if ( memTable[memID].ptr == ptrold ) break;
      memID++;
    }

  if ( memID == memTableSize )
    {
      if ( ptrold != NULL )
	memInternalProblem(__func__, "Item at %p not found.", ptrold);
    }
  else
    {
      item = memTable[memID].item;

      sizeold = memTable[memID].size*memTable[memID].nobj;

      memTable[memID].ptr   = ptr;
      memTable[memID].size  = size;
      memTable[memID].nobj  = 1;
      memTable[memID].mtype = REALLOC_FUNC;
      memTable[memID].line  = line;

      if ( file )
	{
	  size_t len = strlen(file);
	  if ( len > MEM_MAXNAME-1 ) len = MEM_MAXNAME-1;

	  (void) memcpy(memTable[memID].file, file, len);
	  memTable[memID].file[len] = '\0';
	}
      else
	{
	  (void) strcpy(memTable[memID].file, "unknown");
	}

      if ( caller )
	{
	  size_t len = strlen(caller);
	  if ( len > MEM_MAXNAME-1 ) len = MEM_MAXNAME-1;

	  (void) memcpy(memTable[memID].caller, caller, len);
	  memTable[memID].caller[len] = '\0';
	}
      else
	{
	  (void) strcpy(memTable[memID].caller, "unknown");
	}

      MemUsed -= sizeold;
      MemUsed += size;
      if ( MemUsed > MaxMemUsed ) MaxMemUsed = MemUsed;
    }

  return item;
}


void *Calloc(const char *caller, const char *file, int line, size_t nobjs, size_t size)
{
  void *ptr = NULL;
  int item = MEM_UNDEFID;

  memInit();

  if ( nobjs*size > 0 )
    {
      ptr = protected_calloc(nobjs, size);

      if ( MEM_Info )
	{
	  memAccess++;

	  if ( ptr ) item = memListNewEntry(CALLOC_FUNC, ptr, size, nobjs, caller, file, line);

	  if ( MEM_Debug ) memListPrintEntry(CALLOC_FUNC, item, size*nobjs, ptr, caller, file, line);
	}

      if ( ptr == NULL && dmemory_ExitOnError )
	memError(caller, file, line, size*nobjs);
    }
  else
    fprintf(stderr, "Warning (%s) : Allocation of 0 bytes! [ line %d file %s ]\n", caller, line, file);

  return(ptr);
}


void *Malloc(const char *caller, const char *file, int line, size_t size)
{
  void *ptr = NULL;
  int item = MEM_UNDEFID;

  memInit();

  if ( size > 0 )
    {
      ptr = protected_malloc(size);

      if ( MEM_Info )
	{
	  memAccess++;

	  if ( ptr ) item = memListNewEntry(MALLOC_FUNC, ptr, size, 1, caller, file, line);

	  if ( MEM_Debug ) memListPrintEntry(MALLOC_FUNC, item, size, ptr, caller, file, line);
	}

      if ( ptr == NULL && dmemory_ExitOnError )
	memError(caller, file, line, size);
    }
  else
    fprintf(stderr, "Warning (%s) : Allocation of 0 bytes! [ line %d file %s ]\n", caller, line, file);

  return ptr;
}


void *Realloc(const char *caller, const char *file, int line, void *ptrold, size_t size)
{
  void *ptr = NULL;
  int item = MEM_UNDEFID;

  memInit();

  if ( size > 0 )
    {
      ptr = protected_realloc(ptrold, size);

      if ( MEM_Info )
	{
	  memAccess++;

	  if ( ptr )
	    {
	      item = memListChangeEntry(ptrold, ptr, size, caller, file, line);

	      if ( item == MEM_UNDEFID ) item = memListNewEntry(REALLOC_FUNC, ptr, size, 1, caller, file, line);
	    }

	  if ( MEM_Debug ) memListPrintEntry(REALLOC_FUNC, item, size, ptr, caller, file, line);
	}

      if ( ptr == NULL && dmemory_ExitOnError )
	memError(caller, file, line, size);
    }
  else
    fprintf(stderr, "Warning (%s) : Allocation of 0 bytes! [ line %d file %s ]\n", caller, line, file);

  return (ptr);
}


void Free(const char *caller, const char *file, int line, void *ptr)
{
  int item;
  size_t size;

  memInit();

  if ( MEM_Info )
    {
      if ( (item = memListDeleteEntry(ptr, &size)) >= 0 )
	{
	  if ( MEM_Debug ) memListPrintEntry(FREE_FUNC, item, size, ptr, caller, file, line);
	}
      else
	{
	  if ( ptr && MEM_Debug  )
	    fprintf(stderr, "%s info: memory entry at %p not found. [line %4d file %s (%s)]\n",
		    __func__, ptr, line, file, caller);
	}
    }

  protected_free(ptr);
}

#include "error.h" // cdiAbort()

void *cdiXmalloc(size_t size, const char *filename, const char *functionname,
                 int line)
{
  void *value = protected_malloc(size);
  if (size == 0 || value != NULL) ; else
    cdiAbort(filename, functionname, line, "malloc failed: %s", strerror(errno));
  return value;
}

void *cdiXcalloc(size_t nmemb, size_t size, const char *filename,
                 const char *functionname, int line)
{
  void *value = protected_calloc(nmemb, size);
  if (size == 0 || value != NULL) ; else
    cdiAbort(filename, functionname, line, "calloc failed: %s", strerror(errno));
  return value;
}

void *cdiXrealloc(void *p, size_t size, const char *functionname,
                  const char *filename, int line)
{
  void *value = protected_realloc(p, size);
  if (size == 0 || value != NULL) ; else
    cdiAbort(filename, functionname, line, "realloc failed: %s", strerror(errno));
  return value;
}

size_t memTotal(void)
{
  size_t memtotal = 0;
#if  defined  (HAVE_MALLINFO)
  struct mallinfo meminfo = mallinfo();
  if ( MEM_Debug )
    {
      fprintf(stderr, "arena      %8zu (non-mmapped space allocated from system)\n", (size_t)meminfo.arena);
      fprintf(stderr, "ordblks    %8zu (number of free chunks)\n", (size_t)meminfo.ordblks);
      fprintf(stderr, "smblks     %8zu (number of fastbin blocks)\n", (size_t) meminfo.smblks);
      fprintf(stderr, "hblks      %8zu (number of mmapped regions)\n", (size_t) meminfo.hblks);
      fprintf(stderr, "hblkhd     %8zu (space in mmapped regions)\n", (size_t) meminfo.hblkhd);
      fprintf(stderr, "usmblks    %8zu (maximum total allocated space)\n", (size_t) meminfo.usmblks);
      fprintf(stderr, "fsmblks    %8zu (maximum total allocated space)\n", (size_t) meminfo.fsmblks);
      fprintf(stderr, "uordblks   %8zu (total allocated space)\n", (size_t) meminfo.uordblks);
      fprintf(stderr, "fordblks   %8zu (total free space)\n", (size_t) meminfo.fordblks);
      fprintf(stderr, "Memory in use:   %8zu bytes\n", (size_t) meminfo.usmblks + (size_t)meminfo.uordblks);
      fprintf(stderr, "Total heap size: %8zu bytes\n", (size_t) meminfo.arena);

      /* malloc_stats(); */
    }
  memtotal = (size_t)meminfo.arena;
#endif

  return memtotal;
}


void memExitOnError(void)
{
  dmemory_ExitOnError = 1;
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
