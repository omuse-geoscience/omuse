#ifndef _DMEMORY_H
#define _DMEMORY_H

#include <stdlib.h>

/*
 * if DEBUG_MEMORY is defined setenv MEMORY_DEBUG to debug memory
 */

#define  DEBUG_MEMORY

#ifndef  WITH_CALLER_NAME
#define  WITH_CALLER_NAME
#endif

extern size_t  memTotal(void);
extern void    memDebug(int debug);
extern void    memExitOnError(void);

#if  defined  DEBUG_MEMORY

extern void   *Realloc(const char *caller, const char *file, int line, void *ptr, size_t size);
extern void   *Calloc (const char *caller, const char *file, int line, size_t nmemb, size_t size);
extern void   *Malloc (const char *caller, const char *file, int line, size_t size);
extern void    Free   (const char *caller, const char *file, int line, void *ptr);

#if  defined  calloc
#  undef  calloc
#endif

#if  defined  WITH_CALLER_NAME
#  define  realloc(p, s)  Realloc(__func__, __FILE__, __LINE__, (p), (s))
#  define   calloc(n, s)   Calloc(__func__, __FILE__, __LINE__, (n), (s))
#  define   malloc(s)      Malloc(__func__, __FILE__, __LINE__, (s))
#  define     free(p)        Free(__func__, __FILE__, __LINE__, (p))
#else
#  define  realloc(p, s)  Realloc((void *) NULL, __FILE__, __LINE__, (p), (s))
#  define   calloc(n, s)   Calloc((void *) NULL, __FILE__, __LINE__, (n), (s))
#  define   malloc(s)      Malloc((void *) NULL, __FILE__, __LINE__, (s))
#  define     free(p)        Free((void *) NULL, __FILE__, __LINE__, (p))
#endif

#endif /* DEBUG_MEMORY */

#endif /* _DMEMORY_H */
