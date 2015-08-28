#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <errno.h>

#if !defined (NAMESPACE_H)
#include "namespace.h"
#endif

int _ExitOnError   = 1;	/* If set to 1, exit on error       */
int _Verbose = 1;	/* If set to 1, errors are reported */
int _Debug   = 0;       /* If set to 1, debugging           */


void SysError_(const char *caller, const char *fmt, ...)
{
  va_list args;

  va_start(args, fmt);

  printf("\n");
   fprintf(stderr, "Error (%s) : ", caller);
  vfprintf(stderr, fmt, args);
   fprintf(stderr, "\n");

  va_end(args);

  if ( errno )
    perror("System error message ");

  exit(EXIT_FAILURE);
}


void Error_(const char *caller, const char *fmt, ...)
{
  va_list args;

  va_start(args, fmt);

  printf("\n");
   fprintf(stderr, "Error (%s) : ", caller);
  vfprintf(stderr, fmt, args);
   fprintf(stderr, "\n");

  va_end(args);

  if ( _ExitOnError ) exit(EXIT_FAILURE);
}

typedef void (*cdiAbortCFunc)(const char * caller, const char * filename,
                              const char *functionname, int line,
                              const char * errorString, va_list ap)
#ifdef __GNUC__
  __attribute__((noreturn))
#endif
;

void cdiAbortC(const char * caller, const char * filename,
               const char *functionname, int line,
               const char * errorString, ... )
{
  va_list ap;
  va_start(ap, errorString);
  cdiAbortCFunc cdiAbortC_p
    = (cdiAbortCFunc)namespaceSwitchGet(NSSWITCH_ABORT).func;
  cdiAbortC_p(caller, filename, functionname, line, errorString, ap);
  va_end(ap);
}

void
cdiAbortC_serial(const char *caller, const char *filename,
                 const char *functionname, int line,
                 const char *errorString, va_list ap)
{
  fprintf(stderr, "ERROR, %s, %s, line %d%s%s\nerrorString: \"",
          functionname, filename, line, caller?", called from ":"",
          caller?caller:"");
  vfprintf(stderr, errorString, ap);
  fputs("\"\n", stderr);
  exit(EXIT_FAILURE);
}

typedef void (*cdiWarningFunc)(const char * caller, const char * fmt,
                               va_list ap);

void Warning_(const char *caller, const char *fmt, ...)
{
  va_list args;

  va_start(args, fmt);

  if ( _Verbose )
    {
      cdiWarningFunc cdiWarning_p
        = (cdiWarningFunc)namespaceSwitchGet(NSSWITCH_WARNING).func;
      cdiWarning_p(caller, fmt, args);
    }

  va_end(args);
}

void cdiWarning(const char *caller, const char *fmt, va_list ap)
{
  fprintf(stderr, "Warning (%s) : ", caller);
  vfprintf(stderr, fmt, ap);
  fputc('\n', stderr);
}


void Message_(const char *caller, const char *fmt, ...)
{
  va_list args;

  va_start(args, fmt);

   fprintf(stdout, "%-18s : ", caller);
  vfprintf(stdout, fmt, args);
   fprintf(stdout, "\n");

  va_end(args);
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
