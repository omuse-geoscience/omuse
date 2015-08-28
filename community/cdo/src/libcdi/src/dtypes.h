#ifndef _DTYPES_H
#define _DTYPES_H

#include <stdio.h>
#include <limits.h>

/* INT32 */

#if ! defined (INT_MAX)
#  error INT_MAX undefined
#endif

#undef  INT32
#if  INT_MAX == 2147483647L
#  define  INT32  int
#elif LONG_MAX == 2147483647L
#  define  INT32  long
#endif

/* INT64 */

#if ! defined (LONG_MAX)
#  error LONG_MAX undefined
#endif

#undef  INT64
#if  LONG_MAX > 2147483647L
#  define  INT64  long
#else
#  define  INT64  long long
#endif

/* FLT32 */

#undef   FLT32
#define  FLT32  float

/* FLT64 */

#undef   FLT64
#define  FLT64  double

/* UINT32 and UINT64 */

#define  UINT32   unsigned INT32
#define  UINT64   unsigned INT64

#endif  /* _DTYPES_H */
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
