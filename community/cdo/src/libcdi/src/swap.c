#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <inttypes.h>
#include <stdio.h>

#include "error.h"
#include "binary.h"

void swap4byte(void *ptr, size_t size)
{
  int32_t *ptrtmp = ptr;

  if (sizeof (int32_t) == 4)
    for (size_t i = 0; i < size; ++i)
      ptrtmp[i] = (((ptrtmp[i] >> 24) & 0x00ff) | ((ptrtmp[i] & 0x00ff) << 24) |
                   ((ptrtmp[i] >>  8) & 0xff00) | ((ptrtmp[i] & 0xff00) <<  8));
  else
    Error("not implemented for %d byte data", sizeof(int32_t));
}

void swap8byte(void *ptr, size_t size)
{
  int64_t *ptrtmp = ptr;

  if (sizeof (int64_t) == 8)
    for (size_t i = 0; i < size; ++i)
      ptrtmp[i] = (((ptrtmp[i] >> 56) & 0x000000ff) | ((ptrtmp[i] & 0x000000ff) << 56) |
                   ((ptrtmp[i] >> 40) & 0x0000ff00) | ((ptrtmp[i] & 0x0000ff00) << 40) |
                   ((ptrtmp[i] >> 24) & 0x00ff0000) | ((ptrtmp[i] & 0x00ff0000) << 24) |
                   ((ptrtmp[i] >>  8) & 0xff000000) | ((ptrtmp[i] & 0xff000000) <<  8));
  else
    Error("not implemented for %d byte data", sizeof(int64_t));
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
