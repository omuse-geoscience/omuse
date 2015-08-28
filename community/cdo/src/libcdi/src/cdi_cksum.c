#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <inttypes.h>
#include <sys/types.h>
#include <stdlib.h>

#include "cdi_cksum.h"
#include "cksum.h"
#include "error.h"
#include "serialize.h"

uint32_t cdiCheckSum(int type, int count, const void *buffer)
{
  uint32_t s = 0U;
  xassert(count >= 0);
  size_t elemSize = (size_t)serializeGetSizeInCore(1, type, NULL);
  memcrc_r_eswap(&s, (const unsigned char *)buffer, (size_t)count, elemSize);
  s = memcrc_finish(&s, (off_t)(elemSize * (size_t)count));
  return s;
}
