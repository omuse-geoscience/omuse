#ifndef CDI_CKSUM_H_
#define CDI_CKSUM_H_

#include <inttypes.h>

uint32_t cdiCheckSum(int type, int count, const void *data);

#endif
