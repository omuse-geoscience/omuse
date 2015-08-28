#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#ifndef SERIALIZE_H
#define SERIALIZE_H

#include <string.h>

#include "cdi.h"
#ifndef  CDI_CKSUM_H_
#include "cdi_cksum.h"
#endif
#ifndef  _ERROR_H
#include "error.h"
#endif

/*
 * Generic interfaces for (de-)marshalling
 */
int serializeGetSize(int count, int datatype, void *context);
void serializePack(const void *data, int count, int datatype,
                   void *buf, int buf_size, int *position, void *context);
void serializeUnpack(const void *buf, int buf_size, int *position,
                     void *data, int count, int datatype, void *context);

/*
 * (de-)marshalling function for common data structures
 */
static inline int
serializeStrTabGetPackSize(const char **strTab, int numStr,
                           void *context)
{
  xassert(numStr >= 0);
  int packBuffSize = 0;
  for (size_t i = 0; i < (size_t)numStr; ++i)
  {
    size_t len = strlen(strTab[i]);
    packBuffSize +=
      serializeGetSize(1, DATATYPE_INT, context)
      + serializeGetSize((int)len, DATATYPE_TXT, context);
  }
  packBuffSize +=
    serializeGetSize(1, DATATYPE_UINT32, context);
  return packBuffSize;
}

static inline void
serializeStrTabPack(const char **strTab, int numStr,
                    void *buf, int buf_size, int *position, void *context)
{
  uint32_t d = 0;
  xassert(numStr >= 0);
  for (size_t i = 0; i < (size_t)numStr; ++i)
  {
    size_t len = strlen(strTab[i]);
    serializePack(&(int){(int)len}, 1, DATATYPE_INT,
                  buf, buf_size, position, context);
    serializePack(strTab[i], (int)len, DATATYPE_TXT,
                  buf, buf_size, position, context);
    d ^= cdiCheckSum(DATATYPE_TXT, (int)len, strTab[i]);
  }
  serializePack(&d, 1, DATATYPE_UINT32,
                buf, buf_size, position, context);
}

static inline void
serializeStrTabUnpack(const void *buf, int buf_size, int *position,
                      char **strTab, int numStr, void *context)
{
  uint32_t d, d2 = 0;
  xassert(numStr >= 0);
  for (size_t i = 0; i < (size_t)numStr; ++i)
    {
      int len;
      serializeUnpack(buf, buf_size, position,
                      &len, 1, DATATYPE_INT, context);
      serializeUnpack(buf, buf_size, position,
                      strTab[i], len, DATATYPE_TXT, context);
      strTab[i][len] = '\0';
      d2 ^= cdiCheckSum(DATATYPE_TXT, (size_t)len, strTab[i]);
    }
  serializeUnpack(buf, buf_size, position,
                  &d, 1, DATATYPE_UINT32, context);
  xassert(d == d2);
}

/*
 * Interfaces for marshalling within a single memory domain
 */
int serializeGetSizeInCore(int count, int datatype, void *context);
void serializePackInCore(const void *data, int count, int datatype,
                          void *buf, int buf_size, int *position, void *context);
void serializeUnpackInCore(const void *buf, int buf_size, int *position,
                           void *data, int count, int datatype, void *context);

#endif
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
