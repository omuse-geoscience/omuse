/*
 * An implementation of the iterator interface for GRIB files.
 * Since GRIB files do not contain an index, this avoids scanning the entire file to generate an in-memory index as streamOpenRead() does.
 * Consequently, using this interface is much more efficient for GRIB files.
 */

#ifndef INCLUDE_GUARD_CDI_ITERATOR_GRIB_H
#define INCLUDE_GUARD_CDI_ITERATOR_GRIB_H

#include "iterator.h"
#include "input_file.h"

#ifdef HAVE_LIBGRIB_API
#include <grib_api.h>
#endif

typedef struct recordList recordList;
struct CdiGribIterator {
  CdiIterator super;

  CdiInputFile* file;
  off_t fileOffset;
  unsigned char* gribBuffer;
  size_t bufferSize, curRecordSize;
#ifdef HAVE_LIBGRIB_API
  grib_handle* gribHandle;
#else
  void* gribHandle;
#endif
};

CdiIterator* cdiGribIterator_new(const char* path, int filetype);
CdiGribIterator* cdiGribIterator_makeClone(CdiIterator* me);
char* cdiGribIterator_serialize(CdiIterator* me);
CdiGribIterator* cdiGribIterator_deserialize(const char* me);

int cdiGribIterator_nextField(CdiIterator* me);

char* cdiGribIterator_inqTime(CdiIterator* me, bool getEndTime);
int cdiGribIterator_levelType(CdiIterator* me, int levelSelector, char** outName, char** outLongName, char** outStdName, char** outUnit);
int cdiGribIterator_level(CdiIterator* me, int levelSelector, double* outValue1, double* outValue2);
int cdiGribIterator_zaxisUuid(CdiIterator* me, int* outVgridNumber, int* outLevelCount, unsigned char (*outUuid)[16]);
char* cdiGribIterator_copyVariableName(CdiIterator* me);

void cdiGribIterator_readField(CdiIterator* me, double* buffer, size_t* nmiss);
void cdiGribIterator_readFieldF(CdiIterator* me, float* buffer, size_t* nmiss);

#endif
