#ifndef _ZAXIS_H
#define _ZAXIS_H

void zaxisGetTypeDescription(int zaxisType, int* outPositive, const char** outName, const char** outLongName, const char** outStdName, const char** outUnit);  //The returned const char* point to static storage. Don't free or modify them.

unsigned cdiZaxisCount(void);

void cdiZaxisGetIndexList(unsigned numIDs, int IDs[numIDs]);

void
zaxisUnpack(char * unpackBuffer, int unpackBufferSize,
            int * unpackBufferPos, int originNamespace, void *context,
            int force_id);

void zaxisDefLtype2(int zaxisID, int ltype2);

extern const resOps zaxisOps;

#endif
