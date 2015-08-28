#ifndef _EXTRA_H
#define _EXTRA_H

#define  EXT_REAL   1
#define  EXT_COMP   2


typedef struct {
  int    checked;
  int    byteswap;
  int    header[4];
  int    prec;      /* single or double precison */
  int    number;    /* real or complex */
  size_t datasize;
  size_t buffersize;
  void  *buffer;
}
extrec_t;


const char *extLibraryVersion(void);

void extDebug(int debug);

int  extCheckFiletype(int fileID, int *swap);

void *extNew(void);
void extDelete(void *ext);

int  extRead(int fileID, void *ext);
int  extWrite(int fileID, void *ext);

int  extInqHeader(void *ext, int *header);
int  extInqDataSP(void *ext, float *data);
int  extInqDataDP(void *ext, double *data);

int  extDefHeader(void *ext, const int *header);
int  extDefDataSP(void *ext, const float *data);
int  extDefDataDP(void *ext, const double *data);

#endif  /* _EXTRA_H */
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
