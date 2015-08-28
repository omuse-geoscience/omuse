#ifndef _SERVICE_H
#define _SERVICE_H


typedef struct {
  int    checked;
  int    byteswap;
  int    header[8];
  int    hprec;      /* header precision */
  int    dprec;      /* data   precision */
  size_t datasize;
  size_t buffersize;
  void  *buffer;
}
srvrec_t;


const char *srvLibraryVersion(void);

void srvDebug(int debug);

int  srvCheckFiletype(int fileID, int *swap);

srvrec_t *srvNew(void);
void srvDelete(srvrec_t *srvp);

int  srvRead(int fileID, srvrec_t *srvp);
void srvWrite(int fileID, srvrec_t *srvp);

int  srvInqHeader(srvrec_t *srvp, int *header);
int  srvInqDataSP(srvrec_t *srvp, float *data);
int  srvInqDataDP(srvrec_t *srvp, double *data);

int  srvDefHeader(srvrec_t *srvp, const int *header);
int  srvDefDataSP(srvrec_t *srvp, const float *data);
int  srvDefDataDP(srvrec_t *srvp, const double *data);


#endif  /* _SERVICE_H */
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
