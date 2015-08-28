#ifndef _CDF_H
#define _CDF_H

void cdfDebug(int debug);

const char *cdfLibraryVersion(void);
const char *hdfLibraryVersion(void);

int  cdfOpen(const char *filename, const char *mode);
int  cdfOpen64(const char *filename, const char *mode);
int  cdf4Open(const char *filename, const char *mode, int *filetype);
void cdfClose(int fileID);

#endif  /* _CDF_H */
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
