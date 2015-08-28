#ifndef _STREAM_CGRIBEX_H
#define _STREAM_CGRIBEX_H

int cgribexScanTimestep1(stream_t * streamptr);
int cgribexScanTimestep2(stream_t * streamptr);
int cgribexScanTimestep(stream_t * streamptr);

int cgribexDecode(unsigned char *gribbuffer, int gribsize, double *data, int gridsize,
		  int unreduced, int *nmiss, double missval);

size_t cgribexEncode(int memtype, int varID, int levelID, int vlistID, int gridID, int zaxisID,
		     int vdate, int vtime, int tsteptype, int numavg, 
		     long datasize, const double *data, int nmiss, unsigned char *gribbuffer, size_t gribbuffersize);

#endif  /* _STREAM_CGRIBEX_H */
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
