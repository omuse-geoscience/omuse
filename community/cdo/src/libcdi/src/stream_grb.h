#ifndef _STREAM_GRB_H
#define _STREAM_GRB_H

int   grbBitsPerValue(int datatype);

int   grbInqContents(stream_t * streamptr);
int   grbInqTimestep(stream_t * streamptr, int tsID);

int   grbInqRecord(stream_t * streamptr, int *varID, int *levelID);
void  grbDefRecord(stream_t * streamptr);
void  grbReadRecord(stream_t * streamptr, double *data, int *nmiss);
void  grb_write_record(stream_t * streamptr, int memtype, const void *data, int nmiss);
void  grbCopyRecord(stream_t * streamptr2, stream_t * streamptr1);

void  grbReadVarDP(stream_t * streamptr, int varID, double *data, int *nmiss);
void  grb_write_var(stream_t * streamptr, int varID, int memtype, const void *data, int nmiss);

void  grbReadVarSliceDP(stream_t * streamptr, int varID, int levelID, double *data, int *nmiss);
void  grb_write_var_slice(stream_t *streamptr, int varID, int levelID, int memtype, const void *data, int nmiss);

int   grib1ltypeToZaxisType(int grib_ltype);
int   grib2ltypeToZaxisType(int grib_ltype);

int   zaxisTypeToGrib1ltype(int zaxistype);
int   zaxisTypeToGrib2ltype(int zaxistype);

#endif  /* _STREAM_GRB_H */
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
