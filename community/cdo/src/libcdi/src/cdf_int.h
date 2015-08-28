#ifndef _CDF_INT_H
#define _CDF_INT_H

#if  defined  (HAVE_LIBNETCDF)

#include <stdlib.h>
#include <netcdf.h>

void cdf_create (const char *path, int cmode, int *idp);
int  cdf_open   (const char *path, int omode, int *idp);
void cdf_close  (int ncid);

void cdf_redef(int ncid);
void cdf_enddef(int ncid);
void cdf__enddef(const int ncid, const size_t hdr_pad);
void cdf_sync(int ncid);

void cdf_inq (int ncid, int *ndimsp, int *nvarsp, int *ngattsp, int *unlimdimidp);

void cdf_def_dim (int ncid, const char *name, size_t len, int *idp);
void cdf_inq_dimid (int ncid, const char *name, int *dimidp);
void cdf_inq_dim (int ncid, int dimid, char *name, size_t * lengthp);
void cdf_inq_dimname (int ncid, int dimid, char *name);
void cdf_inq_dimlen (int ncid, int dimid, size_t * lengthp);
void cdf_def_var (int ncid, const char *name, nc_type xtype, int ndims, const int dimids[], int *varidp);
void cdf_def_var_serial(int ncid, const char *name, nc_type xtype, int ndims, const int dimids[], int *varidp);
void cdf_inq_varid(int ncid, const char *name, int *varidp);
void cdf_inq_nvars(int ncid, int *nvarsp);
void cdf_inq_var(int ncid, int varid, char *name, nc_type *xtypep, int *ndimsp, int dimids[], int *nattsp);
void cdf_inq_varname (int ncid, int varid, char *name);
void cdf_inq_vartype (int ncid, int varid, nc_type *xtypep);
void cdf_inq_varndims (int ncid, int varid, int *ndimsp);
void cdf_inq_vardimid (int ncid, int varid, int dimids[]);
void cdf_inq_varnatts (int ncid, int varid, int *nattsp);

void cdf_copy_att (int ncid_in, int varid_in, const char *name, int ncid_out, int varid_out);
void cdf_put_var_text   (int ncid, int varid, const char *tp);
void cdf_put_var_uchar  (int ncid, int varid, const unsigned char *up);
void cdf_put_var_schar  (int ncid, int varid, const signed char *cp);
void cdf_put_var_short  (int ncid, int varid, const short *sp);
void cdf_put_var_int    (int ncid, int varid, const int *ip);
void cdf_put_var_long   (int ncid, int varid, const long *lp);
void cdf_put_var_float  (int ncid, int varid, const float *fp);
void cdf_put_var_double (int ncid, int varid, const double *dp);

void cdf_get_var_text   (int ncid, int varid, char *tp);
void cdf_get_var_uchar  (int ncid, int varid, unsigned char *up);
void cdf_get_var_schar  (int ncid, int varid, signed char *cp);
void cdf_get_var_short  (int ncid, int varid, short *sp);
void cdf_get_var_int    (int ncid, int varid, int *ip);
void cdf_get_var_long   (int ncid, int varid, long *lp);
void cdf_get_var_float  (int ncid, int varid, float *fp);
void cdf_get_var_double (int ncid, int varid, double *dp);

void cdf_get_var1_text(int ncid, int varid, const size_t index[], char *tp);

void cdf_get_var1_double(int ncid, int varid, const size_t index[], double *dp);
void cdf_put_var1_double(int ncid, int varid, const size_t index[], const double *dp);

void cdf_get_vara_uchar(int ncid, int varid, const size_t start[], const size_t count[], unsigned char *tp);
void cdf_get_vara_text(int ncid, int varid, const size_t start[], const size_t count[], char *tp);

void cdf_get_vara_double(int ncid, int varid, const size_t start[], const size_t count[], double *dp);
void cdf_put_vara_double(int ncid, int varid, const size_t start[], const size_t count[], const double *dp);

void cdf_get_vara_float(int ncid, int varid, const size_t start[], const size_t count[], float *fp);
void cdf_put_vara_float(int ncid, int varid, const size_t start[], const size_t count[], const float *fp);

void cdf_put_att_text(int ncid, int varid, const char *name, size_t len, const char *tp);
void cdf_put_att_int(int ncid, int varid, const char *name, nc_type xtype, size_t len, const int *ip);
void cdf_put_att_double(int ncid, int varid, const char *name, nc_type xtype, size_t len, const double *dp);

void cdf_get_att_string(int ncid, int varid, const char *name, char **tp);
void cdf_get_att_text  (int ncid, int varid, const char *name, char *tp);
void cdf_get_att_int   (int ncid, int varid, const char *name, int *ip);
void cdf_get_att_double(int ncid, int varid, const char *name, double *dp);

void cdf_inq_att    (int ncid, int varid, const char *name, nc_type * xtypep, size_t * lenp);
void cdf_inq_atttype(int ncid, int varid, const char *name, nc_type *xtypep);
void cdf_inq_attlen (int ncid, int varid, const char *name, size_t *lenp);
void cdf_inq_attname(int ncid, int varid, int attnum, char *name);
void cdf_inq_attid  (int ncid, int varid, const char *name, int *attnump);

typedef int (*cdi_nc__create_funcp)(const char *path, int cmode,
                                    size_t initialsz, size_t *chunksizehintp,
                                    int *ncidp);

typedef void (*cdi_cdf_def_var_funcp)(int ncid, const char *name,
                                      nc_type xtype, int ndims,
                                      const int dimids[], int *varidp);

#endif

#endif  /* _CDF_INT_H */
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
