#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>

#include "cdi.h"
#include "cdi_int.h"
#include "cdf_int.h"
#include "namespace.h"

extern int CDF_Debug;

#if  defined  (HAVE_LIBNETCDF)
/*
#if ! defined (MIN_BUF_SIZE)
#  define  MIN_BUF_SIZE  131072L
#endif

static size_t ChunkSizeMin = MIN_BUF_SIZE;
*/
void cdf_create(const char *path, int cmode, int *ncidp)
{
  int status;
  int oldfill;
  size_t initialsz = 0, chunksizehint = 0;
  /*
#if defined (HAVE_STRUCT_STAT_ST_BLKSIZE)
  struct stat filestat;
  char basename[1024];
  char *pend;

  pend = strrchr(path, '/');
  if ( pend == 0 )
    strcpy(basename, "./");
  else
    {
      memcpy(basename, path, pend-path);
      basename[pend-path] = 0;
    }

  if ( stat(basename, &filestat) != 0 )
    SysError(basename);

  chunksizehint = (size_t) filestat.st_blksize * 4;
#endif

  if ( chunksizehint < ChunkSizeMin ) chunksizehint = ChunkSizeMin;
  */
#if defined(__SX__) || defined(ES)
  chunksizehint = 16777216; /* 16 MB */
#endif

  if ( cdiNcChunksizehint != CDI_UNDEFID )
    chunksizehint = (size_t)cdiNcChunksizehint;

  cdi_nc__create_funcp my_nc__create =
    (cdi_nc__create_funcp)namespaceSwitchGet(NSSWITCH_NC__CREATE).func;
  status = my_nc__create(path, cmode, initialsz, &chunksizehint, ncidp);

  if ( CDF_Debug || status != NC_NOERR )
    Message("ncid = %d  mode = %d  file = %s", *ncidp, cmode, path);

  if ( CDF_Debug || status != NC_NOERR )
    Message("chunksizehint %d", chunksizehint);

  if ( status != NC_NOERR ) Error("%s: %s", path, nc_strerror(status));

  status = nc_set_fill(*ncidp, NC_NOFILL, &oldfill);

  if ( status != NC_NOERR ) Error("%s: %s", path, nc_strerror(status));
}


int cdf_open(const char *path, int omode, int *ncidp)
{
  int status = 0;
  int dapfile = FALSE;
  struct stat filestat;
  size_t chunksizehint = 0;

#if  defined  (HAVE_LIBNC_DAP)
  if ( strncmp(path, "http:", 5) == 0 || strncmp(path, "https:", 6) == 0 ) dapfile = TRUE;
#endif

  if ( dapfile )
    {
      status = nc_open(path, omode, ncidp);
    }
  else
    {
      if ( stat(path, &filestat) != 0 ) SysError(path);

#if defined (HAVE_STRUCT_STAT_ST_BLKSIZE)
      chunksizehint = (size_t) filestat.st_blksize * 4;
#endif
      /*
      if ( chunksizehint < ChunkSizeMin ) chunksizehint = ChunkSizeMin;
      */
      if ( cdiNcChunksizehint != CDI_UNDEFID )
        chunksizehint = (size_t)cdiNcChunksizehint;

      /* FIXME: parallel part missing */
      status = nc__open(path, omode, &chunksizehint, ncidp);

      if ( CDF_Debug ) Message("chunksizehint %d", chunksizehint);
    }

  if ( CDF_Debug )
    Message("ncid = %d  mode = %d  file = %s", *ncidp, omode, path);

  if ( CDF_Debug && status != NC_NOERR ) Message("%s", nc_strerror(status));

  return (status);
}


void cdf_close(int ncid)
{
  int status;

  status = nc_close(ncid);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void cdf_redef(int ncid)
{
  int status;

  status = nc_redef(ncid);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void cdf_enddef(int ncid)
{
  int status;

  status = nc_enddef(ncid);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void cdf__enddef(const int ncid, const size_t hdr_pad)
{
  int status;
  const size_t v_align   = 4UL; /* [B] Alignment of beginning of data section for fixed variables */
  const size_t v_minfree = 0UL; /* [B] Pad at end of data section for fixed size variables */
  const size_t r_align   = 4UL; /* [B] Alignment of beginning of data section for record variables */

  /* nc_enddef(ncid) is equivalent to nc__enddef(ncid, 0, 4, 0, 4) */
  status = nc__enddef(ncid, hdr_pad, v_align, v_minfree, r_align);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void cdf_sync(int ncid)
{
  int status;

  status = nc_sync(ncid);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void cdf_inq(int ncid, int *ndimsp, int *nvarsp, int *ngattsp, int *unlimdimidp)
{
  int status;

  status = nc_inq(ncid, ndimsp, nvarsp, ngattsp, unlimdimidp);

  if ( CDF_Debug || status != NC_NOERR )
    Message("ncid = %d ndims = %d nvars = %d ngatts = %d unlimid = %d",
	    ncid, *ndimsp, *nvarsp, *ngattsp, *unlimdimidp);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void cdf_def_dim(int ncid, const char *name, size_t len, int *dimidp)
{
  int status;

  status = nc_def_dim(ncid, name, len, dimidp);

  if ( CDF_Debug || status != NC_NOERR )
    Message("ncid = %d  name = %s  len = %d", ncid, name, len);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void cdf_inq_dimid(int ncid, const char *name, int *dimidp)
{
  int status;

  status = nc_inq_dimid(ncid, name, dimidp);

  if ( CDF_Debug || status != NC_NOERR )
    Message("ncid = %d  name = %s  dimid= %d", ncid, name, *dimidp);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void cdf_inq_dim(int ncid, int dimid, char *name, size_t * lengthp)
{
  int status;

  status = nc_inq_dim(ncid, dimid, name, lengthp);

  if ( CDF_Debug || status != NC_NOERR )
    Message("ncid = %d  dimid = %d  length = %d name = %s", ncid, dimid, *lengthp, name);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void cdf_inq_dimname(int ncid, int dimid, char *name)
{
  int status;

  status = nc_inq_dimname(ncid, dimid, name);

  if ( CDF_Debug || status != NC_NOERR )
    Message("ncid = %d  dimid = %d  name = %s", ncid, dimid, name);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void cdf_inq_dimlen(int ncid, int dimid, size_t * lengthp)
{
  int status;

  status = nc_inq_dimlen(ncid, dimid, lengthp);

  if ( CDF_Debug || status != NC_NOERR )
    Message("ncid = %d dimid = %d length = %d", ncid, dimid, *lengthp);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void cdf_def_var(int ncid, const char *name, nc_type xtype, int ndims,
                 const int dimids[], int *varidp)
{
  cdi_cdf_def_var_funcp my_cdf_def_var
    = (cdi_cdf_def_var_funcp)namespaceSwitchGet(NSSWITCH_CDF_DEF_VAR).func;
  my_cdf_def_var(ncid, name, xtype, ndims, dimids, varidp);
}

void
cdf_def_var_serial(int ncid, const char *name, nc_type xtype, int ndims,
                   const int dimids[], int *varidp)
{
  int status = nc_def_var(ncid, name, xtype, ndims, dimids, varidp);

  if ( CDF_Debug || status != NC_NOERR )
    Message("ncid = %d  name = %s  xtype = %d  ndims = %d  varid = %d",
	    ncid, name, xtype, ndims, *varidp);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}



void cdf_inq_varid(int ncid, const char *name, int *varidp)
{
  int status;

  status = nc_inq_varid(ncid, name, varidp);

  if ( CDF_Debug || status != NC_NOERR )
    Message("ncid = %d  name = %s  varid = %d ", ncid, name, *varidp);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void cdf_inq_nvars(int ncid, int *nvarsp)
{
  int status;

  status = nc_inq_nvars(ncid, nvarsp);

  if ( CDF_Debug || status != NC_NOERR )
    Message("ncid = %d  nvars = %d", ncid, *nvarsp);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void cdf_inq_var(int ncid, int varid, char *name, nc_type *xtypep, int *ndimsp,
		 int dimids[], int *nattsp)
{
  int status;

  status = nc_inq_var(ncid, varid, name, xtypep, ndimsp, dimids, nattsp);

  if ( CDF_Debug || status != NC_NOERR )
    Message("ncid = %d varid = %d ndims = %d xtype = %d natts = %d name = %s",
	    ncid, varid, *ndimsp, *xtypep, *nattsp, name);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void cdf_inq_varname(int ncid, int varid, char *name)
{
  int status;

  status = nc_inq_varname(ncid, varid, name);

  if ( CDF_Debug || status != NC_NOERR )
    Message("ncid = %d varid = %d name = %s", ncid, varid, name);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void cdf_inq_vartype(int ncid, int varid, nc_type *xtypep)
{
  int status;

  status = nc_inq_vartype(ncid, varid, xtypep);

  if ( CDF_Debug || status != NC_NOERR )
    Message("ncid = %d varid = %d xtype = %s", ncid, varid, *xtypep);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void cdf_inq_varndims(int ncid, int varid, int *ndimsp)
{
  int status;

  status = nc_inq_varndims(ncid, varid, ndimsp);

  if ( CDF_Debug || status != NC_NOERR )
    Message("ncid = %d varid = %d", ncid, varid);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void cdf_inq_vardimid(int ncid, int varid, int dimids[])
{
  int status;

  status = nc_inq_vardimid(ncid, varid, dimids);

  if ( CDF_Debug || status != NC_NOERR )
    Message("ncid = %d varid = %d", ncid, varid);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void cdf_inq_varnatts(int ncid, int varid, int *nattsp)
{
  int status;

  status = nc_inq_varnatts(ncid, varid, nattsp);

  if ( CDF_Debug || status != NC_NOERR )
    Message("ncid = %d varid = %d nattsp = %d", ncid, varid, *nattsp);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void cdf_put_var_text(int ncid, int varid, const char *tp)
{
  int status;

  status = nc_put_var_text(ncid, varid, tp);

  if ( CDF_Debug || status != NC_NOERR )
    Message("%d %d %s", ncid, varid, tp);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void cdf_put_var_short(int ncid, int varid, const short *sp)
{
  int status;

  status = nc_put_var_short(ncid, varid, sp);

  if ( CDF_Debug || status != NC_NOERR )
    Message("%d %d %hd", ncid, varid, *sp);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void cdf_put_var_int(int ncid, int varid, const int *ip)
{
  int status;

  status = nc_put_var_int(ncid, varid, ip);

  if ( CDF_Debug || status != NC_NOERR )
    Message("%d %d %d", ncid, varid, *ip);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void cdf_put_var_long(int ncid, int varid, const long *lp)
{
  int status;

  status = nc_put_var_long(ncid, varid, lp);

  if ( CDF_Debug || status != NC_NOERR )
    Message("%d %d %ld", ncid, varid, *lp);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void cdf_put_var_float(int ncid, int varid, const float *fp)
{
  int status;

  status = nc_put_var_float(ncid, varid, fp);

  if ( CDF_Debug || status != NC_NOERR )
    Message("%d %d %f", ncid, varid, *fp);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void cdf_put_vara_double(int ncid, int varid, const size_t start[],
                         const size_t count[], const double *dp)
{
  int status;

  status = nc_put_vara_double(ncid, varid, start, count, dp);

  if ( CDF_Debug || status != NC_NOERR )
    Message("ncid = %d varid = %d val0 = %f", ncid, varid, *dp);

  if ( status != NC_NOERR )
    {
      char name[256];
      nc_inq_varname(ncid, varid, name);
      Message("varname = %s", name);
    }

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void  cdf_put_vara_float(int ncid, int varid, const size_t start[],
                         const size_t count[], const float *fp)
{
  int status;

  status = nc_put_vara_float(ncid, varid, start, count, fp);

  if ( CDF_Debug || status != NC_NOERR )
    Message("ncid = %d varid = %d val0 = %f", ncid, varid, *fp);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void  cdf_get_vara_int(int ncid, int varid, const size_t start[],
                       const size_t count[], int *dp)
{
  int status;

  status = nc_get_vara_int(ncid, varid, start, count, dp);

  if ( CDF_Debug || status != NC_NOERR )
    Message("ncid = %d varid = %d", ncid, varid);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void  cdf_get_vara_double(int ncid, int varid, const size_t start[],
                          const size_t count[], double *dp)
{
  int status;

  status = nc_get_vara_double(ncid, varid, start, count, dp);

  if ( CDF_Debug || status != NC_NOERR )
    Message("ncid = %d varid = %d", ncid, varid);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void  cdf_get_vara_float(int ncid, int varid, const size_t start[],
                         const size_t count[], float *fp)
{
  int status;

  status = nc_get_vara_float(ncid, varid, start, count, fp);

  if ( CDF_Debug || status != NC_NOERR )
    Message("ncid = %d varid = %d", ncid, varid);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void  cdf_get_vara_text(int ncid, int varid, const size_t start[],
			const size_t count[], char *tp)
{
  int status;

  status = nc_get_vara_text(ncid, varid, start, count, tp);

  if ( CDF_Debug || status != NC_NOERR )
    Message("ncid = %d varid = %d", ncid, varid);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void  cdf_get_vara_uchar(int ncid, int varid, const size_t start[], const size_t count[], unsigned char *tp)
{
  int status;

  status = nc_get_vara_uchar(ncid, varid, start, count, tp);

  if ( CDF_Debug || status != NC_NOERR )
    Message("ncid = %d varid = %d", ncid, varid);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void cdf_put_var_double(int ncid, int varid, const double *dp)
{
  int status;

  status = nc_put_var_double(ncid, varid, dp);

  if ( CDF_Debug || status != NC_NOERR )
    Message("ncid = %d varid = %d val0 = %f", ncid, varid, *dp);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void cdf_get_var1_text(int ncid, int varid, const size_t index[], char *tp)
{
  int status;

  status = nc_get_var1_text(ncid, varid, index, tp);

  if ( CDF_Debug || status != NC_NOERR )
    Message("ncid = %d varid = %d", ncid, varid);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void cdf_get_var1_double(int ncid, int varid, const size_t index[], double *dp)
{
  int status;

  status = nc_get_var1_double(ncid, varid, index, dp);

  if ( CDF_Debug || status != NC_NOERR )
    Message("ncid = %d varid = %d", ncid, varid);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void cdf_put_var1_double(int ncid, int varid, const size_t index[], const double *dp)
{
  int status;

  status = nc_put_var1_double(ncid, varid, index, dp);

  if ( CDF_Debug || status != NC_NOERR )
    Message("ncid = %d varid = %d val = %f", ncid, varid, *dp);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void cdf_get_var_text(int ncid, int varid, char *tp)
{
  int status;

  status = nc_get_var_text(ncid, varid, tp);

  if ( CDF_Debug || status != NC_NOERR )
    Message("ncid = %d varid = %d", ncid, varid);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void cdf_get_var_short(int ncid, int varid, short *sp)
{
  int status;

  status = nc_get_var_short(ncid, varid, sp);

  if ( CDF_Debug || status != NC_NOERR )
    Message("ncid = %d varid = %d", ncid, varid);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void cdf_get_var_int(int ncid, int varid, int *ip)
{
  int status;

  status = nc_get_var_int(ncid, varid, ip);

  if ( CDF_Debug || status != NC_NOERR )
    Message("ncid = %d varid = %d", ncid, varid);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void cdf_get_var_long(int ncid, int varid, long *lp)
{
  int status;

  status = nc_get_var_long(ncid, varid, lp);

  if ( CDF_Debug || status != NC_NOERR )
    Message("ncid = %d varid = %d", ncid, varid);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void cdf_get_var_float(int ncid, int varid, float *fp)
{
  int status;

  status = nc_get_var_float(ncid, varid, fp);

  if ( CDF_Debug || status != NC_NOERR )
    Message("ncid = %d varid = %d", ncid, varid);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void cdf_get_var_double(int ncid, int varid, double *dp)
{
  int status;

  status = nc_get_var_double(ncid, varid, dp);

  if ( CDF_Debug || status != NC_NOERR )
    Message("ncid = %d varid = %d val[0] = %f", ncid, varid, *dp);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void cdf_copy_att(int ncid_in, int varid_in, const char *name, int ncid_out,
		  int varid_out)
{
  int status;

  status = nc_copy_att(ncid_in, varid_in, name, ncid_out, varid_out);

  if ( CDF_Debug || status != NC_NOERR )
    Message("%d %d %s %d %d", ncid_in, varid_out, name, ncid_out, varid_out);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void cdf_put_att_text(int ncid, int varid, const char *name, size_t len,
		      const char *tp)
{
  int status;

  status = nc_put_att_text(ncid, varid, name, len, tp);

  if ( CDF_Debug || status != NC_NOERR )
    Message("ncid = %d varid = %d att = %s text = %.*s", ncid, varid, name, (int)len, tp);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void cdf_put_att_int(int ncid, int varid, const char *name, nc_type xtype,
		     size_t len, const int *ip)
{
  int status;

  status = nc_put_att_int(ncid, varid, name, xtype, len, ip);

  if ( CDF_Debug || status != NC_NOERR )
    Message("ncid = %d varid = %d att = %s val = %d", ncid, varid, name, *ip);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void cdf_put_att_double(int ncid, int varid, const char *name, nc_type xtype,
			size_t len, const double *dp)
{
  int status;

  status = nc_put_att_double(ncid, varid, name, xtype, len, dp);

  if ( CDF_Debug || status != NC_NOERR )
    Message("%d %d %f", ncid, varid, *dp);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void cdf_get_att_text(int ncid, int varid, const char *name, char *tp)
{
  int status;

  status = nc_get_att_text(ncid, varid, name, tp);

  if ( CDF_Debug || status != NC_NOERR )
    Message("ncid = %d varid = %d name = %s", ncid, varid, name);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}

void cdf_get_att_string(int ncid, int varid, const char *name, char **tp)
{
#if  defined  (HAVE_NETCDF4)
  int status;

  status = nc_get_att_string(ncid, varid, name, tp);

  if ( CDF_Debug || status != NC_NOERR )
    Message("ncid = %d varid = %d name = %s", ncid, varid, name);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
#endif
}


void cdf_get_att_int(int ncid, int varid, const char *name, int *ip)
{
  int status;

  status = nc_get_att_int(ncid, varid, name, ip);

  if ( CDF_Debug || status != NC_NOERR )
    Message("ncid = %d varid = %d att = %s val = %d", ncid, varid, name, *ip);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void cdf_get_att_double(int ncid, int varid, const char *name, double *dp)
{
  int status;

  status = nc_get_att_double(ncid, varid, name, dp);

  if ( CDF_Debug || status != NC_NOERR )
    Message("ncid = %d varid = %d att = %s val = %.9g",
	    ncid, varid, name, *dp);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void cdf_inq_att(int ncid, int varid, const char *name, nc_type *xtypep,
		 size_t *lenp)
{
  int status;

  status = nc_inq_att(ncid, varid, name, xtypep, lenp);

  if ( CDF_Debug || status != NC_NOERR )
    Message("ncid = %d varid = %d", ncid, varid);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void cdf_inq_atttype(int ncid, int varid, const char *name, nc_type * xtypep)
{
  int status;

  status = nc_inq_atttype(ncid, varid, name, xtypep);

  if ( CDF_Debug || status != NC_NOERR )
    Message("ncid = %d varid = %d", ncid, varid);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void cdf_inq_attlen(int ncid, int varid, const char *name, size_t * lenp)
{
  int status;

  status = nc_inq_attlen(ncid, varid, name, lenp);

  if ( CDF_Debug || status != NC_NOERR )
    Message("ncid = %d varid = %d name = %s len = %d", ncid, varid, name, *lenp);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void cdf_inq_attname(int ncid, int varid, int attnum, char *name)
{
  int status;

  status = nc_inq_attname(ncid, varid, attnum, name);

  if ( CDF_Debug || status != NC_NOERR )
    Message("ncid = %d varid = %d attnum = %d name = %s", ncid, varid, attnum, name);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}


void cdf_inq_attid(int ncid, int varid, const char *name, int *attnump)
{
  int status;

  status = nc_inq_attid(ncid, varid, name, attnump);

  if ( CDF_Debug || status != NC_NOERR )
    Message("ncid = %d varid = %d", ncid, varid);

  if ( status != NC_NOERR ) Error("%s", nc_strerror(status));
}

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
