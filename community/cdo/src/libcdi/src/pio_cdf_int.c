#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#if defined (HAVE_NETCDF4) && defined (HAVE_PARALLEL_NC4)

#include <netcdf.h>
#ifdef HAVE_NETCDF_PAR_H
#include <netcdf_par.h>
#endif

#include "namespace.h"
#include "pio.h"
#include "cdipio.h"
#include "pio_comm.h"
#include "pio_cdf_int.h"
#include "pio_util.h"

#include "pio_cdf_int.h"

static int
cdiPio_nc__create(const char *path, int cmode,
                  size_t initialsz, size_t *chunksizehintp,
                  int *ncidp)
{
  int status;
  if (cmode & NC_NETCDF4
      && commInqIOMode() != PIO_NONE)
    {
      cmode |= NC_MPIPOSIX;
      status = nc_create_par(path, cmode, commInqCommColl(),
                             MPI_INFO_NULL, ncidp);
    }
  else if (cmode & (NC_64BIT_OFFSET | NC_CLASSIC_MODEL)
           && commInqIOMode() != PIO_NONE)
    {
      /* FIXME: improve handling of pnetcdf here */
      abort();
    }
  else
    status = nc__create(path, cmode, initialsz, chunksizehintp, ncidp);
  return status;
}

static void
cdiPioCdfDefVar(int ncid, const char *name, nc_type xtype, int ndims,
                const int dimids[], int *varidp)
{
  cdf_def_var_serial(ncid, name, xtype, ndims, dimids, varidp);
  if (commInqIOMode() != PIO_NONE)
    {
      xdebug("%s", "calling nc_var_par_access");
      int status = nc_var_par_access(ncid, *varidp, NC_COLLECTIVE);
      if ( status != NC_NOERR )
        Error("%s", nc_strerror(status));
    }
}

void
cdiPioEnableNetCDFParAccess()
{
  namespaceSwitchSet(NSSWITCH_NC__CREATE, NSSW_FUNC(cdiPio_nc__create));
  namespaceSwitchSet(NSSWITCH_CDF_DEF_VAR, NSSW_FUNC(cdiPioCdfDefVar));
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
