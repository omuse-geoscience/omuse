      module amuse_qgcm
      
      USE parameters
      USE atconst
      USE occonst
      USE radiate
      USE timinfo
      USE monitor
#ifndef ocean_only
      USE athomog
      USE amlsubs
      USE qgasubs
      USE atisubs
#endif
#if defined tau_udiff || !defined ocean_only
      USE atstate
#endif
#ifndef atmos_only
      USE ochomog
      USE ocstate
      USE omlsubs
      USE qgosubs
      USE ocisubs
#endif
      USE intrfac
      USE conhoms
      USE vorsubs
      USE xfosubs, ONLY : xforc
      USE eigmode
      USE valsubs
      USE radsubs
      USE topsubs
#ifdef use_netcdf
      USE nc_subs
#endif /* use_netcdf */
*     Diagnostics modules
      USE mondiag
#ifdef get_areav
      USE areasubs
#endif
#ifndef atmos_only
#  ifdef qoc_diag
      USE qocdiag
#  endif /* qoc_diag */
#endif /* atmos_only */
#ifdef get_covar
      USE covaria
#endif /* get_covar */
      USE timavge

      contains



      end module
