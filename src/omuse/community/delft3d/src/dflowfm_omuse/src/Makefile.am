lib_LTLIBRARIES=libdflowfm_omuse.la
#AM_CXXFLAGS = $(INTI_CFLAGS)

# Import optimization and OPENMP (in AM_FCFLAGS) from common.am
include $(top_srcdir)/common.am

INCLUDES = \
    -I$(builddir)/../../dflowfm_kernel/include \
    -I$(builddir)/../../dflowfm_kernel/src


# these are generated using:
# cd $(srcdir)/src
# $(top_srcdir)/scripts/generate.py modules.f90 network_data.f90 --template-dir=templates

EXTRA_DIST = 

# default to empty so we can conditionally expand it
include_HEADERS=

EXTRALIBS=

# For source distributions: the tar.gz will contain the correct version numbers in the generated *version.F90
# So, just before packaging the .tar.gz dist, disable version number REgeneration parts of this Makefile.am.
dist-hook:
	sed -i -e 's/^REGENERATE_VN[ \t]*=[ \t]*yes/REGENERATE_VN=no/' $(distdir)/Makefile.am
	sed -i -e 's/^REGENERATE_VN[ \t]*=[ \t]*yes/REGENERATE_VN=no/' $(distdir)/Makefile.in

libdflowfm_omuse_la_SOURCES= dflowfm_omuse_lib.F90

# This line would prevent warnings from gfortran about assigned goto's (but introduce warning from ifort about unknown flag)
#saadf.o : FFLAGS += -std=legacy

# Version number generation (only when in SVN source dir, not for tar.gz dist)

.FORCE: 

.PHONY: .FORCE

MOSTLYCLEANFILES = *.mod


#-------------------------------------------------------------------------------
#   Fortran flags

FINCLUDES += \
	-I$(builddir)/../include
#~     -I$(builddir)/../../dflowfm_kernel/include \
#~     -I$(builddir)/../../dflowfm_kernel/src

FMODULES += \
    -I$(top_builddir)/utils_lgpl/ec_module/packages/ec_module/src \
    -I$(top_builddir)/utils_lgpl/io_netcdf/packages/io_netcdf/src \
    -I$(top_builddir)/utils_lgpl/gridgeom/packages/gridgeom/src \
    -I$(top_builddir)/utils_lgpl/deltares_common/packages/deltares_common/src \
    -I$(top_builddir)/utils_lgpl/deltares_common/packages/deltares_common_mpi/src \
    -I$(top_builddir)/utils_gpl/morphology/packages/morphology_data/src \
    -I$(top_builddir)/utils_gpl/morphology/packages/morphology_io/src \
    -I$(top_builddir)/utils_gpl/morphology/packages/morphology_kernel/src \
    -I$(top_builddir)/utils_gpl/trachytopes/packages/trachytopes_io/src \
    -I$(top_builddir)/utils_gpl/trachytopes/packages/trachytopes_kernel/src \
    -I$(top_builddir)/utils_gpl/flow1d/packages/flow1d_core/src \
    -I$(top_builddir)/utils_gpl/flow1d/packages/flow1d_io/src \
    -I$(top_builddir)/third_party_open/kdtree2/src-f90 \
    -I$(top_builddir)/third_party_open/libsigwatch \
    -I$(top_builddir)/utils_lgpl/kdtree_wrapper/packages/kdtree_wrapper/src \
    -I$(top_builddir)/engines_gpl/waq/include \
    -I$(top_builddir)/engines_gpl/waq/packages/wq_processes/src \
    -I$(top_builddir)/engines_gpl/waq/packages/waq_process/src \
    -I$(top_builddir)/engines_gpl/waq/packages/waq_process/src/bloom \
    -I$(top_builddir)/engines_gpl/waq/packages/waq_process/src/charon \
    -I$(top_builddir)/engines_gpl/waq/packages/waq_process/src/proc_preprocess \
    -I$(top_builddir)/engines_gpl/waq/packages/waq_process/src/waq_process \
    -I$(top_builddir)/engines_gpl/waq/packages/waq_utils_f/src/mod_timings \
    -I$(top_builddir)/engines_gpl/waq/packages/waq_utils_f/src/mod_couplib \
    -I$(top_builddir)/engines_gpl/waq/packages/waq_utils_f/src/waq_utils \
    -I$(top_builddir)/utils_lgpl/deltares_common/packages/deltares_common/src/rd_token \
    -I$(top_builddir)/utils_lgpl/deltares_common/packages/deltares_common/src

# NOTE:
# includes, modules, optimization and openmp are all added to AM_FCFLAGS by common.am

#-------------------------------------------------------------------------------
#   Fortran 90 flags

AM_FCFLAGS+= \
    -DHAVE_DISPLAY=0 \
    -DNetCDF4

if USENETCDF
AM_FCFLAGS+=${NETCDF_CFLAGS} ${NETCDF_FORTRAN_CFLAGS}
# Manualy add the netcdff (fortran dll), bug in pkg-config file is not fixed until in NetCDF 4.4 HEAD
LIBS+=${NETCDF_LIBS} ${NETCDF_FORTRAN_LIBS} -lnetcdff
endif

if USEMPI
AM_FCFLAGS+=${MPI_CFLAGS}
AM_FCFLAGS+=${MPI_FCFLAGS}
# Manualy add the MPI (fortran dll)
LIBS+=${MPILIBS}
endif

if USEPETSC
AM_FCFLAGS+=${PETSc_CFLAGS}
LIBS+=${PETSc_LIBS}
endif

if USEPROJ
AM_FCFLAGS+=${PROJ_CFLAGS}
LIBS+=${PROJ_LIBS}
USEFORTRANGIS=TRUE
endif

if USESHAPELIB
AM_FCFLAGS+=${SHAPELIB_CFLAGS}
LIBS+=${SHAPELIB_LIBS}
USEFORTRANGIS=TRUE
endif

if USEFORTRANGIS
FMODULES += \
    -I$(top_builddir)/third_party_open/fortrangis/libfortrangis
EXTRALIBS += \
    $(top_builddir)/third_party_open/fortrangis/libfortranc/libfortranc.la \
    $(top_builddir)/third_party_open/fortrangis/libfortrangis/libfortrangis.la
endif

if USEMETIS
AM_FCFLAGS+=${METIS_CFLAGS}
LIBS+=${METIS_LIBS}
endif

#if USETECIO
#AM_FCFLAGS+=${TECIO_FCFLAGS}
#LIBS+=${TECIO_LIBS}
#endif

if USEPARMS
LIBS+=-L../third_party_open/pARMS/pARMS_3.2/lib -lparms
LIBS+=-L${blasdir}/lib -lrefblas
endif

if USEMKL
LIBS+=${MKL_LIBS}
AM_FCFLAGS+=${MKL_FCFLAGS}
endif

# if USEESMF
# ESMFSOURCES=unstruc_esmf.F90
# else
# ESMFSOURCES=
# endif



EXTRALIBS+=\
        $(top_builddir)/engines_gpl/dflowfm/packages/dflowfm_kernel/src/libdflowfm.la \
        $(top_builddir)/utils_gpl/flow1d/packages/flow1d/src/libflow1d.la \
        $(top_builddir)/utils_gpl/morphology/packages/morphology_plugins_c/src/libMorphologyPluginsC.la \
        $(top_builddir)/utils_gpl/morphology/packages/morphology_kernel/src/libMorphologyKernel.la \
        $(top_builddir)/utils_gpl/morphology/packages/morphology_io/src/libMorphologyIO.la \
        $(top_builddir)/utils_gpl/morphology/packages/morphology_data/src/libMorphologyData.la \
        $(top_builddir)/utils_gpl/trachytopes/packages/trachytopes_kernel/src/libTrachytopesKernel.la \
        $(top_builddir)/utils_gpl/trachytopes/packages/trachytopes_io/src/libTrachytopesIO.la \
        $(top_builddir)/engines_gpl/waq/packages/waq_process/src/proc_preprocess/libproc_preprocess.la \
        $(top_builddir)/engines_gpl/waq/packages/waq_process/src/bloom/libbloom.la \
        $(top_builddir)/engines_gpl/waq/packages/waq_process/src/waq_process/libwaq_process.la \
        $(top_builddir)/engines_gpl/waq/packages/waq_plugin_wasteload/src/libwaq_plugin_wasteload.la \
        $(top_builddir)/engines_gpl/waq/packages/waq_utils_c/src/libwaq_utils_c.la \
        $(top_builddir)/engines_gpl/waq/packages/waq_utils_f/src/interact/libinteract.la \
        $(top_builddir)/engines_gpl/waq/packages/waq_utils_f/src/waq_utils/libwaq_utils.la \
        $(top_builddir)/engines_gpl/waq/packages/wq_processes/src/libwqprocesses.la \
        $(top_builddir)/utils_lgpl/nefis/packages/nefis/src/libNefis.la \
        $(top_builddir)/utils_lgpl/deltares_common/packages/deltares_common/src/libDeltaresCommon.la \
        $(top_builddir)/utils_lgpl/deltares_common/packages/deltares_common_mpi/src/libDeltaresCommonMpi.la \
        $(top_builddir)/utils_lgpl/kdtree_wrapper/packages/kdtree_wrapper/src/libkdtreewrapper.la \
        $(top_builddir)/utils_lgpl/ec_module/packages/ec_module/src/libECModule.la \
        $(top_builddir)/utils_lgpl/io_netcdf/packages/io_netcdf/src/libIONetcdf.la \
        $(top_builddir)/utils_lgpl/gridgeom/packages/gridgeom/src/libgridgeom.la \
        $(top_builddir)/utils_lgpl/metistools/packages/metisoptions/src/libmetisoptions.la \
        $(top_builddir)/third_party_open/triangle/libtriangle_c.la \
        $(top_builddir)/third_party/interacter_stub/src/libintstub.la \
        $(top_builddir)/third_party_open/kdtree2/src-f90/libkdtree2f.la \
        $(top_builddir)/third_party_open/polypack/src/libpolypack.la \
        $(top_builddir)/third_party_open/libsigwatch/libsigwatch.la \
        $(CCLIBS) -luuid

# dflowfm_LDADD = \
#         libdflowfm.la		#
dflowfm_omuse_LDFLAGS=-static
libdflowfm_omuse_la_LIBADD=\
        $(EXTRALIBS)
