ifeq ($(origin AMUSE_DIR), undefined)
  AMUSE_DIR := $(shell amusifier --get-amuse-dir)
endif
-include ${AMUSE_DIR}/config.mk

MPIFC ?= mpif90

#MPIFC = ifort -I/cm/shared/apps/mpich/ge/open64/64/3.1/include -L/cm/shared/apps/mpich/ge/open64/64/3.1/lib -lmpich 

# NOTE on DALES versions
# OMUSE previously relied on a fork of DALES containing additions necessary for interfacing with OMUSE.
# These additions are part of the official DALES from version 4.3 on.
# 
# When using the previous fork of DALES, the following settings apply:
# INCLUDEDIR = $(BUILDDIR)/src/ 
#
# DALES_GIT_TAG ?= master
# DALES_GIT_REPO ?= https://github.com/CloudResolvingClimateModeling/dales.git
#
# Also, when using the older fork, replace "make -j 4" by "make" in the "dales:" rule below,
# to avoid a parallel build issue which was fixed in version 4.3 

DALESDIR = dales-repo
BUILDDIR = $(CURDIR)/dales-build
INCLUDEDIR = $(BUILDDIR)/src/library_modules/

CODE_LIB = $(BUILDDIR)/libdales.a

#DALESFLAGS := $(POPFFLAGS) $(NETCDF_FLAGS) $(NETCDFF_FLAGS) $(NETCDF_LIBS) $(NETCDFF_LIBS)
FS_LIBS += $(NETCDF_LIBS) $(NETCDFF_LIBS) $(SC_MPI_FCLIBS)

CODE_GENERATOR = $(AMUSE_DIR)/build.py


# DALES_FCFLAGS are extra compilation flags specific to Dales. Currently used for the size of 'real'
# FC is set in the Amuse configuration.  Amuse supplies -O2

# on ECMWF Cray, all compilers are invoked as ftn, so this test does not work. Specify  DALES_FCFLAGS outside this script on Cray.

ifeq ($(findstring ifort, $(notdir $(FC))), ifort)
# using ifort
  DALES_FCFLAGS ?= -r8 -xHost -fpp
  SYST ?= localpc_ifort
  export SYST
endif

ifeq ($(findstring gfortran, $(notdir $(FC))), gfortran)
# using gfortran
  DALES_FCFLAGS ?= -fdefault-real-8 -cpp
  SYST ?= gfortran
  export SYST
endif

DEBUG ?= 0
CMAKE_BLD = Release
ifeq ($(DEBUG), 1)
    CMAKE_BLD = Debug
endif
ifeq ($(DEBUG), yes)
    CMAKE_BLD = Debug
endif

# DOWNLOAD_CODES=all will clone the entire repo and check out $DALES_GIT_TAG, intended for developers of the components.
# DOWNLOAD_CODES=<anything else> will (shallow) checkout $DALES_GIT_TAG.
#   $DALES_GIT_TAG defaults to the release tag "4.3-rc.1"

DALES_GIT_TAG ?= 4.3-rc.1
DALES_GIT_REPO ?= https://github.com/dalesteam/dales.git
GIT_CLONE_CMD = git clone -b $(DALES_GIT_TAG) --single-branch --depth=1

ifeq ($(DOWNLOAD_CODES), all)
    GIT_CLONE_CMD = git clone -b $(DALES_GIT_TAG)
endif

all: dales dales_worker test/run_dales data

update: $(DALESDIR) 
	cd $(DALESDIR); git pull

dales: $(DALESDIR)
	mkdir -p $(BUILDDIR)
	cd $(BUILDDIR); cmake ../$(DALESDIR) -DCMAKE_BUILD_TYPE=$(CMAKE_BLD); make -j 4; make install

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py DalesInterface -o $@

dales_worker: worker_code.f90 $(OBJS) interface.f90 dales
	$(MPIFC) $(FCFLAGS) $(DALES_FCFLAGS) -I$(INCLUDEDIR) -c interface.f90 -o interface.o
	$(MPIFC) $(FCFLAGS) $(DALES_FCFLAGS) -I$(INCLUDEDIR) $(SC_FLAGS) $(FS_FLAGS) $< interface.o $(CODE_LIB) $(FS_LIBS) -o $@

test/run_dales: test/run_dales.f90 $(OBJS) interface.f90 dales
	$(MPIFC) $(FCFLAGS) $(DALES_FCFLAGS) -I$(INCLUDEDIR) -c interface.f90 -o interface.o
	$(MPIFC) $(FCFLAGS) $(DALES_FCFLAGS) -I$(INCLUDEDIR) $(SC_FLAGS) $(FS_FLAGS) $< interface.o $(CODE_LIB) $(FS_LIBS) -o $@

clean:
	$(RM) -f *.so *.o *.pyc *.mod
	$(RM) -f test/*.so test/*.o test/*.pyc test/*.mod
	$(RM) -rf $(BUILDDIR)
	$(RM) -f worker_code.* dales_worker test/fortran_test

distclean: clean
	$(RM) -rf $(DALESDIR)

download:
	$(GIT_CLONE_CMD) $(DALES_GIT_REPO) $(DALESDIR)

ifdef DOWNLOAD_CODES
$(DALESDIR):
	make -C . download
else
$(DALESDIR):
	@echo "ERROR: need to download Dales, but DOWNLOAD_CODES is not set" 
	@false 
endif

data:
	mkdir -p data
	cp -r dales-repo/cases data/
