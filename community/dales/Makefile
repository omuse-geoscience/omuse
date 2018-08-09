# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR ?= ../../../..
-include $(AMUSE_DIR)/config.mk

MPIFC ?= mpif90

#MPIFC = ifort -I/cm/shared/apps/mpich/ge/open64/64/3.1/include -L/cm/shared/apps/mpich/ge/open64/64/3.1/lib -lmpich 

#LDFLAGS  += -lm $(MUSE_LD_FLAGS)
# MUSE is probably a typo - guess its not needed
# try without lm, ifort warns it  may supersede the Intel(R) math library and affect performance

DALESDIR = dales-repo
BUILDDIR = $(CURDIR)/dales-build
INCLUDEDIR = $(BUILDDIR)/src

CODE_LIB = $(BUILDDIR)/libdales.a

#DALESFLAGS := $(POPFFLAGS) $(NETCDF_FLAGS) $(NETCDFF_FLAGS) $(NETCDF_LIBS) $(NETCDFF_LIBS)
FS_LIBS += $(NETCDF_LIBS) $(NETCDFF_LIBS) $(SC_MPI_FCLIBS)

CODE_GENERATOR = $(AMUSE_DIR)/build.py


# DALES_FCFLAGS are extra compilation flags specific to Dales. Currently used for the size of 'real'
# FC is set in the Amuse configuration.  Amuse supplies -O2

# on ECMWF Cray, all compilers are invoked as ftn, so this test does not work. Specify  DALES_FCFLAGS outside this script on Cray.

ifeq ($(FC), ifort)
  DALES_FCFLAGS ?= -r8 -xHost -fpp
  SYST ?= localpc_ifort
  export SYST
endif

ifeq ($(FC), gfortran)
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

# DOWNLOAD_CODES=all will checkout entire repo with ssh, intended for developers of the components.
# DOWNLOAD_CODES=latest will (shallow) checkout latest revision only
# DOWNLOAD_CODES=<anythin else> will (shallow) checkout release tag spifs_v1.0.0
GIT_DEFAULT_TAG = spifs_v1.0.0
GIT_CLONE_CMD = git clone -b $(GIT_DEFAULT_TAG) --single-branch --depth=1
GIT_REPO = https://github.com/CloudResolvingClimateModeling/dales.git
ifeq ($(DOWNLOAD_CODES), all)
    GIT_CLONE_CMD = git clone
    GIT_REPO = git@github.com:CloudResolvingClimateModeling/dales.git
endif
ifeq ($(DOWNLOAD_CODES), latest)
    GIT_CLONE_CMD = git clone --depth=1
endif

all: dales dales_worker test/fortran_test

update: $(DALESDIR) 
	cd $(DALESDIR); git pull

dales: $(DALESDIR)
	mkdir -p $(BUILDDIR)
	cd $(BUILDDIR); cmake ../$(DALESDIR) -DCMAKE_BUILD_TYPE=$(CMAKE_BLD); make; make install

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py DalesInterface -o $@

dales_worker: worker_code.f90 $(OBJS) interface.f90 dales
	$(MPIFC) $(FCFLAGS) $(DALES_FCFLAGS) -I$(INCLUDEDIR) -c interface.f90 -o interface.o
	$(MPIFC) $(FCFLAGS) $(DALES_FCFLAGS) -I$(INCLUDEDIR) $(SC_FLAGS) $(FS_FLAGS) $< interface.o $(CODE_LIB) $(FS_LIBS) -o $@

test/fortran_test: test/fortran_test.f90 $(OBJS) interface.f90 dales
	$(MPIFC) $(FCFLAGS) $(DALES_FCFLAGS) -I$(INCLUDEDIR) -c interface.f90 -o interface.o
	$(MPIFC) $(FCFLAGS) $(DALES_FCFLAGS) -I$(INCLUDEDIR) $(SC_FLAGS) $(FS_FLAGS) $< interface.o $(CODE_LIB) $(FS_LIBS) -o $@

clean:
	$(RM) -f *.so *.o *.pyc *.mod
	$(RM) -f test/*.so test/*.o test/*.pyc test/*.mod
	$(RM) -rf $(BUILDDIR)
	$(RM) -f worker_code.* dales_worker test/fortran_test

distclean: clean

download:
	$(GIT_CLONE_CMD) $(GIT_REPO) $(DALESDIR)

ifdef DOWNLOAD_CODES
$(DALESDIR):
	make -C . download
else
$(DALESDIR):
	@echo "ERROR: need to download Dales, but DOWNLOAD_CODES is not set" 
	@false 
endif
