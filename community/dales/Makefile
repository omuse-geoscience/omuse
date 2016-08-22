# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR ?= ../../../..
-include $(AMUSE_DIR)/config.mk

MPIFC ?= mpif90

#MPIFC = ifort -I/cm/shared/apps/mpich/ge/open64/64/3.1/include -L/cm/shared/apps/mpich/ge/open64/64/3.1/lib -lmpich 

LDFLAGS  += -lm $(MUSE_LD_FLAGS)

SOURCEDIR = $(CURDIR)/dales-repo
BUILDDIR = $(CURDIR)/dales-build
INCLUDEDIR = $(BUILDDIR)/src

CODE_LIB = $(BUILDDIR)/libdales.a

#DALESFLAGS := $(POPFFLAGS) $(NETCDF_FLAGS) $(NETCDFF_FLAGS) $(NETCDF_LIBS) $(NETCDFF_LIBS)
FS_LIBS += $(NETCDF_LIBS) $(NETCDFF_LIBS) $(SC_MPI_FCLIBS)

CODE_GENERATOR = $(AMUSE_DIR)/build.py

all: dales dales_worker test/fortran_test

update:
	cd dales-repo; git pull

dales:
	mkdir -p $(BUILDDIR)
#	cd $(BUILDDIR); cmake ../dales-repo; make; make install
	cd $(BUILDDIR); cmake ../dales-repo -DCMAKE_BUILD_TYPE=Debug; make; make install

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py DalesInterface -o $@

dales_worker: worker_code.f90 $(OBJS) interface.f90
	$(MPIFC) $(FCFLAGS) -I$(INCLUDEDIR) -c interface.f90 -o interface.o
	$(MPIFC) $(FCFLAGS) -I$(INCLUDEDIR) $(SC_FLAGS) $(FS_FLAGS) $< interface.o $(CODE_LIB) $(FS_LIBS) -o $@

test/fortran_test: test/fortran_test.f90 $(OBJS) interface.f90
	$(MPIFC) $(FCFLAGS) -I$(INCLUDEDIR) -c interface.f90 -o interface.o
	$(MPIFC) $(FCFLAGS) -I$(INCLUDEDIR) $(SC_FLAGS) $(FS_FLAGS) $< interface.o $(CODE_LIB) $(FS_LIBS) -o $@

clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h *.mod
	$(RM) -rf $(BUILDDIR)

distclean: clean

