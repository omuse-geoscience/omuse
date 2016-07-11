# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR ?= ../../../..
-include $(AMUSE_DIR)/config.mk

MPIFC ?= mpif90

#MPIFC = ifort -I/cm/shared/apps/mpich/ge/open64/64/3.1/include -L/cm/shared/apps/mpich/ge/open64/64/3.1/lib -lmpich 

LDFLAGS  += -lm $(MUSE_LD_FLAGS)

SOURCEDIR = $(CURDIR)/dales-repo
BUILDDIR = $(CURDIR)/dales-build

#DALESFLAGS := $(POPFFLAGS) $(NETCDF_FLAGS) $(NETCDFF_FLAGS) $(NETCDF_LIBS) $(NETCDFF_LIBS)
#FS_FLAGS += $(NETCDF_LIBS) $(NETCDFF_LIBS)

CODE_GENERATOR = $(AMUSE_DIR)/build.py

all: dales

update:
	cd dales-repo; git pull

dales:
	mkdir -p $(BUILDDIR)
	cd $(BUILDDIR); cmake -DCMAKE_Fortran_FLAGS=$(DALESFLAGS) ../dales-repo; make

clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h *.mod
	$(RM) -rf $(BUILDDIR)

distclean: clean

