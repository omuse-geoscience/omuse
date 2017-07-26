# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include $(AMUSE_DIR)/config.mk

MPIFC ?= mpif90

FCFLAGS += -fopenmp -DBLAS -DLINUX -DINTEGER_IS_INT -DF90 -DPARAL -DNONCRAYF -m64 -fconvert=big-endian

LDFLAGS  += -lm $(MUSE_LD_FLAGS)

OIFS_GRIB_API_DIR  ?= /usr/local/
export OIFS_GRIB_API_DIR 
export OIFS_GRIB_API_LIB ?= -L$(OIFS_GRIB_API_DIR)/lib -lgrib_api_f90 -lgrib_api
export OIFS_LAPACK_LIB ?= -L/usr/lib -llapack -lblas

FS_LIBS +=  $(OIFS_GRIB_API_LIB)  $(OIFS_LAPACK_LIB)

OIFSDIR=oifslib

CODE_GENERATOR = $(AMUSE_DIR)/build.py

# choose compiler - if not specified already using environment variables
ifeq ($(FC),ifort)
export OIFS_COMP ?= intel
else
export OIFS_COMP ?= gnu
endif

export OIFS_BUILD ?= noopt



all: oifslib oifssrc openifs_worker

update: oifslib
	cd $(OIFSDIR); git pull

oifslib:
	git clone -b oifslib https://software.ecmwf.int/stash/scm/~g.vandenoord_esciencecenter.nl/oifs40r1-lib.git $(OIFSDIR)

oifssrc:
	./$(OIFSDIR)/fcm/bin/fcm make -v -j4 -f $(OIFSDIR)/make/oifslib.cfg
	mkdir -p $(OIFSDIR)/make/$(OIFS_COMP)-$(OIFS_BUILD)/oifs/lib
	ar rcs $(OIFSDIR)/make/$(OIFS_COMP)-$(OIFS_BUILD)/oifs/lib/liboifs.a $(OIFSDIR)/make/$(OIFS_COMP)-$(OIFS_BUILD)/oifs/o/*.o 

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py OpenIFSInterface -o $@

openifs_worker: worker_code.f90 $(OBJS) interface.f90
	$(MPIFC) $(FCFLAGS) -O0 -g -I$(OIFSDIR)/make/$(OIFS_COMP)-$(OIFS_BUILD)/oifs/include -c interface.f90 -o interface.o
	$(MPIFC) $(FCFLAGS) -O0 -g -I$(OIFSDIR)/make/$(OIFS_COMP)-$(OIFS_BUILD)/oifs/include $(SC_FLAGS) $(FS_FLAGS) $< interface.o $(OIFSDIR)/make/$(OIFS_COMP)-$(OIFS_BUILD)/oifs/lib/liboifs.a $(FS_LIBS) $(SC_FCLIBS) -o $@

#TODO: Clean oifs directory
clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h *.mod
	$(RM) *~ worker_code worker_code.f90 
	rm -f oifs_worker_*
	rm oifslib/make/$(OIFS_COMP)-$(OIFS_BUILD) -rf

distclean: clean
	rm -rf $(OIFSDIR)


