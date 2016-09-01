# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include $(AMUSE_DIR)/config.mk

MPIFC ?= mpif90

LDFLAGS  += -lm $(MUSE_LD_FLAGS)

OIFSDIR=oifslib

OIFSPATH=$(CURDIR)/$(OIFSDIR)

CODE_GENERATOR = $(AMUSE_DIR)/build.py
	
export OIFS_COMP = gnu
export OIFS_BUILD = opt
export OIFS_GRIB_API_DIR = /usr/local/
export OIFS_GRIB_API_LIB = -L/usr/local/lib -lgrib_api_f90 -lgrib_api
export OIFS_LAPACK_LIB = -L/usr/lib -llapack -lblas

all: oifslib oifssrc

update: oifslib
	cd $(OIFSDIR); git pull

#       For some reason ssh access to the oifs repo doesn't work for me...
#	git clone ssh://git@software.ecmwf.int:7999/~g.vandenoord_esciencecenter.nl/oifs40r1-lib.git $(OIFSDIR)
oifslib:
	git clone -b feature/crm https://software.ecmwf.int/stash/scm/~g.vandenoord_esciencecenter.nl/oifs40r1-lib.git $(OIFSDIR)

oifssrc:
	./$(OIFSDIR)/fcm/bin/fcm make -v -j4 -f $(OIFSDIR)/make/oifs.cfg

#TODO: Clean oifs directory
clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h *.mod
	$(RM) *~ worker_code worker_code.f90 
	rm -f oifs_worker_*

distclean: clean
	rm -rf $(OIFSDIR)


