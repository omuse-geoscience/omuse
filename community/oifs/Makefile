# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include $(AMUSE_DIR)/config.mk

MPIFC ?= mpif90

LDFLAGS  += -lm $(MUSE_LD_FLAGS)

OIFSDIR=oifs40r1

OIFSPATH=$(CURDIR)/$(OIFSDIR)

CODE_GENERATOR = $(AMUSE_DIR)/build.py
OIFS_COMP = gnu
OIFS_BUILD = opt
OIFS_GRIB_API_DIR = /usr/local/
OIFS_GRIB_API_LIB = -L/usr/local/lib -lgrib_api_f90 -lgrib_api
OIFS_LAPACK_LIB = -L/usr/lib -llapack -lblas

all: oifs

update: oifs
	cd $(OIFSDIR); git pull

#       For some reason ssh access to the oifs repo doesn't work for me...
#	git clone ssh://git@software.ecmwf.int:7999/~g.vandenoord_esciencecenter.nl/oifs40r1-lib.git $(OIFSDIR)
oifs:
	git clone https://software.ecmwf.int/stash/scm/~g.vandenoord_esciencecenter.nl/oifs40r1-lib.git $(OIFSDIR)
	./$(OIFSDIR)/fcm/bin/fcm make -v -j4 -f $(OIFSDIR)/make/oifs.cfg

#TODO: Clean oifs directory
clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h *.mod
	$(RM) *~ worker_code worker_code.f90 
	rm -f oifs_worker_*

distclean: clean
	rm -rf $(OIFSDIR)

$(CODELIB1): src_320x384x40/POP_DomainSizeMod.F90
	make -C $(POPEXEDIR1) POPARCH=$(POPARCH) POPEXEDIR=$(POPEXEDIR1) POPDIR=$(POPDIR) AMUSE_DIR="$(realpath $(AMUSE_DIR))" POPFFLAGS="$(POPFFLAGS)" MPIFC="$(MPIFC)"


