# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include $(AMUSE_DIR)/config.mk

MPIFC ?= mpif90

FCFLAGS += -DBLAS -DLINUX -DINTEGER_IS_INT -DF90 -DPARAL -DNONCRAYF -m64 -fconvert=big-endian

#LDFLAGS  += -lm $(MUSE_LD_FLAGS)

# change this (not yet in amuse config.mk)
OIFS_GRIB_API_DIR  ?= /usr/local/
export OIFS_GRIB_API_DIR 
export OIFS_GRIB_API_LIB ?= -L$(OIFS_GRIB_API_DIR)/lib -lgrib_api_f90 -lgrib_api
#export OIFS_LAPACK_LIB ?= -L/usr/lib -llapack -lblas # the -L path causes problems. 
export OIFS_LAPACK_LIB ?= -llapack -lblas

FS_LIBS +=  $(OIFS_GRIB_API_LIB)  $(OIFS_LAPACK_LIB)

OIFSDIR = oifs-repo

CODE_GENERATOR = $(AMUSE_DIR)/build.py

# choose compiler - if not specified already using environment variables
CPPFLAGS = -cpp
ifeq ($(FC),ifort)
export OIFS_COMP ?= intel
FCFLAGS += -traceback
endif
ifeq ($(FC),gfortran)
export OIFS_COMP ?= gnu
CPPFLAGS = -x f95-cpp-input
endif

DEBUG ?= 0
export OIFS_BUILD ?= nomp
ifeq ($(DEBUG), 1)
    export OIFS_BUILD ?= noopt
endif
ifeq ($(DEBUG), yes)
    export OIFS_BUILD ?= noopt
endif

GIT_DEFAULT_TAG = spifs_v1.0.0
GIT_CLONE_CMD = git clone -b $(GIT_DEFAULT_TAG) --single-branch --depth=1
ifeq ($(DOWNLOAD_CODES), all)
    GIT_CLONE_CMD = git clone
endif
ifeq ($(DOWNLOAD_CODES), latest)
    GIT_CLONE_CMD = git clone --depth=1
endif

all: oifslib openifs_worker

update: $(OIFSDIR)
	cd $(OIFSDIR); git pull

oifslib: $(OIFSDIR)
	./$(OIFSDIR)/fcm/bin/fcm make -v -j4 -f $(OIFSDIR)/make/oifslib.cfg
	mkdir -p $(OIFSDIR)/make/$(OIFS_COMP)-$(OIFS_BUILD)/oifs/lib
	ar rcs $(OIFSDIR)/make/$(OIFS_COMP)-$(OIFS_BUILD)/oifs/lib/liboifs.a $(OIFSDIR)/make/$(OIFS_COMP)-$(OIFS_BUILD)/oifs/o/*.o 

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py OpenIFSInterface -o $@

openifs_worker: worker_code.f90 $(OBJS) interface.f90
	$(MPIFC) $(FCFLAGS) $(CPPFLAGS) -O2 -g -I$(OIFSDIR)/make/$(OIFS_COMP)-$(OIFS_BUILD)/oifs/include -c interface.f90 -o interface.o
	$(MPIFC) $(FCFLAGS) -O2 -g -I$(OIFSDIR)/make/$(OIFS_COMP)-$(OIFS_BUILD)/oifs/include $(SC_FLAGS) $(FS_FLAGS) $< interface.o $(OIFSDIR)/make/$(OIFS_COMP)-$(OIFS_BUILD)/oifs/lib/liboifs.a $(FS_LIBS) $(SC_FCLIBS) -o $@

clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h *.mod
	$(RM) *~ worker_code worker_code.f90 
	rm -f openifs_worker
	rm $(OIFSDIR)/make/$(OIFS_COMP)-$(OIFS_BUILD) -rf

distclean: clean

download:
	$(GIT_CLONE_CMD) https://git.ecmwf.int//scm/~g.vandenoord_esciencecenter.nl/oifs40r1-lib.git $(OIFSDIR)	


ifdef DOWNLOAD_CODES
$(OIFSDIR):
	make -C . download
else
$(OIFSDIR):
	@echo "ERROR: need to download OpenIFS, but DOWNLOAD_CODES is not set" 
	@false 
endif
