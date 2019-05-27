# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=$(shell amusifier --get-amuse-dir)
-include $(AMUSE_DIR)/config.mk

MPIFC ?= mpif90

#MPIFC = ifort -I/cm/shared/apps/mpich/ge/open64/64/3.1/include -L/cm/shared/apps/mpich/ge/open64/64/3.1/lib -lmpich 

TARGETS := $(patsubst src_%, pop_worker_%, ${shell ls -d src_*})

ifeq ($(findstring gfortran, $(notdir $(FC))), gfortran)
FCFLAGS+=-fconvert=swap
POPFFLAGS=-O2 -Wall -fdefault-double-8 -fdefault-real-8 -fconvert=swap -fimplicit-none -fbounds-check
endif

ifeq ($(findstring ifort, $(notdir $(FC))), ifort)
FCFLAGS+=-O2 -r8 -g -convert big_endian -assume byterecl
POPFFLAGS=-O2 -r8 -g -convert big_endian -assume byterecl
endif

LDFLAGS  += -lm $(MUSE_LD_FLAGS)

POPARCH=amuse
POPDIR=$(CURDIR)/eSalsa-POP

POPFFLAGS := $(POPFFLAGS) $(NETCDF_FLAGS) $(NETCDFF_FLAGS) $(NETCDF_LIBS) $(NETCDFF_LIBS)
FS_FLAGS += $(NETCDF_LIBS) $(NETCDFF_LIBS)

CODE_GENERATOR ?= $(AMUSE_DIR)/build.py

POPREPO=https://github.com/nlesc-smcm/eSalsa-POP
POPREPOBRANCH=master


.PHONY: all download src_%/lib

all: eSalsa-POP $(TARGETS)

update: eSalsa-POP
	cd eSalsa-POP; git pull

ifdef DOWNLOAD_CODES
eSalsa-POP:
	make -C . download
else
eSalsa-POP:
	@echo "ERROR: need to download eSalsa-POP, but DOWNLOAD_CODES is not set" 
	@false 
endif

download:
	git clone -b $(POPREPOBRANCH) $(POPREPO)

clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h *.mod
	$(RM) *~ worker_code worker_code.f90 
	$(RM) -rf src_*/compile
	$(RM) -f src_*/libpop.a
	rm -f pop_worker_*

distclean: clean
	rm -rf eSalsa-POP

src_%/libpop.a: src_%/POP_DomainSizeMod.F90 eSalsa-POP
	make -C src_$* POPARCH=$(POPARCH) POPEXEDIR="$(CURDIR)/src_$*" POPDIR="$(POPDIR)" AMUSE_DIR="$(realpath $(AMUSE_DIR))" POPFFLAGS="$(POPFFLAGS)" MPIFC="$(MPIFC)"

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py POPInterface -o $@

$(TARGETS): pop_worker_%: worker_code.f90 src_%/libpop.a interface.f90
	$(MPIFC) -Isrc_$*/compile/ $(FCFLAGS) -c interface.f90 -o interface_$*.o 
	$(MPIFC) -Isrc_$*/compile/ $(FCFLAGS) $(FS_FLAGS) $< interface_$*.o src_$*/libpop.a $(FS_LIBS) -o $@
