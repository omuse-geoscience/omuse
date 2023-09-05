ifeq ($(origin AMUSE_DIR), undefined)
  AMUSE_DIR := $(shell amusifier --get-amuse-dir)
endif
-include ${AMUSE_DIR}/config.mk
.NOTPARALLEL:

export PYTHONPATH := $(PYTHONPATH):$(AMUSE_DIR)/src:$(AMUSE_DIR)/test
MPIFC ?= mpifort #/usr/local/bin/mpif90

TARGETS := $(patsubst src_%, pop_worker_%, ${wildcard src_*})

ifeq ($(findstring gfortran, $(notdir $(FC))), gfortran)
FCFLAGS+=-fconvert=swap -std=legacy -fmax-stack-var-size=1048576
POPFFLAGS=-O2 -Wall -fdefault-double-8 -fdefault-real-8 -std=legacy -fconvert=swap -fimplicit-none -fbounds-check  -fmax-stack-var-size=1048576
endif

ifeq ($(findstring ifort, $(notdir $(FC))), ifort)
FCFLAGS+=-O2 -r8 -g -convert big_endian -assume byterecl
POPFFLAGS=-O2 -r8 -g -convert big_endian -assume byterecl
endif

LDFLAGS  += -lm $(MUSE_LD_FLAGS)

POPARCH=amuse
POPDIR=$(CURDIR)/eSalsa-POP

POPFFLAGS := $(POPFFLAGS) $(NETCDF_FLAGS) $(NETCDFF_FLAGS) $(NETCDF_LIBS) $(NETCDFF_LIBS)
FS_LIBS += $(NETCDF_LIBS) $(NETCDFF_LIBS)

CODE_GENERATOR ?= $(AMUSE_DIR)/build.py

POPREPO=https://github.com/nlesc-smcm/eSalsa-POP
POPREPOBRANCH=master

.PHONY: all download src_%

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
	make -j1 -C src_$* POPARCH=$(POPARCH) POPEXEDIR="$(CURDIR)/src_$*" POPDIR="$(POPDIR)" AMUSE_DIR="$(realpath $(AMUSE_DIR))" POPFFLAGS="$(POPFFLAGS)" MPIFC="$(MPIFC)"

src_%:
	mkdir src_$*

src_%/POP_DomainSizeMod.F90: | src_%
	cp _src_template/GNUmakefile src_$*
	cp _src_template/domain_size.F90 src_$*
	cd src_$*; ${PYTHON} ../_src_template/POP_DomainSizeMod.py $*

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py POPInterface -o $@

#~ $(TARGETS): pop_worker_%: worker_code.f90 src_%/libpop.a interface.f90
pop_worker_%: worker_code.f90 src_%/libpop.a interface.f90
	$(MPIFC) -Isrc_$*/compile/ $(FCFLAGS) -c interface.f90 -o interface_$*.o 
	$(MPIFC) -Isrc_$*/compile/ $(FCFLAGS) $(FS_FLAGS) $< interface_$*.o src_$*/libpop.a $(FS_LIBS) -o $@

.PRECIOUS: src_% src_%/libpop.a src_%/POP_DomainSizeMod.F90
