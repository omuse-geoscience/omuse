# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../../..
-include $(AMUSE_DIR)/config.mk

MPIFC ?= mpif90

#MPIFC = ifort -I/cm/shared/apps/mpich/ge/open64/64/3.1/include -L/cm/shared/apps/mpich/ge/open64/64/3.1/lib -lmpich 

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
POPEXEDIR1=$(CURDIR)/src_320x384x40
POPEXEDIR2=$(CURDIR)/src_3600x2400x42
POPDIR=$(CURDIR)/eSalsa-POP

#-------------- NETCDF inserts
#NETCDFPATH = /cm/shared/apps/netcdf/open64/64/4.3.1.1/
NETCDFPATH = /home/ben
NETCDFINC = -I$(NETCDFPATH)/include
NETCDFLIB = -L/cm/shared/apps/netcdf/open64/64/4.3.1.1/lib
LDFLAGS := $(LDFLAGS) $(NETCDFINC) $(NETCDFLIB)
FS_LIBS := $(FS_LIBS) $(NETCDFLIB) -lnetcdf -lnetcdff
POPFFLAGS := $(POPFFLAGS) $(NETCDFINC) $(NETCDFLIB)
#-------------- End of NETCDF inserts

#FS_LIBS+=-lnetcdf -lnetcdff



POPINCPATH1 = -Isrc_320x384x40/compile/
POPINCPATH2 = -Isrc_3600x2400x42/compile/

#OBJS = interface.o 
CODELIB1 = src_320x384x40/libpop.a
CODELIB2 = src_3600x2400x42/libpop.a
CODE_GENERATOR = $(AMUSE_DIR)/build.py

all: pop_worker_320x384x40 pop_worker_3600x2400x42

update: eSalsa-POP
	cd eSalsa-POP; git pull

eSalsa-POP:
	git clone -b amuse https://github.com/NLeSC/eSalsa-POP.git

clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h *.mod
	$(RM) *~ worker_code worker_code.f90 
	make -C $(POPEXEDIR1) clean POPARCH=$(POPARCH) POPEXEDIR=$(POPEXEDIR1) POPDIR=$(POPDIR)
	make -C $(POPEXEDIR2) clean POPARCH=$(POPARCH) POPEXEDIR=$(POPEXEDIR2) POPDIR=$(POPDIR)
	rm -f pop_worker_*

distclean: clean
	rm -rf eSalsa-POP

$(CODELIB1):
	make -C $(POPEXEDIR1) POPARCH=$(POPARCH) POPEXEDIR=$(POPEXEDIR1) POPDIR=$(POPDIR) AMUSE_DIR="$(realpath $(AMUSE_DIR))" POPFFLAGS="$(POPFFLAGS)" MPIFC="$(MPIFC)"

$(CODELIB2):
	make -C $(POPEXEDIR2) POPARCH=$(POPARCH) POPEXEDIR=$(POPEXEDIR2) POPDIR=$(POPDIR) AMUSE_DIR="$(realpath $(AMUSE_DIR))" POPFFLAGS="$(POPFFLAGS)" MPIFC="$(MPIFC)"

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py POPInterface -o $@

pop_worker_320x384x40: worker_code.f90 $(CODELIB1) $(OBJS) interface.f90
	$(MPIFC) $(POPINCPATH1) $(FCFLAGS) -c interface.f90 -o interface.o
	$(MPIFC) $(POPINCPATH1) $(FCFLAGS) $(FS_FLAGS) $< interface.o $(CODELIB1) $(FS_LIBS) -o $@

pop_worker_3600x2400x42: worker_code.f90 $(CODELIB2) $(OBJS) interface.f90
	$(MPIFC) $(POPINCPATH2) $(FCFLAGS) -c interface.f90 -o interface.o 
	$(MPIFC) $(POPINCPATH2) $(FCFLAGS) $(FS_FLAGS) $< interface.o $(CODELIB2) $(FS_LIBS) -o $@



