# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../../..
-include $(AMUSE_DIR)/config.mk

MPIFC ?= mpif90

ifeq ($(findstring gfortran, $(notdir $(FC))), gfortran)
FCFLAGS+=-fconvert=swap
POPFFLAGS=-O2 -Wall -fdefault-double-8 -fdefault-real-8 -fconvert=swap -fimplicit-none -fbounds-check
endif

ifeq ($(findstring ifort, $(notdir $(FC))), ifort)
FCFLAGS+=-convert big_endian 
POPFFLAGS=-O2 -r8 -convert big_endian
endif

LDFLAGS  += -lm $(MUSE_LD_FLAGS)

POPARCH=amuse
POPEXEDIR=$(CURDIR)/src
POPDIR=$(CURDIR)/eSalsa-POP

#-------------- NETCDF inserts
NETCDFPATH = /opt/local
NETCDFINC = -I$(NETCDFPATH)/include
NETCDFLIB = -L$(NETCDFPATH)/lib
#~ LDFLAGS := $(LDFLAGS) $(NETCDFINC)
#~ FS_LIBS := $(FS_LIBS) $(NETCDFLIB) -lnetcdf -lnetcdff
#-------------- End of NETCDF inserts

FS_LIBS+=$(NETCDFLIB) -lnetcdf -lnetcdff


POPINCPATH = -Isrc/compile/

OBJS = interface.o
CODELIB = src/libpop.a
CODE_GENERATOR = $(AMUSE_DIR)/build.py

all: pop_worker 

update: eSalsa-POP
	cd eSalsa-POP; git pull

eSalsa-POP:
	git clone -b amuse https://github.com/NLeSC/eSalsa-POP.git

clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h *.mod
	$(RM) *~ worker_code worker_code.f90 
	make -C src clean POPARCH=$(POPARCH) POPEXEDIR=$(POPEXEDIR) POPDIR=$(POPDIR)
	rm -f pop_worker

distclean: clean
	rm -rf eSalsa-POP

$(CODELIB):
	make -C src POPARCH=$(POPARCH) POPEXEDIR=$(POPEXEDIR) POPDIR=$(POPDIR) AMUSE_DIR="$(realpath $(AMUSE_DIR))" POPFFLAGS="$(POPFFLAGS)"

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py POPInterface -o $@

pop_worker: worker_code.f90 $(CODELIB) $(OBJS)
	$(MPIFC) $(POPINCPATH) $(FCFLAGS) $(FS_FLAGS) $< $(OBJS) $(CODELIB) $(FS_LIBS) -o $@

%.o: %.f90
	$(MPIFC) $(POPINCPATH) $(FCFLAGS) -c -o $@ $<
