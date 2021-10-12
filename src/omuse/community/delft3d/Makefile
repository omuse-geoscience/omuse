# standard amuse configuration include
# config.mk will be made after ./configure has run
ifeq ($(origin AMUSE_DIR), undefined)
  AMUSE_DIR := $(shell amusifier --get-amuse-dir)
endif
AMUSE_DIR?=$(shell amusifier --get-amuse-dir)
-include $(AMUSE_DIR)/config.mk

MPIFC ?= mpif90
FC    ?= $(MPIFC)

FCFLAGS += -I/usr/include -Isrc/include

OBJS = interface.o

CODELIB = src/libdflowfm_omuse.a src/libdflowfm.a

LIBS = -ldl -lnetcdf -lnetcdff -lmetis -lpetsc -lstdc++

all: dflowfm_worker 

clean:
	$(RM) -rf __pycache__ 
	$(RM) -f *.mod *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ worker_code worker_code.f90
	make -C src clean

distclean: clean
	$(RM) -f dflowfm_worker partitioner_worker
	make -C src distclean

$(CODELIB): BUILDLIB

BUILDLIB:
	make -C src all

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py DFlowFMInterface -o $@

partitioner_worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py PartitionerInterface -o $@

dflowfm_worker: worker_code.f90 $(CODELIB) $(OBJS)
	$(MPIFC) $(FCFLAGS) $(FS_FLAGS) $< $(OBJS) $(CODELIB) $(FS_LIBS) -o $@ $(LIBS)

partitioner_worker: partitioner_interface.o partitioner_worker_code.f90 $(CODELIB) 
	$(MPIFC) $(FCFLAGS) $(FS_FLAGS) partitioner_worker_code.f90 partitioner_interface.o $(CODELIB) $(FS_LIBS) -o $@ $(LIBS)

partitioner_getter_setters.f90: interface.py
	$(PYTHON) -c "import interface; interface.code_generator_partitioner.print_getters_setters()" > $@

partitioner_interface.o: partitioner_getter_setters.f90

%.o: %.f90
	$(FC) $(FCFLAGS) -c -o $@ $<
