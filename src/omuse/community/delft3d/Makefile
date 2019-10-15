# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=$(shell amusifier --get-amuse-dir)
-include $(AMUSE_DIR)/config.mk

MPIFC ?= mpif90
FC      = $(MPIFC)

FCFLAGS += -Isrc/include
LDFLAGS += -lm $(MUSE_LD_FLAGS)

OBJS = interface.o

CODELIB = src/libdflowfm_omuse.a src/libdflowfm.a

LIBS = -ldl -lnetcdf -lnetcdff -lmetis -lpetsc -lstdc++

all: dflowfm_worker 

clean:
	$(RM) -rf __pycache__
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ worker_code worker_code.f90
	make -C src clean

distclean: clean
	make -C src distclean

$(CODELIB): BUILDLIB

BUILDLIB:
	make -C src all

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py DFlowFMInterface -o $@

dflowfm_worker: worker_code.f90 $(CODELIB) $(OBJS)
	$(MPIFC) $(FCFLAGS) $(FS_FLAGS) $< $(OBJS) $(CODELIB) $(FS_LIBS) -o $@ $(LIBS)

%.o: %.f90
	$(FC) $(FCFLAGS) -c -o $@ $<
