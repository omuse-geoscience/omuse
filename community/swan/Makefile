AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

CODE_GENERATOR = $(AMUSE_DIR)/build.py

CLASSNAME=SwanInterface

FCFLAGS += -I./src/ -I./src/swan/

ifeq ($(findstring gfortran, $(notdir $(FC))), gfortran)
compiler="gnu"
endif

ifeq ($(findstring ifort, $(notdir $(FC))), ifort)
compiler="intel"
endif

all: src src/swan src/libamuse_swan.a src/amuse_swan.o swan_worker

src/swan:
	cp -r swan/swan4101 src/swan
	
src/libamuse_swan.a: 
	make -C src/ libamuse_swan.a 
#BUILDTYPE=amuse FC="$(FC)" PFC="$(MPIFC)" compiler=$(compiler) FSO=-FI

src/amuse_swan.o: src/amuse_swan.f90
	make -C src/ amuse_swan.o
#BUILDTYPE=amuse FC="$(FC)" PFC="$(MPIFC)" compiler=$(compiler) FSO=-FI

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 $< $(CLASSNAME) -o $@

swan_worker: worker_code.f90 interface.o src/amuse_swan.o src/libamuse_swan.a 
	$(MPIFC) $(FCFLAGS) $(SC_FLAGS) $(FS_FLAGS) $^ -o $@  $(LIBS) $(SC_FCLIBS) $(FS_LIBS)

%.o: %.f90
	$(FC) $(FCFLAGS) -c -o $@ $<

clean:
	make -C src/ clean
	rm -f *.pyc *.mod
	rm -f interface.o swan_worker.f90 worker_code.f90
	rm -f swan_worker 

distclean: clean
	make -C src/ distclean
