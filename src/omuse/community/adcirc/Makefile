AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

CODE_GENERATOR = $(AMUSE_DIR)/build.py

CLASSNAME=AdcircInterface

FCFLAGS += -I./src/ -I./src/odir3/

ifeq ($(findstring gfortran, $(notdir $(FC))), gfortran)
compiler="gnu"
endif

ifeq ($(findstring ifort, $(notdir $(FC))), ifort)
compiler="intel"
FCFLAGS+=-traceback 
endif


all: src prep src/libamuse_adcirc.a src/amuse_adcirc.o adcirc_worker

src:
	$(error "please link src to correct src_v[version] directory")  

src/adcirc_src/src/*.F: prep

src/libamuse_adcirc.a: src/adcirc_src/src/*.F
	make -C src/ libamuse_adcirc.a BUILDTYPE=amuse FC="$(FC)" PFC="$(MPIFC)" compiler=$(compiler) FSO=-FI

src/amuse_adcirc.o: src/amuse_adcirc.F90
	make -C src/ amuse_adcirc.o BUILDTYPE=amuse FC="$(FC)" PFC="$(MPIFC)" compiler=$(compiler) FSO=-FI

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 $< $(CLASSNAME) -o $@

adcirc_worker: worker_code.f90 interface.o src/amuse_adcirc.o src/libamuse_adcirc.a 
	$(MPIFC) $(FCFLAGS) $(SC_FLAGS) $(FS_FLAGS) $^ -o $@  $(LIBS) $(SC_FCLIBS) $(FS_LIBS)

prep:
	make -C src/ -f makefile_prep patched_source

%.o: %.f90
	$(FC) $(FCFLAGS) -c -o $@ $<

clean:
	-make -C src/ clean
	rm -f *.pyc *.mod
	rm -f interface.o adcirc_worker.f90 worker_code.f90
	rm -f adcirc_worker 

distclean: clean
	-make -C src/ LIBADC=libadcirc.a clobber
	make -C src/ -f makefile_prep distclean
	
download:
	make -C src/ -f makefile_prep download
	
update:
	make -C src/ -f makefile_prep update
