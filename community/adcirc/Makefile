AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

CODE_GENERATOR = $(AMUSE_DIR)/build.py

CLASSNAME=AdcircInterface

FFLAGS += -I./src/ -I./src/odir3/

all: src/libadcirc.a src/amuse_adcirc.o adcirc_worker

src/libadcirc.a: adcirc/v50release_121024/src/*.F
	make -C src/ libadcirc.a LIBADC=libadcirc.a BUILDTYPE=libadcirc.a FC="$(FC)" PFC="$(MPIFC)" compiler="amuse"

src/amuse_adcirc.o: src/amuse_adcirc.F90
	make -C src/ amuse_adcirc.o LIBADC=libadcirc.a BUILDTYPE=libadcirc.a FC="$(FC)" PFC="$(MPIFC)" compiler="amuse"

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 $< $(CLASSNAME) -o $@

adcirc_worker: worker_code.f90 interface.o src/amuse_adcirc.o src/libadcirc.a 
	$(MPIFC) $(FFLAGS) $(SC_FLAGS) $(FS_FLAGS) $^ -o $@  $(LIBS) $(SC_FCLIBS) $(FS_LIBS)

%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $<

clean:
	make -C src/ clean
	rm -f *.pyc *.mod
	rm -f interface.o adcirc_worker.f90 worker_code.f90
	rm -f adcirc_worker 

distclean: clean
	make -C src/ LIBADC=libadcirc.a clobber
