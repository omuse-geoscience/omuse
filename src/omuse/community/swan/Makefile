ifeq ($(origin AMUSE_DIR), undefined)
  AMUSE_DIR := $(shell amusifier --get-amuse-dir)
endif
-include ${AMUSE_DIR}/config.mk

CODE_GENERATOR = $(AMUSE_DIR)/build.py

CLASSNAME=SwanInterface

FCFLAGS+= -O2

FCFLAGS+=$(OPENMP_FCFLAGS)

FCINCS = -I./src/ -I./src/work/

all: src/libamuse_swan.a src/amuse_swan.o swan_worker
	
src/libamuse_swan.a: 
	make -C src/ libamuse_swan.a MPIFC="$(MPIFC)"

src/amuse_swan.o: src/libamuse_swan.a src/amuse_swan.ftn90
	make -C src/ amuse_swan.o MPIFC="$(MPIFC)"

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 $< $(CLASSNAME) -o $@

swan_worker: worker_code.f90 interface.o src/amuse_swan.o src/libamuse_swan.a 
	$(MPIFC) $(FCFLAGS) $(FCINCS) $(SC_FLAGS) $(FS_FLAGS) $^ -o $@  $(LIBS) $(SC_FCLIBS) $(FS_LIBS)

%.o: %.f90
	$(FC) $(FCFLAGS) $(FCINCS) -c -o $@ $<

interface.o: interface.f90 getter_setters.f90 src/amuse_swan.o src/libamuse_swan.a
	$(FC) $(FCFLAGS) $(FCINCS) -c -o $@ interface.f90

getter_setters.f90: interface.py
	python -c "import interface; interface.generate_getters_setters()"

clean:
	make -C src/ clean
	rm -f *.pyc *.mod
	rm -f interface.o swan_worker.f90 worker_code.f90
	rm -f swan_worker 

distclean: clean
	make -C src/ distclean
