AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

TARGETS := $(patsubst %.config, q-gcm_worker_%, ${shell ls *.config})

BUILDDIR = build

SRCDIR=./src/q-gcm/src

CODE_GENERATOR = $(AMUSE_DIR)/build.py

ifneq ($(MAKEFILE_OPTIONS_FILE), )
include $(MAKEFILE_OPTIONS_FILE)
endif

.SECONDARY:

.PHONY: all $(BUILDDIR)_%/lib download

all: download q-gcm_worker

q-gcm_worker: ${TARGETS}

q-gcm_%: % ;

FFLAGS ?= -g -O2
FFLAGS2=$(FFLAGS) $(NETCDFF_FLAGS) -80 -fpp -std90
LIBS+=$(NETCDFF_LIBS)

q-gcm_worker_%: worker_code.f90 $(BUILDDIR)_%/libq-gcm.a $(BUILDDIR)_%/q-gcm_utility.o $(BUILDDIR)_%/amuse_q-gcm.o interface_%.o 
	$(MPIFC) $(FFLAGS) $(SC_FLAGS) $(FS_FLAGS) $(LDFLAGS) -I$(BUILDDIR)_$* -o $@ $< interface_$*.o $(BUILDDIR)_$*/amuse_q-gcm.o $(BUILDDIR)_$*/q-gcm_utility.o $(BUILDDIR)_$*/libq-gcm.a $(LIBS) $(SC_FCLIBS) $(FS_LIBS) $(LIBS)
 
interface_%.o: interface.f90 getter_setters.f90	default_param.f90 Makefile $(BUILDDIR)_%/amuse_q-gcm.o
	$(MPIFC) $(FFLAGS) -I$(BUILDDIR)_$*  -c -o $@ $< 
 
$(BUILDDIR)_%:
	-mkdir $@
	cp -r $(SRCDIR)/* $@/

$(BUILDDIR)_%/Makefile: $(SRCDIR)/Makefile $(BUILDDIR)_%/make.config $(BUILDDIR)_%/make.macro | $(BUILDDIR)_%
	cp $(SRCDIR)/Makefile $@

$(BUILDDIR)_%/make.config: %.config | $(BUILDDIR)_%
	cp $*.config $@

$(BUILDDIR)_%/make.macro: make.macro | $(BUILDDIR)_%
	cp make.macro $@

$(BUILDDIR)_%/parameters_data.F: %_parameters_data.F | $(BUILDDIR)_%
	cp $*_parameters_data.F $@

$(BUILDDIR)_%/q-gcm: $(BUILDDIR)_%/Makefile $(BUILDDIR)_%/parameters_data.F 
	make -C $(BUILDDIR)_$* q-gcm FC="$(FC)" FFLAGS="$(FFLAGS2)" NETCDFF_PREFIX="$(NETCDFF_PREFIX)"

$(BUILDDIR)_%/libq-gcm.a: $(BUILDDIR)_%/q-gcm $(BUILDDIR)_%/lib ;

$(BUILDDIR)_%/lib: OBJS=$(filter-out  $(BUILDDIR)_$*/q-gcm.o,$(shell ls $(BUILDDIR)_$*/*.o ))
$(BUILDDIR)_%/lib: 
#~ $(BUILDDIR)_%/libq-gcm.a: OBJS=$(shell ls $(BUILDDIR)_%/*.o )
	ar crs $(BUILDDIR)_$*/libq-gcm.a $(OBJS)

$(BUILDDIR)_%/amuse_q-gcm.o: src/amuse_q-gcm.F $(BUILDDIR)_%/make.config $(BUILDDIR)_%/libq-gcm.a
ifeq ($(MAKEFILE_OPTIONS_FILE), )
	make -C . $@ MAKEFILE_OPTIONS_FILE=$(BUILDDIR)_$*/make.config
else
	$(FC) $(FFLAGS2) $(QGOPTS) -I$(BUILDDIR)_$* $(NETCDFF_FLAGS) -c -o $@ $<
endif  

$(BUILDDIR)_%/q-gcm_utility.o: src/q-gcm_utility.F $(BUILDDIR)_%/make.config $(BUILDDIR)_%/libq-gcm.a
ifeq ($(MAKEFILE_OPTIONS_FILE), )
	make -C . $@ MAKEFILE_OPTIONS_FILE=$(BUILDDIR)_$*/make.config
else
	$(FC) $(FFLAGS2) $(QGOPTS) -I$(BUILDDIR)_$* $(NETCDFF_FLAGS) -c -o $@ $<
endif

worker_code.f90: interface.py $(AMUSE_INTERFACE_DEPS)
	$(CODE_GENERATOR) --type=f90 --needs-mpi=false $< QGCMInterface -o $@

%.o: %.f90 Makefile
	$(MPIFC) $(FFLAGS)  -c -o $@ $< 

getter_setters.f90: interface.py
	python -c "import interface; interface.code_generator.generate_getters_setters()"

update_q-gcm: src/q-gcm
	cd src/q-gcm; git pull

src/q-gcm:
	cd src;git clone -b omuse https://github.com/ipelupessy/q-gcm.git

download: src/q-gcm

clean:
	-rm -Rf $(BUILDDIR)_*
	-rm -f worker_code.f90 getter_setters.f90
	-rm -f interface*.o q-gcm_worker* *.mod

distclean: clean
	-rm -Rf src/q-gcm

