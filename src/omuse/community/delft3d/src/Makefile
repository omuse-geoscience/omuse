MPIFC ?= mpif90
FC      = $(MPIFC)

LDFLAGS  += -lm $(MUSE_LD_FLAGS)

CODELIB = libdelf3dfm.a

CODEOBJS = test.o

INSTALLDIR = $(PWD)

AR = ar ruv
RANLIB = ranlib
RM = rm

all: libdflowfm_omuse.a include

clean:
	$(RM) -f include libdflowfm.a libdflowfm_omuse.a
	$(RM) -f *.o *.a

distclean: clean

$(CODELIB): $(CODEOBJS)
	$(RM) -f $@
	$(AR) $@ $(CODEOBJS)
	$(RANLIB) $@

%.o: %.f90
	$(FC) $(FCFLAGS) -c -o $@ $<

delft3dfm:
	$(error delft3dfm source not found )

delft3dfm/src/engines_gpl/dflowfm/packages/dflowfm_omuse:
	cp -r dflowfm_omuse delft3dfm/src/engines_gpl/dflowfm/packages/

delft3dfm/src/Makefile:  | delft3dfm delft3dfm/src/engines_gpl/dflowfm/packages/dflowfm_omuse
	cd delft3dfm/src; ./autogen.sh
	cd delft3dfm/src/third_party_open/kdtree2; ./autogen.sh
	cd delft3dfm/src; ./configure --prefix=$(INSTALLDIR) --with-netcdf --with-mpi --with-metis --with-petsc

libdflowfm_omuse.a: delft3dfm/src/Makefile
#~ 	make -C delft3dfm/src/
	make -C delft3dfm/src/engines_gpl/dflowfm/
	ln -s delft3dfm/src/engines_gpl/dflowfm/packages/dflowfm_kernel/src/.libs/libdflowfm.a libdflowfm.a
	ln -s delft3dfm/src/engines_gpl/dflowfm/packages/dflowfm_omuse/src/.libs/libdflowfm_omuse.a libdflowfm_omuse.a

include: libdflowfm_omuse.a
	ln -s delft3dfm/src/engines_gpl/dflowfm/packages/dflowfm_omuse/src include
