FC ?= ifort
MPIFC ?= mpif90

.PHONY: work config code objects

libamuse_swan.a: work config objects
	ar crs $@ work/*.o

amuse_swan.o: amuse_swan.f90
	$(FC) -Iwork $(FCFLAGS) -c -o $@ $<

amuse_swan.f90: work amuse_swan.ftn90
	$(eval $(shell grep 'swch =' macros.inc))
	perl work/switch.pl $(swch) amuse_swan.ftn90

config: work
	cp macros.inc work/macros.inc

objects: config
	-cd work; patch -N -p1 < ../swan.patch
	make -C work amuse MPIFC="$(MPIFC)"

code: config
	-cd work; patch -N -p1 -R < ../swan.patch
	make -C work/ ser MPIFC="$(MPIFC)"

work:
	cp -r swan work

clean:
	rm -f *.mod *.o libamuse_swan.a
	make -C work/ clean

distclean: clean
	rm -rf work
