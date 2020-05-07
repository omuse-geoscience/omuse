FC ?= gfortran
MPIFC ?= $(FC)

.PHONY: patch apply_patch apply_reverse_patch config code objects

libamuse_swan.a: work config objects
	ar crs $@ work/*.o

amuse_swan.o: amuse_swan.f90
	$(FC) -Iwork $(FCFLAGS) -c -o $@ $<

amuse_swan.f90: work amuse_swan.ftn90
	$(eval $(shell grep 'swch =' macros.inc))
	perl work/switch.pl $(swch) amuse_swan.ftn90

config: work
	cp macros.inc work/macros.inc

objects: config apply_patch
	make -C work amuse MPIFC="$(MPIFC)"

code: config apply_reverse_patch
	make -C work/ ser MPIFC="$(MPIFC)"

work:
	cp -r swan work

apply_patch:
	-cd work; patch -N -p1 < ../swan.patch

apply_reverse_patch:
	-cd work; patch -N -p1 -R < ../swan.patch

patch: clean
	diff -uN swan work > new.patch || exit 0

clean:
	rm -f *.mod *.o libamuse_swan.a
	make -C work/ clobber

distclean: clean
	rm -rf work
