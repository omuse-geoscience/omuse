FC ?= ifort
MPIFC ?= mpif90

libamuse_swan.a: work config code 
	ar crs $@ work/*.o

amuse_swan.o: amuse_swan.f90
	$(FC) -Iwork $(FCFLAGS) -c -o $@ $<

config: work
	cp macros.inc work/macros.inc

code:
	make -C work/ ser MPIFC="$(MPIFC)"

work:
	cp -r swan work
#	cd work; patch -p1 < ../swan.patch

clean:
	rm -f *.mod *.o libamuse_swan.a
	make -C work/ clean

distclean: clean
	rm -rf work
