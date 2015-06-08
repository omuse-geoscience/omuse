libamuse_swan.a: $(filter-out swanmain.o, swan/*.o)
	ar crs $@ $^

amuse_swan.o: amuse_swan.f90
	$(FC) -Iswan $(FCFLAGS) -c -o $@ $<

code:
	make -C swan/ config ser 

clean:
	make -C swan/ clean
	rm -f *.mod *.o libamuse_swan.a

distclean: clean
	make -C swan/ clobber
