SRC=blktri.f cblktri.f cmgnbn.f comf.f fftpack.f genbun.f gnbnaux.f \
    hstcrt.f hstcsp.f hstcyl.f hstplr.f hstssp.f hw3crt.f hwscrt.f \
    hwscsp.f hwscyl.f hwsplr.f hwsssp.f pois3d.f poistg.f sepaux.f \
    sepeli.f sepx4.f

OBJ=$(subst .f,.o,$(SRC))
OBJS=$(addprefix ../objs/,$(OBJ))

../lib/libfishpack.a: $(OBJS)
	$(AR) -rv $@ $? 

../objs/%.o : %.f
	$(FC) $(FCFLAGS) -c $< -o $@

clean:
	rm -f ../lib/libfishpack.a $(OBJS)
