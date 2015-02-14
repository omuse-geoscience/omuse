AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

all: src/libadcirc.a

src/libadcirc.a:
	make -C src/ libadcirc.a LIBADC=libadcirc.a BUILDTYPE=libadcirc.a FC="$(FC)" PFC="$(MPIFC)" compiler="amuse"

clean:

distclean: clean
	make -C src/ LIBADC=libadcirc.a clobber
