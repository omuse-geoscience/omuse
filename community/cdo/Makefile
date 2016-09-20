# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include $(AMUSE_DIR)/config.mk

CFLAGS   += -Wall -g
CXXFLAGS += $(CFLAGS)
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

OBJS = interface.o

CODE = cdo-1.7.0rc5

CODELIB = $(CODE)/src/libcdo-amuse.a

CDILIB = $(CODE)/libcdi/src/.libs/libcdi.a

CODE_GENERATOR = $(AMUSE_DIR)/build.py

all: remapper_worker

patches:
	$(RM) -rf $(CODE)
	tar zxf cdo-current.tar.gz
	cd $(CODE); ./configure --with-netcdf=$(NETCDF_PREFIX)
	cd .; ./make_patches.py diff
	$(RM) -rf $(CODE)

clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h
	$(RM) *~ remapper_worker worker_code.cc
#	make -C $(CODE) clean

distclean: clean
	rm -rf $(CODE) src

$(CODE):
	tar zxf cdo-current.tar.gz
	ln -s `tar -tf cdo-current.tar.gz | awk -F/ '{print $1;exit}'` src
	cd $(CODE); ./configure --disable-openmp --with-netcdf=$(NETCDF_PREFIX) FC="$(FC)" FCFLAGS="$(FCFLAGS)"
	cd .; ./make_patches.py patch
	cp script.ar $(CODE)/src/script.ar

$(CODELIB): $(CODE)
	make -C $(CODE) 
	make -C $(CODE)/src libcdo-cdolib.a
	cd $(CODE)/src/; ar rvsT ./libcdo-amuse.a libcdo-cdolib.a .libs/libcdo.a
	cd $(CODE)/src/; ranlib libcdo-amuse.a

worker_code.cc: interface.py worker_code.h
	$(CODE_GENERATOR) --type=c interface.py CDOInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H interface.py CDOInterface -o $@

remapper_worker: worker_code.o $(CODELIB) $(OBJS)
	$(eval $(shell grep '\-luuid' cdo-1.7.0rc5/libcdi/Makefile))
	$(MPICXX) $(CXXFLAGS) -z muldefs $< $(OBJS) $(CODELIB) $(CDILIB)  -lm -lz $(NETCDF_LIBS) $(UUID_C_LIB) $(OPENMP_CFLAGS) -o $@

.cc.o: $<
	$(MPICXX) $(CXXFLAGS) -c -o $@ $<

.PHONY: patches 
