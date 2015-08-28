# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include $(AMUSE_DIR)/config.mk

MPICXX   ?= mpicxx

CFLAGS   += -Wall -g
CXXFLAGS += $(CFLAGS)
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

OBJS = interface.o

CODE = cdo-1.7.0rc5

CODELIB = src/src/cdo-biglib.a
#CODELIB = src/src/libcdo-amuse.a

CDOLIB = src/src/.libs/libcdo.a
CDILIB = src/libcdi/src/.libs/libcdi.a

CODE_GENERATOR = $(AMUSE_DIR)/build.py

all: remapper_worker

clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h
	$(RM) *~ remapper_worker worker_code.cc
#	make -C $(CODE) clean

$(CODE):
	tar zxf cdo-current.tar.gz
	cd $(CODE)
	./configure
	cd -

$(CODELIB):
	make -C $(CODE) 
	make -C $(CODE)/src libcdo-amuse.a

worker_code.cc: interface.py worker_code.h
	$(CODE_GENERATOR) --type=c interface.py CDOInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H interface.py CDOInterface -o $@

remapper_worker: worker_code.o $(CODELIB) $(OBJS)
#	$(LD) $< $(OBJS) $(CDOLIB) $(CODELIB) $(CDILIB) -lm -lz -lnetcdf -o $@
	$(MPICXX) $(CXXFLAGS) -z muldefs $< $(OBJS) $(CDOLIB) $(CODELIB) $(CDILIB) -lm -lz -lnetcdf -fopenmp -o $@

.cc.o: $<
	$(MPICXX) $(CXXFLAGS) -c -o $@ $<

