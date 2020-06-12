# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=$(shell amusifier --get-amuse-dir)
-include $(AMUSE_DIR)/config.mk

-include src/opt/include/Config.mk

.PHONY: all update clean distclean

export PYTHONPATH := $(PYTHONPATH):$(AMUSE_DIR)/src:$(AMUSE_DIR)/test
CODE_GENERATOR:=amusifier
MPICXX   ?= mpicxx

CFLAGS   += -Wall -g

IEMICLIBPATH=$(realpath src/opt/lib)

UNAME := $(shell uname -s)
ifeq ($(UNAME),Darwin)
CODELIB=$(IEMICLIBPATH)/libiemic.dylib
else
CODELIB=$(IEMICLIBPATH)/libiemic.so

# explicitly locate parallel hdf5
#~ LDFLAGS  += -L/usr/lib64/mpi/gcc/openmpi2/lib64 -lhdf5
LDFLAGS += $(HDF5_LIBS)
endif

export LDFLAGS

LIBS += -L$(IEMICLIBPATH) -liemic -lcoupledmodel -ljdqz_tools -lutils -lglobaldefs -lifpack_mrilu -lmrilucpp

OBJS = interface.o paramset.o

all: iemic_worker

update: | src/jdqzpp/ src/i-emic/
	cd src/jdqzpp; git pull
	cd src/i-emic; git pull

src/:
	mkdir -p $@

src/jdqzpp/: | src/
ifdef DOWNLOAD_CODES
	cd src; git clone https://github.com/erik808/jdqzpp.git
else
	@echo "ERROR: need to download JDQZPP, but DOWNLOAD_CODES is not set"
	@false
endif

src/jdqzpp/build/: | src/jdqzpp/
	mkdir -p $@
	cd $@; cmake .. -DCMAKE_INSTALL_RPATH="$(realpath src)/opt/" -DCMAKE_INSTALL_PREFIX:PATH="$(realpath src)/opt/"

src/opt/lib/cmake/JDQZPP/JDQZPPConfig.cmake: | src/jdqzpp/build/
	cd $|; make install

src/i-emic/: | src/
ifdef DOWNLOAD_CODES
	cd src; git clone https://github.com/nlesc-smcm/i-emic.git
else
	@echo "ERROR: need to download i-emic, but DOWNLOAD_CODES is not set"
	@false
endif

src/i-emic/build/: src/opt/lib/cmake/JDQZPP/JDQZPPConfig.cmake | src/i-emic/
	mkdir -p $@
	cd $@; cmake ..  -DCMAKE_INSTALL_RPATH="$(realpath src)/opt/" -DCMAKE_INSTALL_PREFIX:PATH="$(realpath src)/opt/" \
	    -DJDQZPP_DIR="$(realpath src)/opt/lib/cmake/JDQZPP/"

$(CODELIB) src/opt/include/Config.mk: | src/i-emic/build/
	cd $|; make install

clean:
	$(RM) -rf __pycache__
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h
	$(RM) -f *~ iemic_worker worker_code.cc

distclean: clean
	$(RM) -rf src/opt src/jdqzpp src/i-emic

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py iemicInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H interface.py iemicInterface -o $@

iemic_worker: worker_code.o $(OBJS) $(CODELIB)
	$(MPICXX) -Wl,-rpath $(IEMICLIBPATH) $^ -o $@ $(LDFLAGS) $(LIBS)

%.o: %.cc $(CODELIB) src/opt/include/Config.mk
	$(MPICXX) $(CXXFLAGS) -std=c++11 -Isrc/opt/include/ -c -o $@ $<

worker_code.o: worker_code.h
interface.o: interface.hpp
