# standard amuse configuration include
# config.mk will be made after ./configure has run
ifeq ($(origin AMUSE_DIR), undefined)
  AMUSE_DIR := $(shell amusifier --get-amuse-dir)
endif
-include $(AMUSE_DIR)/config.mk

ifneq ($(MAKECMDGOALS),download)
-include src/opt/include/Config.mk
else
DOWNLOAD_CODES:=true
endif

.PHONY: all download update clean distclean FORCE

export PYTHONPATH := $(PYTHONPATH):$(AMUSE_DIR)/src:$(AMUSE_DIR)/test
CODE_GENERATOR:=amusifier
MPICXX   ?= mpicxx

CFLAGS   += -Wall -g
LIBS:=iemic coupledmodel globaldefs utils ifpack_mrilu mrilucpp

UNAME := $(shell uname -s)
ifeq ($(UNAME),Darwin)
LIBEXT:=dylib
WORKER_RPATH:=@executable_path
else
LIBEXT=so
WORKER_RPATH:=$$ORIGIN

# explicitly locate parallel hdf5
LDFLAGS += $(HDF5_LIBS)
export LDFLAGS
endif

all: iemic_worker

download: update

update: | src/jdqzpp/ src/i-emic/
	cd src/jdqzpp; git pull
	cd src/i-emic; git pull

src/:
	mkdir -p $@

data/:
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
	cd $@; cmake .. -DCMAKE_INSTALL_RPATH="$(realpath src)/opt/lib" -DCMAKE_INSTALL_PREFIX:PATH="$(realpath src)/opt/"

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
	cd $@; cmake ..  -DCMAKE_INSTALL_RPATH='$(WORKER_RPATH)' -DCMAKE_INSTALL_PREFIX:PATH="$(realpath src)/opt/" \
	    -DJDQZPP_DIR="$(realpath src)/opt/lib/cmake/JDQZPP/"

IEMICLIBS:=$(LIBS:%=src/opt/lib/lib%.$(LIBEXT))
$(IEMICLIBS) src/opt/include/Config.mk: FORCE | src/i-emic/build/
	cd $|; make install

clean:
	$(RM) -rf __pycache__ data
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h
	$(RM) -f *~ iemic_worker worker_code.cc

distclean: clean
	$(RM) -rf src/opt src/jdqzpp src/i-emic

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py iemicInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H interface.py iemicInterface -o $@

iemic_worker: worker_code.o interface.o paramset.o $(LIBS:%=data/lib%.$(LIBEXT))
	$(MPICXX) -Wl,-rpath,'$(WORKER_RPATH)/data/' $(filter %.o,$^) -o $@ $(LDFLAGS) -Ldata/ $(LIBS:%=-l%)

data/lib%.$(LIBEXT): src/opt/lib/lib%.$(LIBEXT) | data/
	cp $< $@

%.o: %.cc $(CODELIB) src/opt/include/Config.mk
	$(MPICXX) $(CXXFLAGS) -std=c++11 -Isrc/opt/include/ -c -o $@ $<

paramset.o: paramset.hpp
worker_code.o: worker_code.h
interface.o: interface.hpp paramset.hpp

FORCE:
