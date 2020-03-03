ifeq ($(origin AMUSE_DIR), undefined)
  AMUSE_DIR := $(shell amusifier --get-amuse-dir)
endif
-include ${AMUSE_DIR}/config.mk

MPIFC ?= mpif90
FORTRAN=$(FC)

ifneq (,$(findstring gfortran, $(notdir $(FORTRAN))))
FISHPACK_FLAGS = -fdefault-real-8 
export FISHPACK_FLAGS
LDFLAGS  += -L./src/fishpack4.1/lib/ -lfishpack
FCFLAGS += -fno-automatic
endif

ifeq ($(findstring ifort, $(notdir $(FORTRAN))), ifort)
# ifort flags
LDFLAGS  += -O2 -lm -mkl -I./src/include
FCFLAGS += -O2 -mkl -I./src/include  -heap-arrays 1024
endif

OBJS = interface.o

CODELIB = src/libqgmodel.a

CODE_GENERATOR = $(AMUSE_DIR)/build.py

all: qgmodel_worker 

clean:
	$(RM) -f *.so *.o *.pyc *.mod worker_code.cc worker_code.h 
	$(RM) *~ worker_code worker_code.f90 qgmodel_worker
	$(RM) -rf qgmodel_worker.dSYM 
	make -C src clean

$(CODELIB):
	make -C src all

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py QGmodelInterface -o $@

qgmodel_worker: worker_code.f90 $(CODELIB) $(OBJS)
	$(MPIFC) $(FCFLAGS) $(FS_FLAGS) $< $(OBJS) $(CODELIB) $(FS_LIBS) $(LDFLAGS) -o $@

%.o: %.f90
	$(MPIFC) $(FCFLAGS) -c -o $@ $<
