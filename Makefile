# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=/Users/brunnabend/amuse-svn/amuse-svn/
-include $(AMUSE_DIR)/config.mk

MPIFC ?= mpif90
FORTRAN=$(FC)

ifneq (,$(findstring gfortran, $(notdir $(FORTRAN))))
FISHPACK_FLAGS = -fdefault-real-8 
export FISHPACK_FLAGS
LDFLAGS  += -L./src/fishpack4.1/lib/ -lfishpack
endif

ifeq ($(findstring ifort, $(notdir $(FORTRAN))), ifort)
# ifort flags
LDFLAGS  += -lm -mkl -I./src/include
FCFLAGS += -mkl -I./src/include
endif



OBJS = interface.o

CODELIB = src/libqgmodel.a

CODE_GENERATOR = $(AMUSE_DIR)/build.py

all: qgmodel_worker 

clean:
	$(RM) -f *.so *.o *.pyc *.mod worker_code.cc worker_code.h 
	$(RM) *~ worker_code worker_code.f90 qgmodel_worker
	make -C src clean

$(CODELIB):
	make -C src all

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py QGmodelInterface -o $@

qgmodel_worker: worker_code.f90 $(CODELIB) $(OBJS)
	$(MPIFC) $(FCFLAGS) $(FS_FLAGS) $< $(OBJS) $(CODELIB) $(FS_LIBS) $(LDFLAGS) -o $@

%.o: %.f90
	$(MPIFC) $(FCFLAGS) -c -o $@ $<
