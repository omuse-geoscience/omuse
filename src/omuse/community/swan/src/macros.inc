# template for amuse

ifeq ($(findstring gfortran, $(notdir $(FC))), gfortran)
FCFLAGS += -fno-second-underscore
OMPFLAG = -fopenmp
F90FLAGS = -ffree-line-length-none
endif

ifeq ($(findstring ifort, $(notdir $(FC))), ifort)
FCFLAGS += -W0 -assume byterecl -traceback -diag-disable 8290 -diag-disable 8291 -diag-disable 8293 
OMPFLAG = -openmp
F90FLAGS = 
endif

F90_SER = $(FC)
F90_OMP = $(FC)
F90_MPI = $(MPIFC)
FLAGS_OPT = -O2
FLAGS_MSC = $(FCFLAGS)
FLAGS90_MSC = $(FLAGS_MSC) $(F90FLAGS)
FLAGS_SER =
FLAGS_OMP = $(OMPFLAG)
FLAGS_MPI =
NETCDFROOT =
INCS_SER =
INCS_OMP =
INCS_MPI =
LIBS_SER =
LIBS_OMP =
LIBS_MPI =
NCF_OBJS =
O_DIR = ../work/odir4/
OUT = -o 
EXTO = o
MAKE = make
RM = rm -f
swch = -unix
