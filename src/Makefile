MPIFC ?= mpif90
FORTRAN= $(FC)
ifneq (,$(findstring gfortran, $(notdir $(FORTRAN))))
FISHPACK_FLAGS += -fdefault-real-8 
export FISHPACK_FLAGS
CODEOBJS = beta.o  d_Poisson_2D_f_fishpack.o  main.o  vis_bot.o  chi.o \
  jacobian.o vis_lat.o
endif

ifeq ($(findstring ifort, $(notdir $(FORTRAN))), ifort)
# ifort flags
LDFLAGS  += -lm -mkl 
FCFLAGS += -mkl -I./include
CODEOBJS = beta.o d_Poisson_2D_f.o  main.o  vis_bot.o  chi.o \
  jacobian.o vis_lat.o
endif


CODELIB = libqgmodel.a

AR = ar ruv
RANLIB = ranlib
RM = rm

all: mkl $(CODELIB)


clean:
	$(RM) -f *.o *.a *.mod a.out
	make -C include clean
	make -C fishpack4.1 clean

$(CODELIB): $(CODEOBJS)
	$(RM) -f $@
	$(AR) $@ $(CODEOBJS)
	$(RANLIB) $@

%.o: %.f90
	$(FC) $(FCFLAGS) -c -o $@ $<

mkl:
ifeq ($(findstring ifort, $(notdir $(FORTRAN))), ifort)
	make -C include all
else
	make -C fishpack4.1 all FC="$(FC) $(FISHPACK_FLAGS)"
endif
