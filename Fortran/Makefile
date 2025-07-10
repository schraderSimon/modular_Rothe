# --------- Makefile -------------------------------------------------
FC       := gfortran
STD_INC  := /usr/local/include/fortran_stdlib/GNU-13.3.0
LIBS     := -lfortran_stdlib -lopenblas

OPT      := -ffree-line-length-none -Ofast -march=native -funroll-loops -ffast-math
WARN     := -std=f2018 -pedantic

# Object files
OBJS     := lbfgs.o ho_kernels.o RotheSolver_mod.o optimizers.o propagate_Gaussians.o

# Final target
all: propagate_Gaussians

propagate_Gaussians: $(OBJS)
	$(FC) $(OBJS) -I$(STD_INC) $(LIBS) -o $@

# Explicit dependencies
ho_kernels.o: ho_kernels.f90
lbfgs.o: lbfgs.f
RotheSolver_mod.o: RotheSolver_mod.f90
optimizers.o: optimizers.f90 RotheSolver_mod.o lbfgs.o
propagate_Gaussians.o: propagate_Gaussians.f90 optimizers.o ho_kernels.o RotheSolver_mod.o

%.o: %.f90
	$(FC) -c $(OPT) $(WARN) -I$(STD_INC) $< -o $@

%.o: %.f
	$(FC) -c $(OPT)       -I$(STD_INC) $< -o $@

clean:
	rm -f *.o *.mod propagate_Gaussians
.PHONY: all clean
# --------------------------------------------------------------------
