PLAT=intel+nvidia

FF=mpifort
FFLAGS=-O2 -r8 -acc -g -Mcuda

MPI_ROOT=${I_MPI_ROOT}

MKL_LIB=${MKLROOT}/lib/intel64
MKL_INC=${MKLROOT}/mkl/include

MKLROOT=/opt/intel/compilers_and_libraries_2020.4.304/linux/mkl
MKLPATH=${MKLROOT}/lib/intel64
MKLINCLUDE=${MKLROOT}/include

LDFLAGS=-L$(MKLROOT) -L$(MKLROOT)/lib/intel64/ -I$(MKLROOT)/include -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_blas95_lp64.a $(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a -Wl,--end-group -lpthread -lgfortran -lmkl_rt -lcudart -lcublas -Mcuda

