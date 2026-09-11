# fortran compiler
FC = mpif90
FFLAGS := -O3 -ffixed-line-length-none -Wall
BIG :=  -mcmodel=large

interpolate_MPI.x : config_param.o mpi_utilities.o interpolate_MPI.o 
	$(FC) $(FFLAGS) $(BIG) config_param.o mpi_utilities.o interpolate_MPI.o -o interpolate_MPI.x

interpolate_MPI.o : interpolate_MPI.f90 
	$(FC) $(FFLAGS) $(BIG) $(INCLUDE) -c interpolate_MPI.f90 

config_param.o : config_param.f90 
	$(FC) $(FFLAGS) -c config_param.f90 -o config_param.o 

mpi_utilities.o : mpi_utilities.f90 
	$(FC) $(FFLAGS) -c mpi_utilities.f90 -o mpi_utilities.o

.PHONY: clean
clean:
	rm -f *.o *.x *.mod *.out ._*
