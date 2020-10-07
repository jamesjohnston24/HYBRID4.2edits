F90 = mpif90 -I/usr/local/Cluster-Apps/intel/2017.4/compilers_and_libraries_2017.4.196/linux/mpi/intel64/include

FFLAGS = -O3
LDFLAGS = -lnetcdff -lnetcdf

SRC = \
	hybrid4_2.f90

OBJ = $(SRC:.f90=.o)

hybrid4_2.exe : $(OBJ)
	$(F90) $(FFLAGS) -o hybrid4_2.exe $(OBJ) $(LDFLAGS)

# Main routine.
hybrid4_2.o : hybrid4_2.f90
	$(F90) $(FFLAGS) -c hybrid4_2.f90
