F90 = mpif90 -I/usr/local/Cluster-Apps/intel/2017.4/compilers_and_libraries_2017.4.196/linux/mpi/intel64/include

# FFLAGS = -O3 -check -traceback
FFLAGS = -O3
LDFLAGS = -lnetcdff -lnetcdf

SRC = \
	hybrid4_2.f90		\
	initl.f90		\
	init_yr.f90		\
	init_grid_box.f90	\
	getclmn.f90		\
	hydrology.f90		\
	carbon.f90		\
	nitrogen.f90		\
	grass.f90		\
	phenology.f90		\
	soil.f90		\
	allocate.f90		\
	mortal.f90		\
	regen.f90		\
	radiation.f90		\
	annual_diagnostics.f90	\
	run_finish.f90		\
	shared.f90

OBJ = $(SRC:.f90=.o)

hybrid4_2.exe : $(OBJ)
	$(F90) $(FFLAGS) -o hybrid4_2.exe $(OBJ) $(LDFLAGS)

# Main routine.
hybrid4_2.o : shared.o initl.o init_yr.o init_grid_box.o getclmn.o hydrology.o carbon.o nitrogen.o grass.o phenology.o soil.o allocate.o mortal.o regen.o radiation.o annual_diagnostics.o run_finish.o hybrid4_2.f90
	$(F90) $(FFLAGS) -c hybrid4_2.f90
	
# Subroutines.
initl.o : shared.o initl.f90
	$(F90) $(FFLAGS) -c initl.f90

init_yr.o : shared.o init_yr.f90
	$(F90) $(FFLAGS) -c init_yr.f90

init_grid_box.o : shared.o init_grid_box.f90
	$(F90) $(FFLAGS) -c init_grid_box.f90
	
getclmn.o : shared.o getclmn.f90
	$(F90) $(FFLAGS) -c getclmn.f90
	
hydrology.o : shared.o hydrology.f90
	$(F90) $(FFLAGS) -c hydrology.f90
	
carbon.o : shared.o carbon.f90
	$(F90) $(FFLAGS) -c carbon.f90
	
nitrogen.o : shared.o nitrogen.f90
	$(F90) $(FFLAGS) -c nitrogen.f90
	
grass.o : shared.o grass.f90
	$(F90) $(FFLAGS) -c grass.f90
	
phenology.o : shared.o phenology.f90
	$(F90) $(FFLAGS) -c phenology.f90
	
soil.o : shared.o soil.f90
	$(F90) $(FFLAGS) -c soil.f90
	
allocate.o : shared.o allocate.f90
	$(F90) $(FFLAGS) -c allocate.f90
	
mortal.o : shared.o mortal.f90
	$(F90) $(FFLAGS) -c mortal.f90
	
regen.o : shared.o regen.f90
	$(F90) $(FFLAGS) -c regen.f90

radiation.o : shared.o radiation.f90
	$(F90) $(FFLAGS) -c radiation.f90
	
annual_diagnostics.o : shared.o annual_diagnostics.f90
	$(F90) $(FFLAGS) -c annual_diagnostics.f90
	
run_finish.o : shared.o run_finish.f90
	$(F90) $(FFLAGS) -c run_finish.f90
	
# Modules.
shared.o : shared.f90
	$(F90) $(FFLAGS) -c shared.f90
