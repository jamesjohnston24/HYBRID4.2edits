!======================================================================!
program hybrid4_2
!----------------------------------------------------------------------!
! just to make change
!----------------------------------------------------------------------!
! Learn from /home/adf10/MODELS/HYBRID10 on CSD3 for how to use
! TRENDY forcings. FW00 is Friend AD, White A. 2000. Eval...
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Enable Network Common Data Forum (netCDF) libraries.
!----------------------------------------------------------------------!
use netcdf
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Enable library Message Passing Interface (MPI) library.
!----------------------------------------------------------------------!
use mpi
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Declare all parameters and variables shared between subroutines.
!----------------------------------------------------------------------!
use shared ! checked
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Initialise and allocate global arrays.
! Read run control parameters from driver.in.
! Open diagnostic output files.
! Read annual CO2 values.
! Set derived GPT-level parameters.
! Read initial state from restart file if requested, otherwise
! Initialise all state variables to standard values.
! Set variables derived from state variables.
! Read grid-box ice/water fractions and areas.
! Read phenology parameters.
! Set daylengths.
!----------------------------------------------------------------------!
call initl
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Loop through gridboxes and integrate state variables at land points.
! Climate files start at 1901, run through 2019.
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
do kyr = syr, eyr

 !---------------------------------------------------------------------!
 ! Zero annual diagnostics and variables.
 ! Read global annual climate.
 ! Set stomatal CO2 scalar, fc.
 !---------------------------------------------------------------------!
 !****adf
 if(kyr==syr)call init_yr
 !---------------------------------------------------------------------!

 !---------------------------------------------------------------------!
 ! Loop over all half-degree boxes on Earth.
 !---------------------------------------------------------------------!
 do j = j1, j2 ! Latitude count (starts at SP).
  do i = i1, i2 ! Longitude count (starts at IDL and works eastwards).
   !-------------------------------------------------------------------!
   ! Climate for this box (and hence some land)?
   !-------------------------------------------------------------------!
   if (tmp (i,j,1) /= fillvalue) then
   !-------------------------------------------------------------------!
   
    !------------------------------------------------------------------!
    ! Set annual phenology and diagnostic values of grid-box to zero.
    !------------------------------------------------------------------!
    call init_grid_box
    !------------------------------------------------------------------!
    
    !------------------------------------------------------------------!
    do it = 1, ntimes ! Loop (1460 * 6-hr) timepoints in year.
    !------------------------------------------------------------------!
     
     !-----------------------------------------------------------------!
     ! Compute kday if start of day (day of year).
     !-----------------------------------------------------------------!
     if (mod (it+3, 4) == zero) then ! Start of day.
      kday = (it + 3) / 4
     end if ! mod start of day.
     !-----------------------------------------------------------------!
     ! Set degree-days and chilling-days to zero if Feb. 1st and Nov.
     ! 1st resp. in NH, if July 1st and March 1st resp. in SH.
     !-----------------------------------------------------------------!
     if (lat (j) >= zero) then
      if (kday <=  32) dd (i,j) = zero
      if (kday == 305) cd (i,j) = 0
     else
      !****adf should be 215 and 123? changed from 182 and 60
      if (kday <= 215) dd (i,j) = zero
      if (kday == 123) cd (i,j) = 0
     end if
     !----------------------------------------------------------------!
     
     !----------------------------------------------------------------!
     ! Assign climate variables from global arrays.
     ! Compute derived climate values.
     ! Accumulate degree-days.
     ! Compute temperature factors.
     !----------------------------------------------------------------!
     call getclmn
     !----------------------------------------------------------------!
     
     !----------------------------------------------------------------!
     ! Compute N uptake temperature effect (why here?).
     ! Compute effect soil water potential for each individual.
     ! Compute stomatal conductance for each individual.
     ! Compute canopy temperature for each plot.
     ! Update soil water in each layer and snowpack.
     ! Compute proportion of water lost as outflow.
     !----------------------------------------------------------------!
     call hydrology
     !----------------------------------------------------------------!
     
     !----------------------------------------------------------------!
     ! Compute carbon balance of each individual.
     ! Compute carbon balance of lowest foliage layers.
     !----------------------------------------------------------------!
     call carbon
     !----------------------------------------------------------------!
     
     !----------------------------------------------------------------!
     if (mod (it, 4) == zero) then ! End of day.
     !----------------------------------------------------------------!
      
      !---------------------------------------------------------------!
      ! Compute N uptake by each individual.
      ! Subtract N uptake from soil mineral N.
      !---------------------------------------------------------------!
      call nitrogen
      !---------------------------------------------------------------!
      
      !---------------------------------------------------------------!
      ! Update C and N pools in grass layer.
      !---------------------------------------------------------------!
      call grass
      !---------------------------------------------------------------!
      
      !---------------------------------------------------------------!
      ! Accumulate chilling-days.
      ! Compute if CD trees have foliage on.
      ! Compute daily litter and compensation growth from trees.
      ! DD trees can lose leaves.
      ! Compute drought and cold damage.
      ! Compute new plot LAI.
      ! Compute new radiation in lowest layer and related
      ! physiological variables if no tree foliage area change.
      !---------------------------------------------------------------!
      call phenology
      !---------------------------------------------------------------!
      
      !---------------------------------------------------------------!
      ! Update soil C and N pools.
      ! Compute water holding capacity in each layer.
      !---------------------------------------------------------------!
      call soil
      !---------------------------------------------------------------!
      
      !---------------------------------------------------------------!
      ! If tree foliage area changed in phenology.f90, compute new
      ! radiation profiles and related physiological variables for all
      ! individuals.
      !---------------------------------------------------------------!
      call radiation
      !---------------------------------------------------------------!
      
     !----------------------------------------------------------------!
     end if ! (mod (it, 4) == zero)
     !----------------------------------------------------------------!
     
    !-----------------------------------------------------------------!
    end do ! it = 1, ntimes
    !-----------------------------------------------------------------!
    
    !-----------------------------------------------------------------!
    ! Allocate tree C and N.
    ! Compute new plot LAI.
    !-----------------------------------------------------------------!
    call allocate
    !-----------------------------------------------------------------!
    
    !-----------------------------------------------------------------!
    ! Assign live = 0 if tree dead and move C and N to litter pools.
    !-----------------------------------------------------------------!
    call mortal
    !-----------------------------------------------------------------!
    
    !-----------------------------------------------------------------!
    ! Plant trees if spaces available.
    !-----------------------------------------------------------------!
    call regen
    !-----------------------------------------------------------------!
    
   !------------------------------------------------------------------!
   end if ! tmp (i,j,k) /= fillvalue
   !------------------------------------------------------------------!
   
  end do ! i = i1, i2 ! Next longitude.
 end do ! j = j1, j2 ! Next latitude.
 ! End of Earth loop.
 !---------------------------------------------------------------------!
 
 !---------------------------------------------------------------------!
 ! Compute annual grid diagnostics and global totals.
 ! Report all if local simulation.
 ! Report global totals to unit 20 (annual_global.dat) if global
 ! simulation.
 !---------------------------------------------------------------------!
 call annual_diagnostics
 !---------------------------------------------------------------------!
 
end do ! kyr = syr, eyr
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! If not local simulation output global grid diagnostics.
! If required, write state variables to restart file.
! Close unit 20 (annual_global.dat).
!----------------------------------------------------------------------!
call run_finish
!----------------------------------------------------------------------!

end program hybrid4_2
!======================================================================!
