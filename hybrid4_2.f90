!======================================================================!
program hybrid4_2
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Learn from /home/adf10/MODELS/HYBRID10 on CSD3 for how to use
! TRENDY forcings. FW00 is Friend AD, White A. 2000. Eval...
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
use netcdf
use mpi
use shared
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------!

character (len = 250) :: file_name ! Generic filename
integer :: ncid
integer :: varid
integer :: lon_dimid, lat_dimid
integer :: lon_varid, lat_varid
integer, dimension (2) :: dimids_two

!----------------------------------------------------------------------!
integer, parameter :: ntimes      =  1460
integer, parameter :: fillvalue_int = -99999
character (len = 100) :: var_name
character (len =   4) :: char_year
integer :: varid_Cm,varid_Cu,varid_Cn,varid_Cv
integer :: varid_Ca,varid_Cs,varid_Cpa
integer :: varid_Nm,varid_Nu,varid_Nn,varid_Nv
integer :: varid_Na,varid_Ns,varid_Npa
integer :: varid_snmin
integer :: varid_nind
integer :: varid_snow,varid_soilw1,varid_soilw2,varid_soilw3
integer :: varid_ind
integer :: varid_land_index
integer :: varid_cd
integer :: varid_dd
integer :: varid_bgs
integer :: varid_egs
integer :: varid_k_ind
integer :: varid_kspp
integer :: varid_dbh
integer :: varid_hdbh
integer :: varid_cfoliage
integer :: varid_cfiner
integer :: varid_cstore
integer :: varid_nfoliage
integer :: varid_nbswood
integer :: varid_nheart
integer :: varid_navail
integer :: varid_nfiner
integer :: varid_lasa
integer :: varid_nitf
integer :: varid_hbc
integer :: plot_dimid, ind_p_dimid, ind_t_dimid
integer :: plot_varid, ind_p_varid, ind_t_varid
integer :: land_dimid,land_varid
integer, dimension (1) :: dimids_one
integer, dimension (3) :: dimids_three
integer :: ii,jj
integer :: k,k1,k2
integer :: ki
integer :: kp,p
integer :: in,kil
integer :: m
integer :: ksp
integer :: il
integer :: nlayers
!integer :: cd_flag
real, allocatable, dimension (:,:) :: icwtr_qd ! Ice/water fraction
real, allocatable, dimension (:,:) :: isc ! Initial soil C (kg[C] m-2)
real, allocatable, dimension (:,:) :: larea_qd ! Grid-box area (km2)
!----------------------------------------------------------------------!
integer :: nplots_total
real :: isc_total ! Total global soil C (Pg[C])
!real :: wlitterc
!real :: wlittern
!real :: flitterc
!real :: flittern
!real :: rlitterc
!real :: rlittern
!real :: h2o_30
!real :: pleach
!real :: pmf
!real :: puf
!real :: pnr
!real :: pvr
!real :: psa
!real :: temp,C0,N0,Cleach,Nmin,Usoil
!real :: sresp ! Soil respiration (kg[C] m-1 day-1)
!real, dimension (nspp) :: mbiomsp
real :: fd
real :: kf,kSWf
real :: fracs
real :: frac
real :: lati,loni,rltdli,ddreqi
real :: phi,rn,delta,tn,gn
!----------------------------------------------------------------------!

allocate (icwtr_qd (2*nlon,2*nlat))
allocate (larea_qd (2*nlon,2*nlat))

!----------------------------------------------------------------------!

call initl

!----------------------------------------------------------------------!
! Read in ice/water fractions for each grid-box, and areas (km2).
! Water is ocean and freshwater bodies from map.
!----------------------------------------------------------------------!
file_name = '/home/adf10/rds/rds-mb425-geogscratch/adf10/FORCINGS/&
&LUH2_new/staticData_quarterdeg.nc'
write (*, *) 'Reading from ', trim (file_name)
call check (nf90_open (trim (file_name), nf90_nowrite, ncid))
varid = 7
call check (nf90_get_var (ncid, varid, icwtr_qd))
varid = 9
call check (nf90_get_var (ncid, varid, larea_qd))
call check (nf90_close (ncid))
!----------------------------------------------------------------------!
! Input file is 1/4 degree, so gridded to 1/2 degree. Need to invert.
! Assume no missing values.
!----------------------------------------------------------------------!
jj = 1
do j = 1, nlat
 ii = 1
 do i = 1, nlon
  icwtr (i, nlat-j+1) = sum (icwtr_qd (ii:ii+1, jj:jj+1)) / 4.0
  larea (i, nlat-j+1) = sum (larea_qd (ii:ii+1, jj:jj+1))
  ii = ii + 2
 end do
 jj = jj + 2
end do
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Read phenology parameters for each grid-box.
! Should re-create these based on TRENDY grid and climate.
!----------------------------------------------------------------------!
rltdl (:,:) = zero
ddreq (:,:) = zero
open (10,file='spin_up.phen',status='old')
do k = 1, 67420
 read (10,*) lati,loni,rltdli,ddreqi
 i = nint ((loni - 0.25 + 180.0) * 2.0) + 1
 j = nint ((lati - 0.25 +  90.0) * 2.0) + 1
 rltdl (i,j) = rltdli
 ddreq (i,j) = ddreqi
end do
close (10)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Daylengths (s). Try to only do local when local.
!----------------------------------------------------------------------!
do j = 1, nlat
 !...phi is latitude (radians)
 phi = pi * abs (lat (j)) / 180.0
 ! Calculate daylength for each yearday.
 do jdl = 1, nd
  ! rn is climatological day number
  rn = float (jdl)
  ! delta is solar declination (radians), depends on time of year.
  ! Calculation taken from Spitters et al., 1986 (Eqn. 16, p. 226).
  delta = asin (-1.0 * sin (23.45 * pi / 180.0) * &
          cos (360.0 * pi / 180.0 * (rn + 10.0) / float (nd)))
  ! Daylength calculation based on France & Thornley (1984).
  ! tn is used in daylength calculation.
  tn = - 1.0 * tan (phi) * tan (delta)
  ! Check if daylength between 0 and 24 hr.
  if ((tn >= -1.0) .and. (tn <= 1.0)) then
   ! gn is decimal parts of a day.
   gn = (180.0 / pi) * 2.0 * acos (tn) / 360.0
  else
   if (tn > one) then
    gn = zero
   else
    gn = one
   end if
  end if
  dl (j,jdl) = gn * sday
 end do ! jdl
 ! Invert if in southern hemisphere.
 if (lat (j) < zero) then
  do kday = 0, nd
   dl (j,kday) = sday - dl (j,kday)
  end do
 end if
 !dl (j,0) = dl (j,nd) ! needed?
end do ! j
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
allocate (ball (nind_total))
allocate (foff (nind_total))
allocate (fon  (nind_total))
NPP_grid     (:,:) = fillvalue
Cv_grid      (:,:) = fillvalue
Cs_grid      (:,:) = fillvalue
Cv_C3GR_grid (:,:) = fillvalue
Cv_C4GR_grid (:,:) = fillvalue
Cv_BRCD_grid (:,:) = fillvalue
Cv_NLEV_grid (:,:) = fillvalue
!----------------------------------------------------------------------!

summer = .FALSE.

!----------------------------------------------------------------------!
! Loop through gridboxes and integrate state variables at land points.
! Climate files start at 1901, run through 2019.
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
do kyr = syr, eyr

 !---------------------------------------------------------------------!
 GPP_global = zero
 NPP_global = zero
 Cv_global  = zero
 Cs_global  = zero
 !---------------------------------------------------------------------!

 !---------------------------------------------------------------------!
 ! C balance of lowest foliage layer (kg[C] ind-1 yr-1).
 !---------------------------------------------------------------------!
 ball (:) = zero
 !---------------------------------------------------------------------!
 
 !---------------------------------------------------------------------!
 ! Flag for cold-deciduous diagnostic.
 !---------------------------------------------------------------------!
 !cd_flag = 0
 !---------------------------------------------------------------------!
 
 !---------------------------------------------------------------------!
 ! Flag for dry-deciduous diagnostic.
 !---------------------------------------------------------------------!
 dd_flag (:) = 0
 !---------------------------------------------------------------------!
 
 !---------------------------------------------------------------------!
 ! Flags for deciduousness in year.
 !---------------------------------------------------------------------!
 foff (:) = 0
 fon  (:) = 0
 !---------------------------------------------------------------------!

 !---------------------------------------------------------------------!
 write (char_year, '(i4)') kyr
 !---------------------------------------------------------------------!

 !---------------------------------------------------------------------!
 ! Read global temperature fields for year kyr into tmp (K).
 !---------------------------------------------------------------------!
 var_name = 'tmp'
 file_name = '/rds/user/adf10/rds-mb425-geogscratch/adf10/FORCINGS/&
 &CRUJRA_2.1/CRUJRA2020/'//TRIM(var_name)//'/crujra.v2.1.5d.'&
 &//TRIM(var_name)//'.'//char_year//'.365d.noc.nc'
 write (*,*) 'Reading from ',trim(file_name)
 call check (nf90_open (trim (file_name), nf90_nowrite, ncid))
 ! Read temperatures (K).
 varid = 4
 call check (nf90_get_var (ncid, varid, tmp))
 call check (nf90_close (ncid))
 !---------------------------------------------------------------------!

 !---------------------------------------------------------------------!
 ! Read precipitation fields for year kyr intp pre (mm 6hr-1).
 !---------------------------------------------------------------------!
 var_name = 'pre'
 file_name = '/rds/user/adf10/rds-mb425-geogscratch/adf10/FORCINGS/&
 &CRUJRA_2.1/CRUJRA2020/'//TRIM(var_name)//'/crujra.v2.1.5d.'&
 &//TRIM(var_name)//'.'//char_year//'.365d.noc.nc'
 write (*,*) 'Reading from ',trim(file_name)
 call check (nf90_open (trim (file_name), nf90_nowrite, ncid))
 ! Read precipitation (mm 6hr-1).
 varid = 4
 call check (nf90_get_var (ncid, varid, pre))
 call check (nf90_close (ncid))
 !---------------------------------------------------------------------!
 
 !---------------------------------------------------------------------!
 ! Read global downward solar radiation flux fields for year kyr into
 ! dswrf (J m-2).
 !---------------------------------------------------------------------!
 var_name = 'dswrf'
 file_name = '/rds/user/adf10/rds-mb425-geogscratch/adf10/FORCINGS/&
 &CRUJRA_2.1/CRUJRA2020/'//TRIM(var_name)//'/crujra.v2.1.5d.'&
 &//TRIM(var_name)//'.'//char_year//'.365d.noc.nc'
 write (*,*) 'Reading from ',trim(file_name)
 call check (nf90_open (trim (file_name), nf90_nowrite, ncid))
 ! Read solar radation (J m-2).
 varid = 4
 call check (nf90_get_var (ncid, varid, dswrf))
 call check (nf90_close (ncid))
 !---------------------------------------------------------------------!
 
 !---------------------------------------------------------------------!
 ! Read specific humidity fields for year kyr into spfh (kg kg-1).
 !---------------------------------------------------------------------!
 var_name = 'spfh'
 file_name = '/rds/user/adf10/rds-mb425-geogscratch/adf10/FORCINGS/&
 &CRUJRA_2.1/CRUJRA2020/'//TRIM(var_name)//'/crujra.v2.1.5d.'&
 &//TRIM(var_name)//'.'//char_year//'.365d.noc.nc'
 write (*,*) 'Reading from ',trim(file_name)
 call check (nf90_open (trim (file_name), nf90_nowrite, ncid))
 ! Read specific humidities (kg kg-1).
 varid = 4
 call check (nf90_get_var (ncid, varid, spfh))
 call check (nf90_close (ncid))
 !---------------------------------------------------------------------!

 !---------------------------------------------------------------------!
 ! Read pressure fields for year kyr into pres (Pa).
 !---------------------------------------------------------------------!
 var_name = 'pres'
 file_name = '/rds/user/adf10/rds-mb425-geogscratch/adf10/FORCINGS/&
 &CRUJRA_2.1/CRUJRA2020/'//TRIM(var_name)//'/crujra.v2.1.5d.'&
 &//TRIM(var_name)//'.'//char_year//'.365d.noc.nc'
 write (*,*) 'Reading from ',trim(file_name)
 call check (nf90_open (trim (file_name), nf90_nowrite, ncid))
 ! Read pressures (Pa).
 varid = 4
 call check (nf90_get_var (ncid, varid, pres))
 call check (nf90_close (ncid))
 !---------------------------------------------------------------------!
 
 !---------------------------------------------------------------------!
 ! CO2 stomatal response (assuming closing above 80 Pa). Should give
 ! curve based on Fig. 1 of M&G. kyr must be 1700-1900.
 !---------------------------------------------------------------------!
 fc = one - 0.008333 * cao (kyr-1699)
 fc = max (0.33, fc)
 !---------------------------------------------------------------------!

 !---------------------------------------------------------------------!
 nland = 0 ! Just for diagnostic.
 !---------------------------------------------------------------------!
 do j = j1, j2 ! Latitude count (starts at SP).
  do i = i1, i2 ! Longitude count (starts at IDL and works eastwards).
   !-------------------------------------------------------------------!
   ! Land?
   !-------------------------------------------------------------------!
   if (tmp (i,j,1) /= fillvalue) then
    nland = nland + 1 ! For diagnostic.
    
    !------------------------------------------------------------------!
    ! Set initial annual values to zero.
    !------------------------------------------------------------------!
    ! Annual growth respiration (kg[C] m-2 yr-1).
    !------------------------------------------------------------------!
    rgs = zero
    !------------------------------------------------------------------!
    ! Degree-days.
    !------------------------------------------------------------------!
    dd (i,j) = zero !****adf OK over SH?
    !------------------------------------------------------------------!
    ! Chilling days.
    !------------------------------------------------------------------!
    cd (i,j) = zero !****adf OK over SH?
    !------------------------------------------------------------------!
    ! Flag for tree foliage area change.
    !------------------------------------------------------------------!
    phenf (:) = 1
    !------------------------------------------------------------------!
    ! CD flag.
    !------------------------------------------------------------------!
    ! cd_flag = 0
    !------------------------------------------------------------------!
    ! Flag for dry deciduousness.
    !------------------------------------------------------------------!
    dd_flag (:) = 0
    !------------------------------------------------------------------!
    ! Reset fturn as may have been altered in phen for litter calc. in
    ! cflux.
    !------------------------------------------------------------------!
    do kp = 1, nplots
     do ksp = 1, nspp
      fturn_plot (kp,ksp) = fturn_save (ksp)
     end do ! ksp
     !-----------------------------------------------------------------!
     do ki = 1, nind (i,j,kp)
      k = k_ind (land_index(i,j),kp,ki)
      ball (k) = zero
     end do ! ki
    end do ! kp
    !------------------------------------------------------------------!
    foff (:) = 0
    fon  (:) = 0
    !------------------------------------------------------------------!
    NPP_grid (i,j) = zero
    !------------------------------------------------------------------!
    if (local) mnppsp (:) = zero
    !------------------------------------------------------------------!
    ! End set zero.
    !------------------------------------------------------------------!
    
    !------------------------------------------------------------------!
    kday = 1
    !------------------------------------------------------------------!
    do it = 1, ntimes ! Loop timepoints in year.
     
     !-----------------------------------------------------------------!
     if (mod (it+3, 4) == zero) then ! Start of day.
      kday = (it + 3) / 4
      ! Set degree-days and chilling gdays to zero (oC).
      if (lat (j) >= zero) then
       if (kday <=  32) dd (i,j) = zero
       if (kday == 305) cd (i,j) = 0
      else
       if (kday <= 182) dd (i,j) = zero
       if (kday ==  60) cd (i,j) = 0
      end if
      !----------------------------------------------------------------!
     end if ! mod start of day.
     !-----------------------------------------------------------------!
     
     call getclmn
     call hydrology
     call carbon
     
     ! End of day?
     if (mod (it, 4) == zero) then
      
      pptod = sum (pre (i, j, it-3:it)) / 1000.0
      call nitrogen
      call grass
      call phenology
      call soil
      call radiation
      
     end if ! (mod (it, 4) == zero)
    end do ! it = 1, ntimes
    
    call allocate
    call mortal
    call regen
    
   end if ! tmp (i,j,k) /= fillvalue
   
  end do ! i = i1, i2
 end do ! j = j1, j2
 
 call annual_diagnostics
 
end do ! kyr = syr, eyr
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
call run_finish
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
contains
 subroutine check (status)

 integer, intent (in) :: status
 if (status /= nf90_noerr) then
  print *, trim (nf90_strerror(status))
  stop  "Stopped"
 end if
 end subroutine check
!----------------------------------------------------------------------!

end program hybrid4_2
!======================================================================!
