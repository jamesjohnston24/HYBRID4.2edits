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
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Options if run locally.
!----------------------------------------------------------------------!
logical :: local  = .TRUE. ! Run only local site?
!logical :: local  = .FALSE. ! Run only local site?
logical :: wrclim = .FALSE. ! Write local climate?
!----------------------------------------------------------------------!
! Uggla site.
!----------------------------------------------------------------------!
!real, parameter :: lon_w = 19.0 + 46.0 / 60.0
!real, parameter :: lat_w = 64.0 + 21.0 / 60.0
!----------------------------------------------------------------------!
! Grudd site.
!----------------------------------------------------------------------!
!real, parameter :: lon_w = 19.75
!real, parameter :: lat_w = 68.25
!----------------------------------------------------------------------!
! Cambridge.
!----------------------------------------------------------------------!
real, parameter :: lon_w =  0.108
real, parameter :: lat_w = 52.198
!----------------------------------------------------------------------!
! Indices of local site, if used.
!----------------------------------------------------------------------!
integer :: i_w
integer :: j_w
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Noise on initial soil C across plots?
logical, parameter :: ran_s = .TRUE.
! Start from restart file?
!logical, parameter :: rsf_in = .FALSE.
logical, parameter :: rsf_in = .TRUE.
! Write to restart file?
logical, parameter :: rsf_out = .FALSE.
!logical, parameter :: rsf_out = .TRUE.
integer, parameter :: nlon = 720
integer, parameter :: nlat = 360
integer, parameter :: ntimes = 1460
integer, parameter :: nplots = 10
real, parameter :: fillvalue = 1.0e20
real, parameter :: tf = 273.15
real, parameter :: eps = 1.0e-8
real, parameter :: zero = 0.0
real, parameter :: one = 1.0
! Mol. weight of dry air (kg [air] mol-1).
real, parameter :: m_air = 28.9647
! Mol. weight of water (kg [water] mol-1).
real, parameter :: m_water = 18.01528
real, parameter :: pi = 3.14159265359
real, parameter :: iscv = 0.2 ! Range of initial soil C
! N:C ratio of surface structural litter.
real, parameter :: vu =  1.0 / 150.0
! N:C ratio of soil structural litter.
real, parameter :: vv =  1.0 / 150.0
! N:C ratio of active soil organic matter pool.
real, parameter :: va =  1.0 /  15.0
! N:C ratio of slow soil organic matter pool.
real, parameter :: vs =  1.0 /  20.0
! N:C ratio of slow soil passive organic matter pool.
real, parameter :: vpa = 1.0 / 10.0
real, parameter :: area = 200.0 ! Plot area (m2)
character (len = 250) :: file_name ! Generic filename
character (len = 100) :: var_name
character (len =     4) :: char_year
integer :: kyr
integer :: it
integer :: syr
integer :: eyr
integer :: varid
integer :: varid_Cm,varid_Cu,varid_Cn,varid_Cv,varid_Ca,varid_Cs,varid_Cpa
integer :: varid_Nm,varid_Nu,varid_Nn,varid_Nv,varid_Na,varid_Ns,varid_Npa
integer :: varid_snmin
integer :: varid_nind
integer :: varid_soilw1,varid_soilw2,varid_soilw3
integer :: ncid
integer :: lon_dimid
integer :: lat_dimid
integer :: plot_dimid
integer :: lon_varid
integer :: lat_varid
integer :: plot_varid
integer, dimension (2) :: dimids_two
integer, dimension (3) :: dimids_three
integer :: i,ii,i1,i2
integer :: j,jj,j1,j2
integer :: k
integer :: kp
integer :: ki
integer :: nland
integer, allocatable, dimension (:) :: plot ! No. of plot
real, allocatable, dimension (:,:) :: icwtr_qd ! Ice/water fraction
real, allocatable, dimension (:,:) :: icwtr ! Ice/water fraction
real, allocatable, dimension (:,:,:) :: tmp ! K
real, allocatable, dimension (:,:,:) :: pre ! mm 6hr-1
real, allocatable, dimension (:,:,:) :: spfh ! kg kg-1
real, allocatable, dimension (:,:,:) :: pres ! Pa
real, allocatable, dimension (:) :: lat ! degrees north
real, allocatable, dimension (:) :: lon ! degrees east
real, allocatable, dimension (:,:) :: isc ! Initial soil C (kg[C] m-2)
real, allocatable, dimension (:,:) :: isc_grid ! Initial soil C (kg[C] m-2)
real, allocatable, dimension (:,:) :: larea_qd ! Grid-box area (km2)
real, allocatable, dimension (:,:) :: larea ! Grid-box area (km2)
real, allocatable, dimension (:,:,:) :: Cm
real, allocatable, dimension (:,:,:) :: Cu
real, allocatable, dimension (:,:,:) :: Cn
real, allocatable, dimension (:,:,:) :: Cv
real, allocatable, dimension (:,:,:) :: Ca
real, allocatable, dimension (:,:,:) :: Cs
real, allocatable, dimension (:,:,:) :: Cpa
real, allocatable, dimension (:,:,:) :: Nu
real, allocatable, dimension (:,:,:) :: Nm
real, allocatable, dimension (:,:,:) :: Nv
real, allocatable, dimension (:,:,:) :: Nn
real, allocatable, dimension (:,:,:) :: Na
real, allocatable, dimension (:,:,:) :: Ns
real, allocatable, dimension (:,:,:) :: Npa
real, allocatable, dimension (:,:,:) :: snmin
integer, allocatable, dimension (:,:,:) :: nind
real, allocatable, dimension (:,:,:) :: soilw1
real, allocatable, dimension (:,:,:) :: soilw2
real, allocatable, dimension (:,:,:) :: soilw3
real, allocatable, dimension (:) :: wlittc
real, allocatable, dimension (:) :: wlittn
real, allocatable, dimension (:) :: flittc
real, allocatable, dimension (:) :: flittn
real, allocatable, dimension (:) :: rlittc
real, allocatable, dimension (:) :: rlittn
real :: t ! Temperature (degree C)
real :: t_d ! Dew-point (degree C)
real :: eo ! Penman evaporation from lake (mm day-1)
real :: pt ! Annual precipiation (m yr-1)
real :: vap ! Vapour pressure (Pa)
real :: e ! Vapour pressure (mbar)
real :: isc_total ! Total global soil C (Pg[C])
real :: ran ! Random number (0-1)
real :: Noise ! Noise on initial soil C (ratio)
real :: vm ! N:C ratio of surface metabolic litter
real :: vn ! N:C ratio of soil metabolic litter
real :: tsoil ! Soil temperature (oC)
real :: et_soil ! Soil decomposition temperature modifier (fraction)
real :: em_soil ! Soil decomposition water modifer (fraction)
real :: ev ! Overall soil decomposition modifier (fraction)
real :: wlitterc
real :: wlittern
real :: flitterc
real :: flittern
real :: rlitterc
real :: rlittern
real :: swct ! Soil water holding capacity (m)
real :: wfps ! Water-filled pore space (%)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
allocate (lon  (nlon))
allocate (lat  (nlat))
allocate (plot (nplots))
allocate (icwtr_qd (2*nlon,2*nlat))
allocate (icwtr (nlon,nlat))
! Initial total soil C on non-ice/water fraction (kg[C] m-2)
allocate (isc  (nlon,nlat))
! Initial total soil C grid-box mean (kg[C] m-2)
allocate (isc_grid  (nlon,nlat))
allocate (larea_qd (2*nlon,2*nlat))
allocate (larea (nlon,nlat))
allocate (tmp  (nlon,nlat,ntimes))
allocate (pre  (nlon,nlat,ntimes))
allocate (spfh (nlon,nlat,ntimes))
allocate (pres (nlon,nlat,ntimes))
allocate (Cm   (nlon,nlat,nplots))
allocate (Cu   (nlon,nlat,nplots))
allocate (Cn   (nlon,nlat,nplots))
allocate (Cv   (nlon,nlat,nplots))
allocate (Ca   (nlon,nlat,nplots))
allocate (Cs   (nlon,nlat,nplots))
allocate (Cpa  (nlon,nlat,nplots))
allocate (Nu   (nlon,nlat,nplots))
allocate (Nm   (nlon,nlat,nplots))
allocate (Nv   (nlon,nlat,nplots))
allocate (Nn   (nlon,nlat,nplots))
allocate (Na   (nlon,nlat,nplots))
allocate (Ns   (nlon,nlat,nplots))
allocate (Npa  (nlon,nlat,nplots))
allocate (snmin (nlon,nlat,nplots))
allocate (nind (nlon,nlat,nplots))
allocate (soilw1 (nlon,nlat,nplots))
allocate (soilw2 (nlon,nlat,nplots))
allocate (soilw3 (nlon,nlat,nplots))
allocate (wlittc (nplots))
allocate (wlittn (nplots))
allocate (flittc (nplots))
allocate (flittn (nplots))
allocate (rlittc (nplots))
allocate (rlittn (nplots))
!----------------------------------------------------------------------!
Cn  (:,:,:) = zero
Cu  (:,:,:) = zero
Cm  (:,:,:) = zero
Cv  (:,:,:) = zero
Ca  (:,:,:) = zero
Cs  (:,:,:) = zero
Cpa (:,:,:) = zero
Nn  (:,:,:) = zero
Nu  (:,:,:) = zero
Nm  (:,:,:) = zero
Nv  (:,:,:) = zero
Na  (:,:,:) = zero
Ns  (:,:,:) = zero
Npa (:,:,:) = zero
snmin (:,:,:) = zero
nind (:,:,:) = 0
soilw1 (:,:,:) = fillvalue
soilw2 (:,:,:) = fillvalue
soilw3 (:,:,:) = fillvalue
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
if (local) then
 !---------------------------------------------------------------------!
 ! Gridbox indices for local site.
 !---------------------------------------------------------------------!
 i_w = nint ((719.0 / 2.0) * (1.0 + lon_w / 179.75)) + 1
 j_w = nint ((359.0 / 2.0) * (1.0 + lat_w /   89.75)) + 1
 write(*,*) 'i_w j_w ',i_w,j_w
 if (lat_w == (64.0 + 21.0 / 60.0)) open (20,file='uggla.clm',&
  status='unknown')
 if (lat_w == 68.25) open (20,file='grudd.clm',status='unknown')
 !---------------------------------------------------------------------!
 i1 = i_w
 i2 = i_w
 j1 = j_w
 j2 = j_w
else
 i1 = 1
 i2 = nlon
 j1 = 1
 j2 = nlat
end if
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
if (rsf_in) then
 file_name = '/home/adf10/rds/rds-mb425-geogscratch/adf10/RSF/rsf_in.nc'
 write (*,*) 'Reading from ',trim(file_name)
 call check (nf90_open (trim (file_name), nf90_nowrite, ncid))
 varid = 1
 call check (nf90_get_var (ncid, varid, lon))
 varid = 2
 call check (nf90_get_var (ncid, varid, lat))
 varid = 3
 call check (nf90_get_var (ncid, varid, plot))
 varid = 4
 call check (nf90_get_var (ncid, varid, Cm))
 varid = 5
 call check (nf90_get_var (ncid, varid, Cu))
 varid = 6
 call check (nf90_get_var (ncid, varid, Cn))
 varid = 7
 call check (nf90_get_var (ncid, varid, Cv))
 varid = 8
 call check (nf90_get_var (ncid, varid, Ca))
 varid = 9
 call check (nf90_get_var (ncid, varid, Cs))
 varid = 10
 call check (nf90_get_var (ncid, varid, Cpa))
 varid = 11
 call check (nf90_get_var (ncid, varid, Nm))
 varid = 12
 call check (nf90_get_var (ncid, varid, Nu))
 varid = 13
 call check (nf90_get_var (ncid, varid, Nn))
 varid = 14
 call check (nf90_get_var (ncid, varid, Nv))
 varid = 15
 call check (nf90_get_var (ncid, varid, Na))
 varid = 16
 call check (nf90_get_var (ncid, varid, Ns))
 varid = 17
 call check (nf90_get_var (ncid, varid, Npa))
 varid = 18
 call check (nf90_get_var (ncid, varid, snmin))
 varid = 19
 call check (nf90_get_var (ncid, varid, nind))
 call check (nf90_close (ncid))
else
 do kp = 1, nplots
  plot (kp) = kp
 end do
end if
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
syr = 1901
eyr = 1901
!eyr = 2019
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
if (.NOT. (rsf_in)) then ! Set soil C, N, water;  plot nind.
do kyr = syr, syr

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
 ! If first year, get lats (start at SP; degrees north) and lons.
 if (kyr == syr) then
  varid = 2
  call check (nf90_get_var (ncid, varid, lon))
  varid = 3
  call check (nf90_get_var (ncid, varid, lat))
 end if
 ! Read temperatures (K).
 varid = 4
 call check (nf90_get_var (ncid, varid, tmp))
 call check (nf90_close (ncid))
 !---------------------------------------------------------------------!
 ! Read precipitation fields for year kyr intp pre (mm 6hr-1.
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
 ! Read specific humidity fields for year kyr intp spfh (kg kg-1).
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
 ! Read pressure fields for year kyr intp pres (Pa).
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
 if ((local) .and. (wrclim)) then
  !--------------------------------------------------------------------!
  ! Climate output of local climate.
  !--------------------------------------------------------------------!
  do it = 1, ntimes
   write (20,'(2i5,f12.4,f12.5)') kyr,it,tmp(i_w,j_w,it),pre(i_w,j_w,it)
  end do
  !--------------------------------------------------------------------!
 end if
 !---------------------------------------------------------------------!

 !---------------------------------------------------------------------!
 ! Set total soil C using relationship in Fig. 2 of
 ! FW00, given by Eqn. 3.
 ! Annual potential evaporation, from Linacre (1977) (m yr-1).
 ! Could use 30-yr climatology, but this is simpler for now.
 !---------------------------------------------------------------------!
  isc_total = 0.0
  isc (:, :) = 0.0
  isc_grid (:,:) = 0.0
  !--------------------------------------------------------------------!
  do j = j1, j2
   !-------------------------------------------------------------------!
   do i = i1, i2
    if (tmp (i, j, 1) /= fillvalue) then
     eo = 0.0 ! Penman lake evaporation (m yr-1)
     do it = 1, ntimes
      if (mod (it,4) == 0.0) then
       ! Temperature (degree C).
       t = sum (tmp (i, j, it-3:it)) / 4.0 - tf
       ! Vapour pressure (Pa).
       vap = sum (pres (i, j, it-3:it)) * sum (spfh (i, j, it-3:it)) * &
             m_air / m_water
       ! Vapour pressure (mbar).
       e = vap / 100.0
       e = max (eps, e)
       ! Dew-point, based on https://archive.eol.ucar.edu/projects/
       !ceop/dm/documents/refdata_report/eqns.html (degree C).
       t_d = log (e / 6.112) * 243.5 / (17.67 - log (e / 6.112))
       t_d = min (t, t_d)
       ! Sum Penman lake evaporation (m yr-1).
       eo = eo + 0.001 * (700.0 * (t + 0.6) / &
            (100.0 - abs (lat (j))) + 15.0 * (t - t_d)) / (80.0 - t)
       eo = max (zero, eo)
       !if ((local) .and. (i == i_w) .and. (j == j_w)) &
       ! write (*,*)i_w,j_w,eo,t,t_d
      end if
     end do ! it
     ! Annual precipitation (m yr-1).
     pt = max (eps, sum (pre (i, j, 1:ntimes))) / 1000.0
     ! Local soil C on vegetated fraction (kg[C] / m-2).
     if ((eo / pt) < 0.5) then
      isc (i, j) = 36.7 - 53.3 * eo / pt
     else
      isc (i, j) = 10.8 - 1.6 * eo / pt
     end if ! eo / pt
     isc (i, j) = max (isc (i, j), zero)
     isc_grid (i, j) = (1.0 - icwtr (i, j)) * isc (i, j)
     isc_total = isc_total + isc_grid (i, j) * larea (i, j) * 1.0e6
     !if ((local) .and. (i == i_w) .and. (j == j_w)) &
     ! write (*,*)i_w,j_w,isc (i_w, j_w),eo,pt
     !-----------------------------------------------------------------!
     ! Now assign values to all soil C and N pools.
     !-----------------------------------------------------------------!
     do kp = 1, nplots
      !----------------------------------------------------------------!
      ! Noise on initial soil C across plots?
      !----------------------------------------------------------------!
      if (ran_s) then
       !---------------------------------------------------------------!
       call random_number (ran)
       !---------------------------------------------------------------!
       ! iscv is noise in initial soil C.
       !---------------------------------------------------------------!
       noise = 2.0 * iscv * ran + (1.0 - iscv)
       !---------------------------------------------------------------!
      else
       !---------------------------------------------------------------!
       noise = 1.0
       !---------------------------------------------------------------!
      end if
      !----------------------------------------------------------------!
      ! Initialise C and N pools. These values are different from the
      ! PLANT.f90 of HYBRID4, conforming to the F97 paper.
      !----------------------------------------------------------------!
      ! Surface structural litter C and N (kg m-2).
      !----------------------------------------------------------------!
      Cu  (i, j, kp) = isc (i, j) * noise * 0.5  / 14.34
      Nu  (i, j, kp) = Cu (i, j, kp) * vu
      !----------------------------------------------------------------!
      ! Surface metabolic litter C and N (kg m-2).
      !----------------------------------------------------------------!
      Cm  (i, j, kp) = isc (i, j) * noise * 0.0  / 14.34
      vm = 0.07 ! Initial N:C.
      Nm  (i, j, kp) = Cm (i, j, kp) * vm
      !----------------------------------------------------------------!
      ! Soil structural C and N (kg m-2).
      !----------------------------------------------------------------!
      Cv  (i, j, kp) = isc (i, j) * noise * 0.5  / 14.34
      Nv  (i, j, kp) = Cv (i, j, kp) * vv
      !----------------------------------------------------------------!
      ! Soil metabolic litter C and N (kg m-2).
      !----------------------------------------------------------------!
      Cn  (i, j, kp) = isc (i, j) * noise * 0.0  / 14.34
      vn = 0.07 ! Initial N:C.
      Nn  (i, j, kp) = Cn (i, j, kp) * vn
      !----------------------------------------------------------------!
      ! Surface + soil active (microbe) C and N (kg m-2).
      !----------------------------------------------------------------!
      Ca  (i, j, kp) = isc (i, j) * noise * 0.34 / 14.34
      Na  (i, j, kp) = Ca (i, j, kp) * va
      !----------------------------------------------------------------!
      ! Slow C and N (kg m-2).
      !----------------------------------------------------------------!
      Cs  (i, j, kp) = isc (i, j) * noise * 5.0  / 14.34
      Ns  (i, j, kp) = Cs (i, j, kp) * vs
      !----------------------------------------------------------------!
      ! Passive C and N (kg m-2).
      !----------------------------------------------------------------!
      Cpa (i, j, kp) = isc (i, j) * noise * 8.0  / 14.34
      Npa  (i, j, kp) = Cpa (i, j, kp) * vpa
      !----------------------------------------------------------------!
      ! Mineral N (kg[N] m-2).
      !----------------------------------------------------------------!
      snmin (i, j, kp) = 0.004
      !----------------------------------------------------------------!
      ! No. individual trees + 2 (grass types) in plot.
      !----------------------------------------------------------------!
      nind (i, j, kp) = 3
      !----------------------------------------------------------------!
      ! Soil water (m).
      !----------------------------------------------------------------!
      soilw1 (i, j, kp) = 0.1
      soilw2 (i, j, kp) = 0.1
      soilw3 (i, j, kp) = 0.1
      !----------------------------------------------------------------!
      if (local) write (*, *) kp, isc (i, j), Cm (i, j, kp), &
       Cv (i, j, kp), Cn (i, j, kp), Ca (i, j, kp), Cs (i, j, kp), &
       Cpa (i, j, kp)
      if (local) write (*, *) kp, isc (i, j), Nm (i, j, kp), &
       Nv (i, j, kp), Nn (i, j, kp), Na (i, j, kp), Ns (i, j, kp), &
       Npa (i, j, kp)
     end do ! kp
    end if ! fillvalue
   end do ! i
  end do ! j
  if (.NOT. (local)) then
   write (*,*) 'isc_total = ', isc_total / 1.0e12, 'Pg[C]'
   ! Output the global isc field.
   file_name = "isc_grid.nc"
   write (*, *) 'Writing to ', trim (file_name)
   ! Create netCDF dataset and enter define mode.
   call check (nf90_create (trim (file_name), cmode = nf90_clobber, &
               ncid = ncid))
   ! Define the dimensions.
   call check (nf90_def_dim (ncid, "longitude", nlon, lon_dimid))
   call check (nf90_def_dim (ncid, "latitude" , nlat, lat_dimid))
   ! Define coordinate variables.
   call check (nf90_def_var (ncid, "longitude", nf90_float, lon_dimid, &
               lon_varid))
   call check (nf90_def_var (ncid, "latitude" , nf90_float, lat_dimid, &
               lat_varid))
   dimids_two = (/ lon_dimid, lat_dimid /)
   ! Assign units attributes to coordinate data.
   call check (nf90_put_att (ncid, lon_varid, "units", "degrees_east"))
   call check (nf90_put_att (ncid, lat_varid, "units", "degrees_north"))
   ! Define variable.
   call check (nf90_def_var (ncid, "Soil_C", nf90_float, dimids_two, varid))
   call check (nf90_put_att (ncid, varid, "units", "kg[C] m-2"))
   call check (nf90_put_att (ncid, varid, "_FillValue", fillvalue))
   ! End definitions.
   call check (nf90_enddef (ncid))
   ! Write data.
   call check (nf90_put_var (ncid, lon_varid, lon))
   call check (nf90_put_var (ncid, lat_varid, lat))
   call check (nf90_put_var (ncid,     varid, isc_grid))
   ! Close file.
   call check (nf90_close (ncid))
  end if
 end do ! kyr = syr, syr
end if ! .NOT. (rsf_in)
  !--------------------------------------------------------------------!
 !---------------------------------------------------------------------!
 
 !---------------------------------------------------------------------!
 ! Loop through gridboxes and integrate state variables at land points.
 ! Climate files start at 
 !---------------------------------------------------------------------!
do kyr = syr, eyr

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
 nland = 0
 do j = j1, j2
  do i = i1, i2
   ! Land?
   if (tmp (i, j, 1) /= fillvalue) then
    nland = nland + 1
    do it = 1, ntimes
     ! End of day?
     if (mod (it, 4) == zero) then
      do kp = 1, nplots
       wlittc (kp) = zero
       wlittn (kp) = zero
       flittc (kp) = zero
       flittn (kp) = zero
       rlittc (kp) = zero
       rlittn (kp) = zero
       do ki = 1, nind (i, j, kp)
        !****adf
        wlitterc = 0.0
        wlittern = 0.0
        flitterc = 0.0
        flittern = 0.0
        rlitterc = 0.0
        rlittern = 0.0
        !****adf
        wlittc (kp) = wlittc (kp) + wlitterc
        wlittn (kp) = wlittn (kp) + wlittern
        flittc (kp) = flittc (kp) + flitterc
        flittn (kp) = flittn (kp) + flittern
        rlittc (kp) = rlittc (kp) + rlitterc
        rlittn (kp) = rlittn (kp) + rlittern
       end do
      end do
      ! et_soil is cited in F97 as from Comins & McMurtrie (1993).
      ! Equation from HYBRID4 code.
      tsoil = sum (tmp (i, j, it-3:it)) / 4.0 - tf
      if (tsoil > eps) then
       et_soil = 0.0326 + 0.00351 * tsoil ** 1.652 - &
                (0.023953 * tsoil) ** 7.19
      else
       et_soil = 0.0326
      end if
      et_soil = max (zero, et_soil)
      et_soil = min (one , et_soil)
      !write (*,*) it/4,Cpa(i,j,1),et_soil
      do kp = 1, nplots
       wlittc (kp) = wlittc (kp) / area
       wlittn (kp) = wlittn (kp) / area
       flittc (kp) = flittc (kp) / area
       flittn (kp) = flittn (kp) / area
       rlittc (kp) = rlittc (kp) / area
       rlittn (kp) = rlittn (kp) / area
       ! Soil water holding capacity (Eqn. (1) of FW00; m).
       swct = 0.213 + 0.00227 * (Cm (i, j, kp) + Cu (i, j, kp) + &
        Cn (i, j, kp) + Cv (i, j, kp) + Ca (i, j, kp) + &
        Cs (i, j, kp) + Cpa (i, j, kp))
       ! Convert soil water to water-filled pore space.
       ! Assumes micro-pore space = swc and macro-pore space = 42% saturation
       ! content (from TEM, for loam; Raich et al., 1991).
       if (swct > eps) then
        wfps = 100.0 * (soilw1 (i, j, kp) + soilw2 (i, j, kp) + &
                        soilw3 (i, j, kp)) / (1.7241 * swct)
       else
        wfps = 0.00001
       end if
       ! Soil-water decomposition modifer (fraction). Equations are 
       ! fitted to Fig. 8 of Williams et al. (1992).
       if (wfps < 60.0) then
        em_soil = exp ((wfps - 60.0) ** 2 / (-800.0))
       else
        em_soil = 0.000371 * wfps ** 2 - 0.0748 * wfps + 4.13
       endif
       em_soil = max (zero, em_soil)
       em_soil = min (one , em_soil)
       ! Overall decomposition modifier.
       ev = et_soil * em_soil
      write (*,*) it/4,Cpa(i,j,1),et_soil,em_soil,ev
      end do
     end if ! 
    end do ! it = 1, ntimes
!    stop
   end if ! tmp (i,j,k) /= fillvalue
  end do ! j
 end do ! i
 !---------------------------------------------------------------------!
end do ! kyr = syr, eyr
!----------------------------------------------------------------------!

! Output state variables to restart file if required.
if (rsf_out) then
   file_name = "/home/adf10/rds/rds-mb425-geogscratch/adf10/RSF/rsf_out.nc"
   write (*, *) 'Writing to ', trim (file_name)
   ! Create netCDF dataset and enter define mode.
   call check (nf90_create (trim (file_name), cmode = nf90_clobber, &
               ncid = ncid))
   ! Define the dimensions.
   call check (nf90_def_dim (ncid, "longitude", nlon  , lon_dimid))
   call check (nf90_def_dim (ncid, "latitude" , nlat  , lat_dimid))
   call check (nf90_def_dim (ncid, "plot"     , nplots, plot_dimid))
   ! Define coordinate variables.
   call check (nf90_def_var (ncid, "longitude", nf90_float, lon_dimid, &
               lon_varid))
   call check (nf90_def_var (ncid, "latitude" , nf90_float, lat_dimid, &
               lat_varid))
   call check (nf90_def_var (ncid, "plot"     , nf90_int  , plot_dimid, &
               plot_varid))
   dimids_three = (/ lon_dimid, lat_dimid, plot_dimid /)
   ! Assign units attributes to coordinate data.
   call check (nf90_put_att (ncid, lon_varid, "units", "degrees_east"))
   call check (nf90_put_att (ncid, lat_varid, "units", "degrees_north"))
   ! Define variables.
   call check (nf90_def_var (ncid, "Surface_metabolic_organic_C", &
    nf90_float, dimids_three, varid_Cm))
   call check (nf90_def_var (ncid, "Surface_structural_organic_C", &
    nf90_float, dimids_three, varid_Cu))
   call check (nf90_def_var (ncid, "Soil_metabolic_organic_C", &
    nf90_float, dimids_three, varid_Cn))
   call check (nf90_def_var (ncid, "Soil_structural_organic_C", &
    nf90_float, dimids_three, varid_Cv))
   call check (nf90_def_var (ncid, "Active_organic_C", &
    nf90_float, dimids_three, varid_Ca))
   call check (nf90_def_var (ncid, "Slow_organic_C", &
    nf90_float, dimids_three, varid_Cs))
   call check (nf90_def_var (ncid, "Passive_soil_organic_C", &
    nf90_float, dimids_three, varid_Cpa))
   call check (nf90_def_var (ncid, "Surface_metabolic_organic_N", &
    nf90_float, dimids_three, varid_Nm))
   call check (nf90_def_var (ncid, "Surface_structural_organic_N", &
    nf90_float, dimids_three, varid_Nu))
   call check (nf90_def_var (ncid, "Soil_metabolic_organic_N", &
    nf90_float, dimids_three, varid_Nn))
   call check (nf90_def_var (ncid, "Soil_structural_organic_N", &
    nf90_float, dimids_three, varid_Nv))
   call check (nf90_def_var (ncid, "Active_organic_N", &
    nf90_float, dimids_three, varid_Na))
   call check (nf90_def_var (ncid, "Slow_organic_N", &
    nf90_float, dimids_three, varid_Ns))
   call check (nf90_def_var (ncid, "Passive_soil_organic_N", &
    nf90_float, dimids_three, varid_Npa))
   call check (nf90_def_var (ncid, "Mineral_N", &
    nf90_float, dimids_three, varid_snmin))
   call check (nf90_def_var (ncid, "Number_individuals", &
    nf90_int, dimids_three, varid_nind))
   call check (nf90_def_var (ncid, "Soil_water_layer_1", &
    nf90_float, dimids_three, varid_soilw1))
   call check (nf90_def_var (ncid, "Soil_water_layer_2", &
    nf90_float, dimids_three, varid_soilw2))
   call check (nf90_def_var (ncid, "Soil_water_layer_3", &
    nf90_float, dimids_three, varid_soilw3))
   call check (nf90_put_att (ncid, varid_Cm , "units", "kg[C] m-2"))
   call check (nf90_put_att (ncid, varid_Cu , "units", "kg[C] m-2"))
   call check (nf90_put_att (ncid, varid_Cn , "units", "kg[C] m-2"))
   call check (nf90_put_att (ncid, varid_Cv , "units", "kg[C] m-2"))
   call check (nf90_put_att (ncid, varid_Ca , "units", "kg[C] m-2"))
   call check (nf90_put_att (ncid, varid_Cs , "units", "kg[C] m-2"))
   call check (nf90_put_att (ncid, varid_Cpa, "units", "kg[C] m-2"))
   call check (nf90_put_att (ncid, varid_Nm , "units", "kg[N] m-2"))
   call check (nf90_put_att (ncid, varid_Nu , "units", "kg[N] m-2"))
   call check (nf90_put_att (ncid, varid_Nn , "units", "kg[N] m-2"))
   call check (nf90_put_att (ncid, varid_Nv , "units", "kg[N] m-2"))
   call check (nf90_put_att (ncid, varid_Na , "units", "kg[N] m-2"))
   call check (nf90_put_att (ncid, varid_Ns , "units", "kg[N] m-2"))
   call check (nf90_put_att (ncid, varid_Npa, "units", "kg[N] m-2"))
   call check (nf90_put_att (ncid, varid_snmin, "units", "kg[N] m-2"))
   call check (nf90_put_att (ncid, varid_nind, "units", "n"))
   call check (nf90_put_att (ncid, varid_soilw1, "units", "m"))
   call check (nf90_put_att (ncid, varid_soilw2, "units", "m"))
   call check (nf90_put_att (ncid, varid_soilw3, "units", "m"))
   ! End definitions.
   call check (nf90_enddef (ncid))
   ! Write data.
   call check (nf90_put_var (ncid, lon_varid, lon))
   call check (nf90_put_var (ncid, lat_varid, lat))
   call check (nf90_put_var (ncid, plot_varid, plot))
   call check (nf90_put_var (ncid, varid_Cm , Cm))
   call check (nf90_put_var (ncid, varid_Cu , Cu))
   call check (nf90_put_var (ncid, varid_Cn , Cn))
   call check (nf90_put_var (ncid, varid_Cv , Cv))
   call check (nf90_put_var (ncid, varid_Ca , Ca))
   call check (nf90_put_var (ncid, varid_Cs , Cs))
   call check (nf90_put_var (ncid, varid_Cpa, Cpa))
   call check (nf90_put_var (ncid, varid_Nm , Nm))
   call check (nf90_put_var (ncid, varid_Nu , Nu))
   call check (nf90_put_var (ncid, varid_Nn , Nn))
   call check (nf90_put_var (ncid, varid_Nv , Nv))
   call check (nf90_put_var (ncid, varid_Na , Na))
   call check (nf90_put_var (ncid, varid_Ns , Ns))
   call check (nf90_put_var (ncid, varid_Npa, Npa))
   call check (nf90_put_var (ncid, varid_snmin, snmin))
   call check (nf90_put_var (ncid, varid_nind, nind))
   call check (nf90_put_var (ncid, varid_soilw1, soilw1))
   call check (nf90_put_var (ncid, varid_soilw2, soilw2))
   call check (nf90_put_var (ncid, varid_soilw3, soilw3))
   ! Close file.
   call check (nf90_close (ncid))
end if

if ((local) .and. (wrclim)) close (20) ! Local climate output.

!----------------------------------------------------------------------!
write (*,*) 'nland = ', nland
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

