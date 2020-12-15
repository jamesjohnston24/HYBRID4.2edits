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
logical :: local = .FALSE. ! Run only local site?
logical :: wrclim = .FALSE. ! Write local climate?
!----------------------------------------------------------------------!
! Uggla site.
!----------------------------------------------------------------------!
!real, parameter :: lon_w = 19.0 + 46.0 / 60.0
!real, parameter :: lat_w = 64.0 + 21.0 / 60.0
!----------------------------------------------------------------------!
! Grudd site.
!----------------------------------------------------------------------!
real, parameter :: lon_w = 19.75
real, parameter :: lat_w = 68.25
!----------------------------------------------------------------------!
! Indices of local site, if used.
!----------------------------------------------------------------------!
integer :: i_w
integer :: j_w
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
integer, parameter :: nlon = 720
integer, parameter :: nlat = 360
integer, parameter :: ntimes = 1460
real, parameter :: fillvalue = 1.0e20
real, parameter :: tf = 273.15
real, parameter :: eps = 1.0e-8
real, parameter :: zero = 0.0
! Mol. weight of dry air (kg [air] mol-1)
real, parameter :: m_air = 28.9647
! Mol. weight of water (kg [water] mol-1)
real, parameter :: m_water = 18.01528
! Earth's circumference at equator (m2).
real, parameter :: circ = 40075.0e3
real, parameter :: pi = 3.14159265359
character (len = 250) :: file_name ! Generic filename
character (len = 100) :: var_name
character (len =     4) :: char_year
integer :: kyr
integer :: it
integer :: syr
integer :: eyr
integer :: varid
integer :: ncid
integer :: lon_dimid
integer :: lat_dimid
integer :: lon_varid
integer :: lat_varid
integer, dimension (2) :: dimids
integer :: i,ii
integer :: j,jj
integer :: k
integer :: nland
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
real, allocatable, dimension (:,:) :: larea_qd ! Grid-box area (km2).
real, allocatable, dimension (:,:) :: larea ! Grid-box area (km2).
real :: t ! Temperature (degree C)
real :: t_d ! Dew-point (degree C)
real :: eo ! Penman evaporation from lake (mm day-1)
real :: pt ! Annual precipiation (m yr-1)
real :: vap ! Vapour pressure (Pa)
real :: e ! Vapour pressure (mbar)
real :: isc_total ! Total global soil C (Pg[C])
real :: barea ! Area of grid-box at equation (m2)
real :: rlat ! Latitude (radians)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
allocate (lon  (nlon))
allocate (lat  (nlat))
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
!----------------------------------------------------------------------!

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
end if
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
syr = 1901
eyr = 2019
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
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
 ! If first year, set total soil C using relationship in Fig. 2 of
 ! FW00, given by Eqn. 3.
 ! Annual potential evaporation, from Linacre (1977)? (m yr-1).
 !---------------------------------------------------------------------!
 if (kyr == syr) then
  isc_total = 0.0
  isc (:, :) = 0.0
  isc_grid (:,:) = 0.0
  !--------------------------------------------------------------------!
  ! Area of grid-box at equator (m2).
  !--------------------------------------------------------------------!
  barea = (circ / float (nlon)) ** 2
  !--------------------------------------------------------------------!
  do j = 1, nlat
   !-------------------------------------------------------------------!
   ! Latitude (radians).
   !-------------------------------------------------------------------!
   rlat = lat (j) * pi / 180.0
   !-------------------------------------------------------------------!
   do i = 1, nlon
    if (tmp (i, j, 1) /= fillvalue) then
     eo = 0.0 ! Penman lake evaporation (m yr-1)
     do it = 1, ntimes
      ! Temperature (degree C).
      t = tmp (i, j, it) - tf
      ! Vapour pressure (Pa).
      vap = pres (i, j, it) * spfh (i, j, it) * m_air / m_water
      ! Vapour pressure (mbar).
      e = vap / 100.0
      e = max (eps, e)
      ! Dew-point, based on https://archive.eol.ucar.edu/projects/
      !ceop/dm/documents/refdata_report/eqns.html (degree C).
      t_d =  log (e / 6.112) * 243.5 / (17.67 - log (e / 6.112))
      t_d = min (t, t_d)
      ! Sum Penman lake evaporation (m yr-1).
      eo = eo + (0.001 / 4.0) * (700.0 * (t + 0.6) / &
           (100.0 - abs (lat (j))) + 15.0 * (t - t_d)) / (80.0 - t)
      eo = max (zero, eo)
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
     !write (*,*)i,j,isc (i, j),eo,pt
     ! Area of grid-box (m2).
     !larea = cos (rlat) * barea
     isc_grid (i, j) = (1.0 - icwtr (i, j)) * isc (i, j)
     isc_total = isc_total + isc_grid (i, j) * larea (i, j) * 1.0e6
    end if ! fillvalue
   end do ! i
  end do ! j
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
  dimids = (/ lon_dimid, lat_dimid /)
  ! Assign units attributes to coordinate data.
  call check (nf90_put_att (ncid, lon_varid, "units", "degrees_east"))
  call check (nf90_put_att (ncid, lat_varid, "units", "degrees_north"))
  ! Define variable.
  call check (nf90_def_var (ncid, "Soil C", nf90_float, dimids, varid))
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
 end if ! kyr == syr
 !---------------------------------------------------------------------!
 
 !---------------------------------------------------------------------!
 ! Loop through gridboxes and integrate state variables at land points.
 ! Climate files start at 
 !---------------------------------------------------------------------!
 nland = 0
 do i = 1, nlon
  do j = 1, nlat
   k = 1
   !write (*,*) i,j,tmp(i,j,k)
   if (tmp (i,j,k) /= fillvalue) nland = nland + 1
   do it = 1, ntimes
   end do ! it
  end do ! j
 end do ! i
 !---------------------------------------------------------------------!
stop
end do ! kyr = syr, eyr
!----------------------------------------------------------------------!
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
