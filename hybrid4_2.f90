!==========================================================!
program hybrid4_2
!---------------------------------------------------------------------------------------------------------------------------------------!

!---------------------------------------------------------------------------------------------------------------------------------------!
! Learn from /home/adf10/MODELS/HYBRID10 on CSD3 for how to use TRENDY forcings.
! FW00 is Friend AD, White A. 2000. Eval...
!---------------------------------------------------------------------------------------------------------------------------------------!

!---------------------------------------------------------------------------------------------------------------------------------------!
use netcdf
use mpi
!---------------------------------------------------------------------------------------------------------------------------------------!

!---------------------------------------------------------------------------------------------------------------------------------------!
implicit none
!---------------------------------------------------------------------------------------------------------------------------------------!

!---------------------------------------------------------------------------------------------------------------------------------------!
integer, parameter :: nlon = 720
integer, parameter :: nlat = 360
integer, parameter :: ntimes = 1460
real, parameter :: fillvalue = 1.0e20
real, parameter :: tf = 273.15
real, parameter :: eps = 1.0e-8
real, parameter :: zero = 0.0
! Uggla site.
real, parameter :: lon_w = 19.0 + 46.0 / 60.0
real, parameter :: lat_w = 64.0 + 21.0 / 60.0
! Grudd site.
!real, parameter :: lon_w = 19.75
!real, parameter :: lat_w = 68.25
character (len = 250) :: file_name ! Generic filename
character (len = 100) :: var_name
character (len =     4) :: char_year
integer :: kyr
integer :: it
integer :: syr
integer :: eyr
integer :: varid
integer :: ncid
integer :: i
integer :: j
integer :: k
integer :: nland
integer :: i_w
integer :: j_w
real, allocatable, dimension (:,:,:) :: tmp ! K
real, allocatable, dimension (:,:,:) :: pre ! mm 6hr-1
real, allocatable, dimension (:,:,:) :: spfh ! kg kg-1
real, allocatable, dimension (:,:,:) :: pres ! Pa
real, allocatable, dimension (: ) :: lat ! degrees north
real, allocatable, dimension (:) :: lon ! degrees east
real :: isc ! Initial total soil C (kg[C] m-2)
real :: t ! Temperature (degree C)
real :: t_d ! Dew-point (degree C)
real :: eo ! Penman evaporation from lake (mm day-1)
real :: pt ! Annual precipiation (m yr-1)
real :: vap ! Vapour pressure (Pa)
real :: e ! Vapour pressure (mbar)
real :: m_air = 28.9647 ! Molecular weight of dry air (kg [air] mol-1)
real :: m_water = 18.01528 ! Molecular weight of water (kg [water] mol-1)
!---------------------------------------------------------------------------------------------------------------------------------------!

!---------------------------------------------------------------------------------------------------------------------------------------!
allocate (tmp (nlon,nlat,ntimes))
allocate (pre (nlon,nlat,ntimes))
allocate (spfh (nlon,nlat,ntimes))
allocate (pres (nlon,nlat,ntimes))
allocate (lon(nlon))
allocate (lat (nlat))
!---------------------------------------------------------------------------------------------------------------------------------------!

!---------------------------------------------------------------------------------------------------------------------------------------!
! Gridbox indices for containing wanted site in climate data.
!---------------------------------------------------------------------------------------------------------------------------------------!
i_w = nint ((719.0 / 2.0) * (1.0 + lon_w / 179.75)) + 1
j_w = nint ((359.0 / 2.0) * (1.0 + lat_w /   89.75)) + 1
write(*,*) 'i_w j_w ',i_w,j_w
if (lat_w == (64.0 + 21.0 / 60.0)) open (20,file='uggla.clm',status='unknown')
if (lat_w == 68.25) open (20,file='grudd.clm',status='unknown')
!---------------------------------------------------------------------------------------------------------------------------------------!

!---------------------------------------------------------------------------------------------------------------------------------------!
syr = 1901
eyr = 2019
!---------------------------------------------------------------------------------------------------------------------------------------!

!---------------------------------------------------------------------------------------------------------------------------------------!
do kyr = syr, eyr

 !--------------------------------------------------------------------------------------------------------------------------------------!
 write (char_year, '(i4)') kyr
 !--------------------------------------------------------------------------------------------------------------------------------------!

 !--------------------------------------------------------------------------------------------------------------------------------------!
 ! Read global temperature fields for year kyr into tmp (K).
 !--------------------------------------------------------------------------------------------------------------------------------------!
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
 !--------------------------------------------------------------------------------------------------------------------------------------!
 ! Read precipitation fields for year kyr intp pre (mm 6hr-1.
 !--------------------------------------------------------------------------------------------------------------------------------------!
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
 !--------------------------------------------------------------------------------------------------------------------------------------!
 ! Read specific humidity fields for year kyr intp spfh (kg kg-1).
 !--------------------------------------------------------------------------------------------------------------------------------------!
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
 !--------------------------------------------------------------------------------------------------------------------------------------!
 ! Read pressure fields for year kyr intp pres (Pa).
 !--------------------------------------------------------------------------------------------------------------------------------------!
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
 !--------------------------------------------------------------------------------------------------------------------------------------!
 
 !--------------------------------------------------------------------------------------------------------------------------------------!
 ! Climate output if needed.
 do it = 1, ntimes
  write (20,'(2i5,f12.4,f12.5)') kyr,it,tmp(i_w,j_w,it),pre(i_w,j_w,it)
 end do
 !--------------------------------------------------------------------------------------------------------------------------------------!

 ! If first year, set total soil C using relationship in Fig. 2 of FW00, given by Eqn. 3.
 ! Annual potential evaporation, from Linacre (1977)? (m yr-1).
 if (kyr == syr) then
  do j = 1, nlat
   do i = 1, nlon
    if (tmp (i, j, 1) /= fillvalue) then
     eo = 0.0 ! Penman lake evaporation (m yr-1)
     do it = 1, ntimes
      ! Temperature (degree C)>
      t = tmp (i, j, it) - tf
      ! Vapour pressure (Pa).
      vap = pres (i, j, it) * spfh (i, j, it) * m_air / m_water
      ! Vapour pressure (mbar).
      e = vap / 100.0
      e = max (eps, e)
      ! Dew-point, based on https://archive.eol.ucar.edu/projects/ceop/dm/documents/
      ! refdata_report/eqns.html (degree C).
      t_d =  log (e / 6.112) * 243.5 / (17.67 - log (e / 6.112))
      t_d = min (t, t_d)
      ! Sum Penman lake evaporation (m yr-1).
      eo = eo + (0.001 / 4.0) * (700.0 * (t + 0.6) / (100.0 - abs (lat (j))) + 15.0 * (t - t_d)) / &
              (80.0 - t)
      eo = max (zero, eo)
     end do ! it
     ! Annual precipitation (m yr-1).
     pt = max (eps, sum (pre (i, j, 1:ntimes))) / 1000.0
     if ((eo / pt) < 0.5) then
      isc = 36.7 - 53.3 * eo / pt
     else
      isc = 10.8 - 1.6 * eo / pt
     end if ! eo / pt
     !write(*,*)lat(j),pt,eo,isc
    end if ! fillvalue
   end do ! i
  end do ! j
 end if ! kyr == syr
 
 !--------------------------------------------------------------------------------------------------------------------------------------!
 ! Loop through gridboxes and integrate state variables at land points.
 ! Climate files start at 
 !--------------------------------------------------------------------------------------------------------------------------------------!
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
 !--------------------------------------------------------------------------------------------------------------------------------------!

end do ! kyr = syr, eyr
!---------------------------------------------------------------------------------------------------------------------------------------!
 close (20) ! Climate output if needed.

!---------------------------------------------------------------------------------------------------------------------------------------!
write (*,*) 'nland = ', nland
!---------------------------------------------------------------------------------------------------------------------------------------!

!---------------------------------------------------------------------------------------------------------------------------------------!
contains
 subroutine check (status)

 integer, intent (in) :: status
 if (status /= nf90_noerr) then
  print *, trim (nf90_strerror(status))
  stop  "Stopped"
 end if
 end subroutine check
!---------------------------------------------------------------------------------------------------------------------------------------!

end program hybrid4_2
!==========================================================!
