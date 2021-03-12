subroutine init_yr

use netcdf
use mpi
use shared

implicit none

!----------------------------------------------------------------------!
! Some variables used locally for netCDF operations.
!----------------------------------------------------------------------!
character (len = 250) :: file_name ! Generic filename
character (len = 100) :: var_name  ! Generic variable name
character (len =   4) :: char_year ! CE year as text
integer :: ncid
integer :: varid
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
GPP_global = zero
NPP_global = zero
Cv_global  = zero
Cs_global  = zero
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Flag for cold-deciduous diagnostic.
!----------------------------------------------------------------------!
!cd_flag = 0
!----------------------------------------------------------------------!
 
!----------------------------------------------------------------------!
! Flag for dry-deciduous diagnostic.
!----------------------------------------------------------------------!
dd_flag (:) = 0
!----------------------------------------------------------------------!
 
!----------------------------------------------------------------------!
write (char_year, '(i4)') kyr
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Read global temperature fields for year kyr into tmp (K).
!----------------------------------------------------------------------!
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
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Read precipitation fields for year kyr intp pre (mm 6hr-1).
!----------------------------------------------------------------------!
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
!----------------------------------------------------------------------!
 
!----------------------------------------------------------------------!
! Read global downward solar radiation flux fields for year kyr into
! dswrf (J m-2).
!----------------------------------------------------------------------!
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
!----------------------------------------------------------------------!
 
!----------------------------------------------------------------------!
! Read specific humidity fields for year kyr into spfh (kg kg-1).
!----------------------------------------------------------------------!
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
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Read pressure fields for year kyr into pres (Pa).
!----------------------------------------------------------------------!
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
!----------------------------------------------------------------------!
 
!----------------------------------------------------------------------!
! CO2 stomatal response (assuming closing above 80 Pa). Should give
! curve based on Fig. 1 of M&G. kyr must be 1700-1900.
!----------------------------------------------------------------------!
fc = one - 0.008333 * cao (kyr-1699)
fc = max (0.33, fc)
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

end subroutine init_yr
