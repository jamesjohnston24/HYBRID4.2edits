subroutine initl

use netcdf
use mpi
use shared
implicit none

integer, parameter :: ntimes      =  1460
character (len = 250) :: file_name ! Generic filename
integer, parameter :: fillvalue_int = -99999
character (len = 100) :: var_name
character (len =   4) :: char_year
integer :: ncid
integer :: varid
integer :: plot_dimid, ind_p_dimid, ind_t_dimid
integer :: plot_varid, ind_p_varid, ind_t_varid
integer :: lon_dimid, lat_dimid
integer :: lon_varid, lat_varid
integer, dimension (2) :: dimids_two

real, parameter :: swpmax = -0.033 ! MPa
real, parameter :: bsoil = 5.0
integer :: kp
integer :: ki
integer :: ksp
integer :: k
integer :: p
integer :: il
integer :: m
integer :: nplots_total
integer :: iregl,ireglc
integer :: nlayers
real :: warea,harea,tarea,saparea
real :: wheight,hheight,theight
real :: wwood,hwood,swood
real :: diamw
real :: ht
real :: nfp
real :: sow1,sc1,sw1,sw2,tr_rat
real :: fd,kf,kSWf,frac,fracs
real :: t_d ! Dew-point (degree C)
real :: t ! Temperature (degree C)
real :: eo ! Penman evaporation from lake (mm day-1)
real :: e ! Vapour pressure (mbar)
real :: pt ! Annual precipitation (m yr-1)
real :: isc_total ! Total global soil C (Pg[C])
real :: ran ! Random number (0-1)
real :: noise ! Noise on initial soil C (ratio)
real :: caod ! ppm
real, dimension (mh+1) :: fPARt
real, dimension (mh+1) :: fSWt
real, dimension (mh+1) :: fPARi
real, dimension (mh+1) :: fSWi
real, dimension (mh+1) :: skz
real, dimension (mh+1) :: skSWz
real, dimension (mh+1) :: eskz
real, dimension (mh+1) :: eskSWz
real, allocatable, dimension (:,:) :: isc ! Initial soil C (kg[C] m-2)
real, allocatable, dimension (:,:) :: isc_grid ! Init soil C (kg[C] m-2)

!----------------------------------------------------------------------!
! Variables from 'driver.in'.
!----------------------------------------------------------------------!
real :: lon_w
real :: lat_w
!****adflogical :: rsf_in  ! Start from restart file?
!****adflogical :: ran_s   ! Noise on initial soil C across plots?
!****adfreal?integer :: iscv    ! Range of initial soil C (+/-fraction)
!****adfinteger :: syr     ! Start year of simulation (>= 1901)
!****adfinteger :: eyr     ! End year of simulation (<= 2019)
integer :: imax    ! Max. tree density (tree/ha)
!----------------------------------------------------------------------!
!****adfinteger :: limp     ! Max. trees/plot
!----------------------------------------------------------------------!
! Indices of local site, if used.
!----------------------------------------------------------------------!
!****adfinteger :: i_w
!****adfinteger :: j_w
!----------------------------------------------------------------------!
! Initial total soil C on non-ice/water fraction (kg[C] m-2)
allocate (isc  (nlon,nlat))
! Initial total soil C grid-box mean (kg[C] m-2)
allocate (isc_grid  (nlon,nlat))

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
land_index (:,:) = 0
p_plot (:,:) = 0
k_ind  (:,:,:) = 0
snow   (:,:,:) = fillvalue
soilw1 (:,:,:) = fillvalue
soilw2 (:,:,:) = fillvalue
soilw3 (:,:,:) = fillvalue

uwind = 5.0

!----------------------------------------------------------------------!
! Read simulation control values from 'driver.in'.
!----------------------------------------------------------------------!
open (10,file='driver.in',status='old')
read (10,*) local
read (10,*) lon_w
read (10,*) lat_w
read (10,*) rsf_in
read (10,*) rsf_out
read (10,*) ran_s
read (10,*) iscv
read (10,*) syr
read (10,*) eyr
read (10,*) nplots
read (10,*) area
read (10,*) imax
read (10,*) idbh
read (10,*) idbhv
read (10,*) ndepi
! Convert from kg[N] ha-1 yr-1 to kg[N] m-2 d-1.
ndepi = ndepi / (10000.0 * float (nd))
close (10)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
limp = imax * NINT (area) / 10000 ! Max. trees/plotdo kp = 1, nplots
nind_max = limp + 2                   ! Max. individuals/plot (>1).
allocate (p_plot     (nland_max,nplots))
allocate (k_ind      (nland_max,nplots,nind_max))
allocate (dd_flag    (nplots))
allocate (phenf      (nplots))
allocate (swp1       (nplots))
allocate (swp2       (nplots))
allocate (fturn_plot (nplots,nspp))
allocate (lon  (nlon))
allocate (lat  (nlat))
allocate (plot (nplots))
allocate (icwtr (nlon,nlat))
allocate (larea (nlon,nlat))
allocate (tmp   (nlon,nlat,ntimes))
allocate (dswrf (nlon,nlat,ntimes))
allocate (pre   (nlon,nlat,ntimes))
allocate (spfh (nlon,nlat,ntimes))
allocate (pres (nlon,nlat,ntimes))
allocate (Cm   (nlon,nlat,nplots))
allocate (Cu   (nlon,nlat,nplots))
allocate (Cn   (nlon,nlat,nplots))
allocate (Cv   (nlon,nlat,nplots))
allocate (Ca   (nlon,nlat,nplots))
allocate (Cs   (nlon,nlat,nplots))
allocate (Cpa  (nlon,nlat,nplots))
allocate (soilC (nlon,nlat,nplots))
allocate (Nu   (nlon,nlat,nplots))
allocate (Nm   (nlon,nlat,nplots))
allocate (Nv   (nlon,nlat,nplots))
allocate (Nn   (nlon,nlat,nplots))
allocate (Na   (nlon,nlat,nplots))
allocate (Ns   (nlon,nlat,nplots))
allocate (Npa  (nlon,nlat,nplots))
allocate (snmin (nlon,nlat,nplots))
allocate (soilN (nlon,nlat,nplots))
allocate (nind (nlon,nlat,nplots))
allocate (laip (nlon,nlat,nplots))
allocate (swct   (nlon,nlat,nplots))
allocate (snow   (nlon,nlat,nplots))
allocate (soilw1 (nlon,nlat,nplots))
allocate (soilw2 (nlon,nlat,nplots))
allocate (soilw3 (nlon,nlat,nplots))
allocate (swc1 (nlon,nlat,nplots))
allocate (swc2 (nlon,nlat,nplots))
allocate (swc3 (nlon,nlat,nplots))
allocate (outflow(nplots))
allocate (wfps   (nplots))
allocate (tfolK  (nplots))
allocate (rcstom (nind_max))
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Open diagnostic output files.
!----------------------------------------------------------------------!
open (20,file='annual_global.dat',status='unknown')
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Read annual global CO2 mixing ratios (ppm -> Pa).
! Do everything here that is in co2.f of original code.
!----------------------------------------------------------------------!
open (10,file='/home/adf10/rds/rds-mb425-geogscratch/adf10/FORCINGS/&
&CO2field/global_co2_ann_1700_2019.txt',status='old')
do kyr = 1700, 2019
 read (10,*) i, caod
 cao   (kyr-1699) = caod / 10.0 ! Pa
end do
close (10)
! N deposition at 0 m precipitation (kg[N] m-2 day-).
nf = ndepi
!----------------------------------------------------------------------!

do kp = 1, nplots
 do ksp = 1, nspp
  fturn_plot (kp,ksp) = fturn (ksp)
 end do ! ksp
end do ! kp
! Grass coarse root to fine root biomass ratio.
crratio (1) = bark (1)
crratio (2) = bark (2)
do ksp = 1, nspp
 fturn_save (ksp) = fturn (ksp)
 ! Maximum stomatal conductance factor (m s-1).
 ngr   (ksp) = ngr (ksp) / 41.0
 pruba (ksp) = (one - ao (ksp)) / (one + 12.5 / nrc (ksp))
 slope (ksp) = 71.4 / (one + 12.5 / nrc (ksp))
 ! Convert to integers for tests.
 weff  (ksp) = nint (weffi  (ksp))
 ptype (ksp) = nint (ptypei (ksp))
 stf (ksp) = one / (one - stf (ksp))
 rah (ksp) = one / (1.5 * 6.62e-03 * (uwind / d_leaf (ksp)) ** 0.5)
 rhr (ksp) = one / (one / rah (ksp) + one / 213.21)
 rac (ksp) = rah (ksp) / 0.76
 raw (ksp) = rah (ksp) / 1.08
 f1  (ksp) = stf (ksp) * formf (ksp) * woodd (ksp) * &
             pi * ah (ksp) / 4.0
 f2  (ksp) = one / (bh (ksp) + 2.0)
 f3  (ksp) = stf (ksp) * formf (ksp) * woodd (ksp)
end do ! ksp

! Initial value of N:C ratio of surface metabolic litter.
vm = 0.07
! Initial value of N:C ratio of soil metabolic litter.
vn = 0.07
!----------------------------------------------------------------------!
! Boundary layer resistance factor for gammas.
!----------------------------------------------------------------------!
raf = (rah (1) + raw (1)) / rhr (1)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
if (local) then
 !---------------------------------------------------------------------!
 ! Gridbox indices for local site.
 !---------------------------------------------------------------------!
 i_w = nint ((719.0 / 2.0) * (1.0 + lon_w / 179.75)) + 1
 j_w = nint ((359.0 / 2.0) * (1.0 + lat_w /   89.75)) + 1
 write (*,*) 'lon_w lat_w',lon_w,lat_w
 write (*,*) 'i_w   j_w  ',i_w  ,j_w
 if (lat_w == (64.0 + 21.0 / 60.0)) open (20,file='uggla.clm',&
  status='unknown')
 if (lat_w == 68.25) open (20,file='grudd.clm',status='unknown')
 !---------------------------------------------------------------------!
 open (21,file='soil_carbon.txt',status='unknown')
 open (22,file= 'soil_water.txt',status='unknown')
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
end if ! local
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
if (rsf_in) then
 file_name = '/home/adf10/rds/rds-mb425-geogscratch/adf10/RSF/rsf_in.nc'
 write (*,*) 'Reading from ',trim(file_name)
 call check (nf90_open (trim (file_name), nf90_nowrite, ncid))
 ! Get number of individuals.
 call check (nf90_inq_dimid(ncid,"ind_total",ind_t_dimid))
 call check (nf90_inquire_dimension(ncid,ind_t_dimid,len=nind_total))
 write(*,*) 'nind_total = ',nind_total
 allocate (ind     (nind_total)) !****adf what is this?
 allocate (kspp    (nind_total))
 allocate (dbh     (nind_total))
 allocate (hdbh    (nind_total))
 allocate (cfoliage(nind_total))
 allocate (cfiner  (nind_total))
 allocate (cstore  (nind_total))
 allocate (nfoliage(nind_total))
 allocate (nbswood (nind_total))
 allocate (nheart  (nind_total))
 allocate (navail  (nind_total))
 allocate (nfiner  (nind_total))
 allocate (lasa    (nind_total))
 allocate (nitf    (nind_total))
 allocate (hbc     (nind_total))
 varid = 1
 call check (nf90_get_var (ncid, varid, lon))
 varid = varid + 1
 call check (nf90_get_var (ncid, varid, lat))
 !varid = varid + 1
 !call check (nf90_get_var (ncid, varid, plot))
 varid = varid + 5
 call check (nf90_get_var (ncid, varid, land_index))
 varid = varid + 1
 call check (nf90_get_var (ncid, varid, cd))
 varid = varid + 1
 call check (nf90_get_var (ncid, varid, dd))
 varid = varid + 1
 call check (nf90_get_var (ncid, varid, bgs))
 varid = varid + 1
 call check (nf90_get_var (ncid, varid, egs))
 !***p_plot and wlittc as well, etc. litter pools, wlittn, alive
 ! flittc, flittn, rlittc, rlittn
 varid = varid + 1
 call check (nf90_get_var (ncid, varid, k_ind))
 varid = varid + 1
 call check (nf90_get_var (ncid, varid, Cm))
 varid = varid + 1
 call check (nf90_get_var (ncid, varid, Cu))
 varid = varid + 1
 call check (nf90_get_var (ncid, varid, Cn))
 varid = varid + 1
 call check (nf90_get_var (ncid, varid, Cv))
 varid = varid + 1
 call check (nf90_get_var (ncid, varid, Ca))
 varid = varid + 1
 call check (nf90_get_var (ncid, varid, Cs))
 varid = varid + 1
 call check (nf90_get_var (ncid, varid, Cpa))
 varid = varid + 1
 call check (nf90_get_var (ncid, varid, Nm))
 varid = varid + 1
 call check (nf90_get_var (ncid, varid, Nu))
 varid = varid + 1
 call check (nf90_get_var (ncid, varid, Nn))
 varid = varid + 1
 call check (nf90_get_var (ncid, varid, Nv))
 varid = varid + 1
 call check (nf90_get_var (ncid, varid, Na))
 varid = varid + 1
 call check (nf90_get_var (ncid, varid, Ns))
 varid = varid + 1
 call check (nf90_get_var (ncid, varid, Npa))
 varid = varid + 1
 call check (nf90_get_var (ncid, varid, snmin))
 varid = varid + 1
 call check (nf90_get_var (ncid, varid, nind))
 varid = varid + 1
 call check (nf90_get_var (ncid, varid, snow))
 varid = varid + 1
 call check (nf90_get_var (ncid, varid, soilw1))
 varid = varid + 1
 call check (nf90_get_var (ncid, varid, soilw2))
 varid = varid + 1
 call check (nf90_get_var (ncid, varid, soilw3))
 varid = varid + 1
 call check (nf90_get_var (ncid, varid, kspp))
 varid = varid + 1
 call check (nf90_get_var (ncid, varid, dbh))
 varid = varid + 1
 call check (nf90_get_var (ncid, varid, hdbh))
 varid = varid + 1
 call check (nf90_get_var (ncid, varid, cfoliage))
 varid = varid + 1
 call check (nf90_get_var (ncid, varid, cfiner))
 varid = varid + 1
 call check (nf90_get_var (ncid, varid, cstore))
 varid = varid + 1
 call check (nf90_get_var (ncid, varid, nfoliage))
 varid = varid + 1
 call check (nf90_get_var (ncid, varid, nbswood))
 varid = varid + 1
 call check (nf90_get_var (ncid, varid, nheart))
 varid = varid + 1
 call check (nf90_get_var (ncid, varid, navail))
 varid = varid + 1
 call check (nf90_get_var (ncid, varid, nfiner))
 varid = varid + 1
 call check (nf90_get_var (ncid, varid, lasa))
 varid = varid + 1
 call check (nf90_get_var (ncid, varid, nitf))
 varid = varid + 1
 call check (nf90_get_var (ncid, varid, hbc))
 call check (nf90_close (ncid))
 ! Do not think this needs to be in RSF because is not a state
 ! variable.
 allocate (lsap  (nind_total))
 allocate (cwood (nind_total))
 do k = 1, nind_total
  ksp = kspp (k)
  !****adf crratio for grass?
  lsap (k) = bark (ksp) * cfoliage (k)
  if (ksp <= 2) then
   !****adf not sure if this is not a state var.
   cwood (k) = lsap (k)
  else
   theight = ah (ksp) * dbh (k) ** bh (ksp)
   tarea = pi * (dbh (k) / 2.0) ** 2
   cwood (k) = stf (ksp) * formf (ksp) * theight * tarea * woodd (ksp)
   cwood (k) = cwood (k) - cstore (k)
  end if
 end do ! k
else
 do kp = 1, nplots
  plot (kp) = kp
 end do
end if
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
if (.NOT. (rsf_in)) then

 ! Set soil C, N, water; plot nind; kspp; dbh; hdbh; cfoliage, cfiner,
 ! cstore, nfoliage, nbswood, nheart, navail, nfiner; lasa; nitf; hbc.
 nind_total = nland_max * nplots * nind_max ! Decrease when have regen.
 allocate (kspp    (nind_total))
 allocate (dbh     (nind_total))
 allocate (hdbh    (nind_total))
 allocate (cfoliage(nind_total))
 allocate (cfiner  (nind_total))
 allocate (cstore  (nind_total))
 allocate (nfoliage(nind_total))
 allocate (nbswood (nind_total))
 allocate (nheart  (nind_total))
 allocate (navail  (nind_total))
 allocate (nfiner  (nind_total))
 allocate (lasa    (nind_total))
 allocate (nitf    (nind_total))
 allocate (hbc     (nind_total))
 allocate (lsap    (nind_total))
 allocate (cwood   (nind_total))
 allocate (alive   (nind_total))
 nplots_total = nland_max * nplots
 allocate (wlittc (nplots_total))
 allocate (wlittn (nplots_total))
 allocate (flittc (nplots_total))
 allocate (flittn (nplots_total))
 allocate (rlittc (nplots_total))
 allocate (rlittn (nplots_total))
 allocate (SWf    (nplots_total))
 alive (:) = 0
 wlittc (:) = zero
 wlittn (:) = zero
 flittc (:) = zero
 flittn (:) = zero
 rlittc (:) = zero
 rlittn (:) = zero
 cd (:,:) = fillvalue_int
 dd (:,:) = fillvalue
 bgs (:,:) = fillvalue_int
 egs (:,:) = fillvalue_int
 
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
  p = 0
  k = 0
  il = 0
  !--------------------------------------------------------------------!
  do j = j1, j2
   !-------------------------------------------------------------------!
   do i = i1, i2
    if (tmp (i, j, 1) /= fillvalue) then
     il = il + 1
     land_index (i,j) = il
     cd (i,j) = 0
     dd (i,j) = zero
     bgs (i,j) = nd + 1
     egs (i,j) = nd + 2
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
      p = p + 1
      p_plot (land_index(i,j),kp) = p
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
      nind (i, j, kp) = 0
      !----------------------------------------------------------------!
      ! Plant grass in plots.
      !----------------------------------------------------------------!
      do ki = 1, 2
       nind (i,j,kp) = nind (i,j,kp) + 1
       k = k + 1
       alive (k) = 1
       ksp = ki
       !ksp = 2 !****adf
       dbh (k) = 0.01 ! LAI of grass
       hbc (k) = 0
       k_ind (land_index(i,j),kp,ki) = k
       kspp (k) = ksp
       cfoliage (k) = dbh (k) * area / sla (ksp)
       !---------------------------------------------------------------!
       ! Grass structure C (kg[C] m-1).
       !---------------------------------------------------------------!
       lsap (k) = crratio (ksp) * cfoliage (k)
       !---------------------------------------------------------------!
       cwood (k) = lsap (k)
       nheart (k) = zero
       navail (k) = 3.0 * 0.04 * cfoliage (k)
       cfiner (k) = rlratio (ksp) * cfoliage (k)
       cstore (k) = zero
       nfoliage (k) = navail (k) * cfoliage (k) / (cfoliage (k) +      &
	              fsr (ksp) * lsap (k) + frr (ksp) * cfiner (k))
       nbswood (k) = navail (k) * fsr (ksp) * lsap (k) / &
                    (cfoliage (k) + fsr (ksp) * lsap (k) +  &
                    frr (ksp) * cfiner (k))
       nfiner (k) = navail (k) * frr (ksp) * cfiner (k) / &
                   (cfoliage (k) + fsr (ksp) * lsap (k) + frr (ksp) * &
		    cfiner (k))
       navail (k) = zero !****adf
       nitf (k) = one
      end do ! ki = 1, 2
      !----------------------------------------------------------------!
      ! Plant trees in plots.
      !----------------------------------------------------------------!
      iregl = ((limp - (nind (i,j,kp) - 2)) / (nspp - 2))
      !iregl = min (iregl, 10)
      if (((nind (i,j,kp) - 2) + iregl * (nspp - 2)) > limp) &
       iregl = iregl - 1
      if (iregl > 0) then
      ki = 2
      !****adfdo ksp = 3, nspp
      do ksp = 3, nspp
       do ireglc = 1, iregl
       nind (i,j,kp) = nind (i,j,kp) + 1
       ki = ki + 1
       !do ki = 3, nind (i,j,kp)
       !---------------------------------------------------------------!
       k = k + 1
       k_ind (land_index(i,j),kp,ki) = k
       alive (k) = 1
       !---------------------------------------------------------------!
       !ksp = 3
       kspp (k) = ksp
       call random_number (ran)
       dbh (k) = idbh * (2.0 * idbhv * ran + (1.0 - idbhv))
       hdbh (k) = zero
       nfp = 0.02
       diamw = dbh (k) * (one - 2.0 * bark (ksp))
       warea = pi * (0.5 * diamw) ** 2
       harea = pi * (0.5 * hdbh (k)) ** 2
       tarea = pi * (dbh (k) / 2.0) ** 2
       wheight = ah (ksp) * diamw ** bh (ksp)
       hheight = ah (ksp) * hdbh (k) ** bh (ksp)
       theight = ah (ksp) * dbh (k) ** bh (ksp)
       wwood = stf (ksp) * formf (ksp) * wheight * &
               warea * woodd (ksp)
       hwood = stf (ksp) * formf (ksp) * hheight * &
               harea * woodd (ksp)
       cwood (k) = stf (ksp) * formf (ksp) * theight * tarea * &
                   woodd (ksp)
       swood = wwood - hwood
       saparea = warea - harea
       !---------------------------------------------------------------!
       ! Foliage area to sapwood area ratio (-).
       !---------------------------------------------------------------!
       lasa (k) = lsave (ksp)
       !---------------------------------------------------------------!
       ! Foliage C of each individual (kg[C] ind-1).
       !---------------------------------------------------------------!
       cfoliage (k) = lasa (k) * saparea / sla (ksp)
       !---------------------------------------------------------------!
       ! Fine root C (kg[C] ind-1).
       !---------------------------------------------------------------!
       cfiner (k) = rlratio (ksp) * cfoliage (k)
       !---------------------------------------------------------------!
       lsap (k) = live (ksp) * swood
       cstore (k) = storef (ksp) * lsap (k)
       !---------------------------------------------------------------!
       ! Foliage N (kg[N] ind-1).
       !---------------------------------------------------------------!
       nfoliage (k) = nfp * 2.0 * cfoliage (k)
       !---------------------------------------------------------------!
       ! Stem and bark N (kg[N] ind-1).
       !---------------------------------------------------------------!
       nbswood (k) = swood / ((cfoliage (k) / nfoliage (k)) / fsr (ksp))
       !---------------------------------------------------------------!
       ! Heartwood N (kg[C] ind-1).
       !---------------------------------------------------------------!
       nheart (k)  = hwood / ((cfoliage (k) / nfoliage (k)) / fsr (ksp))
       !---------------------------------------------------------------!
       ! Fine root N (kg[N] ind-1).
       !---------------------------------------------------------------!
       nfiner (k) = cfiner (k) / ((cfoliage (k) / nfoliage (k)) / &
                    frr (ksp))
       !---------------------------------------------------------------!
       ! Frost effect (fraction).
       !---------------------------------------------------------------!
       nitf (k) = one
       !---------------------------------------------------------------!
       ! N storage (kg[N] ind-1).
       !---------------------------------------------------------------!
       navail (k) = zero
       !---------------------------------------------------------------!
      !end do ! ki
       end do ! il
      end do ! ksp
      end if
      !----------------------------------------------------------------!
      ! Total soil C (kg[C] m-2).
      !----------------------------------------------------------------!
      soilC (i,j,kp) = Cu (i,j,kp) + Cm (i,j,kp) + Cv (i,j,kp) + &
                       Cn (i,j,kp) + Ca (i,j,kp) + Cs (i,j,kp) + &
                       Cpa (i,j,kp)
      !----------------------------------------------------------------!
      ! Total soil water holding capacity (Eqn. (1) of FW00; m).
      !----------------------------------------------------------------!
      swct (i,j,kp) = 0.213 + 0.00227 * soilC (i,j,kp)
      !----------------------------------------------------------------!
      ! Soil water holding capacities of layers (m).
      ! d_one is max. swc of top layer for roots.
      !----------------------------------------------------------------!
      if (swct (i,j,kp) <= 0.05) then
       swc1 (i,j,kp) = swct (i,j,kp)
       swc2 (i,j,kp) = 0.0
       swc3 (i,j,kp) = 0.0
      else
       if (swct (i,j,kp) <= d_one) then
        swc1 (i,j,kp) = 0.05
        swc2 (i,j,kp) = swct (i,j,kp) - 0.05
        swc3 (i,j,kp) = 0.0
       else
        swc1 (i,j,kp) = 0.05
        swc2 (i,j,kp) = d_one - 0.05
        swc3 (i,j,kp) = swct (i,j,kp) - d_one
       end if
      end if
      !----------------------------------------------------------------!
      ! Snowpack (m).
      !----------------------------------------------------------------!
      snow (i,j,kp) = 0.0
      !----------------------------------------------------------------!
      ! Soil water (m).
      !----------------------------------------------------------------!
      soilw1 (i,j,kp) = 0.5 * swc1 (i,j,kp)
      soilw2 (i,j,kp) = 0.5 * swc2 (i,j,kp)
      soilw3 (i,j,kp) = 0.5 * swc3 (i,j,kp)
      !----------------------------------------------------------------!
      !if (local) write (*, *) kp, isc (i, j), Cm (i, j, kp), &
      ! Cv (i, j, kp), Cn (i, j, kp), Ca (i, j, kp), Cs (i, j, kp), &
      ! Cpa (i, j, kp)
      !if (local) write (*, *) kp, isc (i, j), Nm (i, j, kp), &
      ! Nv (i, j, kp), Nn (i, j, kp), Na (i, j, kp), Ns (i, j, kp), &
      ! Npa (i, j, kp)
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
   call check (nf90_def_var (ncid, "Soil_C", nf90_float, &
               dimids_two, varid))
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
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Initialised either from rsf or standard values. Now set some
! derived initial quantities.
!----------------------------------------------------------------------!
allocate (skzg   (nlon,nlat,nplots))
allocate (fPARtg (nlon,nlat,nplots))
allocate (skSWzg (nlon,nlat,nplots))
allocate (fSWtg  (nlon,nlat,nplots))
allocate (gmax   (nind_total))
allocate (cfact  (nind_total))
allocate (isfact (nind_total))
allocate (et_cat (nind_total))
allocate (gc     (nind_total)) ! really only needed within grid-box
allocate (ratiol (nind_total))
allocate (ipfact (nind_total))
allocate (jmaxfn (nind_total))
allocate (height (nind_total))
allocate (kzg    (nind_total))
allocate (kSWzg  (nind_total))
allocate (fPARiig(nind_total))
allocate (fSWiig (nind_total))
allocate (fPAR   (nind_total))
allocate (fSW    (nind_total))
allocate (rcfoliage(nind_total))
allocate (rcwood   (nind_total))
allocate (rcfiner  (nind_total))
allocate (rnfoliage(nind_total))
allocate (rnbswood (nind_total))
allocate (rnfiner  (nind_total))
allocate (rlold    (nind_total))
allocate (farea    (nind_total))
allocate (kz       (mh+1,nind_max))
allocate (kSWz     (mh+1,nind_max))
!----------------------------------------------------------------------!
do j = j1, j2
 !---------------------------------------------------------------------!
 do i = i1, i2
  if (tmp (i, j, 1) /= fillvalue) then
   do kp = 1, nplots
    !------------------------------------------------------------------!
    p = p_plot (land_index(i,j),kp)
    !------------------------------------------------------------------!
    ! Total soil C (kg[C] m-2).
    !------------------------------------------------------------------!
    soilC (i,j,kp) = Cu (i,j,kp) + Cm (i,j,kp) + Cv (i,j,kp) + &
                     Cn (i,j,kp) + Ca (i,j,kp) + Cs (i,j,kp) + &
                     Cpa (i,j,kp)
    !------------------------------------------------------------------!
    ! Total soil N (kg[N] m-2).
    !------------------------------------------------------------------!
    soilN (i,j,kp) = Nu (i,j,kp) + Nm (i,j,kp) + Nv (i,j,kp) + &
                     Nn (i,j,kp) + Na (i,j,kp) + Ns (i,j,kp) + &
                     Npa (i,j,kp)
    !****adf include snmin?
    !------------------------------------------------------------------!
    ! Total soil water holding capacity (Eqn. (1) of FW00; m).
    !------------------------------------------------------------------!
    swct (i,j,kp) = 0.213 + 0.00227 * soilC (i, j, kp)
    !------------------------------------------------------------------!
    ! Soil water holding capacities of layers (m).
    ! d_one is max. swc of top layer for roots.
    !------------------------------------------------------------------!
    if (swct (i,j,kp) <= 0.05) then
     swc1 (i,j,kp) = swct (i,j,kp)
     swc2 (i,j,kp) = zero
     swc3 (i,j,kp) = zero
    else
     if (swct (i,j,kp) <= d_one) then
       swc1 (i,j,kp) = 0.05
       swc2 (i,j,kp) = swct (i,j,kp) - 0.05
       swc3 (i,j,kp) = zero
     else
       swc1 (i,j,kp) = 0.05
       swc2 (i,j,kp) = d_one - 0.05
       swc3 (i,j,kp) = swct (i,j,kp) - d_one
     end if
    end if
    !****following not needed as loops over sites afterwards
    !****but is same as orig, so maybe needed for first site!
    !------------------------------------------------------------------!
    ! Soil water potential top-layer roots (MPa).
    !------------------------------------------------------------------!
    sow1 = soilw1 (i,j,kp) + soilw2 (i,j,kp)
    sc1  = swc1   (i,j,kp) + swc2   (i,j,kp)
    if (sow1 > eps) then
     if (sc1 > eps) then
      sw1 = swpmax * (1.0 / ((sow1 / sc1) ** bsoil))
     else
      sw1 = -100.0
     end if
    else
     sw1 = -100.0
    end if
    !------------------------------------------------------------------!
    ! Soil water potential for bottom-layer roots (MPa).
    !-----------------------------------------------------------------!
    if (soilw3 (i,j,kp) > eps) then
     if (swc3 (i,j,kp) > eps) then
      sw2 = swpmax * (1.0 / ((soilw3 (i,j,kp) / &
	    swc3 (i,j,kp)) ** bsoil))
     else
      sw2 = -100.0
     end if
    else
     sw2 = -100.0
    end if
    !------------------------------------------------------------------!
    ! Soil water potential for grass (MPa).
    !------------------------------------------------------------------!
    swp1 (kp) = sw1
    !------------------------------------------------------------------!
    ! Soil water potential for trees (MPa).
    !------------------------------------------------------------------!
    if (sw1 > sw2) then
     ! All transpiration from top two layers.
     swp2 (kp) = sw1
     tr_rat = one
    else
     ! All transpiration from bottom layer.
     swp2 (kp) = sw2
     tr_rat = zero
    end if
    !****above not needed? as loops over sites afterwards
    !****but is same as orig, so maybe needed for first site!
    !------------------------------------------------------------------!
    ! Factor to calculate absorbed SW by whole plot (ratio).
    !------------------------------------------------------------------!
    SWf (p) = zero
    !------------------------------------------------------------------!
    ! Set sums of foliage area * extinction coefficient in each layer
    ! to zero.
    !------------------------------------------------------------------!
    do m = (mh + 1), 1, -1
     skz   (m) = zero
     skSWz (m) = zero
    end do
    ! Foliage area in each layer.
    do ki = 1, nind (i,j,kp)
     k = k_ind (land_index(i,j),kp,ki)
     ksp = kspp (k)
     ! Set fraction of PAR on top layer absorbed by individual to zero.
     fPAR (k) = zero
     ! Set fraction of SW on top layer absorbed by individual to zero.
     fSW (k) = zero
     do m = (mh + 1), 1, -1
      kz   (m,ki) = zero
      kSWz (m,ki) = zero
     end do
     farea (k) = cfoliage (k) * sla (ksp)
     rlold (k) = farea (k)
     !-----------------------------------------------------------------!
     ! Height of crown (m) allometrically from dbh (m).
     !----------------------------------------------------------------!
     if (dbh (k) > eps) then
      ht = ah (ksp) * dbh (k) ** bh (ksp)
     else
      ht = zero
     end if
     height (k) = nint (ht + 0.5)
     if (height (k) < 1) height (k) = 1
     if (height (k) > mh) height (k) = mh
     !-----------------------------------------------------------------!
     ! Height to base of crown (m).
     !-----------------------------------------------------------------!
     hbc (k) = min (hbc (k), (height (k) - 1))
     hbc (k) = max (hbc (k), 0)
     !-----------------------------------------------------------------!
     ! Number of 1 m layers in crown (n).
     nlayers = height (k) - hbc (k)
     nlayers = max (1, nlayers)
     if (farea (k) > eps) then
      ! fd is LAI of individual in each layer.
      fd = farea (k) / (float (nlayers) * area)
      ! Factors to make following calculation easier.
      kf   = kpar (ksp) * fd
      kSWf = ksw  (ksp) * fd
      ! Loop down through canopy, calculate LAI in each layer.
      ! m refers to number of layer below height m (m).
      do m = height (k), (hbc (k) + 1), -1
       ! Individual KI,j * Zf,i,j in layer.
       kz   (m,ki) = kf
       kSWz (m,ki) = kSWf
       ! Sum of KI,j * Zf,i,j in layer.
       skz   (m) = skz   (m) + kz   (m,ki)
       skSWz (m) = skSWz (m) + kSWz (m,ki)
      end do ! m
     end if
    end do ! ki
    !------------------------------------------------------------------!
    ! Fraction of PAR penetrating layer above top of plot.
    !------------------------------------------------------------------!
    eskz  (mh+1) = one
    fPARt (mh+1) = one
    !------------------------------------------------------------------!
    ! Fraction of SW radiation penetrating layer above top of plot
    ! (fraction).
    !------------------------------------------------------------------!
    eskSWz (mh+1) = one
    fSWt   (mh+1) = one
    !------------------------------------------------------------------!
    ! Proportion of light at top of each layer absorbed by that layer.
    !------------------------------------------------------------------!
    do m = mh, 1, -1
     ! Fraction of PAR penetrating layer.
     eskz (m) = exp (-1.0 * skz (m))
     ! Fraction of SW penetrating layer.
     eskSWz (m) = exp (-1.0 * skSWz (m))
     ! Assume tiny bit of radiation intercepted by structure (this stops
     ! eskz and eskSWz being 1 at very low foliage areas, and causing
     ! later divide by zero).
     eskz   (m) = min (0.99999, eskz (m))
     eskSWz (m) = min (0.99999, eskSWz (m))
     ! Total fraction of PAR incident on plot incident on layer,
     ! equal to fraction penetrating layer above multiplied by fraction
     ! of PAR at top of plot incident on layer above.
     fPARt (m) = eskz (m+1) * fPARt (m+1)
     ! Total fraction of PAR incident on plot incident on layer,
     ! equal to fraction penetrating layer above multiplied by fraction
     ! of PAR at top of plot incident on layer above.
     fSWt (m) = eskSWz (m+1) * fSWt (m+1)
     ! Total fraction of radiation at top of plot absorbed in layer.
     fSWi (m) = (1.0 - eskSWz (m)) * fSWt (m)
     ! Total fraction of radiation at top of plot absorbed in layer.
     fPARi (m) = (1.0 - eskz (m)) * fPARt (m)
     ! Individual absorptions.
     do ki = 1, nind (i,j,kp)
      k = k_ind (land_index(i,j),kp,ki)
      ! Save species number for easier coding.
      ksp = kspp (k)
      ! Check if foliage in layer.
      if (skz (m) > eps) then
       ! Fraction of absorption in layer due to individual.
       frac  = kz   (m,ki) / skz   (m)
       fracs = kSWz (m,ki) / skSWz (m)
       ! Sum total fraction absorbed by individual relative to top of
       ! plot.
       fPAR (k) = fPAR (k) + frac * fPARi (m)
       ! Sum total fraction absorbed by individual relative to top of
       ! plot.
       fSW (k) = fSW (k) + fracs * fSWi (m)
       ! Sum factor to calculate absorbed SW by foliage from incident
       ! SW on plot.
       SWf (p) = SWf (p) + (one - rhos (ksp)) * fracs * fSWi (m)
      end if
     end do ! ki
    end do ! m
    !stop
    !------------------------------------------------------------------!
    laip (i,j,kp) = zero
    !------------------------------------------------------------------!
    ! Canopy net photosynthesis and conductance.
    !------------------------------------------------------------------!
    do ki = 1, nind (i,j,kp)
     k = k_ind (land_index(i,j),kp,ki)
     ksp = kspp (k)
     !-----------------------------------------------------------------!
     ! Save required C and N tree compartment sizes.
     !-----------------------------------------------------------------!
     if (ksp > 2) then ! tree
      rcfoliage (k) = cfoliage (k)
      rcwood    (k) = cwood    (k)
      rcfiner   (k) = cfiner   (k)
      rnfoliage (k) = nfoliage (k)
      rnbswood  (k) = nbswood  (k)
      rnfiner   (k) = nfiner   (k)
     end if
     !-----------------------------------------------------------------!
     laip (i,j,kp) = laip (i,j,kp) + farea (k)
     !-----------------------------------------------------------------!
    end do ! ki
    laip (i,j,kp) = laip (i,j,kp) / area
   end do ! kp
   phenf (:) = 1
   call radiation !****replaces some of above?
  end if ! tmp (i,j,1) /= fillvalue
 end do ! i
end do ! j
outflow (:) = zero
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

end subroutine initl
