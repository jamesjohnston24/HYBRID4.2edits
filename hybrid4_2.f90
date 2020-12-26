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
integer, parameter :: nlon     =  720
integer, parameter :: nlat     =  360
integer, parameter :: ntimes   = 1460
integer, parameter :: nplots   =   10
integer, parameter :: nind_max =  100
integer, parameter :: nspp     =    4
integer, parameter :: mh       =  110
real, parameter :: dt = 6.0 * 60. * 60.0 ! s timestep-1
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
real, parameter :: iscv = 0.2 ! Range of initial soil C (+/-fraction)
real, parameter :: area = 200.0 ! Plot area (m2)
!                                   C3GR C4GR BRCD NLEV
real, dimension (nspp) :: rhos  = (/0.20,0.20,  0.20,  0.11/)
real, dimension (nspp) :: ksw   = (/0.48,0.48,  0.48,  0.37/)
real, dimension (nspp) :: lsave = (/0.00,0.00,4167.0,3333.0/)
real, dimension (nspp) :: sla   = (/36.0,36.0,  36.0,  12.0/)
real, parameter :: Tc_local = 0.3 ! Clay fraction?
real, parameter :: Ts_local = 0.3 ! Sand fraction?
real, parameter :: plf = 0.01 + 0.04 * Ts_local
real, parameter :: swpmax = -0.033 ! MPa
real, parameter :: bsoil = 5.0
! Biomass C content (fraction).
real, parameter :: w = 0.45
! Lignin to biomass ratio in leaf litter.
real, parameter :: lfl = 0.20
! Lignin to biomass ratio in root litter.
real, parameter :: lrl = 0.16
real, parameter :: texture = 0.5 ! Soil texture (units?)
! N:C ratio of surface structural litter.
real, parameter :: vu =  1.0 / 150.0
! N:C ratio of soil structural litter.
real, parameter :: vv =  1.0 / 150.0
! N:C ratio of active soil organic matter pool.
real, parameter :: va =  1.0 /  15.0
! N:C ratio of slow soil organic matter pool.
real, parameter :: vs =  1.0 /  20.0
! N:C ratio of soil passive organic matter pool.
real, parameter :: vpa = 1.0 / 10.0
real, parameter :: emf = 0.05 ! N emission (fraction).
! Partition coefficients.
real, parameter :: pau = 0.55 * (one - lfl)
real, parameter :: psu = 0.7 * lfl
real, parameter :: pav = 0.45 * (one - lfl)
real, parameter :: psv = 0.7 * lfl
real, parameter :: pam = 0.45
real, parameter :: pan = 0.45
real, parameter :: psab = 0.996 - (0.85 - 0.68 * texture)
real, parameter :: ppa = 0.004
real, parameter :: pas = 0.42
real, parameter :: pps = 0.03
real, parameter :: pap = 0.45
! Decay rates base values from Comins & McMurtrie (1993), but
! doubled to allow for inclusion of the soil water factor (/day).
real, parameter :: d1pb = 2.0 * 0.076 * exp (-3.0 * lfl)      / 7.0
real, parameter :: d2pb = 2.0 * 0.28                          / 7.0
real, parameter :: d3pb = 2.0 * 0.094 * exp (-3.0 * lrl)      / 7.0
real, parameter :: d4pb = 2.0 * 0.35                          / 7.0
real, parameter :: d5pb = 2.0 * 0.14 * (one - 0.75 * texture) / 7.0
real, parameter :: d6pb = 2.0 * 0.0038                        / 7.0
real, parameter :: d7pb = 2.0 * 0.00013                       / 7.0
real, parameter :: Nad = 0.0 ! N addition (kg[N] m-2 day-1).
real, parameter :: nfx = 0.0 ! N fixation (kg[N] m-2 day-1).
! Pre-industrial N deposition at 0 m ppt (kg[N] ha-2 yr-1 ->
! kg[N] m-2 day-1).
real, parameter :: ndepi = 10.0 / (10000.0 * 365.0)
! d_one is max. swc of top layer for roots (m).
real, parameter :: d_one = 0.2
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
integer :: varid_kspp
integer :: varid_cfoliage
integer :: varid_lasa
integer :: ncid
integer :: lon_dimid, lat_dimid, plot_dimid, ind_dimid
integer :: lon_varid, lat_varid, plot_varid, ind_varid
integer, dimension (2) :: dimids_two
integer, dimension (3) :: dimids_three
integer, dimension (4) :: dimids_four
integer :: i,ii,i1,i2
integer :: j,jj,j1,j2
integer :: k
integer :: kp
integer :: ki
integer :: ksp
integer :: m
integer :: nland
integer, allocatable, dimension (:) :: plot ! No. of plot
real, allocatable, dimension (:,:) :: icwtr_qd ! Ice/water fraction
real, allocatable, dimension (:,:) :: icwtr ! Ice/water fraction
real, allocatable, dimension (:,:,:) :: tmp ! K
real, allocatable, dimension (:,:,:) :: dswrf ! J m-2
real, allocatable, dimension (:,:,:) :: pre ! mm 6hr-1
real, allocatable, dimension (:,:,:) :: spfh ! kg kg-1
real, allocatable, dimension (:,:,:) :: pres ! Pa
real, allocatable, dimension (:) :: lat ! degrees north
real, allocatable, dimension (:) :: lon ! degrees east
real, allocatable, dimension (:,:) :: isc ! Initial soil C (kg[C] m-2)
real, allocatable, dimension (:,:) :: isc_grid ! Init soil C (kg[C] m-2)
real, allocatable, dimension (:,:) :: larea_qd ! Grid-box area (km2)
real, allocatable, dimension (:,:) :: larea ! Grid-box area (km2)
real, allocatable, dimension (:,:,:) :: Cm
real, allocatable, dimension (:,:,:) :: Cu
real, allocatable, dimension (:,:,:) :: Cn
real, allocatable, dimension (:,:,:) :: Cv
real, allocatable, dimension (:,:,:) :: Ca
real, allocatable, dimension (:,:,:) :: Cs
real, allocatable, dimension (:,:,:) :: Cpa
real, allocatable, dimension (:,:,:) :: soilC
real, allocatable, dimension (:,:,:) :: Nu
real, allocatable, dimension (:,:,:) :: Nm
real, allocatable, dimension (:,:,:) :: Nv
real, allocatable, dimension (:,:,:) :: Nn
real, allocatable, dimension (:,:,:) :: Na
real, allocatable, dimension (:,:,:) :: Ns
real, allocatable, dimension (:,:,:) :: Npa
real, allocatable, dimension (:,:,:) :: snmin
real, allocatable, dimension (:,:,:) :: soilN
integer, allocatable, dimension (:,:,:) :: nind
! Soil water holding capacity (m).
real, allocatable, dimension (:,:,:) :: swct
real, allocatable, dimension (:,:,:) :: soilw1
real, allocatable, dimension (:,:,:) :: soilw2
real, allocatable, dimension (:,:,:) :: soilw3
integer, allocatable, dimension (:,:,:,:) :: kspp
real   , allocatable, dimension (:,:,:,:) :: cfoliage
real   , allocatable, dimension (:,:,:,:) :: lasa
real   , allocatable, dimension (:,:,:) :: swc1
real   , allocatable, dimension (:,:,:) :: swc2
real   , allocatable, dimension (:,:,:) :: swc3
real   , allocatable, dimension (:,:,:,:) :: isfact
integer, allocatable, dimension (:,:,:,:) :: height
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
real :: radsw ! 
real :: isc_total ! Total global soil C (Pg[C])
real :: ran ! Random number (0-1)
real :: Noise ! Noise on initial soil C (ratio)
real :: vm = 0.07 ! N:C ratio of surface metabolic litter
real :: vn = 0.07 ! N:C ratio of soil metabolic litter
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
real :: wfps ! Water-filled pore space (%)
real :: outflow
real :: h2o_30
real :: pleach
real :: pmf
real :: puf
real :: pnr
real :: pvr
real :: psa
real :: d1p,d2p,d3p,d4p,d5p,d6p,d7p
real :: d1C,d2C,d3C,d4C,d5C,d6C,d7C
real :: dCu,dCm,dCv,dCn,dCa,dCs,dCpa,dCle
real :: Ju,Jv,Ja,Js,Jpa
real :: NmfNmw,Nnr
real :: dNu,dNm,dNv,dNn,dNa,dNs,dNpa
real :: temp,C0,N0,Cleach,Nmin,U,Ndepo,nf
real :: sresp ! Soil respiration (kg[C] m-1 day-1)
real :: st ! SW radiation absorbed in upper crown layer (W m-2)
real, dimension (mh+1) :: fSWt
real :: fst,fdt
real :: swp1,swp2,swp3
real :: sow1
real :: sc1
real :: sw1,sw2
real :: tr_rat
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
allocate (swct   (nlon,nlat,nplots))
allocate (soilw1 (nlon,nlat,nplots))
allocate (soilw2 (nlon,nlat,nplots))
allocate (soilw3 (nlon,nlat,nplots))
allocate (kspp     (nlon,nlat,nplots,nind_max))
allocate (cfoliage (nlon,nlat,nplots,nind_max))
allocate (lasa     (nlon,nlat,nplots,nind_max))
allocate (swc1 (nlon,nlat,nplots))
allocate (swc2 (nlon,nlat,nplots))
allocate (swc3 (nlon,nlat,nplots))
allocate (isfact (nlon,nlat,nplots,nind_max))
allocate (height (nlon,nlat,nplots,nind_max))
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
kspp     (:,:,:,:) = zero
cfoliage (:,:,:,:) = zero
lasa     (:,:,:,:) = zero
height   (:,:,:,:) = 0
isfact   (:,:,:,:) = zero
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
 open (21,file='soil.txt',status='unknown')
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
 varid = varid + 1
 call check (nf90_get_var (ncid, varid, lat))
 varid = varid + 1
 call check (nf90_get_var (ncid, varid, plot))
 varid = varid + 2
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
 call check (nf90_get_var (ncid, varid, soilw1))
 varid = varid + 1
 call check (nf90_get_var (ncid, varid, soilw2))
 varid = varid + 1
 call check (nf90_get_var (ncid, varid, soilw3))
 varid = varid + 1
 call check (nf90_get_var (ncid, varid, kspp))
 varid = varid + 1
 call check (nf90_get_var (ncid, varid, cfoliage))
 varid = varid + 1
 call check (nf90_get_var (ncid, varid, lasa))
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
if (.NOT. (rsf_in)) then
! Set soil C, N, water; plot nind; cfoliage; lasa.
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
      kspp (i,j,kp,1) = 1
      kspp (i,j,kp,2) = 2
      !****adf
      kspp (i,j,kp,3) = 3
      !****adf
      !----------------------------------------------------------------!
      do ki = 1, nind (i,j,kp)
       !****adf
       !---------------------------------------------------------------!
       ksp = kspp (i,j,kp,ki)
       !---------------------------------------------------------------!
       ! Foliage area to sapwood area ratio (-).
       !---------------------------------------------------------------!
       lasa (i,j,kp,ki) = lsave (ksp)
       !---------------------------------------------------------------!
       ! Foliage C of each individual (kg[C] ind-1).
       !---------------------------------------------------------------!
       cfoliage (i,j,kp,ki) = (area * 3.0 / sla (ksp)) / 3.0
       !---------------------------------------------------------------!
       !****adf
      end do ! ki
      !----------------------------------------------------------------!
      ! Total soil C (kg[C] m-2).
      !----------------------------------------------------------------!
      soilC (i,j,kp) = Cu (i,j,kp) + Cm (i,j,kp) + Cv (i,j,kp) + &
                       Cn (i,j,kp) + Ca (i,j,kp) + Cs (i,j,kp) + &
                       Cpa (i,j,kp)
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
      ! Soil water (m).
      !----------------------------------------------------------------!
      soilw1 (i,j,kp) = 0.5 * swc1 (i,j,kp)
      soilw2 (i,j,kp) = 0.5 * swc2 (i,j,kp)
      soilw3 (i,j,kp) = 0.5 * swc3 (i,j,kp)
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
do j = j1, j2
 !---------------------------------------------------------------------!
 do i = i1, i2
  if (tmp (i, j, 1) /= fillvalue) then
   do kp = 1, nplots
    ! Total soil C (kg[C] m-2).
    soilC (i,j,kp) = Cu (i,j,kp) + Cm (i,j,kp) + Cv (i,j,kp) + &
                     Cn (i,j,kp) + Ca (i,j,kp) + Cs (i,j,kp) + &
                     Cpa (i,j,kp)
    ! Total soil N (kg[N] m-2).
    soilN (i,j,kp) = Nu (i,j,kp) + Nm (i,j,kp) + Nv (i,j,kp) + &
                     Nn (i,j,kp) + Na (i,j,kp) + Ns (i,j,kp) + &
                     Npa (i,j,kp)
    !------------------------------------------------------------------!
    ! Soil water holding capacity (Eqn. (1) of FW00; m).
    !------------------------------------------------------------------!
    swct (i,j,kp) = 0.213 + 0.00227 * soilC (i, j, kp)
    !------------------------------------------------------------------!
    ! Soil water holding capacities of layers (m).
    ! d_one is max. swc of top layer for roots.
    !------------------------------------------------------------------!
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
    do ki = 1, nind (i,j,kp)
     ksp = kspp (i,j,kp,ki)
     !****adf
     !-----------------------------------------------------------------!
     ! Plant height (m).
     !-----------------------------------------------------------------!
     height (i,j,kp,ki) = 1
     !-----------------------------------------------------------------!
     ! Fraction of SW radiation penetrating layer above (fraction).
     !-----------------------------------------------------------------!
     fSWt (mh+1) = 1.0
     !-----------------------------------------------------------------!
     do m = mh, 1, -1
      fSWt (m) = 1.0
     end do ! m
     !-----------------------------------------------------------------!
     !****adf
     !-----------------------------------------------------------------!
     ! Factor to calculate absorbed solar radiation.
     !-----------------------------------------------------------------!
     isfact (i,j,kp,ki) = (1.0 - rhos (ksp)) * &
                           fSWt (height (i,j,kp,ki)) * ksw (ksp)
     !-----------------------------------------------------------------!
    end do ! ki
   end do ! kp
  end if ! tmp (i,j,1) /= fillvalue
 end do ! i
end do ! j
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Loop through gridboxes and integrate state variables at land points.
! Climate files start at 1901, run through 2019.
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
 ! Read temperatures (K).
 varid = 4
 call check (nf90_get_var (ncid, varid, tmp))
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
 nland = 0
 do j = j1, j2
  do i = i1, i2
   ! Land?
   if (tmp (i, j, 1) /= fillvalue) then
    nland = nland + 1
    do it = 1, ntimes
     !-----------------------------------------------------------------!
     ! Downwelling solar radiation (W m-2).
     !-----------------------------------------------------------------!
     radsw = dswrf (i,j,it) / dt
     !-----------------------------------------------------------------!
     ! Soil temperature (oC).
     !-----------------------------------------------------------------!
     tsoil = sum (tmp (i,j,it-3:it)) / 4.0 - tf
     !-----------------------------------------------------------------!
     do kp = 1, nplots
      !----------------------------------------------------------------!
      ! Soil water potentials (after Johnson et al., 1991) (MPa).
      !----------------------------------------------------------------!
      if (tsoil > zero) then
       ! Soil water potential top-layer roots (MPa).
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
       !---------------------------------------------------------------!
       ! Soil water potential for bottom-layer roots (MPa).
       !---------------------------------------------------------------!
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
       !---------------------------------------------------------------!
       ! Soil water potential for grass (MPa).
       !---------------------------------------------------------------!
       swp1 = sw1
       !---------------------------------------------------------------!
       ! Soil water potential for trees (MPa).
       !---------------------------------------------------------------!
       if (sw1 > sw2) then
        ! All transpiration from top two layers.
        swp2 = sw1
	tr_rat = one
       else
        ! All transpiration from bottom layer.
        swp2 = sw2
	tr_rat = zero
       end if
       !---------------------------------------------------------------!
      else
       !---------------------------------------------------------------!
       ! Soil frozen.
       !---------------------------------------------------------------!
       swp1 = -1.48
       swp2 = -1.48
       tr_rat = zero
       !---------------------------------------------------------------!
      end if ! tsoil > zero
      !----------------------------------------------------------------!
      do ki = 1, nind (i,j,kp)
       !---------------------------------------------------------------!
       ! SW radiation absorbed in upper crown layer (W m-2).
       !---------------------------------------------------------------!
       st = radsw * isfact (i,j,kp,ki)
       !---------------------------------------------------------------!
       ! Stomatal conductance solar-radiation factor (fraction).
       !---------------------------------------------------------------!
       fst = 1.1044 * st / (st + 104.4)
       !---------------------------------------------------------------!
       ! Stomatal conductance soil-water factor (fraction).
       !---------------------------------------------------------------!
       if ((ki == 1) .or. (ki == 2)) then
        ! Soil water factor for grass (fraction).
        if (swp1 > (-1.5)) fdt = 1.155 + 0.77 * &
	 (swp1 - 0.015 * height (i,j,kp,ki))
         if (swp1 >  (-0.2)) fdt = one
         if (swp1 <= (-1.5)) fdt = zero
       else
        if (swp2 > (-1.5)) fdt = 1.155 + 0.77 * &
	 (swp2 - (lsave (ksp) / lasa (i,j,kp,ki)) * 0.005 * &
	 height (i,j,kp,ki) - 0.01 * height (i,j,kp,ki))
        if (swp2 >  (-0.2)) fdt = one
        if (swp2 <= (-1.5)) fdt = zero
       end if
       !---------------------------------------------------------------!
      end do
     end do
     !****adf
     outflow = 0.001 ! Daily outflow (m).
     !****adf
     ! End of day?
     if (mod (it, 4) == zero) then
      !****adf
      nf = ndepi
      !****adf
      Ndepo = nf + 4.0e-3 * sum (pre (i, j, it-3:it)) / 1000.0
      do kp = 1, nplots
       wlittc (kp) = zero
       wlittn (kp) = zero
       flittc (kp) = zero
       flittn (kp) = zero
       rlittc (kp) = zero
       rlittn (kp) = zero
       do ki = 1, nind (i,j,kp)
        !****adf
	! Assume woody mass is 5 kg[C] m-2, 10% annual turnover;
	! 2% N in foliage and fine roots; 1% N in wood.
	! Fine root and foliage equal.
	ksp = kspp (i,j,kp,ki)
	if (ksp >= 3) then
         wlitterc = area * 0.1 * 5.0 / 365.0
         wlittern = 2.0 * 0.01 * wlitterc
	else
	 wlitterc = zero
	 wlittern = zero
	end if
        flitterc = cfoliage (i,j,kp,ki) / 365.0
        flittern = 2.0 * 0.02 * flitterc
        rlitterc = flitterc
        rlittern = flittern
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
       !---------------------------------------------------------------!
       ! Soil water holding capacity (Eqn. (1) of FW00; m).
       !---------------------------------------------------------------!
       swct (i,j,kp) = 0.213 + 0.00227 * soilC (i, j, kp)
       !---------------------------------------------------------------!
       ! Soil water holding capacities of layers (m).
       ! d_one is max. swc of top layer for roots.
       !---------------------------------------------------------------!
       if (swct (i,j,kp) <= 0.05) then
        swc1 (i,j,kp) = swct (i,j,kp)
        swc2 (i,j,kp) = 0.0
        swc3 (i,j,kp) = 0.0
       else
        if (swct (i,j,kp) <= d_one) THEN
          swc1 (i,j,kp) = 0.05
          swc2 (i,j,kp) = swct (i,j,kp) - 0.05
          swc3 (i,j,kp) = 0.0
        else
          swc1 (i,j,kp) = 0.05
          swc2 (i,j,kp) = d_one - 0.05
          swc3 (i,j,kp) = swct (i,j,kp) - d_one
        end if
       end if
       !---------------------------------------------------------------!
       ! Convert soil water to water-filled pore space.
       ! Assumes micro-pore space = swc and macro-pore space = 42%
       ! saturation content (from TEM, for loam; Raich et al., 1991).
       !---------------------------------------------------------------!
       if (swct (i,j,kp) > eps) then
        wfps = 100.0 * (soilw1 (i, j, kp) + soilw2 (i, j, kp) + &
                        soilw3 (i, j, kp)) / (1.7241 * swct (i,j,kp))
	wfps = min (100.0, wfps)
       else
        wfps = 0.00001
       end if
       !---------------------------------------------------------------!
       ! Soil-water decomposition modifer (fraction). Equations are 
       ! fitted to Fig. 8 of Williams et al. (1992).
       !---------------------------------------------------------------!
       if (wfps < 60.0) then
        em_soil = exp ((wfps - 60.0) ** 2 / (-800.0))
       else
        em_soil = 0.000371 * wfps ** 2 - 0.0748 * wfps + 4.13
       endif
       !---------------------------------------------------------------!
       em_soil = max (zero, em_soil)
       em_soil = min (one , em_soil)
       !---------------------------------------------------------------!
       ! Overall decomposition modifier.
       !---------------------------------------------------------------!
       ev = et_soil * em_soil
       !---------------------------------------------------------------!
       ! Soil microbe leach fraction.
       !---------------------------------------------------------------!
       h2o_30 = outflow
       pleach = 0.05555 * h2o_30 * plf
       ! Partitioning coefficients.
       if (flittn (kp) > eps) then
        pmf = 0.85 - 0.018 * lfl / (w * flittn (kp) / flittc (kp))
	pmf = max (zero, pmf)
       else
        pmf = zero
       end if
       puf = one - pmf
       if (rlittn (kp) > eps) then
        pnr = 0.85 - 0.018 * lrl / (w * rlittn (kp) / rlittc (kp))
	pnr = max (zero, pnr)
       else
        pnr = zero
       end if
       pvr = one - pnr
       psa = psab - pleach
       ! Decay rates (/day).
       d1p = ev * d1pb
       d2p = ev * d2pb
       d3p = ev * d3pb
       d4p = ev * d4pb
       d5p = ev * d5pb
       d6p = ev * d6pb
       d7p = ev * d7pb
       ! C decays (kg[C] m-2 day-1).
       d1C = d1p * Cu  (i,j,kp)
       d2C = d2p * Cm  (i,j,kp)
       d3C = d3p * Cv  (i,j,kp)
       d4C = d4p * Cn  (i,j,kp)
       d5C = d5p * Ca  (i,j,kp)
       d6C = d6p * Cs  (i,j,kp)
       d7C = d7p * Cpa (i,j,kp)
       ! C fluxes between litter pools.
       dCu = puf * flittc (kp) + wlittc (kp) - d1C
       dCm = pmf * flittc (kp)               - d2C
       dCv = pvr * rlittc (kp)               - d3C
       dCn = pnr * rlittc (kp)               - d4C
       ! C fluxes to soil pools.
       Ja  = pau * d1C + pam * d2C + pav * d3C + pan * d4C + &
             pas * d6C + pap * d7C
       Js  = psu * d1C + psv * d3C + psa * d5C
       Jpa = ppa * d5C + pps * d6C
       ! C fluxes between soil pools.
       dCa  = Ja  - d5C
       dCs  = Js  - d6C
       dCpa = Jpa - d7C
       dCle = pleach * d5C
       ! C fluxes to structural litter pools.
       Ju = puf * flittc (kp) + wlittc (kp)
       Jv = pvr * rlittc (kp)
       ! N flux from surface litter to surface metabolic.
       NmfNmw = flittn (kp) + wlittn (kp) - vu * Ju
       ! N flux from root litter to soil metabolic.
       Nnr = rlittn (kp) - vv * Jv
       ! N fluxes between litter pools.
       dNu = vu * Ju - d1p * Nu (i,j,kp)
       dNm = NmfNmw  - d2p * Nm (i,j,kp)
       dNv = vv * Jv - d3p * Nv (i,j,kp)
       dNn = Nnr     - d4p * Nn (i,j,kp)
       ! Constrained forms, if required.
       if ((Cm (i,j,kp) + dCm) > eps) then
        temp = (Nm (i,j,kp) + dNm) / (Cm (i,j,kp) + dCm)
       else
	temp = one
       end if
       if (temp < (one / 25.0)) then
        vm = one / 25.0
        NmfNmw = vm * pmf * flittc (kp)
        dNm    = NmfNmw - d2p * Nm (i,j,kp)
       end if
       if (temp > (one / 10.0)) then
        vm = one / 10.0
        NmfNmw = vm * pmf * flittc (kp)
        dNm    = NmfNmw - d2p * Nm (i,j,kp)
       end if
       if ((Cn (i,j,kp) + dCn) > eps) then
        temp = (Nn (i,j,kp) + dNn) / (Cn (i,j,kp) + dCn)
       else
        temp = one
       end if
       if (temp < (one / 25.0)) then
        vn = one / 25.0
        Nnr    = vn * pnr * rlittc (kp)
        dNn    = Nnr    - d4p * Nn (i,j,kp)
       end if
       if (temp > (one / 10.0)) then
        vn = one / 10.0
        Nnr    = vn * pnr * rlittc (kp)
        dNn    = Nnr    - d4p * Nn (i,j,kp)
       end if
       ! N fluxes between soil pools.
       dNa  = va  * Ja  - d5p * Na  (i,j,kp)
       dNs  = vs  * Js  - d6p * Ns  (i,j,kp)
       dNpa = vpa * Jpa - d7p * Npa (i,j,kp)
       ! Inital total C.
       C0 = flittc (kp) + wlittc (kp) + rlittc (kp) + &
            Cu (i,j,kp) + Cm (i,j,kp) + Cv (i,j,kp) + Cn (i,j,kp) + &
            Ca (i,j,kp) + Cs (i,j,kp) + Cpa (i,j,kp)
       ! Initial total N.
       N0 = flittn (kp) + wlittn (kp) + rlittn (kp) + &
            Nu (i,j,kp) + Nm (i,j,kp) + Nv (i,j,kp) + Nn (i,j,kp) + &
            Na (i,j,kp) + Ns (i,j,kp) + Npa (i,j,kp)
       ! Update C pools.
       ! Surface structural litter.
       Cu (i,j,kp) = Cu (i,j,kp) + dCu
       ! Surface metabolic litter.
       Cm (i,j,kp) = Cm (i,j,kp) + dCm
       ! Soil structural litter.
       Cv (i,j,kp) = Cv (i,j,kp) + dCv
       ! Soil metabolic litter.
       Cn (i,j,kp) = Cn (i,j,kp) + dCn
       ! Active.
       Ca (i,j,kp) = Ca (i,j,kp) + dCa
       ! Slow.
       Cs (i,j,kp) = Cs (i,j,kp) + dCs
       ! Passive.
       Cpa (i,j,kp) = Cpa (i,j,kp) + dCpa
       ! Following does not appear to be used.
       Cleach = dCle
       ! Update N pools.
       Nu (i,j,kp) = Nu (i,j,kp) + dNu
       Nm (i,j,kp) = Nm (i,j,kp) + dNm
       Nv (i,j,kp) = Nv (i,j,kp) + dNv
       Nn (i,j,kp) = Nn (i,j,kp) + dNn
       Na (i,j,kp) = Na (i,j,kp) + dNa
       Ns (i,j,kp) = Ns (i,j,kp) + dNs
       Npa (i,j,kp) = Npa (i,j,kp) + dNpa
       ! Total soil C (kg[C] m-2).
       soilC (i,j,kp) = Cu (i,j,kp) + Cm (i,j,kp) + Cv (i,j,kp) + &
                        Cn (i,j,kp) + Ca (i,j,kp) + Cs (i,j,kp) + &
			Cpa (i,j,kp)
       soilN (i,j,kp) = Nu (i,j,kp) + Nm (i,j,kp) + Nv (i,j,kp) + &
                        Nn (i,j,kp) + Na (i,j,kp) + Ns (i,j,kp) + &
			Npa (i,j,kp)
       ! N mineralization rate (kg[N] m-2 day-1).
       Nmin = N0 - soilN (i,j,kp)
       ! Production rate of plant-available N (kg[N] m-2 day-1).
       U = (one - emf) * (Nmin + Nad + Nfx + Ndepo)
       ! Soil respiration (kg[C] m-2 day-1).
       sresp = C0 - (soilC (i,j,kp) + dCle)
       snmin (i,j,kp) = snmin (i,j,kp) + U
      end do
      if (local) then
       kp = 1
       !write (*,*) 
       write (*,*) kyr,it/4,cfoliage(i,j,1,1)
       !write (*,'(2i5,9f8.4)') kyr,it/4,Cu(i,j,kp),Cm(i,j,kp),&
       !                     Cv(i,j,kp), &
       !                     Cn(i,j,kp),Ca(i,j,kp),Cs(i,j,kp),Cpa(i,j,kp),&
       !		    soilC(i,j,kp),1.0e3*sresp
       write (21,'(2i5,9f8.4)') kyr,it/4,Cu(i,j,kp),Cm(i,j,kp),&
                            Cv(i,j,kp), &
                            Cn(i,j,kp),Ca(i,j,kp),Cs(i,j,kp),Cpa(i,j,kp),&
			    soilC(i,j,kp),1.0e3*sresp
       !write (*,'(2i5,8f8.4)') kyr,it/4,Nu(i,j,kp),Nm(i,j,kp),&
       !                     Nv(i,j,kp), &
       !                     Nn(i,j,kp),Na(i,j,kp),Ns(i,j,kp),Npa(i,j,kp),&
       !		    soilN(i,j,kp)
       !write (*,*) kyr,it/4,1.0e3*snmin(i,j,kp),&
       !    soilN(i,j,kp)/soilC(i,j,kp),soilC(i,j,kp)/soilN(i,j,kp)
      end if
     end if ! (mod (it, 4) == zero)
    end do ! it = 1, ntimes
!    stop
   end if ! tmp (i,j,k) /= fillvalue
  end do ! j
 end do ! i
 !---------------------------------------------------------------------!
end do ! kyr = syr, eyr
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Output state variables to restart file if required.
!----------------------------------------------------------------------!
if (rsf_out) then
 file_name = &
            "/home/adf10/rds/rds-mb425-geogscratch/adf10/RSF/rsf_out.nc"
  write (*, *) 'Writing to ', trim (file_name)
  !--------------------------------------------------------------------!
  ! Create netCDF dataset and enter define mode.
  ! nf90_64bit_offset required because file is large. Found rather
  ! randomly. Allows it to work when added lasa to output.
  !--------------------------------------------------------------------!
  call check (nf90_create (trim (file_name), &
              cmode = nf90_64bit_offset, &
              ncid = ncid))
  !--------------------------------------------------------------------!
  ! Define the dimensions.
  !--------------------------------------------------------------------!
  call check (nf90_def_dim (ncid, "longitude", nlon    , lon_dimid ))
  call check (nf90_def_dim (ncid, "latitude" , nlat    , lat_dimid ))
  call check (nf90_def_dim (ncid, "plot"     , nplots  , plot_dimid))
  call check (nf90_def_dim (ncid, "ind"      , nind_max, ind_dimid ))
  !--------------------------------------------------------------------!
  ! Define coordinate variables.
  !--------------------------------------------------------------------!
  call check (nf90_def_var (ncid, "longitude", nf90_float, lon_dimid , &
              lon_varid ))
  call check (nf90_def_var (ncid, "latitude" , nf90_float, lat_dimid , &
              lat_varid ))
  call check (nf90_def_var (ncid, "plot"     , nf90_int  , plot_dimid, &
              plot_varid))
  call check (nf90_def_var (ncid, "ind"      , nf90_int  , ind_dimid , &
              ind_varid ))
  !--------------------------------------------------------------------!
  dimids_three = (/ lon_dimid, lat_dimid, plot_dimid /)
  dimids_four  = (/ lon_dimid, lat_dimid, plot_dimid, ind_dimid /)
  !--------------------------------------------------------------------!
  ! Assign units attributes to coordinate data.
  !--------------------------------------------------------------------!
  call check (nf90_put_att (ncid, lon_varid, "units", "degrees_east"))
  call check (nf90_put_att (ncid, lat_varid, "units", "degrees_north"))
  !--------------------------------------------------------------------!
  ! Define variables.
  !--------------------------------------------------------------------!
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
  call check (nf90_def_var (ncid, "Species_number", &
   nf90_int, dimids_four, varid_kspp))
  call check (nf90_def_var (ncid, "Foliage_C", &
   nf90_float, dimids_four, varid_cfoliage))
  call check (nf90_def_var (ncid, "Foliage_sapwood_area_ratio", &
   nf90_float, dimids_four, varid_lasa))
  !--------------------------------------------------------------------!
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
  call check (nf90_put_att (ncid, varid_kspp  , "units", "n"))
  call check (nf90_put_att (ncid, varid_cfoliage, "units", "kg[C] ind-1"))
  call check (nf90_put_att (ncid, varid_lasa, "units", "m2 m-2"))
  !--------------------------------------------------------------------!
  ! End definitions.
  !--------------------------------------------------------------------!
  call check (nf90_enddef (ncid))
  !--------------------------------------------------------------------!
  ! Write data.
  !--------------------------------------------------------------------!
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
  call check (nf90_put_var (ncid, varid_kspp  , kspp  ))
  call check (nf90_put_var (ncid, varid_cfoliage, cfoliage))
  !call check (nf90_put_var (ncid, varid_lasa, lasa))
  !--------------------------------------------------------------------!
  ! Close file.
  !--------------------------------------------------------------------!
  call check (nf90_close (ncid))
  !--------------------------------------------------------------------!
end if

if ((local) .and. (wrclim)) close (20) ! Local climate output.
if (local) close (21)

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

