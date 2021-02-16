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
!logical :: local  = .TRUE. ! Run only local site?
logical :: local  = .FALSE. ! Run only local site?
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
! Manaus.
!----------------------------------------------------------------------!
real, parameter :: lon_w = -60.035
real, parameter :: lat_w =  -3.053
!----------------------------------------------------------------------!
! Cambridge.
!----------------------------------------------------------------------!
!real, parameter :: lon_w =  0.108
!real, parameter :: lat_w = 52.198
!----------------------------------------------------------------------!
! Mali.
!----------------------------------------------------------------------!
!real, parameter :: lon_w =  0.0
!real, parameter :: lat_w = 21.0
!----------------------------------------------------------------------!
! Indices of local site, if used.
!----------------------------------------------------------------------!
integer :: i_w
integer :: j_w
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Noise on initial soil C across plots?
logical, parameter :: ran_s = .TRUE.
!logical, parameter :: ran_s = .FALSE.
! Start from restart file?
logical, parameter :: rsf_in = .FALSE.
!logical, parameter :: rsf_in = .TRUE.
! Write to restart file?
!logical, parameter :: rsf_out = .FALSE.
logical, parameter :: rsf_out = .TRUE.
integer, parameter :: syr         =  1901 ! Min. 1901
integer, parameter :: eyr         =  2000 ! Max. 2016
integer, parameter :: nland_max   = 67420
integer, parameter :: nlon        =   720
integer, parameter :: nlat        =   360
integer, parameter :: nd          =   365
integer, parameter :: ntimes      =  1460
integer, parameter :: nplots      =     1
integer, parameter :: nind_max    =    12 ! Max. individuals/plot (>1).
integer, parameter :: nspp     =    4
integer, parameter :: mh       =  110
integer, parameter :: sday = 86400 ! s day-1.
real, parameter :: dt = 6.0 * 60. * 60.0 ! s timestep-1
integer, parameter :: fillvalue_int = -99999
real, parameter :: fillvalue = 1.0e20
real, parameter :: tf = 273.15
real, parameter :: eps = 1.0e-8 !****adf
real, parameter :: zero = 0.0
real, parameter :: one = 1.0
! Mol. weight of dry air (kg [air] mol-1).
real, parameter :: m_air = 28.9647
! Mol. weight of water (kg [water] mol-1).
real, parameter :: m_water = 18.01528
real, parameter :: pi = 3.14159265359
real, parameter :: R = 8.31446261815324 ! J K-1 mol-1
real, parameter :: cp = 1012.0
real, parameter :: iscv = 0.2 ! Range of initial soil C (+/-fraction)
real, parameter :: area = 200.0 ! Plot area (m2)
real, parameter :: uwind = 5.0
! Initial slope of photosynthetic light response (mol m-1).
real, parameter :: alphas = 0.04
! Chilling-day upper limit (oC).
real, parameter :: thold = 5.0
! Lower limit of pruba for no frost damage to foliage.
real, parameter :: prubal = 0.13
!                                       C3GR   C4GR   BRCD   NLEV
real, dimension (nspp) :: rhos    = (/  0.20,  0.20,  0.20,  0.11/)
real, dimension (nspp) :: rhop    = (/  0.05,  0.05,  0.05,  0.03/)
real, dimension (nspp) :: ksw     = (/  0.48,  0.48,  0.48,  0.37/)
real, dimension (nspp) :: lsave   = (/  0.00,  0.00,4167.0,3333.0/)
real, dimension (nspp) :: sla     = (/  36.0,  36.0,  36.0,  12.0/)
real, dimension (nspp) :: ao      = (/  0.67,  0.99,  0.67,  0.83/)
real, dimension (nspp) :: fsr     = (/ 0.145, 0.145, 0.145, 0.145/)
real, dimension (nspp) :: frr     = (/  0.86,  0.86,  0.86,  0.86/)
real, dimension (nspp) :: bark    = (/   1.0,   1.0, 0.033,  0.01/)
real, dimension (nspp) :: rlratio = (/   1.0,   1.0,   1.0,   1.0/)
real, dimension (nspp) :: kpar    = (/  0.65,  0.65,  0.65,  0.50/)
real, dimension (nspp) :: d_leaf  = (/  0.04,  0.04,  0.04,  0.04/)
real, dimension (nspp) :: nrc     = (/  9.31,  9.31,  9.31,  9.31/)
real, dimension (nspp) :: ah      = (/  0.00,  0.00, 28.51, 32.95/)
real, dimension (nspp) :: bh      = (/  0.00,  0.00,0.4667,0.5882/)
real, dimension (nspp) :: stf     = (/  0.00,  0.00,  0.22, 0.222/)
real, dimension (nspp) :: formf   = (/  0.00,  0.00,   0.6,  0.56/)
real, dimension (nspp) :: woodd   = (/  0.00,  0.00, 305.0, 205.0/)
real, dimension (nspp) :: live    = (/  0.00,  0.00,  0.17,0.0708/)
real, dimension (nspp) :: storef  = (/  0.00,  0.00,  0.67,  0.67/)
real, dimension (nspp) :: nupc    = (/ 0.036, 0.036, 0.036, 0.036/)
real, dimension (nspp) :: fturn   = (/ 0.031, 0.031,  0.33,  0.33/)
real, dimension (nspp) :: wturn   = (/ 0.031, 0.031,  0.01,  0.01/)
real, dimension (nspp) :: rturn   = (/ 0.031, 0.031,  2.00,  2.00/)
real, dimension (nspp) :: frcoeff = (/   0.5,   0.5,   0.5,   0.5/)
real, dimension (nspp) :: rrcoeff = (/   1.0,   1.0,   1.0,   1.0/)
real, dimension (nspp) :: rgf     = (/  0.75,  0.75,  0.75,  0.75/)
real, dimension (nspp) :: wmf     = (/   0.1,   0.1,   0.1,   0.1/)
integer, dimension (nspp) :: weff  = (/     2,     2,     2,     2/)
integer, dimension (nspp) :: ptype = (/     1,     1,     3,     2/)
! Ratio between gmax and Rubisco N (mol[H2O] m-2 s-1 / (kg[N] m-2)).
real, dimension (nspp) :: ngr     = (/1359.0,4652.0,1672.0,2223.0/)
! Cuticular conductance to CO2 (m s-1).
real, dimension (nspp) :: gmin    = (/4.81e-5,4.81e-5,4.81e-5,4.81e-5/)
real, parameter :: Tc_local = 0.3 ! Clay fraction?
real, parameter :: Ts_local = 0.3 ! Sand fraction?
real, parameter :: plf = 0.01 + 0.04 * Ts_local
real, parameter :: swpmax = -0.033 ! MPa
real, parameter :: bsoil = 5.0
real, parameter :: idbh = 0.001 ! Mean initial tree dbh (m)
!real, parameter :: idbh = 0.1 ! Mean initial tree dbh (m)
real, parameter :: idbhv = 0.9 ! Fraction.
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
real, parameter :: ftmin = 0.0
real, parameter :: ftmax = 40.0
real, parameter :: fkfive = 18.35
real, parameter :: fkeight = (ftmax - fkfive) / (fkfive - ftmin)
real, parameter :: intc   = 0.0005 * dt / float (sday) ! m LAI-1 dt-1
real, parameter :: smeltc = 0.0007 * dt / float (sday) ! m K-1 dt-1.
real, parameter :: fout   = 0.9 ! fraction
character (len = 250) :: file_name ! Generic filename
character (len = 100) :: var_name
character (len =     4) :: char_year
integer :: kyr
integer :: it
integer :: varid,varid1,varid2,varid3
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
integer :: ncid
integer :: lon_dimid, lat_dimid, plot_dimid, ind_p_dimid, ind_t_dimid
integer :: lon_varid, lat_varid, plot_varid, ind_p_varid, ind_t_varid
integer :: land_dimid,land_varid
integer, dimension (1) :: dimids_one
integer, dimension (2) :: dimids_two
integer, dimension (3) :: dimids_three
integer, dimension (3) :: dimids_land
integer :: i,ii,i1,i2
integer :: j,jj,j1,j2
integer :: k,k1,k2
integer :: kp
integer :: ki
integer :: iregl,ireglc
integer :: ksp
integer :: kday,ktemp,jdl
integer :: m
integer :: nland ! For diagnostic.s
integer :: il
integer :: nlayers
integer :: phenf
integer :: cd_flag
integer, dimension (nplots) :: dd_flag
logical :: summer
integer, dimension (nlon,nlat) :: land_index
! Chilling days (days).
integer, dimension (nlon,nlat) :: cd
integer, dimension (nlon,nlat) :: bgs
integer, dimension (nlon,nlat) :: egs
real, dimension (nlat,nd) :: dl
! Degree-days (oC).
real, dimension (nlon,nlat) :: dd
real, dimension (nlon,nlat) :: rltdl
real, dimension (nlon,nlat) :: ddreq
integer, dimension (nland_max,nplots,nind_max) :: k_ind
real, dimension (nlon,nlat) :: LAI_grid
! Grid-box NPP (kg[C] m-1 yr-1).
real, dimension (nlon,nlat) :: NPP_grid
integer, allocatable, dimension (:) :: plot ! No. of plot
real, allocatable, dimension (:,:) :: icwtr_qd ! Ice/water fraction
real, allocatable, dimension (:,:) :: icwtr ! Ice/water fraction
real, allocatable, dimension (:,:,:) :: tmp   ! K
real, allocatable, dimension (:,:,:) :: dswrf ! J m-2
real, allocatable, dimension (:,:,:) :: pre   ! mm 6hr-1
real, allocatable, dimension (:,:,:) :: spfh  ! kg kg-1
real, allocatable, dimension (:,:,:) :: pres  ! Pa
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
real, allocatable, dimension (:,:,:) :: laip
! Soil water holding capacity (m).
real, allocatable, dimension (:,:,:) :: swct
real, allocatable, dimension (:,:,:) :: snow
real, allocatable, dimension (:,:,:) :: soilw1
real, allocatable, dimension (:,:,:) :: soilw2
real, allocatable, dimension (:,:,:) :: soilw3
!----------------------------------------------------------------------!
integer :: nind_total
integer, allocatable, dimension (:) :: ind
integer, allocatable, dimension (:) :: kspp
real, allocatable, dimension (:) :: dbh
real, allocatable, dimension (:) :: hdbh
real, allocatable, dimension (:) :: cfoliage
real, allocatable, dimension (:) :: cfiner
real, allocatable, dimension (:) :: cstore
real, allocatable, dimension (:) :: nfoliage
real, allocatable, dimension (:) :: nbswood
real, allocatable, dimension (:) :: nheart
real, allocatable, dimension (:) :: navail
real, allocatable, dimension (:) :: nfiner
real, allocatable, dimension (:) :: lasa
real, allocatable, dimension (:) :: nitf
integer, allocatable, dimension (:) :: hbc
!----------------------------------------------------------------------!
real, allocatable, dimension (:) :: gmax
real, allocatable, dimension (:) :: cfact
real, allocatable, dimension (:) :: isfact
real, allocatable, dimension (:) :: et_cat
real, allocatable, dimension (:) :: ratiol
real, allocatable, dimension (:) :: ipfact
real, allocatable, dimension (:) :: jmaxfn
integer, allocatable, dimension (:) :: height
real, allocatable, dimension (:) :: kzg
real, allocatable, dimension (:) :: kSWzg
real, allocatable, dimension (:) :: fPARiig
real, allocatable, dimension (:) :: fSWiig
real, allocatable, dimension (:) :: fPAR
real, allocatable, dimension (:) :: fSW
real, allocatable, dimension (:) :: rcfoliage
real, allocatable, dimension (:) :: rcwood
real, allocatable, dimension (:) :: rcfiner
real, allocatable, dimension (:) :: rnfoliage
real, allocatable, dimension (:) :: rnbswood
real, allocatable, dimension (:) :: rnfiner
real, allocatable, dimension (:) :: lsap
real, allocatable, dimension (:) :: cwood
real, allocatable, dimension (:) :: ball
real, allocatable, dimension (:) :: foff
real, allocatable, dimension (:) :: fon
real, allocatable, dimension (:) :: rlold
!----------------------------------------------------------------------!
real, allocatable, dimension (:,:,:) :: swc1
real, allocatable, dimension (:,:,:) :: swc2
real, allocatable, dimension (:,:,:) :: swc3
real, allocatable, dimension (:,:,:) :: SWf
real, allocatable, dimension (:,:,:) :: skzg
real, allocatable, dimension (:,:,:) :: fPARtg
real, allocatable, dimension (:,:,:) :: skSWzg
real, allocatable, dimension (:,:,:) :: fSWtg
real, allocatable, dimension (:) :: wlittc
real, allocatable, dimension (:) :: wlittn
real, allocatable, dimension (:) :: flittc
real, allocatable, dimension (:) :: flittn
real, allocatable, dimension (:) :: rlittc
real, allocatable, dimension (:) :: rlittn
real, allocatable, dimension (:) :: outflow
real, allocatable, dimension (:) :: wfps
real :: t ! Temperature (degree C)
real :: t_d ! Dew-point (degree C)
real :: eo ! Penman evaporation from lake (mm day-1)
real :: pt ! Annual precipitation (m yr-1)
real :: vap ! Vapour pressure (Pa)
real :: e ! Vapour pressure (mbar)
real :: radsw ! W m-2
real :: radpar ! mol[PAR] m-2 s-1
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
real :: temp,C0,N0,Cleach,Nmin,Usoil,Ndepo,nf
real :: sresp ! Soil respiration (kg[C] m-1 day-1)
real :: st ! SW radiation absorbed in upper crown layer (W m-2)
real, dimension (mh+1) :: fPARt
real, dimension (mh+1) :: fSWt
real, dimension (mh+1) :: fPARi
real, dimension (mh+1) :: fSWi
real, dimension (mh+1) :: skz
real, dimension (mh+1) :: skSWz
real, dimension (mh+1) :: eskz
real, dimension (mh+1) :: eskSWz
real, dimension (nind_max) :: farea
real, dimension (nind_max) :: rcstom
real, dimension (nind_max) :: ngains
real, dimension (mh+1,nind_max) :: kz
real, dimension (mh+1,nind_max) :: kSWz
real :: fst,fdt,fdq,ft,fc,fdqftfc
real :: swp1,swp3
real :: swp2 (nplots)
real :: sow1
real :: sc1
real :: sw1,sw2
real, dimension (nspp) :: pruba
real, dimension (nspp) :: slope
real, dimension (nspp) :: rac
real, dimension (nspp) :: rah
real, dimension (nspp) :: rhr
real, dimension (nspp) :: raw
real, dimension (nspp) :: fwsoil
real, dimension (nspp) :: fturn_save
real, dimension (nspp) :: f1,f2,f3
real, dimension (nplots,nspp) :: fturn_plot
real :: prub
real :: nit
real :: gs
real :: caod ! ppm
real, dimension (2019-1700+1) :: cao ! Pa
real :: tairK ! K
real :: tairC ! oC
real :: ppt   ! m dt-1
real :: tt
real :: cwa
real :: cwsta
real :: cwstl
real :: dwo
real :: ET,Ek
real :: cwstf
real :: tfolK
real :: veg_frac
real :: cgc
real :: rwi
real :: rw
real :: gc
real :: raws
real :: rsoil
real :: psin
real :: ef1
real :: lambda
real :: gamma
real :: gammas
real :: s
real :: vpd
real :: betae
real :: fd
real :: kf,kSWf
real :: fracs
real :: psinf
real :: ts
real :: frac
real :: tnfact
real :: pchl
real :: fparb
real :: total
real :: top
real :: rhoa
real :: raf
real :: ccair
real :: tr_grass,tr_tot,tr_rat
real :: rcif,rci,rcn
real :: lmbae
real :: transp
real :: pint,aint
real :: pevap,pevap1,aevap,aevapr
real :: srain
real :: snowi
real :: smelt
real :: transp1,transp2,transp3
real :: out1,out2,outt
real :: ploss
real :: kci,koi,klco,klcc,c,d,jlcmaxn,v,Y
real :: par
real :: efact
real :: vmaxo
real :: tfac
real :: num,den,vmax
real :: kt,cair,rch,rd
real :: ci
real :: ai,ac,ae,A
real :: jmaxn,jn,vcmax,vomax
real :: pcpi
real :: X,Z
real :: alc,bcarb,brubp,ccarb,crubp
real :: dum,q,ccicarb,ccirubp
real :: acarb,arubp
real :: tmind
real :: cffol ! kg[C] ind-1 s-1.
real :: cfwood
real :: cffiner
real :: cflux
real :: tfact,tfaca1,tfaca2,tfacs,ftsoil
real :: diamw
real :: warea,harea,tarea,saparea
real :: wheight,hheight,theight
real :: wwood,hwood,swood
real :: nfp
real :: NPP_global
real :: tup
real :: totalC,totalN
real :: rnfrac
real :: rat
real :: cfmin,maxcf
real :: cso,csn,cto,ctn
real :: cavail,potf,max_folg,potfg
real :: grg,rnmax,rnlose,rnsave
real :: kzog1,kzog2
real :: skzo
real :: kSWzog1,kSWzog2
real :: skSWzo
real :: fSWi_phen
real :: fPARiio,fPARi_ll,fSWiio
real :: lati,loni,rltdli,ddreqi
real :: phi,rn,delta,tn,gn
real :: fcneed,rcneed,wcneed,fnneed,rnneed,wnneed
real :: rat_fol,rat_finc
real :: finc,winc,rinc,tinc
real :: fgr,wgr,rgr
real :: won,woff
real :: hwoodo
real :: rnlayers
real :: ratio
real :: nhinc
real :: winc_min
real :: cfol_max,cfro_max
real :: TOL,X1,X2,wincn,ZBRENT,dw,wood
real :: max_la,sneed,sget
real :: swoodo
real :: ht
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
allocate (laip (nlon,nlat,nplots))
allocate (swct   (nlon,nlat,nplots))
allocate (snow   (nlon,nlat,nplots))
allocate (soilw1 (nlon,nlat,nplots))
allocate (soilw2 (nlon,nlat,nplots))
allocate (soilw3 (nlon,nlat,nplots))
allocate (swc1 (nlon,nlat,nplots))
allocate (swc2 (nlon,nlat,nplots))
allocate (swc3 (nlon,nlat,nplots))
allocate (wlittc (nplots))
allocate (wlittn (nplots))
allocate (flittc (nplots))
allocate (flittn (nplots))
allocate (rlittc (nplots))
allocate (rlittn (nplots))
allocate (outflow(nplots))
allocate (wfps   (nplots))
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
land_index (:,:) = 0
k_ind  (:,:,:) = 0
snow   (:,:,:) = fillvalue
soilw1 (:,:,:) = fillvalue
soilw2 (:,:,:) = fillvalue
soilw3 (:,:,:) = fillvalue
do ksp = 1, nspp
 pruba (ksp) = (one - ao (ksp)) / (one + 12.5 / nrc (ksp))
 slope (ksp) = 71.4 / (one + 12.5 / nrc (ksp))
 rah (ksp) = one / (1.5 * 6.62e-03 * (uwind / d_leaf (ksp)) ** 0.5)
 rhr (ksp) = one / (one / rah (ksp) + one / 213.21)
 rac (ksp) = rah (ksp) / 0.76
 raw (ksp) = rah (ksp) / 1.08
 stf (ksp) = one / (one - stf (ksp))
 f1  (ksp) = stf (ksp) * formf (ksp) * woodd (ksp) * &
             pi * ah (ksp) / 4.0
 f2  (ksp) = one / (bh (ksp) + 2.0)
 f3 (ksp)  = stf (ksp) * formf (ksp) * woodd (ksp)
 fturn_save (ksp) = fturn (ksp) ! Not sure why needed !****adf
end do ! ksp
!----------------------------------------------------------------------!
! Boundary layer resistance factor for gammas.
!----------------------------------------------------------------------!
raf = (rah (1) + raw (1)) / rhr (1)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Read annual global CO2 mixing ratios (ppm -> Pa).
!----------------------------------------------------------------------!
open (10,file='/home/adf10/rds/rds-mb425-geogscratch/adf10/FORCINGS/&
&CO2field/global_co2_ann_1700_2019.txt',status='old')
do kyr = 1700, 2019
 read (10,*) i, caod
 cao   (kyr-1699) = caod / 10.0 ! Pa
end do
close (10)
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
 allocate (ind     (nind_total))
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
 nind_total = nland_max*nplots*nind_max ! Decrease when have regen.
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
      nind (i, j, kp) = 0
      !----------------------------------------------------------------!
      ! Plant grass in plots.
      !----------------------------------------------------------------!
      do ki = 1, 2
       nind (i,j,kp) = nind (i,j,kp) + 1
       ksp = ki
       k = k + 1
       dbh (k) = 0.01 ! LAI of grass
       hbc (k) = 0
       k_ind (land_index(i,j),kp,ki) = k
       kspp (k) = ksp
       cfoliage (k) = dbh (k) * area / sla (ksp)
       !---------------------------------------------------------------!
       ! Grass structure C (kg[C] m-1).
       !---------------------------------------------------------------!
       lsap (k) = bark (ksp) * cfoliage (k)
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
      iregl = 5
      ki = 2
      do ksp = 3, nspp
       do ireglc = 1, iregl
       nind (i,j,kp) = nind (i,j,kp) + 1
       ki = ki + 1
       !do ki = 3, nind (i,j,kp)
       !---------------------------------------------------------------!
       k = k + 1
       k_ind (land_index(i,j),kp,ki) = k
       !---------------------------------------------------------------!
       !ksp = 3
       kspp (k) = 3
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
allocate (SWf    (nlon,nlat,nplots))
allocate (skzg   (nlon,nlat,nplots))
allocate (fPARtg (nlon,nlat,nplots))
allocate (skSWzg (nlon,nlat,nplots))
allocate (fSWtg  (nlon,nlat,nplots))
allocate (gmax   (nind_total))
allocate (cfact  (nind_total))
allocate (isfact (nind_total))
allocate (et_cat (nind_total))
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
!----------------------------------------------------------------------!
do j = j1, j2
 !---------------------------------------------------------------------!
 do i = i1, i2
  if (tmp (i, j, 1) /= fillvalue) then
   do kp = 1, nplots
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
    !------------------------------------------------------------------!
    ! Factor to calculate absorbed SW by whole plot (ratio).
    !------------------------------------------------------------------!
    SWf (i,j,kp) = zero
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
     farea (ki) = cfoliage (k) * sla (ksp)
     rlold (k) = farea (ki)
     !****adf
     !-----------------------------------------------------------------!
     ! Plant height (m).
     !-----------------------------------------------------------------!
     height (k) = 1
     hbc (k) = 0
     !****adf
     ! Number of 1 m layers in crown (n).
     nlayers = height (k) - hbc (k)
     nlayers = max (1, nlayers)
     if (farea (ki) > eps) then
      ! fd is LAI of individual in each layer.
      fd = farea (ki) / (float (nlayers) * area)
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
       SWf (i,j,kp) = SWf (i,j,kp) + (one - rhos (ksp)) * fracs * &
                      fSWi (m)
      end if
     end do ! ki
    end do ! m
    !------------------------------------------------------------------!
    laip (i,j,kp) = zero
    !------------------------------------------------------------------!
    ! Canopy net photosynthesis and conductance.
    !------------------------------------------------------------------!
    do ki = 1, nind (i,j,kp)
     k = k_ind (land_index(i,j),kp,ki)
     ksp = kspp (k)
     !-----------------------------------------------------------------!
     if (farea (ki) > eps) then
      ! Fraction of radiation incident at top of crown absorbed by
      ! individual.
      fPAR (k) = fPAR (k) / fPARt (height (k))
      ! Fraction of radiation incident at top of crown absorbed by
      ! individual.
      fSW (k)  = fSW (k)  / fSWt  (height (k))
      ! Stop FPE overflow from Bob.
      fPAR (k) = max (fPAR (k), eps)
      ! Nitrogen content in upper layer of each crown, scaled with
      ! relative fPAR (kg N m-2) from ratio between radiation absorbed
      ! by top leaf of crown and mean.
      tnfact = farea (ki) * kpar (ksp) / (fPAR (k) * area)
      tnfact = max (one, tnfact)
      ! If foliage in top leaf is too high, reduce foliage N. Assume
      ! extra N just not used for photosynthesis.
      ! IF (tnfact . 2.0) tnfact = 2.0
      nit = tnfact * nfoliage (k) / farea (ki)
      ! nit = min (4.0e-3, nit)
      ! Proportion of N bound in Rubisco in upper layer.
      prub = pruba (ksp) + slope (ksp) * nit
      ! Proportion of N bound in chlorophyll in upper layer.
      pchl = prub / nrc (ksp)
      ! Leaf catalytic site content (mol m-2) is calculated from the
      ! amount of leaf N bound in Rubisco.
      et_cat (k) = prub * nit / 11.0
      ! Factor for calculating jmax.
      jmaxfn (k) = pchl * nit / 0.056
      ! Maximum stomatal conductance in upper layer (m s-1).
      gmax (k) = ngr (ksp) * nit * prub
      ! Fraction of radiation incident at top of individual absorbed
      ! by bottom layer of individual.
      fparb = fPARi (hbc (k) + 1) * kz (hbc (k) + 1, ki) / &
              skz (hbc (k) + 1) / fPARt (height (k))
      ! Factor to calculate net photosynthesis in basal crown layer
      ! from total crown rate.
      ratiol (k) = fparb / fPAR (k)
      ! Factor to calculate canopy net photosynthesis from mean rate
      ! in upper layer.
      cfact (k) = area * fPAR (k) / kpar (ksp)
      ! Factor to calculate absorbed PAR in top layer from plot PAR.
      ipfact (k) = (one - rhop (ksp)) * fPARt (height (k)) * kpar (ksp)
      ! Factor to calculate absorbed SW in top layer from plot SW.
      isfact (k) = (one - rhos (ksp)) * fSWt (height (k)) * ksw (ksp)
      ! Need to set cfact as used for test in energy.f.
     else
      cfact (k) = zero
     end if
     kzg (k) = kz (1,ki)
     fPARiig (k) = fPARi (1) * kz (1,ki) / skz (1)
     kSWzg (k) = kSWz (1,ki)
     fSWiig (k) = fSWi (1) * kz (1,ki) / skz (1)
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
     laip (i,j,kp) = laip (i,j,kp) + farea (ki)
     !-----------------------------------------------------------------!
    end do ! ki
    skzg   (i,j,kp) = skz   (1)
    fPARtg (i,j,kp) = fPARt (1)
    skSWzg (i,j,kp) = skSWz (1)
    fSWtg  (i,j,kp) = fSWt  (1)
    !...Grass lowest layer values.
    !...Total fPAR of lowest layer.
    total = one - eskz (1)
    !...Top 90 % of lowest layer.
    top = one - exp (-0.9 * skz (1))
    !...Grass basal layer ratio.
    k = k_ind (land_index(i,j),kp,1)
    ratiol (k) = (total - top) / total
    k = k_ind (land_index(i,j),kp,2)
    ratiol (k) = (total - top) / total
    laip (i,j,kp) = laip (i,j,kp) / area
   end do ! kp
  end if ! tmp (i,j,1) /= fillvalue
 end do ! i
end do ! j
outflow (:) = zero
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
NPP_grid (:,:) = fillvalue
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Loop through gridboxes and integrate state variables at land points.
! Climate files start at 1901, run through 2019.
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
do kyr = syr, eyr

 !---------------------------------------------------------------------!
 NPP_global = zero
 !---------------------------------------------------------------------!

 !---------------------------------------------------------------------!
 ! C balance of lowest foliage layer (kg[C] ind-1 yr-1).
 !---------------------------------------------------------------------!
 ball (:) = zero
 !---------------------------------------------------------------------!
 
 !---------------------------------------------------------------------!
 ! Flag for cold-deciduous diagnostic.
 !---------------------------------------------------------------------!
 cd_flag = 0
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
    NPP_grid (i,j) = zero
    !------------------------------------------------------------------!
    ! Set GPT-level daily (for grass) or annual (for trees) foliage
    ! turnover rates to default values at start of year.
    !------------------------------------------------------------------!
    do kp = 1, nplots
     !-----------------------------------------------------------------!
     fturn_plot (kp,:) = fturn_save (:)
     !-----------------------------------------------------------------!
    end do ! kp
    !------------------------------------------------------------------!
    do it = 1, ntimes ! Loop timepoints in year.
     !-----------------------------------------------------------------!
     ! Air temperature (K).
     !-----------------------------------------------------------------!
     tairK = tmp (i,j,it)
     !-----------------------------------------------------------------!
     ! Air temperature (oC).
     !-----------------------------------------------------------------!
     tairC = tairK - tf
     !-----------------------------------------------------------------!
     ! Precipitation (m dt-1).
     !-----------------------------------------------------------------!
     ppt = pre (i,j,it) / 1.0e3
     !-----------------------------------------------------------------!
     ! Downwelling solar radiation (W m-2).
     !-----------------------------------------------------------------!
     radsw = dswrf (i,j,it) / dt
     !-----------------------------------------------------------------!
     ! Downwelling PAR (mol[photons] m-2 s-1).
     !-----------------------------------------------------------------!
     radpar = 2.3e-6 * radsw
     !-----------------------------------------------------------------!
     ! Soil temperature, 24-hr minimum temperature (oC).
     ! Slight trick on first day to avoid carry-over requirement
     ! from previous year!****adf.
     !-----------------------------------------------------------------!
     if (it > 3) then
      tsoil = sum (tmp (i,j,it-3:it)) / 4.0 - tf
      tmind = minval (tmp (i,j,it-3:it)) - tf
     else
      tsoil = sum (tmp (i,j,1:it)) / float (it) - tf
      tmind = minval (tmp (i,j,1:it)) - tf
     end if
     !-----------------------------------------------------------------!
     ! Accumulate degree-days (oC).
     !-----------------------------------------------------------------!
     if (tairC > thold) dd (i,j) = dd (i,j) + (tairC - thold) * &
                                   dt / sday
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
      ! Accumulate chilling-days (days).
      !----------------------------------------------------------------!
      if (tmind < thold) cd (i,j) = cd (i,j) + 1
      !----------------------------------------------------------------!
      ! Yearday on hemisphere-basis.
      if (lat (j) >= zero) then
       ktemp = kday
      else
       ktemp = kday - 183
       if (ktemp <= 0) ktemp = ktemp + nd
      end if
      if ((dd (i,j) >= ddreq (i,j)) .and. (bgs (i,j) == (nd + 1))) then
       bgs (i,j) = ktemp + 10
       ! Following ensures DECC act as evergreen if no chance of frosts.
       if (ddreq (i,j) == zero) then
        bgs = ktemp
       else
        ! Set cd_flag because cold-deciduousness must occur.
	! Appears to be only a diagnostic.
        cd_flag = 1
       end if ! ddreq (i,j) == zero
      end if ! dd >= ddreq
      !----------------------------------------------------------------!
      if (ktemp >= bgs (i,j)) then
       summer = .TRUE.
      else
       summer = .FALSE.
      end if
      !----------------------------------------------------------------!
      if (kday > 1) then
       if (summer .and. (dl (j,kday) < dl (j,kday - 1) )) then
        if ((dl (j,kday) <= rltdl (i,j)) .and. (egs (i,j) == (nd + 2)))&
         egs (i,j)  = ktemp
       end if
      end if
      !----------------------------------------------------------------!
      if (ktemp >= egs (i,j)) summer = .FALSE.
      !----------------------------------------------------------------!
      ! Reset bgs and egs to flag values.
      if (ktemp == nd) then
       bgs (i,j) = nd + 1
       egs (i,j) = nd + 2
      end if
      !----------------------------------------------------------------!
     end if ! mod start of day.
     !-----------------------------------------------------------------!
     ! Vapour pressure (Pa).
     !-----------------------------------------------------------------!
     vap = pres (i,j,it) * spfh (i,j,it) * m_air / m_water
     !-----------------------------------------------------------------!
     ! Density of dry air (kg m-3).
     !-----------------------------------------------------------------!
     rhoa = 2.42 - 4.12e-3 * tairK
     !-----------------------------------------------------------------!
     ! Latent heat of vapourisation of H2O (increases by 3% between
     ! 30 oC and 0oC; J kg-1).
     !-----------------------------------------------------------------!
     lambda = 3.15e+6 - 2.375e+3 * tairK
     !-----------------------------------------------------------------!
     ! Psychrometer 'constant' (mol m-3 K-1).
     !-----------------------------------------------------------------!
     gamma = 1.983e7 / (lambda * tairK)
     !-----------------------------------------------------------------!
     ! Adjusted psychrometer `constant' (mol m-3 K-1).
     !-----------------------------------------------------------------!
     gammas = gamma * raf
     !-----------------------------------------------------------------!
     ! Leaf temperature calculation factor (mol J-1).
     !-----------------------------------------------------------------!
     ef1 = gamma / (rhoa * cp)
     !-----------------------------------------------------------------!
     ! Apparent radiative temperature of atmosphere (K)
     !-----------------------------------------------------------------!
     ts = tairK - 0.825 * exp (3.54e-3 * radsw)
     !-----------------------------------------------------------------!
     ! Isothermal net radiation factor.
     !-----------------------------------------------------------------!
     psinf = 5.6703e-8 * (ts ** 4.0 - tairK ** 4.0)
     !-----------------------------------------------------------------!
     ! Atmospheric concentration of CO2 (mol m-3).
     !-----------------------------------------------------------------!
     ccair = cao (kyr-1699) / (R * tairK)
     !-----------------------------------------------------------------!
     ! Saturation concentration of water at air temperature;
     ! after Jones 1992 (mol m-3).
     !-----------------------------------------------------------------!
     cwsta = 613.78 * exp (17.502 * (tairK- tf) / &
            (tairK - 32.18)) / (R * tairK)
     !-----------------------------------------------------------------!
     ! Saturation concentration of water close to air temperature
     ! (mol m-3).
     !-----------------------------------------------------------------!
     cwstl = 613.78 * exp (17.502 * ((tairK + 0.1) - tf) / &
             ((tairK + 0.1) - 32.18)) / (R * (tairK + 0.1))
     !-----------------------------------------------------------------!
     ! Concentration of water vapour in air outside leaf boundary
     ! layer (mol m-3).
     !-----------------------------------------------------------------!
     cwa = vap / (R * tairK)
     cwa = max (zero , cwa)
     cwa = min (cwsta, cwa)
     !-----------------------------------------------------------------!
     ! Water vapour concentration deficit of air outside boundary
     ! layer (mol m-3).
     !-----------------------------------------------------------------!
     dwo = cwsta - cwa
     !-----------------------------------------------------------------!
     ! Slope of saturation vapour pressure curve (mol m-3 K-1).
     !-----------------------------------------------------------------!
     s = 10.0 * (cwstl - cwsta)
     !-----------------------------------------------------------------!
     ! Humidity stomatal conductance response function (fraction).
     ! Full closure at 3000 Pa and passing through 0.2348 at 1500 Pa.
     !-----------------------------------------------------------------!
     fdq = 0.46955 - 0.3815 * dwo
     fdq = max (zero, fdq)
     !-----------------------------------------------------------------!
     ! Temperature stomatal conductance response (fraction).
     !-----------------------------------------------------------------!
     tt = tairC
     tt = max (ftmin, tt)
     tt = min (ftmax, tt)
     if (tt /= ftmax) then
      ft = ((tt - ftmin) * (ftmax - tt) ** fkeight) / &
           ((fkfive - ftmin) * (ftmax - fkfive) ** fkeight)
     else
      ft = zero
     end if
     ft = max (zero, ft)
     ft = min (one , ft)
     !-----------------------------------------------------------------!
     ! Temperature response function.
     !-----------------------------------------------------------------!
     tfact = exp (-6595.0 / tairK)
     !-----------------------------------------------------------------!
     ! Dark respiration temperature effect.
     !-----------------------------------------------------------------!
     if (radpar <= eps) then
      tfaca1 = -42.6e3 * tfact
     else
      tfaca1 = zero
     end if
     !-----------------------------------------------------------------!
     ! Structure respiration temperature effect.
     !-----------------------------------------------------------------!
     tfaca2 = -83.14 * tfact
     !-----------------------------------------------------------------!
     ! Fine root respiration temperature effect.
     !-----------------------------------------------------------------!
     tfacs  = -42.6e3 * exp (-6595.0 / (tsoil + tf))
     !-----------------------------------------------------------------!
     ! Effect of soil temperature on N uptake (m2 kg[C]-1 s-1).
     !-----------------------------------------------------------------!
     ftsoil = (tsoil * (60.0 - tsoil)) / (800.0 * sday)
     ftsoil = max (zero, ftsoil)
     ftsoil = min (one , ftsoil)
     !-----------------------------------------------------------------!
     ! Overall stomatal factor to allow for vpd, temperature, and CO2.
     !-----------------------------------------------------------------!
     fdqftfc = fdq * ft * fc
     !-----------------------------------------------------------------!
     do kp = 1, nplots
      !----------------------------------------------------------------!
      phenf = 0 ! I hope this is correct here?
      !----------------------------------------------------------------!
      ! Convert soil water to water-filled pore space.
      ! Assumes micro-pore space = swc and macro-pore space = 42%
      ! saturation content (from TEM, for loam; Raich et al., 1991).
      !----------------------------------------------------------------!
      if (swct (i,j,kp) > eps) then
       wfps (kp) = 100.0 * (soilw1 (i,j,kp) + soilw2 (i,j,kp) + &
                   soilw3 (i,j,kp)) / (1.7241 * swct (i,j,kp))
       wfps (kp) = min (100.0, wfps (kp))
      else
       wfps (kp) = 0.00001
      end if
      !----------------------------------------------------------------!
      ! Effect of soil water on N uptake through anaerobic conditions at
      ! high soil water and lack of solution medium at low. Function
      ! for decomposition used. weff is GPT parameter for effect of
      ! soil water saturation on N uptake.
      !----------------------------------------------------------------!
      if (snmin (i,j,kp) > eps) then
       if (wfps (kp) < 60.0) then
        fwsoil (:) = exp ((wfps (kp) - 60.0) ** 2 * (-0.00125))
       else
        ! Soil water effects of N uptake by GPT.
        do ksp = 1, nspp
	 if (weff (ksp) == 1) then
          !------------------------------------------------------------!
	  ! Sensitive (class 1) - no uptake at 100% of saturation.
          !------------------------------------------------------------!
          fwsoil (ksp) = 0.000611 * wfps (kp) ** 2 - &
                         0.1222   * wfps (kp) + 6.11
          !------------------------------------------------------------!
          if (wfps (kp) > 100.0) fwsoil (ksp) = zero
          !------------------------------------------------------------!
	 else if (weff (ksp) == 2) then
          !------------------------------------------------------------!
	  ! Average/intermediate (class 2) - same as Nmin. multiplier.
          !------------------------------------------------------------!
          fwsoil (ksp) = 0.000371 * wfps (kp) ** 2 - &
                         0.0748   * wfps (kp) + 4.13
          !------------------------------------------------------------!
          if (wfps (kp) > 100.0) fwsoil (ksp) = 0.36
          !------------------------------------------------------------!
	 else if (weff (ksp) == 3) then
          !------------------------------------------------------------!
	  ! Tolerant (class 3) - 2/3 of peak uptake when soil is
	  ! saturated.
          !------------------------------------------------------------!
          fwsoil (ksp) = 0.0001943 * wfps (kp) ** 2 - &
                         0.0388625 * wfps (kp) + 2.61
          !------------------------------------------------------------!
          if (wfps (kp) > 100.0) fwsoil (ksp) = 0.6667
          !------------------------------------------------------------!
	 end if
         !-------------------------------------------------------------!
	 ! Limit range of fwsoil.
         !-------------------------------------------------------------!
	 fwsoil (ksp) = min (one , fwsoil (ksp))
         fwsoil (ksp) = max (zero, fwsoil (ksp))
         !-------------------------------------------------------------!
	end do ! ksp = 1, nspp
       end if ! wfps < 60
      end if ! snmin > eps
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
        swp2 (kp) = sw1
	tr_rat = one
       else
        ! All transpiration from bottom layer.
        swp2 (kp) = sw2
	tr_rat = zero
       end if
       !---------------------------------------------------------------!
      else
       !---------------------------------------------------------------!
       ! Soil frozen.
       !---------------------------------------------------------------!
       swp1      = -1.48
       swp2 (kp) = -1.48
       tr_rat = zero
       !---------------------------------------------------------------!
      end if ! tsoil > zero
      !----------------------------------------------------------------!
      ! Total plot isothermal net radiation on ground area basis
      ! (W m-2).
      !----------------------------------------------------------------!
      psin = SWf (i,j,kp) * (radsw + psinf)
      !----------------------------------------------------------------!
      ! Total canopy conductance (m[H2O] s-1).
      !----------------------------------------------------------------!
      cgc = zero
      !----------------------------------------------------------------!
      ! Transpiration factors to determine tree vs. grass.
      !----------------------------------------------------------------!
      tr_grass = zero
      tr_tot   = zero
      !----------------------------------------------------------------!
      ! Total N uptake rate in plot (kg[N] s-1).
      !----------------------------------------------------------------!
      tup = zero
      !----------------------------------------------------------------!
      do ki = 1, nind (i,j,kp)
       !---------------------------------------------------------------!
       ! Get index of individual.
       !---------------------------------------------------------------!
       k = k_ind (land_index(i,j),kp,ki)
       ksp = kspp (k)
       !---------------------------------------------------------------!
       if (ksp > 2) then
        harea = pi * (0.5 * hdbh (k)) ** 2
        hheight = ah (ksp) * hdbh (k) ** bh (ksp)
        hwood = stf (ksp) * formf (ksp) * hheight * &
                harea * woodd (ksp)
        totalC = cfoliage (k) + cwood (k) + cstore (k) + &
	         cfiner  (k) - hwood
       else
        totalC = cfoliage (k) + lsap (k) + cfiner (k) + cstore (k)
       end if
       !---------------------------------------------------------------!
       totalN = nfoliage (k) + nbswood (k) + navail (k) + nfiner (k)
       !---------------------------------------------------------------!
       if (totalC > eps) then
        rnfrac = totalN / totalC
       else
        rnfrac = one
       end if
       !---------------------------------------------------------------!
       ! N uptake by individual (kg[N] ind-1 s-1). Denominator
       ! makes it saturating with respect to fine root N. Should add
       ! check to make sure cannot take up more than is there!
       !---------------------------------------------------------------!
       if ((rnfrac < 0.1) .and. (rnfrac > eps)) then
        ngains (ki) = cfiner (k) * nupc (ksp) * snmin (i,j,kp) * &
	              fwsoil (ksp) * ftsoil / rnfrac
        !--------------------------------------------------------------!
        ! For grass, only N available from top layer.
        !--------------------------------------------------------------!
        if (ksp <= 2) then
	 if (swct (i,j,kp) > eps) then
          ngains (ki) = ngains (ki) * &
	                (swc1 (i,j,kp) + swc2 (i,j,kp)) / swct (i,j,kp)
         end if
        end if
       else
        ngains (ki) = zero
       end if
       !---------------------------------------------------------------!
       ! Total N uptake rate in plot (kg[N] s-1).
       !---------------------------------------------------------------!
       tup = tup + ngains (ki)
       !---------------------------------------------------------------!
       ! SW radiation absorbed in upper crown layer (W m-2).
       !---------------------------------------------------------------!
       st = radsw * isfact (k)
       !---------------------------------------------------------------!
       ! Stomatal conductance solar-radiation factor (fraction).
       !---------------------------------------------------------------!
       fst = 1.1044 * st / (st + 104.4)
       !---------------------------------------------------------------!
       ! Stomatal conductance soil-water factor (fraction).
       !---------------------------------------------------------------!
       if ((ki == 1) .or. (ki == 2)) then
        ! Soil water factor for grass.
        if (swp1 > (-1.5)) fdt = 1.155 + 0.77 * &
	 (swp1 - 0.015 * height (k))
         if (swp1 >  (-0.2)) fdt = one
         if (swp1 <= (-1.5)) fdt = zero
       else
        ! Soil water factor for trees.
        if (swp2 (kp) > (-1.5)) fdt = 1.155 + 0.77 * &
	 (swp2 (kp) - (lsave (ksp) / lasa (k)) * 0.005 * &
	 height (k) - 0.01 * height (k))
        if (swp2 (kp) >  (-0.2)) fdt = one
        if (swp2 (kp) <= (-1.5)) fdt = zero
       end if
       !---------------------------------------------------------------!
       ! Stomatal conductance to H2O (mol[H2O] m-2 s-1).
       !---------------------------------------------------------------!
       gs = gmax (k) * fst * fdt * fdqftfc
       gs = max (gmin (ksp), gs)
       !---------------------------------------------------------------!
       ! Stomatal resistance to CO2 (m s-1).
       !---------------------------------------------------------------!
       rcstom (ki) = 1.42 / gs
       !---------------------------------------------------------------!
       ! Total resistance to water transfer across leaf boundary layer
       ! and leaf surface (s m[H2O]-1).
       !---------------------------------------------------------------!
       if (cfact (k) > eps) then
        !--------------------------------------------------------------!
        rwi = (0.607 * rac (ksp) + one / gs) * area / cfact (k)
        !--------------------------------------------------------------!
        ! Sum canopy conductance across individuals (m[H2O s-1).
        !--------------------------------------------------------------!
        cgc = cgc + one / rwi
        !--------------------------------------------------------------!
        ! Sum relative amount of transpiration due to trees.
        !--------------------------------------------------------------!
        tr_tot = tr_tot + one / rwi
        if (ksp <= 2) tr_grass = tr_grass + one / rwi
        !--------------------------------------------------------------!
       end if
      end do ! ki
      !----------------------------------------------------------------!
      if ((dt * tup) > eps) then
       rat = area * snmin (i,j,kp) / (dt * tup)
      else
       rat = zero
      end if
      !----------------------------------------------------------------!
      ! Limit N uptake if more than available in plot.
      !----------------------------------------------------------------!
      if (rat < one) then
       do ki = 1, nind (i,j,kp)
        ngains (ki) = rat * ngains (ki)
       end do ! ki
       tup = rat * tup
      end if
      !----------------------------------------------------------------!
      ! Subtract plant N uptake from soil (kg[N] m-2).
      !----------------------------------------------------------------!
      snmin (i,j,kp) = snmin (i,j,kp) - dt * tup / area
      !----------------------------------------------------------------!
      do ki = 1, nind (i,j,kp)
       !---------------------------------------------------------------!
       k = k_ind (land_index(i,j),kp,ki)
       !---------------------------------------------------------------!
       ! Update individual N store (kg[N] ind-1).
       !---------------------------------------------------------------!
       navail (k) = navail (k) + dt * ngains (ki)
       !---------------------------------------------------------------!
      end do ! ki
      !----------------------------------------------------------------!
      ! tr_rat is converted from flag to the fraction of plot
      ! transpiration from top 2 layers (adjusted if trees taking water
      ! from bottom layer, i.e. tr_rat is set from previous day in this
      ! routine).
      !----------------------------------------------------------------!
      if (tr_rat < 0.1) then
       ! Tree transpiration from bottom (third) layer only.
       if (tr_tot > eps) then
        tr_rat = tr_grass / tr_tot
       else
        tr_rat = 0.5
       end if
      end if
      !----------------------------------------------------------------!
      ! Total resistance to water transfer across leaf boundary layer
      ! and leaf surface on ground area basis (m[H2O] s-1).
      !----------------------------------------------------------------!
      if (cgc > eps) then
       rw = one / cgc
      else
       rw = 0.607 * rac (ksp) + one / gmin (1)
      end if
      !----------------------------------------------------------------!
      ! raws is boundary layer conductance to water of bare ground
      ! calculated from p. 68 of Jones.
      !----------------------------------------------------------------!
      raws = 170.5
      if (snow (i,j,kp) > zero) then
       rsoil = raws + 150.0
      else
       if (swc1 (i,j,kp) > eps) then
        betae = (soilw1 (i,j,kp) / swc1 (i,j,kp) - 0.19) / 0.016
       else
        betae = 0.01
       end if
       if (tsoil <= zero) betae = 0.01
       betae = min (one, betae)
       betae = max (eps, betae)
       rsoil = raws * (one + (one - betae) / betae)
      end if
      rsoil = max (eps, rsoil)
      !----------------------------------------------------------------!
      ! Vegetation canopy fraction of plot (fraction).
      !----------------------------------------------------------------!
      veg_frac = 0.6 + 0.0875 * laip (i,j,kp)
      if (laip (i,j,kp) == zero) veg_frac = zero
      veg_frac = min (0.95, veg_frac)
      !----------------------------------------------------------------!
      ! Update rw to include soil.
      !----------------------------------------------------------------!
      rw = one / (veg_frac / rw + (one - veg_frac) / rsoil)
      !----------------------------------------------------------------!
      ! Mean plot canopy temperature (K) (ef1 is gamma / (rhoa * cp)).
      !----------------------------------------------------------------!
      ksp = 1 ! Use grass for now for boundary conductance.
      !****adf I cannot see how psin is correct for here.
      tfolK = tairK + (rw * psin * ef1 - dwo) / &
              (gamma * (rw + rah (ksp)) / rhr (ksp) + s)
      ! Resistance to CO2 transfer from leaf surface to outside
      ! mesophyll liquid phase factor (83.0e3 * 300.0e-6). 300.0e-6 m
      ! is leaf thickness. 
      rcif = 24.9
      ! Resistance to CO2 transfer from leaf surface to outside
      ! mesophyll liquid phase (s m-1).
      rci = rcif * (tfolK / 293.15) ** (-1.75)
      ! Kinetic parameters in air space equivalents (mol[CO2] m-3).
      kci = 1.9925e15 * exp (- 10127.01  / tfolK) / tfolK
      koi = 2.4239e6  * exp (- 1828.1536 / tfolK) / tfolK
      ! Rubisco oxygenation turnover number (klco) is temperature
      ! dependent.
      klco = 4.397e07 * exp (- 5292.02 / tfolK)
      ! Rubisco carboxylation turnover number (klcc) is temperature
      ! dependent.
      klcc = 2.897e14 * exp (- 9862.41 / tfolK)
      ! Light saturated rate of electron transport
      ! (mol mol(chl)-1 s-1).
      c = 3.486e13 * exp (- 9561.7 / tfolK)
      d = one + exp (78.178 - 23934.4 / tfolK)
      jlcmaxn = c / d
      ! Factors for Ci solution in PGEN.
      v = kci * (one + 8.471 / koi)
      Y = 1.231e-7 * tfolK
      !----------------------------------------------------------------!
      ! efact allows for effect of transpiration on movement of CO2
      ! molecules.
      efact = 2.462e-7 * tfolK
      !----------------------------------------------------------------!
      ! Individual stomatal conductance to CO2 of uppermost leaf.
      !----------------------------------------------------------------!
      do ki = 1, nind (i,j,kp)
       !---------------------------------------------------------------!
       ! Get index of individual.
       !---------------------------------------------------------------!
       k = k_ind (land_index(i,j),kp,ki)
       ksp = kspp (k)
       ! Total resistance to CO2 transfer from outside leaf boundary
       ! layer to outside mesophyll liquid phase (s m-1).
       rcn = rac (ksp) + rcstom (ki) + rci
       rcn = max (eps, rcn)
       ! Total conductance to CO2 transfer from outside leaf boundary
       ! layer to outside mesophyll liquid phase of upper layer
       ! (m s-1), used in PGEN.
       gc = one / rcn
       !---------------------------------------------------------------!
       if ((cfoliage (k) > eps) .and. (et_cat (k) > eps)) then
        if (radpar > zero) then
	 par = radpar * ipfact (k)
	 if (ksp == 2) then ! C4 (uses Collatz et al., 1992).
	  ! Vcmax at 25 oC (mol m-2 s-1).
          vmaxo = 1.248 * et_cat (k)
	  ! Q10 2 response.
          tfac = 2 ** ((tfolK - 298.15) / 10.0)
	  num = vmaxo * tfac
          den = (one + exp (0.3 * (286.15 - tfolK))) * &
                (one + exp (0.3 * (tfolK - 309.15)))
          vmax = num / den
	  kt = 18000.0 * vmax
	  cair = pres (i,j,it) / (R * tfolK)
          rch = one / gc
          rd = 0.021 * vmax
	  num = 2.0 * rch * rd - rch * efact * ccair + 2.0 * ccair
          den = 2.0 + rch * efact + 2.0 * rch * kt / cair
          ci = num / den
	  ai = alphas * par - rd
          ac = kt * ci / cair - rd
          ae = vmax - rd
	  ! Find limiting rate.
          A = ai
          IF (ac < A) A = ac
          IF (ae < A) A = ae
	 else ! C3.
	  ! Light saturated rate of electron transport (mol m-2 s-1).
          jmaxn = jmaxfn (k) * jlcmaxn
	  ! Rate of electron transport (mol m-2 s-1).
          jn = jmaxn * radpar / (par + 2.1 * jmaxn)
	  ! Maximum rate of carboxylation by Rubisco depends on klcc
	  ! and leaf catalytic site content.
          vcmax = klcc * et_cat (k)
	  ! Dark respiration after Collatz et al. (1991)
	  ! (mol CO2 m-2 s-1).
          rd = 0.015 * vcmax
	  ! Maximum rate of oxygenation by Rubisco depends on klco and
	  ! leaf catalytic site content content
	  vomax = klco * et_cat (k)
	  ! Photorespiratory compensation point (mol m-3, air space).
          pcpi = 0.5 * vomax * kci * 8.471 / (vcmax * koi)
	  ! Now analytically solve for the two Cc,i values.
          ! Calculate initial variables.
          X = jn / 4.5
          Z = 2.3333 * pcpi
	  ! alc is the same for all solutions.
          alc = gc + Y
	  bcarb = gc * (v - ccair) + Y * (v + ccair) + vcmax - rd
          brubp = gc * (Z - ccair) + Y * (Z + ccair) + X - rd
          ccarb = v * ccair * (y - gc) - vcmax * pcpi - rd * v
          crubp = Z * ccair * (Y - gc) - pcpi * (X + 2.3333 * rd)
	  ! Stable root solution from Numerical Recipes.
          dum = sqrt (bcarb * bcarb - 4.0 * alc * ccarb)
	  ! 'sign (x, y)' returns abs. of x times sign of y.
          dum = sign (dum, bcarb)
          q = (-0.5) * (bcarb + dum)
          ccicarb = max (q / alc, ccarb / q)
	  dum = sqrt (brubp * brubp - 4.0 * alc * crubp)
          dum = sign (dum, brubp)
          q = (-0.5) * (brubp + dum)
          ccirubp = max (q / alc, crubp / q)
	  ! Carboxyation-limited net photosynthesis.
          acarb = gc * (ccair - ccicarb) - &
	          ((ccair + ccicarb) / 2.0) * efact
          ! RuBP-regeneration-limited net photosynthesis.
          arubp = gc * (ccair - ccirubp) - &
                  ((ccair + ccirubp) / 2.0) * efact
          ! Non co-limited rate of photosynthesis (mol[CO2] m-2 -s1).
          A = min (acarb, arubp)
	  ! Set net photosynthesis to zero if frost.
	  if (tmind <= zero) A = zero
	 end if
	 ! Allow for minimum temperature factor.
         A = nitf (k) * A
	 ! Change units from mol[CO2] to kg[C].
	 A = A * 0.012
	 ! Canopy carbon balance with light (kg[C] ind-1 s-1).
	 cffol = A * cfact (k)
	else
	 ! Canopy carbon balance with no light (kg[C] ind-1 s-1).
	 cffol = tfaca1 * nfoliage (k)
	end if
       else
       	! Canopy carbon balance if no foliage (kg[C] ind-1 s-1).
	cffol = zero
       end if
       !---------------------------------------------------------------!
       ! Structural respiration (kg[C] ind-1 s-1).
       !---------------------------------------------------------------!
       if (ksp <= 2) then
        cfwood = tfaca2 * (lsap (k) + cstore (k))
       else
        cfwood = tfaca2 * lsap (k)
       end if
       !---------------------------------------------------------------!
       ! Fine root respiration (kg[C] ind-1 s-1).
       !---------------------------------------------------------------!
       cffiner = tfacs * nfiner (k)
       !---------------------------------------------------------------!
       ! Individual carbon change (kg[C] ind-1 s-1).
       !---------------------------------------------------------------!
       cflux = cffol + cfwood + cffiner
       !---------------------------------------------------------------!
       ! Update individual carbon store (kg[C] ind-1).
       !---------------------------------------------------------------!
       cstore (k) = cstore (k) + dt * cflux
       !---------------------------------------------------------------!
       ! Annual C balance of lowest foliage layer (kg C ind-2 year-1).
       !---------------------------------------------------------------!
       ball (k) = ball (k) + dt * ratiol (k) * cffol
       !---------------------------------------------------------------!
       ! Grid-box NPP (kg[C] m-2 yr-1).
       !---------------------------------------------------------------!
       NPP_grid (i,j) = NPP_grid (i,j) + cflux
       !---------------------------------------------------------------!
      end do ! ki
      !----------------------------------------------------------------!
      ! Saturation concentration of water at leaf temperature
      ! (mol[H2O] m-3).
      !----------------------------------------------------------------!
      cwstf = 613.78 * exp (17.502 * (tfolK - tf) / &
              (tfolK - 32.18)) / (R * tfolK)
      !----------------------------------------------------------------!
      ! Canopy-to-air water vapour deficit (mol[H2O] m-3).
      !----------------------------------------------------------------!
      vpd = cwstf - cwa
      vpd = max (zero, vpd)
      !----------------------------------------------------------------!
      ! Rate of transpiration from plot (mol[H2O] m-2 s-1).
      !----------------------------------------------------------------!
      ET = vpd / rw
      !----------------------------------------------------------------!
      ! Rate of transpiration from plot (kg[H2O] m-2 dt-1).
      !----------------------------------------------------------------!
      Ek = ET * m_water * dt / 1.0e3
      !----------------------------------------------------------------!
      ! Mean energy used in transpiration (W m-2).
      !----------------------------------------------------------------!
      lmbae = ET * m_water * lambda / 1.0e3
      !----------------------------------------------------------------!
      ! Plot transpiration, assuming 1 g = 1 cm-3 (m[H2O] dt-1).
      transp = Ek / 1000.0
      transp = min (transp, (soilw1 (i,j,kp) + soilw2 (i,j,kp) + &
               soilw3 (i,j,kp)))
      !----------------------------------------------------------------!
      ! Soil water balance.
      if (tairC > eps) then ! Rain.
       ! Potential interception (m dt-1).
       pint = laip (i,j,kp) * intc
       ! Actual interception (m dt-1).
       aint = min (pint, ppt)
       ! Potential intercepted rain (m d-1), calculated from isothermal
       ! net radiation of whole canopy.
       ! Potential evaporation (mol m-2 s-1). Assuming species 1.
       ! Reduce net radiation by that used for evapotranspiration
       ! (not the most correct way, but conserves energy).
       psin = psin - lmbae
       ! mol[H2O] m-2 s-1.
       pevap1 = ((s * psin + rhoa * cp * dwo / rhr (1)) / &
                (s + gammas)) / (m_water * lambda / 1.0e3)
       ! Convert from mol m-2 s-1 to m dt-1.
       pevap = dt * pevap1 * m_water / 1.0e6
       pevap = max (zero, pevap)
       ! Actual evaporated intercepted rain (m dt-1).
       aevap = min (pevap, aint)
       ! Total rain reaching soil (m dt-1)
       srain = ppt - aevap
       ! Evaporation converted to kg[H2O] m-2 dt-1.
       aevapr = aevap * 1.0e3
       ! No snow
       snowi = 0.0
      else ! Snow.
       ! Have snow (m dt-1).
       snowi = ppt
       ! No rain (m dt-1).
       srain = zero
       ! Evaporation (kg[H2O] m-2 dt-1).
       aevapr = zero
      end if
      ! Potential snowmelt (m dt-1).
      smelt = max (zero, (smeltc * tairC))
      ! Actual snowmelt (m dt-1).
      smelt = min (smelt, snow (i,j,kp))
      ! Transpiration from first two soil layers is weighted by soil
      ! water content (assumes root distribution follows soil
      ! water). tr_rat is the fraction of total transpiration from
      ! the first two layers.
      if ((soilw1 (i,j,kp) + soilw2 (i,j,kp)) > eps) then
       transp1 = tr_rat * transp * soilw1 (i,j,kp) / &
                 (soilw1 (i,j,kp) + soilw2 (i,j,kp))
      else
       transp1 = zero
      end if
      transp2 = tr_rat * transp - transp1
      transp3 = transp - (transp1 + transp2)
      ! Update soil water in top layer (m).
      soilw1 (i,j,kp) = soilw1 (i,j,kp) + srain + smelt - transp1
      ! If too much transpiration for top layer to cope with, make
      ! top layer water zero, and subtract required amount from
      ! second layer.
      if (soilw1 (i,j,kp) < zero) then
       transp2 = transp2 - soilw1 (i,j,kp)
       soilw1 (i,j,kp) = zero
      end if
      ! Outflow from first to second layer if top layer has more water
      ! than capacity (m dt-1).
      out1 = max (zero, (soilw1 (i,j,kp) - swc1 (i,j,kp)))
      ! Restrict outflow from top layer to stop second layer having
      ! more than saturation.
      out1 = min (out1, swc2 (i,j,kp) / 0.58 - soilw2 (i,j,kp))
      out1 = max (out1, zero)
      ! Update first (top) layer.
      soilw1 (i,j,kp) = soilw1 (i,j,kp) - out1
      ! Update second (middle) layer.
      soilw2 (i,j,kp) = soilw2 (i,j,kp) + out1 - transp2
      ! If too much transpiration for second layer to cope with,
      ! subtract from third layer.
      if (soilw2 (i,j,kp) < zero) then
       transp3 = transp3 - soilw2 (i,j,kp)
       soilw2 (i,j,kp) = zero
      end if
      ! Outflow from second to third layer if second layer has more
      ! water than capacity (m dt-1).
      out2 = max (zero, (soilw2 (i,j,kp) - swc2 (i,j,kp)))
      ! Restrict outflow from second layer to stop third layer having
      ! more than saturation.
      out2 = min (out2, swc3 (i,j,kp) / 0.58 - soilw3 (i,j,kp))
      out2 = max (out2, zero)
      ! Update second (middle) layer.
      soilw2 (i,j,kp) = soilw2 (i,j,kp) - out2
      ! Update third (bottom) layer.
      soilw3 (i,j,kp) = soilw3 (i,j,kp) + out2 - transp3
      ! Pos. of not conserving water if too much for bottom layer,
      ! but probably very minor problem (should test properly).
      soilw3 (i,j,kp) = max (zero, soilw3 (i,j,kp))
      ! Water above holding capacity of 3rd layer goes to outflow
      ! (m dt-1).
      outt = max (zero, (soilw3 (i,j,kp) - swc3 (i,j,kp)))
      ! Reduce outflow due to drainage impedence.
      outt = outt * fout
      ! Sum outflow for daily use by soil calculations (m day-1).
      outflow (kp) = outflow (kp) + outt
      ! Subtract outflow from soil water in third layer.
      soilw3 (i,j,kp) = soilw3 (i,j,kp) - outt
      !...Proportion of water lost as outflow (fraction)
      if ((soilw1 (i,j,kp) + soilw2 (i,j,kp) + soilw3 (i,j,kp) + &
       outt) > eps) then
       ploss = outt / (soilw1 (i,j,kp) + soilw2 (i,j,kp) + &
               soilw3 (i,j,kp) + outt)
      else
       ploss = zero
      end if
      !-------------------------------------------------------------!
      ! Subtract mineral N in outflow from soil mineral N (kg[N] m-2),
      ! allowing for possibility of negative snmin.
      if (snmin (i,j,kp) > eps) snmin (i,j,kp) = snmin (i,j,kp) * &
       (one - ploss)
      !-------------------------------------------------------------!
      ! Update snowpack (m).
      !-------------------------------------------------------------!
      snow (i,j,kp) = snow (i,j,kp) + snowi - smelt
      !-------------------------------------------------------------!
     end do ! kp
     ! End of day?
     if (mod (it, 4) == zero) then
      !****adf
      nf = ndepi
      !****adf
      !-------------------------------------------------------------!
      ! Mineral N addition from atmosphere (kg[C] m-1 d-1).
      !-------------------------------------------------------------!
      Ndepo = nf + 4.0e-3 * sum (pre (i, j, it-3:it)) / 1000.0
      !-------------------------------------------------------------!
      do kp = 1, nplots
       wlittc (kp) = zero
       wlittn (kp) = zero
       flittc (kp) = zero
       flittn (kp) = zero
       rlittc (kp) = zero
       rlittn (kp) = zero
       laip (i,j,kp) = zero
       do ki = 1, nind (i,j,kp)
        k = k_ind (land_index(i,j),kp,ki)
	ksp = kspp (k)
        !-----------------------------------------------------------!
	if (ksp <= 2) then
         !----------------------------------------------------------!
	 ! Min. amount of grass foliage mass allowed (kg[C] plot-1).
         !----------------------------------------------------------!
	 cfmin = dbh (k) * area / sla (ksp)
         !----------------------------------------------------------!
	 ! Reduce maximum leaf area if ball less than zero.
	 ! Could put in temperature limitation to growth here.
	 ! Also, could make this reduction proportional to LAI.
         !----------------------------------------------------------!
         if (ball (k) < zero) then
          maxcf = 0.9 * cfoliage (k)
         else
          maxcf = cfoliage (k) * 10.0
         end if
         maxcf = max (cfmin, maxcf)
         !----------------------------------------------------------!
	 ! Grass daily C and N  litter (kg day-1).
         !----------------------------------------------------------!
         flitterc = fturn_plot (kp,ksp) * cfoliage (k)
         rlitterc = rturn (ksp) * cfiner (k)
         wlitterc = wturn (ksp) * lsap (k)
         flittern = frcoeff (ksp) * fturn_plot (kp,ksp) * &
                    nfoliage (k)
         rlittern = rrcoeff (ksp) * rturn (ksp) * nfiner (k)
         wlittern = wturn (ksp) * nbswood (k)
         !----------------------------------------------------------!
	 ! Subtract daily grass C and N litter from tissue.
         !----------------------------------------------------------!
         cfoliage (k) = cfoliage (k) - flitterc
         cfiner   (k) = cfiner   (k) - rlitterc
         lsap     (k) = lsap     (k) - wlitterc
         nfoliage (k) = nfoliage (k) - flittern
         nfiner   (k) = nfiner   (k) - rlittern
         nbswood  (k) = nbswood  (k) - wlittern
         !----------------------------------------------------------!
	 ! Increment litter pools (kg plot-1).
         !----------------------------------------------------------!
         wlittc (kp) = wlittc (kp) + wlitterc
         wlittn (kp) = wlittn (kp) + wlittern
         flittc (kp) = flittc (kp) + flitterc
         flittn (kp) = flittn (kp) + flittern
         rlittc (kp) = rlittc (kp) + rlitterc
         rlittn (kp) = rlittn (kp) + rlittern
         !----------------------------------------------------------!
	 ! Initial grass structural C (kg[C] plot-1). Used to calculate
	 ! growth respiration.
         !----------------------------------------------------------!
         cso = cfoliage (k) + lsap (k) + cfiner (k)
         !----------------------------------------------------------!
	 ! Initial total grass C (kg[C] plot-1). Used to calculate C
	 ! balance.
         !----------------------------------------------------------!
         cto = cso + cstore (k)
         !----------------------------------------------------------!
	 ! Total grass C available for partitioning (kg[C] plot-1).
         !----------------------------------------------------------!
         cavail = cfoliage (k) + lsap (k) + cfiner (k) + cstore (k)
         cavail = max (zero, cavail)
         !----------------------------------------------------------!
	 ! Maximum possible grass foliage C based on available C
	 ! (kg[C] plot-1).
         !----------------------------------------------------------!
         potf = cavail / (one + (one + bark (ksp)) * rlratio (ksp))
         !----------------------------------------------------------!
	 ! Grass foliage growth limited by soil water potential.
	 !****adf swp2 should be mean over day.
         !----------------------------------------------------------!
         max_folg = 0.1 * cfoliage (k) * (swp2 (kp) + 1.5) / 1.5
         max_folg = max (zero, max_folg)
         potfg    = cfoliage (k) + max_folg
         !----------------------------------------------------------!
	 ! Actual foliage C lower than potential if ball was -ve.
         !----------------------------------------------------------!
         cfoliage (k) = min (maxcf, potf)
         cfoliage (k) = min (cfoliage (k), potfg)
         !----------------------------------------------------------!
         ! Do not allow foliage C to fall below a minimum (set by
	 ! 'dbh').
         !----------------------------------------------------------!
         cfoliage (k) = max (cfoliage (k), cfmin)
         !----------------------------------------------------------!
	 lsap (k) = bark (ksp) * cfoliage (k)
         cfiner (k) = rlratio (ksp) * cfoliage (k)
         !----------------------------------------------------------!
	 ! New grass structural C (kg[C] plot-1). Used to calculate
	 ! growth respiration and store.
         !----------------------------------------------------------!
         csn = cfoliage (k) + lsap (k) + cfiner (k)
         cstore (k) = cavail - csn
         !----------------------------------------------------------!
	 ! cstore may go negative if cfoliage set to minimum, hence:
         !----------------------------------------------------------!
         cstore (k) = max (cstore (k), zero)
         !----------------------------------------------------------!
	 ! Grass growth respiration, if positive growth.
         !----------------------------------------------------------!
         grg = (one - rgf (ksp)) * max (zero, (csn - cso))
         grg = min (lsap (k), grg)
         !rgs = rgs + grg
         !----------------------------------------------------------!
	 farea (ki) = cfoliage (k) * sla (ksp)
         !----------------------------------------------------------!
	 ! The following allows for growth respiration, but is really
	 ! a cheat as should maintain correct relative proportions.
         !----------------------------------------------------------!
         lsap (k) = lsap (k) - grg
         cwood (k) = lsap (k) + cstore (k)
         !----------------------------------------------------------!
	 ! Grass N available for partitioning (kg[N] plot-1).
         !----------------------------------------------------------!
         navail (k) = navail (k) + nfoliage (k) + &
                      nfiner (k) + nbswood (k)
         rnmax = 0.1 * csn
         rnlose = max (zero, (navail (k) - rnmax))
         navail (k) = navail (k) - rnlose
         rnsave = rnlose
         !----------------------------------------------------------!
	 ! Grass N allocation.
         !----------------------------------------------------------!
         nfoliage (k) = navail (k) * cfoliage (k) / (cfoliage (k) + &
                        fsr (ksp) * lsap (k) + &
                        frr (ksp) * cfiner (k))
         nbswood (k) = navail (k) * fsr (ksp) * lsap (k) / &
	               (cfoliage (k) + fsr (ksp) * lsap (k) + &
                       frr (ksp) * cfiner (k))
         nfiner (k) = navail (k) * frr (ksp) * cfiner (k) / &
                      (cfoliage (k) + fsr (ksp) * lsap (k) + &
                      frr (ksp) * cfiner (k))
         navail (k) = rnsave
         nitf (k) = one
         !----------------------------------------------------------!
	 ! Set C balance of lowest 10 % of grass canopy to zero at
	 ! the end of each day.
         !----------------------------------------------------------! 
         ball (k) = zero
         !----------------------------------------------------------!
	 ! New grass structural C (kg[C] plot=1). Used to calculate
	 ! growth respiration.
         !----------------------------------------------------------!
         ctn = cfoliage (k) + lsap (k) + cfiner (k) + cstore (k) + &
	       grg
         if (cto /= ctn) then
          Cpa (i,j,kp) = Cpa (i,j,kp) - (ctn - cto) / area
          Cpa (i,j,kp) = max (Cpa (i,j,kp), eps)
         end if
         !----------------------------------------------------------!
	else
         !----------------------------------------------------------!
	 ! Tree GPT phenology.
         !----------------------------------------------------------!
	 !****adf
         !wlitterc = area * 0.1 * 5.0 / 365.0
         !wlittern = 2.0 * 0.01 * wlitterc
	 !phenf = 0
	 !****adf
	 ! Calculate and subtract daily litter.
         flitterc = cfoliage (k) * fturn_plot (kp,ksp) / float (nd)
         cfoliage (k) = cfoliage (k) - flitterc
         wlitterc = cwood (k) * wturn (ksp) / float (nd)
         cwood (k) = cwood (k) - wlitterc
         rlitterc = cfiner (k) * rturn (ksp) / float (nd)
         cfiner (k) = cfiner (k) - rlitterc
         !----------------------------------------------------------!
	 flittern = frcoeff (ksp) * nfoliage (k) * &
                    fturn_plot (kp,ksp) / float (nd)
         nfoliage (k) = nfoliage (k) - flittern
         wlittern = nbswood (k) * wturn (ksp) / float (nd)
         nbswood (k) = nbswood (k) - wlittern
         rlittern = rrcoeff (ksp) * nfiner (k) * rturn (ksp) / &
                    float (nd)
         nfiner (k) = nfiner (k) - rlittern
         wlittern = wlittern + navail (k) * wturn (ksp) / float (nd)
         navail (k) = navail (k) - navail (k) * wturn (ksp) / &
	              float (nd)
         !----------------------------------------------------------!
         flittc (kp) = flittc (kp) + flitterc
         flittn (kp) = flittn (kp) + flittern
         wlittc (kp) = wlittc (kp) + wlitterc
         wlittn (kp) = wlittn (kp) + wlittern
         rlittc (kp) = rlittc (kp) + rlitterc
         rlittn (kp) = rlittn (kp) + rlittern
         !----------------------------------------------------------!
	 ! Replenish compartments from stores if possible. Growth
         ! respiration must be allowed for.
         !----------------------------------------------------------!
         fcneed = (rcfoliage (k) - cfoliage (k)) / rgf (ksp)
         rcneed = (rcfiner   (k) - cfiner   (k)) / rgf (ksp)
         wcneed = (rcwood    (k) - cwood    (k)) / rgf (ksp)
         fnneed = rnfoliage (k) - nfoliage (k)
         rnneed = rnfiner   (k) - nfiner   (k)
         wnneed = rnbswood  (k) - nbswood  (k)
         !----------------------------------------------------------!
	 fcneed = max (zero, fcneed)
         rcneed = max (zero, rcneed)
         wcneed = max (zero, wcneed)
         fnneed = max (zero, fnneed)
         rnneed = max (zero, rnneed)
         wnneed = max (zero, wnneed)
         !----------------------------------------------------------!
	 if ((rcfoliage (k) + rcfiner (k)) > eps) then
          rat_fol = rcfoliage (k) / (rcfoliage (k) + rcfiner (k))
         else
          rat_fol = 0.5
         end if
         !----------------------------------------------------------!
	 finc = zero
         winc = zero
         rinc = zero
         !----------------------------------------------------------!
	 if (cstore (k) > eps) then
          if ((fcneed + rcneed) > cstore (k)) then
           finc = rat_fol * cstore (k)
           rinc = (one - rat_fol) * cstore (k)
           winc = zero
           cstore (k) = eps
          else
           finc = fcneed
           rinc = rcneed
           cstore (k) = cstore (k) - (finc + rinc)
           if (wcneed > cstore (k)) then
            winc = cstore (kp)
            cstore (k) = eps
           else
            winc = wcneed
            cstore (k) = cstore (k) - wcneed
           end if
          end if
	 endif
         !----------------------------------------------------------!
	 ! Calculate growth respiration of each compartment (kg C).
         !----------------------------------------------------------!
         fgr = finc * (one - rgf (ksp))
         wgr = winc * (one - rgf (ksp))
         rgr = rinc * (one - rgf (ksp))
         !----------------------------------------------------------!
	 ! Add increments to each compartment, and subtract growth
         ! respiration (kg C).
         !----------------------------------------------------------!
         cfoliage (k) = cfoliage (k) + finc - fgr
         cwood    (k) = cwood    (k) + winc - wgr
         cfiner   (k) = cfiner   (k) + rinc - rgr
         !----------------------------------------------------------!
	 ! Sum site growth respiration (kg C).
         !----------------------------------------------------------!
         !rgs = rgs + fgr + wgr + rgr
         !----------------------------------------------------------!
	 if ((fnneed + rnneed) > navail (k)) then
          nfoliage (k) = nfoliage (k) + rat_fol * navail (k)
          nfiner   (k) = nfiner   (k) + (one - rat_fol) * navail (k)
          navail   (k) = zero
         else
          nfoliage (k) = nfoliage (k) + fnneed
          nfiner   (k) = nfiner   (k) + rnneed
          navail   (k) = navail (k) - (fnneed + rnneed)
          if (wnneed > navail (k)) then
           nbswood (k) = nbswood (k) + navail (k)
           navail (k)  = zero
          else
           nbswood (k) = nbswood (k) + wnneed
           navail  (k) = navail  (k) - wnneed
          end if
         end if
         !----------------------------------------------------------!
	 ! Ensure only lose and/or gain leaves once per year. At
	 ! beginning of year assume dd tree is broadleaf evergreen.
	 ! Will have to improve this for southern hemisphere sites.
         !----------------------------------------------------------!
         if (kday == 1) then
          foff (k) = 0
          fon  (k) = 0
         end if
         !----------------------------------------------------------!
	 ! Evergreen tree, so set leaf area.
         !----------------------------------------------------------!
         if (ptype (ksp) == 2) then
          farea (ki) = rcfoliage (k) * sla (ksp)
         end if
	 ! If this species is cold deciduous, then test for leaf
	 ! area change today.
         if (ptype (ksp) == 3) then
          if (.NOT. summer) then
           farea (ki) = zero
          else
           farea (ki) = rcfoliage (k) * sla (ksp)
          end if
	  ! If leaves fell off then recalculate fturn for cd
	  ! (leaves may not go on and off in ptype 3 if the signals
	  ! and/or phenology parameter values are not appropriate).
          ! rlold is the actual leaf area on the previous day.
	  ! farea is not the actual between allocate and phen; it
	  ! is the actual between phen and allocate.
          if ((farea (ki) == zero) .AND. (rlold (k) > eps) .AND. &
	   (foff (k) == 0)) then
           ! Set leaf-off flag for rest of year.
           foff (k) = 1
           ! Re-calculated fturn to make all foliage enter litter
	   ! pools by the end of the year.
           fturn_plot (kp,ksp) = (float (nd) * (one - &
            (fturn_save (ksp) / float (nd)) * &
            float (kday - 1))) / float (nd + 1 - kday)
          end if
          ! If leaves came on then set to no cold or drought damage.
          if ((farea (ki) > eps) .AND. (rlold (k) == zero)) then
           ! Remove any cold damage.
           nitf (k) = one
           ! Remove any drought damage.
           lasa (k) = lsave (ksp)
          end if
	 ! End of ptype 3 (CD) conditional.
         end if
         !----------------------------------------------------------!
	end if ! grass/tree
	! If this species is potentially dry deciduous then test for
	! leaf area change. Made to lose leaves (once per year) if
	! soil gets below a critical water potential (woff; MPa).
	! Back on when above another value (won; MPa).
        if (ptype (ksp) == 4) then
         woff = -1.49
         won  = -0.5
         ! Set leaf area to previous day's as it may be the first
	 ! day of the year and would have been reset in allocate.
         farea (ki) = rlold (k)
         ! Take leaves off if swp falls low enough. Only do once,
	 ! but could make this happen more often if works OK in
	 ! ALLOCATE with fturn>1.
         if (farea (ki) > eps) then
          if (swp2 (kp) < woff) then
           if (foff (k) == 0) then
            farea (ki) = zero
            foff (k) = 1
            ! Plot-level flag for dry deciduousness.
            dd_flag (kp) = 1
            ! Re-calculated fturn to make all foliage enter litter
	    ! pools by the end of the year.
            fturn_plot (kp,ksp) = (float (nd) * &
	     (one - (fturn_save (ksp) / float (nd)) * &
             float (kday - 1))) / float (nd + 1 - kday)
           end if
          end if
         else
          if (swp2 (kp) > won) then
           if (fon (k) == 0) then
            farea (ki) = rcfoliage (k) * sla (ksp)
            fon (k) = 1
            ! Remove any cold damage.
            nitf (k) = one
            ! Remove any drought damage.
            lasa (k) = lsave (ksp)
            !write (*, *) 'Leaf on because soil wet'
           end if
          end if
         end if
        ! End of ptype 4 (DD) conditional.
        end if
	! If the leaf area changed then set tree leaf area change
	! flag.
        if (farea (ki) /= rlold (k)) then
         ! Flag for tree leaf area change.
         phenf = 1
        end if ! End of leaf area change conditional.
	! See if leaves are present for drought and cold damage.
        if (farea (ki) > (zero + eps)) then
         ! Embolism in tree if soil dry.
         if (swp2 (kp) <= -1.5) then
          lasa (k) = max (one, lasa (k) * 0.9)
         end if
         ! Cold damage if pruba (intercept) high and cold night.
         if (pruba (ksp) > prubal) then
          if (tmind < (-2.3)) then
           nitf (k) = max (0.0001, 0.5 * nitf (k))
          end if
         end if
        ! End of leaf area conditional.
        end if
	! Re-filling of xylem if soil wet enough.
        if (swp2 (kp) >= (-0.5)) then
         lasa (k) = min (lsave (ksp), lasa (k) * 1.1)
        end if
        !-----------------------------------------------------------!
	!laip (i,j,kp) = laip (i,j,kp) + farea (ki)
        !-----------------------------------------------------------!
	! Save foliage area.
        rlold (k) = farea (ki)
       end do ! ki = 1, nind (i,j,kp)
       !------------------------------------------------------------!
       laip (i,j,kp) = laip (i,j,kp) / area
       !------------------------------------------------------------!
       ! Calculate radiation and physiological variables.
       ! Section from phenology.f.
       !------------------------------------------------------------!
       if (phenf == 0) then ! Only grass foliage area changed.
        k1 = k_ind (land_index(i,j),kp,1)
	kzog1 = kzg (k1)
        kzg (k1) = kpar (1) * farea (1) / area
        k2 = k_ind (land_index(i,j),kp,2)
        kzog2 = kzg (k2)
        kzg (k2) = kpar (2) * farea (2) / area
        skzo = skzg (i,j,kp)
        skzg (i,j,kp) = skzo - (kzog1 + kzog2) + &
                        (kzg (k1) + kzg (k2))
	! New fPAR of lowest layer.
        fPARi_ll = one - exp (-skzg (i,j,kp))
	kSWzog1 = kSWzg (k1)
        kSWzog2 = kSWzg (k2)
        kSWzg (k1) = ksw (1) * farea (1) / area
        kSWzg (k2) = ksw (2) * farea (2) / area
        skSWzo = skSWzg (i,j,kp)
        skSWzg (i,j,kp) = skSWzo - (kSWzog1 + kSWzog2) + &
                          (kSWzg (k1) + kSWzg (k2))
	! New fSW of lowest layer.
        fSWi_phen = one - exp (-skSWzg (i,j,kp))
	do ki = 1, nind (i,j,kp)
         k = k_ind (land_index(i,j),kp,ki)
	 ksp = kspp (k)
	 ! Individual's fPAR in lowest layer yesterday.
         fPARiio = fPARiig (k)
	 ! Individual's fPAR in lowest layer today.
	 fPARiig (k) = fPARi_ll * kzg (k) / skzg (i,j,kp)
         ! Total fPAR of individual today.
         fPAR (k) = fPAR (k) + fPARtg (i,j,kp) * &
                    (fPARiig (k) - fPARiio)
         ! Individual's fSW in lowest layer yesterday.
         fSWiio  = fSWiig (k)
	 ! Individual's fSW in lowest layer today.
         fSWiig (k) = fSWi_phen * kSWzg (k) / skSWzg (i,j,kp)
	 ! Total fSW of individual today.
         fSW (k) = fSW (k) + fSWtg (i,j,kp) * &
	           (fSWiig (kp) - fSWiio)
	 ! New SWf.
         SWf (i,j,kp) = SWf (i,j,kp) + fSWtg (i,j,kp) * &
	                (one - rhos (ksp)) * (fSWiig (k) - fSWiio)
	end do ! ki
        !--------------------------------------------------------------!
	! Grass N partitioning.
        !--------------------------------------------------------------!
	do ki = 1, 2
         !-------------------------------------------------------------!
         k = k_ind (land_index(i,j),kp,ki)
	 ksp = kspp (k)
         !-------------------------------------------------------------!
	 ! N content in upper layer of grass, scaled with
         ! relative fPAR (kg[N] m-2) from ratio between radiation
	 ! absorbed by top leaf of crown and mean.
         !-------------------------------------------------------------!
         tnfact = farea (ki) * kpar (ksp) / (fPAR (k) * area)
         tnfact = max (one, tnfact)
	 nit = tnfact * nfoliage (k) / farea (ki)
	 ! Proportion of grass N bound in Rubisco in top foliage.
         prub = pruba (ksp) + slope (ksp) * nit
	 ! Proportion of grass N bound in chlorophyll in top foliage.
         pchl = prub / nrc (ksp)
	 ! Grass foliage catalytic site content (mol m-2) is
	 ! calculated from the amount of leaf N bound in Rubisco.
         et_cat (k) = prub * nit / 11.0
	 ! Factor for calculating jmax.
         jmaxfn (k) = pchl * nit / 0.056
	 ! Maximum stomatal conductance in upper layer (m s-1).
         gmax (k) = ngr (ksp) * nit * prub
	 ! Top 90 % of lowest layer.
         top = one - exp (-0.9 * skzg (i,j,kp))
	 ! Grass basal layer ratio.
         ratiol (k) = (fPARi_ll - top) / fPARi_ll
	 ! Factor to calculate grass net photosynthesis from mean
	 ! rate in upper layer.
         cfact (k) = area * fPAR (k) / kpar (ksp)
         ipfact (k) = (one - rhop (ksp)) * fPARtg (i,j,kp) * &
	              kpar (ksp)
	 ! Factor to calculate absorbed SW in top layer from plot SW.
         isfact (k) = (one - rhos (ksp)) * fSWtg (i,j,kp) * ksw (ksp)
         !-------------------------------------------------------------!
	end do ! ki
       else ! phenf = 1, so tree foliage changed.
	!--------------------------------------------------------------!
        phenf = 0 ! Not sure this needed here.
    !------------------------------------------------------------------!
    ! Factor to calculate absorbed SW by whole plot (ratio).
    !------------------------------------------------------------------!
    SWf (i,j,kp) = zero
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
     nlayers = height (k) - hbc (k)
     nlayers = max (1, nlayers)
     if (farea (ki) > eps) then
      ! fd is LAI of individual in each layer.
      fd = farea (ki) / (float (nlayers) * area)
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
       SWf (i,j,kp) = SWf (i,j,kp) + (one - rhos (ksp)) * fracs * &
                      fSWi (m)
      end if
     end do ! ki
    end do ! m
    !------------------------------------------------------------------!
    ! Canopy net photosynthesis and conductance.
    !------------------------------------------------------------------!
    do ki = 1, nind (i,j,kp)
     k = k_ind (land_index(i,j),kp,ki)
     ksp = kspp (k)
     !-----------------------------------------------------------------!
     if (farea (ki) > eps) then
      ! Fraction of radiation incident at top of crown absorbed by
      ! individual.
      fPAR (k) = fPAR (k) / fPARt (height (k))
      ! Fraction of radiation incident at top of crown absorbed by
      ! individual.
      fSW (k)  = fSW (k)  / fSWt  (height (k))
      ! Stop FPE overflow from Bob.
      fPAR (k) = max (fPAR (k), eps)
      ! Nitrogen content in upper layer of each crown, scaled with
      ! relative fPAR (kg N m-2) from ratio between radiation absorbed
      ! by top leaf of crown and mean.
      tnfact = farea (ki) * kpar (ksp) / (fPAR (k) * area)
      tnfact = max (one, tnfact)
      ! If foliage in top leaf is too high, reduce foliage N. Assume
      ! extra N just not used for photosynthesis.
      ! IF (tnfact . 2.0) tnfact = 2.0
      nit = tnfact * nfoliage (k) / farea (ki)
      ! nit = min (4.0e-3, nit)
      ! Proportion of N bound in Rubisco in upper layer.
      prub = pruba (ksp) + slope (ksp) * nit
      ! Proportion of N bound in chlorophyll in upper layer.
      pchl = prub / nrc (ksp)
      ! Leaf catalytic site content (mol m-2) is calculated from the
      ! amount of leaf N bound in Rubisco.
      et_cat (k) = prub * nit / 11.0
      ! Factor for calculating jmax.
      jmaxfn (k) = pchl * nit / 0.056
      ! Maximum stomatal conductance in upper layer (m s-1).
      gmax (k) = ngr (ksp) * nit * prub
      ! Fraction of radiation incident at top of individual absorbed
      ! by bottom layer of individual.
      fparb = fPARi (hbc (k) + 1) * kz (hbc (k) + 1, ki) / &
              skz (hbc (k) + 1) / fPARt (height (k))
      ! Factor to calculate net photosynthesis in basal crown layer
      ! from total crown rate.
      ratiol (k) = fparb / fPAR (k)
      ! Factor to calculate canopy net photosynthesis from mean rate
      ! in upper layer.
      cfact (k) = area * fPAR (k) / kpar (ksp)
      ! Factor to calculate absorbed PAR in top layer from plot PAR.
      ipfact (k) = (one - rhop (ksp)) * fPARt (height (k)) * kpar (ksp)
      ! Factor to calculate absorbed SW in top layer from plot SW.
      isfact (k) = (one - rhos (ksp)) * fSWt (height (k)) * ksw (ksp)
      ! Need to set cfact as used for test in energy.f.
     else
      cfact (k) = zero
     end if
     kzg (k) = kz (1,ki)
     fPARiig (k) = fPARi (1) * kz (1,ki) / skz (1)
     kSWzg (k) = kSWz (1,ki)
     fSWiig (k) = fSWi (1) * kz (1,ki) / skz (1)
     !-----------------------------------------------------------------!
    end do ! ki
    skzg   (i,j,kp) = skz   (1)
    fPARtg (i,j,kp) = fPARt (1)
    skSWzg (i,j,kp) = skSWz (1)
    fSWtg  (i,j,kp) = fSWt  (1)
    !...Grass lowest layer values.
    !...Total fPAR of lowest layer.
    total = one - eskz (1)
    !...Top 90 % of lowest layer.
    top = one - exp (-0.9 * skz (1))
    !...Grass basal layer ratio.
    k = k_ind (land_index(i,j),kp,1)
    ratiol (k) = (total - top) / total
    k = k_ind (land_index(i,j),kp,2)
    ratiol (k) = (total - top) / total
        !--------------------------------------------------------------!
       end if ! phenf = 1
      end do ! kp
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
      do kp = 1, nplots
       wlittc (kp) = wlittc (kp) / area
       wlittn (kp) = wlittn (kp) / area
       flittc (kp) = flittc (kp) / area
       flittn (kp) = flittn (kp) / area
       rlittc (kp) = rlittc (kp) / area
       rlittn (kp) = rlittn (kp) / area
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
       ! Soil-water decomposition modifer (fraction). Equations are 
       ! fitted to Fig. 8 of Williams et al. (1992).
       ! wfps should really be mean over day.
       !---------------------------------------------------------------!
       if (wfps (kp) < 60.0) then
        em_soil = exp ((wfps (kp) - 60.0) ** 2 / (-800.0))
       else
        em_soil = 0.000371 * wfps (kp) ** 2 - 0.0748 * wfps (kp) + 4.13
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
       h2o_30 = outflow (kp)
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
       !---------------------------------------------------------------!
       ! Soil water holding capacity (Eqn. (1) of FW00; m).
       !---------------------------------------------------------------!
       swct (i,j,kp) = 0.213 + 0.00227 * soilC (i,j,kp)
       !---------------------------------------------------------------!
       ! N mineralization rate (kg[N] m-2 day-1).
       Nmin = N0 - soilN (i,j,kp)
       ! Production rate of plant-available N (kg[N] m-2 day-1).
       Usoil = (one - emf) * (Nmin + Nad + Nfx + Ndepo)
       ! Soil respiration (kg[C] m-2 day-1).
       sresp = C0 - (soilC (i,j,kp) + dCle)
       snmin (i,j,kp) = snmin (i,j,kp) + Usoil
      end do ! kp
      if (local) then
       kp = 1
       !write (*,*)
       !write (*,'(2i5,9f8.4)') kyr,it/4,Cu(i,j,kp),Cm(i,j,kp),&
       ! Cv(i,j,kp),Cn(i,j,kp),Ca(i,j,kp),Cs(i,j,kp),Cpa(i,j,kp),&
       !  soilC(i,j,kp),1.0e3*sresp
       write (21,'(2i5,9f8.4)') kyr,it/4,Cu(i,j,kp),Cm(i,j,kp),   &
                                Cv(i,j,kp),Cn(i,j,kp),Ca(i,j,kp), &
				Cs(i,j,kp),Cpa(i,j,kp),&
			        soilC(i,j,kp),1.0e3*sresp
       write (22,'(2i5,3f8.4)') kyr,it/4,&
                                soilw1(i,j,kp)/swc1(i,j,kp),&
                                soilw2(i,j,kp)/swc2(i,j,kp),&
                                soilw3(i,j,kp)/swc3(i,j,kp)
       !write (*,'(2i5,8f8.4)') kyr,it/4,Nu(i,j,kp),Nm(i,j,kp),&
       !                     Nv(i,j,kp), &
      !              Nn(i,j,kp),Na(i,j,kp),Ns(i,j,kp),Npa(i,j,kp),&
       !		    soilN(i,j,kp)
       !write (*,*) kyr,it/4,1.0e3*snmin(i,j,kp),&
       !    soilN(i,j,kp)/soilC(i,j,kp),soilC(i,j,kp)/soilN(i,j,kp)
      end if
      outflow (:) = zero
     end if ! (mod (it, 4) == zero)
    end do ! it = 1, ntimes
    !------------------------------------------------------------------!
    ! Tree C and N allocation.
    !------------------------------------------------------------------!
    do kp = 1, nplots
     ! Calculate laip.
     laip (i,j,kp) = farea (1) + farea (2)
     do ki = 3, nind (i,j,kp)
      k = k_ind (land_index(i,j),kp,ki)
      ksp = kspp (k)
      !--------------------------C ALLOCATION--------------------------!
      ! Increase heartwood area, etc if base of crown C balance is
      ! negative.
      if (ball (k) < zero) then
       ! Save old heartwood mass (kg C) to calculate proportion of
       ! sapwood converted into heartwood.
       harea = pi * (0.5 * hdbh (k)) ** 2
       hheight = ah (ksp) * hdbh (k) ** bh (ksp)
       hwood = stf (ksp) * formf (ksp) * hheight * &
               harea * woodd (ksp)
       hwoodo = hwood
       ! Real of number of canopy layers
       nlayers = height (k) - hbc (k)
       nlayers = max (1, nlayers)
       rnlayers = float (nlayers)
       ! Increase height to base of crown by 1 m.
       hbc (k) = hbc (k) + 1
       ! Increase heartwood area
       diamw = dbh (k) * (one - 2.0 * bark (ksp))
       warea = pi * (0.5 * diamw) ** 2
       saparea = warea - harea
       harea = harea + saparea / rnlayers
       ! Heartwood area cannot be greater than wood area.
       harea = min (harea, warea)
       ! New sapwood area (m2).
       saparea = warea - harea
       ! New heartwood biomass (kg C).
       if (harea > eps) then
        hdbh (k) = 2.0 * sqrt (harea / pi)
        hheight = ah (ksp) * hdbh (k) ** bh (ksp)
        wheight = ah (ksp) * diamw ** bh (ksp)
        hwood = stf (ksp) * formf (ksp) * hheight * harea * woodd (ksp)
        wwood = stf (ksp) * formf (ksp) * wheight * warea * woodd (ksp)
        hwood = min (wwood, hwood)
       else
        hdbh (k) = zero
        hheight  = zero
        hwood    = zero
       end if
       ! Ratio between heartwood mass increment and old sapwood plus
       ! barkwood mass (because do not separate bark and sapwood N).
       if ((cwood (k) + cstore (k) - hwoodo) > eps) then
        ratio = (hwood - hwoodo) / (cwood (k) + cstore (k) - hwoodo)
       else
        ratio = one
       end if
       ! Make sure ratio positive but not above 1.
       ratio = max (zero, ratio)
       ratio = min (one , ratio)
       !...Heartwood nitrogen increment and sapwood reduction
       nhinc = ratio * nbswood (k)
       nheart (k) = nheart (k) + nhinc
       nbswood (k) = nbswood (k) - nhinc
       ! End of (ball (k) < zero) conditional.
      end if
      ! Make sure cstore is not negative. If was, then subtract
      ! C from wood.
      if (cstore (k) <= eps) then
       cwood (k) = cwood (k) + cstore (k)
       cstore (k) = zero
       winc_min = zero
       if (cwood (k) < eps) then
        ! All cwood used in respiration (N.B. C will not balance now!).
        ! Options are to use foliage and fine root C also.
        ! Or, balance cstore on each day. But, problem of not allowing
        ! for respiration components, favouring EVGR at Penn.
        ! Maintain carbon balance by subtracting deficit from soil.
        Cpa (i,j,kp) = Cpa (i,j,kp) + cwood (k) / area
        cwood (k) = zero
        Cpa (i,j,kp) = max (eps, Cpa (i,j,kp))
        ! Make sure tree dies in mortality, but litter contains all of
        ! the net C and N, allowing for the respiration of C as far as
        ! possible. For convenience, all of C and N put into fine roots.
        ! Setting cfoliage and cstore to 0 ensures mortality.
        cfiner (k) = cfoliage (k) + cwood  (k) + &
	             cfiner   (k) + cstore (k)
        cfoliage (k) = zero
        cwood    (k) = zero
        cfiner   (k) = max (zero, cfiner (k))
        cstore   (k) = zero
        ! Jump over remaining allocation for this individual.
        go to 200 !****adf
       !...End of (cwood (k) <= eps) conditional.
       end if
      else
       ! Enough cstore for some cwood growth.
       ! Minimum wood mass increment (kg[C]).
       winc_min = wmf (ksp) * cstore (k)
       ! Update C store (kg[C]).
       cstore (k) = cstore (k) - winc_min
       ! New sapwood area, if have only have winc_min (m2).
       cwood (kp) = cwood (k) + winc_min
      ! End (cstore (k) <= eps) conditional.
      end if
      ! New dbh and warea.
      ! Added lots of calcs. here as not state variables.
      dbh (k) = (cwood (k) / f1 (ksp)) ** f2 (ksp)
      diamw = dbh (k) * (one - 2.0 * bark (ksp))
      warea = pi * ((0.5 * diamw) ** 2)
      ! warea may have shrunk due to respiration, so keep harea and
      ! hwood sensible.
      harea = pi * (0.5 * hdbh (k)) ** 2
      harea = min (warea, harea)
      hheight = ah (ksp) * hdbh (k) ** bh (ksp)
      wheight = ah (ksp) * diamw ** bh (ksp)
      hwood = stf (ksp) * formf (ksp) * hheight * harea * woodd (ksp)
      wwood = stf (ksp) * formf (ksp) * wheight * warea * woodd (ksp)
      hwood = min (wwood, hwood)
      ! New sapwood area (m2).
      saparea = warea - harea
      saparea = max (zero, saparea)
      ! Maximum foliage mass (kg C).
      cfol_max = lsave (ksp) * saparea / sla (ksp)
      ! Maximum fine root mass (kg C).
      cfro_max = rlratio (ksp) * cfol_max
      ! Required foliage increment (kg C). Can be negative.
      finc = cfol_max - cfoliage (k)
      ! Required fine root increment (kg C). Can be negative.
      rinc = cfro_max - cfiner (k)
      ! Ratio of foliage to foliage plus fine root increments.
      if (abs (finc + rinc) > eps) then
       rat_finc = finc / (finc + rinc)
      else
       rat_finc = zero
      end if
      ! Actual total increment (kg C).
      tinc = min ((finc + rinc), cstore (k))
      ! Subtract total increment from store (kg C).
      cstore (k) = cstore (k) - tinc
      ! Actual foliage increment (kg C).
      finc = rat_finc * tinc
      ! Actual fine root increment (kg C).
      rinc = tinc - finc
      ! Add increments to foliage and fine root compartments (kg C).
      cfoliage (k) = cfoliage (k) + finc
      cfiner   (k) = cfiner   (k) + rinc
      ! See if any cstore left for additional allocation.
      !****adf ZBRENT keeps stopping, so increase this to 1 g as an
      ! experiment.
      if (cstore (k) > 1.0e-3) then
       ! Function ZBRENT solves the pipe model for wincn. Should make
       ! it a subroutine, then could return many useful variables
       ! (save calc. then again).
       TOL = 0.0000000001
       X1 = 0.0
       X2 = cstore (k)
       wincn = ZBRENT(X1,X2,TOL,dw,cwood(k),cfoliage(k), &
               cfiner(k),cstore(k),harea, &
               bark(ksp),ah(ksp),bh(ksp),f1(ksp), &
               f2(ksp),f3(ksp),lsave(ksp),sla(ksp),rlratio (ksp))
       ! New total woody C, including store (subtracted from wincn
       ! later) (kg C). ZBRENT may have had to guess at the root,
       ! hence the following check.
       wincn = min (cstore (k), wincn)
       wincn = max (zero, wincn)
       wood = cwood (k) + wincn
       ! New dbh (m).
       dbh (k) = (wood / f1 (ksp)) ** f2 (ksp)
       ! New inside bark wood diameter (m).
       diamw = dbh (k) * (one - 2.0 * bark (ksp))
       ! New inside bark area (m2).
       warea = pi * ((0.5 * diamw) ** 2)
       ! New inside bark height (m).
       wheight = ah (ksp) * diamw ** bh (ksp)
       ! New inside bark mass (kg C).
       wwood = f3 (ksp) * wheight * warea
       ! New sapwood area (m2).
       saparea = warea - harea
       !...Maximum foliage area from sapwood area (m2).
       max_la   = lsave (ksp) * saparea
       ! Maximum foliage mass (kg C).
       cfol_max = max_la / sla (ksp)
       ! Increment required by foliage (kg C).
       finc = cfol_max - cfoliage (k)
       ! Increment required by fine roots (kg C).
       rinc = rlratio (ksp) * cfol_max - cfiner (k)
       ! New live sapwood mass (kg C).
       lsap (k) = (wwood - hwood) * live (ksp)
       ! Maximum store mass (kg C).
       sneed = storef (ksp) * lsap (k)
       ! Mass available for store, up to maximum (kg C).
       sget = min (sneed, cstore (k) - finc - rinc)
       ! Take store increment out of increment to wood.
       wincn = wincn - sget
       ! Following ensures conservation of C. This is necessary
       ! because the pipe model solution is always an approximation
       ! to TOL.
       wincn = cstore (k) - (finc + rinc + sget)
       ! Add increments to each compartment.
       cfoliage (k) = cfoliage (k) + finc
       cwood (k)    = cwood    (k) + wincn
       cfiner (k)   = cfiner   (k) + rinc
       cstore (k)   = sget
       hwood = min (wwood, hwood)
      else
       ! Cstore is LE epsiln following minimum increment, so need to
       ! adjust harea, etc.
       ! Save old heartwood mass (kg C) for heartwood N calculation.
       hwoodo = hwood
       swoodo = wwood - hwoodo
       ! New inside-bark mass (kg C).
       wheight = ah (ksp) * diamw ** bh (ksp)
       wwood  = f3 (ksp) * wheight * warea
       ! New sapwood area (m2) required to satisfy pipe model.
       saparea = sla (ksp) * cfoliage (k) / lsave (ksp)
       ! New heartwood area (m2).
       harea = warea - saparea
       harea = MAX (zero, harea)
       ! New heartwood biomass (kg C).
       if ( harea > eps) then
        hdbh (k) = 2.0 * sqrt (harea / pi)
        hheight  = ah (ksp) * hdbh (k) ** bh (ksp)
        hwood    = f3 (ksp) * hheight * harea
       else
        hdbh (k) = zero
        hheight  = zero
        hwood    = zero
       end if
       ! Ratio between old heartwood mass increment and old sapwood
       ! plus barkwood mass (because do not separate bark and
       ! sapwood N).
       if (swoodo > eps) then
        ratio = (hwood - hwoodo) / swoodo
       else
        ratio = one
       end if
       ! Make sure ratio OK, hwood may have shrunk.
       ratio = max (-1.0, ratio)
       ratio = min ( 1.0, ratio)
       ! Heartwood nitrogen increment and sapwood reduction.
       nhinc       = ratio * nbswood (k)
       nheart (k)  = nheart (k)  + nhinc
       nbswood (k) = nbswood (k) - nhinc
       ! All of store used for foliage and fine roots.
       cstore (k) = zero
      end if ! cstore (k) > eps
      ! Finished C allocation, now calculate some useful variables.
      ! New sapwood mass (kg C).
      swood = wwood - hwood
      ! New live sapwood mass (kg C).
      lsap  (k) = swood * live (ksp)
      ! Calculate laip.
      laip (i,j,kp) = laip (i,j,kp) + cfoliage (k) * sla (ksp)
      ! Set foliage area to zero (kg C).
      farea (ki) = zero
      ! Height of crown (m) allometrically from dbh (m).
      if (dbh (k) > eps) then
       ht = ah (ksp) * dbh (k) ** bh (ksp)
      else
       ht = zero
      end if
      height (k) = nint (ht + 0.5)
      if (height (k) < 1) height (k) = 1
      if (height (k) > mh) height (k) = mh
      ! Height to base of crown (m).
      hbc (k) = min (hbc (k), (height (k) - 1))
      hbc (k) = max (hbc (k), 0)
      ! Number of 1 m layers in crown. !****adfNeeded?
      nlayers = height (k) - hbc (k)
      nlayers = max (1, nlayers)
      ! Total annual growth respiration (kg C site-1).
      ! Growth respiration of each compartment (kg C).
      finc = cfoliage (k) - rcfoliage (k)
      winc = cwood    (k) - rcwood    (k)
      rinc = cfiner   (k) - rcfiner   (k)
      finc = max (zero, finc)
      winc = max (zero, winc)
      rinc = max (zero, rinc)
      fgr = finc * (one - rgf (ksp))
      wgr = winc * (one - rgf (ksp))
      rgr = rinc * (one - rgf (ksp))
      ! Save required C tree compartment sizes for daily litter
      ! calculations. Done here so that during year can make up to
      ! pipe model solution using store.
      rcfoliage (k) = cfoliage (k)
      rcwood    (k) = cwood    (k)
      rcfiner   (k) = cfiner   (k)
      !Subtract growth respiration from each compartment.
      cfoliage (k) = cfoliage (k) - fgr
      cwood    (k) = cwood    (k) - wgr
      cfiner   (k) = cfiner   (k) - rgr
      ! Sum plot growth respiration (kg C).
      !rgp (kp) = rgp (kp) + fgr + wgr + rgr
      ! Sum plot growth respiration (kg C).
      !rgp (kp) = rgp (kp) + fgr + wgr + rgr
      !-----------------------END OF C ALLOCATION-----------------------
      !---------------------------N ALLOCATION--------------------------
      ! Pool all nitrogen in plant (except heartwood) for allocation
      ! betweenfoliage, sapwood and fine roots (kg).
      navail (k) = navail (k) + nfoliage (k) + nbswood (k) + nfiner (k)
      if (cfoliage (k) > eps) then
       nfoliage (k) = navail (k) * cfoliage (k) / (cfoliage (k) + &
                      fsr (ksp) * swood + frr (ksp) * cfiner (k))
       nbswood (k) = navail (k) * fsr (ksp) * swood / &
                     (cfoliage (k) + fsr (ksp) * swood + &
                     frr (ksp) * cfiner (k))
       nfiner (k) = navail (k) * frr (ksp) * cfiner (k) / &
                    (cfoliage (k) + fsr (ksp) * swood + &
                    frr (ksp) * cfiner (k))
      else
       nfoliage (k) = zero
       nbswood  (k) = navail (k)
       nfiner   (k) = zero
      end if
      navail (k) = zero
      ! Save required N tree compartment sizes for daily litter
      ! calculations.
      rnfoliage (k) = nfoliage (k)
      rnbswood  (k) = nbswood  (k)
      rnfiner   (k) = nfiner   (k)
      !-----------------------END OF N ALLOCATION-----------------------
  200 continue !****adf
     end do ! ki = 3, nind (i,j,kp)
     laip (i,j,kp) = laip (i,j,kp) / area
    end do ! kp
    !------------------------------------------------------------------!
    ! End of tree C and N allocation.
    !------------------------------------------------------------------!
    ! Mortality routines (kill individual trees if low growth.
    !------------------------------------------------------------------!
    do kp = 1, nplots
    end do
    !------------------------------------------------------------------!
    ! End of mortality routines.
    !------------------------------------------------------------------!
    ! m2 m-2.
    LAI_grid (i,j) = sum (laip (i,j,:)) / float (nplots)
    ! kg[C] m-2 yr-1
    NPP_grid (i,j) = (one - icwtr (i,j)) * dt * NPP_grid (i,j) / &
                     (float (nplots) * area)
    ! Pg[C] yr-1.
    !write (99,*)i,j,NPP_grid(i,j)
    NPP_global = NPP_global + larea (i,j) * NPP_grid (i,j) * &
                 1.0e6 / 1.0e12
		 if(isNaN(NPP_global))then
		  write (*,*)i,j,NPP_grid(i,j)
		  write (*,*)i,j,lon(i),lat(j)
		  stop
		 end if
   end if ! tmp (i,j,k) /= fillvalue
  end do ! i = i1, i2
 end do ! j = j1, i2
 !---------------------------------------------------------------------!
 if (local) then
  write (*,*) 'kyr NPP_grid ',kyr,NPP_grid(i1,j1),'kg[C] m-2 yr-1'
  write (*,*) 'kyr LAI_grid ',kyr,LAI_grid(i1,j1),'m2 m-2'
 else
  write (*,*) 'kyr NPP_global ',kyr,NPP_global,'Pg[C] yr-1'
 end if
end do ! kyr = syr, eyr
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
!do j = j1, j2
! do i = i1, i2
!  if (tmp (i,j,1) /= fillvalue) then
!   LAI_grid (i,j) = (one - icwtr (i,j)) * 3.0
!  else
!   LAI_grid (i,j) = fillvalue
!  end if
! end do
!end do
!----------------------------------------------------------------------!
! Output the global LAI field.
!----------------------------------------------------------------------!
file_name = "LAI_grid.nc"
write (*, *) 'Writing to ', trim (file_name)
!----------------------------------------------------------------------!
! Create netCDF dataset and enter define mode.
!----------------------------------------------------------------------!
call check (nf90_create (trim (file_name), cmode = nf90_clobber, &
            ncid = ncid))
!----------------------------------------------------------------------!
! Define the dimensions.
!----------------------------------------------------------------------!
call check (nf90_def_dim (ncid, "longitude", nlon, lon_dimid))
call check (nf90_def_dim (ncid, "latitude" , nlat, lat_dimid))
!----------------------------------------------------------------------!
! Define coordinate variables.
!----------------------------------------------------------------------!
call check (nf90_def_var (ncid, "longitude", nf90_float, lon_dimid, &
            lon_varid))
call check (nf90_def_var (ncid, "latitude" , nf90_float, lat_dimid, &
            lat_varid))
dimids_two = (/ lon_dimid, lat_dimid /)
!----------------------------------------------------------------------!
! Assign units attributes to coordinate data.
!----------------------------------------------------------------------!
call check (nf90_put_att (ncid, lon_varid, "units", "degrees_east"))
call check (nf90_put_att (ncid, lat_varid, "units", "degrees_north"))
!----------------------------------------------------------------------!
! Define variable.
!----------------------------------------------------------------------!
call check (nf90_def_var (ncid, "LAI", nf90_float, &
            dimids_two, varid))
call check (nf90_put_att (ncid, varid, "units", "m2 m-2"))
call check (nf90_put_att (ncid, varid, "_FillValue", fillvalue))
!----------------------------------------------------------------------!
! End definitions.
!----------------------------------------------------------------------!
call check (nf90_enddef (ncid))
!----------------------------------------------------------------------!
! Write data.
!----------------------------------------------------------------------!
call check (nf90_put_var (ncid, lon_varid, lon))
call check (nf90_put_var (ncid, lat_varid, lat))
call check (nf90_put_var (ncid,     varid, LAI_grid))
!----------------------------------------------------------------------!
! Close file.
!----------------------------------------------------------------------!
call check (nf90_close (ncid))
!----------------------------------------------------------------------!
! Output global soil water fields.
!----------------------------------------------------------------------!
file_name = "Soil_water_grid.nc"
write (*, *) 'Writing to ', trim (file_name)
!----------------------------------------------------------------------!
! Create netCDF dataset and enter define mode.
!----------------------------------------------------------------------!
call check (nf90_create (trim (file_name), cmode = nf90_clobber, &
            ncid = ncid))
!----------------------------------------------------------------------!
! Define the dimensions.
!----------------------------------------------------------------------!
call check (nf90_def_dim (ncid, "longitude", nlon, lon_dimid))
call check (nf90_def_dim (ncid, "latitude" , nlat, lat_dimid))
!----------------------------------------------------------------------!
! Define coordinate variables.
!----------------------------------------------------------------------!
call check (nf90_def_var (ncid, "longitude", nf90_float, lon_dimid, &
            lon_varid))
call check (nf90_def_var (ncid, "latitude" , nf90_float, lat_dimid, &
            lat_varid))
dimids_two = (/ lon_dimid, lat_dimid /)
!----------------------------------------------------------------------!
! Assign units attributes to coordinate data.
!----------------------------------------------------------------------!
call check (nf90_put_att (ncid, lon_varid, "units", "degrees_east"))
call check (nf90_put_att (ncid, lat_varid, "units", "degrees_north"))
!----------------------------------------------------------------------!
! Define variables.
!----------------------------------------------------------------------!
call check (nf90_def_var (ncid, "Soil_water1", nf90_float, &
            dimids_two, varid1))
call check (nf90_put_att (ncid, varid1, "units", "m"))
call check (nf90_put_att (ncid, varid1, "_FillValue", fillvalue))
call check (nf90_def_var (ncid, "Soil_water2", nf90_float, &
            dimids_two, varid2))
call check (nf90_put_att (ncid, varid2, "units", "m"))
call check (nf90_put_att (ncid, varid2, "_FillValue", fillvalue))
call check (nf90_def_var (ncid, "Soil_water3", nf90_float, &
            dimids_two, varid3))
call check (nf90_put_att (ncid, varid3, "units", "m"))
call check (nf90_put_att (ncid, varid3, "_FillValue", fillvalue))
!----------------------------------------------------------------------!
! End definitions.
!----------------------------------------------------------------------!
call check (nf90_enddef (ncid))
!----------------------------------------------------------------------!
! Write data.
!----------------------------------------------------------------------!
call check (nf90_put_var (ncid, lon_varid, lon))
call check (nf90_put_var (ncid, lat_varid, lat))
call check (nf90_put_var (ncid,    varid1, soilw1))
call check (nf90_put_var (ncid,    varid2, soilw2))
call check (nf90_put_var (ncid,    varid3, soilw3))
!----------------------------------------------------------------------!
! Close file.
!----------------------------------------------------------------------!
call check (nf90_close (ncid))
!----------------------------------------------------------------------!
! Output the global NPP field.
!----------------------------------------------------------------------!
file_name = "NPP_grid.nc"
write (*, *) 'Writing to ', trim (file_name)
!----------------------------------------------------------------------!
! Create netCDF dataset and enter define mode.
!----------------------------------------------------------------------!
call check (nf90_create (trim (file_name), cmode = nf90_clobber, &
            ncid = ncid))
!----------------------------------------------------------------------!
! Define the dimensions.
!----------------------------------------------------------------------!
call check (nf90_def_dim (ncid, "longitude", nlon, lon_dimid))
call check (nf90_def_dim (ncid, "latitude" , nlat, lat_dimid))
!----------------------------------------------------------------------!
! Define coordinate variables.
!----------------------------------------------------------------------!
call check (nf90_def_var (ncid, "longitude", nf90_float, lon_dimid, &
            lon_varid))
call check (nf90_def_var (ncid, "latitude" , nf90_float, lat_dimid, &
            lat_varid))
dimids_two = (/ lon_dimid, lat_dimid /)
!----------------------------------------------------------------------!
! Assign units attributes to coordinate data.
!----------------------------------------------------------------------!
call check (nf90_put_att (ncid, lon_varid, "units", "degrees_east"))
call check (nf90_put_att (ncid, lat_varid, "units", "degrees_north"))
!----------------------------------------------------------------------!
! Define variable.
!----------------------------------------------------------------------!
call check (nf90_def_var (ncid, "NPP", nf90_float, &
            dimids_two, varid))
call check (nf90_put_att (ncid, varid, "units", "kg[C] m-2 yr-1"))
call check (nf90_put_att (ncid, varid, "_FillValue", fillvalue))
!----------------------------------------------------------------------!
! End definitions.
!----------------------------------------------------------------------!
call check (nf90_enddef (ncid))
!----------------------------------------------------------------------!
! Write data.
!----------------------------------------------------------------------!
call check (nf90_put_var (ncid, lon_varid, lon))
call check (nf90_put_var (ncid, lat_varid, lat))
call check (nf90_put_var (ncid,     varid, NPP_grid))
!----------------------------------------------------------------------!
! Close file.
!----------------------------------------------------------------------!
call check (nf90_close (ncid))
!----------------------------------------------------------------------!
! Output global soil water fields.
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Output state variables to restart file if required.
!----------------------------------------------------------------------!
if (rsf_out) then
 !---------------------------------------------------------------------!
 ! For coordinates in restart file.
 ! 4 byte integer can go up to 2,147,483,647 , I think.
 ! 10 plots over 67420 grid-boxes would allow 3185 individuals/plot.
 !---------------------------------------------------------------------!
 allocate (ind(nind_total))
 do k = 1, nind_total
  ind (k) = k
 end do
 !----------------------------------------------------------------------!

 !----------------------------------------------------------------------!
 file_name = &
            "/home/adf10/rds/rds-mb425-geogscratch/adf10/RSF/rsf_out.nc"
 write (*, *) 'Writing to ', trim (file_name)
 !---------------------------------------------------------------------!
 ! Create netCDF dataset and enter define mode.
 ! nf90_64bit_offset required because file is large. Found rather
 ! randomly. Allows it to work when added lasa to output.
 !---------------------------------------------------------------------!
 call check (nf90_create (trim (file_name), cmode = nf90_64bit_offset, &
             ncid = ncid))
 !---------------------------------------------------------------------!
 ! Define the dimensions.
 !---------------------------------------------------------------------!
 call check (nf90_def_dim (ncid, "longitude", nlon      , lon_dimid  ))
 call check (nf90_def_dim (ncid, "latitude" , nlat      , lat_dimid  ))
 call check (nf90_def_dim (ncid, "plot"     , nplots    , plot_dimid ))
 call check (nf90_def_dim (ncid, "ind_plot" , nind_max  , ind_p_dimid))
 call check (nf90_def_dim (ncid, "ind_total", nind_total, ind_t_dimid))
 call check (nf90_def_dim (ncid, "iland"    , nland_max , land_dimid ))
 !---------------------------------------------------------------------!
 ! Define coordinate variables.
 !---------------------------------------------------------------------!
 call check (nf90_def_var (ncid, "longitude", nf90_float, lon_dimid , &
             lon_varid ))
 call check (nf90_def_var (ncid, "latitude" , nf90_float, lat_dimid , &
             lat_varid ))
 call check (nf90_def_var (ncid, "plot"     , nf90_int  , plot_dimid, &
             plot_varid))
 call check (nf90_def_var (ncid, "ind_plot" , nf90_int  , ind_p_dimid, &
             ind_p_varid))
 call check (nf90_def_var (ncid, "ind_total", nf90_int  , ind_t_dimid, &
             ind_t_varid))
 call check (nf90_def_var (ncid, "iland"    , nf90_int  , land_dimid , &
             land_varid))
 !---------------------------------------------------------------------!
 dimids_one   = (/ ind_t_dimid /)
 dimids_two   = (/ lon_dimid, lat_dimid /)
 dimids_three = (/ lon_dimid, lat_dimid, plot_dimid /)
 dimids_land  = (/ land_dimid, plot_dimid, ind_p_dimid /)
 !---------------------------------------------------------------------!
 ! Assign units attributes to coordinate data.
 !--------------------------------------------------------------------!
 call check (nf90_put_att (ncid, lon_varid, "units", "degrees_east"))
 call check (nf90_put_att (ncid, lat_varid, "units", "degrees_north"))
 !---------------------------------------------------------------------!
 ! Define variables.
 !---------------------------------------------------------------------!
 call check (nf90_def_var (ncid, "Land_index", &
  nf90_int, dimids_two, varid_land_index))
 call check (nf90_def_var (ncid, "Chilling-days", &
  nf90_int, dimids_two, varid_cd))
 call check (nf90_def_var (ncid, "Degree-days", &
  nf90_float, dimids_two, varid_dd))
 call check (nf90_def_var (ncid, "Beginning_of_growing_season", &
  nf90_int, dimids_two, varid_bgs))
 call check (nf90_def_var (ncid, "End_of_growing_season", &
  nf90_int, dimids_two, varid_egs))
 call check (nf90_def_var (ncid, "Individual_index", &
  nf90_int, dimids_land, varid_k_ind))
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
 call check (nf90_def_var (ncid, "Number_individuals_per_plot", &
  nf90_int, dimids_three, varid_nind))
 call check (nf90_def_var (ncid, "Snowpack", &
  nf90_float, dimids_three, varid_snow))
 call check (nf90_def_var (ncid, "Soil_water_layer_1", &
  nf90_float, dimids_three, varid_soilw1))
 call check (nf90_def_var (ncid, "Soil_water_layer_2", &
  nf90_float, dimids_three, varid_soilw2))
 call check (nf90_def_var (ncid, "Soil_water_layer_3", &
  nf90_float, dimids_three, varid_soilw3))
 call check (nf90_def_var (ncid, "Species_number", &
  nf90_int, dimids_one, varid_kspp))
 call check (nf90_def_var (ncid, "DBH", &
  nf90_float, dimids_one, varid_dbh))
 call check (nf90_def_var (ncid, "Heartwood_DBH", &
  nf90_float, dimids_one, varid_hdbh))
 call check (nf90_def_var (ncid, "Foliage_C", &
  nf90_float, dimids_one, varid_cfoliage))
 call check (nf90_def_var (ncid, "Fine_root_C", &
  nf90_float, dimids_one, varid_cfiner))
 call check (nf90_def_var (ncid, "Storage_C", &
  nf90_float, dimids_one, varid_cstore))
 call check (nf90_def_var (ncid, "Foliage_N", &
  nf90_float, dimids_one, varid_nfoliage))
 call check (nf90_def_var (ncid, "Non-heartwood_structure_N", &
  nf90_float, dimids_one, varid_nbswood))
 call check (nf90_def_var (ncid, "Heartwood_N", &
  nf90_float, dimids_one, varid_nheart))
 call check (nf90_def_var (ncid, "Storage_N", &
  nf90_float, dimids_one, varid_navail))
 call check (nf90_def_var (ncid, "Fine_root_N", &
  nf90_float, dimids_one, varid_nfiner))
 call check (nf90_def_var (ncid, "Foliage_sapwood_area_ratio", &
  nf90_float, dimids_one, varid_lasa))
 call check (nf90_def_var (ncid, "Foliage_cold_damage", &
  nf90_float, dimids_one, varid_nitf))
 call check (nf90_def_var (ncid, "Height_to_base_of_crown", &
  nf90_int  , dimids_one, varid_hbc))
 !---------------------------------------------------------------------!
 call check (nf90_put_att (ncid, varid_cd , "units", "days"))
 call check (nf90_put_att (ncid, varid_dd , "units", "oC"))
 call check (nf90_put_att (ncid, varid_bgs, "units", "DOY"))
 call check (nf90_put_att (ncid, varid_egs, "units", "DOY"))
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
 call check (nf90_put_att (ncid, varid_nind , "units", "n"))
 call check (nf90_put_att (ncid, varid_snow , "units", "m"))
 call check (nf90_put_att (ncid, varid_soilw1  , "units", "m"))
 call check (nf90_put_att (ncid, varid_soilw2  , "units", "m"))
 call check (nf90_put_att (ncid, varid_soilw3  , "units", "m"))
 call check (nf90_put_att (ncid, varid_kspp    , "units", "n"))
 call check (nf90_put_att (ncid, varid_dbh     , "units", "m"))
 call check (nf90_put_att (ncid, varid_hdbh    , "units", "m"))
 call check (nf90_put_att (ncid, varid_cfoliage, "units", "kg[C] ind-1"))
 call check (nf90_put_att (ncid, varid_cfiner  , "units", "kg[C] ind-1"))
 call check (nf90_put_att (ncid, varid_cstore  , "units", "kg[C] ind-1"))
 call check (nf90_put_att (ncid, varid_nfoliage, "units", "kg[N] ind-1"))
 call check (nf90_put_att (ncid, varid_nbswood , "units", "kg[N] ind-1"))
 call check (nf90_put_att (ncid, varid_nheart  , "units", "kg[N] ind-1"))
 call check (nf90_put_att (ncid, varid_navail  , "units", "kg[N] ind-1"))
 call check (nf90_put_att (ncid, varid_nfiner  , "units", "kg[N] ind-1"))
 call check (nf90_put_att (ncid, varid_lasa    , "units", "m2 m-2"))
 call check (nf90_put_att (ncid, varid_nitf    , "units", "1-fraction"))
 call check (nf90_put_att (ncid, varid_hbc     , "units", "m"))
 !---------------------------------------------------------------------!
 ! End definitions.
 !---------------------------------------------------------------------!
 call check (nf90_enddef (ncid))
 !---------------------------------------------------------------------!
 ! Write data.
 !---------------------------------------------------------------------!
 call check (nf90_put_var (ncid, lon_varid  ,   lon))
 call check (nf90_put_var (ncid, lat_varid  ,   lat))
 call check (nf90_put_var (ncid, plot_varid ,  plot))
 call check (nf90_put_var (ncid, ind_t_varid,   ind))
 call check (nf90_put_var (ncid, varid_land_index, land_index))
 call check (nf90_put_var (ncid, varid_cd, cd))
 call check (nf90_put_var (ncid, varid_dd, dd))
 call check (nf90_put_var (ncid, varid_bgs, bgs))
 call check (nf90_put_var (ncid, varid_egs, egs))
 call check (nf90_put_var (ncid, varid_k_ind, k_ind))
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
 call check (nf90_put_var (ncid, varid_snow, snow))
 call check (nf90_put_var (ncid, varid_soilw1, soilw1))
 call check (nf90_put_var (ncid, varid_soilw2, soilw2))
 call check (nf90_put_var (ncid, varid_soilw3, soilw3))
 call check (nf90_put_var (ncid, varid_kspp  , kspp  ))
 call check (nf90_put_var (ncid, varid_dbh     , dbh))
 call check (nf90_put_var (ncid, varid_hdbh    , hdbh))
 call check (nf90_put_var (ncid, varid_cfoliage, cfoliage))
 call check (nf90_put_var (ncid, varid_cfiner  , cfiner))
 call check (nf90_put_var (ncid, varid_cstore  , cstore))
 call check (nf90_put_var (ncid, varid_nfoliage, nfoliage))
 call check (nf90_put_var (ncid, varid_nbswood , nbswood))
 call check (nf90_put_var (ncid, varid_nheart  , nheart))
 call check (nf90_put_var (ncid, varid_navail  , navail))
 call check (nf90_put_var (ncid, varid_nfiner  , nfiner))
 call check (nf90_put_var (ncid, varid_lasa    , lasa))
 call check (nf90_put_var (ncid, varid_nitf    , nitf))
 call check (nf90_put_var (ncid, varid_hbc     , hbc))
 !---------------------------------------------------------------------!
 ! Close file.
 !---------------------------------------------------------------------!
 call check (nf90_close (ncid))
 !---------------------------------------------------------------------!
end if

if ((local) .and. (wrclim)) close (20) ! Local climate output.
if (local) close (21) ! soil carbon output.
if (local) close (22) ! soil water output.

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

! -------------------------------------------------------------------- !
      FUNCTION ZBRENT(X1,X2,TOL, &
                     dw, cwood, cfoliage, cfiner, cstore, harea, &
                     bark, ah, bh, f1, f2, f3, lasa, sla, rlratio)
! -------------------------------------------------------------------- !
! This code is taken directly from p. 253 of Press et al. 1989.
! Using Brent's method, find the root of a function Root known to lie
! between X1 and X2. The root, returned as ZBRENT, will be refined
! until its accuracy is TOL.
! -------------------------------------------------------------------- !
      REAL dw, cwood, cfoliage, cfiner, cstore, harea, &
          bark, ah, bh, f3, lasa, sla, rlratio
      REAL wood, dbh, diamw, warea, wheight, wwood, saparea
      REAL max_la, cfol_max, finc, rinc
!      PARAMETER (ITMAX=100,EPS=3.E-8,pi = 3.14159)
!      PARAMETER (ITMAX=10,EPS=0.12E-06,pi = 3.14159)
      PARAMETER (ITMAX=10)
      EPS=0.12E-06
      pi = 3.14159

      A=X1
      B=X2
      
      dw=A
      wood = cwood + dw
      dbh = (wood / f1) ** f2
      diamw = dbh * (1.0 - 2.0 * bark)
      warea = pi * ((0.5 * diamw) ** 2)
      wheight = ah * diamw * bh
      wwood = f3 * wheight * warea
      saparea = max (EPS, (warea - harea))
      max_la = lasa * saparea
      cfol_max = max_la / sla
      finc = cfol_max - cfoliage
      rinc = rlratio * cfol_max - cfiner
      FA = cstore - dw - finc - rinc
      
      dw=B
      wood = cwood + dw
      dbh = (wood / f1) ** f2
      diamw = dbh * (1.0 - 2.0 * bark)
      warea = pi * ((0.5 * diamw) ** 2)
      wheight = ah * diamw * bh
      wwood = f3 * wheight * warea
      saparea = max (EPS, (warea - harea))
      max_la = lasa * saparea
      cfol_max = max_la / sla
      finc = cfol_max - cfoliage
      rinc = rlratio * cfol_max - cfiner
      FB = cstore - dw - finc - rinc
      
      IF(FB*FA.GT.0.) THEN
          WRITE (*, *) 'A, B', A, B
          WRITE (*, *) 'FA, FB', FA, FB
          WRITE (*, *) 'cstore, cwood', cstore, cwood
          WRITE (*, *) 'wood', wood
          WRITE (*, *) dbh
          WRITE (*, *) diamw
          WRITE (*, *) 'warea=',warea
	  write (*,*) 'harea=',harea
          WRITE (*, *) wheight
          WRITE (*, *) wwood
          WRITE (*, *) saparea
          WRITE (*, *) max_la
          WRITE (*, *) cfol_max
          WRITE (*, *) finc
          WRITE (*, *) rinc
          PAUSE 'Root must be bracketed for ZBRENT.'
      ENDIF
      FC=FB
      DO 100 ITER=1,ITMAX
          IF(FB*FC.GT.0.) THEN
              C=A
              FC=FA
              D=B-A
              E=D
          ENDIF
          IF(ABS(FC).LT.ABS(FB)) THEN
              A=B
              B=C
              C=A
              FA=FB
              FB=FC
              FC=FA
          ENDIF
          TOL1=2.*EPS*ABS(B)+0.5*TOL
          XM=.5*(C-B)
          IF(ABS(XM).LE.TOL1 .OR. FB.EQ.0.) THEN
              ZBRENT=B
              RETURN
          ENDIF
          IF(ABS(E).GE.TOL1 .AND. ABS(FA).GT.ABS(FB)) THEN
              S=FB/FA
              IF(A.EQ.C) THEN
                  P=2.*XM*S
                  Q=1.-S
              ELSE
                  Q=FA/FC
                  R=FB/FC
                  P=S*(2.*XM*Q*(Q-R)-(B-A)*(R-1.))
                  Q=(Q-1.)*(R-1.)*(S-1.)
              ENDIF
              IF(P.GT.0.) Q=-Q
              P=ABS(P)
              IF(2.*P .LT. MIN(3.*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
                  E=D
                  D=P/Q
              ELSE
                  D=XM
                  E=D
              ENDIF
          ELSE
              D=XM 
              E=D
          ENDIF
          A=B
          FA=FB
          IF(ABS(D) .GT. TOL1) THEN
              B=B+D
          ELSE
              B=B+SIGN(TOL1,XM)
          ENDIF
      
          dw=B
          wood = cwood + dw
          dbh = (wood / f1) ** f2
          diamw = dbh * (1.0 - 2.0 * bark)
          warea = pi * ((0.5 * diamw) ** 2)
          wheight = ah * diamw * bh
          wwood = f3 * wheight * warea
          saparea = warea - harea
          max_la = lasa * saparea
          cfol_max = max_la / sla
          finc = cfol_max - cfoliage
          rinc = rlratio * cfol_max - cfiner
          FB = cstore - dw - finc - rinc
      
  100 CONTINUE
!....ADF_new      PAUSE 'ZBRENT exceeded maxumum iterations.'
      ZBRENT=B
      RETURN
      END
! -------------------------------------------------------------------- !

