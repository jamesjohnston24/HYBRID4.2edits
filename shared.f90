module shared
implicit none
integer, parameter :: nd        =    365
integer, parameter :: nland_max =  67420
integer, parameter :: nlon      =    720
integer, parameter :: nlat      =    360
integer, parameter :: mh        =    110
integer, parameter :: nspp      =      4
real   , parameter :: zero      =    0.0
real   , parameter :: one       =    1.0
real   , parameter :: eps       = 1.0e-8 !****adf
real   , parameter :: sday      = 86400.0 ! s day-1
real   , parameter :: dt        = 6.0 * 60.0 * 60.0 ! s timestep-1
real   , parameter :: pi        = 3.14159265359
real   , parameter :: R         = 8.31446261815324 ! J K-1 mol-1
real   , parameter :: cp = 1012.0
real   , parameter :: tf        = 273.15
! Mol. weight of dry air (g [air] mol-1).
real, parameter :: m_air   = 28.9647
! Mol. weight of water (g [water] mol-1).
real, parameter :: m_water = 18.01528
! d_one is max. swc of top layer for roots (m).
real, parameter :: d_one = 0.2
logical, parameter :: wrclim = .FALSE. ! Write local climate?
!                                       C3GR   C4GR   BRCD   NLEV
real, dimension (nspp) :: kpar    = (/  0.65,  0.65,  0.65,  0.50/)
real, dimension (nspp) :: kSW     = (/  0.48,  0.48,  0.48,  0.37/)
real, dimension (nspp) :: rhos    = (/  0.20,  0.20,  0.20,  0.11/)
real, dimension (nspp) :: nrc     = (/  9.31,  9.31,  9.31,  9.31/)
real, dimension (nspp) :: rhop    = (/  0.05,  0.05,  0.05,  0.03/)
real, dimension (nspp) :: sla     = (/  36.0,  36.0,  36.0,  12.0/)
real, dimension (nspp) :: bark    = (/   1.0,   1.0, 0.033,  0.01/)
real, dimension (nspp) :: ah      = (/  0.00,  0.00, 28.51, 32.95/)
real, dimension (nspp) :: bh      = (/  0.00,  0.00,0.4667,0.5882/)
real, dimension (nspp) :: stf     = (/  0.00,  0.00,  0.22, 0.222/)
real, dimension (nspp) :: formf   = (/  0.00,  0.00,   0.6,  0.56/)
real, dimension (nspp) :: woodd   = (/  0.00,  0.00, 305.0, 205.0/)
real, dimension (nspp) :: lsave   = (/  0.00,  0.00,4167.0,3333.0/)
real, dimension (nspp) :: rlratio = (/   1.0,   1.0,   1.0,   1.0/)
real, dimension (nspp) :: live    = (/  0.00,  0.00,  0.17,0.0708/)
real, dimension (nspp) :: storef  = (/  0.00,  0.00,  0.67,  0.67/)
real, dimension (nspp) :: fsr     = (/ 0.145, 0.145, 0.145, 0.145/)
real, dimension (nspp) :: frr     = (/  0.86,  0.86,  0.86,  0.86/)
real, dimension (nspp) :: wmf     = (/   0.1,   0.1,   0.1,   0.1/)
real, dimension (nspp) :: rgf     = (/  0.75,  0.75,  0.75,  0.75/)
real, dimension (nspp) :: wturn   = (/ 0.031, 0.031,  0.01,  0.01/)
real, dimension (nspp) :: rturn   = (/ 0.031, 0.031,  2.00,  2.00/)
real, dimension (nspp) :: frcoeff = (/   0.5,   0.5,   0.5,   0.5/)
real, dimension (nspp) :: rrcoeff = (/   1.0,   1.0,   1.0,   1.0/)
real, dimension (nspp) :: nupc    = (/ 0.036, 0.036, 0.036, 0.036/)
integer, dimension (nspp) :: ptype = (/    1,     1,     3,     2/)
integer, dimension (nspp) :: weff  = (/    2,     2,     2,     2/)
! Ratio between gmax and Rubisco N (mol[H2O] m-2 s-1 / (kg[N] m-2)).
real, dimension (nspp) :: ngr     = (/1359.0,4652.0,1672.0,2223.0/)
! Cuticular conductance to CO2 (m s-1).
real, dimension (nspp) :: gmin    = (/4.81e-5,4.81e-5,4.81e-5,4.81e-5/)
! Initial slope of photosynthetic light response (mol m-1).
real, parameter :: alphas  = 0.04
! Chilling-day upper limit (oC).
real, parameter :: thold   = 5.0
real, allocatable, dimension (:) :: et_cat
real, dimension (nspp) :: pruba
real, dimension (nspp) :: slope
real, dimension (nspp) :: f1,f2,f3
real, dimension (nspp) :: fturn_save
real, dimension (nspp) :: rac
real, dimension (nspp) :: rah
real, dimension (nspp) :: rhr
real, dimension (nspp) :: mnppsp
logical :: local   ! Run for only local site?
logical :: rsf_out ! Write restart file (if not local)?
real    :: area    ! Plot area (m2)
real    :: idbh    ! Mean initial tree dbh (m).
real    :: idbhv   ! Fraction var. idbh
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
!----------------------------------------------------------------------!
logical :: summer
real    :: rgs ! Growth respiration (kg[C] plot-1 yr-1)
real    :: GPP_global,NPP_global,Cv_global,Cs_global
real    :: tsoil ! Soil temperature (oC)
real    :: ftsoil
real    :: tmind
real    :: Ndepo
real    :: nf
real    :: pptod
real    :: vm ! N:C ratio of surface metabolic litter
real    :: vn ! N:C ratio of soil metabolic litter
real    :: radpar ! mol[PAR] m-2 s-1
real    :: P_Pa
real    :: tfacs
real    :: ccair
real    :: tairK ! K
real    :: raf
real    :: radsw ! W m-2
real    :: vap ! Vapour pressure (Pa)
real    :: tairC ! oC
real    :: ppt   ! m dt-1
real    :: fc
real    :: tfaca1,tfaca2
real, dimension (2019-1700+1) :: cao ! Pa
integer :: i,i1
integer :: j,j1
integer :: kyr
integer :: kday,jdl
integer :: it
!----------------------------------------------------------------------!
integer :: nplots     ! No. plots/grid-box.
integer :: nind_total ! Total no. plants in simulation.
integer :: nind_max   ! Max. individuals/plot.
integer :: nland
!----------------------------------------------------------------------!
! netCDF parameters and variables.
!----------------------------------------------------------------------!
real, parameter :: fillvalue = 1.0e20
integer :: varid1,varid2,varid3
integer :: varid_C3GR,varid_C4GR,varid_BRCD,varid_NLEV
!----------------------------------------------------------------------!
real   , dimension (nlat,nd) :: dl
integer, dimension (nlon,nlat) :: land_index
real   , dimension (nlon,nlat) :: ddreq
real   , dimension (nlon,nlat) :: rltdl
real   , dimension (nlon,nlat) :: LAI_grid
real   , dimension (nlon,nlat) :: NPP_grid     ! NPP (kg[C] m-1 yr-1)
real   , dimension (nlon,nlat) :: Cv_grid
real   , dimension (nlon,nlat) :: Cs_grid
real   , dimension (nlon,nlat) :: Cv_C3GR_grid
real   , dimension (nlon,nlat) :: Cv_C4GR_grid
real   , dimension (nlon,nlat) :: Cv_BRCD_grid
real   , dimension (nlon,nlat) :: Cv_NLEV_grid
integer, dimension (nlon,nlat) :: cd           ! Chilling days (days).
real   , dimension (nlon,nlat) :: dd           ! Degree-days (oC).
integer, dimension (nlon,nlat) :: bgs
integer, dimension (nlon,nlat) :: egs
real, allocatable, dimension (:,:,:) :: tmp   ! K
real, allocatable, dimension (:,:,:) :: dswrf ! J m-2
real, allocatable, dimension (:,:,:) :: pre   ! mm 6hr-1
real, allocatable, dimension (:,:,:) :: pres  ! Pa
real, allocatable, dimension (:,:,:) :: spfh  ! kg kg-1
real, allocatable, dimension (:) :: lat ! degrees north
real, allocatable, dimension (:) :: lon ! degrees east
integer, allocatable, dimension (:,:,:) :: nind !****adf change to p
integer, allocatable, dimension (:) :: plot ! No. of plot
real   , allocatable, dimension (:) :: phenf
real   , allocatable, dimension (:) :: SWf
real   , allocatable, dimension (:) :: wlittc
real   , allocatable, dimension (:) :: wlittn
real   , allocatable, dimension (:) :: flittc
real   , allocatable, dimension (:) :: flittn
real   , allocatable, dimension (:) :: rlittc
real   , allocatable, dimension (:) :: rlittn
real   , allocatable, dimension (:) :: wfps
real   , allocatable, dimension (:) :: swp2
integer, allocatable, dimension (:) :: dd_flag
real   , allocatable, dimension (:) :: outflow
real   , allocatable, dimension (:) :: tfolK
real   , allocatable, dimension (:) :: gc
real   , allocatable, dimension (:) :: rcstom
integer, allocatable, dimension (:,:) :: p_plot
integer, allocatable, dimension (:,:,:) :: k_ind
integer, allocatable, dimension (:) :: ind
integer, allocatable, dimension (:) :: kspp
integer, allocatable, dimension (:) :: height
integer, allocatable, dimension (:) :: hbc
real   , allocatable, dimension (:) :: hdbh
integer, allocatable, dimension (:) :: alive
real   , allocatable, dimension (:) :: fPAR
real   , allocatable, dimension (:) :: fSW
real   , allocatable, dimension (:) :: farea
real   , allocatable, dimension (:) :: dbh
real   , allocatable, dimension (:) :: nfoliage
real   , allocatable, dimension (:) :: cfoliage
real   , allocatable, dimension (:) :: cfiner
real   , allocatable, dimension (:) :: cstore
real   , allocatable, dimension (:) :: cwood
real   , allocatable, dimension (:) :: lsap
real   , allocatable, dimension (:) :: nbswood
real   , allocatable, dimension (:) :: nheart
real   , allocatable, dimension (:) :: navail
real   , allocatable, dimension (:) :: nfiner
real   , allocatable, dimension (:) :: rcfoliage
real   , allocatable, dimension (:) :: rcwood
real   , allocatable, dimension (:) :: rcfiner
real   , allocatable, dimension (:) :: rnfoliage
real   , allocatable, dimension (:) :: rnbswood
real   , allocatable, dimension (:) :: rnfiner
real   , allocatable, dimension (:) :: rlold
real   , allocatable, dimension (:) :: lasa
real   , allocatable, dimension (:) :: nitf
real   , allocatable, dimension (:) :: jmaxfn
real   , allocatable, dimension (:) :: gmax
real   , allocatable, dimension (:) :: ratiol
real   , allocatable, dimension (:) :: cfact
real   , allocatable, dimension (:) :: ipfact
real   , allocatable, dimension (:) :: isfact
real   , allocatable, dimension (:) :: kzg
real   , allocatable, dimension (:) :: kSWzg
real   , allocatable, dimension (:) :: fPARiig
real   , allocatable, dimension (:) :: fSWiig
real   , allocatable, dimension (:) :: ball
real   , allocatable, dimension (:) :: foff
real   , allocatable, dimension (:) :: fon
real   , allocatable, dimension (:,:) :: larea ! Grid-box area (km2)
real   , allocatable, dimension (:,:) :: icwtr ! Ice/water fraction
real   , allocatable, dimension (:,:) :: kz
real   , allocatable, dimension (:,:) :: kSWz
real   , allocatable, dimension (:,:) :: fturn_plot
real   , allocatable, dimension (:,:,:) :: skzg
real   , allocatable, dimension (:,:,:) :: fPARtg
real   , allocatable, dimension (:,:,:) :: skSWzg
real   , allocatable, dimension (:,:,:) :: fSWtg
real   , allocatable, dimension (:,:,:) :: soilw1
real   , allocatable, dimension (:,:,:) :: soilw2
real   , allocatable, dimension (:,:,:) :: soilw3
real   , allocatable, dimension (:,:,:) :: snow
real   , allocatable, dimension (:,:,:) :: Cm
real   , allocatable, dimension (:,:,:) :: Cu
real   , allocatable, dimension (:,:,:) :: Cn
real   , allocatable, dimension (:,:,:) :: Cv
real   , allocatable, dimension (:,:,:) :: Ca
real   , allocatable, dimension (:,:,:) :: Cs
real   , allocatable, dimension (:,:,:) :: Cpa
real   , allocatable, dimension (:,:,:) :: Nu
real   , allocatable, dimension (:,:,:) :: Nm
real   , allocatable, dimension (:,:,:) :: Nv
real   , allocatable, dimension (:,:,:) :: Nn
real   , allocatable, dimension (:,:,:) :: Na
real   , allocatable, dimension (:,:,:) :: Ns
real   , allocatable, dimension (:,:,:) :: Npa
real   , allocatable, dimension (:,:,:) :: snmin
real   , allocatable, dimension (:,:,:) :: laip
real   , allocatable, dimension (:,:,:) :: soilC
real   , allocatable, dimension (:,:,:) :: soilN
real   , allocatable, dimension (:,:,:) :: swct ! Soil water hold. capacity (m)
real   , allocatable, dimension (:,:,:) :: swc1
real   , allocatable, dimension (:,:,:) :: swc2
real   , allocatable, dimension (:,:,:) :: swc3
end module shared
