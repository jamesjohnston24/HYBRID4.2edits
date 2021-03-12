!======================================================================!
module shared

!----------------------------------------------------------------------!
! Declares parameters and variables shared between subroutines.
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Control parameters from driver.in.
!----------------------------------------------------------------------!
logical :: local   ! Local (single grid-box) only?                (flag)
logical :: rsf_out ! Output final state to restart file?          (flag)
integer :: syr     ! Start year of simulation                       (CE)
integer :: eyr     ! End year of simulation                         (CE)
integer :: nplots  ! Number of plots per grid-box                    (n)
real :: area  ! Plot area                                           (m2)
real :: idbh  ! Initial tree DBH mean                                (m)
real :: idbhv ! Variation around idbh                         (fraction)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Other control and dimension parameters.
!----------------------------------------------------------------------!
logical, parameter :: wrclim = .FALSE. ! Write local climate?
integer, parameter :: nd        =   365 ! Days yr-1
integer, parameter :: ntimes    =  1460 ! Climate timepoints yr-1
integer, parameter :: nlon      =   720 ! No. longitudes
integer, parameter :: nlat      =   360 ! No. latitudes
integer, parameter :: nland_max = 67420 ! Nax. number of grid-boxes
integer, parameter :: mh        =   110 ! No. 1 m height layers in plot
integer, parameter :: nspp      =     8 ! No. GPTs
real, parameter :: fillvalue = 1.0e20
real, parameter :: eps       = 1.0e-8        ! Limit on accuracy!****adf
real, parameter :: dt        = 6.0*60.0*60.0 ! s 6hr-1
integer :: i,i1,i2    ! Longitude index              (index)
integer :: j,j1,j2    ! Latitude index               (index)
integer :: kyr        ! Year                            (CE)
integer :: kday       ! Day of year                    (day)
integer :: it         ! Timepoint in year            (index)
integer :: nind_total ! Total no. plants in simulation (ind)
integer :: limp       ! Max. trees/plot for regen.         (ind plot-1)
integer :: nind_max   ! Max. individuals/plot   (ind plot-1)
integer, dimension (nlon,nlat) :: land_index ! Index of grid-box (index)
integer, allocatable, dimension (:,:) :: p_plot ! Plot index (index)
integer, allocatable, dimension (:,:,:) :: k_ind ! Ind. index (index)
integer, allocatable, dimension (:) :: nind  ! No. individuals plot-1
integer, allocatable, dimension (:) :: plot  ! No. of plot
integer, allocatable, dimension (:) :: phenf ! Flag for plot foliage
real, allocatable, dimension (:) :: lat ! degrees north
real, allocatable, dimension (:) :: lon ! degrees east
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Physical constants.
!----------------------------------------------------------------------!
real, parameter :: zero = 0.0
real, parameter :: one  = 1.0
real, parameter :: sday = 86400.0 ! s day-1
real, parameter :: pi   = 3.14159265359    ! Ratio of circle circ/rad
real, parameter :: R    = 8.31446261815324 ! Gas constant (J K-1 mol-1)
real, parameter :: tf   = 273.15           ! Freezing point of water (K)
real, parameter :: m_water = 18.01528 ! MW water (g[water] mol-1)
real, parameter :: r_air_water   = 28.9647 / m_water ! MW[air]/MW[water]
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Physical parameters.
!----------------------------------------------------------------------!
real, dimension (nlat,0:nd) :: dl ! Daylength (s)
real, allocatable, dimension (:,:) :: larea ! Grid-box area        (km2)
real, allocatable, dimension (:,:) :: icwtr ! Ice/water frac. (fraction)
real, parameter :: swpmax = -0.033 ! Max. SWP                      (MPa)
real, parameter :: bsoil  =    5.0 ! Texture SPW parameter           (-)
real, parameter :: uwind =     5.0 ! Wind speed                  (m s-1)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Biological parameters.
!----------------------------------------------------------------------!
! Specified.
!----------------------------------------------------------------------!
!                                        C3GR     C4GR     BREV     BRCD     BRDD    NLEV      NLCD     NLDD
real, dimension (nspp) :: kpar    = (/   0.65,    0.65,    0.65,    0.65,    0.65,    0.50,    0.50,    0.50/)
real, dimension (nspp) :: ksw     = (/   0.48,    0.48,    0.48,    0.48,    0.48,    0.37,    0.37,    0.37/)
real, dimension (nspp) :: rhop    = (/   0.05,    0.05,    0.05,    0.05,    0.05,    0.03,    0.03,    0.03/)
real, dimension (nspp) :: rhos    = (/   0.20,    0.20,    0.20,    0.20,    0.20,    0.11,    0.11,    0.11/)
real, dimension (nspp) :: fturn   = (/ 0.0031,  0.0031,    0.33,    0.33,    0.33,    0.33,    0.33,    0.33/)
real, dimension (nspp) :: wturn   = (/ 0.0031,  0.0031,    0.01,    0.01,    0.01,    0.01,    0.01,    0.01/)
real, dimension (nspp) :: rturn   = (/ 0.0031,  0.0031,    2.00,    2.00,    2.00,    2.00,    2.00,    2.00/)
real, dimension (nspp) :: sla     = (/   36.0,    36.0,    36.0,    36.0,    36.0,    12.0,    12.0,    12.0/)
real, dimension (nspp) :: bark    = (/    1.0,     1.0,   0.033,   0.033,   0.033,    0.01,    0.01,    0.01/)
real, dimension (nspp) :: lsave   = (/   0.00,    0.00,  4167.0,  4167.0,  4167.0,  3333.0,  3333.0,  3333.0/)
real, dimension (nspp) :: ah      = (/   0.00,    0.00,   28.51,   28.51,   28.51,   32.95,   32.95,   32.95/)
real, dimension (nspp) :: bh      = (/   0.00,    0.00,  0.4667,  0.4667,  0.4667,  0.5882,  0.5882,  0.5882/)
real, dimension (nspp) :: stf     = (/   0.00,    0.00,    0.22,    0.22,    0.22,   0.222,   0.222,   0.222/)
real, dimension (nspp) :: rlratio = (/    1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0/)
real, dimension (nspp) :: frcoeff = (/    0.5,     0.5,     0.5,     0.5,     0.5,     0.5,     0.5,     0.5/)
real, dimension (nspp) :: rrcoeff = (/    1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0/)
real, dimension (nspp) :: woodd   = (/   0.00,    0.00,   305.0,   305.0,   305.0,   205.0,   205.0,   205.0/)
real, dimension (nspp) :: formf   = (/   0.00,    0.00,     0.6,     0.6,     0.6,    0.56,    0.56,    0.56/)
real, dimension (nspp) :: fsr     = (/  0.145,   0.145,   0.145,   0.145,   0.145,   0.145,   0.145,   0.145/)
real, dimension (nspp) :: frr     = (/   0.86,    0.86,    0.86,    0.86,    0.86,    0.86,    0.86,    0.86/)
real, dimension (nspp) :: live    = (/   0.00,    0.00,    0.17,    0.17,    0.17,  0.0708,  0.0708,  0.0708/)
real, dimension (nspp) :: storef  = (/   0.00,    0.00,    0.67,    0.67,    0.67,    0.67,    0.67,    0.67/)
real, dimension (nspp) :: nupc    = (/  0.036,   0.036,   0.036,   0.036,   0.036,   0.036,   0.036,   0.036/)
! Ratio between gmax and Rubisco N (mol[H2O] m-2 s-1 / (kg[N] m-2))
real, dimension (nspp) :: ngr     = (/ 1359.0,  4652.0,  1672.0,  1672.0,  1672.0,  2223.0,  2223.0,  2223.0/)
! Cuticular conductance to CO2 (m s-1)
real, dimension (nspp) :: gmin    = (/4.81e-5, 4.81e-5, 4.81e-5, 4.81e-5, 4.81e-5, 4.81e-5, 4.81e-5, 4.81e-5/)
real, dimension (nspp) :: pruba   = (/   0.67,    0.99,    0.67,    0.67,    0.67,    0.83,    0.67,    0.83/)
real, dimension (nspp) :: nrc     = (/   9.31,    9.31,    9.31,    9.31,    9.31,    9.31,    9.31,    9.31/)
real, dimension (nspp) :: d_leaf  = (/   0.04,    0.04,    0.04,    0.04,    0.04,    0.04,    0.04,    0.04/)
real, dimension (nspp) :: weffi   = (/    2.0,     2.0,     2.0,     2.0,     2.0,     2.0,     2.0,     2.0/)
real, dimension (nspp) :: ptypei = (/     1.0,     1.0,     2.0,     3.0,     4.0,     2.0,     3.0,     4.0/)
real, dimension (nspp) :: rgf     = (/   0.75,    0.75,    0.75,    0.75,    0.75,    0.75,    0.75,    0.75/)
real, dimension (nspp) :: wmf     = (/    0.1,     0.1,     0.1,     0.1,     0.1,     0.1,     0.1,     0.1/)
real, parameter :: thold = 5.0 ! Chilling-day/degree-day limit (oC)
real, parameter :: nfp   = 0.02 ! Initial tree foliage N (fraction)
real, dimension (nlon,nlat) :: ddreq ! Degree-day requirement (oC)
real, dimension (nlon,nlat) :: rltdl ! Daylength requirement   (s)
!----------------------------------------------------------------------!
! Derived.
!----------------------------------------------------------------------!
integer, dimension (nspp) :: weff  ! Water effect on N uptake    (index)
integer, dimension (nspp) :: ptype ! Phenology type              (index)
real, dimension (nspp) :: crratio ! Grass stucture/foliage (ratio)
! Slope in fraction non-photosynthetic N / foliage N (fraction)
real, dimension (nspp) :: slope
! Allometric factors for trees.
real, dimension (nspp) :: f1,f2,f3
! Baseline foliage turnover rate (fraction day-1)
real, dimension (nspp) :: fturn_save
!Canopy boundary layer resistance to CO2 (s m-1)
real, dimension (nspp) :: rac
! Canopy boundary layer resistance to sensible heat (s m-1)
real, dimension (nspp) :: rah
! Canopy boundary layer resistance to convective and radiative heat
! (s m-1)
real, dimension (nspp) :: rhr
real :: raf ! Adjustment to psychrometric constant (unitless)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------! Soil parameters.
!----------------------------------------------------------------------!
real, parameter :: d_one = 0.2 ! Max. SWC of top layer for roots (m)
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

!----------------------------------------------------------------------!
! State variables.
!----------------------------------------------------------------------!
integer, dimension (nlon,nlat) :: cd  ! Chilling days (days)
integer, dimension (nlon,nlat) :: bgs ! Begining of summer (day)
integer, dimension (nlon,nlat) :: egs ! End of summer (day)
real, dimension (nlon,nlat) :: dd ! Degree-days (oC)
!----------------------------------------------------------------------!
! Plot (indexed by p).
!----------------------------------------------------------------------!
real, allocatable, dimension (:) :: soilw1 ! Layer 1 soil water      (m)
real, allocatable, dimension (:) :: soilw2 ! Layer 2 soil water      (m)
real, allocatable, dimension (:) :: soilw3 ! Layer 3 soil water      (m)
real, allocatable, dimension (:) :: snow   ! Snow water eq.          (m)
real, allocatable, dimension (:) :: Cm ! Surf. metab. soil C (kg[C] m-2)
real, allocatable, dimension (:) :: Cu ! Surf. struc. soil C (kg[C] m-2)
real, allocatable, dimension (:) :: Cn ! Soil metab. soil C  (kg[C] m-2)
real, allocatable, dimension (:) :: Cv ! Soil struc. soil C  (kg[C] m-2)
real, allocatable, dimension (:) :: Ca ! Active soil C       (kg[C] m-2)
real, allocatable, dimension (:) :: Cs ! Slow soil C         (kg[C] m-2)
real, allocatable, dimension (:) :: Cpa ! Passive soil C     (kg[C] m-2)
real, allocatable, dimension (:) :: Nm ! Surf. metab. soil N (kg[N] m-2)
real, allocatable, dimension (:) :: Nu ! Surf. struc. soil N (kg[N] m-2)
real, allocatable, dimension (:) :: Nv ! Soil struc. soil N  (kg[N] m-2)
real, allocatable, dimension (:) :: Nn ! Soil metab. soil N  (kg[N] m-2)
real, allocatable, dimension (:) :: Na ! Active soil N       (kg[N] m-2)
real, allocatable, dimension (:) :: Ns ! Slow soil N         (kg[N] m-2)
real, allocatable, dimension (:) :: Npa ! Passive soil N     (kg[N] m-2)
real, allocatable, dimension (:) :: snmin ! Mineral N        (kg[N] m-2)
!----------------------------------------------------------------------!
! Individual (indexed by k).
!----------------------------------------------------------------------!
integer, allocatable, dimension (:) :: ind   ! Individual index  (index)
integer, allocatable, dimension (:) :: kspp  ! GPT-index         (index)
integer, allocatable, dimension (:) :: alive ! Alive?             (flag)
integer, allocatable, dimension (:) :: hbc ! Height to base of crown (m)
real, allocatable, dimension (:) :: dbh  ! Tree stem DBH/min. gr LAI (m)
real, allocatable, dimension (:) :: hdbh ! Heartwood DBH             (m)
real, allocatable, dimension (:) :: cstore   ! Stored C    (kg[C] ind-1)
real, allocatable, dimension (:) :: navail   ! Available N (kg[N] ind-1)
real, allocatable, dimension (:) :: cfoliage ! Foliage C   (kg[C] ind-1)
real, allocatable, dimension (:) :: nfoliage ! Foliage N   (kg[N] ind-1)
real, allocatable, dimension (:) :: nbswood  ! Structure N (kg[N] ind-1)
real, allocatable, dimension (:) :: nheart   ! Heartwood N (kg[N] ind-1)
real, allocatable, dimension (:) :: cfiner   ! Fine root C (kg[C] ind-1)
real, allocatable, dimension (:) :: nfiner   ! Fine root N (kg[N] ind-1)
real, allocatable, dimension (:) :: lasa     ! FA/SA            (m2 m-2)
real, allocatable, dimension (:) :: nitf     ! T factor for A (fraction)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Plot variables (indexed by p; not state variables).
!----------------------------------------------------------------------!
real, allocatable, dimension (:) :: SWf    ! Plot abs. SW factor (ratio)
real, allocatable, dimension (:) :: wlittc ! Structure lit C (kg[C] m-2)
real, allocatable, dimension (:) :: wlittn ! Structure lit N (kg[C] m-2)
real, allocatable, dimension (:) :: flittc ! Foliage lit C   (kg[C] m-2)
real, allocatable, dimension (:) :: flittn ! Foliage lit N   (kg[C] m-2)
real, allocatable, dimension (:) :: rlittc ! Fine root lit C (kg[C] m-2)
real, allocatable, dimension (:) :: rlittn ! Fine root lit N (kg[C] m-2)
real, allocatable, dimension (:) :: laip   ! Plot LAI           (m2 m-2)
real, allocatable, dimension (:) :: soilC  ! Total soil C    (kg[C] m-2)
real, allocatable, dimension (:) :: soilN  ! Total soil N    (kg[N] m-2)
real, allocatable, dimension (:) :: swct   ! Total SWHC              (m)
real, allocatable, dimension (:) :: swc1   ! SWHC of layer 1         (m)
real, allocatable, dimension (:) :: swc2   ! SWHC of layer 2         (m)
real, allocatable, dimension (:) :: swc3   ! SWHC of layer 3         (m)
!----------------------------------------------------------------------!
! Plot variables (indexed by kp).
!----------------------------------------------------------------------!
integer, allocatable, dimension (:) :: dd_flag ! DD diagnostic    (flag)
real, allocatable, dimension (:) :: wfps ! Water-filled pore space   (%)
real, allocatable, dimension (:) :: swp1 ! Grass SWP               (MPa)
real, allocatable, dimension (:) :: swp2 ! Tree SWP                (MPa)
real, allocatable, dimension (:) :: outflow ! Water outflow    (m day-1)
real, allocatable, dimension (:) :: tfolK   ! Foliage temperature    (K)
real, allocatable, dimension (:) :: tr_rat  ! Top transp.     (fraction)
!----------------------------------------------------------------------!
! Plot variable (indexed by kp,ksp).
!----------------------------------------------------------------------!
real, allocatable, dimension (:,:) :: fturn_plot ! Fol. tover (fraction)
!----------------------------------------------------------------------!
! Plot variables (indexed by i,j,kp).
!----------------------------------------------------------------------!
real, allocatable, dimension (:,:,:) :: skzg   ! APAR of grass   (ratio)
real, allocatable, dimension (:,:,:) :: skSWzg ! ASW of grass    (ratio)
real, allocatable, dimension (:,:,:) :: fPARtg ! fPAR of grass   (ratio)
real, allocatable, dimension (:,:,:) :: fSWtg  ! fSW of grass    (ratio)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Individual variables (indexed by k; not state variables).
!----------------------------------------------------------------------!
integer, allocatable, dimension (:) :: foff ! Flag, foliage on    (flag)
integer, allocatable, dimension (:) :: fon  ! Flag, foliage off   (flag)
integer, allocatable, dimension (:) :: height ! Total height         (m)
real, allocatable, dimension (:) :: et_cat ! Rub cat site c't  (mol m-2)
real, allocatable, dimension (:) :: fPAR   ! fPAR of ind.     (fraction)
real, allocatable, dimension (:) :: fSW    ! fSW of ind.      (fraction)
real, allocatable, dimension (:) :: farea  ! Foliage area     (m2 ind-1)
real, allocatable, dimension (:) :: lsap   ! Living struc. (kg[C] ind-1)
real, allocatable, dimension (:) :: cwood  ! Structural C  (kg[C] ind-1)
real, allocatable, dimension (:) :: jmaxfn ! Jmax factor  (mol[chl] m-2)
real, allocatable, dimension (:) :: gmax   ! Max.st.c.(mol[H2O] m-2 s-1)
real, allocatable, dimension (:) :: ratiol ! fPARb/fPAR          (ratio)
real, allocatable, dimension (:) :: cfact  ! Ind/m2 A scalar  (m2 ind-1)
real, allocatable, dimension (:) :: ipfact ! APAR scalar         (ratio)
real, allocatable, dimension (:) :: isfact ! ISW scalar          (ratio)
real, allocatable, dimension (:) :: kzg    ! Grass layer APAR    (ratio)
real, allocatable, dimension (:) :: kSWzg  ! Grass layer ASW     (ratio)
real, allocatable, dimension (:) :: fPARiig ! fPARb of grass     (ratio)
real, allocatable, dimension (:) :: fSWiig ! fSWb of grass       (ratio)
real, allocatable, dimension (:) :: ball   ! Basal layer C bal.  (kg[C])
!----------------------------------------------------------------------!
! Individual variables (indexed by kp,ki).
!----------------------------------------------------------------------!
real, allocatable, dimension (:,:) :: gc ! Canopy cond.     (m[CO2] s-1)
!----------------------------------------------------------------------!
! Individual variables (indexed by m,ki).
!----------------------------------------------------------------------!
real, allocatable, dimension (:,:) :: kz   ! APAR (ratio)
real, allocatable, dimension (:,:) :: kSWz ! ASW  (ratio)
!----------------------------------------------------------------------!
! Previous values.
!----------------------------------------------------------------------!
real, allocatable, dimension (:) :: rcfoliage ! cfoliage   (kg[C] ind-1)
real, allocatable, dimension (:) :: rnfoliage ! nfoliage   (kg[N] ind-1)
real, allocatable, dimension (:) :: rcwood    ! cwood      (kg[C] ind-1)
real, allocatable, dimension (:) :: rnbswood  ! nbswood    (kg[N] ind-1)
real, allocatable, dimension (:) :: rcfiner   ! cfiner     (kg[C] ind-1)
real, allocatable, dimension (:) :: rnfiner   ! nfiner     (kg[N] ind-1)
real, allocatable, dimension (:) :: rlold     ! farea         (m2 ind-1)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Other dynamic variables.
!----------------------------------------------------------------------!
logical :: summer ! Foliage display if cold-deciduous (flag)
real :: radsw  ! Downwelling shortwave radiation           (W m-2)
real :: radpar ! Downwelling PAR                (mol[PAR] m-2 s-1)
real :: tairK  ! Air temperature                               (K)
real :: tair   ! Air temperature                              (oC)
real :: ppt    ! Precipitation                           (m 6hr-1)
real :: P_Pa   ! Atmospheric pressure                         (Pa)
real :: vap    ! Atmospheric water vapour pressure            (Pa)
real :: vapr   ! Atmospheric water vapour pressure           (hPa)
real :: ccair  ! Atmospheric CO2                         (mol m-3)
real :: tsoil  ! Soil temperature                             (oC)
real :: tfacs  ! Fine root resp. T factor      (kg[C] kg[N]-1 s-1)
real :: tfaca1 ! Foliage respiration T factor  (kg[C] kg[N]-1 s-1)
real :: tfaca2 ! Structure respiration T factor  (kg[C] kg[C] s-1)
real :: tmind  ! Min. temperature of previous 24 hr           (oC)
real :: Ndepo  ! N deposition                    (kg[N] m-2 day-1)
real :: nf     ! Baseline N deposition rate      (kg[N] m-1 day-1)
real :: vm     ! N:C ratio of surface metabolic litter     (ratio)
real :: vn     ! N:C ratio of soil metabolic litter        (ratio)
real :: fc     ! CO2 effect on stomatal conductance       (scalar)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Input physical variables.
!----------------------------------------------------------------------!
real, dimension (2019-1700+1) :: cao ! Atmospheric CO2 (Pa)
real, allocatable, dimension (:,:,:) :: tmp   ! Air temperature      (K)
real, allocatable, dimension (:,:,:) :: dswrf ! Downwelling SW   (J m-2)
real, allocatable, dimension (:,:,:) :: pre   ! Precipitation (mm 6hr-1)
real, allocatable, dimension (:,:,:) :: pres  ! Atm. pressure       (Pa)
real, allocatable, dimension (:,:,:) :: spfh  ! Atm. sp. hum.  (kg kg-1)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Diagnostic variables.
!----------------------------------------------------------------------!
real, dimension (nspp) :: mgppsp! GPP by GPT if local (kg[C] m-2 yr-1)
real, dimension (nspp) :: mnppsp! NPP by GPT if local (kg[C] m-2 yr-1)
real :: rgs ! Growth respiration (kg[C] grid-box-1 yr-1)
real :: GPP_global ! Global GPP          (kg[C] yr-1)
real :: NPP_global ! Global NPP          (kg[C] yr-1)
real :: Cv_global  ! Global vegetation C (kg[C] yr-1)
real :: Cs_global  ! Global soil C       (kg[C] yr-1)
real, dimension (nlon,nlat) :: LAI_grid     ! LAI         (m2 m-2)
real, dimension (nlon,nlat) :: GPP_grid     ! GPP (kg[C] m-1 yr-1)
real, dimension (nlon,nlat) :: NPP_grid     ! NPP (kg[C] m-1 yr-1)
real, dimension (nlon,nlat) :: Cv_grid      ! Veg. C   (kg[C] m-2)
real, dimension (nlon,nlat) :: Cs_grid      ! Soil C   (kg[C] m-2)
real, dimension (nlon,nlat) :: Cv_C3GR_grid ! C3GR C   (kg[C] m-2)
real, dimension (nlon,nlat) :: Cv_C4GR_grid ! C4GR C   (kg[C] m-2)
real, dimension (nlon,nlat) :: Cv_BREV_grid ! BREV C   (kg[C] m-2)
real, dimension (nlon,nlat) :: Cv_BRCD_grid ! BRCD C   (kg[C] m-2)
real, dimension (nlon,nlat) :: Cv_BRDD_grid ! BRDD C   (kg[C] m-2)
real, dimension (nlon,nlat) :: Cv_NLEV_grid ! NLEV C   (kg[C] m-2)
real, dimension (nlon,nlat) :: Cv_NLCD_grid ! NLCD C   (kg[C] m-2)
real, dimension (nlon,nlat) :: Cv_NLDD_grid ! NLDD C   (kg[C] m-2)
!----------------------------------------------------------------------!

end module shared
!======================================================================!

