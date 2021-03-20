!======================================================================!
subroutine hydrology
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Computes N uptake temperature effect (why here?).
! Computes effect soil water potential for each individual.
! Computes stomatal conductance for each individual.
! Computes canopy temperature for each plot.
! Updates soil water in each layer and snowpack.
! Computes proportion of water lost as outflow.
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
use shared
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Parameters only used in this subroutine.
!----------------------------------------------------------------------!
real, parameter :: cp = 1012.0 ! Sp. heat capacity of air   (J kg-1 K-1)
real, parameter :: ftmin  = 0.0   ! Parameter for ft                (oC)
real, parameter :: ftmax  = 40.0  ! Parameter for ft                (oC)
real, parameter :: fkfive = 18.35 ! Parameter for ft                (oC)
real, parameter :: fkeight = (ftmax - fkfive) / (fkfive - ftmin)
real, parameter :: intc   = 0.0005 * dt / sday ! m LAI-1 dt-1
real, parameter :: smeltc = 0.0007 * dt / sday ! m oC-1 dt-1
real, parameter :: fout   = 0.9    * dt / sday ! fraction dt-1
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Miscellaneous variables used in this subroutine.
!----------------------------------------------------------------------!
integer :: kp  ! Plot number in grid-box                          (plot)
integer :: p   ! Plot number in simulation                        (plot)
integer :: ki  ! Individual number in plot                         (ind)
integer :: k   ! Individual number in simulation                   (ind)
integer :: ksp ! GPT number                                        (GPT)
!----------------------------------------------------------------------!
real, dimension (nind_max) :: rcstom ! Stomatal res. to CO2 (m[CO2] s-1)
!----------------------------------------------------------------------!
real :: rhoa   ! Density of dry air                        (kg[air] m-3)
real :: lambda ! Latent heat of vapourisation of H2O            (J kg-1)
real :: gamma  ! Psychrometer 'constant'                   (mol m-3 K-1)
real :: gammas ! Adjusted psychrometer `constant'          (mol m-3 K-1)
real :: ef1    ! Leaf temperature calculation factor           (mol J-1)
real :: ts     ! Apparent radiative temperature of atmosphere        (K)
real :: psinf  ! Isothermal net radiation factor                 (W m-2)
real :: cwsta  ! Saturation conc'n of water at air temperature (mol m-3)
real :: cwstl  ! Saturation conc'n of water close to air temp. (mol m-3)
real :: cwa    ! Conc'n of water vapour in air outside LBL     (mol m-3)
real :: dwo    ! Water vapour conc'n deficit of air outside BL (mol m-3)
real :: s      ! Slope of saturation vapour pressure curve (mol m-3 K-1)
real :: fdq    ! Humidity stomatal cond. response function    (fraction)
real :: ft     ! Temperature stomatal conductance response    (fraction)
real :: fdqftfc ! Stom. factor to allow for vpd, temp., CO2   (fraction)
real :: tt     ! Factor to calculate fdqftfc
real :: sow1   ! Soil water in top two layers                        (m)
real :: sc1    ! Soil water holding capacity in top two layers       (m)
real :: sw1    ! Soil water potential of top two layers            (MPa)
real :: sw2    ! Soil water potential of lowest layer              (MPa)
real :: tr_grass ! Plot conductance for grass transpiration (m[H2O] s-1)
real :: tr_tot ! Plot conductance for total transpiration   (m[H2O] s-1)
real :: psin   ! Total plot isoth. net rad. on ground area basis (W m-2)
real :: cgc    ! Total canopy conductance                   (m[H2O] s-1)
real :: st     ! SW radiation absorbed in upper crown layer      (W m-2)
real :: fst    ! Stomatal conductance solar-radiation factor  (fraction)
real :: fdt    ! Stomatal conductance soil-water factor       (fraction)
real :: gs     ! Stomatal conductance to H2O          (mol[H2O] m-2 s-1)
real :: rwi    ! Total res. across LBL and leaf surface     (s m[H2O]-1)
real :: rw     ! rwi on ground area basis                   (m[H2O] s-1)
real :: raws   ! BL conductance to water of bare ground     (m[H2O] s-1)
real :: rsoil  ! Resistance to water loss of bare ground    (m[H2O] s-1)
real :: betae  ! Soil water control on rsoil                  (fraction)
real :: veg_frac ! Vegetation canopy fraction of plot         (fraction)
real :: rcif   ! Baseline res. from leaf surface to mesophyll (s m[CO2])
real :: rci    ! Resistance from leaf surface to mesophyll    (s m[CO2])
real :: rcn    ! Resistance from outside LBL to mesophyll     (s m[CO2])
real :: cwstf  ! Saturation conc'n of water at leaf temp' (mol[H2O] m-3)
real :: vpd    ! Canopy-to-air water vapour deficit       (mol[H2O] m-3)
real :: ET     ! Rate of transpiration from plot      (mol[H2O] m-2 s-1)
real :: Ek     ! Rate of transpiration from plot      (kg[H2O] m-2 dt-1)
real :: lmbae  ! Mean energy used in transpiration               (W m-2)
real :: transp ! Plot transpiration, assuming 1 g = 1 cm-3 (m[H2O] dt-1)
real :: pint   ! Potential interception                         (m dt-1)
real :: aint   ! Actual interception                            (m dt-1)
real :: pevap  ! Potential evaporation                          (m dt-1)
real :: pevap1 ! Potential evaporation                (mol[H2O] m-2 s-1)
real :: aevap  ! Actual evaporated intercepted rain             (m dt-1)
real :: aevapr ! Actual evaporated intercepted rain   (kg[H2O] m-2 dt-1)
real :: srain  ! Total rain reaching soil                       (m dt-1)
real :: snowi  ! Snow                                           (m dt-1)
real :: smelt  ! Snowmelt                                       (m dt-1)
real :: transp1 ! Transpiration from first soil layer           (m dt-1)
real :: transp2 ! Transpiration from second soil layer          (m dt-1)
real :: transp3 ! Transpiration from third soil layer           (m dt-1)
real :: out1   ! Outflow from first to second layer             (m dt-1)
real :: out2   ! Outflow from second to third layer             (m dt-1)
real :: outt   ! Outflow from plot                              (m dt-1)
real :: ploss  ! Outflow as proportion of plot water          (fraction)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Density of dry air (kg m-3).
!----------------------------------------------------------------------!
rhoa = 2.42 - 4.12e-3 * tairK
!-----------------------------------------------------------------!
! Latent heat of vapourisation of H2O (increases by 3% between
! 30 oC and 0oC; J kg-1).
!-----------------------------------------------------------------!
lambda = 3.15e+6 - 2.375e+3 * tairK
!-----------------------------------------------------------------!
! Psychrometeric 'constant' (mol m-3 K-1).
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
!----------------------------------------------------------------------!
! Atmospheric concentration of CO2 (mol m-3).
!----------------------------------------------------------------------!
ccair = cao (kyr-1699) / (R * tairK)
!----------------------------------------------------------------------!
! Saturation concentration of water at air temperature;
! after Jones 1992 (mol m-3).
!----------------------------------------------------------------------!
cwsta = 613.78 * exp (17.502 * (tairK- tf) / &
        (tairK - 32.18)) / (R * tairK)
!----------------------------------------------------------------------!
! Saturation concentration of water close to air temperature
! (mol m-3).
!----------------------------------------------------------------------!
cwstl = 613.78 * exp (17.502 * ((tairK + 0.1) - tf) / &
        ((tairK + 0.1) - 32.18)) / (R * (tairK + 0.1))
!----------------------------------------------------------------------!
! Concentration of water vapour in air outside leaf boundary
! layer (mol m-3).
!----------------------------------------------------------------------!
cwa = vapr * 100.0 / (R * tairK)
cwa = max (zero , cwa)
cwa = min (cwsta, cwa)
!----------------------------------------------------------------------!
! Water vapour concentration deficit of air outside boundary
! layer (mol m-3).
!----------------------------------------------------------------------!
dwo = cwsta - cwa
!----------------------------------------------------------------------!
! Slope of saturation vapour pressure curve (mol m-3 K-1).
!----------------------------------------------------------------------!
s = 10.0 * (cwstl - cwsta)
!----------------------------------------------------------------------!
! Humidity stomatal conductance response function (fraction).
! Full closure at 3000 Pa and passing through 0.2348 at 1500 Pa.
!----------------------------------------------------------------------!
fdq = 0.46955 - 0.3815 * dwo
fdq = max (zero, fdq)
!----------------------------------------------------------------------!
! Temperature stomatal conductance response (fraction).
!----------------------------------------------------------------------!
tt = tair
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
!----------------------------------------------------------------------!
! Overall stomatal factor to allow for vpd, temperature, and CO2.
!----------------------------------------------------------------------!
fdqftfc = fdq * ft * fc
!----------------------------------------------------------------------!
do kp = 1, nplots
 p = p_plot (land_index(i,j),kp)
 !---------------------------------------------------------------------!
 ! Total plot isothermal net radiation on ground area basis
 ! (W m-2).
 !---------------------------------------------------------------------!
 psin = SWf (p) * (radsw + psinf)
 !---------------------------------------------------------------------!
 ! Total canopy conductance (m[H2O] s-1).
 !---------------------------------------------------------------------!
 cgc = zero
 !---------------------------------------------------------------------!
 ! Transpiration factors to determine tree vs. grass.
 !---------------------------------------------------------------------!
 tr_grass = zero
 tr_tot   = zero
 do ki = 1, nind (p)
  !--------------------------------------------------------------------!
  ! Get index of individual.
  !--------------------------------------------------------------------!
  k = k_ind (land_index(i,j),kp,ki)
  if (alive (k) == 1) then
  ksp = kspp (k)
  !--------------------------------------------------------------------!
  ! SW radiation absorbed in upper crown layer (W m-2).
  !--------------------------------------------------------------------!
  st = radsw * isfact (k)
  !--------------------------------------------------------------------!
  ! Stomatal conductance solar-radiation factor (fraction).
  !--------------------------------------------------------------------!
  fst = 1.1044 * st / (st + 104.4)
  !--------------------------------------------------------------------!
  ! Stomatal conductance soil-water factor (fraction).
  !--------------------------------------------------------------------!
  if ((ki == 1) .or. (ki == 2)) then
   ! Soil water factor for grass.
   if (swp1 (kp) > (-1.5)) fdt = 1.155 + 0.77 * &
	 (swp1 (kp) - 0.015 * float (height (k)))
    if (swp1 (kp) >  (-0.2)) fdt = one
    if (swp1 (kp) <= (-1.5)) fdt = zero
  else
   ! Soil water factor for trees.
   if (swp2 (kp) > (-1.5)) fdt = 1.155 + 0.77 * &
	 (swp2 (kp) - (lsave (ksp) / lasa (k)) * 0.005 * &
	 float (height (k)) - 0.01 * float (height (k)))
   if (swp2 (kp) >  (-0.2)) fdt = one
   if (swp2 (kp) <= (-1.5)) fdt = zero
  end if
  !--------------------------------------------------------------------!
  ! Stomatal conductance to H2O (mol[H2O] m-2 s-1).
  !--------------------------------------------------------------------!
  gs = gmax (k) * fst * fdt * fdqftfc
  gs = max (gmin (ksp), gs)
  !--------------------------------------------------------------------!
  ! Stomatal resistance to CO2 (m s-1).
  !--------------------------------------------------------------------!
  rcstom (ki) = 1.42 / gs
  !--------------------------------------------------------------------!
  ! Total resistance to water transfer across leaf boundary layer
  ! and leaf surface (s m[H2O]-1).
  !--------------------------------------------------------------------!
  if (cfact (k) > eps) then
   !-------------------------------------------------------------------!
   rwi = (0.607 * rac (ksp) + one / gs) * area / cfact (k)
   !-------------------------------------------------------------------!
   ! Sum canopy conductance across individuals (m[H2O] s-1).
   !-------------------------------------------------------------------!
   cgc = cgc + one / rwi
   !-------------------------------------------------------------------!
   ! Sum relative amount of transpiration due to trees.
   !-------------------------------------------------------------------!
   tr_tot = tr_tot + one / rwi
   if (ki <= 2) tr_grass = tr_grass + one / rwi
   !-------------------------------------------------------------------!
  end if
  else
   rcstom (ki) = one / eps
  end if ! alive
 end do ! ki
 !---------------------------------------------------------------------!
 ! tr_rat is converted from flag to the fraction of plot
 ! transpiration from top 2 layers (adjusted if trees taking water
 ! from bottom layer, i.e. tr_rat is set from previous day in this
 ! routine).
 !***adf I think tr_rat needs to carry over sites.
 !---------------------------------------------------------------------!
 if (tr_rat (kp) < 0.1) then
  ! Tree transpiration from bottom (third) layer only.
  if (tr_tot > eps) then
   tr_rat (kp) = tr_grass / tr_tot
  else
   tr_rat (kp) = 0.5
  end if
 end if
 !---------------------------------------------------------------------!
 ! Total resistance to water transfer across leaf boundary layer
 ! and leaf surface on ground area basis (m[H2O] s-1).
 !---------------------------------------------------------------------!
 if (cgc > eps) then
  rw = one / cgc
 else
  rw = 0.607 * rac (ksp) + one / gmin (1)
 end if
 !---------------------------------------------------------------------!
 ! raws is boundary layer conductance to water of bare ground
 ! calculated from p. 68 of Jones.
 !---------------------------------------------------------------------!
 raws = 170.5
 if (snow (p) > zero) then
  rsoil = raws + 150.0
 else
  if (swc1 (p) > eps) then
   betae = (soilw1 (p) / swc1 (p) - 0.19) / 0.016
  else
   betae = 0.01
  end if
  if (tsoil <= zero) betae = 0.01
  betae = min (one, betae)
  betae = max (eps, betae)
  rsoil = raws * (one + (one - betae) / betae)
 end if
 rsoil = max (eps, rsoil)
 !---------------------------------------------------------------------!
 ! Vegetation canopy fraction of plot (fraction).
 !---------------------------------------------------------------------!
 veg_frac = 0.6 + 0.0875 * laip (p)
 if (laip (p) == zero) veg_frac = zero
 veg_frac = min (0.95, veg_frac)
 !---------------------------------------------------------------------!
 ! Update rw to include soil.
 !---------------------------------------------------------------------!
 rw = one / (veg_frac / rw + (one - veg_frac) / rsoil)
 !---------------------------------------------------------------------!
 ! Mean plot canopy temperature (K) (ef1 is gamma / (rhoa * cp)).
 !---------------------------------------------------------------------!
 ksp = 1 ! Use grass for now for boundary conductance.
 !****adf I cannot see how psin is correct for here.
 tfolK (kp) = tairK + (rw * psin * ef1 - dwo) / &
    (gamma * (rw + rah (ksp)) / rhr (ksp) + s)
 ! Resistance to CO2 transfer from leaf surface to outside
 ! mesophyll liquid phase factor (83.0e3 * 300.0e-6). 300.0e-6 m
 ! is leaf thickness. 
 rcif = 24.9
 ! Resistance to CO2 transfer from leaf surface to outside
 ! mesophyll liquid phase (s m-1).
 rci = rcif * (tfolK (kp) / 293.15) ** (-1.75)
 !---------------------------------------------------------------------!
 ! Individual stomatal conductance to CO2 of uppermost leaf.
 !---------------------------------------------------------------------!
 do ki = 1, nind (p)
  !--------------------------------------------------------------------!
  ! Get index of individual.
  !--------------------------------------------------------------------!
  k = k_ind (land_index(i,j),kp,ki)
  if (alive (k) == 1) then
  ksp = kspp (k)
  ! Total resistance to CO2 transfer from outside leaf boundary
  ! layer to outside mesophyll liquid phase (s m-1).
  rcn = rac (ksp) + rcstom (ki) + rci
  rcn = max (eps, rcn)
  ! Total conductance to CO2 transfer from outside leaf boundary
  ! layer to outside mesophyll liquid phase of upper layer
  ! (m s-1), used in PGEN.
  gc (kp,ki) = one / rcn
  !--------------------------------------------------------------------!
  else
   gc (kp,ki) = eps
  end if ! alive
 end do ! ki
 !---------------------------------------------------------------------!
 ! Saturation concentration of water at leaf temperature
 ! (mol[H2O] m-3).
 !---------------------------------------------------------------------!
 cwstf = 613.78 * exp (17.502 * (tfolK (kp) - tf) / &
    (tfolK (kp) - 32.18)) / (R * tfolK (kp))
 !---------------------------------------------------------------------!
 ! Canopy-to-air water vapour deficit (mol[H2O] m-3).
 !---------------------------------------------------------------------!
 vpd = cwstf - cwa
 vpd = max (zero, vpd)
 !---------------------------------------------------------------------!
 ! Rate of transpiration from plot (mol[H2O] m-2 s-1).
 !---------------------------------------------------------------------!
 ET = vpd / rw
 !---------------------------------------------------------------------!
 ! Rate of transpiration from plot (kg[H2O] m-2 dt-1).
 !---------------------------------------------------------------------!
 Ek = ET * m_water * dt / 1.0e3
 !---------------------------------------------------------------------!
 ! Mean energy used in transpiration (W m-2).
 !---------------------------------------------------------------------!
 lmbae = ET * m_water * lambda / 1.0e3
 !---------------------------------------------------------------------!
 ! Plot transpiration, assuming 1 g = 1 cm-3 (m[H2O] dt-1).
 transp = Ek / 1000.0
 transp = min (transp, (soilw1 (p) + soilw2 (p) + &
          soilw3 (p)))
 !---------------------------------------------------------------------!
 ! Soil water balance.
 if (tair > eps) then ! Rain.
  ! Potential interception (m dt-1).
  pint = laip (p) * intc
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
  pevap = 0.001 * dt * pevap1 * m_water / 1.0e3
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
 smelt = max (zero, (smeltc * tair))
 ! Actual snowmelt (m dt-1).
 smelt = min (smelt, snow (p))
 ! Transpiration from first two soil layers is weighted by soil
 ! water content (assumes root distribution follows soil
 ! water). tr_rat is the fraction of total transpiration from
 ! the first two layers.
 if ((soilw1 (p) + soilw2 (p)) > eps) then
  transp1 = tr_rat (kp) * transp * soilw1 (p) / &
  (soilw1 (p) + soilw2 (p))
 else
  transp1 = zero
 end if
 transp2 = tr_rat (kp) * transp - transp1
 transp3 = transp - (transp1 + transp2)
 ! Update soil water in top layer (m).
 soilw1 (p) = soilw1 (p) + srain + smelt - transp1
 ! If too much transpiration for top layer to cope with, make
 ! top layer water zero, and subtract required amount from
 ! second layer.
 if (soilw1 (p) < zero) then
  transp2 = transp2 - soilw1 (p)
  soilw1 (p) = zero
 end if
 ! Outflow from first to second layer if top layer has more water
 ! than capacity (m dt-1).
 out1 = max (zero, (soilw1 (p) - swc1 (p)))
 ! Restrict outflow from top layer to stop second layer having
 ! more than saturation.
 out1 = min (out1, swc2 (p) / 0.58 - soilw2 (p))
 out1 = max (out1, zero)
 ! Update first (top) layer.
 soilw1 (p) = soilw1 (p) - out1
 ! Update second (middle) layer.
 soilw2 (p) = soilw2 (p) + out1 - transp2
 ! If too much transpiration for second layer to cope with,
 ! subtract from third layer.
 if (soilw2 (p) < zero) then
  transp3 = transp3 - soilw2 (p)
  soilw2 (p) = zero
 end if
 ! Outflow from second to third layer if second layer has more
 ! water than capacity (m dt-1).
 out2 = max (zero, (soilw2 (p) - swc2 (p)))
 ! Restrict outflow from second layer to stop third layer having
 ! more than saturation.
 out2 = min (out2, swc3 (p) / 0.58 - soilw3 (p))
 out2 = max (out2, zero)
 ! Update second (middle) layer.
 soilw2 (p) = soilw2 (p) - out2
 ! Update third (bottom) layer.
 soilw3 (p) = soilw3 (p) + out2 - transp3
 ! Pos. of not conserving water if too much for bottom layer,
 ! but probably very minor problem (should test properly).
 soilw3 (p) = max (zero, soilw3 (p))
 ! Water above holding capacity of 3rd layer goes to outflow
 ! (m dt-1).
 outt = max (zero, (soilw3 (p) - swc3 (p)))
 ! Reduce outflow due to drainage impedence.
 outt = outt * fout
 ! Sum outflow for daily use by soil calculations (m day-1).
 outflow (kp) = outflow (kp) + outt
 ! Subtract outflow from soil water in third layer.
 soilw3 (p) = soilw3 (p) - outt
 ! Proportion of water lost as outflow (fraction)
 if ((soilw1 (p) + soilw2 (p) + soilw3 (p) + &
  outt) > eps) then
  ploss = outt / (soilw1 (p) + soilw2 (p) + &
          soilw3 (p) + outt)
 else
  ploss = zero
 end if
 !------------------------------------------------------------------!
 ! Subtract mineral N in outflow from soil mineral N (kg[N] m-2),
 ! allowing for possibility of negative snmin.
 if (snmin (p) > eps) snmin (p) = snmin (p) * (one - ploss)
 !------------------------------------------------------------------!
 ! Update snowpack (m).
 !------------------------------------------------------------------!
 snow (p) = snow (p) + snowi - smelt
 !------------------------------------------------------------------!
 ! Soil water holding capacity.
 !---------------------------------------------------------------------!
 swct (p) = swc1 (p) + swc2 (p) + swc3 (p)
 !---------------------------------------------------------------------!
 ! Soil water potentials (after Johnson et al., 1991) (MPa).
 !---------------------------------------------------------------------!
 if (tsoil > zero) then
  ! Soil water potential top-layer roots (MPa).
  sow1 = soilw1 (p) + soilw2 (p)
  sc1  = swc1   (p) + swc2   (p)
  if (sow1 > eps) then
   if (sc1 > eps) then
    sw1 = swpmax * (1.0 / ((sow1 / sc1) ** bsoil))
   else
    sw1 = -100.0
   end if
  else
   sw1 = -100.0
  end if
  !--------------------------------------------------------------------!
  ! Soil water potential for bottom-layer roots (MPa).
  !--------------------------------------------------------------------!
  if (soilw3 (p) > eps) then
   if (swc3 (p) > eps) then
    sw2 = swpmax * (1.0 / ((soilw3 (p) / &
	  swc3 (p)) ** bsoil))
   else
    sw2 = -100.0
   end if
  else
   sw2 = -100.0
  end if
  !--------------------------------------------------------------------!
  ! Soil water potential for grass (MPa).
  !--------------------------------------------------------------------!
  swp1 (kp) = sw1
  !--------------------------------------------------------------------!
  ! Soil water potential for trees (MPa).
  !--------------------------------------------------------------------!
  if (sw1 > sw2) then
   ! All transpiration from top two layers.
   swp2   (kp) = sw1
   tr_rat (kp) = one
  else
   ! All transpiration from bottom layer.
   swp2   (kp) = sw2
   tr_rat (kp) = zero
  end if
  !--------------------------------------------------------------------!
 else
  !--------------------------------------------------------------------!
  ! Soil frozen.
  !--------------------------------------------------------------------!
  swp1   (kp) = -1.48
  swp2   (kp) = -1.48
  tr_rat (kp) = zero
  !--------------------------------------------------------------------!
 end if ! tsoil > zero
end do ! kp
end subroutine hydrology
