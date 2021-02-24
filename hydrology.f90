subroutine hydrology

use shared
implicit none

integer :: kp
integer :: ki
integer :: p
integer :: k
integer :: ksp
real, parameter :: ftmin = 0.0
real, parameter :: ftmax = 40.0
real, parameter :: fkfive = 18.35
real, parameter :: fkeight = (ftmax - fkfive) / (fkfive - ftmin)
real, parameter :: swpmax = -0.033 ! MPa
real, parameter :: bsoil = 5.0
real, parameter :: intc   = 0.0005 * dt / sday ! m LAI-1 dt-1
real, parameter :: smeltc = 0.0007 * dt / sday ! m K-1 dt-1
real, parameter :: fout   = 0.9 ! fraction
real :: rhoa
real :: lambda
real :: gamma
real :: gammas
real :: ef1
real :: ts
real :: psinf
real :: cwsta
real :: cwstl
real :: cwa
real :: dwo
real :: s
real :: fdq,ft,fdqftfc
real :: tt
real :: sow1
real :: sc1
real :: sw1,sw2
real :: swp1,swp3
real :: tr_grass,tr_tot,tr_rat
real :: psin
real :: cgc
real :: st ! SW radiation absorbed in upper crown layer (W m-2)
real :: fst,fdt
real :: gs
real :: rwi,rw,raws,rsoil
real :: betae
real :: veg_frac
real :: rcif,rci,rcn
real :: cwstf
real :: vpd
real :: ET,Ek
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
cwa = vap / (R * tairK)
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
!----------------------------------------------------------------------!
! Effect of soil temperature on N uptake (m2 kg[C]-1 s-1).
!----------------------------------------------------------------------!
ftsoil = (tsoil * (60.0 - tsoil)) / (800.0 * sday)
ftsoil = max (zero, ftsoil)
ftsoil = min (one , ftsoil)
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
 do ki = 1, nind (i,j,kp)
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
   if (swp1 > (-1.5)) fdt = 1.155 + 0.77 * &
	 (swp1 - 0.015 * float (height (k)))
    if (swp1 >  (-0.2)) fdt = one
    if (swp1 <= (-1.5)) fdt = zero
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
   !rwi = (0.607 * rac (ksp) + one / gs) * area / cfact (k)
   !-------------------------------------------------------------------!
   ! Sum canopy conductance across individuals (m[H2O s-1).
   !-------------------------------------------------------------------!
   !cgc = cgc + one / rwi
   !-------------------------------------------------------------------!
   ! Sum relative amount of transpiration due to trees.
   !-------------------------------------------------------------------!
   !tr_tot = tr_tot + one / rwi
   !if (ksp <= 2) tr_grass = tr_grass + one / rwi
   !-------------------------------------------------------------------!
  end if
  end if ! alive
 end do ! ki
 !---------------------------------------------------------------------!
 ! tr_rat is converted from flag to the fraction of plot
 ! transpiration from top 2 layers (adjusted if trees taking water
 ! from bottom layer, i.e. tr_rat is set from previous day in this
 ! routine).
 !---------------------------------------------------------------------!
 if (tr_rat < 0.1) then
  ! Tree transpiration from bottom (third) layer only.
  if (tr_tot > eps) then
   tr_rat = tr_grass / tr_tot
  else
   tr_rat = 0.5
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
 !---------------------------------------------------------------------!
 ! Vegetation canopy fraction of plot (fraction).
 !---------------------------------------------------------------------!
 veg_frac = 0.6 + 0.0875 * laip (i,j,kp)
 if (laip (i,j,kp) == zero) veg_frac = zero
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
 do ki = 1, nind (i,j,kp)
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
  gc (ki) = one / rcn
  !--------------------------------------------------------------------!
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
 transp = min (transp, (soilw1 (i,j,kp) + soilw2 (i,j,kp) + &
          soilw3 (i,j,kp)))
 !---------------------------------------------------------------------!
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
 !------------------------------------------------------------------!
 ! Subtract mineral N in outflow from soil mineral N (kg[N] m-2),
 ! allowing for possibility of negative snmin.
 if (snmin (i,j,kp) > eps) snmin (i,j,kp) = snmin (i,j,kp) * &
  (one - ploss)
 !------------------------------------------------------------------!
 ! Update snowpack (m).
 !------------------------------------------------------------------!
 snow (i,j,kp) = snow (i,j,kp) + snowi - smelt
 !------------------------------------------------------------------!
 ! Soil water holding capacity.
 !---------------------------------------------------------------------!
 swct (i,j,kp) = swc1 (i,j,kp) + swc2 (i,j,kp) + swc3 (i,j,kp)
 !---------------------------------------------------------------------!
 ! Soil water potentials (after Johnson et al., 1991) (MPa).
 !---------------------------------------------------------------------!
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
  !--------------------------------------------------------------------!
  ! Soil water potential for bottom-layer roots (MPa).
  !--------------------------------------------------------------------!
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
  !--------------------------------------------------------------------!
  ! Soil water potential for grass (MPa).
  !--------------------------------------------------------------------!
  swp1 = sw1
  !--------------------------------------------------------------------!
  ! Soil water potential for trees (MPa).
  !--------------------------------------------------------------------!
  if (sw1 > sw2) then
   ! All transpiration from top two layers.
   swp2 (kp) = sw1
	tr_rat = one
  else
   ! All transpiration from bottom layer.
   swp2 (kp) = sw2
	tr_rat = zero
  end if
  !--------------------------------------------------------------------!
 else
  !--------------------------------------------------------------------!
  ! Soil frozen.
  !--------------------------------------------------------------------!
  swp1 = -1.48
  swp2 (kp) = -1.48
  tr_rat = zero
  !--------------------------------------------------------------------!
 end if ! tsoil > zero
end do ! kp
end subroutine hydrology
