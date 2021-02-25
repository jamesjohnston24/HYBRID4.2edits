subroutine carbon

use shared
implicit none

integer :: kp
integer :: ki
integer :: k
integer :: ksp
real :: kci
real :: koi
real :: klco
real :: klcc
real :: c,d
real :: jlcmaxn
real :: v,Y
real :: par
real :: tfac
real :: num,den
real :: kt
real :: cair
real :: vmax
real :: rch
real :: A,E
real :: efact
real :: vmaxo
real :: cffol,cfwood,cffiner,cflux
real :: rd
real :: ci
real :: ai
real :: ac,ae
real :: jmaxn,jn
real :: vcmax,vomax,pcpi,X,Z,alc,bcarb,brubp,ccarb,crubp
real :: dum,q,ccicarb,ccirubp,acarb,arubp

do kp = 1, nplots
 ! Kinetic parameters in air space equivalents (mol[CO2] m-3).
 kci = 1.9925e15 * exp (- 10127.01  / tfolK (kp)) / tfolK (kp)
 koi = 2.4239e6  * exp (- 1828.1536 / tfolK (kp)) / tfolK (kp)
 ! Rubisco oxygenation turnover number (klco) is temperature
 ! dependent.
 klco = 4.397e07 * exp (- 5292.02 / tfolK (kp))
 ! Rubisco carboxylation turnover number (klcc) is temperature
 ! dependent.
 klcc = 2.897e14 * exp (- 9862.41 / tfolK (kp))
 ! Light saturated rate of electron transport (mol mol(chl)-1 s-1).
 c = 3.486e13 * exp (- 9561.7 / tfolK (kp))
 d = one + exp (78.178 - 23934.4 / tfolK (kp))
 jlcmaxn = c / d
 ! Factors for Ci solution in PGEN.
 v = kci * (one + 8.471 / koi)
 Y = 1.231e-7 * tfolK (kp)
 do ki = 1, nind (i,j,kp)
  !---------------------------------------------------------------!
  ! Get index of individual.
  !---------------------------------------------------------------!
  k = k_ind (land_index(i,j),kp,ki)
  if (alive (k) == 1) then
   ksp = kspp (k)
   if ((cfoliage (k) > eps) .and. (et_cat (k) > eps)) then
    if (radpar > zero) then
     par = radpar * ipfact (k)
     !------------------------------------------------------------!
     ! efact allows for effect of transpiration on movement of CO2
     ! molecules.
     efact = 2.462e-7 * tfolK (kp)
     if (ksp == 2) then ! C4 (uses Collatz et al., 1992).
      ! Vcmax at 25 oC (mol m-2 s-1).
      vmaxo = 1.248 * et_cat (k)
      ! Q10 2 response.
      tfac = 2 ** ((tfolK (kp) - 298.15) / 10.0)
      num = vmaxo * tfac
      den = (one + exp (0.3 * (286.15 - tfolK (kp)))) * &
            (one + exp (0.3 * (tfolK (kp) - 309.15)))
      vmax = num / den
	  kt = 18000.0 * vmax
	  cair = P_Pa / (R * tfolK (kp))
          rch = one / gc (k)
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
          alc = gc (k) + Y
	  bcarb = gc (k) * (v - ccair) + Y * (v + ccair) + vcmax - rd
          brubp = gc (k) * (Z - ccair) + Y * (Z + ccair) + X - rd
          ccarb = v * ccair * (y - gc (k)) - vcmax * pcpi - rd * v
          crubp = Z * ccair * (Y - gc (k)) - pcpi * (X + 2.3333 * rd)
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
          acarb = gc (k) * (ccair - ccicarb) - &
	          ((ccair + ccicarb) / 2.0) * efact
          ! RuBP-regeneration-limited net photosynthesis.
          arubp = gc (k) * (ccair - ccirubp) - &
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
       ! Grid-box NPP (kg[C] s-1).
       !---------------------------------------------------------------!
       NPP_grid (i,j) = NPP_grid (i,j) + cflux
       !---------------------------------------------------------------!
       if (local) then
        mnppsp (ksp) = mnppsp (ksp) + cflux
       end if
       end if ! alive
      end do ! ki
end do ! kp
end subroutine carbon
