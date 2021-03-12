!======================================================================!
subroutine radiation
!----------------------------------------------------------------------!
use shared
! ---------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------!
integer :: kp
integer :: ki
integer :: p
integer :: k
integer :: ksp
integer :: m
integer :: nlayers
real :: fd
real :: kf
real :: kSWf
real :: frac
real :: fracs
real :: tnfact
real :: nit
real :: prub
real :: pchl
real :: fparb
real :: total
real :: top
real, dimension (mh+1) :: skz
real, dimension (mh+1) :: skSWz
real, dimension (mh+1) :: eskz
real, dimension (mh+1) :: fPARt
real, dimension (mh+1) :: eskSWz
real, dimension (mh+1) :: fSWt
real, dimension (mh+1) :: fSWi
real, dimension (mh+1) :: fPARi
do kp = 1, nplots
 p = p_plot (land_index(i,j),kp)
 if (phenf (kp) == 1) then ! Foliage area changed.
  phenf (kp) = 0
  !---------------------------------------------------------------------!
  ! Factor to calculate absorbed SW by whole plot (ratio).
  !---------------------------------------------------------------------!
  SWf (p) = zero
  !---------------------------------------------------------------------!
  ! Set sums of foliage area * extinction coefficient in each layer
  ! to zero.
  !---------------------------------------------------------------------!
  do m = (mh + 1), 1, -1
   skz   (m) = zero
   skSWz (m) = zero
  end do
  ! Foliage area in each layer.
  do ki = 1, nind (p)
   k = k_ind (land_index(i,j),kp,ki)
   if (alive (k)) then
    ksp = kspp (k)
    ! Set fraction of PAR on top layer absorbed by individual to zero.
    fPAR (k) = zero
    ! Set fraction of SW on top layer absorbed by individual to zero.
    fSW (k) = zero
    do m = (mh + 1), 1, -1
     kz   (m,ki) = zero
     kSWz (m,ki) = zero
    end do
    if (farea (k) > eps) then
     nlayers = height (k) - hbc (k)
     nlayers = max (1, nlayers)
     ! fd is LAI of individual in each layer.
     fd = farea (k) / (float (nlayers) * area)
     ! Factors to make following calculation easier.
     kf   = kpar (ksp) * fd
     kSWf = kSW  (ksp) * fd
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
    end if ! farea > eps
   end if ! alive
  end do ! ki
  !--------------------------------------------------------------------!
  ! Fraction of PAR penetrating layer above top of plot.
  !--------------------------------------------------------------------!
  eskz  (mh+1) = one
  fPARt (mh+1) = one
  !------------------------------------------------------------------!
  ! Fraction of SW radiation penetrating layer above top of plot
  ! (fraction).
  !------------------------------------------------------------------!
  eskSWz (mh+1) = one
  fSWt   (mh+1) = one
  !--------------------------------------------------------------------!
  ! Proportion of light at top of each layer absorbed by that layer.
  !--------------------------------------------------------------------!
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
   do ki = 1, nind (p)
    k = k_ind (land_index(i,j),kp,ki)
    if (alive (k)) then
     ! Save species number for easier coding.
     ksp = kspp (k)
     ! Check if foliage in layer.
     if (skz (m) > eps) then
      ! Fraction of absorption in layer due to individual.
      frac  = kz   (m,ki) / skz   (m)
      fracs = kSWz (m,ki) / skSWz (m)
      ! Sum total fraction absorbed by individual relative to top of
      ! plot.
      ! needed outside radiation?
      fPAR (k) = fPAR (k) + frac * fPARi (m)
      ! Sum total fraction absorbed by individual relative to top of
      ! plot.
      fSW (k) = fSW (k) + fracs * fSWi (m)
      ! Sum factor to calculate absorbed SW by foliage from incident
      ! SW on plot.
      SWf (p) = SWf (p) + (one - rhos (ksp)) * fracs * fSWi (m)
     end if ! skz > eps
    end if ! alive
   end do ! ki
  end do ! m
  !--------------------------------------------------------------------!
  ! Calculate canopy net photosynthesis and conductance parameters.
  !--------------------------------------------------------------------!
  do ki = 1, nind (p)
   k = k_ind (land_index(i,j),kp,ki)
   if (alive (k)) then
    ksp = kspp (k)
    !------------------------------------------------------------------!
    if (farea (k) > eps) then
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
     tnfact = farea (k) * kpar (ksp) / (fPAR (k) * area)
     tnfact = max (one, tnfact)
     ! If foliage in top leaf is too high, reduce foliage N. Assume
     ! extra N just not used for photosynthesis.
     ! IF (tnfact . 2.0) tnfact = 2.0
     nit = tnfact * nfoliage (k) / farea (k)
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
    !------------------------------------------------------------------!
   end if ! alive
  end do! ki
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
  !--------------------------------------------------------------------!
 end if ! phenf (kp) = 1
end do ! kp
!----------------------------------------------------------------------!
end subroutine radiation
!======================================================================!
