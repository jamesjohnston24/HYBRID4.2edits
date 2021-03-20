subroutine nitrogen

use shared
implicit none

integer :: kp,p
integer :: ki
integer :: ksp
integer :: k
real :: ftsoil ! N uptake temperature effect            (fraction)
real :: pptod
real :: tup
real :: harea
real :: hheight
real :: hwood
real :: totalC,totalN
real :: rnfrac
real :: rat
real, dimension (nind_max) :: ngains
real, dimension (nspp) :: fwsoil

pptod = sum (pre (i, j, it-3:it)) / 1000.0

      !-------------------------------------------------------------!
      ! Mineral N addition from atmosphere (kg[C] m-1 d-1).
      !-------------------------------------------------------------!
      Ndepo = nf + 4.0e-3 * pptod
      !-------------------------------------------------------------!
!----------------------------------------------------------------!

!----------------------------------------------------------------------!
! Effect of soil temperature on N uptake (m2 kg[C]-1 s-1).
!----------------------------------------------------------------------!
ftsoil = (tsoil * (60.0 - tsoil)) / (800.0 * sday)
ftsoil = max (zero, ftsoil)
ftsoil = min (one , ftsoil)
!----------------------------------------------------------------------!

do kp = 1, nplots
 p = p_plot (land_index(i,j),kp)
  !---------------------------------------------------------------------!
 ! Convert soil water to water-filled pore space.
 ! Assumes micro-pore space = swc and macro-pore space = 42%
 ! saturation content (from TEM, for loam; Raich et al., 1991).
 !---------------------------------------------------------------------!
 if (swct (p) > eps) then
  wfps (kp) = 100.0 * (soilw1 (p) + soilw2 (p) + &
    soilw3 (p)) / (1.7241 * swct (p))
  wfps (kp) = min (100.0, wfps (kp))
 else
  wfps (kp) = 0.00001
 end if
      ! Effect of soil water on N uptake through anaerobic conditions at
      ! high soil water and lack of solution medium at low. Function
      ! for decomposition used. weff is GPT parameter for effect of
      ! soil water saturation on N uptake.
      !----------------------------------------------------------------!
      if (snmin (p) > eps) then
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
       !----------------------------------------------------------------!
      ! Total N uptake rate in plot (kg[N] s-1).
      !----------------------------------------------------------------!
      tup = zero
      !----------------------------------------------------------------!
      do ki = 1, nind (p)
       !---------------------------------------------------------------!
       ! Get index of individual.
       !---------------------------------------------------------------!
       k = k_ind (land_index(i,j),kp,ki)
       if (alive (k) == 1) then
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
        ngains (ki) = cfiner (k) * nupc (ksp) * snmin (p) * &
	             fwsoil (ksp) * ftsoil / rnfrac
        !--------------------------------------------------------------!
        ! For grass, only N available from top layer.
        !--------------------------------------------------------------!
        if (ksp <= 2) then
	 if (swct (p) > eps) then
          ngains (ki) = ngains (ki) * &
		       (swc1 (p) + swc2 (p)) / swct (p)
         end if
        end if
       else
        ngains (ki) = zero
       end if
       !---------------------------------------------------------------!
       ! Total N uptake rate in plot (kg[N] s-1).
       !---------------------------------------------------------------!
       tup = tup + ngains (ki)
       end if ! alive
       end do ! ki
       !----------------------------------------------------------------!
      if ((sday * tup) > eps) then
       rat = area * snmin (p) / (sday * tup)
      else
       rat = zero
      end if
      !----------------------------------------------------------------!
      ! Limit N uptake if more than available in plot.
      !----------------------------------------------------------------!
      if (rat < one) then
       do ki = 1, nind (p)
       k = k_ind (land_index(i,j),kp,ki)
       if (alive (k) == 1) then
        ngains (ki) = rat * ngains (ki)
       end if ! alive
       end do ! ki
       tup = rat * tup
      end if
      !----------------------------------------------------------------!
      ! Subtract plant N uptake from soil (kg[N] m-2).
      !----------------------------------------------------------------!
      snmin (p) = snmin (p) - sday * tup / area
      !snmin (p) = 0.004
      !----------------------------------------------------------------!
      do ki = 1, nind (p)
       !---------------------------------------------------------------!
       k = k_ind (land_index(i,j),kp,ki)
       if (alive (k) == 1) then
       !---------------------------------------------------------------!
       ! Update individual N store (kg[N] ind-1).
       !---------------------------------------------------------------!
       navail (k) = navail (k) + sday * ngains (ki)
       !---------------------------------------------------------------!
       end if ! alive
      end do ! ki
      end if ! snmin > eps
end do ! kp
end subroutine nitrogen
