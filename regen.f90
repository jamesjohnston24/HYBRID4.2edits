subroutine regen

use shared
implicit none

integer :: kp
integer :: ki
integer :: ksp
integer :: k
real :: ran
real :: nfp
real :: diamt,diamw,diamh
real :: tarea  ,warea  ,harea
real :: theight,wheight,hheight
real :: hwood,swood,wwood
real :: saparea
real :: cnf
real :: nit
real :: prub,pchl
real :: ht

! Regeneration.
    !------------------------------------------------------------------!
    do kp = 1, nplots
     !****adf for now just fill gaps
     ! Number of each spp. to add to get max. density.
     !iregl = ((limp - (nind (i,j,kp) - 2)) / (nspp - 2))
     !iregl = min (iregl, 10)
     !if (((nind (i,j,kp) - 2) + iregl * (nspp - 2)) > limp) &
     !     iregl = iregl - 1
     !write(*,*)'iregl nind =',iregl,nind(i,j,kp)
     !if (iregl > 0) then
     ksp = 3 !****adf should alter each yr
     do ki = 1, nind (i,j,kp)
      k = k_ind (land_index(i,j),kp,ki)
      if (alive (k) == 0) then ! plant a tree
       alive (k) = 1
       call random_number (ran)
       dbh (k) = idbh * (2.0 * idbhv * ran + (1.0 - idbhv))
       hdbh (k) = zero
       hbc (k) = 0
       nfp = 0.02
       kspp (k) = ksp
       diamt = dbh (k)
       tarea = pi * (0.5 * diamt) ** 2
       diamw = dbh (k) * (one - 2.0 * bark (ksp))
       warea = pi * (0.5 * diamw) ** 2
       diamh = hdbh (k)
       harea = pi * (0.5 * diamh) ** 2
       harea = pi * (0.5 * hdbh (k)) ** 2
       theight = ah (ksp) * diamt ** bh (ksp)
       wheight = ah (ksp) * diamw ** bh (ksp)
       hheight = ah (ksp) * diamh ** bh (ksp)
       cwood (k) = stf (ksp) * formf (ksp) * theight * tarea * &
                   woodd (ksp)
       wwood = stf (ksp) * formf (ksp) * wheight * &
               warea * woodd (ksp)
       hwood = stf (ksp) * formf (ksp) * hheight * &
               harea * woodd (ksp)
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
       farea (k) = zero !****adf set in phen
       !---------------------------------------------------------------!
       ! Fine root C (kg[C] ind-1).
       !---------------------------------------------------------------!
       cfiner (k) = rlratio (ksp) * cfoliage (k)
       !---------------------------------------------------------------!
       lsap (k) = live (ksp) * swood
       cstore (k) = storef (ksp) * lsap (k)
       cwood (k) = cwood (k) - cstore (k)
       !---------------------------------------------------------------!
       ! Foliage N (kg[N] ind-1).
       !---------------------------------------------------------------!
       nfoliage (k) = nfp * 2.0 * cfoliage (k)
       cnf = cfoliage (k) / nfoliage (k)
       nit = nfoliage (k) / (sla (ksp) * cfoliage (k))
       prub = pruba (ksp) + slope (ksp) * nit
       pchl = prub / nrc (ksp)
       gmax (k) = ngr (ksp) * nit * prub
       !---------------------------------------------------------------!
       ! Stem and bark N (kg[N] ind-1).
       !---------------------------------------------------------------!
       nbswood (k) = swood / (cnf / fsr (ksp))
       !---------------------------------------------------------------!
       ! Heartwood N (kg[C] ind-1).
       !---------------------------------------------------------------!
       nheart (k)  = hwood / (cnf / fsr (ksp))
       !---------------------------------------------------------------!
       ! Fine root N (kg[N] ind-1).
       !---------------------------------------------------------------!
       nfiner (k) = cfiner (k) / (cnf / frr (ksp))
       !---------------------------------------------------------------!
       !ngains (k) = zero !****adf
       !---------------------------------------------------------------!
       ! Height of crown (m) allometrically from dbh (m).
       !---------------------------------------------------------------!
       if (dbh (k) > eps) then
        ht = ah (ksp) * dbh (k) ** bh (ksp)
       else
        ht = zero
       end if
       height (k) = nint (ht + 0.5)
       if (height (k) < 1) height (k) = 1
       if (height (k) > mh) height (k) = mh
       !---------------------------------------------------------------!
       ! Frost effect (fraction).
       !---------------------------------------------------------------!
       nitf (k) = one
       !---------------------------------------------------------------!
       ! N storage (kg[N] ind-1).
       !---------------------------------------------------------------!
       navail (k) = zero
       ! Save required C and N tree compartment sizes.
       rcfoliage (k) = cfoliage (k)
       rcwood    (k) = cwood    (k)
       rcfiner   (k) = cfiner   (k)
       rnfoliage (k) = nfoliage (k)
       rnbswood  (k) = nbswood  (k)
       rnfiner   (k) = nfiner   (k)
       ! Initialise previous day's leaf area.
       rlold (k) = farea (k)
       ! Need to conserve C and N in system.
       ! Subtract new individual's C  & N from soil.
       Cpa (i,j,kp) = Cpa (i,j,kp) - (cfoliage (k) + cwood (k) + &
                      cfiner (k) + cstore (k)) / area
       Cpa (i,j,kp) = max (Cpa (i,j,kp), eps)
       Npa (i,j,kp) = Npa (i,j,kp) - (nfoliage (k) + nbswood (k) + &
                      nfiner (k) + navail (k)) / area
       Npa (i,j,kp) = max (Npa (i,j,kp), eps)
       
       !****adf will need to subtract new C and N from somewhere!
       if (ksp == nspp) then
        ksp = 3
       else
        ksp = ksp + 1
       end if
      end if
     end do ! ki
     !****adfphenf (kp) = 1
    end do ! kp
    !------------------------------------------------------------------!
    ! End of regeneration routines.
    !------------------------------------------------------------------!
end subroutine regen
