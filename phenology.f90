subroutine phenology

use shared
implicit none

! Lower limit of pruba for no frost damage to foliage.
real, parameter :: prubal = 0.13
integer :: kp
integer :: ki
integer :: p
integer :: k,k1,k2
integer :: ksp
real :: wlitterc
real :: wlittern
real :: flitterc
real :: flittern
real :: rlitterc
real :: rlittern
real :: fcneed,rcneed,wcneed,fnneed,rnneed,wnneed
real :: rat_fol,finc,winc,rinc
real :: fgr,wgr,rgr
real :: woff,won
real :: kzog1,kzog2
real :: fPARi_ll
real :: skzo
real :: kSWzog1,kSWzog2
real :: skSWzo
real :: fSWi_phen
real :: fPARiio,fSWiio
real :: tnfact,nit,prub,pchl,top
integer :: ktemp
!----------------------------------------------------------------!
! Accumulate chilling-days (days).
!----------------------------------------------------------------!
      if (tmind < thold) cd (i,j) = cd (i,j) + 1
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
        !cd_flag = 1
       end if ! ddreq (i,j) == zero
      end if ! dd >= ddreq
     
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
      end if !----------------------------------------------------------------!
      if (ktemp >= bgs (i,j)) then
       summer = .TRUE.
      else
       summer = .FALSE.
      end if
!----------------------------------------------------------!
do kp = 1, nplots
 p = p_plot (land_index(i,j),kp)
 do ki = 3, nind (i,j,kp)
  k = k_ind (land_index(i,j),kp,ki)
  if (alive (k) == 1) then
	ksp = kspp (k)
	 ! Tree GPT phenology.
         !----------------------------------------------------------!
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
         flittc (p) = flittc (p) + flitterc
         flittn (p) = flittn (p) + flittern
         wlittc (p) = wlittc (p) + wlitterc
         wlittn (p) = wlittn (p) + wlittern
         rlittc (p) = rlittc (p) + rlitterc
         rlittn (p) = rlittn (p) + rlittern
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
         rgs = rgs + fgr + wgr + rgr
	 if (local) mnppsp (ksp) = mnppsp (ksp) - &
	                           (fgr + wgr + rgr) / dt
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
          farea (k) = rcfoliage (k) * sla (ksp)
         end if
	 ! If this species is cold deciduous, then test for leaf
	 ! area change today.
         if (ptype (ksp) == 3) then
          if (.NOT. summer) then
           farea (k) = zero
          else
           farea (k) = rcfoliage (k) * sla (ksp)
          end if
	  ! If leaves fell off then recalculate fturn for cd
	  ! (leaves may not go on and off in ptype 3 if the signals
	  ! and/or phenology parameter values are not appropriate).
          ! rlold is the actual leaf area on the previous day.
	  ! farea is not the actual between allocate and phen; it
	  ! is the actual between phen and allocate.
          if ((farea (k) == zero) .AND. (rlold (k) > eps) .AND. &
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
          if ((farea (k) > eps) .AND. (rlold (k) == zero)) then
           ! Remove any cold damage.
           nitf (k) = one
           ! Remove any drought damage.
           lasa (k) = lsave (ksp)
          end if
	 ! End of ptype 3 (CD) conditional.
         end if
         !----------------------------------------------------------!
	
	! If this species is potentially dry deciduous then test for
	! leaf area change. Made to lose leaves (once per year) if
	! soil gets below a critical water potential (woff; MPa).
	! Back on when above another value (won; MPa).
        if (ptype (ksp) == 4) then
         woff = -1.49
         won  = -0.5
         ! Set leaf area to previous day's as it may be the first
	 ! day of the year and would have been reset in allocate.
         farea (k) = rlold (k)
         ! Take leaves off if swp falls low enough. Only do once,
	 ! but could make this happen more often if works OK in
	 ! ALLOCATE with fturn>1.
         if (farea (k) > eps) then
          if (swp2 (kp) < woff) then
           if (foff (k) == 0) then
            farea (k) = zero
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
            farea (k) = rcfoliage (k) * sla (ksp)
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
        if (farea (k) /= rlold (k)) then
         ! Flag for tree leaf area change.
         phenf (kp) = 1
        end if ! End of leaf area change conditional.
	! See if leaves are present for drought and cold damage.
        if (farea (k) > (zero + eps)) then
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
	laip (i,j,kp) = laip (i,j,kp) + farea (k)
        !-----------------------------------------------------------!
	! Save foliage area.
        rlold (k) = farea (k)
  end if ! alive
 end do ! ki = 3, nind (i,j,kp)
       !------------------------------------------------------------!
       laip (i,j,kp) = laip (i,j,kp) / area
       !------------------------------------------------------------!
       ! Calculate radiation and physiological variables.
       ! Section from phenology.f.
       !------------------------------------------------------------!
       if (phenf (kp) == 0) then ! Only grass foliage area changed.
        k1 = k_ind (land_index(i,j),kp,1)
	kzog1 = kzg (k1)
        kzg (k1) = kpar (1) * farea (k1) / area
        k2 = k_ind (land_index(i,j),kp,2)
        kzog2 = kzg (k2)
        kzg (k2) = kpar (2) * farea (k2) / area
        skzo = skzg (i,j,kp)
        skzg (i,j,kp) = skzo - (kzog1 + kzog2) + &
                        (kzg (k1) + kzg (k2))
	! New fPAR of lowest layer.
        fPARi_ll = one - exp (-skzg (i,j,kp))
	kSWzog1 = kSWzg (k1)
        kSWzog2 = kSWzg (k2)
        kSWzg (k1) = ksw (1) * farea (k1) / area
        kSWzg (k2) = ksw (2) * farea (k2) / area
        skSWzo = skSWzg (i,j,kp)
        skSWzg (i,j,kp) = skSWzo - (kSWzog1 + kSWzog2) + &
                          (kSWzg (k1) + kSWzg (k2))
	! New fSW of lowest layer.
        fSWi_phen = one - exp (-skSWzg (i,j,kp))
	do ki = 1, nind (i,j,kp)
         k = k_ind (land_index(i,j),kp,ki)
	 if (alive (k) == 1) then
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
         SWf (p) = SWf (p) + fSWtg (i,j,kp) * &
	           (one - rhos (ksp)) * (fSWiig (k) - fSWiio)
	 end if ! alive
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
         tnfact = farea (k) * kpar (ksp) / (fPAR (k) * area)
         tnfact = max (one, tnfact)
	 nit = tnfact * nfoliage (k) / farea (k)
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
	end do ! ki = 1, 2
       !else ! phenf (kp) = 1, so tree foliage changed.
	!--------------------------------------------------------------!
        !phenf (kp) = 0 ! Not sure this needed here.
    !------------------------------------------------------------------!
    ! Factor to calculate absorbed SW by whole plot (ratio).
    !------------------------------------------------------------------!
    !SWf (p) = zero
    end if ! phenf = 1
    end do ! kp
    !------------------------------------------------------------------!
    
end subroutine phenology
