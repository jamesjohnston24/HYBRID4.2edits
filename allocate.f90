subroutine allocate

use shared
implicit none

integer :: kp
integer :: ki,k1,k2
integer :: k
integer :: ksp
integer :: nlayers,rnlayers
real :: harea,warea,saparea
real :: hheight,wheight
real :: wwood,wood,swood
real :: diamw
real :: hwood,hwoodo,swoodo
real :: ratio
real :: nhinc
real :: winc_min,wincn
real :: cfol_max,cfro_max
real :: finc,rinc,tinc,winc
real :: rat_finc
real :: ZBRENT,TOL,X1,X2,dw
real :: max_la
real :: sneed,sget
real :: ht
real :: fgr,wgr,rgr

!------------------------------------------------------------------!
    ! Tree C and N allocation.
    !------------------------------------------------------------------!
    do kp = 1, nplots
     ! Calculate laip.
     k1 = k_ind (land_index(i,j),kp,1)
     k2 = k_ind (land_index(i,j),kp,2)
     laip (i,j,kp) = farea (k1) + farea (k2)
     do ki = 3, nind (i,j,kp)
      k = k_ind (land_index(i,j),kp,ki)
      if (alive (k) == 1) then
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
      farea (k) = zero
      !----------------------------------------------------------------!
      ! Height of crown (m) allometrically from dbh (m).
      !----------------------------------------------------------------!
      if (dbh (k) > eps) then
       ht = ah (ksp) * dbh (k) ** bh (ksp)
      else
       ht = zero
      end if
      height (k) = nint (ht + 0.5)
      if (height (k) < 1) height (k) = 1
      if (height (k) > mh) height (k) = mh
      !----------------------------------------------------------------!
      ! Height to base of crown (m).
      !----------------------------------------------------------------!
      hbc (k) = min (hbc (k), (height (k) - 1))
      hbc (k) = max (hbc (k), 0)
      !----------------------------------------------------------------!
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
      rgs = rgs + fgr + wgr + rgr
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
      end if ! alive
     end do ! ki = 3, nind (i,j,kp)
     laip (i,j,kp) = laip (i,j,kp) / area
    end do ! kp
    !------------------------------------------------------------------!
    ! End of tree C and N allocation.
    !------------------------------------------------------------------!
    
end subroutine allocate
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
