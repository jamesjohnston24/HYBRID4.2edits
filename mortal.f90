subroutine mortal

use shared
implicit none

integer :: kp
integer :: ki
integer :: p
integer :: k
integer :: in
real :: carbonm,carbon1,carbon2
real :: prob
real :: ran

! Mortality routines (kill individual trees if low growth.
    !------------------------------------------------------------------!
    do kp = 1, nplots
     p = p_plot (land_index(i,j),kp)
     k = k_ind (land_index(i,j),kp,1)
     alive (k) = 1
     k = k_ind (land_index(i,j),kp,2)
     alive (k) = 1
     ! Set number of live individuals in plot to two (for grass).
     in = 2
     do ki = 3, nind (i,j,kp)
      k = k_ind (land_index(i,j),kp,ki)
      if (alive(k) == 1) then
       !ksp = kspp (k)
       ! Sum test value for tree death, fine roots used since evergreen
       ! trees may hang onto leaves after lost all fine roots through
       ! turnover, foliage for inverse reason.
       carbon1 = cfiner   (k) + cstore (k)
       carbon2 = cfoliage (k) + cstore (k)
       carbonm = min (carbon1, carbon2)
       ! Stochastic mortality, probability of mortality is prob.
       prob = one / 400.0
       call random_number (ran)
       ! Test for stochastic mortality.
       if (ran > (one - prob)) carbonm = zero
       ! Individual is dead.
       if (carbonm <= 0.00001) then
        ! Add C and N to plot litter.
        wlittc (p) = wlittc (p) + cwood (k)
        wlittn (p) = wlittn (p) + nbswood (k) + nheart (k)
        flittn (p) = flittn (p) + nfoliage (k)
        rlittn (p) = rlittn (p) + nfiner   (k)
        flittc (p) = flittc (p) + cfoliage (k) + cstore (k)
        rlittc (p) = rlittc (p) + cfiner   (k)
        alive (k) = 0
        ! Go to next individual
        GOTO 300
       end if
       ! Alive, so add 1 to number alive in plot.
       alive (k) = 1
       !in = in + 1
   300 CONTINUE
      end if
     end do ! ki
     ! Number of individuals left in plot.
     !nind (i,j,kp) = in
     !write(*,*)nind(i,j,kp)
    end do ! kp
    !------------------------------------------------------------------!
    ! End of mortality routines.
    !------------------------------------------------------------------!
end subroutine mortal
