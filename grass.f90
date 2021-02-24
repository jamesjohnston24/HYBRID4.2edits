subroutine grass

use shared
implicit none

integer :: kp
integer :: ki
integer :: p
integer :: k
integer :: ksp
real :: cfmin,maxcf
real :: flitterc,wlitterc,rlitterc
real :: flittern,wlittern,rlittern
real :: cso,cto
real :: cavail,potf,max_folg,potfg
real :: csn,ctn
real :: grg
real :: rnmax,rnlose,rnsave

do kp = 1, nplots
       p = p_plot (land_index(i,j),kp)
       laip (i,j,kp) = zero
       do ki = 1, nind (i,j,kp)
        k = k_ind (land_index(i,j),kp,ki)
	if (alive (k) == 1) then
	ksp = kspp (k)
        !-----------------------------------------------------------!
	if (ksp <= 2) then
         !----------------------------------------------------------!
	 ! Min. amount of grass foliage mass allowed (kg[C] plot-1).
         !----------------------------------------------------------!
	 cfmin = dbh (k) * area / sla (ksp)
         !----------------------------------------------------------!
	 ! Reduce maximum leaf area if ball less than zero.
	 ! Could put in temperature limitation to growth here.
	 ! Also, could make this reduction proportional to LAI.
         !----------------------------------------------------------!
         if (ball (k) < zero) then
          maxcf = 0.9 * cfoliage (k)
         else
          maxcf = cfoliage (k) * 10.0
         end if
         maxcf = max (cfmin, maxcf)
         !----------------------------------------------------------!
	 ! Grass daily C and N  litter (kg day-1).
         !----------------------------------------------------------!
         flitterc = fturn_plot (kp,ksp) * cfoliage (k)
         rlitterc = rturn (ksp) * cfiner (k)
         wlitterc = wturn (ksp) * lsap (k)
         flittern = frcoeff (ksp) * fturn_plot (kp,ksp) * &
                    nfoliage (k)
         rlittern = rrcoeff (ksp) * rturn (ksp) * nfiner (k)
         wlittern = wturn (ksp) * nbswood (k)
         !----------------------------------------------------------!
	 ! Subtract daily grass C and N litter from tissue.
         !----------------------------------------------------------!
         cfoliage (k) = cfoliage (k) - flitterc
         cfiner   (k) = cfiner   (k) - rlitterc
         lsap     (k) = lsap     (k) - wlitterc
         nfoliage (k) = nfoliage (k) - flittern
         nfiner   (k) = nfiner   (k) - rlittern
         nbswood  (k) = nbswood  (k) - wlittern
         !----------------------------------------------------------!
	 ! Increment litter pools (kg plot-1).
         !----------------------------------------------------------!
	 wlittc (p) = wlittc (p) + wlitterc
         wlittn (p) = wlittn (p) + wlittern
         flittc (p) = flittc (p) + flitterc
         flittn (p) = flittn (p) + flittern
         rlittc (p) = rlittc (p) + rlitterc
         rlittn (p) = rlittn (p) + rlittern
         !----------------------------------------------------------!
	 ! Initial grass structural C (kg[C] plot-1). Used to calculate
	 ! growth respiration.
         !----------------------------------------------------------!
         cso = cfoliage (k) + lsap (k) + cfiner (k)
         !----------------------------------------------------------!
	 ! Initial total grass C (kg[C] plot-1). Used to calculate C
	 ! balance.
         !----------------------------------------------------------!
         cto = cso + cstore (k)
         !----------------------------------------------------------!
	 ! Total grass C available for partitioning (kg[C] plot-1).
         !----------------------------------------------------------!
         cavail = cfoliage (k) + lsap (k) + cfiner (k) + cstore (k)
         cavail = max (zero, cavail)
         !----------------------------------------------------------!
	 ! Maximum possible grass foliage C based on available C
	 ! (kg[C] plot-1).
         !----------------------------------------------------------!
         potf = cavail / (one + (one + bark (ksp)) * rlratio (ksp))
         !----------------------------------------------------------!
	 ! Grass foliage growth limited by soil water potential.
	 !****adf swp2 should be mean over day.
         !----------------------------------------------------------!
         max_folg = 0.1 * cfoliage (k) * (swp2 (kp) + 1.5) / 1.5
         max_folg = max (zero, max_folg)
         potfg    = cfoliage (k) + max_folg
         !----------------------------------------------------------!
	 ! Actual foliage C lower than potential if ball was -ve.
         !----------------------------------------------------------!
         cfoliage (k) = min (maxcf, potf)
         cfoliage (k) = min (cfoliage (k), potfg)
         !----------------------------------------------------------!
         ! Do not allow foliage C to fall below a minimum (set by
	 ! 'dbh').
         !----------------------------------------------------------!
         cfoliage (k) = max (cfoliage (k), cfmin)
         !----------------------------------------------------------!
	 lsap (k) = bark (ksp) * cfoliage (k)
         cfiner (k) = rlratio (ksp) * cfoliage (k)
         !----------------------------------------------------------!
	 ! New grass structural C (kg[C] plot-1). Used to calculate
	 ! growth respiration and store.
         !----------------------------------------------------------!
         csn = cfoliage (k) + lsap (k) + cfiner (k)
         cstore (k) = cavail - csn
         !----------------------------------------------------------!
	 ! cstore may go negative if cfoliage set to minimum, hence:
         !----------------------------------------------------------!
         cstore (k) = max (cstore (k), zero)
         !----------------------------------------------------------!
	 ! Grass growth respiration, if positive growth.
         !----------------------------------------------------------!
         grg = (one - rgf (ksp)) * max (zero, (csn - cso))
         grg = min (lsap (k), grg)
         rgs = rgs + grg
	 if (local) mnppsp (ksp) = mnppsp (ksp) - grg / dt
         !----------------------------------------------------------!
	 farea (k) = cfoliage (k) * sla (ksp)
         !----------------------------------------------------------!
	 ! The following allows for growth respiration, but is really
	 ! a cheat as should maintain correct relative proportions.
         !----------------------------------------------------------!
         lsap (k) = lsap (k) - grg
         cwood (k) = lsap (k) + cstore (k)
         !----------------------------------------------------------!
	 ! Grass N available for partitioning (kg[N] plot-1).
         !----------------------------------------------------------!
         navail (k) = navail (k) + nfoliage (k) + &
                      nfiner (k) + nbswood (k)
         rnmax = 0.1 * csn
         rnlose = max (zero, (navail (k) - rnmax))
         navail (k) = navail (k) - rnlose
         rnsave = rnlose
         !----------------------------------------------------------!
	 ! Grass N allocation.
         !----------------------------------------------------------!
         nfoliage (k) = navail (k) * cfoliage (k) / (cfoliage (k) + &
                        fsr (ksp) * lsap (k) + &
                        frr (ksp) * cfiner (k))
         nbswood (k) = navail (k) * fsr (ksp) * lsap (k) / &
	               (cfoliage (k) + fsr (ksp) * lsap (k) + &
                       frr (ksp) * cfiner (k))
         nfiner (k) = navail (k) * frr (ksp) * cfiner (k) / &
                      (cfoliage (k) + fsr (ksp) * lsap (k) + &
                      frr (ksp) * cfiner (k))
         navail (k) = rnsave
         nitf (k) = one
         !----------------------------------------------------------!
	 ! Set C balance of lowest 10 % of grass canopy to zero at
	 ! the end of each day.
         !----------------------------------------------------------! 
         ball (k) = zero
         !----------------------------------------------------------!
	 ! New grass structural C (kg[C] plot=1). Used to calculate
	 ! growth respiration.
         !----------------------------------------------------------!
         ctn = cfoliage (k) + lsap (k) + cfiner (k) + cstore (k) + &
	       grg
         if (cto /= ctn) then
          Cpa (i,j,kp) = Cpa (i,j,kp) - (ctn - cto) / area
          Cpa (i,j,kp) = max (Cpa (i,j,kp), eps)
         end if
         !----------------------------------------------------------!
	end if
	end if ! alive
       end do ! ki
      end do ! kp
      
end subroutine grass
