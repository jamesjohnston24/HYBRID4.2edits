subroutine annual_diagnostics

use shared
implicit none

integer :: kp
integer :: ki
integer :: k
integer :: ksp
integer, dimension (nspp) :: ntreessp
real   , dimension (nspp) :: mlaisp
real   , dimension (nspp) :: mheightsp

!----------------------------------------------------------------------!
do j = j1, j2
 do i = i1, i2
  !--------------------------------------------------------------------!
  ! m2 m-2.
  !--------------------------------------------------------------------!
  LAI_grid (i,j) = sum (laip (i,j,:)) / float (nplots)
  !--------------------------------------------------------------------!
  ! kg[C] m-2 yr-1
  !--------------------------------------------------------------------!
  NPP_grid (i,j) = (one - icwtr (i,j)) * (dt * NPP_grid (i,j) - rgs) / &
                   (float (nplots) * area)
  !--------------------------------------------------------------------!
  Cv_grid      (i,j) = zero
  Cs_grid      (i,j) = zero
  Cv_C3GR_grid (i,j) = zero
  Cv_C4GR_grid (i,j) = zero
  Cv_BRCD_grid (i,j) = zero
  Cv_NLEV_grid (i,j) = zero
  !--------------------------------------------------------------------!
  do kp = 1, nplots
   Cs_grid (i,j) = Cs_grid (i,j) + soilC (i,j,kp)
   do ki = 1, nind (i,j,kp)
    k = k_ind (land_index(i,j),kp,ki)
    if (alive (k) == 1) then
     Cv_grid (i,j) = Cv_grid (i,j) + cfoliage (k) + cwood (k) + &
	             cfiner (k)
     ksp = kspp (k)
     select case (ksp)
      Case (1)
       Cv_C3GR_grid (i,j) = Cv_C3GR_grid (i,j) + &
	                    cfoliage (k) + cwood (k) + cfiner (k)
      Case (2)
       Cv_C4GR_grid (i,j) = Cv_C4GR_grid (i,j) + &
	                    cfoliage (k) + cwood (k) + cfiner (k)
      Case (3)
       Cv_BRCD_grid (i,j) = Cv_BRCD_grid (i,j) + &
	                    cfoliage (k) + cwood (k) + cfiner (k)
      Case (4)
       Cv_NLEV_grid (i,j) = Cv_NLEV_grid (i,j) + &
	                    cfoliage (k) + cwood (k) + cfiner (k)
     end select
    end if ! alive
   end do ! ki
  end do ! kp
  !--------------------------------------------------------------------!
  Cv_grid (i,j) = (one - icwtr (i,j)) * Cv_grid (i,j) / &
                  (float (nplots) * area)
  Cs_grid (i,j) = (one - icwtr (i,j)) * Cs_grid (i,j) / float (nplots)
  Cv_C3GR_grid (i,j) = (one - icwtr (i,j)) * Cv_C3GR_grid (i,j) / &
                       (float (nplots) * area)
  Cv_C4GR_grid (i,j) = (one - icwtr (i,j)) * Cv_C4GR_grid (i,j) / &
                       (float (nplots) * area)
  Cv_BRCD_grid (i,j) = (one - icwtr (i,j)) * Cv_BRCD_grid (i,j) / &
                       (float (nplots) * area)
  Cv_NLEV_grid (i,j) = (one - icwtr (i,j)) * Cv_NLEV_grid (i,j) / &
                       (float (nplots) * area)
  !--------------------------------------------------------------------!
  if (.NOT. local) then
   ! Pg[C] yr-1.
   NPP_global = NPP_global + larea (i,j) * NPP_grid (i,j) * &
                1.0e6 / 1.0e12
   if(isNaN(NPP_global))then
    write (*,*)i,j,NPP_grid(i,j)
    write (*,*)i,j,lon(i),lat(j)
    stop
   end if
   Cv_global = Cv_global + larea (i,j) * Cv_grid (i,j) * &
               1.0e6 / 1.0e12
   Cs_global = Cs_global + larea (i,j) * Cs_grid (i,j) * &
               1.0e6 / 1.0e12
  end if ! .NOT. local
  !--------------------------------------------------------------------!
 end do ! i1, i2
end do ! j1, j2
    
!----------------------------------------------------------------------!
if (local) then
 write (*,*) Cv_C3GR_grid(i,j),Cv_C3GR_grid(i1,j1)
 write (*,*) 'kyr NPP_grid      ',kyr,NPP_grid(i1,j1),'kg[C] m-2 yr-1'
 write (*,*) 'kyr Cv_grid       ',kyr,Cv_grid (i1,j1),'kg[C] m-2'
 write (*,*) 'kyr Cs_grid       ',kyr,Cs_grid (i1,j1),'kg[C] m-2'
 write (*,*) 'kyr Cv_C3GR_grid  ',kyr,Cv_C3GR_grid (i1,j1),'kg[C] m-2'
 write (*,*) 'kyr Cv_C4GR_grid  ',kyr,Cv_C4GR_grid (i1,j1),'kg[C] m-2'
 write (*,*) 'kyr Cv_BRCD_grid  ',kyr,Cv_BRCD_grid (i1,j1),'kg[C] m-2'
 write (*,*) 'kyr Cv_NLEV_grid  ',kyr,Cv_NLEV_grid (i1,j1),'kg[C] m-2'
 write (*,*) 'kyr LAI_grid      ',kyr,LAI_grid(i1,j1),'m2 m-2'
 ntreessp (:) = 0
 mlaisp (:) = zero
 mheightsp (:) = zero
 do kp = 1, nplots
  do ki = 1, nind (i1,j1,kp)
   k = k_ind (land_index(i1,j1),kp,ki)
   if (alive (k) == 1) then
    ksp = kspp (k)
    ntreessp  (ksp) = ntreessp (ksp) + 1
    mlaisp    (ksp) = mlaisp (ksp) + cfoliage (k) * sla (ksp)
    mheightsp (ksp) = mheightsp (ksp) + height (k)
   end if ! alive
  end do ! ki
 end do ! kp
 write (*,*) '        ksp   ntrees/plot      LAI         height&
 &           NPP'
 do ksp = 1, nspp
  if (ntreessp (ksp) > 0) then
   write (*,*) ksp,ntreessp(ksp)/float(nplots),&
               mlaisp(ksp)/(area*float(nplots)), &
 	       mheightsp(ksp)/float(ntreessp(ksp)), &
	       dt*mnppsp(ksp)/(area*float(nplots))
  else
   write (*,*) ksp,ntreessp(ksp)/float(nplots),&
               mlaisp(ksp)/(area*float(nplots)), &
 	       0.0,dt*mnppsp(ksp)/(area*float(nplots))
  end if
 end do
else
 write (20,*) kyr,GPP_global,NPP_global,Cv_global,Cs_global
end if
!----------------------------------------------------------------------!
  
end subroutine annual_diagnostics
