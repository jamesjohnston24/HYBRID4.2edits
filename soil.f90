subroutine soil

use shared
implicit none

real, parameter :: texture = 0.5 ! Soil texture (units?)
real   , parameter :: Tc_local  = 0.3 ! Clay fraction?
real   , parameter :: Ts_local  = 0.3 ! Sand fraction?
real   , parameter :: plf       = 0.01 + 0.04 * Ts_local
! Lignin to biomass ratio in leaf litter.
real, parameter :: lfl = 0.20
! Lignin to biomass ratio in root litter.
real, parameter :: lrl = 0.16
! Biomass C content (fraction).
real, parameter :: w = 0.45
real, parameter :: psab = 0.486
! Decay rates base values from Comins & McMurtrie (1993), but
! doubled to allow for inclusion of the soil water factor (/day).
real, parameter :: d1pb = 2.0 * 0.076 * exp (-3.0 * lfl)      / 7.0
real, parameter :: d2pb = 2.0 * 0.28                          / 7.0
real, parameter :: d3pb = 2.0 * 0.094 * exp (-3.0 * lrl)      / 7.0
real, parameter :: d4pb = 2.0 * 0.35                          / 7.0
real, parameter :: d5pb = 2.0 * 0.14 * (one - 0.75 * texture) / 7.0
real, parameter :: d6pb = 2.0 * 0.0038                        / 7.0
real, parameter :: d7pb = 2.0 * 0.00013                       / 7.0
! Partition coefficients.
real, parameter :: pau = 0.55 * (one - lfl)
real, parameter :: psu = 0.7 * lfl
real, parameter :: pav = 0.45 * (one - lfl)
real, parameter :: psv = 0.7 * lfl
real, parameter :: pam = 0.45
real, parameter :: pan = 0.45
real, parameter :: ppa = 0.004
real, parameter :: pas = 0.42
real, parameter :: pps = 0.03
real, parameter :: pap = 0.45
real, parameter :: emf = 0.05 ! N emission (fraction)
real, parameter :: Nad = 0.0 ! N addition (kg[N] m-2 day-1)
real, parameter :: nfx = 0.0 ! N fixation (kg[N] m-2 day-1)
integer :: kp
integer :: p
real :: et_soil ! Soil decomposition temperature modifier (fraction)
real :: em_soil ! Soil decomposition water modifer (fraction)
real :: ev      ! Overall soil decomposition modifier (fraction)
real :: h2o_30
real :: pleach
real :: pmf,puf,pnr,pvr,psa
real :: d1p,d2p,d3p,d4p,d5p,d6p,d7p
real :: d1C,d2C,d3C,d4C,d5C,d6C,d7C
real :: dCu,dCm,dCv,dCn,dCa,dCs,dCpa,dCle
real :: dNu,dNm,dNv,dNn,dNa,dNs,dNpa
real :: Ju,Jv,Ja,Js,Jpa
real :: NmfNmw,Nnr
real :: temp
real :: C0,N0
real :: Cleach
real :: Nmin
real :: Usoil

! Soil routine.
      ! et_soil is cited in F97 as from Comins & McMurtrie (1993).
      ! Equation from HYBRID4 code.
      if (tsoil > eps) then
       et_soil = 0.0326 + 0.00351 * tsoil ** 1.652 - &
                (0.023953 * tsoil) ** 7.19
      else
       et_soil = 0.0326
      end if
      et_soil = max (zero, et_soil)
      et_soil = min (one , et_soil)
      do kp = 1, nplots
       p = p_plot (land_index(i,j),kp)
       wlittc (p) = wlittc (p) / area
       wlittn (p) = wlittn (p) / area
       flittc (p) = flittc (p) / area
       flittn (p) = flittn (p) / area
       rlittc (p) = rlittc (p) / area
       rlittn (p) = rlittn (p) / area
 !---------------------------------------------------------------------!
 ! Convert soil water to water-filled pore space.
 ! Assumes micro-pore space = swc and macro-pore space = 42%
 ! saturation content (from TEM, for loam; Raich et al., 1991).
 !---------------------------------------------------------------------!
 ! Moved to nitrogen.f90 as needed there first.
 !if (swct (p) > eps) then
 ! wfps (kp) = 100.0 * (soilw1 (p) + soilw2 (p) + &
 !   soilw3 (p)) / (1.7241 * swct (p))
 ! wfps (kp) = min (100.0, wfps (kp))
 !else
 ! wfps (kp) = 0.00001
 !end if
       !---------------------------------------------------------------!
       ! Soil-water decomposition modifer (fraction). Equations are 
       ! fitted to Fig. 8 of Williams et al. (1992).
       ! wfps should really be mean over day.
       !---------------------------------------------------------------!
       if (wfps (kp) < 60.0) then
        em_soil = exp ((wfps (kp) - 60.0) ** 2 / (-800.0))
       else
        em_soil = 0.000371 * wfps (kp) ** 2 - 0.0748 * wfps (kp) + 4.13
       endif
       !---------------------------------------------------------------!
       em_soil = max (zero, em_soil)
       em_soil = min (one , em_soil)
       !---------------------------------------------------------------!
       ! Overall decomposition modifier.
       !---------------------------------------------------------------!
       ev = et_soil * em_soil
       !---------------------------------------------------------------!
       ! Soil microbe leach fraction.
       !---------------------------------------------------------------!
       h2o_30 = outflow (kp)
       outflow (kp) = zero
       pleach = 0.05555 * h2o_30 * plf
       ! Partitioning coefficients.
       if (flittn (p) > eps) then
        pmf = 0.85 - 0.018 * lfl / (w * flittn (p) / flittc (p))
	pmf = max (zero, pmf)
       else
        pmf = zero
       end if
       puf = one - pmf
       if (rlittn (p) > eps) then
        pnr = 0.85 - 0.018 * lrl / (w * rlittn (p) / rlittc (p))
	pnr = max (zero, pnr)
       else
        pnr = zero
       end if
       pvr = one - pnr
       psa = psab - pleach
       ! Decay rates (/day).
       d1p = ev * d1pb
       d2p = ev * d2pb
       d3p = ev * d3pb
       d4p = ev * d4pb
       d5p = ev * d5pb
       d6p = ev * d6pb
       d7p = ev * d7pb
       ! C decays (kg[C] m-2 day-1).
       d1C = d1p * Cu  (p)
       d2C = d2p * Cm  (kp)
       d3C = d3p * Cv  (p)
       d4C = d4p * Cn  (p)
       d5C = d5p * Ca  (p)
       d6C = d6p * Cs  (p)
       d7C = d7p * Cpa (p)
       ! C fluxes between litter pools.
       dCu = puf * flittc (p) + wlittc (p) - d1C
       dCm = pmf * flittc (p)               - d2C
       dCv = pvr * rlittc (p)               - d3C
       dCn = pnr * rlittc (p)               - d4C
       ! C fluxes to soil pools.
       Ja  = pau * d1C + pam * d2C + pav * d3C + pan * d4C + &
             pas * d6C + pap * d7C
       Js  = psu * d1C + psv * d3C + psa * d5C
       Jpa = ppa * d5C + pps * d6C
       ! C fluxes between soil pools.
       dCa  = Ja  - d5C
       dCs  = Js  - d6C
       dCpa = Jpa - d7C
       dCle = pleach * d5C
       ! C fluxes to structural litter pools.
       Ju = puf * flittc (p) + wlittc (p)
       Jv = pvr * rlittc (p)
       ! N flux from surface litter to surface metabolic.
       NmfNmw = flittn (p) + wlittn (p) - vu * Ju
       ! N flux from root litter to soil metabolic.
       Nnr = rlittn (p) - vv * Jv
       ! N fluxes between litter pools.
       dNu = vu * Ju - d1p * Nu (p)
       dNm = NmfNmw  - d2p * Nm (p)
       dNv = vv * Jv - d3p * Nv (p)
       dNn = Nnr     - d4p * Nn (p)
       ! Constrained forms, if required.
       if ((Cm (p) + dCm) > eps) then
        temp = (Nm (p) + dNm) / (Cm (p) + dCm)
       else
	temp = one
       end if
       if (temp < (one / 25.0)) then
        vm = one / 25.0
        NmfNmw = vm * pmf * flittc (p)
        dNm    = NmfNmw - d2p * Nm (p)
       end if
       if (temp > (one / 10.0)) then
        vm = one / 10.0
        NmfNmw = vm * pmf * flittc (p)
        dNm    = NmfNmw - d2p * Nm (p)
       end if
       if ((Cn (p) + dCn) > eps) then
        temp = (Nn (p) + dNn) / (Cn (p) + dCn)
       else
        temp = one
       end if
       if (temp < (one / 25.0)) then
        vn = one / 25.0
        Nnr    = vn * pnr * rlittc (p)
        dNn    = Nnr    - d4p * Nn (p)
       end if
       if (temp > (one / 10.0)) then
        vn = one / 10.0
        Nnr    = vn * pnr * rlittc (p)
        dNn    = Nnr    - d4p * Nn (p)
       end if
       ! N fluxes between soil pools.
       dNa  = va  * Ja  - d5p * Na  (p)
       dNs  = vs  * Js  - d6p * Ns  (p)
       dNpa = vpa * Jpa - d7p * Npa (p)
       ! Inital total C.
       C0 = flittc (p) + wlittc (p) + rlittc (p) + &
            Cu (p) + Cm (p) + Cv (p) + Cn (p) + &
            Ca (p) + Cs (p) + Cpa (p)
       ! Initial total N.
       N0 = flittn (p) + wlittn (p) + rlittn (p) + &
            Nu (p) + Nm (p) + Nv (p) + Nn (p) + &
            Na (p) + Ns (p) + Npa (p)
       ! Update C pools.
       ! Surface structural litter.
       Cu (p) = Cu (p) + dCu
       ! Surface metabolic litter.
       Cm (p) = Cm (p) + dCm
       ! Soil structural litter.
       Cv (p) = Cv (p) + dCv
       ! Soil metabolic litter.
       Cn (p) = Cn (p) + dCn
       ! Active.
       Ca (p) = Ca (p) + dCa
       ! Slow.
       Cs (p) = Cs (p) + dCs
       ! Passive.
       Cpa (p) = Cpa (p) + dCpa
       ! Following does not appear to be used.
       Cleach = dCle
       ! Update N pools.
       Nu (p) = Nu (p) + dNu
       Nm (p) = Nm (p) + dNm
       Nv (p) = Nv (p) + dNv
       Nn (p) = Nn (p) + dNn
       Na (p) = Na (p) + dNa
       Ns (p) = Ns (p) + dNs
       Npa (p) = Npa (p) + dNpa
       ! Total soil C (kg[C] m-2).
       soilC (p) = Cu (p) + Cm (p) + Cv (p) + &
                        Cn (p) + Ca (p) + Cs (p) + &
			Cpa (p)
       soilN (p) = Nu (p) + Nm (p) + Nv (p) + &
                        Nn (p) + Na (p) + Ns (p) + &
			Npa (p)
			
       ! N mineralization rate (kg[N] m-2 day-1).
       Nmin = N0 - soilN (p)
       ! Production rate of plant-available N (kg[N] m-2 day-1).
       Usoil = (one - emf) * (Nmin + Nad + Nfx + Ndepo)
       ! Soil respiration (kg[C] m-2 day-1).
       !****adfsresp = C0 - (soilC (p) + dCle)
       snmin (p) = snmin (p) + Usoil
       !---------------------------------------------------------------!
       ! Soil water holding capacity (Eqn. (1) of FW00; m).
       !---------------------------------------------------------------!
       swct (p) = 0.213 + 0.00227 * soilC (p)
       !---------------------------------------------------------------!
       
       !---------------------------------------------------------------!
       ! Soil water holding capacities of layers (m).
       ! d_one is max. swc of top layer for roots.
       !---------------------------------------------------------------!
       if (swct (p) <= 0.05) then
        swc1 (p) = swct (p)
        swc2 (p) = 0.0
        swc3 (p) = 0.0
       else
        if (swct (p) <= d_one) THEN
          swc1 (p) = 0.05
          swc2 (p) = swct (p) - 0.05
          swc3 (p) = 0.0
        else
          swc1 (p) = 0.05
          swc2 (p) = d_one - 0.05
          swc3 (p) = swct (p) - d_one
        end if
       end if
       wlittc (p) = zero
       wlittn (p) = zero
       flittc (p) = zero
       flittn (p) = zero
       rlittc (p) = zero
       rlittn (p) = zero
      end do ! kp
      !----------------------------------------------------------------!
      ! End of soil routine.
      !----------------------------------------------------------------!
      
end subroutine soil
