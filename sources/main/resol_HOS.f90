MODULE resol_HOS
!
! This module contains necessary routine to solve the Free Surface Boundary Conditions
! This includes the actual High-Order Spectral to compute the vertical velocity at FS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    Copyright (C) 2014 - LHEEA Lab., Ecole Centrale de Nantes, UMR CNRS 6598
!
!    This program is part of HOS-NWT
!
!    HOS-NWT is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
USE type
USE variables
USE Fourier_FFTW
USE filtering
USE dealiasing
USE resol_wmkr
USE variablechange
!
CONTAINS
!
! start solveHOS ******************************************************
!      
!     ==========================================================
subroutine solveHOS(phisrk,dphisrk,etark,detark,phis_add_rk,dphis_add_rk,isubsteprk,time_loc,c)
!     ==========================================================
!     
! This routine solves the Free Surface Boundary Conditions with 1st order wavemaker
!
IMPLICIT NONE
!
!% INPUT VARIABLES
!	free surface elevation
REAL(RP),DIMENSION(m1,m2),INTENT(IN)  :: etark
REAL(RP),DIMENSION(m1,m2),INTENT(OUT) :: detark
!	free surface potential
REAL(RP),DIMENSION(m1,m2),INTENT(IN)  :: phisrk
REAL(RP),DIMENSION(m1,m2),INTENT(OUT) :: dphisrk
!	wavemaker additional potential
REAL(RP),DIMENSION(m3_add,m2),INTENT(IN)  :: phis_add_rk
REAL(RP),DIMENSION(m3_add,m2),INTENT(OUT) :: dphis_add_rk
!
REAL(RP),INTENT(IN) :: time_loc
INTEGER ,INTENT(IN) :: isubsteprk
REAL(RP), OPTIONAL  :: c
!
!% LOCAL VARIABLES
!	extended FSBCs terms
REAL(RP),DIMENSION(md1,md2) :: detarkext,dphisrkext, phiz_etax_ext,phiz_etay_ext
REAL(RP),DIMENSION(m1,m2)   :: phiz_etay, phiz_etax
!
INTEGER                   :: i1,i2,i3
REAL(RP),DIMENSION(m1,m2) :: grad_phi
REAL(RP)                  :: dphin,dphix,dphiy,dphiz
!
REAL(RP) :: ti,tf
!
REAL(RP),DIMENSION(md1,md2) :: phix_add_ext, phiy_add_ext, phiz_add_ext
REAL(RP),DIMENSION(md1,md2) :: gradphis_gradphi_ext, gradeta_gradphi_ext
REAL(RP),DIMENSION(m1,m2)   :: gradphis_gradphi, gradeta_gradphi
! 
!	CPU times inlet
if(iCPUtime.eq.1) then
   print*,'entering subroutine solveHOS'
   call CPU_TIME(ti)
endif
!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! resolution of the different order of the FS potential boundary value problem to get modes
!
call phisxy_etaxy(phisrk,etark)
!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! calculation of the horizontal derivatives of phis and eta
!
call HOSphis_modes(phisrk)
!
! RHS of FSBCs correctly dealiased order p
do i1 = 1,nd1
   do i2=1,nd2
      detarkext(i1,i2) = - phisx(i1,i2) * etax(i1,i2) - phisy(i1,i2) * etay(i1,i2)
      dphisrkext(i1,i2) = - phisx(i1,i2) * phisx(i1,i2) - phisy(i1,i2) * phisy(i1,i2)
   end do
end do
CALL dealias(2,detarkext,'cos','cos')
CALL dealias(2,dphisrkext,'cos','cos')
!
do i1 = 1,nd1
   do i2=1,nd2
      detarkext(i1,i2)  = detarkext(i1,i2)  + geta2phiz(i1,i2)
      dphisrkext(i1,i2) = dphisrkext(i1,i2) + geta2phiz2(i1,i2) + phiz2(i1,i2)
   enddo
end do
!
call filter_ext(detarkext,detark,'cos','cos')
call filter_ext(dphisrkext,dphisrk,'cos','cos')
!
! Wavemaker part (only 1st order here)
IF ((icase==2.OR.&
	icase==3.OR.icase==31.OR.icase==32.OR.icase==33.OR.icase==4.OR.icase==41.OR.icase==42).AND.(igeom==1.OR.igeom==2)) THEN
   if (isubsteprk.eq.0) then
      call ramp
      call pos_1st
   endif
   !
   !	RHS of the wavemaker boundary condition
   do i3=1,n3
      do i2=1,n2
         dphis_add_rk(i3,i2) =  - alpha * d2pos1stdt2(i3,i2) - 2.d0*alphat*dpos1stdt(i3,i2) - alphatt*pos1st(i3,i2)
      end do
   end do
   !
   ! Matching surface : polynom to ensure the C3 character of phiadd + limited amplitude
   CALL match_surface_wmk(dphis_add_rk(1:n3_add,1:n2),n2,n3,n3_add,l_add,delz)
   !
   !   wavemaker additional modes resolution
   dphis_add_rk = space_2_Fourier_add(dphis_add_rk,'cos','cos')
   !
   do i3=1,n3_add
      do i2=1,n2
         dphis_add_rk(i3,i2) = dphis_add_rk(i3,i2) / k_add_thk_add(i3,i2)
      end do
   end do
   !
   ! Initial derivatives calculation
   !
   IF(ABS(time_loc) < tiny) THEN
      !
      !	re-construction of all 1st-order additional quantities needed
      call phiadd_1st_bis(phis_add_rk, dphis_add_rk, phi_add, phit_add,&
           phix_add, phiy_add, phiz_add, etark)
   ELSE
      !
      !	re-construction of all 1st-order additional quantities needed
      call phiadd_1st_bis(phis_add_rk, dphis_add_rk, phi_add, phit_add,&
           phix_add, phiy_add, phiz_add, etark)
   ENDIF
   !
ELSE IF ((icase==2.OR.icase==3.OR.icase==31.OR.icase==32.OR.icase==33.OR.icase==4.OR.icase==41.OR.icase==42).AND.(igeom==3)) THEN
   !
   call ramp
   call phi_dipoles(phi_add, phit_add, phix_add, phiz_add, etark)
ELSE IF (icase==1) THEN
   phi_add  = 0.0_rp
   phit_add = 0.0_rp
   phix_add = 0.0_rp
   phiz_add = 0.0_rp
END IF
!
do i1 = 1, nd1
   do i2=1,nd2
      phiz_etax_ext(i1,i2) = etax(i1,i2) * phizMm1(i1,i2)
      phiz_etay_ext(i1,i2) = etay(i1,i2) * phizMm1(i1,i2)
   end do
end do
!
call filter_ext(phiz_etax_ext,phiz_etax,'sin','cos')
call filter_ext(phiz_etay_ext,phiz_etay,'cos','sin')
!
! grad eta * grad phiadd
phix_add = space_2_Fourier(phix_add,'cos','cos')
phiy_add = space_2_Fourier(phiy_add,'cos','cos')
phiz_add = space_2_Fourier(phiz_add,'cos','cos')
!
phix_add_ext = extend(phix_add)
phiy_add_ext = extend(phiy_add)
phiz_add_ext = extend(phiz_add)
!
phix_add = Fourier_2_space(phix_add,'cos','cos')
phiy_add = Fourier_2_space(phiy_add,'cos','cos')
phiz_add = Fourier_2_space(phiz_add,'cos','cos')
phix_add_ext = Fourier_2_space_big(phix_add_ext,'cos','cos')
phiy_add_ext = Fourier_2_space_big(phiy_add_ext,'cos','cos')
phiz_add_ext = Fourier_2_space_big(phiz_add_ext,'cos','cos')
!
gradeta_gradphi_ext(1:nd1,1:nd2)  = etax(1:nd1,1:nd2)  * phix_add_ext(1:nd1,1:nd2) + etay(1:nd1,1:nd2)  * phiy_add_ext(1:nd1,1:nd2)
gradphis_gradphi_ext(1:nd1,1:nd2) = phisx(1:nd1,1:nd2) * phix_add_ext(1:nd1,1:nd2) + phisy(1:nd1,1:nd2) * phiy_add_ext(1:nd1,1:nd2)&
     + 0.5_rp * (phix_add_ext(1:nd1,1:nd2)*phix_add_ext(1:nd1,1:nd2) + phiy_add_ext(1:nd1,1:nd2)*phiy_add_ext(1:nd1,1:nd2) &
     + phiz_add_ext(1:nd1,1:nd2)*phiz_add_ext(1:nd1,1:nd2))
!
CALL filter_ext(gradeta_gradphi_ext,gradeta_gradphi,'cos','cos')
CALL filter_ext(gradphis_gradphi_ext,gradphis_gradphi,'cos','cos')
!
IF (mHOS ==1) THEN
   do i1=1,n1
      do i2=1,n2
         dphin          = phiz(i1,i2) + phiz_add(i1,i2) 
         detark(i1,i2)  = dphin
         dphisrk(i1,i2) = - etark(i1,i2) - phit_add(i1,i2) - nu(i1,i2)*dphin
      end do
   enddo
   IF (PRESENT(c)) THEN
      CALL filter_ext(phisx,gradeta_gradphi,'sin','cos')
      CALL filter_ext(phisy,gradphis_gradphi,'cos','sin')
      do i1=1,n1
         do i2=1,n2
            dphix           = gradeta_gradphi(i1,i2)  + phix_add(i1,i2)
            dphiy           = gradphis_gradphi(i1,i2) + phiy_add(i1,i2)
            dphiz           = phiz(i1,i2)      + phiz_add(i1,i2) 
            grad_phi(i1,i2) = SQRT(dphix*dphix + dphiy*dphiy + dphiz*dphiz)
         end do
      enddo
   END IF
ELSE
   do i1=1,n1
      do i2=1,n2
         dphin           = detark(i1,i2) + phiz(i1,i2) + phiz_add(i1,i2) &
              				- gradeta_gradphi(i1,i2)
         detark(i1,i2)   = dphin
         grad_phi(i1,i2) = 0.50_rp * dphisrk(i1,i2)
         dphisrk(i1,i2)  = 0.50_rp * dphisrk(i1,i2) - phit_add(i1,i2) - etark(i1,i2) &
              				- nu(i1,i2) * (dphin) - gradphis_gradphi(i1,i2)
      end do
   enddo
   !
   IF (PRESENT(c)) THEN
      CALL filter_ext(phisx,gradeta_gradphi,'sin','cos')
      CALL filter_ext(phisy,gradphis_gradphi,'cos','sin')
      do i1=1,n1
         do i2=1,n2
            dphix           = gradeta_gradphi(i1,i2)  - phiz_etax(i1,i2) + phix_add(i1,i2)
            dphiy           = gradphis_gradphi(i1,i2) - phiz_etay(i1,i2) + phiy_add(i1,i2)
            dphiz           = phiz(i1,i2) + phiz_add(i1,i2) 
            grad_phi(i1,i2) = SQRT(dphix*dphix + dphiy*dphiy + dphiz*dphiz)
         end do
      enddo
   END IF
END IF
!
IF (PRESENT(c)) THEN
   c  = MAXVAL(grad_phi)
END IF
!
!	CPU times outlet
if(iCPUtime.eq.1) then
   call CPU_TIME(tf)
   write(*,910)'quitting subroutine solveHOS, total CPU time: ',tf-ti,'s'
endif
!
910 format(a,1ES11.4,a)
!
return      
!
end subroutine solveHOS
!
! end solveHOS ********************************************************
!
!
!
! start solveHOS2 ******************************************************
!      
!     ==========================================================
      subroutine solveHOS2(a1strk,da1strk,eta1strk,deta1strk,a2ndrk,da2ndrk,eta2ndrk,deta2ndrk,G2,dG2,G1,dG1,phis_add_rk_1,&
           phis_add_rk_2,dphis_add_rk_1,dphis_add_rk_2,dphis_add_rkramp,phis_add_rk_3,dphis_add_rk_3, &
           isubsteprk, delta_t, c)
!     ==========================================================
!     
! This routine solves the Free Surface Boundary Conditions with different wavemaker nonlinearities
! Analytical integration of linear part
!
IMPLICIT NONE
!
!% INPUT VARIABLES
!
!	new variable G
REAL(RP),DIMENSION(m1,m2),INTENT(IN)  :: G1,G2
REAL(RP),DIMENSION(m1,m2),INTENT(OUT) :: dG1,dG2
!
!	wavemaker additional potential
REAL(RP),DIMENSION(m3_add,m2),INTENT(IN)  :: phis_add_rk_1,phis_add_rk_2,phis_add_rk_3
REAL(RP),DIMENSION(m3_add,m2),INTENT(OUT) :: dphis_add_rk_1,dphis_add_rk_2,dphis_add_rkramp,dphis_add_rk_3
!
! 1st and 2nd order free surface quantities
REAL(RP),DIMENSION(m1,m2),INTENT(IN)  :: a1strk,a2ndrk,eta1strk,eta2ndrk
REAL(RP),DIMENSION(m1,m2),INTENT(OUT) :: da1strk,da2ndrk,deta1strk,deta2ndrk
!
REAL(RP),INTENT(IN) :: delta_t
INTEGER ,INTENT(IN) :: isubsteprk
REAL(RP), OPTIONAL  :: c
!
!% LOCAL VARIABLES
!	reconstructed free surface elevation
!
REAL(RP),DIMENSION(m1,m2) :: etark,detark
!
!	reconstructed free surface potential
REAL(RP),DIMENSION(m1,m2) :: phisrk,dphisrk
!
!	Nonlinear terms in derived equation
!
REAL(RP),DIMENSION(m1,m2) :: NL1, NL2
!
!	extended FSBCs terms
REAL(RP),DIMENSION(md1,md2) :: detarkext,dphisrkext, phiz_etax_ext, phiz_etay_ext
REAL(RP),DIMENSION(m1,m2)   :: phiz_etax,phiz_etay
!
REAL(RP),DIMENSION(m3_add,m2) :: a_add_tot, da_add_tot
!
INTEGER :: i1,i2
REAL(RP),DIMENSION(m1,m2) :: grad_phi
REAL(RP) :: dphin, dphix,dphiy,dphiz
!
REAL(RP) :: ti,tf,dphisrk_1
!
REAL(RP),DIMENSION(md1,md2) :: phix_add_ext, phiy_add_ext, phiz_add_ext
REAL(RP),DIMENSION(md1,md2) :: gradeta_gradphi_ext, gradphis_gradphi_ext
REAL(RP),DIMENSION(m1,m2)   :: gradphis_gradphi, gradeta_gradphi
!
!	CPU times inlet
!
if(iCPUtime.eq.1) then
   print*,'entering subroutine solveHOS2'
   call CPU_TIME(ti)
endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! reconstruction of surface elevation and velocity potential
!
call GtoF(G1,G2,etark,phisrk,delta_t,iCPUtime)
!
!om2phis(1) = phisrk(1)*omega(1)
!
do i1=2,n1
   phisrk(i1,1) = phisrk(i1,1)/omega(i1,1)
end do
do i1=1,n1
   do i2=2,n2
      phisrk(i1,i2) = phisrk(i1,i2)/omega(i1,i2)
   enddo
enddo
!
! Correction of mean volume
!
! NL1 is the FS in space
NL1  = Fourier_2_space(etark,'cos','cos')
!
! Correct the volume, i.e. first mode on eta: set the value to theoretical value
! This is useful for very long time simulations
! FIXME: check the influence on Results
if(igeom.EQ.2) then
    etark(1,1) = 0.5_rp*alpha*(1-d_hinge)*pos1st(n3,1) / xlen + NL1(1,1)*alpha*pos1st(n3,1) / xlen
elseif(igeom.EQ.1) then
    etark(1,1) = alpha*pos1st(n3,1) / xlen + NL1(1,1)*alpha*pos1st(n3,1) / xlen
else
    etark(1,1) = 0.0_rp
endif
!
NL1 = 0.0_rp
!
etark  = Fourier_2_space(etark,'cos','cos')
phisrk = Fourier_2_space(phisrk,'cos','cos')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! resolution of the different order of the FS potential boundary value problem to get modes
call phisxy_etaxy(phisrk,etark)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! calculation of the horizontal derivatives of phis and eta
call HOSphis_modes(phisrk)
!
! RHS of FSBCs correctly dealiased order p
do i1 = 1,nd1
   do i2=1,nd2
      detarkext(i1,i2)  = - phisx(i1,i2) * etax(i1,i2)  - phisy(i1,i2) * etay(i1,i2)
      dphisrkext(i1,i2) = - phisx(i1,i2) * phisx(i1,i2) - phisy(i1,i2) * phisy(i1,i2)
   end do
end do
CALL dealias(2,detarkext,'cos','cos')
CALL dealias(2,dphisrkext,'cos','cos')
!
do i1 = 1,nd1
   do i2=1,nd2
      detarkext(i1,i2)  = detarkext(i1,i2)  + geta2phiz(i1,i2)
      dphisrkext(i1,i2) = dphisrkext(i1,i2) + geta2phiz2(i1,i2) + phiz2(i1,i2)
   end do
end do
!
call filter_ext(detarkext,detark,'cos','cos')
call filter_ext(dphisrkext,dphisrk,'cos','cos')
!
! Wavemaker part (1st,2nd or 3rd order)
!
IF ((icase==2.OR.icase==3.OR.icase==31.OR.icase==32.OR.icase==33.OR.icase==4.OR.icase==41.OR.icase==42).AND.&
	(igeom==1.OR.igeom==2)) THEN
   !
   IF(i_wmk.EQ.3) THEN ! 3rd order wavemaker condition
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! resolution HOS_SWEET_1: evaluates 1st order solution (FS + wmk)
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      CALL solveHOS_SWEET_1(a1strk,da1strk,eta1strk,deta1strk,phis_add_rk_1,&
           dphis_add_rk_1,dphis_add_rkramp,isubsteprk)
      !   
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! resolution HOS_SWEET_2: evaluates 2nd order solution (FS + wmk)
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      CALL solveHOS_SWEET_2(a2ndrk,da2ndrk,eta2ndrk,deta2ndrk,eta1strk,&
           phis_add_rk_2,dphis_add_rk_2,isubsteprk)
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! resolution HOS_SWEET_3: evaluates 3rd order solution (wmk)
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      CALL solveHOS_SWEET_3(a1strk,da1strk,a2ndrk,da2ndrk,phis_add_rk_1,dphis_add_rk_1,phis_add_rk_2,dphis_add_rk_2, &
           phis_add_rk_3,dphis_add_rk_3,isubsteprk)
      !
   ELSEIF(i_wmk.EQ.2) THEN ! 2nd order wavemaker condition
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! resolution HOS_SWEET_1: evaluates 1st order solution (FS + wmk)
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      CALL solveHOS_SWEET_1(a1strk,da1strk,eta1strk,deta1strk,phis_add_rk_1,&
           dphis_add_rk_1,dphis_add_rkramp,isubsteprk)
      !   
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! resolution HOS_SWEET_2: evaluates 2nd order solution (FS + wmk)
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      CALL solveHOS_SWEET_2(a2ndrk,da2ndrk,eta2ndrk,deta2ndrk,eta1strk,&
           phis_add_rk_2,dphis_add_rk_2,isubsteprk)
      !
      dphis_add_rk_3 = 0.d0
      phi_add_3   = 0.d0
      phit_add_3  = 0.d0
      phix_add_3  = 0.d0
      phiy_add_3  = 0.d0
      phiz_add_3  = 0.d0
      !
   ELSEIF(i_wmk.EQ.1) THEN ! 1st order wavemaker condition
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! resolution HOS_SWEET_1: evaluates 1st order solution (FS + wmk)
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      CALL solveHOS_SWEET_1(a1strk,da1strk,eta1strk,deta1strk,phis_add_rk_1,&
           dphis_add_rk_1,dphis_add_rkramp,isubsteprk)
      !
      dphis_add_rk_2 = 0.d0
      phi_add_2   = 0.d0
      phit_add_2  = 0.d0
      phix_add_2  = 0.d0
      phiy_add_2  = 0.d0
      phiz_add_2  = 0.d0
      !
      dphis_add_rk_3 = 0.d0
      phi_add_3   = 0.d0
      phit_add_3  = 0.d0
      phix_add_3  = 0.d0
      phiy_add_3  = 0.d0
      phiz_add_3  = 0.d0
      !
   ENDIF
   !
   a_add_tot(1:n3_add,1:n2)  = phis_add_rk_1(1:n3_add,1:n2)+phis_add_rk_2(1:n3_add,1:n2)+phis_add_rk_3(1:n3_add,1:n2)
   da_add_tot(1:n3_add,1:n2) = dphis_add_rk_1(1:n3_add,1:n2)+dphis_add_rk_2(1:n3_add,1:n2)+dphis_add_rk_3(1:n3_add,1:n2)
   !
   CALL phiadd_SL(a_add_tot,da_add_tot,phi_add,phix_add,phiy_add,phiz_add,phit_add, etark)
   !
ELSE IF ((icase==2.OR.icase==3.OR.icase==31.OR.icase==32.OR.icase==33.OR.icase==4.OR.icase==41.OR.icase==42).AND.(igeom==3)) THEN
   call ramp
   call phi_dipoles(phi_add, phit_add, phix_add, phiz_add, etark)
ELSE IF (icase==1) THEN
   phi_add  = 0.0_rp
   phit_add = 0.0_rp
   phix_add = 0.0_rp
   phiy_add = 0.0_rp
   phiz_add = 0.0_rp
   !
   dphis_add_rkramp = 0.d0
   dphis_add_rk_1   = 0.d0
   !
   dphis_add_rk_2 = 0.d0
   phi_add_2   = 0.d0
   phit_add_2  = 0.d0
   phix_add_2  = 0.d0
   phiy_add_2  = 0.d0
   phiz_add_2  = 0.d0
   !
   dphis_add_rk_3 = 0.d0
   phi_add_3   = 0.d0
   phit_add_3  = 0.d0
   phix_add_3  = 0.d0
   phiy_add_3  = 0.d0
   phiz_add_3  = 0.d0
END IF
!
do i1 = 1,nd1
   do i2=1,nd2
      phiz_etax_ext(i1,i2) = etax(i1,i2) * phizMm1(i1,i2)
      phiz_etay_ext(i1,i2) = etay(i1,i2) * phizMm1(i1,i2)
   end do
end do
!
call filter_ext(phiz_etax_ext,phiz_etax,'sin','cos')
call filter_ext(phiz_etay_ext,phiz_etay,'cos','sin')
!
! grad eta * grad phiadd
!
phix_add = space_2_Fourier(phix_add,'cos','cos')
phiy_add = space_2_Fourier(phiy_add,'cos','cos')
phiz_add = space_2_Fourier(phiz_add,'cos','cos')
!
phix_add_ext = extend(phix_add)
phiy_add_ext = extend(phiy_add)
phiz_add_ext = extend(phiz_add)
!
phix_add = Fourier_2_space(phix_add,'cos','cos')
phiy_add = Fourier_2_space(phiy_add,'cos','cos')
phiz_add = Fourier_2_space(phiz_add,'cos','cos')
phix_add_ext = Fourier_2_space_big(phix_add_ext,'cos','cos')
phiy_add_ext = Fourier_2_space_big(phiy_add_ext,'cos','cos')
phiz_add_ext = Fourier_2_space_big(phiz_add_ext,'cos','cos')
!
gradeta_gradphi_ext(1:nd1,1:nd2)  = etax(1:nd1,1:nd2)  * phix_add_ext(1:nd1,1:nd2) + etay(1:nd1,1:nd2)  * phiy_add_ext(1:nd1,1:nd2)
gradphis_gradphi_ext(1:nd1,1:nd2) = phisx(1:nd1,1:nd2) * phix_add_ext(1:nd1,1:nd2) + phisy(1:nd1,1:nd2) * phiy_add_ext(1:nd1,1:nd2)&
     + 0.5_rp * (phix_add_ext(1:nd1,1:nd2)*phix_add_ext(1:nd1,1:nd2) + phiy_add_ext(1:nd1,1:nd2)*phiy_add_ext(1:nd1,1:nd2) &
     + phiz_add_ext(1:nd1,1:nd2)*phiz_add_ext(1:nd1,1:nd2))
!
CALL filter_ext(gradeta_gradphi_ext,gradeta_gradphi,'cos','cos')
CALL filter_ext(gradphis_gradphi_ext,gradphis_gradphi,'cos','cos')
!
IF (mHOS ==1) THEN
   do i1=1,n1
      do i2=1,n2
         dphin          = phiz(i1,i2) + phiz_add(i1,i2) 
         detark(i1,i2)  = dphin - W1(i1,i2)
         dphisrk(i1,i2) = - phit_add(i1,i2) - nu(i1,i2)*dphin
      end do
   enddo
   IF (PRESENT(c)) THEN
      gradeta_gradphi_ext(:,:)  = phisx(:,:)
      gradphis_gradphi_ext(:,:) = phisy(:,:)
      !
      CALL filter_ext(gradeta_gradphi_ext,gradeta_gradphi,'sin','cos')
      CALL filter_ext(gradphis_gradphi_ext,gradphis_gradphi,'cos','sin')
      !
      do i1=1,n1
         do i2=1,n2
            dphix        = gradeta_gradphi(i1,i2) + phix_add(i1,i2)
            dphiy        = gradphis_gradphi(i1,i2) + phiy_add(i1,i2)
            dphiz        = phiz(i1,i2)   + phiz_add(i1,i2) 
            grad_phi(i1,i2) = SQRT(dphix*dphix + dphiy*dphiy + dphiz*dphiz)
         end do
      enddo
   END IF
ELSE
   do i1=1,n1
      do i2=1,n2
         dphin           = detark(i1,i2) + phiz(i1,i2) + phiz_add(i1,i2) &
             		 		- gradeta_gradphi(i1,i2)
         detark(i1,i2)   = dphin  - W1(i1,i2)
         grad_phi(i1,i2) = 0.50_rp * dphisrk(i1,i2)
         dphisrk(i1,i2)  = 0.50_rp * dphisrk(i1,i2) - phit_add(i1,i2) &
              				- nu(i1,i2) * (dphin) - gradphis_gradphi(i1,i2)
      end do
   enddo
   !
   IF (PRESENT(c)) THEN
      gradeta_gradphi_ext(:,:) = phisx(:,:)
      gradphis_gradphi_ext(:,:) = phisy(:,:)
      !
      CALL filter_ext(gradeta_gradphi_ext,gradeta_gradphi,'sin','cos')
      CALL filter_ext(gradphis_gradphi_ext,gradphis_gradphi,'cos','sin')
      do i1=1,n1
         do i2=1,n2
            dphix        = gradeta_gradphi(i1,i2) - phiz_etax(i1,i2) + phix_add(i1,i2)
            dphiy        = gradphis_gradphi(i1,i2) - phiz_etay(i1,i2) + phiy_add(i1,i2)
            dphiz        = phiz(i1,i2) + phiz_add(i1,i2) 
            grad_phi(i1,i2) = SQRT(dphix*dphix + dphiy*dphiy + dphiz*dphiz)
         end do
      enddo
   END IF
END IF
IF (PRESENT(c)) THEN
   c  = MAXVAL(grad_phi)
END IF
!
! Filtering if prescribed
IF (coeffilt(1) < 1.0_rp) THEN
   CALL filtering_x(detark,coeffilt(1))
   CALL filtering_x(dphisrk,coeffilt(1))
END IF
IF (n2 /= 1 .AND. coeffilt(2) < 1.0_rp) THEN
   CALL filtering_y(detark,coeffilt(2))
   CALL filtering_y(dphisrk,coeffilt(2))
END IF
!
! Construction of new variable G
do i1=1,n1
   do i2=1,n2
      NL1(i1,i2) = detark(i1,i2)
      NL2(i1,i2) = dphisrk(i1,i2)
   end do
end do
!
NL1 = space_2_Fourier(NL1,'cos','cos')
NL2 = space_2_Fourier(NL2,'cos','cos')
!
dphisrk_1 = NL2(1,1)
!
do i1=1,n1
   do i2=1,n2
      NL2(i1,i2) = NL2(i1,i2)*omega(i1,i2)
   end do
end do
!
call expon(NL1,NL2,dG1,dG2,delta_t,iCPUtime)
!
dG2(1,1) = dphisrk_1
!
!	CPU times outlet
if(iCPUtime.eq.1) then
   call CPU_TIME(tf)
   write(*,910)'quitting subroutine solveHOS2, total CPU time: ',tf-ti,'s'
endif
!
910 format(a,1ES11.4,a)
!
return      
!
end subroutine solveHOS2
!
! end solveHOS2 ********************************************************
!
!
!
! start HOSphis_modes *************************************************
!
!     ======================================
      subroutine HOSphis_modes(phisrk)
!     ======================================
!   
IMPLICIT NONE
!
!% INPUT VARIABLES
!
!	free surface potential
REAL(RP),DIMENSION(m1,m2),INTENT(IN) :: phisrk
!
!% LOCAL VARIABLES
!
!	order m of phi and i^th derivatives at order m
REAL(RP),DIMENSION(m1,m2)          :: phim
REAL(RP),DIMENSION(md1,md2,maxHOS) :: phizm,phizm2
REAL(RP),DIMENSION(md1,md2)        :: phizi,phimext
REAL(RP),DIMENSION(md1,md2)        :: temp1, temp2
REAL(RP),DIMENSION(m1,m2)          :: temp3, temp4
!	  
INTEGER :: i1, j, iHOS, jm1, j2, i2, jp1, mHOSm2mj, q
!
REAL(RP) :: ti,tf
!
!	CPU times inlet
!
if(iCPUtime.eq.1) then
   print*,'entering subroutine HOSphis_modes'
   call CPU_TIME(ti)
endif
!
! initializations
do i1 = 1, n1
   do i2=1,n2
      phim(i1,i2) = phisrk(i1,i2)
   end do
end do
!
do i1 = 1, nd1
   do i2=1, nd2
      phiz2(i1,i2)    = 0.d0
      phizMm1(i1,i2)  = 0.d0
      phizMm2(i1,i2)  = 0.d0
      phiz2Mm2(i1,i2) = 0.d0
      do j = 1, mHOS
         phizm(i1,i2,j) = 0.d0
      end do
   end do
end do
!
phim = space_2_Fourier(phim,'cos','cos')
!
aHOS_ext(:,:,1) = extend(phim)
!
phim = Fourier_2_space(phim,'cos','cos')
!
! calculation at each order iHOS, problem at next order: phim and phiz in orders of eta
do iHOS=1,mHOS
   do i1 = 1, nd1
      do i2=1, nd2
         phimext(i1,i2) = 0.d0
      end do
   end do
   !
   do j = iHOS,2,-1
      jm1 = j - 1
      j2  = iHOS - jm1  
      !
      DO i1 = 1, nd1
         DO i2 = 1, nd2
            phizi(i1,i2) = aHOS_ext(i1,i2,j2) * kth_ext(i1,i2,j)
         ENDDO
      ENDDO
      !
      phizi = Fourier_2_space_big(phizi,'cos','cos')
      DO i1 = 1, nd1
         DO i2 = 1, nd2
            !
            !	construction of phi_iHOS+1
            phimext(i1,i2) = phimext(i1,i2) - phizi(i1,i2) * etapmext(i1,i2,j)
            !
            !	construction of phiz components at each order in eta
            phizm(i1,i2,iHOS) = phizm(i1,i2,iHOS) + phizi(i1,i2) * etapmext(i1,i2,jm1)
         ENDDO
      end do
   end do
   !
   !	direct FFT: iHOS^th z-derivatives of phi_1
   DO i1 = 1, nd1
      DO i2 = 1, nd2
         phizi(i1,i2) = aHOS_ext(i1,i2,iHOS) * kth_ext(i1,i2,1)
      ENDDO
   ENDDO
   !
   !call FFTdir_csxcsy_big_IMSL(phizi)
   phizi = Fourier_2_space_big(phizi,'cos','cos')
   DO i1 = 1, nd1
      DO i2 = 1, nd2
         !
         !	last component of phi_iHOS+1
         phimext(i1,i2)=phimext(i1,i2)-phizi(i1,i2)*etapmext(i1,i2,1)
         !
         !	construction of phiz components at each order in eta
         phizm(i1,i2,iHOS) = phizm(i1,i2,iHOS) + phizi(i1,i2)
      ENDDO
   ENDDO
   !
   ! De-aliasing of phizm at order i_m
   IF (iHOS > 1) THEN
      CALL dealias(iHOS, phizm(:,:,iHOS),'cos','cos')
   END IF
   ! De-aliasing of phim at order i_m+1
   IF (iHOS < mHOS) THEN
      CALL dealias(iHOS+1, phimext,'cos','cos')
      ! Storing it for next i_m
      aHOS_ext(:,:,iHOS+1) = phimext(:,:)
      aHOS_ext(:,:,iHOS+1) = space_2_Fourier_big(aHOS_ext(:,:,iHOS+1),'cos','cos')
   END IF
end do
!
! assembling of phiz at orders M and M-2 (no dealiasing needed)
if (mHOS.ge.2) then
   do j=1,mHOS-2
      do i1 = 1,nd1
         do i2 = 1,nd2
            phizMm2(i1,i2)  = phizMm2(i1,i2) + phizm(i1,i2,j)
            phizm2(i1,i2,j) = phizm(i1,i2,j)
         end do
      end do
   end do
   do i1 = 1,nd1
      do i2 = 1,nd2
         phizMm1(i1,i2)  = phizMm2(i1,i2) + phizm(i1,i2,mHOS-1)
      end do
   end do
   !
   temp1(:,:) = phizMm1(:,:)
   temp2(:,:) = phizm(:,:,mHOS)
	CALL filter_ext(temp1,temp3,'cos','cos')
   	CALL filter_ext(temp2,temp4,'cos','cos')

   do i1=1,n1
      do i2=1,n2
         phiz(i1,i2)=temp3(i1,i2) + temp4(i1,i2)
      end do
   end do
else
   temp1(:,:) = phizm(:,:,1)
   CALL filter_ext(temp1,temp4,'cos','cos')
   do i1=1,n1
      do i2=1,n2
         phiz(i1,i2)=temp4(i1,i2)
      end do
   end do
endif
!
! Linear vertical velocity W1
temp1(1:nd1,1:nd2) = phizm(1:nd1,1:nd2,1)
CALL filter_ext(temp1,W1,'cos','cos')
!
if (mHOS.ge.4) then
   do j=1,mHOSm3div2
      mHOSm2mj=mHOS-2-j
      jp1=j+1
      do i1=1,nd1
         do i2=1,nd2
            phiz2Mm2(i1,i2)   = phiz2Mm2(i1,i2)   + 2.d0*phizm2(i1,i2,j)*phizm2(i1,i2,mHOSm2mj)
            phizm2(i1,i2,jp1) = phizm2(i1,i2,jp1) + phizm2(i1,i2,j)
         end do
      end do
   end do
   do i1=1,nd1
      do i2=1,nd2
         phiz2Mm2(i1,i2) = phiz2Mm2(i1,i2)+phizm2(i1,i2,mHOSdiv2m1)**2
      end do
   end do
endif
!
! Assembling of gradeta2*phiz at order M (i.e. with phiz at order M-2)
!  with correct p dealiasing
geta2phiz(1:nd1,1:nd2) = 0.0_rp
DO j=1,mHOS-3
   temp1(1:nd1,1:nd2) = gradeta2(1:nd1,1:nd2) * phizm(1:nd1,1:nd2,j)
   CALL dealias(j+2,temp1,'cos','cos')
   geta2phiz(1:nd1,1:nd2) = geta2phiz(1:nd1,1:nd2) + temp1(1:nd1,1:nd2)
END DO
IF (mHOS >= 3) THEN
   temp1(1:nd1,1:nd2) = gradeta2(1:nd1,1:nd2) * phizm(1:nd1,1:nd2,MAX(mHOS-2,1))
   CALL dealias(mHOS,temp1,'cos','cos')
   geta2phiz(1:nd1,1:nd2) = geta2phiz(1:nd1,1:nd2) + temp1(1:nd1,1:nd2)
END IF
!
! Assembling of phiz2 at order M and 
!  gradeta2*phiz at order M (i.e. with phiz at order M-2) 
!  with correct p dealiasing
!  in two steps
phiz2(1:nd1,1:nd2)      = 0.0_rp
geta2phiz2(1:nd1,1:nd2) = 0.0_rp
temp1(1:nd1,1:nd2)  = 0.0_rp
temp2(1:nd1,1:nd2)  = 0.0_rp
!
! Step 1: phiz2 up to order M-2 and gradeta2*phiz2 at order M
! squares
DO j = 1, (mHOS-2)/2
   temp2(1:nd1,1:nd2) = phizm(1:nd1,1:nd2,j) * phizm(1:nd1,1:nd2,j)
   CALL dealias(2*j,temp2,'cos','cos')
   temp1(1:nd1,1:nd2)  = temp1(1:nd1,1:nd2) + temp2(1:nd1,1:nd2)
   temp2(1:nd1,1:nd2) =  gradeta2(1:nd1,1:nd2) * temp2(1:nd1,1:nd2)
   CALL dealias(2*j+2,temp2,'cos','cos')
   geta2phiz2(1:nd1,1:nd2) = geta2phiz2(1:nd1,1:nd2) + temp2(1:nd1,1:nd2)
END DO
! double products
DO j = 2, mHOS-3
   DO q = 1, MIN(mHOS-2-j,j-1)
      temp2(1:nd1,1:nd2) = 2.0_rp * phizm(1:nd1,1:nd2,j) * phizm(1:nd1,1:nd2,q)
      CALL dealias(j+q,temp2,'cos','cos')
      temp1(1:nd1,1:nd2)  = temp1(1:nd1,1:nd2) + temp2(1:nd1,1:nd2)
      temp2(1:nd1,1:nd2) =  gradeta2(1:nd1,1:nd2) * temp2(1:nd1,1:nd2)
      CALL dealias(j+q+2,temp2,'cos','cos')
      geta2phiz2(1:nd1,1:nd2) = geta2phiz2(1:nd1,1:nd2) + temp2(1:nd1,1:nd2)
   END DO
END DO
! Step 2: phiz2 up to order M
! squares
IF (mHOS >= 2) THEN
   temp2(1:nd1,1:nd2) = phizm(1:nd1,1:nd2,mHOS/2) * phizm(1:nd1,1:nd2,mHOS/2)
   CALL dealias(2*(mHOS/2),temp2,'cos','cos')
   temp1(1:nd1,1:nd2)  = temp1(1:nd1,1:nd2) + temp2(1:nd1,1:nd2)
END IF
! double products
DO j = 2, mHOS-1
   DO q = MAX(0,MIN(mHOS-2-j,j-1))+1, MIN(mHOS-j,j-1)
      temp2(1:nd1,1:nd2) = 2.0_rp * phizm(1:nd1,1:nd2,j) * phizm(1:nd1,1:nd2,q)
      CALL dealias(j+q,temp2,'cos','cos')
      temp1(1:nd1,1:nd2)  = temp1(1:nd1,1:nd2) + temp2(1:nd1,1:nd2)
   END DO
END DO
!
phiz2(1:nd1,1:nd2) = temp1(1:nd1,1:nd2)
!
! this looks surprising but the cases M=2, 3 and 4 are correctly treated
! (cf work done in december 2007)
!
!	CPU times outlet
if (iCPUtime.eq.1) then
   call CPU_TIME(tf)
   write(*,910)'quitting subroutine HOSphis_modes, total CPU time:',tf-ti,'s'
endif
!
910 format(a,1ES11.4,a)
!
end subroutine HOSphis_modes
!
! end HOSphis_modes ***************************************************
!
!
! 
! start phisxy_etaxy **************************************************
!
!     =====================================
subroutine phisxy_etaxy(phisrk,etark)
!     =====================================
!
IMPLICIT NONE
!
!% INPUT VARIABLES
!	free surface potential
REAL(RP),DIMENSION(m1,m2),INTENT(IN) :: phisrk
!	free surface elevation
REAL(RP),DIMENSION(m1,m2),INTENT(IN) :: etark
!	  
!% LOCAL VARIABLES
!	modes of phisrk or etark
REAL(RP),DIMENSION(m1,m2) :: as !etalg(4*m12t4)
INTEGER  :: i1, i2, j, jp1
REAL(RP) :: ti, tf
!
!	CPU times inlet
if (iCPUtime.eq.1) then
   print*,'entering subroutine phisxy_etaxy'
   call CPU_TIME(ti)
endif
!
! inverse FFT: determination of the modes of eta
do i1=1,n1
   do i2=1,n2
      as(i1,i2)=etark(i1,i2) 
   end do
end do
!
as = space_2_Fourier(as,'cos','cos')
!
do i1=1,n1
   do i2=1,n2
      etapmext(i1,i2,1)=as(i1,i2)
   end do
   do i2=n2+1,nd2
      etapmext(i1,i2,1)=0.d0
   end do
end do
!
do i1=n1+1,nd1
   do i2=1,nd2
      etapmext(i1,i2,1)=0.d0
   end do
end do
!
etapmext(:,:,1) = Fourier_2_space_big(etapmext(:,:,1),'cos','cos')
!
!	de-aliasing of required powers of eta
do j=1,mHOS-1
   jp1 = j+1
   do i1=1,nd1
      do i2=1,nd2
         etapmext(i1,i2,jp1) = etapmext(i1,i2,j) / jp1 * etapmext(i1,i2,1)
      end do
   end do
   call dealias(jp1,etapmext(:,:,jp1),'cos','cos')
end do
!
!   calculation of etax by direct FFT
do i1=1,n1
   do i2=1,n2
      etax(i1,i2)=-kx(i1)*as(i1,i2)
   end do
   do i2=n2+1,nd2
      etax(i1,i2)=0.d0
   end do
end do
!
do i1=n1+1,nd1
   do i2=1,nd2
      etax(i1,i2)=0.d0
   end do
end do
!
etax = Fourier_2_space_big(etax,'sin','cos')
!
do i1=1,n1
   do i2=1,n2
      etay(i1,i2)=-ky(i2)*as(i1,i2)
   end do
   do i2=n2+1,nd2
      etay(i1,i2)=0.d0
   end do
end do
!
do i1=n1+1,nd1
   do i2=1,nd2
      etay(i1,i2)=0.d0
   end do
end do
!
etay = Fourier_2_space_big(etay,'cos','sin')
! 
! grad(eta)^2 (on a (p+1)/2 N grid)
gradeta2(1:nd1,1:nd2) = etax(1:nd1,1:nd2) * etax(1:nd1,1:nd2) + etay(1:nd1,1:nd2) * etay(1:nd1,1:nd2)
CALL dealias(2,gradeta2,'cos','cos')
!
! inverse FFT: determination of the modes of phis
do i1=1,n1
   do i2=1,n2
      as(i1,i2)=phisrk(i1,i2)
   end do
end do
!
as = space_2_Fourier(as,'cos','cos')
!
!   calculation of phisx by direct FFT
do i1=1,n1
   do i2=1,n2
      phisx(i1,i2)=-kx(i1)*as(i1,i2)
   end do
   do i2=n2+1,nd2
      phisx(i1,i2)=0.d0
   end do
end do
!
do i1=n1+1,nd1
   do i2=1,nd2
      phisx(i1,i2)=0.d0
   end do
end do
!
phisx = Fourier_2_space_big(phisx,'sin','cos')
!
!   calculation of phisy by direct FFT
do i1=1,n1
   do i2=1,n2
      phisy(i1,i2)=-ky(i2)*as(i1,i2)
   end do
   do i2=n2+1,nd2
      phisy(i1,i2)=0.d0
   end do
end do
!
do i1=n1+1,nd1
   do i2=1,nd2
      phisy(i1,i2)=0.d0
   end do
end do
!
phisy = Fourier_2_space_big(phisy,'cos','sin')
!
!	CPU times outlet
if(iCPUtime.eq.1) then
   call CPU_TIME(tf)
   write(*,910)'quitting subroutine phisxy_etaxy, total CPU time:',tf-ti,'s'
endif
!
910 format(a,1ES11.4,a)
!
return
!
end subroutine phisxy_etaxy
!
! end phisxy_etaxy ****************************************************
!
!
! end module
!
END MODULE resol_HOS
