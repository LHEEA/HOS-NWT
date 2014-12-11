MODULE resol_wmkr
!
! This module contains the different subroutines for wavemaker resolution in HOS-NWT
! This is done through an additional potential solved at different orders (up to 3rd)
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
USE wavemaker
USE Fourier_FFTW
USE filtering
!
CONTAINS
!
!
! start phiadd_1st ****************************************************
!      
!     ==================================================================
subroutine phiadd_1st(a1st_add_l,da1st_add_l,phi2_1st,phi2t_1st,phi2x_1st,phi2y_1st,phi2z_1st, eta_l)
!     ==================================================================
!     
implicit none
!
!% INPUT VARIABLES
!	first-order additional time modal amplitudes
REAL(RP), DIMENSION(m3_add,m2), INTENT(IN) :: a1st_add_l,da1st_add_l
!	first-order additional quantities needed
!	  on the free surface
REAL(RP), DIMENSION(m1,m2), INTENT(IN)  :: eta_l
REAL(RP), DIMENSION(m1,m2), INTENT(OUT) :: phi2_1st,phi2t_1st,phi2x_1st,phi2y_1st,phi2z_1st
!
!% LOCAL VARIABLES
REAL(RP) :: ti, tf
INTEGER  :: i1, i2, i3
REAL(RP) :: kx_add_etam1, coskx_add_eta, sinkx_add_eta
!
!	CPU times inlet
if (iCPUtime.eq.1) then
   print*,'entering subroutine phiadd_1st'
   call CPU_TIME(ti)
endif
!
! re-construction on the free surface (z=eta)	
do i2=1,n2
   do i1=1,n1
      phi2_1st(i1,i2)   = 0.0_rp
      phi2t_1st(i1,i2)  = 0.0_rp
      phi2x_1st(i1,i2)  = 0.0_rp
      phi2y_1st(i1,i2)  = 0.0_rp
      phi2z_1st(i1,i2)  = 0.0_rp
      do i3=1,n3_add
         kx_add_etam1  = kx_add(i3) * (eta_l(i1,i2)+1.0_rp)
         coskx_add_eta = cos(kx_add_etam1)
         sinkx_add_eta = sin(kx_add_etam1)
         phi2_1st(i1,i2)   = phi2_1st(i1,i2)   +  a1st_add_l(i3,i2) * coskx_add_eta * csh_add_x(i1,i2,i3)
         phi2t_1st(i1,i2)  = phi2t_1st(i1,i2)  + da1st_add_l(i3,i2) * coskx_add_eta * csh_add_x(i1,i2,i3)
         phi2x_1st(i1,i2)  = phi2x_1st(i1,i2)  -  a1st_add_l(i3,i2) * coskx_add_eta * k_add_sh_add_x(i1,i2,i3)
         phi2y_1st(i1,i2)  = phi2y_1st(i1,i2)  -  a1st_add_l(i3,i2) * coskx_add_eta * kycsh_add_x(i1,i2,i3)
         phi2z_1st(i1,i2)  = phi2z_1st(i1,i2)  -  a1st_add_l(i3,i2) * sinkx_add_eta * kx_add_csh_add_x(i1,i2,i3)
      end do
   end do
end do
!
do i1=1,n1
   phi2_1st(i1,1:n2)  = Fourier_2_space_y(phi2_1st(i1,1:n2), 'cos')
   phi2t_1st(i1,1:n2) = Fourier_2_space_y(phi2t_1st(i1,1:n2),'cos')
   phi2x_1st(i1,1:n2) = Fourier_2_space_y(phi2x_1st(i1,1:n2),'cos')
   phi2y_1st(i1,1:n2) = Fourier_2_space_y(phi2y_1st(i1,1:n2),'sin')
   phi2z_1st(i1,1:n2) = Fourier_2_space_y(phi2z_1st(i1,1:n2),'cos')
end do
!
!	CPU times outlet
if (iCPUtime.eq.1) then
   call CPU_TIME(tf)
   write(*,910)'quitting subroutine phiadd_1st, total CPU time: ',tf-ti,'s'
endif
!
910 format(a,1ES11.4,a)
!
end subroutine phiadd_1st
!
! end phiadd_1st ******************************************************
!
!
!
! start solveHOS_SWEET_1 ******************************************************
!      
!     ==========================================================
      subroutine solveHOS_SWEET_1(a1strk_l,da1strk_l,eta1strk_l,deta1strk_l,phis_add_rk_1_loc,&
     					   dphis_add_rk_1_loc,dphis_add_rkramp_loc, isubsteprk)
!     ==========================================================S
!
IMPLICIT NONE  
!
!% INPUT VARIABLES
!	free surface elevation
REAL(RP), DIMENSION(m1,m2), INTENT(IN)  :: eta1strk_l,a1strk_l
REAL(RP), DIMENSION(m1,m2), INTENT(OUT) :: da1strk_l,deta1strk_l
!
!	wavemaker additional potential
REAL(RP), DIMENSION(m3_add,m2), INTENT(IN)  :: phis_add_rk_1_loc
REAL(RP), DIMENSION(m3_add,m2), INTENT(OUT) :: dphis_add_rk_1_loc ,dphis_add_rkramp_loc
!
INTEGER, INTENT(IN)  :: isubsteprk
!
! Local variables
REAL(RP) :: ti, tf
INTEGER  :: i1,i2,i3
! 
!	CPU times inlet
if(iCPUtime.eq.1) then
   print*,'entering subroutine solveHOS'
   call CPU_TIME(ti)
endif
!
! Wavemaker part
if (isubsteprk.eq.0) then
   call ramp
   call pos_1st
endif
!
!	RHS of the wavemaker boundary condition
do i3=1,n3
   do i2=1,n2
      dphis_add_rk_1_loc(i3,i2)   = - alpha * d2pos1stdt2(i3,i2) - 2.d0*alphat*dpos1stdt(i3,i2) - alphatt*pos1st(i3,i2)
      dphis_add_rkramp_loc(i3,i2) = - alphat*dpos1stdt(i3,i2) - alphatt*pos1st(i3,i2)
   end do
end do
!
! Matching surface : polynom to ensure the C3 character of additional quantities + limited amplitude
CALL match_surface_wmk(dphis_add_rkramp_loc(1:n3_add,1:n2),n2,n3,n3_add,l_add,delz)
CALL match_surface_wmk(dphis_add_rk_1_loc(1:n3_add,1:n2),n2,n3,n3_add,l_add,delz)
!
!   wavemaker additional modes resolution
dphis_add_rk_1_loc   = space_2_Fourier_add(dphis_add_rk_1_loc, 'cos','cos')
dphis_add_rkramp_loc = space_2_Fourier_add(dphis_add_rkramp_loc, 'cos','cos')
!
IF (coeffilt(3) < 1.0_rp) THEN
   CALL filtering_z_modes(dphis_add_rk_1_loc,coeffilt(3))
   CALL filtering_z_modes(dphis_add_rkramp_loc,coeffilt(3))
END IF
!
do i3=1,n3_add
   do i2=1,n2
      dphis_add_rk_1_loc(i3,i2) = dphis_add_rk_1_loc(i3,i2) / k_add_thk_add(i3,i2)
      dphis_add_rkramp_loc(i3,i2) = dphis_add_rkramp_loc(i3,i2)  / k_add_thk_add(i3,i2)
   end do
end do
!
!	re-construction of all 1st-order additional quantities needed
call phiadd_1st_SWEET(phis_add_rk_1_loc,dphis_add_rk_1_loc,phi_add_1,phit_add_1,phix_add_1,&
	phiy_add_1,phiz_add_1,phizz_add_1,phixt_add_1,phiyt_add_1,&
	phizt_add_1,phixx_add_1,phixz_add_1,phiyz_add_1,&
	phi2xx_add_1st,phi2xxt_add_1st,phi2y_add_1st,phi2yt_add_1st,&
	phi2z_add_1st,phi2zt_add_1st)
!
!   1st-order free surface potential and its vertical derivative
call phidphiz(a1strk_l,phi1st,phiz1st,phi_add_1,phiz_add_1)
!
!   RHS of the FS boundary conditions
do i1=1,n1
   do i2=1,n2
      da1strk_l(i1,i2)   = - eta1strk_l(i1,i2) - phit_add_1(i1,i2) - iabsnb * nu(i1,i2) * phiz1st(i1,i2)
      deta1strk_l(i1,i2) =   phiz1st(i1,i2)
   end do
end do
!
!   free surface modes resolution
da1strk_l = space_2_Fourier(da1strk_l, 'cos','cos')
!
IF((i_wmk.EQ.2).OR.(i_wmk.EQ.3)) THEN
   !
   !	other first order derivatives required at second order
   call d1st_for_2nd(a1strk_l,da1strk_l,phix_add_1,phiy_add_1,phizz_add_1,&
        phixt_add_1,phiyt_add_1,phizt_add_1,phixx_add_1,&
        phixz_add_1,phiyz_add_1,phi2xx_add_1st,phi2xxt_add_1st,phi2y_add_1st&
        ,phi2yt_add_1st,phi2z_add_1st,phi2zt_add_1st)
ENDIF
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
end subroutine solveHOS_SWEET_1
!
! end solveHOS_SWEET_1 ********************************************************
!
!
!
! start solveHOS_SWEET_2 ******************************************************
!      
!     ==========================================================
      subroutine solveHOS_SWEET_2(a2ndrk,da2ndrk_loc,eta2ndrk,deta2ndrk,eta1strk,&
     					a2nd_add_rk,da2nd_add_rk,isubsteprk)
!     ==========================================================
! 
IMPLICIT NONE
!
!	second-order rk4 spectral time modal amplitudes
REAL(RP), DIMENSION(m1,m2), INTENT(IN)  :: a2ndrk
REAL(RP), DIMENSION(m1,m2), INTENT(OUT) :: da2ndrk_loc
!
!	second-order rk4 additional time modal amplitudes
REAL(RP), DIMENSION(m3_add,m2), INTENT(IN)  :: a2nd_add_rk
REAL(RP), DIMENSION(m3_add,m2), INTENT(OUT) :: da2nd_add_rk
!
!	second and first-order rk4 water elevation
REAL(RP), DIMENSION(m1,m2), INTENT(IN)  :: eta2ndrk,eta1strk
REAL(RP), DIMENSION(m1,m2), INTENT(OUT) :: deta2ndrk
!
INTEGER,INTENT(IN)  :: isubsteprk
!
!% LOCAL VARIABLES
REAL(RP) :: ti, tf, a2nd_add_noramp
INTEGER  :: i1, i2, i3
!
!	CPU times inlet
!
if(iCPUtime.eq.1) then
   print*,'entering subroutine solveHOS_SWEET_2'
   call CPU_TIME(ti)
endif
!
! Wavemaker part
!
if (isubsteprk.eq.0) then
   call dposdy_1st
   call dposdz_1st
   if ((icase.eq.2).AND.(ifree.eq.1)) then
      call pos_2nd
   end if
end if
!
!	RHS of the wavemaker boundary condition
!
do i3=1,n3
   do i2=1,n2
      a2nd_add_noramp  = - pos1st(i3,i2)    * phixx1st_add(i3,i2)&
           + dpos1stdy(i3,i2) * phiy1st_add(i3,i2)&
           + dpos1stdz(i3,i2) * phiz1st_add(i3,i2)

      da2nd_add_rk(i3,i2) = - alpha * (- pos1st(i3,i2)       * phixxt1st_add(i3,i2)&
           - dpos1stdt(i3,i2)    * phixx1st_add(i3,i2)&
           + dpos1stdy(i3,i2)    * phiyt1st_add(i3,i2)&
           + d2pos1stdydt(i3,i2) * phiy1st_add(i3,i2)&
           + dpos1stdz(i3,i2)    * phizt1st_add(i3,i2)&
           + d2pos1stdzdt(i3,i2) * phiz1st_add(i3,i2))&
           - alphat * a2nd_add_noramp&
           - ifree   * (alpha           * d2pos2nddt2(i3,i2)&
           + 2.d0 * alphat * dpos2nddt(i3,i2)&
           + alphatt       * pos2nd(i3,i2))
   end do
end do
!
! Matching surface : polynom to ensure the C3 character of additional quantities + limited amplitude
CALL match_surface_wmk(da2nd_add_rk(1:n3_add,1:n2),n2,n3,n3_add,l_add,delz)
!
!	wavemaker additional modes resolution
da2nd_add_rk = space_2_Fourier_add(da2nd_add_rk, 'cos','cos')
!
IF (coeffilt(3) < 1.0_rp) THEN
    CALL filtering_z_modes(da2nd_add_rk,coeffilt(3))
END IF
!
do i3=1,n3_add
   do i2=1,n2
      da2nd_add_rk(i3,i2)=da2nd_add_rk(i3,i2)/k_add_thk_add(i3,i2)
   end do
end do
!
!	re-construction of all 2nd-order additional quantities needed
call phiadd_2nd(a2nd_add_rk,da2nd_add_rk,phi_add_2,phix_add_2,phiy_add_2,phit_add_2,phiz_add_2)
!
!   2nd-order free surface potential and its vertical derivative
call phidphiz(a2ndrk,phi2nd,phiz2nd,phi_add_2,phiz_add_2)
!
!   RHS of the FS boundary conditions
do i1=1,n1
   do i2=1,n2
      deta2ndrk(i1,i2)=phiz2nd(i1,i2)&
           -phix1st(i1,i2)*etax1st(i1,i2)-phiy1st(i1,i2)*etay1st(i1,i2)&
           +eta1strk(i1,i2)*phizz1st(i1,i2)		
      da2ndrk_loc(i1,i2)=-eta2ndrk(i1,i2)-phit_add_2(i1,i2)&
           -5.d-1*(phix1st(i1,i2)**2+phiy1st(i1,i2)**2+phiz1st(i1,i2)**2)&
           -eta1strk(i1,i2)*phizt1st(i1,i2)&
           -iabsnb*nu(i1,i2)*deta2ndrk(i1,i2)
   end do
end do
!
!   free surface modes resolution
da2ndrk_loc = space_2_Fourier(da2ndrk_loc, 'cos','cos')
!
!	CPU times outlet
if(iCPUtime.eq.1) then
   call CPU_TIME(tf)
   write(*,910)'quitting subroutine solveHOS_SWEET_2, total CPU time: ',tf- ti,'s'
endif
!
910 format(a,1ES11.4,a)
!
return      
!
end subroutine solveHOS_SWEET_2
!
! end solveHOS_SWEET_2 ********************************************************
!
!
!
! start solveHOS_SWEET_3 ******************************************************
!      
!     ==========================================================
      subroutine solveHOS_SWEET_3(a1strk,da1strk,a2ndrk,da2ndrk,a1st_add_rk,da1st_add_rk,a2nd_add_rk,da2nd_add_rk, &
           a3rd_add_rk,da3rd_add_rk,isubsteprk)
!     ==========================================================
!
IMPLICIT NONE
!
!	First & second-order rk4 spectral time modal amplitudes
REAL(RP), DIMENSION(m1,m2), INTENT(IN) :: a1strk, da1strk, a2ndrk, da2ndrk
!
! First-Second-Third order wavemaker time modal amplitudes
REAL(RP), DIMENSION(m3_add,m2), INTENT(IN)  :: a1st_add_rk,da1st_add_rk,a2nd_add_rk,da2nd_add_rk,a3rd_add_rk
REAL(RP), DIMENSION(m3_add,m2), INTENT(OUT) :: da3rd_add_rk
!
INTEGER,INTENT(IN)  :: isubsteprk
!
!% LOCAL VARIABLES
REAL(RP) :: a3rd_add_alpha, a3rd_add_alphat, a3rd_add_alpha2, a3rd_add_alphaalphat
!
REAL(RP) :: ti, tf
INTEGER  :: i2, i3
!
!	CPU times inlet
if(iCPUtime.eq.1) then
   print*,'entering subroutine solveHOS_SWEET_3'
   call CPU_TIME(ti)
endif
!
! Wavemaker part
if (isubsteprk.eq.0) then
   if ((icase.eq.2).AND.(ifree.eq.1)) then
      write(*,*) 'check signs of 2nd order'
      call dposdy_2nd
      call dposdz_2nd
   end if
end if
!
! First and second order quantities needed to compute third order wavemaker
CALL wmk_3rdorder(a1strk,da1strk,a2ndrk,da2ndrk,a1st_add_rk,da1st_add_rk,a2nd_add_rk,da2nd_add_rk)
!
!	RHS of the wavemaker boundary condition
do i3=1,n3
   do i2=1,n2
      a3rd_add_alpha   =  - alpha * (- dpos1stdt(i3,i2) * phixx2nd_add(i3,i2) &
           - pos1st(i3,i2)    * phixxt2nd_add(i3,i2) &
           + d2pos1stdydt(i3,i2) * phiy2nd_add(i3,i2)&
           + dpos1stdy(i3,i2) * phiyt2nd_add(i3,i2)&
           + d2pos1stdzdt(i3,i2) * phiz2nd_add(i3,i2)&
           + dpos1stdz(i3,i2) * phizt2nd_add(i3,i2))
      !
      a3rd_add_alphat  =  - alphat * (- pos1st(i3,i2)    * phixx2nd_add(i3,i2)&
           + dpos1stdy(i3,i2) * phiy2nd_add(i3,i2)&
           + dpos1stdz(i3,i2) * phiz2nd_add(i3,i2))
      !
      a3rd_add_alpha2  =  - alpha**2*(- (pos1st(i3,i2)**2)/2.0_rp   * phixxxt1st_add(i3,i2)&
           - dpos1stdt(i3,i2) * pos1st(i3,i2)* phixxx1st_add(i3,i2)&
           + dpos1stdt(i3,i2) *dpos1stdy(i3,i2) * phixy1st_add(i3,i2)&
           + pos1st(i3,i2)*dpos1stdy(i3,i2)*phixyt1st_add(i3,i2)&
           + pos1st(i3,i2)*d2pos1stdydt(i3,i2)*phixy1st_add(i3,i2)&
           + dpos1stdt(i3,i2) *dpos1stdz(i3,i2) * phixz1st_add(i3,i2)&
           + pos1st(i3,i2)*dpos1stdz(i3,i2)*phixzt1st_add(i3,i2)&
           + pos1st(i3,i2)*d2pos1stdzdt(i3,i2)*phixz1st_add(i3,i2))
      !
      a3rd_add_alphaalphat = - alpha*alphat*(-pos1st(i3,i2)**2 * phixxx1st_add(i3,i2)&
           + 2.0_rp*pos1st(i3,i2)*dpos1stdy(i3,i2)*phixy1st_add(i3,i2)&
           + 2.0_rp*pos1st(i3,i2)*dpos1stdz(i3,i2)*phixz1st_add(i3,i2))
      !
      da3rd_add_rk(i3,i2) = a3rd_add_alpha + a3rd_add_alphat + a3rd_add_alpha2 + a3rd_add_alphaalphat
   end do
end do
!
! Matching surface : polynom to ensure the C3 character of phiadd + limited amplitude
CALL match_surface_wmk(da3rd_add_rk(1:n3_add,1:n2),n2,n3,n3_add,l_add,delz)
!
!	wavemaker additional modes resolution
da3rd_add_rk = space_2_Fourier_add(da3rd_add_rk, 'cos','cos')
!
IF (coeffilt(3) < 1.0_rp) THEN
    CALL filtering_z_modes(da3rd_add_rk,coeffilt(3))
END IF
!
do i3=1,n3_add
   do i2=1,n2
      da3rd_add_rk(i3,i2)=da3rd_add_rk(i3,i2)/k_add_thk_add(i3,i2)
   end do
end do
!
!	re-construction of all 3rd-order additional quantities needed
call phiadd_3rd(a3rd_add_rk,da3rd_add_rk,phi_add_3,phix_add_3,phiy_add_3,phit_add_3,phiz_add_3)
!
!	CPU times outlet
if(iCPUtime.eq.1) then
   call CPU_TIME(tf)
   write(*,910)'quitting subroutine solveHOS_SWEET_3, total CPU time: ',tf- ti,'s'
endif
!
910 format(a,1ES11.4,a)
!
return      
!
end subroutine solveHOS_SWEET_3
!
! end solveHOS_SWEET_3 ********************************************************
!
!
!
! start phiadd_1st_SWEET ****************************************************
!      
!     ==================================================================
      subroutine phiadd_1st_SWEET(a1st_add_l,da1st_add_l,phi2_1st,phi2t_1st,phi2x_1st,&
     			phi2y_1st,phi2z_1st,phi2zz_1st,phi2xt_1st,phi2yt_1st,&
     			phi2zt_1st,phi2xx_1st,phi2xz_1st,phi2yz_1st&
     			,phi2xx_add_1st_l,phi2xxt_add_1st_l,phi2y_add_1st_l,phi2yt_add_1st_l,&
     			phi2z_add_1st_l,phi2zt_add_1st_l)
!     ==================================================================
!
IMPLICIT NONE
!
!% INPUT VARIABLES
!
!	first-order additional time modal amplitudes
REAL(RP),DIMENSION(m3_add,m2),INTENT(IN) :: a1st_add_l,da1st_add_l
!
!	first-order additional quantities needed on the free surface
REAL(RP),DIMENSION(m1,m2),INTENT(OUT) :: phi2_1st,phi2t_1st,phi2x_1st,phi2y_1st,phi2z_1st,phi2zz_1st
REAL(RP),DIMENSION(m1,m2),INTENT(OUT) :: phi2xt_1st,phi2yt_1st,phi2zt_1st,phi2xx_1st,phi2xz_1st,phi2yz_1st
!
!	  on the wavemaker
REAL(RP),DIMENSION(m3,m2),INTENT(OUT) :: phi2xx_add_1st_l,phi2xxt_add_1st_l,phi2y_add_1st_l
REAL(RP),DIMENSION(m3,m2),INTENT(OUT) :: phi2yt_add_1st_l,phi2z_add_1st_l,phi2zt_add_1st_l
!
!% LOCAL VARIABLES
REAL(RP),DIMENSION(m3_add,m2) :: phi2xx_add_1st2,phi2xxt_add_1st2
REAL(RP),DIMENSION(m3_add,m2) :: phi2y_add_1st2,phi2yt_add_1st2
REAL(RP),DIMENSION(m3_add,m2) :: phi2z_add_1st2,phi2zt_add_1st2
!    
REAL(RP) :: ti, tf
INTEGER  :: i1, i2, i3, j2
!
REAL(RP) :: coskx_add_eta, sinkx_add_eta
!
!	CPU times inlet
if(iCPUtime.ge.1) then
   print*,'entering subroutine phiadd_1st_SWEET'
   call CPU_TIME(ti)
endif
!
! re-construction  on the free surface (z=0 because used for SWEET additional computation)
do i2=1,n2
   do i1=1,n1     
      phi2_1st(i1,i2)   = 0.0_rp
      phi2t_1st(i1,i2)  = 0.0_rp
      phi2x_1st(i1,i2)  = 0.0_rp
      phi2y_1st(i1,i2)  = 0.0_rp
      phi2z_1st(i1,i2)  = 0.0_rp
      phi2zz_1st(i1,i2) = 0.0_rp
      phi2xt_1st(i1,i2) = 0.0_rp
      phi2yt_1st(i1,i2) = 0.0_rp
      phi2zt_1st(i1,i2) = 0.0_rp
      phi2xx_1st(i1,i2) = 0.0_rp
      phi2xz_1st(i1,i2) = 0.0_rp
      phi2yz_1st(i1,i2) = 0.0_rp
      !
      do i3=1,n3_add
         do j2=1,n2
            coskx_add_eta = cos(kx_add(i3))
            sinkx_add_eta = sin(kx_add(i3))
            !
            phi2_1st(i1,i2)   = phi2_1st(i1,i2)   +  a1st_add_l(i3,j2) * coskx_add_eta * csh_add_x(i1,j2,i3) &
            	* cos(ky(j2)*cy(i2)) 
            phi2t_1st(i1,i2)  = phi2t_1st(i1,i2)  + da1st_add_l(i3,j2) * coskx_add_eta * csh_add_x(i1,j2,i3) &
            	* cos(ky(j2)*cy(i2)) 
            phi2x_1st(i1,i2)  = phi2x_1st(i1,i2)  -  a1st_add_l(i3,j2) * coskx_add_eta * k_add_sh_add_x(i1,j2,i3) &
            	* cos(ky(j2)*cy(i2)) 
            phi2y_1st(i1,i2)  = phi2y_1st(i1,i2)  -  a1st_add_l(i3,j2) * coskx_add_eta * kycsh_add_x(i1,j2,i3) &
            	* sin(ky(j2)*cy(i2))
            phi2z_1st(i1,i2)  = phi2z_1st(i1,i2)  -  a1st_add_l(i3,j2) * sinkx_add_eta * kx_add_csh_add_x(i1,j2,i3) &
            	* cos(ky(j2)*cy(i2)) 
            phi2zz_1st(i1,i2) = phi2zz_1st(i1,i2) -  a1st_add_l(i3,j2) * coskx_add_eta * kx_add_csh_add_x(i1,j2,i3) * kx_add(i3) &
                * cos(ky(j2)*cy(i2))
            phi2xt_1st(i1,i2) = phi2xt_1st(i1,i2) - da1st_add_l(i3,j2) * coskx_add_eta*  k_add_sh_add_x(i1,j2,i3) &
            	* cos(ky(j2)*cy(i2)) 
            phi2yt_1st(i1,i2) = phi2yt_1st(i1,i2) - da1st_add_l(i3,j2) * coskx_add_eta*  kycsh_add_x(i1,j2,i3) &
            	* sin(ky(j2)*cy(i2))
            phi2zt_1st(i1,i2) = phi2zt_1st(i1,i2) - da1st_add_l(i3,j2) * sinkx_add_eta*  kx_add_csh_add_x(i1,j2,i3) &
            	* cos(ky(j2)*cy(i2)) 
            phi2xx_1st(i1,i2) = phi2xx_1st(i1,i2) +  a1st_add_l(i3,j2) * coskx_add_eta * csh_add_x(i1,j2,i3)*k_add(i3,j2) &
            	* k_add(i3,j2) * cos(ky(j2)*cy(i2))
            phi2xz_1st(i1,i2) = phi2xz_1st(i1,i2) +  a1st_add_l(i3,j2) * sinkx_add_eta* k_add_sh_add_x(i1,j2,i3)* kx_add(i3) &
            	* cos(ky(j2)*cy(i2))
            phi2yz_1st(i1,i2) = phi2yz_1st(i1,i2) +  a1st_add_l(i3,j2) * sinkx_add_eta* ky(j2) * csh_add_x(i1,j2,i3) &
            	* sin(ky(j2)*cy(i2))
         end do
      enddo
   end do
end do
!
!     on the wavemaker
do i3=1,n3_add
   do i2=1,n2
      phi2xx_add_1st2(i3,i2)  =   a1st_add_l(i3,i2) * k_add(i3,i2)*k_add(i3,i2)
      phi2xxt_add_1st2(i3,i2) =  da1st_add_l(i3,i2) * k_add(i3,i2)*k_add(i3,i2)
      phi2y_add_1st2(i3,i2)   = - a1st_add_l(i3,i2) * ky(i2)
      phi2yt_add_1st2(i3,i2)  = -da1st_add_l(i3,i2) * ky(i2)
      phi2z_add_1st2(i3,i2)   = - a1st_add_l(i3,i2) * kx_add(i3)
      phi2zt_add_1st2(i3,i2)  = -da1st_add_l(i3,i2) * kx_add(i3)
   end do
end do
!
phi2xx_add_1st2  = Fourier_2_space_add(phi2xx_add_1st2,  'cos','cos')
phi2xxt_add_1st2 = Fourier_2_space_add(phi2xxt_add_1st2, 'cos','cos')
phi2y_add_1st2   = Fourier_2_space_add(phi2y_add_1st2,   'cos','sin')
phi2yt_add_1st2  = Fourier_2_space_add(phi2yt_add_1st2,  'cos','sin')
phi2z_add_1st2   = Fourier_2_space_add(phi2z_add_1st2,   'sin','cos')
phi2zt_add_1st2  = Fourier_2_space_add(phi2zt_add_1st2,  'sin','cos')
!
do i3=1,n3
   do i2=1,n2
      phi2xx_add_1st_l(i3,i2)  = phi2xx_add_1st2(i3,i2)
      phi2xxt_add_1st_l(i3,i2) = phi2xxt_add_1st2(i3,i2)
      phi2y_add_1st_l(i3,i2)   = phi2y_add_1st2(i3,i2)
      phi2yt_add_1st_l(i3,i2)  = phi2yt_add_1st2(i3,i2)
      phi2z_add_1st_l(i3,i2)   = phi2z_add_1st2(i3,i2)
      phi2zt_add_1st_l(i3,i2)  = phi2zt_add_1st2(i3,i2)
   end do
end do
!
!	CPU times outlet
if(iCPUtime.ge.1) then
   call CPU_TIME(tf)
   write(*,910)'quitting subroutine phiadd_1st_SWEET, total CPU time: ',tf-ti,'s'
endif
!
910 format(a,1ES11.4,a)
!
return
!
end subroutine phiadd_1st_SWEET
!
! end phiadd_1st ******************************************************
!
!
!
! start phidphiz ******************************************************
!
!     ============================================
      subroutine phidphiz(ark,phi_l,phiz_l,phi2_l,phi2z_l)
!     ============================================
!     
implicit none
!
!% INPUT VARIABLES				(at any of the two orders)
!
!	rk4 spectral time modal amplitudes
REAL(RP),DIMENSION(m1,m2),INTENT(IN) :: ark
!
!	spectral potential and its vertical derivative
REAL(RP),DIMENSION(m1,m2),INTENT(OUT) :: phi_l,phiz_l
!
!	additional potential and its vertical derivative
REAL(RP),DIMENSION(m1,m2),INTENT(IN) :: phi2_l,phi2z_l
!
REAL(RP) :: ti, tf
INTEGER  :: i1, i2
!
!	CPU times inlet
if(iCPUtime.ge.1) then
   print*,'entering subroutine phidphidz'
   call CPU_TIME(ti)
endif
!
!	re-construction on the free surface	
do i1=1,n1
   do i2=1,n2
      phi_l(i1,i2) = ark(i1,i2)
      phiz_l(i1,i2)= ark(i1,i2)*kth(i1,i2,1)
   end do
end do
!
phi_l  = Fourier_2_space(phi_l,'cos','cos')
phiz_l = Fourier_2_space(phiz_l,'cos','cos')
!
do i1=1,n1
   do i2=1,n2
      phi_l(i1,i2)  = phi_l(i1,i2)  + phi2_l(i1,i2)
      phiz_l(i1,i2) = phiz_l(i1,i2) + phi2z_l(i1,i2)
   end do
end do
!
!	CPU times outlet
if(iCPUtime.ge.1) then
   call CPU_TIME(tf)
   write(*,910)'quitting subroutine phidphiz, total CPU time: ',tf-ti,'s'
endif
!
910 format(a,1ES11.4,a)
!
return
!
end subroutine phidphiz
!
! end phidphiz ********************************************************
!
!
!
! start d1st_for_2nd **************************************************
!      
!     ==================================================================
      subroutine d1st_for_2nd(a1strk_loc,da1strk_loc,phi2x_1st_l,phi2y_1st_l,&
     			phi2zz_1st_l,phi2xt_1st_l,phi2yt_1st_l,phi2zt_1st_l,phi2xx_1st_l&
     			,phi2xz_1st_l,phi2yz_1st_l,phi2xx_add_1st_l,&
     			phi2xxt_add_1st_l,phi2y_add_1st_l,phi2yt_add_1st_l,phi2z_add_1st_l&
                        ,phi2zt_add_1st_l)
!     ==================================================================
!
IMPLICIT NONE
!
!% INPUT VARIABLES
!
!	first-order additional time modal amplitudes
!
REAL(RP),DIMENSION(m1,m2),INTENT(IN) :: a1strk_loc,da1strk_loc
!
!	first-order additional quantities needed
!	  on the free surface
REAL(RP),DIMENSION(m1,m2),INTENT(IN) :: phi2x_1st_l,phi2y_1st_l,phi2zz_1st_l,phi2xt_1st_l,phi2yt_1st_l
REAL(RP),DIMENSION(m1,m2),INTENT(IN) :: phi2xx_1st_l,phi2xz_1st_l,phi2yz_1st_l,phi2zt_1st_l
!
!	  on the wavemaker
REAL(RP),DIMENSION(m3,m2),INTENT(IN) :: phi2xx_add_1st_l,phi2xxt_add_1st_l,phi2y_add_1st_l,phi2yt_add_1st_l
REAL(RP),DIMENSION(m3,m2),INTENT(IN) :: phi2z_add_1st_l,phi2zt_add_1st_l
!
!% LOCAL VARIABLES
!
!	first-order additional quantities needed for spatial derivatives of eta
!	  on the free surface
!
REAL(RP),DIMENSION(m1,m2) :: phixz1st,phiyz1st
!
REAL(RP) :: ti, tf
INTEGER  :: i1, i2, i3
!
!	CPU times inlet
if(iCPUtime.ge.1) then
   print*,'entering subroutine d1st_for_2nd'
   call CPU_TIME(ti)
endif
!
IF(coeffilt(1).LT.1.d0) THEN
   CALL filtering_x_modes(a1strk_loc,coeffilt(1))
   CALL filtering_x_modes(da1strk_loc,coeffilt(1))
ENDIF
IF (n2 /= 1 .AND. coeffilt(2) < 1.0_rp) THEN
   CALL filtering_y_modes(a1strk_loc,coeffilt(2))
   CALL filtering_y_modes(da1strk_loc,coeffilt(2))
END IF
!
!	re-construction on the wavemaker
!
do i2=1,n2
   do i3=1,n3
      phixx1st_add(i3,i2)  = 0.0_rp
      phixxt1st_add(i3,i2) = 0.0_rp
      phiy1st_add(i3,i2)   = 0.0_rp
      phiyt1st_add(i3,i2)  = 0.0_rp
      phiz1st_add(i3,i2)   = 0.0_rp
      phizt1st_add(i3,i2)  = 0.0_rp
      do i1=1,n1
         phixx1st_add(i3,i2)  = phixx1st_add(i3,i2)  -  a1strk_loc(i1,i2) * kx2cshx_add(i1,i2,i3)
         phixxt1st_add(i3,i2) = phixxt1st_add(i3,i2) - da1strk_loc(i1,i2) * kx2cshx_add(i1,i2,i3)
         phiy1st_add(i3,i2)   = phiy1st_add(i3,i2)   -  a1strk_loc(i1,i2) * kycshx_add(i1,i2,i3)
         phiyt1st_add(i3,i2)  = phiyt1st_add(i3,i2)  - da1strk_loc(i1,i2) * kycshx_add(i1,i2,i3)
         phiz1st_add(i3,i2)   = phiz1st_add(i3,i2)   +  a1strk_loc(i1,i2) * kshx_add(i1,i2,i3)
         phizt1st_add(i3,i2)  = phizt1st_add(i3,i2)  + da1strk_loc(i1,i2) * kshx_add(i1,i2,i3)
      end do
   end do
end do
!
do i3=1,n3
   phixx1st_add(i3,1:n2)  = Fourier_2_space_y(phixx1st_add(i3,1:n2), 'cos')
   phixxt1st_add(i3,1:n2) = Fourier_2_space_y(phixxt1st_add(i3,1:n2),'cos')
   phiy1st_add(i3,1:n2)   = Fourier_2_space_y(phiy1st_add(i3,1:n2),  'sin')
   phiyt1st_add(i3,1:n2)  = Fourier_2_space_y(phiyt1st_add(i3,1:n2), 'sin')
   phiz1st_add(i3,1:n2)   = Fourier_2_space_y(phiz1st_add(i3,1:n2),  'cos')
   phizt1st_add(i3,1:n2)  = Fourier_2_space_y(phizt1st_add(i3,1:n2), 'cos')
enddo
!
do i3=1,n3
   do i2=1,n2
      phixx1st_add(i3,i2)  = phixx1st_add(i3,i2)  + phi2xx_add_1st_l(i3,i2)
      phixxt1st_add(i3,i2) = phixxt1st_add(i3,i2) + phi2xxt_add_1st_l(i3,i2)
      phiy1st_add(i3,i2)   = phiy1st_add(i3,i2)   + phi2y_add_1st_l(i3,i2)
      phiyt1st_add(i3,i2)  = phiyt1st_add(i3,i2)  + phi2yt_add_1st_l(i3,i2)
      phiz1st_add(i3,i2)   = phiz1st_add(i3,i2)   + phi2z_add_1st_l(i3,i2)
      phizt1st_add(i3,i2)  = phizt1st_add(i3,i2)  + phi2zt_add_1st_l(i3,i2)
   end do
end do
!
!	re-construction on the free surface
do i1=1,n1
   do i2=1,n2
      phix1st(i1,i2)  = - a1strk_loc(i1,i2) * kx(i1)
      phiy1st(i1,i2)  = - a1strk_loc(i1,i2) * ky(i2)
      etax1st(i1,i2)  =  da1strk_loc(i1,i2) * kx(i1)
      etay1st(i1,i2)  =  da1strk_loc(i1,i2) * ky(i2)
      phizt1st(i1,i2) =  da1strk_loc(i1,i2) * kth(i1,i2,1)
      phixx1st(i1,i2) = - a1strk_loc(i1,i2) * kx2(i1)
      phizz1st(i1,i2) =   a1strk_loc(i1,i2) * k(i1,i2)*k(i1,i2)
      phixz1st(i1,i2) = - a1strk_loc(i1,i2) * kx(i1) * kth(i1,i2,1)
      phiyz1st(i1,i2) = phiy1st(i1,i2) * kth(i1,i2,1)
   end do
end do
!
phix1st = Fourier_2_space(phix1st,   'sin','cos')
phiy1st = Fourier_2_space(phiy1st,   'cos','sin')
etax1st = Fourier_2_space(etax1st,   'sin','cos')
etay1st = Fourier_2_space(etay1st,   'cos','sin')
phizt1st = Fourier_2_space(phizt1st, 'cos','cos')
phixx1st = Fourier_2_space(phixx1st, 'cos','cos')
phizz1st = Fourier_2_space(phizz1st, 'cos','cos')
phixz1st = Fourier_2_space(phixz1st, 'sin','cos')
phiyz1st = Fourier_2_space(phiyz1st, 'cos','sin')
!
do i1=1,n1
   do i2=1,n2
      phix1st(i1,i2)  = phix1st(i1,i2)  + phi2x_1st_l(i1,i2)
      phiy1st(i1,i2)  = phiy1st(i1,i2)  + phi2y_1st_l(i1,i2)
      phizt1st(i1,i2) = phizt1st(i1,i2) + phi2zt_1st_l(i1,i2)
      phixx1st(i1,i2) = phixx1st(i1,i2) + phi2xx_1st_l(i1,i2)
      phizz1st(i1,i2) = phizz1st(i1,i2) + phi2zz_1st_l(i1,i2)
      phixz1st(i1,i2) = phixz1st(i1,i2) + phi2xz_1st_l(i1,i2)
      phiyz1st(i1,i2) = phiyz1st(i1,i2) + phi2yz_1st_l(i1,i2)
   end do
end do
!
!	re-construction of FS elevation derivatives including the absorption
do i1=1,n1
   do i2=1,n2
      etax1st(i1,i2)=etax1st(i1,i2)-phi2xt_1st_l(i1,i2)&
           -iabsnb*(dnux(i1,i2)*phiz1st(i1,i2)+nu(i1,i2)*phixz1st(i1,i2)) 
      !
      etay1st(i1,i2)=etay1st(i1,i2)-phi2yt_1st_l(i1,i2)&
           -iabsnb*(dnuy(i1,i2)*phiz1st(i1,i2)+nu(i1,i2)*phiyz1st(i1,i2))   								
   end do
end do
!
!	CPU times outlet
if(iCPUtime.ge.1) then
   call CPU_TIME(tf)
   write(*,910)'quitting subroutine d1st_for_2nd, total CPU time: ',tf-ti,'s'
endif
!
910 format(a,1ES11.4,a)
!
return
!
end subroutine d1st_for_2nd
!
! end d1st_for_2nd ****************************************************
!
!
!
! start phiadd_2nd ****************************************************
!      
!     ============================================================
      subroutine phiadd_2nd(a2nd_add_rk_l,da2nd_add_rk_l, phi2_2nd, phi2x_2nd, phi2y_2nd, phi2t_2nd, phi2z_2nd)
!     ============================================================
!
IMPLICIT NONE
!
!% INPUT VARIABLES
!
!	second-order additional time modal amplitudes
REAL(RP),DIMENSION(m3_add,m2),INTENT(IN)  :: a2nd_add_rk_l,da2nd_add_rk_l
!
!	second-order additional quantities needed (on the free surface)
REAL(RP),DIMENSION(m1,m2),INTENT(OUT) :: phi2_2nd,phi2x_2nd,phi2y_2nd,phi2t_2nd,phi2z_2nd
!
REAL(RP) :: ti, tf
INTEGER  :: i1, i2, i3, j2
!
REAL(RP) ::  coskx_add_eta, sinkx_add_eta
!
!	CPU times inlet
!
if (iCPUtime.ge.1) then
   print*,'entering subroutine phiadd_2nd'
   call CPU_TIME(ti)
endif
!
!	re-construction on the free surface (z=0, part of additional computation i.e. SWEET equivalent)	
do i2=1,n2
   do i1=1,n1
      phi2_2nd(i1,i2)  = 0.d0
      phi2t_2nd(i1,i2) = 0.d0
      phi2z_2nd(i1,i2) = 0.d0
      phi2x_2nd(i1,i2) = 0.d0
      phi2y_2nd(i1,i2) = 0.d0
      do i3=1,n3_add
         do j2=1,n2
            coskx_add_eta = cos(kx_add(i3))
            sinkx_add_eta = sin(kx_add(i3))
            phi2_2nd(i1,i2)  = phi2_2nd(i1,i2)  + a2nd_add_rk_l(i3,j2)* coskx_add_eta * csh_add_x(i1,j2,i3) &
            	* cos(ky(j2)*cy(i2))  
            phi2t_2nd(i1,i2) = phi2t_2nd(i1,i2) + da2nd_add_rk_l(i3,j2)* coskx_add_eta * csh_add_x(i1,j2,i3) &
            	* cos(ky(j2)*cy(i2)) 
            phi2z_2nd(i1,i2) = phi2z_2nd(i1,i2) - a2nd_add_rk_l(i3,j2) * kx_add(i3) *sinkx_add_eta * csh_add_x(i1,j2,i3) &
            	* cos(ky(j2)*cy(i2))
            phi2x_2nd(i1,i2) = phi2x_2nd(i1,i2) - a2nd_add_rk_l(i3,j2)* coskx_add_eta * k_add_sh_add_x(i1,j2,i3) &
            	* cos(ky(j2)*cy(i2)) 
            phi2y_2nd(i1,i2) = phi2y_2nd(i1,i2) - a2nd_add_rk_l(i3,j2) * coskx_add_eta * kycsh_add_x(i1,j2,i3) &
            	* sin(ky(j2)*cy(i2)) 
         end do
      end do
   end do
end do
!
!	CPU times outlet
if (iCPUtime.ge.1) then
   call CPU_TIME(tf)
   write(*,910)'quitting subroutine phiadd_2nd, total CPU time: ',tf-ti,'s'
endif
!
910 format(a,1ES11.4,a)
!
return
!
end subroutine phiadd_2nd
!
! end phiadd_2nd ******************************************************
!
!
! start phiadd_3rd ****************************************************
!      
!     ============================================================
      subroutine phiadd_3rd(a3rd_add_rk_l,da3rd_add_rk_l, phi2_3rd, phi2x_3rd, phi2y_3rd, phi2t_3rd, phi2z_3rd)
!     ============================================================
!
IMPLICIT NONE
!
!% INPUT VARIABLES
!
!	second-order additional time modal amplitudes
REAL(RP),DIMENSION(m3_add,m2),INTENT(IN) :: a3rd_add_rk_l,da3rd_add_rk_l
!
!	second-order additional quantities needed (on the free surface)
REAL(RP),DIMENSION(m1,m2),INTENT(OUT) :: phi2_3rd,phi2x_3rd,phi2y_3rd,phi2t_3rd,phi2z_3rd
!
REAL(RP) :: ti, tf
INTEGER  :: i1, i2, i3, j2
!
REAL(RP) ::  coskx_add_eta, sinkx_add_eta
!
!	CPU times inlet
if (iCPUtime.ge.1) then
   print*,'entering subroutine phiadd_3rd'
   call CPU_TIME(ti)
endif
!
!	re-construction on the free surface (z=0, part of additional computation i.e. SWEET equivalent)		
do i2=1,n2
   do i1=1,n1
      phi2_3rd(i1,i2)  = 0.d0
      phi2t_3rd(i1,i2) = 0.d0
      phi2z_3rd(i1,i2) = 0.d0
      phi2x_3rd(i1,i2) = 0.d0
      phi2y_3rd(i1,i2) = 0.d0
      do i3=1,n3_add
         do j2=1,n2
            coskx_add_eta = cos(kx_add(i3))
            sinkx_add_eta = sin(kx_add(i3))
            phi2_3rd(i1,i2)  = phi2_3rd(i1,i2)  + a3rd_add_rk_l(i3,j2)* coskx_add_eta * csh_add_x(i1,j2,i3) &
            	* cos(ky(j2)*cy(i2))  
            phi2t_3rd(i1,i2) = phi2t_3rd(i1,i2) + da3rd_add_rk_l(i3,j2)* coskx_add_eta * csh_add_x(i1,j2,i3) &
            	* cos(ky(j2)*cy(i2)) 
            phi2z_3rd(i1,i2) = phi2z_3rd(i1,i2) - a3rd_add_rk_l(i3,j2) * kx_add(i3) *sinkx_add_eta * csh_add_x(i1,j2,i3) &
            	* cos(ky(j2)*cy(i2))
            phi2x_3rd(i1,i2) = phi2x_3rd(i1,i2) - a3rd_add_rk_l(i3,j2) * coskx_add_eta * k_add_sh_add_x(i1,j2,i3) &
            	* cos(ky(j2)*cy(i2)) 
            phi2y_3rd(i1,i2) = phi2y_3rd(i1,i2) - a3rd_add_rk_l(i3,j2) * coskx_add_eta * kycsh_add_x(i1,j2,i3) &
            	* sin(ky(j2)*cy(i2))
         end do
      end do
   end do
end do
!
!	CPU times outlet
if (iCPUtime.ge.1) then
   call CPU_TIME(tf)
   write(*,910)'quitting subroutine phiadd_3rd, total CPU time: ',tf-ti,'s'
endif
!
910 format(a,1ES11.4,a)
!
return
!
end subroutine phiadd_3rd
!
! end phiadd_3rd ******************************************************
!
!
!
! start wmk_3rdorder ****************************************************
!      
!     ==================================================================
      subroutine wmk_3rdorder(a1strk_loc,da1strk_loc,a2ndrk_loc,da2ndrk_loc,a1st_add_l,da1st_add_l,a2nd_add_l,da2nd_add_l)
!     ==================================================================
!     
implicit none
!
!% INPUT VARIABLES
!
!	first-order additional time modal amplitudes
REAL(RP),DIMENSION(m3_add,m2),INTENT(IN) :: a1st_add_l,da1st_add_l,a2nd_add_l,da2nd_add_l
!
!       first order amplitudes
REAL(RP),DIMENSION(m1,m2),INTENT(IN) :: a1strk_loc,da1strk_loc,a2ndrk_loc,da2ndrk_loc
!
!% LOCAL VARIABLES
!
!% Local variables needed by third order
REAL(RP),DIMENSION(m3_add,m2) :: phi2xxx_add_1st2,phi2xxxt_add_1st2,phi2xy_add_1st2
REAL(RP),DIMENSION(m3_add,m2) :: phi2xyt_add_1st2,phi2xz_add_1st2,phi2xzt_add_1st2
!
REAL(RP),DIMENSION(m3_add,m2) :: phi2xx_add_2nd2,phi2xxt_add_2nd2,phi2y_add_2nd2
REAL(RP),DIMENSION(m3_add,m2) :: phi2yt_add_2nd2,phi2z_add_2nd2,phi2zt_add_2nd2
!
REAL(RP),DIMENSION(m3_add,m2) :: phi_add_1st,phit_add_1st,phi_add_2nd,phit_add_2nd
! 
REAL(RP) :: ti, tf
INTEGER  :: i1, i2, i3
!
!	CPU times inlet
if(iCPUtime.ge.1) then
   print*,'entering subroutine wmk_3rdorder'
   call CPU_TIME(ti)
endif
!
!     on the wavemaker
do i3=1,n3_add
   do i2=1,n2
      !
      ! First order
      phi2xxx_add_1st2(i3,i2) =  -a1st_add_l(i3,i2) * k_add(i3,i2)*k_add(i3,i2)*k_add_thk_add(i3,i2)
      phi2xxxt_add_1st2(i3,i2)= -da1st_add_l(i3,i2) * k_add(i3,i2)*k_add(i3,i2)*k_add_thk_add(i3,i2)
      phi2xy_add_1st2(i3,i2)  =   a1st_add_l(i3,i2) * ky(i2)*k_add_thk_add(i3,i2)
      phi2xyt_add_1st2(i3,i2) =  da1st_add_l(i3,i2) * ky(i2)*k_add_thk_add(i3,i2)
      phi2xz_add_1st2(i3,i2)  =   a1st_add_l(i3,i2) * kx_add(i3)*k_add_thk_add(i3,i2)
      phi2xzt_add_1st2(i3,i2) =  da1st_add_l(i3,i2) * kx_add(i3)*k_add_thk_add(i3,i2)
      !
      ! Second order
      phi2xx_add_2nd2(i3,i2) =   a2nd_add_l(i3,i2) * k_add(i3,i2)*k_add(i3,i2)
      phi2xxt_add_2nd2(i3,i2)=  da2nd_add_l(i3,i2) * k_add(i3,i2)*k_add(i3,i2)
      phi2y_add_2nd2(i3,i2)  =  -a2nd_add_l(i3,i2) * ky(i2)
      phi2yt_add_2nd2(i3,i2) = -da2nd_add_l(i3,i2) * ky(i2)
      phi2z_add_2nd2(i3,i2)  =  -a2nd_add_l(i3,i2) * kx_add(i3)
      phi2zt_add_2nd2(i3,i2) = -da2nd_add_l(i3,i2) * kx_add(i3)
      !
      phi_add_1st(i3,i2) =  a1st_add_l(i3,i2) 
      phit_add_1st(i3,i2)= da1st_add_l(i3,i2)
      phi_add_2nd(i3,i2) =  a2nd_add_l(i3,i2)
      phit_add_2nd(i3,i2)= da2nd_add_l(i3,i2)
   end do
end do
!
! First order
phi2xxx_add_1st2  = Fourier_2_space_add(phi2xxx_add_1st2, 'cos','cos')
phi2xxxt_add_1st2 = Fourier_2_space_add(phi2xxxt_add_1st2,'cos','cos')
phi2xy_add_1st2   = Fourier_2_space_add(phi2xy_add_1st2,  'cos','sin')
phi2xyt_add_1st2  = Fourier_2_space_add(phi2xyt_add_1st2, 'cos','sin')
phi2xz_add_1st2   = Fourier_2_space_add(phi2xz_add_1st2,  'sin','cos')
phi2xzt_add_1st2  = Fourier_2_space_add(phi2xzt_add_1st2, 'sin','cos')
!
! Second order
phi2xx_add_2nd2 = Fourier_2_space_add(phi2xx_add_2nd2, 'cos','cos')
phi2xxt_add_2nd2 = Fourier_2_space_add(phi2xxt_add_2nd2,'cos','cos')
phi2y_add_2nd2 = Fourier_2_space_add(phi2y_add_2nd2,  'cos','sin')
phi2yt_add_2nd2 = Fourier_2_space_add(phi2yt_add_2nd2, 'cos','sin')
phi2z_add_2nd2 = Fourier_2_space_add(phi2z_add_2nd2,  'sin','cos')
phi2zt_add_2nd2 = Fourier_2_space_add(phi2zt_add_2nd2, 'sin','cos')
!
phi_add_1st = Fourier_2_space_add(phi_add_1st, 'cos','cos')
phit_add_1st = Fourier_2_space_add(phit_add_1st,'cos','cos')
phi_add_2nd = Fourier_2_space_add(phi_add_2nd, 'cos','cos')
phit_add_2nd = Fourier_2_space_add(phit_add_2nd,'cos','cos')
!
! Filtering if necessary
IF(coeffilt(1).LT.1.d0) THEN
   CALL filtering_x_modes(a1strk_loc,coeffilt(1))
   CALL filtering_x_modes(da1strk_loc,coeffilt(1))
   CALL filtering_x_modes(a2ndrk_loc,coeffilt(1))
   CALL filtering_x_modes(da2ndrk_loc,coeffilt(1))
ENDIF
IF (n2 /= 1 .AND. coeffilt(2) < 1.0_rp) THEN
   CALL filtering_y_modes(a1strk_loc,coeffilt(2))
   CALL filtering_y_modes(da1strk_loc,coeffilt(2))
   CALL filtering_y_modes(a2ndrk_loc,coeffilt(2))
   CALL filtering_y_modes(da2ndrk_loc,coeffilt(2))
END IF
!
!	re-construction on the wavemaker
do i2=1,n2
   do i3=1,n3
      !
      ! First order
      phixxx1st_add(i3,i2)  = 0.0_rp
      phixxxt1st_add(i3,i2) = 0.0_rp
      phixy1st_add(i3,i2)   = 0.0_rp
      phixyt1st_add(i3,i2)  = 0.0_rp
      phixz1st_add(i3,i2)   = 0.0_rp
      phixzt1st_add(i3,i2)  = 0.0_rp
      !
      ! Second order
      phixx2nd_add(i3,i2)  = 0.0_rp
      phixxt2nd_add(i3,i2) = 0.0_rp
      phiy2nd_add(i3,i2)   = 0.0_rp
      phiyt2nd_add(i3,i2)  = 0.0_rp
      phiz2nd_add(i3,i2)   = 0.0_rp
      phizt2nd_add(i3,i2)  = 0.0_rp
      !
      do i1=1,n1
         !
         ! First order equal to zero on x=0 (sinus(x) serie)
         !
         ! Second order
         phixx2nd_add(i3,i2)  = phixx2nd_add(i3,i2)  -  a2ndrk_loc(i1,i2) * kx2cshx_add(i1,i2,i3)
         phixxt2nd_add(i3,i2) = phixxt2nd_add(i3,i2) - da2ndrk_loc(i1,i2) * kx2cshx_add(i1,i2,i3)
         phiy2nd_add(i3,i2)   = phiy2nd_add(i3,i2)   -  a2ndrk_loc(i1,i2) * kycshx_add(i1,i2,i3)
         phiyt2nd_add(i3,i2)  = phiyt2nd_add(i3,i2)  - da2ndrk_loc(i1,i2) * kycshx_add(i1,i2,i3)
         phiz2nd_add(i3,i2)   = phiz2nd_add(i3,i2)   +  a2ndrk_loc(i1,i2) * kshx_add(i1,i2,i3)
         phizt2nd_add(i3,i2)  = phizt2nd_add(i3,i2)  + da2ndrk_loc(i1,i2) * kshx_add(i1,i2,i3)
      end do
   end do
end do
!
! First order equal to zero on x=0 (sinus(x) serie)
!
! Second order
do i3=1,n3
   phixx2nd_add(i3,1:n2)  = Fourier_2_space_y(phixx2nd_add(i3,1:n2), 'cos')
   phixxt2nd_add(i3,1:n2) = Fourier_2_space_y(phixxt2nd_add(i3,1:n2),'cos')
   phiy2nd_add(i3,1:n2)   = Fourier_2_space_y(phiy2nd_add(i3,1:n2),  'sin')
   phiyt2nd_add(i3,1:n2)  = Fourier_2_space_y(phiyt2nd_add(i3,1:n2), 'sin')
   phiz2nd_add(i3,1:n2)   = Fourier_2_space_y(phiz2nd_add(i3,1:n2),  'cos')
   phizt2nd_add(i3,1:n2)  = Fourier_2_space_y(phizt2nd_add(i3,1:n2), 'cos')
enddo
!
do i3=1,n3
   do i2=1,n2
      !
      ! First order
      phixxx1st_add(i3,i2)  = phixxx1st_add(i3,i2)  + phi2xxx_add_1st2(i3,i2)
      phixxxt1st_add(i3,i2) = phixxxt1st_add(i3,i2) + phi2xxxt_add_1st2(i3,i2)
      phixy1st_add(i3,i2)   = phixy1st_add(i3,i2)   + phi2xy_add_1st2(i3,i2)
      phixyt1st_add(i3,i2)  = phixyt1st_add(i3,i2)  + phi2xyt_add_1st2(i3,i2)
      phixz1st_add(i3,i2)   = phixz1st_add(i3,i2)   + phi2xz_add_1st2(i3,i2)
      phixzt1st_add(i3,i2)  = phixzt1st_add(i3,i2)  + phi2xzt_add_1st2(i3,i2)
      !
      ! Second order
      phixx2nd_add(i3,i2)  = phixx2nd_add(i3,i2)  + phi2xx_add_2nd2(i3,i2)
      phixxt2nd_add(i3,i2) = phixxt2nd_add(i3,i2) + phi2xxt_add_2nd2(i3,i2)
      phiy2nd_add(i3,i2)   = phiy2nd_add(i3,i2)   + phi2y_add_2nd2(i3,i2)
      phiyt2nd_add(i3,i2)  = phiyt2nd_add(i3,i2)  + phi2yt_add_2nd2(i3,i2)
      phiz2nd_add(i3,i2)   = phiz2nd_add(i3,i2)   + phi2z_add_2nd2(i3,i2)
      phizt2nd_add(i3,i2)  = phizt2nd_add(i3,i2)  + phi2zt_add_2nd2(i3,i2)
   end do
end do
!
!	CPU times outlet
if(iCPUtime.ge.1) then
   call CPU_TIME(tf)
   write(*,910)'quitting subroutine wmk_3rdorder, total CPU time: ',tf-ti,'s'
endif
!
910 format(a,1ES11.4,a)
!
return
!
end subroutine wmk_3rdorder
!
! end wmk_3rdorder ******************************************************
!
!
!
! start phiadd_SL ****************************************************
!      
!     ==================================================================
      subroutine phiadd_SL(a_add_l,da_add_l,phi_add_l,phix_add_l,phiy_add_l,phiz_add_l,phit_add_l, eta_l)
!     ==================================================================
!      
implicit none
!
!% INPUT VARIABLES
!
!	additional time modal amplitudes
REAL(RP),DIMENSION(m3_add,m2),INTENT(IN)  :: a_add_l,da_add_l
REAL(RP),DIMENSION(m1,m2),INTENT(IN)      :: eta_l
!
!	additional quantities needed on the free surface
REAL(RP),DIMENSION(m1,m2),INTENT(OUT) :: phi_add_l,phix_add_l,phiy_add_l,phiz_add_l,phit_add_l
!    
REAL(RP) :: ti, tf
INTEGER  :: i1, i2, i3, j2
REAL(RP) :: kx_add_etam1, coskx_add_eta, sinkx_add_eta
!
!	CPU times inlet
if(iCPUtime.ge.1) then
   print*,'entering subroutine phiadd_SL'
   call CPU_TIME(ti)
endif
!
!	re-construction  on the free surface (z=eta for HOS-NWT)
do i2=1,n2
   do i1=1,n1     
      phi_add_l(i1,i2)   = 0.0_rp
      phit_add_l(i1,i2)  = 0.0_rp
      phix_add_l(i1,i2)  = 0.0_rp
      phiy_add_l(i1,i2)  = 0.0_rp
      phiz_add_l(i1,i2)  = 0.0_rp
      do i3=1,n3_add
         do j2=1,n2
            kx_add_etam1  = kx_add(i3) * (eta_l(i1,i2)+1.0_rp)
            coskx_add_eta = cos(kx_add_etam1)
            sinkx_add_eta = sin(kx_add_etam1)
            !
            phi_add_l(i1,i2)   = phi_add_l(i1,i2)   +  a_add_l(i3,j2) * coskx_add_eta * csh_add_x(i1,j2,i3) &
            	* cos(ky(j2)*cy(i2)) 
            phit_add_l(i1,i2)  = phit_add_l(i1,i2)  + da_add_l(i3,j2) * coskx_add_eta * csh_add_x(i1,j2,i3) &
            	* cos(ky(j2)*cy(i2)) 
            phix_add_l(i1,i2)  = phix_add_l(i1,i2)  -  a_add_l(i3,j2) * coskx_add_eta * k_add_sh_add_x(i1,j2,i3) &
            	* cos(ky(j2)*cy(i2)) 
            phiy_add_l(i1,i2)  = phiy_add_l(i1,i2)  -  a_add_l(i3,j2) * coskx_add_eta * kycsh_add_x(i1,j2,i3) &
            	* sin(ky(j2)*cy(i2))
            phiz_add_l(i1,i2)  = phiz_add_l(i1,i2)  -  a_add_l(i3,j2) * sinkx_add_eta * kx_add_csh_add_x(i1,j2,i3) &
            	* cos(ky(j2)*cy(i2)) 
         end do
      enddo
   end do
end do
!
!	CPU times outlet
if(iCPUtime.ge.1) then
   call CPU_TIME(tf)
   write(*,910)'quitting subroutine phiadd_SL, total CPU time: ',tf-ti,'s'
endif
!
910 format(a,1ES11.4,a)
!
return
!
end subroutine phiadd_SL
!
! end phiadd_SL ******************************************************
!
!
!
! start phiadd_1st_bis ****************************************************
!      
!     ==================================================================
      subroutine phiadd_1st_bis(a1st_add_l,da1st_add_l,phi2_1st,phi2t_1st,phi2x_1st,&
     			phi2y_1st,phi2z_1st, eta_l)
!     ==================================================================
!      
implicit none
!
!% INPUT VARIABLES
!	first-order additional time modal amplitudes
REAL(RP),DIMENSION(m3_add,m2),INTENT(IN) :: a1st_add_l,da1st_add_l
!
!	first-order additional quantities needed
!	  on the free surface
REAL(RP),DIMENSION(m1,m2),INTENT(IN) :: eta_l
REAL(RP),DIMENSION(m1,m2),INTENT(OUT) :: phi2_1st,phi2t_1st,phi2x_1st,phi2y_1st,phi2z_1st
!
!% LOCAL VARIABLES
REAL(RP) :: ti, tf
INTEGER  :: i1, i2, i3, j2
REAL(RP) :: kx_add_etam1, coskx_add_eta, sinkx_add_eta 
!
!	CPU times inlet
if (iCPUtime.eq.1) then
   print*,'entering subroutine phiadd_1st'
   call CPU_TIME(ti)
endif
!
!	re-construction on the free surface (z=eta)
do i2=1,n2
   do i1=1,n1
      phi2_1st(i1,i2)   = 0.0_rp
      phi2t_1st(i1,i2)  = 0.0_rp
      phi2x_1st(i1,i2)  = 0.0_rp
      phi2y_1st(i1,i2)  = 0.0_rp
      phi2z_1st(i1,i2)  = 0.0_rp
      do i3=1,n3_add
         do j2=1,n2
            kx_add_etam1  = kx_add(i3) * (eta_l(i1,i2)+1.0_rp)
            coskx_add_eta = cos(kx_add_etam1)
            sinkx_add_eta = sin(kx_add_etam1)
            !
            phi2_1st(i1,i2)   = phi2_1st(i1,i2)   +  a1st_add_l(i3,j2) * coskx_add_eta * csh_add_x(i1,j2,i3) &
            	* cos(ky(j2)*cy(i2))
            phi2t_1st(i1,i2)  = phi2t_1st(i1,i2)  + da1st_add_l(i3,j2) * coskx_add_eta * csh_add_x(i1,j2,i3) &
            	* cos(ky(j2)*cy(i2)) 
            phi2x_1st(i1,i2)  = phi2x_1st(i1,i2)  -  a1st_add_l(i3,j2) * coskx_add_eta * k_add_sh_add_x(i1,j2,i3) &
            	* cos(ky(j2)*cy(i2))
            phi2y_1st(i1,i2)  = phi2y_1st(i1,i2)  -  a1st_add_l(i3,j2) * coskx_add_eta * kycsh_add_x(i1,j2,i3) &
            	* sin(ky(j2)*cy(i2))
            phi2z_1st(i1,i2)  = phi2z_1st(i1,i2)  -  a1st_add_l(i3,j2) * sinkx_add_eta * kx_add_csh_add_x(i1,j2,i3) &
            	* cos(ky(j2)*cy(i2))
         end do
      enddo
   end do
end do
!
!	CPU times outlet
if (iCPUtime.eq.1) then
   call CPU_TIME(tf)
   write(*,910)'quitting subroutine phiadd_1st, total CPU time: ',tf-ti,'s'
endif
!
910 format(a,1ES11.4,a)
!
end subroutine phiadd_1st_bis
!
! end phiadd_1st_bis ******************************************************
!
!
!
! start ramp **********************************************************
!           
!     ===============
subroutine ramp
!     ===============
!      
IMPLICIT NONE
!
REAL(RP) :: T_swtt
REAL(RP) :: tmp, T2, tim2, timsup, tmp2, rap, rap2, rap22, t_inter
!
!  ramp duration checkings
if (Tramp > T_stop) write(*,*) 'overlaping ramps'
!
!	ramp time coefficients
IF (iramp == 11) THEN
   T_swtt = 3.0_rp / sqrt(5.0_rp / 9.81_rp)
   IF (time <= T_swtt) THEN
      alpha   = time   / T_swtt
      alphat  = 1.0_rp / T_swtt
      alphatt = 0.0_rp
   ELSE
      alpha   = 1.0_rp
      alphat  = 0.0_rp
      alphatt = 0.0_rp
   END IF
ELSE IF (iramp == 33) THEN
   T_swtt = - log(0.1d0) / (10.d0 * T_mono)
   alpha   = 1.d0 - exp( - T_swtt * time)
   alphat  = (1.d0 - alpha) * T_swtt
   alphatt = (1.d0 - alpha) * T_swtt * T_swtt
ELSE
   if (time.le.Tramp) then
      if (iramp.eq.0) then
         alpha   = 1.0d0
         alphat  = 0
         alphatt = 0
      elseif (iramp.eq.1) then
         alpha=time/Tramp
         alphat=1.d0/Tramp
         alphatt=0.d0
      elseif (iramp.eq.2) then
         tmp     = time/Tramp 
         alpha   = tmp * tmp
         alphat  = 2.0d0 * tmp / Tramp
         alphatt = 2.0d0 / (Tramp * Tramp)
      elseif (iramp.eq.21) then
         T2=Tramp**2
         tim2=time**2
         timsup=time/Tramp-1.d0
         if(time.le.(Tramp/2.d0)) then
            alpha=2.d0*tim2/T2
            alphat=4.d0*time/T2
            alphatt=4.d0/T2
         else
            alpha=1.d0-2.d0*timsup**2
            alphat=-4.d0*timsup/Tramp
            alphatt=-4.d0/T2
         endif
      elseif (iramp.eq.4) then
         tmp     = time/Tramp
         tmp2    = tmp * tmp
         alpha   = tmp2 * tmp2
         alphat  = 4.0d0 * tmp2 * tmp / Tramp
         alphatt = 12.0d0 * tmp2 / (Tramp * Tramp)
      elseif (iramp.eq.41) then
         T2=Tramp**2
         tim2=time**2
         rap=8.d0*time/Tramp
         rap2=tim2/T2
         rap22=3.d0*rap2
         alpha=(rap22-rap+6.d0)*rap2
         alphat=(4.d0*rap22-3.d0*rap+12.d0)*time/T2
         alphatt=(12.d0*rap22-6.d0*rap+12.d0)/T2
      else
         write(*,*) 'Unknown ramp'
      endif
   else 
	 if (time .ge. T_stop + Tramp) then
		alpha=0.d0
		alphat=0.d0
		alphatt=0.d0
	 else if (time .ge. T_stop) then
		t_inter = T_stop + Tramp - time
		if (iramp.eq.1) then
		   alpha=t_inter/Tramp
		   alphat=-1.d0/Tramp
		   alphatt=0.d0
		elseif (iramp.eq.2) then
		   tmp     = t_inter/Tramp 
		   alpha   = tmp * tmp
		   alphat  = -2.0d0 * tmp / Tramp
		   alphatt = 2.0d0 / (Tramp * Tramp)
		elseif (iramp.eq.21) then
		   T2=Tramp**2
		   tim2=t_inter**2
		   timsup=t_inter/Tramp-1.d0
		   if(t_inter.le.(Tramp/2.d0)) then
			  alpha=2.d0*tim2/T2
			  alphat=-4.d0*t_inter/T2
			  alphatt=4.d0/T2
		   else
			  alpha=1.d0-2.d0*timsup**2
			  alphat=4.d0*timsup/Tramp
			  alphatt=-4.d0/T2
		   endif
		elseif (iramp.eq.4) then
		   tmp     = t_inter/Tramp
		   tmp2    = tmp * tmp
		   alpha   = tmp2 * tmp2
		   alphat  = -4.0d0 * tmp2 * tmp / Tramp
		   alphatt = 12.0d0 * tmp2 / (Tramp * Tramp)
		elseif (iramp.eq.41) then
		   T2=Tramp**2
		   tim2=t_inter**2
		   rap=8.d0*t_inter/Tramp
		   rap2=tim2/T2
		   rap22=3.d0*rap2
		   alpha=(rap22-rap+6.d0)*rap2
		   alphat=-(4.d0*rap22-3.d0*rap+12.d0)*t_inter/T2
		   alphatt=(12.d0*rap22-6.d0*rap+12.d0)/T2
		endif
	 else
		alpha=1.d0
		alphat=0.d0
		alphatt=0.d0
	 endif
   endif
END IF
!
return
!
end subroutine ramp
!
! end ramp ************************************************************
!
!
! start dipoles_ini **********************************************************
!
SUBROUTINE dipoles_ini
!
IMPLICIT NONE
!
INTEGER                             :: n_f, iloop
REAL(RP)                            :: x_dip_left, x_dip_right
REAL(RP), ALLOCATABLE, DIMENSION(:) :: z_dip
!
REAL(RP), ALLOCATABLE, DIMENSION(:) :: u
!
n_beta = n1
n_u    = 0
!
ALLOCATE(beta(n_beta,n2), beta_x(n_beta,n2))
ALLOCATE(u(n_u))
!
beta   = 1.0_rp
beta_x = 0.0_rp
n_f = SIZE(wav%ampli)
!
ALLOCATE(z_dip(n_f))
!
WHERE (wav%freq > 1.0_rp)
  z_dip = - 0.35_rp
ELSEWHERE (wav%freq > 0.6_rp)
  z_dip = - 0.2_rp
ELSEWHERE (wav%freq > 0.5_rp)
   z_dip = - 0.4_rp
ELSEWHERE
   z_dip = - 0.6_rp
END WHERE
!
IF (icase == 2) THEN
   z_dip = z_dipmono
END IF
!
DO iloop=1,n_f
   wav%TF(iloop) = TF_adim_new(wav%config%typ_wmk, wav%k(iloop), z_dip(iloop))
ENDDO
!
WRITE(*,*) wav%TF(1)
!
!   mu(freq,x==1 ou z==2)
ALLOCATE(mu(n_f,4))
!
FORALL (iloop=1:n_f)
   mu(iloop,3) = z_dip(iloop)
   mu(iloop,4) = -2.0_rp - z_dip(iloop)
   mu(iloop,2) = wav%ampli(iloop) * wav%TF(iloop) * EXP(i * (wav%k(iloop)*x_dip &
        + wav%phase(iloop)))
   mu(iloop,1) = - i * mu(iloop,2) * TANH(wav%k(iloop) * (z_dip(iloop) +1))
END FORALL
!
WRITE(*,*) ABS(mu(1,2)/TWOPI)
!
! position along x-axis
ALLOCATE(delta_x(n1,n2, -n_side:n_side), delta_x2(n1,n2, -n_side:n_side))
!
!   main dipole
delta_x(:,:,0)  = x - x_dip
delta_x2(:,:,0) = delta_x(:,:,0) * delta_x(:,:,0)
!
!   images
DO iloop=1, n_side
   IF (MOD(iloop, 2) == 0) THEN
      x_dip_right = x_dip + xlen * iloop
      x_dip_left  = x_dip - xlen * iloop
   ELSE
      x_dip_right = - x_dip + xlen * (iloop + 1)
      x_dip_left  = - x_dip - xlen * (iloop - 1)
   END IF
   delta_x(:,:,iloop)   = x - x_dip_right
   delta_x(:,:,-iloop)  = x - x_dip_left
   delta_x2(:,:,iloop)  = delta_x(:,:,iloop)  * delta_x(:,:,iloop)
   delta_x2(:,:,-iloop) = delta_x(:,:,-iloop) * delta_x(:,:,-iloop)
END DO
!
! for z<-h, phase-shift of pi on mu_z
!
! for image in x, phase-shift of pi on mu_x
!
END SUBROUTINE dipoles_ini
!
! end dipoles_ini ****************************************
!
!
! debut phi_dipoles ****************************************
!
!     ===============================================     
      subroutine phi_dipoles(phi_dip, dphi_dip_dt, dphi_dip_dx, dphi_dip_dz, etatd)
!     ===============================================     
!
IMPLICIT NONE
!
! Input variables
REAL(RP), DIMENSION(:,:) :: etatd
!
! Output variables
REAL(RP), DIMENSION(:,:) :: phi_dip, dphi_dip_dt, dphi_dip_dx, dphi_dip_dz
!
! Local variables
INTEGER     :: i1, j1, n_f, i2
REAL(RP)    :: ti,tf
COMPLEX(CP) :: eiwt
REAL(RP), ALLOCATABLE, DIMENSION(:,:)    :: eta_local
REAL(RP), ALLOCATABLE, DIMENSION(:,:)    :: r2, dz_up, dz_down, dz_up2, dz_down2
REAL(RP), ALLOCATABLE, DIMENSION(:,:)    :: dx_left, dx_right, dx_left2, dx_right2
COMPLEX(CP), ALLOCATABLE, DIMENSION(:,:) :: muxor2, muzor2
COMPLEX(CP), ALLOCATABLE, DIMENSION(:,:) :: resu, resu_dx, resu_dz, tmp
!
COMPLEX(CP), ALLOCATABLE, DIMENSION(:,:) :: phi, dphi_dt, dphi_dx, dphi_dz
!            
!	CPU times inlet
if (iCPUtime.eq.1) then
   print*,'entering subroutine phi_dipoles'
   call CPU_TIME(ti)
endif
!
ALLOCATE(eta_local(n_beta,n2))
!
DO i1 = 1, n_beta
   DO i2 = 1, n2
	  eta_local(i1,i2) = etatd(i1,i2)
   END DO
END DO
!
ALLOCATE(r2(n_beta,n2), dz_up(n_beta,n2), dz_down(n_beta,n2), dz_up2(n_beta,n2), dz_down2(n_beta,n2))
ALLOCATE(dx_left(n_beta,n2), dx_right(n_beta,n2), dx_left2(n_beta,n2), dx_right2(n_beta,n2))
ALLOCATE(resu(n_beta,n2), resu_dx(n_beta,n2), resu_dz(n_beta,n2), tmp(n_beta,n2))
ALLOCATE(phi(n_beta,n2), dphi_dt(n_beta,n2), dphi_dx(n_beta,n2), dphi_dz(n_beta,n2))
ALLOCATE(muxor2(n_beta,n2), muzor2(n_beta,n2))
!
phi     = 0.0_cp
dphi_dt = 0.0_cp
dphi_dx = 0.0_cp
dphi_dz = 0.0_cp
!
n_f    = SIZE(wav%ampli)
!
DO i1 = 1, n_f
   dz_up     = eta_local - REAL(mu(i1,3)) ! mu(:,3) is real
   dz_up2    = dz_up * dz_up
   dz_down   = eta_local - REAL(mu(i1,4)) ! mu(:,4) is real
   dz_down2  = dz_down * dz_down
   !
   ! main dipole:
   r2      = 1.0_rp / (delta_x2(1:n_beta,1:n2,0) + dz_up2)
   muxor2  = mu(i1,1) * r2
   muzor2  = mu(i1,2) * r2
   tmp     = muxor2 * delta_x(1:n_beta,1:n2,0) + muzor2 * dz_up
   resu    = tmp
   resu_dx = muxor2 - tmp * 2.0_rp * delta_x(1:n_beta,1:n2,0) * r2
   resu_dz = muzor2 - tmp * 2.0_rp * dz_up * r2
   !
   ! images
   ! main z-image:
   r2      = 1.0_rp / (delta_x2(1:n_beta,1:n2,0) + dz_down2)
   muxor2  = mu(i1,1) * r2
   muzor2  = mu(i1,2) * r2
   tmp     = muxor2 * delta_x(1:n_beta,1:n2,0) - muzor2 * dz_down
   resu    = resu    + tmp
   resu_dx = resu_dx + muxor2 - tmp * 2.0_rp * delta_x(1:n_beta,1:n2,0) * r2
   resu_dz = resu_dz - muzor2 - tmp * 2.0_rp * dz_down * r2
   !
   ! auxiliary images:
   DO j1 = 1, n_side ! a double pair of images
      dx_left   = delta_x(:,:,-j1)
      dx_left2  = delta_x2(:,:,-j1)
      dx_right  = delta_x(:,:,j1)
      dx_right2 = delta_x2(:,:,j1)
      !
      ! left x-image
      r2      = 1.0_rp / (dx_left2 + dz_up2)
      muxor2  = mu(i1,1) * r2
      muzor2  = mu(i1,2) * r2
      IF (MOD(j1,2) == 1) muxor2 = - muxor2
      tmp     = muxor2 * dx_left + muzor2 * dz_up
      resu    = resu    + tmp
      resu_dx = resu_dx + muxor2 - tmp * 2.0_rp * dx_left * r2
      resu_dz = resu_dz + muzor2 - tmp * 2.0_rp * dz_up * r2
      !
      ! right x-image 
      r2      = 1.0_rp / (dx_right2 + dz_up2)
      muxor2  = mu(i1,1) * r2
      muzor2  = mu(i1,2) * r2
      IF (MOD(j1,2) == 1) muxor2 = - muxor2
      tmp     = muxor2 * dx_right + muzor2 * dz_up
      resu    = resu    + tmp
      resu_dx = resu_dx + muxor2 - tmp * 2.0_rp * dx_right * r2
      resu_dz = resu_dz + muzor2 - tmp * 2.0_rp * dz_up * r2
      !
      ! right xz-image
      r2      = 1.0_rp / (dx_right2 + dz_down2)
      muxor2  = mu(i1,1) * r2
      muzor2  = mu(i1,2) * r2
      IF (MOD(j1,2) == 1) muxor2 = - muxor2
      tmp     = muxor2 * dx_right - muzor2 * dz_down
      resu    = resu    + tmp
      resu_dx = resu_dx + muxor2 - tmp * 2.0_rp * dx_right * r2
      resu_dz = resu_dz - muzor2 - tmp * 2.0_rp * dz_down * r2
      !
      ! left xz-image
      r2      = 1.0_rp / (dx_left2 + dz_down2)
      muxor2  = mu(i1,1) * r2
      muzor2  = mu(i1,2) * r2
      IF (MOD(j1,2) == 1) muxor2 = - muxor2
      tmp     = muxor2 * dx_left - muzor2 * dz_down
      resu    = resu    + tmp
      resu_dx = resu_dx + muxor2 - tmp * 2.0_rp * dx_left * r2
      resu_dz = resu_dz - muzor2 - tmp * 2.0_rp * dz_down * r2
   END DO
   eiwt    = EXP(i * TWOPI * wav%freq(i1) * time)
   phi     = phi     + resu    * eiwt
   dphi_dt = dphi_dt + resu    * i * TWOPI * wav%freq(i1) * eiwt
   dphi_dx = dphi_dx + resu_dx * eiwt
   dphi_dz = dphi_dz + resu_dz * eiwt
END DO
!
phi     = phi     / TWOPI
dphi_dt = dphi_dt / TWOPI
dphi_dx = dphi_dx / TWOPI
dphi_dz = dphi_dz / TWOPI
!
phi_dip(1:n1,1:n2)     = 0.0_rp
dphi_dip_dt(1:n1,1:n2) = 0.0_rp
dphi_dip_dx(1:n1,1:n2) = 0.0_rp
dphi_dip_dz(1:n1,1:n2) = 0.0_rp

phi_dip(1:n_beta,1:n2)     = REAL(beta * phi * alpha)
dphi_dip_dt(1:n_beta,1:n2) = REAL(beta * (phi * alphat + dphi_dt * alpha))
dphi_dip_dx(1:n_beta,1:n2) = REAL(beta * dphi_dx * alpha + beta_x * phi * alpha)
dphi_dip_dz(1:n_beta,1:n2) = REAL(beta * dphi_dz * alpha)
!
!	CPU times outlet
if (iCPUtime.eq.1) then
   call CPU_TIME(tf)
   write(*,910)'quitting subroutine phi_dipoles, total CPU time: ',tf-ti,'s'
endif
! 
910 format(a,1ES11.4,a)
!
return
!
END SUBROUTINE phi_dipoles
!
! fin phi_dipoles **********************************************
!
! 
! start match_surface_wmk **************************************************
!
!     =====================================
SUBROUTINE match_surface_wmk(data_wmk,n2,n3,n3_add,l_add,delz)
!     =====================================
!
IMPLICIT NONE
!
! Input variables
REAL(RP), DIMENSION(n3_add,n2), INTENT(INOUT) :: data_wmk
INTEGER, INTENT(IN)                           :: n2,n3,n3_add,l_add
REAL(RP), INTENT(IN)                          :: delz
!
! Local variables
INTEGER :: i2,i3
REAL(RP) :: x_med, f, df, d2f, alpha_add, beta_add, inv_det, xx, xx2, pxx
REAL(RP), DIMENSION(4) :: c_poly_add_2
!
x_med = float(l_add)-1.0d0
do i2=1,n2
  if (l_add.NE.1) then
	 df = (data_wmk(n3-2,i2)-4.0d0*data_wmk(n3-1,i2)&
		   +3.0d0*data_wmk(n3,i2))/(2.0d0*delz)
	 f  = data_wmk(n3,i2)
	 d2f= (-data_wmk(n3-3,i2) + 4.d0*data_wmk(n3-2,i2) - 5.d0*data_wmk(n3-1,i2)&
		  +2.d0*data_wmk(n3,i2))/(delz**2)
	 !
	 alpha_add = 6.d0
	 beta_add  = 1.d0/3.d0
	 !
	 inv_det = 1.d0/(-320*x_med**13*exp(-alpha_add*beta_add)*beta_add**4+128*x_med**13*exp(-alpha_add*beta_add)*beta_add**3 &
		  +288*x_med**13*exp(-alpha_add*beta_add)*beta_add**5-112*x_med**13*exp(-alpha_add*beta_add)*beta_add**6&
		  +16*x_med**13*exp(-alpha_add*beta_add)*beta_add**7)
	 !
	 c_poly_add_2(1)= inv_det*(2*(3*x_med*beta_add**4*df+2*x_med**2*beta_add**4*df*alpha_add+f*alpha_add**2*x_med**2*beta_add**4 &
		  +3*f*alpha_add*x_med*beta_add**4+x_med**2*beta_add**4*d2f+3*f*beta_add**4-12*x_med*beta_add**3*df &
		  -8*x_med**2*beta_add**3*df*alpha_add-4*f*alpha_add**2*x_med**2*beta_add**3-12*f*alpha_add*x_med*beta_add**3 &
		  -4*x_med**2*beta_add**3*d2f-12*f*beta_add**3+4*x_med**2*beta_add**2*d2f+8*x_med**2*beta_add**2*df*alpha_add &
		  +8*f*beta_add**2+8*f*alpha_add*x_med*beta_add**2+4*f*alpha_add**2*x_med**2*beta_add**2+8*x_med*beta_add**2*df &
		  +8*f*beta_add+8*x_med*beta_add*df+8*f*alpha_add*x_med*beta_add+8*f)*(-1+beta_add)*x_med**6*exp(-alpha_add*beta_add))
	 !
	 c_poly_add_2(2) = inv_det*(-2*(24*f-20*f*alpha_add**2*x_med**2*beta_add**3+24*f*alpha_add*x_med*beta_add &
		  +12*x_med**2*beta_add**2*d2f-6*f*alpha_add**2*x_med**2*beta_add**5+3*f*x_med*beta_add**6*alpha_add &
		  +f*x_med**2*beta_add**6*alpha_add**2-60*f*alpha_add*x_med*beta_add**3+12*f*alpha_add**2*x_med**2*beta_add**2 &
		  +15*f*alpha_add**2*x_med**2*beta_add**4+45*f*alpha_add*x_med*beta_add**4-18*f*alpha_add*x_med*beta_add**5 &
		  +2*x_med**2*beta_add**6*df*alpha_add-12*x_med**2*beta_add**5*df*alpha_add+30*x_med**2*beta_add**4*df*alpha_add &
		  -40*x_med**2*beta_add**3*df*alpha_add+24*x_med**2*beta_add**2*df*alpha_add+24*f*beta_add**2-60*f*beta_add**3 &
		  +24*f*beta_add+45*f*beta_add**4-18*f*beta_add**5+3*f*beta_add**6+24*f*alpha_add*x_med*beta_add**2+3*x_med*beta_add**6*df &
		  -18*x_med*beta_add**5*df+45*x_med*beta_add**4*df-60*x_med*beta_add**3*df+24*x_med*beta_add**2*df &
		  +x_med**2*beta_add**6*d2f-6*x_med**2*beta_add**5*d2f+15*x_med**2*beta_add**4*d2f-20*x_med**2*beta_add**3*d2f &
		  +24*x_med*beta_add*df)*(-1+beta_add)*x_med**8*exp(-alpha_add*beta_add))
	 !
	 c_poly_add_2(3) = inv_det*(2*(24*f-28*f*alpha_add**2*x_med**2*beta_add**3+24*f*alpha_add*x_med*beta_add &
		  +12*x_med**2*beta_add**2*d2f-12*f*alpha_add**2*x_med**2*beta_add**5+10*f*x_med*beta_add**6*alpha_add &
		  +2*f*x_med**2*beta_add**6*alpha_add**2-116*f*alpha_add*x_med*beta_add**3+12*f*alpha_add**2*x_med**2*beta_add**2 &
		  +27*f*alpha_add**2*x_med**2*beta_add**4+129*f*alpha_add*x_med*beta_add**4-60*f*alpha_add*x_med*beta_add**5 &
		  +4*x_med**2*beta_add**6*df*alpha_add-24*x_med**2*beta_add**5*df*alpha_add &
		  +54*x_med**2*beta_add**4*df*alpha_add-56*x_med**2*beta_add**3*df*alpha_add+24*x_med**2*beta_add**2*df*alpha_add &
		  +24*f*beta_add**2-116*f*beta_add**3+24*f*beta_add+129*f*beta_add**4-60*f*beta_add**5+10*f*beta_add**6 &
		  +24*f*alpha_add*x_med*beta_add**2+10*x_med*beta_add**6*df-60*x_med*beta_add**5*df+129*x_med*beta_add**4*df &
		  -116*x_med*beta_add**3*df+24*x_med*beta_add**2*df+2*x_med**2*beta_add**6*d2f-12*x_med**2*beta_add**5*d2f &
		  +27*x_med**2*beta_add**4*d2f-28*x_med**2*beta_add**3*d2f+24*x_med*beta_add*df)*(-1+beta_add) &
		  *x_med**10*exp(-alpha_add*beta_add))
	 !
	 c_poly_add_2(4) = inv_det*(-2*(8*f-12*f*alpha_add**2*x_med**2*beta_add**3+8*f*alpha_add*x_med*beta_add &
		  +4*x_med**2*beta_add**2*d2f-6*f*alpha_add**2*x_med**2*beta_add**5+7*f*x_med*beta_add**6*alpha_add &
		  +f*x_med**2*beta_add**6*alpha_add**2-68*f*alpha_add*x_med*beta_add**3+4*f*alpha_add**2*x_med**2*beta_add**2 &
		  +13*f*alpha_add**2*x_med**2*beta_add**4+87*f*alpha_add*x_med*beta_add**4-42*f*alpha_add*x_med*beta_add**5 &
		  +2*x_med**2*beta_add**6*df*alpha_add-12*x_med**2*beta_add**5*df*alpha_add+26*x_med**2*beta_add**4*df*alpha_add &
		  -24*x_med**2*beta_add**3*df*alpha_add+8*x_med**2*beta_add**2*df*alpha_add+8*f*beta_add**2-132*f*beta_add**3 &
		  +8*f*beta_add+183*f*beta_add**4-90*f*beta_add**5+15*f*beta_add**6+8*f*alpha_add*x_med*beta_add**2+7*x_med*beta_add**6*df &
		  -42*x_med*beta_add**5*df+87*x_med*beta_add**4*df-68*x_med*beta_add**3*df+8*x_med*beta_add**2*df &
		  +x_med**2*beta_add**6*d2f-6*x_med**2*beta_add**5*d2f+13*x_med**2*beta_add**4*d2f-12*x_med**2*beta_add**3*d2f &
		  +8*x_med*beta_add*df)*(-1+beta_add)*x_med**12*exp(-alpha_add*beta_add))
  end if
  !
  do i3=n3+1,n3_add
	 xx=dfloat(i3-1)/dfloat(n3-1)-dfloat(l_add)
	 xx2=xx*xx
	 !
	 !  seventh order polynom
	 !
	 pxx=(((c_poly_add_2(1)*xx2+c_poly_add_2(2))*xx2+c_poly_add_2(3))*xx2+c_poly_add_2(4))*xx*exp(-alpha_add*(xx+x_med))
	 data_wmk(i3,i2) = pxx
  end do
end do
!
RETURN
!
END SUBROUTINE match_surface_wmk
!
! end match_surface_wmk ****************************************************
!
!
!end module
!
END MODULE resol_wmkr
