MODULE runge_kutta
!
! This module defines the time-stepping procedure based on Runge-Kutta for HOS-NWT
! This is with adaptative time-stepping (Runge-Kutta Fehlberg scheme) 
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
USE resol_HOS
USE resol_wmkr
USE Fourier_FFTW
USE velocities
USE variablechange
!
CONTAINS
!
! start fill_butcher_array ******************************************************
!
!     ==============================
      SUBROUTINE fill_butcher_array()
!     ==============================
!
IMPLICIT NONE
!
WRITE(*,*) '4(5) Runge Kutta Fehlberg scheme'
RK_A      = 0.0_rp
RK_A(2,1) = 1.0_rp / 2.0_rp
RK_A(3,2) = 1.0_rp / 2.0_rp
RK_A(4,3) = 1.0_rp
RK_A(5,1) =-1.0_rp
RK_A(5,2) = 2.0_rp
RK_b(1)   = 1.0_rp / 6.0_rp
RK_b(2)   = 1.0_rp / 3.0_rp
RK_b(3)   = RK_b(2)
RK_b(4)   = RK_b(1)
RK_c(1)   = 0.0_rp
RK_c(2)   = 1.0_rp / 2.0_rp
RK_c(3)   = 1.0_rp / 2.0_rp
RK_c(4)   = 1.0_rp
RK_c(5)   = 1.0_rp
RK_e      = 0.0_rp
RK_e(1)   = 1.0_rp / 6.0_rp
RK_e(2)   = 2.0_rp / 3.0_rp
RK_e(5)   = 1.0_rp / 6.0_rp
RK_p      = 4
RK_substep(1) = 0
RK_substep(2) = 0
RK_substep(3) = 0
RK_substep(4) = 0
RK_substep(5) = 0
!
END SUBROUTINE fill_butcher_array
!
! end fill_butcher_array ******************************************************
!
!
!
! debut runge4_adapt2
!
!     ===================================================
subroutine runge4_adapt2(t_o, etark, phisrk, a_add_rk_1, a_add_rk_2 ,a_add_rkramp, a_add_rk_3, eta1st_l, eta2nd_l, &
     a1st_l, a2nd_l, h_loc, err, tiny)
!     ==================================================
!
IMPLICIT NONE
!
! Input variables
REAL(RP), INTENT(IN)                          :: h_loc, t_o
REAL(RP), DIMENSION(m1,m2), INTENT(INOUT)     :: etark, phisrk, eta1st_l, eta2nd_l,a1st_l, a2nd_l
REAL(RP), DIMENSION(m3_add,m2), INTENT(INOUT) :: a_add_rk_1,a_add_rk_2,a_add_rkramp,a_add_rk_3
!
REAL(RP),INTENT(IN)           :: tiny
REAL(RP),INTENT(OUT)          :: err
!
!% LOCAL VARIABLES
INTEGER                        :: i1,i2,jloop,ns
!
! Free surface variables with variable change
REAL(RP), DIMENSION(m1,m2)       :: G1m, G2m, dG1m, dG2m , G1, G2
REAL(RP), DIMENSION(m1,m2,2)     :: erm
REAL(RP), DIMENSION(m1,m2,RK_s)  :: k_G1, k_G2
!
! Wavemaking quantities
REAL(RP), DIMENSION(m3_add,m2)      :: da_add_m_1, a_add_m_1, da_add_m_2, a_add_m_2, da_add_mramp
REAL(RP), DIMENSION(m3_add,m2)      :: a_add_mramp, da_add_m_3, a_add_m_3
REAL(RP), DIMENSION(m3_add,m2,3)    :: erm_add
REAL(RP), DIMENSION(m3_add,m2,RK_s) :: k_a_add_1, k_a_add_2, k_a_add_ramp, k_a_add_3
!
! 1st and 2nd order free surface quantities
REAL(RP), DIMENSION(m1,m2)      :: a1stm, da1stm, eta1stm, deta1stm
REAL(RP), DIMENSION(m1,m2,RK_s) :: k_a1st, k_eta1st
REAL(RP), DIMENSION(m1,m2)      :: a2ndm, da2ndm, eta2ndm, deta2ndm
REAL(RP), DIMENSION(m1,m2,RK_s) :: k_a2nd, k_eta2nd
!
REAL(RP) :: ti, tf, phisrk0, error
!
REAL(RP), DIMENSION(m1,m2) :: dphisrk
!
!	CPU times inlet
if(iCPUtime.eq.1) then
   print*,'entering subroutine runge4_adapt2'
   call CPU_TIME(ti)
endif
!
! Construction of new variable G
do i1=1,n1
   do i2=1,n2
      G1(i1,i2) = etark(i1,i2)
      G2(i1,i2) = phisrk(i1,i2)
   end do
end do
!
G1 = space_2_Fourier(G1,'cos','cos')
G2 = space_2_Fourier(G2,'cos','cos')
!
phisrk0 = G2(1,1)
do i1=1,n1
   do i2=1,n2
      G2(i1,i2) = G2(i1,i2)*omega(i1,i2)
   end do
end do
!
time = t_o
!
! Derivatives at t=time
call solveHOS2(a1st_l, k_a1st(:,:,1), eta1st_l,k_eta1st(:,:,1), a2nd_l,k_a2nd(:,:,1), eta2nd_l,k_eta2nd(:,:,1),&
     G2,k_G2(:,:,1),G1,k_G1(:,:,1),a_add_rk_1,a_add_rk_2,k_a_add_1(:,:,1),k_a_add_2(:,:,1), &
     k_a_add_ramp(:,:,1),a_add_rk_3,k_a_add_3(:,:,1), RK_substep(1),h_loc*RK_c(1))
!
k_G2(1,1,1)=0.0_rp
!
! Runge-Kutta loop
DO jloop = 2, RK_s
   dG1m(1:n1,1:n2)             = 0.0_rp
   dG2m(1:n1,1:n2)             = 0.0_rp
   da_add_m_1(1:n3_add,1:n2)   = 0.0_rp
   da_add_m_2(1:n3_add,1:n2)   = 0.0_rp
   da_add_mramp(1:n3_add,1:n2) = 0.0_rp
   da_add_m_3(1:n3_add,1:n2)   = 0.0_rp
   !
   da1stm(1:n1,1:n2)   = 0.0_rp
   deta1stm(1:n1,1:n2) = 0.0_rp
   da2ndm(1:n1,1:n2)   = 0.0_rp
   deta2ndm(1:n1,1:n2) = 0.0_rp
   !
   DO ns = 1, jloop-1
      dG1m(1:n1,1:n2)              = dG1m(1:n1,1:n2) + k_G1(1:n1,1:n2,ns) * RK_A(jloop, ns)
      dG2m(1:n1,1:n2)              = dG2m(1:n1,1:n2) + k_G2(1:n1,1:n2,ns) * RK_A(jloop, ns)
      da_add_m_1(1:n3_add,1:n2)    = da_add_m_1(1:n3_add,1:n2) + k_a_add_1(1:n3_add,1:n2,ns) * RK_A(jloop, ns)
      da_add_m_2(1:n3_add,1:n2)    = da_add_m_2(1:n3_add,1:n2) + k_a_add_2(1:n3_add,1:n2,ns) * RK_A(jloop, ns)
      da_add_mramp(1:n3_add,1:n2)  = da_add_mramp(1:n3_add,1:n2) + k_a_add_ramp(1:n3_add,1:n2,ns) * RK_A(jloop, ns)
      da_add_m_3(1:n3_add,1:n2)    = da_add_m_3(1:n3_add,1:n2) + k_a_add_3(1:n3_add,1:n2,ns) * RK_A(jloop, ns)
      !
      da1stm(1:n1,1:n2)   = da1stm(1:n1,1:n2) + k_a1st(1:n1,1:n2,ns) * RK_A(jloop, ns)
      deta1stm(1:n1,1:n2) = deta1stm(1:n1,1:n2) + k_eta1st(1:n1,1:n2,ns) * RK_A(jloop, ns)
      da2ndm(1:n1,1:n2)   = da2ndm(1:n1,1:n2) + k_a2nd(1:n1,1:n2,ns) * RK_A(jloop, ns)
      deta2ndm(1:n1,1:n2) = deta2ndm(1:n1,1:n2) + k_eta2nd(1:n1,1:n2,ns) * RK_A(jloop, ns)
   END DO
   !
   G1m(1:n1,1:n2)             = G1(1:n1,1:n2) + h_loc * dG1m(1:n1,1:n2)
   G2m(1:n1,1:n2)             = G2(1:n1,1:n2) + h_loc * dG2m(1:n1,1:n2)
   a_add_m_1(1:n3_add,1:n2)   = a_add_rk_1(1:n3_add,1:n2) + h_loc * da_add_m_1(1:n3_add,1:n2)
   a_add_m_2(1:n3_add,1:n2)   = a_add_rk_2(1:n3_add,1:n2) + h_loc * da_add_m_2(1:n3_add,1:n2)
   a_add_mramp(1:n3_add,1:n2) = a_add_rkramp(1:n3_add,1:n2) + h_loc * da_add_mramp(1:n3_add,1:n2)
   a_add_m_3(1:n3_add,1:n2)   = a_add_rk_3(1:n3_add,1:n2) + h_loc * da_add_m_3(1:n3_add,1:n2)
   !
   a1stm(1:n1,1:n2)   = a1st_l(1:n1,1:n2)   + h_loc*da1stm(1:n1,1:n2)
   eta1stm(1:n1,1:n2) = eta1st_l(1:n1,1:n2) + h_loc*deta1stm(1:n1,1:n2)
   a2ndm(1:n1,1:n2)   = a2nd_l(1:n1,1:n2)   + h_loc*da2ndm(1:n1,1:n2)
   eta2ndm(1:n1,1:n2) = eta2nd_l(1:n1,1:n2) + h_loc*deta2ndm(1:n1,1:n2)
   !
   time = t_o + h_loc * RK_c(jloop)
   !
   IF (jloop == RK_s) THEN
      call solveHOS2(a1stm, k_a1st(:,:,jloop),eta1stm, k_eta1st(:,:,jloop),a2ndm, k_a2nd(:,:,jloop), &
           eta2ndm, k_eta2nd(:,:,jloop), G2m, k_G2(:,:,jloop), G1m, k_G1(:,:,jloop), &
           a_add_m_1,a_add_m_2,k_a_add_1(:,:,jloop),k_a_add_2(:,:,jloop),k_a_add_ramp(:,:,jloop),&
           a_add_m_3,k_a_add_3(:,:,jloop) , RK_substep(jloop), h_loc * RK_c(jloop))
      k_G2(1,1,jloop)=0.0_rp
   ELSE
      call solveHOS2(a1stm, k_a1st(:,:,jloop), eta1stm, k_eta1st(:,:,jloop),a2ndm, k_a2nd(:,:,jloop),&
           eta2ndm, k_eta2nd(:,:,jloop),G2m, k_G2(:,:,jloop), G1m, k_G1(:,:,jloop),&
           a_add_m_1,a_add_m_2,k_a_add_1(:,:,jloop),k_a_add_2(:,:,jloop),k_a_add_ramp(:,:,jloop),&
           a_add_m_3,k_a_add_3(:,:,jloop),  RK_substep(jloop), h_loc * RK_c(jloop))
      k_G2(1,1,jloop)=0.0_rp
   END IF
END DO
!
! Initialization to 0
erm(1:n1,1:n2,:)            = 0.0_rp
erm_add(1:n3_add,1:n2,:)    = 0.0_rp
dG1m(1:n1,1:n2)             = 0.0_rp
dG2m(1:n1,1:n2)             = 0.0_rp
da_add_m_1(1:n3_add,1:n2)   = 0.0_rp
da_add_m_2(1:n3_add,1:n2)   = 0.0_rp
da_add_mramp(1:n3_add,1:n2) = 0.0_rp
da_add_m_3(1:n3_add,1:n2)   = 0.0_rp
!
da1stm(1:n1,1:n2)   = 0.0_rp
deta1stm(1:n1,1:n2) = 0.0_rp
da2ndm(1:n1,1:n2)   = 0.0_rp
deta2ndm(1:n1,1:n2) = 0.0_rp
!
! Construct next time-step
DO jloop = 1, RK_s
   dG1m(1:n1,1:n2)             = dG1m(1:n1,1:n2) + k_G1(1:n1,1:n2,jloop) * RK_b(jloop)
   dG2m(1:n1,1:n2)             = dG2m(1:n1,1:n2) + k_G2(1:n1,1:n2,jloop) * RK_b(jloop)
   da_add_m_1(1:n3_add,1:n2)   = da_add_m_1(1:n3_add,1:n2) + k_a_add_1(1:n3_add,1:n2,jloop) * RK_b(jloop)
   da_add_m_2(1:n3_add,1:n2)   = da_add_m_2(1:n3_add,1:n2) + k_a_add_2(1:n3_add,1:n2,jloop) * RK_b(jloop)
   da_add_m_3(1:n3_add,1:n2)   = da_add_m_3(1:n3_add,1:n2) + k_a_add_3(1:n3_add,1:n2,jloop) * RK_b(jloop)
   da_add_mramp(1:n3_add,1:n2) = da_add_mramp(1:n3_add,1:n2) + k_a_add_ramp(1:n3_add,1:n2,jloop) * RK_b(jloop)
   !
   da1stm(1:n1,1:n2)   = da1stm(1:n1,1:n2) + k_a1st(1:n1,1:n2,jloop) * RK_b(jloop)
   deta1stm(1:n1,1:n2) = deta1stm(1:n1,1:n2) + k_eta1st(1:n1,1:n2,jloop) * RK_b(jloop)
   da2ndm(1:n1,1:n2)   = da2ndm(1:n1,1:n2) + k_a2nd(1:n1,1:n2,jloop) * RK_b(jloop)
   deta2ndm(1:n1,1:n2) = deta2ndm(1:n1,1:n2) + k_eta2nd(1:n1,1:n2,jloop) * RK_b(jloop)
   !
   erm(1:n1,1:n2,1)         = erm(1:n1,1:n2,1) + k_G1(1:n1,1:n2,jloop) * RK_e(jloop)
   erm(1:n1,1:n2,2)         = erm(1:n1,1:n2,2) + k_G2(1:n1,1:n2,jloop) * RK_e(jloop)
   erm_add(1:n3_add,1:n2,1) = erm_add(1:n3_add,1:n2,1) + k_a_add_1(1:n3_add,1:n2,jloop) * RK_e(jloop)
   erm_add(1:n3_add,1:n2,2) = erm_add(1:n3_add,1:n2,2) + k_a_add_2(1:n3_add,1:n2,jloop) * RK_e(jloop)
   erm_add(1:n3_add,1:n2,3) = erm_add(1:n3_add,1:n2,3) + k_a_add_3(1:n3_add,1:n2,jloop) * RK_e(jloop)
END DO
!
! Time-stepping
G1(1:n1,1:n2)               = G1(1:n1,1:n2) + h_loc * dG1m(1:n1,1:n2)
G2(1:n1,1:n2)               = G2(1:n1,1:n2) + h_loc * dG2m(1:n1,1:n2)
a_add_rk_1(1:n3_add,1:n2)   = a_add_rk_1(1:n3_add,1:n2) + h_loc*da_add_m_1(1:n3_add,1:n2)
a_add_rk_2(1:n3_add,1:n2)   = a_add_rk_2(1:n3_add,1:n2) + h_loc*da_add_m_2(1:n3_add,1:n2)
a_add_rk_3(1:n3_add,1:n2)   = a_add_rk_3(1:n3_add,1:n2) + h_loc*da_add_m_3(1:n3_add,1:n2)
a_add_rkramp(1:n3_add,1:n2) = a_add_rkramp(1:n3_add,1:n2) + h_loc*da_add_mramp(1:n3_add,1:n2)
!
a1st_l(1:n1,1:n2)   = a1st_l(1:n1,1:n2)   + h_loc*da1stm(1:n1,1:n2)
eta1st_l(1:n1,1:n2) = eta1st_l(1:n1,1:n2) + h_loc*deta1stm(1:n1,1:n2)
a2nd_l(1:n1,1:n2)   = a2nd_l(1:n1,1:n2)   + h_loc*da2ndm(1:n1,1:n2)
eta2nd_l(1:n1,1:n2) = eta2nd_l(1:n1,1:n2) + h_loc*deta2ndm(1:n1,1:n2)
!
da1st(1:n1,1:n2)         = k_a1st(1:n1,1:n2,RK_s)
da2nd(1:n1,1:n2)         = k_a2nd(1:n1,1:n2,RK_s)
da1st_add(1:n3_add,1:n2) = k_a_add_1(1:n3_add,1:n2,RK_s)
da2nd_add(1:n3_add,1:n2) = k_a_add_2(1:n3_add,1:n2,RK_s)
da3rd_add(1:n3_add,1:n2) = k_a_add_3(1:n3_add,1:n2,RK_s)
deta1st(1:n1,1:n2)       = k_eta1st(1:n1,1:n2,RK_s)
deta2nd(1:n1,1:n2)       = k_eta2nd(1:n1,1:n2,RK_s)
!
! Error evaluation for adaptative time-stepping
erm(1:n1,1:n2,1)         = (erm(1:n1,1:n2,1) - dG1m(1:n1,1:n2)) * h_loc
erm(1:n1,1:n2,2)         = (erm(1:n1,1:n2,2) - dG2m(1:n1,1:n2)) * h_loc
erm_add(1:n3_add,1:n2,1) = (erm_add(1:n3_add,1:n2,1) - da_add_m_1(1:n3_add,1:n2)) * h_loc
erm_add(1:n3_add,1:n2,2) = (erm_add(1:n3_add,1:n2,2) - da_add_m_2(1:n3_add,1:n2)) * h_loc
erm_add(1:n3_add,1:n2,3) = (erm_add(1:n3_add,1:n2,3) - da_add_m_3(1:n3_add,1:n2)) * h_loc
erm(1:n1,1:n2,1)         = ABS(erm(1:n1,1:n2,1)) / (MAXVAL(ABS(G1(1:n1,1:n2))) + tiny)
erm(1:n1,1:n2,2)         = ABS(erm(1:n1,1:n2,2)) / (MAXVAL(ABS(G2(1:n1,1:n2))) + tiny)
erm_add(1:n3_add,1:n2,1) = ABS(erm_add(1:n3_add,1:n2,1)) / (MAXVAL(ABS(a_add_rk_1(1:n3_add,1:n2))) + tiny)
erm_add(1:n3_add,1:n2,2) = ABS(erm_add(1:n3_add,1:n2,2)) / (MAXVAL(ABS(a_add_rk_2(1:n3_add,1:n2))) + tiny)
erm_add(1:n3_add,1:n2,3) = ABS(erm_add(1:n3_add,1:n2,3)) / (MAXVAL(ABS(a_add_rk_3(1:n3_add,1:n2))) + tiny)
!
err = MAXVAL(erm(1:n1,1:n2,1:2))
err = MAX(err,MAXVAL(erm_add(1:n3_add,1:n2,1:3)))
!
call GtoF(G1,G2,etark,phisrk,h_loc,iCPUtime)
! For time-derivative of eta
call GtoF(k_G1(:,:,RK_s),k_G2(:,:,RK_s),deta,dphisrk,h_loc,iCPUtime)
!
phisrk(1,1) = phisrk0
!
do i1=2,n1
   phisrk(i1,1) = phisrk(i1,1)/omega(i1,1)
   deta(i1,1)   = deta(i1,1) + phisrk(i1,1)*omega(i1,1)**2
end do
do i1=1,n1
   do i2=2,n2
      phisrk(i1,i2) = phisrk(i1,i2)/omega(i1,i2)
      deta(i1,i2)   = deta(i1,i2) + phisrk(i1,i2)*omega(i1,i2)**2
   end do
end do
!
! Correction of mean volume
!
! modesFS is the Free Surface in space
modesFS  = Fourier_2_space(etark,'cos','cos')
!
! Correct the volume, i.e. first mode on eta: set the value to theoretical value
! This is useful for very long time simulations
! FIXME: check the influence on Results
if(igeom.EQ.2) then
    etark(1,1) = 0.5_rp*alpha*(1-d_hinge)*pos1st(n3,1) / xlen + modesFS(1,1)*alpha*pos1st(n3,1) / xlen
elseif(igeom.EQ.1) then
    etark(1,1) = alpha*pos1st(n3,1) / xlen + modesFS(1,1)*alpha*pos1st(n3,1) / xlen
else
    etark(1,1) = 0.0_rp
endif
!
! modesFS are the modes of 'corrected' FS
do i1=1,n1
   do i2=1,n2
      modesFS(i1,i2) = etark(i1,i2)
   end do
enddo
!
etark  = Fourier_2_space(etark,'cos','cos')
phisrk = Fourier_2_space(phisrk,'cos','cos')
deta   = Fourier_2_space(deta, 'cos','cos')
!
! Output of volumic information of the wavefield computed with HOS-NWT
! This may be used for coupling, e.g. SWENSE method or with velocity/pressure cards 
! (Available in post processing)
IF(i_sw.EQ.1) THEN
   IF(ABS(time_cur + h_loc-time_next).LT.tiny) then
		print*, 'H2 operator'
		CALL HOSvel2_SWEET(10,etark,phisrk,a_add_rk_1,da1st_add,a_add_rk_2,da2nd_add&
           ,a_add_rk_3,da3rd_add,a1st_l,da1st, eta1st_l, deta1st, a2nd_l, da2nd, eta2nd_l, deta2nd, &
           G2, k_G2(:,:,RK_s),G1, k_G1(:,:,RK_s), k_a_add_ramp(:,:,RK_s),h_loc)

    	CALL error_FS(modesspecx,modesspecy,modesspecz,modesspect,modesadd,modesaddt, modesFS, modesFSt, &
            vitxref_FS,vityref_FS, vitzref_FS, phitref_FS, error)

	IF (error.GT.0.1_rp) THEN
		print*, 'Max relative error = ', error, '... direct method'    
		  CALL HOSvel2_direct(etark,phisrk,a_add_rk_1,da1st_add,a_add_rk_2,da2nd_add&
			   ,a_add_rk_3,da3rd_add,a1st_l,da1st, eta1st_l, deta1st, a2nd_l, da2nd, eta2nd_l, deta2nd, &
			   G2, k_G2(:,:,RK_s),G1, k_G1(:,:,RK_s), k_a_add_ramp(:,:,RK_s),h_loc)
			   
		  CALL error_FS(modesspecx,modesspecy,modesspecz,modesspect,modesadd,modesaddt, modesFS, modesFSt, &
            vitxref_FS,vityref_FS, vitzref_FS, phitref_FS, error)
	ENDIF        
   ENDIF
ENDIF  
!
!	CPU times outlet
if(iCPUtime.eq.1) then
   call CPU_TIME(tf)
   write(*,910)'quitting subroutine runge4bis, total CPU time: ',tf-ti,'s'
endif
!   
910 format(a,1ES11.4,a)
!
return
!      
end subroutine runge4_adapt2
! 
! end Runge-Kutta 4 adapt2 ***************************************************
!
END MODULE runge_kutta
