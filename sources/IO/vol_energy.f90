MODULE vol_energy
!
! This module contains routines for volume and energy evaluation 
! during an HOS-NWT computation
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
USE variables
!
!
!
CONTAINS
!
! start volume ********************************************************
!
!     =======================
SUBROUTINE volume
!     =======================
!
!     volume conservation
!     >>>>>>>>>>>>>>>>>>>
!
!     on the free surface
!
IMPLICIT NONE
!
! LOCAL VARIABLES
!
INTEGER  :: iloop
REAL(RP), DIMENSION(m1) :: sumvolx
!
do iloop=1,n1
   call int2(eta,sumvolx,iloop) ! warning: without dely
end do
!
call int1(sumvolx,volfs,n1,delx) 
!
volfs = volfs / xlen
if (n2.ne.1) volfs = volfs / ylen * dely
!
voltot = volfs
!
END SUBROUTINE volume
!
! end volume **********************************************************
!
! start energy ********************************************************
!
!     =======================
SUBROUTINE energy(eta, phis, deta)
!     =======================
!	        
IMPLICIT NONE
!
! INPUT VARIABLES
REAL(RP), DIMENSION(m1,m2), INTENT(IN) :: eta, phis, deta
! LOCAL VARIABLES
INTEGER                    :: iloop
REAL(RP), DIMENSION(m1)    :: sumenex
REAL(RP), DIMENSION(m1,m2) :: integxy
!          
!     energy conservation
!     >>>>>>>>>>>>>>>>>>>
!
!	  potential energy
!
!	  on the free surface
!
integxy = eta * eta
do iloop=1,n1
   call int2(integxy,sumenex,iloop) ! warning: without dely
end do
!
call int1(sumenex,enepot,n1,delx)
!
enepot = enepot / 2.0_rp
!
if (n2.ne.1) enepot = enepot * dely
!
!	  kinetic energy
!
!	  on the free surface
!
integxy = phis * deta
!
do iloop=1,n1
   call int2(integxy,sumenex,iloop) ! warning: without dely
end do
!
call int1(sumenex,enekin,n1,delx)      
!
enekin = enekin / 2.0_rp
!
if (n2.ne.1) enekin = enekin * dely
!
!	  total energy
!
enetot = enepot + enekin
!
END SUBROUTINE energy
!
! end energy **********************************************************
!
!     ===============================================     
SUBROUTINE phi_dip_chk_flux(n_z_chk, dVx0, dVxLx, dVz)
!     ===============================================     
!
IMPLICIT NONE
!
! Input variables
INTEGER     :: n_z_chk
! Output variables
REAL(RP)    :: dVx0, dVxLx, dVz
! Local variables
INTEGER, PARAMETER :: m_z_chk = 10
INTEGER     :: i1, j1, n_f, iloop
REAL(RP)    :: ti,tf
COMPLEX(CP) :: eiwt
REAL(RP), DIMENSION(2) :: x_base
REAL(RP), DIMENSION(n_z_chk+1)      :: r2, dz_up, dz_down, dz_up2, dz_down2, dz
REAL(RP), DIMENSION(n_z_chk+1)      :: dx_left, dx_right, dx_left2, dx_right2
COMPLEX(CP), DIMENSION(n_z_chk+1)   :: muxor2, muzor2
COMPLEX(CP), DIMENSION(n_z_chk+1,2) :: dphi_dx, resu_dx
COMPLEX(CP), DIMENSION(n_z_chk+1)   :: tmp
!
REAL(RP), DIMENSION(n1,n2)    :: r2n, dz_upn, dz_downn, dz_up2n, dz_down2n
REAL(RP), DIMENSION(n1,n2)    :: dx_leftn, dx_rightn, dx_left2n, dx_right2n
COMPLEX(CP), DIMENSION(n1,n2) :: muxor2n, muzor2n
COMPLEX(CP), DIMENSION(n1,n2) :: tmpn
COMPLEX(CP), DIMENSION(n1,n2) :: dphi_dzn, resu_dzn
!            
!	CPU times inlet
!
if (iCPUtime.eq.1) then
   print*,'entering subroutine phi_dip_chk_flux'
   call CPU_TIME(ti)
endif
!
IF (n_z_chk > m_z_chk) THEN
   WRITE(*,'(A)') 'm_z_chk too small in phi_dip_chk_flux'
   STOP 
END IF
!
x_base = (/-x_dip, x_dip/)
!
dz = - (/(iloop,iloop=0,n_z_chk)/) / n_z_chk
dphi_dx = 0.0_cp
n_f    = SIZE(wav%ampli)
!
DO i1 = 1, n_f
   dz_up     = dz - REAL(mu(i1,3)) ! mu(:,3) is real
   dz_up2    = dz_up * dz_up
   dz_down   = dz - REAL(mu(i1,4)) ! mu(:,4) is real
   dz_down2  = dz_down * dz_down
   !
   ! main dipole on x=0 wall
   !
   dx_right  = 0.0_rp - x_dip
   dx_right2 = dx_right * dx_right
   r2      = 1.0_rp / (dx_right2 + dz_up2)
   muxor2  = mu(i1,1) * r2
   muzor2  = mu(i1,2) * r2
   tmp     = muxor2 * dx_right + muzor2 * dz_up
   resu_dx(:,1) = muxor2 - tmp * 2.0_rp * dx_right * r2
   ! main dipole z-image
   r2      = 1.0_rp / (dx_right2 + dz_down2)
   muxor2  = mu(i1,1) * r2
   muzor2  = mu(i1,2) * r2
   tmp     = muxor2 * dx_right - muzor2 * dz_down
   resu_dx(:,1) = resu_dx(:,1) + muxor2 - tmp * 2.0_rp * dx_right * r2
   !
   ! main dipole on x=Lx wall
   !
   dx_right  = xlen - x_dip
   dx_right2 = dx_right * dx_right
   r2      = 1.0_rp / (dx_right2 + dz_up2)
   muxor2  = mu(i1,1) * r2
   muzor2  = mu(i1,2) * r2
   tmp     = muxor2 * dx_right + muzor2 * dz_up
   resu_dx(:,2) = muxor2 - tmp * 2.0_rp * dx_right * r2
   ! main dipole z-image
   r2      = 1.0_rp / (dx_right2 + dz_down2)
   muxor2  = mu(i1,1) * r2
   muzor2  = mu(i1,2) * r2
   tmp     = muxor2 * dx_right - muzor2 * dz_down
   resu_dx(:,2) = resu_dx(:,2) + muxor2 - tmp * 2.0_rp * dx_right * r2
   !
   DO j1 = 1, n_side ! a double pair of images
      IF (MOD(j1,2) == 0) THEN
         dx_left   = 0.0_rp - (x_dip - xlen * j1)
         dx_right  = 0.0_rp - (x_dip + xlen * j1)
      ELSE
         dx_left   = 0.0_rp - (- x_dip - xlen * (j1-1))
         dx_right  = 0.0_rp - (- x_dip + xlen * (j1+1))
      END IF
      dx_left2  = dx_left  * dx_left
      dx_right2 = dx_right * dx_right
      ! right x-image on x = 0 wall
      r2      = 1.0_rp / (dx_right2 + dz_up2)
      muxor2  = mu(i1,1) * r2
      muzor2  = mu(i1,2) * r2
      IF (MOD(j1,2) == 1) muxor2 = - muxor2
      tmp     = muxor2 * dx_right + muzor2 * dz_up
      resu_dx(:,1) = resu_dx(:,1) + muxor2 - tmp * 2.0_rp * dx_right * r2
      ! left x-image on x = 0 wall
      r2      = 1.0_rp / (dx_left2 + dz_up2)
      muxor2  = mu(i1,1) * r2
      muzor2  = mu(i1,2) * r2
      IF (MOD(j1,2) == 1) muxor2 = - muxor2
      tmp     = muxor2 * dx_left + muzor2 * dz_up
      resu_dx(:,1) = resu_dx(:,1) + muxor2 - tmp * 2.0_rp * dx_left * r2
      ! right xz-image on x = 0 wall
      r2      = 1.0_rp / (dx_right2 + dz_down2)
      muxor2  = mu(i1,1) * r2
      muzor2  = mu(i1,2) * r2
      IF (MOD(j1,2) == 1) muxor2 = - muxor2
      tmp     = muxor2 * dx_right - muzor2 * dz_down
      resu_dx(:,1) = resu_dx(:,1) + muxor2 - tmp * 2.0_rp * dx_right * r2
      ! left xz-image on x = 0 wall
      r2      = 1.0_rp / (dx_left2 + dz_down2)
      muxor2  = mu(i1,1) * r2
      muzor2  = mu(i1,2) * r2
      IF (MOD(j1,2) == 1) muxor2 = - muxor2
      tmp     = muxor2 * dx_left - muzor2 * dz_down
      resu_dx(:,1) = resu_dx(:,1) + muxor2 - tmp * 2.0_rp * dx_left * r2
      !
      ! x = Lx wall
      IF (MOD(j1,2) == 0) THEN
         dx_left   = xlen - (x_dip - xlen * j1)
         dx_right  = xlen - (x_dip + xlen * j1)
      ELSE
         dx_left   = xlen - (- x_dip - xlen * (j1-1))
         dx_right  = xlen - (- x_dip + xlen * (j1+1))
      END IF
      dx_left2  = dx_left * dx_left
      dx_right2 = dx_right * dx_right
      ! right x-image on x = Lx wall
      r2      = 1.0_rp / (dx_right2 + dz_up2)
      muxor2  = mu(i1,1) * r2
      muzor2  = mu(i1,2) * r2
      IF (MOD(j1,2) == 1) muxor2 = - muxor2
      tmp     = muxor2 * dx_right + muzor2 * dz_up
      resu_dx(:,2) = resu_dx(:,2) + muxor2 - tmp * 2.0_rp * dx_right * r2
      ! left x-image on x = Lx wall
      r2      = 1.0_rp / (dx_left2 + dz_up2)
      muxor2  = mu(i1,1) * r2
      muzor2  = mu(i1,2) * r2
      IF (MOD(j1,2) == 1) muxor2 = - muxor2
      tmp     = muxor2 * dx_left + muzor2 * dz_up
      resu_dx(:,2) = resu_dx(:,2) + muxor2 - tmp * 2.0_rp * dx_left * r2
      ! right xz-image on x = Lx wall
      r2      = 1.0_rp / (dx_right2 + dz_down2)
      muxor2  = mu(i1,1) * r2
      muzor2  = mu(i1,2) * r2
      IF (MOD(j1,2) == 1) muxor2 = - muxor2
      tmp     = muxor2 * dx_right - muzor2 * dz_down
      resu_dx(:,2) = resu_dx(:,2) + muxor2 - tmp * 2.0_rp * dx_right * r2
      ! left xz-image on x = Lx wall
      r2      = 1.0_rp / (dx_left2 + dz_down2)
      muxor2  = mu(i1,1) * r2
      muzor2  = mu(i1,2) * r2
      IF (MOD(j1,2) == 1) muxor2 = - muxor2
      tmp     = muxor2 * dx_left - muzor2 * dz_down
      resu_dx(:,2) = resu_dx(:,2) + muxor2 - tmp * 2.0_rp * dx_left * r2
   END DO
   !
   eiwt    = EXP(i * TWOPI * wav%freq(i1) * time)
   dphi_dx = dphi_dx + resu_dx * eiwt / TWOPI
END DO
!
dVx0  = REAL(SUM(dphi_dx(:,1)) * alpha) / n_z_chk
dVx0 = dVx0 * delt / xlen
IF (n2 /= 1) dVx0 = dVx0 / ylen
!
dVxLx = - REAL(SUM(dphi_dx(:,2)) * alpha) / n_z_chk
dVxLx = dVxLx * delt / xlen
IF (n2 /= 1) dVxLx = dVxLx / ylen
!
dphi_dzn = 0.0_cp
DO i1 = 1, n_f
   dz_upn     = -1.0_rp - REAL(mu(i1,3)) ! mu(:,3) is real
   dz_up2n    = dz_upn * dz_upn
   dz_downn   = -1.0_rp - REAL(mu(i1,4)) ! mu(:,4) is real
   dz_down2n  = dz_downn * dz_downn
   !
   ! main dipole on z=-1 wall
   dx_rightn  = delta_x(:,:,0)
   dx_right2n = delta_x2(:,:,0)
   r2n      = 1.0_rp / (dx_right2n + dz_up2n)
   muxor2n  = mu(i1,1) * r2n
   muzor2n  = mu(i1,2) * r2n
   tmpn     = muxor2n * dx_rightn + muzor2n * dz_upn
   resu_dzn = muzor2n - tmpn * 2.0_rp * dz_upn * r2n
   ! main dipole z-image
   r2n      = 1.0_rp / (dx_right2n + dz_down2n)
   muxor2n  = mu(i1,1) * r2n
   muzor2n  = mu(i1,2) * r2n
   tmpn     = muxor2n * dx_rightn - muzor2n * dz_downn
   resu_dzn = resu_dzn - muzor2n - tmpn * 2.0_rp * dz_downn * r2n
   !
   DO j1 = 1, n_side ! a double pair of images
      dx_leftn   = delta_x(:,:,-j1)
      dx_left2n  = delta_x2(:,:,-j1)
      dx_rightn  = delta_x(:,:,j1)
      dx_right2n = delta_x2(:,:,j1)
      ! right x-image on z=-1 wall
      r2n      = 1.0_rp / (dx_right2n + dz_up2n)
      muxor2n  = mu(i1,1) * r2n
      muzor2n  = mu(i1,2) * r2n
      IF (MOD(j1,2) == 1) muxor2n = - muxor2n
      tmpn     = muxor2n * dx_rightn + muzor2n * dz_upn
      resu_dzn = resu_dzn + muzor2n - tmpn * 2.0_rp * dz_upn * r2n
      ! left x-image on x = Lx wall
      r2n      = 1.0_rp / (dx_left2n + dz_up2n)
      muxor2n  = mu(i1,1) * r2n
      muzor2n  = mu(i1,2) * r2n
      IF (MOD(j1,2) == 1) muxor2n = - muxor2n
      tmpn     = muxor2n * dx_leftn + muzor2n * dz_upn
      resu_dzn = resu_dzn + muzor2n - tmpn * 2.0_rp * dz_upn * r2n
      ! right xz-image on x = Lx wall
      r2n      = 1.0_rp / (dx_right2n + dz_down2n)
      muxor2n  = mu(i1,1) * r2n
      muzor2n  = mu(i1,2) * r2n
      IF (MOD(j1,2) == 1) muxor2n = - muxor2n
      tmpn     = muxor2n * dx_rightn - muzor2n * dz_downn
      resu_dzn = resu_dzn - muzor2n - tmpn * 2.0_rp * dz_downn * r2n
      ! left xz-image on x = Lx wall
      r2n      = 1.0_rp / (dx_left2n + dz_down2n)
      muxor2n  = mu(i1,1) * r2n
      muzor2n  = mu(i1,2) * r2n
      IF (MOD(j1,2) == 1) muxor2n = - muxor2n
      tmpn     = muxor2n * dx_leftn - muzor2n * dz_downn
      resu_dzn = resu_dzn - muzor2n - tmpn * 2.0_rp * dz_downn * r2n
   END DO
   !
   eiwt    = EXP(i * TWOPI * wav%freq(i1) * time)
   dphi_dzn = dphi_dzn + resu_dzn * eiwt / TWOPI
END DO
dVz = REAL(SUM(dphi_dzn) * alpha) * delx
dVz = dVz * delt / xlen
IF (n2 /= 1) dVz = dVz / ylen
!
!	CPU times outlet
!
if (iCPUtime.eq.1) then
   call CPU_TIME(tf)
   write(*,910)'quitting subroutine phi_dip_chk_flux, total CPU time: ',tf-ti,'s'
endif
910 format(a,1ES11.4,a)
!
END SUBROUTINE phi_dip_chk_flux
!
! start int2 **********************************************************
!      
!     ==========================
      subroutine int2(fy,sumf,i1)
!     ==========================
!
IMPLICIT NONE
!
! INPUT VARIABLES
REAL(RP),DIMENSION(:,:) :: fy
REAL(RP), DIMENSION(:)  :: sumf
INTEGER  :: i1,i2
!
REAL(RP) :: sumy
!
sumy=5.d-1*(fy(i1,1)+fy(i1,n2))
do i2=2,n2-1
   sumy=sumy+fy(i1,i2)
end do
sumf(i1)=sumy
!      
end subroutine int2
!
! end int2 ************************************************************
!
!
!
! start int1 **********************************************************
!           
!     ====================================
      subroutine int1(sumf,sumvol,npt,del)
!     ====================================
!      
IMPLICIT NONE
!
! INPUT VARIABLES
REAL(RP) :: sumf(*),del,sumvol
INTEGER  :: npt
!
INTEGER  :: ii
REAL(RP) :: sum
!            
sum=5.d-1*(sumf(1)+sumf(npt))
do ii=2,npt-1
   sum=sum+sumf(ii)
end do
sumvol=sum*del
!
end subroutine int1
!
! end int1 ************************************************************
!	
END MODULE vol_energy
