MODULE variablechange
!
! This module defines the change of variable 
! for analytical time-integration of linear part of free-surface boundary conditions
!
! Here is defined the exponential of a matrix
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
!
CONTAINS
!
! start expon *********************************************
!
!     ===========================================     
subroutine expon(NL1,NL2,F1,F2,delta_t,iCPUtime)
!     ===========================================
!
! Computation of the exponential of a matrix: (exp(A(t-t0)) cf. p37 PhD thesis G. Ducrozet (ECN)
!
!
IMPLICIT NONE     
!
INTEGER  :: iCPUtime
REAL(RP) :: delta_t
REAL(RP), DIMENSION(m1,m2), INTENT(IN)  :: NL1,NL2
REAL(RP), DIMENSION(m1,m2), INTENT(OUT) :: F1,F2
!
! Local variables
!
INTEGER  :: i1,i2
REAL(RP) :: cosomt,sinomt,ti,tf
!
!	CPU times inlet
!
if (iCPUtime.eq.1) then
	print*,'entering subroutine expon'
	call CPU_TIME(ti)
endif
!
do i1 = 1, n1
   do i2 = 1, n2
      cosomt = cos(omega(i1,i2)*delta_t)
      sinomt = sin(omega(i1,i2)*delta_t)
      F1(i1,i2) = NL1(i1,i2) * cosomt - NL2(i1,i2) * sinomt
      F2(i1,i2) = NL2(i1,i2) * cosomt + NL1(i1,i2) * sinomt
   end do
end do
!      
!	CPU times outlet
!
if (iCPUtime.eq.1) then
	call CPU_TIME(tf)
	write(*,910)'quitting subroutine expon, total CPU time: ',tf-ti,'s'
endif
!
910 format(a,1d12.5,a)
!
end subroutine expon
!
! end expon *********************************************
!
!
! start GtoF *********************************************
!
!     ===========================================     
subroutine GtoF(G1,G2,F1,F2,delta_t,iCPUtime)
!     ===========================================
!
! Transformation G to F cf. p37 PhD thesis G. Ducrozet (ECN)
!
IMPLICIT NONE   
!
INTEGER  :: iCPUtime
REAL(RP) :: delta_t
REAL(RP), DIMENSION(m1,m2), INTENT(IN)  :: G1,G2
REAL(RP), DIMENSION(m1,m2), INTENT(OUT) :: F1,F2
!
! Local variables
!
INTEGER  :: i1,i2
REAL(RP) :: cosomt,sinomt,ti,tf
!
!	CPU times inlet
!
if (iCPUtime.eq.1) then
   print*,'entering subroutine GtoF'
   call CPU_TIME(ti)
endif
!
do i1 = 1, n1
   do i2 = 1, n2
      cosomt = cos(-omega(i1,i2)*delta_t)
      sinomt = sin(-omega(i1,i2)*delta_t)
      F1(i1,i2) = G1(i1,i2) * cosomt - G2(i1,i2) * sinomt
      F2(i1,i2) = G2(i1,i2) * cosomt + G1(i1,i2) * sinomt
   end do
end do
!  
!	CPU times outlet
!
if (iCPUtime.eq.1) then
	call CPU_TIME(tf)
	write(*,910)'quitting subroutine GtoF, total CPU time: ',tf-ti,'s'
endif
!
910 format(a,1d12.5,a)
!
end subroutine GtoF
!
! end GtoF *********************************************
!
! end module
END MODULE variablechange

