MODULE dealiasing
!
! This module defines all necessary subroutines for correct dealiasing procedure
! Dealiasing may be total or partial w.r.t. to choice of p1 and p2 in computation
!
! Different cases depending if quantity is cos or sin serie for x and y
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
USE fourier_FFTW
!
!
!
CONTAINS
!
!
!
SUBROUTINE dealias(order,todealias,type_x,type_y)
!
IMPLICIT NONE
!
! INPUT VARIABLES
INTEGER, INTENT(IN)          :: order
CHARACTER(LEN=*), INTENT(IN) :: type_x, type_y
! Product quantity to dealias and resulting dealiased quantity
REAL(RP), DIMENSION(md1,md2), INTENT(INOUT) :: todealias
!
! Local variables
REAL(RP) :: ti, tf
!
!	CPU times inlet
IF(iCPUtime.eq.1) then
   print*,'entering subroutine dealias'
   call CPU_TIME(ti)
ENDIF
!
! To prevent from multiples FFTs...
IF (((p1 /= 1).AND.((MOD(order-1,p1-1) == 0))).OR.((p2 /= 1).AND.((MOD(order-1,p2-1) == 0)))) THEN
	!
	todealias = space_2_Fourier_big(todealias,type_x,type_y)
ENDIF
!
! analysis of the quantity to dealias on (Nd1,Nd2) modes and keep only (n1,n2) non-zeros
IF (p1 /= 1) THEN
	IF ((MOD(order-1,p1-1) == 0)) THEN ! partial dealiasing along x-direction
		!
		todealias(n1+1:Nd1,1:Nd2) = 0.0_rp
	END IF
ENDIF
IF ((n2 /= 1).AND.(p2 /= 1)) THEN
	IF ((MOD(order-1,p2-1) == 0)) THEN! partial dealiasing along y-direction
		!
		todealias(1:Nd1,n2+1:Nd2) = 0.0_rp
	ENDIF
ENDIF
IF (((p1 /= 1).AND.((MOD(order-1,p1-1) == 0))).OR.((p2 /= 1).AND.((MOD(order-1,p2-1) == 0)))) THEN
	!
	todealias = Fourier_2_space_big(todealias,type_x,type_y)
ENDIF
!
!	CPU times outlet
IF(iCPUtime.eq.1) then
   call CPU_TIME(tf)
   write(*,910)'quitting subroutine dealias, total CPU time:',tf-ti,'s'
ENDIF
!
910 format(a,1ES11.4,a)
!
RETURN
!
END SUBROUTINE dealias
!
!
!
SUBROUTINE filter_ext(todealias,dealiased,type_x,type_y)
!      
IMPLICIT NONE
!
! INPUT VARIABLES
CHARACTER(LEN=*), INTENT(IN) :: type_x, type_y
! Product quantity to dealias and resulting dealiased quantity
REAL(RP),DIMENSION(md1,md2), INTENT(INOUT) :: todealias
REAL(RP),DIMENSION(m1,m2), INTENT(OUT)     :: dealiased
!
! Local variables
REAL(RP) :: ti, tf
!
!	CPU times inlet
IF(iCPUtime.eq.1) then
   print*,'entering subroutine filter_ext'
   call CPU_TIME(ti)
ENDIF
!
! Direct FFT: determination of the (Nd1,Nd2) 'aliased' modes
todealias = space_2_Fourier_big(todealias,type_x,type_y)
!
! calculation of deliased quantity
dealiased(1:n1,1:n2) = todealias(1:n1,1:n2)
!
! Inverse FFT of dealiased on (n1,n2) points
dealiased = Fourier_2_space(dealiased,type_x,type_y)
!
!	CPU times outlet
IF(iCPUtime.eq.1) then
   call CPU_TIME(tf)
   write(*,910)'quitting subroutine filter_ext, total CPU time: ',tf-ti,'s'
ENDIF
!
910 format(a,1ES11.4,a)
!
RETURN
!
END SUBROUTINE filter_ext
!
!
!
END MODULE dealiasing