MODULE Fourier_FFTW
!
! This module contains all necessary features for the use of FFTW 3.3.4 in HOS-NWT
! with comprehensive interface and normalization
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
USE, INTRINSIC :: iso_c_binding
!
IMPLICIT NONE
!
include 'fftw3.f03'
!
INTERFACE SPACE_2_FOURIER
   MODULE PROCEDURE DCT_FFTW
END INTERFACE SPACE_2_FOURIER
!
INTERFACE FOURIER_2_SPACE
   MODULE PROCEDURE IDCT_FFTW
END INTERFACE FOURIER_2_SPACE
!
INTERFACE FOURIER_2_SPACE_y
   MODULE PROCEDURE IDCT_FFTW_y
END INTERFACE FOURIER_2_SPACE_y
!
INTERFACE SPACE_2_FOURIER_add
   MODULE PROCEDURE DCT_FFTW_add
END INTERFACE SPACE_2_FOURIER_add
! !
INTERFACE FOURIER_2_SPACE_add
   MODULE PROCEDURE IDCT_FFTW_add
END INTERFACE FOURIER_2_SPACE_add
!
! GD : add fourier transforms on 'dealiased' domain i.e. extended
!
INTERFACE SPACE_2_FOURIER_BIG
   MODULE PROCEDURE DCT_FFTW_BIG
END INTERFACE SPACE_2_FOURIER_BIG
!
INTERFACE FOURIER_2_SPACE_BIG
   MODULE PROCEDURE IDCT_FFTW_BIG
END INTERFACE FOURIER_2_SPACE_BIG
!
type(C_PTR)                           :: plan_CC, plan_SC, plan_CS, plan_SS
REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: s_2_f, f_2_s
REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: in, out, in_sin, out_sin
!
type(C_PTR)                          :: plan_Cy, plan_Sy
REAL(RP), ALLOCATABLE, DIMENSION(:)   :: s_2_f_y, f_2_s_y
REAL(RP), ALLOCATABLE, DIMENSION(:)   :: in_y, out_y
!
type(C_PTR)                          :: plan_CC_add, plan_SC_add, plan_CS_add, plan_SS_add
type(C_PTR)                           :: plan_CC_add_I, plan_SC_add_I, plan_CS_add_I, plan_SS_add_I
REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: s_2_f_add, f_2_s_add
REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: in_add, out_add
!
! Fourier transforms on 'dealiased' domain i.e. extended
!
type(C_PTR)                           :: plan_CC_big, plan_SC_big, plan_CS_big, plan_SS_big
REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: s_2_f_big, f_2_s_big
REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: in_big, out_big, in_sin_big, out_sin_big
!
!
CONTAINS
!
!
!
SUBROUTINE Fourier_ini(library)
!
! Initializes the Fourier Transforms
! default library: FFTW-3.x (library=3)
! 
!
IMPLICIT NONE
!
! Input variables
INTEGER, INTENT(IN), OPTIONAL :: library
! Local variables
INTEGER :: lib, flag
!
! Input argument
IF (PRESENT(library)) THEN
   lib = library
ELSE
   ! default library is FFTW3
   lib = 3
END IF
!
! Memory allocation
ALLOCATE(in(n1,n2), out(n1,n2), s_2_f(n1,n2), f_2_s(n1,n2))
ALLOCATE(in_sin(n1-2,n2), out_sin(n1-2,n2))
!
ALLOCATE(in_y(n2), out_y(n2), s_2_f_y(n2), f_2_s_y(n2))
!
ALLOCATE(in_add(n3_add,n2), out_add(n3_add,n2), s_2_f_add(n3_add,n2), f_2_s_add(n3_add,n2))
!
! Fourier transforms on 'dealiased' domain i.e. extended
ALLOCATE(in_big(nd1,nd2), out_big(nd1,nd2), s_2_f_big(nd1,nd2), f_2_s_big(nd1,nd2))
ALLOCATE(in_sin_big(nd1-2,nd2), out_sin_big(nd1-2,nd2))
!
! Initializing the FFT library
SELECT CASE (lib)
   CASE (3)
      ! FFTW-3.x library
      WRITE(*,'(A)') 'Initializing FFTW-3.3.4'
      flag = FFTW_PATIENT ! FFTW_MEASURE ! FFTW_PATIENT ! FFTW_EXHAUSTIVE
      IF (n2 == 1) THEN
         ! 1D FFTs
         CALL dfftw_plan_r2r_1d_(plan_CC, n1,   in(1,1),     out(1,1),     FFTW_REDFT00, flag)
         CALL dfftw_plan_r2r_1d_(plan_SC, n1-2, in_sin(1,1), out_sin(1,1), FFTW_RODFT00, flag)
         !
         CALL dfftw_plan_r2r_1d_(plan_CC_add,   n3_add, in_add(1,1), out_add(1,1), FFTW_REDFT01, flag)
         CALL dfftw_plan_r2r_1d_(plan_CC_add_I, n3_add, in_add(1,1), out_add(1,1), FFTW_REDFT10, flag)
         CALL dfftw_plan_r2r_1d_(plan_SC_add,   n3_add, in_add(1,1), out_add(1,1), FFTW_RODFT01, flag)
         CALL dfftw_plan_r2r_1d_(plan_SC_add_I, n3_add, in_add(1,1), out_add(1,1), FFTW_RODFT10, flag)
	 	 !
         ! Fourier transforms on 'dealiased' domain i.e. extended
         CALL dfftw_plan_r2r_1d(plan_CC_big, nd1,   in_big(1,1),     out_big(1,1),     FFTW_REDFT00, flag)
         CALL dfftw_plan_r2r_1d(plan_SC_big, nd1-2, in_sin_big(1,1), out_sin_big(1,1), FFTW_RODFT00, flag)
      ELSE
         ! 2D FFTs
         CALL dfftw_plan_r2r_2d_(plan_CC, n1,   n2,   in(1,1),     out(1,1),     FFTW_REDFT00, FFTW_REDFT00, flag)
         CALL dfftw_plan_r2r_2d_(plan_SC, n1-2, n2,   in_sin(1,1), out_sin(1,1), FFTW_RODFT00, FFTW_REDFT00, flag)
         CALL dfftw_plan_r2r_2d_(plan_CS, n1,   n2-2, in(1,2),     out(1,2),     FFTW_REDFT00, FFTW_RODFT00, flag)
         CALL dfftw_plan_r2r_2d_(plan_SS, n1-2, n2-2, in_sin(1,2), out_sin(1,2), FFTW_RODFT00, FFTW_RODFT00, flag)
         ! 2D FFTs y-direction only
         CALL dfftw_plan_r2r_1d_(plan_Cy, n2,   in_y(1), out_y(1), FFTW_REDFT00, flag)
         CALL dfftw_plan_r2r_1d_(plan_Sy, n2-2, in_y(2), out_y(2), FFTW_RODFT00, flag)
         ! 2D FFTs additionnal potential
         CALL dfftw_plan_r2r_2d_(plan_CC_add,   n3_add, n2,   in_add(1,1), out_add(1,1), FFTW_REDFT01, FFTW_REDFT00, flag)
         CALL dfftw_plan_r2r_2d_(plan_CC_add_I, n3_add, n2,   in_add(1,1), out_add(1,1), FFTW_REDFT10, FFTW_REDFT00, flag)
         CALL dfftw_plan_r2r_2d_(plan_SC_add,   n3_add, n2,   in_add(1,1), out_add(1,1), FFTW_RODFT01, FFTW_REDFT00, flag)
         CALL dfftw_plan_r2r_2d_(plan_SC_add_I, n3_add, n2,   in_add(1,1), out_add(1,1), FFTW_RODFT10, FFTW_REDFT00, flag)
         CALL dfftw_plan_r2r_2d_(plan_CS_add,   n3_add, n2-2, in_add(1,2), out_add(1,2), FFTW_REDFT01, FFTW_RODFT00, flag)
         CALL dfftw_plan_r2r_2d_(plan_CS_add_I, n3_add, n2-2, in_add(1,2), out_add(1,2), FFTW_REDFT10, FFTW_RODFT00, flag)
         CALL dfftw_plan_r2r_2d_(plan_SS_add,   n3_add, n2-2, in_add(1,2), out_add(1,2), FFTW_RODFT01, FFTW_RODFT00, flag)
         CALL dfftw_plan_r2r_2d_(plan_SS_add_I, n3_add, n2-2, in_add(1,2), out_add(1,2), FFTW_RODFT10, FFTW_RODFT00, flag)
	 	!
	 	! 2D FFTs in 'dealiased' domain (i.e. extended) 
         CALL dfftw_plan_r2r_2d(plan_CC_big, nd1,   nd2,   in_big(1,1),     out_big(1,1),     FFTW_REDFT00, FFTW_REDFT00, flag)
         CALL dfftw_plan_r2r_2d(plan_SC_big, nd1-2, nd2,   in_sin_big(1,1), out_sin_big(1,1), FFTW_RODFT00, FFTW_REDFT00, flag)
         CALL dfftw_plan_r2r_2d(plan_CS_big, nd1,   nd2-2, in_big(1,2),     out_big(1,2),     FFTW_REDFT00, FFTW_RODFT00, flag)
         CALL dfftw_plan_r2r_2d(plan_SS_big, nd1-2, nd2-2, in_sin_big(1,2), out_sin_big(1,2), FFTW_RODFT00, FFTW_RODFT00, flag)
      END IF
      WRITE(*,'(A)') 'Done'
   CASE DEFAULT
      STOP 'Unknown DFT library in Fourier_ini'
END SELECT
!
! Evaluating conversion coefficients
!  when computing a space to Fourier Transform
!  on the free surface
s_2_f           = 1.0_rp / REAL(MAX(1,2*(n1-1)), RP)
s_2_f(2:n1-1,:) = 2.0_rp * s_2_f(2:n1-1,:)
IF (n2 /= 1) THEN
   s_2_f           = s_2_f / REAL(MAX(1,2*(n2-1)), RP)
   s_2_f(:,2:n2-1) = 2.0_rp * s_2_f(:,2:n2-1)
END IF
!
IF (n2 /= 1) THEN
   s_2_f_y         = 1.0_rp / REAL(MAX(1,2*(n2-1)), RP)
   s_2_f_y(2:n2-1) = 2.0_rp * s_2_f_y(2:n2-1)
END IF
!
!  on the wavemaker
s_2_f_add           = 1.0_rp / REAL(n3_add, RP) ! it is 2/(2*n)=2/(2*n3_add)=1/n3_add in fact
IF (n2 /= 1) THEN
   s_2_f_add           = s_2_f_add / REAL(MAX(1,2*(n2-1)), RP)
   s_2_f_add(:,2:n2-1) = 2.0_rp * s_2_f_add(:,2:n2-1)
END IF
!
! on the free surface 'dealiased' domain
s_2_f_big            = 1.0_rp / REAL(MAX(1,2*(nd1-1)), RP)
s_2_f_big(2:nd1-1,:) = 2.0_rp * s_2_f_big(2:nd1-1,:)
IF (n2 /= 1) THEN
   s_2_f_big            = s_2_f_big / REAL(MAX(1,2*(nd2-1)), RP)
   s_2_f_big(:,2:nd2-1) = 2.0_rp * s_2_f_big(:,2:nd2-1)
END IF
!
!  when computing a Fourier to space Transform
!  on the free surface
f_2_s           = 1.0_rp
f_2_s(2:n1-1,:) = 0.5_rp * f_2_s(2:n1-1,:)
IF (n2 /= 1) THEN
   f_2_s(:,2:n2-1) = 0.5_rp * f_2_s(:,2:n2-1)
END IF
!
IF (n2 /= 1) THEN
   f_2_s_y         = 1.0_rp
   f_2_s_y(2:n2-1) = 0.5_rp * f_2_s_y(2:n2-1)
END IF
!
!  on the wavemaker
f_2_s_add = 0.5_rp
IF (n2 /= 1) THEN
   f_2_s_add(:,2:n2-1) = 0.5_rp * f_2_s_add(:,2:n2-1)
END IF
!
! on the free surface 'dealiased' domain
f_2_s_big            = 1.0_rp
f_2_s_big(2:nd1-1,:) = 0.5_rp * f_2_s_big(2:nd1-1,:)
IF (n2 /= 1) THEN
   f_2_s_big(:,2:nd2-1) = 0.5_rp * f_2_s_big(:,2:nd2-1)
END IF
!
END SUBROUTINE Fourier_ini
!
!
!
FUNCTION DCT_FFTW(f, CS_x, CS_y)
!
IMPLICIT NONE
!
! Input variables, output variable
REAL(RP), DIMENSION(m1,m2) :: f, DCT_FFTW
CHARACTER(LEN=*)           :: CS_x, CS_y
!
in   = f(1:n1,1:n2)
!
CALL execute(CS_x, CS_y)
!
DCT_FFTW(1:n1,1:n2) = out * s_2_f
!
END FUNCTION DCT_FFTW
!
!
!
FUNCTION IDCT_FFTW(f, CS_x, CS_y)
!
IMPLICIT NONE
!
! Input variables, output variable
REAL(RP), DIMENSION(m1,m2) :: f, IDCT_FFTW
CHARACTER(LEN=*)           :: CS_x, CS_y
!
in   = f(1:n1,1:n2) * f_2_s
!
CALL execute(CS_x, CS_y)
!
IDCT_FFTW(1:n1,1:n2) = out
!
END FUNCTION IDCT_FFTW
!
!
!
SUBROUTINE execute(CS_x, CS_y)
!
IMPLICIT NONE
!
! Input variables
CHARACTER(LEN=*) :: CS_x, CS_y
! Local variables
type(C_PTR)      :: plan
!
plan = select_plan(CS_x, CS_y)
!
IF ((n2 /= 1) .AND. (CS_y == 'sin')) THEN
   in(:,1)  = 0.0_rp
   in(:,n2) = 0.0_rp
END IF
!
IF (CS_x == 'sin') in_sin = in(2:n1-1,:)
!
call dfftw_execute_(plan)
!
IF (CS_x == 'sin') THEN
   out(1,:)      = 0.0_rp
   out(2:n1-1,:) = out_sin
   out(n1,:)     = 0.0_rp
END IF
!
IF ((n2 /= 1) .AND. (CS_y == 'sin')) THEN
   out(:,1)  = 0.0_rp
   out(:,n2) = 0.0_rp
END IF
!
END SUBROUTINE execute
!
!
!
FUNCTION select_plan(CS_x, CS_y)
!
IMPLICIT NONE
!
! Input variables
CHARACTER(LEN=*)           :: CS_x, CS_y
! Output variable
type(C_PTR)                :: select_plan
!
IF ((CS_x /= 'cos') .AND. (CS_x /= 'sin')) STOP 'Unknown x-transform in Fourier_FFTW:select_plan'
IF ((CS_y /= 'cos') .AND. (CS_y /= 'sin')) STOP 'Unknown y-transform in Fourier_FFTW:select_plan'
!
IF (n2 == 1) THEN
   IF (CS_x == 'cos') THEN
      select_plan = plan_CC
   ELSE IF (CS_x == 'sin') THEN
      select_plan = plan_SC
   END IF
ELSE
   IF ((CS_x == 'cos') .AND. (CS_y == 'cos')) THEN
      select_plan = plan_CC
   ELSE IF ((CS_x == 'sin') .AND. (CS_y == 'cos')) THEN
      select_plan = plan_SC
   ELSE IF ((CS_x == 'cos') .AND. (CS_y == 'sin')) THEN
      select_plan = plan_CS
   ELSE IF ((CS_x == 'sin') .AND. (CS_y == 'sin')) THEN
      select_plan = plan_SS
   END IF
END IF
!
END FUNCTION select_plan
!
!
!
FUNCTION DCT_FFTW_y(f, CS_y)
!
IMPLICIT NONE
!
! Input variables, output variable
REAL(RP), DIMENSION(m2) :: f, DCT_FFTW_y
CHARACTER(LEN=*)        :: CS_y
!
IF (n2 /= 1) THEN
   in_y = f(1:n2)
   !
   CALL execute_y(CS_y)
   !
   DCT_FFTW_y(1:n2) = out_y * s_2_f_y
ELSE
   DCT_FFTW_y = f
END IF
!
END FUNCTION DCT_FFTW_y
!
!
!
FUNCTION IDCT_FFTW_y(f, CS_y)
!
IMPLICIT NONE
!
! Input variables, output variable
REAL(RP), DIMENSION(m2)    :: f, IDCT_FFTW_y
CHARACTER(LEN=*)           :: CS_y
!
IF (n2 /= 1) THEN
   in_y = f(1:n2) * f_2_s_y
   !
   CALL execute_y(CS_y)
   !
   IDCT_FFTW_y(1:n2) = out_y
ELSE
   IDCT_FFTW_y = f
END IF
!
END FUNCTION IDCT_FFTW_y
!
!
!
SUBROUTINE execute_y(CS_y)
!
IMPLICIT NONE
!
! Input variables
CHARACTER(LEN=*) :: CS_y
! Local variables
type(C_PTR)      :: plan
!
plan = select_plan_y(CS_y)
!
IF (CS_y == 'sin') THEN
   in_y(1)  = 0.0_rp
   in_y(n2) = 0.0_rp
END IF
!
CALL dfftw_execute_(plan)
!
IF (CS_y == 'sin') THEN
   out_y(1)  = 0.0_rp
   out_y(n2) = 0.0_rp
END IF
!
END SUBROUTINE execute_y
!
!
!
FUNCTION select_plan_y(CS_y)
!
IMPLICIT NONE
!
! Input variables
CHARACTER(LEN=*) :: CS_y
! Output variable
type(C_PTR)      :: select_plan_y
!
IF ((CS_y /= 'cos') .AND. (CS_y /= 'sin')) STOP 'Unknown y-transform in Fourier_FFTW:select_plan'
!
IF (CS_y == 'cos') THEN
   select_plan_y = plan_Cy
ELSE IF (CS_y == 'sin') THEN
   select_plan_y = plan_Sy
END IF
!
END FUNCTION select_plan_y
!
!
!
FUNCTION DCT_FFTW_add(f, CS_z, CS_y)
!
IMPLICIT NONE
!
! Input variables, output variable
REAL(RP), DIMENSION(m3_add,m2) :: f, DCT_FFTW_add
CHARACTER(LEN=*)               :: CS_z, CS_y
!
in_add = f(1:n3_add,1:n2)
!
CALL execute_add(CS_z, CS_y, 1)
!
DCT_FFTW_add(1:n3_add,1:n2) = out_add * s_2_f_add
!
END FUNCTION DCT_FFTW_add
!
!
!
FUNCTION IDCT_FFTW_add(f, CS_z, CS_y)
!
IMPLICIT NONE
!
! Input variables, output variable
REAL(RP), DIMENSION(m3_add,m2) :: f, IDCT_FFTW_add
CHARACTER(LEN=*)               :: CS_z, CS_y
!
in_add = f(1:n3_add,1:n2) * f_2_s_add
!
CALL execute_add(CS_z, CS_y, -1)
!
IDCT_FFTW_add(1:n3_add,1:n2) = out_add
!
END FUNCTION IDCT_FFTW_add
!
!
!
SUBROUTINE execute_add(CS_z, CS_y, FB)
!
IMPLICIT NONE
!
! Input variables
CHARACTER(LEN=*) :: CS_z, CS_y
INTEGER          :: FB
! Local variables
type(C_PTR)      :: plan
!
plan = select_plan_add(CS_z, CS_y, FB)
!
IF ((n2 /= 1) .AND. (CS_y == 'sin')) THEN
   in_add(:,1)  = 0.0_rp
   in_add(:,n2) = 0.0_rp
END IF
!
CALL dfftw_execute_(plan)
!
IF ((n2 /= 1) .AND. (CS_y == 'sin')) THEN
   out_add(:,1)  = 0.0_rp
   out_add(:,n2) = 0.0_rp
END IF
!
END SUBROUTINE execute_add
!
!
!
FUNCTION select_plan_add(CS_z, CS_y, FB)
!
IMPLICIT NONE
!
! Input variables
CHARACTER(LEN=*) :: CS_z, CS_y
INTEGER          :: FB
! Output variable
type(C_PTR)      :: select_plan_add
!
IF ((CS_z /= 'cos') .AND. (CS_z /= 'sin')) STOP 'Unknown z-transform in Fourier_FFTW:select_plan_add'
IF ((CS_y /= 'cos') .AND. (CS_y /= 'sin')) STOP 'Unknown y-transform in Fourier_FFTW:select_plan_add'
IF ((FB   /= 1)     .AND. (FB   /= -1))    STOP 'Unknown direction in Fourier_FFTW:select_plan_add'
!
IF (n2 == 1) THEN
   IF (FB == 1) THEN
      IF (CS_z == 'cos') THEN
         select_plan_add = plan_CC_add
      ELSE IF (CS_z == 'sin') THEN
         select_plan_add = plan_SC_add
      END IF
   ELSE IF (FB == -1) THEN
      IF (CS_z == 'cos') THEN
         select_plan_add = plan_CC_add_I
      ELSE IF (CS_z == 'sin') THEN
         select_plan_add = plan_SC_add_I
      END IF
   END IF
ELSE
   IF (FB == 1) THEN
      IF ((CS_z == 'cos') .AND. (CS_y == 'cos')) THEN
         select_plan_add = plan_CC_add
      ELSE IF ((CS_z == 'sin') .AND. (CS_y == 'cos')) THEN
         select_plan_add = plan_SC_add
      ELSE IF ((CS_z == 'cos') .AND. (CS_y == 'sin')) THEN
         select_plan_add = plan_CS_add
      ELSE IF ((CS_z == 'sin') .AND. (CS_y == 'sin')) THEN
         select_plan_add = plan_SS_add
      END IF
   ELSE IF (FB == -1) THEN
      IF ((CS_z == 'cos') .AND. (CS_y == 'cos')) THEN
         select_plan_add = plan_CC_add_I
      ELSE IF ((CS_z == 'sin') .AND. (CS_y == 'cos')) THEN
         select_plan_add = plan_SC_add_I
      ELSE IF ((CS_z == 'cos') .AND. (CS_y == 'sin')) THEN
         select_plan_add = plan_CS_add_I
      ELSE IF ((CS_z == 'sin') .AND. (CS_y == 'sin')) THEN
         select_plan_add = plan_SS_add_I
      END IF
   END IF
END IF
!
END FUNCTION select_plan_add
!
! Fourier transforms on 'dealiased' domain i.e. extended
!
FUNCTION DCT_FFTW_big(f, CS_x, CS_y)
!
IMPLICIT NONE
!
! Input variables, output variable
REAL(RP), DIMENSION(md1,md2) :: f, DCT_FFTW_big
CHARACTER(LEN=*)             :: CS_x, CS_y
!
in_big   = f(1:nd1,1:nd2)
!
CALL execute_big(CS_x, CS_y)
!
DCT_FFTW_big(1:nd1,1:nd2) = out_big * s_2_f_big
!
END FUNCTION DCT_FFTW_big
!
!
!
FUNCTION IDCT_FFTW_big(f, CS_x, CS_y)
!
IMPLICIT NONE
!
! Input variables, output variable
REAL(RP), DIMENSION(md1,md2) :: f, IDCT_FFTW_big
CHARACTER(LEN=*)             :: CS_x, CS_y
!
in_big   = f(1:nd1,1:nd2) * f_2_s_big
!
CALL execute_big(CS_x, CS_y)
!
IDCT_FFTW_big(1:nd1,1:nd2) = out_big
!
END FUNCTION IDCT_FFTW_big
!
!
!
SUBROUTINE execute_big(CS_x, CS_y)
!
IMPLICIT NONE
!
! Input variables
CHARACTER(LEN=*) :: CS_x, CS_y
! Local variables
type(C_PTR)      :: plan
!
plan = select_plan_big(CS_x, CS_y)
!
IF ((n2 /= 1) .AND. (CS_y == 'sin')) THEN
   in_big(:,1)  = 0.0_rp
   in_big(:,nd2) = 0.0_rp
END IF
!
IF (CS_x == 'sin') in_sin_big = in_big(2:nd1-1,:)
!
call dfftw_execute_(plan)
!
IF (CS_x == 'sin') THEN
   out_big(1,:)       = 0.0_rp
   out_big(2:nd1-1,:) = out_sin_big
   out_big(nd1,:)     = 0.0_rp
END IF
!
IF ((n2 /= 1) .AND. (CS_y == 'sin')) THEN
   out_big(:,1)   = 0.0_rp
   out_big(:,nd2) = 0.0_rp
END IF
!
END SUBROUTINE execute_big
!
!
!
FUNCTION select_plan_big(CS_x, CS_y)
!
IMPLICIT NONE
!
! Input variables
CHARACTER(LEN=*)           :: CS_x, CS_y
! Output variable
type(C_PTR)                :: select_plan_big
!
IF ((CS_x /= 'cos') .AND. (CS_x /= 'sin')) STOP 'Unknown x-transform in Fourier_FFTW:select_plan'
IF ((CS_y /= 'cos') .AND. (CS_y /= 'sin')) STOP 'Unknown y-transform in Fourier_FFTW:select_plan'
!
IF (n2 == 1) THEN
   IF (CS_x == 'cos') THEN
      select_plan_big = plan_CC_big
   ELSE IF (CS_x == 'sin') THEN
      select_plan_big = plan_SC_big
   END IF
ELSE
   IF ((CS_x == 'cos') .AND. (CS_y == 'cos')) THEN
      select_plan_big = plan_CC_big
   ELSE IF ((CS_x == 'sin') .AND. (CS_y == 'cos')) THEN
      select_plan_big = plan_SC_big
   ELSE IF ((CS_x == 'cos') .AND. (CS_y == 'sin')) THEN
      select_plan_big = plan_CS_big
   ELSE IF ((CS_x == 'sin') .AND. (CS_y == 'sin')) THEN
      select_plan_big = plan_SS_big
   END IF
END IF
!
END FUNCTION select_plan_big
!
!
!
END MODULE Fourier_FFTW
