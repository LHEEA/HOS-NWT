MODULE variables
!
! This module defines the different common variables for post-processing
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
!
INTEGER :: n1,n2,n3
REAL(RP), ALLOCATABLE, DIMENSION(:)     :: x,y
REAL(RP), ALLOCATABLE, DIMENSION(:,:)   :: eta
!
REAL(RP), ALLOCATABLE, DIMENSION(:)     :: kx,ky
REAL(RP), ALLOCATABLE, DIMENSION(:,:)   :: k,kth
!
REAL(RP), ALLOCATABLE, DIMENSION(:)     :: kx_add,kx2_add
REAL(RP), ALLOCATABLE, DIMENSION(:,:)   :: k_add,k_add_2,k_add_thk_add
REAL(RP), ALLOCATABLE, DIMENSION(:,:,:) :: csh_add_x,k_add_sh_add_x,kycsh_add_x,kx_add_csh_add_x
!
! Input file
INTEGER               :: i_card, tecplot, i_zvect
REAL(RP)              :: T_start, T_stop, x_min, x_max, y_min, y_max, z_min, z_max
CHARACTER(LEN=100)    :: file_mod
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
INTEGER :: i_unit
INTEGER :: l_add, n3_add
! test for fourier
INTEGER :: m1,m2,m3,m3_add,Nd1,Nd2,md1,md2
!
CONTAINS
!
LOGICAL FUNCTION iseven(n)
!
IMPLICIT NONE
!
INTEGER :: n
!
iseven = (MOD(n,2) == 0)
!
END FUNCTION iseven
!
END MODULE variables