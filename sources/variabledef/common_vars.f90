MODULE common_vars
!
! This module defines different common variables
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
IMPLICIT NONE
!
! Number of modes
INTEGER, PARAMETER :: n1 = 513
INTEGER, PARAMETER :: n2 = 1
INTEGER, PARAMETER :: n3 = 33
! HOS nonlinearity order
INTEGER, PARAMETER :: mHOS = 3
! Dealiasing parameters
INTEGER, PARAMETER :: p1 = mHOS
INTEGER, PARAMETER :: p2 = 1
! Number of modes with the extended grid
INTEGER, PARAMETER :: Nd1 = ((p1+1) * n1) / 2 ! required modes
INTEGER, PARAMETER :: Nd2 = ((p2+1) * n2) / 2 ! required modes
! Possible filtering in each direction
REAL(RP), DIMENSION(3), PARAMETER :: coeffilt = (/1.0_rp, 1.0_rp, 1.0_rp/)
! Max number of probes in files
INTEGER, PARAMETER :: maxprobes = 10
! Wavemaker definition
INTEGER, PARAMETER :: l_add=2, n3_add = l_add*(n3-1)
! 
INTEGER, PARAMETER :: mHOSdiv2m1 = MAX(1,mHOS / 2 - 1), mHOSm3div2 = MAX(1,(mHOS-3) / 2)
!
! Variables definition for possible specification of n1,n2,n3,mHOS in input files
INTEGER, PARAMETER :: m1=n1, m2=n2, m3=n3, maxm1m2 = max(m1,m2)
INTEGER, PARAMETER :: m3_add = l_add*(m3-1)
INTEGER, PARAMETER :: maxHOS=mHOS, p1max=p1, p2max=p2
INTEGER, PARAMETER :: md1 = (m1*(p1max+1))/2+1, md2 = (m2*(p2max+1))/2+1 !dealiased
INTEGER, PARAMETER :: tecplot = 11, iCPUtime = 0 ! For possible time output of each routine 
!
!
!
END MODULE common_vars