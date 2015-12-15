MODULE velocities
!
! This module evaluates the volumic information (potential, velocities, pressure...)
! from the free surface informations
! It uses the H2 operator developed in DNO solver (see Bateman et al.)
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
USE dealiasing
USE resol_HOS
USE variablechange
USE Fourier_FFTW
!
CONTAINS
!
! start HOSvel *************************************************
!
!     ======================================
      subroutine HOSvel_SWEET(vitx_int,vity_int,vitz_int,phit_int)
!     ======================================
!
IMPLICIT NONE
!
REAL(RP),DIMENSION(md1,md2) :: vitx_int, vity_int, vitz_int, phit_int
!
!% LOCAL VARIABLES
!	order m of phi and i^th derivatives at order m
REAL(RP),DIMENSION(md1,md2)        :: phimx_ext,phimy_ext,phimz_ext,phimt_ext
REAL(RP),DIMENSION(m1,m2)          :: phimx, phimy,phimz,phimt
REAL(RP),DIMENSION(md1,md2)        :: phizix,phiziy,phiziz, phizit
REAL(RP),DIMENSION(md1,md2,maxHOS) :: aHOSx_ext, aHOSy_ext, aHOSz_ext, aHOSt_ext
REAL(RP),DIMENSION(md1,md2)        :: vitx_ext, vity_ext, vitz_ext, dphit_ext
!	  
INTEGER  :: i1,i2,iHOS,j,jm1,j2
REAL(RP) :: ti,tf
!
!	CPU times inlet
if(iCPUtime.eq.1) then
   print*,'entering subroutine HOSvel'
   call CPU_TIME(ti)
endif
!
! initializations
do i1=1,nd1
   do i2=1,nd2
      vitx_ext(i1,i2)    = vitx_int(i1,i2)
      vity_ext(i1,i2)    = vity_int(i1,i2)
      vitz_ext(i1,i2)    = vitz_int(i1,i2)
      dphit_ext(i1,i2)   = phit_int(i1,i2)
      aHOSx_ext(i1,i2,1) = vitx_int(i1,i2)
      aHOSy_ext(i1,i2,1) = vity_int(i1,i2)
      aHOSz_ext(i1,i2,1) = vitz_int(i1,i2)
      aHOSt_ext(i1,i2,1) = phit_int(i1,i2)
   end do
enddo
!
call filter_ext(vitx_int,phimx,'sin','cos')
call filter_ext(vity_int,phimy,'cos','sin')
call filter_ext(vitz_int,phimz,'cos','cos')
call filter_ext(phit_int,phimt,'cos','cos')
!
aHOSx_ext(:,:,1) = space_2_Fourier_big(aHOSx_ext(:,:,1),'sin','cos')
aHOSy_ext(:,:,1) = space_2_Fourier_big(aHOSy_ext(:,:,1),'cos','sin')
aHOSz_ext(:,:,1) = space_2_Fourier_big(aHOSz_ext(:,:,1),'cos','cos')
aHOSt_ext(:,:,1) = space_2_Fourier_big(aHOSt_ext(:,:,1),'cos','cos')
!
! calculation at each order iHOS
do iHOS=1,mHOS
   do i1 = 1, nd1
      do i2=1,nd2
         phimx_ext(i1,i2) = 0.0_rp
         phimy_ext(i1,i2) = 0.0_rp
         phimz_ext(i1,i2) = 0.0_rp
         phimt_ext(i1,i2) = 0.0_rp
      end do
   end do
   do j = iHOS,2,-1
      jm1 = j - 1
      j2  = iHOS - jm1  
      do i1=1,nd1
         do i2=1,nd2
            phizix(i1,i2) = aHOSx_ext(i1,i2,j2) * kth_ext(i1,i2,j)
            phiziy(i1,i2) = aHOSy_ext(i1,i2,j2) * kth_ext(i1,i2,j)
            phiziz(i1,i2) = aHOSz_ext(i1,i2,j2) * kth2_ext(i1,i2,j)
            phizit(i1,i2) = aHOSt_ext(i1,i2,j2) * kth_ext(i1,i2,j)
         end do
      enddo
      !
      phizix = Fourier_2_space_big(phizix,'sin','cos')
      phiziy = Fourier_2_space_big(phiziy,'cos','sin')
      phiziz = Fourier_2_space_big(phiziz,'cos','cos')
      phizit = Fourier_2_space_big(phizit,'cos','cos')
      !
      do i1=1,nd1
         do i2=1,nd2
            phimx_ext(i1,i2) = phimx_ext(i1,i2) - phizix(i1,i2) * etapmext(i1,i2,j)
            phimy_ext(i1,i2) = phimy_ext(i1,i2) - phiziy(i1,i2) * etapmext(i1,i2,j)
            phimz_ext(i1,i2) = phimz_ext(i1,i2) - phiziz(i1,i2) * etapmext(i1,i2,j)
            phimt_ext(i1,i2) = phimt_ext(i1,i2) - phizit(i1,i2) * etapmext(i1,i2,j)
         end do
      end do
   end do
   !
   do i1=1,nd1
      do i2=1,nd2
         phizix(i1,i2) = aHOSx_ext(i1,i2,iHOS) * kth_ext(i1,i2,1)
         phiziy(i1,i2) = aHOSy_ext(i1,i2,iHOS) * kth_ext(i1,i2,1)
         phiziz(i1,i2) = aHOSz_ext(i1,i2,iHOS) * kth2_ext(i1,i2,1)
         phizit(i1,i2) = aHOSt_ext(i1,i2,iHOS) * kth_ext(i1,i2,1)
      end do
   enddo
   !
   phizix = Fourier_2_space_big(phizix,'sin','cos')
   phiziy = Fourier_2_space_big(phiziy,'cos','sin')
   phiziz = Fourier_2_space_big(phiziz,'cos','cos')
   phizit = Fourier_2_space_big(phizit,'cos','cos')
   !
   do i1=1,nd1
      do i2=1,nd2
         phimx_ext(i1,i2)= phimx_ext(i1,i2) - phizix(i1,i2) * etapmext(i1,i2,1)
         phimy_ext(i1,i2)= phimy_ext(i1,i2) - phiziy(i1,i2) * etapmext(i1,i2,1)
         phimz_ext(i1,i2)= phimz_ext(i1,i2) - phiziz(i1,i2) * etapmext(i1,i2,1)
         phimt_ext(i1,i2)= phimt_ext(i1,i2) - phizit(i1,i2) * etapmext(i1,i2,1)
         !
         !	construction of velocities
         vitx_ext(i1,i2)  = vitx_ext(i1,i2)  + phimx_ext(i1,i2)
         vity_ext(i1,i2)  = vity_ext(i1,i2)  + phimy_ext(i1,i2)
         vitz_ext(i1,i2)  = vitz_ext(i1,i2)  + phimz_ext(i1,i2)
         dphit_ext(i1,i2) = dphit_ext(i1,i2) + phimt_ext(i1,i2)
      end do
   end do
   !
   CALL dealias(iHOS+1, vitx_ext, 'sin', 'cos')
   CALL dealias(iHOS+1, vity_ext,'cos','sin')
   CALL dealias(iHOS+1, vitz_ext, 'cos', 'cos')
   CALL dealias(iHOS+1, dphit_ext, 'cos', 'cos')
   !
   IF (iHOS < mHOS) THEN
	  CALL dealias(iHOS+1, phimx_ext,'sin','cos')
      CALL dealias(iHOS+1, phimy_ext,'cos','sin')
      CALL dealias(iHOS+1, phimz_ext, 'cos', 'cos')
      CALL dealias(iHOS+1, phimt_ext, 'cos', 'cos')
      ! Storing it for next i_m
      aHOSx_ext(:,:,iHOS+1) = phimx_ext(:,:)
      aHOSy_ext(:,:,iHOS+1) = phimy_ext(:,:)
      aHOSz_ext(:,:,iHOS+1) = phimz_ext(:,:)
      aHOSt_ext(:,:,iHOS+1) = phimt_ext(:,:)
      !
      aHOSx_ext(:,:,iHOS+1) = space_2_Fourier_big(aHOSx_ext(:,:,iHOS+1),'sin','cos')
      aHOSy_ext(:,:,iHOS+1) = space_2_Fourier_big(aHOSy_ext(:,:,iHOS+1),'cos','sin')
      aHOSz_ext(:,:,iHOS+1) = space_2_Fourier_big(aHOSz_ext(:,:,iHOS+1),'cos','cos')
      aHOSt_ext(:,:,iHOS+1) = space_2_Fourier_big(aHOSt_ext(:,:,iHOS+1),'cos','cos')
   END IF
end do
!
call filter_ext(vitx_ext,modesspecx,'sin','cos')
call filter_ext(vity_ext,modesspecy,'cos','sin')
call filter_ext(vitz_ext,modesspecz,'cos','cos')
call filter_ext(dphit_ext,modesspect,'cos','cos')
!
modesspecx = space_2_Fourier(modesspecx,'sin','cos')
modesspecy = space_2_Fourier(modesspecy,'cos','sin')
modesspecz = space_2_Fourier(modesspecz,'cos','cos')
modesspect = space_2_Fourier(modesspect,'cos','cos')
!
! Filtering done on the modal description of volumic informations
!
call filtering_x_modes(modesspecx,coeffilt(1))
call filtering_z_modes(modesspecx,coeffilt(3))
!
call filtering_x_modes(modesspecy,coeffilt(1))
call filtering_z_modes(modesspecy,coeffilt(3))
!
call filtering_x_modes(modesspecz,coeffilt(1))
call filtering_z_modes(modesspecz,coeffilt(3))
!
call filtering_x_modes(modesspect,coeffilt(1))
call filtering_z_modes(modesspect,coeffilt(3))
!
IF (n2 /= 1) THEN
	call filtering_y_modes(modesspecx,coeffilt(2))
	call filtering_y_modes(modesspecy,coeffilt(2))
	call filtering_y_modes(modesspecz,coeffilt(2))
	call filtering_y_modes(modesspect,coeffilt(2))
ENDIF
!
!	CPU times outlet
if (iCPUtime.eq.1) then
   call CPU_TIME(tf)
   write(*,910)'quitting subroutine HOSvel, total CPU time:',tf-ti,'s'
endif
!
910 format(a,1ES11.4,a)
!
end subroutine HOSvel_SWEET
!
! end HOSvel ***************************************************
!
! start HOSvel2_SWEET *************************************************
!
!     ======================================
      subroutine HOSvel2_SWEET(meta2,eta_l,phis_l,a1st_add_l,da1st_add_l,a2nd_add_l,da2nd_add_l,a3rd_add_l,da3rd_add_l, &
           a1st_loc,da1st_loc, eta1st_loc, deta1st_loc, a2nd_loc, da2nd_loc, eta2nd_loc, deta2nd_loc, &
           G2_loc, dG2_loc, G1_loc, dG1_loc, da_add_ramp_loc, h_loc)
!     ======================================
!
!USE filtering
IMPLICIT NONE
!
!% INPUT VARIABLES
!	Number of levels for H2 operator
INTEGER	:: meta2
!
REAL(RP),DIMENSION(m1,m2)     :: eta_l,phis_l,a1st_loc,da1st_loc,eta1st_loc,deta1st_loc,a2nd_loc,da2nd_loc,eta2nd_loc
REAL(RP),DIMENSION(m3_add,m2) :: a1st_add_l,da1st_add_l,a2nd_add_l,da2nd_add_l,a3rd_add_l,da3rd_add_l,da_add_ramp_loc
REAL(RP),DIMENSION(m1,m2)     :: deta2nd_loc,G2_loc,dG2_loc,G1_loc,dG1_loc
REAL(RP) :: h_loc
!
!% LOCAL VARIABLES
REAL(RP),DIMENSION(m1,m2)          :: eta2,as2
REAL(RP),DIMENSION(m1,m2)          :: phimx,phimy,phimz,phimt
REAL(RP),DIMENSION(md1,md2,maxHOS) :: etapmext2, aHOSx_ext, aHOSy_ext, aHOSz_ext, aHOSt_ext
!
REAL(RP),DIMENSION(md1,md2) :: phizix,phiziy,phiziz,phizit
REAL(RP),DIMENSION(md1,md2) :: phimx_ext,phimy_ext,phimz_ext,phimt_ext
REAL(RP),DIMENSION(md1,md2) :: phimx2_ext, phimy2_ext, phimz2_ext,phimt2_ext
!
!
REAL(RP),DIMENSION(md1,md2)   :: vitx_int, vity_int, vitz_int, phit_int
!	  
INTEGER  :: i1, i2, iHOS, j, jm1, jp1, j2, ii1
REAL(RP) :: ti,tf
!
!	CPU times inlet
if(iCPUtime.eq.1) then
   print*,'entering subroutine HOSvel2_SWEET'
   call CPU_TIME(ti)
endif
!
! Initialisation of velocity computations
CALL init_velocity_comp(eta_l,phis_l,da1st_add_l,a1st_add_l,da2nd_add_l,a2nd_add_l,da3rd_add_l,a3rd_add_l, &
     	a1st_loc,da1st_loc,eta1st_loc,deta1st_loc,a2nd_loc,da2nd_loc,eta2nd_loc,deta2nd_loc,&
     	G2_loc,dG2_loc,G1_loc,dG1_loc,da_add_ramp_loc,h_loc,phimx,phimy,phimz,phimt)
! FIXME: CHECK where it is computed
do i1 = 1,nd1
   do i2 = 1,nd2
      etapmext2(i1,i2,1:mHOS)=etapmext(i1,i2,1:mHOS)
   end do
enddo
!
! Extension of quantities on dealiased dimension
phimx = space_2_Fourier(phimx,'sin','cos')
phimy = space_2_Fourier(phimy,'cos','sin')
phimz = space_2_Fourier(phimz,'cos','cos')
phimt = space_2_Fourier(phimt,'cos','cos')
!
vitx_int = extend(phimx)
vity_int = extend(phimy)
vitz_int = extend(phimz)
phit_int = extend(phimt)
!
aHOSx_ext(:,:,1) = vitx_int
aHOSy_ext(:,:,1) = vity_int
aHOSz_ext(:,:,1) = vitz_int
aHOSt_ext(:,:,1) = phit_int
!
phimx = Fourier_2_space(phimx,'sin','cos')
phimy = Fourier_2_space(phimy,'cos','sin')
phimz = Fourier_2_space(phimz,'cos','cos')
phimt = Fourier_2_space(phimt,'cos','cos')
!
vitx_int = Fourier_2_space_big(vitx_int,'sin','cos')
vity_int = Fourier_2_space_big(vity_int,'cos','sin')
vitz_int = Fourier_2_space_big(vitz_int,'cos','cos')
phit_int = Fourier_2_space_big(phit_int,'cos','cos')
!
! H2 operator (transforms with meta2 steps quantities from z=eta to z=0)
do ii1=1,meta2-1
   !
   do i1=1,nd1
      do i2=1,nd2
         etapmext(i1,i2,:)=etapmext2(i1,i2,:)
      enddo
   enddo
   !
   do i1=1,n1
      do i2=1,n2
         eta2(i1,i2)=eta_l(i1,i2)*(meta2-1-(ii1-1))/(meta2)
      enddo
   enddo
   !   
   ! inverse FFT: determination of the modes of eta2
   do i1=1,n1
      do i2=1,n2
         as2(i1,i2)=eta2(i1,i2) 
      enddo
   end do
   !
   as2 = space_2_Fourier(as2,'cos','cos')
   !
   do i1=1,n1
      do i2=1,n2
         etapmext2(i1,i2,1)=as2(i1,i2)
      end do
      do i2=n2+1,nd2
         etapmext2(i1,i2,1)=0.0_rp
      end do
   end do
   do i1=n1+1,nd1
      do i2=1,nd2
         etapmext2(i1,i2,1)=0.0_rp
      end do
   end do
   !
   etapmext2(:,:,1) = Fourier_2_space_big(etapmext2(:,:,1),'cos','cos')
   !
   !	de-aliasing of required powers of eta2
   do j=1,mHOS-1
      jp1 = j+1
      do i1=1,nd1
         do i2=1,nd2
            etapmext2(i1,i2,jp1) = etapmext2(i1,i2,j) / jp1 * etapmext2(i1,i2,1)
         end do
      end do
      call dealias(j,etapmext2(:,:,jp1),'cos','cos')
   end do
   !
   ! calculation at each order iHOS
   do iHOS=1,mHOS
      !
      do i1 = 1, nd1
         do i2=1,nd2
            phimx_ext(i1,i2)  = 0.0_rp
            phimy_ext(i1,i2)  = 0.0_rp
            phimz_ext(i1,i2)  = 0.0_rp
            phimt_ext(i1,i2)  = 0.0_rp
            phimx2_ext(i1,i2) = 0.0_rp
            phimy2_ext(i1,i2) = 0.0_rp
            phimz2_ext(i1,i2) = 0.0_rp
            phimt2_ext(i1,i2) = 0.0_rp
         end do
      end do
      !
      do j = iHOS,2,-1
         jm1 = j - 1
         j2  = iHOS - jm1  
         !
         do i1=1,nd1
            do i2=1,nd2
               phizix(i1,i2) = aHOSx_ext(i1,i2,j2) * kth_ext(i1,i2,j)
               phiziy(i1,i2) = aHOSy_ext(i1,i2,j2) * kth_ext(i1,i2,j)
               phiziz(i1,i2) = aHOSz_ext(i1,i2,j2) * kth2_ext(i1,i2,j)
               phizit(i1,i2) = aHOSt_ext(i1,i2,j2) * kth_ext(i1,i2,j)
            end do
         enddo
         !
         phizix = Fourier_2_space_big(phizix,'sin','cos')
         phiziy = Fourier_2_space_big(phiziy,'cos','sin')
         phiziz = Fourier_2_space_big(phiziz,'cos','cos')
         phizit = Fourier_2_space_big(phizit,'cos','cos')
         !
         do i1=1,nd1
            do i2=1,nd2
               phimx_ext(i1,i2)  = phimx_ext(i1,i2)  - phizix(i1,i2) * etapmext(i1,i2,j)
               phimy_ext(i1,i2)  = phimy_ext(i1,i2)  - phiziy(i1,i2) * etapmext(i1,i2,j)
               phimz_ext(i1,i2)  = phimz_ext(i1,i2)  - phiziz(i1,i2) * etapmext(i1,i2,j)
               phimt_ext(i1,i2)  = phimt_ext(i1,i2)  - phizit(i1,i2) * etapmext(i1,i2,j)
               phimx2_ext(i1,i2) = phimx2_ext(i1,i2) + phizix(i1,i2) * etapmext2(i1,i2,j)
               phimy2_ext(i1,i2) = phimy2_ext(i1,i2) + phiziy(i1,i2) * etapmext2(i1,i2,j)
               phimz2_ext(i1,i2) = phimz2_ext(i1,i2) + phiziz(i1,i2) * etapmext2(i1,i2,j)
               phimt2_ext(i1,i2) = phimt2_ext(i1,i2) + phizit(i1,i2) * etapmext2(i1,i2,j)
            end do
         end do
      end do
      !
      do i1=1,nd1
         do i2=1,nd2
            phizix(i1,i2) = aHOSx_ext(i1,i2,iHOS) * kth_ext(i1,i2,1)
            phiziy(i1,i2) = aHOSy_ext(i1,i2,iHOS) * kth_ext(i1,i2,1)
            phiziz(i1,i2) = aHOSz_ext(i1,i2,iHOS) * kth2_ext(i1,i2,1)
            phizit(i1,i2) = aHOSt_ext(i1,i2,iHOS) * kth_ext(i1,i2,1)
         end do
      enddo
      !
      phizix = Fourier_2_space_big(phizix,'sin','cos')
      phiziy = Fourier_2_space_big(phiziy,'cos','sin')
      phiziz = Fourier_2_space_big(phiziz,'cos','cos')
      phizit = Fourier_2_space_big(phizit,'cos','cos')
      do i1=1,nd1
         do i2=1,nd2
            phimx_ext(i1,i2) = phimx_ext(i1,i2) - phizix(i1,i2) * etapmext(i1,i2,1) 
            phimy_ext(i1,i2) = phimy_ext(i1,i2) - phiziy(i1,i2) * etapmext(i1,i2,1)
            phimz_ext(i1,i2) = phimz_ext(i1,i2) - phiziz(i1,i2) * etapmext(i1,i2,1)
            phimt_ext(i1,i2) = phimt_ext(i1,i2) - phizit(i1,i2) * etapmext(i1,i2,1)
            phimx2_ext(i1,i2) = phimx2_ext(i1,i2) + phizix(i1,i2) * etapmext2(i1,i2,1) 
            phimy2_ext(i1,i2) = phimy2_ext(i1,i2) + phiziy(i1,i2) * etapmext2(i1,i2,1)
            phimz2_ext(i1,i2) = phimz2_ext(i1,i2) + phiziz(i1,i2) * etapmext2(i1,i2,1)
            phimt2_ext(i1,i2) = phimt2_ext(i1,i2) + phizit(i1,i2) * etapmext2(i1,i2,1)
            !
            !	construction of velocities
            vitx_int(i1,i2) = vitx_int(i1,i2) + phimx_ext(i1,i2) + phimx2_ext(i1,i2)
            vity_int(i1,i2) = vity_int(i1,i2) + phimy_ext(i1,i2) + phimy2_ext(i1,i2)
            vitz_int(i1,i2) = vitz_int(i1,i2) + phimz_ext(i1,i2) + phimz2_ext(i1,i2)
            phit_int(i1,i2) = phit_int(i1,i2) + phimt_ext(i1,i2) + phimt2_ext(i1,i2)
         end do
      end do
      !
      CALL dealias(iHOS, vitx_int, 'sin', 'cos')
      CALL dealias(iHOS, vity_int, 'cos', 'sin')
      CALL dealias(iHOS, vitz_int, 'cos', 'cos')
      CALL dealias(iHOS, phit_int, 'cos', 'cos')
      !
      IF (iHOS < mHOS) THEN
         !
         CALL dealias(iHOS+1, phimx_ext, 'sin', 'cos')
         CALL dealias(iHOS+1, phimy_ext, 'cos', 'sin')
         CALL dealias(iHOS+1, phimz_ext, 'cos', 'cos')
         CALL dealias(iHOS+1, phimt_ext, 'cos', 'cos')
         !
         ! Storing it for next i_m
         aHOSx_ext(:,:,iHOS+1) = phimx_ext(:,:)
         aHOSy_ext(:,:,iHOS+1) = phimy_ext(:,:)
         aHOSz_ext(:,:,iHOS+1) = phimz_ext(:,:)
         aHOSt_ext(:,:,iHOS+1) = phimt_ext(:,:)
         !
         aHOSx_ext(:,:,iHOS+1) = space_2_Fourier_big(aHOSx_ext(:,:,iHOS+1),'sin','cos')
         aHOSy_ext(:,:,iHOS+1) = space_2_Fourier_big(aHOSy_ext(:,:,iHOS+1),'cos','sin')
         aHOSz_ext(:,:,iHOS+1) = space_2_Fourier_big(aHOSz_ext(:,:,iHOS+1),'cos','cos')
         aHOSt_ext(:,:,iHOS+1) = space_2_Fourier_big(aHOSt_ext(:,:,iHOS+1),'cos','cos')
      END IF
   end do
   !
   do i1=1,nd1
      do i2=1,nd2
         phimx_ext(i1,i2) = vitx_int(i1,i2)
         phimy_ext(i1,i2) = vity_int(i1,i2)
         phimz_ext(i1,i2) = vitz_int(i1,i2)
         phimt_ext(i1,i2) = phit_int(i1,i2)
      end do
   enddo
   !
   ! Storing it for next eta2
   aHOSx_ext(:,:,1) = phimx_ext(:,:)
   aHOSy_ext(:,:,1) = phimy_ext(:,:)
   aHOSz_ext(:,:,1) = phimz_ext(:,:)
   aHOSt_ext(:,:,1) = phimt_ext(:,:)
   !
   aHOSx_ext(:,:,1) = space_2_Fourier_big(aHOSx_ext(:,:,1),'sin','cos')
   aHOSy_ext(:,:,1) = space_2_Fourier_big(aHOSy_ext(:,:,1),'cos','sin')
   aHOSz_ext(:,:,1) = space_2_Fourier_big(aHOSz_ext(:,:,1),'cos','cos')
   aHOSt_ext(:,:,1) = space_2_Fourier_big(aHOSt_ext(:,:,1),'cos','cos')
end do
!
do i1=1,nd1
   do i2=1,nd2
      etapmext(i1,i2,1:mHOS)=etapmext2(i1,i2,1:mHOS)
   end do
enddo
!
! Last step to transform to z=0
CALL HOSvel_SWEET(vitx_int,vity_int,vitz_int,phit_int)
!
!	CPU times outlet
if (iCPUtime.eq.1) then
   call CPU_TIME(tf)
   write(*,910)'quitting subroutine HOSvel2_SWEET, total CPU time:',tf-ti,'s'
endif
!
910 format(a,1ES11.4,a)
!
end subroutine HOSvel2_SWEET
!
! end HOSvel2_bis ***************************************************
!
!
!
! start init_velocity_comp ***************************************************
!
!    ======================================
     subroutine init_velocity_comp(eta_l,phis_l,da1st_add_l,a1st_add_l,da2nd_add_l,a2nd_add_l,da3rd_add_l,a3rd_add_l, &
     	a1st_loc,da1st_loc,eta1st_loc,deta1st_loc,a2nd_loc,da2nd_loc,eta2nd_loc,deta2nd_loc,&
     	G2_loc,dG2_loc,G1_loc,dG1_loc,da_add_ramp_loc,h_loc,phimx,phimy,phimz,phimt)
!    ======================================
!
IMPLICIT NONE
!
! Input variables
REAL(RP), DIMENSION(m1,m2), INTENT(IN)     :: eta_l,phis_l
REAL(RP), DIMENSION(m3_add,m2), INTENT(IN) :: da1st_add_l,a1st_add_l,da2nd_add_l,a2nd_add_l,da3rd_add_l,a3rd_add_l,da_add_ramp_loc
REAL(RP), DIMENSION(m1,m2), INTENT(IN)     :: a1st_loc,da1st_loc,eta1st_loc,deta1st_loc,a2nd_loc,da2nd_loc,eta2nd_loc,deta2nd_loc
REAL(RP), DIMENSION(m1,m2), INTENT(IN)     :: G2_loc,dG2_loc,G1_loc,dG1_loc
REAL(RP), INTENT(IN) :: h_loc
!
! Output variables
REAL(RP), DIMENSION(m1,m2), INTENT(OUT)    :: phimx,phimy,phimz,phimt
!
! Local variables
!INTEGER :: i1,i2
REAL(RP), DIMENSION(m3_add,m2) :: da1st_add_int,a1st_add_int,da2nd_add_int,a2nd_add_int,da3rd_add_int,a3rd_add_int,da_add_ramp_int
REAL(RP), DIMENSION(m1,m2)     :: a1st_int,da1st_int,eta1st_int,deta1st_int,a2nd_int,da2nd_int,eta2nd_int,deta2nd_int
REAL(RP), DIMENSION(m1,m2)     :: G2_int,dG2_int,G1_int,dG1_int
REAL(RP), DIMENSION(m1,m2)     :: eta3,phis3,deta3,dphis3
REAL(RP) :: dphis3_1
!
REAL(RP),DIMENSION(md1,md2) :: phimx_ext,phimy_ext
REAL(RP),DIMENSION(md1,md2) :: phimx2_ext, phimy2_ext, phimz2_ext
!
! Verification of FSBCs
REAL(RP),DIMENSION(md1,md2)  :: detaphiz,phimxetax_ext,etax2, etay2,phiz_etax_ext,phiz_etay_ext,phizMm2_etax_ext,phizMm2_etay_ext
REAL(RP),DIMENSION(m1,m2)    :: phimx2, phimy2, phimz2,detaphiz_d,phimxetax, phiz_etax, phiz_etay
REAL(RP),DIMENSION(m1,m2)    :: CCSL, CCSL_bis, Press_bis, Press
!
! temporary variables for etax and etay on initial grid (n1,n2)
REAL(RP),DIMENSION(m1,m2)     :: etax_tmp, etay_tmp
!	  
INTEGER  :: i1, i2, i3, isubsteprk
REAL(RP) :: c_end
! dealiasing add. part
REAL(RP),DIMENSION(md1,md2) :: phix_add_ext, phiy_add_ext, phiz_add_ext, phiz_ext, tmp1_add_ext, tmp2_add_ext
REAL(RP),DIMENSION(m1,m2)   :: tmp1_add, tmp2_add
!
! Initialization
do i1=1,n1
   do i2=1,n2
      eta3(i1,i2) = eta_l(i1,i2)
      phis3(i1,i2)= phis_l(i1,i2)
   end do
enddo
!
do i1=1,n3_add
   do i2=1,n2
      da1st_add_int(i1,i2) = da1st_add_l(i1,i2)
      a1st_add_int(i1,i2)  = a1st_add_l(i1,i2)
      da2nd_add_int(i1,i2) = da2nd_add_l(i1,i2)
      a2nd_add_int(i1,i2)  = a2nd_add_l(i1,i2)
      da3rd_add_int(i1,i2) = da3rd_add_l(i1,i2)
      a3rd_add_int(i1,i2)  = a3rd_add_l(i1,i2)
   end do
enddo
!
a1st_int(1:n1,1:n2)    = a1st_loc(1:n1,1:n2)
da1st_int(1:n1,1:n2)   = da1st_loc(1:n1,1:n2)
eta1st_int(1:n1,1:n2)  = eta1st_loc(1:n1,1:n2)
deta1st_int(1:n1,1:n2) = deta1st_loc(1:n1,1:n2)
a2nd_int(1:n1,1:n2)    = a2nd_loc(1:n1,1:n2)
da2nd_int(1:n1,1:n2)   = da2nd_loc(1:n1,1:n2)
eta2nd_int(1:n1,1:n2)  = eta2nd_loc(1:n1,1:n2)
deta2nd_int(1:n1,1:n2) = deta2nd_loc(1:n1,1:n2)
G2_int(1:n1,1:n2)      = G2_loc(1:n1,1:n2)
dG2_int(1:n1,1:n2)     = dG2_loc(1:n1,1:n2)
G1_int(1:n1,1:n2)      = G1_loc(1:n1,1:n2)
dG1_int(1:n1,1:n2)     = dG1_loc(1:n1,1:n2)
da_add_ramp_int(1:n3_add,1:n2) = da_add_ramp_loc(1:n3_add,1:n2)
!
! Evaluate local time derivatives
IF ((igeom==1.OR.igeom==2).AND.&
	(icase==2.OR.icase==3.OR.icase==31.OR.icase==32.OR.icase==33.OR.icase==4.OR.icase==41.OR.icase==42)) THEN
   isubsteprk = 0
   call solveHOS2(a1st_int, da1st_int, eta1st_int, deta1st_int,a2nd_int,da2nd_int, &
        eta2nd_int, deta2nd_int,G2_int, dG2_int, G1_int, dG1_int, &
        a1st_add_int,a2nd_add_int,da1st_add_int,da2nd_add_int,da_add_ramp_int, &
        a3rd_add_int,da3rd_add_int, isubsteprk, h_loc)
   dphis3_1 = dG2_int(1,1)
   dG2_int(1,1) = 0.0_rp
   !
   eta3  = space_2_Fourier(eta3,'cos','cos')
   phis3 = space_2_Fourier(phis3,'cos','cos')
   !
   call GtoF(dG1_int,dG2_int,deta3,dphis3,h_loc,iCPUtime)
   !
   dphis3(1,1)   = dphis3_1 - eta3(1,1)
   modesFSt(1,1) = deta3(1,1) 
   !
   ! FIXME: missing 
   do i2=2,n2
     dphis3(1,i2)   = dphis3(1,i2)/omega(1,i2) - eta3(1,i2)
     deta3(1,i2)    = deta3(1,i2) + phis3(1,i2)*omega(1,i2)**2
     modesFSt(1,i2) = deta3(1,i2)
   end do
   do i1=2,n1
      do i2=1,n2
         dphis3(i1,i2)   = dphis3(i1,i2)/omega(i1,i2) - eta3(i1,i2)
         deta3(i1,i2)    = deta3(i1,i2) + phis3(i1,i2)*omega(i1,i2)**2
         modesFSt(i1,i2) = deta3(i1,i2)
      end do
   end do
   !
   eta3   = Fourier_2_space(eta3,  'cos','cos')
   phis3  = Fourier_2_space(phis3, 'cos','cos')
   deta3  = Fourier_2_space(deta3, 'cos','cos')
   dphis3 = Fourier_2_space(dphis3,'cos','cos')
   !
ELSE
   isubsteprk = 0
   call solveHOS(phis3,dphis3,eta3,deta3,a1st_add_int,da1st_add_int,isubsteprk,c_end)
   !
   modesFSt = deta3
   modesFSt = space_2_Fourier(modesFSt,'cos','cos')
END IF
!
! Additional modes
do i3=1,n3_add
   do i2=1,n2
      modesadd(i3,i2)  = a1st_add_l(i3,i2) + a2nd_add_l(i3,i2) + a3rd_add_l(i3,i2)
      modesaddt(i3,i2) = da1st_add_int(i3,i2) + da2nd_add_int(i3,i2) + da3rd_add_int(i3,i2)
   end do
enddo
!
! de-aliasing of terms involved in multiple products
do i1 = 1,nd1
   do i2 = 1,nd2
      etax2(i1,i2) = etax(i1,i2) * etax(i1,i2)
      etay2(i1,i2) = etay(i1,i2) * etay(i1,i2)
   end do
end do
!
call dealias(2,etax2,'cos','cos')
call dealias(2,etay2,'cos','cos')
!
! initialization of H2 calculations
do i1 = 1,nd1
   do i2 = 1,nd2
      phiz_etax_ext(i1,i2) = etax(i1,i2) * phizMm1(i1,i2)
      phiz_etay_ext(i1,i2) = etay(i1,i2) * phizMm1(i1,i2)
      !
      phizMm2_etax_ext(i1,i2) = etax(i1,i2) * phizMm2(i1,i2)
      phizMm2_etay_ext(i1,i2) = etay(i1,i2) * phizMm2(i1,i2)
   end do
end do
!
call dealias(2,phiz_etax_ext,'sin','cos')
call dealias(2,phiz_etay_ext,'cos','sin')
!
call dealias(2,phizMm2_etax_ext,'sin','cos')
call dealias(2,phizMm2_etay_ext,'cos','sin')
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Test dealiasing add part
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
phix_add = space_2_Fourier(phix_add,'cos','cos')
phiy_add = space_2_Fourier(phiy_add,'cos','cos')
phiz_add = space_2_Fourier(phiz_add,'cos','cos')
phimz    = space_2_Fourier(phiz,'cos','cos')
!
phix_add_ext = extend(phix_add)
phiy_add_ext = extend(phiy_add)
phiz_add_ext = extend(phiz_add)
phiz_ext     = extend(phimz)
!
phix_add = Fourier_2_space(phix_add,'cos','cos')
phiy_add = Fourier_2_space(phiy_add,'cos','cos')
phiz_add = Fourier_2_space(phiz_add,'cos','cos')
phix_add_ext = Fourier_2_space_big(phix_add_ext,'cos','cos')
phiy_add_ext = Fourier_2_space_big(phiy_add_ext,'cos','cos')
phiz_add_ext = Fourier_2_space_big(phiz_add_ext,'cos','cos')
!
phimz    = Fourier_2_space(phimz,'cos','cos')
phiz_ext = Fourier_2_space_big(phiz_ext,'cos','cos')
!
do i1 = 1,nd1
   do i2 = 1,nd2
      phimx_ext(i1,i2)     = phisx(i1,i2) - etax(i1,i2) * phizMm1(i1,i2)
      phimy_ext(i1,i2)     = phisy(i1,i2) - etay(i1,i2) * phizMm1(i1,i2)
      phimx2_ext(i1,i2)    = phisx(i1,i2)*phisx(i1,i2) + etax2(i1,i2)*phiz2Mm2(i1,i2) &
                                - 2.0_rp*phisx(i1,i2)*phizMm2_etax_ext(i1,i2) 
      phimy2_ext(i1,i2)    = phisy(i1,i2)*phisy(i1,i2) + etay2(i1,i2)*phiz2Mm2(i1,i2) &
                                - 2.0_rp*phisy(i1,i2)*phizMm2_etay_ext(i1,i2)
      phimz2_ext(i1,i2)    = phiz2(i1,i2)
      phimxetax_ext(i1,i2) = geta2phiz(i1,i2) - phisx(i1,i2) * etax(i1,i2) - phisy(i1,i2) * etay(i1,i2)
      ! dealiasing add_part
      detaphiz(i1,i2)      = geta2phiz2(i1,i2) + phiz2(i1,i2) - phisx(i1,i2) * phizMm2_etax_ext(i1,i2) &
           - phisy(i1,i2) * phizMm2_etay_ext(i1,i2) - phix_add_ext(i1,i2)*phiz_etax_ext(i1,i2) &
            - phiy_add_ext(i1,i2)*phiz_etay_ext(i1,i2) + phiz_add_ext(i1,i2)*phiz_ext(i1,i2)
      tmp1_add_ext(i1,i2)  = - 0.5_rp*(phix_add_ext(i1,i2)**2 + phiy_add_ext(i1,i2)**2 &
           +phiz_add_ext(i1,i2)**2 + 2.0_rp*phiz_ext(i1,i2)*phiz_add_ext(i1,i2) + 2.0_rp*phimx_ext(i1,i2)*phix_add_ext(i1,i2) &
           + 2.0_rp*phimy_ext(i1,i2)*phiy_add_ext(i1,i2)) !warning dealiasing phimx_ext... 3 multiplication...
      tmp2_add_ext(i1,i2)  = phix_add_ext(i1,i2)*etax(i1,i2)+ phiy_add_ext(i1,i2)*etay(i1,i2)
   end do
enddo
!
call dealias(2,phimx_ext, 'sin', 'cos')
call dealias(2,phimy_ext, 'cos', 'sin')
call dealias(2,phimx2_ext,'cos', 'cos')
call dealias(2,phimy2_ext,'cos', 'cos')
!
call filter_ext(detaphiz,detaphiz_d,    'cos', 'cos')
call filter_ext(phimxetax_ext,phimxetax,'cos', 'cos')
!
call filter_ext(phiz_etax_ext,phiz_etax,'sin', 'cos')
call filter_ext(phiz_etay_ext,phiz_etay,'cos', 'sin')
! dealiasing add_part
call filter_ext(tmp1_add_ext,tmp1_add,'cos', 'cos')
call filter_ext(tmp2_add_ext,tmp2_add,'cos', 'cos')
!
do i1=1,n1
   do i2=1,n2
      ! dealiasing add_part
      phimt(i1,i2) = dphis3(i1,i2) - (detaphiz_d(i1,i2)) 
   enddo
enddo
!
! Modif for partial dealiasing
call filter_ext(phimx_ext,phimx,  'sin', 'cos')
call filter_ext(phimy_ext,phimy,  'cos', 'sin')
!
call filter_ext(phimx2_ext,phimx2,'cos', 'cos')
call filter_ext(phimy2_ext,phimy2,'cos', 'cos')
call filter_ext(phimz2_ext,phimz2,'cos', 'cos')
!
! Transform etax and etay on the grid (1:n1,1:n2)
call filter_ext(etax,etax_tmp,'sin', 'cos')
call filter_ext(etay,etay_tmp,'cos', 'sin')
!
! Verification of FSBCs...
do i1=1,n1
   do i2=1,n2
      vitxref_FS(i1,i2)  = phimx(i1,i2)+phix_add(i1,i2)
      vityref_FS(i1,i2)  = phimy(i1,i2)+phiy_add(i1,i2)
      vitzref_FS(i1,i2)  = phimz(i1,i2)+phiz_add(i1,i2)
      phitref_FS(i1,i2)  = phimt(i1,i2)+phit_add(i1,i2)
      vitx2ref_FS(i1,i2) = phimx2(i1,i2)
      vitz2ref_FS(i1,i2) = phimz2(i1,i2)
      !
      ! CDSL OK (2D) - no dealiasing on add. part
      !
      !Pressref_FS(i1,i2) = -eta3(i1,i2) - 0.5_rp*(phimx2(i1,i2)+phix_add(i1,i2)**2 + phimy2(i1,i2)+phiy_add(i1,i2)**2 &
      !     + phimz2(i1,i2)+phiz_add(i1,i2)**2 + 2.0_rp*phimz(i1,i2)*phiz_add(i1,i2) + 2.0_rp*phimx(i1,i2)*phix_add(i1,i2) &
      !     + 2.0_rp*phimy(i1,i2)*phiy_add(i1,i2)) - (phimt(i1,i2)+phit_add(i1,i2))
      ! dealiasing add_part
      Press(i1,i2) = -eta3(i1,i2) - 0.5_rp*(phimx2(i1,i2) + phimy2(i1,i2) &
           + phimz2(i1,i2)) - (phimt(i1,i2)+phit_add(i1,i2)) + tmp1_add(i1,i2)
      !
      ! CCSL OK
      !
      !CCSL(i1,i2) = deta3(i1,i2) - phimxetax(i1,i2) + phix_add(i1,i2)*etax_tmp(i1,i2)+ phiy_add(i1,i2)*etay_tmp(i1,i2) - (phimz(i1,i2)+phiz_add(i1,i2)) 
      ! dealiasing add_part
      CCSL(i1,i2) = deta3(i1,i2) - phimxetax(i1,i2) - (phimz(i1,i2)+phiz_add(i1,i2)) + tmp2_add(i1,i2)
      ! Approximation (no dealiasing)
      Press_bis(i1,i2) = -eta3(i1,i2) - 0.5_rp*((vitxref_FS(i1,i2))**2+(vityref_FS(i1,i2))**2+(vitzref_FS(i1,i2))**2) &
            - (phitref_FS(i1,i2))
      CCSL_bis(i1,i2) = deta3(i1,i2) - (vitzref_FS(i1,i2)) + (vitxref_FS(i1,i2))*etax_tmp(i1,i2) &
            + (vityref_FS(i1,i2))*etay_tmp(i1,i2)
   enddo
enddo
!
! Look at error only outside absorption zones...
DO i1=1,n1
    DO i2=1,n2
        IF(ABS(nu(i1,i2)) >  tiny) THEN
            Press(i1,i2)       = 0.0_rp
            CCSL(i1,i2)        = 0.0_rp
            Press_bis(i1,i2)   = 0.0_rp
            CCSL_bis(i1,i2)    = 0.0_rp
        ENDIF
    ENDDO
ENDDO
!
!write(*,*) 'KFSBC =', MAXVAL(ABS(CCSL(1:n1,1:n2))), 'approx. (no dealiasing)=', MAXVAL(ABS(CCSL_bis(1:n1,1:n2)))
!write(*,*) 'DFSBC =', MAXVAL(ABS(Press(1:n1,1:n2))), 'approx. (no dealiasing)=', MAXVAL(ABS(Press_bis(1:n1,1:n2)))
!
end subroutine init_velocity_comp
!
! end init_velocity_comp ***************************************************
!
!
!
! start error_FS ****************************************************
!      
!     ============================================================
subroutine error_FS(modesspecx_l,modesspecy_l,modesspecz_l,modesspect_l,modesadd_l,modesaddt_l, modesFS_l, modesFSt_l, &
            vitxrefFS, vityrefFS, vitzrefFS, phitrefFS, error)
!     ============================================================
!
IMPLICIT NONE
!
!% INPUT VARIABLES
REAL(RP), DIMENSION(m1,m2), INTENT(IN)     :: modesspecx_l,modesspecy_l,modesspecz_l,modesspect_l,modesFS_l,modesFSt_l
REAL(RP), DIMENSION(m1,m2), INTENT(IN)     :: vitxrefFS, vityrefFS, vitzrefFS, phitrefFS
REAL(RP), DIMENSION(m3_add,m2), INTENT(IN) :: modesadd_l, modesaddt_l
!
REAL(RP), INTENT(OUT) :: error
!
! Local variables
REAL(RP), DIMENSION(m1,m2) :: vitx_l, vity_l, vitz_l, phitd_l, eta_l, deta_l, etax_l, etay_l
REAL(RP), DIMENSION(m1,m2) :: phiadd,phixadd,phizadd,phiyadd,phitadd
REAL(RP), DIMENSION(m1,m2) :: vitx,vity,vitz,Press,phitd,vitx2,vity2,vitz2
REAL(RP), DIMENSION(m1,m2) :: mat_coscosvit, mat_sincosvit, mat_cossinvit
!
REAL(RP) :: coeff
INTEGER  :: i1,i2,ii1,ii2
! Test FSBC
REAL(RP), DIMENSION(m1,m2) :: CCSL, Press_bis, CCSL_bis
!
! Reconstruct from modesFS and modesFSt
do i1=1,n1
   do i2=1,n2
      eta_l(i1,i2)  = modesFS_l(i1,i2) 
      deta_l(i1,i2) = modesFSt_l(i1,i2)
      etax_l(i1,i2) = -modesFS_l(i1,i2)*kx(i1)
      etay_l(i1,i2) = -modesFS_l(i1,i2)*ky(i2)
   end do
end do
!
eta_l  = Fourier_2_space(eta_l, 'cos','cos')
deta_l = Fourier_2_space(deta_l,'cos','cos')
etax_l = Fourier_2_space(etax_l,'sin','cos')
etay_l = Fourier_2_space(etay_l,'cos','sin')
!
! Add. part
CALL reconstruction_add_FS(modesadd_l,modesaddt_l,phiadd,phixadd,phizadd,&
     			phiyadd,phitadd,eta_l)
!
! Spec. part
do ii1=1,n1
   do ii2=1,n2
      do i2=1,n2
         do i1=1,n1
            mat_coscosvit(i1,i2)     = cos(kx(i1) * x(ii1,ii2))*cos(ky(i2) * y(ii1,ii2))
            mat_sincosvit(i1,i2)     = sin(kx(i1) * x(ii1,ii2))*cos(ky(i2) * y(ii1,ii2))
            mat_cossinvit(i1,i2)     = cos(kx(i1) * x(ii1,ii2))*sin(ky(i2) * y(ii1,ii2))
         end do
      end do
      vitx_l(1:n1,1:n2)  = 0.0_rp
      vity_l(1:n1,1:n2)  = 0.0_rp
      vitz_l(1:n1,1:n2)  = 0.0_rp
      phitd_l(1:n1,1:n2) = 0.0_rp
      do i1=1,n1
         do i2=1,n2
               IF(ABS(k(i1,i2)).GE.50.) THEN
                  coeff=EXP(k(i1,i2)*(eta_l(ii1,ii2)))
                  vitx_l(i1,i2)  = modesspecx_l(i1,i2)*coeff 
                  vity_l(i1,i2)  = modesspecy_l(i1,i2)*coeff
                  vitz_l(i1,i2)  = modesspecz_l(i1,i2)*coeff
                  phitd_l(i1,i2) = modesspect_l(i1,i2)*coeff
               ELSE
                  coeff=COSH(k(i1,i2)*(eta_l(ii1,ii2)+1.0_rp))/COSH(k(i1,i2))
                  vitx_l(i1,i2) = modesspecx_l(i1,i2)*coeff 
                  vity_l(i1,i2) = modesspecy_l(i1,i2)*coeff
                  phitd_l(i1,i2) = modesspect_l(i1,i2)*coeff    
                  IF (ABS(k(i1,i2)) > tiny) THEN
                     vitz_l(i1,i2) = modesspecz_l(i1,i2)*SINH(k(i1,i2)*(eta_l(ii1,ii2)+1.0_rp))/SINH(k(i1,i2))
                  ELSE
                     vitz_l(i1,i2) = modesspecz_l(i1,i2)! GD: change 07/2014... 0.0_rp !modesspecz_l(i1,i2) !modesspecz(j+jj) 
               ENDIF
            ENDIF
         end do
      end do
      !
      vitx(ii1,ii2)  = 0.0_rp
      vity(ii1,ii2)  = 0.0_rp
      vitz(ii1,ii2)  = 0.0_rp
      phitd(ii1,ii2) = 0.0_rp
      !
      DO i2=1,n2
         vitx(ii1,ii2) = vitx(ii1,ii2)  + DOT_PRODUCT(vitx_l(1:n1,i2) ,mat_sincosvit(1:n1,i2))
         vity(ii1,ii2) = vity(ii1,ii2)  + DOT_PRODUCT(vity_l(1:n1,i2) ,mat_cossinvit(1:n1,i2))
         vitz(ii1,ii2) = vitz(ii1,ii2)  + DOT_PRODUCT(vitz_l(1:n1,i2) ,mat_coscosvit(1:n1,i2))
         phitd(ii1,ii2)= phitd(ii1,ii2) + DOT_PRODUCT(phitd_l(1:n1,i2),mat_coscosvit(1:n1,i2))
      ENDDO      
      vitx2(ii1,ii2)= vitx(ii1,ii2)**2
      vity2(ii1,ii2)= vity(ii1,ii2)**2
      vitz2(ii1,ii2)= vitz(ii1,ii2)**2
      !
      Press(ii1,ii2) = -(eta_l(ii1,ii2))-0.5_rp*((vitx(ii1,ii2)+phixadd(ii1,ii2))**2+&
		  (vity(ii1,ii2)+phiyadd(ii1,ii2))**2+(vitz(ii1,ii2)+phizadd(ii1,ii2))**2)-(phitd(ii1,ii2)+phitadd(ii1,ii2))
      CCSL(ii1,ii2) = deta_l(ii1,ii2)- (vitz(ii1,ii2)+phizadd(ii1,ii2))+(vitx(ii1,ii2)+phixadd(ii1,ii2))*etax_l(ii1,ii2)&
            +(vity(ii1,ii2)+phiyadd(ii1,ii2))*etay_l(ii1,ii2)
      !
      vitx(ii1,ii2) = vitx(ii1,ii2) + phixadd(ii1,ii2)
      vity(ii1,ii2) = vity(ii1,ii2) + phiyadd(ii1,ii2)
      vitz(ii1,ii2) = vitz(ii1,ii2) + phizadd(ii1,ii2)
      phitd(ii1,ii2)= phitd(ii1,ii2)+ phitadd(ii1,ii2)
   enddo
enddo
!
! Check input values
DO ii1=1,n1
    DO ii2=1,n2
        Press_bis(ii1,ii2) = -(eta_l(ii1,ii2))-0.5_rp*((vitxrefFS(ii1,ii2))**2+&
              (vityrefFS(ii1,ii2))**2+(vitzrefFS(ii1,ii2))**2)-(phitrefFS(ii1,ii2))
        CCSL_bis(ii1,ii2) = deta_l(ii1,ii2)- (vitzrefFS(ii1,ii2))+(vitxrefFS(ii1,ii2))*etax_l(ii1,ii2) &
                +(vityrefFS(ii1,ii2))*etay_l(ii1,ii2)
    ENDDO
ENDDO
!
! Look at error only outside absorption zones...
DO ii1=1,n1
    DO ii2=1,n2
        IF(ABS(nu(ii1,ii2)) > tiny) THEN
            Press(ii1,ii2)     = 0.0_rp
            CCSL(ii1,ii2)      = 0.0_rp
            Press_bis(ii1,ii2) = 0.0_rp
            CCSL_bis(ii1,ii2)  = 0.0_rp
        ENDIF
    ENDDO
ENDDO
!
! Evaluate error on free surface boundary conditions (initial with aliasing errors and reconstructed)
write(*,*) 'KFSBC/DFSBC max. initial = ', MAX(MAXVAL(ABS(CCSL_bis(1:n1,1:n2))),MAXVAL(ABS(Press_bis(1:n1,1:n2)))), &
	', reconstructed = ', MAX(MAXVAL(ABS(CCSL(1:n1,1:n2))),MAXVAL(ABS(Press(1:n1,1:n2))))
!
! Compute max error on velocities and phit
error = MAXVAL(ABS(vitx-vitxrefFS))/MAXVAL(ABS(vitxrefFS))
error = MAX(error,MAXVAL(ABS(vitz-vitzrefFS))/MAXVAL(ABS(vitzrefFS)))
error = MAX(error,MAXVAL(ABS(phitd-phitrefFS))/MAXVAL(ABS(phitrefFS)))
IF (n2 /= 1) THEN
	error = MAX(error,MAXVAL(ABS(vity-vityrefFS))/MAXVAL(ABS(vityrefFS)))
ENDIF
!
end subroutine error_FS
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Reconstruction additional part on Free Surface
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE reconstruction_add_FS(modesadd_l,modesaddt_l,phiadd,phixadd,phizadd,phiyadd,phitadd,eta_l)
!
IMPLICIT NONE
!
DOUBLE PRECISION, DIMENSION(m3_add,m2) :: modesadd_l,modesaddt_l
DOUBLE PRECISION, DIMENSION(m1,m2)     :: eta_l
DOUBLE PRECISION, DIMENSION(m1,m2)     :: phiadd,phixadd,phizadd,phiyadd,phitadd
!
REAL(RP) :: kx_add_z, coskx_add_eta, sinkx_add_eta, cosky, sinky, coeff1, coeff2, coeff3
!
!% Local variables
INTEGER :: ii1,ii2,i2,i3
!
if(n2.EQ.1) then
   do ii1=1,n1
		phiadd(ii1,1)   = 0.0d0
		phitadd(ii1,1)  = 0.0d0
		phixadd(ii1,1)  = 0.0d0
		phiyadd(ii1,1)  = 0.0d0
		phizadd(ii1,1)  = 0.0d0
		do i3=1,n3_add
		   kx_add_z  = kx_add(i3) * (eta_l(ii1,1)+1.d0) !kx_add(i3) * (zvect(i1)+1.d0)
		   phiadd(ii1,1)   = phiadd(ii1,1)   +  modesadd_l(i3,1) * cos(kx_add_z) * csh_add_x(ii1,1,i3)
		   phitadd(ii1,1)  = phitadd(ii1,1)  + modesaddt_l(i3,1) * cos(kx_add_z) * csh_add_x(ii1,1,i3)
		   phixadd(ii1,1)  = phixadd(ii1,1)  -  modesadd_l(i3,1) * cos(kx_add_z) * k_add_sh_add_x(ii1,1,i3)
		   phizadd(ii1,1)  = phizadd(ii1,1)  -  modesadd_l(i3,1) * sin(kx_add_z) * kx_add_csh_add_x(ii1,1,i3)
		end do
   end do
else
    do ii1=1,n1
        do ii2=1,n2
            phiadd(ii1,ii2)   = 0.0d0
            phitadd(ii1,ii2)  = 0.0d0
            phixadd(ii1,ii2)  = 0.0d0
            phiyadd(ii1,ii2)  = 0.0d0
            phizadd(ii1,ii2)  = 0.0d0
            do i3=1,n3_add
                kx_add_z  = kx_add(i3) * (eta_l(ii1,1)+1.d0) !kx_add(i3) * (zvect(i1)+1.d0)
                coskx_add_eta = cos(kx_add_z)
                sinkx_add_eta = sin(kx_add_z)
                do i2=1,n2
                    cosky =  cos(ky(i2)*y(ii1,ii2))
                    sinky =  sin(ky(i2)*y(ii1,ii2))
                    coeff1 = coskx_add_eta * csh_add_x(ii1,i2,i3)* cosky
                    coeff2 = sinkx_add_eta * kx_add_csh_add_x(ii1,i2,i3)
                    coeff3 = k_add_sh_add_x(ii1,i2,i3)
                    phiadd(ii1,ii2)     =  phiadd(ii1,ii2)   +  modesadd(i3,i2) * coeff1 
                    phitadd(ii1,ii2)    =  phitadd(ii1,ii2)  +  modesaddt(i3,i2)* coeff1 
                    phixadd(ii1,ii2)    =  phixadd(ii1,ii2)  -  modesadd(i3,i2) * coskx_add_eta * coeff3 * cosky 
                    phiyadd(ii1,ii2)    =  phiyadd(ii1,ii2)  -  modesadd(i3,i2) * coskx_add_eta * kycsh_add_x(ii1,i2,i3) * sinky
                    phizadd(ii1,ii2)    =  phizadd(ii1,ii2)  -  modesadd(i3,i2) * coeff2 * cosky
                enddo
            enddo
        enddo
    enddo
endif
!
END SUBROUTINE reconstruction_add_FS
!
! end reconstruction_add_FS
!
!
!
! start velpress
!
!     ============================================================
subroutine velpress(modesx,modesy,modesz,modest,modeseta,modes_add,modes_addt,n_press,zpress,ypress, &
           vitx_press,vity_press,vitz_press,Press_press,eta_press)
!     ============================================================
!
IMPLICIT NONE
!
REAL(RP), DIMENSION(m1,m2), INTENT(IN)      :: modesx,modesy,modesz,modest,modeseta
REAL(RP), DIMENSION(m3_add,m2), INTENT(IN)  :: modes_add,modes_addt
REAL(RP), DIMENSION(maxprobes), INTENT(OUT) :: vitx_press,vity_press,vitz_press,Press_press,eta_press
!
INTEGER  :: n_press
REAL(RP), DIMENSION(maxprobes) :: ypress,zpress
!
REAL(RP), DIMENSION(maxprobes) :: phiadd_press, phixadd_press, phizadd_press, phiyadd_press, phitadd_press
REAL(RP), DIMENSION(maxprobes) :: vitx2_press, vity2_press, vitz2_press, phitd_press
!
! Local variables
REAL(RP), DIMENSION(m1,m2) :: vitx_l, vity_l, vitz_l, phitd_l
REAL(RP) :: coeff
INTEGER  :: i1, i2, ii
!
! Additional part
CALL phiadd_1st_press(modes_add,modes_addt,phiadd_press,phixadd_press,phizadd_press,phiyadd_press,phitadd_press,&
						zpress, ypress, n_press)
!
! Spectral part
do ii=1,n_press
   vitx_l(1:n1,1:n2)   = 0.0_rp
   vitz_l(1:n1,1:n2)   = 0.0_rp
   phitd_l(1:n1,1:n2)  = 0.0_rp
   do i1=1,n1 
      do i2=1,n2
         IF(ABS(k(i1,i2)*(zpress(ii)+1.0_rp)).GE.50.) THEN
            IF(ABS(k(i1,i2)).GE.50.) THEN
               coeff=EXP(k(i1,i2)*(zpress(ii)+1.0_rp-1.0_rp))
               vitx_l(i1,i2)  = modesx(i1,i2)*coeff 
               vity_l(i1,i2)  = modesy(i1,i2)*coeff 
               vitz_l(i1,i2)  = modesz(i1,i2)*coeff 
               phitd_l(i1,i2) = modest(i1,i2)*coeff	
            ELSE
               coeff=EXP(k(i1,i2)*(zpress(ii)+1.0_rp))/(2.0_rp*COSH(k(i1,i2)))
               vitx_l(i1,i2)  = modesx(i1,i2)*coeff 
               vity_l(i1,i2)  = modesy(i1,i2)*coeff
               vitz_l(i1,i2)  = modesz(i1,i2)*EXP(k(i1,i2)*(zpress(ii)+1.0_rp))/(2.0_rp*SINH(k(i1,i2)))
               phitd_l(i1,i2) = modest(i1,i2)*coeff
            ENDIF
         ELSE
            IF(ABS(k(i1,i2)).GE.50.) THEN
               coeff=COSH(k(i1,i2)*(zpress(ii)+1.0_rp))*EXP(-k(i1,i2))*2.0_rp
               vitx_l(i1,i2)  = modesx(i1,i2)*coeff 
               vity_l(i1,i2)  = modesy(i1,i2)*coeff
               vitz_l(i1,i2)  = modesz(i1,i2)*SINH(k(i1,i2)*(zpress(ii)+1.0_rp))*EXP(-k(i1,i2))*2.0_rp
               phitd_l(i1,i2) = modest(i1,i2)*coeff
            ELSE
               coeff=COSH(k(i1,i2)*(zpress(ii)+1.0_rp))/COSH(k(i1,i2))
               vitx_l(i1,i2)=modesx(i1,i2)*coeff 
               vity_l(i1,i2)=modesy(i1,i2)*coeff
               phitd_l(i1,i2)=modest(i1,i2)*coeff
               !
               IF (ABS(k(i1,i2)) > tiny) THEN
                  vitz_l(i1,i2) = modesz(i1,i2)*SINH(k(i1,i2)*(zpress(ii)+1.0_rp))/SINH(k(i1,i2)) 
               ELSE
                  vitz_l(i1,i2) = 0.0_rp 
               ENDIF
            ENDIF
         ENDIF
      end do
   end do
   !
   vitx_press(ii)  = 0.0_rp
   vity_press(ii)  = 0.0_rp
   vitz_press(ii)  = 0.0_rp
   phitd_press(ii) = 0.0_rp
   eta_press(ii)   = 0.0_rp
   !
   vitx2_press(ii) = 0.0_rp
   vity2_press(ii) = 0.0_rp
   vitz2_press(ii) = 0.0_rp
   Press_press(ii) = 0.0_rp 
   !
   DO i2=1,n2
      vitx_press(ii)  = vitx_press(ii)   + DOT_PRODUCT(vitx_l(1:n1,i2),mat_sincos(1:n1,i2,ii))
      vity_press(ii)  = vity_press(ii)   + DOT_PRODUCT(vity_l(1:n1,i2),mat_cossin(1:n1,i2,ii))
      vitz_press(ii)  = vitz_press(ii)   + DOT_PRODUCT(vitz_l(1:n1,i2),mat_coscos(1:n1,i2,ii))
      phitd_press(ii) = phitd_press(ii)  + DOT_PRODUCT(phitd_l(1:n1,i2),mat_coscos(1:n1,i2,ii))
      eta_press(ii)   = eta_press(ii)    + DOT_PRODUCT(modeseta(1:n1,i2),mat_coscos(1:n1,i2,ii))
   ENDDO
   !
   vitx2_press(ii)= vitx_press(ii)**2
   vity2_press(ii)= vity_press(ii)**2
   vitz2_press(ii)= vitz_press(ii)**2
   !
   ! Computation of pressure only under the free surface
   if(zpress(ii).GT.eta_press(ii) - (voltot-0.5_rp*alpha*(1-d_hinge)*pos1st(n3,1) / xlen  &
        - eta(1,1)*alpha*pos1st(n3,1) / xlen)) then
      Press_press(ii) = 0.0_rp
   else
      Press_press(ii) = 0.0_rp-0.0_rp*(zpress(ii))-0.5_rp*(vitx2_press(ii)+(phixadd_press(ii))**2 &
           +vity2_press(ii)+(phiyadd_press(ii))**2+vitz2_press(ii)+(phizadd_press(ii))**2  &
           + 2.0_rp*vitx_press(ii)*phixadd_press(ii) + 2.0_rp*vity_press(ii)*phiyadd_press(ii) &
           + 2.0_rp*vitz_press(ii)*phizadd_press(ii))-(phitd_press(ii)+phitadd_press(ii))
   endif
   !
   vitx_press(ii) = vitx_press(ii) + phixadd_press(ii)
   vity_press(ii) = vity_press(ii) + phiyadd_press(ii)
   vitz_press(ii) = vitz_press(ii) + phizadd_press(ii)
enddo
!
end subroutine velpress
!
! end velpress
!
!
!
! start phiadd_1st_press ****************************************************
     
!     ==================================================================
      subroutine phiadd_1st_press(a1st_add_l,da1st_add_l,phiadd,phixadd,phizadd,phiyadd,phitadd, &
           zpress_l, ypress_l, n_press)
!     ==================================================================
!      
IMPLICIT NONE
!
!% INPUT VARIABLES
!
!	definition of horizontal and vertical position of pressure probes
REAL(RP) :: ypress_l(*), zpress_l(*)
INTEGER  :: n_press
!
!	first-order additional time modal amplitudes
REAL(RP),DIMENSION(m3_add,m2) :: a1st_add_l,da1st_add_l
!
! Quantities needed for velocity and pressure (at probe position)
REAL(RP), DIMENSION(maxprobes) :: phiadd,phixadd,phizadd,phiyadd,phitadd
!
!% LOCAL VARIABLES
REAL(RP) :: ti, tf
INTEGER  :: j1, i2, i3
REAL(RP) :: kx_add_zpress, coskx_add_zpress, sinkx_add_zpress, cosky, sinky
!
!	CPU times inlet
if (iCPUtime.eq.1) then
   print*,'entering subroutine phiadd_1st_press'
   call CPU_TIME(ti)
endif
!
!	re-construction on the position of pressure probes
if(n2.EQ.1) then
   do j1 = 1,n_press
      phiadd(j1)   = 0.0_rp
      phitadd(j1)  = 0.0_rp
      phixadd(j1)  = 0.0_rp
      phiyadd(j1)  = 0.0_rp
      phizadd(j1)  = 0.0_rp
      do i3=1,n3_add
         kx_add_zpress    = kx_add(i3) * (zpress_l(j1)+1.d0)
         coskx_add_zpress = cos(kx_add_zpress)
         sinkx_add_zpress = sin(kx_add_zpress)
         phiadd(j1)    = phiadd(j1)   +  a1st_add_l(i3,1) * coskx_add_zpress * csh_add_x_probe(i3,1,j1)
         phitadd(j1)   = phitadd(j1)  + da1st_add_l(i3,1) * coskx_add_zpress * csh_add_x_probe(i3,1,j1)
         phixadd(j1)   = phixadd(j1)  -  a1st_add_l(i3,1) * coskx_add_zpress * k_add_sh_add_x_probe(i3,1,j1)
         phiyadd(j1)   = 0.0_rp
         phizadd(j1)   = phizadd(j1)  -  a1st_add_l(i3,1) * sinkx_add_zpress * kx_add_csh_add_x_probe(i3,1,j1)
      end do
   end do
else
   do j1 = 1,n_press
      phiadd(j1)   = 0.0_rp
      phitadd(j1)  = 0.0_rp
      phixadd(j1)  = 0.0_rp
      phiyadd(j1)  = 0.0_rp
      phizadd(j1)  = 0.0_rp
      do i2=1,n2
         cosky = COS(ky(i2)*ypress_l(j1))
         sinky = SIN(ky(i2)*ypress_l(j1))
         do i3=1,n3_add
            kx_add_zpress  = kx_add(i3) * (zpress_l(j1)+1.d0)
            coskx_add_zpress = cos(kx_add_zpress)
            sinkx_add_zpress = sin(kx_add_zpress)
            phiadd(j1)   = phiadd(j1)   +  a1st_add_l(i3,i2) * coskx_add_zpress * csh_add_x_probe(i3,i2,j1)*cosky
            phitadd(j1)  = phitadd(j1)  + da1st_add_l(i3,i2) * coskx_add_zpress * csh_add_x_probe(i3,i2,j1)*cosky
            phixadd(j1)  = phixadd(j1)  -  a1st_add_l(i3,i2) * coskx_add_zpress * k_add_sh_add_x_probe(i3,i2,j1)*cosky
            phiyadd(j1)  = phiyadd(j1)  -  a1st_add_l(i3,i2) * coskx_add_zpress * kycsh_add_x_probe(i3,i2,j1)*sinky
            phizadd(j1)  = phizadd(j1)  -  a1st_add_l(i3,i2) * sinkx_add_zpress * kx_add_csh_add_x_probe(i3,i2,j1)*cosky
         end do
      end do
   end do
endif
!
!	CPU times outlet
if (iCPUtime.eq.1) then
   call CPU_TIME(tf)
   write(*,910)'quitting subroutine phiadd_1st_press, total CPU time: ',tf-ti,'s'
endif
!
910 format(a,1ES11.4,a)
!
end subroutine phiadd_1st_press
!
! end phiadd_1st_press ******************************************************
!
!
! start HOSvel2_direct *************************************************
!
!     ======================================
      subroutine HOSvel2_direct(eta_l,phis_l,a1st_add_l,da1st_add_l,a2nd_add_l,da2nd_add_l,a3rd_add_l,da3rd_add_l, &
           a1st_loc,da1st_loc, eta1st_loc, deta1st_loc, a2nd_loc, da2nd_loc, eta2nd_loc, deta2nd_loc, &
           G2_loc, dG2_loc, G1_loc, dG1_loc, da_add_ramp_loc, h_loc)
!     ======================================
!
IMPLICIT NONE
!
!% INPUT VARIABLES
REAL(RP),DIMENSION(m1,m2)     :: eta_l,phis_l,a1st_loc,da1st_loc,eta1st_loc,deta1st_loc,a2nd_loc,da2nd_loc,eta2nd_loc
REAL(RP),DIMENSION(m3_add,m2) :: a1st_add_l,da1st_add_l,a2nd_add_l,da2nd_add_l,a3rd_add_l,da3rd_add_l,da_add_ramp_loc
REAL(RP),DIMENSION(m1,m2)     :: deta2nd_loc,G2_loc,dG2_loc,G1_loc,dG1_loc
REAL(RP) :: h_loc
!
!% LOCAL VARIABLES
REAL(RP),DIMENSION(m1,m2)          :: phimx,phimy,phimz,phimt
REAL(RP),DIMENSION(md1,md2,maxHOS) ::  etapmext2
!	  
INTEGER  :: i1, i2
REAL(RP) :: ti,tf
!
!	CPU times inlet
!
if(iCPUtime.eq.1) then
   print*,'entering subroutine HOSvel2_direct'
   call CPU_TIME(ti)
endif
!
! Initialisation of velocity computations
CALL init_velocity_comp(eta_l,phis_l,da1st_add_l,a1st_add_l,da2nd_add_l,a2nd_add_l,da3rd_add_l,a3rd_add_l, &
     	a1st_loc,da1st_loc,eta1st_loc,deta1st_loc,a2nd_loc,da2nd_loc,eta2nd_loc,deta2nd_loc,&
     	G2_loc,dG2_loc,G1_loc,dG1_loc,da_add_ramp_loc,h_loc,phimx,phimy,phimz,phimt)
!
do i1 = 1,nd1
   do i2 = 1,nd2
      etapmext2(i1,i2,1:mHOS)=etapmext(i1,i2,1:mHOS)
   end do
enddo
!
!
IF (n2.eq.1) THEN
CALL Velocity_modes_direct(eta_l(:,1),phimx(:,1),phimy(:,1),phimz(:,1),phimt(:,1), &
	modesspecx(:,1),modesspecy(:,1),modesspecz(:,1),modesspect(:,1))
ELSE
	print*, 'Direct method should not be used in 3D due to RAM memory use'
	stop
ENDIF
!
!	CPU times outlet
if(iCPUtime.eq.1) then
   call CPU_TIME(tf)
   print*,'quitting subroutine HOSvel2_direct', tf-ti
endif
!
end subroutine HOSvel2_direct
!
! end HOSvel2_direct ***************************************************
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE Velocity_modes_direct(eta,phix,phiy,phiz,phit,&
	modesphix,modesphiy,modesphiz,modesphit)
!
USE variables, ONLY : m1,m2,k,kx,ky,x,y
!
IMPLICIT NONE
!
REAL(RP), DIMENSION(m1,1), INTENT(IN)  :: eta,phix,phiy,phiz,phit
REAL(RP), DIMENSION(m1,1), INTENT(OUT) :: modesphix, modesphiy, modesphiz, modesphit
!
! Local variable 
REAL(RP), DIMENSION(m1*1,m1*1) :: M
!
! Check m2=1
IF (m2.NE.1) THEN
	print*,'direct inversion only working (RAM memory use) in 2D'
	stop
ENDIF
!
! Use velocities on FS
!
! Create Matrix for dphi_dx and solve AX=B
CALL CreateMatrix_sincos_cosh(kx,ky(1),k(:,1),x(:,1),y(:,1),eta,M)
CALL solve_system(M,phix,modesphix)
!
! Create Matrix for dphi_dz and solve AX=B
CALL CreateMatrix_coscos_sinh(kx,ky(1),k(:,1),x(:,1),y(:,1),eta,M)
CALL solve_system(M,phiz,modesphiz)
!
!modesphiz(1,1) = 0.0_rp
!
IF(n2.EQ.1)THEN
	modesphiy=0.0_rp
ELSE
	! Create Matrix for dphi_dy and solve AX=B
	CALL CreateMatrix_cossin_cosh(kx,ky(1),k(:,1),x(:,1),y(:,1),eta,M)
	CALL solve_system(M,phiy,modesphiy)
ENDIF
!
! Create the matrix for dphisdt and solve AX=B
CALL CreateMatrix_coscos_cosh(kx,ky(1),k(:,1),x(:,1),y(:,1),eta,M)
CALL Solve_system(M,phit,modesphit)
!
END SUBROUTINE Velocity_modes_direct
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE Solve_system(A,B,X)
!
IMPLICIT NONE
!
REAL(RP), DIMENSION(m1*1,m1*1), INTENT(IN) :: A
REAL(RP), DIMENSION(m1*1,m1*1) :: A_save
REAL(RP), DIMENSION(m1,1), INTENT(IN) :: B
REAL(RP), DIMENSION(m1*1) :: work,work2
REAL(RP), DIMENSION(m1,1), INTENT(OUT) :: X
! 
INTEGER, DIMENSION(m1*1)  :: ipiv
!
INTEGER :: i1,i2,index1,info
!
! Solve AX=B
DO i2=1,1
	DO i1=1,n1
		index1 = i1+n1*(i2-1)
		work(index1) = B(i1,i2)
		work2(index1)= work(index1)
	ENDDO
ENDDO
!
A_save(1:n1*1,1:n1*1)=A(1:n1*1,1:n1*1)
CALL dgesv(n1*1,1,A_save(1:n1*1,1:n1*1),n1*1,ipiv(1:n1*1), work(1:n1*1), n1*1, info)
!
!CALL ZGETRF(n1*n2,n1*n2,AF,n1*n2,IPIV,INFO)
!CALL zgesvx('E','N',n1*n2,1,A,n1*n2,A,n1*n2,ipiv,EQUED,R,C,work,n1*n2,work,n1*n2,RCOND,FERR,BERR
!
IF(info/=0)THEN
  PRINT *, 'Problems with inversion of the matrix, info=',info
  STOP
END IF
!!
!! Evaluate error of inversion
!!
!DO index1=1,n1*n2
!	error(index1) = SUM(A(index1,1:n1*n2)*work(1:n1*n2))-work2(index1)
!ENDDO
!!
!print*,'error max=', MAXVAL(ABS(error)),MAXVAL(ABS(work2))
!!
DO i1=1,n1
	DO i2=1,1
		index1=i1+n1*(i2-1)
		X(i1,i2) = work(index1)
	ENDDO
ENDDO
!
END SUBROUTINE solve_system
!
!
!
SUBROUTINE CreateMatrix_coscos_cosh(kx,ky,k,x,y,eta,M)
! By G. Ducrozet.
USE variables, ONLY : m1
IMPLICIT NONE
!
INTEGER :: i1,i2,ii1,ii2,index1,index2
REAL(RP), DIMENSION(m1,1), INTENT(IN) :: x, y, k
REAL(RP), DIMENSION(m1), INTENT(IN)   :: kx
REAL(RP), DIMENSION(1), INTENT(IN)    :: ky
REAL(RP), DIMENSION(m1,1), INTENT(IN) :: eta
REAL(RP), DIMENSION(m1*1,m1*1), INTENT(OUT) :: M
!
M = 0.0_rp
!
DO i2=1,1
	DO i1=1,n1
		index1=i1+n1*(i2-1)
		DO ii2=1,1
			DO ii1=1,n1
				index2=ii1+n1*(ii2-1)
				IF(k(ii1,ii2).LT.50.0_rp) THEN
					M(index1,index2) = COS(kx(ii1)*x(i1,i2))*COS(ky(ii2)*y(i1,i2))*COSH(k(ii1,ii2)*(eta(i1,i2)+1.0_rp))/COSH(k(ii1,ii2))
				ELSE
					M(index1,index2) = COS(kx(ii1)*x(i1,i2))*COS(ky(ii2)*y(i1,i2))*EXP(k(ii1,ii2)*eta(i1,i2))
				ENDIF
			ENDDO
		ENDDO
	ENDDO
ENDDO
!
END SUBROUTINE CreateMatrix_coscos_cosh
!
!
!
SUBROUTINE CreateMatrix_sincos_cosh(kx,ky,k,x,y,eta,M)
!
! By G. Ducrozet.
USE variables, ONLY : m1
!
IMPLICIT NONE
!
INTEGER :: i1,i2,ii1,ii2,index1,index2
REAL(RP), DIMENSION(m1,1), INTENT(IN) :: x, y, k
REAL(RP), DIMENSION(m1), INTENT(IN)    :: kx
REAL(RP), DIMENSION(1), INTENT(IN)    :: ky
REAL(RP), DIMENSION(m1,1), INTENT(IN) :: eta
REAL(RP), DIMENSION(m1*1,m1*1), INTENT(OUT) :: M
!
M = 0.0_rp
!
DO i2=1,1
	DO i1=1,n1
		index1=i1+n1*(i2-1)
		DO ii2=1,1
			ii1=1
			index2=ii1+n1*(ii2-1)
			M(index1,index2) = 1.0_rp !to enable the inversion of the matrix... 
			!(it has no influence on result since to reconstruct this mode is multiplied by zero afterward)
			DO ii1=2,n1
				index2=ii1+n1*(ii2-1)
				IF(k(ii1,ii2).LT.50.0_rp) THEN
					M(index1,index2) = SIN(kx(ii1)*x(i1,i2))*COS(ky(ii2)*y(i1,i2))*COSH(k(ii1,ii2)*(eta(i1,i2)+1.0_rp))/COSH(k(ii1,ii2))
				ELSE
					M(index1,index2) = SIN(kx(ii1)*x(i1,i2))*COS(ky(ii2)*y(i1,i2))*EXP(k(ii1,ii2)*eta(i1,i2))
				ENDIF
			ENDDO
		ENDDO
	ENDDO
ENDDO
!
END SUBROUTINE CreateMatrix_sincos_cosh
!
!
!
SUBROUTINE CreateMatrix_cossin_cosh(kx,ky,k,x,y,eta,M)
!
! By G. Ducrozet.
USE variables, ONLY : m1
!
IMPLICIT NONE
!
INTEGER :: i1,i2,ii1,ii2,index1,index2
REAL(RP), DIMENSION(m1,1), INTENT(IN)  :: x, y, k
REAL(RP), DIMENSION(m1), INTENT(IN)    :: kx
REAL(RP), DIMENSION(1), INTENT(IN)     :: ky
REAL(RP), DIMENSION(m1,1), INTENT(IN)  :: eta
REAL(RP), DIMENSION(m1*1,m1*1), INTENT(OUT) :: M
!
M = 0.0_rp
!
DO i2=1,1
	DO i1=1,n1
		index1=i1+n1*(i2-1)
		ii2=1
		DO ii1=1,n1
			index2=ii1+n1*(ii2-1)
			M(index1,index2) = 1.0_rp !to enable the inversion of the matrix... 
			!(it has no influence on result since to reconstruct this mode is multiplied by zero afterward)
			IF(k(ii1,ii2).LT.50.0_rp) THEN
				M(index1,index2) = COS(kx(ii1)*x(i1,i2))*SIN(ky(ii2)*y(i1,i2))*COSH(k(ii1,ii2)*(eta(i1,i2)+1.0_rp))/COSH(k(ii1,ii2))
			ELSE
				M(index1,index2) = COS(kx(ii1)*x(i1,i2))*SIN(ky(ii2)*y(i1,i2))*EXP(k(ii1,ii2)*eta(i1,i2))
			ENDIF
		ENDDO
		DO ii2=2,1
			DO ii1=1,n1
				index2=ii1+n1*(ii2-1)
				IF(k(ii1,ii2).LT.50.0_rp) THEN
					M(index1,index2) = COS(kx(ii1)*x(i1,i2))*SIN(ky(ii2)*y(i1,i2))*COSH(k(ii1,ii2)*(eta(i1,i2)+1.0_rp))/COSH(k(ii1,ii2))
				ELSE
					M(index1,index2) = COS(kx(ii1)*x(i1,i2))*SIN(ky(ii2)*y(i1,i2))*EXP(k(ii1,ii2)*eta(i1,i2))
				ENDIF
			ENDDO
		ENDDO
	ENDDO
ENDDO
!
END SUBROUTINE CreateMatrix_cossin_cosh
!
!
!
SUBROUTINE CreateMatrix_coscos_sinh(kx,ky,k,x,y,eta,M)
!
! By G. Ducrozet.
USE variables, ONLY : m1
!
IMPLICIT NONE
!
INTEGER :: i1,i2,ii1,ii2,index1,index2
REAL(RP), DIMENSION(m1,1), INTENT(IN) :: x, y, k
REAL(RP), DIMENSION(m1), INTENT(IN)   :: kx
REAL(RP), DIMENSION(1), INTENT(IN)    :: ky
REAL(RP), DIMENSION(m1,1), INTENT(IN) :: eta
REAL(RP), DIMENSION(m1*1,m1*1), INTENT(OUT) :: M
!
M = 0.0_rp
!
DO i2=1,1
	DO i1=1,n1
		index1=i1+n1*(i2-1)
		DO ii2=1,1
			DO ii1=1,n1
				index2=ii1+n1*(ii2-1)
				IF(k(ii1,ii2).LT.50.0_rp) THEN
					M(index1,index2) = COS(kx(ii1)*x(i1,i2))*COS(ky(ii2)*y(i1,i2))*SINH(k(ii1,ii2)*(eta(i1,i2)+1.0_rp))/SINH(k(ii1,ii2))
				ELSE
					M(index1,index2) = COS(kx(ii1)*x(i1,i2))*COS(ky(ii2)*y(i1,i2))*EXP(k(ii1,ii2)*eta(i1,i2))
				ENDIF
			ENDDO
		ENDDO
	ENDDO
ENDDO
!
! Specific treatment of k=0
M(:,1) = 1.0_rp
!
END SUBROUTINE CreateMatrix_coscos_sinh
!
!
!
SUBROUTINE InvertMatrix(M,rank)
!
! Inverse real matrix with LAPACK
IMPLICIT NONE
!
INTEGER :: rank, lwork, info
REAL(RP), DIMENSION(rank,rank) :: M
REAL(RP), DIMENSION(5*rank**2) :: work
INTEGER, DIMENSION(rank)       :: ipiv
!
rank = SIZE(M,1)
CALL dgetrf(rank,rank,M,rank,ipiv,info)
IF(info/=0)THEN
   PRINT *, 'Problems with L-U of the matrix',info
   STOP
END IF
lwork=rank
CALL dgetri(rank,M,rank,ipiv,work,lwork,info)
IF(info/=0)THEN
  PRINT *, 'Problems with inversion of the matrix',info
  STOP
END IF
!
END SUBROUTINE InvertMatrix
!
! end module
!
END MODULE velocities
