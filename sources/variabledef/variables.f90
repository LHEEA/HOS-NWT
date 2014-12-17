MODULE variables
!
! This module defines all global variables used in HOS-NWT
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
USE wave_def
USE common_vars
!
IMPLICIT NONE
!
! Adaptative Time Step Runge Kutta scheme
INTEGER, PARAMETER             :: RK_s = 5
INTEGER                        :: RK_p
INTEGER, DIMENSION(RK_s)       :: RK_substep
REAL(RP), DIMENSION(RK_s,RK_s) :: RK_A ! Butcher array
REAL(RP), DIMENSION(RK_s)      :: RK_b ! Butcher array: increment factor for slope
REAL(RP), DIMENSION(RK_s)      :: RK_c ! Butcher array: increment factor for abscissa
REAL(RP), DIMENSION(RK_s)      :: RK_e ! Butcher array: weighting factors for the error estimate
REAL(RP)                       :: toler, f_out, time_cur, time_next
!
! Variable change
REAL(RP), DIMENSION(m1,m2) :: omega
!
! Free Surface elevation and its derivatives
REAL(RP), DIMENSION(m1,m2)          :: eta, deta
REAL(RP), DIMENSION(md1,md2)        :: etax, etay
REAL(RP), DIMENSION(md1,md2,maxHOS) :: etapmext
!
! Free Surface potential and its derivatives
REAL(RP), DIMENSION(m1,m2)          :: phis
REAL(RP), DIMENSION(md1,md2)        :: phisx, phisy
!
! Ordered components of phis derivatives for de-aliasing
REAL(RP), DIMENSION(m1,m2)          :: phiz, W1
REAL(RP), DIMENSION(md1,md2)        :: phiz2, phizMm1, phizMm2, phiz2Mm2
REAL(RP), DIMENSION(md1,md2)        :: gradeta2, geta2phiz2, geta2phiz
!
! Modes of the different orders of phi
REAL(RP), DIMENSION(m1,m2,maxHOS)   :: aHOS
REAL(RP), DIMENSION(md1,md2,maxHOS) :: aHOS_ext
!
! Vertical derivatives modal coefficients
REAL(RP), DIMENSION(m1,m2,maxHOS)   :: kth, kth2
REAL(RP), DIMENSION(md1,md2,maxHOS) :: kth_ext, kth2_ext
!
! Useful mathematical functions
REAL(RP), DIMENSION(maxHOS) :: fact
!
! Pre-calculated wave number vectors, and trigonometric vectors and matrix
REAL(RP), DIMENSION(m3_add,m2)    :: k_add_thk_add 
REAL(RP), DIMENSION(m1,m2,m3_add) :: csh_add_x,k_add_sh_add_x, kycsh_add_x, kx_add_csh_add_x
!
! Wave-numbers
REAL(RP), DIMENSION(m1)        :: kx
REAL(RP), DIMENSION(m2)        :: ky
REAL(RP), DIMENSION(m3_add,m2) :: k2_add,k_add_2,k_add
REAL(RP), DIMENSION(m3_add)    :: kx_add
REAL(RP), DIMENSION(m1,m2)     :: k
REAL(RP), DIMENSION(m1)        :: kx2			
!
! Wave basin parameters
REAL(RP)          :: h,xlen,ylen
!
! Test case studied during computation
INTEGER           :: icase
!
! Sloshing case
INTEGER           :: islosh
REAL(RP)          :: aslosh
!
! Monochromatic wave case
INTEGER           :: ibat, ifree
REAL(RP)          :: amp_mono, nu_mono, theta_mono, ph_mono, T_mono, k_mono, lambda_mono, w_mono
REAL(RP)          :: xd_mono, z_dipmono
! 
! Frequency file case
CHARACTER(LEN=100) :: file_name
INTEGER            :: i_cut
REAL(RP)           :: nuc_low, nuc_high
!
! Wavemaker
INTEGER           :: iwidth, iramp, i_wmk, igeom, n_side
REAL(RP)          :: T_stop, Tramp, ywmk, yl, yr, y_rmp, d_hinge, x_dip
!
! Numerical beach
INTEGER                    :: iabsnb
REAL(RP)                   :: xabsf, xabsb, yabsl, yabsr, coeffabsf, coeffabsb, coeffabsl, coeffabsr
REAL(RP), DIMENSION(m1,m2) :: nu, dnux, dnuy
!
! Time ramp
REAL(RP)          :: alpha,alphat,alphatt
!
! Width ramp
REAL(RP), DIMENSION(m2) :: rampy, drampy
!
! Numerical probes (FS elevation and pressure)
INTEGER           :: iprobes, nprobes
CHARACTER(LEN=8)  :: pro_file
!
! Time stepping
REAL(RP)          :: time, delt
!
! Output files
INTEGER           :: idim, i3d, imodes, iwmk, i_sw, i_card
!
! Wave object
TYPE(wave) :: wav
!
! Space discretization
REAL(RP) :: delx,dely,delz
!
! Mesh coordinate
REAL(RP), DIMENSION(m1)    :: cx 
REAL(RP), DIMENSION(m2)    :: cy 
REAL(RP), DIMENSION(m3)    :: cz
REAL(RP), DIMENSION(m1,m2) :: x, y
!
! Spinning dipoles
COMPLEX(CP), ALLOCATABLE, DIMENSION(:,:) :: mu
REAL(RP), ALLOCATABLE, DIMENSION(:,:)    :: beta, beta_x, beta_y
REAL(RP), ALLOCATABLE, DIMENSION(:,:,:)  :: delta_x, delta_x2
INTEGER                                  :: n_beta, n_u, n_z_chk
!
! Energies and volume
REAL(RP) :: volfs, voltot, enepot, enekin, enetot
!
! Wavemaker position: 1st order and dervatives
REAL(RP), DIMENSION(m3,m2) :: pos1st,dpos1stdt,d2pos1stdt2
REAL(RP), DIMENSION(m3,m2) :: dpos1stdy,dpos1stdz,d2pos1stdydt,d2pos1stdzdt
!
! Wavemaker position: 2nd order and derivatives
REAL(RP), DIMENSION(m3,m2)     :: dpos2nddy,dpos2nddz,d2pos2nddydt,d2pos2nddzdt
COMPLEX(CP), DIMENSION(m3,m2)  :: pos2nd_ini
REAL(RP), DIMENSION(m3,m2)     :: pos2nd,dpos2nddt,d2pos2nddt2
!
! Interpolation polynom
REAL(RP) :: c_poly_add(2)		
!
! All position routines adapted for spectrum case
COMPLEX(CP), ALLOCATABLE, dimension(:,:) :: pos_ini, dposdy_ini
COMPLEX(CP), ALLOCATABLE, dimension(:)   :: pos_gzed, dpos_gzeddz
COMPLEX(CP), ALLOCATABLE, dimension(:)   :: pos_free_ini, dposdy2nd_ini
!
! Additional potential and derivatives on free surface (total, 1st, 2nd and 3rd order)
REAL(RP), DIMENSION(m1,m2) :: phiz_add, phit_add, phix_add, phiy_add, phi_add
REAL(RP), DIMENSION(m1,m2) :: phiz_add_1,phit_add_1,phix_add_1,phiy_add_1,phi_add_1
REAL(RP), DIMENSION(m1,m2) :: phiz_add_2,phit_add_2,phix_add_2,phiy_add_2,phi_add_2
REAL(RP), DIMENSION(m1,m2) :: phiz_add_3,phit_add_3,phix_add_3,phiy_add_3,phi_add_3
REAL(RP), DIMENSION(m1,m2) :: phizz_add_1,phixt_add_1,phiyt_add_1, phizt_add_1,phixx_add_1
REAL(RP), DIMENSION(m1,m2) :: phixz_add_1,phiyz_add_1
!
! Additional potential and derivatives on wavemaker
REAL(RP), DIMENSION(m3_add,m2) :: phi2xx_add_1st,phi2xxt_add_1st,phi2y_add_1st,phi2yt_add_1st, phi2z_add_1st, phi2zt_add_1st
!
! Additional potential - modal amplitudes (total, 1st, 2nd and 3rd order)+ ramp adjustment
REAL(RP) , DIMENSION(m3_add,m2):: a1st_add, a2nd_add, a3rd_add, a_add_ramp
!
! 1st order quantities 
REAL(RP), DIMENSION(m1,m2)  :: a1st, phi1st, phiz1st, eta1st 
!
! 2nd order quantities
REAL(RP), DIMENSION(m1,m2)  :: a2nd, phi2nd, phiz2nd, eta2nd
!
! derivatives of first-order quantities needed at second order
!
! free surface	
REAL(RP), DIMENSION(m1,m2)  :: etax1st,etay1st,phix1st,phiy1st
REAL(RP), DIMENSION(m1,m2)  :: phixx1st,phizt1st,phizz1st
!
! wavemaker
REAL(RP), DIMENSION(m3,m2)  :: phixx1st_add,phixxt1st_add
REAL(RP), DIMENSION(m3,m2)  :: phiy1st_add,phiyt1st_add,phiz1st_add,phizt1st_add
!
! derivatives needed to compute 3rd order wavemaker
REAL(RP), DIMENSION(m3,m2) :: phixx2nd_add,phixxt2nd_add,phixxx1st_add,phixxxt1st_add
REAL(RP), DIMENSION(m3,m2) :: phiy2nd_add,phiyt2nd_add,phixy1st_add,phixyt1st_add
REAL(RP), DIMENSION(m3,m2) :: phiz2nd_add,phizt2nd_add,phixz1st_add,phixzt1st_add
!
! Trigonometric functions on wavemaker
REAL(RP), DIMENSION(m1,m2,m3) :: kx2cshx_add, kycshx_add, kshx_add
!
! Temporal derivatives
REAL(RP), DIMENSION(m1,m2)      :: da1st, da2nd, deta1st, deta2nd
REAL(RP) , DIMENSION(m3_add,m2) :: da1st_add, da2nd_add, da3rd_add
!
! Pressure and velocities calculation
REAL(RP), DIMENSION(m1,m2,maxprobes) :: phix_addm, phiy_addm, phiz_addm
REAL(RP), DIMENSION(m1,m2,maxprobes) :: phit_addm, phi_addm
!
! Trigonometric functions for pressure calculation on probes
REAL(RP), DIMENSION(m3_add,m2,maxprobes) :: csh_add_x_probe,k_add_sh_add_x_probe
REAL(RP), DIMENSION(m3_add,m2,maxprobes) :: kycsh_add_x_probe, kx_add_csh_add_x_probe
!
! Pressure/velocities validation
REAL(RP), DIMENSION(m1,m2) :: vitxref_FS,vityref_FS, vitzref_FS,vitx2ref_FS, vitz2ref_FS, phitref_FS
!
! SWENSE outputs
REAL(RP), DIMENSION(m1,m2)     :: modesspecx, modesspecz, modesspecy, modesspect, modesFS, modesFSt
REAL(RP), DIMENSION(m3_add,m2) :: modesadd, modesaddt
INTEGER                        :: it !(increment to write the modes_HOS_swense.dat)
!
! Non-dimensional parameters
REAL(RP) :: x_adim, t_adim
!
! Trigonometric matrix for probes calculation
REAL(RP), DIMENSION(m1,m2,maxprobes) :: mat_coscos, mat_sincos, mat_cossin, mat_sinsin
!
! probe locations
REAL(RP), DIMENSION(maxprobes) :: xprobe,yprobe,zprobe
integer, DIMENSION(maxprobes)  :: i_xprobe,i_yprobe,i_zprobe
!
! Irregular wave generation in input file
REAL(RP) :: Hs, Tp, gamma
INTEGER  :: iseed
!
!
!
CONTAINS
!
!
!
FUNCTION extend(a)
!
IMPLICIT NONE
!
REAL(RP), DIMENSION(m1,m2)   :: a
REAL(RP), DIMENSION(md1,md2) :: extend
!
INTEGER :: i1,i2
!
do i1 = 1, n1
   do i2 = 1,n2
      extend(i1,i2) = a(i1,i2)
   end do
   do i2 = n2+1,nd2
      extend(i1,i2) = 0.0_rp
   end do
end do
do i1 = n1+1, nd1
   do i2 = 1,nd2
      extend(i1,i2) = 0.0_rp
   end do
end do
!
END FUNCTION extend
!
!
!
FUNCTION reduce(a_big)
!
IMPLICIT NONE
!
REAL(RP), DIMENSION(md1,md2) :: a_big
REAL(RP), DIMENSION(m1,m2)   :: reduce
!
reduce(1:n1,1:n2) = a_big(1:n1,1:n2)
!
END FUNCTION reduce
!
!
!
END MODULE variables