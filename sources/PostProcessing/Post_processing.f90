PROGRAM Post_processing
!
! This is the main program for post-processing HOS-NWT computations
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
USE input_post_process
USE read_files
USE output_post_process
USE reconstruction
USE fourier_FFTW
!
IMPLICIT NONE
!
REAL(RP) :: time, time_prev
!
! For SWENSE-type outputs + velocity/pressure cards
!
!INTEGER :: nz
REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: modesspecx,modesspecy,modesspecz,modesspect,modesFS,modesFSt,modesadd,modesaddt
REAL(RP), ALLOCATABLE, DIMENSION(:)   :: xvect,yvect,zvect
REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: vitx,vity,vitz,phit,dudt,dvdt,dwdt,zlocal
REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: phixadd,phizadd,phiyadd,phitadd,phixtadd,phiytadd,phiztadd
!
REAL(RP) :: dt_out_star,T_stop_star,xlen_star,ylen_star,depth
INTEGER  :: imin,imax,jmin,jmax,i_xvect,i_yvect
INTEGER  :: i1,i2,i3,i_test
!
REAL(RP) :: T_adim, L_adim
REAL(RP) :: tiny_sp
!
! tiny_sp is single precision: useful for inequalities check with values read from files
tiny_sp = epsilon(1.0)
!
! Read input file to define parameters of post-processing
CALL read_input('input_post_process.dat')
!
! Initialize outputs
CALL init_output_post_process(i_card)
!
! Velocities and pressure inside domain
IF (i_card /= 0) THEN
	!
	i_unit = 201
	!
	! Initialize computations of volumic informations
	! Everything is non-dimensional in file_mod
	!
	CALL recons_HOS_init(file_mod,i_unit,n1,n2,n3_add,dt_out_star,T_stop_star,xlen_star,ylen_star,depth, &
		modesspecx,modesspecy,modesspecz,modesspect,modesFS,modesFSt,modesadd,modesaddt)
	!
	! Initialize Fourier (global variables have to be defined)
	m1      = n1
	m2      = n2
	m3_add  = n3_add
	Nd1     = 1
	Nd2     = 1
	md1     = 1
	md2     = 1
	!
	! Initialize Fourier transforms (FFTW)
	CALL fourier_ini(3)
	!
	! Define the scales to rebuild at HOS-NWT scale: may be changed to fit the scale for coupling for instance...
	L_adim = depth
	T_adim = SQRT(depth/g)
	!
	! Check (x_min, x_max, y_min, y_max) w.r.t. domain size
	! + time window (t_min, t_max)
	CALL check_sizes(n2,x_min,x_max,y_min,y_max,z_min,T_start,T_stop,xlen_star,ylen_star,T_stop_star,L_adim,T_adim)
	!
	ALLOCATE(x(n1),y(n2),ky(n2),kx(n1),k(n1,n2),kth(n1,n2),eta(n1,n2))
	!
	ALLOCATE(kx_add(n3_add),kx2_add(n3_add),k_add(n3_add,n2),k_add_2(n3_add,n2),k_add_thk_add(n3_add,n2))
	ALLOCATE(csh_add_x(n1,n2,n3_add),k_add_sh_add_x(n1,n2,n3_add),kycsh_add_x(n1,n2,n3_add),kx_add_csh_add_x(n1,n2,n3_add))
	!
	! Initialize mesh in physical and modal space (whole domain in HOS-NWT)
	CALL build_mesh_global(xlen_star,ylen_star,n1,n2,n3,n3_add,x,y,kx,ky,k,kth,&
		kx_add, kx2_add,k_add, k_add_2, k_add_thk_add,k_add_sh_add_x, kycsh_add_x, kx_add_csh_add_x)
	!
	! Define local meshes for zone of study
	CALL build_mesh_local(x_min,x_max,y_min,y_max,z_min,z_max,xlen_star,ylen_star,L_adim,n1,n2,i_zvect, &
		xvect,yvect,zvect,imin,imax,jmin,jmax)
	!
	! Reconstruction of fields
	! First ALLOCATE the matrices
	i_xvect = imax-imin+1
	i_yvect = jmax-jmin+1
	!
	ALLOCATE(zlocal(i_xvect, i_yvect), vitx(i_xvect, i_yvect), vity(i_xvect, i_yvect), vitz(i_xvect, i_yvect), &
		phit(i_xvect, i_yvect), dudt(i_xvect, i_yvect), dvdt(i_xvect, i_yvect), dwdt(i_xvect, i_yvect))
	!
	ALLOCATE(phixadd(i_xvect, i_yvect), phizadd(i_xvect, i_yvect),phiyadd(i_xvect, i_yvect), &
     phitadd(i_xvect, i_yvect), phixtadd(i_xvect, i_yvect), phiytadd(i_xvect, i_yvect), phiztadd(i_xvect, i_yvect))
	!
	! Define first time as the closest to T_start (input file)
	time      = NINT(T_start/T_adim/dt_out_star)*dt_out_star 
	time_prev = 0.0_rp
	!
	DO WHILE (time*T_adim <= T_stop+tiny_sp)
		!
		write(*,'(A,ES8.1)') 'time = ',time*T_adim
		! It reads the corresponding time in file_mod (closest to time)
		IF (time >= dt_out_star/2) THEN
			CALL read_mod(file_mod,i_unit,time,dt_out_star,n1,n2,n3_add, &
				modesspecx,modesspecy,modesspecz,modesspect,modesFS,modesFSt,modesadd,modesaddt)
		ENDIF
		!
		! For outputs
		eta = fourier_2_space(modesFS,'cos','cos')
		!
		! Make a loop over all the elements in z
		DO i3 = 1, i_zvect
			IF (i_card == 1) THEN
				! Construct the field at each zvect...
				CALL reconstruction_FFTs(modesspecx,modesspecy,modesspecz,modesspect,eta,&
						imin,imax,jmin,jmax,zvect(i3),vitx,vity,vitz,phit,dudt,dvdt,dwdt)
				! Additional contribution of wavemaker
				CALL reconstruction_add_FFTs(modesadd,modesaddt,eta,imin,imax,jmin,jmax,zvect(i3),&
						phixadd,phiyadd,phizadd,phitadd,phixtadd,phiytadd,phiztadd)
				!
				! vitx, vity, vitz, phit : total quantities
				vitx = vitx + phixadd
				vity = vity + phiyadd
				vitz = vitz + phizadd
				phit = phit + phitadd
				! dudt, dvdt, dwdt : total quantities
				dudt = dudt + phixtadd
				dvdt = dvdt + phiytadd
				dwdt = dwdt + phiztadd
				! Necessary for output
				DO i1=1, imax-imin+1
					DO i2=1, jmax-jmin+1
						zlocal(i1,i2) = zvect(i3)
					ENDDO
				ENDDO
			ELSEIF (i_card == 2) THEN
				! Construct the field at each zvect...
				CALL reconstruction_direct(modesspecx,modesspecy,modesspecz,modesspect,eta, &
						imin,imax,jmin,jmax,z_min/L_adim,i3,i_zvect,vitx,vity,vitz,phit,dudt,dvdt,dwdt,zlocal)
				! Additional contribution of wavemaker
				CALL reconstruction_add_direct(modesadd,modesaddt,imin,imax,jmin,jmax,zlocal,&
						phixadd,phiyadd,phizadd,phitadd,phixtadd,phiytadd,phiztadd)
				!
				! vitx, vity, vitz, phit : total quantities
				vitx = vitx + phixadd
				vity = vity + phiyadd
				vitz = vitz + phizadd
				phit = phit + phitadd
				! dudt, dvdt, dwdt : total quantities
				dudt = dudt + phixtadd
				dvdt = dvdt + phiytadd
				dwdt = dwdt + phiztadd
			ENDIF
			! Test to know if it is first z-element to write corresponding header
			IF (i3 == 1) THEN
				i_test = 1
			ELSE
				i_test = 0
			ENDIF
			! Output time step
			CALL output_time_step_card(i_card,tecplot,time,dt_out_star,zlocal,z_min,z_max,T_start,L_adim,T_adim,i_test,&
				imin,imax,jmin,jmax,i_zvect,vitx,vity,vitz,phit,eta)
		ENDDO
		!
		! next time-step
		time_prev = time
		time      = time + dt_out_star
		!
	ENDDO
	!
	! Close all files (including those open in output...
	CLOSE(i_unit)
	CLOSE(30)
	CLOSE(31)
	CLOSE(32)
ENDIF
!
! End of main program 
!
END PROGRAM Post_processing