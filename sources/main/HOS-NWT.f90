!     =============
program HOS_NWT ! Numerical Wave Tank based on High-Order Spectral
!     =============
!
!
! This is the main program for HOS-NWT computations
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    Copyright (C) 2014 - LHEEA Lab., Ecole Centrale de Nantes, UMR CNRS 6598
!
!    This program is free software: you can redistribute it and/or modify
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
! start Main ***********************************************************
!
USE type
USE variables
USE vol_energy
USE wavemaker
USE fourier_FFTW
USE variablechange
USE runge_kutta
USE filtering
USE resol_HOS
USE resol_wmkr
USE dealiasing
USE input
USE initial_condition
USE output
!
!% VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     __________________________________________________________________
!
IMPLICIT NONE
!
!% LOCAL VARIABLES
! coordinate system on wavemaker
!
REAL(RP), DIMENSION(m3) :: cx_add
!
! wave numbers
REAL(RP), DIMENSION(m2)     :: ky2
REAL(RP), DIMENSION(m3_add) :: kx2_add
!
REAL(RP)         :: t_i, t_f
CHARACTER(LEN=1) :: t_unit
INTEGER          :: i1, i2, i3, j, ii
REAL(RP)         :: Tscale, xlen_add, pixlen, piylen, thk, pixlen_add
REAL(RP)         :: k_m, t_tot
INTEGER          :: endramp
INTEGER          :: n_hour, n_min, n_sec
!
TYPE(config)     :: cfg
!
REAL(RP) :: expon3, expon3p1, expon1, expon2, expon12, csh_add, sh_add, k_add_x, k_add_x1, k_add_x2, k_add_x_max
!
! Time-stepping
REAL(RP) :: dt_choice(3), dt_rk4(3), dt_ch, dt_stab, c_rk4, c_CFL
REAL(RP) :: h_rk, h_loc, error, delt_out
!
REAL(RP), DIMENSION(m1,m2)     :: eta_rk, phis_rk
REAL(RP)                       :: t_i_indiv, t_f_indiv
INTEGER                        :: n_rk, n_rk_tot, n_error, n_er_tot
REAL(RP)                       :: dt_correc
REAL(RP), DIMENSION(m3_add,m2) :: a_add_rk_1, a_add_rk_2, a_add_rkramp, a_add_rk_3
REAL(RP), DIMENSION(m1,m2)     :: a1st_rk, a2nd_rk, eta1st_rk, eta2nd_rk
!
! Trigonometric functions
REAL(RP) ::  cx_add_k ,csh , cshx_add_a 
REAL(RP), DIMENSION(m3_add) ::  cs_add, k_add_s_add, s_add
!
! Read input file
CALL read_input('input_HOS-NWT.dat')
!
! Initiate parameters
CALL initiate_parameters
!
! Initiate output
CALL init_output(i3d,imodes,iprobes,iwmk,igeom,i_sw,idim)
!     __________________________________________________________________
!
!% NONDIMENTIONALIZATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     __________________________________________________________________
!
!	general constants
Tscale = 1.d0 / dsqrt(g/h)
!
!	nondimensionalizations [L]=h ;[T]=1/sqrt(g/h)
xlen     = xlen / h
ylen     = ylen / h
ywmk     = ywmk / h
yl       = yl   / h
yr       = yr   / h
amp_mono = amp_mono / h
xlen_add = 2.d0 * REAL(l_add)		! influence of additional tank beam MUST be INVESTIGATED
xd_mono  = xd_mono / h
d_hinge  = d_hinge / h
x_dip    = x_dip   / h
aslosh   = aslosh  / h
!
Tramp   = Tramp  / Tscale
T_stop  = T_stop / Tscale
f_out   = f_out * Tscale
nu_mono = nu_mono  * Tscale
!
w_mono = nu_mono * TWOPI
T_mono = 1.0d0 / nu_mono
!     __________________________________________________________________
!
!% DISCRETIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     __________________________________________________________________
!
!     __________________________________________________________________
!
! Automatic choice of time step
!
WRITE(*,'(A)') 'RK4-requested time step: '  
!
k_m          = PI / xlen
dt_choice(1) = SQRT(k_m * xlen**2 / TANH(k_m)) / REAL(n1-1)
dt_rk4(1)    = SQRT(8 / PI**2) * xlen / REAL(n1-1)
!
IF (n2 /= 1) THEN
   k_m          = PI / ylen
   dt_choice(2) = SQRT(k_m * ylen**2 / TANH(k_m)) / REAL(n2-1)
   dt_rk4(2)    = SQRT(8 / PI**2) * ylen / REAL(n2-1)
ELSE
   dt_choice(2) = dt_choice(1) + 1.0_rp
   dt_rk4(2)    = dt_rk4(1)    + 1.0_rp
END IF
!
IF (igeom /= 3) THEN
   k_m          = PI / (2.0_rp * REAL(l_add))
   dt_choice(3) = SQRT(k_m * (2.0_rp * REAL(l_add))**2 / TANH(k_m)) / REAL( REAL(l_add) * (n3-1))
   dt_rk4(3)    = SQRT(8 / PI**2) * (2.0_rp * REAL(l_add)) / REAL( l_add * (n3-1))
ELSE
   dt_choice(3) = dt_choice(1) + 1.0_rp
   dt_rk4(3)    = dt_rk4(1)    + 1.0_rp
END IF
!
if (igeom==3 .OR. icase==1) then
   IF (n2 == 1) THEN
      dt_ch   = dt_choice(1)
      dt_stab = dt_rk4(1)
      write(*,'(A)') 'CFL x-modes limiting'
   ELSE
      dt_ch   = MIN(dt_choice(1), dt_choice(2))
      dt_stab = MIN(dt_rk4(1), dt_rk4(2))
      if (ABS(dt_choice(1)-dt_ch) < tiny) then
         write(*,'(A)') 'CFL x-modes limiting'
      else
         write(*,'(A)') 'CFL y-modes limiting'
      end if
   END IF
else if (n2.eq.1) then
   dt_ch   = MIN(dt_choice(1), dt_choice(3))
   dt_stab = MIN(dt_rk4(1), dt_rk4(3))
   if (ABS(dt_choice(1)-dt_ch) < tiny) then
      write(*,'(A)') 'CFL x-modes limiting'
   else
      write(*,'(A)') 'CFL z-modes limiting'
   end if
else
   dt_ch   = MIN(dt_choice(1), dt_choice(2), dt_choice(3))
   dt_stab = MIN(dt_rk4(1), dt_rk4(2), dt_rk4(3))
   if (ABS(dt_choice(1)-dt_ch) < tiny) then
      write(*,'(A)') 'CFL x-modes limiting'
   else if (ABS(dt_choice(2)-dt_ch) < tiny) then
      write(*,'(A)') 'CFL y-modes limiting'
   else
      write(*,'(A)') 'CFL z-modes limiting'
   end if
end if
!
write(*,'(A,4F10.6)') 'RK4 stability condition     => delt = ',dt_stab, dt_rk4
write(*,'(A,F10.6)')  'West et al.'' CFL condition => delt = ',dt_ch
!
! security factor
delt = 0.95_rp * dt_stab
!
delt_out = 1.0_rp / f_out
!
delt = MIN(delt, delt_out)
!
IF (ABS(delt-delt_out) < tiny) THEN
   WRITE(*,*) 'DECREASE f_out'
   STOP
END IF
!
if (ABS(dt_rk4(1)-dt_stab) < tiny) then
   c_rk4 = xlen / (REAL(n1-1) * dt_stab)
else if (ABS(dt_rk4(2)-dt_stab) < tiny) then
   c_rk4 = ylen / (REAL(n2-1) * dt_stab)
else
   c_rk4 = 2.0_rp * REAL(l_add) / (REAL(n3-1) * dt_stab)
end if
!
if (ABS(dt_choice(1)-dt_ch) < tiny) then
   c_CFL = xlen / (REAL(n1-1) * dt_ch)
else if (ABS(dt_choice(2)-dt_ch) < tiny) then
   c_CFL = ylen / (REAL(n2-1) * dt_ch)
else
   c_CFL = 2.0_rp * REAL(l_add) / (REAL(n3-1) * dt_ch)
end if
!
!   mesh generation
delx = xlen / (n1 - 1)
delz = 1.d0 / (n3 - 1)
!
if (n2.ne.1) then
   dely = ylen / (n2 - 1)
else
   dely = 0.d0
endif
!
do i2 = 1, n2
   cy(i2) = (i2 - 1) * dely
end do
!
do i1 = 1,n1
   cx(i1) = (i1 - 1) * delx
   do i2 = 1, n2
      x(i1,i2) = cx(i1)
      y(i1,i2) = cy(i2)
   end do
end do
!
do i3=1,n3
   cx_add(i3) = delz*(i3-1)
   cz(i3)   = cx_add(i3)-1.d0
end do
!
!% WAVEMAKER INITIALIZATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     __________________________________________________________________
!
CALL initialize_wav_cfg(icase,igeom,ibat,Tscale,cfg,wav)
!     __________________________________________________________________
!
!% ABSORBING SYSTEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     __________________________________________________________________     
!
IF(iabsnb /= 0) THEN
    call absorb
ELSE ! ensure it is initialized
    nu   = 0.0_rp
    dnux = 0.0_rp
    dnuy = 0.0_rp
ENDIF
!     __________________________________________________________________
!
!% SPECTRAL RESOLUTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     __________________________________________________________________
!
!	wave numbers
!
pixlen=pi/xlen
!
IF (n2 /= 1) THEN
   piylen=pi/ylen
ELSE
   piylen=0.d0
ENDIF
!
DO i2=1,n2
   ky(i2)  = (i2-1) * piylen
   ky2(i2) = ky(i2) * ky(i2)
ENDDO
!
DO i1 = 1, n1
   kx(i1)  = (i1-1) * pixlen
   kx2(i1) = kx(i1) * kx(i1)
   DO i2 = 1, n2
      k(i1,i2)     = SQRT(kx2(i1) + ky2(i2))
      omega(i1,i2) = SQRT(k(i1,i2)*TANH(k(i1,i2)))
   ENDDO
ENDDO
!
! HOS modal coefficients of the vertical derivatives
DO i1 = 1, n1
   DO i2 = 1,n2
      thk = TANH(k(i1,i2))
      DO j = 1,mHOS,2
         kth(i1,i2,j) = k(i1,i2)**j * thk
         IF((i1.EQ.1).AND.(i2.EQ.1)) THEN
            kth2(i1,i2,j) = 0.0_rp
         ELSE
            kth2(i1,i2,j) = k(i1,i2)**j / thk
         ENDIF
      ENDDO
      DO j=2,mHOS,2
         kth(i1,i2,j)  = k(i1,i2)**j
         kth2(i1,i2,j) = k(i1,i2)**j
      ENDDO
   ENDDO
ENDDO
!
! kth_ext
DO i1 = 1, Nd1
   DO i2 = 1, Nd2
      thk = TANH(SQRT(((i1 - 1) * pixlen)**2 + ((i2 - 1) * piylen)**2))
      DO ii = 1,mHOS,2
         kth_ext(i1,i2,ii) = (SQRT(((i1 - 1) * pixlen)**2 + ((i2 - 1) * piylen)**2))**ii * thk
         IF((i1.EQ.1).AND.(i2.EQ.1)) THEN
            kth2_ext(i1,i2,ii) = 0.0_rp
         ELSE
            kth2_ext(i1,i2,ii) = (SQRT(((i1 - 1) * pixlen)**2 + ((i2 - 1) * piylen)**2))**ii / thk
         ENDIF
      ENDDO
      DO ii=2,mHOS,2
         kth_ext(i1,i2,ii)  = (SQRT(((i1 - 1) * pixlen)**2 + ((i2 - 1) * piylen)**2))**ii
         kth2_ext(i1,i2,ii) = kth_ext(i1,i2,ii)
      ENDDO
   ENDDO
ENDDO
!
! on the wavemaker surface
pixlen_add = PI / xlen_add
!
! wavemaker modes part
DO i3 = 1, n3_add
   kx_add(i3)  = (2*i3-1) * pixlen_add
   kx2_add(i3) = kx_add(i3) * kx_add(i3)
   DO i2 = 1, n2
      k_add(i3,i2)  = SQRT(kx2_add(i3) + ky2(i2))
      k_add_2(i3,i2) = 2.0_rp * k_add(i3,i2) ! GD : to check
      k_add_thk_add(i3,i2) = k_add(i3,i2) * dtanh(k_add(i3,i2) * xlen) ! on the wavemaker
   ENDDO
ENDDO
!
! known spectral matrix or vectors needed
! free surface modes part
DO i1= 1, n1
   DO i2 = 1, n2
      DO i3=1,n3				   ! on the wavemaker
         cx_add_k   = k(i1,i2) * cx_add(i3)
         csh        = COSH(k(i1,i2))
         cshx_add_a = COSH(cx_add_k) / csh
         kx2cshx_add(i1,i2,i3) = cshx_add_a * kx2(i1)
         kycshx_add(i1,i2,i3)  = cshx_add_a * ky(i2)
         kshx_add(i1,i2,i3)    = k(i1,i2) * SINH(cx_add_k) / csh
      ENDDO
   ENDDO
ENDDO
!	    
! wavemaker modes part
k_add_X_max = 700.0_rp
!
DO i3 = 1, n3_add
   cs_add(i3)   = COS(kx_add(i3))
   s_add(i3)    = SIN(kx_add(i3))
   k_add_s_add(i3) = kx_add(i3) * SIN(kx_add(i3))
   !
   DO i2 = 1, n2
      k_add_x = k_add_2(i3,i2) * xlen
      !
      IF (k_add_x <= k_add_X_max) THEN
         expon3 = EXP(-k_add_2(i3,i2) * xlen)
      ELSE
         expon3 = 0.0_rp  
      END IF
      !
      expon3p1  = expon3 + 1.d0
      DO i1 = 1, n1	                             ! on the free surface
         k_add_x1 = k_add(i3,i2) * cx(i1)
         !
         IF (k_add_x1 <= k_add_X_max) THEN
            expon1    = EXP(-k_add_x1) / expon3p1
         ELSE
            expon1 = 0.0_rp  
         END IF
         !
         k_add_x2 = k_add_2(i3,i2) * (xlen - cx(i1))
         IF (k_add_x2 <= k_add_X_max) THEN
            expon2    = EXP(- k_add_x2)
         ELSE
            expon2 = 0.0_rp  
         END IF
         !
         IF (k_add_x1 + k_add_x2 <= k_add_x_max) THEN
            expon12   = expon1 * expon2
         ELSE
            expon12   = 0.0_rp
         END IF
         !
         csh_add     = expon1 + expon12
         sh_add      = expon1 - expon12
         csh_add_x(i1,i2,i3)     = csh_add
         k_add_sh_add_x(i1,i2,i3)   = sh_add * k_add(i3,i2)
         kycsh_add_x(i1,i2,i3)   = csh_add * ky(i2)
         kx_add_csh_add_x(i1,i2,i3) = csh_add * kx_add(i3)
      ENDDO
   ENDDO
ENDDO
!     __________________________________________________________________
!      
!% INITIAL CONDITIONS : FS + time derivatives %%%%%%%%%%%%%%%%%%%%%%%%%%
!     __________________________________________________________________
!
CALL initialize
!
! Initial volume and energy
CALL volume
CALL energy(eta, phis+phi_add, deta)
!     __________________________________________________________________
!
!% OUTPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     __________________________________________________________________      
!
IF (idim == 0) THEN
   t_adim = 1.0_rp
   x_adim = 1.0_rp
   t_unit = ''
ELSE IF (idim == 1) THEN
   t_adim = Tscale
   x_adim = h
   t_unit = 's'
END IF
!     __________________________________________________________________
!
!     ***************************************************************
!      
!     BIG TIME LOOP
!      
!     ***************************************************************
!
CALL CPU_TIME(t_i)
!
n_rk     = 0
n_rk_tot = 0
n_error  = 0
n_er_tot = 0
!
DO WHILE (time_cur < T_stop + delt_out)        ! instead of DO WHILE (time_cur < T_stop) 
                                               ! because output stuff is before RK stuff
   n_er_tot = n_error + n_er_tot
   n_rk_tot = n_rk_tot + n_rk
   it =it+1
   !
   CALL CPU_TIME(t_f)
   !
   IF (time_cur <= 100 * delt_out) THEN
      if ((t_f-t_i).gt.2.0d-2) WRITE(*,'(A,F6.2,2(X,I3))') 'CPU time for Output time step ',t_f-t_i, n_rk, n_error
   END IF
   !
   ! CPU estimation at beginning of computation
   if (time_cur <= 10*delt_out) then
      t_tot  = (t_f - t_i) * T_stop * f_out
      n_hour = FLOOR(t_tot / 3600)
      n_min  = FLOOR((t_tot - 3600 * n_hour) / 60)
      n_sec  = FLOOR(t_tot - 60* (n_min + 60 * n_hour))
      write(*,'(A,I3,A,I2,A,I2,A)') 'Expected time ',n_hour,' h ',n_min, ' min ',n_sec,' sec.'
   end if
   t_i = t_f
   !
   !     ________________________________________________________
   ! ccccccccccc WRITING ZONE cccccccccccccccccccccccccccccccccccc!      
   !     ________________________________________________________
   !
   ! Write output files
   CALL output_time_step(i3d,imodes,iprobes,iwmk,igeom,i_sw,time_cur,t_unit)
   !
   ! Display on screen
	IF (n2 /= 1) THEN
	  write(*,200) voltot*x_adim*x_adim, time_cur*t_adim, delt*t_adim,n_rk_tot,n_er_tot
	ELSE
	  write(*,200) voltot*x_adim, time_cur*t_adim, delt*t_adim,n_rk_tot,n_er_tot
	END IF
	200 format('    volume = ',1ES12.4,' t= ',1ES11.4,' timestep:',ES11.4,2(X,I5)) 
   !
   !ccccccccccc END OF WRITING ZONE ccccccccccccccccccccccccccccc!
   !
   IF (ABS(time_cur - T_stop) <= tiny) EXIT ! output of the last zone is done
   !     ________________________________________________________
   !
   !      =======================================
   !      =  TIME INTEGRATION (Runge-Kutta 4th) =
   !      =======================================
   !
   time_next = time_cur + delt_out
   IF (time_next > T_stop + delt_out) time_next = T_stop + delt_out
   h_rk    = delt
   n_rk    = 0
   n_error = 0
   !
   DO WHILE(time_cur < time_next)
      !
      CALL CPU_TIME(t_i_indiv)
      h_loc = h_rk
      !
      IF (time_next - time_cur < h_rk) h_loc = time_next - time_cur
      !
      ! For time marching
      eta_rk(1:n1,1:n2)  = eta(1:n1,1:n2)
      phis_rk(1:n1,1:n2) = phis(1:n1,1:n2)
      a_add_rk_1(1:n3_add,1:n2)   = a1st_add(1:n3_add,1:n2)
      a_add_rk_2(1:n3_add,1:n2)   = a2nd_add(1:n3_add,1:n2)
      a_add_rkramp(1:n3_add,1:n2) = a_add_ramp(1:n3_add,1:n2)
      a_add_rk_3(1:n3_add,1:n2)   = a3rd_add(1:n3_add,1:n2)
      a1st_rk(1:n1,1:n2)   = a1st(1:n1,1:n2)
      a2nd_rk(1:n1,1:n2)   = a2nd(1:n1,1:n2)
      eta1st_rk(1:n1,1:n2) = eta1st(1:n1,1:n2)
      eta2nd_rk(1:n1,1:n2) = eta2nd(1:n1,1:n2)
      !
      ! adaptative time step
      call runge4_adapt2(time_cur, eta_rk, phis_rk, a_add_rk_1, a_add_rk_2, &
           a_add_rkramp,a_add_rk_3, eta1st_rk,eta2nd_rk, a1st_rk, a2nd_rk, &
           h_loc, error, tiny)
      !
      IF (error > toler .AND. h_rk > 0.001_rp * 0.95_rp * dt_stab) THEN
         !
         ! Reduction of time step 
         dt_correc = 0.9*(error/toler)**(-(1.0/RK_p))
         n_error = n_error + 1
         !
         if((ABS(Tramp-(time_cur+h_loc)).LE.tiny).AND.(Tramp.GT.time_cur)) then
            endramp = 0
         endif
      ELSE
         !
         ! Quantities are time marched
         eta(1:n1,1:n2)  = eta_rk(1:n1,1:n2)
         phis(1:n1,1:n2) = phis_rk(1:n1,1:n2)
         a1st_add(1:n3_add,1:n2) = a_add_rk_1(1:n3_add,1:n2)
         a2nd_add(1:n3_add,1:n2) = a_add_rk_2(1:n3_add,1:n2)
         a_add_ramp(1:n3_add,1:n2) = a_add_rkramp(1:n3_add,1:n2)
         a3rd_add(1:n3_add,1:n2) = a_add_rk_3(1:n3_add,1:n2)
         a1st(1:n1,1:n2)   = a1st_rk(1:n1,1:n2)
         a2nd(1:n1,1:n2)   = a2nd_rk(1:n1,1:n2)
         eta1st(1:n1,1:n2) = eta1st_rk(1:n1,1:n2)
         eta2nd(1:n1,1:n2) = eta2nd_rk(1:n1,1:n2)
         !
         ! Correction wavemaker ramp
         if((iramp == 1).AND.(Tramp-(time_cur+h_loc).LE.tiny).AND.(Tramp.GT.time_cur).AND.(endramp.NE.1)) then
              a1st_add(1:n3_add,1:n2) = a1st_add(1:n3_add,1:n2) - a_add_ramp(1:n3_add,1:n2)
              endramp=1
              write(*,*) '********************* end of time-ramp *********************'
         endif
         time_cur = time_cur + h_loc
         !
         ! Evaluation next time-step
         IF (ABS(error) < tiny) THEN
         	dt_correc = 2.0_rp ! Give a value greater than 1 so that it uses maximum possible time-step
         ELSE
         	dt_correc = 0.9_rp*(error/toler)**(-(1.0_rp/REAL(RK_p)))
         ENDIF
         n_rk = n_rk + 1
      END IF
      h_rk = h_rk * dt_correc
      h_rk = MIN(h_rk, 0.95_rp * dt_stab, delt_out)
      h_rk = MAX(h_rk, 0.0001_rp * 0.95_rp * dt_stab, 0.0001_rp * delt_out)
      !
      CALL CPU_TIME(t_f_indiv)
      !
      IF (time_cur <= 10 * delt_out) THEN
         if ((t_f_indiv-t_i_indiv).gt.2.0d-2) WRITE(*,'(A,F6.2,2(X,ES11.4))') 'CPU time for individual time step ' &
              ,t_f_indiv-t_i_indiv, h_rk, error
      END IF
   END DO
   delt = h_rk
   !
   ! Check the slope of FS elevation: possible filtering if above threshold
   CALL check_slope
   !
   ! Volume and energy computation
   CALL volume
   CALL energy(eta, phis+phi_add, deta)
   !
END DO ! while	
!				
!     *****************************************************************
!	END OF BIG TIME LOOP
!     *****************************************************************
!
end program HOS_NWT
!
! end Main ************************************************************
!