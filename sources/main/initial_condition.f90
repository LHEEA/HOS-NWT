MODULE initial_condition
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
USE runge_kutta
USE linear_wave
USE read_ocean_txt
!
CONTAINS
!
SUBROUTINE initiate_parameters()
!
IMPLICIT NONE
!
REAL(RP) :: kp
!
! Checking dimensions
IF (n1 .gt. m1) THEN
   write(*,'(A,I4,A)') 'Error in HOS: m1 = ',m1,' too small'
   STOP 'm1'
ELSE IF (n2 .gt. m2) THEN
   write(*,'(A,I4,A)') 'Error in HOS: m2 = ',m2,' too small'
   STOP 'm2'
ELSE IF (n3 .gt. m3) THEN
   write(*,'(A,I4,A)') 'Error in HOS: m3 = ',m3,' too small'
   STOP 'm3'
ELSE IF (mHOS .gt. maxHOS) THEN
   write(*,'(A,I4,A)') 'Error in HOS: maxHOS = ',maxHOS,' too small'
   STOP 'maxHOS'
ELSE IF (nprobes .gt. maxprobes) THEN
   write(*,'(A,I4,A)') 'Error in HOS: maxprobes = ',maxprobes,' too small'
   STOP 'maxprobes'
ELSEIF (p1 .gt. p1max) THEN
   write(*,'(A,I4,A)') 'Error in p1: p1max = ',p1max,' too small'
   STOP 'p1'
ELSE IF (p2 .gt. p2max) THEN
   write(*,'(A,I4,A)') 'Error in p2: p2max = ',p1max,' too small'
   STOP 'p2'
END IF
!
! RK time stepping
CALL fill_butcher_array()
!
! Number of vertical points for dipole images checking
IF (igeom == 3) THEN
   n_z_chk = 10
END IF
!
WRITE(*,'(A)',ADVANCE='NO') 'Choosen sea test : '
SELECT CASE (icase)
	CASE (1)
   		WRITE(*,'(A,I2,A,F6.3,A)') 'sloshing case with mode ', islosh-1, ' and amplitude ', aslosh, ' m.'
   		k_mono = (islosh - 1) * PI / xlen
   		w_mono = SQRT(g*k_mono * tanh(k_mono * h))
   		nu_mono = w_mono / TWOPI
   		WRITE(*,'(A,F7.4,A,F7.4,A)') 'Frequency =', nu_mono, ' Hz. (', nu_mono*dsqrt(h/g), ' adim)'
   		WRITE(*,'(A,F5.2,A,F7.4,A)') 'Linear steepness =', 100.0_rp * k_mono * aslosh / PI, ' % or ',k_mono * aslosh, ' %'
	CASE (2)
   		WRITE(*,'(A)') 'monochromatic case'
   	CASE (3)
   		WRITE(*,'(2A)') 'Wavemaker movement from .dat & .cfg files: ', TRIM(file_name)
	CASE (31)
   		WRITE(*,'(2A)') 'OCEAN file case - ECN Wave Tank, with file: ', TRIM(file_name)
	CASE (32)
   		WRITE(*,'(2A)') 'OCEAN file case - ECN Towing Tank, with file: ', TRIM(file_name)
   	CASE (33)
   		WRITE(*,'(2A)') 'OCEAN file case - other wave tank, with file: ', TRIM(file_name)
   	CASE (4)
   		WRITE(*,'(A)') 'Long-crested irregular wave with JONSWAP spectrum'
   	CASE (41)
   		WRITE(*,'(A)') 'Long-crested irregular wave with Bretschneider spectrum'
   	CASE (42) !FIXME
   		WRITE(*,'(A)') 'Short-crested irregular wave with JOSWAP spectrum'
   		WRITE(*,'(A)') 'Not done yet!'
   		STOP
   	CASE DEFAULT
   		WRITE(*,'(A)') 'Test case unknown... check icase'
   		STOP
END SELECT
!
! Initializations of FFTW
CALL Fourier_ini
!
! Checking numerical parameters
CALL check_range(mHOS,'M ')
CALL check_range(p1,'p1')
CALL check_range(p2,'p2')
!
IF (p1 > mHOS) STOP 'initiate_parameters: p1 is higher than M'
IF (p2 > mHOS) STOP 'initiate_parameters: p2 is higher than M'
!
! Test steepness of sea state
!
!	generated wavelength and wave number calculation
IF ((icase == 4).OR.(icase == 41).OR.(icase == 42)) THEN
	kp = wave_number_r(1.0_rp/Tp, h, 1.d-12)
	IF ((kp*Hs/(2.0_rp*sqrt(2.0_rp))).GT.0.15_rp) THEN
		print*, 'Mean steepness kp*Hs/(2*sqrt(2)) is greater than 15%, possible breaking'
		print*, 'steepness =', kp*Hs/(2.0_rp*sqrt(2.0_rp))
	ENDIF
ENDIF
!
END SUBROUTINE initiate_parameters
!
!
!
SUBROUTINE initialize_wav_cfg(icase,igeom,ibat,Tscale,cfg,wav)
!
IMPLICIT NONE
!
! Inputs/ outputs
INTEGER, INTENT(INOUT)        :: icase
INTEGER, INTENT(IN)           :: igeom,ibat
REAL(RP), INTENT(IN)          :: Tscale
TYPE(config), INTENT(OUT)     :: cfg
TYPE(wave), INTENT(OUT)       :: wav
!
! Local variables
REAL(RP)     :: Hs_int, f_base, sigma
INTEGER      :: i1, iloop, n_pts
TYPE(wave)   :: wav_tmp1, wav_tmp2
!
INTEGER      :: rnum, harmo
REAL(RP)     :: clock
!
IF (icase == 2) THEN
   ! For monochromatic wave, define similar objects to icase=3/31/32/...
   !
   ! Config object
   rnum  = 1
   clock = 2.0_rp*nu_mono ! Already non-dimensional
   harmo = 1
   !
   IF (igeom == 2) THEN
      cfg = build_config(rnum, clock = clock, depth = h, typ_wmk= 'hinged', hinge = d_hinge, Ly = ylen)
   ELSE IF (igeom == 1) THEN
      cfg = build_config(rnum, clock = clock, depth = h, typ_wmk= 'piston', Ly = ylen)
   ELSE IF (igeom == 3) THEN
      cfg = build_config(rnum, clock = clock, depth = h, typ_wmk= 'dipole', Ly = ylen)
   END IF
   ! Wave object
   IF (ibat == 2) THEN
      wav = front(cfg, harmo, amp_mono, theta_mono, ph_mono)
   ELSE IF (ibat == 3) THEN
      wav = dalrymple(cfg, xd_mono, harmo, amp_mono, theta_mono, ph_mono)
   END IF
   IF (igeom == 3) THEN
      CALL dipoles_ini
   ELSE
      CALL wmk_ini
   END IF
   CALL check_mods(nu_mono, theta_mono)
   !	generated wavelength and wave number calculation
   k_mono      = wav%k(1) 
   lambda_mono = TWOPI / k_mono
ELSE IF (icase == 3 .OR. icase == 31 .OR. icase == 32 .OR. icase == 33) THEN
   !
   ! Read .txt file from Ocean output (Edimburgh Design wavemaker control software) 
   IF (icase == 31 .OR. icase == 32) THEN
      CALL ocean_txt2cfg_dat(file_name, i_cut, nuc_low, nuc_high, icase)
   ELSEIF (icase == 33) THEN ! Specified characteristics of wave tank
      ! FIXME, one has to define clock... More comprehensive variable?
      clock = xd_mono
      CALL ocean_txt2cfg_dat(file_name, i_cut, nuc_low, nuc_high, icase, clock, h, ylen, d_hinge, igeom)
   ENDIF
   wav       = file2wave(file_name)
   IF (ABS(wav%config%depth-h) > tiny)        WRITE(*,*) 'Warning : Different basin depth'
   IF (ABS(wav%config%Ly-ylen) > tiny)        WRITE(*,*) 'Warning : Different basin width'
   IF (ABS(wav%config%hinge-d_hinge) > tiny)  WRITE(*,*) 'Warning : Different hinge heights'
   !
   IF (ibat == 3.AND.n2/=1) THEN
      DO iloop = 1, SIZE(wav%ampli)
         wav_tmp1 = dalrymple(wav%config, xd_mono, wav%harmo(iloop)&
              , wav%ampli(iloop), wav%angle(iloop), wav%phase(iloop))
         wav_tmp2 = wav_tmp2 + wav_tmp1
      END DO
      wav = wav_tmp2
   END IF
   CALL display(wav)
   !
   IF ((igeom == 2).OR.(igeom == 1)) THEN
      CALL wmk_ini
   ELSE
      CALL dipoles_ini
   END IF
   call check_mods(MAXVAL(wav%freq))
   !	generated wavelength and wave number calculation
   k_mono      = wav%k(1) 
   IF (ABS(k_mono) < tiny) k_mono = wav%k(2)
   lambda_mono = TWOPI / k_mono
ELSE IF (icase == 4.OR.icase==41.OR.icase==42) THEN !Irregular wave
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! equivalent to file .cfg
    !
    wav%config%clock       = 32.0_rp
    wav%config%rnum        = CEILING(LOG10(T_stop*Tscale*wav%config%clock)/LOG10(2.0_rp)) !superior round up
    wav%config%depth       = h
    wav%config%ramp        = 3.0_rp*wav%config%clock !3 seconds
    wav%config%typ_ramp    = 'lin'
    wav%config%cutoff_low  = 0.0_rp
    wav%config%cutoff_high = 3.0_rp/Tp ! dimensional
    wav%config%Ly          = ylen*h    ! dimensional
    IF (igeom ==1) THEN
        wav%config%typ_wmk = 'piston'
    ELSEIF (igeom == 2) THEN
        wav%config%typ_wmk = 'hinged'
    ELSE
        print*, 'igeom should be 1 or 2'
        stop
    ENDIF
    wav%config%hinge = d_hinge*h !dimensional
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! equivalent to file .dat
    ! construction of spectrum
    !
    n_pts  = NINT(2*2**wav%config%rnum/wav%config%clock+1)
    f_base = REAL(wav%config%clock)/2.0_rp**wav%config%rnum
    ALLOCATE(wav%harmo(n_pts), wav%freq(n_pts), wav%ampli(n_pts), &
        wav%angle(n_pts), wav%phase(n_pts))
    !
    ! Initialize random number generation
    !
    CALL SRAND(iseed)
    Hs_int = 0.0_rp
    !
    DO i1=1,n_pts
        wav%harmo(i1) = i1-1
        wav%freq(i1)  = i1*f_base
        wav%angle(i1) = 0.0_rp
        wav%phase(i1) = TWOPI*RAND() !TWOPI*RAN(iseed)
        !
        IF (wav%freq(i1).LT.1.0_rp/Tp) THEN
            sigma = 0.07
        ELSE
            sigma = 0.09
        ENDIF
        ! S(f) df
        IF (icase == 4) THEN ! JONSWAP
            wav%ampli(i1) = 1.0_rp/(TWOPI**4)*g**2/wav%freq(i1)**5 &
                * exp(-5.0_rp/4.0_rp*(wav%freq(i1)*Tp)**(-4)) &
                * gamma**(exp((-(wav%freq(i1)-1.0_rp/Tp)**2)/(2.0_rp*(sigma/Tp)**2))) &
                * f_base
        ELSEIF (icase == 41) THEN !Bretschneider
            wav%ampli(i1) = 5.0_rp/16.0_rp/(Tp**4*wav%freq(i1)**5)*Hs*exp(-5.0_rp/4.0_rp*(wav%freq(i1)*Tp)**(-4))*f_base
        ELSEIF (icase == 42) THEN ! Directional spectrum
            print*,'case not available'            
            stop
        ENDIF
        !
        ! Compute integral to adjust Hs
        !
        Hs_int = Hs_int + wav%ampli(i1)
        wav%ampli(i1) = SQRT(2.0_rp*wav%ampli(i1))
    ENDDO
    Hs_int = 4.0_rp*SQRT(Hs_int)
    print*,'Hs target =',Hs ,'Tp =', Tp, 'Hs_calc =', Hs_int
    !
    ! Adjustment of amplitude to have correct Hs
    !
    wav%ampli = wav%ampli*Hs/Hs_int
    !
    ! Making wave data nondimensionalised
    wav%freq  = wav%freq * Tscale
    wav%ampli = wav%ampli / h
    wav%angle = wav%angle * pi / 180.0_rp
    !
    ! Non-dimensionalization
    !
    wav%config%clock       = wav%config%clock  * Tscale
    wav%config%cutoff_low  = wav%config%cutoff_low * Tscale
    wav%config%cutoff_high = wav%config%cutoff_high * Tscale
    wav%config%Ly          = wav%config%Ly / h
    wav%config%hinge       = wav%config%hinge / h
    !
    ! Building wavenumber and Wavemaker Transfer Function at each frequency
    wav = wave_init(wav)
    !
    CALL wmk_ini
    !call check_mods(MAXVAL(wav%freq))
    call check_mods(2.0_rp/Tp * Tscale)
    !	generated wavelength and wave number calculation
    k_mono      = wav%k(1) 
    IF (ABS(k_mono) < tiny) k_mono = wav%k(2)
    lambda_mono = TWOPI / k_mono
END IF
!
END SUBROUTINE initialize_wav_cfg
!
!
!
SUBROUTINE initialize
!
IMPLICIT NONE
!
! Local variables
INTEGER                        :: i1, i2, i3, j, isubsteprk
REAL(RP), DIMENSION(m1,m2)     :: dphis_rk
!
time_cur = 0.0_rp
time     = 0.0_rp
!
! number of time steps for modes_HOS_SWENSE.dat
it=1
!
!     Initial water elevations (sloshing)
IF (icase == 1) THEN
   DO i1=1,n1
      DO i2=1,n2
         eta(i1,i2) = aslosh * cos(kx(islosh)*cx(i1)) ! * cos(ky(islosh)*cy(i2))
      ENDDO
   ENDDO
ELSE
   DO i1=1,n1
      DO i2=1,n2
         eta(i1,i2) = 0.0_rp
      ENDDO
   ENDDO
END IF
!
!     Initial potentials
DO i1=1,n1
   DO i2=1,n2
      phis(i1,i2) = 0.0_rp
   ENDDO
ENDDO
!
! Initialize matrices needed for probes
CALL initialize_probes_matrices(iprobes,nprobes)
!
! useful mathematical functions for de-aliasing
fact(1)=1.d0
DO j=2,mHOS
   fact(j)=fact(j-1)*j
ENDDO
!
!     __________________________________________________________________
!
!% INITIAL TIME DERIVATIVES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     __________________________________________________________________
!
IF ((igeom==1.OR.igeom==2).AND. &
	(icase==2.OR.icase==3.OR.icase==31.OR.icase==32.OR.icase==33.OR.icase==4.OR.icase==41.OR.icase==42)) THEN
   isubsteprk = 0
   !
   a1st    = 0.0_rp
   eta1st  = 0.0_rp
   a1st_add = 0.0_rp
   !
   a2nd    = 0.0_rp
   eta2nd  = 0.0_rp
   a2nd_add  = 0.0_rp
   !
   a3rd_add  = 0.0_rp
   ! Starting from rest with wavemaker
   deta = 0.0_rp
   !
   ! Initialize wavemaker position... (is useful if no time ramp and initial position /= 0)
   ! FIXME: should we initialize a2nd_add, a3rd_add in same way?
   CALL ramp
   CALL pos_1st
	!	RHS of the wavemaker boundary condition
	do i3=1,n3
	   do i2=1,n2
		  a1st_add(i3,i2) = -alpha*dpos1stdt(i3,i2) - alphat*pos1st(i3,i2)
	   end do
	end do
	!
	! Matching surface : polynom to ensure the C3 character of additional quantities + limited amplitude
	CALL match_surface_wmk(a1st_add(1:n3_add,1:n2),n2,n3,n3_add,l_add,delz)
	!
	!   wavemaker additional modes resolution
	a1st_add   = space_2_Fourier_add(a1st_add, 'cos','cos')
	!
	IF (coeffilt(3) < 1.0_rp) THEN
	   CALL filtering_z_modes(a1st_add,coeffilt(3))
	END IF
	!
	do i3=1,n3_add
	   do i2=1,n2
		  a1st_add(i3,i2) = a1st_add(i3,i2) / k_add_thk_add(i3,i2)
	   end do
	end do
	!
   IF(i_wmk.EQ.3) THEN
      !     __________________________________________________________________
      !
      !%  Resolution HOS_SWEET_1 : wavemaker 1st order
      !     __________________________________________________________________
      !
      CALL solveHOS_SWEET_1(a1st,da1st,eta1st,deta1st,a1st_add,da1st_add,a_add_ramp,isubsteprk)
      !
      !     __________________________________________________________________
      !
      !%  Resolution HOS_SWEET_2 : wavemaker 2nd order
      !     __________________________________________________________________
      !
      CALL solveHOS_SWEET_2(a2nd,da2nd,eta2nd,deta2nd,eta1st,a2nd_add,da2nd_add,isubsteprk)
      !
      !     __________________________________________________________________
      !
      !%  Resolution HOS_SWEET_3 : wavemaker 3rd order
      !     __________________________________________________________________
      !
      CALL solveHOS_SWEET_3(a1st,da1st,a2nd,da2nd,a1st_add,da1st_add,a2nd_add,da2nd_add,a3rd_add,da3rd_add,isubsteprk)
      !
   ELSEIF(i_wmk.EQ.2) THEN
      !
      !     __________________________________________________________________
      !
      !%  Resolution HOS_SWEET_1 : wavemaker 1st order
      !     __________________________________________________________________
      !
      CALL solveHOS_SWEET_1(a1st,da1st,eta1st,deta1st,a1st_add,da1st_add,a_add_ramp,isubsteprk)
      !
      !     __________________________________________________________________
      !
      !%  Resolution HOS_SWEET_2 : wavemaker 2nd order
      !     __________________________________________________________________
      !
      CALL solveHOS_SWEET_2(a2nd,da2nd,eta2nd,deta2nd,eta1st,a2nd_add,da2nd_add,isubsteprk)
      !
      da3rd_add = 0.0_rp
      !
   ELSEIF(i_wmk.EQ.1) THEN
      !
      !     __________________________________________________________________
      !
      !%  Resolution HOS_SWEET_1 : wavemaker 1st order
      !     __________________________________________________________________
      !
      CALL solveHOS_SWEET_1(a1st,da1st,eta1st,deta1st,a1st_add,da1st_add,a_add_ramp,isubsteprk)
      !
      da2nd_add = 0.0_rp  
      da3rd_add = 0.0_rp 
      !
   ENDIF
   !
   ! a_add_ramp has to be evaluated since in solveHOS_SWEET_1 call this is taken as a temp variable
   a_add_ramp = a1st_add
   !
   ! Initialize additional potential on FS
   CALL phiadd_SL(a1st_add+a2nd_add+a3rd_add,da1st_add+da2nd_add+da3rd_add, &
   	phi_add,phix_add,phiy_add,phiz_add,phit_add, eta)
   !
ELSE
   isubsteprk = 0
   !
   a1st_add   = 0.0_rp
   a2nd_add   = 0.0_rp
   da2nd_add = 0.0_rp  
   a3rd_add   = 0.0_rp
   da3rd_add = 0.0_rp
   !
   call solveHOS(phis,dphis_rk,eta,deta,a1st_add,da1st_add,isubsteprk,time)
   !
   a_add_ramp   = a1st_add
END IF
!
END SUBROUTINE initialize
!
!
!
SUBROUTINE initialize_probes_matrices(iprobes,nprobes)
!
IMPLICIT NONE
!
! Input variables
INTEGER, INTENT(IN) :: iprobes, nprobes
!
! Local variables
INTEGER  :: i1, i2, i3, j
REAL(RP) :: expon3, expon3p1, expon1, expon2, expon12, csh_add, sh_add, k_add_x, k_add_x1, k_add_x2, k_add_X_max
REAL(RP), DIMENSION(m3_add) :: cs_add, s_add, k_add_s_add
!    ___________________________________________________________________
!
!% WAVE PROBES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     __________________________________________________________________
!
IF (iprobes.eq.1) THEN
   DO j=1,nprobes
      DO i1=1,n1
         DO i2=1,n2
             mat_coscos(i1,i2,j) = cos(kx(i1) * xprobe(j)) * cos(ky(i2) * yprobe(j))
         ENDDO
      ENDDO
   ENDDO     
ENDIF
!    ___________________________________________________________________
!
!% PRESSURE PROBES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     __________________________________________________________________
!
!
! Matrices needed for pressure gradients calculation
k_add_X_max = 700.0_rp
!
IF (iprobes.EQ.2) THEN
   DO i3 = 1, n3_add
      cs_add(i3)   = COS(kx_add(i3))
      s_add(i3)    = SIN(kx_add(i3))
      k_add_s_add(i3) = kx_add(i3) * SIN(kx_add(i3))
      !
      do i2 = 1, n2
         k_add_x = k_add_2(i3,i2) * xlen
         IF (k_add_x <= k_add_X_max) THEN
            expon3 = dexp(-k_add_2(i3,i2) * xlen)
         ELSE
            expon3 = 0.0_rp  
         END IF
         !
         expon3p1  = expon3 + 1.d0
         DO j = 1, nprobes							   ! on the free surface
            k_add_x1 = k_add(i3,i2) * xprobe(j)
            IF (k_add_x1 <= k_add_X_max) THEN
               expon1    = dexp(-k_add_x1) / expon3p1
            ELSE
               expon1 = 0.0_rp  
            END IF
            !
            k_add_x2 = k_add_2(i3,i2) * (xlen - xprobe(j))
            IF (k_add_x2 <= k_add_X_max) THEN
               expon2    = dexp(- k_add_x2)
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
            csh_add_x_probe(i3,i2,j)     = csh_add
            k_add_sh_add_x_probe(i3,i2,j)   = sh_add * k_add(i3,i2)
            kycsh_add_x_probe(i3,i2,j)   = csh_add * ky(i2)
            kx_add_csh_add_x_probe(i3,i2,j) = csh_add * kx_add(i3)
         ENDDO
      ENDDO
   ENDDO
ENDIF
!
IF (iprobes.eq.2) THEN
   DO j=1, nprobes
	 DO i1=1,n1
		DO i2=1,n2
		   mat_coscos(i1,i2,j)     = cos(kx(i1) * xprobe(j)) * cos(ky(i2) * yprobe(j))
		   mat_sincos(i1,i2,j)     = sin(kx(i1) * xprobe(j)) * cos(ky(i2) * yprobe(j))
		   mat_cossin(i1,i2,j)     = cos(kx(i1) * xprobe(j)) * sin(ky(i2) * yprobe(j))
		   mat_sinsin(i1,i2,j)     = sin(kx(i1) * xprobe(j)) * sin(ky(i2) * yprobe(j))
		ENDDO
	 ENDDO
   ENDDO
ENDIF

END SUBROUTINE initialize_probes_matrices
!
!
!
SUBROUTINE check_range(n,name)
!
! This subroutine check the values of parameters p1/p2 or M
! FIXME: With FFTW, maybe not necessary: test the loss of efficiency for large prime numbers
!
IMPLICIT NONE
!
INTEGER :: n
CHARACTER(LEN=2) :: name
!
SELECT CASE (n)
   CASE (1,2,3,4,5,7,8,9,11,14,15,17,19,23,29)
      WRITE(*,'(A,A,A)') 'initiate_parameters: value of ',name,' in correct range'
   CASE (6,10,12,13,16,18,20,21,22,24,25,26,27,28)
      WRITE(*,'(A,A,A)') 'initiate_parameters: forbidden value of ',name,'. cf readme.txt for instructions'
      STOP
   CASE DEFAULT ! n above 30
      WRITE(*,'(A,A,A)') 'initiate_parameters: forbidden value of ',name,'. cf readme.txt for instructions'
      STOP
END SELECT
!
END SUBROUTINE check_range
!
!
!
!
!
! start check_mods ****************************************************
!           
!     ===========================
subroutine check_mods(freq, theta)
!     ===========================
!
IMPLICIT NONE
!
REAL(RP)           :: freq
REAL(RP), OPTIONAL :: theta
REAL(RP) :: ti, tf
REAL(RP) :: theta1, prec, k_1, k_2
INTEGER  :: nt
!
!	CPU times inlet
!
if (iCPUtime.eq.1) then
   print*,'entering subroutine check_mods'
   call CPU_TIME(ti)
endif
!
IF (PRESENT(theta)) THEN
   theta1 = theta
ELSE
   theta1 = 0.0_rp
END IF
!
prec=1.d-12
!
!	generated wavelength and wave number calculation
!
k_1 = wave_number_adim_r(freq, prec)
!
!	free wave wavelength and wave number calculation
!
k_2 = wave_number_adim_r(2.d0*freq, prec)
!
!	first order checking
!
WRITE(*,*) 'At first order'
nt = FLOOR(k_1*xlen*DCOS(theta1)/PI)
IF (nt.GT.n1) WRITE(*,*) 'WARNING'
WRITE(*,'(A,I4,A,I4)') 'required x-pic mode ',nt,' enabled',n1
nt = FLOOR(k_1*ylen*DSIN(theta1)/PI)
IF (nt.GT.n2) WRITE(*,*) 'WARNING'
WRITE(*,'(A,I4,A,I4)') 'required y-pic mode ',nt,' enabled',n2
!
!	second order checking
!
WRITE(*,*) 'At second order'
nt = FLOOR(k_2*xlen*DCOS(theta1)/PI)
IF (nt.GT.n1) WRITE(*,*) 'WARNING'
WRITE(*,'(A,I4,A,I4)') 'required x-pic mode ',nt,' enabled',n1
nt = FLOOR(k_2*ylen*DSIN(theta1)/PI)
IF (nt.GT.n2) WRITE(*,*) 'WARNING'
WRITE(*,'(A,I4,A,I4)') 'required y-pic mode ',nt,' enabled',n2
!      
!	CPU times outlet
!
if (iCPUtime.eq.1) then
   call CPU_TIME(tf)
   write(*,910)'quitting subroutine check_mods, total CPU time: ',tf-ti,'s'
endif
!
910 format(a,1ES11.4,a)
!
return 
!
end subroutine check_mods
!
! end check_mods ******************************************************
!
!
! start check_slope ******************************************************
!
!     ===========================
      SUBROUTINE check_slope
!     ===========================
!
IMPLICIT NONE
!
REAL(RP) :: slopex_max, slopey_max, threshold
!
threshold = 8.0_rp
slopex_max = MAXVAL(ABS(etax(1:nd1,1:nd2)))
slopey_max = MAXVAL(ABS(etay(1:nd1,1:nd2)))
!
IF (slopex_max > threshold) THEN
   print*,'etax',slopex_max,MAXLOC(ABS(etax(1:nd1,1:nd2)))
   WRITE(*,'(A)') '******************** Slope of FS elevation too high ******************** '
   WRITE(*,'(A)') 'Trying to filter eta and phis...'
   call filtering_x(phis,0.7_rp)
   call filtering_x(eta,0.7_rp)
END IF
!
IF (slopey_max > threshold) THEN
   print*,'etay',slopey_max,MAXLOC(ABS(etay(1:nd1,1:nd2)))
   WRITE(*,'(A)') '******************** Slope of FS elevation too high ******************** '
   WRITE(*,'(A)') 'Trying to filter eta and phis...'
   IF (n2 /= 1) THEN
      call filtering_y(phis,0.7_rp)
      call filtering_y(eta,0.7_rp)
   END IF
END IF
!
END SUBROUTINE check_slope
!
! end check_slope ******************************************************
!
!
! beginning absorb ************************************************
!
!     =================     
subroutine absorb
!     =================     
!
IMPLICIT NONE
!
INTEGER  :: i1,i2
REAL(RP) :: x0, x1, x2, x3, y1, y2, y3, y4
LOGICAL  :: front_abs, back_abs
!
! parabolic damping coefficient      
!
x0 = xlen * xabsf
x1 = 0.d0
!
IF (ABS(xabsb) > tiny) THEN
   back_abs = .TRUE.
   x1 = xlen * xabsb
ELSE
   back_abs= .FALSE.
END IF
!
IF (ABS(xabsf-1.0_rp) > tiny) THEN
   front_abs = .TRUE.
   x2 = xlen - x0
ELSE
   front_abs = .FALSE.
END IF
!
!  B600	specificities included
if (iwidth.eq.2) then ! B600: width limited wavemaker
   y1 = yr * (1.0_rp - yabsr)
   y2 = yr + ywmk + yl * yabsl
   y3 = ylen - y2
else ! full width wavemaker
   y1 = 0.0_rp
   y2 = ylen
   y3 = 1.0_rp
end if
!
do i1=1,n1
   DO i2=1,n2
      if ((x(i1,i2) >= x0) .AND. front_abs) then
         x3     = (x(i1,i2) - x0) / x2
         nu(i1,i2)   = coeffabsf * x3 * x3  ! * rampy(..)
         dnux(i1,i2) = coeffabsf  * 2 * x3
         dnuy(i1,i2) = 0.0_rp
      else if ((x(i1,i2) <= x1) .AND. back_abs) then
         x3     = (x(i1,i2) - x1) / x1
         nu(i1,i2)   = coeffabsb * x3 * x3  ! * rampy(..)
         dnux(i1,i2) = coeffabsf * 2 * x3
         dnuy(i1,i2) = 0.0_rp
      else 
         nu(i1,i2)   = 0.0_rp
         dnux(i1,i2) = 0.0_rp
         dnuy(i1,i2) = 0.0_rp
      endif
      if (y(i1,i2).lt.y1) then ! right beach
         y4     = (y1 - y(1,i2)) / y1
         nu(i1,i2)   = nu(i1,i2) + coeffabsr * y4 * y4 * (3.0_rp - 2.0_rp * y4)
         dnuy(i1,i2) = dnuy(i1,i2)-6.0_rp*coeffabsr*(1.0_rp-y4)*y4/y1
      else if (y(i1,i2).gt.y2) then ! left beach
         y4     = (y(i1,i2) - y2) / y3
         nu(i1,i2)   = nu(i1,i2) + coeffabsl * y4 * y4 * (3.0_rp - 2.0_rp * y4)
         dnuy(i1,i2) = dnuy(i1,i2)+6.0_rp*coeffabsl*(1.0_rp-y4)*y4/y3
      end if
   end do
ENDDO
!
return
!
END SUBROUTINE absorb
!
! end absorb ****************************************************
! 
!
END MODULE initial_condition