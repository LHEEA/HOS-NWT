MODULE wavemaker
!
! This module defines necessary subroutines for wavemaker movement in HOS-NWT
! This includes the definition of geometry, first and second order movements
! Free waves correction of wavemaker motion is also included
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
USE variables
!
REAL(RP), ALLOCATABLE, DIMENSION(:)  :: freq
!
CONTAINS
!
!     ==================
subroutine wmk_ini
!     ==================
!
IMPLICIT NONE
!
COMPLEX(CP) :: pos_tmp
REAL(RP)    :: t_i,t_f
REAL(RP)    :: df,y1,dy,ylimr,yliml,y_t,y2,x_med
INTEGER     :: j1,i1, n_wave
INTEGER     :: harm, n_max, n_indiv, n_harm, n_harm_max
REAL(RP), ALLOCATABLE, DIMENSION(:,:)  :: decomp
COMPLEX(CP), ALLOCATABLE, DIMENSION(:) :: TF
REAL(RP), ALLOCATABLE, DIMENSION(:)    :: k, f
!
!	CPU times inlet
if(iCPUtime.eq.1) then
   print*,'entering subroutine wmk_ini'
   call CPU_TIME(t_i)
endif
!
if ( igeom.ne.1 .AND. igeom.ne.2 ) then
   WRITE(*,'(A,I2)') 'Unknown wavemaker geometry, igeom=',igeom
   STOP
end if
!
WRITE(*,'(A)') 'Wave first order'
!
!	wavemaker initialization
if ( igeom.eq.1 ) then
   df=0.d0
else if ( igeom.eq.2 ) then
   df=1.d0/(1.d0-d_hinge)
else
   WRITE(*,'(A,I2)') 'Unknown wavemaker geometry, igeom=',igeom
end if
!
!  ramp against y for B600 geometry
if (iwidth.eq.2) then
   y1=yr+ywmk
   dy=y_rmp/100.0d0*ywmk
   ylimr=yr+dy
   yliml=y1-dy
   do j1=1,n2
      rampy(j1)=1.d0
      drampy(j1)=0.d0
      y_t=cy(j1)
      if ((y_t.ge.yr).and.(y_t.le.ylimr)) then
         y2=(y_t-yr)/dy
         rampy(j1)=y2**2*(3.d0-2.d0*y2)
         drampy(j1)=6.d0*y2*(1.d0-y2)/dy
      else if ((y_t.ge.yliml).and.(y_t.le.y1)) then
         y2=(y1-y_t)/dy
         rampy(j1)=y2**2*(3.d0-2.d0*y2)
         drampy(j1)=-6.d0*y2*(1.d0-y2)/dy
      else if ((y_t.le.yr).or.(y_t.ge.y1)) then
         rampy(j1)=0.d0
      end if
   end do
else
   rampy(1:n2)=1.d0
   drampy(1:n2)=0.d0
end if
!
! initialisation of the control law for spectrum
!
! number of frequency components between 0 and 1 Hz
!
n_max  = MAX(1,FLOOR((wav%config%clock * SQRT(g / wav%config%depth))))
n_max  = 2**wav%config%rnum / n_max
! total number of components
n_wave = size(wav%ampli)
! number of individual frequencies
n_indiv    = 0
n_harm_max = 0
DO harm = 0, 4 * n_max ! maximum frequency = 4 Hz (cf n_max def.)
   IF (ANY(wav%harmo == harm)) THEN
      n_indiv = n_indiv + 1
      n_harm = COUNT(wav%harmo == harm)
      IF (n_harm > n_harm_max) n_harm_max = n_harm
   END IF
END DO
!
allocate(pos_ini(n2,n_indiv), dposdy_ini(n2,n_indiv), freq(n_indiv))
allocate(pos_gzed(n3),dpos_gzeddz(n3))
!
do i1=1,n3
   pos_gzed(i1)    = CMPLX(gzed(cz(i1)), 0.0_rp, KIND = CP)
   dpos_gzeddz(i1) = CMPLX(dgzeddz(cz(i1)),0.0_rp, KIND = CP)
end do
!
n_indiv = 0
ALLOCATE(decomp(n_harm_max,3), TF(n_harm_max), k(n_harm_max), f(n_harm_max))
!
DO harm = 0, 4 * n_max ! maximum frequency = 4 Hz (cf n_max def.)
   IF (ANY(wav%harmo == harm)) THEN
      n_harm = COUNT(wav%harmo == harm)
      n_indiv = n_indiv + 1
      decomp(1:n_harm,1)    = PACK(wav%ampli,wav%harmo == harm)
      decomp(1:n_harm,2)    = PACK(wav%angle,wav%harmo == harm)
      decomp(1:n_harm,3)    = PACK(wav%phase,wav%harmo == harm)
      TF(1:n_harm)          = PACK(wav%TF,wav%harmo == harm)
      k(1:n_harm)           = PACK(wav%k,wav%harmo == harm)
      f(1:n_harm)           = PACK(wav%freq,wav%harmo == harm)
      freq(n_indiv)         = f(1)
      pos_ini(:,n_indiv)    = 0.0_cp
      dposdy_ini(:,n_indiv) = 0.0_cp
      !
      DO j1 = 1, n2
         DO i1 = 1, n_harm
            pos_tmp                = TF(i1) * COS(decomp(i1,2)) * decomp(i1,1) * &
                                      EXP(i * (- k(i1) * SIN(decomp(i1,2)) * cy(j1) + &
                                      decomp(i1,3) ) )
            pos_ini(j1,n_indiv)    = pos_ini(j1,n_indiv) + rampy(j1) * pos_tmp
            dposdy_ini(j1,n_indiv) = dposdy_ini(j1,n_indiv) + (drampy(j1) -&
                                          i * k(i1)*sin(decomp(i1,2)) * rampy(j1)) * pos_tmp
         END DO
      END DO
   END IF
END DO
DEALLOCATE(decomp, TF, k, f)
!
!   calculation of first-order wavemaker condition polynom
IF (l_add.NE.1) THEN
   x_med         = float(l_add)-1.0d0
   c_poly_add(1) = (1.d0+df*x_med)/(2.d0*x_med**3)
   c_poly_add(2) = -(3.d0+df*x_med)/(2.d0*x_med)
END IF
!
!   initialization for 2nd-order free waves correction
if ((icase.eq.2).AND.(ifree.eq.1)) call free_ini
!
!	CPU times outlet
if (iCPUtime.eq.1) then
   call CPU_TIME(t_f)
   write(*,910)'quitting subroutine wmk_ini, total CPU time: ',t_f-t_i,'s'
endif
910 format(a,1ES11.4,a)
!
end subroutine wmk_ini
!
! end wmk_ini *********************************************************
!
!
!
! start gzed **********************************************************
!
!     ================
FUNCTION gzed(z)
!     ================
!
!	Vertical wavemaker shape
!
IMPLICIT NONE
!
REAL(RP) :: z, gzed, delta, U
!
! To remove discontinuity at d_hinge for hinged-flap wmk
! FIXME: can be adjusted, what is best choice?
delta = MIN(d_hinge,(1.0_rp-d_hinge)/10.0_rp)
!
if (igeom.eq.1) then
   gzed = 1.0D0
else if (igeom.eq.2) then
   IF ( z.GE.(d_hinge - 1.0D0+delta) ) THEN
      gzed = (z + 1.0D0 - d_hinge) / (1.0D0 - d_hinge)
   ELSEIF ( z.LE.(d_hinge - 1.0D0-delta) ) THEN
      gzed = 0.0D0
   ELSE ! Polynomial interpolation
      U = z-(d_hinge - 1.d0)
      gzed = delta/(2.d0*(1.0D0 - d_hinge))*(5.d0/16.d0 + (U/delta) &
           + 15.d0/16.d0*((U/delta)**2 -1.d0/3.d0*(U/delta)**4+1.d0/15.d0*(U/delta)**6))
   END IF
else
   write(*,*) 'Unknown type of wavemaker geometry in gzed'
end if
!
END FUNCTION gzed
!
! end gzed ************************************************************
!
!
!
! start dgzeddz *******************************************************
!
!     ===================
FUNCTION dgzeddz(z)
!     ===================
!
!	Vertical derivative of the wavemaker shape
!
IMPLICIT NONE
!
REAL(RP) :: z, dgzeddz, delta, U
!
! To remove discontinuity at d_hinge for hinged-flap wmk
! FIXME: can be adjusted, what is best choice?
delta = MIN(d_hinge,(1.0_rp-d_hinge)/10.0_rp) !(equal to previous delta)...
!
if (igeom.eq.1) then
   dgzeddz = 0.0D0
else if (igeom.eq.2) then
   IF ( z.GE.(d_hinge - 1.0D0+delta) ) THEN
      dgzeddz = 1.0D0  / (1.0D0 - d_hinge)
   ELSEIF ( z.LE.(d_hinge - 1.0D0-delta) ) THEN
      dgzeddz = 0.0D0
   ELSE ! Polynomial interpolation
      U = z-(d_hinge - 1.d0)
      dgzeddz = 1.0D0  / (1.0D0 - d_hinge)*(1.d0/2.d0+15.d0/16.d0*((U/delta)-2.d0/3.d0*(U/delta)**3+1.d0/5.d0*(U/delta)**5))
   END IF
else
   write(*,*) 'Unknown type of wavemaker geometry in dgzeddz'
end if
!
END FUNCTION dgzeddz
!
! end dgzeddz *********************************************************
!
!
!
! start pos_1st *******************************************************
!
!     ==================
SUBROUTINE pos_1st
!     ==================
!
IMPLICIT NONE
!
! LOCAL VARIABLES
COMPLEX(CP) :: resu, eiwt, tmp_c
INTEGER     :: n_w, i2, i3
REAL(RP)    :: w
REAL(RP)    :: ti, tf
!
!	CPU times inlet
if (iCPUtime.eq.1) then
   print*,'entering subroutine pos_1st'
   call CPU_TIME(ti)
endif
!
!	 1st-order wavemaker position
pos1st      = 0.d0
dpos1stdt   = 0.d0
d2pos1stdt2 = 0.d0
do n_w = 1, size(freq)
   w     = TWOPI * freq(n_w)
   eiwt  = EXP(i * w * time)
   tmp_c = i * w
   do i2=1,n2
      do i3=1,n3
         resu               = pos_ini(i2,n_w) * pos_gzed(i3)
         resu               = resu * eiwt
         pos1st(i3,i2)      = pos1st(i3,i2) + DREAL(resu)
         resu               = resu * tmp_c
         dpos1stdt(i3,i2)   = dpos1stdt(i3,i2) + DREAL(resu)
         resu               = resu * tmp_c
         d2pos1stdt2(i3,i2) = d2pos1stdt2(i3,i2) + DREAL(resu)
      end do
   end do
end do
!
!	CPU times outlet
if (iCPUtime.eq.1) then
   call CPU_TIME(tf)
   write(*,910)'quitting subroutine pos_1st, total CPU time: ',tf-ti,'s'
endif
910 format(a,1ES11.4,a)
!
END subroutine pos_1st
!
! end pos_1st *********************************************************
!
!
!
! start dposdy_1st ****************************************************
!
!     =====================
SUBROUTINE dposdy_1st
!     =====================
!
IMPLICIT NONE
!
! LOCAL VARIABLES
COMPLEX(CP) :: resu, eiwt, tmp_c
INTEGER     :: n_w, i2,i3
REAL(RP)    :: w
REAL(RP)    :: ti, tf
!
!	CPU times inlet
if (iCPUtime.eq.1) then
   print*,'entering subroutine dposdy_1st'
   call CPU_TIME(ti)
endif
!
!	y-derivative of 1st-order wavemaker position
dpos1stdy    = 0.d0
d2pos1stdydt = 0.d0
do n_w = 1, size(freq)
   w = TWOPI * freq(n_w)
   eiwt  = EXP(i * w * time)
   tmp_c = DCMPLX(0.0d0, w)
   do i2=1,n2
      do i3=1,n3
         resu                = dposdy_ini(i2,n_w) * pos_gzed(i3)
         resu                = resu * eiwt
         dpos1stdy(i3,i2)    = dpos1stdy(i3,i2) + DREAL(resu)
         resu                = resu * tmp_c
         d2pos1stdydt(i3,i2) = d2pos1stdydt(i3,i2) + DREAL(resu)
      end do
   end do
end do
!
!	CPU times outlet
if (iCPUtime.eq.1) then
   call CPU_TIME(tf)
   write(*,910)'quitting subroutine dposdy_1st, total CPU time: ',tf-ti,'s'
endif
910 format(a,1ES11.4,a)
!
END subroutine dposdy_1st
!
! end dposdy_1st *****************************************************
!
!
!
! start dposdz_1st ****************************************************
!
!     =====================
SUBROUTINE dposdz_1st
!     =====================
!
IMPLICIT NONE
!
! LOCAL VARIABLES
REAL(RP)    :: ti, tf
COMPLEX(CP) ::	resu, eiwt, tmp_c
INTEGER     ::  n_w, i2, i3
REAL(RP)    :: w
!
!	CPU times inlet
if(iCPUtime.eq.1) then
   print*,'entering subroutine dposdz_1st'
   call CPU_TIME(ti)
endif
!
!	z-derivative of first-order wavemaker position
dpos1stdz    = 0.d0
d2pos1stdzdt = 0.d0
do n_w = 1, size(freq)
   w = TWOPI * freq(n_w)
   eiwt  = EXP(i * w * time)
   tmp_c = DCMPLX(0.0d0, w)
   do i2=1,n2
      do i3=1,n3
         resu                = pos_ini(i2,n_w)*dpos_gzeddz(i3)
         resu                = resu * eiwt
         dpos1stdz(i3,i2)    = dpos1stdz(i3,i2) + DREAL(resu)
         resu                = resu * tmp_c
         d2pos1stdzdt(i3,i2) = d2pos1stdzdt(i3,i2) + DREAL(resu)
      end do
   end do
end do
!
!	CPU times outlet
if(iCPUtime.eq.1) then
   call CPU_TIME(tf)
   write(*,910)'quitting subroutine dposdz_1st, total CPU time: ',tf-ti,'s'
endif
910 format(a,1ES11.4,a)
!
END subroutine dposdz_1st
!
! end dposdz_1st ******************************************************
!
!
!
! start pos_2nd *******************************************************
!     ==================
SUBROUTINE pos_2nd
!     ==================
!
IMPLICIT NONE
!
! LOCAL VARIABLES
COMPLEX(CP) :: resu, eiwt, tmp_c
REAL(RP)    :: w
INTEGER     :: i2, i3
REAL(RP)    :: ti, tf
!
!	CPU times inlet
if (iCPUtime.eq.1) then
   print*,'entering subroutine pos_2nd'
   call CPU_TIME(ti)
endif
!
! 2nd-order wavemaker postion
w     = TWOPI * 2.0_rp * freq(1)
eiwt  = EXP(i * w * time)
tmp_c = i * w
DO i2=1,n2
   DO i3=1,n3
      resu               = pos_free_ini(i2) * pos_gzed(i3)
      resu               = resu * eiwt
      pos2nd(i3,i2)      = REAL(resu)
      resu               = resu * tmp_c
      dpos2nddt(i3,i2)   = REAL(resu)
      resu               = resu * tmp_c
      d2pos2nddt2(i3,i2) = REAL(resu)
   END DO
END DO
!
!	CPU times outlet
if (iCPUtime.eq.1) then
   call CPU_TIME(tf)
   write(*,910)'quitting subroutine pos_2nd, total CPU time: ',tf-ti,'s'
endif
910 format(a,1ES11.4,a)
!
END SUBROUTINE pos_2nd
!
! end pos_2nd *********************************************************
!
!
!
! start dposdy_1st ****************************************************
!
!     =====================
SUBROUTINE dposdy_2nd
!     =====================
!
IMPLICIT NONE
!
! LOCAL VARIABLES
COMPLEX(CP) :: resu, eiwt, tmp_c
INTEGER     :: n_w, i2, i3
REAL(RP)    :: w
REAL(RP)    :: ti, tf
!
!	CPU times inlet
if (iCPUtime.eq.1) then
   print*,'entering subroutine dposdy_2nd'
   call CPU_TIME(ti)
endif
!
!	y-derivative of 2nd-order wavemaker position
dpos2nddy    = 0.d0
d2pos2nddydt = 0.d0
do n_w = 1, size(freq)
   w = TWOPI * freq(n_w)
   eiwt  = EXP(i * w * time)
   tmp_c = DCMPLX(0.0d0, w)
   do i2=1,n2
      do i3=1,n3
         resu                = dposdy2nd_ini(i2) * pos_gzed(i3)
         resu                = resu * eiwt
         dpos2nddy(i3,i2)    = dpos2nddy(i3,i2) + DREAL(resu)
         resu                = resu * tmp_c
         d2pos2nddydt(i3,i2) = d2pos2nddydt(i3,i2) + DREAL(resu)
      end do
   end do
end do
!
!	CPU times outlet
if (iCPUtime.eq.1) then
   call CPU_TIME(tf)
   write(*,910)'quitting subroutine dposdy_2nd, total CPU time: ',tf-ti,'s'
endif
910 format(a,1ES11.4,a)
!
END subroutine dposdy_2nd
!
! end dposdy_2nd *****************************************************
!
!
!
! start dposdz_2nd ****************************************************
!
!     =====================
SUBROUTINE dposdz_2nd
!     =====================
!
IMPLICIT NONE
!
! LOCAL VARIABLES
REAL(RP)    :: ti, tf
COMPLEX(CP) :: resu, eiwt, tmp_c
INTEGER     :: n_w, i2 ,i3
REAL(RP)    :: w
!
!	CPU times inlet
if(iCPUtime.eq.1) then
   print*,'entering subroutine dposdz_1st'
   call CPU_TIME(ti)
endif
!
!	z-derivative of first-order wavemaker position
dpos2nddz    = 0.d0
d2pos2nddzdt = 0.d0
do n_w = 1, size(freq)
   w = TWOPI * freq(n_w)
   eiwt  = EXP(i * w * time)
   tmp_c = DCMPLX(0.0d0, w)
   do i2=1,n2
      do i3=1,n3
         resu                = pos_free_ini(i2)*dpos_gzeddz(i3)
         resu                = resu * eiwt
         dpos2nddz(i3,i2)    = dpos2nddz(i3,i2) + DREAL(resu)
         resu                = resu * tmp_c
         d2pos2nddzdt(i3,i2) = d2pos2nddzdt(i3,i2) + DREAL(resu)
      end do
   end do
end do
!
!	CPU times outlet
if(iCPUtime.eq.1) then
   call CPU_TIME(tf)
   write(*,910)'quitting subroutine dposdz_1st, total CPU time: ',tf-ti,'s'
endif
910 format(a,1ES11.4,a)
!
END subroutine dposdz_2nd
!
! end dposdz_2nd ******************************************************
!
!
!
! start free_ini ******************************************************
!     ======================================
SUBROUTINE free_ini !(wave)
!     ======================================
!
USE definition2
!
IMPLICIT NONE
!
! LOCAL VARIABLES
INTEGER       :: j,N_free,nbou
REAL(RP)      :: k_free
REAL(RP)      :: sumprec = 1.0E-4_rp
REAL(RP)      :: ti,tf
!
! Old definition2 stuff
TYPE(size_basin) :: basin
TYPE(wave_12)    :: wave12
INTEGER          :: Nf
REAL(RP)         :: prec
COMPLEX(CP), DIMENSION(:), ALLOCATABLE :: coeff
!
!	CPU times inlet
if (iCPUtime.eq.1) then
   print*,'entering subroutine free_ini'
   call CPU_TIME(ti)
endif
!
WRITE(*,'(A)') 'Wave second order'
!
! Compatibility checkings
!     already done before calling free_ini
!
! To be done...
!
! calculate a_on^free up to a given N_2 (given as in Dalrymple method)
!  rmk : i)   snake motion is imposed (1st order) without cutting for large n>N_1
!               so 2nd order correction must be calculated with n to inf.
!        ii)  Dalrymple's motion cutted at n=N_1
!               so 2nd order correction must be calculated with n <= N_1 !
!        iii) maybe calculated above n = N_2 to remove evanescent free modes as well.
!
! check it for snake principle: does such a N_2 exists ?
!
! separate each a_on^free in two (+ or minus direction)
!
! impose a snake motion for each part at 2w frequency.
!
k_free = wave_number_adim_r(2.0_rp * wav%freq(1),1.0E-15_rp)
IF (n2 == 1) THEN
   N_free = 0
ELSE
   N_free = FLOOR(k_free * wav%config%Ly / PI)
END IF
!
ALLOCATE(pos_free_ini(n2), dposdy2nd_ini(n2))
!
! Old definition2 stuff
if ( igeom.eq.1 ) then
   call initialize_basin(basin,1.0d0,0.0d0,ylen,'piston',d_hinge,.FALSE.,48,yr,yl,'',0)
else if ( igeom.eq.2 ) then
   call initialize_basin(basin,1.0d0,0.0d0,ylen,'hinged',d_hinge,.FALSE.,48,yr,yl,'',0)
else
   WRITE(*,'(A,I2)') 'Unknown wavemaker geometry, igeom=',igeom
end if
!
WRITE(*,'(A)') 'Wave first order'
prec  = 1.d-12
!
if (n2.ne.1) then
   call initialize_wave1(wave12, basin, nu_mono, theta_mono, amp_mono, ph_mono, prec, '3D')
else
   call initialize_wave1(wave12, basin, nu_mono, theta_mono, amp_mono, ph_mono, prec, '2D')
end if
!
IF (ibat == 2) THEN
   CALL select_wavemaker(wave12, 'snake')
ELSE IF (ibat == 3) THEN
   CALL select_wavemaker(wave12, 'Dalrymple', xd_mono)
ELSE IF (ibat == 4) THEN
   CALL select_wavemaker(wave12, 'B600')
END IF
!
WRITE(*,'(A)') 'Wave second order'
!
CALL initialize_wave2(wave12,'free',sumprec)
if (n2.ne.1 .AND. ABS(wave12%angle) > tiny) then
   Nf = wave12%N(2)
else
   Nf = 0
end if
!
ALLOCATE(coeff(0:Nf))
OPEN(222,FILE='coeff_sweet.dat')
DO nbou = 0, Nf
   coeff(nbou) = al_mn(wave12,0,nbou,'direct','direct','direct')
   WRITE(222,*) nbou, coeff(nbou)
END DO
CLOSE(222)
!
! initialisation of the control law
DO j = 1, n2
   pos_free_ini(j)  = DCMPLX(0.0D0,0.0D0)
   dposdy2nd_ini(j) = DCMPLX(0.0D0,0.0D0)
   DO nbou = 0, Nf
      pos_free_ini(j) = pos_free_ini(j) + &
           DCMPLX( rampy(j) * DCOS(mu_n(wave12,nbou)*cy(j)), 0.0_rp) &
           * coeff(nbou) &
           * gamma_mn(wave12,0,nbou)
      dposdy2nd_ini(j) =  dposdy2nd_ini(j) + &
           DCMPLX( drampy(j) * DCOS(mu_n(wave12,nbou)*cy(j)), 0.0_rp) &
           * coeff(nbou) &
           * gamma_mn(wave12,0,nbou) + &
           DCMPLX( rampy(j) * (-mu_n(wave12,nbou)*DSIN(mu_n(wave12,nbou)*cy(j))), 0.0_rp) &
           * coeff(nbou) &
           * gamma_mn(wave12,0,nbou)
   END DO
END DO
!
pos_free_ini(:n2) = pos_free_ini(:n2) * wave12%ampli * wave12%ampli * wave12%TF_f / beta_m(wave12,0)
pos_free_ini(:n2) = - pos_free_ini(:n2)
!
dposdy2nd_ini(:n2) = dposdy2nd_ini(:n2) * wave12%ampli * wave12%ampli * wave12%TF_f / beta_m(wave12,0)
dposdy2nd_ini(:n2) = - dposdy2nd_ini(:n2)
!
!	CPU times outlet
if (iCPUtime.eq.1) then
   call CPU_TIME(tf)
   write(*,910)'quitting subroutine free_ini, total CPU time: ',tf-ti,'s'
endif
910 format(a,1ES11.4,a)
!
END subroutine free_ini
!
! end free_ini *****************************************************
!
END MODULE wavemaker
