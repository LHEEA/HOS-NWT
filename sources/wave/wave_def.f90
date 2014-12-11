MODULE wave_def
!
! This module defines the wave objects and the corresponding actions associated
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
USE config_def
USE linear_wave
!
IMPLICIT NONE
!
TYPE wave
   TYPE(config)                       :: config
   INTEGER, POINTER, DIMENSION(:)     :: harmo => NULL( )
   REAL(RP), POINTER, DIMENSION(:)    :: freq  => NULL( )
   REAL(RP), POINTER, DIMENSION(:)    :: ampli => NULL( )
   REAL(RP), POINTER, DIMENSION(:)    :: angle => NULL( )
   REAL(RP), POINTER, DIMENSION(:)    :: phase => NULL( )
   LOGICAL                            :: isinit = .FALSE.
   COMPLEX(CP), POINTER, DIMENSION(:) :: TF => NULL( )
   REAL(RP), POINTER, DIMENSION(:)    :: k  => NULL( )
END TYPE wave
!
INTERFACE OPERATOR (+)
   MODULE PROCEDURE add_wave
END INTERFACE
!
!
!
CONTAINS
!
!
!
FUNCTION front(cfg, harmo, ampli, angle, phase)
!
IMPLICIT NONE
!
! Input variables
TYPE(config) :: cfg
INTEGER      :: harmo
REAL(RP)     :: ampli, angle, phase
TYPE(wave)   :: front
! Routine variables
INTEGER      :: n_f
!
front%config = cfg
!
n_f = 1
!
ALLOCATE(front%harmo(n_f), front%freq(n_f), front%ampli(n_f), front%angle(n_f), front%phase(n_f))
!
front%harmo(n_f) = harmo
front%freq(n_f)  = harmo * front%config%clock / 2**front%config%rnum
front%ampli(n_f) = ampli
front%angle(n_f) = angle
front%phase(n_f) = phase
!
front = wave_init(front)
!
END FUNCTION front
!
FUNCTION single(cfg, freq, ampli, angle, phase)
!
IMPLICIT NONE
!
! Input variables
TYPE(config) :: cfg
REAL(RP)     :: freq, ampli, angle, phase
TYPE(wave)   :: single
! Routine variables
INTEGER      :: n_f
!
single%config = cfg
!
n_f = 1
!
ALLOCATE(single%harmo(n_f), single%freq(n_f), single%ampli(n_f), single%angle(n_f), single%phase(n_f))
!
single%harmo(n_f) = 0
single%freq(n_f)  = freq
single%ampli(n_f) = ampli
single%angle(n_f) = angle
single%phase(n_f) = phase
!
single = wave_init(single)
!
END FUNCTION single
!
!
!
FUNCTION dalrymple(cfg, X_d, harmo, ampli, angle, phas)
!
IMPLICIT NONE
!
! Input variables
TYPE(config) :: cfg
REAL(RP)     :: X_d, ampli, angle, phas
INTEGER      :: harmo
TYPE(wave)   :: dalrymple
! Routine variables
INTEGER      :: n_f, N, iloop
REAL(RP)     :: freq, k, beta
COMPLEX(CP)  :: tmp
REAL(RP), ALLOCATABLE, DIMENSION(:)    :: mu_n, theta_n, delta_on
COMPLEX(CP), ALLOCATABLE, DIMENSION(:) :: I_n, k_n, a_on
!
dalrymple%config = cfg
!
freq    = harmo * dalrymple%config%clock / 2**(dalrymple%config%rnum)
k       = wave_number_adim_r(freq,1.0E-12_rp)
N       = FLOOR(k * dalrymple%config%Ly / PI)
ALLOCATE(mu_n(0:N), theta_n(0:N), delta_on(0:N), I_n(0:N), a_on(0:N), k_n(0:N))
mu_n    = (/(iloop,iloop=0,N)/) * pi / dalrymple%config%Ly
theta_n = ASIN(mu_n / k)
beta    = k * SIN(angle)
delta_on(1:N) = 0.0_rp
delta_on(0)   = 1.0_rp
!
IF (ABS(angle) < tiny) THEN
    I_n    = delta_on
ELSE
    tmp = - exp(-i * beta * dalrymple%config%Ly)
    DO iloop = 0,N
        tmp = - tmp
        IF (ABS(mu_n(iloop)-abs(beta)) < tiny) THEN
            I_n(iloop) = 1
        ELSE
            I_n(iloop) = 2 * i * beta * (tmp - 1) / (beta**2 - mu_n(iloop)**2) / &
                ((1 + delta_on(iloop)) * dalrymple%config%Ly)
        END IF
    END DO
END IF
k_n   = i * SQRT(k**2 - mu_n**2)
a_on  = ampli * exp(i * phas) * I_n * exp((k_n - i * k * COS(angle)) * X_d)
!
n_f = 2 * (N+1)
!
ALLOCATE(dalrymple%harmo(n_f), dalrymple%freq(n_f), dalrymple%ampli(n_f), dalrymple%angle(n_f), dalrymple%phase(n_f))
!
dalrymple%harmo = harmo
!
dalrymple%freq  = freq
!
dalrymple%ampli(1:N+1)       = abs(a_on / 2.0_rp)
dalrymple%ampli(N+2:2*(N+1)) = abs(a_on / 2.0_rp)
!
dalrymple%angle(1:N+1)       = theta_n
dalrymple%angle(N+2:2*(N+1)) = - theta_n
!
dalrymple%phase(1:N+1)       = phase(a_on)
dalrymple%phase(N+2:2*(N+1)) = phase(a_on)
!
dalrymple = wave_init(dalrymple)
!
DEALLOCATE(mu_n, theta_n, delta_on, I_n, a_on, k_n)
!
END FUNCTION dalrymple
!
!
!
FUNCTION file2wave(file_name)
!
IMPLICIT NONE
!
! Input variables
CHARACTER(LEN=*) :: file_name
TYPE(wave)       :: file2wave
! Routine variables
INTEGER          :: n_f, iloop, jloop, n_start, n_stop, n_good
REAL(RP)         :: Tscale, h
INTEGER,  ALLOCATABLE, DIMENSION(:) :: harmo
REAL(RP), ALLOCATABLE, DIMENSION(:) :: freq, ampli, angle, phas
!
! Config file
OPEN(999,FILE=TRIM(file_name)//'.cfg')
!
READ(999,'(I3)') file2wave%config%rnum
READ(999,*) file2wave%config%clock
READ(999,*) file2wave%config%depth
h      = file2wave%config%depth
Tscale = sqrt(h/g)
READ(999,*) file2wave%config%ramp
READ(999,*) file2wave%config%typ_ramp
READ(999,*) file2wave%config%cutoff_low
READ(999,*) file2wave%config%cutoff_high
READ(999,*) file2wave%config%Ly
READ(999,*) file2wave%config%typ_wmk
READ(999,*) file2wave%config%hinge

file2wave%config%clock  = file2wave%config%clock  * Tscale
file2wave%config%cutoff_low  = file2wave%config%cutoff_low * Tscale
file2wave%config%cutoff_high = file2wave%config%cutoff_high * Tscale
file2wave%config%Ly     = file2wave%config%Ly / h
file2wave%config%hinge  = file2wave%config%hinge / h
!
CLOSE(999)
!
! Frequency data file
OPEN(999,FILE=TRIM(file_name)//'.dat')
!
n_f = 0
!
DO
   READ(999,*,END=98)
   n_f = n_f + 1
   CYCLE
98 EXIT
END DO
!
WRITE(*,'(A,A,A,I5,A)') 'File ',TRIM(file_name),' contains ',n_f,' samples'
!
REWIND(999)
!
ALLOCATE(harmo(n_f), freq(n_f), ampli(n_f), angle(n_f), phas(n_f))
!
n_start = 1
n_stop  = 0
DO iloop = 1, n_f
   READ(999,*) harmo(iloop), freq(iloop), ampli(iloop), angle(iloop), phas(iloop)
   IF (freq(iloop) <= file2wave%config%cutoff_high / Tscale) THEN
      n_stop = iloop
   END IF
   IF (freq(iloop) < file2wave%config%cutoff_low / Tscale) THEN
      n_start = n_start + 1
   END IF
END DO
!
n_good = n_stop - n_start + 1
!
WRITE(*,'(A,A,A,I5,A)') 'File ',TRIM(file_name),' contains ', n_good,' good samples'
!
IF (n_good == 0) STOP 'n_good == 0'
!
ALLOCATE(file2wave%harmo(n_good), file2wave%freq(n_good), file2wave%ampli(n_good), &
file2wave%angle(n_good), file2wave%phase(n_good))
!
DO iloop = 1, n_good
   jloop = n_start + (iloop - 1)
   file2wave%harmo(iloop) = harmo(jloop)
   file2wave%freq(iloop)  = freq(jloop)
   file2wave%ampli(iloop) = ampli(jloop)
   file2wave%angle(iloop) = angle(jloop)
   file2wave%phase(iloop) = phas(jloop)
END DO
!
DEALLOCATE(harmo, freq, ampli, angle, phas)
!
CLOSE(999)
!
! Making wave data nondimensionalised
file2wave%freq  = file2wave%freq * Tscale
file2wave%ampli = file2wave%ampli / h
file2wave%angle = file2wave%angle * pi / 180.0_rp
!
! Building wavenumber and Wavemaker Transfer Function at each frequency
file2wave = wave_init(file2wave)
!
END FUNCTION file2wave
!
!
!
FUNCTION wave_init(wv)
!
TYPE(wave) :: wv, wave_init
!
INTEGER    :: n_f, iloop
REAL(RP)   :: k_tmp
REAL(RP), ALLOCATABLE, DIMENSION(:)    :: k
COMPLEX(CP), ALLOCATABLE, DIMENSION(:) :: TF
!
IF (.NOT.wv%isinit) THEN
   n_f    = SIZE(wv%harmo)
   ALLOCATE(TF(n_f), k(n_f))
   !
   IF (wv%config%typ_wmk == 'dipole') THEN
      DO iloop = 1, n_f
         wv%k(iloop)  = wave_number_adim_r(wv%freq(iloop),1.0E-15_rp)
      END DO
   ELSE
      DO iloop = 1, n_f
         k_tmp     = wave_number_adim_r(wv%freq(iloop),1.0E-15_rp)
         k(iloop)  = k_tmp
         IF (ABS(k_tmp)>tiny) THEN
            TF(iloop) = TF_adim_new(wv%config%typ_wmk, k_tmp, wv%config%hinge)
         ELSE
            TF(iloop) = 1.0_rp
         ENDIF
      END DO
   END IF
!
   WHERE (wv%freq < 1.0E-8_rp)
      k  = 0.0_rp
      !TF = CMPLX(0.0_rp,0.0_rp, KIND=CP)
      ! FIXME: To take into account possible non-zero mean (not fully tested)
	  TF = CMPLX(1.0_rp,0.0_rp, KIND=CP)
   END WHERE
   !
   wave_init = wv
   !
   ALLOCATE(wave_init%TF(n_f), wave_init%k(n_f))
   wave_init%TF     = TF
   wave_init%k      = k
   wave_init%isinit = .TRUE.
   DEALLOCATE(TF, k)
ELSE
   wave_init = wv
END IF
!
END FUNCTION wave_init
!
!
!
FUNCTION add_wave(wave1,wave2)
!
IMPLICIT NONE
!
! Input variables
TYPE(wave), INTENT(IN)   :: wave1, wave2
TYPE(wave)               :: add_wave
! Local variables
INTEGER :: n_w, n_w1, n_w2, wloop, wtmp
!
IF ( ASSOCIATED(wave1%ampli) ) THEN
   n_w1 = SIZE(wave1%ampli)
ELSE
   n_w1 = 0
END IF
!
IF ( ASSOCIATED(wave2%ampli) ) THEN
   n_w2 = SIZE(wave2%ampli)
ELSE
   n_w2 = 0
END IF
!
n_w = n_w1 + n_w2
IF (n_w /= 0) THEN
   add_wave%config = wave2%config
   ALLOCATE(add_wave%harmo(n_w), add_wave%freq(n_w), add_wave%ampli(n_w), add_wave%angle(n_w), add_wave%phase(n_w))
   ALLOCATE(add_wave%k(n_w), add_wave%TF(n_w))
   DO wloop=1, n_w1
      add_wave%harmo(wloop) = wave1%harmo(wloop)
      add_wave%freq(wloop)  = wave1%freq(wloop)
      add_wave%ampli(wloop) = wave1%ampli(wloop)
      add_wave%angle(wloop) = wave1%angle(wloop)
      add_wave%phase(wloop) = wave1%phase(wloop)
      add_wave%k(wloop)     = wave1%k(wloop)
      add_wave%TF(wloop)    = wave1%TF(wloop)
   END DO
   !
   DO wloop=n_w1+1, n_w
      wtmp = wloop - n_w1
      add_wave%harmo(wloop) = wave2%harmo(wtmp)
      add_wave%freq(wloop)  = wave2%freq(wtmp)
      add_wave%ampli(wloop) = wave2%ampli(wtmp)
      add_wave%angle(wloop) = wave2%angle(wtmp)
      add_wave%phase(wloop) = wave2%phase(wtmp)
      add_wave%k(wloop)     = wave2%k(wtmp)
      add_wave%TF(wloop)    = wave2%TF(wtmp)
   END DO
   !
   add_wave%isinit = .TRUE.
END IF
!
END FUNCTION add_wave
!
!
!
SUBROUTINE kill_wave(wav)
!
IMPLICIT NONE
!
! Input variables
TYPE(wave), INTENT(INOUT) :: wav
!
IF (ASSOCIATED(wav%harmo)) THEN
   NULLIFY(wav%harmo)
!    wav%harmo => NULL( )
END IF
!
IF (ASSOCIATED(wav%freq)) THEN
   NULLIFY(wav%freq)
!    wav%freq => NULL( )
END IF
!
IF (ASSOCIATED(wav%ampli)) THEN
   NULLIFY(wav%ampli)
!    wav%ampli => NULL( )
END IF
!
IF (ASSOCIATED(wav%angle)) THEN
   NULLIFY(wav%angle)
!    wav%angle => NULL( )
END IF
!
IF (ASSOCIATED(wav%phase)) THEN
   NULLIFY(wav%phase)
!    wav%phase => NULL( )
END IF
!
IF (ASSOCIATED(wav%TF)) THEN
   NULLIFY(wav%TF)
!    wav%TF => NULL( )
END IF
!
IF (ASSOCIATED(wav%k)) THEN
   NULLIFY(wav%k)
!    wav%k => NULL( )
END IF
!
wav%isinit = .FALSE.
!
END SUBROUTINE kill_wave
!
!
!
SUBROUTINE display(wv)
!
IMPLICIT NONE
!
TYPE(wave) :: wv
INTEGER :: n_w, n_max
!
CALL display_config(wv%config)
IF (ASSOCIATED(wv%ampli)) THEN
   n_max = 20
   n_w   = SIZE(wv%ampli)
   IF (n_w <= n_max) THEN
      WRITE(*,'(A,20(I4,X))') 'harmo = ',wv%harmo
      WRITE(*,'(A,20(ES11.4,X))') 'freq  = ',wv%freq
      WRITE(*,'(A,20(ES11.4,X))') 'ampli = ',wv%ampli
      WRITE(*,'(A,20(ES11.4,X))') 'angle = ',wv%angle
      WRITE(*,'(A,20(ES11.4,X))') 'phase = ',wv%phase
   ELSE
      WRITE(*,'(A,20(I4,X))') 'harmo = ',wv%harmo(1:20)
      WRITE(*,'(A,20(ES11.4,X))') 'freq  = ',wv%freq(1:20)
      WRITE(*,'(A,20(ES11.4,X))') 'ampli = ',wv%ampli(1:20)
      WRITE(*,'(A,20(ES11.4,X))') 'angle = ',wv%angle(1:20)
      WRITE(*,'(A,20(ES11.4,X))') 'phase = ',wv%phase(1:20)
   END IF
END IF

END SUBROUTINE display
!
!
!
SUBROUTINE wave2file(wv, filename)
!
IMPLICIT NONE
!
TYPE(wave) :: wv
CHARACTER(LEN=*) :: filename
INTEGER :: n_w, n
!
CALL display_config(wv%config)
OPEN(UNIT=111, FILE=filename)
!
IF (ASSOCIATED(wv%ampli)) THEN
   n_w   = SIZE(wv%ampli)
   DO n = 1, n_w
      WRITE(111,'(I4,A,4(ES11.4,A))') wv%harmo(n),' , ',wv%freq(n),' , ',wv%ampli(n),' , ' &
      ,wv%angle(n),' , ',wv%phase(n),' , '
   END DO
END IF
!
CLOSE(111)
!
END SUBROUTINE wave2file
!
!
!
ELEMENTAL FUNCTION phase(z)
!
! gives the phase of the complex number z (between 0 and 2*PI)
!
USE type
!
IMPLICIT NONE
! Input variable
COMPLEX(CP), INTENT(IN) :: z
! Local variables
COMPLEX(CP)             :: ztemp
REAL(RP)                :: phase, ptemp
!
IF ( ABS(z) > tiny ) THEN
   phase = 0.0_rp
   RETURN
ELSE
   ztemp = z / ABS(z)
   ptemp = ASIN( AIMAG( ztemp ) )
   IF ( REAL(ztemp) >= 0.0_rp ) THEN
      IF ( AIMAG(ztemp) >= 0.0_rp ) THEN
         phase = ptemp
         RETURN
      ELSE
         phase = TWOPI + ptemp
         RETURN
      END IF
   ELSE
      phase = PI - ptemp
      RETURN
   END IF
END IF
!
END FUNCTION phase
!
!
!
END MODULE wave_def
