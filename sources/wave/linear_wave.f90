MODULE linear_wave
!
! This module contains functions allowing the use of linear (1st order) dispersion relations
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
INTERFACE alpha_adim
   MODULE PROCEDURE alpha_adim_r, alpha_adim_v
END INTERFACE
INTERFACE TE_adim
   MODULE PROCEDURE TE_adim_r, TE_adim_v
END INTERFACE
!
CONTAINS
!
FUNCTION alpha_r(f, N, prec, h)
!
IMPLICIT NONE
REAL(RP)       :: f, prec, h
INTEGER        :: N
REAL(RP)       :: alpha_r
REAL(RP)       :: asympt
!
! Computation of alpha_N
asympt = (N-1.0_rp)*PI
IF (xtanx(asympt-PI/4.0_rp) <= - (TWOPI*f)**2 * h / g) THEN
   alpha_r = dichoto( xtanx, - (TWOPI*f)**2 * h / g, asympt-PI/4.0_rp, asympt, prec )
   alpha_r = alpha_r / h
ELSE IF (xtanx(asympt-PI/4.0_rp) > - (TWOPI*f)**2 * h / g) THEN
   alpha_r = dichoto( xtanx, - (TWOPI*f)**2 * h / g, asympt-PI/2.0_rp+prec, asympt-PI/4.0_rp, prec )
   alpha_r = alpha_r / h
END IF
!
END FUNCTION alpha_r
!
!
!
FUNCTION alpha_adim_r(f, N, prec)
!
IMPLICIT NONE
!
REAL(RP)     :: f, prec
INTEGER      :: N
REAL(RP)     :: alpha_adim_r
REAL(RP)     :: alpha_min, alpha_max, tmp, f0, Npi
!
! Computation of alpha_N
Npi       = N * PI
alpha_max = Npi
tmp       = prec
alpha_min = Npi - PI/2.0_rp + tmp
f0        = - (TWOPI*f)**2
DO WHILE(xtanx(alpha_min) > f0)
	tmp       = tmp / 2.0_rp
    WRITE(*,*) N,tmp, f0, xtanx(alpha_min)
	READ(*,*)
	alpha_min = alpha_min - tmp
END DO
alpha_adim_r = dichoto( xtanx, f0, alpha_min, alpha_max, prec )
!
END FUNCTION alpha_adim_r
!
!
!
FUNCTION alpha_adim_v(f, N, prec)
!
IMPLICIT NONE
!
REAL(RP)                     :: f, prec
INTEGER, DIMENSION(:)        :: N
REAL(RP), DIMENSION(SIZE(N)) :: alpha_adim_v
INTEGER                      :: ibou
REAL(RP)                     :: alpha_min, alpha_max, tmp, f0, Npi
!
! Computation of alpha_N
f0 = - (TWOPI*f)**2
DO ibou = 1, SIZE(N)
   !
   Npi       = N(ibou) * PI
   alpha_max = Npi
   tmp       = prec
   alpha_min = Npi - PI/2.0_rp + tmp
   DO WHILE(xtanx(alpha_min) > f0)
      tmp       = tmp / 2.0_rp
      WRITE(*,*) N(ibou),tmp, f0, xtanx(alpha_min)
      alpha_min = alpha_min - tmp
   END DO
   alpha_adim_v(ibou) = dichoto( xtanx, f0, alpha_min, alpha_max, prec )
   !
END DO
!
END FUNCTION alpha_adim_v
!
!
!
FUNCTION wave_number_r(f,h,prec)
!
IMPLICIT NONE
REAL(RP)       :: f,h,prec
REAL(RP)       :: wave_number_r
!
! Computation of k
IF ( (TWOPI*f)**2 / g * h < log(HUGE(1.0_rp))-1.0_rp ) THEN
   wave_number_r = dichoto(xthx,(TWOPI*f)**2 * h / g,0.0_rp,10000000.0_rp,prec)
   wave_number_r = wave_number_r / h
ELSE
   wave_number_r = (TWOPI*f)**2 / g
END IF
!
END FUNCTION wave_number_r
!
!
!
FUNCTION wave_number_adim_r(f,prec)
!
IMPLICIT NONE
REAL(RP)       :: f,prec
REAL(RP)       :: wave_number_adim_r
!
! Computation of k
IF ( (TWOPI*f)**2 < log(HUGE(1.0_rp))-1.0_rp ) THEN
   wave_number_adim_r = dichoto(xthx,(TWOPI*f)**2,0.0_rp,1.0e6_rp,prec)
ELSE
   WRITE(*,*) 'big frequency'
   wave_number_adim_r = (TWOPI*f)**2
END IF
!
END FUNCTION wave_number_adim_r
!
!
!
FUNCTION wave_number_v(f,h,prec)
!
IMPLICIT NONE
REAL(RP), DIMENSION(:)       :: f
REAL(RP)                     :: h,prec
REAL(RP), DIMENSION(SIZE(f)) :: wave_number_v
INTEGER                      :: ibou
!
! Computation of k
DO ibou = 1, SIZE(f)
   wave_number_v(ibou) = wave_number_r(f(ibou), h, prec)
END DO
!
END FUNCTION wave_number_v
!
!
!
FUNCTION xthx(x)
!
IMPLICIT NONE
!
REAL(RP) :: x,xthx
!
xthx = x*tanh(x)
!
END FUNCTION xthx
!
!
!
FUNCTION xtanx(x)
!
IMPLICIT NONE
!
REAL(RP) :: x,xtanx
!
xtanx = x*tan(x)
!
END FUNCTION xtanx
!
!
!
FUNCTION phase_velocity(f,h)
!
IMPLICIT NONE
REAL(RP)   :: f, h
REAL(RP)   :: phase_velocity
REAL(RP)   :: k, prec
!
prec          = 1E-7_rp
k             = wave_number_r(f, h, prec)
phase_velocity = TWOPI * f / k
!
END FUNCTION phase_velocity
!
!
!
FUNCTION group_velocity(f,h)
!
IMPLICIT NONE
REAL(RP)   :: f, h
REAL(RP)   :: group_velocity
REAL(RP)   :: k, prec
!
prec = EPSILON(f)
k    = wave_number_r(f, h, prec)
IF ( 2.0_rp*k*h < log(HUGE(1.0_rp))-1.0_rp ) THEN
   group_velocity = phase_velocity(f,h) / 2.0_rp * (1.0_rp + (2.0_rp*k*h) / SINH(2.0_rp*k*h))
ELSE
   group_velocity = phase_velocity(f,h) / 2.0_rp
END IF
!
END FUNCTION group_velocity
!
!
!
FUNCTION omega_second(f,h)
!
IMPLICIT NONE
REAL(RP)   :: f, h
REAL(RP)   :: omega_second
REAL(RP)   :: k, prec, v_g
!
prec = EPSILON(f)
k    = wave_number_r(f, h, prec)
v_g  = group_velocity(f,h)
IF ( 2.0_rp*k*h < log(HUGE(1.0_rp))-1.0_rp ) THEN
   omega_second = v_g/k*(1.0_rp-2.0_rp*k*h*TANH(k*h))-v_g**2/(TWOPI*f)+phase_velocity(f,h) / (2.0_rp*k) * (2.0_rp*k*h / &
   TANH(2.0_rp*k*h)-1.0_rp)
ELSE
   WRITE(*,*) 'not finished'
   omega_second = - v_g**2/(TWOPI*f)
END IF
!
END FUNCTION omega_second
!
!
!
FUNCTION TF_adim_r(nu, wmk, d, prec)
!
! adimensional Transfer Function, inspired from ECN report potentiel.tex
!   b = TF * a
!
IMPLICIT NONE
!
REAL(RP)            :: nu, d, prec
CHARACTER(LEN=*)    :: wmk
COMPLEX(CP)         :: TF_adim_r
REAL(RP)            :: k
!
k = wave_number_adim_r(nu, prec)
!
IF ( wmk == 'ECN' ) THEN
   TF_adim_r = k * (1.0_rp - d) * (2.0_rp * k + SINH(2.0_rp * k)) / (4.0_rp * i * SINH(k))
   TF_adim_r = TF_adim_r / (k * (1.0_rp - d) * SINH(k) + COSH(k * d) - COSH(k))
ELSE
   WRITE(*,*) 'Unknown wavemaker type for transfer function'
END IF
!
END FUNCTION TF_adim_r
!
!
!
ELEMENTAL FUNCTION TF_adim_new(typ, k, hinge)
!
! adimensional Transfer Function, inspired from ECN report potentiel.tex
!   b = TF * a
!
IMPLICIT NONE
!
! Input variables
CHARACTER(LEN=*), INTENT(IN)   :: typ
REAL(RP), INTENT(IN)           :: k
REAL(RP), INTENT(IN), OPTIONAL :: hinge
! Output variable
COMPLEX(CP)        :: TF_adim_new
! Local variables
COMPLEX(CP)        :: resu
REAL(RP)           :: k2, w2, w4, d
!
IF ( typ == 'hinged' .OR. typ == 'monofl' ) THEN
   d  = hinge
   IF ((k.LT.50).AND.(k*d.LT.50)) THEN
       resu = - i * k * (1.0_rp - d) * (2.0_rp * k + &
                    SINH(2.0_rp * k)) / (4.0_rp * SINH(k))
       resu = resu / (k * (1.0_rp - d) * SINH(k) + COSH(k * d) - COSH(k))
   ELSEIF ((k.GE.50).AND.(k*d.LT.50)) THEN
       resu = - i * k * (1.0_rp - d) / (2.0_rp * k * (1.0_rp - d) -2.0_rp + 4.0_rp*COSH(k * d)*EXP(-k))
   ELSE ! both are greater than 50
       IF (k*(1.0_rp-d).GE.20) THEN
           resu = - i * 0.5_rp
       ELSE
           resu = - i * k * (1.0_rp - d) / (2.0_rp * k * (1.0_rp - d) -2.0_rp + 2.0_rp*EXP(k*(d-1.0_rp)))
       ENDIF
   ENDIF
ELSE IF ( typ == 'piston' ) THEN
   IF (k.LT.50.) THEN
       resu = CMPLX(0.0_rp, - (2.0_rp * k + &
                    SINH(2.0_rp * k)) / (4.0_rp * SINH(k) * SINH(k)), KIND=CP)
   ELSE
       resu = CMPLX(0.0_rp, - 0.5_rp, KIND=CP)
   ENDIF
ELSE IF ( typ == 'dipole' ) THEN
   d = hinge
   IF ( k <= 20.0_rp ) THEN
      k2   = k * k
      w2   = k * TANH(k)
      w4   = w2 * w2
      resu = (k2 - w4 + w2) / (2.0_rp * SQRT(w2) * (k2 - w4) * COSH(k) * SINH(k*(d+1)))
   ELSE
      k2   = k * k
      w2   = k * TANH(k)
      resu = SQRT(w2) / (2.0_rp * EXP(k*d) * k2)
   END IF
ELSE
   resu = CMPLX(1.0_rp, 0.0_rp, KIND=CP)
END IF
!
TF_adim_new = resu
!
END FUNCTION TF_adim_new
!
!
!
FUNCTION TE_adim_r(nu, N, wmk, d, prec)
!
! adimensional Transfer Function, inspired from ECN report potentiel.tex,
! for evanescent modes
!   a_n = TE_n * a
!
IMPLICIT NONE
!
REAL(RP)            :: nu, d, prec
INTEGER             :: N
CHARACTER(LEN=*)    :: wmk
COMPLEX(CP)         :: TE_adim_r, resu
REAL(RP)            :: alpha
!
alpha = alpha_adim(nu, N, prec)
!
IF ( wmk == 'ECN' ) THEN
   resu = TF_adim_r(nu, wmk, d, prec)
   resu = resu * 4.0_rp * SIN(alpha) / (alpha * (1.0_rp - d) * (2.0_rp * alpha + &
   SIN(2.0_rp * alpha)))
   resu = resu * (alpha * (1.0_rp - d) * SIN(alpha) + COS(alpha) - COS(alpha * d))
ELSE
   WRITE(*,*) 'Unknown wavemaker type for transfer function of evanescent modes' 
   resu = CMPLX(0.0_rp, 0.0_rp, KIND = CP)
END IF
!
TE_adim_r = resu
!
END FUNCTION TE_adim_r
!
!
!
FUNCTION TE_adim_v(nu, N, wmk, d, prec)
!
! adimensional Transfer Function, inspired from ECN report potentiel.tex,
! for evanescent modes
!   a_n = TE_n * a
!
IMPLICIT NONE
!
REAL(RP)                        :: nu, d, prec
INTEGER, DIMENSION(:)           :: N
CHARACTER(LEN=*)                :: wmk
COMPLEX(CP), DIMENSION(SIZE(N)) :: TE_adim_v
INTEGER                         :: ibou
!
!
DO ibou = 1, SIZE(N)
   !
   TE_adim_v(ibou) = TE_adim_r(nu, N(ibou), wmk, d, prec)
   !
END DO
!
END FUNCTION TE_adim_v
!
!
!
FUNCTION ampli_free_adim(x, y, w, z, a_1, a_2)
!
! Evaluation of the contribution to amplitude of second order free mode (prop. or evan.)
! of two first order modes (prop. or evan. as well)
! Computation is done in complex notation, see report potentiel.tex
!
IMPLICIT NONE
!
REAL(RP)    :: w        ! angular frequency of wave at first order
COMPLEX(CP) :: x, y     ! wave numbers at first order (propagative mode : i k,
!                                                       evanescent mode : alpha_n)
COMPLEX(CP) :: z        ! free wave numbers at second order (propagative mode : i k_l,
!                                                             evanescent mode : beta_m)
COMPLEX(CP) :: a_1, a_2 ! amplitudes of first order modes
!
COMPLEX(CP) :: ampli_free_adim, one
COMPLEX(CP), DIMENSION(3) :: tmp
!
one = CMPLX(1.0_rp, 0.0_rp, KIND = CP)
!
tmp(1) = EXP(i * z)
tmp(2) = EXP(2.0_rp * i * z)
tmp(3) = x + y

ampli_free_adim = tmp(3) * a_1 * a_2 * (5.0_rp * w**4 + x**2 + x * y + y**2)&
* (tmp(1) + one / tmp(1))**2 / (tmp(3)**2 - z**2) &
/ (2.0_rp * z + (tmp(2) - one / tmp(2)) / (2.0_rp * i))
!
END FUNCTION ampli_free_adim
!
!
!
FUNCTION dichoto(func,f0,a,b,prec)
!
!     Dichotomy to solve func(x)=f0 between a and b, with a relative accuracy prec on x
!
! inputs :
!   func => function used (declared as EXTERNAL in the calling program)
!   f0   => target value of the function
!   a,b  => limits of the interval of search
!   prec => accuracy wanted on x
!
IMPLICIT NONE
REAL(RP) :: dichoto,f0,a,b,prec,func
REAL(RP) :: x1,x2,x3
!
IF ( ((func(a)-f0)*(func(b)-f0)) >= 0.0_rp ) THEN
   IF ( ABS(func(a)-f0) < tiny ) THEN
      dichoto = a
   RETURN
   ELSE IF (ABS(func(b)-f0) < tiny ) THEN
      dichoto = b
   RETURN
   END IF
   WRITE(*,*) 'a and b are ill positioned in dichoto'
END IF
x1 = a
x2 = b
DO WHILE (ABS(x1-x2) > prec*ABS(x2))
   x3=(x1+x2)/2.0_rp
   IF (ABS(func(x3)-f0) < tiny) THEN
      dichoto = x3
      RETURN
   ELSE IF ((func(a)-f0)*(func(x3)-f0) < 0.0_rp) THEN
      x2 = x3
   ELSE
      x1 = x3
   END IF
END DO
dichoto = x3
!
END FUNCTION dichoto
!
!
!
END MODULE linear_wave
