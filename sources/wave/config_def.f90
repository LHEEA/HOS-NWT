MODULE config_def
!
! This module allows the definition of the wave-tank configuration (especially wavemaker)
! ECN ocean wave basin taken as default values
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
TYPE config
   INTEGER          :: rnum
   REAL(RP)         :: clock
   REAL(RP)         :: depth
   REAL(RP)         :: ramp
   CHARACTER(LEN=3) :: typ_ramp
   REAL(RP)         :: cutoff_low
   REAL(RP)         :: cutoff_high
   REAL(RP)         :: Ly
   CHARACTER(LEN=6) :: typ_wmk
   REAL(RP)         :: hinge
END TYPE config
!
!
!
CONTAINS
!
!
!
FUNCTION build_config(rnum, clock, depth, ramp, typ_ramp, cutoff_low, cutoff_high, Ly, typ_wmk, hinge)
!
INTEGER                    :: rnum
REAL(RP), OPTIONAL         :: clock
REAL(RP), OPTIONAL         :: depth
REAL(RP), OPTIONAL         :: ramp
CHARACTER(LEN=3), OPTIONAL :: typ_ramp
REAL(RP), OPTIONAL         :: cutoff_high, cutoff_low
REAL(RP), OPTIONAL         :: Ly
CHARACTER(LEN=6), OPTIONAL :: typ_wmk
REAL(RP), OPTIONAL         :: hinge
!
TYPE(config)     :: build_config
!
build_config%rnum = rnum
!
IF (PRESENT(clock)) THEN
   build_config%clock = clock
ELSE
   build_config%clock = 32.0_rp
END IF
!
IF (PRESENT(depth)) THEN
   build_config%depth = depth
ELSE
   build_config%depth = 5.0_rp
END IF
!
IF (PRESENT(ramp)) THEN
   build_config%ramp = ramp
ELSE
   build_config%ramp = 32.0_rp
END IF
!
IF (PRESENT(typ_ramp)) THEN
   build_config%typ_ramp = typ_ramp
ELSE
   build_config%typ_ramp = 'lin'
END IF
!
IF (PRESENT(cutoff_high)) THEN
   build_config%cutoff_high = cutoff_high
ELSE
   build_config%cutoff_high = 10.0_rp
END IF
!
IF (PRESENT(cutoff_low)) THEN
   build_config%cutoff_low = cutoff_low
ELSE
   build_config%cutoff_low = 0.0_rp
END IF
!
IF (PRESENT(Ly)) THEN
   build_config%Ly = Ly
ELSE
   build_config%Ly = 29.74_rp
END IF
!
IF (PRESENT(typ_wmk)) THEN
   build_config%typ_wmk = typ_wmk
ELSE
   build_config%typ_wmk = 'hinged'
END IF
!
IF (PRESENT(hinge)) THEN
   build_config%hinge = hinge
ELSE
   build_config%hinge = 2.147_rp
END IF
!
END FUNCTION build_config
!
!
!
SUBROUTINE display_config(cfg)
!
IMPLICIT NONE
!
TYPE(config) :: cfg
!
WRITE(*,*) 'rnum        ', cfg%rnum
WRITE(*,*) 'clock       ', cfg%clock
WRITE(*,*) 'depth       ', cfg%depth
WRITE(*,*) 'ramp        ', cfg%ramp
WRITE(*,*) 'typ_ramp    ', cfg%typ_ramp
WRITE(*,*) 'cutoff low  ', cfg%cutoff_low
WRITE(*,*) 'cutoff high ', cfg%cutoff_high
WRITE(*,*) 'Ly          ', cfg%Ly
WRITE(*,*) 'typ_wmk     ', cfg%typ_wmk
WRITE(*,*) 'hinge       ', cfg%hinge
!
END SUBROUTINE display_config
!
!
!
END MODULE config_def
