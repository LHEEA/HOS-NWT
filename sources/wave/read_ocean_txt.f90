MODULE read_ocean_txt
!
! This module allows the use of .txt files from wave generation software used at ECN wave basins
! This software, named ocean, may generate a .txt file containing wavemaker movement informations
! The file is read and transformed  into classical .dat and .cfg files used in HOS-NWT
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
USE linear_wave
!
CONTAINS
!
SUBROUTINE ocean_txt2cfg_dat(filename, i_cut, nuc_low, nuc_high, i_case, clock, h, ylen, d_hinge, igeom)
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN) :: filename
INTEGER, INTENT(IN)          :: i_cut, i_case
REAL(RP), INTENT(IN)         :: nuc_low, nuc_high
! Optional
INTEGER, INTENT(IN), OPTIONAL  :: igeom
REAL(RP), INTENT(IN), OPTIONAL :: clock, h, ylen, d_hinge
!
INTEGER                      :: nlines, iloop
INTEGER                      :: i_scan
CHARACTER(LEN=200)           :: temp
CHARACTER(LEN=200), ALLOCATABLE, DIMENSION(:) :: temp_array
INTEGER                      :: rnum, ramp, rcomp, n_front, n_start
LOGICAL                      :: comp
REAL(RP)                     :: clock_tmp, f_base, f_c_high, f_c_low, k
INTEGER,  ALLOCATABLE, DIMENSION(:) :: harmo
REAL(RP), ALLOCATABLE, DIMENSION(:) :: ampli, freq, angle, phase
!
WRITE(*,'(A)') 'Reading a .txt file from ocean...'
!
OPEN(999,FILE=TRIM(filename)//'.txt')
!
! Number of lines...
nlines = 0
DO
   READ(999,*,END=98)
   nlines = nlines + 1
   CYCLE
98 EXIT
END DO
REWIND(999)
!
WRITE(*,'(A,i5,A)') 'File is containing ',nlines,' lines.'
!
ALLOCATE(temp_array(nlines))
!
comp = .FALSE.
DO iloop = 1, nlines
   READ(999,'(200A)') temp
   temp_array(iloop) = temp
   ! run rnumber
   i_scan = INDEX(temp, 'rnumber  =')
   IF (i_scan /= 0) THEN
      rnum = convert(temp, 18)
      WRITE(*,'(A,I2)') 'rnum = ',rnum
   END IF
   ! run ramp length
   i_scan = INDEX(temp, 'ramp     =')
   IF (i_scan /= 0) THEN
      ramp = convert(temp, 18)
      WRITE(*,'(A,I4)') 'ramp = ',ramp
   END IF
   ! run rnum compression factor
   i_scan = INDEX(temp, 'rnumber compression')
   IF (i_scan /= 0) THEN
      comp  = .TRUE.
      rcomp = convert(temp, 62)
      WRITE(*,'(A,I2)') 'rcomp = ',rcomp
   END IF
   ! run number of fronts
   i_scan = INDEX(temp, 'Wave has ')
   IF (i_scan /= 0) THEN
      n_front = convert(temp, 10, -1)
      n_start = iloop + 2
      WRITE(*,'(A,I4)') 'n_front = ',n_front
   END IF
END DO
!
REWIND(999)
!
DO iloop = 1, n_start-1
   READ(999,'(200A)')
END DO
!
ALLOCATE(harmo(n_front), ampli(n_front), freq(n_front), angle(n_front), phase(n_front))
!
DO iloop = 1, n_front
   READ(999,*) harmo(iloop), freq(iloop), ampli(iloop), angle(iloop), phase(iloop)
END DO
!
CLOSE(999)
!
IF (comp) THEN
   harmo = harmo * rcomp
END IF
!
IF (PRESENT(clock)) THEN
	clock_tmp = clock
ELSE
	clock_tmp = 32.0_rp ! For ECN wave basins
ENDIF
f_base = clock_tmp / 2**rnum
freq   = harmo * f_base
!
! FIXME: where does it come from?
! Rk: act only in 3D (angle/=0)
IF (i_case == 31) THEN ! ECN wave basin
	DO iloop = 1, n_front
   		k = wave_number_r(freq(iloop), 5.0_rp, 1.0e-12_rp)
   		phase(iloop) = phase(iloop) + k * 29.74_rp / 2.0_rp * SIN(angle(iloop))
	END DO
ELSEIF (i_case == 32) THEN ! ECN towing tank
	DO iloop = 1, n_front
   		k = wave_number_r(freq(iloop), 2.88_rp, 1.0e-12_rp)
	END DO
ELSE
	!FIXME: test if this is necessary
	IF(PRESENT(h) .AND. PRESENT(ylen)) THEN
		DO iloop = 1, n_front
			k = wave_number_r(freq(iloop), h, 1.0e-12_rp)
			phase(iloop) = phase(iloop) + k * ylen / 2.0_rp * SIN(angle(iloop))
		END DO	
	ELSE
		print*,'Unknown type for .txt files'
		stop
	ENDIF
ENDIF
!
angle = 180.0_rp * angle / PI
!
! frequency bandpass filter
IF (i_cut == 1) THEN
   f_c_low  = nuc_low
   f_c_high = nuc_high
ELSE
   f_c_low  = 0.0_rp
   f_c_high = 100.0_rp
END IF

OPEN(999,FILE=TRIM(filename)//'.cfg')
!
IF (i_case == 31) THEN ! ECN wave basin
	WRITE(999,'(I2)')   rnum
	WRITE(999,'(F5.2)') clock_tmp
	WRITE(999,'(F5.2)') 5.0_rp
	WRITE(999,'(I4)')   ramp
	WRITE(999,'(A)')    'lin'
	WRITE(999,'(F5.2)') f_c_low
	WRITE(999,'(F6.2)') f_c_high
	WRITE(999,'(F5.2)') 30.00_rp
	WRITE(999,'(A)')    'hinged'
	WRITE(999,'(F6.3)') 2.147_rp
ELSEIF (i_case == 32) THEN ! ECN towing tank
	WRITE(999,'(I2)')   rnum
	WRITE(999,'(F5.2)') clock_tmp
	WRITE(999,'(F5.2)') 2.88_rp
	WRITE(999,'(I4)')   ramp
	WRITE(999,'(A)')    'lin'
	WRITE(999,'(F5.2)') f_c_low
	WRITE(999,'(F6.2)') f_c_high
	WRITE(999,'(F5.2)') 5.00_rp
	WRITE(999,'(A)')    'hinged'
	WRITE(999,'(F6.3)') 0.0_rp
ELSE
	IF (PRESENT(h).AND.PRESENT(d_hinge).AND.PRESENT(igeom)) THEN
		WRITE(999,'(I2)')   rnum
		WRITE(999,'(F5.2)') clock_tmp
		WRITE(999,'(F5.2)') h
		WRITE(999,'(I4)')   ramp
		WRITE(999,'(A)')    'lin'
		WRITE(999,'(F5.2)') f_c_low
		WRITE(999,'(F6.2)') f_c_high
		WRITE(999,'(F5.2)') ylen*h
		IF (igeom.EQ.2) THEN
			WRITE(999,'(A)')    'piston'
		ELSE
			WRITE(999,'(A)')    'hinged'
		ENDIF
		WRITE(999,'(F6.3)') d_hinge*h
	ELSE
		print*,'unknown types in read_ocean_txt'
		stop
	ENDIF
ENDIF
!
CLOSE(999)
!
OPEN(999,FILE=TRIM(filename)//'.dat')
!
DO iloop=1, n_front
   IF (ampli(iloop) < 1.0e-98_rp) THEN
      ampli(iloop) = 0.0_rp
   END IF
END DO
DO iloop=1, n_front
   IF (freq(iloop) >= f_c_low .AND. freq(iloop) <= f_c_high) THEN
      WRITE(999,'(I5,4(X,G20.7))') harmo(iloop), freq(iloop), ampli(iloop), angle(iloop), phase(iloop)
   END IF
END DO
!
CLOSE(999)
!
END SUBROUTINE ocean_txt2cfg_dat
!
!
!
FUNCTION convert(temp, i_unit, i_way, i_char)
!
IMPLICIT NONE
!
CHARACTER(LEN=200) :: temp
INTEGER :: i_unit
INTEGER, OPTIONAL :: i_way, i_char
INTEGER :: i_tmp, i_size, i_start, i_stop
INTEGER :: convert, i_test
!
IF (PRESENT(i_way)) THEN
   IF (present(i_char)) THEN
      i_test = i_char
   ELSE
      i_test = 32 !assume you look for blank in general case
   ENDIF
   i_tmp  = i_unit
   i_size = 0
   DO WHILE (IACHAR(temp(i_tmp:i_tmp)) /= i_test)
      i_size = i_size + i_way
      i_tmp  = i_tmp  - i_way
   END DO
ELSE
   i_tmp  = i_unit
   i_size = 0
   DO WHILE (IACHAR(temp(i_tmp:i_tmp)) /= 32)
      i_size = i_size + 1
      i_tmp  = i_tmp  - 1
   END DO
END IF
!
IF (i_size < 0) THEN
   i_size  = - i_size
   i_start = i_unit
   i_stop  = i_unit + i_size - 1
ELSE
   i_start = i_unit - i_size + 1
   i_stop  = i_unit
END IF
!
convert = strtoint(temp(i_start:i_stop))
!
END FUNCTION convert
!
!
!
END MODULE read_ocean_txt