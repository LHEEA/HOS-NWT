MODULE input
!
! This module contains the input related routines
!  Subroutines :  read_input
!                 read_datum
!                 read_blank_line
!                 build_format
!                 error_message
!                 write_input
!                 write_datum
!                 write_blank_line
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
!
IMPLICIT NONE
!
INTERFACE read_datum
     MODULE PROCEDURE read_datum_i, read_datum_r, read_datum_c
END INTERFACE
!
INTERFACE write_datum
     MODULE PROCEDURE write_datum_i, write_datum_r, write_datum_c
END INTERFACE
!
INTEGER, PARAMETER            :: N_descr = 33
INTEGER, PARAMETER            :: N_tot   = 52
INTEGER, PARAMETER            :: len_form_read = 7
INTEGER, PARAMETER            :: len_form_write = 25
CHARACTER(LEN=len_form_read)  :: format_read(0:4)
CHARACTER(LEN=len_form_write) :: format_write(0:2)
INTEGER                       :: line_counter
CHARACTER(LEN=N_tot)          :: description
!
!
!
CONTAINS
!
!
!
SUBROUTINE read_input(filename)
!
IMPLICIT NONE
! Input variables
CHARACTER(LEN=*), INTENT(IN)  :: filename
! Local variables
INTEGER  :: unit
!
unit = 100
OPEN(unit, FILE=filename)
line_counter = 0
!
CALL read_blank_line(unit)
WRITE(*,*)
CALL read_datum(unit, xlen)             ! Length(m) of the wavetank
CALL read_datum(unit, ylen)             ! Beam(m) of the wavetank
CALL read_datum(unit, h)                ! Depth(m) of the wavetank
CALL read_blank_line(unit)
WRITE(*,*)
CALL read_datum(unit, icase)            ! Choice of computed case
CALL read_blank_line(unit)
WRITE(*,*)
CALL read_datum(unit, islosh)           ! Number of the natural mode
CALL read_datum(unit, aslosh)           ! Amplitude (m) of the mode
CALL read_blank_line(unit)
WRITE(*,*)
CALL read_datum(unit, amp_mono)         ! Amplitude (m) of wave train
CALL read_datum(unit, nu_mono)          ! Frequency (Hz) of wave train
CALL read_datum(unit, theta_mono)       ! Angle (deg) from x-axis
CALL read_datum(unit, ph_mono)          ! Phasis (rad) of wave train
CALL read_datum(unit, ibat)             ! Directional wmk type
CALL read_datum(unit, xd_mono)          ! Wave target distance (ibat=3)
CALL read_blank_line(unit)
WRITE(*,*)
CALL read_datum(unit, file_name)        ! Name of frequency file
CALL read_datum(unit, i_cut)            ! Frequency cut-off
CALL read_datum(unit, nuc_low)          ! Low cut_off frequency (Hz)
CALL read_datum(unit, nuc_high)         ! High cut_off frequency (Hz)
CALL read_blank_line(unit)
WRITE(*,*)
CALL read_datum(unit, Hs)               ! Significant wave height (m)
CALL read_datum(unit, Tp)               ! Peak period (s)
CALL read_datum(unit, gamma)            ! Shape factor (Jonswap)
CALL read_datum(unit, iseed)            ! Seed number for random numb.
CALL read_blank_line(unit)
WRITE(*,*)
CALL read_datum(unit, i_wmk)            ! Nonlinear order of wavemaker
CALL read_datum(unit, igeom)            ! Type (1: piston, 2: hinged)
CALL read_datum(unit, d_hinge)          ! Rotation axis distance
CALL read_datum(unit, iramp)            ! Time ramp
CALL read_datum(unit, Tramp)            ! Time ramp duration
CALL read_datum(unit, T_stop)           ! Time before stopping wmk 
CALL read_blank_line(unit)
WRITE(*,*)
CALL read_datum(unit, iabsnb)           ! Absorption numerical beach
CALL read_datum(unit, xabsf)            ! Beginning front num. beach
CALL read_datum(unit, coeffabsf)        ! Absorption strength front 
CALL read_blank_line(unit)
WRITE(*,*)
CALL read_datum(unit, iprobes)           ! Use of probes
CALL read_datum(unit, itypprobes)        ! 1: on mesh, 2: on specif pts
CALL read_datum(unit, pro_file)          ! Filename of probe positions
CALL read_blank_line(unit)
WRITE(*,*)
CALL read_datum(unit, toler)             ! Time tolerance: RK 4(5)
CALL read_datum(unit, f_out)             ! Output frequency
CALL read_blank_line(unit)
WRITE(*,*)
CALL read_datum(unit, idim)              ! Output: 1-dim. ; 0-nondim.
CALL read_datum(unit, i3d)               ! Free surface plot
CALL read_datum(unit, imodes)            ! Modes plot
CALL read_datum(unit, iwmk)              ! Wavemaker motion plot
CALL read_datum(unit, i_sw)              ! Swense output 1='yes',0='no'
!
! Spinning dipoles removed for now... 
! Useful for someone?
z_dipmono = -0.25_rp
x_dip     = 2.5_rp
n_side    = 6
! No back absorbing zone with wavemaker only!
xabsb     = 0.0_rp 
coeffabsb = 0.0_rp
!
! Specific developments for B600 removed for now
! Keep this here in case limited wavemaker length + possible absorbing zone in y direction
! Tested in december, 2014: working
iwidth    = 1
ywmk      = 0.0_rp
yl        = 0.0_rp
yr        = 0.0_rp
y_rmp     = 0.0_rp
yabsl     = 0.0_rp
yabsr     = 0.0_rp
coeffabsl = 0.0_rp
coeffabsr = 0.0_rp
!
! No 2nd order free-wave correction
! FIXME: not fully validated
! FIXME: Check if OK for Dalrymple regular wave in HOS-NWT (validated in SWEET, F.B. thesis)
! Note that available only if wavemaker 2nd or 3rd order --> change this also?
ifree = 0
!
! Check the different values are correct
! Mode number is mode+1 (1 is constant mode after FFTs)
islosh = islosh +1
!
if (n2.eq.1) theta_mono=0.d0
theta_mono = theta_mono * PI / 180.d0
!
IF (iramp == 0) THEN
   Tramp = 0.0d0
END IF
!
IF (iwidth == 2) THEN
   IF(dabs(ylen-ywmk-yl-yr).ge.1.d-5) THEN
      print*, 'partial width, iwidth=2'
      print*, 'ylen should be equal to ywmk+yl+yr'
      stop 'error: y-lengths not properly defined'
   ENDIF
ELSE
   ywmk  = ylen
   yl    = 0.0d0
   yr    = 0.0d0
   y_rmp = 0.0d0
END IF
!
IF (igeom == 3 .AND. ABS(coeffabsb) > tiny .AND. (x_dip < xlen * xabsb) .AND. (iabsnb==1)) THEN
   STOP 'Spinning dipole(s) in the back numerical beach'
END IF
!
IF (igeom == 2 .AND. ABS(coeffabsb) > tiny .AND. (iabsnb==1)) THEN
   STOP 'Back numerical beach while generating with a flux'
END IF
!
if (iprobes.eq.2) then !Velocity/pressure computation
   i_sw = 1
endif
!
CLOSE(unit)
!
END SUBROUTINE read_input
!
!
!
SUBROUTINE read_datum_i(unit, input)
!
IMPLICIT NONE
! Input variables
INTEGER, INTENT(IN)  :: unit
! Output variables
INTEGER, INTENT(OUT) :: input
!
line_counter = line_counter + 1
CALL build_read_format('I')
READ(unit,format_read,ERR=101,END=101) description, input
WRITE(*,'(A," : ",I5)') description(1:N_descr), input
RETURN
101 CALL error_message(description)
!
END SUBROUTINE read_datum_i
!
!
!
SUBROUTINE read_datum_r(unit, input)
!
IMPLICIT NONE
! Input variables
INTEGER, INTENT(IN)   :: unit
! Output variables
REAL(RP), INTENT(OUT) :: input
!
line_counter = line_counter + 1
CALL build_read_format('R')
READ(unit,format_read,ERR=101,END=101) description, input
WRITE(*,'(A," : ",ES25.16)') description(1:N_descr), input
IF (ABS(input) > tiny .AND. ABS(input) * 1.0E+16_rp < 1.0E+5_rp) THEN
   WRITE(*,'(A,A)') 'Numeric point is probably missing in current input ',description
   STOP
END IF
RETURN
101 CALL error_message(description)
!
END SUBROUTINE read_datum_r
!
!
!
SUBROUTINE read_datum_c(unit, input)
!
IMPLICIT NONE
! Input variables
INTEGER, INTENT(IN)   :: unit
! Output variables
CHARACTER(LEN=*), INTENT(OUT) :: input
!
line_counter = line_counter + 1
CALL build_read_format('A')
READ(unit,format_read,ERR=101,END=101) description, input
WRITE(*,'(A," : ",A)') description(1:N_descr), input
RETURN
101 CALL error_message(description)
!
END SUBROUTINE read_datum_c
!
!
!
SUBROUTINE read_blank_line(unit)
!
IMPLICIT NONE
! Input variables
INTEGER, INTENT(IN)   :: unit
!
line_counter = line_counter + 1
READ(unit,*) 
!
END SUBROUTINE read_blank_line
!
!
!
SUBROUTINE build_read_format(code)
!
IMPLICIT NONE
!
CHARACTER(LEN=*) :: code
!
format_read(0) = '(A'
format_read(1) = TRIM(int2str(N_tot))
format_read(2) = ','
SELECT CASE (code)
   CASE('I')          ! Integer
      format_read(3) = 'I5'
   CASE('F','R')      ! Real number
      format_read(3) = 'ES25.16'
   CASE('S','C','A')  ! Character string
      format_read(3) = 'A'
END SELECT
format_read(4) = ')'
! WRITE(*,*) format_read
!
END SUBROUTINE build_read_format
!
!
!
SUBROUTINE error_message(description)
!
IMPLICIT NONE
! Input variables
CHARACTER(LEN=N_tot) :: description
!
WRITE(*,'(A,I2)') 'Error while reading the input file on line: ', line_counter 
WRITE(*,'(A)') description
STOP

!
END SUBROUTINE error_message
!
!
!
SUBROUTINE write_input(unit)
!
IMPLICIT NONE
! Input variables
INTEGER, INTENT(IN) :: unit
!
CALL write_blank_line(unit,'----------------------- Tank dimensions -----------------------')
CALL write_datum(unit, xlen,              'xlen',           'Length(m) of the wavetank')
CALL write_datum(unit, ylen,              'ylen',           'Beam(m) of the wavetank')
CALL write_datum(unit, h,                 'h',              'Depth(m) of the wavetank')
!
CALL write_blank_line(unit,'------------------------ Select case --------------------------')
CALL write_datum(unit, icase,             'icase',          'Choice of computed case')
!
CALL write_blank_line(unit,'------------------ Sloshing case (icase = 1) ------------------')
CALL write_datum(unit, islosh,            'islosh',         'Number of the natural mode')
CALL write_datum(unit, aslosh,            'aslosh',         'Amplitude (m) of the mode')
!
CALL write_blank_line(unit,'--------------- Monochromatic case (icase = 2) ----------------')
CALL write_datum(unit, amp_mono,          'amp_mono',       'Amplitude (m) of wave train')
CALL write_datum(unit, nu_mono,           'nu_mono',        'Frequency (Hz) of wave train')
CALL write_datum(unit, theta_mono,        'theta_mono',     'Angle (deg) from x-axis')
CALL write_datum(unit, ph_mono,           'ph_mono',        'Phasis (rad) of wave train')
CALL write_datum(unit, ibat,              'ibat',           'Directional wmk type')
CALL write_datum(unit, xd_mono,           'xd_mono',        'Wave target distance (ibat=3)')
!
CALL write_blank_line(unit,'-------------------- File case (icase = 3) --------------------')
CALL write_datum(unit, file_name,         'file_name',      'Name of frequency file')
CALL write_datum(unit, i_cut,             'i_cut',          'Frequency cut-off')
CALL write_datum(unit, nuc_low,           'nuc_low',        'Low cut_off frequency (Hz)')
CALL write_datum(unit, nuc_high,          'nuc_high',       'High cut_off frequency (Hz) ')
!
CALL write_blank_line(unit,'------------------ Irregular wave (icase=4) -------------------')
CALL write_datum(unit, Hs,                'Hs',             'Significant wave height (m)')
CALL write_datum(unit, Tp,                'Tp',             'Peak period (s)')
CALL write_datum(unit, gamma,             'gamma',          'Shape factor (Jonswap)')
CALL write_datum(unit, iseed,             'iseed',          'Seed number for random numb. ')
!
CALL write_blank_line(unit,'-------------------- Wavemaker definition ---------------------')
CALL write_datum(unit, i_wmk,             'i_wmk',          'Nonlinear order of wavemaker')
CALL write_datum(unit, igeom,             'igeom',          'Type (1: piston, 2: hinged)')
CALL write_datum(unit, d_hinge,           'd_hinge',        'Rotation axis distance')
CALL write_datum(unit, iramp,             'iramp',          'Time ramp')
CALL write_datum(unit, Tramp,             'Tramp',          'Time ramp duration')
CALL write_datum(unit, T_stop,            'T_stop',         'Time before stopping wmk')
!
CALL write_blank_line(unit,'----------------------- Numerical beach -----------------------')
CALL write_datum(unit, iabsnb,            'iabsnb',         'Absorption numerical beach')
CALL write_datum(unit, xabsf,             'xabsf',          'Beginning front num. beach')
CALL write_datum(unit, coeffabsf,         'coeffabsf',      'Absorption strength front')
!
CALL write_blank_line(unit,'------------- Elevation/Velocity-Pressure probes --------------')
CALL write_datum(unit, iprobes,           'iprobes',        'Use of probes')
CALL write_datum(unit, itypprobes,        'itypprobes',     '1: on mesh, 2: on specif pts')
CALL write_datum(unit, pro_file,          'pro_file',       'Filename of probe positions')
!
CALL write_blank_line(unit,'---------------------- Time-integration -----------------------')
CALL write_datum(unit, toler,             'toler',          'Time tolerance: RK 4(5)')
CALL write_datum(unit, f_out,             'f_out',          'Output frequency')
!
CALL write_blank_line(unit,'--------------------------- Output ----------------------------')
CALL write_datum(unit, idim,              'idim',           'Output: 1-dim. ; 0-nondim.')
CALL write_datum(unit, i3d,               'i3d',            'Free surface plot')
CALL write_datum(unit, imodes,            'imodes',         'Modes plot')
CALL write_datum(unit, iwmk,              'iwmk',           'Wavemaker motion plot')
CALL write_datum(unit, i_sw,              'i_sw',           'Swense output 1=''yes'',0=''no''')
!
CALL write_blank_line(unit,'--------------------- Numerical parameters ---------------------')
CALL write_datum(unit, n1,                'n1',             'Number of nodes/modes in x')
CALL write_datum(unit, n2,                'n2',             'Number of nodes/modes in y')
CALL write_datum(unit, n3,                'n3',             'Number of nodes/modes on wmk')
CALL write_datum(unit, mHOS,              'mHOS',           'HOS order')
CALL write_datum(unit, p1,                'p1',             'Dealiasing in x-direction')
CALL write_datum(unit, p2,                'p2',             'Dealiasing in y-direction')
CALL write_datum(unit, coeffilt(1),       'coeffilt(1)',    'Filtering ratio in x-modes')
CALL write_datum(unit, coeffilt(2),       'coeffilt(2)',    'Filtering ratio in y-modes')
CALL write_datum(unit, coeffilt(3),       'coeffilt(3)',    'Filtering ratio in z-modes')
!
END SUBROUTINE write_input
!
!
!
SUBROUTINE write_datum_i(unit, input, variable, text)
!
IMPLICIT NONE
! Input variables
INTEGER, INTENT(IN)  :: unit
CHARACTER(LEN=*), INTENT(IN) :: text, variable
INTEGER, INTENT(IN) :: input
!
CALL build_write_format('I')
WRITE(unit,format_write) text, variable, input
!
END SUBROUTINE write_datum_i
!
!
!
SUBROUTINE write_datum_r(unit, input, variable, text)
!
IMPLICIT NONE
! Input variables
INTEGER, INTENT(IN)   :: unit
CHARACTER(LEN=*), INTENT(IN) :: text, variable
REAL(RP), INTENT(IN) :: input
!
CALL build_write_format('R')
WRITE(unit,format_write) text, variable, input
!
END SUBROUTINE write_datum_r
!
!
!
SUBROUTINE write_datum_c(unit, input, variable, text)
!
IMPLICIT NONE
! Input variables
INTEGER, INTENT(IN)   :: unit
CHARACTER(LEN=*), INTENT(IN) :: text, variable
CHARACTER(LEN=*), INTENT(IN) :: input
!
CALL build_write_format('A')
WRITE(unit,format_write) text, variable, input
!
END SUBROUTINE write_datum_c
!
!
!
SUBROUTINE write_blank_line(unit, text)
!
IMPLICIT NONE
! Input variables
INTEGER, INTENT(IN)          :: unit
CHARACTER(LEN=*), INTENT(IN) :: text
!
WRITE(unit,'("#",A)') text
!
END SUBROUTINE write_blank_line
!
!
!
SUBROUTINE build_write_format(code)
!
IMPLICIT NONE
!
CHARACTER(LEN=*) :: code
!
format_write(0) = '("#",A'//TRIM(int2str(N_descr))//',":: ",A'//TRIM(int2str(N_tot-(N_descr+3)-3))//',":: ",'
SELECT CASE (code)
   CASE('I')          ! Integer
      format_write(1) = 'I5'
   CASE('F','R')      ! Real number
      format_write(1) = 'ES25.16'
   CASE('S','C','A')  ! Character string
      format_write(1) = 'A'
END SELECT
format_write(2) = ')'
!
END SUBROUTINE build_write_format
!
!
!
END MODULE input
