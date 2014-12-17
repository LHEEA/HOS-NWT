MODULE output_post_process
!
! This module contains the input related routines
!  Subroutines :  read_input
!                 read_datum
!                 read_blank_line
!                 build_format
!                 error_message
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
!
!
CONTAINS
!
!
!
SUBROUTINE init_output_post_process(i_card)
!
IMPLICIT NONE
!
! Input variables
INTEGER, INTENT(IN) :: i_card
! Local variables
!
!INTEGER :: i1, i2
!
!
IF (i_card /= 0) THEN
	!
	! Velocity and pressure card
	IF (i_card == 1) THEN
    	! Tecplot output
    	OPEN(31,FILE='Results/VP_card.dat')
		WRITE(31,'(A)') 'TITLE =" Velocity and pressure field "'
    	WRITE(31,'(A)') 'VARIABLES="x","y","z","vitx","vity","vitz","Press"'
    ELSEIF (i_card == 2) THEN ! Boundary fitted coordinates
    	! Tecplot output
    	OPEN(32,FILE='Results/VP_card_fitted.dat')
    	WRITE(32,'(A)') 'TITLE =" Velocity and pressure field "'
    	WRITE(32,'(A)') 'VARIABLES="x","y","z","vitx","vity","vitz","Press"'
    ENDIF
    !
ENDIF
!
END SUBROUTINE init_output_post_process
!
!
!
SUBROUTINE output_time_step_card(i_card,tecplot,time,dt_out_star,zlocal,z_min,z_max,T_start,L_out,T_out,i_test,&
		imin,imax,jmin,jmax,i_zvect,vitx,vity,vitz,phit)
!
! This subroutine performs the output of velocity/pressure card
!
IMPLICIT NONE
!
! Input variables
INTEGER, INTENT(IN)                  :: i_card, tecplot
REAL(RP), INTENT(IN)                 :: time, dt_out_star
!
REAL(RP), INTENT(IN)                 :: T_start,L_out,T_out,z_min,z_max
INTEGER, INTENT(IN)                  :: i_test,imin,imax,jmin,jmax,i_zvect
REAL(RP), DIMENSION(:,:), INTENT(IN) :: zlocal,vitx,vity,vitz,phit
! Local variables
!
INTEGER  :: i1,i2
REAL(RP) :: Press
!
REAL(RP)             :: tiny_sp
!
! tiny_sp is single precision: useful for inequalities check with values read from files
tiny_sp = epsilon(1.0)
!
IF (i_card /= 0) THEN
	!
	! Velocity and pressure card
	IF (i_card == 1) THEN
		!
   		! These are informations useful for eventual coupling using files VP_card
    	IF (time*T_out <= T_start+tiny_sp) THEN ! First time-step
    		OPEN(30,file='Results/data_VP_card.dat',status='unknown')
			WRITE(30,'(2(ES16.9,X))') x(imin)*L_out, x(imax)*L_out
			WRITE(30,'(2(ES16.9,X))') y(jmin)*L_out, y(jmax)*L_out
			WRITE(30,'(2(ES16.9,X))') z_min, z_max ! should be used only for (i_card == 1)
			WRITE(30,'(3(I4,X))') imax-imin+1, jmax-jmin+1, i_zvect
			WRITE(30,'(3(ES16.9,X),I6)') T_start, T_stop, dt_out_star*T_out, FLOOR((T_stop-T_start)/T_out/dt_out_star) !T_stop, T_start are dimensional quantities
			CLOSE(30)
		ENDIF
		!
		IF (i_test == 1) THEN ! First element in the z-loop
			IF (time*T_out <= T_start+tiny_sp) THEN ! First time-step
				IF (tecplot == 11) THEN
					WRITE(31,104)'ZONE SOLUTIONTIME = ',time*T_out,', I=', imax-imin+1,' J=', jmax-jmin+1,' K=', i_zvect
				ELSE
					WRITE(31,104)'ZONE T = "',time*T_out,'", I=', imax-imin+1,' J=', jmax-jmin+1,' K=', i_zvect
				END IF
			ELSE ! Following time-steps
				IF (tecplot == 11) THEN
					WRITE(31,104)'ZONE SOLUTIONTIME = ',time*T_out,', D=(1,2,3), I=', imax-imin+1,' J=', jmax-jmin+1,' K=', i_zvect
				ELSE
					WRITE(31,104)'ZONE T = "',time*T_out,'", D=(1,2,3), I=', imax-imin+1,' J=', jmax-jmin+1,' K=', i_zvect
				END IF
			ENDIF
		ENDIF
		!
		DO i2=1,jmax-jmin+1
			DO i1=1,imax-imin+1
				Press = - zlocal(i1,i2) - 0.5_rp*(vitx(i1,i2)**2+vity(i1,i2)**2+vitz(i1,i2)**2)-phit(i1,i2)
				if(time*T_out <= T_start+tiny_sp) then
					WRITE(31,'(7(ES16.9,X))') x(i1+imin-1)*L_out, y(i2+jmin-1)*L_out, zlocal(i1,i2)*L_out, &
							vitx(i1,i2)*L_out/T_out, vity(i1,i2)*L_out/T_out, vitz(i1,i2)*L_out/T_out, &
							Press*L_out**2/T_out**2
				else
					WRITE(31,'(4(ES16.9,X))') vitx(i1,i2)*L_out/T_out, vity(i1,i2)*L_out/T_out, &
							vitz(i1,i2)*L_out/T_out, Press*L_out**2/T_out**2
				endif
			ENDDO
		ENDDO
    ELSEIF (i_card == 2) THEN
    	IF (i_test == 1) THEN ! First element in the z-loop
			IF (tecplot == 11) THEN
				WRITE(32,104)'ZONE SOLUTIONTIME = ',time*T_out,', I=', imax-imin+1,' J=', jmax-jmin+1,' K=', i_zvect
			ELSE
				WRITE(32,104)'ZONE T = "',time*T_out,', I=', imax-imin+1,' J=', jmax-jmin+1,' K=', i_zvect
			END IF    	
    	ENDIF
		!
		DO i2=1,jmax-jmin+1
        	DO i1=1,imax-imin+1
				Press = - zlocal(i1,i2) - 0.5_rp*(vitx(i1,i2)**2+vity(i1,i2)**2+vitz(i1,i2)**2)-phit(i1,i2)
				WRITE(32,'(7(ES16.9,X))') x(i1+imin-1)*L_out, y(i2+jmin-1)*L_out, zlocal(i1,i2)*L_out, vitx(i1,i2)*L_out/T_out, &
					vity(i1,i2)*L_out/T_out, vitz(i1,i2)*L_out/T_out, Press*L_out**2/T_out**2 
			ENDDO
		ENDDO
    ENDIF
ENDIF
!
104 FORMAT(A,F9.2,A,I5,A,I5,A,I5)
                  
END SUBROUTINE output_time_step_card
!
END MODULE output_post_process