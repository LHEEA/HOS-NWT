MODULE output
!
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
USE input
USE Fourier_FFTW
USE velocities
!
!
!
CONTAINS
!
!
!
SUBROUTINE init_output(i3d,imodes,iprobes,iwmk,igeom,i_sw,idim)
!
IMPLICIT NONE
!
! Input variables
INTEGER, INTENT(IN) :: i3d,imodes,iprobes,iwmk,igeom,i_sw,idim
! Local variables
INTEGER :: i1, ii
!
! Free-surface elevation
if (i3d.eq.1)  then
   open(1,file='Results/3d.dat',status='unknown')
   call write_input(1)
   !
   write(1,*)'TITLE=" 3D free surface elevation "'
   IF (idim == 0) THEN
      write(1,'(A)') 'VARIABLES="x/h","y/h","eta/h","phis*sqrt(h/g)/h**2"'
   ELSE
      write(1,'(A)') 'VARIABLES="x","y","eta","phis"'
   END IF
end if
!
! Volume and energy
open(2,file='Results/vol_energy.dat',status='unknown')
call write_input(2)
!
write(2,'(A)')'TITLE=" volume and energy evolution "'
!
IF (idim == 0) THEN
   write(2,'(A)') 'VARIABLES="time*sqrt(g/h)", "voltot", "enepot", "enekin", "enetot"'
ELSE
   write(2,'(A)') 'VARIABLES="time", "voltot", "enepot", "enekin", "enetot"'
END IF
!
! Modal description
if (imodes.eq.1) then
   open(5,file='Results/HOS_modes.dat',status='unknown')
   call write_input(5)
   !
   write(5,'(A)') 'TITLE = " 3D spectral potential modes"'
   IF (idim == 0) THEN
      write(5,500) 'VARIABLES = "kx*h", "ky*h","a_eta/h","a_phis*sqrt(h/g)/h**2","LOG10-a-eta","LOG10-a-phis"'
   ELSE
      write(5,500) 'VARIABLES = "kx", "ky","a_eta","a_phis","LOG10-a-eta","LOG10-a-phis"'
   ENDIF
end if
500 format(15A)
!
! Wavemaker motion
if (iwmk.ne.0 .AND. igeom/=3) then
   open(8,file='Results/wmk_motion.dat',status='unknown')
   call write_input(8)
   IF ( n2 == 1) THEN
      write(8,*)'TITLE="Wavemaker motion against time"'
      IF (idim == 0) THEN
         write(8,800)'VARIABLES="time*sqrt(g/h)", "X-top/h", "Vx-top/sqrt(g*h)"'
      ELSE
         write(8,800)'VARIABLES="time", "X-top", "Vx-top"'
      END IF
   ELSE
      write(8,*)'TITLE="Wavemaker motion against time and transverse direction"'
      IF (idim == 0) THEN
         write(8,800)'VARIABLES="y/h", "X-top/h", "Vx-top/sqrt(g*h)"'
      ELSE
         write(8,800)'VARIABLES="y", "X-top", "Vx-top"'
      END IF
   END IF
800 format(a)
endif
!
! Wave probes/Pressure probes location
if (iprobes.eq.1) then
   WRITE(*,'(A)') 'Probes:'
   OPEN(55,FILE=pro_file)
   nprobes=0
   DO
      READ(55,*,END=98) 
      nprobes = nprobes + 1
      CYCLE
98    EXIT
   END DO
   REWIND(55)
   !
   write(*,'(I3,2A)') nprobes, ' probes found in the file ',pro_file
   !
   if (nprobes.gt.maxprobes) stop 'error: increase maxprobes in MAIN parameters'
   DO i1=1,nprobes
      IF (n2 == 1) THEN
         READ(55,*) xprobe(i1)
         yprobe(i1) = 0.0d0
      ELSE
         READ(55,*) xprobe(i1), yprobe(i1)
      END IF
   END DO
   CLOSE(55)
   !
   WRITE(*,'(A)') 'Probes position:'
   !
   DO i1=1,nprobes
      write(*,902)'xprobe(',i1,')=',xprobe(i1),' m, yprobe(',i1,')=',yprobe(i1),' m'
   END DO
   !
   xprobe = xprobe / h
   yprobe = yprobe / h
END IF
!
902 format(a,2(i2,a,1es10.3,a))
!
if (iprobes.eq.2) then
   WRITE(*,'(A)') 'Press:'
   OPEN(55,FILE=pro_file)
   nprobes=0
   DO
      READ(55,*,END=998) 
      nprobes = nprobes + 1
      CYCLE
998   EXIT
   END DO
   REWIND(55)
   !
   write(*,'(I3,2A)') nprobes, ' pressure probes found in the file ',pro_file
   !
   if (nprobes.gt.maxprobes) stop 'error: increase maxprobes in MAIN parameters'
   DO i1=1,nprobes
      IF (n2 == 1) THEN
         READ(55,*) xprobe(i1), zprobe(i1)
         yprobe(i1) = 0.0_rp
      ELSE
         READ(55,*) xprobe(i1), yprobe(i1), zprobe(i1)
      END IF
      IF (zprobe(i1).LT.-h) THEN
         write(*,*) 'error on vertical position of pressure probe'
         stop
      ENDIF
   END DO
   CLOSE(55)
   !
   WRITE(*,'(A)') 'Pressure probes position:'
   !
   DO i1=1,nprobes
      write(*,903)'xprobe(',i1,')=',xprobe(i1),' m, yprobe(',i1,')=',yprobe(i1),' m and zprobe( ',i1,')=',zprobe(i1),' m'
   END DO
   xprobe = xprobe / h
   yprobe = yprobe / h
   zprobe = zprobe / h
END IF
!
903 format(a,3(i2,a,1es10.3,a))
!
!  WAVE PROBES
if (iprobes.eq.1) then
   open(9,file='Results/probes.dat',status='unknown')
   call write_input(9)
   write(9,'(A)') 'TITLE="Probes elevation records versus time"'		!
   IF (idim == 0) THEN
      write(9,'(200A)') 'VARIABLES = "time*sqrt(g/h)" ',('"p'//TRIM(int2str(ii))//'/h" ',ii=1,nprobes)
   ELSE
      write(9,'(200A)') 'VARIABLES = "time" ',          ('"p'//TRIM(int2str(ii))//'" ',ii=1,nprobes)
   END IF
end if
!
!  VELOCITIES AND PRESSURE ON PROBES
if (iprobes.eq.2) then
   open(99,file='Results/probes.dat',status='unknown')
   call write_input(99)
   write(99,'(A)') 'TITLE="Velocities/pressure on probes records versus time"'
   IF (idim == 0) THEN
      write(99,'(2000A)') 'VARIABLES = "time*sqrt(g/h)" ',('"etap'//TRIM(int2str(ii))//'/h" ',ii=1,nprobes), &
           ('"vitxp'//TRIM(int2str(ii))//'*sqrt(h/g)/h" ',ii=1,nprobes), &
           ('"vityp'//TRIM(int2str(ii))//'*sqrt(h/g)/h" ',ii=1,nprobes), &
           ('"vitzp'//TRIM(int2str(ii))//'*sqrt(h/g)/h" ',ii=1,nprobes), &
           ('"Pressp'//TRIM(int2str(ii))//'/(rho*g)" ',ii=1,nprobes)
   ELSE
      write(99,'(2000A)') 'VARIABLES = "time" ',          ('"etap'//TRIM(int2str(ii))//'" ',ii=1,nprobes),&
           ('"vitxp'//TRIM(int2str(ii))//'" ',ii=1,nprobes),&
           ('"vityp'//TRIM(int2str(ii))//'" ',ii=1,nprobes),&
           ('"vitzp'//TRIM(int2str(ii))//'" ',ii=1,nprobes),&
           ('"Pressp'//TRIM(int2str(ii))//'" ',ii=1,nprobes)
   END IF
end if
!
! Description of volumic informations in a specific file (for coupling with e.g. SWENS approach or velocity/pressure cards)
IF(i_sw.eq.1) then
   OPEN(123,file='Results/modes_HOS_SWENSE.dat',status='REPLACE', FORM='FORMATTED', ACCESS='DIRECT',RECL=18*(n1)) 
   !
   ! We assume that n1 is greater n3_add in the writing
   IF (n3_add.GT.n1) THEN
      PRINT*, 'Change the value of n1 and/or n3 so that n1>=n3_add'
      STOP
   ENDIF
ENDIF
!
END SUBROUTINE init_output
!
!
SUBROUTINE output_time_step(i3d,imodes,iprobes,iwmk,igeom,i_sw,time_cur,t_unit)
!
IMPLICIT NONE
!
! Input variables
INTEGER, INTENT(IN)          :: i3d,imodes,iprobes,iwmk,igeom,i_sw
REAL(RP), INTENT(IN)         :: time_cur
CHARACTER(LEN=1), INTENT(IN) :: t_unit
!
! Local variables
REAL(RP), DIMENSION(m1,m2)     :: FT_tmp, temp1, temp2
REAL(RP)                       :: eta_tmp, eta_target
INTEGER                        :: i1,i2,j,i3
!
! Pressure probes velocities/pressure
REAL(RP), DIMENSION(maxprobes) :: vitx_probe, vity_probe, vitz_probe, Press_probe, eta_probe
!
!	Volume (correct up to 2nd order) 
!   Volume is now corrected in runge_kutta to be exactly the volume variation of wavemaker...
if(igeom.EQ.2) then
  	write(2,201) time_cur * t_adim, voltot - 0.5_rp*alpha*(1-d_hinge)*pos1st(n3,1) / xlen - eta(1,1)*alpha*pos1st(n3,1) / xlen, &
	   enepot, enekin, enetot
elseif(igeom.EQ.1) then
  	write(2,201) time_cur * t_adim, voltot - alpha*pos1st(n3,1) / xlen - eta(1,1)*alpha*pos1st(n3,1) / xlen, &
	   enepot, enekin, enetot
else
  	write(2,201) time_cur * t_adim, voltot, &
	   enepot, enekin, enetot
endif
!
201 format(5(ES12.5,X))
!
! Probes elevation
IF (iprobes.EQ.1) THEN
  DO j = 1, nprobes
	 FT_tmp(1:n1,1:n2) = space_2_Fourier(eta,'cos','cos')
	 eta_tmp = 0.0_rp
	 DO i2=1,n2
	    eta_tmp = eta_tmp + DOT_PRODUCT(FT_tmp(1:n1,i2),mat_coscos(1:n1,i2,j))
	 ENDDO
	 ! Volume is now corrected in runge_kutta to be exactly the volume variation of wavemaker... No adjustment necessary
	 eta_probe(j) = eta_tmp
  ENDDO
  WRITE(9,109) time_cur*t_adim,(eta_probe(j)*x_adim,j=1,nprobes)
  !
  IF (nprobes.GT.60) then
  	print*,'writing format of wave probes has to be changed...'
  	stop
  ENDIF
109   format(61d13.5)
ENDIF
!
! Velocities, pressure and pressure gradient calculation (OMAE2006)
IF (iprobes.eq.2) THEN
  CALL velpress(modesspecx,modesspecy,modesspecz,modesspect,modesFS,modesadd,modesaddt,nprobes,zprobe,yprobe, &
	   vitx_probe,vity_probe,vitz_probe,Press_probe,eta_probe)
ENDIF
!
! Free surface elevations
if (i3d.eq.1) then
  if (ABS(time_cur) < tiny) then
	 IF (icase == 2) THEN
		IF (tecplot == 11) THEN
			WRITE(1,105)'ZONE SOLUTIONTIME = ',time_cur*t_adim,', I=',n1,', J=',n2
		ELSE
		 	WRITE(1,103)'ZONE T = "t = ',time_cur*t_adim, t_unit,', t/T = ',time_cur/T_mono,'", I=',n1,', J=',n2
		END IF
	 ELSE
		IF (tecplot == 11) THEN
			WRITE(1,105)'ZONE SOLUTIONTIME = ',time_cur*t_adim,', I=',n1,', J=',n2
		ELSE
		 	WRITE(1,104)'ZONE T = "t = ',time_cur*t_adim, t_unit, '", I=', n1, ', J=', n2
		END IF
	 END IF
	 DO i2=1,n2
		DO i1=1,n1
		   write(1,101) x(i1,i2)*x_adim,y(i1,i2)*x_adim,eta(i1,i2)*x_adim,phis(i1,i2)*x_adim**2/t_adim
		ENDDO
	 ENDDO
  ELSE
	 IF (icase == 2) THEN
	 	IF (tecplot == 11) THEN
			WRITE(1,105)'ZONE SOLUTIONTIME = ',time_cur*t_adim,', D=(1,2), I=',n1,', J=',n2
		ELSE
		 	WRITE(1,103)'ZONE T = "t = ',time_cur*t_adim, t_unit,', t/T = ',time_cur/T_mono,'", D=(1,2), I=',n1,', J=',n2
		END IF
	 ELSE
	 	IF (tecplot == 11) THEN
			WRITE(1,105)'ZONE SOLUTIONTIME = ',time_cur*t_adim,', D=(1,2), I=',n1,', J=',n2
		ELSE
		 	WRITE(1,104)'ZONE T = "t = ',time_cur*t_adim, t_unit, '", D=(1,2), I=', n1, ', J=', n2
		END IF
	 END IF
	 DO i2=1,n2
		DO i1=1,n1
		   IF (icase == 2) THEN
			  IF (igeom == 1.OR.igeom==2) THEN
				 eta_target = amp_mono * cos(w_mono*time_cur - k_mono*(x(i1,i2)*cos(theta_mono)+y(i1,i2)*sin(theta_mono)) + ph_mono)
			  ELSE IF (igeom==3) THEN
				 eta_target = amp_mono * cos(w_mono*time_cur - k_mono*((x(i1,i2)-x_dip)*cos(theta_mono)+y(i1,i2)*sin(theta_mono)) + ph_mono)
			  END IF
		   ELSE
			  eta_target = 0.0_rp
		   END IF
		   WRITE(1,102) eta(i1,i2)*x_adim,(eta(i1,i2)-eta_target)*x_adim
		ENDDO
	 ENDDO
  ENDIF
ENDIF
101 format(2(F10.5,X),2(ES12.5,X))
102 format(2(ES12.5,X))
103 format(A,F6.2,2A,F5.1,A,I4,A,I4)
104 format(A,F6.2,2A,I4,A,I4)
105 format(A,F6.2,A,I4,A,I4)
!	
! Wavemaker motion
IF (iwmk.EQ.1 .AND. igeom/=3) THEN
  IF (n2 == 1) THEN
     WRITE(8,801) time_cur*t_adim, alpha*pos1st(n3,1)*x_adim , &
     	(alphat*pos1st(n3,1)+alpha*dpos1stdt(n3,1))*x_adim/t_adim
  ELSE
  	IF (tecplot == 11) THEN
		WRITE(8,'(A,F8.2,A,I3)') 'ZONE SOLUTIONTIME = ',time_cur*t_adim,', I=',n2
	ELSE
		 WRITE(8,'(A,F8.2,A,I3)')'ZONE T = "',time_cur*t_adim,'" I = ',n2
	END IF
	DO i2=1,n2
		WRITE(8,801) cy(i2)*x_adim, alpha*pos1st(n3,i2)*x_adim, &
			(alphat*pos1st(n3,i2)+alpha*dpos1stdt(n3,i2))*x_adim/t_adim
	ENDDO
  END IF
801   format(3d13.5)  	
ENDIF
!
! Potential modes
!
if (imodes.eq.1) then
  temp1 = space_2_Fourier(eta*x_adim,'cos','cos')
  temp2 = space_2_Fourier(phis*x_adim**2/t_adim,'cos','cos')
  if (ABS(time_cur) < tiny) then
	 IF (icase == 2) THEN
		IF (tecplot == 11) THEN
			WRITE(5,505)'ZONE SOLUTIONTIME = ',time_cur*t_adim,', I=',n1,', J=',n2
		ELSE
		 	write(5,503)'ZONE T = "t = ',time_cur*t_adim, t_unit,', t/T = ',time_cur/T_mono,'", I=',n1,', J=',n2
		END IF
	 ELSE
		IF (tecplot == 11) THEN
			WRITE(5,505)'ZONE SOLUTIONTIME = ',time_cur*t_adim,', I=',n1,', J=',n2
		ELSE
		 	write(5,504)'ZONE T = "t = ',time_cur*t_adim, t_unit, '", I=', n1, ', J=', n2
		END IF
	 END IF
	 !
	 do i2=1,n2
		do i1=1,n1
		   write(5,501) kx(i1)/x_adim,ky(i2)/x_adim,temp1(i1,i2),temp2(i1,i2), &
		   	LOG10(MAX(tiny,ABS(temp1(i1,i2)))), LOG10(MAX(tiny,ABS(temp2(i1,i2))))
		end do
	 end do
  else
	 IF (icase == 2) THEN
	 	IF (tecplot == 11) THEN
			WRITE(5,505)'ZONE SOLUTIONTIME = ',time_cur*t_adim,', D=(1,2), I=',n1,', J=',n2
		ELSE
		 	write(5,503)'ZONE T = "t = ',time_cur*t_adim, t_unit,', t/T = ',time_cur/T_mono,'", D=(1,2), I=',n1,', J=',n2
		END IF
	 ELSE
	 	IF (tecplot == 11) THEN
			WRITE(5,505)'ZONE SOLUTIONTIME = ',time_cur*t_adim,', D=(1,2), I=',n1,', J=',n2
		ELSE
		 	write(5,504)'ZONE T = "t = ',time_cur*t_adim, t_unit, '", D=(1,2), I=', n1, ', J=', n2
		END IF
	 END IF
	 do i2=1,n2
		do i1=1,n1
	   		write(5,502) temp1(i1,i2),temp2(i1,i2), &
		   	LOG10(MAX(tiny,ABS(temp1(i1,i2)))), LOG10(MAX(tiny,ABS(temp2(i1,i2))))
		end do
	 end do
  end if
end if
!
! Spectral Modes Needed in SWENSE (JC, Romain, PM)
!
IF(i_sw.EQ.1) THEN
! 10/2013... remove to test influence FIXME
!    modesFS = space_2_Fourier(eta,'cos','cos')
!      ! GD: modif niveau moyen (07/2012)
!      if(igeom.EQ.2) then
!         modesFS = space_2_Fourier(eta - (voltot-0.5_rp*alpha*(1-d_hinge)*pos1st(n3,1) / xlen - eta(1,1)*alpha*pos1st(n3,1) / xlen),'cos','cos')
!      elseif(igeom.EQ.1) then
!         modesFS = space_2_Fourier(eta - (voltot-alpha*pos1st(n3,1) / xlen - eta(1,1)*alpha*pos1st(n3,1) / xlen),'cos','cos')
!      else
!         modesFS = space_2_Fourier(eta - (voltot),'cos','cos')
!      endif
  if (ABS(time_cur) < tiny) then
	 WRITE(123,'(5000(ES17.10,1X))',REC=1) 1.0_rp*n1, 1.0_rp*n3_add , 1.0_rp*n2, 1.0_rp/f_out,T_stop &
		  ,xlen , ylen, h, (modesspecx(i1,1), i1=9,n1)
	 
	 WRITE(123,'(5000(ES17.10,1X))',REC=2) (modesspecy(i1,1), i1=1,n1)
	 WRITE(123,'(5000(ES17.10,1X))',REC=3) (modesspecz(i1,1), i1=1,n1)
	 WRITE(123,'(5000(ES17.10,1X))',REC=4) (modesspect(i1,1), i1=1,n1)
	 WRITE(123,'(5000(ES17.10,1X))',REC=5) (modesFS(i1,1)   , i1=1,n1)
	 WRITE(123,'(5000(ES17.10,1X))',REC=6) (modesFSt(i1,1)  , i1=1,n1)
	 WRITE(123,'(5000(ES17.10,1X))',REC=7) (modesadd(i3,1)  , i3=1,n3_add)
	 WRITE(123,'(5000(ES17.10,1X))',REC=8) (modesaddt(i3,1) , i3=1,n3_add)
	 do i2=2,n2
		WRITE(123,'(5000(ES17.10,1X))',REC=1+8*(i2-1)) (modesspecx(i1,i2), i1=1,n1)
		WRITE(123,'(5000(ES17.10,1X))',REC=2+8*(i2-1)) (modesspecy(i1,i2), i1=1,n1)
		WRITE(123,'(5000(ES17.10,1X))',REC=3+8*(i2-1)) (modesspecz(i1,i2), i1=1,n1)
		WRITE(123,'(5000(ES17.10,1X))',REC=4+8*(i2-1)) (modesspect(i1,i2), i1=1,n1)
		WRITE(123,'(5000(ES17.10,1X))',REC=5+8*(i2-1)) (modesFS(i1,i2)   , i1=1,n1)
		WRITE(123,'(5000(ES17.10,1X))',REC=6+8*(i2-1)) (modesFSt(i1,i2)  , i1=1,n1)
		WRITE(123,'(5000(ES17.10,1X))',REC=7+8*(i2-1)) (modesadd(i3,i2)  , i3=1,n3_add)
		WRITE(123,'(5000(ES17.10,1X))',REC=8+8*(i2-1)) (modesaddt(i3,i2) , i3=1,n3_add)
	 enddo
  endif
  do i2=1,n2
	 WRITE(123,'(5000(ES17.10,1X))',REC=((it-1)*n2*8)+1+8*(i2-1)) (modesspecx(i1,i2), i1=1,n1)
	 WRITE(123,'(5000(ES17.10,1X))',REC=((it-1)*n2*8)+2+8*(i2-1)) (modesspecy(i1,i2), i1=1,n1)
	 WRITE(123,'(5000(ES17.10,1X))',REC=((it-1)*n2*8)+3+8*(i2-1)) (modesspecz(i1,i2), i1=1,n1)
	 WRITE(123,'(5000(ES17.10,1X))',REC=((it-1)*n2*8)+4+8*(i2-1)) (modesspect(i1,i2), i1=1,n1)
	 WRITE(123,'(5000(ES17.10,1X))',REC=((it-1)*n2*8)+5+8*(i2-1)) (modesFS(i1,i2)   , i1=1,n1)
	 WRITE(123,'(5000(ES17.10,1X))',REC=((it-1)*n2*8)+6+8*(i2-1)) (modesFSt(i1,i2)  , i1=1,n1)
	 WRITE(123,'(5000(ES17.10,1X))',REC=((it-1)*n2*8)+7+8*(i2-1)) (modesadd(i3,i2)  , i3=1,n3_add)
	 WRITE(123,'(5000(ES17.10,1X))',REC=((it-1)*n2*8)+8+8*(i2-1)) (modesaddt(i3,i2) , i3=1,n3_add)
  enddo
  !
501   format(6(ES12.5,X))
502   format(4(ES12.5,X))
503   format(A,F6.2,2A,F5.1,A,I4,A,I4)
504   format(A,F6.2,2A,I4,A,I4)
505   format(A,F6.2,A,I4,A,I4)
  !
ENDIF
!
!	Velocities and pressure on probes
if (iprobes.eq.2) THEN
  if (nprobes.GT.60) then
  	print*,'writing format of pressure probes has to be changed...'
  	stop
  endif
  ! writing format should be at least 5*maxprobes+1
  write(99,'(301(ES16.9,X))') time_cur*t_adim,(eta_probe(j)*x_adim,j=1,nprobes),&
       (vitx_probe(j)*x_adim/t_adim,j=1,nprobes),&
	   (vity_probe(j)*x_adim/t_adim,j=1,nprobes), &
	   (vitz_probe(j)*x_adim/t_adim,j=1,nprobes), (Press_probe(j)*(x_adim**2)/(t_adim**2),j=1,nprobes)
end if
!
END SUBROUTINE output_time_step
!
!
!
END MODULE output