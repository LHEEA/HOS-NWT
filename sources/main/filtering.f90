MODULE filtering
!
! This module gives the spatial filters that may be applied in HOS-NWT simulations
! Filtering may be applied in x,y and/or z directions independently
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
USE Fourier_FFTW
!
CONTAINS
!
!
! start filtering_x ******************************************************
!           
!     ====================================
SUBROUTINE filtering_x(filtered,coeff_cut)
!     ====================================
!
IMPLICIT NONE
!
!% INPUT VARIABLES
!	quantity to filter before and after
REAL(RP), DIMENSION(m1,m2) :: filtered
REAL(RP) :: coeff_cut
!
!% LOCAL VARIABLES
!	temporary quantity to be filtered
REAL(RP), DIMENSION(m1,m2) :: as
INTEGER :: n1filt, i1, i2
!
! This subroutines filters quantity filtered with filtering defined by coeff_cut
n1filt=ABS(int(n1*coeff_cut))
!
if(coeff_cut.ge.0.9999d0) n1filt=n1
!     
! inverse FFT: determination of the modes of 'filtered'
do i1=1,n1
   do i2=1,n2
      as(i1,i2)=filtered(i1,i2) 
   end do
end do
!
as = space_2_Fourier(as,'cos','cos')
!
! Filtering along first direction
CALL filtering_base(as(1:n1,1:n2),n1,n2,n1filt,1,coeff_cut)
!
! FFT: determination of the 'filtered' 
as = Fourier_2_space(as,'cos','cos')
!
do i1=1,n1
   do i2=1,n2
      filtered(i1,i2)=as(i1,i2) 
   end do
end do
!
!
END SUBROUTINE  filtering_x
!
! end filtering_x *******************************************************
!
!
!
! beginning filtering_y ******************************************************
!           
!     ====================================
SUBROUTINE filtering_y(filtered,coeff_cut)
!     ====================================
!
IMPLICIT NONE
!
!% INPUT VARIABLES
!	quantity to filter before and after
!
REAL(RP), DIMENSION(m1,m2) :: filtered
REAL(RP) :: coeff_cut
!
!% LOCAL VARIABLES
!	temporary quantity to be filtered
!
REAL(RP), DIMENSION(m1,m2) :: as
INTEGER :: n2filt, i1, i2
!
! This subroutines filters quantity filtered with filtering defined by coeff_cut
if (n2.ne.1) then
   n2filt=ABS(int(n2*coeff_cut))
   if(coeff_cut.ge.0.9999d0) n2filt=n2
else
   n2filt=1
endif
!      
! inverse FFT: determination of the modes of 'filtered'
do i1=1,n1
   do i2=1,n2
      as(i1,i2)=filtered(i1,i2) 
   end do
end do
!
as = space_2_Fourier(as,'cos','cos')
!
! Filtering along second direction
CALL filtering_base(as(1:n1,1:n2),n1,n2,n2filt,2,coeff_cut)
!
! FFT: determination of the 'filtered' 
as = Fourier_2_space(as,'cos','cos')
!
do i1=1,n1
   do i2=1,n2
      filtered(i1,i2)=as(i1,i2) 
   end do
end do

END SUBROUTINE filtering_y 
!
!  end filtering_y *******************************************************
!
!
!
! beginning filtering_z ******************************************************
!           
!     ====================================
SUBROUTINE filtering_z(filtered,coeff_cut)
!     ====================================

IMPLICIT NONE
!
!% INPUT VARIABLES
!	quantity to filter before and after
!
REAL(RP), DIMENSION(m3_add,m2) :: filtered
REAL(RP) :: coeff_cut
!
!% LOCAL VARIABLES
!	temporary quantity to be filtered
!
REAL(RP), DIMENSION(m3_add,m2) :: as
INTEGER :: n3_add_filt, i2, i3
!
! This subroutines filters quantity filtered with filtering defined by coeff_cut
n3_add_filt=ABS(int(n3_add*coeff_cut))
if(coeff_cut.ge.0.9999d0) n3_add_filt=n3_add
!      
! inverse FFT: determination of the modes of 'filtered'
do i3=1,n3_add
   do i2=1,n2
      as(i3,i2) = filtered(i3,i2)
   end do
end do
!
as = space_2_Fourier_add(as,'cos','cos')
!
! Filtering along first direction
CALL filtering_base(as(1:n3_add,1:n2),n3_add,n2,n3_add_filt,1,coeff_cut)
!
! FFT: determination of the 'filtered' 
as = Fourier_2_space_add(as,'cos','cos')
!
do i3=1,n3_add
   do i2=1,n2
      filtered(i3,i2) = as(i3,i2)
   end do
end do
!
END SUBROUTINE filtering_z 
!
!  end filtering_z *******************************************************
!
!
!
! beginning filtering_x_modes ******************************************************
!           
!     ====================================
SUBROUTINE filtering_x_modes(filtered,coeff_cut)
!     ====================================
!
! Input and output are modal description
IMPLICIT NONE
!
!% INPUT VARIABLES
!	quantity to filter before and after
REAL(RP), DIMENSION(m1,m2) :: filtered
REAL(RP) :: coeff_cut
!
!% LOCAL VARIABLES
!	temporary quantity to be filtered
REAL(RP), DIMENSION(m1,m2) :: as
INTEGER :: n1filt, i1, i2
!
! This subroutines filters quantity filtered with filtering defined by coeff_cut
n1filt=ABS(int(n1*coeff_cut))
if(coeff_cut.ge.0.9999d0) n1filt=n1
!     
do i1=1,n1
   do i2=1,n2
      as(i1,i2)=filtered(i1,i2) 
   end do
end do
!
! Filtering along first direction
CALL filtering_base(as(1:n1,1:n2),n1,n2,n1filt,1,coeff_cut)
!
do i1=1,n1
   do i2=1,n2
      filtered(i1,i2)=as(i1,i2) 
   end do
end do
!
!
END SUBROUTINE  filtering_x_modes
!
! end filtering_x_modes *******************************************************
!
!
!
! beginning filtering_y_modes ******************************************************
!           
!     ====================================
SUBROUTINE filtering_y_modes(filtered,coeff_cut)
!     ====================================
!
! Input and output are modal description
IMPLICIT NONE
!
!% INPUT VARIABLES
!	quantity to filter before and after
!
REAL(RP), DIMENSION(m1,m2) :: filtered
REAL(RP) :: coeff_cut
!
!% LOCAL VARIABLES
!	temporary quantity to be filtered
!
REAL(RP), DIMENSION(m1,m2) :: as
INTEGER :: n2filt, i1, i2
!
! This subroutines filters quantity filtered with filtering defined by coeff_cut
if (n2.ne.1) then
   n2filt=ABS(int(n2*coeff_cut))
   if(coeff_cut.ge.0.9999d0) n2filt=n2
else
   n2filt=1
endif
!      
do i1=1,n1
   do i2=1,n2
      as(i1,i2)=filtered(i1,i2) 
   end do
end do
!
! Filtering along second direction
CALL filtering_base(as(1:n1,1:n2),n1,n2,n2filt,2,coeff_cut)
!
do i1=1,n1
   do i2=1,n2
      filtered(i1,i2)=as(i1,i2) 
   end do
end do

END SUBROUTINE filtering_y_modes
!
!  end filtering_y_modes *******************************************************
!
!
!
! beginning filtering_z_modes ******************************************************
!           
!     ====================================
SUBROUTINE filtering_z_modes(filtered,coeff_cut)
!     ====================================
!
! Input and output are modal description
IMPLICIT NONE
!
!% INPUT VARIABLES
!	quantity to filter before and after
!
REAL(RP), DIMENSION(m3_add,m2) :: filtered
REAL(RP) :: coeff_cut
!
!% LOCAL VARIABLES
!	temporary quantity to be filtered
!
REAL(RP), DIMENSION(m3_add,m2) :: as
INTEGER :: n3_add_filt, i2, i3
!
! This subroutines filters quantity filtered with filtering defined by coeff_cut
n3_add_filt=ABS(int(n3_add*coeff_cut))
if(coeff_cut.ge.0.9999d0) n3_add_filt=n3_add
!      
do i3=1,n3_add
   do i2=1,n2
      as(i3,i2) = filtered(i3,i2)
   end do
end do
!
! Filtering along first direction
CALL filtering_base(as(1:n3_add,1:n2),n3_add,n2,n3_add_filt,1,coeff_cut)
!
do i3=1,n3_add
   do i2=1,n2
      filtered(i3,i2) = as(i3,i2)
   end do
end do
!
END SUBROUTINE filtering_z_modes
!
!  end filtering_z_modes *******************************************************
!
!
! start filtering_base ******************************************************
!           
!     ====================================
SUBROUTINE filtering_base(as,n1,n2,nfilt,filt_1or2,coeff_cut)
!     ====================================
!
IMPLICIT NONE
!
!% INPUT VARIABLES
!	modal description of quantity to filter
INTEGER, INTENT(IN)                       :: n1,n2,nfilt,filt_1or2
REAL(RP), DIMENSION(n1,n2), INTENT(INOUT) :: as
REAL(RP) :: coeff_cut
!
!% LOCAL VARIABLES
INTEGER  :: i1, i2, p
REAL(RP) :: sigma0,alpha_cut,ERFC
!
! This subroutines defines the different filtering types and apply filter on input
! which has to be a modal description
!
alpha_cut= 2.d0
IF(coeff_cut.GE.0.d0) THEN
   ! 
	!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Initial filter
	!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! modes over n1filt or n2filt set to zero
   !
   IF (filt_1or2 == 1) THEN
	   do i1=nfilt+1,n1
		  do i2=1,n2
			 as(i1,i2)=0.d0
		  end do
	   end do
	ELSEIF (filt_1or2 == 2) THEN
		do i1=1,n1
		  do i2=nfilt+1,n2
			 as(i1,i2)=0.d0
		  end do
	   end do
	ELSE
		print*, 'Problem in filtering choice'
		STOP
	ENDIF
   !
ELSEIF(ABS(coeff_cut+1.0_rp) < tiny) THEN
   !
	!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Cesaro filter
	!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   IF (filt_1or2 == 1) THEN
	   do i1=1,n1
		  do i2=1,n2
			 as(i1,i2)=as(i1,i2)*(1.d0-REAL((i1)/(n1)))
		  end do
	   end do
   ELSEIF (filt_1or2 == 2) THEN
	   do i1=1,n1
		  do i2=1,n2
			 as(i1,i2)=as(i1,i2)*(1.d0-REAL((i2)/(n2)))
		  end do
	   end do   
   ELSE
		print*, 'Problem in filtering choice'
		STOP
	ENDIF
   !
ELSEIF(ABS(coeff_cut+2.0_rp) < tiny) THEN
   !
	!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Lanczos filter
	!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   IF (filt_1or2 == 1) THEN
	   do i1=1,n1
		  do i2=1,n2
			 as(i1,i2)=as(i1,i2)*SIN(TWOPI*i1/n1)/(TWOPI*i1/n1)
		  end do
	   end do
   ELSEIF (filt_1or2 == 2) THEN
	   do i1=1,n1
		  do i2=1,n2
			 as(i1,i2)=as(i1,i2)*SIN(TWOPI*i2/n2)/(TWOPI*i2/n2)
		  end do
	   end do
   ELSE
		print*, 'Problem in filtering choice'
		STOP
	ENDIF
   !
ELSEIF(ABS(coeff_cut+3.0_rp) < tiny) THEN
   !
	!!!!!!!!!!!!!!!!!!!!!!!!!
   ! raised cosine filter
	!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   IF (filt_1or2 == 1) THEN
	   do i1=1,n1
		  do i2=1,n2
			 as(i1,i2)=as(i1,i2)*(1.d0+COS(TWOPI*i1/n1))/2.d0
		  end do
	   end do
   ELSEIF (filt_1or2 == 2) THEN
	   do i1=1,n1
		  do i2=1,n2
			 as(i1,i2)=as(i1,i2)*(1.d0+COS(TWOPI*i2/n2))/2.d0
		  end do
	   end do
   ELSE
		print*, 'Problem in filtering choice'
		STOP
	ENDIF
   !
ELSEIF(ABS(coeff_cut+4.0_rp) < tiny) THEN
   !
	!!!!!!!!!!!!!!!!!!!!!!!!!
   ! sharpened raised cosine filter
	!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   IF (filt_1or2 == 1) THEN
	   do i1=1,n1
		  do i2=1,n2
			 sigma0 = (1.d0+COS(TWOPI*i1/n1))/2.d0
			 as(i1,i2)=as(i1,i2)*sigma0**4*(35.d0-84.d0*sigma0+70.d0*sigma0**2-20.d0*sigma0**3)
		  end do
	   end do
   ELSEIF (filt_1or2 == 2) THEN
	   do i1=1,n1
		  do i2=1,n2
			 sigma0 = (1.d0+COS(TWOPI*i2/n2))/2.d0
			 as(i1,i2)=as(i1,i2)*sigma0**4*(35.d0-84.d0*sigma0+70.d0*sigma0**2-20.d0*sigma0**3)
		  end do
	   end do   
   ELSE
		print*, 'Problem in filtering choice'
		STOP
	ENDIF
   !
ELSEIF(ABS(coeff_cut+5.0_rp) < tiny) THEN
   !
	!!!!!!!!!!!!!!!!!!!!!!!!!
   ! erfc-log filter
	!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   IF (filt_1or2 == 1) THEN
	   do i1=1,n1
		  do i2=1,n2
			 p = 10
			 sigma0 = (REAL(i1)-3.d0/4.d0*REAL(n1))/(REAL(n1)+1.d0)
			 as(i1,i2)=as(i1,i2)*1.d0/2.d0*ERFC(2.d0*REAL(p)**0.5d0*sigma0 &
				  *SQRT((-LOG(1.d0-(4.d0/3.d0*sigma0)**2))/((4.d0/3.d0*sigma0)**2)))
		  end do
	   end do
   ELSEIF (filt_1or2 == 2) THEN
	   do i1=1,n1
		  do i2=1,n2
			 p = 10
			 sigma0 = (REAL(i2)-3.d0/4.d0*REAL(n2))/(REAL(n2)+1.d0)
			 as(i1,i2)=as(i1,i2)*1.d0/2.d0*ERFC(2.d0*REAL(p)**0.5d0*sigma0 &
				  *SQRT((-LOG(1.d0-(4.d0/3.d0*sigma0)**2))/((4.d0/3.d0*sigma0)**2)))
		  end do
	   end do   
   ELSE
		print*, 'Problem in filtering choice'
		STOP
	ENDIF
ELSE
   !
	!!!!!!!!!!!!!!!!!!!!!!!!!
   ! exponential cut-off filter
	!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   IF (filt_1or2 == 1) THEN
	   do i1=nfilt+1,n1
		  do i2=1,n2
			 as(i1,i2)=as(i1,i2)*exp(-alpha_cut*(TWOPI*REAL(i1/n1)-TWOPI*REAL(nfilt/n1))**4)
		  end do
	   end do
   ELSEIF (filt_1or2 == 2) THEN
	   do i1=n1,n1
		  do i2=nfilt+1,n2
			 as(i1,i2)=as(i1,i2)*exp(-alpha_cut*(TWOPI*REAL(i2/n2)-TWOPI*REAL(nfilt/n2))**4)
		  end do
	   end do
   ELSE
		print*, 'Problem in filtering choice'
		STOP
	ENDIF
ENDIF
!
!
END SUBROUTINE  filtering_base
!
! end filtering_base *******************************************************
!
!
! end module
END MODULE filtering
