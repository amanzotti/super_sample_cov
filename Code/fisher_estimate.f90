program fisher_estimate

  USE healpix_types 

	implicit none
	integer(i4b)::i,jlo,l
	integer,parameter :: lmax=3000
	real(dp), parameter ::fsky=1d0,thetaFWHM=1d0*0.000290888209d0
	real(dp), allocatable, dimension(:,:)::Cl
	real(dp), allocatable, dimension(:,:,:)::CovMatrix,CovMatrixInv

	real(dp), allocatable, dimension(:):: xa,ya,y2a,deltaCl,deltaCl_pol,y_pol,y_pol2,x_pol !! variable for spline the TT spectrum and EE
	real(dp), allocatable, dimension(:):: xTE,yTE,yTE2,deltaClTE!! variable for spline the TE spectrum
    LOGICAL ::OK_FLAG	
	character(len=128) :: filename
	integer  :: ios,err,istat
	real(dp)::logCl,logCl_pol,fisher,fisher_pol,fisherTE, variance,variance_pol
	filename='/home/manzotti/Software/camb/test_scalCls.dat'
        



write(*,*),filename 
	open(unit=10, file=filename, iostat=ios, &
     status="old", action="read", form='formatted')
	if ( ios /= 0 ) stop "Error opening file name"
  write(*,*) 'here'

	allocate(CovMatrix(3,3,lmax),CovMatrixInv(3,3,lmax), stat=err)
	if (err /= 0) print *, "array: Allocation request denied"
	allocate(cl(lmax,4), stat=err)
	if (err /= 0) print *, "array: Allocation request denied"
	allocate(xa(lmax),ya(lmax),y2a(lmax),deltaCl(lmax),y_pol(lmax)&
             ,y_pol2(lmax),deltaCl_pol(lmax),yTE(lmax),yTE2(lmax),deltaClTE(lmax),stat=err)
	if (err /= 0) print *, "array: Allocation request denied"
	!! read the Cl from camb


	do i = 1, lmax, 1

		read(10, *) cl(i,1),cl(i,2),cl(i,3),cl(i,4)
		if ( istat /= 0 ) stop "Read error in file 10"				
	end do


!	write(*,*) cl(10,1),cl(10,2)

!! normalize CAMB
	do i = 1, lmax, 1
		l=cl(i,1)
		cl(i,2)=cl(i,2)*(twopi/(l*(l+1)))!/(7.4311e12)
		cl(i,3)=cl(i,3)*(twopi/(l*(l+1)))!/(7.4311e12)
		cl(i,4)=cl(i,4)*(twopi/(l*(l+1)))!/(7.4311e12)

	end do	
	write(*,*) cl(10,1),cl(10,2),cl(10,3),cl(10,4)

	
	!!! spline the TT spectrum

	xa=log(cl(:,1)) !log(l)
	ya=log(cl(:,1)**2*cl(:,2)) !log (l^2 Cl)
	call spline(xa,ya,lmax,1.d32,1.d32,y2a)

!!! spline the EE spectrum

	xa=log(cl(:,1)) !log(l)
	y_pol=log(cl(:,1)**2*cl(:,3)) !log (l^2 CEE)
 !       write(*,*) shape(y_pol),shape(xa),shape(y_pol2)
	call spline(xa,y_pol,lmax,1.d32,1.d32,y_pol2)


!!! spline the TE spectrum
!!! note: numerical problematic cause it change SIGN. This is not implemented directly in a logaritmic fashion

	xa=log(cl(:,1)) !log(l)
	yTE=(cl(:,1))**2*(cl(:,4)) !log (l^2 CEE)
!        write(*,*) shape(yTE),shape(xa),shape(yTE2)
	call spline(xa,yTE,lmax,1.d32,1.d32,yTE2)

!!write(*,*) yTE,y_pol2

!stop

!	write(*,*) shape(y2a),y2a(10)
	write(*,*) 'power spectrum TT&EE&TE splined'
!	logCl=logl2Cl(10._dp)
!	write(*,*) 'power spectrum splined'!,logl2Cl(11._dp),dlnl2Cldlogl(11d0)

	!!! prepare for fishere elements, compute deltaCl

	!!! deltaCl=sqrt(2/(2l+1)fsky) (Cl+Nl) we need delta cl squared


        !! COMPUTE OUTPUT ANS PRINT

        !cl(i,2) C^T
        !cl(i,3) C^E
        !cl(i,4) C^TE

do i=2,lmax,1
  CovMatrix(1,1,i)=cl(i,2)**2 !! TT TT variance look 9807130v2
  CovMatrix(1,2,i)=cl(i,2)*cl(i,4) !! TT TE variance
  CovMatrix(2,1,i)=CovMatrix(1,2,i)
  CovMatrix(1,3,i)=cl(i,4)**2 !! TT EE variance
  CovMatrix(3,1,i)=CovMatrix(1,3,i)
  CovMatrix(2,2,i)=(cl(i,4)**2+cl(i,2)*cl(i,3)) !! TE TE variance
  CovMatrix(2,3,i)=(cl(i,3)*cl(i,4)) !! TE EE 
  CovMatrix(3,2,i)=CovMatrix(2,3,i)
  CovMatrix(3,3,i)=cl(i,3)**2 !! EE EE
  CovMatrix(:,:,i)=2/(2*cl(i,1)+1)*CovMatrix(:,:,i)
enddo
 !! EE EE 
open(unit=30, file='variancenonoisenofsky_fisher.txt')

!! fisher elements  sum_l 1/deltaC (deriv)^2
	fisher=0d0
        fisher_pol=0d0
        fisherTE=0d0
            write(*,*) 'fisher', fisher

	do i = 100, lmax, 1

		!! invert the 3 by 3 covariance matrix
write(*,*) 'i', i
write(*,*) CovMatrix(1,:,i)
write(*,*) CovMatrix(2,:,i)
write(*,*) CovMatrix(3,:,i)
write(*,*) ''
		call m33inv(CovMatrix(:,:,i),CovMatrixInv(:,:,i),OK_FLAG)
if ( OK_FLAG .neqv. .true. ) stop "Something wrong while inverting the matrix"


write(*,*) CovMatrixInv(1,:,i)
write(*,*) CovMatrixInv(2,:,i)
write(*,*) CovMatrixInv(3,:,i)
write(*,*) ''



	fisher = fisher+cl(i,2)*cl(i,2)*(dlnl2Cldlogl2(DBLE(i),ya,y2a))*(dlnl2Cldlogl2(DBLE(i),ya,y2a))*CovMatrixInv(1,1,i)& !TTTT
	+(cl(i,3)*cl(i,3)*(dlnl2Cldlogl2(DBLE(i),y_pol,y_pol2))*(dlnl2Cldlogl2(DBLE(i),y_pol,y_pol2))*CovMatrixInv(3,3,i))& !EEEE
	+(cl(i,4)*cl(i,4)*(dlnl2Cldlogl2(DBLE(i),yTE,yTE2))/((cl(i-1,1))**2*(cl(i-1,4)))*&
    (dlnl2Cldlogl2(DBLE(i),yTE,yTE2))/((cl(i-1,1))**2*(cl(i-1,4)))*CovMatrixInv(2,2,i))& !TETE
	+(cl(i,4)*cl(i,2)*(dlnl2Cldlogl2(DBLE(i),yTE,yTE2))/((cl(i-1,1))**2*(cl(i-1,4)))&
    *(dlnl2Cldlogl2(DBLE(i),ya,y2a))*CovMatrixInv(1,2,i))& !TETT
	+(cl(i,3)*cl(i,2)*(dlnl2Cldlogl2(DBLE(i),y_pol,y_pol2))*(dlnl2Cldlogl2(DBLE(i),ya,y2a))*CovMatrixInv(1,3,i))& !EETT
	+(cl(i,3)*cl(i,4)*(dlnl2Cldlogl2(DBLE(i),y_pol,y_pol2))*(dlnl2Cldlogl2(DBLE(i),yTE,yTE2))&
    /((cl(i-1,1))**2*(cl(i-1,4)))*CovMatrixInv(2,3,i)) !EETE

write(*,*) 'fisher', fisher                

		fisher_pol=fisher_pol+(1/(2d0/((2*cl(i,1)+1))))*(dlnl2Cldlogl2(DBLE(i),y_pol,y_pol2))**2!*(cl(i,3))**2
		fisherTE=fisherTE+(1/(2d0/((2*cl(i,1)+1))))*(dlnl2Cldlogl2(DBLE(i),ya,y2a))**2!*(cl(i,4))**2

		write(30,*) cl(i,1),1/fisher,1/fisher_pol,1/fisherTE

	end do
        close(30)
!! now get variance of the parameter from  fisher



	open(unit=20, file='derivatives2.txt')

	do i = 1, lmax, 1

		write(20,*) cl(i,1),(dlnl2Cldlogl2(DBLE(i),ya,y2a)),logl2Cl(dble(i)),&
        (dlnl2Cldlogl2(DBLE(i),yTE,yTE2))/(cl(i-1,1)**2*cl(i-1,4)),logl2Cl_pol(dble(i)),(dlnl2Cldlogl2(DBLE(i),y_pol,y_pol2))
	end do
	close(unit=20)


	variance=(1/fisher)

	write(*,*) 'your variance is sigma=',variance

	if (allocated(cl)) deallocate(cl, stat=err)
	if (err /= 0) print *, "array: Deallocation request denied"

	

	Contains

		!!! the two functions to get lnl2Cl from spline

real(DP) FUNCTION logl2Cl(ak)
  IMPLICIT none
  real(dp) :: a,b,h,x,y,ak
  x  = log(ak)  
  CALL hunt(xa,lmax,x,jlo)
  h=xa(jlo+1)-xa(jlo)
  a=(xa(jlo+1)-x)/h
  b=(x-xa(jlo))/h
  y=a*ya(jlo)+b*ya(jlo+1)+((a**3-a)*y2a(jlo)+(b**3-b)*y2a(jlo+1))*(h**2)/6.
  logl2Cl = y
  100 return
END FUNCTION logl2Cl


	real(DP) FUNCTION logl2Cl_pol(ak)
  IMPLICIT none
  real(dp) :: a,b,h,x,y,ak
  x  = log(ak)  
  CALL hunt(xa,lmax,x,jlo)
  h=xa(jlo+1)-xa(jlo)
  a=(xa(jlo+1)-x)/h
  b=(x-xa(jlo))/h
  y=a*y_pol(jlo)+b*y_pol(jlo+1)+((a**3-a)*y_pol2(jlo)+(b**3-b)*y_pol2(jlo+1))*(h**2)/6.
  logl2Cl_pol = y
  100 return
END FUNCTION logl2Cl_pol

DOUBLE PRECISION FUNCTION dlnl2Cldlogl(ak)
  IMPLICIT none
  DOUBLE PRECISION :: a,b,h,x,y,ak
  x  = log(ak)
  CALL hunt(xa,lmax,x,jlo)
  h=xa(jlo+1)-xa(jlo)
  a=(xa(jlo+1)-x)/h
  b=(x-xa(jlo))/h
  y=(ya(jlo+1)-ya(jlo))/h+(-(3.*a**2-1.)*y2a(jlo)+(3.*b**2-1.)*y2a(jlo+1))*h/6.
  dlnl2Cldlogl = y
  return
	END FUNCTION dlnl2Cldlogl



DOUBLE PRECISION FUNCTION dlnl2Cldlogl2(ak,ya,ya2)
  IMPLICIT none
  DOUBLE PRECISION :: a,b,h,x,y,ak
  Double precision, dimension(:)::ya2,ya
  x  =  log(ak)
  CALL hunt(xa,lmax,x,jlo)
  h=xa(jlo+1)-xa(jlo)
  a=(xa(jlo+1)-x)/h
  b=(x-xa(jlo))/h
  y=(ya(jlo+1)-ya(jlo))/h+(-(3.*a**2-1.)*ya2(jlo)+(3.*b**2-1.)*ya2(jlo+1))*h/6.
  dlnl2Cldlogl2 = y
  return
END FUNCTION dlnl2Cldlogl2


 FUNCTION M33DET (A) RESULT (DET)

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN)  :: A(:,:)

      DOUBLE PRECISION :: DET

if (size(A).ne.9) then
         write(*,*) 'Dimension of the matrix- function do not agree'
         STOP  
      END IF

      DET =   A(1,1)*A(2,2)*A(3,3)  &
            - A(1,1)*A(2,3)*A(3,2)  &
            - A(1,2)*A(2,1)*A(3,3)  &
            + A(1,2)*A(2,3)*A(3,1)  &
            + A(1,3)*A(2,1)*A(3,2)  &
            - A(1,3)*A(2,2)*A(3,1)

      RETURN

      END FUNCTION M33DET


SUBROUTINE M33INV (A, AINV, OK_FLAG)

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN)  :: A(:,:)
      DOUBLE PRECISION, INTENT(OUT)  :: AINV(:,:)
      LOGICAL, INTENT(OUT) :: OK_FLAG

      DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-50
      DOUBLE PRECISION :: DET
      DOUBLE PRECISION, DIMENSION(3,3) :: COFACTOR
      
      if (size(A).ne.9) then
         write(*,*) 'Dimension of the matrix- function do not agree'
           STOP  
        END IF

      DET =   A(1,1)*A(2,2)*A(3,3)  &
            - A(1,1)*A(2,3)*A(3,2)  &
            - A(1,2)*A(2,1)*A(3,3)  &
            + A(1,2)*A(2,3)*A(3,1)  &
            + A(1,3)*A(2,1)*A(3,2)  &
            - A(1,3)*A(2,2)*A(3,1)
write(*,*) 'det', det
      IF (ABS(DET) .LE. EPS) THEN
         AINV = 0.0D0
         OK_FLAG = .FALSE.
         RETURN
      END IF

      COFACTOR(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
      COFACTOR(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
      COFACTOR(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
      COFACTOR(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
      COFACTOR(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
      COFACTOR(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
      COFACTOR(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
      COFACTOR(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
      COFACTOR(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))

      AINV = TRANSPOSE(COFACTOR) / DET

      OK_FLAG = .TRUE.

      RETURN

      END SUBROUTINE M33INV



	end program fisher_estimate



