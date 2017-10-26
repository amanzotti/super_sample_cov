program fisher_estimate

  USE healpix_types 

	implicit none
	integer(i4b)::i,jlo,l
	integer,parameter :: lmax=3500
	real(dp), parameter ::fsky=1d0,thetaFWHM=1d0*0.000290888209d0
	real(dp), allocatable, dimension(:,:)::Cl
	real(dp), allocatable, dimension(:):: xa,ya,y2a,deltaCl,deltaCl_pol,y_pol,y_pol2,x_pol !! variable for spline the TT spectrum and EE
	character(len=128) :: filename
	integer  :: ios,err,istat
	real(dp)::logCl,logCl_pol,fisher,fisher_pol, variance,variance_pol
	filename='/Users/alessandromanzotti/Downloads/camb/test_scal2Cls.dat'





	open(unit=10, file=filename, iostat=ios, &
             status="old", action="read", form='formatted')
        
	if ( ios /= 0 ) stop "Error opening file name"
	allocate(cl(lmax,4), stat=err)
	if (err /= 0) print *, "array: Allocation request denied"
	allocate(xa(lmax),ya(lmax),y2a(lmax),deltaCl(lmax),y_pol(lmax),y_pol2(lmax),deltaCl_pol(lmax) ,stat=err)
	if (err /= 0) print *, "array: Allocation request denied"
	!! read the Cl from camb


	do i = 1, lmax, 1

		read(10, *) cl(i,1),cl(i,2),cl(i,3)
		if ( istat /= 0 ) stop "Read error in file 10"
		
		
	end do



!	write(*,*) cl(10,1),cl(10,2)

!! normalize CAMB
	do i = 1, lmax, 1
		l=cl(i,1)
		cl(i,2)=cl(i,2)*(twopi/(l*(l+1)))/(7.4311e12)
		cl(i,3)=cl(i,3)*(twopi/(l*(l+1)))/(7.4311e12)

	end do	
	write(*,*) cl(10,1),cl(10,2),cl(10,3)

	

	!!! spline the TT spectrum

	xa=log(cl(:,1)) !log(l)
	ya=log(cl(:,1)**2*cl(:,2)) !log (l^2 Cl)
	call spline(xa,ya,lmax,1.d30,1.d30,y2a)

!!! spline the EE spectrum

	xa=log(cl(:,1)) !log(l)
	y_pol=log(cl(:,1)**2*cl(:,3)) !log (l^2 CEE)
        write(*,*) shape(y_pol),shape(xa),shape(y_pol2)
	call spline(xa,y_pol,lmax,1.d30,1.d30,y_pol2)

!	write(*,*) shape(y2a),y2a(10)
	write(*,*) 'power spectrum TT&EE splined'
!	logCl=logl2Cl(10._dp)
!	write(*,*) 'power spectrum splined'!,logl2Cl(11._dp),dlnl2Cldlogl(11d0)

	!!! prepare for fishere elements, compute deltaCl

	!!! deltaCl=sqrt(2/(2l+1)fsky) (Cl+Nl) we need delta cl squared

	do i = 1, lmax, 1
           l=cl(i,1)
           !deltaCl(i)=(2d0/((2*cl(i,1)+1)*fsky))*(cl(i,2)+((4d-6)*(3.1415926535d0/10800))**2*exp(l*(l+1)*((thetaFWHM)**2)/(8d0*log(2d0))))**2 !! add noise
           deltaCl(i)=(2d0/((2*cl(i,1)+1)))*(cl(i,2))**2 !! no noise
           deltaCl_pol(i)=(2d0/((2*cl(i,1)+1)))*(cl(i,3))**2 !! no noise POL

	end do
		

        !! COMPUTE OUTPUT ANS PRINT
        open(unit=30, file='variancenonoisenofsky_fisher.txt')

!! fisher elements  sum_l 1/deltaC (deriv)^2
	fisher=0d0
	do i = 1, lmax, 1
		write(*,*) fisher
		fisher=fisher+(1/deltaCl(i))*(dlnl2Cldlogl(DBLE(i)))**2*(cl(i,2))**2
		write(30,*) cl(i,1),1/fisher

	end do
        close(30)
!! now get variance of the parameter from  fisher

        open(unit=60, file='variancenonoisenofsky_fisher_pol.txt')
        fisher_pol=0d0
	do i = 1, lmax, 1
	!	write(*,*) fisher
		fisher_pol=fisher_pol+(1/deltaCl_pol(i))*(dlnl2Cldlogl_pol(DBLE(i)))**2*(cl(i,3))**2
		write(60,*) cl(i,1),1/fisher_pol

	end do
        close(60)

	open(unit=20, file='derivatives.txt')

	do i = 1, lmax, 1

		write(20,*) cl(i,1),dlnl2Cldlogl(DBLE(i)),logl2Cl(dble(i))
		
	end do
	close(unit=20)

	open(unit=20, file='derivatives_pol.txt')

	do i = 1, lmax, 1

		write(20,*) cl(i,1),dlnl2Cldlogl_pol(DBLE(i)),logl2Cl_pol(dble(i))
		
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



DOUBLE PRECISION FUNCTION dlnl2Cldlogl_pol(ak)
  IMPLICIT none
  DOUBLE PRECISION :: a,b,h,x,y,ak
  x  = log(ak)
  CALL hunt(xa,lmax,x,jlo)
  h=xa(jlo+1)-xa(jlo)
  a=(xa(jlo+1)-x)/h
  b=(x-xa(jlo))/h
  y=(y_pol(jlo+1)-y_pol(jlo))/h+(-(3.*a**2-1.)*y_pol2(jlo)+(3.*b**2-1.)*y_pol2(jlo+1))*h/6.
  dlnl2Cldlogl_pol = y
  return
END FUNCTION dlnl2Cldlogl_pol
	 
	end program fisher_estimate



