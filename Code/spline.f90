


SUBROUTINE spline(x,y,n,yp1,ypn,y2)
  USE healpix_types 
  implicit none
      INTEGER, intent(in) :: n
      real(dp), intent(in) :: x(n), y(n), yp1, ypn
      real(dp), intent(out) :: y2(n)
      INTEGER i,k
      real(dp) p,qn,sig,un
      real(dp), dimension(:), allocatable :: u

      write(*,*) 'you are using this one'
       
      Allocate(u(1:n))
      if (yp1.gt..1d30) then
        y2(1)=0._dp
        u(1)=0._dp
      else
        y2(1)=-0.5d0
        u(1)=(3._dp/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      
      do i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2._dp 
   
        y2(i)=(sig-1._dp)/p
      
         u(i)=(6._dp*((y(i+1)-y(i))/(x(i+ &
         1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig* &
         u(i-1))/p
      end do
      if (ypn.gt..1d30) then
        qn=0._dp
        un=0._dp
      else
        qn=0.5d0
        un=(3._dp/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1._dp)
      do k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
      end do
      
      Deallocate(u)
  
!  (C) Copr. 1986-92 Numerical Recipes Software =$j*m,).
      END SUBROUTINE spline
