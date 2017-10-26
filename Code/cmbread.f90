PROGRAM flatlens
  use random
  USE healpix_types
  USE fitstools, only: fits2cl
  use alm_tools, only: create_alm
  use rngmod, only: rand_init, planck_rng 
  use AMLutils, only: spline_double
  implicit none
  type(planck_rng) :: rng_handle
  
  character(len=128)  :: clfile
  integer(I4b):: lmax,ncl
  integer(I8b)::i
  real(dp), dimension(:,:), allocatable :: cl
  real(dp), dimension(:), allocatable:: CTT2,ell
  complex, dimension(:,:,:), allocatable :: almloaded
  character(len=80), dimension(1:70) :: header,units
  logical :: RandInited = .false.
  real(dp) ::cl3,tempCl

  !! grid parameters

  complex, dimension(:,:,:),allocatable::PHYXY
  integer(I8b),parameter::xmax=2**3,ymax=2**3
  real(dp),parameter::delx=0.5d0,dely=0.5d0
  real(dp)::tamp,x,y,lenght
  integer(i8b)::xcounter,ycounter,xcenter,ycenter
  allocate(PHYXY(0:(xmax-1),0:(ymax-1),1:1))
  write(*,*) size(PHYXY)

!!!!


  !character(len=80):: units
  ncl=1
  lmax=1000
  allocate(cl(0:lmax,1:ncl))
  allocate(almloaded(1:1,0:lmax,0:128))
  clfile='manzotti_scalCls.fits'
  call fits2cl(clfile,cl,lmax,ncl,header,units)
  !! slpine Cls to fill non integer k (same as l ) on the 2d grid
  allocate(CTT2(0:lmax))
  allocate(ell(0:lmax))
  
  DO i=0,lmax
     ell(i)=i
  END DO
  
  call spline_double(ell,cl(0:lmax,1),lmax+1,cTT2(0:lmax))

  stop
  
do l = 1, 1000, 1
  call splint(ell,cl(0:lmax,1),cTT2(0:lmax),cTT2(0:lmax),tempCl)
!write(unit=iounit, fmt="(format string)", iostat=ios, advance='NO') cl(l,1),tempCl
!if ( ios /= 0 ) stop "Write error in file unit iounit"


end do

  
!!! now create lx ly coefficient in 2d Flat sky approximation
  !! Prepere for gaussian generation
  
  call InitRandom
  RandInited = .true.

  !Write(*,*) gaussian1()
!! then you can use Gaussian(1) to extract


!! grid

x=0d0

!!look notes---> thanks to simmetry we can loop on only the half of the grid.

xcenter=xmax/2+1
ycenter=ymax/2+1
!write(*,*) ymax/2+1,xmax/2+1

  DO xcounter=0,xmax-1
    
     y=0d0

     
     DO ycounter=0,ymax/2+1
        
        lenght=sqrt(((xcounter-xcenter)*delx)**2+((ycounter-ycenter)*dely)**2)
       ! write(*,*) xcounter,ycounter
       ! write(*,*)dble((xcounter-xcenter)*delx),dble((xcounter-xcenter))*dely,lenght
       ! write(*,*)


        call splint(ell,cl(0:lmax,1),cTT2(0:lmax),lmax+1,lenght,cl3)
       
        tamp = sqrt(cl3/2)

        ! note in the center line you are doing this 2 times, improve this for performance      
  
        PHYXY(xcounter,ycounter,1)= cmplx(Gaussian1(),Gaussian1())*tamp
        PHYXY(-xcounter,-ycounter,1)= cmplx(Gaussian1(),-Gaussian1())*tamp

        write(*,*) PHYXY(xcounter,ycounter,1), PHYXY(-xcounter,-ycounter,1),xcounter,ycounter
        
       
     END DO



      !! At this point you shoulf FFT
  

  END DO










CONTAINS 





   SUBROUTINE splint(xa, ya, y2a, n, x, y)
!   USE nrtype
!
! Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a function
! (with the xa(i) in order), and given the array y2a(1:n), which is the output
! from the subroutine spline, and given a value of x, this routine returns a
! cubic spline interpolated value y.
! (adopted from Numerical Recipes in FORTRAN 77)
!
   INTEGER, PARAMETER :: DP = KIND(1.0D0)
   INTEGER:: n
   REAL(DP):: x, y, xa(n), y2a(n), ya(n)
   INTEGER:: k, khi, klo
   REAL(DP):: a, b, h

     klo=1
     khi=n
1   if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if (xa(k).gt.x) then
           khi=k
        else
           klo=k
        endif
        goto 1
     endif

     h=xa(khi)-xa(klo)
     if (h.eq.0.) pause 'bad xa input in splint'

     a=(xa(khi)-x)/h
     b=(x-xa(klo))/h
     y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.

     return
   
   END SUBROUTINE splint



END PROGRAM flatlens
