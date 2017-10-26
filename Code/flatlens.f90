PROGRAM flatlens

  use random
  USE healpix_types
  USE fitstools, only: fits2cl
  use alm_tools, only: create_alm
  use rngmod, only: rand_init, planck_rng 
  use AMLutils, only: spline_double
  USE sub_flatlens

  implicit none


  include "fftw3.f"



!!!!========================================================!!!!!!!
!! ALLOCATE
!!!!========================================================!!!!!!!
  type(planck_rng) :: rng_handle
  
  character(len=128)  :: clfile
  integer(I4b):: lmax,ncl,ios
  integer(I8b)::i
  real(dp), dimension(:,:), allocatable :: ps1d,check,cl
  real(dp), dimension(:), allocatable:: CTT2,ell
  complex, dimension(:,:,:), allocatable :: almloaded
  character(len=80), dimension(1:70) :: header,units
  logical :: RandInited = .false.
  real(dp) ::cl3,tempCl,t1,t2


  INTEGER OMP_GET_MAX_THREADS
  INTEGER OMP_GET_NUM_THREADS
  INTEGER OMP_GET_THREAD_NUM


!!!!========================================================!!!!!!!
!! FFTW3 PARAMETERS
!!!!========================================================!!!!!!!
  integer(I8b)::plan_backward,plan_forward
  complex(kind=DPC), dimension(:,:), allocatable :: in,out

!!!!========================================================!!!!!!!


  !! grid parameters
!!! NOTE WE ARE STARTING FROM (0,0) indeces yu can always do that in FORTRAN
!!! initializing the array as below

!!!!========================================================!!!!!!!

  complex(kind=DPC), dimension(:,:,:),allocatable::GRIDHARMTPHI
  complex(kind=DPC), dimension(:,:,:),allocatable::GRIDTPHI
  integer(I4b),parameter::imax=2**10,jmax=2**10
  real(dp)::delx,dely
  real(dp)::tamp,phiamp, x,y,lenght,l
  integer(i4b)::icounter,jcounter,icenter,jcenter
  allocate(GRIDHARMTPHI(0:(imax-1),0:(jmax-1),1:2))
  allocate(GRIDTPHI(0:(imax-1),0:(jmax-1),1:2))
  allocate(check(0:(imax-1),0:(jmax-1)))
  check(:,:)=0
  !write(*,*) size(GRIDHARMTPHI)

!!!!========================================================!!!!!!!

!! set sizes load CL and allocate all the stuff we need
!!!!========================================================!!!!!!!


  !character(len=80):: units
  ncl=2
  lmax=3000

  !! find a clever way to automatically set the uvcell dimension
  !! to get a good sample between l=0 and l<lmax
  delx=max(1 , floor((SQRT2*lmax/imax)-1))
  dely=delx

  write(*,*) 'delx ',delx
  !!!!========================================================!!!!!!!

  allocate(cl(0:lmax,1:ncl))
  allocate(ps1d(0:lmax,1:1))
  allocate(almloaded(1:1,0:lmax,0:128))
!!!!========================================================!!!!!!!

  clfile='test_scal2Cls.dat'
  cl=0
  call read_camb(clfile, cl, lmax)

  !call fits2cl(clfile,cl,lmax,ncl,header,units)
  !! slpine Cls to fill non integer k (same as l ) on the 2d grid
  
  allocate(CTT2(0:lmax))
  allocate(ell(0:lmax))
  
!!!!========================================================!!!!!!!

!! this is only to get a double version of l in Cl maybee not needed

  DO i=0,lmax
     ell(i)=i
  END DO
 ! Write(*,*) cl(120,1)
  call spline_double(ell,cl(0:lmax,1),lmax+1,cTT2(0:lmax))
!!!!========================================================!!!!!!!

  !!!!========================================================!!!!!!!
  ! check spline
  !!!!========================================================!!!!!!!

  
! do l = 1, 1000, 1

! write(*,*) cl(l,1),tempCl
! !write(unit=iounit, fmt="(format string)", iostat=ios, advance='NO') cl(l,1),tempCl
! !if ( ios /= 0 ) stop "Write error in file unit iounit"


! end do
!!!!========================================================!!!!!!!

!!! now create lx ly coefficient in 2d Flat sky approximation
  !! Prepere for gaussian generation
  

  call InitRandom
  RandInited = .true.
  !Write(*,*) gaussian1()
!! then you can use Gaussian(1) to extract

!!!!========================================================!!!!!!!
!! grid
!!!!========================================================!!!!!!!


!!look notes---> thanks to simmetry we can loop on only the half of the grid.

!GOTO 100

icenter=0
jcenter=0
!write(*,*) jmax/2+1,imax/2+1

!!!!========================================================!!!!!!!
!!!! FILL THE GRID WITH CL OF TT IN THIS CASE
!!!! TODO GENERALIZE FOR PHI
!!!!========================================================!!!!!!!

  GRIDHARMTPHI(0,0,:)=0

  DO icounter=1,imax/2-1
    jcounter=0
    lenght=sqrt(((icounter-imax/2)*delx)**2+((jcounter-jmax/2)*dely)**2)
    call splint(dble(ell),cl(0:lmax,1),cTT2(0:lmax),lmax+1,lenght,tempCl)
    tamp = sqrt(tempCl/2d0)
    phiamp=sqrt(cl(int(lenght),2)/2)
    !write(*,*) icounter,jcounter,imax-icounter


    GRIDHARMTPHI(icounter,jcounter,1)= cmplx(Gaussian1(),Gaussian1())*tamp
    GRIDHARMTPHI(imax-icounter,jcounter,1)=  CONJG(GRIDHARMTPHI(icounter,jcounter,1))
    GRIDHARMTPHI(icounter,jcounter,2)= cmplx(Gaussian1(),Gaussian1())*phiamp
    GRIDHARMTPHI(imax-icounter,jcounter,2)=  CONJG(GRIDHARMTPHI(icounter,jcounter,2))
   ! check(icounter,jcounter)=1
  !  check(imax-icounter,jcounter)=1
 enddo


 do jcounter=1,jmax/2-1
  icounter=0
  lenght=sqrt(((icounter-imax/2)*delx)**2+((jcounter-jmax/2)*dely)**2)
  call splint(dble(ell),cl(0:lmax,1),cTT2(0:lmax),lmax+1,lenght,tempCl)
  tamp = sqrt(tempCl/2d0)
  phiamp=sqrt(cl(int(lenght),2)/2)

  GRIDHARMTPHI(icounter,jcounter,1)= cmplx(Gaussian1(),Gaussian1())*tamp
    GRIDHARMTPHI(icounter,jmax-jcounter,1)= CONJG(GRIDHARMTPHI(icounter,jcounter,1))
    GRIDHARMTPHI(icounter,jcounter,2)= cmplx(Gaussian1(),Gaussian1())*phiamp
    GRIDHARMTPHI(imax-icounter,jcounter,2)=  CONJG(GRIDHARMTPHI(icounter,jcounter,2))
     
   !  check(icounter,jcounter)=1
    
  !  check(icounter,jmax-jcounter)=1
 end do


!OMP PARALLEL DO PRIVATE(icounter,jcounter)
    write(*,*) 'here',imax/2-1

do icounter=1,imax/2-1
  do jcounter=1,jmax-1
    

      lenght=sqrt(((icounter-imax/2)*delx)**2+((jcounter-jmax/2)*dely)**2)
      call splint(dble(ell),cl(0:lmax,1),cTT2(0:lmax),lmax+1,lenght,tempCl)
      tamp = sqrt(tempCl/2d0)
      phiamp=sqrt(cl(int(lenght),2)/2)

    !write(*,*) icounter,jcounter,imax-icounter,jmax-jcounter

    GRIDHARMTPHI(icounter,jcounter,1)= cmplx(Gaussian1(),Gaussian1())*tamp
    GRIDHARMTPHI(imax-icounter,jmax-jcounter,1)= Conjg(GRIDHARMTPHI(icounter,jcounter,1))
    GRIDHARMTPHI(icounter,jcounter,2)= cmplx(Gaussian1(),Gaussian1())*phiamp
    GRIDHARMTPHI(imax-icounter,jcounter,2)=  CONJG(GRIDHARMTPHI(icounter,jcounter,2))
    !         check(icounter,jcounter)=1
   ! check(imax-icounter,jmax-jcounter)=1

  end do
!OMP END PARALLEL DO NOWAIT



  
end do

do jcounter=1,jmax/2-1
  icounter=imax/2

  lenght=sqrt(((icounter-icenter)*delx)**2+((jcounter-jcenter)*dely)**2)
  call splint(dble(ell),cl(0:lmax,1),cTT2(0:lmax),lmax+1,lenght,tempCl)
  tamp = sqrt(tempCl/2d0)
  phiamp=sqrt(cl(int(lenght),2)/2)


  GRIDHARMTPHI(icounter,jcounter,1)= cmplx(Gaussian1(),Gaussian1())*tamp
  GRIDHARMTPHI(imax-icounter,jmax-jcounter,1)= conjg(GRIDHARMTPHI(icounter,jcounter,1))
  GRIDHARMTPHI(icounter,jcounter,2)= cmplx(Gaussian1(),Gaussian1())*phiamp
  GRIDHARMTPHI(imax-icounter,jcounter,2)=  CONJG(GRIDHARMTPHI(icounter,jcounter,2))
   !check(icounter,jcounter)=1
   ! check(imax-icounter,jmax-jcounter)=1
   end do
 
 icounter=imax/2
 jcounter=0
  lenght=sqrt(((icounter-icenter)*delx)**2+((jcounter-jcenter)*dely)**2)
call splint(dble(ell),cl(0:lmax,1),cTT2(0:lmax),lmax+1,lenght,tempCl)
  tamp = sqrt(tempCl)
  phiamp=sqrt(cl(int(lenght),2))

  GRIDHARMTPHI(icounter,jcounter,1)= cmplx(Gaussian1(),0)*tamp
  GRIDHARMTPHI(jcounter,icounter,1)= cmplx(Gaussian1(),0)*tamp
  GRIDHARMTPHI(icounter,jcounter,1)= cmplx(Gaussian1(),0)*phiamp
  GRIDHARMTPHI(jcounter,icounter,1)= cmplx(Gaussian1(),0)*phiamp

    ! check(icounter,jcounter)=1
   ! check(jcounter,icounter)=1


  jcounter=jmax/2
  call splint(dble(ell),cl(0:lmax,1),cTT2(0:lmax),lmax+1,lenght,tempCl)
  tamp = sqrt(tempCl)
  phiamp=sqrt(cl(int(lenght),2))

  GRIDHARMTPHI(icounter,jcounter,1)= cmplx(Gaussian1(),0)*tamp  
  GRIDHARMTPHI(icounter,jcounter,1)= cmplx(Gaussian1(),0)*phiamp  

    !check(icounter,jcounter)=1



!Write(*,*) GRIDHARMTPHI(0,:,1),'\'
!Write(*,*) GRIDHARMTPHI(1,:,1)
!Write(*,*) GRIDHARMTPHI(2,:,1)
!Write(*,*) GRIDHARMTPHI(3,:,1)
!Write(*,*) GRIDHARMTPHI(4,:,1)
!Write(*,*) GRIDHARMTPHI(5,:,1)
!Write(*,*) GRIDHARMTPHI(6,:,1)
!Write(*,*) GRIDHARMTPHI(7,:,1)


!!!!========================================================!!!!!!!
!!!! CHECK IF I FILL THE GRID CORRECTLY
!!!!========================================================!!!!!!!

write(*,*) GRIDHARMTPHI(0,0,1)
call ps_azim(imax,jmax,GRIDHARMTPHI(:,:,2),icenter,jcenter,lmax,delx,ps1d)

open(unit=10, file='sanitycheckCl.dat', iostat=ios, status="REPLACE", action="write")
if ( ios /= 0 ) stop "Error opening file name"

do l = 1, lmax, 1
!write(*,"(E14.7)") ps1d(l,1)
write(10,'(I5.2,E14.7)') int(l), ps1d(l,1)

end do

close(unit=10, iostat=ios)
if ( ios /= 0 ) stop "Error closing file unit 10"


!!!!========================================!!!
  ! INITIALIZE THE FOURIER PLAN 
!!!!========================================!!!


200 call dfftw_plan_dft_2d_ ( plan_backward, imax, jmax, GRIDHARMTPHI(:,:,2), GRIDTPHI(:,:,2), FFTW_BACKWARD, &
    FFTW_ESTIMATE )


  call dfftw_execute_ (plan_backward)

!Write(*,*) GRIDHARMTPHI(0,:,1)
!Write(*,*) GRIDHARMTPHI(1,:,1)
!Write(*,*) GRIDHARMTPHI(2,:,1)
!Write(*,*) GRIDHARMTPHI(3,:,1)
!write(*,*) '\n'

!Write(*,*) GRIDTPHI(0,:,1)
!Write(*,*) GRIDTPHI(1,:,1)
!Write(*,*) GRIDTPHI(2,:,1)
!Write(*,*) GRIDTPHI(3,:,1)
  

  call dfftw_plan_dft_2d_ ( plan_forward,imax, jmax, GRIDTPHI(:,:,2), GRIDHARMTPHI(:,:,2), FFTW_FORWARD, &
    FFTW_ESTIMATE)

  call dfftw_execute_ ( plan_forward )

  

  call ps_azim(imax,jmax,GRIDHARMTPHI(:,:,2),icenter,jcenter,lmax,delx,ps1d)

open(unit=10, file='sanitycheck2Cl.dat', iostat=ios, status="REPLACE", action="write")
if ( ios /= 0 ) stop "Error opening file name"

do l = 1, lmax, 1
!write(*,"(E14.7)") ps1d(l,1)
write(10,'(I5.2,E14.7)') int(l), ps1d(l,1)/(real(imax*jmax,DP)**2)

end do

close(unit=10, iostat=ios)
if ( ios /= 0 ) stop "Error closing file unit 10"

  write(*,*)  shape(GRIDHARMTPHI(:,:,1)),shape(GRIDTPHI(:,:,1))!,GRIDTPHI(:,:,1)

  call dfftw_destroy_plan_ ( plan_forward )
  call dfftw_destroy_plan_ ( plan_backward )


stop




 !-----------------------------------------------------------------------
  !                       end of routine
  !-----------------------------------------------------------------------

  WRITE(*,9000) " "
  WRITE(*,9000) " flatlens > normal completion"

9000 format(a)



















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
