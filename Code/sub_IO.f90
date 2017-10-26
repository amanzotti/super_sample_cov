Module sub_output

Implicit none

Contains

!**************************************************
! Subroutine to read in power spectrum, either from a fits file
! using Healpix subroutine or from a file which is of the form
! l-value, C_l-value 
! also computes total power 'sum'
!**************************************************
Subroutine read_powerspec(infile,lmin, lmax, cls, sum)

Use healpix_types
USE fitstools, ONLY : read_asctab
use utilities, only : die_alloc

Implicit none

  Character(len = *), Intent(In) :: infile
  Integer, Intent(In)       :: lmin
  Integer, Intent(InOut)    :: lmax
  Real(SP), dimension(0:lmax), Intent(Out) :: cls
  Real , Intent(Out) :: sum

  Logical fitscl
  Integer(I4B) :: nlheader
  Integer(I4B) :: l_max, l
  CHARACTER(LEN=80), DIMENSION(1:180)          :: header_file
  Real(SP) :: clread(0:lmax,1)

!-----------------------------------------------------------------------

  cls = 0.0
  nlheader = SIZE(header_file)
  l_max = lmax

  ! Find out whether input cl file is a fits file
  fitscl = (index(infile, '.fits') /= 0)
  If (fitscl) Then
     clread = 0.0
     call read_asctab(infile, clread, lmax, 1, header_file, nlheader)
     cls(0:lmax) =  clread(0:lmax,1)
  Else
         call read_camb(infile, cls, lmax)
     !Check to see if enough values were present in file
     If (lmax .LT. l_max) then
        Write (*,*) "Maximum l-value in cl file is ",lmax
        Write (*,*) "Power spectrum above this value set to zero"
!        l_max = lmax
     end if
  End If

sum = 0.0
  Do l = 2, lmax
     If (l .lt. lmin .or. l .gt. lmax) then
        sum = sum + 0.0
     else
        sum = sum +(real(2.0*l+1.0)*cls(l))
     endif     
  End Do
sum = sum/FOURPI

End Subroutine read_powerspec

! Subroutine called by read_powerspec to read
! in .dat file of the form output by CAMB 

 Subroutine read_camb(infile, cl_array, lmax)
   
    Use healpix_types

    Implicit none

    Character (len = *), intent(IN) :: infile
    Integer(I4B), intent(INOUT) :: lmax
    Real(SP), Intent(OUT), Dimension(0:lmax) :: cl_array

    Real(DP) :: clval
    Integer :: lval = 0
    Real(DP) :: cmbT = 2.726d0

    ! Read in the cl values and convert to units of uK^2
    Open (Unit = 23, file = infile, status = 'old', action = 'read', err = 9980)
    Do while (lval .LT. lmax)
       Read (23, *, end = 800) lval, clval
       cl_array(lval) = TWOPI*clval/(lval*(lval+1))*cmbT*cmbT*1.d12
    End Do
800 lmax = lval

    Close (23)
    Return

9980 Stop 'Error opening cl file to read'

  End Subroutine read_camb

!*************************************************
! Subroutine to write out values in Fourier space
!*************************************************

Subroutine write_uv(map, nx, ny, cellx, celly)

Complex, Intent(In) :: map(nx,ny)
Integer, Intent(In) :: nx, ny
Real, Intent(In)    :: cellx, celly
!Real, parameter     :: cmbtemp = 2.726e6

Integer i,j, icentre, jcentre, count, unit
Character(len = 50) :: fname
Real :: error = 0.0

unit = 15
Write (*,*) 'Input filename for output in uv plane'
Read (*,*) fname

Open(unit=unit, action = 'write', file = fname)

icentre=nx/2+1
jcentre=ny/2+1 

!count = 1
Do i = 1, nx
   Do j = 1, ny
      Write (unit, '(I5, 2F11.4, 2E14.5, F5.1)') count, (i-icentre)*cellx, (j-jcentre)*celly,&
& real(map(i,j)), imag(map(i,j)), error
      count = count + 1
      If (count .gt. 99999) count = 1
   End Do
EndDo

Close(unit)

End Subroutine write_uv

End Module sub_output
