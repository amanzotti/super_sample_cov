Module sub_flatlens
  USE healpix_types
!USE flatlens_types

Implicit None

type grid

integer(i4b) :: nx,ny
logical isT,isPhi
	
end type grid



Contains

subroutine getcenter(mygrid,centercoord)
	type (grid), intent(in) :: mygrid
		integer(i8b), intent(out),dimension(0:1) :: centercoord
		if ( (mod(mygrid%nx,2).ne.0) .or. (mod(mygrid%ny,2) .ne. 0) ) then
			stop "only implemented for even grid dimensions"	
		end if
		centercoord(0)=mygrid%nx/2
		centercoord(1)=mygrid%ny/2
	
end subroutine getcenter

subroutine real2matcoord(mygrid,realcoor,mat)
	type (grid), intent(in) :: mygrid

	integer(i8b), intent(in),dimension(0:1) :: realcoor 

	integer(i8b), intent(out),dimension(0:1) :: mat
	integer(i8b), dimension(0:1) :: center
	call getcenter(mygrid,center)
	mat(0)=2*center(0)+realcoor(0)
	mat(1)=2*center(1)-realcoor(1)


end subroutine real2matcoord

subroutine mat2realcoord(mygrid,realcoor,mat)
	type (grid), intent(in) :: mygrid

	integer(i8b), intent(in),dimension(0:1) :: realcoor 

	integer(i8b), intent(out),dimension(0:1) :: mat
	integer(i8b), dimension(0:1) :: center
	call getcenter(mygrid,center)
	mat(0)=2*center(0)-realcoor(0)
	mat(1)=2*center(1)-realcoor(1)


end subroutine mat2realcoord



!subroutine grid2cl(mygrid,Cl)
!	complex(DPC), intent(in), dimension(:,:) :: mygrid
!	real(DP),intent(out),dimension(:,:)::Cl
!end subroutine grid2cl


!------------------------------------------------------------------------         

! subroutine to compute the ps1d(l) by azimuthal averaging the power in each fourier mode

        SUBROUTINE ps_azim(nx,ny,map,icentre,jcentre,lmax,uvcell,ps1d)
        	USE healpix_types

          Implicit None

        Integer(i4b), Intent(In) :: nx, ny, icentre,jcentre,lmax

        Complex(DPC), Intent(In) ::  map(0:nx-1,0:ny-1)
        real(DP), Intent(Out)   :: ps1d(0:lmax)

        Real(DP)  :: uvdist,sumps,uvcell
        Integer :: i,j,k,l, count(0:lmax)      
       write(*,*) lmax
        write(*,*) 'center',icentre,jcentre
! initialize variables

        count = 0
        ps1d = 0d0

! Work out which value of l each pixel corresponds to and add up
! the power for each l

!write(*,*) size(map)

     Do i = 0,nx-1
        Do j = 0,ny-1
           uvdist = uvcell * Sqrt(real(i-nx/2)**2+real(j-ny/2)**2)
          
           l = uvdist


           If((l .le. lmax)) Then
           	  !write(*,*)'l',l,i,j

              ps1d(l) = ps1d(l)+Abs(map(i,j))**2
              count(l) = count(l) + 1
             ! write(*,*)'count(l)', count(l)
          End if
         End Do
      End Do

! the azimuthal average is scaled by uvcell**2 

      !ps1d = ps1d * (2.0*pi*uvcell)**2

      !sumps=0.
      !Do k=lmin,lmax
      !   If(ps1d(k).ne.0) Then
       !     sumps = sumps+ps1d(k)
       do l = 1, lmax, 1
           	
           	if ( count(l)/=0 ) then !! prevent by 0 division
           		
            ps1d(l) = ps1d(l)/count(l)

           ! write(*,*) l,'  ',ps1d(l)


        end if

       
        end do    


       !  Endif
     ! End Do
      
      ! Normalise sumps by map size
      !sumps = sumps * uvcell**2
      
     ! Write(*,*) 'Sqrt(total power) = ', sqrt(sumps)

        END SUBROUTINE ps_azim



 Subroutine read_camb(infile, cl_array, lmax)
   
    Use healpix_types

    Implicit none

    Character (len = *), intent(IN) :: infile
    Integer(I4B), intent(INOUT) :: lmax
    Real(DP), Intent(OUT), Dimension(0:lmax,1:2) :: cl_array

    Real(DP) :: clT,clE,clET,clphi
    Integer :: lval = 0
    Real(DP) :: cmbT = 2.726d0
    cl_array(0:1,:)=0
    ! Read in the cl values and convert to units of uK^2
    Open (Unit = 23, file = infile, status = 'old', action = 'read', err = 9980)
    Do while (lval .LT. lmax)
       Read (23, *, end = 800) lval, clT,clE,clET,clphi
       cl_array(lval,1) = clT!/(lval*(lval+1))*cmbT*cmbT*1.d12*TWOPI
       cl_array(lval,2) = clphi!/(lval**4)*cmbT*cmbT*1.d12*TWOPI

    End Do
800 lmax = lval

    Close (23)
    Return

9980 Stop 'Error opening cl file to read'

  End Subroutine read_camb



end module sub_flatlens















