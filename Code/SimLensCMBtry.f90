program SimLensCMBtry
 use HealpixObj
 use HealpixVis
 use Random
 use spinalm_tools
 use IniFile
 use AMLUtils
 implicit none
 Type(HealpixInfo)  :: H
 Type(HealpixMap)   :: M, GradPhi
 Type(HealpixPower) :: P
 Type(HealpixAlm)   :: A
 
 integer            :: nside, lmax
 integer(I_NPIX)    :: npix
 character(LEN=256)  :: w8name = '../Healpix_2.00/data/'
 character(LEN=256)  :: file_stem, cls_file, out_file_root, cls_lensed_file
 integer, parameter :: lens_interp =1, lens_exact = 2
 integer :: lens_method = lens_interp
 integer :: mpi_division_method = division_equalrows
 integer ::  interp_method,  rand_seed
 logical :: err, want_pol
 real :: interp_factor
#ifdef MPIPIX
 integer i
 
 call mpi_init(i)
#endif

 Ini_Fail_On_Not_Found = .true.
 call Ini_Open(GetParam(1), 3,err)
 if (err) then
#ifdef MPIPIX
    call mpi_finalize(i)
#endif
   stop 'No ini'
 end if
 nside  = Ini_Read_Int('nside')
 npix = nside2npix(nside)

 lmax   = Ini_Read_Int('lmax')  
 cls_file = Ini_Read_String('cls_file')
 x out_file_root = Ini_Read_String('out_file_root')

 want_pol = Ini_Read_Logical('want_pol')
 rand_seed = Ini_Read_Int('rand_seed')

 interp_method = Ini_read_int('interp_method')
 
 Ini_Fail_On_Not_Found = .false.
 
 
 w8name = Ini_Read_String('w8dir')
 interp_factor=0
 if (lens_method == lens_interp) interp_factor = Ini_Read_Real('interp_factor',3.)
#ifdef MPIPIX
 mpi_division_method = Ini_Read_Int('mpi_division_method',division_balanced);
#endif 

 call Ini_Close

 file_stem =  trim(out_file_root)//'_lmax'//trim(IntToStr(lmax))//'_nside'//trim(IntTOStr(nside))// &
              '_interp'//trim(RealToStr(interp_factor,3))//'_method'//trim(IntToStr(interp_method))//'_'

 if (want_pol) file_stem=trim(file_stem)//'pol_'
 file_stem = trim(file_stem)//trim(IntToStr(lens_method)) 
 
 cls_lensed_file  = trim(file_stem)//'.dat'
  
 call SetIdlePriority();

!Healpix 2.01 has up to 8192 
! if (nside > 1024) then
 if (w8name=='') then
  write (*,*) 'Warning: using unit weights as no w8dir specified'
  call HealpixInit(H,nside, lmax,.true., w8dir='', method= mpi_division_method) 
 else
  call HealpixInit(H,nside, lmax,.true., w8dir=w8name,method=mpi_division_method) 
 end if 

 if (H%MpiID ==0) then !if we are main thread
  !All but main thread stay in HealpixInit
 
  call HealpixPower_nullify(P)
  call HealpixAlm_nullify(A)
  call HealpixMap_nullify(GradPhi)
  call HealpixMap_nullify(M)

  call HealpixPower_ReadFromTextFile(P,cls_file,lmax,pol=.true.,dolens = .true.)
  !Reads in unlensed C_l text files as produced by CAMB (or CMBFAST if you aren't doing lensing)
  call HealpixPower_ReadFromTextFile(P,cls_file,lmax,pol=.true.,dolens = .true.)
  !Reads in unlensed C_l text files as produced by CAMB (or CMBFAST if you aren't doing lensing)

  call HealpixAlm_Sim(A, P, rand_seed,HasPhi=.true., dopol = want_pol)
  call HealpixAlm2Power(A,P)
  call HealpixPower_Write(P,trim(file_stem)//'_unlensed_simulated.dat')

END PROGRAM SimLensCMBtry
