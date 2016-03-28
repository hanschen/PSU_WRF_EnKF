
      module module_wrfchem_lib

        use netcdf
      implicit none

!     include files
      include 'netcdf.inc'

!     public  procedures
      public  :: init_wrfchem_lib, exit_wrfchem_lib
      public  :: wrfchem_readscalar, wrfchem_read2d,  wrfchem_read3d
      public  :: wrfchem_write4d

!     private procedures
      private :: handle_error, define_variable

!     public variables
      integer, public  :: nx, ny, nz, nw, ntime      ! data dimension length
      integer, private :: nscalar, nstring           ! data dimension length
      integer, public  :: year_start, month_start, day_start
      integer, public  :: hour_start, minute_start, second_start
      integer, public  :: year_end, month_end, day_end
      integer, public  :: hour_end, minute_end, second_end
      integer, public  :: dt                        ! temporal resolution

!     variables for reading netCDF data
      integer          :: ncid, ncidi                          ! netCDF file IDs
      integer, private :: nstt(4), ncnt(4), ncntx(4), ncnty(4) ! start index and counter array

      contains

!---------------------------------------------------------------------------------------------------
!     initialize netcdf file
!---------------------------------------------------------------------------------------------------

      subroutine init_wrfchem_lib( wrfchem_bdy_fn, wrfchem_input_fn, specname, nspec, ntime )

      implicit none

!     input arguments
      integer, intent(in) :: nspec
      integer, intent(out) :: ntime
      character(len=*), intent(in) ::  wrfchem_bdy_fn, wrfchem_input_fn
!     declaration of specname changed to match what is in main_bc_wrfchem.f90  2014/09/29      
      !character(len=8), intent(in), DIMENSION(nspec) :: specname
      character(len=100), intent(in), DIMENSION(nspec) :: specname

!     local arguments
      integer :: status
      integer :: xdimid,  ydimid, zdimid, wdimid, tdimid, sdimid
      integer :: ndims, dimids(4)
      integer :: vid
      integer :: ip, n, ns, i
      integer :: year_first, month_first, day_first
      integer :: hour_first, minute_first, second_first
      character(len=100) :: dtstrings, dtstringe, dtstringf
      character(len=100)  :: spec1
      character(len=3)   :: order

      ! open of wrfinput changed to nowrite 2014/09/29
      write(*,*) 'init_wrfchem_lib: specname, nspec ',specname, nspec  !debugging display 2014/09/29

      write(*,*) 'open the bdy NETCDF file : ',wrfchem_bdy_fn
      status = nf_open( wrfchem_bdy_fn, nf_write, ncid )
      if( status /= nf_noerr ) call handle_error( status )

      write(*,*) 'open the input NETCDF file : ', wrfchem_input_fn
      status = nf_open( wrfchem_input_fn, nf_nowrite, ncidi )
      if( status /= nf_noerr ) call handle_error( status )

      write(*,*) 'get spacial dimension lengths'
      status = nf_inq_dimid( ncid, 'west_east', xdimid )
      if( status /= nf_noerr ) call handle_error( status )

      status = nf_inq_dimlen( ncid, xdimid, nx )
      if( status /= nf_noerr ) call handle_error( status )

      status = nf_inq_dimid( ncid, 'south_north', ydimid )
      if( status /= nf_noerr ) call handle_error( status )

      status = nf_inq_dimlen( ncid, ydimid, ny )
      if( status /= nf_noerr ) call handle_error( status )

      status = nf_inq_dimid( ncid, 'bottom_top', zdimid )
      if( status /= nf_noerr ) call handle_error( status )

      status = nf_inq_dimlen( ncid, zdimid, nz )
      if( status /= nf_noerr ) call handle_error( status )

      status = nf_inq_dimid( ncid, 'bdy_width', wdimid )
      if( status /= nf_noerr ) call handle_error( status )

      status = nf_inq_dimlen( ncid, wdimid, nw )
      if( status /= nf_noerr ) call handle_error( status )

      write(*,*) 'get the time dimension length'
      status = nf_inq_dimid( ncid, 'Time', tdimid )
      if( status /= nf_noerr )  call handle_error( status )

      status = nf_inq_dimlen( ncid, tdimid, ntime )
      if( status /= nf_noerr ) call handle_error( status )

!     get the scalar dimension lenth
!      status = nf_inq_dimid( ncid, 'ext_scalar', sdimid )
!      if( status /= nf_noerr )  call handle_error( status )
!
!      status = nf_inq_dimlen( ncid, sdimid, nscalar )
!      if( status /= nf_noerr ) call handle_error( status )
!      print*,'here',nscalar

      write(*,*) 'get the string dimension length'
      status = nf_inq_dimid( ncid, 'DateStrLen', sdimid )
      if( status /= nf_noerr )  call handle_error( status )

      status = nf_inq_dimlen( ncid, sdimid, nstring )
      if( status /= nf_noerr ) call handle_error( status )

      write(*,*) 'enter redefine mode'
      status = nf_redef( ncid )
      if( status /= nf_noerr ) call handle_error( status )

      ! Note that these variables should already exist! 

!     define variables      
      ndims = 4
      
      ! west-east
      dimids(1:4) = (/ ydimid, zdimid, wdimid, tdimid /)
      ! west
      do i = 1, 2
         do ns = 1, nspec
            if(i==1) spec1 = trim(specname(ns))//'_BXS'
            if(i==2) spec1 = trim(specname(ns))//'_BTXS'
            status = nf_inq_varid( ncid, trim(spec1), vid )
            !write(*,*) 'ncid, spec1, vid, status: ', ncid, trim(spec1), vid, status !temp debugging 2014/09/29
            if( status /= nf_noerr ) then 
               print*, 'Define new species : ', trim(spec1)
               order = 'YSZ'
               call define_variable( trim(spec1), nf_real, ndims, dimids, order )
            endif
         enddo
      end do

      ! east
      do i = 1, 2 
         do ns = 1, nspec
         if(i==1) spec1 = trim(specname(ns))//'_BXE'
         if(i==2) spec1 = trim(specname(ns))//'_BTXE'
            status = nf_inq_varid( ncid, trim(spec1), vid )
            if( status /= nf_noerr ) then 
               print*, 'Define new species : ', trim(spec1)
               order = 'YEZ'
               call define_variable( trim(spec1), nf_real, ndims, dimids, order )
            endif
         enddo
      end do
      
      ! north-south
      dimids(1:4) = (/ xdimid, zdimid, wdimid, tdimid /)
      ! north
      do i = 1, 2 
         do ns = 1, nspec
         if(i==1) spec1 = trim(specname(ns))//'_BYS'
         if(i==2) spec1 = trim(specname(ns))//'_BTYS'
            status = nf_inq_varid( ncid, trim(spec1), vid )
            if( status /= nf_noerr ) then 
               print*, 'Define new species : ', trim(spec1)
               order = 'XEZ'
               call define_variable( trim(spec1), nf_real, ndims, dimids, order )
            endif
         enddo
      end do

      ! south
      do i = 1, 2 
         do ns = 1, nspec
         if(i==1) spec1 = trim(specname(ns))//'_BYE'
         if(i==2) spec1 = trim(specname(ns))//'_BTYE'
            status = nf_inq_varid( ncid, trim(spec1), vid )
            if( status /= nf_noerr ) then 
               print*, 'Define new species : ', trim(spec1)
               order = 'XSZ'
               call define_variable( trim(spec1), nf_real, ndims, dimids, order )
            endif
         enddo
      end do
 

!     end define mode
      status = nf_enddef( ncid )
      if( status /= nf_noerr ) call handle_error( status )

      write(*,*) 'read start and end date and second time stamp'
      status = nf_inq_varid( ncid, 'Times', vid )
      if( status /= nf_noerr )  call handle_error( status )

      nstt(1:2) = (/ 1, 1 /)
      ncnt(1:2) = (/ nstring, 1 /)
      status = nf90_get_var( ncid, vid, dtstrings, nstt(1:2), ncnt(1:2) )

      nstt(1:2) = (/ 1, ntime /)
      status = nf90_get_var( ncid, vid, dtstringe, nstt(1:2), ncnt(1:2) )

      nstt(1:2) = (/ 1, 2 /)
      status = nf90_get_var( ncid, vid, dtstringf, nstt(1:2), ncnt(1:2) )

      write(*,*) 'get the location of _ '
      ip = index( dtstrings, '_' )

      write(*,*) 'delete -, _ and : '
      do n = 1, nstring 
        if( (dtstrings(n:n)=='-') .or. (dtstrings(n:n)=='_') .or.  &
             dtstrings(n:n)==':' ) then
          dtstrings(n:n) = ' '
        endif
        if( (dtstringe(n:n)=='-') .or. (dtstringe(n:n)=='_') .or.  &
             dtstringe(n:n)==':' ) then
          dtstringe(n:n) = ' '
        endif
        if( (dtstringf(n:n)=='-') .or. (dtstringf(n:n)=='_') .or.  &
             dtstringf(n:n)==':' ) then
          dtstringf(n:n) = ' '
        endif
      enddo

      write(*,*) 'read year, month and day'
      read( dtstrings(:ip), * ) year_start, month_start, day_start
      read( dtstrings(ip:nstring), * ) hour_start, minute_start, second_start
      !read( dtstringe(:ip), * ) year_end, month_end, day_end
      !read( dtstringe(ip:nstring), * ) hour_end, minute_end, second_end
      !read( dtstringf(:ip), * ) year_first, month_first, day_first
      !read( dtstringf(ip:nstring), * ) hour_first, minute_first, second_first
      year_end=year_start
      if ((month_start.eq.1).or.(month_start.eq.3).or.(month_start.eq.5).or.(month_start.eq.7).or.(month_start.eq.8).or.(month_start.eq.10).or.(month_start.eq.12)) then
         if ((day_start.eq.31).and.(hour_start.eq.18)) then
            month_end=month_start+1
            day_end=1
            hour_end=0
         else
            month_end=month_start
            if (hour_start.ne.18) then
               day_end=day_start
               hour_end=hour_start+6.
            else
               day_end=day_start+1
               hour_end=0
            endif
         endif
      elseif ((month_start.eq.4).or.(month_start.eq.6).or.(month_start.eq.9).or.(month_start.eq.11)) then
         if ((day_start.eq.30).and.(hour_start.eq.18)) then
            month_end=month_start+1
            day_end=1
            hour_end=0
         else
            month_end=month_start
            if (hour_start.ne.18) then
               day_end=day_start
               hour_end=hour_start+6.
            else
               day_end=day_start+1
               hour_end=0
            endif
         endif
      else
         if ((day_start.eq.28).and.(hour_start.eq.18)) then
            month_end=month_start+1
            day_end=1
            hour_end=0
         else
            month_end=month_start
            if (hour_start.ne.18) then
               day_end=day_start
               hour_end=hour_start+6.
            else
               day_end=day_start+1
               hour_end=0
            endif
         endif
      endif

      minute_end=minute_start
      second_end=second_start

      year_first=year_start
      month_first=month_start
      day_first=day_start
      hour_first=hour_start
      minute_first=minute_start 
      second_first=second_start

      write(*,*) 'start',year_start, month_start, day_start, hour_start, minute_start, second_start
      write(*,*) 'end', year_end, month_end, day_end, hour_end, minute_end, second_end

      ! changed declaration to integer (as defined) 2014/09/29
      write(*,*) 'calculate temporal resolution' 
           ! dt = 6.
            dt = 6
      
      print*, 'succesfully init module_wrfchem_lib ...'

      end subroutine init_wrfchem_lib

!---------------------------------------------------------------------------------------------------
!     define variable
!     1. Define variable according input arguments
!     2. Define the attributes of the variable
!---------------------------------------------------------------------------------------------------

      subroutine define_variable( vname, nf_type, ndims, dimids, order )

      implicit none

!     input arguments
      character (len=*), intent(in) :: vname
      integer, intent(in)           :: nf_type
      integer, intent(in)           :: ndims
      integer, intent(in)           :: dimids(ndims)
      character (len=*), intent(in) :: order

!     local arguments
      integer :: status
      character(len=20) :: string 
      integer :: vid

!     define variables
      status = nf_def_var( ncid, vname, nf_type, ndims, dimids, vid )
      if( status /= nf_noerr )  call handle_error( status )

!     define attributes
      status = nf_put_att_int( ncid, vid, 'FieldType', nf_int, 1, 104 )
      if( status /= nf_noerr )  call handle_error( status )

      string = trim(order) 
      status = nf_put_att_text( ncid, vid, 'MemoryOrder', len_trim(string), trim(string) )
      if( status /= nf_noerr )  call handle_error( status )

      string = '-' 
      status = nf_put_att_text( ncid, vid, 'description', len_trim(string), trim(string) ) 
      if( status /= nf_noerr )  call handle_error( status )

      string = '-'
      status = nf_put_att_text( ncid, vid, 'units', len_trim(string), trim(string) )
      if( status /= nf_noerr )  call handle_error( status )

      string = ' ' 
      status = nf_put_att_text( ncid, vid, 'stagger', len_trim(string), trim(string) )
      if( status /= nf_noerr )  call handle_error( status )

      end subroutine define_variable

!---------------------------------------------------------------------------------------------------
!     read scalar data 
!---------------------------------------------------------------------------------------------------

      subroutine wrfchem_readscalar( vname, val )   ! read scalar variable 

      implicit none

!     input arguments
      character (len=*), intent(in) :: vname

!     output arugments
      real, intent(out)        :: val

!     local arguments
      integer :: status
      integer :: varid

!     get variable ID, if no such a variable, return
      status = nf90_inq_varid( ncidi, vname, varid )
      if( status /= nf_noerr ) call handle_error( status )
      write(*,*) ' wrfchem_readscalar: get ',vname,varid

!     read value
      nstt(1:2) = (/ 1, 1 /)
      ncnt(1:2) = (/ 1, 1 /)
      status = nf90_get_var( ncidi, varid, val)
      if( status /= nf_noerr ) call handle_error( status )

      write(*,*)' wrfchem_readscalar: value ',val,status

      end subroutine wrfchem_readscalar

!---------------------------------------------------------------------------------------------------
!     read integer data 
!---------------------------------------------------------------------------------------------------

      subroutine wrfchem_readinteger( vname, val )   ! read integer variable 

      implicit none

!     input arguments
      character (len=*), intent(in) :: vname

!     output arugments
      integer, intent(out)        :: val

!     local arguments
      integer :: status
      integer :: varid

!     get variable ID, if no such a variable, return
      status = nf90_inq_varid( ncidi, vname, varid )
      if( status /= nf_noerr ) call handle_error( status )
      print *,' wrfchem_readint: get ',vname,varid,ncidi

!     read value
      nstt(1:2) = (/ 1, 1 /)
      ncnt(1:2) = (/ 1, 1 /)
      status = nf90_get_var( ncidi, varid, val)
      if( status /= nf_noerr ) call handle_error( status )

      print *,' wrfchem_readint: value ',val,status

      end subroutine wrfchem_readinteger

!---------------------------------------------------------------------------------------------------
!     read one-dimensional data
!---------------------------------------------------------------------------------------------------

      subroutine wrfchem_read2d( vname, val )   ! read 2d variable( a slice of data )

      implicit none

!     input arguments
      character (len=*), intent(in) :: vname

!     output arugments
      real, intent(out)        :: val(nz)

!     local arguments
      integer :: status
      integer :: varid

!     get variable ID, if no such a variable, return
      status = nf_inq_varid( ncidi, vname, varid )
      if( status /= nf_noerr ) call handle_error( status )

!     read value
      nstt(1:2) = (/ 1, 1 /)
      ncnt(1:2) = (/ nz, 1 /)
      status = nf90_get_var( ncidi, varid, val(:), nstt(1:2), ncnt(1:2) )
      if( status /= nf_noerr ) call handle_error( status )

      end subroutine wrfchem_read2d

!---------------------------------------------------------------------------------------------------
!     read three-dimensional data
!---------------------------------------------------------------------------------------------------

      subroutine wrfchem_read3d( vname, val )   ! read 3d variable( a slice of data )

      implicit none

!     input arguments
      character (len=*), intent(in) :: vname

!     output arugments
      real, intent(out)        :: val(nx,ny)

!     local arguments
      integer :: status
      integer :: varid

!     get variable ID, if no such a variable, return
      status = nf_inq_varid( ncidi, vname, varid )
      if( status /= nf_noerr ) call handle_error( status )

!     read value
      nstt(1:3) = (/ 1, 1, 1 /)
      ncnt(1:3) = (/ nx, ny, 1 /)
      status = nf_get_vara_real( ncidi, varid, nstt(1:3), ncnt(1:3), val(:,:) )
      if( status /= nf_noerr ) call handle_error( status )

      end subroutine wrfchem_read3d

!---------------------------------------------------------------------------------------------------
!     write four-dimensional data
!---------------------------------------------------------------------------------------------------

      subroutine wrfchem_write4d( vname, valxs, valxe, valys, valye, vatxs, vatxe, vatys, vatye,it )

      implicit none

!     input arguments
      character (len=*), intent(in) :: vname
      real, intent(in)         :: valxs(ny,nz,nw)
      real, intent(in)         :: valxe(ny,nz,nw)
      real, intent(in)         :: valys(nx,nz,nw)
      real, intent(in)         :: valye(nx,nz,nw)
      real, intent(in)         :: vatxs(ny,nz,nw)
      real, intent(in)         :: vatxe(ny,nz,nw)
      real, intent(in)         :: vatys(nx,nz,nw)
      real, intent(in)         :: vatye(nx,nz,nw)
      integer, intent(in)      :: it

!     output arugments

!     local arguments
      integer :: status
      integer :: varid
      integer :: i
      character(len=100) :: bcname
      character(len=5),DIMENSION(8)  :: varext

      varext = (/'_BXS ', '_BTXS', '_BXE ', '_BTXE', '_BYS ', '_BTYS',&
                    '_BYE ', '_BTYE'/)

      do i = 1, 8
         bcname = trim(vname)//trim(varext(i))
         !write(*,*) 'get variable id in the file : ',bcname
         status = nf90_inq_varid( ncid, trim(bcname), varid )
         if( status /= nf_noerr )  call handle_error( status )
         
         !write(*,*) 'set start and count arrays'
         nstt(1:4) = (/ 1, 1, 1, it /)
        
         if(i <= 4) ncnt(1:4) = (/ ny, nz, nw, 1 /)
         if(i >  4) ncnt(1:4) = (/ nx, nz, nw, 1 /)
         
         !write(*,*) 'write BC values'
         ! west
         if ( i == 1 ) then
            status = nf90_put_var( ncid, varid, valxs(:,:,:), nstt(:), ncnt(:) )
            if( status /= nf_noerr ) call handle_error( status )     
            ! east       
         elseif ( i == 3 ) then
            status = nf90_put_var( ncid, varid, valxe(:,:,:), nstt(:), ncnt(:) )
            if( status /= nf_noerr ) call handle_error( status )  
            ! north
         elseif ( i == 5 ) then
            status = nf90_put_var( ncid, varid, valys(:,:,:), nstt(:), ncnt(:) )
            if( status /= nf_noerr ) call handle_error( status ) 
            ! south
         elseif ( i == 7 ) then
            status = nf90_put_var( ncid, varid, valye(:,:,:), nstt(:), ncnt(:) )
            if( status /= nf_noerr ) call handle_error( status )
         else
            if ( it >= 2 ) then
               nstt(1:4) = (/ 1, 1, 1, it-1 /)
               ! Tendencies 
               ! west  
               if ( i == 2 ) then
                  status = nf90_put_var( ncid, varid, vatxs(:,:,:), nstt(:), ncnt(:) )
                  if( status /= nf_noerr ) call handle_error( status )  
                  ! east        
               elseif ( i == 4 ) then
                  status = nf90_put_var( ncid, varid, vatxe(:,:,:), nstt(:), ncnt(:) )
                  if( status /= nf_noerr ) call handle_error( status ) 
                  ! north
               elseif ( i == 6 ) then
                  status = nf90_put_var( ncid, varid, vatys(:,:,:), nstt(:), ncnt(:) )
                  if( status /= nf_noerr ) call handle_error( status )  
                  ! south
               elseif ( i == 8 ) then
                  status = nf90_put_var( ncid, varid, vatye(:,:,:), nstt(:), ncnt(:) )
                  if( status /= nf_noerr ) call handle_error( status )
               end if
            end if
         end if
         
      end do
      end subroutine wrfchem_write4d 


!---------------------------------------------------------------------------------------------------
!     handle errors produced by calling netCDF functions
!---------------------------------------------------------------------------------------------------

      subroutine handle_error( status )

      implicit none

!     input arguments :
      integer, intent(in) :: status

!     print the error information from processing NETcdf file
      print*, nf_strerror( status )

!     exit from the bconLib
      call exit_wrfchem_lib( flag=1 )

      end subroutine handle_error

!---------------------------------------------------------------------------------------------------
!     exit from bconLib
!---------------------------------------------------------------------------------------------------

      subroutine exit_wrfchem_lib( flag )

!     input arguments
      integer, optional, intent(in) :: flag

!     local arguments
      integer :: status

!     close netCDF file
      if( ncid /= 0 ) status = nf_close( ncid )

!     output information
      if( present(flag) ) then
        select case( flag )
          case( 1 ); print*, 'module_wrfchem_lib: fail to process netCDF file...'
          case( 2 ); print*, 'no such a species to save ...' 
          case default; print*, 'unknown error(s) occurred ...'
        endselect
        stop ' in module_wrfchem_lib ...'
      else
        print*, 'successfully exit from module_wrfchem_lib ...'
      endif 

      end subroutine exit_wrfchem_lib

      end module module_wrfchem_lib

!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
