
      module module_GEOS_lib

      use netcdf      
      implicit none

!     include files
      include 'netcdf.inc'

!     public readprocedures
      public  :: init_GEOS_lib, mozart_interpolate4d, exit_mozart_lib

!     private procedures
      ! [PSU]: removed translate_radm2_mozart in process of removing refs to raqms processing 2014/09/29
      private :: handle_error_moz!, translate_radm2_mozart

!     public variables
      integer, public :: nlon, nlat, nlev                    ! data dimension length
      integer, public  :: year_start_m, month_start_m, day_start_m
      integer, public  :: hour_start_m, minute_start_m, second_start_m
      integer, public  :: year_end_m, month_end_m, day_end_m
      integer, public  :: hour_end_m
      integer, public  :: dt_m, ntime_m

!     variables for reading netCDF data
      integer, private :: ncid                               ! netCDF file ID
      integer, private :: nstt(4), ncnt(4)                   ! start index and counter array
      integer, private :: start_index = 1                    ! monthly mean data
      integer, private  :: ntime      ! data dimension length
      integer, private :: nscalar, nstring           ! data dimension length
      integer, private, allocatable :: date(:), datesec(:) 

!     variables used by interpolation
!     [PSU]:  note that the only variables currently used are ix, jy, kz, ax, by,
!             and the cz.. series (czxs, czxe, czys, czye),
!             others are calculated but not used
      integer, private, allocatable :: ix(:,:,:), jy(:,:,:), kz(:,:,:,:)  ! index used by interpolation
      real, private, allocatable :: ax(:,:,:), by(:,:,:), cz(:,:,:,:)           ! weight coef. all domain
      real, private, allocatable :: axxs(:,:,:), byxs(:,:,:), czxs(:,:,:,:)     ! weight coef. west
      real, private, allocatable :: axxe(:,:,:), byxe(:,:,:), czxe(:,:,:,:)     ! weight coef. east
      real, private, allocatable :: axys(:,:,:), byys(:,:,:), czys(:,:,:,:)     ! weight coef. south
      real, private, allocatable :: axye(:,:,:), byye(:,:,:), czye(:,:,:,:)     ! weight coef. north

      real, private, allocatable :: mozval(:,:,:), mozval_temp(:,:,:), tmpval(:,:,:), mozval2(:,:,:,:)

      contains

!---------------------------------------------------------------------------------------------------
!     initialize netcdf file
!---------------------------------------------------------------------------------------------------

      subroutine init_GEOS_lib( mozart_fn, ptop_wrf, x2d, y2d, znu, nx, ny, nz )

      implicit none

!     input arguments
      integer, intent(in)      :: nx, ny, nz
      character(*), intent(in) :: mozart_fn 
      real, intent(in)         :: ptop_wrf
      real, intent(in)         :: x2d(nx,ny), y2d(nx,ny)
      real, intent(in)         :: znu(nz)

!     local argumets
      integer :: status
      integer :: dimid, varid, timelen
      integer :: i, j, k, n, l, nou
      integer :: vid
      real :: ptop_wrf3   !ptop_wrf in units of Pa*1e3 - added March 2015
      real, allocatable :: x1d(:), y1d(:)
      real, allocatable :: hyam(:), hybm(:)
      real, allocatable :: p_moz(:), p_wrf(:)
      real, allocatable :: ps_moz(:,:), ps_wrf(:,:)
      real, allocatable :: dateCT(:,:)
      character(100)    :: dtstrings, dtstringe, dtstringf

      write(*,*) 'open the input netCDF file ',mozart_fn
      status = nf_open( mozart_fn, nf_nowrite, ncid )
      if( status /= nf_noerr ) call handle_error_moz( status )

      write(*,*) 'get the longitude dimension length of GEOS-CHEM'
      status = nf_inq_dimid( ncid, 'lon', dimid )
      if( status /= nf_noerr )  call handle_error_moz( status )

      status = nf_inq_dimlen( ncid, dimid, nlon )
      if( status /= nf_noerr )  call handle_error_moz( status )

      write(*,*) 'get the latitude dimension length of the GEOS-CHEM'
      status = nf_inq_dimid( ncid, 'lat', dimid )
      if( status /= nf_noerr )  call handle_error_moz( status )

      status = nf_inq_dimlen( ncid, dimid, nlat )
      if( status /= nf_noerr )  call handle_error_moz( status )

!     get the vertical dimension length of the GEOS-Chem data
      status = nf_inq_dimid( ncid, 'lev', dimid )
      if( status /= nf_noerr )  call handle_error_moz( status )

      status = nf_inq_dimlen( ncid, dimid, nlev )
      if( status /= nf_noerr )  call handle_error_moz( status )

      write(*,*) 'allocate sigma pressure variable hyam'
      allocate( hyam(nlev), hybm(nlev) )
      write(*,*) 'read sigma levels',nlev   ![PSU]: not really, but we will treat as such

      status = nf_inq_varid( ncid, 'lev', varid )
      if( status /= nf_noerr )  call handle_error_moz( status )

      nstt(1:4) = (/ 1, 1, 1, 1 /)
      ncnt(1:4) = (/ 1, 1, nlev, 1/)
      
!      status = nf_get_vara_real( ncid, varid, nstt(1:4), ncnt(1:4), hyam(1:nlev) )

      status = nf_get_var_real( ncid, varid, hyam(1:nlev) )
      if( status /= nf_noerr )  call handle_error_moz( status )

      hyam(:)= (hyam(:)*(101325.-1.)+1.)*1e3

!      write(*,*) 'Test sigma pressure level : ',hyam(1:nlev)  !commented to avoid excessive displays 2014/09/29
!      status = nf_inq_varid( ncid, 'hybm', varid )
!      if( status /= nf_noerr )  call handle_error_moz( status )

!      status = nf_get_var_real( ncid, varid, hybm )
!      if( status /= nf_noerr )  call handle_error_moz( status )

      write(*,*) 'read longitudes'
      allocate( x1d(nlon), y1d(nlat) )

      status = nf_inq_varid( ncid, 'lon', varid )
      if( status /= nf_noerr )  call handle_error_moz( status )

      status = nf_get_var_real( ncid, varid, x1d )
      if( status /= nf_noerr )  call handle_error_moz( status )
!     [PSU]: changed 2014/09/25 for CMS longitude vector
!            in input file, it is [180.0, 185.0,...,355.0, 0.0, 5.0,...,175.0]
      where( x1d > 179.0 )
        x1d = x1d - 360.0
      endwhere

!     read latitudes
      status = nf_inq_varid( ncid, 'lat', varid )
      if( status /= nf_noerr )  call handle_error_moz( status )

      status = nf_get_var_real( ncid, varid, y1d )
      if( status /= nf_noerr )  call handle_error_moz( status )

      write(*,*) 'read number of dates in the file'
      status = nf_inq_dimid( ncid, 'time', varid )
      if( status /= nf_noerr )  call handle_error_moz( status )

      status = nf_inq_dimlen( ncid, varid, nstring )
      if( status /= nf_noerr ) call handle_error_moz( status )

      ntime_m = nstring
      ! hours between 1985 and 2010: 219144

      allocate ( dateCT(4,nstring) )

      write(*,*) 'check time dimensions : ',timelen, nstring

      ! [PSU]: note that these are not used!
      !    The variables day_start and hour_start in the main module are used to
      !    find the starting index into the GEOS/CMS concatenated files
      year_start_m   = 2010! dateCT(1,1)       !floor(float(date(1)/10000))
      month_start_m  = 1! dateCT(2,1)       !floor(float(date(1)/100-year_start_m*100))
      day_start_m    = 0! dateCT(3,1)       !float(date(1)-year_start_m*10000-month_start_m*100)
      hour_start_m   = 0! dateCT(4,1)       !datesec(1)/3600
      
      year_end_m   = 2010! dateCT(1,ntime)    !floor(float(date(nstring)/10000))
      month_end_m  = 12! dateCT(2,ntime)    ! floor(float(date(nstring)/100-year_end_m*100))
      day_end_m    = 31! dateCT(3,ntime)    !float(date(nstring)-year_end_m*10000-month_end_m*100)
      hour_end_m   = 21! dateCT(4,ntime)    !datesec(nstring)/3600
      
      !write(*,*) 'Dates in GEOS-CHEM file - Start : ',dateCT(:,1)
      !write(*,*) 'Dates in GEOS-CHEM file - End : ',dateCT(:,ntime)

      ! temporal resolution of mozart file
      ! [PSU]: changed to integer (as defined) 2014/09/29
      !dt_m         = 3.             !  abs(datesec(2)/3600 - datesec(1)/3600)
      dt_m         = 3

!     allocate memory space to store interpolation coef.
!     [PSU]: again, note that only ax, by, ix, jy, kz, and cz.. are used
      allocate( ax(0:1,nx,ny), by(0:1,nx,ny), cz(0:1,nx,ny,nz) )
      allocate( axxs(0:1,nx,ny), byxs(0:1,nx,ny), czxs(0:1,nx,ny,nz) )
      allocate( axxe(0:1,nx,ny), byxe(0:1,nx,ny), czxe(0:1,nx,ny,nz) )
      allocate( axys(0:1,nx,ny), byys(0:1,nx,ny), czys(0:1,nx,ny,nz) )
      allocate( axye(0:1,nx,ny), byye(0:1,nx,ny), czye(0:1,nx,ny,nz) )
      allocate( ix(0:1,nx,ny), jy(0:1,nx,ny), kz(0:1,nx,ny,nz) )

!     horizontal interpolation coefs.

      write(*,*) 'x2d : ',x2d(1:10,23)

      ! all domain
      do i = 1, nx
      do j = 1, ny 
         do n = 1, nlon
            if( x2d(i,j) < x1d(n) ) exit
         enddo
         
         ix(0,i,j) = min( nlon-1, max(n-1,1) )
         ix(1,i,j) = ix(0,i,j) + 1
         ax(0,i,j) = ( x2d(i,j) - x1d(ix(0,i,j)) )/( x1d(ix(1,i,j)) - x1d(ix(0,i,j)) )
         ax(1,i,j) = 1.0 - ax(0,i,j)
         
         do n = 1, nlat
            if( y2d(i,j) < y1d(n) ) exit
         enddo
         
         jy(0,i,j) = min( nlat-1, max(n-1,1) )
         jy(1,i,j) = jy(0,i,j) + 1
         by(0,i,j) = ( y2d(i,j) - y1d(jy(0,i,j)) )/( y1d(jy(1,i,j)) - y1d(jy(0,i,j)) )
         by(1,i,j) = 1.0 - by(0,i,j)
         
      enddo
   end do

      ![PSU]: note that none of the following are used in the interpolation and could be omitted

      ! west
      i = 1
      do j = 1, ny 
         do n = 1, nlon
            if( x2d(i,j) < x1d(n) ) exit
         enddo
         
         axxs(0,i,j) = ( x2d(i,j) - x1d(ix(0,i,j)) )/( x1d(ix(1,i,j)) - x1d(ix(0,i,j)) )
         axxs(1,i,j) = 1.0 - axxs(0,i,j)
         
         do n = 1, nlat
            if( y2d(i,j) < y1d(n) ) exit
         enddo
         
         byxs(0,i,j) = ( y2d(i,j) - y1d(jy(0,i,j)) )/( y1d(jy(1,i,j)) - y1d(jy(0,i,j)) )
         byxs(1,i,j) = 1.0 - byxs(0,i,j)
         
      enddo

      ! east
      i = nx
      do j = 1, ny 
         
         do n = 1, nlon
            if( x2d(i,j) < x1d(n) ) exit
         enddo
         
         axxe(0,i,j) = ( x2d(i,j) - x1d(ix(0,i,j)) )/( x1d(ix(1,i,j)) - x1d(ix(0,i,j)) )
         axxe(1,i,j) = 1.0 - axxe(0,i,j)
         
         do n = 1, nlat
            if( y2d(i,j) < y1d(n) ) exit
         enddo
         
         byxe(0,i,j) = ( y2d(i,j) - y1d(jy(0,i,j)) )/( y1d(jy(1,i,j)) - y1d(jy(0,i,j)) )
         byxe(1,i,j) = 1.0 - byxe(0,i,j)
         
      enddo

      ! north
      j = ny
      do i = 1, nx
         
         do n = 1, nlon
            if( x2d(i,j) < x1d(n) ) exit
         enddo
         
         axye(0,i,j) = ( x2d(i,j) - x1d(ix(0,i,j)) )/( x1d(ix(1,i,j)) - x1d(ix(0,i,j)) )
         axye(1,i,j) = 1.0 - axye(0,i,j)
         
         do n = 1, nlat
            if( y2d(i,j) < y1d(n) ) exit
         enddo
         
         byye(0,i,j) = ( y2d(i,j) - y1d(jy(0,i,j)) )/( y1d(jy(1,i,j)) - y1d(jy(0,i,j)) )
         byye(1,i,j) = 1.0 - byye(0,i,j)
         
      enddo
 
      ! south
      j = 1 
      do i = 1, nx

        do n = 1, nlon
          if( x2d(i,j) < x1d(n) ) exit
        enddo

        axys(0,i,j) = ( x2d(i,j) - x1d(ix(0,i,j)) )/( x1d(ix(1,i,j)) - x1d(ix(0,i,j)) )
        axys(1,i,j) = 1.0 - axys(0,i,j)

        do n = 1, nlat
          if( y2d(i,j) < y1d(n) ) exit
        enddo

        byys(0,i,j) = ( y2d(i,j) - y1d(jy(0,i,j)) )/( y1d(jy(1,i,j)) - y1d(jy(0,i,j)) )
        byys(1,i,j) = 1.0 - byys(0,i,j)

      enddo
      ![PSU] end of derivation of horizontal interpolation coefficients not used

      write(*,*) 'read surface pressure and interpolate ...'
      allocate( ps_moz(nlon,nlat), ps_wrf(nx,ny) )

      status = nf_inq_varid( ncid, 'PEDGE_S__PSURF', varid )
      if( status /= nf_noerr )  call handle_error_moz( status )

      nstt(1:3) = (/ 1, 1, start_index /)  ![PSU]:  start_index=1 is first time step of composite CMS CO2 file
      ncnt(1:3) = (/ nlon, nlat, 1/)
      status = nf_get_vara_real( ncid, varid, nstt(1:3), ncnt(1:3), ps_moz )
      if( status /= nf_noerr ) call handle_error_moz( status )

      write(*,*) 'Test horizontal interpolation coefficients : ',ix(0,23,23),ax(0,23,23),by(0,23,23),jy(0,23,23)
   
      ![PSU]: for clarification, ps_moz is CMS surface pressure and
      !                          ps_wrf is CMS surface pressure interpolated to the WRF grid
      ps_moz(:,:)=ps_moz(:,:)*1e2
      ! Constant pressure at the surface (for now)  
      write(*,*) start_index,ps_moz(23,23)

      do j = 1, ny 
      do i = 1, nx 
        ps_wrf(i,j) = 1e3*(ps_moz(ix(0,i,j),jy(0,i,j))*ax(1,i,j)*by(1,i,j) +   &
                      ps_moz(ix(0,i,j),jy(1,i,j))*ax(1,i,j)*by(0,i,j) +   &
                      ps_moz(ix(1,i,j),jy(0,i,j))*ax(0,i,j)*by(1,i,j) +   &
                      ps_moz(ix(1,i,j),jy(1,i,j))*ax(0,i,j)*by(0,i,j))
      enddo
      enddo

      write(*,*) 'Examples of interpolated surface pressure : ',ps_wrf(1:5,2)

      write(*,*) 'vertical interpolation coefs'
      allocate( p_moz(nlev), p_wrf(nz) )
      p_moz(:)=0.
      p_wrf(:)=0.

      ! make ptop_wrf same units as ps_wrf [Pa * 1e3]  - March 2015 bug fix
      ! variable used in p_wrf calculations below changed from ptop_wrf to ptop_wrf3
      ! get_var in main module changed to keep ptop in Pa units
      ptop_wrf3 = ptop_wrf * 1e3
      write(*,*) 'ptop_wrf for p_wrf calculation ', ptop_wrf3
      
      ![PSU]: note that before September 2015, only kz[0,i,j,k] is used in the interpolate subroutine
      !   This has been revised to use kz[0:1,i,j,k] and the corresponding czxs,czxe,czys,czye vertical
      !   interpolation coefficients
      !west
      i = 1
      ! write(*,*) hyam  !commented 2014/09/29 to avoid excessive displays
      do j = 1, ny
         !write(*,*) hyam ps_wrf(i,j)
         !p_moz(:) = log( hybm(:)*ps_wrf(i,j) + hyam(:)*1.0e5 )
         do l=1,nlev
            if (hyam(l)+ps_wrf(i,j)-hyam(1).gt.1.) then
               p_moz(l) = log( hyam(l)+ps_wrf(i,j)-hyam(1) )
               !p_moz(l) = log( hyam(l)+1013.25-hyam(1) )
            else
               p_moz(l)=0.
            endif
         enddo
         p_wrf(:) = log( znu(:)*ps_wrf(i,j) + (1.0-znu(:))*ptop_wrf3 )

         
         do k = 1, nz
            
            do n = 1, nlev 
               if( p_wrf(k) > p_moz(n) ) exit
            enddo
            
            kz(0,i,j,k) = min( nlev-1, max(n-1,1) )
            kz(1,i,j,k) = kz(0,i,j,k) + 1
            czxs(0,i,j,k) = ( p_wrf(k) - p_moz(kz(0,i,j,k)) )/( p_moz(kz(1,i,j,k)) - p_moz(kz(0,i,j,k)) )
            czxs(1,i,j,k) = 1.0 - czxs(0,i,j,k)
            !write(*,*) 'chouchou : ',p_moz(kz(1,i,j,k)), p_moz(kz(0,i,j,k))
            !write(*,*) 'chouchou2 : ',czxs(0,i,j,k),czxs(1,i,j,k)
         enddo
            
      end do
       ! write(*,*) 'vertical interpolation vector for WRF : ',p_wrf(:)
        !write(*,*) 'vertical interpolation vector for CMS : ',p_moz(:)
        !write(*,*) 'kz levels : ',kz(0,1,25,:)  !commented to reduce excessive displays 2014/09/29
 

      ! east
      i = nx
      do j = 1, ny
         
         !p_moz(:) = log( hybm(:)*ps_wrf(i,j) + hyam(:)*1.0e5 )
         p_moz(:) = log( hyam(:)+ps_wrf(i,j)-hyam(1) )
         p_wrf(:) = log( znu(:)*ps_wrf(i,j) + (1.0-znu(:))*ptop_wrf3 )
         
         do k = 1, nz
            
            do n = 1, nlev 
               if( p_wrf(k) > p_moz(n) ) exit
            enddo

            kz(0,i,j,k) = min( nlev-1, max(n-1,1) )
            kz(1,i,j,k) = kz(0,i,j,k) + 1
            czxe(0,i,j,k) = ( p_wrf(k) - p_moz(kz(0,i,j,k)) )/( p_moz(kz(1,i,j,k)) - p_moz(kz(0,i,j,k)) )
            czxe(1,i,j,k) = 1.0 - czxe(0,i,j,k)
            
         enddo
      enddo

      ! north
      j = ny
      do i = 1, nx
         
         p_moz(:) = log( hyam(:)+ps_wrf(i,j)-hyam(1) )
         !p_moz(:) = log( hybm(:)*ps_wrf(i,j) + hyam(:)*1.0e5 )
         p_wrf(:) = log( znu(:)*ps_wrf(i,j) + (1.0-znu(:))*ptop_wrf3 )
         
         do k = 1, nz
            
            do n = 1, nlev 
               if( p_wrf(k) > p_moz(n) ) exit
            enddo
            
            kz(0,i,j,k) = min( nlev-1, max(n-1,1) )
            kz(1,i,j,k) = kz(0,i,j,k) + 1
            czye(0,i,j,k) = ( p_wrf(k) - p_moz(kz(0,i,j,k)) )/( p_moz(kz(1,i,j,k)) - p_moz(kz(0,i,j,k)) )
            czye(1,i,j,k) = 1.0 - czye(0,i,j,k)
            
         enddo
         
      enddo
      
      ! south
      j = 1
      do i = 1, nx

        !p_moz(:) = log( hybm(:)*ps_wrf(i,j) + hyam(:)*1.0e5 )
        p_moz(:) = log( hyam(:)+ps_wrf(i,j)-hyam(1) )
        p_wrf(:) = log( znu(:)*ps_wrf(i,j) + (1.0-znu(:))*ptop_wrf3 )

        do k = 1, nz

          do n = 1, nlev 
            if( p_wrf(k) > p_moz(n) ) exit
          enddo

          kz(0,i,j,k) = min( nlev-1, max(n-1,1) )
          kz(1,i,j,k) = kz(0,i,j,k) + 1
          czys(0,i,j,k) = ( p_wrf(k) - p_moz(kz(0,i,j,k)) )/( p_moz(kz(1,i,j,k)) - p_moz(kz(0,i,j,k)) )
          if( czys(0,i,j,k) > 1.)  czys(0,i,j,k) = 1.
          czys(1,i,j,k) = 1.0 - czys(0,i,j,k)

        enddo

      enddo

      write(*,*) 'release memory space not used'
      deallocate( dateCT, x1d, y1d, hyam, hybm, ps_moz, ps_wrf, p_moz, p_wrf )

      write(*,*) 'allocate memory space for reading'
      allocate( mozval(nlon,nlat,nlev), mozval_temp(nlon,nlat,nlev), tmpval(nlon,nlat,nlev), mozval2(nlon,nlat,nlev,2) )

      write(*,*) 'End of subroutine init_CT_lib'
      end subroutine init_GEOS_lib

!---------------------------------------------------------------------------------------------------
!     interpolate four-dimensional field ... 
!---------------------------------------------------------------------------------------------------

      subroutine mozart_interpolate4d( wrfspn, wrfxs, wrfxe, wrfys, wrfye, nx, ny, nz, nw, i_con, i_carbtrac  )
      use netcdf
      implicit none

!     input arguments
      integer, intent(in)      :: nx, ny, nz, nw     ! dimensions
      integer, intent(in)      :: i_con, i_carbtrac
     ! [PSU]: character length changed to agree with specification in main_nc_wrfchem.f90  2014/09/29
      !character(50), intent(in) :: wrfspn         ! wrfchem species' name
      character(100), intent(in) :: wrfspn         ! wrfchem species' name

!     output arguments, [PSU]: but note that they are passed as t1..(ny,nx,nz,1) variables from the main module
!                              where they are allocated
      ! real, intent(out) :: wrfval(nx,ny,nz)         ! wrfchem vmr(ppm)
      real, intent(out) :: wrfxs(ny,nz,nw)         ! wrfchem vmr(ppm)
      real, intent(out) :: wrfxe(ny,nz,nw)         ! wrfchem vmr(ppm)
      real, intent(out) :: wrfys(nx,nz,nw)         ! wrfchem vmr(ppm)
      real, intent(out) :: wrfye(nx,nz,nw)         ! wrfchem vmr(ppm)

!     local arguments
      integer :: status, varid, you, youmax, new_i_carbtrac
     
      integer :: nspec
      !character(20) :: mozspn(20)   ![PSU]: changed to be non-dimensioned for CMS (may be different for CT)
      character(20) :: mozspn
      real :: scale, chouchou

      integer :: i, j, k, n, p,q,r
      integer, dimension(4) :: tchouk, tchak

      nspec=1  ! Bug in the code (terrible coding... )
      write(*,*) 'Process wrfchem species : ', trim(wrfspn)
      write(*,*) 'Number of species : ', nspec

      wrfxs(:,:,:) = 385.
      wrfxe(:,:,:) = 385.
      wrfys(:,:,:) = 385.
      wrfye(:,:,:) = 385.

!     read MOZART data  [PSU]: in this case the GEOS-Chem/CMS CO2 mole fractions
!     [PSU]:  Modifications made September 2015 to bring time steps in line with those used 
!             for CarbonTracker (something lost in earlier translation?)
!             'it' (wrfbdy time step) and 'it_CT' (increment in GEOS-Chem file) in main_bc_wrfchem
!             become known as 'i_con' and 'i_carbtrac' here in the GEOS module
!             Modifications made March 2015 to simplify getvar for CMS only
!             some of the loops below could be removed, as well... [someday...]
!             Some logic retained by commented from earlier CT implementations.
!             

      mozval(:,:,:) = 0.
      !mozspn(1) = 'IJ_AVG_S__NOx'
      ! [PSU]: changed to be non-dimensioned for GEOS
      mozspn = 'IJ_AVG_S__NOx'
   

      write(*,*) 'number of species (should be one with CO2)',nspec
      ! [PSU]: youmax condition specifically for CMS:
      if (i_con.eq.1) then    ![PSU]: this is the variable called it in main (youmax needed for CT)
        youmax=1
      else
        youmax=2
      endif
    
      do you=1,youmax
        do n = 1, nspec
          !do i=1,1   ![PSU]:  This loop for CT where multiple tracers used to build CO2 bdy
            tmpval(:,:,:)=0.
            new_i_carbtrac = i_carbtrac+you-1
            !write(*,*) 'Read and GEOS/CMS CO2 : ', trim(mozspn(i)), new_i_carbtrac
            write(*,*) 'Read GEOS/CMS CO2 : ', trim(mozspn), new_i_carbtrac
            !status = nf90_inq_varid( ncid, mozspn(i), varid )   ![PSU]: for CT
            status = nf90_inq_varid( ncid, mozspn, varid )
            if( status /= nf_noerr )  call handle_error_moz( status )
            !do p=1,nlon
            !   do q=1,nlat
            !      do r=1,nlev
            !         !write(*,*) p,q,r,new_i_carbtrac
            !          tchouk(1:4) = (/ p, q, r, new_i_carbtrac /)
            !          tchak(1:4) = (/ 1, 1, 1, 1 /)
            !          status = nf_get_vara_real( ncid, varid, tchouk(:), tchak(:), chouchou )
            !          tmpval(p,q,r)=real(chouchou,kind=8)
            nstt(1:4) = (/ 1, 1, 1, new_i_carbtrac /)
            ncnt(1:4) = (/ nlon, nlat, nlev, 1 /)
            status = nf_get_vara_real( ncid, varid, nstt(:), ncnt(:), tmpval(:,:,:) )
                      if( status /= nf_noerr ) call handle_error_moz( status )
            !      enddo
            !    enddo
            !enddo
            mozval(:,:,:) = mozval(:,:,:) + tmpval(:,:,:)*1e-3   ! change unit to ppm
          !enddo
          write(*,*) 'Some values just after reading the CMS CO2 : ', mozval(23:26,23,1)
        enddo
      enddo
      ! [PSU]: following should acknowledge that youmax is an integer
      !mozval(:,:,:)=mozval(:,:,:)/youmax
      mozval(:,:,:)=mozval(:,:,:)/float(youmax)

      ! [PSU]: following is for use with CarbonTracker (I think)
      !do p=1,nlon
      !   do q=1,nlat
      !       do r=1,nlev-2
      !          mozval(p,q,r)=(mozval(p,q,r+1)+mozval(p,q,r+2)+mozval(p,q,r))/3.
      !       enddo
      !       mozval(p,q,nlev-1)=(mozval(p,q,nlev-2)+mozval(p,q,nlev-1)+mozval(p,q,nlev))/3.
      !       mozval(p,q,nlev)=mozval(p,q,nlev-1)
      !   enddo
      !enddo
      !do k=1,nlev
      !    mozval(:,:,k)=mozval_temp(:,:,nlev-k+1)
      !enddo

      write(*,*) 'Examples final :: ', mozval(2,3,1),mozval(3,4,5)
      !if( scale >= 0.0 ) mozval(:,:,:) = mozval(:,:,:)*scale

!     interpolate in three dimensions
!     [PSU]: note that only ix(0,i,j), jy(0,i,j), are used here horizontally
!     [PSU]: vertical interpolation using kz and the cz.. coefficients added September 2015 for CMS

   
      ! west
      i= 1
      do k = 1, nz
         do j = 1, ny
            wrfxs(j,k,:)  = mozval(ix(0,i,j),jy(0,i,j),kz(0,i,j,k))*czxs(1,i,j,k) + &
                            mozval(ix(0,i,j),jy(0,i,j),kz(1,i,j,k))*czxs(0,i,j,k)
            !wrfxs(j,k,:)  = mozval(ix(0,i,j),jy(0,i,j),kz(0,i,j,k))!*axxs(1,i,j)*byxs(1,i,j)*czxs(1,i,j,k) + &
                            !mozval(ix(0,i,j),jy(0,i,j),kz(1,i,j,k))*axxs(1,i,j)*byxs(1,i,j)*czxs(0,i,j,k) + &
                            !mozval(ix(0,i,j),jy(1,i,j),kz(0,i,j,k))*axxs(1,i,j)*byxs(0,i,j)*czxs(1,i,j,k) + &
                            !mozval(ix(0,i,j),jy(1,i,j),kz(1,i,j,k))*axxs(1,i,j)*byxs(0,i,j)*czxs(0,i,j,k) + &
                            !mozval(ix(1,i,j),jy(0,i,j),kz(0,i,j,k))*axxs(0,i,j)*byxs(1,i,j)*czxs(1,i,j,k) + &
                            !mozval(ix(1,i,j),jy(0,i,j),kz(1,i,j,k))*axxs(0,i,j)*byxs(1,i,j)*czxs(0,i,j,k) + &
                            !mozval(ix(1,i,j),jy(1,i,j),kz(0,i,j,k))*axxs(0,i,j)*byxs(0,i,j)*czxs(1,i,j,k) + &
                            !mozval(ix(1,i,j),jy(1,i,j),kz(1,i,j,k))*axxs(0,i,j)*byxs(0,i,j)*czxs(0,i,j,k)
            !write(*,*) mozval(ix(0,i,j),jy(0,i,j),kz(0,i,j,k)),axxs(1,i,j),byxs(1,i,j),czxs(1,i,j,k)
         enddo
      enddo
      

      ! east
      i = nx
      do k = 1, nz
         do j = 1, ny
            wrfxe(j,k,:)  = mozval(ix(0,i,j),jy(0,i,j),kz(0,i,j,k))*czxe(1,i,j,k) + &
                            mozval(ix(0,i,j),jy(0,i,j),kz(1,i,j,k))*czxe(0,i,j,k)
            !wrfxe(j,k,:)  = mozval(ix(0,i,j),jy(0,i,j),kz(0,i,j,k))!*axxe(1,i,j)*byxe(1,i,j)*czxe(1,i,j,k) + &
                            !mozval(ix(0,i,j),jy(0,i,j),kz(1,i,j,k))*axxe(1,i,j)*byxe(1,i,j)*czxe(0,i,j,k) + &
                            !mozval(ix(0,i,j),jy(1,i,j),kz(0,i,j,k))*axxe(1,i,j)*byxe(0,i,j)*czxe(1,i,j,k) + &
                            !mozval(ix(0,i,j),jy(1,i,j),kz(1,i,j,k))*axxe(1,i,j)*byxe(0,i,j)*czxe(0,i,j,k) + &
                            !mozval(ix(1,i,j),jy(0,i,j),kz(0,i,j,k))*axxe(0,i,j)*byxe(1,i,j)*czxe(1,i,j,k) + &
                            !mozval(ix(1,i,j),jy(0,i,j),kz(1,i,j,k))*axxe(0,i,j)*byxe(1,i,j)*czxe(0,i,j,k) + &
                            !mozval(ix(1,i,j),jy(1,i,j),kz(0,i,j,k))*axxe(0,i,j)*byxe(0,i,j)*czxe(1,i,j,k) + &
                            !mozval(ix(1,i,j),jy(1,i,j),kz(1,i,j,k))*axxe(0,i,j)*byxe(0,i,j)*czxe(0,i,j,k)

         enddo
      enddo

      ! north
      j = ny
      do k = 1, nz
         do i = 1, nx
            wrfye(i,k,:)  = mozval(ix(0,i,j),jy(0,i,j),kz(0,i,j,k))*czye(1,i,j,k) + &
                            mozval(ix(0,i,j),jy(0,i,j),kz(1,i,j,k))*czye(0,i,j,k) 
            !wrfye(i,k,:)  = mozval(ix(0,i,j),jy(0,i,j),kz(0,i,j,k))!*axye(1,i,j)*byye(1,i,j)*czye(1,i,j,k) + &
                            !mozval(ix(0,i,j),jy(0,i,j),kz(1,i,j,k))*axye(1,i,j)*byye(1,i,j)*czye(0,i,j,k) + &
                            !mozval(ix(0,i,j),jy(1,i,j),kz(0,i,j,k))*axye(1,i,j)*byye(0,i,j)*czye(1,i,j,k) + &
                            !mozval(ix(0,i,j),jy(1,i,j),kz(1,i,j,k))*axye(1,i,j)*byye(0,i,j)*czye(0,i,j,k) + &
                            !mozval(ix(1,i,j),jy(0,i,j),kz(0,i,j,k))*axye(0,i,j)*byye(1,i,j)*czye(1,i,j,k) + &
                            !mozval(ix(1,i,j),jy(0,i,j),kz(1,i,j,k))*axye(0,i,j)*byye(1,i,j)*czye(0,i,j,k) + &
                            !mozval(ix(1,i,j),jy(1,i,j),kz(0,i,j,k))*axye(0,i,j)*byye(0,i,j)*czye(1,i,j,k) + &
                            !mozval(ix(1,i,j),jy(1,i,j),kz(1,i,j,k))*axye(0,i,j)*byye(0,i,j)*czye(0,i,j,k)
            
         enddo
      enddo

      ! south
      j = 1
      do k = 1, nz
         do i = 1, nx
            wrfys(i,k,:)  = mozval(ix(0,i,j),jy(0,i,j),kz(0,i,j,k))*czys(1,i,j,k) + &
                            mozval(ix(0,i,j),jy(0,i,j),kz(1,i,j,k))*czys(0,i,j,k) 
            !wrfys(i,k,:)  = mozval(ix(0,i,j),jy(0,i,j),kz(0,i,j,k))!*axys(1,i,j)*byys(1,i,j)*czys(1,i,j,k) + &
                            !mozval(ix(0,i,j),jy(0,i,j),kz(1,i,j,k))*axys(1,i,j)*byys(1,i,j)*czys(0,i,j,k) + &
                            !mozval(ix(0,i,j),jy(1,i,j),kz(0,i,j,k))*axys(1,i,j)*byys(0,i,j)*czys(1,i,j,k) + &
                            !mozval(ix(0,i,j),jy(1,i,j),kz(1,i,j,k))*axys(1,i,j)*byys(0,i,j)*czys(0,i,j,k) + &
                            !mozval(ix(1,i,j),jy(0,i,j),kz(0,i,j,k))*axys(0,i,j)*byys(1,i,j)*czys(1,i,j,k) + &
                            !mozval(ix(1,i,j),jy(0,i,j),kz(1,i,j,k))*axys(0,i,j)*byys(1,i,j)*czys(0,i,j,k) + &
                            !mozval(ix(1,i,j),jy(1,i,j),kz(0,i,j,k))*axys(0,i,j)*byys(0,i,j)*czys(1,i,j,k) + &
                            !mozval(ix(1,i,j),jy(1,i,j),kz(1,i,j,k))*axys(0,i,j)*byys(0,i,j)*czys(0,i,j,k)

         enddo
      enddo

      end subroutine mozart_interpolate4d

!---------------------------------------------------------------------------------------------------
!     handle errors produced by calling netCDF functions
!---------------------------------------------------------------------------------------------------

      subroutine handle_error_moz( status )

      implicit none

!     input arguments :
      integer, intent(in) :: status

!     print the error information from processing NETcdf file
      print*, nf_strerror( status )

!     exit from the bconLib
      call exit_mozart_lib( flag=1 )

      end subroutine handle_error_moz

!---------------------------------------------------------------------------------------------------
!     exit from module_mozart_lib 
!---------------------------------------------------------------------------------------------------

      subroutine exit_mozart_lib( flag )

!     input arguments
      integer, optional, intent(in) :: flag

!     local arguments
      integer :: status

!     release memory space
      if( allocated(ix) ) deallocate( ix )
      if( allocated(jy) ) deallocate( jy )
      if( allocated(kz) ) deallocate( kz )
      if( allocated(ax) ) deallocate( ax )
      if( allocated(by) ) deallocate( by )
      if( allocated(cz) ) deallocate( cz )
      if( allocated(axxs) ) deallocate( axxs )
      if( allocated(byxs) ) deallocate( byxs )
      if( allocated(czxs) ) deallocate( czxs )
      if( allocated(axxe) ) deallocate( axxe )
      if( allocated(byxe) ) deallocate( byxe )
      if( allocated(czxe) ) deallocate( czxe )
      if( allocated(axys) ) deallocate( axys )
      if( allocated(byys) ) deallocate( byys )
      if( allocated(czys) ) deallocate( czys )
      if( allocated(axye) ) deallocate( axye )
      if( allocated(byye) ) deallocate( byye )
      if( allocated(czye) ) deallocate( czye )

      if( allocated(mozval) ) deallocate( mozval )
      if( allocated(tmpval) ) deallocate( tmpval )

!     close netCDF file
      if( ncid /= 0 ) status = nf_close( ncid )

!     output information
      if( present(flag) ) then
        select case( flag )
          case( 1 ); print*, 'fail to process netCDF file...'
          case default; print*, 'unknown error(s) occurred ...'
        endselect
        stop ' in module_GEOS_lib ...'
      else
        print*, 'successfully exit from module_GEOS_lib ...'
      endif

      end subroutine exit_mozart_lib

!  [PSU]: large amount of dead code removed 2014/09/29 -- see original code if interested


      end module module_GEOS_lib
