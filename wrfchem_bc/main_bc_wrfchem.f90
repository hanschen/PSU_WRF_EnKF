program main_bc_wrfchem

  !<DESCRIPTION>
  !
  ! Program wrfchembc.
  !
  ! Authors:  Rainer Schmitz (University of Chile - Santiago, Chile)
  !           Steven Peckham (NOAA/ESRL/GSD - Boulder, CO)
  !
  ! Main program of WRF/Chem global boundary condition generation code.  Responsible for 
  ! starting up the model, reading in configuration data, defining and initializing
  ! the top-level domain, either from initial or restart data, setting up time-keeping, and
  ! then calling the proper routine to get the global boundary condition data.
  ! After the boundary condition data construction is completed and the data is added to the
  ! wrfbdy_d01 data file, the model is properly shut down.
  !
  ! This code has been liberally adapted at Penn State to support the 'nesting' of 
  ! a WRF domain within a global model (to date: CarbonTracker and GEOS-Chem/CMS).
  ! Differences are noted below with a [PSU] annotate comment.
  !
  ! Other files needed:
  !  module_mozart_lib.f90  [PSU]: use module specific to global model
  !                                currently either module_GEOS_lib.f90 or module_CT_lib.f90 
  !                                BUT note that subroutines inside those modules will still refer to 'mozart'
  !  module_raqms_lib.f90   [PSU]: This module is not used: removed from calls and makefile for GEOS version
  !  module_wrfchem_lib.f90
  !  wrfchembc_namelist.input
  !  Makefile
  !
  ! The namelist input data file for the program is structured like the following
  ! &control
  ! 
  ! dir_wrf = '/data/wrfchem_data/wrfchem_raqmsbc/'   ! Run directory, wrf data directory
  ! fnb_wrf  = 'wrfbdy_d01'                           ! WRF boundary condition data file
  ! fni_wrf  = 'wrfinput_d01'                         ! WRF initial condition data file

  ! chem_bc_opt = 2                                   ! Global Model data used: 1=MOZART, 2=RAQMS
  !                                                   ! [PSU]: 1=global model, 2=fixed constant, 3=zero
  ! dir_global = '/data/RAQMS_DATA/'                  ! Global model data directory
  ! fn_global  = 'uwhyb_09_02_2006_00Z.chem.assim.nc' ! Global model data file name
  ! 
  ! nspec = 10                                        ! number of species listed by name in specname
  !                                                   ! [PSU]: nspec = 1; nspec must be 1 in this version
  ! 
  ! specname = 'o3',                                  ! name of chemical species to be added to boundary conditions
  !            'co',                                  ! [PSU]: specname = 'tracer_n' where n is the tracer number
  !            'no',                                  !        tracer name must match what is in wrfbdy file
  !            'no2',
  !            'n2o5',
  !            'h2o2',
  !            'ho',
  !            'ho2',
  !            'eth',
  !            'iso'
  ! /
  !
  ! To compile, type "make". The make file is set-up to compile the code for a 
  ! single processor using the Intel compiler (We assume you are using a 
  ! linux cluster).  The Makefile will need to be modified to work with other 
  ! computer systems.
  !
  ! To execute the program, type 
  !
  !  wrfchembc < wrfchembc_namelist.input
  !
  !</DESCRIPTION>

  use netcdf
  use module_wrfchem_lib
  use module_GEOS_lib     ! [PSU]: This version supports the GEOS-Chem/CMS global model
  !use module_raqms_lib   ! [PSU]: removed references to module_raqms_lib 2014/09/29


  implicit none

  !     parameters
  integer, parameter :: maxsize = 100
  integer, parameter :: kin = 15
  !     control variables
  character(len=maxsize) :: dir_global
  character(len=maxsize) :: fnb_wrf, fni_wrf, dir_wrf
  integer :: nspec, chem_bc_opt
  character(len=maxsize),DIMENSION(50) :: specname  ![PSU]: should dimension be removed?

  namelist /control/ dir_global, fnb_wrf, fni_wrf, dir_wrf, &
       chem_bc_opt, fn_global, nspec, specname

  !     wrf variables 
  !integer :: iptop_wrf   ![PSU]: WRF variable P_TOP is already defined as real/float
  real    :: ptop_wrf
  real, allocatable :: znu(:)
  real, allocatable :: xlon(:,:), xlat(:,:)

  !     global file name
  character(len=maxsize) :: fn_global

  !     species vmr
  real, allocatable :: bcxs(:,:,:)
  real, allocatable :: bcxe(:,:,:)
  real, allocatable :: bcys(:,:,:)
  real, allocatable :: bcye(:,:,:)
  !     species vmr tendencies
  real, allocatable :: btxs(:,:,:)
  real, allocatable :: btxe(:,:,:)
  real, allocatable :: btys(:,:,:)
  real, allocatable :: btye(:,:,:)
  !     species vmr old
  real, allocatable :: bcxso(:,:,:)
  real, allocatable :: bcxeo(:,:,:)
  real, allocatable :: bcyso(:,:,:)
  real, allocatable :: bcyeo(:,:,:)
  !     species vmr temp1
  real, allocatable :: t1xs(:,:,:,:)
  real, allocatable :: t1xe(:,:,:,:)
  real, allocatable :: t1ys(:,:,:,:)
  real, allocatable :: t1ye(:,:,:,:)
  !     species vmr temp2
  real, allocatable :: t2xs(:,:,:,:)
  real, allocatable :: t2xe(:,:,:,:)
  real, allocatable :: t2ys(:,:,:,:)
  real, allocatable :: t2ye(:,:,:,:)

  !     other working variables
  ! [PSU]:    fng added for 'global' file name (fni was reused for this) 2014/09/29
  !           dt_r added as integer variable (was defined in raqms module now removed) 2014/09/29
  !           Defined here even though the code it is in really does not matter
  !           it_CT is used to find the comparable time step in the global model data file
  !           i_con is used in the called model interpolation subroutine, and should not be here
  !           What is ios?
  integer            :: ns, ios, it, it_CT, dt_r  !,i_con
  character(len=maxsize) :: fnb, fni, fng
  character(len=maxsize)      :: spec1
  logical            :: time                       ![PSU]: not needed?
  real               :: t_ratio, it_ratio, f_int   ![PSU]: not needed?

  print *,"******************************************************"
  print *,"*   PROGRAM: WRF/CHEM GLOBAL BOUNDARIES              *"
  print *,"*              FOR WRF/CHEM VERSION 2.2              *"
  print *,"*                                                    *"
  print *,"*    PLEASE REPORT ANY BUGS TO WRF/CHEM HELP at      *"
  print *,"*                                                    *"
  print *,"*              wrfchemhelp.gsd@noaa.gov              *"
  print *,"*                                                    *"
  print *,"******************************************************"

  !     read control variables
  read( 5, nml=control )

  !     intialize module_wrfchem_lib
  fnb = trim(dir_wrf)//adjustl(fnb_wrf)
  fni = trim(dir_wrf)//adjustl(fni_wrf)

  !     print *,fnb
  print *,fn_global

  write(*,*) 'call init_wrfchem_lib: specname, nspec ', specname, nspec
  call init_wrfchem_lib( fnb, fni, specname, nspec, ntime )

  !     read top ref. pressure from wrf input data
  !  [PSU]: call changed to access model top pressure directly
  !         and NOT adjust the units from native [Pa]  2015-03-17
  call wrfchem_readscalar( 'P_TOP', ptop_wrf )
  !call wrfchem_readinteger( 'P_TOP', iptop_wrf )
  !ptop_wrf = float( iptop_wrf )*.01        ! convert to mb 
  print*, 'read  top ref. pressure from wrf input data'

  !     read eta values on half (mass) levels from wrfchem
  allocate( znu(nz) )
  call wrfchem_read2d( 'ZNU', znu )
  print*, 'read eta values on half (mass) levels',znu(1)

  !     read longitudes and latitudes from wrfchem input data
  allocate( xlon(nx,ny), xlat(nx,ny) )
  call wrfchem_read3d( 'XLONG', xlon )
  call wrfchem_read3d( 'XLAT',  xlat )

!  where( xlon < 0.0 )
!     xlon = xlon + 360.0
!  endwhere
!  print*, xlon,xlat   !commented out excessive displays 2014/09/29
  print*,'read longitudes and latitudes from wrfchem input data'

  !     initialize global data module
  ! [PSU]: fng used here instead of reuse of variable name fni  2014/09/29
  fng = trim(dir_global)//adjustl(fn_global)
  if(chem_bc_opt == 1) then
     call init_GEOS_lib( fng, ptop_wrf, xlon, xlat, znu, nx, ny, nz )
     print*,'Read CO2 boundary data from GEOS-CHEM'
  endif
  

  !     read global boundary data, interpolate and then write as BC
  allocate( bcxs(ny,nz,nw) )
  allocate( bcxe(ny,nz,nw) )
  allocate( bcys(nx,nz,nw) )
  allocate( bcye(nx,nz,nw) )
  allocate( bcxso(ny,nz,nw) )
  allocate( bcxeo(ny,nz,nw) )
  allocate( bcyso(nx,nz,nw) )
  allocate( bcyeo(nx,nz,nw) )
  allocate( btxs(ny,nz,nw) )
  allocate( btxe(ny,nz,nw) )
  allocate( btys(nx,nz,nw) )
  allocate( btye(nx,nz,nw) )
  allocate( t1xs(ny,nz,nw,nspec) )
  allocate( t1xe(ny,nz,nw,nspec) )
  allocate( t1ys(nx,nz,nw,nspec) )
  allocate( t1ye(nx,nz,nw,nspec) )
  allocate( t2xs(ny,nz,nw,nspec) )
  allocate( t2xe(ny,nz,nw,nspec) )
  allocate( t2ys(nx,nz,nw,nspec) )
  allocate( t2ye(nx,nz,nw,nspec) )

  !---------------------------------------------------------------------------------------------------
  !     time check of global data and WRF data
  !---------------------------------------------------------------------------------------------------

  ! [PSU]: review the following time check logic, not needed?
  !        main processing could begin at the do it=1, ntime loop
  !        For GEOS-Chem/CMS global model, dt=6 is set unconditionally in module_wrfchem_lib,
  !        and dt_m=3 is set unconditionally in module_GEOS_lib,
  !        and dt_r is not used,
  !        and time interpolation has been removed from this module (f_int, etc. commented out below).
  time = .true.

  it_ratio = 1.
  if ( time .and. chem_bc_opt == 1) then
     t_ratio  = float(dt)/float(dt_m)
     print*, 'Times are ok.'
     if ( t_ratio /= 1. ) then
        it_ratio = float(dt_m)/float(dt)
        t_ratio  = float(dt)/float(dt_m)
     endif
  else if ( time .and. chem_bc_opt == 2) then
     dt_r = 3   !arbitrary value assigned to variable not initialized 2014/09/29
     t_ratio  = float(dt)/float(dt_r)
     print *, 'Times are ok.' 
     if ( t_ratio /= 1. ) then
        it_ratio = float(dt_r)/float(dt)
        t_ratio  = float(dt)/float(dt_r)
     endif
  end if
  !write(*,*) 'We have to interpolate in time with factor: ', t_ratio, it_ratio  ![PSU]: not used

  !i_con = 2     ![PSU]: note this unconditional setting: This is not the same i_con as in the GEOS module!
                ![PSU]: This is confusing: variables it and it_CT here are passed to the GEOS module
                ![PSU]: and received there as i_con and i_carbtrac; This statement is therefore commented out.
  !f_int=.5      ![PSU]: not used
  do it = 1, ntime
     write(*,*) 'Time step number : ',it
     write(*,*) 'Time interval : ',dt
     ![PSU]: Correction to setting of it_CT for it==1 made September 2015
     ![PSU]: At same time hour_start of 12 UTC is added
     ![PSU]: Variable 'it' represents the time increment in wrfbdy; 'it_CT' represents the time increment 
     ![PSU]:   in GEOS-Chem file; these become i_con and i_carbtrac
     if (it.eq.1) then
        it_CT=(day_start - 1)*8 + 1
     else
        it_CT=(day_start - 1)*8+(it-1)*2
     endif
     if (hour_start == 12) then
        it_CT = it_CT + 4
     endif
     !write(*,*) 'Increment de lecture de CT/PCTM : ',it_CT
     write(*,*) 'Time step of GEOS/CMS CO2 at call to interpolate subroutine : ',it_CT
     do ns = 1, nspec
        ! [PSU]: intent here is to save the last time step for use in interpolation between,
        ! the t1.. variables are initialized in called module (but not before this point for it=1)
        ! for chem_bc_opt == 1
        ! and the t2...variable use is commented out below
        t2xs(:,:,:,ns) = t1xs(:,:,:,ns)
        t2xe(:,:,:,ns) = t1xe(:,:,:,ns)
        t2ys(:,:,:,ns) = t1ys(:,:,:,ns)
        t2ye(:,:,:,ns) = t1ye(:,:,:,ns)
        write(*,*) 'Name of the specie : ',specname(ns)  
      
        if(chem_bc_opt == 1) then 
           call mozart_interpolate4d( specname(ns), t1xs(:,:,:,ns), t1xe(:,:,:,ns), t1ys(:,:,:,ns), &
                t1ye(:,:,:,ns), nx, ny, nz, nw, it, it_CT )
        endif
! [PSU]: changed from 100. to 200. 2014/09/25 and to 300. in September 2015         
        if(chem_bc_opt == 2) then 
           t1xs(:,:,:,:)=300.
           t1xe(:,:,:,:)=300.
           t1ys(:,:,:,:)=300.
           t1ye(:,:,:,:)=300.
        endif
! [PSU]: use option 3 only for positive definite flux components      
        if(chem_bc_opt == 3) then 
           t1xs(:,:,:,:)=0.
           t1xe(:,:,:,:)=0.
           t1ys(:,:,:,:)=0.
           t1ye(:,:,:,:)=0.
        endif
      
        if (it >=2) then
           bcxso = bcxs
           bcxeo = bcxe
           bcyso = bcys
           bcyeo = bcye
        end if

        ! BC
        bcxs = t1xs(:,:,:,ns)!f_int * t2xs(:,:,:,ns) + ( 1. - f_int ) * t1xs(:,:,:,ns)
        bcxe = t1xe(:,:,:,ns)!f_int * t2xe(:,:,:,ns) + ( 1. - f_int ) * t1xe(:,:,:,ns)
        bcys = t1ys(:,:,:,ns)!f_int * t2ys(:,:,:,ns) + ( 1. - f_int ) * t1ys(:,:,:,ns)
        bcye = t1ye(:,:,:,ns)!f_int * t2ye(:,:,:,ns) + ( 1. - f_int ) * t1ye(:,:,:,ns)

        ! [PSU] tendencies for options 2 and 3 updated November 2015
        !       These were set at 1e-5 which is not nearly small enough,
        !       and introduced 0.036 ppm 'leaks' per tracer in waves on inflows on the walls,
        !       eventually amounting to approximately 0.1 ppm per tracer (option 2 and 3) in the interior.
        !       Tendencies changed to 0.
        if (it >=2) then
           ! Tendencies
           if(chem_bc_opt == 1) then 
              btxs = ( bcxs - bcxso )/(dt*3600.)
              btxe = ( bcxe - bcxeo )/(dt*3600.)
              btys = ( bcys - bcyso )/(dt*3600.)
              btye = ( bcye - bcyeo )/(dt*3600.)
           endif
           if(chem_bc_opt == 2) then 
              btxs = 0.
              btxe = 0.
              btys = 0.
              btye = 0.
           endif
           if(chem_bc_opt == 3) then 
              btxs = 0.
              btxe = 0.
              btys = 0.
              btye = 0.
           endif
        endif
     
        !print *,'f_int ',1. - f_int, f_int   ![PSU]: not used, so why print it?
      
        write(*,*) 'Write the variable in the netcdf file'
        call wrfchem_write4d( specname(ns), bcxs, bcxe, bcys, bcye, btxs, btxe, btys, btye, it )
     enddo
   
  end do

  !     exit from libs
  ! [PSU]: call to exit_raqms_lib removed 2014/09/29 - and module removed from Makefile
  call exit_wrfchem_lib( )
  if(chem_bc_opt == 1) then
     call exit_mozart_lib( )
  !elseif(chem_bc_opt == 2) then
  !   call exit_raqms_lib( )
  endif

  !     deallocate memory space
  deallocate( znu, xlon, xlat, bcxs, bcxe, bcys, bcye,    &
       bcxso, bcxeo, bcyso, bcyeo, btxs, btxe, btys, btye,&
       t1xs, t1xe, t1ys, t1ye, t2xs, t2xe, t2ys, t2ye  )

  print*,'bc_wrfchem completed successfully'
end program main_bc_wrfchem
