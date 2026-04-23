
module read_climate_nudge_data_nlm_mod

use fms_mod, only: check_nml_error, &
                   stdlog, mpp_pe, mpp_root_pe, write_version_number, &
                   string, error_mesg, FATAL, NOTE
use mpp_mod, only: input_nml_file
use fms2_io_mod,   only: FmsNetcdfFile_t, open_file, close_file, read_data, &
                        get_num_dimensions, get_dimension_names, get_dimension_size, &
                        get_num_variables, get_variable_names, get_variable_num_dimensions, &
                        get_variable_size, get_variable_attribute, variable_exists
use constants_mod, only: PI, GRAV, RDGAS, RVGAS
use fv_arrays_nlm_mod,  only: REAL4, REAL8, FVPRC

implicit none
private

public :: read_climate_nudge_data_init, read_time, read_grid,  &
          read_climate_nudge_data, read_climate_nudge_data_end
public :: read_sub_domain_init

interface read_climate_nudge_data
   module procedure read_climate_nudge_data_2d
   module procedure read_climate_nudge_data_3d
end interface

  real(FVPRC), parameter :: P0 = 1.e5
  real(FVPRC), parameter :: D608 = RVGAS/RDGAS - 1.

  integer, parameter :: NUM_REQ_AXES = 3
  integer, parameter :: INDEX_LON = 1, INDEX_LAT = 2, INDEX_LEV = 3
  character(len=8), dimension(NUM_REQ_AXES) :: required_axis_names = &
                                     (/ 'lon', 'lat', 'lev' /)

  integer, parameter :: NUM_REQ_FLDS = 9
  integer, parameter :: INDEX_P0 = 1, INDEX_AK = 2, INDEX_BK = 3, &
                        INDEX_ZS = 4, INDEX_PS = 5,               &
                        INDEX_T  = 6, INDEX_Q  = 7,               &
                        INDEX_U  = 8, INDEX_V  = 9
  character(len=8), dimension(NUM_REQ_FLDS) :: required_field_names = &
       (/ 'P0  ', 'hyai', 'hybi', 'PHI ', 'PS  ', 'T   ', 'Q   ', 'U   ', 'V   ' /)
 
  integer, parameter :: MAXFILES = 53
  character(len=256) :: filenames(MAXFILES)
  character(len=256) :: filename_tails(MAXFILES)
  character(len=256) :: filename_head
  integer :: read_buffer_size
  integer :: nfiles = 0
  logical :: module_is_initialized = .false.

  namelist /read_climate_nudge_data_nml/ filename_tails, read_buffer_size, &
                                         filename_head

! dimensions for checking
  integer :: global_axis_size(NUM_REQ_AXES), numtime, sub_domain_latitude_size
  integer, allocatable :: file_index(:)

type filedata_type
  type(FmsNetcdfFile_t) :: fileobj
  integer, allocatable :: length_axes(:)
  integer :: ndim, nvar, ntim
  integer :: time_offset
  integer, dimension(NUM_REQ_FLDS) :: field_index
  integer, dimension(NUM_REQ_AXES) :: axis_index
  character(len=64), dimension(NUM_REQ_AXES) :: axis_dim_names
  character(len=64), dimension(NUM_REQ_FLDS) :: field_var_names
end type

  type(filedata_type), allocatable :: Files(:)

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine read_climate_nudge_data_init (nlon, nlat, nlev, ntime)
integer, intent(out) :: nlon, nlat, nlev, ntime
! returns dimension lengths of input data set
!   nlon, nlat  lat/lon grid size
!   nlev        number of levels
!   ntime       number of time levels

  integer :: iunit, ierr, io
  character(len=128) :: name
  integer :: istat, i, j, k, n, nd, siz(4), i1, i2
  character(len=64), allocatable :: dim_names(:)
  character(len=64), allocatable :: var_names(:)
  integer :: dim_size

  if (module_is_initialized) return
  ! initial file names to blanks
  do n = 1, MAXFILES
     do i = 1, len(filename_tails(n))
        filename_tails(n)(i:i) = ' '
        filenames(n)(i:i) = ' '
     enddo
  enddo

!----- read namelist -----
  read (input_nml_file, nml=read_climate_nudge_data_nml, iostat=io)
  ierr = check_nml_error (io, 'read_climate_nudge_data_nml')

!----- write version and namelist to log file -----

  iunit = stdlog()
  call write_version_number ( '0.0', 'fv3-jedi-lm' )
  if (mpp_pe() == mpp_root_pe()) write (iunit, nml=read_climate_nudge_data_nml)

  ! determine the number of files
  do n = 1, MAXFILES
     if (filename_tails(n)(1:1) .eq. ' ') exit
     nfiles = n
  enddo
  do n=1,nfiles
    filenames(n) = trim(filename_head)//trim(filename_tails(n))
  end do

  allocate(Files(nfiles))
  numtime = 0

! open input file(s)
  do n = 1, nfiles
      if (.not. open_file(Files(n)%fileobj, trim(filenames(n)), 'read')) then
         call error_mesg ('read_climate_nudge_data_nlm_mod', 'cannot open file '//trim(filenames(n)), FATAL)
      endif

      Files(n)%ndim = get_num_dimensions(Files(n)%fileobj)
      Files(n)%nvar = get_num_variables(Files(n)%fileobj)
      call get_dimension_size(Files(n)%fileobj, 'time', Files(n)%ntim)

      allocate (Files(n)%length_axes(Files(n)%ndim))
      allocate (dim_names(Files(n)%ndim))
      call get_dimension_names(Files(n)%fileobj, dim_names)

      ! inquire dimension sizes
      do i = 1, Files(n)%ndim
         call get_dimension_size(Files(n)%fileobj, dim_names(i), Files(n)%length_axes(i))
         name = dim_names(i)
         do j = 1, NUM_REQ_AXES
            if (trim(name) .eq. trim(required_axis_names(j))) then
               call check_axis_size (j,Files(n)%length_axes(i))
               Files(n)%axis_dim_names(j) = dim_names(i)
               Files(n)%axis_index(j) = i
               exit
            endif
         enddo
      enddo
      deallocate(dim_names)
      ! time axis indexing
      Files(n)%time_offset = numtime
      numtime = numtime + Files(n)%ntim

      allocate(var_names(Files(n)%nvar))
      call get_variable_names(Files(n)%fileobj, var_names)
      Files(n)%field_index = 0
      do i = 1, Files(n)%nvar
         name = var_names(i)
         nd = get_variable_num_dimensions(Files(n)%fileobj, trim(name))
         if (nd > 0) then
            call get_variable_size(Files(n)%fileobj, trim(name), siz(1:nd))
         endif
         do j = 1, NUM_REQ_FLDS
            if (trim(name) .eq. trim(required_field_names(j))) then
               Files(n)%field_index(j) = i
               Files(n)%field_var_names(j) = name
               if (j .gt. 3) then
                  call check_resolution (siz(1:nd))
               endif
               exit
            endif
         enddo
      enddo
      deallocate(var_names)

  enddo ! "n" files loop

  ! setup file indexing
  allocate(file_index(numtime))
  i2 = 0
  do n = 1, nfiles
     i1 = i2+1
     i2 = i2+Files(n)%ntim
     file_index(i1:i2) = n
  enddo

      sub_domain_latitude_size = global_axis_size(INDEX_LAT)

    ! output arguments
      nlon = global_axis_size(INDEX_LON)
      nlat = global_axis_size(INDEX_LAT)
      nlev = global_axis_size(INDEX_LEV)
      ntime = numtime

      module_is_initialized = .true.

end subroutine read_climate_nudge_data_init

!###############################################################################

subroutine read_time ( times, units, calendar )
real(FVPRC),          intent(out) :: times(:)
character(len=*), intent(out) :: units, calendar
integer :: istat, i1, i2, n

   if (.not.module_is_initialized) then
     call error_mesg ('read_climate_nudge_data_nlm_mod/read_time',  &
                                        'module not initialized', FATAL)
   endif

   if (size(times(:)) < numtime) then
      call error_mesg ('read_climate_nudge_data_nlm_mod', 'argument times too small in read_time', FATAL)
   endif

 ! data
   i2 = 0
   do n = 1, nfiles
      i1 = i2+1
      i2 = i2+Files(n)%ntim
      if( n == 1) then
         if (variable_exists(Files(n)%fileobj, 'time')) then
            call get_variable_attribute(Files(n)%fileobj, 'time', 'units', units)
            call get_variable_attribute(Files(n)%fileobj, 'time', 'calendar', calendar)
         else
            units = 'days since 0001-01-01 00:00:00'
            calendar = 'gregorian'
         endif
      endif
      call read_data(Files(n)%fileobj, 'time', times(i1:i2))
   enddo

! NOTE: need to do the conversion to time_type in this routine
!       this will allow different units and calendars for each file

end subroutine read_time

!###############################################################################

subroutine read_grid ( lon, lat, ak, bk )
real(FVPRC), intent(out), dimension(:) :: lon, lat, ak, bk

 real(FVPRC) :: pref
 integer :: istat

   if (.not.module_is_initialized) then
     call error_mesg ('read_climate_nudge_data_nlm_mod/read_grid',  &
                                        'module not initialized', FATAL)
   endif


    ! static fields from first file only
      call read_data(Files(1)%fileobj, Files(1)%axis_dim_names(INDEX_LON), lon)
      call read_data(Files(1)%fileobj, Files(1)%axis_dim_names(INDEX_LAT), lat)

    ! units are assumed to be degrees east and north
    ! convert to radians
      lon = lon * PI/180.
      lat = lat * PI/180.

    ! vertical coodinate
      if (Files(1)%field_index(INDEX_AK) .gt. 0) then
         call read_data(Files(1)%fileobj, Files(1)%field_var_names(INDEX_AK), ak)
         if (Files(1)%field_index(INDEX_P0) .gt. 0) then
            call read_data(Files(1)%fileobj, Files(1)%field_var_names(INDEX_P0), pref)
         else
            pref = P0
         endif
         ak = ak*pref
      else
         ak = 0.
      endif
 
      call read_data(Files(1)%fileobj, Files(1)%field_var_names(INDEX_BK), bk)


end subroutine read_grid

!###############################################################################

subroutine read_sub_domain_init ( ylo, yhi, ydat, js, je )
 real(FVPRC),    intent(in)  :: ylo, yhi, ydat(:)
 integer, intent(out) :: js, je
 integer :: j

   if (.not.module_is_initialized) then
     call error_mesg ('read_climate_nudge_data_nlm_mod/read_sub_domain_init',  &
                                        'module not initialized', FATAL)
   endif
   ! increasing data
   if (ydat(1) < ydat(2)) then
      js = 1
      do j = 1, size(ydat(:))-1
         if (ylo >= ydat(j) .and. ylo <= ydat(j+1)) then
            js = j
            exit
         endif
      enddo

      if (ylo < -1.5) then
         print *, 'ylo=',ylo
         print *, 'js,ydat=',js,ydat(js)
         print *, 'ydat=',ydat(:js+2)
      endif

      je = size(ydat(:))
      do j = js, size(ydat(:))-1
         if (yhi >= ydat(j) .and. yhi <= ydat(j+1)) then
            je = j+1
            exit
         endif
      enddo

      if (yhi > 1.5) then
         print *, 'yhi=',yhi
         print *, 'je,ydat=',je,ydat(je)
         print *, 'ydat=',ydat(je-2:)
      endif

   ! decreasing data (may not work)
   else
      call error_mesg ('read_climate_nudge_data_nlm_mod', 'latitude values for observational data decrease with increasing index', NOTE)
      je = size(ydat(:))-1
      do j = 1, size(ydat(:))-1
         if (ylo >= ydat(j+1) .and. ylo <= ydat(j)) then
            je = j+1
            exit
         endif
      enddo

      js = 1
      do j = 1, je
         if (yhi >= ydat(j+1) .and. yhi <= ydat(j)) then
            js = j
            exit
         endif
      enddo

   endif

   sub_domain_latitude_size = je-js+1

 end subroutine read_sub_domain_init

!###############################################################################

subroutine read_climate_nudge_data_2d (itime, field, dat, is, js)
integer,          intent(in) :: itime
character(len=4), intent(in) :: field
real(FVPRC),             intent(out), dimension(:,:) :: dat
integer,          intent(in),  optional       :: is, js
integer :: istat, atime, n, this_index
integer :: nread(4), start(4)

   if (.not.module_is_initialized) then
     call error_mesg ('read_climate_nudge_data_nlm_mod',  &
                                        'module not initialized', FATAL)
   endif
     ! time index check
      if (itime < 1 .or. itime > numtime) then
         call error_mesg ('read_climate_nudge_data_nlm_mod', 'itime out of range', FATAL)
      endif

     ! check dimensions 
     if (present(js)) then
        if (size(dat,1) .ne. global_axis_size(INDEX_LON) .or. &
            size(dat,2) .ne. sub_domain_latitude_size) then
            !write (*,'(a)') 'climate_nudge_data_mod: size dat2d = '//trim(string(size(dat,1)))//' x '//trim(string(size(dat,2)))// &
            !         '  <-vs->  '//trim(string(global_axis_size(INDEX_LON)))//' x '//trim(string(sub_domain_latitude_size))
            call error_mesg ('read_climate_nudge_data_nlm_mod', 'incorrect 2d array dimensions', FATAL)
        endif
     else
        if (size(dat,1) .ne. global_axis_size(INDEX_LON) .or. &
            size(dat,2) .ne. global_axis_size(INDEX_LAT))     &
            call error_mesg ('read_climate_nudge_data_nlm_mod', 'incorrect 2d array dimensions', FATAL)
     endif

     ! check field
     if (field .eq. 'phis') then
        this_index = INDEX_ZS
     else if (field .eq. 'psrf') then
        this_index = INDEX_PS
     else
         call error_mesg ('read_climate_nudge_data_nlm_mod', 'incorrect field requested in read_climate_nudge_data_2d', FATAL)
     endif
     
     ! file index and actual time index in file
     n = file_index(itime)
     atime = itime - Files(n)%time_offset

     start = 1
     if (present(is)) start(1) = is
     if (present(js)) start(2) = js
     start(3) = atime

     nread = 1
     nread(1) = size(dat,1)
     nread(2) = size(dat,2)
     
     call read_data(Files(n)%fileobj, Files(n)%field_var_names(this_index), dat, corner=start, edge_lengths=nread)
  
      ! geopotential height (convert to m2/s2 if necessary)
     if (field .eq. 'phis') then
        if (maxval(dat) > 1000.*GRAV) then
          ! do nothing
        else
          dat = dat * GRAV
        endif
     endif

end subroutine read_climate_nudge_data_2d

!###############################################################################

subroutine read_climate_nudge_data_3d (itime, field, dat, is, js)
integer,          intent(in) :: itime
character(len=4), intent(in) :: field
real(FVPRC),             intent(out), dimension(:,:,:) :: dat
integer,          intent(in),  optional         :: is, js
integer :: istat, atime, n, this_index, start(4), nread(4)
!logical :: convert_virt_temp = .false.

   if (.not.module_is_initialized) then
     call error_mesg ('read_climate_nudge_data_nlm_mod',  &
                                        'module not initialized', FATAL)
   endif

     ! time index check
     if (itime < 1 .or. itime > numtime) then
        call error_mesg ('read_climate_nudge_data_nlm_mod', 'itime out of range', FATAL)
     endif

     ! check dimensions
     if (present(js)) then
        if (size(dat,1) .ne. global_axis_size(INDEX_LON) .or. &
            size(dat,2) .ne. sub_domain_latitude_size    .or. &
            size(dat,3) .ne. global_axis_size(INDEX_LEV)) then
            !write (*,'(a)') 'climate_nudge_data_mod: size dat3d = '//trim(string(size(dat,1)))//' x '//trim(string(size(dat,2)))// &
            !                                        ' x '//trim(string(size(dat,3)))
            call error_mesg ('read_climate_nudge_data_nlm_mod', 'incorrect 3d array dimensions', FATAL)
        endif
     else
        if (size(dat,1) .ne. global_axis_size(INDEX_LON) .or. &
            size(dat,2) .ne. global_axis_size(INDEX_LAT) .or. &
            size(dat,3) .ne. global_axis_size(INDEX_LEV))     &
            call error_mesg ('read_climate_nudge_mod', 'incorrect 3d array dimensions', FATAL)
     endif

     ! check field
     if (field .eq. 'temp') then
        this_index = INDEX_T
     else if (field .eq. 'qhum') then
        this_index = INDEX_Q
     else if (field .eq. 'uwnd') then
        this_index = INDEX_U
     else if (field .eq. 'vwnd') then
        this_index = INDEX_V
     else
        call error_mesg ('read_climate_nudge_data_nlm_mod', 'incorrect field requested in read_climate_nudge_data_3d', FATAL)
     endif
     

     ! file index and actual time index in file
     n = file_index(itime)
     atime = itime - Files(n)%time_offset

     start = 1
     if (present(is)) start(1) = is
     if (present(js)) start(2) = js
     start(4) = atime

     nread = 1
     nread(1) = size(dat,1)
     nread(2) = size(dat,2)
     nread(3) = size(dat,3)

     call read_data(Files(n)%fileobj, Files(n)%field_var_names(this_index), dat, corner=start, edge_lengths=nread)

     ! convert virtual temp to temp
     ! necessary for some of the high resol AVN analyses
     !if (convert_virt_temp) then
     !   temp = temp/(1.+D608*qhum)
     !endif

end subroutine read_climate_nudge_data_3d

!###############################################################################

subroutine read_climate_nudge_data_end
integer :: istat, n

  if ( .not.module_is_initialized) return
  do n = 1, nfiles
     call close_file(Files(n)%fileobj)
  enddo
  deallocate (Files)
  module_is_initialized = .false.

end subroutine read_climate_nudge_data_end

!###############################################################################

 subroutine check_axis_size (ind,lendim)
 integer, intent(in) :: ind,lendim

   ! once the axis size is set all subsuquent axes must be the same
   if (global_axis_size(ind) .gt. 0) then
      if (global_axis_size(ind) .ne. lendim) then
         call error_mesg ('read_climate_nudge_data_nlm_mod', 'incorrect axis size for axis = '//trim(required_axis_names(ind)), FATAL)
      endif
   else
      global_axis_size(ind) = lendim
   endif

 end subroutine check_axis_size

!------------------------------------

 subroutine check_resolution (axis_len)
 integer, intent(in) :: axis_len(:)

   if (size(axis_len(:)) .lt. 2) then
      call error_mesg ('read_climate_nudge_data_nlm_mod', 'incorrect number of array dimensions', FATAL)
   endif
   if (axis_len(1) .ne. global_axis_size(INDEX_LON)) then
      call error_mesg ('read_climate_nudge_data_nlm_mod', 'incorrect array dimension one', FATAL)
   endif
   if (axis_len(2) .ne. global_axis_size(INDEX_LAT)) then
      call error_mesg ('read_climate_nudge_data_nlm_mod', 'incorrect array dimension two', FATAL)
   endif
   if (size(axis_len(:)) .gt. 3) then
      if (axis_len(3) .ne. global_axis_size(INDEX_LEV)) then
         call error_mesg ('read_climate_nudge_data_nlm_mod', 'incorrect array dimension three', FATAL)
      endif
   endif

 end subroutine check_resolution

!###############################################################################

end module read_climate_nudge_data_nlm_mod

