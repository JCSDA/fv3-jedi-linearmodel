!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of fvGFS.                                       *
!*                                                                     *
!* fvGFS is free software; you can redistribute it and/or modify it    *
!* and are expected to follow the terms of the GNU General Public      *
!* License as published by the Free Software Foundation; either        *
!* version 2 of the License, or (at your option) any later version.    *
!*                                                                     *
!* fvGFS is distributed in the hope that it will be useful, but        *
!* WITHOUT ANY WARRANTY; without even the implied warranty of          *
!* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU   *
!* General Public License for more details.                            *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************
!!!NOTE: Merging in the seasonal forecast initialization code
!!!!     has proven problematic in the past, since many conflicts
!!!!     occur. Leaving this for now --- lmh 10aug15

module fv_io_nlm_mod

  !<OVERVIEW>
  ! Restart facilities for FV core
  !</OVERVIEW>
  !<DESCRIPTION>
  ! This module writes and reads restart files for the FV core. Additionally
  ! it provides setup and calls routines necessary to provide a complete restart
  ! for the model.
  !</DESCRIPTION>

  use fms2_io_mod,              only: FmsNetcdfFile_t, FmsNetcdfDomainFile_t, &
                                     open_file, close_file, &
                                     register_restart_field, read_restart, write_restart, &
                                     file_exists, variable_exists, set_filename_appendix
  use fms_mod,                 only: set_domain, nullify_domain
  use mpp_mod,                 only: mpp_error, FATAL, NOTE, WARNING, mpp_root_pe, &
                                     mpp_sync, mpp_pe, mpp_declare_pelist
  use mpp_domains_mod,         only: domain2d, EAST, WEST, NORTH, CENTER, SOUTH, CORNER, &
                                     mpp_get_compute_domain, mpp_get_data_domain, & 
                                     mpp_get_layout, mpp_get_ntile_count, &
                                     mpp_get_global_domain
  use tracer_manager_mod,      only: tr_get_tracer_names=>get_tracer_names, &
                                     get_tracer_names, get_number_tracers, &
                                     set_tracer_profile, &
                                     get_tracer_index
  use field_manager_mod,       only: MODEL_ATMOS  
  use external_sst_nlm_mod,        only: sst_ncep, sst_anom, use_ncep_sst
  use fv_arrays_nlm_mod,           only: fv_atmos_type, fv_nest_BC_type_3D
  use fv_eta_nlm_mod,              only: set_eta

  use fv_mp_nlm_mod,               only: ng, mp_gather, is_master
  implicit none
  private

  public :: fv_io_init, fv_io_exit, fv_io_read_restart, remap_restart, fv_io_write_restart
  public :: fv_io_read_tracers, fv_io_register_restart, fv_io_register_nudge_restart
  public :: fv_io_register_restart_BCs, fv_io_register_restart_BCs_NH
  public :: fv_io_write_BCs, fv_io_read_BCs

  logical                       :: module_is_initialized = .FALSE.

contains 

  !#####################################################################
  ! <SUBROUTINE NAME="fv_io_init">
  !
  ! <DESCRIPTION>
  ! Initialize the fv core restart facilities
  ! </DESCRIPTION>
  !
  subroutine fv_io_init()
    module_is_initialized = .TRUE.
  end subroutine fv_io_init
  ! </SUBROUTINE> NAME="fv_io_init"


  !#####################################################################
  ! <SUBROUTINE NAME="fv_io_exit">
  !
  ! <DESCRIPTION>
  ! Close the fv core restart facilities
  ! </DESCRIPTION>
  !
  subroutine fv_io_exit
    module_is_initialized = .FALSE.
  end subroutine fv_io_exit
  ! </SUBROUTINE> NAME="fv_io_exit"



  !#####################################################################
  ! <SUBROUTINE NAME="fv_io_read_restart">
  !
  ! <DESCRIPTION>
  ! Write the fv core restart quantities 
  ! </DESCRIPTION>
  subroutine  fv_io_read_restart(fv_domain,Atm)
    type(domain2d),      intent(inout) :: fv_domain
    type(fv_atmos_type), intent(inout) :: Atm(:)

    character(len=64)    :: fname, tracer_name
    character(len=6)  :: stile_name
    integer              :: isc, iec, jsc, jec, n, nt, nk, ntracers
    integer              :: ntileMe
    integer              :: ks, ntiles
    real                 :: ptop

    character(len=128)           :: tracer_longname, tracer_units

    ntileMe = size(Atm(:))  ! This will need mods for more than 1 tile per pe

    call read_restart(Atm(1)%Fv_restart)
    call close_file(Atm(1)%Fv_restart)

    if ( use_ncep_sst .or. Atm(1)%flagstruct%nudge .or. Atm(1)%flagstruct%ncep_ic ) then
       call mpp_error(NOTE, 'READING FROM SST_RESTART DISABLED')
       !call restore_state(Atm(1)%SST_restart)
    endif

! fix for single tile runs where you need fv_core.res.nc and fv_core.res.tile1.nc
    ntiles = mpp_get_ntile_count(fv_domain)
    if(ntiles == 1 .and. .not. Atm(1)%neststruct%nested) then
       stile_name = '.tile1'
    else
       stile_name = ''
    endif
 
    do n = 1, ntileMe
       call read_restart(Atm(n)%Fv_tile_restart)
       call close_file(Atm(n)%Fv_tile_restart)

!--- restore data for fv_tracer - if it exists
       fname = 'INPUT/fv_tracer.res'//trim(stile_name)//'.nc'
       if (file_exists(fname)) then
         call read_restart(Atm(n)%Tra_restart)
         call close_file(Atm(n)%Tra_restart)
       else
         call mpp_error(NOTE,'==> Warning from fv_read_restart: Expected file '//trim(fname)//' does not exist')
       endif

!--- restore data for surface winds - if it exists
       fname = 'INPUT/fv_srf_wnd.res'//trim(stile_name)//'.nc'
       if (file_exists(fname)) then
         call read_restart(Atm(n)%Rsf_restart)
         call close_file(Atm(n)%Rsf_restart)
         Atm(n)%flagstruct%srf_init = .true.
       else
         call mpp_error(NOTE,'==> Warning from fv_read_restart: Expected file '//trim(fname)//' does not exist')
         Atm(n)%flagstruct%srf_init = .false.
       endif

       if ( Atm(n)%flagstruct%fv_land ) then
!--- restore data for mg_drag - if it exists
         fname = 'INPUT/mg_drag.res'//trim(stile_name)//'.nc'
         if (file_exists(fname)) then
           call read_restart(Atm(n)%Mg_restart)
           call close_file(Atm(n)%Mg_restart)
         else
           call mpp_error(NOTE,'==> Warning from fv_read_restart: Expected file '//trim(fname)//' does not exist')
         endif
!--- restore data for fv_land - if it exists
         fname = 'INPUT/fv_land.res'//trim(stile_name)//'.nc'
         if (file_exists(fname)) then
           call read_restart(Atm(n)%Lnd_restart)
           call close_file(Atm(n)%Lnd_restart)
         else
           call mpp_error(NOTE,'==> Warning from fv_read_restart: Expected file '//trim(fname)//' does not exist')
         endif
       endif

    end do

    return

  end subroutine  fv_io_read_restart
  ! </SUBROUTINE> NAME="fv_io_read_restart"
  !#####################################################################


  subroutine fv_io_read_tracers(fv_domain,Atm)
    type(domain2d),      intent(inout) :: fv_domain
    type(fv_atmos_type), intent(inout) :: Atm(:)
    integer :: n, ntracers, ntprog, nt, isc, iec, jsc, jec
    character(len=6) :: stile_name
    character(len=64):: fname, tracer_name
    type(FmsNetcdfDomainFile_t) :: Tra_restart_r
    integer :: ntiles

    n = 1
    isc = Atm(n)%bd%isc
    iec = Atm(n)%bd%iec
    jsc = Atm(n)%bd%jsc
    jec = Atm(n)%bd%jec
    call get_number_tracers(MODEL_ATMOS, num_tracers=ntracers, num_prog=ntprog)

! fix for single tile runs where you need fv_core.res.nc and fv_core.res.tile1.nc
    ntiles = mpp_get_ntile_count(fv_domain)
    if(ntiles == 1 .and. .not. Atm(1)%neststruct%nested) then
       stile_name = '.tile1'
    else
       stile_name = ''
    endif

    fname = 'fv_tracer.res'//trim(stile_name)//'.nc'
    if (file_exists('INPUT/'//trim(fname))) then
      if (open_file(Tra_restart_r, fname, 'read', fv_domain, is_restart=.true.)) then
        do nt = 2, ntprog
           call get_tracer_names(MODEL_ATMOS, nt, tracer_name)
           call set_tracer_profile (MODEL_ATMOS, nt, Atm(n)%q(isc:iec,jsc:jec,:,nt)  )
           if (variable_exists(Tra_restart_r, tracer_name)) then
             call register_restart_field(Tra_restart_r, tracer_name, Atm(n)%q(:,:,:,nt))
           endif
        enddo
        do nt = ntprog+1, ntracers
           call get_tracer_names(MODEL_ATMOS, nt, tracer_name)
           call set_tracer_profile (MODEL_ATMOS, nt, Atm(n)%qdiag(isc:iec,jsc:jec,:,nt)  )
           if (variable_exists(Tra_restart_r, tracer_name)) then
             call register_restart_field(Tra_restart_r, tracer_name, Atm(n)%qdiag(:,:,:,nt))
           endif
        enddo
        call read_restart(Tra_restart_r)
        call close_file(Tra_restart_r)
      endif
    else
      call mpp_error(NOTE,'==> Warning from fv_io_read_tracers: Expected file '//trim(fname)//' does not exist')
    endif

    return

  end subroutine  fv_io_read_tracers


  subroutine  remap_restart(fv_domain,Atm)
  use fv_mapz_nlm_mod,       only: rst_remap

    type(domain2d),      intent(inout) :: fv_domain
    type(fv_atmos_type), intent(inout) :: Atm(:)

    character(len=64)    :: fname, tracer_name
    character(len=6)     :: stile_name
    integer              :: isc, iec, jsc, jec, n, nt, nk, ntracers, ntprog, ntdiag
    integer              :: isd, ied, jsd, jed
    integer              :: ntiles
    type(FmsNetcdfFile_t) :: Fv_restart_r
    type(FmsNetcdfDomainFile_t) :: FV_tile_restart_r, Tra_restart_r

!
!-------------------------------------------------------------------------
    real, allocatable:: ak_r(:), bk_r(:)
    real, allocatable:: u_r(:,:,:), v_r(:,:,:), pt_r(:,:,:), delp_r(:,:,:)
    real, allocatable:: w_r(:,:,:), delz_r(:,:,:), ze0_r(:,:,:)
    real, allocatable:: q_r(:,:,:,:), qdiag_r(:,:,:,:)
!-------------------------------------------------------------------------
    integer npz, npz_rst, ng

    npz     = Atm(1)%npz       ! run time z dimension
    npz_rst = Atm(1)%flagstruct%npz_rst   ! restart z dimension
    isc = Atm(1)%bd%isc; iec = Atm(1)%bd%iec; jsc = Atm(1)%bd%jsc; jec = Atm(1)%bd%jec
    ng = Atm(1)%ng

    isd = isc - ng;  ied = iec + ng
    jsd = jsc - ng;  jed = jec + ng


!   call get_number_tracers(MODEL_ATMOS, num_tracers=ntracers)
    ntprog = size(Atm(1)%q,4)  ! Temporary until we get tracer manager integrated
    ntdiag = size(Atm(1)%qdiag,4)
    ntracers = ntprog+ntdiag

!    ntileMe = size(Atm(:))  ! This will have to be modified for mult tiles per PE


! Allocate arrays for reading old restart file:
    allocate ( ak_r(npz_rst+1) )
    allocate ( bk_r(npz_rst+1) )

    allocate ( u_r(isc:iec,  jsc:jec+1,npz_rst) )
    allocate ( v_r(isc:iec+1,jsc:jec  ,npz_rst) )

    allocate (   pt_r(isc:iec, jsc:jec,  npz_rst) )
    allocate ( delp_r(isc:iec, jsc:jec,  npz_rst) )
    allocate (    q_r(isc:iec, jsc:jec,  npz_rst, ntprog) )
    allocate (qdiag_r(isc:iec, jsc:jec,  npz_rst, ntprog+1:ntracers) )

    if ( (.not.Atm(1)%flagstruct%hydrostatic) .and. (.not.Atm(1)%flagstruct%make_nh) ) then
           allocate (    w_r(isc:iec, jsc:jec,  npz_rst) )
           allocate ( delz_r(isc:iec, jsc:jec,  npz_rst) )
           if ( Atm(1)%flagstruct%hybrid_z )   &
           allocate ( ze0_r(isc:iec, jsc:jec,  npz_rst+1) )
    endif

    fname = 'fv_core.res.nc'
    if (open_file(Fv_restart_r, fname, 'read', is_restart=.true.)) then
      call register_restart_field(Fv_restart_r, 'ak', ak_r(:))
      call register_restart_field(Fv_restart_r, 'bk', bk_r(:))
      call read_restart(Fv_restart_r)
      call close_file(Fv_restart_r)
    endif

! fix for single tile runs where you need fv_core.res.nc and fv_core.res.tile1.nc
    ntiles = mpp_get_ntile_count(fv_domain)
    if(ntiles == 1 .and. .not. Atm(1)%neststruct%nested) then
       stile_name = '.tile1'
    else
       stile_name = ''
    endif

!    do n = 1, ntileMe
    n = 1
       fname = 'fv_core.res'//trim(stile_name)//'.nc'
       if (open_file(Fv_tile_restart_r, fname, 'read', fv_domain, is_restart=.true.)) then
         call register_restart_field(Fv_tile_restart_r, 'u', u_r)
         call register_restart_field(Fv_tile_restart_r, 'v', v_r)
         if (.not.Atm(n)%flagstruct%hydrostatic) then
            if (variable_exists(Fv_tile_restart_r, 'W')) &
              call register_restart_field(Fv_tile_restart_r, 'W', w_r)
            if (variable_exists(Fv_tile_restart_r, 'DZ')) &
              call register_restart_field(Fv_tile_restart_r, 'DZ', delz_r)
            if ( Atm(n)%flagstruct%hybrid_z ) then
               if (variable_exists(Fv_tile_restart_r, 'ZE0')) &
                 call register_restart_field(Fv_tile_restart_r, 'ZE0', ze0_r)
            endif
         endif
         call register_restart_field(Fv_tile_restart_r, 'T', pt_r)
         call register_restart_field(Fv_tile_restart_r, 'delp', delp_r)
         call register_restart_field(Fv_tile_restart_r, 'phis', Atm(n)%phis)
         call read_restart(FV_tile_restart_r)
         call close_file(FV_tile_restart_r)
       endif
       fname = 'INPUT/fv_srf_wnd.res'//trim(stile_name)//'.nc'
       if (file_exists(fname)) then
         call read_restart(Atm(n)%Rsf_restart)
         call close_file(Atm(n)%Rsf_restart)
         Atm(n)%flagstruct%srf_init = .true.
       else
         call mpp_error(NOTE,'==> Warning from remap_restart: Expected file '//trim(fname)//' does not exist')
         Atm(n)%flagstruct%srf_init = .false.
       endif

       if ( Atm(n)%flagstruct%fv_land ) then
!--- restore data for mg_drag - if it exists
         fname = 'INPUT/mg_drag.res'//trim(stile_name)//'.nc'
         if (file_exists(fname)) then
           call read_restart(Atm(n)%Mg_restart)
           call close_file(Atm(n)%Mg_restart)
         else
           call mpp_error(NOTE,'==> Warning from remap_restart: Expected file '//trim(fname)//' does not exist')
         endif
!--- restore data for fv_land - if it exists
         fname = 'INPUT/fv_land.res'//trim(stile_name)//'.nc'
         if (file_exists(fname)) then
           call read_restart(Atm(n)%Lnd_restart)
           call close_file(Atm(n)%Lnd_restart)
         else
           call mpp_error(NOTE,'==> Warning from remap_restart: Expected file '//trim(fname)//' does not exist')
         endif
       endif

       fname = 'fv_tracer.res'//trim(stile_name)//'.nc'
       if (file_exists('INPUT/'//trim(fname))) then
         if (open_file(Tra_restart_r, fname, 'read', fv_domain, is_restart=.true.)) then
           do nt = 1, ntprog
              call get_tracer_names(MODEL_ATMOS, nt, tracer_name)
              call set_tracer_profile (MODEL_ATMOS, nt, q_r(isc:iec,jsc:jec,:,nt)  )
              if (variable_exists(Tra_restart_r, tracer_name)) &
                call register_restart_field(Tra_restart_r, tracer_name, q_r(:,:,:,nt))
           enddo
           do nt = ntprog+1, ntracers
              call get_tracer_names(MODEL_ATMOS, nt, tracer_name)
              call set_tracer_profile (MODEL_ATMOS, nt, qdiag_r(isc:iec,jsc:jec,:,nt)  )
              if (variable_exists(Tra_restart_r, tracer_name)) &
                call register_restart_field(Tra_restart_r, tracer_name, qdiag_r(:,:,:,nt))
           enddo
           call read_restart(Tra_restart_r)
           call close_file(Tra_restart_r)
         endif
       else
         call mpp_error(NOTE,'==> Warning from remap_restart: Expected file '//trim(fname)//' does not exist')
       endif

       call rst_remap(npz_rst, npz, isc, iec, jsc, jec, isd, ied, jsd, jed, ntracers, ntprog,      &
                      delp_r,      u_r,      v_r,      w_r,      delz_r,      pt_r,  q_r,  qdiag_r,&
                      Atm(n)%delp, Atm(n)%u, Atm(n)%v, Atm(n)%w, Atm(n)%delz, Atm(n)%pt, Atm(n)%q, &
                      Atm(n)%qdiag, ak_r,  bk_r, Atm(n)%ptop, Atm(n)%ak, Atm(n)%bk,                &
                      Atm(n)%flagstruct%hydrostatic, Atm(n)%flagstruct%make_nh, Atm(n)%domain,     &
                      Atm(n)%gridstruct%square_domain)
    !end do

    deallocate( ak_r )
    deallocate( bk_r )
    deallocate( u_r )
    deallocate( v_r )
    deallocate( pt_r )
    deallocate( delp_r )
    deallocate( q_r )
    deallocate( qdiag_r )

    if ( (.not.Atm(1)%flagstruct%hydrostatic) .and. (.not.Atm(1)%flagstruct%make_nh) ) then
         deallocate ( w_r )
         deallocate ( delz_r )
         if ( Atm(1)%flagstruct%hybrid_z ) deallocate ( ze0_r )
    endif

  end subroutine  remap_restart


  !#####################################################################
  ! <SUBROUTINE NAME="fv_io_register_nudge_restart">
  !
  ! <DESCRIPTION>
  !   register restart nudge field to be written out to restart file. 
  ! </DESCRIPTION>
  subroutine  fv_io_register_nudge_restart(Atm)
    type(fv_atmos_type), intent(inout) :: Atm(:)
    character(len=64) :: fname

! use_ncep_sst may not be initialized at this point?
    call mpp_error(NOTE, 'READING FROM SST_restart DISABLED')
!!$    if ( use_ncep_sst .or. Atm(1)%nudge .or. Atm(1)%ncep_ic ) then
!!$       fname = 'sst_ncep.res.nc'
!!$       id_restart = register_restart_field(Atm(1)%SST_restart, fname, 'sst_ncep', sst_ncep)
!!$       id_restart = register_restart_field(Atm(1)%SST_restart, fname, 'sst_anom', sst_anom)
!!$    endif

  end subroutine  fv_io_register_nudge_restart
  ! </SUBROUTINE> NAME="fv_io_register_nudge_restart"


  !#####################################################################
  ! <SUBROUTINE NAME="fv_io_register_restart">
  !
  ! <DESCRIPTION>
  !   register restart field to be written out to restart file. 
  ! </DESCRIPTION>
  subroutine  fv_io_register_restart(fv_domain,Atm)
    type(domain2d),      intent(inout) :: fv_domain
    type(fv_atmos_type), intent(inout) :: Atm(:)

    character(len=64) :: fname, tracer_name
    character(len=6)  :: gn, stile_name
    integer           :: n, nt, ntracers, ntprog, ntdiag, ntileMe, ntiles

    ntileMe = size(Atm(:)) 
    ntprog = size(Atm(1)%q,4) 
    ntdiag = size(Atm(1)%qdiag,4) 
    ntracers = ntprog+ntdiag

!--- set the 'nestXX' appendix for all files using fms_io
    if (Atm(1)%grid_number > 1) then
       write(gn,'(A4, I2.2)') "nest", Atm(1)%grid_number
    else
       gn = ''
    end if
    call set_filename_appendix(gn)

!--- fix for single tile runs where you need fv_core.res.nc and fv_core.res.tile1.nc
    ntiles = mpp_get_ntile_count(fv_domain)
    if(ntiles == 1 .and. .not. Atm(1)%neststruct%nested) then
       stile_name = '.tile1'
    else
       stile_name = ''
    endif

! use_ncep_sst may not be initialized at this point?
#ifndef DYCORE_SOLO
    call mpp_error(NOTE, 'READING FROM SST_RESTART DISABLED')
!!$   if ( use_ncep_sst .or. Atm(1)%flagstruct%nudge .or. Atm(1)%flagstruct%ncep_ic ) then
!!$       fname = 'sst_ncep'//trim(gn)//'.res.nc'
!!$       id_restart = register_restart_field(Atm(1)%SST_restart, fname, 'sst_ncep', sst_ncep)
!!$       id_restart = register_restart_field(Atm(1)%SST_restart, fname, 'sst_anom', sst_anom)
!!$   endif
#endif

    fname = 'fv_core.res.nc'
    if (open_file(Atm(1)%Fv_restart, fname, 'read', is_restart=.true.)) then
      call register_restart_field(Atm(1)%Fv_restart, 'ak', Atm(1)%ak(:))
      call register_restart_field(Atm(1)%Fv_restart, 'bk', Atm(1)%bk(:))
    endif

    do n = 1, ntileMe
       fname = 'fv_core.res'//trim(stile_name)//'.nc'
       if (open_file(Atm(n)%Fv_tile_restart, fname, 'read', fv_domain, is_restart=.true.)) then
         call register_restart_field(Atm(n)%Fv_tile_restart, 'u', Atm(n)%u)
         call register_restart_field(Atm(n)%Fv_tile_restart, 'v', Atm(n)%v)
         if (.not.Atm(n)%flagstruct%hydrostatic) then
            if (variable_exists(Atm(n)%Fv_tile_restart, 'W')) &
              call register_restart_field(Atm(n)%Fv_tile_restart, 'W', Atm(n)%w)
            if (variable_exists(Atm(n)%Fv_tile_restart, 'DZ')) &
              call register_restart_field(Atm(n)%Fv_tile_restart, 'DZ', Atm(n)%delz)
            if ( Atm(n)%flagstruct%hybrid_z ) then
               if (variable_exists(Atm(n)%Fv_tile_restart, 'ZE0')) &
                 call register_restart_field(Atm(n)%Fv_tile_restart, 'ZE0', Atm(n)%ze0)
            endif
         endif
         call register_restart_field(Atm(n)%Fv_tile_restart, 'T', Atm(n)%pt)
         call register_restart_field(Atm(n)%Fv_tile_restart, 'delp', Atm(n)%delp)
         call register_restart_field(Atm(n)%Fv_tile_restart, 'phis', Atm(n)%phis)

         !--- include agrid winds in restarts for use in data assimilation 
         if (Atm(n)%flagstruct%agrid_vel_rst) then
           if (variable_exists(Atm(n)%Fv_tile_restart, 'ua')) &
             call register_restart_field(Atm(n)%Fv_tile_restart, 'ua', Atm(n)%ua)
           if (variable_exists(Atm(n)%Fv_tile_restart, 'va')) &
             call register_restart_field(Atm(n)%Fv_tile_restart, 'va', Atm(n)%va)
         endif
       endif

       fname = 'fv_srf_wnd.res'//trim(stile_name)//'.nc'
       if (open_file(Atm(n)%Rsf_restart, fname, 'read', fv_domain, is_restart=.true.)) then
         call register_restart_field(Atm(n)%Rsf_restart, 'u_srf', Atm(n)%u_srf)
         call register_restart_field(Atm(n)%Rsf_restart, 'v_srf', Atm(n)%v_srf)
       endif
#ifdef SIM_PHYS
       call register_restart_field(Atm(n)%Rsf_restart, 'ts', Atm(n)%ts)
#endif

       if ( Atm(n)%flagstruct%fv_land ) then
          !-------------------------------------------------------------------------------------------------
          ! Optional terrain deviation (sgh) and land fraction (oro)
          fname = 'mg_drag.res'//trim(stile_name)//'.nc'
          if (open_file(Atm(n)%Mg_restart, fname, 'read', fv_domain, is_restart=.true.)) then
            call register_restart_field(Atm(n)%Mg_restart, 'ghprime', Atm(n)%sgh)
          endif

          fname = 'fv_land.res'//trim(stile_name)//'.nc'
          if (open_file(Atm(n)%Lnd_restart, fname, 'read', fv_domain, is_restart=.true.)) then
            call register_restart_field(Atm(n)%Lnd_restart, 'oro', Atm(n)%oro)
          endif
       endif

       fname = 'fv_tracer.res'//trim(stile_name)//'.nc'
       if (open_file(Atm(n)%Tra_restart, fname, 'read', fv_domain, is_restart=.true.)) then
         do nt = 1, ntprog
            call get_tracer_names(MODEL_ATMOS, nt, tracer_name)
            ! set all tracers to an initial profile value
            call set_tracer_profile (MODEL_ATMOS, nt, Atm(n)%q(:,:,:,nt)  )
            if (variable_exists(Atm(n)%Tra_restart, tracer_name)) &
              call register_restart_field(Atm(n)%Tra_restart, tracer_name, Atm(n)%q(:,:,:,nt))
         enddo
         do nt = ntprog+1, ntracers
            call get_tracer_names(MODEL_ATMOS, nt, tracer_name)
            ! set all tracers to an initial profile value
            call set_tracer_profile (MODEL_ATMOS, nt, Atm(n)%qdiag(:,:,:,nt)  )
            if (variable_exists(Atm(n)%Tra_restart, tracer_name)) &
              call register_restart_field(Atm(n)%Tra_restart, tracer_name, Atm(n)%qdiag(:,:,:,nt))
         enddo
       endif

    enddo

  end subroutine  fv_io_register_restart
  ! </SUBROUTINE> NAME="fv_io_register_restart"



  !#####################################################################
  ! <SUBROUTINE NAME="fv_io_write_restart">
  !
  ! <DESCRIPTION>
  ! Write the fv core restart quantities 
  ! </DESCRIPTION>
  subroutine  fv_io_write_restart(Atm, grids_on_this_pe, timestamp)

    type(fv_atmos_type),        intent(inout) :: Atm(:)
    logical, intent(IN) :: grids_on_this_pe(:)
    character(len=*), optional, intent(in) :: timestamp
    integer                                :: n, ntileMe

    ntileMe = size(Atm(:))  ! This will need mods for more than 1 tile per pe

    if ( use_ncep_sst .or. Atm(1)%flagstruct%nudge .or. Atm(1)%flagstruct%ncep_ic ) then
       call mpp_error(NOTE, 'READING FROM SST_RESTART DISABLED')
       !call save_restart(Atm(1)%SST_restart, timestamp)
    endif
 
    do n = 1, ntileMe
       if (.not. grids_on_this_pe(n)) cycle

       if ( (use_ncep_sst .or. Atm(n)%flagstruct%nudge) .and. .not. Atm(n)%gridstruct%nested ) then
          call write_restart(Atm(n)%SST_restart)
          call close_file(Atm(n)%SST_restart)
       endif
 
       call write_restart(Atm(n)%Fv_restart)
       call close_file(Atm(n)%Fv_restart)
       call write_restart(Atm(n)%Fv_tile_restart)
       call close_file(Atm(n)%Fv_tile_restart)
       call write_restart(Atm(n)%Rsf_restart)
       call close_file(Atm(n)%Rsf_restart)

       if ( Atm(n)%flagstruct%fv_land ) then
          call write_restart(Atm(n)%Mg_restart)
          call close_file(Atm(n)%Mg_restart)
          call write_restart(Atm(n)%Lnd_restart)
          call close_file(Atm(n)%Lnd_restart)
       endif

       call write_restart(Atm(n)%Tra_restart)
       call close_file(Atm(n)%Tra_restart)

    end do

  end subroutine  fv_io_write_restart

  subroutine register_bcs_2d(Atm, BCfile_ne, BCfile_sw, fname_ne, fname_sw, &
                             var_name, var, var_bc, istag, jstag)
    type(fv_atmos_type),           intent(in)    :: Atm
    type(FmsNetcdfDomainFile_t),   intent(inout) :: BCfile_ne, BCfile_sw
    character(len=120),            intent(in)    :: fname_ne, fname_sw
    character(len=*),              intent(in)    :: var_name
    real, dimension(:,:),          intent(in), optional :: var
    type(fv_nest_BC_type_3D),      intent(in), optional :: var_bc
    integer,                       intent(in), optional :: istag, jstag

!register west halo data in t1
    if (present(var_bc)) call register_restart_field(BCfile_sw, trim(var_name)//'_west_t1', var_bc%west_t1)
!register west prognostic halo data
    if (present(var)) call register_restart_field(BCfile_sw, trim(var_name)//'_west', var)

!register east halo data in t1
    if (present(var_bc)) call register_restart_field(BCfile_ne, trim(var_name)//'_east_t1', var_bc%east_t1)
!register east prognostic halo data
    if (present(var)) call register_restart_field(BCfile_ne, trim(var_name)//'_east', var)

!register south halo data in t1
    if (present(var_bc)) call register_restart_field(BCfile_sw, trim(var_name)//'_south_t1', var_bc%south_t1)
!register south prognostic halo data
    if (present(var)) call register_restart_field(BCfile_sw, trim(var_name)//'_south', var)

!register north halo data in t1
    if (present(var_bc)) call register_restart_field(BCfile_ne, trim(var_name)//'_north_t1', var_bc%north_t1)
!register north prognostic halo data
    if (present(var)) call register_restart_field(BCfile_ne, trim(var_name)//'_north', var)

  end subroutine register_bcs_2d


  subroutine register_bcs_3d(Atm, BCfile_ne, BCfile_sw, fname_ne, fname_sw, &
                             var_name, var, var_bc, istag, jstag, mandatory)
    type(fv_atmos_type),           intent(in)    :: Atm
    type(FmsNetcdfDomainFile_t),   intent(inout) :: BCfile_ne, BCfile_sw
    character(len=120),            intent(in)    :: fname_ne, fname_sw
    character(len=*),              intent(in)    :: var_name
    real, dimension(:,:,:),        intent(in), optional :: var
    type(fv_nest_BC_type_3D),      intent(in), optional :: var_bc
    integer,                       intent(in), optional :: istag, jstag
    logical,                       intent(IN), optional :: mandatory

!register west halo data in t1
    if (present(var_bc)) call register_restart_field(BCfile_sw, trim(var_name)//'_west_t1', var_bc%west_t1)
!register west prognostic halo data
    if (present(var)) call register_restart_field(BCfile_sw, trim(var_name)//'_west', var)

!register east halo data in t1
    if (present(var_bc)) call register_restart_field(BCfile_ne, trim(var_name)//'_east_t1', var_bc%east_t1)
!register east prognostic halo data
    if (present(var)) call register_restart_field(BCfile_ne, trim(var_name)//'_east', var)

!register south halo data in t1
    if (present(var_bc)) call register_restart_field(BCfile_sw, trim(var_name)//'_south_t1', var_bc%south_t1)
!register south prognostic halo data
    if (present(var)) call register_restart_field(BCfile_sw, trim(var_name)//'_south', var)

!register north halo data in t1
    if (present(var_bc)) call register_restart_field(BCfile_ne, trim(var_name)//'_north_t1', var_bc%north_t1)
!register north prognostic halo data
    if (present(var)) call register_restart_field(BCfile_ne, trim(var_name)//'_north', var)

  end subroutine register_bcs_3d


  ! </SUBROUTINE> NAME="fv_io_regsiter_restart_BCs"
  !#####################################################################

  subroutine fv_io_register_restart_BCs(Atm)
    type(fv_atmos_type),        intent(inout) :: Atm

    integer :: n, ntracers, ntprog, ntdiag
    character(len=120) :: tname, fname_ne, fname_sw
    type(FmsNetcdfDomainFile_t) :: fileobj_tmp

    fname_ne = 'fv_BC_ne.res.nc'
    fname_sw = 'fv_BC_sw.res.nc'

    ntprog=size(Atm%q,4)
    ntdiag=size(Atm%qdiag,4)
    ntracers=ntprog+ntdiag

    call set_domain(Atm%domain)

    if (open_file(Atm%neststruct%BCfile_ne, fname_ne, 'read', Atm%domain, is_restart=.true.)) then
    endif
    if (open_file(Atm%neststruct%BCfile_sw, fname_sw, 'read', Atm%domain, is_restart=.true.)) then
    endif

    call register_bcs_2d(Atm, Atm%neststruct%BCfile_ne, Atm%neststruct%BCfile_sw, &
                         fname_ne, fname_sw, 'phis', var=Atm%phis)
    call register_bcs_3d(Atm, Atm%neststruct%BCfile_ne, Atm%neststruct%BCfile_sw, &
                         fname_ne, fname_sw, 'delp', Atm%delp, Atm%neststruct%delp_BC)
    do n=1,ntprog
       call get_tracer_names(MODEL_ATMOS, n, tname)
       call register_bcs_3d(Atm, Atm%neststruct%BCfile_ne, Atm%neststruct%BCfile_sw, &
                            fname_ne, fname_sw, trim(tname), Atm%q(:,:,:,n), Atm%neststruct%q_BC(n), mandatory=.false.)
    enddo
    do n=ntprog+1,ntracers
       call get_tracer_names(MODEL_ATMOS, n, tname)
       call register_bcs_3d(Atm, Atm%neststruct%BCfile_ne, Atm%neststruct%BCfile_sw, &
                            fname_ne, fname_sw, trim(tname), var=Atm%qdiag(:,:,:,n), mandatory=.false.)
    enddo
#ifndef SW_DYNAMICS
    call register_bcs_3d(Atm, Atm%neststruct%BCfile_ne, Atm%neststruct%BCfile_sw, &
                         fname_ne, fname_sw, 'pt', Atm%pt, Atm%neststruct%pt_BC)
    if ((.not.Atm%flagstruct%hydrostatic) .and. (.not.Atm%flagstruct%make_nh)) then
       if (is_master()) print*, 'fv_io_register_restart_BCs: REGISTERING NH BCs', Atm%flagstruct%hydrostatic, Atm%flagstruct%make_nh
      call register_bcs_3d(Atm, Atm%neststruct%BCfile_ne, Atm%neststruct%BCfile_sw, &
                           fname_ne, fname_sw, 'w', Atm%w, Atm%neststruct%w_BC, mandatory=.false.)
      call register_bcs_3d(Atm, Atm%neststruct%BCfile_ne, Atm%neststruct%BCfile_sw, &
                           fname_ne, fname_sw, 'delz', Atm%delz, Atm%neststruct%delz_BC, mandatory=.false.)
    endif
#ifdef USE_COND
       call register_bcs_3d(Atm, Atm%neststruct%BCfile_ne, Atm%neststruct%BCfile_sw, &
                            fname_ne, fname_sw,'q_con', var_bc=Atm%neststruct%q_con_BC, mandatory=.false.)
#ifdef MOIST_CAPPA
       call register_bcs_3d(Atm, Atm%neststruct%BCfile_ne, Atm%neststruct%BCfile_sw, &
            fname_ne, fname_sw, 'cappa', var_bc=Atm%neststruct%cappa_BC, mandatory=.false.)
#endif
#endif
#endif
    call register_bcs_3d(Atm, Atm%neststruct%BCfile_ne, Atm%neststruct%BCfile_sw, &
                         fname_ne, fname_sw, 'u', Atm%u, Atm%neststruct%u_BC, jstag=1)
    call register_bcs_3d(Atm, Atm%neststruct%BCfile_ne, Atm%neststruct%BCfile_sw, &
                         fname_ne, fname_sw, 'v', Atm%v, Atm%neststruct%v_BC, istag=1)
    call register_bcs_3d(Atm, Atm%neststruct%BCfile_ne, Atm%neststruct%BCfile_sw, &
                         fname_ne, fname_sw, 'uc', var_bc=Atm%neststruct%uc_BC, istag=1)
    call register_bcs_3d(Atm, Atm%neststruct%BCfile_ne, Atm%neststruct%BCfile_sw, &
                         fname_ne, fname_sw, 'vc', var_bc=Atm%neststruct%vc_BC, jstag=1)
    call register_bcs_3d(Atm, Atm%neststruct%BCfile_ne, Atm%neststruct%BCfile_sw, &
                         fname_ne, fname_sw, 'divg', var_bc=Atm%neststruct%divg_BC, istag=1,jstag=1, mandatory=.false.)
    Atm%neststruct%divg_BC%initialized = variable_exists(Atm%neststruct%BCfile_ne, 'divg_north_t1')


    return
  end subroutine fv_io_register_restart_BCs


  subroutine fv_io_register_restart_BCs_NH(Atm)
    type(fv_atmos_type),        intent(inout) :: Atm

    integer :: n
    character(len=120) :: tname, fname_ne, fname_sw

    fname_ne = 'fv_BC_ne.res.nc'
    fname_sw = 'fv_BC_sw.res.nc'

    call set_domain(Atm%domain)

    if (is_master()) print*, 'fv_io_register_restart_BCs_NH: REGISTERING NH BCs', Atm%flagstruct%hydrostatic, Atm%flagstruct%make_nh
#ifndef SW_DYNAMICS
    call register_bcs_3d(Atm, Atm%neststruct%BCfile_ne, Atm%neststruct%BCfile_sw, &
         fname_ne, fname_sw, 'w', Atm%w, Atm%neststruct%w_BC)
    call register_bcs_3d(Atm, Atm%neststruct%BCfile_ne, Atm%neststruct%BCfile_sw, &
         fname_ne, fname_sw, 'delz', Atm%delz, Atm%neststruct%delz_BC)
#endif

    return
  end subroutine fv_io_register_restart_BCs_NH


  subroutine fv_io_write_BCs(Atm, timestamp)
    type(fv_atmos_type), intent(inout) :: Atm
    character(len=*),    intent(in), optional :: timestamp

    call write_restart(Atm%neststruct%BCfile_ne)
    call close_file(Atm%neststruct%BCfile_ne)
    call write_restart(Atm%neststruct%BCfile_sw)
    call close_file(Atm%neststruct%BCfile_sw)

    return
  end subroutine fv_io_write_BCs


  subroutine fv_io_read_BCs(Atm)
    type(fv_atmos_type), intent(inout) :: Atm

    call read_restart(Atm%neststruct%BCfile_ne)
    call close_file(Atm%neststruct%BCfile_ne)
    call read_restart(Atm%neststruct%BCfile_sw)
    call close_file(Atm%neststruct%BCfile_sw)

    return
  end subroutine fv_io_read_BCs

end module fv_io_nlm_mod
