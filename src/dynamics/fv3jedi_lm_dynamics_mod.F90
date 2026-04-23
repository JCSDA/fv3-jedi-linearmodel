module fv3jedi_lm_dynamics_mod

use fv3jedi_lm_utils_mod
use fv3jedi_lm_kinds_mod
use fv3jedi_lm_const_mod

use fms_mod,         only: set_domain, nullify_domain
use fms2_io_mod,     only: FmsNetcdfDomainFile_t, open_file, close_file, &
                           register_restart_field, register_axis, &
                           read_restart, write_restart, unlimited, &
                           register_variable_attribute
use mpp_mod,         only: mpp_pe, mpp_root_pe, mpp_error, FATAL
use mpp_domains_mod, only: mpp_update_domains, mpp_get_boundary, DGRID_NE, mpp_get_boundary_ad, &
                           east, north, center

use fv_control_nlm_mod,     only: fv_init, pelist_all
use fv_control_tlmadm_mod,  only: fv_init_pert
use fv_arrays_nlm_mod,      only: fv_atmos_type, deallocate_fv_atmos_type
use fv_arrays_tlmadm_mod,   only: fv_atmos_pert_type, deallocate_fv_atmos_pert_type
use fv_dynamics_nlm_mod,    only: fv_dynamics
use fv_dynamics_tlm_mod,    only: fv_dynamics_tlm, fv_dynamics_nlm => fv_dynamics
use fv_dynamics_adm_mod,    only: fv_dynamics_fwd, fv_dynamics_bwd
use fv_pressure_mod,        only: compute_fv3_pressures, compute_fv3_pressures_tlm, compute_fv3_pressures_bwd



use tapenade_iter, only: cp_iter, cp_iter_controls, initialize_cp_iter, finalize_cp_iter
use tapenade_iter, only: cp_mod_ini, cp_mod_mid, cp_mod_end, pushrealarray, poprealarray

use fv3jedi_lm_dynutils_mod, only: init_fv_diag_type

!> Top level for fv3jedi linearized dynamical core

implicit none
private
public :: fv3jedi_lm_dynamics_type

! Precision of dyncore
#ifdef SINGLE_FV
  integer, parameter :: fvprec = 4
#else
  integer, parameter :: fvprec = 8
#endif

type fv3jedi_lm_dynamics_type
 type(fv_atmos_type),      allocatable :: FV_Atm(:)       !<Traj FV3 structure
 type(fv_atmos_pert_type), allocatable :: FV_AtmP(:)      !<Pert FV3 structure
 real(fvprec), allocatable, dimension(:,:) :: ebuffery    !<Halo holder
 real(fvprec), allocatable, dimension(:,:) :: nbufferx    !<Halo holder
 real(fvprec), allocatable, dimension(:,:) :: wbuffery    !<Halo holder
 real(fvprec), allocatable, dimension(:,:) :: sbufferx    !<Halo holder
 integer         :: isc,iec,jsc,jec     !<Grid, compute region
 integer         :: isd,ied,jsd,jed     !<Grid, with halo
 integer         :: npz                 !<Number of vertical levels
 integer :: cp_dyn_ind
 integer :: linmodtest = 0
 contains
  procedure :: create
  procedure :: init_nl
  procedure :: init_tl
  procedure :: init_ad
  procedure :: step_nl
  procedure :: step_tl
  procedure :: step_ad
  procedure :: delete
  procedure :: traj_to_fv3
  procedure :: fv3_to_traj
  procedure :: pert_to_fv3
  procedure :: fv3_to_pert
  procedure :: read_d_grid_winds
  procedure :: write_d_grid_winds
  procedure :: a_to_d_nl
end type fv3jedi_lm_dynamics_type

contains

! ------------------------------------------------------------------------------

subroutine create(self,conf)

 implicit none

 class(fv3jedi_lm_dynamics_type), target, intent(inout) :: self
 type(fv3jedi_lm_conf), intent(inout)    :: conf

 logical, allocatable :: grids_on_this_pe(:)
 integer :: p_split = 1
 integer :: i,j,tmp
 real(kind_real) :: f_coriolis_angle

 type(fv_atmos_type), pointer :: FV_Atm(:)

  call fv_init(self%FV_Atm, real(conf%dt, fvprec), grids_on_this_pe, p_split)

  FV_Atm => self%FV_Atm

  if (allocated(grids_on_this_pe)) deallocate(grids_on_this_pe)
  if (allocated(pelist_all)) deallocate(pelist_all)

  !Halo holders for domain grid
  allocate(self%wbuffery(FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,FV_Atm(1)%npz))
  allocate(self%sbufferx(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%npz))
  allocate(self%ebuffery(FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,FV_Atm(1)%npz))
  allocate(self%nbufferx(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%npz))

  !Set ptop, ak, bk in fv3 structure
  FV_Atm(1)%ak = conf%ak
  FV_Atm(1)%bk = conf%bk
  FV_Atm(1)%ptop = conf%ptop

  !Always allocate w, delz, q_con for now
  deallocate(FV_Atm(1)%w)
  deallocate(FV_Atm(1)%delz)
  deallocate(FV_Atm(1)%q_con)
  allocate  ( FV_Atm(1)%w (FV_Atm(1)%bd%isd:FV_Atm(1)%bd%ied,FV_Atm(1)%bd%jsd:FV_Atm(1)%bd%jed,&
                                 FV_Atm(1)%flagstruct%npz) )
  allocate  ( FV_Atm(1)%delz (FV_Atm(1)%bd%isd:FV_Atm(1)%bd%ied,FV_Atm(1)%bd%jsd:FV_Atm(1)%bd%jed,&
                                 FV_Atm(1)%flagstruct%npz) )
  allocate  ( FV_Atm(1)%q_con(FV_Atm(1)%bd%isd:FV_Atm(1)%bd%ied,FV_Atm(1)%bd%jsd:FV_Atm(1)%bd%jed,&
                                 FV_Atm(1)%flagstruct%npz) )
  FV_Atm(1)%w = 0.0_kind_real
  FV_Atm(1)%delz = 0.0_kind_real
  FV_Atm(1)%q_con = 0.0_kind_real

  f_coriolis_angle = 0.0_kind_real

  !fC and f0
  if (FV_Atm(1)%flagstruct%grid_type == 4) then
     FV_Atm(1)%gridstruct%fC(:,:) = 2.0_kind_real*omega*sin(FV_Atm(1)%flagstruct%deglat/180.0_kind_real*pi)
     FV_Atm(1)%gridstruct%f0(:,:) = 2.0_kind_real*omega*sin(FV_Atm(1)%flagstruct%deglat/180.0_kind_real*pi)
  else
     if (f_coriolis_angle == -999.0_kind_real) then
        FV_Atm(1)%gridstruct%fC(:,:) = 0.0_kind_real
        FV_Atm(1)%gridstruct%f0(:,:) = 0.0_kind_real
     else
        do j=FV_Atm(1)%bd%jsd,FV_Atm(1)%bd%jed+1
           do i=FV_Atm(1)%bd%isd,FV_Atm(1)%bd%ied+1
              FV_Atm(1)%gridstruct%fC(i,j) = 2.0_kind_real*omega*( -COS(FV_Atm(1)%gridstruct%grid(i,j,1))*&
                                             COS(FV_Atm(1)%gridstruct%grid(i,j,2))*SIN(f_coriolis_angle) + &
                                             SIN(FV_Atm(1)%gridstruct%grid(i,j,2))*COS(f_coriolis_angle) )
           enddo
        enddo
        do j=FV_Atm(1)%bd%jsd,FV_Atm(1)%bd%jed
           do i=FV_Atm(1)%bd%isd,FV_Atm(1)%bd%ied
              FV_Atm(1)%gridstruct%f0(i,j) = 2.0_kind_real*omega*( -COS(FV_Atm(1)%gridstruct%agrid(i,j,1))*&
                                             COS(FV_Atm(1)%gridstruct%agrid(i,j,2))*SIN(f_coriolis_angle) + &
                                             SIN(FV_Atm(1)%gridstruct%agrid(i,j,2))*COS(f_coriolis_angle) )
           enddo
        enddo
     endif
  endif

  !Pointer to self when not nested
  if (.not. FV_Atm(1)%gridstruct%nested) FV_Atm(1)%parent_grid => FV_Atm(1)

  !Hardwire some flags
  FV_Atm(1)%flagstruct%reproduce_sum = .false.
  FV_Atm(1)%flagstruct%fill = .false.
  FV_Atm(1)%flagstruct%fv_debug = .false.
  FV_Atm(1)%flagstruct%adiabatic = .false.
  FV_Atm(1)%flagstruct%do_sat_adj = .false.
  FV_Atm(1)%flagstruct%breed_vortex_inline = .false.

  !Initialze the perturbation fv3 structure
  call fv_init_pert(self%FV_Atm,self%FV_AtmP,conf%inputpert_filename)

  !Global
  cp_iter_controls%cp_i  = 0
  cp_iter_controls%cp_nt = 4
  cp_iter_controls%cp_gb = -0.1
  cp_iter_controls%cp_nm = 1
  call initialize_cp_iter

  if (cp_iter_controls%cp_i .ne. 0) then

     !Dynamics
     self%cp_dyn_ind = 1
     cp_iter(self%cp_dyn_ind)%my_name(1:3) = 'dyn'

     cp_iter(self%cp_dyn_ind)%cp_test = .false.
     tmp = 0
     if (tmp==1) cp_iter(self%cp_dyn_ind)%cp_test = .true.

     cp_iter(self%cp_dyn_ind)%cp_rep = .false.
     tmp = 0
     if (tmp==1) cp_iter(self%cp_dyn_ind)%cp_test = .true.

     !Hardwire these for now
     cp_iter(self%cp_dyn_ind)%check_st_control = .false.
     cp_iter(self%cp_dyn_ind)%check_st_integer = .false.
     cp_iter(self%cp_dyn_ind)%check_st_real_r4 = .false.
     cp_iter(self%cp_dyn_ind)%check_st_real_r8 = .false.

     cp_iter(self%cp_dyn_ind)%test_dim_st_control = 0
     cp_iter(self%cp_dyn_ind)%test_dim_st_integer = 0
     cp_iter(self%cp_dyn_ind)%test_dim_st_real_r4 = 0
     cp_iter(self%cp_dyn_ind)%test_dim_st_real_r8 = 0

     cp_iter(self%cp_dyn_ind)%test_dim_cp_control = 0
     cp_iter(self%cp_dyn_ind)%test_dim_cp_integer = 0
     cp_iter(self%cp_dyn_ind)%test_dim_cp_real_r4 = 0
     cp_iter(self%cp_dyn_ind)%test_dim_cp_real_r8 = 0

  endif

  !Convenience
  self%isc = FV_Atm(1)%bd%isc
  self%iec = FV_Atm(1)%bd%iec
  self%jsc = FV_Atm(1)%bd%jsc
  self%jec = FV_Atm(1)%bd%jec
  self%isd = FV_Atm(1)%bd%isd
  self%ied = FV_Atm(1)%bd%ied
  self%jsd = FV_Atm(1)%bd%jsd
  self%jed = FV_Atm(1)%bd%jed
  self%npz = FV_Atm(1)%npz

  conf%rpe = .false.
  if (mpp_pe() == mpp_root_pe()) conf%rpe = .true.

  ! Initialize the idiag structure to zero for safety
  call init_fv_diag_type(FV_Atm(1)%idiag)

endsubroutine create

! ------------------------------------------------------------------------------

subroutine init_nl(self,conf,pert,traj)

 implicit none

 class(fv3jedi_lm_dynamics_type), intent(inout), target :: self
 type(fv3jedi_lm_conf), intent(in) :: conf
 type(fv3jedi_lm_pert), intent(inout) :: pert
 type(fv3jedi_lm_traj), intent(in) :: traj
 type(fv_atmos_type), pointer :: FV_Atm(:)

 !Convenience pointer to the main FV_Atm structure
 !------------------------------------------------
 FV_Atm => self%FV_Atm

 FV_Atm(1)%flagstruct%ncnst = size(traj%tracers,4)
 FV_Atm(1)%ncnst = size(traj%tracers,4)

 call allocate_tracers(FV_Atm)

endsubroutine init_nl

! ------------------------------------------------------------------------------

subroutine read_d_grid_winds(self, dpath, fname, traj)

  ! Args
  class(fv3jedi_lm_dynamics_type), intent(inout) :: self
  character(len=*), intent(in) :: dpath, fname
  type(fv3jedi_lm_traj), intent(inout) :: traj

  ! Locals
  integer :: isc, iec, jsc, jec, npz
  type(FmsNetcdfDomainFile_t) :: rst
  character(len=512) :: fpath
  real(kind=kind_real), allocatable, dimension(:,:,:) :: u, v

  ! Convenience
  isc = self%isc
  iec = self%iec
  jsc = self%jsc
  jec = self%jec
  npz = self%npz

  allocate ( u(isc:iec,  jsc:jec+1,npz) )
  allocate ( v(isc:iec+1,jsc:jec  ,npz) )
  u = 0.0_kind_real
  v = 0.0_kind_real

  ! Build the full file path
  fpath = trim(adjustl(dpath))//'/'//trim(fname)

  ! Open, register, read, close
  if (open_file(rst, trim(fpath), 'read', self%FV_Atm(1)%domain, is_restart=.true.)) then
    call register_axis(rst, 'xaxis_1', 'x', domain_position=center)
    call register_axis(rst, 'xaxis_2', 'x', domain_position=east)
    call register_axis(rst, 'yaxis_1', 'y', domain_position=north)
    call register_axis(rst, 'yaxis_2', 'y', domain_position=center)
    call register_axis(rst, 'zaxis_1', npz)
    call register_axis(rst, 'Time', unlimited)
     
    call register_restart_field(rst, 'u', u, (/'xaxis_1', 'yaxis_1', 'zaxis_1', 'Time   '/))
    call register_restart_field(rst, 'v', v, (/'xaxis_2', 'yaxis_2', 'zaxis_1', 'Time   '/))

    call read_restart(rst)
    call close_file(rst)
  endif

  ! Now copy the data to the traj structure
  traj%u(isc:iec, jsc:jec, 1:npz) = u(isc:iec, jsc:jec, 1:npz)
  traj%v(isc:iec, jsc:jec, 1:npz) = v(isc:iec, jsc:jec, 1:npz)

endsubroutine read_d_grid_winds

! ------------------------------------------------------------------------------

subroutine write_d_grid_winds(self, dpath, fname, traj)

  ! Args
  class(fv3jedi_lm_dynamics_type), intent(inout) :: self
  character(len=*), intent(in) :: dpath, fname
  type(fv3jedi_lm_traj), intent(inout) :: traj

  ! Locals
  integer :: isc, iec, jsc, jec, npz
  type(FmsNetcdfDomainFile_t) :: rst
  real(kind=kind_real), allocatable, dimension(:,:,:) :: u, v
  character(len=512) :: fpath

  ! Convenience
  isc = self%isc
  iec = self%iec
  jsc = self%jsc
  jec = self%jec
  npz = self%npz

  ! Allocate D-Grid winds with the correct array shapes and copy from traj
  allocate ( u(isc:iec,  jsc:jec+1,npz) )
  allocate ( v(isc:iec+1,jsc:jec  ,npz) )
  u = 0.0_kind_real
  v = 0.0_kind_real
  u(isc:iec, jsc:jec, 1:npz) = traj%u(isc:iec, jsc:jec, 1:npz)
  v(isc:iec, jsc:jec, 1:npz) = traj%v(isc:iec, jsc:jec, 1:npz)

  ! Build full file path
  fpath = trim(adjustl(dpath))//'/'//trim(fname)

  ! Open, register axes and fields, write, close
  if (open_file(rst, trim(fpath), 'overwrite', self%FV_Atm(1)%domain, is_restart=.true.)) then
    call register_axis(rst, 'xaxis_1', 'x', domain_position=center)
    call register_axis(rst, 'xaxis_2', 'x', domain_position=east)
    call register_axis(rst, 'yaxis_1', 'y', domain_position=north)
    call register_axis(rst, 'yaxis_2', 'y', domain_position=center)
    call register_axis(rst, 'zaxis_1', npz)
    call register_axis(rst, 'Time', unlimited)

    call register_restart_field(rst, 'u', u, (/'xaxis_1', 'yaxis_1', 'zaxis_1', 'Time   '/))
    call register_variable_attribute(rst, 'u', 'long_name', 'u_component_of_native_D_grid_wind', str_len=len('u_component_of_native_D_grid_wind'))
    call register_variable_attribute(rst, 'u', 'units', 'ms-1', str_len=len('ms-1'))
 
    call register_restart_field(rst, 'v', v, (/'xaxis_2', 'yaxis_2', 'zaxis_1', 'Time   '/))
    call register_variable_attribute(rst, 'v', 'units', 'ms-1', str_len=len('ms-1'))
    call register_variable_attribute(rst, 'v', 'long_name', 'v_component_of_native_D_grid_wind', str_len=len('v_component_of_native_D_grid_wind'))

    call write_restart(rst)
    call close_file(rst)
  endif

endsubroutine write_d_grid_winds

! ------------------------------------------------------------------------------

subroutine init_tl(self,conf,pert,traj)

 implicit none

 class(fv3jedi_lm_dynamics_type), intent(inout), target :: self
 type(fv3jedi_lm_conf), intent(in) :: conf
 type(fv3jedi_lm_pert), intent(inout) :: pert
 type(fv3jedi_lm_traj), intent(in) :: traj
 type(fv_atmos_type), pointer :: FV_Atm(:)
 type(fv_atmos_pert_type), pointer :: FV_AtmP(:)

 !Convenience pointer to the main FV_Atm structure
 !------------------------------------------------
 FV_Atm  => self%FV_Atm
 FV_AtmP => self%FV_AtmP

 FV_Atm(1)%flagstruct%ncnst = size(traj%tracers,4)
 FV_Atm(1)%ncnst = size(traj%tracers,4)

 ! Check for equivalency between traj and pert tracers
 call check_tracers(traj, pert)

 call allocate_tracers(FV_Atm)
 call allocate_tracersP(FV_AtmP,FV_Atm)

endsubroutine init_tl

! ------------------------------------------------------------------------------

subroutine init_ad(self,conf,pert,traj)

 implicit none

 class(fv3jedi_lm_dynamics_type), intent(inout), target :: self
 type(fv3jedi_lm_conf), intent(in) :: conf
 type(fv3jedi_lm_pert), intent(inout) :: pert
 type(fv3jedi_lm_traj), intent(in) :: traj
 type(fv_atmos_type), pointer :: FV_Atm(:)
 type(fv_atmos_pert_type), pointer :: FV_AtmP(:)

 !Convenience pointer to the main FV_Atm structure
 !------------------------------------------------
 FV_Atm  => self%FV_Atm
 FV_AtmP => self%FV_AtmP

 FV_Atm(1)%flagstruct%ncnst = size(traj%tracers,4)
 FV_Atm(1)%ncnst = size(traj%tracers,4)

 ! Check for equivalency between traj and pert tracers
 call check_tracers(traj, pert)

 call allocate_tracers(FV_Atm)
 call allocate_tracersP(FV_AtmP,FV_Atm)

endsubroutine init_ad

! ------------------------------------------------------------------------------

subroutine step_nl(self,conf,traj)

 implicit none

 class(fv3jedi_lm_dynamics_type), intent(inout), target :: self
 type(fv3jedi_lm_traj), intent(inout) :: traj
 type(fv3jedi_lm_conf), intent(in) :: conf

 type(fv_atmos_type), pointer :: FV_Atm(:)
 integer :: i,j,k

 !Convenience pointer to the main FV_Atm structure
 !------------------------------------------------
 FV_Atm => self%FV_Atm

 !Copy from traj to the fv3 structure
 !-----------------------------------
 call traj_to_fv3(self,conf,traj)

 ! MPP set domain
 ! --------------
 call set_domain(FV_Atm(1)%domain)

 !Propagate FV3 one time step
 !---------------------------
 if (self%linmodtest == 0) then
    call fv_dynamics( FV_Atm(1)%npx, FV_Atm(1)%npy, FV_Atm(1)%npz, FV_Atm(1)%ncnst, FV_Atm(1)%ng,  &
                      real(conf%dt, fvprec), FV_Atm(1)%flagstruct%consv_te, FV_Atm(1)%flagstruct%fill,           &
                      FV_Atm(1)%flagstruct%reproduce_sum, real(kappa, fvprec),                                   &
                      real(cp, fvprec), real(zvir, fvprec), FV_Atm(1)%ptop, FV_Atm(1)%ks, FV_Atm(1)%flagstruct%ncnst,          &
                      FV_Atm(1)%flagstruct%n_split, FV_Atm(1)%flagstruct%q_split,                  &
                      FV_Atm(1)%u, FV_Atm(1)%v, FV_Atm(1)%w, FV_Atm(1)%delz,                       &
                      FV_Atm(1)%flagstruct%hydrostatic, FV_Atm(1)%pt, FV_Atm(1)%delp, FV_Atm(1)%q, &
                      FV_Atm(1)%ps, FV_Atm(1)%pe, FV_Atm(1)%pk, FV_Atm(1)%peln, FV_Atm(1)%pkz,     &
                      FV_Atm(1)%phis, FV_Atm(1)%q_con, FV_Atm(1)%omga,                             &
                      FV_Atm(1)%ua, FV_Atm(1)%va, FV_Atm(1)%uc, FV_Atm(1)%vc,                      &
                      FV_Atm(1)%ak, FV_Atm(1)%bk,                                                  &
                      FV_Atm(1)%mfx, FV_Atm(1)%mfy, FV_Atm(1)%cx, FV_Atm(1)%cy, FV_Atm(1)%ze0,     &
                      FV_Atm(1)%flagstruct%hybrid_z, FV_Atm(1)%gridstruct, FV_Atm(1)%flagstruct,   &
                      FV_Atm(1)%neststruct, FV_Atm(1)%idiag, FV_Atm(1)%bd, FV_Atm(1)%parent_grid,  &
                      FV_Atm(1)%domain )
 else
    call fv_dynamics_nlm( FV_Atm(1)%npx, FV_Atm(1)%npy, FV_Atm(1)%npz, FV_Atm(1)%ncnst, FV_Atm(1)%ng,  &
                          real(conf%dt, fvprec), FV_Atm(1)%flagstruct%consv_te, FV_Atm(1)%flagstruct%fill,           &
                          FV_Atm(1)%flagstruct%reproduce_sum, real(kappa, fvprec),                                   &
                          real(cp, fvprec), real(zvir, fvprec), FV_Atm(1)%ptop, FV_Atm(1)%ks, FV_Atm(1)%flagstruct%ncnst,          &
                          FV_Atm(1)%flagstruct%n_split, FV_Atm(1)%flagstruct%q_split,                  &
                          FV_Atm(1)%u, FV_Atm(1)%v, FV_Atm(1)%w, FV_Atm(1)%delz,                       &
                          FV_Atm(1)%flagstruct%hydrostatic, FV_Atm(1)%pt, FV_Atm(1)%delp, FV_Atm(1)%q, &
                          FV_Atm(1)%ps, FV_Atm(1)%pe, FV_Atm(1)%pk, FV_Atm(1)%peln, FV_Atm(1)%pkz,     &
                          FV_Atm(1)%phis, FV_Atm(1)%q_con, FV_Atm(1)%omga,                             &
                          FV_Atm(1)%ua, FV_Atm(1)%va, FV_Atm(1)%uc, FV_Atm(1)%vc,                      &
                          FV_Atm(1)%ak, FV_Atm(1)%bk,                                                  &
                          FV_Atm(1)%mfx, FV_Atm(1)%mfy, FV_Atm(1)%cx, FV_Atm(1)%cy, FV_Atm(1)%ze0,     &
                          FV_Atm(1)%flagstruct%hybrid_z, FV_Atm(1)%gridstruct, FV_Atm(1)%flagstruct,   &
                          self%FV_AtmP(1)%flagstruct,                                                       &
                          FV_Atm(1)%neststruct, FV_Atm(1)%idiag, FV_Atm(1)%bd, FV_Atm(1)%parent_grid,  &
                          FV_Atm(1)%domain )
 endif


 ! MPP nulify
 ! ----------
 call nullify_domain()

 !Copy from fv3 back to traj structure
 !------------------------------------
 call fv3_to_traj(self,conf,traj)


endsubroutine step_nl

! ------------------------------------------------------------------------------

subroutine step_tl(self,conf,traj,pert)

 implicit none

 class(fv3jedi_lm_dynamics_type), target, intent(inout) :: self
 type(fv3jedi_lm_conf), intent(in)    :: conf
 type(fv3jedi_lm_traj), intent(in)    :: traj
 type(fv3jedi_lm_pert), intent(inout) :: pert

 type(fv_atmos_type), pointer :: FV_Atm(:)
 type(fv_atmos_pert_type), pointer :: FV_AtmP(:)
 integer :: i,j,k


 !Convenience pointer to the main FV_Atm structure
 !------------------------------------------------
 FV_Atm => self%FV_Atm
 FV_AtmP => self%FV_AtmP


 ! Set diagnostics to zeros
 ! ------------------------
 call zero_pert_vars(FV_AtmP(1))


 !Copy from traj/pert to the fv3 structures
 !-----------------------------------------
 call traj_to_fv3(self,conf,traj)
 call pert_to_fv3(self,conf,pert)


 !A-grid winds are diagnostic
 !---------------------------
 FV_AtmP(1)%uap = 0.0
 FV_AtmP(1)%vap = 0.0


 !Edge of pert always needs to be filled
 !--------------------------------------
 call mpp_get_boundary( FV_AtmP(1)%up, FV_AtmP(1)%vp, FV_Atm(1)%domain, &
                        wbuffery=self%wbuffery, ebuffery=self%ebuffery, &
                        sbufferx=self%sbufferx, nbufferx=self%nbufferx, &
                        gridtype=DGRID_NE, complete=.true. )
 do k=1,self%npz
    do i=self%isc,self%iec
       FV_AtmP(1)%up(i,self%jec+1,k) = self%nbufferx(i,k)
    enddo
 enddo
 do k=1,self%npz
    do j=self%jsc,self%jec
       FV_AtmP(1)%vp(self%iec+1,j,k) = self%ebuffery(j,k)
    enddo
 enddo


 !Compute the other pressure variables needed by FV3
 !--------------------------------------------------
 call compute_fv3_pressures_tlm( self%isc, self%iec, self%jsc, self%jec, &
                                 self%isd, self%ied, self%jsd, self%jed, &
                                 self%npz, real(kappa, fvprec), FV_Atm(1)%ptop, &
                                 FV_Atm(1)%delp, FV_AtmP(1)%delpp, &
                                 FV_Atm(1)%pe, FV_AtmP(1)%pep, &
                                 FV_Atm(1)%pk, FV_AtmP(1)%pkp, &
                                 FV_Atm(1)%pkz, FV_AtmP(1)%pkzp, &
                                 FV_Atm(1)%peln, FV_AtmP(1)%pelnp )


 ! MPP set domain
 ! --------------
 call set_domain(FV_Atm(1)%domain)


 !Propagate TLM one time step
 !---------------------------
 call fv_dynamics_tlm(FV_Atm(1)%npx, FV_Atm(1)%npy, FV_Atm(1)%npz, FV_Atm(1)%ncnst, FV_Atm(1)%ng,                      &
                      real(conf%DT, fvprec), FV_Atm(1)%flagstruct%consv_te, FV_Atm(1)%flagstruct%fill,                               &
                      FV_Atm(1)%flagstruct%reproduce_sum, real(kappa, fvprec),                                                       &
                      real(cp, fvprec), real(zvir, fvprec), FV_Atm(1)%ptop, FV_Atm(1)%ks, FV_Atm(1)%flagstruct%ncnst,                              &
                      FV_Atm(1)%flagstruct%n_split, FV_Atm(1)%flagstruct%q_split,                                      &
                      FV_Atm(1)%u, FV_AtmP(1)%up, FV_Atm(1)%v, FV_AtmP(1)%vp, FV_Atm(1)%w, FV_AtmP(1)%wp,                 &
                      FV_Atm(1)%delz, FV_AtmP(1)%delzp, FV_Atm(1)%flagstruct%hydrostatic,                               &
                      FV_Atm(1)%pt, FV_AtmP(1)%ptp, FV_Atm(1)%delp, FV_AtmP(1)%delpp,                                    &
                      FV_Atm(1)%q, FV_AtmP(1)%qp, FV_Atm(1)%ps, FV_AtmP(1)%psp, FV_Atm(1)%pe, FV_AtmP(1)%pep,             &
                      FV_Atm(1)%pk, FV_AtmP(1)%pkp, FV_Atm(1)%peln, FV_AtmP(1)%pelnp, FV_Atm(1)%pkz, FV_AtmP(1)%pkzp,     &
                      FV_Atm(1)%phis, FV_Atm(1)%q_con, FV_Atm(1)%omga, FV_AtmP(1)%omgap,                                &
                      FV_Atm(1)%ua, FV_AtmP(1)%uap, FV_Atm(1)%va, FV_AtmP(1)%vap,                                        &
                      FV_Atm(1)%uc, FV_AtmP(1)%ucp, FV_Atm(1)%vc, FV_AtmP(1)%vcp,                                        &
                      FV_Atm(1)%ak, FV_Atm(1)%bk,                                                                      &
                      FV_Atm(1)%mfx, FV_AtmP(1)%mfxp, FV_Atm(1)%mfy, FV_AtmP(1)%mfyp,                                    &
                      FV_Atm(1)%cx, FV_AtmP(1)%cxp, FV_Atm(1)%cy,  FV_AtmP(1)%cyp, FV_Atm(1)%ze0,                       &
                      FV_Atm(1)%flagstruct%hybrid_z, FV_Atm(1)%gridstruct, FV_Atm(1)%flagstruct, FV_AtmP(1)%flagstruct, &
                      FV_Atm(1)%neststruct, FV_Atm(1)%idiag, FV_Atm(1)%bd, FV_Atm(1)%parent_grid,FV_Atm(1)%domain      )


 ! MPP nulify
 ! ----------
 call nullify_domain()


 !Copy from fv3 back to pert structure
 !------------------------------------
 call fv3_to_pert(self,conf,pert)


 ! Set diagnostics to zeros
 ! ------------------------
 call zero_pert_vars(FV_AtmP(1))


endsubroutine step_tl

! ------------------------------------------------------------------------------

subroutine step_ad(self,conf,traj,pert)

 implicit none

 class(fv3jedi_lm_dynamics_type), target, intent(inout) :: self
 type(fv3jedi_lm_conf), intent(in)    :: conf
 type(fv3jedi_lm_traj), intent(in)    :: traj
 type(fv3jedi_lm_pert), intent(inout) :: pert

 type(fv_atmos_type), pointer :: FV_Atm(:)
 type(fv_atmos_pert_type), pointer :: FV_AtmP(:)
 integer :: i,j,k


 !Convenience pointer to the main FV_Atm structure
 !------------------------------------------------
 FV_Atm  => self%FV_Atm
 FV_AtmP => self%FV_AtmP


 ! Set diagnostics to zeros
 ! ------------------------
 call zero_pert_vars(FV_AtmP(1))


 !Copy from traj/pert to the fv3 structures
 !-----------------------------------------
 call traj_to_fv3(self,conf,traj)
 call pert_to_fv3(self,conf,pert)


 ! MPP set domain
 ! --------------
 call set_domain(FV_Atm(1)%domain)


 ! Initilize the module level checkpointing
 ! ----------------------------------------
 if (cp_iter_controls%cp_i .ne. 0) then
    call cp_mod_ini(self%cp_dyn_ind)
 endif


 ! Forward sweep of the dynamics with saving of checkpoints for use in backward sweep
 ! ----------------------------------------------------------------------------------
 if (cp_iter_controls%cp_i <= 3) then

    call fv_dynamics_fwd(FV_Atm(1)%npx, FV_Atm(1)%npy, FV_Atm(1)%npz, FV_Atm(1)%ncnst, FV_Atm(1)%ng,                      &
                         real(conf%DT, fvprec), FV_Atm(1)%flagstruct%consv_te, FV_Atm(1)%flagstruct%fill,                               &
                         FV_Atm(1)%flagstruct%reproduce_sum, real(kappa, fvprec),                                                       &
                         real(cp, fvprec), real(zvir, fvprec), FV_Atm(1)%ptop, FV_Atm(1)%ks, FV_Atm(1)%flagstruct%ncnst,                              &
                         FV_Atm(1)%flagstruct%n_split, FV_Atm(1)%flagstruct%q_split,                                      &
                         FV_Atm(1)%u, FV_Atm(1)%v, FV_Atm(1)%w,                                                           &
                         FV_Atm(1)%delz, FV_Atm(1)%flagstruct%hydrostatic,                                                &
                         FV_Atm(1)%pt, FV_Atm(1)%delp,                                                                    &
                         FV_Atm(1)%q, FV_Atm(1)%ps, FV_Atm(1)%pe,                                                         &
                         FV_Atm(1)%pk, FV_Atm(1)%peln, FV_Atm(1)%pkz,                                                     &
                         FV_Atm(1)%phis, FV_Atm(1)%q_con, FV_Atm(1)%omga,                                                 &
                         FV_Atm(1)%ua, FV_Atm(1)%va,                                                                      &
                         FV_Atm(1)%uc, FV_Atm(1)%vc,                                                                      &
                         FV_Atm(1)%ak, FV_Atm(1)%bk,                                                                      &
                         FV_Atm(1)%mfx, FV_Atm(1)%mfy,                                                                    &
                         FV_Atm(1)%cx, FV_Atm(1)%cy, FV_Atm(1)%ze0,                                                       &
                         FV_Atm(1)%flagstruct%hybrid_z, FV_Atm(1)%gridstruct, FV_Atm(1)%flagstruct,                       &
                         FV_AtmP(1)%flagstruct,                                                                           &
                         FV_Atm(1)%neststruct, FV_Atm(1)%idiag, FV_Atm(1)%bd, FV_Atm(1)%parent_grid,FV_Atm(1)%domain      )

    if (cp_iter_controls%cp_i .ne. 0) then
       !Push end of timestep trajectory to stack
       call PUSHREALARRAY(FV_Atm(1)%u   ,(self%ied-self%isd+1)*(self%jed-self%jsd+2)*self%npz)
       call PUSHREALARRAY(FV_Atm(1)%v   ,(self%ied-self%isd+2)*(self%jed-self%jsd+1)*self%npz)
       call PUSHREALARRAY(FV_Atm(1)%w   ,(self%ied-self%isd+1)*(self%jed-self%jsd+1)*self%npz)
       call PUSHREALARRAY(FV_Atm(1)%delz,(self%ied-self%isd+1)*(self%jed-self%jsd+1)*self%npz)
       call PUSHREALARRAY(FV_Atm(1)%pt  ,(self%ied-self%isd+1)*(self%jed-self%jsd+1)*self%npz)
       call PUSHREALARRAY(FV_Atm(1)%delp,(self%ied-self%isd+1)*(self%jed-self%jsd+1)*self%npz)
       call PUSHREALARRAY(FV_Atm(1)%q   ,(self%ied-self%isd+1)*(self%jed-self%jsd+1)*self%npz*FV_Atm(1)%ncnst)
       call PUSHREALARRAY(FV_Atm(1)%ps  ,(self%ied-self%isd+1)*(self%jed-self%jsd+1))
       call PUSHREALARRAY(FV_Atm(1)%pe  ,(self%iec-self%isc+3)*(self%jec-self%jsc+3)*(self%npz+1))
       call PUSHREALARRAY(FV_Atm(1)%pk  ,(self%iec-self%isc+1)*(self%jec-self%jsc+1)*(self%npz+1))
       call PUSHREALARRAY(FV_Atm(1)%peln,(self%iec-self%isc+1)*(self%jec-self%jsc+1)*(self%npz+1))
       call PUSHREALARRAY(FV_Atm(1)%pkz ,(self%iec-self%isc+1)*(self%jec-self%jsc+1)*self%npz)
       call PUSHREALARRAY(FV_Atm(1)%phis,(self%ied-self%isd+1)*(self%jed-self%jsd+1))
       call PUSHREALARRAY(FV_Atm(1)%omga,(self%ied-self%isd+1)*(self%jed-self%jsd+1)*self%npz)
       call PUSHREALARRAY(FV_Atm(1)%ua  ,(self%ied-self%isd+1)*(self%jed-self%jsd+1)*self%npz)
       call PUSHREALARRAY(FV_Atm(1)%va  ,(self%ied-self%isd+1)*(self%jed-self%jsd+1)*self%npz)
       call PUSHREALARRAY(FV_Atm(1)%uc  ,(self%ied-self%isd+2)*(self%jed-self%jsd+1)*self%npz)
       call PUSHREALARRAY(FV_Atm(1)%vc  ,(self%ied-self%isd+1)*(self%jed-self%jsd+2)*self%npz)
       call PUSHREALARRAY(FV_Atm(1)%mfx ,(self%iec-self%isc+2)*(self%jec-self%jsc+1)*self%npz)
       call PUSHREALARRAY(FV_Atm(1)%mfy ,(self%iec-self%isc+1)*(self%jec-self%jsc+2)*self%npz)
       call PUSHREALARRAY(FV_Atm(1)%cx  ,(self%iec-self%isc+2)*(self%jed-self%jsd+1)*self%npz)
       call PUSHREALARRAY(FV_Atm(1)%cy  ,(self%ied-self%isd+1)*(self%jec-self%jsc+2)*self%npz)
       !Trick checkpoint schemes into not considering these superfluous checkpoints,
       !about to recover with the pop anyway.
       FV_Atm(1)%u    = 2.0_kind_real*FV_Atm(1)%u
       FV_Atm(1)%v    = 2.0_kind_real*FV_Atm(1)%v
       FV_Atm(1)%w    = 2.0_kind_real*FV_Atm(1)%w
       FV_Atm(1)%delz = 2.0_kind_real*FV_Atm(1)%delz
       FV_Atm(1)%pt   = 2.0_kind_real*FV_Atm(1)%pt
       FV_Atm(1)%delp = 2.0_kind_real*FV_Atm(1)%delp
       FV_Atm(1)%q    = 2.0_kind_real*FV_Atm(1)%q
       FV_Atm(1)%ps   = 2.0_kind_real*FV_Atm(1)%ps
       FV_Atm(1)%pe   = 2.0_kind_real*FV_Atm(1)%pe
       FV_Atm(1)%pk   = 2.0_kind_real*FV_Atm(1)%pk
       FV_Atm(1)%peln = 2.0_kind_real*FV_Atm(1)%peln
       FV_Atm(1)%pkz  = 2.0_kind_real*FV_Atm(1)%pkz
       FV_Atm(1)%phis = 2.0_kind_real*FV_Atm(1)%phis
       FV_Atm(1)%omga = 2.0_kind_real*FV_Atm(1)%omga
       FV_Atm(1)%ua   = 2.0_kind_real*FV_Atm(1)%ua
       FV_Atm(1)%va   = 2.0_kind_real*FV_Atm(1)%va
       FV_Atm(1)%uc   = 2.0_kind_real*FV_Atm(1)%uc
       FV_Atm(1)%vc   = 2.0_kind_real*FV_Atm(1)%vc
       FV_Atm(1)%mfx  = 2.0_kind_real*FV_Atm(1)%mfx
       FV_Atm(1)%mfy  = 2.0_kind_real*FV_Atm(1)%mfy
       FV_Atm(1)%cx   = 2.0_kind_real*FV_Atm(1)%cx
       FV_Atm(1)%cy   = 2.0_kind_real*FV_Atm(1)%cy
    endif

 endif


 ! Checkpoint mid point, reset counters etc
 ! ----------------------------------------
 if (cp_iter_controls%cp_i .ne. 0) then
    call cp_mod_mid
 endif

 if (cp_iter_controls%cp_i .ne. 0) then
    !Populate end of timestep trajectory from stack
    call POPREALARRAY(FV_Atm(1)%cy  ,(self%ied-self%isd+1)*(self%jec-self%jsc+2)*self%npz)
    call POPREALARRAY(FV_Atm(1)%cx  ,(self%iec-self%isc+2)*(self%jed-self%jsd+1)*self%npz)
    call POPREALARRAY(FV_Atm(1)%mfy ,(self%iec-self%isc+1)*(self%jec-self%jsc+2)*self%npz)
    call POPREALARRAY(FV_Atm(1)%mfx ,(self%iec-self%isc+2)*(self%jec-self%jsc+1)*self%npz)
    call POPREALARRAY(FV_Atm(1)%vc  ,(self%ied-self%isd+1)*(self%jed-self%jsd+2)*self%npz)
    call POPREALARRAY(FV_Atm(1)%uc  ,(self%ied-self%isd+2)*(self%jed-self%jsd+1)*self%npz)
    call POPREALARRAY(FV_Atm(1)%va  ,(self%ied-self%isd+1)*(self%jed-self%jsd+1)*self%npz)
    call POPREALARRAY(FV_Atm(1)%ua  ,(self%ied-self%isd+1)*(self%jed-self%jsd+1)*self%npz)
    call POPREALARRAY(FV_Atm(1)%omga,(self%ied-self%isd+1)*(self%jed-self%jsd+1)*self%npz)
    call POPREALARRAY(FV_Atm(1)%phis,(self%ied-self%isd+1)*(self%jed-self%jsd+1))
    call POPREALARRAY(FV_Atm(1)%pkz ,(self%iec-self%isc+1)*(self%jec-self%jsc+1)*self%npz)
    call POPREALARRAY(FV_Atm(1)%peln,(self%iec-self%isc+1)*(self%jec-self%jsc+1)*(self%npz+1))
    call POPREALARRAY(FV_Atm(1)%pk  ,(self%iec-self%isc+1)*(self%jec-self%jsc+1)*(self%npz+1))
    call POPREALARRAY(FV_Atm(1)%pe  ,(self%iec-self%isc+3)*(self%jec-self%jsc+3)*(self%npz+1))
    call POPREALARRAY(FV_Atm(1)%ps  ,(self%ied-self%isd+1)*(self%jed-self%jsd+1))
    call POPREALARRAY(FV_Atm(1)%q   ,(self%ied-self%isd+1)*(self%jed-self%jsd+1)*self%npz*FV_Atm(1)%ncnst)
    call POPREALARRAY(FV_Atm(1)%delp,(self%ied-self%isd+1)*(self%jed-self%jsd+1)*self%npz)
    call POPREALARRAY(FV_Atm(1)%pt  ,(self%ied-self%isd+1)*(self%jed-self%jsd+1)*self%npz)
    call POPREALARRAY(FV_Atm(1)%delz,(self%ied-self%isd+1)*(self%jed-self%jsd+1)*self%npz)
    call POPREALARRAY(FV_Atm(1)%w   ,(self%ied-self%isd+1)*(self%jed-self%jsd+1)*self%npz)
    call POPREALARRAY(FV_Atm(1)%v   ,(self%ied-self%isd+2)*(self%jed-self%jsd+1)*self%npz)
    call POPREALARRAY(FV_Atm(1)%u   ,(self%ied-self%isd+1)*(self%jed-self%jsd+2)*self%npz)
 endif


 ! Backward adjoint sweep of the dynamics
 ! --------------------------------------
 call fv_dynamics_bwd(FV_Atm(1)%npx, FV_Atm(1)%npy, FV_Atm(1)%npz, FV_Atm(1)%ncnst, FV_Atm(1)%ng,                      &
                      real(conf%DT, fvprec), FV_Atm(1)%flagstruct%consv_te, FV_Atm(1)%flagstruct%fill,                               &
                      FV_Atm(1)%flagstruct%reproduce_sum, real(kappa, fvprec),                                                       &
                      real(cp, fvprec), real(zvir, fvprec), FV_Atm(1)%ptop, FV_Atm(1)%ks, FV_Atm(1)%flagstruct%ncnst,                              &
                      FV_Atm(1)%flagstruct%n_split, FV_Atm(1)%flagstruct%q_split,                                      &
                      FV_Atm(1)%u, FV_AtmP(1)%up, FV_Atm(1)%v, FV_AtmP(1)%vp, FV_Atm(1)%w, FV_AtmP(1)%wp,                 &
                      FV_Atm(1)%delz, FV_AtmP(1)%delzp, FV_Atm(1)%flagstruct%hydrostatic,                               &
                      FV_Atm(1)%pt, FV_AtmP(1)%ptp, FV_Atm(1)%delp, FV_AtmP(1)%delpp,                                    &
                      FV_Atm(1)%q, FV_AtmP(1)%qp, FV_Atm(1)%ps, FV_AtmP(1)%psp, FV_Atm(1)%pe, FV_AtmP(1)%pep,             &
                      FV_Atm(1)%pk, FV_AtmP(1)%pkp, FV_Atm(1)%peln, FV_AtmP(1)%pelnp, FV_Atm(1)%pkz, FV_AtmP(1)%pkzp,     &
                      FV_Atm(1)%phis, FV_Atm(1)%q_con, FV_Atm(1)%omga, FV_AtmP(1)%omgap,                                &
                      FV_Atm(1)%ua, FV_AtmP(1)%uap, FV_Atm(1)%va, FV_AtmP(1)%vap,                                        &
                      FV_Atm(1)%uc, FV_AtmP(1)%ucp, FV_Atm(1)%vc, FV_AtmP(1)%vcp,                                        &
                      FV_Atm(1)%ak, FV_Atm(1)%bk,                                                                      &
                      FV_Atm(1)%mfx, FV_AtmP(1)%mfxp, FV_Atm(1)%mfy, FV_AtmP(1)%mfyp,                                    &
                      FV_Atm(1)%cx, FV_AtmP(1)%cxp, FV_Atm(1)%cy,  FV_AtmP(1)%cyp, FV_Atm(1)%ze0,                      &
                      FV_Atm(1)%flagstruct%hybrid_z, FV_Atm(1)%gridstruct, FV_Atm(1)%flagstruct,                       &
                      FV_AtmP(1)%flagstruct,                                                                           &
                      FV_Atm(1)%neststruct, FV_Atm(1)%idiag, FV_Atm(1)%bd, FV_Atm(1)%parent_grid,FV_Atm(1)%domain      )


 !Adjoint of compute the other pressure variables needed by FV3
 !-------------------------------------------------------------
 call compute_fv3_pressures_bwd( self%isc, self%iec, self%jsc, self%jec, &
                                 self%isd, self%ied, self%jsd, self%jed, &
                                 self%npz, real(kappa, fvprec), FV_Atm(1)%ptop, &
                                 FV_Atm(1)%delp, FV_AtmP(1)%delpp, &
                                 FV_Atm(1)%pe, FV_AtmP(1)%pep, &
                                 FV_Atm(1)%pk, FV_AtmP(1)%pkp, &
                                 FV_Atm(1)%pkz, FV_AtmP(1)%pkzp, &
                                 FV_Atm(1)%peln, FV_AtmP(1)%pelnp )


 !Edge of pert always needs to be filled
 !--------------------------------------
 self%nbufferx = 0.0_kind_real
 do k=1,self%npz
    do i=self%isc,self%iec
       self%nbufferx(i,k) = FV_AtmP(1)%up(i,self%jec+1,k)
    enddo
 enddo
 self%ebuffery = 0.0_kind_real
 do k=1,self%npz
    do j=self%jsc,self%jec
       self%ebuffery(j,k) = FV_AtmP(1)%vp(self%iec+1,j,k)
    enddo
 enddo

 call mpp_get_boundary_ad( FV_AtmP(1)%up, FV_AtmP(1)%vp, FV_Atm(1)%domain, &
                           wbuffery=self%wbuffery, ebuffery=self%ebuffery, sbufferx=self%sbufferx, nbufferx=self%nbufferx, &
                           gridtype=DGRID_NE, complete=.true. )


 ! MPP nulify
 ! ----------
 call nullify_domain()


 !A-grid winds are diagnostic
 !---------------------------
 FV_AtmP(1)%uap = 0.0
 FV_AtmP(1)%vap = 0.0


 !Copy from fv3 back to pert structure
 !------------------------------------
 call fv3_to_pert(self,conf,pert)


 ! Set diagnostics to zeros
 ! ------------------------
 call zero_pert_vars(FV_AtmP(1))


endsubroutine step_ad

! ------------------------------------------------------------------------------

subroutine delete(self,conf)

 implicit none

 class(fv3jedi_lm_dynamics_type), intent(inout) :: self
 type(fv3jedi_lm_conf), intent(in) :: conf

 deallocate(self%ebuffery)
 deallocate(self%wbuffery)
 deallocate(self%nbufferx)
 deallocate(self%sbufferx)

 call deallocate_fv_atmos_type(self%FV_Atm(1))
 deallocate(self%FV_Atm)

 call deallocate_fv_atmos_pert_type(self%FV_AtmP(1))
 deallocate(self%FV_AtmP)

 if (cp_iter_controls%cp_i .ne. 0) call finalize_cp_iter

endsubroutine delete

! ------------------------------------------------------------------------------

subroutine traj_to_fv3(self,conf,traj)

 implicit none

 class(fv3jedi_lm_dynamics_type), intent(inout) :: self
 type(fv3jedi_lm_conf), intent(in) :: conf
 type(fv3jedi_lm_traj), intent(in) :: traj

 integer :: i,j,k

 !Zero the halos
 !--------------
 self%FV_Atm(1)%u     = 0.0_kind_real
 self%FV_Atm(1)%v     = 0.0_kind_real
 self%FV_Atm(1)%pt    = 0.0_kind_real
 self%FV_Atm(1)%delp  = 0.0_kind_real
 self%FV_Atm(1)%q     = 0.0_kind_real
 self%FV_Atm(1)%w     = 0.0_kind_real
 self%FV_Atm(1)%delz  = 0.0_kind_real
 self%FV_Atm(1)%phis  = 0.0_kind_real
 self%FV_Atm(1)%pe    = 0.0_kind_real
 self%FV_Atm(1)%peln  = 0.0_kind_real
 self%FV_Atm(1)%pk    = 0.0_kind_real
 self%FV_Atm(1)%pkz   = 0.0_kind_real
 self%FV_Atm(1)%ua    = 0.0_kind_real
 self%FV_Atm(1)%va    = 0.0_kind_real
 self%FV_Atm(1)%uc    = 0.0_kind_real
 self%FV_Atm(1)%vc    = 0.0_kind_real
 self%FV_Atm(1)%omga  = 0.0_kind_real
 self%FV_Atm(1)%mfx   = 0.0_kind_real
 self%FV_Atm(1)%mfy   = 0.0_kind_real
 self%FV_Atm(1)%cx    = 0.0_kind_real
 self%FV_Atm(1)%cy    = 0.0_kind_real
 self%FV_Atm(1)%ze0   = 0.0_kind_real
 self%FV_Atm(1)%q_con = 0.0_kind_real
 self%FV_Atm(1)%ps    = 0.0_kind_real


 !Copy from traj
 !--------------
 self%FV_Atm(1)%u   (self%isc:self%iec,self%jsc:self%jec,:) = traj%u   (self%isc:self%iec,self%jsc:self%jec,:)
 self%FV_Atm(1)%v   (self%isc:self%iec,self%jsc:self%jec,:) = traj%v   (self%isc:self%iec,self%jsc:self%jec,:)
 self%FV_Atm(1)%pt  (self%isc:self%iec,self%jsc:self%jec,:) = traj%t   (self%isc:self%iec,self%jsc:self%jec,:)
 self%FV_Atm(1)%delp(self%isc:self%iec,self%jsc:self%jec,:) = traj%delp(self%isc:self%iec,self%jsc:self%jec,:)
 self%FV_Atm(1)%q(self%isc:self%iec,self%jsc:self%jec,:,:) = traj%tracers(self%isc:self%iec,self%jsc:self%jec,:,:)

 if (.not. self%FV_Atm(1)%flagstruct%hydrostatic) then
    self%FV_Atm(1)%delz(self%isc:self%iec  ,self%jsc:self%jec  ,:  ) = traj%delz(self%isc:self%iec  ,self%jsc:self%jec  ,:  )
    self%FV_Atm(1)%w   (self%isc:self%iec  ,self%jsc:self%jec  ,:  ) = traj%w   (self%isc:self%iec  ,self%jsc:self%jec  ,:  )
 endif

 self%FV_Atm(1)%phis(self%isc:self%iec,self%jsc:self%jec) = traj%phis(self%isc:self%iec,self%jsc:self%jec)

 !Update edges of d-grid winds
 !----------------------------
 call mpp_get_boundary(self%FV_Atm(1)%u, self%FV_Atm(1)%v, self%FV_Atm(1)%domain, &
                       wbuffery=self%wbuffery, ebuffery=self%ebuffery, &
                       sbufferx=self%sbufferx, nbufferx=self%nbufferx, &
                       gridtype=DGRID_NE, complete=.true. )
 do k=1,self%npz
    do i=self%isc,self%iec
       self%FV_Atm(1)%u(i,self%jec+1,k) = self%nbufferx(i,k)
    enddo
 enddo
 do k=1,self%npz
    do j=self%jsc,self%jec
       self%FV_Atm(1)%v(self%iec+1,j,k) = self%ebuffery(j,k)
    enddo
 enddo


 ! Fill phi halos
 ! --------------
 call mpp_update_domains(self%FV_Atm(1)%phis, self%FV_Atm(1)%domain, complete=.true.)


 !Compute the other pressure variables needed by FV3
 !--------------------------------------------------
 call compute_fv3_pressures( self%isc, self%iec, self%jsc, self%jec, self%isd, self%ied, self%jsd, self%jed, &
                             self%npz, real(kappa, fvprec), self%FV_Atm(1)%ptop, &
                             self%FV_Atm(1)%delp, self%FV_Atm(1)%pe, self%FV_Atm(1)%pk, self%FV_Atm(1)%pkz, self%FV_Atm(1)%peln )

endsubroutine traj_to_fv3

! ------------------------------------------------------------------------------

subroutine fv3_to_traj(self,conf,traj)

 implicit none

 class(fv3jedi_lm_dynamics_type), intent(in) :: self
 type(fv3jedi_lm_conf), intent(in)    :: conf
 type(fv3jedi_lm_traj), intent(inout) :: traj

 traj%u   (self%isc:self%iec,self%jsc:self%jec,:) = self%FV_Atm(1)%u   (self%isc:self%iec,self%jsc:self%jec,:)
 traj%v   (self%isc:self%iec,self%jsc:self%jec,:) = self%FV_Atm(1)%v   (self%isc:self%iec,self%jsc:self%jec,:)
 traj%t   (self%isc:self%iec,self%jsc:self%jec,:) = self%FV_Atm(1)%pt  (self%isc:self%iec,self%jsc:self%jec,:)
 traj%delp(self%isc:self%iec,self%jsc:self%jec,:) = self%FV_Atm(1)%delp(self%isc:self%iec,self%jsc:self%jec,:)
 traj%tracers(self%isc:self%iec,self%jsc:self%jec,:,:) = self%FV_Atm(1)%q(self%isc:self%iec,self%jsc:self%jec,:,:)

 if (.not. self%FV_Atm(1)%flagstruct%hydrostatic) then
    traj%delz(self%isc:self%iec,self%jsc:self%jec,:) = self%FV_Atm(1)%delz(self%isc:self%iec,self%jsc:self%jec,:)
    traj%w   (self%isc:self%iec,self%jsc:self%jec,:) = self%FV_Atm(1)%w   (self%isc:self%iec,self%jsc:self%jec,:)
 endif

 traj%ua(self%isc:self%iec,self%jsc:self%jec,:) = self%FV_Atm(1)%ua(self%isc:self%iec,self%jsc:self%jec,:)
 traj%va(self%isc:self%iec,self%jsc:self%jec,:) = self%FV_Atm(1)%va(self%isc:self%iec,self%jsc:self%jec,:)

 traj%phis(self%isc:self%iec,self%jsc:self%jec) = self%FV_Atm(1)%phis(self%isc:self%iec,self%jsc:self%jec)

endsubroutine fv3_to_traj

! ------------------------------------------------------------------------------

subroutine pert_to_fv3(self,conf,pert)

 implicit none

 class(fv3jedi_lm_dynamics_type), intent(inout) :: self
 type(fv3jedi_lm_conf), intent(in) :: conf
 type(fv3jedi_lm_pert), intent(in) :: pert

 !To zero the halos
 self%FV_AtmP(1)%up    = 0.0
 self%FV_AtmP(1)%vp    = 0.0
 self%FV_AtmP(1)%ptp   = 0.0
 self%FV_AtmP(1)%delpp = 0.0
 self%FV_AtmP(1)%qp    = 0.0
 self%FV_AtmP(1)%wp    = 0.0
 self%FV_AtmP(1)%delzp = 0.0
 self%FV_AtmP(1)%uap   = 0.0
 self%FV_AtmP(1)%vap   = 0.0

 self%FV_AtmP(1)%up   (self%isc:self%iec,self%jsc:self%jec,:) = pert%u   (self%isc:self%iec,self%jsc:self%jec,:)
 self%FV_AtmP(1)%vp   (self%isc:self%iec,self%jsc:self%jec,:) = pert%v   (self%isc:self%iec,self%jsc:self%jec,:)
 self%FV_AtmP(1)%ptp  (self%isc:self%iec,self%jsc:self%jec,:) = pert%T   (self%isc:self%iec,self%jsc:self%jec,:)
 self%FV_AtmP(1)%delpp(self%isc:self%iec,self%jsc:self%jec,:) = pert%delp(self%isc:self%iec,self%jsc:self%jec,:)
 self%FV_AtmP(1)%qp(self%isc:self%iec,self%jsc:self%jec,:,:) = pert%tracers(self%isc:self%iec,self%jsc:self%jec,:,:)

 if (.not. self%FV_Atm(1)%flagstruct%hydrostatic) then
    self%FV_AtmP(1)%delzp(self%isc:self%iec,self%jsc:self%jec,:) = pert%delz(self%isc:self%iec,self%jsc:self%jec,:)
    self%FV_AtmP(1)%wp   (self%isc:self%iec,self%jsc:self%jec,:) = pert%w   (self%isc:self%iec,self%jsc:self%jec,:)
 endif

 self%FV_AtmP(1)%uap(self%isc:self%iec,self%jsc:self%jec,:) = 0.0!pert%ua(self%isc:self%iec,self%jsc:self%jec,:)
 self%FV_AtmP(1)%vap(self%isc:self%iec,self%jsc:self%jec,:) = 0.0!pert%va(self%isc:self%iec,self%jsc:self%jec,:)

endsubroutine pert_to_fv3

! ------------------------------------------------------------------------------

subroutine fv3_to_pert(self,conf,pert)

 implicit none

 class(fv3jedi_lm_dynamics_type), intent(inout) :: self
 type(fv3jedi_lm_conf), intent(in)    :: conf
 type(fv3jedi_lm_pert), intent(inout) :: pert

 pert%u   (self%isc:self%iec,self%jsc:self%jec,:) = self%FV_AtmP(1)%up   (self%isc:self%iec,self%jsc:self%jec,:)
 pert%v   (self%isc:self%iec,self%jsc:self%jec,:) = self%FV_AtmP(1)%vp   (self%isc:self%iec,self%jsc:self%jec,:)
 pert%T   (self%isc:self%iec,self%jsc:self%jec,:) = self%FV_AtmP(1)%ptp  (self%isc:self%iec,self%jsc:self%jec,:)
 pert%delp(self%isc:self%iec,self%jsc:self%jec,:) = self%FV_AtmP(1)%delpp(self%isc:self%iec,self%jsc:self%jec,:)
 pert%tracers(self%isc:self%iec,self%jsc:self%jec,:,:) = self%FV_AtmP(1)%qp(self%isc:self%iec,self%jsc:self%jec,:,:)

 if (.not. self%FV_Atm(1)%flagstruct%hydrostatic) then
   pert%delz(self%isc:self%iec,self%jsc:self%jec,:) = self%FV_AtmP(1)%delzp(self%isc:self%iec,self%jsc:self%jec,:)
   pert%w   (self%isc:self%iec,self%jsc:self%jec,:) = self%FV_AtmP(1)%wp   (self%isc:self%iec,self%jsc:self%jec,:)
 endif

 pert%ua(self%isc:self%iec,self%jsc:self%jec,:) = 0.0!self%FV_AtmP(1)%uap(self%isc:self%iec,self%jsc:self%jec,:)
 pert%va(self%isc:self%iec,self%jsc:self%jec,:) = 0.0!self%FV_AtmP(1)%vap(self%isc:self%iec,self%jsc:self%jec,:)

 self%FV_AtmP(1)%up    = 0.0
 self%FV_AtmP(1)%vp    = 0.0
 self%FV_AtmP(1)%ptp   = 0.0
 self%FV_AtmP(1)%delpp = 0.0
 self%FV_AtmP(1)%qp    = 0.0
 self%FV_AtmP(1)%wp    = 0.0
 self%FV_AtmP(1)%delzp = 0.0
 self%FV_AtmP(1)%uap   = 0.0
 self%FV_AtmP(1)%vap   = 0.0

endsubroutine fv3_to_pert

! ------------------------------------------------------------------------------

subroutine zero_pert_vars(FV_AtmP)

implicit none
type(fv_atmos_pert_type), intent(inout) :: FV_AtmP

!Prognostic
FV_AtmP%up = 0.0
FV_AtmP%vp = 0.0
FV_AtmP%ptp = 0.0
FV_AtmP%delpp = 0.0
FV_AtmP%qp = 0.0
FV_AtmP%wp = 0.0
FV_AtmP%delzP = 0.0
FV_AtmP%qp = 0.0

!Outputs
FV_AtmP%ze0p = 0.0
FV_AtmP%q_conp = 0.0
FV_AtmP%psp = 0.0
FV_AtmP%pep = 0.0
FV_AtmP%pkp = 0.0
FV_AtmP%pelnp = 0.0
FV_AtmP%pkzp = 0.0
FV_AtmP%omgap = 0.0
FV_AtmP%uap = 0.0
FV_AtmP%vap = 0.0
FV_AtmP%ucp = 0.0
FV_AtmP%vcp = 0.0
FV_AtmP%mfxp = 0.0
FV_AtmP%mfyp = 0.0
FV_AtmP%cxp = 0.0
FV_AtmP%cyp = 0.0

end subroutine zero_pert_vars

! ------------------------------------------------------------------------------

subroutine check_tracers(traj, pert)

 type(fv3jedi_lm_traj), intent(in) :: traj
 type(fv3jedi_lm_pert), intent(in) :: pert

 integer :: i

 ! Check that tracers are allocated
 if (.not. allocated(traj%tracers)) &
   call mpp_error(FATAL, 'check_tracers: .not. allocated(traj%tracers)')

 if (.not. allocated(pert%tracers)) &
   call mpp_error(FATAL, 'check_tracers: .not. allocated(pert%tracers)')

 ! Check that size of traj%tracers is the same as pert%tracers
 if (size(traj%tracers,4) /= size(pert%tracers,4)) then
   if (mpp_pe() == mpp_root_pe()) then
     print*, 'Size of traj%tracers = ', size(traj%tracers,4), &
             'Size of pert%tracers = ', size(pert%tracers,4)
     do i = 1, size(traj%tracers,4)
       print*, 'Traj tracer at index ', i, ' = ', trim(traj%tracer_names(i))
     enddo
     do i = 1, size(pert%tracers,4)
       print*, 'Pert tracer at index ', i, ' = ', trim(pert%tracer_names(i))
     enddo
   endif
   call mpp_error(FATAL, 'check_tracers: traj and pert tracer arrays have different size')
 endif

 ! Loop over traj%tracer_names and pert%tracer_names and check they are all the same
 do i = 1, size(traj%tracer_names)
   if (trim(traj%tracer_names(i)) /= trim(pert%tracer_names(i))) then
     if (mpp_pe() == mpp_root_pe()) then
       print*, 'Traj tracer with index ', i, ' has name ', trim(traj%tracer_names(i))
       print*, 'Pert tracer with index ', i, ' has name ', trim(pert%tracer_names(i))
     endif
     call mpp_error(FATAL, 'check_tracers: traj and pert tracers have different names in the same position')
   endif
 enddo

endsubroutine check_tracers

! ------------------------------------------------------------------------------

subroutine allocate_tracers(FV_Atm)

  type(fv_atmos_type), intent(inout), target :: FV_Atm(:)

  ! If the size is correct then return
  if (size(FV_Atm(1)%q,4) == FV_Atm(1)%ncnst) return

  ! If internal tracers are allocated then deallocate
  if (allocated(FV_Atm(1)%q)) deallocate(FV_Atm(1)%q)

  ! Allocate internal tracers
  allocate(FV_Atm(1)%q(FV_Atm(1)%bd%isd:FV_Atm(1)%bd%ied, &
                       FV_Atm(1)%bd%jsd:FV_Atm(1)%bd%jed, &
                       FV_Atm(1)%flagstruct%npz, FV_Atm(1)%ncnst))

endsubroutine allocate_tracers

! ------------------------------------------------------------------------------
subroutine allocate_tracersP(FV_AtmP, FV_Atm)

  type(fv_atmos_pert_type), intent(inout), target :: FV_AtmP(:)
  type(fv_atmos_type), intent(in), target :: FV_Atm(:)

  ! If the size is correct then return
  if (size(FV_AtmP(1)%qp,4) == FV_Atm(1)%ncnst) return

  ! If internal perturabtion tracers are allocated then deallocate
  if (allocated(FV_AtmP(1)%qp)) deallocate(FV_AtmP(1)%qp)

  ! Allocate internal perturbation tracers
  allocate(FV_AtmP(1)%qp(FV_Atm(1)%bd%isd:FV_Atm(1)%bd%ied, &
                         FV_Atm(1)%bd%jsd:FV_Atm(1)%bd%jed, &
                         FV_Atm(1)%flagstruct%npz, &
                         FV_Atm(1)%ncnst))

endsubroutine allocate_tracersP

! ------------------------------------------------------------------------------

subroutine a_to_d_nl(self, isc, iec, jsc, jec, npz, ua, va, ud, vd)

! Args
class(fv3jedi_lm_dynamics_type), target, intent(inout) :: self
integer, intent(in) :: isc, iec, jsc, jec, npz
real(kind=kind_real), intent(in)  :: ua(isc:iec, jsc:jec, npz)
real(kind=kind_real), intent(in)  :: va(isc:iec, jsc:jec, npz)
real(kind=kind_real), intent(out) :: ud(isc:iec, jsc:jec, npz)
real(kind=kind_real), intent(out) :: vd(isc:iec, jsc:jec, npz)

! Locals
type(fv_atmos_type), pointer :: Atm(:)

integer :: npx, npy
integer :: i,j,k, im2,jm2

real(kind_real), allocatable, dimension(:,:,:)   :: uatemp, vatemp

real(kind=kind_real) :: v3(isc-1:iec+1,jsc-1:jec+1,3)
real(kind=kind_real) :: ue(isc-1:iec+1,jsc  :jec+1,3)    ! 3D winds at edges
real(kind=kind_real) :: ve(isc  :iec+1,jsc-1:jec+1,3)    ! 3D winds at edges
real(kind=kind_real), dimension(isc:iec):: ut1, ut2, ut3
real(kind=kind_real), dimension(jsc:jec):: vt1, vt2, vt3

! Convenience to access the FV_Atm structure
Atm => self%FV_Atm

! Halfway points of the tiles
im2 = (Atm(1)%npx-1)/2
jm2 = (Atm(1)%npy-1)/2

! Allocate uatemp and vatemp
allocate(uatemp(Atm(1)%bd%isd:Atm(1)%bd%ied,Atm(1)%bd%jsd:Atm(1)%bd%jed,npz))
allocate(vatemp(Atm(1)%bd%isd:Atm(1)%bd%ied,Atm(1)%bd%jsd:Atm(1)%bd%jed,npz))

uatemp(:,:,:) = 0.0
vatemp(:,:,:) = 0.0

uatemp(isc:iec,jsc:jec,:) = ua(isc:iec,jsc:jec,:)
vatemp(isc:iec,jsc:jec,:) = va(isc:iec,jsc:jec,:)

call mpp_update_domains(uatemp, Atm(1)%domain, complete=.true.)
call mpp_update_domains(vatemp, Atm(1)%domain, complete=.true.)

do k=1, npz

  do j=jsc-1,jec+1
    do i=isc-1,iec+1
      v3(i,j,1) = uatemp(i,j,k)*Atm(1)%gridstruct%vlon(i,j,1) + vatemp(i,j,k)*Atm(1)%gridstruct%vlat(i,j,1)
      v3(i,j,2) = uatemp(i,j,k)*Atm(1)%gridstruct%vlon(i,j,2) + vatemp(i,j,k)*Atm(1)%gridstruct%vlat(i,j,2)
      v3(i,j,3) = uatemp(i,j,k)*Atm(1)%gridstruct%vlon(i,j,3) + vatemp(i,j,k)*Atm(1)%gridstruct%vlat(i,j,3)
    enddo
  enddo

  do j=jsc,jec+1
    do i=isc-1,iec+1
      ue(i,j,1) = v3(i,j-1,1) + v3(i,j,1)
      ue(i,j,2) = v3(i,j-1,2) + v3(i,j,2)
      ue(i,j,3) = v3(i,j-1,3) + v3(i,j,3)
    enddo
  enddo

  do j=jsc-1,jec+1
    do i=isc,iec+1
      ve(i,j,1) = v3(i-1,j,1) + v3(i,j,1)
      ve(i,j,2) = v3(i-1,j,2) + v3(i,j,2)
      ve(i,j,3) = v3(i-1,j,3) + v3(i,j,3)
    enddo
  enddo


  if ( isc==1 ) then
    i = 1
    do j=jsc,jec
      if ( j>jm2 ) then
        vt1(j) = Atm(1)%gridstruct%edge_vect_w(j)*ve(i,j-1,1)+(1.-Atm(1)%gridstruct%edge_vect_w(j))*ve(i,j,1)
        vt2(j) = Atm(1)%gridstruct%edge_vect_w(j)*ve(i,j-1,2)+(1.-Atm(1)%gridstruct%edge_vect_w(j))*ve(i,j,2)
        vt3(j) = Atm(1)%gridstruct%edge_vect_w(j)*ve(i,j-1,3)+(1.-Atm(1)%gridstruct%edge_vect_w(j))*ve(i,j,3)
      else
        vt1(j) = Atm(1)%gridstruct%edge_vect_w(j)*ve(i,j+1,1)+(1.-Atm(1)%gridstruct%edge_vect_w(j))*ve(i,j,1)
        vt2(j) = Atm(1)%gridstruct%edge_vect_w(j)*ve(i,j+1,2)+(1.-Atm(1)%gridstruct%edge_vect_w(j))*ve(i,j,2)
        vt3(j) = Atm(1)%gridstruct%edge_vect_w(j)*ve(i,j+1,3)+(1.-Atm(1)%gridstruct%edge_vect_w(j))*ve(i,j,3)
      endif
    enddo
    do j=jsc,jec
      ve(i,j,1) = vt1(j)
      ve(i,j,2) = vt2(j)
      ve(i,j,3) = vt3(j)
    enddo
  endif

  if ( (iec+1)==Atm(1)%npx ) then
    i = Atm(1)%npx
    do j=jsc,jec
      if ( j>jm2 ) then
        vt1(j) = Atm(1)%gridstruct%edge_vect_e(j)*ve(i,j-1,1)+(1.-Atm(1)%gridstruct%edge_vect_e(j))*ve(i,j,1)
        vt2(j) = Atm(1)%gridstruct%edge_vect_e(j)*ve(i,j-1,2)+(1.-Atm(1)%gridstruct%edge_vect_e(j))*ve(i,j,2)
        vt3(j) = Atm(1)%gridstruct%edge_vect_e(j)*ve(i,j-1,3)+(1.-Atm(1)%gridstruct%edge_vect_e(j))*ve(i,j,3)
      else
        vt1(j) = Atm(1)%gridstruct%edge_vect_e(j)*ve(i,j+1,1)+(1.-Atm(1)%gridstruct%edge_vect_e(j))*ve(i,j,1)
        vt2(j) = Atm(1)%gridstruct%edge_vect_e(j)*ve(i,j+1,2)+(1.-Atm(1)%gridstruct%edge_vect_e(j))*ve(i,j,2)
        vt3(j) = Atm(1)%gridstruct%edge_vect_e(j)*ve(i,j+1,3)+(1.-Atm(1)%gridstruct%edge_vect_e(j))*ve(i,j,3)
      endif
    enddo
    do j=jsc,jec
      ve(i,j,1) = vt1(j)
      ve(i,j,2) = vt2(j)
      ve(i,j,3) = vt3(j)
    enddo
  endif

  if ( jsc==1 ) then
    j = 1
    do i=isc,iec
      if ( i>im2 ) then
        ut1(i) = Atm(1)%gridstruct%edge_vect_s(i)*ue(i-1,j,1)+(1.-Atm(1)%gridstruct%edge_vect_s(i))*ue(i,j,1)
        ut2(i) = Atm(1)%gridstruct%edge_vect_s(i)*ue(i-1,j,2)+(1.-Atm(1)%gridstruct%edge_vect_s(i))*ue(i,j,2)
        ut3(i) = Atm(1)%gridstruct%edge_vect_s(i)*ue(i-1,j,3)+(1.-Atm(1)%gridstruct%edge_vect_s(i))*ue(i,j,3)
      else
        ut1(i) = Atm(1)%gridstruct%edge_vect_s(i)*ue(i+1,j,1)+(1.-Atm(1)%gridstruct%edge_vect_s(i))*ue(i,j,1)
        ut2(i) = Atm(1)%gridstruct%edge_vect_s(i)*ue(i+1,j,2)+(1.-Atm(1)%gridstruct%edge_vect_s(i))*ue(i,j,2)
        ut3(i) = Atm(1)%gridstruct%edge_vect_s(i)*ue(i+1,j,3)+(1.-Atm(1)%gridstruct%edge_vect_s(i))*ue(i,j,3)
      endif
    enddo
    do i=isc,iec
      ue(i,j,1) = ut1(i)
      ue(i,j,2) = ut2(i)
      ue(i,j,3) = ut3(i)
    enddo
  endif

  if ( (jec+1)==Atm(1)%npy ) then
    j = Atm(1)%npy
    do i=isc,iec
      if ( i>im2 ) then
        ut1(i) = Atm(1)%gridstruct%edge_vect_n(i)*ue(i-1,j,1)+(1.-Atm(1)%gridstruct%edge_vect_n(i))*ue(i,j,1)
        ut2(i) = Atm(1)%gridstruct%edge_vect_n(i)*ue(i-1,j,2)+(1.-Atm(1)%gridstruct%edge_vect_n(i))*ue(i,j,2)
        ut3(i) = Atm(1)%gridstruct%edge_vect_n(i)*ue(i-1,j,3)+(1.-Atm(1)%gridstruct%edge_vect_n(i))*ue(i,j,3)
      else
        ut1(i) = Atm(1)%gridstruct%edge_vect_n(i)*ue(i+1,j,1)+(1.-Atm(1)%gridstruct%edge_vect_n(i))*ue(i,j,1)
        ut2(i) = Atm(1)%gridstruct%edge_vect_n(i)*ue(i+1,j,2)+(1.-Atm(1)%gridstruct%edge_vect_n(i))*ue(i,j,2)
        ut3(i) = Atm(1)%gridstruct%edge_vect_n(i)*ue(i+1,j,3)+(1.-Atm(1)%gridstruct%edge_vect_n(i))*ue(i,j,3)
      endif
    enddo
    do i=isc,iec
      ue(i,j,1) = ut1(i)
      ue(i,j,2) = ut2(i)
      ue(i,j,3) = ut3(i)
    enddo
  endif

  do j=jsc,jec
    do i=isc,iec
      ud(i,j,k) = 0.5*( ue(i,j,1)*Atm(1)%gridstruct%es(1,i,j,1) +  &
                        ue(i,j,2)*Atm(1)%gridstruct%es(2,i,j,1) +  &
                        ue(i,j,3)*Atm(1)%gridstruct%es(3,i,j,1) )
    enddo
  enddo

  do j=jsc,jec
    do i=isc,iec
      vd(i,j,k) = 0.5*( ve(i,j,1)*Atm(1)%gridstruct%ew(1,i,j,2) +  &
                        ve(i,j,2)*Atm(1)%gridstruct%ew(2,i,j,2) +  &
                        ve(i,j,3)*Atm(1)%gridstruct%ew(3,i,j,2) )
    enddo
  enddo

enddo

end subroutine a_to_d_nl

! ------------------------------------------------------------------------------

end module fv3jedi_lm_dynamics_mod
