module fv3jedi_lm_mod

use fv3jedi_lm_kinds_mod
use fv3jedi_lm_utils_mod

use fv3jedi_lm_dynamics_mod, only: fv3jedi_lm_dynamics_type
use fv3jedi_lm_physics_mod,  only: fv3jedi_lm_physics_type

!> Top level for fv3jedi linearized model

!> Combines fv3 tlm/adm with the GEOS tlm/adm physics
!> All developed by NASA's Global Modeling and Assimilation Office
!> daniel.holdaway@nasa.gov, Code 610.1 Goddard Space Flight Center,
!> Greenbelt, MD 20771 USA

implicit none
private
public :: fv3jedi_lm_type

type fv3jedi_lm_type
 type(fv3jedi_lm_conf) :: conf
 type(fv3jedi_lm_pert) :: pert
 type(fv3jedi_lm_traj) :: traj
 type(fv3jedi_lm_dynamics_type) :: fv3jedi_lm_dynamics
 type(fv3jedi_lm_physics_type)  :: fv3jedi_lm_physics
 contains
  procedure :: create
  procedure :: allocate_tracers
  procedure :: read_d_grid_winds
  procedure :: write_d_grid_winds
  procedure :: store_winds
  procedure :: reinitialize_winds
  procedure :: initialize_dwinds_from_awinds
  procedure :: init_nl
  procedure :: init_tl
  procedure :: init_ad
  procedure :: step_nl
  procedure :: step_tl
  procedure :: step_ad
  procedure :: final_nl
  procedure :: final_tl
  procedure :: final_ad
  procedure :: delete
end type

contains

! ------------------------------------------------------------------------------

subroutine create(self,dt,npx,npy,npz,ptop,ak,bk)

 class(fv3jedi_lm_type), intent(inout) :: self

 real(kind=kind_real), intent(in) :: dt, ptop
 integer, intent(in) :: npx,npy,npz
 real(kind=kind_real), intent(in) :: ak(npz+1), bk(npz+1)

 self%conf%dt = dt
 self%conf%ptop = ptop

 allocate(self%conf%ak(npz+1))
 allocate(self%conf%bk(npz+1))
 self%conf%ak = ak
 self%conf%bk = bk

 !Call dynamics create, provides the model grid
 call self%fv3jedi_lm_dynamics%create(self%conf)

 !Make grid available to all components
 self%conf%isc = self%fv3jedi_lm_dynamics%FV_Atm(1)%bd%isc
 self%conf%iec = self%fv3jedi_lm_dynamics%FV_Atm(1)%bd%iec
 self%conf%jsc = self%fv3jedi_lm_dynamics%FV_Atm(1)%bd%jsc
 self%conf%jec = self%fv3jedi_lm_dynamics%FV_Atm(1)%bd%jec
 self%conf%isd = self%fv3jedi_lm_dynamics%FV_Atm(1)%bd%isd
 self%conf%ied = self%fv3jedi_lm_dynamics%FV_Atm(1)%bd%ied
 self%conf%jsd = self%fv3jedi_lm_dynamics%FV_Atm(1)%bd%jsd
 self%conf%jed = self%fv3jedi_lm_dynamics%FV_Atm(1)%bd%jed
 self%conf%npx = self%fv3jedi_lm_dynamics%FV_Atm(1)%npx
 self%conf%npy = self%fv3jedi_lm_dynamics%FV_Atm(1)%npy
 self%conf%npz = self%fv3jedi_lm_dynamics%FV_Atm(1)%npz
 self%conf%hydrostatic = self%fv3jedi_lm_dynamics%FV_Atm(1)%flagstruct%hydrostatic

 !Physics grid
 self%conf%im = (self%conf%iec-self%conf%isc+1)
 self%conf%jm = (self%conf%jec-self%conf%jsc+1)
 self%conf%lm = self%conf%npz

 !All physics switch
 if (self%conf%do_phy_trb == 0 .and. self%conf%do_phy_mst == 0) self%conf%do_phy = 0

 !Call the physics create
 if (self%conf%do_phy == 1) call self%fv3jedi_lm_physics%create(self%conf)

 !Sanity checks
 if ((self%conf%npx .ne. npx) .or. (self%conf%npy .ne. npy) .or. (self%conf%npz .ne. npz) ) then
  if (self%conf%rpe) print*, 'fv3jedi tlm/adm problem: dynamics creating different grid than inputs created for'
  call exit(1)
 endif

 !Allocate main traj and pert structures
 call allocate_traj(self%traj,self%conf%isc,self%conf%iec,self%conf%jsc,self%conf%jec,&
                    self%conf%npz,self%conf%hydrostatic,self%conf%do_phy_mst)
 call allocate_pert(self%pert,self%conf%isc,self%conf%iec,self%conf%jsc,self%conf%jec,self%conf%npz,self%conf%hydrostatic)

endsubroutine create

! ------------------------------------------------------------------------------

subroutine allocate_tracers(self, ntracers, include_pert_tracers)

 class(fv3jedi_lm_type), intent(inout) :: self
 integer, intent(in) :: ntracers
 logical, intent(in) :: include_pert_tracers

 call allocate_traj_tracers(self%traj, self%conf%isc, self%conf%iec, self%conf%jsc, self%conf%jec, &
                            self%conf%npz, ntracers)
 if (include_pert_tracers) then
   call allocate_pert_tracers(self%pert, self%conf%isc, self%conf%iec, self%conf%jsc, &
                              self%conf%jec, self%conf%npz, ntracers)
 end if

endsubroutine allocate_tracers

! ------------------------------------------------------------------------------

subroutine read_d_grid_winds(self, datapath, filename)

  ! Args
  class(fv3jedi_lm_type), intent(inout) :: self
  character(len=*), intent(in) :: datapath, filename

  ! Call dynamics restart init
  call self%fv3jedi_lm_dynamics%read_d_grid_winds(datapath, filename, self%traj)

endsubroutine read_d_grid_winds

! ------------------------------------------------------------------------------

subroutine write_d_grid_winds(self, datapath, filename)

  ! Args
  class(fv3jedi_lm_type), intent(inout) :: self
  character(len=*), intent(in) :: datapath, filename

  ! Call dynamics restart init
  call self%fv3jedi_lm_dynamics%write_d_grid_winds(datapath, filename, self%traj)

endsubroutine write_d_grid_winds

! ------------------------------------------------------------------------------

subroutine store_winds(self)

  ! Args
  class(fv3jedi_lm_type), intent(inout) :: self

  ! Winds that are stored for use in reinitialize_winds
  if (.not.allocated(self%traj%u_stored )) &
    allocate(self%traj%u_stored (self%conf%isc:self%conf%iec, self%conf%jsc:self%conf%jec, self%conf%npz))
  if (.not.allocated(self%traj%v_stored )) &
    allocate(self%traj%v_stored (self%conf%isc:self%conf%iec, self%conf%jsc:self%conf%jec, self%conf%npz))
  if (.not.allocated(self%traj%ua_stored)) &
    allocate(self%traj%ua_stored(self%conf%isc:self%conf%iec, self%conf%jsc:self%conf%jec, self%conf%npz))
  if (.not.allocated(self%traj%va_stored)) &
    allocate(self%traj%va_stored(self%conf%isc:self%conf%iec, self%conf%jsc:self%conf%jec, self%conf%npz))

  self%traj%u_stored = self%traj%u
  self%traj%v_stored = self%traj%v
  self%traj%ua_stored = self%traj%ua
  self%traj%va_stored = self%traj%va

endsubroutine store_winds

! ------------------------------------------------------------------------------

subroutine reinitialize_winds(self)

  ! Args
  class(fv3jedi_lm_type), intent(inout) :: self

  ! Locals
  real(kind=kind_real), allocatable, dimension(:,:,:) :: ua_pert, va_pert
  real(kind=kind_real), allocatable, dimension(:,:,:) :: ud_pert, vd_pert

  ! Compute the difference between the current and stored A-Grid winds
  allocate(ua_pert(self%conf%isc:self%conf%iec, self%conf%jsc:self%conf%jec, self%conf%npz))
  allocate(va_pert(self%conf%isc:self%conf%iec, self%conf%jsc:self%conf%jec, self%conf%npz))
  ua_pert = self%traj%ua - self%traj%ua_stored
  va_pert = self%traj%va - self%traj%va_stored

  ! Convert the perturbations to D-Grid
  allocate(ud_pert(self%conf%isc:self%conf%iec, self%conf%jsc:self%conf%jec, self%conf%npz))
  allocate(vd_pert(self%conf%isc:self%conf%iec, self%conf%jsc:self%conf%jec, self%conf%npz))
  call self%fv3jedi_lm_dynamics%a_to_d_nl(self%conf%isc, self%conf%iec, self%conf%jsc, &
                                       self%conf%jec, self%conf%npz, ua_pert, va_pert, ud_pert, &
                                       vd_pert)

  ! Add the D-Grid wind perturbation to the D-Grid winds in the trajectory
  self%traj%u = self%traj%u_stored + ud_pert
  self%traj%v = self%traj%v_stored + vd_pert

endsubroutine reinitialize_winds

! ------------------------------------------------------------------------------

subroutine initialize_dwinds_from_awinds(self, ua, va)

  ! Args
  class(fv3jedi_lm_type), intent(inout) :: self
  real(kind=kind_real),   intent(in)    :: ua(:,:,:)
  real(kind=kind_real),   intent(in)    :: va(:,:,:)

  ! Locals
  real(kind=kind_real), allocatable, dimension(:,:,:) :: ud, vd

  ! Convert the A-Grid winds to D-Grid
  allocate(ud(self%conf%isc:self%conf%iec, self%conf%jsc:self%conf%jec, self%conf%npz))
  allocate(vd(self%conf%isc:self%conf%iec, self%conf%jsc:self%conf%jec, self%conf%npz))
  call self%fv3jedi_lm_dynamics%a_to_d_nl(self%conf%isc, self%conf%iec, self%conf%jsc, &
                                          self%conf%jec, self%conf%npz, ua, va, ud, vd)

  ! Initialize the D-Grid winds in the trajectory
  self%traj%u = ud
  self%traj%v = vd

endsubroutine initialize_dwinds_from_awinds

! ------------------------------------------------------------------------------

subroutine init_nl(self)

 class(fv3jedi_lm_type), intent(inout) :: self

 if (self%conf%do_dyn == 1) call self%fv3jedi_lm_dynamics%init_nl(self%conf,self%pert,self%traj)
 if (self%conf%do_phy == 1) call self%fv3jedi_lm_physics%init_nl(self%conf,self%pert,self%traj)

endsubroutine init_nl

! ------------------------------------------------------------------------------

subroutine init_tl(self)

 class(fv3jedi_lm_type), intent(inout) :: self

 call ipert_to_zero(self%pert)

 if (self%conf%do_dyn == 1) call self%fv3jedi_lm_dynamics%init_tl(self%conf,self%pert,self%traj)
 if (self%conf%do_phy == 1) call self%fv3jedi_lm_physics%init_tl(self%conf,self%pert,self%traj)

endsubroutine init_tl

! ------------------------------------------------------------------------------

subroutine init_ad(self)

 class(fv3jedi_lm_type), intent(inout) :: self

 call ipert_to_zero(self%pert)

 if (self%conf%do_dyn == 1) call self%fv3jedi_lm_dynamics%init_ad(self%conf,self%pert,self%traj)
 if (self%conf%do_phy == 1) call self%fv3jedi_lm_physics%init_ad(self%conf,self%pert,self%traj)

endsubroutine init_ad

! ------------------------------------------------------------------------------

subroutine step_nl(self)

 class(fv3jedi_lm_type), intent(inout) :: self

 if (self%conf%do_dyn == 1) call self%fv3jedi_lm_dynamics%step_nl(self%conf,self%traj)
 if (self%conf%do_phy == 1) call self%fv3jedi_lm_physics%step_nl(self%conf,self%traj)

endsubroutine step_nl

! ------------------------------------------------------------------------------

subroutine step_tl(self)

 class(fv3jedi_lm_type), intent(inout) :: self

 call ipert_to_zero(self%pert)
 if (self%conf%do_dyn == 1) call self%fv3jedi_lm_dynamics%step_tl(self%conf,self%traj,self%pert)
 if (self%conf%do_phy == 1) call self%fv3jedi_lm_physics%step_tl(self%conf,self%traj,self%pert)
 call ipert_to_zero(self%pert)

endsubroutine step_tl

! ------------------------------------------------------------------------------

subroutine step_ad(self)

 class(fv3jedi_lm_type), intent(inout) :: self

 call ipert_to_zero(self%pert)
 if (self%conf%do_phy == 1) call self%fv3jedi_lm_physics%step_ad(self%conf,self%traj,self%pert)
 if (self%conf%do_dyn == 1) call self%fv3jedi_lm_dynamics%step_ad(self%conf,self%traj,self%pert)
 call ipert_to_zero(self%pert)

endsubroutine step_ad

! ------------------------------------------------------------------------------

subroutine final_nl(self)

 class(fv3jedi_lm_type), intent(inout) :: self

endsubroutine final_nl

! ------------------------------------------------------------------------------

subroutine final_tl(self)

 class(fv3jedi_lm_type), intent(inout) :: self

 call ipert_to_zero(self%pert)

endsubroutine final_tl

! ------------------------------------------------------------------------------

subroutine final_ad(self)

 class(fv3jedi_lm_type), intent(inout) :: self

 call ipert_to_zero(self%pert)

endsubroutine final_ad

! ------------------------------------------------------------------------------

subroutine delete(self)

 class(fv3jedi_lm_type), intent(inout) :: self

 if (self%conf%do_dyn == 1) call self%fv3jedi_lm_dynamics%delete(self%conf)
 if (self%conf%do_phy == 1) call self%fv3jedi_lm_physics%delete(self%conf)

 deallocate(self%conf%ak)
 deallocate(self%conf%bk)

 call deallocate_traj(self%traj)
 call deallocate_pert(self%pert)

endsubroutine delete

! ------------------------------------------------------------------------------

subroutine ipert_to_zero(pert)

 !> Intenal part of pert to zero

 type(fv3jedi_lm_pert), intent(inout) :: pert

 pert%ua = 0.0_kind_real
 pert%va = 0.0_kind_real
 pert%cfcn = 0.0_kind_real

endsubroutine ipert_to_zero

! ------------------------------------------------------------------------------

end module fv3jedi_lm_mod
