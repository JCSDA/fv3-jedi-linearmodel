module fv3jedi_lm_moist_mod

use fv3jedi_lm_utils_mod
use fv3jedi_lm_kinds_mod
use fv3jedi_lm_const_mod

use MAPL_ConstantsMod

! Moist Schemes
use convection
use convection_ad
use convection_tl
use cloud
use cloud_ad
use cloud_tl

! Saturation table
use qsat_util

use tapenade_iter, only: cp_iter, cp_iter_controls, initialize_cp_iter, finalize_cp_iter
use tapenade_iter, only: cp_mod_ini, cp_mod_mid, cp_mod_end, pushrealarray, poprealarray

!> Top level for fv3jedi linearized model

implicit none
private
public :: fv3jedi_lm_moist_type

!> Local trajectory objects
type local_traj_moist
  real(8), allocatable, dimension(:,:,:) :: UT, VT, PTT, QVT
  real(8), allocatable, dimension(:,:)   :: TS, FRLAND
  integer, allocatable, dimension(:,:)   :: KCBL, KHu, KHl
  real(8), allocatable, dimension(:,:,:) :: PLE, CNV_PLE, PK
  real(8), allocatable, dimension(:,:,:) :: CNV_DQLDTT, CNV_MFDT, CNV_PRC3T, CNV_UPDFT
  real(8), allocatable, dimension(:,:,:) :: PTT_C, QVT_C
  real(8), allocatable, dimension(:,:,:) :: CNV_DQLDTT_C, CNV_MFDT_C, CNV_PRC3T_C, CNV_UPDFT_C
  integer, allocatable, dimension(:,:)   :: SEEDRAS
  real(8), allocatable, dimension(:,:)   :: CO_AUTO
  integer, allocatable, dimension(:,:)   :: DOCONVEC
  real(8), allocatable, dimension(:,:,:) :: WGT0, WGT1
  real(8), allocatable, dimension(:,:,:) :: QILST, QLLST, QICNT, QLCNT
  real(8), allocatable, dimension(:,:,:) :: CFLST, CFCNT
  real(8), allocatable, dimension(:,:,:) :: ILSF, ICNF, LLSF, LCNF
  logical :: set = .false.
endtype local_traj_moist

!> Local perturbation objects
type local_pert_moist
  real(8), allocatable, dimension(:,:,:) :: UP, VP, PTP, QVP
  real(8), allocatable, dimension(:,:,:) :: CNV_DQLDTP, CNV_MFDP, CNV_PRC3P, CNV_UPDFP
  real(8), allocatable, dimension(:,:,:) :: QILSP, QLLSP, QICNP, QLCNP
  real(8), allocatable, dimension(:,:,:) :: CFLSP, CFCNP
endtype local_pert_moist

!> Local constants objects
type local_cnst_moist
  integer :: ICMIN
  real(8), allocatable, dimension(:) :: ESTBLX, SIGE
  real(8) :: RASPARAMS(25),  CLOUDPARAMS (57)
  real(8) :: MAPL8_CP, MAPL8_ALHL, MAPL8_GRAV, MAPL8_P00, MAPL8_KAPPA
  real(8) :: MAPL8_RGAS, MAPL8_H2OMW, MAPL8_AIRMW, MAPL8_VIREPS
  real(8) :: MAPL8_RUNIV, MAPL8_ALHF, MAPL8_PI, MAPL8_ALHS
  real(8) :: MAPL8_TICE, MAPL8_RVAP
endtype local_cnst_moist

!> Moist class (self)
type fv3jedi_lm_moist_type
 type(local_traj_moist), allocatable :: ltraj(:)
 type(local_pert_moist) :: lpert
 type(local_cnst_moist) :: lcnst
 contains
  procedure :: create
  procedure :: init_nl
  procedure :: init_tl
  procedure :: init_ad
  procedure :: step_nl
  procedure :: step_tl
  procedure :: step_ad
  procedure :: delete
end type fv3jedi_lm_moist_type

contains

! ------------------------------------------------------------------------------

subroutine create(self,conf)

 implicit none

 class(fv3jedi_lm_moist_type), target, intent(inout) :: self
 type(fv3jedi_lm_conf), intent(in)    :: conf

 integer :: imsize, L, n
 real(8), allocatable, dimension(:) :: pref
 real(8), parameter :: PMIN_DET = 3000.0
 integer, parameter :: DEGSUBS   = 100
 real(8), parameter :: TMINTBL   = 150.0, TMAXTBL = 333.0
 integer, parameter :: TABLESIZE = nint(TMAXTBL-TMINTBL)*DEGSUBS + 1

 self%lcnst%MAPL8_CP     = dble(MAPL_CP)
 self%lcnst%MAPL8_ALHL   = dble(MAPL_ALHL)
 self%lcnst%MAPL8_GRAV   = dble(MAPL_GRAV)
 self%lcnst%MAPL8_P00    = dble(MAPL_P00)
 self%lcnst%MAPL8_KAPPA  = dble(MAPL_KAPPA)
 self%lcnst%MAPL8_RGAS   = dble(MAPL_RGAS)
 self%lcnst%MAPL8_H2OMW  = dble(MAPL_H2OMW)
 self%lcnst%MAPL8_AIRMW  = dble(MAPL_AIRMW)
 self%lcnst%MAPL8_VIREPS = dble(MAPL_VIREPS)
 self%lcnst%MAPL8_RUNIV  = dble(MAPL_RUNIV)
 self%lcnst%MAPL8_ALHF   = dble(MAPL_ALHF)
 self%lcnst%MAPL8_PI     = dble(MAPL_PI)
 self%lcnst%MAPL8_ALHS   = dble(MAPL_ALHS)
 self%lcnst%MAPL8_TICE   = dble(MAPL_TICE)
 self%lcnst%MAPL8_RVAP   = dble(MAPL_RVAP)

 imsize = conf%im*4

 !RAS Parameters
 self%lcnst%RASPARAMS( 1) = 1.000
 self%lcnst%RASPARAMS( 2) = 0.05
 self%lcnst%RASPARAMS( 3) = 0.0   ! NOW IN CO_AUTO (CONTAINED WITHIN SUBROUTINE)
 self%lcnst%RASPARAMS( 4) = 8.0e-4
 self%lcnst%RASPARAMS( 5) = 1800.
 self%lcnst%RASPARAMS( 6) = 43200.0
 self%lcnst%RASPARAMS( 7) = -300. !RASNCL, CONTROLS FINDDTLS, USE OF RANDOM NUMBER
 self%lcnst%RASPARAMS( 8) = 4.0
 self%lcnst%RASPARAMS( 9) = 0.0
 self%lcnst%RASPARAMS(10) = 200.
 self%lcnst%RASPARAMS(11) = 7.5e-4
 self%lcnst%RASPARAMS(12) = 1.0
 self%lcnst%RASPARAMS(13) =-1.0
 self%lcnst%RASPARAMS(14) = 1.3
 self%lcnst%RASPARAMS(15) = 1.3
 self%lcnst%RASPARAMS(16) = 263.
 self%lcnst%RASPARAMS(17) = 0.5
 self%lcnst%RASPARAMS(18) = 1.0
 self%lcnst%RASPARAMS(19) = 0.0
 self%lcnst%RASPARAMS(20) = 0.1
 self%lcnst%RASPARAMS(21) = 0.8
 self%lcnst%RASPARAMS(22) = 1.0
 if( imsize .le. 200                      ) self%lcnst%RASPARAMS(23) = 4000.0
 if( imsize .gt. 200 .and. imsize.le.400  ) self%lcnst%RASPARAMS(23) = 2000.0
 if( imsize .gt. 400 .and. imsize.le.800  ) self%lcnst%RASPARAMS(23) = 700.0
 if( imsize .gt. 800 .and. imsize.le.1600 ) self%lcnst%RASPARAMS(23) = 450.0
 if( imsize .gt. 1600                     ) self%lcnst%RASPARAMS(23) = 450.0
 self%lcnst%RASPARAMS(24) = 0.5
 self%lcnst%RASPARAMS(25) = 0.65

 !SET SBAC PARAMETERS
 self%lcnst%CLOUDPARAMS( 1) = 10.0
 self%lcnst%CLOUDPARAMS( 2) = 4.0
 self%lcnst%CLOUDPARAMS( 3) = 4.0
 self%lcnst%CLOUDPARAMS( 4) = 1.0
 self%lcnst%CLOUDPARAMS( 5) = 2.0e-3
 self%lcnst%CLOUDPARAMS( 6) = 8.0e-4
 self%lcnst%CLOUDPARAMS( 7) = 2.0
 self%lcnst%CLOUDPARAMS( 8) = 1.0
 self%lcnst%CLOUDPARAMS( 9) = -1.0
 self%lcnst%CLOUDPARAMS(10) = 0.0
 self%lcnst%CLOUDPARAMS(11) = 1.3
 self%lcnst%CLOUDPARAMS(12) = 1.0e-9
 self%lcnst%CLOUDPARAMS(13) = 3.3e-4
 self%lcnst%CLOUDPARAMS(14) = 20.
 self%lcnst%CLOUDPARAMS(15) = 4.8
 self%lcnst%CLOUDPARAMS(16) = 4.8
 self%lcnst%CLOUDPARAMS(17) = 230.
 self%lcnst%CLOUDPARAMS(18) = 1.0
 self%lcnst%CLOUDPARAMS(19) = 1.0
 self%lcnst%CLOUDPARAMS(20) = 230.
 self%lcnst%CLOUDPARAMS(21) = 14400.
 self%lcnst%CLOUDPARAMS(22) = 50.
 self%lcnst%CLOUDPARAMS(23) = 0.01
 self%lcnst%CLOUDPARAMS(24) = 0.1
 self%lcnst%CLOUDPARAMS(25) = 200.
 self%lcnst%CLOUDPARAMS(26) = 0.
 self%lcnst%CLOUDPARAMS(27) = 0.
 self%lcnst%CLOUDPARAMS(28) = 0.5
 self%lcnst%CLOUDPARAMS(29) = 0.5
 self%lcnst%CLOUDPARAMS(30) = 2000.
 self%lcnst%CLOUDPARAMS(31) = 0.8
 self%lcnst%CLOUDPARAMS(32) = 0.5
 self%lcnst%CLOUDPARAMS(33) = -40.0
 self%lcnst%CLOUDPARAMS(34) = 1.0
 self%lcnst%CLOUDPARAMS(35) = 4.0
 self%lcnst%CLOUDPARAMS(36) = 0.0
 self%lcnst%CLOUDPARAMS(37) = 0.0
 self%lcnst%CLOUDPARAMS(38) = 0.0
 self%lcnst%CLOUDPARAMS(39) = 1.0e-3
 self%lcnst%CLOUDPARAMS(40) = 8.0e-4
 self%lcnst%CLOUDPARAMS(41) = 1.0
 if( imsize .le. 200                      ) self%lcnst%CLOUDPARAMS(42) = 0.80
 if( imsize .gt. 200 .and. imsize.le.400  ) self%lcnst%CLOUDPARAMS(42) = 0.90
 if( imsize .gt. 400 .and. imsize.le.800  ) self%lcnst%CLOUDPARAMS(42) = 0.93
 if( imsize .gt. 800 .and. imsize.le.1600 ) self%lcnst%CLOUDPARAMS(42) = 0.95
 if( imsize .gt. 1600                     ) self%lcnst%CLOUDPARAMS(42) = 0.97
 self%lcnst%CLOUDPARAMS(43) = 1.0
 self%lcnst%CLOUDPARAMS(44) = 0.0
 self%lcnst%CLOUDPARAMS(45) = 750.0
 self%lcnst%CLOUDPARAMS(46) = self%lcnst%CLOUDPARAMS(42)+0.01
 self%lcnst%CLOUDPARAMS(47) = 1.0
 self%lcnst%CLOUDPARAMS(48) = 1.0
 self%lcnst%CLOUDPARAMS(49) = 0.0
 self%lcnst%CLOUDPARAMS(50) = 0.0
 self%lcnst%CLOUDPARAMS(51) = 10.e-6
 self%lcnst%CLOUDPARAMS(52) = 20.e-6
 self%lcnst%CLOUDPARAMS(53) = 21.e-6
 self%lcnst%CLOUDPARAMS(54) = 40.e-6
 self%lcnst%CLOUDPARAMS(55) = 30.e-6
 self%lcnst%CLOUDPARAMS(56) = 1.0
 self%lcnst%CLOUDPARAMS(57) = 1.0

 if (conf%saveltraj) then
   allocate(self%ltraj(conf%nt))
   do n = 1,conf%nt
     call allocate_ltraj(conf%im,conf%jm,conf%lm,self%ltraj(n))
   enddo
 else
   allocate(self%ltraj(1))
   call allocate_ltraj(conf%im,conf%jm,conf%lm,self%ltraj(1))
 endif

 call allocate_lpert(conf%im,conf%jm,conf%lm,self%lpert)

 !Allocate self
 allocate(self%lcnst%ESTBLX(TABLESIZE))
 allocate(self%lcnst%SIGE(0:conf%LM))

 !ESTBLX
 call ESINIT(self%lcnst%ESTBLX)

 allocate(PREF(0:conf%lm))
 do L = 0,conf%lm
   PREF(L) = conf%ak(L+1) + conf%bk(L+1)*self%lcnst%MAPL8_P00
 enddo

 self%lcnst%ICMIN = max(1,count(PREF < PMIN_DET))
 self%lcnst%SIGE = PREF/PREF(conf%LM)

 deallocate(pref)

endsubroutine create

! ------------------------------------------------------------------------------

subroutine init_nl(self,pert,traj)

 implicit none

 class(fv3jedi_lm_moist_type), intent(inout) :: self
 type(fv3jedi_lm_pert), intent(inout) :: pert
 type(fv3jedi_lm_traj), intent(in) :: traj

endsubroutine init_nl

! ------------------------------------------------------------------------------

subroutine init_tl(self,pert,traj)

 implicit none

 class(fv3jedi_lm_moist_type), intent(inout) :: self
 type(fv3jedi_lm_pert), intent(inout) :: pert
 type(fv3jedi_lm_traj), intent(in) :: traj

endsubroutine init_tl

! ------------------------------------------------------------------------------

subroutine init_ad(self,pert,traj)

 implicit none

 class(fv3jedi_lm_moist_type), intent(inout) :: self
 type(fv3jedi_lm_pert), intent(inout) :: pert
 type(fv3jedi_lm_traj), intent(in) :: traj

endsubroutine init_ad

! ------------------------------------------------------------------------------

subroutine step_nl(self,conf,traj)

 implicit none

 class(fv3jedi_lm_moist_type), target, intent(inout) :: self
 type(fv3jedi_lm_traj), target,        intent(inout) :: traj
 type(fv3jedi_lm_conf),                 intent(in)   :: conf

 integer :: i,j,im,jm,lm,isc,iec,jsc,jec
 integer :: it_qv,it_qi,it_ql
 type(local_traj_moist), pointer :: ltraj
 type(local_pert_moist), pointer :: lpert
 type(local_cnst_moist), pointer :: lcnst

 im = conf%im
 jm = conf%jm
 lm = conf%lm
 isc = conf%isc
 iec = conf%iec
 jsc = conf%jsc
 jec = conf%jec

 !Convenience pointers
 lpert => self%lpert
 lcnst => self%lcnst
 if (conf%saveltraj) then
   ltraj => self%ltraj(conf%n)
 else
   ltraj => self%ltraj(1)
 endif

 !Set up the local trajectory
 if (.not. ltraj%set) call set_ltraj(conf,lcnst,traj,ltraj)

 !Local pert (not really needed)
 lpert%up    = 0.0_8
 lpert%vp    = 0.0_8
 lpert%ptp   = 0.0_8
 lpert%qvp   = 0.0_8
 lpert%cflsp = 0.0_8
 lpert%cfcnp = 0.0_8
 lpert%qilsp = 0.0_8
 lpert%qicnp = 0.0_8
 lpert%qllsp = 0.0_8
 lpert%qlcnp = 0.0_8

 ltraj%cnv_dqldtt  = 0.0_8
 ltraj%cnv_mfdt    = 0.0_8
 ltraj%cnv_prc3t   = 0.0_8
 ltraj%cnv_updft   = 0.0_8
 lpert%cnv_dqldtp  = 0.0_8
 lpert%cnv_mfdp    = 0.0_8
 lpert%cnv_prc3p   = 0.0_8
 lpert%cnv_updfp   = 0.0_8

 !Call the tangent linear convection scheme
 do i = 1,conf%im
   do j = 1,conf%jm

     if ( ltraj%doconvec(i,j) == 1) then
       call rase_d( 1, 1, conf%lm, lcnst%icmin, conf%dt,              &
                    lcnst%mapl8_cp, lcnst%mapl8_alhl,                 &
                    lcnst%mapl8_grav, lcnst%mapl8_rgas,               &
                    lcnst%mapl8_h2omw, lcnst%mapl8_airmw,             &
                    lcnst%mapl8_vireps,                               &
                    ltraj%seedras(i,j), lcnst%sige,                   &
                    ltraj%kcbl(i,j),                                  &
                    ltraj%wgt0(i,j,:), ltraj%wgt1(i,j,:),             &
                    ltraj%frland(i,j), ltraj%ts(i,j),                 &
                    ltraj%ptt(i,j,:), lpert%ptp(i,j,:),               &
                    ltraj%qvt(i,j,:), lpert%qvp(i,j,:),               &
                    ltraj%ut(i,j,:), lpert%up(i,j,:),                 &
                    ltraj%vt(i,j,:), lpert%vp(i,j,:),                 &
                    ltraj%co_auto(i,j), ltraj%cnv_ple(i,j,:),         &
                    ltraj%cnv_dqldtt(i,j,:), lpert%cnv_dqldtp(i,j,:), &
                    ltraj%cnv_mfdt(i,j,:),   lpert%cnv_mfdp(i,j,:),   &
                    ltraj%cnv_prc3t(i,j,:),  lpert%cnv_prc3p(i,j,:),  &
                    ltraj%cnv_updft(i,j,:),  lpert%cnv_updfp(i,j,:),  &
                    lcnst%rasparams, lcnst%estblx                     )
     endif

   enddo
 enddo

 !Call the tangent linear cloud scheme.
 call cloud_driver_d ( conf%dt, conf%im, conf%jm, conf%lm,                                         &
                       ltraj%ptt_c, lpert%ptp,                                                     &
                       ltraj%qvt_c, lpert%qvp,                                                     &
                       ltraj%ple,                                                                  &
                       ltraj%cnv_dqldtt_c, lpert%cnv_dqldtp, ltraj%cnv_mfdt_c,  lpert%cnv_mfdp,    &
                       ltraj%cnv_prc3t_c,  lpert%cnv_prc3p,  ltraj%cnv_updft_c, lpert%cnv_updfp,   &
                       ltraj%qilst, lpert%qilsp, ltraj%qllst, lpert%qllsp,                         &
                       ltraj%qicnt, lpert%qicnp, ltraj%qlcnt, lpert%qlcnp,                         &
                       ltraj%cflst, lpert%cflsp, ltraj%cfcnt, lpert%cfcnp,                         &
                       ltraj%frland, lcnst%cloudparams, lcnst%estblx, ltraj%khu, ltraj%khl,        &
                       lcnst%mapl8_runiv, lcnst%mapl8_kappa, lcnst%mapl8_airmw, lcnst%mapl8_h2omw, &
                       lcnst%mapl8_grav, lcnst%mapl8_alhl, lcnst%mapl8_alhf,   lcnst%mapl8_pi,     &
                       lcnst%mapl8_rgas, lcnst%mapl8_cp,   lcnst%mapl8_vireps, lcnst%mapl8_alhs,   &
                       lcnst%mapl8_tice, lcnst%mapl8_rvap, lcnst%mapl8_p00, conf%do_phy_mst        )

 !Back to traj
 traj%u(isc:iec,jsc:jec,:) = real(ltraj%ut(1:im,1:jm,:),kind_real)
 traj%v(isc:iec,jsc:jec,:) = real(ltraj%vt(1:im,1:jm,:),kind_real)
 traj%t(isc:iec,jsc:jec,:) = real(ltraj%pk(1:im,1:jm,:) * ltraj%PTT(1:im,1:jm,:) / p00**kappa,kind_real)
 traj%cfcn(isc:iec,jsc:jec,:) = real(ltraj%cfcnt(1:im,1:jm,:),kind_real)

 call get_tracer_index(traj, 'water_vapor_mixing_ratio_wrt_moist_air', it_qv)
 call get_tracer_index(traj, 'cloud_liquid_ice', it_qi)
 call get_tracer_index(traj, 'cloud_liquid_water', it_ql)
 traj%tracers(isc:iec,jsc:jec,:,it_qv) = real(ltraj%qvt(1:im,1:jm,:),kind_real)
 traj%tracers(isc:iec,jsc:jec,:,it_qi) = real(ltraj%qilst(1:im,1:jm,:) + ltraj%qicnt(1:im,1:jm,:),kind_real)
 traj%tracers(isc:iec,jsc:jec,:,it_ql) = real(ltraj%qllst(1:im,1:jm,:) + ltraj%qlcnt(1:im,1:jm,:),kind_real)

endsubroutine step_nl

! ------------------------------------------------------------------------------

subroutine step_tl(self,conf,traj,pert)

 implicit none

 class(fv3jedi_lm_moist_type), target, intent(inout) :: self
 type(fv3jedi_lm_conf), intent(in)    :: conf
 type(fv3jedi_lm_traj), intent(in)    :: traj
 type(fv3jedi_lm_pert), intent(inout) :: pert

 integer :: i,j,im,jm,lm,isc,iec,jsc,jec
 integer :: it_qv,it_qi,it_ql
 type(local_traj_moist), pointer :: ltraj
 type(local_pert_moist), pointer :: lpert
 type(local_cnst_moist), pointer :: lcnst

 im = conf%im
 jm = conf%jm
 lm = conf%lm
 isc = conf%isc
 iec = conf%iec
 jsc = conf%jsc
 jec = conf%jec

 !Convenience pointers
 lpert => self%lpert
 lcnst => self%lcnst
 if (conf%saveltraj) then
   ltraj => self%ltraj(conf%n)
 else
   ltraj => self%ltraj(1)
 endif

 !Set up the local trajectory
 if (.not. ltraj%set) call set_ltraj(conf,lcnst,traj,ltraj)

 !Local pert
 lpert%up(1:im,1:jm,:)    = dble(pert%u (isc:iec,jsc:jec,:))
 lpert%vp(1:im,1:jm,:)    = dble(pert%v (isc:iec,jsc:jec,:))
 lpert%ptp(1:im,1:jm,:)   = dble(pert%t (isc:iec,jsc:jec,:)) * p00**kappa / ltraj%pk(1:im,1:jm,:)

 call get_tracer_index(traj, 'water_vapor_mixing_ratio_wrt_moist_air', it_qv)
 call get_tracer_index(traj, 'cloud_liquid_ice', it_qi)
 call get_tracer_index(traj, 'cloud_liquid_water', it_ql)

 lpert%qvp(1:im,1:jm,:)   = dble(pert%tracers(isc:iec,jsc:jec,:,it_qv))
 lpert%qilsp(1:im,1:jm,:) = dble(pert%tracers(isc:iec,jsc:jec,:,it_qi)) * ltraj%ilsf(1:im,1:jm,:)
 lpert%qicnp(1:im,1:jm,:) = dble(pert%tracers(isc:iec,jsc:jec,:,it_qi)) * ltraj%icnf(1:im,1:jm,:)
 lpert%qllsp(1:im,1:jm,:) = dble(pert%tracers(isc:iec,jsc:jec,:,it_ql)) * ltraj%llsf(1:im,1:jm,:)
 lpert%qlcnp(1:im,1:jm,:) = dble(pert%tracers(isc:iec,jsc:jec,:,it_ql)) * ltraj%lcnf(1:im,1:jm,:)
 lpert%cflsp(1:im,1:jm,:) = 0.0_8
 lpert%cfcnp(1:im,1:jm,:) = dble(pert%cfcn(isc:iec,jsc:jec,:))

 ltraj%cnv_dqldtt  = 0.0_8
 ltraj%cnv_mfdt    = 0.0_8
 ltraj%cnv_prc3t   = 0.0_8
 ltraj%cnv_updft   = 0.0_8
 lpert%cnv_dqldtp  = 0.0_8
 lpert%cnv_mfdp    = 0.0_8
 lpert%cnv_prc3p   = 0.0_8
 lpert%cnv_updfp   = 0.0_8

 !Call the tangent linear convection scheme
 do i = 1,conf%im
   do j = 1,conf%jm

     if ( ltraj%doconvec(i,j) == 1) then
       call rase_d( 1, 1, conf%lm, lcnst%icmin, conf%dt,              &
                    lcnst%mapl8_cp, lcnst%mapl8_alhl,                 &
                    lcnst%mapl8_grav, lcnst%mapl8_rgas,               &
                    lcnst%mapl8_h2omw, lcnst%mapl8_airmw,             &
                    lcnst%mapl8_vireps,                               &
                    ltraj%seedras(i,j), lcnst%sige,                   &
                    ltraj%kcbl(i,j),                                  &
                    ltraj%wgt0(i,j,:), ltraj%wgt1(i,j,:),             &
                    ltraj%frland(i,j), ltraj%ts(i,j),                 &
                    ltraj%ptt(i,j,:), lpert%ptp(i,j,:),               &
                    ltraj%qvt(i,j,:), lpert%qvp(i,j,:),               &
                    ltraj%ut(i,j,:), lpert%up(i,j,:),                 &
                    ltraj%vt(i,j,:), lpert%vp(i,j,:),                 &
                    ltraj%co_auto(i,j), ltraj%cnv_ple(i,j,:),         &
                    ltraj%cnv_dqldtt(i,j,:), lpert%cnv_dqldtp(i,j,:), &
                    ltraj%cnv_mfdt(i,j,:),   lpert%cnv_mfdp(i,j,:),   &
                    ltraj%cnv_prc3t(i,j,:),  lpert%cnv_prc3p(i,j,:),  &
                    ltraj%cnv_updft(i,j,:),  lpert%cnv_updfp(i,j,:),  &
                    lcnst%rasparams, lcnst%estblx                     )
     endif

   enddo
 enddo

 !Call the tangent linear cloud scheme.
 call cloud_driver_d ( conf%dt, conf%im, conf%jm, conf%lm,                                         &
                       ltraj%ptt_c, lpert%ptp,                                                     &
                       ltraj%qvt_c, lpert%qvp,                                                     &
                       ltraj%ple,                                                                  &
                       ltraj%cnv_dqldtt_c, lpert%cnv_dqldtp, ltraj%cnv_mfdt_c,  lpert%cnv_mfdp,    &
                       ltraj%cnv_prc3t_c,  lpert%cnv_prc3p,  ltraj%cnv_updft_c, lpert%cnv_updfp,   &
                       ltraj%qilst, lpert%qilsp, ltraj%qllst, lpert%qllsp,                         &
                       ltraj%qicnt, lpert%qicnp, ltraj%qlcnt, lpert%qlcnp,                         &
                       ltraj%cflst, lpert%cflsp, ltraj%cfcnt, lpert%cfcnp,                         &
                       ltraj%frland, lcnst%cloudparams, lcnst%estblx, ltraj%khu, ltraj%khl,        &
                       lcnst%mapl8_runiv, lcnst%mapl8_kappa, lcnst%mapl8_airmw, lcnst%mapl8_h2omw, &
                       lcnst%mapl8_grav, lcnst%mapl8_alhl, lcnst%mapl8_alhf,   lcnst%mapl8_pi,     &
                       lcnst%mapl8_rgas, lcnst%mapl8_cp,   lcnst%mapl8_vireps, lcnst%mapl8_alhs,   &
                       lcnst%mapl8_tice, lcnst%mapl8_rvap, lcnst%mapl8_p00, conf%do_phy_mst        )

 !Back to pert
 pert%u(isc:iec,jsc:jec,:)    = real(lpert%up (1:im,1:jm,:),kind_real)
 pert%v(isc:iec,jsc:jec,:)    = real(lpert%vp (1:im,1:jm,:),kind_real)
 pert%t(isc:iec,jsc:jec,:)    = real(lpert%ptp(1:im,1:jm,:) * ltraj%pk(1:im,1:jm,:) / p00**kappa,kind_real)
 pert%cfcn(isc:iec,jsc:jec,:) = real(lpert%cfcnp(1:im,1:jm,:),kind_real)

 ! order of tracers is the same for traj and pert
 call get_tracer_index(traj, 'water_vapor_mixing_ratio_wrt_moist_air',  it_qv)
 call get_tracer_index(traj, 'cloud_liquid_ice',   it_qi)
 call get_tracer_index(traj, 'cloud_liquid_water', it_ql)

 pert%tracers(isc:iec,jsc:jec,:,it_qv) = real(lpert%qvp(1:im,1:jm,:),kind_real)
 pert%tracers(isc:iec,jsc:jec,:,it_qi) = real(lpert%qilsp(1:im,1:jm,:) + lpert%qicnp(1:im,1:jm,:),kind_real)
 pert%tracers(isc:iec,jsc:jec,:,it_ql) = real(lpert%qllsp(1:im,1:jm,:) + lpert%qlcnp(1:im,1:jm,:),kind_real)

endsubroutine step_tl

! ------------------------------------------------------------------------------

subroutine step_ad(self,conf,traj,pert)

 implicit none

 class(fv3jedi_lm_moist_type), target, intent(inout) :: self
 type(fv3jedi_lm_conf), intent(in)    :: conf
 type(fv3jedi_lm_traj), intent(in)    :: traj
 type(fv3jedi_lm_pert), intent(inout) :: pert

 integer :: i,j,im,jm,lm,isc,iec,jsc,jec
 integer :: it_qv,it_qi,it_ql
 type(local_traj_moist), pointer :: ltraj
 type(local_pert_moist), pointer :: lpert
 type(local_cnst_moist), pointer :: lcnst

 im = conf%im
 jm = conf%jm
 lm = conf%lm
 isc = conf%isc
 iec = conf%iec
 jsc = conf%jsc
 jec = conf%jec

 !Convenience pointers
 lpert => self%lpert
 lcnst => self%lcnst
 if (conf%saveltraj) then
   ltraj => self%ltraj(conf%n)
 else
   ltraj => self%ltraj(1)
 endif

 !Set up the local trajectory
 if (.not. ltraj%set) call set_ltraj(conf,lcnst,traj,ltraj)

 !Local pert
 lpert%up(1:im,1:jm,:)    = dble(pert%u (isc:iec,jsc:jec,:))
 lpert%vp(1:im,1:jm,:)    = dble(pert%v (isc:iec,jsc:jec,:))
 lpert%ptp(1:im,1:jm,:)   = dble(pert%t (isc:iec,jsc:jec,:)) * ltraj%pk(1:im,1:jm,:) / p00**kappa
 lpert%cflsp(1:im,1:jm,:) = 0.0_8
 lpert%cfcnp(1:im,1:jm,:) = dble(pert%cfcn(isc:iec,jsc:jec,:))

 call get_tracer_index(traj, 'water_vapor_mixing_ratio_wrt_moist_air', it_qv)
 call get_tracer_index(traj, 'cloud_liquid_ice', it_qi)
 call get_tracer_index(traj, 'cloud_liquid_water', it_ql)
 lpert%qvp(1:im,1:jm,:)   = dble(pert%tracers(isc:iec,jsc:jec,:,it_qv))
 lpert%qilsp(1:im,1:jm,:) = dble(pert%tracers(isc:iec,jsc:jec,:,it_qi))
 lpert%qicnp(1:im,1:jm,:) = dble(pert%tracers(isc:iec,jsc:jec,:,it_qi))
 lpert%qllsp(1:im,1:jm,:) = dble(pert%tracers(isc:iec,jsc:jec,:, it_ql))
 lpert%qlcnp(1:im,1:jm,:) = dble(pert%tracers(isc:iec,jsc:jec,:, it_ql))

 ltraj%cnv_dqldtt  = 0.0_8
 ltraj%cnv_mfdt    = 0.0_8
 ltraj%cnv_prc3t   = 0.0_8
 ltraj%cnv_updft   = 0.0_8
 lpert%cnv_dqldtp  = 0.0_8
 lpert%cnv_mfdp    = 0.0_8
 lpert%cnv_prc3p   = 0.0_8
 lpert%cnv_updfp   = 0.0_8

 !Call the tangent linear cloud scheme.
 call cloud_driver_b ( conf%dt, conf%im, conf%jm, conf%lm,                                         &
                       ltraj%ptt_c, lpert%ptp,                                                     &
                       ltraj%qvt_c, lpert%qvp,                                                     &
                       ltraj%ple,                                                                  &
                       ltraj%cnv_dqldtt_c, lpert%cnv_dqldtp, ltraj%cnv_mfdt_c,  lpert%cnv_mfdp,    &
                       ltraj%cnv_prc3t_c,  lpert%cnv_prc3p,  ltraj%cnv_updft_c, lpert%cnv_updfp,   &
                       ltraj%qilst, lpert%qilsp, ltraj%qllst, lpert%qllsp,                         &
                       ltraj%qicnt, lpert%qicnp, ltraj%qlcnt, lpert%qlcnp,                         &
                       ltraj%cflst, lpert%cflsp, ltraj%cfcnt, lpert%cfcnp,                         &
                       ltraj%frland, lcnst%cloudparams, lcnst%estblx, ltraj%khu, ltraj%khl,        &
                       lcnst%mapl8_runiv, lcnst%mapl8_kappa, lcnst%mapl8_airmw, lcnst%mapl8_h2omw, &
                       lcnst%mapl8_grav, lcnst%mapl8_alhl, lcnst%mapl8_alhf,   lcnst%mapl8_pi,     &
                       lcnst%mapl8_rgas, lcnst%mapl8_cp,   lcnst%mapl8_vireps, lcnst%mapl8_alhs,   &
                       lcnst%mapl8_tice, lcnst%mapl8_rvap, lcnst%mapl8_p00, conf%do_phy_mst        )

 !Call the tangent linear convection scheme
 do i = 1,conf%im
   do j = 1,conf%jm

     if ( ltraj%doconvec(i,j) == 1) then
       call rase_b( 1, 1, conf%lm, lcnst%icmin, conf%dt,              &
                    lcnst%mapl8_cp, lcnst%mapl8_alhl,                 &
                    lcnst%mapl8_grav, lcnst%mapl8_rgas,               &
                    lcnst%mapl8_h2omw, lcnst%mapl8_airmw,             &
                    lcnst%mapl8_vireps,                               &
                    ltraj%seedras(i,j), lcnst%sige,                   &
                    ltraj%kcbl(i,j),                                  &
                    ltraj%wgt0(i,j,:), ltraj%wgt1(i,j,:),             &
                    ltraj%frland(i,j), ltraj%ts(i,j),                 &
                    ltraj%ptt(i,j,:), lpert%ptp(i,j,:),               &
                    ltraj%qvt(i,j,:), lpert%qvp(i,j,:),               &
                    ltraj%ut(i,j,:), lpert%up(i,j,:),                 &
                    ltraj%vt(i,j,:), lpert%vp(i,j,:),                 &
                    ltraj%co_auto(i,j), ltraj%cnv_ple(i,j,:),         &
                    ltraj%cnv_dqldtt(i,j,:), lpert%cnv_dqldtp(i,j,:), &
                    ltraj%cnv_mfdt(i,j,:),   lpert%cnv_mfdp(i,j,:),   &
                    ltraj%cnv_prc3t(i,j,:),  lpert%cnv_prc3p(i,j,:),  &
                    ltraj%cnv_updft(i,j,:),  lpert%cnv_updfp(i,j,:),  &
                    lcnst%rasparams, lcnst%estblx                     )
     endif

   enddo
 enddo

 !Back to pert
 pert%u   (isc:iec,jsc:jec,:) = real(lpert%up (1:im,1:jm,:),kind_real)
 pert%v   (isc:iec,jsc:jec,:) = real(lpert%vp (1:im,1:jm,:),kind_real)
 pert%t   (isc:iec,jsc:jec,:) = real(lpert%ptp(1:im,1:jm,:) * p00**kappa,kind_real) / ltraj%pk(1:im,1:jm,:)

 call get_tracer_index(traj, 'water_vapor_mixing_ratio_wrt_moist_air', it_qv)
 call get_tracer_index(traj, 'cloud_liquid_ice', it_qi)
 call get_tracer_index(traj, 'cloud_liquid_water', it_ql)
 pert%tracers  (isc:iec,jsc:jec,:,it_qv) = real(lpert%qvp(1:im,1:jm,:),kind_real)
 pert%tracers  (isc:iec,jsc:jec,:,it_qi) = real(lpert%qilsp(1:im,1:jm,:)*ltraj%ilsf(1:im,1:jm,:) + &
                                                lpert%qicnp(1:im,1:jm,:)*ltraj%icnf(1:im,1:jm,:),kind_real)
 pert%tracers  (isc:iec,jsc:jec,:,it_ql) = real(lpert%qllsp(1:im,1:jm,:)*ltraj%llsf(1:im,1:jm,:) + &
                                                lpert%qlcnp(1:im,1:jm,:)*ltraj%lcnf(1:im,1:jm,:),kind_real)
 pert%cfcn(isc:iec,jsc:jec,:) = real(lpert%cfcnp(1:im,1:jm,:),kind_real)

endsubroutine step_ad

! ------------------------------------------------------------------------------

subroutine delete(self,conf)

 implicit none
 class(fv3jedi_lm_moist_type), intent(inout) :: self
 type(fv3jedi_lm_conf), intent(in)    :: conf

 integer :: n

 deallocate(self%lcnst%ESTBLX)
 deallocate(self%lcnst%SIGE)

 call deallocate_lpert(self%lpert)

 if (conf%saveltraj) then
   do n = 1,conf%nt
     call deallocate_ltraj(self%ltraj(n))
   enddo
 else
   call deallocate_ltraj(self%ltraj(1))
 endif

 deallocate(self%ltraj)

endsubroutine delete

! ------------------------------------------------------------------------------

subroutine set_ltraj(conf,lcnst,traj,ltraj)

 type(fv3jedi_lm_conf),         intent(in)    :: conf
 type(local_cnst_moist),        intent(in)    :: lcnst
 type(fv3jedi_lm_traj), target, intent(in)    :: traj
 type(local_traj_moist),        intent(inout) :: ltraj

 real(8), allocatable, dimension(:,:,:) :: PLO, PK
 real(8), allocatable, dimension(:,:,:) :: TEMP
 real(8), allocatable, dimension(:,:,:) :: PTT_F, QVT_F
 real(8), allocatable, dimension(:,:,:) :: PTT_L, QVT_L
 real(8), allocatable, dimension(:,:,:) :: HEAT
 integer, allocatable, dimension(:,:)   :: CTOP
 real(8), allocatable, dimension(:,:)   :: sumHEAT
 real(8), allocatable, dimension(:,:)   :: JACOBIAN
 real(8), allocatable, dimension(:)     :: PT_pert, QV_pert
 real(8), allocatable, dimension(:)     :: PT_pert_in, QV_pert_in
 real(8), allocatable, dimension(:)     :: H_pert, M_pert
 real(8), allocatable, dimension(:,:,:) :: fQi

 integer :: i,j,l,im,jm,lm,isc,iec,jsc,jec,i_qv
 real(8), parameter :: PMIN_DET = 3000.0, AUTOC_CN_OCN  = 2.5e-3, AUTOC_CN_LAND = AUTOC_CN_OCN
 integer :: maxcondep

 im  = conf%im
 jm  = conf%jm
 lm  = conf%lm
 isc = conf%isc
 iec = conf%iec
 jsc = conf%jsc
 jec = conf%jec

 allocate(plo(im,jm,lm)   )
 allocate(pk(im,jm,lm)    )
 allocate(temp(im,jm,lm)  )
 allocate(ptt_f(im,jm,lm) )
 allocate(qvt_f(im,jm,lm) )
 allocate(ptt_l(im,jm,lm) )
 allocate(qvt_l(im,jm,lm) )
 allocate(heat(im,jm,lm)  )
 allocate(ctop(im,jm)     )
 allocate(sumheat(im,jm)  )
 allocate(jacobian(2*lm,2))
 allocate(h_pert(lm)      )
 allocate(m_pert(lm)      )
 allocate(pt_pert(lm)     )
 allocate(qv_pert(lm)     )
 allocate(pt_pert_in(lm)  )
 allocate(qv_pert_in(lm)  )
 allocate(fqi(im,jm,lm)   )

 !!Still on the D-Grid!!
 ltraj%ut(1:im,1:jm,:) = dble(traj%u(isc:iec,jsc:jec,:))
 ltraj%vt(1:im,1:jm,:) = dble(traj%v(isc:iec,jsc:jec,:))

 !Pressure hPa, half levels and p^kappa
 call compute_pressures(im,jm,lm,conf%ptop,traj%delp(isc:iec,jsc:jec,:),ltraj%ple,plo,ltraj%pk)

 !Potential temperature from temperature (use proper pk form to get pt)
 ltraj%PTT(1:im,1:jm,:) = p00**kappa * dble(traj%t(isc:iec,jsc:jec,:)) / ltraj%pk(1:im,1:jm,:)

 !Pressures in form used be GEOS moist physics
 ltraj%cnv_ple     = 0.01_8*ltraj%ple
 PLO               = 0.5_8*(ltraj%CNV_PLE(:,:,0:LM-1) +  ltraj%CNV_PLE(:,:,1:LM  ) )
 PK                = (PLO/1000.0_8)**(lcnst%MAPL8_RGAS/lcnst%MAPL8_CP) !Formulation in GEOS moist
 TEMP              = ltraj%PTT*PK

 !Some moist vars
 call get_tracer_index(traj, 'water_vapor_mixing_ratio_wrt_moist_air', i_qv)


 ltraj%qvt(1:im,1:jm,:) = dble(traj%tracers(isc:iec,jsc:jec,:,i_qv))
 ltraj%cflst  = 0.0_8

 ltraj%cfcnt(1:im,1:jm,:)  = dble(traj%cfcn(isc:iec,jsc:jec,:))

 ltraj%ts    (1:im,1:jm) = dble(traj%ts    (isc:iec,jsc:jec))
 ltraj%frland(1:im,1:jm) = dble(traj%frland(isc:iec,jsc:jec))
 ltraj%kcbl  (1:im,1:jm) = nint(traj%kcbl  (isc:iec,jsc:jec))
 ltraj%khl   (1:im,1:jm) = nint(traj%khl   (isc:iec,jsc:jec))
 ltraj%khu   (1:im,1:jm) = nint(traj%khu   (isc:iec,jsc:jec))

 ltraj%PTT_C = ltraj%PTT
 ltraj%QVT_C = ltraj%QVT

 ltraj%CNV_DQLDTT_C  = 0.0_8
 ltraj%CNV_MFDT_C    = 0.0_8
 ltraj%CNV_PRC3T_C   = 0.0_8
 ltraj%CNV_UPDFT_C   = 0.0_8

 PTT_F = ltraj%PTT
 QVT_F = ltraj%QVT

 PTT_L = ltraj%PTT
 QVT_L = ltraj%QVT

 !Not linearised as could produce unpredictable behaviour (and very sensitive).
 ltraj%SEEDRAS(:,:) = 1000000 * ( 100*TEMP(:,:,LM) - INT( 100*TEMP(:,:,LM) ) )

 !Strapping levels
 DO I = 1,IM
    DO J = 1,JM

       ltraj%WGT0(I,J,:)                    = 0.0_8
       ltraj%WGT0(I,J,ltraj%KCBL(I,J):LM)   = 1.0_8
       ltraj%WGT1(I,J,:)                    = 0.0_8
       ltraj%WGT1(I,J,ltraj%KCBL(I,J):LM)   = 1.0_8
    ENDDO
 ENDDO

 where (ltraj%FRLAND<0.1_8)
    ltraj%CO_AUTO = AUTOC_CN_OCN   ! ocean value
 elsewhere
    ltraj%CO_AUTO = AUTOC_CN_LAND  ! land value
 end where

 !Call nonlinear convection scheme. We need to do this because the
 !cloud adjoint scheme needs the outputs from the convection.
 !We will also only call the convective linearizations for profiles
 !where convection is occuring for efficiency.
 CALL RASE0(IM*JM, IM*JM, LM, lcnst%ICMIN, conf%DT,                   &
            lcnst%MAPL8_CP, lcnst%MAPL8_ALHL, lcnst%MAPL8_GRAV, lcnst%MAPL8_RGAS,  &
            lcnst%MAPL8_H2OMW, lcnst%MAPL8_AIRMW, lcnst%MAPL8_VIREPS,        &
            ltraj%SEEDRAS, lcnst%SIGE,                                 &
            ltraj%KCBL, ltraj%WGT0, ltraj%WGT1, ltraj%FRLAND, ltraj%TS,                  &
            ltraj%PTT_C, ltraj%QVT_C,                                  &
            ltraj%CO_AUTO, ltraj%CNV_PLE,                              &
            ltraj%CNV_DQLDTT_C, ltraj%CNV_MFDT_C,                      &
            ltraj%CNV_PRC3T_C, ltraj%CNV_UPDFT_C,                      &
            lcnst%RASPARAMS, lcnst%ESTBLX                              )

 ! Do the filtering to determine whether linear convection should be called
 ! ------------------------------------------------------------------------
 !Figure out whether or not each convective profile should be linearized.
 !Conditions  - Convection is happening (heating rate is nonzero)
 !            - Convection is deep enough (>= 10 levels)
 !            - The heating rate profile is not just a spike at one level
 !            - Gradients are not too steep (Jacobian filtering).
 ltraj%DOCONVEC = 0
 HEAT = 0.0
 CTOP = LM
 sumHEAT = 0.0

 !Be more lenient on profiles let in for 4DVAR, less likely to encounter problems
 !at shorter lead times and gives more realistic low level cloud perturbations.
 if (conf%do_phy_mst == 1) then
    MAXCONDEP = 1
 elseif (conf%do_phy_mst == 2) then
    MAXCONDEP = 10
 endif

 DO I = 1,IM
    DO J = 1,JM

       !Compute the heating rate.
       HEAT(I,J,:) = (ltraj%PTT_C(I,J,:) - ltraj%PTT(I,J,:))/conf%DT

       !Starting at the top scan downwards to look for nonzero heating rate,
       !record index of highest level convection reaches (ignoring v.small heating rate)
       DO L = 1,LM
          IF ( abs(HEAT(I,J,L)) .gt. 0.01*maxval(abs(HEAT(I,J,:)),1) ) THEN
             CTOP(I,J) = L
             exit
          ENDIF
       ENDDO

       !Compute sort of integral of the heating rate.
       if ( (CTOP(I,J) .ne. LM) .and. (ltraj%KCBL(I,J) - CTOP(I,J) > 0) ) then
          sumHEAT(I,J) = (   sum(abs(HEAT(I,J,CTOP(I,J):ltraj%KCBL(I,J)-1))) &
                                 - maxval(abs(HEAT(I,J,CTOP(I,J):ltraj%KCBL(I,J)-1)),1) ) &
                                 / ( ltraj%KCBL(I,J) - CTOP(I,J)  )
       endif

       !Compare `integral` to maximum absolute heating rate
       IF ( ltraj%KCBL(I,J) - CTOP(I,J) >= MAXCONDEP ) THEN
          IF (sumHEAT(I,J) / maxval(abs(HEAT(I,J,1:ltraj%KCBL(I,J)-1)),1) > 0.125 ) then
             ltraj%DOCONVEC(I,J) = 1
          endif
       endif

       !Compute two columns of the Jacobian to check for steep gradients
       !This prevents instability and floating point issues that can cause failure of the TLM/ADJ dot prod test
       IF ( ltraj%DOCONVEC(I,J) == 1 ) THEN
          call jacobian_filter_tlm
       ENDIF

    enddo
 enddo

 !Prepare the inputs for the cloud scheme
 !---------------------------------------
 !Compute the ice fraction for each grid box
 DO i = 1,IM
    DO j = 1,JM
       DO l = 1,LM
          call IceFraction( TEMP(i,j,l), fQi(i,j,l) )
       enddo
    enddo
 enddo

 !Split the input large scale and convective cloud into ice and liquid parts.
 ltraj%QILST(1:im,1:jm,:) = dble(traj%QLS(isc:iec,jsc:jec,:)) * fQi(1:im,1:jm,:)
 ltraj%QLLST(1:im,1:jm,:) = dble(traj%QLS(isc:iec,jsc:jec,:)) * (1-fQi(1:im,1:jm,:))
 ltraj%QICNT(1:im,1:jm,:) = dble(traj%QCN(isc:iec,jsc:jec,:)) * fQi(1:im,1:jm,:)
 ltraj%QLCNT(1:im,1:jm,:) = dble(traj%QCN(isc:iec,jsc:jec,:)) * (1-fQi(1:im,1:jm,:))

 !Split the perturbations for total cloud water and ice into the convective and large scale parts.
 !Spilitting is based on fraction of total cloud ice/water that is attributed to large scale and anvil
 !in the trajectory at this time step. Note that we don't really need to protect for small values of
 !total cloud sum, if the denominator is small then so is the numerator, though compiler may complain.
 ltraj%ILSF = 0.0
 ltraj%ICNF = 0.0
 ltraj%LLSF = 0.0
 ltraj%LCNF = 0.0
 DO I = 1,IM
    DO J = 1,JM
       DO L = 1,LM

          if ( ltraj%QILST(i,j,l) + ltraj%QICNT(i,j,l) .gt. 0.0_8 ) then
             ltraj%ILSF(i,j,l) = ltraj%QILST(i,j,l) / ( ltraj%QILST(i,j,l) + ltraj%QICNT(i,j,l) )
             ltraj%ICNF(i,j,l) = ltraj%QICNT(i,j,l) / ( ltraj%QILST(i,j,l) + ltraj%QICNT(i,j,l) )
          endif
          if ( ltraj%QLLST(i,j,l) + ltraj%QLCNT(i,j,l) .gt. 0.0_8 ) then
             ltraj%LLSF(i,j,l) = ltraj%QLLST(i,j,l) / ( ltraj%QLLST(i,j,l) + ltraj%QLCNT(i,j,l) )
             ltraj%LCNF(i,j,l) = ltraj%QLCNT(i,j,l) / ( ltraj%QLLST(i,j,l) + ltraj%QLCNT(i,j,l) )
          endif

       enddo
    enddo
 enddo

 deallocate(plo)
 deallocate(pk)
 deallocate(temp)
 deallocate(ptt_f)
 deallocate(qvt_f)
 deallocate(ptt_l)
 deallocate(qvt_l)
 deallocate(heat)
 deallocate(ctop)
 deallocate(sumheat)
 deallocate(jacobian)
 deallocate(h_pert)
 deallocate(m_pert)
 deallocate(pt_pert)
 deallocate(qv_pert)
 deallocate(pt_pert_in)
 deallocate(qv_pert_in)
 deallocate(fqi)

 contains

 subroutine jacobian_filter_tlm

  !Compute two specific columns of the Jacobian for a single profile. Then filter if values in those
  !columns are over a certain value, implying too steep gradients.

  JACOBIAN = 0.0_8

  DO L = 1,1 !Perturb level by level

     PT_pert = 0.0_8
     QV_pert = 0.0_8

     if (L == 1) then

        PT_pert(ltraj%KCBL(I,J)) = 1.0_8

     elseif (L == 2) then

        if (ltraj%KCBL(I,J) == lm) then
          QV_pert(ltraj%KCBL(I,J)) = 1.0_8
        else
           QV_pert(ltraj%KCBL(I,J) + 1) = 1.0_8
        endif

     endif

     !Save precall prognostic variables
     PT_pert_in = PT_pert
     QV_pert_in = QV_pert

     !Do nonlinear convection to create inputs trajectory for large scale adjoint
     call rase0_d ( 1, 1, lm, lcnst%icmin, conf%dt,          &
                    lcnst%mapl8_cp, lcnst%mapl8_alhl,       &
                    lcnst%mapl8_grav, lcnst%mapl8_rgas,     &
                    lcnst%mapl8_h2omw, lcnst%mapl8_airmw,   &
                    lcnst%mapl8_vireps,                      &
                    ltraj%seedras(i,j), lcnst%sige,          &
                    ltraj%kcbl(i,j),                          &
                    ltraj%wgt0(i,j,:), ltraj%wgt1(i,j,:),     &
                    ltraj%frland(i,j), ltraj%ts(i,j),         &
                    ptt_f(i,j,:), pt_pert,                    &
                    qvt_f(i,j,:), qv_pert,                    &
                    ltraj%co_auto(i,j), ltraj%cnv_ple(i,j,:), &
                    lcnst%rasparams, lcnst%estblx           )

     !Compute perturbation heating and moistening rates
     H_pert = (PT_pert - PT_pert_in)/conf%DT
     M_pert = (QV_pert - QV_pert_in)/conf%DT

     !Uncomment here if just doing two columns of the Jacobian
     if (L == 1) then
        Jacobian(0*LM+1:1*LM,1) = H_pert
        Jacobian(1*LM+1:2*LM,1) = M_pert
     elseif (L == 2) then
        Jacobian(0*LM+1:1*LM,2) = H_pert
        Jacobian(1*LM+1:2*LM,2) = M_pert
     endif

  endDO

  !Constants here determined so as to remove as many of the problematic points as possible.
  !The constants used in this if loop are tuned from looking at many Jacobians for many time steps. Values are choosen
  !so as to balance between keeping the natural behahiour for as many points as possible without
  if ( (maxval(abs(Jacobian(1:lm     ,1))) .gt. 0.00010  ) .or. &
       (maxval(abs(Jacobian(1:lm     ,2))) .gt. 0.25     ) .or. &
       (maxval(abs(Jacobian(lm+1:2*lm,1))) .gt. 1.0e-07  ) .or. &
       (maxval(abs(Jacobian(lm+1:2*lm,2))) .gt. 0.000250 ) ) then

     ltraj%doconvec(i,j) = 0

  else

     ltraj%doconvec(i,j) = 1

  endif

 endsubroutine jacobian_filter_tlm

endsubroutine set_ltraj

! ------------------------------------------------------------------------------

subroutine allocate_ltraj(im,jm,lm,ltraj)

 integer, intent(in) :: im,jm,lm
 type(local_traj_moist), intent(inout) :: ltraj

 !double precision trajectories
 allocate(ltraj%ut(im,jm,lm)          )
 allocate(ltraj%vt(im,jm,lm)          )
 allocate(ltraj%ptt(im,jm,lm)         )
 allocate(ltraj%qvt(im,jm,lm)         )
 allocate(ltraj%ts(im,jm)             )
 allocate(ltraj%frland(im,jm)         )
 allocate(ltraj%kcbl(im,jm)           )
 allocate(ltraj%khu(im,jm)            )
 allocate(ltraj%khl(im,jm)            )

 !pressure and temperature variables
 allocate(ltraj%ple(im,jm,0:lm)       )
 allocate(ltraj%cnv_ple(im,jm,0:lm)   )
 allocate(ltraj%pk(im,jm,lm)          )

 !convection varaibles
 allocate(ltraj%cnv_dqldtt(im,jm,lm)  )
 allocate(ltraj%cnv_mfdt(im,jm,lm)    )
 allocate(ltraj%cnv_prc3t(im,jm,lm)   )
 allocate(ltraj%cnv_updft(im,jm,lm)   )
 allocate(ltraj%ptt_c(im,jm,lm)       )
 allocate(ltraj%qvt_c(im,jm,lm)       )
 allocate(ltraj%cnv_dqldtt_c(im,jm,lm))
 allocate(ltraj%cnv_mfdt_c(im,jm,lm)  )
 allocate(ltraj%cnv_prc3t_c(im,jm,lm) )
 allocate(ltraj%cnv_updft_c(im,jm,lm) )
 allocate(ltraj%seedras(im,jm)        )
 allocate(ltraj%co_auto(im,jm)        )
 allocate(ltraj%doconvec(im,jm)       )
 allocate(ltraj%wgt0(im,jm,lm)        )
 allocate(ltraj%wgt1(im,jm,lm)        )

 !Cloud variables
 allocate(ltraj%qilst(im,jm,lm)       )
 allocate(ltraj%qllst(im,jm,lm)       )
 allocate(ltraj%qicnt(im,jm,lm)       )
 allocate(ltraj%qlcnt(im,jm,lm)       )
 allocate(ltraj%cflst(im,jm,lm)       )
 allocate(ltraj%cfcnt(im,jm,lm)       )
 allocate(ltraj%ilsf(im,jm,lm)        )
 allocate(ltraj%icnf(im,jm,lm)        )
 allocate(ltraj%llsf(im,jm,lm)        )
 allocate(ltraj%lcnf(im,jm,lm)        )

endsubroutine allocate_ltraj

! ------------------------------------------------------------------------------

subroutine allocate_lpert(im,jm,lm,lpert)

 integer, intent(in) :: im,jm,lm
 type(local_pert_moist), intent(inout) :: lpert

 allocate(lpert%up(im,jm,lm)          )
 allocate(lpert%vp(im,jm,lm)          )
 allocate(lpert%ptp(im,jm,lm)         )
 allocate(lpert%qvp(im,jm,lm)         )
 allocate(lpert%cnv_dqldtp(im,jm,lm)  )
 allocate(lpert%cnv_mfdp(im,jm,lm)    )
 allocate(lpert%cnv_prc3p(im,jm,lm)   )
 allocate(lpert%cnv_updfp(im,jm,lm)   )
 allocate(lpert%qilsp(im,jm,lm)       )
 allocate(lpert%qllsp(im,jm,lm)       )
 allocate(lpert%qicnp(im,jm,lm)       )
 allocate(lpert%qlcnp(im,jm,lm)       )
 allocate(lpert%cflsp(im,jm,lm)       )
 allocate(lpert%cfcnp(im,jm,lm)       )

endsubroutine allocate_lpert

! ------------------------------------------------------------------------------

subroutine deallocate_ltraj(ltraj)

 type(local_traj_moist), intent(inout) :: ltraj

 deallocate(ltraj%ut)
 deallocate(ltraj%vt)
 deallocate(ltraj%ptt)
 deallocate(ltraj%qvt)
 deallocate(ltraj%ts)
 deallocate(ltraj%frland)
 deallocate(ltraj%kcbl)
 deallocate(ltraj%khu)
 deallocate(ltraj%khl)
 deallocate(ltraj%ple)
 deallocate(ltraj%cnv_ple)
 deallocate(ltraj%pk)
 deallocate(ltraj%cnv_dqldtt)
 deallocate(ltraj%cnv_mfdt)
 deallocate(ltraj%cnv_prc3t)
 deallocate(ltraj%cnv_updft)
 deallocate(ltraj%ptt_c)
 deallocate(ltraj%qvt_c)
 deallocate(ltraj%cnv_dqldtt_c)
 deallocate(ltraj%cnv_mfdt_c)
 deallocate(ltraj%cnv_prc3t_c)
 deallocate(ltraj%cnv_updft_c)
 deallocate(ltraj%seedras)
 deallocate(ltraj%co_auto)
 deallocate(ltraj%doconvec)
 deallocate(ltraj%wgt0)
 deallocate(ltraj%wgt1)
 deallocate(ltraj%qilst)
 deallocate(ltraj%qllst)
 deallocate(ltraj%qicnt)
 deallocate(ltraj%qlcnt)
 deallocate(ltraj%cflst)
 deallocate(ltraj%cfcnt)
 deallocate(ltraj%ilsf)
 deallocate(ltraj%icnf)
 deallocate(ltraj%llsf)
 deallocate(ltraj%lcnf)

endsubroutine deallocate_ltraj

! ------------------------------------------------------------------------------

subroutine deallocate_lpert(lpert)

 type(local_pert_moist), intent(inout) :: lpert

 deallocate(lpert%up)
 deallocate(lpert%vp)
 deallocate(lpert%ptp)
 deallocate(lpert%qvp)
 deallocate(lpert%cnv_dqldtp)
 deallocate(lpert%cnv_mfdp)
 deallocate(lpert%cnv_prc3p)
 deallocate(lpert%cnv_updfp)
 deallocate(lpert%qilsp)
 deallocate(lpert%qllsp)
 deallocate(lpert%qicnp)
 deallocate(lpert%qlcnp)
 deallocate(lpert%cflsp)
 deallocate(lpert%cfcnp)

endsubroutine deallocate_lpert

! ------------------------------------------------------------------------------


end module fv3jedi_lm_moist_mod
