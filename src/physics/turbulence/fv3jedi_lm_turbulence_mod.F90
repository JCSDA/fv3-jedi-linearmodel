module fv3jedi_lm_turbulence_mod

use fv3jedi_lm_utils_mod
use fv3jedi_lm_kinds_mod
use fv3jedi_lm_const_mod

! BL Schemes
use bldriver

! Saturation table
use qsat_util

use tapenade_iter, only: cp_iter, cp_iter_controls, initialize_cp_iter, finalize_cp_iter
use tapenade_iter, only: cp_mod_ini, cp_mod_mid, cp_mod_end, pushrealarray, poprealarray

!> Turbulence module

implicit none
private
public :: fv3jedi_lm_turbulence_type

!> Local trajectory object
type local_traj_turbulence
  real(kind_real), allocatable, dimension(:,:,:) :: akq, bkq, ckq
  real(kind_real), allocatable, dimension(:,:,:) :: aks, bks, cks
  real(kind_real), allocatable, dimension(:,:,:) :: akv, bkv, ckv
  real(kind_real), allocatable, dimension(:,:,:) :: pk
  logical :: set = .false.
endtype local_traj_turbulence

!> Local constants object
type local_cnst_turbulence
  real(kind_real), dimension(22) :: TURBPARAMS
  integer,         dimension(4)  :: TURBPARAMSI
endtype local_cnst_turbulence

!> Turbulence class (self)
type fv3jedi_lm_turbulence_type
 type(local_traj_turbulence), allocatable :: ltraj(:)
 type(local_cnst_turbulence) :: lcnst
 contains
  procedure :: create
  procedure :: init_nl
  procedure :: init_tl
  procedure :: init_ad
  procedure :: step_nl
  procedure :: step_tl
  procedure :: step_ad
  procedure :: delete
end type fv3jedi_lm_turbulence_type

contains

! ------------------------------------------------------------------------------

subroutine create(self,conf)

 implicit none

 class(fv3jedi_lm_turbulence_type), target, intent(inout) :: self
 type(fv3jedi_lm_conf), intent(in)    :: conf

 real(kind_real), allocatable, dimension(:) :: pref
 integer :: l, n

 if (conf%saveltraj) then
   allocate(self%ltraj(conf%nt))
   do n = 1,conf%nt
     call allocate_ltraj(conf%im,conf%jm,conf%lm,self%ltraj(n))
   enddo
 else
   allocate(self%ltraj(1))
   call allocate_ltraj(conf%im,conf%jm,conf%lm,self%ltraj(1))
 endif

 allocate(pref(0:conf%lm))
 DO l = 0,conf%lm
   pref(l) = conf%AK(l+1) + conf%BK(l+1)*p00
 enddo

 !Turbulence Parameters
 self%lcnst%TurbParams(1)  = 5.0_kind_real                    !LOUIS
 self%lcnst%TurbParams(2)  = 160.0_kind_real                  !LAMBDAM
 self%lcnst%TurbParams(3)  = 1.0_kind_real                    !LAMBDAM2
 self%lcnst%TurbParams(4)  = 160.0_kind_real                  !LAMBDAH
 self%lcnst%TurbParams(5)  = 1.0_kind_real                    !LAMBDAH2
 self%lcnst%TurbParams(6)  = 3000.0_kind_real                 !ZKMENV
 self%lcnst%TurbParams(7)  = 3000.0_kind_real                 !ZKHENV
 self%lcnst%TurbParams(8)  = 0.1_kind_real                    !MINTHICK
 self%lcnst%TurbParams(9)  = 0.0030_kind_real                 !MINSHEAR
 self%lcnst%TurbParams(10) = 2.5101471e-8_kind_real           !C_B
 self%lcnst%TurbParams(11) = 1500.0_kind_real                 !LAMBDA_B
 self%lcnst%TurbParams(12) = 500.0_kind_real                  !AKHMMAX
 self%lcnst%TurbParams(13) = 1.0_kind_real                    !PRANDTLSFC
 self%lcnst%TurbParams(14) = 0.75_kind_real                   !PRANDTLRAD
 self%lcnst%TurbParams(15) = 0.50_kind_real                   !BETA_RAD
 self%lcnst%TurbParams(16) = 0.25_kind_real                   !BETA_SURF
 self%lcnst%TurbParams(17) = 0.85_kind_real                   !KHRADFAC
 self%lcnst%TurbParams(18) = 0.45_kind_real                   !KHSFCFAC
 self%lcnst%TurbParams(19) = 20.0_kind_real                   !TPFAC_SURF
 self%lcnst%TurbParams(20) = 1.5e-3_kind_real                 !ENTRATE_SURF
 self%lcnst%TurbParams(21) = 0.5_kind_real                    !PCEFF_SURF
 self%lcnst%TurbParams(22) = -999.0_kind_real                 !LOUIS_MEMORY
 self%lcnst%TurbParamsI(1) = count(pref < 50000.0_kind_real)  !KPBLMIN
 self%lcnst%TurbParamsI(2) = 1                                !LOCK_ON
 self%lcnst%TurbParamsI(3) = 1                                !PBLHT_OPTION
 self%lcnst%TurbParamsI(4) = 0                                !RADLW_DEP

 deallocate(pref)

endsubroutine create

! ------------------------------------------------------------------------------

subroutine init_nl(self,pert,traj)

 implicit none

 class(fv3jedi_lm_turbulence_type), intent(inout) :: self
 type(fv3jedi_lm_pert), intent(inout) :: pert
 type(fv3jedi_lm_traj), intent(in) :: traj

endsubroutine init_nl

! ------------------------------------------------------------------------------

subroutine init_tl(self,pert,traj)

 implicit none

 class(fv3jedi_lm_turbulence_type), intent(inout) :: self
 type(fv3jedi_lm_pert), intent(inout) :: pert
 type(fv3jedi_lm_traj), intent(in) :: traj

endsubroutine init_tl

! ------------------------------------------------------------------------------

subroutine init_ad(self,pert,traj)

 implicit none

 class(fv3jedi_lm_turbulence_type), intent(inout) :: self
 type(fv3jedi_lm_pert), intent(inout) :: pert
 type(fv3jedi_lm_traj), intent(in) :: traj

endsubroutine init_ad

! ------------------------------------------------------------------------------

subroutine step_nl(self,conf,traj)

 implicit none

 integer :: f

 class(fv3jedi_lm_turbulence_type), target, intent(inout) :: self
 type(fv3jedi_lm_traj), target,             intent(inout) :: traj
 type(fv3jedi_lm_conf),                     intent(in)    :: conf

 type(local_traj_turbulence), pointer :: ltraj
 real(kind_real), pointer, dimension(:,:,:) :: p_u
 real(kind_real), pointer, dimension(:,:,:) :: p_v
 real(kind_real), pointer, dimension(:,:,:) :: p_t
 real(kind_real), pointer, dimension(:,:,:) :: p_delp

 !Pointers with ind starting at 1
 p_u   (1:,1:,1:) => traj%u
 p_v   (1:,1:,1:) => traj%v
 p_t   (1:,1:,1:) => traj%t
 p_delp(1:,1:,1:) => traj%delp

 !Convenience pointers
 if (conf%saveltraj) then
   ltraj => self%ltraj(conf%n)
 else
   ltraj => self%ltraj(1)
 endif

 !Set up the local trajectory
 if (.not. ltraj%set) call set_ltraj(conf,self%lcnst,traj,ltraj)

 !t2pt
 p_t = p00**kappa * p_t / ltraj%pk

 call vtrisolvepert(conf%im,conf%jm,conf%lm,ltraj%akv,ltraj%bkv,ltraj%ckv,p_u ,1,1)
 call vtrisolvepert(conf%im,conf%jm,conf%lm,ltraj%akv,ltraj%bkv,ltraj%ckv,p_v ,1,1)
 call vtrisolvepert(conf%im,conf%jm,conf%lm,ltraj%aks,ltraj%bks,ltraj%cks,p_t ,1,1)
 do f = 1, size(traj%tracer_names)
   call vtrisolvepert(conf%im,conf%jm,conf%lm,ltraj%akq,ltraj%bkq,ltraj%ckq,traj%tracers(:,:,:,f),1,0)
 end do

 !pt2t
 p_t = ltraj%pk * p_t / p00**kappa

 !Nullify
 nullify(p_u)
 nullify(p_v)
 nullify(p_t)
 nullify(p_delp)
 nullify(ltraj)

endsubroutine step_nl

! ------------------------------------------------------------------------------

subroutine step_tl(self,conf,traj,pert)

 implicit none

 integer :: f

 class(fv3jedi_lm_turbulence_type), target, intent(inout) :: self
 type(fv3jedi_lm_conf),                     intent(in)    :: conf
 type(fv3jedi_lm_traj),                     intent(in)    :: traj
 type(fv3jedi_lm_pert), target,             intent(inout) :: pert

 type(local_traj_turbulence), pointer :: ltraj
 real(kind_real), pointer, dimension(:,:,:) :: p_u
 real(kind_real), pointer, dimension(:,:,:) :: p_v
 real(kind_real), pointer, dimension(:,:,:) :: p_t
 real(kind_real), pointer, dimension(:,:,:) :: p_delp

 !Pointers with ind starting at 1
 p_u   (1:,1:,1:) => pert%u
 p_v   (1:,1:,1:) => pert%v
 p_t   (1:,1:,1:) => pert%t
 p_delp(1:,1:,1:) => pert%delp

 !Convenience pointers
 if (conf%saveltraj) then
   ltraj => self%ltraj(conf%n)
 else
   ltraj => self%ltraj(1)
 endif

 !Set up the local trajectory
 if (.not. ltraj%set) call set_ltraj(conf,self%lcnst,traj,ltraj)

 !t2pt
 p_t = p00**kappa * p_t / ltraj%pk

 call vtrisolvepert(conf%im,conf%jm,conf%lm,ltraj%akv,ltraj%bkv,ltraj%ckv,p_u ,1,1)
 call vtrisolvepert(conf%im,conf%jm,conf%lm,ltraj%akv,ltraj%bkv,ltraj%ckv,p_v ,1,1)
 call vtrisolvepert(conf%im,conf%jm,conf%lm,ltraj%aks,ltraj%bks,ltraj%cks,p_t ,1,1)
 call vtrisolvepert(conf%im,conf%jm,conf%lm,ltraj%akq,ltraj%bkq,ltraj%ckq,pert%tracers(:,:,:,1),1,1)

 do f = 2, size(traj%tracer_names)
   call vtrisolvepert(conf%im,conf%jm,conf%lm,ltraj%akq,ltraj%bkq,ltraj%ckq,pert%tracers(:,:,:,f),1,0)
 end do

 !pt2t
 p_t = ltraj%pk * p_t / p00**kappa

 !Nullify
 nullify(p_u)
 nullify(p_v)
 nullify(p_t)
 nullify(p_delp)
 nullify(ltraj)

endsubroutine step_tl

! ------------------------------------------------------------------------------

subroutine step_ad(self,conf,traj,pert)

 implicit none

 integer :: f

 class(fv3jedi_lm_turbulence_type), target, intent(inout) :: self
 type(fv3jedi_lm_conf),                     intent(in)    :: conf
 type(fv3jedi_lm_traj),                     intent(in)    :: traj
 type(fv3jedi_lm_pert), target,             intent(inout) :: pert

 type(local_traj_turbulence), pointer :: ltraj
 real(kind_real), pointer, dimension(:,:,:) :: p_u
 real(kind_real), pointer, dimension(:,:,:) :: p_v
 real(kind_real), pointer, dimension(:,:,:) :: p_t
 real(kind_real), pointer, dimension(:,:,:) :: p_delp

 !Pointers with ind starting at 1
 p_u   (1:,1:,1:) => pert%u
 p_v   (1:,1:,1:) => pert%v
 p_t   (1:,1:,1:) => pert%t
 p_delp(1:,1:,1:) => pert%delp

 !Convenience pointers
 if (conf%saveltraj) then
   ltraj => self%ltraj(conf%n)
 else
   ltraj => self%ltraj(1)
 endif

 !Set up the local trajectory
 if (.not. ltraj%set) call set_ltraj(conf,self%lcnst,traj,ltraj)

 !t2pt adjoint
 p_t = ltraj%pk * p_t / p00**kappa

 call vtrisolvepert(conf%im,conf%jm,conf%lm,ltraj%akv,ltraj%bkv,ltraj%ckv,p_u ,2,1)
 call vtrisolvepert(conf%im,conf%jm,conf%lm,ltraj%akv,ltraj%bkv,ltraj%ckv,p_v ,2,1)
 call vtrisolvepert(conf%im,conf%jm,conf%lm,ltraj%aks,ltraj%bks,ltraj%cks,p_t ,2,1)
 call vtrisolvepert(conf%im,conf%jm,conf%lm,ltraj%akq,ltraj%bkq,ltraj%ckq,pert%tracers(:,:,:,1),2,1)

 do f = 2, size(traj%tracer_names)
   call vtrisolvepert(conf%im,conf%jm,conf%lm,ltraj%akq,ltraj%bkq,ltraj%ckq,pert%tracers(:,:,:,f),2,0)
 end do

 !pt2t adjoint
 p_t = p00**kappa * p_t / ltraj%pk

 !Nullify
 nullify(p_u)
 nullify(p_v)
 nullify(p_t)
 nullify(p_delp)
 nullify(ltraj)

endsubroutine step_ad

! ------------------------------------------------------------------------------

subroutine delete(self,conf)

 implicit none
 class(fv3jedi_lm_turbulence_type), intent(inout) :: self
 type(fv3jedi_lm_conf),             intent(in)    :: conf

 integer :: n

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
 type(local_cnst_turbulence),   intent(in)    :: lcnst
 type(fv3jedi_lm_traj), target, intent(in)    :: traj
 type(local_traj_turbulence),   intent(inout) :: ltraj

 integer :: i,j,l,im,jm,lm,isc,iec,jsc,jec
 integer :: it_qv,it_qi,it_ql

 real(kind_real), allocatable, dimension(:,:,:) :: PTT1
 real(kind_real), allocatable, dimension(:,:)   :: ZPBL1, CT1
 real(kind_real), allocatable, dimension(:,:,:) :: EKV, FKV
 real(kind_real), allocatable, dimension(:,:,:) :: pet, pmt
 real(kind_real), allocatable, dimension(:,:,:) :: QIT1, QLT1
 real(kind_real), allocatable, dimension(:,:,:) :: fQi

 real(kind_real), pointer, dimension(:,:,:) :: p_u
 real(kind_real), pointer, dimension(:,:,:) :: p_v
 real(kind_real), pointer, dimension(:,:,:) :: p_t
 real(kind_real), pointer, dimension(:,:,:) :: p_delp
 real(kind_real), pointer, dimension(:,:,:) :: p_qv
 real(kind_real), pointer, dimension(:,:)   :: p_frland
 real(kind_real), pointer, dimension(:,:)   :: p_frocean
 real(kind_real), pointer, dimension(:,:)   :: p_varflt
 real(kind_real), pointer, dimension(:,:)   :: p_cm
 real(kind_real), pointer, dimension(:,:)   :: p_cq
 real(kind_real), pointer, dimension(:,:)   :: p_ustar
 real(kind_real), pointer, dimension(:,:)   :: p_bstar

 im = conf%im
 jm = conf%jm
 lm = conf%lm

 isc = conf%isc
 iec = conf%iec
 jsc = conf%jsc
 jec = conf%jec

 allocate(PTT1 (1:im,1:jm,1:lm))
 allocate(EKV  (1:im,1:jm,1:lm))
 allocate(FKV  (1:im,1:jm,1:lm))
 allocate(pet  (1:im,1:jm,0:lm))
 allocate(pmt  (1:im,1:jm,1:lm))
 allocate(QIT1 (1:im,1:jm,1:lm))
 allocate(QLT1 (1:im,1:jm,1:lm))
 allocate(fQi  (1:im,1:jm,1:lm))
 allocate(ZPBL1(1:im,1:jm))
 allocate(CT1  (1:im,1:jm))

 p_u      (1:,1:,1:) => traj%u
 p_v      (1:,1:,1:) => traj%v
 p_t      (1:,1:,1:) => traj%t
 p_delp   (1:,1:,1:) => traj%delp

 p_frland (1:,1:) => traj%FRLAND
 p_frocean(1:,1:) => traj%FROCEAN
 p_varflt (1:,1:) => traj%VARFLT
 p_cm     (1:,1:) => traj%CM
 p_cq     (1:,1:) => traj%CQ
 p_ustar  (1:,1:) => traj%USTAR
 p_bstar  (1:,1:) => traj%BSTAR

 call get_tracer_index(traj, 'water_vapor_mixing_ratio_wrt_moist_air', it_qv)
 p_qv     (1:,1:,1:) => traj%tracers(:,:,:,it_qv)

 !Compute pressures from delp
 call compute_pressures(im,jm,lm,conf%ptop,p_delp,pet,pmt,ltraj%pk)

 !Use local copied to avoid overwrite
 ZPBL1(1:im,1:jm) = traj%ZPBL(isc:iec,jsc:jec)
 CT1  (1:im,1:jm) = traj%CT  (isc:iec,jsc:jec)

 ptt1 = p00**kappa * p_t / ltraj%pk

 !Calculate total cloud ice and liquid trajectory
 if (conf%do_phy_mst == 0) then
    call get_tracer_index(traj, 'cloud_liquid_ice', it_qi )
    call get_tracer_index(traj, 'cloud_liquid_water', it_ql )
    QIT1(1:im,1:jm,:) = traj%tracers(isc:iec,jsc:jec,:,it_qi)
    QLT1(1:im,1:jm,:) = traj%tracers(isc:iec,jsc:jec,:,it_ql)
 else
   do l = 1,lm
     do j = 1,jm
       do i = 1,im
         call IceFraction( p_t(i,j,l), fQi(i,j,l) )
       enddo
     enddo
   enddo

   QIT1(1:im,1:jm,:) = (traj%QLS(isc:iec,jsc:jec,:) + traj%QCN(isc:iec,jsc:jec,:)) * fQi(1:im,1:jm,:)
   QLT1(1:im,1:jm,:) = (traj%QLS(isc:iec,jsc:jec,:) + traj%QCN(isc:iec,jsc:jec,:)) * (1-fQi(1:im,1:jm,:))
 endif

 !Initialize tri-diagonal arrays
 ltraj%AKV = 0.0
 ltraj%BKV = 0.0
 ltraj%CKV = 0.0
 ltraj%AKS = 0.0
 ltraj%BKS = 0.0
 ltraj%CKS = 0.0
 ltraj%AKQ = 0.0
 ltraj%BKQ = 0.0
 ltraj%CKQ = 0.0

 !Call boundary layer routine. These routines will return the lower (AK*), main (BK*),
 !and upper (CK*) diagonals.

 call BL_DRIVER( IM                             , &
                 JM                             , &
                 LM                             , &
                 conf%DT                        , &
                 p_u                            , &
                 p_v                            , &
                 PTT1                           , &
                 p_QV                           , &
                 PET                            , &
                 QIT1                           , &
                 QLT1                           , &
                 p_FRLAND                       , &
                 p_FROCEAN                      , &
                 p_VARFLT                       , &
                 ZPBL1                          , &
                 p_CM                           , &
                 CT1                            , &
                 p_CQ                           , &
                 lcnst%TURBPARAMS               , &
                 lcnst%TURBPARAMSI              , &
                 p_USTAR                        , &
                 p_BSTAR                        , &
                 ltraj%AKS, ltraj%BKS, ltraj%CKS, &
                 ltraj%AKQ, ltraj%BKQ, ltraj%CKQ, &
                 ltraj%AKV, ltraj%BKV, ltraj%CKV, &
                 EKV, FKV                         )

 !Solver part 1, perform LU decomposition
 call VTRILUPERT(IM,JM,LM,ltraj%AKV,ltraj%BKV,ltraj%CKV)
 call VTRILUPERT(IM,JM,LM,ltraj%AKS,ltraj%BKS,ltraj%CKS)
 call VTRILUPERT(IM,JM,LM,ltraj%AKQ,ltraj%BKQ,ltraj%CKQ)

 deallocate(PTT1 )
 deallocate(ZPBL1)
 deallocate(CT1  )
 deallocate(EKV  )
 deallocate(FKV  )
 deallocate(pet  )
 deallocate(pmt  )
 deallocate(QIT1 )
 deallocate(QLT1 )
 deallocate(fQi  )

 nullify(p_u)
 nullify(p_v)
 nullify(p_t)
 nullify(p_delp)
 nullify(p_qv)
 nullify(p_frland)
 nullify(p_frocean)
 nullify(p_varflt)
 nullify(p_cm)
 nullify(p_cq)
 nullify(p_ustar)
 nullify(p_bstar)

 if (conf%saveltraj) ltraj%set = .true.

endsubroutine set_ltraj

! ------------------------------------------------------------------------------

subroutine allocate_ltraj(im,jm,lm,ltraj)

 integer, intent(in) :: im,jm,lm
 type(local_traj_turbulence), intent(inout) :: ltraj

 allocate(ltraj%pk(im,jm,lm))
 allocate(ltraj%akv(im,jm,lm))
 allocate(ltraj%bkv(im,jm,lm))
 allocate(ltraj%ckv(im,jm,lm))
 allocate(ltraj%aks(im,jm,lm))
 allocate(ltraj%bks(im,jm,lm))
 allocate(ltraj%cks(im,jm,lm))
 allocate(ltraj%akq(im,jm,lm))
 allocate(ltraj%bkq(im,jm,lm))
 allocate(ltraj%ckq(im,jm,lm))

endsubroutine allocate_ltraj

! ------------------------------------------------------------------------------

subroutine deallocate_ltraj(ltraj)

 type(local_traj_turbulence), intent(inout) :: ltraj

 deallocate(ltraj%pk )
 deallocate(ltraj%akv)
 deallocate(ltraj%bkv)
 deallocate(ltraj%ckv)
 deallocate(ltraj%aks)
 deallocate(ltraj%bks)
 deallocate(ltraj%cks)
 deallocate(ltraj%akq)
 deallocate(ltraj%bkq)
 deallocate(ltraj%ckq)

endsubroutine deallocate_ltraj

! ------------------------------------------------------------------------------

subroutine vtrilupert(im,jm,lm,a,b,c)

 !perform lu decomposition of tridiagonal system

 implicit none

 integer,                              intent(in   ) :: im, jm, lm
 real(kind_real), dimension(im,jm,lm), intent(in   ) :: c
 real(kind_real), dimension(im,jm,lm), intent(inout) :: a, b

 integer :: l

 b(:,:,1) = 1. / b(:,:,1)
 do l = 2,lm
    a(:,:,l) = a(:,:,l) * b(:,:,l-1)
    b(:,:,l) = 1. / ( b(:,:,l) - c(:,:,l-1) * a(:,:,l) )
 end do

end subroutine vtrilupert

! ------------------------------------------------------------------------------

subroutine vtrisolvepert(im,jm,lm,a,b,c,y,phase,ygswitch)

 !solve lu decomposed tridiagonal system

 implicit none

 !arguments
 integer,                              intent(in   ) :: im, jm, lm, phase, ygswitch
 real(kind_real), dimension(im,jm,lm), intent(in   ) :: a, b, c
 real(kind_real), dimension(im,jm,lm), intent(inout) :: y

 !locals
 integer :: l

 if (phase == 1) then !TLM

    !solve (lu)ynew = yold

    !sweep down, modifying rhs with multiplier a
    do l=2,lm
       y(:,:,l) = y(:,:,l) - a(:,:,l)*y(:,:,l-1)
    enddo

    !sweep up, solving for updated value. note b has the inverse of the main diagonal
    if (ygswitch == 1) then !winds, temperature and q
       y(:,:,lm) = y(:,:,lm)*b(:,:,lm) !ygprime = 0
    else !tracers
       y(:,:,lm) = y(:,:,lm)*b(:,:,lm-1)/(b(:,:,lm-1) - a(:,:,lm)*(1.0+c(:,:,lm-1)*b(:,:,lm-1) ))
    endif

    do l=lm-1,1,-1
       y(:,:,l) = b(:,:,l)*(y(:,:,l)-c(:,:,l)*y(:,:,l+1) )
    enddo

 elseif (phase == 2) then !ADM

    if (ygswitch == 1) then !can solve (lu)'ynew = u'l'ynew = yold

       !sweep down but with u' instead of l
       y(:,:,1) = y(:,:,1)*b(:,:,1)
       do l=2,lm
          y(:,:,l) = b(:,:,l) * (y(:,:,l) - c(:,:,l-1)*y(:,:,l-1))
       enddo

       !sweep up but with l' instead of u
       do l=lm-1,1,-1
          y(:,:,l) = y(:,:,l) - a(:,:,l+1)*y(:,:,l+1)
       enddo

    else !change for surface means transpose doesnt work so use line-by-line adjoint

       !adjoint of sweep up, solving for updated value. note b has the inverse of the main diagonal
       do l=1,lm-1,1
          y(:,:,l+1) = y(:,:,l+1) - c(:,:,l)*b(:,:,l)*y(:,:,l)
          y(:,:,l) = b(:,:,l)*y(:,:,l)
       end do

       !adjoint of surface fix
       y(:,:,lm) = b(:,:,lm-1)*y(:,:,lm)/(b(:,:,lm-1)-a(:,:,lm)*(c(:,:,lm-1)*b(:,:,lm-1)+1.0))

       !adjoint of sweep down, modifying rhs with multiplier a
       do l=lm,2,-1
          y(:,:,l-1) = y(:,:,l-1) - a(:,:,l)*y(:,:,l)
       end do

    endif

 endif

end subroutine vtrisolvepert


end module fv3jedi_lm_turbulence_mod
