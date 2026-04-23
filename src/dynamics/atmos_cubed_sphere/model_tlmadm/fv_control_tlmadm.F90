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
!
!----------------
! FV contro panel
!----------------

!Prepare the FV_AtmP derived type that holds the perturbation variables and the
!coefficients used for advection, remapping and damping in the tangent linear
!and adjoint versions of FV3.

module fv_control_tlmadm_mod

 use fv_arrays_nlm_mod,     only: fv_atmos_type
 use fv_arrays_tlmadm_mod, only: fv_atmos_pert_type, allocate_fv_atmos_pert_type, deallocate_fv_atmos_pert_type
 use fms_mod,           only: check_nml_error
 use mpp_mod,           only: stdlog, mpp_pe, mpp_root_pe

 implicit none
 private

 integer, public :: ngrids = 1

!Convenience pointer to AtmP(n)
 logical, pointer :: split_hord
 integer, pointer :: hord_mt_pert
 integer, pointer :: hord_vt_pert
 integer, pointer :: hord_tm_pert
 integer, pointer :: hord_dp_pert
 integer, pointer :: hord_tr_pert
 logical, pointer :: split_kord
 integer, pointer :: kord_mt_pert
 integer, pointer :: kord_wz_pert
 integer, pointer :: kord_tm_pert
 integer, pointer :: kord_tr_pert
 logical, pointer :: split_damp
 logical, pointer :: do_vort_damp_pert
 integer, pointer :: nord_pert
 real,    pointer :: dddmp_pert
 real,    pointer :: d2_bg_pert
 real,    pointer :: d4_bg_pert
 real,    pointer :: vtdm4_pert
 real,    pointer :: d2_bg_k1_pert
 real,    pointer :: d2_bg_k2_pert
 real,    pointer :: d2_bg_ks_pert
 logical, pointer :: split_damp_tr
 integer, pointer :: nord_tr_pert
 real,    pointer :: trdm2_pert
 integer, pointer :: n_sponge_pert
 logical, pointer :: hord_ks_pert
 integer, pointer :: hord_mt_ks_pert
 integer, pointer :: hord_vt_ks_pert
 integer, pointer :: hord_tm_ks_pert
 integer, pointer :: hord_dp_ks_pert
 integer, pointer :: hord_tr_ks_pert
 logical, pointer :: hord_ks_traj
 integer, pointer :: hord_mt_ks_traj
 integer, pointer :: hord_vt_ks_traj
 integer, pointer :: hord_tm_ks_traj
 integer, pointer :: hord_dp_ks_traj
 integer, pointer :: hord_tr_ks_traj


 public fv_init_pert, fv_end_pert

 contains

!-------------------------------------------------------------------------------

 subroutine fv_init_pert(Atm, AtmP, inputpert_filename)

 type(fv_atmos_type), allocatable, intent(inout), target :: Atm(:)
 type(fv_atmos_pert_type), allocatable, intent(inout), target :: AtmP(:)
 character(len=*), intent(in) :: inputpert_filename

 integer :: n, ntilesMe

  allocate(AtmP(ngrids))

  call run_setup_pert(AtmP, Atm, inputpert_filename)

  ntilesMe = size(AtmP(:))

  do n=1,ntilesMe

     call allocate_fv_atmos_pert_type( AtmP(n), Atm(n)%bd%isd, Atm(n)%bd%ied, Atm(n)%bd%jsd, Atm(n)%bd%jed, &
                                                Atm(n)%bd%isc, Atm(n)%bd%iec, Atm(n)%bd%jsc, Atm(n)%bd%jec, &
                                                Atm(n)%flagstruct%npz,  Atm(n)%flagstruct%ncnst )

  enddo

 end subroutine fv_init_pert

!-------------------------------------------------------------------------------

 subroutine fv_end_pert(AtmP)

  type(fv_atmos_pert_type), intent(inout) :: AtmP(:)

  integer :: n, ntilesMe

  ntilesMe = size(AtmP(:))

  do n=1,ntilesMe

   call deallocate_fv_atmos_pert_type(AtmP(n))

  enddo

 end subroutine fv_end_pert

!-------------------------------------------------------------------------------

 subroutine setup_pointers_pert(AtmP)

  type(fv_atmos_pert_type), intent(INOUT), target :: AtmP

   !Linearized model pointers
   split_hord        => AtmP%flagstruct%split_hord
   hord_mt_pert      => AtmP%flagstruct%hord_mt_pert
   hord_vt_pert      => AtmP%flagstruct%hord_vt_pert
   hord_tm_pert      => AtmP%flagstruct%hord_tm_pert
   hord_dp_pert      => AtmP%flagstruct%hord_dp_pert
   hord_tr_pert      => AtmP%flagstruct%hord_tr_pert
   split_kord        => AtmP%flagstruct%split_kord
   kord_mt_pert      => AtmP%flagstruct%kord_mt_pert
   kord_wz_pert      => AtmP%flagstruct%kord_wz_pert
   kord_tm_pert      => AtmP%flagstruct%kord_tm_pert
   kord_tr_pert      => AtmP%flagstruct%kord_tr_pert
   split_damp        => AtmP%flagstruct%split_damp
   nord_pert         => AtmP%flagstruct%nord_pert
   dddmp_pert        => AtmP%flagstruct%dddmp_pert
   d2_bg_pert        => AtmP%flagstruct%d2_bg_pert
   d4_bg_pert        => AtmP%flagstruct%d4_bg_pert
   do_vort_damp_pert => AtmP%flagstruct%do_vort_damp_pert
   d2_bg_k1_pert     => AtmP%flagstruct%d2_bg_k1_pert
   d2_bg_k2_pert     => AtmP%flagstruct%d2_bg_k2_pert
   d2_bg_ks_pert     => AtmP%flagstruct%d2_bg_ks_pert
   vtdm4_pert        => AtmP%flagstruct%vtdm4_pert
   split_damp_tr     => AtmP%flagstruct%split_damp_tr
   nord_tr_pert      => AtmP%flagstruct%nord_tr_pert
   trdm2_pert        => AtmP%flagstruct%trdm2_pert
   n_sponge_pert     => AtmP%flagstruct%n_sponge_pert
   hord_ks_pert      => AtmP%flagstruct%hord_ks_pert
   hord_mt_ks_pert   => AtmP%flagstruct%hord_mt_ks_pert
   hord_vt_ks_pert   => AtmP%flagstruct%hord_vt_ks_pert
   hord_tm_ks_pert   => AtmP%flagstruct%hord_tm_ks_pert
   hord_dp_ks_pert   => AtmP%flagstruct%hord_dp_ks_pert
   hord_tr_ks_pert   => AtmP%flagstruct%hord_tr_ks_pert
   hord_ks_traj      => AtmP%flagstruct%hord_ks_traj
   hord_mt_ks_traj   => AtmP%flagstruct%hord_mt_ks_traj
   hord_vt_ks_traj   => AtmP%flagstruct%hord_vt_ks_traj
   hord_tm_ks_traj   => AtmP%flagstruct%hord_tm_ks_traj
   hord_dp_ks_traj   => AtmP%flagstruct%hord_dp_ks_traj
   hord_tr_ks_traj   => AtmP%flagstruct%hord_tr_ks_traj

  end subroutine setup_pointers_pert

!-------------------------------------------------------------------------------

 subroutine run_setup_pert(AtmP, Atm, inputpert_filename)

  type(fv_atmos_pert_type), intent(inout), target :: AtmP(:)
  type(fv_atmos_type), intent(inout), target :: Atm(:)
  character(len=*), intent(in) :: inputpert_filename

  integer :: f_unit, n, ierr, ios, unit
  logical :: file_exists

  namelist /fv_core_pert_nml/split_hord, hord_mt_pert, hord_vt_pert, hord_tm_pert, hord_dp_pert, hord_tr_pert, &
                             split_kord, kord_mt_pert, kord_wz_pert, kord_tm_pert, kord_tr_pert, split_damp, do_vort_damp_pert, &
                             nord_pert, dddmp_pert, d2_bg_pert, d4_bg_pert, vtdm4_pert, d2_bg_k1_pert, d2_bg_k2_pert, d2_bg_ks_pert, &
                             split_damp_tr, nord_tr_pert, trdm2_pert, &
                             n_sponge_pert, &
                             hord_ks_traj, hord_mt_ks_traj, hord_vt_ks_traj, hord_tm_ks_traj, hord_dp_ks_traj, hord_tr_ks_traj, &
                             hord_ks_pert, hord_mt_ks_pert, hord_vt_ks_pert, hord_tm_ks_pert, hord_dp_ks_pert, hord_tr_ks_pert

  do n=1,size(AtmP)

     call setup_pointers_pert(AtmP(n))

     ! Check file exists
     inquire( file=trim(inputpert_filename), exist = file_exists )

     if (file_exists) then

       ! Open the file
       open(newunit=f_unit, file=trim(inputpert_filename), status='old', form='formatted', action='read')

       !Read linearized FVCORE namelist
       read (f_unit,fv_core_pert_nml,iostat=ios)
       ierr = check_nml_error(ios,'fv_core_pert_nml')

       close(f_unit)

       unit = stdlog()
       write(unit, nml=fv_core_pert_nml)

     endif

     !Unless specfied the trajectory uses the coeffs suitable for the perts
     if (.not. AtmP(n)%flagstruct%split_damp) then
        Atm(n)%flagstruct%nord         = AtmP(n)%flagstruct%nord_pert
        Atm(n)%flagstruct%dddmp        = AtmP(n)%flagstruct%dddmp_pert
        Atm(n)%flagstruct%d2_bg        = AtmP(n)%flagstruct%d2_bg_pert
        Atm(n)%flagstruct%d4_bg        = AtmP(n)%flagstruct%d4_bg_pert
        Atm(n)%flagstruct%do_vort_damp = AtmP(n)%flagstruct%do_vort_damp_pert
        Atm(n)%flagstruct%vtdm4        = AtmP(n)%flagstruct%vtdm4_pert
        Atm(n)%flagstruct%d2_bg_k1     = AtmP(n)%flagstruct%d2_bg_k1_pert
        Atm(n)%flagstruct%d2_bg_k2     = AtmP(n)%flagstruct%d2_bg_k2_pert
     endif

     !Unless specfied the trajectory uses the coeffs suitable for the perts
     if (.not. AtmP(n)%flagstruct%split_damp_tr) then
        Atm(n)%flagstruct%nord_tr = AtmP(n)%flagstruct%nord_tr_pert
        Atm(n)%flagstruct%trdm2   = AtmP(n)%flagstruct%trdm2_pert
     endif

     !Unless specfied the trajectory uses hord suitable for the perts
     if (.not. AtmP(n)%flagstruct%split_hord) then
        Atm(n)%flagstruct%hord_mt = AtmP(n)%flagstruct%hord_mt_pert
        Atm(n)%flagstruct%hord_vt = AtmP(n)%flagstruct%hord_vt_pert
        Atm(n)%flagstruct%hord_tm = AtmP(n)%flagstruct%hord_tm_pert
        Atm(n)%flagstruct%hord_dp = AtmP(n)%flagstruct%hord_dp_pert
        Atm(n)%flagstruct%hord_tr = AtmP(n)%flagstruct%hord_tr_pert
     endif

     !Unless specfied the trajectory uses hord suitable for the perts
     if (.not. AtmP(n)%flagstruct%split_kord) then
        Atm(n)%flagstruct%kord_mt = AtmP(n)%flagstruct%kord_mt_pert
        Atm(n)%flagstruct%kord_wz = AtmP(n)%flagstruct%kord_wz_pert
        Atm(n)%flagstruct%kord_tm = AtmP(n)%flagstruct%kord_tm_pert
        Atm(n)%flagstruct%kord_tr = AtmP(n)%flagstruct%kord_tr_pert
     endif

     if (mpp_pe() == mpp_root_pe()) then

        print*, ''
        print*, '|-----------------------------------------------|'
        print*, '| Advection, remapping and damping coefficients |'
        print*, '|-----------------------------------------------|'
        print*, ''
        print*, ' Splitting (off for speed, on for accuracy)'
        print*, '  split_hord = ', split_hord
        print*, '  split_kord = ', split_kord
        print*, '  split_damp = ', split_damp
        print*, '  split_damp_tr = ', split_damp_tr
        print*, ''
        print*, ' Advection of the trajectory'
        print*, '  hord_mt = ', Atm(n)%flagstruct%hord_mt
        print*, '  hord_vt = ', Atm(n)%flagstruct%hord_vt
        print*, '  hord_tm = ', Atm(n)%flagstruct%hord_tm
        print*, '  hord_dp = ', Atm(n)%flagstruct%hord_dp
        print*, '  hord_tr = ', Atm(n)%flagstruct%hord_tr
        print*, ''
        print*, ' Advection of the perturbations'
        print*, '  hord_mt_pert = ', AtmP(n)%flagstruct%hord_mt_pert
        print*, '  hord_vt_pert = ', AtmP(n)%flagstruct%hord_vt_pert
        print*, '  hord_tm_pert = ', AtmP(n)%flagstruct%hord_tm_pert
        print*, '  hord_dp_pert = ', AtmP(n)%flagstruct%hord_dp_pert
        print*, '  hord_tr_pert = ', AtmP(n)%flagstruct%hord_tr_pert
        print*, ''
        print*, ' Number of sponge layers for the perturbations'
        print*, '  n_sponge_pert = ', AtmP(n)%flagstruct%n_sponge_pert
        print*, ''
        print*, ' Sponge layer advection of the trajecotry'
        print*, '  hord_ks_traj = '   , AtmP(n)%flagstruct%hord_ks_traj
        print*, '  hord_mt_ks_traj = ', AtmP(n)%flagstruct%hord_mt_ks_traj
        print*, '  hord_vt_ks_traj = ', AtmP(n)%flagstruct%hord_vt_ks_traj
        print*, '  hord_tm_ks_traj = ', AtmP(n)%flagstruct%hord_tm_ks_traj
        print*, '  hord_dp_ks_traj = ', AtmP(n)%flagstruct%hord_dp_ks_traj
        print*, '  hord_tr_ks_traj = ', AtmP(n)%flagstruct%hord_tr_ks_traj
        print*, ' '
        print*, ' Sponge layer advection of the perturbations'
        print*, '  hord_ks_pert = '   , AtmP(n)%flagstruct%hord_ks_pert
        print*, '  hord_mt_ks_pert = ', AtmP(n)%flagstruct%hord_mt_ks_pert
        print*, '  hord_vt_ks_pert = ', AtmP(n)%flagstruct%hord_vt_ks_pert
        print*, '  hord_tm_ks_pert = ', AtmP(n)%flagstruct%hord_tm_ks_pert
        print*, '  hord_dp_ks_pert = ', AtmP(n)%flagstruct%hord_dp_ks_pert
        print*, '  hord_tr_ks_pert = ', AtmP(n)%flagstruct%hord_tr_ks_pert
        print*, ''
        print*, ' Remapping of the trajectory'
        print*, '  kord_mt = ', Atm(n)%flagstruct%kord_mt
        print*, '  kord_wz = ', Atm(n)%flagstruct%kord_wz
        print*, '  kord_tm = ', Atm(n)%flagstruct%kord_tm
        print*, '  kord_tr = ', Atm(n)%flagstruct%kord_tr
        print*, ''
        print*, ' Remapping of the perturbations'
        print*, '  kord_mt_pert = ', AtmP(n)%flagstruct%kord_mt_pert
        print*, '  kord_wz_pert = ', AtmP(n)%flagstruct%kord_wz_pert
        print*, '  kord_tm_pert = ', AtmP(n)%flagstruct%kord_tm_pert
        print*, '  kord_tr_pert = ', AtmP(n)%flagstruct%kord_tr_pert
        print*, ''
        print*, ' Dynamics damping, trajectory'
        print*, '  nord         = ', Atm(n)%flagstruct%nord
        print*, '  dddmp        = ', Atm(n)%flagstruct%dddmp
        print*, '  d2_bg        = ', Atm(n)%flagstruct%d2_bg
        print*, '  d4_bg        = ', Atm(n)%flagstruct%d4_bg
        print*, '  do_vort_damp = ', Atm(n)%flagstruct%do_vort_damp
        print*, '  vtdm4        = ', Atm(n)%flagstruct%vtdm4
        print*, '  d2_bg_k1     = ', Atm(n)%flagstruct%d2_bg_k1
        print*, '  d2_bg_k2     = ', Atm(n)%flagstruct%d2_bg_k2

        print*, ''
        print*, ' Dynamics damping, perturbations'
        print*, '  nord_pert         = ', AtmP(n)%flagstruct%nord_pert
        print*, '  dddmp_pert        = ', AtmP(n)%flagstruct%dddmp_pert
        print*, '  d2_bg_pert        = ', AtmP(n)%flagstruct%d2_bg_pert
        print*, '  d4_bg_pert        = ', AtmP(n)%flagstruct%d4_bg_pert
        print*, '  do_vort_damp_pert = ', AtmP(n)%flagstruct%do_vort_damp_pert
        print*, '  vtdm4_pert        = ', AtmP(n)%flagstruct%vtdm4_pert
        print*, '  d2_bg_k1_pert     = ', AtmP(n)%flagstruct%d2_bg_k1_pert
        print*, '  d2_bg_k2_pert     = ', AtmP(n)%flagstruct%d2_bg_k2_pert
        print*, '  d2_bg_ks_pert     = ', AtmP(n)%flagstruct%d2_bg_ks_pert
        print*, ''
        print*, ' Tracer damping, trajectory'
        print*, '  nord_tr = ', Atm(n)%flagstruct%nord_tr
        print*, '  trdm2   = ', Atm(n)%flagstruct%trdm2
        print*, ''
        print*, ' Tracer damping, perturbations'
        print*, '  nord_tr_pert = ', AtmP(n)%flagstruct%nord_tr_pert
        print*, '  trdm2_pert   = ', AtmP(n)%flagstruct%trdm2_pert
        print*, ''
        print*, '|-----------------------------------------------|'
        print*, ''

     endif

   enddo

  end subroutine run_setup_pert

end module fv_control_tlmadm_mod
