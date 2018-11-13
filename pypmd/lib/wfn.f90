subroutine wfn_driver(epscuttz, rmaxatom, nmo, nprim, icen, ityp, &
                      oexp, mo_coeff, mo_occ, natm, coords, charges, &
                      newnprim, maxgrp, numshells, npc, ngroup, &
                      icenat, nzexp, nuexp, rcutte) bind(c)

  use mod_prec, only: rp, ip
  use mod_param, only: init_param
  use mod_wfn, only: nmo_, nprims_, epscuttz_, rmaxatom_, ncent_, &
                     allocate_space_for_wfn, init_wfn, coords_, &
                     coeff_, oexp_, occ_, charges_, icen_, ityp_, &
                     deallocate_space_for_wfn, setup_wfn, filtergto, &
                     maxgrp_, numshells_, npc_, ngroup_, mgrp, ngtoh, &
                     icenat_, nzexp_, nuexp_, rcutte_, &
                     allocate_space_for_shells, deallocate_space_for_shells
  implicit none

  real(kind=rp), intent(in), value :: epscuttz
  real(kind=rp), intent(in), value :: rmaxatom
  integer(kind=ip), intent(in), value :: nmo
  integer(kind=ip), intent(in), value :: nprim
  integer(kind=ip), intent(in), value :: natm
  
  integer(kind=ip), intent(in), dimension(nprim) :: icen
  integer(kind=ip), intent(in), dimension(nprim) :: ityp
  real(kind=rp), intent(in), dimension(nprim) :: oexp
  real(kind=rp), intent(in), dimension(nmo,nprim) :: mo_coeff
  real(kind=rp), intent(in), dimension(nmo) :: mo_occ
  real(kind=rp), intent(in), dimension(3,natm) :: coords
  integer(kind=ip), intent(in), dimension(natm) :: charges
  
  integer(kind=ip), intent(out) :: newnprim
  integer(kind=ip), intent(out) :: maxgrp
  integer(kind=ip), intent(out) :: numshells
  integer(kind=ip), intent(out), dimension(natm) :: npc
  integer(kind=ip), intent(out), dimension(natm) :: ngroup
  integer(kind=ip), intent(out), dimension(natm,nprim) :: icenat
  integer(kind=ip), intent(out), dimension(natm,mgrp) :: nzexp
  integer(kind=ip), intent(out), dimension(natm,mgrp,ngtoh) :: nuexp
  real(kind=rp), intent(out), dimension(natm,mgrp) :: rcutte

  call init_param ()
  call init_wfn ()

  nprims_ = nprim
  nmo_ = nmo
  ncent_ = natm
  epscuttz_ = epscuttz
  rmaxatom_ = rmaxatom

  call allocate_space_for_wfn ()
  
  coords_ = coords
  coeff_ = mo_coeff
  oexp_ = oexp
  occ_ = mo_occ
  charges_ = charges
  icen_ = icen
  ityp_ = ityp

  call setup_wfn ()
  call allocate_space_for_shells ()
  call filtergto ()
  
  newnprim = nprims_
  maxgrp = maxgrp_
  numshells = numshells_
  npc = npc_
  ngroup = ngroup_
  icenat = transpose(icenat_)
  nzexp = nzexp_
  nuexp = nuexp_
  rcutte = rcutte_

  call deallocate_space_for_wfn ()
  call deallocate_space_for_shells ()

end subroutine wfn_driver
