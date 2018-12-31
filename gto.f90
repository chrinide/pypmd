subroutine overlap(nmo, nprims, icen, ityp, oexp, ngroup, nzexp, &
                   nuexp, mo_coeff, mo_occ, natm, coords, output) bind(c)

  use mod_prec, only: rp, ip
  use mod_param, only: init_param
  use mod_mole, only: ncent_, coords_, allocate_space_for_mole, &
                                       deallocate_space_for_mole
  use mod_basis, only: nmo_, nprims_, allocate_space_for_basis, init_basis, &
                       coeff_, oexp_, occ_, icen_, ityp_, nzexp_, nuexp_, & 
                       deallocate_space_for_basis, ngroup_, mgrp, ngtoh
  use mod_gto, only: gtogto
  implicit none

  integer(kind=ip), intent(in), value :: natm
  real(kind=rp), intent(in), dimension(3,natm) :: coords
  
  integer(kind=ip), intent(in), value :: nmo
  integer(kind=ip), intent(in), value :: nprims
  integer(kind=ip), intent(in), dimension(nprims) :: icen
  integer(kind=ip), intent(in), dimension(nprims) :: ityp
  real(kind=rp), intent(in), dimension(nprims) :: oexp
  real(kind=rp), intent(in), dimension(nmo,nprims) :: mo_coeff
  real(kind=rp), intent(in), dimension(nmo) :: mo_occ
  
  integer(kind=ip), intent(in), dimension(natm) :: ngroup
  integer(kind=ip), intent(in), dimension(natm,mgrp) :: nzexp
  integer(kind=ip), intent(in), dimension(natm,ngtoh,mgrp) :: nuexp

  real(kind=rp), intent(out), dimension(nprims,nprims) :: output

  call init_param ()
  call init_basis ()

  nprims_ = nprims
  nmo_ = nmo
  ncent_ = natm

  call allocate_space_for_mole ()
  call allocate_space_for_basis ()
  
  coords_ = coords

  coeff_ = mo_coeff
  oexp_ = oexp
  occ_ = mo_occ
  icen_ = icen
  ityp_ = ityp
  
  ngroup_ = ngroup
  nzexp_ = nzexp
  nuexp_ = nuexp

  call gtogto (output)

  call deallocate_space_for_basis ()
  call deallocate_space_for_mole ()

end subroutine overlap
