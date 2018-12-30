subroutine becke_driver(nmo, nprims, icen, ityp, oexp, ngroup, nzexp, &
                        nuexp, rcutte, mo_coeff, mo_occ, natm, coords, &
                        npoints, points, weights) bind(c)

  use mod_prec, only: rp, ip
  use iso_fortran_env, only: uout=>output_unit
  use mod_param, only: init_param
  use mod_mole, only: ncent_, coords_, allocate_space_for_mole, &
                                       deallocate_space_for_mole
  use mod_basis, only: nmo_, nprims_, allocate_space_for_basis, init_basis, &
                       coeff_, oexp_, occ_, icen_, ityp_, nzexp_, nuexp_, & 
                       deallocate_space_for_basis, ngroup_, mgrp, ngtoh, &
                       rcutte_
  use mod_fields, only: rho_grad=>density_grad_shell
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
  real(kind=rp), intent(in), dimension(natm,mgrp) :: rcutte

  integer(kind=ip), intent(in), value :: npoints
  real(kind=rp), intent(in), dimension(3,npoints) :: points
  real(kind=rp), intent(in), dimension(npoints) :: weights

  real(kind=rp) :: rho, trho, grad(3), gradmod
  integer(kind=ip) :: i

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
  rcutte_ = rcutte
  nuexp_ = nuexp

  rho = 0.0_rp
  do i = 1,npoints
    call rho_grad (points(:,i),trho,grad,gradmod)
    rho = rho + trho*weights(i)
  end do
  write (*,*) "Rho",rho

  call deallocate_space_for_basis ()
  call deallocate_space_for_mole ()

end subroutine becke_driver
