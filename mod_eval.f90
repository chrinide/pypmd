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

subroutine eval_basin(nmo, nprims, icen, ityp, oexp, ngroup, nzexp, &
                      nuexp, rcutte, mo_coeff, mo_occ, natm, coords, &
                      npoints, points, output) bind(c)

  use mod_prec, only: rp, ip
  use mod_param, only: init_param
  use mod_mole, only: ncent_, coords_, allocate_space_for_mole, &
                                       deallocate_space_for_mole
  use mod_basis, only: nmo_, nprims_, allocate_space_for_basis, init_basis, &
                       coeff_, oexp_, occ_, icen_, ityp_, nzexp_, nuexp_, & 
                       deallocate_space_for_basis, ngroup_, mgrp, ngtoh, &
                       rcutte_
  use mod_fields, only: basin
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
  real(kind=rp), intent(out), dimension(3,npoints) :: output

  integer(kind=ip) :: i
  real(kind=rp) :: rho, kin, lap

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

  !$omp parallel do &
  !$omp default(none) &
  !$omp private(i,rho,kin,lap) & 
  !$omp shared(npoints,output,points) &
  !$omp schedule(dynamic) 
  do i = 1,npoints
    call basin(points(:,i),rho,kin,lap)
    output(1,i) = rho
    output(2,i) = kin
    output(3,i) = lap
  end do
  !$omp end parallel do

  call deallocate_space_for_basis ()
  call deallocate_space_for_mole ()

end subroutine eval_basin

subroutine eval_rho(nmo, nprims, icen, ityp, oexp, ngroup, nzexp, &
                        nuexp, rcutte, mo_coeff, mo_occ, natm, coords, &
                        npoints, points, output) bind(c)

  use mod_prec, only: rp, ip
  use mod_param, only: init_param
  use mod_mole, only: ncent_, coords_, allocate_space_for_mole, &
                                       deallocate_space_for_mole
  use mod_basis, only: nmo_, nprims_, allocate_space_for_basis, init_basis, &
                       coeff_, oexp_, occ_, icen_, ityp_, nzexp_, nuexp_, & 
                       deallocate_space_for_basis, ngroup_, mgrp, ngtoh, &
                       rcutte_
  use mod_fields, only: rho_shell=>density_shell
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
  real(kind=rp), intent(out), dimension(npoints) :: output

  integer(kind=ip) :: i
  real(kind=rp) :: rho

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

  !$omp parallel do &
  !$omp default(none) &
  !$omp private(i,rho) & 
  !$omp shared(npoints,output,points) &
  !$omp schedule(dynamic) 
  do i = 1,npoints
    call rho_shell(points(:,i),rho)
    output(i) = rho
  end do
  !$omp end parallel do

  call deallocate_space_for_basis ()
  call deallocate_space_for_mole ()

end subroutine eval_rho

subroutine eval_rho_gradmod(nmo, nprims, icen, ityp, oexp, ngroup, nzexp, &
                            nuexp, rcutte, mo_coeff, mo_occ, natm, coords, &
                            npoints, points, output) bind(c)

  use mod_prec, only: rp, ip
  use mod_param, only: init_param
  use mod_mole, only: ncent_, coords_, allocate_space_for_mole, &
                                       deallocate_space_for_mole
  use mod_basis, only: nmo_, nprims_, allocate_space_for_basis, init_basis, &
                       coeff_, oexp_, occ_, icen_, ityp_, nzexp_, nuexp_, & 
                       deallocate_space_for_basis, ngroup_, mgrp, ngtoh, &
                       rcutte_
  use mod_fields, only: rho_grad_shell=>density_grad_shell
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
  real(kind=rp), intent(out), dimension(2,npoints) :: output

  integer(kind=ip) :: i
  real(kind=rp) :: rho, grad(3), gradmod

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

  !$omp parallel do &
  !$omp default(none) &
  !$omp private(i,rho,grad,gradmod) & 
  !$omp shared(npoints,output,points) &
  !$omp schedule(dynamic) 
  do i = 1,npoints
    call rho_grad_shell(points(:,i),rho,grad,gradmod)
    output(1,i) = rho
    output(2,i) = gradmod
  end do
  !$omp end parallel do

  call deallocate_space_for_basis ()
  call deallocate_space_for_mole ()

end subroutine eval_rho_gradmod

