subroutine surf_driver(nmo, nprims, icen, ityp, oexp, ngroup, nzexp, &
                       nuexp, rcutte, mo_coeff, mo_occ, natm, coords, &
                       npang, inuc, xyzrho, filename, ct, st, cp, sp, &
                       angw, backend, ntrial, epsiscp, epsroot, rmaxsurf, &
                       epsilon, rprimer, step, mstep, nlimsurf, rlimsurf) bind(c)

  use mod_prec, only: rp, ip
  use iso_c_binding, only: c_null_char, c_char
  use iso_fortran_env, only: uout=>output_unit
  use mod_io, only: string
  use mod_param, only: init_param
  use mod_mole, only: ncent_, coords_, allocate_space_for_mole, &
                                       deallocate_space_for_mole
  use mod_basis, only: nmo_, nprims_, allocate_space_for_basis, init_basis, &
                       coeff_, oexp_, occ_, icen_, ityp_, nzexp_, nuexp_, & 
                       deallocate_space_for_basis, ngroup_, mgrp, ngtoh, &
                       rcutte_
  use mod_surface, only: init_surf, rlimsurf_, nlimsurf_, minter, epsilon_, &
                      xnuc_, inuc_, xyzrho_, npang_, rmaxsurf_, surf, epsroot_, &
                      allocate_space_for_surface, deallocate_space_for_surface, &
                      epsiscp_, mstep_, ntrial_, steeper_, step_, rprimer_
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

  integer(kind=ip), intent(in), value :: npang, inuc
  real(kind=rp), intent(in), dimension(3,natm) :: xyzrho
  character(kind=c_char,len=1), intent(in) :: filename(*)

  integer(kind=ip), intent(in), value :: backend, ntrial, mstep
  real(kind=rp), intent(in), value :: epsiscp, epsilon, rmaxsurf
  real(kind=rp), intent(in), value :: step, epsroot, rprimer
  real(kind=rp), intent(in), dimension(npang) :: ct, st, cp, sp, angw
  real(kind=rp), intent(out) :: rlimsurf(ntrial,npang) 
  integer(kind=ip), intent(out) :: nlimsurf(npang) 

  real(kind=rp) :: rsurf(minter,2)
  integer(kind=ip) :: j, k, nsurf

  call init_param ()
  call init_basis ()
  call init_surf ()

  nprims_ = nprims
  nmo_ = nmo
  ncent_ = natm
  npang_ = npang
  inuc_ = inuc
  xnuc_(:) = xyzrho(:,inuc_)
  steeper_ = backend
  ntrial_ = ntrial
  epsiscp_ = epsiscp
  epsilon_ = epsilon
  epsroot_ = epsroot ! Change in mod_surf
  rmaxsurf_ = rmaxsurf
  step_ = step
  mstep_ = mstep
  rprimer_ = rprimer

  call allocate_space_for_mole ()
  call allocate_space_for_basis ()
  call allocate_space_for_surface (ncent_,npang_,minter)
  
  coords_ = coords
  xyzrho_ = xyzrho

  coeff_ = mo_coeff
  oexp_ = oexp
  occ_ = mo_occ
  icen_ = icen
  ityp_ = ityp
  
  ngroup_ = ngroup
  nzexp_ = nzexp
  rcutte_ = rcutte
  nuexp_ = nuexp

  !$omp parallel default(none) &
  !$omp private(j,nsurf,rsurf) &
  !$omp shared(npang_,ct,st,cp,sp,rlimsurf,nlimsurf)
  !$omp do schedule(dynamic)
  do j = 1,npang_
    call surf (ct(j),st(j),cp(j),sp(j),rsurf,nsurf)
    do k = 1,nsurf
      rlimsurf(k,j) = rsurf(k,2)
    end do
    nlimsurf(j) = nsurf
  end do
  !$omp end do nowait
  !$omp end parallel
  
  call deallocate_space_for_basis ()
  call deallocate_space_for_mole ()
  call deallocate_space_for_surface ()

end subroutine surf_driver
