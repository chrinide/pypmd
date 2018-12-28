subroutine surf_driver(nmo, nprims, icen, ityp, oexp, ngroup, nzexp, &
                       nuexp, rcutte, mo_coeff, mo_occ, natm, coords, &
                       npang, inuc, xyzrho, filename, ct, st, cp, sp, &
                       angw, backend, ntrial, epsiscp, epsilon, rmaxsurf, &
                       step, mstep, nlimsurf, rlimsurf, verbose) bind(c)

  use mod_prec, only: rp, ip
  use iso_c_binding, only: c_null_char, c_char
  use iso_fortran_env, only: uout=>output_unit
  use mod_io, only: string, fourchar, mline
  use mod_param, only: init_param
  use mod_mole, only: ncent_, coords_, allocate_space_for_mole, &
                                       deallocate_space_for_mole
  use mod_basis, only: nmo_, nprims_, allocate_space_for_basis, init_basis, &
                       coeff_, oexp_, occ_, icen_, ityp_, nzexp_, nuexp_, & 
                       deallocate_space_for_basis, ngroup_, mgrp, ngtoh, &
                       rcutte_
  use mod_surface, only: init_surf, rlimsurf_, nlimsurf_, minter, epsilon_, &
                      xnuc_, inuc_, xyzrho_, npang_, rmaxsurf_, surf, &
                      allocate_space_for_surface, deallocate_space_for_surface, &
                      epsiscp_, mstep_, ntrial_, steeper_, step_, rprimer_
  use mod_fields, only: density_grad_shell
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
  
  integer(kind=ip), intent(out), dimension(natm) :: ngroup
  integer(kind=ip), intent(out), dimension(natm,mgrp) :: nzexp
  integer(kind=ip), intent(out), dimension(natm,ngtoh,mgrp) :: nuexp
  real(kind=rp), intent(out), dimension(natm,mgrp) :: rcutte

  integer(kind=ip), intent(in), value :: npang, inuc, verbose
  real(kind=rp), intent(in), dimension(3,natm) :: xyzrho
  character(kind=c_char,len=1), intent(in) :: filename(*)

  integer(kind=ip), intent(in), value :: backend, ntrial, mstep
  real(kind=rp), intent(in), value :: epsiscp, epsilon , rmaxsurf, step
  real(kind=rp), intent(in), dimension(npang) :: ct, st, cp, sp, angw
  real(kind=rp), intent(out) :: rlimsurf(ntrial,npang) 
  integer(kind=ip), intent(out) :: nlimsurf(npang) 

  real(kind=rp) :: rmaxs, rmins, rsurf(minter,2)
  integer(kind=ip) :: i, lsu, ltxt, j, nchars, nsurf, k
  character(len=mline) :: files
  character(len=:), allocatable :: filename_
  character(len=4) :: d4
  integer(kind=ip) :: suma
  real(kind=rp) :: xpoint(3)
  real(kind=rp) :: rho, grad(3), gradmod

  i = 1
  do
    if (filename(i) == c_null_char) exit
    i = i + 1
  end do
  nchars = i - 1
  allocate(character(len=nchars) :: filename_)
  filename_ = transfer(filename(1:nchars), filename_)

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
  rmaxsurf_ = rmaxsurf
  step_ = step
  mstep_ = mstep

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

  ! Init
  ltxt = 999
  lsu = 998
  d4 = fourchar(inuc_)
  files = trim(filename_)//".surf"

  if (verbose == 1) then
    write (*,*) "nmp",nmo_
    write (*,*) "nprims",nprims_
    write (*,*) "ncent",ncent_
    write (*,*) "npang",npang_
    write (*,*) "inuc",inuc_
    write (*,*) "ntrial",ntrial_
    write (*,*) "epsiscp",epsiscp_
    write (*,*) "epsilon",epsilon_
    write (*,*) "rmaxsurf", rmaxsurf_
    write (*,*) "step", step_
    write (*,*) "mstep",mstep_
    write (*,*) "xnuc",xnuc_
    write (*,*) "rprimer",rprimer_
    do i = 1,ncent_
      write (*,*) xyzrho_(:,i)
    end do
  end if
  xpoint(1) = 0.0_rp
  xpoint(2) = 0.0_rp
  xpoint(3) = 0.3_rp
  call density_grad_shell (xpoint,rho,grad,gradmod)
  write (uout,'(1x,a)') string('# --------------------------------------')
  write (uout,'(1x,a)') string('# --- Follow values for test fshells ---')
  write (uout,'(1x,a,1x,3(1x,f0.8))') string('# Coordinates of point'), xpoint
  write (uout,'(1x,a,1x,f0.8)')  string('# Density value'), rho
  write (uout,'(1x,a,1x,3(1x,f0.8))') string('# Gradient value'), grad(:)
  write (uout,'(1x,a,1x,f0.8)') string('# Gradmod value'), gradmod
  write (uout,'(1x,a)') string('# --------------------------------------')

  !$omp parallel default(none) &
  !$omp private(j,nsurf,rsurf) &
  !$omp shared(npang_,ct,st,cp,sp,rlimsurf,nlimsurf)
  !$omp do schedule(dynamic)
  do j = 1,npang_
    !write (*,*) ct(j), st(j), cp(j), sp(j), angw(j)
    call surf (ct(j),st(j),cp(j),sp(j),rsurf,nsurf)
    do k = 1,nsurf
      rlimsurf(k,j) = rsurf(k,2)
    end do
    nlimsurf(j) = nsurf
  end do
  !$omp end do nowait
  !$omp end parallel

  open (ltxt,file=trim(files)//"-txt"//d4)
  open (lsu,file=trim(files)//d4,form='unformatted')
  write (ltxt,1111) npang_, inuc_
  write (lsu) npang_, inuc_
  write (ltxt,3333) (nlimsurf(j),j=1,npang)
  write (ltxt,1090)
  write (lsu) (nlimsurf(j),j=1,npang)
  rmins = 1000_rp
  rmaxs = 0.0_rp
  do j = 1,npang_
    nsurf = nlimsurf(j)
    rmins = min(rmins,rlimsurf(1,j))
    rmaxs = max(rmaxs,rlimsurf(nsurf,j))
    write (ltxt,2222) ct(j),st(j),cp(j),sp(j),angw(j),(rlimsurf(k,j),k=1,nsurf)
    write (lsu) ct(j),st(j),cp(j),sp(j),angw(j),(rlimsurf(k,j),k=1,nsurf)
  end do
  close (ltxt)
  close (lsu)

1090 format (9x,'cos(theta)',13x,'sin(theta)',13x,'cos(phi)',15x,'sin(phi)',15x,'weight')
1111 format (2(1x,i5),' <--- (Angular points & Atom)')
3333 format (20(1x,i2),4x,'(Surface intersections)')
2222 format (15(1x,F22.15))

  call deallocate_space_for_basis ()
  call deallocate_space_for_mole ()
  call deallocate_space_for_surface ()

end subroutine surf_driver
