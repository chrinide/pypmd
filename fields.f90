subroutine driver(nmo, nprims, icen, ityp, oexp, ngroup, nzexp, &
                  nuexp, rcutte, mo_coeff, mo_occ, natm, coords) bind(c)

  use mod_prec, only: rp, ip
  use iso_fortran_env, only: uout=>output_unit
  use mod_io, only: string
  use mod_param, only: init_param
  use mod_mole, only: ncent_, coords_, allocate_space_for_mole, &
                                       deallocate_space_for_mole
  use mod_basis, only: nmo_, nprims_, allocate_space_for_basis, init_basis, &
                       coeff_, oexp_, occ_, icen_, ityp_, nzexp_, nuexp_, & 
                       deallocate_space_for_basis, ngroup_, mgrp, ngtoh, &
                       rcutte_
  use mod_fields, only: density_grad, density_grad_shell
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

  integer(kind=ip) :: i, j, k, suma
  real(kind=rp) :: xpoint(3)
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

  do i = 1,ncent_
    write (uout,'(1x,a,1x,i0,1x,3(1x,f0.8))') &
    string("# Atomic coordinates for center"), i, coords_(:,i)
  end do
  do i = 1,ncent_
    suma = 0_ip
    write (uout,'(1x,a,1x,i0,1x,i0)') string("# ngroup "), i, ngroup_(i)
    do j = 1,ngroup_(i)
      write (uout,'(1x,a,1x,i0,1x,i0)') string("# nzexp"), j, nzexp_(i,j)
      write (uout,'(1x,a,1x,i0,1x,f0.8)') string("# rcutte"), j, rcutte_(i,j)
      suma = suma + nzexp_(i,j)
      do k = 1,nzexp_(i,j)
        write (uout,'(1x,a,1x,i0)') string("# nuexp"), nuexp_(i,k,j)
      end do
    end do
    write (uout,'(1x,a,1x,i0)') string("# suma"), suma
  end do

  xpoint(1) = 0.0_rp
  xpoint(2) = 0.0_rp
  xpoint(3) = 0.3_rp
  call density_grad (xpoint,rho,grad,gradmod)
  write (uout,'(1x,a)') string('# --------------------------------------')
  write (uout,'(1x,a)') string('# --- Follow values for test fsimple ---')
  write (uout,'(1x,a,1x,3(1x,f0.8))') string('# Coordinates of point'), xpoint
  write (uout,'(1x,a,1x,f0.8)')  string('# Density value'), rho
  write (uout,'(1x,a,1x,3(1x,f0.8))') string('# Gradient value'), grad(:)
  write (uout,'(1x,a,1x,f0.8)') string('# Gradmod value'), gradmod
  write (uout,'(1x,a)') string('# --------------------------------------')

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

  call deallocate_space_for_basis ()
  call deallocate_space_for_mole ()

end subroutine driver
