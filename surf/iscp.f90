logical function iscp (p,nuc)
 
  use mod_prec, only: rp, ip
  use mod_wfn, only: ncent
  use mod_surf, only: xyzrho, epsiscp
  implicit none

  integer(kind=ip), intent(out) :: nuc
  real(kind=rp), intent(in) :: p(3)

  real(kind=rp) :: gradmod, grad(3), rho
  real(kind=rp) :: x(3)
  integer(kind=ip) :: i

  interface
    subroutine pointr1 (p,rho,grad,gradmod)
      import rp
      real(kind=rp), intent(in) :: p(3)
      real(kind=rp), intent(out) :: grad(3)
      real(kind=rp), intent(out) :: rho
      real(kind=rp), intent(out) :: gradmod
    end subroutine
  end interface

  iscp = .false.
  nuc = 0
  x(:) = p(:)
  call pointr1 (x,rho,grad,gradmod)
  do i = 1,ncent
    if (abs(p(1)-xyzrho(i,1)).lt.epsiscp .and. &
        abs(p(2)-xyzrho(i,2)).lt.epsiscp .and. &
        abs(p(3)-xyzrho(i,3)).lt.epsiscp) then
      iscp = .true.
      nuc = i
    end if
  end do

  if (gradmod.le.1d-10) then
    iscp = .true.
    if (rho.le.1d-10) nuc = -1
  end if

end function
