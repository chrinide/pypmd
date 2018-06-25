! Copyright (c) 2018 
! Jose Luis Casals Sainz <jluiscasalssainz@gmail.com>,
! Ángel Martín Pendás <angel@fluor.quimica.uniovi.es> and 
! Evelio Francisco Miguelez <evelio@uniovi.es>. 
!
! promolden is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or (at
! your option) any later version.
! 
! promolden is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
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
