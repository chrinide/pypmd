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
! Abscissas and weights for the theta quadrature. 
!
subroutine weightheta (iq,x,w,n)

  use mod_prec, only: rp, ip
  use mod_io, only: ferror, faterr
  implicit none

  ! Arguments
  integer(kind=ip), intent(in) :: iq, n
  real(kind=rp), intent(out) :: x(n), w(n)

  if (iq.eq.1) then ! Gauss-Legendre quadrature 
    call gauleg (-1.0_rp,+1.0_rp,x,w,n)
  else if (iq.eq.2) then ! Clenshaw-Curtis quadrature
    call genclcu (-1.0_rp,+1.0_rp,x,w,n)
  else if (iq.eq.3) then ! Gauss-Chebychev of first kind 
    call gaucheb1 (-1.0_rp,+1.0_rp,x,w,n)
  else if (iq.eq.4) then ! Gauss-Chebychev of second kind 
    call gaucheb2 (x,w,n)
  else if (iq.eq.5) then ! Perez-Jorda (Gauss-Chebychev, 2nd kind) quadrature. 
    call pjt (x,w,n)
  else
    call ferror ('weightheta', 'bad theta quadrature', faterr)
  end if
 
end subroutine
