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
! get the number of xmat in the symmetry operator list.
! A -1 value will be returned if xmat is not in the list.
!
integer(kind=ip) function numbersym (xmat)

  use mod_prec, only: rp, ip
  use mod_sym, only: nopsym, opsym, toleqvm
  implicit none
  real(kind=rp), dimension(3,3), intent(in) :: xmat

  integer(kind=ip) :: i, j, k
  real(kind=rp) :: diff

  numbersym = -1_ip
  do i = 1,nopsym
    diff = 0.0_rp
    do j = 1,3
      do k = 1,3
        diff = diff + abs(xmat(j,k)-opsym(i,j,k))
      end do
    end do
    if (diff .le. TOLeqvm) then
      numbersym = i
      return
    end if
  end do

end function numbersym
