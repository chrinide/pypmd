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
! apply the molecular symmetry operations on the
! (xa1,ya1,za1) point and print all the images.
!
subroutine replicatesym (xa1, ya1, za1)

  use iso_fortran_env, only: uout=>output_unit
  use mod_prec, only: rp, ip
  use mod_param, only: verbose
  use mod_sym, only: nopsym, toldist, nopsym, opsym, &
                     opproper, mopsym
  implicit none
  real(kind=rp), intent(in) :: xa1, ya1, za1

  real(kind=rp) :: xx, yy, zz, xcpy(MOPSYM,3)
  integer(kind=ip) :: i, j, ncpy
  logical :: newcpy

  if (verbose) write (uout,1) xa1, ya1, za1
  ncpy = 0_ip
  do i = 1,nopsym
    xx = opsym(i,1,1)*xa1 + opsym(i,1,2)*ya1 + opsym(i,1,3)*za1
    yy = opsym(i,2,1)*xa1 + opsym(i,2,2)*ya1 + opsym(i,2,3)*za1
    zz = opsym(i,3,1)*xa1 + opsym(i,3,2)*ya1 + opsym(i,3,3)*za1
    newcpy = .true.
    j = 1_ip
    do while (newcpy .and. j.le.ncpy)
      if (abs(xcpy(j,1)-xx) .le. TOLdist .and. &
          abs(xcpy(j,2)-yy) .le. TOLdist .and. &
          abs(xcpy(j,3)-zz) .le. TOLdist) then
        newcpy = .false.
      else
        j = j + 1_ip
      end if
    end do
    if (newcpy) then
      ncpy = ncpy + 1
      xcpy(ncpy,1) = xx
      xcpy(ncpy,2) = yy
      xcpy(ncpy,3) = zz
      if (verbose) then
        if (opproper(i)) then
          write (uout,5)  i, xx, yy, zz
        else
          write (uout,5) -i, xx, yy, zz
        end if
      end if
    end if
  end do

1 format (/                                                          &
  1x, '++SYMREPLICATE: Different symmetrical images of the point ',  &
  3f12.8)
5 format (1x, 'Opsym: ', i5, '  Image: ', 3f18.12) 

end subroutine replicatesym
