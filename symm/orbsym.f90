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
! classify all atoms into orbits. All atoms in an orbit
! have the same atomic number and the same distance to the center.
!
subroutine orbsym()

  use iso_fortran_env, only: uout=>output_unit
  use mod_io, only: ferror, faterr
  use mod_prec, only: ip, rp
  use mod_param, only: debug
  use mod_sym, only: matorb, norbit, morbit, toldist, molradius, &
                     orbdis, orbz, orbmol, natorb, iatorb, sdim, &
                     atzmol, ax, ay, az
  implicit none
 
  real(kind=rp) :: dis2
  integer(kind=ip) :: i, j, iorb, n

  n = sdim
  ! Classify all atoms into orbits. All atoms in an orbit have the
  ! same atomic number and the same distance to the center of mass:
  norbit = 0
  molradius = 0.0_rp
  do i = 1,n
    dis2 = sqrt(ax(i)*ax(i) + ay(i)*ay(i) + az(i)*az(i))
    if (dis2 .gt. molradius) molradius = dis2
    iorb = -1
    j = 1
    do while (j.le.norbit .and. iorb.lt.0)
      if (abs(dis2-orbdis(j)) .le. TOLdist .and. atZmol(i) .eq. orbZ(j)) then
        iorb = j
      end if
      j = j + 1
    end do
    if (iorb .lt. 0) then
      norbit = norbit + 1
      if (norbit .gt. MORBIT) then
        call ferror ('orbsym','too many orbits', faterr)
      end if
      orbdis(norbit) = dis2
      orbZ(norbit) = atZmol(i)
      natorb(norbit) = 0
      iorb = norbit
    end if
    orbmol(i) = iorb
    natorb(iorb) = natorb(iorb) + 1
    if (natorb(iorb) .gt. MATORB) then
      call ferror ('orbsym', 'too many atoms in an orbit', faterr)
    end if
    iatorb(iorb,natorb(iorb)) = i
  end do

  if (debug) then
    write (uout,550) norbit
    do i = 1,norbit
      write (uout,555) i, natorb(i), orbZ(i), orbdis(i), &
                        (iatorb(i,j),j=1,natorb(i))
    end do
  end if

550 format (/                                                    &
    1x, '++SYMORB: Determining the atomic orbits.'/              &
    1x, 'ALG(symorb) Number of orbits: ', i5)

555 format (                                                     &
    1x, 'ALG(symorb) Orbit ', i5/                                &
    1x, 'ALG(symorb) Number of atoms: ', i5/                     &
    1x, 'ALG(symorb) Atomic number and distance: ', i5, f18.9/   &
    1x, 'ALG(symorb) List of atoms in this orbit:'/              &
    (1x, 10i6))

end subroutine orbsym
