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
subroutine infowfn()

  use iso_fortran_env, only: uout=>output_unit
  use mod_prec, only: ip
  use mod_io, only: string
  use mod_wfn, only: ncent, charge, atnam, xyz, nvirtual, noccupied, &
                     occupied, epsocc, epsortho
  implicit none

  integer(kind=ip) :: i

  write (uout,'(1x,a,1x,i0)') string('# Number of centers :'), ncent
  atnam = adjustl(atnam)
  do i = 1,ncent
    write (uout,'(1x,a,1x,i0,1x,a,1x,f4.1,1x,3f12.6)') '#', i, &
                  atnam(i)(1:2), charge(i), xyz(i,:)
  end do
  write (uout,'(1x,a,1x,e13.6)') string('# Occupied eps ='), epsocc
  write (uout,'(1x,a,1x,i0)') string('# Number of occupied orbitals ='), noccupied
  write (uout,*) string('#'), occupied(:)
  write (uout,'(1x,a,1x,i0)') string('# Number of virtual orbitals ='), nvirtual
  write (uout,'(1x,a,1x,e13.6)') string('# epsortho = '), epsortho

end subroutine infowfn
