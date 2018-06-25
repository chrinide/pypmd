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
subroutine init_param ()

  use iso_fortran_env, only: uout=>output_unit
  use mod_io, only: string
  use mod_prec, only: rp
  use mod_param, only: scratch, verbose, debug, vbig, &
                       isdata, d1mach, d1mach_, debug, &
                       vsmall, epsreal

  implicit none
  integer :: i, n, clock
  integer, dimension(:), allocatable :: seed

  isdata = .false.
  scratch = './'
  verbose = .false.
  debug = .false.

  ! random seed
  call random_seed(size=n)
  allocate(seed(n))
  call system_clock(count=clock)
  seed = clock + 37*(/(i-1, i =1,n)/)
  call random_seed(put=seed)
  deallocate(seed)

  ! machine constants
  do i = 1,5
    d1mach_(i) = d1mach(i)
  end do
  vsmall = d1mach_(1)
  vbig = d1mach_(2)
  epsreal = epsilon(1.0_rp)

  if (debug) then
    write (uout,'(1x,a,1x,e13.6)') string('# (DEBUG) smallest number :'), vsmall
    write (uout,'(1x,a,1x,e13.6)') string('# (DEBUG) biggest number :'), vbig
    write (uout,'(1x,a,1x,e13.6)') string('# (DEBUG) machine eps :'), epsreal
  end if

end subroutine init_param
