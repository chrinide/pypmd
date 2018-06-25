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
subroutine optsparam(var,val)

  use iso_fortran_env, only: uout=>output_unit
  use mod_io, only: equal, faterr, ferror, string
  use mod_param, only: verbose, debug
  implicit none

  character(len=*), intent(in) :: var
  logical :: val

  if (equal(var,'verbose')) then
    verbose = val
    write (uout,'(1x,a)') string('# *** Verbose mode is enabled')
  else if (equal(var,'debug')) then
    debug = val
    verbose = val
    write (uout,'(1x,a)') string('# !! WARNING: DEBUG MODE ENABLED !!')
    write (uout,'(1x,a)') string('# !! Lot of info will be printed !!')
    write (uout,'(1x,a)') string('# !! WARNING: DEBUG MODE ENABLED !!')
    write (uout,*)
  else
    call ferror ('optsparam', 'unknown option', faterr)
  end if

end subroutine optsparam
