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
subroutine optssurf(var,rval,ival)

  use iso_fortran_env, only: uout=>output_unit
  use mod_prec, only: rp, ip
  use mod_io, only: equal, faterr, ferror, string
  use mod_surf, only: steeper, ntrial, rprimer, epsilon, &
                      epsiscp, rmaxsurf
  use mod_param, only: verbose
  implicit none

  character(len=*), intent(in) :: var
  real(kind=rp), optional :: rval
  integer(kind=ip), optional :: ival

  if (equal(var,'steeper')) then
    steeper = abs(ival)
    if (verbose) then
      write (uout,'(1x,a,1x,i0)') string('# *** Variable steeper changed to :'), steeper
    end if
  else if (equal(var,'ntrial')) then
    ntrial = abs(ival)
    if (verbose) then
      write (uout,'(1x,a,1x,i0)') string('# *** Variable ntrial changed to :'), ntrial
    end if
  else if (equal(var,'rmaxsurf')) then
    rmaxsurf = abs(rval)
    if (verbose) then
      write (uout,'(1x,a,1x,e13.6)') string('# *** Variable rmaxsurf changed to :'), rmaxsurf
    end if
  else if (equal(var,'epsilon')) then
    epsilon = abs(rval)
    if (verbose) then
      write (uout,'(1x,a,1x,e13.6)') string('# *** Variable epsilon changed to :'), epsilon
    end if
  else if (equal(var,'epsiscp')) then
    epsiscp = abs(rval)
    if (verbose) then
      write (uout,'(1x,a,1x,e13.6)') string('# *** Variable epsiscp changed to :'), epsiscp
    end if
  else if (equal(var,'rprimer')) then
    rprimer = abs(rval)
    if (verbose) then
      write (uout,'(1x,a,1x,e13.6)') string('# *** Variable rprimer changed to :'), rprimer
    end if
  else
    call ferror ('optssurf', 'unknown option', faterr)
  end if

end subroutine optssurf
