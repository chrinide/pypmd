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
subroutine optswfn(var,val)

  use iso_fortran_env, only: uout=>output_unit
  use mod_prec, only: rp, ip
  use mod_wfn, only: cuttz, epsocc, rmaxatom, epsortho 
  use mod_io, only: equal, faterr, ferror, string
  use mod_param, only: verbose
  implicit none

  character(len=*), intent(in) :: var
  real(kind=rp) :: val

  if (equal(var,'cuttz')) then
    cuttz = abs(val)
    if (verbose) then
      write (uout,'(1x,a,1x,e13.6)') string('# *** Variable cuttz changed to :'), cuttz
    end if
  else if (equal(var,'epsocc')) then
    epsocc = abs(val)
    if (verbose) then
      write (uout,'(1x,a,1x,e13.6)') string('# *** Variable epsocc changed to :'), epsocc
    end if
  else if (equal(var,'epsortho')) then
    epsortho = abs(val)
    if (verbose) then
      write (uout,'(1x,a,1x,e13.6)') string('# *** Variable epsortho changed to :'), epsortho
    end if
  else if (equal(var,'rmaxatom')) then
    rmaxatom = abs(val)
    if (verbose) then
      write (uout,'(1x,a,1x,e13.6)') string('# *** Variable rmaxatom changed to :'), rmaxatom
    end if
  else
    call ferror ('optswfn', 'unknown option '//string(var), faterr)
  end if

end subroutine optswfn
