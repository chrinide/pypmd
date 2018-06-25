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
subroutine findnuc ()

  use iso_fortran_env, only: uout=>output_unit
  use mod_prec, only: rp, ip
  use mod_io, only: string, faterr, ferror
  use mod_surf, only: xyzrho
  use mod_wfn, only: xyz, charge, ncent
  use mod_param, only: verbose
  implicit none

  integer(kind=ip) :: i
  real(kind=rp) :: rho, grad(3), gradmod, p(3)
  logical :: inf

  interface
    subroutine gradrho (xpoint,hi,iup,inf)
      import rp, ip
      logical, intent(inout) :: inf
      integer(kind=ip), intent(in) :: iup
      real(kind=rp), intent(in) :: hi
      real(kind=rp), intent(inout) :: xpoint(3)
    end subroutine
    subroutine pointr1 (p,rho,grad,gradmod)
      import rp
      real(kind=rp), intent(in) :: p(3)
      real(kind=rp), intent(out) :: grad(3)
      real(kind=rp), intent(out) :: rho
      real(kind=rp), intent(out) :: gradmod
    end subroutine
  end interface

  do i = 1,ncent
    p(:) = xyz(i,:)
    call gradrho (p,0.05_rp,1,inf)
    call pointr1 (p,rho,grad,gradmod)
    if (gradmod.gt.1d-4) then
      if (charge(i).gt.2.0_rp) then
        if (verbose) write (uout,321) i
        xyzrho(i,:) = xyz(i,:)
      else
        call ferror('findnuc', 'failed finding nucleus '//string(i), faterr)
      end if
    else
      xyzrho(i,:) = xyz(i,:)
    end if
  end do
 
 321 format (1x,'# Assuming nuclei ',i0,' position: Check!')

end subroutine
