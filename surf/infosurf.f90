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
subroutine info ()

  use iso_fortran_env, only: uout=>output_unit
  use mod_prec, only: ip
  use mod_io, only: string, mline
  use mod_surf, only: rmaxsurf, epsiscp, epsilon, inuc, &
                      ntrial, rprimer, nangleb, steeper
  implicit none

  character(len=mline), dimension(5) :: rqudstr
  character(len=mline), dimension(5) :: ssteeper

  rqudstr(1) = 'Gauss-Legendre'
  rqudstr(2) = 'Clenshaw-Curtis'
  rqudstr(3) = 'Gauss-Chebychev 1st kind'
  rqudstr(4) = 'Gauss-Chebychev 2nd kind'
  rqudstr(5) = 'Perez-Jorda (Gauss-Chebychev) 2nd kind'

  ssteeper(1) = 'Runge-Kutta-Cash-Karp'
  ssteeper(2) = 'Calvo-Montijano-Randez'
  ssteeper(3) = 'Dormand-Prince method'

  write (uout,'(1x,a,1x,e13.6)') string('# Rmaxsur ='), rmaxsurf 
  write (uout,'(1x,a,1x,a)') string('# Steeper ='), string(ssteeper(steeper))
  write (uout,'(1x,a,1x,e13.6)') string('# Surface precision ='), epsilon
  write (uout,'(1x,a,1x,e13.6)') string('# EPSISCP parameter ='), epsiscp
  write (uout,'(1x,a,1x,i0)') string('# Ntrial ='), ntrial
  write (uout,'(1x,a,1x,e13.6)') string('# Rprimer ='), rprimer
! logical, allocatable, dimension(:), public :: lstart
! integer(kind=ip), allocatable, dimension(:), public :: nrsearch
! real(kind=rp), allocatable, dimension(:,:), public :: rstart
  write (uout,'(1x,a)') string('# Angular Quadratures')
  if (nangleb(1).eq.1) then
    write (uout,'(1x,a,1x,i0,1x,a,1x,i0)') &
    string('# Atom'), inuc, string('lebedev points'), nangleb(2)
  else if (nangleb(1).eq.0) then
    write (uout,'(1x,a)') string('# Phi quadrature is always trapezoidal')
    write (uout,'(1x,a,1x,i0,1x,a,1x,i0,1x,i0,1x,a)') &
    string('# Atom'), inuc, string('(ntheta,nphi,iqudt'), nangleb(2), nangleb(3), &
    string(rqudstr(nangleb(4)))
  end if

end subroutine info
