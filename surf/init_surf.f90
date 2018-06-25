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
subroutine init_surf ()

  use mod_prec, only: rp, ip
  use mod_wfn, only: ncent
  use mod_surf, only: steeper, ntrial, rprimer, inuc, &
                      epsilon, epsiscp, nangleb, &
                      allocate_space_for_cp, rmaxsurf 
  implicit none

  steeper = 1_ip
  ntrial = 11_ip
  rprimer = 0.4_rp
  inuc = 0_ip
  epsilon = 1d-5 
  epsiscp = 0.08_rp
  nangleb(1) = 1
  nangleb(2) = 434
  nangleb(3) = 0
  nangleb(4) = 0
  rmaxsurf = 10.0_rp

  call allocate_space_for_cp (ncent)

end subroutine
