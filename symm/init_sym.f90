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
subroutine init_sym()

  use mod_prec, only: rp, ip
  use mod_sym, only: toldist, toleqvm, tolsng, tolisint, &
                     tolnull, toleigen, toldirty, toltriplet, &
                     errsng, erreigen, inf_order, mol_linear, &
                     mol_planar
  implicit none
    
  ! default values  
  toldist = 3e-5_rp
  toleqvm = 2.0_rp*toldist
  tolsng = 1e-5_rp
  tolisint = 3e-5_rp
  tolnull = 3e-6_rp
  toleigen = 3e-5_rp
  toldirty = .false.
  toltriplet = 0.1_rp
  errsng = 0.0_rp
  erreigen = 0.0_rp
  inf_order = 16_ip
  mol_linear = .false.
  mol_planar = .false.


end subroutine init_sym
