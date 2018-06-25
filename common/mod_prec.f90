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
module mod_prec

  use iso_c_binding, only: c_int, c_double, c_size_t
  implicit none
  save
  private

  public :: rp, rp_size
  public :: ip, ip_size
  public :: size_t, size_t_size
  public :: smallest
  public :: print_kind_info

  integer, parameter :: ip = c_int
  integer, parameter :: ip_size = bit_size(int(0, ip))/8
  integer, parameter :: size_t = c_size_t
  integer, parameter :: size_t_size = bit_size(int(0, size_t))/8 
  integer, parameter :: rp = c_double
  integer, parameter :: rp_size = bit_size(int(0, rp))/8  
  
  real(kind=rp), parameter :: smallest = epsilon(0.0_rp)

contains

  subroutine print_kind_info(uout)

    implicit none
    integer(kind=ip), intent(in) :: uout
 
    write (uout,'(t2,a)') '# '
    write (uout,'(t2,a)') '# Data type information '
    write (uout,'(t2,a,/,t2,a,2(/,t2,a,i2),3(/,t2,a,E15.8))') &
      '# REAL ', &
      '#  Data type name : rp', &
      '#  Kind value : ', kind(0.0_rp), &
      '#  Precision : ', precision(0.0_rp), &
      '#  Smallest nonnegligible quantity relative to 1 : ', epsilon(0.0_rp), &
      '#  Smallest positive number : ', tiny(0.0_rp), &
      '#  Largest representable number : ', huge(0.0_rp)
    write (uout,'(t2,a,/,t2,a,2(/,t2,a,i2),t2,a,i20)') &
      '# INTEGER ', &
      '#  Data type name : ip', &
      '#  Kind value : ', kind(0_ip), &
      '#  Bit size : ', bit_size(0_ip), &
      '#  Largest representable number : ', huge(0_ip)
    write (uout,'(t2,a,/,t2,a,2(/,t2,a,i2),t2,a,i20)') &
      '# INTEGER ', &
      '#  Data type name : size_t', &
      '#  Kind value : ', kind(0_size_t), &
      '#  Bit size : ', bit_size(0_size_t), &
      '#  Largest representable number : ', huge(0_size_t)
    write (uout,'(t2,a,/,t2,a,/,t2,a,i2)') &
      '# LOGICAL ', &
      '#  Data type name : default', &
      '#  Kind value : ', kind(.true.)
    write (uout,'(t2,a,/,t2,a,/,t2,a,i2)') &
      '# CHARACTER ', &
      '#  Data type name : default', &
      '#  Kind value : ', kind('c')
 
  end subroutine print_kind_info

end module mod_prec
