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
module mod_param

  use mod_io, only: mline
  use mod_prec, only: rp, ip
  implicit none
  public

  ! global parameters
  real(kind=rp), parameter :: pi = 3.14159265358979323846264338328_rp
  real(kind=rp), parameter :: h = 6.62606896e-34 !< Planck ct. [J.s] (nist2006)
  real(kind=rp), parameter :: na = 6.02214179d23 !< Avogadro ct. [1/mol] (nist2006)
  real(kind=rp), parameter :: bohrtoa = 0.52917720859_rp !< bohr to angstrom conversion factor
  real(kind=rp), parameter :: c = 2.99792458e10 !< light speed [cm/s] (nist2006)
  real(kind=rp), parameter :: rad2deg = 180.0_rp/pi
  real(kind=rp), parameter :: deg2rad = pi/180.0_rp

  ! Enumerate for structure formats
  integer(kind=ip), parameter :: isformat_unknown = 0
  integer(kind=ip), parameter :: isformat_cube = 1
  integer(kind=ip), parameter :: isformat_crystal = 2 
  integer(kind=ip), parameter :: isformat_xyz = 3 
  integer(kind=ip), parameter :: isformat_wfn = 4 
  integer(kind=ip), parameter :: isformat_wfx = 5 

  ! some options
  character(len=132) :: scratch
  logical :: verbose, debug
  
  real(kind=rp) :: d1mach_(5)
  private :: d1mach
  real(kind=rp) :: vbig
  real(kind=rp) :: vsmall 
  real(kind=rp) :: epsreal
  
  logical :: isdata
  character(len=mline) :: filename

contains

  subroutine optsparam (var,val)

    use iso_fortran_env, only: uout=>output_unit
    use mod_io, only: equal, faterr, ferror, string
    implicit none

    character(len=*), intent(in) :: var
    logical :: val

    if (equal(var,'verbose')) then
      verbose = val
      write (uout,'(1x,a)') string('# *** Verbose mode is enabled')
    else if (equal(var,'debug')) then
      debug = val
      verbose = val
      write (uout,'(1x,a)') string('# *** Debug mode is enabled')
    else
      call ferror ('optsparam', trim(var)//' unknown option', faterr)
    end if

  end subroutine optsparam

  subroutine init_param ()

    use iso_fortran_env, only: uout=>output_unit
    use mod_io, only: string

    implicit none
    integer :: i, n, clock
    integer, dimension(:), allocatable :: seed

    scratch = './'
    verbose = .false.
    debug = .false.
    isdata = .false.

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
      write (uout,'(1x,a,1x,e13.6)') string('# (DEBUG) Smallest number :'), vsmall
      write (uout,'(1x,a,1x,e13.6)') string('# (DEBUG) Biggest number :'), vbig
      write (uout,'(1x,a,1x,e13.6)') string('# (DEBUG) Machine eps :'), epsreal
    end if

  end subroutine init_param

  real(kind=rp) function d1mach (i)
  
    use iso_fortran_env, only: uout=>output_unit
    use mod_io, only: faterr, ferror
    implicit none
    integer(kind=ip) :: i
    real(kind=rp) :: b, x

    x = 1.0_rp
    b = radix(x)
    select case (i)
      case (1)
        d1mach = b**(minexponent(x)-1) ! the smallest positive magnitude.
      case (2)
        d1mach = huge(x)               ! the largest magnitude.
      case (3)
        d1mach = b**(-digits(x))       ! the smallest relative spacing.
      case (4)
        d1mach = b**(1-digits(x))      ! the largest relative spacing.
      case (5)
        d1mach = log10(b)
      case default
        call ferror ('d1mach', 'i out of bounds', faterr)
    end select

    return

  end function

end module mod_param
