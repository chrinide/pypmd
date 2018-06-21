module mod_param

  use mod_prec, only: rp, ip
  implicit none
  public

  ! global parameters
  real(kind=rp), parameter :: VBIG = huge(1) !< a very big number
  real(kind=rp), parameter :: VSMALL = 1d-80   !< a very small number
  real(kind=rp), parameter :: pi = 3.14159265358979323846264338328_rp
  real(kind=rp), parameter :: h = 6.62606896d-34 !< Planck ct. [J.s] (nist2006)
  real(kind=rp), parameter :: na = 6.02214179d23 !< Avogadro ct. [1/mol] (nist2006)
  real(kind=rp), parameter :: bohrtoa = 0.52917720859d0 !< bohr to angstrom conversion factor
  real(kind=rp), parameter :: c = 2.99792458d10 !< light speed [cm/s] (nist2006)
  real(kind=rp), parameter :: rad2deg = 180.0_rp/pi
  real(kind=rp), parameter :: deg2rad = pi/180.0_rp

  ! Enumerate for structure formats
  logical :: isdata
  integer(kind=ip), parameter :: isformat_unknown = 0
  integer(kind=ip), parameter :: isformat_cube = 1
  integer(kind=ip), parameter :: isformat_crystal = 2 
  integer(kind=ip), parameter :: isformat_xyz = 3 
  integer(kind=ip), parameter :: isformat_wfn = 4 
  integer(kind=ip), parameter :: isformat_wfx = 5 
  integer(kind=ip), parameter :: isformat_fchk = 6 
  integer(kind=ip), parameter :: isformat_molden = 7 

  ! some options
  character(len=132) :: scratch
  logical :: verbose, debug
  
  real(kind=rp) :: d1mach_(5)
  public :: d1mach

contains

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
