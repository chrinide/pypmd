module mod_param

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

  real(kind=rp) :: d1mach_(5)
  private :: d1mach
  real(kind=rp) :: vbig
  real(kind=rp) :: vsmall 
  real(kind=rp) :: epsreal
  
contains

  subroutine init_param ()

    implicit none
    integer :: i, n, clock
    integer, dimension(:), allocatable :: seed

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

  end subroutine init_param

  real(kind=rp) function d1mach (i)
  
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
