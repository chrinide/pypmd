module mod_param

  use mod_prec, only: rp, ip, rp_size
  implicit none
  public

  ! global parameters
  real(kind=rp), parameter :: pi = 3.14159265358979323846264338328_rp

  ! some options
  logical :: largwr, debug, bzip
  integer(kind=ip) :: reclength, bytespow

contains

  subroutine init_param ()

    implicit none

    reclength = rp_size
    bytespow = 63
    bzip = .false.
    largwr = .false.
    debug = .false.

  end subroutine init_param

end module mod_param
