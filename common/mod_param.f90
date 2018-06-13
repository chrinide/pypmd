module mod_param

  use mod_prec, only: rp, ip
  implicit none
  public

  ! global parameters
  real(kind=rp), parameter :: pi = 3.14159265358979323846264338328_rp

  ! some options
  character(len=132) :: scratch
  logical :: largwr, debug

contains

  subroutine init_param ()

    implicit none

    scratch = './'
    largwr = .false.
    debug = .false.

  end subroutine init_param

end module mod_param
