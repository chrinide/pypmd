module mod_param

  use mod_prec, only: rp, ip
  implicit none
  public

  ! global parameters
  real(kind=rp), parameter :: pi = 3.14159265358979323846264338328_rp

  ! some options
  character(len=132) :: scratch
  logical :: verbose, debug

  logical :: isdata

end module mod_param
