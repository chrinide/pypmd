module mod_geom

  use mod_prec, only: rp, ip
  implicit none
  public

  real(kind=rp) :: covx
  real(kind=rp) :: cq(3), cm(3)
  real(kind=rp) :: moi(3,3), emoi(3), evmoi(3,3)

end module mod_geom
