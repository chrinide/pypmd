module mod_odeint      
  
  use mod_prec, only: rp, ip
  use mod_io, only: ferror, faterr, string
  implicit none
  private

  public :: odeint

contains

  subroutine odeint (ystart,h1,iup,inf,eps,xnuc,steeper)
 
    implicit none

    ! Parameters
    integer(kind=ip), parameter :: maxstp = 5000
    real(kind=rp), parameter :: tiny = 1d-40
    real(kind=rp), parameter :: epsg = 1d-15
    real(kind=rp), parameter :: epsg1 = 1d-15

    ! Arguments
    integer(kind=ip), intent(in) :: steeper
    real(kind=rp), intent(inout) :: ystart(3)
    real(kind=rp), intent(in) :: h1
    real(kind=rp), intent(in) :: iup
    logical, intent(inout) :: inf
    real(kind=rp), intent(in) :: eps
    real(kind=rp), intent(in) :: xnuc(3)

  end subroutine

end module mod_odeint
