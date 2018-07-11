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
module mod_odeint      
  
  use mod_prec, only: rp, ip
  implicit none
  private

  public :: odeint

contains

  subroutine odeint (ystart,h1,iup,inf,eps,xnuc,steeper,field,check)
 
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

    interface
      subroutine field (p,rho,grad,gradmod)
        import rp
        real(kind=rp), intent(in) :: p(3)
        real(kind=rp), intent(out) :: grad(3)
        real(kind=rp), intent(out) :: rho
        real(kind=rp), intent(out) :: gradmod
      end subroutine
      logical function  check (p,nuc)
        import :: rp, ip
        integer(kind=ip), intent(out) :: nuc
        real(kind=rp), intent(in) :: p(3)
      end function  
    end interface

  end subroutine

  ! Butcher table for Dormand-Prince method (ode45)
  subroutine dop45 (y,dydx,h,yout,yerr,field)

    implicit none
  
    ! Arguments
    real(kind=rp), intent(in) :: dydx(3)
    real(kind=rp), intent(in) :: y(3)
    real(kind=rp), intent(out) :: yerr(3)
    real(kind=rp), intent(out) :: yout(3)
    real(kind=rp), intent(in) :: h

    real(kind=rp), parameter :: b21 = 1.0_rp/5.0_rp
    real(kind=rp), parameter :: b31 = 3.0_rp/40.0_rp
    real(kind=rp), parameter :: b32 = 9.0_rp/40.0_rp
    real(kind=rp), parameter :: b41 = 44.0_rp/45.0_rp
    real(kind=rp), parameter :: b42 = -56.0_rp/15.0_rp
    real(kind=rp), parameter :: b43 = 32.0_rp/9.0_rp
    real(kind=rp), parameter :: b51 = 19372.0_rp/6561.0_rp
    real(kind=rp), parameter :: b52 = -25360.0_rp/2187.0_rp
    real(kind=rp), parameter :: b53 = 64448.0_rp/6561.0_rp
    real(kind=rp), parameter :: b54 = -212.0_rp/729.0_rp
    real(kind=rp), parameter :: b61 = 9017.0_rp/3168.0_rp
    real(kind=rp), parameter :: b62 = -355.0_rp/33.0_rp
    real(kind=rp), parameter :: b63 = 46732.0_rp/5247.0_rp
    real(kind=rp), parameter :: b64 = 49.0_rp/176.0_rp
    real(kind=rp), parameter :: b65 = -5103.0_rp/18656.0_rp
    real(kind=rp), parameter :: b71 = 35.0_rp/384.0_rp
    real(kind=rp), parameter :: b73 = 500.0_rp/1113.0_rp
    real(kind=rp), parameter :: b74 = 125.0_rp/192.0_rp
    real(kind=rp), parameter :: b75 = -2187.0_rp/6784.0_rp
    real(kind=rp), parameter :: b76 = 11.0_rp/84.0_rp

    real(kind=rp), parameter :: c1 = 5179.0_rp/57600.0_rp
    real(kind=rp), parameter :: c3 = 7571.0_rp/16695.0_rp
    real(kind=rp), parameter :: c4 = 393.0_rp/640.0_rp
    real(kind=rp), parameter :: c5 = -92097.0_rp/339200.0_rp
    real(kind=rp), parameter :: c6 = 187.0_rp/2100.0_rp
    real(kind=rp), parameter :: c7 = 1.0_rp/40.0_rp

    real(kind=rp), parameter :: d1 = 35.0_rp/384.0_rp
    real(kind=rp), parameter :: d3 = 500.0_rp/1113.0_rp
    real(kind=rp), parameter :: d4 = 125.0_rp/192.0_rp
    real(kind=rp), parameter :: d5 = -2187.0_rp/6784.0_rp
    real(kind=rp), parameter :: d6 = 11.0_rp/84.0_rp

    real(kind=rp) :: rho, grad(3), gradmod
    real(kind=rp) :: ak1(3), ak2(3), ak3(3), ak4(3), ak5(3), ak6(3), ak7(3)
    real(kind=rp) :: y4(3), y5(3)
   
    interface
      subroutine field (p,rho,grad,gradmod)
        import rp
        real(kind=rp), intent(in) :: p(3)
        real(kind=rp), intent(out) :: grad(3)
        real(kind=rp), intent(out) :: rho
        real(kind=rp), intent(out) :: gradmod
      end subroutine
    end interface

    call field (y,rho,grad,gradmod)
    ak1(:) = dydx(:)
    yout = y + h*b21*ak1

    call field (yout,rho,grad,gradmod)
    ak2(:) = grad(:)
    yout = y + h*(b31*ak1+b32*ak2)

    call field (yout,rho,grad,gradmod)
    ak3(:) = grad(:)
    yout = y + h*(b41*ak1+b42*ak2+b43*ak3)

    call field (yout,rho,grad,gradmod)
    ak4(:) = grad(:)
    yout = y + h*(b51*ak1+b52*ak2+b53*ak3+b54*ak4)

    call field (yout,rho,grad,gradmod)
    ak5(:) = grad(:)
    yout = y + h*(b61*ak1+b62*ak2+b63*ak3+b64*ak4+b65*ak5)

    call field (yout,rho,grad,gradmod)
    ak6(:) = grad(:)
    yout = y + h*(b71*ak1+b73*ak3+b74*ak4+b75*ak5+b76*ak6)

    call field (yout,rho,grad,gradmod)
    ak7(:) = grad(:)

    ! compute forth-order solution
    y4 = y + h*(c1*ak1+c3*ak3+c4*ak4+c5*ak5+c6*ak6+c7*ak7)
    ! compute fifth-order solution 
    y5 = y + h*(d1*ak1+d3*ak3+d4*ak4+d5*ak5+d6*ak6)
    ! estimated error
    yerr = abs(y5-y4)

  end subroutine

end module mod_odeint
