module mod_steeper      
  
  use mod_prec, only: rp, ip
  implicit none
  private

  public :: rkqs

contains

  subroutine rkqs (y,dydx,x,htry,eps,yscal,hnext)
      
    use mod_io, only: ferror, faterr
    implicit none

    ! Parameters
    real(kind=rp), parameter :: safety = 0.9_rp
    real(kind=rp), parameter :: pgrow = -0.2_rp
    real(kind=rp), parameter :: pshrnk = -0.25_rp
    real(kind=rp), parameter :: errcon = 1.89d-4
    
    ! Arguments
    real(kind=rp), dimension(3), intent(in) :: dydx
    real(kind=rp), dimension(3), intent(inout) :: y
    real(kind=rp), dimension(3), intent(in) :: yscal
    real(kind=rp), intent(inout) :: x
    real(kind=rp), intent(in) :: eps
    real(kind=rp), intent(in) :: htry
    real(kind=rp), intent(out) :: hnext

    ! Local vars
    integer(kind=ip) :: i
    real(kind=rp), dimension(3) :: yerr, ytemp
    real(kind=rp) :: h, errmax, htemp, xnew
      
    h = htry
    hnext = 0.0_rp
    errmax = 0.0_rp
    yerr = 0.0_rp

    do
      call rkck (y,dydx,h,ytemp,yerr)
      ! adaptive and error estimation 
      errmax = 0.0_rp
      do i = 1,3
        errmax = max(errmax,abs(yerr(i)/yscal(i)))
      end do
      errmax = errmax/eps
      if (errmax.gt.1.0_rp) then
        htemp = safety*h*(errmax**pshrnk)
        !h = sign(max(abs(htemp),0.1_rp*abs(h)),h)
        h = min(max(abs(htemp),0.1_rp*abs(h)),h)
        xnew = x + h
        if (xnew.eq.x) then
          call ferror ('mod_odeint/rkqs', 'stepsize underflow', faterr)
          return
        end if
        cycle
      else
        if (errmax.gt.errcon) then
          hnext = safety*h*(errmax**pgrow)
        else
          hnext = 5.0_rp*h
        end if
        x = x + h
        y = ytemp
        return
      end if
    end do
 
  end subroutine rkqs

  ! Runge-Kutta-Cash-Karp
  subroutine rkck (y,dydx,h,yout,yerr)
      
    use mod_fields, only: pointr1
    implicit none
  
    ! Arguments
    real(kind=rp), intent(in) :: dydx(3)
    real(kind=rp), intent(in) :: y(3)
    real(kind=rp), intent(out) :: yerr(3)
    real(kind=rp), intent(out) :: yout(3)
    real(kind=rp), intent(in) :: h
  
    ! Local vars
    real(kind=rp) :: rho, grad(3), gradmod
    real(kind=rp) :: ak2(3), ak3(3), ak4(3), ak5(3), ak6(3)
    real(kind=rp), parameter :: b21=0.2d0,                           &
                                b31=3.0d0/40.0d0,                    &
                                b32=9.0d0/40.0d0,                    &
                                b41=0.3d0,b42=-0.9d0,b43=1.2d0,      &
                                b51=-11.0d0/54.0d0,b52=2.5d0,        &
                                b53=-70.0d0/27.0d0,b54=35.d0/27.0d0, &
                                b61=1631.0d0/55296.0d0,              &
                                b62=175.0d0/512.0d0,                 &
                                b63=575.0d0/13824.0d0,               &
                                b64=44275.0d0/110592.0d0,            &
                                b65=253.0d0/4096.0d0                  
    real(kind=rp), parameter :: c1=37.0d0/378.0d0,c3=250.0d0/621.0d0,&
                                c4=125.0d0/594.0d0,                  &
                                c6=512.0d0/1771.0d0                   
    real(kind=rp), parameter :: dc1=c1-2825.0d0/27648.0d0,           &
                                dc3=c3-18575.0d0/48384.0d0,          &
                                dc4=c4-13525.0d0/55296.0d0,          &
                                dc5=-277.0d0/14336.0d0,              &
                                dc6=c6-0.25d0                         
   
    yout = y + b21*h*dydx
   
    call pointr1 (yout,rho,grad,gradmod)
    ak2(:) = grad(:)
    yout = y + h*(b31*dydx+b32*ak2)

    call pointr1 (yout,rho,grad,gradmod)
    ak3(:) = grad(:)
    yout = y + h*(b41*dydx+b42*ak2+b43*ak3)
  
    call pointr1 (yout,rho,grad,gradmod)
    ak4(:) = grad(:)
    yout = y + h*(b51*dydx+b52*ak2+b53*ak3+b54*ak4)
  
    call pointr1 (yout,rho,grad,gradmod)
    ak5(:) = grad(:)
    yout = y + h*(b61*dydx+b62*ak2+b63*ak3+b64*ak4+b65*ak5)

    call pointr1 (yout,rho,grad,gradmod)
    ak6(:) = grad(:)
    yout = y + h*(c1*dydx+c3*ak3+c4*ak4+c6*ak6)
   
    yerr = h*(dc1*dydx+dc3*ak3+dc4*ak4+dc5*ak5+dc6*ak6)
 
  end subroutine

end module mod_steeper
