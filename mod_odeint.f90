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
 
    use mod_io, only: faterr, ferror, string
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
    
    ! Local vars
    integer(kind=ip) :: nstp, nuc
    real(kind=rp), parameter :: hmin = 0.0_rp ! minimal step size
    real(kind=rp) :: x1 ! intial point
    real(kind=rp) :: x2 ! final point
    real(kind=rp) :: p(3) ! initial solution point
    real(kind=rp) :: h ! initial step size
    real(kind=rp) :: x ! update point
    real(kind=rp) :: hnext ! next steep size
    real(kind=rp) :: dydx(3), y(3), yscal(3)
    real(kind=rp) :: a1, a2, a3
    real(kind=rp) :: rho, grad(3), gradmod

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
 
    inf = .false.
    p(1) = ystart(1)
    p(2) = ystart(2)
    p(3) = ystart(3)
    call field (p,rho,grad,gradmod)
    if (gradmod.lt.epsg .and. rho.lt.epsg1) then
      inf = .true.
      return
    end if
 
    x1 = 0.0_rp
    x2 = 1d40*iup
    x = x1 ! initial point
    !h = sign(h1,x2-x1) ! initial steep size
    h = min(h1,x2-x1) ! initial steep size
    y(:) = ystart(:) ! initial point 
 
    do nstp = 1,maxstp
      call field (y,rho,grad,gradmod)
      dydx(:) = grad(:)
      yscal(:) = max(abs(y(:))+abs(h*dydx(:))+tiny,eps)
      if ((x+h-x2)*(x+h-x1).gt.0.0_rp) h = x2 - x
      call rkqs (y,dydx,x,h,eps,yscal,hnext,steeper,field)
      if ((x-x2)*(x2-x1).ge.0.0_rp .or. check(y,nuc)) then
        ystart(:) = y(:)
        return
      end if
      if (abs(hnext).lt.hmin) then
        call ferror ('mod_odeint/odeint', 'stepsize small than minimum', faterr)
      end if
      if (nstp.eq.maxstp) then
        call ferror ('mod_odeint/odeint', 'reached maxstp', faterr)
      end if 
      h = hnext
    end do

    ! Test if the point is far from RMAXSURF from current atom. 
    a1 = y(1) - xnuc(1)
    a2 = y(2) - xnuc(2)
    a3 = y(3) - xnuc(3)
    if ((a1*a1+a2*a2+a3*a3).ge.5d0*5d0) then
      inf = .true.
      return
    else
      call ferror ('mod_odeint/odeint', 'Non nuclear maximum at : ' &
                                         //string(y(1),'e')//' '    &  
                                         //string(y(2),'e')//' '    &  
                                         //string(y(3),'e'), faterr) 
    end if

  end subroutine

  subroutine rkqs (y,dydx,x,htry,eps,yscal,hnext,steeper,field)
      
    use mod_io, only: faterr, ferror
    implicit none

    ! Parameters
    !real(kind=rp), parameter :: safety = 0.9_rp
    !real(kind=rp), parameter :: pgrow = -0.2_rp
    !real(kind=rp), parameter :: pshrnk = -0.25_rp
    real(kind=rp), parameter :: errcon = 1.89d-4
    
    ! Arguments
    integer(kind=ip), intent(in) :: steeper
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
    real(kind=rp) :: pshrnk, safety, pgrow

    interface
      subroutine field (p,rho,grad,gradmod)
        import rp
        real(kind=rp), intent(in) :: p(3)
        real(kind=rp), intent(out) :: grad(3)
        real(kind=rp), intent(out) :: rho
        real(kind=rp), intent(out) :: gradmod
      end subroutine
    end interface
      
    h = htry
    hnext = 0.0_rp
    errmax = 0.0_rp
    yerr = 0.0_rp
    if (steeper.eq.1 .or. steeper.eq.3) then
      pshrnk = -1.0_rp/4.0_rp
      safety = 0.9_rp
      pgrow = -1.0_rp/5.0_rp
    else if (steeper.eq.2) then
      pshrnk = -1.0_rp/6.0_rp
      safety = 0.8_rp
      pgrow = -1.0_rp/9.0_rp
    !else if (steeper.eq.3) then
    !  pshrnk = -1.0_rp/10.0_rp
    !  safety = 0.8_rp
    !  pgrow = -1.0_rp/17.0_rp
    end if

    do
      if (steeper.eq.1) then
        call rkck (y,dydx,h,ytemp,yerr,field)
      else if (steeper.eq.2) then
        !call cmr (y,dydx,h,ytemp,yerr)
      else if (steeper.eq.3) then
        call dop45 (y,dydx,h,ytemp,yerr,field)
      end if
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

  ! Runge-Kutta-Cash-Karp embedded 4(5)-order, with local extrapolation.
  ! Runge-Kutta embedded 4th order, Cash-Karp parametrization.
  ! ##### 6 stages, 5th order
  ! This scheme is due to Cash and Karp, see [1].
  !    
  !   0    | 0
  !   1/5	 | 1/5
  !   3/10 | 3/40	         9/40
  !   3/5	 | 3/10	         -9/10	      6/5
  !   1	   | -11/54	       5/2	        -70/27	    35/27
  !   7/8	 | 1631/55296    175/512      575/13824   44275/110592     253/4096     0
  !  ----------------------------------------------------------------------------------------
  !        | 37/378        0           250/621      125/594          0            512/1771
  !        | 2825/27648    0           18575/48384  13525/55296      277/14336    1/4
  ! 
  ! [1] *A variable order Runge-Kutta method for initial value problems with rapidly varying right-hand sides*, J. R. Cash,
  ! A. H. Karp, ACM Transactions on Mathematical Software, vol. 16,  pp. 201--222, 1990, doi:10.1145/79505.79507.
  subroutine rkck (y,dydx,h,yout,yerr, field)
      
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
    real(kind=rp), parameter :: b21=0.2_rp,                           &
                                b31=3.0_rp/40.0_rp,                    &
                                b32=9.0_rp/40.0_rp,                    &
                                b41=0.3_rp,b42=-0.9_rp,b43=1.2_rp,      &
                                b51=-11.0_rp/54.0_rp,b52=2.5_rp,        &
                                b53=-70.0_rp/27.0_rp,b54=35._rp/27.0_rp, &
                                b61=1631.0_rp/55296.0_rp,              &
                                b62=175.0_rp/512.0_rp,                 &
                                b63=575.0_rp/13824.0_rp,               &
                                b64=44275.0_rp/110592.0_rp,            &
                                b65=253.0_rp/4096.0_rp                  
    real(kind=rp), parameter :: c1=37.0_rp/378.0_rp,c3=250.0_rp/621.0_rp,&
                                c4=125.0_rp/594.0_rp,                  &
                                c6=512.0_rp/1771.0_rp                   
    real(kind=rp), parameter :: dc1=c1-2825.0_rp/27648.0_rp,           &
                                dc3=c3-18575.0_rp/48384.0_rp,          &
                                dc4=c4-13525.0_rp/55296.0_rp,          &
                                dc5=-277.0_rp/14336.0_rp,              &
                                dc6=c6-0.25_rp                         
   
    interface
      subroutine pointr1 (p,rho,grad,gradmod)
        import rp
        real(kind=rp), intent(in) :: p(3)
        real(kind=rp), intent(out) :: grad(3)
        real(kind=rp), intent(out) :: rho
        real(kind=rp), intent(out) :: gradmod
      end subroutine
    end interface

    yout = y + b21*h*dydx
   
    call field (yout,rho,grad,gradmod)
    ak2(:) = grad(:)
    yout = y + h*(b31*dydx+b32*ak2)

    call field (yout,rho,grad,gradmod)
    ak3(:) = grad(:)
    yout = y + h*(b41*dydx+b42*ak2+b43*ak3)
  
    call field (yout,rho,grad,gradmod)
    ak4(:) = grad(:)
    yout = y + h*(b51*dydx+b52*ak2+b53*ak3+b54*ak4)
  
    call field (yout,rho,grad,gradmod)
    ak5(:) = grad(:)
    yout = y + h*(b61*dydx+b62*ak2+b63*ak3+b64*ak4+b65*ak5)

    call field (yout,rho,grad,gradmod)
    ak6(:) = grad(:)
    yout = y + h*(c1*dydx+c3*ak3+c4*ak4+c6*ak6)
   
    yerr = h*(dc1*dydx+dc3*ak3+dc4*ak4+dc5*ak5+dc6*ak6)
 
  end subroutine

end module mod_odeint
