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

    interface
      logical function  iscp (p,nuc)
        import :: rp, ip
        integer(kind=ip), intent(out) :: nuc
        real(kind=rp), intent(in) :: p(3)
      end function  iscp
      subroutine pointr1 (p,rho,grad,gradmod)
        import rp
        real(kind=rp), intent(in) :: p(3)
        real(kind=rp), intent(out) :: grad(3)
        real(kind=rp), intent(out) :: rho
        real(kind=rp), intent(out) :: gradmod
      end subroutine
    end interface

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
  subroutine rkck (y,dydx,h,yout,yerr)
      
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

  !<##### 9 stages, 6th order
  !  This scheme is due to Calvo et al., see [1].
  !    
  !   0                 | 0
  !   2/15              | 2/15
  !   1/5               | 1/20                  3/20
  !   3/10              | 3/40                  0                      9/40
  !   14/25             | 86727015/196851553    -60129073/52624712     957436434/1378352377    83886832/147842441
  !   19/25             | -86860849/45628967    111022885/25716487     108046682/101167669     -141756746/36005461
  !   35226607/35688279 | 77759591/16096467     -49252809/6452555      -381680111/51572984     879269579/66788831
  !   1                 | 237564263/39280295    -100523239/10677940    -265574846/27330247     317978411/18988713
  !   1                 | 17572349/289262523    0                      57513011/201864250      15587306/354501571
  !  --------------------------------------------------------------------------------------------------------------
  !                     | 17572349/289262523    0                      57513011/201864250      15587306/354501571
  !                     | 15231665/510830334    0                      59452991/116050448      -28398517/122437738
  !  ...continued...
  !   0                 |
  !   2/15              |
  !   1/5               |
  !   3/10              |
  !   14/25             |
  !   19/25             | 73139862/60170633
  !   35226607/35688279 | -90453121/33722162     111179552/157155827
  !   1                 | -124494385/35453627    86822444/100138635     -12873523/724232625
  !   1                 | 71783021/234982865     29672000/180480167     65567621/127060952     -79074570/210557597    0
  !  -----------------------------------------------------------------------------------------------------------------------
  !                     | 71783021/234982865     29672000/180480167     65567621/127060952     -79074570/210557597    0
  !                     | 56673824/137010559     68003849/426673583     7097631/37564021       -71226429/583093742    1/20
  !
  ! [1] *A New Embedded Pair of Runge-Kutta Formulas of orders 5 and 6*, M. Calvo, J.I. Montijano, L. Randez, Computers & Mathematics
  ! wwith Applications, Volume 20, Issue 1, 1990, Pages 15--24, ISSN 0898-1221, http://dx.doi.org/10.1016/0898-1221(90)90064-Q.
  subroutine cmr (y,dydx,h,yout,yerr)

    implicit none
  
    ! Arguments
    real(kind=rp), intent(in) :: dydx(3)
    real(kind=rp), intent(in) :: y(3)
    real(kind=rp), intent(out) :: yerr(3)
    real(kind=rp), intent(out) :: yout(3)
    real(kind=rp), intent(in) :: h
  
    ! Local vars
    real(kind=rp), parameter :: c1 = 17572349._rp/289262523._rp  
    real(kind=rp), parameter :: c3 = 57513011._rp/201864250._rp  
    real(kind=rp), parameter :: c4 = 15587306._rp/354501571._rp  
    real(kind=rp), parameter :: c5 = 71783021._rp/234982865._rp  
    real(kind=rp), parameter :: c6 = 29672000._rp/180480167._rp  
    real(kind=rp), parameter :: c7 = 65567621._rp/127060952._rp  
    real(kind=rp), parameter :: c8 = -79074570._rp/210557597._rp 
    real(kind=rp), parameter :: dc1 = c1-15231665._rp/510830334._rp 
    real(kind=rp), parameter :: dc3 = c3-59452991._rp/116050448._rp 
    real(kind=rp), parameter :: dc4 = c4+28398517._rp/122437738._rp 
    real(kind=rp), parameter :: dc5 = c5-56673824._rp/137010559._rp 
    real(kind=rp), parameter :: dc6 = c6-68003849._rp/426673583._rp 
    real(kind=rp), parameter :: dc7 = c7-7097631._rp/37564021._rp   
    real(kind=rp), parameter :: dc8 = c8+71226429._rp/583093742._rp 
    real(kind=rp), parameter :: dc9 = -1._rp/20._rp               

    real(kind=rp), parameter :: b21 = 2._rp/15._rp
    real(kind=rp), parameter :: b31 = 1._rp/20._rp                
    real(kind=rp), parameter :: b41 = 3._rp/40._rp                
    real(kind=rp), parameter :: b51 = 86727015._rp/196851553._rp  
    real(kind=rp), parameter :: b61 = -86860849._rp/45628967._rp  
    real(kind=rp), parameter :: b71 = 77759591._rp/16096467._rp   
    real(kind=rp), parameter :: b81 = 237564263._rp/39280295._rp  
    real(kind=rp), parameter :: b91 = 17572349._rp/289262523._rp  
    real(kind=rp), parameter :: b32 = 3._rp/20._rp                
    real(kind=rp), parameter :: b52 = -60129073._rp/52624712._rp  
    real(kind=rp), parameter :: b62 = 111022885._rp/25716487._rp  
    real(kind=rp), parameter :: b72 = -49252809._rp/6452555._rp   
    real(kind=rp), parameter :: b82 = -100523239._rp/10677940._rp 
    real(kind=rp), parameter :: b43 = 9._rp/40._rp
    real(kind=rp), parameter :: b53 = 957436434._rp/1378352377._rp  
    real(kind=rp), parameter :: b63 = 108046682._rp/101167669._rp   
    real(kind=rp), parameter :: b73 = -381680111._rp/51572984._rp   
    real(kind=rp), parameter :: b83 = -265574846._rp/27330247._rp   
    real(kind=rp), parameter :: b93 = 57513011._rp/201864250._rp    
    real(kind=rp), parameter :: b54 = 83886832._rp/147842441._rp  
    real(kind=rp), parameter :: b64 = -141756746._rp/36005461._rp 
    real(kind=rp), parameter :: b74 = 879269579._rp/66788831._rp  
    real(kind=rp), parameter :: b84 = 317978411._rp/18988713._rp  
    real(kind=rp), parameter :: b94 = 15587306._rp/354501571._rp  
    real(kind=rp), parameter :: b65 = 73139862._rp/60170633._rp
    real(kind=rp), parameter :: b75 = -90453121._rp/33722162._rp  
    real(kind=rp), parameter :: b85 = -124494385._rp/35453627._rp 
    real(kind=rp), parameter :: b95 = 71783021._rp/234982865._rp  
    real(kind=rp), parameter :: b76 = 111179552._rp/157155827._rp 
    real(kind=rp), parameter :: b86 = 86822444._rp/100138635._rp  
    real(kind=rp), parameter :: b96 = 29672000._rp/180480167._rp  
    real(kind=rp), parameter :: b87 = -12873523._rp/724232625._rp
    real(kind=rp), parameter :: b97 = 65567621._rp/127060952._rp 
    real(kind=rp), parameter :: b98 = -79074570._rp/210557597._rp

    real(kind=rp) :: rho, grad(3), gradmod
    real(kind=rp) :: ak2(3), ak3(3), ak4(3), ak5(3), ak6(3), ak7(3), ak8(3), ak9(3)
   
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
   
    call pointr1 (yout,rho,grad,gradmod)
    ak2(:) = grad(:)
    yout = y + h*(b31*dydx+b32*ak2)

    call pointr1 (yout,rho,grad,gradmod)
    ak3(:) = grad(:)
    yout = y + h*(b41*dydx+b43*ak3)
  
    call pointr1 (yout,rho,grad,gradmod)
    ak4(:) = grad(:)
    yout = y + h*(b51*dydx+b52*ak2+b53*ak3+b54*ak4)
  
    call pointr1 (yout,rho,grad,gradmod)
    ak5(:) = grad(:)
    yout = y + h*(b61*dydx+b62*ak2+b63*ak3+b64*ak4+b65*ak5)

    call pointr1 (yout,rho,grad,gradmod)
    ak6(:) = grad(:)
    yout = y + h*(b71*dydx+b72*ak2+b73*ak3+b74*ak4+b75*ak5+b76*ak6)

    call pointr1 (yout,rho,grad,gradmod)
    ak7(:) = grad(:)
    yout = y + h*(b81*dydx+b82*ak2+b83*ak3+b84*ak4+b85*ak5+b86*ak6+b87*ak7)

    call pointr1 (yout,rho,grad,gradmod)
    ak8(:) = grad(:)
    yout = y + h*(b91*dydx+b93*ak3+b94*ak4+b95*ak5+b96*ak6+b97*ak7+b98*ak8)

    ! Solution and error
    call pointr1 (yout,rho,grad,gradmod)
    ak9(:) = grad(:)
    yout = y + h*(c1*dydx+c3*ak3+c4*ak4+c5*ak5+c6*ak6+c7*ak7+c8*ak8)
   
    yerr = h*(dc1*dydx+dc3*ak3+dc4*ak4+dc5*ak5+dc6*ak6+dc7*ak7+dc8*ak8+dc9*ak9)

  end subroutine

  ! Butcher table for Dormand-Prince method (ode45)
  subroutine dop45 (y,dydx,h,yout,yerr)

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
      subroutine pointr1 (p,rho,grad,gradmod)
        import rp
        real(kind=rp), intent(in) :: p(3)
        real(kind=rp), intent(out) :: grad(3)
        real(kind=rp), intent(out) :: rho
        real(kind=rp), intent(out) :: gradmod
      end subroutine
    end interface

    call pointr1 (y,rho,grad,gradmod)
    ak1(:) = dydx(:)
    yout = y + h*b21*ak1

    call pointr1 (yout,rho,grad,gradmod)
    ak2(:) = grad(:)
    yout = y + h*(b31*ak1+b32*ak2)

    call pointr1 (yout,rho,grad,gradmod)
    ak3(:) = grad(:)
    yout = y + h*(b41*ak1+b42*ak2+b43*ak3)

    call pointr1 (yout,rho,grad,gradmod)
    ak4(:) = grad(:)
    yout = y + h*(b51*ak1+b52*ak2+b53*ak3+b54*ak4)

    call pointr1 (yout,rho,grad,gradmod)
    ak5(:) = grad(:)
    yout = y + h*(b61*ak1+b62*ak2+b63*ak3+b64*ak4+b65*ak5)

    call pointr1 (yout,rho,grad,gradmod)
    ak6(:) = grad(:)
    yout = y + h*(b71*ak1+b73*ak3+b74*ak4+b75*ak5+b76*ak6)

    call pointr1 (yout,rho,grad,gradmod)
    ak7(:) = grad(:)

    ! compute forth-order solution
    y4 = y + h*(c1*ak1+c3*ak3+c4*ak4+c5*ak5+c6*ak6+c7*ak7)
    ! compute fifth-order solution 
    y5 = y + h*(d1*ak1+d3*ak3+d4*ak4+d5*ak5+d6*ak6)
    ! estimated error
    yerr = abs(y5-y4)

  end subroutine

end module mod_odeint
