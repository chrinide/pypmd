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
! A tanh-sinh quadrature rule.
subroutine tanh_sinh (x,w,n)

  use mod_prec, only: rp, ip
  use mod_param, only: pi
  implicit none

  integer(kind=ip) :: n, i
  real(kind=rp) :: ct, ct2, h, st, t, w(n), w_sum, x(n)

! This choice of H is only one of many.
! For our choice, the family ORDER = 1, 3, 7, 15, 31, 63, ... is nested.
! h = 16.0_rp/real(n+1,rp)
! h =  8.0_rp/real(n+1,rp)
  h =  4.0_rp/real(n+1,rp)

  do i = 1,n 
    t = real(2*i-n-1,rp)*h/2.0_rp
    ct = cosh(t)
    st = sinh(t)
    ct2 = cosh(0.5_rp*pi*st)
    x(i) = tanh(0.5_rp*pi*st)
    w(i) = 0.5_rp*pi*h*ct/ct2/ct2
  end do

  ! Normalize the weights so that they sum to 2.0.
  w_sum = sum(w(1:n))
  w(1:n) = 2.0_rp*w(1:n)/w_sum

  ! Reverse the data.
  x(1:n) = x(n:1:-1)
  w(1:n) = w(n:1:-1)

end subroutine

subroutine chebyshev1 (x,w,n)
  
  use mod_prec, only: rp, ip
  use mod_param, only: pi
  implicit none
  real(kind=rp), parameter :: half = 0.5_rp

  integer(kind=ip) :: n, i
  real(kind=rp) :: x(n),w(n)
  real(kind=rp) :: z

  do i = 1,n
    z = pi*(real(i,rp)-half)/real(n,rp)
    x(i) = cos(z)
    w(i) = (pi/real(n,rp))*sin(z)
  end do

end subroutine

subroutine chebyshev2 (x,w,n)
  
  use mod_prec, only: rp, ip
  use mod_param, only: pi
  implicit none

  integer(kind=ip) :: n, i
  real(kind=rp) :: x(n),w(n)
  real(kind=rp) :: xad, xa

  xad = pi/real(n+1,rp)
  do i = 1,n
    xa = i*xad
    w(i) = xad*sin(xa)
    x(i) = cos(xa)
  end do

end subroutine

! This routine computes the points and weights of the Gauss-Chebyshev
! quadrature of the first kind using the Perez-Jorda transformation.
subroutine pjt (x,w,n)

  use mod_prec, only: rp, ip
  use mod_param, only: pi
  implicit none
  real(kind=rp), parameter :: twothird = 2.0_rp/3.0_rp

  integer(kind=ip) :: n, i
  real(kind=rp) :: x(n),w(n)
  real(kind=rp) :: z,st,ct,dn

  do i = 1,n
    dn = real(n+1_ip,rp)
    z = i*pi/dn
    st = sin(z)
    ct = cos(z)
    w(i) = 16.0_rp/3.0_rp/dn*st*st*st*st
    x(i) = 1.0_rp-real(2*i,rp)/dn+2.0_rp/pi*(1.0_rp+twothird*st*st)*ct*st
  end do

end subroutine

! The abscissas for the rule of order N can be regarded as the cosines of
! equally spaced angles between 180 and 0 degrees:
!   X(I) = cos((N-I)*PI/(N-1))
! except for the basic case N = 1, when
!   X(1) = 0.
! A Clenshaw-Curtis rule that uses N points will integrate exactly all 
! polynomials of degrees 0 through N-1. If N is odd, then by symmetry the 
! polynomial of degree N will also be integrated exactly. If the value of N 
! is increased in a sensible way, then the new set of abscissas will include 
! the old ones. One such sequence would be N(K) = 2*K+1 for K = 0, 1, 2, ...
subroutine clenshaw_curtis (x,w,n)

  use mod_prec, only: rp, ip
  use mod_param, only: pi
  implicit none

  integer(kind=ip) :: n
  integer(kind=ip) :: i, j
  real(kind=rp) :: b
  real(kind=rp) :: w(n)
  real(kind=rp) :: x(n)
  real(kind=rp) :: theta(n)

  if (n == 1) then
    x(1) = 0.0_rp
    w(1) = 2.0_rp
    return
  end if

  do i = 1,n
    theta(i) = real(n-i,rp)*pi/real(n-1,rp)
  end do
  x(1:n) = cos(theta(1:n))
  do i = 1,n
    w(i) = 1.0_rp
    do j = 1,(n-1)/2
      if (2*j == (n-1)) then
        b = 1.0_rp
      else
        b = 2.0_rp
      end if
      w(i) = w(i)-b*cos(2.0_rp*real(j,rp)*theta(i))/real(4*j*j-1,rp)
    end do
  end do

  w(1) = w(1)/real(n-1,rp)
  w(2:n-1) = 2.0_rp*w(2:n-1)/real(n-1,rp)
  w(n) = w(n)/real(n-1,rp)

  ! Reverse the data.
  x(1:n) = x(n:1:-1)
  w(1:n) = w(n:1:-1)

end subroutine

subroutine fejer1 (x,w,n)

  use mod_prec, only: rp, ip
  use mod_param, only: pi
  implicit none

  integer(kind=ip) :: n
  integer(kind=ip) :: i, j
  real(kind=rp) :: w(n)
  real(kind=rp) :: x(n)
  real(kind=rp) :: theta(n)

  if (n == 1) then
    x(1) = 0.0_rp
    w(1) = 2.0_rp
    return
  end if

  do i = 1,n
    theta(i) = real(2*(n+1-i)-1,rp)*pi/real(2*n,rp)
  end do
  x(1:n) = cos(theta(1:n))
  do i = 1,n
    w(i) = 1.0_rp
    do j = 1,(n/2)
      w(i) = w(i) - 2.0_rp*cos(2.0_rp*real(j,rp)*theta(i))/real(4*j*j-1,rp)
    end do
  end do
  w(1:n) = 2.0_rp*w(1:n)/real(n,rp)

  ! Reverse the data.
  x(1:n) = x(n:1:-1)
  w(1:n) = w(n:1:-1)

end subroutine

subroutine fejer2 (x,w,n)

  use mod_prec, only: rp, ip
  use mod_param, only: pi
  implicit none

  integer(kind=ip) :: n
  integer(kind=ip) :: i, j
  real(kind=rp) :: p
  real(kind=rp) :: w(n)
  real(kind=rp) :: x(n)
  real(kind=rp) :: theta(n)

  if (n == 1) then
    x(1) = 0.0_rp
    w(1) = 2.0_rp
    return
  else if (n == 2) then
    x(1) = -0.5_rp
    x(2) = 0.5_rp
    w(1:2) = 1.0_rp
    return
  end if

  do i = 1,n
    theta(i) = real(n+1-i,rp)*pi/real(n+1,rp)
  end do
  x(1:n) = cos(theta(1:n))
  do i = 1,n
    w(i) = 1.0_rp
    do j = 1,((n-1)/2)
      w(i) = w(i) - 2.0_rp*cos(2.0_rp*real(j,rp)*theta(i))/real(4*j*j-1,rp)
    end do
    if (2 < n) then
      p = 2.0_rp*real(((n+1)/2),rp)-1.0_rp
      w(i) = w(i) - cos((p+1.0_rp)*theta(i))/p
    end if
  end do
  w(1:n) = 2.0_rp*w(1:n)/real(n+1,rp)

  ! Reverse the data.
  x(1:n) = x(n:1:-1)
  w(1:n) = w(n:1:-1)

end subroutine

! The quadrature rule will integrate exactly all polynomials up to X**(2*N-3).
! The Lobatto rule is distinguished by the fact that both endpoints (-1 and 1) 
! are always abscissas.
subroutine lobatto (x,w,n)

  use mod_prec, only: rp, ip
  use mod_param, only: pi
  implicit none

  integer(kind=ip) :: n
  integer(kind=ip) :: i, j
  real(kind=rp) :: w(n)
  real(kind=rp) :: x(n)
  real(kind=rp) :: p(n,n)
  real(kind=rp) :: tolerance
  real(kind=rp) :: xold(n)

  if (n < 2) then
    stop 'lobatto N must be at least 2' 
  end if

  tolerance = 100.0_rp*epsilon(tolerance)
  ! Initial estimate for the abscissas is the Chebyshev-Gauss-Lobatto nodes.
  do i = 1,n
    x(i) = cos(pi*real(i-1,rp)/real(n-1,rp))
  end do
  xold(1:n) = 2.0_rp

  do while (tolerance < maxval(abs(x(1:n)-xold(1:n))))
    xold(1:n) = x(1:n)
    p(1:n,1) = 1.0_rp
    p(1:n,2) = x(1:n)
    do j = 2,n-1
      p(1:n,j+1) = (real(2*j-1,rp)*x(1:n)*p(1:n,j)+real(-j+1,rp)*p(1:n,j-1))/real(j,rp)
    end do
    x(1:n) = xold(1:n)-(x(1:n)*p(1:n,n)-p(1:n,n-1))/(real(n,rp)*p(1:n,n))
  end do
  x(1:n) = x(n:1:-1)
  w(1:n) = 2.0_rp/(real((n-1)*n,rp)*p(1:n,n)**2)

  ! Reverse the data.
  x(1:n) = x(n:1:-1)
  w(1:n) = w(n:1:-1)
        
end subroutine

! The Radau rule is distinguished by the fact that the left endpoint (-1) 
! is always an abscissa. The quadrature rule will integrate exactly all 
! polynomials up to X**(2*N-2).
subroutine radau (x,w,n)

  use mod_prec, only: rp, ip
  use mod_param, only: pi
  implicit none

  integer(kind=ip) :: n
  integer(kind=ip) :: i, j
  real(kind=rp) :: w(n)
  real(kind=rp) :: x(n)
  real(kind=rp) :: p(n,n+1)
  real(kind=rp) :: tolerance
  real(kind=rp) :: xold(n)

  tolerance = 100.0_rp*epsilon(tolerance)
  ! Initial estimate for the abscissas is the Chebyshev-Gauss-Radau nodes.
  do i = 1,n
    x(i) = -cos(2.0_rp*pi*real(i-1,rp)/real(2*n-1,rp))
  end do
  xold(1:n) = 2.0_rp

  do while (tolerance < maxval(abs(x(1:n)-xold(1:n))))
    xold(1:n) = x(1:n)
    do j = 1,n+1
      p(1,j) = (-1.0_rp)**(j-1)
    end do
    p(2:n,1) = 1.0_rp
    p(2:n,2) = x(2:n)
    do j = 2, n
      p(2:n,j+1) = (real(2*j-1,rp)*x(2:n)*p(2:n,j)+real(-j+1,rp)*p(2:n,j-1))/real(j,rp)
    end do
    x(2:n) = xold(2:n)-((1.0_rp-xold(2:n))/real(n,rp))*(p(2:n,n)+p(2:n,n+1))/(p(2:n,n)-p(2:n,n+1))
  end do
  w(1) = 2.0_rp/real(n*n,rp)
  w(2:n) = (1.0_rp-x(2:n))/(real(n,rp)*p(2:n,n))**2

  ! Reverse the data.
  x(1:n) = x(n:1:-1)
  w(1:n) = w(n:1:-1)

end subroutine

! Gauss-Legendre quadrature by Stroud-Secrest method.
! The Stroud and Secrest reference did not print a specific computer program
! for the Gauss-Legendre case. Therefore, this code is based on the
! printed code for the Gauss-Jacobi case, with ALPHA and BETA set to 0.
! This means that the LEGENDRE_SS_ROOT and LEGENDRE_SS_RECUR routines,
! while appropriate for this computation, do not use the standard
! normalization for the Legendre polynomials in which Pn(1) = 1.
! The unusual scaling does not, however, affect the location of the
! roots, which is the primary thing of interest.
! The integral:
!   integral (-1 <= x <= 1) f(x) dx
subroutine legendre_ss (x,w,n)

  use mod_prec, only: rp, ip
  implicit none

  integer(kind=ip) :: n, i
  real(kind=rp) :: c(n), cc, dp2, p1, r, w(n), x(n), xtemp

  ! Set the recursion coefficients.
  do i = 1,n
    c(i) = real((i-1)*(i-1),rp)/real((2*i-1)*(2*i-3),rp)
  end do
  cc = 2.0_rp*product(c(2:n))

  do i = 1,n
    if (i == 1) then
      r = 2.78_rp/(4.0_rp+real(n*n,rp))
      xtemp = 1.0_rp - r
    else if (i == 2) then
      r = 1.0_rp + 0.06_rp*real(n-8,rp)/real(n,rp)
      xtemp = xtemp - 4.1_rp*r*(1.0_rp-xtemp)
    else if (i == 3) then
      r = 1.0_rp + 0.22_rp*real(n-8,rp)/real(n,rp)
      xtemp = xtemp - 1.67_rp*r*(x(1)-xtemp)
    else if (i < n-1) then
      xtemp = 3.0_rp*x(i-1) - 3.0_rp*x(i-2) + x(i-3)
    else if (i == n-1) then
      r = 1.0_rp/(1.0_rp+0.639_rp*real(n-4,rp)/(1.0_rp+0.71_rp*real(n-4,rp)))
      xtemp = xtemp + r*(xtemp-x(i-2))/0.766_rp
    else if (i == n) then
      r = 1.0_rp/(1.0_rp+0.22_rp*real(n-8,rp)/real(n,rp))
      xtemp = xtemp + r*(xtemp-x(i-2))/1.67_rp
    end if
    call legendre_ss_root (xtemp, n, dp2, p1, c)
    x(i) = xtemp
    w(i) = cc/dp2/p1
  end do

end subroutine

! Given the lower and upper limits of integration x1 and x2, 
! and given n, this routine returns arrays x and w of lenght n, 
! containing the abscissas and weights of the Gauss-Legendre 
! n-point quadrature formula.
! High precision is a good idea for this routine
subroutine gauleg (x1,x2,x,w,n)

  use mod_prec, only: rp, ip
  use mod_param, only: pi
  implicit none

  ! Parameters
  real(kind=rp), parameter :: eps = 3.0d-16

  ! Arguments
  integer(kind=ip), intent(in) :: n
  real(kind=rp), intent(in) :: x1, x2
  real(kind=rp), intent(out) :: x(n), w(n)

  ! Local vars
  real(kind=rp) :: p1, p2, p3, pp, xl, xm, z, z1
  integer(kind=ip) :: i, j, m

  ! The roots are symmetric in the interval, so we only have to find
  ! half of them.
  m = (n+1_ip)/2_ip
  xm = 0.5_rp*(x2+x1)
  xl = 0.5_rp*(x2-x1)

  ! Loop over the desired roots.
  do i = 1,m
    z = cos(pi*(i-0.25_rp)/(n+0.5_rp))
1   continue
      p1 = 1.0_rp
      p2 = 0.0_rp
      ! Loop up the recurrence relation to get the Legendre polynomial 
      ! evaluated at z.
      do j = 1,n
        p3 = p2
        p2 = p1
        p1 = ((2.0_rp*j-1.0_rp)*z*p2-(j-1.0_rp)*p3)/j
      end do
      ! p1 is now the desired Legendre polynomial. We next compute pp,
      ! derivative , by a standard relation involving p2, the polyn-
      ! omial of one lower order.
      pp = n*(z*p1-p2)/(z*z-1.0_rp)
      z1 = z
      ! Newton method.
      z = z1 - p1/pp
    if (abs(z-z1).gt.eps) go to 1
    ! Scale the root to the desired interval.
    x(i) = xm - xl*z
    ! and put in its symmetric counterpart.
    x(n+1-i) = xm + xl*z
    ! compute the weight.
    w(i) = 2.0_rp*xl/((1.0_rp-z*z)*pp*pp)
    ! and its symmetric counterpart.
    w(n+1-i)=w(i)
  end do

  ! Reverse the data.
  x(1:n) = x(n:1:-1)
  w(1:n) = w(n:1:-1)

end subroutine
      
! This routine computes the points and weights of the Gen-Clurtis quadrature
subroutine genclcu (a,b,x,w,n)

  use mod_prec, only: rp, ip
  use mod_param, only: pi
  implicit none
  real(kind=rp), parameter :: half = 0.5_rp
  real(kind=rp), parameter :: one = 1.0_rp
  real(kind=rp), parameter :: two = 2.0_rp

  ! Arguments
  integer(kind=ip), intent(in) :: n
  real(kind=rp), intent(in) :: a, b
  real(kind=rp), intent(out) :: x(n), w(n) 

  ! Local vars
  integer(kind=ip) :: i, j, nn, nn11, nn12, nnh
  real(kind=rp) :: z, factor, piovernn, term, x1, x2

  x1 = (b-a)*half
  x2 = (b+a)*half
  nn = n - 1_ip
  nn11 = nn*nn
  nn12 = nn*(nn-1_ip)
  nnh = (nn-1_ip)/2_ip
  piovernn = pi/real(nn,rp)
  factor = two*(b-a)/real(nn,rp)
  x(1) = b
  x(n) = a
  do i = 1,nn-1
    z = piovernn*i
    x(i+1) = x1*cos(z) + x2
    w(i+1) = 0.0_rp
    do j = 0,nnh
      if (j.eq.0_ip) then
        term = half
      else
        term = one/(1.0_rp-4.0_rp*j*j)*cos(2.0_rp*j*z)
      end if
      w(i+1) = w(i+1) + term
    end do
    w(i+1) = w(i+1)*factor
  end do
  if (mod(n,2).eq.0) then
    w(1) = x1/nn11
    w(n) = w(1)
  else
    w(1) = x1/nn12
    w(n) = w(1)
  end if
 
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Computes a tanh quadrature rule.
subroutine tanh_rule (x,w,n,h)
 
  use mod_prec, only: rp, ip
  use mod_param, only: pi
  implicit none

  integer(kind=ip) :: n, i, ij
  real(kind=rp) :: ct, h, t, w(2*n+1), x(2*n+1)

  ij = 0_ip
  do i = -n,n
    ij = ij + 1_ip
    t = real(i,rp)*h/2.0_rp
    ct = cosh(t)
    x(ij) = tanh(t)
    w(ij) = 0.5_rp*h/ct/ct
  end do

  ! Reverse the data.
  !x(1:2*n+1) = x(2*n+1:1:-1)
  !w(1:2*n+1) = w(2*n+1:1:-1)

end subroutine

! Computes a tanh-sinh quadrature rule.
subroutine tanh_sinh_rule (x,w,n,h)

  use mod_prec, only: rp, ip
  use mod_param, only: pi
  implicit none

  integer(kind=ip) :: n, i, ij
  real(kind=rp) :: ct, h, t, w(2*n+1), x(2*n+1), ct2, st

  ij = 0_ip
  do i = -n,n
    ij = ij + 1_ip
    t = real(i,rp)*h
    ct = cosh(t)
    st = sinh(t)
    ct2 = cosh(0.5_rp*pi*st)
    x(ij) = tanh(0.5_rp*pi*st)
    w(ij) = 0.5_rp*pi*h*ct/ct2/ct2
  end do

  ! Reverse the data.
  !x(1:2*n+1) = x(2*n+1:1:-1)
  !w(1:2*n+1) = w(2*n+1:1:-1)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine legendre_ss_recur (p2, dp2, p1, x, n, c)

  use mod_prec, only: rp, ip
  implicit none

  integer(kind=ip) :: n, i
  real(kind=rp) :: c(n), dp0, dp1, dp2, p0, p1, p2, x

  p1 = 1.0_rp
  dp1 = 0.0_rp
  p2 = x
  dp2 = 1.0_rp

  do i = 2,n
    p0 = p1
    dp0 = dp1
    p1 = p2
    dp1 = dp2
    p2 = x*p1 - c(i)*p0
    dp2 = x*dp1 + p1 - c(i)*dp0
  end do

end subroutine

subroutine legendre_ss_root (x, n, dp2, p1, c)

  use mod_prec, only: rp, ip
  implicit none
  integer(kind=ip), parameter :: step_max = 30

  integer(kind=ip) :: n, step
  real(kind =rp) :: c(n), d, dp2, eps, p1, p2, x

  eps = epsilon(eps)

  do step = 1,step_max
    call legendre_ss_recur (p2, dp2, p1, x, n, c)
    d = p2/dp2
    x = x - d
    if (abs(d) <= eps*(abs(x)+1.0_rp)) then
      return
    end if
  end do

end subroutine

! Computes N as a function of H and TOL.
subroutine tanh_h_to_n (h, tol, n)

  use mod_prec, only: rp, ip
  use mod_param, only: pi
  implicit none

  integer(kind=ip) :: n
  real(kind=rp) :: ct,h,t,tol,w

  n = 0
  do
    t = real(n,rp)*h/2.0_rp
    ct = cosh(t)
    w = 0.5_rp*h/ct/ct
    if (w <= tol) then
      exit
    end if
    n = n + 1
  end do

end subroutine

! Computes H as a function of M.
subroutine tanh_m_to_h (m, h)

  use mod_prec, only: rp, ip
  implicit none

  real(kind=rp) :: h
  integer(kind=ip) :: i, m

  h = 1.0_rp
  do i = -1,m,-1
    h = h*2.0_rp
  end do
  do i = 1,m
    h = h/2.0_rp
  end do

end subroutine 

! Computes N as a function of H and TOL for tanh-sinh rule.
subroutine tanh_sinh_h_to_n (h, tol, n)

  use mod_prec, only: rp, ip
  use mod_param, only: pi
  implicit none

  integer(kind=ip) :: n
  real(kind=rp) :: ct,ct2,h,st,t,tol,w

  n = 0_ip
  do
    t = real(n,rp)*h
    ct = cosh(t)
    st = sinh(t)
    ct2 = cosh(0.5_rp*pi*st)
    w = 0.5_rp*pi*h*ct/ct2/ct2
    if (w <= tol) then
      exit
    end if
    n = n + 1_ip
  end do

end subroutine
