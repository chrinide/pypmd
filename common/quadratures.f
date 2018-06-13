      subroutine gauleg (x1,x2,x,w,n)
c
c.....Given the lower and upper limits of integration x1 and x2, 
c     and given n, this routine returns arrays x and w of lenght n, 
c     containing the abscissas and weights of the Gauss-Legendre 
c     n-point quadrature formula.
c.....High precision is a good idea for this routine
c     
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
1      continue
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
          
      end subroutine
c
c-----------------------------------------------------------------------
c
      subroutine genclcu (a,b,x,w,n)
c
c-----This routine computes the points and weights of the Gen-Clurtis quadrature
c
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
c
c-----------------------------------------------------------------------
c
      subroutine gaucheb1 (a,b,x,w,n)
c
c-----This routine computes the points and weights of the Gauss-Chebyshev
c     quadrature of the first kind, neccessary to obtain the integral 
c     F(x) dx = sum(i=1,n) W_i F(x_i)
c
      use mod_prec, only: rp, ip
      use mod_param, only: pi
      implicit none
      real(kind=rp), parameter :: half = 0.5_rp
c
      integer(kind=ip) :: n, i
      real(kind=rp) :: x(n),w(n)
      real(kind=rp) :: a,b,z,x1,x2
c
      x1 = (b-a)*half
      x2 = (b+a)*half
      do i = 1,n
        z = pi*(real(i,rp)-half)/real(n,rp)
        x(i) = x1*cos(z)+x2
        w(i) = (pi/real(n,rp))*x1*sin(z)
      end do
c
      end subroutine
c
c-----------------------------------------------------------------------
c
      subroutine gaucheb2 (x,w,n)
c
c-----This routine computes the points and weights of the Gauss-Chebyshev
c     quadrature of the first kind, neccessary to obtain the integral 
c     F(x) dx = sum(i=1,n) W_i F(x_i)
c
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
c
      end subroutine
c
c-----------------------------------------------------------------------
c
      subroutine pjt (x,w,n)
c
c-----This routine computes the points and weights of the Gauss-Chebyshev
c     quadrature of the first kind using the Perez-Jorda transformation.
c
      use mod_prec, only: rp, ip
      use mod_param, only: pi
      implicit none
      real(kind=rp), parameter :: twothird = 2.0_rp/3.0_rp
c
      integer(kind=ip) :: n, i
      real(kind=rp) :: x(n),w(n)
      real(kind=rp) :: z,st,ct,dn
c
      do i = 1,n
        dn = real(n+1_ip,rp)
        z = i*pi/dn
        st = sin(z)
        ct = cos(z)
        w(i) = 16.0_rp/3.0_rp/dn*st*st*st*st
        x(i) = 1.0_rp-real(2*i,rp)/dn+2.0_rp
     &                            /pi*(1.0_rp+twothird*st*st)*ct*st
      end do
c
      end subroutine
