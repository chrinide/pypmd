! Setup radial quadratures
subroutine rquad (r0,rfar,rad,nr,rmesh,dvol,dvoln,rwei,mapr,iqudr)

  use mod_prec, only: rp, ip
  use mod_io, only: faterr, ferror
  use mod_integ, only: hlevel, nlevel
  implicit none
  real(kind=rp), parameter :: epsden = 1d-7
 
  integer(kind=ip), intent(in) :: nr, mapr, iqudr
  real(kind=rp), intent(in) :: r0, rfar, rad 
  real(kind=rp), intent(out) :: rmesh(nr), dvol(nr), dvoln(nr), rwei(nr)
 
  integer(kind=ip) :: n
  real(kind=rp) :: u, aa, bb, r, rfarc, den, eta
  real(kind=rp), dimension(nr) :: xr
 
  ! Make checks and init
  if (rfar-r0.le.0.001_rp) then
    call ferror ('rquad', 'rmax < rmin ???', faterr)
  end if
  if (iqudr.lt.1 .or. iqudr.gt.9) then
    call ferror ('rquad', 'not allowed radial quadrature', faterr)
  end if
  if (mapr.lt.0 .or. mapr.gt.3) then
    call ferror ('rquad', 'not allowed radial mapping', faterr)
  end if
  rmesh = 0.0_rp
  dvoln = 0.0_rp
  dvol = 0.0_rp
  rwei = 0.0_rp
  rfarc = 0.0_rp
  eta = 0.0_rp

  ! Determine eta parameter in case of radial mapping
  rfarc = rfar - r0
  if (mapr.eq.1) then
    eta = 2.0_rp*rad/rfarc
  else if (mapr.eq.2) then
    eta = 2.0_rp*exp(-rfarc/rad)/(1.0_rp-exp(-rfarc/rad))
  end if

  ! Determine abscissas and weights of the quadrature.
  if (iqudr.eq.1) then 
    call legendre_ss (xr,rwei,nr)
  else if (iqudr.eq.2) then 
    call chebyshev1 (xr,rwei,nr)
  else if (iqudr.eq.3) then 
    call chebyshev2 (xr,rwei,nr)
  else if (iqudr.eq.4) then  
    call fejer1 (xr,rwei,nr)
  else if (iqudr.eq.5) then  
    call fejer2 (xr,rwei,nr)
  else if (iqudr.eq.6) then  
    call lobatto (xr,rwei,nr)
  else if (iqudr.eq.7) then  
    call radau (xr,rwei,nr)
  else if (iqudr.eq.8) then  
    call pjt (xr,rwei,nr)
  else if (iqudr.eq.9) then 
    call clenshaw_curtis (xr,rwei,nr)
  else
    call ferror ('rquad', 'not allowed radial quadrature', faterr)
  end if

  ! Determine abscissas and volume elements.
  if (mapr.eq.0) then
    ! for finite range (a..b) the transformation is y = (b-a)*x/2+(b+a)/2
    ! x = (b-a)*0.5_rp*x+(b+a)*0.5_rp
    ! w = w*(b-a)*0.5_rp
    do n = 1,nr
      aa = (rfar-r0)/2.0_rp
      bb = (rfar+r0)/2.0_rp
      u = xr(n)
      r = aa*u+bb
      rmesh(n) = r
      dvoln(n) = r*aa
      dvol(n) = dvoln(n)*r
    end do  
  else if (mapr.eq.1) then
    do n = 1,nr
      u = xr(n)
      den = (1-u+eta)
      r = rad*(1+u)/den + r0
      rmesh(n) = r
      if (abs(den).gt.epsden) then !this
        dvoln(n) = rad*(2.0_rp+eta)/den/den*r
      else
        dvoln(n) = 0.0_rp
      endif
      dvol(n) = dvoln(n)*r
    enddo
  else if (mapr.eq.2) then
    do n = 1,nr
      u = xr(n)
      den = (1.0_rp-u+eta)
      r = rad*log((2.0_rp+eta)/den) + r0
      rmesh(n) = r
      if (abs(den).gt.epsden) then !this
        dvoln(n) = r*rad/den
      else
        dvoln(n) = 0.0_rp
      endif
      dvol(n) = dvoln(n)*r
    end do
  else if (mapr.eq.3) then
    do n = 1,nr
      u = xr(n)
      r = (1.0_rp+u)*rfarc/2.0_rp + r0
      rmesh(n) = r
      dvoln(n) = r*rfarc/2.0_rp
      dvol(n) = dvoln(n)*r
    end do
  end if
 
  end subroutine
