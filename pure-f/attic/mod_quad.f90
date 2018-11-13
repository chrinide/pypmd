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
module mod_quad

  use mod_prec, only: rp, ip
  implicit none
  private

  public :: rotagrid, rquad

contains

  ! Setup radial quadratures
  subroutine rquad (r0,rfar,rad,nr,rmesh,dvol,dvoln,rwei,mapr,iqudr)

    use mod_io, only: faterr, ferror
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
    if (iqudr.lt.1 .or. iqudr.gt.10) then
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
    else if (iqudr.eq.10) then 
      call tanh_sinh (xr,rwei,nr)
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
        end if
        dvol(n) = dvoln(n)*r
      end do
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
        end if
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

  ! This routine performs three successive rotations R_z, R_y, R_x
  ! about the Z, Y, and X axes, respectively, of NANG points placed on 
  ! an unit radius sphere, and defined by the polar coordinates 
  ! cos(theta)=ct, sin(theta)=st, cos(phi)=cp, and sin(phi)=sp.
  ! input
  ! angx, angy, angz = Rotation angle about the X, Y, and Z axes.
  ! ct,st,cp,sp = cos(theta),sin(theta),cos(phi),sin(phi). 
  ! nang = number of points on the unit radius sphere.
  ! maxang = maximum number of points.
  ! output
  ! ct,st,cp,sp = New polar coordinates after rotations.
  subroutine rotagrid (ct,st,cp,sp,npang,angx,angy,angz)

    implicit none
      
    real(kind=rp), parameter :: eps = 1.0e-07
    real(kind=rp), parameter :: zero = 0.0_rp
    real(kind=rp), parameter :: one = 1.0_rp
  
    integer(kind=ip), intent(in) :: npang
    real(kind=rp), intent(inout) :: ct(npang),st(npang),cp(npang),sp(npang)
    real(kind=rp), intent(in) :: angx, angy, angz
  
    real(kind=rp) :: xmat(3,3),ymat(3,3),zmat(3,3),rot(3,3)
    real(kind=rp) :: prod,x,y,z,xp,yp,zp,r,rxy,cangx,cangy,cangz,sangx,sangy,sangz
    integer(kind=ip) :: i,j,k

    ! Rotacion matrix about X axis.
    xmat = zero
    xmat(1,1) = one
    cangx = cos(angx)
    sangx = sin(angx)
    xmat(2,2) = +cangx
    xmat(3,2) = +sangx
    xmat(2,3) = -sangx
    xmat(3,3) = +cangx

    ! Rotacion matrix about Y axis.
    ymat = zero
    ymat(2,2) = one
    cangy = cos(angy)
    sangy = sin(angy)
    ymat(1,1) = +cangy
    ymat(3,1) = +sangy
    ymat(1,3) = -sangy
    ymat(3,3) = +cangy

    ! Rotacion matrix about Z axis.
    zmat = zero
    zmat(3,3) = one
    cangz = cos(angz)
    sangz = sin(angz)
    zmat(1,1) = +cangz
    zmat(2,1) = +sangz
    zmat(1,2) = -sangz
    zmat(2,2) = +cangz

    ! Full rotacion matrix. R = R_X * R_Y * R_Z
    do i = 1,3
      do j = 1,3
        prod = zero
        do k = 1,3
          prod = prod + ymat(i,k)*zmat(k,j)
        end do
        rot(i,j) = prod
      end do
    end do
    zmat(1:3,1:3) = rot(1:3,1:3)
    do i = 1,3
      do j = 1,3
        prod = zero
        do k = 1,3
          prod = prod + xmat(i,k)*zmat(k,j)
        end do
        rot(i,j) = prod
      end do
    end do
  
    ! Rotate angles.
    do i = 1,npang
      x = st(i)*cp(i)
      y = st(i)*sp(i)
      z = ct(i)
      xp = rot(1,1)*x + rot(1,2)*y + rot(1,3)*z
      yp = rot(2,1)*x + rot(2,2)*y + rot(2,3)*z
      zp = rot(3,1)*x + rot(3,2)*y + rot(3,3)*z
      ! Rebuild CT,ST,CP, and SP
      rxy = xp*xp + yp*yp
      r = sqrt(rxy+zp*zp)
      if (rxy.lt.eps) then
        if (zp.ge.zero) then
          ct(i) = +one
        else
          ct(i) = -one
        endif
        st(i) = zero
        sp(i) = zero
        cp(i) = one
      else
        rxy = sqrt(rxy)
        ct(i) = zp/r
        st(i) = sqrt((one-ct(i))*(one+ct(i)))
        cp(i) = xp/rxy
        sp(i) = yp/rxy
      end if
    end do
            
  end subroutine

end module mod_quad
