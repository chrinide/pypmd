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
subroutine rotagrid (ct,st,cp,sp,npang)

  use mod_prec, only: rp, ip
  use mod_integ, only: angx, angy, angz
  implicit none
      
  real(kind=rp), parameter :: eps = 1.0d-07
  real(kind=rp), parameter :: zero = 0.0_rp
  real(kind=rp), parameter :: one = 1.0_rp

  integer(kind=ip), intent(in) :: npang
  real(kind=rp), intent(inout) :: ct(npang),st(npang),cp(npang),sp(npang)

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
