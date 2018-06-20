subroutine cmcq ()
    
  use mod_prec, only: rp, ip
  use mod_io, only: faterr, ferror
  use mod_linalg, only: jacobi
  use mod_datatm, only: wgatm
  use mod_wfn, only: ncent, charge, xyz
  use mod_geom, only: cm, cq, moi, emoi, evmoi

  implicit none
  real(kind=rp) :: wt, zt, twt, tzt, wi, xx, yy, zz
  integer(kind=ip) :: i, nrot

  cm = 0.0_rp
  cq = 0.0_rp
  wt = 0.0_rp
  zt = 0.0_rp
  moi = 0.0_rp

  ! compute the center of mass and charge
  do i = 1,ncent
    twt = wgatm(int(charge(i)))
    wt = wt + twt
    cm(1) = cm(1) + twt*xyz(i,1)
    cm(2) = cm(2) + twt*xyz(i,2)
    cm(3) = cm(3) + twt*xyz(i,3)
    tzt = charge(i)
    zt = zt + tzt
    cq(1) = cq(1) + tzt*xyz(i,1)
    cq(2) = cq(2) + tzt*xyz(i,2)
    cq(3) = cq(3) + tzt*xyz(i,3)
  end do
  cm = cm / wt
  cq = cq / zt

  ! compute the inertia matrix and diagonalize it
  do i = 1,ncent
    wi = wgatm(int(charge(i)))
    xx = xyz(i,1) - cm(1)
    yy = xyz(i,2) - cm(2)
    zz = xyz(i,3) - cm(3)
    moi(1,1) = moi(1,1) + wi*(yy*yy + zz*zz)
    moi(2,2) = moi(2,2) + wi*(zz*zz + xx*xx)
    moi(3,3) = moi(3,3) + wi*(xx*xx + yy*yy)
    moi(1,2) = moi(1,2) - wi*xx*yy
    moi(2,3) = moi(2,3) - wi*yy*zz
    moi(3,1) = moi(3,1) - wi*zz*xx
  end do
  moi(2,1) = moi(1,2)
  moi(3,2) = moi(2,3)
  moi(1,3) = moi(3,1)

  call jacobi (moi,emoi,evmoi,nrot)
  if (nrot .lt. 0) call ferror ('cmcq','fail diag inertia', faterr)

end subroutine
