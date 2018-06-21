! obtain the eigenvalues and eigenvectors of the
! orthogonal matrix xmat(,) by calling to the dgeev() routine in
! the lapack package. If xdet is +/-1, the routine uses this
! value as the determinant of xmat(,). Otherwise the determinant
! is obtained.
!
! The eigenvalues are sorted with xval(3) being the equivalent
! to the matrix determinant and xvect(:,3) being the normal vector
! corresponding to the rotation axis.
! The two other eigenvalues are meaningless (the imaginary part is
! not returned) but their eigenvectors define the plane normal to
! the rotation axis.
! The three eigenvectors are real and orthonormal.
subroutine eigensym (xmat, xval, xvect, xdet, ierr)

  use mod_prec, only: rp, ip
  use mod_sym, only: toleigen, erreigen
  implicit none

  integer(kind=ip), intent(out) :: ierr
  real(kind=rp), intent(in), dimension(3,3) :: xmat
  real(kind=rp), intent(out), dimension(3,3) :: xvect
  real(kind=rp), intent(out), dimension(3) :: xval
  real(kind=rp), intent(out) :: xdet
  
  integer(kind=ip), parameter :: lwork = 102
  real(kind=rp), dimension(lwork) ::  work

  integer(kind=ip) :: i, ii, ij, ik
  real(kind=rp) :: dif, difmax, xnorm
  real(kind=rp), dimension(3) :: wr, wi
  real(kind=rp), dimension(3,3) :: ar, vr, vl

  ! Get the determinant if it is needed:
  if (abs(abs(xdet)-1.0_rp) .ge. TOLeigen) then
    xdet = xmat(1,1)*(xmat(2,2)*xmat(3,3) - xmat(2,3)*xmat(3,2)) &
         + xmat(1,2)*(xmat(2,3)*xmat(3,1) - xmat(2,1)*xmat(3,3)) &
         + xmat(1,3)*(xmat(2,1)*xmat(3,2) - xmat(2,2)*xmat(3,1))
  end if
  ar = xmat
  
  call dgeev ("N", "V", 3, ar, 3, wr, wi, vl, 3, vr, 3, work, lwork, ierr)

  ! Get the eigenvalue identical to the determinant and return it as
  ! the third one:
  ik = 0
  difmax = 1d40
  do i = 1,3
    dif = abs(wr(i)-xdet) + abs(wi(i))
    if (dif .lt. difmax) then
      ik = i
      difmax = dif
    end if
  end do
  ii = mod(ik,3) + 1
  ij = mod(ii,3) + 1
  ERReigen = max(ERReigen, difmax)
  xval(3) = wr(ik)
  xvect(1,3) = vr(1,ik)
  xvect(2,3) = vr(2,ik)
  xvect(3,3) = vr(3,ik)
  xval(1) = wr(ii)
  xvect(1,1) = vr(1,ii)
  xvect(2,1) = vr(2,ii)
  xvect(3,1) = vr(3,ii)
  xval(2) = wr(ij)

  ! Get explicitely the v2 = v3 x v1 eigenvector:
  xvect(1,2) = xvect(2,3)*xvect(3,1) - xvect(3,3)*xvect(2,1)
  xvect(2,2) = xvect(3,3)*xvect(1,1) - xvect(1,3)*xvect(3,1)
  xvect(3,2) = xvect(1,3)*xvect(2,1) - xvect(2,3)*xvect(1,1)

  ! Enforce the normalization:
  do i = 1,3
    xnorm = xvect(1,i)**2 + xvect(2,i)**2 + xvect(3,i)**2
    xnorm = 1.0_rp/sqrt(xnorm)
    xvect(1,i) = xvect(1,i)*xnorm
    xvect(2,i) = xvect(2,i)*xnorm
    xvect(3,i) = xvect(3,i)*xnorm
  end do
 
end subroutine eigensym
