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
module mod_math

  use mod_prec, only: rp, ip
  implicit none
  public

contains
  
  ! check if this 3x3 matrix is singular
  logical function singular(xmat)

    implicit none
    real(kind=rp), dimension(3,3), intent(in) :: xmat
   
    real(kind=rp), parameter :: tolsng = 1e-5_rp
    real(kind=rp) :: det
    real(kind=rp), dimension(3,3) :: adj

    adj(1,1) = xmat(2,2)*xmat(3,3) - xmat(2,3)*xmat(3,2)
    adj(1,2) = xmat(2,3)*xmat(3,1) - xmat(2,1)*xmat(3,3)
    adj(1,3) = xmat(3,2)*xmat(2,1) - xmat(2,2)*xmat(3,1)

    det = xmat(1,1)*adj(1,1) + xmat(1,2)*adj(1,2) + xmat(1,3)*adj(1,3)

    if (abs(det) .le. tolsng) then
      singular = .true.
    else
      singular = .false.
    end if

    return

  end function singular

  ! Get the inverse of a 3x3 nonsingular matrix.
  ! The input matrix is not modified.
  subroutine inv (xmat, xinv, xdet, ierror)
  
    implicit none

    real(kind=rp), parameter :: tolsng = 1e-5_rp
    real(kind=rp), intent(in) :: xmat(3,3)
    real(kind=rp), intent(out) :: xinv(3,3), xdet
    integer(kind=ip), intent(out) :: ierror

    integer(kind=ip) :: i, j
    real(kind=rp) :: adj(3,3), det, ddet

    ierror = 0_ip
    adj(1,1) = xmat(2,2)*xmat(3,3) - xmat(2,3)*xmat(3,2)
    adj(1,2) = xmat(2,3)*xmat(3,1) - xmat(2,1)*xmat(3,3)
    adj(1,3) = xmat(3,2)*xmat(2,1) - xmat(2,2)*xmat(3,1)
    det = xmat(1,1)*adj(1,1) + xmat(1,2)*adj(1,2) + xmat(1,3)*adj(1,3)
    if (abs(det) .le. tolsng) then
      ierror = -1
      return
    end if
    adj(2,1) = xmat(1,3)*xmat(3,2) - xmat(1,2)*xmat(3,3)
    adj(2,2) = xmat(1,1)*xmat(3,3) - xmat(1,3)*xmat(3,1)
    adj(2,3) = xmat(1,2)*xmat(3,1) - xmat(1,1)*xmat(3,2)
    adj(3,1) = xmat(1,2)*xmat(2,3) - xmat(1,3)*xmat(2,2)
    adj(3,2) = xmat(1,3)*xmat(2,1) - xmat(1,1)*xmat(2,3)
    adj(3,3) = xmat(1,1)*xmat(2,2) - xmat(1,2)*xmat(2,1)
    ddet = 1.0_rp/det
    do i = 1,3
      do j = 1,3
        xinv(i,j) = adj(j,i)*ddet
      end do
    end do
    xdet = det
 
  end subroutine inv

  ! Fill in the xmat matrix.
  subroutine matfill (xmat, a11, a12, a13, a21, a22, a23, a31, a32, a33)
 
    implicit none

    real(kind=rp), dimension(3,3), intent(out) :: xmat
    real(kind=rp), intent(in) :: a11, a12, a13, a21
    real(kind=rp), intent(in) :: a22, a23, a31, a32, a33

    xmat(1,1) = a11
    xmat(1,2) = a12
    xmat(1,3) = a13
    xmat(2,1) = a21
    xmat(2,2) = a22
    xmat(2,3) = a23
    xmat(3,1) = a31
    xmat(3,2) = a32
    xmat(3,3) = a33

  end subroutine matfill

  ! perform the "C = A B" matrix multiplication between
  ! 3x3 matrices. The final matrix "C" can be equal to either "A" or
  ! "B", so the product will be stored on a secondary matrix before
  ! being moved to "C".
  subroutine matprod (a, b, c)

    implicit none

    real(kind=rp), dimension(3,3), intent(in) :: a, b
    real(kind=rp), dimension(3,3), intent(out) :: c

    c = matmul(a,b)
 
  end subroutine matprod

  ! get the (z,x,z) counterclockwise Euler angles that
  ! correspond to the rotation input "matrix".
  subroutine euler (matrix, eulera)

    use mod_param, only: pi
    implicit none

    real(kind=rp), dimension(3,3), intent(in) :: matrix
    real(kind=rp), dimension(3), intent(out) :: eulera

    if (abs(matrix(3,3)) .gt. 1.0_rp) then
      eulera(2) = acos(sign(1.0_rp,matrix(3,3)))
    else
      eulera(2) = acos(matrix(3,3))
    end if

    if (abs(abs(matrix(3,3))-1.0_rp) .le. 1d-6) then
      eulera(3) = 0.0_rp
      if (matrix(3,3) .gt. 0.0_rp) then
        eulera(1) = atan2(matrix(2,1), matrix(1,1))
      else
        eulera(1) = atan2(-matrix(2,1), matrix(1,1))
      end if
    else
      eulera(1) = atan2(matrix(3,1), matrix(3,2))
      eulera(3) = atan2(matrix(1,3), -matrix(2,3))
    end if

    if (eulera(1) .lt. 0.0_rp) eulera(1) = eulera(1) + pi + pi
    if (eulera(3) .lt. 0.0_rp) eulera(3) = eulera(3) + pi + pi

  end subroutine euler

  real(kind=rp) function dist (xyza, xyzb)

    implicit none
    real(kind=rp), intent(in), dimension(3) :: xyza, xyzb

    real(kind=rp) :: a, b, c

    a = (xyzb(1) - xyza(1))*(xyzb(1) - xyza(1)) 
    b = (xyzb(2) - xyza(2))*(xyzb(2) - xyza(2)) 
    c = (xyzb(3) - xyza(3))*(xyzb(3) - xyza(3)) 

    dist = sqrt(a+b+c)
    return

  end function dist

  ! calculate unit vector between to 3-d cartesian coordinates
  subroutine unitv (xyza, xyzb, xyzc)

    implicit none
    real(kind=rp), intent(in), dimension(3) :: xyza, xyzb
    real(kind=rp), intent(out), dimension(3) :: xyzc

    real(kind=rp) :: d
    real(kind=rp) :: a, b, c

    d = dist(xyza,xyzb)
    a = (xyzb(1) - xyza(1))*(xyzb(1) - xyza(1)) 
    b = (xyzb(2) - xyza(2))*(xyzb(2) - xyza(2)) 
    c = (xyzb(3) - xyza(3))*(xyzb(3) - xyza(3)) 
    xyzc(1) = a/d
    xyzc(2) = b/d
    xyzc(3) = c/d
  
  end subroutine unitv 

  ! calculate dot product between two unit vectors
  real(kind=rp) function unitdp (xyza, xyzb)

    implicit none
    real(kind=rp), intent(in), dimension(3) :: xyza, xyzb
    
    real(kind=rp) :: a, b, c

    unitdp = 0.0_rp
    a = xyzb(1)*xyza(1)
    b = xyzb(2)*xyza(2)
    c = xyzb(3)*xyza(3)
    unitdp = a + b + c
    unitdp = max(min(unitdp, 1.0_rp), -1.0_rp)

    return
  
  end function unitdp 

  ! calculate unit cross product between two unit vectors
  subroutine unitcp (xyza, xyzb, xyzc)

    implicit none
    real(kind=rp), intent(in), dimension(3) :: xyza, xyzb
    real(kind=rp), intent(out), dimension(3) :: xyzc

    real(kind=rp) :: cos_12, sin_12

    cos_12 = unitdp(xyza, xyzb)
    sin_12 = sqrt(1_rp - cos_12*cos_12)
    xyzc(1) = (xyza(2)*xyzb(3) - xyza(3)*xyzb(2))/sin_12
    xyzc(2) = (xyza(3)*xyzb(1) - xyza(1)*xyzb(3))/sin_12
    xyzc(3) = (xyza(1)*xyzb(2) - xyza(2)*xyzb(1))/sin_12

  end subroutine unitcp

  ! calculate angle between three 3-d cartesian coordinates
  real(kind=rp) function a123 (xyza, xyzb, xyzc)

    use mod_param, only: rad2deg
    implicit none
    real(kind=rp), intent(in), dimension(3) :: xyza, xyzb, xyzc

    real(kind=rp), dimension(3) :: u21, u23
    real(kind=rp) :: dp2123

    call unitv(xyzb, xyza, u21)
    call unitv(xyzb, xyzc, u23)
    dp2123 = unitdp(u21, u23)
    a123 = rad2deg*acos(dp2123)
    return 

  end function a123

end module mod_math
