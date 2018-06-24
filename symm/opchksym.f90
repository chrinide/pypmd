! check if xmat is a sym operator.
logical function opchksym(xmat)
  
  use iso_fortran_env, only: uout=>output_unit
  use mod_prec, only: ip, rp
  use mod_sym, only: nmol=>sdim, ax, ay, az, atgroup=>atzmol, &
                     nopsym, toldist, opsym, toleqvm
  implicit none

  real(kind=rp), dimension(3,3), intent(in) :: xmat

  integer :: i, j, k
  real(kind=rp) :: ximg, yimg, zimg, diff, xii, xij, xji 
  logical :: found

  ! First check: The matrix must be orthogonal:
  opchksym = .false.
  do i = 1,3
    xii = xmat(1,i)*xmat(1,i) + xmat(2,i)*xmat(2,i) + xmat(3,i)*xmat(3,i)
    if (abs(xii-1d0) .gt. TOLeqvm) return
    do j = i+1,3
      xij = xmat(1,i)*xmat(1,j) + xmat(2,i)*xmat(2,j) + xmat(3,i)*xmat(3,j)
      if (abs(xij) .gt. TOLeqvm) return
      xji = xmat(1,j)*xmat(1,i) + xmat(2,j)*xmat(2,i) + xmat(3,j)*xmat(3,i)
      if (abs(xji) .gt. TOLeqvm) return
    end do
  end do

  ! If the number of atoms is large, perhaps is better to check
  ! for repeated operations of symmetry before testing if it
  ! really transforms the molecule into itself:
  if (nmol .gt. nopsym+nopsym) then
    do i = 1,nopsym
      diff = 0.0_rp
      do j = 1,3
        do k = 1,3
          diff = diff + abs(xmat(j,k) - opsym(i,j,k))
        end do
      end do
      if(diff .le. TOLeqvm) return ! This is not a new operator
    end do
  end if

  ! Transform every atom and check for the symmetric image:
  do i = 1,nmol
    ximg = xmat(1,1)*ax(i) + xmat(1,2)*ay(i) + xmat(1,3)*az(i)
    yimg = xmat(2,1)*ax(i) + xmat(2,2)*ay(i) + xmat(2,3)*az(i)
    zimg = xmat(3,1)*ax(i) + xmat(3,2)*ay(i) + xmat(3,3)*az(i)
    found = .false.
    j = 1
    do while (.not.found .and. j.le.nmol)
      if (atgroup(i).eq.atgroup(j)) then
        found = abs(ximg-ax(j))+abs(yimg-ay(j))+abs(zimg-az(j)) .le. TOLdist
      end if
      j = j + 1
    end do
    if (.not.found) return
  end do
  opchksym = .true.
 
end function opchksym

