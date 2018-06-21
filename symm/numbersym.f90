! get the number of xmat in the symmetry operator list.
! A -1 value will be returned if xmat is not in the list.
integer(kind=ip) function numbersym (xmat)

  use mod_prec, only: rp, ip
  use mod_sym, only: nopsym, opsym, toleqvm
  implicit none
  real(kind=rp), dimension(3,3), intent(in) :: xmat

  integer(kind=ip) :: i, j, k
  real(kind=rp) :: diff

  numbersym = -1_ip
  do i = 1,nopsym
    diff = 0.0_rp
    do j = 1,3
      do k = 1,3
        diff = diff + abs(xmat(j,k)-opsym(i,j,k))
      end do
    end do
    if (diff .le. TOLeqvm) then
      numbersym = i
      return
    end if
  end do

end function numbersym
