! put back the molecule to the original orientation.
! The routine uses the orientation matrices obtained by symorient()
! to reverse its action.
subroutine backsym (ax, ay, az, n)

  use iso_fortran_env, only: uout=>output_unit
  use mod_prec, only: rp, ip
  use mod_io, only: string
  use mod_param, only: verbose
  use mod_sym, only: or_mat, or_imat
  implicit none

  integer(kind=ip), intent(in)  :: n
  real(kind=rp), intent(inout) :: ax(n), ay(n), az(n)

  integer(kind=ip) :: i, j, iat
  real(kind=rp) :: xx, yy, zz, v(3,3), vinv(3,3)

  if (verbose) then
    write (uout,'(1x,a,1x,a)') string('# +++ SYMBACK: Return the molecule')// &
                               string('to the original orientation.')
  end if

  ! Get the transformation matrices from the symmetry database:
  do i = 1,3
    do j = 1,3
      v(i,j) = or_mat(i,j)
      vinv(i,j) = or_imat(i,j)
    end do
  end do

  ! Transform the local copy of the coordinates:
  do iat = 1,n
    xx = v(1,1)*ax(iat) + v(1,2)*ay(iat) + v(1,3)*az(iat)
    yy = v(2,1)*ax(iat) + v(2,2)*ay(iat) + v(2,3)*az(iat)
    zz = v(3,1)*ax(iat) + v(3,2)*ay(iat) + v(3,3)*az(iat)
    ax(iat) = xx
    ay(iat) = yy
    az(iat) = zz
  end do

  ! Transform the symmetry matrices too:
  call transformsym (v, vinv)

  ! Print out the transformed coordinates?
  if (verbose) then
    write (uout,'(1x,a)') string('# Transformed molecular coordiantes')
    write (uout,'(1x,i0,1x,3(1x,f12.8))') (i, ax(i), ay(i), az(i), i=1,n)
    write (uout,'(1x,a)') string('# Transformed molecular coordiantes')
  end if

end subroutine backsym
