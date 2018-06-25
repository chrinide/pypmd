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
module mod_linalg

  use mod_prec, only : i4=>ip, rp
  implicit none
  private

  integer(kind=i4), parameter :: maxrot = 500 ! max jacobi rotations

  ! overlad
  interface assert_shape
    module procedure dassert_shape_2d 
    module procedure dassert_shape_1d 
  end interface assert_shape

  !> public routines
  public :: jacobi

contains

  ! Resolution of a symmetric eigenvalue equation 
  subroutine jacobi(a, d, v, nrot)

    use mod_io, only : string, ferror, faterr
    implicit none

    integer(kind=i4), intent(out) :: nrot
    real(kind=rp), dimension(:), intent(out) :: d
    real(kind=rp), dimension(:,:), intent(inout) :: a
    real(kind=rp), dimension(:,:), intent(out) :: v

    integer(kind=i4) :: i, ip, iq, n
    real(kind=rp) :: c, g, h, s, sm, t, tau, theta, tresh
    real(kind=rp), dimension(size(d)) :: b, z

    ! get sizes
    n = size(a(1,:))
    
    ! check shape of matrix and vector
    call assert_shape(a, [n, n], "jacobi", "a")
    call assert_shape(d, [n, n], "jacobi", "d")
    call assert_shape(v, [n, n], "jacobi", "v")
  
    ! Init data
    nrot = 0_i4
    v(:,:) = 0.0_rp
    do ip = 1,n
      v(ip,ip) = 1.0_rp
      b(ip) = a(ip,ip)
    end do
    d(:) = b(:)
    z(:) = 0.0_rp

    ! Begin rotations
    do i = 1,maxrot
      sm = 0.0_rp
      do ip = 1,n-1 ! to parallelize collapse loops
        do iq = ip+1,n
          sm = sm + abs(a(ip,iq))
        end do
      end do
      if (sm.eq.0.0_rp) return
      tresh = merge(0.2_rp*sm/(n*n), 0.0_rp, i < 4)
      do ip = 1,n-1 ! to parallelize collapse loops 
        do iq = ip+1,n
          g = 100_rp*abs(a(ip,iq))
          if ((i.gt.4) .and. (abs(d(ip))+g.eq.abs(d(ip))) &
                       .and. (abs(d(iq))+g.eq.abs(d(iq)))) then
            a(ip,iq) = 0.0_rp
          else if (abs(a(ip,iq)).gt.tresh) then
            h = d(iq) - d(ip)
            if (abs(h)+g.eq.abs(h)) then
              t = a(ip,iq)/h
            else
              theta = 0.5_rp*h/a(ip,iq)
              t = 1.0_rp/(abs(theta)+sqrt(1.0_rp+(theta*theta)))
              if (theta.lt.0.0_rp) t = -t
            end if
            c = 1.0_rp/sqrt(1.0_rp+(t*t))
            s = t*c
            tau = s/(1.0_rp+c)
            h = t*a(ip,iq)
            z(ip) = z(ip) - h
            z(iq) = z(iq) + h
            d(ip) = d(ip) - h
            d(iq) = d(iq) + h
            a(ip,iq) = 0.0_rp
            call jrotate(a(1:ip-1,ip), a(1:ip-1, iq))
            call jrotate(a(ip, ip+1:iq-1), a(ip+1:iq-1, iq))
            call jrotate(a(ip, iq+1:n), a(iq,iq+1:n))
            call jrotate(v(:,ip), v(:,iq))
            nrot = nrot + 1_i4
          end if
        end do
      end do
      b(:) = b(:) + z(:)
      d(:) = b(:)
      z(:) = 0.0_rp
    end do

    call ferror('mod_linalg/jacobi', string(maxrot)//' iterations should never happen', faterr)

    contains

      subroutine jrotate(a1, a2)
    
        use mod_prec, only: rp
        implicit none
        real(kind=rp), dimension(:), intent(inout) :: a1, a2
        real(kind=rp), dimension(size(a1)) :: wk1

        wk1(:) = a1(:)
        a1(:) = a1(:) - s*(a2(:)+a1(:)*tau)
        a2(:) = a2(:) + s*(wk1(:)-a2(:)*tau)

      end subroutine jrotate

  end subroutine jacobi

  ! Diagonalize a real symmetric matrix, on ouput
  ! a(i,i) - eigenvalues
  ! x(i,j) - eigenvectors
  subroutine new_jacobi(a, x, abserr)

    use mod_io, only : string, ferror, faterr
    implicit none

    real(kind=rp), intent(inout), dimension(:,:) :: a
    real(kind=rp), intent(out), dimension(:,:) :: x
    real(kind=rp), intent(in) :: abserr

    integer(kind=i4) :: i, j, k, n, niter
    real(kind=rp) :: b2, bar, beta, coeff, c, s, cs, sc

    ! get sizes
    n = size(a(1,:))
    
    ! check shape of matrix and vector
    call assert_shape(a, [n, n], "new_jacobi", "a")
    call assert_shape(x, [n, n], "new_jacobi", "x")

    ! initialize
    x = 0.0_rp
    do i = 1,n
      x(i,i) = 1.0_rp
    end do

    ! find the sum of all off-diagonal elements (squared)
    b2 = 0.0_rp
    do i = 1,n
      do j = 1,n
        if (i.ne.j) b2 = b2 + (a(i,j)*a(i,j))
      end do
    end do

    ! alreadi diagonal ?
    if (b2 <= abserr) return

    ! average for off-diagonal elements /2
    bar = 0.5_rp*b2/real(n*n,rp)

    ! do the rotations, to parallel
    niter = 0_i4
    do while (b2.gt.abserr)
      do i = 1,n-1
        do j = i+1,n
          if (a(j,i)*a(j,i) <= bar) cycle ! do not touch small elements
          b2 = b2 - 2.0_rp*(a(j,i)*a(j,i))
          bar = 0.5_rp*b2/real(n*n,rp)
          ! calculate coefficient c and s for givens matrix
          beta = (a(j,j)-a(i,i))/(2.0_rp*a(j,i))
          coeff = 0.5_rp*beta/sqrt(1.0_rp+(beta*beta))
          s = sqrt(max(0.5_rp+coeff,0.0_rp))
          c = sqrt(max(0.5_rp-coeff,0.0_rp))
          ! recalculate rows i and j
          do k = 1,n
            cs = c*a(i,k)+s*a(j,k)
            sc = -s*a(i,k)+c*a(j,k)
            a(i,k) = cs
            a(j,k) = sc
          end do
          ! new matrix a_{k+1} from a_{k}, and eigenvectors 
          do k = 1,n
            cs = c*a(k,i)+s*a(k,j)
            sc = -s*a(k,i)+c*a(k,j)
            a(k,i) = cs
            a(k,j) = sc
            cs = c*x(k,i)+s*x(k,j)
            sc = -s*x(k,i)+c*x(k,j)
            x(k,i) = cs
            x(k,j) = sc
          end do
        end do
      end do
      if (niter.gt.maxrot) then
        call ferror('mod_linalg/new_jacobi', string(maxrot)//' iterations should never happen', faterr)
      end if
      niter = niter + 1_i4
    end do

  end subroutine

  !> make sure a given real matrix has a given shape
  subroutine dassert_shape_2d(a, shap, routine, matname)

    use mod_io, only: uout, faterr, ferror
    implicit none

    real(kind=rp), intent(in) :: a(:,:)
    integer(kind=i4), intent(in) :: shap(:)

    character(len=*) :: routine, matname

    if (any(shape(a) /= shap)) then
      write (uout,*) "# In routine "//routine//" matrix "//matname//" has illegal shape ", shape(a)
      write (uout,*) "# shape should be ", shap
      call ferror("mod_linalg/dassert", "aborting due to illegal operation", faterr)
    end if

  end subroutine dassert_shape_2d

  !> make sure a given real vector has a given shape
  subroutine dassert_shape_1d(a, shap, routine, matname)

    use mod_io, only: uout, faterr, ferror
    implicit none

    real(kind=rp), intent(in) :: a(:)
    integer(kind=i4), intent(in) :: shap(:)

    character(len=*) :: routine, matname

    if (any(shape(a) /= shap)) then
      write (uout,*) "# In routine "//routine//" vector "//matname//" has illegal shape ", shape(a)
      write (uout,*) "# shape should be ", shap
      call ferror("mod_linalg/dassert", "aborting due to illegal operation", faterr)
    end if

  end subroutine dassert_shape_1d

end module
