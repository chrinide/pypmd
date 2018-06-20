      subroutine polich (x, n, np, c)
!
!-----Let x(i), i=1,n, the n eigenvalues of a real symmetric matrix A
!     of logical size n x n and physical size np x np; i.e. the n solu-
!     tions of the secular equation of A
!           
!     det | x I - A | = SUM (k=0,N) c(k) x^k = 0
!             -   -
!     are x(1), x(2),..., x(n). c(k) are the cofficients of the charac-
!     teristic polynomial and are the objective of this routine.
!
      use mod_prec, only: rp, ip
      use mod_io, only: ferror, faterr
      implicit none
 
      integer(kind=ip), intent(in) :: n, np
      real(kind=rp) :: x(np), c(0:np)

      integer(kind=ip) :: i, k
 
      if (n.gt.0) then
        c(0) = -x(1)
        c(1) = 1.0_rp
        do i = 2,n
          c(i) = c(i-1)
          do k = i-1,1,-1
            c(k) = c(k-1)-x(i)*c(k)
          end do
          c(0) = -x(i)*c(0)
        end do
      else
         call ferror('polich','improper dimension', faterr)
      end if

      end subroutine
