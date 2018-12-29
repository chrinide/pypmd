      subroutine gradrho (xpoint,hi,iup,inf)
c
c.....Integration of a trajectory in the vector field of the
c     electron density.
c
c.....Input data:
c
c     xpoint() .... starting point of the trajectory
c     step ........ integration step. Enter negative value if default
c                value is wanted.
c     iup ......... +1 tells the routine to climb up in the gradient field
c                   (usually going to end in a nucleus).
c                   -1 tells the routine to go down in the field.
c
c.....Output data:
c
c     xpoint() .... end point of the trajectory 
c
      use mod_prec, only: rp, ip
      implicit none
c
      real(kind=rp), parameter :: epsg = 1d-10
      real(kind=rp), parameter :: eps = 1d-8
      real(kind=rp), parameter :: epsg1 = 1d-10
      real(kind=rp), parameter :: hminimal = 1d-40
      integer(kind=ip), parameter :: mstep = 500
c
      logical, intent(inout) :: inf
      integer(kind=ip), intent(in) :: iup
      real(kind=rp), intent(in) :: hi
      real(kind=rp), intent(inout) :: xpoint(3)
c
      real(kind=rp) :: gradmod, rho, grad(3)
      integer(kind=ip) :: ier, niter, npoints, i
      real(kind=rp) :: xtemp(3), grdt(3), h0, escalar, grdmodule
c
      interface
        subroutine pointr1 (p,rho,grad,gradmod)
          import rp
          real(kind=rp), intent(in) :: p(3)
          real(kind=rp), intent(out) :: grad(3)
          real(kind=rp), intent(out) :: rho
          real(kind=rp), intent(out) :: gradmod
        end subroutine
        subroutine pointshells (p,rho,grad,gradmod)
          import rp
          real(kind=rp), intent(in) :: p(3)
          real(kind=rp), intent(out) :: grad(3)
          real(kind=rp), intent(out) :: rho
          real(kind=rp), intent(out) :: gradmod
        end subroutine
      end interface
c
      h0 = hi
      npoints = 1_ip
      inf = .false.
*     call pointr1 (xpoint,rho,grad,gradmod)
      call pointshells (xpoint,rho,grad,gradmod)
      if (gradmod.lt.epsg .and. rho.lt.epsg1) then
        inf = .true.
        return
      end if
      grdt(1) = grad(1)/gradmod
      grdt(2) = grad(2)/gradmod
      grdt(3) = grad(3)/gradmod
      grdmodule = gradmod
c
      escalar = 1_ip
      niter = 1_ip
      do while (grdmodule.gt.eps .and. niter.lt.mstep)
        niter = niter + 1_ip
        ier = 1_ip
        do while (ier.ne.0)
          xtemp(1) = xpoint(1) + h0*iup*grdt(1)
          xtemp(2) = xpoint(2) + h0*iup*grdt(2)
          xtemp(3) = xpoint(3) + h0*iup*grdt(3)
*         call pointr1 (xtemp,rho,grad,gradmod)
          call pointshells (xtemp,rho,grad,gradmod)
          escalar = 0.0_rp
          do i = 1,3
            escalar = escalar + grdt(i)*(grad(i)/(gradmod+1d-80))
          end do
c
c.........Good direction
c
          if (escalar.lt.0.707_rp) then
c
c...........It should'nt happen that h0 goes to zero, except if there
c           are problems with the gradient, for instance the gradient
c           has discontinuities or large rounding errors. Anyway, we
c           introduce a safety check to avoid nasty infinite loops if
c           h0 = 0.
c
            if (h0.ge.hminimal) then
              h0 = h0/2.0_rp
              ier = 1_ip
            else
              ier = 0_ip
            end if
          else
            if (escalar.gt.0.9_rp) h0 = min(hi,h0*1.6_rp)
            ier = 0_ip
            ! Yes
            do i = 1,3
              xpoint(i) = xtemp(i)
              grdt(i) = grad(i)/gradmod
            enddo
            grdmodule = gradmod
          end if
        end do
      end do
c
      end subroutine
