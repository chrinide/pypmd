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
      subroutine gradrho (xpoint,hi,iup,inf)
!
!.....Integration of a trajectory in the vector field of the
!     electron density.
!
!.....Input data:
!
!     xpoint() .... starting point of the trajectory
!     step ........ integration step. Enter negative value if default
!                value is wanted.
!     iup ......... +1 tells the routine to climb up in the gradient field
!                   (usually going to end in a nucleus).
!                   -1 tells the routine to go down in the field.
!
!.....Output data:
!
!     xpoint() .... end point of the trajectory 
!
      use mod_prec, only: rp, ip
      implicit none
!
      real(kind=rp), parameter :: epsg = 1d-10
      real(kind=rp), parameter :: eps = 1d-8
      real(kind=rp), parameter :: epsg1 = 1d-10
      real(kind=rp), parameter :: hminimal = 1d-40
      integer(kind=ip), parameter :: mstep = 500
!
      logical, intent(inout) :: inf
      integer(kind=ip), intent(in) :: iup
      real(kind=rp), intent(in) :: hi
      real(kind=rp), intent(inout) :: xpoint(3)
!
      real(kind=rp) :: gradmod, rho, grad(3)
      integer(kind=ip) :: ier, niter, npoints, i
      real(kind=rp) :: xtemp(3), grdt(3), h0, escalar, grdmodule
!
      interface
        subroutine pointr1 (p,rho,grad,gradmod)
          import rp
          real(kind=rp), intent(in) :: p(3)
          real(kind=rp), intent(out) :: grad(3)
          real(kind=rp), intent(out) :: rho
          real(kind=rp), intent(out) :: gradmod
        end subroutine
      end interface
!
      h0 = hi
      npoints = 1_ip
      inf = .false.
      call pointr1 (xpoint,rho,grad,gradmod)
      if (gradmod.lt.epsg .and. rho.lt.epsg1) then
        inf = .true.
        return
      end if
      grdt(1) = grad(1)/gradmod
      grdt(2) = grad(2)/gradmod
      grdt(3) = grad(3)/gradmod
      grdmodule = gradmod
!
      escalar = 1_ip
      niter = 1_ip
      do while (grdmodule.gt.eps .and. niter.lt.mstep)
        niter = niter + 1_ip
        ier = 1_ip
        do while (ier.ne.0)
          xtemp(1) = xpoint(1) + h0*iup*grdt(1)
          xtemp(2) = xpoint(2) + h0*iup*grdt(2)
          xtemp(3) = xpoint(3) + h0*iup*grdt(3)
          call pointr1 (xtemp,rho,grad,gradmod)
          escalar = 0.0_rp
          do i = 1,3
            escalar = escalar + grdt(i)*(grad(i)/(gradmod+1d-80))
          end do
!
!.........Good direction
!
          if (escalar.lt.0.707_rp) then
!
!...........It should'nt happen that h0 goes to zero, except if there
!           are problems with the gradient, for instance the gradient
!           has discontinuities or large rounding errors. Anyway, we
!           introduce a safety check to avoid nasty infinite loops if
!           h0 = 0.
!
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
!
      end subroutine
