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
module mod_fields

  use mod_prec, only: rp, ip
  implicit none
  private

  public :: pointr1, testpointr1

contains
        
  subroutine testpointr1 (xpoint)

    use iso_fortran_env, only: uout=>output_unit
    use mod_param, only: verbose
    use mod_io, only: string
    implicit none

    real(kind=rp), intent(in) :: xpoint(3)

    real(kind=rp) :: rho, grad(3), gradmod

    call pointr1 (xpoint,rho,grad,gradmod)
   
    if (verbose) then
      write (uout,'(1x,a)') string('# --- Follow values for test pointr1 ---')
      write (uout,'(1x,a,1x,3(1x,f0.8))') string('# Coordinates of point'), xpoint
      write (uout,'(1x,a,1x,f0.8)')  string('# Density value'), rho
      write (uout,'(1x,a,1x,3(1x,f0.8))') string('# Gradient value'), grad(:)
      write (uout,'(1x,a,1x,f0.8)') string('# Gradmod value'), gradmod
      write (uout,'(1x,a)') string('# --------------------------------------')
    end if
  
  end subroutine

  ! Computes the spatial density and gradient at point p
  subroutine pointr1 (p,rho,grad,gradmod)

    use mod_wfn, only: ncent, nlm, xyz, ityp, rcutte, occv, occupied, &
                       oexp, ngroup, nuexp, nzexp, coefnat, noccupied 
    implicit none
 
    real(kind=rp), intent(in) :: p(3)
    real(kind=rp), intent(out) :: grad(3)
    real(kind=rp), intent(out) :: rho
    real(kind=rp), intent(out) :: gradmod
  
    integer(kind=ip) :: it(3), i, ic, itip, j, jj, k, m, n, idx
    real(kind=rp) :: aexp, cfj, dis2, dp2, f12, fc, x2
    real(kind=rp) :: f123, fa, fb , ori, x
    real(kind=rp) :: xcoor(3), fun(3), fun1(3), dp2x, dp2x2
    real(kind=rp), dimension(noccupied,3) :: gun1
    real(kind=rp), dimension(noccupied) :: gun
 
    it = 0_ip
    fun = 0.0_rp
    fun1 = 0.0_rp
    rho = 0.0_rp 
    grad = 0.0_rp
    gradmod = 0.0_rp
    gun = 0.0_rp
    gun1 = 0.0_rp
  
    do ic = 1,ncent
      xcoor(:) = p(:) - xyz(ic,:)
      dis2 = xcoor(1)*xcoor(1) + xcoor(2)*xcoor(2) + xcoor(3)*xcoor(3)
      do m = 1,ngroup(ic)
        k = nuexp(ic,m,1)
        !if (dis2.gt.rcutte(ic,m)*rcutte(ic,m)) goto 2
        if (dis2.gt.rcutte(ic,m)) goto 2
        ori = -oexp(k)
        dp2 = ori + ori
        aexp = exp(ori*dis2)
        do jj = 1,nzexp(ic,m)
          i = nuexp(ic,m,jj)
          itip = ityp(i)
          it(:) = nlm(itip,:)
          do j = 1,3
            n = it(j)
            x = xcoor(j)
            if (n.eq.0) then
              dp2x = dp2*x
              fun1(j) = dp2x
              fun(j) = 1.0_rp
            else if (n.eq.1) then
              x2 = x*x
              dp2x2 = dp2*x2
              fun1(j) = 1.0_rp+dp2x2
              fun(j) = x
            else if (n.eq.2) then
              x2 = x*x
              dp2x2 = dp2*x2
              fun1(j) = x*(2.0_rp+dp2x2)
              fun(j) = x2
            else if (n.eq.3) then
              x2 = x*x
              dp2x2 = dp2*x2
              fun1(j) = x2*(3.0_rp+dp2x2)
              fun(j) = x*x2
            else if (n.eq.4) then
              x2 = x*x
              dp2x2 = dp2*x2
              fun1(j) = x2*x*(4.0_rp+dp2x2)
              fun(j) = x2*x2
            else if (n.eq.5) then
              x2 = x*x
              dp2x2 = dp2*x2
              fun1(j) = x2*x2*(5.0_rp+dp2x2)
              fun(j) = x2*x2*x
            end if
          end do
          f12 = fun(1)*fun(2)*aexp
          f123 = f12*fun(3)
          fa = fun1(1)*fun(2)*fun(3)*aexp
          fb = fun1(2)*fun(1)*fun(3)*aexp
          fc = fun1(3)*f12
          do j = 1,noccupied
            idx = occupied(j) 
            cfj = coefnat(idx,i)
            gun(idx) = gun(idx) + cfj*f123
            gun1(idx,1) = gun1(idx,1) + cfj*fa
            gun1(idx,2) = gun1(idx,2) + cfj*fb
            gun1(idx,3) = gun1(idx,3) + cfj*fc
          end do
        end do
2      continue
      end do
    end do
  
    rho = dot_product(occv,gun*gun)
    grad(1) = dot_product(occv,gun*gun1(:,1))
    grad(2) = dot_product(occv,gun*gun1(:,2))
    grad(3) = dot_product(occv,gun*gun1(:,3))
 
    grad(:) = grad(:) + grad(:)
    !gradmod = sqrt(grad(1)*grad(1)+grad(2)*grad(2)+grad(3)*grad(3))
    gradmod = norm2(grad)
    !grad = grad/(gradmod+1e-16) !TODO:change by small machine number
 
  end subroutine

end module mod_fields
