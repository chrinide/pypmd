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
! Given coeff return orbitals 
!
subroutine orbfun (p,gun,coef)

  use mod_prec, only: rp, ip
  use mod_wfn, only: ncent, nlm, xyz, ityp, rcutte, nprims, &
                     oexp, ngroup, nuexp, nzexp, nmo
  implicit none
 
  real(kind=rp), intent(in) :: p(3)
  real(kind=rp), dimension(nmo), intent(out) :: gun
  real(kind=rp), dimension(nmo,nprims), intent(in) :: coef
 
  integer(kind=ip) :: it(3), i, ic, itip, j, jj, k, m, n
  real(kind=rp) :: aexp, cfj, dis2, dp2, x2, f123, ori, x
  real(kind=rp) :: xcoor(3), fun(3)
 
  it = 0_ip
  fun = 0.0_rp
  gun = 0.0_rp
 
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
            fun(j) = 1.0_rp
          else if (n.eq.1) then
            fun(j) = x
          else if (n.eq.2) then
            x2 = x*x
            fun(j) = x2
          else if (n.eq.3) then
            x2 = x*x
            fun(j) = x*x2
          else if (n.eq.4) then
            x2 = x*x
            fun(j) = x2*x2
          else if (n.eq.5) then
            x2 = x*x
            fun(j) = x2*x2*x
          end if
        end do
        f123 = fun(1)*fun(2)*fun(3)*aexp
        do j = 1,nmo
          cfj = coef(j,i)
          gun(j) = gun(j) + cfj*f123
        end do
      end do
 2    continue
    end do
  end do
 
end subroutine
