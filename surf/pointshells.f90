   subroutine pointshells (p,rho,grad,gradmod)
!
!..computes the spatial density and gradient at point p().
!..When this routine is called it is always true that we are computing
!  the electron density and its gradient at a point p() associated
!  the nucleus 'INUC'
! 
   use mod_prec, only: rp, ip
   use mod_surf, only: inuc
   use mod_wfn, only: noccupied, nlm, xyz, oexp, occv, ityp, rcutte, &
                      nzexp, coef, atcenter, ishell, nshell, nuexp,  &
                      occupied 
   implicit none
!
   real(kind=rp), intent(in) :: p(3)
   real(kind=rp), intent(out) :: rho, grad(3), gradmod
!
   integer(kind=ip) :: it(3), i, ic, itip, j, jj, is, k, m, n, idx
   real(kind=rp) :: aexp, cfj, dis2, dp2, dp2x, dp2x2
   real(kind=rp) :: f12, f123, fa, fb, fc
   real(kind=rp) :: xcoor(3), fun(3), fun1(3), ori, x, x2
   real(kind=rp), dimension(noccupied,3) :: gun1
   real(kind=rp), dimension(noccupied) :: gun
!
   it = 0_ip
   fun = 0.0_rp
   fun1 = 0.0_rp
   rho = 0.0_rp 
   grad = 0.0_rp
   gradmod = 0.0_rp
   gun = 0.0_rp
   gun1 = 0.0_rp
!
!  Run over centers
!
   do is = 1,nshell(inuc)
     ic = atcenter(inuc,is)
     xcoor(:) = p(:) - xyz(ic,:)
     dis2 = xcoor(1)*xcoor(1) + xcoor(2)*xcoor(2) + xcoor(3)*xcoor(3)
     m = ishell(inuc,is)
     k = nuexp(ic,m,1)
!
!....skip to compute this primitive if distance is too big.
!
!    if (dis2.gt.rcutte(ic,m)*rcutte(ic,m)) goto 2
     if (dis2.gt.rcutte(ic,m)) goto 2
!
     ori = -oexp(k)
     dp2 = ori + ori
!
!....All primitives in a shell share the same exponent.
!
     aexp = exp(ori*dis2)
!
!....Loop over the different primitives in this shell.
!
     do jj = 1,nzexp(ic,m)
!
!......"i" is the original index of the primitive in the WFN.
!
       i = nuexp(ic,m,jj)
       itip = ityp(i)
!
!......Integer coefficients.
!
       it(:) = nlm(itip,:)
       do j = 1,3
         n = it(j)
         x = xcoor(j)
         if (n.eq.0) then
           dp2x = dp2*x
           fun1(j) = dp2x
           fun(j) = 1d0
         else if (n.eq.1) then
           x2 = x*x
           dp2x2 = dp2*x2
           fun1(j) = 1d0+dp2x2
           fun(j) = x
         else if (n.eq.2) then
           x2 = x*x
           dp2x2 = dp2*x2
           fun1(j) = x*(2d0+dp2x2)
           fun(j) = x2
         else if (n.eq.3) then
           x2 = x*x
           dp2x2 = dp2*x2
           fun1(j) = x2*(3d0+dp2x2)
           fun(j) = x*x2
         else if (n.eq.4) then
           x2 = x*x
           dp2x2 = dp2*x2
           fun1(j) = x2*x*(4d0+dp2x2)
           fun(j) = x2*x2
         else if (n.eq.5) then
           x2 = x*x
           dp2x2 = dp2*x2
           fun1(j) = x2*x2*(5d0+dp2x2)
           fun(j) = x2*x2*x
         end if
       end do
       f123 = fun(1)*fun(2)*fun(3)*aexp
       fa = fun1(1)*fun(2)*fun(3)*aexp
       fb = fun1(2)*fun(1)*fun(3)*aexp
       fc = fun1(3)*fun(1)*fun(2)*aexp 
!
!......run over orbitals.
!
       do j = 1,noccupied
         idx = occupied(j)
         cfj = coef(idx,i)
         gun(idx) = gun(idx) + cfj*f123
         gun1(idx,1)= gun1(idx,1) + cfj*fa
         gun1(idx,2)= gun1(idx,2) + cfj*fb
         gun1(idx,3)= gun1(idx,3) + cfj*fc
       end do
     end do
 2 continue
   end do
!
!..run again over orbitals.
!
   rho = dot_product(occv,gun(:)*gun(:))
   grad(1) = dot_product(occv,gun(:)*gun1(:,1))
   grad(2) = dot_product(occv,gun(:)*gun1(:,2))
   grad(3) = dot_product(occv,gun(:)*gun1(:,3))
!
   grad(:) = grad(:) + grad(:)
   gradmod = sqrt(grad(1)*grad(1)+grad(2)*grad(2)+grad(3)*grad(3))
!
   end subroutine
