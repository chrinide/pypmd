      subroutine pointcrystal(p)
c 
      use mod_crystal
      use space_for_wfnbasis
      use space_for_wfncoef
      implicit none
      include 'param.inc'    
      include 'wfn.inc'    
      include 'point.inc'    
c
c----->>> To Kpoints: resize f123, fa, fb, fc, gun, gun1, factor to store
c----->>> the value of each point, then take the mediam for equally weight
c----->>> kpoints
c
      real(kind=8), intent(in) :: p(3)
c
      integer :: it(3), i, ic, itip, j, n, t
      real(kind=8) :: xcoor(3), aexp, cfj, dis2, fun(3), dp2
      real(kind=8) :: fac, x, ori, x2, fun1(3), dp2x, dp2x2
      real(kind=8) :: f123, gun(nmo), factor, gun1(nmo,3), fa, fb, fc
c
      it=0
      gun=0d0
      gun1=0d0
      xcoor=0d0
      fun=0d0
      fun1=0d0
      rho=0d0
      grad=0d0
      gradmod=0d0
c
      do i=1,nprims
        ic=icen(i)
        itip=ityp(i)
        it(1)=nlm(itip,1)
        it(2)=nlm(itip,2)
        it(3)=nlm(itip,3)
        ori=-oexp(i)
        dp2=ori+ori
c
c------>>Loop over T vectors to get SALCS
        fa=0d0
        fb=0d0
        fc=0d0
        f123=0d0
        do t=1,ntvectors
          xcoor(1)=p(1)-(xyz(ic,1)+tvectors(t,1))
          xcoor(2)=p(2)-(xyz(ic,2)+tvectors(t,2))
          xcoor(3)=p(3)-(xyz(ic,3)+tvectors(t,3))
          dis2=xcoor(1)*xcoor(1)+xcoor(2)*xcoor(2)+xcoor(3)*xcoor(3)
          aexp=exp(ori*dis2)
          do j=1,3
            n=it(j)
            x=xcoor(j)
            if (n.eq.0) then
              dp2x=dp2*x
              fun1(j)=dp2x
              fun(j)=1d0
            else if (n.eq.1) then
              x2=x*x
              dp2x2=dp2*x2
              fun1(j)=1d0+dp2x2
              fun(j)=x
            else if (n.eq.2) then
              x2=x*x
              dp2x2=dp2*x2
              fun1(j)=x*(2d0+dp2x2)
              fun(j)=x2
            else if (n.eq.3) then
              x2=x*x
              dp2x2=dp2*x2
              fun1(j)=x2*(3d0+dp2x2)
              fun(j)=x*x2
            else if (n.eq.4) then
              x2=x*x
              dp2x2=dp2*x2
              fun1(j)=x2*x*(4d0+dp2x2)
              fun(j)=x2*x2
            else if (n.eq.5) then
              x2=x*x
              dp2x2=dp2*x2
              fun1(j)=x2*x2*(5d0+dp2x2)
              fun(j)=x2*x2*x
            end if
          end do
c
c------>>Loop over K points and store for each k-point
!         factor=exp(dot_product(xyz(ic,:)+tvectors(t,:),
!    &                           imag*kpoints(1,:)))
!         GAMMA POINT only real
          factor=1.0
          f123=f123+fun(1)*fun(2)*fun(3)*aexp*factor 
          fa=fa+fun1(1)*fun(2)*fun(3)*aexp*factor
          fb=fb+fun1(2)*fun(1)*fun(3)*aexp*factor
          fc=fc+fun1(3)*fun(1)*fun(2)*aexp*factor 
c------>>Loop over K points and store for each k-point
c
        enddo
c------>>Loop over T vectors to get SALCS
c
c------>>Loop over K points and store for each k-point
        do j=1,nmo
          cfj=coef(j,i)
          gun(j)=gun(j)+cfj*f123
          gun1(j,1)=gun1(j,1)+cfj*fa
          gun1(j,2)=gun1(j,2)+cfj*fb
          gun1(j,3)=gun1(j,3)+cfj*fc
        enddo
c------>>Loop over K points and store for each k-point
c
      enddo
c
c------>>Loop over K points and store for each k-point
!     do i=1,nmo
!       fac=occ(i)
!       rho=rho+fac*conjg(gun(i))*gun(i) 
!       grad(1)=grad(1)+fac*gun(i)*conjg(gun1(i,1))
!       grad(2)=grad(2)+fac*gun(i)*conjg(gun1(i,2))
!       grad(3)=grad(3)+fac*gun(i)*conjg(gun1(i,3))
!     enddo
      rho=dot_product(occ,gun*gun)
      grad(1)=dot_product(occ,gun*gun1(:,1))
      grad(2)=dot_product(occ,gun*gun1(:,2))
      grad(3)=dot_product(occ,gun*gun1(:,3))
c
      grad(1:3)=grad(1:3)+grad(1:3)
!     gradmod=sqrt(grad(1)*grad(1)+grad(2)*grad(2)+grad(3)*grad(3))
      gradmod=norm2(grad)
c------>>Loop over K points and store for each k-point
c
      end subroutine
