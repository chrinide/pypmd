module mod_fields

  use mod_prec, only: rp, ip
  implicit none
  private

  public :: density_grad, density_grad_shell
  public :: density_shell, mo_shell, basin

contains
  
  subroutine mo_shell (p,idx,mo)

    use mod_mole, only: coords_, ncent_
    use mod_basis, only: nlm, ityp_, oexp_, coeff_, &
                         ngroup_, nuexp_, nzexp_, rcutte_
    implicit none
 
    real(kind=rp), intent(in) :: p(3)
    real(kind=rp), intent(out) :: mo
    integer(kind=ip), intent(in) :: idx ! orbital to eval

    integer(kind=ip) :: it(3), i, ic, itip, j, jj, k, m, n
    real(kind=rp) :: aexp, cfj, dis2, f12, x2, f123, ori, x
    real(kind=rp) :: xcoor(3), fun(3), gun
 
    mo = 0.0_rp 
    gun = 0.0_rp

    do ic = 1,ncent_
      xcoor = p - coords_(:,ic)
      dis2 = xcoor(1)*xcoor(1) + xcoor(2)*xcoor(2) + xcoor(3)*xcoor(3)
      do m = 1,ngroup_(ic)
        k = nuexp_(ic,1,m)
        if (dis2.gt.rcutte_(ic,m)) goto 2
        ori = -oexp_(k)
        aexp = exp(ori*dis2)
        do jj = 1,nzexp_(ic,m)
          i = nuexp_(ic,jj,m)
          itip = ityp_(i)
          it(:) = nlm(itip,:)
          do j = 1,3
            n = it(j)
            x = xcoor(j)
            if (n.eq.0) then
              fun(j) = 1.0_rp
            else if (n.eq.1) then
              x2 = x*x
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
          f12 = fun(1)*fun(2)*aexp
          f123 = f12*fun(3)
          cfj = coeff_(idx,i)
          gun = gun + cfj*f123
        end do
 2      continue
      end do
    end do
 
    mo = gun
 
  end subroutine

  subroutine density_shell (p,rho)

    use mod_mole, only: coords_, ncent_
    use mod_basis, only: nlm, ityp_, occ_, nmo_, oexp_, coeff_, &
                         ngroup_, nuexp_, nzexp_, rcutte_
    implicit none
 
    real(kind=rp), intent(in) :: p(3)
    real(kind=rp), intent(out) :: rho

    integer(kind=ip) :: it(3), i, ic, itip, j, jj, k, m, n
    real(kind=rp) :: aexp, cfj, dis2, f12, x2, f123, ori, x
    real(kind=rp) :: xcoor(3), fun(3)
    real(kind=rp), dimension(nmo_) :: gun
 
    rho = 0.0_rp 
    gun = 0.0_rp

    do ic = 1,ncent_
      xcoor = p - coords_(:,ic)
      dis2 = xcoor(1)*xcoor(1) + xcoor(2)*xcoor(2) + xcoor(3)*xcoor(3)
      do m = 1,ngroup_(ic)
        k = nuexp_(ic,1,m)
        if (dis2.gt.rcutte_(ic,m)) goto 2
        ori = -oexp_(k)
        aexp = exp(ori*dis2)
        do jj = 1,nzexp_(ic,m)
          i = nuexp_(ic,jj,m)
          itip = ityp_(i)
          it(:) = nlm(itip,:)
          do j = 1,3
            n = it(j)
            x = xcoor(j)
            if (n.eq.0) then
              fun(j) = 1.0_rp
            else if (n.eq.1) then
              x2 = x*x
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
          f12 = fun(1)*fun(2)*aexp
          f123 = f12*fun(3)
          do j = 1,nmo_
            cfj = coeff_(j,i)
            gun(j) = gun(j) + cfj*f123
          end do
        end do
 2      continue
      end do
    end do
 
    rho = dot_product(occ_,gun*gun)
 
  end subroutine

  subroutine density_grad_shell (p,rho,grad,gradmod)

    use mod_mole, only: coords_, ncent_
    use mod_basis, only: nlm, ityp_, occ_, nmo_, oexp_, coeff_, &
                         ngroup_, nuexp_, nzexp_, rcutte_
    implicit none
 
    real(kind=rp), intent(in) :: p(3)
    real(kind=rp), intent(out) :: grad(3)
    real(kind=rp), intent(out) :: rho
    real(kind=rp), intent(out) :: gradmod

    integer(kind=ip) :: it(3), i, ic, itip, j, jj, k, m, n
    real(kind=rp) :: aexp, cfj, dis2, dp2, f12, fc, x2
    real(kind=rp) :: f123, fa, fb , ori, x
    real(kind=rp) :: xcoor(3), fun(3), fun1(3), dp2x, dp2x2
    real(kind=rp), dimension(nmo_,3) :: gun1
    real(kind=rp), dimension(nmo_) :: gun
 
    rho = 0.0_rp 
    grad = 0.0_rp
    gradmod = 0.0_rp
    gun = 0.0_rp
    gun1 = 0.0_rp

    do ic = 1,ncent_
      xcoor = p - coords_(:,ic)
      dis2 = xcoor(1)*xcoor(1) + xcoor(2)*xcoor(2) + xcoor(3)*xcoor(3)
      do m = 1,ngroup_(ic)
        k = nuexp_(ic,1,m)
        if (dis2.gt.rcutte_(ic,m)) goto 2
        ori = -oexp_(k)
        dp2 = ori + ori
        aexp = exp(ori*dis2)
        do jj = 1,nzexp_(ic,m)
          i = nuexp_(ic,jj,m)
          itip = ityp_(i)
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
          do j = 1,nmo_
            cfj = coeff_(j,i)
            gun(j) = gun(j) + cfj*f123
            gun1(j,1) = gun1(j,1) + cfj*fa
            gun1(j,2) = gun1(j,2) + cfj*fb
            gun1(j,3) = gun1(j,3) + cfj*fc
          end do
        end do
 2      continue
      end do
    end do
 
    rho = dot_product(occ_,gun*gun)
    grad(1) = dot_product(occ_,gun*gun1(:,1))
    grad(2) = dot_product(occ_,gun*gun1(:,2))
    grad(3) = dot_product(occ_,gun*gun1(:,3))
    grad(:) = grad(:) + grad(:)
    gradmod = norm2(grad)
 
  end subroutine
        
  ! Computes the spatial density and gradient at point p
  subroutine density_grad (p,rho,grad,gradmod)

    use mod_mole, only: coords_
    use mod_basis, only: nlm, nprims_, ityp_, occ_, nmo_, oexp_, coeff_, icen_
    implicit none
 
    real(kind=rp), intent(in) :: p(3)
    real(kind=rp), intent(out) :: grad(3)
    real(kind=rp), intent(out) :: rho
    real(kind=rp), intent(out) :: gradmod
  
    integer(kind=ip) :: it(3), i, ic, itip, j, n
    real(kind=rp) :: aexp, cfj, dis2, dp2, f12, fc, x2
    real(kind=rp) :: f123, fa, fb , ori, x
    real(kind=rp) :: xcoor(3), fun(3), fun1(3), dp2x, dp2x2
    real(kind=rp), dimension(nmo_,3) :: gun1
    real(kind=rp), dimension(nmo_) :: gun
 
    rho = 0.0_rp 
    grad = 0.0_rp
    gradmod = 0.0_rp
    gun = 0.0_rp
    gun1 = 0.0_rp

    do i = 1,nprims_
      ic = icen_(i)
      itip = ityp_(i)
      it(:) = nlm(itip,:)
      ori = -oexp_(i)
      dp2 = ori+ori
      xcoor = p - coords_(:,ic)
      dis2 = xcoor(1)*xcoor(1)+xcoor(2)*xcoor(2)+xcoor(3)*xcoor(3)
      aexp = exp(ori*dis2)
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
      do j = 1,nmo_
        cfj = coeff_(j,i)
        gun(j) = gun(j) + cfj*f123
        gun1(j,1) = gun1(j,1) + cfj*fa
        gun1(j,2) = gun1(j,2) + cfj*fb
        gun1(j,3) = gun1(j,3) + cfj*fc
      end do
    end do
  
    rho = dot_product(occ_,gun*gun)
    grad(1) = dot_product(occ_,gun*gun1(:,1))
    grad(2) = dot_product(occ_,gun*gun1(:,2))
    grad(3) = dot_product(occ_,gun*gun1(:,3))
    grad(:) = grad(:) + grad(:)
    gradmod = norm2(grad)
 
  end subroutine

  subroutine basin (p,rho,kin,lap)

    use mod_mole, only: coords_, ncent_
    use mod_basis, only: nlm, ityp_, occ_, nmo_, oexp_, coeff_, &
                         ngroup_, nuexp_, nzexp_, rcutte_
    implicit none
 
    real(kind=rp), intent(in) :: p(3)
    real(kind=rp), intent(out) :: rho
    real(kind=rp), intent(out) :: kin
    real(kind=rp), intent(out) :: lap

    integer(kind=ip) :: it(3), i, ic, itip, j, jj, k, m, n
    real(kind=rp) :: hess11, hess22, hess33
    real(kind=rp) :: aexp, cfj, dis2, dp2, f12, x2
    real(kind=rp) :: f123, f23 , ori, x, f13, fun2(3)
    real(kind=rp) :: xcoor(3), fun(3), fun1(3), dp2x, dp2x2
    real(kind=rp), dimension(nmo_,3) :: gun1
    real(kind=rp), dimension(nmo_,6) :: gun2
    real(kind=rp), dimension(nmo_) :: gun
 
    rho = 0.0_rp 
    kin = 0.0_rp
    lap = 0.0_rp
    gun = 0.0_rp
    gun1 = 0.0_rp
    gun2 = 0.0_rp

    do ic = 1,ncent_
      xcoor = p - coords_(:,ic)
      dis2 = xcoor(1)*xcoor(1) + xcoor(2)*xcoor(2) + xcoor(3)*xcoor(3)
      do m = 1,ngroup_(ic)
        k = nuexp_(ic,1,m)
        if (dis2.gt.rcutte_(ic,m)) goto 2
        ori = -oexp_(k)
        dp2 = ori + ori
        aexp = exp(ori*dis2)
        do jj = 1,nzexp_(ic,m)
          i = nuexp_(ic,jj,m)
          itip = ityp_(i)
          it(:) = nlm(itip,:)
          do j = 1,3
            n = it(j)
            x = xcoor(j)
            if (n.eq.0) then
              dp2x = dp2*x
              dp2x2 = dp2*x*x
              fun2(j) = dp2*(1.0_rp+dp2x2)
              fun1(j) = dp2x
              fun(j) = 1.0_rp
            else if (n.eq.1) then
              x2 = x*x
              dp2x2 = dp2*x2
              fun2(j) = dp2*x*(3.0_rp+dp2x2)
              fun1(j) = 1.0_rp+dp2x2
              fun(j) = x
            else if (n.eq.2) then
              x2 = x*x
              dp2x2 = dp2*x2
              fun2(j) = 2.0_rp+dp2x2*(5.0_rp+dp2x2)
              fun1(j) = x*(2.0_rp+dp2x2)
              fun(j) = x2
            else if (n.eq.3) then
              x2 = x*x
              dp2x2 = dp2*x2
              fun2(j) = x*(6.0_rp+dp2x2*(7.0_rp+dp2x2))
              fun1(j) = x2*(3.0_rp+dp2x2)
              fun(j) = x*x2
            else if (n.eq.4) then
              x2 = x*x
              dp2x2 = dp2*x2
              fun2(j) = x2*(x2*(dp2*(dp2x2+9.0_rp))+12.0_rp)
              fun1(j) = x2*x*(4.0_rp+dp2x2)
              fun(j) = x2*x2
            else if (n.eq.5) then
              x2 = x*x
              dp2x2 = dp2*x2
              fun2(j) = x2*x*(x2*(dp2*(dp2x2+11.0_rp))+20.0_rp) 
              fun1(j) = x2*x2*(5.0_rp+dp2x2)
              fun(j) = x2*x2*x
            end if
          end do
          f23 = fun(2)*fun(3)*aexp
          f13 = fun(1)*fun(3)*aexp
          f12 = fun(1)*fun(2)*aexp
          f123 = f12*fun(3)
          do j = 1,nmo_
            cfj = coeff_(j,i)
            gun(j) = gun(j) + cfj*f123
            gun1(j,1) = gun1(j,1) + cfj*(fun1(1)*f23)
            gun1(j,2) = gun1(j,2) + cfj*(fun1(2)*f13)
            gun1(j,3) = gun1(j,3) + cfj*(fun1(3)*f12)
            gun2(j,1) = gun2(j,1) + cfj*(fun2(1)*f23)
            gun2(j,3) = gun2(j,3) + cfj*(fun2(2)*f13)
            gun2(j,6) = gun2(j,6) + cfj*(fun2(3)*f12)
          end do
        end do
 2      continue
      end do
    end do
 
    rho = dot_product(occ_,gun*gun)
    kin = dot_product(occ_,(gun1(:,1)*gun1(:,1)))
    kin = kin + dot_product(occ_,(gun1(:,2)*gun1(:,2)))
    kin = kin + dot_product(occ_,(gun1(:,3)*gun1(:,3)))
    kin = kin*0.5_rp
    hess11 = dot_product(occ_,gun2(:,1)*gun(:))
    hess11 = hess11 + dot_product(occ_,gun1(:,1)*gun1(:,1)) 
    hess11 = hess11 + hess11
    hess22 = dot_product(occ_,gun2(:,3)*gun(:))
    hess22 = hess22 + dot_product(occ_,gun1(:,2)*gun1(:,2)) 
    hess22 = hess22 + hess22
    hess33 = dot_product(occ_,gun2(:,6)*gun(:))
    hess33 = hess33 + dot_product(occ_,gun1(:,3)*gun1(:,3)) 
    hess33 = hess33 + hess33
    lap = hess11 + hess22 + hess33

  end subroutine

end module mod_fields
