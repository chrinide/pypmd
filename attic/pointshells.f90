  subroutine pointshells (point,rho,grad,gradmod,norm)

    use mod_prec, only: rp, ip
    use mod_wfn, only: ncent, rcutte, xyz, nshell, atcenter, &
                       ishell, nuexp, nzexp, nuexp, ityp, nlm, &
                       coefnat, occ, nmo, oexp 
    implicit none
  
    real(kind=rp), dimension(3), intent(in) :: point
    real(kind=rp), intent(out) :: rho, grad(3), gradmod
    logical, intent(in), optional :: norm

    real(kind=rp) :: xcoor(3), fun(3), fun1(3), dp2x, dp2x2
    real(kind=rp) :: aexp, cfj, dis2, dp2, f12, fc, x2, f123, fa, fb , ori, x 
    integer(kind=ip) :: it(3), i, ic, itip, j, jj, k, m, n, inuc, is
    real(kind=rp), dimension(nmo,3) :: gun1
    real(kind=rp), dimension(nmo) :: gun

    rho = 0.0_rp
    grad = 0.0_rp
    gradmod = 0.0_rp
    gun = 0.0_rp
    gun1 = 0.0_rp

    inuc = 1
    !do inuc = 1,ncent --> change by shells in each center
      do is = 1,nshell(inuc)
        ic = atcenter(inuc,is)
        xcoor(:) = point(:) - xyz(ic,:)
        dis2 = xcoor(1)*xcoor(1) + xcoor(2)*xcoor(2) + xcoor(3)*xcoor(3)
        m = ishell(inuc,is)
        k = nuexp(ic,m,1)
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
          do j = 1,nmo
            cfj = coefnat(j,i)
            gun(j) = gun(j) + cfj*f123
            gun1(j,1) = gun1(j,1) + cfj*fa
            gun1(j,2) = gun1(j,2) + cfj*fb
            gun1(j,3) = gun1(j,3) + cfj*fc
          end do
        end do
 2     continue
      end do
    !end do

    rho = dot_product(occ,gun*gun)
    grad(1) = dot_product(occ,gun*gun1(:,1))
    grad(2) = dot_product(occ,gun*gun1(:,2))
    grad(3) = dot_product(occ,gun*gun1(:,3))

    grad(:) = grad(:) + grad(:)
    gradmod = norm2(grad)

    !if (present(norm)) then
    !  if (norm) grad = grad/gradmod
    !end if
 
  end subroutine

