module mod_gaunt

  use mod_prec, only: rp, ip
  implicit none
  private

  real(kind=rp), dimension(:,:), allocatable :: app
  real(kind=rp), dimension(:,:), allocatable :: bpp
  integer(kind=ip), dimension(:,:), allocatable :: islm

  public :: eval_gaunt

contains

  subroutine eval_gaunt (lmax, rint, coeff, output) bind(c)

    use mod_slm, only: init_slm, il0, jlm, deltam, &
                       allocate_space_for_slm, deallocate_space_for_slm
    implicit none
    integer(kind=ip), intent(in), value :: lmax
    real(kind=rp), intent(in) :: rint(3)
    real(kind=rp), intent(out) :: coeff(0:lmax*(lmax+2),0:lmax*(lmax+2))
    integer(kind=ip), intent(out) :: output(0:lmax*(lmax+2))

    real(kind=rp) :: tmp, term, r3
    integer(kind=ip) :: lmaxgaunt, lmax21
    integer(kind=ip) :: lmax2, lm1, lm2, l1, il1, l2, ll, ilm, l, lm, ier
    real(kind=rp), dimension(:), allocatable :: sql
    real(kind=rp), dimension(:,:), allocatable :: sgaunt

    lmax2 = lmax*(lmax+2)
    lmax21 = lmax2 + 1
    lmaxgaunt = lmax21*lmax21
    r3 = rint(1)*rint(1) + rint(2)*rint(2) + rint(3)*rint(3)
    r3 = sqrt(r3)

    call allocate_space_for_slm (lmax)
    call allocate_space_for_gaunt (lmax)
    allocate (sql(0:lmax2),stat=ier)
    if (ier.ne.0) stop 'cannot allocate sql'
    allocate (sgaunt(0:lmax,lmaxgaunt),stat=ier)
    sgaunt = 0.0
    call init_slm (lmax)
    call init_gaunt (lmax)
    ! Obtain coefficients of the multipolar expansion
    l = 0
    do lm = 0,lmax2
      sql(lm) = sqrt(real(l+l+1,rp))
      if (lm.eq.il0(l)+l) l = l + 1
    end do
    call sumgaunt (lmax,rint,sgaunt,lmaxgaunt)
    do lm1 = 0,lmax2
      l1 = jlm(lm1,1)
      il1 = (-1)**l1
      do lm2 = 0,lmax2
        l2 = jlm(lm2,1)
        ll = min(l1,l2)
        ilm = islm(lm1,lm2)
        tmp = 2.0*il1*deltam(l1,l2)*sgaunt(ll,ilm)
        term = tmp*sql(lm1)*sql(lm2)
        !if (okxcmp) corig(lm1,lm2)=term
        coeff(lm1,lm2) = term*(r3**(-1-l1-l2))
      end do
    end do
    output(:) = jlm(:,1)
    deallocate (sql, sgaunt)
    call deallocate_space_for_gaunt ()
    call deallocate_space_for_slm ()

  end subroutine eval_gaunt

  subroutine sumgaunt (lmax,rint,sgaunt,lmaxgaunt)

    use mod_slm, only: rsh, eps, il0
    use mod_param, only: pi
    implicit none
    integer(kind=ip), intent(in) :: lmax
    integer(kind=ip), intent(in) :: lmaxgaunt
    real(kind=rp), intent(out) :: sgaunt(0:lmax,lmaxgaunt)
    real(kind=rp), intent(in) :: rint(3)

    real(kind=rp) :: slm(0:4*lmax*(lmax+1))
    real(kind=rp) :: cp, ct, fac, r, rxy, sp, st
    integer(kind=ip) :: ii, l1, l2, l2max, l3, lm1, lm2, m1
    integer(kind=ip) :: lmax2, m2, m3, n, n3

    l2max = lmax+lmax
    lmax2 = lmax*(lmax+2)
 
    r = rint(1)*rint(1) + rint(2)*rint(2) + rint(3)*rint(3)
    r = sqrt(r)
    ct = rint(3)/r
    st = sqrt(1.0-ct*ct)
    rxy = sqrt(rint(1)*rint(1) + rint(2)*rint(2))
    if (rxy.lt.eps) then
      cp = 1.0
      sp = 0.0
    else
      cp = rint(1)/rxy
      sp = rint(2)/rxy
    end if
    call rsh (ct,st,cp,sp,slm,l2max)
 
    n = 0
    l1 = 0
    m1 = -l1
    do lm1 = 0,lmax2
      l2 = 0
      m2 = -l2 
      do lm2 = 0,lmax2
        n = n + 1
        islm(lm1,lm2) = n
        n3 = 0
        do l3 = abs(l1-l2),l1+l2,2
          fac = sqrt((l3+l3+1)/(4.0*pi)) !To transform to true slm
          ii = il0(l3)
          sgaunt(n3,n) = 0.0
          do m3 = -l3,l3
            sgaunt(n3,n) = sgaunt(n3,n)+slm(ii+m3)*dsupk(l1,m1,l2,m2,l3,m3)*fac
          end do
          n3 = n3 + 1
        end do
        if (m2.eq.l2) then
          l2 = l2 + 1
          m2 = -l2
        else
          m2 = m2 + 1
        end if
      end do
      if (m1.eq.l1) then
        l1 = l1 + 1
        m1 = -l1
      else
        m1 = m1 + 1
      end if
    end do

  end subroutine

  subroutine init_gaunt (lmax)

    implicit none
    integer(kind=ip), intent(in) :: lmax

    real(kind=rp) :: aux, bux, t1, t2
    integer(kind=ip) :: k1, k2, k20, k200, kk, kk0, kk00
    integer(kind=ip) :: mmp, n, mp, l, lp, m, m1 

    ! Computation of the app matrix 
    ! starting elements app(lm,00)(n) = delta(l,n)
    k1 = 0
    do l = 0,lmax
      do m = 0,l
        kk = k1*(k1+1)/2
        app(kk,l) = 1.0
        k1 = k1 + 1
      end do
    end do
    ! elements type app(lm,m'm')(n)
    do mp = 1,lmax
      k2 = mp*(mp+1)/2 + mp
      k20 = (mp-1)*mp/2 + mp-1
      do l = mp,lmax
        if (l.eq.mp) then
           m1 = mp
        else
           m1 = 0
        end if
        do m = m1,l
          k1 = l*(l+1)/2 + m
          kk = k1*(k1+1)/2 + k2
          kk0 = k1*(k1+1)/2 + k20
          do n = l-mp,l+mp,2
            if ( n.ge.m+mp) then
              app(kk,n) = (2*mp-1)*(app(kk0,n-1)/real(n+n-1,rp) - app(kk0,n+1)/real(n+n+3,rp))
            end if
          end do
        end do
      end do
    end do
    ! elements type app(lm,l'm')(n)
    do mp = 0,lmax
      do lp = mp+1,lmax
        k2 = lp*(lp+1)/2 + mp
        k20 = (lp-1)*lp/2 + mp
        k200 = (lp-2)*(lp-1)/2 + mp
        do l = lp,lmax
          if (l.eq.lp) then
            m1 = mp
          else
            m1 = 0
          end if
          do m = m1,l
            k1 = l*(l+1)/2 + m
            kk = k1*(k1+1)/2 + k2
            kk0 = k1*(k1+1)/2 + k20
            kk00 = k1*(k1+1)/2 + k200
            do n = l-lp,l+lp,2
              if (n.ge.m+mp) then
                aux = app(kk0,n+1)*(n+m+mp+1)/real(n+n+3,rp)
                if (n.gt.m+mp) aux = aux + app(kk0,n-1)*(n-m-mp)/real(n+n-1,rp)
                aux = aux*(lp+lp-1)
                if (lp.gt.mp+1) aux = aux - (lp+mp-1)*app(kk00,n)
                app(kk,n) = aux/real(lp-mp,rp)
              end if
            end do
          end do
        end do
      end do
    end do

    ! Computation bpp Matrix
    ! starting elements bpp(lm,00)(n) = delta(l,n)
    k1 = 0
    do l = 0,lmax
      do m = 0,l
        kk = k1*(k1+1) / 2
        bpp(kk,l) = 1.0
        k1 = k1 + 1
      end do
    end do
    ! elements type bpp(lm,m'm')(n)
    do mp = 1,lmax
      k2 = mp*(mp+1)/2 + mp
      k20 = (mp-1)*mp/2 + mp-1
      do l = mp,lmax
        if ( l.eq.mp ) then
          m1 = mp
        else
          m1 = 0
        end if
        do m = m1,l
          k1 = l*(l+1)/2 + m
          kk = k1*(k1+1)/2 + k2
          kk0 = k1*(k1+1)/2 + k20
          do n = l-mp,l+mp,2
            if (mp.gt.m) then
              t1 = 1.0
              t2 = 1.0
            else
              t1 = -(n-(m-mp+1))*(n-(m-mp+1)+1)
              t2 = -(n+(m-mp+1))*(n+(m-mp+1)+1)
            end if
            if (n.ge.abs(m-mp)) then
              if (n.eq.0) then
                bux = 0.0
              else
                bux = t1*bpp(kk0,n-1)/dfloat(n+n-1)
              end if
              bpp(kk,n) = (2*mp-1)*(bux-t2*bpp(kk0,n+1)/real(n+n+3,rp))
            end if
          end do
        end do
      end do
    end do
    ! elements type bpp(lm,l'm')(n)
    do mp = 0,lmax
      do lp = mp+1,lmax
        k2 = lp*(lp+1)/2 + mp
        k20 = (lp-1)*lp/2 + mp
        k200 = (lp-2)*(lp-1)/2 + mp
        do l = lp,lmax
          if (l.eq.lp) then
            m1 = mp
          else
            m1 = 0
          end if
          do m = m1,l
            k1 = l*(l+1)/2 + m
            kk = k1*(k1+1)/2 + k2
            kk0 = k1*(k1+1)/2 + k20
            kk00 = k1*(k1+1)/2 + k200
            do n = l-lp,l+lp,2
              mmp = abs(m-mp)
              if (n.ge.mmp) then
                aux = bpp(kk0,n+1)*(n+mmp+1)/real(n+n+3,rp)
                if (n.gt.mmp) aux = aux + bpp(kk0,n-1)*(n-mmp)/real(n+n-1,rp)
                aux = aux*(lp+lp-1)
                if (lp.gt.mp+1) aux = aux - (lp+mp-1)*bpp(kk00,n)
                bpp(kk,n) = aux/real(lp-mp,rp)
              end if
            end do
          end do
        end do
      end do
    end do

  end subroutine

  ! Evaluates the integral between three real spherical
  ! harmonics:
  !  l1 m1
  ! d            = < S      S      S     >
  !  l2 m2 l3 m3      l1 m1  l2 m2  l3 m3
  !
  ! This integral is symmetric with respect to the interchange of
  ! functions, i.e. (li,mi) pairs. It is sometimes called Gaunt
  ! integral or Gaunt coefficient. This routine calls csupk to
  ! calculate the \theta part, expecting to receive a
  ! Slater-Condon-Shortley coefficient between complex spherical
  ! harmonics, defined as:
  !
  !  l3                                
  ! C  (l1 m1, l2 m2) = sqrt(4*pi/2*l3+1) < Y     | Y     | Y     >
  !                                          l1 m1   l2 m2   l3 m3 
  !
  ! This version assumes that it will be called with (m1,m2,m3)
  ! values with non-zero I_{m1,m2,m3} (\phi integral);
  real(kind=rp) function dsupk (l1,m1,l2,m2,l3,m3)

    use mod_param, only: pi
    implicit none
    integer(kind=ip), intent(in) :: l1, m1, l2, m2, l3, m3

    integer(kind=ip) :: m1a, m2a, m3a, masum, mprod, msum
    real(kind=rp) :: dvalue
 
    ! Initialize variables
    mprod = m1*m2*m3
    if (mprod.lt.0) then
      dsupk = 0.0
      return
    end if
    msum = m1 + m2 + m3
    m1a = abs(m1)
    m2a = abs(m2)
    m3a = abs(m3)
    masum = m1a + m2a + m3a
    if (mprod.eq.0 .and. msum.eq.0 .and. masum.ne.0) then
      dsupk = 0.0
      return
    end if

    ! All m's are zero
    if (masum.eq.0) then
      dsupk = sqrt((l3+l3+1)/(4.0*pi))*csupk(l3,l1,0,l2,0)
      return
    ! |m1|=|m2|+|m3|
    else if (m1a.eq.m2a+m3a) then
      dvalue = sqrt((l3+l3+1)/(4.0*pi))*csupk(l3,l1,m1a,l2,m2a)
    ! |m2|=|m1|+|m3|
    else if (m2a.eq.m1a+m3a) then
      dvalue = sqrt((l3+l3+1)/(4.0*pi))*csupk(l3,l2,m2a,l1,m1a)
    ! |m3|=|m1|+|m2| 
    else if (m3a.eq.m1a+m2a) then
      dvalue = sqrt((l1+l1+1)/(4.0*pi))*csupk(l1,l3,m3a,l2,m2a)
    else
      dsupk = 0.0
      return
    end if

    ! All positive or two negative and one positive
    if (mprod.gt.0) then
      dvalue = dvalue/sqrt(2.0)
      ! If the positive value equals the negatives' sum
      if (msum.eq.0) dvalue = -dvalue
    end if

    ! Everything done
    dsupk = dvalue
    return

  end function

  ! Computes Slater-Condon csuperk coefficients using Gaunt expansion 
  ! for 3j symbols. The routine returns:
  !    k
  !   C (l m ,l m ) = sqrt(4*pi/2*j3+1) <j1m1|j3(m1-m2)|j2m2>
  !       1 1  2 2     
  !
  ! ji and mi are the associated quantum numbers
  real(kind=rp) function csupk(j3,j1,m1,j2,m2)
 
    use mod_slm, only: fact
    implicit none
    integer(kind=ip), intent(in) :: j3, j1, m1, j2, m2

    integer(kind=ip) :: m3, j, kmin, kmax, ik, k
    real(kind=rp) :: coef1, coef2, dumy
 
    ! Local variables
    m3 = m1-m2
    j = (j1+j2+j3)/2

    ! Gaunt Algorythm
    coef1 = fact(j1+j2-j3)*fact(j1-j2+j3)*fact(j2+j3-j1)*fact(j)
    coef2 = fact(j-j1)*fact(j-j2)*fact(j-j3)*fact(j+j+1)
    coef1 = coef1/coef2
    coef2 = fact(j1+m1)*fact(j2-m2)*fact(j2+m2)*fact(j3-m3)*fact(j3+m3)*fact(j1-m1)
    coef2 = sqrt(coef2*(2*j1+1)*(2*j2+1))
    coef1 = coef1*coef2
    coef2 = 0.0
    kmin = max(0,j2-j3+m1,j1-j3+m2)
    kmax = min(j1+j2-j3,j1+m1,j2+m2)
    coef2 = 0.0

    ik = (-1)**kmin
    do k = kmin,kmax
      dumy = fact(k)*fact(j1+j2-j3-k)*fact(j1+m1-k)*fact(j2+m2-k)*fact(j3-j2-m1+k)*fact(j3-j1-m2+k)
      coef2 = coef2+ik/dumy
      ik = -ik
    end do
    csupk = coef1*coef2*(-1)**(m2+j3+j)
    return

  end function

  subroutine allocate_space_for_gaunt (maxl)
  integer(kind=ip), intent(in) :: maxl
  integer(kind=ip) :: mxlcof, mxkcof
  integer(kind=ip) :: ier
  mxlcof = maxl*(maxl+3)/2
  mxkcof = mxlcof*(mxlcof+3)/2 
  allocate (app(0:mxkcof,0:2*maxl+1) ,stat=ier) 
  if (ier.ne.0) stop "cannot alloc memory"
  allocate (bpp(0:mxkcof,0:2*maxl+1) ,stat=ier) 
  if (ier.ne.0) stop "cannot alloc memory"
  allocate (islm(0:maxl*(maxl+2),0:maxl*(maxl+2)),stat=ier)
  if (ier.ne.0) stop "cannot alloc memory"
  end subroutine allocate_space_for_gaunt
  
  subroutine deallocate_space_for_gaunt ()
  integer(kind=ip) :: ier
  deallocate (app,stat=ier) 
  deallocate (bpp,stat=ier) 
  deallocate (islm,stat=ier) 
  end subroutine deallocate_space_for_gaunt

end module mod_gaunt 
