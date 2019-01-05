! Real spherical harmonic (RSH) S_{l,m} will be stored as
! element (l*(l+1)+m) of a supervector slm(0:maxl(maxl+2))
!
! rshnorm() ....... normalization constant of the RSH
! il0() ........... index of the S_{l,0} RSH (l*(l+1))
! pil0() .......... pil0(k) = Number of RSH with l < k
!
module mod_slm

  use mod_prec, only: rp, ip
  implicit none
  private

  integer(kind=ip), parameter :: mxfact = 100
  integer(kind=ip), parameter :: mcomb = 12

  real(kind=rp), allocatable, dimension(:) :: rshnorm, sqrin
  real(kind=rp), allocatable, dimension(:,:) :: deltam
  integer(kind=ip), allocatable, dimension(:) :: il0, pil0
  integer(kind=ip), allocatable, dimension(:,:) :: jlm, islm
  
  real(kind=rp) :: fact(0:mxfact)
  real(kind=rp) :: comb(0:mcomb,0:mcomb)

  public :: init_slm, eval_rsh, allocate_space_for_slm, deallocate_space_for_slm

contains

  ! ct - cos(theta), st-sin(theta)
  ! cf - cos(phi),   sf-sin(phi)
  ! rsh - computes real spherical harmonics up to l=lmax for a
  ! relative position (x,y,z), distance r; results are stored as
  ! a supervector in slm, where S_{l,m} (\hat{r}) is stored in the
  ! slm(l*(l+1)+m) element
  ! The routine WILL NOT test if r is zero; if it is, it will return
  ! the same value given for \theta=0
  subroutine eval_rsh (ct,st,cf,sf,slm,lmax)
 
    implicit none
    real(kind=rp), intent(in) :: ct, st, cf, sf
    integer(kind=ip), intent(in) :: lmax
    real(kind=rp), intent(out) :: slm (0:lmax*(lmax+2))
 
    real(kind=rp), parameter :: eps = 1d-7

    integer(kind=ip) :: i, ifac, il, im, l, l2, ll, m, m2, ipp
    real(kind=rp) :: a1, a2, cc, facsqrt, fifac, fifac1, ss, tmp

    slm = 0.0
    ! test if (x,y,z) points in the z direction
    if (abs(abs(ct)-1.0).le.eps) then
      ! positive z => \theta=0
      if (ct.ge.0.0) then
        do l = 0,lmax
          slm(il0(l)) = rshnorm(il0(l))
        end do
      ! negative z => \theta=\pi
      else 
        do l = 0,lmax-1,2
          slm(il0(l)) = rshnorm(il0(l))
          slm(il0(l+1)) = - rshnorm(il0(l+1))
        end do
        if (mod(lmax,2).eq.0) then
          slm(il0(lmax)) = rshnorm(il0(lmax))
        end if
      end if
      ! RSH with m<>0 are 0
      i = -1 
      do l = 0,lmax-1
        i = i + 1 
        do m = 0,l+l
          i = i + 1
          slm(i) = 0.0
        end do
      end do 
      i = i + 1 
      do m = 1,lmax
        i = i + 1
        slm(i) = 0.0
      end do
    ! (x,y,z) pointing elsewhere
    else
      ! sines and cosines
      slm(0) = 1.0
      ifac = 1
      ! compute associate Legendre polinomials
      do l = 1, lmax
        m = l - 1
        im = il0(m) + m
        il = il0(l) + l
        fifac = sqrin(ifac)
        fifac1 = sqrin(ifac+1)
        facsqrt = fifac/fifac1
        slm(il) = facsqrt*st*slm(im)
        slm(il-1) = fifac*ct*slm(im)
        l2 = ifac
        ifac = ifac + 2
        do ll = l+1,lmax
          l2 = l2 + 2
          a1 = sqrin(ll-m)/sqrin(ll+m)
          a2 = sqrin(ll-m-1)/sqrin(ll+m-1)
          slm(il0(ll)+m) = (ct*l2*a1*slm(il0(ll-1)+m)-(ll+m-1)*a1*a2*slm(il0(ll-2)+m))/(ll-m)
        end do
      end do
      ! compute \Phi functions and unnormalized RSH
      cc = 1.0
      ss = 0.0
      do m = 1, lmax
        tmp = cc*cf - ss*sf
        ss = ss*cf + cc*sf
        cc = tmp
        m2 = m+m
        ipp = il0(m-1)+m
        do l2 = m2,lmax+lmax,2
          ipp = ipp + l2
          slm(ipp-m2) = slm(ipp)*ss
          slm(ipp) = slm(ipp)*cc
        end do
      end do
      ! normalize RSH
      do i = 0,lmax*(lmax+2)
        slm(i) = slm(i)*rshnorm(i)
      end do
    end if

  end subroutine

  ! Compute the sqrt[(2l+1)/(4 pi)] factor of the normalization 
  ! constants of  RSH.
  ! In this version the slm stored are already multiplied by
  ! the sqrt(4pi/(2l+1)) factor
  subroutine init_slm (lmaxi)

    use mod_param, only: pi
    implicit none
    integer(kind=ip), intent(in) :: lmaxi

    real(kind=rp), parameter :: twopi = 2.0*pi

    real(kind=rp) :: temp
    integer(kind=ip) :: l, m, i, il, l1, l12, l2, l22, lmax, ls, j, ij

    fact(0) = 1.0
    do i = 1,mxfact
      fact(i) = fact(i-1)*real(i,rp)
    end do
    do i = 0,mcomb
      do j = 0,i/2
        ij = i - j
        comb(i,j) = fact(i)/fact(j)/fact(ij)
        if (j.ne.ij) comb(i,ij) = comb(i,j)
      end do
    end do

    l = 0
    m = 0
    do i = 0,lmaxi*(lmaxi+2)
      jlm(i,1) = l
      jlm(i,2) = m
      if (m.eq.l) then
         l = l + 1
         m = -l
      else
         m = m + 1
      end if
    end do
 
    lmax = lmaxi + lmaxi 
    il0(0) = 0
    pil0(0) = 0
    rshnorm(0) = 1.0
    il = 0
    do l = 1,lmax
      il = il + l + l
      il0(l) = il
      pil0(l) = pil0(l-1) + l + l - 1
      rshnorm(il) = 1.0
      do m = 1,l
        rshnorm(il+m) = rshnorm(il) * sqrt(2.0)
        rshnorm(il-m) = rshnorm(il+m)
      end do
    end do

    ! Compute square root of integer from 0 to (2*lmax+1).
    do i = 0,2*lmax+1
      sqrin(i) = sqrt(dble(i))
    end do

    ! Compute deltam
    lmax = lmaxi
    do l1 = 0,lmax
       l12 = l1 + l1
       do l2 = 0,lmax
          l22 = l2 + l2
          ls = l1 + l2
          temp = (-1)**ls*twopi
          deltam(l1,l2) = temp*fact(l1)*fact(l2)*fact(l12+l22)/(fact(ls)*fact(l12+1)*fact(l22+1))
       end do
    end do
 
  end subroutine

  subroutine allocate_space_for_slm (maxl)
  integer(kind=ip), intent(in) :: maxl
  integer(kind=ip) :: maxl2
  integer :: ier
  maxl2=maxl+maxl
  allocate (rshnorm(0:maxl2*(maxl2+2)),stat=ier) 
  allocate (sqrin(0:2*maxl2+1),stat=ier) 
  allocate (deltam(0:maxl,0:maxl),stat=ier) 
  allocate (il0(0:maxl2),stat=ier) 
  allocate (pil0(0:maxl2),stat=ier)
  allocate (jlm(0:maxl*(maxl+2),2),stat=ier) 
  allocate (islm(0:maxl*(maxl+2),0:maxl*(maxl+2)),stat=ier) 
  end subroutine allocate_space_for_slm

  subroutine deallocate_space_for_slm
  integer :: ier
  deallocate (rshnorm,stat=ier) 
  deallocate (sqrin,stat=ier) 
  deallocate (deltam,stat=ier) 
  deallocate (il0,stat=ier) 
  deallocate (pil0,stat=ier)
  deallocate (jlm,stat=ier) 
  deallocate (islm,stat=ier) 
  end subroutine deallocate_space_for_slm

end module mod_slm

!allocate (slm(0:lmax*(lmax+2),npang),stat=ier)
!call initslm (lmax) !for bicentric
!do j = 1,npang
!  cost = ct(j)
!  sint = st(j)
!  cosp = cp(j)
!  sinp = sp(j)
!  call rsh (cost,sint,cosp,sinp,slm(0,j),lmax)
!end do
