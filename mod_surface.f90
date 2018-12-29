module mod_surface

  use mod_prec, only: rp, ip
  implicit none
  private

  integer(kind=ip), parameter, public :: minter = 11 
  integer(kind=ip), parameter, public :: maxtrial = 21
  ! The MAXSTART parameter has to do with the RSEARCH order, 
  ! RSEARCH nsearch nrstart, (rstart(i),i=1,nrstart)
  ! This order can be used when one wants to avoid using the
  ! default starting values of the radial coordinate in
  ! the initial search of interatomic surfaces. If nsearch is 0
  ! the starting values (rstart(i),i=1,nrstart) are used for all
  ! the atoms. If nsearch is different from 0, the rstart(i)
  ! (i=1,nrstart) values are used only for the atom 'nsearch'
  ! and all its symmetry equivalent atoms.
  integer(kind=ip), parameter, public :: maxstart = 40

  integer(kind=ip), public :: inuc_
  real(kind=rp), public :: xnuc_(3)
  real(kind=rp), allocatable, dimension(:,:), public :: xyzrho_
  real(kind=rp), public :: rmaxsurf_
  integer(kind=ip), public :: npang_
  real(kind=rp), allocatable, dimension(:), public :: rstart_
  real(kind=rp), allocatable, dimension(:,:), public :: rlimsurf_
  integer(kind=ip), allocatable, dimension(:), public :: nlimsurf_

  ! options
  real(kind=rp), public :: step_
  integer(kind=ip), public :: mstep_
  ! Precision in interatomic surfaces. 
  real(kind=rp), public :: epsilon_
  real(kind=rp), public :: epsroot_
  ! If a non nuclear maximum is found whose distance to any nucleus 
  ! is larger than EPSISCP, promolden STOPS. This control of non nuclear
  ! maxima is performed in the 'iscp.f' routine. When Hydrogen atoms
  ! are involved it is very convenient to change the default value
  ! for EPSISCP (0.08) to a smaller value by using the 'EPSISCP value'
  ! order in the input file.
  real(kind=rp), public :: epsiscp_
  ! NTRIAL  = number of sub-intervals in the search of the surface.
  ! RPRIMER = First point in the in the search of the surface.
  integer(kind=ip), public :: ntrial_
  real(kind=rp), public :: rprimer_
  logical, public :: rotgrid_
  real(kind=rp), public :: angx_, angy_, angz_
  integer(kind=ip), public :: steeper_

  public :: init_surf, surf, rotagrid, odeint, rkck
  public :: allocate_space_for_surface, deallocate_space_for_surface

contains

! Determine the limit of the zero flux surface from xpoint along
! the theta phi direction. 
subroutine surf (cot,sit,cop,sip,rsurf,nsurf)

  use mod_io, only: ferror, faterr
  use mod_mole, only: ncent_
  implicit none
 
  real(kind=rp), parameter :: half = 0.5_rp
 
  integer(kind=ip), intent(out) :: nsurf
  real(kind=rp), intent(in) :: cot, sit, cop, sip
  real(kind=rp), intent(out) :: rsurf(minter,2) 
 
  logical :: inf, good
  integer(kind=ip) :: nintersec, nt
  integer(kind=ip) :: isurf(minter,2), i, ia, ib, im, ncount
  real(kind=rp) :: sintcosp, sintsinp
  real(kind=rp) :: xin(3), xfin(3), xpoint(3), xmed(3), ra, ract
  real(kind=rp) :: xdeltain(3), xsurf(0:minter,3), rb, rm, cost
  ! GEOFAC  = ((rmaxsurf-0.1d0)/rprimer)**(1d0/(ntrial-1))
  !           (RMAXSURF is defined in passed in integ.inc)
  real(kind=rp) :: geofac 
 
  cost = cot
  sintcosp = sit*cop
  sintsinp = sit*sip
 
  if (ncent_.eq.1) then
    nsurf = 1
    rsurf(1,1) = 0.0_rp
    rsurf(1,2) = rmaxsurf_
    return
  end if
 
  inf = .false.
  ncount = 0_ip
  nintersec = 0_ip
  ia = inuc_
  ra = 0.0_rp
  geofac = ((rmaxsurf_-0.1_rp)/rprimer_)**(1.0_rp/(ntrial_-1))
  do i = 1,ntrial_
    ract = rprimer_*geofac**(i-1)
    xdeltain(1) = ract*sintcosp
    xdeltain(2) = ract*sintsinp
    xdeltain(3) = ract*cost    
    xpoint(:) = xnuc_(:) + xdeltain(:)
    call odeint (xpoint,0.1_rp,1.0_rp,inf,epsilon_)
    good = iscp(xpoint,ib)
    rb = ract
    if (ib.ne.ia .and. (ia.eq.inuc_ .or. ib.eq.inuc_)) then
      if (ia.ne.inuc_ .or. ib.ne.-1) then
        nintersec = nintersec + 1_ip
        if (nintersec.gt.minter) then
          call ferror('mod_surface', 'increase minter in mod_surf', faterr)
        end if
        xsurf(nintersec,1) = ra
        xsurf(nintersec,2) = rb
        isurf(nintersec,1) = ia
        isurf(nintersec,2) = ib
      end if
    end if
    ia = ib
    ra = rb
  end do        ! looking for intersections

  ! We have now a consistent set of trial points. 
  ! Let us refine by bipartition, consistency check added. 
  ! No other nuclei basins can be found other than those explicit in isurf
  do i = 1,nintersec
    ia = isurf(i,1)
    ib = isurf(i,2)
    ra = xsurf(i,1)
    rb = xsurf(i,2)
    xin(1) = xnuc_(1) + ra*sintcosp
    xin(2) = xnuc_(2) + ra*sintsinp
    xin(3) = xnuc_(3) + ra*cost
    xfin(1) = xnuc_(1) + rb*sintcosp
    xfin(2) = xnuc_(2) + rb*sintsinp
    xfin(3) = xnuc_(3) + rb*cost
    do while (abs(ra-rb).gt.epsilon_)
      ! Mean Value 
      xmed(:) = half*(xfin(:)+xin(:))    
      rm = (ra+rb)*half
      xpoint(:) = xmed(:)
      call odeint (xpoint,0.1_rp,1.0_rp,inf,epsilon_)
      good = iscp(xpoint,im)
      ! bipartition
      if (im.eq.ia) then
        xin(:) = xmed(:)
        ra = rm
      else if (im.eq.ib) then
        xfin(:) = xmed(:)
        rb = rm
      else
        if (ia.eq.inuc_) then
          xfin(:) = xmed(:)
          rb = rm
        else
          xin(:) = xmed(:)
          ra = rm
        end if
      end if
    end do
    ! found. Mean value
    xpoint(:) = half*(xfin(:)+xin(:))    
    xsurf(i,3) = (ra+rb)*half
  end do

  ! organize pairs
  nsurf = 0_ip
  xsurf(0,3) = 0.0_rp
  ia = inuc_
  nsurf = nintersec
  do i = 1,nsurf
    rsurf(i,2) = xsurf(i,3)
  end do
  nt = mod(nintersec,2)
  if (nt.eq.0) then
    nsurf = nsurf + 1_ip
    rsurf(nsurf,2) = rmaxsurf_
  end if
  write(*,*) cot,sit,cop,sip,nsurf,rsurf(1:nsurf,2)
 
end subroutine

  subroutine odeint (ystart,h1,iup,inf,eps)
 
    use mod_io, only: faterr, ferror, string
    use mod_fields, only: pointr1=>density_grad_shell
    implicit none

    ! Parameters
    real(kind=rp), parameter :: tiny = 1d-40
    real(kind=rp), parameter :: epsg = 1d-15
    real(kind=rp), parameter :: epsg1 = 1d-15

    ! Arguments
    real(kind=rp), intent(inout) :: ystart(3)
    real(kind=rp), intent(in) :: h1
    real(kind=rp), intent(in) :: iup
    logical, intent(inout) :: inf
    real(kind=rp), intent(in) :: eps

    ! Local vars
    integer(kind=ip) :: nstp, nuc
    real(kind=rp), parameter :: hmin = 0.0_rp ! minimal step size
    real(kind=rp) :: x1 ! intial point
    real(kind=rp) :: x2 ! final point
    real(kind=rp) :: p(3) ! initial solution point
    real(kind=rp) :: h ! initial step size
    real(kind=rp) :: x ! update point
    real(kind=rp) :: hnext ! next steep size
    real(kind=rp) :: dydx(3), y(3), yscal(3)
    real(kind=rp) :: a1, a2, a3
    real(kind=rp) :: rho, grad(3), gradmod
 
    inf = .false.
    p(1) = ystart(1)
    p(2) = ystart(2)
    p(3) = ystart(3)
    call pointr1 (p,rho,grad,gradmod)
    if (gradmod.lt.epsg .and. rho.lt.epsg1) then
      inf = .true.
      return
    end if
 
    x1 = 0.0_rp
    x2 = 1d40*iup
    x = x1 ! initial point
    !h = sign(h1,x2-x1) ! initial steep size
    h = min(h1,x2-x1) ! initial steep size
    y(:) = ystart(:) ! initial point 
 
    do nstp = 1,mstep_
      call pointr1 (y,rho,grad,gradmod)
      dydx(:) = grad(:)
      yscal(:) = max(abs(y(:))+abs(h*dydx(:))+tiny,eps)
      if ((x+h-x2)*(x+h-x1).gt.0.0_rp) h = x2 - x
      call rkqs (y,dydx,x,h,eps,yscal,hnext,steeper_)
      if ((x-x2)*(x2-x1).ge.0.0_rp .or. iscp(y,nuc)) then
        ystart(:) = y(:)
        return
      end if
      if (abs(hnext).lt.hmin) then
        call ferror ('mod_surface/odeint', 'stepsize small than minimum', faterr)
      end if
      if (nstp.eq.mstep_) then
        call ferror ('mod_surface/odeint', 'reached maxstp', faterr)
      end if 
      h = hnext
    end do

    ! Test if the point is far from RMAXSURF from current atom. 
    a1 = y(1) - xnuc_(1)
    a2 = y(2) - xnuc_(2)
    a3 = y(3) - xnuc_(3)
    if ((a1*a1+a2*a2+a3*a3).ge.5d0*5d0) then
      inf = .true.
      return
    else
      call ferror ('mod_surface/odeint', 'Non nuclear maximum at : ' &
                                         //string(y(1),'e')//' '    &  
                                         //string(y(2),'e')//' '    &  
                                         //string(y(3),'e'), faterr) 
    end if
 
  end subroutine 

  subroutine rkqs (y,dydx,x,htry,eps,yscal,hnext,steeper)
      
    use mod_io, only: ferror, faterr
    implicit none

    ! Parameters
    real(kind=rp), parameter :: safety = 0.9_rp
    real(kind=rp), parameter :: pgrow = -0.2_rp
    real(kind=rp), parameter :: pshrnk = -0.25_rp
    real(kind=rp), parameter :: errcon = 1.89d-4
    
    ! Arguments
    integer(kind=ip), intent(in) :: steeper
    real(kind=rp), dimension(3), intent(in) :: dydx
    real(kind=rp), dimension(3), intent(inout) :: y
    real(kind=rp), dimension(3), intent(in) :: yscal
    real(kind=rp), intent(inout) :: x
    real(kind=rp), intent(in) :: eps
    real(kind=rp), intent(in) :: htry
    real(kind=rp), intent(out) :: hnext

    ! Local vars
    integer(kind=ip) :: i
    real(kind=rp), dimension(3) :: yerr, ytemp
    real(kind=rp) :: h, errmax, htemp, xnew
      
    h = htry
    hnext = 0.0_rp
    errmax = 0.0_rp
    yerr = 0.0_rp

    do
      call rkck (y,dydx,h,ytemp,yerr)
      ! adaptive and error estimation 
      errmax = 0.0_rp
      do i = 1,3
        errmax = max(errmax,abs(yerr(i)/yscal(i)))
      end do
      errmax = errmax/eps
      if (errmax.gt.1.0_rp) then
        htemp = safety*h*(errmax**pshrnk)
        !h = sign(max(abs(htemp),0.1_rp*abs(h)),h)
        h = min(max(abs(htemp),0.1_rp*abs(h)),h)
        xnew = x + h
        if (xnew.eq.x) then
          call ferror ('mod_odeint/rkqs', 'stepsize underflow', faterr)
          return
        end if
        cycle
      else
        if (errmax.gt.errcon) then
          hnext = safety*h*(errmax**pgrow)
        else
          hnext = 5.0_rp*h
        end if
        x = x + h
        y = ytemp
        return
      end if
    end do
 
  end subroutine rkqs

  ! Runge-Kutta-Cash-Karp embedded 4(5)-order, with local extrapolation.
  ! Runge-Kutta embedded 4th order, Cash-Karp parametrization.
  ! ##### 6 stages, 5th order
  ! This scheme is due to Cash and Karp, see [1].
  !    
  !   0    | 0
  !   1/5	 | 1/5
  !   3/10 | 3/40	         9/40
  !   3/5	 | 3/10	         -9/10	      6/5
  !   1	   | -11/54	       5/2	        -70/27	    35/27
  !   7/8	 | 1631/55296    175/512      575/13824   44275/110592     253/4096     0
  !  ----------------------------------------------------------------------------------------
  !        | 37/378        0           250/621      125/594          0            512/1771
  !        | 2825/27648    0           18575/48384  13525/55296      277/14336    1/4
  ! 
  ! [1] *A variable order Runge-Kutta method for initial value problems with rapidly varying right-hand sides*, J. R. Cash,
  ! A. H. Karp, ACM Transactions on Mathematical Software, vol. 16,  pp. 201--222, 1990, doi:10.1145/79505.79507.
  subroutine rkck (y,dydx,h,yout,yerr)
      
    use mod_fields, only: pointr1=>density_grad_shell
    implicit none
  
    ! Arguments
    real(kind=rp), intent(in) :: dydx(3)
    real(kind=rp), intent(in) :: y(3)
    real(kind=rp), intent(out) :: yerr(3)
    real(kind=rp), intent(out) :: yout(3)
    real(kind=rp), intent(in) :: h
  
    ! Local vars
    real(kind=rp) :: rho, grad(3), gradmod
    real(kind=rp) :: ak2(3), ak3(3), ak4(3), ak5(3), ak6(3)
    real(kind=rp), parameter :: b21=0.2d0,                           &
                                b31=3.0d0/40.0d0,                    &
                                b32=9.0d0/40.0d0,                    &
                                b41=0.3d0,b42=-0.9d0,b43=1.2d0,      &
                                b51=-11.0d0/54.0d0,b52=2.5d0,        &
                                b53=-70.0d0/27.0d0,b54=35.d0/27.0d0, &
                                b61=1631.0d0/55296.0d0,              &
                                b62=175.0d0/512.0d0,                 &
                                b63=575.0d0/13824.0d0,               &
                                b64=44275.0d0/110592.0d0,            &
                                b65=253.0d0/4096.0d0                  
    real(kind=rp), parameter :: c1=37.0d0/378.0d0,c3=250.0d0/621.0d0,&
                                c4=125.0d0/594.0d0,                  &
                                c6=512.0d0/1771.0d0                   
    real(kind=rp), parameter :: dc1=c1-2825.0d0/27648.0d0,           &
                                dc3=c3-18575.0d0/48384.0d0,          &
                                dc4=c4-13525.0d0/55296.0d0,          &
                                dc5=-277.0d0/14336.0d0,              &
                                dc6=c6-0.25d0                         
   
    yout = y + b21*h*dydx
   
    call pointr1 (yout,rho,grad,gradmod)
    ak2(:) = grad(:)
    yout = y + h*(b31*dydx+b32*ak2)

    call pointr1 (yout,rho,grad,gradmod)
    ak3(:) = grad(:)
    yout = y + h*(b41*dydx+b42*ak2+b43*ak3)
  
    call pointr1 (yout,rho,grad,gradmod)
    ak4(:) = grad(:)
    yout = y + h*(b51*dydx+b52*ak2+b53*ak3+b54*ak4)
  
    call pointr1 (yout,rho,grad,gradmod)
    ak5(:) = grad(:)
    yout = y + h*(b61*dydx+b62*ak2+b63*ak3+b64*ak4+b65*ak5)

    call pointr1 (yout,rho,grad,gradmod)
    ak6(:) = grad(:)
    yout = y + h*(c1*dydx+c3*ak3+c4*ak4+c6*ak6)
   
    yerr = h*(dc1*dydx+dc3*ak3+dc4*ak4+dc5*ak5+dc6*ak6)
 
  end subroutine

  logical function  iscp (p,nuc)
 
    use mod_mole, only: ncent_
    use mod_fields, only: pointr1=>density_grad_shell
    implicit none
 
    integer(kind=ip), intent(out) :: nuc
    real(kind=rp), intent(in) :: p(3)
 
    real(kind=rp) :: grad(3), gradmod, rho
    real(kind=rp) :: x(3)
    integer(kind=ip) :: i
 
    iscp = .false.
    nuc = 0
    x(:) = p(:)
    call pointr1 (x,rho,grad,gradmod)
    do i = 1,ncent_
      if (abs(p(1)-xyzrho_(1,i)).lt.epsiscp_ .and. &
          abs(p(2)-xyzrho_(2,i)).lt.epsiscp_ .and. &
          abs(p(3)-xyzrho_(3,i)).lt.epsiscp_) then
          iscp = .true.
          nuc = i
      end if
    end do
 
    if (gradmod.le.1d-10) then
      iscp = .true.
      if (rho.le.1d-10) nuc = -1
    end if
 
  end function

! This routine performs three successive rotations R_z, R_y, R_x
! about the Z, Y, and X axes, respectively, of NANG points placed on 
! an unit radius sphere, and defined by the polar coordinates 
! cos(theta)=ct, sin(theta)=st, cos(phi)=cp, and sin(phi)=sp.
! input
! angx, angy, angz = Rotation angle about the X, Y, and Z axes.
! ct,st,cp,sp = cos(theta),sin(theta),cos(phi),sin(phi). 
! nang = number of points on the unit radius sphere.
! maxang = maximum number of points.
! output
! ct,st,cp,sp = New polar coordinates after rotations.
subroutine rotagrid (ct,st,cp,sp,npang)

  implicit none
      
  real(kind=rp), parameter :: eps = 1.0d-07
  real(kind=rp), parameter :: zero = 0.0_rp
  real(kind=rp), parameter :: one = 1.0_rp

  integer(kind=ip), intent(in) :: npang
  real(kind=rp), intent(inout) :: ct(npang),st(npang),cp(npang),sp(npang)

  real(kind=rp) :: xmat(3,3),ymat(3,3),zmat(3,3),rot(3,3)
  real(kind=rp) :: prod,x,y,z,xp,yp,zp,r,rxy,cangx,cangy,cangz,sangx,sangy,sangz
  integer(kind=ip) :: i,j,k

  ! Rotacion matrix about X axis.
  xmat = zero
  xmat(1,1) = one
  cangx = cos(angx_)
  sangx = sin(angx_)
  xmat(2,2) = +cangx
  xmat(3,2) = +sangx
  xmat(2,3) = -sangx
  xmat(3,3) = +cangx

  ! Rotacion matrix about Y axis.
  ymat = zero
  ymat(2,2) = one
  cangy = cos(angy_)
  sangy = sin(angy_)
  ymat(1,1) = +cangy
  ymat(3,1) = +sangy
  ymat(1,3) = -sangy
  ymat(3,3) = +cangy

  ! Rotacion matrix about Z axis.
  zmat = zero
  zmat(3,3) = one
  cangz = cos(angz_)
  sangz = sin(angz_)
  zmat(1,1) = +cangz
  zmat(2,1) = +sangz
  zmat(1,2) = -sangz
  zmat(2,2) = +cangz

  ! Full rotacion matrix. R = R_X * R_Y * R_Z
  do i = 1,3
    do j = 1,3
      prod = zero
      do k = 1,3
        prod = prod + ymat(i,k)*zmat(k,j)
      end do
      rot(i,j) = prod
    end do
  end do
  zmat(1:3,1:3) = rot(1:3,1:3)
  do i = 1,3
    do j = 1,3
      prod = zero
      do k = 1,3
        prod = prod + xmat(i,k)*zmat(k,j)
      end do
      rot(i,j) = prod
    end do
  end do

  ! Rotate angles.
  do i = 1,npang
    x = st(i)*cp(i)
    y = st(i)*sp(i)
    z = ct(i)
    xp = rot(1,1)*x + rot(1,2)*y + rot(1,3)*z
    yp = rot(2,1)*x + rot(2,2)*y + rot(2,3)*z
    zp = rot(3,1)*x + rot(3,2)*y + rot(3,3)*z
    ! Rebuild CT,ST,CP, and SP
    rxy = xp*xp + yp*yp
    r = sqrt(rxy+zp*zp)
    if (rxy.lt.eps) then
      if (zp.ge.zero) then
        ct(i) = +one
      else
        ct(i) = -one
      endif
      st(i) = zero
      sp(i) = zero
      cp(i) = one
    else
      rxy = sqrt(rxy)
      ct(i) = zp/r
      st(i) = sqrt((one-ct(i))*(one+ct(i)))
      cp(i) = xp/rxy
      sp(i) = yp/rxy
    end if
  end do
            
end subroutine

  subroutine init_surf ()

    implicit none

    steeper_ = 1_ip
    rotgrid_ = .false.
    ntrial_ = 11_ip
    rprimer_ = 0.4_rp
    inuc_ = 0_ip
    epsilon_ = 1d-5 
    epsroot_ = 1d-4
    epsiscp_ = 0.08_rp
    angx_ = 0.0_rp
    angy_ = 0.0_rp
    angz_ = 0.0_rp
    npang_ = 5810
    rmaxsurf_ = 10.0_rp
    step_ = 0.1_rp
    mstep_ = 100_ip

  end subroutine
                                                                        
  subroutine allocate_space_for_surface (ncent,nangular,ncutsurf)

    use mod_memory, only: alloc
    implicit none
    integer(kind=ip), intent(in) :: nangular, ncutsurf, ncent

    call alloc ('mod_surface', 'xyzrho', xyzrho_, 3, ncent)
    call alloc ('mod_surface', 'rstart', rstart_, maxstart)
    call alloc ('mod_surface', 'rlimsurf', rlimsurf_, nangular, ncutsurf)
    call alloc ('mod_surface', 'nlimsurf', nlimsurf_, nangular)

  end subroutine allocate_space_for_surface

  subroutine deallocate_space_for_surface

    use mod_memory, only: free
    implicit none

    call free ('mod_surface', 'xyzrho', xyzrho_)
    call free ('mod_surface', 'rstart', rstart_)
    call free ('mod_surface', 'rlimsurf', rlimsurf_)
    call free ('mod_surface', 'nlimsurf', nlimsurf_)

  end subroutine deallocate_space_for_surface

end module mod_surface
