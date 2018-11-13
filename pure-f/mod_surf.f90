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
module mod_surf

  use mod_prec, only: rp, ip
  implicit none
  private

  real(kind=rp), parameter :: hminimal = 1e-6
  real(kind=rp), parameter :: epsrho = 1e-10
  real(kind=rp), parameter :: epsgrad = 1e-10
  integer(kind=ip), parameter :: maxstp = 500
  integer(kind=ip), parameter :: minter = 10 
  integer(kind=ip), parameter :: maxtrial = 13

  integer(kind=ip) :: inuc
  real(kind=rp) :: xnuc(3)

  real(kind=rp), allocatable, dimension(:,:), public :: xyzrho
  real(kind=rp), allocatable, dimension(:,:), public :: rlimsurf
  integer(kind=ip), allocatable, dimension(:), public :: nlimsurf

  ! options
  integer(kind=ip), dimension(4), public :: nangleb
  real(kind=rp), public :: rmaxsurf
  real(kind=rp), public :: epssurf
  ! If a non nuclear maximum is found whose distance to any nucleus 
  ! is larger than EPSISCP, promolden STOPS.
  real(kind=rp), public :: epsiscp
  ! NTRIAL  = number of sub-intervals in the search of the surface.
  ! RPRIMER = First point in the in the search of the surface.
  integer(kind=ip), public :: ntrial
  real(kind=rp), public :: rprimer
  integer(kind=ip), public :: steeper

  public :: optssurf, init_surf, surface

contains

  subroutine surface(idx)

    use iso_fortran_env, only: uout=>output_unit
    !$ use omp_lib, only: omp_get_wtime
    use mod_param, only: verbose, filename
    use mod_memory, only: free, alloc
    use mod_io, only: mline, string
    use mod_wfn, only: ncent
    implicit none

    integer(kind=ip), intent(in) :: idx

    character(len=4) :: d4
    character(len=mline) :: filesurf
    integer(kind=ip) :: npang, nsurf, ltxt, lsu, j, k
    real(kind=rp) :: time1, time2, rmins, rmaxs, rsurf(minter,2)
    real(kind=rp), allocatable, dimension(:) :: ct, st, cp, sp, angw

    interface
      subroutine lebgrid (ct,st,cp,sp,w,npoints)
        import rp, ip
        integer(kind=ip), intent(in) :: npoints
        real(kind=rp), intent(out) :: ct(npoints)
        real(kind=rp), intent(out) :: st(npoints)
        real(kind=rp), intent(out) :: cp(npoints)
        real(kind=rp), intent(out) :: sp(npoints)
        real(kind=rp), intent(out) :: w(npoints)
      end subroutine
    end interface

    ! Init
    filesurf = trim(filename)//".surf"
    call allocate_space_for_cp (ncent)
    call findnuc ()
    inuc = idx
    d4 = string(inuc,length=4,pad0=.true.)
    xnuc = xyzrho(inuc,:)
    if (verbose) call info_surf ()

    ! Begin
    npang = nangleb(2)
    call alloc ('surface', 'ct', ct, npang)
    call alloc ('surface', 'st', st, npang)
    call alloc ('surface', 'cp', cp, npang)
    call alloc ('surface', 'sp', sp, npang)
    call alloc ('surface', 'angw', angw, npang)
    call lebgrid (ct,st,cp,sp,angw,npang)
    call allocate_space_for_surface (npang,minter)

    if (verbose) write (uout,'(1x,a)') string('# Computing surface')
    call cpu_time (time1)
    !$ time1 = omp_get_wtime()
    !$omp parallel default(none) &
    !$omp private(j,nsurf,rsurf) &
    !$omp shared(npang,ct,st,cp,sp,rlimsurf,nlimsurf)
    !$omp do schedule(dynamic)
    do j = 1,npang
      call surf (ct(j),st(j),cp(j),sp(j),rsurf,nsurf)
      do k = 1,nsurf
        rlimsurf(j,k) = rsurf(k,2)
      end do
      nlimsurf(j) = nsurf
    end do
    !$omp end do nowait
    !$omp end parallel
    call cpu_time (time2)
    !$ time2 = omp_get_wtime()
    if (verbose) then
      write (uout,'(1x,a,1x,f12.5)') string('# Elapsed seconds :'), time2-time1
    end if

    ltxt = 999
    lsu = 998
    open (ltxt,file=trim(filesurf)//"-txt"//d4)
    open (lsu,file=trim(filesurf)//"-"//d4,form='unformatted')
    write (ltxt,111) npang, inuc
    write (ltxt,333) (nlimsurf(j),j=1,npang)
    write (ltxt,090)
    write (lsu) npang, inuc
    write (lsu) (nlimsurf(j),j=1,npang)
    rmins = 1000_rp
    rmaxs = 0.0_rp
    do j = 1,npang
      nsurf = nlimsurf(j)
      rmins = min(rmins,rlimsurf(j,1))
      rmaxs = max(rmaxs,rlimsurf(j,nsurf))
      write (ltxt,222) ct(j),st(j),cp(j),sp(j),angw(j),(rlimsurf(j,k),k=1,nsurf)
      write (lsu) ct(j),st(j),cp(j),sp(j),angw(j),(rlimsurf(j,k),k=1,nsurf)
    end do
    write (ltxt,222) rmins,rmaxs
    write (lsu) rmins,rmaxs
    close (ltxt)
    close (lsu)
 
    ! End
    call free ('surface', 'ct', ct)
    call free ('surface', 'ct', st)
    call free ('surface', 'ct', cp)
    call free ('surface', 'ct', sp)
    call free ('surface', 'angw', angw)
    call deallocate_space_for_surface ()
    call deallocate_space_for_cp ()

090 format ('# cos(theta)',1x,'sin(theta)',1x,'cos(phi)',1x,'sin(phi)',2x,'weight')
111 format (2(1x,i5),' <--- (Angular points & Atom)')
333 format (20(1x,i2),4x,'(Surface intersections)')
222 format (15(1x,f12.6))

  end subroutine surface

  ! Determine the limit of the zero flux surface from xpoint along
  ! the theta phi direction. 
  subroutine surf (cot,sit,cop,sip,rsurf,nsurf)

    use mod_io, only: ferror, faterr
    use mod_wfn, only: ncent
    implicit none

    integer(kind=ip), intent(out) :: nsurf
    real(kind=rp), intent(in) :: cot, sit, cop, sip
    real(kind=rp), intent(out) :: rsurf(minter,2) 
 
    logical :: good
    integer(kind=ip) :: nintersec, nt
    integer(kind=ip) :: isurf(minter,2), i, ia, ib, im, ncount
    real(kind=rp) :: sintcosp, sintsinp, geofac
    real(kind=rp) :: xin(3), xfin(3), xpoint(3), xmed(3), ra, ract
    real(kind=rp) :: xdeltain(3), xsurf(0:minter,3), rb, rm, cost
 
    cost = cot
    sintcosp = sit*cop
    sintsinp = sit*sip
 
    if (ncent.eq.1) then
      nsurf = 1
      rsurf(1,1) = 0.0_rp
      rsurf(1,2) = rmaxsurf
      return
    end if
 
    ncount = 0_ip
    nintersec = 0_ip
    ia = inuc
    ra = 0.0_rp
    geofac = ((rmaxsurf-0.1_rp)/rprimer)**(1.0_rp/(ntrial-1))
    do i = 1,ntrial
      ract = rprimer*geofac**(i-1)
      xdeltain(1) = ract*sintcosp
      xdeltain(2) = ract*sintsinp
      xdeltain(3) = ract*cost    
      xpoint(:) = xnuc(:) + xdeltain(:)
      call odeint (xpoint,0.1_rp,1.0_rp)
      good = iscp(xpoint,ib)
      rb = ract
      if (ib.ne.ia .and. (ia.eq.inuc .or. ib.eq.inuc)) then
        if (ia.ne.inuc .or. ib.ne.-1) then
          nintersec = nintersec + 1_ip
          if (nintersec.gt.minter) then
            call ferror('surf', 'increase minter in mod_surf', faterr)
          end if
          xsurf(nintersec,1) = ra
          xsurf(nintersec,2) = rb
          isurf(nintersec,1) = ia
          isurf(nintersec,2) = ib
        end if
      end if
      ia = ib
      ra = rb
    end do ! looking for possible intersections
  
    ! We have now a consistent set of trial points. 
    ! Let us refine by bipartition, consistency check added. 
    ! No other nuclei basins can be found other than those explicit in isurf
    do i = 1,nintersec
      ia = isurf(i,1)
      ib = isurf(i,2)
      ra = xsurf(i,1)
      rb = xsurf(i,2)
      xin(1) = xnuc(1) + ra*sintcosp
      xin(2) = xnuc(2) + ra*sintsinp
      xin(3) = xnuc(3) + ra*cost
      xfin(1) = xnuc(1) + rb*sintcosp
      xfin(2) = xnuc(2) + rb*sintsinp
      xfin(3) = xnuc(3) + rb*cost
      do while (abs(ra-rb).gt.epssurf)
        ! Mean Value 
        xmed(:) = 0.5_rp*(xfin(:)+xin(:))    
        rm = (ra+rb)*0.5_rp
        xpoint(:) = xmed(:)
        call odeint (xpoint,0.1_rp,1.0_rp)
        good = iscp(xpoint,im)
        ! bipartition
        if (im.eq.ia) then
          xin(:) = xmed(:)
          ra = rm
        else if (im.eq.ib) then
          xfin(:) = xmed(:)
          rb = rm
        else
          if (ia.eq.inuc) then
            xfin(:) = xmed(:)
            rb = rm
          else
            xin(:) = xmed(:)
            ra = rm
          end if
        end if
      end do
      ! found. Mean value
      xpoint(:) = 0.5_rp*(xfin(:)+xin(:))    
      xsurf(i,3) = (ra+rb)*0.5_rp
    end do

    ! organize pairs
    nsurf = 0_ip
    xsurf(0,3) = 0.0_rp
    ia = inuc
    nsurf = nintersec
    do i = 1,nsurf
      rsurf(i,2) = xsurf(i,3)
    end do
    nt = mod(nintersec,2)
    if (nt.eq.0) then
      nsurf = nsurf + 1_ip
      rsurf(nsurf,2) = rmaxsurf
    end if
 
  end subroutine

  subroutine odeint (ystart,h1,iup)
 
    use mod_io, only: faterr, ferror, string
    use mod_param, only: vsmall
    use mod_fields, only: pointr1
    use mod_steeper, only: rkqs
    implicit none

    ! Parameters
    real(kind=rp), parameter :: hmin = 0.0_rp ! minimal step size

    ! Arguments
    real(kind=rp), intent(inout) :: ystart(3)
    real(kind=rp), intent(in) :: h1
    real(kind=rp), intent(in) :: iup

    ! Local vars
    integer(kind=ip) :: nstp, nuc
    real(kind=rp) :: x1 ! intial point
    real(kind=rp) :: x2 ! final point
    real(kind=rp) :: p(3) ! initial solution point
    real(kind=rp) :: h ! initial step size
    real(kind=rp) :: x ! update point
    real(kind=rp) :: hnext ! next steep size
    real(kind=rp) :: dydx(3), y(3), yscal(3)
    real(kind=rp) :: a1, a2, a3
    real(kind=rp) :: rho, grad(3), gradmod
 
    p(1) = ystart(1)
    p(2) = ystart(2)
    p(3) = ystart(3)
    call pointr1 (p,rho,grad,gradmod)
    if (gradmod.lt.epsgrad .and. rho.lt.epsrho) then
      return
    end if
 
    x1 = 0.0_rp
    x2 = 1d40*iup
    x = x1 ! initial point
    !h = sign(h1,x2-x1) ! initial steep size
    h = min(h1,x2-x1) ! initial steep size
    y(:) = ystart(:) ! initial point 
 
    do nstp = 1,maxstp
      call pointr1 (y,rho,grad,gradmod)
      dydx(:) = grad(:)
      yscal(:) = max(abs(y(:))+abs(h*dydx(:))+vsmall,epssurf)
      if ((x+h-x2)*(x+h-x1).gt.0.0_rp) h = x2 - x
      call rkqs (y,dydx,x,h,epssurf,yscal,hnext)
      if ((x-x2)*(x2-x1).ge.0.0_rp .or. iscp(y,nuc)) then
        ystart(:) = y(:)
        return
      end if
      if (abs(hnext).lt.hmin) then
        call ferror ('odeint', 'odeint stepsize small than minimum', faterr)
      end if
      if (nstp.eq.maxstp) then
        call ferror ('odeint', 'odeint reached maxstp', faterr)
      end if 
      h = hnext
    end do

    ! Test if the point is far from RMAXSURF from current atom. 
    a1 = y(1) - xnuc(1)
    a2 = y(2) - xnuc(2)
    a3 = y(3) - xnuc(3)
    if ((a1*a1+a2*a2+a3*a3).ge.5d0*5d0) then
      return
    else
      call ferror ('odeint', 'Non nuclear maximum at : ' &
                              //string(y(1),'e')//' '    &  
                              //string(y(2),'e')//' '    &  
                              //string(y(3),'e'), faterr) 
    end if
 
  end subroutine 


  subroutine findnuc ()

    use mod_fields, only: pointr1
    use mod_io, only: faterr, ferror, string
    use mod_wfn, only: xyz, charge, ncent
    implicit none

    integer(kind=ip) :: i
    real(kind=rp) :: rho, grad(3), gradmod, p(3)
    logical :: inf

    do i = 1,ncent
      p(:) = xyz(i,:)
      call gradrho (p,0.05_rp,1,inf)
      call pointr1 (p,rho,grad,gradmod)
      if (gradmod.gt.1d-4) then
        if (charge(i).gt.2.0_rp) then
          xyzrho(i,:) = xyz(i,:)
        else
          call ferror('findnuc', 'failed finding nucleus '//string(i), faterr)
        end if
      else
        xyzrho(i,:) = xyz(i,:)
      end if
    end do
 
  end subroutine

  ! Integration of a trajectory in the vector field of the electron density.
  ! xpoint() .... starting point of the trajectory
  ! iup ......... +1 tells the routine to climb up in the gradient field
  !               (usually going to end in a nucleus).
  !               -1 tells the routine to go down in the field.
  !
  ! Output data:
  ! xpoint() .... end point of the trajectory 
  subroutine gradrho (xpoint,hi,iup,inf)

    use mod_param, only: vsmall
    use mod_fields, only: pointr1
    implicit none
 
    logical, intent(inout) :: inf
    integer(kind=ip), intent(in) :: iup
    real(kind=rp), intent(in) :: hi
    real(kind=rp), intent(inout) :: xpoint(3)
 
    real(kind=rp) :: gradmod, rho, grad(3)
    integer(kind=ip) :: ier, niter, npoints, i
    real(kind=rp) :: xtemp(3), grdt(3), h0, escalar, grdmodule

    h0 = hi
    npoints = 1_ip
    inf = .false.
    call pointr1 (xpoint,rho,grad,gradmod)
    if (gradmod.lt.epsgrad .and. rho.lt.epsrho) then
      inf = .true.
      return
    end if
    grdt(1) = grad(1)/(gradmod+vsmall)
    grdt(2) = grad(2)/(gradmod+vsmall)
    grdt(3) = grad(3)/(gradmod+vsmall)
    grdmodule = gradmod
 
    escalar = 1_ip
    niter = 1_ip
    do while (grdmodule.gt.epsgrad .and. niter.lt.maxstp)
      niter = niter + 1_ip
      ier = 1_ip
      do while (ier.ne.0)
        xtemp(1) = xpoint(1) + h0*iup*grdt(1)
        xtemp(2) = xpoint(2) + h0*iup*grdt(2)
        xtemp(3) = xpoint(3) + h0*iup*grdt(3)
        call pointr1 (xtemp,rho,grad,gradmod)
        escalar = 0.0_rp
        do i = 1,3
          escalar = escalar + grdt(i)*(grad(i)/(gradmod+vsmall))
        end do
        ! It should'nt happen that h0 goes to zero, except if there
        ! are problems with the gradient, for instance the gradient
        ! has discontinuities or large rounding errors. Anyway, we
        ! introduce a safety check to avoid nasty infinite loops if
        ! h0 = 0.
        if (escalar.lt.0.707_rp) then
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
            grdt(i) = grad(i)/(gradmod+vsmall)
          enddo
          grdmodule = gradmod
        end if
      end do
    end do
 
  end subroutine

  logical function iscp (p,nuc)
 
    use mod_fields, only: pointr1
    use mod_wfn, only: ncent
    implicit none

    integer(kind=ip), intent(out) :: nuc
    real(kind=rp), intent(in) :: p(3)

    real(kind=rp) :: gradmod, grad(3), rho
    real(kind=rp) :: x(3)
    integer(kind=ip) :: i

    iscp = .false.
    nuc = 0
    x(:) = p(:)
    call pointr1 (x,rho,grad,gradmod)
    do i = 1,ncent
      if (abs(p(1)-xyzrho(i,1)).lt.epsiscp .and. &
          abs(p(2)-xyzrho(i,2)).lt.epsiscp .and. &
          abs(p(3)-xyzrho(i,3)).lt.epsiscp) then
        iscp = .true.
        nuc = i
      end if
    end do

    if (gradmod.le.epsgrad) then
      iscp = .true.
      if (rho.le.epsrho) nuc = -1
    end if

  end function

  subroutine optssurf(var,rval,ival)

    use iso_fortran_env, only: uout=>output_unit
    use mod_io, only: equal, faterr, ferror, string
    use mod_param, only: verbose
    implicit none

    character(len=*), intent(in) :: var
    real(kind=rp), optional :: rval
    integer(kind=ip), optional :: ival

    if (equal(var,'steeper')) then
      steeper = abs(ival)
      if (verbose) then
        write (uout,'(1x,a,1x,i0)') string('# *** Variable steeper changed to :'), steeper
      end if
    else if (equal(var,'ntrial')) then
      ntrial = abs(ival)
      if (verbose) then
        write (uout,'(1x,a,1x,i0)') string('# *** Variable ntrial changed to :'), ntrial
      end if
    else if (equal(var,'rmaxsurf')) then
      rmaxsurf = abs(rval)
      if (verbose) then
        write (uout,'(1x,a,1x,e13.6)') string('# *** Variable rmaxsurf changed to :'), rmaxsurf
      end if
    else if (equal(var,'epssurf')) then
      epssurf = abs(rval)
      if (verbose) then
        write (uout,'(1x,a,1x,e13.6)') string('# *** Variable epssurf changed to :'), epssurf
      end if
    else if (equal(var,'epsiscp')) then
      epsiscp = abs(rval)
      if (verbose) then
        write (uout,'(1x,a,1x,e13.6)') string('# *** Variable epsiscp changed to :'), epsiscp
      end if
    else if (equal(var,'rprimer')) then
      rprimer = abs(rval)
      if (verbose) then
        write (uout,'(1x,a,1x,e13.6)') string('# *** Variable rprimer changed to :'), rprimer
      end if
    else
      call ferror ('optssurf', 'unknown option', faterr)
    end if

  end subroutine optssurf

  subroutine init_surf ()

    implicit none

    steeper = 1_ip
    ntrial = 11_ip
    rprimer = 0.4_rp
    inuc = 0_ip
    epssurf = 1d-5 
    epsiscp = 0.08_rp
    nangleb(1) = 1
    nangleb(2) = 434
    nangleb(3) = 0
    nangleb(4) = 0
    rmaxsurf = 10.0_rp

  end subroutine

  subroutine allocate_space_for_cp (ncent)
  
    use mod_memory, only: alloc
    implicit none
    integer(kind=ip), intent(in) :: ncent

    call alloc ('mod_surf', 'xyzrho', xyzrho, ncent, 3)
 
  end subroutine allocate_space_for_cp

  subroutine deallocate_space_for_cp ()

    use mod_memory, only: free
    implicit none

    call free ('mod_surf', 'xyzrho', xyzrho)

  end subroutine deallocate_space_for_cp

  subroutine allocate_space_for_surface (nangular,ncutsurf)

    use mod_memory, only: alloc
    implicit none
    integer(kind=ip), intent(in) :: nangular, ncutsurf

    call alloc ('mod_surf', 'rlimsurf', rlimsurf, nangular, ncutsurf)
    call alloc ('mod_surf', 'nlimsurf', nlimsurf, nangular)

  end subroutine allocate_space_for_surface

  subroutine deallocate_space_for_surface ()

    use mod_memory, only: free
    implicit none

    call free ('mod_surf', 'rlimsurf', rlimsurf)
    call free ('mod_surf', 'nlimsurf', nlimsurf)

  end subroutine deallocate_space_for_surface

  subroutine info_surf ()

    use iso_fortran_env, only: uout=>output_unit
    use mod_io, only: string, mline
    use mod_wfn, only: ncent, atnam, charge
    implicit none

    character(len=mline), dimension(5) :: rqudstr
    character(len=mline), dimension(5) :: ssteeper
    integer(kind=ip) :: i

    rqudstr(1) = 'Gauss-Legendre'
    ssteeper(1) = 'Runge-Kutta-Cash-Karp'

    atnam = adjustl(atnam)
    write (uout,'(1x,a,1x)') string('# Positions of nuclear maxima')
    do i = 1,ncent
      write (uout,'(1x,a,1x,i0,1x,a,1x,f4.1,1x,3f12.6)') '#', i, &
                    atnam(i)(1:2), charge(i), xyzrho(i,:)
    end do
    write (uout,'(1x,a,1x,i0)') string('# Computing SURFACE for atom :'), inuc
    write (uout,'(1x,a,1x,3f12.6)') string('# Nuclear center for atom :'), xnuc(:)
    write (uout,'(1x,a,1x,e13.6)') string('# Rmaxsurf :'), rmaxsurf 
    write (uout,'(1x,a,1x,a)') string('# Steeper :'), string(ssteeper(steeper))
    write (uout,'(1x,a,1x,e13.6)') string('# Surface precision :'), epssurf
    write (uout,'(1x,a,1x,e13.6)') string('# EPSISCP parameter :'), epsiscp
    write (uout,'(1x,a,1x,i0)') string('# Ntrial :'), ntrial
    write (uout,'(1x,a,1x,e13.6)') string('# Rprimer :'), rprimer
    write (uout,'(1x,a)') string('# Angular Quadrature')
    if (nangleb(1).eq.1) then
      write (uout,'(1x,a,1x,i0,1x,a,1x,i0)') &
      string('# Atom'), inuc, string('lebedev points'), nangleb(2)
    else if (nangleb(1).eq.0) then
      write (uout,'(1x,a)') string('# Phi quadrature is always trapezoidal')
      write (uout,'(1x,a,1x,i0,1x,a,1x,i0,1x,i0,1x,a)') &
      string('# Atom'), inuc, string('(ntheta,nphi,iqudt'), nangleb(2), nangleb(3), &
      string(rqudstr(nangleb(4)))
    end if

  end subroutine info_surf

end module mod_surf
