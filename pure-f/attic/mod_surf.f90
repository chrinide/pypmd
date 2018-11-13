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

  integer(kind=ip), parameter :: minter = 10 
  integer(kind=ip), parameter :: maxtrial = 13

  integer(kind=ip), public :: inuc
  real(kind=rp), public :: xnuc(3)
  real(kind=rp), allocatable, dimension(:,:), public :: rlimsurf
  integer(kind=ip), allocatable, dimension(:), public :: nlimsurf

  ! nuclear critical points
  real(kind=rp), allocatable, dimension(:,:), public :: xyzrho

  ! options
  integer(kind=ip), dimension(4), public :: nangleb
  ! Precision in interatomic surfaces. 
  real(kind=rp), public :: rmaxsurf
  ! Precision in interatomic surfaces. 
  real(kind=rp), public :: epsilon
  ! If a non nuclear maximum is found whose distance to any nucleus 
  ! is larger than EPSISCP, promolden STOPS. This control of non nuclear
  ! maxima is performed in the 'iscp.f' routine. When Hydrogen atoms
  ! are involved it is very convenient to change the default value
  ! for EPSISCP (0.08) to a smaller value by using the 'EPSISCP value'
  ! order in the input file.
  real(kind=rp), public :: epsiscp
  ! NTRIAL  = number of sub-intervals in the search of the surface.
  ! RPRIMER = First point in the in the search of the surface.
  ! GEOFAC  = ((rmaxsurf-0.1)/rprimer)**(1.0/(ntrial-1))
  integer(kind=ip), public :: ntrial
  real(kind=rp), public :: rprimer
  real(kind=rp), public :: geofac
  integer(kind=ip), public :: steeper

  public :: optssurf, init_surf, surface

contains

  subroutine surface (nuc)

    use iso_fortran_env, only: uout=>output_unit
    !$ use omp_lib, only: omp_get_wtime
    use mod_io, only: mline, fourchar, flush_unit, string
    use mod_wfn, only: ncent, filename
    use mod_param, only: verbose, pi
    use mod_memory, only: free, alloc
    implicit none
 
    integer(kind=ip), intent(in) :: nuc
 
    character(len=4) :: d4
    character(len=mline) :: filesurf
    real(kind=rp) :: time1, time2, phi, delphi, rsurf(minter,2)
    real(kind=rp) :: rmins, rmaxs, thang 
    integer(kind=ip) :: lsu, ltxt, j, k, npang, nsurf, i, iqudt
    integer(kind=ip) :: ntheta, nphi, it, ipp
    real(kind=rp), allocatable, dimension(:) :: ct, st, cp, sp, angw
    real(kind=rp), allocatable, dimension(:) :: tp, tw

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
    ltxt = 999
    lsu = 998
    filesurf = trim(filename)//".surf"
    call allocate_space_for_cp (ncent)
    call findnuc ()
    inuc = nuc
    xnuc(:) = xyzrho(inuc,:)
    if (verbose) call info_surf ()
      
    ! Begin
    call flush_unit (uout)
    d4 = fourchar(inuc)
    if (nangleb(1).eq.1_ip) then
      npang = nangleb(2)
      call good_lebedev (npang)
      call alloc ('surface', 'ct', ct, npang)
      call alloc ('surface', 'st', st, npang)
      call alloc ('surface', 'cp', cp, npang)
      call alloc ('surface', 'sp', sp, npang)
      call alloc ('surface', 'angw', angw, npang)
      call lebgrid (ct,st,cp,sp,angw,npang)
    else
      ntheta = nangleb(2)
      nphi = nangleb(3)
      iqudt = nangleb(4)
      npang = ntheta*nphi  
      call alloc ('surface', 'tp', tp, ntheta)
      call alloc ('surface', 'tw', tw, ntheta)
      call alloc ('surface', 'ct', ct, npang)
      call alloc ('surface', 'st', st, npang)
      call alloc ('surface', 'cp', cp, npang)
      call alloc ('surface', 'sp', sp, npang)
      call alloc ('surface', 'angw', angw, npang)
      call weightheta (iqudt,tp,tw,ntheta)
      delphi = 2.0_rp*pi/nphi
      i = 0_ip
      do ipp = 0,nphi-1
        phi = ipp*delphi
        do it = 1,ntheta
          i = i + 1_ip
          thang = tp(it)
          ct(i) = thang
          st(i) = sqrt(1.0_rp-thang*thang)
          cp(i) = cos(phi)
          sp(i) = sin(phi)
          angw(i) = tw(it)*delphi
        end do
      end do
      call free ('surface', 'tp', tp)
      call free ('surface', 'tw', tw)
    end if

    call deallocate_space_for_surface ()
    call free ('surface', 'ct', ct)
    call free ('surface', 'ct', st)
    call free ('surface', 'ct', cp)
    call free ('surface', 'ct', sp)
    call free ('surface', 'angw', angw)

    ! Deallocate all arrays
    call deallocate_space_for_cp ()

109 format (9x,'cos(theta)',13x,'sin(theta)',13x,'cos(phi)',15x,'sin(phi)',15x,'weight')
111 format (2(1x,i0),' <--- (Angular points & Atom)')
333 format (20(1x,i0),4x,'(Surface intersections)')
222 format (15(1x,f22.15))

  end subroutine surface

  ! Integration of a trajectory in the vector field of the electron density.
  ! xpoint() .... starting point of the trajectory
  ! iup ......... +1 tells the routine to climb up in the gradient field
  !               (usually going to end in a nucleus).
  !               -1 tells the routine to go down in the field.
  !
  ! Output data:
  ! xpoint() .... end point of the trajectory 
  subroutine gradrho (xpoint,hi,iup,inf)

    use mod_fields, only: pointr1
    implicit none
 
    real(kind=rp), parameter :: epsg = 1d-10
    real(kind=rp), parameter :: eps = 1d-8
    real(kind=rp), parameter :: epsg1 = 1d-10
    real(kind=rp), parameter :: hminimal = 1d-40
    integer(kind=ip), parameter :: mstep = 500
 
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
    if (gradmod.lt.epsg .and. rho.lt.epsg1) then
      inf = .true.
      return
    end if
    grdt(1) = grad(1)/gradmod
    grdt(2) = grad(2)/gradmod
    grdt(3) = grad(3)/gradmod
    grdmodule = gradmod
 
    escalar = 1_ip
    niter = 1_ip
    do while (grdmodule.gt.eps .and. niter.lt.mstep)
      niter = niter + 1_ip
      ier = 1_ip
      do while (ier.ne.0)
        xtemp(1) = xpoint(1) + h0*iup*grdt(1)
        xtemp(2) = xpoint(2) + h0*iup*grdt(2)
        xtemp(3) = xpoint(3) + h0*iup*grdt(3)
        call pointr1 (xtemp,rho,grad,gradmod)
        escalar = 0.0_rp
        do i = 1,3
          escalar = escalar + grdt(i)*(grad(i)/(gradmod+1d-80))
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
            grdt(i) = grad(i)/gradmod
          enddo
          grdmodule = gradmod
        end if
      end do
    end do
 
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

    if (gradmod.le.1d-10) then
      iscp = .true.
      if (rho.le.1d-10) nuc = -1
    end if

  end function

  ! Abscissas and weights for the theta quadrature. 
  subroutine weightheta (iq,x,w,n)

    use mod_io, only: ferror, faterr
    implicit none

    ! Arguments
    integer(kind=ip), intent(in) :: iq, n
    real(kind=rp), intent(out) :: x(n), w(n)

    if (iq.eq.1) then ! Gauss-Legendre quadrature 
      call gauleg (-1.0_rp,+1.0_rp,x,w,n)
    else if (iq.eq.2) then ! Clenshaw-Curtis quadrature
      call genclcu (-1.0_rp,+1.0_rp,x,w,n)
    else if (iq.eq.3) then ! Gauss-Chebychev of first kind 
      call chebyshev1 (x,w,n)
    else if (iq.eq.4) then ! Gauss-Chebychev of second kind 
      call chebyshev2 (x,w,n)
    else if (iq.eq.5) then ! Perez-Jorda (Gauss-Chebychev, 2nd kind) quadrature. 
      call pjt (x,w,n)
    else
      call ferror ('weightheta', 'bad theta quadrature', faterr)
    end if
 
  end subroutine

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
    else if (equal(var,'epsilon')) then
      epsilon = abs(rval)
      if (verbose) then
        write (uout,'(1x,a,1x,e13.6)') string('# *** Variable epsilon changed to :'), epsilon
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
    epsilon = 1d-5 
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
    rqudstr(2) = 'Clenshaw-Curtis'
    rqudstr(3) = 'Gauss-Chebychev 1st kind'
    rqudstr(4) = 'Gauss-Chebychev 2nd kind'
    rqudstr(5) = 'Perez-Jorda (Gauss-Chebychev) 2nd kind'

    ssteeper(1) = 'Runge-Kutta-Cash-Karp'
    ssteeper(2) = 'Calvo-Montijano-Randez'
    ssteeper(3) = 'Dormand-Prince method'

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
    write (uout,'(1x,a,1x,e13.6)') string('# Surface precision :'), epsilon
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
