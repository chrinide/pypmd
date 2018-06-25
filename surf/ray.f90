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
! Determine the limit of the zero flux surface from xpoint along
! the theta phi direction. 
!
subroutine ray (cot,sit,cop,sip,rsurf,nsurf)

  use mod_prec, only: rp, ip
  use mod_io, only: ferror, faterr
  use mod_wfn, only: ncent
  use mod_odeint, only: odeint
  use mod_surf, only: geofac, rprimer, ntrial, epsilon, steeper, &
                      minter, inuc, xnuc, rmaxsurf 
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
 
  interface
    logical function  iscp (p,nuc)
      import :: rp, ip
      integer(kind=ip), intent(out) :: nuc
      real(kind=rp), intent(in) :: p(3)
    end function  iscp
  end interface
 
  cost = cot
  sintcosp = sit*cop
  sintsinp = sit*sip
 
  if (ncent.eq.1) then
    nsurf = 1
    rsurf(1,1) = 0.0_rp
    rsurf(1,2) = rmaxsurf
    return
  end if
 
  inf = .false.
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
    call odeint (xpoint,0.1_rp,1.0_rp,inf,epsilon,xnuc,steeper)
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
  end do        ! looking for intersections

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
    do while (abs(ra-rb).gt.epsilon)
      ! Mean Value 
      xmed(:) = half*(xfin(:)+xin(:))    
      rm = (ra+rb)*half
      xpoint(:) = xmed(:)
      call odeint (xpoint,0.1_rp,1.0_rp,inf,epsilon,xnuc,steeper)
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
    xpoint(:) = half*(xfin(:)+xin(:))    
    xsurf(i,3) = (ra+rb)*half
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
