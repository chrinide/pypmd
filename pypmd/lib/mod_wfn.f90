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
module mod_wfn
  
  use mod_prec, only: ip, rp
  implicit none
  public

  integer(kind=ip), parameter :: mgrp = 500_ip
  integer(kind=ip), parameter :: ngtoh = 21_ip
  integer(kind=ip), parameter, private :: maxtype = 56_ip

  ! nlm keeps the nlm values of x^n y^l z^m gaussian
  integer(kind=ip) :: nlm(maxtype,3)
  integer(kind=ip), private :: nprimsold
  integer(kind=ip) :: nprims_
  integer(kind=ip) :: nmo_
  integer(kind=ip) :: ncent_
  real(kind=rp) :: epscuttz_
  real(kind=rp) :: rmaxatom_
  integer(kind=ip) :: maxgrp_
  integer(kind=ip) :: numshells_

  real(kind=rp), allocatable, dimension(:,:) :: coeff_
  real(kind=rp), allocatable, dimension(:,:) :: coords_
  real(kind=rp), allocatable, dimension(:) :: oexp_
  real(kind=rp), allocatable, dimension(:) :: occ_
  real(kind=rp), allocatable, dimension(:) :: charges_
  integer(kind=ip), allocatable, dimension(:) :: icen_
  integer(kind=ip), allocatable, dimension(:) :: ityp_

  integer(kind=ip), allocatable, dimension(:) :: npc_
  integer(kind=ip), allocatable, dimension(:) :: ngroup_
  integer(kind=ip), allocatable, dimension(:,:) :: icenat_
  integer(kind=ip), allocatable, dimension(:,:) :: nzexp_
  integer(kind=ip), allocatable, dimension(:,:,:) :: nuexp_
  real(kind=rp), allocatable, dimension(:,:) :: rcutte_

  integer(kind=ip), allocatable, dimension(:,:):: iden_
  integer(kind=ip), allocatable, dimension(:) :: nden_
  integer(kind=ip), allocatable, dimension(:,:) :: ishell_
  integer(kind=ip), allocatable, dimension(:) :: nshell_
  integer(kind=ip), allocatable, dimension(:,:) :: atcenter_

  private :: etijcalc

contains

  subroutine setup_wfn ()

    use mod_memory, only: alloc, free
    use mod_io, only: ferror, faterr
    implicit none

    real(kind=rp) :: alph
    integer(kind=ip) :: npa, i, ic, icentro, inda, isud
    integer(kind=ip) :: isuk, itd, itip, itk, j, k, m, npcant
    real(kind=rp), allocatable, dimension(:) :: oexpa
    real(kind=rp), allocatable, dimension(:,:) :: coefa
    integer(kind=ip), allocatable, dimension(:) :: icena
    integer(kind=ip), allocatable, dimension(:) :: itypa
    
    ! Init
    call alloc ('mod_wfn/setupwfn', 'oexpa', oexpa, nprims_)
    call alloc ('mod_wfn/setupwfn', 'coefa', coefa, nmo_, nprims_)
    call alloc ('mod_wfn/setupwfn', 'icena', icena, nprims_)
    call alloc ('mod_wfn/setupwfn', 'itypa', itypa, nprims_)
    npa = 1_ip
    icena(1) = icen_(1)
    oexpa(1) = oexp_(1)
    itypa(1) = ityp_(1)
    coefa(1:nmo_,1) = coeff_(1:nmo_,1)
    cyclej: do j= 2,nprims_
      do m = 1,npa
        if (icen_(j).eq.icena(m) .and. ityp_(j).eq.itypa(m) .and. abs(oexp_(j)-oexpa(m)).le.1d-10) then
          coefa(1:nmo_,m) = coefa(1:nmo_,m) + coeff_(1:nmo_,j)
          cycle cyclej
        end if
      end do
      npa = npa + 1_ip
      icena(npa) = icen_(j)
      oexpa(npa) = oexp_(j)
      itypa(npa) = ityp_(j)
      coefa(1:nmo_,npa) = coeff_(1:nmo_,j)
    end do cyclej

    ! Recompute the original variables
    nprimsold = nprims_
    nprims_ = npa
    do j = 1,nprims_
      icen_(j) = icena(j)
      oexp_(j) = oexpa(j)
      ityp_(j) = itypa(j)
      coeff_(1:nmo_,j) = coefa(1:nmo_,j)
    end do

    ! Determine primitives corresponding to each center.
    do ic = 1,ncent_
      npc_(ic) = 0_ip
    end do
    do j = 1,nprims_
      ic = icen_(j)
      npc_(ic) = npc_(ic) + 1_ip
      inda = npc_(ic)
      icenat_(inda,ic) = j
    end do

    ! Classify primitives in each center by types and similar exponents.
    do ic = 1,ncent_
      scyclej: do j = 1,npc_(ic)
        k = icenat_(j,ic)
        itk = ityp_(k)
        isuk = nlm(itk,1) + nlm(itk,2) + nlm(itk,3)
        if (j.eq.1) then
          ngroup_(ic) = 1_ip
          nzexp_(ic,1) = 1_ip
          nuexp_(ic,1,1) = k
        else
          do m = 1,ngroup_(ic)
            inda = nuexp_(ic,m,1)
            itd = ityp_(inda)
            isud = nlm(itd,1) + nlm(itd,2) + nlm(itd,3)
            if (abs(oexp_(k)-oexp_(inda)).lt.1d-8) then
              if (itk.eq.1.and.itd.eq.1) then
                call ferror ('mod_wfn/setupwfn', 'two s primitives with equal exponents', faterr)
              else
                if (isuk.eq.isud) then
                  nzexp_(ic,m) = nzexp_(ic,m) + 1_ip
                  nuexp_(ic,m,nzexp_(ic,m)) = k
                  cycle scyclej
                end if
              end if
            end if
          end do
          ngroup_(ic) = ngroup_(ic) + 1_ip
          if (ngroup_(ic).gt.mgrp) then
            call ferror ('mod_wfn/setupwfn', 'increase mgrp in mod_wfn file', faterr)
          end if
          nzexp_(ic,ngroup_(ic)) = 1_ip
          nuexp_(ic,ngroup_(ic),1) = k
        end if   
      end do scyclej
    end do
  
    ! Reconstruct the values of ityp(),icen(),oexp(), and coef()
    i = 0_ip
    do ic = 1,ncent_
      do m = 1,ngroup_(ic)
        do k = 1,nzexp_(ic,m)
          j = nuexp_(ic,m,k)
          alph = oexp_(j)
          itip = ityp_(j)
          icentro = icen_(j)
          i = i + 1_ip
          itypa(i) = itip
          oexpa(i) = alph
          icena(i) = icentro
          coefa(1:nmo_,i) = coeff_(1:nmo_,j)
        enddo
      end do
    end do
  
    call free ('mod_wfn/setupwfn', 'coeff', coeff_)
    call alloc ('mod_wfn/setupwfn', 'coeff', coeff_, nmo_, nprims_)
    do i = 1,nprims_
      ityp_(i) = itypa(i)
      oexp_(i) = oexpa(i)
      icen_(i) = icena(i)
      coeff_(1:nmo_,i) = coefa(1:nmo_,i)
    end do
   
    npcant = 0_ip
    do ic = 1,ncent_
      do k = 1,npc_(ic)
        icenat_(k,ic) = k + npcant
      end do
      npcant = npcant + npc_(ic)
    end do

    ! Reconstruct the values of nuexp(). Now, they are ordered.
    ! Determine also which is the maximum value of ngroup(ic)
    ! Determine the total number of shells.
    i = 0_ip
    numshells_ = 0_ip
    maxgrp_ = 0_ip
    do ic = 1,ncent_
      numshells_ = numshells_ + ngroup_(ic)
      if (ngroup_(ic).gt.maxgrp_) maxgrp_ = ngroup_(ic)
      do m = 1,ngroup_(ic)
        do k = 1,nzexp_(ic,m)
          i = i + 1_ip
          nuexp_(ic,m,k) = i
        end do
      end do
    end do

    ! Deallocate arrays.
    call free ('mod_wfn/setupwfn', 'oexpa', oexpa)
    call free ('mod_wfn/setupwfn', 'coefa', coefa)
    call free ('mod_wfn/setupwfn', 'icena', icena)
    call free ('mod_wfn/setupwfn', 'itypa', itypa)

  end subroutine

  subroutine filtergto ()

    implicit none

    logical :: okcen
    integer(kind=ip) :: j, ic, i, jc, lsum, m
    real(kind=rp) :: dis, x1, x2, y1, y2, z1, z2, xmax, zz  
    real(kind=rp), dimension(ncent_,ncent_) :: rint

    ! Evaluate internuclear distances
    do i = 1,ncent_
      x1 = coords_(1,i)
      y1 = coords_(2,i)
      z1 = coords_(3,i)
      do j = 1,i
        x2 = coords_(1,j)
        y2 = coords_(2,j)
        z2 = coords_(3,j)
        dis = sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
        rint(i,j) = dis
        rint(j,i) = dis
      end do
    end do

    ! Maximum distance at which it is necessary to compute a shell.
    do ic = 1,ncent_
      do m = 1,ngroup_(ic)
        i = nuexp_(ic,m,1)
        lsum = nlm(ityp_(i),1)+nlm(ityp_(i),2)+nlm(ityp_(i),3)
        zz = oexp_(i)
        x1 = 0.1_rp
        do 
          if (x1**lsum*exp(-zz*x1*x1).le.abs(epscuttz_)) exit
          x1 = x1 + 0.1_rp
        end do
        rcutte_(ic,m) = x1
      end do
    end do

    nden_(1:ncent_) = 0_ip
    nshell_(1:ncent_) = 0_ip
    do ic = 1,ncent_
      xmax = rmaxatom_
      do jc = 1,ncent_
        dis = rint(ic,jc)
        okcen = .false.
        ! The shell 'm' of center 'jc' contributes to the density, orbitals, 
        ! orbital products, etc, at any point inside the center 'ic' if the 
        ! following condition holds.
        do m = 1,ngroup_(jc)
          if (dis.lt.xmax+rcutte_(jc,m)) then
            nshell_(ic) = nshell_(ic) + 1_ip
            atcenter_(ic,nshell_(ic)) = jc
            ishell_(ic,nshell_(ic)) = m
            okcen = .true.
          end if
        end do
        if (okcen) then
          nden_(ic) = nden_(ic) + 1_ip
          iden_(ic,nden_(ic)) = jc
        end if
      end do
    end do

    do ic = 1,ncent_
      do m = 1,ngroup_(ic)
        rcutte_(ic,m) = rcutte_(ic,m)*rcutte_(ic,m)
      end do
    end do

  end subroutine

  subroutine init_wfn()
                
    implicit none

    !.p's
    nlm(2,1)=1      !px
    nlm(2,2)=0      !px
    nlm(2,3)=0      !px
    nlm(3,1)=0      !py
    nlm(3,2)=1      !py
    nlm(3,3)=0      !py
    nlm(4,1)=0      !pz
    nlm(4,2)=0      !pz
    nlm(4,3)=1      !pz
    !.d's
    nlm(5,1)=2      !xx
    nlm(6,2)=2      !yy
    nlm(7,3)=2      !zz
    nlm(8,1)=1      !xy
    nlm(8,2)=1
    nlm(9,1)=1      !xz
    nlm(9,3)=1
    nlm(10,2)=1     !yz
    nlm(10,3)=1
    !.f's
    nlm(11,1)=3     !xxx
    nlm(12,2)=3     !yyy
    nlm(13,3)=3     !zzz    
    nlm(14,1)=2     !xxy
    nlm(14,2)=1
    nlm(15,1)=2     !xxz
    nlm(15,3)=1
    nlm(16,2)=2     !yyz
    nlm(16,3)=1 
    nlm(17,1)=1     !xyy
    nlm(17,2)=2 
    nlm(18,1)=1     !xzz
    nlm(18,3)=2     
    nlm(19,2)=1     !yzz
    nlm(19,3)=2
    nlm(20,1)=1     !xyz
    nlm(20,2)=1 
    nlm(20,3)=1 
    !.g's
    nlm(21,1)=4     !xxxx
    nlm(22,2)=4     !yyyy
    nlm(23,3)=4     !zzzz
    nlm(24,1)=3     !xxxy
    nlm(24,2)=1
    nlm(25,1)=3     !xxxz
    nlm(25,3)=1 
    nlm(26,1)=1     !xyyy
    nlm(26,2)=3 
    nlm(27,2)=3     !yyyz
    nlm(27,3)=1 
    nlm(28,1)=1     !xzzz 
    nlm(28,3)=3 
    nlm(29,2)=1     !yzzz
    nlm(29,3)=3
    nlm(30,1)=2     !xxyy
    nlm(30,2)=2 
    nlm(31,1)=2     !xxzz
    nlm(31,3)=2
    nlm(32,2)=2     !yyzz 
    nlm(32,3)=2 
    nlm(33,1)=2     !xxyz 
    nlm(33,2)=1
    nlm(33,3)=1
    nlm(34,1)=1     !xyyz
    nlm(34,2)=2 
    nlm(34,3)=1
    nlm(35,1)=1     !xyzz
    nlm(35,2)=1
    nlm(35,3)=2
    !.h's
    nlm(36,1)=0
    nlm(36,2)=0
    nlm(36,3)=5

    nlm(37,1)=0
    nlm(37,2)=1
    nlm(37,3)=4

    nlm(38,1)=0
    nlm(38,2)=2
    nlm(38,3)=3

    nlm(39,1)=0
    nlm(39,2)=3
    nlm(39,3)=2

    nlm(40,1)=0
    nlm(40,2)=4
    nlm(40,3)=1

    nlm(41,1)=0
    nlm(41,2)=5
    nlm(41,3)=0

    nlm(42,1)=1
    nlm(42,2)=0
    nlm(42,3)=4

    nlm(43,1)=1
    nlm(43,2)=1
    nlm(43,3)=3

    nlm(44,1)=1
    nlm(44,2)=2
    nlm(44,3)=2

    nlm(45,1)=1
    nlm(45,2)=3
    nlm(45,3)=1

    nlm(46,1)=1
    nlm(46,2)=4
    nlm(46,3)=0

    nlm(47,1)=2
    nlm(47,2)=0
    nlm(47,3)=3
  
    nlm(48,1)=2
    nlm(48,2)=1
    nlm(48,3)=2
  
    nlm(49,1)=2
    nlm(49,2)=2
    nlm(49,3)=1

    nlm(50,1)=2
    nlm(50,2)=3
    nlm(50,3)=0

    nlm(51,1)=3
    nlm(51,2)=0
    nlm(51,3)=2

    nlm(52,1)=3
    nlm(52,2)=1
    nlm(52,3)=1

    nlm(53,1)=3
    nlm(53,2)=2
    nlm(53,3)=0

    nlm(54,1)=4
    nlm(54,2)=0
    nlm(54,3)=1

    nlm(55,1)=4
    nlm(55,2)=1
    nlm(55,3)=0

    nlm(56,1)=5
    nlm(56,2)=0
    nlm(56,3)=0
  
  end subroutine init_wfn

  subroutine allocate_space_for_wfn ()
  
    use mod_memory, only: alloc
    implicit none

    call alloc ('mod_wfn', 'coeff', coeff_, nmo_, nprims_)
    call alloc ('mod_wfn', 'coords', coords_, 3, ncent_)
    call alloc ('mod_wfn', 'oexp', oexp_, nprims_)
    call alloc ('mod_wfn', 'occ', occ_, nmo_)
    call alloc ('mod_wfn', 'charge', charges_, ncent_)
    call alloc ('mod_wfn', 'icen', icen_, nprims_)
    call alloc ('mod_wfn', 'ityp', ityp_, nprims_)

    call alloc ('mod_wfn', 'npc', npc_, ncent_)
    call alloc ('mod_wfn', 'ngroup', ngroup_, ncent_)
    call alloc ('mod_wfn', 'icenat', icenat_, nprims_, ncent_)
    call alloc ('mod_wfn', 'nzexp', nzexp_, ncent_, mgrp)
    call alloc ('mod_wfn', 'nuexp', nuexp_, ncent_, mgrp, ngtoH)

    call alloc ('mod_wfn', 'rcutte', rcutte_, ncent_, mgrp)
 
  end subroutine allocate_space_for_wfn

  subroutine deallocate_space_for_wfn ()

    use mod_memory, only: free
    implicit none

    call free ('mod_wfn', 'coeff', coeff_)
    call free ('mod_wfn', 'coords', coords_)
    call free ('mod_wfn', 'oexp', oexp_)
    call free ('mod_wfn', 'occ', occ_)
    call free ('mod_wfn', 'charge', charges_)
    call free ('mod_wfn', 'icen', icen_)
    call free ('mod_wfn', 'ityp', ityp_)

    call free ('mod_wfn', 'npc', npc_)
    call free ('mod_wfn', 'ngroup', ngroup_)
    call free ('mod_wfn', 'icenat', icenat_)
    call free ('mod_wfn', 'nzexp', nzexp_)
    call free ('mod_wfn', 'nuexp', nuexp_)

    call free ('mod_wfn', 'rcutte', rcutte_)

  end subroutine deallocate_space_for_wfn

  subroutine allocate_space_for_shells ()

    use mod_memory, only: alloc
    implicit none
  
    call alloc ('mod_wfn', 'iden', iden_, ncent_, ncent_)
    call alloc ('mod_wfn', 'nden', nden_, ncent_)
    call alloc ('mod_wfn', 'nshell', nshell_, ncent_)
    call alloc ('mod_wfn', 'ishell', ishell_, ncent_, numshells_)
    call alloc ('mod_wfn', 'atcenter', atcenter_, ncent_, numshells_)

  end subroutine allocate_space_for_shells
                                                                        
  subroutine deallocate_space_for_shells ()

    use mod_memory, only: free
    implicit none

    call free ('mod_wfn', 'iden', iden_)
    call free ('mod_wfn', 'nden', nden_)
    call free ('mod_wfn', 'nshell', nshell_)
    call free ('mod_wfn', 'ishell', ishell_)
    call free ('mod_wfn', 'atcenter', atcenter_)

  end subroutine deallocate_space_for_shells

  subroutine etijcalc (m,lamx,la,lb,ce,a,b,ax,bx)
 
    use mod_io, only: faterr, ferror

    implicit none
    integer(kind=ip) :: la,lb,lab,i,j,t,i1,j1,t1,m,lamx
    real(kind=rp) :: ce(-1:2*lamx,-1:lamx,-1:lamx,3)
    real(kind=rp) :: a,b,ax,bx,p,ab,pa,pb,tp

    if (la.lt.0) call ferror ('etijcalc', 'fatal error, la < 0', faterr)
    if (lb.lt.0) call ferror ('etijcalc', 'fatal error, lb < 0', faterr)
    if (la.gt.lamx) call ferror ('etijcalc', 'fatal error, la > lamx', faterr)
    if (lb.gt.lamx) call ferror ('etijcalc', 'fatal error, lb > lamx', faterr)
  
    lab = la + lb
    ce(-1:lab,-1:la,-1:lb,m) = 0.0_rp
    ce(0,0,0,m) = 1.0_rp
    if (lab.eq.0) return
    p  = a + b
    ab = ax - bx
    pa = -b*ab/p
    pb = +a*ab/p
    tp = 1.0_rp/(2.0_rp*p)
    do i = 0,la
      i1 = i-1
      do j = 0,lb
        j1 = j-1
        do t = 1,i+j
          t1 = t-1
          ce(t,i,j,m) = tp*(i*ce(t1,i1,j,m) + j*ce(t1,i,j1,m))/real(t,rp)
        end do
        if (i.lt.la) ce(0,i+1,j,m) = pa*ce(0,i,j,m) + ce(1,i,j,m)
        if (j.lt.lb) ce(0,i,j+1,m) = pb*ce(0,i,j,m) + ce(1,i,j,m)
      end do
    end do 
 
  end subroutine


end module mod_wfn
