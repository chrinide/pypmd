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
  use mod_io, only: mline
  implicit none
  public

  integer(kind=ip), parameter, private :: mgrp = 500_ip
  integer(kind=ip), parameter, private :: ngtoh = 21_ip
  integer(kind=ip), parameter, private :: maxtype = 56_ip

  ! nlm keeps the nlm values of x^n y^l z^m gaussian
  integer(kind=ip) :: nlm(maxtype,3)
  integer(kind=ip) :: nprims
  integer(kind=ip), private :: nprimsold
  integer(kind=ip) :: nmo
  integer(kind=ip) :: ncent
  integer(kind=ip) :: maxgrp
  integer(kind=ip) :: numshells

  real(kind=rp), allocatable, dimension(:,:) :: coeff
  integer(kind=ip), allocatable, dimension(:) :: npc
  integer(kind=ip), allocatable, dimension(:) :: ngroup
  integer(kind=ip), allocatable, dimension(:,:) :: icenat
  integer(kind=ip), allocatable, dimension(:,:) :: nzexp
  integer(kind=ip), allocatable, dimension(:,:,:) :: nuexp
  real(kind=rp), allocatable, dimension(:,:) :: xyz
  real(kind=rp), allocatable, dimension(:) :: oexp
  real(kind=rp), allocatable, dimension(:,:) :: rcutte
  real(kind=rp), allocatable, dimension(:,:) :: rint
  real(kind=rp), allocatable, dimension(:) :: occ
  real(kind=rp), allocatable, dimension(:) :: charge
  character(len=8), allocatable, dimension(:) :: atnam
  integer(kind=ip), allocatable, dimension(:) :: icen
  integer(kind=ip), allocatable, dimension(:) :: ityp
  integer(kind=ip), allocatable, dimension(:,:):: iden
  integer(kind=ip), allocatable, dimension(:) :: nden
  integer(kind=ip), allocatable, dimension(:,:) :: ishell
  integer(kind=ip), allocatable, dimension(:) :: nshell
  integer(kind=ip), allocatable, dimension(:,:) :: atcenter

  ! OPTIONS:
  real(kind=rp) :: epscuttz
  real(kind=rp) :: rmaxatom
  real(kind=rp) :: epsortho

  private :: etijcalc, gtogto
  
contains

  subroutine loadwfn(filedat)

    use mod_param, only: verbose, filename
    implicit none

    character(len=*), intent(in) :: filedat

    filename = filedat

    call rdwfn (filename)
    call setupwfn ()
    call filtergto ()
    if (verbose) call info_wfn ()
    !call isorthowfn ()

  end subroutine loadwfn

  subroutine rdwfn (wfnfile)

    use mod_memory, only: alloc, free
    use mod_io, only: ferror, faterr, udat, string, warning
    implicit none
 
    ! Arguments
    character(len=*), intent(in) :: wfnfile
 
    ! Local vars
    real(kind=rp) :: tmp
    integer(kind=ip) :: i, iwfn, j, k
    character(len=80) :: wfnttl
    character(len=4) :: mode
    character(len=8) :: check

    ! Init data
    open (udat,file=wfnfile,status='old') 
    iwfn = udat
    read (iwfn,101) wfnttl
    read (iwfn,102) mode, nmo, nprims, ncent
    call allocate_space_for_wfn ()
    do i = 1,ncent
      read (iwfn,103) atnam(i),j,(xyz(j,k),k=1,3),charge(j)
    end do

    ! Read 
    read (iwfn,104) (icen(i),i=1,nprims)
    read (iwfn,104) (ityp(i),i=1,nprims)
    read (iwfn,105) (oexp(i),i=1,nprims)
    do i = 1,nprims
      if (ityp(i).gt.maxtype) then
        call ferror('rdwfn', 'cannot work with , i- or higher primitives', faterr)
      end if
    end do
    do i = 1,nmo
      read (iwfn,106) occ(i), tmp
      if (occ(i).lt.0.0_rp) then
        call ferror('rdwfn', 'nmo with negative occupation : '//string(i), warning)
      end if  
      read (iwfn,107) (coeff(i,j),j=1,nprims)
    end do
    read (iwfn,108) check
    if (check .ne. 'END DATA') then
      call ferror('rdwfn', 'end card 1 not found', faterr)
    endif

    close (iwfn)

    ! formats
101 format (a80)
102 format (4x,a4,10x,3(i5,15x))
103 format (a8,11x,i3,2x,3f12.8,10x,f5.1)
104 format (20x,20i3)
105 format (10x,5e14.7)
106 format (35x,f12.8,15x,f12.8)
107 format (5e16.8)
108 format (a8)
 
  end subroutine

  subroutine setupwfn ()

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
    call alloc ('setupwfn', 'oexpa', oexpa, nprims)
    call alloc ('setupwfn', 'coefa', coefa, nmo, nprims)
    call alloc ('setupwfn', 'icena', icena, nprims)
    call alloc ('setupwfn', 'itypa', itypa, nprims)
    npa = 1_ip
    icena(1) = icen(1)
    oexpa(1) = oexp(1)
    itypa(1) = ityp(1)
    coefa(1:nmo,1) = coeff(1:nmo,1)
    cyclej: do j= 2,nprims
      do m = 1,npa
        if (icen(j).eq.icena(m) .and. ityp(j).eq.itypa(m) .and. abs(oexp(j)-oexpa(m)).le.1d-10) then
          coefa(1:nmo,m) = coefa(1:nmo,m)+coeff(1:nmo,j)
          cycle cyclej
        end if
      end do
      npa = npa + 1_ip
      icena(npa) = icen(j)
      oexpa(npa) = oexp(j)
      itypa(npa) = ityp(j)
      coefa(1:nmo,npa) = coeff(1:nmo,j)
    end do cyclej

    ! Recompute the original variables
    nprimsold = nprims
    nprims = npa
    do j = 1,nprims
      icen(j) = icena(j)
      oexp(j) = oexpa(j)
      ityp(j) = itypa(j)
      coeff(1:nmo,j) = coefa(1:nmo,j)
    end do

    ! Determine primitives corresponding to each center.
    do ic = 1,ncent
      npc(ic) = 0_ip
    end do
    do j = 1,nprims
      ic = icen(j)
      npc(ic) = npc(ic) + 1_ip
      inda = npc(ic)
      icenat(inda,ic) = j
    end do

    ! Classify primitives in each center by types and similar exponents.
    do ic = 1,ncent
      scyclej: do j = 1,npc(ic)
        k = icenat(j,ic)
        itk = ityp(k)
        isuk = nlm(itk,1) + nlm(itk,2) + nlm(itk,3)
        if (j.eq.1) then
          ngroup(ic) = 1_ip
          nzexp(ic,1) = 1_ip
          nuexp(ic,1,1) = k
        else
          do m = 1,ngroup(ic)
            inda = nuexp(ic,m,1)
            itd = ityp(inda)
            isud = nlm(itd,1) + nlm(itd,2) + nlm(itd,3)
            if (abs(oexp(k)-oexp(inda)).lt.1d-8) then
              if (itk.eq.1.and.itd.eq.1) then
                call ferror ('setupwfn', 'two s primitives with equal exponents', faterr)
              else
                if (isuk.eq.isud) then
                  nzexp(ic,m) = nzexp(ic,m) + 1_ip
                  nuexp(ic,m,nzexp(ic,m)) = k
                  cycle scyclej
                end if
              end if
            end if
          end do
          ngroup(ic) = ngroup(ic) + 1_ip
          if (ngroup(ic).gt.mgrp) then
            call ferror ('setupwfn', 'increase mgrp in mod_wfn file', faterr)
          end if
          nzexp(ic,ngroup(ic)) = 1_ip
          nuexp(ic,ngroup(ic),1) = k
        end if   
      end do scyclej
    end do
  
    ! Reconstruct the values of ityp(),icen(),oexp(), and coef()
    i = 0_ip
    do ic = 1,ncent
      do m = 1,ngroup(ic)
        do k = 1,nzexp(ic,m)
          j = nuexp(ic,m,k)
          alph = oexp(j)
          itip = ityp(j)
          icentro = icen(j)
          i = i + 1_ip
          itypa(i) = itip
          oexpa(i) = alph
          icena(i) = icentro
          coefa(1:nmo,i) = coeff(1:nmo,j)
        enddo
      end do
    end do
  
    call free ('setupwfn', 'coeff', coeff)
    call alloc ('setupwfn', 'coeff', coeff, nmo, nprims)
    do i = 1,nprims
      ityp(i) = itypa(i)
      oexp(i) = oexpa(i)
      icen(i) = icena(i)
      coeff(1:nmo,i) = coefa(1:nmo,i)
    end do
   
    npcant = 0_ip
    do ic = 1,ncent
      do k = 1,npc(ic)
        icenat(k,ic) = k + npcant
      end do
      npcant = npcant + npc(ic)
    end do

    ! Reconstruct the values of nuexp(). Now, they are ordered.
    ! Determine also which is the maximum value of ngroup(ic)
    ! Determine the total number of shells.
    i = 0_ip
    numshells = 0_ip
    maxgrp = 0_ip
    do ic = 1,ncent
      numshells = numshells + ngroup(ic)
      if (ngroup(ic).gt.maxgrp) maxgrp = ngroup(ic)
      do m = 1,ngroup(ic)
        do k = 1,nzexp(ic,m)
          i = i + 1_ip
          nuexp(ic,m,k) = i
        end do
      end do
    end do

    ! Allocate space for nshell(), ishell(), and atcenter() arrays.
    call allocate_space_for_shells ()

    ! Deallocate arrays.
    call free ('setupwfn', 'oexpa', oexpa)
    call free ('setupwfn', 'coefa', coefa)
    call free ('setupwfn', 'icena', icena)
    call free ('setupwfn', 'itypa', itypa)

  end subroutine

  subroutine filtergto ()

    implicit none

    logical :: okcen
    integer(kind=ip) :: j, ic, i, jc, lsum, m
    real(kind=rp) :: dis, x1, x2, y1, y2, z1, z2, xmax, zz  

    ! Evaluate internuclear distances
    do i = 1,ncent
      x1 = xyz(i,1)
      y1 = xyz(i,2)
      z1 = xyz(i,3)
      do j = 1,i
        x2 = xyz(j,1)
        y2 = xyz(j,2)
        z2 = xyz(j,3)
        dis = sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
        rint(i,j) = dis
        rint(j,i) = dis
      end do
    end do

    ! Maximum distance at which it is necessary to compute a shell.
    do ic = 1,ncent
      do m = 1,ngroup(ic)
        i = nuexp(ic,m,1)
        lsum = nlm(ityp(i),1)+nlm(ityp(i),2)+nlm(ityp(i),3)
        zz = oexp(i)
        x1 = 0.1_rp
        do 
          if (x1**lsum*exp(-zz*x1*x1).le.abs(epscuttz)) exit
          x1 = x1 + 0.1_rp
        end do
        rcutte(ic,m) = x1
      end do
    end do

    nden(1:ncent) = 0_ip
    nshell(1:ncent) = 0_ip
    do ic = 1,ncent
      xmax = rmaxatom
      do jc = 1,ncent
        dis = rint(ic,jc)
        okcen = .false.
        ! The shell 'm' of center 'jc' contributes to the density, orbitals, 
        ! orbital products, etc, at any point inside the center 'ic' if the 
        ! following condition holds.
        do m = 1,ngroup(jc)
          if (dis.lt.xmax+rcutte(jc,m)) then
            nshell(ic) = nshell(ic) + 1_ip
            atcenter(ic,nshell(ic)) = jc
            ishell(ic,nshell(ic)) = m
            okcen = .true.
          end if
        end do
        if (okcen) then
          nden(ic) = nden(ic) + 1_ip
          iden(ic,nden(ic)) = jc
        end if
      end do
    end do

    do ic = 1,ncent
      do m = 1,ngroup(ic)
        rcutte(ic,m) = rcutte(ic,m)*rcutte(ic,m)
      end do
    end do

  end subroutine

  subroutine optswfn(var,val)

    use iso_fortran_env, only: uout=>output_unit
    use mod_io, only: equal, faterr, ferror, string
    use mod_param, only: verbose
    implicit none

    character(len=*), intent(in) :: var
    real(kind=rp) :: val

    if (equal(var,'epscuttz')) then
      epscuttz = abs(val)
      if (verbose) then
        write (uout,'(1x,a,1x,e13.6)') string('# *** Variable epscuttz changed to :'), epscuttz
      end if
    else if (equal(var,'epsortho')) then
      epsortho = abs(val)
      if (verbose) then
        write (uout,'(1x,a,1x,e13.6)') string('# *** Variable epsortho changed to :'), epsortho
      end if
    else if (equal(var,'rmaxatom')) then
      rmaxatom = abs(val)
      if (verbose) then
        write (uout,'(1x,a,1x,e13.6)') string('# *** Variable rmaxatom changed to :'), rmaxatom
      end if
    else
      call ferror ('optswfn', 'unknown option '//string(var), faterr)
    end if

  end subroutine optswfn

  subroutine init_wfn()
                
    implicit none

    epsortho = 1d-5
    epscuttz = 1d-14
    rmaxatom = 12.0

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

  subroutine end_wfn()

    implicit none

    call deallocate_space_for_wfn ()
    call deallocate_space_for_shells ()

  end subroutine end_wfn

  subroutine allocate_space_for_wfn ()
  
    use mod_memory, only: alloc
    use mod_io, only: faterr, ferror
    implicit none

    integer(kind=ip) :: ier

    call alloc ('mod_wfn', 'coeff', coeff, nmo, nprims)
    call alloc ('mod_wfn', 'npc', npc, ncent)
    call alloc ('mod_wfn', 'ngroup', ngroup, ncent)
    call alloc ('mod_wfn', 'icenat', icenat, nprims, ncent)
    call alloc ('mod_wfn', 'nzexp', nzexp, ncent, mgrp)
    call alloc ('mod_wfn', 'nuexp', nuexp, ncent, mgrp, ngtoH)
    call alloc ('mod_wfn', 'xyz', xyz, ncent, 3)
    call alloc ('mod_wfn', 'oexp', oexp, nprims)
    call alloc ('mod_wfn', 'rcutte', rcutte, ncent, mgrp)
    call alloc ('mod_wfn', 'rint', rint, ncent, ncent)
    call alloc ('mod_wfn', 'iden', iden, ncent, ncent)
    call alloc ('mod_wfn', 'nden', nden, ncent)
    call alloc ('mod_wfn', 'occ', occ, nmo)
    call alloc ('mod_wfn', 'charge', charge, ncent)
    call alloc ('mod_wfn', 'icen', icen, nprims)
    call alloc ('mod_wfn', 'ityp', ityp, nprims)
    if (.not.allocated(atnam)) then
      allocate (atnam(ncent),stat=ier) 
      if (ier.ne.0) then
        call ferror('mod_wfn', 'cannot allocate atnam', faterr)
      end if
    end if
 
  end subroutine allocate_space_for_wfn

  subroutine deallocate_space_for_wfn ()

    use mod_io, only: faterr, ferror
    use mod_memory, only: free
    implicit none

    integer(kind=ip) :: ier

    call free ('mod_wfn', 'coeff', coeff)
    call free ('mod_wfn', 'npc', npc)
    call free ('mod_wfn', 'ngroup', ngroup)
    call free ('mod_wfn', 'icenat', icenat)
    call free ('mod_wfn', 'nzexp', nzexp)
    call free ('mod_wfn', 'nuexp', nuexp)
    call free ('mod_wfn', 'xyz', xyz)
    call free ('mod_wfn', 'oexp', oexp)
    call free ('mod_wfn', 'rcutte', rcutte)
    call free ('mod_wfn', 'rint', rint)
    call free ('mod_wfn', 'iden', iden)
    call free ('mod_wfn', 'nden', nden)
    call free ('mod_wfn', 'occ', occ)
    call free ('mod_wfn', 'charge', charge)
    call free ('mod_wfn', 'icen', icen)
    call free ('mod_wfn', 'ityp', ityp)
    if (allocated(atnam)) then
      deallocate (atnam,stat=ier) 
      if (ier.ne.0) then
        call ferror('mod_wfn', 'cannot deallocate atnam', faterr)
      end if
    end if

  end subroutine deallocate_space_for_wfn
                                                                        
  subroutine allocate_space_for_shells ()

    use mod_memory, only: alloc
    implicit none
  
    call alloc ('mod_wfn', 'nshell', nshell, ncent)
    call alloc ('mod_wfn', 'ishell', ishell, ncent, numshells)
    call alloc ('mod_wfn', 'atcenter', atcenter, ncent, numshells)

  end subroutine allocate_space_for_shells
                                                                        
  subroutine deallocate_space_for_shells ()

    use mod_memory, only: free
    implicit none

    call free ('mod_wfn', 'nshell', nshell)
    call free ('mod_wfn', 'ishell', ishell)
    call free ('mod_wfn', 'atcenter', atcenter)

  end subroutine deallocate_space_for_shells

  subroutine isorthowfn ()

    use iso_fortran_env, only: uout=>output_unit
    use mod_memory, only: alloc, free
    use mod_param, only: verbose, debug
    use mod_io, only: string, ferror, warning
    implicit none
  
    logical :: ortho 
    real(kind=rp) :: solap
    integer(kind=ip) :: i, j
    real(kind=rp), allocatable, dimension(:,:) :: sprim, overlap, coeftmp

    call alloc ('isorthowfn', 'sprim', sprim , nprims, nprims)
    call alloc ('isorthowfn', 'overlap', overlap , nmo, nmo)
    call alloc ('isorthowfn', 'coeftmp', coeftmp , nmo, nprims)
    call gtogto (sprim)
 
    if (verbose) then
      write (uout,'(1x,a)') string('# Testing orthogonality of canonical MOs')
    end if
    ortho = .true.
    coeftmp = matmul(sprim,transpose(coeff))
    overlap = matmul(coeff,coeftmp)
    do i = 1,nmo
      do j = 1,i
        solap = overlap(i,j)
        if (i.eq.j) then
          if (abs(abs(solap)-1.0_rp) .gt. epsortho) then
            ortho = .false.
            if (debug) write (uout,222) i,solap
          end if
        else
          if (abs(solap) .gt. epsortho) then
            ortho = .false.
            if (debug) write (uout,223) i,j,solap
          end if
        end if
      end do 
    end do
    if (.not.ortho) then
      call ferror ('isorthowfn', 'the set of canonical mos are not orthonormal', warning)
    else
      if (verbose) then
        write (uout,'(1x,a)') string('# Orthogonal canonical MOs')
      end if
    end if

    call free ('isorthowfn', 'sprim', sprim)
    call free ('isorthowfn', 'overlap', overlap)
    call free ('isorthowfn', 'coeftmp', coeftmp)

222 format (1x,'# MO number ',i0, ' is not exactly normalized, NORM = ', e19.12)
223 format (1x,'# MOs ',i0,' and ',i0, ' are not exactly orthogonal, S = ', e19.12)

  end subroutine isorthowfn

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

  ! Overlap matrix between primitive Cartesian Gaussian Functions
  subroutine gtogto (sprim)

    use mod_param, only: pi
    implicit none
    integer(kind=ip), parameter :: lamx = 12
 
    real(kind=rp), intent(out) :: sprim(nprims,nprims)
 
    real(kind=rp), dimension(ncent,ncent) :: ab2
    real(kind=rp) :: ax(1:3), bx(1:3), za, zb, p, pioverp, abaux
    real(kind=rp) :: prefactor, prod, xmu
    integer(kind=ip) :: i, j, k, l, m, ica, icb, nua, nub, n
    integer(kind=ip) :: itipa, itipb, la, lb, ka, kb, ma, mb
    ! ceabx() are the coefficients, except for the factor 
    ! EXP(-XMU*R_AB^2), where XMU=a*b/(a+b), that result from the 
    ! expansion of the product of two primitive cartesian Gaussian
    real(kind=rp) :: ceabx(-1:2*lamx,-1:lamx,-1:lamx,3)
 
    do ica = 1,ncent
      do icb = 1,ica
        abaux = 0.0_rp
        do j = 1,3
          abaux = abaux + (xyz(ica,j)-xyz(icb,j))**2
        end do
        ab2(ica,icb) = abaux
        ab2(icb,ica) = abaux
      end do
    end do
   
    ! Compute the electronic molecular electrostatic potential.
    do ica = 1,ncent
      ax(1:3) = xyz(ica,1:3)
      do ma = 1,ngroup(ica)   
        nua = nuexp(ica,ma,1)
        itipa = ityp(nua)
        la = nlm(itipa,1) + nlm(itipa,2) + nlm(itipa,3)
        za = oexp(nua)
        do icb = 1,ica
          bx(1:3) = xyz(icb,1:3)
          do mb = 1,ngroup(icb)
            nub = nuexp(icb,mb,1)
            itipb = ityp(nub)
            lb = nlm(itipb,1) + nlm(itipb,2) + nlm(itipb,3)
            zb = oexp(nub)
            p = za + zb
            xmu = za*zb/p
            prefactor = exp(-xmu*ab2(ica,icb))
            pioverp = pi/p
            pioverp = sqrt(pioverp*pioverp*pioverp)
            do j = 1,3
              call etijcalc (j,lamx,la,lb,ceabx,za,zb,ax(j),bx(j))
            end do
            ! Compute the target functions for all the products of
            ! of Gaussian primitives.
            do ka = 1,nzexp(ica,ma)
              nua = nuexp(ica,ma,ka)
              itipa = ityp(nua)
              i = nlm(itipa,1)
              k = nlm(itipa,2)
              m = nlm(itipa,3)
              do kb = 1,nzexp(icb,mb)
                nub = nuexp(icb,mb,kb)
                if (nua.ge.nub) then
                  itipb = ityp(nub)
                  j = nlm(itipb,1)
                  l = nlm(itipb,2)
                  n = nlm(itipb,3)
                  prod = ceabx(0,i,j,1)*ceabx(0,k,l,2)*ceabx(0,m,n,3)
                  sprim(nua,nub) = prod*pioverp*prefactor
                  sprim(nub,nua) = sprim(nua,nub)
                end if
              end do
            end do
          end do
        end do
      end do
    end do
   
  end subroutine

  subroutine info_wfn()

    use iso_fortran_env, only: uout=>output_unit
    use mod_io, only: string
    implicit none
    integer(kind=ip), parameter :: ncentw = 100

    integer(kind=ip) :: i, j, m, ic, k, lsum
    character(len=1) :: lb(0:5), lbl
    logical :: wrout
    data (lb(i),i=0,5) /'S','P','D','F','G','H'/

    if (ncent.gt.ncentw) then
      wrout = .false.
    else
      wrout = .true.
    end if

    write (uout,'(1x,a,1x,i0)') string('# Number of centers :'), ncent
    atnam = adjustl(atnam)
    do i = 1,ncent
      write (uout,'(1x,a,1x,i0,1x,a,1x,f4.1,1x,3f12.6)') '#', i, &
                    atnam(i)(1:2), charge(i), xyz(i,:)
    end do

    write (uout,'(1x,a,1x,e13.6)') string('# GTO eps :'), epscuttz
    write (uout,'(1x,a,1x,i0)') string('# Original number of primitives :'), nprimsold
    write (uout,'(1x,a,1x,i0)') string('# Actual number of primitives :'), nprims
    write (uout,'(1x,a,1x,i0)') string('# Total number of shells :'), numshells
    write (uout,'(1x,a,1x,e13.6)') string('# Rmaxatom for shells :'), rmaxatom
    if (wrout) then
      do ic = 1,ncent
        write (uout,210) ic
        write (uout,300) nshell(ic),ic
        write (uout,301) (ishell(ic,j),atcenter(ic,j),j=1,nshell(ic))
        do m = 1,ngroup(ic)
          i = nuexp(ic,m,1)
          lsum = nlm(ityp(i),1)+nlm(ityp(i),2)+nlm(ityp(i),3)
          lbl = lb(lsum)
          write (uout,613) lbl,oexp(i),sqrt(rcutte(ic,m)),(nuexp(ic,m,k),k=1,nzexp(ic,m))
        end do
      end do
    end if

    write (uout,'(1x,a,1x,e13.6)') string('# Orthogonality eps :'), epsortho
    write (uout,'(1x,a,1x,i0)') string('# Number of molecular orbitals :'), nmo

300 format (1x,'# ',i0,' shells contribute to the basin of center ',i0, &
     /,' # [ shell(atom) means shell number "shell" of atom "atom" ]')
301 format (1000(1x,8(I6,'(',I4,')'),/))
210 format (1x,'# CENTER ',i0)
613 format (1x,'# ',a,' Shell   Exp = ',e16.8,4x,'Cutoff = ',f13.6, &
            4x,'Primitives : ',21(1x,i0))

  end subroutine info_wfn

end module mod_wfn
