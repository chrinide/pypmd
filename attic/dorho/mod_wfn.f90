module mod_wfn
  
  use mod_prec, only: ip, rp
  implicit none
  public
  
  integer(kind=ip), parameter :: mgrp = 500_ip
  integer(kind=ip), parameter :: ngtoh = 21_ip
  integer(kind=ip), parameter :: maxtype = 56_ip

  ! nlm keeps the nlm values of x^n y^l z^m gaussian
  integer(kind=ip) :: nlm(maxtype,3)

  integer(kind=ip) :: nprims
  integer(kind=ip) :: nmo
  integer(kind=ip) :: ncent
  integer(kind=ip) :: maxgrp
  integer(kind=ip) :: numshells

  real(kind=rp), allocatable, dimension(:,:) :: coef
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
  real(kind=rp), allocatable, dimension(:) :: eorb
  real(kind=rp), allocatable, dimension(:) :: charge
  character(len=8), allocatable, dimension(:) :: atnam
  integer(kind=ip), allocatable, dimension(:) :: icen
  integer(kind=ip), allocatable, dimension(:) :: ityp
  integer(kind=ip), allocatable, dimension(:,:):: iden
  integer(kind=ip), allocatable, dimension(:) :: nden

  integer(kind=ip), allocatable, dimension(:,:) :: ishell
  integer(kind=ip), allocatable, dimension(:) :: nshell
  integer(kind=ip), allocatable, dimension(:,:) :: atcenter

  real(kind=rp), allocatable, dimension(:,:) :: c1et
  real(kind=rp), allocatable, dimension(:,:) :: c1ea, c1eb
  real(kind=rp), allocatable, dimension(:,:) :: c2e

  logical :: corr
  logical :: wfnat
  logical :: rhf
  logical :: rohf
  logical :: uhf
  logical :: cciqa
  logical :: mpiqa
  logical :: dmftxc

  integer(kind=ip), allocatable, dimension(:) :: ialpha, ibeta
  real(kind=rp) :: multiplicity, totel
  integer(kind=ip) :: ndets, ncore, nact, nelact
  integer(kind=ip) :: nel, ndou, nsin, nalpha, nbeta
  integer(kind=ip) :: norbp

  integer(kind=ip), allocatable, dimension(:,:) :: ipak
  real(kind=rp), allocatable, dimension(:) :: cidet
  !allocate (cidet(0:nmo+nmo),stat=ier) 

  ! Options
  ! Integer parameter that defines the type of Density Matrix
  ! Functional Theory (DMFT) exchange-correlation that is used.
  integer(kind=ip) :: idmft
  ! Descard elements if 
  real(kind=rp) :: epsrho, epsdet

contains

  subroutine init_wfn

    implicit none

    idmft = 2_ip ! Default XC functional when dmftxc becomes .true.
    epsdet = 1d-12
    epsrho = 1d-14

    corr = .false.
    wfnat = .false.
    rhf = .false.
    rohf = .false.
    uhf = .false.
    mpiqa = .false.
    dmftxc = .false.
    cciqa = .false.

    multiplicity = 0.0_rp
    totel = 0.0_rp
    ndets = 0_ip
    ncore = 0_ip
    nact = 0_ip
    nelact = 0_ip
    norbp = 0_ip
    nel = 0_ip
    ndou = 0_ip
    nsin = 0_ip
    nalpha = 0_ip
    nbeta = 0_ip

!...p's
    nlm(2,1)=1      !px
    nlm(2,2)=0      !px
    nlm(2,3)=0      !px
    nlm(3,1)=0      !py
    nlm(3,2)=1      !py
    nlm(3,3)=0      !py
    nlm(4,1)=0      !pz
    nlm(4,2)=0      !pz
    nlm(4,3)=1      !pz
!...d's
    nlm(5,1)=2      !xx
    nlm(6,2)=2      !yy
    nlm(7,3)=2      !zz
    nlm(8,1)=1      !xy
    nlm(8,2)=1
    nlm(9,1)=1      !xz
    nlm(9,3)=1
    nlm(10,2)=1     !yz
    nlm(10,3)=1
!...f's
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
!...g's
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
!...h's
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

end module mod_wfn
!-----------------------------------------------------------------------
subroutine allocate_space_for_wfn (ncent,nprims,nmo,mgrp,ngtoH)
use mod_prec, only: ip
use mod_io, only: faterr, ferror
use mod_wfn, only: coef, npc, ngroup, icenat, nzexp, nuexp
use mod_wfn, only: xyz, oexp, rcutte, rint, occ, eorb, &
                   charge, icen, ityp, iden, nden, atnam, &
                   ialpha, ibeta 
use mod_memory, only: alloc
implicit none
integer(kind=ip) :: ier
integer(kind=ip), intent(in) :: nprims, ncent, mgrp, ngtoH, nmo
call alloc ('mod_wfn', 'coef', coef, nmo+nmo, nprims)
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
call alloc ('mod_wfn', 'eorb', eorb, nmo)
call alloc ('mod_wfn', 'charge', charge, ncent)
call alloc ('mod_wfn', 'icen', icen, nprims)
call alloc ('mod_wfn', 'ityp', ityp, nprims)
call alloc ('mod_wfn', 'ialpha', ialpha, nmo)
call alloc ('mod_wfn', 'ibeta', ibeta, nmo)
if (.not.allocated(atnam)) then
  allocate (atnam(ncent),stat=ier) 
  if (ier.ne.0) then
    call ferror('mod_wfn', 'cannot allocate atnam', faterr)
  end if
end if
end subroutine allocate_space_for_wfn
!-----------------------------------------------------------------------
subroutine deallocate_space_for_wfn
use mod_prec, only: ip
use mod_io, only: faterr, ferror
use mod_wfn, only: coef, npc, ngroup, icenat, nzexp, nuexp
use mod_wfn, only: xyz, oexp, rcutte, rint, occ, eorb, &
                   charge, icen, ityp, iden, nden, atnam, &
                   ialpha, ibeta 
use mod_memory, only: free
implicit none
integer(kind=ip) :: ier
call free ('mod_wfn', 'coef', coef)
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
call free ('mod_wfn', 'eorb', eorb)
call free ('mod_wfn', 'charge', charge)
call free ('mod_wfn', 'icen', icen)
call free ('mod_wfn', 'ityp', ityp)
call free ('mod_wfn', 'ialpha', ialpha)
call free ('mod_wfn', 'ibeta', ibeta)
if (allocated(atnam)) then
  deallocate (atnam,stat=ier) 
  if (ier.ne.0) then
    call ferror('mod_wfn', 'cannot deallocate atnam', faterr)
  end if
end if
end subroutine deallocate_space_for_wfn
!-----------------------------------------------------------------------
subroutine allocate_space_for_shells (ncent,numshells)
use mod_prec, only: ip
use mod_wfn, only: nshell, ishell, atcenter
use mod_memory, only: alloc
implicit none
integer(kind=ip), intent(in) :: ncent, numshells
call alloc ('mod_wfn', 'nshell', nshell, ncent)
call alloc ('mod_wfn', 'ishell', ishell, ncent, numshells)
call alloc ('mod_wfn', 'atcenter', atcenter, ncent, numshells)
end subroutine allocate_space_for_shells
!-----------------------------------------------------------------------
subroutine deallocate_space_for_shells
use mod_wfn, only: nshell, ishell, atcenter
use mod_memory, only: free
implicit none
call free ('mod_wfn', 'nshell', nshell)
call free ('mod_wfn', 'ishell', ishell)
call free ('mod_wfn', 'atcenter', atcenter)
end subroutine deallocate_space_for_shells
!-----------------------------------------------------------------------
subroutine allocate_space_for_rdm (nmo)
use mod_prec, only: ip
use mod_wfn, only: c1ea, c1eb, c1et, c2e, ipak
use mod_memory, only: alloc
implicit none
integer(kind=ip), intent(in) :: nmo
integer(kind=ip) :: nmopair, nsize
nmopair = (nmo*(nmo+1))/2
nsize = nmo*(nmo+2)
call alloc ('mod_wfn', 'c1et', c1et, nmo, nmo)
call alloc ('mod_wfn', 'c1ea', c1ea, nmo, nmo)
call alloc ('mod_wfn', 'c1eb', c1eb, nmo, nmo)
call alloc ('mod_wfn', 'c2e', c2e, nmopair, nmopair)
call alloc ('mod_wfn', 'ipak', ipak, nsize, 2)
end subroutine allocate_space_for_rdm
!-----------------------------------------------------------------------
subroutine deallocate_space_for_rdm
use mod_wfn, only: c1ea, c1eb, c1et, c2e, ipak
use mod_memory, only: free
implicit none
call free ('mod_wfn', 'c1et', c1et)
call free ('mod_wfn', 'c1ea', c1ea)
call free ('mod_wfn', 'c1eb', c1eb)
call free ('mod_wfn', 'c2e', c2e)
call free ('mod_wfn', 'ipak', ipak)
end subroutine deallocate_space_for_rdm
