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
  integer(kind=ip), parameter :: maxtype = 56_ip

  ! nlm keeps the nlm values of x^n y^l z^m gaussian
  integer(kind=ip) :: nlm(maxtype,3)

  integer(kind=ip) :: nprims
  integer(kind=ip) :: ncore
  integer(kind=ip) :: nmo
  integer(kind=ip) :: ncent
  integer(kind=ip) :: maxgrp
  integer(kind=ip) :: numshells

  real(kind=rp), allocatable, dimension(:,:) :: coefcan
  real(kind=rp), allocatable, dimension(:,:) :: coefnat
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

  integer(kind=ip) :: nvirtual, noccupied
  integer(kind=ip), allocatable, dimension(:) :: occupied
  integer(kind=ip), allocatable, dimension(:) :: virtual
  real(kind=rp), allocatable, dimension(:) :: occv

  integer(kind=ip), allocatable, dimension(:,:) :: ishell
  integer(kind=ip), allocatable, dimension(:) :: nshell
  integer(kind=ip), allocatable, dimension(:,:) :: atcenter

  real(kind=rp), allocatable, dimension(:,:) :: c1et

  logical :: rdm
  ! OPTIONS:
  ! Any gaussian primitive and gaussian derivative will be assumed
  ! to be zero if it is smaller than cuttz.
  real(kind=rp) :: cuttz
  real(kind=rp) :: epsocc
  real(kind=rp) :: rmaxatom
  real(kind=rp) :: epsortho
  
contains

  subroutine allocate_space_for_wfn ()
  
    use mod_memory, only: alloc
    use mod_io, only: faterr, ferror
    implicit none

    integer(kind=ip) :: ier

    call alloc ('mod_wfn', 'coefcan', coefcan, nmo, nprims)
    call alloc ('mod_wfn', 'coefnat', coefnat, nmo, nprims)
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

    call free ('mod_wfn', 'coefcan', coefcan)
    call free ('mod_wfn', 'coefnat', coefnat)
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

  subroutine allocate_space_for_rdm ()

    use mod_memory, only: alloc
    implicit none

    call alloc ('mod_wfn', 'c1et', c1et, nmo, nmo)

  end subroutine allocate_space_for_rdm
                                                                        
  subroutine deallocate_space_for_rdm ()

    use mod_memory, only: free
    implicit none

    call free ('mod_wfn', 'c1et', c1et)

  end subroutine deallocate_space_for_rdm

  subroutine allocate_space_for_rho ()

    use mod_memory, only: alloc
    implicit none

    call alloc ('mod_wfn', 'occv', occv, noccupied)
    call alloc ('mod_wfn', 'occupied', occupied, noccupied)

  end subroutine allocate_space_for_rho
                                                                        
  subroutine deallocate_space_for_rho ()

    use mod_memory, only: free
    implicit none

    call free ('mod_wfn', 'occupied', occupied)
    call free ('mod_wfn', 'occv', occv)

  end subroutine deallocate_space_for_rho

end module mod_wfn
