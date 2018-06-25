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
! check the closure of the symmetry operators set
! and get the properties of the symmetry operators.
!
subroutine closuresym ()

  use iso_fortran_env, only: uout=>output_unit
  use mod_io, only: ferror, faterr
  use mod_prec, only: rp, ip
  use mod_param, only: verbose
  use mod_sym, only: opsym, mopsym, mclas, optable, opclas, &
                     opinv, nopsym, opisgener, opt_identity, &
                     minclas, optype, nclas, ninclas, iinclas
  implicit none

  real(kind=rp) :: xmat(3,3)
  integer(kind=ip) :: i, j, ii, jj, kk, inum, old_nopsym, ini_nopsym
  integer(kind=ip) :: iop, jinv, iclas
  logical :: newop, found, generated(MOPSYM)

  interface
    integer(kind=ip) function numbersym (xmat)
      import rp, ip
      implicit none
      real(kind=rp), dimension(3,3), intent(in) :: xmat
    end function
  end interface
                   
  ! Initialize the Cayley table:
  if (verbose) write (uout,600)
  do i = 1,MOPSYM
    do j = 1,MOPSYM
      optable(i,j) = 0
    end do
  end do

  ! Build up the Cayley table and look for new operators at the same time
  ini_nopsym = nopsym
  newop = .true.
  do while (newop)
    newop = .false.
    old_nopsym = nopsym
    do i = 1,old_nopsym
      do j = 1,old_nopsym
        if (optable(i,j) .le. 0) then
          ! Get the product of operators i x j and determine if
          ! this is a new or an old symmetry operator:
          do ii = 1,3
            do jj = 1,3
              xmat(ii,jj) = 0.0_rp
              do kk = 1,3
                xmat(ii,jj) = xmat(ii,jj) + opsym(i,ii,kk)*opsym(j,kk,jj)
              end do
            end do
          end do
          inum = numbersym (xmat)
          if (inum .le. 0) then
            newop = .true.
            optable(i,j) = nopsym + 1
            call opaddsym (xmat)
          else
            optable(i,j) = inum
          end if
        end if
      end do
    end do
  end do

  ! Determine the operation inverses:
  do i = 1,nopsym
    j = 1
    found = .false.
    do while (.not.found .and. j.le.nopsym)
      kk = optable(i,j)
      if (optype(kk) .eq. opt_identity) then
        opinv(i) = j
        found = .true.
      else
        j = j + 1
      end if
    end do
    if (.not.found .and. verbose) write (uout,605) i
  end do

  ! Check for a set of generators:
  do i = 1,nopsym
    generated(i) = .false.
    opisgener(i) = .false.
  end do
  do i = 1,nopsym
    if (.not.generated(i) .and. optype(i).ne.opt_identity) then
      opisgener(i) = .true.
    end if
    do j = 1,i
      kk = optable(i,j)
      generated(kk) = .true.
      kk = optable(j,i)
      generated(kk) = .true.
    end do
  end do
 
  ! Determine the classes of the symmetry operations:
  ! (Non-essential)
  ! A and B belong in the same class if there exist any operation X in
  ! the group such that: X^{-1} A X = B.
  do i = 1,nopsym
    opclas(i) = -1
  end do
  nclas = 0
  do i = 1,nopsym
    if (opclas(i) .le. 0) then
      ! This op belongs to a new class:
      nclas = nclas + 1
      if (nclas.gt.MCLAS) then
        call ferror ('closuresym', 'too many classes', faterr)
      end if
      opclas(i) = nclas
      iclas = nclas
      ninclas(nclas) = 1
      iinclas(nclas,1) = i
    else
      ! The op belongs to a class already found:
      iclas = opclas(i)
    end if
    ! Apply similarity transforms to get all ops in the same class:
    do j = 1,nopsym
      jinv = opinv(j)
      iop = optable(optable(jinv, i), j)
      if (opclas(iop) .le. 0) then
        opclas(iop) = iclas
        ninclas(nclas) = ninclas(nclas) + 1
        if (ninclas(nclas).gt.MINCLAS) then
          call ferror ('symclosure', 'too many ops in a class', faterr)
        end if
        iinclas(nclas,ninclas(nclas)) = iop
      else if (opclas(iop) .ne. iclas) then
        ! This is an error: iop and i should be in the same class,
        ! but they are not.
        if (verbose) write (uout,620) i, iop
      end if
    end do
  end do
 
  if (nopsym.gt.ini_nopsym .and. verbose) then
    write (uout,600)
    write (uout,602) nopsym-ini_nopsym
  end if
 
600 format (/                                                     &
    1x, '++SYMCLOSURE: Complete the group by multiplying the',    &
    1x, 'matrices already known.')

602 format (                                                      &
    1x, 'ALG(symclosure) ', i3,                                   &
    ' new operations found because of the group closure')

605 format (                                                      &
    1x, 'DBG(symclosure): Inverse not found for opsym ', i3)

620 format (                                                      &
    1x, 'WARNING(symclosure): The ops', 2i5, ' should be in the', &
    1x, 'same class but they are not. We continue anyway!')
 
end subroutine closuresym
