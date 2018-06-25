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
! check symmetry of nonlinear molecules.
!
subroutine d3sym ( )

  use iso_fortran_env, only: uout=>output_unit
  use mod_io, only: ferror, faterr
  use mod_prec, only: ip, rp
  use mod_param, only: debug, verbose
  use mod_math, only: matfill, inv
  use mod_sym, only: nmol=>sdim, toldirty, orbmol, toltriplet, &
                     ax, ay, az, linear_mol, errsng, nopsym
                     
  implicit none

  real(kind=rp) :: dbest, xdet
  real(kind=rp), dimension(3,3) :: xmat, xm, xv, xop
  integer(kind=ip) :: k, ii, jj, i1, i2, i3, j1, j2
  integer(kind=ip) :: ntest1, ntest2, j3, ierr, ii1, ii2, ii3 

  interface
    logical function opchksym(xmat)
      import ip, rp
      implicit none
      real(kind=rp), dimension(3,3), intent(in) :: xmat
    end function
  end interface

  ! Classify all atoms into orbits. All atoms in an orbit have the
  ! same atomic number and the same distance to the center of mass:
  call orbsym ()

  linear_mol = .false.
  nopsym = 0

  ! The identity is always a sym operator:
  call matfill (xmat, 1.0_rp,0.0_rp,0.0_rp, &
                      0.0_rp,1.0_rp,0.0_rp, &
                      0.0_rp,0.0_rp,1.0_rp)
  call opaddsym (xmat)

  ! Is the inversion a sym op?
  call matfill (xmat, -1.0_rp,0.0_rp,0.0_rp, &
                       0.0_rp,-1.0_rp,0.0_rp, &
                       0.0_rp,0.0_rp,-1.0_rp)
  if (opchksym(xmat)) then
    call opaddsym (xmat)
  end if

  ! Find a linear independent triplet of atoms.
  ! Try first to find a very good triplet, or use the best available
  ! one anyway.
  dbest = 0.0_rp
  do i1 = 1,nmol
    xmat(1,1) = ax(i1)
    xmat(2,1) = ay(i1)
    xmat(3,1) = az(i1)
    do i2 = i1+1,nmol
      xmat(1,2) = ax(i2)
      xmat(2,2) = ay(i2)
      xmat(3,2) = az(i2)
      do i3 = i2+1,nmol
        xmat(1,3) = ax(i3)
        xmat(2,3) = ay(i3)
        xmat(3,3) = az(i3)
        call inv (xmat, xv, xdet, ierr)
        if (ierr .eq. 0) then
          if (abs(xdet) .gt. TOLtriplet) goto 1001
          if (abs(xdet) .gt. dbest) then
            dbest = abs(xdet)
            ii1 = i1
            ii2 = i2
            ii3 = i3
          end if
        end if
      end do ! i3
    end do ! i2
  end do ! i1
  xmat(1,1) = ax(ii1)
  xmat(2,1) = ay(ii1)
  xmat(3,1) = az(ii1)
  xmat(1,2) = ax(ii2)
  xmat(2,2) = ay(ii2)
  xmat(3,2) = az(ii2)
  xmat(1,3) = ax(ii3)
  xmat(2,3) = ay(ii3)
  xmat(3,3) = az(ii3)
  call inv (xmat, xv, xdet, ierr)
  if (ierr .ne. 0) then
    call ferror('sym3d','singular triplet matrix',faterr)
  end if

! Run over all triplets that are compatible with the linear
! independent triplet already found. Each compatible triplet migth
! produce a new symmetry operator. To be compatible, both triplets
! must be formed by the same atoms in the same order.
1001 continue

  ERRsng = abs(xdet)
  ntest1 = 0
  ntest2 = 0
  do j1 = 1,nmol
    if (orbmol(i1).ne.orbmol(j1)) goto 1002
    xm(1,1) = ax(j1)
    xm(2,1) = ay(j1)
    xm(3,1) = az(j1)
    do j2 = 1,nmol
      if (orbmol(i2).ne.orbmol(j2)) goto 1003
      if (j1.eq.j2) goto 1003
      xm(1,2) = ax(j2)
      xm(2,2) = ay(j2)
      xm(3,2) = az(j2)
      do j3 = 1,nmol
        if (orbmol(i3).ne.orbmol(j3)) goto 1004
        if (j1.eq.j3 .or. j2.eq.j3) goto 1004
        ntest1 = ntest1 + 1
        xm(1,3) = ax(j3)
        xm(2,3) = ay(j3)
        xm(3,3) = az(j3)
        do ii = 1,3
          do jj = 1,3
            xop(ii,jj) = 0.0_rp
            do k = 1, 3
              xop(ii,jj) = xop(ii,jj) + xm(ii,k)*xv(k,jj)
            end do
          end do
        end do
        if (debug) then
          write (uout,520) j1,j2,j3
          write (uout,525) 'xmati', ((xmat(ii,jj), jj=1,3), ii=1,3)
          write (uout,525) 'xinvi', ((xv(ii,jj), jj=1,3), ii=1,3)
          write (uout,525) 'xmatj', ((xm(ii,jj), jj=1,3), ii=1,3)
          write (uout,525) 'xop', ((xop(ii,jj), jj=1,3), ii=1,3)
        end if
        ! Check if this is a new sym operator:
        if (opchksym(xop)) then
          ntest2 = ntest2 + 1
          call opaddsym (xop)
          if (TOLdirty) call closuresym (uout)
          if (debug) then
            write (uout,535) 'xmati', ((xmat(ii,jj), jj=1,3), ii=1,3)
            write (uout,535) 'xinvi', ((xv(ii,jj), jj=1,3), ii=1,3)
            write (uout,535) 'xmatj', ((xm(ii,jj), jj=1,3), ii=1,3)
            write (uout,535) 'Operator xop', ((xop(ii,jj), jj=1,3), ii=1,3)
          end if
        end if
1004    continue
      end do  !j3
1003  continue
    end do  !j2
1002 continue
  end do  !j1
 
  if (verbose) write (uout,60) ntest1, ntest2
  
  ! Check the closure of the symmetry operators set:
  call closuresym ()
  call getgroup ()

520 format(1x, 'DBG(sym3d) Triplet: ', 6i5)
525 format(1x, 'DBG(sym3d) ', a, ' matrix: '/ (3f15.6))
535 format(1x, 'DBG(sym3d) ', a, ' matrix: '/ (3f15.6))
60 format (                                           &
   1x, 'ALG(sym3d) Triplet pairs tested:     ', i12/  &
   1x, 'ALG(sym3d) Possible operators found: ', i12)

end subroutine 
