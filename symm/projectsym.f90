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
! Clean the atomic coordinates by projecting the atoms
! quite close to the symmetry elements onto those symmetry elements.
! This routine is suppossed to work after sympurify() has cleaned
! the symmetry matrices, but the algorithm do not depends on that.
!
subroutine projectsym ()

  use iso_fortran_env, only: uout=>output_unit
  use mod_io, only: ferror, faterr
  use mod_prec, only: ip, rp
  use mod_param, only: verbose
  use mod_sym, only: toldist, n=>sdim, opt_inversion, ax, ay, az, nopsym, &
                     norbit, opt_imp_rotation, opt_rotation, iatorb, natorb, &
                     opaxis, opt_sigma, optype, opsym
  implicit none

  integer(kind=ip)  :: iop, iorb, i, iat, k, kmin, ntrouble
  real(kind=rp) :: dd, xproj, aproj(3), dmin, dtrouble, dproj
  real(kind=rp) :: ximg, yimg, zimg, xx, yy, zz
  real(kind=rp) :: axx(n), ayy(n), azz(n)
  logical  :: cleaned(n), found, trouble

  ! Keep a copy of the original coordinates:
  axx = ax
  ayy = ay
  azz = az

 
  ! Project first the first atom in each atomic orbit.
  if (verbose) write (uout,600)
  do iop = 1, nopsym
    if (optype(iop).eq.opt_inversion) then
      ! Move atoms quite close to the inversion center:
      do iorb = 1,norbit
        do i = 1,natorb(iorb)
          iat = iatorb(iorb,i)
          dd = sqrt(ax(iat)**2 + ay(iat)**2 + az(iat)**2)
          if (dd .le. TOLdist) then
            if (verbose) write (uout,605) iat, dd, ax(iat), ay(iat), az(iat)
            ax(iat) = 0.0_rp
            ay(iat) = 0.0_rp
            az(iat) = 0.0_rp
          end if
        end do
      end do
    ! Atoms too close to a rotation axis will be projected onto the axis:
    else if (optype(iop).eq.opt_rotation .or. optype(iop).eq.opt_imp_rotation) then
      do iorb = 1,norbit
        do i = 1,natorb(iorb)
          iat = iatorb(iorb,i)
          xproj = ax(iat)*opaxis(iop,1)+ ay(iat)*opaxis(iop,2) + az(iat)*opaxis(iop,3)
          aproj(1) = xproj*opaxis(iop,1)
          aproj(2) = xproj*opaxis(iop,2)
          aproj(3) = xproj*opaxis(iop,3)
          xx = aproj(1) - ax(iat)
          yy = aproj(2) - ay(iat)
          zz = aproj(3) - az(iat)
          dd = sqrt(xx*xx + yy*yy + zz*zz)
          if (dd .le. TOLdist) then
            if (verbose) write (uout,610) iat, iop, dd, ax(iat), ay(iat), &
                                          az(iat), aproj(1), aproj(2), aproj(3)
            ax(iat) = aproj(1)
            ay(iat) = aproj(2)
            az(iat) = aproj(3)
          end if
        end do
      end do
    ! Atoms too close to a symmetry plane will be projected onto the plane:
    else if (optype(iop).eq.opt_sigma) then
      do iorb = 1,norbit
        do i = 1,natorb(iorb)
          iat = iatorb(iorb,i)
          xproj = ax(iat)*opaxis(iop,1) + ay(iat)*opaxis(iop,2) + az(iat)*opaxis(iop,3)
          if (abs(xproj) .le. TOLdist) then
            aproj(1) = ax(iat) - xproj*opaxis(iop,1)
            aproj(2) = ay(iat) - xproj*opaxis(iop,2)
            aproj(3) = az(iat) - xproj*opaxis(iop,3)
            if (verbose) write (uout,615) iat, iop, xproj, ax(iat), ay(iat), &
                                          az(iat), aproj(1), aproj(2), aproj(3)
            ax(iat) = aproj(1)
            ay(iat) = aproj(2)
            az(iat) = aproj(3)
          end if
        end do
      end do
    end if
  end do
 
  ! Regenerate each orbit by applying the symmetry matrices to the
  ! first orbit atom.
  do i = 1, n
    cleaned(i) = .false.
  enddo
  ntrouble = 0
  dtrouble = 0.0_rp
  trouble = .false.
  do iorb = 1,norbit
    do i = 1,natorb(iorb)
      iat = iatorb(iorb,i)
      cleaned(iat) = .true.
      do iop = 2,nopsym
        ! Get the image of the atom by the symmetry operation:
        ximg = opsym(iop,1,1)*ax(iat) + opsym(iop,1,2)*ay(iat) + opsym(iop,1,3)*az(iat)
        yimg = opsym(iop,2,1)*ax(iat) + opsym(iop,2,2)*ay(iat) + opsym(iop,2,3)*az(iat)
        zimg = opsym(iop,3,1)*ax(iat) + opsym(iop,3,2)*ay(iat) + opsym(iop,3,3)*az(iat)
        ! Find this image atom or the closest approximation:
        found = .false.
        dmin = 1d30
        kmin = -1
        k = 1
        do while (.not.found .and. k.le.n)
          dd = (ximg-ax(k))**2+(yimg-ay(k))**2+(zimg-az(k))**2
          dd = sqrt(dd)
          if (dd .le. TOLdist) then
            found = .true.
            dmin = dd
            kmin = k
          else if (dd .lt. dmin) then
            dmin = dd
            kmin = k
          end if
          k = k + 1
        end do
        if (.not.found) then
          ntrouble = ntrouble + 1
          dtrouble = max(dtrouble, dmin)
          if (.not.cleaned(kmin)) trouble = .true.
        else if (.not.cleaned(kmin)) then
          cleaned(kmin) = .true.
          dtrouble = max(dtrouble, dmin)
          ax(kmin) = ximg
          ay(kmin) = yimg
          az(kmin) = zimg
        end if
      end do
    end do
  end do

  if (verbose .or. ntrouble.gt.0) write (uout,630) ntrouble, dtrouble

  ! Check that all atoms have been cleaned:
  if (trouble) then
    ntrouble = 0
    do i = 1,n
      if (.not.cleaned(i)) then
        ntrouble = ntrouble + 1
        if (verbose) write (uout,635) i
      end if
    end do
    if (ntrouble.gt.0) then
      write (uout,636) ntrouble
      call ferror ('symproject', 'error cleaning atoms', faterr)
    end if
  end if
 
  ! Write down the cleaned coordinates:
  if (verbose) write (uout,640)
  dproj = 0.0_rp
  do i = 1,n
    xx = ax(i)-axx(i)
    yy = ay(i)-ayy(i)
    zz = az(i)-azz(i)
    dproj = dproj + sqrt(xx*xx + yy*yy + zz*zz)
    if (verbose) write (uout,645) i, ax(i), ay(i), az(i), ax(i)-axx(i), ay(i)-ayy(i), az(i)-azz(i)
  end do
  dproj = dproj / n
  if (verbose) write (uout,650) dproj
    
600 format (/                                                           &
    1x, '++ SYMPROJECT: Project onto the symmetry elements atoms',      &
    1x, 'too close to them.')
605 format (                                                            &
    1x, 'ALG(symproject) Move atom', i4, ' to the origin.',             &
    1x, 'Distance:', 1p, e12.4, 0p, 'Original coordinates: ', 3f15.9)
610 format (                                                            &
    1x, 'ALG(symproject) Move atom', i4, ' to the rotation axis', i4,   &
    1x, 'Distance:', 1p, e12.4, 0p,                                     &
    1x, 'Old and new coordinates: '/ (3f15.9))
615 format (                                                            &
    1x, 'ALG(symproject) Move atom', i4, ' to the mirror plane ', i4,   &
    1x, 'Distance:', 1p, e12.4, 0p,                                     &
    1x, 'Old and new coordinates: '/ (3f15.9))
630 format (                                                            &
    1x, 'ALG(symproject): Atomic images not found and max error ->',    &
    1x, i6, 1p, e12.4)
635 format (                                                            &
    1x, 'PROBLEM(symproject): Atom ', i4, ' was not cleaned!')
636 format (                                                            &
    1x, 'PANIC(symproject): Number of atoms not cleaned: ', i4)
640 format (/                                                           &
    1x, 'ALG(symproject): Cleaned molecular coordinates.'/              &
    1x, '-i------new-x----------new-y----------new-z------',            &
            '---dif-x-----dif-y-----dif-z--')
645 format (i4, 3f15.9, 1p, 3e10.2)
650 format (/                                                           &
    1x, 'ALG(symproject): Average displacement of atoms: ', 1p,e15.6/   &
    1x, 'ALG(symproject): end of cleaning ok!')
 
end subroutine
