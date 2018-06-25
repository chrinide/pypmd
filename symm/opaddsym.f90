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
! if xmat is yet unknown add it to the list of sym ops.
!
subroutine opaddsym (xmat)

  use iso_fortran_env, only: uout=>output_unit
  use mod_io, only: ferror, faterr, warning
  use mod_prec, only: ip, rp
  use mod_param, only: debug, pi
  use mod_sym, only: nopsym, opsym, toleqvm, mopsym, opeuler, &
                     tolisint, opproper, opm, opsymbol, optype, &
                     oporder, opangle, opaxis, opt_identity, &
                     opt_inversion, opt_imp_rotation, opt_rotation, &
                     opt_sigma, tolnull
  use mod_math, only: inv , euler
  implicit none
 
  real(kind=rp), intent(in) :: xmat(3,3)

  integer(kind=ip), parameter :: MOP = 20
  integer(kind=ip) :: i, j, k, ierr, ii, ibest 
  real(kind=rp) :: diff, xtrace, xdet, xnm, ang1, angle 
  real(kind=rp) :: err, besterr, xinv(3,3), xval(3), xvect(3,3), xeul(3)
  character(len=40) :: fmt
  logical :: found
 
  do i = 1,nopsym
    diff = 0.0_rp
    do j = 1,3
      do k = 1,3
        diff = diff + abs(xmat(j,k) - opsym(i,j,k))
      end do
    end do
    if (diff .le. TOLeqvm) return  ! This is not a new operator
  end do

  nopsym = nopsym + 1
  if (nopsym .gt. MOPSYM) then
    call ferror ('opaddsym','too many sym operators', faterr)
  end if

  opsym(nopsym,:,:) = xmat
  if (debug) write (uout,500) ((xmat(i,j),j=1,3),i=1,3)

  ! Get the properties of the symmetry matrix:
  xtrace = 0.0_rp
  do i = 1,3
    xtrace = xtrace + xmat(i,i)
  end do
  call inv (xmat, xinv, xdet, ierr)
  if (ierr .lt. 0) then
    call ferror ('opaddsym','operator has a singular matrix', faterr)
  end if

  ! Get the Euler angles for the rotation:
  call euler (xmat, xeul)
  opeuler(nopsym,:) = xeul(:)
 
  ! Is this a proper or improper rotation?
  if (abs(abs(xdet)-1d0) .gt. TOLisint) then
    call ferror ('opaddsym', 'determinant is not +-1', faterr)
  else if (xdet .gt. 0.0_rp) then
    opproper(nopsym) = .true.
  else
    opproper(nopsym) = .false.
  end if

  ! Get the rotation axis except for E and i, for which it is
  ! fully degenerated, else Get the rotation angle and the n 
  ! and m of the C_n^m or S_n^m symmetry operator:
  if (abs(abs(xtrace)-3.0_rp) .le. TOLisint) then
    if (xtrace .gt. 0.0_rp) then
      opsymbol(nopsym) = 'E'
      opaxis(nopsym,1) = 0d0
      opaxis(nopsym,2) = 0d0
      opaxis(nopsym,3) = 1d0
      oporder(nopsym) = 1
      opm(nopsym) = 1
      opangle(nopsym) = pi + pi
      optype(nopsym) = opt_identity
    else
      opsymbol(nopsym) = 'i'
      opaxis(nopsym,1) = 0d0
      opaxis(nopsym,2) = 0d0
      opaxis(nopsym,3) = 1d0
      oporder(nopsym) = 2
      opm(nopsym) = 1
      opangle(nopsym) = pi
      optype(nopsym) = opt_inversion
    end if
  else
    ang1 = 0.5_rp*(xtrace-xdet)
    if (abs(ang1) .gt. 1.0_rp) ang1 = sign(1.0_rp, ang1)
    angle = acos(ang1)
    if (abs(angle) .le. TOLnull) then
      xnm = 1.0_rp
      angle = pi + pi
    else
      xnm = 2.0_rp*pi/angle
    end if
    opangle(nopsym) = angle
    ii = 0
    found = .false.
    do while (.not. found .and. ii.le.MOP)
      ii = ii + 1
      found = abs(xnm*ii - nint(xnm*ii)) .le. TOLisint
    end do
    if (found) then
      oporder(nopsym) = nint(xnm*ii)
      opm(nopsym) = ii
    else
      call ferror ('opaddsym', 'using best approx order', warning)
      ibest = 1
      besterr = abs(xnm - nint(xnm))
      do ii = 2, MOP
        err = abs(xnm*ii - nint(xnm*ii))
        if (err .le. besterr) then
          besterr = err
          ibest = ii
        end if
      end do
      oporder(nopsym) = nint(xnm*ibest)
      opm(nopsym) = ibest
    end if
    write (fmt,200) 1+int(log10(1d0*oporder(nopsym))), &
                    1+int(log10(1d0*opm(nopsym)))
    if (debug) write (uout,501) fmt
    call eigensym (xmat, xval, xvect, xdet, ierr)
    if (ierr .ne. 0) then
      call ferror ('opaddsym', 'trouble finding the rotation axis', faterr)
    endif
    opaxis(nopsym,1) = xvect(1,3)
    opaxis(nopsym,2) = xvect(2,3)
    opaxis(nopsym,3) = xvect(3,3)
    if (xdet .gt. 0.0_rp) then
      write (opsymbol(nopsym),fmt) 'C', oporder(nopsym), opm(nopsym)
      optype(nopsym) = opt_rotation
    else if (abs(xtrace - 1d0) .le. TOLisint) then
      write (opsymbol(nopsym),'(a)') 'sigma'
      optype(nopsym) = opt_sigma
    else if (oporder(nopsym).eq.1 .and. opm(nopsym).eq.1) then
      write (opsymbol(nopsym),'(a)') 'sigma'
      optype(nopsym) = opt_sigma
    else
      write (opsymbol(nopsym),fmt) 'S', oporder(nopsym), opm(nopsym)
      optype(nopsym) = opt_imp_rotation
    end if
  end if
  
  if (debug) then
    write (uout,700) nopsym, opsymbol(nopsym), fmt, optype(nopsym), xdet, xtrace, &
                     oporder(nopsym), opm(nopsym), ierr, ((xmat(i,j),j=1,3),i=1,3)
  end if

700 format (/                                                  &
    1x, 'DBG(symopadd): Trial symmetry matrix number ', i6/    &
    1x, '   Symbol  ', a/                                      &
    1x, '   Format  ', a/                                      &
    1x, '   Op type ', i6/                                     &
    1x, '   Det and trace ', 2f15.9/                           &
    1x, '   Order   ', 2i6/                                    &
    1x, '   Ierror  ', i6/                                     &
    1x, '   Matrix: '/                                         &
    (1x, 3f15.9))
 
501 format (1x, 'DBG(symopadd) fmt: ', a)
200 format ('(a, "_", i', i2, ', "^", i', i2, ')')
500 format (1x, 'DBG(symopadd) matrix: '/ (3f15.6))

end subroutine opaddsym
