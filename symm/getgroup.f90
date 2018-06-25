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
! determine the point group from the set of symmetry
! operators.
!
subroutine getgroup ()
 
  use iso_fortran_env, only: uout=>output_unit
  use mod_io, only: ferror, faterr
  use mod_prec, only: ip, rp
  use mod_param, only: debug
  use mod_sym, only: opt_rotation, nopsym, optype, opm, opt_sigma, &
                     opsymbol, oporder, opmainord, opt_imp_rotation, &
                     point_group, point_g2, tolisint, tolnull, opaxis, &
                     opt_inversion 
  implicit none

  real(kind=rp) :: xprod
  integer(kind=ip), parameter :: MOP = 20
  integer(kind=ip) :: mainorder, mainmult, imain(MOP), mainimproper
  integer(kind=ip) :: binmult, ibin(MOP), nsigma, isigma(MOP), isigmah
  integer(kind=ip) :: nbinperp, i, j, j1
  character(len=4) :: chmain, chimp
  logical :: inversion, sigmah
 
  ! Get the main axis, if any, and their multiplicity:
  ! Only the C_n^1 operations will be considered.
  ! Unfortunately, C_n^1 and C_n^(n-1) operations cannot be
  ! distinguished.
  ! Get also the list of binary axes, the sigma planes, and the
  ! order of the main improper axis (n > 2).
  mainorder = 0
  mainmult = 0
  mainimproper = 0
  binmult = 0
  inversion = .false.
  nsigma = 0
  do i = 1,nopsym
    if (optype(i) .eq. opt_rotation .and. opm(i) .eq. 1) then
      if (oporder(i) .gt. mainorder) then
        mainorder = oporder(i)
        mainmult = 1
        imain(mainmult) = i
      else if (oporder(i) .eq. mainorder) then
        mainmult = mainmult + 1
        if (mainmult .gt. MOP) call ferror ('getgroup','main axis mult too high', faterr)
        imain(mainmult) = i
      end if
      if (oporder(i) .eq. 2) then
        binmult = binmult + 1
        if (binmult .gt. MOP) call ferror ('getgroup','binary axis mult too high', faterr)
        ibin(binmult) = i
      end if
    else if (optype(i) .eq. opt_inversion) then
      inversion = .true.
    else if (optype(i) .eq. opt_sigma) then
      nsigma = nsigma + 1
      if (nsigma .gt. mop) call ferror ('getgroup', 'too many sigma planes', faterr)
      isigma(nsigma) = i
    else if (optype(i) .eq. opt_imp_rotation) then
      if (oporder(i) .gt. mainimproper) mainimproper = oporder(i)
    end if
  end do
  if (mainorder .gt. 2) mainmult = mainmult/2

  ! If there is a single main order group, look for sigma_h planes:
  sigmah = .false.
  if (mainmult .eq. 1 .or. mainorder .eq. 2) then
    i = imain(1)
    j = 1
    do while (.not.sigmah .and. j .le. nsigma)
      j1 = isigma(j)
      xprod = opaxis(i,1)*opaxis(j1,1)  &
            + opaxis(i,2)*opaxis(j1,2)  &
            + opaxis(i,3)*opaxis(j1,3)
      if (abs(abs(xprod)-1d0) .le. TOLisint) then
        sigmah = .true.
        isigmah = j1
      end if
      j = j + 1
    end do
  end if
  if (sigmah) then
    opsymbol(isigmah) = 'sigma_h'
  end if
 
  ! If there is a single main order group, look for n C2 perpendicular
  ! axis:
  nbinperp = 0
  if (mainmult .eq. 1 .or. mainorder .eq. 2) then
    i = imain(1)
    do j = 1,binmult
      j1 = ibin(j)
      xprod = opaxis(i,1)*opaxis(j1,1) &
            + opaxis(i,2)*opaxis(j1,2) &
            + opaxis(i,3)*opaxis(j1,3)
      if (abs(xprod) .le. TOLnull) then
        nbinperp = nbinperp + 1
      end if
    end do
  end if
  if (debug) then
    write (uout,700) mainorder, mainmult, inversion, binmult, &
                     nbinperp, nsigma, sigmah, mainimproper
  end if

  ! Store as a character the main proper and improper orders:
  if (mainorder .gt. 999) then
    write (chmain,'(i4)') mainorder
  else if (mainorder .gt. 99) then
    write (chmain,'(i3)') mainorder
  else if (mainorder .gt. 9) then
    write (chmain,'(i2)') mainorder
  else
    write (chmain,'(i1)') mainorder
  end if
  if (mainimproper .gt. 999) then
    write (chimp,'(i4)') mainimproper
  else if (mainimproper .gt. 99) then
    write (chimp,'(i3)') mainimproper
  else if (mainimproper .gt. 9) then
    write (chimp,'(i2)') mainimproper
  else
    write (chimp,'(i1)') mainimproper
  end if
  opmainord = mainorder
 
  ! Decision tree to get the point group. Start with the cubic-like
  ! groups:
  if (mainmult .gt. 1 .and. mainorder .gt. 2) then
    if (mainorder .eq. 5 .and. mainmult .ge. 6) then
      if (inversion) then
        point_group = 'Ih'
      else
        point_group = 'I'
      endif
      point_g2 = point_group
      opmainord = 5
    else if (mainorder .eq. 4 .and. mainmult .ge. 3) then
      if (inversion) then
        point_group = 'Oh'
      else
        point_group = 'O'
      end if
      point_g2 = point_group
      opmainord = 4
    else if (mainorder .eq. 3 .and. mainmult .ge. 4) then
      if (inversion) then
        point_group = 'Th'
      else if (nsigma .eq. 6) then
        point_group = 'Td'
      else
        point_group = 'T'
      end if
      point_g2 = point_group
      opmainord = 3
    else
      if (debug) then
        write (uout,700) mainorder, mainmult, inversion, binmult, &
                         nbinperp, nsigma, sigmah, mainimproper
      end if
      call ferror ('getgroup', 'unknown cubic-like group', faterr)
    end if
  else if (mainorder .ge. 2) then
    if (mainorder .eq. nbinperp) then
      if (sigmah) then
        !write (point_group,500) 'Dnh', mainorder
        point_group = "D" // chmain(1:len_trim(chmain)) // "h"
        !point_group = "D" // chmain // "h"
        point_g2 = 'Dnh'
      else if (mainorder .eq. nsigma) then
        !write (point_group,500) 'Dnd', mainorder
        point_group = "D" // chmain(1:len_trim(chmain)) // "d"
        !point_group = "D" // chmain // "d"
        point_g2 = 'Dnd'
      else
        !write (point_group,500) 'Dn', mainorder
        point_group = "D" // chmain(1:len_trim(chmain))
        !point_group = "D" // chmain
        point_g2 = 'Dn'
      end if
    else
      if (sigmah) then
        !write (point_group,500) 'Cnh', mainorder
        point_group = "C" // chmain(1:len_trim(chmain)) // "h"
        !point_group = "C" // chmain // "h"
        point_g2 = 'Cnh'
      else if (mainorder .eq. nsigma) then
        !write (point_group,500) 'Cnv', mainorder
        point_group = "C" // chmain(1:len_trim(chmain)) // "v"
        !point_group = "C" // chmain // "v"
        point_g2 = 'Cnv'
      else if (2*mainorder .eq. mainimproper) then
        !write (point_group,500) 'S2n', mainorder
        point_group = "S" // chimp(1:len_trim(chimp))
        !point_group = "S" // chmain 
        point_g2 = 'S2n'
        opmainord = mainorder
        opmainord = mainimproper
      else
        !write (point_group,500) 'Cn', mainorder
        point_group = "C" // chmain(1:len_trim(chmain))
        !point_group = "C" // chmain 
        point_g2 = 'Cn'
      end if
    end if
  else if (nopsym .eq. 2) then
    if (nsigma .eq. 1) then
          point_group = 'Cs'
          point_g2 = 'Cs'
    else if (inversion) then
          point_group = 'Ci'
          point_g2 = 'Ci'
    else
      if (debug) then          
        write (uout,700) mainorder, mainmult, inversion, binmult, &
                         nbinperp, nsigma, sigmah, mainimproper
      end if
      call ferror ('getgroup', 'unknown group', faterr)
    end if
  else if (nopsym .eq. 1) then
    point_group = 'C1'
    point_g2 = 'C1'
  else
    if (debug) then          
      write (uout,700) mainorder, mainmult, inversion, binmult, &
                       nbinperp, nsigma, sigmah, mainimproper
    end if
    call ferror ('getgroup', 'unknown group', faterr)
  end if
 
!500 format (a, ', n = ', i2)
700 format (/                                               &
    1x, '++SYMGETGROUP:'/                                   &
    1x, 'DBG(symgetgroup) main order & mult.......:', 2i5/  &
    1x, 'DBG(symgetgroup) inversion...............:', l5/   &
    1x, 'DBG(symgetgroup) binary axes.............:', i5/   &
    1x, 'DBG(symgetgroup) binary axes perp to main:', i5/   &
    1x, 'DBG(symgetgroup) symmetry planes.........:', i5/   &
    1x, 'DBG(symgetgroup) sigma_h plane...........:', l5/   &
    1x, 'DBG(symgetgroup) main improper axis......:', i5)
 
end subroutine
