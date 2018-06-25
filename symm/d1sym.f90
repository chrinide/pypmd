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
! check symmetry of linear molecules.
!
subroutine d1sym()

  use iso_fortran_env, only: uout=>output_unit
  use mod_io, only: ferror, faterr
  use mod_prec, only: ip, rp
  use mod_param, only: debug, verbose, pi
  use mod_math, only: matfill, inv
  use mod_sym, only: nmol=>sdim, ax, ay, az, inf_order, opmainord, &
                     linear_mol, nopsym, point_group, point_g2
  implicit none

  real(kind=rp), dimension(3,3) :: xmat
  real(kind=rp) :: xx, yy, zz, xnorm, alfa, theta, phi
  real(kind=rp) ::  sa, ca, st, ct, sp, cp, s2t, c2t, s2p, c2p
  integer(kind=ip) :: i, j, norder

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
 
  linear_mol = .true.
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
    point_group = "D_inf_h"
  else
    point_group = "C_inf_v"
  end if
  point_g2 = point_group
  opmainord = 999

  ! I don't know how to manage an infinity number of rotation
  ! operations. I will enforce a Dnh or Cnv group instead with n=4,
  ! for instance.
  norder = inf_order
  
  ! Determine the molecular axis and transform the direction into
  ! spherical polar angles:
  xnorm = 0.0_rp
  do i = 1,nmol
    if (ax(i)*ax(i) + ay(i)*ay(i) + az(i)*az(i) .gt. xnorm) then
      xx = ax(i)
      yy = ay(i)
      zz = az(i)
      xnorm = xx*xx + yy*yy + zz*zz
    end if
  end do
  if (xnorm .le. 0.0_rp) then
    call ferror('sym1d', 'Molecular axis not found!', faterr)
  end if
  xnorm = sqrt(xnorm)
  theta = acos(zz/xnorm)
  phi = atan2(yy, xx)
  if (verbose) then
    write(uout,600) xx/xnorm, yy/xnorm, zz/xnorm, &
                    theta*180.0_rp/pi, phi*180.0_rp/pi
    !write (uout,605) point_group(1:leng(point_group)), norder
    write (uout,605) point_group, norder
  end if
  
  ! Compose the C_n^1 rotation matrix around the molecular axis:
  alfa = (pi+pi)/norder
  ca = cos(alfa)
  sa = sin(alfa)
  ct = cos(theta)
  st = sin(theta)
  xmat(1,1) = ca*ct*ct + st*st
  xmat(1,2) = -sa*ct
  xmat(1,3) = (1-ca)*st*ct
  xmat(2,1) = sa*ct
  xmat(2,2) = ca
  xmat(2,3) = -sa*st
  xmat(3,1) = (1-ca)*st*ct
  xmat(3,2) = sa*st
  xmat(3,3) = ct*ct + ca*st*st
  if (opchksym(xmat)) then
    call symopadd (xmat)
  else
    if (debug) write (uout,610) ((xmat(i,j),j=1,3),i=1,3)
    call ferror ('sym1d', 'Algorithm error: C_n^1 not sym op!', faterr)
  end if
   
  ! Compose a sigma_v operation and test it:
  cp = cos(phi)
  sp = sin(phi)
  c2t = cos(theta+theta)
  s2t = sin(theta+theta)
  c2p = cos(phi+phi)
  s2p = sin(phi+phi)
  xmat(1,1) = ct*ct*c2p + st*st
  xmat(1,2) = ct*s2p
  xmat(1,3) = s2t*sp*sp
  xmat(2,1) = ct*s2p
  xmat(2,2) = -c2p
  xmat(2,3) = -2.0_rp*st*sp*cp
  xmat(3,1) = s2t*sp*sp
  xmat(3,2) = -2.0_rp*st*sp*cp
  xmat(3,3) = ct*ct + st*st*c2p
  if (opchksym(xmat)) then
    call symopadd (xmat)
  end if
 
  ! Check the closure of the symmetry operators set:
  call closuresym ()
  call getgroup ()
 
600 format (/                                                    &
    1x, '++SYM1D: Symmetry of a linear molecule'/                &
    1x, 'Molecular axis direction:       ', 3f15.6/              &
    1x, 'Direction in polar coordinates: ', 2f15.6, ' degrees')
605 format (                                                     &
    1x, 'True symmetry group: ', a/                              &
    1x, 'The infinite order axis will be changed to order: ', i4)
610 format (/                                                    &
    1x, 'ALG(sym1d) Rotation around the molecular axis not',     &
    1x, 'recognized as a symmetry operation:'/                   &
    (1x, 3f15.6))
 
end subroutine 
