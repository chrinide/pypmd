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
! print the symmetry operators.
! 
subroutine infosym ()

  use iso_fortran_env, only: uout=>output_unit
  use mod_prec, only: ip, rp
  use mod_param, only: pi
  use mod_math, only: matfill
  use mod_sym, only: nopsym, opaxis, opclas, opeuler, opinv, &
                     point_group, opsymbol, opisgener, opsym, &
                     optable, MOPSYM, nclas, ninclas, iinclas, &
                     opangle
  implicit none

  integer(kind=ip) :: i, j, k, itrouble, is, js 
  real(kind=rp) :: ang(MOPSYM), xprod, xtrouble, xmat(3,3)
  real(kind=rp) :: ovr(MOPSYM), ovlp, idn(3,3) 

  write (uout,700) point_group, nopsym
  do i = 1,nopsym
    write (uout,705) i, opsymbol(i), opangle(i)*180.0_rp/pi, &  
                       (opaxis(i,k), k=1,3), opisgener(i), opinv(i), &
                       opclas(i), (opeuler(i,k)*180.0_rp/pi, k=1,3), &
                       ((opsym(i,j,k), k=1,3), j=1,3)
  end do
  write (uout,780) nclas
  do i = 1,nclas
    write (uout,785) i, ninclas(i), opsymbol(iinclas(i,1)), &
                      (iinclas(i,j), j=1,ninclas(i))
  end do
  write (uout,750)
  do i = 1,nopsym
    write (uout,755) (optable(i,j), j=1,nopsym)
  end do

  itrouble = 0
  xtrouble = 0.0_rp
  write (uout,760) (i, opsymbol(i), i=1,nopsym)
  do i = 1,nopsym
    do j = 1,nopsym
      xprod = opaxis(i,1)*opaxis(j,1) + opaxis(i,2)*opaxis(j,2) &
            + opaxis(i,3)*opaxis(j,3)
      if (abs(xprod) .gt. 1.0_rp) then
        itrouble = itrouble + 1
        xtrouble = max(xtrouble, abs(xprod)-1.0_rp)
        xprod = sign(1.0_rp, xprod)
      end if
      ang(j) = acos(xprod)*180.0_rp/pi
    end do
    write (uout,765) i, opsymbol(i), (ang(j), j=1,nopsym)
  end do
  if (itrouble .gt. 0) write (uout,770) itrouble, xtrouble


  call matfill (idn, 1.0_rp,0.0_rp,0.0_rp, &
                     0.0_rp,1.0_rp,0.0_rp, &
                     0.0_rp,0.0_rp,1.0_rp)
  do is = 1,nopsym
    do js = 1,nopsym
      ovlp = 0.0_rp
      do i = 1,3
        do j = 1,3
          xmat(i,j) = 0.0_rp
          do k = 1,3
            xmat(i,j) = xmat(i,j) + opsym(is,k,i)*opsym(js,k,j)
          end do !k
          ovlp = ovlp + abs(xmat(i,j) - idn(i,j))
        end do !j
      end do !i
      write (uout,900) is, js, ((xmat(i,j),j=1,3),i=1,3)
      ovr(js) = ovlp
    end do !js
    write (uout,910) is, (ovr(j), j=1,nopsym)
  end do !is   _

760 format (/                                   &
    1x, 'Angles (deg) among symop directions'/  &
    9x, 120(i6, 3x, a5, 1x))
765 format (i4, a5, 120f9.3)
770 format (/                                                  &
    1x, 'Deviations from the [-1,+1] range (number, max dev)', i9, 1p, e20.12)
900 format (/                                            &
    1x, 'ALG(symprint) Transp ', i3, ' por matriz ', i3/ &
    (1x, 3f15.9))
910 format (/                                                  &
    1x, 'ALG(symprint) Overlap between sym matrices: is=', i3/ &
    (1x, 6f12.6) )
700 format (/                                        &
    1x, '++SYMPRINT: List of symmetry operators:'/   &
    1x, 'Point group: ', a/                          &
    1x, 'Number of operators found: ', i5)
705 format (/                                        &
    1x, 'Operator #', i4, ' :'/                      &
    1x, 'Symbol: ', a/                               &
    1x, 'Rotation angle:   ', f15.9, ' degrees'/     &
    1x, 'Rotation axis:    ', 3f15.9/                &
    1x, 'Is a generator?   ', l6/                    &
    1x, 'Inverse operation:', i6/                    &
    1x, 'Operations class: ', i6/                    &
    1x, 'Euler zxz angles: ', 3f15.9, 3x, 'degrees'/ &
    1x, 'Operator 3x3 matrix (xyz representation):'/ &
    (1x, 3f15.9))
750 format (/                                        &
    1x, 'Multiplication (Cayley) table:')
755 format (120i4)
780 format (/                                        &
    1x, 'Number of symmetry operation classes:', i5/ &
    1x, 'Class...#op....Repres..Operations...')      
785 format (2i6, 2x, a10, (20i4))

end subroutine 
