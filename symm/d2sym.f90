! check symmetry of planar molecules.
subroutine d2sym()

  use iso_fortran_env, only: uout=>output_unit
  use mod_io, only: ferror, faterr
  use mod_prec, only: ip, rp
  use mod_param, only: debug, verbose
  use mod_math, only: matfill, inv
  use mod_sym, only: mol_linear, mol_planar, toldist, nmol=>sdim, &
                     ax, ay, az, atzmol, linear_mol, errsng, nopsym, &
                     toldirty, tolsng, orbmol
  implicit none
 
  real(kind=rp), dimension(3,3) :: v, vinv, xmat, xm, xv, xop
  real(kind=rp) :: vdet, xx, yy, zz, xnorm, ztest, xdet
  integer :: i, j, ierr, ii, jj, i1, i2, j1, j2

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

  ! Get the new coordinates with Z being perpendicular to the
  ! molecular plane:
  call matfill (v, 1.0_rp,0.0_rp,0.0_rp, &
                   0.0_rp,1.0_rp,0.0_rp, &
                   0.0_rp,0.0_rp,1.0_rp)
  call matfill (vinv, 1.0_rp,0.0_rp,0.0_rp, &
                      0.0_rp,1.0_rp,0.0_rp, &
                      0.0_rp,0.0_rp,1.0_rp)
  do i = 1,nmol
    do j = i+1,nmol
      xx = ay(i)*az(j) - az(i)*ay(j)
      yy = az(i)*ax(j) - ax(i)*az(j)
      zz = ax(i)*ay(j) - ay(i)*ax(j)
      xnorm = xx*xx + yy*yy + zz*zz
      if (xnorm .gt. 1d-1) goto 1100
    end do
  end do
  call ferror ('2dsym', 'Normal to the molecular plane not found!', faterr)
1100 continue

  xnorm = 1.0/sqrt(xnorm)
  v(1,3) = xx*xnorm
  v(2,3) = yy*xnorm
  v(3,3) = zz*xnorm

  if (debug) write (uout,560) v(1,3), v(2,3), v(3,3)

  ! Check if the molecule is already in the XY plane and no
  ! reorientation is needed:
  if (abs(v(1,3)).le.TOLdist .and. abs(v(2,3)).le.TOLdist .and. &
      abs(v(3,3)-1_rp).le.TOLdist) goto 1103

  ! Otherwise get the two  directions of the molecular plane:
  ! (they will be converted into the new XY axes)
  do i = 1,nmol
    xx = ax(i)
    yy = ay(i)
    zz = az(i)
    xnorm = xx*xx + yy*yy + zz*zz
    if (xnorm .gt. 1d-2) goto 1101
  end do
  call ferror ('sym2d', 'New x axis not found!', faterr)
1101 continue

  xnorm = 1.0/sqrt(xnorm)
  v(1,1) = xx*xnorm
  v(2,1) = yy*xnorm
  v(3,1) = zz*xnorm

  ! The vectorial product v3 x v1 = v2 should be normalized:
  xx = v(2,3)*v(3,1) - v(3,3)*v(2,1)
  yy = v(3,3)*v(1,1) - v(1,3)*v(3,1)
  zz = v(1,3)*v(2,1) - v(2,3)*v(1,1)
  xnorm = xx*xx + yy*yy + zz*zz
  xnorm = 1.0/sqrt(xnorm)
  v(1,2) = xx*xnorm
  v(2,2) = yy*xnorm
  v(3,2) = zz*xnorm

  ! The V and V^{-1} matrices produce the change of coordinates and
  ! transform the symmetry matrices from the old to the new
  ! coordinates:
  call inv (v, vinv, vdet, ierr)
  if (ierr.ne.0) then
    call ferror ('sym2d','error in the inversion of the V matrix!', faterr)
  end if
  if (debug) then 
    write (uout,500) ((v(i,j),j=1,3),i=1,3), ((vinv(i,j),j=1,3),i=1,3)
  end if

  ztest = 0.0_rp
  do i = 1,nmol
    xx = vinv(1,1)*ax(i) + vinv(1,2)*ay(i) + vinv(1,3)*az(i)
    yy = vinv(2,1)*ax(i) + vinv(2,2)*ay(i) + vinv(2,3)*az(i)
    zz = vinv(3,1)*ax(i) + vinv(3,2)*ay(i) + vinv(3,3)*az(i)
    ztest = ztest + abs(zz)
    ax(i) = xx
    ay(i) = yy
    az(i) = zz
  end do
  if(ztest .gt. 1d-4) then
    write (uout,500) ((v(i,j),j=1,3),i=1,3),((vinv(i,j),j=1,3),i=1,3)
    write (uout,505) (i,ax(i),ay(i),az(i),i=1,nmol)
    call ferror ('sym2d','Transformed atoms do not fulfill z=0 test!', faterr)
  end if

  ! Start looking for symmetry operations:
1103 continue

  if (debug) then
    write (uout,500) ((v(i,j),j=1,3),i=1,3),((vinv(i,j),j=1,3),i=1,3)
  end if

  ! The identity is always a sym operator:
  linear_mol = .false.
  nopsym = 0
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
 
  ! Reflection in the molecular plane should be a symmetry operation:
  call matfill (xmat, 1.0_rp,0.0_rp,0.0_rp, &
                         0.0_rp,1.0_rp,0.0_rp, &
                         0.0_rp,0.0_rp,-1.0_rp)
  if (opchksym(xmat)) then
    call opaddsym (xmat)
  end if
 
  ! Run over all doblet pairs. Each pair migth produce a new symmetry
  ! operator. Only doblets having the same atoms and in the same
  ! order are compatible.
  ! The position matrix of the first doblet must be invertible for
  ! the algorithm to work.
  !
  do i1 = 1,nmol
    xmat(1,1) = ax(i1)
    xmat(2,1) = ay(i1)
    do i2 = i1+1,nmol
      xmat(1,2) = ax(i2)
      xmat(2,2) = ay(i2)
      xdet = xmat(1,1)*xmat(2,2) - xmat(1,2)*xmat(2,1)
      if (abs(xdet) .le. TOLsng) goto 1001
      xv(1,1) = xmat(2,2)/xdet
      xv(1,2) = -xmat(1,2)/xdet
      xv(2,1) = -xmat(2,1)/xdet
      xv(2,2) = xmat(1,1)/xdet
      ERRsng = max(ERRsng, abs(xdet))
      if (debug) then
        write (uout,510) i1, i2
        write (uout,515) 'xmat', ((xmat(ii,jj),jj=1,2),ii=1,2)
        write (uout,515) 'xinv', ((xv(ii,jj),jj=1,2),ii=1,2)
      end if
      do j1 = 1,nmol
        if (orbmol(i1).ne.orbmol(j1)) goto 1002
        xm(1,1) = ax(j1)
        xm(2,1) = ay(j1)
        do j2 = 1,nmol
          if (orbmol(i2).ne.orbmol(j2)) goto 1003
          xm(1,2) = ax(j2)
          xm(2,2) = ay(j2)
          xop(1,1) = xm(1,1)*xv(1,1) + xm(1,2)*xv(2,1)
          xop(1,2) = xm(1,1)*xv(1,2) + xm(1,2)*xv(2,2)
          xop(2,1) = xm(2,1)*xv(1,1) + xm(2,2)*xv(2,1)
          xop(2,2) = xm(2,1)*xv(1,2) + xm(2,2)*xv(2,2)
          xop(1,3) = 0_rp
          xop(2,3) = 0_rp
          xop(3,1) = 0_rp
          xop(3,2) = 0_rp
          xop(3,3) = 1_rp
          if (debug) then
            write (uout,520) i1,i2, j1,j2
            write (uout,524) 'xmati', ((xmat(ii,jj),jj=1,2),ii=1,2)
            write (uout,524) 'xinvi', ((xv(ii,jj),jj=1,2),ii=1,2)
            write (uout,524) 'xmatj', ((xm(ii,jj),jj=1,2),ii=1,2)
            write (uout,525) 'xop', ((xop(ii,jj),jj=1,3),ii=1,3)
          end if
          ! Check if this is a new sym operator:
          if (opchksym(xmat)) then
            call opaddsym (xop)
            if (TOLdirty) call closuresym ()
            if (debug) then
              write (uout,534) 'xmati', ((xmat(ii,jj),jj=1,2),ii=1,2)
              write (uout,534) 'xinvi', ((xv(ii,jj),jj=1,2),ii=1,2)
              write (uout,534) 'xmatj', ((xm(ii,jj),jj=1,2),ii=1,2)
              write (uout,535) 'xop', ((xop(ii,jj),jj=1,3),ii=1,3)
            end if
        end if
1003      continue
        end do
1002    continue
      end do
1001  continue
    end do
  end do

  ! Transform back the coordinates and the symmetry operations to the
  ! original orientation:
  do i = 1,nmol
    xx = v(1,1)*ax(i) + v(1,2)*ay(i) + v(1,3)*az(i)
    yy = v(2,1)*ax(i) + v(2,2)*ay(i) + v(2,3)*az(i)
    zz = v(3,1)*ax(i) + v(3,2)*ay(i) + v(3,3)*az(i)
    ax(i) = xx
    ay(i) = yy
    az(i) = zz
  end do
  call transformsym (v, vinv)

  ! Check the closure of the symmetry operators set:
  call closuresym ()
  call getgroup ()


560 format (1x, 'DBG(sym2d) New z direction: ', 3f15.9)
500 format (/                                                   &
    1x, 'DBG(sym2d) Planar molecule transformation problem: '/  &
    1x, 'Transformation matrix and inverse: '/                  &
    (1x, 3f16.9))
520 format (1x, 'DBG(sym2d) Doblet pair: ', 6i5)
524 format (1x, 'DBG(sym2d) ', a, ' matrix: '/ (2f15.6))
525 format (1x, 'DBG(sym2d) ', a, ' matrix: '/ (3f15.6))
510 format (1x, 'DBG(sym2d) Doblet: ', 3i5)
515 format (1x, 'DBG(sym2d) ', a, ' matrix: '/ (2f15.6))
534 format (1x, 'DBG(sym2d) ', a, ' matrix: '/ (2f15.6))
535 format (1x, 'DBG(sym2d) ', a, ' matrix: '/ (3f15.6))
505 format (                                                   &
    1x, 'DBG(sym2d) Transformed molecular coordinates:'/       &
    1x, i5, 3f16.9)

end subroutine d2sym
