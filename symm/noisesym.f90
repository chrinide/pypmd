! Estimate the average deviation of atoms with respect
! to the perfectly symmetric structure.
subroutine noisesym ()

  use iso_fortran_env, only: uout=>output_unit
  use mod_io, only: ferror, faterr
  use mod_prec, only: ip, rp
  use mod_param, only: verbose, pi
  use mod_sym, only: nopsym, norbit, opsym, opt_identity, opm, &
                     natorb, toldist, orbdis, oporder, opangle, &
                     optype, iatorb, ax, ay, az, n=>sdim
  implicit none
 
  integer(kind=ip) :: i, j, j1, k, iat, iorb, iop, kbest, kat
  real(kind=rp) :: R1, R2, R3
  real(kind=rp) :: dorb, dj, alpha, xxx, yyy, zzz, dbest, tol2
  logical :: done

  ! All atoms in a orbit should have the same distance to the
  ! center of mass. Get a simple noise measurement from checking
  ! this assumption:
  R1 = 0.0_rp
  R2 = 0.0_rp
  R3 = 0.0_rp
  k = 0
  do i = 1,norbit
    dorb = orbdis(i)
    do j = 1,natorb(i)
      k = k + 1
      j1 = iatorb(i,j)
      dj = sqrt(ax(j1)*ax(j1) + ay(j1)*ay(j1) + az(j1)*az(j1))
      R1 = R1 + (dj - dorb)
      R2 = R2 + (dj - dorb)**2
      R3 = R3 + abs(dj - dorb)
    end do
  end do
  R1 = R1 / k
  R2 = R2 / k
  R3 = R3 / k
  if (verbose) write (uout,605) R3, sqrt(R2 - R1*R1)

  ! Are the rotation angles of symmetry matrices exact?
  R1 = 0.0_rp
  R2 = 0.0_rp
  R3 = 0.0_rp
  do i = 1,nopsym
    alpha = opm(i)*2.0_rp*pi/oporder(i)
    if (opangle(i) .lt. 0.0_rp) alpha = -alpha
    R1 = R1 + (alpha - opangle(i))
    R2 = R2 + (alpha - opangle(i))**2
    R3 = R3 + abs(alpha - opangle(i))
  end do
  R1 = R1 / nopsym
  R2 = R2 / nopsym
  R3 = R3 / nopsym
  if (verbose) write (uout,610) R3, sqrt(R2 - R1*R1)

  ! A desperate movement:
  !   Produce all copies of all atoms operated by all symmetry
  !   matrices. Get the average distance between every atom and all
  !   of its copies:
  TOL2 = TOLdist * TOLdist
  R1 = 0.0_rp
  kat = 0
  do iop = 1,nopsym
    if (optype(iop) .ne. opt_identity) then
      do iorb = 1,norbit
        do i = 1,natorb(iorb)
          iat = iatorb(iorb,i)
          xxx = opsym(iop,1,1)*ax(iat) + opsym(iop,1,2)*ay(iat) &
                + opsym(iop,1,3)*az(iat)
          yyy = opsym(iop,2,1)*ax(iat) + opsym(iop,2,2)*ay(iat) &
                + opsym(iop,2,3)*az(iat)
          zzz = opsym(iop,3,1)*ax(iat) + opsym(iop,3,2)*ay(iat) &
                + opsym(iop,3,3)*az(iat)
          kbest = 0
          dbest = 1d30
          j = 1
          done = .false.
          do while (.not.done .and. j.le.natorb(iorb))
            k = iatorb(iorb,j)
            dj = (xxx-ax(k))**2+(yyy-ay(k))**2+(zzz-az(k))**2
            if (dj.le.dbest) then
              kbest = k
              dbest = dj
            end if
            if (dj .le. TOL2) then
              done = .true.
            else
              j = j + 1
            end if
          end do
          if (done) then
            kat = kat + 1
            R1 = R1 + sqrt(dj)
          else
            if (verbose) write (uout,650) iat, iop, kbest, dbest
          end if
        end do
      end do
    end if
  end do
  R1 = R1 / kat
  if (verbose) write (uout,615) R1, n*(nopsym-1)-kat

605 format (/                &
    1x, '++SYMNOISE:'/       &
    1x, 'Error on the orbit radii. Average and sdev: ', 1p, 2e15.6)

610 format (1x, 'Error on the rot. angles. Average and sdev: ', 1p, 2e15.6)

615 format (1x, 'Distance between sym/at images. Average reject: ', 1p, e15.6,0p, i8)

650 format (1x, 'DBG(symnoise): sym/at image not found', 3i5, 1p, e15.6)

end subroutine noisesym
