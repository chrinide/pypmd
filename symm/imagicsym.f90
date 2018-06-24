! Check the elements in the symmetry matrices against
! the icosahedral magic numbers internally known by the routine.
subroutine imagicsym ()

  use iso_fortran_env, only: uout=>output_unit
  use mod_io, only: ferror, faterr
  use mod_prec, only: i4=>ip, rp
  use mod_param, only: debug, verbose
  use mod_sym, only: optable, opsym, MOPSYM, nopsym, opisgener, &
                     toleqvm, opt_inversion, opt_identity, optype
  implicit none
 
  real(kind=rp) :: x, ax, xbest, xdif, TOLxxx, xmat(3,3) 
  integer(kind=i4) :: iop, i, j, k, imin, imax, imed
  integer(kind=i4) :: kk, kbest, iter, ip, iq, ir
  logical :: found, doagain, generated(MOPSYM)
 
  integer(kind=i4) :: nmagic
  real(kind=rp) :: magic(24)

  data nmagic /24/
  data (magic(i), i = 1,24) /0.00000000000000000000_rp,   &
  0.05278640450004206072_rp, 0.13819660112501051518_rp,   &
  0.16245984811645316308_rp, 0.26286555605956680302_rp,   &
  0.27639320225002103036_rp, 0.30901699437494742410_rp,   &
  0.36180339887498948482_rp, 0.42532540417601996609_rp,   &
  0.44721359549995793928_rp, 0.50000000000000000000_rp,   &
  0.52573111211913360603_rp, 0.58778525229247312917_rp,   &
  0.63819660112501051518_rp, 0.67082039324993690892_rp,   &
  0.68819096023558676910_rp, 0.72360679774997896964_rp,   &
  0.80901699437494742410_rp, 0.85065080835203993218_rp,   &
  0.86180339887498948482_rp, 0.89442719099991587856_rp,   &
  0.94721359549995793928_rp, 0.95105651629515357212_rp,   &
  1.00000000000000000000_rp /

  ! 0.003842921 is the smallest difference between two magic numbers.
  ! Use a larger tolerance for the comparisons.
  TOLxxx = min(4d-3, TOLeqvm)

  ! Apply the magic cleaning to the group generators. Identity and
  ! inversion matrices are clean by construction.
  do iop = 1,nopsym
    if (optype(iop).eq.opt_identity .or. optype(iop).eq.opt_inversion) then
      generated(iop) = .true.
    else if (opisgener(iop)) then
      generated(iop) = .true.
      do i = 1,3
        do j = 1,3
          x = opsym(iop,i,j)
          ax = abs(x)
          imin = 1
          imax = nmagic
          k = 0
          found = .false.
          do while (.not.found)
            k = k + 1
            imed = (imin+imax)/2
            if (abs(ax-magic(imed)) .le. TOLxxx) then
              opsym(iop,i,j) = sign(magic(imed), x)
              found = .true.
            else if ((imax-imin).lt.1 .or. k.gt.20) then
              if (verbose) then
                write (uout,700) iop, i, j, ax, imed,imin,imax,k
              end if
              kbest = -1
              xbest = 1d30
              do kk = max(1,imed-3), min(nmagic,imed+3)
                xdif = abs(ax-magic(kk))
                if (xdif .lt. xbest) then
                  xbest = xdif
                  kbest = kk
                end if
              end do
              if (kbest .lt. 0) then
                write (uout,705) iop, i, j, x
                call ferror ('symimagic', 'Magic number not found!', faterr)
              end if
              opsym(iop,i,j) = sign(magic(kbest), x)
              if (verbose) then 
                write (uout,710) iop,i,j, x,opsym(iop,i,j),xbest
              end if
              found = .true.
            else if (magic(imed) .gt. ax) then
              imax = imed - 1
            else
              imin = imed + 1
            end if
          end do
        end do
      end do
    else
      generated(iop) = .false.
    end if
  end do

  ! Regenerate the symmetry operations from the cleaned generators:
  doagain = .true.
  iter = 0
  do while (doagain .and. iter.lt.20)
    doagain = .false.
    iter = iter + 1
    do i = 1,nopsym
      if (generated(i)) then
        do j = 1,nopsym
          k = optable(i,j)
          if (generated(j) .and. (.not.generated(k))) then
            ! This is a non regenerated opsym:
            doagain = .true.
            do ip = 1,3
              do iq = 1,3
                xmat(ip,iq) = 0.0_rp
                do ir = 1,3
                  xmat(ip,iq) = xmat(ip,iq) + opsym(i,ip,ir)*opsym(j,ir,iq)
                end do
              end do
            end do
            if (debug) then
              write (uout,910) k, ((opsym(k,ip,iq), ip=1,3), iq=1,3), &
                               ((xmat(ip,iq), ip=1,3), iq=1,3)
            end if
            do ip = 1,3
              do iq = 1,3
                opsym(k,ip,iq) = xmat(ip,iq)
              end do
            end do
            generated(k) = .true.
          end if
        end do
      end if
    end do
  end do
  if (doagain) then
    call error ('symimagic','the group was not regenerated in 20 iterations!', faterr)
  end if

  !.....Regenerate the symmetry matrices properties:
  call matfill (xmat,1.0_rp,0.0_rp,0.0_rp, &
                     0.0_rp,1.0_rp,0.0_rp, &
                     0.0_rp,0.0_rp,1.0_rp)
  call transformsym (xmat, xmat)
 
700 format (/                                                           &
    1x, '**SYMIMAGIC: Problem finding an icosahedral magic number'/     &
    1x, 'Debug data: ', 3i5, 1p, e17.9, 0p, 4i5/                        &
    1x, 'TROUBLE! The binary search has failed. We will use the',       &
    1x, 'magic number closest to the dirty unknown!')
705 format (                                                            &
    1x, 'PANIC! The best approach search has failed too:',              &
    1x, 'op(', i3, ',', i1, ',', i1, ') dirty value is ', 1p, e17.9)
710 format (                                                            &
    1x, 'Best approach result: op(', i3, ',', i1, ',', i1, ')',         &
    1x, 'dirty, clean, and diff. values ', 1p, 2e17.9, e11.3)
910 format (/                                                           &
    1x, 'DBG(symimagic): Regenerating sym matrix ', i5/                 &
    (1x, 3f15.9))
 
end subroutine imagicsym
