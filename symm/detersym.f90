! Determine linear, planar, 3d
subroutine detersym()

  use iso_fortran_env, only: uout=>output_unit
  use mod_io, only: ferror, faterr
  use mod_prec, only: ip, rp
  use mod_param, only: debug, verbose
  use mod_sym, only: mol_linear, mol_planar, toldist, sdim, &
                     ax, ay, az
  implicit none

  integer(kind=ip)  :: i, j, k, n
  real(kind=rp) :: an(sdim), nor, sum1, sum2, p_esc, p_mixto, contri

  n = sdim

  if (verbose) write (uout,800)
  do i = 1,n
    nor = sqrt(ax(i)*ax(i) + ay(i)*ay(i) + az(i)*az(i))
    if (nor > 1d-9) then
      an(i) = 1.0_rp/nor
    else
      an(i) = 1.0_rp
    end if
  end do

  ! Test for linear molecules:
  ! The normalized scalar product of the position vectors of each pair
  ! of atoms must be +1 or -1.
  sum1 = 0.0_rp
  do i = 1,n
    do j = i+1,n
      p_esc = (ax(i)*ax(j)+ay(i)*ay(j)+az(i)*az(j))*an(i)*an(j)
      contri = abs(abs(p_esc) - 1.0_rp)
      sum1 = sum1 + contri
      if (debug) then
        if (abs(contri) .ge. TOLdist) then
          write (uout,300) i, j, p_mixto
        end if
      end if
    end do
  end do
  if (verbose) write (uout,820) sum1, TOLdist
  if (sum1 .le. TOLdist) then
    mol_linear = .true.
    if (verbose) write (uout,810)
  end if

  ! Test for planar molecules:
  ! The triple (volume) product of the position vectors of each group
  ! of three different atoms must be zero.
  sum2 = 0.0_rp
  do i = 1,n
    do j = i+1,n
      do k = j+1,n
        p_mixto = ax(i)*(ay(j)*az(k) - az(j)*ay(k)) &
                + ay(i)*(az(j)*ax(k) - ax(j)*az(k)) &
                + az(i)*(ax(j)*ay(k) - ay(j)*ax(k))
        sum2 = sum2 + p_mixto
        if (debug) then
          if (abs(p_mixto) .ge. TOLdist) then
            write (uout,310) i, j, k, p_mixto
          end if
        end if
      end do
    end do
  end do
  if (verbose) write (uout,821) sum2, TOLdist
  if (sum2 .le. TOLdist) then
    mol_planar = .true.
    if (verbose) write (uout,815)
  end if

 300 format (' DBG(molsym): non-linear pair ', 2i4, ' contri ', e12.5)
 310 format (' DBG(molsym): non-planar trio ', 3i4, ' p_mixto ', e12.5)
 800 format (/                                             &
     1x, '++MOLSYM: Determining the point group symmetry') 
 810 format (                                              &
     1x, 'Linear molecule detected')                       
 815 format (                                              &
     1x, 'Planar molecule detected')
 820 format (                                              &
     1x, 'Linearity test value (TOL): ', 1p, e10.3, '   (',e10.3,')')
 821 format (                                              &
     1x, 'Planarity test value (TOL): ', 1p, e10.3, '   (',e10.3,')')

end subroutine detersym
