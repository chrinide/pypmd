! This routine test if the input MOs are orthogonal
subroutine isortho ()

  use iso_fortran_env, only: uout=>output_unit
  use mod_prec, only: rp, ip
  use mod_memory, only: alloc, free
  use mod_param, only: verbose, debug
  use mod_io, only: string, ferror, warning
  use mod_wfn, only: epsortho, nprims, nmo
  implicit none

  logical :: ortho 
  real(kind=rp) :: solap
  integer(kind=ip) :: i, j
  real(kind=rp), allocatable, dimension(:,:) :: sprim

  call alloc ('isorthowfn', 'sprim', sprim , nprims, nprims)
  call gtogto (sprim)
 
  if (verbose) write (uout,'(1x,a)') string('# Testing orthogonality of natural MOs')
  ortho = .true.
  do i = 1,nmo
    do j = 1,i
      if (i.eq.j) then
        if (abs(abs(solap)-1.0_rp) .gt. epsortho) then
          ortho = .false.
          if (debug) write (uout,222) i,solap
        end if
      else
        if (abs(solap) .gt. epsortho) then
          ortho = .false.
          if (debug) write (uout,223) i,j,solap
        end if
      end if
    end do 
  end do
  if (.not.ortho) then
    call ferror ('isorthowfn', 'the set of natural mos are not orthonormal', warning)
  end if

  if (verbose) write (uout,'(1x,a)') string('# Testing orthogonality of canonical MOs')
  ortho = .true.
  do i = 1,nmo
    do j = 1,i
      if (i.eq.j) then
        if (abs(abs(solap)-1.0_rp) .gt. epsortho) then
          ortho = .false.
          if (debug) write (uout,222) i,solap
        end if
      else
        if (abs(solap) .gt. epsortho) then
          ortho = .false.
          if (debug) write (uout,223) i,j,solap
        end if
      end if
    end do 
  end do
  if (.not.ortho) then
    call ferror ('isorthowfn', 'the set of canonical mos are not orthonormal', warning)
  end if

  call free ('isorthowfn', 'sprim', sprim)
 
 222  format (1x,'# MO number ',i0, ' is not exactly normalized, NORM = ', e22.16)
 223  format (1x,'# MOs ',i0,' and ',i0, ' are not exactly orthogonal, S = ', e22.16)

end subroutine isortho
