! This routine test if the input MOs are orthogonal
subroutine isorthowfn ()

  use iso_fortran_env, only: uout=>output_unit
  use mod_prec, only: rp, ip
  use mod_memory, only: alloc, free
  use mod_param, only: verbose, debug
  use mod_io, only: string, ferror, warning
  use mod_wfn, only: epsortho, nprims, nmo, coefcan, coefnat
  implicit none

  logical :: ortho 
  real(kind=rp) :: solap
  integer(kind=ip) :: i, j
  real(kind=rp), allocatable, dimension(:,:) :: sprim, overlap, coeftmp

  call alloc ('isorthowfn', 'sprim', sprim , nprims, nprims)
  call alloc ('isorthowfn', 'overlap', overlap , nmo, nmo)
  call alloc ('isorthowfn', 'coeftmp', coeftmp , nmo, nprims)
  call gtogto (sprim)
 
  if (verbose) write (uout,'(1x,a)') string('# Testing orthogonality of natural MOs')
  coeftmp = matmul(sprim,transpose(coefnat))
  overlap = matmul(coefnat,coeftmp)
  ortho = .true.
  do i = 1,nmo
    do j = 1,i
      solap = overlap(i,j)
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
  coeftmp = matmul(sprim,transpose(coefcan))
  overlap = matmul(coefcan,coeftmp)
  do i = 1,nmo
    do j = 1,i
      solap = overlap(i,j)
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
  call free ('isorthowfn', 'overlap', overlap)
  call free ('isorthowfn', 'coeftmp', coeftmp)
 
 222  format (1x,'# MO number ',i0, ' is not exactly normalized, NORM = ', e22.16)
 223  format (1x,'# MOs ',i0,' and ',i0, ' are not exactly orthogonal, S = ', e22.16)

end subroutine isorthowfn
