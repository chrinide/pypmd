subroutine infogeom ()

  use iso_fortran_env, only: uout=>output_unit
  use mod_prec, only: ip, rp
  use mod_io, only: string, mline
  use mod_geom, only: cm, cq, moi, emoi
  use mod_wfn, only: ncent, atnam, xyz, charge
  use mod_param, only: h, pi, na, c
  implicit none

  ! threshold for degenerate principal moments of inertia
  real(kind=rp), parameter :: mom_thresh = 1e-3

  real(kind=rp) :: mhz, cm1
  integer(kind=ip) :: i
  character(len=mline) :: moltype
  logical :: same12, same13, same23, onezero, allzero
  logical :: iszero, degen
  
  !cm = cm*bohrtoa
  !cq = cq*bohrtoa

  write (uout,'(1x,a,1x,3(1x,f0.6))') '# Center of mass', cm(1), cm(2), cm(3)
  write (uout,'(1x,a)') string('# Center of mass translated geometry')
  do i = 1,ncent
    write (uout,'(1x,a,1x,i0,1x,a,1x,f4.1,1x,3f12.6)') &
           '#', i, atnam(i)(1:2), charge(i), xyz(i,:) - cm(:)
  end do
  write (uout,'(1x,a,1x,3(1x,f0.6))') '# Center of charge', cq(1), cq(2), cq(3)
  write (uout,'(1x,a)') string('# Center of charge translated geometry')
  do i = 1,ncent
    write (uout,'(1x,a,1x,i0,1x,a,1x,f4.1,1x,3f12.6)') &
           '#', i, atnam(i)(1:2), charge(i), xyz(i,:) - cq(:)
  end do
  write (uout,'(1x,a)') '# Moment of inertia tensor (amu * A^2)'
  do i = 1,3
    write (uout,'(1x,a,1x,3(1x,f0.6))') '#', moi(i,:)
  end do
  write (uout,'(1x,a)') '# Principal moments of inertia (amu * A^2)'
  write (uout,'(1x,a,1x,3(1x,f0.6))') '#', emoi(:)

  same12 = are_same(emoi(1), emoi(2), mom_thresh)
  same13 = are_same(emoi(1), emoi(3), mom_thresh)
  same23 = are_same(emoi(2), emoi(3), mom_thresh)
  onezero = are_same(emoi(1), 0.0_rp, mom_thresh)
  allzero = are_same(emoi(3), 0.0_rp, mom_thresh)
  if (allzero) then
    moltype = 'monatomic'
  else if (onezero) then
    moltype = 'linear'
  else if (same13) then
    moltype = 'a spherical top'
  else if (same12 .or. same23) then
    moltype = 'a symmetric top'
  else
    moltype = 'an asymmetric top'
  end if
  write (uout,'(1x,a,1x,a)') '# The molecule is', string(moltype)

  !Calculate rotational frequencies in MHz and cm-1
  do i = 1,3
    iszero = are_same(emoi(i), 0.0_rp, mom_thresh)
    if (i>1) then
      degen = are_same(emoi(i-1), emoi(i), mom_thresh)
    end if
    if (iszero .or. degen) then
      continue
    end if
    mhz = h/(8.0_rp*pi**2*emoi(i))
    mhz = mhz * (1e10)**2 *na*1e-3
    cm1 = mhz / c*1e6
    write (uout, *) mhz, cm1
  end do

  call connect (1.2_rp)

contains

  logical function are_same(n1, n2, tol)

    use mod_prec, only: rp, ip
    implicit none
  
    real(kind=rp), intent(in) :: tol
    real(kind=rp), intent(in) :: n1, n2

    real(kind=rp) :: comp

    are_same = .false.
    comp = abs((n2 - n1))

    if (comp <= tol) then
      are_same = .true.
    end if

    return 

  end function

end subroutine infogeom
