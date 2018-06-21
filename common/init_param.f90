subroutine init_param ()

  use iso_fortran_env, only: uout=>output_unit
  use mod_io, only: string
  use mod_prec, only: rp
  use mod_param, only: scratch, verbose, debug, vbig, &
                       isdata, d1mach, d1mach_, debug, &
                       vsmall, epsreal

  implicit none
  integer :: i, n, clock
  integer, dimension(:), allocatable :: seed

  isdata = .false.
  scratch = './'
  verbose = .false.
  debug = .false.

  ! random seed
  call random_seed(size=n)
  allocate(seed(n))
  call system_clock(count=clock)
  seed = clock + 37*(/(i-1, i =1,n)/)
  call random_seed(put=seed)
  deallocate(seed)

  ! machine constants
  do i = 1,5
    d1mach_(i) = d1mach(i)
  end do
  vsmall = d1mach_(1)
  vbig = d1mach_(2)
  epsreal = epsilon(1.0_rp)

  if (debug) then
    write (uout,'(1x,a,1x,e13.6)') string('# (DEBUG) smallest number :'), vsmall
    write (uout,'(1x,a,1x,e13.6)') string('# (DEBUG) biggest number :'), vbig
    write (uout,'(1x,a,1x,e13.6)') string('# (DEBUG) machine eps :'), epsreal
  end if

end subroutine init_param
