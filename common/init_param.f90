  subroutine init_param ()

    use mod_param, only: scratch, verbose, debug, &
                         isdata, d1mach, d1mach_

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

  end subroutine init_param
