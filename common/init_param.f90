  subroutine init_param ()

    use mod_param, only: scratch, verbose, debug, isdata
    implicit none

    isdata = .false.
    scratch = './'
    verbose = .false.
    debug = .false.

  end subroutine init_param
