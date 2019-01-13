module mod_parallel

  !$ use omp_lib
  implicit none
  private
      
  character(len=132) :: scratch
  integer, public :: ichunk, nthreads, itile

  public :: init_parallel, end_parallel
  public :: info_parallel

contains

  subroutine init_parallel() 
    implicit none
    itile = 1
    nthreads = 1
    ichunk = 1
    scratch = './'
!$omp parallel
!$omp master
!$  nthreads=omp_get_num_threads()
!$omp end master
!$omp end parallel
  end subroutine init_parallel

  subroutine end_parallel()
    implicit none
  end subroutine end_parallel

  subroutine info_parallel(lu)
    implicit none
    integer, intent(in) :: lu
    write (lu,'(a,i0)') " # Number of OMP threads : ", nthreads
  end subroutine info_parallel

end module mod_parallel      
