module mod_pbc

  use mod_prec, only: rp, ip
  implicit none
  private

  real(kind=rp), dimension(3,3), public :: lattice
  integer(kind=ip), public :: nkpoints
  integer(kind=ip), public :: ntvectors
  real(kind=rp), allocatable, dimension(:,:), public :: kpoints
  real(kind=rp), allocatable, dimension(:), public :: kweights
  real(kind=rp), allocatable, dimension(:,:), public :: tvectors

end module mod_pbc
