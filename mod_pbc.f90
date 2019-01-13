module mod_pbc

  use mod_prec, only: rp, ip
  implicit none
  private

  real(kind=rp), dimension(3,3), public :: lattice
  integer(kind=ip), parameter, public :: nkpoints = 1
  real(kind=rp), dimension(3), parameter, public :: kpoints = 0.0
  integer(kind=ip), public :: ntvectors
  real(kind=rp), allocatable, dimension(:,:), public :: tvectors

end module mod_pbc
