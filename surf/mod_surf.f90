module mod_surf

  use mod_prec, only: rp, ip
  implicit none
  private

  integer(kind=ip), parameter, public :: minter = 10 
  integer(kind=ip), parameter, public :: maxtrial = 13

  integer(kind=ip), public :: inuc
  real(kind=rp), public :: xnuc(3)
  real(kind=rp), allocatable, dimension(:,:), public :: rlimsurf
  integer(kind=ip), allocatable, dimension(:), public :: nlimsurf

  ! nuclear critical points
  real(kind=rp), allocatable, dimension(:,:), public :: xyzrho

  ! options
  integer(kind=ip), dimension(4), public :: nangleb
  ! Precision in interatomic surfaces. 
  real(kind=rp), public :: rmaxsurf
  ! Precision in interatomic surfaces. 
  real(kind=rp), public :: epsilon
  ! If a non nuclear maximum is found whose distance to any nucleus 
  ! is larger than EPSISCP, promolden STOPS. This control of non nuclear
  ! maxima is performed in the 'iscp.f' routine. When Hydrogen atoms
  ! are involved it is very convenient to change the default value
  ! for EPSISCP (0.08) to a smaller value by using the 'EPSISCP value'
  ! order in the input file.
  real(kind=rp), public :: epsiscp
  ! NTRIAL  = number of sub-intervals in the search of the surface.
  ! RPRIMER = First point in the in the search of the surface.
  ! GEOFAC  = ((rmaxsurf-0.1d0)/rprimer)**(1d0/(ntrial-1))
  !           (RMAXSURF is defined in passed in integ.inc)
  integer(kind=ip), public :: ntrial
  real(kind=rp), public :: rprimer
  real(kind=rp), public :: geofac
  integer(kind=ip), public :: steeper

  public :: allocate_space_for_surface, deallocate_space_for_surface
  public :: allocate_space_for_cp, deallocate_space_for_cp

contains

  subroutine allocate_space_for_cp (ncent)
  
    use mod_memory, only: alloc
    implicit none
    integer(kind=ip), intent(in) :: ncent

    call alloc ('mod_surf', 'xyzrho', xyzrho, ncent, 3)
 
  end subroutine allocate_space_for_cp

  subroutine deallocate_space_for_cp ()

    use mod_memory, only: free
    implicit none

    call free ('mod_surf', 'xyzrho', xyzrho)

  end subroutine deallocate_space_for_cp

  subroutine allocate_space_for_surface (nangular,ncutsurf)

    use mod_memory, only: alloc
    implicit none
    integer(kind=ip), intent(in) :: nangular, ncutsurf

    call alloc ('mod_surf', 'rlimsurf', rlimsurf, nangular, ncutsurf)
    call alloc ('mod_surf', 'nlimsurf', nlimsurf, nangular)

  end subroutine allocate_space_for_surface

  subroutine deallocate_space_for_surface ()

    use mod_memory, only: free
    implicit none

    call free ('mod_surf', 'rlimsurf', rlimsurf)
    call free ('mod_surf', 'nlimsurf', nlimsurf)

  end subroutine deallocate_space_for_surface

end module mod_surf
