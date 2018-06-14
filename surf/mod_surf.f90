module mod_surf

  use mod_prec, only: rp, ip
  implicit none
  private

  integer(kind=ip), parameter, public :: minter = 10 
  integer(kind=ip), parameter, public :: maxtrial = 13
  ! The MAXSTART parameter has to do with the RSEARCH order, 
  ! RSEARCH nsearch nrstart, (rstart(i),i=1,nrstart)
  ! This order can be used when one wants to avoid using the
  ! default starting values of the radial coordinate in
  ! the initial search of interatomic surfaces. If nsearch is 0
  ! the starting values (rstart(i),i=1,nrstart) are used for all
  ! the atoms. If nsearch is different from 0, the rstart(i)
  ! (i=1,nrstart) values are used only for the atom 'nsearch'
  ! and all its symmetry equivalent atoms.
  integer(kind=ip), parameter, public :: maxstart = 40

  !integer(kind=ip), public :: inuc
  !real(kind=rp), public :: xnuc(3)
  real(kind=rp), allocatable, dimension(:,:), public :: xyzrho
  real(kind=rp), allocatable, dimension(:), public :: rmaxsurf
  integer(kind=ip), allocatable, dimension(:,:), public :: nangleb
  integer(kind=ip), public :: neqsurf
  integer(kind=ip), allocatable, dimension(:), public :: insurf
  real(kind=rp), allocatable, dimension(:,:), public :: rstart
  integer(kind=ip), allocatable, dimension(:), public :: nrsearch
  logical, allocatable, dimension(:), public :: lstart
  real(kind=rp), allocatable, dimension(:,:), public :: rlimsurf
  integer(kind=ip), allocatable, dimension(:), public :: nlimsurf

  ! options
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
  logical, public :: rotgrid
  real(kind=rp), public :: angx, angy, angz
  logical, public :: sphere
  integer(kind=ip), public :: steeper

  public :: init_surf
  public :: allocate_space_for_surface, deallocate_space_for_surface
  public :: allocate_space_for_integ, deallocate_space_for_integ

contains

  subroutine init_surf ()

    implicit none

    sphere = .false.
    steeper = 1_ip
    rotgrid = .false.
    ntrial = 11_ip
    rprimer = 0.4_rp
    inuc = 0_ip
    epsilon = 1d-5 
    epsiscp = 0.08_rp
    angx = 0.0_rp
    angy = 0.0_rp
    angz = 0.0_rp

  end subroutine
                                                                        
  subroutine allocate_space_for_integ (ncent)

    use mod_memory, only: alloc
    implicit none

    integer(kind=ip) :: i
    integer(kind=ip), intent(in) :: ncent

    neqsurf = ncent
    call alloc ('mod_surf', 'insurf', insurf, ncent)
    forall (i=1:neqsurf) insurf(i) = i
    call alloc ('mod_surf', 'xyzrho', xyzrho, ncent, 3)
    call alloc ('mod_surf', 'nangleb', nangleb, ncent, 4)
    nangleb(:,1) = 1
    nangleb(:,2) = 434
    nangleb(:,3) = 0
    nangleb(:,4) = 0
    call alloc ('mod_surf', 'rmaxsurf', rmaxsurf, ncent)
    rmaxsurf = 10.0_rp
    call alloc ('mod_surf', 'rstart', rstart, ncent, maxstart)
    call alloc ('mod_surf', 'nrsearch', nrsearch, ncent)
    ! By default no explicit values of the radial coordinate are given to 
    ! the atoms in the initial search of their interatomic surfaces.
    call alloc ('mod_surf', 'lstart', lstart, ncent)
    lstart = .false.

  end subroutine allocate_space_for_integ

  subroutine deallocate_space_for_integ ()

    use mod_memory, only: free
    implicit none

    call free ('mod_surf', 'xyzrho', xyzrho)
    call free ('mod_surf', 'rmaxsurf', rmaxsurf)
    call free ('mod_surf', 'nangleb', nangleb)
    call free ('mod_surf', 'rstart', rstart)
    call free ('mod_surf', 'nrsearch', nrsearch)
    call free ('mod_surf', 'lstart', lstart)
    call free ('mod_surf', 'insurf', insurf)

  end subroutine deallocate_space_for_integ

  subroutine allocate_space_for_surface (nangular,ncutsurf)

    use mod_memory, only: alloc
    implicit none
    integer(kind=ip), intent(in) :: nangular, ncutsurf

    call alloc ('mod_surf', 'rlimsurf', rlimsurf, nangular, ncutsurf)
    call alloc ('mod_surf', 'nlimsurf', nlimsurf, nangular)

  end subroutine allocate_space_for_surface

  subroutine deallocate_space_for_surface

    use mod_memory, only: free
    implicit none

    call free ('mod_surf', 'rlimsurf', rlimsurf)
    call free ('mod_surf', 'nlimsurf', nlimsurf)

  end subroutine deallocate_space_for_surface

end module mod_surf
