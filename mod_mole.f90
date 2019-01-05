module mod_mole
  
  use mod_prec, only: ip, rp
  implicit none
  public
  
  integer(kind=ip) :: ncent_
  real(kind=rp), allocatable, dimension(:,:) :: coords_
  integer(kind=ip), allocatable, dimension(:) :: charges_
  
contains

  subroutine allocate_space_for_mole ()
  
    use mod_memory, only: alloc
    implicit none

    call alloc ('mod_mole', 'coords_', coords_, 3, ncent_)
    call alloc ('mod_mole', 'charges_', charges_, ncent_)
 
  end subroutine allocate_space_for_mole

  subroutine deallocate_space_for_mole

    use mod_memory, only: free
    implicit none

    call free ('mod_mole', 'coords_', coords_)
    call free ('mod_mole', 'charges_', charges_)

  end subroutine deallocate_space_for_mole
                                                                        
end module mod_mole
