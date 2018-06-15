subroutine end_surf ()

  use mod_surf, only:  deallocate_space_for_cp    
  implicit none

  call deallocate_space_for_cp ()

end subroutine
