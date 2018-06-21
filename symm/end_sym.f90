subroutine end_sym()

  use mod_memory, only: free
  use mod_sym, only: xyzcom
  implicit none
    
  call free ('mod_sym', 'xyzcom', xyzcom)

end subroutine end_sym
