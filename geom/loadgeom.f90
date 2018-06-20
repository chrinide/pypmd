subroutine loadgeom()

  use mod_param, only: verbose
  implicit none

  call cmcq ()
  if (verbose) call infogeom ()

end subroutine loadgeom
