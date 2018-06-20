subroutine loadgeom()

  use mod_param, only: verbose
  implicit none

  call cmcq ()
  if (verbose) call infogeom ()
  call connect ()

end subroutine loadgeom
