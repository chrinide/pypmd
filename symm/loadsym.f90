subroutine loadsym ()

  use mod_param, only: verbose
  implicit none

  call rwcoord (.true., .true.)
  call detersym ( )
  call d2sym ()
  if (verbose) call infosym ()
  call noisesym ()
  if (verbose) call infosym ()

end subroutine loadsym

