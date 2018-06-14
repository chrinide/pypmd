subroutine end_wfn()

  use mod_wfn, only: deallocate_space_for_wfn, &
                     deallocate_space_for_shells, &
                     deallocate_space_for_rdm, &
                     deallocate_space_for_rho, rdm
  implicit none

  call deallocate_space_for_wfn ()
  call deallocate_space_for_shells ()
  if (rdm) call deallocate_space_for_rdm ()
  call deallocate_space_for_rho ()

end subroutine end_wfn
