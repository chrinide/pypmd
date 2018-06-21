subroutine end_wfn()

  use mod_wfn, only: deallocate_space_for_wfn, &
                     deallocate_space_for_shells, &
                     deallocate_space_for_rho, rdm, &
                     deallocate_space_for_rdm
  implicit none

  call deallocate_space_for_wfn ()
  call deallocate_space_for_shells ()
  call deallocate_space_for_rho ()
  if (rdm) call deallocate_space_for_rdm ()

end subroutine end_wfn
