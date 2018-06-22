! Determine the surface
subroutine surface (filename)

  !$ use omp_lib, only: omp_get_wtime
  use mod_prec, only: rp, i4=>ip
  use mod_io, only: mline
  use mod_surf, only: minter
  implicit none
 
  character(len=*), intent(in) :: filename
 
  integer(kind=i4) :: lsu, ltxt
  character(len=mline) :: files

  interface
    subroutine ray (cot,sit,cop,sip,rsurf,nsurf)
      import rp, i4, minter
      integer(kind=i4), intent(out) :: nsurf
      real(kind=rp), intent(in) :: cot, sit, cop, sip
      real(kind=rp), intent(out) :: rsurf(minter,2) 
    end subroutine
    subroutine weightheta (iq,x,w,n)
      import rp, i4
      integer(kind=i4), intent(in) :: iq, n
      real(kind=rp), intent(out) :: x(n), w(n)
    end subroutine
    subroutine lebgrid (ct,st,cp,sp,w,npoints)
      import rp, i4
      integer(kind=i4), intent(in) :: npoints
      real(kind=rp), intent(out) :: ct(npoints)
      real(kind=rp), intent(out) :: st(npoints)
      real(kind=rp), intent(out) :: cp(npoints)
      real(kind=rp), intent(out) :: sp(npoints)
      real(kind=rp), intent(out) :: w(npoints)
    end subroutine
  end interface
      
  ! Init
  ltxt = 999
  lsu = 998
  files = trim(filename)//".surf"
  
  ! Find nucleus
  call findnuc ()

end subroutine surface
