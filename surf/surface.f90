! Copyright (c) 2018 
! Jose Luis Casals Sainz <jluiscasalssainz@gmail.com>,
! Ángel Martín Pendás <angel@fluor.quimica.uniovi.es> and 
! Evelio Francisco Miguelez <evelio@uniovi.es>. 
!
! promolden is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or (at
! your option) any later version.
! 
! promolden is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! Determine the surface
!
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
